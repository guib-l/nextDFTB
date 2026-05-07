"""Interface principale du module charge_ml (objet ChargeML)."""

from pathlib import Path

import numpy as np
import torch
from ase.data import chemical_symbols

from .data import load_xyz_dataset, numpy_to_torch, split_data_ttv
from .default import (
    DEFAULT_CHARGES,
    DEFAULT_DESC,
    DEFAULT_EPOCH,
    DEFAULT_HUBBARD,
)
from .display import print_epoch
from .features import build_descriptor, featurize
from .model import ChargeModel, ElectronegativityNet
from .train import _evaluate, _train_one_epoch


class ChargeML:
    """Interface haut-niveau de la librairie charge_ml."""

    def __init__(
        self,
        elements: list[str] | None = None,
        hubbard: dict[int, float] | None = None,
        hardness: dict[int, float] | None = None,
        desc_kind: str | None = None,
        desc_kwargs: dict | None = None,
        target: str = "q",
        reference: dict[int, float] | None = None,
        device: torch.device | str | None = None,
        hidden: int = 64,
        seed: int = 0,
    ):
        if target not in ("q", "dq"):
            raise ValueError(f"target inconnu: {target!r}")

        self.elements = list(elements) if elements else []
        self.hubbard = dict(hubbard) if hubbard else dict(DEFAULT_HUBBARD)
        self.hardness = dict(hardness) if hardness else dict(self.hubbard)
        self.desc_kind = desc_kind if desc_kind else DEFAULT_DESC["kind"]
        self.desc_kwargs = (
            dict(desc_kwargs) if desc_kwargs is not None
            else {k: v for k, v in DEFAULT_DESC.items() if k != "kind"}
        )
        self.target = target
        self.reference = (
            dict(reference) if reference is not None
            else (dict(DEFAULT_CHARGES) if target == "dq" else None)
        )
        if device is None:
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.device = torch.device(device)
        self.hidden = hidden
        self.seed = seed

        self.desc = None
        self.model: ChargeModel | None = None
        self.samples: list[dict] = []
        self.feats: list[np.ndarray] = []
        self.splits: dict | None = None
        self.history: dict[str, list[float]] = {}

    # ----------------------- I/O données -----------------------

    def load_data(self, folder: str | Path) -> list[dict]:
        """Charge les fichiers .xyz d'un dossier."""
        self.samples = load_xyz_dataset(
            folder, target=self.target, reference=self.reference,
        )
        if not self.elements:
            seen = sorted({
                chemical_symbols[int(z)]
                for s in self.samples for z in s["Z"]
            })
            self.elements = seen
        return self.samples

    def save_data(self, folder: str | Path) -> dict[str, Path]:
        """Sauve les splits courants en 3 fichiers xyz train/test/valid."""
        if self.splits is None:
            raise ValueError("aucun split: appelez build_dataset d'abord")
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        paths: dict[str, Path] = {}
        for name in ("train", "test", "valid"):
            path = folder / f"{name}.xyz"
            self._write_xyz(path, self.splits[name]["samples"])
            paths[name] = path
        return paths

    @staticmethod
    def _write_xyz(path: Path, samples: list[dict]) -> None:
        with open(path, "w") as f:
            for s in samples:
                Z = s["Z"]
                R = s["R"]
                q = s["q_ref"]
                n = len(Z)
                f.write(f"{n}\n\n")
                for k in range(n):
                    sym = chemical_symbols[int(Z[k])]
                    f.write(
                        f"{sym} {float(R[k, 0]):.6f} "
                        f"{float(R[k, 1]):.6f} {float(R[k, 2]):.6f} "
                        f"{float(q[k]):.6f}\n"
                    )

    # ----------------------- Dataset -----------------------

    def build_dataset(
        self,
        percentages: tuple[float, float, float] = (80, 10, 10),
        seed: int | None = None,
    ) -> dict:
        """Calcule descripteurs et split en train/test/valid."""
        if not self.samples:
            raise ValueError("aucune donnée chargée: appelez load_data d'abord")
        if not self.elements:
            raise ValueError("self.elements doit être défini")
        if self.desc is None:
            self.desc = build_descriptor(
                self.desc_kind, self.elements, **self.desc_kwargs,
            )
        self.feats = [
            featurize(self.desc, s["Z"], s["R"]) for s in self.samples
        ]
        self.splits = split_data_ttv(
            {"samples": self.samples, "feats": self.feats},
            percentages=percentages,
            seed=seed if seed is not None else self.seed,
        )
        return self.splits

    def save_dataset(self, folder: str | Path) -> dict[str, Path]:
        """Sauve le dataset (samples + descripteurs) en 3 fichiers npz."""
        if self.splits is None:
            raise ValueError("aucun split: appelez build_dataset d'abord")
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        paths: dict[str, Path] = {}
        for name in ("train", "test", "valid"):
            path = folder / f"{name}.npz"
            samples = self.splits[name]["samples"]
            feats = self.splits[name]["feats"]
            np.savez(
                path,
                Z=np.array([s["Z"] for s in samples], dtype=object),
                R=np.array([s["R"] for s in samples], dtype=object),
                q=np.array([s["q_ref"] for s in samples], dtype=object),
                Q_tot=np.array([s["Q_tot"] for s in samples], dtype=np.float32),
                feats=np.array(feats, dtype=object),
            )
            paths[name] = path
        return paths

    # ----------------------- Modèle -----------------------

    def save_model(
        self,
        folder: str | Path = "model",
        filename: str = "checkpoint.pt",
    ) -> Path:
        """Sauve le modèle dans un répertoire (créé par défaut: model/)."""
        if self.model is None:
            raise ValueError("aucun modèle à sauvegarder")
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / filename
        ckpt = dict(
            state_dict=self.model.state_dict(),
            d_in=self.model.chi_net.d_in,
            hidden=self.model.chi_net.hidden,
            elements=list(self.elements),
            desc_kind=self.desc_kind,
            desc_kwargs=self.desc_kwargs,
            hardness=self.hardness,
            hubbard=self.hubbard,
            target=self.target,
            reference=self.reference,
            history=self.history,
        )
        torch.save(ckpt, path)
        return path

    def load_model(self, path: str | Path) -> ChargeModel:
        """Charge un modèle depuis un checkpoint."""
        ckpt = torch.load(path, map_location="cpu", weights_only=False)
        self.elements = list(ckpt["elements"])
        self.desc_kind = ckpt["desc_kind"]
        self.desc_kwargs = dict(ckpt["desc_kwargs"])
        self.hubbard = dict(ckpt["hubbard"])
        self.hardness = dict(ckpt["hardness"])
        self.target = ckpt.get("target", "q")
        self.reference = ckpt.get("reference", None)
        self.history = ckpt.get("history", {})
        self.desc = build_descriptor(
            self.desc_kind, self.elements, **self.desc_kwargs,
        )
        chi_net = ElectronegativityNet(
            d_in=ckpt["d_in"], elements=self.elements, hidden=ckpt["hidden"],
        )
        self.model = ChargeModel(
            chi_net, hardness=self.hardness, hubbard=self.hubbard,
        )
        self.model.load_state_dict(ckpt["state_dict"])
        self.model.to(self.device)
        return self.model

    # ----------------------- Entraînement (3 méthodes) -----------------------

    def train(
        self,
        epochs: int = DEFAULT_EPOCH,
        lr: float = 1e-3,
        patience: int = 20,
        verbose: bool = True,
    ) -> dict:
        """Entraîne sur le set train; arrêt anticipé sur la métrique test."""
        if self.splits is None:
            raise ValueError("dataset non préparé: appelez build_dataset d'abord")
        torch.manual_seed(self.seed)
        rng = np.random.default_rng(self.seed)

        if self.model is None:
            d_in = self.feats[0].shape[1]
            chi_net = ElectronegativityNet(
                d_in=d_in, elements=self.elements, hidden=self.hidden,
            )
            self.model = ChargeModel(
                chi_net, hardness=self.hardness, hubbard=self.hubbard,
            )
        self.model.to(self.device)
        opt = torch.optim.Adam(self.model.parameters(), lr=lr)

        train_s = self.splits["train"]["samples"]
        train_f = self.splits["train"]["feats"]

        history = {"train_mse": [], "test_mse": [], "test_mae": []}
        best_test = float("inf")
        best_state = {
            k: v.detach().clone() for k, v in self.model.state_dict().items()
        }
        bad = 0

        for epoch in range(epochs):
            train_mse = _train_one_epoch(
                self.model, train_s, train_f, opt, rng, self.device,
            )
            test_mse, test_mae = self.test()
            history["train_mse"].append(train_mse)
            history["test_mse"].append(test_mse)
            history["test_mae"].append(test_mae)
            if verbose:
                print_epoch(epoch, dict(
                    train_mse=train_mse,
                    test_mse=test_mse,
                    test_mae=test_mae,
                ))
            if test_mse < best_test - 1e-8:
                best_test = test_mse
                best_state = {
                    k: v.detach().clone()
                    for k, v in self.model.state_dict().items()
                }
                bad = 0
            else:
                bad += 1
                if bad >= patience:
                    break

        self.model.load_state_dict(best_state)
        self.history = history
        return history

    def test(self) -> tuple[float, float]:
        """Évalue le modèle sur le set test (MSE, MAE)."""
        if self.splits is None or self.model is None:
            raise ValueError("dataset/model non prêt")
        return _evaluate(
            self.model,
            self.splits["test"]["samples"],
            self.splits["test"]["feats"],
            self.device,
        )

    def validate(self) -> tuple[float, float]:
        """Évalue le modèle sur le set valid (MSE, MAE)."""
        if self.splits is None or self.model is None:
            raise ValueError("dataset/model non prêt")
        return _evaluate(
            self.model,
            self.splits["valid"]["samples"],
            self.splits["valid"]["feats"],
            self.device,
        )

    # ----------------------- Inférence -----------------------

    def predict(self, Z, R, Q_tot: float = 0.0) -> np.ndarray:
        """Inférence pour une géométrie unique."""
        if self.model is None:
            raise ValueError("aucun modèle chargé")
        if self.desc is None:
            self.desc = build_descriptor(
                self.desc_kind, self.elements, **self.desc_kwargs,
            )
        Z = np.asarray(Z, dtype=np.int64)
        R = np.asarray(R, dtype=np.float32)
        X = featurize(self.desc, Z, R)
        X_t = numpy_to_torch(X).to(self.device)
        R_t = numpy_to_torch(R).to(self.device)
        Z_t = numpy_to_torch(Z).to(self.device)
        self.model.eval()
        with torch.no_grad():
            q = self.model(X_t, Z_t, R_t, float(Q_tot))
        return q.cpu().numpy()
