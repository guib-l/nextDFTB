from pathlib import Path

import numpy as np
import torch
from ase.data import chemical_symbols, atomic_masses

from .features import build_descriptor, featurize
from .model import ChargeModel, ElectronegativityNet
from . import electrostatic as es


class ChargePredictor:
    def __init__(self, model: ChargeModel, desc, elements: list[str]):
        self.model = model.eval()
        self.desc = desc
        self.elements = list(elements)

    @classmethod
    def load(cls, path: str | Path) -> "ChargePredictor":
        ckpt = torch.load(path, map_location="cpu", weights_only=False)
        desc = build_descriptor(ckpt["desc_kind"], ckpt["elements"], **ckpt["desc_kwargs"])
        chi_net = ElectronegativityNet(
            d_in=ckpt["d_in"], elements=ckpt["elements"], hidden=ckpt["hidden"],
        ).double()
        model = ChargeModel(chi_net, hardness=ckpt["hardness"], hubbard=ckpt["hubbard"])
        model.load_state_dict(ckpt["state_dict"])
        return cls(model, desc, ckpt["elements"])

    def _check(self, Z: np.ndarray, R: np.ndarray) -> None:
        symbols = {chemical_symbols[int(z)] for z in Z}
        unknown = symbols - set(self.elements)
        if unknown:
            raise ValueError(f"unknown species: {unknown}")
        if len(Z) >= 2:
            d = np.linalg.norm(R[:, None, :] - R[None, :, :], axis=-1)
            i, j = np.triu_indices(len(Z), k=1)
            if d[i, j].min() < 0.5:
                import warnings
                warnings.warn("min interatomic distance < 0.5 Å — geometry may be broken")

    def predict(self, Z: np.ndarray, R: np.ndarray, Q_tot: float = 0.0) -> np.ndarray:
        Z = np.asarray(Z)
        R = np.asarray(R, dtype=np.float64)
        self._check(Z, R)
        X = featurize(self.desc, Z, R)
        X_t = torch.from_numpy(np.asarray(X, dtype=np.float64))
        R_t = torch.from_numpy(R)
        Z_t = torch.from_numpy(Z)
        with torch.no_grad():
            q = self.model(X_t, Z_t, R_t, float(Q_tot))
        return q.numpy()

    def dipole(self, Z: np.ndarray, R: np.ndarray, Q_tot: float = 0.0) -> np.ndarray:
        q = self.predict(Z, R, Q_tot)
        Z = np.asarray(Z)
        R_t = torch.from_numpy(np.asarray(R, dtype=np.float64))
        q_t = torch.from_numpy(np.asarray(q, dtype=np.float64))
        masses = torch.tensor([atomic_masses[int(z)] for z in Z], dtype=torch.float64)
        return es.dipole(q_t, R_t, masses).numpy()

    def coulomb_energy(self, Z: np.ndarray, R: np.ndarray, Q_tot: float = 0.0) -> float:
        q = self.predict(Z, R, Q_tot)
        Z_t = torch.from_numpy(np.asarray(Z))
        R_t = torch.from_numpy(np.asarray(R, dtype=np.float64))
        q_t = torch.from_numpy(np.asarray(q, dtype=np.float64))
        return float(es.coulomb_energy(q_t, R_t, Z_t, self.model.hubbard))
