from pathlib import Path

import numpy as np
import torch
from ase.data import chemical_symbols

from .data import numpy_to_torch, split_data_ttv
from .default import DEFAULT_DESC, DEFAULT_EPOCH, DEFAULT_HUBBARD
from .display import print_epoch
from .features import build_descriptor, featurize
from .metric import metric_mae, metric_mse
from .model import ChargeModel, ElectronegativityNet


def _evaluate(model, samples, feats_by_idx, device):
    """Renvoie (mse, mae) moyens sur la liste de samples."""
    model.eval()
    if not samples:
        return float("nan"), float("nan")
    mse_sum = 0.0
    mae_sum = 0.0
    with torch.no_grad():
        for i, s in enumerate(samples):
            X = numpy_to_torch(feats_by_idx[i]).to(device)
            R = numpy_to_torch(s["R"]).to(device)
            Z = numpy_to_torch(s["Z"]).to(device)
            q_ref = numpy_to_torch(s["q_ref"]).to(device)
            q_pred = model(X, Z, R, float(s["Q_tot"]))
            mse_sum += float(metric_mse(q_pred, q_ref))
            mae_sum += float(metric_mae(q_pred, q_ref))
    return mse_sum / len(samples), mae_sum / len(samples)


def _train_one_epoch(model, samples, feats_by_idx, opt, rng, device):
    model.train()
    order = list(range(len(samples)))
    rng.shuffle(order)
    loss_sum = 0.0
    for i in order:
        s = samples[i]
        X = numpy_to_torch(feats_by_idx[i]).to(device)
        R = numpy_to_torch(s["R"]).to(device)
        Z = numpy_to_torch(s["Z"]).to(device)
        q_ref = numpy_to_torch(s["q_ref"]).to(device)
        q_pred = model(X, Z, R, float(s["Q_tot"]))
        loss = metric_mse(q_pred, q_ref)
        opt.zero_grad()
        loss.backward()
        opt.step()
        loss_sum += float(loss.detach())
    return loss_sum / max(1, len(samples))


def train(
    samples: list[dict],
    elements: list[str],
    hardness: dict[int, float] | None = None,
    hubbard: dict[int, float] | None = None,
    desc_kind: str | None = None,
    desc_kwargs: dict | None = None,
    hidden: int = 64,
    epochs: int = DEFAULT_EPOCH,
    lr: float = 1e-3,
    percentages: tuple[float, float, float] = (80, 10, 10),
    patience: int = 20,
    seed: int = 0,
    device: torch.device | str | None = None,
    verbose: bool = True,
    out: str | Path | None = None,
) -> dict:
    """Entraîne ChargeModel et renvoie un checkpoint dict.

    Sépare les données en train/test/valid via split_data_ttv. Le set
    'test' sert de critère d'arrêt anticipé ; le set 'valid' est évalué
    une fois en fin d'entraînement.
    """
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)

    if desc_kwargs is None:
        desc_kwargs = {k: v for k, v in DEFAULT_DESC.items() if k != "kind"}
    if desc_kind is None:
        desc_kind = DEFAULT_DESC["kind"]
    desc = build_descriptor(desc_kind, elements, **desc_kwargs)

    hubbard = dict(hubbard) if hubbard is not None else dict(DEFAULT_HUBBARD)
    hardness = dict(hardness) if hardness is not None else dict(hubbard)

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    device = torch.device(device)

    feats: list[np.ndarray] = []
    for s in samples:
        symbols = [chemical_symbols[int(z)] for z in s["Z"]]
        unknown = set(symbols) - set(elements)
        if unknown:
            raise ValueError(f"unknown species in sample: {unknown}")
        feats.append(featurize(desc, s["Z"], s["R"]))
    d_in = feats[0].shape[1]

    splits = split_data_ttv(
        {"samples": samples, "feats": feats},
        percentages=percentages,
        seed=seed,
    )
    train_s, train_f = splits["train"]["samples"], splits["train"]["feats"]
    test_s, test_f = splits["test"]["samples"], splits["test"]["feats"]
    valid_s, valid_f = splits["valid"]["samples"], splits["valid"]["feats"]

    chi_net = ElectronegativityNet(d_in=d_in, elements=elements, hidden=hidden)
    model = ChargeModel(chi_net, hardness=hardness, hubbard=hubbard).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr)

    history: dict[str, list[float]] = {
        "train_mse": [], "test_mse": [], "test_mae": [],
    }
    best_test = float("inf")
    best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    bad = 0

    for epoch in range(epochs):
        train_mse = _train_one_epoch(model, train_s, train_f, opt, rng, device)
        test_mse, test_mae = _evaluate(model, test_s, test_f, device)
        history["train_mse"].append(train_mse)
        history["test_mse"].append(test_mse)
        history["test_mae"].append(test_mae)
        if verbose:
            print_epoch(epoch, dict(train_mse=train_mse,
                                    test_mse=test_mse, test_mae=test_mae))

        if test_mse < best_test - 1e-8:
            best_test = test_mse
            best_state = {k: v.detach().clone()
                          for k, v in model.state_dict().items()}
            bad = 0
        else:
            bad += 1
            if bad >= patience:
                break

    model.load_state_dict(best_state)
    valid_mse, valid_mae = _evaluate(model, valid_s, valid_f, device)
    if verbose:
        print(f"validation | mse={valid_mse:.4e}  mae={valid_mae:.4e}",
              flush=True)

    ckpt = dict(
        state_dict=model.state_dict(),
        d_in=d_in,
        hidden=hidden,
        elements=list(elements),
        desc_kind=desc_kind,
        desc_kwargs=desc_kwargs,
        hardness=hardness,
        hubbard=hubbard,
        history=history,
        best_test=best_test,
        valid_mse=valid_mse,
        valid_mae=valid_mae,
    )
    if out is not None:
        torch.save(ckpt, out)
    return ckpt
