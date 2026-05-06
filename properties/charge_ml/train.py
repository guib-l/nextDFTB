from pathlib import Path

import numpy as np
import torch
from ase.data import chemical_symbols

from .features import build_descriptor, featurize
from .model import ChargeModel, ElectronegativityNet
from .electrostatic import HUBBARD_DEFAULT


def train(
    samples: list[dict],
    elements: list[str],
    hardness: dict[int, float] | None = None,
    hubbard: dict[int, float] | None = None,
    desc_kind: str = "acsf",
    desc_kwargs: dict | None = None,
    hidden: int = 64,
    epochs: int = 200,
    lr: float = 1e-3,
    val_frac: float = 0.1,
    patience: int = 20,
    seed: int = 0,
    out: str | Path | None = None,
) -> dict:
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)

    desc_kwargs = desc_kwargs or dict(
        r_cut=5.0,
        g2_params=[[1.0, 1.0], [1.0, 2.0], [1.0, 3.0], [1.0, 4.0]],
        g4_params=[[1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, -1.0]],
    )
    desc = build_descriptor(desc_kind, elements, **desc_kwargs)

    hubbard = dict(hubbard) if hubbard is not None else dict(HUBBARD_DEFAULT)
    hardness = dict(hardness) if hardness is not None else dict(hubbard)

    feats = []
    for s in samples:
        symbols = [chemical_symbols[int(z)] for z in s["Z"]]
        if any(sym not in elements for sym in symbols):
            raise ValueError(f"unknown species in sample: {set(symbols) - set(elements)}")
        feats.append(featurize(desc, s["Z"], s["R"]))
    d_in = feats[0].shape[1]

    indices = np.arange(len(samples))
    rng.shuffle(indices)
    n_val = max(1, int(val_frac * len(indices)))
    val_idx = set(indices[:n_val].tolist())

    chi_net = ElectronegativityNet(d_in=d_in, elements=elements, hidden=hidden).double()
    model = ChargeModel(chi_net, hardness=hardness, hubbard=hubbard)
    opt = torch.optim.Adam(model.parameters(), lr=lr)

    best_val = float("inf")
    best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    bad = 0

    for epoch in range(epochs):
        model.train()
        order = [i for i in indices if i not in val_idx]
        rng.shuffle(order)
        for i in order:
            s = samples[i]
            X = torch.from_numpy(np.asarray(feats[i], dtype=np.float64))
            R = torch.from_numpy(np.asarray(s["R"], dtype=np.float64))
            Z = torch.from_numpy(np.asarray(s["Z"]))
            q_ref = torch.from_numpy(np.asarray(s["q_ref"], dtype=np.float64))
            q_pred = model(X, Z, R, float(s["Q_tot"]))
            loss = ((q_pred - q_ref) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()

        model.eval()
        with torch.no_grad():
            val_loss = 0.0
            for i in val_idx:
                s = samples[i]
                X = torch.from_numpy(np.asarray(feats[i], dtype=np.float64))
                R = torch.from_numpy(np.asarray(s["R"], dtype=np.float64))
                Z = torch.from_numpy(np.asarray(s["Z"]))
                q_ref = torch.from_numpy(np.asarray(s["q_ref"], dtype=np.float64))
                q_pred = model(X, Z, R, float(s["Q_tot"]))
                val_loss += float(((q_pred - q_ref) ** 2).mean())
            val_loss /= max(1, len(val_idx))

        if val_loss < best_val - 1e-8:
            best_val = val_loss
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            bad = 0
        else:
            bad += 1
            if bad >= patience:
                break

    model.load_state_dict(best_state)

    ckpt = dict(
        state_dict=model.state_dict(),
        d_in=d_in,
        hidden=hidden,
        elements=list(elements),
        desc_kind=desc_kind,
        desc_kwargs=desc_kwargs,
        hardness=hardness,
        hubbard=hubbard,
        best_val=best_val,
    )
    if out is not None:
        torch.save(ckpt, out)
    return ckpt
