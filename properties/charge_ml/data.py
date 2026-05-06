from pathlib import Path

import numpy as np
from ase.io import read


def load_xyz_dataset(folder: str | Path) -> list[dict]:
    folder = Path(folder)
    samples = []
    for path in sorted(folder.glob("*.xyz")):
        for atoms in read(path, index=":"):
            info = atoms.info
            Q_tot = int(info.get("charge", 0))
            q_str = info.get("charges", None)
            q_ref = (
                np.fromstring(q_str.strip('"'), sep=" ")
                if isinstance(q_str, str)
                else np.array(info.get("charges", []), dtype=float)
            )
            samples.append(
                dict(
                    Z=atoms.numbers.copy(),
                    R=atoms.positions.copy(),
                    Q_tot=Q_tot,
                    q_ref=q_ref if q_ref.size else None,
                )
            )
    return samples


def load_npz_dataset(path: str | Path) -> list[dict]:
    data = np.load(path, allow_pickle=True)
    Z, R, q, Q_tot = data["Z"], data["R"], data["q"], data["Q_tot"]
    return [
        dict(Z=np.asarray(Z[i]), R=np.asarray(R[i]),
             Q_tot=int(Q_tot[i]), q_ref=np.asarray(q[i]))
        for i in range(len(Z))
    ]
