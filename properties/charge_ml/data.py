from pathlib import Path

import numpy as np
import torch
from ase.data import atomic_numbers


def numpy_to_torch(x) -> torch.Tensor:
    """Convertit un tableau numpy/list en torch.Tensor.

    Les tableaux flottants sont coulés en float32, les autres types
    (entiers en particulier) sont conservés tels quels.
    """
    arr = np.asarray(x)
    if np.issubdtype(arr.dtype, np.floating):
        arr = arr.astype(np.float32, copy=False)
    return torch.from_numpy(np.ascontiguousarray(arr))


def _parse_xyz_frames(path: Path):
    """Itère sur les frames (Z, R, q_ref|None) d'un xyz multi-frame.

    Les charges atomiques sont attendues en 5e colonne. Si elle est absente
    ou vide pour un atome, q_ref est None pour la frame. La ligne de
    commentaire n'est pas exploitée.
    """
    with open(path, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        n = int(line)

        Z = np.empty(n, dtype=np.int64)
        R = np.empty((n, 3), dtype=np.float32)
        q = np.empty(n, dtype=np.float32)
        has_q = True

        for k in range(n):
            parts = lines[i + 2 + k].split()
            if len(parts) < 4:
                raise ValueError(
                    f"{path}: ligne {i + 2 + k + 1} mal formée"
                )
            Z[k] = atomic_numbers[parts[0]]
            R[k] = (float(parts[1]), float(parts[2]), float(parts[3]))
            if len(parts) >= 5 and parts[4] != "":
                try:
                    q[k] = float(parts[4])
                except ValueError:
                    has_q = False
            else:
                has_q = False

        yield Z, R, (q if has_q else None)
        i += 2 + n


def load_xyz_dataset(folder: str | Path,
                     target: str = "q",
                     reference: dict[int, float] | None = None) -> list[dict]:
    """Charge un dataset xyz (charges en 5e colonne).

    Paramètres
    ----------
    folder : dossier contenant des fichiers .xyz (mono ou multi-frame).
    target : 'q' (par défaut) pour entraîner sur les charges, 'dq' pour
        entraîner sur la variation q - reference[Z].
    reference : dict {Z: charge_ref} requis si target='dq'.
    """
    if target not in ("q", "dq"):
        raise ValueError(f"target inconnu: {target!r}")
    if target == "dq" and reference is None:
        raise ValueError("target='dq' requiert un dict 'reference'")

    folder = Path(folder)
    samples: list[dict] = []
    for path in sorted(folder.glob("*.xyz")):
        for Z, R, q in _parse_xyz_frames(path):
            if q is None:
                raise ValueError(
                    f"{path}: charges manquantes en 5e colonne (q_ref=None)"
                )
            if target == "dq":
                ref = np.array([reference[int(z)] for z in Z],
                               dtype=np.float32)
                q = q - ref
            Q_tot = float(q.sum())
            samples.append(dict(Z=Z, R=R, Q_tot=Q_tot, q_ref=q))
    return samples


def split_data_ttv(arrays: dict,
                   percentages: tuple[float, float, float] = (80, 10, 10),
                   seed: int | None = None) -> dict:
    """Sépare un ou plusieurs tableaux en train/test/validation.

    Le découpage se fait toujours sur la première dimension. Tous les
    tableaux fournis doivent avoir la même longueur sur cet axe et sont
    mélangés selon la même permutation.

    Paramètres
    ----------
    arrays : dict {nom: tableau} (numpy.ndarray, list, ...).
    percentages : (p_train, p_test, p_valid). Somme <= 100.
    seed : graine pour le shuffle (optionnel).

    Renvoie
    -------
    dict {'train': {nom: ...}, 'test': {nom: ...}, 'valid': {nom: ...}}.
    """
    if not arrays:
        raise ValueError("aucun tableau fourni")
    if sum(percentages) > 100 + 1e-9:
        raise ValueError(f"somme des pourcentages > 100: {sum(percentages)}")
    sizes = {name: len(arr) for name, arr in arrays.items()}
    n = next(iter(sizes.values()))
    if any(s != n for s in sizes.values()):
        raise ValueError(f"tableaux de longueurs différentes: {sizes}")

    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)

    p1, p2, p3 = percentages
    n1 = int(round(p1 * n / 100.0))
    n2 = int(round(p2 * n / 100.0))
    n3 = int(round(p3 * n / 100.0))
    if n1 + n2 + n3 > n:
        n3 = max(0, n - n1 - n2)

    idx_train = perm[:n1]
    idx_test = perm[n1:n1 + n2]
    idx_valid = perm[n1 + n2:n1 + n2 + n3]

    def _take(arr, idx):
        if isinstance(arr, np.ndarray):
            return arr[idx]
        return [arr[i] for i in idx]

    out = {"train": {}, "test": {}, "valid": {}}
    for name, arr in arrays.items():
        out["train"][name] = _take(arr, idx_train)
        out["test"][name] = _take(arr, idx_test)
        out["valid"][name] = _take(arr, idx_valid)
    return out


