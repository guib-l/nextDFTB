"""Exemple de prédiction de charges atomiques avec un modèle entraîné.

Charge le checkpoint produit par exemple_water.py, prend une frame du
fichier h2o_data_01.xyz, et compare la prédiction à la référence.
"""

from pathlib import Path

import numpy as np

from charge_ml import ChargeML
from charge_ml.data import _parse_xyz_frames,load_xyz_dataset
from charge_ml.default import DEFAULT_CHARGES


HERE = Path(__file__).resolve().parent
MODEL_PATH = HERE / "model" / "checkpoint.pt"
XYZ_PATH = HERE / "data" 

FRAME_INDEX = 1587


def main() -> None:
    if not MODEL_PATH.exists():
        raise FileNotFoundError(
            f"aucun modèle à {MODEL_PATH}: lancez exemple_water.py d'abord"
        )

    ml = ChargeML()
    ml.load_model(MODEL_PATH)
    print(f"modèle chargé (target={ml.target}, elements={ml.elements})")

    samples = load_xyz_dataset(XYZ_PATH, target="q", reference=DEFAULT_CHARGES)
    s = samples[FRAME_INDEX]
    Z, R, q_ref = s["Z"], s["R"], s["q_ref"]

    Q_tot = float(q_ref.sum())
    if ml.target == "dq" and ml.reference is not None:
        ref = np.array([ml.reference[int(z)] for z in Z], dtype=np.float32)
        q_ref = q_ref - ref
        Q_tot = float(q_ref.sum())

    q_pred = ml.predict(Z, R, Q_tot=Q_tot)

    print(f"frame : {len(Z)} atomes, Q_tot={Q_tot:+.4f}")
    print(f"{'i':>3} {'sym':>3} {'q_ref':>10} {'q_pred':>10} {'err':>10}")
    for i in range(len(Z)):
        sym = "H" if int(Z[i]) == 1 else "O"
        print(f"{i:3d} {sym:>3s} "
              f"{float(q_ref[i]):>+10.4f} "
              f"{float(q_pred[i]):>+10.4f} "
              f"{float(q_pred[i] - q_ref[i]):>+10.4f}")
    mae = float(np.mean(np.abs(q_pred - q_ref)))
    print(f"MAE = {mae:.4e}")


if __name__ == "__main__":
    main()
