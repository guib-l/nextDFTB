"""Entraînement de charge_ml sur des agrégats d'eau (h2o_data_*.xyz).

Étapes :
  - chargement des données h2o_data_*.xyz à la racine du projet
  - sous-échantillonnage à 35 % des frames disponibles
  - découpage 80/10/10 train/test/valid
  - entraînement sur les charges partielles (target='dq')
  - validation finale
  - sauvegarde du modèle dans exemples/model/
"""

from pathlib import Path

import numpy as np

from charge_ml import ChargeML
from charge_ml.default import DEFAULT_CHARGES


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA_DIR = HERE / "data"
MODEL_DIR = HERE / "model"

SUBSAMPLE = 1.0
SEED = 0


def main() -> None:
    ml = ChargeML(
        elements=["H", "O"],
        target="dq",
        seed=SEED,
        hidden=32,
    )

    samples = ml.load_data(DATA_DIR)
    n_total = len(samples)
    rng = np.random.default_rng(SEED)
    n_keep = int(round(SUBSAMPLE * n_total))
    keep_idx = rng.choice(n_total, size=n_keep, replace=False)
    ml.samples = [samples[i] for i in keep_idx]

    print(f"{n_total} frames disponibles, {n_keep} conservées ({SUBSAMPLE:.0%})")

    splits = ml.build_dataset(percentages=(80, 10, 10), seed=SEED)
    print({k: len(v["samples"]) for k, v in splits.items()})

    ml.train(epochs=400, lr=1e-3, patience=5, verbose=True)

    val_mse, val_mae = ml.validate()
    print(f"validation finale | mse={val_mse:.4e}  mae={val_mae:.4e}")

    path = ml.save_model(MODEL_DIR)
    print(f"modèle sauvegardé : {path}")




if __name__ == "__main__":
    main()
