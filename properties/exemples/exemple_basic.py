"""Utilisation basique de la librairie charge_ml.

Démontre :
  - création des descripteurs (via build_dataset)
  - chargement des données et découpage train/test/validation
  - chargement d'un modèle quelconque depuis un checkpoint
  - sauvegarde des données et du dataset après traitement
"""

from pathlib import Path

from charge_ml import ChargeML


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]               # racine du projet nextDFTB
DATA_DIR = ROOT                       # h2o_data_*.xyz à la racine
OUT_DIR = HERE / "data" / "_out"      # sortie de save_data / save_dataset
MODEL_DIR = HERE / "model"


def main() -> None:
    # 1. Configuration de l'objet ChargeML
    ml = ChargeML(
        elements=["H", "O"],
        target="dq",
        seed=0,
        hidden=32,
    )

    # 2. Chargement des fichiers .xyz (charges en 5e colonne)
    samples = ml.load_data(DATA_DIR)
    print(f"frames chargées : {len(samples)}")

    # 3. Création des descripteurs ACSF + split 80/10/10
    splits = ml.build_dataset(percentages=(80, 10, 10), seed=0)
    print({k: len(v["samples"]) for k, v in splits.items()})

    # 4. Sauvegarde des données (xyz) et du dataset (npz) découpés
    ml.save_data(OUT_DIR / "xyz")
    ml.save_dataset(OUT_DIR / "npz")
    print("save_data / save_dataset OK")

    # 5. Chargement d'un modèle existant si présent
    ckpt = MODEL_DIR / "checkpoint.pt"
    if ckpt.exists():
        ml.load_model(ckpt)
        print(f"modèle chargé depuis {ckpt}")
    else:
        print(f"aucun modèle trouvé à {ckpt} (étape ignorée)")


if __name__ == "__main__":
    main()
