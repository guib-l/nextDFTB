"""Exemple : prédire les charges via QEq sans réseau de neurones.

L'utilisateur fournit une électronégativité chi par espèce ; on résout
le système QEq (matrice gamma DFTB + dureté Hubbard) pour obtenir les
charges atomiques sous contrainte Q_tot.
"""

from pathlib import Path

import numpy as np
import torch

from charge_ml.default import DEFAULT_CHARGES
from charge_ml import load_xyz_dataset, qeq_charges


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
FRAME_INDEX = 1587


def main() -> None:
    # 1. Chargement complet
    samples = load_xyz_dataset(DATA_DIR, target="dq", reference=DEFAULT_CHARGES)
    s = samples[FRAME_INDEX]
    Z, R, qref = s["Z"], s["R"], s["q_ref"]
    n = len(Z)
    Q_tot = 0.0  

    # 2. Chi fourni par l'utilisateur 
    chi_user = {
        1: 0.2,   # H
        8: 0.5,   # O
    }

    # 3. Résolution QEq 
    q_pred = qeq_charges(R, Z, Q_tot, chi=chi_user)
    q_pred_np = q_pred.detach().cpu().numpy()

    # 4. Affichage 
    print(f"frame : {n} atomes  Q_tot imposé = {Q_tot:+.4f}")
    print(f"sum(q_pred) = {q_pred_np.sum():+.6e}")
    print(f"{'i':>3} {'sym':>3} {'chi':>8} {'q_pred':>10} {'q_real':>10}")

    for i in range(n):
        sym = "H" if int(Z[i]) == 1 else "O"
        print(f"{i:3d} {sym:>3s} {chi_user[int(Z[i])]:>8.3f} ",
              f"{float(q_pred_np[i]):>+10.4f}",
              f"{float(qref[i]):>+10.4f}")

    # 5. Vérification 
    q_O = q_pred_np[np.asarray(Z) == 8]
    q_H = q_pred_np[np.asarray(Z) == 1]
    print(f"\n<q_O> = {q_O.mean():+.4f}   <q_H> = {q_H.mean():+.4f}")


if __name__ == "__main__":
    main()
