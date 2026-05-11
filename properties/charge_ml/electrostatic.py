import numpy as np
import torch
from torch import Tensor

from ase.data import atomic_masses

from .default import DEFAULT_HUBBARD


HUBBARD_DEFAULT: dict[int, float] = dict(DEFAULT_HUBBARD)


def Gamma(rab: float, Ua: float, Ub: float) -> float:
    """Noyau gamma DFTB (variante stdr) entre deux atomes A et B.

    Calcule la valeur scalaire gamma_AB pour une distance rab et des
    paramètres de Hubbard Ua, Ub. Les valeurs renvoyées sont positives.
    """
    def _gamma_e(r: float, a: float, b: float) -> float:
        d = a * a - b * b
        return np.exp(-a * r) * (
            (0.5 * b ** 4 * a) / (d ** 2)
            - (b ** 6 - 3.0 * b ** 4 * a ** 2) / (r * d ** 3)
        )

    eq = abs(Ua - Ub) < 1e-6
    tauA, tauB = 3.2 * Ua, 3.2 * Ub

    if rab < 1e-6:
        if eq:
            return 0.5 * (Ua + Ub)
        return 0.5 * (
            (tauA * tauB) / (tauA + tauB)
            + (tauA * tauB) ** 2 / (tauA + tauB) ** 3
        )
    if eq:
        tauM = 0.5 * (tauA + tauB)
        return float(
            np.exp(-tauM * rab) * (
                1.0 / rab
                + 0.6875 * tauM
                + 0.1875 * rab * tauM ** 2
                + 0.0208333333333 * rab ** 2 * tauM ** 3
            )
        )
    return float(_gamma_e(rab, tauA, tauB) + _gamma_e(rab, tauB, tauA))


def calc_gamma(R: Tensor, Z, hubbard: dict[int, float]) -> Tensor:
    """Construit la matrice gamma (n, n) symétrique via Gamma(rab, Ua, Ub)."""
    R_t = R if torch.is_tensor(R) else torch.as_tensor(R)
    Z_t = torch.as_tensor(Z)
    n = R_t.shape[0]
    R_np = R_t.detach().cpu().numpy()
    u_list = [hubbard[int(z)] for z in Z_t]
    G = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i, n):
            rij = float(np.linalg.norm(R_np[i] - R_np[j]))
            if rij > 1e-7:
                v = 1/rij - Gamma(rij, u_list[i], u_list[j])
            else:
                v = Gamma(rij, u_list[i], u_list[j])
            G[i, j] = v
            G[j, i] = v
    return torch.from_numpy(G).to(dtype=R_t.dtype, device=R_t.device)


def coulomb_energy(q: Tensor, R: Tensor, Z: Tensor,
                   hubbard: dict[int, float]) -> Tensor:
    g = calc_gamma(R, Z, hubbard).to(q.dtype)
    n = q.shape[0]
    mask = 1.0 - torch.eye(n, dtype=g.dtype, device=g.device)
    qq = q.unsqueeze(0) * q.unsqueeze(1)
    return 0.5 * (qq * g * mask).sum()


def dipole(q: Tensor, R: Tensor, masses: Tensor | None = None) -> Tensor:
    if masses is None:
        Rcm = R.mean(dim=0)
    else:
        w = masses / masses.sum()
        Rcm = (w.unsqueeze(-1) * R).sum(dim=0)
    return (q.unsqueeze(-1) * (R - Rcm)).sum(dim=0)


def masses_from_Z(Z: Tensor) -> Tensor:
    return torch.tensor([atomic_masses[int(z)] for z in Z], dtype=torch.float64)


def qeq_charges(R, Z, Q_tot: float,
                chi: "dict[int, float] | Tensor",
                hardness: "dict[int, float] | None" = None,
                hubbard: "dict[int, float] | None" = None) -> Tensor:
    """Charges QEq pour un chi fourni par l'utilisateur (sans NN).

    Paramètres
    ----------
    R, Z : géométrie (positions et numéros atomiques).
    Q_tot : charge totale imposée.
    chi : électronégativité par espèce ({Z: chi}) ou par atome (tensor (n,)).
    hardness, hubbard : valeurs par défaut DEFAULT_HUBBARD si None.
    """
    Z_t = torch.as_tensor(Z, dtype=torch.long)
    R_t = R if torch.is_tensor(R) else torch.as_tensor(R, dtype=torch.float32)
    if isinstance(chi, dict):
        chi_t = torch.tensor([chi[int(z)] for z in Z_t], dtype=R_t.dtype)
    else:
        chi_t = chi.to(R_t.dtype) if torch.is_tensor(chi) \
                else torch.as_tensor(chi, dtype=R_t.dtype)
    if hardness is None:
        hardness = dict(HUBBARD_DEFAULT)
    if hubbard is None:
        hubbard = dict(HUBBARD_DEFAULT)
    return qeq_solve(R_t, Z_t, float(Q_tot), chi_t, hardness, hubbard)


def qeq_solve(R: Tensor, Z, Q_tot: float, chi: Tensor,
              hardness: dict[int, float],
              hubbard: dict[int, float]) -> Tensor:
    Z_t = torch.as_tensor(Z)
    n = R.shape[0]
    g = calc_gamma(R, Z_t, hubbard).to(chi.dtype)
    #diag = torch.tensor([hardness[int(z)] for z in Z_t],
    #                    dtype=chi.dtype, device=chi.device)
    A = g #- torch.diag(torch.diagonal(g)) + torch.diag(diag)
    M = torch.zeros((n + 1, n + 1), dtype=chi.dtype, device=chi.device)
    M[:n, :n] = A
    M[:n, n] = 1.0
    M[n, :n] = 1.0
    rhs = torch.zeros(n + 1, dtype=chi.dtype, device=chi.device)
    rhs[:n] = -chi
    rhs[n] = float(Q_tot)
    sol = torch.linalg.solve(M, rhs)
    return sol[:n]
