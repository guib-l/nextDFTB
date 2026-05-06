import torch
from torch import Tensor

from ase.data import atomic_masses


HUBBARD_DEFAULT: dict[int, float] = {
    1: 0.4196, 6: 0.3647, 7: 0.4309, 8: 0.4954,
}


def _hubbard_pair(Z: Tensor, hubbard: dict[int, float]) -> Tensor:
    U = torch.tensor([hubbard[int(z)] for z in Z], dtype=torch.float64)
    Ui = U.unsqueeze(0)
    Uj = U.unsqueeze(1)
    return 2.0 * Ui * Uj / (Ui + Uj)


def gamma_klopman(R: Tensor, Z: Tensor, hubbard: dict[int, float]) -> Tensor:
    R = R.to(torch.float64)
    Z = torch.as_tensor(Z)
    rij = torch.cdist(R, R)
    U_ij = _hubbard_pair(Z, hubbard).to(R.dtype)
    return 1.0 / torch.sqrt(rij ** 2 + (1.0 / U_ij) ** 2)


def coulomb_energy(q: Tensor, R: Tensor, Z: Tensor,
                   hubbard: dict[int, float]) -> Tensor:
    g = gamma_klopman(R, Z, hubbard)
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


def qeq_solve(R: Tensor, Z, Q_tot: float, chi: Tensor,
              hardness: dict[int, float],
              hubbard: dict[int, float]) -> Tensor:
    Z_t = torch.as_tensor(Z)
    n = R.shape[0]
    g = gamma_klopman(R, Z_t, hubbard).to(chi.dtype)
    diag = torch.tensor([hardness[int(z)] for z in Z_t],
                        dtype=chi.dtype, device=chi.device)
    A = g - torch.diag(torch.diagonal(g)) + torch.diag(diag)
    M = torch.zeros((n + 1, n + 1), dtype=chi.dtype, device=chi.device)
    M[:n, :n] = A
    M[:n, n] = 1.0
    M[n, :n] = 1.0
    rhs = torch.zeros(n + 1, dtype=chi.dtype, device=chi.device)
    rhs[:n] = -chi
    rhs[n] = float(Q_tot)
    sol = torch.linalg.solve(M, rhs)
    return sol[:n]
