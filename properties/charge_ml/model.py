import torch
from torch import nn, Tensor
from ase.data import atomic_numbers

from . import electrostatic as es


class ElectronegativityNet(nn.Module):
    def __init__(self, d_in: int, elements: list[str], hidden: int = 64):
        super().__init__()
        self.d_in = d_in
        self.hidden = hidden
        self.elements = list(elements)
        self.nets = nn.ModuleDict({
            sym: nn.Sequential(
                nn.Linear(d_in, hidden), nn.SiLU(),
                nn.Linear(hidden, hidden), nn.SiLU(),
                nn.Linear(hidden, 1),
            )
            for sym in elements
        })

    def forward(self, X: Tensor, Z: Tensor) -> Tensor:
        Z = torch.as_tensor(Z)
        chi = torch.zeros(X.shape[0], dtype=X.dtype, device=X.device)
        for sym in self.elements:
            zi = atomic_numbers[sym]
            idx = (Z == zi).nonzero(as_tuple=True)[0]
            if idx.numel():
                chi = chi.index_copy(0, idx, self.nets[sym](X[idx]).squeeze(-1))
        return chi


class ChargeModel(nn.Module):
    def __init__(self, chi_net: ElectronegativityNet,
                 hardness: dict[int, float],
                 hubbard: dict[int, float]):
        super().__init__()
        self.chi_net = chi_net
        self.hardness = dict(hardness)
        self.hubbard = dict(hubbard)

    def forward(self, X: Tensor, Z: Tensor, R: Tensor, Q_tot: float) -> Tensor:
        chi = self.chi_net(X, Z)
        return es.qeq_solve(R, Z, float(Q_tot), chi, self.hardness, self.hubbard)
