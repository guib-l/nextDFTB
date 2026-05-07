import numpy as np
import torch

from charge_ml.electrostatic import (
    coulomb_energy,
    dipole,
    calc_gamma,
    masses_from_Z,
)


HUBBARD = {1: 0.4196, 6: 0.3647, 7: 0.4309, 8: 0.4954}


def test_long_range_limit():
    Z = torch.tensor([1, 1])
    r = 50.0
    R = torch.tensor([[0.0, 0.0, 0.0], [r, 0.0, 0.0]], dtype=torch.float64)
    g = calc_gamma(R, Z, HUBBARD)
    assert abs(float(g[0, 1])) < 1e-3


def test_zero_distance_finite():
    Z = torch.tensor([8, 8])
    R = torch.zeros((2, 3), dtype=torch.float64)
    g = calc_gamma(R, Z, HUBBARD)
    assert torch.isfinite(g).all()
    assert abs(float(g[0, 0]) - HUBBARD[8]) < 1e-6


def test_translation_rotation_invariance():
    Z = torch.tensor([8, 1, 1])
    R = torch.tensor([
        [0.0, 0.0, 0.117],
        [0.0, 0.757, -0.469],
        [0.0, -0.757, -0.469],
    ], dtype=torch.float64)
    q = torch.tensor([-0.834, 0.417, 0.417], dtype=torch.float64)

    E0 = coulomb_energy(q, R, Z, HUBBARD)
    masses = masses_from_Z(Z)
    mu0 = dipole(q, R, masses)

    R_t = R + torch.tensor([1.5, -2.0, 3.7], dtype=torch.float64)
    E_t = coulomb_energy(q, R_t, Z, HUBBARD)
    mu_t = dipole(q, R_t, masses)
    assert abs(float(E_t - E0)) < 1e-10
    assert torch.allclose(mu_t, mu0, atol=1e-10)

    theta = 0.7
    Rot = torch.tensor([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta), np.cos(theta), 0.0],
        [0.0, 0.0, 1.0],
    ], dtype=torch.float64)
    R_r = R @ Rot.T
    E_r = coulomb_energy(q, R_r, Z, HUBBARD)
    mu_r = dipole(q, R_r, masses)
    assert abs(float(E_r - E0)) < 1e-10
    assert torch.allclose(mu_r, Rot @ mu0, atol=1e-10)
