import numpy as np
import torch

from charge_ml.electrostatic import HUBBARD_DEFAULT
from charge_ml.features import build_descriptor, featurize
from charge_ml.model import ChargeModel, ElectronegativityNet
from charge_ml.predict import ChargePredictor


ELEMENTS = ["H", "C", "N", "O"]
DESC_KW = dict(
    r_cut=5.0,
    g2_params=[[1.0, 1.0], [1.0, 2.0]],
    g4_params=[[1.0, 1.0, 1.0]],
)


def _water():
    Z = np.array([8, 1, 1])
    R = np.array([
        [0.0, 0.0, 0.117],
        [0.0, 0.757, -0.469],
        [0.0, -0.757, -0.469],
    ])
    return Z, R


def _make_predictor():
    desc = build_descriptor("acsf", ELEMENTS, **DESC_KW)
    Z, R = _water()
    d_in = featurize(desc, Z, R).shape[1]
    chi_net = ElectronegativityNet(d_in=d_in, elements=ELEMENTS, hidden=8)
    model = ChargeModel(chi_net, hardness=HUBBARD_DEFAULT, hubbard=HUBBARD_DEFAULT)
    return ChargePredictor(model, desc, ELEMENTS)


def test_conservation():
    pred = _make_predictor()
    Z, R = _water()
    for Q in (0.0, 1.0, -1.0):
        q = pred.predict(Z, R, Q_tot=Q)
        assert abs(q.sum() - Q) < 1e-5


def test_invariance_translation_rotation():
    pred = _make_predictor()
    Z, R = _water()
    q0 = pred.predict(Z, R, Q_tot=0.0)

    R_t = R + np.array([1.5, -2.0, 3.7])
    q_t = pred.predict(Z, R_t, Q_tot=0.0)
    assert np.allclose(q0, q_t, atol=1e-5)

    theta = 0.7
    Rot = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta), np.cos(theta), 0.0],
        [0.0, 0.0, 1.0],
    ])
    R_r = R @ Rot.T
    q_r = pred.predict(Z, R_r, Q_tot=0.0)
    assert np.allclose(q0, q_r, atol=1e-5)


def test_permutation_equivariance():
    pred = _make_predictor()
    Z, R = _water()
    q0 = pred.predict(Z, R, Q_tot=0.0)

    perm = np.array([2, 0, 1])
    q_p = pred.predict(Z[perm], R[perm], Q_tot=0.0)
    assert np.allclose(q0[perm], q_p, atol=1e-5)


def test_gradient_flows_through_qeq():
    desc = build_descriptor("acsf", ELEMENTS, **DESC_KW)
    Z_np, R_np = _water()
    X = featurize(desc, Z_np, R_np)
    d_in = X.shape[1]
    chi_net = ElectronegativityNet(d_in=d_in, elements=ELEMENTS, hidden=8)
    model = ChargeModel(chi_net, hardness=HUBBARD_DEFAULT, hubbard=HUBBARD_DEFAULT)

    X_t = torch.from_numpy(np.asarray(X, dtype=np.float32))
    R_t = torch.from_numpy(R_np.astype(np.float32))
    Z_t = torch.from_numpy(Z_np)
    q_ref = torch.tensor([-0.834, 0.417, 0.417], dtype=torch.float32)

    q = model(X_t, Z_t, R_t, 0.0)
    loss = ((q - q_ref) ** 2).mean()
    loss.backward()

    grads = [p.grad for p in model.parameters() if p.grad is not None]
    assert len(grads) > 0
    assert any(g.abs().sum() > 0 for g in grads)
