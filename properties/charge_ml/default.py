"""Valeurs par défaut centralisées du module charge_ml."""


DEFAULT_DESC: dict = dict(
    kind="soap",
    r_cut=8.0,
    n_max = 6,
    l_max = 5,
)

DEFAULT_HUBBARD: dict[int, float] = {
    1: 0.4196,
    6: 0.3647,
    7: 0.4309,
    8: 0.4954,
}

DEFAULT_CHARGES: dict[int, float] = {
    1: 1.0,
    6: 4.0,
    7: 5.0,
    8: 6.0,
}

DEFAULT_EPOCH: int = 200

DEFAULT_BATCH: int = 8
