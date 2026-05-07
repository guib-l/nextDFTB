"""Métriques utilisées pour l'entraînement et l'évaluation."""

import torch
from torch import Tensor


def metric_mae(pred: Tensor, ref: Tensor) -> Tensor:
    """Mean Absolute Error."""
    return torch.mean(torch.abs(pred - ref))


def metric_mse(pred: Tensor, ref: Tensor) -> Tensor:
    """Mean Squared Error."""
    return torch.mean((pred - ref) ** 2)
