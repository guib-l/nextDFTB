"""Affichage textuel de la progression d'un entraînement."""


def print_epoch(epoch: int, metrics: dict[str, float]) -> None:
    """Affiche une ligne récapitulative pour une époque."""
    items = "  ".join(f"{k}={v:.4e}" for k, v in metrics.items())
    print(f"epoch {epoch:4d} | {items}", flush=True)


def print_history(history: dict[str, list[float]]) -> None:
    """Affiche un tableau récapitulatif d'un historique d'entraînement."""
    if not history:
        return
    keys = list(history.keys())
    n = len(history[keys[0]])
    header = "epoch  " + "  ".join(f"{k:>12s}" for k in keys)
    print(header)
    for i in range(n):
        row = f"{i:5d}  " + "  ".join(f"{history[k][i]:12.4e}" for k in keys)
        print(row)
