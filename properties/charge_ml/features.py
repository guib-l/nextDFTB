import numpy as np
from ase import Atoms
from dscribe.descriptors import ACSF, SOAP


def build_descriptor(kind: str, species: list[str], **kwargs):
    if kind == "acsf":
        return ACSF(species=species, **kwargs)
    if kind == "soap":
        return SOAP(species=species, **kwargs)
    raise ValueError(kind)


def featurize(desc, Z: np.ndarray, R: np.ndarray) -> np.ndarray:
    return desc.create(Atoms(numbers=Z, positions=R))
