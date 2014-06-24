"""
Utilities for simple strucural analysis.
"""
import math
from numpy import cross, ndarray
from numpy.linalg import norm

from .objects import Atom

__all__ = ['angle', 'dihedral', 'distance', 'improper']


def _vector_angle(a, b):
    """
    Returns angle (in degrees) between two vectors.
    """
    assert isinstance(a, ndarray)
    assert isinstance(b, ndarray)
    cross_prod = cross(a, b)
    sine = math.sqrt(cross_prod.dot(cross_prod))
    cosine = a.dot(b)
    return math.degrees(math.atan2(sine, cosine))


def distance(a, b):
    """
    Returns distance between two atoms.
    """
    assert isinstance(a, Atom)
    assert isinstance(b, Atom)
    assert a != b
    return norm(a.coords - b.coords)


def angle(a, b, c):
    """
    Returns angle between three atoms a--b--c in degrees.
    """
    assert isinstance(a, Atom)
    assert isinstance(b, Atom)
    assert isinstance(c, Atom)
    assert a != b and a != c and b != c
    coords_1 = a.coords
    coords_2 = b.coords
    coords_3 = c.coords
    vec_1 = coords_1 - coords_2
    vec_2 = coords_3 - coords_2
    return _vector_angle(vec_1, vec_2)


def dihedral(a, b, c, d):
    """
    Returns dihedral angle of four atoms a--b--c--d in degrees.
    """
    assert isinstance(a, Atom)
    assert isinstance(b, Atom)
    assert isinstance(c, Atom)
    assert isinstance(d, Atom)
    assert a != b and a != c and a != d and b != c and b != d and c != d
    coords_1 = a.coords
    coords_2 = b.coords
    coords_3 = c.coords
    coords_4 = d.coords
    vec_1 = coords_2 - coords_1
    vec_2 = coords_3 - coords_2
    vec_3 = coords_4 - coords_3

    norm_1 = cross(vec_1, vec_2)
    norm_2 = cross(vec_2, vec_3)
    sine = norm_1.dot(vec_3) * norm(vec_2)
    cosine = norm_1.dot(norm_2)
    return math.degrees(math.atan2(sine, cosine))


def improper(a, b, c, d):
    """
    Returns improper dihedral angle of four atoms in degrees.

    a--b--c
       |
       d
    """
    return dihedral(a, b, c, d)
