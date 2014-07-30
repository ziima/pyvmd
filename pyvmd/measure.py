"""
Utilities for simple strucural analysis.
"""
import math

from numpy import cross
from numpy.linalg import norm

from .atoms import Atom

__all__ = ['angle', 'dihedral', 'distance']


def coords_distance(a, b):
    """
    Returns distance between two coordinates.
    """
    assert a.shape == (3, )
    assert b.shape == (3, )
    return norm(a - b)


def distance(a, b):
    """
    Returns distance between two atoms.
    """
    assert isinstance(a, Atom)
    assert isinstance(b, Atom)
    return coords_distance(a.coords, b.coords)


def coords_angle(a, b, c):
    """
    Returns angle between three coordinates a--b--c in degrees.
    """
    assert a.shape == (3, )
    assert b.shape == (3, )
    assert c.shape == (3, )
    # Get vectors b-->a and b-->c
    vec_1 = a - b
    vec_2 = c - b
    # Compute angle between the two vectors
    cross_prod = cross(vec_1, vec_2)
    sine = math.sqrt(cross_prod.dot(cross_prod))
    cosine = vec_1.dot(vec_2)
    return math.degrees(math.atan2(sine, cosine))


def angle(a, b, c):
    """
    Returns angle between three atoms a--b--c in degrees.
    """
    assert isinstance(a, Atom)
    assert isinstance(b, Atom)
    assert isinstance(c, Atom)
    return coords_angle(a.coords, b.coords, c.coords)


def coords_dihedral(a, b, c, d):
    """
    Returns dihedral angle of four coordinates a--b--c--d in degrees.
    """
    assert a.shape == (3, )
    assert b.shape == (3, )
    assert c.shape == (3, )
    assert d.shape == (3, )
    # Get vectors a-->b, b-->c and c-->d
    vec_1 = b - a
    vec_2 = c - b
    vec_3 = d - c
    # Compute the dihedral of the three vectors
    norm_1 = cross(vec_1, vec_2)
    norm_2 = cross(vec_2, vec_3)
    sine = norm_1.dot(vec_3) * norm(vec_2)
    cosine = norm_1.dot(norm_2)
    return math.degrees(math.atan2(sine, cosine))


def dihedral(a, b, c, d):
    """
    Returns dihedral or improper dihedral angle of four atoms in degrees.

    Atoms are considered in structure a--b--c--d for dihedral angle or

    a--b--c
       |
       d

    for improper dihedral angle.
    """
    assert isinstance(a, Atom)
    assert isinstance(b, Atom)
    assert isinstance(c, Atom)
    assert isinstance(d, Atom)
    return coords_dihedral(a.coords, b.coords, c.coords, d.coords)
