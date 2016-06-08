"""
Utilities for structure analysis.
"""
from . import measure
from .atoms import Selection

__all__ = ['hydrogen_bonds', 'HydrogenBond']


class HydrogenBond(object):
    """
    Represents hydrogen bond.

    @ivar: Hydrogen donor atom
    @ivar: Hydrogen atom
    @ivar: Hydrogen acceptor atom
    """
    def __init__(self, donor, hydrogen, acceptor):
        self.donor = donor
        self.hydrogen = hydrogen
        self.acceptor = acceptor

    def __repr__(self):
        return '<%s: %s--%s..%s>' % (type(self).__name__, self.donor, self.hydrogen, self.acceptor)

    def __eq__(self, other):
        return type(self) == type(other) and self.donor == other.donor and self.hydrogen == other.hydrogen and \
            self.acceptor == other.acceptor

    def __ne__(self, other):
        return not self.__eq__(other)


def _get_bonds(donor, acceptor, donor_hydrogens, angle):
    # Utility function which finds hydrogen bonds between donor atom and acceptor atom
    hydrogens = (a for a in donor.bonded if a in donor_hydrogens)
    for hydrogen in hydrogens:
        # Check the angle. If it's big enough, then it is a hydrogen bond
        if measure.angle(donor, hydrogen, acceptor) >= angle:
            yield HydrogenBond(donor, hydrogen, acceptor)


def hydrogen_bonds(donors, acceptors=None, distance=3.0, angle=135):
    """
    Returns iterator of hydrogen bonds between the selections.

    @param donors: Hydrogen donors selection
    @type donors: Selection
    @param acceptors: Hydrogen acceptors selection
    @type donors: Selection or None
    @param distance: Maximal distance between donor and acceptor
    @type distance: Non-negative number
    @param angle: Minimal angle in degrees between donor, hydrogen and acceptor
    @type angle: Number between 0 and 180
    @rtype: Generator of HydrogenBond objects
    """
    assert isinstance(donors, Selection)
    assert acceptors is None or isinstance(acceptors, Selection)
    assert distance >= 0
    assert 0 <= angle <= 180

    # Remove hydrogen atoms from selection. This can be done safely, hydrogens are never donors.
    donor_heavy = Selection('(%s) and noh' % donors.selection, donors.molecule, donors.frame)
    # Create selections of hydrogens for donor molecule. This will be used to find the hydrogen involved in the bond.
    donor_hydrogens = Selection('hydrogen', donors.molecule, donors.frame)
    if acceptors is None:
        # Acceptor is same as donor, just copy
        acceptor_heavy = donor_heavy
        acceptor_hydrogens = donor_hydrogens
    else:
        # Acceptor is not the same as donor. Make same selections as for donor.
        acceptor_heavy = Selection('(%s) and noh' % acceptors.selection, acceptors.molecule, acceptors.frame)
        acceptor_hydrogens = Selection('hydrogen', acceptors.molecule, acceptors.frame)

    for donor, acceptor in donor_heavy.contacts(acceptor_heavy, distance):
        for hbond in _get_bonds(donor, acceptor, donor_hydrogens, angle):
            yield hbond
        # If acceptors and donors share atoms, contacts return pair only once.
        # Check if donor and acceptors can have opposite roles.
        if donor in acceptor_heavy and acceptor in donor_heavy:
            # Donor can be acceptor and acceptor can be donor, check the hydrogen bonds.
            for hbond in _get_bonds(acceptor, donor, acceptor_hydrogens, angle):
                yield hbond
