"""
Tests for drawing methods of molecule representations.
"""
from VMD import molecule as _molecule, molrep as _molrep

from pyvmd.representations import (DRAW_BEADS, DRAW_BONDS, DRAW_CARTOON, DRAW_CPK, DRAW_DOTTED, DRAW_DYNAMIC_BONDS,
                                   DRAW_FIELD_LINES, DRAW_HBONDS, DRAW_ISOSURFACE, DRAW_LICORICE, DRAW_LINES, DRAW_MSMS,
                                   DRAW_NEW_CARTOON, DRAW_NEW_RIBBONS, DRAW_ORBITAL, DRAW_PAPER_CHAIN, DRAW_POINTS,
                                   DRAW_POLYHEDRA, DRAW_QUICKSURF, DRAW_RIBBONS, DRAW_SOLVENT, DRAW_SURF, DRAW_TRACE,
                                   DRAW_TUBE, DRAW_TWISTER, DRAW_VDW, DRAW_VOLUME_SLICE, Representation)

from .utils import data, PyvmdTestCase


class TestDrawingMethods(PyvmdTestCase):
    """
    Test drawing methods are defined correctly.
    """
    def setUp(self):
        self.molid = _molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        self.rep = Representation('rep0')
        self.style = self.rep.style

    def test_lines(self):
        self.style.method = DRAW_LINES
        self.assertEqual(self.style.method, DRAW_LINES)
        self.assertEqual(self.style.get_parameters(), {'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Lines')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_LINES)
        self.assertEqual(self.style.get_parameters(), {'size': 3})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Lines 3.0')

    def test_bonds(self):
        self.style.method = DRAW_BONDS
        self.assertEqual(self.style.method, DRAW_BONDS)
        self.assertEqual(self.style.get_parameters(), {'size': 0.3, 'resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Bonds')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_BONDS)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Bonds 3.0 10')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_BONDS)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 20})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Bonds 3.0 20.0')

    def test_dynamic_bonds(self):
        self.style.method = DRAW_DYNAMIC_BONDS
        self.assertEqual(self.style.method, DRAW_DYNAMIC_BONDS)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 3, 'size': 0.3, 'resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'DynamicBonds')

        self.style.set_parameters(size=2.)
        self.assertEqual(self.style.method, DRAW_DYNAMIC_BONDS)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 3, 'size': 2, 'resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'DynamicBonds 3.0 2.0 10')

        self.style.set_parameters(cutoff=7.)
        self.assertEqual(self.style.method, DRAW_DYNAMIC_BONDS)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 7, 'size': 2, 'resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'DynamicBonds 7.0 2.0 10.0')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_DYNAMIC_BONDS)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 7, 'size': 2, 'resolution': 20})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'DynamicBonds 7.0 2.0 20.0')

    def test_hbonds(self):
        self.style.method = DRAW_HBONDS
        self.assertEqual(self.style.method, DRAW_HBONDS)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 3, 'angle_cutoff': 20, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'HBonds')

        self.style.set_parameters(cutoff=3.5)
        self.assertEqual(self.style.method, DRAW_HBONDS)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 3.5, 'angle_cutoff': 20, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'HBonds 3.5 20.0 1.0')

        self.style.set_parameters(angle_cutoff=30.)
        self.assertEqual(self.style.method, DRAW_HBONDS)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 3.5, 'angle_cutoff': 30, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'HBonds 3.5 30.0 1.0')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_HBONDS)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 3.5, 'angle_cutoff': 30, 'size': 3})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'HBonds 3.5 30.0 3.0')

    def test_points(self):
        self.style.method = DRAW_POINTS
        self.assertEqual(self.style.method, DRAW_POINTS)
        self.assertEqual(self.style.get_parameters(), {'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Points')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_POINTS)
        self.assertEqual(self.style.get_parameters(), {'size': 3})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Points 3.0')

    def test_vdw(self):
        self.style.method = DRAW_VDW
        self.assertEqual(self.style.method, DRAW_VDW)
        self.assertEqual(self.style.get_parameters(), {'size': 1, 'resolution': 12})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'VDW')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_VDW)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 12})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'VDW 3.0 12')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_VDW)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 20})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'VDW 3.0 20.0')

    def test_cpk(self):
        self.style.method = DRAW_CPK
        self.assertEqual(self.style.method, DRAW_CPK)
        self.assertEqual(self.style.get_parameters(),
                         {'size': 1, 'bond_radius': 0.3, 'resolution': 12, 'bond_resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'CPK')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_CPK)
        self.assertEqual(self.style.get_parameters(),
                         {'size': 3, 'bond_radius': 0.3, 'resolution': 12, 'bond_resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'CPK 3.0 0.3 12 10')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_CPK)
        self.assertEqual(self.style.get_parameters(),
                         {'size': 3, 'bond_radius': 0.3, 'resolution': 20, 'bond_resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'CPK 3.0 0.3 20.0 10.0')

        self.style.set_parameters(bond_radius=0.7)
        self.assertEqual(self.style.method, DRAW_CPK)
        self.assertEqual(self.style.get_parameters(),
                         {'size': 3, 'bond_radius': 0.7, 'resolution': 20, 'bond_resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'CPK 3.0 0.7 20.0 10.0')

        self.style.set_parameters(bond_resolution=30.)
        self.assertEqual(self.style.method, DRAW_CPK)
        self.assertEqual(self.style.get_parameters(),
                         {'size': 3, 'bond_radius': 0.7, 'resolution': 20, 'bond_resolution': 30})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'CPK 3.0 0.7 20.0 30.0')

    def test_licorice(self):
        self.style.method = DRAW_LICORICE
        self.assertEqual(self.style.method, DRAW_LICORICE)
        self.assertEqual(self.style.get_parameters(), {'size': 0.3, 'resolution': 10, 'bond_resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Licorice')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_LICORICE)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 10, 'bond_resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Licorice 3.0 10 10')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_LICORICE)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 20, 'bond_resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Licorice 3.0 20.0 10.0')

        self.style.set_parameters(bond_resolution=30.)
        self.assertEqual(self.style.method, DRAW_LICORICE)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 20, 'bond_resolution': 30})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Licorice 3.0 20.0 30.0')

    def test_polyhedra(self):
        self.style.method = DRAW_POLYHEDRA
        self.assertEqual(self.style.method, DRAW_POLYHEDRA)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 3.})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Polyhedra')

        self.style.set_parameters(cutoff=5.)
        self.assertEqual(self.style.method, DRAW_POLYHEDRA)
        self.assertEqual(self.style.get_parameters(), {'cutoff': 5.})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Polyhedra 5.0')

    def test_trace(self):
        self.style.method = DRAW_TRACE
        self.assertEqual(self.style.method, DRAW_TRACE)
        self.assertEqual(self.style.get_parameters(), {'size': 0.3, 'resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Trace')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_TRACE)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Trace 3.0 10')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_TRACE)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 20})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Trace 3.0 20.0')

    def test_tube(self):
        self.style.method = DRAW_TUBE
        self.assertEqual(self.style.method, DRAW_TUBE)
        self.assertEqual(self.style.get_parameters(), {'size': 0.3, 'resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Tube')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_TUBE)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Tube 3.0 10')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_TUBE)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 20})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Tube 3.0 20.0')

    def test_ribbons(self):
        self.style.method = DRAW_RIBBONS
        self.assertEqual(self.style.method, DRAW_RIBBONS)
        self.assertEqual(self.style.get_parameters(), {'margin_size': 0.3, 'resolution': 10, 'size': 2})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Ribbons')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_RIBBONS)
        self.assertEqual(self.style.get_parameters(), {'margin_size': 0.3, 'resolution': 10, 'size': 3})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Ribbons 0.3 10 3.0')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_RIBBONS)
        self.assertEqual(self.style.get_parameters(), {'margin_size': 0.3, 'resolution': 20, 'size': 3})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Ribbons 0.3 20.0 3.0')

        self.style.set_parameters(margin_size=0.7)
        self.assertEqual(self.style.method, DRAW_RIBBONS)
        self.assertEqual(self.style.get_parameters(), {'margin_size': 0.7, 'resolution': 20, 'size': 3})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Ribbons 0.7 20.0 3.0')

    def test_new_ribbons(self):
        self.style.method = DRAW_NEW_RIBBONS
        self.assertEqual(self.style.method, DRAW_NEW_RIBBONS)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.3, 'resolution': 10, 'width': 3, 'spline': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewRibbons')

        self.style.set_parameters(width=5.)
        self.assertEqual(self.style.method, DRAW_NEW_RIBBONS)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.3, 'resolution': 10, 'width': 5, 'spline': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewRibbons 0.3 10 5.0 0')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_NEW_RIBBONS)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.3, 'resolution': 20, 'width': 5, 'spline': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewRibbons 0.3 20.0 5.0 0.0')

        self.style.set_parameters(spline=1)
        self.assertEqual(self.style.method, DRAW_NEW_RIBBONS)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.3, 'resolution': 20, 'width': 5, 'spline': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewRibbons 0.3 20.0 5.0 1')

        self.style.set_parameters(thickness=0.7)
        self.assertEqual(self.style.method, DRAW_NEW_RIBBONS)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.7, 'resolution': 20, 'width': 5, 'spline': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewRibbons 0.7 20.0 5.0 1.0')

    def test_cartoon(self):
        self.style.method = DRAW_CARTOON
        self.assertEqual(self.style.method, DRAW_CARTOON)
        self.assertEqual(self.style.get_parameters(), {'helix_size': 2.1, 'resolution': 12, 'sheet_size': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Cartoon')

        self.style.set_parameters(helix_size=3.)
        self.assertEqual(self.style.method, DRAW_CARTOON)
        self.assertEqual(self.style.get_parameters(), {'helix_size': 3, 'resolution': 12, 'sheet_size': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Cartoon 3.0 12 5.0')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_CARTOON)
        self.assertEqual(self.style.get_parameters(), {'helix_size': 3, 'resolution': 20, 'sheet_size': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Cartoon 3.0 20.0 5.0')

        self.style.set_parameters(sheet_size=0.7)
        self.assertEqual(self.style.method, DRAW_CARTOON)
        self.assertEqual(self.style.get_parameters(), {'helix_size': 3, 'resolution': 20, 'sheet_size': 0.7})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Cartoon 3.0 20.0 0.7')

    def test_new_cartoon(self):
        self.style.method = DRAW_NEW_CARTOON
        self.assertEqual(self.style.method, DRAW_NEW_CARTOON)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.3, 'resolution': 10, 'width': 4.5, 'spline': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewCartoon')

        self.style.set_parameters(width=5.)
        self.assertEqual(self.style.method, DRAW_NEW_CARTOON)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.3, 'resolution': 10, 'width': 5, 'spline': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewCartoon 0.3 10 5.0 0')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_NEW_CARTOON)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.3, 'resolution': 20, 'width': 5, 'spline': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewCartoon 0.3 20.0 5.0 0.0')

        self.style.set_parameters(spline=1)
        self.assertEqual(self.style.method, DRAW_NEW_CARTOON)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.3, 'resolution': 20, 'width': 5, 'spline': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewCartoon 0.3 20.0 5.0 1')

        self.style.set_parameters(thickness=0.7)
        self.assertEqual(self.style.method, DRAW_NEW_CARTOON)
        self.assertEqual(self.style.get_parameters(), {'thickness': 0.7, 'resolution': 20, 'width': 5, 'spline': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'NewCartoon 0.7 20.0 5.0 1.0')

    def test_paper_chain(self):
        self.style.method = DRAW_PAPER_CHAIN
        self.assertEqual(self.style.method, DRAW_PAPER_CHAIN)
        self.assertEqual(self.style.get_parameters(), {'height': 1, 'max_ring_size': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'PaperChain')

        self.style.set_parameters(height=3.)
        self.assertEqual(self.style.method, DRAW_PAPER_CHAIN)
        self.assertEqual(self.style.get_parameters(), {'height': 3, 'max_ring_size': 10})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'PaperChain 3.0 10')

        self.style.set_parameters(max_ring_size=12.)
        self.assertEqual(self.style.method, DRAW_PAPER_CHAIN)
        self.assertEqual(self.style.get_parameters(), {'height': 3, 'max_ring_size': 12})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'PaperChain 3.0 12.0')

    def test_twister(self):
        self.style.method = DRAW_TWISTER
        self.assertEqual(self.style.method, DRAW_TWISTER)
        self.assertEqual(self.style.get_parameters(),
                         {'start': 1, 'hide_shared_links': 0, 'steps': 10, 'width': 0.3, 'height': 0.05,
                          'max_ring_size': 10, 'linking_distance': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Twister')

        self.style.set_parameters(width=0.5)
        self.assertEqual(self.style.method, DRAW_TWISTER)
        self.assertEqual(self.style.get_parameters(),
                         {'start': 1, 'hide_shared_links': 0, 'steps': 10, 'width': 0.5, 'height': 0.05,
                          'max_ring_size': 10, 'linking_distance': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Twister 1 0 10 0.5 0.05 10 5')

        self.style.set_parameters(height=0.1)
        self.assertEqual(self.style.method, DRAW_TWISTER)
        self.assertEqual(self.style.get_parameters(),
                         {'start': 1, 'hide_shared_links': 0, 'steps': 10, 'width': 0.5, 'height': 0.1,
                          'max_ring_size': 10, 'linking_distance': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Twister 1.0 0.0 10.0 0.5 0.1 10.0 5.0')

        self.style.set_parameters(start=0.)
        self.assertEqual(self.style.method, DRAW_TWISTER)
        self.assertEqual(self.style.get_parameters(),
                         {'start': 0, 'hide_shared_links': 0, 'steps': 10, 'width': 0.5, 'height': 0.1,
                          'max_ring_size': 10, 'linking_distance': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Twister 0.0 0.0 10.0 0.5 0.1 10.0 5.0')

        self.style.set_parameters(hide_shared_links=1.)
        self.assertEqual(self.style.method, DRAW_TWISTER)
        self.assertEqual(self.style.get_parameters(),
                         {'start': 0, 'hide_shared_links': 1, 'steps': 10, 'width': 0.5, 'height': 0.1,
                          'max_ring_size': 10, 'linking_distance': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Twister 0.0 1.0 10.0 0.5 0.1 10.0 5.0')

        self.style.set_parameters(steps=12.)
        self.assertEqual(self.style.method, DRAW_TWISTER)
        self.assertEqual(self.style.get_parameters(),
                         {'start': 0, 'hide_shared_links': 1, 'steps': 12, 'width': 0.5, 'height': 0.1,
                          'max_ring_size': 10, 'linking_distance': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Twister 0.0 1.0 12.0 0.5 0.1 10.0 5.0')

        self.style.set_parameters(max_ring_size=17.)
        self.assertEqual(self.style.method, DRAW_TWISTER)
        self.assertEqual(self.style.get_parameters(),
                         {'start': 0, 'hide_shared_links': 1, 'steps': 12, 'width': 0.5, 'height': 0.1,
                          'max_ring_size': 17, 'linking_distance': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Twister 0.0 1.0 12.0 0.5 0.1 17.0 5.0')

        self.style.set_parameters(linking_distance=7.)
        self.assertEqual(self.style.method, DRAW_TWISTER)
        self.assertEqual(self.style.get_parameters(),
                         {'start': 0, 'hide_shared_links': 1, 'steps': 12, 'width': 0.5, 'height': 0.1,
                          'max_ring_size': 17, 'linking_distance': 7})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Twister 0.0 1.0 12.0 0.5 0.1 17.0 7.0')

    def test_quicksurf(self):
        self.style.method = DRAW_QUICKSURF
        self.assertEqual(self.style.method, DRAW_QUICKSURF)
        self.assertEqual(self.style.get_parameters(),
                         {'sphere_scale': 1, 'density_isovalue': 0.5, 'grid_spacing': 1, 'resolution': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'QuickSurf')

        self.style.set_parameters(density_isovalue=0.7)
        self.assertEqual(self.style.method, DRAW_QUICKSURF)
        self.assertEqual(self.style.get_parameters(),
                         {'sphere_scale': 1, 'density_isovalue': 0.7, 'grid_spacing': 1, 'resolution': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'QuickSurf 1.0 0.7 1.0 0')

        self.style.set_parameters(sphere_scale=1.3)
        self.assertEqual(self.style.method, DRAW_QUICKSURF)
        self.assertEqual(self.style.get_parameters(),
                         {'sphere_scale': 1.3, 'density_isovalue': 0.7, 'grid_spacing': 1, 'resolution': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'QuickSurf 1.3 0.7 1.0 0.0')

        self.style.set_parameters(grid_spacing=1.6)
        self.assertEqual(self.style.method, DRAW_QUICKSURF)
        self.assertEqual(self.style.get_parameters(),
                         {'sphere_scale': 1.3, 'density_isovalue': 0.7, 'grid_spacing': 1.6, 'resolution': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'QuickSurf 1.3 0.7 1.6 0.0')

        self.style.set_parameters(resolution=2.3)
        self.assertEqual(self.style.method, DRAW_QUICKSURF)
        self.assertEqual(self.style.get_parameters(),
                         {'sphere_scale': 1.3, 'density_isovalue': 0.7, 'grid_spacing': 1.6, 'resolution': 2.3})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'QuickSurf 1.3 0.7 1.6 2.3')

    def test_msms(self):
        self.style.method = DRAW_MSMS
        self.assertEqual(self.style.method, DRAW_MSMS)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1.5, 'density': 1.5, 'atoms': 0, 'method': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'MSMS')

        self.style.set_parameters(density=1.7)
        self.assertEqual(self.style.method, DRAW_MSMS)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1.5, 'density': 1.7, 'atoms': 0, 'method': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'MSMS 1.5 1.7 0 0')

        self.style.set_parameters(probe_radius=1.2)
        self.assertEqual(self.style.method, DRAW_MSMS)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1.2, 'density': 1.7, 'atoms': 0, 'method': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'MSMS 1.2 1.7 0.0 0.0')

        self.style.set_parameters(atoms=1.)
        self.assertEqual(self.style.method, DRAW_MSMS)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1.2, 'density': 1.7, 'atoms': 1, 'method': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'MSMS 1.2 1.7 1.0 0.0')

        self.style.set_parameters(method=2.)
        self.assertEqual(self.style.method, DRAW_MSMS)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1.2, 'density': 1.7, 'atoms': 1, 'method': 2})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'MSMS 1.2 1.7 1.0 2.0')

    def test_surf(self):
        self.style.method = DRAW_SURF
        self.assertEqual(self.style.method, DRAW_SURF)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1.4, 'method': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Surf')

        self.style.set_parameters(probe_radius=1.2)
        self.assertEqual(self.style.method, DRAW_SURF)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1.2, 'method': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Surf 1.2 0')

        self.style.set_parameters(method=2.)
        self.assertEqual(self.style.method, DRAW_SURF)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1.2, 'method': 2})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Surf 1.2 2.0')

    def test_volume_slice(self):
        self.style.method = DRAW_VOLUME_SLICE
        self.assertEqual(self.style.method, DRAW_VOLUME_SLICE)
        self.assertEqual(self.style.get_parameters(), {'slice': 0.5, 'volume_id': 0, 'axis': 0, 'quality': 2})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'VolumeSlice')

        self.style.set_parameters(slice=1.2)
        self.assertEqual(self.style.method, DRAW_VOLUME_SLICE)
        self.assertEqual(self.style.get_parameters(), {'slice': 1.2, 'volume_id': 0, 'axis': 0, 'quality': 2})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'VolumeSlice 1.2 0 0 2')

        self.style.set_parameters(volume_id=3.)
        self.assertEqual(self.style.method, DRAW_VOLUME_SLICE)
        self.assertEqual(self.style.get_parameters(), {'slice': 1.2, 'volume_id': 3, 'axis': 0, 'quality': 2})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'VolumeSlice 1.2 3.0 0.0 2.0')

        self.style.set_parameters(axis=2.)
        self.assertEqual(self.style.method, DRAW_VOLUME_SLICE)
        self.assertEqual(self.style.get_parameters(), {'slice': 1.2, 'volume_id': 3, 'axis': 2, 'quality': 2})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'VolumeSlice 1.2 3.0 2.0 2.0')

        self.style.set_parameters(quality=5.)
        self.assertEqual(self.style.method, DRAW_VOLUME_SLICE)
        self.assertEqual(self.style.get_parameters(), {'slice': 1.2, 'volume_id': 3, 'axis': 2, 'quality': 5})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'VolumeSlice 1.2 3.0 2.0 5.0')

    def test_isosurface(self):
        self.style.method = DRAW_ISOSURFACE
        self.assertEqual(self.style.method, DRAW_ISOSURFACE)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.5, 'volume_id': 0, 'display': 2, 'method': 2, 'step': 1, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Isosurface')

        self.style.set_parameters(isovalue=0.7)
        self.assertEqual(self.style.method, DRAW_ISOSURFACE)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.7, 'volume_id': 0, 'display': 2, 'method': 2, 'step': 1, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Isosurface 0.7 0 2 2 1 1')

        self.style.set_parameters(volume_id=7.)
        self.assertEqual(self.style.method, DRAW_ISOSURFACE)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.7, 'volume_id': 7, 'display': 2, 'method': 2, 'step': 1, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Isosurface 0.7 7.0 2.0 2.0 1.0 1.0')

        self.style.set_parameters(display=1.)
        self.assertEqual(self.style.method, DRAW_ISOSURFACE)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.7, 'volume_id': 7, 'display': 1, 'method': 2, 'step': 1, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Isosurface 0.7 7.0 1.0 2.0 1.0 1.0')

        self.style.set_parameters(method=4.)
        self.assertEqual(self.style.method, DRAW_ISOSURFACE)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.7, 'volume_id': 7, 'display': 1, 'method': 4, 'step': 1, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Isosurface 0.7 7.0 1.0 4.0 1.0 1.0')

        self.style.set_parameters(step=6.)
        self.assertEqual(self.style.method, DRAW_ISOSURFACE)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.7, 'volume_id': 7, 'display': 1, 'method': 4, 'step': 6, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Isosurface 0.7 7.0 1.0 4.0 6.0 1.0')

        self.style.set_parameters(size=12.)
        self.assertEqual(self.style.method, DRAW_ISOSURFACE)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.7, 'volume_id': 7, 'display': 1, 'method': 4, 'step': 6, 'size': 12})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Isosurface 0.7 7.0 1.0 4.0 6.0 12.0')

    def test_field_lines(self):
        self.style.method = DRAW_FIELD_LINES
        self.assertEqual(self.style.method, DRAW_FIELD_LINES)
        self.assertEqual(self.style.get_parameters(),
                         {'volume_id': 0, 'gradient': 1.8, 'min_length': 10, 'max_length': 50, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'FieldLines')

        self.style.set_parameters(gradient=0.7)
        self.assertEqual(self.style.method, DRAW_FIELD_LINES)
        self.assertEqual(self.style.get_parameters(),
                         {'volume_id': 0, 'gradient': 0.7, 'min_length': 10, 'max_length': 50, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'FieldLines 0 0.7 10 50 1')

        self.style.set_parameters(volume_id=7.)
        self.assertEqual(self.style.method, DRAW_FIELD_LINES)
        self.assertEqual(self.style.get_parameters(),
                         {'volume_id': 7, 'gradient': 0.7, 'min_length': 10, 'max_length': 50, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'FieldLines 7.0 0.7 10.0 50.0 1.0')

        self.style.set_parameters(min_length=12.)
        self.assertEqual(self.style.method, DRAW_FIELD_LINES)
        self.assertEqual(self.style.get_parameters(),
                         {'volume_id': 7, 'gradient': 0.7, 'min_length': 12, 'max_length': 50, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'FieldLines 7.0 0.7 12.0 50.0 1.0')

        self.style.set_parameters(max_length=42.)
        self.assertEqual(self.style.method, DRAW_FIELD_LINES)
        self.assertEqual(self.style.get_parameters(),
                         {'volume_id': 7, 'gradient': 0.7, 'min_length': 12, 'max_length': 42, 'size': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'FieldLines 7.0 0.7 12.0 42.0 1.0')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_FIELD_LINES)
        self.assertEqual(self.style.get_parameters(),
                         {'volume_id': 7, 'gradient': 0.7, 'min_length': 12, 'max_length': 42, 'size': 3})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'FieldLines 7.0 0.7 12.0 42.0 3.0')

    def test_orbital(self):
        self.style.method = DRAW_ORBITAL
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.05, 'orbital_id': 0, 'display': 0, 'method': 0, 'grid_spacing': 0.075,
                          'size': 1, 'wavefunction': 0, 'spin': 0, 'excitation': 0, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital')

        self.style.set_parameters(grid_spacing=0.1)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.05, 'orbital_id': 0, 'display': 0, 'method': 0, 'grid_spacing': 0.1,
                          'size': 1, 'wavefunction': 0, 'spin': 0, 'excitation': 0, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.05 0 0 0 0.1 1 0 0 0 1')

        self.style.set_parameters(isovalue=0.2)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.2, 'orbital_id': 0, 'display': 0, 'method': 0, 'grid_spacing': 0.1,
                          'size': 1, 'wavefunction': 0, 'spin': 0, 'excitation': 0, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.2 0.0 0.0 0.0 0.1 1.0 0.0 0.0 0.0 1.0')

        self.style.set_parameters(orbital_id=2.)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.2, 'orbital_id': 2, 'display': 0, 'method': 0, 'grid_spacing': 0.1,
                          'size': 1, 'wavefunction': 0, 'spin': 0, 'excitation': 0, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.2 2.0 0.0 0.0 0.1 1.0 0.0 0.0 0.0 1.0')

        self.style.set_parameters(display=1.)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.2, 'orbital_id': 2, 'display': 1, 'method': 0, 'grid_spacing': 0.1,
                          'size': 1, 'wavefunction': 0, 'spin': 0, 'excitation': 0, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.2 2.0 1.0 0.0 0.1 1.0 0.0 0.0 0.0 1.0')

        self.style.set_parameters(method=4.)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.2, 'orbital_id': 2, 'display': 1, 'method': 4, 'grid_spacing': 0.1,
                          'size': 1, 'wavefunction': 0, 'spin': 0, 'excitation': 0, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.2 2.0 1.0 4.0 0.1 1.0 0.0 0.0 0.0 1.0')

        self.style.set_parameters(size=2.5)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.2, 'orbital_id': 2, 'display': 1, 'method': 4, 'grid_spacing': 0.1,
                          'size': 2.5, 'wavefunction': 0, 'spin': 0, 'excitation': 0, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.2 2.0 1.0 4.0 0.1 2.5 0.0 0.0 0.0 1.0')

        self.style.set_parameters(wavefunction=7.)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.2, 'orbital_id': 2, 'display': 1, 'method': 4, 'grid_spacing': 0.1,
                          'size': 2.5, 'wavefunction': 7, 'spin': 0, 'excitation': 0, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.2 2.0 1.0 4.0 0.1 2.5 7.0 0.0 0.0 1.0')

        self.style.set_parameters(spin=3.)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.2, 'orbital_id': 2, 'display': 1, 'method': 4, 'grid_spacing': 0.1,
                          'size': 2.5, 'wavefunction': 7, 'spin': 3, 'excitation': 0, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.2 2.0 1.0 4.0 0.1 2.5 7.0 3.0 0.0 1.0')

        self.style.set_parameters(excitation=5.)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.2, 'orbital_id': 2, 'display': 1, 'method': 4, 'grid_spacing': 0.1,
                          'size': 2.5, 'wavefunction': 7, 'spin': 3, 'excitation': 5, 'step': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.2 2.0 1.0 4.0 0.1 2.5 7.0 3.0 5.0 1.0')

        self.style.set_parameters(step=12.)
        self.assertEqual(self.style.method, DRAW_ORBITAL)
        self.assertEqual(self.style.get_parameters(),
                         {'isovalue': 0.2, 'orbital_id': 2, 'display': 1, 'method': 4, 'grid_spacing': 0.1,
                          'size': 2.5, 'wavefunction': 7, 'spin': 3, 'excitation': 5, 'step': 12})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Orbital 0.2 2.0 1.0 4.0 0.1 2.5 7.0 3.0 5.0 12.0')

    def test_beads(self):
        self.style.method = DRAW_BEADS
        self.assertEqual(self.style.method, DRAW_BEADS)
        self.assertEqual(self.style.get_parameters(), {'size': 1, 'resolution': 12})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Beads')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_BEADS)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 12})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Beads 3.0 12')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_BEADS)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 20})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Beads 3.0 20.0')

    def test_dotted(self):
        self.style.method = DRAW_DOTTED
        self.assertEqual(self.style.method, DRAW_DOTTED)
        self.assertEqual(self.style.get_parameters(), {'size': 1, 'resolution': 12})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Dotted')

        self.style.set_parameters(size=3.)
        self.assertEqual(self.style.method, DRAW_DOTTED)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 12})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Dotted 3.0 12')

        self.style.set_parameters(resolution=20.)
        self.assertEqual(self.style.method, DRAW_DOTTED)
        self.assertEqual(self.style.get_parameters(), {'size': 3, 'resolution': 20})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Dotted 3.0 20.0')

    def test_solvent(self):
        self.style.method = DRAW_SOLVENT
        self.assertEqual(self.style.method, DRAW_SOLVENT)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 0, 'detail': 7, 'method': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Solvent')

        self.style.set_parameters(detail=3.)
        self.assertEqual(self.style.method, DRAW_SOLVENT)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 0, 'detail': 3, 'method': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Solvent 0 3.0 1')

        self.style.set_parameters(probe_radius=1.)
        self.assertEqual(self.style.method, DRAW_SOLVENT)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1, 'detail': 3, 'method': 1})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Solvent 1.0 3.0 1.0')

        self.style.set_parameters(method=0.)
        self.assertEqual(self.style.method, DRAW_SOLVENT)
        self.assertEqual(self.style.get_parameters(), {'probe_radius': 1, 'detail': 3, 'method': 0})
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Solvent 1.0 3.0 0.0')
