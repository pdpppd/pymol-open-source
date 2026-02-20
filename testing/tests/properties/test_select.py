'''
Testing selection by properties
'''

import unittest
from pymol import cmd, testing

class TestProperties(testing.PyMOLTestCase):
    def testNumeric(self):
        cmd.fragment('ala')
        cmd.alter_state(1, 'all', 'p.x, p.y, p.index = x, y, index')
        cmd.select('sele_p_x', 'p.x < 0')
        cmd.select('sele_p_y', 'p.y > 0')
        cmd.select('sele_p_i', 'p.index = 3')
        cmd.select('sele_x', 'x < 0.0')
        cmd.select('sele_y', 'y > 0.0')
        cmd.select('sele_i', 'index 3')
        counts = [
            cmd.count_atoms('sele_p_x'),
            cmd.count_atoms('sele_x'),
            cmd.count_atoms('sele_x & sele_p_x'),
        ]
        self.assertEqual(counts[0], 8)
        self.assertEqual(counts[0], counts[1])
        self.assertEqual(counts[0], counts[2])
        counts = [
            cmd.count_atoms('sele_p_y'),
            cmd.count_atoms('sele_y'),
            cmd.count_atoms('sele_y & sele_p_y'),
        ]
        self.assertEqual(counts[0], 6)
        self.assertEqual(counts[0], counts[1])
        self.assertEqual(counts[0], counts[2])
        counts = [
            cmd.count_atoms('sele_i'),
            cmd.count_atoms('sele_p_i'),
            cmd.count_atoms('sele_i & sele_p_i'),
        ]
        self.assertEqual(counts[0], 1)
        self.assertEqual(counts[0], counts[1])
        self.assertEqual(counts[0], counts[2])

    def testString(self):
        cmd.fragment('ala')
        cmd.alter('all', 'p.name = name')
        cmd.select('sele_p_x', 'p.name in C+N+O')
        cmd.select('sele_x', 'name C+N+O')
        counts = [
            cmd.count_atoms('sele_p_x'),
            cmd.count_atoms('sele_x'),
            cmd.count_atoms('sele_x & sele_p_x'),
        ]
        self.assertEqual(counts[0], 3)
        self.assertEqual(counts[0], counts[1])
        self.assertEqual(counts[0], counts[2])

    def testCast(self):
        cmd.fragment('ala')
        cmd.alter('all', 'p.s, p.index = "10e0", index')
        cmd.alter('index 2+3', 'p.s = "30e0"')
        cmd.select('sele_p_s', 'p.s > 20')
        cmd.select('sele_p_i', 'p.index in 2+3')
        counts = [
            cmd.count_atoms('sele_p_s'),
            cmd.count_atoms('sele_p_i'),
        ]
        self.assertEqual(counts[0], 2)
        self.assertEqual(counts[1], 2)

    def testObjectPropertyNumeric(self):
        """Object-level properties (set_property) should be queryable in
        selections, falling back from atom-level to CoordSet-level."""
        cmd.pseudoatom('obj1')
        cmd.pseudoatom('obj2')
        cmd.set_property('CD', 7.9, 'obj1')
        cmd.set_property('CD', 12.0, 'obj2')

        # numeric comparisons
        self.assertEqual(cmd.count_atoms('p.CD < 8'), 1)
        self.assertEqual(cmd.count_atoms('p.CD > 10'), 1)
        self.assertEqual(cmd.get_object_list('p.CD < 8'), ['obj1'])
        self.assertEqual(cmd.get_object_list('p.CD > 10'), ['obj2'])

    def testObjectPropertyString(self):
        """Object-level string properties should be queryable via p. selections."""
        cmd.pseudoatom('obj1')
        cmd.pseudoatom('obj2')
        cmd.set_property('tag', 'alpha', 'obj1')
        cmd.set_property('tag', 'beta', 'obj2')

        self.assertEqual(cmd.get_object_list('p.tag in alpha'), ['obj1'])
        self.assertEqual(cmd.get_object_list('p.tag in beta'), ['obj2'])

    def testAtomPropertyTakesPrecedence(self):
        """Atom-level properties should take precedence over object-level."""
        cmd.pseudoatom('obj1')
        cmd.set_property('val', 100.0, 'obj1')
        cmd.set_atom_property('val', 5.0, 'obj1')

        # Atom property (5.0) should be used, not object property (100.0)
        self.assertEqual(cmd.count_atoms('p.val < 10'), 1)
        self.assertEqual(cmd.count_atoms('p.val > 50'), 0)

    def testObjectPropertyFallbackWithOtherAtomProps(self):
        """When an atom has some properties but not the queried one,
        fall back to the CoordSet for the queried property."""
        cmd.pseudoatom('obj1')
        cmd.set_atom_property('other', 42, 'obj1')
        cmd.set_property('CD', 7.9, 'obj1')

        # Atom has 'other' but not 'CD', should fall back to CoordSet for 'CD'
        self.assertEqual(cmd.count_atoms('p.CD < 8'), 1)
        # 'other' is on the atom, should still work directly
        self.assertEqual(cmd.count_atoms('p.other = 42'), 1)
