# Test of Wall class and intersection methods
# Copyright (C) 2019  Robin Scheibler, Cyril Cadoux
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# You should have received a copy of the MIT License along with this program. If
# not, see <https://opensource.org/licenses/MIT>.
from __future__ import division

import unittest

import numpy as np
import pyroomacoustics as pra

corners = np.asarray([[0., 1., 1., 0., 0.], [0., 0., 1., 1., 0.], [0., 0., 0., 0., 0.]])
hole = np.asarray([[0.25, 0.75, 0.75, 0.25, 0.25], [0.25, 0.25, 0.75, 0.75, 0.25], [0., 0., 0., 0., 0.]])
holes = np.expand_dims(hole, 0)


class TestWallWithHoles(unittest.TestCase):

    def test_same_wall_true(self):
        w1 = pra.Wall(
            corners, holes=holes
        )
        w2 = pra.Wall(
            corners, holes=holes
        )
        self.assertTrue(w1.same_as(w2))

    def test_same_wall_false_hole_no_hole(self):
        w1 = pra.Wall(
            corners, holes=holes
        )
        w2 = pra.Wall(
            corners
        )
        self.assertTrue(not w1.same_as(w2))

    def test_same_wall_false_no_hole_hole(self):
        w1 = pra.Wall(
            corners
        )
        w2 = pra.Wall(
            corners, holes=holes
        )
        self.assertTrue(not w1.same_as(w2))

    def test_same_wall_false_more_corners(self):

        c2 = np.asarray([[0., 1., 1., 0., 0.1, 0.], [0., 0., 1., 1., 0.5, 0.], [0., 0., 0., 0., 0., 0.]])
        w1 = pra.Wall(
            corners, holes=holes
        )
        w2 = pra.Wall(
            c2, holes=holes
        )
        self.assertTrue(not w1.same_as(w2))

    def test_same_wall_false_more_holes(self):

        h2 = np.asarray([
            [0.25, 0.45, 0.45, 0.25, 0.25],
            [0.25, 0.25, 0.75, 0.75, 0.25],
            [0., 0., 0., 0., 0.]
        ])
        h3 = np.asarray([
            [0.55, 0.75, 0.75, 0.55, 0.25],
            [0.25, 0.25, 0.75, 0.75, 0.25],
            [0., 0., 0., 0., 0.]
        ])
        w1 = pra.Wall(
            corners, holes=np.expand_dims(h2, 0)
        )
        w2 = pra.Wall(
            corners, holes=np.stack((h2, h3))
        )
        self.assertTrue(not w1.same_as(w2))

    def test_same_wall_false_more_hole_corners(self):

        h2 = np.asarray([
            [0.25, 0.75, 0.75, 0.25, 0.3, 0.25],
            [0.25, 0.25, 0.75, 0.75, 0.5, 0.25],
            [0., 0., 0., 0., 0., 0.]
        ])
        w1 = pra.Wall(
            corners, holes=holes
        )
        w2 = pra.Wall(
            corners, holes=np.expand_dims(h2, 0)
        )
        self.assertTrue(not w1.same_as(w2))


if __name__ == "__main__":
    unittest.main()
    # TODO intersection tests
