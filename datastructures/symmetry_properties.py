# Copyright (c) 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Software available at https://github.com/sandialabs/ISMAGS
# (POC) Mark DeBonis (mjdebon@sandia.gov)

class SymmetryProperties:
    """Groups all information on the symmetric properties of the motif.

    Attributes:
        number_of_nodes(int): Number of motif nodes.
        smaller(dict): Smaller set of nodes.
        larger(dict): Larger set of nodes.
        permutations(list(list(int))): List of permutations.

    Methods:
        add_permutation(permutation): Add a permutation.
        fix(motif_node_id, orbits): Extracts symmetry-breaking constraints from the orbit
            of the specified motif node.
        add_constraint(lower_id, higher_id): Adds a constraint of the form lowerID higherID to the symmetric
            properties. Constraints are transitively propagated.

    TODO: Determine if number_of_nodes is needed. Doesn't SEEM to be used.

    """
    def __init__(self, number_of_nodes=0, smaller=None, larger=None, permutations=None):
        """Initialize symmetric properties of the motif.

        Keyword Args:
            number_of_nodes (int): Number of motif nodes. Defaults to 0.
            smaller (dict): Smaller set of nodes. Defaults to None.
            larger (dict): Larger set of nodes. Defaults to None.
            permutations (list(list(int))): List of permutations. Defaults to None.
        """
        self.number_of_nodes = number_of_nodes

        if smaller is None:
            self.smaller = dict()
        else:
            self.smaller = smaller

        if larger is None:
            self.larger = dict()
        else:
            self.larger = larger

        if permutations is None:
            self.permutations = []
        else:
            self.permutations = permutations

    def add_permutation(self, permutation):
        """Add a permutation.

        Args:
            permutation(list(int)): Permutation to add.
        """
        self.permutations.append(permutation)

    def fix(self, motif_node_id, orbits):
        """ Extracts symmetry-breaking constraints from the orbit
            of the specified motif node.

        Args:
            motif_node_id(int): Motif node to generate constraints for.
            orbits(list(int)): Orbit partition.
        """
        if motif_node_id >= len(orbits):
            return

        orbit = orbits[motif_node_id]
        for i in range(motif_node_id+1, len(orbits)):
            if orbit == orbits[i]:
                self.add_constraint(motif_node_id, i)

    def add_constraint(self, lower_id, higher_id):
        """Adds a constraint of the form lowerID higherID to the symmetric
            properties. Constraints are transitively propagated.

        Args:
            lower_id(int): ID of the lower motif node.
            higher_id(int): ID of the higher motif node.
        """
        if lower_id in self.smaller:
            smaller_set_a = self.smaller[lower_id]
        else:
            smaller_set_a = set()
            self.smaller[lower_id] = smaller_set_a

        if higher_id in self.smaller:
            smaller_set_b = self.smaller[higher_id]
        else:
            smaller_set_b = set()
            self.smaller[higher_id] = smaller_set_b

        if lower_id in self.larger:
            larger_set_a = self.larger[lower_id]
        else:
            larger_set_a = set()
            self.larger[lower_id] = larger_set_a

        if higher_id in self.larger:
            larger_set_b = self.larger[higher_id]
        else:
            larger_set_b = set()
            self.larger[higher_id] = larger_set_b

        smaller_set_a.add(higher_id)
        larger_set_b.add(lower_id)

        smaller_set_a.update(smaller_set_b)
        larger_set_b.update(larger_set_a)

        for i in smaller_set_b:
            self.larger[i].add(lower_id)

        for i in larger_set_a:
            self.smaller[i].add(higher_id)
