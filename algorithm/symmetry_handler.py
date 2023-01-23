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

import copy
import math
import sys

from datastructures.priority_queue import (PriorityObject, PriorityQueueMap)
from datastructures.symmetry_graph import SymmetryGraph
from datastructures.symmetry_properties import SymmetryProperties
from motifs.motif import Motif


class SymmetryHandler:
    """Creates a new SymmetryHandler. This class is responsible
        for analysing the motif and providing the constraints
        to the MotifFinder class.

    Attributes:
        mapping (list[NodeIterator]): Handle to NodeIterators containing constraining neighbor lists.
        mapped_nodes (list[Node]): Handle to partial node mapping.
        motif (Motif): Motif to be analysed.
        priority_queue_map (PriorityQueueMap): Priotity queue mapping of the motif nodes.
        mapped_positions (set(int)): Mapped positions of motif nodes.
        smaller (dict): Dictionary of smaller motif nodes for a given motif node ID.
        larger (dict): Dictionary of larger motif nodes for a given motif node ID.
        number_of_orbits (int): Number of orbits in the motif
        symmetric_properties (SymmetryProperties): All information on the symmetric properties of the motif

    Methods:
        get_next_best_iterator(unmapped_motif_nodes): Determines the next motif nodes and candidates to be mapped.
        map_node(motif_node, graph_node): Maps a graph node to a motif node and updates the neighbor lists used for intersecting.
        remove_node_mapping(motif_node, graph_node): Un-maps a graph node previously mapped to a motif node, ensuring
            consistency in constraining neighbour lists.

    """

    def __init__(self, mapping=None, motif=Motif(), mapped_nodes=None):
        """Initialize a SymmetryHandler object to deal with the specified motif.

        Keywords Args:
            mapping (list[NodeIterator]): Handle to NodeIterators containing constraining neighbor lists.
            motif (Motif): Motif to be analyzed.
            mapped_nodes (list[Node]): Handle to partial node mapping.
        """
        if mapping is None:
            self.mapping = []
        else:
            self.mapping = mapping

        if mapped_nodes is None:
            self.mapped_nodes = []
        else:
            self.mapped_nodes = mapped_nodes

        self.motif = motif
        self.priority_queue_map = PriorityQueueMap(len(self.mapping))
        self.mapped_positions = set()
        self.smaller = dict()
        self.larger = dict()
        self.number_of_orbits = 0

        self.symmetric_properties = self._analyze_motif(motif)

    def get_next_best_iterator(self, unmapped_motif_nodes):
        """Determines the next motif nodes and candidates to be mapped.

        Args:
            unmapped_motif_nodes (set(int)): Unmapped motif nodes from which to select the next node.

        Returns:
            Next motif node and graph node candidates.
        """
        poll = self.priority_queue_map.poll(unmapped_motif_nodes)
        motif_node_id = poll.to_position
        node_iterator = self.mapping[motif_node_id]

        # Determine lower bound for graph node candidates
        if motif_node_id in self.larger:
            min_set = self.larger[motif_node_id]
        else:
            min_set = None
        min_node = None
        min_value = -math.inf

        if min_set is not None:
            for integer in min_set:
                if integer in self.mapped_positions and min_value < self.mapped_nodes[integer].id:
                    min_value = self.mapped_nodes[integer].id
                    min_node = self.mapped_nodes[integer]

        # Determine upper bound for graph node candidates
        if motif_node_id in self.larger:
            max_set = self.smaller[motif_node_id]
        else:
            max_set = None

        max_node = None
        max_value = sys.maxsize

        if max_set is not None:
            for integer in max_set:
                if integer in self.mapped_positions and max_value > self.mapped_nodes[integer].id:
                    max_value = self.mapped_nodes[integer].id
                    max_node = self.mapped_nodes[integer]
                    # Abort when bounds conflict
                    if min_value > max_value:
                        return None

        # Determine nodes by intersecting using the bounds
        return node_iterator.intersect(min_node, max_node)

    def map_node(self, motif_node, graph_node):
        """Maps a graph node to a motif node and updates the neighbor lists used for intersecting.

        Args:
            motif_node (int): Motif node to be mapped on.
            graph_node (Node): Graph node to be mapped.

        Returns:
            True if graph_node is suitable for mapping on the motif node, otherwise False.
        """

        # Prior to this step the motif will be initialized on the MotifFinder side of things
        connections = self.motif.final_connections[motif_node]
        restrictions = self.motif.links[motif_node]

        for (connection, motif_link) in zip(connections, restrictions):
            if self.mapped_nodes[connection] is not None:
                continue

            links = graph_node.neighbours_per_type[motif_link.motif_link_id]
            if links is None:
                return False
            else:
                self.mapping[connection].add_restriction_list(links, node=graph_node)
                self.priority_queue_map.add(PriorityObject(start_node=graph_node,
                                                           from_position=motif_node,
                                                           to_position=connection,
                                                           num_neighbors=len(links)))
        return True

    def remove_node_mapping(self, motif_node, graph_node):
        """Un-maps a graph node previously mapped to a motif node, ensuring
            consistency in constraining neighbour lists.

        Args:
            motif_node (int): Motif node mapped to.
            graph_node (Node): Graph node mapped to motif node.
        """
        neighbors = self.motif.final_connections[motif_node]
        for i in neighbors:
            self.mapping[i].remove_restriction_list(graph_node)
            self.priority_queue_map.remove_motif_node(motif_node, i)

    def _merge_orbits(self, a, b, orbits):
        """Updates the orbit partitioning by merging the orbit partition
            cells of the specified motif nodes.

        Args:
            a (int): First cell to be merged.
            b (int): Second cell to be merged.
            orbits (list[int]): Current orbit partition.
        """
        orbit_a = orbits[a]
        orbit_b = orbits[b]
        if orbit_a == -1 and orbit_b == -1:
            self.number_of_orbits+=1
            orbits[a] = self.number_of_orbits
            orbits[b] = self.number_of_orbits
        elif orbit_b == -1:
            orbits[b] = orbit_a
        elif orbit_a == -1:
            orbits[a] = orbit_b
        else:
            for i in range(len(orbits)):
                if orbits[i] == orbit_a:
                    orbits[i] = orbit_b

    def _analyze_motif(self, motif):
        """Initializes and returns full motif analysis.

        Args:
            motif (Motif): Motif to be analyzed.
        Return:
            The symmetric properties of the given motif.
        """
        number_of_motif_nodes = motif.number_of_motif_nodes
        symmetry_graph = SymmetryGraph(motif=motif)
        symmetric_properties = SymmetryProperties(number_of_nodes=number_of_motif_nodes, smaller=self.smaller, larger=self.larger)
        orbits = [-1] * number_of_motif_nodes
        self._map_nodes(symmetric_properties, orbits, symmetry_graph)
        return symmetric_properties

    def _map_nodes(self, symmetric_properties, orbits, symmetry_graph, main=True):
        """Recursive motif analysis

        Args:
            symmetric_properties (SymmetryProperties): Stores all permutations and symmetry-breaking constraints for the motif.
            orbits (list(int)): Orbit partitioning of the motif nodes.
            symmetry_graph (SymmetryGraph): Current state in motif analysis.

        Keyword Args:
            main (bool): True if all previously coupled motif nodes are coupled to themselves.
        """
        all_one = True
        split_color = -1
        lowest_unassigned_motif_node = sys.maxsize

        # Determing the next motif node to map
        for i in range(len(symmetry_graph.color_to_top_motif_node)):
            list_i = symmetry_graph.color_to_top_motif_node[i]
            size_i = len(list_i)
            if size_i != 1:
                all_one = False
                try:
                    for motif_node_id in list_i:
                        if motif_node_id < lowest_unassigned_motif_node:
                            split_color = i
                            lowest_unassigned_motif_node = motif_node_id
                            raise ValueError
                except ValueError:
                    continue

        # If all nodes are mapped, export permutation
        if all_one:
            permutation = [0] * symmetry_graph.motif.number_of_motif_nodes
            for j in range(len(permutation)):
                bottom_color = symmetry_graph.color_to_bottom_motif_node[j][0]
                top_color = symmetry_graph.color_to_top_motif_node[j][0]
                permutation[top_color] = bottom_color
                self._merge_orbits(bottom_color, top_color, orbits)

            symmetric_properties.add_permutation(permutation)
            return

        # Map lowest uncoupled node in the upper partition on all motif nodes in the lower partition
        top_split = symmetry_graph.color_to_top_motif_node[split_color]
        top = lowest_unassigned_motif_node
        bottom_split = symmetry_graph.color_to_bottom_motif_node[split_color]
        new_symmetry_graph = copy.deepcopy(symmetry_graph)

        # Deal with the initial OPP and the cases for which the OPP has identical subsets
        if len(top_split) != symmetry_graph.motif.number_of_motif_nodes and all(item in top_split for item in bottom_split):
            permutation = [0] * symmetry_graph.motif.number_of_motif_nodes
            identity_permutation = True
            for k in range(len(new_symmetry_graph.color_to_bottom_motif_node)):
                bottom_nodes = new_symmetry_graph.color_to_bottom_motif_node[k]
                top_nodes = new_symmetry_graph.color_to_top_motif_node[k]
                bottom_node = bottom_nodes[0]
                top_node = top_nodes[0]
                if bottom_node != top_node:
                    self._merge_orbits(bottom_node, top_node, orbits)

            # TODO (mjfadem): The original Java code for this loop makes no sense, seems like it could easily infinitely loop
            for j in range(len(new_symmetry_graph.color_to_bottom_motif_node)):
                bottom_nodes = new_symmetry_graph.color_to_bottom_motif_node[j]
                top_nodes = new_symmetry_graph.color_to_top_motif_node[j]
                if len(bottom_nodes) > 1:
                    bottom_node = bottom_nodes[0]
                    top_node = top_nodes[0]
                    new_symmetry_graph = new_symmetry_graph.map_node_between_partitions(top_node, bottom_node, j)
                    if len(new_symmetry_graph.color_to_bottom_motif_node) != symmetry_graph.motif.number_of_motif_nodes:
                        j = -1 # TODO (mjfadem): This is all kinds of wrong from a logic standpoint
                        continue
                    else:
                        break

            for j in range(len(new_symmetry_graph.color_to_bottom_motif_node)):
                bottom_nodes = new_symmetry_graph.color_to_bottom_motif_node[k]
                top_nodes = new_symmetry_graph.color_to_top_motif_node[k]
                bottom_node = bottom_nodes[0]
                top_node = top_nodes[0]
                permutation[top_node] = bottom_node
                if bottom_node != top_node:
                    identity_permutation = False

            if not identity_permutation:
                symmetric_properties.add_permutation(permutation)
                return

        # Iterate over all possible couplings
        for motif_node in bottom_split:
            if orbits[top] != -1 and orbits[top] == orbits[motif_node]:
                continue
            new_symmetry_graph = symmetry_graph.map_node_between_partitions(top, motif_node, split_color)
            # Keep track of couplings and if the first are coupled to themselves
            new_main = (main and (motif_node == top))
            if new_symmetry_graph is not None:
                self._map_nodes(symmetric_properties, orbits, new_symmetry_graph, main=new_main)

        # Export partial orbit cells as symmetry-breaking constraints
        if main:
            symmetric_properties.fix(top, orbits)
