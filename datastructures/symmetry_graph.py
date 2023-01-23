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

from motifs.motif import Motif


class SymmetryGraph:
    """Keeps track of the ordered pair partition (OPP) state during the
       motif analysis.

    Attributes:
        motif (Motif): Subgraph/Motif being analyzed.
        colors_to_recheck (set(int)): Used for rechecking during refinement.
        top_motif_node_to_color (list(int)): Top row of OPP. Keeps track of
            node partition using different integers.
        color_to_bottom_motif_node (dict(int,list(int))): Bottom row of OPP.
            For each branch in the search tree, keeps track of order of nodes
        color_to_top_motif_node (dict(int,list(int))): Top row of OPP.
            For each branch in the search tree, keeps track of order of nodes.

    Methods:
        refine_colors(color): Performs a motif refinement, starting with a specific
            motif partition cell/color that needs refinement.
        map_node_between_partitions(top_id, bottom_id, split_color): Maps a specific node in the top partition
            to a node in the bottom partition.

    """

    def __init__(self, motif=None):
        """Create an initial ordered pair partition (OPP).

        Keyword Args:
            motif (Motif): Subgraph/Motif being analyzed. Defaults to None.
        """

        if motif is None:
            self.motif = Motif()
        else:
            self.motif = motif

        self.colors_to_recheck = set()
        self.top_motif_node_to_color = [0] * self.motif.number_of_motif_nodes
        self.color_to_bottom_motif_node = {}
        self.color_to_top_motif_node = {}

        list1 = []
        list2 = []
        for i in range(self.motif.number_of_motif_nodes):
            list1.append(i)
            list2.append(i)

        self.color_to_bottom_motif_node[0] = list2.copy()
        self.color_to_top_motif_node[0] = list1.copy()

    def refine_colors(self, color):
        """Performs a motif refinement, starting with a specific
            motif partition cell/color that needs refinement

        Args:
            color (int): Color starting cell for refinement

        Return:
            False if resulting OPP is invalid, True otherwise
        """
        ok = self._refine(color)

        while ok and len(self.colors_to_recheck) > 0:

            # In order to grab the first element with minimal overhead
            # pop and re-add the element as sets can't be directly indexed
            color_to_check = self.colors_to_recheck.pop()
            self.colors_to_recheck.add(color_to_check)

            ok = self._refine(color_to_check)
            self.colors_to_recheck.remove(color_to_check)

        return ok

    def map_node_between_partitions(self, top_id, bottom_id, split_color):
        """Maps a specific node in the top partition to a node in the bottom
            partition.

        Args:
            top_id (int): Motif node in top partition.
            bottom_id (int): Motif node in bottom partition.
            split_color (int): ID of partition cell.

        Return:
            New OPP state.
        """
        new_sym = SymmetryGraph(self.motif)

        new_sym.top_motif_node_to_color = self.top_motif_node_to_color.copy()

        new_sym.color_to_bottom_motif_node = copy.deepcopy(self.color_to_bottom_motif_node)
        new_sym.color_to_top_motif_node = copy.deepcopy(self.color_to_top_motif_node)

        new_sym.color_to_bottom_motif_node[split_color].remove(bottom_id)
        new_sym.color_to_top_motif_node[split_color].remove(top_id)

        new_color = len(self.color_to_top_motif_node)

        new_sym.color_to_bottom_motif_node[new_color] = [bottom_id]
        new_sym.color_to_top_motif_node[new_color] = [top_id]

        if new_sym.refine_colors(new_color):
            return new_sym
        return None

    def _refine(self, color):
        """Performs one step of the refinement procedure.

        Args:
            color (int): Color/Cell to be refined.

        Return:
            True if refinement was successful, False otherwise.
        """
        number_of_motif_nodes = self.motif.number_of_motif_nodes

        degrees_top = [list() for _ in range(number_of_motif_nodes)]
        degrees_bottom = [list() for _ in range(number_of_motif_nodes)]

        top_nodes = self.color_to_top_motif_node[color]
        bottom_nodes = self.color_to_bottom_motif_node[color]

        number_of_link_types = len(self.motif.link_types)

        for i in range(number_of_motif_nodes):
            degrees_bottom[i] = [0] * (number_of_link_types * 2)
            degrees_top[i] = [0] * (number_of_link_types * 2)

        reached_colors = set()

        for node in top_nodes:
            links = self.motif.links[node]
            links_degrees = self.motif.final_connections[node]

            for j in range(len(links)):
                i = links_degrees[j]
                motif_link = links[j]

                if motif_link.link_type.motif_link.motif_link_id == motif_link.motif_link_id:
                    degrees_top[i][motif_link.link_type.link_type_id] += 1
                    degrees_top[node][number_of_link_types + motif_link.link_type.link_type_id] += 1
                else:
                    degrees_top[node][motif_link.link_type.link_type_id] += 1
                    degrees_top[i][number_of_link_types + motif_link.link_type.link_type_id] += 1

                reached_colors.add(self.top_motif_node_to_color[i])

        for node in bottom_nodes:
            links = self.motif.links[node]
            links_degrees = self.motif.final_connections[node]

            for j in range(len(links)):
                i = links_degrees[j]
                motif_link = links[j]

                # TODO: In order to stay closer to the Java code we should probably create a __eq__ function for MotifLink()
                if motif_link.link_type.motif_link.motif_link_id == motif_link.motif_link_id:
                    degrees_bottom[i][motif_link.link_type.link_type_id] += 1
                    degrees_bottom[node][number_of_link_types + motif_link.link_type.link_type_id] += 1
                else:
                    degrees_bottom[node][motif_link.link_type.link_type_id] += 1
                    degrees_bottom[i][number_of_link_types + motif_link.link_type.link_type_id] += 1


        for integer in reached_colors:
            nodes_in_color = self.color_to_top_motif_node[integer]
            current_color_mapping = {}
            current_color_mapping[integer] = degrees_top[nodes_in_color[0]]
            start_set = []
            start_set.append(nodes_in_color[0])
            self.color_to_top_motif_node[integer] = start_set

            for i in range(1, len(nodes_in_color)):
                node = nodes_in_color[i]
                i_s = degrees_top[node]
                added = False

                for entry in current_color_mapping.items():
                    connections_color = entry[1]

                    if self._compare_rows(connections_color, i_s):
                        self.color_to_top_motif_node[entry[0]].append(node)
                        self.top_motif_node_to_color[node] = entry[0]
                        added = True
                        break

                if not added:
                    new_color = len(self.color_to_top_motif_node)
                    self.colors_to_recheck.add(new_color)
                    self.colors_to_recheck.add(color)

                    new_set = []
                    new_set.append(node)
                    current_color_mapping[new_color] = i_s
                    self.color_to_top_motif_node[new_color] = new_set
                    self.top_motif_node_to_color[node] = new_color

            nodes_in_bottom_color = self.color_to_bottom_motif_node[integer]
            self.color_to_bottom_motif_node.pop(integer)

            for i in range(len(nodes_in_bottom_color)):
                node_id = nodes_in_bottom_color[i]
                self._refine_bottom(node_id, degrees_bottom, current_color_mapping)

            for integer1 in self.color_to_top_motif_node.keys():
                bottom_set = self.color_to_bottom_motif_node.get(integer1)
                top_set = self.color_to_top_motif_node.get(integer1)
                if bottom_set is None or len(top_set) != len(bottom_set):
                    return False

        return True

    def _refine_bottom(self, node_id, degrees_bottom, current_color_mapping):
        """Refine the bottom row of the motif.

        Args:
            node_id (int): ID of the given node.
            degrees_bottom (list(list(int))): Number of degrees of the bottom row.
            current_color_mapping (list(list(int))): The current mapping of colors.
        """
        i_s = degrees_bottom[node_id]
        for entry in current_color_mapping.items():
            if self._compare_rows(entry[1], i_s):
                color_for_bottom = entry[0]
                get = self.color_to_bottom_motif_node.get(color_for_bottom)
                if get is None:
                    get = []
                    self.color_to_bottom_motif_node[color_for_bottom] = get
                get.append(node_id)
                return
        return

    def _compare_rows(self, a, b):
        """Compare rows of the motif.

        Args:
            a (int): First row to compare.
            b (int): Second row to compare.

        Returns:
            True if the rows are identical, False otherwise.
        """
        for i in range(len(b)):
            if(a[i] != b[i]):
                return False
        return True
