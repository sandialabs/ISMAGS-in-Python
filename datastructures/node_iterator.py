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

import sys
from collections import deque

class NodeIterator:
    """Class that Keeps track of all lists that need to be intersected to obtain candidates graph nodes.

    Attributes:
        nodes (list): Node candidates.
        parent (NodeIterator): Parent NodeIterator.
        description (int): Node description.
        motif_node_id (int): Motif node for which the NodeIterator will determine candidates.
        min_set_size (int): Size of the smallest set of nodes
        neighbor_lists (deque): Stack of nodes that constrained by a given graph node
        node_causing_restriction (deque): Stack of constraining nodes
        initial_lists (list[list]): List of lists containing the intial nodes that aren't constrained

    Methods:
        add_restriction_list(node_list, node=None): Adds a new list of candidate nodes for the
            motif node during the initialization or if initialization has occured
            adds a new list of candidate nodes for the motif node to the set of constraining lists.
        remove_restriction_list(node): Removes the constraining list of candidate
            nodes induced by the graph node
        get_node_set(): Retrieve the candiate graph nodes.
        node_index(node_list, target): Perform a binary search to find the given node target.
        intersect(minimum, maximum): Creates a NodeIterator based on the constraint lists.

    """

    def __init__(self, nodes=None, parent=None, motif_node_id=0):
        """initialize a NodeIterator with or without an intital set of node candidates.

        Keyword Args:
            nodes (list): Node candidates. Defaults to None.
            parent (NodeIterator): Parent NodeIterator. Defaults to None.
            motif_node_id (int): Motif node for which the NodeIterator will determine
                candidates. Defaults to 0.
        """
        self.nodes = nodes
        self.parent = parent

        if nodes is None and parent is None:
            self.motif_node_id = motif_node_id
            self.min_set_size = sys.maxsize
        else:
            self.motif_node_id = parent.motif_node_id
            self.min_set_size = len(nodes)

        self.neighbor_lists = deque() # Stack
        self.node_causing_restriction = deque() # Stack
        self.initial_lists = []

    def add_restriction_list(self, node_list, node=None):
        """Adds a new list of candidate nodes for the motif node during the
            initialization or if initialization has occured adds a new list of
            candidate nodes for the motif node to the set of constraining lists.

        Args:
            node_list (list): List of candidate nodes.
            node (Node): Graph node causing the constraint. Defaults to None.
        """
        if node is None:
            if len(node_list) < self.min_set_size:
                self.initial_lists.insert(0, node_list)
                self.min_set_size = len(node_list)
            else:
                self.initial_lists.append(node_list)
        else:
            if node_list not in self.neighbor_lists:
                self.neighbor_lists.append(node_list)
                if len(node_list) < self.min_set_size:
                    self.min_set_size = len(node_list)
                self.node_causing_restriction.append(node)

    def remove_restriction_list(self, node):
        """Removes the constraining list of candidate nodes induced by the graph node

        Args:
            node (Node): Graph node inducing the constraining list
        """
        # TODO (mjfadem): Is ID the best thing to check here for the nodes?
        while len(self.node_causing_restriction) > 0 and self.node_causing_restriction[-1].id == node.id:
            self.neighbor_lists.pop()
            self.node_causing_restriction.pop()

    def get_node_set(self):
        """Retrieve the candiate graph nodes.

        Returns:
            The candidate graph nodes. If no set has been calculated
            (=initially) a linear sweep is used
        """
        if self.nodes is not None:
            return self.nodes
        else:
            initial_lists_length = len(self.initial_lists)
            if initial_lists_length == 1:
                return self.initial_lists[0]
            iter_list = []
            for i in range(initial_lists_length):
                iter_list.append(iter(self.initial_lists[i]))

            # All of the following exception logic is mimicking the Java label logic
            try:
                node_i = next(iter_list[0])
            except StopIteration:
                return []
            j = 1
            sim = 1
            node_set = []
            while True:
                try:
                    node_j = next(iter_list[j])
                    try:
                        while node_i - node_j > 0:
                            try:
                                node_j = next(iter_list[j])
                                sim = 1
                            except StopIteration:
                                raise ValueError
                    except ValueError:
                        break
                    try:
                        while node_i - node_j < 0:
                            try:
                                node_i = next(iter_list[0])
                                sim = 1
                            except StopIteration:
                                raise ValueError
                    except ValueError:
                        break
                    if node_j.id == node_i.id:
                        sim += 1
                        if sim == initial_lists_length:
                            node_set.append(node_j)
                            sim = 1
                    else:
                        sim = 1
                    j = (j % (initial_lists_length - 1)) + 1
                except StopIteration:
                    break

            return list({n for nodes in self.initial_lists for n in nodes})

    def node_index(self, node_list, target):
        """Perform a binary search to find the given node target.

        Args:
            node_list (list): List of nodes to search
            target (node): Target node to find in the node_list

        Returns:
            middle: Either the position of the given target or where the target would be inserted.
                This "should" be inline with the original algorithm and Java's implementation of
                binary search.
        """
        left = 0
        right = len(node_list)-1
        while left <= right:
            middle = left + (right - left) // 2

            if node_list[middle].id == target.id:
                return middle
            elif target.id > node_list[middle].id:
                left = middle + 1
            elif target.id < node_list[middle].id:
                right = middle - 1

        return -middle-1

    def intersect(self, minimum, maximum):
        """Creates a NodeIterator based on the constraint lists.

        Args:
            minimum (Node): lower bound on the nodes (ID-based).
            maximum (Node): upper bound on the nodes (ID-based).

        Returns:
            child NodeIterator object or None if no candidates were found
        """
        if len(self.neighbor_lists) == 0:
            return None

        result = []
        stack_size = len(self.neighbor_lists)
        smallest_set = stack_size - 1
        node_set = self.neighbor_lists[smallest_set]
        if len(node_set) > self.min_set_size:
            for i in range(stack_size - 1):
                neighbor_node_set = self.neighbor_lists[i]
                if len(neighbor_node_set) < self.min_set_size:
                    node_set = neighbor_node_set
                    self.min_set_size = len(node_set)
                    smallest_set = i

        nodes = node_set
        start_index = 0
        if minimum is not None:
            p = self.node_index(nodes, minimum)
            if p >= 0:
                start_index = p + 1
            else:
                start_index = -p - 1

        end_index = len(nodes)
        if maximum is not None:
            p = self.node_index(nodes, maximum)
            if p >= 0:
                end_index = p
            else:
                end_index = -p - 1

        # All of the following exception logic is mimicking the Java label logic
        for node in nodes[start_index:end_index]:
            if node.used:
                continue
            try:
                for i in range(len(self.neighbor_lists)):
                    if i != smallest_set and node not in self.neighbor_lists[i]:
                        raise ValueError
            except ValueError:
                continue
            result.append(node)
        if len(result) == 0:
            return None
        return NodeIterator(nodes=result, parent=self)
