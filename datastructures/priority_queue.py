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

import heapq
import sys

class PriorityObject:
    """Base object used in our priority queues.

    Object used to determine priority and track its
    own place in the queue.

    Attributes:
        start_node(Node): The node asociated with this object
        from_position(int): motif/network position identification
                            based on the position the node is transitioned from
        to_position(int): motif/network position identification
                            based on the position the node is transitioned to
        num_neighbors(int): number of neighbors to this neighbor

    Args:
        start_node(Node): The node asociated with this object
        from_position(int): motif/network position identification
                            based on the position the node is transitioned from
        to_position(int): motif/network position identification
                            based on the position the node is transitioned to
        num_neighbors(int): number of neighbors to this neighbor

    """
    def __init__(self, start_node, from_position, to_position, num_neighbors):
        self.start_node = start_node
        self.from_position = from_position
        self.to_position = to_position
        self.num_neighbors = num_neighbors

    def __repr__(self):
        repr = "<" + str(self.start_node) + "," + \
                str(self.from_position) + "," + \
                str(self.to_position) + "," + \
                str(self.num_neighbors) + ">"
        return repr

    def __str__(self):
        return repr(self)

    def __eq__(self, other):
        if isinstance(self, other.__class__):
            return self.start_node == other.start_node and\
                self.from_position == other.from_position and\
                self.to_position == other.to_position and\
                self.num_neighbors == other.num_neighbors
        return False

    def __lt__(self, other):
        return self.num_neighbors < other.num_neighbors

class PriorityQueue:
    """Priority Queue Wrapper

    This wraps priority queue functionality.
    Mostly implemented in heapq. The removal of items
    and reheapify-ing might be a little inefficient.

    Attributes:
        pq(list): The priority queue. Sorted to a priorirty queue with heapq
        motif_map(dictionary(int, Node)): Map of elements used for arbitrary
                        deletion based on motif.

    TODO: Add checks of passed in lists and dicts. Current module implementation
        won't use it this way but may be best to be sure and safe.

    """
    def __init__(self, pq=None, motif_map=None):
        if pq is None:
            self.pq = []
        else:
            self.pq = pq

        if motif_map is None:
            self.motif_map = {}
        else:
            self.motif_map = motif_map

    def __len__(self):
        return len(self.pq)

    def add(self, element):
        heapq.heappush(self.pq, element)
        self.motif_map[element.from_position] = element

    def peek(self):
        if len(self.pq):
            return self.pq[0]
        return None

    def remove_object(self, priority_obj):
        obj_index = self.pq.index(priority_obj)
        self.pq.pop(obj_index)
        heapq.heapify(self.pq)

    def remove_motif_node(self, motif_node):
        if motif_node in self.motif_map:
            object = self.motif_map.pop(motif_node)
            self.remove_object(object)


class PriorityQueueMap:
    """List of priority queues

    Maps priority queues according to indices.

    Args:
        size(int): size of the pq_map

    TODO: Add checks of passed in lists. Current module implementation
        won't use it this way but may be best to be sure and safe.

    """

    def __init__(self, size, pq_map=None):
        if pq_map is None:
            self.pq_map = []
        else:
            self.pq_map = pq_map
        for _ in range(size):
            self.pq_map.append(PriorityQueue())

    def add(self, priority_obj):
        self.pq_map[priority_obj.to_position].add(priority_obj)

    def poll(self, indices):
        index_iter = iter(indices)
        min_pq = self.pq_map[next(index_iter)]
        priority_obj = min_pq.peek()
        if priority_obj is None:
            min_score = sys.maxsize
        else:
            min_score = priority_obj.num_neighbors

        for element in index_iter:
            next_pq = self.pq_map[element]
            if next_pq is not None:
                if len(next_pq):
                    score = next_pq.peek().num_neighbors
                    if score < min_score:
                        min_score = score
                        min_pq = self.pq_map[element]

        return min_pq.peek()

    def remove_object(self, priority_obj):
        self.pq_map[priority_obj.to_position].remove_object(priority_obj)

    def remove_motif_node(self, motif_node, i):
        self.pq_map[i].remove_motif_node(motif_node)

