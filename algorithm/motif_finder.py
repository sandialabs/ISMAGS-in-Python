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

import logging
import sys
import time

from algorithm.symmetry_handler import SymmetryHandler
from datastructures.node_iterator import NodeIterator
from motifs.motif_instance import MotifInstance
from motifs.motif_link import MotifLink
from network.network import Network

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MotifFinder:
    """Creates a new MotifFinder. This class is responsible
        for finding all motif instances.

    Attributes:
        network (Network): Network to be searched.
        symmetry_handler (SymmetryHandler): Symmetry handler to use when finding the motif.
        unmapped_nodes (set(int)): Nodes that need to be mapped to motif.
        used_links (set(set(Node))): Set of links that have been used while finding the motif.
        cancelled (bool): Flag to terminate parallel proccessing. **Unused Currently**.

    Methods:
        find_motif(motif, save_links): Finds and returns all instances of specified motif in the network.

    """

    def __init__(self, network=Network()):
        """Initialize a MotifFinder object to find all motifs in a network.

        Keywords Args:
            network (Network): Network to be searched.
        """
        logger.info("Initializing ISMAGS MotifFinder...")
        self.network = network
        self.symmetry_handler = None # No need to init symmetry_handler ahead of time
        self.unmapped_nodes = set()
        self.used_links = set() # set(set(Node))
        
        self.cancelled = False


    def find_motif(self, motif, save_links=False):
        """Finds and returns all instances of specified motif in the network.

        Args:
            motif (Motif): Motif of which instances need to be found.

        Keywords Args:
            save_links (bool): Keep a set of links used in the result set. Defaults to False.

        Returns:
            All instances of the given motif within the network.
        """

        logger.info("Performing motif search...")
        timer = time.perf_counter()

        for i in range(motif.number_of_motif_nodes):
            self.unmapped_nodes.add(i)

        number_of_motif_nodes = motif.number_of_motif_nodes

        # Determining first motif node to be investigated based on number of
		# edges in network
        mapping = [None for _ in range(number_of_motif_nodes)]
        best_motif_node = -1
        size_of_list_of_best_node = sys.maxsize

        for i in range(number_of_motif_nodes):
            # Determine nodes mappable on node i
            number_of_links = [0] * MotifLink.NUMBER_OF_LINK_IDS
            links_from_i = motif.links[i]
            nodes_connected_to_i = motif.final_connections[i]
            number_of_connections = len(nodes_connected_to_i)
            node_iterator = NodeIterator(motif_node_id=i)
            size_of_smallest_list_node_i = sys.maxsize

            # For each outgoing link, add the list of nodes in the network
			# having that edge type
            for k in range(number_of_connections):
                link = links_from_i[k]
                number_of_links[link.motif_link_id] += 1

                if number_of_links[link.motif_link_id] == 1:
                    nodes_of_type = self.network.get_nodes_of_type(link)
                    node_iterator.add_restriction_list(nodes_of_type)

                    if size_of_smallest_list_node_i > len(nodes_of_type):
                        size_of_smallest_list_node_i = len(nodes_of_type)


            # First node to be mapped is the node with the smallest candidate sublist
            if size_of_smallest_list_node_i < size_of_list_of_best_node:
                size_of_list_of_best_node = size_of_smallest_list_node_i
                best_motif_node = i

            mapping[i] = node_iterator

        instances = set()
        mapped_nodes = [None for _ in range(number_of_motif_nodes)]

        # Initialize symmetry handler to analyze the motif
        self.symmetry_handler = SymmetryHandler(mapping=mapping, motif=motif, mapped_nodes=mapped_nodes)

        if save_links:
            self.used_links = set()

        self._map_next(motif, instances, best_motif_node, mapped_nodes, save_links, number_of_mapped=0)
        logger.info(f"Completed motif search in {time.perf_counter()-timer:.6f} seconds")
        logger.info(f"Found {len(instances)} instances of {motif.description} motif")
        return instances

    def _map_next(self, motif, instances, motif_node, mapped_nodes, save_links, number_of_mapped=0):
        """Recursively called to map graph nodes to the next motif node.

        Args:
            motif (Motif): Subgraph to be searched for.
            instances (set(MotifInstance)): Set to store motif instances in.
            motif_node (int): Next node to be mapped.
            mapped_nodes (list[Node]): Current partial node mapping.
            save_links (bool): Keep a set of links used in the result set.

        Keyword Args:
            number_of_mapped (int): Number of nodes already in the partial mapping. Defaults to 0.
        """
        nodes = self.symmetry_handler.mapping[motif_node].get_node_set()

        # if the current node mapping will complete the mapping, export the instances
        if number_of_mapped == motif.number_of_motif_nodes-1:
            if save_links and len(nodes) > 0:
                for i in range(motif.number_of_motif_nodes):
                    if mapped_nodes[i] == None:
                        continue
                    links = motif.final_connections[i]
                    for j in range(len(links)):
                        if mapped_nodes[links[j]] == None:
                            continue
                        elif links[j] > i:
                            break
                        link = set()
                        link.add(mapped_nodes[i])
                        link.add(mapped_nodes[links[j]])
                        self.used_links.add(link)

            for node in nodes:
                mapped_nodes[motif_node] = node
                instances.add(MotifInstance(mapping=mapped_nodes))
                if save_links:
                    links = motif.final_connections[motif_node]
                    for j in range(len(links)):
                        link = set()
                        link.add(mapped_nodes[links[j]])
                        link.add(node)
                        self.used_links.add(link)

            mapped_nodes[motif_node] = None
        else:
            # For each possible node, map
            node_iterator = iter(nodes)
            self.symmetry_handler.mapped_positions.add(motif_node)
            self.unmapped_nodes.remove(motif_node)

            for node in node_iterator: # This should function the same as the Java while loop...
                mapped_nodes[motif_node] = node
                node.used = True

                # Map graph node to motif node, early termination if graph node
                # does not support all edges of motif node
                success_mapping = self.symmetry_handler.map_node(motif_node, node)

                if success_mapping:
                    # Determine next node to be mapped
                    next_iterator = self.symmetry_handler.get_next_best_iterator(self.unmapped_nodes)
                    if next_iterator is not None:
                        self.symmetry_handler.mapping[next_iterator.motif_node_id] = next_iterator
                        # Recursively call _map_next
                        self._map_next(motif, instances, next_iterator.motif_node_id, mapped_nodes, save_links, number_of_mapped=number_of_mapped+1)
                        # Backtracking
                        self.symmetry_handler.mapping[next_iterator.motif_node_id] = next_iterator.parent

                # Backtracking
                self.symmetry_handler.remove_node_mapping(motif_node, node)
                node.used = False
                mapped_nodes[motif_node] = None

            self.symmetry_handler.mapped_positions.remove(motif_node)
            self.unmapped_nodes.add(motif_node)
