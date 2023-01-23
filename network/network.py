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
from functools import cmp_to_key

from network.link import Link
from network.node import Node

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Network:
    """Network of nodes and their associated links.

    Attributes:
        nodes_by_id (dict): Dict of nodes within the network based on their IDs.
        nodes_by_description (dict): Dict of nodes within the network based on their descriptions.
        node_sets_departing_from_link (dict{MotifLink: set}): Sets of nodes related to a MotifLink that depart from a given link.
        nodes_with_link (dict): All nodes that share a common/specific link.
        number_of_nodes (int): Number of nodes in the network.
        number_of_links (int): Number of links in the network.

    Methods:
        get_node_by_id(id=None): Retrieve all nodes from network with a specific ID.
        get_nodes_by_description(description=None): Retrieve all nodes from network with a specific description.
        get_nodes_of_type(node_type): Retrieve nodes from network with a specific motif type.
        add_node(node): Add node to the network.
        add_link(link): Add link to the network.
        finalize_network_construction(): Optimize network structure for further processing.
        read_network_from_files(filenames, link_types): Read in network structures from file(s).

    """

    def __init__(self):

        """Initializes a network of nodes and their associated links.
        """

        self.nodes_by_id = {}
        self.nodes_by_description = {}
        self.node_sets_departing_from_link = {}
        self.nodes_with_link = {}
        self.number_of_links = 0

    @property
    def number_of_nodes(self):
        return len(self.nodes_by_id)

    def get_node_by_id(self, id):
        """Retrieve all nodes from network with a specific ID.

            TODO (mjfadem): Provided that all nodes in a network are suppose to be
                unique then this will always return a single node which seems to
                disagree with the original comments. Investigate further and consider
                removing this function if it would be easier/simpler to directly access
                nodes_by_id.

        Args:
            id (int): Node ID to find in network. Defaults to None.

        Returns:
            All nodes with specified ID.
        """
        return self.nodes_by_id[id]

    def get_nodes_by_description(self, description=None):
        """Retrieve all nodes from network with a specific description.

        Args:
            description (str): Node description to find in network. Defaults to None.

        Returns:
            All nodes with specified description, if no description is provided all nodes
            with descriptions are returned.
        """
        if description:
            return self.nodes_by_description[description]
        else:
            return self.nodes_by_description

    def _get_set_of_type(self, set_type):
        """Retrieve a set of nodes that are departing from a link within the network with a specific set type.

        Args:
            set_type (MotifLink): Type of set to retrieve from node_sets_departing_from_link.

        Returns:
            Set of nodes matching the set type, if the set type doesn't exist then
            an empty set is returned.
        """
        if set_type in self.node_sets_departing_from_link:
            node_set = self.node_sets_departing_from_link[set_type]
        else:
            node_set = set(())
            self.node_sets_departing_from_link[set_type] = node_set
        return node_set

    def get_nodes_of_type(self, node_type):
        """Retrieve all nodes with links within the network with a specific node type.

        Args:
            node_type (MotifLink): Type to retrieve from nodes_with_link.

        Returns:
            List of sorted nodes matching the node type, if the type doesn't exist then
            an empty list is returned.
        """
        node_list = self.nodes_with_link[node_type]
        if node_list:
            return node_list
        else:
            return []

    def add_node(self, node):
        """Add node to the network.

        Args:
            node (Node): Node to add to the network.
        """
        self.nodes_by_id[node.id] = node
        self.nodes_by_description[node.description] = node

    def add_link(self, link):
        """Add link to the network.

        Args:
            link (Link): Link to add to the network.
        """
        self.number_of_links += 1
        type_id = link.type.motif_link.motif_link_id

        set_of_link = self._get_set_of_type(link.type.motif_link)
        set_of_link.add(link.start)

        # Add the link for the start node
        if link.type.directed:
            reverse_set = self._get_set_of_type(link.type.inverse_motif_link)
            reverse_set.add(link.end)
        else:
            set_of_link.add(link.end)

        if type_id < len(link.start.neighbours_per_type):
            node_list = link.start.neighbours_per_type[type_id]
        else:
            node_list = []
            link.start.neighbours_per_type[type_id] = node_list

        if link.end not in node_list:
            link.start.neighbours_per_type[type_id].append(link.end)

        # Add the link for the end node
        if link.type.directed:
            type_id = link.type.inverse_motif_link.motif_link_id

        if type_id < len(link.end.neighbours_per_type):
            node_list = link.end.neighbours_per_type[type_id]
        else:
            node_list = []
            link.end.neighbours_per_type[type_id] = node_list

        if link.start not in node_list:
            link.end.neighbours_per_type[type_id].append(link.start)

    def finalize_network_construction(self):
        """Optimize network structure for further processing.
        """
        key_set = self.node_sets_departing_from_link.keys()
        for motif_link in key_set:
            nodes = list(self.node_sets_departing_from_link[motif_link])
            # This _should_ be equivalent to the overloaded `CompareTo()` in the Java
            # Node class when sorting
            self.nodes_with_link[motif_link] = sorted(nodes, key=cmp_to_key(node_id_compare))
        self.node_sets_departing_from_link = None

def node_id_compare(n1, n2):
    return n1.id - n2.id

def read_network_from_files(filenames, link_types):
    """Read in network structures from file(s).

    Example:
        A tab seperated text file containing with each file of the file being of
        the format "source_node\tdestination_node" e.g.
        1	2
        1	3
        1	4
        1	5
        1	6
        1	7
        1	8
        1	9
        1	10

    Args:
        filenames (list[str]): List of filenames to read in and create network from.
        link_types (list[LinkType]): Link types to assign to the links in the networks being generated.
    """
    network = Network()
    for i, filename in enumerate(filenames):
        link_type = link_types[i]
        links = 0
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                tab = line.index('\t')
                if tab <= 0 or '#' in line:
                    continue
                node_1 = line.strip('\n').split('\t')[0] + link_type.source_network
                node_2 = line.strip('\n').split('\t')[1] + link_type.destination_network

                if node_1 == node_2:
                    continue

                if node_1 in network.nodes_by_description:
                    origin = network.get_nodes_by_description(node_1)
                else:
                    origin = Node(description=node_1)
                    network.add_node(origin)

                if node_2 in network.nodes_by_description:
                    destination = network.get_nodes_by_description(node_2)
                else:
                    destination = Node(description=node_2)
                    network.add_node(destination)

                nodes = origin.neighbours_per_type[link_type.motif_link.motif_link_id]
                if nodes is not None and destination in nodes:
                    continue

                link = Link(origin, destination, link_type)

                links += 1
                network.add_link(link)

        for node in network.nodes_by_description.values():
            for node_set in range(len(node.neighbours_per_type)):
                if node.neighbours_per_type[node_set] is not None:
                    # This _should_ be equivalent to the overloaded `CompareTo()` in the Java
                    # Node class when sorting
                    node.neighbours_per_type[node_set].sort(key=cmp_to_key(node_id_compare))

        logger.info(f"Read: {filename} | Links: {links}")

    network.finalize_network_construction()

    logger.info(f"Number of Nodes: {len(network.nodes_by_description)}")
    logger.info(f"Number of Links: {network.number_of_links}")

    return network
