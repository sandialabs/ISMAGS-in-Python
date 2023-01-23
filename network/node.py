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

from motifs.motif_link import MotifLink

class Node:
    """Class representing a node in graph.

    Attributes:
        used (bool): Weather or not the node has been used.
        id (int): Node ID.
        description (int): Node description.
        NEXT_AVAILABLE_ID (int): Next ID that can be assigned to a Node

    """

    NEXT_AVAILABLE_ID = 0

    def __init__(self, used=False, description=""):
        """Initalize a node in graph

        Keyword Args:
            used (bool): Weather or not the node has been used. Defaults to False.
            description (str): Node description. Defaults to "".
        """
        self.used = used
        Node.NEXT_AVAILABLE_ID += 1
        self.id = Node.NEXT_AVAILABLE_ID
        self.description = description

        self.number_of_motif_link_types = MotifLink.NUMBER_OF_LINK_IDS
        self.neighbours_per_type = []
        for _ in range(self.number_of_motif_link_types):
            self.neighbours_per_type.append([])

    def __str__(self):
        return self.description

    def __lt__(self, node):
        """Compare the ID numbers of two nodes.

        Args:
            node (Node): Node to compare against.

        Raises:
            TypeError: If the node being compared against isn't a Node() instance.

        Returns:
            bool: True or false depending on if the given Node ID is less than another Node's ID.
        """
        if isinstance(self, node.__class__):
            return self.id < node.id
        else:
            raise TypeError

    def __sub__(self, node):
        """Calculate the difference between node IDs

        Args:
            node (Node): Node to compare against.

        Raises:
            TypeError: If the node being compared against isn't a Node() instance.

        Returns:
            int: The difference between the node's ID numbers.
        """
        if isinstance(self, node.__class__):
            return self.id - node.id
        else:
            raise TypeError