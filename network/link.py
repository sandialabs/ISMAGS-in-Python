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

class LinkType:
    """Class representing the the type of link in the graph.

    Attributes:
        LINK_TYPES (dict): Contains a mapping between link type's IDs and their type.
        directed (bool): Weather the link between nodes is directed or not.
        link_type_id (int): ID of the given link type.
        source_network (str): Source network description.
        destination_network (str): Destination network description.
        motif_link (MotifLink): Link to corresponding motif/subgraph and its attributes.
        inverse_motif_link (MotifLink): Link to corresponding inversed motif/subgraph and its attributes.

    """

    LINK_TYPES = {}

    def __init__(self, directed=False,
                        link_type_id=0,
                        source_network="",
                        destination_network=""):
        """Initializes a graph edge's (link's) attributes/type.

        Keyword Args:
            directed (bool): Weather the link between nodes is directed or not. Defaults to False.
            link_type_id (int): ID of the given link type. Defaults to 0.
            source_network (str): Source network description. Defaults to "".
            destination_network (str): Destination network description. Defaults to "".
        """

        self.directed = directed
        self.link_type_id = link_type_id
        self.source_network = source_network
        self.destination_network = destination_network
        self.motif_link = None
        self.inverse_motif_link = None

        if directed:
            self.motif_link = MotifLink(link_type=self, directed=True)
            self.inverse_motif_link = MotifLink(link_type=self, directed=False)
        else:
            self.motif_link = MotifLink(link_type=self, directed=True)
            self.inverse_motif_link = self.motif_link

        LinkType.LINK_TYPES[self.link_type_id] = self

    def __len__(self):
        return len(LinkType.LINK_TYPES)

class Link:
    """Class representsing graph edge between two nodes.

    Attributes:
        id (int): ID of the given graph edge link.
        start (Node): Starting node of the graph edge link.
        end (Node): Ending node the graph edge link.
        type (LinkType): Type of link between nodes.
        NEXT_AVAILABLE_ID (int): Next ID that can be assigned to a Link

    """

    NEXT_AVAILABLE_ID = -1

    def __init__(self, start, end, type):
        """Initializes a graph edge (link) between two nodes.

        Keyword Args:
            start (Node): Starting node of the graph edge link. Defaults to a default Node object.
            end (Node): Ending node the graph edge link. Defaults to a default Node object.
            type (LinkType): Type of link between nodes. Defaults to a default LinkType object.
        """
        Link.NEXT_AVAILABLE_ID += 1
        self.id = Link.NEXT_AVAILABLE_ID
        self.start = start
        self.end = end
        self.type = type