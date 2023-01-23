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
import math

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Motif:
    """Creates a new motif without any edges.

    Attributes:
        number_of_motif_nodes (int): Number of nodes in a given motif.
        links (list[list[MotifLink()]]): The links that the motif is comprised of.
        initial_connections (dict): Initial unoptimized connections between nodes in the motif.
        final_connections (list[list]): Final unoptimized connections between nodes in the motif.

    Methods:
        add_motif_link(start_node, end_node, link_type): Add a motif link between two nodes.
        finalize_motif(): Optimize motif structure for further processing.

    """

    def __init__(self, number_of_motif_nodes=0):
        """Initialize a new motif without any edges.

        Keyword Args:
            number_of_motif_nodes (int): Number of nodes in a given motif. Defaults to 0.
        """
        self.number_of_motif_nodes = number_of_motif_nodes
        self.links = [self.__init_links(self.number_of_motif_nodes) for _ in range(self.number_of_motif_nodes)]
        self.initial_connections = {}
        self.final_connections = None

        self.link_types = set()

        self.description = ""

    def __str__(self):
        return self.description

    def __init_links(self, count):
        return [None for _ in range(count)]

    def add_motif_link(self, start_node, end_node, link_type):
        """Add a motif link between two nodes.

        Args:
            start_node (int): Starting node of motif link.
            end_node (int): Ending node of motif link.
            link_type (LinkType): Type of link between nodes.
        """
        self.links[start_node][end_node] = link_type.motif_link
        self.links[end_node][start_node] = link_type.inverse_motif_link

        if start_node not in self.initial_connections:
            self.initial_connections[start_node] = []

        self.initial_connections[start_node].append(end_node)

        if end_node not in self.initial_connections:
            self.initial_connections[end_node] = []

        self.initial_connections[end_node].append(start_node)
        self.link_types.add(link_type.link_type_id)

    def finalize_motif(self):
        """Finalizes motif construction by optimizing internal data structures for
            further use in the algorithm
        """
        self.final_connections = [list() for _ in range(self.number_of_motif_nodes)]
        for i in range(self.number_of_motif_nodes):
            arr = self.initial_connections[i]
            self.final_connections[i] = [0 for _ in arr]
            final_links = [0 for _ in range(len(arr))]
            for j in range(len(arr)):
                self.final_connections[i][j] = arr[j]
                final_links[j] = self.links[i][arr[j]]
            self.links[i] = final_links
        self.initial_connections = None

def create_motif(motif_description, link_type_translation):
    """Given a motif description and link type translation value
        generate a motif object to find in the graph.

    Args:
        motif_description (str): String description of the motif. For example, "AA00AA".
        link_type_translation (dict): Dictionary of motif description characters and their corresponding
            LinkType(). For example, {'A':LinkType()}.

    Raises:
        ValueError: If the len of the motif description doesn't equal the number of nodes required to form motif.

    Returns:
        Motif(): Generated motif to find in the network.
    """
    motif_length = len(motif_description)
    number_of_nodes = int(math.ceil(math.sqrt(2 * motif_length)))
    if motif_length != number_of_nodes * (number_of_nodes - 1) / 2:
        raise ValueError(f'Motif description `{motif_description}` has invalid length!')
    counter = 0
    motif = Motif(number_of_motif_nodes=number_of_nodes)
    motif.description = motif_description

    for i in range(1, number_of_nodes):
        for j in range(i):
            character = motif_description[counter]
            counter += 1
            if character == '0':
                continue
            link_type = link_type_translation[character]
            if character.isupper():
                motif.add_motif_link(j, i, link_type)
            else:
                motif.add_motif_link(i, j, link_type)

    motif.finalize_motif()
    return motif

def write_motifs(motifs, output):
    """Write instances of the motif found in a given network(s) to a specified file.

    Args:
        motifs (set): Motif instances to write to a file.
        output (str): File to write motif instances to.
    """
    logger.info(f'Writing motif instances to `{output}`')
    with open(output, 'w') as f:
        for motif in motifs:
            f.write(f"{str(motif)}\n")
