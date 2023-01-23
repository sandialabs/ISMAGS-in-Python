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

class MotifLink:
    """Represents a motif link with all related attributes.

    Attributes:
        link_type (LinkType): Type of link for the given motif.
        directed (bool): Whether the link is directed or not.
        number_of_link_ids (int): The total number of link ids in the motif link
        motif_link_id (int): The id of the link
        link_id_to_link (dict): Mapping of the link id and the link

    """

    NUMBER_OF_LINK_IDS = 0

    def __init__(self, link_type=None,
                        directed=False):
        """Initialize a link within a motif

        Keyword Args:
            link_type (LinkType): Type of link for the given motif. Defaults to None.
            directed (bool): Whether the link is directed or not. Defaults to False.
        """
        self.link_type = link_type
        self.directed = directed
        self.motif_link_id = MotifLink.NUMBER_OF_LINK_IDS
        MotifLink.NUMBER_OF_LINK_IDS += 1
        self.link_id_to_motif_link = {self.motif_link_id: self}

