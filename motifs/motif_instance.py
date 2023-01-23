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

class MotifInstance:
    """Represents a motif instance in the graph.

    Attributes:
        mapping (list): Mapping of all nodes in the given motif.

    """

    def __init__(self, mapping=None):
        """Initialize a motif instance.

        If no mapping is provided then the mapping is an empty list
        otherwise the mapping is a copy of the provided mapping.

        Args:
            mapping (list): Mapping of all nodes in the given motif.  Defaults to None.
        """
        if mapping is None:
            self.mapping = []
        else:
            self.mapping = mapping.copy()

    def __str__(self):
        """String representation of the motif.

        Returns:
            A string representation of the motif of the form "1;2;3;4".
        """
        descriptions = [str(node.description) for node in self.mapping]
        return ';'.join(descriptions)
