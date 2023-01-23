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

import argparse
import logging
import os
import sys
sys.path.insert(1, '../')
import textwrap
from time import perf_counter
from algorithm.motif_finder import MotifFinder
from motifs.motif import create_motif, write_motifs
from network.link import LinkType
from network.network import read_network_from_files

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class CLI():
    """Creates a new CLI instance that will read in and parse command line arguments
        for generating networks and motifs.

    Attributes:
        folder (str): String containing the path where the input files(s) are located.
        link_types (str): String containing comma seperated link types.
            Note: The first link type will also correspond to the first network file.
        networks (str): String containing comma seperated network file names.
            Note: The first link type will also correspond to the first network file.
        motif_description (str): String containing the description of the motif in the network.
        output (str): String containing the path to desired output file.
        motif (Motif): Motif object that was generated from the given command line descriptors.

    Methods:
        write_motifs(motifs): Write instances of the motif to find in the network to a specified file.

    """

    def __init__(self):
        """Initialize a new CLI instance and parse the command line arguments.
        """
        parser_desciption = textwrap.dedent('''\
                            The Index-based Subgraph Matching Algorithm with General Symmetries
                            -------------------------------------------------------------------
                            Original Copyright: Copyright (c) 2013-2014 Maarten Houbraken
                            Sandia Copyright: Copyright (c) 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS)''')
        parser = argparse.ArgumentParser(description=parser_desciption, formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("-f", "--folder", dest="folder", help="Folder path containing network files", default=os.getcwd())
        required = parser.add_argument_group('required arguments')
        required.add_argument("-l", "--link-types", dest="link_types", help="Link types seperated by commas e.g. \"A u P P\" or \"A u P P,A d P P\"", required=True)
        required.add_argument("-n", "--networks", dest="networks", help="Network files seperated by commas e.g. file1.txt or file1.txt,file2.txt", required=True)
        required.add_argument("-m", "--motif", dest="motif_description", help="Motif description e.g. AA0A00", required=True)
        required.add_argument("-o", "--output", dest="output", help="Output file name", required=True)
        args = parser.parse_args()

        # Grab all of the CLI args and do some light preprocessing
        self.folder = args.folder
        self.link_types = args.link_types.split(',')
        self.link_types = [x.strip() for x in self.link_types]
        self.networks = args.networks.split(',')
        if self.folder is not None:
            if self.folder[-1] != os.path.sep:
                self.folder = self.folder + os.path.sep
            self.networks = [self.folder + x.strip() for x in self.networks]
        else:
            self.networks = [x.strip() for x in self.networks]
        self.motif_description = args.motif_description
        self.output = args.output

        # Walk through the link_types and networks and generate the corresponding internal structures
        link_types_list = []
        link_type_translation = {}
        logger.info('Creating link types...')
        for i, link_type in enumerate(self.link_types):
            link_type_char_list = link_type.split(' ')
            if len(link_type_char_list) < 4 or len(link_type_char_list) > 4:
                logger.error(f'link type `{link_type}` doesn\'t meet specification, ignoring given link type.')
                if len(self.link_types) == 1:
                    raise ValueError('No valid link types to process, exiting.')
                else:
                    continue
            directed = link_type_char_list[1] == 'd'
            if link_type_char_list[0] not in link_type_translation:
                t = LinkType(directed=directed, link_type_id=i, source_network=link_type_char_list[2], destination_network=link_type_char_list[3])
            else:
                t = link_type_translation[link_type_char_list[0]]
            link_types_list.append(t)
            link_type_translation[link_type_char_list[0]] = t

        logger.info('Reading in networks...')
        self.network = read_network_from_files(self.networks, link_types_list)

        logger.info('Creating motif data structure...')
        self.motif = create_motif(self.motif_description, link_type_translation)

def main():
    """Simple main to run CLI parsing and the ISMAGS algorithm.
    """
    cli = CLI()

    motif_finder = MotifFinder(cli.network)
    motifs = motif_finder.find_motif(cli.motif, False)
    write_motifs(motifs, cli.output)

if __name__ == "__main__":
    main()
