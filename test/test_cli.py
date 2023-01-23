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

import os
import pytest

def test_cli():
    # Test 1
    os.system("python ../cli/cli.py -f '../data/' -l \"A d t t\" -n 'graph1_Ad.txt' -m AA0A00 -o 'graph1.out'")
    os.remove("graph1.out")
    # Test 2
    os.system("python ../cli/cli.py -f '../data/' -l \"A d t t,B u t t\" -n 'graph2_Ad.txt,graph2_Bu.txt' -m AB0B00 -o 'graph2.out'")
    os.remove("graph2.out")
    # Test 3
    os.system("python ../cli/cli.py -f '../data/' -l \"A d t t,B d t t,C d t t,D d t t,E d t t,F d t t\" \
-n 'graph3_Ad.txt,graph3_Bd.txt,graph3_Cd.txt,graph3_Dd.txt,graph3_Ed.txt,graph3_Fd.txt' \
-m AB00C00F0000E0000D000 -o 'graph3.out'")
    os.remove("graph3.out")
    # Test 4
    os.system("python ../cli/cli.py -f '../data/' -l \"A d t t,B d t t,C d t t,D u t t\" \
-n 'graph4_Ad.txt,graph4_Bd.txt,graph4_Cd.txt,graph4_Du.txt' \
-m ABDC00 -o 'graph4.out'")
    os.remove("graph4.out")
    # Test 5
    os.system("python ../cli/cli.py -f '../data/' -l \"A d t t,B d t t,C d t t\" \
-n 'graph5_Ad.txt,graph5_Bd.txt,graph5_Cd.txt' \
-m AB0C00 -o 'graph5.out'")
    os.remove("graph5.out")