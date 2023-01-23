# Indexed-based Subgraph Matching Algorithm with General Symmetries (ISMAGS)

## Overview

ISMAGS is a python software package based on the original [java software](https://github.com/biointec/ismags) written by the authors of this algorithm and described in this [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0097896).

ISMAGS finds all instances - up to symmetry - of a subgraph (called a **motif**) embedded in a larger graph (called a **network**).

## Installation

pip install -r requirements.txt

## Usage

If you are already familiar with the java-based version of ISMAGS, then this python-based version should be very easy to use. Besides a slight change in command line syntax the operation of this software is identical to its java-based version.

```python
python cli.py [-h] [-f FOLDER] -l LINK_TYPES -n NETWORKS -m MOTIF_DESCRIPTION -o OUTPUT
```

- **Folder** is the path containing network files.
- **LINK_TYPES** are separated by commas e.g. \"A u P P\" or \"A u P P,B d P P\"
- **NETWORKS** are separated by commas e.g. file1.txt or file1.txt,file2.txt
- **MOTIF_DESCRIPTION** e.g. AA0A00 is a 3-star where all the connections are of type 'A'
- **OUTPUT** is the output file name.
- Note: at the command line type 'python cli.py -h' to obtain the above information.

A brief description of syntax usage can be found in the 'doc' folder.

For more information contact Mark DeBonis - mjdebon-at-sandia-dot-gov