#!/usr/bin/env bash
# Authors:
#   Šimon Petrov (primary implementation, testing)
#   Kamila Riedlová (implementation review, extensions, testing)
#   Parts of the implementation were assisted by ChatGPT5.

# requirements: biopython, wsl,
# program runs in wsl. All paths must be in linux path format (e.g. /mnt/c/Users/...)

python3 /path/to/Docking_and_docking_preparation/01.1_fetch_structures-script.py \
  --dataset /path/to/Docking_and_docking_preparation/cryptobench-main-structures.json \
  --holo-dir /path/to/Docking_and_docking_preparation/Docking/holo_cif \
  --apo-dir /path/to/Docking_and_docking_preparation/Docking/apo_cif \
  #--write-holo-noligand \
  #--holo-noligand-dir /path/to/Docking_and_docking_preparation/Docking/holo_no_ligand_pdb
