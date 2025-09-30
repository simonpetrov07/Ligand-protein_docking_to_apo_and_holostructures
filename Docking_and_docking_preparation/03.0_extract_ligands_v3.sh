#!/usr/bin/env bash
# Authors:
#   Šimon Petrov (primary implementation, testing)
#   Kamila Riedlová (implementation review, adjustmen and extensions, testing)
#   Parts of the implementation were assisted by ChatGPT5.

# program runs in wsl. All paths must be in linux path format (e.g. /mnt/c/Users/...)

python3 /path/to/Docking_and_docking_preparation/03.1_extract_ligands-script_v3.py \
  --dataset /path/to/Docking_and_docking_preparation/cryptobench-main-structures.json \
  --holo-dir /path/to/Docking_and_docking_preparation/Docking/holo_cif \
  --out-dir  /path/to/Docking_and_docking_preparation/Docking/ligands \
  --receptor-ext .cif
