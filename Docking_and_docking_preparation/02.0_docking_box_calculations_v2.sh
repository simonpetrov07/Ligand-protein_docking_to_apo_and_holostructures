#!/usr/bin/env bash
# Authors:
#   Šimon Petrov (primary implementation, testing)
#   Kamila Riedlová (implementation review, adjustmen and extensions, testing)
#   Parts of the implementation were assisted by ChatGPT5.

# program runs in wsl. All paths must be in linux path format (e.g. /mnt/c/Users/...)

python3 /path/to/Docking_and_docking_preparation/02.1_docking_box_calculations_v2-script.py \
  --dataset /path/to/Docking_and_docking_preparation/cryptobench-main-structures.json \
  --holo-dir /path/to/Docking_and_docking_preparation/Docking/holo_cif \
  --apo-dir /path/to/Docking_and_docking_preparation/Docking/apo_cif \
  --out-dir  /path/to/Docking_and_docking_preparation/Docking/docking_parameters \
  --receptor-ext .cif \
  --ligand-ext .pdb \
  --padding 5.0 \
  --min-size 10.0 \
  --write-holo --write-apo \
  --chain-select dataset --min-found-frac-chain 0.5 \
  --resnum-mode auto --collapse-icodes \
  --ligand-name-mode apo \
 #--lig-dir ./Docking/ligands \
 #--require-ligand
 # --exhaustiveness 32 --num-modes 9 --energy-range 7 --seed -123456
