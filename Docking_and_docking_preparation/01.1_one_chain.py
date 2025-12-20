#!/usr/bin/env python3
#
# Authors:
#   Šimon Petrov (writing, primary implementation, testing)
#   Parts of the implementation were assisted by ChatGPT5.
#
# Change the paths to your own directories before running.
#
# This code can download only 100 pairs at once. To download more, run multiple times or modify the batch_size variable. 
# The batch size can be adjusted though it is recommended to keep it at 100 to avoid exceeding programms limits.
#
# File: 01.0_one_chain.py

from pymol import cmd
import json
import os

with open(r"C:/Users/path/to/Docking_and_docking_preparation/cryptobench-main-structures.json", "r") as a: # change path here
    dict1 = json.load(a)

structures = list(dict1.keys())

batch_size = 100

BASE_OUT = r"C:/Users/path/to/Docking_and_docking_preparation/Docking" # change path here

HOLO_DIR = os.path.join(BASE_OUT, "holo_cif")
APO_DIR = os.path.join(BASE_OUT, "apo_cif")

os.makedirs(HOLO_DIR, exist_ok=True)
os.makedirs(APO_DIR, exist_ok=True)

def main(structures, dict1):
    processed = 0
    for structure in structures:
        if processed >= batch_size:
          cmd.quit()

        holo_structure = dict1[structure]["holo_pdb_id"]
        holo_chain = dict1[structure]["holo_chain"]
        apo_chain = dict1[structure]["apo_chain"]
        apo_structure = structure

        holo_out_path = os.path.join(HOLO_DIR, f"{holo_structure}.cif")
        apo_out_path =  os.path.join(APO_DIR, f"{apo_structure}.cif")

        did_holo = holo_download(holo_structure, holo_chain, holo_out_path)
        did_apo = apo_download(apo_structure, apo_chain, apo_out_path)

        if did_holo or did_apo:
          processed += 1

def holo_download(holo_structure, holo_chain, holo_out_path):
      if os.path.exists(holo_out_path):
        print(f"Skipping holo {holo_structure} (exists)")
        return False

      try:
        cmd.reinitialize()
        cmd.fetch(holo_structure)
        cmd.remove(f"not chain {holo_chain} and model {holo_structure}")
        cmd.save(holo_out_path)
        return True
      except Exception as e:
        print(f"Warning: Could not fetch {holo_structure}: {e}")
        return False
      
def apo_download(apo_structure, apo_chain, apo_out_path):
      if os.path.exists(apo_out_path):
        print(f"Skipping apo {apo_structure} (exists)")
        return False
      
      try:
        cmd.reinitialize()
        cmd.fetch(apo_structure)
        cmd.remove(f"not chain {apo_chain} and model {apo_structure}")
        cmd.save(apo_out_path)
        return True
      
      except Exception as e:
        print(f"Warning: Could not fetch {apo_structure}: {e}")
        return False

main(structures, dict1)