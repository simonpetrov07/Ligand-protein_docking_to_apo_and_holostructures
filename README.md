# Ligand-protein_docking_to_apo_and_holostructures

## Usage
Before the docking you need to prepare files first. After that, you can carry out the docking and further analysis.
At the end of the preparation you should have a folder structure like this:

```text
Docking/
 ├─ holo_cif/                           # HOLO structures (E.g. 3rj9.cif)
 ├─ apo_cif/                            # APO structures (E.g. 1a4u.cif)
 ├─ ligands/                            # ligands (E.g. NAD_1a4u.pdb)
 ├─ docking_parameters/                 # JSON parameters for Vina (docking_parameters_*.json)
 ├─ cryptobench-main-structures.json    # dataset 
 └─ Docking_output/                     # This folder doesn't need to be in the structure. It will be created automaticaly after the start of docking. 
```
### Files preparation
If you are using Windows OS, it is need to first download WSL and Biopython. Then you need to copy the set of commands from  
`01.0_fetch_structures_run.sh` and paste them in wsl console. 
