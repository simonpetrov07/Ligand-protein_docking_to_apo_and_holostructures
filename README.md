# Ligand-protein docking to apo and holo protein structures
Proteins are fundamental components of all living organisms. Their main functions include participating in cellular processes and forming larger cellular structures. In order to understand how this occurs, we first need to understand how proteins interact with small molecules, known as ligands. These interactions depend on the 3D structure of a protein. However, a protein can change its form, which in turn alters its interaction with the ligand. Understanding how changes in a protein's form affect ligand interactions is crucial for precisely understanding how proteins influence cellular processes and how we can use this knowledge to improve drug design.
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
If you are using Windows OS, it's needed to first download WSL and Biopython.

 Firstly run the set of commands from  `01.0_fetch_structures_run.sh` in wsl console. This will create folders holo_cif and apo_cif, that will contain all the .cif structures from the dataset. These structures will be needed further for extraction of ligands and for docking itself.

Secondly run the set of commands from `02.0_dockbox_calculations_run_v2`. This will create the docking_parameters folder that contains all the docking parameters needed for the further docking.

Thirdly run the set of commands from `03.0_extract_ligands_run_v3`. This will extract all the ligands from the holostructures and save them to the ligands folder.
