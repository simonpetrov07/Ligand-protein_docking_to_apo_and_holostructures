# Ligand-protein docking to apo and holo protein structures
Proteins are fundamental components of all living organisms. Their main functions include participating in cellular processes and forming larger cellular structures. In order to understand how this occurs, we first need to understand how proteins interact with small molecules, known as ligands. These interactions depend on the 3D structure of a protein. However, a protein can change its form, which in turn alters its interaction with the ligand. Understanding how changes in a protein's form affect ligand interactions is crucial for precisely understanding how proteins influence cellular processes and how we can use this knowledge to improve drug design.
## Usage
Before docking, you need to prepare the necessary files first. After that, you can carry out the docking and further analysis.
At the end of the preparation you should have a folder structure like this:

```text

Docking/
 ├─ holo_cif/                           # HOLO structures (E.g. 3rj9.cif)
 ├─ apo_cif/                            # APO structures (E.g. 1a4u.cif)
 ├─ ligands/                            # ligands (E.g. NAD_1a4u.pdb)
 ├─ docking_parameters/                 # JSON parameters for Vina (docking_parameters_*.json)
 ├─ cryptobench-main-structures.json    # dataset - dataset will not be provided, as it has not been published yet.
 ├─ Dockerfile                          # file needed for docking
 ├─ run_all.sh                          # file needed for docking
 ├─ run_docking.py                      # docking program 
 └─ Docking_output/                     # This folder doesn't need to be in the structure. It will be created automatically after the start of docking. 
```
### Files preparation
If you are using Windows OS, you first need to download WSL and Biopython..

First, run the set of commands from  [`01.0_fetch_structures_run.sh`](Docking_and_docking_preparation/01.0_fetch_structures_run.sh) in WSL console. You need to specify the path to the [`Docking`](Docking_and_docking_preparation/Docking) and to the [`Docking_and_docking_preparation`](Docking_and_docking_preparation) folder. This will create folders holo_cif and apo_cif, that will contain all the .cif structures from the dataset. These structures will be needed further for extraction of ligands and for docking itself.

Second, run the set of commands from [`02.0_dockbox_calculations_run_v2`](Docking_and_docking_preparation/02.0_dockbox_calculations_run_v2).You need to specify the path to the [`Docking`](Docking_and_docking_preparation/Docking) and to the [`Docking_and_docking_preparation`](Docking_and_docking_preparation) folder. This will create the docking_parameters folder that contains all the docking parameters needed for the further docking.

Third, run the set of commands from [`03.0_extract_ligands_run_v3`](Docking_and_docking_preparation/03.0_extract_ligands_run_v3).You need to specify the path to the [`Docking`](Docking_and_docking_preparation/Docking) and to the [`Docking_and_docking_preparation`](Docking_and_docking_preparation) folder. This will extract all the ligands from the holostructures and save them to the ligands folder. These ligands will be further used in docking.

### Docking
After the preparation you should have all the files in the structure above.
For the start of the docking, open WSL terminal in the Docking folder. 

To run the docking, you first need to build the image.
```console
docker build -t apo-holo-docking .
```
After the image is ready, you can start the docking using this command.
```console
docker run --rm -u "$(id -u):$(id -g)" -v "/path/to/Docking:/data" apo-holo-docking 
```
Docker will process all the .JSON docking_parameters in a row and conduct the docking.

### Output
At the end, you should have in `Docking_output` two subfolders: `APO` with results of docking to apo structures and `HOLO` with results of docking to holo structures. An output is a folder e.g. `/1a4u` that would contain `1a4u_output.pdbqt`. This file contains 9 models showing the binding affinity between the protein's binding site and the ligand (Example of outputs are in folder [`Output_examples`](Output_examples))

