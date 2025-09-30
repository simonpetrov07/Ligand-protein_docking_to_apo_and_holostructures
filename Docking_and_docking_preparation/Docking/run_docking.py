#!/usr/bin/env python3
# Created by: Kamila Riedlová (parts of the implementation were taken over from DODO), some edits were assisted by ChatGPT5
# File: run_docking.py
"""
Docking runner with robust receptor/ligand preparation and resume markers.

Features:
- Reads ONE docking-parameter JSON (receptor path, ligand filename, box center/size, etc.).
- Finds receptor either in /data/holo_cif or /data/apo_cif.
- Filters receptor to protein-only and prepares receptor PDBQT via MGLTools.
- Prepares the ligand to PDBQT (MGLTools → Open Babel fallback → RDKit-from-SMILES as last resort).
- Accepts ligand fallback naming (e.g., NAD_1a4u.pdb not found → try NAD_<holo>.pdb based on dataset).
- Writes outputs into /data/Docking_output/{HOLO|APO}/{ID}/ and creates a .ok marker for resume.
- If output and .ok marker already exist, the job is skipped (useful for interrupted batches).

This script is intended to be called by run_all.sh for many JSONs, or directly for one JSON.
"""

import os, sys, json, shutil, subprocess, time
from pathlib import Path
from typing import Optional, List

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select
from rdkit import Chem
from rdkit.Chem import AllChem

# ---- Container paths (mounted from the host) ----
DATA_ROOT   = Path("/data")
HOLO_DIR    = DATA_ROOT / "holo_cif"      # HOLO receptors (.cif or .pdb)
APO_DIR     = DATA_ROOT / "apo_cif"       # APO receptors  (.cif or .pdb)
LIG_DIR     = DATA_ROOT / "ligands"       # ligands (NAD_1a4u.pdb, etc.)
OUT_ROOT    = DATA_ROOT / "Docking_output"
TMP_ROOT    = Path("/tmp")                

# MGLTools utilities used for receptor/ligand preparation
MGLTOOLS_BIN_PATH = "/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/"

# Soft timeouts; in Docker we often run with use_timeout=False to avoid unwanted kills
PREPARATION_TIMEOUT_SECONDS = 600
DOCKING_TIMEOUT_SECONDS     = 900

def log(msg: str):
    """Timestamped console log."""
    print(time.strftime("[%H:%M:%S]"), msg, flush=True)

def is_valid_pdbqt(pdbqt_file: Path) -> bool:
    """
    Quick sanity check that a PDBQT is not empty/corrupt:
    Find at least one ATOM/HETATM with a parseable (non-zero-ish) X coordinate.
    """
    try:
        with open(pdbqt_file, "r") as f:
            for line in f:
                if line.startswith(("ATOM","HETATM")):
                    try:
                        x = float(line[30:38].strip())
                        if abs(x) > 1e-3:
                            return True
                    except Exception:
                        continue
    except Exception:
        return False
    return False

class ProteinOnlySelect(Select):
    """
    Biopython selector: keep standard amino acids (+ common alt histidines/cystines).
    All HETATMs get converted to ATOM later in filter_pdb; here we just filter residues.
    """
    def accept_residue(self, residue):
        return residue.get_resname().strip().upper() in {
            "ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS",
            "LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL",
            "TRP","TYR","HID","HSP","HIE","HIP","CYX","CSS"
        }

def filter_cif_file_biopython(input_cif: Path, output_cif: Path):
    """
    Write a protein-only copy of the input mmCIF.
    We also preserve a few mmCIF header categories (_entry.id, _cell, _symmetry)
    so that downstream MGLTools doesn't complain about missing blocks.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("ID", str(input_cif))
    io = MMCIFIO(); io.set_structure(structure)
    tmp = output_cif.with_suffix(output_cif.suffix + ".tmp")
    io.save(str(tmp), select=ProteinOnlySelect())

    needed = ["_entry.id", "_cell", "_symmetry"]; saved = {}
    with open(input_cif, "r") as f:
        for line in f:
            for k in needed:
                if line.startswith(k):
                    saved.setdefault(k, []).append(line.rstrip("\n"))

    with open(tmp, "r") as infile, open(output_cif, "w") as out:
        # mmCIF writer emits a header; skip first two lines defensively
        try:
            next(infile); next(infile)
        except StopIteration:
            pass
        out.write("data_ID\n# \n")
        for k, vals in saved.items():
            for ln in vals: out.write(ln + "\n")
            out.write("# \n")
        out.writelines(infile)
    tmp.unlink(missing_ok=True)

def filter_pdb(input_pdb: Path, output_pdb: Path):
    """
    Very simple PDB filter: keep ATOM; once ATOMs seen, convert later HETATMs to ATOM
    (so protein-like modified residues remain). Stop on first TER/END.
    """
    with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
        writing_hetatm = False
        for line in fin:
            if line.startswith("ATOM"):
                fout.write(line); writing_hetatm = True
            elif line.startswith("HETATM") and writing_hetatm:
                fout.write("ATOM  " + line[6:])
            elif line.startswith(("TER","END")):
                break

def prepare_receptor(in_filtered: Path, out_pdbqt: Path, use_timeout=True):
    """Run MGLTools prepare_receptor4.py on the filtered receptor."""
    log(f"Preparing receptor (MGLTools) → {out_pdbqt.name}")
    cmd = [
        "/usr/local/bin/python2.7",
        os.path.join(MGLTOOLS_BIN_PATH, "prepare_receptor4.py"),
        "-A","hydrogens","-U","nphs _lps_waters_deleteAltB",
        "-r", str(in_filtered), "-o", str(out_pdbqt),
    ]
    subprocess.run(cmd, check=True,
                   timeout=PREPARATION_TIMEOUT_SECONDS if use_timeout else None)

def prepare_ligand_mgltools(ligand_pdb_local: Path, out_pdbqt: Path, work_dir: Path, use_timeout=True):
    """
    Prepare ligand via MGLTools. 
    """
    cmd = [
        "/usr/local/bin/python2.7",
        os.path.join(MGLTOOLS_BIN_PATH, "prepare_ligand4.py"),
        "-l", ligand_pdb_local.name,
        "-o", str(out_pdbqt),
    ]
    cwd = os.getcwd(); os.chdir(work_dir)
    try:
        log(f"Preparing ligand (MGLTools) → {out_pdbqt.name}")
        subprocess.run(cmd, check=True,
                       timeout=PREPARATION_TIMEOUT_SECONDS if use_timeout else None)
    finally:
        os.chdir(cwd)

def prepare_ligand_obabel(ligand_pdb: Path, out_pdbqt: Path, use_timeout=True):
    """
    Fallback via Open Babel: PDB → PDBQT with hydrogens and Gasteiger charges.
    """
    log(f"Preparing ligand (Open Babel fallback) → {out_pdbqt.name}")
    subprocess.run(
        ["obabel","-ipdb",str(ligand_pdb),
         "-opdbqt","-O",str(out_pdbqt),
         "-h","--partialcharge","gasteiger"],
        check=True,
        timeout=PREPARATION_TIMEOUT_SECONDS if use_timeout else None
    )

def rdkit_fallback_from_smiles(smiles_file: Path, out_pdbqt: Path):
    """
    Last-resort fallback for .smi ligands:
    - build 3D with ETKDG
    - UFF optimization
    - Gasteiger partial charges
    - write a minimal PDBQT
    """
    smi = smiles_file.read_text(encoding="utf-8").strip().splitlines()[0]
    mol = Chem.MolFromSmiles(smi)
    if mol is None: raise ValueError(f"Invalid SMILES: {smi}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDG(); params.randomSeed = 42
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params)
    if not cids: raise RuntimeError("EmbedMultipleConfs failed")
    energies = AllChem.UFFOptimizeMoleculeConfs(mol)
    best = min(range(len(energies)), key=lambda i: energies[i][1])
    conf = mol.GetConformer(best)
    mol.RemoveAllConformers(); mol.AddConformer(conf, assignId=True)
    AllChem.UFFOptimizeMolecule(mol); AllChem.ComputeGasteigerCharges(mol)
    lines = ["REMARK   Generated by RDKit fallback", "ROOT"]
    for i, atom in enumerate(mol.GetAtoms(), start=1):
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        charge = float(atom.GetProp("_GasteigerCharge")) if atom.HasProp("_GasteigerCharge") else 0.0
        el = atom.GetSymbol()
        lines.append(f"ATOM  {i:5d}  {el:<4s}{'LIG':>3s} A{1:4d}    {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00 {charge:6.3f}          {el:>2s}")
    lines += ["ENDROOT","TORSDOF 0","END"]
    out_pdbqt.write_text("\n".join(lines), encoding="utf-8")

def read_json(json_file: Path):
    """
    Load docking parameters JSON:
      {
        "receptor": "1a4u.cif",
        "ligand": "NAD_1a4u.pdb",
        "output": "1a4u_output.pdbqt",
        "center": {"x":..., "y":..., "z":...},
        "size":   {"x":..., "y":..., "z":...},
        ... optional: exhaustiveness, num_modes, energy_range ...
      }
    """
    data = json.loads(json_file.read_text(encoding="utf-8"))
    receptor_name = data["receptor"]
    ligand_name   = data["ligand"]
    out_name      = data["output"]
    center        = data["center"]; size = data["size"]
    extra = {
        "exhaustiveness": str(data.get("exhaustiveness", 32)),
        "num_modes":      str(data.get("num_modes", 9)),
        "energy_range":   str(data.get("energy_range", 7)),
    }
    return receptor_name, ligand_name, out_name, center, size, extra

def locate_receptor(receptor_name: str):
    """
    Try HOLO first, then APO. Return (path, kind, id),
    where kind ∈ {"HOLO","APO"} and id is the PDB code without extension.
    """
    path_h = HOLO_DIR / receptor_name
    path_a = APO_DIR / receptor_name
    if path_h.exists(): return path_h, "HOLO", Path(receptor_name).stem
    if path_a.exists(): return path_a, "APO",  Path(receptor_name).stem
    raise FileNotFoundError(f"Receptor {receptor_name} not found in {HOLO_DIR} or {APO_DIR}")

def _dataset_record_for_apo(apo_id: str):
    """
    Optional dataset lookup (cryptobench-main-structures.json at /data):
    used to map APO id → corresponding HOLO id for ligand-name fallback.
    """
    ds_path = DATA_ROOT / "cryptobench-main-structures.json"
    if not ds_path.exists():
        return None
    try:
        ds = json.loads(ds_path.read_text(encoding="utf-8"))
        return ds.get(apo_id)
    except Exception:
        return None

def locate_ligand(ligand_name: str) -> Optional[Path]:
    """
    Locate ligand file inside /data/ligands.
    - First try exact match.
    - If name looks like LIG_<APOID>.ext and exact file is missing, try LIG_<HOLOID>.<same/ext alternatives>
      based on the dataset mapping.
    - Finally, try any LIG_*.{pdb,sdf,mol2,cif,smi} with same LIG prefix (WARN).
    Return a Path or None.
    """
    # 1) exact
    p = LIG_DIR / ligand_name
    if p.exists():
        return p

    # 2) fallback via APO→HOLO mapping
    lig_code = None; apo_id = None; ext = None
    try:
        lig_code, rest = ligand_name.split("_", 1)
        apo_id, ext = rest.rsplit(".", 1)
        ext = "." + ext
    except ValueError:
        lig_code = ligand_name.split("_", 1)[0] if "_" in ligand_name else ligand_name

    if lig_code and apo_id:
        rec = _dataset_record_for_apo(apo_id)
        if rec:
            holo_id = rec.get("holo_pdb_id", "").strip()
            if holo_id:
                exts: List[str] = [ext] if ext else []
                for alt in [".pdb", ".sdf", ".mol2", ".cif", ".smi"]:
                    if alt and alt not in exts:
                        exts.append(alt)
                for e in exts:
                    cand = LIG_DIR / f"{lig_code}_{holo_id}{e}"
                    if cand.exists():
                        print(f"[WARN] Using ligand fallback: {cand.name} for expected {ligand_name}")
                        return cand

    # 3) last resort: any candidate with same ligand code prefix
    if not lig_code:
        lig_code = ligand_name.split("_", 1)[0] if "_" in ligand_name else ligand_name
    candidates = []
    for e in (".pdb",".sdf",".mol2",".cif",".smi"):
        candidates += list(LIG_DIR.glob(f"{lig_code}_*{e}"))
    if candidates:
        cand = sorted(candidates)[0]
        print(f"[WARN] Using ligand candidate: {cand.name} for expected {ligand_name}")
        return cand

    return None

def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)

def run_vina(receptor_pdbqt: Path, ligand_pdbqt: Path, out_pdbqt: Path,
             center, size, extra, log_file: Path, use_timeout=True):
    """Run AutoDock Vina with given receptor/ligand PDBQT and box."""
    cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand",   str(ligand_pdbqt),
        "--out",      str(out_pdbqt),
        "--center_x", str(center["x"]),
        "--center_y", str(center["y"]),
        "--center_z", str(center["z"]),
        "--size_x",   str(size["x"]),
        "--size_y",   str(size["y"]),
        "--size_z",   str(size["z"]),
        "--exhaustiveness", extra["exhaustiveness"],
        "--num_modes",      extra["num_modes"],
        "--energy_range",   extra["energy_range"],
    ]
    log(f"Running Vina → {out_pdbqt.name}")
    with open(log_file, "w") as lf:
        subprocess.run(
            cmd, check=True, stdout=lf, stderr=subprocess.STDOUT,
            timeout=DOCKING_TIMEOUT_SECONDS if use_timeout else None
        )

def run_one(json_path: Path, use_timeout=False):
    """
    Main worker for a single JSON:
    - read JSON,
    - locate receptor & ligand,
    - prepare receptor/ligand,
    - run Vina,
    - write .ok marker on success,
    - copy inputs (receptor/ligand PDBQT) into the output folder.
    """
    receptor_name, ligand_name, out_name, center, size, extra = read_json(json_path)

    # Locate receptor (HOLO or APO)
    try:
        receptor_in, kind, rec_id = locate_receptor(receptor_name)
    except FileNotFoundError as e:
        print(f"[SKIP] {json_path.name}: {e}")
        return 0

    # Locate ligand (with fallbacks)
    ligand_in = locate_ligand(ligand_name)
    if ligand_in is None:
        print(f"[SKIP] {json_path.name}: ligand {ligand_name} not found (tried HOLO fallback + candidates).")
        return 0

    # Output dir + resume marker
    out_dir = OUT_ROOT / kind / rec_id; ensure_dir(out_dir)
    out_pdbqt = out_dir / out_name
    ok_marker = out_pdbqt.with_suffix(out_pdbqt.suffix + ".ok")
    log_file  = out_dir / (Path(out_name).stem + ".log")

    # Resume: when both output and .ok marker are present → skip
    if out_pdbqt.exists() and ok_marker.exists():
        print(f"[SKIP] {json_path.name}: already completed ({out_pdbqt.name}).")
        return 0

    # Scratch per receptor id
    tmp = TMP_ROOT / rec_id; ensure_dir(tmp)
    rec_base = tmp / Path(receptor_in).stem

    # Filter receptor (mmCIF vs PDB)
    if receptor_in.suffix.lower() == ".cif":
        rec_filtered = rec_base.with_name(rec_base.name + "_filtered.cif")
        log(f"Filtering receptor CIF → {rec_filtered.name}")
        filter_cif_file_biopython(receptor_in, rec_filtered)
    else:
        rec_filtered = rec_base.with_name(rec_base.name + "_filtered.pdb")
        log(f"Filtering receptor PDB → {rec_filtered.name}")
        filter_pdb(receptor_in, rec_filtered)

    # Receptor PDBQT via MGLTools
    receptor_pdbqt = rec_base.with_suffix(".pdbqt")
    try:
        prepare_receptor(rec_filtered, receptor_pdbqt, use_timeout)
    except subprocess.CalledProcessError as e:
        print(f"[SKIP] {json_path.name}: receptor prep failed ({e}).")
        return 0

    if not is_valid_pdbqt(receptor_pdbqt):
        print(f"[SKIP] {json_path.name}: receptor PDBQT invalid.")
        return 0
    log(f"Receptor prepared: {receptor_pdbqt.name}")

    # Ligand (convert to PDB if needed)
    lig_base = TMP_ROOT / Path(ligand_in).stem
    if ligand_in.suffix.lower() == ".pdb":
        ligand_pdb = ligand_in
    else:
        ligand_pdb = lig_base.with_suffix(".pdb")
        fmt = ligand_in.suffix.lower().lstrip(".")
        log(f"Converting ligand → PDB via Open Babel: {ligand_in.name}")
        try:
            subprocess.run(
                ["obabel","-i",fmt,str(ligand_in),"-o","pdb","-O",str(ligand_pdb),"--gen3d"],
                check=True, timeout=PREPARATION_TIMEOUT_SECONDS if use_timeout else None
            )
        except subprocess.CalledProcessError as e:
            # If original was .smi, try RDKit path and dock from generated PDBQT
            if ligand_in.suffix.lower() == ".smi":
                log(f"OBabel PDB conversion failed: {e} → RDKit from SMILES")
                ligand_pdbqt = lig_base.with_suffix(".pdbqt")
                rdkit_fallback_from_smiles(ligand_in, ligand_pdbqt)
                if not is_valid_pdbqt(ligand_pdbqt):
                    print(f"[SKIP] {json_path.name}: ligand PDBQT invalid after RDKit fallback.")
                    return 0
                try:
                    run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt, center, size, extra, log_file, use_timeout)
                except subprocess.CalledProcessError as e:
                    print(f"[SKIP] {json_path.name}: vina failed (code {e.returncode}). See log: {log_file.name}")
                    return 0
                for src in (receptor_pdbqt, ligand_pdbqt):
                    dst = out_dir / src.name
                    if not dst.exists(): shutil.copy2(src, dst)
                ok_marker.write_text(time.strftime("OK %Y-%m-%d %H:%M:%S\n"), encoding="utf-8")
                log(f"[{kind}] {rec_id}: done → {out_pdbqt.name} (log: {log_file.name})")
                return 0
            else:
                print(f"[SKIP] {json_path.name}: ligand conversion to PDB failed ({e}).")
                return 0

    # Ligand PDBQT (MGLTools → OBabel fallback)
    ligand_pdbqt = lig_base.with_suffix(".pdbqt")
    local_lig = TMP_ROOT / ligand_pdb.name
    if ligand_pdb.resolve() != local_lig.resolve():
        shutil.copy2(ligand_pdb, local_lig)

    mgltools_ok = True
    try:
        prepare_ligand_mgltools(local_lig, ligand_pdbqt, TMP_ROOT, use_timeout)
    except Exception as e:
        log(f"MGLTools ligand prep failed: {e}")
        mgltools_ok = False

    if mgltools_ok and not is_valid_pdbqt(ligand_pdbqt):
        log("MGLTools produced invalid ligand PDBQT → switching to Open Babel fallback")
        mgltools_ok = False

    if not mgltools_ok:
        try:
            prepare_ligand_obabel(local_lig, ligand_pdbqt, use_timeout)
        except Exception as e2:
            if ligand_in.suffix.lower() == ".smi":
                log(f"Open Babel fallback failed: {e2} → RDKit from SMILES")
                rdkit_fallback_from_smiles(ligand_in, ligand_pdbqt)
            else:
                print(f"[SKIP] {json_path.name}: ligand prep failed ({e2}).")
                return 0

    if not is_valid_pdbqt(ligand_pdbqt):
        print(f"[SKIP] {json_path.name}: ligand PDBQT invalid after fallbacks.")
        return 0
    log(f"Ligand prepared: {ligand_pdbqt.name}")

    # Docking (AutoDock Vina)
    try:
        run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt, center, size, extra, log_file, use_timeout)
    except subprocess.CalledProcessError as e:
        print(f"[SKIP] {json_path.name}: vina failed (code {e.returncode}). See log: {log_file.name}")
        return 0

    # Archive inputs (for audit / reproducibility)
    for src in (receptor_pdbqt, ligand_pdbqt):
        dst = out_dir / src.name
        if not dst.exists(): shutil.copy2(src, dst)

    # Success → write resume marker
    ok_marker.write_text(time.strftime("OK %Y-%m-%d %H:%M:%S\n"), encoding="utf-8")
    log(f"[{kind}] {rec_id}: done → {out_pdbqt.name} (log: {log_file.name})")
    return 0

if __name__ == "__main__":
    # Script expects exactly ONE JSON path argument (usually provided by run_all.sh)
    if len(sys.argv) != 2:
        print("Usage: python3 run_docking.py /data/docking_parameters/<file>.json")
        sys.exit(1)
    run_one(Path(sys.argv[1]), use_timeout=False)
