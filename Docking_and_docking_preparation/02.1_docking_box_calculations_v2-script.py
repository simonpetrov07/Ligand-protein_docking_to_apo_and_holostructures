#!/usr/bin/env python3
#
# Authors:
#   Student Name <student.email@example.com> (primary implementation, testing)
#   Kamila Riedlová <kamila.riedlova@matfyz.cuni.cz> (implementation review, adjustmen and extensions, testing)
#   Parts of the implementation were assisted by ChatGPT5.
#
# File: 02.1_dockbox_from_structures-script.py

"""
Compute AutoDock Vina box (center + size) directly from PDB/mmCIF using CryptoBench selections (apo/holo pocket residue lists + chains).
Output structure JSON parameter files consumable by DODO.


Usage
-----
python3 02.1_dockbox_from_structures-script.py \
  --dataset cryptobench-main-structures.json \
  --holo-dir ./Docking/holo_cif \
  --apo-dir  ./Docking/apo_cif \
  --out-dir  ./Docking/docking_parameters \
  --receptor-ext .cif \
  --ligand-ext .pdb \
  --padding 5.0 \
  --min-size 10.0 \
  --write-holo --write-apo \
  --chain-select dataset --min-found-frac-chain 0.5 \
  --resnum-mode auto --collapse-icodes \
  --ligand-name-mode apo \
  --lig-dir ./Docking/ligands \
  --require-ligand
 # --exhaustiveness 32 --num-modes 10 --energy-range 7 --seed -131415392

"""

from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Structure import Structure

# -------------------- collects finnished structures ------------------
finished_list = []
def collect_finnished(pdb_id):
    finished_list.append(pdb_id)

# ---------------------------- small utils ----------------------------
def log(msg: str) -> None:
    print(msg, flush=True)

# ---------------------------- geometry box ----------------------------
class Box:
    def __init__(self):
        self.xmin = math.inf; self.xmax = -math.inf
        self.ymin = math.inf; self.ymax = -math.inf
        self.zmin = math.inf; self.zmax = -math.inf
    def update(self, x: float, y: float, z: float):
        if x < self.xmin: self.xmin = x
        if x > self.xmax: self.xmax = x
        if y < self.ymin: self.ymin = y
        if y > self.ymax: self.ymax = y
        if z < self.zmin: self.zmin = z
        if z > self.zmax: self.zmax = z
    def valid(self) -> bool:
        return all(math.isfinite(v) for v in (self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax))
    def center(self) -> Tuple[float, float, float]:
        return (
            (self.xmax + self.xmin) / 2.0,
            (self.ymax + self.ymin) / 2.0,
            (self.zmax + self.zmin) / 2.0,
        )
    def size(self, padding: float, min_size: float, cubic: bool) -> Tuple[float, float, float]:
        sx = max((self.xmax - self.xmin) + padding, min_size)
        sy = max((self.ymax - self.ymin) + padding, min_size)
        sz = max((self.zmax - self.zmin) + padding, min_size)
        if cubic:
            m = max(sx, sy, sz)
            return (m, m, m)
        return (sx, sy, sz)

# ---------------------------- parsing helpers ----------------------------
def parse_selection_ids(ids: Iterable[str]) -> List[Tuple[Optional[str], str, Optional[str]]]:
    """Return list of (sel_chain, resseq_str, icode) tuples."""
    out: List[Tuple[Optional[str], str, Optional[str]]] = []
    for item in ids:
        sel_chain: Optional[str] = None
        res = item
        if "_" in item:
            sel_chain, res = item.split("_", 1)
        if res and res[-1:].isalpha():
            out.append((sel_chain, res[:-1], res[-1]))
        else:
            out.append((sel_chain, res, None))
    return out

def load_structure(path: Path, struct_id: str) -> Structure:
    if path.suffix.lower() == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(struct_id, str(path))

def residue_matches(residue, resseq: str, icode: Optional[str]) -> bool:
    """Match by (author) resseq + optional insertion code. Includes modified residues (MSE, CME, ...)."""
    het, ridx, ic = residue.get_id()
    try:
        target = int(resseq)
    except ValueError:
        return False
    if ridx != target:
        return False
    if icode is not None and (ic or "") != icode:
        return False
    return True

def collect_box(structure: Structure, chain_id: str, selection: Sequence[Tuple[Optional[str], str, Optional[str]]]) -> Tuple[Box, List[Tuple[str, Optional[str]]]]:
    """Aggregate box over atoms for the given chain and (resseq, icode) list.
    Returns (box, missing_residues_list_as[(resseq, icode)]).
    """
    wanted = {(r, i) for (_sch, r, i) in selection}
    missing = set(wanted)
    box = Box()

    # Locate target chain
    model = next(structure.get_models())
    chain = None
    for ch in model:
        if ch.id == chain_id:
            chain = ch; break
    if chain is None:
        raise ValueError(f"Chain '{chain_id}' not found in structure")

    for residue in chain.get_residues():
        for (resseq, icode) in list(missing):
            if residue_matches(residue, resseq, icode):
                for atom in residue.get_atoms():
                    x, y, z = atom.coord
                    box.update(float(x), float(y), float(z))
                if (resseq, icode) in missing:
                    missing.remove((resseq, icode))
                break

    def _sort_key(t):
        r, i = t
        try: rr = int(r)
        except: rr = 0
        return (rr, i or "")
    return box, sorted(list(missing), key=_sort_key)

# ---- numbering maps (author <-> label) & chain resolver from mmCIF ----
def build_resnum_maps_from_cif(cif_path: Path, chain_id: str):
    d = MMCIF2Dict(str(cif_path))
    auth_asym = d.get("_atom_site.auth_asym_id", [])
    auth_seq  = d.get("_atom_site.auth_seq_id", [])
    label_asym= d.get("_atom_site.label_asym_id", [])
    label_seq = d.get("_atom_site.label_seq_id", [])
    icode     = d.get("_atom_site.pdbx_PDB_ins_code", [])
    n = min(len(auth_asym), len(auth_seq), len(label_asym), len(label_seq), len(icode))

    auth2label, label2auth = {}, {}
    for i in range(n):
        if str(label_asym[i]) != str(chain_id) and str(auth_asym[i]) != str(chain_id):
            continue
        aseq = auth_seq[i]; lseq = label_seq[i]; ic = icode[i]
        if aseq == "?" or lseq == "?":
            continue
        ic = ("" if ic in ("?", ".") else str(ic))
        key_auth = (int(aseq), ic)   # (author_seq, icode)
        val_lab  = int(lseq)         # label_seq
        auth2label[key_auth] = val_lab
        if val_lab not in label2auth:
            label2auth[val_lab] = key_auth
    return auth2label, label2auth

def build_chain_maps_from_cif(cif_path: Path):
    d = MMCIF2Dict(str(cif_path))
    label_asym = d.get("_atom_site.label_asym_id", [])
    auth_asym  = d.get("_atom_site.auth_asym_id", [])
    n = min(len(label_asym), len(auth_asym))
    pairs, seen = [], set()
    for i in range(n):
        la = str(label_asym[i]); aa = str(auth_asym[i])
        key = (la, aa)
        if key not in seen:
            seen.add(key); pairs.append(key)
    return pairs  # list[(label_asym_id, auth_asym_id)]

def resolve_chain_id(structure: Structure, cif_path: Path, target_chain: str) -> str:
    # direct match
    model = next(structure.get_models())
    for ch in model:
        if ch.id == target_chain:
            return target_chain
    # CIF mapping
    if cif_path.suffix.lower() == ".cif":
        pairs = build_chain_maps_from_cif(cif_path)
        label_ids = {la for (la, aa) in pairs}
        auth_ids  = {aa for (la, aa) in pairs}
        if target_chain in label_ids:
            return target_chain
        if target_chain in auth_ids:
            for (la, aa) in pairs:
                if aa == target_chain:
                    return la
    # fallback
    return target_chain

def list_chain_ids(structure: Structure) -> List[str]:
    model = next(structure.get_models())
    return [ch.id for ch in model]

def expand_missing_with_icodes(selection: Sequence[Tuple[Optional[str], str, Optional[str]]], missing: List[Tuple[str, Optional[str]]]):
    tried = ["A","B","C","D","E"]
    expanded = list(selection)
    missing_set = set(missing)
    for (_sch, res, ic) in selection:
        if (res, ic) in missing_set and ic is None:
            for t in tried:
                expanded.append((None, res, t))
    return expanded

def collapse_icodes(miss_disp_list: List[str]) -> List[str]:
    out, seen = [], set()
    for s in miss_disp_list:
        base = s[:-1] if s and s[-1].isalpha() and s[:-1].isdigit() else s
        if base not in seen:
            seen.add(base); out.append(base)
    return out

def try_collect(struct: Structure,
                chain: str,
                selection: Sequence[Tuple[Optional[str], str, Optional[str]]],
                mode: str,
                cif_path: Path,
                min_found_frac: float,
                try_icode: bool,
                collapse: bool):
    original_len = max(1, len(selection))

    def _collect(sel):
        box, missing = collect_box(struct, chain, sel)
        frac = 1.0 - (len(missing) / original_len)
        miss_disp = [f"{r}{(i or '')}" for (r, i) in missing]
        if collapse:
            miss_disp = collapse_icodes(miss_disp)
        return box, missing, miss_disp, frac

    # maps for LABEL→AUTHOR (if CIF)
    a2l, l2a = {}, {}
    if cif_path.suffix.lower() == ".cif":
        try:
            a2l, l2a = build_resnum_maps_from_cif(cif_path, chain)
        except Exception:
            a2l, l2a = {}, {}

    def map_label_to_author(sel):
        mapped = []
        for (_sch, res, ic) in sel:
            try:
                l = int(res)
            except Exception:
                mapped.append((None, res, ic)); continue
            ai = l2a.get(l, None)  # (author_seq, icode)
            if ai is not None:
                aseq, aic = ai
                mapped.append((None, str(aseq), aic if aic != "" else None))
            else:
                mapped.append((None, res, ic))
        return mapped

    # explicit modes
    if mode == "auth":
        sel = selection
        box, missing, miss_disp, frac = _collect(sel)
        if try_icode and missing:
            sel2 = expand_missing_with_icodes(selection, missing)
            box, missing2, miss_disp, frac = _collect(sel2)
        return box, miss_disp, frac, "auth"

    if mode == "label":
        sel = map_label_to_author(selection)
        box, missing, miss_disp, frac = _collect(sel)
        if try_icode and missing:
            sel2 = expand_missing_with_icodes(sel, missing)
            box, missing2, miss_disp, frac = _collect(sel2)
        return box, miss_disp, frac, "label"

    # AUTO: compute BOTH paths and pick the better one
    sel_a = selection
    box_a, miss_a, disp_a, frac_a = _collect(sel_a)
    if try_icode and miss_a:
        sel_a2 = expand_missing_with_icodes(sel_a, [(m, None) for m in miss_a])
        box_a, miss_a, disp_a, frac_a = _collect(sel_a2)

    best = (box_a, disp_a, frac_a, "auth")

    if cif_path.suffix.lower() == ".cif":
        sel_l = map_label_to_author(selection)
        box_l, miss_l, disp_l, frac_l = _collect(sel_l)
        if try_icode and miss_l:
            sel_l2 = expand_missing_with_icodes(sel_l, miss_l)
            box_l, miss_l, disp_l, frac_l = _collect(sel_l2)
        if frac_l > frac_a:
            best = (box_l, disp_l, frac_l, "label")

    return best

# ---------------------------- ligand naming & validation ----------------------------
def choose_ligand_basename(
    lig_code: str,
    apo_id: str,
    holo_id: str,
    mode: str,
    lig_dir: Optional[Path],
    lig_ext: str,
    require: bool,
    context_tag: str
) -> Optional[str]:
    """
    Returns the chosen ligand basename (without extension) OR None when require=True and file is missing.
    Validates existence if lig_dir is provided. Logs what is chosen and why.
    """
    candidates = []
    if mode == "apo":
        candidates = [(f"{lig_code}_{apo_id}", "apo")]
    elif mode == "holo":
        candidates = [(f"{lig_code}_{holo_id}", "holo")]
    else:  # auto → prefer apo name, fall back to holo name
        candidates = [(f"{lig_code}_{apo_id}", "apo"), (f"{lig_code}_{holo_id}", "holo")]

    if lig_dir is None:
        # No validation possible
        name = candidates[0][0]
        log(f"[{context_tag}] ligand-name-mode={mode} → {name}{lig_ext} (no --lig-dir, not validated)")
        return name

    for base, tag in candidates:
        p = lig_dir / f"{base}{lig_ext}"
        if p.exists():
            log(f"[{context_tag}] ligand-name-mode={mode} → chose {base}{lig_ext} (found in {lig_dir})")
            return base

    # none found
    expected = ", ".join([f"{base}{lig_ext}" for base, _ in candidates])
    log(f"[{context_tag}] WARNING: no ligand found in {lig_dir}. Tried: {expected}")
    if require:
        log(f"[{context_tag}] --require-ligand is set → skipping JSON emission.")
        return None
    # still return the first naming pattern (for downstream tools that might provide it later)
    name = candidates[0][0]
    log(f"[{context_tag}] proceeding with non-validated name: {name}{lig_ext}")
    return name

# ---------------------------- JSON builder ----------------------------
def make_param_json(
    receptor_basename: str,
    ligand_basename: str,
    box: Box,
    padding: float,
    min_size: float,
    receptor_ext: str,
    ligand_ext: str,
    cubic: bool,
    extras: Optional[Dict] = None,
   # pdb_id: str = "",
) -> Dict:
    #if pdb_id in finished_list:
    #    count = finished_list.count(pdb_id) + 1 
    cx, cy, cz = box.center()
    sx, sy, sz = box.size(padding=padding, min_size=min_size, cubic=cubic)
    data = {
        "receptor": f"{receptor_basename}{receptor_ext}",
        "ligand": f"{ligand_basename}{ligand_ext}",
        "output": f"{receptor_basename}_output.pdbqt",
        #"output": f"{receptor_basename}_{count}_output.pdbqt",
        "center": {"x": round(cx, 4), "y": round(cy, 4), "z": round(cz, 4)},
        "size": {"x": round(sx, 3), "y": round(sy, 3), "z": round(sz, 3)},
    }
    if extras:
        data.update(extras)
    return data

# ---------------------------- main ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="Path to cryptobench JSON")
    ap.add_argument("--holo-dir", required=True, help="Folder with HOLO structures (.cif or .pdb)")
    ap.add_argument("--apo-dir", required=True, help="Folder with APO structures (.cif or .pdb)")
    ap.add_argument("--out-dir", required=True, help="Output folder for docking parameter JSONs")
    ap.add_argument("--receptor-ext", default=".cif", choices=[".cif", ".pdb"], help="Receptor file extension")
    ap.add_argument("--ligand-ext", default=".pdb", choices=[".pdb", ".sdf", ".smi", ".pdbqt"], help="Ligand file extension for the JSON field")
    ap.add_argument("--padding", type=float, default=5.0, help="Extra Å added to each axis range before min-size clamp")
    ap.add_argument("--min-size", type=float, default=15.0, help="Minimum per-axis size in Å")
    ap.add_argument("--cubic", action="store_true", help="Force cubic box (max axis applied to all)")
    ap.add_argument("--write-holo", action="store_true", help="Emit HOLO docking JSONs")
    ap.add_argument("--write-apo", action="store_true", help="Emit APO docking JSONs")
    ap.add_argument("--exhaustiveness", type=int, default=None)
    ap.add_argument("--num-modes", type=int, default=None)
    ap.add_argument("--energy-range", type=int, default=None)
    ap.add_argument("--seed", type=int, default=None)
    # Numbering & logs
    ap.add_argument("--resnum-mode", default="auto", choices=["auto","auth","label"],
                    help="Residue numbering: author, label, or auto-choose the better one")
    ap.add_argument("--min-found-frac", type=float, default=0.5,
                    help="Compatibility threshold (auto anyway compares both)")
    ap.add_argument("--try-icode", action="store_true",
                    help="If residue N is missing, try N with insertion codes (A,B,...)")
    ap.add_argument("--collapse-icodes", action="store_true",
                    help="Collapse A..E variants in missing-residue logs")
    # Chain selection
    ap.add_argument("--chain-select", default="dataset", choices=["dataset","auto"],
                    help="Use dataset chain as-is or automatically pick the chain with best coverage")
    ap.add_argument("--min-found-frac-chain", type=float, default=0.5,
                    help="Minimum coverage to accept an auto-selected chain")
    # Ligand naming + validation
    ap.add_argument("--ligand-name-mode", default="apo", choices=["apo","holo","auto"],
                    help="How to name ligand files in JSONs: LIG_<apo_id>.* or LIG_<holo_id>.* (auto prefers apo)")
    ap.add_argument("--lig-dir", type=str, default=None,
                    help="Directory with ligand files to validate existence (e.g., ./Docking/ligands)")
    ap.add_argument("--require-ligand", action="store_true",
                    help="If set, skip JSON emission when selected ligand file is missing in --lig-dir")
    args = ap.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    holo_dir = Path(args.holo_dir)
    apo_dir = Path(args.apo_dir)
    lig_dir = Path(args.lig_dir).resolve() if args.lig_dir else None

    dataset: Dict[str, Dict] = json.loads(Path(args.dataset).read_text(encoding="utf-8"))

    extras = {}
    if args.exhaustiveness is not None: extras["exhaustiveness"] = args.exhaustiveness
    if args.num_modes is not None: extras["num_modes"] = args.num_modes
    if args.energy_range is not None: extras["energy_range"] = args.energy_range
    if args.seed is not None: extras["seed"] = args.seed

    def process_one(kind: str, pdb_id: str, chain_hint: str, sel_ids: Sequence[str], struct_dir: Path,
                    lig_code: str, apo_id: str, holo_id: str):
        context = f"{kind} {pdb_id}"
        struct_path = struct_dir / f"{pdb_id}{args.receptor_ext}"
        if not struct_path.exists():
            log(f"[{kind}] {pdb_id}: missing file {struct_path}")
            return

        try:
            struct = load_structure(struct_path, pdb_id)
            # parse selection
            sel = parse_selection_ids(sel_ids)

            # chain choice
            if args.chain_select == "auto":
                candidate_chains = list_chain_ids(struct)
                best = (None, -1.0, None)
                for ch_id in candidate_chains:
                    box, missing_disp, frac, used = try_collect(
                        struct, ch_id, sel, args.resnum_mode, struct_path,
                        args.min_found_frac, args.try_icode, args.collapse_icodes
                    )
                    if frac > best[1]:
                        best = (ch_id, frac, (box, missing_disp, used))
                chosen_chain, best_frac, res = best
                if chosen_chain is None or best_frac < args.min_found_frac_chain:
                    chosen_chain = resolve_chain_id(struct, struct_path, chain_hint)
                    box, missing_disp, frac, used = try_collect(
                        struct, chosen_chain, sel, args.resnum_mode, struct_path,
                        args.min_found_frac, args.try_icode, args.collapse_icodes
                    )
                else:
                    box, missing_disp, used = res
                final_chain = chosen_chain
            else:
                final_chain = resolve_chain_id(struct, struct_path, chain_hint)
                box, missing_disp, frac, used = try_collect(
                    struct, final_chain, sel, args.resnum_mode, struct_path,
                    args.min_found_frac, args.try_icode, args.collapse_icodes
                )

            # verbose logs
            if not box.valid():
                raise ValueError("no atoms aggregated for pocket")
            if missing_disp:
                log(f"[{kind}] {pdb_id}: missing residues in chain {final_chain} ({used}): {','.join(missing_disp)}")
            cx, cy, cz = box.center(); sx, sy, sz = box.size(args.padding, args.min_size, args.cubic)
            log(f"[{kind}] {pdb_id}: chain={final_chain} mode={used} cover={frac:.2%} "
                f"center=({cx:.3f},{cy:.3f},{cz:.3f}) size=({sx:.3f},{sy:.3f},{sz:.3f})")

            # ligand basename selection + validation
            lig_base = choose_ligand_basename(
                lig_code=lig_code,
                apo_id=apo_id,
                holo_id=holo_id,
                mode=args.ligand_name_mode,
                lig_dir=lig_dir,
                lig_ext=args.ligand_ext,
                require=args.require_ligand,
                context_tag=f"{kind} {pdb_id}"
            )
            if lig_base is None:
                return  # required but missing

            # write JSON
            param = make_param_json(
                receptor_basename=pdb_id,
                ligand_basename=lig_base,
                box=box,
                padding=args.padding,
                min_size=args.min_size,
                receptor_ext=args.receptor_ext,
                ligand_ext=args.ligand_ext,
                cubic=args.cubic,
                extras=extras if extras else None,
            )
            
            if pdb_id in finished_list: # avoid overwrites
                count = finished_list.count(pdb_id) + 1
                out_path = out_dir / f"docking_parameters_{pdb_id}_{count}.json"
                out_path.write_text(json.dumps(param, indent=2), encoding="utf-8")
            else:
                out_path = out_dir / f"docking_parameters_{pdb_id}.json"
                out_path.write_text(json.dumps(param, indent=2), encoding="utf-8")

            collect_finnished(pdb_id)
        except Exception as e:
            log(f"[{kind}] {pdb_id}: {e}")

    # -------- main loop over dataset records --------
    for apo_pdb_id, rec in dataset.items():
        lig_code = rec.get("ligand", "LIG").strip()
        holo_pdb_id = rec.get("holo_pdb_id", "").strip()
        holo_chain = rec.get("holo_chain", "").strip()
        apo_chain  = rec.get("apo_chain", "").strip()

        if args.write_holo and holo_pdb_id:
            process_one(
                kind="HOLO",
                pdb_id=holo_pdb_id,
                chain_hint=holo_chain,
                sel_ids=rec.get("holo_pocket_selection", []),
                struct_dir=holo_dir,
                lig_code=lig_code,
                apo_id=apo_pdb_id,
                holo_id=holo_pdb_id,
            )

        if args.write_apo and apo_pdb_id:
            process_one(
                kind="APO",
                pdb_id=apo_pdb_id,
                chain_hint=apo_chain,
                sel_ids=rec.get("apo_pocket_selection", []),
                struct_dir=apo_dir,
                lig_code=lig_code,
                apo_id=apo_pdb_id,
                holo_id=holo_pdb_id,
            )

if __name__ == "__main__":
    main()
