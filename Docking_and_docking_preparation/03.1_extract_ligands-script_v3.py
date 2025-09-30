#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Authors:
#   Šimon Petrov (primary implementation, testing)
#   Kamila Riedlová (design/review, adjustments & extensions, testing)
#   Parts of the implementation were assisted by ChatGPT5.
#
# File: 03.1_extract_ligands-script.py
"""
Extract HOLO ligands from mmCIF/PDB into standalone PDB files named <LIG>_<apo_id>.pdb.

Usage:
python3 03.1_extract_ligands-script_v3.py \
  --dataset cryptobench-main-structures.json \
  --holo-dir ./Docking/holo_cif \
  --out-dir ./Docking/ligands \
  --receptor-ext .cif
"""


from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Structure import Structure

# ---------------- helpers ----------------
def load_structure(path: Path, struct_id: str) -> Structure:
    if path.suffix.lower() == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(struct_id, str(path))

def list_chain_ids(structure: Structure) -> List[str]:
    model = next(structure.get_models())
    return [ch.id for ch in model]

def build_chain_maps_from_cif(cif_path: Path):
    """Return list of (label_asym_id, auth_asym_id) unique pairs."""
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
    return pairs

def resolve_chain_id(structure: Structure, cif_path: Path, target_chain: str) -> str:
    """Try exact match in Biopython chain ids; else map label/auth ids via CIF."""
    for ch in next(structure.get_models()):
        if ch.id == target_chain:
            return target_chain
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
    return target_chain  # keep; we’ll remap on write if needed

def build_resnum_maps_from_cif(cif_path: Path, chain_id: str):
    """
    Map label_seq_id → (auth_seq_id, icode) for a given chain.
    Skip '.'/'?' entries.
    """
    d = MMCIF2Dict(str(cif_path))
    auth_asym = d.get("_atom_site.auth_asym_id", [])
    auth_seq  = d.get("_atom_site.auth_seq_id", [])
    label_asym= d.get("_atom_site.label_asym_id", [])
    label_seq = d.get("_atom_site.label_seq_id", [])
    icode     = d.get("_atom_site.pdbx_PDB_ins_code", [])
    n = min(len(auth_asym), len(auth_seq), len(label_asym), len(label_seq), len(icode))
    label2auth: Dict[int, Tuple[int, str]] = {}
    for i in range(n):
        if str(label_asym[i]) != str(chain_id) and str(auth_asym[i]) != str(chain_id):
            continue
        aseq = str(auth_seq[i]); lseq = str(label_seq[i]); ic = str(icode[i])
        if aseq in ("?", ".") or lseq in ("?", "."):
            continue
        ic = ("" if ic in ("?", ".") else ic)
        try:
            l = int(lseq); a = int(aseq)
        except ValueError:
            continue
        if l not in label2auth:
            label2auth[l] = (a, ic)
    return label2auth

def parse_selection_ids(ids: Sequence[str]) -> List[Tuple[str, Optional[str]]]:
    """Return [(resseq_str, icode)] for pocket residues; ignores chain part."""
    out: List[Tuple[str, Optional[str]]] = []
    for item in ids:
        res = item.split("_", 1)[1] if "_" in item else item
        if res and res[-1:].isalpha():
            out.append((res[:-1], res[-1]))
        else:
            out.append((res, None))
    return out

def pocket_centroid(structure: Structure, chain_id: str, selection: Sequence[Tuple[str, Optional[str]]]) -> Optional[Tuple[float,float,float]]:
    """Compute centroid of listed residues (auth numbering + optional icode) on a chain."""
    model = next(structure.get_models())
    chain = None
    for ch in model:
        if ch.id == chain_id:
            chain = ch; break
    if chain is None: return None
    total = [0.0, 0.0, 0.0]; count = 0
    wanted = {(r, i) for (r, i) in selection}
    for res in chain.get_residues():
        het, ridx, ic = res.get_id()
        key = (str(ridx), (ic or "") or None)
        if key in wanted:
            for atom in res.get_atoms():
                x, y, z = atom.coord
                total[0] += float(x); total[1] += float(y); total[2] += float(z)
                count += 1
    if count == 0: return None
    return (total[0]/count, total[1]/count, total[2]/count)

def residue_center(res) -> Tuple[float,float,float]:
    total = [0.0, 0.0, 0.0]; n = 0
    for a in res.get_atoms():
        x, y, z = a.coord
        total[0] += float(x); total[1] += float(y); total[2] += float(z); n += 1
    return (total[0]/n, total[1]/n, total[2]/n) if n else (math.nan, math.nan, math.nan)

class OnlyChosenLigand(Select):
    """Select exactly one residue by chain+resseq+icode+resname (resname may be '*ANY*' to ignore)."""
    def __init__(self, chain_id: str, resseq: int, icode: Optional[str], resname: Optional[str]):
        self.chain_id = chain_id
        self.resseq   = resseq
        self.icode    = (icode or "")
        self.resname  = (resname.strip().upper() if resname else None)
    def accept_chain(self, chain):
        return 1 if chain.id == self.chain_id else 0
    def accept_residue(self, residue):
        het, ridx, ic = residue.get_id()
        if ridx != self.resseq: return 0
        if (ic or "") != self.icode: return 0
        if self.resname is not None:
            if residue.get_resname().strip().upper() != self.resname:
                return 0
        return 1
    def accept_atom(self, atom):
        return 1

def pick_free_chain_id(occupied: List[str]) -> str:
    """Pick a free single-char chain id not in occupied."""
    candidates = [chr(c) for c in range(ord('A'), ord('Z')+1)] + [str(d) for d in range(10)]
    for c in candidates:
        if c not in occupied:
            return c
    # last resort
    return 'Z'

def save_residue_as_pdb(structure: Structure, chain, resseq: int, icode: str, resname: Optional[str], out_path: Path) -> None:
    """
    Save single residue to PDB via PDBIO, ensuring chain id is 1-char and unique.
    Temporarily remap chain.id if needed; restore afterwards.
    """
    io = PDBIO(); io.set_structure(structure)
    model = next(structure.get_models())
    occupied = [ch.id for ch in model]
    orig_chain_id = chain.id
    temp_id = orig_chain_id
    try:
        if len(orig_chain_id) != 1 or sum(1 for ch in model if ch.id == orig_chain_id) != 1:
            # Need a unique 1-char id
            temp_id = pick_free_chain_id(occupied)
            chain.id = temp_id
        selector = OnlyChosenLigand(temp_id, resseq, icode, resname)
        io.save(str(out_path), select=selector)
    finally:
        chain.id = orig_chain_id  # restore

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--holo-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--receptor-ext", default=".cif", choices=[".cif", ".pdb"])
    args = ap.parse_args()

    ds = json.loads(Path(args.dataset).read_text(encoding="utf-8"))
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    for apo_id, rec in ds.items():
        lig_code  = rec.get("ligand", "").strip().upper()
        lig_idx_s = rec.get("ligand_index", "").strip()  # author index (string)
        lig_chain = rec.get("ligand_chain", "").strip()
        holo_id   = rec.get("holo_pdb_id")
        sel_holo  = rec.get("holo_pocket_selection", [])

        if not (lig_code and lig_chain and holo_id):
            print(f"[LIG] {apo_id}: missing ligand metadata, skipping")
            continue

        holo_path = Path(args.holo_dir) / f"{holo_id}{args.receptor_ext}"
        if not holo_path.exists():
            print(f"[LIG] {apo_id}: missing HOLO file {holo_path}")
            continue

        try:
            struct = load_structure(holo_path, holo_id)
            model = next(struct.get_models())

            # resolve intended chain id (accept label/auth)
            resolved_chain_id = resolve_chain_id(struct, holo_path, lig_chain)

            def get_chain(cid: str):
                for ch in model:
                    if ch.id == cid:
                        return ch
                return None

            chain = get_chain(resolved_chain_id)
            chain_pool: List = [chain] if chain is not None else list(model)

            # 1) candidates by resname (strict match)
            candidates: List[Tuple[object,int,str,object]] = []  # (res, ridx, icode, chain_obj)
            for ch in chain_pool:
                for res in ch.get_residues():
                    het, ridx, ic = res.get_id()
                    if res.get_resname().strip().upper() == lig_code:
                        candidates.append((res, ridx, (ic or ""), ch))

            chosen: Optional[Tuple[object,int,str,object]] = None

            # 2) if given, try by ligand_index (author numbering) — even if resname didn’t match
            lig_idx: Optional[int] = None
            if lig_idx_s:
                try:
                    lig_idx = int(lig_idx_s)
                except Exception:
                    lig_idx = None

            if lig_idx is not None and chosen is None:
                # prefer declared chain; else any chain
                try_first = chain_pool
                for ch in try_first:
                    for res in ch.get_residues():
                        het, ridx, ic = res.get_id()
                        if ridx == lig_idx:
                            chosen = (res, ridx, (ic or ""), ch); break
                    if chosen: break

                # If nothing yet and CIF, try label→author map across chains (when dataset uses label ids)
                if chosen is None and holo_path.suffix.lower() == ".cif":
                    for ch in model:
                        l2a = build_resnum_maps_from_cif(holo_path, ch.id)
                        if lig_idx in l2a:
                            aseq, aic = l2a[lig_idx]
                            for res in ch.get_residues():
                                het, ridx, ic = res.get_id()
                                if ridx == aseq and (ic or "") == (aic or ""):
                                    chosen = (res, ridx, (ic or ""), ch); break
                        if chosen: break

            # 3) if still no pick, but have candidates by resname → choose nearest to pocket centroid
            if chosen is None and candidates:
                ref_chain_id = resolved_chain_id if get_chain(resolved_chain_id) else candidates[0][3].id
                pc = pocket_centroid(struct, ref_chain_id, parse_selection_ids(sel_holo))
                if pc is not None:
                    px, py, pz = pc
                    best = (None, float("inf"))
                    for (res, ridx, ic, ch) in candidates:
                        cx, cy, cz = residue_center(res)
                        d2 = (cx - px)**2 + (cy - py)**2 + (cz - pz)**2
                        if d2 < best[1]:
                            best = ((res, ridx, ic, ch), d2)
                    chosen, _ = best
                else:
                    chosen = candidates[0]

            # 4) if still nothing → report and continue
            if chosen is None:
                print(f"[LIG] {apo_id}: ligand {lig_code} not found in chain {resolved_chain_id}")
                continue

            res, ridx, ic, ch = chosen
            out_path = out_dir / f"{lig_code}_{apo_id}.pdb"
            # When chosen by index (no reliable resname), allow any resname for writing
            resname_for_write = lig_code if candidates else None
            save_residue_as_pdb(struct, ch, ridx, ic, resname_for_write, out_path)
            print(f"[LIG] {apo_id}: wrote {out_path.name} (chain {ch.id}, res {ridx}{ic or ''})")

        except Exception as e:
            print(f"[LIG] {apo_id}: ERROR {e}")

if __name__ == "__main__":
    main()
