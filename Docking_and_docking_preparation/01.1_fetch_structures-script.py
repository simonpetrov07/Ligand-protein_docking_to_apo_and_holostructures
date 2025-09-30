#!/usr/bin/env python3
#
# Authors:
#   Šimon Petrov (primary implementation, testing)
#   Kamila Riedlová (implementation review,  extensions, testing)
#   Parts of the implementation were assisted by ChatGPT5.
#
# File: 01.1_fetch_structures-script.py

"""
Download CIFs structures (APO + HOLO) from RCSB according to a CryptoBench dataset JSON  and optionally export HOLO-without-ligand PDBs.

Usage
-----
python3 01.1_fetch_structures-script.py \
  --dataset cryptobench-main-structures.json \
  --apo-dir Docking/apo_cif \
  --holo-dir Docking/holo_cif \
  --holo-noligand-dir Docking/holo_no_ligand_pdb \
  --write-holo-noligand
"""

from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
from typing import Dict, Optional

import requests
from Bio.PDB import MMCIFParser, PDBIO, Select

RCSB_FILE_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
UA = {"User-Agent": "cryptobench-downloader/1.0"}


def http_get_file(url: str, out_path: Path, retries: int = 3, timeout: int = 30) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    last_err: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, headers=UA, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                tmp = out_path.with_suffix(out_path.suffix + ".part")
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1024 * 128):
                        if chunk:
                            fh.write(chunk)
                tmp.replace(out_path)
                return
        except Exception as e:
            last_err = e
            if attempt < retries:
                time.sleep(2 * attempt)
            else:
                raise last_err


def download_cif(pdb_id: str, out_dir: Path) -> Path:
    url = RCSB_FILE_URL.format(pdb_id=pdb_id)
    out_path = out_dir / f"{pdb_id}.cif"
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    http_get_file(url, out_path)
    return out_path


class ExcludeLigand(Select):
    def __init__(self, chain_id: str, resseq: int, icode: str = ""):
        super().__init__()
        self.chain_id = chain_id
        self.resseq = resseq
        self.icode = icode
    def accept_residue(self, residue):
        het, ridx, ic = residue.get_id()  # residue id tuple: (hetflag, resseq, icode)
        # Keep all residues that are not the target ligand residue
        if residue.parent.id != self.chain_id:
            return True
        if ridx != self.resseq:
            return True
        if (ic or "") != (self.icode or ""):
            return True
        # This is the ligand residue to drop
        return False


def write_holo_without_ligand(cif_path: Path, out_pdb: Path, ligand_chain: str, ligand_index: str) -> None:
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure(cif_path.stem, str(cif_path))
    try:
        resseq = int(ligand_index.rstrip().rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"))
        icode = ligand_index[len(str(resseq)):]  # supports optional insertion code suffix
    except Exception:
        # Fallback: treat whole as int, no icode
        resseq = int(ligand_index)
        icode = ""
    io = PDBIO()
    io.set_structure(struct)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(out_pdb), select=ExcludeLigand(ligand_chain, resseq, icode))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="Path to cryptobench JSON (structures dict)")
    ap.add_argument("--apo-dir", required=True, help="Output dir for APO CIFs")
    ap.add_argument("--holo-dir", required=True, help="Output dir for HOLO CIFs")
    ap.add_argument("--write-holo-noligand", action="store_true", help="Also write HOLO-without-ligand as PDB")
    ap.add_argument("--holo-noligand-dir", default="holo_no_ligand_pdb", help="Output dir for HOLO-no-ligand PDBs")
    args = ap.parse_args()

    apo_dir = Path(args.apo_dir); apo_dir.mkdir(parents=True, exist_ok=True)
    holo_dir = Path(args.holo_dir); holo_dir.mkdir(parents=True, exist_ok=True)
    nolig_dir = Path(args.holo_noligand_dir)

    with open(args.dataset, "r", encoding="utf-8") as fh:
        dataset: Dict[str, Dict] = json.load(fh)

    for apo_pdb_id, rec in dataset.items():
        # APO
        try:
            apo_path = download_cif(apo_pdb_id, apo_dir)
            print(f"[APO] {apo_pdb_id}: downloaded -> {apo_path}")
        except Exception as e:
            print(f"[APO] {apo_pdb_id}: download failed: {e}")

        # HOLO
        holo_id = rec.get("holo_pdb_id")
        if not holo_id:
            print(f"[HOLO] {apo_pdb_id}: missing holo_pdb_id in dataset")
            continue
        try:
            holo_path = download_cif(holo_id, holo_dir)
            print(f"[HOLO] {holo_id}: downloaded -> {holo_path}")
        except Exception as e:
            print(f"[HOLO] {holo_id}: download failed: {e}")
            continue

        if args.write_holo_noligand:
            ligand_chain = rec.get("ligand_chain")
            ligand_index = rec.get("ligand_index")
            if not ligand_chain or not ligand_index:
                print(f"[HOLO-noLIG] {holo_id}: missing ligand_chain/index; skipping")
            else:
                try:
                    out_pdb = nolig_dir / f"{holo_id}_no_ligand.pdb"
                    write_holo_without_ligand(holo_path, out_pdb, ligand_chain, ligand_index)
                    print(f"[HOLO-noLIG] {holo_id}: wrote {out_pdb}")
                except Exception as e:
                    print(f"[HOLO-noLIG] {holo_id}: failed to create no-ligand PDB: {e}")

if __name__ == "__main__":
    main()
