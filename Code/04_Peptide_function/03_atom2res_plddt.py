#!/usr/bin/env python3
import argparse, os, sys, math, glob
from collections import defaultdict, namedtuple

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
       "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE"}

Rec = namedtuple("Rec","model chain resseq icode resname atom altloc b")

def parse_pdb(pdb_path):
    recs = []
    with open(pdb_path) as f:
        model = 1
        for line in f:
            if line.startswith("MODEL"):
                try: model = int(line.split()[1])
                except: model += 1
            if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
            resname = line[17:20].strip()
            if resname not in AA3: continue
            atom   = line[12:16].strip()
            altloc = line[16].strip()
            chain  = line[21].strip() or "_"
            try:
                resseq = int(line[22:26])
            except:
                continue
            icode = line[26].strip()
            try:
                b = float(line[60:66])  # B-factor = pLDDT
            except:
                continue
            recs.append(Rec(model, chain, resseq, icode, resname, atom, altloc, b))
    return recs

def parse_mmcif(cif_path):
    """
    Light-weight mmCIF parser for per-atom B-factors.
    Requires that the file has a standard _atom_site loop with these tags:
    label_asym_id, label_comp_id, label_seq_id, auth_asym_id, auth_comp_id,
    auth_seq_id, B_iso_or_equiv, label_alt_id, group_PDB, type_symbol, pdbx_PDB_model_num, ins_code
    We’ll use auth_* when available; fall back to label_*.
    """
    # Minimal mmCIF reader
    def read_loop(lines, start):
        tags = []
        i = start
        while i < len(lines) and lines[i].startswith("_"):
            tags.append(lines[i].strip())
            i += 1
        rows = []
        if i < len(lines) and lines[i].lower().startswith("loop_"):
            i += 1
        # Already consumed loop_ and tags; now rows until blank or next loop_/data_
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("data_") or s.startswith("loop_") or s.startswith("_"):
                break
            rows.append(s)
            i += 1
        return tags, rows, i

    with open(cif_path) as fh:
        lines = fh.readlines()

    recs = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.lower() == "loop_":
            # Look ahead to tags
            i += 1
            start = i
            while i < len(lines) and lines[i].startswith("_"):
                i += 1
            tags = [t.strip() for t in lines[start:i]]
            if "_atom_site.B_iso_or_equiv" in tags and "_atom_site.type_symbol" in tags:
                # Collect rows
                rows = []
                while i < len(lines):
                    s = lines[i].strip()
                    if not s or s.startswith("data_") or s.startswith("loop_") or s.startswith("_"):
                        break
                    rows.append(s)
                    i += 1
                # Tokenize with CIF rules (very simplified: honor quoted fields)
                def split_cif_row(row):
                    out, tok, in_sq, in_dq = [], "", False, False
                    for c in row:
                        if c=="'" and not in_dq: in_sq = not in_sq; continue
                        if c=='"' and not in_sq: in_dq = not in_dq; continue
                        if c.isspace() and not in_sq and not in_dq:
                            if tok: out.append(tok); tok=""
                        else:
                            tok += c
                    if tok: out.append(tok)
                    return out
                tag_idx = {t:i for i,t in enumerate(tags)}
                def get(row, key, default=""):
                    j = tag_idx.get(key, None)
                    if j is None: return default
                    v = row[j] if j < len(row) else default
                    return v if v!="." and v!="?" else default

                for r in rows:
                    cols = split_cif_row(r)
                    group = get(cols, "_atom_site.group_PDB")
                    if group not in ("ATOM","HETATM"): continue
                    resn = get(cols, "_atom_site.auth_comp_id") or get(cols, "_atom_site.label_comp_id")
                    if resn not in AA3: continue
                    chain = get(cols, "_atom_site.auth_asym_id") or get(cols, "_atom_site.label_asym_id") or "_"
                    try:
                        resseq = int(get(cols, "_atom_site.auth_seq_id") or get(cols, "_atom_site.label_seq_id"))
                    except:
                        continue
                    icode = get(cols, "_atom_site.pdbx_PDB_ins_code") or ""
                    atom  = get(cols, "_atom_site.label_atom_id") or get(cols, "_atom_site.auth_atom_id") or ""
                    alt   = get(cols, "_atom_site.label_alt_id")
                    try:
                        b = float(get(cols, "_atom_site.B_iso_or_equiv"))
                    except:
                        continue
                    try:
                        model = int(get(cols, "_atom_site.pdbx_PDB_model_num") or "1")
                    except:
                        model = 1
                    recs.append(Rec(model, chain, resseq, icode, resn, atom, alt, b))
            else:
                # Not the atom loop; continue
                continue
        else:
            i += 1
    return recs

def per_residue_stats(recs, mode="mean"):
    """
    Aggregate per (model, chain, resseq, icode, resname).
    mode: "mean" (default), "median", or "CA"
    """
    # filter altLoc: prefer altloc == '' or 'A'
    filtered = [r for r in recs if (r.altloc in ("", "A"))]
    key = lambda r: (r.model, r.chain, r.resseq, r.icode, r.resname)

    by_res = defaultdict(list)
    for r in filtered:
        by_res[key(r)].append(r)

    import statistics as stats
    out = []
    for k, atoms in sorted(by_res.items()):
        if mode.upper()=="CA":
            ca = [a for a in atoms if a.atom == "CA"]
            if not ca: 
                # fallback to mean if no CA present
                val = stats.mean(a.b for a in atoms)
                n = len(atoms)
            else:
                val = ca[0].b
                n = 1
        elif mode.lower()=="median":
            val = stats.median(a.b for a in atoms)
            n = len(atoms)
        else:
            vals = [a.b for a in atoms]
            val = sum(vals)/len(vals)
            n = len(vals)
        out.append((*k, round(val,2), n))
    return out  # list of tuples

def write_tsv(records, out_path):
    with open(out_path, "w") as w:
        w.write("model\tchain\tresnum\ticode\tresname\tplddt\tatom_count\n")
        for model, chain, resseq, icode, resname, plddt, n in records:
            resid = f"{resseq}{icode}" if icode else str(resseq)
            w.write(f"{model}\t{chain}\t{resid}\t{icode}\t{resname}\t{plddt:.2f}\t{n}\n")

def main():
    ap = argparse.ArgumentParser(description="Convert per-atom pLDDT to per-residue pLDDT (PDB or mmCIF).")
    ap.add_argument("input", help="PDB/mmCIF file, directory, or glob (e.g. 'out/*.pdb')")
    ap.add_argument("-o","--outdir", default="plddt_residue")
    ap.add_argument("--mode", choices=["mean","median","CA"], default="mean",
                    help="Aggregation per residue (default: mean).")
    args = ap.parse_args()

    paths = []
    if os.path.isdir(args.input):
        for ext in ("*.pdb","*.cif","*.mmcif"):
            paths.extend(glob.glob(os.path.join(args.input, ext)))
    else:
        globbed = glob.glob(args.input)
        paths = globbed if globbed else [args.input]

    if not paths:
        sys.exit("No input files found.")

    os.makedirs(args.outdir, exist_ok=True)

    for p in sorted(paths):
        ext = os.path.splitext(p)[1].lower()
        if ext == ".pdb":
            recs = parse_pdb(p)
        elif ext in (".cif",".mmcif"):
            recs = parse_mmcif(p)
        else:
            print(f"Skip {p} (unknown extension)", file=sys.stderr); continue
        stats = per_residue_stats(recs, mode=args.mode)
        base = os.path.basename(p)
        out  = os.path.join(args.outdir, base.rsplit(".",1)[0] + f".residue_plddt.{args.mode}.tsv")
        write_tsv(stats, out)
        print(f"Wrote {out} ({len(stats)} residues)")

if __name__ == "__main__":
    main()
