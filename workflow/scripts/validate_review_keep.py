#!/usr/bin/env python3
import argparse, pathlib
p=argparse.ArgumentParser(); p.add_argument("--base",required=True); p.add_argument("--reviewed",required=True); p.add_argument("--output",required=True); a=p.parse_args()
def ids(path): return [tuple(x.split()[:2]) for x in open(path) if x.strip() and not x.startswith("#")]
base=set(ids(a.base)); reviewed=ids(a.reviewed)
if not reviewed: raise SystemExit(f"Review keep file is empty: {a.reviewed}")
bad=[x for x in reviewed if x not in base]
if bad: raise SystemExit(f"Review keep contains {len(bad)} IDs absent from previous stage; first: {bad[0]}")
pathlib.Path(a.output).parent.mkdir(parents=True,exist_ok=True)
with open(a.output,"w") as out:
    for fid,iid in reviewed: out.write(f"{fid}\t{iid}\n")

