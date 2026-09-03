#!/usr/bin/env python3
import argparse, pathlib
p=argparse.ArgumentParser(); p.add_argument("--output",required=True); p.add_argument("stages",nargs="+"); a=p.parse_args()
pathlib.Path(a.output).parent.mkdir(parents=True,exist_ok=True)
with open(a.output,"w") as out:
    out.write("stage\tn_samples\n")
    for item in a.stages:
        stage,path=item.split("=",1)
        n=sum(1 for x in open(path) if x.strip() and not x.startswith("#"))
        out.write(f"{stage}\t{n}\n")

