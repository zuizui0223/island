"""Restore the archived, publicly tracked own/mixed files of one FD snapshot.

Only materializes the local archive. External candidates and excluded files
remain explicit omissions. Does not run historical workflows or fetch data.
"""
import argparse,hashlib,json
from pathlib import Path,PurePosixPath
def main():
    p=argparse.ArgumentParser()
    p.add_argument("snapshot")
    p.add_argument("--archive",type=Path,default=Path("research/fd_own_db"))
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if a.out.exists():
        raise ValueError("choose a new empty path")
    if a.snapshot!="main" and not (a.snapshot.startswith("pr") and a.snapshot[2:].isdigit()):
        raise ValueError("invalid snapshot")
    d=json.loads((a.archive/"snapshots"/(a.snapshot+".json")).read_text())
    n=0
    for e in d["files"]:
        if not e["archived"]:continue
        rel=PurePosixPath(e["path"])
        if rel.is_absolute() or ".." in rel.parts:raise ValueError("unsafe path")
        raw=(a.archive/e["object"]).read_bytes()
        if hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()!=e["git_blob"]:
            raise ValueError("hash mismatch")
        dest=a.out/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(raw)
        if e["mode"]=="100755":dest.chmod(0o755)
        n+=1
    (a.out/"MIGRATION_OMISSIONS.json").write_text(json.dumps([e for e in d["files"] if not e["archived"]],indent=2)+"\n")
    print("restored",n,"public files; no source input recovered or old CI activated")
if __name__=="__main__":main()
