"""Archive historical FD research without touching raw data or activating legacy CI.

Public fd Git objects only. Exact pinned main + 70 PR head snapshots.
External candidates are marked separately, not silently called independent.
All permitted public research blobs stored byte-identically and deduplicated;
snapshot indexes preserve path, blob, mode, provenance and disposition.
A legacy archival copy of external/mixed context is not an active analysis.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess

def git(cwd,*args):
    return subprocess.check_output(["git","-C",str(cwd),*args])

def put(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists() and path.read_bytes()!=data:
        raise ValueError("archive collision "+str(path))
    path.write_bytes(data)

def safe_path(path,mode,size):
    p=PurePosixPath(path)
    if mode not in ("100644","100755") or p.is_absolute() or ".." in p.parts:
        return False
    if size>8_000_000 or any(t in p.parts for t in ("data","private","raw","secrets","credentials",".git")):
        return False
    if len(p.parts)==1:
        return path in ("README.md",".gitignore","pyproject.toml","LICENSE","LICENSE.md")
    if p.parts[0] not in ("docs","results","src","scripts","tests","config",".github"):
        return False
    if p.parts[0]==".github" and p.parts[1]!="workflows":
        return False
    return p.suffix.lower() in {".py",".r",".md",".json",".csv",".tsv",".yml",".yaml",".toml",".txt",".mjs",".js",".svg",".png",".pdf",".rmd"}

def blob_hash(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--config",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    c=json.loads(a.config.read_text())
    if len(c["pins"])!=70 or set(r["pr"] for r in c["pins"])!=set(range(1,71)):
        raise ValueError("70 frozen PR snapshots required")
    refs=[{"pr":0,"head":c["source_public_main_sha"],"state":"main","title":"fd main at migration freeze"}]+c["pins"]
    if git(a.source,"rev-parse",c["source_public_main_sha"]+"^{tree}").decode().strip()!=c["source_main_tree"]:
        raise ValueError("main tree hash mismatch")
    ext=set(c["external_candidates_main"])
    for row in c["external_candidates_extra"]:
        ext.update(row["paths"])
    shared=set(c["shared_scaffold"])
    manifest=[]; objects={}; summary=[]; excluded=[]
    for r in refs:
        tag="main" if r["pr"]==0 else f"pr{r['pr']:03d}"
        sha=r["head"]
        assert git(a.source,"rev-parse",sha+"^{commit}").decode().strip()==sha
        entries=[]
        for item in git(a.source,"ls-tree","-rlz",sha).split(b"\0"):
            if not item: continue
            meta,rawpath=item.split(b"\t",1)
            mode,kind,oid,size=meta.decode().split()
            path=rawpath.decode("utf-8")
            if kind!="blob":
                excluded.append({"snapshot":tag,"path":path,"blob":oid,"reason":"nonregular"})
                continue
            size=int(size)
            ownership="external_candidate_fd" if path in ext else "shared_scaffold" if path in shared else "island_legacy_own_mixed_or_review"
            entry={"path":path,"mode":mode,"git_blob":oid,"bytes":size,"ownership":ownership}
            if not safe_path(path,mode,size):
                entry["archived"]=False
                excluded.append({"snapshot":tag,"path":path,"blob":oid,"reason":"not permitted public research type"})
                entries.append(entry)
                continue
            # The public external source-only candidates are preserved as exact
            # references in manifest, not copied as island-owned research.
            if ownership=="external_candidate_fd":
                entry["archived"]=False
                entry["reason"]="external-only candidate retained in original fd public history"
                entries.append(entry)
                continue
            if oid not in objects:
                raw=git(a.source,"cat-file","blob",oid)
                if len(raw)!=size or blob_hash(raw)!=oid:
                    raise ValueError("blob identity failure "+path)
                # Don't duplicate an already public raw dataset disguised as
                # LFS, an opaque archive, or credential-bearing file.
                if raw.startswith(b"version https://git-lfs.github.com/spec/"):
                    excluded.append({"snapshot":tag,"path":path,"blob":oid,"reason":"LFS pointer not copied"})
                    entry["archived"]=False
                    entries.append(entry)
                    continue
                if b"-----BEGIN PRIVATE KEY-----" in raw or re.search(rb"gh[pousr]_[A-Za-z0-9]{30,}",raw):
                    raise ValueError("credential-like content; stop without publishing "+path)
                suffix=PurePosixPath(path).suffix or ".txt"
                obj=f"objects/{oid[:2]}/{oid}{suffix}"
                put(a.out/obj,raw)
                objects[oid]=obj
            entry["archived"]=True
            entry["object"]=objects[oid]
            entries.append(entry)
        idx={"snapshot":tag,"source_commit":sha,"source_pr":r["pr"],"state":r["state"],"title":r["title"],"files":entries}
        put(a.out/"snapshots"/f"{tag}.json",(json.dumps(idx,indent=2,ensure_ascii=False)+"\n").encode())
        summary.append({"snapshot":tag,"commit":sha,"paths":len(entries),"archived":sum(x["archived"] for x in entries),"fd_external_candidates":sum(x["ownership"]=="external_candidate_fd" for x in entries)})
        manifest.extend((tag,sha,e["path"],e["git_blob"],e["ownership"],e.get("object",""),e["archived"]) for e in entries)
    # Verify all materialized bytes independently before committing.
    for oid,rel in objects.items():
        if blob_hash((a.out/rel).read_bytes())!=oid:
            raise ValueError("archive verification failure")
    with (a.out/"path_index.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["snapshot","source_commit","path","git_blob","ownership","local_object","archived"]);w.writerows(manifest)
    receipt={"source":c["source_repo"],"source_main_sha":c["source_public_main_sha"],"snapshots":len(summary),"path_occurrences":len(manifest),"unique_archived_blobs":len(objects),"archived_bytes":sum((a.out/v).stat().st_size for v in objects.values()),"excluded":excluded,"snapshot_support":summary,
      "limitations":"Conservative archival routing, not per-file scientific provenance certification. Verified external candidates excluded by path; unreviewed historical variants remain source-pinned. Old workflows are archival, not activated. No private source artifacts."}
    put(a.out/"receipt.json",(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n").encode())
    # Human-readable directory with relative links to all exact indexed worktrees.
    lines=["# FD legacy own/mixed research archive","",
    "All preserved bytes are addressable by their original Git blob SHA. Inspect the snapshot indexes or restore a selected worktree with the restore script. No legacy CI is active here. Historical narratives may have stale GIFT labels. Active ownership is defined by the repository source contract, not the archived wording.","",
    "| Snapshot | Pinned commit | Indexed paths | Local archived paths | External candidates retained in fd |","|---|---|---:|---:|---:|"]
    for x in summary:
        lines.append(f"| [{x['snapshot']}](snapshots/{x['snapshot']}.json) | [{x['commit'][:12]}](https://github.com/zuizui0223/fd/tree/{x['commit']}) | {x['paths']} | {x['archived']} | {x['fd_external_candidates']} |")
    put(a.out/"README.md",("\n".join(lines)+"\n").encode())
    print("FD_ARCHIVE_RECEIPT "+json.dumps({k:receipt[k] for k in ("snapshots","path_occurrences","unique_archived_blobs","archived_bytes")}|{"excluded":len(excluded)}))
    for x in excluded[:50]:
        print("FD_ARCHIVE_EXCLUDED "+json.dumps(x))

if __name__=="__main__":
    main()
