#!/usr/bin/env python3
"""Read-only archive verification, or exact-byte backup into explicitly authorized storage.

Does not upload, publish, change retention settings, or extract raw data. Output
receipts contain only public artifact metadata, byte hashes and pin checks.
"""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, tempfile, zipfile
from pathlib import Path, PurePosixPath

HERE=Path(__file__).resolve().parent
MANIFEST=HERE/"SOURCE_PRESERVATION_MANIFEST.json"
REPO=HERE.parents[2]
def stream_hash(f):
    h=hashlib.sha256()
    for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def file_hash(path):
    with path.open("rb") as f:return stream_hash(f)
def ensure_archive(path,spec):
    if path.stat().st_size!=spec["bytes"] or file_hash(path)!=spec["digest"].split(":")[1]:
        raise ValueError("archive size or digest mismatch: "+str(spec["id"]))
def inspect(path,spec,pins):
    ensure_archive(path,spec)
    found=[]
    with zipfile.ZipFile(path) as z:
        names=z.namelist()
        if len(names)!=len(set(names)):raise ValueError("duplicate ZIP entry")
        if any(PurePosixPath(n).is_absolute() or ".." in PurePosixPath(n).parts or "\\" in n for n in names):
            raise ValueError("unsafe ZIP member")
        for p in pins:
            if p["archive_id"] not in (None,spec["id"]):continue
            cand=[n for n in names if PurePosixPath(n).name==p["basename"]]
            hashes=[]
            for n in cand:
                with z.open(n) as f:hashes.append(stream_hash(f))
            matches=hashes.count(p["sha256"])
            found.append({"label":p["label"],"archive_id":spec["id"],"basename_matches":len(cand),"exact_hash_matches":matches})
            if p["archive_id"]==spec["id"] and matches<1:
                raise ValueError("required inner pin absent or different: "+p["label"])
        if z.testzip() is not None:raise ValueError("ZIP CRC verification failed")
    return {"id":spec["id"],"run":spec["run"],"bytes":path.stat().st_size,"archive_sha256":spec["digest"],"crc_verified":True,"pins":found}
def download(repo,spec,dest):
    name=str(spec["id"])+"-"+spec["name"]+".zip"
    p=dest/name
    if p.exists():
        if p.is_symlink():raise ValueError("refusing symlink destination")
        ensure_archive(p,spec);return p
    fd,tmp=tempfile.mkstemp(prefix=".source-",suffix=".partial",dir=dest);os.close(fd);tmp=Path(tmp)
    try:
        with tmp.open("wb") as f:
            subprocess.run(["gh","api",f"repos/{repo}/actions/artifacts/{spec['id']}/zip"],stdout=f,check=True)
        ensure_archive(tmp,spec)
        # No silent overwrite if another custodian populated the destination.
        try:os.link(tmp,p)
        except FileExistsError:ensure_archive(p,spec)
        return p
    finally:tmp.unlink(missing_ok=True)
def run(manifest,dest,all_snapshots):
    rows=[]
    for a in manifest["archive_metadata"]:
        if not all_snapshots and not a["required_for_current_lineage"]:continue
        path=download(manifest["owner"],a,dest)
        rows.append(inspect(path,a,manifest["inner_pins"]))
    pins=[]
    for p in manifest["inner_pins"]:
        n=sum(x["exact_hash_matches"] for row in rows for x in row["pins"] if x["label"]==p["label"])
        pins.append({"label":p["label"],"exact_hash_matches":n,"status":"verified" if n else "not_recovered"})
    repo=[]
    for p in manifest["repo_pins"]:
        path=REPO/p["path"]
        v=file_hash(path)
        item={"path":p["path"],"sha256":v,"expected_sha256":p["sha256"],"status":"exact_bytes" if v==p["sha256"] else "unresolved_mismatch"}
        # Only the historically Windows-generated YAML may have a checkout
        # line-ending discrepancy. Test exact CRLF reconstruction, not arbitrary
        # YAML reserialization or semantic equivalence. Keep both hashes.
        if v!=p["sha256"] and p.get("allow_exact_crlf_reconstruction"):
            raw=path.read_bytes()
            crlf=raw.replace(b"\r\n",b"\n").replace(b"\n",b"\r\n")
            reconstructed=hashlib.sha256(crlf).hexdigest()
            git_blob=subprocess.check_output(["git","hash-object",str(path)],text=True).strip()
            item.update({"git_blob":git_blob,"crlf_sha256":reconstructed})
            if reconstructed==p["sha256"] and git_blob==p["git_blob"]:
                item["status"]="exact_historical_bytes_reconstructed_by_crlf"
                item["transformation"]="LF to CRLF only; original checkout unchanged"
        repo.append(item)
    unresolved=[x for x in repo if x["status"]=="unresolved_mismatch"]
    return {"owner":manifest["owner"],"archives":rows,"inner_pins":pins,"repo_pins":repo,
        "no_upload":True,"unavailable_originals_not_restored":manifest["unavailable"],"unresolved_repo_pins":unresolved,"all_required_repository_pins_verified":not unresolved}
def self_test():
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/"x.zip";b=b"synthetic only"
        with zipfile.ZipFile(p,"w") as z:z.writestr("a/test.dat",b)
        spec={"id":1,"run":2,"bytes":p.stat().st_size,"digest":"sha256:"+file_hash(p)}
        pin={"label":"synthetic","archive_id":1,"basename":"test.dat","sha256":hashlib.sha256(b).hexdigest()}
        assert inspect(p,spec,[pin])["pins"][0]["exact_hash_matches"]==1
        try:inspect(p,spec,[dict(pin,sha256="0"*64)])
        except ValueError:pass
        else:raise AssertionError("bad inner hash accepted")
    print("PRESERVATION_SELF_TEST_PASS")
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--verify-only",action="store_true",help="temporary validation; does not preserve raw bytes")
    a.add_argument("--destination",type=Path,help="existing authorized local/institutional storage")
    a.add_argument("--ack-authorized-storage",action="store_true")
    a.add_argument("--include-alpha1",action="store_true")
    a.add_argument("--receipt",type=Path)
    a.add_argument("--self-test",action="store_true")
    p=a.parse_args()
    if p.self_test:self_test();return
    m=json.loads(MANIFEST.read_text())
    os.umask(0o077)
    if p.verify_only:
        if p.destination is not None:raise ValueError("choose one mode")
        with tempfile.TemporaryDirectory(prefix="island-preserve-") as d:result=run(m,Path(d),p.include_alpha1)
        result["durable_custodian_backup_created"]=False
        result["mode"]="temporary_verify_only"
    else:
        if p.destination is None or not p.ack_authorized_storage:raise ValueError("explicit authorized existing destination required")
        dest=p.destination.expanduser().resolve(strict=True)
        if not dest.is_dir():raise ValueError("destination must be existing directory")
        result=run(m,dest,p.include_alpha1)
        result["durable_custodian_backup_created"]=True
        result["mode"]="explicit_custodian_local_backup"
    if p.receipt:
        p.receipt.parent.mkdir(parents=True,exist_ok=True)
        p.receipt.write_text(json.dumps(result,indent=2)+"\n")
    print("PRESERVATION_VERIFIED "+json.dumps(result,allow_nan=False),flush=True)
    if not result["all_required_repository_pins_verified"]:raise SystemExit("Unresolved repository pin mismatch; receipt retained")
if __name__=="__main__":main()
