#!/usr/bin/env python3
"""Read an authenticated, sealed Europe PMC metadata receipt locally.

Use only the intended private custodian key and an authorized directory.
No network calls, no output of private metadata to logs, no key export.
"""
import argparse, base64, hashlib, json, zipfile
from pathlib import Path
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import padding
from cryptography.hazmat.primitives.ciphers.aead import AESGCM
AAD=b"island-epmc-rights-statement-audit-v1"
def run(keyfile,archive,expected_zip_sha,aggregate,out):
    actual=hashlib.sha256(archive.read_bytes()).hexdigest()
    if actual!=expected_zip_sha.removeprefix("sha256:"):raise ValueError("Transport archive digest mismatch")
    with zipfile.ZipFile(archive) as z:
        if z.namelist()!=["metadata.sealed.json"] or z.testzip() is not None:raise ValueError("Unexpected transport ZIP")
        obj=json.loads(z.read("metadata.sealed.json"))
    if obj["schema"]!=1 or obj["aad"]!=AAD.decode() or obj["algorithm"]!="RSA-OAEP-SHA256 + AES-256-GCM":raise ValueError("Unexpected sealed schema")
    key=serialization.load_pem_private_key(keyfile.read_bytes(),password=None)
    pub=key.public_key().public_bytes(serialization.Encoding.PEM,serialization.PublicFormat.SubjectPublicKeyInfo)
    if hashlib.sha256(pub).hexdigest()!=obj["public_key_sha256"]:raise ValueError("Wrong recipient key")
    sym=key.decrypt(base64.b64decode(obj["wrapped_key"],validate=True),padding.OAEP(mgf=padding.MGF1(hashes.SHA256()),algorithm=hashes.SHA256(),label=AAD))
    plaintext=AESGCM(sym).decrypt(base64.b64decode(obj["nonce"],validate=True),base64.b64decode(obj["ciphertext"],validate=True),AAD)
    obj=json.loads(plaintext)
    if obj["summary"]!=json.loads(aggregate.read_text()):raise ValueError("Public aggregate differs from sealed summary")
    records=obj["records"];seen={(r["historical_source"],r["historical_id"]) for r in records}
    if len(records)!=len(seen) or len(records)!=obj["summary"]["target_article_ids"]:raise ValueError("Article roster mismatch")
    if not all(set(r)=={"historical_source","historical_id","search_metadata","statement"} for r in records):raise ValueError("Unexpected fields")
    import os
    os.umask(0o077)
    with out.open("xb") as f:f.write(plaintext+b"\n")
    out.chmod(0o600)
    print(json.dumps({"private_records":len(records),"receipt_sha256":hashlib.sha256(out.read_bytes()).hexdigest(),"verified_public_aggregate":True}))
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--key",type=Path,required=True)
    p.add_argument("--zip",type=Path,required=True)
    p.add_argument("--sha",required=True)
    p.add_argument("--aggregate",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.key,a.zip,a.sha,a.aggregate,a.out)
