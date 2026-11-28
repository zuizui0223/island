#!/usr/bin/env python3
"""Verify article-license statements for the sole-Europe-PMC-blocker cohort.

The official OAI pmc_fm record is front matter, not a full article. Only
already-OA, identity-matched search candidates are requested. The permitted
REST fullTextXML route is a fallback if the front-matter permissions are
missing. XML is processed transiently in memory; no full texts, abstracts,
species identities, or trait values are exported. The public result is an
aggregate. Optional sealed transport is encrypted for the private custodian.
Explicit article license metadata still does not establish the historical
acquisition version or rights for third-party data within the article.
"""
from __future__ import annotations
import argparse, base64, collections as C, hashlib, json, os, re, time
from pathlib import Path
from datetime import datetime, timezone
import xml.etree.ElementTree as ET
import requests
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import padding
from cryptography.hazmat.primitives.ciphers.aead import AESGCM
from europe_pmc_license_audit import load_inputs, get_metadata, normalize_license

OAI="https://europepmc.org/oai.cgi"
REST="https://www.ebi.ac.uk/europepmc/webservices/rest"
URI=re.compile(r'https?://(?:www\.)?creativecommons\.org/(?:licenses/(?:by(?:-sa|-nc|-nd|-nc-sa|-nc-nd)?)/[1-4]\.0|publicdomain/zero/1\.0)(?:/[^\s<>"\']*)?',re.I)
ALLOW={"cc-by","cc-by-sa","cc0"}
AAD=b"island-epmc-rights-statement-audit-v1"
def name(e):return e.tag.rsplit("}",1)[-1]
def text(e):return " ".join(" ".join(e.itertext()).split())
def canonical(v):
    k,n=normalize_license(v)
    if not n:return None
    if k=="cc0":return "https://creativecommons.org/publicdomain/zero/1.0/"
    if k.startswith("cc-by"):return "https://creativecommons.org/licenses/"+k[3:]+"/"+n+"/"
    return None
def ids_equal(ids,rec,key):
    pmc=str(rec.get("pmcid","")).upper()
    doi=str(rec.get("doi","")).lower().strip()
    got={str(x).upper().removeprefix("PMC") for x in ids.get("pmc",[])+ids.get("pmcid",[])}
    if not re.fullmatch(r"PMC[0-9]+",pmc) or pmc[3:] not in got:return False
    if not doi or doi not in [x.lower().strip() for x in ids.get("doi",[])]:return False
    if key[0]=="MED" and key[1] not in ids.get("pmid",[]):return False
    if key[0]=="PMC" and key[1]!=pmc:return False
    return True
def parse_xml(payload,rec,key):
    # Never resolve external XML entities. ElementTree is used only on the
    # official OA metadata response with a hard response-size ceiling.
    root=ET.fromstring(payload)
    if any(name(e)=="error" for e in root.iter()):return {"status":"provider_error"}
    metas=[e for e in root.iter() if name(e)=="article-meta"]
    if len(metas)!=1:return {"status":"article_meta_missing_or_ambiguous"}
    m=metas[0];ids=C.defaultdict(list)
    for e in m:
        if name(e)=="article-id":ids[e.attrib.get("pub-id-type","").lower()].append(text(e))
    if not ids_equal(ids,rec,key):return {"status":"identity_mismatch"}
    perms=[e for e in m if name(e)=="permissions"]
    if len(perms)!=1:return {"status":"article_permissions_missing_or_ambiguous"}
    p=perms[0]
    licenses=[e for e in p if name(e)=="license"]
    if not licenses:return {"status":"license_statement_missing"}
    raw=[];hashes0=[];uri_set=set();unrecognized=False
    for e in licenses:
        # Include explicit URI in license link attrs and license-p text, never
        # search references or body for an unrelated CC hyperlink.
        st=text(e);hashes0.append(hashlib.sha256(st.encode()).hexdigest())
        vals=[st]+[str(v) for node in e.iter() for k,v in node.attrib.items() if k.rsplit("}",1)[-1]=="href"]
        urls=[u.rstrip(".,;:)]}") for val in vals for u in URI.findall(val)]
        norm={canonical(u) for u in urls}
        norm.discard(None)
        if not norm:unrecognized=True
        uri_set.update(norm)
        raw.append(st)
    out={"license_statement_sha256":sorted(set(hashes0)),"license_uris":sorted(uri_set),
         "identity_verified":True,
         "article_license_type_labels":sorted(set(e.attrib.get("license-type","") for e in licenses)),
         "possible_third_party_exception_mentioned":bool(re.search(r"third.party|except where|unless otherwise|not included|credit line|separate permission"," ".join(raw),re.I))}
    if unrecognized or len(uri_set)!=1:
        return {"status":"version_uri_missing_or_multiple",**out}
    uri=next(iter(uri_set));lic,version=normalize_license(uri);meta_lic,_=normalize_license(rec.get("license"))
    if lic!=meta_lic or lic not in ALLOW:return {"status":"license_statement_conflict",**out}
    # Article-level statement only: retain actual version, do not infer 4.0.
    out["normalized_license"]=lic;out["version"]=version
    # Capture minimal citation metadata for private follow-up, not article text.
    title_groups=[e for e in m if name(e)=="title-group"]
    title=next((text(e) for t in title_groups for e in t if name(e)=="article-title"),"")
    out["title"]=title[:500]
    out["article_ids"]={k:v for k,v in ids.items() if k in {"pmid","pmc","pmcid","doi"}}
    out["copyright_years"]=[text(e) for e in p if name(e)=="copyright-year"]
    out["copyright_holders"]=[text(e)[:300] for e in p if name(e)=="copyright-holder"]
    return {"status":"explicit_permissive_article_license",**out}
def request_xml(session,url,params=None):
    for attempt in range(2):
        try:
            with session.get(url,params=params,timeout=(10,25),stream=True) as r:
                if r.status_code in (400,403,404,410):return None
                r.raise_for_status()
                b=bytearray()
                for chunk in r.iter_content(65536):
                    b.extend(chunk)
                    if len(b)>20*1024*1024:raise ValueError("response too large")
                return bytes(b)
        except (requests.RequestException,ValueError):
            if attempt==0:time.sleep(1)
    return None
def verify_statement(session,key,rec):
    pmc=str(rec.get("pmcid","")).upper()
    if not re.fullmatch(r"PMC[0-9]+",pmc) or rec.get("isOpenAccess")!="Y":
        return {"status":"oa_pmc_scope_missing"}
    params={"verb":"GetRecord","metadataPrefix":"pmc_fm","identifier":"oai:europepmc.org:"+pmc[3:]}
    x=request_xml(session,OAI,params)
    oai_status="request_failed"
    if x:
        try:
            y=parse_xml(x,rec,key);oai_status=y["status"]
            if y["status"]=="explicit_permissive_article_license":
                return {"retrieval":"OAI_pmc_fm","response_sha256":hashlib.sha256(x).hexdigest(),**y}
            # An explicit identity conflict or conflicting license is never
            # overridden by a second representation.
            if y["status"] in {"identity_mismatch","license_statement_conflict"}:
                return {"retrieval":"OAI_pmc_fm","response_sha256":hashlib.sha256(x).hexdigest(),**y}
        except ET.ParseError:oai_status="xml_parse_failure"
    # Only an OA article explicitly identified by the first-pass metadata may
    # use this official fallback endpoint, solely to inspect permissions.
    x=request_xml(session,REST+"/"+pmc+"/fullTextXML")
    if x:
        try:return {"retrieval":"REST_fullTextXML_permissions_only","oai_status":oai_status,"response_sha256":hashlib.sha256(x).hexdigest(),**parse_xml(x,rec,key)}
        except ET.ParseError:return {"status":"xml_parse_failure","oai_status":oai_status}
    return {"status":"statement_unavailable","oai_status":oai_status}
def summary(cells):
    return {"cells":len(cells),"by_axis":dict(C.Counter(c[1] for c in cells)),"by_quality":dict(C.Counter(c[2] for c in cells)),"derived_cells":sum(c[4] for c in cells)}
def seal(obj,public_key,out):
    pk=serialization.load_pem_public_key(public_key.read_bytes())
    key=AESGCM.generate_key(bit_length=256);nonce=os.urandom(12)
    payload=json.dumps(obj,sort_keys=True,allow_nan=False).encode()
    wrapped=pk.encrypt(key,padding.OAEP(mgf=padding.MGF1(hashes.SHA256()),algorithm=hashes.SHA256(),label=AAD))
    ciphertext=AESGCM(key).encrypt(nonce,payload,AAD)
    out.write_text(json.dumps({"schema":1,"algorithm":"RSA-OAEP-SHA256 + AES-256-GCM","aad":AAD.decode(),"public_key_sha256":hashlib.sha256(public_key.read_bytes()).hexdigest(),"wrapped_key":base64.b64encode(wrapped).decode(),"nonce":base64.b64encode(nonce).decode(),"ciphertext":base64.b64encode(ciphertext).decode()},sort_keys=True)+"\n")
def run(source,out,pub,sealed):
    keys,cells=load_inputs(source)
    eligible=[c for c in cells if c[3] and all(k[0] in {"MED","PMC"} for k in c[0])]
    roster=set().union(*(c[0] for c in eligible))
    assert len(eligible)==451 and len(roster)==368
    records,req=get_metadata(roster)
    candidates={k for k,v in records.items() if v["status"]=="permissive_metadata_candidate"}
    targeted=[c for c in eligible if c[0]<=candidates]
    target_keys=set().union(*(c[0] for c in targeted)) if targeted else set()
    print("EPMC_STATEMENT_TARGET "+json.dumps({"eligible_cells":len(eligible),"metadata_ids":len(roster),"candidate_cells":len(targeted),"target_article_ids":len(target_keys),"metadata_requests":req}),flush=True)
    session=requests.Session()
    session.headers["User-Agent"]="island-chapter1-license-audit/1.1 (permissions audit only)"
    checked={}
    for i,k in enumerate(sorted(target_keys),1):
        checked[k]=verify_statement(session,k,records[k]["record"])
        if i%20==0:print("EPMC_STATEMENT_PROGRESS "+str(i)+"/"+str(len(target_keys)),flush=True)
        time.sleep(.35)
    verified={k for k,v in checked.items() if v["status"]=="explicit_permissive_article_license"}
    by_uri=C.Counter(v["license_uris"][0] for k,v in checked.items() if k in verified)
    keep=[c for c in targeted if c[0]<=verified]
    exception_ids={k for k in verified if checked[k]["possible_third_party_exception_mentioned"]}
    keep_with_notice=[c for c in keep if c[0]&exception_ids]
    result={"timestamp_utc":datetime.now(timezone.utc).isoformat(),"eligible_sole_family_med_pmc":summary(eligible),"metadata_ids":len(roster),"metadata_statuses":dict(C.Counter(v["status"] for v in records.values())),"metadata_requests":req,"metadata_candidate_cells":summary(targeted),"previous_full_audit_candidate_cells":87,"target_article_ids":len(target_keys),"statement_statuses":dict(C.Counter(v["status"] for v in checked.values())),"verified_article_license_uri_counts":dict(by_uri),"article_statement_retrieval":dict(C.Counter(v.get("retrieval","unavailable") for v in checked.values())),"all_required_article_statements_verified_cells":summary(keep),"verified_article_cells_with_third_party_notice":summary(keep_with_notice),"new_cells_authorized":0,"policy_changed":False,"scope":"Explicit article license URI/identity verification, not original evidence version/third-party-material clearance or attribution approval. No rights promotion."}
    out.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
    if pub and sealed:
        private={"summary":result,"records":[{"historical_source":k[0],"historical_id":k[1],"search_metadata":records[k],"statement":checked[k]} for k in sorted(target_keys)]}
        seal(private,pub,sealed)
    print("EPMC_STATEMENT_RESULT "+json.dumps(result,sort_keys=True,allow_nan=False),flush=True)
    if not req["successful_batches"]:raise SystemExit("Metadata request failure")
def self_test():
    rec={"pmcid":"PMC123","doi":"10.1/x","license":"cc by","isOpenAccess":"Y"}
    x=b'<article><front><article-meta><article-id pub-id-type="pmc">123</article-id><article-id pub-id-type="pmid">456</article-id><article-id pub-id-type="doi">10.1/x</article-id><permissions><license xmlns:xlink="http://www.w3.org/1999/xlink" xlink:href="https://creativecommons.org/licenses/by/4.0/"><license-p>Licensed under CC BY.</license-p></license></permissions></article-meta></front></article>'
    r=parse_xml(x,rec,("MED","456"))
    assert r["status"]=="explicit_permissive_article_license" and r["version"]=="4.0"
    assert parse_xml(x,rec,("MED","999"))["status"]=="identity_mismatch"
    assert parse_xml(x.replace(b"/by/4.0",b"/by-nc/4.0"),rec,("MED","456"))["status"]=="license_statement_conflict"
    assert parse_xml(x.replace(b"https://creativecommons.org/licenses/by/4.0/",b"https://example.org/"),rec,("MED","456"))["status"]=="version_uri_missing_or_multiple"
    print("EPMC_STATEMENT_SELF_TEST_PASS")
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--source-dir",type=Path);p.add_argument("--out",type=Path);p.add_argument("--public-key",type=Path);p.add_argument("--sealed-out",type=Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test()
    else:run(a.source_dir,a.out,a.public_key,a.sealed_out)
