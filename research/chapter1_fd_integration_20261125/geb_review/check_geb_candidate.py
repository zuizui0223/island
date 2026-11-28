#!/usr/bin/env python3
from pathlib import Path
import re,json,hashlib
import xml.etree.ElementTree as ET
R=Path(__file__).resolve().parent
ROOT=R.parents[2]
M=(R/"MANUSCRIPT_BLINDED_REVIEW.md").read_text()
S=(R/"SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW.md").read_text()
C=(R/"COVER_LETTER_REVIEW.md").read_text()
P=(ROOT/"submission/geb_current/MANUSCRIPT_BLINDED.md").read_text()
a=M.split("## Abstract",1)[1].split("**Keywords:**",1)[0]
b="## Introduction"+M.split("## Introduction",1)[1].split("## References",1)[0]
aw=len(a.split());bw=len(b.split())
assert aw<=300 and bw<=5400
for h in ("Aim","Location","Time period","Major taxa studied","Methods","Results","Main conclusions"):
    assert a.count("**"+h+":**")==1,h
k=M.split("**Keywords:**",1)[1].splitlines()[0].strip().split(",")
k=[x.strip() for x in k];assert 6<=len(k)<=10 and k==sorted(k,key=str.casefold)
old=P.split("## Results\n",1)[1].split("## Discussion\n",1)[0]
new=M.split("## Results\n",1)[1].split("## Discussion\n",1)[0]
assert all(x in new for x in old.strip().split("\n\n") if x.strip()),"original result changed"
for x in (M,S,(R/"FIGURE_CAPTIONS_REVIEW.md").read_text(),(R/"Figure_S8.svg").read_text()):
    assert not re.search(r"(?i)zuizui0223|zhang ruiqi|@gmail\.com|@tohoku|10\.5281/zenodo\.22704973",x),"identity leak"
for x in ("0.02762","0.01967","0.02463","0.02089","0.000221","32 global","eight tropical"):
    assert x in M,x
assert S.count("# Appendix S8.")==1
assert all(("Table_S8"+x) in S for x in "abcde")
assert "PENDING" in M and "176,414" in M
cw=len(C.split("**Why this paper matters to GEB readers (under 250 words):**",1)[1].split("The analyses reuse",1)[0].split())
assert cw<=250
svg=ET.parse(R/"Figure_S8.svg").getroot()
st=" ".join(svg.itertext());assert "(a)" in st and "(b)" in st
print(json.dumps({"status":"GEB_INTEGRATED_TEXT_REVIEW_PASS","abstract_words":aw,"body_words":bw,"cover_interest_words":cw,"original_result_paragraphs_preserved":True,"new_H4_fit":False,"anonymous_text_checks":True,"full_visual_proof":False,"ready_to_submit":False}))
