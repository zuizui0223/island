from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from lxml import etree as E
import json
w=Path.cwd();repo=w/'work/island-traitwise-20261004';out=repo/'results/poster_q1_20261004'
p=w/'outputs/Island_Biology_2026_A0_parallel_methods_20261001_v2.pptx'
with ZipFile(p) as z:d={k:z.read(k) for k in z.namelist()}
n={'a':'http://schemas.openxmlformats.org/drawingml/2006/main','p':'http://schemas.openxmlformats.org/presentationml/2006/main'}
r=E.fromstring(d['ppt/slides/slide1.xml'])
for i,name in zip([6,7,8,9],['poster_h1','poster_h2','poster_h3','poster_h4']):
 d[f'ppt/media/q1_preserved_image{i}.png']=(out/f'{name}.png').read_bytes()
 d[f'ppt/media/q1_preserved_image{i-3}.svg']=(out/f'{name}.svg').read_bytes()
replacements={
 'Joint Wald test of distance slopes.':'Spatial-block t tests (G−1 df).',
 'Spatial-block robust CIs. Joint q<.005 in all regions. Dots = islands.':'Spatial-block 95% t CIs. Dots = islands. No pooled syndrome test.',
 'Supplemental-only: selfing p=.290; accessibility p=.0397. Post-hoc.':'Supplemental-only: selfing p=.291; accessibility p=.0424. Post-hoc.',
 '(p=.00037).':'(p=.00112).',
 'Tropical Direct access: q=.120.':'Tropical Direct access: q=.127.',
 'H4 atomic traits: autonomous selfing p=2.9×10⁻⁸; radial symmetry p=1.2×10⁻⁵; generalized form p=.0443; self-compatibility p=.150.':'H4 atomic traits: autonomous selfing p=4.85×10⁻⁸; radial symmetry p=1.47×10⁻⁵; generalized form p=.0454; self-compatibility p=.150.',
 'Publication-clustered SE.':'Publication-clustered SE; t inference.',
}
changes=[]
for node in r.xpath('//a:t',namespaces=n):
 old=node.text or '';new=old
 for a,b in replacements.items():new=new.replace(a,b)
 if new!=old:changes.append({'before':old,'after':new});node.text=new
# Add concise WCVP receipt to notes; preserve layout and Q2 exactly.
notes='ppt/notesSlides/notesSlide1.xml'
if notes in d:
 nr=E.fromstring(d[notes]);ts=nr.xpath('//a:t',namespaces=n)
 if ts:ts[0].text=(ts[0].text or '')+' Q1 UPDATE 2026-10-04: H1 uses final seven separate All traits, spatial-block t(G-1) intervals and unadjusted two-sided p, no Holm and no pooled score. WCVP sensitivity is supplied as a separate companion figure: 513320 records / 2372 islands before filtering, regional compatibility not exact island nativity. H2-H4 retain current finite-cluster/publication inference; H2 retains its existing BH q family. H4 figure shows model contrasts with pointwise t intervals, not raw scatter. Q2 retained unchanged; not reaudited by this Q1 update.'
 d[notes]=E.tostring(nr,xml_declaration=True,encoding='UTF-8',standalone=True)
d['ppt/slides/slide1.xml']=E.tostring(r,xml_declaration=True,encoding='UTF-8',standalone=True)
final=w/'outputs/Island_Biology_2026_A0_Q1_final_20261004.pptx'
with ZipFile(final,'w',ZIP_DEFLATED) as z:
 for k,v in d.items():z.writestr(k,v)
(out/'poster_text_changes.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2),encoding='utf-8')
print(final)
