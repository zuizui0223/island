import fs from "node:fs";
import path from "node:path";
import {fileURLToPath} from "node:url";
const dir=path.resolve(path.dirname(fileURLToPath(import.meta.url)),"../docs/manuscript_20261124");
const data=JSON.parse(fs.readFileSync(path.join(dir,"supplement_data.json"),"utf8"));
const render=function render(D){
 const C={ink:"#172B3A",sub:"#546675",teal:"#007C78",blue:"#3C6FA3",orange:"#B96A34",gray:"#A0ADB5",pale:"#E7EDF1",line:"#DCE3E7",white:"#FFFFFF"};
 const e=s=>String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;");
 const t=(x,y,s,z=14,w=400,a="start",c=C.ink)=>`<text x="${x}" y="${y}" font-size="${z}" font-weight="${w}" text-anchor="${a}" fill="${c}">${e(s)}</text>`;
 const l=(x,y,x2,y2,c=C.line,w=1,d=false)=>`<line x1="${x}" y1="${y}" x2="${x2}" y2="${y2}" stroke="${c}" stroke-width="${w}"${d?' stroke-dasharray="5 4"':""}/>`;
 const r=(x,y,w,h,c,rx=2)=>`<rect x="${x}" y="${y}" width="${Math.max(w,0)}" height="${h}" rx="${rx}" fill="${c}"/>`;
 const p=(x,y,c,open=false)=>`<circle cx="${x}" cy="${y}" r="5" stroke="${c}" stroke-width="1.7" fill="${open?C.white:c}"/>`;
 const f=(lo,hi,x,w)=>v=>x+(v-lo)/(hi-lo)*w;
 const signed=(v,d=3)=>(v>0?"+":"")+v.toFixed(d);
 const ax=(x,y,w,h,lo,hi,ticks,label)=>{const z=f(lo,hi,x,w);let a="";for(const v of ticks){a+=l(z(v),y,z(v),y+h,v===0?C.sub:C.line,v===0?1.2:1,v===0)+t(z(v),y+h+21,v,12,400,"middle",C.sub)}return [a+l(x,y+h,x+w,y+h,C.sub)+t(x+w/2,y+h+44,label,13,400,"middle"),z]};
 const range=(f,y,v,b,c=C.teal)=>l(f(b[0]),y,f(b[1]),y,c,2)+l(f(b[0]),y-5,f(b[0]),y+5,c)+l(f(b[1]),y-5,f(b[1]),y+5,c)+p(f(v),y,c);
 const panel=(x,y,a,s)=>t(x,y,a,20,700)+t(x+30,y,s,16,700);
 const svg=(w,h,ttl,sub,body,foot)=>`<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" role="img" aria-label="${e(ttl)}"><rect width="100%" height="100%" fill="#FFFFFF"/><g font-family="Arial, Helvetica, sans-serif">${t(24,32,ttl,23,700)}${t(24,55,sub,13,400,"start",C.sub)}${l(24,69,w-24,69)}${body}${l(24,h-34,w-24,h-34)}${t(24,h-14,foot,11,400,"start",C.sub)}</g></svg>\n`;
 const out={};
 {let b=panel(24,102,"A","Separate source-relative contrasts")+panel(650,102,"B","Joint geological/LGM support");
 const [a,z]=ax(255,130,335,170,0,12,[0,3,6,9,12],"Island − source joint-frequency contrast (pp)");b+=a;
 D.geology.contrasts.forEach((v,i)=>{let y=174+82*i;b+=t(242,y-6,v[0],13,400,"end")+t(242,y+13,`n = ${v[4]} vs ${v[5]}`,12,400,"end",C.sub)+range(z,y,v[1],[v[2],v[3]])+t(z(v[1]),y-14,signed(v[1],2),13,700,"middle")});
 b+=t(865,139,"LGM connected",13,700,"middle")+t(1000,139,"LGM unconnected",13,700,"middle");
 D.geology.overlap.forEach((a,i)=>{b+=t(782,184+i*70,["Continental","Oceanic"][i],13,700,"end");a.forEach((v,j)=>{b+=r(803+135*j,151+i*70,124,61,v===0?"#F5E4DA":C.pale)+t(865+135*j,189+i*70,v,24,700,"middle",v===0?C.orange:C.ink)})});
 b+=t(650,322,"No jointly labelled oceanic-connected islands.",13,700,"start",C.orange)+t(650,343,"Geological origin and past connection cannot be isolated.",12,400,"start",C.sub);
 out["figS01_geology_connectivity.svg"]=svg(1100,400,"Figure S1 | Geology and past land connectivity","Tropical actinomorphic × deep-tube phenotype; separate comparison cohorts.",b,"Source: verified original geological/LGM evidence ledger. Whiskers: archived exploratory grouped intervals.");
 }
 {const es=[...new Set(D.q.map(v=>v.endpoint))].sort(),rs=["northern_midlatitude","tropical","southern_extratropical"],sc=[...new Set(D.q.map(v=>v.scope))];if(es.length!==13||sc.length!==2)throw Error("Q layout");
 let b="";const colors={positive:C.teal,negative:C.orange,ambiguous:C.gray,unsupported:"#E6EBEF"};const sym={positive:"+",negative:"−",ambiguous:"?",unsupported:"·"};
 for(let si=0;si<2;si++){let x=260+si*466;b+=panel(x,103,["A","B"][si],sc[si].replaceAll("_"," "));rs.forEach((k,j)=>{let xx=x+j*146;b+=t(xx+64,134,["North mid","Tropics","South extra"][j],13,700,"middle")+t(xx+32,153,"Category",11,400,"middle")+t(xx+100,153,"Q",11,400,"middle")})}
 es.forEach((k,i)=>{let y=169+i*31;b+=t(24,y+17,k.replaceAll("_"," "),12.5);for(let si=0;si<2;si++)for(let j=0;j<3;j++){let z=D.q.find(v=>v.scope===sc[si]&&v.region===rs[j]&&v.endpoint===k);if(!z)throw Error("missing Q cell");["category","Q"].forEach((key,a)=>{let x=260+si*466+j*146+a*68;let v=z[key];if(!colors[v])throw Error("unknown Q state "+v);b+=r(x,y,62,24,colors[v],2)+t(x+31,y+17,sym[v],14,700,"middle",v==="unsupported"?C.sub:C.white)})}});
 let y=596;["positive","negative","ambiguous","unsupported"].forEach((v,i)=>{let x=24+i*270;b+=r(x,y,20,18,colors[v])+t(x+28,y+14,v,13)});
 b+=t(24,642,"All-recorded: 9/39 category gains, 13/39 corrected-Q gains, 6/39 both; 3/39 unsupported.",14,700)+t(24,666,"Native-only: 1/39 category gains, 4/39 corrected-Q gains, 1/39 both; 16/39 unsupported.",14);
 out["figS02_all13_nominal_Q.svg"]=svg(1210,718,"Figure S2 | Categorical prediction and nominal functional diversity","All 78 archived trait × macroregion × provenance cells; gain relative to strong local category prior.",b,"Category = dual loss/dual grouping gate; Q = finite-sample corrected nominal Rao Q dual grouping gate. Not continuous trait volume.");
 }
 {let b=panel(24,102,"A","Wing shape and within-archipelago occupancy")+panel(24,365,"B","Butterfly recorded wing size and larval host breadth");
 const [a,z]=ax(300,127,425,170,0,1.2,[0,.3,.6,.9,1.2],"Occupancy-fraction slope per +10 HWI");b+=a;
 D.flight.forEach((v,i)=>{let y=146+i*33;b+=t(285,y+4,`${v[0]} (species=${v[1]})`,13,400,"end")+p(z(v[4]),y,i===3?C.orange:C.blue)+t(754,y+4,`${v[2]}/${v[3]} positive group correlations`,13)});
 b+=t(24,342,"≥50% diet estimate is sparse and exceeds the occupancy scale; not a credible comparative effect size.",12,700,"start",C.orange);
 const masks=["all","precision_le1km","no_gbifid_with_reference"],labels=["All","Precision ≤1 km","No GBIF ID + reference"];const tt=["log_wingspan","log_host_families"];
 const [a2,z2]=ax(360,404,340,252,-.06,.04,[-.06,-.04,-.02,0,.02,.04],"Log-trait slope / isolation doubling");b+=a2;
 const [a3,z3]=ax(806,404,340,252,-.06,.04,[-.06,-.04,-.02,0,.02,.04],"Log-trait slope / isolation doubling");b+=a3;
 b+=t(530,395,"Wingspan",14,700,"middle")+t(976,395,"Documented host-family count",14,700,"middle");
 masks.forEach((m,i)=>[10,50].forEach((n,j)=>{let y=424+(i*2+j)*40;b+=t(346,y+4,`${labels[i]} / ≥${n} records`,12.5,400,"end");tt.forEach((v,k)=>{const name=`${m}_${v}_${n}`,rs=D.butter.filter(q=>q.panel===name);const zz=k===0?z2:z3;if(!rs.length){b+=t(zz(0)+6,y+4,"unsupported",11,400,"start",C.sub);return}rs.forEach((q,a)=>{const col=q.model==="raw"?C.gray:C.teal;b+=p(zz(q.slope),y+(q.model==="raw"?-6:6),col,q.model==="raw")});b+=t(zz(-.06)+5,y+17,`n=${rs[0].n}; g=${rs[0].groups}`,10,400,"start",C.sub)})}));
 b+=p(25,724,C.gray,true)+t(40,729,"Area adjusted",12)+p(205,724,C.teal)+t(220,729,"Also record-count adjusted",12)+t(535,729,"Source masks and thresholds overlap; no independent ecological replicates.",12,400,"start",C.sub);
 b+=t(24,762,"Host-family count and host documentation are related; wing shape and wing size do not measure realized dispersal.",13,700);
 out["figS03_animal_trait_sensitivity.svg"]=svg(1190,810,"Figure S3 | Flight morphology versus recorded resource breadth","Bird HWI and butterfly traits: separate populations, different observation processes.",b,"Sources: HWI run 45341369223 and butterfly trait filtering run 45358140791. Complete butterfly table provided.");
 }
 {let b=panel(24,105,"A","All qualified original permutation comparisons")+panel(24,390,"B","Adjusted evidence gate");
 const x0=310,w=500;const names={source_relative:"Source-relative joint phenotype",form_depth:"Direct form and depth",dual:"Two-sided taxonomic null"};
 const z=f(0,500,x0,w);
 Object.entries(D.nullcounts).forEach(([k,v],i)=>{let y=153+i*74;b+=t(x0-15,y-3,names[k],14,700,"end")+t(x0-15,y+16,`${v.n} qualified tests`,12,400,"end",C.sub)+r(x0,y-14,z(v.n)-x0,26,C.pale)+r(x0,y-14,z(v.raw)-x0,26,C.teal)+t(830,y+3,`${v.raw} raw p < 0.05`,14,700)+t(830,y+22,`${v.adjusted} adjusted p < 0.05`,12,400,"start",C.sub)});
 b+=r(25,342,20,16,C.pale)+t(53,355,"All qualified",12)+r(200,342,20,16,C.teal)+t(228,355,"Nominal passes (not discoveries)",12);
 Object.entries(D.nullcounts).forEach(([k,v],i)=>{let x=24+i*345;b+=r(x,413,322,90,"#F1F5F7",5)+t(x+15,437,names[k],13,700)+t(x+15,467,`0 / ${v.n}`,25,700,"start",C.orange)+t(x+15,490,`Minimum adjusted p = ${v.minq.toFixed(3)}`,12,400,"start",C.sub)});
 b+=t(24,532,"996 tests across three correlated analysis families; no adjusted positive result.",14,700);
 out["figS04_original_null_families.svg"]=svg(1080,582,"Figure S4 | Original-data taxonomic nulls","Nominal versus corrected passes: all recorded endpoints and terms, not just favourable diet contrasts.",b,"BH for 336 and 500 tests; Holm IUT for 160. Exact p values and unsupported panel status tables accompany this figure.");
 }
 {let b=panel(24,101,"A","Change in deep-class residual coefficient")+panel(24,386,"B","Unresolved-origin partial identification");
 const [a,z]=ax(335,129,720,186,-.3,.35,[-.3,-.2,-.1,0,.1,.2,.3],"Coefficient / exact change component");b+=a;
 const labels=["All 30 islands","Island selection","Scored-count control change","Remove known introduced","Remove unresolved origin","Final literal native (20)"];
 labels.forEach((v,i)=>{let y=144+i*31;b+=t(320,y+3,v,13,400,"end");[D.origin.geo5,D.origin.climate10].forEach((a,j)=>{const val=a[i];b+=p(z(val),y+(j===0?-5:6),j===0?C.teal:C.blue,j===1)})});
 const [aa,g]=ax(335,414,720,135,-1.5,1.5,[-1.5,-1,-.5,0,.5,1,1.5],"Fixed-design residual identification bound");b+=aa;
 D.origin.bounds.forEach((a,i)=>{let y=433+i*35;b+=t(320,y+4,`≤${a[0]*100}% unresolved native`,13,400,"end")+l(g(a[1]),y,g(a[2]),y,C.orange,2)+l(g(a[1]),y-5,g(a[1]),y+5,C.orange)+l(g(a[2]),y-5,g(a[2]),y+5,C.orange)});
 b+=p(25,357,C.teal)+t(40,362,"Geographic sources",12)+p(220,357,C.blue,true)+t(235,362,"Climate-matched sources",12);
 b+=t(24,621,"0 / 460 feasible endpoint-wise bounds identify coefficient sign.",17,700,"start",C.orange)+t(24,644,"Bounds are not confidence intervals. Panel A is arithmetic accounting, not an intervention effect.",12,400,"start",C.sub);
 out["figS05_origin_uncertainty.svg"]=svg(1100,692,"Figure S5 | Native status and cohort selection","Same-species form/depth cohort, all-recorded min5 versus literal native; fixed-design sensitivity.",b,"Sources: original origin decomposition run 45385821139; native origin bounds run 45386194167.");
 }
 return out
};
for(const [name,svg] of Object.entries(render(data))){const p=path.join(dir,"supplement",name);if(process.argv.includes("--check")){if(fs.readFileSync(p,"utf8")!==svg)throw Error("Replay mismatch: "+name)}else{fs.mkdirSync(path.dirname(p),{recursive:true});fs.writeFileSync(p,svg)}console.log(name)}
