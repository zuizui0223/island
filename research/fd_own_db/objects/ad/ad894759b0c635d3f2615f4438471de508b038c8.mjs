import fs from "node:fs";
import path from "node:path";
import {fileURLToPath} from "node:url";
const here=path.dirname(fileURLToPath(import.meta.url));
const dir=path.resolve(here,"../docs/manuscript_20261124");
const data=JSON.parse(fs.readFileSync(path.join(dir,"figure_data.json"),"utf8"));
function renderRegion(d){
 const C={ink:"#172B3A",sub:"#51616D",grid:"#DCE3E7",teal:"#007C78",blue:"#3C6FA3",orange:"#B96A34",gray:"#9BA8B1",white:"#FFFFFF"};
 const e=s=>String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;");
 const T=(x,y,s,o={})=>`<text x="${x}" y="${y}" font-size="${o.sz??15}" font-weight="${o.bold?700:400}" text-anchor="${o.anchor??"start"}" fill="${o.color??C.ink}">${e(s)}</text>`;
 const L=(x,y,x2,y2,o={})=>`<line x1="${x}" y1="${y}" x2="${x2}" y2="${y2}" stroke="${o.color??C.grid}" stroke-width="${o.w??1}"${o.dash?' stroke-dasharray="5 4"':""}/>`;
 const P=(x,y,c,r=5,open=false)=>`<circle cx="${x}" cy="${y}" r="${r}" fill="${open?C.white:c}" stroke="${c}" stroke-width="1.6"/>`;
 const fmt=(v,n=3)=>(v>0?"+":"")+v.toFixed(n);
 const scale=(lo,hi,x,w)=>v=>x+(v-lo)/(hi-lo)*w;
 const panel=(x,y,k,s)=>T(x,y,k,{sz:20,bold:true})+T(x+30,y,s,{sz:17,bold:true});
 function ax(x,y,w,h,lo,hi,ticks,lab){const f=scale(lo,hi,x,w);let s="";for(const t of ticks)s+=L(f(t),y,f(t),y+h,{color:t===0?C.sub:C.grid,w:t===0?1.2:1,dash:t===0})+T(f(t),y+h+21,t,{anchor:"middle",sz:12,color:C.sub});return [s+L(x,y+h,x+w,y+h,{color:C.sub})+T(x+w/2,y+h+42,lab,{anchor:"middle",sz:13}),f]}
 function point(f,y,v,ci,color,open=false){return (ci?L(f(ci[0]),y,f(ci[1]),y,{color,w:1.8})+L(f(ci[0]),y-4,f(ci[0]),y+4,{color})+L(f(ci[1]),y-4,f(ci[1]),y+4,{color}):"")+P(f(v),y,color,5,open)}
 const w=1280,h=1028;let b=panel(24,105,"A","Radial × deep across geographical regions")+panel(660,105,"B","Tropical deep-tube combinations");
 let [a,f]=ax(244,133,355,202,-5,15,[-5,0,5,10,15],"Island − candidate mainland joint frequency (pp)");b+=a;
 d.rows.forEach((v,i)=>{const y=154+i*48;b+=T(230,y-2,v[0],{anchor:"end",sz:13})+T(230,y+15,`n = ${v[1].toLocaleString("en-US")}`,{anchor:"end",sz:12,color:C.sub})+point(f,y,v[2],[v[3],v[4]],v[2]>0?C.teal:C.orange)+T(f(v[4])+6,y-10,fmt(v[2],2),{sz:11,bold:true})});
 let [a2,g]=ax(875,133,362,202,-5,15,[-5,0,5,10,15],"Joint frequency difference (pp)");b+=a2;
 d.tropical_pairs.forEach((v,i)=>{const y=169+i*67;b+=T(860,y-2,v.label,{anchor:"end",sz:13})+T(860,y+15,`n = ${v.n}`,{anchor:"end",sz:12,color:C.sub})+point(g,y,v.joint,v.joint_interval,v.joint>0?C.teal:C.orange)+T(v.joint>=0?g(v.joint)+7:g(v.joint)-7,y-12,fmt(v.joint,2),{sz:12,bold:true,anchor:v.joint>=0?"start":"end"})});
 b+=L(24,399,1256,399);
 b+=panel(24,432,"C","Tropical covariance and genus accounting")+panel(660,432,"D","Known-native recorded subset (96 islands)");
 const [a3,k]=ax(230,460,369,198,-5,6,[-4,-2,0,2,4,6],"Island − source covariance difference (pp)");b+=a3;
 d.tropical_pairs.forEach((v,i)=>{const y=488+i*63;b+=T(216,y+3,v.label,{anchor:"end",sz:12})+point(k,y-12,v.covariance,v.cov_interval,C.ink)+point(k,y+2,v.between,v.between_interval,C.teal)+point(k,y+16,v.within,v.within_interval,C.blue,true)});
 b+=P(29,724,C.ink)+T(42,728,"Total ΔC",{sz:12})+P(178,724,C.teal)+T(191,728,"Between genera",{sz:12})+P(366,724,C.blue,5,true)+T(379,728,"Within genera",{sz:12});
 const [a4,j]=ax(875,460,362,198,-10,20,[-10,-5,0,5,10,15,20],"Source-relative difference (pp)");b+=a4;
 d.native.rows.forEach((v,i)=>{const y=506+i*91;b+=T(860,y+2,v.label,{anchor:"end",sz:13})+point(j,y-13,v.joint,v.joint_interval,C.teal)+point(j,y+13,v.covariance,v.cov_interval,C.blue,true)+T(j(v.joint),y-25,fmt(v.joint,2),{anchor:"middle",sz:12,bold:true})+T(j(v.covariance),y+32,fmt(v.covariance,2),{anchor:"middle",sz:12,bold:true})});
 b+=P(680,724,C.teal)+T(693,728,"Joint ΔJ",{sz:12})+P(817,724,C.blue,5,true)+T(830,728,"Covariance ΔC",{sz:12})+T(660,752,"Nonrandom provenance-documented subset; not an independent replicate.",{sz:12,color:C.sub});
 b+=L(24,774,1256,774)+panel(24,806,"E","Source and documentation sensitivity of radial × deep");
 b+=T(24,833,"Source comparison",{sz:12,bold:true})+T(418,833,"Joint ΔJ (pp)",{sz:12,bold:true})+T(700,833,"Covariance ΔC (pp)",{sz:12,bold:true});
 const ir=(v,ci)=>`${fmt(v,2)} [${fmt(ci[0],2)}, ${fmt(ci[1],2)}]`;
 d.radial_sensitivity.forEach((v,i)=>{const y=857+i*24;b+=T(24,y,`${v.label} (n=${v.n})`,{sz:12})+T(418,y,ir(v.joint,v.joint_interval),{sz:12})+T(700,y,ir(v.covariance,v.cov_interval),{sz:12})});
 b+=T(24,966,"Animal-matched tropical subset (123 islands): geo5 +9.52 pp; climate10 +9.41 pp. All sensitivities use overlapping cohorts.",{sz:12,color:C.sub});
 const ttl="Figure 2 | Reassembly depends on geography and floral trait combinations";
 return `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" role="img" aria-label="${e(ttl)}"><rect width="100%" height="100%" fill="#FFFFFF"/><g font-family="Arial, Helvetica, sans-serif">${T(24,32,ttl,{sz:23,bold:true})}${T(24,55,"Species-Direct paired categories relative to candidate mainland pools. Panels A–D use geo5.",{sz:13,color:C.sub})}${L(24,69,w-24,69)}${b}${L(24,h-34,w-24,h-34)}${T(24,h-14,"Whiskers: archived conditional 95% geographical-block intervals. Panel E includes alternative source-group intervals. No causal effects implied.",{sz:11,color:C.sub})}</g></svg>\n`
}
const buildFigures=function buildFigures(D){
 const C={ink:"#172B3A",sub:"#51616D",grid:"#DCE3E7",teal:"#007C78",blue:"#3C6FA3",orange:"#B96A34",light:"#EAF3F2",gray:"#9BA8B1",black:"#172B3A",white:"#FFFFFF"};
 const e=s=>String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;");
 const T=(x,y,s,o={})=>`<text x="${x}" y="${y}" font-size="${o.sz??15}" font-weight="${o.bold?700:400}" text-anchor="${o.anchor??"start"}" fill="${o.color??C.ink}"${o.italic?' font-style="italic"':""}>${e(s)}</text>`;
 const L=(x,y,x2,y2,o={})=>`<line x1="${x}" y1="${y}" x2="${x2}" y2="${y2}" stroke="${o.color??C.grid}" stroke-width="${o.w??1}"${o.dash?' stroke-dasharray="5 4"':""}/>`;
 const R=(x,y,w,h,fill,rx=0)=>`<rect x="${x}" y="${y}" width="${Math.max(0,w)}" height="${h}" rx="${rx}" fill="${fill}"/>`;
 const P=(x,y,c,r=5,open=false)=>`<circle cx="${x}" cy="${y}" r="${r}" fill="${open?C.white:c}" stroke="${c}" stroke-width="1.6"/>`;
 const fmt=(v,d=3)=>(v>0?"+":"")+v.toFixed(d);
 const N=v=>v.toLocaleString("en-US");
 const scale=(lo,hi,x,w)=>v=>x+(v-lo)/(hi-lo)*w;
 function chartAxis(x,y,w,h,lo,hi,ticks,lab){
   const f=scale(lo,hi,x,w);let s="";
   for(const t of ticks){s+=L(f(t),y,f(t),y+h,{color:t===0?C.sub:C.grid,w:t===0?1.2:1,dash:t===0})+T(f(t),y+h+22,t.toString(),{anchor:"middle",sz:13,color:C.sub})}
   s+=L(x,y+h,x+w,y+h,{color:C.sub})+T(x+w/2,y+h+45,lab,{anchor:"middle",sz:14});return [s,f];
 }
 function dotRange(f,y,v,range,color=C.teal,open=false){let s="";if(range){const[a,b]=range;s+=L(f(a),y,f(b),y,{color,w:1.8})+L(f(a),y-4,f(a),y+4,{color,w:1.5})+L(f(b),y-4,f(b),y+4,{color,w:1.5})}return s+P(f(v),y,color,5,open)}
 function panel(x,y,k,s){return T(x,y,k,{sz:20,bold:true})+T(x+30,y,s,{sz:17,bold:true})}
 function svg(w,h,title,subtitle,body,foot){
   return `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" role="img" aria-label="${e(title)}"><rect width="100%" height="100%" fill="#FFFFFF"/><g font-family="Arial, Helvetica, sans-serif">${T(24,32,title,{sz:23,bold:true})}${T(24,55,subtitle,{sz:13,color:C.sub})}${L(24,69,w-24,69)}${body}${L(24,h-34,w-24,h-34)}${T(24,h-14,foot,{sz:11,color:C.sub})}</g></svg>\n`
 }
 const out={};
 {const d=D.recurrence;let b=panel(24,101,"A","Botanical support")+panel(390,101,"B","Recurrence of recorded colour × form");
  const cards=[["GIFT island floras",N(d.parent)],["Eligible ≥5 trait names",N(d.eligible)],["Island pairs",N(d.rows[0][1])]];
  cards.forEach((v,i)=>{const y=124+i*85;b+=R(24,y,322,69,i===1?C.light:"#F2F5F7",5)+T(38,y+24,v[0],{sz:14,color:C.sub})+T(38,y+52,v[1],{sz:25,bold:true})});
  let [a,f]=chartAxis(654,128,350,245,0,.5,[0,.1,.2,.3,.4,.5],"Mean total-variation overlap");b+=a;
  d.rows.forEach((r,i)=>{const y=165+i*73;b+=T(635,y-2,r[0],{anchor:"end",sz:15})+T(635,y+18,`${N(r[1])} pairs`,{anchor:"end",sz:12,color:C.sub})+R(f(0),y-13,f(r[2])-f(0),25,i===0?C.blue:C.teal,2)+T(f(r[2])+8,y+4,r[2].toFixed(3),{sz:14,bold:true})});
  b+=T(24,442,"Strict genus-disjoint pairs: observed squared distance 0.281 vs shuffled-trait reference 0.257.",{sz:14,color:C.sub})+T(24,465,"Nonzero overlap is recurrence; it does not establish excess convergence.",{sz:14,bold:true});
  out["fig01_recurrence.svg"]=svg(1060,510,"Figure 1 | Floral repertoires recur despite taxonomic turnover","Observed high/medium species-Direct colour × form; full-flora taxonomic disjointness.",b,"Source: fd verified positive evidence ledger; 24 Nov 2026 synthesis. Pair cohorts overlap; no independence assumed.");
 }
 out["fig02_geographical_reassembly.svg"]=renderRegion(D.region);
 {const d=D.geography;let b=panel(24,101,"A","Independent source support")+panel(450,101,"B","Isolation-distance sensitivity");
  d.source_counts.forEach((r,i)=>{const y=142+i*52;b+=T(24,y,r[0],{sz:14})+R(24,y+9,310*r[1]/305,12,C.blue,2)+T(345,y+20,N(r[1]),{sz:14,bold:true,anchor:"end"})});
  b+=T(24,377,"38 tropical bird–bee candidates",{sz:17,bold:true})+T(24,399,"7 archipelagos; 6 informative",{sz:13,color:C.sub})+T(24,421,"9 tropical three-guild candidates",{sz:13,color:C.sub});
  let [a,f]=chartAxis(660,130,365,229,-.12,.32,[-.1,0,.1,.2,.3],"Slope difference per doubling of 1 + distance");b+=a;
  d.distances.forEach((r,i)=>{const y=168+i*70;b+=T(644,y+4,r[0],{anchor:"end",sz:13})+dotRange(f,y-7,r[1],null,C.teal)+dotRange(f,y+10,r[2],null,C.blue,true)});
  b+=P(470,418,C.teal)+T(482,423,"Nectar >0 minus bees",{sz:13})+P(716,418,C.blue,5,true)+T(728,423,"Nectar >0 minus nonnectar birds",{sz:13});
  b+=L(24,445,1036,445)+T(24,472,"Mainland bootstrap ranges include zero: bee contrast [−0.005, 2.127]; nonnectar [−0.129, 1.595].",{sz:14,color:C.sub})+T(24,495,"For ≥50% diet, the mainland contrast against bees is −0.0015. No robust guild replacement result.",{sz:14,bold:true});
  out["fig03_partner_geography.svg"]=svg(1060,540,"Figure 3 | Geographical differentiation of potential partners","Provisional bird island matches; positive butterfly records do not measure complete native richness.",b,"Distance coefficients are from separate correlated-distance models. Only mainland exploratory resampling ranges are available here.");
 }
 {const d=D.diet;let b=panel(24,101,"A","Diet-class associations (32 islands)")+panel(540,101,"B","Tubular contrast across nested panels");
  let [a,f]=chartAxis(219,127,268,150,-.10,.16,[-.1,0,.1],"Tubular-form coefficient");b+=a;
  d.coefs.forEach((r,i)=>{const y=153+i*50;b+=T(205,y+4,r[0],{anchor:"end",sz:13})+R(Math.min(f(0),f(r[1])),y-10,Math.abs(f(r[1])-f(0)),20,i===0?C.teal:i===1?C.orange:C.gray,2)+T(r[1]>0?f(r[1])+6:f(r[1])-5,y-15,fmt(r[1]),{anchor:r[1]>0?"start":"end",sz:12,bold:true})});
  b+=T(24,343,"Specialist − partial = +0.187",{sz:17,bold:true,color:C.teal})+T(24,364,"Six informative groups; 5-name threshold.",{sz:13,color:C.sub});
  let [a2,g]=chartAxis(820,125,395,422,-.18,.9,[-.1,0,.2,.4,.6,.8],"Specialist − partial coefficient");b+=a2;
  const labels={all:"All matches",tight:"Tight matches",exclude_shared_boundary:"Exclude boundaries"};
  d.contrasts.forEach((r,i)=>{let y=142+i*34;let m=r.id.match(/^(direct_single_state|direct)_n(5|50)_(.*)$/);if(!m)throw Error(r.id);let prefix=m[1]==="direct"?"Direct":"Single";let title=prefix+" ≥"+m[2]+" · "+labels[m[3]];b+=T(806,y+4,title+` (n=${r.n})`,{anchor:"end",sz:11.5,color:r.estimate===null?C.gray:C.ink});if(r.estimate===null)b+=T(g(0)+6,y+4,"Unsupported",{sz:12,color:C.sub});else b+=dotRange(g,y,r.estimate,r.range,C.teal,r.range===null)});
  b+=T(540,610,"11/11 estimable contrasts > 0; one panel unsupported.",{sz:14,bold:true})+T(540,630,"Whiskers: group-resampling percentile ranges; open = range unavailable.",{sz:12,color:C.sub});
  b+=panel(24,405,"C","Bird dietary null: tubular form");
  b+=T(24,433,"Plant scope",{sz:12,bold:true})+T(282,433,"Family q",{sz:12,bold:true,anchor:"end"})+T(475,433,"Genus q",{sz:12,bold:true,anchor:"end"})+L(24,441,487,441);
  d.tests.forEach((r,i)=>{let y=465+i*31;b+=T(24,y,`${r[0]} (n=${r[1]})`,{sz:13})+T(282,y,r[2].toFixed(3),{anchor:"end",sz:13,bold:r[2]<.05,color:r[2]<.05?C.teal:C.ink})+T(475,y,r[3].toFixed(3),{anchor:"end",sz:13})});
  b+=T(24,599,"6 / 24 tests pass BH in the family null",{sz:14,bold:true})+T(24,620,"(3 tubular, 3 extended); 0 in the genus null.",{sz:13,color:C.sub});
  b+=L(24,650,1232,650)+panel(24,681,"D","Held-out within-archipelago prediction (tubular form)");
  let headers=["Scope","Broad nectar MSE","Split diet MSE","Groups improved"];[24,385,665,960].forEach((x,i)=>b+=T(x,709,headers[i],{sz:13,bold:true}));
  d.tests.forEach((r,i)=>{let y=734+i*25;b+=T(24,y,`${r[0]} (n=${r[1]})`,{sz:13})+T(385,y,r[4].toFixed(5),{sz:13})+T(665,y,r[5].toFixed(5),{sz:13,bold:r[5]<r[4],color:r[5]<r[4]?C.teal:C.orange})+T(960,y,`${r[7]}/${r[6]}`,{sz:13})});
  out["fig04_dietary_differentiation.svg"]=svg(1260,865,"Figure 4 | Bird dietary classes differ in observed floral correspondence","Historical public Direct name-unit aggregates. Dietary fraction refers to nectar, pollen and exudates.",b,"Nested exploratory panels and provisional geography; group-resampling ranges are not calibrated confidence intervals; no causal matching inferred.");
 }
 {const d=D.dual;let b=panel(24,101,"A","Same plants and same 30 islands")+panel(625,101,"B","Residual sensitivity to candidate source");
  let [a,f]=chartAxis(238,133,335,222,-.12,.38,[-.1,0,.1,.2,.3],"Specialist − partial coefficient");b+=a;
  d.rows.forEach((r,i)=>{const y=167+i*75;b+=T(224,y+4,r[0],{anchor:"end",sz:13})+dotRange(f,y-9,r[1],null,C.teal)+dotRange(f,y+9,r[2],null,C.gray,true)});
  b+=P(26,419,C.teal)+T(38,423,"Observed",{sz:13})+P(151,419,C.gray,5,true)+T(163,423,"Genus × grade expected",{sz:13});
  let [a2,g]=chartAxis(818,133,385,222,-.12,.24,[-.1,0,.1,.2],"Observed − taxonomy expectation");b+=a2;
  d.rows.forEach((r,i)=>{let y=167+i*75;b+=T(803,y+4,r[0],{anchor:"end",sz:13})+dotRange(g,y-9,r[3],null,C.teal)+dotRange(g,y+9,r[4],null,C.blue,true)});
  b+=P(650,419,C.teal)+T(662,423,"Geographic sources",{sz:13})+P(860,419,C.blue,5,true)+T(872,423,"Climate-matched sources",{sz:13});
  b+=L(24,446,1236,446)+panel(24,478,"C","Deep-class residual: native-status support")+panel(760,478,"D","Two-sided evidence gate");
  let [a3,k]=chartAxis(339,502,350,176,-.2,.25,[-.2,-.1,0,.1,.2],"Source-relative residual coefficient");b+=a3;
  d.native.forEach((r,i)=>{const y=530+i*55;b+=T(325,y+4,r[0],{anchor:"end",sz:12})+dotRange(k,y-6,r[1],null,C.teal)+dotRange(k,y+8,r[2],null,C.blue,true)});
  b+=T(760,523,"160 directional intersection-union tests",{sz:15,bold:true})+T(760,553,"4 unadjusted passes; same depth cohort",{sz:14})+T(760,586,"0 after Holm correction",{sz:20,bold:true,color:C.orange})+T(760,619,"0 / 460 feasible origin bounds identify sign",{sz:14,bold:true})+T(760,645,"Positive depth residual disappears before",{sz:13,color:C.sub})+T(760,663,"restricting plant origin on common islands.",{sz:13,color:C.sub});
  out["fig05_two_sided_specificity.svg"]=svg(1260,765,"Figure 5 | Two-sided taxonomy and origin constrain floral specificity","Original protected species-Direct form/depth cohorts; points are coefficient diagnostics, not fitness effects.",b,"Genus-grade expectation is statistical accounting. No adjusted dual-null support or stable native effect. No independent replication claimed.");
 }
 return out
};
for(const [name,svg] of Object.entries(buildFigures(data))){const p=path.join(dir,"figures",name);if(process.argv.includes("--check")){if(fs.readFileSync(p,"utf8")!==svg)throw Error("Non-reproducible: "+name);}else{fs.mkdirSync(path.dirname(p),{recursive:true});fs.writeFileSync(p,svg);}console.log(name);}
