"""Regenerate synthetic input, CSV/JSON results and self-contained SVG figures."""
import csv
import json
from pathlib import Path
from generate import make_system
from waternet import analyze

ROOT=Path(__file__).parent


def svg(name, title, body):
    text='<svg xmlns="http://www.w3.org/2000/svg" width="720" height="440" viewBox="0 0 720 440">'
    text+='<rect width="720" height="440" fill="#f8fafc"/><g font-family="Arial,sans-serif" fill="#0f172a">'
    text+=f'<text x="36" y="38" font-size="22">{title}</text>'+body+'</g></svg>'
    (ROOT/'results'/name).write_text(text,encoding='utf-8',newline='\n')


def main():
    for folder in ['data','results']: (ROOT/folder).mkdir(exist_ok=True)
    s=make_system()
    (ROOT/'data'/'synthetic_water.json').write_text(json.dumps(s,indent=2)+'\n',encoding='utf-8',newline='\n')
    runs=[]
    for cutoff in [2.5,2.8,3.2,3.5,4.5]:
        r=analyze(s['oxygen'],s['hydrogen'],s['box'],cutoff)
        runs.append({'cutoff_angstrom':cutoff,**r})
    (ROOT/'results'/'analysis.json').write_text(json.dumps(runs,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'results'/'cutoff_sweep.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f, lineterminator='\n'); w.writerow(['cutoff_angstrom','hbond_pair_edges','distance_only_edges','mean_degree'])
        w.writerows((r['cutoff_angstrom'],r['edges'],r['distance_edges'],sum(r['degrees'])/4) for r in runs)
    base=runs[3]
    with (ROOT/'results'/'adjacency.csv').open('w',newline='',encoding='utf-8') as f:
        csv.writer(f, lineterminator='\n').writerows(base['adjacency'])
    # Adjacency matrix: hand-checkable four-node path.
    body='<text x="36" y="70" font-size="15">3.5 A, 30 degrees | molecule pairs, undirected</text>'
    for i,row in enumerate(base['adjacency']):
        body+=f'<text x="176" y="{140+i*60}" font-size="18">{i}</text><text x="{213+i*60}" y="100" font-size="18">{i}</text>'
        for j,x in enumerate(row):
            body+=f'<rect x="{194+j*60}" y="{110+i*60}" width="56" height="56" fill="{ "#0f766e" if x else "#e2e8f0"}"/><text x="{216+j*60}" y="{144+i*60}" fill="{ "white" if x else "#334155"}" font-size="20">{x}</text>'
    body+='<text x="36" y="402" font-size="16">Expected edges: 0-1 (across boundary), 1-2, 2-3</text>'
    svg('adjacency.svg','Synthetic water | adjacency matrix',body)
    body='<text x="36" y="72" font-size="15">3.5 A, 30 degrees | degrees = [1, 2, 2, 1]</text>'
    body+='<path d="M90 100 V340 H650" fill="none" stroke="#64748b"/>'
    for k in range(4):
        p=base['degree_distribution'].get(str(k),0); x=135+k*125; height=360*p
        body+=f'<rect x="{x}" y="{340-height}" width="65" height="{height}" fill="#0f766e"/><text x="{x+23}" y="365">{k}</text><text x="{x+15}" y="{325-height}">{p:.2f}</text>'
    body+='<text x="300" y="405" font-size="16">Degree k</text><text x="20" y="98" font-size="16">P(k)</text>'
    svg('degree_distribution.svg','Synthetic water | degree distribution',body)
    body='<text x="36" y="72" font-size="15">Teal: geometric H bonds | orange: distance-only control</text>'
    body+='<path d="M90 100 V340 H650" fill="none" stroke="#64748b"/>'
    for y in range(5): body+=f'<text x="63" y="{345-y*50}">{y}</text><path d="M90 {340-y*50} H650" stroke="#e2e8f0"/>'
    for key,color in [('distance_edges','#ea580c'),('edges','#0f766e')]:
        points=' '.join(f'{120+(r["cutoff_angstrom"]-2.5)*240},{340-r[key]*50}' for r in runs)
        body+=f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="3"/>'
        for r in runs: body+=f'<circle cx="{120+(r["cutoff_angstrom"]-2.5)*240}" cy="{340-r[key]*50}" r="5" fill="{color}"/>'
    for r in runs: body+=f'<text x="{108+(r["cutoff_angstrom"]-2.5)*240}" y="367">{r["cutoff_angstrom"]}</text>'
    body+='<text x="245" y="410">Cutoff (angstrom; sampled settings)</text><text x="20" y="98">Edges</text>'
    svg('cutoff_sweep.svg','Synthetic water | cutoff sensitivity',body)
    print('Synthetic demo complete. At 3.5 A:',base['edges'],'edges; degrees',base['degrees'])


if __name__ == '__main__': main()
