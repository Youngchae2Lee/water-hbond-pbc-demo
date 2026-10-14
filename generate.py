"""Construct four hand-designed water-shaped molecules, without research data."""
import json
import math
from pathlib import Path


def make_system():
    box=[12.0,12.0,12.0]
    oxygen=[[10.5,1.,1.],[1.2,1.,1.],[3.9,1.,1.],[3.9,4.4,1.]]
    directions=[([1,0,0],[0,1,0]),([1,0,0],[0,1,0]),
                ([0,1,0],[0,0,1]),([0,0,1],[1,0,0])]
    hydrogen=[]
    theta=math.radians(104.5); length=0.9572
    for o,(u,v) in zip(oxygen,directions):
        hydrogen.append([[ (o[k]+length*u[k]) % box[k] for k in range(3)],
                         [ (o[k]+length*(math.cos(theta)*u[k]+math.sin(theta)*v[k])) % box[k] for k in range(3)]])
    return {'description':'Hand-designed geometry; not an MD snapshot or force-field validation',
            'units':'angstrom','box':box,'oxygen':oxygen,'hydrogen':hydrogen}


if __name__ == '__main__':
    out=Path(__file__).parent/'data'; out.mkdir(exist_ok=True)
    (out/'synthetic_water.json').write_text(json.dumps(make_system(),indent=2)+'\n',encoding='utf-8',newline='\n')
