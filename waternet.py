"""Independent synthetic-water geometry demo; Python standard library only."""
import math
from collections import Counter


def vector3(value, name):
    if len(value) != 3 or not all(math.isfinite(x) for x in value):
        raise ValueError(name + ' must contain three finite numbers')
    return tuple(value)


def minimum_image(delta, box):
    """Orthorhombic MIC in [-L/2,L/2); ties map to the negative image."""
    vector3(delta, 'delta')
    vector3(box, 'box')
    if any(x <= 0 for x in box):
        raise ValueError('Box lengths must be positive')
    return tuple(d - L * math.floor(d / L + 0.5) for d, L in zip(delta, box))


def displacement(a, b, box):
    """Return the minimum-image vector pointing from a to b."""
    return minimum_image(tuple(y-x for x, y in zip(a, b)), box)


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def analyze(oxygen, hydrogen, box, cutoff=3.5, angle_degrees=30.0):
    """Coordinates/cutoff in angstrom; angle H-D-A at donor (0 is linear).

    directed[i][j]: i donates through either H to j.
    adjacency: logical OR of directed edges; counts molecule pairs, not H atoms.
    distance_adjacency: distance-only control, not a hydrogen-bond assignment.
    """
    vector3(box, 'box')
    if any(x <= 0 for x in box):
        raise ValueError('Nonpositive box')
    half_shortest_box = min(box) / 2
    if not math.isfinite(cutoff) or not 0 < cutoff < half_shortest_box:
        raise ValueError('Require 0 < cutoff < half the shortest box length')
    if not math.isfinite(angle_degrees) or not 0 <= angle_degrees <= 180:
        raise ValueError('Angle must be between 0 and 180 degrees')
    n=len(oxygen)
    if n == 0 or len(hydrogen) != n:
        raise ValueError('Require N oxygens and N hydrogen pairs')
    # Each molecule explicitly owns two H atoms; build their O-H unit vectors.
    bonds=[]
    for o, pair in zip(oxygen, hydrogen):
        vector3(o, 'oxygen')
        if len(pair) != 2:
            raise ValueError('Exactly two hydrogens per molecule')
        bh=[]
        for h in pair:
            vector3(h, 'hydrogen')
            v=displacement(o,h,box)
            length=norm(v)
            if not 0 < length < half_shortest_box:
                raise ValueError('Invalid O-H vector')
            bh.append(tuple(x/length for x in v))
        bonds.append(bh)
    directed=[[0]*n for _ in range(n)]
    distance=[[0]*n for _ in range(n)]
    cos_limit=math.cos(math.radians(angle_degrees))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            # i is the donor, j the acceptor; zero angle means O-H points at j.
            v=displacement(oxygen[i],oxygen[j],box)
            r=norm(v)
            if r <= 1e-12:
                raise ValueError('Distinct oxygens overlap under PBC')
            # Numerical tolerance must not admit an ambiguous half-box pair.
            if r >= half_shortest_box:
                continue
            if r <= cutoff + 1e-12:
                distance[i][j]=1
                if any(sum(a*b for a,b in zip(u,v))/r >= cos_limit-1e-12 for u in bonds[i]):
                    directed[i][j]=1
    # Merge donor directions into one undirected molecular-pair edge.
    adjacency=[[int(directed[i][j] or directed[j][i]) for j in range(n)] for i in range(n)]
    degrees=[sum(row) for row in adjacency]
    counts=Counter(degrees)
    return {'directed':directed,'adjacency':adjacency,'distance_adjacency':distance,
            'degrees':degrees,'degree_distribution':{str(k):counts[k]/n for k in sorted(counts)},
            'edges':sum(degrees)//2,'distance_edges':sum(map(sum,distance))//2}
