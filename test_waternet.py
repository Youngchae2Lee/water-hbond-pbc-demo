import copy
import itertools
import math
import random
import unittest
from generate import make_system
from waternet import analyze, minimum_image, norm


class WaterTests(unittest.TestCase):
    def setUp(self): self.s=make_system()
    def run_s(self, s=None, **kw):
        s=self.s if s is None else s
        return analyze(s['oxygen'],s['hydrogen'],s['box'],**kw)

    def test_known_path(self):
        r=self.run_s()
        self.assertEqual(r['adjacency'],[[0,1,0,0],[1,0,1,0],[0,1,0,1],[0,0,1,0]])
        self.assertEqual(r['directed'],[[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0]])
        self.assertEqual(r['degrees'],[1,2,2,1])
        self.assertEqual(r['degree_distribution'],{'1':.5,'2':.5})

    def test_boundary_pair(self):
        self.assertAlmostEqual(norm(minimum_image([-9.3,0,0],[12,12,12])),2.7)
        self.assertGreater(9.3,3.5)
        self.assertEqual(self.run_s()['adjacency'][0][1],1)

    def test_mic_against_image_enumeration(self):
        rng=random.Random(572)
        for _ in range(100):
            box=[8.,11.,14.]; d=[rng.uniform(-L,L) for L in box]
            brute=min(norm([d[k]+image[k]*box[k] for k in range(3)]) for image in itertools.product([-1,0,1],repeat=3))
            self.assertAlmostEqual(norm(minimum_image(d,box)),brute,places=11)

    def test_translation_and_independent_wrapping(self):
        s=copy.deepcopy(self.s)
        s['oxygen']=[[x+17 for x in o] for o in s['oxygen']]
        s['hydrogen']=[[[x+17 for x in h] for h in pair] for pair in s['hydrogen']]
        for o in s['oxygen']: o[0]+=24
        for pair in s['hydrogen']:
            for h in pair: h[1]-=36
        self.assertEqual(self.run_s(s),self.run_s())

    def test_permutation(self):
        order=[2,0,3,1]; s=copy.deepcopy(self.s)
        for key in ['oxygen','hydrogen']: s[key]=[s[key][i] for i in order]
        a=self.run_s()['adjacency']
        self.assertEqual(self.run_s(s)['adjacency'],[[a[i][j] for j in order] for i in order])

    def test_cutoff_controls(self):
        edges=[self.run_s(cutoff=c)['edges'] for c in [2.5,2.8,3.2,3.5,4.5]]
        self.assertEqual(edges,[0,2,2,3,3])
        self.assertEqual(self.run_s(cutoff=4.5)['distance_edges'],4)

    def test_distance_boundary_inclusive(self):
        self.assertEqual(self.run_s(cutoff=2.7)['edges'],2)
        self.assertEqual(self.run_s(cutoff=2.7-1e-8)['edges'],0)

    def test_cutoff_tolerance_does_not_admit_half_box_tie(self):
        s={'box':[12.,12.,12.], 'oxygen':[[0.,1.,1.],[6.,1.,1.]],
           'hydrogen':[[[1.,1.,1.],[0.,2.,1.]],
                       [[5.,1.,1.],[6.,2.,1.]]]}
        r=self.run_s(s,cutoff=math.nextafter(6.0,0.0))
        self.assertEqual(r['edges'],0)
        self.assertEqual(r['distance_edges'],0)

    def test_angle_boundary_and_rejection(self):
        s={'box':[12.,12.,12.],'oxygen':[[1.,1.,1.],[4.,1.,1.]],
           'hydrogen':[[[1+math.cos(math.pi/6),1+math.sin(math.pi/6),1],[1,1,2]],[[4,1,2],[4,2,1]]]}
        self.assertEqual(self.run_s(s,angle_degrees=30)['edges'],1)
        self.assertEqual(self.run_s(s,angle_degrees=29.9)['edges'],0)

    def test_simple_graph_invariants(self):
        r=self.run_s(); a=r['adjacency']
        self.assertTrue(all(a[i][i]==0 for i in range(4)))
        self.assertTrue(all(a[i][j]==a[j][i] for i in range(4) for j in range(4)))
        self.assertEqual(sum(r['degrees']),2*r['edges'])
        self.assertEqual(sum(r['degree_distribution'].values()),1)

    def test_invalid_inputs(self):
        for cutoff in [0,-1,6,float('nan')]:
            with self.assertRaises(ValueError): self.run_s(cutoff=cutoff)
        for angle in [-1,181,float('inf')]:
            with self.assertRaises(ValueError): self.run_s(angle_degrees=angle)
        for box in [[0,12,12],[12,12],[12,float('nan'),12]]:
            s=copy.deepcopy(self.s); s['box']=box
            with self.assertRaises(ValueError): self.run_s(s)

    def test_overlap_and_zero_oh(self):
        s=copy.deepcopy(self.s); s['oxygen'][1]=s['oxygen'][0][:]
        with self.assertRaises(ValueError): self.run_s(s)
        s=copy.deepcopy(self.s); s['hydrogen'][0][0]=s['oxygen'][0][:]
        with self.assertRaises(ValueError): self.run_s(s)

    def test_second_hydrogen_can_donate(self):
        # Only the second hydrogen of molecule 0 points at molecule 1.
        s={'box':[12.,12.,12.], 'oxygen':[[1.,1.,1.],[4.,1.,1.]],
           'hydrogen':[[[1.,2.,1.],[2.,1.,1.]],
                       [[4.,2.,1.],[4.,1.,2.]]]}
        r=self.run_s(s)
        self.assertEqual(r['directed'], [[0,1],[0,0]])
        self.assertEqual(r['degrees'], [1,1])

    def test_bidirectional_donation_is_one_pair_edge(self):
        # Two directed donations still represent only one molecular pair.
        s={'box':[12.,12.,12.], 'oxygen':[[1.,1.,1.],[4.,1.,1.]],
           'hydrogen':[[[2.,1.,1.],[1.,2.,1.]],
                       [[3.,1.,1.],[4.,2.,1.]]]}
        r=self.run_s(s)
        self.assertEqual(r['directed'], [[0,1],[1,0]])
        self.assertEqual(r['edges'], 1)

    def test_single_molecule_and_empty_system(self):
        s=copy.deepcopy(self.s)
        s['oxygen']=s['oxygen'][:1]; s['hydrogen']=s['hydrogen'][:1]
        self.assertEqual(self.run_s(s)['degree_distribution'], {'0':1.0})
        s['oxygen']=[]; s['hydrogen']=[]
        with self.assertRaises(ValueError): self.run_s(s)


if __name__ == '__main__': unittest.main()
