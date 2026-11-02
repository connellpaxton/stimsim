"""Check lossless skeleton cleanup and geometry-dependent cable construction."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

TOOLS=Path(__file__).resolve().parents[1]/'tools'
sys.path.insert(0,str(TOOLS))
from prepare_minnie_skeletons import collapse_identical_vertices, repair_zero_radii
from reconstructed_cell import Curve


class SkeletonTests(unittest.TestCase):
    def test_duplicate_vertex_collapse_preserves_branch_endpoint_and_length(self):
        rows=np.array([[0,1,0,0,0,5,-1],[1,3,10,0,0,1,0],[2,3,10,0,0,1,1],[3,3,20,0,0,.8,2]],dtype=float)
        clean,parent,order,changes=collapse_identical_vertices(rows,np.array([-1,0,1,2]),np.arange(4))
        self.assertEqual(clean[:,0].tolist(),[0,1,3])
        self.assertEqual(clean[-1,6],1)
        self.assertEqual(parent.tolist(),[-1,0,1])
        self.assertEqual(len(changes),1)
        self.assertEqual(sum(np.linalg.norm(clean[i,2:5]-clean[p,2:5]) for i,p in enumerate(parent) if p>=0),20.)

    def test_zero_radius_repair_leaves_measured_positive_values_unchanged(self):
        rows=np.array([[0,1,0,0,0,5,-1],[1,3,10,0,0,.2,0],[2,3,20,0,0,0,1],[3,3,30,0,0,.1,2]],dtype=float)
        radii,repairs=repair_zero_radii(rows,np.array([-1,0,1,2]),np.arange(4))
        np.testing.assert_array_equal(radii,[5,.2,.2,.1])
        self.assertEqual(repairs[0]['measured_ancestor_node_id'],1)
        self.assertEqual(rows[2,5],0)

    def test_reference_curve_combines_same_distance_and_clamps_endpoints(self):
        curve=Curve([0,0,10],[1,3,5])
        np.testing.assert_allclose(curve.at([-1,0,5,10,11]),[2,2,3.5,5,5])


@unittest.skipUnless(importlib.util.find_spec('neuron'),'Requires NEURON')
class CableTests(unittest.TestCase):
    def test_actual_branched_swc_uses_distinct_compartments_and_stays_at_rest(self):
        script='''
import os
from neuron import h
from reconstructed_cell import ReconstructedCell, ElectricalProfile
h.load_file('stdrun.hoc');h.load_file('import3d.hoc')
points=[]
for region,g in [(1,.0001),(2,.0001),(3,.0002),(4,.0004)]:
 points.append({'region':region,'distance_um':0.,'Ra_ohm_cm':100.,'cm_uf_cm2':1.,'v_init_mv':-65.,'mechanisms':{'pas':{'g':g,'e':-65.}},'ions':{},'ion_styles':{}})
profile=ElectricalProfile({'points':points});pc=h.ParallelContext()
node={'root_id':'864691135999999999','profile':'passive_probe','cell_model':'passive','cell_type':'probe','manual_cell_type':'23P','ei_class':'E'}
cell=ReconstructedCell(os.environ['NEURONSIM_TEST_SWC'],node,0,profile,pc)
assert len(cell.secLists['spiny'])==2
assert abs(cell.secs['dend_0']['hObj'](.5).g_pas-.0002)<1e-12
assert abs(cell.secs['apic_0']['hObj'](.5).g_pas-.0004)<1e-12
assert cell.secs['apic_0']['hObj'].L>90
assert all(info['hObj'].nseg%2==1 for info in cell.secs.values())
pc.set_maxstep(10);h.finitialize(-65);pc.psolve(.1)
assert all(abs(seg.v+65)<1e-8 for info in cell.secs.values() for seg in info['hObj'])
'''
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'branched.swc'
            path.write_text('0 1 0 0 0 5 -1\n1 3 0 5 0 1 0\n2 3 0 100 0 .5 1\n3 4 5 0 0 1 0\n4 4 100 0 0 .5 3\n5 2 0 -5 0 .5 0\n6 2 0 -40 0 .5 5\n')
            env=dict(os.environ,NEURON_MODULE_OPTIONS='-nogui',PYTHONPATH=str(TOOLS),NEURONSIM_TEST_SWC=str(path))
            result=subprocess.run([sys.executable,'-c',script],env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)


if __name__=='__main__':unittest.main()
