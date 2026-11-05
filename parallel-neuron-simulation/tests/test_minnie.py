"""Test anatomical identity, typed receptors, and sparse graph accounting."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from prepare_minnie import column_nodes
from minnie_network import graph_summary, receptor_components


class MinnieTests(unittest.TestCase):
    def test_annotation_selection_preserves_identity_and_converts_voxels(self):
        row={'pt_root_id':'864691135999999999','target_id':'12','classification_system':'aibs_coarse_excitatory','cell_type':'23P','pt_position_x':'250','pt_position_y':'500','pt_position_z':'25'}
        nodes,excluded=column_nodes([row])
        self.assertEqual(nodes[0]['root_id'],'864691135999999999')
        self.assertEqual([nodes[0][a+'_um'] for a in 'xyz'],[1.,2.,1.])
        self.assertEqual(excluded,[])

    def test_ambiguous_and_merged_identities_are_not_silently_used(self):
        def row(root,kind):return {'pt_root_id':root,'target_id':root,'classification_system':kind,'cell_type':'23P','pt_position_x':'1','pt_position_y':'1','pt_position_z':'1'}
        selected,excluded=column_nodes([row('1','aibs_coarse_unclear'),row('2','aibs_coarse_excitatory'),row('2','aibs_coarse_inhibitory')])
        self.assertEqual(selected,[])
        self.assertEqual(len(excluded),3)

    def test_sparse_reciprocal_component_and_isolated_cell(self):
        nodes=[{'root_id':root} for root in ['a','b','c']]
        edges=[{'pre_pt_root_id':'a','post_pt_root_id':'b','n_syn':'2'},{'pre_pt_root_id':'b','post_pt_root_id':'a','n_syn':'3'}]
        result=graph_summary(nodes,edges)
        self.assertEqual(result['anatomical_contacts'],5)
        self.assertEqual(result['weak_component_sizes'],[2,1])
        self.assertEqual(result['strong_component_sizes'],[2,1])
        self.assertEqual(result['isolated_neurons'],['c'])

    def test_observed_inhibitory_class_selects_inhibitory_receptors(self):
        exc={'ei_class':'E','manual_cell_type':'23P'}
        basket={'ei_class':'I','manual_cell_type':'BC'}
        martinotti={'ei_class':'I','manual_cell_type':'MC'}
        self.assertEqual(receptor_components(exc,exc),[('AMPA',.5),('NMDA',.5)])
        self.assertEqual(receptor_components(basket,exc),[('GABAA',1.)])
        self.assertEqual(receptor_components(martinotti,exc),[('GABAASlow',.9),('GABAB',.1)])
        self.assertEqual(receptor_components(martinotti,basket),[('GABAA',1.)])


if __name__=='__main__':unittest.main()
