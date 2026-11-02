"""Data-only utilities for an observed Minnie column and approximate cell fits."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from neuronsim.graph import graph_summary

PROFILES = {
    '23P': ('IT2_reduced','IT','HH_reduced'),
    '4P': ('IT4_reduced','IT','HH_reduced'),
    '5P-IT': ('IT5B_reduced','IT','HH_reduced'),
    '5P-NP': ('IT5B_reduced','IT','HH_reduced'),
    '5P-PT': ('PT5B_full','PT','HH_full'),
    '6P-IT': ('IT6_reduced','IT','HH_reduced'),
    '6P-CT': ('CT6_reduced','CT','HH_reduced'),
    '6P-U': ('IT6_reduced','IT','HH_reduced'),
    'WM-P': ('IT6_reduced','IT','HH_reduced'),
    'Unsure E': ('IT5B_reduced','IT','HH_reduced'),
    'BC': ('PV_simple','PV','HH_simple'),
    'MC': ('SOM_simple','SOM','HH_simple'),
    'BPC': ('SOM_simple','SOM','HH_simple'),
    'NGC': ('SOM_simple','SOM','HH_simple'),
    'Unsure I': ('SOM_simple','SOM','HH_simple'),
}
APPROXIMATE_TYPES = {'5P-NP','6P-U','WM-P','Unsure E','BPC','NGC','Unsure I'}


def receptor_components(pre, post):
    """Use the published M1 receptor families and its E/SOM mixing fractions."""
    if pre['ei_class']=='E':
        return [('AMPA',.5),('NMDA',.5)]
    if pre['manual_cell_type']!='BC' and post['ei_class']=='E':
        return [('GABAASlow',.9),('GABAB',.1)]
    return [('GABAA',1.)]
