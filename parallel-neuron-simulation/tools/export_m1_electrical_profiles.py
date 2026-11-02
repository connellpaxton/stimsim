"""Sample published M1 density parameters for transfer onto recorded skeletons."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import runpy
import sys

from minnie_network import PROFILES

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'vendor/m1_260015'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=args.output.resolve();output.parent.mkdir(parents=True,exist_ok=True)
    os.environ['NEURON_MODULE_OPTIONS']=os.environ.get('NEURON_MODULE_OPTIONS','')+' -nogui'
    os.environ['MPLBACKEND']='Agg';os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
    import neuron
    from neuron import h
    from netpyne import sim
    manifest=json.loads((SOURCE/'SOURCE_MANIFEST.json').read_text())
    for item in manifest['files']:
        if hashlib.sha256((SOURCE/item['path']).read_bytes()).hexdigest()!=item['sha256']:raise RuntimeError('Published source changed')
    directory=SOURCE/'sim';neuron.load_mechanisms(str(directory));os.chdir(directory);sys.path.insert(0,str(directory))
    global cfg
    cfg=runpy.run_path(str(directory/'cfg.py'))['cfg']
    cfg.addConn=cfg.addLongConn=cfg.addSubConn=cfg.addIClamp=0
    cfg.recordLFP=[];cfg.recordDipole=False;cfg.recordDipoles=False;cfg.analysis={};cfg.verbose=False;cfg.printRunTime=False
    parameters=runpy.run_path(str(directory/'netParams.py'))['netParams']
    parameters.popParams.clear();parameters.connParams.clear();parameters.subConnParams.clear()
    parameters.stimSourceParams.clear();parameters.stimTargetParams.clear()
    definitions={label:(kind,model) for label,kind,model in PROFILES.values()}
    for label in list(parameters.cellParams):
        if label not in definitions:del parameters.cellParams[label]
        else:parameters.cellParams[label]['conds']={'pop':label}
    for label,(kind,model) in sorted(definitions.items()):
        parameters.popParams[label]={'cellType':kind,'cellModel':model,'cellsList':[{'x':0.,'y':0.,'z':0.}]}
    sim.initialize(simConfig=cfg,netParams=parameters);sim.pc.nthread(1);sim.net.createPops();sim.net.createCells();h.celsius=34
    profiles={}
    for cell in sim.net.cells:
        label=cell.tags['pop'];soma=cell.secs['soma']['hObj'];h.distance(0,.5,sec=soma)
        apical=set(cell.secLists.get('apicdend',[]));points=[];geometry={}
        for name,info in cell.secs.items():
            sec=info['hObj']
            region=1 if name.startswith('soma') else 2 if name.startswith('axon') else 4 if name in apical or name.startswith(('apic','Adend')) else 3
            geometry[name]=[[float(sec.x3d(i)),float(sec.y3d(i)),float(sec.z3d(i)),float(sec.diam3d(i))] for i in range(sec.n3d())]
            for segment in sec:
                mechanisms={}
                for mechanism,specification in info.get('mechs',{}).items():
                    obj=getattr(segment,mechanism)
                    mechanisms[mechanism]={parameter:float(getattr(obj,parameter)) for parameter in specification}
                ions={}
                for ion,specification in info.get('ions',{}).items():
                    names={'e':'e'+ion,'i':ion+'i','o':ion+'o'}
                    ions[ion]={parameter:float(getattr(segment,names[parameter])) for parameter in specification if parameter in names}
                styles={ion:int(h.ion_style(ion+'_ion',sec=sec)) for ion in ions}
                points.append({'region':region,'reference_section':name,'distance_um':float(h.distance(segment.x,sec=sec)),
                               'Ra_ohm_cm':float(sec.Ra),'cm_uf_cm2':float(segment.cm),
                               'v_init_mv':float(info.get('vinit',cfg.hParams['v_init'])),
                               'mechanisms':mechanisms,'ions':ions,'ion_styles':styles})
        profiles[label]={'points':points,'reference_geometry':geometry}
    result={'source_commit':manifest['commit'],'neuron_version':neuron.__version__,'temperature_c':34,
            'method':'Published input parameters sampled by compartment class and soma path distance. These are transferred priors, not measured Minnie conductances.',
            'profiles':profiles}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({label:len(profile['points']) for label,profile in profiles.items()},indent=2),flush=True)


if __name__=='__main__':main()
