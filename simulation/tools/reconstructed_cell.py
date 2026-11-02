"""Recorded SWC geometry with explicitly transferred electrical-density priors."""

from collections import Counter, defaultdict
import itertools

import numpy as np


class Curve:
    def __init__(self, distances, values):
        grouped=defaultdict(list)
        for distance,value in zip(distances,values):grouped[round(float(distance),6)].append(float(value))
        self.x=np.array(sorted(grouped));self.y=np.array([np.median(grouped[x]) for x in self.x])
        if not np.isfinite(self.y).all():raise ValueError('Nonfinite electrical reference parameter')
        self.constant=bool(np.all(self.y==self.y[0]))

    def at(self, distances):
        return np.interp(distances,self.x,self.y)


class ElectricalProfile:
    """Interpolate parameter densities by compartment and physical path distance."""
    def __init__(self, reference):
        self.regions={}
        for region in sorted({point['region'] for point in reference['points']}):
            points=[p for p in reference['points'] if p['region']==region]
            combinations={tuple(sorted(p['mechanisms'])) for p in points}
            if len(combinations)!=1:raise ValueError('Reference mechanism complements vary within a compartment class')
            flat=[]
            for point in points:
                row={('Ra',):point['Ra_ohm_cm'],('cm',):point['cm_uf_cm2'],('vinit',):point['v_init_mv']}
                for mechanism,parameters in point['mechanisms'].items():
                    for parameter,value in parameters.items():row[('mech',mechanism,parameter)]=value
                for ion,parameters in point['ions'].items():
                    for parameter,value in parameters.items():row[('ion',ion,parameter)]=value
                flat.append(row)
            keys=set().union(*(row.keys() for row in flat));fields={}
            for key in keys:
                present=[i for i,row in enumerate(flat) if key in row]
                fields[key]=Curve([points[i]['distance_um'] for i in present],[flat[i][key] for i in present])
            styles={}
            for ion in set().union(*(point['ion_styles'].keys() for point in points)):
                values={point['ion_styles'][ion] for point in points if ion in point['ion_styles']}
                if len(values)!=1:raise ValueError('Reference ion handling varies within a compartment class')
                styles[ion]=values.pop()
            self.regions[region]={'fields':fields,'mechanisms':set(points[0]['mechanisms']),'ion_styles':styles}

    def region(self, requested):
        if requested in self.regions:return self.regions[requested],requested
        fallback=3 if 3 in self.regions else 1
        return self.regions[fallback],fallback


class ReconstructedCell:
    _instances=itertools.count()

    def __init__(self,path,node,gid,profile,parallel_context,d_lambda=.1,frequency=100.):
        from neuron import h
        self.h=h;self.root_id=node['root_id'];self.gid=gid;self.instance=next(self._instances);self.stims=[]
        self.tags={'root_id':self.root_id,'pop':node['profile'],'cellModel':node['cell_model'],'cellType':node['cell_type'],
                   'manual_cell_type':node['manual_cell_type'],'ei_class':node['ei_class']}
        reader=h.Import3d_SWC_read();reader.input(str(path));h.Import3d_GUI(reader,False).instantiate(self)
        if not getattr(self,'soma',[]):raise ValueError(f'No reconstructed soma: {self.root_id}')
        self.secs={};self.section_regions={};self.section_areas={};self.field_fallbacks=Counter()
        for region,prefix in [(1,'soma'),(2,'axon'),(3,'dend'),(4,'apic')]:
            for i,section in enumerate(list(getattr(self,prefix,[]))):
                name='soma' if prefix=='soma' and i==0 else f'{prefix}_{i}'
                self.secs[name]={'hObj':section};self.section_regions[name]=region
        sections=[item['hObj'] for item in self.secs.values()]
        if len(sections)!=len(list(self.all)):raise RuntimeError('Imported sections lack compartment labels')
        h.distance(0,.5,sec=self.secs['soma']['hObj'])
        self.secLists={'perisom':['soma'],'spiny':[],'alldend':[],'apicdend':[],'basal':[],'axonal':[]}
        for name,info in self.secs.items():
            section=info['hObj'];region=self.section_regions[name];reference,used=profile.region(region)
            if used!=region:self.field_fallbacks[f'{region}->{used}']+=1
            if section.L<=0:raise ValueError('Nonpositive imported section length')
            midpoint=float(h.distance(.5,sec=section));fields=reference['fields']
            section.Ra=float(fields[('Ra',)].at(midpoint));section.cm=float(fields[('cm',)].at(midpoint))
            wavelength=float(h.lambda_f(frequency,sec=section))
            if not np.isfinite(wavelength) or wavelength<=0:raise ValueError('Invalid electrotonic length')
            section.nseg=int((section.L/(d_lambda*wavelength)+.9)/2)*2+1
            info['mechs']={mechanism:{} for mechanism in reference['mechanisms']}
            for mechanism in reference['mechanisms']:section.insert(mechanism)
            segments=list(section);distances=np.array([h.distance(segment.x,sec=section) for segment in segments])
            for key,curve in fields.items():
                if key[0]=='Ra':continue
                values=curve.at(distances)
                if key[0]=='vinit':info['vinit']=float(values[0]) if curve.constant else values;continue
                if key[0]=='cm':attribute='cm'
                elif key[0]=='mech':attribute=f'{key[2]}_{key[1]}'
                elif key[0]=='ion':attribute=({'e':'e'+key[1],'i':key[1]+'i','o':key[1]+'o'})[key[2]]
                else:raise RuntimeError('Unknown reference field')
                if curve.constant:setattr(section,attribute,float(values[0]))
                else:
                    for segment,value in zip(segments,values):setattr(segment,attribute,float(value))
            for ion,style in reference['ion_styles'].items():
                h.ion_style(ion+'_ion',style%4,(style//8)%4,(style//32)%2,(style//64)%2,(style//4)%2,sec=section)
            self.section_areas[name]=sum(float(h.area(segment.x,sec=section)) for segment in segments)
            if region in {3,4}:
                self.secLists['spiny'].append(name);self.secLists['alldend'].append(name)
                self.secLists['apicdend' if region==4 else 'basal'].append(name)
            if region==2:self.secLists['axonal'].append(name)
        parallel_context.set_gid2node(gid,0)
        self.spike_detector=h.NetCon(self.secs['soma']['hObj'](.5)._ref_v,None,sec=self.secs['soma']['hObj'])
        self.spike_detector.threshold=0.
        parallel_context.cell(gid,self.spike_detector)

    def initV(self):
        for info in self.secs.values():
            values=info['vinit'];section=info['hObj']
            if np.isscalar(values):section.v=float(values)
            else:
                for segment,value in zip(section,values):segment.v=float(value)

    def __str__(self):
        return f'minnie_skeleton_{self.root_id}_{self.instance}'
