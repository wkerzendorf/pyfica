#old version of dalek before update with spectral class. Basically the first version of Dalek

machine_conf={'mithrandir-local':
    {'processors':2,
    'type':'local', 
    'path':'/Volumes/manana1/sn_rad_trans/', 
    'fica_script':'~/scripts/python/pyfica/bin/fica.exe'}, 
    'myriad':
        {'processors':8, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':'/usr/bin/rsh myriad %slaunch_fica_solarisx86.sh myriad'%default_script_path, 
        'prep_proc':\
        '/usr/bin/rsh myriad %sprep_dir_myriad.sh' %default_script_path}, 
    'moron-local':
        {'processors':2, 
        'type':'local', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '%slaunch_fica_moron-local.sh'%default_script_path }, 
    'maggot':
        {'processors':4, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh maggot %slaunch_fica_linuxx86_64.sh maggot'%default_script_path}, 
    'miner':
        {'processors':32,
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh miner %slaunch_fica_linuxx86_64.sh miner'%default_script_path}, 

    'moron':
        {'processors':4, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh moron %slaunch_fica_linuxx86_64.sh moron'%default_script_path}, 
    'miami':
        {'processors':2, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh miami %slaunch_fica_solarisx86.sh miami'%default_script_path}, 
    'merino':
        {'processors':2, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh merino %slaunch_fica_solarisx86.sh merino'%default_script_path}, 
    'minotaur':
        {'processors':2, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh minotaur %slaunch_fica_solarisx86.sh minotaur'%default_script_path},
    'munch':
        {'processors':2, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh munch %slaunch_fica_solarisx86.sh munch'%default_script_path},
    'macerate':
        {'processors':2, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh macerate %slaunch_fica_solarisx86_32.sh macerate'%default_script_path}, 
    'magpie':
        {'processors':2, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh magpie %slaunch_fica_solarissparc.sh magpie'%default_script_path},
    'manana':
        {'processors':2, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh manana %slaunch_fica_solarissparc.sh manana'%default_script_path},
    'marvin':
        {'processors':2, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh marvin %slaunch_fica_solarissparc.sh marvin'%default_script_path},
    'matrix':
        {'processors':2, 
        'type':'remote', 
        'path':'/priv/manana1/wkerzend/sn_rad_trans/', 
        'fica_script':\
        '/usr/bin/rsh matrix %slaunch_fica_solarissparc.sh matrix'%default_script_path}
    }
from scipy.integrate import trapz
from scipy import polyfit, ndimage
#from scipy.interpolate import splrep, splev
#BAD THINGS CHECK YOUR DATA FOR THESE THINGS:
# it must be evenly spaced , otherwise rewrite si fitting routine and maybe others

#Bad solution for temperature fitting. if it doesnt have 5739 SI III line in model its disregarded as bad. need to change this
from launcher import launch, multi_dist_launch

from .fitelem import siMerit, siSection, getCaMerit, getSMerit, getMgMerit
from .launcher import launch, multi_dist_launch
from .util import findLast, dictBuilder, calcTemp, median_spectra, gaussian_spectra, multi_fit, getBinMerits, syn_div_orig
from scipy import interpolate, polyfit, polyval
from bisect import bisect
from numpy import exp, power, array, matrix, mean, std, argmax, argmin, argsort, abs, arange, median, average,  vstack, linspace, diff
from math import pi
import math
from pyspec.spectrum import spectrum
from pyphot import datafile
import sys
from pickle import dump, load
from random import random
from glob import glob
from os import path
import os
import pickle
from .fileio import spctfile, dicafile, compfile, sbibfile
import shutil
import pygsl.statistics
initial_lum=(10.0, 9.0)
initial_vph=(12000.0, 10000.0)
funcDict={'xisquared':lambda data: mean((data[:,1]-1)**2), 
                'mean':lambda data: mean(data[:,1]-1)}



bin_weights=[1, 1, 2, 2]
seedno=250819801106

def calc_integral(spec, start=None, end=None):
    if not (start==None and end==None):
        lower=bisect(spec[:, 0], start)
        higher=bisect(spec[:, 0], end)
        spec=spec[lower:higher, :]

    #Integrating the spectra
    return trapz(spec[:, 1], spec[:, 0])

def mod2dist(mod):
    return 3.08568025e18*10**(0.2*mod+1)

def calc_spec_grad(spec, sample=10):
    grad_spec=[]
    for i in range(len(spec)-sample):
        j=i+sample
        grad=(spec[j, 1]-spec[i, 1])/(spec[j, 0]-spec[i, 0])
        grad_spec.append([spec[i, 0], grad])
    return grad_spec



def bbody(wl, T):
    wl=wl*1e-8
    h=6.6262e-27
    c=2.9979e10
    k=1.3807e-16
    #print "Am I over 600?:", (h*c)/(wl*k*T)
    i=((8*pi*h*c**2)/(wl**5))*(1/(exp((h*c)/(wl*k*T))-1))
    return i/1e8


def calc_grid_lum_vph(lum_bound, vph_bound, no_cells=10):
    lums=[]
    vphs=[]
    lums=linspace(lum_bound[0], lum_bound[1], no_cells)
    vphs=linspace(vph_bound[0], vph_bound[1], no_cells)
    #for i in range(no_cells):
     #   lums=linspace.append(lum_bound[0]+(lum_bound[1]-lum_bound[0])*i/no_cells)
     #   vphs.append(vph_bound[0]+(vph_bound[1]-vph_bound[0])*i/no_cells)
    dicas=[]
    comps=[]
    for lum in lums:
        for vph in vphs:
            dicas.append({'log_lbol':lum, 'v_ph':vph})
            comps.append({})
    return dicas, comps
    #multi_launch(dicas, comps, 'myriad', run_dir=run_dir)

def calc_jacobian(lum, vph, grad, intercept):
    print "lum %s vph %s grad %s intercept %s"%(lum, vph, grad, intercept)

    lumderiv=2*((grad[1]-grad[0])/(lum[1]-lum[0])+
                                              (intercept[1]-1)*(intercept[1]-intercept[0])/(lum[1]-lum[0]))
    vphderiv=2*((grad[1]-grad[0])/(vph[1]-vph[0])+
                                              (intercept[1]-1)*(intercept[1]-intercept[0])/(vph[1]-vph[0]))
    return lumderiv, vphderiv
def mkdivspect(slower=0, upper=20000):
    snspect=snspect[bisect(snspect[:, 0], lower):bisect(snspect[:, 0], upper)]
    return syn_div_orig(aspect, snspect)
def slope(data):
    return polyfit(data[:, 0], data[:, 1], 1)[0]

def offset(data):
    return mean(data[:, 1])
def advisor():
    print "test"
def referee(accused, evidence, penalty_score):
    ids=range(len(penalty_score))
    ids.pop(accused)
    before, after=evidence
    #before=before[ids]
    #after=after[ids]
    if mean(after) > mean(before):
        penalty_score[accused]=penalty_score[accused]+1
    else:
        penalty_score[accused]=1
    #penalty system implementation
    return penalty_score



def newton_lum_vph(init_dica, init_comp, init_cond=None):
    merit = lambda fits: (fits[0]**2+(fits[1]-1)**2)
    snspect=array(datafile('snspect/origspect.dat').read_data(sel_columns=[0,1]))
    lums=[]
    vphs=[]
    merits=[]
    penalty_score=[1]*len(freq_bins)
    if init_cond==None:
        slopes=[]
        offsets=[]
        for i in range(2):
            #initializing values
            old_slopes=slopes
            old_offsets=offsets
            dica=init_dica
            dica_param={'log_lbol':initial_lum[i], 'v_ph':initial_vph[i]}            
            #launching fica
            aspect=launch(dica_param,{})
            #fitting functions and eval
            divspect=syn_div_orig(aspect, snspect)

            slopes, offsets=multi_fit(divspect, freq_bins, [slope, offset])
            fits_merits=[merit(ifit) for ifit in zip(slopes, offsets)]
            lums.append(initial_lum[i])
            vphs.append(initial_vph[i])
            merits.append(fits_merits)
            print "input values:", zip(slopes, offsets)
            print "merits:", fits_merits
            max_id=argmax(fits_merits)
            print "max_id:", max_id
        js=[]
        deltaxs=[]
        for i in range(len(slopes)):
            ij=array(calc_jacobian(initial_lum, initial_vph, (old_slopes[i], slopes[i]), (old_offsets[i],offsets[i])))
            ideltax=-merit([slopes[i], offsets[i]])/ij
            js.append(ij)
            deltaxs.append(ideltax)

        deltax=mean(deltaxs, axis=0)
        #init_cond=[initial_lum, initial_vph, (oldgrad, grad), (oldinter, inter)]
        #dump(init_cond, file('init_cond.pkl', 'w'))
        #print init_cond
    else:
        lums, vphs, grads, inters=load(file(init_cond))
        j=array(calc_jacobian(initial_lum, initial_vph, (oldgrad, grad), (oldinter, inter)))
    i=0

    head_banging_lum=0
    head_banging_vph=0
    while i<15:

        old_lum=dica_param['log_lbol']
        old_vph=dica_param['v_ph']

        print "new proposed luminosities %s and photospheric velocities:\n %s"%(array(deltaxs)[:, 0]+old_lum, array(deltaxs)[:, 1]+old_vph)
        new_lum=dica_param['log_lbol']+deltax[0]
        new_vph=dica_param['v_ph']+deltax[1]
        if old_lum==new_lum or old_vph==new_vph or str(new_vph)=='nan' or str(new_vph)=='nan':
            print "Ran into trouble !!! Exiting ....."
            break
        #print "1/j %s and deltax %s"%(1/j, deltax)
        print "Current luminosity %s current photospheric velocity %s"%(dica_param['log_lbol'], dica_param['v_ph'])
        #print "determined grad %s and intercept%s"%(slopes[max_id], offsets[max_id])
        #print "new luminosity: %s new photospheric velocity: %s"%(new_lum, new_vph)
        #if raw_input("Want to continue [Y/n]").lower().startswith('n'):
        i=i+1    
        print "----------------------------------------------------------------------------------------"
        #    break
        #keeping fica from crashing, I need to address this issue some time:

        if new_lum<8.9: new_lum=8.9+0.2*random()
        if new_lum>10: new_lum = 10-0.2*random()

        if new_vph<6000: new_vph=6000+1000*random()
        if new_vph>15000: new_vph = 15000-1000*random()

        dica_param['log_lbol']=new_lum
        dica_param['v_ph']=new_vph
        aspect=launch(dica_param,{})
        old_slopes=slopes
        old_offsets=offsets
        divspect=syn_div_orig(aspect, snspect)
        slopes, offsets=multi_fit(divspect, freq_bins, [slope, offset])
        old_fits_merits=fits_merits
        fits_merits=[merit(ifit) for ifit in zip(slopes, offsets)]
        #print "input values:", zip(slopes, offsets)
        #print "merits:", fits_merits

        js=[]
        deltaxs=[]
        for i in range(len(slopes)):
            ij=array(calc_jacobian((new_lum, old_lum),(new_vph, old_vph), (old_slopes[i], slopes[i]), (old_offsets[i],offsets[i])))
            ideltax=-merit([slopes[i], offsets[i]])/ij
            js.append(ij)
            deltaxs.append(ideltax)
        deltax=deltaxs[max_id]
        penalty_score=referee(max_id, [old_fits_merits, fits_merits], penalty_score)
        print penalty_score
        max_id=argmax(array(fits_merits)/array(penalty_score))
        print "max_id:", max_id
        lums.append(new_lum)
        vphs.append(new_vph)
        merits.append(fits_merits)
    return lums, vphs, merits
def checkElemLines( divspect, elem, llist, sampleSize=0, selWL=None):    
    sel_list=[item for item in llist if item[3].lower()==elem.lower() and min(divspect[:, 0])<item[1]<max(divspect[:, 0])]
    result=[]
    for line in sel_list:
        if selWL!=None and line[2] not in selWL: continue
        if not divspect[0, 0]<line[1]<divspect[-1, 0]: continue
        else:
            if sampleSize==0:
                lindex=bisect(divspect[:, 0], line[1])
                result.append(divspect[lindex, 1])
            else:
                lindex=bisect(divspect[:, 0], line[1])
                result.append(mean(divspect[lindex-sampleSize:lindex+sampleSize, 1]))
    return result

def getSelLlist(divspect, llist):
    ndivspect=abs(divspect[:, 1]-1)
    max_spect=argmax(ndivspect)
    if divspect[max_spect, 1]<1: sign=-1
    else: sign=1
    print "Found Maximum at %s Angstrom"%divspect[max_spect, 0]
    i=0
    while True:
        i+=1
        if ndivspect[i+max_spect]<ndivspect[max_spect]/2: break
    lower=divspect[max_spect, 0]
    upper=divspect[max_spect+2*i, 0]
    print "Lower %s Upper %s"%(lower, upper)
    sel_list=[item for item in llist if lower<item[1]<upper]
    if sel_list==[]: raise Exception('Found no lines in the region, giving up')
    return sel_list, [lower, upper], sign

def calc_grid_elem(element,bounds, samples=10):
    comps=[]
    dicas=[]

    for abund in linspace(bounds[0], bounds[1], samples):
        comps.append({element:abund})
        dicas.append({})
    return dicas, comps

def getLineMerit(models, elem, selWl=None, offsets=None, func='mean', returnMean=True, sample=10):
    abund=[]
    merits=[]
    for i, model in enumerate(models):
        if selWl!=None:
            
            elemWl=[[item[1]-sample, item[1]+sample] 
                                          for item in model['sbib']['llist'] if item[3].lower()==elem.lower() and int(item[2]) in selWl]
        else:
            elemWl=[[item[1]-sample, item[1]+sample] for item in model['sbib']['llist'] if item[3].lower()==elem.lower()]
        if offsets:
            for offset in offsets:
                elemWl+=[[low+offset, high+offset] for low, high in elemWl]
        if returnMean:
            merits.append(mean(multi_fit(model['divspect'], elemWl,[funcDict[func]])))
        else:
            merits.append(multi_fit(model['divspect'], elemWl,[funcDict[func]]))
        abund.append(model['comp'][elem])
    return abund, merits

def getLineMeritMin(abund, merits, sampling=1000, thresh=0.1, retMerit=False):
    delta=(abund[-1]-abund[0])/sampling
    abund_new=arange(abund[0], abund[-1]+delta, delta)
    rep=interpolate.splrep(abund, merits)
    merits_new=interpolate.splev(abund_new, rep)
    #Differentiating
    diffMerits=[(merits_new[i+1]-merits_new[i])/(abund_new[i+1]-abund_new[i]) for i in range(len(abund_new)-1)]
    diffMerits.append(diffMerits[-1])
    if retMerit:
        return diffMerits
    else:
        for i in range(len(diffMerits)):
            if abs(diffMerits[i])<thresh: return abund_new[i]



def singleFit(steps=['lumvph1', 'ca1'],  initModel=None, dalekStep=None, initStep=1):
    selDica=None
    selModel=None
    initSample=10
    origspect=array(datafile('origspect.dat').read_data())
    origspect[:, 0]/=(1+dicafile('init_cond/dica.dat').read_data()['z'])
    if not os.path.exists('dalek/'): os.mkdir('dalek')
    use_machines=['miner','myriad','maggot','moron','merino','minotaur','miami','munch', 'macerate']

    def initDicaComp():
        if initModel!=None:
            if initModel=='init':
                return pickle.load(file('init_cond/init_model.pkl'))
            return initModel
        elif dalekStep!=None:
            return pickle.load(file('dalek/model.pkl.step'+dalekStep))
        else:
            step=findLast('dica.dat.step', 'dalek/')
            return readModel(step, origspect, 'dalek/', '.step')
    def updateModel(model, origspect, bounds):
        print
    def writeModel2Dalek(model, step, runDir):
        if type(step)==int: step="%03d"%step

        dicafile('dalek/dica.dat.step'+step, 'w').write_data(model['dica'])
        compfile('dalek/comp.ind.step'+step, 'w').write_data(model['comp'])
        datafile('dalek/divspect.dat.step'+step, 'w').write_data(model['divspect'])
        shutil.copy(runDir+'sbib.dat.bak%s'%model['id'], 'dalek/sbib.dat.step'+step)
        shutil.copy(runDir+'spct.dat.bak%s'%model['id'], 'dalek/spct.dat.step'+step)
        pickle.dump(model, file('dalek/model.pkl.step'+step, 'w'))
    
    def fitS(selModel, step, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='srun%s/'%step
        dicas, comps=calc_grid_elem('S', [selComp['Si']/5,selComp['Si']/3])
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        model=readModels(origspect, run_dir, backup=True)
        #getting the limits to the Ca line by differentiating 2x
        abund, merits=getSMerit(model)
        wAbund=pygsl.statistics.wmean(abs(merits[:,0]), abund)
        wAbundSd=pygsl.statistics.wsd(abs(merits[:,0]), abund)
        selModelNo=argmin(abs(merits))
        selModel.update(model[selModelNo])
        print "Selected Model %s with S Abundance %s"%(selModelNo, selModel['comp']['S'])
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    def fitMg(selModel, step, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='mgrun%s/'%step
        dicas, comps=calc_grid_elem('Mg', [0, .3])
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        model=readModels(origspect, run_dir, backup=True)
        #getting the limits to the Ca line by differentiating 2x
        abund, merits=getMgMerit(model)
        wAbund=pygsl.statistics.wmean(abs(merits[:,0]), abund)
        wAbundSd=pygsl.statistics.wsd(abs(merits[:,0]), abund)
        selModelNo=argmin(abs(merits))
        selModel.update(model[selModelNo])
        print "Selected Model %s with Mg Abundance %s"%(selModelNo, selModel['comp']['Mg'])
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    def fitCa(selModel, step, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='carun%s/'%step
        caTester=mean(checkElemLines(selDivspect, 'ca', selSbib['llist']))
        if caTester>1: 
            print "Adding Ca"
            dicas, comps=calc_grid_elem('Ca', [selComp['Ca'], 0.1])
        else: 
            print "Reducing Ca"        
            dicas, comps=calc_grid_elem('Ca', [0, selComp['Ca']])
        multi_dist_launch(dicas,comps,use_machines,run_dir, init_dica=selDica, init_comp=selComp)
        
        model=readModels(origspect, run_dir, backup=True)
        #getting the limits to the Ca line by differentiating 2x
        abund, merits=getCaMerit(model)
        
        wAbund=pygsl.statistics.wmean(abs(merits[:,0]), abund)
        wAbundSd=pygsl.statistics.wsd(abs(merits[:,0]), abund)
        selModelNo=argmin(abs(merits[:,0]))
        
        selModel.update(model[selModelNo])
        print "Selected Model %s with Ca Abundance %s"%(selModelNo, selModel['comp']['Ca'])
        print "Calculated weighted mean and weighted standard deviation for abundance: mean %s sd %s"%(wAbund, wAbundSd)
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    def fitSi(selModel, step, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='sirun1/'
        if siMerit([selModel])[1][0]>0: 
            print "Adding Si"
            if selComp['Si']+selComp['O']>0.6: upperLimit=0.6
            else: upperLimit=selComp['Si']+selComp['O']
            dicas, comps=calc_grid_elem('Si', [selComp['Si'], upperLimit])
        else: 
            print "Reducing Si"        
            dicas, comps=calc_grid_elem('Si', [0, selComp['Si']])
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        #Launching the fit
        model=readModels(origspect, run_dir, backup=True)
        #abund=[item['comp']['Si'] for item in model]
        abund, merits=siMerit(model)
        wAbund=pygsl.statistics.wmean(abs(merits[:, 0]), abund)
        wAbundSd=pygsl.statistics.wsd(abs(merits[:, 0]), abund)
        
        selModelNo=argmin(abs(merits))
        selModel.update(model[selModelNo])
        print "Selected Model %s"%selModelNo
        print "Si abundance %s"%selModel['comp']['Si']
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
        
    #Grid Fitters
    def initLumVphFit(step=1):
        init_lum=[8.8,10]
        init_vph=[6000,15000]
        init_rundir='init_lumvph_run/'
        run_dir=init_rundir
        dicas,comps=calc_grid_lum_vph(init_lum,init_vph,initSample)
        print "Doing first grid: Lum: %s vph: %s"%(init_lum, init_vph)
        multi_dist_launch(dicas,comps,use_machines,run_dir)
        for i in range(2, 5):
            model=readModels(origspect, run_dir, backup='True')
            lums,vphs=zip(*[[item['dica']['log_lbol'],item['dica']['v_ph']] for item in model])
            merits=getBinMerits(model)
            avgmerit=average(merits, axis=1)
            lumavg, lumstd=pygsl.statistics.wmean(1/avgmerit**2, lums), pygsl.statistics.wsd(1/avgmerit**2,lums)
            vphavg, vphstd=pygsl.statistics.wmean(1/avgmerit**2, vphs), pygsl.statistics.wsd(1/avgmerit**2,vphs)
            lumBounds=[lumavg-lumstd, lumavg+lumstd]
            vphBounds=[vphavg-vphstd, vphavg+vphstd]
            dicas,comps=calc_grid_lum_vph(lumBounds,vphBounds,initSample)
            print "Doing %s grid: Lum: %s vph: %s"%(i,lumBounds, vphBounds)
            multi_dist_launch(dicas,comps,use_machines,run_dir)
        model=readModels(origspect, run_dir, backup='True')
        merits=getBinMerits(model)
        avgmerit=average(merits, axis=1)
        selModelNo=argmin(avgmerit)
        print "First selection: %s"%selModelNo
        writeModel2Dalek(model[selModelNo], step, run_dir)
        return model[selModelNo]
    def fitGridLumVph(selModel, step=1, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        lumBounds, vphBounds=selModel['params']['lumvph']['bounds']
        run_dir='lumvph_run%s/'%step
        print "Doing grid: Lum: %s vph: %s"%(lumBounds, vphBounds)
        dicas,comps=calc_grid_lum_vph(lumBounds,vphBounds,initSample)
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        model=readModels(origspect, run_dir, backup='True')
        lums,vphs=zip(*[[item['dica']['log_lbol'],item['dica']['v_ph']] for item in model])
        if params['fitMode']=='normGrid':
            merits=getBinMerits(model)
            avgmerit=average(merits, axis=1)
        elif params['fitMode']=='tempGrid':
            #tempBins=[[4380, 4480], [5500, 5600]]
            abund, merits=getLineMerit(model,'Si',selWl=[5739],offsets=[100],sample=50,returnMean=False)
            merits=array([item[0] if item!=[[]] else [10,20] for item in merits])
            #Checking if both lines have the same ratio as in he origspect
            tempTest1=abs(merits[:, 0]-merits[:, 1])
            #Checking if both lines are closeto origspect
            tempTest2=abs(average(merits, axis=1))
            avgmerit=tempTest1+tempTest2
        lumavg, lumstd=pygsl.statistics.wmean(1/avgmerit**2, lums), pygsl.statistics.wsd(1/avgmerit**2,lums)
        vphavg, vphstd=pygsl.statistics.wmean(1/avgmerit**2, vphs), pygsl.statistics.wsd(1/avgmerit**2,vphs)
        lumBounds=[lumavg-lumstd, lumavg+lumstd]
        vphBounds=[vphavg-vphstd, vphavg+vphstd]
        if params['fitMode']=='tempGrid':
            t=model[0]['dica']['t']
            tempMax=calcTemp(t, lumBounds[1], vphBounds[0])
            tempMin=calcTemp(t, lumBounds[0], vphBounds[1])
            print "Temperature min %s and max %s"%(tempMin, tempMax)
        selModelNo=argmin(abs(avgmerit))
        selModel.update(model[selModelNo])
        print "Selected model %s"%selModelNo
        print "Parameters for selected Model: lum %s vph %s"%(selModel['dica']['log_lbol'], selModel['dica']['v_ph'])
        print "New Bounds are: lum %s and vph %s"%(lumBounds, vphBounds)
        print "Calculated lum %s and vph %s. (this is just the weighted mean and will not be used further)"%(lumavg, vphavg)
        
        if params['saveBounds']:
            selModel['params']['lumvph']['bounds']=[lumBounds, vphBounds]
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    
    def fitGridMetals(selModel, step, params=None):
        if params==None:
            fit='log'
        else:
            fit=params['fitMode']
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='metalrun%s/'%step
        comps=[]
        dicas=[]
        TiCrBounds=selModel['params']['TiCr']['bounds']
        NiBounds=selModel['params']['Ni']['bounds']
        if fit=='log':
            TiCrRange=10**linspace(math.log(TiCrBounds[0], 10),math.log(TiCrBounds[1], 10), 5)
        elif fit=='linear':TiCrRange=linspace(TiCrBounds[0],TiCrBounds[1], 5)
        
        for TiCrAbund in TiCrRange:
            for NiAbund in linspace(0.01, 0.1, 5):
                for FeAbund in linspace(0.001,0.05, 5):
                    comps.append({'Ti':TiCrAbund, 'Cr':TiCrAbund, 'Ni':NiAbund, 'Fe':FeAbund})
                    dicas.append({})
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        model=readModels(origspect, run_dir, backup=True)
        meritUV=getBinMerits(model,[[0, 3500]],func='mean',doProc=None)
        meritAll=getBinMerits(model,[[0,10000]],func='mean',doProc=None)
        CoAbund, CoMerit=getLineMerit(model,'Co')
        selModelNo=argmin(abs(CoMerit)+abs(meritUV[:, 0])+abs(meritAll[:, 0]))
        print "Selected Model %s"%selModelNo
        print "Abundances for selected Model: TI Cr %s Ni %s Fe %s"%\
            (selModel['comp']['Ti'], selModel['comp']['Ni'], selModel['comp']['Fe'])
        selModel.update(model[selModelNo])
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    
   
    #MAINLOOP
    fitterDict={'lumvph1':[initLumVphFit, {'saveBounds':True}],
        'ca1':[fitCa, None],
        'si1':[fitSi, None],
        's':[fitS, None], 
        'temp':[fitGridLumVph, {'saveBounds':True,'fitMode':'tempGrid'}],
        'metal1':[fitGridMetals, None],
        'glumvph':[fitGridLumVph, {'saveBounds':True, 'fitMode':'normGrid'}],
        'glumvphns':[fitGridLumVph, {'saveBounds':False,'fitMode':'normGrid'}], 
        'mg':[fitMg, None]}
    for i, step in enumerate(steps):
        if step=='lumvph1':
            fitterDict[step]()
            continue
        else:
            selModel=fitterDict[step][0](selModel, i+initStep, params=fitterDict[step][1])


 
