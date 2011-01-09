#The fitting steps for the dalek fitter
initialLum=[8.8,10]
initialVph=[6000,15000]
from .launcher import multi_dist_launch
from .read import readModels
import os,pickle
import numpy as np
import util,config,initialize,abund


use_machines=['miner','myriad','maggot','merino','minotaur','miami','munch']
#use_machines=['miner']
def initLumVphGrid(step=1,gridSize=10,doPickle=True):
    t=config.getTimeFromExplosion()
    vph=initialize.time2vph(t)
    mainConfigDir=config.getMainConfigDir()
    
    modelW7=initialize.readW7Data(os.path.join(mainConfigDir,'w7.combined.dat'))
    initDica=config.getVanillaDica()
    initDica['t']=t
    #Preparing the normalization
    initComp=config.getVanillaComp()
    initComp.update(initialize.getW7Comp(modelW7,t))
    initComp=abund.setNiDecay(initComp,t)
    initComp=abund.setCONe(initComp)
    initComp=abund.normAbundances(initComp)
    initRunDir='init_lumvph_run/'
    initialVph=[0.6*vph,1.4*vph]
    return lumVphGrid(initialLum,initialVph,1,initRunDir,initComp=initComp,initDica=initDica,gridSize=gridSize,doPickle=doPickle)
def singleElement(element,abundances,step,runDir,initDica=None,initComp=None):
    comps=[]
    dicas=[{}]*len(abundances)
    #Creating the dicas and comps
    for abund in abundances:
        comps.append({element:abund})
    #Launching with different properties
    multi_dist_launch(dicas,comps,use_machines,runDir,init_dica=initDica,init_comp=initComp)
    return readModels(runDir)
    
def lumVphGrid(lumLimits,vphLimits,step,runDir,gridSize=10,initDica={},initComp={},doPickle=True):
    lums,vphs=util.makeGrid(lumLimits,vphLimits,gridSize)
    print "Doing grid: Lum: %s vph: %s"%(lumLimits,vphLimits)
    dicas=[]
    for lum,vph in zip(lums,vphs):
        tmpDica=initDica.copy()
        tmpDica.update({'log_lbol':lum,'v_ph':vph})
        #print lum,vph
        #print tmpDica
        dicas.append(tmpDica.copy())
    comps=[initComp]*len(dicas)
    print "Using Machines %s"%use_machines
    #Launching with different properties
    multi_dist_launch(dicas,comps,use_machines,runDir,init_dica=initDica,init_comp=initComp)
    #Reading model
    if doPickle:
        pickle.dump(readModels(runDir),file('tmp.pkl','w'))
    else:
        return readModels(runDir)