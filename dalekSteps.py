#The fitting steps for the dalek fitter
initialLum=[8.8,10]
initialVph=[6000,15000]
if __name__!='__main__':
    from .launcher import multi_dist_launch
    from .read import readModels
import read
import fit
import os,pickle
import numpy as np
import util,config,initialize,abund,param,fitelem

##### FOR DEBUG PURPOSES ONLY SHOULD NOT BE USE LATER ON
remoteBasePath='/priv/manana1/wkerzend/sn_rad_trans'
localBasePath='/Users/wkerzend/wdata'
##### FOR DEBUG PURPOSES ONLY SHOULD NOT BE USE LATER ON


use_machines=['miner','myriad','maggot','merino','minotaur','miami','munch']
#use_machines=['miner']

def getInitLumParams(runDirs,params,doPickle=True):
    ##### FOR DEBUG PURPOSES ONLY SHOULD NOT BE USE LATER ON
    #params=pickle.load(file('params.pkl'))
    #runDirs=pickle.load(file('rundirs.pkl'))
    #runDirs=[item.replace(remoteBasePath,localBasePath) for item in runDirs]
    ##### FOR DEBUG PURPOSES ONLY SHOULD NOT BE USE LATER ON
    
    
    modelGrid=read.modelGrid(runDirs,params)
    if doPickle:
        pickle.dump(modelGrid,file('mg.pkl','w'))
def checkElementBounds(intervals,bounds):
    for element in intervals.keys():
        boundCheck=zip(intervals[element],bounds[element])
        checkedInterval=[np.max(boundCheck[0]),np.min(boundCheck[1])]
        if checkedInterval!=intervals[element]:
            print "%s interval out of Bounds. Changing from %s to %s"%(element,intervals[element],checkedInterval)
            intervals[element]=checkedInterval
    return intervals

def getNextElementParams(element,MG,initParam,bounds,sampleSize=3):
    elementAbundance=MG[element]
    curValue=initParam[element]
    lineMerits,eweights,dweights=fitelem.getSetMerits(element,MG)
    emask=np.all(eweights==0,axis=1)
    eweights[emask]=1
    dweights[emask]=1
    merits=np.average(lineMerits,weights=dweights*eweights,axis=1)
    merits[emask]=10.
    paramSortID=np.argsort(np.abs(merits))
    bestParam=paramSortID[:sampleSize]
    m,t=np.polyfit(elementAbundance[bestParam],merits[bestParam],1)
    bestValue=(1-t)/m
    newValue=np.mean([bestValue,curValue])
    newDev=np.abs(bestValue-curValue)
    newInterval=[newValue-newDev,newValue+newDev]
    boundCheck=zip(newInterval,bounds)
    checkedInterval=[np.max(boundCheck[0]),np.min(boundCheck[1])]
    steady=checkSteadiness(merits)
    print "Element %s steadiness %s fraction %s"%(element,steady,newDev/newValue)
    if steady>7 or newDev/newValue>0.1:
        finished=False
    else:
        finished=True
    return checkedInterval,finished

def checkSteadiness(merits):
    dMerits=np.diff(merits)
    steadParam=np.abs(np.sum(dMerits>0)-np.sum(dMerits<0))
    return steadParam
def getNextLumInterval(params,sampleSize=3):
    #lumMerits=read.getGridInt(params['divSpec'])
    lumMerits=2.*read.getGridInt(params['subSpec'])/read.getGridInt(params['addSpec'])
    paramSortID=np.argsort(np.abs(lumMerits))
    bestParam=paramSortID[:sampleSize]
    m,t=np.polyfit(params['lum'][bestParam],lumMerits[bestParam],1)
    bestLum=(-t)/m
    lumDev=min(np.abs(params['lum']-bestLum))
    interval=[bestLum-lumDev,bestLum+lumDev]
    minIntervalCheck=False
    
    steadParam=checkSteadiness(lumMerits)
    if steadParam<4:
        print "Warning merits not steady suggestions might be problematic (steadParam=%s)"%steadParam
    #WARNING. CHECK SAMPLES THIS WILL BE BAD LATER ON
    minInterval=0.02
    if 2*lumDev<minInterval:
        minIntervalCheck=True
        print "Interval too small"
        intervalIDs=[0,1]
        i=2
        while True:
            intervalPoints=params['lum'][paramSortID][intervalIDs]
            minIntervalPoint,maxIntervalPoint=min(intervalPoints),max(intervalPoints)
            if np.abs(minIntervalPoint-maxIntervalPoint)>minInterval:
                interval=[minIntervalPoint,maxIntervalPoint]
                break
            
            if len(intervalIDs)>len(lumMerits)-3:
                interval=[bestLum-2*lumDev,bestLum+2*lumDev]
                break
            
            intervalIDs.append(i)
            i+=1

    evalDic={'suggestValue':bestLum,
             'bestFitID':bestParam[0],
             'merit':abs(lumMerits[bestParam[0]]),
             'merits':lumMerits,
             'interval':interval,
             'dev':lumDev,
            'sortedModelIDX':paramSortID,
            'fitKey':'lum',
            'mininterval':minIntervalCheck,
            'steady':steadParam}
    return evalDic
    return bestLum,bestParam[0],abs(lumMerits[bestParam[0]]),[bestLum-lumDev,bestLum+lumDev]
    
def getNextVphInterval(params,sampleSize=4):
    contin=[item.smoothg(75) for item in params['subspec']]
    #contin=[item.fitContinuum(func='poly1') for item in params['subspec']]
    vphMerits=read.getGridSlope(params['subspec'])
    curVph=initialize.getCurVph()
    #vphMerits=vphMerits[(params['vph']>curVph-4000.)*(params['vph']<curVph+4000.)]
    paramSortID=np.argsort(np.abs(vphMerits))
    bestParam=paramSortID[:sampleSize]
    m,t=np.polyfit(params['vph'][bestParam],vphMerits[bestParam],1)
    bestVph=-t/m
    vphDev=min(np.abs(params['vph']-bestVph))
    minIntervalCheck=False
    minInterval=100
    interval=[bestVph-vphDev,bestVph+vphDev]
    steadParam=checkSteadiness(vphMerits)
    if 2*vphDev<minInterval:
        print "Interval too small"
        minIntervalCheck=True
        intervalIDs=[0,1]
        i=2
        while True:
            intervalPoints=params['vph'][paramSortID][intervalIDs]
            minIntervalPoint,maxIntervalPoint=min(intervalPoints),max(intervalPoints)
            if np.abs(minIntervalPoint-maxIntervalPoint)>minInterval:
                interval=[minIntervalPoint,maxIntervalPoint]
                break
            if len(intervalIDs)>len(vphMerits)-3:
                interval=[bestVph-2*vphDev,bestVph+2*vphDev]
                break
            intervalIDs.append(i)
            i+=1
    evalDic={'suggestValue':bestVph,
             'bestFitID':bestParam[0],
             'merit':abs(vphMerits[bestParam[0]]),
             'merits':vphMerits,
             'dev':vphDev,
             'interval':interval,
            'sortedModelIDX':paramSortID,
            'fitKey':'vph',
            'mininterval':minIntervalCheck,
            'steady':steadParam}
    return evalDic
    return bestVph,bestParam[0],abs(vphMerits[bestParam[0]]),[bestVph-vphDev,bestVph+vphDev]

def getNextIGEInterval(params,IGEElement,sampleSize=3):
    IGEMeritsOptical=read.getGridOptical(params['subspec'])
    IGEMeritsUV=read.getGridUV(params['subspec'])
    #VERY IMPORTANT NEED TO CHANGE WEIGHTS. THIS IS ONLY TESTING ATM
    IGEMerits=IGEMeritsUV-IGEMeritsOptical
    #IGEMerits=read.getGridUV(params['subSpec'])/read.getGridUV(params.origSpec)
    paramSortID=np.argsort(np.abs(IGEMerits))
    bestParam=paramSortID[:sampleSize]
    m,t=np.polyfit(params[IGEElement][bestParam],IGEMerits[bestParam],1)
    bestIGE=-t/m
    IGEDev=min(np.abs(params[IGEElement]-bestIGE))
    minInterval=1e-8
    interval=[bestIGE-IGEDev,bestIGE+IGEDev]
    steadParam=checkSteadiness(IGEMerits)
    minIntervalCheck=False
    if 2*IGEDev<minInterval:
        print "Interval too small"
        minIntervalCheck=True
        intervalIDs=[0,1]
        i=2
        while True:
            intervalPoints=params[IGEElement][intervalIDs]
            minIntervalPoint,maxIntervalPoint=min(intervalPoints),max(intervalPoints)
            if np.abs(minIntervalPoint-maxIntervalPoint)>minInterval:
                interval=[minIntervalPoint,maxIntervalPoint]
                break
            if len(intervalIDs)>len(IGEMerits)-3:
                interval=[bestIGE-2*IGEDev,bestIGE+2*IGEDev]
                break
            intervalIDs.append(i)
            i+=1
    evalDic={'suggestValue':bestIGE,
             'bestFitID':bestParam[0],
             'merit':abs(IGEMerits[bestParam[0]]),
             'merits':IGEMerits,
             'interval':interval,
             'dev':IGEDev,
            'sortedModelIDX':paramSortID,
            'fitKey':IGEElement,
            'mininterval':minIntervalCheck,
            'steady':steadParam}
    return evalDic
    return bestIGE,bestParam[0],abs(IGEMerits[bestParam[0]]),[bestIGE-IGEDev,bestIGE+IGEDev]
    
def readTriModel(params,runDirs,doPickle=True):
    lumModelGrid=read.modelGrid(runDirs[0],params[0])
    vphModelGrid=read.modelGrid(runDirs[1],params[1])
    IGEModelGrid=read.modelGrid(runDirs[2],params[2])
    if doPickle:
        pickle.dump(lumModelGrid,file('lmg.pkl','w'))
        pickle.dump(vphModelGrid,file('vmg.pkl','w'))
        pickle.dump(IGEModelGrid,file('img.pkl','w'))
    return lumModelGrid,vphModelGrid,IGEModelGrid
    
def getLumScale(params,noForFit=2):
    lums=params['lum']
    lumMerits=read.getGridInt(params['divSpec'])
    selLumMerits=np.argsort(lumMerits)[:noForFit]
    return np.polyfit(lums[selLumMerits],lumMerits[selLumMerits],1)
def getNextLumGuess(params,noForFit=2):
    lums=params['lum']
    lumMerits=read.getGridInt(params['subspec'])
    selLumMerits=np.argsort(lumMerits)[:noForFit]
    print lums[selLumMerits],lumMerits[selLumMerits]
    return np.polyfit(lums[selLumMerits],lumMerits[selLumMerits],1)
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

#def checkLu
