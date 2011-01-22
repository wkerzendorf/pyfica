#Algorithms for the geneticAlgorithm
import param
import config
import initialize
from numpy import random
import numpy as np


def createRandomLogNormalValueW7(value,element,bounds):
    if value>bounds[1]: value=bounds[1]-np.abs(bounds[1]-value)
    elif value<bounds[1]: value=bounds[0]+np.abs(bounds[0]-value)
    if random.random()>0.5:
        while True:
            factor=10**abs(random.normal(loc=0,scale=-np.log10(value/bounds[1])/2))
            newUpperValue=factor*value
            if newUpperValue<bounds[1]: return newUpperValue
    else:
        while True:    
            factor=10**-abs(random.normal(loc=0,scale=-np.log10(bounds[0]/value)/2))
            newLowerValue=factor*value
            if newLowerValue>bounds[0]: return newLowerValue
            
            
def createRandomLogNormalValueBestFit(value,element,bounds):
    fname=glob('*.bf.pkl')[0]
    bestFitParam=cPickle.load(file(fname))
    value=bestFitParam[element]
    if not lockExperiment:
        if value>bounds[1] or value<bounds[0]: value=np.mean(bounds)
        if random.random()>0.5:
            while True:    
                factor=10**abs(random.normal(loc=0,scale=-np.log10(value/bounds[1])/2))
                newUpperValue=factor*value
                if newUpperValue<bounds[1]: return newUpperValue
        else:
            while True:    
                factor=10**-abs(random.normal(loc=0,scale=-np.log10(bounds[0]/value)/2))
                newLowerValue=factor*value
                if newLowerValue>bounds[0]: return newLowerValue
    else:
        return value


def createRandomParam(randomParamFunc=None):
    #creating new paramObject
    lumLimits = config.GAConfDict['lumLimits']
    vphLimits = config.GAConfDict['vphLimits']
    randomParam=param.param()
    if config.GAConfDict['lockTiCr']:
        randomParam.comp.lockTiCr=True
    if config.GAConfDict['lockScTi']:
        randomParam.comp.lockScTi=True
    if config.GAConfDict['lockVCr']:
        randomParam.comp.lockVCr=True
    if randomParamFunc==createRandomLogNormalValueBestFit:
        fname=glob('*.bf.pkl')[0]
        bestFitParam=cPickle.load(file(fname))
        #Limiting luminosity and photospheric velocity near best fit value, commented out atm
        #randomParam['lum']=random.uniform(bestFitParam['lum']-0.1,bestFitParam['lum']+0.1)
        #no limit for lum parameter
        #randomParam['lum']=bestFitParam['lum']
        randomParam['lum']=random.uniform(*lumLimits)
    else:
        randomParam['lum']=random.uniform(*lumLimits)
    
    randomParam['vph']=random.uniform(*vphLimits)
    
    
    rChoice=random.permutation(config.GAConfDict['selElements'])
    for element in rChoice:
        bounds=initialize.getElementBounds(element,randomParam.comp)
        curAbundance=randomParam[element]
        #newAbundance=random.uniform(bounds[0],bounds[1])
        newAbundance=randomParamFunc(curAbundance,element,bounds)
        if any(np.array(bounds)<0): raise Exception('Negative Bounds')
        randomParam[element]=newAbundance
        randomParam.comp.resetOxygen()
        if randomParam['O']<0: pdb.set_trace()
        if randomParam[element]<0: pdb.set_trace()
        if randomParam.comp.data.has_key('element'): raise Exception()
    #randomParam['Ca']=0.01
    #randomParam['vph']=11700
    return randomParam
    #Selecting random
    
def createRandomParamSet(n, randomParamFunc=createRandomLogNormalValueW7):
    randomParamSet = param.multiParam()
    randomParamList = []
    for i in range(n):
        randomParam = createRandomParam(randomParamFunc=randomParamFunc)
        randomParam.targetDir = 'gen%d' % i
        randomParamList.append(randomParam)
    randomParamSet.paramGrid = np.array(randomParamList)
    return randomParamSet

def selectRoulette(fitness):
    #normFitness=(np.max(fitness)*1.1-fitness)
    normFitness=fitness/np.sum(fitness)
    randNum = random.random()
    partSum=0
    for i in range(len(normFitness)):
        partSum += normFitness[i]
        if partSum >= randNum:
            return i
    else:
        return i
