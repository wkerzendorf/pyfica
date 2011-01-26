#Algorithms for the geneticAlgorithm
import param
import config
import initialize
from numpy import random
import numpy as np
import dalekExceptions


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
    print "Creating %d random initial parameters..." % n
    randomParamSet = param.multiParam()
    randomParamList = []
    for i in range(n):
        randomParam = createRandomParam(randomParamFunc=randomParamFunc)
        randomParam.targetDir = 'gen%d' % i
        randomParamList.append(randomParam)
    randomParamSet.paramGrid = np.array(randomParamList)
    print "Finished creating random parameters"
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


def crossArith(parentA,parentB,mutationRate=0.2):
    childParam=param.param()
    rChoice=random.permutation(selElements)
    for element in rChoice:
        childParam[element]=np.mean([parentA[element],parentB[element]])
        
    if childParam['O']<0: raise Exception('Crossover: Child has negative oxygen')
    
    if any(np.array(childParam.comp.data.values())<0):
        raise dalekExceptions.geneticException('Negative values: %s' % childParam.comp.data.values())
    childParam['lum']=np.mean([parentA['lum'],parentB['lum']])
    childParam['vph']=np.mean([parentA['vph'],parentB['vph']])
    return childParam

def crossSingle(parentA,parentB):
    #Creating Child Param
    childParam=param.param()
    #Locking abundance ratios
    if config.GAConfDict['lockTiCr']:
        childParam.comp.lockTiCr=True
    if config.GAConfDict['lockScTi']:
        childParam.comp.lockScTi=True
    if config.GAConfDict['lockVCr']:
        childParam.comp.lockVCr=True
    
    #Initialising mutation counter
    childParam.mutations=0
    
    #selParameters=selElements+['lum','vph']
    
    selParameters=config.GAConfDict['selParameters']
    
    #closedGAParameters = list(set(config.GAConfDict['openGAParametersDefault'])
    #                          - set(config.GAConfDict['openGAParameters']))
    
    closedGAParameters = []
    selParamLen = len(selParameters)
    
    #Checking crossOverProbability
    if random.random()<config.GAConfDict['crossOverProbability']:
        crossOverPoint=random.uniform(1,selParamLen-1)
        paramChoice=random.permutation(selParameters)
        
        #Crossover Selecting from both parents
        for paramName in paramChoice[:crossOverPoint]:
            childParam[paramName]=parentA[paramName]
            
        for paramName in paramChoice[crossOverPoint:]:
            childParam[paramName]=parentB[paramName]
    else:
        #No crossover, selecting parameters from ParentA
        for paramName in selParameters:
            childParam[paramName]=parentA[paramName]
            
    if closedGAParameters != []:
        for paramName in selParameters:
            childParam[paramName] = parentA[paramName]
            
    childParam=mutateUniform(childParam)
    childParam.comp.resetOxygen()
    if childParam['O']<0: raise dalekExceptions.geneticException('Crossover: Child has negative oxygen')
    
    return childParam


def crossPGA(parentA,parentB,children=2):
    #Creating Child Param
    childParam=param.param()
    #Locking abundance ratios
    if config.GAConfDict['lockTiCr']:
        childParam.comp.lockTiCr=True
    if config.GAConfDict['lockScTi']:
        childParam.comp.lockScTi=True
    if config.GAConfDict['lockVCr']:
        childParam.comp.lockVCr=True
    #Initialising mutation counter
    childParam.mutations=0
    
    selParameters=selElements+['lum','vph']
    selParamLen=len(selParameters)
    
    #Checking crossOverProbability
    crossOverPoint=random.uniform(1,selParamLen-1)
    paramChoice=random.permutation(selParameters)
    #Crossover Selecting from both parents
    for paramName in paramChoice[:crossOverPoint]:
        childParam[paramName]=parentA[paramName]
        
    for paramName in paramChoice[crossOverPoint:]:
        childParam[paramName]=parentB[paramName]
            
    childParam=mutateUniform(childParam)
    childParam.comp.resetOxygen()
    if childParam['O']<0: raise dalekExceptions.geneticException('Crossover: Child has negative oxygen')
    return childParam


def mutateGauss(childParam,mutateScale=0.01,mutationRate=0.2):
    for element in selElements:
        mutationRate=mutationParams['elements']
        if random.random()<mutationRate:
            childParam[element]*=random.normal(1,mutationParams['elements'])
            childParam.mutations+=1
    if childParam['O']<0: raise dalekExceptions.geneticException('Crossover: Child has negative oxygen')
    if random.random()<mutationParams['lum'][0]:
        childParam['lum']*=random.normal(1,mutationParams['lum'][1])
        childParam.mutations+=1
    if random.random()<mutationParam['vph'][0]:
        childParam['vph']*=random.normal(1,mutationParams['vph'][1])
        childParam.mutations+=1
    return childParam


def mutateUniform(childParam):
    #HACK WARNING HACK WARNING
    if False:
    #if not lockExperiment:
        for element in selElements:
            #if element=='Ca': continue
            if random.random()<mutationParams['elements'][0]:
                childParam[element]*=random.uniform(1-config.GAConfDict['mutationParams']['elements'][1],
                                                    1+config.GAConfDict['mutationParams']['elements'][1])
                childParam.mutations+=1
    if childParam['O']<0: raise dalekExceptions.geneticException('Crossover: Child has negative oxygen')
    if random.random() < config.GAConfDict['mutationParams']['lum'][0]:
        #childParam['lum']*=random.uniform(1-mutationParams['lum'][1],
        #                                  1+mutationParams['lum'][1])
        childParam['lum']+=random.uniform(-config.GAConfDict['mutationParams']['lum'][1],
                                          +config.GAConfDict['mutationParams']['lum'][1])
        childParam.mutations+=1
    if random.random() < config.GAConfDict['mutationParams']['vph'][0]:
        #childParam['vph']*=random.uniform(1-mutationParams['vph'][1],
        #                                  1+mutationParams['vph'][1])
        childParam['vph']+=random.uniform(-config.GAConfDict['mutationParams']['vph'][1],
                                          +config.GAConfDict['mutationParams']['vph'][1])
        childParam.mutations+=1
    return childParam

def checkRatio(childParam):
    #checkedRatios=[]
    for element in config.GAConfDict['selElements']:
        if childParam[element]<1e-7:  return False
    if childParam['C']>0.125:
        #print "killed by C"
        return False
    if childParam['Si']<=0.01:
        #print "killed by Si"
        return False
    if childParam['Ca']>=0.05:
        #print "killed by Ca"
        return False
    if childParam['Ti']+childParam['Cr']>0.1:
        #print "killed by Ti Cr"
        return False
    
    if childParam['Ni0']>0.8:
        #print "killed by Ni0"
        return False
    
    CORatio=childParam['C']/childParam['O']
    if CORatio>=1:
        #print "killed by CO"
        return False
    
    if childParam['Mg']>=childParam['Si']:
        #print "killed by Mg"
        return False
    
    SiSRatio=childParam['Si']/childParam['S']
    if not (SiSRatio<=10 and SiSRatio>=1):
        #print "killed by SiS"
        return False
    
    Fe0Ni0Ratio=childParam['Fe0']/childParam['Ni0']
    if Fe0Ni0Ratio>=10:
        #print "killed by Fe0 Ni0"
        return False
    
    CrRatio=childParam['Ni0']*10
    if childParam['Cr']>CrRatio:
        #print "killed by Cr"
        return False
    
    
    return True

