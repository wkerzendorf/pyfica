#implementing genetic algorithm for fica:
import numpy as np
import shutil
from numpy import random
#from ffnet import savenet, loadnet, exportnet
import elauncher
import param
import initialize
#import neuralDalek
import fit
import config
import time
import read
import cPickle
import pdb
import copy
import progressbar
import os
import plotSet
import dalekExceptions
import re
import numpy
import subprocess
import weakref
from glob import glob
import genFitness
#Constants


fitConf=config.getCurFitConfig()
lumLimits=fitConf['limit']['lum']
vphLimits=[8000,15000]
selElements=['C','Mg','Si','S','Ca','Cr','Ti','Ni0','Fe0']
selParameters=['lum','vph']+selElements
seed=25081106
#current modes are subpopulation or elitism
mode='elitism'
generationGapFraction=1.
subPopulationFraction=1.
#or
elitism=0.1

scaleFitness=True
Cmult=2.

mutationParams={'lum':[0.05,0.1],
                'vph':[0.2,1000],
                'elements':[0.07,0.1]}
crossOverProbability=0.9

#Lock abundances
lockTiCr = False
lockScTi = True
lockVCr = True

#open for mutation
#model

openGAParametersDefault = ['lum','vph']+selElements
openGAParameters = openGAParametersDefault


#locking certain values might help us figure out the fitting process
lockExperiment=False
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


def createRandomParam(randomFunc=createRandomLogNormalValueW7):
    #creating new paramObject
    randomParam=param.param()
    if lockTiCr:
        randomParam.comp.lockTiCr=True
    if lockScTi:
        randomParam.comp.lockScTi=True
    if lockVCr:
        randomParam.comp.lockVCr=True
    if randomFunc==createRandomLogNormalValueBestFit:
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
    
    
    rChoice=random.permutation(selElements)
    for element in rChoice:
        bounds=initialize.getElementBounds(element,randomParam.comp)
        curAbundance=randomParam[element]
        #newAbundance=random.uniform(bounds[0],bounds[1])
        newAbundance=randomFunc(curAbundance,element,bounds)
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
    
def createRandomParamSet(n):
    randomParamSet=param.multiParam()
    randomParamList=[]
    for i in range(n):
        randomParam=createRandomParam()
        randomParam.targetDir='gen%d'%i
        randomParamList.append(randomParam)
    randomParamSet.paramGrid=np.array(randomParamList)
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


def crossArith(parentA,parentB,mutationRate=0.2):
    childParam=param.param()
    rChoice=random.permutation(selElements)
    for element in rChoice:
        childParam[element]=np.mean([parentA[element],parentB[element]])
    if childParam['O']<0: raise Exception('Crossover: Child has negative oxygen')
    if any(np.array(childParam.comp.data.values())<0):
        raise dalekExceptions.geneticException('Negative values: %s'%childParam.comp.data.values())
    childParam['lum']=np.mean([parentA['lum'],parentB['lum']])
    childParam['vph']=np.mean([parentA['vph'],parentB['vph']])
    return childParam

def crossSingle(parentA,parentB):
    #Creating Child Param
    childParam=param.param()
    #Locking abundance ratios
    if lockTiCr:
        childParam.comp.lockTiCr=True
    if lockScTi:
        childParam.comp.lockScTi=True
    if lockVCr:
        childParam.comp.lockVCr=True
    #Initialising mutation counter
    childParam.mutations=0
    
    #selParameters=selElements+['lum','vph']
    
    selParameters=openGAParameters
    
    closedGAParameters = list(set(openGAParametersDefault) - set(openGAParameters))
    selParamLen=len(selParameters)
    
    #Checking crossOverProbability
    if random.random()<crossOverProbability:
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
    if lockTiCr:
        childParam.comp.lockTiCr=True
    if lockScTi:
        childParam.comp.lockScTi=True
    if lockVCr:
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
                childParam[element]*=random.uniform(1-mutationParams['elements'][1],
                                                    1+mutationParams['elements'][1])
                childParam.mutations+=1
    if childParam['O']<0: raise dalekExceptions.geneticException('Crossover: Child has negative oxygen')
    if random.random()<mutationParams['lum'][0]:
        #childParam['lum']*=random.uniform(1-mutationParams['lum'][1],
        #                                  1+mutationParams['lum'][1])
        childParam['lum']+=random.uniform(-mutationParams['lum'][1],
                                          +mutationParams['lum'][1])
        childParam.mutations+=1
    if random.random()<mutationParams['vph'][0]:
        #childParam['vph']*=random.uniform(1-mutationParams['vph'][1],
        #                                  1+mutationParams['vph'][1])
        childParam['vph']+=random.uniform(-mutationParams['vph'][1],
                                          +mutationParams['vph'][1])
        childParam.mutations+=1
    return childParam

def checkRatio(childParam):
    #checkedRatios=[]
    for element in selElements:
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

def breed(randomModelSet,popNum=None,select=selectRoulette,cross=crossSingle):

    if popNum==None: popNum=len(randomModelSet.grid)
    
    fitness = randomModelSet['fitness']
    w=randomModelSet['w']
    k=0
    loopNo=0
    m=0
    mutated=0
    
    populationList=[]
    
    if scaleFitness:
        fitness=genFitness.fitnessScale(fitness,Cmult)
    
    k=1
    #breeding until enough new children
    while k<popNum:
        loopNo+=1
        if loopNo>100*popNum: pdb.set_trace()
        parentAID=select(fitness)
        parentBID=select(fitness)
        
        try:
            child=cross(randomModelSet.grid[parentAID].param,randomModelSet.grid[parentBID].param)
        except dalekExceptions.geneticException:
            continue
        if not checkRatio(child):
            #print "Child died due to ratio problems\n%s"%child.comp.data
            m+=1
            continue
        child.targetDir='gen%d'%k
        populationList.append(child)
        k+=1
    population=param.multiParam()
    population.paramGrid=np.array(populationList)
    print "%s children died because of abundance ratios"%m
    return population
        

def launchRandomParamSet(randomParamSet):
    launcherSteps.getModel()
    
def readFitness(searchPattern='generation???.dat'):
    i=[]
    fitness=[]
    fitmax=[]
    for fname in np.sort(glob(searchPattern)):
        i.append(int(re.search('\d\d\d',fname).group()))
        fitness.append(np.mean(np.loadtxt(fname)))
        fitmax.append(np.max(np.loadtxt(fname)))
    return i,fitness,fitmax


def evolve2(generations=200,populationSize=150,savePath='.',continueGA=False,mutationRate=0.2,mutationScale=0.05):

#Initializing random seed
    random.seed(seed)

    evolveHeader='#gen specno fitness %s\n'%' '.join(['lum','vph']+selElements)
    evolveHeaderFMT='%d %d %s '+' '.join(['%s']*len(['lum','vph']+selElements))
    evolveDBPath=os.path.join(savePath,'evolution.dat')
    
    file(evolveDBPath,'w').write(evolveHeader)
    
    generationGapNo=int(populationSize*generationGapFraction)
    
    subPopulationNo=int(populationSize*generationGapFraction*subPopulationFraction)
    
    gws=elauncher.gateways()
    
    basePath=config.getAutoDir()
    
    #Smoothing UV
    origSpec=config.getOrigSpec(preProcess=True)
    try:
            os.remove('break_after_loop')
    except:
        pass
    try:
        os.remove('debug_after_loop')
    except:
        pass
    #continue GA if stuff there if not make new one
    
    if continueGA:
        x=[]
        y=[]
        y2=[]
        yerr=[]
        for fname in np.sort(glob(os.path.join(savePath,'generation*.pkl'))):
            x.append(int(re.search('\d+',os.path.basename(fname)).group()))
            fitness=np.loadtxt(fname.replace('.pkl','.dat'))
            y.append(np.mean(fitness))
            y2.append(np.max(fitness))
            yerr.append(np.std(fitness))
        plotSet.genPlotGenVSFitness(x,y,y2,yerr,outName='gen_vs_fitness_log.png')
        plotSet.genPlotGenVSFitness(x,y,y2,yerr,logPlot=False,outName='gen_vs_fitness_linear.png')
        curGenerationModel=cPickle.load(file(os.path.join(savePath,'generation%04d.pkl'%np.max(x))))
        firstGeneration=np.max(x)+1
    else:
        x=[]#generation Number
        y=[]#fitness medium
        y2=[]#fitness max
        yerr=[]#fitness std
        curGenerationSet=createRandomParamSet(populationSize)
        firstGeneration=0
        gws.checkAvailability()
        print gws.availability
        curGenerationModel,tmp=elauncher.cloudLaunch(curGenerationSet.paramGrid,gws.getAvailGateWays(),origSpec=origSpec)        
        keepChildren=param.multiParam()
        keepChildren.grid=np.array([])
        
    for i in range(firstGeneration,generations):
        
        curTime=time.time()
        #Skipping for continueGA
        if not continueGA:
            curGenerationModel.grid=np.concatenate((keepChildren.grid,curGenerationModel.grid))
            curGenerationModel._initSpecs()
            if i%5 == 0:
                fname="generation%04d.pkl"%i
                print "Saving current generation to %s"%fname
                cPickle.dump(curGenerationModel,file(os.path.join(savePath,fname),'w'))
                print "Saved"
        else:
            continueGA=False
        fitness=curGenerationModel['fitness']
        fitnessIDX=np.argsort(fitness)[::-1]
        #TEST CASE FOR VARYING OPEN PARAMETERS ________ TEST _____ TEST
        if i == 1:
            openGAParameters=['lum','vph']
        
        #TEST CASE FOR VARYING OPEN PARAMETERS ________ TEST _____ TEST
        #Saving the fitness, making plots
        np.savetxt(os.path.join(savePath,'generation%04d.dat'%i),fitness)
        x.append(i)
        y.append(np.mean(fitness))
        y2.append(np.max(fitness))
        yerr.append(np.std(fitness))
        
        #Saving to evolution DB
        evolveDB=zip([i]*populationSize,range(populationSize),
            *[curGenerationModel[item] for item in ['fitness','lum','vph']+selElements])
        np.savetxt(file(evolveDBPath,'a'),evolveDB,fmt=evolveHeaderFMT)
        
        #Plotting the generation
        plotSet.genPlotSingleGeneration(curGenerationModel,
                                        fitness,
                                        generation=i,
                                        no=10,
                                        pdfName=os.path.join(savePath,'generation%04d.pdf'%i))
        
        plotSet.genPlotStatus(curGenerationModel,
                              pdfName=os.path.join(savePath,'gen_status%04d.pdf'%i))
        
        plotSet.genPlotHistogram(curGenerationModel,
                                 fitness=fitness,
                                 pdfName=os.path.join(savePath,'gen_histogram%04d.pdf'%i))
        
        plotSet.genPlotGenVSFitness(x,y,y2,yerr,
                                    outName=os.path.join(savePath,'gen_vs_fitness_log.png'))
        plotSet.genPlotGenVSFitness(x,y,y2,yerr,
                                    logPlot=False,
                                    outName=os.path.join(savePath,'gen_vs_fitness_linear.png'))
        #Checking for break conditions
        if os.path.exists('break_after_loop') or os.path.exists(os.path.join(savePath,'break_after_loop')): break
        if os.path.exists('debug_after_loop'):
            try:
                os.remove('debug_after_loop')
            except:
                pass
            pdb.set_trace()
        
        #made plots
        
        #start selective breeding
        #First the children that are kept for the subpopulation are put into a seperate variable
        if mode=='subpopulation':
            if populationSize==subPopulationNo:
                keepChildren=param.multiParam()
                keepChildren.grid=np.array([])
            else:
                keepChildren=read.modelGrid(paramList=curGenerationModel.
                             grid[fitnessIDX[:(populationSize-subPopulationNo)]])
                keepChildren._initSpecs()
        elif mode=='elitism':
            keepChildren=read.modelGrid(paramList=curGenerationModel.
                             grid[fitnessIDX[:int(populationSize*elitism)]])
            keepChildren._initSpecs()
        #Now we get the population that is used for breeding and submit it to the breed function
        breedPopulation=read.modelGrid(paramList=curGenerationModel.
                                       grid[fitnessIDX[:generationGapNo]])
        
        if mode=='subpopulation':
            curGenerationSet=breed(breedPopulation,
                  popNum=subPopulationNo)
        elif mode=='elitism':
            curGenerationSet=breed(breedPopulation,
                  popNum=int(populationSize*(1-elitism)))
        del curGenerationModel
        
        #Time is kept
        
        ficaTime=time.time()
    
        #Network check for available nodes
        gws.checkAvailability()
        print gws.availability
        
        #Calculating the new generation with elauncher.
        curGenerationModel,tmp=elauncher.cloudLaunch(curGenerationSet.paramGrid,gws.getAvailGateWays(),origSpec=origSpec)
        
        
        #Printing time statements
        print "Took %s for the fica runs"%(time.time()-ficaTime)
        print "Took %s seconds for last loop"%(time.time()-curTime)
        
def getGenModelGrid(nodeDict,params,basePath,origSpec):
    js=launcher2.getOptimalJS(len(params),nodeDict)
    imports=('pyfica.read','pyfica.launcher2','os','shutil','subprocess','numpy')
    workerTemplate=pp.Template(js,genRunModel,(fitInvSSE,fitFunc),imports)
    workers=[]
    for param in params:
        workers.append(workerTemplate.submit(param,origSpec,basePath))
    modelGrid=getWorkerResult(workers)
#    js.destroy()
    return modelGrid
def getWorkerResult(workers):
    pbar=progressbar.ProgressBar(maxval=len(workers)).start()
    while True:
        currentStatus=[item.finished for item in workers]
        if all(currentStatus): break
        pbar.update(np.sum(currentStatus))
        time.sleep(5)
    return [item() for item in workers]
def genRunModel(param,origSpec,basePath):
    currentNode=pyfica.launcher2.getNodeName()
    fica_bin=pyfica.launcher2.machineConf[currentNode]['fica_bin']
    ficaWorkDir=os.path.join(basePath,param.targetDir)
    try:
        os.makedirs(ficaWorkDir)
    except:
        shutil.rmtree(ficaWorkDir)
        os.makedirs(ficaWorkDir)
    param.write2file(ficaWorkDir)
    #print "Running %s on %s in dir %s"%(fica_bin,currentNode,ficaWorkDir)
    proc=subprocess.Popen([fica_bin],stdout=-1,cwd=ficaWorkDir,close_fds=True)
    #Saving Log File:
    file(os.path.join(ficaWorkDir,'fica.log'),'w').write(proc.stdout.read())
    proc.stdout.close()
    #print "tested %s"%proc.stdout.read()
    #return "I have run %s on %s"%(basePath,getNodeName())
    model=pyfica.read.model(param,ficaWorkDir,origSpec)
    model.fitness=fitFunc(model)
    return model
def prepGATry(gaNo,description):
    
    wikiRoot='/home/wkerzend/public_html/dalek/data/pages'
    wikiMediaRoot='/home/wkerzend/public_html/dalek/data/media'
    gaStore='/priv/miner3/skymap/wkerzend/ga_store2'
    SNName=config.getName().lower()
    t=config.getTimeFromExplosion()
    GeneralOverviewFile=os.path.join(wikiRoot,'sn_gas.txt')
    SNOverviewFile=os.path.join(wikiRoot,'%s_gas.txt'%SNName)
    SNTimeFile=os.path.join(wikiRoot,'%s_t%.2f.txt'%(SNName,t))
    GATryFile=os.path.join(wikiRoot,'%s_t%.2f_ga%d.txt'%(SNName,t,gaNo))
    mediaNamespace=":%s:%s:%s:"%(SNName,'%s_t%.2f'%(SNName,t),'ga_try%d'%gaNo)
    try:
        os.makedirs(os.path.join(gaStore,SNName,'%s_t%.2f'%(SNName,t),'ga_try%d'%gaNo))
    except:
        pass
    os.system('ln -s %s .'%os.path.join(gaStore,SNName,'%s_t%.2f'%(SNName,t),'ga_try%d'%gaNo))
    paramString="""
lumLimits=%s

vphLimits=%s

selElements=%s

selParameters=%s

seed=%s

mode=%s

generationGapFraction=%.3f

subPopulationFraction=%.3f

elitism=%.3f

scaleFitness=%s

Cmult=%s

mutationParams=%s

crossOverProbability=%s
    """%(lumLimits,vphLimits,selElements,selParameters,seed,mode,generationGapFraction,subPopulationFraction,elitism,scaleFitness,Cmult,mutationParams,crossOverProbability)
    
    templatePage="""{{%sgen_vs_fitness_linear.png?400}}
{{%sgen_vs_fitness_log.png?400}}

x=generation
y=mean(fitness) in this case sum square error
yerr=std(fitness)

===Current Parameter===
%s

===Description===
%s


^ Spectra      ^ Status       ^ Histograms          ^
| {{filelist>%sgener*.pdf&style=table&tableheader=1&tableshowdate=1&tableshowsize=1&order=desc}}   | {{filelist>%sgen_st*.pdf&style=table&tableheader=1&tableshowdate=1&tableshowsize=1&order=desc}}    | {{filelist>%sgen_his*.pdf&style=table&tableheader=1&tableshowdate=1&tableshowsize=1&order=desc}}        |

    """%(mediaNamespace,mediaNamespace,paramString,description,mediaNamespace,mediaNamespace,mediaNamespace)
    #Create page structure and create media structure
    if os.path.exists(SNOverviewFile):
        if os.path.exists(SNTimeFile):
            file(SNTimeFile,'a').write("[[%s_t%.2f_ga%d|%s t=%.2f ga try=%d]]\n\n"%(SNName,t,gaNo,SNName,t,gaNo))
            file(GATryFile,'w').write(templatePage)
        else:
            file(SNOverviewFile,'a').write('[[%s_t%.2f|%s t=%.3f]]\n\n'%(SNName,t,SNName,t))
            file(SNTimeFile,'a').write("[[%s_t%.2f_ga%d|%s t=%.2f ga try=%d]]\n\n"%(SNName,t,gaNo,SNName,t,gaNo))
            file(GATryFile,'w').write(templatePage)
        
    else:
        file(GeneralOverviewFile,'a').write("[[%s_gas|%s GA overview]]\n\n"%(SNName,SNName))
        file(SNOverviewFile,'a').write('[[%s_t%.2f|%s t=%.3f]]\n\n'%(SNName,t,SNName,t))
        file(SNTimeFile,'a').write("[[%s_t%.2f_ga%d|%s t=%.2f ga try=%d]]\n\n"%(SNName,t,gaNo,SNName,t,gaNo))
        file(GATryFile,'w').write(templatePage)
        
def makeEvolutionFile():
    parameters=['lum','vph']+selElements
    data=[]
    header='# gen mean_fit std_fit max_fit '+' '.join(['mean_%s'%param for param in parameters])+" "+' '.join(['max_%s'%param for param in parameters])
    for fname in np.sort(glob("*.pkl")):
        curGeneration=int(re.search("\d+",fname))
        modelGrid=cPickle.load(file(fname))
        fitness=loadtxt(fname.replace('.pkl','.dat'))
        fitMaxID=np.argmax(fitness)
        tmpRow=[curGeneration,meanFit,stdFit,maxFit]+[np.mean(modelGrid[param]) for param in parameters]+[modelGrid[param][fitMaxID] for param in parameters]
        data.append(tmpRow)