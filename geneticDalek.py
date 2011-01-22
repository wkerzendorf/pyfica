#implementing genetic algorithm for fica:
import numpy as np
import shutil
from numpy import random
import datetime

import elauncher
import param
import initialize
import genFitness
import geneticDalekAlgo
import dalekDB
import inspect
import config
#Constants




openGAParametersDefault = ['lum','vph']+selElements
openGAParameters = openGAParametersDefault


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
        





def evolve(conn, SNSpecID, GARunID=None, description=None, generations=200, populationSize=150,
           breedFunc = breed, select=selectRoulette, cross=crossSingle,
           randomParamFunc = geneticDalekAlgo.createRandomLogNormalValueW7):
    
    
    
    
    #Run the GA

    startTime = datetime.datetime.now()
    #Initializing cursor
    curs = conn.cursor()
    
    
    #trying to remove existing break and debug switches
    try:
        os.remove('break_after_loop')
    except:
        pass
    try:
        os.remove('debug_after_loop')
    except:
        pass
    
    #Initializing random seed
    random.seed(GAConfDict['seed'])
    
    #Initializing mode dependent constants:
    generationGapNo=int(populationSize*generationGapFraction)
    
    subPopulationNo=int(populationSize*generationGapFraction*subPopulationFraction)
    
    #Launching the gateways
    gws=elauncher.gateways()
    
    #getting origSpec and preparing it
    rawOrigSpec = curs.execute('select SPECTRUM from SN_SPECTRA where ID=%d' % SNSpecID).fetchall()[0]
    origSpec = initialize.preProcessOrigSpec(rawOrigSpec)
    
    
    breedSource = dalekDB.makeZipPickle(inspect.getsource(breedFunc))
    GAConfSource = dalekDB.makeZipPickle(GAConfDict)
    crossSource = dalekDB.makeZipPickle(inspect.getsource(crossFunc))
    selectSource = dalekDB.makeZipPickle(inspect.getsource(selectFunc))
    
    #Getting time
    t = config.getTimeFromExplosion(conn, SNSpecID)
    #setting time
    param.snT = t
    
    #Checking for continue or new
    if GARunID!=None:
        raise NotImplementedError('soon.....')
        
    
    #What IS GAID???????
    
    if GARunID == None: #or GA_CUR_GEN=0    
        #creating new GA_RUN entry and inserting the code of various functions
        
        curs.execute('insert into GA_RUN(DESCRIPTION, SN_ID, START_TIME,'
                     'SN_SPECTRUM, GA_POP_SIZE, GA_CONF_DICT, GA_BREED_FUNC,'
                     'GA_CROSS_FUNC, GA_SELECT_FUNC, GA_FITNESS_FUNC)'
                     ' VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                     (description, SNSpecID, startTime, 
                      origSpec, populationSize, GAConfSource, breedSource,
                      crossSource, selectSource, None))
        
        GARunID = curs.lastrowid
        
        curs.execute('insert GA_GENERATIONS(GA_RUN_ID) VALUES(?)', GARunID)
        
        generationID = curs.lastrowid
        
        curGenerationSet = geneticDalekAlgo.createRandomParamSet(populationSize)
        
        #Inserting the dica into the database, this dica will be used throughout the run
        dicaID = dalekDB.insertDica(conn, curGenerationSet.paramGrid[0].dica)
        
        firstGeneration = 0
        
        # Checking availability on gateways
        gws.checkAvailability()
        print gws.availability
        
        #Launching 
        curGenerationModel = elauncher.cloudLaunch(curGenerationSet.paramGrid, gws.getAvailGateWays(), origSpec=origSpec)
        
        keepChildren = param.multiParam()
        keepChildren.grid = np.array([])
        keepModelIDs = np.array([])
        keepFitness = np.array([])
        
    
    
    for i in range(firstGeneration,generations):
    
        #getting current time
        curTime=time.time()
        
        #getting fitnesses from the model
        fitness=curGenerationModel['fitness']
        fitnessIDX=np.argsort(fitness)[::-1]
        
        #saving model to db
        modelIDs = curGenerationModel.toDB(conn, GARunID=GARunID, dicaID=dicaID, storeLList=True, storeWParam=True)
        
        
        #main link between models and the GA for the keeping
        dalekDB.insertGAIndividual(conn, GARunID, keepModelIDs, keepFitness)
        
        
        #main link between models and the GA 
        dalekDB.insertGAIndividual(conn, GARunID, modelIDs, fitness)

        
        
        curGenerationModel.grid=np.concatenate((keepChildren.grid,curGenerationModel.grid))
        curGenerationModel._initSpecs()
        
        
        
        
        #Checking for break conditions
        if os.path.exists('break_after_loop') or os.path.exists(os.path.join(savePath,'break_after_loop')): break
        if os.path.exists('debug_after_loop'):
            try:
                os.remove('debug_after_loop')
            except:
                pass
            pdb.set_trace()
        
        
        #start selective breeding
        #First the children that are kept for the subpopulation are put into a seperate variable
        if mode == config.GAConfDict['subpopulation']:
            if populationSize==subPopulationNo:
                keepChildren=param.multiParam()
                keepChildren.grid=np.array([])
            else:
                keepChildren = model.modelGrid(paramList=curGenerationModel.
                                grid[fitnessIDX[:(populationSize-subPopulationNo)]])
                keepModelIDs = modelIDs[fitnessIDX[:(populationSize - subPopulationNo)]]
                keepFitness = fitness[fitnessIDX[:(populationSize - subPopulationNo)]]
                keepChildren._initSpecs()
        
        elif mode == config.GAConfDict['elitism']:
            keepChildren = modelmodelGrid(paramList=curGenerationModel.
                            grid[fitnessIDX[:int(populationSize*elitism)]])
            keepModelIDs = modelIDs[fitnessIDX[:(populationSize*elitism)]]
            keepFitness = fitness[fitnessIDX[:(populationSize*elitism)]]
            keepChildren._initSpecs()
        
        #Now we get the population that is used for breeding and submit it to the breed function
        breedPopulation=model.modelGrid(paramList=curGenerationModel.
                                       grid[fitnessIDX[:generationGapNo]])
        
        if mode == config.GAConfDict['subpopulation']:
            curGenerationSet=breed(breedPopulation,
                  popNum=subPopulationNo, select=selectFunc, cross=crossFunc)
        elif mode == config.GAConfDict['elitism']:
            curGenerationSet=breed(breedPopulation,
                  popNum=int(populationSize*(1-elitism)), select=selectFunc, cross=crossFunc)
        
        del curGenerationModel
        
        #Time is kept
        ficaTime=time.time()
    
        #Network check for available nodes
        gws.checkAvailability()
        print gws.availability
        
        #Calculating the new generation with elauncher
        curGenerationModel = elauncher.cloudLaunch(curGenerationSet.paramGrid, gws.getAvailGateWays(), origSpec=origSpec)
        
        #Printing time statements
        print "Took %s for the fica runs"%(time.time()-ficaTime)
        print "Took %s seconds for last loop"%(time.time()-curTime)


