#implementing genetic algorithm for fica:
import numpy as np
import shutil
from numpy import random
import datetime
import time
import os
import elauncher
import param
import initialize
import genFitness
import geneticDalekAlgo
import dalekDB
import inspect
import config
import model
import dalekExceptions
#Constants
import pprint

pp = pprint.PrettyPrinter(indent=4)



openGAParametersDefault = ['lum','vph']+config.GAConfDict['selElements']
openGAParameters = openGAParametersDefault



def breed(randomModelSet,popNum=None,select=geneticDalekAlgo.selectRoulette,
          cross=geneticDalekAlgo.crossSingle):

    if popNum==None: popNum=len(randomModelSet.grid)
    
    fitness = randomModelSet['fitness']
    w=randomModelSet['w']
    k=0
    loopNo=0
    m=0
    mutated=0
    
    populationList=[]
    
    if config.GAConfDict['scaleFitness']:
        fitness=genFitness.fitnessScale(fitness,config.GAConfDict['Cmult'])
    
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
        if not geneticDalekAlgo.checkRatio(child):
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
        





def evolve(conn, SNSpecID, GARunID=None, description=None, generations=20, populationSize=150,
           breedFunc = breed, selectFunc=geneticDalekAlgo.selectRoulette,
           crossFunc=geneticDalekAlgo.crossSingle,
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
    random.seed(config.GAConfDict['seed'])
    
    #Initializing mode dependent constants:
    generationGapNo=int(populationSize * config.GAConfDict['generationGapFraction'])
    
    subPopulationNo=int(populationSize * config.GAConfDict['generationGapFraction']
                        * config.GAConfDict['subPopulationFraction'])
    
    #Launching the gateways
    gws=elauncher.gateways()
    
    
    
    #getting origSpec and preparing it
    rawOrigSpec = curs.execute('select SPECTRUM from SN_SPECTRA where ID=%d' % SNSpecID).fetchall()[0]
    origSpec = initialize.preProcessOrigSpec(rawOrigSpec[0])
    
    
    breedSource = dalekDB.makeZipPickle(inspect.getsource(breedFunc))
    GAConfSource = dalekDB.makeZipPickle(config.GAConfDict)
    crossSource = dalekDB.makeZipPickle(inspect.getsource(crossFunc))
    selectSource = dalekDB.makeZipPickle(inspect.getsource(selectFunc))
    fitSource = dalekDB.makeZipPickle(inspect.getsource(genFitness.fitFunc))
    
    
    #Getting SNConfigDict
    SNConfigDict = config.getSNConfigDict(conn)
    SNConfigDict['t'] = config.getTimeFromExplosion(conn, SNSpecID, SNConfigDict)
    
    
    
    #setting time
    param.SNConfigDict = SNConfigDict
    
    #Checking for continue or new
    if GARunID!=None:
        raise NotImplementedError('soon.....')
        
    
    
    if GARunID == None: #or GA_CUR_GEN=0    
        #creating new GA_RUN entry and inserting the code of various functions
        
        curs.execute('insert into GA_RUN(DESCRIPTION, SN_ID, START_TIME,'
                     'SN_SPECTRUM, GA_POP_SIZE, GA_GENERATION_SIZE, GA_CONF_DICT, GA_BREED_FUNC,'
                     'GA_CROSS_FUNC, GA_SELECT_FUNC, GA_FITNESS_FUNC)'
                     ' VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                     (description, SNSpecID, startTime, 
                      origSpec, populationSize, generations, GAConfSource, breedSource,
                      crossSource, selectSource, fitSource))
        
        GARunID = curs.lastrowid
        
        curs.execute('insert into GA_GENERATION(GA_RUN_ID) VALUES(?)', (GARunID,))
        
        generationID = curs.lastrowid
        
        curGenerationSet = geneticDalekAlgo.createRandomParamSet(populationSize)
        
        #Inserting the dica into the database, this dica will be used throughout the run
        dicaID = dalekDB.insertDica(conn, curGenerationSet.paramGrid[0].dica)
        
        firstGeneration = 0
        
        # Checking availability on gateways
        gws.checkAvailability()
        pp.pprint(gws.availability)
        
        #Launching 
        curGenerationModel = elauncher.cloudLaunch(curGenerationSet.paramGrid, gws.getAvailGateWays(), origSpec=origSpec)
        #getting fitnesses from the model
        fitness = curGenerationModel['fitness']
        fitnessIDX = np.argsort(fitness)[::-1]
        
        keepChildren = param.multiParam()
        keepChildren.grid = np.array([])
        keepModelIDs = np.array([])
        keepFitness = np.array([])
        
    
    
    for i in range(firstGeneration,generations):
    
        #getting current time
        curTime=time.time()
        
        
        
        #saving model to db
        modelIDs = curGenerationModel.toDB(conn, dicaID=dicaID, storeLList=False, storeWParam=False)
        
        
        #main link between models and the GA for the keeping
        dalekDB.insertGAIndividual(conn, generationID, keepModelIDs, keepFitness)
        
        #main link between models and the GA 
        dalekDB.insertGAIndividual(conn, generationID, modelIDs, fitness)
        
        curs.execute('update GA_RUN set GA_CUR_GEN=? where ID=?', (generationID, GARunID))
        
        #getting new generation
        curs.execute('insert into GA_GENERATION(GA_RUN_ID) VALUES(?)', (GARunID,))
        
        generationID = curs.lastrowid

        
        #uniting old keepChildren and current generation model
        curGenerationModel.grid=np.concatenate((keepChildren.grid,curGenerationModel.grid))
        curGenerationModel._initSpecs()
        modelIDs = np.concatenate((keepModelIDs, modelIDs))
        
        #getting fitnesses from the model
        fitness = curGenerationModel['fitness']
        fitnessIDX = np.argsort(fitness)[::-1]
        
        #Checking for break conditions
        if os.path.exists('break_after_loop'): break
        if os.path.exists('debug_after_loop'):
            try:
                os.remove('debug_after_loop')
            except:
                pass
            pdb.set_trace()
        
        
        #start selective breeding
        #First the children that are kept for the subpopulation are put into a seperate variable
        if config.GAConfDict['mode'] == 'subpopulation':
            if populationSize==subPopulationNo:
                keepChildren=param.multiParam()
                keepChildren.grid=np.array([])
            else:
                keepChildren = model.modelGrid(paramList=curGenerationModel.
                                grid[fitnessIDX[:(populationSize-subPopulationNo)]],
                                origSpec=origSpec)
                keepModelIDs = modelIDs[fitnessIDX[:(populationSize - subPopulationNo)]]
                keepFitness = fitness[fitnessIDX[:(populationSize - subPopulationNo)]]
                keepChildren._initSpecs()
        
        elif config.GAConfDict['mode'] == 'elitism':
            keepChildren = model.modelGrid(paramList=curGenerationModel.
                            grid[fitnessIDX[:int(populationSize * config.GAConfDict['elitism'])]],
                            origSpec=origSpec)
            keepModelIDs = modelIDs[fitnessIDX[:(populationSize * config.GAConfDict['elitism'])]]
            keepFitness = fitness[fitnessIDX[:(populationSize * config.GAConfDict['elitism'])]]
            keepChildren._initSpecs()
        
        #Now we get the population that is used for breeding and submit it to the breed function
        breedPopulation=model.modelGrid(paramList=curGenerationModel.
                                       grid[fitnessIDX[:generationGapNo]],
                                       origSpec=origSpec)
        
        if config.GAConfDict['mode'] == 'subpopulation':
            curGenerationSet=breed(breedPopulation,
                  popNum=subPopulationNo, select=selectFunc, cross=crossFunc)
        elif config.GAConfDict['mode'] == 'elitism':
            curGenerationSet=breed(breedPopulation,
                  popNum=int(populationSize*(1-config.GAConfDict['elitism'])), select=selectFunc, cross=crossFunc)
        
        del curGenerationModel
        
        #Time is kept
        ficaTime=time.time()
    
        #Network check for available nodes
        gws.checkAvailability()
        pp.pprint(gws.availability)
        
        
        #Calculating the new generation with elauncher
        curGenerationModel = elauncher.cloudLaunch(curGenerationSet.paramGrid, gws.getAvailGateWays(), origSpec=origSpec)
        
        conn.commit()
        
        #Printing time statements
        print "Took %.3f for the fica runs" % (time.time()-ficaTime)
        print "Took %.3f seconds for last loop" % (time.time()-curTime)
        
    curs.execute('update GA_RUN set END_TIME=? where ID=?', (datetime.datetime.now(), GARunID))
