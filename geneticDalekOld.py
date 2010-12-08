#!/usr/bin/env python

def evolve(generations=200,populationSize=150,continueGA=False):
    
    mutationRate=0.1
    mutationScale=0.125
    curMutationRate=mutationRate
    radDuration=0
    generationGapFraction=0.8
    generationGapNo=int(populationSize*generationGapFraction)
    subPopulationFraction=0.9
    subPopulationNo=int(populationSize*generationGapFraction*subPopulationFraction)
    
    if continueGA:
        x=[]
        x=[]
        y=[]
        y2=[]
        yerr=[]
        for fname in np.sort(glob('generation*.pkl')):
            x.append(int(re.search('\d+',fname).group()))
            fitness=np.loadtxt(fname.replace('.pkl','.dat'))
            y.append(np.mean(fitness))
            y2.append(np.max(fitness))
            yerr.append(np.std(fitness))
        plotSet.genPlotGenVSFitness(x,y,y2,yerr,outName='gen_vs_fitness_log.png')
        plotSet.genPlotGenVSFitness(x,y,y2,yerr,logPlot=False,outName='gen_vs_fitness_linear.png')
        initPopulation=cPickle.load(file('generation%03d.pkl'%np.max(x)))
        specMask=initPopulation.specFlag==0
        initPopulation.grid=initPopulation.grid[specMask]
        initPopulation.subSpec=initPopulation.subSpec[specMask]
        #initPopulation.addSpec=initPopulation.addSpec[specMask]
        initPopulation.divSpec=initPopulation.divSpec[specMask]
        initPopulation.specFlag=initPopulation.specFlag[specMask]
        firstGeneration=np.max(x)+1
        #populationSize=len(initPopulation.grid)
        fitness,curGenerationSet=breed(initPopulation,mutationRate=curMutationRate)    
    else:
        curGenerationSet=createRandomParamSet(populationSize)
        firstGeneration=0
        x=[]
        y=[]
        y2=[]
        yerr=[]
    try:
        os.remove('break_after_loop')
    except:
        pass
    try:
        os.remove('debug_after_loop')
    except:
        pass
    for i in range(firstGeneration,generations):
        print "At Generation %s @ mutationRate=%.4f"%(i,curMutationRate)
        curGenerationModel=launcherSteps.getModel(curGenerationSet)[0]
        specMask=curGenerationModel.specFlag==0
        curGenerationModel.grid=curGenerationModel.grid[specMask]
        curGenerationModel.subSpec=curGenerationModel.subSpec[specMask]
        #curGenerationModel.addSpec=curGenerationModel.addSpec[specMask]
        curGenerationModel.divSpec=curGenerationModel.divSpec[specMask]
        curGenerationModel.specFlag=curGenerationModel.specFlag[specMask]
        if 
        fname="generation%03d.pkl"%i
        print "Saving current generation to %s"%fname
        cPickle.dump(curGenerationModel,file(fname,'w'))
        oldGenerationModel=copy.deepcopy(curGenerationModel)
        print "Saved"
        if generations==1: break
        fitness=fitFunc(curGenerationModel)
        fitnessIDX=np.argsort(fitness)[::-1]
        keepChildren=curGenerationModel.grid[fitnessIDX[:(populationSize-subPopulationNo)]]
        breedPopulation=copy.deepcopy(curGenerationModel)
        breedPopulation.grid=curGenerationModel.grid[fitnessIDX[:generationGapNo]]
        curGenerationSet=breed(breedPopulation,popNum=subPopulationNo,mutationRate=curMutationRate,mutationScale=curMutationScale)
        
        
        np.savetxt('generation%03d.dat'%i,fitness)
        x.append(i)
        y.append(np.mean(fitness))
        y2.append(np.max(fitness))
        yerr.append(np.std(fitness))
        
        #Radiation Bursts
        if radDuration==0 and np.sum(np.diff(y2[-5:]))==0:
            curMutationRate=0.5
            radDuration=2
        if radDuration>0:
            if radDuration==1:
                radDuration=0
                curMutationRate=mutationRate
            else:
                radDuration-=1
                
        plotSet.genPlotSingleGeneration(oldGenerationModel,fitness,generation=i,no=10)
        plotSet.genPlotStatus(oldGenerationModel,pdfName='gen_status%03d.pdf'%i)
        plotSet.genPlotHistogram(oldGenerationModel,pdfName='gen_histogram%03d.pdf'%i)
        shutil.move('curGeneration.pdf','generation%03d.pdf'%i)
        plotSet.genPlotGenVSFitness(x,y,y2,yerr,outName='gen_vs_fitness_log.png')
        plotSet.genPlotGenVSFitness(x,y,y2,yerr,logPlot=False,outName='gen_vs_fitness_linear.png')
        del oldGenerationModel
        if os.path.exists('break_after_loop'): break
        if os.path.exists('debug_after_loop'):
            try:
                os.remove('debug_after_loop')
            except:
                pass