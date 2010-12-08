#The neural network part of dalek
import fit
import fitelem
import pdb
import numpy as np
import itertools
import cPickle
import config
import geneticDalek
from ffnet import ffnet, mlgraph
from glob import glob
import os
import re
import launcherSteps
import itertools
import shutil
def readMultiModelParams(pattern='*ts*.pkl',bestFit=None):
#if bestFit==None:
#        bestFit=cPickle.load(file('bestFit.pkl'))
    
    params=[]
    features=[]
    sources=[]
    for modelName in np.sort(glob(pattern)):
        #if bestFit==None:
        #    bestFit=cPickle.load(file(re.sub('\.ts\d+\.','.bf.',modelName)))
        model=cPickle.load(file(modelName))
        params.append([item.param for item in model.grid])
        features.append(getInputFromModel(model))
        sources.append([(modelName,i) for i in range(len(model.grid))])
    params=list(itertools.chain(*params))
    features=np.vstack(features)
    sources=list(itertools.chain(*sources))
    return params,features,sources
def getParamsFromMultiModel(pattern):
    id=[]
    params=[]
    for modelName in glob(pattern):
        model=cPickle.load(file(modelName))
        tmpParams=getParamsFromModel(model)
        params.append(tmpParams)
        lenModel=len(tmpParams)
        id.append(zip([modelName]*lenModel,range(lenModel)))
    return list(itertools.chain(*id)),list(itertools.chain(*params))
    
def getFeaturesFromBestFit(bestFit):
    merits=[]
    i=0
    origSpec=bestFit.origSpec
    aSpec,subSpec,llist,w,t=bestFit['aspec'],bestFit['subspec'],bestFit['llist'],bestFit['w'],bestFit['t']
    addSpec=aSpec+origSpec
    intTotal=fit.getInt(subSpec)/fit.getInt(addSpec)
    intUV=fit.getUVInt(subSpec)/fit.getUVInt(addSpec)
    slope=fit.getSlope(subSpec)
    chiSquared=fit.getChiSquared(subSpec)
    modChiSquared=fit.getModChiSquared(subSpec)
    elementMerits=[]
    for element in ['C','O','Na','Mg','Si','S','Ca','Ti','Cr','Ni','Fe']:
        try:
            lmerits,weights=fitelem.getElementMerit(element,llist,subSpec)
            elementMerits.append(np.average(lmerits,weights=weights))
            #pdb.set_trace()
        except:
            elementMerits.append(-1)
    allMerit=[intTotal,intUV,chiSquared,modChiSquared,slope,w,t]+elementMerits
    print "modincl"
    tmpParam=bestFit.param
    #tmpParam.update(modelGrid.grid[j].param.comp.data)
    merits.append((np.array(allMerit),tmpParam))
    return merits
def getParamsFromModel(modelGrid):
    merits=[]
    i=0
    origSpec=modelGrid.origSpec
    for j,(aSpec,subSpec,llist,w,t) in enumerate(zip(modelGrid['aspec'],modelGrid['subspec'],modelGrid['llist'],modelGrid['w'],modelGrid['t'])):
        print i
        i+=1
        #subSpec=model['subspec']
        addSpec=aSpec+origSpec
        intTotal=fit.getInt(subSpec)/fit.getInt(addSpec)
        intUV=fit.getUVInt(subSpec)/fit.getUVInt(addSpec)
        slope=fit.getSlope(subSpec)
        elementMerits=[]
        for element in ['C','O','Na','Mg','Si','S','Ca','Ti','Cr','Ni','Fe']:
            try:
                lmerits,weights=fitelem.getElementMerit(element,llist,subSpec)
                elementMerits.append(np.average(lmerits,weights=weights))
                #pdb.set_trace()
            except:
                elementMerits.append(-1)
        allMerit=[intTotal,intUV,slope,w,t]+elementMerits
        tmpParam=modelGrid.grid[j].param
        #tmpParam.update(modelGrid.grid[j].param.comp.data)
        merits.append((np.array(allMerit),tmpParam))
    return merits
def getTargetFromSet(params,ids,targetFunc):
    bestFitNames=list(set([re.search('(.+)ts\d+\.pkl',item).groups()[0] for item in set(zip(*ids)[0])]))
    bestFitDict=dict([(item+'bf.pkl',cPickle.load(file(item+'bf.pkl'))) for item in bestFitNames])
    return targetFunc(params,ids,bestFitDict)
    #weights=[C:0.01,Na:0.01,Mg:0.1,Ti:0.01,Cr:0.01,Fe0:0.1,Ni0:0.01]
    return np.nansum(targetFunc(params,ids,bestFitDict),axis=1)/16.

def getFeatureDistance(features,sources,metricFunc=None):
    bestFitDict={}
    featureDist=[]
    for feature,source in zip(features,sources):
        bestFitName=re.sub('\.ts\d+\.','.bf.',source[0])
        if bestFitDict.has_key(bestFitName): bestFit=bestFitDict[bestFitName]
        else:
            bestFitDict[bestFitName]=getFeaturesFromBestFit(cPickle.load(file(bestFitName)))[0][0]
            bestFit=bestFitDict[bestFitName]
        bestFitFeature=bestFit
        targetValues=feature
        if metricFunc==None:
            metricFunc=lambda x,y:x/y-1
        
        featureDist.append([metricFunc(bestFitFeatureValue,targetValue) for bestFitFeatureValue,targetValue in  zip(bestFitFeature,targetValues)])
    return featureDist

def getParamDistance(params,sources,metricFunc=None,selParams=['log_lbol','v_ph','C','O','Mg','Si','S','Ca','Ti','Cr','Ni0','Fe0']):
    bestFitDict={}
    paramDist=[]
    for param,source in zip(params,sources):
        bestFitName=re.sub('\.ts\d+\.','.bf.',source[0])
        if bestFitDict.has_key(bestFitName): bestFit=bestFitDict[bestFitName]
        else:
            bestFitDict[bestFitName]=cPickle.load(file(bestFitName))
            bestFit=bestFitDict[bestFitName]
        paramValues=[param[paramName] for paramName in selParams]
        paramValues[0]=10**paramValues[0]
        targetValues=[bestFit[paramName] for paramName in selParams]
        targetValues[0]=10**targetValues[0]
        if metricFunc==None:
            metricFunc=lambda x,y:x/y-1
        
        paramDist.append([metricFunc(paramValue,targetValue) for paramValue,targetValue in  zip(paramValues,targetValues)])
    return paramDist
        
def calcTargetParamRelDistance(params,ids,bestFitDict,selParams=['log_lbol','v_ph','C','O','Mg','Si','S','Ca','Ti','Cr','Ni0','Fe0']):
    #Calculates the target with relative distance to features
    bestFitTarget={}
    targets=[]
    for bestFitName,bestFitModel in bestFitDict.items():
        bestFitTarget[bestFitName]=np.array([bestFitModel[key] for key in selParams])
    for id,param in zip(ids,params):
        bestFitName=re.search('(.+)ts\d+\.pkl',id[0]).groups()[0]+'bf.pkl'
#        pdb.set_trace()
        tmpParam=np.array([param[1][key] for key in selParams])
        tmpTarget=abs(tmpParam-bestFitTarget[bestFitName])/(tmpParam+bestFitTarget[bestFitName])
        
        tmpTarget=abs(tmpParam/bestFitTarget[bestFitName]-1)
        #if any([np.isnan(item) for item in tmpTarget]): pdb.set_trace()
        targets.append(tmpTarget)
    return np.array(targets)

def calcTargetFeatureRelDistance(params,ids,bestFitDict):
    #Calculates the target with relative distance to features
    bestFitTarget={}
    targets=[]
    for bestFitName,bestFitModel in bestFitDict.items():
        origSpec=bestFitModel.origSpec
        aSpec=bestFitModel.aSpec
        subSpec=fit.getSubSpec(aSpec,origSpec)
        addSpec=aSpec+origSpec
        intTotal=fit.getInt(subSpec)/fit.getInt(addSpec)
        intUV=fit.getUVInt(subSpec)/fit.getUVInt(addSpec)
        slope=fit.getSlope(subSpec)
        elementMerits=[]
        for element in ['C','O','Na','Mg','Si','S','Ca','Ti','Cr','Ni','Fe']:
            try:
                lmerits,weights=fitelem.getElementMerit(element,bestFitModel.llist,subSpec)
                elementMerits.append(np.average(lmerits,weights=weights))
                #pdb.set_trace()
            except:
                elementMerits.append(-1)
        allMerit=[intTotal,intUV,slope,bestFitModel['w'],bestFitModel['t']]+elementMerits
        bestFitTarget[bestFitName]=allMerit
    for id,param in zip(ids,params):
        bestFitName=re.search('(.+)ts\d+\.pkl',id[0]).groups()[0]+'bf.pkl'
        #tmpTarget=[abs(iparam-itarget)/abs(iparam+itarget) for iparam,itarget in zip(param[0],bestFitTarget[bestFitName])]
        tmpTarget=[abs(iparam/itarget)-1 for iparam,itarget in zip(param[0],bestFitTarget[bestFitName]) if itarget!=0]
        targets.append(tmpTarget)
    return np.array(targets)
def getInputFromModel(modelGrid):
    merits=[]
    i=0
    origSpec=modelGrid.origSpec
    for aSpec,subSpec,llist,w,t in zip(modelGrid['aspec'],modelGrid['subspec'],modelGrid['llist'],modelGrid['w'],modelGrid['t']):
        print i
        i+=1
        #subSpec=model['subspec']
        addSpec=aSpec+origSpec
        intTotal=fit.getInt(subSpec)/fit.getInt(addSpec)
        intUV=fit.getUVInt(subSpec)/fit.getUVInt(addSpec)
        slope=fit.getSlope(subSpec)
        chiSquared=fit.getChiSquared(subSpec)
        modChiSquared=fit.getModChiSquared(subSpec)
        MAD=fit.getMAD(aSpec,origSpec)
        elementMerits=[]
        for element in ['C','O','Na','Mg','Si','S','Ca','Ti','Cr','Ni','Fe']:
            try:
                lmerits,weights=fitelem.getElementMerit(element,llist,subSpec)
                elementMerits.append(np.average(lmerits,weights=weights))
                #pdb.set_trace()
            except:
                elementMerits.append(-1)
        allMerit=[intTotal,intUV,chiSquared,modChiSquared,MAD,slope,w,t]+elementMerits
        print 'MAD'
        merits.append(allMerit)
    return np.array(merits)
def addProcSpec2BestFit(pattern='*bf.pkl'):
    for fname in glob(pattern):
        bestFit=cPickle.load(file(fname))
        bestFit.subSpec=fit.getSubSpec(bestFit['aspec'],bestFit.origSpec)
        bestFit.addSpec=fit.getSubSpec(bestFit['aspec'],bestFit.origSpec)
        tmpAspec=bestFit['aSpec']
        tmpAspec=tmpAspec.interpolate(xref=bestFit.origSpec.x)
        bestFit.divSpec=tmpAspec/bestFit.origSpec
        cPickle.dump(bestFit,file(fname,'w'))
    
def getDistanceOutput(bestFit,modelGrid,params=['lum','vph','C','O','Na','Mg','Si','S','Ca','Ti','Cr','Ni0','Fe0']):
    bestFitParam=np.array([bestFit[item] for item in params])
    bestFitParam[0]/=9.5
    bestFitParam[1]/=1e4
    distance=[]
    for model in modelGrid.grid:
        modelParam=np.array([model[item] for item in params])
        modelParam[0]/=9.5
        modelParam[1]/=1e4
        distance.append(modelParam-bestFitParam)
    return np.array(distance)
    
def getSingleDistance(distance):
    return np.sum(np.power(distance,2.),axis=1)
    
def createNNetwork(design):
    conec = mlgraph(design)
    return ffnet(conec)

def trainNetwork(net,inputs,targets,maxfun=0):
    net.train_genetic(inputs,targets,individuals=20,generations=500)
    net.train_tnc(inputs,targets,messages=1)
    return net

def createSingleTrainingSet(setSize,iter,SNName='2002bo'):
    t=config.getTimeFromExplosion()
    outnames=[]
    for i in range(iter):
        curGenerationSet=geneticDalek.createRandomParamSet(setSize)
        print "At Generation %s"%i
        curGenerationModel=launcherSteps.getModel(curGenerationSet)[0]
        outname='%st%.2f.ts%0d.pkl'%(SNName,t,i)
        print "Saving current generation to %s"%outname
        cPickle.dump(curGenerationModel,file(outname,'w'))
        outnames.append(outname)
    return outnames
def createTrainingSet(TSetName,pattern,TSetSize=250,iter=4):
    TSetDir=os.path.join('neuroTrain','trainingSet',TSetName)
    try:os.makedirs(TSetDir)
    except: print "WARNING: Training Set dir already exists. Does the TrainingSet exist?"
    for bestFit in glob(os.path.join('neuroTrain','trainingSet','*.bf.pkl')):
        shutil.copy(bestFit,TSetDir)
    for dir in glob(pattern):
        os.chdir(dir)
        outnames=createSingleTrainingSet(TSetSize,iter)
        for filePKL in outnames:
            shutil.move(filePKL,os.path.join('..',TSetDir))
        os.chdir('..')
    os.chdir(TSetDir)
    ids,targets=readMultiModelParams('*.ts*.pkl')
    cPickle.dump([ids,targets],file('%s.model.pkl'%TSetName,'w'))
    