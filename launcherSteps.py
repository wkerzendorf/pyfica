#functions to launch parameter sets:

import os
import shutil
import pickle
import itertools

import numpy as np
import config,initialize,launcher,param,read

#Pickling is done to work on two different systems in the initial stage of the autofitter
def initLaunchLum(doPickle=True):
    initLumRange=np.linspace(9,10,10)
    return launchLum(initLumRange,step=1,dirName='initLum')
def launchLum(lumRange,step=None,dirName=None):
    autoPath=config.getAutoDir()
    if dirName==None: dirName='lum'
    path=os.path.join(autoPath,dirName)
    if os.path.exists(path): shutil.rmtree(path)
    os.makedirs(path)
    lumParam=param.multiParam()
    #lumParam.initParam['lum']=9.0
    #lumParam.initParam.dica.lockTemp=True
    lumParam['lum']=lumRange
    runDirs=initialize.initProcDirs(lumParam,path)
    launcher.multiDistLaunch(runDirs)
    return runDirs,lumParam

def launchTriCycle(lumRange,vphRange,IGERange,IGEElement,initParam=None,procPath='tri0',step=None):
    lumParam=param.multiParam(initParam)
    lumParam['lum']=lumRange
    
    vphParam=param.multiParam(initParam)
    vphParam['vph']=vphRange
    
    IGEParam=param.multiParam(initParam)
    IGEParam.initParam.comp.lockIGE=True
    IGEParam[IGEElement]=IGERange
    
    
    if initParam!=None:
        curParam=param.multiParam(initParam)
        curParam.initParam.targetDir='%s_param'%procPath
        curParam.paramGrid=np.array([curParam.initParam])
        
    basePath=os.path.join(config.getAutoDir(),procPath)
    if os.path.exists(basePath): shutil.rmtree(basePath)
    os.makedirs(basePath)
    if initParam!=None:
        return getModel(lumParam,vphParam,IGEParam,curParam,basepath=basePath)
    else:
        return getModel(lumParam,vphParam,IGEParam,basepath=basePath)
    
def launchElementCycle(rangeDict,initParam=None,procPath='elemCycle'):
    elementParams=[]
    for key in rangeDict:
        elementParams.append(param.multiParam(initParam))
        elementParams[-1][key]=rangeDict[key]
    
    param.multiParam(initParam)
    basePath=os.path.join(config.getAutoDir(),procPath)
    
    if initParam!=None:
        curParam=param.multiParam(initParam)
        curParam.initParam.targetDir='%s_param'%procPath
        curParam.paramGrid=np.array([curParam.initParam])
    
    if os.path.exists(basePath): shutil.rmtree(basePath)
    os.makedirs(basePath)
    
    return getModel(*(elementParams+[curParam]),**dict(basepath=basePath))
    
def getModel(*params,**kwargs):
    if kwargs.has_key('basepath'):
        basePath=kwargs['basepath']
    else:
        basePath=config.getAutoDir()
    runDirs=[]
    for param in params:
        runDirs.append(initialize.initProcDirs(param,basePath))
    launcher.multiDistLaunch(list(itertools.chain(*runDirs)))
    modelGrids=[]
    for runDir,param in zip(runDirs,params):
        modelGrids.append(read.modelGrid(runDir,param))
    return modelGrids
    