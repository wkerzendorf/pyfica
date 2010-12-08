#!/usr/bin/env python
import numpy as np
import pp
import subprocess
import config
import pyfica.read
#import pyfica.geneticDalek
import time
import progressbar
import shutil
import pdb
import os
import geneticDalek
machineConf=config.getMachineConfigPP()
coala_machines=('mosura01','mosura02','mosura03','mosura04',
                'mosura05','mosura06','mosura07','mosura08',
                'mosura09','mosura10','mosura11','mosura12',
                'mosura13','mosura14','mosura15','mosura16',)

nodes=('mosura01','mosura02','mosura03','mosura04','mosura05','mosura06',
       'mosura09','mosura10','mosura11','mosura12','miner')

def __jobServerDestroy(js):
    """
    Destroy the parallel python job server
    NOTE: this is a bug in the pp code resulting in open FIFO pipes
    """
    # access job server methods for shutting down cleanly
    js._Server__exiting = True
    js._Server__queue_lock.acquire()
    js._Server__queue = []
    js._Server__queue_lock.release()
    #pdb.set_trace()
    for worker in js._Server__workers:
        worker.t.exiting = True
        
        try:
            print "trying to close the file"
            # add worker close()
            worker.t.close()
            os.kill(worker.pid, 0)
            os.waitpid(worker.pid, 0)
            print "closed the file"
        except:
            #pdb.set_trace()
            # something nasty happened
            pass
    for rworker in js._Server__rworkers:
        rworker.close()
    for rworker in js._Server__rworkers_reserved:
        rworker.close()
    for rworker in js._Server__rworkers_reserved4:
        rworker.close()
def launchServers(startPort=60000):
    nodes=[node for node in machineConf if machineConf[node]['use']]
    processList=[]
    ppNodeDict={}
    ppArgs='-n2'
    
    for i,node in enumerate(nodes):
        #killSysCMD='ssh %s killall ppserver.py'%node
        sysCMD='ssh -f %s "killall ppserver.py;%s %s -p %d"'%(node,machineConf[node]['ppserver_bin'],ppArgs,i+startPort)
        print sysCMD
        #subprocess.Popen(killSysCMD,shell=True,stdout=-1)
        processList.append(subprocess.Popen(sysCMD,shell=True,stdout=-1))
        ppNodeDict[node]=i+startPort
    return processList,ppNodeDict
def getModelGrid(modelParams,nodeDict,basePath=None,applyFunc=None):
    noModel=len(modelParams)
    if basePath==None:
        basePath=config.getAutoDir()
    origSpec=config.getOrigSpec()
    nodes=[node for node in machineConf if machineConf[node]['use']]
    def keyer(x): return (not x.startswith('mosura')), x
    nodes=sorted(nodes, key=keyer)
    rload=getRelativeLoad(nodes)
    print "%s"%(zip(nodes,rload))
    nodeList=["%s:%d"%(item,nodeDict[item]) for item in getOptimalNodes(noModel,nodes,rload)]
    
    if len(nodeList)==0: raise Exception('No free nodes')
    
    print "Using these nodes: %s "%','.join(nodeList)
    js=pp.Server(ppservers=tuple(nodeList))
    js.set_ncpus(0)
    modelGrid=[]
    for param in modelParams:
        modelGrid.append(js.submit(runModel,(param,basePath,machineConf,origSpec),(getNodeName,),('subprocess','os','pyfica.read','shutil','numpy','pyfica.geneticDalek')))
    
    pbar=progressbar.ProgressBar(maxval=noModel).start()
    while True:
        currentStatus=[item.finished for item in modelGrid]
        if all(currentStatus): break
        pbar.update(np.sum(currentStatus))
        time.sleep(5)
    del js
    return [item() for item in modelGrid]
    
def getGenModelGrid(js,modelParams,fitness,basePath=None):
    if basePath==None:
        basePath=config.getAutoDir()
    js.set_ncpus(0)
    modelGrid=[]
    

def runModel(modelParam,basePath,machineConf,origSpec):
    currentNode=getNodeName()
    fica_bin=machineConf[currentNode]['fica_bin']
    ficaWorkDir=os.path.join(basePath,modelParam.targetDir)
    try:
        os.makedirs(ficaWorkDir)
    except:
        shutil.rmtree(ficaWorkDir)
        os.makedirs(ficaWorkDir)
    modelParam.write2file(ficaWorkDir)
    #print "Running %s on %s in dir %s"%(fica_bin,currentNode,ficaWorkDir)
    proc=subprocess.Popen([fica_bin],stdout=-1,cwd=ficaWorkDir)
    #Saving Log File:
    file(os.path.join(ficaWorkDir,'fica.log'),'w').write(proc.stdout.read())
    #print "tested %s"%proc.stdout.read()
    #return "I have run %s on %s"%(basePath,getNodeName())
    proc.stdout.close()
    model=pyfica.read.model(modelParam,ficaWorkDir)
    pyfica.geneticDalek.calcModelFitness(model,origSpec)
    return model

class genRunModel(object):
    pass
def getNodeName():
    p=subprocess.Popen(['uname -n'],stdout=-1,shell=True)
    return p.stdout.read().strip()
    
def getRelativeLoad(nodes):
    rLoads=[]
    processes=[]
    for node in nodes:
        p=subprocess.Popen(['ssh',node,'vmstat'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
        processes.append(p)
    print "Waiting for results"
    for i,p in enumerate(processes):
        r,b=map(int,p.stdout.readlines()[2].split()[:2])
        p.stdout.close()
        p.stderr.close()
        rLoads.append(r/float(machineConf[nodes[i]]['processors']))
    return rLoads
def getOptimalJS(noJobs,nodeDict):
    nodes=[node for node in machineConf if machineConf[node]['use']]
    def keyer(x): return (not x.startswith('mosura')), x
    nodes=sorted(nodes, key=keyer)
    rload=getRelativeLoad(nodes)
    print "%s"%(zip(nodes,rload))
    nodeList=["%s:%d"%(item,nodeDict[item]) for item in getOptimalNodes(noJobs,nodes,rload)]
    assert nodeList>0
    #if len(nodeList)==0: raise Exception('No free nodes')
    print "Using these nodes: %s "%','.join(nodeList)
    js=pp.Server(ppservers=tuple(nodeList))
    js.set_ncpus(0)
    return js
def getOptimalNodes(noJobs,nodes,rLoad,exludeList=None,includeList=None):
    
    sortLoadID=np.argsort(rLoad)
    curNoProcessors=0
    useNodes=[]
    
    
    for nodeID in sortLoadID:
        if rLoad[nodeID]>=0.6: continue
        #if nodes[nodeID].startswith('mosura'): continue
        useNodes.append(nodes[nodeID])
        
        curNoProcessors+=machineConf[nodes[nodeID]]['processors']
        if curNoProcessors>noJobs: break
    return useNodes