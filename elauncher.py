import sys
import time
import cPickle as pickle
import pdb
import model
import config
import execnet
import genWorker
import weakref
import shutil
import os
import genFitness
import numpy as np

machineConfig = config.getMachineConfigExecNet()

checkUpTimeStep = 2 #seconds
#currentGWs=openGateWays()

class gateways(object):
    def __init__(self):
        self.availability = {}
        self.availThresh = 0.6
        self.loadavg = {}
        #self.machines=['mosura05','mosura02','mosura10','mutant']
        self.exclude='mosura'
        self.machines=[item for item in machineConfig.keys() if (machineConfig[item]['use'] and (not item.startswith(self.exclude)))]
        #self.machines += ['mosura16']
        self.machines+=['mosura14','mosura15','mosura10','mosura11','mosura12','mosura13','mosura16']
        self.machines+=['mosura01','mosura02','mosura03',
                        'mosura04','mosura05','mosura06',
                        'mosura07','mosura08','mosura09']
        print self.machines
        self.speed=[machineConfig[machine]['speed'] for machine in self.machines]
        self.checkProcNames=['ppns']
        self._openGateWays()
        self.checkAvailability()
    def _openGateWays(self):
        self.gateways=[]
        self.checkGWs={}
        self.machineGateways={}
        j=0
        for machine in self.machines:
            if machineConfig[machine]['use']:
                procs=machineConfig[machine]['processors']
                print "Opening %d to machine %s"%(procs,machine)
                checkGW=True
                for i in range(procs):
                    #implement nice level 
                    GWCommand='ssh=%s//python=%s'%\
                              (machine,machineConfig[machine]['python_bin'])
                    #one True is for self managed the next is for open
                    gw=execnet.makegateway(GWCommand)
                    if checkGW:
                        self.checkGWs[machine]=gw
                        checkGW=False
                    if not self.machineGateways.has_key(machine):
                        self.machineGateways[machine]=[]
                    self.machineGateways[machine].append([machine,gw,True,True])
                    self.gateways.append([machine,gw,True,True])
                    j+=1
                    
    def getAvailGateWays(self):
        availGateWays=[]
        for i in np.argsort(self.speed)[::-1]:
            machine=self.machines[i]
            cores,threads=self.availability[machine]
            for j in range(cores):
                #print "%d machine: %s threads %d"%(j,machine,threads)
                availGateWays.append([machine,self.machineGateways[machine][j][1],threads])
        return availGateWays
    
    
    def _checkVmstat(self,machine):
        remoteExec="import subprocess;\
                   vmstat=subprocess.Popen('vmstat',stdout=subprocess.PIPE)\
                   .stdout.readlines()[2];\
                   channel.send(map(int,vmstat.split()[:2]))"
        gw=self.checkGWs[machine]
        ch=gw.remote_exec(remoteExec)
        r,b=ch.receive()
        return r+b
    def _checkProcesses(self,machine,processName):
        processCheckCommand='ps -a |grep %s'%processName
        processCheck="import subprocess;\
                   procCheck=subprocess.Popen('%s',stdout=subprocess.PIPE,shell=True)\
                   .stdout.readlines();\
                   channel.send(len(procCheck))"%processCheckCommand
        gw=self.checkGWs[machine]
        ch=gw.remote_exec(processCheck)
        noProcs=ch.receive()
        return noProcs
    def checkAvailability(self):
        print "Checking Availability... ",
        for machine in self.machines:
            if machineConfig[machine]['use']:
                noCores=machineConfig[machine]['processors']
                usedCores=0
                for procName in self.checkProcNames:
                    usedCores+=self._checkProcesses(machine,procName)
                useCores=max((0,noCores-usedCores))
                self.availability[machine]=[useCores,1]
                if usedCores==0:
                    usedCores+=self._checkVmstat(machine)
                    useCores=max((0,noCores-usedCores))
                    self.availability[machine]=[useCores,machineConfig[machine]['threads']]
        print "done"
    def printAvailability():
        for machine in self.machines:
            pass
class modelList(list):
    pass
def cloudLaunch(params,gateways,origSpec=None,baseDir=None):
    if origSpec==None: origSpec=config.getOrigSpec(preProcess=True)
    if baseDir==None: baseDir = os.getcwd()
    noModels=len(params)
    models=modelList()
    proxyModels=weakref.proxy(models)
    ##CONSTANTS
    params=list(params)
    machines=zip(*gateways)[0]
    machineSlots=dict([[machine,0] for machine in machines])
    startTime=time.time()
    
    #Configuring which machines to use
    
    def callbackFica(confTuple):
        #execution of fica after the program has run
        protocol, i = confTuple[:2]
        gwConfig = gateways[i]
        machineName = gwConfig[0]
        
        #Reading model and appending
        if protocol == 'centralRead':
            ficaWorkDir=confTuple[2]
            model=getModel(ficaWorkDir,origSpec)
            model.machine = machineName
            proxyModels.append(model)
        elif protocol == 'localRead':
            model = pickle.loads(confTuple[2])
            model.machine = machineName
            proxyModels.append(model)
        
        
        threads = gwConfig[2]
        gw=gwConfig[1]
        machineSlots[machineName] -= 1
        
        try: param=params.pop(0)
        except IndexError: return
        machineSlots[machineName] += 1
        ch=launch(param, machineConfig[machine], threads,baseDir, origSpec, gw, i, genWorker)
        ch.setcallback(callbackFica)
    
    channels=[]
    
    for i,gwConfig in enumerate(gateways):
        try: param=params.pop(0)
        except IndexError: break
        machineName=gwConfig[0]
        gw=gwConfig[1]
        threads=gwConfig[2]
        machineSlots[machineName]+=1
        ch=launch(param,machineConfig[machine],threads,baseDir,origSpec,gw,i,genWorker)
        ch.setcallback(callbackFica)
        
    ########
    #Printing the current progress (works only on VT-100 compliant terminals)
    while True:
        for machineName in machineSlots:
            print "%s: %02d"%(machineName,machineSlots[machineName])
        lines=len(machineSlots)+1
        if len(models)==noModels:break
        if params!=[]:
            print "Next process to launch is %s/%s"%\
            (noModels-len(params)+1,noModels)
            sys.stdout.write("\x1b[%sA"%lines)
        else:
            print "Launched all parameter sets, Still %s items in queue"%\
            (noModels-len(models))
            sys.stdout.write("\x1b[%sA"%lines)
        
        time.sleep(checkUpTimeStep)
    wrefModels=weakref.ref(models)
    modelGrid=model.modelGrid(paramList=models,origSpec=origSpec)
    
    del models
    
    timedelta = time.time()-startTime
    
    print "Fica run took %.4f min and %.4f s/model" % (timedelta/60, timedelta/noModels)
    
    return modelGrid


def launch(param,configParam,threads,basePath,origSpec,gateway,gwid,remoteModule, protocol = 'localRead'):
    channel=gateway.remote_exec(remoteModule)
    python_path=configParam['python_path']
    ficaBin="%s %d"%(configParam['fica_bin'],threads)
    channel.send(python_path)
    #Establishing protocol
    channel.send(protocol)
    
    dica=param.dica.data
    for key in dica:
        if isinstance(dica[key],float):
            dica[key]=float(dica[key])
            
    comp=param.comp.data
    for key in comp:
        if isinstance(comp[key],float):
            comp[key]=float(comp[key])
            
    ficaWorkDir=os.path.join(basePath,param.targetDir)
    
    #sending for different protocols
    if protocol == 'centralRead':
        channel.send((dica,comp,ficaWorkDir,ficaBin,gwid))
    elif protocol == 'localRead':
        pickledOrigSpec = pickle.dumps(origSpec)
        channel.send((dica,comp,ficaWorkDir,ficaBin,gwid, pickledOrigSpec))
    return channel
def getModel(ficaWorkDir,origSpec):
    model=model.model(ficaWorkDir,origSpec=origSpec)
    model.fitness=genFitness.fitFunc(model)
    aSpec=model.aSpec.interpolate(xref=genFitness.origSpec.x)
    model.contOptical=genFitness.fitOpticalContinuum(aSpec)
    model.contIR=genFitness.fitIRContinuum(aSpec)
    model.contIRDiff=(model.contIR-genFitness.contIROrig)/(model.contIR+genFitness.contIROrig)
    model.contOpticalDiff=(model.contOptical-genFitness.contOpticalOrig)/(model.contOptical+genFitness.contOpticalOrig)
    model.lineFit=aSpec[:]
    irThresh=genFitness.irThresh
    model.lineFit[:irThresh]/=model.contOptical
    model.lineFit[irThresh:]/=model.contIR
    model.lineCheck=(((model.lineFit-genFitness.lineOrig)/(genFitness.lineOrig+model.lineFit)))
    model.lineCheck.y**=2
    shutil.rmtree(ficaWorkDir)
    return model