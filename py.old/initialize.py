#initializing values for fit (taking more or less educated guesses )
import math,string,os,shutil
import numpy as np
import config
time2vph=lambda t:17648.5*math.exp(-0.0920127*t)+5436.47
def readW7Data(dataFile=None):
    if dataFile==None: dataFile=os.path.join(config.getMainConfigDir(),'w7.combined.dat')
    data=np.loadtxt(dataFile)
    colHeaders=file(dataFile).readlines()[1].split()
    selProperties=[5,3]
    selElem=[9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
             20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
             31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42]
    colDescriptor=[]
#   Manually setting property colDescriptors
    colDescriptor.append('vph')
    colDescriptor.append('dens')
#Automatically setting the element colDescriptors
    colDescriptor+=[string.upper(colHeaders[i][0])+colHeaders[i][1:] for i in selElem]
    dataDict=dict([(key,np.array(value)) for key,value in zip(colDescriptor,zip(*data[:,selProperties+selElem]))])
    dataDict['vph']/=1e5
    return dataDict

def getW7Comp(w7Data=None,t=None,vph=None,selElements=['C','Ca','O','Ne','Si','S','Sc','Ti','V','Cr','Fe','Ni56','Mg','Mn','Cu','Zn']):
    if w7Data==None:
        w7Data=readW7Data()
    if t!=None:
        curVph=time2vph(t)
    elif vph!=None:
        curVph=vph
    elif t==None and vph==None:
        raise Exception('Please specify a time or a photospheric velocity')
    #setting the max composition speed to 12000 otherwise elements become 0, which poses problems
    if curVph>12000: curVph=12000.
    print "Evaluating at photospheric velocity %s km/s on day %s from Explosion"%(curVph,t)
    vphIndex=w7Data['vph'].searchsorted(curVph)
    compDict={}
    getIntElem=lambda elem,index: np.trapz(w7Data[elem][vphIndex:]*w7Data['dens'][vphIndex:],w7Data['vph'][vphIndex:])/np.trapz(w7Data['dens'][vphIndex:],w7Data['vph'][vphIndex:])
    for elem in selElements:
        compDict.update({elem:getIntElem(elem,vphIndex)})
    if compDict.has_key('Ni56'):
        compDict['Ni0']=compDict['Ni56']
        compDict.pop('Ni56')
    compDict['Fe0']=compDict['Fe']
    compDict.pop('Fe')
    norm=sum(compDict.values())
    if norm <0.8: raise Exception('Total sum of elements is less than 80%')
    for key,value in compDict.items():
        compDict[key]=value/norm
    return compDict

def initProcDirs(paramGrid,baseDir='.'):
    if paramGrid.paramGrid.size==0: raise Exception('ParamGrid is empty')
    for paramSet in paramGrid.paramGrid.flatten():
        prepProcDir(paramSet,baseDir=baseDir)
    
def prepProcDir(param,baseDir):
    procPath=os.path.join(baseDir,param.targetDir)
    if os.path.exists(procPath):
        try:
            shutil.rmtree(procPath)
        except:
            print "Process Directory %s exists. Having trouble removing it"%procPath
    os.mkdir(procPath,0777)
    param.write2file(procPath)
        