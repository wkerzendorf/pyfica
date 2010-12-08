import numpy as np
import itertools
import pdb
def getElementLList(element,llist):
    selLList=llist[llist['atom']==element.upper()]
    #selLList=selLList[selLList['eqw']>threshold]
    #if selLList.size==0: raise Exception('No lines found for Element %s'%element)
    return selLList
def getSetLList(element,MG):
    return [getElementLList(element,item) for item in MG['llist']]
def getRestWL(element,MG,cutThreshold=1):
    llists=getSetLList(element,MG)
    restWL=list(itertools.chain(*[item['rest'] for item in llists]))
    return [value for value, group in itertools.groupby(sorted(restWL))\
        if (len(list(group)) > cutThreshold) and ((value>MG.origSpec.x.min()) and value<MG.origSpec.x.max())]

def getLineMerit(lines,llist,spec,samples=np.array([-50.,50.]),lowerLimit=None,upperLimit=None,debug=False):
    merits=[]
    eqwWeights=[]
    if lowerLimit==None:
        lowerLimit=spec.x.min()
    if upperLimit==None:
        upperLimit=spec.x.max()
    if debug==True:
        regionList=[]
    for line in lines:
        selLList=llist[llist['rest']==line]
        if selLList.size>0:
            selLine=selLList[np.argmax(selLList['eqw'])]
            
            sampleSlice=slice(*(selLine['shift']+samples))
            if lowerLimit<sampleSlice.start and upperLimit>sampleSlice.stop:
                if debug==True: regionList.append(selLine['shift']+samples)
                merits.append(np.mean(spec[sampleSlice].y))
                eqwWeights.append(selLine['eqw'])
            else:
                merits.append(-1)
                eqwWeights.append(0)
        else:
            merits.append(-1)
            eqwWeights.append(0)
    if debug==True:
        return merits,eqwWeights,regionList
    else:
        return merits,eqwWeights

def interpolateMerits(merits):
    for i in np.arange(merits.shape[1]):
        line=merits[:,i]
        if all(line==-1): continue
        start = None
        end=None
        for i, item in enumerate(line):
            if start is None and item != -1:
                start = i
            if start is not None and item == -1 and end == None:
             end = i
            if start is not None and item !=-1: end=None
        if start!=None: line[:start]=line[start]
        if end!=None: line[end:]=line[end-1]
def getDiffWeights(merits):
    dWeights=[]
    for i in np.arange(merits.shape[1]):
        line=merits[:,i]
        dline=list(np.diff(line))
        dline.insert(0,dline[0])
        dline.insert(-1,dline[-1])
        dWeights.append([np.mean((item,dline[i+1])) for i,item in enumerate(dline[:-1])])
    return np.array(dWeights).transpose()
        
def getSetMerits(element,MG):
    lines=getRestWL(element,MG)
    merits=[]
    eqwWeights=[]
    origSpec=MG.origSpec
    contOrigSpec=origSpec.fitContinuum(func='poly5',lsigma=3)
    #normOrigSpec=origSpec/contOrigSpec
    scale=np.median(np.diff(origSpec.x))
    normOrigSpec=origSpec.smoothg(15./scale)
    cont=np.mean([item.fitContinuum(func='poly5') for item in MG['divspec']])
    for llist,divSpec in zip(MG['llist'],MG['divspec']):
        #contDivSpec=divSpec.fitContinuum(func='poly5')      
        normSpec=divSpec/cont
        tmpMerits,tmpEqwWeights=getLineMerit(lines,llist,normSpec,lowerLimit=4000)
        merits.append(tmpMerits)
        eqwWeights.append(tmpEqwWeights)
    merits=np.array(merits)
    interpolateMerits(merits)
    diffWeights=getDiffWeights(merits)
    weights=np.array(eqwWeights)*np.abs(diffWeights)
    return merits,np.array(eqwWeights),np.abs(diffWeights)
    return np.average(merits,weights=weights,axis=1)
        
def getElementMerit(element,llist,spec,threshold=5,samples=np.array([-50,50.])):
    cont=spec.fitContinuum(func='poly5')
    normSpec=spec/cont
    selLList=llist[llist['atom']==element.upper()]
    selLList=selLList[selLList['eqw']>threshold]
    if selLList.size==0: raise Exception('No lines found to build merit')
    lineMerit=[]
    eqwWeights=[]
    #return selLList
    for line in selLList:
        sampleSlice=slice(*(line[1]+samples))
        #print element.lower()
        if element.lower()=='ca':
            caSample=np.array([-60,0])
#            print "Ca is a special element and will be treated specially:"
#            print "Integrating sample %s"%(line[1]+caSample)
            sampleSlice=slice(*(line[1]+caSample))
        if normSpec.x.min()<sampleSlice.start and normSpec.x.max()>sampleSlice.stop:
            lineMerit.append(np.mean(normSpec[sampleSlice].y))
            eqwWeights.append(line[0])
    if len(lineMerit)==0: raise Exception('All lines outside of the observed spectrum')
    return lineMerit,eqwWeights

def getIGEMerit(llist,spec,threshold=5,samples=np.array([0.,30.])):
    IGEs=['Ni','Fe','Sc','Ti','V','Cr','Mn','Cu','Zn']
    merits=[]
    for element in IGEs:
        try:
            merits.append(getElementMerit(element,llist,spec,threshold=threshold,samples=samples))
        except:
            print "Could not find lines for element %s"%element
    return merits

def showIGEAbundances(param):
    IGEs=['Ni','Fe0','Fe','Sc','Ti','V','Cr','Mn','Cu','Zn']
    for element in IGEs:
        print "%s: %.7f"%(element,param[element])


def getGridElementMerit(element,modelGrid):
    llists=modelGrid['llist']
    divSpecs=modelGrid['divspec']
    merits=[]
    for i,llist in enumerate(llists):
        try:
            merits.append(getElementMerit(element,llist,divSpecs[i],threshold=0))
        except:
            print "Couldn't find lines for %s=%s"%(element,modelGrid[element][i])
            merits.append(-1)
    return np.array(merits).reshape(divSpecs.shape)
    
def calcRelAbundance(element,relAbund,comp,relElement='O'):
    relElementAbund=1/(np.sum([item[1]/comp[relElement] for item in comp.data.items() if not (item[0]==element or item[0][-1]=='0')])+relAbund)
    print relElementAbund
    return relAbund*relElementAbund
    
#getGridElementMerit=np.vectorize(lambda item:getElementMerit('Si',item['llist'],item['divspec']))





#from bisect import bisect
#from scipy import polyval, polyfit
#from numpy import vstack, diff
#from .fileio import dicafile, compfile
#from os import path
#from .util import median_spectra, gaussian_spectra, getSpectrumScale, diffSpectrum, getBinMerits, findLineEdges
#
#
#
#
#
#
#def getMaxIGE(param):
#    print
#def siSection(divspect, sbib):
#    siWl=[item[1] for item in sbib['llist'] if item[3].lower()=='si' and (int(item[2])==6347 or int(item[2])==6371)]
#    indSiW=[bisect(divspect[:,0],siWl[0]+item) for item in [-500, -300, 300, 500]]
#    siSection=divspect[indSiW[0]:indSiW[3]]
#    return siSection
#def siMerit(models, smoothSize=100):
#    bins=[findLineEdges(models[0]['origspect'], 6100, smoothSize)]
#    abund=[model['comp']['Si'] for model in models]
#    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
#    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
#    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
#    
#    siWl=[item[1] for item in sbib['llist'] if item[3].lower()=='si' and (int(item[2])==6347 or int(item[2])==6371)]
#    indSiW=[bisect(divspect[:,0],siWl[0]+item) for item in [-500, -300, 300, 500]]
#    siMargin=vstack((divspect[indSiW[0]:indSiW[1]], divspect[indSiW[2]:indSiW[3]]))
#    polyParam=polyfit(siMargin[:, 0], siMargin[:, 1], 1)
#    siSection=divspect[indSiW[0]:indSiW[3]]
#    normSiSection=siSection[:, 1]/polyval(polyParam,siSection[:,0])
#    integral=sum([(siSection[i+1, 0]-siSection[i, 0])*normSiSection[i] for i in range(len(siSection[:, 0])-1)])/\
#        abs(siSection[-1, 0]-siSection[0, 0])
#    return integral
#
#def getCaMerit(models, smoothSize=200):
#    bins=[findLineEdges(models[0]['origspect'], 3790, 100)]
#    abund=[model['comp']['Ca'] for model in models]
#    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
#
#def getSMerit(models, smoothSize=100):
#    bins=[findLineEdges(models[0]['origspect'], 5350, 100)]
#    abund=[model['comp']['S'] for model in models]
#    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
#def getMgMerit(models, smoothSize):
#    bins=[findLineEdges(models[0]['origspect'], 4250, 50)]
#    abund=[model['comp']['Mg'] for model in models]
#    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
#def doDecay(confPath=''):
#    t1=6.1
#    t2=77.27    
#    t=dicafile(path.join(confPath,'dica.dat')).read_data()['t']
#    comp=compfile(path.join(confPath,'comp.ind')).read_data()
#    initNiAbund=comp['Ni']
#    decayNiAbund=initNiAbund*(pow(2, -(t/t1)))
#    decayCoAbund=(t2/(t1-t2))*initNiAbund*(pow(2, -t/t1)-pow(2, -t/t2))
#    decayFeAbund=initNiAbund-decayCoAbund-decayNiAbund
#    comp['Ni']=decayNiAbund
#    comp['Co']=decayCoAbund
#    comp['Fe']+=decayFeAbund
#    compfile(path.join(confPath,'comp.ind'), 'w').write_data(comp)
