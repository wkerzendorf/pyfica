from bisect import bisect
from scipy import polyval, polyfit
from numpy import vstack, diff
from .fileio import dicafile, compfile
from os import path
from .util import median_spectra, gaussian_spectra, getSpectrumScale, diffSpectrum, getBinMerits, findLineEdges

def siSection(divspect, sbib):
    siWl=[item[1] for item in sbib['llist'] if item[3].lower()=='si' and (int(item[2])==6347 or int(item[2])==6371)]
    indSiW=[bisect(divspect[:,0],siWl[0]+item) for item in [-500, -300, 300, 500]]
    siSection=divspect[indSiW[0]:indSiW[3]]
    return siSection
def siMerit(models, smoothSize=100):
    bins=[findLineEdges(models[0]['origspect'], 6100, smoothSize)]
    abund=[model['comp']['Si'] for model in models]
    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
    
    siWl=[item[1] for item in sbib['llist'] if item[3].lower()=='si' and (int(item[2])==6347 or int(item[2])==6371)]
    indSiW=[bisect(divspect[:,0],siWl[0]+item) for item in [-500, -300, 300, 500]]
    siMargin=vstack((divspect[indSiW[0]:indSiW[1]], divspect[indSiW[2]:indSiW[3]]))
    polyParam=polyfit(siMargin[:, 0], siMargin[:, 1], 1)
    siSection=divspect[indSiW[0]:indSiW[3]]
    normSiSection=siSection[:, 1]/polyval(polyParam,siSection[:,0])
    integral=sum([(siSection[i+1, 0]-siSection[i, 0])*normSiSection[i] for i in range(len(siSection[:, 0])-1)])/\
        abs(siSection[-1, 0]-siSection[0, 0])
    return integral

def getCaMerit(models, smoothSize=200):
    bins=[findLineEdges(models[0]['origspect'], 3790, 100)]
    abund=[model['comp']['Ca'] for model in models]
    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)

def getSMerit(models, smoothSize=100):
    bins=[findLineEdges(models[0]['origspect'], 5350, 100)]
    abund=[model['comp']['S'] for model in models]
    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
def getMgMerit(models, smoothSize):
    bins=[findLineEdges(models[0]['origspect'], 4250, 50)]
    abund=[model['comp']['Mg'] for model in models]
    return abund, getBinMerits(models, fBins=bins, func='mean', doProc=None)
def doDecay(confPath=''):
    t1=6.1
    t2=77.27    
    t=dicafile(path.join(confPath,'dica.dat')).read_data()['t']
    comp=compfile(path.join(confPath,'comp.ind')).read_data()
    initNiAbund=comp['Ni']
    decayNiAbund=initNiAbund*(pow(2, -(t/t1)))
    decayCoAbund=(t2/(t1-t2))*initNiAbund*(pow(2, -t/t1)-pow(2, -t/t2))
    decayFeAbund=initNiAbund-decayCoAbund-decayNiAbund
    comp['Ni']=decayNiAbund
    comp['Co']=decayCoAbund
    comp['Fe']+=decayFeAbund
    compfile(path.join(confPath,'comp.ind'), 'w').write_data(comp)
