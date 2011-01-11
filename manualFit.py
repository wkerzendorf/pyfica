#!/usr/bin/env python
import os,shutil
import fileio,abund,fit,config,param,fit
from pyspec.spectrum import spectrum
from glob import glob
#import pylab,matplotlib
import numpy as np
#from matplotlib.backends.backend_pdf import PdfPages
#origspec=spectrum('origspect.dat')
#ficaBin=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin/bin/fica.exe')
#useMachines=['mithrandir-local']
#fig=pylab.figure(1)

import elauncher

def runManual(param, gateways=None):
    if gateways==None:
        gateways = elauncher.gateways()
    model = elauncher.cloudLaunch([param,], gateways.getAvailGateWays())
    return model[0]
    

def saveOldFit():
    os.system('cp spct.dat spct.old')
    os.system('cp dica.dat dica.old')
    os.system('cp comp.ind comp.old')
def onAtom(event):
#    print event.xdata,event.ydata,event.key
    if event.key=='a': printLines(event.xdata)
def printLines(wl,sample=80):
    id=getLastID()
    llist=getHistLList(id)
    min,max=llist['shift'].searchsorted(wl-sample),llist['shift'].searchsorted(wl+sample)
    llistSample=llist[min:max]
    sortedID=np.argsort(llistSample)
    print "------------------------------------------------"
    for id in sortedID[::-1]:
        print '%.3f %.3f %.3f %s %s'%tuple(llistSample[['eqw','shift','rest','atom','ion']][id])
    print "------------------------------------------------"
def printCurMerits():
    id=getLastID()
    aSpec=getHistSpec(id)
    aSpecOld=getHistSpec(id-1)
    divSpec=origspec/aSpec
    divSpecOld=origspec/aSpecOld
    subSpec=fit.getSubSpec(aSpec,origspec)
    subSpecOld=fit.getSubSpec(aSpecOld,origspec)
    checkElements=['Ni','Co','Fe','Ca','Ti','Cr','Si','S','O','C','Mg']
    #checkIons=[[]]
    print "General integral: %s last: %s"%(fit.getDiffIntBin(subSpec,1)[1],fit.getDiffIntBin(subSpecOld,1)[1])
    print "Integral slope: %s last: %s"%(fit.getIntSlope(subSpec),fit.getIntSlope(subSpecOld))
    print "UV excess: %s last: %s"%(fit.getUVInt(subSpec),fit.getUVInt(subSpecOld))
    print "UV comparison (UV to optical): %s last: %s"%(fit.getUVIntComp(subSpec),fit.getUVIntComp(subSpecOld))
    print "---------------- Element Ratios --------"
    for element in checkElements:
        try:
            print "Element %s: %.5f  last: %.5f "%(element,getElementMerit(element,divSpec,threshold=0),getElementMerit(element,divSpecOld,threshold=0,id=id-1))
        except Exception as inst:
            print "Problem getting %s merit: %s"%(element,inst)
    print "________________ Ion Ratios ------------"
    
def runElementRatio(dica,comp,element,ratios=[0.5,2],ax=None):
    abundance=comp[element]
    ratioSpec=[]
    for ratio in ratios+[1]:
        comp[element]=ratio*abundance
        runManual(dica,comp,origspec,backup=False)    
        ratioSpec.append(spectrum('spct.dat',usecols=(0,2)))
    if ax==None:
        fig=pylab.figure(1)
        fig.clf()
        ax=fig.add_subplot(111)
    ax.plot(origspec.x,origspec.y,lw=3,color='blue')
    for ratio,spec in zip(ratios+[1],ratioSpec):
        ax.plot(spec.x,spec.y,label='%s = %s'%(element,ratio*abundance))
    ax.legend()
        
def runMultiElement(dica,comp,elements=['Ni0','Fe0','Si','S','Ti','C','Mg'],ratios=[0.5,2],output='multi_element.pdf'):
    pdf = PdfPages(output)
    fig=pylab.figure(1)
    for element in elements:
        print '----------------------------------------------------------------'
        print "Doing Element %s"%element
        os.sytem('growlnotify -m Doing Element %s'%element)
        fig.clf()
        ax=fig.add_subplot(111)
        runElementRatio(dica,comp,element,ratios,ax=ax)
        fig.savefig(pdf,format='pdf')
    pdf.close()
    os.sytem('growlnotify --sticky -m Done with multielement')
def getElementMerit(element,divSpec,threshold=5,samples=np.array([0.,30.]),id=None):
    if id==None: id=getLastID()
    llist=getHistLList(id)
    cont=divSpec.fitContinuum(func='poly5')
    normDivSpec=divSpec/cont
    selLList=llist[llist['atom']==element.upper()]
    selLList=selLList[selLList['eqw']>threshold]
    if selLList.size==0: raise Exception('No lines found to build merit')
    lineMerit=[]
    eqwWeights=[]
    for line in selLList:
        sampleSlice=slice(*(line[1]+samples))
        #print element.lower()
        if element.lower()=='ca':
            caSample=np.array([-60,0])
#            print "Ca is a special element and will be treated specially:"
#            print "Integrating sample %s"%(line[1]+caSample)
            sampleSlice=slice(*(line[1]+caSample))
        
        if divSpec.x.min()<sampleSlice.start and divSpec.x.max()>sampleSlice.stop:
            lineMerit.append(np.median(normDivSpec[sampleSlice].y))
            eqwWeights.append(line[0])
    if len(lineMerit)==0: raise Exception('All lines outside of the observed spectrum')
    return np.average(lineMerit,weights=eqwWeights)
def getIonMerit(element,ion,divSpec,threshold=5,samples=np.array([0.,30.]),llist=None):
    if llist==None: llist=getCurLList()
    cont=divSpec.fitContinuum(func='poly5')
    normDivSpec=divSpec/cont
    selLList=llist[llist['atom']==element.upper()]
    selLList=selLList[selLList['ion']==ion.upper()]
    selLList=selLList[selLList['eqw']>threshold]
    if selLList.size==0: raise Exception('No lines found to build merit')
    lineMerit=[]
    eqwWeights=[]
    for line in selLList:
        sampleSlice=slice(*(line[1]+samples))
        
        if divSpec.x.min()<sampleSlice.start and divSpec.x.max()>sampleSlice.stop:
            lineMerit.append(np.median(normDivSpec[sampleSlice].y))
            eqwWeights.append(line[0])
    return np.average(lineMerit,weights=eqwWeights)
def plotHist(id=None):
    if id==None:
        id=getLastID()
    aspec=getHistSpec(id)
    aspecOld=getHistSpec(id-1)
    fig=pylab.figure(1)
    fig.clf()
    ax=fig.add_subplot(111)
    ax.plot(origspec.x,origspec.y,'k')
    ax.plot(aspec.x,aspec.y,'r')
    ax.plot(aspecOld.x,aspecOld.y,'g')
    
def plotElement(element,threshold=5,cmap=None,alpha=0.5):
    llist=getCurLList()
    selLList=llist[llist['atom']==element.upper()]
    ax=pylab.gca()
    normfunc=matplotlib.colors.LogNorm()
    normColors=normfunc(selLList['eqw'])
    for i,line in enumerate(selLList):
        if line['eqw']<threshold: continue
        ax.axvline(line['shift'],alpha=alpha,lw=4,color=cmap(normColors[i]))
    pylab.get_current_fig_manager().canvas.draw()
def plotCur():
    aspec=getCurSpectrum()
    fig=pylab.figure(1)
    fig.clf()
    ax=fig.add_subplot(111)
    ax.plot(origspec.x,origspec.y,'k')
    ax.plot(aspec.x,aspec.y,'r')
    
def setElemAbundNormToO(comp,elem,abund):
    comp=comp.copy()
    oldElemAbundance=comp[elem]
    deltaElemAbundance=oldElemAbundance-abund
    comp['O']+=deltaElemAbundance
    comp[elem]=abund
    return comp
def getCurLList():
    id=getLastID()
    #sbib=fileio.sbibfile('sbib.dat').read_data()
    return getHistLList(id)
def getCurSpectrum():
    return spectrum('spct.dat',usecols=(0,2))

def getHistSpec(id):
    return spectrum('hist/spct.bak%04d'%id,usecols=(0,2))
def printHist(dicaKeys=['log_lbol','v_ph'],compKeys=['Fe0','Fe','Ni0','Si','S','Cr','Ti','Mg'],max=10):
    #dica=fileio.dicafile('dica.dat').read_data()
    #comp=fileio.compfile('comp.ind').read_data()
    #print comp
    dicaStr=['%s=%s'%item for item in dica.items() if item[0] in dicaKeys]
    compStr=['%s=%s'%item for item in comp.items() if item[0] in compKeys]
    #print 'Current: '+' '.join(dicaStr+compStr)
    #print '------------------------------------'
    dicaFiles=glob('hist/dica.bak????')
    compFiles=glob('hist/comp.bak????')
    dicaFiles.sort()
    compFiles.sort()
    for dicaFile,compFile in zip(dicaFiles,compFiles)[::-1][:max]:
        id=int(dicaFile[-4:])
        dica=fileio.dicafile(dicaFile).read_data()
        comp=fileio.compfile(compFile).read_data()
        #print comp
        dicaStr=['%s=%s'%item for item in dica.items() if item[0] in dicaKeys]
        compStr=['%s=%s'%item for item in comp.items() if item[0] in compKeys]
        print 'id=%s |'%id+' '.join(dicaStr+compStr)
def backupFit(id=None):
    #backupFiles=['dica.dat', 'comp.ind', 'spct.dat', 'sbib.dat', 'stst.dat','yhea.dat']
    backupFiles=['dica.dat','comp.ind','stst.dat', 'diagn.dat', 'yhea.dat', 'lhea.dat', 'spcp.dat', 'sica.dat', 'eica.dat', 'spct.dat', 'atmd.oud', 'ptfn.dat', 'taul.dat', 'sbib.dat']
    if id==None:
        files=glob('hist/dica.bak????')
        if files==[]: id=1
        else: id=max([int(item[-4:]) for item in files])+1
    for iFile in backupFiles:
        try:
            print "Copying %s to %s"%(iFile,'hist/%s.bak%04d'%(iFile[:-4],id))
            shutil.move(iFile,'hist/%s.bak%04d'%(iFile[:-4],id))
        except:
            print "Had trouble copying %s"%iFile
    return id
def getLastID():
    files=glob('hist/dica.bak????')
    return max([int(item[-4:]) for item in files])
def getHistComp(id):
    compData=fileio.compfile('hist/comp.bak%04d'%id).read_data()
    return param.comp(compData)
def getHistDica(id):
    dicaData=fileio.dicafile('hist/dica.bak%04d'%id).read_data()
    return param.dica(dicaData)
def getHistLList(id):
    sbib=fileio.sbibfile('hist/sbib.bak%04d'%id).read_data()
    return sbib['llist']
def getLastComp():
    id=getLastID()
    return getHistComp(id)
def getLastDica():
    id=getLastID()
    return getHistDica(id)
def getPrevDica():
    dicaFile=os.path.join(config.getLastMainDir(),'manual','dica.dat')
    dica=fileio.dicafile(dicaFile).read_data()
    print "Time=%s"%config.getTimeFromExplosion()
    dica['t']=config.getTimeFromExplosion()
    return dica
def getPrevComp():
    compFile=os.path.join(config.getLastMainDir(),'manual','comp.ind')
    comp=param.comp(fileio.compfile(compFile).read_data())
    return comp
