from scipy.integrate import trapz
from scipy import polyfit, ndimage

import pdb
import os
import pickle
import copy
import numpy as np
import dalekSteps,launcherSteps,read,param,initialize,fitelem,plot,util
import pprint
import select
import sys
seedno=250819801106
def initLum(doPickle=True):
    runDirs,params=launcherSteps.initLaunchLum()
    lumModelGrid=read.modelGrid(runDirs,params)
    basePath=os.path.split(runDirs[0])[0]
    if doPickle:
        pickle.dump(lumModelGrid,file(os.path.join(basePath,'lum0.pkl'),'w'))
    return lumModelGrid,dalekSteps.getNextLumInterval(lumModelGrid)

def runLumCycle(doPickle=True,maxIter=None,meritThresh=0.01,samples=10,fitHist=None):
    if fitHist==None: fitHist=param.fitHistory()
    lumModelGrid,interval=initLum()
    #fitHist.addHistItem([interval],[lumModelGrid],[None])
    
    #print interval
    i=1
    while True:
        runDirs,params=launcherSteps.launchLum(np.linspace(interval['interval'][0],interval['interval'][1],samples))
        lumModelGrid=read.modelGrid(runDirs,params)
        basePath=os.path.split(runDirs[0])[0]
        if doPickle:
            pickle.dump(lumModelGrid,file(os.path.join(basePath,'lum%d.pkl'%i),'w'))
        i+=1
        interval=dalekSteps.getNextLumInterval(lumModelGrid)
        #fitHist.addHistItem([interval],[lumModelGrid],[None])
        print "Current interval:",interval
        if maxIter!=None:
            if i>maxIter: break
        #fitHist.write2pickle('lumCycle.pkl')
        if interval['merit']<meritThresh: break
        
    return interval

def initTriCycle(IGEElement,doPickle=True,samples=10):
    fitHist=param.fitHistory()
    lumInterval=runLumCycle(maxIter=5)
    curParam=param.param()
    curParam['lum']=lumInterval['suggestValue']
    intervals={'luminterval':lumInterval['interval'],
               'vphinterval':initialize.getVphBounds(),
                #make sure that no other metal really shoots up
                'igeinterval':initialize.getElementBounds(IGEElement,curParam.comp)}
    #pdb.set_trace()
    lumRange=np.linspace(intervals['luminterval'][0],intervals['luminterval'][1],num=samples)
    vphRange=np.linspace(intervals['vphinterval'][0],intervals['vphinterval'][1],num=samples)
    IGERange=np.linspace(intervals['igeinterval'][0],intervals['igeinterval'][1],num=samples)
    lumMG,vphMG,IGEMG,curParamMG=launcherSteps.launchTriCycle(lumRange,vphRange,IGERange,IGEElement,initParam=curParam)
    if doPickle:
        pickle.dump(lumMG,file('lumMG0.pkl','w'))
        pickle.dump(vphMG,file('vphMG0.pkl','w'))
        pickle.dump(IGEMG,file('igeMG0.pkl','w'))
    
    evalTriCycle(lumMG,vphMG,IGEMG,curParamMG,curParam,intervals,fitHist,IGEElement,mode='init')
    return fitHist,curParam,intervals
def runTriCycle(param=None,samples=10,fitHist=None,doPickle=True):
    IGEElement='Fe0'
    fitHist,curParam,intervals=initTriCycle(IGEElement)
    i=1
    while True:
        print "starting cycle %s"%i
        lumRange=np.linspace(intervals['luminterval'][0],intervals['luminterval'][1],num=samples)
        vphRange=np.linspace(intervals['vphinterval'][0],intervals['vphinterval'][1],num=samples)
#        IGERange=np.linspace(intervals['igeinterval'][0],intervals['igeinterval'][1],num=samples)
        IGERange=np.logspace(np.log10(intervals['igeinterval'][0]),np.log10(intervals['igeinterval'][1]),num=samples)
        lumMG,vphMG,IGEMG,curParamMG=launcherSteps.launchTriCycle(lumRange,vphRange,IGERange,IGEElement,initParam=copy.deepcopy(curParam),procPath='tri%02d'%i)
        print IGERange
        print "EVALING Cycle %s"%i
        print
        done=evalTriCycle(lumMG,vphMG,IGEMG,curParamMG,curParam,intervals,fitHist,IGEElement)
        #pdb.set_trace()
        #if done: break
        i+=1
    pickle.dump(curParam,file('curTriParam.pkl','w'))
    return curParam

def evalTriCycle(lumMG,vphMG,IGEMG,curParamMG,curParam,intervals,fitHist,IGEElement,mode='std'):
    #IMPORTANT!!!!!!!! the function DOES NOT need curParamMG, only for adding it to the fithist object
    #REMOVE AS SOON AS DEBUGGING for
    pp=pprint.PrettyPrinter()
    selectSteadyThresh=6
    steadyThresh=4
    pp.pprint(intervals)
    pivots=[0,0,0]
    lumBounds=initialize.getLumBounds()
    vphBounds=initialize.getVphBounds()
    IGEBounds=initialize.getElementBounds(IGEElement,curParam.comp)
    lumInterval=dalekSteps.getNextLumInterval(lumMG)
    vphInterval=dalekSteps.getNextVphInterval(vphMG)
    IGEInterval=dalekSteps.getNextIGEInterval(IGEMG,IGEElement)
    
    lastSuggestLum=np.mean(intervals['luminterval'])
    lastSuggestVph=np.mean(intervals['vphinterval'])
    lastSuggestIGE=np.mean(intervals['igeinterval'])
    print "CurLum=%.4f CurVph=%.2f CurFe0=%.7f"%(curParam['lum'],curParam['vph'],curParam['Fe0'])
    print "---------------------------------"
    
    print "LUM:"
    pp.pprint(lumInterval)
    print
    print "VPH:"
    
    pp.pprint(vphInterval)
    print
    print "IGE:"
    pp.pprint(IGEInterval)
    
    #Updating the intervals if the solution is not steady
    steadyIncFactor=2
    selectSteadyIncFactor=2
    
    curLum=curParam['lum']
    closestLumID=np.argmin(lumMG['lum']-curLum)
    
        
    #intervals=[lumInterval,vphInterval,IGEInterval]
    merits=np.argsort([item[2] for item in intervals])
    #pivots[merits[2]]+=3
    #pivots[merits[1]]+=1
    pivots[0]=3
    pivots[1]=2
    pivots[2]=1
    #making sure that we are close in luminosity
    print IGEMG['Fe0']
    if not (abs((lumMG['lum'][closestLumID]-curLum))<0.05 and abs(lumInterval['merits'][closestLumID])<0.1):
        print "lum doesnt allow other parameters atm"
        print "closestLum %s closestMerit %s"%(lumMG['lum'][closestLumID]-curLum,lumInterval['merits'][closestLumID])
        pivots[1]-3
        pivots[2]-3
        pivots[0]=10
    #Checking if the interval is too small already
    if lumInterval['mininterval']: pivots[0]-=5
    if vphInterval['mininterval']: pivots[1]-=5
    #IGE NOT WORKING properly
    #pivots[2]=-10
    print "PIVOTS %s"%pivots
    
    #Checking w factors:
    print 'vph w %s'%vphMG['w']
    print 'ige w %s'%IGEMG['w']
    if np.sum(np.abs(vphMG['w']-0.5)>0.1)>4: pivots[1]-=1.5
    if np.sum(np.abs(IGEMG['w']-0.5)>0.1)>4: pivots[2]-=1
    
    damping=True
    #the factor does not do that much atm BEWARE
    dampingFactor=2
    if np.argmax(pivots)==0:
        #param=copy.copy(lumMG.grid[lumInterval['bestFitID']].param)
        if lumInterval['steady']>selectSteadyThresh:
            if damping:
                curParam['lum']=np.mean([lumInterval['suggestValue'],lastSuggestLum])
                intervals.update({'luminterval':[curParam['lum']+i*dampingFactor*lumInterval['dev'] for i in [-1,1]]})
            else:
                curParam['lum']=lumInterval['suggestValue']
                intervals.update({'luminterval':lumInterval['interval']})
        else:
            print "Not steady solution will increase last interval and try again"
            dev=selectSteadyIncFactor*np.abs(intervals['luminterval'][0]-curParam['lum'])
            intervals.update({'luminterval':[curParam['lum']-dev,curParam['lum']+dev]})
        print "Chose Lum: updated for next run with suggestValue %s and Interval %s"%(curParam['lum'],intervals['luminterval'])
    else:
        if lumInterval['steady']<steadyThresh:
            print "Updating the lum interval as the current solution is not steady"
            dev=steadyIncFactor*np.abs(intervals['luminterval'][0]-curParam['lum'])
            intervals.update({'luminterval':[curParam['lum']-dev,curParam['lum']+dev]})

    if np.argmax(pivots)==1:
        #param=copy.copy(vphMG.grid[vphInterval['bestFitID']].param)
        if vphInterval['steady']>selectSteadyThresh:
            if damping:
                curParam['vph']=np.mean([vphInterval['suggestValue'],lastSuggestVph])
                intervals.update({'vphinterval':[curParam['vph']+i*dampingFactor*vphInterval['dev'] for i in [-1,1]]})
            else:
                curParam['vph']=vphInterval['suggestValue']
                intervals.update({'vphinterval':vphInterval['interval']})
        else:
            print "Not steady solution will increase last interval and try again"
            dev=selectSteadyIncFactor*np.abs(intervals['vphinterval'][0]-curParam['vph'])
            intervals.update({'vphinterval':[curParam['vph']-dev,curParam['vph']+dev]})
        print "Chose Vph: updated for next run with suggestValue %s and Interval %s"%(curParam['vph'],intervals['vphinterval'])
    else:
        if vphInterval['steady']<steadyThresh:
            print "Updating the vph interval as the current solution is not steady"
            dev=steadyIncFactor*np.abs(intervals['vphinterval'][0]-curParam['vph'])
            intervals.update({'vphinterval':[curParam['vph']-dev,curParam['vph']+dev]})
        

    if np.argmax(pivots)==2:
        if IGEInterval['steady']>selectSteadyThresh:
            if damping:
                curParam[IGEElement]=np.mean([IGEInterval['suggestValue'],lastSuggestIGE])
                intervals.update({'igeinterval':[curParam[IGEElement]+i*dampingFactor*IGEInterval['dev'] for i in [-1,1]]})
            else:
                curParam[IGEElement]=IGEInterval['suggestValue']
                intervals.update({'igeinterval':IGEInterval['interval']})
        else:
            print "Not steady solution will increase last interval and try again"
            dev=selectSteadyIncFactor*np.abs(intervals['igeinterval'][0]-curParam[IGEElement])
            intervals.update({'igeinterval':[curParam[IGEElement]-dev,curParam[IGEElement]+dev]})
        print "Chose IGE: updated for next run with suggestValue %s and Interval %s"%(curParam[IGEElement],intervals['igeinterval'])
    else:
            if IGEInterval['steady']<steadyThresh:
                print "Updating the IGE interval as the current solution is not steady"
                dev=steadyIncFactor*np.abs(intervals['igeinterval'][0]-curParam[IGEElement])
                intervals.update({'igeinterval':[curParam[IGEElement]-dev,curParam[IGEElement]+dev]})


    #Checking for crossing Boundaries:
    boundCheck=zip(lumBounds,intervals['luminterval'])
    lumCheck=[np.max(boundCheck[0]),np.min(boundCheck[1])]
    boundCheck=zip(vphBounds,intervals['vphinterval'])
    vphCheck=[np.max(boundCheck[0]),np.min(boundCheck[1])]
    boundCheck=zip(IGEBounds,intervals['igeinterval'])
    igeCheck=[np.max(boundCheck[0]),np.min(boundCheck[1])]
    intervals.update({'luminterval':lumCheck,'vphinterval':vphCheck,'igeinterval':igeCheck})


    #Readjusting suggest values to be in the middle of the intervals:
    curParam['lum']=np.mean(intervals['luminterval'])
    curParam['vph']=np.mean(intervals['vphinterval'])
    curParam[IGEElement]=np.mean(intervals['igeinterval'])
    #Printing the intervals:
    pp.pprint(intervals)

    #Adding to fitHist
    fitHist.addHistItem([lumInterval,vphInterval,IGEInterval],[lumMG,vphMG,IGEMG,curParamMG],pivots,intervals,reallyAdd=False,saveSingle=True)
    #pdb.set_trace()
    #Break condition
    if np.argmax(pivots)>0:
    #    pdb.set_trace()
        return True
        
    else:
        return False
def initElementCycle(initParam,elements,samples,scaling):
    bounds=initialize.getElementsBounds(initParam.comp)
    divIntervals={'Si':[5e-3,5],
                'S':[1e-3,1],
                'Ca':[0.5e-3,0.5],
                'Mg':[0.5e-3,2],
                'C':[0.5e-3,0.5],
                'Fe0':[1e-3,1]}
    
    intervals=copy.deepcopy(bounds)
    initParam.lockIGE=False
    initParam.lockIGEwNi=False
    initParam.lockIGEwoNi=True
    ranges={}
    #intervalHistory={}
    for element in intervals.keys():
        if scaling=='log':
            ranges[element]=np.logspace(*np.log10(intervals[element]),**dict(num=samples))
            
        elif scaling=='lin':
            ranges[element]=np.linspace(*intervals[element],**dict(num=samples))
        else:
            raise Exception('Unknown scaling option %s. Please use lin or log.'%scaling)
    finished=False
    for element in ['Si','S','Mg','C']:
        if scaling=='log':
            elementRange=np.logspace(*np.log10(intervals[element]),**dict(num=samples))
        elif scaling=='lin':
            elementRange=np.linspace(*intervals[element],**dict(num=samples))
        
        elementMG=launcherSteps.launchElementCycle({element:elementRange},initParam=initParam)
        plot.plotElementModelGrid(element,elementMG[0])    
        newInterval,finished=dalekSteps.getNextElementParams(element,elementMG[0],initParam,bounds[element])
        if not finished: intervals[element]=newInterval
        print "New Interval for %s: %s "%(element,intervals[element])
        initParam[element]=np.mean(intervals[element])
        bounds=initialize.getElementsBounds(initParam.comp)
        intervals=dalekSteps.checkElementBounds(intervals,bounds)
        #pdb.set_trace()
    #allMG=launcherSteps.launchElementCycle(ranges,initParam=initParam)
    #elementMG=dict(zip(ranges.keys(),allMG[:-1]))%debug
    #curParamMG=allMG[-1]
    
    return intervals,initParam
    
def runElementCycle(initParam,elements=['C','Ca','Mg','Si','S','Ti','Ni0','Fe0'],samples=10,scaling='lin'):
    #initParam.autoRelAbund='O'
    intervals,initParam=initElementCycle(initParam,elements,samples,scaling)
    #print "Dumping initial run in initElemCycle.pkl"
    #pickle.dump(elementMG,file('initElemCycle.pkl','w'))
    bounds=initialize.getElementsBounds(initParam.comp)
    intervalHistory={'Si':[],'S':[],'Mg':[],'C':[]}
    i=1
    while True:
        converged=[]
        for element in ['Si','S','Mg','C']:
            if scaling=='log':
                elementRange=np.logspace(*np.log10(intervals[element]),**dict(num=samples))
            elif scaling=='lin':
                elementRange=np.linspace(*intervals[element],**dict(num=samples))
            
            elementMG=launcherSteps.launchElementCycle({element:elementRange},initParam=initParam)
            plot.plotElementModelGrid(element,elementMG[0],suffix=str(i))    
            intervals[element],finished=dalekSteps.getNextElementParams(element,elementMG[0],initParam,bounds[element])
            converged.append(finished)
            intervalHistory[element].append(intervals[element])
            print "New Interval for %s: %s "%(element,intervals[element])
            initParam[element]=np.mean(intervals[element])
            bounds=initialize.getElementsBounds(initParam.comp)
            intervals=dalekSteps.checkElementBounds(intervals,bounds)
        print "Convergence status %s"%converged
        if all(converged): break
        print "Looking for stop.txt in %s"%os.getcwd()
        if os.path.exists('stop.txt'): break
        i+=1
    return intervalHistory