#parameter classes -- especially for comp
#todo: rewrite locking in comp (mutually excluding ticr and igelocking and ni0 hanndling better)
import math,copy,os
import pickle
import numpy as np
import config,fileio,initialize,abund
#import pylab
import pdb
day2sec=3600*24
sigma=5.6705e-8
LSun=3.846e26
joule2ergs=1e7
class comp(object):
    def __init__(self,initComp=None, t=None, mode='initW7'):
        if t == None: self.t=config.getTimeFromExplosion()
        else: self.t = t
        self.autoOxyCorr=True
        self.decayNi=True
        self.lockTiCr=False
        self.lockScTi=False
        self.lockVCr=False
        self.TiCrSet=set(('Ti','Cr'))
        self.ScTiSet=set(('Sc','Ti'))
        self.lockIGEwoNi=False
        self.lockIGEwNi=False
        self.lockIGE=False
        self.oxWarn=True
        self.autoRelAbund=None
        self.keepElementRel=False
        
        self.VCrSet=set(('V','Cr'))
        self.IGEwoNiSet=set(('Sc','Ti','V','Cr','Mn','Cu','Zn'))
        self.IGEwNiSet=self.IGEwoNiSet.union(['Ni0'])
        self.IGESet=self.IGEwNiSet.union(['Fe0'])
        self.data={'Al': 0.0,
                   'Ar': 0.0,
                   'B': 0.0,
                   'Be': 0.0,
                   'C': 0.0,
                   'Ca': 0.0, 
                   'Cl': 0.0,
                   'Co': 0.0,
                   'Cr': 0.0,
                   'Cu': 0.0,
                   'F': 0.0,
                   'Fe': 0.0,
                   'H': 0.0,
                   'He': 0.0,
                   'K': 0.0,
                   'Li': 0.0,
                   'Mg': 0.0,
                   'Mn': 0.0,
                   'N': 0.0,
                   'Na': 0.0,
                   'Ne': 0.0,
                   'Ni': 0.0,
                   'O': 0.0,
                   'P': 0.0,
                   'S': 0.0,
                   'Sc': 0.0,
                   'Si': 0.0,
                   'Ti': 0.0,
                   'V': 0.0,
                   'Zn': 0.0,
                   'Ni0':0.0,
                   'Fe0':0.0}
        self.FeDecayed=0.0
        if initComp!=None:
            self.data.update(initComp)
            if self.data['Ni0']==0.0 and self.data['Fe0']==0.0:
                print "Calculating Ni0 and Fe0"
                curNi,curCo,curFe=self.data['Ni'],self.data['Co'],self.data['Fe']
                curSum=np.sum((curNi,curCo,curFe))
                decNi,decCo,decFe=abund.calcNiDecay(curSum,self.t)
                factorNi=curSum/decNi
                self['Ni0']=curNi*factorNi
                self['Fe0']=curFe-self['Fe']
                if self['Fe0']<0: raise Exception('Problems while calculating Ni0 and Fe0')
                
        elif mode=='initW7':
            self.data.update(initialize.getW7Comp(t=self.t))
            self.Fe0=self.data['Fe']
            curAutoOxyCorr=self.autoOxyCorr
            self.autoOxyCorr=False
            self._setNiDecay()
            self.resetOxygen()
            self.autoOxyCorr=curAutoOxyCorr
            #self.normalize()
    def __getitem__(self, key):
        if isinstance(key,str):
            key=key[0].upper()+key[1:]
            if key in self.data.keys():
                return self.data[key]
        return key
    def __setitem__(self,key,value):
        if isinstance(key,str):
            key=key[0].upper()+key[1:]
            if key in self.TiCrSet and self.lockTiCr:
                self._setElementKeepRatio(key,value,*(self.TiCrSet.difference([key])))
            if key in self.ScTiSet and self.lockScTi:
                self._setElementKeepRatio(key,value,*(self.ScTiSet.difference([key])))
            if key in self.VCrSet and self.lockVCr:
                self._setElementKeepRatio(key,value,*(self.VCrSet.difference([key])))
                
            elif key in self.IGEwoNiSet and self.lockIGEwoNi:
                self._setElementKeepRatio(key,value,*(self.IGEwoNiSet.difference([key])))
            elif key in self.IGEwNiSet and self.lockIGEwNi:
                self._setElementKeepRatio(key,value,*(self.IGEwNiSet.difference([key])))
            elif key in self.IGESet and self.lockIGE:
                self._setElementKeepRatio(key,value,*(self.IGESet.difference([key])))
            elif key=='Co':
                print 'Warning: setting non-decaying Co. NOT RECOMMENDED.'
                self._setElement(key,value)
            elif key in self.data.keys():
                self._setElement(key,value)
    def _setNiDecay(self):
        curOxWarn=self.oxWarn
        self.oxWarn=False
        ni,co,fe=abund.calcNiDecay(self.data['Ni0'],self.t)
        self.FeDecayed=fe
        self._setElement('Ni',ni)
        self._setElement('Co',co)
        self.oxWarn=curOxWarn
        self._setElement('Fe',self.data['Fe0']+fe)
    def _setElement(self,element,abundance):
        if element=='Ni0' and self.decayNi:
            self.data['Ni0']=abundance
            self._setNiDecay()
        if element=='Fe0' and self.decayNi:
            self.data['Fe0']=abundance
            self._setElement('Fe',self.data['Fe0']+self.FeDecayed)
        if self.autoRelAbund!=None:
            relElement=self.autoRelAbund
            if self.keepElementRel:
                relAbundances=dict([(item[0],item[1]/self.data[relElement]) for item in self.data.items() if not (item[0]==element or item[0][-1]=='0')])
                relAbundances[element]=abundance
                newRelElementAbund=1/np.sum(relAbundances.values())
                for item in relAbundances.items(): self.data[item[0]]=item[1]*newRelElementAbund
            else:
                while True:
                    newRelElementAbund=1/(np.sum([item[1]/self.data[relElement] for item in self.data.items() if not (item[0]==element or item[0][-1]=='0')])+abundance)
                    self.data[element]=newRelElementAbund*abundance
                    self.data[relElement]=newRelElementAbund
                    if np.abs(self.getNorm()-1)<1e-10: break
        elif self.autoOxyCorr:
            self.oxyCorrect(element,abundance)
        else:
            self.data[element]=abundance
            self.normalize()
    def _setElementKeepRatio(self,element,abundance,*args):
        ratios=[self.data[item]/self.data[element] for item in args]
        self._setElement(element,abundance)
        for iElement,ratio in zip(args,ratios):
            #print "Setting Element %s from %s to %s"%(iElement,self[iElement],ratio*abundance)
            self._setElement(iElement,ratio*abundance)
    def oxyCorrect(self,element,abundance):
        if element[-1]!='0':
            deltaElem=self.data[element]-abundance
            self.data['O']+=deltaElem
            self.data[element]=abundance
            if self.data['O']<0 and self.oxWarn:
                print "WARNING!!!!: OXYGEN abundance negative %s (setting of Element %s=%s did this)"%(self.data['O'],element,abundance)
        else:
            self.data[element]=abundance
    def write2file(self,fileName='comp.ind'):
        fileio.compfile(fileName,'w').write_data(self.data)
    def getNorm(self):
        return sum([value for key,value in self.data.items() if key[-1]!='0'])
    def resetOxygen(self):
        sumElementsWOOxygen=sum([value for key,value in self.data.items() if (key[-1]!='0' and key!='O')])
        #pdb.set_trace()
        self.data['O']=1-sumElementsWOOxygen
        if self.data['O']<0:
            print "WARNING!!!!: OXYGEN abundance negative. The from pyfica import iinitial parameters were not normalized"
    def normalize(self):
        norm=self.getNorm()
        for key,value in self.data.items():
            self.data[key]=value/norm

#def intervalComp(comp):
#    def __init__(self)
        #super(self).__init__
class dica(object):
    def __init__(self,initDica=None,mode='init'):
        if mode=='init':
            self.data={ 'chl': 25.0,
                        'e_b-v_gal': -1.0,
                        'e_b-v_host': -1.0,
                        'em_high': 6500.0,
                        'em_low': 1000.0,
                        'grid': 1.25,
                        'inthigh': 10000.0,
                        'intlow': 2500.0,
                        'itt': 4.0,
                        'js': 20.0,
                        'kb': 501.0,
                        'kr': 250819801106,
                        'lg_tau': -3.5,
                        'log_l_low_high': 8.9399999999999995,
                        'log_lbol': -1.0,
                        'm-m': -1.0,
                        'mb': 10000.0,
                        'mu': 1500.0,
                        'nc': 5.0,
                        'np5': 4.0,
                        'options': '1   1   1   1   1',
                        't': -1.0,
                        'tb': 10000.0,
                        'v_ph': initialize.getCurVph(),
                        'wl': 0.25,
                        'z': -1.0,
                        'xe1':0.20}
            self.tRise=19.5
            self._lockTemp=False
            self._lockVph=False
            self._lockLum=False
            self.data.update(config.getInitDicaParam(tRise=self.tRise))
        elif mode!='init':
            self.data=dict()
        if initDica!=None:
            self.data.update(initDica)
        else:
            pass
            #pdb.set_trace()
            #print "Warning: Dica initialized with non-sensible default values. This might be problematic."
    def __getitem__(self, key):
        if isinstance(key,str):
            key=key.lower()
            if key=='lum': key='log_lbol'
            if key=='vph': key='v_ph'
            
            if key in self.data.keys():
                return self.data[key]
            else: raise KeyError
        else: raise KeyError
    def __setitem__(self, key,value):
        if isinstance(key,str):
            key=key.lower()
            if key=='log_lbol' or key=='lum':
                if self.lockTemp:
                    T=self.getTemp()
                    self.data['log_lbol']=value
                    self.data['v_ph']=self._TLum2Vph(T,value)
                else:
                    self.data['log_lbol']=value    
            elif key=='v_ph' or key=='vph':
                if self.lockTemp:
                    T=self.getTemp()
                    self.data['v_ph']=value
                    self.data['log_lbol']=self._TVph2Lum(T,value)
                else:
                    self.data['v_ph']=value
            elif key in self.data.keys():
                self.data[key]=value
            else: raise KeyError
        else: raise KeyError
    def getLBol(self):
        return self._logLBol2LBol(self.data['log_lbol'])
    def _logLBol2LBol(self,logLBol):
        return (10**(logLBol))*LSun
    def getLockTemp(self):
        return self._lockTemp
    def setLockTemp(self,value):
        self._lockTemp=value
        if value==True:
            self._lockLum=False
            self._lockVph=False
    def getLockLum(self):
        return self._lockLum
    def setLockLum(self,value):
        self._lockLum=value
        if value==True:
            self._lockT=False
            self._lockVph=False
    def getLockVph(self):
        return self._lockVph
    def setLockVph(self,value):
        self._lockVph=value
        if value==True:
            self._lockT=False
            self._lockLum=False
    def getTemp(self):
        return self._lumVph2T(self.data['log_lbol'],self.data['v_ph'])
    def setTemp(self,temp):
        if self.lockLum:
            self.data['v_ph']=self._TLum2Vph(temp,self.data['log_lbol'])
        elif self.lockVph:
            self.data['log_lbol']=self._TVph2Lum(temp,self.data['v_ph'])
    def _lumVph2T(self,lum,vph):
        t=self.data['t']*day2sec
        r=vph*1e3*t
        return (self._logLBol2LBol(lum)/(4*math.pi*sigma*(r**2)))**0.25
        
    def _TLum2Vph(self,temp,lum):
        t=self.data['t']*day2sec
        return ((self._logLBol2LBol(lum)/(4*math.pi*sigma*(temp**4)*(t**2) ))**0.5)*1e-3
        
    def _TVph2Lum(self,temp,vph):
        t=self.data['t']*day2sec
        r=vph*1e3*t
        LSn=4*math.pi*sigma*(r**2)*(temp**4)
        return math.log(LSn/LSun, 10)
    def getLocks(self):
        print "lockLum %s"%self._lockLum
        print "lockVph %s"%self._lockVph
        print "lockTemp %s"%self._lockTemp
    lBol=property(getLBol)
    T=property(getTemp,setTemp)
    lockLum=property(getLockLum,setLockLum)
    lockVph=property(getLockVph,setLockVph)
    lockTemp=property(getLockTemp,setLockTemp)
    def write2file(self,fileName='dica.dat'):
        fileio.dicafile(fileName,'w').write_data(self.data)

class param(object):
    def __init__(self,initDica=None,initComp=None,targetDir=''):
        if initDica==None: self.dica=dica()
        else: self.dica=initDica
        if initComp==None: self.comp=comp()
        else: self.comp=initComp
        self.targetDir=targetDir
    def __getitem__(self,key):
        if key.lower()=='targetdir':
            return self.targetDir
        else:
            try: return self.dica[key]
            except KeyError: return self.comp[key]
    def __setitem__(self,key,value):
        try: self.dica[key]=value
        except KeyError: self.comp[key]=value
    def write2file(self,baseDir='.'):
        self.dica.write2file(os.path.join(baseDir,'dica.dat'))
        self.comp.write2file(os.path.join(baseDir,'comp.ind'))
    
        
        
        
        
class multiParam(object):
    def __init__(self,initParam=None):
        if initParam==None: self.initParam=param()
        else: self.initParam=copy.deepcopy(initParam)
        self.paramGrid=np.array([])
        self.paramsInGrid=[]
    def __setitem__(self,key,value):
        self._setParamSet(key,value)
    def __getitem__(self,key):
        getFunc=np.vectorize(lambda item:item[key])
        return getFunc(self.paramGrid)
    def _setParamSet(self,key,params):
        tmpList=[]
        if not np.iterable(params): raise Exception('Needs a sequence of items to work on')
        if self.paramGrid.size==0:
            for param in params:
                tmpParam=copy.deepcopy(self.initParam)
                tmpParam[key]=param
                tmpParam.targetDir='%s_%s'%(key,param)
                tmpList.append(tmpParam)
        else:
            for param in params:
                def setItem(item):
                    item[key]=param
                    item.targetDir+='%s_%s'%(key,param)
                    return None
                setterFunc=np.vectorize(setItem)
                tmpParam=copy.deepcopy(self.paramGrid)
                #print [item.targetDir for item in tmpParam.flatten()]
                #print tmpParam
                for item in tmpParam.flatten(): setItem(item)
                #setterFunc(tmpParam.flatten())
                tmpList.append(tmpParam)
        self.paramsInGrid.append(key)
            #okay so you dont forget [a,a,a,a,a,a]
            #then [[a,a,a,a,a,a],[a,a,a,a,a,a],[a,a,a,a,a,a],[a,a,a,a,a,a]]
            #then reshape -- simple
                
                    
        self.paramGrid=np.array(tmpList)

class fitHistory(object):
    def __init__(self):
        self.curStep=0
        self.evalDics=[]
        self.modelGrids=[]
        self.evalPivots=[]
        self.paramHist=[]
        self.specPdf=PdfPages('specFit.pdf')
        self.paramPdf=PdfPages('paramFit.pdf')
        self.writePickle=False
        self.writeSpecPdf=True
        self.writeParamPdf=True
        #self.modelPerStep=3
    def __getstate__(self):
        return
    def __del__(self):
        print "fitHist is being deleted"
        if self.writePickle:
            self.write2pickle('fitHist.pkl')
        if self.writeSpecPdf:
            self.specPdf.close()
        if self.writeParamPdf:
            self.paramPdf.close()
    def __getitem__(self,key):
        return self.modelGrids[key]
    def addHistItem(self,evalDics,modelGrids,pivots,intervals,reallyAdd=False,saveSingle=True,makeSpecPdf=True,makeParamPdf=True):
        self.curStep+=1
        self.evalDics.append(copy.copy(evalDics))
        self.evalPivots.append(copy.copy(pivots))
        if reallyAdd:
            self.modelGrids.append(copy.copy(modelGrids))

        else:
            print "Nothing added to fitHistory"
        if saveSingle:
            modelGridFileName=('modelgrids%02d.pkl'%self.curStep)
            print "Saving Current Model %s"%modelGridFileName
            pickle.dump(modelGrids,file(modelGridFileName,'w'))
        #self.plot
        #self.fitParams.append()
        if makeSpecPdf:
            specPdfName='specFit%02d.pdf'%self.curStep
            print "Saving Plot to %s"%specPdfName
            pdf=PdfPages(specPdfName)
            self.plotSpecItem(modelGrids,evalDics,pivots,pdfOut=pdf)
            pdf.close()
        if makeParamPdf:
            self.paramHist.append([self.curStep,copy.copy(intervals)])
            self.plotParamItem(copy.copy(pivots))
            #pdb.set_trace()
                
    def plotHistory(self,output,plotBestIDs=[0,1,2,-1,-2,-3]):
        pdf = PdfPages(output)
        for modelGrids,evalDics,evalPivots in zip(self.modelGrids,self.evalDics,self.evalPivots):
            fig=pylab.figure(1)
            fig.clf()
            modelGridsNo=len(evalDics)
            for i,(modelGrid,evalDic,evalPivot) in enumerate(zip(modelGrids,evalDics,evalPivots)):
                #ax=fig.add_subplot(modelGridsNo,1,i+1)
                ax=fig.add_subplot(111)
                ax.plot(modelGrid.origSpec.x,modelGrid.origSpec.y,'k')
                for bestID in plotBestIDs:
                    aSpec=modelGrid['aspec'][evalDic['sortedModelIDX'][bestID]]
                    fitParam=modelGrid[evalDic['fitKey']][evalDic['sortedModelIDX'][bestID]]
                    ax.plot(aSpec.x,aSpec.y,label='%s=%.3f fit %d merit %.5f'%(evalDic['fitKey'],fitParam,evalDic['sortedModelIDX'][bestID],evalDic['merits'][evalDic['sortedModelIDX'][bestID]]))
                ax.legend()
                ax.set_title('%s interval: %.5f-%.5f\n curLum=%s,curVph=%s,curTi=%s'%(evalDic['fitKey'],modelGrid[evalDic['fitKey']][0],modelGrid[evalDic['fitKey']][-1],modelGrid['lum'][0],modelGrid['vph'][0],modelGrid['ti'][0]))
                ax.set_xlim(modelGrid.origSpec.x.min(),modelGrid.origSpec.x.max())
                fig.savefig(pdf,format='pdf')
                fig.clf()
                ax=fig.add_subplot(111)
                param=modelGrid[evalDic['fitKey']]
                merits=evalDic['merits']
                ax.plot(param,merits)
                fig.savefig(pdf,format='pdf')
                fig.clf()
        pdf.close() 
    def plotSpecItem(self,modelGrids,evalDics,evalPivots,pdfOut=None,plotBestIDs=[0,1,2,-1,-2,-3]):
        if pdfOut==None: pdfOut=self.specPdf
        mapPivot2Param=['lum','vph','Fe0']
        fig=pylab.figure(1)
        fig.clf()
        for i,(modelGrid,evalDic,evalPivot) in enumerate(zip(modelGrids,evalDics,evalPivots)):
                #ax=fig.add_subplot(modelGridsNo,1,i+1)
                np.argmax(evalPivots)
                if mapPivot2Param[np.argmax(evalPivots)].lower()==evalDic['fitKey'].lower():
                    ax=fig.add_subplot(111,axisbg='pink')
                else:
                    ax=fig.add_subplot(111)
                #pdb.set_trace()
                ax.plot(modelGrid.origSpec.x,modelGrid.origSpec.y,'k')
                for bestID in plotBestIDs:
                    aSpec=modelGrid['aspec'][evalDic['sortedModelIDX'][bestID]]
                    fitParam=modelGrid[evalDic['fitKey']][evalDic['sortedModelIDX'][bestID]]
                    ax.plot(aSpec.x,aSpec.y,label='%s=%.3f fit %d merit %.5f'%(evalDic['fitKey'],fitParam,evalDic['sortedModelIDX'][bestID],evalDic['merits'][evalDic['sortedModelIDX'][bestID]]))
                ax.legend()
                ax.set_title('%s interval: %.5f-%.5f pivot=%s\n curLum=%s,curVph=%s,curTi=%s'%(evalDic['fitKey'],modelGrid[evalDic['fitKey']][0],modelGrid[evalDic['fitKey']][-1],evalPivot,modelGrid['lum'][0],modelGrid['vph'][0],modelGrid['ti'][0]))
                ax.set_xlim(modelGrid.origSpec.x.min(),modelGrid.origSpec.x.max())
                fig.savefig(pdfOut,format='pdf')
                fig.clf()
                ax=fig.add_subplot(111)
                param=modelGrid[evalDic['fitKey']]
                merits=evalDic['merits']
                ax.plot(param,merits)
                fig.savefig(pdfOut,format='pdf')
                fig.clf()
        
        ax=fig.add_subplot(111,axisbg='lightgreen')
        
        ax.plot(modelGrid.origSpec.x,modelGrid.origSpec.y,'k')
        aSpec=modelGrids[-1]['aspec'][0]
        ax.plot(aSpec.x,aSpec.y,'r')
        ax.set_xlim(modelGrids[-1].origSpec.x.min(),modelGrids[-1].origSpec.x.max())
        fig.savefig(pdfOut,format='pdf')
    
    def plotParamItem(self,pivots,pdfName='paramHist.pdf'):
        pdf = PdfPages(pdfName)
        #pdb.set_trace()
        mapPivot2Param=['luminterval','vphinterval','igeinterval']
        for key in  self.paramHist[0][1].keys():
            evalIDx=mapPivot2Param.index(key)
            if evalIDx==np.argmax(pivots):
                color='green'
            else:
                color='blue'
            fig=pylab.figure(1)
            fig.clf()
            ax=fig.add_subplot(111)
            x=np.array([item[0] for item in self.paramHist])
            yrange=[item[1][key] for item in self.paramHist]
            y=[np.mean(item) for item in yrange]
            yerr=[np.abs(item[0]-item[1])/2. for item in yrange]
            ySuggest=[item[evalIDx]['suggestValue'] for item in self.evalDics]
            ySuggestErr=[item[evalIDx]['dev'] for item in self.evalDics]
            ax.errorbar(x,y,yerr,marker='x',color=color)
            
            ax.errorbar(x+0.5,ySuggest,ySuggestErr,marker='x',color='r')
            if key=='igeinterval': ax.set_ylim(0,1)
            fig.savefig(pdf,format='pdf')
        #pdb.set_trace()
        
        pdf.close()
                
    def write2pickle(self,fname):
        pickle.dump(self,file(fname,'w'))
        
#class triCycleHistory(fitHistory)
    