#parameter classes -- especially for comp
#todo: rewrite locking in comp (mutually excluding ticr and igelocking and ni0 hanndling better)
import math,copy,os
import numpy as np
import config,fileio,initialize,abund
day2sec=3600*24
sigma=5.6705e-8
LSun=3.846e26
joule2ergs=1e7
class comp(object):
    def __init__(self,initComp=None,mode='initW7'):
        self.t=config.getTimeFromExplosion()
        self.autoOxyCorr=True
        self.decayNi=False
        self.lockTiCr=False
        self.TiCrSet=set(('Ti','Cr'))
        self.lockIGEwoNi=True
        self.lockIGEwNi=False
        self.IGEwoNiSet=set(('Sc','Ti','V','Cr','Mn','Cu','Zn'))
        self.IGEwNiSet=self.IGEwoNiSet.union(['Ni0'])
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
        if initComp!=None:
            self.data.update(initComp)
        elif mode=='initW7':
            self.data.update(initialize.getW7Comp(t=self.t))
            self.Fe0=self.data['Fe']
            self.normalize()
        #for key in self.data:
         #   self.__dict__[key]=self[key]
            #setattr(self,key,property(lambda self: self.__getitem__(self,key),lambda self,value: self.__setitem__(self,key,value)))
    def __getitem__(self, key):
        if type(key)==str:
            key=key[0].upper()+key[1:]
            if key in self.data.keys():
                return self.data[key]
        return key
    def __setitem__(self,key,value):
        if type(key)==str:
            key=key[0].upper()+key[1:]
            if key in self.TiCrSet and self.lockTiCr:
                print key
                self._setElementKeepRatio(key,value,*(self.TiCrSet.difference([key])))
            elif key in self.IGEwoNiSet and self.lockIGEwoNi:
                self._setElementKeepRatio(key,value,*(self.IGEwoNiSet.difference([key])))
            elif key in self.IGEwNiSet and self.lockIGEwNi:
                self._setElementKeepRatio(key,value,*(self.IGEwNiSet.difference([key])))
            elif key=='Co':
                print 'Warning: setting non-decaying Co. NOT RECOMMENDED.'
                self._setElement(key,value)
            elif key in self.data.keys():
                self._setElement(key,value)
    def _setNiDecay(self):
        ni,co,fe=abund.calcNiDecay(self.data['Ni0'],self.t)
        self._setElement('Ni',ni)
        self._setElement('Co',co)
        self._setElement('Fe',self.data['Fe0']+fe)
    def _setElement(self,element,abundance):
        if element=='Ni0' and self.decayNi:
            self.data['Ni0']=abundance
            self._setNiDecay()
        if element=='Fe0' and self.decayNi:
            self.data['Fe0']=abundance
            self._setNiDecay()
        if self.autoOxyCorr:
            self.oxyCorrect(element,abundance)
        else:
            self.data['element']=abundance
            self.normalize()
    def _setElementKeepRatio(self,element,abundance,*args):
        ratios=[self.data[item]/self.data[element] for item in args]
        self._setElement(element,abundance)
        for iElement,ratio in zip(args,ratios):
            print "Setting Element %s from %s to %s"%(iElement,self[iElement],ratio*abundance)
            self._setElement(iElement,ratio*abundance)
    def oxyCorrect(self,element,abundance):
        if element[-1]!='0':
            deltaElem=self.data[element]-abundance
            self.data['O']+=deltaElem
            self.data[element]=abundance
            if self.data['O']<0: print "WARNING!!!!: OXYGEN abundance negative"
        else:
            self.data[element]=abundance
    def write2file(self,fileName='comp.ind'):
        fileio.compfile(fileName,'w').write_data(self.data)
    def normalize(self):
        norm=sum([value for key,value in self.data.items() if key[-1]!='0'])
        for key,value in self.data.items():
            self.data[key]=value/norm
    
class dica(object):
    def __init__(self,initDica=None,mode='init'):
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
                    'lg_tau': -2.0,
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
                    'v_ph': -1.0,
                    'wl': 0.25,
                    'z': -1.0,
                    'xe1':0.20}
        self.tRise=19.5
        self._lockTemp=False
        self._lockVph=False
        self._lockLum=False
        self.data.update(config.getInitDicaParam(tRise=self.tRise))
        if initDica!=None:
            self.data.update(initDica)
        else:
            print "Warning: Dica initialized with non-sensible default values. This might be problematic."
    def __getitem__(self, key):
        if type(key)==str:
            if key=='lum': key='log_lbol'
            if key=='vph': key='v_ph'
            
            if key in self.data.keys():
                return self.data[key]
            else: raise KeyError
        else: raise KeyError
    def __setitem__(self, key,value):
        if type(key)==str:
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
        else: self.dica=dica
        if initComp==None: self.comp=comp()
        else: self.comp=comp
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
        else: self.initParam=initParam
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
                print "single"
                tmpParam=copy.deepcopy(self.initParam)
                tmpParam[key]=param
                tmpParam.targetDir='%s_%s'%(key,param)
                tmpList.append(tmpParam)
        else:
            for param in params:
                print "multiple"
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
                print [item.targetDir for item in tmpParam.flatten()]
                tmpList.append(tmpParam)
        self.paramsInGrid.append(key)
            #okay so you dont forget [a,a,a,a,a,a]
            #then [[a,a,a,a,a,a],[a,a,a,a,a,a],[a,a,a,a,a,a],[a,a,a,a,a,a]]
            #then reshape -- simple
                
                    
        self.paramGrid=np.array(tmpList)