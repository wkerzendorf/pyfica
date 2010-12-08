#Used to read Models
from pyspec.spectrum import spectrum
#from .fileio import dicafile,compfile,sbibfile
import param as paramMod
import fileio
import config
import fit
from glob import glob
import os,sys,re
import numpy as np
def readModels(modelDir='.', tracerName='dica.dat',backup=False, suffix='.bak'):
    print "Reading in Model from %s"%os.path.join(modelDir, tracerName+'*')
    files=glob(os.path.join(modelDir, tracerName+'*'))    
    model=[0]*len(files)
    for i, ifile in enumerate(files):
        id=int(re.search('\d+',ifile).group())
        filePath=os.path.dirname(ifile)
        print "Reading file set %04d"%id
        print "Progress %s/%s"%(i+1, len(files))
        sys.stdout.write('\x1b[2A')
        model[id]=readModel(id,modelDir=filePath)
#   model[0].update({'origspect':origspect})
    sys.stdout.write('\x1b[2B')
    if backup: 
        print "Backing up model to model.pkl"
        pickle.dump(model,file(path.join(modelDir,'model.pkl'), 'w'))
    return model
def createFakeSpec(dica):
    wl=dica['wl']*1e4,dica['grid']*1e4
    return spectrum(zip(wl,[1,1]))    

def readModel(id,modelDir='.',suffix='.bak',fmt='%04d'):
    suffix=suffix+fmt%id
    dica=dicafile(os.path.join(modelDir,'dica.dat%s'%suffix)).read_data()
    comp=compfile(os.path.join(modelDir,'comp.ind%s'%suffix)).read_data()
    try:
        synSpec=spectrum(os.path.join(modelDir,'spct.dat%s'%suffix),usecols=(0,2))
    except:
        print "Trouble with Model %s"%id
        synSpec=createFakeSpec(dica)
    sbibData=sbibfile(os.path.join(modelDir,'sbib.dat%s'%suffix)).read_data()
    return {'aspect':synSpec,'dica':dica,'comp':comp,'sbib':sbibData,'id':id}
"""    
class modelGrid(object):
    def __init__(self,models,shape=None,origspect=None):
        if shape!=None: self.data=np.array(models).reshape(shape)
        else: self.data=models
        #if origspect==None:
        #    print('You didnt specify a real spectra, please do so. I won\'t crash, but I\'m no fun to work with')
    def getKeyword(self,keyword):
        func=np.vectorize(lambda item:item[keyword])
        return modelGrid(func(self.data))
    def getAspect(self):
        return self.getKeyword('aspect')
    def getDica(self):
        return self.getKeyword('dica')
    def getComp(self):
        return self.getKeyword('comp')
    def getSbib(self):
        return self.getKeyword('sbib')
    def getIntTrapz(self,xref=None):
        integrate=np.vectorize(lambda item:item.intTrapz())
        if xref==None:
            integral=integrate(self.getAspect().data)
        else:
            interpolate=np.vectorize(lambda item:item.interpolate(xref=xref))
            integral=integrate(interpolate(self.getAspect().data))
        return integral
    def getDivSpect(self,origspect):
        divide=vectorize(lambda item:origspect/item)
        return modelGrid(divide(self.getAspect().data))
    def getXiSquared(self,origspect):
        print
"""
class model(object):
    def __init__(self,basePath='.',param=None,origSpec=None):
        if param!=None:
            self.param=param
        else:
            dicaData=fileio.dicafile(os.path.join(basePath,'dica.dat')).read_data()
            compData=fileio.compfile(os.path.join(basePath,'comp.ind')).read_data()
            dica=paramMod.dica(dicaData)
            comp=paramMod.comp(compData)
            self.param=paramMod.param(initDica=dica,initComp=comp)
        
        aSpecPath=os.path.join(basePath,'spct.dat')
        try:
            self.aSpec=spectrum(aSpecPath,usecols=(0,2))
            self.sbib=fileio.sbibfile(os.path.join(basePath,'sbib.dat')).read_data()
            self.llist=self.sbib['llist']
            self.wParams=fileio.ststfile(os.path.join(basePath,'stst.dat')).getWParams()
            self.specFlag=0
        except:
            print "Creating fake Spectrum @%s"%basePath
            self.aSpec=spectrum(zip(np.linspace(2000,10000,10),range(1,11)))
            self.sbib={'llist':[]}
            self.llist=self.sbib['llist']
            self.wParams=[]
            self.specFlag=-1
        """
        if (not os.path.exists(aSpecPath)) or os.stat(aSpecPath).st_size<10:
            print "Creating fake Spectrum @%s"%basePath
            self.aSpec=spectrum(zip(np.linspace(2000,10000,10),range(1,11)))
            self.sbib={'llist':[]}
            self.llist=self.sbib['llist']
            self.wParams=[]
        else:
            self.aSpec=spectrum(aSpecPath,usecols=(0,2))
            self.sbib=fileio.sbibfile(os.path.join(basePath,'sbib.dat')).read_data()
            self.llist=self.sbib['llist']
            self.wParams=fileio.ststfile(os.path.join(basePath,'stst.dat')).getWParams()
        """
        self.log=list(file(os.path.join(basePath,'fica.log')))
        #self.error=list(file(os.path.join(basePath,'error.log')))
        tmpAspec=self.aSpec.interpolate(xref=origSpec.x)
        self.divSpec=tmpAspec/origSpec
        self.subSpec=fit.getSubSpec(tmpAspec,origSpec)
    def __getitem__(self,key):
        if key.lower()=='llist':
            return self.llist
        elif key.lower()=='w':
            if self.wParams!=[]:
                return self.wParams[-1][0]['w']
            else:
                return -1
        elif key.lower().startswith('aspec'):
            return self.aSpec
        elif key.lower().startswith('divpec'):
            return self.divSpec
        elif key.lower().startswith('addspec'):
            return self.addSpec
        elif key.lower().startswith('subspec'):
            return self.subSpec
        elif key.lower().startswith('err'):
            return self.error
        elif key.lower().startswith('divspec'):
            return self.divSpec
        else:
            return self.param[key]

       
class modelGrid(object):
    def __init__(self,runDirs=None,multiParam=None,paramList=None,origSpec=None):
        tmpParam=[]
        if paramList!=None:
            self.grid=np.array(paramList)
            del paramList
        else:
            for dir,param in zip(runDirs,multiParam.paramGrid.flatten()):
                #print "Reading Directory %s"%dir
                tmpParam.append(model(param,dir))
            self.grid=np.array(tmpParam)
            del tmpParam
        #self.grid.reshape(multiParam.paramGrid.shape)
        if origSpec!=None: self.origSpec=origSpec
        else: self.origSpec=config.getOrigSpec(preProcess=True)
        self._initSpecs()
        self.specMask=self.specFlag==0
        self._initSpecs()
    def _filterSpectra(self):
        self.grid=self.grid[self.specMask]
    def _initSpecs(self):
        self._getSubSpec()
        self._getDivSpec()
        #self._getAddSpec()
        self._getSpecFlag()
        self._getFitness()
    def __getitem__(self,key):
        if key=='tobedetermined':
            print
        elif isinstance(key,int):
            return self.grid[key]
        elif key.lower().startswith('fit'):
            return self.fitness
        elif key.lower().startswith('divspec'):
            return self.divSpec
        elif key.lower().startswith('specflag'):
            return self.specFlag
        elif key.lower().startswith('subspec'):
            return self.subSpec
        elif key.lower().startswith('addspec'):
            return self.addSpec
        elif key.lower().startswith('llis'):
            return [item['llist'] for item in self.grid]
        else:
            getFunc=np.vectorize(lambda item:item[key])
            return getFunc(self.grid)
    def _getSubSpec(self):
        #vecFunc=np.vectorize(lambda item:fit.getSubSpec(item['aspec'],self.origSpec))
        #self.subSpec=vecFunc(self.grid)
        self.subSpec=np.array([item.subSpec for item in self.grid])
    def _getAddSpec(self):
        raise Exception('add not implemented yet')
        vecFunc=np.vectorize(lambda item:fit.getAddSpec(item['aspec'],self.origSpec))
        self.addSpec=vecFunc(self.grid)
    def _getSpecFlag(self):
        vecFunc=np.vectorize(lambda item:item.specFlag)
        self.specFlag=vecFunc(self.grid)
    def _getDivSpec(self):
        self.divSpec=np.array([item.divSpec for item in self.grid])
    def _getFitness(self):
        self.fitness=np.array([item.fitness for item in self.grid])
        #def divSpecGetter(item):
        #    tmpAspec=item['aspec']
        #    tmpAspec=tmpAspec.interpolate(xref=self.origSpec.x)
        #    return tmpAspec/self.origSpec
        #vecFunc=np.vectorize(divSpecGetter)
        #self.divSpec=vecFunc(self.grid)
        
class lumVphGrid(modelGrid):
    def getLums(self):
        return self.getDica().getKeyword('log_lbol').data
    def getVphs(self):
        return self.getDica().getKeyword('v_ph').data
    def getLumsVphsLabels(self,fmt='%.2f'):
        labelStr='lum=%s vph=%s'%(fmt,fmt)
        labels=[]
        for lum,vph in zip(self.getLums().data.flatten(),self.getVphs().data.flatten()):
            labels.append(labelStr%(lum,vph))
        return np.array(labels).reshape(self.data.shape)
          
        
#simple Functions to extract merits
#getGridInt=np.vectorize(lambda item: fit.getDiffIntBin(item,1)[1][0])
getGridInt=np.vectorize(lambda item: fit.getInt(item))
getGridUV=np.vectorize(lambda item: fit.getUVInt(item,norm=False))
getGridOptical=np.vectorize(lambda item: fit.getInt(item,lower=3950.,norm=False))
getGridUVComp=np.vectorize(lambda item: fit.getUVIntComp(item))
getGridSlope=np.vectorize(lambda item: fit.getIntSlope(item))
getGridSlope=np.vectorize(lambda item: fit.getSlope(item))