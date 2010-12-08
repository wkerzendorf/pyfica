#Used to read Models
from pyspec.spectrum import spectrum
from .fileio import dicafile,compfile,sbibfile
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