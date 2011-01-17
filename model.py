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

import dalekDB


class model(object):
    #class to store one model
    @classmethod
    def fromDB(cls, conn, modelID, origSpecID):
        curs = conn.cursor()
    
        #retrieving origSpec:
        if origSpecID != None:
            origSpec = curs.execute('select spectrum from sn_spectra where id=%s' % origSpecID)
        else: origSpec = None
        #Retrieving the Model
        (machineName, execTime, wFactor, errorString, ficaLog,
        abundanceID, dicaID,
        lumVphID, spectrumID) = curs.execute('select MACHINE, TIME, W, ERROR, FICA_LOG, '
                     'ABUNDANCE_ID, DICA_ID, LUMVPH_ID, SPECTRUM_ID '
                     'from FICA_MODEL where FICA_MODEL.ID=%s' % modelID).fetchall()[0]
        #getting dica params
        colNames = zip(*curs.execute('PRAGMA table_info(fica_dica)').fetchall())[1]
        colNames = map(str, colNames)
        colValues = curs.execute('select * from fica_dica where id=%s' % dicaID).fetchall()[0]
        dicaDict = dict(zip([dalekDB.convertFields2Dica[item] for item in colNames[1:]], colValues[1:]))
        lum, vph = curs.execute('select LUM, VPH from FICA_LUMVPH where FICA_LUMVPH.ID=%s' % lumVphID).fetchall()[0]
        dicaDict['log_lbol'] = lum
        dicaDict['v_ph'] = vph
        
        dica = paramMod.dica(initDica=dicaDict, mode='fromDict')
        
        
        #getting abundances
        colNames = zip(*curs.execute('PRAGMA table_info(fica_abundance)').fetchall())[1]
        colNames = map(str, colNames)
        colValues = curs.execute('select * from fica_abundance where id=%s' % abundanceID).fetchall()[0]
        compDict = dict(zip(colNames[1:], colValues[1:]))
        comp = paramMod.comp(initComp=compDict, t=dica['t'])
        comp._setNiDecay()
        
        #getting aSpec
        wl = dalekDB.createWLGrid(dicaDict['wl']*1e4, dicaDict['grid']*1e4, dicaDict['mu'])
        intens = curs.execute('select spectrum from fica_spectrum where id=%s' % abundanceID).fetchall()[0][0]
        aSpec = spectrum(wl,intens)
        #getting llist
        colValues = curs.execute('select eqw, shift, rest, atom, ion, param1, param2, param3 '
                                 'from FICA_LLIST where model_id=%d' % modelID).fetchall()
        #checking if llist exists for current model
        if colValues == []: llist = None
        else:
            colNames = zip(*curs.execute('PRAGMA table_info(fica_abundance)').fetchall())[1]
            colNames = [(item.lower(), '|S2') if item=='ATOM' else (item.lower(), float)
                for item in colNames[2:]]
            llist = np.array(colValues, dtype=colNames)
        
        """ Commented out until wParams becomes important, W is safed none the less    
        #getting wParams
        colValues = curs.execute('select XS, VS, LOGRH, TE, TR, W '
                                 'from FICA_WPARAM where FICA_WPARAM.model_id=%d' % model_id).fetchall()
        #checking if WParams exists for current model
        if colValues == []: llist = None
        else:
            colNames = zip(*curs.execute('PRAGMA table_info(fica_WPARAM)').fetchall())[1]
            colNames = [(item.lower(), '|S2') if item=='ATOM' else (item.lower(), float)
                for item in colNames[2:]]
            llist = np.array(colValues, dtype=colNames)
        """
        wParam=None
            
        curParam = paramMod.param(initDica=dica, initComp=comp)
        
        return cls(aSpec, curParam, wFactor, machineName=machineName, execTime=execTime,
                 wParam=wParam, error=errorString, ficaLog = ficaLog,
                 llist=None, origSpec=origSpec)
        
        
    @classmethod
    def fromPath(cls, basePath='.',machineName=None, param=None,origSpec=None):
        if param==None:
            dicaData=fileio.dicafile(os.path.join(basePath,'dica.dat')).read_data()
            compData=fileio.compfile(os.path.join(basePath,'comp.ind')).read_data()
            dica=paramMod.dica(dicaData)
            comp=paramMod.comp(compData)
            param = paramMod.param(initDica=dica,initComp=comp)
        
        aSpecPath=os.path.join(basePath,'spct.dat')
        try:
            aSpec=spectrum(aSpecPath,usecols=(0,2))
            sbib=fileio.sbibfile(os.path.join(basePath,'sbib.dat')).read_data()
            llist=self.sbib['llist']
            wParams=fileio.ststfile(os.path.join(basePath,'stst.dat')).getWParams()
            specFlag=0
        except:
            print "Creating fake Spectrum @%s"%basePath
            aSpec=spectrum(zip(np.linspace(2000,10000,10),range(1,11)))
            sbib={'llist':[]}
            llist=self.sbib['llist']
            wParams=[]
            specFlag=-1

        log=list(file(os.path.join(basePath,'fica.log')))
        error=list(file(os.path.join(basePath,'error.log')))
        
        
        if self.wParams!=[]:
            self.w = self.wParams[-1][0]
        else:
            self.w = -1
        
        return cls(param, w, machine, execTime, wParam, error, log, llist, origSpec)    
    def __init__(self, aSpec, param, w, machineName=None, execTime=None,
                 wParam=None, error=None, ficaLog = None,
                 llist=None, origSpec=None, specFlag = -1):
        self.param = param
        self.w = w
        self.machine = machineName
        self.execTime = execTime
        self.wParam = wParam
        self.log = ficaLog
        self.error = error
        self.llist = llist
        self.origSpec = None
        self.specFlag = specFlag
        self.aSpec = aSpec
        #Initializing subspec, divspec....
        if origSpec!=None:
            tmpAspec=self.aSpec.interpolate(xref=origSpec.x)
            self.divSpec=tmpAspec/origSpec
            self.subSpec=fit.getSubSpec(tmpAspec,origSpec)
        else:
            self.subSpec = None
            self.divSpec = None
    
    def __getitem__(self,key):
        if key.lower()=='llist':
            return self.llist
        elif key.lower()=='w':
            return self.w
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
    def toDB(self, conn, dicaID=None, storeLList=False, storeWParam=False):
        curs = conn.cursor()
        dica = self.param.dica.data.copy()
        comp = self.param.comp.data.copy()
        aSpec = self.aSpec.y
        comp.pop('Ni')
        comp.pop('Co')
        comp.pop('Fe')
        
        lum = dica.pop('log_lbol')
        vph = dica.pop('v_ph')
        
        dicaFields = [convertDica2Fields[item] for item in dica.keys()]
        dicaValues = dica.values()
        
        compFields = comp.keys()
        compValues = comp.values()
        
        #Inserting dica values
        if dicaID == None:
            curs.execute('insert into FICA_DICA (%s) values  (%s)'
                         % (','.join(dicaFields), ','.join('?' * len(dicaValues))),
                        dicaValues)
            dicaID = curs.lastrowid
        
        #Inserting lum, vph values
        curs.execute('insert into FICA_LUMVPH (LUM, VPH) VALUES (?, ?)', (lum,vph) )
        lumVphID = curs.lastrowid
        
        #Inserting comp values
        curs.execute('insert into FICA_ABUNDANCE (%s) VALUES (%s)'
                     % (','.join(compFields), ','.join('?' * len(compValues))),
                        compValues)
        compID = curs.lastrowid
        
            
        #Inserting spectrum wl values
        zASpec = sqlite3.Binary(zlib.compress(cPickle.dumps(aSpec)))
        curs.execute('insert into FICA_SPECTRUM (SPECTRUM) VALUES (?)', (zASpec,))
        specID = curs.lastrowid
        if not hasattr(model,'error'):
            self.error = "None"
            
        if self.log != None:
            ficaLog = sqlite3.Binary(zlib.compress(cPickle.dumps(model.log)))
        else:
            ficaLog = 'None'
        #merging the dataset
        curs.execute('insert into FICA_MODEL'
                     '(MACHINE, TIME, W, ERROR, FICA_LOG, '
                     'ABUNDANCE_ID, DICA_ID, LUMVPH_ID, SPECTRUM_ID)'
                     'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
                     ("None", 0, 0, str(model.error), ficaLog, compID, dicaID, lumVphID, specID))
        #returning FICA_MODEL_ID
        modelID = curs.lastrowid
        
        #Storing Line List
        if storeLList:
            #Inserting line list
            for line in self.llist:
                curLine = list(line)
                curLine[4]=transformIon[curLine[4].lower()]
                curLine = tuple([modelID] + [float(item) if i!=3 else str(item) for i, item in enumerate(curLine)])
                curs.execute('insert into FICA_LLIST (MODEL_ID, EQW, SHIFT, REST, ATOM, ION, PARAM1, PARAM2, PARAM3)'
                             'values  (?, ?, ?, ?, ?, ?, ?, ?, ?)',
                             curLine)
        #Storing whole wParams
        if storeWParam:
            for i, wParam in enumerate(model.wParams):
                for line in wParam:
                    curLine=[modelID, i] + map(float,list(line))
                    curs.execute('insert into FICA_WPARAMS (MODEL_ID, WSET_ID, XS, VS, LOGRH, TE, TR, W)'
                             'values  (?, ?, ?, ?, ?, ?, ?, ?)',
                             curLine)
                     
        #Returning the modelID    
        return modelID
        
       
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
        
       
#simple Functions to extract merits
#getGridInt=np.vectorize(lambda item: fit.getDiffIntBin(item,1)[1][0])
getGridInt=np.vectorize(lambda item: fit.getInt(item))
getGridUV=np.vectorize(lambda item: fit.getUVInt(item,norm=False))
getGridOptical=np.vectorize(lambda item: fit.getInt(item,lower=3950.,norm=False))
getGridUVComp=np.vectorize(lambda item: fit.getUVIntComp(item))
getGridSlope=np.vectorize(lambda item: fit.getIntSlope(item))
getGridSlope=np.vectorize(lambda item: fit.getSlope(item))