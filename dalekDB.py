#module to interact with the sqlite database storing the GA data
from glob import glob
import os
import sqlite3
import cPickle
import zlib
import datetime
import config
import param
paramDir=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'conf.d/')
from pyspec.spectrum import spectrum

import numpy as np

convertDica2Fields = {'e_b-v_host':'EXT_HOST',
                          'e_b-v_gal':'EXT_GAL',
                        'wl':'WL',
                        'js':'JS',
                        'm-m':'DIST_MODULUS',
                        'grid':'GRID',
                        'xe1':'XE1',
                        'log_l_low_high':'LOG_L_LOW_HIGH',
                        'np5':'NP5',
                        'itt':'ITT',
                        'lg_tau':'LOG_TAU',
                        'em_high':'EM_HIGH',
                        'em_low':'EM_LOW',
                        'kb':'KB',
                        'z':'Z',
                        'mb':'MB',
                        'intlow':'INTLOW',
                        'inthigh':'INTHIGH',
                        'nc':'NC',
                        'chl':'CHL',
                        'mu':'MU',
                        'kr':'KR',
                        't':'T',
                        'tb':'TB',
                        'options':'OPTIONS'}

convertFields2Dica = dict([item[::-1] for item in convertDica2Fields.items()])

def importOldDalekDir(path, conn):
    #importing dalek conf and spectra from before the sqlite era
    pass
        
def importOldSNDir(path, conn):
    #importing the old sn directory to sqlite from before the sqlite era
    snSpectraDirs = glob(os.path.join(os.path.abspath(path), '????-??-??T??-??-??'))
    for specDir in snSpectraDirs:
        snSpec = spectrum(os.path.join(os.path.abspath(specDir),
                                'spectra','origspect.dat'))
        snDate = datetime.datetime.strptime(os.path.basename(specDir),
                                            '%Y-%m-%dT%H-%M-%S')
        conn.execute('insert into SN_SPECTRA (DATE, SPECTRUM) '
                     'values (?, ?)', (snDate, makeZipPickle(snSpec)))
        

def importOldConf(dalekDir, conn):
    #importing the old sn directory to sqlite from before the sqlite era
    SNConfig = config.getMainConfig(dalekDir)
    for item in SNConfig.items('snconf'):
        conn.execute('insert into SN_PARAM (NAME, VALUE, VALUE_TYPE) '
                     'values (?, ?, ?)', item + ('float',))

def convertZipPickle(blob):
    return cPickle.loads(zlib.decompress(blob))

def makeZipPickle(object):
    return sqlite3.Binary(zlib.compress(cPickle.dumps(object)))

def createTestDB(dbName=':memory:'):
    #creating tmp database to play with
    schema = file(os.path.join(paramDir, 'dalekDB.schema')).read()
    #deleting old play database
    conn = sqlite3.connect(dbName, detect_types=sqlite3.PARSE_DECLTYPES)
    conn.executescript(schema)
    return conn

def createWLGrid(wlStart, wlEnd, wlSteps):
    #Creating a wavelength grid from the data given in fica
    logWLStep = (np.log10(wlEnd) - np.log10(wlStart))/wlSteps
    wlVacGrid = 10**(np.log10(wlStart) + np.arange(wlSteps)*logWLStep)
    if wlStart < 2000:
        raise NotImplementedError()
        for i in np.arange(wlSteps):
            if wlVacGrid[i] > 2000.0:
                rw2 = 1/(wlVacGrid[i] * 1.0e-4)**2
                wlAirGrid.append(wlVacGrid[i] / (1.0 + 64.328e-6 + 29498.1e-6 /
                                                 (146.0 - rw2) + 255.4e-6/(41.0 - rw2)))
            else:
                wlAirGrid.append(wlVacGrid[i])
    else:
        wlFinalGrid = wlVacGrid
    #return wlAirGrid
    wlFinalGrid = 0.5*(wlFinalGrid[:-1] + wlFinalGrid[1:])
    return wlFinalGrid[3:-4]

def insertDica(conn, dica):
    curs = conn.cursor()
    dicaDict = dica.data.copy()
    
    #removing lum and vph
    dicaDict.pop('log_lbol')
    dicaDict.pop('v_ph')
    
    dicaFields = [convertDica2Fields[item] for item in dicaDict.keys()]
    dicaValues = dicaDict.values()
    
    curs.execute('insert into FICA_DICA (%s) values  (%s)'
                % (','.join(dicaFields), ','.join('?' * len(dicaValues))),
                    dicaValues)
    dicaID = curs.lastrowid
    
    return dicaID


def insertGAIndividual(conn, GARunID, modelIDs, fitness):
    curs = conn.cursor()
    for fit, mID in zip(fitness, modelIDs):
        curs.execute('insert into GA_INDIVIDUAL(GENERATION_ID, MODEL_ID, FITNESS) '
                     'values(?, ?, ?)' % (GARunID, mID, fit))

def insertFicaModel(conn, model, dicaID=None, storeLList=False, storeWParam=False):
    curs = conn.cursor()
    dica = model.param.dica.data.copy()
    comp = model.param.comp.data.copy()
    aSpec = model.aSpec.y
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
        model.error = "None"
        
    if model.log != None:
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
        for line in model.llist:
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

def getFicaModel(conn, modelID, origSpecID=None):
    curs = conn.cursor()
    
    #retrieving origSpec:
    if origSpecID != None:
        origSpec = curs.execute('select spectrum from sn_spectra where id=%s' % origSpecID)
    
    #Retrieving the Model
    (machineName, execTime, wFactor, errorString,
    abundanceID, dicaID,
    lumVphID, spectrumID) = curs.execute('select MACHINE, TIME, W, ERROR, '
                 'ABUNDANCE_ID, DICA_ID, LUMVPH_ID, SPECTRUM_ID '
                 'from FICA_MODEL where FICA_MODEL.ID=%s' % modelID).fetchall()[0]
    #getting dica params
    colNames = zip(*curs.execute('PRAGMA table_info(fica_dica)').fetchall())[1]
    colNames = map(str, colNames)
    colValues = curs.execute('select * from fica_dica where id=%s' % dicaID).fetchall()[0]
    dicaDict = dict(zip([convertFields2Dica[item] for item in colNames[1:]], colValues[1:]))
    lum, vph = curs.execute('select LUM, VPH from FICA_LUMVPH where FICA_LUMVPH.ID=%s' % lumVphID).fetchall()[0]
    dicaDict['log_lbol'] = lum
    dicaDict['v_ph'] = vph
    
    dica = param.dica(initDica=dicaDict, mode='fromDict')
    return dica
    
    #getting abundances
    colNames = zip(*curs.execute('PRAGMA table_info(fica_abundance)').fetchall())[1]
    colNames = map(str, colNames)
    colValues = curs.execute('select * from fica_abundance where id=%s' % abundanceID).fetchall()[0]
    compDict = dict(zip(colNames[1:], colValues[1:]))
    comp = param.comp(initComp=compDict, t=dica['t'])
    comp._setNiDecay()
    
    curParam = param.param(initDica=dica, initComp=comp)
    
    
    
    
    
    #Retrieving 
    

sqlite3.register_converter("PYSPEC_ONED", convertZipPickle)
sqlite3.register_converter("NP_ARRAY", convertZipPickle)
sqlite3.register_converter("ZLOG", convertZipPickle)
transformIon = dict(i=1, ii=2, iii=3, iv=4, v=5)
