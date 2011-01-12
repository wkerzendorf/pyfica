#module to interact with the sqlite database storing the GA data
from glob import glob
import os
import sqlite3
import cPickle
import zlib
import datetime
import config
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
                     'values (?, ?)', (snDate, snSpec))
        

def importOldConf(dalekDir, conn):
    #importing the old sn directory to sqlite from before the sqlite era
    SNConfig = config.getMainConfig(dalekDir)
    for item in SNConfig.items('snconf'):
        conn.execute('insert into SN_PARAM (NAME, VALUE, VALUE_TYPE) '
                     'values (?, ?, ?)', item + ('float',))

def convertZipPickle(blob):
    return cPickle.loads(zlib.decompress(blob))

def createTestDB():
    #creating tmp database to play with
    schema = file(os.path.join(paramDir, 'dalekDB.schema')).read()
    #deleting old play database
    conn = sqlite3.connect(":memory:", detect_types=sqlite3.PARSE_DECLTYPES)
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
    return wlFinalGrid[3:-3]



def insertFicaModel(conn, model, dicaID=None):
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
    
    #merging the dataset
    curs.execute('insert into FICA_MODEL'
                 '(MACHINE, TIME, W, ERROR,'
                 'ABUNDANCE_ID, DICA_ID, LUMVPH_ID, SPECTRUM_ID)'
                 'VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                 ("None", 0, 0, str(model.error), compID, dicaID, lumVphID, specID))
    

sqlite3.register_converter("PYSPEC_ONED", convertZipPickle)
sqlite3.register_converter("NP_ARRAY", convertZipPickle)

