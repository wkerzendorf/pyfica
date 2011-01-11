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

def importOldDalekDir(path, conn):
    #importing dalek conf and spectra from before the sqlite era
    pass
        
def importOldSNDir(path, conn):
    #importing the old sn directory to sqlite from before the sqlite era
    snSpectraDirs = glob(os.path.join(os.path.abspath(path), '????-??-??T??-??-??'))
    for specDir in snSpectraDirs:
        snSpec = spectrum(os.path.join(os.path.abspath(specDir),
                                'spectra','origspect.dat'))
        zSnSpec = sqlite3.Binary(zlib.compress(cPickle.dumps(snSpec)))
        snDate = datetime.datetime.strptime(os.path.basename(specDir),
                                            '%Y-%m-%dT%H-%M-%S')
        conn.execute('insert into SN_SPECTRA (DATE, SPECTRUM) '
                     'values (?, ?)', (snDate, zSnSpec))
        

def importOldConf(dalekDir, conn):
    #importing the old sn directory to sqlite from before the sqlite era
    SNConfig = config.getMainConfig(dalekDir)
    for item in SNConfig.items('snconf'):
        conn.execute('insert into SN_PARAM (NAME, VALUE, VALUE_TYPE) '
                     'values (?, ?, ?)', item + ('float',))
        
def createTestDB():
    #creating tmp database to play with
    schema = file(os.path.join(paramDir, 'dalekDB.schema')).read()
    #deleting old play database
    conn = sqlite3.connect(":memory:")
    conn.executescript(schema)
    return conn
