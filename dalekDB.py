#module to interact with the sqlite database storing the GA data
from glob import glob
import os
import sqlite3

import config
paramDir=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'conf.d/')


def importOldDalekDir(path, dbPath):
    #importing dalek conf and spectra from before the sqlite era
    snSpectraDirs = glob('????-??-??T??-??-??')
    for specDir in snSpectraDirs:
        
        origSpec = os.path.join(os.path.abspath(specDir),
                                'spectra','origspect.dat')
        
def importOldSNDir(dalekDir, dbPath):
    #importing the old sn directory to sqlite from before the sqlite era
    pass

def importOldConf(dalekDir, conn):
    #importing the old sn directory to sqlite from before the sqlite era
    SNConfig = config.getMainConfig(dalekDir)
    conn = sqlite3.connect(dbPath)
    for item in SNCONFIG.items('snconf'):
        conn.execute('insert into SN_PARAM (NAME, VALUE, VALUE_TYPE) '
                     'values (?, ?, ?)', item + ('float',))
        
def createTestDB():
    #creating tmp database to play with
    schema = file(os.path.join(paramDir, 'dalekDB.schema')).read()
    #deleting old play database
    conn = sqlite3.connect(":MEMORY:")
    conn.execute(schema)
    return conn
