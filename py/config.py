import ConfigParser,os,datetime,ephem,string,re
from glob import glob
import fileio,initialize
paramDir=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'conf.d/')

def getMainConfigDir(confDirName='conf'):
    curpath=os.getcwd()
    while True:
        if os.path.exists(os.path.join(curpath,confDirName)):
            return os.path.join(curpath,confDirName)
        if curpath=='/':
            raise Exception('Configuration directory %s not found in path'%confDirName)
        curpath=os.path.split(curpath)[0]
def getMainConfig(confDirName='conf',confFileName='sn.cfg'):
    config=ConfigParser.ConfigParser()
    print "Reading Configuration file %s"%os.path.join(getMainConfigDir(),confFileName)
    config.read(os.path.join(getMainConfigDir(),confFileName))
    return config
def getDir2DateTime(pathName):
    dirName=os.path.basename(getCurMainDir())
    return datetime.datetime(*map(int,dirName.replace('T','-').split('-')))
def getTimeFromMax(config=None):
    if config==None:
        config=getMainConfig()
    return ephem.julian_date(getDir2DateTime(os.getcwd()))-config.getfloat('snconf','tbmax')
def getTimeFromExplosion(config=None,tRise=19.5):
    return getTimeFromMax(config=config)+tRise
def getVanillaDica():
    dicaConf=os.path.join(getMainConfigDir(),'dica.dat')
    return fileio.dicafile(dicaConf).read_data()
    #vanillaDica=
def getInitDicaParam(tRise=19.5):
    snConfig=getMainConfig()
    return {'z':snConfig.getfloat('snconf','z'),
                 'e_b-v_gal': snConfig.getfloat('snconf','extgal'),
                 'e_b-v_host': snConfig.getfloat('snconf','exthost'),
                 'm-m': snConfig.getfloat('snconf','mu'),
                 't':getTimeFromExplosion(tRise=tRise)}
def getInitDica(tRise=19.5):
    dica=getVanillaDica()
    snConfig=getMainConfig()
    dica.update(getInitDicaParam(tRise=tRise))
    return dica
def getCurMainDir():
    curpath=os.getcwd()
    while True:
        dirName=os.path.basename(curpath)
        matchingDate=re.match('\d{4}\-\d{2}\-\d{2}T\d{2}\-\d{2}\-\d{2}',dirName)
        if matchingDate!=None:
            break
        elif curpath=='/':
            raise Exception('Configuration directory %s not found in path'%confDirName)
        curpath=os.path.split(curpath)[0]
    return curpath
def getLastMainDir():
    epochPaths=glob(os.path.dirname(getMainConfigDir())+'/????-??-??T??-??-??')
    epochPaths.sort()
    epochs=[os.path.basename(item) for item in epochPaths]
    curMainDir=getCurMainDir()
    index=epochs.index(os.path.basename(curMainDir))
    if index==0: raise Exception('You are in the earliest observation directory. There is no observation before that!')
    else:
        return epochPaths[index-1]
def getInitComp(model='w7'):
    if model=='w7':
        w7Data=initialize.readW7Data()
        initialize.getW7Comp(w7Data)
def getVanillaComp(compType='comp_0'):
    compConf=os.path.join(getMainConfigDir(),'initialSetup.cfg')
    config=ConfigParser.ConfigParser()
    config.read(compConf)
    tmpOutput=[]
    return dict([(string.upper(elem[0])+elem[1:],config.getfloat(compType,elem)) for elem,val in config.items(compType)])
    
def getMachineConfig():
    config=ConfigParser.ConfigParser()
    config.read(os.path.join(paramDir,'machines.cfg'))
    machineConfig={}
    for section in config.sections():
        machineConfig[section]=dict(config.items(section))
        machineConfig[section]['processors']=int(machineConfig[section]['processors'])
    return machineConfig