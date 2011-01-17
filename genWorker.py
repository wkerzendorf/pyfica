import site
import inspect
import os
import shutil
import subprocess
import cPickle as pickle
import time


#import json
if __name__ == '__channelexec__':
    
    
    #Setting up the proper paths on the end node
    for sitedir in channel.receive():
        site.addsitedir(sitedir)
    #after adding paths import pyfica modules
    from pyfica import model,param,fileio, elauncher
    
    protocol = channel.receive()
    if protocol == 'centralRead':
        dica,comp,ficaWorkDir,ficaBin,i=channel.receive()
    if protocol == 'localRead':
        startTime = time.time()
        origSpecPickle,dica,comp,ficaWorkDir,ficaBin,i=channel.receive()
        origSpec = pickle.loads(origSpecPickle)
    try:
        os.makedirs(ficaWorkDir)
    except:
        shutil.rmtree(ficaWorkDir)
        os.makedirs(ficaWorkDir)
    
    #Writing configuration files
    fileio.dicafile(os.path.join(ficaWorkDir,'dica.dat'),'w').write_data(dica)
    fileio.compfile(os.path.join(ficaWorkDir,'comp.ind'),'w').write_data(comp)
    
    #Launching Fica
    proc=subprocess.Popen([ficaBin],stdout=-1,cwd=ficaWorkDir,shell=True,close_fds=True)
    
    #Saving Log File:
    file(os.path.join(ficaWorkDir,'fica.log'),'w').write(proc.stdout.read())
    proc.stdout.close()
    
    
    
    #reading or sending ficaWorkDir depending on protocol being used
    if protocol == 'centralRead':
        channel.send((protocol, i,ficaWorkDir))    
    elif protocol == 'localRead':
        model = elauncher.getModel(ficaWorkDir, origSpec)
        model.time = time.time() - startTime
        
    


