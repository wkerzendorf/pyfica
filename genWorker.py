import site
import inspect
import os
import shutil
import subprocess
import cPickle as pickle
import time
#import json
if __name__ == '__channelexec__':
    startTime = time.time()
    
    #Setting up the proper paths on the end node
    for sitedir in channel.receive():
        site.addsitedir(sitedir)
    from pyfica import read,param,fileio
    #objectString=channel.receive()
    dica,comp,ficaWorkDir,ficaBin,i=channel.receive()
    try:
        os.makedirs(ficaWorkDir)
    except:
        shutil.rmtree(ficaWorkDir)
        os.makedirs(ficaWorkDir)
    #param.write2file(ficaWorkDir)
    fileio.dicafile(os.path.join(ficaWorkDir,'dica.dat'),'w').write_data(dica)
    fileio.compfile(os.path.join(ficaWorkDir,'comp.ind'),'w').write_data(comp)
    
    #print "Running %s on %s in dir %s"%(fica_bin,currentNode,ficaWorkDir)
    proc=subprocess.Popen([ficaBin],stdout=-1,cwd=ficaWorkDir,shell=True,close_fds=True)
    #Saving Log File:
    file(os.path.join(ficaWorkDir,'fica.log'),'w').write(proc.stdout.read())
    proc.stdout.close()
    
    channel.send((i,ficaWorkDir))


