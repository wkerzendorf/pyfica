import subprocess
import os
import numpy as np
from ConfigParser import ConfigParser



import config

from os import chdir,mkdir, getcwd, path, symlink, system
from shutil import copy, rmtree, move
from copy import copy as objcopy
from .fileio import compfile, dicafile
from .util import setAbundances
from commands import getoutput
from glob import glob
#from .util import find_last
from pyphot import datafile
from numpy import array
import threading
import time, sys
import pprint
from datetime import datetime, timedelta
import logging
import sys
from ConfigParser import ConfigParser

localFicaBin=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin/bin/fica.exe')
#confDir=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'conf.d')
#reading
machineConf=config.getMachineConf()
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
    datefmt='%m-%d %H:%M',
    filename='runlog.log',
    filemode='w')
mdLog=logging.getLogger('pyfica.launcher')
#console = logging.StreamHandler()
#console.setLevel(logging.INFO)
#logging.getLogger('').addHandler(console)

from .fitelem import doDecay
#fica_bin='/Users/wkerzend/wdata/sn_rad_trans/fica/fica.exe'
#local_fica_bin=fica_bin
#no_processor=0
error_pref="----- ERROR ----- "
if path.exists('/Users/wkerzend'):
    conf_dir='/Users/wkerzend/wdata/sn_rad_trans/fica/conf/'
else:
    conf_dir='/home/wkerzend/wdata/sn_rad_trans/fica/conf/'
init_cond='init_cond/'
default_rundir='run_tmp/'
default_script_path='/priv/manana1/wkerzend/sn_rad_trans/fica/scripts/'
#Giving all machines x more processing slots
default_extraproc=0
del_proc=True
nice='nice'
default_machines=['myriad', 'maggot', 'miami', 'merino', 'minotaur', 'munch', 'macerate', 'magpie', 'manana', 'marvin', 'matrix']
checkUpTimeStep=timedelta(0, 50)

    
backup_files=['dica.dat', 'comp.ind', 'spct.dat', 'sbib.dat', 'fica.log', 'error.log','stst.dat','yhea.dat']
backup_suffix='.bak'
def multiDistLaunch(runDirs,baseDir='.',useMachines='default'):
    logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
    datefmt='%m-%d %H:%M',
    filename='%s/runlog.log',
    filemode='w')
    #Configuring which machines to use
    if useMachines=='default': useMachines=default_machines
    #Creating stack for each machine to check whats running
    runnings=dict([(machine, []) for machine in conf_names])
    #Creating process stack
    procStack=[]
    for i,runDir in enumerate(runDirs):
        tmpParam={'popen':None,
                    'runDir':runDir,
                    'mconf':None,
                    'runID':i,
                    'freq':0}
        procStack.append(tmpParam)
    
    while True:
        
        if procStack==[] and sum([len(item) for item in runnings.values()])==0: 
            #all processes launched and all processes returned, end launcher
            break
        elif procStack!=[]:
            #Checking run stacks for empty slots and adding processes            
            for machine in useMachines:
                while len(runnings[machine])<\
                    machine_conf[machine]['processors']+default_extraproc:
                    
                    if procStack!=[]: procParam=procStack.pop(0)
                    else: break

                    mdLog.debug("--------------------------------")
                    mdLog.info("Adding process %s to the machine %s"%(proc_param['run_id'], machine))
                    #Increasing frequency counter
                    freq=proc_param['freq']+1
                    #launching the process
                    proc_handle=launch(proc_param['dica'],
                                       proc_param['comp'],
                                       machine_conf[machine],
                                       run_id=proc_param['run_id'], 
                                       preped=True, 
                                       init_dica=init_dica, 
                                       init_comp=init_comp)
                    #Updating the stats
                    procParam.update({'popen':proc_handle,'freq':freq, 'mconf':machine_conf[machine]})
                    procParam.update({'startTime':datetime.now(),'checkTime':datetime.now()+checkUpTimeStep})
                    mdLog.debug("Added process %s"%proc_param['run_id'])
                    mdLog.debug("--------------------------------")
                    #append to queue
                    runnings[machine].append(procParam)
            
        #checking if processes are finished and evaluating what to do
        for machine in runnings.keys():
            exit_flags={}
            for j, irun in enumerate(runnings[machine]):
                merror_pref='%s: %s'%(machine, error_pref)
                if irun['popen'].poll()!=None:
                    #print merror_pref+"Stderr msg for proc %s is not readable. poll indicates %s"%(irun['run_id'], irun['popen'].poll())
                    #print "Trying to get error code from file"
                    try:
                        error_file=file("proc%s/error.log"%irun['run_id']).readlines()
                        if not error_file[-1].startswith('error'):
                            mdLog.debug(machine + " Found error.log with the following in it %s"%error_file)
                        return_code=int(error_file[-1].split()[1])
                    except:
                        return_code=-999
                        mdLog.error("error.log could not be found for process %s @ %s"%(irun['run_id'], machine))
                    
                    try:
                        return_code=int(return_code)
                    except:
                        mdLog.error("Process %s returned a non standard code %s @ %s"%(irun['run_id'], return_code, machine))
                        return_code=-999
                    
                    mdLog.info("Process %s has returned %s"%(irun['run_id'], return_code))
                    if return_code==0:
                        mdLog.debug('Marking %s for cleanup'%irun['run_id'])
                        exit_flags.update({j:'cleanup'})
                        #cleanups.append()
                    else:
                        if irun['freq']>3: 
                            mdLog.error("Process %s has failed more than 3 times @ %s"%(irun['run_id'], machine))
                            exit_flags.update({j:'multifail'})
                        if return_code==2:
                            mdLog.error('%s has encountered a stale nfs handle @ %s. Resubmitting to queue'%(irun['run_id'], machine))
                            exit_flag.update({j:'resubmit'})
                        elif return_code==-999:
                            mdLog.debug("User set return code on %s"%irun['run_id'])
                            raise Exception ('shit is happening: %s'%irun)
                            exit_flags.update({j:'resubmit'})
                        else: 
                            mdLog.error("Unknown return code %s @ %s"%(return_code, machine))
                else:
                    x=0
            sexit_index=exit_flags.keys()
            sexit_index.sort()
            sexit_index.reverse()
            for iexit_flag in sexit_index:
                exit_proc=runnings[machine].pop(iexit_flag)
                if exit_flags[iexit_flag]=='cleanup':
                    proc_cleanup(exit_proc)                
                elif exit_flags[iexit_flag]=='resubmit':
                    proc_cleanup(exit_proc)
                    exit_proc.update({'popen':None})
                    proc_params.insert(0, exit_proc)
                else:
                    mdLog.critical("%s had a critical error. It won't be resubmitted or cleaned up"%exit_proc['run_id'])
        for machine in runnings.keys():
            print "%s: %02d"%(machine, len(runnings[machine]))
        lines=len(runnings.keys())+1
        if proc_params!=[]:
            mdLog.info("Next process to launch is %s"%proc_params[0]['run_id'])
            print "Next process to launch is %s/%s"%(proc_params[0]['run_id'], len(dicas))
            sys.stdout.write("\x1b[%sA"%lines)
        else:
            print "Launched all parameter sets, Still %s items in queue"%\
            sum([len(item) for item in runnings.values()])
            sys.stdout.write("\x1b[%sA"%lines)
            mdLog.info("Launched all parameter sets, Still %s items in queue"%\
            sum([len(item) for item in runnings.values()]))
            for machine in runnings.keys():
                mdLog.info("%s: %s"%(machine, len(runnings[machine])))
        time.sleep(0.2)
    sys.stdout.write('\x1b[%sB'%lines)
    print "The directory %s took %s minutes "%(getcwd(), (datetime.now()-start).seconds/60.0)
    mdLog.info("The directory %s took %s minutes "%(getcwd(), (datetime.now()-start).seconds/60.0))
    chdir('../')



def multi_dist_launch(dicas, comps,  conf_names=['local'],run_dir=default_rundir,  init_dica=None,  init_comp=None):
    logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
    datefmt='%m-%d %H:%M',
    filename='runlog.log',
    filemode='w')

    if conf_names=='default': conf_names=default_machines
    
    start=datetime.now()
    #creating individual process stacks for each machine
    runnings=dict([(machine, []) for machine in conf_names])
    run_dir=path.join(run_dir, '')
    if not path.exists(run_dir): mkdir(run_dir)
    else:
        print "Deleting existing runDir %s"%run_dir
        rmtree(run_dir)
        mkdir(run_dir)
    chdir(run_dir)
    proc_params=[]
    #Creating the params stack:
    for i, (dica, comp) in enumerate(zip(dicas, comps)):
        tmp_param={'popen':None,
                    'dica':dica, 
                    'comp':comp,
                    'mconf':None, 'run_id':i, 'freq':0}
        proc_params.append(tmp_param)
    #preparing the directories
    #if init_dica==None: dica=dicafile('../'+init_cond+'dica.dat').read_data()
    #else: dica=init_dica
    #if init_comp==None: comp=compfile('../'+init_cond+'comp.ind').read_data()
    #else: comp=init_comp
    print "Preparing directories"
    for i, proc_param in enumerate(proc_params):
        run_id=proc_param['run_id']
        dica_params=proc_param['dica']
        comp_params=proc_param['comp']
        #Changing the parameters parsed to the launcher
        dica_new=dica
        dica_new.update(dica_params)
        comp_new=objcopy(comp)
        #comp_new=setAbundances(comp_new, comp_params)
        proc_dir=prep_proc('proc%s'%run_id, wait=False)
        dicafile(proc_dir+'dica.dat', 'w').write_data(dica_new)
        compfile(proc_dir+'comp.ind', 'w').write_data(comp_new)
        proc_params[i]['comp']=objcopy(comp_new)
        #doDecay(proc_dir)
    time.sleep(1)
    #raise Exception('test')
    print "Starting FICA..."
    while True:
        #
        if proc_params==[] and sum([len(item) for item in runnings.values()])==0: 
            #pp.pprint(runnings)
            break
        elif proc_params!=[]:
            #Checking run stacks for empty slots and adding processes            
            for machine in conf_names:
                while len(runnings[machine])<\
                    machine_conf[machine]['processors']+default_extraproc:
                    if proc_params!=[]: proc_param=proc_params.pop(0)
                    else: break

                    mdLog.debug("--------------------------------")
                    mdLog.info("Adding process %s to the machine %s"%(proc_param['run_id'], machine))
                    #Increasing frequency counter
                    freq=proc_param['freq']+1
                    #launching the process
                    proc_handle=launch(proc_param['dica'],
                                       proc_param['comp'],
                                       machine_conf[machine],
                                       run_id=proc_param['run_id'], 
                                       preped=True, 
                                       init_dica=init_dica, 
                                       init_comp=init_comp)
                    #Updating the stats
                    proc_param.update({'popen':proc_handle,'freq':freq, 'mconf':machine_conf[machine]})
                    proc_param.update({'startTime':datetime.now(),'checkTime':datetime.now()+checkUpTimeStep})
                    mdLog.debug("Added process %s"%proc_param['run_id'])
                    mdLog.debug("--------------------------------")
                    #append to queue
                    runnings[machine].append(proc_param)
            
        #checking if processes are finished and evaluating what to do
        for machine in runnings.keys():
            exit_flags={}
            for j, irun in enumerate(runnings[machine]):
                merror_pref='%s: %s'%(machine, error_pref)
                if irun['popen'].poll()!=None:
                    #print merror_pref+"Stderr msg for proc %s is not readable. poll indicates %s"%(irun['run_id'], irun['popen'].poll())
                    #print "Trying to get error code from file"
                    try:
                        error_file=file("proc%s/error.log"%irun['run_id']).readlines()
                        if not error_file[-1].startswith('error'):
                            mdLog.debug(machine + " Found error.log with the following in it %s"%error_file)
                        return_code=int(error_file[-1].split()[1])
                    except:
                        return_code=-999
                        mdLog.error("error.log could not be found for process %s @ %s"%(irun['run_id'], machine))
                    
                    try:
                        return_code=int(return_code)
                    except:
                        mdLog.error("Process %s returned a non standard code %s @ %s"%(irun['run_id'], return_code, machine))
                        return_code=-999
                    
                    mdLog.info("Process %s has returned %s"%(irun['run_id'], return_code))
                    if return_code==0:
                        mdLog.debug('Marking %s for cleanup'%irun['run_id'])
                        exit_flags.update({j:'cleanup'})
                        #cleanups.append()
                    else:
                        if irun['freq']>3: 
                            mdLog.error("Process %s has failed more than 3 times @ %s"%(irun['run_id'], machine))
                            exit_flags.update({j:'multifail'})
                        if return_code==2:
                            mdLog.error('%s has encountered a stale nfs handle @ %s. Resubmitting to queue'%(irun['run_id'], machine))
                            exit_flag.update({j:'resubmit'})
                        elif return_code==-999:
                            mdLog.debug("User set return code on %s"%irun['run_id'])
                            raise Exception ('shit is happening: %s'%irun)
                            exit_flags.update({j:'resubmit'})
                        else: 
                            mdLog.error("Unknown return code %s @ %s"%(return_code, machine))
                else:
                    x=0
            sexit_index=exit_flags.keys()
            sexit_index.sort()
            sexit_index.reverse()
            for iexit_flag in sexit_index:
                exit_proc=runnings[machine].pop(iexit_flag)
                if exit_flags[iexit_flag]=='cleanup':
                    proc_cleanup(exit_proc)                
                elif exit_flags[iexit_flag]=='resubmit':
                    proc_cleanup(exit_proc)
                    exit_proc.update({'popen':None})
                    proc_params.insert(0, exit_proc)
                else:
                    mdLog.critical("%s had a critical error. It won't be resubmitted or cleaned up"%exit_proc['run_id'])
        for machine in runnings.keys():
            print "%s: %02d"%(machine, len(runnings[machine]))
        lines=len(runnings.keys())+1
        if proc_params!=[]:
            mdLog.info("Next process to launch is %s"%proc_params[0]['run_id'])
            print "Next process to launch is %s/%s"%(proc_params[0]['run_id'], len(dicas))
            sys.stdout.write("\x1b[%sA"%lines)
        else:
            print "Launched all parameter sets, Still %s items in queue"%\
            sum([len(item) for item in runnings.values()])
            sys.stdout.write("\x1b[%sA"%lines)
            mdLog.info("Launched all parameter sets, Still %s items in queue"%\
            sum([len(item) for item in runnings.values()]))
            for machine in runnings.keys():
                mdLog.info("%s: %s"%(machine, len(runnings[machine])))
        time.sleep(0.2)
    sys.stdout.write('\x1b[%sB'%lines)
    print "The directory %s took %s minutes "%(getcwd(), (datetime.now()-start).seconds/60.0)
    mdLog.info("The directory %s took %s minutes "%(getcwd(), (datetime.now()-start).seconds/60.0))
    chdir('../')
            
        
def launch2(runDir,machine,firstTime=True):
    ficaScript=machine['fica_script']
    if machine_conf['type']=='local':
        raise Exception('not implemented yet.')
        return subprocess.Popen(proc_cmd, cwd=proc_dir, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    elif machine_conf['type']=='remote':
        #This only works if the path structure is the same on all machines
        #you can include a basePath
        procCMD=ficaScript.split()+[runDir]
        print procCMD
        return subprocess.Popen(procCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
def launch(dica_params, comp_params,machine_conf,run_id=0, preped=False, init_dica=None, init_comp=None):
    #if backup is 0  search the backupdir and make a backup with the last number, if backup = -1 don't backup at all, if back up is 
    #Reading in initial condition
    fica_script=machine_conf['fica_script']
    if not preped:
        if init_dica==None: dica=dicafile('../'+init_cond+'dica.dat').read_data()
        else: dica=init_dica
        if init_comp==None: comp=compfile('../'+init_cond+'comp.ind').read_data()
        else: comp=init_comp
        #Changing the parameters parsed to the launcher
        dica.update(dica_params)
        comp=setAbundances(comp, comp_params)
        proc_dir=prep_proc('proc%s'%run_id)
        dicafile(proc_dir+'dica.dat', 'w').write_data(dica)
        compfile(proc_dir+'comp.ind', 'w').write_data(comp)
        #doDecay(proc_dir)
    else:
        proc_dir='proc%s/'%run_id
    if machine_conf['type']=='local':
        proc_cmd=fica_script.split()
        return subprocess.Popen(proc_cmd, cwd=proc_dir, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
    elif machine_conf['type']=='remote':
        pathlist=getcwd().split('/')
        sub_path='/'.join(pathlist[pathlist.index('sn_rad_trans')+1:])
        proc_cmd=fica_script+' '+\
            path.join(machine_conf['path'], sub_path, proc_dir)
        #print "this is the launch cmd %s"%proc_cmd
        return subprocess.Popen(proc_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #return spect

def prep_proc(directory=None, wait=True):
    #print "---------------------------------------"
    mdLog.info("Preparing simulation directory")
    if directory==None:
        lastproc=find_last('proc')
        directory='proc%s'%(lastproc)
        #print "Chose number %s"%(lastproc+1)
    mdLog.info("Creating directory %s"%directory)
    system('mkdir %s'%directory)
    mdLog.info("symlinking files ... at %s"%getcwd())
    if conf_dir.startswith('/User'):
        syscmd='rsh moron /home/wkerzend/scripts/symlink.sh %s'%(path.join('/home/wkerzend/wdata/', '/'.join(getcwd().split('/')[-6:]) ,directory)+'/')
        print syscmd
        system(syscmd)
    else:
        system('ln -s %s %s/'%(path.join(conf_dir, '*'), directory))
    #for ifile in glob(conf_dir+'*'):
    #    symlink(ifile, directory+'/'+path.basename(ifile))
    #print "---------------------------------------"
    if wait:
        mdLog.info("waiting for symlinks to finish (1s)")
        time.sleep(1)
    return directory+'/'
    
    
def proc_cleanup(run_dict):
    run_id=run_dict['run_id']
    mdLog.debug("Cleaning up process with ID %s"%run_id)
    proc_dir=path.join('proc%s'%run_id, '')
    for ifile in backup_files:
        #print "moving %s to %s at path %s"%(proc_dir+ifile,'./'+ifile+backup_suffix+'%04d'%run_id,  getcwd())
        #backing up the undecayed version
        if ifile=='comp.ind':
            compfile('./'+ifile+backup_suffix+'%04d'%run_id, 'w').write_data(run_dict['comp'])
            continue
        try:
            copy(proc_dir+ifile,'./'+ifile+backup_suffix+'%04d'%run_id)
        except:
            mdLog.error("Problem moving file %s with run_id %s at path %s"\
                %(ifile, run_id, getcwd()))
    if del_proc:
        try:
            rmtree(proc_dir)
        except:
            mdLog.error("Can't remove proc_dir %s at path %s"%(proc_dir, getcwd()))
    mdLog.info("Cleaned up process with ID %s"%run_id)


def backup_proc(proc_no, backup_path):
    proc_dir='proc%s/'%proc_no
    num=proc_no
    for ifile in backup_files:
        while True:
            try:
                move(proc_dir+ifile,backup_path+ifile+backup_suffix+'%04d'%num)
                break
            except:
                exval=sys.exc_value
                print "couldn't move spct.dat @ %s"%proc_no
                file('error.log','a').write('The error occured %s' %exval)
