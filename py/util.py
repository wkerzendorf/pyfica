from bisect import bisect
from math import pi, log
from numpy import argmax, concatenate, diff, array, mean, std
from glob import glob
from os import path
from scipy import polyfit, ndimage, interpolate
import sys
import numpy as np
freq_bins=[[3500,3600],[3940,4910],[6000,6500], [8500, 10000], [6000,6500], [8500, 10000]]
funcDict={'xisquared':lambda data: mean((data[:,1]-1)**2), 
                'mean':lambda data: mean(data[:,1]-1)}
def playwave(audiofile):
    import pyaudio
    import wave
    import sys

    chunk = 1024
    
    wf = wave.open(audiofile, 'rb')
    
    p = pyaudio.PyAudio()
    
# open stream
    stream = p.open(format =
                    p.get_format_from_width(wf.getsampwidth()),
                    channels = wf.getnchannels(),
                    rate = wf.getframerate(),
                    output = True)

# read data
    data = wf.readframes(chunk)

# play stream
    while data != '':
        stream.write(data)
        data = wf.readframes(chunk)

    stream.close()
    p.terminate()
    
class keyhandler:
    def __init__(self,data, linecat):
        self.linecat=linecat
        self.data=data
        self.calc_velo=""
        self.calc_gauss=""
    def __call__(self,event):
        import numpy
        if event.key.lower()=='w':
            print event.xdata, event.ydata
        elif event.key.lower()=='r':
            linecat=self.linecat
            lindex=bisect(zip(*linecat)[1],event.xdata-80)
            hindex=bisect(zip(*linecat)[1],event.xdata+80)
            print
            print "------------------"
            print "The Reference wavelength is %s"%event.xdata

            for irow in linecat[lindex:hindex]:
                print ' '.join(map(str,irow))
        elif event.key.lower()=='v':
            if self.calc_velo=="":
                print "Please select the second line py pressing v again to calculate velocity"
                self.calc_velo=event.xdata
            else: 
                print "velocity difference is %s"%(3e5*(self.calc_velo-event.xdata)/event.xdata)
                self.calc_velo=""
        elif event.key.lower()=='k':
            if self.calc_gauss=="":
                print "Please select the second point py pressing k again to fit a gaussian"
                self.calc_gauss=event.xdata
            else:
                id1=bisect(self.data[:, 0], self.calc_gauss)
                id2=bisect(self.data[:, 0], event.xdata)
                self.calc_gauss=""
                fitdata=self.data[id1:id2]
                self.fitdata=fitdata
                margin=round(len(fitdata)/100+0.5)
                print fitdata[:margin, 1] + fitdata[-margin:, 1]
                norm=numpy.mean(fitdata[:margin, 1] + fitdata[-margin:, 1])
                fitdata[:, 1]/=norm
                fitdata[:, 1]-=1
                if abs(numpy.min(fitdata[:, 1])-1) > abs(numpy.max(fitdata[:, 1])-1):
                    #absorption
                    gamp=numpy.min(fitdata[:, 1])
                    gmean=fitdata[numpy.argmin(fitdata[:, 1]), 0]
                    
                else:
                    #emission
                    gamp=numpy.max(fitdata[:, 1])
                    gmean=fitdata[numpy.argmax(fitdata[:, 1]), 0]
                print "Estimating mean %s and amp %s "%(gmean, gamp)
                par_out, fitted_data=fit_gauss(self.data,  est=[gamp, gmean, 1])
                print "Fitted a Gauss function: Amplitude= %s Mean= %s Std=%s"%(par_out[0][0], par_out[0][1], par_out[0][2])
                
        else:
            print "Not implemented yet"
            

def fit_gauss(data, est=[1, 0, 1]):
    import scipy,  scipy.optimize
    import numpy
    par=numpy.array([0, 0, 0])
    par[0]=est[0] #amplitude
    par[1]= est[1] # mean
    par[2]= est[2] #sigma
    gauss_func = lambda p,  x : p[0]*scipy.exp(-(x-p[1])**2/(2.0*p[2]**2))
    error_func= lambda p,  x,  y: gauss_func(p, x)-y
    par_out = scipy.optimize.leastsq(error_func, par.copy(), args=(data[:, 0], data[:, 1]))
    if len(par_out[0])==3: output=gauss_func(par_out[0], data[:, 0])
    else: output=[]
    return par_out, output
    
def savehist(logmsg=""):
    import shutil, os, glob
    print "Saving current setup to hist"
    if not os.path.exists('./history/'):
        os.mkdir('./history')
        cur_entry=1
        log=file('./history/history.txt', 'w')
    else:
        if os.path.exists('./history/history.txt'): log=file('./history/history.txt', 'a')
        else: log=file('./history/history.txt', 'w')
        filelist=glob.glob('./history/*.hist.????')
        cur_entry=max(int(item[-4:]) for item in filelist)+1
    print "Creating new history item with number %04d"%cur_entry
    suffix=".hist.%04d"%cur_entry
    bup_files=['comp.ind', 'dica.dat', 'spct.dat']
    for ifile in bup_files:
        shutil.copy(ifile, './history/'+ifile+suffix)
    log.write("%04d log:"%(cur_entry)+logmsg+"\n")
    log.close()

def loadhist(number):
    import shutil, os
    if not os.path.exists('./history/'): raise Exception("The history directory does not exists")
    files=['comp.ind', 'dica.dat']
    suffix=".hist.%04d"%number
    for ifile in files:
        shutil.copy('./history/'+ifile+suffix, './'+ifile)

def plothist(number):
    print "test"

    
def setAbundances(comp, elemAbundances, t=None):
    for elem in elemAbundances:
        if elem!='O':
            oldElemAbundance=comp[elem]
            deltaElemAbundance=oldElemAbundance-elemAbundances[elem]
            comp['O']+=deltaElemAbundance
            if comp['O']<0: raise Exception ("Error: negative oxygen abundance")
            comp[elem]=elemAbundances[elem]            
    return comp


def calcTemp(t, lum=None, vph=None, T=None):
    day2sec=3600*24
    sigma=5.6705e-8
    L_sun=3.846e26
    if lum and vph:
        L_sn=(10**(lum))*L_sun
        t=t*day2sec
        r=vph*1e3*t
        T=(L_sn/(4*pi*sigma*(r**2)))**0.25
        return T
    elif lum and T:
        L_sn=(10**(lum))*L_sun
        t=t*day2sec
        v_ph=((L_sn/(4*pi*sigma*(T**4)*(t**2) ))**0.5)*1e-3
        print "Photospheric velocity is %s km/s"%v_ph
        return v_ph
    elif vph and T:
        t=t*day2sec
        r=vph*1e3*t
        L_sn=4*pi*sigma*(r**2)*(T**4)
        log_lbol=log(L_sn/L_sun, 10)
        print "Log(L/L_sol) is %s"%log_lbol
        return log_lbol
    else:
        raise Exception ("please use two of the three variables")
    
def findLast(pref, filePath=''):
    files=glob(path.join(filePath, pref+'*'))
    if files==[]:
        return 1
    last=argmax([int(ifile[len(filePath+pref):]) for ifile in files])
    return files[last][len(filePath+pref):]
    
def dictBuilder(edict,items):
	if items==[]: 
		return edict
	key=items.pop(0)
	if edict.has_key(key): 
		dictBuilder(edict[key],items)
	else:
		edict[key]={}
		dictBuilder(edict[key],items)

def median_spectra(spectra,smooth_size):
    pix_window=bisect(spectra[:,0],spectra[0,0]+smooth_size)
    
    new_spectra=ndimage.median_filter(spectra[:, 1], pix_window)
    return array(zip(spectra[:,0],new_spectra))
def gaussian_spectra(spectra,smooth_size):
    pix_window=bisect(spectra[:,0],spectra[0,0]+smooth_size)
    
    new_spectra=ndimage.gaussian_filter1d(spectra[:, 1], pix_window)
    return array(zip(spectra[:,0],new_spectra))
def diffSpectrum(spectrum, order=1):
    dy=diff(spectrum[:, 1], order)
    dx=diff(spectrum[:, 0])[order-1:]
    return array(zip(spectrum[:, 0], concatenate((dy/dx, [0]*order))))

def getSpectrumScale(spectrum):
    dx=diff(spectrum[:, 0])
    if mean(dx)/1000<std(dx): raise Exception ("Problem with scale. Spectrum probably not evenly spaced")
    return mean(dx)
    
def syn_div_orig(synspec, origspec):
    #linear interpolating the synthetic spectrum on the 
    new_syn=interpolate.splev(origspec[:, 0], interpolate.splrep(synspec[:, 0], synspec[:, 1], k=1))/origspec[:, 1]
    return array(zip(origspec[:, 0], new_syn))

def getBinMerits(model, fBins=freq_bins, func='xisquared', doProc='gaussian', smoothSize=500):
    #print "Calculating LumVph Merits:"
    if type(func)==str: 
        meritFunc=funcDict[func]
    else:
        meritFunc=func
    merits=[0]*(len(model))
    if doProc=='median': med_origspect=median_spectra(model[0]['origspect'],smoothSize)
    if doProc=='gaussian': gauss_origspect=gaussian_spectra(model[0]['origspect'],smoothSize)
    for no, item in enumerate(model):
        print "At model %s"%no
        sys.stdout.write('\x1b[1A') 
        if doProc=='median': 
            med_aspect=median_spectra(item['aspect'],smoothSize)
            med_divspect=syn_div_orig(med_aspect,med_origspect)
            merits[no]=multi_fit(med_divspect, fBins, [meritFunc])[0]
        if doProc=='gaussian': 
            gauss_aspect=median_spectra(item['aspect'],smoothSize)
            gauss_divspect=syn_div_orig(gauss_aspect,gauss_origspect)
            merits[no]=multi_fit(gauss_divspect, fBins, [meritFunc])[0]
        else:
            merits[no]=multi_fit(item['divspect'],  fBins, [meritFunc])[0]
    return array(merits)
    
def multi_fit(data, bins, functions):
    data_max=max(data[:, 0])
    data_min=min(data[:, 0])
    id_bins=[]
    for bin in bins:
        tmp_bin=[bisect(data[:, 0], bin[0]), bisect(data[:, 0], bin[1])]
        if tmp_bin[0]!=tmp_bin[1]: id_bins.append(tmp_bin)

    #Actual fitting
    results=[]
    for ifunc in functions:
        funcbundle=[]
        for ibin in id_bins:
            funcbundle.append(ifunc(data[ibin[0]:ibin[1]]))
        results.append(funcbundle)
    return results

def findLineEdges(spectrum, point, smoothSize=50):
    gSpectrum=gaussian_spectra(spectrum, smoothSize)
    #CaCenter=bisect(gOrigspect[:, 0], [item[1] for item in model[0]['sbib']['llist'] if int(item[2])==3933 and item[3].lower()=='ca'][0])
    center=bisect(gSpectrum[:, 0], point)
    secondDiff=diffSpectrum(gSpectrum, 2)
    specScale=getSpectrumScale(spectrum)
    i=1
    lowerBound=None
    upperBound=None
    #checking 1000 angstrim left and right of the Center
    for i in range(1, int(1000/specScale)):
        if secondDiff[center-i, 1]<0 and lowerBound==None: lowerBound=center-i
        if secondDiff[center+i, 1]<0 and upperBound==None :upperBound=center+i
        if lowerBound!=None and upperBound!=None: break
    if lowerBound==None or upperBound==None: raise Exception("Couldn't find edges to feature")
    bounds=[secondDiff[lowerBound, 0], secondDiff[upperBound, 0]]
    return bounds

def makeGrid(xLimits,yLimits,cells=10):
    x=np.linspace(xLimits[0],xLimits[1],cells)
    y=np.linspace(yLimits[0],yLimits[1],cells)
    X,Y=np.meshgrid(x,y)
    return X.flatten(),Y.flatten()
