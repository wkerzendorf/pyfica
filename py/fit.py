import scipy as sp
import numpy as np
import util
def getUVInt(subSpec):
    return (subSpec[:3750.]).intTrapz()
def getUVIntComp(subSpec):
    return getUVInt(subSpec)-(subSpec[3950.:6350.]).intTrapz()
def getIntSlope(subSpec,bins=2):
    wlBinMiddle,integrals=getDiffIntBin(subSpec,bins=bins)
    return sp.polyfit(wlBinMiddle,integrals,1)[0]

def getDiffIntBin(subSpec,bins):
    integrals=[]
    limits,step=np.linspace(subSpec.x.min(),subSpec.x.max(),bins+1,retstep=True)
    wlBinMiddle=limits[:-1]+step/2.
    for i in range(limits.size-1):
        integrals.append(subSpec[slice(*limits[i:i+2])].intTrapz())
    return wlBinMiddle,integrals

def getSubSpec(aspec,origspec):
    aspecInterp=aspec.interpolate(xref=origspec.x)
    return aspecInterp-origspec

def setLumWithConstTemp(lum,curLum=None,vph=None,t=None,dica=None):
    if dica!=None:
        curLum=dica['log_lbol']
        vph=dica['v_ph']
        t=dica['t']
    curT=util.calcTemp(t,curLum,vph)
    return util.calcTemp(t,lum=lum,T=curT)