import scipy as sp
import numpy as np
import util
def getUVInt(subSpec,norm=True):
    if norm: return (subSpec[:3750.]).intTrapz()/(3750.-subSpec.x.min())
    return (subSpec[:3750.]).intTrapz()
def getUVIntComp(subSpec,norm=True):
    if norm: return getUVInt(subSpec,norm)-(subSpec[3950.:6350.]).intTrapz()/(6350.-3950.)
    return getUVInt(subSpec,norm)-(subSpec[3950.:6350.]).intTrapz()
def getIntSlope(subSpec,bins=2):
    wlBinMiddle,integrals=getDiffIntBin(subSpec,bins=bins)
    return sp.polyfit(wlBinMiddle,integrals,1)[0]
def getChiSquared(subSpec):
    return np.sum(subSpec.y**2)
def getModChiSquared(subSpec):
    cont=subSpec.fitContinuum(func='poly5',iter=5)
    return np.sum((subSpec/cont).y**2)
def getMAD(aSpec,origSpec):
    aspecInterp=aSpec.interpolate(xref=origSpec.x)
    medOrigSpec=np.median(origSpec.y)
    return np.median(np.abs(aspecInterp.y-medOrigSpec))
    
def getSlope(subSpec):
    return sp.polyfit(subSpec[:7e3].x,subSpec[:7e3].y,1)[0]
def getInt(spec,norm=True,lower=None,upper=7000.):
    if norm:
        specInt=(spec[slice(lower,upper)]).intTrapz()
        if lower==None or lower < spec.x.min(): lower=spec.x.min()
        if upper==None or upper > spec.x.max(): upper=spec.x.max()
        return specInt/(upper-lower)
    return (spec[slice(lower,upper)]).intTrapz()
def getDiffIntBin(subSpec,bins,norm=True):
    integrals=[]
    limits,step=np.linspace(subSpec.x.min(),subSpec.x.max(),bins+1,retstep=True)
    wlBinMiddle=limits[:-1]+step/2.
    for i in range(limits.size-1):
        integrals.append(subSpec[slice(*limits[i:i+2])].intTrapz())
    if norm: return wlBinMiddle,np.array(integrals)/(subSpec.x.max()-subSpec.x.min())
    return wlBinMiddle,integrals

def getSubSpec(aspec,origspec):
    aspecInterp=aspec.interpolate(xref=origspec.x)
    return aspecInterp-origspec

def getAddSpec(aspec,origspec):
    aspecInterp=aspec.interpolate(xref=origspec.x)
    return aspecInterp+origspec

def setLumWithConstTemp(lum,curLum=None,vph=None,t=None,dica=None):
    if dica!=None:
        curLum=dica['log_lbol']
        vph=dica['v_ph']
        t=dica['t']
    curT=util.calcTemp(t,curLum,vph)
    return util.calcTemp(t,lum=lum,T=curT)


