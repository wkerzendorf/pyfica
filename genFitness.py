import numpy as np
import config
from scipy import ndimage
fitConfig=config.getCurFitConfig()
uvWeight=0.5
irWeight=0.02
irThresh=fitConfig['threshold']['ir']
origSpec=config.getOrigSpec(preProcess=True)


contFactor=0.7
lineFactor=1.5
wFactor=1.

def fitInvSSE(model):
    return 1/np.sum(model.subSpec.y**2)
def fitInvSSEW(model):
    fitness=1/np.sum(model.subSpec[:6500.].y**2)
    w=model['w']
    wFactor=np.exp(-(abs(w-0.5)**5)/0.00005)
    fitness*=wFactor
    return fitness
def fitInvSSEW2(model):
    uvRange=4000-model.subSpec.x.min()
    opticalRange=6500-4000
    irRange=model.subSpec.x.max()-6500
    uv=np.sum(model.subSpec[:4000.].y**2)
    optical=np.sum(model.subSpec[4000.:6300.].y**2)
    ir=np.sum(model.subSpec[6300.:].y**2)
    w=model['w']
    wFactor=np.exp(-(abs(w-0.5)**5)/0.00005)
    fitness=(1/(1*optical))
    return fitness

def fitInvSSEW3(model):
    w=model['w']
    wFactor=wFilter(w)
    stretch=30
    sFilter=specFilter(model.subSpec.x,6300,4000,uvWeight,irWeight,stretch=30)
    squareSpec=model.subSpec.y**2
    sFitness=1/np.sum(squareSpec*sFilter)
    return sFitness*w*wFactor


def fitInvSSEW4(model):
    w=model['w']
    wFactor=wFilter(w)
    contOpticalFit=fitOpticalContinuum(model.aSpec)
    contOpticalDiff=contOpticalFit-contOrigOptical
    
    fixSubSpec,contIRDiff=fixIR(model)
    contFitness=1/(np.sum(contOpticalDiff.y**2)+(1/irWeight)*np.sum(contIRDiff.y**2))
    fixSubSpec[:irThresh]-=contOpticalDiff
    lineFitness=1/np.sum(fixSubSpec.y**2)
    
    return contFitness*lineFitness*w*wFactor

def fitInvSSEW5(model):
    w=model['w']
    wFactor=wFilter(w)
    contOpticalFit=fitOpticalContinuum(model.aSpec)
    contIRFit=fitIR(model)

    lineFit=model.aSpec[:]
    lineFit[:irThresh]/=contOpticalFit
    lineFit[irThresh:]/=contIRFit
    #
    lineOrig[:irThresh]/=contOpticalOrig
    lineOrig[irThresh:]/=contIROrig
    
    lineFitness=1/np.sum(((lineFit/lineOrig).y-1)**2)

    contIRDiff=contIRFit-contIROrig
    contOpticalDiff=contOpticalFit-contOpticalOrig
    contFitness=1/(np.sum(contOpticalDiff.y**2)+irWeight*np.sum(contIRDiff.y**2))

def fitInvSSEW5(model):
    w=model['w']
    wFactor=wFilter(w)
    contOpticalFit=fitOpticalContinuum(model.aSpec)
    contIRFit=fitIRContinuum(model.aSpec)

    lineFit=model.aSpec[:]
    lineFit[:irThresh]/=contOpticalFit
    lineFit[irThresh:]/=contIRFit
    #
    lineOrig[:irThresh]/=contOpticalOrig
    lineOrig[irThresh:]/=contIROrig
    
    lineFitness=1/np.sum(((lineFit/lineOrig).y-1)**2)

    contIRDiff=contIRFit-contIROrig
    contOpticalDiff=contOpticalFit-contOpticalOrig
    contFitness=1/(np.sum(contOpticalDiff.y**2)+irWeight*np.sum(contIRDiff.y**2))
    
    return (contFitness**contFactor)*(lineFitness**lineFactor)*(w**wFactor)
    
def fitInvSSEW6(model):
    #fitness function with relative linestrengths,
    aSpec=model.aSpec.interpolate(xref=origSpec.x)
    w=model['w']
    wFactor=wFilter(w)
    
    contOpticalFit=fitOpticalContinuum(aSpec)
    contIRFit=fitIRContinuum(aSpec)

    lineFit=aSpec[:]
    lineFit[:irThresh]/=contOpticalFit
    lineFit[irThresh:]/=contIRFit
    
    lineFitness=1/np.sum((((lineFit-lineOrig)/(lineOrig+lineFit)).y)**2)

    contIRDiff=(contIRFit-contIROrig)/(contIRFit+contIROrig)
    contOpticalDiff=(contOpticalFit-contOpticalOrig)/(contOpticalFit+contOpticalOrig)
    contFitness=1/np.sum(contOpticalDiff.y**2)+irWeight/np.sum(contIRDiff.y**2)
    
    return (contFitness**contFactor)*(lineFitness**lineFactor)*(w**wFactor)
    
def wFilter(w):
    stretch=0.01
    upLimit=0.65
    lowLimit=0.4
    filter=(1/(np.exp((w-upLimit)/stretch)+1))*(1/(np.exp((-w+lowLimit)/stretch)+1))
    return filter

def fitOpticalContinuum(spec):
    return spec[:irThresh].smoothMax(400).fitContinuum(lsigma=2,iter=4,func='poly5')
def fitIRContinuum(spec):
    return spec[irThresh:].smoothg(50).fitContinuum(lsigma=1,iter=4,func='poly3')
    
def fixIR(model):
    smoothKernelSize=50
    
    irFit=model.aSpec[irThresh:].interpolate(xref=contIROrig.x)
    
    contFitIR=irFit.smoothg(50).fitContinuum(lsigma=1,iter=4,func='poly3')
    contDiff=contFitIR-contIROrig
    newSubSpec=model.subSpec[:]
    newSubSpec[irThresh:]-=contDiff
    return newSubSpec,contDiff

def fitIR(model):
    irFit=model.aSpec[irThresh:].interpolate(xref=contIROrig.x)
    contFitIR=irFit.smoothg(50).fitContinuum(lsigma=1,iter=4,func='poly3')
    return contFitIR
def diracFermiDown(x,loc,stretch,before,after):
    return ((before-after)/(np.exp((x-loc)/stretch)+1))+after
    
def diracFermiUp(x,loc,stretch,before,after):
    return ((before-after)/(np.exp((-x+loc)/stretch)+1))+after
    
def specFilter(x,loc1,loc2,amp1,amp2,stretch):
    return diracFermiDown(x,loc1,stretch,1,amp1)*diracFermiUp(x,loc2,stretch,1,amp2)
    
    



def fitInvSSE(model):
    return 1/numpy.sum(model.subSpec.y**2)

def fitInvSSEAdjustIR(modelGrid):
    fitnessPreIR=np.array([np.sum(item[:6500.].y**2) for item in modelGrid['subspec']])
    adjustIR=[]
    for aspec,subspec in zip(modelGrid['aspec'],modelGrid['subspec']):
        adjustedContin=subspec[6500.:].fitContinuum(func='poly5')
        adjustIR.append(subspec.interpolate(xref=adjustedContin.x)-adjustedContin)
    fitnessIR=np.array([np.sum(item.y**2) for item in adjustIR])
    fitnessIROrig=np.array([np.sum(item[6500.:].y**2) for item in modelGrid['subspec']])
    fitness=1/(fitnessPreIR+0.5*fitnessIR+0.5*fitnessIROrig)
    w=modelGrid['w']
    fitness*=-4*(w**2)+4*w
    return fitness


def fitRelativInvSSEAdjustIR(modelGrid):
    fitnessPreIR=np.array([np.sum((item[:6500.].y-1)**2) for item in modelGrid['divspec']])
    adjustIR=[]
    for aspec,divspec in zip(modelGrid['aspec'],modelGrid['divspec']):
        adjustedContin=divspec[6500.:].fitContinuum(func='poly5')
        adjustIR.append(divspec.interpolate(xref=adjustedContin.x)/adjustedContin)
    fitnessIR=np.array([np.sum((item.y-1)**2) for item in adjustIR])
    fitnessIROrig=np.array([np.sum((item[6500.:].y-1)**2) for item in modelGrid['subspec']])
    return 1/(fitnessPreIR+0.5*fitnessIR+0.5*fitnessIROrig+0.2*(modelGrid['w']/0.5-1)**2)
    
def fitSSE(modelGrid):
    fitness=np.array([np.trapz(item.y**2,item.x)/(item.x[-1]-item.x[0]) for item in modelGrid['subspec']])
    return fitness
        
def fitIntegral(modelGrid):
    fitness=np.array([1/abs(sum(item.y)) for item in modelGrid['subspec']])
    return fitness/sum(fitness)
def fitW(modelGrid):
    wValue=1/np.abs(modelGrid['w']-0.5)
    return wValue/sum(wValue)
def fitSlope(modelGrid):
    slopeFunc=np.vectorize(fit.getIntSlope)
    slopeValue=1/np.abs(slopeFunc(modelGrid['subspec']))
    return slopeValue/np.sum(slopeValue)
def fitUV(modelGrid):
    UVFunc=np.vectorize(fit.getUVInt)
    UVValue=1/np.abs(UVFunc(modelGrid['subspec']))
    return UVValue/sum(UVValue)
def fitLinComb(modelGrid,weights=[1,3,1,1]):
    fitIntegralValue=fitIntegral(modelGrid)
    fitWValue=fitW(modelGrid)
    fitUVValue=fitUV(modelGrid)
    fitSlopeValue=fitSlope(modelGrid)
    #linear combination
    fitness=weights[0]*fitIntegralValue+weights[1]*fitWValue+weights[2]*fitUVValue+weights[3]*fitSlopeValue
    return fitness/sum(fitness)
def fitNNetwork(modelGrid):
    inputs=neuralDalek.getInputFromModel(modelGrid)
    net=loadnet('curNet.net')
    outputs=np.array(net(inputs)).flatten()
    fitness=1/outputs
    return fitness/sum(fitness)
    
def fitnessScale(fitness,Cmult=2.):
    fmin=fitness.min()
    fmean=fitness.mean()
    fmax=fitness.max()
    #Trying to find scaling factor m and t
    m=(Cmult*fmax-fmean)/(fmax-fmean)
    t=fmean-m*fmean
    if fmin*m+t>0: return fitness*m+t
    else: #fmin goes negative with chosen scaling
        print "Scaling to avoid negativity"
        m=(fmean)/(fmean-fmin)
        t=fmean-m*fmean
    return fitness*m+t
        
fitFunc=fitInvSSEW6

contIROrig=fitIRContinuum(origSpec)
contOpticalOrig=fitOpticalContinuum(origSpec)

#Calculating line fitnesses:
lineOrig=origSpec[:]
lineOrig[:irThresh]/=contOpticalOrig
lineOrig[irThresh:]/=contIROrig

