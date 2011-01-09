import numpy as np
def updateAbundances(comp,norm=True,normMode='all'):
    print notImplementedYet    


def setCONe(comp):
    comp=comp.copy()
    sumCONe=np.sum([comp['C'],comp['O'],comp['Ne']])
    comp['C']=0.01*sumCONe
    comp['O']=0.99*sumCONe
    comp['Ne']=0.0
    return comp

def setTiCr(comp):
    comp=comp.copy()
    sumTiCr=np.sum([comp['Ti'],comp['Cr']])
    comp['Ti']=0.5*sumTiCr
    comp['Cr']=0.5*sumTiCr
    return comp

def normAbundances(comp):
    comp=comp.copy()
    sumAbund=np.sum(comp.values())
    for key in comp.keys():
        comp[key]/=sumAbund
    return comp
def normAbundancesToO(comp):
    comp=comp.copy()
    for elem in comp:
        if elem!='O':
            oldElemAbundance=comp[elem]
            deltaElemAbundance=oldElemAbundance-elemAbundances[elem]
            comp['O']+=deltaElemAbundance
            if comp['O']<0: raise Exception ("Error: negative oxygen abundance")
            comp[elem]=elemAbundances[elem]            
    return comp
    print "Not yet implemented"

def setNiDecay(comp,t):
    t1=6.1
    t2=77.27    
    comp=comp.copy()
    initNiAbund=comp['Ni']
    decayNiAbund=initNiAbund*(pow(2, -(t/t1)))
    decayCoAbund=(t2/(t1-t2))*initNiAbund*(pow(2, -t/t1)-pow(2, -t/t2))
    decayFeAbund=initNiAbund-decayCoAbund-decayNiAbund
    comp['Ni']=decayNiAbund
    comp['Co']=decayCoAbund
    comp['Fe']+=decayFeAbund
    return comp

def calcNiDecay(initNi,t):
    lambdaNi=6.1
    lambdaCo=77.72
    ni=initNi*pow(2,-(t/lambdaNi))
    co=(lambdaCo/(lambdaNi-lambdaCo))*initNi*(pow(2,-t/lambdaNi)-pow(2,-t/lambdaCo))
    fe=initNi-ni-co
    return ni,co,fe
def calcNi0(decayNi,t):
    lambdaNi=6.1
    initNi=decayNi/pow(2,-(t/lambdaNi))
    retrun
def calcConstFe(initNi,constFe,t):
    ni,co,fe=calcNiDecaY(initNi,t)
    return constFe-fe