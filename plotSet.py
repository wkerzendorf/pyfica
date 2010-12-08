import matplotlib
#from matplotlib import pylab
import pylab
#from mpl_toolkits.axes_grid import AxesGrid
from glob import glob
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import cPickle
import pdb
import re
import genFitness
def plotSelNeuralModels(sources,selection,fdist=None,pdist=None,pdfName='selModels.pdf'):
    pdistHeader=['log_lbol','v_ph','C','O','Mg','Si','S','Ca','Ti','Cr','Ni0','Fe0']
    fdistHeader=['int','intUV','slope','w','t','C','O','Na','Mg','Si','S','Ca','Ti','Cr','Ni','Fe']
    modelDict={}
    bestFitDict={}
    pdf=PdfPages(pdfName)
    
    for select in selection:
        source=sources[select]
        fig=pylab.figure(1)
        ax=fig.add_subplot(111)
        modelName=source[0]
        bestFitName=bestFitName=re.sub('\.ts\d+\.','.bf.',source[0])
        modelNo=source[1]
        if not modelDict.has_key(modelName):
            print "Loading model %s"%modelName
            modelDict[modelName]=cPickle.load(file(modelName))
        if not bestFitDict.has_key(bestFitName):
            print "Loading bestFit %s"%bestFitName
            bestFitDict[bestFitName]=cPickle.load(file(bestFitName))
        model=modelDict[modelName]
        bestFit=bestFitDict[bestFitName]
        origSpec=model.origSpec
        aSpec=model['aSpec'][modelNo]
        ax.set_title('%s'%select)
        ax.plot(origSpec.x,origSpec.y,color='k')
        ax.plot(aSpec.x,aSpec.y,color='r')
        ax.plot(bestFit['aspec'].x,bestFit['aspec'].y,color='b')
        fig.savefig(pdf,format='pdf')
        fig.clf()
        plotParam(model.grid[modelNo],fig)
        if pdist!=None and fdist!=None:
            ax=fig.gca()
            pdistString=''.join(["%s %s\n"%item for item in zip(pdistHeader,pdist[select])])
            fdistString=''.join(["%s %s\n"%item for item in zip(fdistHeader,fdist[select])])
            fdistString+="fdist %s"%np.nansum(np.abs(fdist[select]))
            pdistString+="pdist %s"%np.mean(np.abs(pdist[select]))
            ax.text(0.4,0,pdistString)
            ax.text(0.8,0,fdistString)
        fig.savefig(pdf,format='pdf')
        fig.clf()
    pdf.close()
def plotParam(param,fig):
    fig.clf()
    ax=fig.add_subplot(111)
    ax.text(0.01,0,''.join(['%s %s\n'%(item,param[item]) for item in ['log_lbol','vph','t','C','O','Na','Mg','Si','S','Ca','Ti','Cr','Ni','Fe']]))
    
def plotElementModelGrid(element,modelGrid,suffix=""):
    pdf=PdfPages('%splot%s.pdf'%(element,suffix))
    origSpec=modelGrid.origSpec
    for aspec,abund in zip(modelGrid['aspec'],modelGrid[element]):
        fig=pylab.figure(1)
        ax=fig.add_subplot(111)
        ax.plot(origSpec.x,origSpec.y,color='k')
        ax.plot(aspec.x,aspec.y,color='r',label='%s=%s'%(element,abund))
        fig.savefig(pdf,format='pdf')
        fig.clf()
    pdf.close()
def plotModelGrid(modelGrid,id):
    fig=pylab.figure(1)
    fig.clf()
    ax=fig.add_subplot(111)
    ax.plot(modelGrid.origSpec.x,modelGrid.origSpec.y,'k')
    ax.plot(modelGrid.grid[id]['aspect'].x,modelGrid.grid[id]['aspect'].y,'r')
def plotGridAspect(specGrid,shape=None,origSpect=None):
    fig=pylab.gcf()
    if shape==None: shape=specGrid.shape
    grid = AxesGrid(fig, 111,
                nrows_ncols = shape,
                axes_pad=0.1,
                share_all=True,
                aspect=False)
    for i,ispec in enumerate(specGrid.flatten()):
        grid[i].plot(origSpect.x,origSpect.y,linewidth=3)
        grid[i].plot(ispec.x,ispec.y)
def overPlotGridAspect(specGrid,labels=None,origSpect=None,colors=None,cmap=pylab.cm.hot):
    if origSpect!=None: pylab.plot(origSpect.x,origSpect.y,color='blue',linewidth=3)
    if colors!=None:
        normfunc=matplotlib.colors.normalize()
        colors=cmap(normfunc(colors.flatten()))
    for i,ispec in enumerate(specGrid.flatten()):
        specplot=pylab.plot(ispec.x,ispec.y)
        if labels!=None: specplot[0].set_label(labels.flatten()[i])
        if colors!=None:
            
            specplot[0].set_color(colors[i])
    if labels!=None: pylab.legend()
def overPlotGridAspectSelect(specGrid,mask,labels=None,origSpect=None,colors=None,cmap=pylab.cm.hot):
    newSpecGrid=specGrid[mask]
    if colors!=None: newColors=colors[mask]
    else: newColors=None
    if labels!=None: labels=labels[mask]
    overPlotGridAspect(newSpecGrid,labels=labels,origSpect=origSpect,colors=newColors,cmap=cmap)
def plotGridAspectSelect(specGrid,mask,origSpect=None):
    fig=pylab.gcf()
    newSpecGrid=specGrid[mask]
    gridNo=len(newSpecGrid)
    gridShape=int(np.floor(sqrt(gridNo))),int(np.ceil(sqrt(gridNo)))
    plotGridAspect(newSpecGrid,shape=gridShape,origSpect=origSpect)

def createLabels(*args,**kwargs):
    if not all([item.shape==args[0].shape for item in args]): raise Exception('The shapes are not the same for all objects')
#   labelStr='lum=%s vph=%s'%(fmt,fmt)
    labels=[]
    #return zip(*[item.flatten() for item in args])
    for item in zip(*[item.flatten() for item in args]):
        if kwargs.has_key('fmt'): labels.append(kwargs['fmt']%item)
        else: labels.append(' '.join(map(str,item)))
    
    return np.array(labels).reshape(args[0].shape)

def genPlotSetAspec(modelGrid,pdfName='set.pdf'):
    pdf=PdfPages(pdfName)
    fig=pylab.figure(1)
    fig.clf()
    origSpec=modelGrid.origSpec
    for aspec in modelGrid['aspec']:
        ax=fig.add_subplot(111)
        ax.plot(origSpec.x,origSpec.y,color='k')
        ax.plot(aspec.x,aspec.y)
        fig.savefig(pdf,format='pdf')
        fig.clf()
    pdf.close()
def genPlotSingleGeneration(modelGrid,fitness,generation='unknown',no=None,pdfName='curGeneration.pdf'):
    pdf=PdfPages(pdfName)
    fig=pylab.figure(1)
    fig.clf()
    origSpec=modelGrid.origSpec
    sortID=np.argsort(fitness)[::-1]
    if no!=None:
        sortID=sortID[:no]
    for id in sortID:
        aspec=modelGrid['aspec'][id]
        paramString=    """Lum=%.4f\n\
vph=%d\n\
Fe0=%s\n\
Ni0=%s\n\
Cr=%s\n\
Si=%s\n\
S=%s\n\
Ca=%s\n\
Mg=%s\n\
C=%s\n\
                        """%(modelGrid['lum'][id],
                           modelGrid['vph'][id],
                           modelGrid['Fe0'][id],
                           modelGrid['Ni0'][id],
                           modelGrid['Cr'][id],
                           modelGrid['Si'][id],
                           modelGrid['S'][id],
                           modelGrid['Ca'][id],
                           modelGrid['Mg'][id],
                           modelGrid['C'][id],
                           )
        curFitness=fitness[id]
        ax=fig.add_subplot(111)
        ax.plot(origSpec.x,origSpec.y,color='k')
        ax.plot(aspec.x,aspec.y,color='red')
        curXlim=ax.get_xlim()
        #plotting fit continuum
        ax.plot(modelGrid[id].contOptical.x,modelGrid[id].contOptical.y,color='r',lw=3,alpha=0.3)
        ax.plot(modelGrid[id].contIR.x,modelGrid[id].contIR.y,color='r',lw=3,alpha=0.3)
        #plotting origspec continuum
        ax.plot(genFitness.contOpticalOrig.x,genFitness.contOpticalOrig.y,color='k',lw=3,alpha=0.3)
        ax.plot(genFitness.contIROrig.x,genFitness.contIROrig.y,color='k',lw=3,alpha=0.3)
        ax.text(0.8, 0.6,paramString,
             horizontalalignment='center',
             verticalalignment='center',
             transform = ax.transAxes,
             bbox=dict(edgecolor='black',
                       facecolor='none',
                       alpha=0.5)
             )
        
        ax.set_title('Fitness=%s modelid=%s generation=%s (10 Best)'%(fitness[id],id,generation))
        fig.savefig(pdf,format='pdf')
        fig.clf()
        ax=fig.add_subplot(211)
        contOpticalDiff=(modelGrid[id].contOptical-genFitness.contOpticalOrig)/(modelGrid[id].contOptical+genFitness.contOpticalOrig)
        contIRDiff=(modelGrid[id].contIR-genFitness.contIROrig)/(modelGrid[id].contIR+genFitness.contIROrig)
        
        contOpticalDiff.y**=2
        contIRDiff.y**=2
        contOpticalFitness=np.sum(contOpticalDiff.y)
        contIRFitness=np.sum(contIRDiff.y)
        lineFit=modelGrid[id].lineCheck
        lineFitness=np.sum(lineFit.y)
        ax.set_title("contFitness=%s LineFitness=%s\n (first continuum fitness and then linefitness) w=%s"%(1/lineFitness,(1/(contOpticalFitness+genFitness.irWeight*contIRFitness)),modelGrid[id]['w']))
        ax.plot(contOpticalDiff.x,contOpticalDiff.y,color='b',label="unFitness= %s"%contOpticalFitness)
        ax.plot(contIRDiff.x,contIRDiff.y,color='r',label="unFitness= %s"%contIRFitness)
        ax.legend()
        ax.set_xlim(curXlim)
        ax=fig.add_subplot(212)
        
        
        ax.plot(lineFit.x,lineFit.y,color='k',label="unFitness=%s"%lineFitness)
        ax.legend()
        ax.set_xlim(curXlim)
        fig.savefig(pdf,format='pdf')
        fig.clf()
    if no!=None:
        for id in np.argsort(fitness)[::-1][-no:]:
            aspec=modelGrid['aspec'][id]
            paramString="""Lum=%.4f\n\
vph=%d\n\
Fe0=%s\n\
Ni0=%s\n\
Cr=%s\n\
Si=%s\n\
S=%s\n\
Ca=%s\n\
Mg=%s\n\
C=%s\n\
                        """%(modelGrid['lum'][id],
                           modelGrid['vph'][id],
                           modelGrid['Fe0'][id],
                           modelGrid['Ni0'][id],
                           modelGrid['Cr'][id],
                           modelGrid['Si'][id],
                           modelGrid['S'][id],
                           modelGrid['Ca'][id],
                           modelGrid['Mg'][id],
                           modelGrid['C'][id],
                           )
            curFitness=fitness[id]
            ax=fig.add_subplot(111)
            ax.text(0.8, 0.6,paramString,
             horizontalalignment='center',
             verticalalignment='center',
             transform = ax.transAxes,
             bbox=dict(edgecolor='black',
                       facecolor='none',
                       alpha=0.5)
             )
            ax.plot(origSpec.x,origSpec.y,color='k')
            ax.plot(aspec.x,aspec.y,color='red')
            ax.set_title('Fitness=%s modelid=%s generation=%s (10 worst)'%(fitness[id],id,generation))
            fig.savefig(pdf,format='pdf')
            fig.clf()
            ax=fig.add_subplot(211)
        contOpticalDiff=(modelGrid[id].contOptical-genFitness.contOpticalOrig)/(modelGrid[id].contOptical+genFitness.contOpticalOrig)
        contIRDiff=(modelGrid[id].contIR-genFitness.contIROrig)/(modelGrid[id].contIR+genFitness.contIROrig)
        
        contOpticalDiff.y**=2
        contIRDiff.y**=2
        contOpticalFitness=np.sum(contOpticalDiff.y)
        contIRFitness=np.sum(contIRDiff.y)
        lineFit=modelGrid[id].lineCheck
        lineFitness=np.sum(lineFit.y)
        ax.set_title("contFitness=%s LineFitness=%s\n (first continuum fitness and then linefitness) w=%s"%(1/lineFitness,(1/(contOpticalFitness+genFitness.irWeight*contIRFitness)),modelGrid[id]['w']))
        ax.plot(contOpticalDiff.x,contOpticalDiff.y,color='b',label="unFitness= %s"%contOpticalFitness)
        ax.plot(contIRDiff.x,contIRDiff.y,color='r',label="unFitness= %s"%contIRFitness)
        ax.legend()
        ax.set_xlim(curXlim)
        ax=fig.add_subplot(212)
        
        
        ax.plot(lineFit.x,lineFit.y,color='k',label="unFitness=%s"%lineFitness)
        ax.legend()
        ax.set_xlim(curXlim)
        fig.savefig(pdf,format='pdf')
        fig.clf()
    pdf.close()
def genPlotSetHist(modelGrid,pdfName='hist.pdf',params=['lum','vph','C','Mg','Si','S','Ca','Ti','Cr','Ni0','Fe0']):
    pdf=PdfPages(pdfName)
    fig=pylab.figure(1)
    fig.clf()
    for i,param in enumerate(params):
        #ax=fig.add_subplot(3,4,i+1)
        ax=fig.add_subplot(111)
        ax.hist(modelGrid[param])
        ax.set_title('param %s mean: %s std %s'%(param,np.mean(modelGrid[param]),np.std(modelGrid[param])))
        fig.savefig(pdf,format='pdf')
        fig.clf()
    pdf.close()

def genPlotSetHist(modelGrid,pdfName='hist.pdf',params=['lum','vph','C','Mg','Si','S','Ca','Ti','Cr','Ni0','Fe0']):
    pdf=PdfPages(pdfName)
    fig=pylab.figure(1)
    fig.clf()
    for i,param in enumerate(params):
        #ax=fig.add_subplot(3,4,i+1)
        ax=fig.add_subplot(111)
        ax.hist(modelGrid[param])
        ax.set_title('param %s mean: %s std %s'%(param,np.mean(modelGrid[param]),np.std(modelGrid[param])))
        fig.savefig(pdf,format='pdf')
        fig.clf()
    pdf.close()

def genPlotStatus(modelGrid,pdfName='curStatus.pdf'):
    fname=glob('*.bf.pkl')
    if len(fname)==1:
        plotBF=True
        bf=cPickle.load(file(fname[0]))
    else: plotBF=False
    pdf=PdfPages(pdfName)
    fig=pylab.figure(1)
    fig.clf()
    ax=fig.add_subplot(111)
    x=modelGrid['lum']
    y=modelGrid['vph']
    ax.plot(x,y,'b,')
    if plotBF:
        ax.axvline(bf['lum'],color='r')
        ax.axhline(bf['vph'],color='r')
    ax.set_xlabel('lum')
    ax.set_ylabel('vph')
    fig.savefig(pdf,format='pdf')
    
    fig.clf()
    ax=fig.add_subplot(221)
    x=modelGrid['Si']
    y=modelGrid['S']
    ax.plot(x,y,'b,')
    if plotBF:
        ax.axvline(bf['Si'],color='r')
        ax.axhline(bf['S'],color='r')
    ax.set_xlabel('Si')
    ax.set_ylabel('S')
    
    ax=fig.add_subplot(222)
    x=modelGrid['Fe0']
    y=modelGrid['Ni0']
    ax.plot(x,y,'b,')
    if plotBF:
        ax.axvline(bf['Fe0'],color='r')
        ax.axhline(bf['Ni0'],color='r')
    ax.set_xlabel('Fe0')
    ax.set_ylabel('Ni0')
    
    ax=fig.add_subplot(223)
    x=modelGrid['Ti']
    y=modelGrid['Cr']
    ax.plot(x,y,'b,')
    if plotBF:
        ax.axvline(bf['Ti'],color='r')
        ax.axhline(bf['Cr'],color='r')
    ax.set_xlabel('Ti')
    ax.set_ylabel('Cr')
    
    ax=fig.add_subplot(224)
    x=modelGrid['C']
    y=modelGrid['O']
    ax.plot(x,y,'b,')
    if plotBF:
        ax.axvline(bf['C'],color='r')
        ax.axhline(bf['O'],color='r')
    ax.set_xlabel('C')
    ax.set_ylabel('O')
    fig.savefig(pdf,format='pdf')
    fig.clf()
    pdf.close()
def genPlotHistogram(modelGrid,fitness=None,pdfName='curHistogram.pdf'):
    fname=glob('*.bf.pkl')
    pdf=PdfPages(pdfName)
    fig=pylab.figure(1)
    fig.clf()
    if len(fname)==1:
        plotBF=True
        bf=cPickle.load(file(fname[0]))
    else:
        plotBF=False
    for i,param in enumerate(['log_lbol','v_ph','C','O','Mg','Si','S','Ca','Ti','Cr','Ni0','Fe0']):
        ax=fig.add_subplot(3,1,i%3+1)
        if plotBF:
            ax.hist(modelGrid[param],bins=20,label="%s BestFit=%s"%(param,bf[param]))
            ax.axvline(bf[param],color='r')
        else:
            ax.hist(modelGrid[param],bins=20,label=param)
        if fitness!=None:
            ax.axvline(modelGrid[param][np.argsort(fitness)[-1]],color='black',lw=3,label='current Best=%s'%modelGrid[param][np.argsort(fitness)[-1]])
        ax.set_xlabel(param)
        ax.legend(loc=0)
        #ax.semilogy()
        if i%3==2:
            fig.savefig(pdf,format='pdf')
            fig.clf()
    pdf.close()
        
def genPlotGenVSFitness(x,y,y2,yerr,logPlot=True,outName='gen_vs_fitness.png'):
    fig=pylab.figure(1)
    fig.clf()
    ax=fig.add_subplot(111)
    ax.errorbar(x,y,yerr,marker='x')
    ax.plot(x,y2,'r-')
    if logPlot==True: ax.semilogy()
    fig.savefig(outName)