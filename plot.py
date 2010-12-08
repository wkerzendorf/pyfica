import matplotlib
#from matplotlib import pylab
import pylab
from mpl_toolkits.axes_grid import AxesGrid
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

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