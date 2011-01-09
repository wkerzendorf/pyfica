#A script file to explore the dataset
exploreLumVph=True
if exploreLumVph:
    readIn=False
    if readIn:
        #models=read.readModels()
        z=models[0]['dica']['z']
        origspec=spectrum('origspect.dat')
        origspec=origspec.shiftVelocity(z=z)
        modelgrid=lumVphGrid(models,(10,10))
        specgrid=modelgrid.getAspect().data
        lums=modelgrid.getLums().data
        vphs=modelgrid.getVphs().data
        integral=getIntMeritLumVph(modelgrid)
        lvlabel=modelgrid.getLumsVphsLabels()
        lvilabel=createLabels(lvlabel,integral,fmt='%s %.4f')
        
def getIntMeritLumVph(modelGrid):
    integral=modelGrid.getIntTrapz(xref=origspec.x)
    return abs(integral/origspec.intTrapz()-1)
def getXiMeritLumVph(modelGrid):
    print
def imshowLumVph(model):
    imshowGrid(model,lums[0],vphs[:,0])
def imshowGrid(modelGrid,x,y):
    imshow(modelGrid,aspect='auto',extent=(min(x),max(x),min(y),max(y)),interpolation='nearest')
