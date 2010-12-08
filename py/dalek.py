from scipy.integrate import trapz
from scipy import polyfit, ndimage
#from scipy.interpolate import splrep, splev
#BAD THINGS CHECK YOUR DATA FOR THESE THINGS:
# it must be evenly spaced , otherwise rewrite si fitting routine and maybe others

#Bad solution for temperature fitting. if it doesnt have 5739 SI III line in model its disregarded as bad. need to change this



seedno=250819801106

def singleFit(steps=['lumvph1', 'ca1'],  initModel=None, dalekStep=None, initStep=1):
    
    selDica=None
    selModel=None
    initSample=10
    origspect=array(datafile('origspect.dat').read_data())
    origspect[:, 0]/=(1+dicafile('init_cond/dica.dat').read_data()['z'])
    if not os.path.exists('dalek/'): os.mkdir('dalek')
    use_machines=['miner','myriad','maggot','moron','merino','minotaur','miami','munch', 'macerate']

    def initDicaComp():
        if initModel!=None:
            if initModel=='init':
                return pickle.load(file('init_cond/init_model.pkl'))
            return initModel
        elif dalekStep!=None:
            return pickle.load(file('dalek/model.pkl.step'+dalekStep))
        else:
            step=findLast('dica.dat.step', 'dalek/')
            return readModel(step, origspect, 'dalek/', '.step')
    def updateModel(model, origspect, bounds):
        print
    def writeModel2Dalek(model, step, runDir):
        if type(step)==int: step="%03d"%step

        dicafile('dalek/dica.dat.step'+step, 'w').write_data(model['dica'])
        compfile('dalek/comp.ind.step'+step, 'w').write_data(model['comp'])
        datafile('dalek/divspect.dat.step'+step, 'w').write_data(model['divspect'])
        shutil.copy(runDir+'sbib.dat.bak%s'%model['id'], 'dalek/sbib.dat.step'+step)
        shutil.copy(runDir+'spct.dat.bak%s'%model['id'], 'dalek/spct.dat.step'+step)
        pickle.dump(model, file('dalek/model.pkl.step'+step, 'w'))
    
    def fitS(selModel, step, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='srun%s/'%step
        dicas, comps=calc_grid_elem('S', [selComp['Si']/5,selComp['Si']/3])
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        model=readModels(origspect, run_dir, backup=True)
        #getting the limits to the Ca line by differentiating 2x
        abund, merits=getSMerit(model)
        wAbund=pygsl.statistics.wmean(abs(merits[:,0]), abund)
        wAbundSd=pygsl.statistics.wsd(abs(merits[:,0]), abund)
        selModelNo=argmin(abs(merits))
        selModel.update(model[selModelNo])
        print "Selected Model %s with S Abundance %s"%(selModelNo, selModel['comp']['S'])
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    def fitMg(selModel, step, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='mgrun%s/'%step
        dicas, comps=calc_grid_elem('Mg', [0, .3])
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        model=readModels(origspect, run_dir, backup=True)
        #getting the limits to the Ca line by differentiating 2x
        abund, merits=getMgMerit(model)
        wAbund=pygsl.statistics.wmean(abs(merits[:,0]), abund)
        wAbundSd=pygsl.statistics.wsd(abs(merits[:,0]), abund)
        selModelNo=argmin(abs(merits))
        selModel.update(model[selModelNo])
        print "Selected Model %s with Mg Abundance %s"%(selModelNo, selModel['comp']['Mg'])
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    def fitCa(selModel, step, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='carun%s/'%step
        caTester=mean(checkElemLines(selDivspect, 'ca', selSbib['llist']))
        if caTester>1: 
            print "Adding Ca"
            dicas, comps=calc_grid_elem('Ca', [selComp['Ca'], 0.1])
        else: 
            print "Reducing Ca"        
            dicas, comps=calc_grid_elem('Ca', [0, selComp['Ca']])
        multi_dist_launch(dicas,comps,use_machines,run_dir, init_dica=selDica, init_comp=selComp)
        
        model=readModels(origspect, run_dir, backup=True)
        #getting the limits to the Ca line by differentiating 2x
        abund, merits=getCaMerit(model)
        
        wAbund=pygsl.statistics.wmean(abs(merits[:,0]), abund)
        wAbundSd=pygsl.statistics.wsd(abs(merits[:,0]), abund)
        selModelNo=argmin(abs(merits[:,0]))
        
        selModel.update(model[selModelNo])
        print "Selected Model %s with Ca Abundance %s"%(selModelNo, selModel['comp']['Ca'])
        print "Calculated weighted mean and weighted standard deviation for abundance: mean %s sd %s"%(wAbund, wAbundSd)
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    def fitSi(selModel, step, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='sirun1/'
        if siMerit([selModel])[1][0]>0: 
            print "Adding Si"
            if selComp['Si']+selComp['O']>0.6: upperLimit=0.6
            else: upperLimit=selComp['Si']+selComp['O']
            dicas, comps=calc_grid_elem('Si', [selComp['Si'], upperLimit])
        else: 
            print "Reducing Si"        
            dicas, comps=calc_grid_elem('Si', [0, selComp['Si']])
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        #Launching the fit
        model=readModels(origspect, run_dir, backup=True)
        #abund=[item['comp']['Si'] for item in model]
        abund, merits=siMerit(model)
        wAbund=pygsl.statistics.wmean(abs(merits[:, 0]), abund)
        wAbundSd=pygsl.statistics.wsd(abs(merits[:, 0]), abund)
        
        selModelNo=argmin(abs(merits))
        selModel.update(model[selModelNo])
        print "Selected Model %s"%selModelNo
        print "Si abundance %s"%selModel['comp']['Si']
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
        
    #Grid Fitters

    def fitGridLumVph(selModel, step=1, params=None):
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        lumBounds, vphBounds=selModel['params']['lumvph']['bounds']
        run_dir='lumvph_run%s/'%step
        print "Doing grid: Lum: %s vph: %s"%(lumBounds, vphBounds)
        dicas,comps=calc_grid_lum_vph(lumBounds,vphBounds,initSample)
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        model=readModels(origspect, run_dir, backup='True')
        lums,vphs=zip(*[[item['dica']['log_lbol'],item['dica']['v_ph']] for item in model])
        if params['fitMode']=='normGrid':
            merits=getBinMerits(model)
            avgmerit=average(merits, axis=1)
        elif params['fitMode']=='tempGrid':
            #tempBins=[[4380, 4480], [5500, 5600]]
            abund, merits=getLineMerit(model,'Si',selWl=[5739],offsets=[100],sample=50,returnMean=False)
            merits=array([item[0] if item!=[[]] else [10,20] for item in merits])
            #Checking if both lines have the same ratio as in he origspect
            tempTest1=abs(merits[:, 0]-merits[:, 1])
            #Checking if both lines are closeto origspect
            tempTest2=abs(average(merits, axis=1))
            avgmerit=tempTest1+tempTest2
        lumavg, lumstd=pygsl.statistics.wmean(1/avgmerit**2, lums), pygsl.statistics.wsd(1/avgmerit**2,lums)
        vphavg, vphstd=pygsl.statistics.wmean(1/avgmerit**2, vphs), pygsl.statistics.wsd(1/avgmerit**2,vphs)
        lumBounds=[lumavg-lumstd, lumavg+lumstd]
        vphBounds=[vphavg-vphstd, vphavg+vphstd]
        if params['fitMode']=='tempGrid':
            t=model[0]['dica']['t']
            tempMax=calcTemp(t, lumBounds[1], vphBounds[0])
            tempMin=calcTemp(t, lumBounds[0], vphBounds[1])
            print "Temperature min %s and max %s"%(tempMin, tempMax)
        selModelNo=argmin(abs(avgmerit))
        selModel.update(model[selModelNo])
        print "Selected model %s"%selModelNo
        print "Parameters for selected Model: lum %s vph %s"%(selModel['dica']['log_lbol'], selModel['dica']['v_ph'])
        print "New Bounds are: lum %s and vph %s"%(lumBounds, vphBounds)
        print "Calculated lum %s and vph %s. (this is just the weighted mean and will not be used further)"%(lumavg, vphavg)
        
        if params['saveBounds']:
            selModel['params']['lumvph']['bounds']=[lumBounds, vphBounds]
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    
    def fitGridMetals(selModel, step, params=None):
        if params==None:
            fit='log'
        else:
            fit=params['fitMode']
        if selModel==None: selModel=initDicaComp()
        selDica,selComp= selModel['dica'], selModel['comp']
        selDivspect, selSbib=selModel['divspect'], selModel['sbib']
        run_dir='metalrun%s/'%step
        comps=[]
        dicas=[]
        TiCrBounds=selModel['params']['TiCr']['bounds']
        NiBounds=selModel['params']['Ni']['bounds']
        if fit=='log':
            TiCrRange=10**linspace(math.log(TiCrBounds[0], 10),math.log(TiCrBounds[1], 10), 5)
        elif fit=='linear':TiCrRange=linspace(TiCrBounds[0],TiCrBounds[1], 5)
        
        for TiCrAbund in TiCrRange:
            for NiAbund in linspace(0.01, 0.1, 5):
                for FeAbund in linspace(0.001,0.05, 5):
                    comps.append({'Ti':TiCrAbund, 'Cr':TiCrAbund, 'Ni':NiAbund, 'Fe':FeAbund})
                    dicas.append({})
        multi_dist_launch(dicas,comps,use_machines,run_dir,  init_dica=selDica, init_comp=selComp)
        model=readModels(origspect, run_dir, backup=True)
        meritUV=getBinMerits(model,[[0, 3500]],func='mean',doProc=None)
        meritAll=getBinMerits(model,[[0,10000]],func='mean',doProc=None)
        CoAbund, CoMerit=getLineMerit(model,'Co')
        selModelNo=argmin(abs(CoMerit)+abs(meritUV[:, 0])+abs(meritAll[:, 0]))
        print "Selected Model %s"%selModelNo
        print "Abundances for selected Model: TI Cr %s Ni %s Fe %s"%\
            (selModel['comp']['Ti'], selModel['comp']['Ni'], selModel['comp']['Fe'])
        selModel.update(model[selModelNo])
        writeModel2Dalek(selModel, step, run_dir)
        return selModel
    
   
    #MAINLOOP
    fitterDict={'lumvph1':[initLumVphFit, {'saveBounds':True}],
        'ca1':[fitCa, None],
        'si1':[fitSi, None],
        's':[fitS, None], 
        'temp':[fitGridLumVph, {'saveBounds':True,'fitMode':'tempGrid'}],
        'metal1':[fitGridMetals, None],
        'glumvph':[fitGridLumVph, {'saveBounds':True, 'fitMode':'normGrid'}],
        'glumvphns':[fitGridLumVph, {'saveBounds':False,'fitMode':'normGrid'}], 
        'mg':[fitMg, None]}
    for i, step in enumerate(steps):
        if step=='lumvph1':
            fitterDict[step]()
            continue
        else:
            selModel=fitterDict[step][0](selModel, i+initStep, params=fitterDict[step][1])


 
