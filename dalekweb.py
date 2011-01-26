import cherrypy
import sqlite3
from pyfica import dalekDB, model
cherrypy.config.update({'server.socket_port': 8088,}) 
htmlHeader = "<html><body>"
htmlFooter = "</body></html>"
dbName = '2002bo_test.db3'
import StringIO
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import pprint
import pdb
pp = pprint.PrettyPrinter(indent=4)

SNDBDict = {'2002bo':'2002bo_test.db3'}
tableColumnElement = "<td>%s</td>"

collapseScript = """
<script language="JavaScript" type="text/javascript">
<!--
function sizeTbl(h, elementID) {
  var tbl = document.getElementById(elementID);
  tbl.style.display = h;
}
// -->
</script> """

collapseSection = """
%s
<div id=%s name=%s style="overflow:hidden;display:none">
%s
</div>



<a href="javascript:sizeTbl('none', '%s')">Hide</a>
 
<a href="javascript:sizeTbl('block', '%s')">Expand</a>
<br>
"""


def sn_site(kwargs):
    
    dbPath = SNDBDict[kwargs['name']]
    snName = kwargs['name']
    #prepare string for sn site
    conn = dalekDB.getDBConnection(dbPath)
    curs = conn.cursor()
    
    siteString = htmlHeader + "<table>"
    siteString += "<tr>" + tableColumnElement * 4 + "</tr>\n"
    siteString = siteString % ("GA_RUN_ID", "START DATE", "current gen", "SN DATE")
    dbSelectStmt = """select GA_RUN.ID, GA_RUN.START_TIME,
                    GA_RUN.GA_CUR_GEN, SN_SPECTRA.DATE from 
                    GA_RUN, SN_SPECTRA where GA_RUN.SN_ID=SN_SPECTRA.ID
                    """
    for item in curs.execute(dbSelectStmt):
       siteString += "<tr><td><a href=/gaexplore?name=%s&garunid=%d>%s</a></td>" \
                     % (snName,item[0],item[0])
       siteString += tableColumnElement * 3 + "</tr>"
       siteString = siteString % item[1:]
    siteString += "</table>" + htmlFooter
    conn.close()
    return siteString

def garun_site(kwargs):
    dbPath = SNDBDict[kwargs['name']]
    snName = kwargs['name']
    gaRunID = int(kwargs['garunid'])
    #prepare string for sn site
    conn = dalekDB.getDBConnection(dbPath)
    curs = conn.cursor()
    
    siteString = htmlHeader + collapseScript
    siteString += "<img src=garun_plot?name=%s&garunid=%d><p>" % (snName, gaRunID)
    
    gaConfDict = pp.pformat(dalekDB.convertZipPickle(
                    curs.execute('select GA_CONF_DICT from GA_RUN where ID=%d'
                                 % gaRunID).fetchone()[0]))
    gaConfDictTuple = ('GA_CONF_DICT code', 'gaconfdict', 'gaconfdict', gaConfDict, 'gaconfdict', 'gaconfdict')
    siteString += collapseSection % gaConfDictTuple
    
    gaBreedFunc = dalekDB.convertZipPickle(
                    curs.execute('select GA_BREED_FUNC from GA_RUN where ID=%d'
                                 % gaRunID).fetchone()[0])
    gaBreedTuple = ('breed code', 'gabreed', 'gabreed', gaBreedFunc, 'gabreed', 'gabreed')
    siteString += collapseSection % gaBreedTuple
    
    gaCrossFunc = dalekDB.convertZipPickle(
                    curs.execute('select GA_CROSS_FUNC from GA_RUN where ID=%d'
                                 % gaRunID).fetchone()[0])
    gaCrossTuple = ('cross code', 'gacross', 'gacross', gaCrossFunc, 'gacross', 'gacross')
    siteString += collapseSection % gaCrossTuple
    
    gaSelectFunc = dalekDB.convertZipPickle(
                    curs.execute('select GA_SELECT_FUNC from GA_RUN where ID=%d'
                                 % gaRunID).fetchone()[0])
    gaSelectTuple = ('select code', 'gaselect', 'gaselect', gaSelectFunc, 'gaselect', 'gaselect')
    siteString += collapseSection % gaSelectTuple
    
    
    
    
    
    siteString += "<table>\n"
    tmpString = "<tr>" + tableColumnElement * 2 + "</tr>\n"
    siteString += tmpString % ("ID", "FITNESS")
    dbSelectStmt = """select GA_GENERATION.ID, max(GA_INDIVIDUAL.FITNESS) from GA_GENERATION
                    join GA_INDIVIDUAL on GA_GENERATION.ID=GA_INDIVIDUAL.GENERATION_ID where
                    GA_GENERATION.GA_RUN_ID=%d
                    group by GA_INDIVIDUAL.GENERATION_ID""" % gaRunID
                    
    for item in curs.execute(dbSelectStmt):
       siteString += "<tr><td><a href=/gaexplore?name=%s&garunid=%d&genid=%d>%s</a></td>" \
                     % (snName, gaRunID, item[0], item[0])
       tmpString = tableColumnElement + "</tr>\n"
       siteString += tmpString % item[1]
    siteString += "</table>" + htmlFooter
    conn.close()
    return siteString


def generation_site(kwargs):
    dbPath = SNDBDict[kwargs['name']]
    snName = kwargs['name']
    gaRunID = int(kwargs['garunid'])
    genID = int(kwargs['genid'])
    #prepare string for sn site
    conn = dalekDB.getDBConnection(dbPath)
    curs = conn.cursor()
    
    siteString = htmlHeader + "<table>\n"
    siteString += "<tr>" + tableColumnElement * 16 + "</tr>\n"
    siteString = siteString % ("INDIVIDUAL ID", "LUM",
                               "VPH", "C", "O", "NA", "MG", "SI",
                               "S", "CA", "TI", "CR", "MN", "FE0",
                               "NI0", "FITNESS")
    dbSelectStmt = """select GA_INDIVIDUAL.ID, FICA_MODEL.ID, FICA_LUMVPH.LUM, FICA_LUMVPH.VPH, FICA_ABUNDANCE.C,
                    FICA_ABUNDANCE.O, FICA_ABUNDANCE.NA, FICA_ABUNDANCE.MG, FICA_ABUNDANCE.SI,
                    FICA_ABUNDANCE.S, FICA_ABUNDANCE.CA, FICA_ABUNDANCE.TI, FICA_ABUNDANCE.CR,
                    FICA_ABUNDANCE.MN, FICA_ABUNDANCE.FE0, FICA_ABUNDANCE.NI0, GA_INDIVIDUAL.FITNESS
                    from FICA_DICA, FICA_ABUNDANCE, FICA_LUMVPH, FICA_MODEL, GA_INDIVIDUAL
                    where GA_INDIVIDUAL.GENERATION_ID=%d and GA_INDIVIDUAL.MODEL_ID=FICA_MODEL.ID
                    and FICA_ABUNDANCE.ID=FICA_MODEL.ABUNDANCE_ID
                    and FICA_DICA.ID=FICA_MODEL.DICA_ID
                    and FICA_LUMVPH.ID=FICA_MODEL.LUMVPH_ID
                    ORDER BY GA_INDIVIDUAL.FITNESS DESC
                    """ % (genID)
                    
    for item in curs.execute(dbSelectStmt):
       #siteString += "<tr><td><a href=/gaexplore?name=%s&garunid=%d&genid=%d&modelid=%d>%s</a></td>" \
        #             % (snName, gaRunID, genID, item[0], item[0])
       siteString += "<tr><td><a href=/model_plot?name=%s&modelid=%d&garunid=%d>%s</a></td>" \
                     % (snName, item[1], gaRunID, item[0], )
       siteString += tableColumnElement*15 + "</tr>\n"
       siteString = siteString % item[2:]
    siteString += "</table>" + htmlFooter
    conn.close()
    return siteString

def model_site(kwargs):
    dbPath = SNDBDict[kwargs['name']]
    snName = kwargs['name']
    gaRunID = int(kwargs['garunid'])
    genID = int(kwargs['genid'])
    #prepare string for sn site
    conn = dalekDB.getDBConnection(dbPath)
    curs = conn.cursor()
    
    


class GAExplorer(object):
    
    @cherrypy.expose
    def index(self):
        siteString = htmlHeader+"<table>"
        siteString += "<tr><SN Name></tr>"
        for item in SNDBDict.items():
            siteString += "<tr><td><a href=\"gaexplore?name=%s\">%s</a></td></tr>" % (item[0],item[0])
        siteString += "</table>" + htmlFooter
        return siteString
    
    @cherrypy.expose
    def gaexplore(self, **kwargs):
        snPath = SNDBDict[kwargs['name']]
        if 'genid' in kwargs.keys():
            return generation_site(kwargs)
        if 'garunid' in kwargs.keys():
            return garun_site(kwargs)
        if 'name' in kwargs.keys():
            return sn_site(kwargs)
        return str(dict(kwargs))
    
 
            
    
    @cherrypy.expose
    def garun_plot(self, **kwargs):
        dbPath = SNDBDict[kwargs['name']]
        cherrypy.response.headers['Content-Type'] = "image/png"
        snName = kwargs['name']
        #prepare string for sn site
        conn = dalekDB.getDBConnection(dbPath)
        curs = conn.cursor()
        gaRunID = int(kwargs['garunid'])
        dbSelectStmt = """select GA_GENERATION.ID,avg(GA_INDIVIDUAL.FITNESS), max(GA_INDIVIDUAL.FITNESS) from GA_GENERATION
                    join GA_INDIVIDUAL on GA_GENERATION.ID=GA_INDIVIDUAL.GENERATION_ID where
                    GA_GENERATION.GA_RUN_ID=%d
                    group by GA_INDIVIDUAL.GENERATION_ID""" % gaRunID
    
                      
        generation, fitAvg, fitMax= zip(*curs.execute(dbSelectStmt).fetchall())
        fh = StringIO.StringIO()
        fig = Figure()
        canvas = FigureCanvas(fig)
        ax = canvas.figure.add_subplot(111)
        ax.plot(generation, fitAvg, 'b-', label='Fitness average')
        ax.plot(generation, fitMax, 'r-', label='Fitness max')
        ax.legend()
        canvas.print_figure(fh, fmt='png')
        fh.seek(0)
        conn.close()
        return fh.read()
        

        
  
        
        
    @cherrypy.expose    
    def model_plot(self, **kwargs):
        cherrypy.response.headers['Content-Type'] = "image/png"
        snName = kwargs['name']
        modelID = int(kwargs['modelid'])
        gaRunID = int(kwargs['garunid'])
        dbPath = SNDBDict[kwargs['name']]
        
        conn = dalekDB.getDBConnection(dbName)
        curs = conn.cursor()
        fh = StringIO.StringIO()
        modelind = model.model.fromDB(conn, modelID, gaRunID)
        fig = Figure()
        canvas = FigureCanvas(fig)
        ax = canvas.figure.add_subplot(111)
        ax.plot(modelind.origSpec.x, modelind.origSpec.y, 'k-', label='observed')
        ax.plot(modelind.aSpec.x, modelind.aSpec.y, label='model')
        ax.legend()
        canvas.print_figure(fh, fmt='png')
        fh.seek(0)
        conn.close()
        return fh.read()

    
cherrypy.quickstart(GAExplorer())
 