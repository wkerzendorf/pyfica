import cherrypy
import sqlite3
from pyfica import dalekDB
cherrypy.config.update({'server.socket_port': 8088,}) 

class GAExplorer(object):
    def index(self):
        conn = dalekDB.getDBConnection('tmp.db3')
        curs = conn.cursor()
        dbString = "<table>"
        dbString += "<tr><td>GA RUN ID</td><td>Description</td><td>start time</td></tr>"
        for item in conn.execute('select id, description,start_time from GA_RUN'):
            dbString += "<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % item
        dbString += "</table>"
        conn.close()
        return dbString
    index.exposed = True

cherrypy.quickstart(GAExplorer())
 