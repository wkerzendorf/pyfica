day2sec=3600*24
sigma=5.6705e-8
L_sun=3.846e26

from math import *
import numpy
from numpy import sort, array,loadtxt,float32
from pyphot import datafile
class sbibfile(file):
    def read_data(self):
        data_dict={}
        filedata=self.readlines()
        tmp_row=[float(item) for item in filedata[3].split()]
        data_dict['log_r'], data_dict['log_lbol'],data_dict['v_ph'],data_dict['TB']=tmp_row
        #Reading the element information
        data_dict['elems']=[]
        row_format=[int, str, float, float, float]
        for i, row in enumerate(filedata[5:]):
            if row=='\n': break
            row_elems=[func(item) for func, item in zip(row_format, row.split())]
            data_dict['elems'].append(row_elems)
        #Reading the line list
        row_format={'names':('eqw','shift','rest','atom','ion','param1','param2','param3'),
            'formats':(float32,float32,float32,'|S4','|S8',float32,float32,float32)}
        data_dict['llist']=[]
        tmpdata=[]
        for j, row in enumerate(filedata[8+i:]):
            if row=='\n': break
            tmpdata.append(tuple(row.split()))
        data_dict['llistTuple']=tmpdata
        data_dict['llist']=array(tmpdata,row_format)
        return data_dict
                
class dicafile(file):
    def read_data(self):
        #print "reading data"
        data_dict={}
        filedata=self.readlines()
        data_dict['np5'],data_dict['itt']=map(float,filedata[1].split('|')[0].split())
        line2=map(float,filedata[3].split('|')[0].split())
        if len(line2)==3:
            data_dict['js'],data_dict['mb'],data_dict['kb']=line2
            data_dict['xe1']=0.2
        elif len(line2)==4:
            data_dict['js'],data_dict['mb'],data_dict['kb'],data_dict['xe1']=line2
        else: raise Exception('Problem with line2 input file')
        data_dict['log_lbol'],data_dict['v_ph'],data_dict['tb']=map(float,filedata[5].split('|')[0].split())
        data_dict['t'],data_dict['m-m'],data_dict['z']=map(float,filedata[7].split('|')[0].split())
        data_dict['e_b-v_gal'],data_dict['e_b-v_host']=map(float,filedata[9].split('|')[0].split())
        data_dict['kr'],data_dict['lg_tau']=map(float,filedata[11].split('|')[0].split())
        data_dict['chl'],data_dict['nc']=map(float,filedata[13].split('|')[0].split())
        data_dict['wl'],data_dict['grid'],data_dict['mu']=map(float,filedata[15].split('|')[0].split())
        data_dict['intlow'],data_dict['inthigh'],data_dict['log_l_low_high']=map(float,filedata[17].split('|')[0].split())
        data_dict['em_low'],data_dict['em_high']=map(float,filedata[19].split('|')[0].split())
        data_dict['options']=filedata[21].split('|')[0].strip()
        self.data=data_dict
        return data_dict
    
    def write_data(self,data=""):
        if data=="" and not hasattr(self, 'data'): raise Exception('Please specify data or run read_data. to write data')
        elif data=="": data=self.data
        data_dict={"np5":-1, "itt":-1,
                      "js":-1, "mb": -1, "kb":-1,'xe1':0.2,
                      "log_lbol":-1,"v_ph":-1,'tb':-1,
                      "t":-1,"m-m":-1,"z":-1,
                      "e_b-v_gal":-1,"e_b-v_host":-1,
                      "kr":-1,"lg_tau":-1,
                      "chl":-1,"nc":-1,
                      "wl":-1,"grid":-1,"mu":-1, 
                      "intlow":-1,"inthigh":-1,"log_l_low_high":-1,
                      "em_low":-1,"em_high":-1,
                      "options":"1   1   1   1   1"}
        for key,value in data.items():
            if data_dict.has_key(key): data_dict[key]=value
            else: raise Exception("The key you have specified is not a parameter %s"%key)
        #print "writing data"
        self.write("\n")
        self.write("     %d    %d                            | NP5 	 ITT\n"%(data_dict['np5'],data_dict['itt']))
        self.write("\n")
        self.write("   %02d   %05d   %03d  %4.2f               | JS    MB    KB\n"%(data_dict['js'],data_dict['mb'],data_dict['kb'],data_dict['xe1']))
        self.write("\n")
        self.write("   %05.03f   %06.01f   %07.01f            | Log L(Bol)/Ls  Vph(km/s)  Tb(K)\n"%(data_dict['log_lbol'],data_dict['v_ph'],data_dict['tb']))
        self.write("\n")
        self.write("    %05.02f   %05.02f    %06.04f              | t (dys)   m-M    redshift\n"%(data_dict['t'],data_dict['m-m'],data_dict['z']))
        self.write("\n")
        self.write("    %0.02f    %0.02f                        | E(B-V)Gal 	E(B-V)Host\n"%(data_dict['e_b-v_gal'],data_dict['e_b-v_host']))
        self.write("\n")
        self.write("   %s   % 04.02f                    | KR    Log tau(limit)\n"%(('%09d'%data_dict['kr'])[:9],data_dict['lg_tau']))
        self.write("\n")
        self.write("  %04.01f       %d                          | CHL(eV)   NC (no. of cascades)\n"%(data_dict['chl'],data_dict['nc']))
        self.write("\n")
        self.write("   %06.04f   %06.04f   %04d               | Wavelength grid - mu \n"%(data_dict['wl'],data_dict['grid'],data_dict['mu']))
        self.write("\n")
        self.write("   %07.02f   %07.02f   %05.03f            | Obs interval in A   Log L(wv1,wv2)\n"%(data_dict['intlow'],data_dict['inthigh'],data_dict['log_l_low_high']))
        self.write("\n")
        self.write("   %07.02f   %07.02f                    | Em interval at lwr bdry  - KC(3)=0\n"%(data_dict['em_low'],data_dict['em_high']))
        self.write("\n")
        self.write("   %s                    | Options  -  see FICA.NOTES\n"%data_dict['options'])

class compfile(file):
    def write_data(self,data):
        nucl_dict={'Al': 13,
                   'Ar': 18,
                   'B': 5,
                   'Be': 4,
                   'C': 6,
                   'Ca': 20,
                   'Cl': 17,
                   'Co': 27,
                   'Cr': 24,
                   'Cu': 29,
                   'F': 9,
                   'Fe': 26,
                   'H': 1,
                   'He': 2,
                   'K': 19,
                   'Li': 3,
                   'Mg': 12,
                   'Mn': 25,
                   'N': 7,
                   'Na': 11,
                   'Ne': 10,
                   'Ni': 28,
                   'O': 8,
                   'P': 15,
                   'S': 16,
                   'Sc': 21,
                   'Si': 14,
                   'Ti': 22,
                   'V': 23,
                   'Zn': 30,
                   'Ni0':31,
                   'Fe0':32}
        data_dict={'Al': -1.0,
                   'Ar': -1.0,
                   'B': -1.0,
                   'Be': -1.0,
                   'C': -1.0,
                   'Ca': -1.0, 
                   'Cl': -1.0,
                   'Co': -1.0,
                   'Cr': -1.0,
                   'Cu': -1.0,
                   'F': -1.0,
                   'Fe': -1.0,
                   'H': -1.0,
                   'He': -1.0,
                   'K': -1.0,
                   'Li': -1.0,
                   'Mg': -1.0,
                   'Mn': -1.0,
                   'N': -1.0,
                   'Na': -1.0,
                   'Ne': -1.0,
                   'Ni': -1.0,
                   'O': -1.0,
                   'P': -1.0,
                   'S': -1.0,
                   'Sc': -1.0,
                   'Si': -1.0,
                   'Ti': -1.0,
                   'V': -1.0,
                   'Zn': -1.0,
                   'Ni0':-1.0,
                   'Fe0':-1.0}
        for key,value in data.items():
            if data_dict.has_key(key): data_dict[key]=value
            else: raise Exception("The key you have specified (%s) is not a parameter."%key)
        #building lookup dict
        lookup_dict=dict([(b,a) for a,b in nucl_dict.items()])
        for key,value in data_dict.items():
            lookup_dict.update({nucl_dict[key]:key})
        no_list=sort(lookup_dict.keys())
        for ino in no_list:
            elem=lookup_dict[ino]
            self.write("%2d %-2s    %08.06f\n"%(ino,elem,data_dict[elem]))

    def read_data(self):
        data_dict={}
        for line in self:
            ilist=line.split()
            data_dict.update({ilist[1]:float(ilist[2])})
        return data_dict
        
        
class spctfile(datafile):
    def read_data(self):
        return array(datafile.read_data(self, sel_columns=[0,2]))
class ststfile(file):
    def _readWParam(self):
        print

def get_temp(data):
        L_sn=(10**(data['log_lbol']))*L_sun
        t=data['t']*day2sec
        r=data['v_ph']*1e3*t
        T=(L_sn/(4*pi*sigma*(r**2)))**0.25
        return T
def set_temp(data, T, keep='lum'):
    #keep sets what to keep, the choices are 'lum' or 'v_ph' for photospheric velocity
    print "setting temperature"
    T=float(T)
    if keep=='lum':
        L_sn=(10**(data['log_lbol']))*L_sun
        t=data['t']*day2sec
        v_ph=((L_sn/(4*pi*sigma*(T**4)*(t**2) ))**0.5)*1e-3
        print "Photospheric velocity is %s km/s"%v_ph
        data['v_ph']=v_ph
    elif keep=='v_ph':        
        t=data['t']*day2sec
        r=data['v_ph']*1e3*t
        L_sn=4*pi*sigma*(r**2)*(T**4)
        log_lbol=log(L_sn/L_sun, 10)
        print "Log(L/L_sol) is %s"%log_lbol
        data['log_lbol']=log_lbol
        
    else:
        raise Exception("the second argument 'keep' can only have two options: 'lum' and 'v_ph' ")
    
def set_co(t, comp):
    ni_abundance=comp['ni']
    co_abundance=(ni_abundance*(1-pow(2, -(t/6.1))))
    comp['ni']=(ni_abundance*(pow(2, -(t/6.1))) )
    print " The calculated Co abundance is %s "% co_abundance
    comp['co']+=co_abundance
def set_ni_chain(t, comp):
    t1=6.1
    t2=77.27
    ni_abundance_0=comp['ni']
    ni_abundance_t=ni_abundance_0*(pow(2, -(t/t1)))
    co_abundance_t=(t2/(t1-t2))*ni_abundance_0*(pow(2, -t/t1)-pow(2, -t/t2))
    fe_abundance_t=ni_abundance_0-co_abundance_t-ni_abundance_t
    comp['ni']=ni_abundance_t
    comp['co']=co_abundance_t
    comp['fe']+=fe_abundance_t
    print "Ni %s Co %s Fe %s"%(ni_abundance_t, co_abundance_t, fe_abundance_t)
    
def normalize_comp(comp):
    norm_factor=sum(comp.values())
    for key in comp.keys():
        comp[key]/=norm_factor
