
#
#  python3 flux_vs_record.py 
#

import mysql.connector
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as ml
import math
import sys
import random
import pickle
from numpy import loadtxt
from datetime import date
import matplotlib.dates as mdates

def main():

    Species = "Fe"    # H, He, O, Si, Fe
    Moment='avg'
    minDate = '\"2003-10-29\"'
    maxDate = '\"2003-10-30\"'

# get from arguments instead
    if len(sys.argv)>=2:
    	minDate=sys.argv[1]
    	maxDate=sys.argv[2]
   
# plot settings
    msize=6		# markers in mass lines 
    fsize=22		# font size

    xtextpos = 0.50 	# x-position of text /annotation
    ytextpos = 0.99	# y-position of top text /annotation
    deltaytext=0.05     # y-spacing between text lines


# -- Mysql connection ----
    conn = mysql.connector.connect(host="127.0.0.1", port=3309, user="seh", db="Cluster")    
    cursor = conn.cursor()
 

    sqlfile = '../_fluxes_vs_record.sql'
    fd = open(sqlfile, 'r') 
    query = fd.read()
    fd.close()


    query=query.replace("@mom",Moment)
    query=query.replace("@i",Species)

    query=query.replace("@minX","-25").replace("@maxX","25")
    query=query.replace("@minY","-15").replace("@maxY","15")
    query=query.replace("@minZ","-25").replace("@maxZ","25")
    query=query.replace("@minDate",minDate).replace("@maxDate",maxDate)
    
    print (query)
    cursor.execute(query)
    result = cursor.fetchall()

    num=0;
    intensity=[]
    time=[]
    for col in result:
       print(col[0], col[1])
       time.append(col[0])
       intensity.append(col[1])

 #  --- Energy channels - set log10 middle -----
    Energy_range = '709 - 1080 keV'
    if Species == "Si":
       Energy_range='410-648 keV';

    ESi=[410,648,818,934,1082,1272,1520,1854,2066]
    E = ESi

    for i in range(0,8):
       logE= math.log10(E[i])+math.log10(E[i+1])
       E[i] = 10**(logE/2)

    EnergyCh = E[:8]
   
 
    plt.figure(figsize=(30,8))
    plt.rcParams.update({'font.size': fsize})

 
    Si=plt.subplot(1,1,1)
    Si.set_xlabel("UT time");
    Si.set_xticks(time);
    Si.set_ylabel("Intensity [1/(keV s)]");
#    Si.grid(True,which='major')
#    Si.fill_between(range(len(day)), min(day), max(day), where=(day < 10) & (day > 5), alpha=0.5)
    Si.set_yscale('log');
    Si.set_title('Intensity for '+Species+'$^{N+}$,  Cluster Xgse < 0,   Energy range = '+Energy_range + 'Date = '+minDate+ ' - ' +maxDate)
    Si.plot(time,intensity, marker='o',markersize=10,linewidth=1)
#    Si.tick_params(axis ='x', rotation = 30)

    myFmt = mdates.DateFormatter('%H:%M')
    Si.xaxis.set_major_formatter(myFmt)


    plt.savefig("../../fluxplots/"+minDate.strip("\"")+"___"+maxDate.strip("\"")+"_Yamaplot_intensity_record_by_record"+Species+".png",bbox_inches='tight')   
    plt.show()


main()

