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


def main():

    Species = "Si"    # H, He, O, Si, Fe
    Moment  = "avg"   # min, max, avg, median, stddev..... 
    minDate = '\"2001-01-01\"'
    maxDate = '\"2007-12-31\"'

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
 

    sqlfile = '../matrx_spectrum.sql'
    fd = open(sqlfile, 'r') 
    query = fd.read()
    fd.close()


    query=query.replace("mom",Moment)
    query=query.replace("@i",Species)

    query=query.replace("@minX","-25").replace("@maxX","0")
    query=query.replace("@minY","-25").replace("@maxY","25")
    query=query.replace("@minZ","-25").replace("@maxZ","25")
    query=query.replace("@minDate",minDate).replace("@maxDate",maxDate)
    
    print query
    cursor.execute(query)
    result = cursor.fetchall()

    num=0;
    c0=[]
    for col in result:
        num=num+1
	for i in range(0,9):
	   c0.append(col[i])

    num = c0[8];
    print num;
    c0=c0[:8]

#  --- Energy channels - set log10 middle -----

    ESi=[410,648,818,934,1082,1272,1520,1854,2066]
    E = ESi

    for i in range(0,8):
       logE= math.log10(E[i])+math.log10(E[i+1])
       E[i] = 10**(logE/2)

    EnergyCh = E[:8]




  


    plt.figure(figsize=(15,10))
    plt.rcParams.update({'font.size': fsize})

 
    Si=plt.subplot(1,1,1)
    Si.set_xlabel("Energy [keV]");
    Si.set_ylabel("Counts/s");
    Si.grid(True,which='both')
    Si.set_yscale('log');
    Si.set_xscale('log');
    Si.set_xlim(1e2,4e3)
    Si.set_title('Si, all nightside, 2001-2002')
#    Si.set_ylim(1e-3,1e0)
    Si.plot(EnergyCh,c0, marker='o',markersize=20,linewidth=2)
 
    plt.savefig("../../PapSi/figs/matrx_spectrum_"+Species+".png",bbox_inches='tight')   
    plt.show()


main()

