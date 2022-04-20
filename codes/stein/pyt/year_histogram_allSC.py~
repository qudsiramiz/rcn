import mysql.connector
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import pickle

def main():
    
    minDate = '2001-06-01'
    maxDate = '2019-12-31'

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

    w=0.15		# width of bars

# --- writing text files for later use ---

    conn = mysql.connector.connect(host="127.0.0.1", port=3309, user="seh", db="Cluster")    
    cursor = conn.cursor()

    writeFiles=True;
    if writeFiles :	
    	fd = open('../_year_histogram_py.sql', 'r')
    	query = fd.read()
    	cursor.execute(query)
    	result = cursor.fetchall()
    	outfile='outmatrix/year_histogramresult_all.txt'
    	fd.close()

    	with open(outfile,'w') as fp:
		pickle.dump(result,fp)
    	fp.close()


   	fd = open('../_year_histogram_py.sql', 'r')
    	query = fd.read()
    	query=query.replace("@SCID","2")
    	cursor.execute(query)
    	result = cursor.fetchall()
    	outfile='outmatrix/year_histogramresult_C2.txt'
    	fd.close()

    	with open(outfile,'w') as fp:
		pickle.dump(result,fp)
    	fp.close()

  	fd = open('../_year_histogram_py.sql', 'r')
    	query = fd.read()
    	query=query.replace("@SCID","3")
    	cursor.execute(query)
    	result = cursor.fetchall()
    	outfile='outmatrix/year_histogramresult_C3.txt'
    	fd.close()

    	with open(outfile,'w') as fp:
		pickle.dump(result,fp)
    	fp.close()

  	fd = open('../_year_histogram_py.sql', 'r')
    	query = fd.read()
    	query=query.replace("@SCID","4")
    	cursor.execute(query)
    	result = cursor.fetchall()
    	outfile='outmatrix/year_histogramresult_C4.txt'
    	fd.close()

    	with open(outfile,'w') as fp:
		pickle.dump(result,fp)
    	fp.close()




# == Nwo reading C1 only ===
    
    outfile='outmatrix/year_histogramresult_C1.txt'
    with open (outfile, 'rb') as fp:
          result = pickle.load(fp) 

    year1 = []
    freq1 = []
    for column in result:
        year1.append(float(column[0])-.5)
        freq1.append(float(column[1]))

    numtot1 = int(np.sum(freq1))
    print numtot1


# ==== C2 ===
    
    outfile='outmatrix/year_histogramresult_C2.txt'
    with open (outfile, 'rb') as fp:
          result = pickle.load(fp) 

    year2 = []
    freq2 = []
    for column in result:
        year2.append(float(column[0])-.5+w)
        freq2.append(float(column[1]))

    numtot2 = int(np.sum(freq2))
    print numtot2

# ==== C3 ===
    
    outfile='outmatrix/year_histogramresult_C3.txt'
    with open (outfile, 'rb') as fp:
          result = pickle.load(fp) 

    year3 = []
    freq3 = []
    for column in result:
        year3.append(float(column[0])-0.5+2*w)
        freq3.append(float(column[1]))

    numtot3 = int(np.sum(freq3))
    print numtot3

# ==== C4 ===
    
    outfile='outmatrix/year_histogramresult_C4.txt'
    with open (outfile, 'rb') as fp:
          result = pickle.load(fp) 

    year4 = []
    freq4 = []
    for column in result:
        year4.append(float(column[0])-.5+3*w)
        freq4.append(float(column[1]))

    numtot4 = int(np.sum(freq4))
    print numtot4

    numtotAll=numtot1+numtot2+numtot3+numtot4


# ===== F10.7 ======

    outfile='outmatrix/year_histogram_omni.txt'    
    if writeFiles:
       fd = open('../omni_f107.sql', 'r')
       query = fd.read()
       cursor.execute(query)
       result = cursor.fetchall()
       fd.close()
       with open(outfile, 'wb') as fp:
          pickle.dump(result, fp)
    else:
       with open (outfile, 'rb') as fp:
          result = pickle.load(fp) 

    num=0
    yearfrac = []
    f107 = []
    dst=[]
    for column in result:
        num=num+1
        yearfrac.append(float(column[0]))
        f107.append(float(column[1]))
        dst.append(float(column[2]))

    print "Num OMNI data = ",num

#    fe=np.asarray(freqFe[1:18])
#    so=np.asarray(f107[1:18])
#    ds=np.asarray(dst[1:18])


    plt.rcParams.update({'font.size': 25})

    a=0.6

    if (1==1):
      fig=plt.figure(figsize=(28,9))
      ax1 = fig.add_subplot(111)
      ax1.set_ylabel('Number of RMATRX records',color='black')
#      ax1.set_yscale('log')
      ax1.set_xlim([2000,2020])
      ax1.set_ylim([1.1e4,4.4e4])
      ax1.bar(year1,freq1,width=w,color='black',edgecolor='black',alpha=a)
      ax1.bar(year2,freq2,width=w,color='red',edgecolor='red',alpha=a)
      ax1.bar(year3,freq3,width=w,color='green',edgecolor='green',alpha=a)
      ax1.bar(year4,freq4,width=w,color='blue',edgecolor='blue',alpha=a)
      ax1.tick_params(axis='y', colors='black')
      ax1.legend(['SC1','SC2','SC3','SC4'],loc=1)

      ax2 = ax1.twinx()
      ax2.set_ylabel('F10.7 index',color='black')
      ax2.plot(yearfrac,f107,'-',color='grey',alpha=0.8,linewidth=6)
      ax2.plot(yearfrac,f107,'k-',linewidth=2)
      ax2.set_xlim([2000,2020])
      ax2.set_ylim([61,249])
      ax2.tick_params(axis='y', colors='black')
    
      plt.text(xtextpos,ytextpos-3*deltaytext,'Total number of RMATRX records = '+str(numtotAll),transform=ax1.transAxes)
    else:
      cc=plt.figure(figsize=(24,9))
      cc1 = cc.add_subplot(121)
      cc1.set_xlabel('F10.7 index',color='blue')
      cc1.set_ylabel('Fe counts',color='black')
      cc1.set_yscale('log')
      cc1.set_xlim([61,199])
      cc1.set_ylim([7e3,2e5])
      cc1.scatter(so,fe,color='b',marker="h",s=600)
 
      cc2 = cc.add_subplot(122)
      cc2.set_xlabel('Dst index [nT]',color='blue')
      cc2.set_ylabel('Fe counts',color='black')
      cc2.set_yscale('log')
      cc2.set_xlim([-570,10])
      cc2.set_ylim([5e3,1.3e5])
      cc2.scatter(ds,fe,color='b',marker="h",s=600)

    plt.savefig("../../Pap/figs/hist_year_allSC.png",bbox_inches='tight')
    plt.show()
 



#    print yearFe[1]
#    print yearfrac[1]

#    N=np.size(fe)
#    print "CC Fe vs F10.7 = ",np.corrcoef(fe,so)
#    print "CC Fe vs Dst = ",np.corrcoef(fe,ds)


main()
