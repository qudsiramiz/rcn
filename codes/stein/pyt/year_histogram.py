import mysql.connector
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import pickle

def main():
    
    Writefiles=False;
    minDate = '2001-06-01'
    maxDate = '2021-12-31'

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

    if (Writefiles):
       conn = mysql.connector.connect(host="127.0.0.1", port=3309, user="seh", db="Cluster")
       cursor = conn.cursor()

# ==== Si only
    print (".. Si data ...")
    
    outfile='outmatrix/year_histogramresult_Matrx_Si.txt'

    if (Writefiles):
       sqlfile = '../_year_histogram_Matrx_Si.sql'
       fd = open(sqlfile, 'r')
       query = fd.read()
       fd.close()       

       cursor.execute(query)
       result = cursor.fetchall()
       with open(outfile, 'wb') as fp:
          pickle.dump(result, fp)
    else:
       with open (outfile, 'rb') as fp:
          result = pickle.load(fp)


    yearSi = []
    freqSi = []
    for column in result:
        yearSi.append(float(column[0])-.5)
        freqSi.append(float(column[1]))

    numtotSi = int(np.sum(freqSi))

    print ("Num Si", numtotSi)

# === now Fe
    print (".. Fe data ...")
    
    outfile='outmatrix/year_histogramresult_Matrx_Fe.txt'

    if (Writefiles):
       sqlfile = '../_year_histogram_Matrx_Fe.sql'
       fd = open(sqlfile, 'r')
       query = fd.read()

       fd.close()       

       cursor.execute(query)
       result = cursor.fetchall()
       with open(outfile, 'wb') as fp:
          pickle.dump(result, fp)
    else:
       with open (outfile, 'rb') as fp:
          result = pickle.load(fp)


    yearFe= []
    freqFe = []
    for column in result:
        yearFe.append(float(column[0])-.2)
        freqFe.append(float(column[1]))

    numtotFe = int(np.sum(freqFe))

    print ("NumFe", numtotFe)


## === now OMNI ====

    outfile='outmatrix/year_histogram_Matrx_omni.txt'    

    if (Writefiles):
       fd = open('../omni_f107.sql', 'r')
       query = fd.read()
       fd.close()
       cursor.execute(query)
       result = cursor.fetchall()

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
        yearfrac.append(float(column[0])-1.0)
        f107.append(float(column[1]))
        dst.append(float(column[2]))

    print ("Num OMNI data = ",num)

    fe=np.asarray(freqFe[1:18])
    so=np.asarray(f107[1:18])
    ds=np.asarray(dst[1:18])


    plt.rcParams.update({'font.size': 27})


    if (1==1):
#      plt.title("RAPID DE SC "+str(SC))
      fig=plt.figure(figsize=(28,9))      
      ax1 = fig.add_subplot(111)
      ax1.set_ylabel('Number of Si/Fe records',color='black')
      ax1.set_yscale('log')
      ax1.set_xlim([2000,2022])
      ax1.set_ylim([9e2,1.5e4])
      ax1.locator_params(axis='x', nbins=5)
      ax1.bar(yearSi,freqSi,width=0.3, color='blue',edgecolor='blue', hatch="/", alpha=0.6)
      ax1.bar(yearFe,freqFe,width=0.3, color='black',edgecolor='black', hatch="/", alpha=0.6)
      ax1.tick_params(axis='y', colors='black')
      ax1.text(-0.07, 1.02,"a)",fontsize=36, color='blue',transform=ax1.transAxes);

      if (1==1):
         ax2 = ax1.twinx()
         ax2.set_ylabel('F10.7 index',color='red')
         ax2.plot(yearfrac,f107,color='red',linewidth=9,alpha=0.8)
         ax2.set_xlim([2000,2020])
         ax2.set_ylim([61,249])
         ax2.tick_params(axis='y', colors='red')
         plt.text(xtextpos,ytextpos-2*deltaytext,'Total number of Si records = '+str(numtotSi),color='blue',transform=ax1.transAxes)
         plt.text(xtextpos,ytextpos-4*deltaytext,'Total number of Fe records = '+str(numtotFe),color='black',transform=ax1.transAxes)
         panel="F107"
      else:
         ax2 = ax1.twinx()
         ax2.set_ylabel('Dst index [nT]',color='black')
         ax2.plot(yearfrac,dst,'k-',linewidth=4,alpha=0.7)
         ax2.set_xlim([2000,2020])
         ax2.set_ylim([-450,40])
         ax2.tick_params(axis='y', colors='black')
         panel="Dst"
         ax2.grid(True)
         ax2.locator_params(axis='x', nbins=5)
    else:
      cc=plt.figure(figsize=(24,9))
      cc1 = cc.add_subplot(121)
      cc1.set_xlabel('F10.7 index',color='red')
      cc1.set_ylabel('Si/Fe counts',color='black')
      cc1.set_yscale('log')
      cc1.set_xlim([61,199])
      cc1.set_ylim([5e3,2e5])
      cc1.scatter(so,fe,color='b',marker="h",s=600)
 
      cc2 = cc.add_subplot(122)
      cc2.set_xlabel('Dst index [nT]',color='red')
      cc2.set_ylabel('Si/Fe counts',color='black')
      cc2.set_yscale('log')
      cc2.set_xlim([-570,10])
      cc2.set_ylim([5e3,1.3e5])
      cc2.scatter(ds,fe,color='b',marker="h",s=600)

#    plt.savefig("../../Pap/figs/hist_year_SiFe_"+panel+".png",bbox_inches='tight')
    plt.show()
 



    print (yearFe[1])
    print (yearfrac[1])

    N=np.size(fe)
    print ("CC Fe vs F10.7 = ",np.corrcoef(fe,so))
    print ("CC Fe vs Dst = ",np.corrcoef(fe,ds))




main()
