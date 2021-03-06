import mysql.connector
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as ml
import math
import sys
import random
import pickle
from numpy import loadtxt
from datetime import date


def main():

    param='f107'
    split=''
    writeFiles=False
    minCount=0.15   #(0.15 ~ 10/64  ~ 10 counts)
    maxCount=250
    minSpecies=1

    minDate = '2001-06-01'
    maxDate = '2019-12-31'

    minX=-22
    maxX=22

    minY=-22
    maxY=22

    minZ=-22
    maxZ=22

    minR=0
    maxR=22

    minLat=00
    maxLat=99

    minDst=-999;   
    maxDst=-50;

    minF107=0;
    maxF107=500;

    region='all'
    label=['b)', 'c)']
    csys='GSM'
    ptit="Disturbed conditions; Dst < "+str(maxDst)

# get from arguments instead
    if len(sys.argv)>1:
    	param=sys.argv[1]
        split=sys.argv[2]

    if region=='cusp':
	# Region_cusp_x-20_20_y-20_20_z-20_20_r0_20_dst-2500_300.png
    	minLat=60
    	maxLat=99
        label=['d)', 'e)', 'f)']
	ptit="Cusp and polar cap: "
    elif region=='sw':
        minX=0
	minR=13
        csys='GSE'
        label=['a)', 'b)', 'c)']
        ptit="Solar wind : "
    elif region=='ps':
        maxX=-5
	minY=-15
	maxY=15
        minZ=-4
	maxZ=4
        label=['j)', 'k)', 'l)']
        ptit="Plasma sheet: "
    elif region=='lobe':
	maxX=-5
	minY=-15
	maxY=15
	minZ=4
	label=['g)', 'h)', 'i)']
        ptit="Magnetotail lobes: "
    elif region=='inner':
	maxR=6.6    
        label=['m)', 'n)', 'o)']
        ptit="Inner magnetosphere: "



# -- High activity ----

    outfileHi='outmatrix/disturbance_spectrum_dst_disturbed'+split+'.txt'
    maxDst=-35;
    minDst=-999;
    minF107=160;
    maxF107=999;
 
    if param=='f107':
	outfileHi=outfileHi.replace("_dst_","_f107_");
    if (writeFiles):
       print "Connecting to DB server...."
       conn = mysql.connector.connect(host="127.0.0.1", port=3309, user="seh", db="Cluster")    
       cursor = conn.cursor()
       sqlfile = '../_disturbance_spectrum.sql'
       fd = open(sqlfile, 'r') 
       query = fd.read()
       query=query.replace("@mom","avg");
       query=query.replace("@err","stddev");
       query=query.replace("@minCount",str(minCount)).replace("@maxCount",str(maxCount))
       query=query.replace("@minDate",minDate).replace("@maxDate",maxDate);
       query=query.replace("@minR",str(minR)).replace("@maxR",str(maxR))
       query=query.replace("@minX",str(minX)).replace("@maxX",str(maxX));
       query=query.replace("@minY",str(minY)).replace("@maxY",str(maxY));
       query=query.replace("@minZ",str(minZ)).replace("@maxZ",str(maxZ));
       query=query.replace("@minLat",str(minLat)).replace("@maxLat",str(maxLat))
       if param=='dst':
          query=query.replace("O.F107","# O.F107")   # disable F10.7 filter
          query=query.replace("@minDst",str(minDst)).replace("@maxDst",str(maxDst))
       else:
          query=query.replace("O.Dst","# O.Dst")   # disable F10.7 filter
          query=query.replace("@minF107",str(minF107)).replace("@maxF107",str(maxF107))

       if (split=='Si'):
          query=query.replace(" Fe1+Fe2+Fe3","#")   # ignore Fe in WHERE statements here
       if (split=='Fe'):
          query=query.replace(" Si1+Si2+Si3","#")   # ignore Fe in WHERE statements here
# -- sepcial treatment for lobe:two separate regions
       if region=='lobe':
	 query=query.replace("#1","")
       elif region!= 'lobe':
         query=query.replace("#2","");	
       print query
       cursor.execute(query)    
       result = cursor.fetchall()
       fd.close()
       with open(outfileHi, 'wb') as fp:
          pickle.dump(result, fp)
       print "Writing SPECTRA to file..."
    else:
       with open (outfileHi, 'rb') as fp:
          result = pickle.load(fp) 
       print "Reading SPECTRA from file (!NB No DB connection)"

    numHi=0;
    SiFluxHi=[]
    SiErrHi=[]
    FeFluxHi=[]
    FeErrHi=[]
    for col in result:
        numHi=numHi+1
	for i in range(0,8):
	   SiFluxHi.append(float(col[i]))
        for i in range(8,16):
           SiErrHi.append(float(col[i]))
        for i in range(16,24):
           FeFluxHi.append(float(col[i]))
        for i in range(24,32):
           FeErrHi.append(float(col[i]))
        numHi = col[32]

#  --- Energy channels - set log10 middle -----

    ESi=[410,648,818,934,1082,1272,1520,1854,2066]
    E = ESi

    for i in range(0,8):
       logE= math.log10(E[i])+math.log10(E[i+1])
       E[i] = 10**(logE/2)
       SiFluxHi[i] = SiFluxHi[i]/E[i]
       SiErrHi[i] = SiErrHi[i]/np.sqrt(numHi)

    SiEnergy = E[:8]
#    SiSlope = np.log10(SiFlux[0]/SiFlux[5]) / np.log10(E[0]/E[5]);
    SiSlopeHi=slope(SiFluxHi,SiEnergy)

    print "SiEnergy = ", SiEnergy
    print "SiFlux   = ", SiFluxHi
    print "SiErr    = ", SiErrHi
    print "SiSlope  = ", SiSlopeHi 


    EFe=[709,1080,1344,1524,1689,1817,1961,2126,2218]
    E = EFe

    for i in range(0,8):
       logE= math.log10(E[i])+math.log10(E[i+1])
       E[i] = 10**(logE/2)
       FeFluxHi[i] = FeFluxHi[i]/E[i]
       FeErrHi[i] = FeErrHi[i]/np.sqrt(numHi)


    FeEnergy = E[:8]
#    FeSlope = np.log10(FeFlux[0]/FeFlux[5]) / np.log10(E[0]/E[5]);
    FeSlopeHi=slope(FeFluxHi,FeEnergy)

    print "FeEnergy = ",FeEnergy
    print "FeFlux   = ",FeFluxHi
    print "FeErr    = ",FeErrHi
    print "FeSlope  = ",FeSlopeHi
    print "Num points in spectra = ", numHi

    SiTextHi='Storm: (Dst < '+str(maxDst) +' nT). Slope = %.1f' % SiSlopeHi# +'  ('+str(numHi)+')'
    FeTextHi='Storm: (Dst < '+str(maxDst) +' nT). Slope = %.1f' % FeSlopeHi# +'  ('+str(numHi)+')'
    if param=='f107':
       SiTextHi='High: (F10.7 > '+str(minF107) +'). Slope = %.1f' % SiSlopeHi# +'  ('+str(numHi)+')'
       FeTextHi='High: (F10.7 > '+str(minF107) +'). Slope = %.1f' % FeSlopeHi# +'  ('+str(numHi)+')'

#  ====  NOW QUIET  / Low F10.7 ====



    outfileLo='outmatrix/disturbance_spectrum_dst_quiet'+split+'.txt'
    maxDst=99;
    minDst=0;
    minF107=0;
    maxF107=75;
 
    if param=='f107':
	outfileLo=outfileLo.replace("_dst_","_f107_");
    if (writeFiles):
       print "Connecting to DB server...."
       conn = mysql.connector.connect(host="127.0.0.1", port=3309, user="seh", db="Cluster")    
       cursor = conn.cursor()
       sqlfile = '../_disturbance_spectrum.sql'
       fd = open(sqlfile, 'r') 
       query = fd.read()
       query=query.replace("@mom","avg");
       query=query.replace("@err","stddev");
       query=query.replace("@minCount",str(minCount)).replace("@maxCount",str(maxCount))
       query=query.replace("@minDate",minDate).replace("@maxDate",maxDate);
       query=query.replace("@minR",str(minR)).replace("@maxR",str(maxR))
       query=query.replace("@minX",str(minX)).replace("@maxX",str(maxX));
       query=query.replace("@minY",str(minY)).replace("@maxY",str(maxY));
       query=query.replace("@minZ",str(minZ)).replace("@maxZ",str(maxZ));
       query=query.replace("@minLat",str(minLat)).replace("@maxLat",str(maxLat))
       query=query.replace("@minDst",str(minDst)).replace("@maxDst",str(maxDst))
       if param=='dst':
          query=query.replace("O.F107","# O.F107")   # disable F10.7 filter
          query=query.replace("@minDst",str(minDst)).replace("@maxDst",str(maxDst))
       else:
          query=query.replace("O.Dst","# O.Dst")   # disable F10.7 filter
          query=query.replace("@minF107",str(minF107)).replace("@maxF107",str(maxF107))

# -- sepcial treatment for lobe:two separate regions
       if region=='lobe':
	 query=query.replace("#1","")
       elif region!= 'lobe':
         query=query.replace("#2","");	
       print query
       cursor.execute(query)    
       result = cursor.fetchall()
       fd.close()
       with open(outfileLo, 'wb') as fp:
          pickle.dump(result, fp)
       print "Writing SPECTRA to file..."
    else:
       with open (outfileLo, 'rb') as fp:
          result = pickle.load(fp) 
       print "Reading SPECTRA from file (!NB No DB connection)"



    numLo=0;
    SiFluxLo=[]
    SiErrLo=[]
    FeFluxLo=[]
    FeErrLo=[]
    for col in result:
        numLo=numLo+1
	for i in range(0,8):
	   SiFluxLo.append(float(col[i]))
        for i in range(8,16):
           SiErrLo.append(float(col[i]))
        for i in range(16,24):
           FeFluxLo.append(float(col[i]))
        for i in range(24,32):
           FeErrLo.append(float(col[i]))
        numLo = col[32]




#  --- Energy channels - set log10 middle -----

    ESi=[410,648,818,934,1082,1272,1520,1854,2066]
    E = ESi

    for i in range(0,8):
       logE= math.log10(E[i])+math.log10(E[i+1])
       E[i] = 10**(logE/2)
       SiFluxLo[i] = SiFluxLo[i]/E[i]
       SiErrLo[i] = SiErrLo[i]/np.sqrt(numLo)

    SiEnergy = E[:8]
#    SiSlope = np.log10(SiFlux[0]/SiFlux[5]) / np.log10(E[0]/E[5]);
    SiSlopeLo=slope(SiFluxLo,SiEnergy)

    print "SiEnergy = ", SiEnergy
    print "SiFlux   = ", SiFluxLo
    print "SiErr    = ", SiErrLo
    print "SiSlope  = ", SiSlopeLo 


    EFe=[709,1080,1344,1524,1689,1817,1961,2126,2218]
    E = EFe

    for i in range(0,8):
       logE= math.log10(E[i])+math.log10(E[i+1])
       E[i] = 10**(logE/2)
       FeFluxLo[i] = FeFluxLo[i]/E[i]
       FeErrLo[i] = FeErrLo[i]/np.sqrt(numLo)


    FeEnergy = E[:8]
#    FeSlope = np.log10(FeFlux[0]/FeFlux[5]) / np.log10(E[0]/E[5]);
    FeSlopeLo=slope(FeFluxLo,FeEnergy)

    print "FeEnergy = ",FeEnergy
    print "FeFlux   = ",FeFluxLo
    print "FeErr    = ",FeErrLo
    print "FeSlope  = ",FeSlopeLo
    print "Num points in spectra = ", numLo

    SiTextLo='Quiet (Dst > '+str(minDst)+'nT): slope = %.1f' % SiSlopeLo# +'  ('+str(numLo)+')'
    FeTextLo='Quiet:(Dst > '+str(minDst)+'nT): slope = %.1f' % FeSlopeLo# +'  ('+str(numLo)+')'
    if param=='f107':
       SiTextLo='Quiet (F10.7 < '+str(maxF107)+'): slope = %.1f' % SiSlopeLo# +'  ('+str(numLo)+')'
       FeTextLo='Quiet:(F10.7 < '+str(maxF107)+'): slope = %.1f' % FeSlopeLo# +'  ('+str(numLo)+')'

# plot settings

    msize=6		# markers in mass lines 
    fsize=27		# font size

    xtextpos = 0.03 	# x-position of text /annotation
    ytextpos = 0.99	# y-position of top text /annotation
    deltaytext=0.07     # y-spacing between text lines

 
    plt.figure(figsize=(28,9))
    plt.rcParams.update({'font.size': fsize})
    plt.subplots_adjust(hspace=0.3)

    spectrumSi = plt.subplot(121) 


    spectrumSi.set_xlabel('Energy [keV]')
    spectrumSi.set_ylabel('Si Intensity [1/(keV*s)]')
    spectrumSi.set_yscale('log');
    spectrumSi.set_xscale('log');
    spectrumSi.set_xlim(4e2,2.5e3);
    spectrumSi.set_ylim(20e-7,8e-2);
    spectrumSi.grid(True,which='both')
    spectrumSi.set_xticks([500,1000,2000])
    spectrumSi.get_xaxis().set_major_formatter(ml.ticker.ScalarFormatter())

    spectrumSi.errorbar(SiEnergy,SiFluxHi, marker='o',color='blue',markersize=10,linewidth=2)
    spectrumSi.errorbar(SiEnergy,SiFluxLo, marker='o',color='darkblue',markersize=10,linewidth=2)
 
    if param=='dst':
       spectrumSi.text(0.03,ytextpos-1*deltaytext,'Si, geomagnetic activity:',color='blue',ha='left',transform=spectrumSi.transAxes)#    E_TOF.text(xtextpos,ytextpos-3*deltaytext,'Ygse = ['+str(minY)+" - "+ str(maxY)+"]",transform=E_TOF.transAxes)        
    else:
       spectrumSi.text(0.03,ytextpos-1*deltaytext,'Si, solar activity:',color='blue',ha='left',transform=spectrumSi.transAxes)#    E_TOF.text(xtextpos,ytextpos-3*deltaytext,'Ygse = ['+str(minY)+" - "+ str(maxY)+"]",transform=E_TOF.transAxes)        
    
    spectrumSi.text(0.03,ytextpos-2*deltaytext,SiTextHi,color='blue',ha='left',transform=spectrumSi.transAxes)#    E_TOF.text(xtextpos,ytextpos-3*deltaytext,'Ygse = ['+str(minY)+" - "+ str(maxY)+"]",transform=E_TOF.transAxes)  
    spectrumSi.text(0.03,ytextpos-12*deltaytext,SiTextLo,color='darkblue',ha='left',transform=spectrumSi.transAxes)

    spectrumSi.text(-0.16, 1.02,label[0],fontsize=36, color='blue',transform=spectrumSi.transAxes);
    spectrumSi.text(1.09, 1.02,label[1],fontsize=36, color='blue',transform=spectrumSi.transAxes);

    spectrumFe = plt.subplot(122) 
    spectrumFe.set_xlabel('Energy [keV]')
    spectrumFe.set_ylabel('Fe Intensity [1/(keV*s)]')
    spectrumFe.set_yscale('log');
    spectrumFe.set_xscale('log');
    spectrumFe.set_xlim(4e2,2.5e3);
    spectrumFe.set_ylim(20e-7,8e-2);
    spectrumFe.grid(True,which='both')
    spectrumFe.yaxis.tick_right()
    spectrumFe.yaxis.set_label_position("right")
    spectrumFe.set_xticks([500,1000,2000])
    spectrumFe.get_xaxis().set_major_formatter(ml.ticker.ScalarFormatter())

    spectrumFe.errorbar(FeEnergy,FeFluxHi, marker='o',color='red',markersize=10,linewidth=2)
    spectrumFe.errorbar(FeEnergy,FeFluxLo, marker='o',color='darkred',markersize=10,linewidth=2)
         
    if param=='dst':
       spectrumFe.text(0.03,ytextpos-1*deltaytext,"Fe, geomagnetic activity:",color='red',ha='left',transform=spectrumFe.transAxes)#    E_TOF.text(xtextpos,ytextpos-3*deltaytext,'Ygse = ['+str(minY)+" - "+ str(maxY)+"]",transform=E_TOF.transAxes)
    else:
       spectrumFe.text(0.03,ytextpos-1*deltaytext,"Fe, solar activity:",color='red',ha='left',transform=spectrumFe.transAxes)#    E_TOF.text(xtextpos,ytextpos-3*deltaytext,'Ygse = ['+str(minY)+" - "+ str(maxY)+"]",transform=E_TOF.transAxes)

    spectrumFe.text(0.03,ytextpos-2*deltaytext,FeTextHi,color='red',ha='left',transform=spectrumFe.transAxes)#    E_TOF.text(xtextpos,ytextpos-3*deltaytext,'Ygse = ['+str(minY)+" - "+ str(maxY)+"]",transform=E_TOF.transAxes)
    spectrumFe.text(0.03,ytextpos-13*deltaytext,FeTextLo,color='darkred',ha='left',transform=spectrumFe.transAxes)
 

    plotname="Disturbance_spectrum_@param.png"
    plotname=plotname.replace("@param",param);
    plt.savefig("../../Pap/figs/"+plotname, bbox_inches='tight',dpi=300)
    print plotname


    if len(sys.argv)>=1:
	plt.show()


def slope(y,x):
    s= np.polyfit(np.log10(x[:8]),np.log10(y[:8]),1)
    return (s[-2])
 

def roundup(x):
    return int(math.ceil(x / 1.0)) * 1

main()

