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
from matplotlib.patches import Wedge

def main():

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

		region='all'
		label=['a)', 'b)', 'c)']
		csys='GSM'

# get from arguments instead
		if len(sys.argv)>1:
			region=sys.argv[1]

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



# -- Mysql connection ----

		outfile='outmatrix/spectrumresults_'+region+'.txt'
		if (writeFiles):
			 print ("Connecting to DB server....")
			 conn = mysql.connector.connect(host="127.0.0.1", port=3309, user="seh", db="Cluster")    
			 cursor = conn.cursor()
			 sqlfile = '../_matrx_spectrum.sql'
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
# -- sepcial treatment for lobe:two separate regions
			 if region=='lobe':
					query=query.replace("#1","")
			 elif region!= 'lobe':
					query=query.replace("#2","");  
			 print (query)
			 cursor.execute(query)    
			 result = cursor.fetchall()
			 fd.close()
			 with open(outfile, 'wb') as fp:
					pickle.dump(result, fp)
			 print ("Writing SPECTRA to file...")
		else:
			 with open (outfile, 'rb') as fp:
					result = pickle.load(fp) 
			 print ("Reading SPECTRA from file (!NB No DB connection)")


		num=0;
		SiFlux=[]
		SiErr=[]
		FeFlux=[]
		FeErr=[]
		for col in result:
				num=num+1
				for i in range(0,8):
					SiFlux.append(float(col[i]))
				for i in range(8,16):
						SiErr.append(float(col[i]))
				for i in range(16,24):
						FeFlux.append(float(col[i]))
				for i in range(24,32):
						FeErr.append(float(col[i]))
				num = col[32]



#  --- Energy channels - set log10 middle -----

		ESi=[410,648,818,934,1082,1272,1520,1854,2066]
		E = ESi

		for i in range(0,8):
			 logE= math.log10(E[i])+math.log10(E[i+1])
			 E[i] = 10**(logE/2)
			 SiFlux[i] = SiFlux[i]/E[i]
			 SiErr[i] = SiErr[i]/np.sqrt(num)

		SiErr = np.array(SiErr)
		SiErr = 0.434*SiErr/SiFlux
		SiErr = np.array([SiFlux*(1-10**(-SiErr)),SiFlux*(10**(SiErr)-1)])
	 
		SiEnergy = E[:8]
#    SiSlope = np.log10(SiFlux[0]/SiFlux[5]) / np.log10(E[0]/E[5]);
		SiSlope=slope(SiFlux,SiEnergy)

		print ("SiEnergy = ", SiEnergy)
		print ("SiFlux   = ", SiFlux)
		print ("SiErr    = ", SiErr)
		print ("SiSlope  = ", SiSlope )




# ===== now Fe ----



		EFe=[709,1080,1344,1524,1689,1817,1961,2126,2218]
		E = EFe

		for i in range(0,8):
			 logE= math.log10(E[i])+math.log10(E[i+1])
			 E[i] = 10**(logE/2)
			 FeFlux[i] = FeFlux[i]/E[i]
			 FeErr[i] = FeErr[i]/np.sqrt(num)
 
		FeErr = np.array(FeErr)
		FeErr = 0.434*FeErr/FeFlux
		FeErr = np.array([FeFlux*(1-10**(-FeErr)),FeFlux*(10**(FeErr)-1)])


		FeEnergy = E[:8]
#    FeSlope = np.log10(FeFlux[0]/FeFlux[5]) / np.log10(E[0]/E[5]);
		FeSlope=slope(FeFlux,FeEnergy)

		print ("FeEnergy = ",FeEnergy)
		print ("FeFlux   = ",FeFlux)
		print ("FeErr    = ",FeErr)
		print ("FeSlope  = ",FeSlope)
		print ("Num points in spectra = ", num)

#
#  --- panel b), c) show positions
#


		outfile='outmatrix/positionresults_'+region+'.txt'
		if (writeFiles):
			 print ("Connecting to DB server for positions...")
			 sqlfile = '../_matrx_position.sql'
			 fd = open(sqlfile, 'r') 
			 query = fd.read()
			 query=query.replace("@minCount",str(minCount)).replace("@maxCount",str(maxCount))
			 query=query.replace("@minDate",minDate).replace("@maxDate",maxDate);
			 query=query.replace("@minR",str(minR)).replace("@maxR",str(maxR))
			 query=query.replace("@minX",str(minX)).replace("@maxX",str(maxX));
			 query=query.replace("@minY",str(minY)).replace("@maxY",str(maxY));
			 query=query.replace("@minZ",str(minZ)).replace("@maxZ",str(maxZ));
			 query=query.replace("@minLat",str(minLat)).replace("@maxLat",str(maxLat))
# -- sepcial treatment for lobe:two separate regions
			 if region=='lobe':
					query=query.replace("#1",'')
			 elif region!='lobe':
				 query=query.replace("#2","");  
# -- GSE in SW, GSM else ---
			 if region=='sw':
					query=query.replace("#3","")
			 elif region!= 'sw':
					query=query.replace("#4","");  

			 print (query)
			 cursor.execute(query)    
			 result = cursor.fetchall()
			 fd.close()
			 with open(outfile, 'wb') as fp:
					pickle.dump(result, fp)
			 print ("Writing POS to file....")
			 print (query)
		else:
			 with open (outfile, 'rb') as fp:
					result = pickle.load(fp) 
			 print ("Reading POS from file (no DB connected!)")

		xpos=[]
		ypos=[]
		zpos=[]
		color=[]   # 'black' if bot counts
		for col in result:
			xpos.append(float(col[0]))
			ypos.append(float(col[1]))  
			zpos.append(float(col[2]))
				





 # plot settings
		msize=6   # markers in mass lines 
		fsize=27    # font size

		xtextpos = 0.03   # x-position of text /annotation
		ytextpos = 0.99 # y-position of top text /annotation
		deltaytext=0.07     # y-spacing between text lines

 
		plt.figure(figsize=(36,7))
		plt.rcParams.update({'font.size': fsize})

		spectrum = plt.subplot(131) 
		spectrum.set_xlabel('Energy [keV]')
		spectrum.set_ylabel('Intensity [1/(keV*s)]')
		spectrum.set_yscale('log');
		spectrum.set_xscale('log');
		spectrum.set_xlim(4e2,2.5e3);
		spectrum.set_ylim(11e-7,8e-2);
		spectrum.grid(True,which='both')

		spectrum.set_xticks([500,1000,2000])
		spectrum.get_xaxis().set_major_formatter(ml.ticker.ScalarFormatter())

		spectrum.errorbar(FeEnergy,FeFlux, marker='o',color='red',markersize=20,linewidth=2)
		spectrum.errorbar(SiEnergy,SiFlux, marker='o',color='blue',markersize=20,linewidth=2)

#    ax2 = spectrum.twinx()
#    ax2.errorbar(SiEnergy,SiFlux,yerr=SiErr,marker='o',linewidth=0,markersize=60)
#    ax2.set_xlim(4e2,2.5e3);
#    ax2.set_ylim(7e-7,8e-2);

		fac=20;
		lw = 2
		msize=10
		mpy=np.linspace(-25,20);
		mpx=-pow(mpy,2)/fac
		cm1 = plt.cm.get_cmap('jet')



		XY=plt.subplot(132) 
#    XY.set_title("XY$_{GSE}$ projection")
		XY.set_xlabel('X$_{'+csys+'}$ [Re]')
		XY.set_ylabel('Y$_{'+csys+'}$ [Re]')
		XY.set_xlim(21,-21)
		XY.set_ylim(21,-21)
		XY.plot([-21,21],[0,0],'k-',linewidth=lw)
		XY.plot([0,0],[-21,21],'k-',linewidth=lw)
		XY.plot(mpx+11,mpy,'k-',linewidth=lw)
		if (region=='sw' or region=='lobe' or region=='ps'):
			earth(angle=270,zorder=10)
		XY.scatter(xpos,ypos,edgecolors='none', color='black', marker="h",s=msize,cmap=cm1)
		XY.set_aspect(1)
		XY.grid(True)    


		XZ=plt.subplot(133) 
#    XZ.set_title("XZ$_{GSE}$ projection")
		XZ.set_xlabel('X$_{'+csys+'}$ [Re]')
		XZ.set_ylabel('Z$_{'+csys+'}$ [Re]')
		XZ.set_xlim(21,-21)
		XZ.set_ylim(-21,21)
		XZ.plot([-21,21],[0,0],'k-',linewidth=lw)
		XZ.plot([0,0],[-21,21],'k-',linewidth=lw)
		XZ.plot(mpx+11,mpy,'k-',linewidth=lw)
		if (region=='sw' or region=='cusp' or region=='lobe' or region=='ps'):
			earth(angle=270,zorder=10)
		XZ.scatter(xpos,zpos,edgecolors='none', color='black',marker="h",s=msize,cmap=cm1)
		XZ.set_aspect(1)
		XZ.grid(True);



#    spectrum.text(0.95,ytextpos-2*deltaytext,'Date = ['+str(minDate)+" - "+ str(maxDate)+"]",ha='right',transform=spectrum.transAxes)
		spectrum.text(0.95,ytextpos-deltaytext,ptit+' ('+str(num)+' samples)',ha='right', transform=spectrum.transAxes)
		spectrum.text(0.95,ytextpos-2*deltaytext,'Si, slope = '+"%.1f" % SiSlope,color='blue',ha='right',transform=spectrum.transAxes)
		spectrum.text(0.95,ytextpos-3*deltaytext,'Fe, slope = '+"%.1f" % FeSlope,color='red',ha='right',transform=spectrum.transAxes)#    E_TOF.text(xtextpos,ytextpos-3*deltaytext,'Ygse = ['+str(minY)+" - "+ str(maxY)+"]",transform=E_TOF.transAxes)
#    E_TOF.text(xtextpos,ytextpos-4*deltaytext,'Zgse = ['+str(minZ)+" - "+ str(maxZ)+"]",transform=E_TOF.transAxes)
#    E_TOF.text(xtextpos,ytextpos-5*deltaytext,'Dst = ['+str(minDst)+" - "+ str(maxDst)+"]",transform=E_TOF.transAxes)

		spectrum.text(-0.2, 1.02,label[0],fontsize=36, color='blue',transform=spectrum.transAxes);
		spectrum.text(1.15, 1.02,label[1],fontsize=36, color='blue',transform=spectrum.transAxes);
		spectrum.text(2.35, 1.02,label[2],fontsize=36, color='blue',transform=spectrum.transAxes);

		xposstr="_x"+str(minX)+"_"+str(maxX)
		yposstr="_y"+str(minY)+"_"+str(maxY)
		zposstr="_z"+str(minZ)+"_"+str(maxZ)
		rposstr="_r"+str(minR)+"_"+str(maxR)

		plotname="Region_"+region+xposstr+yposstr+zposstr+rposstr
		plotname=plotname.replace(".","")+".png"
#		plt.savefig("../../Pap/figs/"+plotname, bbox_inches='tight',dpi=300)
		print (plotname)


		if len(sys.argv)>=1:
			plt.show()


def slope(y,x):
		s= np.polyfit(np.log10(x[:8]),np.log10(y[:8]),1)
		return (s[-2])
 

def roundup(x):
		return int(math.ceil(x / 1.0)) * 1

#  -----------------------------------------------------------------------
#
#  Draws the Earth with day/nighside
#
#  -----------------------------------------------------------------------
def earth(angle=0, ax=None, colors=('w','k'), **kwargs):
		"""
		Add two half circles to the axes *ax* (or the current axes) at the lower
		left corner of the axes with the specified facecolors *colors* rotated at
		*angle* (in degrees).
		"""
		radius=1.0;
		if ax is None:
				ax = plt.gca()
		center = (0, 0)
		theta1, theta2 = angle, angle + 180
		w1 = Wedge(center, radius, theta1, theta2, fc=colors[0], **kwargs)
		w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], **kwargs)
		for wedge in [w1, w2]:
				ax.add_artist(wedge)
		return [w1, w2]



main()

