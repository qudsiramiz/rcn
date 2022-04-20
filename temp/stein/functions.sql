#
# User defined functions and definitions
#


# Database to use;
USE MMS;


# -- useful constants: proton and electron masses  ---
set @mp=1.6726219e-27;
set @me=9.10938356e-31;

# ------------------------------------------------------------
#
#
# returns deltaVstar = see Phan et al, 1996b,
# input VL1, VM1 = V components at t1 [km/s]
# input VL1, VM1 = V components at t2 [km/s]# 
# input BL1, BM1 = B-field components at t1 [nT]
# input BL2, BM2 = B-field components at t2 [nT]
# input N1, N2 = number densities at t1 and t2 resp [cm^-3]
#
#  Formula from Phan et al, 1996b
#
# ------------------------------------------------------------
drop function if exists DeltaVstar;

delimiter $$
create function DeltaVstar(V1L float, V1M float, VmaxL float, VmaxM float, V2L float, V2M float, 
                           B1L float, B1M float, BmaxL float, BmaxM float, B2L float, B2M float,
                           N1 float, Nmax float, N2 float)
returns float
begin
  declare VA1L float default VA(B1L, N1);
  declare VA2L float default VA(B2L, N2);
  declare VAmaxL float default VA(BmaxL, Nmax);

  declare VA1M float default VA(B1M, N1);
  declare VA2M float default VA(B2M, N2);
  declare VAmaxM float default VA(BmaxM, Nmax);

#inbound 
  declare dVL float default (VmaxL - V1L);
  declare dVAL float default (VAmaxL - VA1L);

  declare dVM float default (VmaxM - V1M);
  declare dVAM float default (VAmaxM - VA1M);

  declare dv float default 0;   
  declare VAsquared float default 0;

#outbound
  if (B2L > B1L) then
     set dVL = (VmaxL - V2L); 
     set dVAL = (VAmaxL - VA2L);
     set dVM  = (VmaxM - V2M);
     set dVAM = (VAmaxM - VA2M);
  end if;


  set dv = dot(dVL,dVM,0,  dVAL,dVAM,0);   
  set VAsquared = dot(dVAL,dVAM,0,  dVAL,dVAM,0);
  return (dv/VAsquared);

end $$
delimiter ;




# ----------------------------------------
#
# Functions from GP mail, 15 Jun 2017 
#
#
# ----------------------------------------

-- Create syntax for FUNCTION 'beta'
drop function if exists beta;
CREATE  FUNCTION `beta`(N float, Tpara float, Tperp float, BL float, BM float, BN float) RETURNS float
return 0.4021*N*T(Tpara,Tperp)/(BL*BL+BM*BM+BN*BN);

-- Create syntax for FUNCTION 'clock'
drop function if exists clock;
CREATE  FUNCTION `clock`(By_gse float, Bz_gse float) RETURNS float
return 57.3*acos(By_gse/sqrt(By_gse*By_gse+Bz_gse*Bz_gse));

-- Create syntax for FUNCTION 'deltat'
drop function if exists deltat;
CREATE  FUNCTION `deltat`(t1 char(25), t2 char(25)) RETURNS float
return time_to_sec(timediff(t1,t2))+microsecond(timediff(t1,t2))/1e6;

-- Create syntax for FUNCTION 'dot'
drop function if exists dot;
CREATE  FUNCTION `dot`(x1 float, y1 float, z1 float, x2 float, y2 float, z2 float) RETURNS float
return (x1*x2 + y1*y2 + z1*z2);

-- Create syntax for FUNCTION 'DT'
drop function if exists DT;
CREATE  FUNCTION `DT`(t1 char(25), t2 char(25)) RETURNS float
return time_to_sec(timediff(t1,t2))+microsecond(timediff(t1,t2))/1e6;

-- Create syntax for FUNCTION 'evr'
drop function if exists evr;
CREATE  FUNCTION `evr`(l2 float, l1 float) RETURNS float
return l2/l1;

-- Create syntax for FUNCTION 'lambda'
drop function if exists lambda;
CREATE  FUNCTION `lambda`(N float) RETURNS float
return 227/sqrt(N);

-- Create syntax for FUNCTION 'mag'
drop function if exists mag;
CREATE  FUNCTION `mag`(x float, y float, z float) RETURNS float
return sqrt(x*x+y*y+z*z);

-- Create syntax for FUNCTION 'shear'
drop function if exists shear;
CREATE  FUNCTION `shear`(x1 float, y1 float, z1 float, x2 float, y2 float, z2 float) RETURNS float
return 57.3*acos(dot(x1,y1,z1,x2,y2,z2)/(mag(x1,y1,z1)*mag(x2,y2,z2)));

-- Create syntax for FUNCTION 'T'
drop function if exists T;
CREATE  FUNCTION `T`(Tpara float, Tperp float) RETURNS float
return (Tpara + 2*Tperp)/3;

-- Create syntax for FUNCTION 'VA'
drop function if exists VA;
CREATE  FUNCTION `VA`(B float, N float) RETURNS float
return 21.83*B/sqrt(N); 


# ------------------------------------------------------------
#
#
# returns rgp = the gyro radius for protons [km]
# input t = perpendicular temp [eV]
# input x,y,z = B-field components [nT]
#
#  Formula from QSAS Plasmaparams r=1300*sqrt(T)/|B|, where T is in MK
#
#
# ------------------------------------------------------------
drop function if exists rgyrop;

create function rgyrop(t float, x float, y float, z float)
returns float
return 140*sqrt(t)/sqrt(x*x+y*y+z*z)
;





# ---------------------------------
#
# returns the shear angle between two
# vectors in LMN representation. N is ignored
# Positive M value of the 2nd argument gives 
# a positive angle 
#
# ---------------------------------
drop function if exists shear;

CREATE FUNCTION `shear`(x1 float, y1 float, z1 float, x2 float, y2 float, z2 float) 
RETURNS float
return 57.3*acos(dot(x1,y1,z1,x2,y2,z2)/(mag(x1,y1,z1)*mag(x2,y2,z2)))
; 

# ---------------------------------
#
# checks whether a crossing is inbound 
# or outbound on basis of BL 
#
# ---------------------------------
drop function if exists isInbound; 

create function isInbound(BL1 float, BL2 float)
returns boolean
return (BL2 > BL1)
;


# ---------------------------------
#
# fills the FlagStr field
#
# ---------------------------------
drop function if exists classifyXing;

delimiter $$
create function classifyXing(BL1 float,BL2 float, HTcc float, ypos float) 
returns char(25)
begin 
  declare flag char(25) default "un,o,nmon,bf,mult,part";

  if (HTcc > 0.9 and abs(ypos) < 10) then
     set flag = replace(flag,"un,","mp,");  
  end if;

  if (BL2 > BL1) then
     set flag = replace(flag,',o,',',i,');
  end if;

  return flag;

end $$
delimiter ;



# ---------------------------------
#
# Calculates then[IMF] clock angle
#
# ---------------------------------
drop function if exists clockAngle;

delimiter $$
create function clockAngle(y float, z float) returns float
begin 
  declare rad2deg float default 57.296;
  if (y < 0) then 
     set rad2deg =  -57.296;
  end if;

  return rad2deg * acos(z/sqrt(y*y+z*z));
end $$
delimiter ;




# ---------------------------------
#
#
#
# ---------------------------------
drop function if exists deltat; 

create function deltat(t1 char(25), t2 char(25))
returns float
return time_to_sec(timediff(t1,t2))+microsecond(timediff(t1,t2))/1e6;
;


# ---------------------------------
#
#  returns magnitude of a vector;
#
#  inputs : x,y,z componenst of vector
#
# ---------------------------------
drop function if exists mag;

create function mag (x float, y float, z float)
returns float
return sqrt(x*x+y*y+z*z)
;

# ---------------------------------
#
#  returns temperature;
#
#  inputs : Tpara and Tperp
#
# ---------------------------------
drop function if exists T;

create function T (Tpara float, Tperp float)
returns float
return (Tpara + 2*Tperp)/3
;

# ---------------------------------
#
#  returns lambda;
#
#  inputs : N
#
# ---------------------------------
drop function if exists lambda;

create function lambda (N float)
returns float
return 227/sqrt(N)
;

# ---------------------------------
#
#  extracts L,M,N coordinates from a GSE vector
#  i.e., dot vector into a LMN unit vector
#
#  inputs : 
#       x,y,z component of vector 1,
#       LMNx,LMNy,LMNz = components of the L, M or N unit vector
#
# ---------------------------------
drop function if exists LMN_from_GSE;

create function LMN_from_GSE(x float, y float, z float, LMNx float, LMNy float, LMNz float)
returns float
return (x*LMNx + y*LMNy + z*LMNz)
;




# ---------------------------------
#
#  returns Bfield  from VAlven and massdensity 
#  i.e., Bx = Vx * sqrt(N * mass). 
#  NB! check units 
#
#  inputs:
#    valfven = Alfven velocity [km/s]
#    density = number density [cm^-3]
#    mass = mass of ion or electron [kg]
#
# ---------------------------------
drop function if exists B_from_VA;

create function B_from_VA(valfven float, density float, mass float)
returns float
return valfven * sqrt(density * mass) * 1e12      # factor 1e12 to get units nT
; 



