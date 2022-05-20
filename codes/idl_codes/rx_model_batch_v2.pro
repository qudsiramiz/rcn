
function getshear, bx0,by0,bz0,bx1,by1,bz1
  dp = (bx0*bx1+by0*by1+bz0*bz1)
  mag0 = sqrt(bx0^2+by0^2+bz0^2)
  mag1 = sqrt(bx1^2+by1^2+bz1^2)
  angle = acos(dp/(mag0*mag1))*180./!PI
  return, angle
end

function getrxben, bx0,by0,bz0,bx1,by1,bz1
  mag0 = sqrt(bx0^2+by0^2+bz0^2)
  mag1 = sqrt(bx1^2+by1^2+bz1^2)
  dp = (bx0*bx1+by0*by1+bz0*bz1)
  hat0 = [bx0,by0,bz0]/mag0
  hat1 = [bx1,by1,bz1]/mag1
  dtp = dotp(hat0,-hat1)
  ;  angle = acos(dp/(mag0*mag1))*180./!PI
  ;  b0 = [bx0,by0,bz0]
  ;  b0hat = b0/norm(b0)
  ;  bisector = mag0*[bx1,by1,bz1] + mag1*[bx0,by0,bz0]
  ;  u_bisect = bisector / norm(bisector)
  ;  rxbmag0 = dotp(u_bisect,[bx0,by0,bz0])
  ;  rxbmag1 = dotp(u_bisect,[bx1,by1,bz1])
  return, 0.5*(mag0+mag1)*(1+dtp)
  ;return, (rxbmag0^2 + rxbmag1^2)*1.03 ; MJ/RE^3
  ;return, (rxbmag0^2 + rxbmag1^2)*3.98*10^(-4);nPa
end

function getvcs, bx0,by0,bz0,bx1,by1,bz1,n0,n1
  va_pl = 21.812 ; conv. nT, m_P/cm^3 product to km/s cassak-shay
  mag0 = sqrt(bx0^2+by0^2+bz0^2)
  mag1 = sqrt(bx1^2+by1^2+bz1^2)
  hat0 = [bx0,by0,bz0]/mag0
  hat1 = [bx1,by1,bz1]/mag1
  dtp = dotp(hat0,-hat1)
  ;angle = acos(dp/(mag0*mag1))*180./!PI
  bisector = mag0*[bx1,by1,bz1] + mag1*[bx0,by0,bz0]
  u_bisect = bisector / norm(bisector)
  ;rxbmag0 = dotp(u_bisect,[bx0,by0,bz0])
  ;rxbmag1 = dotp(u_bisect,[bx1,by1,bz1])
  rxbmag0 = mag0*(1+dtp)/2
  rxbmag1 = mag1*(1+dtp)/2
  return, va_pl * sqrt(rxbmag0*rxbmag1 * (rxbmag0+rxbmag1)/(rxbmag0*n1 + rxbmag1*n0))
end

function getbis, bx0,by0,bz0,bx1,by1,bz1
  mag0 = sqrt(bx0^2+by0^2+bz0^2)
  mag1 = sqrt(bx1^2+by1^2+bz1^2)
  hat0 = [bx0,by0,bz0]/mag0
  hat1 = [bx1,by1,bz1]/mag1
  dtp = dotp(hat0,-hat1)
  ;  angle = acos(dp/(mag0*mag1))*180./!PI
  ;  b0 = [bx0,by0,bz0]
  ;  b0hat = b0/norm(b0)
  ;  bisector = mag0*[bx1,by1,bz1] + mag1*[bx0,by0,bz0]
  ;  u_bisect = bisector / norm(bisector)
  ;  rxbmag0 = dotp(u_bisect,[bx0,by0,bz0])
  ;  rxbmag1 = dotp(u_bisect,[bx1,by1,bz1])
  return, -ACOS(dtp)*180./!dpi
  ;return, (rxbmag0^2 + rxbmag1^2)*1.03 ; MJ/RE^3
  ;return, (rxbmag0^2 + rxbmag1^2)*3.98*10^(-4);nPa
end

function getca, bx0, by0, bz0
  return, ATAN(by0/bz0)*180./!dpi
end

function getmag, bx0, by0, bz0
  return, sqrt(bx0*bx0 + by0*by0 + bz0+bz0)
end

; Brian Walsh
; 15 Feb 2012
; Made to show shear angle between draped imf from the cooling et al 2001 model
; and the internal geomagnetic field from the sun's point of view.
; OMNI_HRO_1min_BX_GSE
PRO RX_model_batch_v2, probe=probe,maximum_shear=maximum_shear, movie=movie, mmsprobe=mmsprobe, times=['2016-12-24 15:08', '2016-12-24 15:12']


  ;for jeff, probe indicates which THEMIS probe to project on the shear distribution map
  ;maximum_shear means to overlay MSM. If not set then no MSM just the shear distribution
  ;movie produces snapshots of map and complies into one movie. If not set then only one snapshot is produced
  ;    omni_init

  Syear=''
  Smonth=''
  Sday=''
  Shour=''
  Ssmin=''
  Ssec=''
  cdum=''
  ;2C3uw79hYAj3xgDU.

  ;    cd, '~/HSRY1'

  if keyword_set(times) then begin
    t1 = times[0]
    dt = time_double(times[1]) - time_double(times[0])
    timespan, t1, dt, /seconds
    omni_hro_load,trange=times
    t = [time_double(times[0]),time_double(times[1])]

    reads,time_string(t[0]),Syear,cdum,Smonth,cdum,Sday,cdum,Shour,cdum,Ssmin,cdum,Ssec,format='(a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a2)'
  endif
  ;load sw and themis s/c info
  ;omni_load,/init,/median

  if keyword_set(probe) then begin
    thm_load_state,probe=probe,trange=[t[0]-30,t[1]+30]
    for i=0, n_elements(probe)-1 do begin
      cotrans,'th'+probe[i]+'_state_pos','th'+probe[i]+'_state_pos_gse',/gei2gse
      cotrans,'th'+probe[i]+'_state_pos_gse','th'+probe[i]+'_state_pos_gsm',/gse2gsm
      sc_pos = 'th'+probe[i]+'_state_pos_gsm'
      get_data,sc_pos,data=pos1
    endfor
  endif

  if keyword_set(mmsprobe) then begin
    ;mms_init

    mms_load_state,probe=mmsprobe,trange=[t[0]-300,t[1]+300]
    sc_pos = 'mms1_mec_r_gsm'
    get_data,sc_pos,data=pos1
    print, 'got MMS probe data, y=' + string(pos1.Y[0,1]/6371) + ', z= ' + string(pos1.Y[0,2]/6371)
  endif
  omni_prefix = 'OMNI_HRO_1min_'
  if not keyword_set(movie) then begin
    time_clip,omni_prefix+'BX_GSE',time_string(times[0]),time_string(times[1]), newname='BX_GSE_tclip';lag 10 min to account for SW propagation within SH
    time_clip,omni_prefix+'BY_GSM',time_string(times[0]),time_string(times[1]), newname='BY_GSM_tclip'
    time_clip,omni_prefix+'BZ_GSM',time_string(times[0]),time_string(times[1]), newname='BZ_GSM_tclip'
    get_data,'BX_GSE_tclip',data=Bimfx_st
    get_data,'BY_GSM_tclip',data=Bimfy_st
    get_data,'BZ_GSM_tclip',data=Bimfz_st
    get_data,omni_prefix+'Vx',data=Vx
    get_data,omni_prefix+'proton_density',data=Nsw
    get_data,omni_prefix+'SYM_H',data=SYM_H ; symmetric (SYM) disturbance index in horizontal
    ; direction

    bx = median(Bimfx_st.y)
    by = median(Bimfy_st.y)
    bz = median(Bimfz_st.y)
    Bimfx=bx
    Bimfy=by
    Bimfz=bz
    fmt = '(f6.2)'
    title = 'Bimf [x,y,z]: '+string(bx,format=fmt)+', '+string(by,format=fmt)+', '+string(bz,format=fmt)

    time_clip,omni_prefix+'Vx',time_string(t[0]),time_string(t[1]), newname='Vx_tclip'
    time_clip,omni_prefix+'proton_density',time_string(t[0]),time_string(t[1]), newname='Nsw_tclip'
    time_clip,omni_prefix+'SYM_H',time_string(t[0]-30),time_string(t[1]+30),newname='SYM_H_tclip' ; no lag
    get_data,'Vx_tclip',data=Vx
    get_data,'Nsw_tclip',data=Nsw
    get_data,'SYM_H_tclip',data=SYM_H
    vel=median(Vx.y)
    n = median(Nsw.y)
    dst=median(SYM_H.y)

    if not finite(vel) then vel=450.
    if not finite(n) then n=4.
    if not finite(dst) then dst=-1

    mp = 1.67262158e-27 ; mass of proton in kg
    ;rho =; kg/m^3
    rho  = 1.94d-6 * n * vel^2
    pdyn = 1.94d-6 * n * vel^2

    param = [rho,dst,by,bz,0,0,0,0,0,0]
    reads,time_string(t[0]),Syear,cdum,Smonth,cdum,Sday,cdum,Shour,cdum,Ssmin,cdum,Ssec,format='(a4,a1,a2,a1,a2,a1,a2,a1,a2,a1,a2)'
    geopack_recalc,fix(Syear),fix(Smonth),fix(Sday),fix(Shour),fix(Ssmin),0,/date ; Setup Geopack global variables.
    geopack_epoch,epo,fix(Syear),fix(Smonth),fix(Sday),fix(Shour),fix(Ssmin),0,/compute_epoch ; Convert between UTC and CDF Epoch time formats.
    ;geopack_getw,n,vel,bz,ww   ; Calculate G parameters for T01 model. Outputs: w: array containing G parameters for each time step.

    v=150
    shear = fltarr(v,v)
    rx_en = fltarr(v,v)
    va_cs = fltarr(v,v)
    bisec = fltarr(v,v)
    ycoord = fltarr(v,v)
    zcoord = fltarr(v,v)
    bshca = fltarr(v,v)
    bshmg = fltarr(v,v)
    nsh = fltarr(v,v)

    ro=(10.22+1.29*tanh(0.184*(bz+8.14)))*(Pdyn)^(-1.0/6.6) ; Shue et al., 1998
    rostring = strtrim(string(ro,format='(F5.2)'),1)
    rostring = rostring.replace('.','_')

    alpha=(0.58-0.007*bz)*(1+0.024*ALOG(Pdyn))
    rmp = ro*(2/(1+COS(0.0)))^alpha  ; Stand off position of the magnetopause
    dtheta = !PI/100.

    ;;you may want to control the thickness of the generated plots
    ;;to make the lines in the output more visible when exported
    ;;you can do so with these variables
    ;axisthick = 4.0
    ;charthick = 6.0
    ;thick = 10.0
    ;charsize = 1.0
    ;symsize = 1.0
    ;!X.STYLE=1
    ;!Y.STYLE=1
    ;!X.MARGIN = 7
    ;HSIZE = 100
    ;popen,'B_subsolar_shear_10min_'+Syear+Smonth+Sday+Shour+Ssmin,encapsulated=1,land=0,xsize=6,ysize=6
    loadct,39
    ;fmt = '(f6.2)'
    ;title = 'Bimf [x,y,z]: '+string(bx,format=fmt)+', '+string(by,format=fmt)+', '+string(bz,format=fmt)

    ; At each distance down the tail from x=8.75 to x=-60 with 1 Re cadience get the B vector at the equator
    a=2.
    leny = 80
    lenz = 80
    for ynum = 0, leny-1 do begin
      y0 = float(ynum)-40 ;40.-float(ynum)
      for znum = 0, lenz-1 do begin
        ; Find caps at bottom and top for each x value
        z0 = float(znum)-40 ;40.-float(znum)

        rp = sqrt(z0^2+y0^2) ; projection of r into the yz plane

        ; Find the theta value that corresponds with this point on the Shue et al 1998 magnetopause
        for index = 0, 100 do begin
          theta = float(index)*dtheta
          r = ro*(2/(1+COS(theta)))^alpha
          zp = r*sin(theta) ;not really in z direction, but a distance in yz plane
          x0 = r*cos(theta)
          signx = x0/abs(x0)
          signy = y0/abs(y0)
          signz = z0/abs(z0)
          if x0 eq 0.0 then signx = 1.0
          if y0 eq 0.0 then signy = 1.0
          if z0 eq 0.0 then signz = 1.0

          if rp le zp then begin
            GOTO, JUMP1
          endif
        endfor

        JUMP1: if rp le zp then begin
          ycoord[ynum,znum] = y0
          zcoord[ynum,znum] = z0
          x_shu = (r-0.5)*cos(theta) ;is this 0.5 Re the thickness of magnetopause?
          phi = atan(z0/y0)
          z_shu = sqrt((rp-0.5)^2/(1+tan(phi)^(-2))) ;;z_shu vs z0? Seem to use z0 for SH, but z_shu within MS
          y_shu = z_shu/tan(phi)

          if (abs(z0) eq 0.0) then begin
            z_shu = 0.0
            y_shu = (r-0.5)*sin(theta)
          endif
          if (abs(y0) eq 0.0) then begin
            y_shu = 0.0
            z_shu = (r-0.5)*sin(theta)
          endif


          ;also in Cooling+2001: Spreiter+66 density model

          rhosh = rho * (1.509 * exp(x_shu/rmp) + .1285)

          ; Force variables to be negative if they were origionally
          y_shu = abs(y_shu)*signy
          z_shu = abs(z_shu)*signz

          ; Cooling et al JGR 2001 model
          l = 3.*rmp/2.-x0
          Bmsx = -a*(-Bimfx*(1-Rmp/(2*l))+bimfy*(y0/l)+bimfz*(z0/l))
          Bmsy = a*(-Bimfx*(y0/(2*l))+bimfy*(2-y0^2/(l*rmp))-bimfz*(y0*z0/(l*rmp)))
          Bmsz= a*(-Bimfx*(z0/(2*l))-bimfy*(y0*z0/(l*rmp))+bimfz*(2-z0^2/(l*rmp)))

          ;btot = sqrt(Bmsx^2+Bmsy^2+Bmsz^2)*0.50
          ;ARROW, y0, z0, y0+Bmsy/btot, z0+Bmsz/btot, /DATA, COLOR=0, /SOLID, THICK=2.5, HSIZE=HSIZE

          ;geopack_ts04,param,x0,y0,z0,ext_0_x,ext_0_y,ext_0_z,epoch=epo
          geopack_t96,param,x_shu,y_shu,z_shu,ext_0_x,ext_0_y,ext_0_z,epoch=epo
          geopack_igrf_gsm,x_shu,y_shu,z_shu,igrf_0_x,igrf_0_y,igrf_0_z,epoch=epo
          bx = ext_0_x+igrf_0_x
          by = ext_0_y+igrf_0_y
          bz = ext_0_z+igrf_0_z
          ;btot = sqrt(bx^2+by^2+bz^2)*0.70
          ;ARROW, y_shu, z_shu, y_shu+by/btot, z_shu+bz/btot, /DATA, COLOR=0, /SOLID, THICK=2.5, HSIZE=HSIZE

          shear[ynum, znum] = getshear(bx,by,bz,Bmsx,Bmsy,Bmsz)
          rx_en[ynum, znum] = getrxben(bx,by,bz,Bmsx,Bmsy,Bmsz)
          va_cs[ynum, znum] = getvcs(bx,by,bz,Bmsx,Bmsy,Bmsz,0.1,rhosh)
          bisec[ynum, znum] = getbis(bx,by,bz,Bmsx,Bmsy,Bmsz)
          bshca[ynum, znum] = getca(Bmsx, Bmsy, Bmsz)
          bshmg[ynum, znum] = getmag(Bmsx, Bmsy, Bmsz)
          nsh[ynum, znum] = rhosh

          if sqrt(y_shu^2+z_shu^2) gt 31 then shear[ynum,znum] = 255
        endif
      endfor
    endfor
 endif
END
