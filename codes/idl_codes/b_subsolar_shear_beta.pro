@'C:\Users\ahmad\IDLWorkspace\Default\coyote\cgcolorbar.pro'

function getshear, bx0,by0,bz0,bx1,by1,bz1
  dp = (bx0*bx1+by0*by1+bz0*bz1)
  mag0 = sqrt(bx0^2+by0^2+bz0^2)
  mag1 = sqrt(bx1^2+by1^2+bz1^2)
  angle = acos(dp/(mag0*mag1))*180./!PI
  return, angle
end

; Brian Walsh
; 15 Feb 2012
; Made to show shear angle between draped imf from the cooling et al 2001 model
; and the internal geomagnetic field from the sun's point of view.

PRO B_subsolar_shear_beta

    thm_init
    ;you may want to control the thickness of the generated plots
    ;to make the lines in the output more visible when exported
    ;you can do so with these variables
    axisthick = 4.0
    charthick = 6.0
    thick = 10.0
    charsize = 1.0
    symsize = 1.0
    !X.STYLE=1
    !Y.STYLE=1
    !X.MARGIN = 7
    HSIZE = 100
    popen,'Z:\Shared drives\Casella\SpacePlasFest\mms\rcn\figures\shear_brian_v4',encapsulated=1,land=0,xsize=6,ysize=6
    loadct,39
    
    ;;;TS05 magnetic field model SW values 

    Bimfx = -2.27
    Bimfy = -5.79
    Bimfz = 0.37   
   
    BY = Bimfy
    BZ = Bimfz
    VEL=638.
    N=3.4
    mp = 1.67262158e-27 ; mass of proton in kg
    DST=-18
    rho = 1.94d-6 * n * vel^2.    
    pdyn = rho
    
    geopack_recalc,2016,12,24,15,08,0,/date
    geopack_epoch,epo,2016,12,24,15,08,0,/compute_epoch
    param = [rho,dst,by,bz,0,0,0,0,0,0]
    geopack_getw,n,vel,bz,ww
 
    fmt = '(f6.2)'
    title = 'Bimf [x,y,z]: '+string(Bimfx,format=fmt)+', '+string(Bimfy,format=fmt)+', '+string(Bimfz,format=fmt)
 ;   plot, [0,0], [1,1], xrange = [-35,35], yrange=[-35, 35], /isotropic, xtitle = 'Y [Re]', ytitle='Z [Re]', title=title
    v=150
    shear = fltarr(v,v)
    ycoord = MAKE_ARRAY(v,v, /INTEGER, VALUE = -9999)
    zcoord = MAKE_ARRAY(v,v, /INTEGER, VALUE = -9999)
    
    xx_shu = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    yy_shu = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    zz_shu = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    
    bms_x = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    bms_y = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    bms_z = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)  

    b_x = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    b_y = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    b_z = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    
    bx_ext = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    by_ext = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    bz_ext = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    
    bx_igrf = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    by_igrf = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
    bz_igrf = MAKE_ARRAY(v,v, /FLOAT, VALUE = -9999)
  
    ro=(10.22+1.29*tanh(0.184*(bz+8.14)))*(Pdyn)^(-1.0/6.6)
    print, ro, bz
    alpha=(0.58-0.007*bz)*(1+0.024*ALOG(Pdyn))
    rmp = ro*(2/(1+COS(0.0)))^alpha  ; Stand off position of the magnetopause
    dtheta = !PI/100.
    
    ; At each distance down the tail from x=8.75 to x=-60 with 1 Re cadience get the B vector at the equator
    a=2.
    leny = 80
    lenz = 80
    for ynum = 0, leny-1 do begin
      y0 = 40.-float(ynum)
      ;print, y0
      for znum = 0, lenz-1 do begin    
        ; Find caps at bottom and top for each x value        
        z0 = 40.-float(znum)

        rp = sqrt(z0^2+y0^2) ; projection of r into the xy plane
        
        ; Find the theta value that corresponds with this point on the Shue et al 1998 magnetopause             
        for index = 0, 100 do begin
          theta = float(index)*dtheta
          r = ro*(2/(1+COS(theta)))^alpha

          zp = r*sin(theta)
          x0 = r*cos(theta)
          signx = x0/abs(x0)
          signy = y0/abs(y0)
          signz = z0/abs(z0)
          if x0 eq 0.0 then signx = 1.0
          if y0 eq 0.0 then signy = 1.0
          if z0 eq 0.0 then signz = 1.0          
  
          if rp le zp then begin      
            ;print, index, rp, zp
            GOTO, JUMP1
          endif
        endfor 
        
        JUMP1: if rp le zp then begin
          ycoord[ynum,znum] = y0       
          zcoord[ynum,znum] = z0                 
          ;print, ynum, znum, y0, z0, ycoord[ynum,znum], zcoord[ynum,znum]
          x_shu = (r-0.5)*cos(theta)
          print, ynum, znum, theta, x_shu
          xx_shu[ynum,znum] = x_shu
          phi = atan(z0/y0)
          z_shu = sqrt((rp-0.5)^2/(1+tan(phi)^(-2)))
          y_shu = z_shu/tan(phi)
  
          if (abs(z0) eq 0.0) then begin
            z_shu = 0.0
            y_shu = (r-0.5)*sin(theta)
          endif
          if (abs(y0) eq 0.0) then begin
            y_shu = 0.0
            z_shu = (r-0.5)*sin(theta)
          endif          
          ; Force variables to be negative if they were origionally
          y_shu = abs(y_shu)*signy
          z_shu = abs(z_shu)*signz
          yy_shu[ynum, znum] = y_shu
          zz_shu[ynum, znum] = z_shu
          ;print, 
          ; Cooling et al JGR 2001 model
          l = 3.*rmp/2.-x0
          Bmsx = -a*(-Bimfx*(1-Rmp/(2*l))+bimfy*(y0/l)+bimfz*(z0/l))
          Bmsy = a*(-Bimfx*(y0/(2*l))+bimfy*(2-y0^2/(l*rmp))-bimfz*(y0*z0/(l*rmp)))
          Bmsz= a*(-Bimfx*(z0/(2*l))-bimfy*(y0*z0/(l*rmp))+bimfz*(2-z0^2/(l*rmp)))
          
          bms_x[ynum,znum] = Bmsx
          bms_y[ynum,znum] = Bmsy
          bms_z[ynum,znum] = Bmsz
       ;   btot = sqrt(Bmsx^2+Bmsy^2+Bmsz^2)*0.50
      ;    ARROW, y0, z0, y0+Bmsy/btot, z0+Bmsz/btot, /DATA, COLOR=0, /SOLID, THICK=2.5, HSIZE=HSIZE

       ;   geopack_ts04,param,x0,y0,z0,ext_0_x,ext_0_y,ext_0_z,epoch=epo          
          geopack_t96,param,x_shu,y_shu,z_shu,ext_0_x,ext_0_y,ext_0_z,epoch=epo
          geopack_igrf_gsm,x_shu,y_shu,z_shu,igrf_0_x,igrf_0_y,igrf_0_z,epoch=epo
          bx = ext_0_x+igrf_0_x
          by = ext_0_y+igrf_0_y
          bz = ext_0_z+igrf_0_z
          
          b_x[ynum,znum] = Bx
          b_y[ynum,znum] = By
          b_z[ynum,znum] = Bz 
          
          bx_ext[ynum,znum] = ext_0_x
          by_ext[ynum,znum] = ext_0_y
          bz_ext[ynum,znum] = ext_0_z

          bx_igrf[ynum,znum] = igrf_0_x
          by_igrf[ynum,znum] = igrf_0_y
          bz_igrf[ynum,znum] = igrf_0_z

   ;       btot = sqrt(bx^2+by^2+bz^2)*0.70      
   ;       ARROW, y_shu, z_shu, y_shu+by/btot, z_shu+bz/btot, /DATA, COLOR=0, /SOLID, THICK=2.5, HSIZE=HSIZE          

          shear[ynum,znum] = getshear(bx,by,bz,Bmsx,Bmsy,Bmsz)
          if sqrt(y_shu^2+z_shu^2) gt 31 then shear[ynum,znum] = 255          
        endif
      endfor      
    endfor
;    stop
    ;B = WHERE_XYZ(ycoord GT -100, XIND=xind, YIND=yind)
    ;print, size(ycoord)
    ;WRITE_CSV, 'bms_x.csv', bms_x
    ;WRITE_CSV, 'bms_y.csv', bms_y
    ;WRITE_CSV, 'bms_z.csv', bms_z
    
    ;WRITE_CSV, 'b_x.csv', b_x
    ;WRITE_CSV, 'b_y.csv', b_y
    ;WRITE_CSV, 'b_z.csv', b_z 
;
;    file = 'Z:\Shared drives\Casella\SpacePlasFest\mms\rcn\data\all_data_rx_model_0.5re_t96_20211012_v01_idl.h5'
;    fid = h5f_create(file)
;
;    data_xxshu = xx_shu
;    xxshu_type_id = H5T_IDL_CREATE(data_xxshu)
;    xxshu_space_id = H5S_CREATE_SIMPLE(size(data_xxshu,/DIMENSIONS))
;    xxshu_set_id = H5D_CREATE(fid,'x_shu',xxshu_type_id,xxshu_space_id)
;        
;    data_yyshu = yy_shu
;    yyshu_type_id = H5T_IDL_CREATE(data_yyshu)
;    yyshu_space_id = H5S_CREATE_SIMPLE(size(data_yyshu,/DIMENSIONS))
;    yyshu_set_id = H5D_CREATE(fid,'y_shu',yyshu_type_id,yyshu_space_id)
;    
;    data_zzshu = zz_shu
;    zzshu_type_id = H5T_IDL_CREATE(data_zzshu)
;    zzshu_space_id = H5S_CREATE_SIMPLE(size(data_zzshu,/DIMENSIONS))
;    zzshu_set_id = H5D_CREATE(fid,'z_shu',zzshu_type_id,zzshu_space_id)
;    
;
;    data_shear= shear
;    shear_type_id = H5T_IDL_CREATE(data_shear)
;    shear_space_id = H5S_CREATE_SIMPLE(size(data_shear,/DIMENSIONS))
;    shear_set_id = H5D_CREATE(fid,'shear',shear_type_id,shear_space_id)
;
;    data_x = b_x
;    bx_type_id = H5T_IDL_CREATE(data_x)
;    bx_space_id = H5S_CREATE_SIMPLE(size(data_x,/DIMENSIONS))
;    bx_set_id = H5D_CREATE(fid,'bx',bx_type_id,bx_space_id)
;    
;    data_y = b_y
;    by_type_id = H5T_IDL_CREATE(data_y)
;    by_space_id = H5S_CREATE_SIMPLE(size(data_y,/DIMENSIONS))
;    by_set_id = H5D_CREATE(fid,'by',by_type_id,by_space_id)
;    
;    data_z = b_z
;    bz_type_id = H5T_IDL_CREATE(data_z)
;    bz_space_id = H5S_CREATE_SIMPLE(size(data_z,/DIMENSIONS))
;    bz_set_id = H5D_CREATE(fid,'bz',bz_type_id,bz_space_id)
;    
;    data_bxigrf = bx_igrf
;    bxigrf_type_id = H5T_IDL_CREATE(data_bxigrf)
;    bxigrf_space_id = H5S_CREATE_SIMPLE(size(data_bxigrf,/DIMENSIONS))
;    bxigrf_set_id = H5D_CREATE(fid,'bx_igrf',bxigrf_type_id,bxigrf_space_id)
;
;    data_byigrf = by_igrf
;    byigrf_type_id = H5T_IDL_CREATE(data_byigrf)
;    byigrf_space_id = H5S_CREATE_SIMPLE(size(data_byigrf,/DIMENSIONS))
;    byigrf_set_id = H5D_CREATE(fid,'by_igrf',byigrf_type_id,byigrf_space_id)
;    
;    data_bzigrf = bz_igrf
;    bzigrf_type_id = H5T_IDL_CREATE(data_bzigrf)
;    bzigrf_space_id = H5S_CREATE_SIMPLE(size(data_bzigrf,/DIMENSIONS))
;    bzigrf_set_id = H5D_CREATE(fid,'bz_igrf',bzigrf_type_id,bzigrf_space_id)
;
;    data_bx_ext = bx_ext
;    bx_ext_type_id = H5T_IDL_CREATE(data_bx_ext)
;    bx_ext_space_id = H5S_CREATE_SIMPLE(size(data_bx_ext,/DIMENSIONS))
;    bx_ext_set_id = H5D_CREATE(fid,'bx_ext',bx_ext_type_id,bx_ext_space_id)
;
;    data_by_ext = by_ext
;    by_ext_type_id = H5T_IDL_CREATE(data_by_ext)
;    by_ext_space_id = H5S_CREATE_SIMPLE(size(data_by_ext,/DIMENSIONS))
;    by_ext_set_id = H5D_CREATE(fid,'by_ext',by_ext_type_id,by_ext_space_id)
;    
;    data_bz_ext = bz_ext
;    bz_ext_type_id = H5T_IDL_CREATE(data_bz_ext)
;    bz_ext_space_id = H5S_CREATE_SIMPLE(size(data_bz_ext,/DIMENSIONS))
;    bz_ext_set_id = H5D_CREATE(fid,'bz_ext',bz_ext_type_id,bz_ext_space_id)
;
;    H5D_WRITE,xxshu_set_id,data_xxshu        
;    H5D_WRITE,yyshu_set_id,data_yyshu
;    H5D_WRITE,zzshu_set_id,data_zzshu
;    
;    H5D_WRITE,shear_set_id,data_shear
;    
;    H5D_WRITE,bx_set_id,data_x
;    H5D_WRITE,by_set_id,data_y
;    H5D_WRITE,bz_set_id,data_z
;    
;    H5D_WRITE,bx_ext_set_id,data_bx_ext
;    H5D_WRITE,by_ext_set_id,data_by_ext
;    H5D_WRITE,bz_ext_set_id,data_bz_ext
;    
;    H5D_WRITE,bxigrf_set_id,data_bxigrf
;    H5D_WRITE,byigrf_set_id,data_byigrf
;    H5D_WRITE,bzigrf_set_id,data_bzigrf
;    ;; close all open identifiers
;    H5D_CLOSE,yyshu_set_id
;    H5S_CLOSE,yyshu_space_id
;    H5T_CLOSE,yyshu_type_id
;    
;    H5D_CLOSE,zzshu_set_id
;    H5S_CLOSE,zzshu_space_id
;    H5T_CLOSE,zzshu_type_id
;    
;    H5D_CLOSE,shear_set_id
;    H5S_CLOSE,shear_space_id
;    H5T_CLOSE,shear_type_id
;    
;    H5D_CLOSE,bx_set_id
;    H5S_CLOSE,bx_space_id
;    H5T_CLOSE,bx_type_id
;    
;    H5D_CLOSE,by_set_id
;    H5S_CLOSE,by_space_id
;    H5T_CLOSE,by_type_id
;    
;    H5D_CLOSE,bz_set_id
;    H5S_CLOSE,bz_space_id
;    H5T_CLOSE,bz_type_id
;    
;    H5D_CLOSE,bxigrf_set_id
;    H5S_CLOSE,bxigrf_space_id
;    H5T_CLOSE,bxigrf_type_id
;    
;    H5D_CLOSE,byigrf_set_id
;    H5S_CLOSE,byigrf_space_id
;    H5T_CLOSE,byigrf_type_id
;    
;    H5D_CLOSE,bzigrf_set_id
;    H5S_CLOSE,bzigrf_space_id
;    H5T_CLOSE,bzigrf_type_id
;    
;    H5D_CLOSE,bx_ext_set_id
;    H5S_CLOSE,bx_ext_space_id
;    H5T_CLOSE,bx_ext_type_id
;    
;    H5D_CLOSE,by_ext_set_id
;    H5S_CLOSE,by_ext_space_id
;    H5T_CLOSE,by_ext_type_id
;
;    H5D_CLOSE,bz_ext_set_id
;    H5S_CLOSE,bz_ext_space_id
;    H5T_CLOSE,bz_ext_type_id
;
;    H5F_CLOSE,fid
    contour, shear, ycoord,zcoord,$
      /CELL_FILL,MIN_VALUE=0.0,MAX_VALUE=180.0,NLEVELS=400,xrange = [-15,15], yrange=[-15, 15],$
      /isotropic, xtitle = 'Y [Re]', ytitle='Z [Re]', title=title;, /IRREGULAR
    
    contour, shear, ycoord,zcoord,levels=[155],$
      xrange = [-15,15], yrange=[-15, 15],thick = 4,c_linestyle=1,$
      /isotropic, xtitle = 'Y [Re]', ytitle='Z [Re]', title=title, /overplot
    
    
;    contour, shear[0:leny-1,0:lenz-1], ycoord[0:leny-1,0:lenz-1],zcoord[0:leny-1,0:lenz-1],$
;      /CELL_FILL,MIN_VALUE=0.0,MAX_VALUE=180.0,NLEVELS=400,xrange = [-35,35], yrange=[-35, 35],$
;       /isotropic, xtitle = 'Y [Re]', ytitle='Z [Re]', title=title, /IRREGULAR 

    cgCOLORBAR, NColors=!D.N_colors-5, MINRANGE=0, MAXRANGE=180, Divisions=4, $
        Format='(I6)', Color=0,RIGHT=1,VERTICAL=1,MINOR=9, TITLE='Shear Angle [Degrees]',$
        Position=[0.91, 0.10, 0.93, 0.85],charsize = 0.8

    ;Plot Terminator
    radius = ro*(2/(1+COS(!PI/2.)))^alpha 
    points = (2*!PI / 99.0) * FINDGEN(100)
    x = radius * COS(points)
    y = radius * SIN(points)
    oplot, y,x, thick = 7
    
    ;Plot End of magnetopause
    radius = 31.0 
    points = (2*!PI / 99.0) * FINDGEN(100)
    x = radius * COS(points)
    y = radius * SIN(points)
    oplot, y,x, thick = 7    

;      myplt.save, 'Z:\Shared drives\Casella\SpacePlasFest\mms\rcn\figures\shear_brian.pdf', BORDER=10, RESOLUTION=300, /TRANSPARENT
    pclose
    print, 'Done Plotting'
END
