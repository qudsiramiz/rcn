pro ridgefinder, image=image,first_deri=first_deri,second_deri=second_deri,ridge=ridge,tan_clockangle=tan_clockangle

;written by Ying, modified from CV_Vessel_Filter_2D
;only image is input, all the rest are output

 image=FLOAT(image)
 
 ; Smooth image, kernel can be changed
    ; a larger kernel will detect bigger vessels
    smooth_kernel=[[1,3,3,1],$
               [3,9,9,3],$
               [3,9,9,3],$
               [1,3,3,1]]

    smooth_image=CONVOL(image, smooth_kernel, 64)
    image=smooth_image

 ridge=image
 ridge[*]=0
     x_size=(SIZE(image, /DIMENSIONS))[0]
     y_size=(SIZE(image, /DIMENSIONS))[1]

    ; Calculate the Hessian Matrix
     x_deriv_kernel=[-1,0,1]
     y_deriv_kernel=TRANSPOSE([-1,0,1])

     Lx=CONVOL(image, x_deriv_kernel)
     Ly=CONVOL(image, y_deriv_kernel)

     Lxy=CONVOL(Lx, y_deriv_kernel)
     Lxx=CONVOL(Lx, x_deriv_kernel)
     Lyy=CONVOL(Ly, y_deriv_kernel)

     first_deri=Lxy
     second_deri=Lxy
     for i=0, x_size-1 do begin
        for j=0, y_size-1 do begin
            result=eigenql([[Lxx[i,j],Lxy[i,j]],[Lxy[i,j],Lyy[i,j]]],/ascending,eigenvectors=eigenv);,/absolute
            second_deri[i,j]=result[0]
            ;help,eigenv
            step=1
            adjac1=bilinear(image,i+eigenv[0,0]*step,j+eigenv[0,1]*step)
            adjac2=bilinear(image,i-eigenv[0,0]*step,j-eigenv[0,1]*step)
            if (adjac1-image[i,j])*(image[i,j]-adjac2) lt 0 and result[0] lt 0 then ridge[i,j]=(second_deri[i,j])^2
        endfor
     endfor
     
     
     ;for i=1, x_size-2 do begin
        ;for j=2, y_size-2 do begin
            ;if second_deri[i,j] lt 0 then begin
                ;print,i,j,second_deri[i,j],first_deri[i,j],first_deri[i-1,j],first_deri[i,j-1]
                ;if first_deri[i+1,j]*first_deri[i,j+1] lt 0 or first_deri[i+1,j]*first_deri[i,j-1] lt 0 or first_deri[i-1,j]*first_deri[i,j-1] lt 0 or ;first_deri[i-1,j]*first_deri[i,j+1] lt 0 then begin
                    ;ridge[i,j]=(second_deri[i,j])^2
                ;endif
                ;;if first_deri[i,j] lt 5 then ridge[i,j]=(second_deri[i,j])^2
            ;endif
        ;endfor
     ;endfor
     
 end
     
