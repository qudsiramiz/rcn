function getshear, bx0,by0,bz0,bx1,by1,bz1
  dp = (bx0*bx1+by0*by1+bz0*bz1)
  mag0 = sqrt(bx0^2+by0^2+bz0^2)
  mag1 = sqrt(bx1^2+by1^2+bz1^2)
  angle = acos(dp/(mag0*mag1))*180./!PI
  return, angle
end