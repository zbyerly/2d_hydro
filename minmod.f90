function minmod(a,b)
  implicit none
  double precision :: a,b,minmod

  minmod = 0.5d0*(sign(1d0,a)+sign(1d0,b))*min(abs(a),abs(b))

end function minmod
