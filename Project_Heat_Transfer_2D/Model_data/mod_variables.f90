module mod_variables
  use mod_param, only: n1, n2, n1p1, n2p1, n1p2, n2p2
  implicit none
  
  integer :: niter, nprof, iread
  real    :: len, rerror, pt
  real :: qx(n1p1, n2p2), qy(n1p2, n2p1), dx, dy, dx2, dy2, dt, re
  real :: qxold(n1p1, n2p2)
  real :: aaa(3), ggg(3), rrr(3)
  real :: hx(2:n1, 2:n2p1), hy(2:n1p1, 2:n2)
  real :: hxp(2:n1, 2:n2p1), hyp(2:n1p1, 2:n2)
  real :: ax(2:n1, 2:n2p1), ay(2:n1p1, 2:n2)
  real :: pres(n1, n2), phi(n1, n2), phip(n1, n2), presp(n1, n2)
  real :: x(n1p1), y(n2p1), xc(n1), yc(n2)
end module mod_variables
