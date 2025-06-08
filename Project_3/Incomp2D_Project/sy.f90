subroutine sy(il, iu, bb, dd, aa, cc)
!*********************************************************************
! Subroutine SY, for solving a tridiagonal system of equations.  Based
! on the Thomas algorithm, from "Computational Fluid Mechanics and
! Heat Transfer," by Anderson, Tannehill, and Pletcher.
! 27 May 1992   Jim DeSpirito
!
  integer, intent(in) :: il, iu
  real, intent(inout) :: bb(il:iu), dd(il:iu), aa(il:iu), cc(il:iu)
!
!      double precision aa(iu),bb(iu),cc(iu),dd(iu)
!
!.....il = subscript of first equation
!.....iu = subscript of last equation
!.....bb = coefficient behind diagonal
!.....dd = coefficient on diagonal
!.....aa = coefficient ahead of diagonal
!.....cc = element of constant vector
!
!.....establish upper triangular matrix
!
      lp = il + 1
      do 10 ii = lp,iu
       r = bb(ii) / dd(ii-1)
       dd(ii) = dd(ii) - r * aa(ii-1)
       cc(ii) = cc(ii) - r * cc(ii-1)
   10 continue
!
!.....back substitution
!
      cc(iu) = cc(iu) / dd(iu)
      do 20 ii = lp,iu
       jj = iu - ii + il
       cc(jj) = (cc(jj) - aa(jj) * cc(jj+1)) / dd(jj)
   20 continue
!
!.....solution stored in cc
!

end subroutine
