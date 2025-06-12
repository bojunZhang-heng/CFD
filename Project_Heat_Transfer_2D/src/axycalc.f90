subroutine axycalc
  use mod_param
  use mod_variables

!      re = sqrt(pr/re)
      cc  = 1.0 / (2.0*re)
      ccx = cc / dx2
      ccy = cc / dy2

      do jj = 2, n2p1
        do ii = 2, n1
          ax1 = qx(ii+1, jj) + qx(ii-1, jj) - 2.0*qx(ii, jj)
          ax2 = qx(ii, jj+1) + qx(ii, jj-1) - 2.0*qx(ii, jj)
          ax(ii, jj) = ax1*ccx + ax2*ccy
        enddo
      enddo
      do jj = 2, n2
        do ii = 2, n1p1
          ay1 = qy(ii+1, jj) + qy(ii-1, jj) - 2.0*qy(ii, jj)
          ay2 = qy(ii, jj+1) + qy(ii, jj-1) - 2.0*qy(ii, jj)
          ay(ii, jj) = ay1*ccx + ay2*ccy
        enddo
      enddo

end subroutine 
