subroutine qhat(irk)

  use mod_param
  use mod_variables
!note: max dim. on a,b,c,d is max(n1p1,n2p1)...
      dimension a(n1p1), b(n1p1), c(n1p1), d(n1p1)
!      re = sqrt(pr/re)
!
! the qx'
!
      do jj = 2, n2p1
        do ii = 2, n1
          ax(ii, jj) = dt * (ggg(irk) * hx(ii, jj) + rrr(irk) * hxp(ii, jj) &
                     - aaa(irk) * (pres(ii, jj-1)-pres(ii-1, jj-1)) / dx) &
                     + 2.0*aaa(irk) * dt * ax(ii, jj)
! Note: in this current ax in LHS is b !!! in the AX=b
!     :                 ax in RHS is just viscous term calculate by qx and qy
        enddo
      enddo
!
! set up vectors for tri-diag calculation
!
      aa = -dt * aaa(irk) / (2.0 * re*dx2)
      bb = 1.0 + dt * aaa(irk) / (re*dx2)
      do jj = 2, n2p1
        do ii= 2, n1
! note: a,b,c are constant here - but leave for generality
          a(ii) = aa
          b(ii) = bb
          c(ii) = a(ii)
          d(ii) = ax(ii, jj)
        enddo

        call sy(2,n1,a,b,c,d)
        do ii = 2, n1
          ax(ii, jj) = d(ii)
! Note: in this current d(i) is the solution for u^ - u^n
        enddo
      enddo

      aa = -dt * aaa(irk) / (2.0*re*dy2)
      bb = 1.0 + dt * aaa(irk) / (re*dy2)
      do ii = 2, n1
        do jj = 2, n2p1
          a(jj) = aa
          b(jj) = bb
          c(jj) = a(jj)
          d(jj) = ax(ii, jj)
        enddo

! LPW2017: boundary conditions
          b(2) = b(2) - a(2)
          b(n2p1) = b(n2p1) - c(n2p1)

        call sy(2,n2p1,a,b,c,d)
        do jj = 2, n2p1
          ax(ii, jj) = d(jj)
          qx(ii, jj) = ax(ii, jj) + qx(ii, jj)
! Note: in this current, qx is u^, the intermediate velocity
        enddo
      enddo
!
! the qy's
!
        do jj = 2, n2
          do ii = 2, n1p1
            ay(ii, jj) = dt * (ggg(irk) * hy(ii, jj) + rrr(irk) * hyp(ii, jj) &
                       - aaa(irk) * (pres(ii-1, jj) - pres(ii-1, jj-1)) / dy) &
                       + 2.0 * aaa(irk) * dt * ay(ii, jj)
          enddo
        enddo
!
! set up the vectors for tri-diag calculation
      aa = -dt * aaa(irk) / (2.0*re*dx2)
      bb = 1.0 + dt * aaa(irk) / (re*dx2)
      do jj = 2, n2
        do ii = 2, n1p1
          a(ii) = aa
          b(ii) = bb
          c(ii) = a(ii)
          d(ii) = ay(ii, jj)
        enddo

! LPW2017: boundary conditions
          b(2) = b(2) - a(2)
          b(n1p1) = b(n1p1) - c(n1p1)

        call sy(2,n1p1,a,b,c,d)
        do  ii = 2, n1p1
          ay(ii, jj) = d(ii)
        enddo
      enddo

      aa = -dt * aaa(irk) / (2.0*re*dy2)
      bb = 1.0 + dt*aaa(irk)/(re*dy2)
      do ii = 2, n1p1
        do jj = 2, n2
          a(jj) = aa
          b(jj) = bb
          c(jj) = a(jj)
          d(jj) = ay(ii, jj)
        enddo

        call sy(2,n2,a,b,c,d)
        do jj = 2, n2
          ay(ii, jj) = d(jj)
          qy(ii, jj) = ay(ii, jj) + qy(ii, jj)
        enddo
      enddo
!
! save previous hx,hy values in hxp,hyp
!
      do jj = 2, n2p1
        do ii = 2, n1
          hxp(ii, jj) = hx(ii, jj)
        enddo
      enddo
      do jj = 2, n2
        do ii = 2, n1p1
          hyp(ii, jj) = hy(ii, jj)
        enddo
      enddo

end subroutine 
