subroutine trigon(irk)
! use trigonometric series method to calculate the poisson equation
  use mod_fft_data
  use mod_param
  use mod_variables
  use fftcos

      dimension a(n2), b(n2), c(n2), d(n2), qrhs(n1,n2)
!
! save old pressure value first (for pressure calc only)
!
        do jj = 1, n2
          do ii = 1, n1
            presp(ii, jj) = pres(ii, jj)
          enddo
        enddo
!
      do jj = 1, n2
        do ii = 1, n1
          qrhs(ii, jj) = ((qx(ii+1, jj+1) - qx(ii, jj+1)) / dx &
                       + (qy(ii+1, jj+1)  - qy(ii+1, jj)) / dy) / (aaa(irk) * dt)
          qrhs(ii, jj) = dx2*qrhs(ii, jj)
        enddo
      enddo
!
! transform the rhs in x-direction
      pi = 6.*asin(0.5)
      CALL R2FC1(qrhs)

      phi(1,1) = 0.0
      phi(1,2) = qrhs(1,1) / (dx2/dy2)
      do jj = 3, n2
        phi(1, jj) = -phi(1, jj-2) + 2.0 * phi(1, jj-1) &
                   + qrhs(1, jj-1) / (dx2/dy2)
      enddo

      do ii = 2, n1

        waven = 2.*(1.0 - cos(pi*(ii-1) / n1))
        do jj = 1, n2
          a(jj) = dx2/dy2
          b(jj) = -2.0* dx2 / dy2 - waven
          c(jj) = a(jj)
          d(jj) = qrhs(ii, jj)
        enddo

        b(1)  = b(1)  + c(1)
        b(n2) = b(n2) + c(n2)

        call sy(1,n2,a,b,c,d)
        do jj = 1, n2
          phi(ii, jj) = d(jj)
        enddo

      enddo
!
! transform the solution back to physical space
      CALL F2RC1(phi)

end subroutine      
