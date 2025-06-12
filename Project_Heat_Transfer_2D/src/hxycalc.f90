subroutine hxycalc
  use mod_param
  use mod_variables

      do jj = 2, n2p1
        do ii = 2, n1
!-----------------------------------------------------------        
! u and v are not in the same grid
!
!          temp1 = qx(ii, jj) * (qx(ii+1, jj) - qx(ii-1, jj))
!          temp2 = qy(ii, jj) * (qx(ii, jj+1) + qx(ii, jj-1))
!          hx(ii, jj) = -temp1/(2.d0*dx) - temp2/(2.d0*dy)

!------------------------------------------------------------       
! Another interpolation scheme I guess
!
!          temp1 = qx(ii, jj) * (qx(ii+1, jj) - qx(ii-1, jj))
!          tmp1  = (qy(ii, jj)   + qy(ii+1, jj))   / 2.d0
!          tmp2  = (qy(ii, jj-1) + qx(ii+1, jj-1)) / 2.d0
!          tmp3  = (tmp1 + tmp2) / 2.d0
!          temp2 = tmp3 * (qx(ii, jj+1) + qx(ii, jj-1))
!          hx(ii, jj) = -temp1/(2.d0*dx) - temp2/(2.d0*dy)

!------------------------------------------------------------
! Benchmark
!
          temp1 = qx(ii+1, jj)**2.0 - qx(ii-1, jj)**2.0
          temp2 = (qx(ii, jj+1) + qx(ii, jj)) * (qy(ii+1, jj) + qy(ii, jj)) &
                - (qx(ii, jj) + qx(ii, jj-1)) * (qy(ii+1, jj-1) + qy(ii,jj-1))
          hx(ii,jj) = -temp1 / (2.0*dx) - temp2 / (4.0*dy)

        enddo
      enddo

      do jj = 2, n2
        do ii = 2, n1p1
!-----------------------------------------------------------        
! Benchmark
          temp1 = (qx(ii, jj+1) + qx(ii, jj)) * (qy(ii+1, jj) + qy(ii, jj)) &
                - (qx(ii-1, jj+1) + qx(ii-1, jj)) * (qy(ii-1, jj) + qy(ii, jj))
          temp2 = qy(ii,jj+1)**2.0 - qy(ii,jj-1)**2.0
          hy(ii,jj) = -temp1 / (4.0*dx) - temp2 / (2.0*dy)

!-----------------------------------------------------------        
! u and v are not in the same grid
!
!          temp1 = qy(ii, jj) * (qy(ii, jj+1) - qy(ii, jj-1))
!          temp2 = qx(ii, jj) * (qy(ii+1, jj) + qy(ii-1, jj))
!          hy(ii, jj) = -temp1/(2.d0*dy) - temp2/(2.d0*dx)

!------------------------------------------------------------       
! Another interpolation scheme I guess
!          temp1 = qy(ii, jj) * (qy(ii, jj+1) - qy(ii, jj-1))
!          tmp1  = (qx(ii-1, jj)   + qx(ii, jj))   / 2.d0
!          tmp2  = (qy(ii-1, jj+1) + qx(ii, jj+1)) / 2.d0
!          tmp3  = (tmp1 + tmp2) / 2.d0
!          temp2 = tmp3 * (qy(ii, jj+1) + qy(ii, jj-1))
!          hx(ii, jj) = -temp1/(2.d0*dy) - temp2/(2.d0*dx)
        enddo
      enddo
end subroutine 
