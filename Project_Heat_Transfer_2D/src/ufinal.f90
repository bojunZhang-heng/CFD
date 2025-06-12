subroutine ufinal(irk)
  use mod_param
  use mod_variables
  use mod_fft_data

! calculate the vel. at t=n+1...
      do jj = 2, n2p1
        do ii = 2, n1
          qx(ii, jj) = qx(ii, jj) - aaa(irk)*(phi(ii, jj-1) - phi(ii-1, jj-1))*dt / dx
        enddo
      enddo
!
      do jj = 2, n2
        do ii = 2, n1p1
          qy(ii, jj) = qy(ii, jj) - aaa(irk)*(phi(ii-1, jj)-phi(ii-1, jj-1)) * dt / dy
        enddo
      enddo
!
! fill bc for velocities defined outside the domain...
      do ii = 2, n1
        qx(ii,1)     = -qx(ii,2)
        qx(ii, n2p2) = 2.0 - qx(ii, n2p1)
      enddo
      do jj = 2, n2
        qy(1, jj)    = -qy(2, jj)
        qy(n1p2, jj) = -qy(n1p1, jj)
      enddo
!
end subroutine
