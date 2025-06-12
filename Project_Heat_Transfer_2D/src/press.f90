subroutine press(irk)

  use mod_param
  use mod_variables
!
!      re = sqrt(pr/re)     
! LP Wang corrected the above typo of DeSpirito's code
      dt2re = aaa(irk) * dt/(2.0*re)
!
! internal...
      do jj = 2, n2m1
        do ii = 2, n1m1
          pres(ii, jj) = phi(ii, jj) - dt2re * &
                         ((phi(ii-1, jj) - 2.0*phi(ii, jj) + phi(ii+1, jj)) / dx2 + &
                         (phi(ii, jj-1)  - 2.0*phi(ii, jj) + phi(ii, jj+1)) / dy2)
          pres(ii, jj) = presp(ii, jj) + pres(ii, jj)
        enddo
      enddo
! sides...
      do jj = 2, n2m1
        pres(1, jj) = phi(1, jj) - dt2re  &
                    * ((phi(2, jj)  - phi(1, jj)) / dx2  &
                    + (phi(1, jj-1) - 2.0*phi(1, jj) + phi(1, jj+1)) / dy2)

        pres(1, jj) = presp(1, jj) + pres(1, jj)

        pres(n1, jj) = phi(n1, jj) - dt2re &
                     * ((phi(n1m1, jj)-phi(n1, jj)) / dx2  &
                     + (phi(n1, jj-1) - 2.0*phi(n1, jj) + phi(n1, jj+1)) / dy2)

        pres(n1, jj) = presp(n1, jj) + pres(n1, jj)
      enddo
      do  ii = 2, n1m1
        pres(ii, 1) = phi(ii, 1) - dt2re  &
                    * ((phi(ii-1, 1) - 2.0 * phi(ii, 1) + phi(ii+1, 1)) / dx2 &
                    + (phi(ii, 2) - phi(ii, 1)) / dy2)

        pres(ii, 1)  = presp(ii, 1) + pres(ii, 1)

        pres(ii, n2) = phi(ii, n2) - dt2re &
                     * ((phi(ii-1, n2) - 2.0*phi(ii, n2) + phi(ii+1, n2)) / dx2 &
                     + (phi(ii, n2m1)  - phi(ii, n2)) / dy2)

        pres(ii,n2)  = presp(ii, n2) + pres(ii, n2)
      enddo
! corners
        pres(1,1) = phi(1, 1) - dt2re &
                  * ((phi(2, 1) - phi(1, 1)) / dx2  &
                  + (phi(1, 2) -  phi(1, 1)) / dy2)

        pres(1,1)  = presp(1, 1) + pres(1, 1)
        
        pres(n1,1) = phi(n1, 1) - dt2re &
                   * ((phi(n1m1, 1) - phi(n1, 1)) / dx2  &
                   + (phi(n1, 2) - phi(n1, 1)) / dy2)

        pres(n1, 1) = presp(n1, 1) + pres(n1, 1)

        pres(1, n2) = phi(1, n2) - dt2re &
                    * ((phi(2, n2)  - phi(1, n2)) / dx2 &
                    + (phi(1, n2m1) - phi(1, n2)) / dy2)

        pres(1, n2)  = presp(1, n2) + pres(1, n2)
        pres(n1, n2) = phi(n1, n2)  - dt2re &
                     * ((phi(n1m1, n2) - phi(n1, n2)) / dx2  &
                     + (phi(n1, n2m1)  - phi(n1, n2)) / dy2)
        pres(n1, n2) = presp(n1, n2) + pres(n1,n2)
end subroutine
