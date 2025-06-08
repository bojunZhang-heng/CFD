! ----------------------------------------------------------------------
! flag = 1 Downwind slope
! flag = 2 van Leer limiter
! flag = 3 SUPERBEE limiter
! Author: Bojun Zhang
! Date  : 2025-04-14 14:08:04
! File  : find_slope.f90
! ----------------------------------------------------------------------
subroutine find_slope(u_old, sigma_slope, dx, N, flag)
  implicit none
  integer                :: jj
  integer, intent(in)    :: N, flag

  real*8,  intent(inout) :: u_old(0:N+1), sigma_slope(N)
  real*8,  intent(in)    :: dx
  real*8                 :: a1, a2, a3, a4, s1, s2, s3
  real*8, allocatable    :: sigma_L(:), sigma_R(:), sign1(:), sign2(:)

  
  allocate(sigma_L(N), sigma_R(N), sign1(N), sign2(N))
  do jj = 1, N
    if(flag == 1) then
      sigma_slope(jj) = (u_old(jj+1) - u_old(jj)) / dx
    elseif(flag == 2) then
      a1 = 2.d0 * (u_old(jj+1) - u_old(jj)) / dx
      a2 = (u_old(jj+1) - u_old(jj-1)) / (2.d0*dx)
      a3 = 2.d0 * (u_old(jj) - u_old(jj-1)) / dx
  
      s1 = sign(1.0d0, a1)
      s2 = sign(1.0d0, a2)
      s3 = sign(1.0d0, a3)
  
      if(s1 == s2 .and. s2 == s3) then
        sigma_slope(jj) = s1 * dmin1(abs(a1), abs(a2), abs(a3))
      else
        sigma_slope(jj) = 0.d0 
      endif
    elseif(flag == 3) then
      a1 = (u_old(jj) - u_old(jj-1)) / dx
      a2 = (u_old(jj+1) - u_old(jj)) / dx
      !write(*,*) a1

      s1 = sign(1.0d0, a1)
      s2 = sign(1.0d0, a2)
      !write(*,*) s1

      if (s1 * s2 .ge. 0.0d0) then
        sigma_L(jj) = s1 * dmin1(abs(2.d0*a1), abs(a2))
        sigma_R(jj) = s1 * dmin1(abs(a1), abs(2.d0*a2))
      else
        sigma_L(jj) = 0.d0
        sigma_R(jj) = 0.d0
      endif

      sign1(jj) = sign(1.d0, sigma_L(jj))
      sign2(jj) = sign(1.d0, sigma_R(jj))
      if (sign1(jj) * sign2(jj) .ge. 0.d0) then
        sigma_slope(jj) = sign1(jj) * dmax1(abs(sigma_L(jj)), abs(sigma_R(jj)))  
      else
        sigma_slope(jj) = 0.d0
      endif
    endif
  enddo
  deallocate(sigma_L, sigma_R)
  deallocate(sign1, sign2)
end subroutine find_slope
!----------------------------------------------
! a1          the slope var
! a2          the slope var
! s1          the sign of a1
! s2          the sign of a2
!
!----------------------------------------------

