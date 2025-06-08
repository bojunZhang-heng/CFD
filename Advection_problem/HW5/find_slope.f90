! ----------------------------------------------------------------------
! flag = 1 Downwind slope
! flag = 2 van Leer limiter
! flag = 3 SUPERBEE limiter
! Author: Bojun Zhang
! Date  : 2025-04-14 14:08:04
! File  : find_slope.f90
! ----------------------------------------------------------------------
subroutine find_slope(u_old, sigma_slope, dx, dt, CFL, N, flag)
  implicit none
  real*8,  intent(inout) :: u_old(0 : N+1), sigma_slope(0 : N+1)
  real*8,  intent(inout) :: dx, dt, CFL
  integer, intent(inout) :: N, flag

  real*8  :: a1, a2, a3, a4, s1, s2, s3
  integer :: jj

  real*8, allocatable :: sigma_L(:), sigma_R(:), sign1(:), sign2(:)
  
  allocate(sigma_L(0 : N), sigma_R(0 : N), sign1(0 : N), sign2(0 : N))
  do jj = 0, N
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
!      a1 = 2.d0 * (u_old(jj) - u_old(jj-1)) / dx
      a1 = (u_old(jj) - u_old(jj-1)) / dx
      a2 = (u_old(jj+1) - u_old(jj)) / dx
!      a3 = (u_old(jj) - u_old(jj-1)) / dx
!      a4 = (u_old(jj+1) - u_old(jj)) / dx
!      a4 = 2.d0 * (u_old(jj+1) - u_old(jj)) / dx
      s1 = sign(1.0d0, a1)
      s2 = sign(1.0d0, a2)
      
      if(s1 ==  s2) then
        sigma_L(jj) = s1 * dmin1(2.d0*a1, a2)
        sigma_R(jj) = s1 * dmin1(a1, 2.d0*a2)
      else
        sigma_L(jj) = 0.d0
        sigma_R(jj) = 0.d0
      endif

      sign1(jj) = sign(1.d0, sigma_L(jj))
      sign2(jj) = sign(1.d0, sigma_R(jj))
      if(sign1(jj) == sign2(jj)) then
        sigma_slope(jj) = sign1(jj) * dmax1(abs(sigma_L(jj)), abs(sigma_R(jj)))  
      else
        sigma_slope(jj) = 0.d0
      endif
    endif
!elseif(flag == 3) then
!  a1 = (u_old(jj) - u_old(jj-1)) / dx
!  a2 = (u_old(jj+1) - u_old(jj)) / dx
!
!  ! Compute sigma_L = minmod(2*a1, a2)
!  s1 = sign(1.0d0, 2.d0*a1)
!  s2 = sign(1.0d0, a2)
!  if (s1 == s2) then
!    sigma_L(jj) = s1 * min(abs(2.d0*a1), abs(a2))
!  else
!    sigma_L(jj) = 0.d0
!  endif
!
!  ! Compute sigma_R = minmod(a1, 2*a2)
!  s1 = sign(1.0d0, a1)
!  s2 = sign(1.0d0, 2.d0*a2)
!  if (s1 == s2) then
!    sigma_R(jj) = s1 * min(abs(a1), abs(2.d0*a2))
!  else
!    sigma_R(jj) = 0.d0
!  endif
!
!  sign1(jj) = sign(1.0d0, sigma_L(jj))
!  sign2(jj) = sign(1.0d0, sigma_R(jj))
!
!  if(sign1(jj) == sign2(jj)) then
!    sigma_slope(jj) = sign1(jj) * max(abs(sigma_L(jj)), abs(sigma_R(jj)))
!  else
!    sigma_slope(jj) = 0.d0
!  endif
!    endif
  enddo
  deallocate(sigma_L, sigma_R)
end subroutine find_slope


