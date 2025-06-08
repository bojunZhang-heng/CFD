!-----------------------------------------------------------------------
! PDE: u,t = alpha * u
! Algebraic equ. u_n+1 = u_n + alpha * u_n * dt
! The error confused me
!
!-----------------------------------------------------------------------
Program main 
implicit none

  integer, parameter         ::  N = 100
  real*8,  dimension(1:N+1)  ::  u, u_exact, actual_error, estimated_error
  real*8,  parameter         ::  dt = 0.02d0
  real*8,  parameter         ::  alpha = 0.1d0

  integer                    ::  ii
  real*8                     ::  tt
  u = 1.d0
  u_exact = 1.d0
  open(10, file='result.txt', status='unknown')  
  do ii = 1, N
    tt       = ii*dt
    u_exact(ii+1) = u(1) * exp(alpha*tt)
    u(ii+1)       = u(ii) + alpha*u(ii)*dt
    actual_error(ii)    = abs(u_exact(ii) - u(ii))
    estimated_error(ii) = 0.5d0 * dt * alpha**2 * u_exact(ii)
  enddo
  do ii = 1, N
    write(*,*) ii, u_exact(ii), u(ii), actual_error(ii), estimated_error(ii) 
    write(10,*) ii, u_exact(ii), u(ii), actual_error(ii), estimated_error(ii) 
  enddo

endProgram
