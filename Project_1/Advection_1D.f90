!-----------------------------------------------------
! Use Lax-Wendroff Scheme
!
!
!
!-----------------------------------------------------
Program main
  implicit none
  integer :: N, jj, N_step, it, nn, mode
  integer :: flag
  real*8  :: dx, dt, L, tt, b2, b3, b4, a_n, k_n, pp, tolerance
  real*8  :: x_min, x_max, x_left, x_right
  real*8  :: t_min, t_max
  real*8  :: CFL, aa, pi, sigma_n, beta_n, A, B
  real*8  :: L1_error_u_N, L1_error_u_M, L1_error_u_AS
  real*8  :: N_L1
  real*8, allocatable :: u_all(:,:)
  real*8, allocatable :: sigma_slope(:), u_NSL(:), u_oldSL(:)
  real*8, allocatable :: u(:), x_coor(:), u_old(:)
  real*8, allocatable :: u_N(:), u_M(:), u_AS(:)

! Preporcess------------------------------------------------------------------------------
!-----------------------------------------------------------
! Declare inital value
  x_min = -10.0
  x_max =  10.0
  t_min =  0.0
  t_max =  10.0
  L   = x_max - x_min
  CFL = 0.5d0           ! CFL = a * dt / dx
!  CFL = 1.d0           ! CFL = a * dt / dx
  N   = 2499
  pp  = 2.d0
  aa  = 1.d0
  dx  = L / real(N+1)
  dt  = CFL * dx / aa
  pi  = 4.d0*atan(1.d0)
  mode   = 10000
  N_step = int((t_max - t_min) / dt)
  b2  = 0.d0
  b3  = -aa * dx**2.d0 / 6 * (1 - CFL**2.d0)
  b4  = -aa * dx**3.d0 / 8 * (1 - CFL**2.d0) * CFL
  a_n = 0.d0
  k_n = 0.d0
  sigma_n   = 0.d0
  beta_n    = 0.d0
  tolerance = 1.0d-8
  L1_error_u_N  = 0.d0
  L1_error_u_M  = 0.d0
  L1_error_u_AS = 0.d0
  open(20, file='solution.txt',  status='unknown')
  open(30, file='L1_error_u_N.txt',  status='unknown', position='append')
  open(40, file='L1_error_u_M.txt',  status='unknown', position='append')
  open(50, file='L1_error_u_AS.txt', status='unknown', position='append')
  open(60, file='solution_SL.txt', status='unknown')

!-----------------------------------------------------------
  open(10, file='output.dat', status = 'unknown')

allocate(u_all(0:N+1, 0:N_step))
allocate(u_old(0:N+1), u_N(0:N+1), u(0:N+1))
allocate(x_coor(0:N+1), u_M(0:N+1), u_AS(0:N+1))
allocate(sigma_slope(N), u_NSL(0:N+1))
allocate(u_oldSL(0:N+1))

! Driver--------------------------------------------------------------------------------
!-----------------------------------------------------------
  u_N     = 0.d0
  u_old   = 0.d0
  u_oldSL = 0.d0
  u_M     = 0.d0
  u_AS    = 0.d0
  u       = 0.d0

!-----------------------------------------------------------
! Intial B.C 
! tt = 0
! Physical coordinate exclude x(0) and x(N+1)
  do jj = 0, N
    x_coor(jj) = x_min + jj * dx
  enddo
  
  x_left  = -1/2.d0
  x_right =  1/2.d0
  do jj = 1, N
    if(x_coor(jj) > -0.5d0 .and. x_coor(jj) < 0.5d0) then
      u_old(jj)  = 1.d0
      u_oldSL(jj) = 1.d0
    else
      u_old(jj)  = 0.d0
      u_oldSL(jj) = 0.d0
    endif
  enddo
  u_all(:, 0) = u_old
  write(*, 130)
130 format(f14.8, 2x, I4, f14.8, f14.8)

  do it = 1, N_step

!-----------------------------------------
! Apply the periodic boundary condition
    u_old(0)   = u_old(N)  
    u_old(N+1) = u_old(1)  
    tt = it * dt

!-----------------------------------------------------------
! Find the numerical solution u_N()  |  Lax-Wendroff
      do jj = 1, N
          u_N(jj) = u_old(jj) - CFL/2.d0 * (u_old(jj+1) - u_old(jj-1))  &
                  + 1/2.d0 * CFL**2.d0 * (u_old(jj+1) - 2*u_old(jj) + u_old(jj-1))
      enddo

!-----------------------------------------------------------
! Find the numerical solution u_N()  |  Lax-Wendroff + SUPERBEE Limiter
  flag = 3
  call find_slope(u_oldSL, sigma_slope, dx, N, flag)
  do jj = 2, N
      u_NSL(jj) = u_oldSL(jj) - CFL*(u_old(jj) - u_old(jj-1)) &
               - CFL/2.d0 * (dx-aa*dt) * (sigma_slope(jj) - sigma_slope(jj-1))
  enddo

!-----------------------------------------------------------
! override the old velocity for the next time step
      u_old = u_N
      u_oldSL = u_NSL
      u_all(:,it) = u_old

!-----------------------------------------------------------
! Find the modified analytical solution u_M()
    if (abs(tt - 8.d0) .lt. tolerance) then
      write(*,*) tt
      do jj = 1, N
        if (x_coor(jj) .gt. -0.5d0+aa*tt .and. x_coor(jj) .lt. 0.5d0+aa*tt) then
      !    write(*,*) "x_coor = ", jj, x_coor(jj)
          do nn = 1, mode
            if (mod(nn, 2) .ne. 0) then
              a_n     =  2 / (pi*nn) * (-1)**((nn-1)/2)
              k_n     =  2 * pi*nn / pp
              sigma_n =  b4 * k_n**4.d0
              beta_n  = -b3 * k_n**3.d0 
              u_M(jj) =  u_M(jj) + a_n * exp(sigma_n*tt) &
                                 * cos(k_n*(x_coor(jj) - aa*tt) + beta_n*tt)
            endif
          enddo
!          write(*,*) x_coor(jj), u_M(jj)
          u_M(jj) = u_M(jj) + 1/2.d0
        else
          u_M(jj) = 0.d0
        endif
      enddo
    endif

!-----------------------------------------------------------
! Find the analytical solution of Lax_Wendroff scheme u_AS()
    do jj = 1, N   
      if (x_coor(jj) .gt. -0.5d0+aa*tt .and. x_coor(jj) .lt. 0.5d0+aa*tt) then
        do nn = 1, mode
          if (mod(nn, 2) .ne. 0) then
            a_n = 2 / (pi*nn) * (-1)**((nn-1)/2)
            k_n = 2*pi*nn / pp
            A   = sqrt(1 - CFL**2.d0 + CFL**4.d0 + cos(k_n*dx)**2.d0 * (CFL**4.d0 - CFL**2.d0) & 
                       + cos(k_n*dx) * (2*CFL**2.d0 - 2*CFL**4.d0))
            B   = atan(((1-CFL**2.d0) * tan(k_n*aa*dt) + CFL**2.d0 * tan(k_n*aa*dt) * cos(k_n*dx) - CFL*sin(k_n*dx)) &
                       / (1-CFL**2.d0 + CFL**2.d0 * cos(k_n*dx) + CFL*tan(k_n*aa*dt) * sin(k_n*dx)))
!            B   = k_n*aa*dt - atan(CFL * sin(k_n*dx) / (1 - CFL**2.d0 + CFL**2.d0*cos(k_n*dx)))
!            write(*,*) A
            u_AS(jj) = u_AS(jj) + A**nn * a_n * cos(k_n*(x_coor(jj) - aa*nn*dt) + B*nn )
          endif
        enddo
        u_AS(jj) = u_AS(jj) + 1/2.d0
      else
        u_AS(jj) = 0.d0
      endif
    enddo

 
!-----------------------------------------------------------
! Find the analytical solution u()
      do jj = 1, N
        if (abs(tt - 8.d0) .lt. tolerance ) then
          if (x_coor(jj) .gt. -0.5+aa*tt .and. x_coor(jj) .lt. 0.5+aa*tt) then
            u(jj) = 1.d0
          else
            u(jj) = 0.d0
          endif
        endif
      enddo

!-----------------------------------------------------------
! Find L1_error
    if (abs(tt - 8.d0) .lt. tolerance) then
      N_L1 = 1.d0 / dx + 1.d0
      do jj = 1, N
        if (x_coor(jj) .gt. -0.5d0+aa*tt .and. x_coor(jj) .lt. 0.5d0+aa*tt) then
          L1_error_u_N  = L1_error_u_N  + (abs(u_N(jj)  - u(jj)))**2.d0 / N_L1
          L1_error_u_M  = L1_error_u_M  + (abs(u_M(jj)  - u(jj)))**2.d0 / N_L1
          L1_error_u_AS = L1_error_u_AS + (abs(u_AS(jj) - u(jj)))**2.d0 / N_L1
        endif
      enddo
    endif

!-----------------------------------------------------------
! Output the solution and L1_error
    if (abs(tt - 8.d0) .lt. tolerance) then
      do jj = N/2, N
        write(20, 150) x_coor(jj), u(jj), u_M(jj), u_N(jj), u_AS(jj)
        write(60, 170) x_coor(jj), u(jj), u_NSL(jj)
      enddo
    endif
  enddo
150 format(5(2x, f14.8))
160 format(I0, f14.8)
170 format(3(2x, f14.8))

  write(30,160) N, L1_error_u_N
  write(40,160) N, L1_error_u_M
  write(50,160) N, L1_error_u_AS

  write(*,160) N, L1_error_u_N
  write(*,160) N, L1_error_u_M
  write(*,160) N, L1_error_u_AS
deallocate(u, u_old, u_N, x_coor)
deallocate(u_M, u_AS, u_NSL, u_oldSL)
deallocate(u_all)

!-------------------------------------------------------------------------------------------------------
!-------------------------Grid Variables
! N                 = grid resolution
! dx                = space step
! dt                = temporal step
! L                 = spatial domain
! x_coor            = spatial coordinate
! T_end             = 
! u_all(:,:)        = store all solution on time and space
!
!-------------------------Equation Variables
! CFL               = Von Neumann stability condition
! aa                = sound speed
! b2                = numerical viscosity 
! b3                = numerical dispersion error(phase) and dissipation error(amplitude)
! b4                = hyper numerical viscosity
! pp                = Period (-P/2, +P/2) 
!
!-------------------------Loop Variables
! jj                = space iteration count
! tt                = Current time size
! it                = temporal iteration count
! nn                = Number of modes
!
!-------------------------The argu for Solution of Modified Equation 
! u_M               = Numerical solution of ME
! sigma_n           = Amplitude coefficient
! beta_n            = Phase coefficient
!
!-------------------------The argu for Analytical Solution of Scheme
! u_AS              = Analytical solution of the Scheme
! A                 = exp(-sigma*dt)
! B                 = beta*dt
! k_n               = coefficient
! a_n               = Furier coefficient
! mode              = Furier series argu

!-------------------------The argu for numerical solution with SUPERBEE Limiter
! u_NSL              = solution
! u_oldSL           = last time step solution
! sigma_slope       = slope for Limiter
! flag              = choose which Limiter
!
!-------------------------L1_error
! L1_error          
! N_L1              = Number of grid points needed for find L1 error
!
!-------------------------------------------------------------------------------------------------------

end Program
