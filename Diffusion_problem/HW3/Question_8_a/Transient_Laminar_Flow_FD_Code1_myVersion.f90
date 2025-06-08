Program main
implicit none

! This grid is abou the whole channel that is 2*aL, i.e. 2.0d0
! on wall implementation, center difference
! Euler time integration

  integer, dimension(1:5)  :: N                               ! Number of elements in the y-direction(spatial resolution), and N+1
!  integer, dimension(1:2)  :: file                            ! 
  integer, dimension(3)   :: Nstep                           ! The number of time steps
  real*8,  dimension(3)   :: output_time                     ! The output at Non-dimension time    

! Input parameter 
  real*8,  parameter      :: u0 = 1.0d0                      ! The maximum velocity
  real*8,  parameter      :: aL = 1.0d0                      ! The half of the channel width
  real*8,  parameter      :: nu = 0.1d0                     ! The kinematic viscosity
  real*8,  parameter      :: rou = 1.d0                      ! Fluid density

  real*8,  allocatable       :: u_old(:), u(:)
  real*8                     :: lamda_k, C_k
  real*8                     :: u_t, u_t_m                   ! single transient theoretical solution
  real*8                     :: u_theory, u_theory_m         ! Linear add for single transient theoretical solution
  real*8                     :: ycc                          ! Y_coordinate
  real*8                     :: L1_N                         ! L1,N(t) error
  real*8                     :: L2_N                         ! L2,N(t) error
  real*8                     :: CFL                          ! Dimensionaless parameter, alpha
  integer, dimension(1:5)    :: Nsteps                       ! The number of time steps
  integer                    :: it                           ! The loop para about time
  integer                    :: jj                           ! The loop para about space
  integer                    :: ii                           ! The loop para about outpur
  integer                    :: kk                           ! The sequence para about Fourier series
  integer                    :: cc                           ! Iteration for gird resolution
  integer                    :: length_output                ! The length of output_time array
  integer                    :: length_N                ! The length of output_time array
  real*8                     :: pi
  real*8, dimension(1:5)     :: dy, dt
  real*8                     :: tt                           ! Dimensional time
  real*8                     :: Tend                         ! The total simulation time(Dimensional time)
  real*8                     :: t_dimensionless              ! dimensionless time
  real*8                     :: tau_analytical               ! shear stress for analytical solution 
  real*8                     :: tau_numerical                ! shear stress for numerical solution 
  real*8                     :: eps                          ! an appropriate tolerance
  real*8                     :: nu_m                         ! The kinematic viscosity
  real*8                     :: CFL_m                        ! The kinematic viscosity
  real*8                     :: g_x                          ! The kinematic viscosity
  integer                    :: y1, y2

  N = (/ 8, 16, 32, 64, 128 /)
  length_N = size(N)
  length_output = size(output_time)
  pi = 4.0d0 * atan(1.0d0)
  eps = 1.0d-5
  g_x = 2.d0*u0*nu/aL**2.d0
  nu_m = 0.d0

  write(*, '(A,F8.5)') 'nu = ', nu
  write(*, '(A,F8.5)') 'aL  = ', aL
  write(*,'(A,I2)') 'length = ', length_output

  tt = 0d0;
  Tend = 5.0d0 * aL**2 / nu                            ! The factor '5.0' denotes the last dimensionless time                     

  write(*,*) ' '
  write(*, '(A,F8.5)') 'Tend   = ', Tend
  open(10, file='L1_error_0.2.txt', status='unknown')
  open(20, file='L2_error_0.2.txt', status='unknown')
  open(30, file='L1_error_1.0.txt', status='unknown')
  open(40, file='L2_error_1.0.txt', status='unknown')
!-------------------------------------------------------------------------------
! Numerical solution 
do cc = 1, length_N
!-----------------------------------------------------
! Initial Condition
  allocate(u_old(N(cc)))
  allocate(u(N(cc)))
  dy(cc) = 2.0d0 * aL / real(N(cc),8)
  dt(cc) = 0.32d0 * dy(cc)**2.d0 / nu
  CFL = dt(cc) * nu / dy(cc)**2.d0
  Nsteps(cc) = Tend / dt(cc)
  u = 0.0d0
  u_t = 0.d0
  u_theory = 0.d0
  u_old = u
  tt =  0.d0
  t_dimensionless =  0.d0
  write(*, '(A,F8.5)') 'CFL = ', CFL
  write(*, '(A,I0)')   'Nsteps = ', Nsteps(cc)
  write(*, '(A,F8.5)') 'dy  = ', dy(cc)
  write(*, '(A,F8.5)') 'dt  = ', dt(cc)
  write(*, '(A,I4)')   'N   = ', N(cc)
  write(*,*) ' '

  do it = 1, Nsteps(cc)
    do jj = 2, N(cc)
      u(jj) = u_old(jj) + CFL*(u_old(jj-1) - 2.d0*u_old(jj) + u_old(jj+1)) & 
              + dt(cc)*g_x
    enddo

    do jj = 2, N(cc)
      u_old(jj) = u(jj)
    enddo

    tt = tt + dt(cc)
    t_dimensionless = nu * tt / aL**2.d0
!-------------------------------------------------------------------------------
! analytical solution
! nu and nu_m
    if((t_dimensionless .gt. 0.2d0-eps) .and. (t_dimensionless .lt. 0.2d0+eps)) then
      L1_N = 0.d0
      L2_N = 0.d0
      do jj = 1, N(cc)+1
          ycc        = -aL + real(jj-1)*dy(cc) 
          u_theory   = (1.d0 - ycc**2.d0 / aL**2.d0)                                  ! steady solution
          u_theory_m = (1.d0 - ycc**2.d0 / aL**2.d0)                                  ! steady solution

          do kk = 0, 100
              nu_m       = nu*(1 + (CFL/2.d0 - 1.d0/12.d0) * 4*(lamda_k*aL/N(cc))**2.d0 &
                                 + (CFL**2.d0/3.d0 - CFL/12.d0 + 1.d0/360.d0) * 16*(lamda_k*aL/N(cc))**4.d0 )
              lamda_k    = (0.5d0 + real(kk)) * pi / aL
              C_k        = -4.0d0 * (-1.0)**kk / (lamda_k*aL)**3.d0
              u_t        = C_k * exp(-(lamda_k**2.d0 * nu * tt)) * cos(lamda_k * ycc)
              u_t_m      = C_k * exp(-(lamda_k**2.d0 * nu_m * tt)) * cos(lamda_k * ycc)
              u_theory   = u_theory + u_t                                               ! steady solution + transient solution
              u_theory_m = u_theory_m + u_t_m                                               ! steady solution + transient solution
          enddo
          if (jj .ge. 2) then
            L1_N = L1_N + abs(u(jj) - u_theory) / (N(cc)-1)
            L2_N = L2_N + (u(jj) - u_theory)**2.d0 / (N(cc)-1)
          endif
      enddo
      write(10, 100)  N(cc), t_dimensionless, L1_N / u0       ! Find truncation error 
      write(*,  100)  N(cc), t_dimensionless, L1_N / u0       ! Find truncation error 
      write(20, 100)  N(cc), t_dimensionless, sqrt(L2_N) / u0 
      write(*,  100)  N(cc), t_dimensionless, sqrt(L2_N) / u0 
    endif

    if((t_dimensionless .gt. 1.d0-eps) .and. (t_dimensionless .lt. 1.d0+eps)) then
      L1_N = 0.d0
      L2_N = 0.d0
      do jj = 1, N(cc)+1
          ycc        = -aL + real(jj-1)*dy(cc) 
          u_theory   = (1.d0 - ycc**2.d0 / aL**2.d0)                                  ! steady solution
          u_theory_m = (1.d0 - ycc**2.d0 / aL**2.d0)                                  ! steady solution

          do kk = 0, 100
              nu_m       = nu*(1 + (CFL/2.d0 - 1.d0/12.d0) * 4*(lamda_k*aL/N(cc))**2.d0 &
                                 + (CFL**2.d0/3.d0 - CFL/12.d0 + 1.d0/360.d0) * 16*(lamda_k*aL/N(cc))**4.d0 )
              lamda_k    = (0.5d0 + real(kk)) * pi / aL
              C_k        = -4.0d0 * (-1.0)**kk / (lamda_k*aL)**3.d0
              u_t        = C_k * exp(-(lamda_k**2.d0 * nu * tt)) * cos(lamda_k * ycc)
              u_t_m      = C_k * exp(-(lamda_k**2.d0 * nu_m * tt)) * cos(lamda_k * ycc)
              u_theory   = u_theory + u_t                                               ! steady solution + transient solution
              u_theory_m = u_theory_m + u_t_m                                               ! steady solution + transient solution
          enddo
          if (jj .ge. 2) then
            L1_N = L1_N + abs(u(jj) - u_theory) / (N(cc)-1)
            L2_N = L2_N + (u(jj) - u_theory)**2.d0 / (N(cc)-1)
          endif
      enddo
      write(30, 100)  N(cc), t_dimensionless, L1_N / u0       ! Find truncation error 
      write(*,  100)  N(cc), t_dimensionless, L1_N / u0       ! Find truncation error 
      write(40, 100)  N(cc), t_dimensionless, sqrt(L2_N) / u0 
      write(*,  100)  N(cc), t_dimensionless, sqrt(L2_N) / u0 
    endif
  enddo
  deallocate(u_old, u)
enddo 
! 2x        skips 2 spaces at the start of the line
! 4F16.12   single-precision, total 16 characters and 12 decimal places
100     format(2x,I4, 2f16.12)
endprogram
