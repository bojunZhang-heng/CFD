Program main
implicit none

! This grid is abou the whole channel that is 2*aL, i.e. 2.0d0
! on wall implementation, center difference
! Euler time integration

!  integer, parameter      :: N = 8                           ! Number of elements in the y-direction(spatial resolution), and N+1
  integer, parameter      :: N = 16                          ! Number of elements in the y-direction(spatial resolution), and N+1
!  integer, parameter      :: N = 32                          ! Number of elements in the y-direction(spatial resolution), and N+1
!  integer, parameter      :: N = 2**10.d0                          ! Number of elements in the y-direction(spatial resolution), and N+1
  integer, dimension(3)   :: Nstep                           ! The number of time steps
  real*8,  dimension(3)   :: output_time                     ! The output at Non-dimension time    
  real*8,  dimension(4)   :: dt                              ! Time step size for Euler integration

! Input parameter 
  real*8,  parameter      :: u0 = 1.0d0                      ! The maximum velocity
  real*8,  parameter      :: aL = 1.0d0                      ! The half of the channel width
  real*8,  parameter      :: nu = 0.1d0                      ! The kinematic viscosity
  real*8,  parameter      :: rou = 1.d0                      ! Fluid density

  real*8,  dimension(1:N+1)  :: u_old, u
  real*8                     :: lamda_k, C_k
  real*8                     :: u_t, u_t_m                   ! single transient theoretical solution
  real*8                     :: u_theory, u_theory_m         ! Linear add for single transient theoretical solution
  real*8                     :: ycc                          ! Y_coordinate
  real*8,  dimension(1:4)    :: CFL                          ! Dimensionaless parameter, alpha
  integer, dimension(1:4)    :: Nsteps                       ! The number of time steps
  integer, dimension(1:4)    :: file                         ! The number of time steps
  integer                    :: it                           ! The loop para about time
  integer                    :: jj                           ! The loop para about space
  integer                    :: ii                           ! The loop para about outpur
  integer                    :: kk                           ! The sequence para about Fourier series
  integer                    :: cc                           ! Iteration for CFL
  integer                    :: length_output                ! The length of output_time array
  integer                    :: length_CFL                   ! The length of CFL array
  real*8                     :: dy, pi
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

  output_time = (/ 0.2d0, 1.0d0, 5.0d0 /)
  CFL         = (/ 0.50d0, 0.505d0, 0.51d0, 0.52d0 /)
  file        = (/ 10, 20, 30, 40 /)
  length_output = size(output_time)
  length_CFL    = size(CFL)
  pi = 4.0d0 * atan(1.0d0)
  dy = 2.0d0 * aL / real(N,8)
  eps = 1.0d-8
  g_x = 2.d0*u0*nu/aL**2.d0
  nu_m = 0.d0

!  CFL = dt * nu / dy**2.d0

  write(*, '(A,I2)')   'N   = ', N
  write(*, '(A,F8.5)') 'nu = ', nu
  write(*, '(A,F8.5)') 'aL  = ', aL
  write(*, '(A,F8.5)') 'dy  = ', dy
  write(*,'(A,I2)') 'length = ', length_output

  tt = 0d0;
  Tend = 5.0d0 * aL**2 / nu                            ! The factor '5.0' denotes the last dimensionless time                     

  write(*,*) ' '
  write(*, '(A,F8.5)') 'Tend   = ', Tend
  open(10, file='velocity_1.txt', status='unknown')
  open(20, file='velocity_2.txt', status='unknown')
  open(30, file='velocity_3.txt', status='unknown')
  open(40, file='velocity_4.txt', status='unknown')

!-------------------------------------------------------------------------------
! Get the space index at ycc = 0 and ycc = -0.5L
  do jj = 1, N+1
    ycc = -aL + real(jj-1)*dy 
    if(ycc == -0.5d0*aL) y1 = jj
    if(ycc == 0.d0)      y2 = jj
  enddo
  write(*,'(A,I0)') 'y1 = ', y1
  write(*,'(A,I0)') 'y2 = ', y2


!-------------------------------------------------------------------------------
! Numerical solution 
do cc = 1, length_CFL
  tt = 0.d0
  dt(cc) = CFL(cc) * dy**2.d0 / nu
  Nsteps(cc) = Tend / dt(cc)
  write(*, '(A,F8.5)') 'CFL = ', CFL(cc)
  write(*, '(A,F8.5)') 'dt  = ', dt(cc)
  write(*, '(A,I0)')   'Nsteps = ', Nsteps(cc)
!------------------------------------------------------
! Initial Condition
  u = 0.0d0
  u_t = 0.d0
  u_theory = 0.d0
  u_old = u;
  do it = 1, Nsteps(cc)
    do jj = 2, N
      u(jj) = u_old(jj) + CFL(cc)*(u_old(jj-1) - 2.d0*u_old(jj) + u_old(jj+1)) & 
              + dt(cc)*g_x
    enddo

    do jj = 2, N
      u_old(jj) = u(jj)
    enddo

    tt = tt + dt(cc)
    t_dimensionless = nu * tt / aL**2.d0

!-------------------------------------------------------------------------------
! analytical solution
! nu and nu_m
    do jj = 1, N+1
      ycc        = -aL + real(jj-1)*dy 
      u_theory   = (1.d0 - ycc**2.d0 / aL**2.d0)                                  ! steady solution
      u_theory_m = (1.d0 - ycc**2.d0 / aL**2.d0)                                  ! steady solution

      if(jj == y2) then
        do kk = 0, 100
          nu_m       = nu*(1 + (CFL(cc)/2.d0 - 1.d0/12.d0) * 4*(lamda_k*aL/N)**2.d0 &
                             + (CFL(cc)**2.d0/3.d0 - CFL(cc)/12.d0 + 1.d0/360.d0) * 16*(lamda_k*aL/N)**4.d0 )
          lamda_k    = (0.5d0 + real(kk)) * pi / aL
          C_k        = -4.0d0 * (-1.0)**kk / (lamda_k*aL)**3.d0
          u_t        = C_k * exp(-(lamda_k**2.d0 * nu * tt)) * cos(lamda_k * ycc)
          u_t_m      = C_k * exp(-(lamda_k**2.d0 * nu_m * tt)) * cos(lamda_k * ycc)
          u_theory   = u_theory + u_t                                                   ! steady solution + transient solution
          u_theory_m = u_theory_m + u_t_m                                               ! steady solution + transient solution
        enddo

        if(t_dimensionless .lt. 5.0) then
          write(file(cc), 100)  ycc, t_dimensionless,   u(jj)/u0              
          write(* , 100)        ycc, t_dimensionless,   u(jj)/u0              
        endif
      endif
    enddo
  enddo
enddo
100     format(2x, f24.8, 4x, f24.8, 4x, f24.8, 4x)
endprogram
