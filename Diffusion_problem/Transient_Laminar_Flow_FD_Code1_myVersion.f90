Program main
implicit none

! This grid is abou the whole channel that is 2*aL, i.e. 2.0d0
! on wall implementation, center difference
! Euler time integration

  integer, parameter      :: N = 8                           ! Number of elements in the y-direction(spatial resolution), and N+1
!  integer, parameter      :: N = 16                          ! Number of elements in the y-direction(spatial resolution), and N+1
!  integer, parameter      :: N = 32                          ! Number of elements in the y-direction(spatial resolution), and N+1
  integer, dimension(3)   :: mm                              ! The number of time steps
  real*8,  dimension(3)   :: output_time                     ! The output at Non-dimension time    
  real*8                  :: dt                              ! Time step size for Euler integration

! Input parameter 
  real*8,  parameter      :: u0 = 1.0d0                      ! The maximum velocity
  real*8,  parameter      :: aL = 1.0d0                      ! The half of the channel width
  real*8,  parameter      :: anu = 0.1d0                     ! The kinematic viscosity
  real*8,  parameter      :: rou = 1.d0                      ! Fluid density

  real*8,  dimension(1:N+1)  :: u_old, u
  real*8                     :: theta_0, theta_1, theory_1, ss
  real*8                     :: ycc                          ! Y_coordinate
  real*8                     :: CFL                          ! Dimensionaless parameter
  integer                    :: Nsteps                       ! The number of time steps
  integer                    :: it                           ! The loop para about time
  integer                    :: jj                           ! The loop para about space
  integer                    :: ii                           ! The loop para about outpur
  integer                    :: kk                           ! The sequence para about Fourier series
  integer                    :: length_output                ! The length of output_time array
  real*8                     :: dy, pi, tt
  real*8                     :: Tend                         ! The total simulation time
  real*8                     :: t_dimensionless              ! dimensionless time
  real*8                     :: tau_analytical               ! shear stress for analytical solution 
  real*8                     :: tau_numerical                ! shear stress for numerical solution 
  real*8                     :: eps                          ! an appropriate tolerance

  output_time = (/ 0.2d0, 1.0d0, 5.0d0 /)
  length_output = size(output_time)
  pi = 4.0d0 * atan(1.0d0)
  dy = 2.0d0 * aL / real(N,8)
  dt = 0.32d0 * dy**2.d0 / anu
  eps = 1.0d-8

! Initial Condition
  u = 0.0d0
  u_old = u;
  CFL = dt * anu / dy**2.d0

  write(*, '(A,I2)')   'N   = ', N
  write(*, '(A,F8.5)') 'anu = ', anu
  write(*, '(A,F8.5)') 'dt  = ', dt
  write(*, '(A,F8.5)') 'aL  = ', aL
  write(*, '(A,F8.5)') 'dy  = ', dy
  write(*, '(A,F8.5)') 'CFL = ', CFL

  tt = 0d0;
  Tend = 3.0d0 * aL**2 / anu
  Nsteps = Tend / dt

  write(*,*) ' '
  write(*, '(A,F8.5)') 'Tend   = ', Tend
  write(*, '(A,F8.5)') 'dt     = ', dt
  write(*, '(A,I0)')   'Nsteps = ', Nsteps

!-------------------------------------------------------------------------------
! Numerical solution
  do it = 1, Nsteps
    do jj = 2, N
      u(jj) = u_old(jj) + CFL*(u_old(jj-1) - 2.d0*u_old(jj) + u_old(jj+1)) & 
              + dt*2.d0*u0*anu/aL**2.d0
    enddo

    do jj = 2, N
      u_old(jj) = u(jj)
    enddo

    tt = tt + dt
    t_dimensionless = tt * anu / aL**2.d0
!    write(*,*) "t_dimensionless = ", t_dimensionless
      do ii = 1, length_output
        if(t_dimensionless > output_time(ii)-eps .and. t_dimensionless < output_time(ii)+eps) then
          mm(ii) = t_dimensionless * N**2.d0 / 1.28  
          write(*,*)
          write(*,*) "The dimensionless time is ", output_time(ii) 
          write(*,*) "mm = ", mm(ii) 
          write(*,*)
        endif
      enddo

!-------------------------------------------------------------------------------
! Print out solutions about shear stress
    open(10, file='shear_stress_result_N32.txt', status='unknown')
    do ii = 1, length_output
      if(it .eq. mm(ii)) then
! analytical solution
        do jj = 1, N+1
          ycc = real(jj-1)*dy -1.0d0
          theory_1 = (1.d0 - ycc**2.d0 / aL**2.d0)

          do kk = 0, 100
            ss = (0.5d0 + real(kk)) * pi
            theta_0 = 4.0d0 * (-1.0)**kk / ss**3.d0
            theta_1 = theta_0 * exp(-ss**2.d0 * t_dimensionless) * cos(ss * ycc / aL)
            theory_1 = theory_1 - theta_1
          enddo

          if(jj == 2) then
            tau_analytical = 11*theory_1  / (7*dt)

            tau_numerical= 11*u(jj)/ (7*dt)
          endif
          if(jj == 3) then
            tau_analytical = tau_analytical - 2*theory_1 / (7*dt)
            tau_analytical = tau_analytical * rou * anu

            tau_numerical = tau_numerical - 2*u(jj)/ (7*dt)
            tau_numerical = tau_numerical * rou * anu
          endif

!          write(10, 100)  ycc, theory_1, u(jj)/u0, abs(u(jj) - theory_1) / theory_1
          write(*, 100)  ycc, theory_1, u(jj)/u0, abs(u(jj) - theory_1) / theory_1

          if(jj .eq. (1+N+1)/2) then
            if((u(jj) - 0.97*u0 > -eps) .and. (u(jj) - 0.97*u0) < eps) then
              write(*,*) "The dimensionless time = ", mm(ii)
            endif
          endif
        enddo
        write(*,*) "tau_analytical is ", tau_analytical 
        write(*,*) "tau_numerical is  ", tau_numerical
        write(10,*) tau_analytical 
        write(10,*) tau_numerical
      endif
    enddo
! 2x        skips 2 spaces at the start of the line
! 4F16.12   single-precision, total 16 characters and 12 decimal places
100     format(2x, 4f16.12)
  enddo
endprogram
