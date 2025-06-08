Program maig
  implicit none
  integer :: N, jj, N_step, it, end_n 
  real*8  :: dx, dt, L, tt, b2, b3, b4, A1, A2, k1, k2
  real*8  :: CFL, aa, pi, T_end, sigma_1, sigma_2, beta_1, beta_2
  real*8, allocatable :: u(:), x_coor(:), u_old(:), u_N(:), u_M(:)

!-----------------------------------------------------------
! Declare inital value
  L  = 2.d0
  CFL = 0.5d0
  N = 20
  end_n = N+2 
  aa = 1.d0
  dx = L / real(N)
  dt = CFL * dx / aa
  pi = 4.d0*atan(1.d0)
  T_end = 5.0 * L / aa                 ! I don't know how to determine the end time
  N_step = T_end / dt
  b2 = aa*dx/2.d0 * (1-CFL)
  b3 = -aa*dx**2/6.d0 * (1 - 3*CFL + 2*CFL**2.d0)
  b4 = aa*dx**3/24.d0 *(1 - 7*CFL + 12*CFL**2 - 6*CFL**3.d0)
  A1 = 1.d0
  A2 = 0.6d0
  k1 = 1.d0
  k2 = 3.d0
  sigma_1 = -b2*k1**2.d0 + b4*k1**2.d0
  sigma_2 = -b2*k2**2.d0 + b4*k2**2.d0
  beta_1  = -b3*k1**3.d0
  beta_2  = -b3*k2**3.d0
!-----------------------------------------------------------
! Print the variables
  open(10, file='output.dat', status = 'unknown')
  open(20, file='t1.txt', status = 'unknown')
  open(30, file='t2.txt', status = 'unknown')
  write(*, 100) 'L      =', L
  write(10,100) 'L      =', L
  
  write(*, 100) 'CFL    =', CFL
  write(10,100) 'CFL    =', CFL
  
  write(*, 100) 'dx     =', dx
  write(10,100) 'dx     =', dx
  
  write(*, 100) 'dt     =', dt
  write(10,100) 'dt     =', dt
  
  write(*, 100) 'aa     =', aa
  write(10,100) 'aa     =', aa
  
  write(*, 100) 'T_end  =', T_end
  write(10,100) 'T_end  =', T_end
  
  write(*, 110) 'N_step =', N_step
  write(10,110) 'N_step =', N_step

  write(*, 100) 'pi     =', pi
  write(10,100) 'pi     =', pi

  write(*, 100) 'b2     =', b2
  write(10,100) 'b2     =', b2
  write(*, 100) 'b3     =', b3
  write(10,100) 'b3     =', b3
  write(*, 100) 'b4     =', b4
  write(10,100) 'b4     =', b4

  write(*, 130)

  write(*, '(A)') "The first mode:"
  write(10,'(A)') "The first mode:"
  write(*,  100) 'A1 = ', A1
  write(10, 100) 'A1 = ', A1

  write(*,  100) 'beta_1 = ', beta_1
  write(10, 100) 'beta_1 = ', beta_1

  write(*,  100) 'sigma_1 = ', sigma_1
  write(10, 100) 'sigma_1 = ', sigma_1

  write(*,  100) 'k1 = ', k1
  write(10, 100) 'k1 = ', k1

  write(*, 130)

  write(*, '(A)') "The second mode:"
  write(10,'(A)') "The second mode:"
  write(*,  100) 'A2 = ', A2
  write(10, 100) 'A2 = ', A2

  write(*,  100) 'beta_2 = ', beta_2
  write(10, 100) 'beta_2 = ', beta_2

  write(*,  100) 'sigma_2 = ', sigma_2
  write(10, 100) 'sigma_2 = ', sigma_2

  write(*,  100) 'k2 = ', k2
  write(10, 100) 'k2 = ', k2

allocate(u_old(0:N+1), u_N(0:N+2), u(0:N+2))
allocate(x_coor(0:N+1), u_M(0:N+1))
!-----------------------------------------------------------
! Find the numerical solution u_N()
  x_coor(0)   = 0.d0
  x_coor(N+1) = L
  u_N   = 0.d0

!-----------------------------------------------------------
! Intial B.C 
! tt = 0
  do jj = 1, N
    x_coor(jj) = (jj-0.5) * dx
  enddo
  do jj = 0, N+1
    u_old(jj) = dsin(pi*x_coor(jj)) + 0.6*dsin(3*pi*x_coor(jj))
!    write(*, 130) tt, jj, x_coor(jj), u_old(jj)
  enddo
  write(*, 130)

!-----------------------------------------------------------
! Find the numerical solution u_N()
  do it = 1, N_step
    tt = it * dt

    do jj = 1, N+1
      if(jj == 1 .or. jj == N+1) then
        u_N(jj) = u_old(jj) * (1 - 2*CFL) + u_old(jj-1) * 2*CFL
      elseif((jj .gt. 1) .and. (jj .lt. N+1)) then
        u_N(jj) = u_old(jj) * (1 - CFL) + u_old(jj-1) * CFL
      endif
    enddo

!-----------------------------------------------------------
! override the old velocity for the next time step
    do jj = 0, N+1
      if(jj == 0) then
        u_N(jj) = u_N(N+1)             ! Apply periodic B.C.
!        write(*,*) u_N(0), u_N(N+1)
      else
        u_old(jj) = u_N(jj)
      endif 
    enddo

!-----------------------------------------------------------
! Find the modified analytical solution u_M()
    do jj = 0, N+1
      u_M(jj) = dexp(sigma_1 * tt) * A1 * dcos(k1*pi*(x_coor(jj) - aa*tt) + beta_1*tt) &
              + dexp(sigma_2 * tt) * A2 * dcos(k2*pi*(x_coor(jj) - aa*tt) + beta_2*tt) 
    enddo

!-----------------------------------------------------------
! Find the analytical solution u()
      do jj = 0, N+1
        u(jj) = dsin(pi*(x_coor(jj) - aa*tt))  + 0.6*dsin(3*pi*(x_coor(jj) - aa*tt)) 
        if(tt == 0.5d0) then
          write(*,   140) tt, jj, x_coor(jj), u_N(jj), u(jj), u_M(jj)
          write(20,  140) tt, jj, x_coor(jj), u_N(jj), u(jj), u_M(jj)

        elseif(tt == 1.5d0) then
          write(*,   140) tt, jj, x_coor(jj), u_N(jj), u(jj), u_M(jj)
          write(30,  140) tt, jj, x_coor(jj), u_N(jj), u(jj), u_M(jj)
        endif
      enddo
  enddo

deallocate(u, u_old, u_N, x_coor, u_M)
100 format(A9, f14.8)
110 format(A9, I6)
120 format(f14.8, I3, f14.8)
130 format(f14.8, I3, f14.8, f14.8)
140 format(f14.8, I3, 4f14.8)
!-----------------------------------------------------------
! Variables
! N                     = grid resolution
! dx                    = space step
! dt                    = temporal step
! CFL                   = Von Neumann stability condition
! aa                    = sound speed
! L                     = spatial domain
! x_coor                = spatial coordinate
! jj, it                = iteration count
! tt                    = temporal location
! b2                    = numerical viscosity 
! b3                    = numerical dispersion error(phase) and dissipation error(amplitude)
! b4                    = hyper numerical viscosity
! a_m1                  = numerical sound speed for mode1 i.e k = 1 
! a_m2                  = numerical sound speed for mode1 i.e k = 3 
!-----------------------------------------------------------

end Program
