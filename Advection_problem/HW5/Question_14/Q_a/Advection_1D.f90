!------------------------------------------------------------------------------
! t1 = 1.5, N = (/ 20, 40, 80, 160, 1000/)
! output N(ii), L1_error
!------------------------------------------------------------------------------

Program main
  implicit none
  integer :: jj, it, ii, unit_id, length_N
  real*8  :: L, tt, A1, A2, k1, k2, L1_error
  real*8  :: CFL, aa, pi, T_end
  real*8,  allocatable :: u(:), x_coor(:), u_old(:), u_N(:), u_M(:)
  real*8,  allocatable :: dx(:), dt(:), b2(:), b3(:), b4(:) 
  real*8,  allocatable :: sigma_1(:), sigma_2(:), beta_1(:), beta_2(:)
  integer, allocatable :: N(:), N_step(:)
  character(len=20)    :: filename

!-----------------------------------------------------------
! Declare inital value
  allocate(N(5))
  N = (/ 20, 40, 80, 160, 1000/)
  length_N = size(N)

  open(10, file='output.dat', status = 'unknown')

  L  = 2.d0
  CFL = 0.5d0
  L1_error = 0.d0
  aa = 1.d0
  pi = 4.d0*atan(1.d0)
  T_end = 5.0 * L / aa       
  write(*, 100) 'L      =', L
  write(10,100) 'L      =', L

  write(*, 100) 'CFL    =', CFL
  write(10,100) 'CFL    =', CFL

  write(*, 100) 'aa     =', aa
  write(10,100) 'aa     =', aa

  write(*, 100) 'T_end  =', T_end
  write(10,100) 'T_end  =', T_end

  write(*, 110) 'length_N  =', length_N
  write(10,110) 'length_N  =', length_N

  write(*, 100) 'pi     =', pi
  write(10,100) 'pi     =', pi

  allocate(dx(length_N), dt(length_N), N_step(length_N))
  allocate(b2(length_N), b3(length_N), b4(length_N))
  allocate(sigma_1(length_N), sigma_2(length_N), beta_1(length_N), beta_2(length_N))


!--------------------------------------------------------------------
! N loop
do ii = 1, length_N
  write(*,*) '***********************************************'
  write(* , '(A,I0,A,I0)') 'N(', ii, ') = ', N(ii)
  write(10, '(A,I0,A,I0)') 'N(', ii, ') = ', N(ii)
  unit_id = 10 * ii + 10

  write(filename, '(A, I0, A)') 'N', N(ii), '.txt'
  open(unit_id, file=filename, status = 'unknown')

  allocate(u_old(0 : N(ii)+1), u_N(0 : N(ii)+1), u(0 : N(ii)+1))
  allocate(x_coor(0 : N(ii)+1), u_M(0 : N(ii)+1))

  dx(ii) = L / real(N(ii))
  dt(ii) = CFL * dx(ii) / aa
  N_step(ii) = int(T_end / dt(ii))

  b2(ii) = 0.d0
  b3(ii) = (-aa * dx(ii)**2.d0)/6.0 * (1-CFL**2.d0)
  b4(ii) = (-aa * dx(ii)**3.d0)/8.0 * (1-CFL**2.d0)*CFL

  A1 = 1.d0
  A2 = 0.6d0
  k1 = 1.d0
  k2 = 3.d0

  sigma_1(ii)  = b4(ii) * (k1*pi)**4.d0
  sigma_2(ii)  = b4(ii) * (k2*pi)**4.d0
  beta_1(ii)   = k1*pi*aa - k1*pi - b3(ii)  * (k1*pi)**3.d0
  beta_2(ii)   = k2*pi*aa - k2*pi - b3(ii)  * (k2*pi)**3.d0


  write(* , '(A,I0,A,F10.6)') 'dx(', ii, ') = ', dx(ii)
  write(10, '(A,I0,A,F10.6)') 'dx(', ii, ') = ', dx(ii)

  write(* , '(A,I0,A,F10.6)') 'dt(', ii, ') = ', dt(ii)
  write(10, '(A,I0,A,F10.6)') 'dt(', ii, ') = ', dt(ii)

  write(* , '(A,I0,A,I0)') 'N_step(', ii, ') = ', N_step(ii)
  write(10, '(A,I0,A,I0)') 'N_step(', ii, ') = ', N_step(ii)

  write(* , '(A,I0,A,F10.6)') 'b2(', ii, ') = ', b2(ii)
  write(10, '(A,I0,A,F10.6)') 'b2(', ii, ') = ', b2(ii)

  write(* , '(A,I0,A,F10.6)') 'b3(', ii, ') = ', b3(ii)
  write(10, '(A,I0,A,F10.6)') 'b3(', ii, ') = ', b3(ii)

  write(* , '(A,I0,A,F10.6)') 'b4(', ii, ') = ', b4(ii)
  write(10, '(A,I0,A,F10.6)') 'b4(', ii, ') = ', b4(ii)

  write(*, 130)

  write(*, '(A)') "The first mode:"
  write(10,'(A)') "The first mode:"
  write(*,  100) 'A1 = ', A1
  write(10, 100) 'A1 = ', A1

  write(*,  100) 'k1 = ', k1
  write(10, 100) 'k1 = ', k1

  write(* , '(A,I0,A,F10.6)') 'beta_1(', ii, ') = ', beta_1(ii)
  write(10, '(A,I0,A,F10.6)') 'beta_1(', ii, ') = ', beta_1(ii)

  write(* , '(A,I0,A,F10.6)') 'sigma_1(', ii, ') = ', sigma_1(ii)
  write(10, '(A,I0,A,F10.6)') 'sigma_1(', ii, ') = ', sigma_1(ii)



  write(*, 130)

  write(*, '(A)') "The second mode:"
  write(10,'(A)') "The second mode:"
  write(*,  100) 'A2 = ', A2
  write(10, 100) 'A2 = ', A2
  
  write(*,  100) 'k2 = ', k2
  write(10, 100) 'k2 = ', k2

  write(* , '(A,I0,A,F10.6)') 'beta_2(', ii, ') = ', beta_2(ii)
  write(10, '(A,I0,A,F10.6)') 'beta_2(', ii, ') = ', beta_2(ii)

  write(* , '(A,I0,A,F10.6)') 'sigma_2(', ii, ') = ', sigma_2(ii)
  write(10, '(A,I0,A,F10.6)') 'sigma_2(', ii, ') = ', sigma_2(ii)

!-----------------------------------------------------------
! Find the numerical solution u_N()
  x_coor(0)   = 0.d0
  x_coor(N(ii) + 1) = L
  u_N   = 0.d0
  u_old = 0.d0

!-----------------------------------------------------------
! Intial B.C 
! tt = 0
  do jj = 1, N(ii)
    x_coor(jj) = (jj-0.5) * dx(ii)
  enddo
  do jj = 0, N(ii)+1
    u_old(jj) = dsin(pi*x_coor(jj)) + 0.6*dsin(3*pi*x_coor(jj))
!    write(*, 130) tt, jj, x_coor(jj), u_old(jj)
  enddo
  write(*, 130)

!-----------------------------------------------------------
! Find the numerical solution u_N()
  do it = 1, N_step(ii)

!-----------------------------------------
! Apply the periodic boundary condition
    u_old(0) = u_old(N(ii))  
    u_old(N(ii) + 1) = u_old(1)  
    tt = it * dt(ii)

    do jj = 1, N(ii)
        u_N(jj) = u_old(jj) * (1 - CFL**2.d0) + u_old(jj-1) * 1/2.d0*CFL * (1+CFL) + u_old(jj+1) * 1/2.d0*CFL * (-1+CFL)
    enddo

!-----------------------------------------------------------
! override the old velocity for the next time step
    do jj = 1, N(ii)
        u_old(jj) = u_N(jj)
    enddo

!-----------------------------------------------------------
! Find the modified analytical solution u_M()
    do jj = 0, N(ii)+1
      u_M(jj) = A1 * dexp(sigma_1(ii) * tt) * dsin(k1*pi*(x_coor(jj) - aa*tt) + beta_1(ii)*tt) &
              + A2 * dexp(sigma_2(ii) * tt) * dsin(k2*pi*(x_coor(jj) - aa*tt) + beta_2(ii)*tt) 
    enddo

!-----------------------------------------------------------
! Find the analytical solution u()
      do jj = 1, N(ii)
        u(jj)  = A1 * dsin(k1*pi*(x_coor(jj) - aa*tt))  + A2 * dsin(k2*pi*(x_coor(jj) - aa*tt))
        if(tt == 0.5d0) then
!          write(*,   140) tt, jj, x_coor(jj), u_N(jj), u(jj)
!          write(20,  140) tt, jj, x_coor(jj), u_N(jj), u(jj)

        elseif(tt == 1.5d0) then
!         write(*,   140) tt, jj, x_coor(jj), u_N(jj), u(jj), u_M(jj)
          write(unit_id,  140) tt, jj, x_coor(jj), u_N(jj), u(jj), u_M(jj)
          L1_error = L1_error + (abs(u_N(jj) - u(jj)))**2.d0 / real(N(ii))
          write(80, '(I4, f14.8)') N(ii), L1_error 
        endif
      enddo
  enddo
  write(70, 150) N(ii), L1_error



100 format(A9, f14.8)
110 format(A9, I6)
!120 format(f14.8, 2x, I4, f14.8)
130 format(f14.8, 2x, I4, f14.8, f14.8)
140 format(f14.8, 2x, I4, 4f14.8)
150 format(I0, 2x, f14.8)
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
! sigma_1               = (Amplitude) dissipation error 
! beta_1                = (Phase) dispersion error
! L1_error              = error
!-----------------------------------------------------------
deallocate(u, u_old, u_N, x_coor, u_M)

enddo
deallocate(sigma_1, sigma_2, beta_1, beta_2)
deallocate(b2, b3, b4, dx, dt, N_step)
deallocate(N)
end Program
