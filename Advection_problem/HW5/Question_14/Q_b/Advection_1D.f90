!------------------------------------------------------------------------------
! Lax-Wendroff sechme and Upwind scheme
! N = 80  CFL = 0.5, 1.0, 1.2
! 
!------------------------------------------------------------------------------

Program main
  implicit none
  integer :: jj, it, ii, kk, unit_id1, unit_id2, length_N, length_CFL
  real*8  :: L, tt, A1, A2, k1, k2, L1_error1, L1_error2
  real*8  :: aa, pi, T_end
  real*8,  allocatable :: u(:), x_coor(:), u_old(:), u_N1(:), u_N2(:), u_M(:)
  real*8,  allocatable :: dx(:), dt(:), b2(:), b3(:), b4(:), CFL(:) 
  real*8,  allocatable :: sigma_1(:), sigma_2(:), beta_1(:), beta_2(:)
  integer, allocatable :: N(:), N_step(:)
  character(len=20)    :: filename1, filename2

!-----------------------------------------------------------
! Declare inital value
  allocate(N(5), CFL(3))
  N   = (/ 20, 30, 40, 50, 80 /)
  CFL = (/ 0.5, 1.0, 1.2 /)
  length_N = size(N)
  length_CFL = size(CFL)

  open(10, file='output.dat', status = 'unknown')

  L  = 2.d0
  aa = 1.d0
  pi = 4.d0*atan(1.d0)
  T_end = 3.0 * L / aa       
  write(*, 100) 'L      =', L
  write(10,100) 'L      =', L


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
  do kk = 1, length_CFL
    unit_id2 = 70 * kk + 10
    write(filename2, '(A, f5.3, A)') 'CFL', CFL(kk), '.txt'
    open(unit_id2, file=filename2, status = 'unknown')

  do ii = 1, length_N
    write(*,*) '***********************************************'
    write(* , '(A,I0,A,I0)') 'N(', ii, ') = ', N(ii)
    write(10, '(A,I0,A,I0)') 'N(', ii, ') = ', N(ii)
    unit_id1 = 10 * ii + 10
  
    write(filename1, '(A, I0, A)') 'N', N(ii), '.txt'
    open(unit_id1, file=filename1, status = 'unknown')

  
    allocate(u_old(0 : N(ii)+1), u_N1(0 : N(ii)+1), u(0 : N(ii)+1), u_N2(0 : N(ii)+1))
    allocate(x_coor(0 : N(ii)+1), u_M(0 : N(ii)+1))
  
    dx(ii) = L / real(N(ii))
    dt(ii) = CFL(kk) * dx(ii) / aa
    N_step(ii) = int(T_end / dt(ii))
  
    b2(ii) = 0.d0
    b3(ii) = (-aa * dx(ii)**2.d0)/6.0 * (1-CFL(kk)**2.d0)
    b4(ii) = (-aa * dx(ii)**3.d0)/8.01 * (1-CFL(kk)**2.d0)*CFL(kk)
  
    A1 = 1.d0
    A2 = 0.6d0
    k1 = 1.d0
    k2 = 3.d0
  
    sigma_1(ii)  = b4(ii) * (k1*pi)**4.d0
    sigma_2(ii)  = b4(ii) * (k2*pi)**4.d0
    beta_1(ii)   = k1*pi*aa - k1*pi - b3(ii)  * (k1*pi)**3.d0
    beta_2(ii)   = k2*pi*aa - k2*pi - b3(ii)  * (k2*pi)**3.d0
  
    write(* , '(A,I0,A,F10.6)') 'CFL(', kk, ') = ', CFL(kk)
    write(10, '(A,I0,A,F10.6)') 'CFL(', kk, ') = ', CFL(kk)
  
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
    u_N1   = 0.d0
    u_N2   = 0.d0
    u_old = 0.d0
  
  !-----------------------------------------------------------
  ! Intial B.C 
   tt = 0.d0
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
      L1_error1 = 0.d0
      L1_error2 = 0.d0 
  !----------------------------------------- 
  ! Apply the periodic boundary condition
      u_old(0) = u_old(N(ii)) 
      u_old(N(ii) + 1) = u_old(1) 
      tt = it * dt(ii) 

      do jj = 1, N(ii) 
        u_N1(jj) = u_old(jj) * (1 - CFL(kk)**2.d0) + u_old(jj-1) * 1/2.d0*CFL(kk) * (1+CFL(kk)) & 
                 + u_old(jj+1) * 1/2.d0*CFL(kk) * (CFL(kk) - 1) 

        u_N2(jj) = u_old(jj) * (1 - CFL(kk)) + u_old(jj-1) * CFL(kk) 
      enddo
  
  !-----------------------------------------------------------
  ! override the old velocity for the next time step
      do jj = 1, N(ii)
          u_old(jj) = u_N1(jj)
          u_old(jj) = u_N2(jj)
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
!          write(unit_id1,  140) tt, jj, x_coor(jj), u_N1(jj), u(jj), u_M(jj)
          L1_error1 = L1_error1 + (abs(u_N1(jj) - u(jj)))**2.d0 / real(N(ii))
          L1_error2 = L1_error2 + (abs(u_N2(jj) - u(jj)))**2.d0 / real(N(ii))
      enddo
      if(N(ii) == 80) then
        write(unit_id2, 160) tt, N(ii), CFL(kk), L1_error1, L1_error2
      endif
    enddo
!    write(70, 150) N(ii), L1_error

  deallocate(u, u_old, u_N1, u_N2, x_coor, u_M)
  enddo
  enddo
deallocate(sigma_1, sigma_2, beta_1, beta_2)
deallocate(b2, b3, b4, dx, dt, N_step)
deallocate(N, CFL)
100 format(A9, f14.8)
110 format(A9, I6)
!120 format(f14.8, 2x, I4, f14.8)
130 format(f14.8, 2x, I4, f14.8, f14.8)
!140 format(f14.8, 2x, I4, 4f14.8)
!150 format(I0, 2x, f14.8)
160 format(f14.8, 2x, I0, 2x, f14.8, 2x, f24.20, 2x, f24.20)
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
! KK                    = CFL iteration
! tt                    = temporal location
! b2                    = numerical viscosity 
! b3                    = numerical dispersion error(phase) and dissipation error(amplitude)
! b4                    = hyper numerical viscosity
! sigma_1               = (Amplitude) dissipation error 
! beta_1                = (Phase) dispersion error
! L1_error              = error
! u_N1                  = Lax-Wendroff scheme
! u_N2                  = Upwind sceme
!-----------------------------------------------------------

end Program




