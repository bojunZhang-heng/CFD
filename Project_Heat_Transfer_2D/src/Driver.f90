program Driver
  use fftcos
  use mod_param
  use mod_variables
  use mod_fft_data
  NAMELIST/INPUT/niter, nprof, iread, re, pr, len, dt, rerror
!  implicit none

      call FFTFAX(n1, IFAXX, TRIGSX)
      open(10, file='velprof1.dat', form='formatted', &
           status='unknown')
      open(11, file='velprof2.dat', form='formatted', &
           status='unknown')

      open(25,  file='uc.dat',   form='formatted', status='unknown')
      open(100, file='eps.dat',  form='formatted', status='unknown')
      open(5,   file='data.inp', form='formatted', status='old', iostat=ios)
      open(7,   file='data.out', form='formatted', status='unknown')

      open(201, file='restart.flo', form='unformatted', &
                status='unknown') 
      open(202, file='endrun.flo',  form='unformatted', &
                status='unknown')
! DX(),DY(),arrays are not needed for dx,dy uniform
! here, x,y are locations of cell boundarys (this will change
! for dx,dy not uniform
!
      zero   = 1.0e-30
! set up Runge-Kutta coefficients
      ggg(1) =  3./2.
      ggg(2) =  0.0
      ggg(3) =  0.0

      rrr(1) = -0.5
      rrr(2) =  0.0
      rrr(3) =  0.0

      aaa(1) =  1.0
      aaa(2) =  0.0
      aaa(3) =  0.0

! input Reynolds number
      read(5, nml=INPUT, iostat=ios)
      write(7, nml=INPUT)
!---------------------------------------------------------------------------
! Initial SETUP
! if IREAD = 0, new start (else RESTART)
      IF (iread .eq. 0) THEN
        itbeg = 0
        tau   = 0.0
        dx  = len / float(n1)
        dy  = len / float(n2)
        dx2 = dx ** 2.0
        dy2 = dy ** 2.0

! fill x,y vectors (grid line coords)
        do ii = 1, n1p1
          x(ii) = float(ii)*dx - dx
        enddo
        do jj = 1,n2p1
          y(jj) = float(jj)*dy - dy
        enddo
  
        do ii = 1,n1p1
        write(21,*) ii, x(ii), y(ii)
        enddo
        
! fill the cell center x,y vectors
        xc(1) = dx / 2.
        do ii = 2, n1
          xc(ii) = xc(ii-1) + dx
        enddo

        yc(1) = dy / 2.
        do jj = 2,n2
          yc(jj) = yc(jj-1) + dy
        enddo

        do ii = 2, n1
          write(*, *) ii, xc(ii), yc(ii)
        enddo
!
! stability criterion (Peyret & Taylor, 1983)...
! there are two: eqs. 6.2.10 & 6.2.11
! for eq. 6.2.10, assume |u0| & |v0| are equal to 1.0
!
        c1 = dt * re
        write(*,*) 'Stability criterion 1 = ',c1
        if (c1 .gt. 1.0) then
          write(*,*) 'Stability criterion #1 greater than 1.0'
!        write(*,*) 'Program stopped'
!        stop
        endif

        c2 = 4*dt/(re*dx2)
        write(*,*) 'Stability criterion 2 = ',c2
        if (c2 .gt. 1.0) then
          write(*,*) 'Stability criterion #2 greater than 1.0'
!        write(*,*) 'Program stopped'
!        stop
        endif
!
! initialize all arrays at t=zero...
! (can include boundary for velocity, but make general
! for boundary velocity not equal to zero)
! for this case, hxp,hyp are initially zero
!
        do jj = 2,n2p1
          do ii = 2,n1
            qx(ii, jj) = zero
            hxp(ii,jj) = zero
          enddo
        enddo

        do ii = 2,n1p1
          do jj = 2,n2
            qy(ii, jj) = zero
            hyp(ii,jj) = zero
          enddo
        enddo

        do jj = 1,n2
          do ii = 1,n1
            pres(ii, jj) = zero
            phi(ii,jj) = zero
          enddo
        enddo
!
! fill bc for velocities defined along the solid
! boundary within the domain
        do jj = 2,n2p1
          qx(1, jj)   = zero
          qx(n1p1,jj) = zero
        enddo
        do ii = 2,n1p1
          qy(ii,1)    = zero
          qy(ii,n2p1) = zero
        enddo
!
! fill bc for velocities defined outside the domain...
        do ii = 1,n1p1
          qx(ii,1)     = -qx(ii,2)                ! Define the no-slip boundary condition for u1 (bottom)
          qx(ii, n2p2) = 2.0 - qx(ii,n2p1)     ! two terms add / 2 = 1.0(Diri BC) for u1 (top)
                                            ! LHS is ghost point
  
        enddo

        do jj = 1,n2p1
          qy(1, jj)    = -qy(2,jj)                ! u2 is all zeros for boundary node
          qy(n1p2, jj) = -qy(n1p1, jj)
        enddo
 
      ELSE
!        read(201) itbeg
!        read(201) tau,dx,dy,dx2,dy2
!        read(201) x,y,xc,yc
!        read(201) qx,qy,hxp,hyp
!        read(201) pres,phi
      ENDIF
!---------------------------------------------------------------------------
! Enter tht time-stepping loop 
! itime  = 1, 2000
!
      itbeg = itbeg + 1
      niter = niter + itbeg-1
      do itime = itbeg, niter

!     if(mod(itime-1,nprof).eq.0)then
!      write(*,*)'itime-1=', itime-1, 'writing'
!      !call output
!     end if

      write(25,*) tau, 0.5*(qx(n1/2+1, n2/2+1) + qx(n1/2+1, n2/2+2)), &
                  0.5*(qy(n1/2+1, n2/2+1) + qy(n1/2+2, n2/2+1))

      tau = tau + dt
!-----------------------------------------------------
! save old qx,qy values to determine if steady-state
! has been reached.
!
      do jj = 1, n2p2
        do ii = 1, n1p1
          qxold(ii, jj) = qx(ii, jj)
        enddo
      enddo

!-----------------------------------------------------
! Back to Adams-Bashforth for the nonlinear term
!         Crank-Nicoholson for the viscous term
!
      irk=1

!-----------------------------------------------------
! calculate non-linear (hx,hy) and viscous (ax,ay) terms
! (sub-step 1)
!
        call hxycalc
        call axycalc
!
! calculate intermediate velocities..
! (sub-step 2)
        call qhat(irk)
!
! calculate the phi from the poisson equation...
! (sub-step 3)
        call trigon(irk)
!
! calculate the final velocities (at current t-step)...
! (sub-step 4)
        call ufinal(irk)
!
! update the pressure (assume wanted)...
! (sub-step 5)
        call press(irk)
!
! determine if steady-state has been reached...
!
      epsmax = 0.0
      sum    = 0.0
      sumtot = 0.0
      do jj = 1, n2p2
        do ii = 1, n1p1
          sum = sum + (qx(ii, jj) - qxold(ii, jj))**2.0
          sumtot = sumtot + (qx(ii, jj))**2.0
        enddo
      enddo
      epsmax = sqrt(sum/sumtot)
! write steady state parameter to file
      if (mod(itime, 50) .eq. 0.0) write(100,*) tau, epsmax

      if (mod(itime, 10) .eq. 0.0) write(*,*) itime, tau, epsmax

      if (itime .gt. 5 .and. epsmax .le. rerror) then
        write(*,*) 'Steady state reached at tau = ', tau
        write(*,*) 'Iteration Number: ', itime
        goto 999
      endif
!
! here is the end of the time-stepping loop...
      enddo
!
  999 continue
      call output
!
! write restart file:
!
      write(202) niter
      write(202) tau, dx, dy, dx2, dy2
      write(202) x, y, xc, yc
      write(202) qx, qy, hxp, hyp
      write(202) pres, phi
end program driver

!----------------------------------------------------------------------
! Variable
! Len                = length,width of cavity
! itbeg              = first time step
! niter              = end time step
! itime              = current time step
! tau                = current time size
!----------------------------------------------------------------------
  

