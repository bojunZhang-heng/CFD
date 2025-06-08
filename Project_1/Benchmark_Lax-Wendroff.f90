program lax_wendroff_advection
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: Nx = 80
    real(dp), parameter :: xmin = 0.0_dp, xmax = 2.0_dp
    real(dp), parameter :: tmin = 0.0_dp, tmax = 1.0_dp
    real(dp), parameter :: a = 1.0_dp, c = 0.9_dp
    real(dp) :: dx, dt
    integer :: Nt, i, n

    real(dp), dimension(0:Nx) :: x, u, u_new, u_initial
    real(dp), dimension(:,:), allocatable :: u_all
    real(dp) :: x_left, x_right

    ! Spatial and temporal discretization
    dx = (xmax - xmin) / real(Nx, dp)
    dt = c * dx / a
    Nt = int((tmax - tmin) / dt)

    allocate(u_all(0:Nx, 0:Nt))

    ! Initialize spatial grid
    do i = 0, Nx
        x(i) = xmin + i * dx
    end do

    ! Select sin^4 initial condition (init_func = 1)
    x_left = 0.25_dp
    x_right = 0.75_dp
    do i = 0, Nx
        if (x(i) > x_left .and. x(i) < x_right) then
            u(i) = sin(acos(-1.0_dp) * (x(i) - x_left) / (x_right - x_left))**4
        else
            u(i) = 0.0_dp
        end if
    end do
    u_initial = u
    u_all(:, 0) = u

    ! Time stepping using Lax-Wendroff
    do n = 1, Nt
        ! Periodic boundary
        u_new(0)   = c/2.0_dp*(1.0_dp + c)*u(Nx-1) + (1.0_dp - c**2)*u(0) - c/2.0_dp*(1.0_dp - c)*u(1)
        u_new(Nx)  = c/2.0_dp*(1.0_dp + c)*u(Nx-1) + (1.0_dp - c**2)*u(Nx) - c/2.0_dp*(1.0_dp - c)*u(1)

        do i = 1, Nx-1
            u_new(i) = c/2.0_dp*(1.0_dp + c)*u(i-1) + (1.0_dp - c**2)*u(i) - c/2.0_dp*(1.0_dp - c)*u(i+1)
        end do

        u = u_new
        u_all(:, n) = u
    end do

    ! Output results to file
    open(unit=10, file="lax_wendroff_output.dat", status="replace")
    do n = 0, 1
        write(10, '(100(2x, F16.8))') (u_all(i, n), i=0, Nx)
    end do
    close(10)

    print *, 'Simulation complete. Output written to lax_wendroff_output.dat'

end program lax_wendroff_advection

