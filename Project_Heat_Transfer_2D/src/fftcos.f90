module fftcos
  use mod_param, only: n1, n2, n1p1, n1p2, n1h
  use mod_fft_data, only: trigsx, ifaxx

contains

subroutine F2RC1(ersatz)
  PARAMETER ( INCX=1,JUMPX=n1p2 )
  real, intent(inout) :: ersatz(n1, n2)
  
  real :: tempw(n1p2, n2)
  real :: workx(n2*n1p1)
  double precision :: theta, wi, wi1, wpi, wpr, wr, wr1, wtemp, pi
  integer :: i, j, ia
  real :: y1, y2

  pi = 4.0d0 * atan(1.0d0)
  theta = 0.5d0 * pi / n1
  wr = 1.0d0
  wi = 0.0d0
  wr1 = cos(theta)
  wi1 = sin(theta)
  wpr = -2.0d0 * wi1**2
  wpi = sin(2.d0 * theta)

  do j = 1, n2
    tempw(1, j) = ersatz(1, j)
    tempw(2, j) = 0.0
    do i = n1, 4, -2
      tempw(i, j) = ersatz(i-2, j) - ersatz(i, j)
    end do
    tempw(n1p1, j) = ersatz(n1, j)
    tempw(n1p2, j) = 0.0
  end do

  do i = 3, n1, 2
    wtemp = wr
    wr = wr * wpr - wi * wpi + wr
    wi = wi * wpr + wtemp * wpi + wi

    do j = 1, n2
      y1 = ersatz(i, j) * wr + tempw(i+1, j) * wi
      y2 = tempw(i+1, j) * wr - ersatz(i, j) * wi
      tempw(i, j) = y1 / 2.0
      tempw(i+1, j) = -y2 / 2.0
    end do
  end do

  call FFT991(tempw(1,1), workx, trigsx, ifaxx, 1, n1p2, n1, n2, +1)

  do i = 1, n1h
    ia = n1p1 - i
    do j = 1, n2
      y1 = tempw(i, j) + tempw(ia, j)
      y2 = 0.5 / wi1 * (tempw(i, j) - tempw(ia, j))
      ersatz(i, j) = 0.5 * (y1 + y2)
      ersatz(ia, j) = 0.5 * (y1 - y2)
    end do

    wtemp = wr1
    wr1 = wr1 * wpr - wi1 * wpi + wr1
    wi1 = wi1 * wpr + wtemp * wpi + wi1
  end do

end subroutine F2RC1

!-----------------------------------------------------------------------

subroutine R2FC1(ersatz)
  real, intent(inout) :: ersatz(n1, n2)
  
  real :: tempw(n1p2, n2)
  real :: workx(n2*n1p1)
  double precision :: theta, wi, wi1, wpi, wpr, wr, wr1, wtemp, pi
  integer :: i, j, ia
  real :: y1, y2
  real :: sum, sum1

  pi = 4.0d0 * atan(1.0d0)
  theta = 0.5d0 * pi / n1
  wr = 1.0d0
  wi = 0.0d0
  wr1 = cos(theta)
  wi1 = sin(theta)
  wpr = -2.0d0 * wi1**2
  wpi = sin(2.d0 * theta)

  do i = 1, n1h
    ia = n1p1 - i
    do j = 1, n2
      y1 = 0.5 * (ersatz(i, j) + ersatz(ia, j))
      y2 = wi1 * (ersatz(i, j) - ersatz(ia, j))
      tempw(i, j) = y1 + y2
      tempw(ia, j) = y1 - y2
    end do

    wtemp = wr1
    wr1 = wr1 * wpr - wi1 * wpi + wr1
    wi1 = wi1 * wpr + wtemp * wpi + wi1
  end do

  do j = 1, n2
    tempw(n1p1, j) = 0.0
    tempw(n1p2, j) = 0.0
  end do

  call FFT991(tempw(1,1), workx, trigsx, ifaxx, 1, n1p2, n1, n2, -1)

  ! - even modes -
  do j = 1, n2
    ersatz(1, j) = tempw(1, j)
  end do

  do i = 3, n1, 2
    wtemp = wr
    wr = wr * wpr - wi * wpi + wr
    wi = wi * wpr + wtemp * wpi + wi

    do j = 1, n2
      y1 = tempw(i, j) * wr + tempw(i+1, j) * wi
      y2 = -tempw(i+1, j) * wr + tempw(i, j) * wi
      ersatz(i, j) = 2.0 * y1
      ersatz(i+1, j) = y2
    end do
  end do

  ! - odd modes -
  do j = 1, n2
    sum = 0.5 * tempw(n1p1, j)
    do i = n1, 2, -2
      sum1 = sum
      sum = sum + ersatz(i, j)
      ersatz(i, j) = 2.0 * sum1
    end do
  end do

end subroutine R2FC1

end module fftcos
