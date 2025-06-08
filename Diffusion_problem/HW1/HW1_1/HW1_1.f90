Program main
  implicit none

  integer, parameter       :: N = 400        ! Number of toss times
  integer                  :: M              ! Number of having the head times 
  real*8                   :: P              ! The probability of having the head M times

  open(unit=10, file='HW1_1_result.txt')

  do M = 0, N
    P = exp(log_binomial_coeff(N,M) - N*log(2.0d0))
    if(P < 1.0d-6) then
      P = 0.0d0
    endif
  write(10,*) M, P
  enddo

  close(10)

  print *, "Computation complete. Data sabed om HW1_1_result.txt"
  print *, "Generating plot uing Gnuplot..."

  ! Create Gnuplot script file
  open(unit=20, file='HW1_1_plot.gp');

  write(20,*) "Set x-label 'Number of Heads (M)'"
  write(20,*) "Set y-label 'Probability (P)'"
  write(20,*) "Set title 'Binomial Probability Distribution (N=400)'"
  write(20,*) "set grid"
  write(20,*) "plot 'HW1_1_result.txt' using 1:2 with linespoints title 'Probability'" 
  close(20)

  call system("gnuplot -persist HW1_1_plot.gp")


contains
  ! compute log binomail coefficient using Stirling's approximation
  function log_binomial_coeff(N,M) result(log_bc)
    integer, intent(in)  :: N, M
    real*8               :: log_bc

    if (M == 0 .or. M==N ) then
      log_bc = 1.0d0
      return
    endif

    ! binomial_coeff(N,M) = n/1 * (n-1)/2 * (n-2)/3 * ... * (N-M+1)/M
    log_bc = stirling_log_factorial(N) - stirling_log_factorial(M) - stirling_log_factorial(N-M)
  end function log_binomial_coeff

  ! Stirling's Approximation for log(n!)
  function stirling_log_factorial(N) result(log_fact)
    integer, intent(in)                :: N
    real(8)                            :: log_fact
 
    if (N == 0) then
      log_fact = 0.0d0
    else
      log_fact = N * log(real(N,8)) - N + 0.5d0 * log(2.0d0 * 3.141592653589793d0 * real(N,8))
    endif
  end function stirling_log_factorial
end Program main




