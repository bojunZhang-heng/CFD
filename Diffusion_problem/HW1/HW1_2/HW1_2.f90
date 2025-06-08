program ieee754
  implicit none
  real*4    :: x
  integer   :: bx
  do
    read(*,*) x
    write(*,*) 'Output = ', x
    bx = transfer(x, bx)
    write(*,'(B32.32)') bx
  end do
end program ieee754

