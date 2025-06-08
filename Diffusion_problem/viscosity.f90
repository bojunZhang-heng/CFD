Program Viscosity
   real*4 mu1, mu2, b, S, mu0, T, TK, T0
   integer i
       
! Compare two models of viscosity

   b = 1.458e-6
   S = 110.4

   mu0 = 1.71e-5
   T0 = 273.

       do i = -5,10
       T = 10.0*real(i)   ! in C
! Convert to K
       TK = T + 273.
       mu1 = b*TK**1.5/(TK+S)
       mu2 = mu0*(TK/T0)**0.7
       write(10,'(2F7.1,2(1PE15.5))') T,TK,mu1,mu2
! Print out a data table
       end do

End Program Viscosity

