subroutine output

      include 'param.inc'
      common/one/qx(n1p1, n2p2), qy(n1p2, n2p1), &
                 dx, dy, dx2, dy2, &
                 dt, re, pr, aaa(3), ggg(3), rrr(3)

      common/two/hx(2:n1, 2:n2p1),  hy(2:n1p1, 2:n2), &
                 hxp(2:n1, 2:n2p1), hyp(2:n1p1, 2:n2)

      common/three/ax(2:n1, 2:n2p1), ay(2:n1p1, 2:n2)

      common/six/pres(n1, n2),   phi(n1, n2), &
                 phip(n1, n2), presp(n1, n2)

      common/seven/x(n1p1), y(n2p1), xc(n1), yc(n2)
      dimension alam(n1, n2)
      dimension uvel(n1p1, n2p1), vvel(n1p1, n2p1), strmfn(n1p1, n2p1), &
                vortfn(n1p1, n2p1), temp(n1p1)

      do jj = 1, n2p1
          u9 = 0.5*(qx(n1/2+1, jj) + qx(n1/2+1, jj+1))
          v9 = 0.5*(qy(n1/2+1, jj) + qy(n1/2+2, jj))
! velocity at the vertical line
        write(10,102) y(jj),u9,v9
      enddo

      do ii = 1, n1p1
          u9 = 0.5*(qx(ii, n2/2+1) + qx(ii, n2/2+2))
          v9 = 0.5*(qy(ii, n2/2+1) + qy(ii+1, n2/2+1))
! velocity at the vertical line
        write(11,102) x(i),u9,v9
      enddo

      write(10,101)
      write(11,101)
101   format(1x)
102   format(2x,3(1pe15.4))
!
! save values of velocity at mesh vertices...
!
      do ii = 1,n1p1
        uvel(ii, 1)    = 0.5*(qx(ii, 1) + qx(ii,2))
        uvel(ii, n2p1) = 0.5*(qx(ii, n2p2) + qx(ii,n2p1))
        do jj = 2,n2
          uvel(ii, jj) = 0.5*(qx(ii, jj+1) + qx(ii, jj))
        enddo
      enddo

      do jj = 1,n2p1
        vvel(1, jj)   = 0.5*(qy(1, jj) + qy(2, jj))
        vvel(n1p1,jj) = 0.5*(qy(n1p2, jj) + qy(n1p1, jj))
        do ii = 2, n1
          vvel(ii, jj) = 0.5*(qy(ii+1, jj) + qy(ii, jj))
        enddo
      enddo


      do jj = 1, n2
        do ii = 1, n1
          dux = (qx(ii+1, jj+1) - qx(ii, jj+1)) / dx
          dvy = (qy(ii+1, jj+1) - qy(ii+1, jj)) / dx
          duy = ((qx(ii, jj+2)  - qx(ii, jj)) / 2.d0 / dy &
              + (qx(ii+1, jj+2) - qx(ii+1, jj)) /2.d0 / dy) / 2.0
          dvy = ((qy(ii+2, jj)  - qy(ii, jj)) / 2.d0 / dx &
              + (qy(ii+2, jj+1) - qy(ii, jj+1)) / 2.d0 / dx) / 2.0
          alam(ii, jj) = 0.0
          s9 = dux*dvy - duy*dvx
        if(s9 .gt. 0.0) alam(ii, jj) = sqrt(s9)
        enddo
      enddo
!
! calculate streamfunction and vorticity
!
      do ii = 1, n1p1
        temp1 = 0.d0
        strmfn(ii, 1) = 0.d0
        do jj = 2, n2p1
          temp1 = temp1 + qx(ii, jj)*dy
          strmfn(ii, jj) = temp1
        enddo
      enddo
!
      do jj = 1, n2p1
        do ii = 1, n1p1
          vortfn(ii, jj) = (qy(ii+1,jj)  - qy(ii, jj)) / dx  &
                         - (qx(ii, jj+1) - qx(ii, jj)) / dy
        enddo
      enddo
!
!      open(11,file='output.dat',form='formatted',
!     >      status='unknown')
!* this is Tecplot format...
!      write(11,'(''TITLE = "2D Cavity"'')')
!      write(11,'(''VARIABLES = X, Y, U, V, VORT, STRM'')')
!*23456789012345678901234567890123456789012345678901234567890123456789012
!      write(11,'(''ZONE T="Zone-one", I='',i4,''  J='',i4,''  F=BLOCK'')
!     >    ') n1p1,n2p1
!      do j=1,n2p1
!        write(11,905) (x(i),i=1,n1p1)
!      enddo
!      do j=1,n2p1
!        do i=1,n1p1
!          temp(i) = y(j)
!        enddo
!        write(11,905) (temp(i),i=1,n1p1)
!      enddo
!      write(11,905) ((uvel(i,j),i=1,n1p1),j=1,n2p1)
!      write(11,905) ((vvel(i,j),i=1,n1p1),j=1,n2p1)
!      write(11,905) ((vortfn(i,j),i=1,n1p1),j=1,n2p1)
!      write(11,905) ((strmfn(i,j),i=1,n1p1),j=1,n2p1)
      write(19,908) ((vortfn(ii, jj), ii = 1, n1p1), jj = 1, n2p1)
      write(18,907) ((strmfn(ii, jj), ii = 1, n1p1), jj = 1, n2p1)
!
      open(13,file='strmfn.dat',form='formatted',status='unknown')
! location of primary vortex center is the region of mimimum streamfunction
      amin = 0.d0
      do jj = 1, n2p1
        do ii = 1, n1p1
          if (strmfn(ii, jj) .lt. amin)then
            amin = strmfn(ii, jj)
            imin = ii
            jmin = jj
          endif
! Rescale to match the LBM code
        write(13, *) strmfn(ii, jj) * 96.d0 / 10.d0
        enddo
!       write(13,906)
      enddo
!
      open(12,file='vortfn.dat',form='formatted',status='unknown')
      do jj = 1, n2p1
        do ii = 1, n1p1
          write(12,*) vortfn(ii, jj)
        enddo
        write(12,906)
      enddo
      write(*,*)'streamfuncton and vorticity at the center', &
         amin,vortfn(imin,jmin), imin, jmin
      close(12)

!
      open(12,file='alam.dat',form='formatted',status='unknown')
      do jj = 1, n2
        do ii = 1, n1
          write(12,*) alam(ii, jj)
        enddo
      enddo
      close(12)
!
      open(14,file='uvel.dat',form='formatted',status='unknown')
      do jj = 1,n2p1
        do ii = 1,n1p1
          write(14,*) x(ii), y(jj), uvel(ii, jj)
        enddo
        write(14,906)
      enddo
!
      open(15,file='vvel.dat',form='formatted',status='unknown')
      do jj = 1, n2p1
        do ii = 1, n1p1
          write(15,*) x(ii), y(jj), vvel(ii, jj)
        enddo
!       write(15,906)
      enddo
!
! added 10/8/96
      open(16,file='pres.dat',form='formatted',status='unknown')
      do jj = 1, n2
        do ii = 1, n1
          write(16,*) xc(ii), yc(jj), pres(ii, jj)
        enddo
        write(16,906)
      enddo
!
  905 format(6(e12.4))
  906 format(x)
  907 format(2x,8f8.4)
  908 format(2x,8f8.2)
end subroutine
