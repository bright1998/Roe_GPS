c    -----------------
      subroutine grid
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      ddx = (xcmax - xcmin)/(ic - 1)
      ddy = (ycmax - ycmin)/(jc - 1)

      if(ic .ne. i0) then
         i1 = (i0 - ic)/2
         i2 = i1 + ic

         do i=i1+1,i2-1
            dx(i) = ddx
         enddo
         do i=i1,1,-1
            dx(i) = dx(i+1)*rax
         enddo
         do i=i2,i0
            dx(i) = dx(i-1)*rax
         enddo

         x(i1+1) = xcmin
         do i=i1+2,i0
            x(i) = x(i-1) + dx(i-1)
         enddo
         do i=i1,1,-1
            x(i) = x(i+1) - dx(i)
         enddo
      else
         do i=1,ic
            dx(i) = ddx
         enddo
         x(1) = xcmin
         do i=2,ic
            x(i) = x(i-1) + dx(i-1)
         enddo
      endif

      if(jc .ne. j0) then
         j1 = (j0 - jc)/2
         j2 = j1 + jc

         do j=j1+1,j2-1
            dy(j) = ddy
         enddo
         do j=j1,1,-1
            dy(j) = dy(j+1)*ray
         enddo
         do j=j2,j0
            dy(j) = dy(j-1)*ray
         enddo

         y(j1+1) = ycmin
         do j=j1+2,j0
            y(j) = y(j-1) + dy(j-1)
         enddo
         do j=j1,1,-1
            y(j) = y(j+1) - dy(j)
         enddo
      else
         do j=1,jc
            dy(j) = ddy
         enddo
         y(1) = ycmin
         do j=2,jc
            y(j) = y(j-1) + dy(j-1)
         enddo
      endif

      do i=1,i0
         dxi(i) = 1.d0/dx(i)
      enddo

      do i=1,i0-1
         xm(i) = 0.5d0*(x(i) + x(i+1))
      enddo
      do i=1,i0-1
         dxm(i) = 0.5d0*(dx(i) + dx(i+1))
      enddo
      do i=1,i0-1
         dxmi(i) = 1.d0/dxm(i)
      enddo

      do j=1,j0
         dyi(j) = 1.d0/dy(j)
      enddo

      do j=1,j0-1
         ym(j) = 0.5d0*(y(j) + y(j+1))
      enddo
      do j=1,j0-1
         dym(j) = 0.5d0*(dy(j) + dy(j+1))
      enddo
      do j=1,j0-1
         dymi(j) = 1.d0/dym(j)
      enddo

      open(10,file='out/grid.data')
      write(10,*)i0,j0
      close(10)

      open(11,file='out/x-grid.data',form='unformatted')
      write(11)(x(i),i=1,i0)
      close(11)

      open(12,file='out/y-grid.data',form='unformatted')
      write(12)(y(j),j=1,j0)
      close(12)

      return
      end subroutine grid
