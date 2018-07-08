c    -----------------
      subroutine read
c    -----------------
      use common
      implicit double precision(a-h,o-z)

C Read the data for Restart Calculation
      open(22,file='out/restart.data',form='unformatted')
      read(22)tm,ns,nf,nbcast
      read(22)
     &((ro(i,j),i=1,i0),j=1,j0),
     &((vx(i,j),i=1,i0),j=1,j0),
     &((vy(i,j),i=1,i0),j=1,j0),
     &((vz(i,j),i=1,i0),j=1,j0),
     &((bx(i,j),i=1,i0),j=1,j0),
     &((by(i,j),i=1,i0),j=1,j0),
     &((bz(i,j),i=1,i0),j=1,j0),
     &((pr(i,j),i=1,i0),j=1,j0),
     &((te(i,j),i=1,i0),j=1,j0)
      close(22)

      write(6,*)' '
      write(6,*)
     &'read data for restart calculation'

C Calculate conserved variable
      do j=1,j0
      do i=1,i0
         w(i,j,1) = ro(i,j)
         w(i,j,2) = ro(i,j)*vx(i,j)
         w(i,j,3) = ro(i,j)*vy(i,j)
         w(i,j,4) = ro(i,j)*vz(i,j)
         w(i,j,5) = bx(i,j)
         w(i,j,6) = by(i,j)
         w(i,j,7) = bz(i,j)
         w(i,j,8) = 0.5d0*ro(i,j)*(vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2)
     &            + pr(i,j)/(gm - 1.d0)
     &            + 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
      enddo
      enddo

      return
      end subroutine read
