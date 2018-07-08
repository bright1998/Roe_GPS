c    --------------------
      subroutine restart
c    --------------------
      use common
      implicit double precision(a-h,o-z)

C Output the data for Restart Calculation
      open(22,file='out/restart.data',form='unformatted')
      write(22)tm,ns,nf,nbcast
      write(22)
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
     &'wrote data for restart calculation'

      return
      end subroutine restart
