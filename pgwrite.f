c    --------------------
      subroutine pgwrite
c    --------------------
      use common
      implicit double precision(a-h,o-z)

      character(len=3) cnm

C Output the data for PGPLOT
      open(21,file='out/pgplot'//cnm(nf)//'.data',form='unformatted')
      write(21)tm,ns,nf,nbcast
      write(21)
     &((ro(i,j),i=1,i0),j=1,j0),
     &((vx(i,j),i=1,i0),j=1,j0),
     &((vy(i,j),i=1,i0),j=1,j0),
     &((vz(i,j),i=1,i0),j=1,j0),
     &((bx(i,j),i=1,i0),j=1,j0),
     &((by(i,j),i=1,i0),j=1,j0),
     &((bz(i,j),i=1,i0),j=1,j0),
     &((pr(i,j),i=1,i0),j=1,j0),
     &((te(i,j),i=1,i0),j=1,j0)
      close(21)

      write(6,*)' '
      write(6,*)
     &'wrote data for pgplot:( nf=',nf,' ns=',ns,' tm=',tm,
     &') in "','pgplot'//cnm(nf)//'.data','"'

      return
      end subroutine pgwrite
