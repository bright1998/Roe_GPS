c    ------------------
      subroutine write
c    ------------------
      use common
      implicit double precision(a-h,o-z)

      character(len=3) cnm
c      character(len=12) fn
      character(len=13) fn

c      fn='out/'//cnm(nf)//'.data'
c      fn='out/x'//cnm(nf)//'.data'
      fn='out/y'//cnm(nf)//'.data'
      open(99,file=fn,form='formatted')
c      j = (j0 + 1)/2
c      do i=1,i0
c         write(99,990,advance="NO")x(i),ro(i,j),vx(i,j),vy(i,j),vz(i,j),
c     &bx(i,j),by(i,j),bz(i,j),pr(i,j),te(i,j)
c         write(99,*)
cc         write(99,991)x(i),ro(i,j)
c      enddo
      i = (i0 + 1)/2
      do j=1,j0
         write(99,990,advance="NO")y(j),ro(i,j),vx(i,j),vy(i,j),vz(i,j),
     &bx(i,j),by(i,j),bz(i,j),pr(i,j),te(i,j)
         write(99,*)
c         write(99,991)y(j),ro(i,j)
      enddo
      close(99)

c 990  FORMAT(10E24.15)
c 991  FORMAT(2E24.15)
 990  FORMAT(10F24.15)
 991  FORMAT(2F24.15)

      write(6,*)' '
      write(6,*)
     &'wrote data:( nf=',nf,' ns=',ns,' tm=',tm,
     &') in "','*'//cnm(nf)//'.data','"'

      nf = nf + 1

      return
      end subroutine write
