c    -----------------
      subroutine bndc
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine bdfrdx(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdfrdx

         subroutine bdfrex(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdfrex

         subroutine bdperx(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdperx

         subroutine bdsymx(da,margin,mbnd,men,coff)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: coff
            integer :: margin,mbnd,men
         end subroutine bdsymx

         subroutine bdfrdy(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdfrdy

         subroutine bdfrey(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdfrey

         subroutine bdpery(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdpery

         subroutine bdsymy(da,margin,mbnd,men,coff)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: coff
            integer :: margin,mbnd,men
         end subroutine bdsymy
      end interface

      if(nbnd1 .eq. 1) then
         call bdfrdx(ro,ro0,3,0,0)
         call bdfrdx(vx,vx0,3,0,0)
         call bdfrdx(vy,vy0,3,0,0)
         call bdfrdx(vz,vz0,3,0,0)
         call bdfrdx(bx,bx0,3,0,0)
         call bdfrdx(by,by0,3,0,0)
         call bdfrdx(bz,bz0,3,0,0)
         call bdfrdx(pr,pr0,3,0,0)
      elseif(nbnd1 .eq. 2) then
         call bdfrex(ro,3,0,0)
         call bdfrex(vx,3,0,0)
         call bdfrex(vy,3,0,0)
         call bdfrex(vz,3,0,0)
         call bdfrex(bx,3,0,0)
         call bdfrex(by,3,0,0)
         call bdfrex(bz,3,0,0)
         call bdfrex(pr,3,0,0)
      elseif(nbnd1 .eq. 3) then
         call bdperx(ro,3,0,0)
         call bdperx(vx,3,0,0)
         call bdperx(vy,3,0,0)
         call bdperx(vz,3,0,0)
         call bdperx(bx,3,0,0)
         call bdperx(by,3,0,0)
         call bdperx(bz,3,0,0)
         call bdperx(pr,3,0,0)
      elseif(nbnd1 .eq. 4) then
         call bdsymx(ro,3,0,0, 1.d0)
         call bdsymx(vx,3,0,0,-1.d0)
         call bdsymx(vy,3,0,0, 1.d0)
         call bdsymx(vz,3,0,0, 1.d0)
         call bdsymx(bx,3,0,0, 1.d0)
         call bdsymx(by,3,0,0,-1.d0)
         call bdsymx(bz,3,0,0,-1.d0)
         call bdsymx(pr,3,0,0, 1.d0)         
      endif

      if(nbnd2 .eq. 1) then
         call bdfrdx(ro,ro0,3,1,0)
         call bdfrdx(vx,vx0,3,1,0)
         call bdfrdx(vy,vy0,3,1,0)
         call bdfrdx(vz,vz0,3,1,0)
         call bdfrdx(bx,bx0,3,1,0)
         call bdfrdx(by,by0,3,1,0)
         call bdfrdx(bz,bz0,3,1,0)
         call bdfrdx(pr,pr0,3,1,0)
      elseif(nbnd2 .eq. 2) then
         call bdfrex(ro,3,1,0)
         call bdfrex(vx,3,1,0)
         call bdfrex(vy,3,1,0)
         call bdfrex(vz,3,1,0)
         call bdfrex(bx,3,1,0)
         call bdfrex(by,3,1,0)
         call bdfrex(bz,3,1,0)
         call bdfrex(pr,3,1,0)
      elseif(nbnd2 .eq. 3) then
         call bdperx(ro,3,1,0)
         call bdperx(vx,3,1,0)
         call bdperx(vy,3,1,0)
         call bdperx(vz,3,1,0)
         call bdperx(bx,3,1,0)
         call bdperx(by,3,1,0)
         call bdperx(bz,3,1,0)
         call bdperx(pr,3,1,0)
      elseif(nbnd2 .eq. 4) then
         call bdsymx(ro,3,1,0, 1.d0)
         call bdsymx(vx,3,1,0,-1.d0)
         call bdsymx(vy,3,1,0, 1.d0)
         call bdsymx(vz,3,1,0, 1.d0)
         call bdsymx(bx,3,1,0, 1.d0)
         call bdsymx(by,3,1,0,-1.d0)
         call bdsymx(bz,3,1,0,-1.d0)
         call bdsymx(pr,3,1,0, 1.d0)         
      endif

      if(nbnd3 .eq. 1) then
         call bdfrdy(ro,ro0,3,0,0)
         call bdfrdy(vx,vx0,3,0,0)
         call bdfrdy(vy,vy0,3,0,0)
         call bdfrdy(vz,vz0,3,0,0)
         call bdfrdy(bx,bx0,3,0,0)
         call bdfrdy(by,by0,3,0,0)
         call bdfrdy(bz,bz0,3,0,0)
         call bdfrdy(pr,pr0,3,0,0)
      elseif(nbnd3 .eq. 2) then
         call bdfrey(ro,3,0,0)
         call bdfrey(vx,3,0,0)
         call bdfrey(vy,3,0,0)
         call bdfrey(vz,3,0,0)
         call bdfrey(bx,3,0,0)
         call bdfrey(by,3,0,0)
         call bdfrey(bz,3,0,0)
         call bdfrey(pr,3,0,0)
      elseif(nbnd3 .eq. 3) then
         call bdpery(ro,3,0,0)
         call bdpery(vx,3,0,0)
         call bdpery(vy,3,0,0)
         call bdpery(vz,3,0,0)
         call bdpery(bx,3,0,0)
         call bdpery(by,3,0,0)
         call bdpery(bz,3,0,0)
         call bdpery(pr,3,0,0)
      elseif(nbnd3 .eq. 4) then
         call bdsymy(ro,3,0,0, 1.d0)
         call bdsymy(vx,3,0,0, 1.d0)
         call bdsymy(vy,3,0,0,-1.d0)
         call bdsymy(vz,3,0,0, 1.d0)
         call bdsymy(bx,3,0,0,-1.d0)
         call bdsymy(by,3,0,0, 1.d0)
         call bdsymy(bz,3,0,0,-1.d0)
         call bdsymy(pr,3,0,0, 1.d0)
      endif

      if(nbnd4 .eq. 1) then
         call bdfrdy(ro,ro0,3,1,0)
         call bdfrdy(vx,vx0,3,1,0)
         call bdfrdy(vy,vy0,3,1,0)
         call bdfrdy(vz,vz0,3,1,0)
         call bdfrdy(bx,bx0,3,1,0)
         call bdfrdy(by,by0,3,1,0)
         call bdfrdy(bz,bz0,3,1,0)
         call bdfrdy(pr,pr0,3,1,0)
      elseif(nbnd4 .eq. 2) then
         call bdfrey(ro,3,1,0)
         call bdfrey(vx,3,1,0)
         call bdfrey(vy,3,1,0)
         call bdfrey(vz,3,1,0)
         call bdfrey(bx,3,1,0)
         call bdfrey(by,3,1,0)
         call bdfrey(bz,3,1,0)
         call bdfrey(pr,3,1,0)
      elseif(nbnd4 .eq. 3) then
         call bdpery(ro,3,1,0)
         call bdpery(vx,3,1,0)
         call bdpery(vy,3,1,0)
         call bdpery(vz,3,1,0)
         call bdpery(bx,3,1,0)
         call bdpery(by,3,1,0)
         call bdpery(bz,3,1,0)
         call bdpery(pr,3,1,0)
      elseif(nbnd4 .eq. 4) then
         call bdsymy(ro,3,1,0, 1.d0)
         call bdsymy(vx,3,1,0, 1.d0)
         call bdsymy(vy,3,1,0,-1.d0)
         call bdsymy(vz,3,1,0, 1.d0)
         call bdsymy(bx,3,1,0,-1.d0)
         call bdsymy(by,3,1,0, 1.d0)
         call bdsymy(bz,3,1,0,-1.d0)
         call bdsymy(pr,3,1,0, 1.d0)
      endif

      return
      end subroutine bndc

c    ------------------
      subroutine bndcx
c    ------------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine bdfrdx(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdfrdx

         subroutine bdfrex(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdfrex

         subroutine bdperx(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdperx

         subroutine bdsymx(da,margin,mbnd,men,coff)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: coff
            integer :: margin,mbnd,men
         end subroutine bdsymx
      end interface

      if(nbnd1 .eq. 1) then
         call bdfrdx(ro,ro0,3,0,0)
         call bdfrdx(vx,vx0,3,0,0)
         call bdfrdx(vy,vy0,3,0,0)
         call bdfrdx(vz,vz0,3,0,0)
         call bdfrdx(bx,bx0,3,0,0)
         call bdfrdx(by,by0,3,0,0)
         call bdfrdx(bz,bz0,3,0,0)
         call bdfrdx(pr,pr0,3,0,0)
      elseif(nbnd1 .eq. 2) then
         call bdfrex(ro,3,0,0)
         call bdfrex(vx,3,0,0)
         call bdfrex(vy,3,0,0)
         call bdfrex(vz,3,0,0)
         call bdfrex(bx,3,0,0)
         call bdfrex(by,3,0,0)
         call bdfrex(bz,3,0,0)
         call bdfrex(pr,3,0,0)
      elseif(nbnd1 .eq. 3) then
         call bdperx(ro,3,0,0)
         call bdperx(vx,3,0,0)
         call bdperx(vy,3,0,0)
         call bdperx(vz,3,0,0)
         call bdperx(bx,3,0,0)
         call bdperx(by,3,0,0)
         call bdperx(bz,3,0,0)
         call bdperx(pr,3,0,0)
      elseif(nbnd1 .eq. 4) then
         call bdsymx(ro,3,0,0, 1.d0)
         call bdsymx(vx,3,0,0,-1.d0)
         call bdsymx(vy,3,0,0, 1.d0)
         call bdsymx(vz,3,0,0, 1.d0)
         call bdsymx(bx,3,0,0, 1.d0)
         call bdsymx(by,3,0,0,-1.d0)
         call bdsymx(bz,3,0,0,-1.d0)
         call bdsymx(pr,3,0,0, 1.d0)         
      endif

      if(nbnd2 .eq. 1) then
         call bdfrdx(ro,ro0,3,1,0)
         call bdfrdx(vx,vx0,3,1,0)
         call bdfrdx(vy,vy0,3,1,0)
         call bdfrdx(vz,vz0,3,1,0)
         call bdfrdx(bx,bx0,3,1,0)
         call bdfrdx(by,by0,3,1,0)
         call bdfrdx(bz,bz0,3,1,0)
         call bdfrdx(pr,pr0,3,1,0)
      elseif(nbnd2 .eq. 2) then
         call bdfrex(ro,3,1,0)
         call bdfrex(vx,3,1,0)
         call bdfrex(vy,3,1,0)
         call bdfrex(vz,3,1,0)
         call bdfrex(bx,3,1,0)
         call bdfrex(by,3,1,0)
         call bdfrex(bz,3,1,0)
         call bdfrex(pr,3,1,0)
      elseif(nbnd2 .eq. 3) then
         call bdperx(ro,3,1,0)
         call bdperx(vx,3,1,0)
         call bdperx(vy,3,1,0)
         call bdperx(vz,3,1,0)
         call bdperx(bx,3,1,0)
         call bdperx(by,3,1,0)
         call bdperx(bz,3,1,0)
         call bdperx(pr,3,1,0)
      elseif(nbnd2 .eq. 4) then
         call bdsymx(ro,3,1,0, 1.d0)
         call bdsymx(vx,3,1,0,-1.d0)
         call bdsymx(vy,3,1,0, 1.d0)
         call bdsymx(vz,3,1,0, 1.d0)
         call bdsymx(bx,3,1,0, 1.d0)
         call bdsymx(by,3,1,0,-1.d0)
         call bdsymx(bz,3,1,0,-1.d0)
         call bdsymx(pr,3,1,0, 1.d0)         
      endif

      return
      end subroutine bndcx

c    ------------------
      subroutine bndcy
c    ------------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine bdfrdy(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdfrdy

         subroutine bdfrey(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdfrey

         subroutine bdpery(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdpery

         subroutine bdsymy(da,margin,mbnd,men,coff)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: coff
            integer :: margin,mbnd,men
         end subroutine bdsymy
      end interface

      if(nbnd3 .eq. 1) then
         call bdfrdy(ro,ro0,3,0,0)
         call bdfrdy(vx,vx0,3,0,0)
         call bdfrdy(vy,vy0,3,0,0)
         call bdfrdy(vz,vz0,3,0,0)
         call bdfrdy(bx,bx0,3,0,0)
         call bdfrdy(by,by0,3,0,0)
         call bdfrdy(bz,bz0,3,0,0)
         call bdfrdy(pr,pr0,3,0,0)
      elseif(nbnd3 .eq. 2) then
         call bdfrey(ro,3,0,0)
         call bdfrey(vx,3,0,0)
         call bdfrey(vy,3,0,0)
         call bdfrey(vz,3,0,0)
         call bdfrey(bx,3,0,0)
         call bdfrey(by,3,0,0)
         call bdfrey(bz,3,0,0)
         call bdfrey(pr,3,0,0)
      elseif(nbnd3 .eq. 3) then
         call bdpery(ro,3,0,0)
         call bdpery(vx,3,0,0)
         call bdpery(vy,3,0,0)
         call bdpery(vz,3,0,0)
         call bdpery(bx,3,0,0)
         call bdpery(by,3,0,0)
         call bdpery(bz,3,0,0)
         call bdpery(pr,3,0,0)
      elseif(nbnd3 .eq. 4) then
         call bdsymy(ro,3,0,0, 1.d0)
         call bdsymy(vx,3,0,0, 1.d0)
         call bdsymy(vy,3,0,0,-1.d0)
         call bdsymy(vz,3,0,0, 1.d0)
         call bdsymy(bx,3,0,0,-1.d0)
         call bdsymy(by,3,0,0, 1.d0)
         call bdsymy(bz,3,0,0,-1.d0)
         call bdsymy(pr,3,0,0, 1.d0)
      endif

      if(nbnd4 .eq. 1) then
         call bdfrdy(ro,ro0,3,1,0)
         call bdfrdy(vx,vx0,3,1,0)
         call bdfrdy(vy,vy0,3,1,0)
         call bdfrdy(vz,vz0,3,1,0)
         call bdfrdy(bx,bx0,3,1,0)
         call bdfrdy(by,by0,3,1,0)
         call bdfrdy(bz,bz0,3,1,0)
         call bdfrdy(pr,pr0,3,1,0)
      elseif(nbnd4 .eq. 2) then
         call bdfrey(ro,3,1,0)
         call bdfrey(vx,3,1,0)
         call bdfrey(vy,3,1,0)
         call bdfrey(vz,3,1,0)
         call bdfrey(bx,3,1,0)
         call bdfrey(by,3,1,0)
         call bdfrey(bz,3,1,0)
         call bdfrey(pr,3,1,0)
      elseif(nbnd4 .eq. 3) then
         call bdpery(ro,3,1,0)
         call bdpery(vx,3,1,0)
         call bdpery(vy,3,1,0)
         call bdpery(vz,3,1,0)
         call bdpery(bx,3,1,0)
         call bdpery(by,3,1,0)
         call bdpery(bz,3,1,0)
         call bdpery(pr,3,1,0)
      elseif(nbnd4 .eq. 4) then
         call bdsymy(ro,3,1,0, 1.d0)
         call bdsymy(vx,3,1,0, 1.d0)
         call bdsymy(vy,3,1,0,-1.d0)
         call bdsymy(vz,3,1,0, 1.d0)
         call bdsymy(bx,3,1,0,-1.d0)
         call bdsymy(by,3,1,0, 1.d0)
         call bdsymy(bz,3,1,0,-1.d0)
         call bdsymy(pr,3,1,0, 1.d0)
      endif

      return
      end subroutine bndcy

c    -------------------------------------------
      subroutine bdfrdx(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd = 1+margin
         do j=1,j0
         do i=1,margin
            da(ibnd-i,j) = (da0(ibnd-i,j) - da0(ibnd,j)) + da(ibnd,j)
         enddo
         enddo
      else
         ibnd = i0-margin-men
         do j=1,j0
         do i=1,margin
            da(ibnd+i,j) = (da0(ibnd+i,j) - da0(ibnd,j)) + da(ibnd,j)
         enddo
         enddo
      endif

      return
      end subroutine bdfrdx

c    ---------------------------------------
      subroutine bdfrex(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,j0
         do i=1,margin
            da(i,j) = da(margin+1,j)
         enddo
         enddo
      else
         ibnd1 = i0+1-men
         do j=1,j0
         do i=1,margin
            da(ibnd1-i,j) = da(ibnd1-margin-1,j)
         enddo
         enddo
      endif

      return
      end subroutine bdfrex

c    ---------------------------------------
      subroutine bdperx(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd = 1+margin
         do i=1,margin
            do j=1,j0
               da(ibnd-i,j) = da(i0-margin-i,j)
            enddo
         enddo
      else
         ibnd = i0-margin-men
         do i=1,margin
            do j=1,j0
               if(men .eq. 0) then
                  da(ibnd+i,j) = da(margin+1+i,j)
               else
                  da(ibnd+i,j) = da(margin+i,j)
               endif
            enddo
         enddo
      endif

      return
      end subroutine bdperx

c    --------------------------------------------
      subroutine bdsymx(da,margin,mbnd,men,coff)
c    --------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: coff
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd1 = 1+margin !(4 for men=0,1)
         ibnd2 = 1+margin-men !(4 for men=0, 3 for men=1)
         do i=1,margin
            do j=1,j0
               da(ibnd1-i,j) = coff*da(ibnd2+i,j)
            enddo
         enddo
      else
         ibnd1 = i0-margin !(i0-3 for men=0,1)
         ibnd2 = i0-margin+men !(i0-3 for men=0, i0-2 for men=1)
         do i=1,margin
            do j=1,j0
               da(ibnd1+i,j) = coff*da(ibnd2-i,j)
            enddo
         enddo
      endif

      end subroutine bdsymx

c    -------------------------------------------
      subroutine bdfrdy(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd = 1+margin
         do j=1,margin
         do i=1,i0
            da(i,jbnd-j) = (da0(i,jbnd-j) - da0(i,jbnd)) + da(i,jbnd)
         enddo
         enddo
      else
         jbnd = j0-margin-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd+j) = (da0(i,jbnd+j) - da0(i,jbnd)) + da(i,jbnd)
         enddo
         enddo
      endif

      return
      end subroutine bdfrdy

c    ---------------------------------------
      subroutine bdfrey(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,margin
         do i=1,i0
            da(i,j) = da(i,margin+1)
         enddo
         enddo
      else
         jbnd1 = j0+1-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd1-j) = da(i,jbnd1-margin-1)
         enddo
         enddo
      endif

      return
      end subroutine bdfrey

c    ---------------------------------------
      subroutine bdpery(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd = 1+margin
         do j=1,margin
            do i=1,i0
               da(i,jbnd-j) = da(i,j0-margin-j)
            enddo
         enddo
      else
         jbnd = j0-margin-men
         do j=1,margin
            do i=1,i0
               if(men .eq. 0) then
                  da(i,jbnd+j) = da(i,margin+1+j)
               else
                  da(i,jbnd+j) = da(i,margin+j)
               endif
            enddo
         enddo
      endif

      return
      end subroutine bdpery

c    --------------------------------------------
      subroutine bdsymy(da,margin,mbnd,men,coff)
c    --------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: coff
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd1 = 1+margin
         jbnd2 = 1+margin-men
         do j=1,margin
            do i=1,i0
               da(i,jbnd1-j) = coff*da(i,jbnd2+j)
            enddo
         enddo
      else
         jbnd1 = j0-margin
         jbnd2 = j0-margin+men
         do j=1,margin
            do i=1,i0
               da(i,jbnd1+j) = coff*da(i,jbnd2-j)
            enddo
         enddo
      endif

      end subroutine bdsymy
