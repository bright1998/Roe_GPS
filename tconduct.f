c    ----------------------------------
      subroutine tconduct(is,ie,js,je)
c    ----------------------------------
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

      dimension :: s(1:i0,1:j0),c1(1:i0,1:j0),c2(1:i0,1:j0),
     &                          d1(1:i0,1:j0),d2(1:i0,1:j0),e(1:i0,1:j0)
      dimension :: rmu(1:i0,1:j0)
      dimension :: r(1:i0,1:j0),work(1:i0,1:j0)
C Magnetic Field and Temperature at (i+1/2,j)
      dimension :: bxmx(1:i0,1:j0),bymx(1:i0,1:j0),bzmx(1:i0,1:j0),
     &             b2mx(1:i0,1:j0),temx(1:i0,1:j0)
C Magnetic Field and Temperature at (i,j+1/2)
      dimension :: bxmy(1:i0,1:j0),bymy(1:i0,1:j0),bzmy(1:i0,1:j0),
     &             b2my(1:i0,1:j0),temy(1:i0,1:j0)
C Temperature at (i+1/2,j+1/2)
      dimension :: tem(1:i0,1:j0)
      integer :: is,ie,js,je

      dti = 1.d0/dt

      do j=1,j0
      do i=1,i0
         rmu(i,j) = peci/ro(i,j)
      enddo
      enddo

C Calculate Magnetic Field and Temperature at (i+1/2,j)
      do j=1,j0
      do i=1,i0-1
         bxmx(i,j) = 0.5d0*(bx(i,j) + bx(i+1,j))
         bymx(i,j) = 0.5d0*(by(i,j) + by(i+1,j))
         bzmx(i,j) = 0.5d0*(bz(i,j) + bz(i+1,j))
         b2mx(i,j) = bxmx(i,j)**2 + bymx(i,j)**2 + bzmx(i,j)**2
c         temx(i,j) = 0.5d0*(te(i,j) + te(i+1,j))
         temx(i,j) = dsqrt(te(i,j)*te(i+1,j))
      enddo
      enddo

C Calculate Magnetic Field and Temperature at (i,j+1/2)
      do j=1,j0-1
      do i=1,i0
         bxmy(i,j) = 0.5d0*(bx(i,j) + bx(i,j+1))
         bymy(i,j) = 0.5d0*(by(i,j) + by(i,j+1))
         bzmy(i,j) = 0.5d0*(bz(i,j) + bz(i,j+1))
         b2my(i,j) = bxmy(i,j)**2 + bymy(i,j)**2 + bzmy(i,j)**2
c         temy(i,j) = 0.5d0*(te(i,j) + te(i,j+1))
         temy(i,j) = dsqrt(te(i,j)*te(i,j+1))
      enddo
      enddo

C Calculate Temperature at (i+1/2,j+1/2)
      do j=1,j0-1
      do i=1,i0-1
         tem(i,j) = 0.25d0*(te(i,j) + te(i+1,j)
     &                    + te(i,j+1) + te(i+1,j+1))
      enddo
      enddo

C Usually is=2,ie=i0-1,js=2,je=j0-1
      do j=js,je
      do i=is,ie
         rkappam1 = bxmx(i,j)*bymx(i,j)*temx(i,j)**2.5/b2mx(i,j)
         rkappam2 = bxmx(i-1,j)*bymx(i-1,j)*temx(i-1,j)**2.5/b2mx(i-1,j)
         rkappam3 = bxmy(i,j)*bymy(i,j)*temy(i,j)**2.5/b2my(i,j)
         rkappam4 = bxmy(i,j-1)*bymy(i,j-1)*temy(i,j-1)**2.5/b2my(i,j-1)

         rkappama = dym(j-1)*dxi(i)*bxmx(i,j)**2*temx(i,j)**2.5
     &/b2mx(i,j)
         rkappamb = dym(j-1)*dxi(i-1)*bxmx(i-1,j)**2*temx(i-1,j)**2.5
     &/b2mx(i-1,j)
         rkappamc = dxm(i-1)*dyi(j)*bymy(i,j)**2*temy(i,j)**2.5
     &/b2my(i,j)
         rkappamd = dxm(i-1)*dyi(j-1)*bymy(i,j-1)**2*temy(i,j-1)**2.5
     &/b2my(i,j-1)

         s(i,j) = dxm(i-1)*dym(j-1)*dti*te(i,j)
     &          + rmu(i,j)*((rkappam1 + rkappam3)*tem(i,j)
     &                    - (rkappam1 + rkappam4)*tem(i,j-1)
     &                    - (rkappam2 + rkappam3)*tem(i-1,j)
     &                    + (rkappam2 + rkappam4)*tem(i-1,j-1))

         d1(i,j) =-rkappama*rmu(i,j)
         c1(i,j) =-rkappamb*rmu(i,j)

         d2(i,j) =-rkappamc*rmu(i,j)
         c2(i,j) =-rkappamd*rmu(i,j)

         e(i,j) = dxm(i-1)*dym(j-1)*dti
     &          - (d1(i,j) + c1(i,j) + d2(i,j) + c2(i,j))
      enddo
      enddo

      sigmar0 = 0.d0
      do j=js,je
      do i=is,ie
         r(i,j) = e(i,j)*te(i,j) + d1(i,j)*te(i+1,j) + c1(i,j)*te(i-1,j)
     &                           + d2(i,j)*te(i,j+1) + c2(i,j)*te(i,j-1)
     &          - s(i,j)
         sigmar0 = sigmar0 + abs(r(i,j))
      enddo
      enddo

      do j=1,j0
      do i=1,i0
         work(i,j) = te(i,j)
      enddo
      enddo

      do n=1,itemax

         do j=js,je,2
         do i=is,ie,2
            r(i,j) = e(i,j)*work(i,j)
     &             + d1(i,j)*work(i+1,j) + c1(i,j)*work(i-1,j)
     &             + d2(i,j)*work(i,j+1) + c2(i,j)*work(i,j-1)
     &             - s(i,j)
         enddo
         enddo
         do j=js,je,2
         do i=is,ie,2
            work(i,j) = work(i,j) - r(i,j)/e(i,j)
         enddo
         enddo

         do j=js+1,je,2
         do i=is+1,ie,2
            r(i,j) = e(i,j)*work(i,j)
     &             + d1(i,j)*work(i+1,j) + c1(i,j)*work(i-1,j)
     &             + d2(i,j)*work(i,j+1) + c2(i,j)*work(i,j-1)
     &             - s(i,j)
         enddo
         enddo
         do j=js+1,je,2
         do i=is+1,ie,2
            work(i,j) = work(i,j) - r(i,j)/e(i,j)
         enddo
         enddo

         do j=js,je,2
         do i=is+1,ie,2
            r(i,j) = e(i,j)*work(i,j)
     &             + d1(i,j)*work(i+1,j) + c1(i,j)*work(i-1,j)
     &             + d2(i,j)*work(i,j+1) + c2(i,j)*work(i,j-1)
     &             - s(i,j)
         enddo
         enddo
         do j=js,je,2
         do i=is+1,ie,2
            work(i,j) = work(i,j) - r(i,j)/e(i,j)
         enddo
         enddo

         do j=js+1,je,2
         do i=is,ie,2
            r(i,j) = e(i,j)*work(i,j)
     &             + d1(i,j)*work(i+1,j) + c1(i,j)*work(i-1,j)
     &             + d2(i,j)*work(i,j+1) + c2(i,j)*work(i,j-1)
     &             - s(i,j)
         enddo
         enddo
         do j=js+1,je,2
         do i=is,ie,2
            work(i,j) = work(i,j) - r(i,j)/e(i,j)
         enddo
         enddo

C Boundary condition for continuing the iteration
         do i=1,i0
            work(i,1) = work(i,2)
            work(i,j0) = work(i,j0-1)
         enddo
         do j=1,j0
            work(1,j) = work(2,j)
            work(i0,j) = work(i0-1,j)
         enddo

         sigmar = 0.d0
         do j=js,je
         do i=is,ie
            sigmar = sigmar + abs(r(i,j))
         enddo
         enddo

         if(sigmar/sigmar0 .lt. epsit) goto 100
      enddo

 100  continue
      do j=js,je
      do i=is,ie
         te(i,j) = work(i,j)
      enddo
      enddo

C Boundary condition
      if(nbnd1 .eq. 1) then
         call bdfrdx(te,te0,3,0,0)
      elseif(nbnd1 .eq. 2) then
         call bdfrex(te,3,0,0)
      elseif(nbnd1 .eq. 3) then
         call bdperx(te,3,0,0)
      elseif(nbnd1 .eq. 4) then
         call bdsymx(te,3,0,0, 1.d0)
      endif

      if(nbnd2 .eq. 1) then
         call bdfrdx(te,te0,3,1,0)
      elseif(nbnd2 .eq. 2) then
         call bdfrex(te,3,1,0)
      elseif(nbnd2 .eq. 3) then
         call bdperx(te,3,1,0)
      elseif(nbnd2 .eq. 4) then
         call bdsymx(te,3,1,0, 1.d0)
      endif

      if(nbnd3 .eq. 1) then
         call bdfrdy(te,te0,3,0,0)
      elseif(nbnd3 .eq. 2) then
         call bdfrey(te,3,0,0)
      elseif(nbnd3 .eq. 3) then
         call bdpery(te,3,0,0)
      elseif(nbnd3 .eq. 4) then
         call bdsymy(te,3,0,0, 1.d0)
      endif

      if(nbnd4 .eq. 1) then
         call bdfrdy(te,te0,3,1,0)
      elseif(nbnd4 .eq. 2) then
         call bdfrey(te,3,1,0)
      elseif(nbnd4 .eq. 3) then
         call bdpery(te,3,1,0)
      elseif(nbnd4 .eq. 4) then
         call bdsymy(te,3,1,0, 1.d0)
      endif

      do j=1,j0
      do i=1,i0
         pr(i,j) = te(i,j)*ro(i,j)/gm
      enddo
      enddo

      return
      end subroutine tconduct
