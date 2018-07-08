c    -------------------
      subroutine roe_xy
c    -------------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine numfx_roe(is,ie,js,je,nstep)
            integer :: is,ie,js,je,nstep
         end subroutine numfx_roe

         subroutine numfy_roe(is,ie,js,je,nstep)
            integer :: is,ie,js,je,nstep
         end subroutine numfy_roe

         subroutine chck(da,damin,is,ie,js,je)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: damin
            integer :: is,ie,js,je
         end subroutine chck

         subroutine tconduct(is,ie,js,je)
            integer :: is,ie,js,je
         end subroutine tconduct
      end interface

      do j=1,j0
      do i=1,i0-1
         qLx(i,j,1) = ro(i,j)
         qLx(i,j,2) = vx(i,j)
         qLx(i,j,3) = vy(i,j)
         qLx(i,j,4) = vz(i,j)
         qLx(i,j,5) = bx(i,j)
         qLx(i,j,6) = by(i,j)
         qLx(i,j,7) = bz(i,j)
         qLx(i,j,8) = pr(i,j)
c E_in & Enthalpy
         qLx(i,j,9) = qLx(i,j,8)/(gm - 1.d0)/qLx(i,j,1)
         qLx(i,j,10) = 0.5d0*(qLx(i,j,2)**2 + qLx(i,j,3)**2
     &                      + qLx(i,j,4)**2)
     &               + gm*qLx(i,j,8)/(gm - 1.d0)/qLx(i,j,1)
     &               + pi4i*(qLx(i,j,5)**2 + qLx(i,j,6)**2
     &                     + qLx(i,j,7)**2)/qLx(i,j,1)

         qRx(i,j,1) = ro(i+1,j)
         qRx(i,j,2) = vx(i+1,j)
         qRx(i,j,3) = vy(i+1,j)
         qRx(i,j,4) = vz(i+1,j)
         qRx(i,j,5) = bx(i+1,j)
         qRx(i,j,6) = by(i+1,j)
         qRx(i,j,7) = bz(i+1,j)
         qRx(i,j,8) = pr(i+1,j)
c E_in & Enthalpy
         qRx(i,j,9) = qRx(i,j,8)/(gm - 1.d0)/qRx(i,j,1)
         qRx(i,j,10) = 0.5d0*(qRx(i,j,2)**2 + qRx(i,j,3)**2
     &                      + qRx(i,j,4)**2)
     &               + gm*qRx(i,j,8)/(gm - 1.d0)/qRx(i,j,1)
     &               + pi4i*(qRx(i,j,5)**2 + qRx(i,j,6)**2
     &                     + qRx(i,j,7)**2)/qRx(i,j,1)
      enddo
      enddo

      call numfx_roe(1,i0,1,j0,1)

      do j=1,j0
      do i=2,i0-1
         wh(i,j,1) = w(i,j,1)
     &             - 0.5d0*dt*dxmi(i-1)*(-fx(i-1,j,1) + fx(i,j,1))
         wh(i,j,2) = w(i,j,2)
     &             - 0.5d0*dt*(dxmi(i-1)*(-fx(i-1,j,2) + fx(i,j,2))
     &             - 0.5d0*w(i,j,1)*Fexx(i,j))
         wh(i,j,3) = w(i,j,3)
     &             - 0.5d0*dt*(dxmi(i-1)*(-fx(i-1,j,3) + fx(i,j,3))
     &             - 0.5d0*w(i,j,1)*Fexy(i,j))
         wh(i,j,4) = w(i,j,4)
     &             - 0.5d0*dt*(dxmi(i-1)*(-fx(i-1,j,4) + fx(i,j,4))
     &             - 0.5d0*w(i,j,1)*Fexz(i,j))
         wh(i,j,5) = w(i,j,5)
         wh(i,j,6) = w(i,j,6)
     &             - 0.5d0*dt*dxmi(i-1)*(-fx(i-1,j,5) + fx(i,j,5))
         wh(i,j,7) = w(i,j,7)
     &             - 0.5d0*dt*dxmi(i-1)*(-fx(i-1,j,6) + fx(i,j,6))
         wh(i,j,8) = w(i,j,8)
     &             - 0.5d0*dt*dxmi(i-1)*(-fx(i-1,j,7) + fx(i,j,7))
      enddo
      enddo

      do j=1,j0
      do i=2,i0-1
         roh(i,j) = wh(i,j,1)
         vxh(i,j) = wh(i,j,2)/wh(i,j,1)
         vyh(i,j) = wh(i,j,3)/wh(i,j,1)
         vzh(i,j) = wh(i,j,4)/wh(i,j,1)
         bxh(i,j) = wh(i,j,5)
         byh(i,j) = wh(i,j,6)
         bzh(i,j) = wh(i,j,7)
         prh(i,j) = (gm - 1.d0)*(wh(i,j,8) - 0.5d0*roh(i,j)*
     &                        (vxh(i,j)**2 + vyh(i,j)**2 + vzh(i,j)**2)
     &           - 0.5d0*pi4i*(bxh(i,j)**2 + byh(i,j)**2 + bzh(i,j)**2))
      enddo
      enddo

      do j=1,j0
      do i=2,i0-2
         qLx(i,j,1) = roh(i,j)
         qLx(i,j,2) = vxh(i,j)
         qLx(i,j,3) = vyh(i,j)
         qLx(i,j,4) = vzh(i,j)
         qLx(i,j,5) = bxh(i,j)
         qLx(i,j,6) = byh(i,j)
         qLx(i,j,7) = bzh(i,j)
         qLx(i,j,8) = prh(i,j)
c E_in & Enthalpy
         qLx(i,j,9) = qLx(i,j,8)/(gm - 1.d0)/qLx(i,j,1)
         qLx(i,j,10) = 0.5d0*(qLx(i,j,2)**2 + qLx(i,j,3)**2
     &                      + qLx(i,j,4)**2)
     &               + gm*qLx(i,j,8)/(gm - 1.d0)/qLx(i,j,1)
     &               + pi4i*(qLx(i,j,5)**2 + qLx(i,j,6)**2
     &                     + qLx(i,j,7)**2)/qLx(i,j,1)

         qRx(i,j,1) = roh(i+1,j)
         qRx(i,j,2) = vxh(i+1,j)
         qRx(i,j,3) = vyh(i+1,j)
         qRx(i,j,4) = vzh(i+1,j)
         qRx(i,j,5) = bxh(i+1,j)
         qRx(i,j,6) = byh(i+1,j)
         qRx(i,j,7) = bzh(i+1,j)
         qRx(i,j,8) = prh(i+1,j)
c E_in & Enthalpy
         qRx(i,j,9) = qRx(i,j,8)/(gm - 1.d0)/qRx(i,j,1)
         qRx(i,j,10) = 0.5d0*(qRx(i,j,2)**2 + qRx(i,j,3)**2
     &                      + qRx(i,j,4)**2)
     &               + gm*qRx(i,j,8)/(gm - 1.d0)/qRx(i,j,1)
     &               + pi4i*(qRx(i,j,5)**2 + qRx(i,j,6)**2
     &                     + qRx(i,j,7)**2)/qRx(i,j,1)
      enddo
      enddo

      call numfx_roe(3,i0-2,1,j0,2)

      do j=1,j0
      do i=4,i0-3
         w(i,j,1) = w(i,j,1) - dt*dxmi(i-1)*(-fx(i-1,j,1) + fx(i,j,1))
         w(i,j,2) = w(i,j,2) - dt*(dxmi(i-1)*(-fx(i-1,j,2) + fx(i,j,2))
     &            - 0.5d0*w(i,j,1)*Fexx(i,j))
         w(i,j,3) = w(i,j,3) - dt*(dxmi(i-1)*(-fx(i-1,j,3) + fx(i,j,3))
     &            - 0.5d0*w(i,j,1)*Fexy(i,j))
         w(i,j,4) = w(i,j,4) - dt*(dxmi(i-1)*(-fx(i-1,j,4) + fx(i,j,4))
     &            - 0.5d0*w(i,j,1)*Fexz(i,j))
         w(i,j,5) = w(i,j,5)
         w(i,j,6) = w(i,j,6) - dt*dxmi(i-1)*(-fx(i-1,j,5) + fx(i,j,5))
         w(i,j,7) = w(i,j,7) - dt*dxmi(i-1)*(-fx(i-1,j,6) + fx(i,j,6))
         w(i,j,8) = w(i,j,8) - dt*dxmi(i-1)*(-fx(i-1,j,7) + fx(i,j,7))
      enddo
      enddo

      do j=1,j0
      do i=4,i0-3
         ro(i,j) = w(i,j,1)
         vx(i,j) = w(i,j,2)/w(i,j,1)
         vy(i,j) = w(i,j,3)/w(i,j,1)
         vz(i,j) = w(i,j,4)/w(i,j,1)
         bx(i,j) = w(i,j,5)
         by(i,j) = w(i,j,6)
         bz(i,j) = w(i,j,7)
         pr(i,j) = (gm - 1.d0)*(w(i,j,8) - 0.5d0*ro(i,j)
     &                        *(vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2)
     &            - 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2))
      enddo
      enddo

      call bndcx

      call chck(ro,dmin,1,i0,1,j0)
      call chck(pr,pmin,1,i0,1,j0)

      do j=1,j0
      do i=1,i0
         te(i,j) = gm*pr(i,j)/ro(i,j)
      enddo
      enddo

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
     &          + pr(i,j)/(gm - 1.d0)
     &          + 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
      enddo
      enddo

      do j=1,j0-1
      do i=1,i0
         qLy(i,j,1) = ro(i,j)
         qLy(i,j,2) = vy(i,j)
         qLy(i,j,3) = vz(i,j)
         qLy(i,j,4) = vx(i,j)
         qLy(i,j,5) = by(i,j)
         qLy(i,j,6) = bz(i,j)
         qLy(i,j,7) = bx(i,j)
         qLy(i,j,8) = pr(i,j)
c E_in & Enthalpy
         qLy(i,j,9) = qLy(i,j,8)/(gm - 1.d0)/qLy(i,j,1)
         qLy(i,j,10) = 0.5d0*(qLy(i,j,2)**2 + qLy(i,j,3)**2
     &                      + qLy(i,j,4)**2)
     &               + gm*qLy(i,j,8)/(gm - 1.d0)/qLy(i,j,1)
     &               + pi4i*(qLy(i,j,5)**2 + qLy(i,j,6)**2
     &                     + qLy(i,j,7)**2)/qLy(i,j,1)

         qRy(i,j,1) = ro(i,j+1)
         qRy(i,j,2) = vy(i,j+1)
         qRy(i,j,3) = vz(i,j+1)
         qRy(i,j,4) = vx(i,j+1)
         qRy(i,j,5) = by(i,j+1)
         qRy(i,j,6) = bz(i,j+1)
         qRy(i,j,7) = bx(i,j+1)
         qRy(i,j,8) = pr(i,j+1)
c E_in & Enthalpy
         qRy(i,j,9) = qRy(i,j,8)/(gm - 1.d0)/qRy(i,j,1)
         qRy(i,j,10) = 0.5d0*(qRy(i,j,2)**2 + qRy(i,j,3)**2
     &                      + qRy(i,j,4)**2)
     &               + gm*qRy(i,j,8)/(gm - 1.d0)/qRy(i,j,1)
     &               + pi4i*(qRy(i,j,5)**2 + qRy(i,j,6)**2
     &                     + qRy(i,j,7)**2)/qRy(i,j,1)
      enddo
      enddo

      call numfy_roe(1,i0,1,j0,1)

      do j=2,j0-1
      do i=1,i0
         wh(i,j,1) = w(i,j,1)
     &             - 0.5d0*dt*dymi(j-1)*(-fy(i,j-1,1) + fy(i,j,1))
         wh(i,j,2) = w(i,j,2)
     &             - 0.5d0*dt*(dymi(j-1)*(-fy(i,j-1,4) + fy(i,j,4))
     &             - 0.5d0*w(i,j,1)*Fexx(i,j))
         wh(i,j,3) = w(i,j,3)
     &             - 0.5d0*dt*(dymi(j-1)*(-fy(i,j-1,2) + fy(i,j,2))
     &             - 0.5d0*w(i,j,1)*Fexy(i,j))
         wh(i,j,4) = w(i,j,4)
     &             - 0.5d0*dt*(dymi(j-1)*(-fy(i,j-1,3) + fy(i,j,3))
     &             - 0.5d0*w(i,j,1)*Fexz(i,j))
         wh(i,j,5) = w(i,j,5)
     &             - 0.5d0*dt*dymi(j-1)*(-fy(i,j-1,6) + fy(i,j,6))
         wh(i,j,6) = w(i,j,6)
         wh(i,j,7) = w(i,j,7)
     &             - 0.5d0*dt*dymi(j-1)*(-fy(i,j-1,5) + fy(i,j,5))
         wh(i,j,8) = w(i,j,8)
     &             - 0.5d0*dt*dymi(j-1)*(-fy(i,j-1,7) + fy(i,j,7))
      enddo
      enddo

      do j=2,j0-1
      do i=1,i0
         roh(i,j) = wh(i,j,1)
         vxh(i,j) = wh(i,j,2)/wh(i,j,1)
         vyh(i,j) = wh(i,j,3)/wh(i,j,1)
         vzh(i,j) = wh(i,j,4)/wh(i,j,1)
         bxh(i,j) = wh(i,j,5)
         byh(i,j) = wh(i,j,6)
         bzh(i,j) = wh(i,j,7)
         prh(i,j) = (gm - 1.d0)*(wh(i,j,8) - 0.5d0*roh(i,j)*
     &                        (vxh(i,j)**2 + vyh(i,j)**2 + vzh(i,j)**2)
     &           - 0.5d0*pi4i*(bxh(i,j)**2 + byh(i,j)**2 + bzh(i,j)**2))
      enddo
      enddo

      do j=2,j0-2
      do i=1,i0
         qLy(i,j,1) = roh(i,j)
         qLy(i,j,2) = vyh(i,j)
         qLy(i,j,3) = vzh(i,j)
         qLy(i,j,4) = vxh(i,j)
         qLy(i,j,5) = byh(i,j)
         qLy(i,j,6) = bzh(i,j)
         qLy(i,j,7) = bxh(i,j)
         qLy(i,j,8) = prh(i,j)
c E_in & Enthalpy
         qLy(i,j,9) = qLy(i,j,8)/(gm - 1.d0)/qLy(i,j,1)
         qLy(i,j,10) = 0.5d0*(qLy(i,j,2)**2 + qLy(i,j,3)**2
     &                      + qLy(i,j,4)**2)
     &               + gm*qLy(i,j,8)/(gm - 1.d0)/qLy(i,j,1)
     &               + pi4i*(qLy(i,j,5)**2 + qLy(i,j,6)**2
     &                     + qLy(i,j,7)**2)/qLy(i,j,1)

         qRy(i,j,1) = roh(i,j+1)
         qRy(i,j,2) = vyh(i,j+1)
         qRy(i,j,3) = vzh(i,j+1)
         qRy(i,j,4) = vxh(i,j+1)
         qRy(i,j,5) = byh(i,j+1)
         qRy(i,j,6) = bzh(i,j+1)
         qRy(i,j,7) = bxh(i,j+1)
         qRy(i,j,8) = prh(i,j+1)
c E_in & Enthalpy
         qRy(i,j,9) = qRy(i,j,8)/(gm - 1.d0)/qRy(i,j,1)
         qRy(i,j,10) = 0.5d0*(qRy(i,j,2)**2 + qRy(i,j,3)**2
     &                      + qRy(i,j,4)**2)
     &               + gm*qRy(i,j,8)/(gm - 1.d0)/qRy(i,j,1)
     &               + pi4i*(qRy(i,j,5)**2 + qRy(i,j,6)**2
     &                     + qRy(i,j,7)**2)/qRy(i,j,1)
      enddo
      enddo

      call numfy_roe(1,i0,3,j0-2,2)

      do j=4,j0-3
      do i=1,i0
         w(i,j,1) = w(i,j,1) - dt*dymi(j-1)*(-fy(i,j-1,1) + fy(i,j,1))
         w(i,j,2) = w(i,j,2) - dt*(dymi(j-1)*(-fy(i,j-1,4) + fy(i,j,4))
     &            - 0.5d0*w(i,j,1)*Fexx(i,j))
         w(i,j,3) = w(i,j,3) - dt*(dymi(j-1)*(-fy(i,j-1,2) + fy(i,j,2))
     &            - 0.5d0*w(i,j,1)*Fexy(i,j))
         w(i,j,4) = w(i,j,4) - dt*(dymi(j-1)*(-fy(i,j-1,3) + fy(i,j,3))
     &            - 0.5d0*w(i,j,1)*Fexz(i,j))
         w(i,j,5) = w(i,j,5) - dt*dymi(j-1)*(-fy(i,j-1,6) + fy(i,j,6))
         w(i,j,6) = w(i,j,6)
         w(i,j,7) = w(i,j,7) - dt*dymi(j-1)*(-fy(i,j-1,5) + fy(i,j,5))
         w(i,j,8) = w(i,j,8) - dt*dymi(j-1)*(-fy(i,j-1,7) + fy(i,j,7))
      enddo
      enddo

      do j=4,j0-3
      do i=1,i0
         ro(i,j) = w(i,j,1)
         vx(i,j) = w(i,j,2)/w(i,j,1)
         vy(i,j) = w(i,j,3)/w(i,j,1)
         vz(i,j) = w(i,j,4)/w(i,j,1)
         bx(i,j) = w(i,j,5)
         by(i,j) = w(i,j,6)
         bz(i,j) = w(i,j,7)
         pr(i,j) = (gm - 1.d0)*(w(i,j,8) - 0.5d0*ro(i,j)
     &                        *(vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2)
     &            - 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2))
      enddo
      enddo

      call bndcy

      call chck(ro,dmin,1,i0,1,j0)
      call chck(pr,pmin,1,i0,1,j0)

C Calculate the thermal conduction
      do j=1,j0
      do i=1,i0
         te(i,j) = gm*pr(i,j)/ro(i,j)
      enddo
      enddo
      if(ntcon .eq. 1) call tconduct(2,i0-1,2,j0-1)

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
     &          + pr(i,j)/(gm - 1.d0)
     &          + 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
      enddo
      enddo

      return
      end subroutine roe_xy

c    -------------------
      subroutine roe_yx
c    -------------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine numfx_roe(is,ie,js,je,nstep)
            integer :: is,ie,js,je,nstep
         end subroutine numfx_roe

         subroutine numfy_roe(is,ie,js,je,nstep)
            integer :: is,ie,js,je,nstep
         end subroutine numfy_roe

         subroutine chck(da,damin,is,ie,js,je)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: damin
            integer :: is,ie,js,je
         end subroutine chck

         subroutine tconduct(is,ie,js,je)
            integer :: is,ie,js,je
         end subroutine tconduct
      end interface

      do j=1,j0-1
      do i=1,i0
         qLy(i,j,1) = ro(i,j)
         qLy(i,j,2) = vy(i,j)
         qLy(i,j,3) = vz(i,j)
         qLy(i,j,4) = vx(i,j)
         qLy(i,j,5) = by(i,j)
         qLy(i,j,6) = bz(i,j)
         qLy(i,j,7) = bx(i,j)
         qLy(i,j,8) = pr(i,j)
c E_in & Enthalpy
         qLy(i,j,9) = qLy(i,j,8)/(gm - 1.d0)/qLy(i,j,1)
         qLy(i,j,10) = 0.5d0*(qLy(i,j,2)**2 + qLy(i,j,3)**2
     &                      + qLy(i,j,4)**2)
     &               + gm*qLy(i,j,8)/(gm - 1.d0)/qLy(i,j,1)
     &               + pi4i*(qLy(i,j,5)**2 + qLy(i,j,6)**2
     &                     + qLy(i,j,7)**2)/qLy(i,j,1)

         qRy(i,j,1) = ro(i,j+1)
         qRy(i,j,2) = vy(i,j+1)
         qRy(i,j,3) = vz(i,j+1)
         qRy(i,j,4) = vx(i,j+1)
         qRy(i,j,5) = by(i,j+1)
         qRy(i,j,6) = bz(i,j+1)
         qRy(i,j,7) = bx(i,j+1)
         qRy(i,j,8) = pr(i,j+1)
c E_in & Enthalpy
         qRy(i,j,9) = qRy(i,j,8)/(gm - 1.d0)/qRy(i,j,1)
         qRy(i,j,10) = 0.5d0*(qRy(i,j,2)**2 + qRy(i,j,3)**2
     &                      + qRy(i,j,4)**2)
     &               + gm*qRy(i,j,8)/(gm - 1.d0)/qRy(i,j,1)
     &               + pi4i*(qRy(i,j,5)**2 + qRy(i,j,6)**2
     &                     + qRy(i,j,7)**2)/qRy(i,j,1)
      enddo
      enddo

      call numfy_roe(1,i0,1,j0,1)

      do j=2,j0-1
      do i=1,i0
         wh(i,j,1) = w(i,j,1)
     &             - 0.5d0*dt*dymi(j-1)*(-fy(i,j-1,1) + fy(i,j,1))
         wh(i,j,2) = w(i,j,2)
     &             - 0.5d0*dt*(dymi(j-1)*(-fy(i,j-1,4) + fy(i,j,4))
     &             - 0.5d0*w(i,j,1)*Fexx(i,j))
         wh(i,j,3) = w(i,j,3)
     &             - 0.5d0*dt*(dymi(j-1)*(-fy(i,j-1,2) + fy(i,j,2))
     &             - 0.5d0*w(i,j,1)*Fexy(i,j))
         wh(i,j,4) = w(i,j,4)
     &             - 0.5d0*dt*(dymi(j-1)*(-fy(i,j-1,3) + fy(i,j,3))
     &             - 0.5d0*w(i,j,1)*Fexz(i,j))
         wh(i,j,5) = w(i,j,5)
     &             - 0.5d0*dt*dymi(j-1)*(-fy(i,j-1,6) + fy(i,j,6))
         wh(i,j,6) = w(i,j,6)
         wh(i,j,7) = w(i,j,7)
     &             - 0.5d0*dt*dymi(j-1)*(-fy(i,j-1,5) + fy(i,j,5))
         wh(i,j,8) = w(i,j,8)
     &             - 0.5d0*dt*dymi(j-1)*(-fy(i,j-1,7) + fy(i,j,7))
      enddo
      enddo

      do j=2,j0-1
      do i=1,i0
         roh(i,j) = wh(i,j,1)
         vxh(i,j) = wh(i,j,2)/wh(i,j,1)
         vyh(i,j) = wh(i,j,3)/wh(i,j,1)
         vzh(i,j) = wh(i,j,4)/wh(i,j,1)
         bxh(i,j) = wh(i,j,5)
         byh(i,j) = wh(i,j,6)
         bzh(i,j) = wh(i,j,7)
         prh(i,j) = (gm - 1.d0)*(wh(i,j,8) - 0.5d0*roh(i,j)*
     &                        (vxh(i,j)**2 + vyh(i,j)**2 + vzh(i,j)**2)
     &           - 0.5d0*pi4i*(bxh(i,j)**2 + byh(i,j)**2 + bzh(i,j)**2))
      enddo
      enddo

      do j=2,j0-2
      do i=1,i0
         qLy(i,j,1) = roh(i,j)
         qLy(i,j,2) = vyh(i,j)
         qLy(i,j,3) = vzh(i,j)
         qLy(i,j,4) = vxh(i,j)
         qLy(i,j,5) = byh(i,j)
         qLy(i,j,6) = bzh(i,j)
         qLy(i,j,7) = bxh(i,j)
         qLy(i,j,8) = prh(i,j)
c E_in & Enthalpy
         qLy(i,j,9) = qLy(i,j,8)/(gm - 1.d0)/qLy(i,j,1)
         qLy(i,j,10) = 0.5d0*(qLy(i,j,2)**2 + qLy(i,j,3)**2
     &                      + qLy(i,j,4)**2)
     &               + gm*qLy(i,j,8)/(gm - 1.d0)/qLy(i,j,1)
     &               + pi4i*(qLy(i,j,5)**2 + qLy(i,j,6)**2
     &                     + qLy(i,j,7)**2)/qLy(i,j,1)

         qRy(i,j,1) = roh(i,j+1)
         qRy(i,j,2) = vyh(i,j+1)
         qRy(i,j,3) = vzh(i,j+1)
         qRy(i,j,4) = vxh(i,j+1)
         qRy(i,j,5) = byh(i,j+1)
         qRy(i,j,6) = bzh(i,j+1)
         qRy(i,j,7) = bxh(i,j+1)
         qRy(i,j,8) = prh(i,j+1)
c E_in & Enthalpy
         qRy(i,j,9) = qRy(i,j,8)/(gm - 1.d0)/qRy(i,j,1)
         qRy(i,j,10) = 0.5d0*(qRy(i,j,2)**2 + qRy(i,j,3)**2
     &                      + qRy(i,j,4)**2)
     &               + gm*qRy(i,j,8)/(gm - 1.d0)/qRy(i,j,1)
     &               + pi4i*(qRy(i,j,5)**2 + qRy(i,j,6)**2
     &                     + qRy(i,j,7)**2)/qRy(i,j,1)
      enddo
      enddo

      call numfy_roe(1,i0,3,j0-2,2)

      do j=4,j0-3
      do i=1,i0
         w(i,j,1) = w(i,j,1) - dt*dymi(j-1)*(-fy(i,j-1,1) + fy(i,j,1))
         w(i,j,2) = w(i,j,2) - dt*(dymi(j-1)*(-fy(i,j-1,4) + fy(i,j,4))
     &            - 0.5d0*w(i,j,1)*Fexx(i,j))
         w(i,j,3) = w(i,j,3) - dt*(dymi(j-1)*(-fy(i,j-1,2) + fy(i,j,2))
     &            - 0.5d0*w(i,j,1)*Fexy(i,j))
         w(i,j,4) = w(i,j,4) - dt*(dymi(j-1)*(-fy(i,j-1,3) + fy(i,j,3))
     &            - 0.5d0*w(i,j,1)*Fexz(i,j))
         w(i,j,5) = w(i,j,5) - dt*dymi(j-1)*(-fy(i,j-1,6) + fy(i,j,6))
         w(i,j,6) = w(i,j,6)
         w(i,j,7) = w(i,j,7) - dt*dymi(j-1)*(-fy(i,j-1,5) + fy(i,j,5))
         w(i,j,8) = w(i,j,8) - dt*dymi(j-1)*(-fy(i,j-1,7) + fy(i,j,7))
      enddo
      enddo

      do j=4,j0-3
      do i=1,i0
         ro(i,j) = w(i,j,1)
         vx(i,j) = w(i,j,2)/w(i,j,1)
         vy(i,j) = w(i,j,3)/w(i,j,1)
         vz(i,j) = w(i,j,4)/w(i,j,1)
         bx(i,j) = w(i,j,5)
         by(i,j) = w(i,j,6)
         bz(i,j) = w(i,j,7)
         pr(i,j) = (gm - 1.d0)*(w(i,j,8) - 0.5d0*ro(i,j)
     &                        *(vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2)
     &            - 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2))
      enddo
      enddo

      call bndcy

      call chck(ro,dmin,1,i0,1,j0)
      call chck(pr,pmin,1,i0,1,j0)

      do j=1,j0
      do i=1,i0
         te(i,j) = gm*pr(i,j)/ro(i,j)
      enddo
      enddo

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
     &          + pr(i,j)/(gm - 1.d0)
     &          + 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
      enddo
      enddo

      do j=1,j0
      do i=1,i0-1
         qLx(i,j,1) = ro(i,j)
         qLx(i,j,2) = vx(i,j)
         qLx(i,j,3) = vy(i,j)
         qLx(i,j,4) = vz(i,j)
         qLx(i,j,5) = bx(i,j)
         qLx(i,j,6) = by(i,j)
         qLx(i,j,7) = bz(i,j)
         qLx(i,j,8) = pr(i,j)
c E_in & Enthalpy
         qLx(i,j,9) = qLx(i,j,8)/(gm - 1.d0)/qLx(i,j,1)
         qLx(i,j,10) = 0.5d0*(qLx(i,j,2)**2 + qLx(i,j,3)**2
     &                      + qLx(i,j,4)**2)
     &               + gm*qLx(i,j,8)/(gm - 1.d0)/qLx(i,j,1)
     &               + pi4i*(qLx(i,j,5)**2 + qLx(i,j,6)**2
     &                     + qLx(i,j,7)**2)/qLx(i,j,1)

         qRx(i,j,1) = ro(i+1,j)
         qRx(i,j,2) = vx(i+1,j)
         qRx(i,j,3) = vy(i+1,j)
         qRx(i,j,4) = vz(i+1,j)
         qRx(i,j,5) = bx(i+1,j)
         qRx(i,j,6) = by(i+1,j)
         qRx(i,j,7) = bz(i+1,j)
         qRx(i,j,8) = pr(i+1,j)
c E_in & Enthalpy
         qRx(i,j,9) = qRx(i,j,8)/(gm - 1.d0)/qRx(i,j,1)
         qRx(i,j,10) = 0.5d0*(qRx(i,j,2)**2 + qRx(i,j,3)**2
     &                      + qRx(i,j,4)**2)
     &               + gm*qRx(i,j,8)/(gm - 1.d0)/qRx(i,j,1)
     &               + pi4i*(qRx(i,j,5)**2 + qRx(i,j,6)**2
     &                     + qRx(i,j,7)**2)/qRx(i,j,1)
      enddo
      enddo

      call numfx_roe(1,i0,1,j0,1)

      do j=1,j0
      do i=2,i0-1
         wh(i,j,1) = w(i,j,1)
     &             - 0.5d0*dt*dxmi(i-1)*(-fx(i-1,j,1) + fx(i,j,1))
         wh(i,j,2) = w(i,j,2)
     &             - 0.5d0*dt*(dxmi(i-1)*(-fx(i-1,j,2) + fx(i,j,2))
     &             - 0.5d0*w(i,j,1)*Fexx(i,j))
         wh(i,j,3) = w(i,j,3)
     &             - 0.5d0*dt*(dxmi(i-1)*(-fx(i-1,j,3) + fx(i,j,3))
     &             - 0.5d0*w(i,j,1)*Fexy(i,j))
         wh(i,j,4) = w(i,j,4)
     &             - 0.5d0*dt*(dxmi(i-1)*(-fx(i-1,j,4) + fx(i,j,4))
     &             - 0.5d0*w(i,j,1)*Fexz(i,j))
         wh(i,j,5) = w(i,j,5)
         wh(i,j,6) = w(i,j,6)
     &             - 0.5d0*dt*dxmi(i-1)*(-fx(i-1,j,5) + fx(i,j,5))
         wh(i,j,7) = w(i,j,7)
     &             - 0.5d0*dt*dxmi(i-1)*(-fx(i-1,j,6) + fx(i,j,6))
         wh(i,j,8) = w(i,j,8)
     &             - 0.5d0*dt*dxmi(i-1)*(-fx(i-1,j,7) + fx(i,j,7))
      enddo
      enddo

      do j=1,j0
      do i=2,i0-1
         roh(i,j) = wh(i,j,1)
         vxh(i,j) = wh(i,j,2)/wh(i,j,1)
         vyh(i,j) = wh(i,j,3)/wh(i,j,1)
         vzh(i,j) = wh(i,j,4)/wh(i,j,1)
         bxh(i,j) = wh(i,j,5)
         byh(i,j) = wh(i,j,6)
         bzh(i,j) = wh(i,j,7)
         prh(i,j) = (gm - 1.d0)*(wh(i,j,8) - 0.5d0*roh(i,j)*
     &                        (vxh(i,j)**2 + vyh(i,j)**2 + vzh(i,j)**2)
     &           - 0.5d0*pi4i*(bxh(i,j)**2 + byh(i,j)**2 + bzh(i,j)**2))
      enddo
      enddo

      do j=1,j0
      do i=2,i0-2
         qLx(i,j,1) = roh(i,j)
         qLx(i,j,2) = vxh(i,j)
         qLx(i,j,3) = vyh(i,j)
         qLx(i,j,4) = vzh(i,j)
         qLx(i,j,5) = bxh(i,j)
         qLx(i,j,6) = byh(i,j)
         qLx(i,j,7) = bzh(i,j)
         qLx(i,j,8) = prh(i,j)
c E_in & Enthalpy
         qLx(i,j,9) = qLx(i,j,8)/(gm - 1.d0)/qLx(i,j,1)
         qLx(i,j,10) = 0.5d0*(qLx(i,j,2)**2 + qLx(i,j,3)**2
     &                      + qLx(i,j,4)**2)
     &               + gm*qLx(i,j,8)/(gm - 1.d0)/qLx(i,j,1)
     &               + pi4i*(qLx(i,j,5)**2 + qLx(i,j,6)**2
     &                     + qLx(i,j,7)**2)/qLx(i,j,1)

         qRx(i,j,1) = roh(i+1,j)
         qRx(i,j,2) = vxh(i+1,j)
         qRx(i,j,3) = vyh(i+1,j)
         qRx(i,j,4) = vzh(i+1,j)
         qRx(i,j,5) = bxh(i+1,j)
         qRx(i,j,6) = byh(i+1,j)
         qRx(i,j,7) = bzh(i+1,j)
         qRx(i,j,8) = prh(i+1,j)
c E_in & Enthalpy
         qRx(i,j,9) = qRx(i,j,8)/(gm - 1.d0)/qRx(i,j,1)
         qRx(i,j,10) = 0.5d0*(qRx(i,j,2)**2 + qRx(i,j,3)**2
     &                      + qRx(i,j,4)**2)
     &               + gm*qRx(i,j,8)/(gm - 1.d0)/qRx(i,j,1)
     &               + pi4i*(qRx(i,j,5)**2 + qRx(i,j,6)**2
     &                     + qRx(i,j,7)**2)/qRx(i,j,1)
      enddo
      enddo

      call numfx_roe(3,i0-2,1,j0,2)

      do j=1,j0
      do i=4,i0-3
         w(i,j,1) = w(i,j,1) - dt*dxmi(i-1)*(-fx(i-1,j,1) + fx(i,j,1))
         w(i,j,2) = w(i,j,2) - dt*(dxmi(i-1)*(-fx(i-1,j,2) + fx(i,j,2))
     &            - 0.5d0*w(i,j,1)*Fexx(i,j))
         w(i,j,3) = w(i,j,3) - dt*(dxmi(i-1)*(-fx(i-1,j,3) + fx(i,j,3))
     &            - 0.5d0*w(i,j,1)*Fexy(i,j))
         w(i,j,4) = w(i,j,4) - dt*(dxmi(i-1)*(-fx(i-1,j,4) + fx(i,j,4))
     &            - 0.5d0*w(i,j,1)*Fexz(i,j))
         w(i,j,5) = w(i,j,5)
         w(i,j,6) = w(i,j,6) - dt*dxmi(i-1)*(-fx(i-1,j,5) + fx(i,j,5))
         w(i,j,7) = w(i,j,7) - dt*dxmi(i-1)*(-fx(i-1,j,6) + fx(i,j,6))
         w(i,j,8) = w(i,j,8) - dt*dxmi(i-1)*(-fx(i-1,j,7) + fx(i,j,7))
      enddo
      enddo

      do j=1,j0
      do i=4,i0-3
         ro(i,j) = w(i,j,1)
         vx(i,j) = w(i,j,2)/w(i,j,1)
         vy(i,j) = w(i,j,3)/w(i,j,1)
         vz(i,j) = w(i,j,4)/w(i,j,1)
         bx(i,j) = w(i,j,5)
         by(i,j) = w(i,j,6)
         bz(i,j) = w(i,j,7)
         pr(i,j) = (gm - 1.d0)*(w(i,j,8) - 0.5d0*ro(i,j)
     &                        *(vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2)
     &            - 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2))
      enddo
      enddo

      call bndcx

      call chck(ro,dmin,1,i0,1,j0)
      call chck(pr,pmin,1,i0,1,j0)

C Calculate the thermal conduction
      do j=1,j0
      do i=1,i0
         te(i,j) = gm*pr(i,j)/ro(i,j)
      enddo
      enddo
      if(ntcon .eq. 1) call tconduct(2,i0-1,2,j0-1)

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
     &          + pr(i,j)/(gm - 1.d0)
     &          + 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
      enddo
      enddo

      return
      end subroutine roe_yx

