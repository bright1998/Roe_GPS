c    -----------------
      subroutine init
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine chck(da,damin,is,ie,js,je)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: damin
            integer :: is,ie,js,je
         end subroutine chck
      end interface

      x0 = 0.5d0
      y0 = 0.5d0
      do j=1,j0
      do i=1,i0
C Sod Shock Tube Problem
cc         if(x(i) .le. x0) then
c         if(y(j) .le. y0) then
c            ro(i,j) = 1.d0
c            vx(i,j) = 0.d0
c            vy(i,j) = 0.d0
c            vz(i,j) = 0.d0
c            bx(i,j) = 0.d0
c            by(i,j) = 0.d0
c            bz(i,j) = 0.d0
c            pr(i,j) = 1.d0
c         else
c            ro(i,j) = 0.125d0
c            vx(i,j) = 0.d0
c            vy(i,j) = 0.d0
c            vz(i,j) = 0.d0
c            bx(i,j) = 0.d0
c            by(i,j) = 0.d0
c            bz(i,j) = 0.d0
c            pr(i,j) = 0.1d0
c         endif
C Brio & Wu
c         if(x(i) .le. x0) then
         if(y(j) .le. y0) then
            ro(i,j) = 1.d0
            vx(i,j) = 0.d0
            vy(i,j) = 0.d0
            vz(i,j) = 0.d0
C X-direction Test
c            bx(i,j) = 0.75d0*dsqrt(pi4)
c            by(i,j) = dsqrt(pi4)
c            bz(i,j) = 0.d0
C Y-direction Test
            bx(i,j) = 0.d0
            by(i,j) = 0.75d0*dsqrt(pi4)
            bz(i,j) = dsqrt(pi4)
            pr(i,j) = 1.d0
         else
            ro(i,j) = 0.125d0
            vx(i,j) = 0.d0
            vy(i,j) = 0.d0
            vz(i,j) = 0.d0
C X-direction Test
c            bx(i,j) = 0.75d0*dsqrt(pi4)
c            by(i,j) =-dsqrt(pi4)
c            bz(i,j) = 0.d0
C Y-direction Test
            bx(i,j) = 0.d0
            by(i,j) = 0.75d0*dsqrt(pi4)
            bz(i,j) =-dsqrt(pi4) 
            pr(i,j) = 0.1d0
         endif
C Ryu & Jones 1a
cc         if(x(i) .le. x0) then
c         if(y(j) .le. y0) then
c            ro(i,j) = 1.d0
c            vx(i,j) = 10.d0
c            vy(i,j) = 0.d0
c            vz(i,j) = 0.d0
c            bx(i,j) = 5.d0
c            by(i,j) = 5.d0
c            bz(i,j) = 0.d0
c            pr(i,j) = 20.d0
c         else
c            ro(i,j) = 1.d0
c            vx(i,j) =-10.d0
c            vy(i,j) = 0.d0
c            vz(i,j) = 0.d0
c            bx(i,j) = 5.d0
c            by(i,j) = 5.d0
c            bz(i,j) = 0.d0
c            pr(i,j) = 1.d0
c         endif
C Ryu & Jones 2a
cc         if(x(i) .le. x0) then
c         if(y(j) .le. y0) then
c            ro(i,j) = 1.08d0
c            vx(i,j) = 1.2d0
c            vy(i,j) = 0.01d0
c            vz(i,j) = 0.5d0
c            bx(i,j) = 2.d0
c            by(i,j) = 3.6d0
c            bz(i,j) = 2.d0
c            pr(i,j) = 0.95d0
c         else
c            ro(i,j) = 1.d0
c            vx(i,j) = 0.d0
c            vy(i,j) = 0.d0
c            vz(i,j) = 0.d0
c            bx(i,j) = 2.d0
c            by(i,j) = 4.d0
c            bz(i,j) = 2.d0
c            pr(i,j) = 1.d0
c         endif
      enddo
      enddo

      call chck(ro,dmin,1,i0,1,j0)
      call chck(pr,pmin,1,i0,1,j0)

C Calculate temperature [(gm - 1)*C_p*T]
      do j=1,j0
      do i=1,i0
         te(i,j) = gm*pr(i,j)/ro(i,j)
      enddo
      enddo

C Preserve the initial values
      do j=1,j0
      do i=1,i0
         ro0(i,j) = ro(i,j)
         vx0(i,j) = vx(i,j)
         vy0(i,j) = vy(i,j)
         vz0(i,j) = vz(i,j)
         bx0(i,j) = bx(i,j)
         by0(i,j) = by(i,j)
         bz0(i,j) = bz(i,j)
         pr0(i,j) = pr(i,j)
         te0(i,j) = te(i,j)
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
     &            + pr(i,j)/(gm - 1.d0)
     &            + 0.5d0*pi4i*(bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
      enddo
      enddo

      return
      end subroutine init
