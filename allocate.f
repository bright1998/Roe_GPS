c    ---------------------
      subroutine allocate
c    ---------------------
      use common
      implicit double precision(a-h,o-z)

      allocate(x(1:i0),dx(1:i0),dxi(1:i0),
     &         y(1:j0),dy(1:j0),dyi(1:j0))

      allocate(ro(1:i0,1:j0),ro0(1:i0,1:j0),
     &         vx(1:i0,1:j0),vx0(1:i0,1:j0),
     &         vy(1:i0,1:j0),vy0(1:i0,1:j0),
     &         vz(1:i0,1:j0),vz0(1:i0,1:j0),
     &         bx(1:i0,1:j0),bx0(1:i0,1:j0),
     &         by(1:i0,1:j0),by0(1:i0,1:j0),
     &         bz(1:i0,1:j0),bz0(1:i0,1:j0),
     &         pr(1:i0,1:j0),pr0(1:i0,1:j0),
     &         te(1:i0,1:j0),te0(1:i0,1:j0))

      allocate(xm(1:i0),dxm(1:i0),dxmi(1:i0),
     &         ym(1:j0),dym(1:j0),dymi(1:j0))

      allocate( w(1:i0,1:j0,8))
      allocate(qLx(1:i0,1:j0,10),qRx(1:i0,1:j0,10),
     &         qLy(1:i0,1:j0,10),qRy(1:i0,1:j0,10))
      allocate(fx(1:i0,1:j0,7),fxL(1:i0,1:j0,7),fxR(1:i0,1:j0,7),
     &         fy(1:i0,1:j0,7),fyL(1:i0,1:j0,7),fyR(1:i0,1:j0,7))

      allocate(wh(1:i0,1:j0,8))
      allocate(roh(1:i0,1:j0),prh(1:i0,1:j0),
     &         vxh(1:i0,1:j0),vyh(1:i0,1:j0),vzh(1:i0,1:j0),
     &         bxh(1:i0,1:j0),byh(1:i0,1:j0),bzh(1:i0,1:j0))

      allocate(vfastL(1:i0,1:j0),vslowL(1:i0,1:j0),valfvL(1:i0,1:j0))
      allocate(vfastR(1:i0,1:j0),vslowR(1:i0,1:j0),valfvR(1:i0,1:j0))

      allocate(Fexx(1:i0,1:j0),Fexy(1:i0,1:j0),Fexz(1:i0,1:j0))

      do j=1,j0
      do i=1,i0
         Fexx(i,j) = 0.d0
         Fexy(i,j) = 0.d0
         Fexz(i,j) = 0.d0
      enddo
      enddo

      return
      end subroutine allocate
