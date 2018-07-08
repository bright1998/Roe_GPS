c    -----------------------------------------
      subroutine numfx_roe(is,ie,js,je,nstep)
c    -----------------------------------------
      use common
      implicit double precision(a-h,o-z)

      dimension :: eigen(1:7),eps(1:7)
      dimension :: Rev(1:7,1:7)
      dimension :: anpw(1:7),anpwL(1:7),anpwR(1:7)
      integer :: is,ie,js,je,nstep

      do j=js,je
      do i=is,ie-1
         a2 = gm*qLx(i,j,8)/qLx(i,j,1)
         aas2 = (gm - 1.d0)*(qLx(i,j,10) - 0.5d0*(qLx(i,j,2)**2
     &                     + qLx(i,j,3)**2 + qLx(i,j,4)**2))
     &        - (gm - 2.d0)*pi4i*(qLx(i,j,5)**2 + qLx(i,j,6)**2
     &                          + qLx(i,j,7)**2)/qLx(i,j,1)
         if(aas2 .lt. 0.d0) aas2 = 0.d0
         valfvL(i,j) = dabs(qLx(i,j,5))/dsqrt(pi4*qLx(i,j,1))

         vf2 = 0.5d0*(aas2 + dsqrt(aas2**2 - 4.d0*a2*valfvL(i,j)**2))
         vs2 = 0.5d0*(aas2 - dsqrt(aas2**2 - 4.d0*a2*valfvL(i,j)**2))

         if(vf2 .lt. 0.d0) vf2 = 0.d0
         if(vs2 .lt. 0.d0) vs2 = 0.d0

         vfastL(i,j) = dsqrt(vf2)
         vslowL(i,j) = dsqrt(vs2)

C Flux interpolated from the left-hand side
         fxL(i,j,1) = qLx(i,j,1)*qLx(i,j,2)
         fxL(i,j,2) = qLx(i,j,1)*qLx(i,j,2)**2 + qLx(i,j,8)
     &              + 0.5d0*pi4i*(qLx(i,j,6)**2 + qLx(i,j,7)**2
     &                          - qLx(i,j,5)**2)
         fxL(i,j,3) = qLx(i,j,1)*qLx(i,j,2)*qLx(i,j,3)
     &              - pi4i*qLx(i,j,5)*qLx(i,j,6)
         fxL(i,j,4) = qLx(i,j,1)*qLx(i,j,2)*qLx(i,j,4)
     &              - pi4i*qLx(i,j,5)*qLx(i,j,7)
         fxL(i,j,5) = qLx(i,j,2)*qLx(i,j,6) - qLx(i,j,3)*qLx(i,j,5)
         fxL(i,j,6) = qLx(i,j,2)*qLx(i,j,7) - qLx(i,j,4)*qLx(i,j,5)
         fxL(i,j,7) = qLx(i,j,1)*qLx(i,j,10)*qLx(i,j,2)
     &              - pi4i*qLx(i,j,5)*
     &               (qLx(i,j,5)*qLx(i,j,2) + qLx(i,j,6)*qLx(i,j,3)
     &              + qLx(i,j,7)*qLx(i,j,4))

         a2 = gm*qRx(i,j,8)/qRx(i,j,1)
         aas2 = (gm - 1.d0)*(qRx(i,j,10) - 0.5d0*(qRx(i,j,2)**2
     &                     + qRx(i,j,3)**2 + qRx(i,j,4)**2))
     &        - (gm - 2.d0)*pi4i*(qRx(i,j,5)**2 + qRx(i,j,6)**2
     &                          + qRx(i,j,7)**2)/qRx(i,j,1)
         if(aas2 .lt. 0.d0) aas2 = 0.d0
         valfvR(i,j) = dabs(qRx(i,j,5))/dsqrt(pi4*qRx(i,j,1))

         vf2 = 0.5d0*(aas2 + dsqrt(aas2**2 - 4.d0*a2*valfvR(i,j)**2))
         vs2 = 0.5d0*(aas2 - dsqrt(aas2**2 - 4.d0*a2*valfvR(i,j)**2))

         if(vf2 .lt. 0.d0) vf2 = 0.d0
         if(vs2 .lt. 0.d0) vs2 = 0.d0

         vfastR(i,j) = dsqrt(vf2)
         vslowR(i,j) = dsqrt(vs2)

C Flux interpolated from the right-hand side
         fxR(i,j,1) = qRx(i,j,1)*qRx(i,j,2)
         fxR(i,j,2) = qRx(i,j,1)*qRx(i,j,2)**2 + qRx(i,j,8)
     &              + 0.5d0*pi4i*(qRx(i,j,6)**2 + qRx(i,j,7)**2
     &                          - qRx(i,j,5)**2)
         fxR(i,j,3) = qRx(i,j,1)*qRx(i,j,2)*qRx(i,j,3)
     &              - pi4i*qRx(i,j,5)*qRx(i,j,6)
         fxR(i,j,4) = qRx(i,j,1)*qRx(i,j,2)*qRx(i,j,4)
     &              - pi4i*qRx(i,j,5)*qRx(i,j,7)
         fxR(i,j,5) = qRx(i,j,2)*qRx(i,j,6) - qRx(i,j,3)*qRx(i,j,5)
         fxR(i,j,6) = qRx(i,j,2)*qRx(i,j,7) - qRx(i,j,4)*qRx(i,j,5)
         fxR(i,j,7) = qRx(i,j,1)*qRx(i,j,10)*qRx(i,j,2)
     &              - pi4i*qRx(i,j,5)*
     &               (qRx(i,j,5)*qRx(i,j,2) + qRx(i,j,6)*qRx(i,j,3)
     &              + qRx(i,j,7)*qRx(i,j,4))

      enddo
      enddo

      do j=js,je
      do i=is,ie-1
         sqro1 = dsqrt(qLx(i,j,1))
         sqro2 = dsqrt(qRx(i,j,1))
c Arithmetical mean of density
         armroi = 1.d0/(sqro1 + sqro2)
C Roe average
         avero = dsqrt(qLx(i,j,1)*qRx(i,j,1))
         avevx = (sqro1*qLx(i,j,2)  + sqro2*qRx(i,j,2))*armroi
         avevy = (sqro1*qLx(i,j,3)  + sqro2*qRx(i,j,3))*armroi
         avevz = (sqro1*qLx(i,j,4)  + sqro2*qRx(i,j,4))*armroi
         avebx = (sqro1*qRx(i,j,5)  + sqro2*qLx(i,j,5))*armroi
         aveby = (sqro1*qRx(i,j,6)  + sqro2*qLx(i,j,6))*armroi
         avebz = (sqro1*qRx(i,j,7)  + sqro2*qLx(i,j,7))*armroi
         aveei = (sqro1*qLx(i,j,9)  + sqro2*qRx(i,j,9))*armroi
         aveen = (sqro1*qLx(i,j,10) + sqro2*qRx(i,j,10))*armroi
         avepr = (sqro1*qLx(i,j,8)  + sqro2*qRx(i,j,8))*armroi

         aveke = 0.5d0*(avevx**2 + avevy**2 + avevz**2)

C Characteristic velocity
         delb2 = (gm - 2.d0)/(gm - 1.d0)*0.5d0*pi4i*
     &           ((qRx(i,j,6) - qLx(i,j,6))**2
     &          + (qRx(i,j,7) - qLx(i,j,7))**2)*armroi**2
         a2 = (gm - 1.d0)*(aveen - aveke
     &       - pi4i*(avebx**2 + aveby**2 + avebz**2)/avero
     &       - delb2)
         aas2 = (gm - 1.d0)*(aveen - aveke - delb2)
     &        - (gm - 2.d0)*pi4i*(avebx**2 + aveby**2 + avebz**2)/avero

         if(a2 .lt. 0.d0) a2 = 0.d0
         if(aas2 .lt. 0.d0) aas2 = 0.d0

c         va = dabs(avebx)/dsqrt(pi4*avero)
         va2 = avebx**2/(pi4*avero)
c         vf2 = 0.5d0*(aas2 + dsqrt(aas2**2 - 4.d0*a2*va**2))
c         vs2 = 0.5d0*(aas2 - dsqrt(aas2**2 - 4.d0*a2*va**2))
         vf2 = 0.5d0*(aas2 + dsqrt(aas2**2 - 4.d0*a2*va2))
         vs2 = 0.5d0*(aas2 - dsqrt(aas2**2 - 4.d0*a2*va2))

         if(vf2 .lt. 0.d0) vf2 = 0.d0
         if(vs2 .lt. 0.d0) vs2 = 0.d0
         if(vf2 .lt. va2) va2 = vf2
         if(vf2 .lt. a2) a2 = vf2

         a1 = dsqrt(a2)
         aas = dsqrt(aas2)
         vf = dsqrt(vf2)
         vs = dsqrt(vs2)
         va = dsqrt(va2)

C Eigen value
         eigen(1) = avevx + vf
         eigen(2) = avevx + va
         eigen(3) = avevx + vs
         eigen(4) = avevx
         eigen(5) = avevx - vs
         eigen(6) = avevx - va
         eigen(7) = avevx - vf

         eps(1) = dmax1(0.d0,dabs(eigen(1) - (qLx(i,j,2) + vfastL(i,j)))
     &                 ,dabs(qRx(i,j,2) + vfastR(i+1,j) - eigen(1)))
         if((eigen(1) .ge. 0.d0) .and. (eigen(1) .lt. eps(1))) then
            eigen(1) = 0.5d0*(eigen(1)**2/eps(1) + eps(1))
         elseif((eigen(1) .lt. 0.d0) .and. (-eigen(1) .lt. eps(1))) then
            eigen(1) =-0.5d0*(eigen(1)**2/eps(1) + eps(1))
         endif

         eps(2) = dmax1(0.d0,dabs(eigen(2) - (qLx(i,j,2) + valfvL(i,j)))
     &                 ,dabs(qRx(i,j,2) + valfvR(i,j) - eigen(2)))
         if((eigen(2) .ge. 0.d0) .and. (eigen(2) .lt. eps(2))) then
            eigen(2) = 0.5d0*(eigen(2)**2/eps(2) + eps(2))
         elseif((eigen(2) .lt. 0.d0) .and. (-eigen(2) .lt. eps(2))) then
            eigen(2) =-0.5d0*(eigen(2)**2/eps(2) + eps(2))
         endif

         eps(3) = dmax1(0.d0,dabs(eigen(3) - (qLx(i,j,2) + vslowL(i,j)))
     &                 ,dabs(qRx(i,j,2) + vslowR(i+1,j) - eigen(3)))
         if((eigen(3) .ge. 0.d0) .and. (eigen(3) .lt. eps(3))) then
            eigen(3) = 0.5d0*(eigen(3)**2/eps(3) + eps(3))
         elseif((eigen(3) .lt. 0.d0) .and. (-eigen(3) .lt. eps(3))) then
            eigen(3) =-0.5d0*(eigen(3)**2/eps(3) + eps(3))
         endif

         eps(4) = dmax1(0.d0,dabs(eigen(4) - qLx(i,j,2)),
     &                  dabs(qRx(i,j,2) - eigen(4)))
         if((eigen(4) .ge. 0.d0) .and. (eigen(4) .lt. eps(4))) then
            eigen(4) = 0.5d0*(eigen(4)**2/eps(4) + eps(4))
         elseif((eigen(4) .lt. 0.d0) .and. (-eigen(4) .lt. eps(4))) then
            eigen(4) =-0.5d0*(eigen(4)**2/eps(4) + eps(4))
         endif

         eps(5) = dmax1(0.d0,dabs(eigen(5) - (qLx(i,j,2) - vslowL(i,j)))
     &                 ,dabs(qRx(i,j,2) - vslowR(i+1,j) - eigen(5)))
         if((eigen(5) .ge. 0.d0) .and. (eigen(5) .lt. eps(5))) then
            eigen(5) = 0.5d0*(eigen(5)**2/eps(5) + eps(5))
         elseif((eigen(5) .lt. 0.d0) .and. (-eigen(5) .lt. eps(5))) then
            eigen(5) =-0.5d0*(eigen(5)**2/eps(5) + eps(5))
         endif

         eps(6) = dmax1(0.d0,dabs(eigen(6) - (qLx(i,j,2) - valfvL(i,j)))
     &                 ,dabs(qRx(i,j,2) - valfvR(i,j) - eigen(6)))
         if((eigen(6) .ge. 0.d0) .and. (eigen(6) .lt. eps(6))) then
            eigen(6) = 0.5d0*(eigen(6)**2/eps(6) + eps(6))
         elseif((eigen(6) .lt. 0.d0) .and. (-eigen(6) .lt. eps(6))) then
            eigen(6) =-0.5d0*(eigen(6)**2/eps(6) + eps(6))
         endif

         eps(7) = dmax1(0.d0,dabs(eigen(7) - (qLx(i,j,2) - vfastL(i,j)))
     &                 ,dabs(qRx(i,j,2) - vfastR(i+1,j) - eigen(7)))
         if((eigen(7) .ge. 0.d0) .and. (eigen(7) .lt. eps(7))) then
            eigen(7) = 0.5d0*(eigen(7)**2/eps(7) + eps(7))
         elseif((eigen(7) .lt. 0.d0) .and. (-eigen(7) .lt. eps(7))) then
            eigen(7) =-0.5d0*(eigen(7)**2/eps(7) + eps(7))
         endif

C Difference of physical variable
         delro = qRx(i,j,1) - qLx(i,j,1)
         delvx = qRx(i,j,2) - qLx(i,j,2)
         delvy = qRx(i,j,3) - qLx(i,j,3)
         delvz = qRx(i,j,4) - qLx(i,j,4)
         delby = qRx(i,j,6) - qLx(i,j,6)
         delbz = qRx(i,j,7) - qLx(i,j,7)
         delpr = qRx(i,j,8) - qLx(i,j,8)

         if((aveby .ne. 0.d0) .or. (avebz .ne. 0.d0)) then
            betay = aveby/dsqrt(aveby**2 + avebz**2)
            betaz = avebz/dsqrt(aveby**2 + avebz**2)
         else
C Avoid the singular point
C Hanawa's paper
            betay = 1.d0
            betaz = 0.d0
C Brio and Wu's paper
c            betay = 1.d0/dsqrt(2.d0)
c            betaz = 1.d0/dsqrt(2.d0)
         endif

         if(vf .ne. vs) then
c            alphaf = dsqrt(vf2 - va**2)/dsqrt(vf2 - vs2)
            alphaf = dsqrt(vf2 - va2)/dsqrt(vf2 - vs2)
            alphas = dsqrt(vf2 - a2)/dsqrt(vf2 - vs2)
         else
            alphaf = 1.d0
            alphas = 0.d0
C Brio and Wu's paper
c            alphaf = 1.d0
c            alphas = 1.d0
         endif

C Sign of Bx
         if(avebx .ge. 0.d0) then
            sgnbx = 1.d0
         else
            sgnbx =-1.d0
         endif

C Right eigenvector
c --- +Fast wave ---
         Rev(1,1) = alphaf
         Rev(2,1) = alphaf*(avevx + vf)
         Rev(3,1) = alphaf*avevy - alphas*betay*va*sgnbx
         Rev(4,1) = alphaf*avevz - alphas*betaz*va*sgnbx
         Rev(5,1) = alphas*betay*vf*dsqrt(pi4/avero)
         Rev(6,1) = alphas*betaz*vf*dsqrt(pi4/avero)
         Rev(7,1) = alphaf*(aveke + delb2 + vf*avevx + vf2/(gm - 1.d0)
     &            + (gm - 2.d0)/(gm - 1.d0)*(vf2 - a2))
     &            - alphas*va*sgnbx*(betay*avevy + betaz*avevz)
c --- +Alfven wave ---
         Rev(1,2) = 0.d0
         Rev(2,2) = 0.d0
         Rev(3,2) =-betaz*sgnbx
         Rev(4,2) = betay*sgnbx
         Rev(5,2) = betaz*dsqrt(pi4/avero)
         Rev(6,2) =-betay*dsqrt(pi4/avero)
         Rev(7,2) =-(betaz*avevy - betay*avevz)*sgnbx

c --- +Slow wave ---
         Rev(1,3) = alphas
         Rev(2,3) = alphas*(avevx + vs)
         Rev(3,3) = alphas*avevy + alphaf*betay*a1*sgnbx
         Rev(4,3) = alphas*avevz + alphaf*betaz*a1*sgnbx
         Rev(5,3) =-alphaf*betay*a2/vf*dsqrt(pi4/avero)
         Rev(6,3) =-alphaf*betaz*a2/vf*dsqrt(pi4/avero) 
         Rev(7,3) = alphas*(aveke + delb2 + vs*avevx + vs2/(gm - 1.d0)
     &            + (gm - 2.d0)/(gm - 1.d0)*(vs2 - a2))
     &            + alphaf*a1*sgnbx*(betay*avevy + betaz*avevz)

c --- Entropy wave ---
         Rev(1,4) = 1.d0
         Rev(2,4) = avevx
         Rev(3,4) = avevy
         Rev(4,4) = avevz
         Rev(5,4) = 0.d0
         Rev(6,4) = 0.d0
         Rev(7,4) = aveke + delb2

c --- -Slow wave ---
         Rev(1,5) = alphas
         Rev(2,5) = alphas*(avevx - vs)
         Rev(3,5) = alphas*avevy - alphaf*betay*a1*sgnbx
         Rev(4,5) = alphas*avevz - alphaf*betaz*a1*sgnbx
         Rev(5,5) =-alphaf*betay*a2/vf*dsqrt(pi4/avero)
         Rev(6,5) =-alphaf*betaz*a2/vf*dsqrt(pi4/avero)
         Rev(7,5) = alphas*(aveke + delb2 - vs*avevx + vs2/(gm - 1.d0)
     &            + (gm - 2.d0)/(gm - 1.d0)*(vs2 - a2))
     &            - alphaf*a1*sgnbx*(betay*avevy + betaz*avevz)

c --- -Alfven wave ---
         Rev(1,6) = 0.d0
         Rev(2,6) = 0.d0
         Rev(3,6) = betaz*sgnbx
         Rev(4,6) =-betay*sgnbx
         Rev(5,6) = betaz*dsqrt(pi4/avero)
         Rev(6,6) =-betay*dsqrt(pi4/avero)
         Rev(7,6) = (betaz*avevy - betay*avevz)*sgnbx

c --- -Fast wave ---
         Rev(1,7) = alphaf
         Rev(2,7) = alphaf*(avevx - vf)
         Rev(3,7) = alphaf*avevy + alphas*betay*va*sgnbx
         Rev(4,7) = alphaf*avevz + alphas*betaz*va*sgnbx
         Rev(5,7) = alphas*betay*vf*dsqrt(pi4/avero)
         Rev(6,7) = alphas*betaz*vf*dsqrt(pi4/avero)
         Rev(7,7) = alphaf*(aveke + delb2 - vf*avevx + vf2/(gm - 1.d0)
     &            + (gm - 2.d0)/(gm - 1.d0)*(vf2 - a2))
     &            + alphas*va*sgnbx*(betay*avevy + betaz*avevz)

         anpw1p7 = alphaf/vf2*(delpr + pi4i*(aveby*delby + avebz*delbz))
     &           + (alphas/a2/vf*((gm - 1.d0)*vs2 - (gm - 2.d0)*a2)
     &            *dsqrt(pi4*avero) + (gm - 2.d0)
     &            *dsqrt(aveby**2 + avebz**2)*alphaf/vf2)*pi4i
     &           *(betay*delby + betaz*delbz)
         anpw1m7 = alphaf/vf*avero*delvx - alphas*vs/vf/a1*sgnbx*avero
     &           *(betay*delvy + betaz*delvz)
         anpw3p5 = alphas/a2*(delpr + pi4i*(aveby*delby + avebz*delbz))
     &           + (alphaf*((gm - 2.d0)/vf - (gm - 1.d0)*vf/a2)
     &            *dsqrt(pi4*avero) + (gm - 2.d0)
     &            *dsqrt(aveby**2 + avebz**2)*alphas/a2)*pi4i
     &           *(betay*delby + betaz*delbz)
         anpw3m5 = alphas*va/vf/a1*avero*delvx
     &           + alphaf/a1*sgnbx*avero*(betay*delvy + betaz*delvz)
         anpw(1) = 0.5d0*(anpw1p7 + anpw1m7)
         anpw(2) = 0.5d0*(-avero*(betaz*delvy - betay*delvz)*sgnbx
     &           + dsqrt(pi4i*avero)*(betaz*delby - betay*delbz))
         anpw(3) = 0.5d0*(anpw3p5 + anpw3m5)
         anpw(4) = delro - alphaf*(anpw1p7) - alphas*(anpw3p5)
         anpw(5) = 0.5d0*(anpw3p5 - anpw3m5)
         anpw(6) = 0.5d0*( avero*(betaz*delvy - betay*delvz)*sgnbx
     &           + dsqrt(pi4i*avero)*(betaz*delby - betay*delbz))
         anpw(7) = 0.5d0*(anpw1p7 - anpw1m7)

C Numerical Flux
         fx(i,j,1) = 0.5d0*(fxL(i,j,1) + fxR(i,j,1))
         fx(i,j,2) = 0.5d0*(fxL(i,j,2) + fxR(i,j,2))
         fx(i,j,3) = 0.5d0*(fxL(i,j,3) + fxR(i,j,3))
         fx(i,j,4) = 0.5d0*(fxL(i,j,4) + fxR(i,j,4))
         fx(i,j,5) = 0.5d0*(fxL(i,j,5) + fxR(i,j,5))
         fx(i,j,6) = 0.5d0*(fxL(i,j,6) + fxR(i,j,6))
         fx(i,j,7) = 0.5d0*(fxL(i,j,7) + fxR(i,j,7))

         do l=1,7
            fx(i,j,1) = fx(i,j,1)
     &                - 0.5d0*dabs(eigen(l))*anpw(l)*Rev(1,l)
            fx(i,j,2) = fx(i,j,2)
     &                - 0.5d0*dabs(eigen(l))*anpw(l)*Rev(2,l)
            fx(i,j,3) = fx(i,j,3)
     &                - 0.5d0*dabs(eigen(l))*anpw(l)*Rev(3,l)
            fx(i,j,4) = fx(i,j,4)
     &                - 0.5d0*dabs(eigen(l))*anpw(l)*Rev(4,l)
            fx(i,j,5) = fx(i,j,5)
     &                - 0.5d0*dabs(eigen(l))*anpw(l)*Rev(5,l)
            fx(i,j,6) = fx(i,j,6)
     &                - 0.5d0*dabs(eigen(l))*anpw(l)*Rev(6,l)
            fx(i,j,7) = fx(i,j,7)
     &                - 0.5d0*dabs(eigen(l))*anpw(l)*Rev(7,l)
         enddo

         if(nstep .eq. 2) then
            delro = qRx(i+1,j,1) - qLx(i+1,j,1)
            delvx = qRx(i+1,j,2) - qLx(i+1,j,2)
            delvy = qRx(i+1,j,3) - qLx(i+1,j,3)
            delvz = qRx(i+1,j,4) - qLx(i+1,j,4)
            delby = qRx(i+1,j,6) - qLx(i+1,j,6)
            delbz = qRx(i+1,j,7) - qLx(i+1,j,7)
            delpr = qRx(i+1,j,8) - qLx(i+1,j,8)

            anpw1p7 = alphaf/vf2*(delpr + pi4i*(aveby*delby
     &              + avebz*delbz)) + (alphas/a2/vf*((gm - 1.d0)*vs2
     &              - (gm - 2.d0)*a2)*dsqrt(pi4*avero) + (gm - 2.d0)
     &               *dsqrt(aveby**2 + avebz**2)*alphaf/vf2)*pi4i
     &              *(betay*delby + betaz*delbz)
            anpw1m7 = alphaf/vf*avero*delvx - alphas*vs/vf/a1*sgnbx
     &               *avero*(betay*delvy + betaz*delvz)
            anpw3p5 = alphas/a2*(delpr + pi4i*(aveby*delby
     &              + avebz*delbz)) + (alphaf*((gm - 2.d0)/vf
     &              - (gm - 1.d0)*vf/a2)*dsqrt(pi4*avero) + (gm - 2.d0)
     &               *dsqrt(aveby**2 + avebz**2)*alphas/a2)*pi4i
     &              *(betay*delby + betaz*delbz)
            anpw3m5 = alphas*va/vf/a1*avero*delvx
     &              + alphaf/a1*sgnbx*avero*(betay*delvy + betaz*delvz)
            anpwR(1) = 0.5d0*(anpw1p7 + anpw1m7)
            anpwR(2) = 0.5d0*(-avero*(betaz*delvy - betay*delvz)*sgnbx
     &               + dsqrt(pi4i*avero)*(betaz*delby - betay*delbz))
            anpwR(3) = 0.5d0*(anpw3p5 + anpw3m5)
            anpwR(4) = delro - alphaf*(anpw1p7) - alphas*(anpw3p5)
            anpwR(5) = 0.5d0*(anpw3p5 - anpw3m5)
            anpwR(6) = 0.5d0*( avero*(betaz*delvy - betay*delvz)*sgnbx
     &               + dsqrt(pi4i*avero)*(betaz*delby - betay*delbz))
            anpwR(7) = 0.5d0*(anpw1p7 - anpw1m7)

            delro = qRx(i-1,j,1) - qLx(i-1,j,1)
            delvx = qRx(i-1,j,2) - qLx(i-1,j,2)
            delvy = qRx(i-1,j,3) - qLx(i-1,j,3)
            delvz = qRx(i-1,j,4) - qLx(i-1,j,4)
            delby = qRx(i-1,j,6) - qLx(i-1,j,6)
            delbz = qRx(i-1,j,7) - qLx(i-1,j,7)
            delpr = qRx(i-1,j,8) - qLx(i-1,j,8)

            anpw1p7 = alphaf/vf2*(delpr + pi4i*(aveby*delby
     &              + avebz*delbz)) + (alphas/a2/vf*((gm - 1.d0)*vs2
     &              - (gm - 2.d0)*a2)*dsqrt(pi4*avero) + (gm - 2.d0)
     &               *dsqrt(aveby**2 + avebz**2)*alphaf/vf2)*pi4i
     &              *(betay*delby + betaz*delbz)
            anpw1m7 = alphaf/vf*avero*delvx - alphas*vs/vf/a1*sgnbx
     &               *avero*(betay*delvy + betaz*delvz)
            anpw3p5 = alphas/a2*(delpr + pi4i*(aveby*delby
     &              + avebz*delbz)) + (alphaf*((gm - 2.d0)/vf
     &              - (gm - 1.d0)*vf/a2)*dsqrt(pi4*avero) + (gm - 2.d0)
     &               *dsqrt(aveby**2 + avebz**2)*alphas/a2)*pi4i
     &              *(betay*delby + betaz*delbz)
            anpw3m5 = alphas*va/vf/a1*avero*delvx
     &              + alphaf/a1*sgnbx*avero*(betay*delvy + betaz*delvz)
            anpwL(1) = 0.5d0*(anpw1p7 + anpw1m7)
            anpwL(2) = 0.5d0*(-avero*(betaz*delvy - betay*delvz)*sgnbx
     &               + dsqrt(pi4i*avero)*(betaz*delby - betay*delbz))
            anpwL(3) = 0.5d0*(anpw3p5 + anpw3m5)
            anpwL(4) = delro - alphaf*(anpw1p7) - alphas*(anpw3p5)
            anpwL(5) = 0.5d0*(anpw3p5 - anpw3m5)
            anpwL(6) = 0.5d0*( avero*(betaz*delvy - betay*delvz)*sgnbx
     &               + dsqrt(pi4i*avero)*(betaz*delby - betay*delbz))
            anpwL(7) = 0.5d0*(anpw1p7 - anpw1m7)

            do l=1,7
               fx(i,j,1) = fx(i,j,1)
     &                   + 0.25d0*(dabs(eigen(l)) - eigen(l))
     &                   *Rev(1,l)*fminmod(anpwR(l),anpw(l))
               fx(i,j,2) = fx(i,j,2)
     &                   + 0.25d0*(dabs(eigen(l)) - eigen(l))
     &                   *Rev(2,l)*fminmod(anpwR(l),anpw(l))
               fx(i,j,3) = fx(i,j,3)
     &                   + 0.25d0*(dabs(eigen(l)) - eigen(l))
     &                   *Rev(3,l)*fminmod(anpwR(l),anpw(l))
               fx(i,j,4) = fx(i,j,4)
     &                   + 0.25d0*(dabs(eigen(l)) - eigen(l))
     &                   *Rev(4,l)*fminmod(anpwR(l),anpw(l))
               fx(i,j,5) = fx(i,j,5)
     &                   + 0.25d0*(dabs(eigen(l)) - eigen(l))
     &                   *Rev(5,l)*fminmod(anpwR(l),anpw(l))
               fx(i,j,6) = fx(i,j,6)
     &                   + 0.25d0*(dabs(eigen(l)) - eigen(l))
     &                   *Rev(6,l)*fminmod(anpwR(l),anpw(l))
               fx(i,j,7) = fx(i,j,7)
     &                   + 0.25d0*(dabs(eigen(l)) - eigen(l))
     &                   *Rev(7,l)*fminmod(anpwR(l),anpw(l))
            enddo

            do l=1,7
               fx(i,j,1) = fx(i,j,1)
     &                   + 0.25d0*(dabs(eigen(l)) + eigen(l))
     &                   *Rev(1,l)*fminmod(anpwL(l),anpw(l))
               fx(i,j,2) = fx(i,j,2)
     &                   + 0.25d0*(dabs(eigen(l)) + eigen(l))
     &                   *Rev(2,l)*fminmod(anpwL(l),anpw(l))
               fx(i,j,3) = fx(i,j,3)
     &                   + 0.25d0*(dabs(eigen(l)) + eigen(l))
     &                   *Rev(3,l)*fminmod(anpwL(l),anpw(l))
               fx(i,j,4) = fx(i,j,4)
     &                   + 0.25d0*(dabs(eigen(l)) + eigen(l))
     &                   *Rev(4,l)*fminmod(anpwL(l),anpw(l))
               fx(i,j,5) = fx(i,j,5)
     &                   + 0.25d0*(dabs(eigen(l)) + eigen(l))
     &                   *Rev(5,l)*fminmod(anpwL(l),anpw(l))
               fx(i,j,6) = fx(i,j,6)
     &                   + 0.25d0*(dabs(eigen(l)) + eigen(l))
     &                   *Rev(6,l)*fminmod(anpwL(l),anpw(l))
               fx(i,j,7) = fx(i,j,7)
     &                   + 0.25d0*(dabs(eigen(l)) + eigen(l))
     &                   *Rev(7,l)*fminmod(anpwL(l),anpw(l))
            enddo
         endif
      enddo
      enddo

      return
      end subroutine numfx_roe
