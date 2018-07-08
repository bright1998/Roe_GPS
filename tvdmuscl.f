c    -----------------------------
      function fminmod(del1,del2)
c    -----------------------------
      implicit double precision(a-h,o-z)

      fminmod = dmax1(0.0d0,dmin1(del2*dsign(1.d0,del1),dabs(del1)))
     &         *dsign(1.d0,del1)
c      if(del1*del2 .le. 0.d0) then
c         fminmod = 0.d0
c      elseif(del1/del2 .le. 1.d0) then
c         fminmod = del1/del2
c      else
c         fminmod = 1.d0
c      endif

      return
      end function fminmod
