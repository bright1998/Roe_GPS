c    -----------------
      subroutine prms
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      open(2,file='analysis.param',form='formatted')
      read(2,*)nsmax
      read(2,*)tmmax
      read(2,*)ftm
      read(2,*)ncflconst
c cn & dtlim are used if ncflconst = 0
      read(2,*)cn
      read(2,*)dtlim
c dtcon is used if ncflconst != 0
      read(2,*)dtcon
      read(2,*)nrestart

      open(3,file='grid.param',form='formatted')
      read(3,*)ic,jc
      read(3,*)xcmin,xcmax
      read(3,*)ycmin,ycmax
      read(3,*)i0,rax
      read(3,*)j0,ray

      open(4,file='const.param',form='formatted')
      read(4,*)gm
      read(4,*)ntcon
C Peclet Number: Nondimensional number for thermal conduction
      read(4,*)pec

      pecmin = 1.0d-10
      if(pec .lt. pecmin) pec = pecmin

      open(5,file='other.param',form='formatted')
      read(5,*)epsit
      read(5,*)itemax
C Lower Limit of Density and Pressure
      read(5,*)dmin,pmin

      open(7,file='bndc.param',form='formatted')
      read(7,*)nbnd1,nbnd2
      read(7,*)nbnd3,nbnd4

      if(ntcon .eq. 1) then
         peci = 1.d0/pec
      else 
         peci = 0.d0
      endif

      call allocate

      close(2)
      close(3)
      close(4)
      close(5)
      close(7)

      return
      end subroutine prms
