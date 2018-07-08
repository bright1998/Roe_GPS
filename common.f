c    ---------------
      module common
c    ---------------
      implicit double precision(a-h,o-z)

c main.f
      double precision :: tm,etm,pi,pi4,pi4i
      integer :: ns,nf,nbcast

c prms.f
      integer :: nsmax,ic,i0,jc,j0,ncflconst,ntcon,itemax,nrestart
      integer :: nbnd1,nbnd2,nbnd3,nbnd4
      double precision :: xcmin,xcmax,rax
      double precision :: ycmin,ycmax,ray
      double precision :: tmmax,ftm,dtlim,cn,dtcon,gm
      double precision :: pec,peci,epsit
      double precision :: dmin,pmin

c cflc.f
      double precision :: dt

C Arraies on Normal Grids
      double precision,allocatable :: x(:),dx(:),dxi(:)
      double precision,allocatable :: y(:),dy(:),dyi(:)

      double precision,allocatable :: ro(:,:),ro0(:,:)
      double precision,allocatable :: vx(:,:),vx0(:,:)
      double precision,allocatable :: vy(:,:),vy0(:,:)
      double precision,allocatable :: vz(:,:),vz0(:,:)
      double precision,allocatable :: bx(:,:),bx0(:,:)
      double precision,allocatable :: by(:,:),by0(:,:)
      double precision,allocatable :: bz(:,:),bz0(:,:)
      double precision,allocatable :: pr(:,:),pr0(:,:)
      double precision,allocatable :: te(:,:),te0(:,:)

C Arraies on Staggered Grids
      double precision,allocatable :: xm(:),dxm(:),dxmi(:)
      double precision,allocatable :: ym(:),dym(:),dymi(:)

      double precision,allocatable ::  w(:,:,:)
      double precision,allocatable :: qLx(:,:,:),qRx(:,:,:)
      double precision,allocatable :: qLy(:,:,:),qRy(:,:,:)
      double precision,allocatable :: fx(:,:,:),fxL(:,:,:),fxR(:,:,:)
      double precision,allocatable :: fy(:,:,:),fyL(:,:,:),fyR(:,:,:)

      double precision,allocatable :: wh(:,:,:)
      double precision,allocatable :: roh(:,:),prh(:,:),
     &                                vxh(:,:),vyh(:,:),vzh(:,:),
     &                                bxh(:,:),byh(:,:),bzh(:,:)
      double precision,allocatable :: vfastL(:,:),vslowL(:,:),
     &                                valfvL(:,:)
      double precision,allocatable :: vfastR(:,:),vslowR(:,:),
     &                                valfvR(:,:)

C Source Term
      double precision,allocatable :: Fexx(:,:),Fexy(:,:),Fexz(:,:)

      end module common
