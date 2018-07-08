c    ---------------------------------------
      subroutine chck(da,damin,is,ie,js,je)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: damin
      integer :: is,ie,js,je

      do j=js,je
      do i=is,ie
         if(da(i,j) .lt. damin) then
            da(i,j) = damin
         endif
      enddo
      enddo

      return
      end subroutine chck
