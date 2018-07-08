c    -------------------
      subroutine source
c    -------------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine gradient(ain,aout1,aout2,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: ain
            double precision,dimension(:,:),intent(out) :: aout1,aout2
            integer :: is,ie,js,je
         end subroutine gradient
      end interface

C Block for Allocating Working Array --- Please add new array freely ---
      dimension :: warray0(1:i0,1:j0)

c#################################################################
c### You can add a steady source term (external force).        ###
c### If it is unnecessary, you need not edit this source file. ###
c#################################################################

c Fexx(i,j)  : the source term of the momentum equation in the x-direction.
c Fexy(i,j)  :                                          in the y-direction.
c Fexz(i,j)  :                                          in the z-direction.
c x(i)       : the x-coordinate.
c y(j)       : the y-coordinate.
c warray(i,j): working array, which you can use freely.

c Start of Editable Block -------------------------------------------
c###################################################################
c If you want to set up the potential like the gravitational energy, 
c please use the array "warray".
c###################################################################
      do j=1,j0
      do i=1,i0
         warray0(i,j) = 0.d0
      enddo
      enddo
c
c##################################################
c Need not calculate the gradient of the potential.
c Please use the subroutine "gradient".
c##################################################
c
      call gradient(warray0,Fexx,Fexy,1,i0,1,j0)
c
c Don't forget to set up the boundary condition.
c Example
      do j=1,j0
         Fexx(1,j) = Fexx(2,j)
         Fexx(i0,j) = Fexx(i0-1,j)
      enddo
      do i=1,i0
         Fexy(i,1) = Fexy(i,2)
         Fexy(i,j0) = Fexy(i,j0-1)
      enddo

c      do j=1,j0
c      do i=1,i0
c         Fexx(i,j) = 
c         Fexy(i,j) = 
c         Fexz(i,j) =
c      enddo
c      enddo
c End of Editable Block

      return
      end subroutine source
