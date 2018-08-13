      subroutine endrun(msg, line)

      implicit none

      character(len=*), intent(in) :: msg
      integer, intent(in) :: line
      write(*,'(A,A,i0)')msg,' at line ',line
      stop
      end subroutine endrun

