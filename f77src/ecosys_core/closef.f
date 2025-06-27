      subroutine closef

      integer, parameter :: iunit1=312
      integer, parameter :: iunit2=313
      integer, parameter :: iunit3=314
      logical :: isopen

      INQUIRE(UNIT=iunit1, OPENED=ISOPEN)
      if(isopen)close(iunit1)

      INQUIRE(UNIT=iunit2, OPENED=ISOPEN)
      if(isopen)close(iunit2)

      INQUIRE(UNIT=iunit3, OPENED=ISOPEN)
      if(isopen)close(iunit3)

      end subroutine closef
