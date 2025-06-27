      SUBROUTINE dspdom(I,J,NHW,NHE,NVN,NVS)
      include "parameters.h"
      include "blkc.h" 
      include "blk13a.h"
      include "blk8b.h"
      include "blk8a.h"
      integer, parameter :: iunit1=312
      integer, parameter :: iunit2=313
      integer, parameter :: iunit3=314
      character(len=36) :: filenm
      character(len=4) :: yrstr
      dimension DXY(JZ)

      if(I.EQ.1.and.J.EQ.1)THEN

      call closef

      write(yrstr,'(I4)')IYRC  
      filenm='oqc.'//yrstr      
      OPEN(UNIT=iUNIT1, FILE=filenm, STATUS='REPLACE')

      filenm='oqa.'//yrstr
      OPEN(UNIT=iUNIT2, FILE=filenm, STATUS='REPLACE')      

      filenm='mactK.'//yrstr
      OPEN(UNIT=iUNIT3, FILE=filenm, STATUS='REPLACE')      
      ENDIF

      DO  NX=NHW,NHE
      DO NY=NVN,NVS 
      DO L=1,10
      DXY(L)=AREA(3,NU(NY,NX),NY,NX)*DLYR(3,L,NY,NX)
      ENDDO
      write(iunit1,'(F7.3,I2.2,I2.2,10(X,E12.5))')I+J/24.
     2,NY,NX
     3,(SUM(OQC(0:4,L,NY,NX))/DXY(L),L=1,10)

      write(iunit2,'(F7.3,I2.2,I2.2,10(X,E12.5))')I+J/24.
     2,NY,NX
     3,(SUM(OQA(0:4,L,NY,NX))/DXY(L),L=1,10)

      write(iunit3,'(F7.3,I2.2,I2.2,10(X,F12.5))')I+J/24.
     2,NY,NX,(TOQCK(L,NY,NX),L=1,10) 
      ENDDO
      ENDDO
      END subroutine dspdom      
