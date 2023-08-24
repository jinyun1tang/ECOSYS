      subroutine omwriter(I,J,NHW,NHE,NVN,NVS,option)

      include "files.h"
      include "parameters.h"
      include "blk13a.h"
      include "blkc.h"
      include "blk8a.h"

C      1st-dimension, protein(1),nonstructural(2)
C        cellulose(3), lignin(4), charcol(5)
C      2nd-dimension, coarse woody litter (0), fine non-woody litter (1),
C      animal manure (2), POM (3), and humus (4)
C
C      OSC(5,0:4,0:JZ,JY,JX),OSA(5,0:4,0:JZ,JY,JX)
      character(len=*)    :: option
      character(len=1024) :: filename
      if(trim(option)=='openfile')then
        write(filename,'(A,I4,A)')trim(outdir),IYRC,'osc'        
        open(unit=101,FILE=trim(filename),STATUS='UNKNOWN')
      elseif(trim(option)=='closefile')then
        close(101)
      else
C     the example below is writing hourly output     
        D9995: DO NX=NHW,NHE
        D9990: DO NY=NVN,NVS
          write(101,*)'coarse_woody_litter_charcoalC',I,J,NY,NX
     2      ,OSC(5,0,0,NY,NX),(OSC(5,0,L,NY,NX),L=NU(NY,NX),NL(NY,NX))
          write(101,*)'fine_litter_charcolC',I,J,NY,NX
     2      ,OSC(5,1,0,NY,NX),(OSC(5,1,L,NY,NX),L=NU(NY,NX),NL(NY,NX))
          write(101,*)'animal_manure_charcolC',I,J,NY,NX
     2      ,OSC(5,2,0,NY,NX),(OSC(5,2,L,NY,NX),L=NU(NY,NX),NL(NY,NX))
          write(101,*)'POM_charcoalC',I,J,NY,NX
     2      ,OSC(5,3,0,NY,NX),(OSC(5,3,L,NY,NX),L=NU(NY,NX),NL(NY,NX))
          write(101,*)'Humus_charcoalC',I,J,NY,NX
     2      ,OSC(5,4,0,NY,NX),(OSC(5,4,L,NY,NX),L=NU(NY,NX),NL(NY,NX))
        ENDDO D9990
        ENDDO D9995
      endif

      end subroutine omwriter
