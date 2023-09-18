      SUBROUTINE MBWRTR(I,J,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE REINITIALIZES HOURLY VARIABLES USED IN OTHER
C     SUBROUTINES
C     OMC(3,7,0:5,0:JZ,JY,JX)
C     1:3, labile, recalcitrant, reserve biomass
C     1:7, functional group id; 
C     0:5, five organic-microbial complexs, + autotrophic complex
C     For K=0:4, N= 1->Obligate aerobes, 2->Facultative anaerobes, 3->Fungi
C     4->Anaerobic fermenters,  5-> Acetotrophic methanogen,  6->Aerobic dizotrohpic N2 fixers
C     7->Anaerobic dizotrophic N2 fixers
C     for K=5, N= 1-> Ammonia oxidizer, 2-> Nitrite oxidizer, 3-> aerobic methanotrophs,
C     5->Hydrotrophic methanogen
C     add necessary head file to acess relevant variables

      include "parameters.h"
      include "blk8a.h"
      include "blk13a.h"

C     suppose to output methanotrophs
      D9995: DO NX=NHW,NHE
      D9990: DO NY=NVN,NVS      
      write(*,*)'methanotroph carbon',I,J,NY,NX
     2,(sum(OMC(:,3,5,L,NY,NX)),L=0,NL(NY,NX))
      write(*,*)'methanotroph nitrogen',I,J,NY,NX
     2,(sum(OMN(:,3,5,L,NY,NX)),L=0,NL(NY,NX))
      write(*,*)'methanotroph phosphorus',I,J,NY,NX
     2,(sum(OMP(:,3,5,L,NY,NX)),L=0,NL(NY,NX))

      ENDDO D9990
      ENDDO D9995

      END SUBROUTINE MBWRTR