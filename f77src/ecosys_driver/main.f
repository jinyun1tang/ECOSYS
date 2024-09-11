      PROGRAM main
C     THIS SUBROUTINE READS THE RUNSCRIPT AND ENTERS FILENAMES INTO DATA ARRAYS
C     FOR USE IN 'READS' AND 'READQ'. WHEN FINISHED THIS SUBROUTINE CALLS
C     'SOIL' WHICH IS THE MAIN SUBROUTINE FROM WHICH ALL OTHERS ARE CALLED
C
      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"
      DIMENSION NA(250),ND(250)
      CHARACTER*16 DATA(30),DATAC(30,250,250),DATAP(JP,JY,JX)
     2,DATAM(JP,JY,JX),DATAX(JP),DATAY(JP),DATAZ(JP,JY,JX)
     3,OUTS(10),OUTP(10),OUTFILS(10,JY,JX),OUTFILP(10,JP,JY,JX)
      CHARACTER*3 CHOICE(102,20)
      CHARACTER*8 CDATE
      CHARACTER*80 BUF,PREFIX
      CHARACTER*100 LINE
      CALL GETCWD(BUF)
C
C     IDENTIFY OPERATING SYSTEM: DOS OR UNIX
C
      IF((.NOT.(BUF(1:1).EQ.'/'.OR.BUF(1:1).EQ.'~'))
     2.AND.BUF(2:2).EQ.':')THEN
      CALL GETARG(1,BUF)
      OPEN(5,FILE=BUF,STATUS='OLD')
      CALL GETARG(2,BUF)
      CALL CHDIR(BUF)
      PREFIX='.\\'
C     make output directory
      outdir=trim(buf)//'\\outputs\\'
      ELSE
      PREFIX='./'
C     make output directory
      outdir='./outputs/'
      ENDIF
      call system('mkdir -p '//trim(outdir))
C
C     READ INPUT FILES
C
10    FORMAT(A16)
      IGO=0
C
C     NUMBER OF COLUMNS AND ROWS
C     
      READ(5,'(A)')LINE

      les=len(line)
      ll=2
      do while (les>=ll)
      if(line(ll-1:ll-1).ne.' ' .and.
     2line(ll:ll).eq.' ')nv=nv+1
      ll=ll+1
      enddo
      idispq=0
      initro=1
      isolut=1
      if(nvs==4)then
      READ(line,*)NHW,NVN,NHE,NVS,idispq
      elseif(nvs==5)then
      READ(line,*)NHW,NVN,NHE,NVS,idispq,initro
      elseif(nvs==6)then
      READ(line,*)NHW,NVN,NHE,NVS,idispq,initro,isolut
      else
      READ(line,*)NHW,NVN,NHE,NVS
      endif
      
C
C     SITE FILE
C
      READ(5,10)DATA(1)
C
C     TOPOGRAPHY FILE
C
      READ(5,10)DATA(2)
C
C     READ THE NUMBER OF TIMES THE SCENARIOS IN THE MODEL RUN
C     ARE TO BE EXECUTED
C
100   READ(5,*)NAX,NDX
      IF(NAX.EQ.0.AND.NDX.EQ.0)GO TO 1000
C
C     NUMBER OF SCENES IN THE NEXT SCENARIO OF THE MODEL RUN
C     AND THE NUMBER OF TIMES THIS SCENARIO IS TO BE EXECUTED
C
      DO 105 NEX=1,NAX
      READ(5,*)NAY,NDY
      NA(NEX)=NAY
      ND(NEX)=NDY
C
C     FOR EACH SCENE IN THIS SCENARIO:
C
      DO 110 NE=1,NA(NEX)
C
C     WEATHER FILE
C
      READ(5,10,END=1000)DATAC(3,NE,NEX)
C
C     WEATHER OPTIONS
C
      READ(5,10)DATAC(4,NE,NEX)
C
C     LAND MANAGEMENT FILE
C
      READ(5,10)DATAC(9,NE,NEX)
C
C     PLANT MANAGEMENT FILE
C
      READ(5,10)DATAC(10,NE,NEX)
C
C     OUTPUT DATA CONTROL
C
      DO 115 N=21,30
      READ(5,10)DATAC(N,NE,NEX)
115   CONTINUE
110   CONTINUE
105   CONTINUE
C
C     RUN THIS SCRIPT
C
      DO 120 NTX=1,NDX
      DO 120 NEX=1,NAX
      DO 120 NT=1,ND(NEX)
      DO 120 NE=1,NA(NEX)
      CALL SOIL(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
      IGO=IGO+1      
120   CONTINUE
C
C     SCRIPT COMPLETED, START NEXT SCRIPT
C
      GO TO 100
1000  STOP
      END
