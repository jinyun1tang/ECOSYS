
      SUBROUTINE reads(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NTZ
     2,NTZX,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE READS ALL SOIL AND PLANT MANAGEMENT INPUT FILES
C
      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"
      include "blk2a.h"
      include "blk2b.h"
      include "blk2c.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk17.h"
      include "blktest.h"
      DIMENSION NA(10),ND(10)
      CHARACTER*16 DATA(30),DATAC(30,250,250),DATAP(JP,JY,JX)
     2,DATAM(JP,JY,JX),DATAX(JP),DATAY(JP),DATAZ(JP,JY,JX)
     3,OUTS(10),OUTP(10),OUTFILS(10,JY,JX),OUTFILP(10,JP,JY,JX)
      CHARACTER*3 CHOICE(102,20)
      CHARACTER*8 CDATE
      CHARACTER*16 OUTW,OUTI,OUTT,OUTN,OUTF
      CHARACTER*4 CHARY
      CHARACTER*1 TTYPE,CTYPE,IVAR(20),VAR(50),TYP(50)
      CHARACTER*80 PREFIX
      DIMENSION IDAT(20),DAT(50),DATK(50),OUT(50)
      dimension datav(50)
      PARAMETER (TWILGT=0.06976)
      DATA IFLGY,IYRX,IYRD/0,0,0/
      SAVE N1,N2,N1X,N2X,IFLGY,IYRX,IYRD
C
C     OPEN WEATHER (DATAC(3,),OPTIONS (DATAC(4,) AND LAND MANAGEMENT
C        (DATAC(9,) FILES FROM DATA ARRAYS READ IN MAIN.F
C
C     PREFIX=path for files in current or higher level directory 
C        from main.f
C
      OPEN(3,FILE=TRIM(PREFIX)//DATAC(3,NE,NEX),STATUS='OLD')
      OPEN(4,FILE=TRIM(PREFIX)//DATAC(4,NE,NEX),STATUS='OLD')
      IF(DATAC(9,NE,NEX).NE.'NO')THEN
      OPEN(13,FILE=TRIM(PREFIX)//DATAC(9,NE,NEX),STATUS='OLD')
      ENDIF
C
C     ARTIFICIAL SOIL WARMING (USED ONLY WITH SOIL WARMING EXPTS)
C
C     soiltemp=file with hourly soil temperatures copied from 
C        hourly soil temperature file in baseline run
C     OUT=hourly soil temperatures from ‘soiltemp’ (oC)
C        assuming 13 variables before temperature of soil layer 1
C     TKSZ=temperature used to calculate additional heat flux 
C        for warming in watsub.f, assuming warming of 4oC 
C
C     OPEN(6,FILE='soiltemp',STATUS='OLD')
C23   READ(6,'(F8.3,4X,A8,I8,50E16.7E3)',END=27)DOY,CDATE,J
C    2,(OUT(L),L=1,28)
C     IF(J.EQ.24)THEN
C     I=INT(DOY)
C     ELSE
C     I=INT(DOY)+1
C     ENDIF
C     DO 24 L=1,20
C     TKSZ(I,J,L)=OUT(L+13)+4.0+273.15
C     WRITE(*,2222) 'READS',I,J,L,TKSZ(I,J,L),OUT(L+13)
2222  FORMAT(A8,3I4,2E12.4)
24    CONTINUE
C     GO TO 23
27    CONTINUE
C
C     END OF ARTIFICIAL SOIL WARMING
C
C     READ START AND END DATES, WEATHER OPTIONS
C
C     IDATA(1),IDATA(2),IDATA(3)=start date of scene DDMMYYYY
C     IDATA(4),IDATA(5),IDATA(6)=end date of scene DDMMYYYY
C     IDATA(7),IDATA(8),IDATA(9)=start date of run DDMMYYYY
C        used to locate checkpoint files if run is resumed from 
C        an earlier one
C     DATA(18),DATA(19),DATA(20)=options for: 
C        visualization in visual.f (18)
C        generating checkpoint files (19)
C        resuming from earlier checkpoint files (20)
C     DRAD,DTMPX,DTMPN,DHUM,DPREC,DIRRI,DWIND,DCO2E,DCNR4,DCNOR
C        annual changes during (1)Dec-Feb,(2)Mar-May,(3)Jun-Aug,
C        (4)Sep-Nov in radiation,max+min temperature, humidity,
C        precipitation, irrigation, windspeed, atmospheric CO2
C        concentration, NH4,NO3 concentrations in precipitation 
C        from those read in the weather file
C        (temperature changes in oC, other changes relative to 
C        existing values)  
C     NPX=number of cycles per time step for water,heat,solute flux
C        calculations (suggest 12)
C     NPY=number of cycles per NPX for gas flux calculations 
C        (suggest 3)
C     JOUT,IOUT,KOUT=output frequency for hourly,daily,checkpoint data
C     ICLM=changes to weather data 
C        :0=none (weather data are used as read)
C        :1=step (simple change in selected weather data)
C        :2=transient (gradual change incremented daily)
C
      READ(4,'(2I2,I4)')IDATA(1),IDATA(2),IDATA(3)
      READ(4,'(2I2,I4)')IDATA(4),IDATA(5),IDATA(6)
      READ(4,'(2I2,I4)')IDATA(7),IDATA(8),IDATA(9)
      READ(4,'(A3)')DATA(18)
      READ(4,'(A3)')DATA(19)
      READ(4,'(A3)')DATA(20)
      DO 25 N=1,4
      READ(4,*)DRAD(N),DTMPX(N),DTMPN(N),DHUM(N),DPREC(N)
     2,DIRRI(N),DWIND(N),DCO2E(N),DCN4R(N),DCNOR(N)
25    CONTINUE
      DO 26 N=5,12
      DRAD(N)=DRAD(N-1)
      DTMPX(N)=DTMPX(N-1)
      DTMPN(N)=DTMPN(N-1)
      DHUM(N)=DHUM(N-1)
      DPREC(N)=DPREC(N-1)
      DIRRI(N)=DIRRI(N-1)
      DWIND(N)=DWIND(N-1)
      DCO2E(N)=DCO2E(N-1)
      DCN4R(N)=DCN4R(N-1)
      DCNOR(N)=DCNOR(N-1)
26    CONTINUE
      READ(4,*)NPX,NPY,JOUT,IOUT,KOUT,ICLM
C
C     INCREMENTS IN START AND END DATES FOR SUCCESSIVE SCENES
C     FROM LOOPS FOR SCENES, SCENARIOS IN RUNSCRIPT SET IN MAIN.F
C
C     IDATA(3),IDATA(6)=start,end year of current scene
C     IYRC=current year
C     IGO: =0, first scene in first scenario.
C
      NTZX=NTZ
      IF(IGO.EQ.0.OR.IDATA(3).NE.0)THEN
      IDATA(3)=IDATA(3)+(NT-1)*NF+(NTX-1)*NFX-NTZX
      IDATA(6)=IDATA(6)+(NT-1)*NF+(NTX-1)*NFX-NTZX
      IYRC=IDATA(3)
      ELSE
      IF(IDATA(1).EQ.1.AND.IDATA(2).EQ.1)THEN
      IDATA(3)=IYRC+1
      ELSE
      IDATA(3)=IYRC
      ENDIF
      IDATA(6)=IDATA(3)
      IYRC=IDATA(3)
      ENDIF
      IF(NE.EQ.1)THEN
      N1=IDATA(3)
      ENDIF
      IF(NE.EQ.NA(NEX))THEN
      N2=IDATA(6)
      NF=N2-N1+1
      IF(IDATA(4).NE.31.OR.IDATA(5).NE.12)THEN
      NTZ=NTZ+1
      ENDIF
      ENDIF
      IF(NE.EQ.1.AND.NT.EQ.1.AND.NEX.EQ.1)THEN
      N1X=IDATA(3)
      ENDIF
      IF(NE.EQ.NA(NEX).AND.NT.EQ.ND(NEX).AND.NEX.EQ.NAX)THEN
      N2X=IDATA(6)
      NFX=N2X-N1X+1
      IF(NE.NE.NA(NEX))THEN
      IF(IDATA(4).NE.31.OR.IDATA(5).NE.12)THEN
      NTZ=NTZ+1
      ENDIF
      ENDIF
      ENDIF
C     WRITE(*,7766)'IDATA3',IGO,IDATA(3),IDATA(6),IYRR,IYRC 
C    2,NE,NT,NEX,NF,NTX,NFX,NTZ,NTZX,N1,N2,N1X,N2X
C    3,NA(NEX),ND(NEX),NAX
7766  FORMAT(A8,30I8)
C
C     OPEN CHECKPOINT FILES FOR SOIL VARIABLES
C
C     IDATE=year label for checkpoint files
C     DATA(1)=site file name
C     W,N=water+heat,nutrient checkpoint files
C     DATA(20)=’YES’:find checkpoint files
C             =’NO’: start new run 
C
      IF(IGO.EQ.0)THEN
      IF(DATA(20).EQ.'YES')THEN
      IDATE=IDATA(9)
      ELSE
      IDATE=IDATA(3)
      ENDIF
      WRITE(CHARY,'(I4)')IDATE
      OUTW='W'//DATA(1)(1:2)//CHARY(1:4)
      OUTN='N'//DATA(1)(1:2)//CHARY(1:4)
      OPEN(21,FILE=OUTW,STATUS='UNKNOWN')
      OPEN(22,FILE=OUTN,STATUS='UNKNOWN')
      ENDIF
C
C     CALCULATE START AND FINISH DATES IN JULIAN DAYS
C     FROM DATE INPUTS IN OPTION FILE
C
C     ISTART,IBEGIN=start dates for current scene
C     IFIN,ILAST=end dates for current scene
C     LYRC,LYRX=number of days in current,previous year
C
      LPY=0
      LYRC=365
      LYRX=365
      DO 575 N=1,7,3
      IF(MOD(IDATA(N+2),4))520,510,520
510   IF(IDATA(N+1).GT.2)LPY=1
      IF(N.EQ.1)LYRC=366
520   IF(IDATA(N+1).EQ.1)GO TO 525
      IDY=30*(IDATA(N+1)-1)+ICOR(IDATA(N+1)-1)+IDATA(N)+LPY
      GO TO 527
525   IDY=IDATA(N)
527   IF(N.EQ.1)ISTART=IDY
      IF(N.EQ.4)IFIN=IDY
      IF(N.EQ.7)IRUN=IDY
      IF(MOD(IDATA(N+2)-1,4))575,530,575
530   IF(N.EQ.1)LYRX=366
575   CONTINUE
      IF(IGO.EQ.0)THEN
      IF(DATA(20).EQ.'NO')IRUN=ISTART
      L=1
      ILAST=ISTART-1
      ITERM=IFIN
      ELSE
      L=2
      ILAST=MIN(ISTART-1,ITERM,IEND)
      ITERM=IFIN
      ENDIF
C
C     READ WEATHER DATA
C
C     DATAC(3,=weather file
C     TTYPE=time step:S=subhourly,H=hourly,3=3-hourly,D=daily
C     CTYPE=calendar format:J=Julian,C=calendar(day,month)
C     NI,NN=number of time,weather data variables
C     IVAR,VAR=time,weather variable type
C     TYP=weather variable units
C     Z0G=windspeed measurement height
C     IFLGW=flag for adjusting Z0G for canopy height:0=no,1=yes
C     ZNOONG=time of solar noon
C     PHRG,CN4RIG,CNORIG,CPORG,CALRG,CFERG,CCARG,CMGRG,CNARG,CKARG, 
C        CSORG,CCLRG=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl 
C        concentration in precipitation (g m-3)
C     IDAT,DAT=arrays holding time,weather variables
C
      IF(DATAC(3,NE,NEX).NE.'NO')THEN
      IFLG3=0
      READ(3,'(2A1,2I2,50A1)')TTYPE,CTYPE,NI,NN,(IVAR(K),K=1,NI)
     2,(VAR(K),K=1,NN)
      READ(3,'(50A1)')(TYP(K),K=1,NN)
      read(3,*)(datav(K),K=1,3)
      Z0G=datav(1)
      IFLGW=int(datav(2))
      ZNOONG=datav(3)
      READ(3,*)PHRG,CN4RIG,CNORIG,CPORG,CALRG,CFERG,CCARG,CMGRG,CNARG
     2,CKARG,CSORG,CCLRG
      DO 55 K=1,NN
      DATK(K)=0.0
55    CONTINUE
      IH=1
60    read(3,*,END=110)(datav(k),k=1,NI),(DAT(K),K=1,NN)
      do k = 1, ni
        idat(k)=int(datav(k))
      enddo
C     WRITE(*,61)(IDAT(K),K=1,NI),(DAT(K),K=1,NN)
61    FORMAT(3I6,50E12.4)
C
C     READ DAILY WEATHER DATA AND CONVERT TO MODEL UNITS
C
      IF(TTYPE.EQ.'D')THEN
C
C     DERIVE DAY I FROM TIME VARIABLES IVAR
C
C     IWTHR=weather data type:1=daily,2=hourly for first(L=1) 
C        or second(L=2) scene
C
      IWTHR(L)=1
      DO 160 K=1,NI
      IF(IVAR(K).EQ.'M')THEN
      M=IDAT(K)
      ELSEIF(IVAR(K).EQ.'D')THEN
      N=IDAT(K)
      ENDIF
      IF(IVAR(K).EQ.'Y')THEN
      IFLGY=1
      IYRX=IDAT(K)+(NTX-1)*NFX
      IF(MOD(IDAT(K),4))170,175,170
175   IYRD=366
170   IYRD=365
      ENDIF
160   CONTINUE
      IF(IFLGY.EQ.1.AND.IYRX.LT.IYRC)GO TO 60
      IF(CTYPE.EQ.'J')THEN
      I=N
      ELSE
      LPY=0
      IF(MOD(IDATA(3),4))70,75,70
75    IF(M.GT.2)LPY=1
70    IF(M.EQ.1)THEN
      I=N
      ELSE
      I=30*(M-1)+ICOR(M-1)+N+LPY
      ENDIF
      ENDIF
C
C     DERIVE START DATE FROM TIME VARIABLES
C
      IF(IFLG3.EQ.0)THEN
      IBEGIN=I
      ISTART=MAX(ISTART,IBEGIN)
      IFLG3=1
      ENDIF
      IF(L.NE.1)THEN
      IF(I.LE.ILAST)GO TO 60
      ENDIF
C
C     CONVERT DAILY WEATHER VARIABLES TO MODEL UNITS 
C     AND ENTER INTO MODEL ARRAYS
C
C     TMPX,TMPN=maximum,minimum temperature (OC) 
C     SRAD=solar radiation (MJ m-2 d-1)
C     WIND=windspeed (m h-1)
C     DWPT=vapor pressure (kPa)
C     RAIN=precipitation (mm d-1)
C
      DO 65 K=1,NN
C
C     MAX,MIN REMPERATURE
C
      IF(VAR(K).EQ.'M')THEN
      IF(TYP(K).EQ.'F')THEN
      TMPX(I)=(DAT(K)-32.0)*0.556
      ELSEIF(TYP(K).EQ.'K')THEN
      TMPX(I)=DAT(K)-273.16
      ELSE
      TMPX(I)=DAT(K)
      ENDIF
      ELSEIF(VAR(K).EQ.'N')THEN
      IF(TYP(K).EQ.'F')THEN
      TMPN(I)=(DAT(K)-32.0)*0.556
      ELSEIF(TYP(K).EQ.'K')THEN
      TMPN(I)=DAT(K)-273.16
      ELSE
      TMPN(I)=DAT(K)
      ENDIF
C
C     SOLAR RADIATION
C
      ELSEIF(VAR(K).EQ.'R')THEN
      IF(TYP(K).EQ.'L')THEN
      SRAD(I)=AMAX1(0.0,DAT(K)/23.87)
      ELSEIF(TYP(K).EQ.'J')THEN
      SRAD(I)=AMAX1(0.0,DAT(K)*0.01)
      ELSEIF(TYP(K).EQ.'W')THEN
      SRAD(I)=AMAX1(0.0,DAT(K)*0.0864)
      ELSE
      SRAD(I)=AMAX1(0.0,DAT(K))
      ENDIF
C
C     WIND SPEED
C
      ELSEIF(VAR(K).EQ.'W')THEN
      IF(TYP(K).EQ.'S')THEN
      WIND(I)=ABS(DAT(K))*3600.0
      ELSEIF(TYP(K).EQ.'H')THEN
      WIND(I)=ABS(DAT(K))*1000.0
      ELSEIF(TYP(K).EQ.'D')THEN
      WIND(I)=ABS(DAT(K))*1000.0/24.0
      ELSEIF(TYP(K).EQ.'M')THEN
      WIND(I)=ABS(DAT(K))*1600.0
      ELSE
      WIND(I)=ABS(DAT(K))
      ENDIF
C
C     VAPOR PRESSURE
C
      ELSEIF(VAR(K).EQ.'H')THEN
      IF(TYP(K).EQ.'D')THEN
      DWPT(1,I)=0.61*EXP(5360.0*(3.661E-03-1.0
     2/(273.15+DAT(K))))
      DWPT(2,I)=0.61*EXP(5360.0*(3.661E-03-1.0
     2/(273.15+DAT(K))))
      ELSEIF(TYP(K).EQ.'F')THEN
      DAT(K)=(DAT(K)-32.0)*0.556
      DWPT(1,I)=0.61*EXP(5360.0*(3.661E-03-1.0
     2/(273.15+DAT(K))))
      DWPT(2,I)=0.61*EXP(5360.0*(3.661E-03-1.0
     2/(273.15+DAT(K))))
      ELSEIF(TYP(K).EQ.'H')THEN
      DAT(K)=AMAX1(0.0,AMIN1(1.0,DAT(K)))
      DWPT(1,I)=0.61*EXP(5360.0*(3.661E-03-1.0
     2/(273.15+(TMPN(I)+TMPX(I))/2)))*DAT(K)
      DWPT(2,I)=0.61*EXP(5360.0*(3.661E-03-1.0
     2/(273.15+TMPN(I))))
      ELSEIF(TYP(K).EQ.'R')THEN
      DAT(K)=AMAX1(0.0,AMIN1(100.0,DAT(K)))
      DWPT(1,I)=0.61*EXP(5360.0*(3.661E-03-1.0
     2/(273.15+(TMPN(I)+TMPX(I))/2)))*DAT(K)*0.01
      DWPT(2,I)=0.61*EXP(5360.0*(3.661E-03-1.0
     2/(273.15+TMPN(I))))
      ELSEIF(TYP(K).EQ.'S')THEN
      DWPT(1,I)=AMAX1(0.0,DAT(K))*0.0289/18.0*101.325
     2*EXP(-ALTIG/7272.0)*288.15/(273.15+(TMPN(I)+TMPX(I))/2)
      DWPT(2,I)=AMAX1(0.0,DAT(K))*0.0289/18.0*101.325
     2*EXP(-ALTIG/7272.0)*288.15/(273.15+TMPN(I))
      ELSEIF(TYP(K).EQ.'G')THEN
      DWPT(1,I)=AMAX1(0.0,DAT(K))*28.9/18.0*101.325
     2*EXP(-ALTIG/7272.0)*288.15/(273.15+(TMPN(I)+TMPX(I))/2)
      DWPT(2,I)=AMAX1(0.0,DAT(K))*28.9/18.0*101.325
     2*EXP(-ALTIG/7272.0)*288.15/(273.15+TMPN(I))
      ELSEIF(TYP(K).EQ.'M')THEN
      DWPT(1,I)=AMAX1(0.0,DAT(K)*0.1)
      DWPT(2,I)=AMAX1(0.0,DAT(K)*0.1)
      ELSE
      DWPT(1,I)=AMAX1(0.0,DAT(K))
      DWPT(2,I)=AMAX1(0.0,DAT(K))
      ENDIF
C
C     PRECIPITATION
C
      ELSEIF(VAR(K).EQ.'P')THEN
      IF(TYP(K).EQ.'M')THEN
      RAIN(I)=AMAX1(0.0,DAT(K))/1.0E+03
      ELSEIF(TYP(K).EQ.'C')THEN
      RAIN(I)=AMAX1(0.0,DAT(K))/1.0E+02
      ELSEIF(TYP(K).EQ.'I')THEN
      RAIN(I)=AMAX1(0.0,DAT(K))*0.0254
      ELSE
      RAIN(I)=AMAX1(0.0,DAT(K))
      ENDIF
      ENDIF
65    CONTINUE
      IX=I
      IF(IFLGY.EQ.1.AND.I.EQ.IYRD)THEN
      GO TO 110
      ENDIF
      GO TO 60
C
C     READ HOURLY WEATHER DATA AND CONVERT TO MODEL UNITS
C
      ELSE
C
C     DERIVE DAY I AND HOUR J FROM TIME VARIABLES IVAR
C
      IWTHR(L)=2
      DO 190 K=1,NI
      IF(IVAR(K).EQ.'M')THEN
      M=IDAT(K)
      ELSEIF(IVAR(K).EQ.'D')THEN
      N=IDAT(K)
      ELSEIF(IVAR(K).EQ.'H')THEN
      J=IDAT(K)
      ENDIF
      IF(IVAR(K).EQ.'Y')THEN
      IFLGY=1
      ENDIF
190   CONTINUE
      IF(IFLGY.EQ.1.AND.IYRX.LT.IYRC)GO TO 60
      IF(CTYPE.EQ.'J')THEN
      I=N
      ELSE
      LPY=0
      IF(MOD(IDATA(3),4))100,115,100
115   IF(M.GT.2)LPY=1
100   IF(M.EQ.1)THEN
      I=N
      ELSE
      I=30*(M-1)+ICOR(M-1)+N+LPY
      ENDIF
      ENDIF
      IF(J.GT.24.AND.(J/100)*100.NE.J)THEN
      DO 80 K=1,NN
      DATK(K)=DATK(K)+DAT(K)
80    CONTINUE
      IH=IH+1
      GO TO 60
      ENDIF
      IF(J.GT.24)J=INT(J/100)
      IF(J.EQ.0)THEN
      J=24
      I=I-1
      IF(I.LT.1)GO TO 60
      ENDIF
C
C     DERIVE START DATE FROM TIME VARIABLES
C
      IF(IFLG3.EQ.0)THEN
      IBEGIN=N
      ISTART=MAX(ISTART,IBEGIN)
      IFLG3=1
      ENDIF
      IF(L.NE.1)THEN
      IF(I.LE.ILAST)GO TO 60
      ENDIF
      XRADH(J,I)=0.0
C
C     CONVERT HOURLY WEATHER VARIABLES TO MODEL UNITS 
C     AND ENTER INTO MODEL ARRAYS
C
C     TMPH=temperature (oC)
C     SRADH=solar radiation (MJ m-2 h-1)
C     WINDH=windspeed (m h-1)
C     DWPTH=vapor pressure (kPa)
C     RAINH=precipitation (mm h-1)
C     XRADH=longwave radiation (MJ m-2 h-1)
C
      DO 95 K=1,NN
C
C     TEMPERATURE
C
      IF(VAR(K).EQ.'T')THEN
      IF(TYP(K).EQ.'F')THEN
      TMPH(J,I)=((DAT(K)+DATK(K))/IH-32.0)*0.556
      ELSEIF(TYP(K).EQ.'K')THEN
      TMPH(J,I)=(DAT(K)+DATK(K))/IH-273.16
      ELSE
      TMPH(J,I)=(DAT(K)+DATK(K))/IH
      ENDIF
C
C     SOLAR RADIATION
C
      ELSEIF(VAR(K).EQ.'R')THEN
      IF(TYP(K).EQ.'W')THEN
      SRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH*0.0036)
      ELSEIF(TYP(K).EQ.'J')THEN
      SRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH*0.01)
      ELSEIF(TYP(K).EQ.'K')THEN
      SRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH*0.001)
      ELSEIF(TYP(K).EQ.'P')THEN
      SRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH*0.0036*0.457)
C     ELSEIF(TYP(K).EQ.'M')THEN
C     SRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH*3.6*0.457)
      ELSE
      SRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH)
      ENDIF
C
C     WIND SPEED
C
      ELSEIF(VAR(K).EQ.'W')THEN
      IF(TYP(K).EQ.'S')THEN
      WINDH(J,I)=(DAT(K)+DATK(K))/IH*3600.0
      ELSEIF(TYP(K).EQ.'H')THEN
      WINDH(J,I)=(DAT(K)+DATK(K))/IH*1000.0
      ELSEIF(TYP(K).EQ.'M')THEN
      WINDH(J,I)=(DAT(K)+DATK(K))/IH*1600.0
      ELSE
      WINDH(J,I)=(DAT(K)+DATK(K))/IH
      ENDIF
C
C     VAPOR PRESSURE
C
      ELSEIF(VAR(K).EQ.'H')THEN
      IF(TYP(K).EQ.'D')THEN
      DWPTH(J,I)=0.61*EXP(5360.0*(3.661E-03-1.0/(273.15
     2+(DAT(K)+DATK(K))/IH)))
      ELSEIF(TYP(K).EQ.'F')THEN
      DAT(K)=(DAT(K)-32.0)*0.556
      DWPTH(J,I)=0.61*EXP(5360.0*(3.661E-03-1.0/(273.15
     2+(DAT(K)+DATK(K))/IH)))
      ELSEIF(TYP(K).EQ.'H')THEN
      DWPTH(J,I)=0.61*EXP(5360.0*(3.661E-03-1.0/(273.15+TMPH(J,I))))
     2*AMAX1(0.0,AMIN1(1.0,(DAT(K)+DATK(K))/IH))
      ELSEIF(TYP(K).EQ.'R')THEN
      DWPTH(J,I)=0.61*EXP(5360.0*(3.661E-03-1.0/(273.15+TMPH(J,I))))
     2*AMAX1(0.0,AMIN1(100.0,(DAT(K)+DATK(K))/IH))*0.01
      ELSEIF(TYP(K).EQ.'S')THEN
      DWPTH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH)*0.0289/18.0*101.325
     2*EXP(-ALTIG/7272.0)*288.15/(273.15+TMPH(J,I))
      ELSEIF(TYP(K).EQ.'G')THEN
      DWPTH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH)*28.9/18.0*101.325
     2*EXP(-ALTIG/7272.0)*288.15/(273.15+TMPH(J,I))
      ELSEIF(TYP(K).EQ.'M')THEN
      DWPTH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH*0.1)
      ELSE
      DWPTH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH)
      ENDIF
C
C     PRECIPITATION
C
      ELSEIF(VAR(K).EQ.'P')THEN
      IF(TYP(K).EQ.'M')THEN
      RAINH(J,I)=AMAX1(0.0,DAT(K)+DATK(K))/1.0E+03
      ELSEIF(TYP(K).EQ.'C')THEN
      RAINH(J,I)=AMAX1(0.0,DAT(K)+DATK(K))/1.0E+02
      ELSEIF(TYP(K).EQ.'I')THEN
      RAINH(J,I)=AMAX1(0.0,DAT(K)+DATK(K))*0.0254
      ELSEIF(TYP(K).EQ.'S')THEN
      IF(TTYPE.EQ.'H')THEN
      RAINH(J,I)=AMAX1(0.0,DAT(K)+DATK(K))*3.6
      ELSE
      RAINH(J,I)=AMAX1(0.0,DAT(K)+DATK(K))*1.8
      ENDIF
      ELSE
      RAINH(J,I)=AMAX1(0.0,DAT(K)+DATK(K))
      ENDIF
C
C     LONGWAVE RADIATION (OPTIONAL)
C
      ELSEIF(VAR(K).EQ.'L')THEN
      IF(TYP(K).EQ.'W')THEN
      XRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH*0.0036)
      ELSEIF(TYP(K).EQ.'J')THEN
      XRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH*0.01)
      ELSEIF(TYP(K).EQ.'K')THEN
      XRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH*0.001)
      ELSE
      XRADH(J,I)=AMAX1(0.0,(DAT(K)+DATK(K))/IH)
      ENDIF
      ENDIF
      DATK(K)=0.0
95    CONTINUE
      IH=1
      IX=I
C
C     INFILL 3-HOURLY WEATHER DATA
C
      IF(TTYPE.EQ.'3')THEN
      JJ=J-3
      II=I
      IF(JJ.EQ.0)THEN
      JJ=24
      II=I-1
      ENDIF
      IF(II.LT.1)THEN
      TMPH(J-2,I)=TMPH(J,I)
      TMPH(J-1,I)=TMPH(J,I)
      SRADH(J-2,I)=SRADH(J,I)
      SRADH(J-1,I)=SRADH(J,I)
      WINDH(J-2,I)=WINDH(J,I)
      WINDH(J-1,I)=WINDH(J,I)
      DWPTH(J-2,I)=DWPTH(J,I)
      DWPTH(J-1,I)=DWPTH(J,I)
      RAINH(J,I)=RAINH(J,I)/3.0
      RAINH(J-2,I)=RAINH(J,I)
      RAINH(J-1,I)=RAINH(J,I)
      XRADH(J-2,I)=XRADH(J,I)
      XRADH(J-1,I)=XRADH(J,I)
      ELSE
      TMPH(J-2,I)=0.667*TMPH(JJ,II)+0.333*TMPH(J,I)
      TMPH(J-1,I)=0.333*TMPH(JJ,II)+0.667*TMPH(J,I)
      SRADH(J-2,I)=0.667*SRADH(JJ,II)+0.333*SRADH(J,I)
      SRADH(J-1,I)=0.333*SRADH(JJ,II)+0.667*SRADH(J,I)
      WINDH(J-2,I)=0.667*WINDH(JJ,II)+0.333*WINDH(J,I)
      WINDH(J-1,I)=0.333*WINDH(JJ,II)+0.667*WINDH(J,I)
      DWPTH(J-2,I)=0.667*DWPTH(JJ,II)+0.333*DWPTH(J,I)
      DWPTH(J-1,I)=0.333*DWPTH(JJ,II)+0.667*DWPTH(J,I)
      RAINH(J,I)=RAINH(J,I)/3.0
      RAINH(J-2,I)=RAINH(J,I)
      RAINH(J-1,I)=RAINH(J,I)
      XRADH(J-2,I)=0.667*XRADH(JJ,II)+0.333*XRADH(J,I)
      XRADH(J-1,I)=0.333*XRADH(JJ,II)+0.667*XRADH(J,I)
      ENDIF
      ENDIF
      IF(IFLGY.EQ.1.AND.I.EQ.IYRD.AND.J.EQ.24)THEN
      GO TO 110
      ENDIF
      GO TO 60
      ENDIF
110   CONTINUE
C
C     ACCOUNT FOR LEAP YEAR
C
      IF(I.EQ.365)THEN
      IF(TTYPE.EQ.'D')THEN
      TMPX(I+1)=TMPX(I)
      TMPN(I+1)=TMPN(I)
      SRAD(I+1)=SRAD(I)
      WIND(I+1)=WIND(I)
      DWPT(1,I+1)=DWPT(1,I)
      DWPT(2,I+1)=DWPT(2,I)
      RAIN(I+1)=RAIN(I)
      ELSE
      DO 130 J=1,24
      TMPH(J,I+1)=TMPH(J,I)
      SRADH(J,I+1)=SRADH(J,I)
      WINDH(J,I+1)=WINDH(J,I)
      DWPTH(J,I+1)=DWPTH(J,I)
      RAINH(J,I+1)=RAINH(J,I)
      XRADH(J,I+1)=XRADH(J,I)
130   CONTINUE
      ENDIF
      IX=I+1
      ENDIF
      ELSE
      IFLGW=1
      Z0G=2.0
      ZNOONG=12.0
      PHRG=7.0
      CN4RIG=0.0
      CNORIG=0.0
      CN4RG=CN4RIG
      CNORG=CNORIG
      CPORG=0.0
      CALRG=0.0
      CFERG=0.0
      CCARG=0.0
      CMGRG=0.0
      CNARG=0.0
      CKARG=0.0
      CSORG=0.0
      CCLRG=0.0
      IX=365
      ENDIF
C
C     CALCULATE PRECIPITATION CONCENTRATIONS IN MOLE UNITS
C
      CN4RIG=CN4RIG/14.0
      CNORIG=CNORIG/14.0
      CN4RG=CN4RIG
      CNORG=CNORIG
      CPORG=CPORG/31.0
      CALRG=CALRG/27.0
      CFERG=CFERG/55.8
      CCARG=CCARG/40.0
      CMGRG=CMGRG/24.3
      CNARG=CNARG/23.0
      CKARG=CKARG/39.1
      CSORG=CSORG/32.0
      CCLRG=CCLRG/35.5
C
C     ALLOCATE PECIITATION ION CONCENTRATIONS TO LANDSCAPE UNITS
C
      DO 8970 NX=NHW,NHE
      DO 8975 NY=NVN,NVS
      Z0(NY,NX)=Z0G
      ZNOON(NY,NX)=ZNOONG
      PHR(NY,NX)=PHRG
      CN4RI(NY,NX)=CN4RIG
      CNORI(NY,NX)=CNORIG
      CN4R(NY,NX)=CN4RIG
      CNOR(NY,NX)=CNORIG
      CPOR(NY,NX)=CPORG
      CALR(NY,NX)=CALRG
      CFER(NY,NX)=CFERG
      CCAR(NY,NX)=CCARG
      CMGR(NY,NX)=CMGRG
      CNAR(NY,NX)=CNARG
      CKAR(NY,NX)=CKARG
      CSOR(NY,NX)=CSORG
      CCLR(NY,NX)=CCLRG
8975  CONTINUE
8970  CONTINUE
C
C     DERIVE END DATES FROM TIME VARIABLES
C
      ICHECK=0
      IF(TTYPE.EQ.'H'.AND.J.NE.24)ICHECK=1
      IEND=IX-ICHECK
      IFIN=MIN(IFIN,IEND)
      IDAYR=MIN(ISTART-1,ILAST)
      IYRR=IDATA(3)
      NYR=0
      IF(IDAYR.EQ.0)THEN
      IDAYR=LYRX
      IYRR=IDATA(3)-1
      NYR=1
      ENDIF
      IFLGY=0
      CLOSE(3)
      CLOSE(4)
C
C     READ LAND MANAGEMENT FILE NAMES FOR EACH GRID CELL
C
C     ROWN,ROWO,ROWP=row spacing for band NH4, NO3 and PO4
C        fertilizer applications
C     ITILL,DCORP=disturbance type,intensity
C     FERT,IYTYP,FDPTH=fertilizer amount,type,application depth
C     RRIG=irrigation timing,depth
C     PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX, 
C     CSOQX,CCLQX=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl 
C        concentration in irrigation water (g m-3)
C
      DO 9980 NX=NHW,NHE
      DO 9985 NY=NVN,NVS
      ROWN(NY,NX)=0.0
      ROWO(NY,NX)=0.0
      ROWP(NY,NX)=0.0
      DO 325 I=1,366
      ITILL(I,NY,NX)=0
      DCORP(I,NY,NX)=0.0
325   CONTINUE
      DO 40 I=1,366
      DO 45 N=1,20
      FERT(N,I,NY,NX)=0.0
45    CONTINUE
      DO 35 N=0,2
      IYTYP(N,I,NY,NX)=0
35    CONTINUE
      FDPTH(I,NY,NX)=0.0
40    CONTINUE
      DO 125 I=1,366
      DO 120 J=1,24
      RRIG(J,I,NY,NX)=0.0
120   CONTINUE
      PHQ(I,NY,NX)=7.0
      CN4Q(I,NY,NX)=0.0
      CNOQ(I,NY,NX)=0.0
      CPOQ(I,NY,NX)=0.0
      CALQ(I,NY,NX)=0.0
      CFEQ(I,NY,NX)=0.0
      CCAQ(I,NY,NX)=0.0
      CMGQ(I,NY,NX)=0.0
      CNAQ(I,NY,NX)=0.0
      CKAQ(I,NY,NX)=0.0
      CSOQ(I,NY,NX)=0.0
      CCLQ(I,NY,NX)=0.0
      WDPTH(I,NY,NX)=0.0
      ROWI(I,NY,NX)=0.0
125   CONTINUE
9985  CONTINUE
9980  CONTINUE
C
C     READ LAND MANAGEMENT FILE DATAC(9 LOADED IN 'MAIN'.
C     THIS FILE CONTAINS NAMES OF DISTURBANCE DATA(8),
C     FERTILIZER DATA(5) AND IRRIGATION DATA(6) FILES
C
      IF(DATAC(9,NE,NEX).NE.'NO')THEN
C
C     NH1,NV1,NH2,NV2=N,W and S,E corners of landscape unit
C     DATA(8),DATA(5),DATA(6)=disturbance,fertilizer,irrigation files
C     PREFIX=path for files in current or higher level directory
C
150   READ(13,*,END=200)NH1,NV1,NH2,NV2
      READ(13,*)DATA(8),DATA(5),DATA(6)
      IF(DATA(8).NE.'NO')THEN
      OPEN(10,FILE=TRIM(PREFIX)//DATA(8),STATUS='OLD')
      ENDIF
      IF(DATA(5).NE.'NO')THEN
      OPEN(8,FILE=TRIM(PREFIX)//DATA(5),STATUS='OLD')
      ENDIF
      IF(DATA(6).NE.'NO')THEN
      OPEN(2,FILE=TRIM(PREFIX)//DATA(6),STATUS='OLD')
      ENDIF
C
C     READ TILLAGE INPUT FILE
C
      IF(DATA(8).NE.'NO')THEN
C
C     DY=date DDMMYYYY
C     IDIST=soil disturbance type 1-10:tillage including crop
C                                 11-20:tillage not including crop 
C                                 21:surface litter removal
C                                 22:fire
C                                 23-24:natural, artificial drainage
C     DDIST=energy (kW m-2) (fire, likely 25 - 75) 
C           or depth (m) (tillage,drainage) of disturbance 
C     RCHGNAG,RCHGEAG,RCHGSAG,RCHGWAG=distance to N,E,S,W external
C        artificial water table (m)
C     RCHGNBG,RCHGEBG,RCHGSBG,RCHGWBG=boundary conditions for 
C        N,E,S,W subsurface exchange with external artificial water 
C        table :0=no flow, 1=unimpeded flow 
C
295   CONTINUE
      READ(10,*,END=305)DY,IDIST,DDIST
      IF(IDIST.EQ.24)THEN
      READ(10,*)RCHGNAG,RCHGEAG,RCHGSAG,RCHGWAG
     2,RCHGNBG,RCHGEBG,RCHGSBG,RCHGWBG
      ENDIF
      LPY=0
      IDY1=INT(DY/1.0E+06)
      IDY2=INT(DY/1.0E+04-IDY1*1.0E+02)
      IDY3=INT(DY-(IDY1*1.0E+06+IDY2*1.0E+04))
      IF(MOD(IDY3,4))3520,3510,3520
3510  IF(IDY2.GT.2)LPY=1
3520  IF(IDY2.EQ.1)GO TO 3535
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
C     IF(IDY2.LE.6)IDY=IDY-0.5*(NTX-1)
C     IF(IDY2.GE.7)IDY=IDY+0.5*(NTX-1)
      GO TO 3530
3535  IDY=IDY1
3530  CONTINUE
C
C     ENTER DISTURBANCE PROPERTIES TO GRID CELLS 
C     IN EACH LANDSCAPE UNIT
C
      DO 8995 NX=NH1,NH2
      DO 8990 NY=NV1,NV2
      ITILL(IDY,NY,NX)=IDIST
      DCORP(IDY,NY,NX)=DDIST
      IF(IDIST.EQ.24)THEN
      RCHGNA(NY,NX)=RCHGNAG
      RCHGEA(NY,NX)=RCHGEAG
      RCHGSA(NY,NX)=RCHGSAG
      RCHGWA(NY,NX)=RCHGWAG
      RCHGNB(NY,NX)=RCHGNBG
      RCHGEB(NY,NX)=RCHGEBG
      RCHGSB(NY,NX)=RCHGSBG
      RCHGWB(NY,NX)=RCHGWBG
      ENDIF
8990  CONTINUE
8995  CONTINUE
      GO TO 295
305   CONTINUE
      CLOSE(10)
      ENDIF
C
C     READ FERTLIZER INPUT FILE
C
      IF(DATA(5).NE.'NO')THEN
C
C     DY=date DDMMYYYY
C     *A,*B=broadcast,banded fertilizer application
C     Z4,Z3,ZU,ZO=NH4,NH3,urea,NO3
C     PM*,PH*=Ca(H2PO4)2,apatite
C     CAC,CAS=CaCO3,CaSO4
C     *1,*2=litter,manure amendments
C     RSC,RSN,RSC=amendment C,N,P content
C     FDPTHI=fertilizer application depth
C     ROWX=row width of band application
C     IRO,IR1,IR2=fertilizer,litter,manure type
C        IR0:0=fast (urea as urine),1=normal,2=slow release
C        IR1,IR2:used to set decomposition parameters in ‘hour1.f’ 
C
1500  CONTINUE
      READ(8,*,END=85)DY,Z4A,Z3A,ZUA,ZOA,Z4B,Z3B,ZUB,ZOB
     2,PMA,PMB,PHA,CAC,CAS,RSC1,RSN1,RSP1,RSC2,RSN2,RSP2,FDPTHI
     3,ROWX,IR0,IR1,IR2
      LPY=0
      IDY1=INT(DY/1.0E+06)
      IDY2=INT(DY/1.0E+04-IDY1*1.0E+02)
      IDY3=INT(DY-(IDY1*1.0E+06+IDY2*1.0E+04))
      IF(MOD(IDY3,4))1520,1510,1520
1510  IF(IDY2.GT.2)LPY=1
1520  IF(IDY2.EQ.1)GO TO 1525
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
C     IF(IDY2.LE.6)IDY=IDY-0.5*(NTX-1)
C     IF(IDY2.GE.7)IDY=IDY+0.5*(NTX-1)
      GO TO 1530
1525  IDY=IDY1
1530  CONTINUE
      DO 8985 NX=NH1,NH2
      DO 8980 NY=NV1,NV2
C
C     ENTER AMENDMENTS INTO GRID CELLS FOR EACH LANDSCAPE UNIT 
C
C     NH4,NH3,UREA,NO3 BROADCAST (A) AND BANDED (B)
C
      FERT(1,IDY,NY,NX)=Z4A
      FERT(2,IDY,NY,NX)=Z3A
      FERT(3,IDY,NY,NX)=ZUA
      FERT(4,IDY,NY,NX)=ZOA
      FERT(5,IDY,NY,NX)=Z4B
      FERT(6,IDY,NY,NX)=Z3B
      FERT(7,IDY,NY,NX)=ZUB
      FERT(8,IDY,NY,NX)=ZOB
C
C     MONOCALCIUM PHOSPHATE OR HYDROXYAPATITE BROADCAST (A)
C     AND BANDED (B)
C
      FERT(9,IDY,NY,NX)=PMA
      FERT(10,IDY,NY,NX)=PMB
      FERT(11,IDY,NY,NX)=PHA
C
C     LIME AND GYPSUM
C
      FERT(12,IDY,NY,NX)=CAC
      FERT(13,IDY,NY,NX)=CAS
C
C     PLANT AND ANIMAL RESIDUE C, N AND P
C
      FERT(14,IDY,NY,NX)=RSC1
      FERT(15,IDY,NY,NX)=RSN1
      FERT(16,IDY,NY,NX)=RSP1
      FERT(17,IDY,NY,NX)=RSC2
      FERT(18,IDY,NY,NX)=RSN2
      FERT(19,IDY,NY,NX)=RSP2
C
C     DEPTH AND WIDTH OF APPLICATION
C
      FDPTH(IDY,NY,NX)=FDPTHI
      ROWI(IDY,NY,NX)=ROWX
C
C     TYPE OF FERTILIZER,PLANT OR ANIMAL RESIDUE
C
      IYTYP(0,IDY,NY,NX)=IR0
      IYTYP(1,IDY,NY,NX)=IR1
      IYTYP(2,IDY,NY,NX)=IR2
8980  CONTINUE
8985  CONTINUE
      GO TO 1500
85    CONTINUE
      CLOSE(8)
      ENDIF
C
C     READ IRRIGATION INPUT FILE
C
      IF(DATA(6).NE.'NO')THEN
C
C     AUTOMATED IRRIGATION
C
C     auto=first 4 characters of the irrigation file name
C
      IF(DATA(6)(1:4).EQ.'auto')THEN
C
C     DST,DEN=start,end dates,hours DDMMHHHH
C     IFLGVX=flag for irrigation criterion,0=SWC,1=canopy water potl 
C     FIRRX=remaining SWC calculated from CIRRX to WP(IFLGV=0) 
C        or minimum canopy water potential(IFLGV=1), to trigger
C        irrigation
C     CIRRX= fraction of FC to which irrigation will raise SWC
C     DIRRX= depth to which water depletion and rewatering is
C        calculated
C     WDPTHI=depth at which irrigation is applied (0=surface) 
C     PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX, 
C        CSOQX,CCLQX=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl 
C        concentration in irrigation water (g m-3)
C
      READ(2,*,END=105)DST,DEN,IFLGVX,FIRRX,CIRRX,DIRRX,WDPTHI
     2,PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX
     3,CSOQX,CCLQX
      LPY=0
      IDY1=INT(DST/1.0E+06)
      IDY2=INT(DST/1.0E+04-IDY1*1.0E+02)
      IDY3=INT(DST-(IDY1*1.0E+06+IDY2*1.0E+04))
      IF(MOD(IDY3,4))4520,4510,4520
4510  IF(IDY2.GT.2)LPY=1
4520  IF(IDY2.EQ.1)GO TO 4535
      IDYS=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
      GO TO 4530
4535  IDYS=IDY1
4530  CONTINUE
      IHRS=IDY3
      LPY=0
      IDY1=INT(DEN/1.0E+06)
      IDY2=INT(DEN/1.0E+04-IDY1*1.0E+02)
      IDY3=INT(DEN-(IDY1*1.0E+06+IDY2*1.0E+04))
      IF(MOD(IDY3,4))5520,5510,5520
5510  IF(IDY2.GT.2)LPY=1
5520  IF(IDY2.EQ.1)GO TO 5535
      IDYE=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
      GO TO 5530
5535  IDYE=IDY1
5530  CONTINUE
      IHRE=IDY3
C
C     TRANSFER AUTOMATED IRRIGATION INPUTS TO MODEL ARRAYS
C     FOR EACH LANDSCAPE UNIT 
C
      DO 7965 NX=NH1,NH2
      DO 7960 NY=NV1,NV2
      IFLGV(NY,NX)=IFLGVX
      IIRRA(1,NY,NX)=IDYS
      IIRRA(2,NY,NX)=IDYE
      IIRRA(3,NY,NX)=INT(IHRS/100)
      IIRRA(4,NY,NX)=INT(IHRE/100)
      FIRRA(NY,NX)=FIRRX
      CIRRA(NY,NX)=CIRRX
      DIRRA(1,NY,NX)=DIRRX
      DIRRA(2,NY,NX)=WDPTHI
C
C     EXPRESS ION CONCENTRATIONS IN MOLAR UNITS
C
      DO 220 I=1,366
      PHQ(IDY,NY,NX)=PHQX
      CN4Q(IDY,NY,NX)=CN4QX/14.0
      CNOQ(IDY,NY,NX)=CNOQX/14.0
      CPOQ(IDY,NY,NX)=CPOQX/31.0
      CALQ(IDY,NY,NX)=CALQX/27.0
      CFEQ(IDY,NY,NX)=CFEQX/55.8
      CCAQ(IDY,NY,NX)=CCAQX/40.0
      CMGQ(IDY,NY,NX)=CMGQX/24.3
      CNAQ(IDY,NY,NX)=CNAQX/23.0
      CKAQ(IDY,NY,NX)=CKAQX/39.1
      CSOQ(IDY,NY,NX)=CSOQX/32.0
      CCLQ(IDY,NY,NX)=CCLQX/35.5
220   CONTINUE
7960  CONTINUE
7965  CONTINUE
      ELSE
C
C     SCHEDULED IRRIGATION
C
2500  CONTINUE
C
C     DY,RR,JST,JEN=date DDMMYYYY,amount (mm),start and end hours
C     WDPTHI=depth at which irrigation is applied (0=surface) 
C     PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX, 
C        CSOQX,CCLQX=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl 
C        concentration in irrigation water
C
      READ(2,*,END=105)DY,RR,JST,JEN,WDPTHI,PHQX,CN4QX,CNOQX,CPOQX
     2,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,CSOQX,CCLQX
      LPY=0
      IDY1=INT(DY/1.0E+06)
      IDY2=INT(DY/1.0E+04-IDY1*1.0E+02)
      IDY3=INT(DY-(IDY1*1.0E+06+IDY2*1.0E+04))
      IF(MOD(IDY3,4))2520,2510,2520
2510  IF(IDY2.GT.2)LPY=1
2520  IF(IDY2.EQ.1)GO TO 2525
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
      GO TO 2530
2525  IDY=IDY1
2530  CONTINUE
C
C     ENTER IRRIGATION AMOUNTS INTO MODEL ARRAY
C     FOR EACH LANDSCAPE UNIT
C
      RRH=RR/(JEN-(JST-1))
      DO 8965 NX=NH1,NH2
      DO 8960 NY=NV1,NV2
      DO 2535 J=1,24
      IF(J.GE.JST.AND.J.LE.JEN)RRIG(J,IDY,NY,NX)=RRH/1000.0
2535  CONTINUE
C
C     TRANSFER ION CONCENTRATIONS IN MOLAR UNITS TO MODEL ARRAYS 
C     FOR EACH LANDSCAPE UNIT
C
      PHQ(IDY,NY,NX)=PHQX
      CN4Q(IDY,NY,NX)=CN4QX/14.0
      CNOQ(IDY,NY,NX)=CNOQX/14.0
      CPOQ(IDY,NY,NX)=CPOQX/31.0
      CALQ(IDY,NY,NX)=CALQX/27.0
      CFEQ(IDY,NY,NX)=CFEQX/55.8
      CCAQ(IDY,NY,NX)=CCAQX/40.0
      CMGQ(IDY,NY,NX)=CMGQX/24.3
      CNAQ(IDY,NY,NX)=CNAQX/23.0
      CKAQ(IDY,NY,NX)=CKAQX/39.1
      CSOQ(IDY,NY,NX)=CSOQX/32.0
      CCLQ(IDY,NY,NX)=CCLQX/35.5
      WDPTH(IDY,NY,NX)=WDPTHI
8960  CONTINUE
8965  CONTINUE
      GO TO 2500
      ENDIF
105   CONTINUE
      ENDIF
      CLOSE(2)
      GO TO 150
200   CONTINUE
      CLOSE(13)
      ENDIF
      IMNG=1
      RETURN
      END


