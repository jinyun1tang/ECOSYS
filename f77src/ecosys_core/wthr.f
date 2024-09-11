      SUBROUTINE wthr(I,J,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE REINITIALIZES WEATHER VARIABLES USED IN OTHER
C     SUBROUTINES
C
      include "parameters.h"
      include "blkc.h"
      include "blk1g.h"
      include "blk2a.h"
      include "blk2b.h"
      include "blk2c.h"
      include "blk3.h"
      include "blk5.h"
      include "blk6.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk9a.h"
      include "blk9b.h"
      include "blk9c.h"
      include "blk10.h"
      include "blk11a.h"
      include "blk11b.h"
      include "blk13a.h"
      include "blk13b.h"
      include "blk13c.h"
      include "blk15a.h"
      include "blk15b.h"
      include "blk16.h"
      include "blk18a.h"
      include "blk18b.h"
      include "blk19a.h"
      include "blk19b.h"
      include "blk19c.h"
      include "blk19d.h"
      include "blk20a.h"
      include "blk20b.h"
      include "blk20c.h"
      include "blk20d.h"
      include "blk20e.h"
      include "blk20f.h"
      include "blk21a.h"
      include "blk21b.h"
      DIMENSION PRECRI(JY,JX),PRECWI(JY,JX),PRECII(JY,JX),PRECUI(JY,JX)
     2,RADN(JY,JX),VPS(JY,JX)
C
C     CDIR,CDIF=fraction of solar SW,sky diffuse radiation in visible
C     PDIR,PDIF=PAR:SW ratio (umol s-1/(MJ h-1)) 
C     TSNOW=temperature below which precipitation is snow (oC)
C
      PARAMETER (CDIR=0.42,CDIF=0.58,PDIR=1269.4,PDIF=1269.4)
      PARAMETER (TSNOW=-0.25,TWILGT=0.06976)
      XJ=J
      DOY=I-1+XJ/24
C
C     SWITCH OUT ECOSYS WEATHER HERE IF CESM WEATHER IS READ IN
C
C     IF 
C
C     IWTHR=weather data type in first(1) or second(2) scene
C     ITYPE=weather data type:1=daily,2=hourly
C
      IF(IGO.EQ.0.OR.I.LE.ILAST)THEN
      ITYPE=IWTHR(1)
      ELSE
      ITYPE=IWTHR(2)
      ENDIF
C
C     CALCULATE HOURLY TEMPERATURE, RADIATION, WINDSPEED, VAPOR PRESSURE
C     AND PRECIPITATION FROM DAILY WEATHER ARRAYS LOADED IN 'READS'
C
      IF(ITYPE.EQ.1)THEN
      DO 9915 NX=NHW,NHE
      DO 9910 NY=NVN,NVS
C
C     IETYP=Koppen climate zone:-2=phytotron
C     RADN=hourky SW radiation 
C     RMAX=maximum hourly radiation
C     ZNOON=time of solar noon
C     DYLN=daylength
C
      IF(IETYP(NY,NX).NE.-2)THEN
      IF(DYLN(NY,NX).GT.ZERO)THEN
      RADN(NY,NX)=AMAX1(0.0,RMAX*SIN((J-(ZNOON(NY,NX)
     2-DYLN(NY,NX)/2.0))*3.1416/DYLN(NY,NX)))
      ELSE
      RADN(NY,NX)=0.0
      ENDIF
      ELSE
      RADN(NY,NX)=RMAX/24.0
      ENDIF
C
C     TCA,TKA=air temperature (oC,K)
C     TAVG*,AMP*=daily averages, amplitudes from day.f
C
      IF(J.LT.(ZNOON(NY,NX)-DYLN(NY,NX)/2))THEN
      TCA(NY,NX)=TAVG1+AMP1*SIN(((J+ZNOON(NY,NX)-3.0)*3.1416
     2/(ZNOON(NY,NX)+9.0-DYLN(NY,NX)/2.0))+1.5708)
      ELSEIF(J.GT.ZNOON(NY,NX)+3)THEN
      TCA(NY,NX)=TAVG3+AMP3*SIN(((J-ZNOON(NY,NX)-3.0)*3.1416
     2/(ZNOON(NY,NX)+9.0-DYLN(NY,NX)/2.0))+1.5708)
      ELSE
      TCA(NY,NX)=TAVG2+AMP2*SIN(((J-(ZNOON(NY,NX)
     2-DYLN(NY,NX)/2.0))*3.1416/(3.0+DYLN(NY,NX)/2.0))-1.5708)
      ENDIF
      TKA(NY,NX)=TCA(NY,NX)+273.15
C
C     VPK,VPS=ambient,saturated vapor pressure
C     VAVG*,VAMP*=daily averages, amplitudes from day.f
C     ALTI=altitude
C
      IF(J.LT.(ZNOON(NY,NX)-DYLN(NY,NX)/2))THEN
      VPK(NY,NX)=VAVG1+VMP1*SIN(((J+ZNOON(NY,NX)-3.0)*3.1416
     2/(ZNOON(NY,NX)+9.0-DYLN(NY,NX)/2.0))+1.5708)
      ELSEIF(J.GT.ZNOON(NY,NX)+3)THEN
      VPK(NY,NX)=VAVG3+VMP3*SIN(((J-ZNOON(NY,NX)-3.0)*3.1416
     2/(ZNOON(NY,NX)+9.0-DYLN(NY,NX)/2.0))+1.5708)
      ELSE
      VPK(NY,NX)=VAVG2+VMP2*SIN(((J-(ZNOON(NY,NX)
     2-DYLN(NY,NX)/2.0))*3.1416 /(3.0+DYLN(NY,NX)/2.0))-1.5708)
      ENDIF
      VPS(NY,NX)=0.61*EXP(5360.0*(3.661E-03-1.0/TKA(NY,NX)))
     2*EXP(-ALTI(NY,NX)/7272.0)
      VPK(NY,NX)=AMIN1(VPS(NY,NX),VPK(NY,NX))
C
C     UA=wind speed
C
      UA(NY,NX)=AMAX1(3600.0,WIND(I))
C
C     TSNOW=temperature below which precipitation is snow (oC)
C     PRECRI,PRECWI=rainfall,snowfall
C
      IF(J.GE.13.AND.J.LE.24)THEN
      IF(TCA(NY,NX).GT.TSNOW)THEN
      PRECRI(NY,NX)=RAIN(I)/12.0
      PRECWI(NY,NX)=0.0
      ELSE
      PRECRI(NY,NX)=0.0
      PRECWI(NY,NX)=RAIN(I)/12.0
C     IF(PRECWI(NY,NX).LT.0.25E-03)PRECWI(NY,NX)=0.0
      ENDIF
      ELSE
      PRECRI(NY,NX)=0.0
      PRECWI(NY,NX)=0.0
      ENDIF
9910  CONTINUE
9915  CONTINUE
C
C     CALCULATE HOURLY TEMPERATURE, RADIATION, WINDSPEED, VAPOR PRESSURE
C     AND PRECIPITATION FROM HOURLY WEATHER ARRAYS LOADED IN 'READS'
C
      ELSE
      DO 9935 NX=NHW,NHE
      DO 9930 NY=NVN,NVS
C
C     RADN=SW radiation at horizontal surface
C     TCA,TKA=air temperature (oC,K)
C     VPK,VPS=ambient,saturated vapor pressure
C     UA=wind speed
C     TSNOW=temperature below which precipitation is snow (oC)
C     PRECRI,PRECWI=rainfall,snowfall
C
      
      RADN(NY,NX)=SRADH(J,I)
      TCA(NY,NX)=TMPH(J,I)
      TKA(NY,NX)=TCA(NY,NX)+273.15
      VPS(NY,NX)=0.61*EXP(5360.0*(3.661E-03-1.0/TKA(NY,NX)))
     2*EXP(-ALTI(NY,NX)/7272.0)
      VPK(NY,NX)=AMIN1(DWPTH(J,I),VPS(NY,NX))
      UA(NY,NX)=AMAX1(3600.0,WINDH(J,I))
      IF(TCA(NY,NX).GT.TSNOW)THEN
      PRECRI(NY,NX)=RAINH(J,I)
      PRECWI(NY,NX)=0.0
      ELSE
      PRECRI(NY,NX)=0.0
      PRECWI(NY,NX)=RAINH(J,I)
      ENDIF
9930  CONTINUE
9935  CONTINUE
      ENDIF
C
C     CALCULATE DIRECT, DIFFUSE AND LONGWAVE RADIATION FROM
C     INCOMING RADIATION READ IN 'READS', SOLAR ANGLE, HUMIDITY,
C     TEMPERATURE AND CLOUDINESS
C
      DO 9965 NX=NHW,NHE
      DO 9960 NY=NVN,NVS
C
C     IF OUTDOORS
C
C     SSIN,SSINN=sine solar angle of current,next hour
C     RADX=solar constant at horizontal surface
C     RADN=SW radiation at horizontal surface
C
      IF(IETYP(NY,NX).GE.-1)THEN
      SSIN(NY,NX)=AMAX1(0.0,AZI+DEC*COS(.2618*(ZNOON(NY,NX)-(J-0.5))))
      SSINN(NY,NX)=AMAX1(0.0,AZI+DEC*COS(.2618*(ZNOON(NY,NX)-(J+0.5))))
C     IF(SSIN(NY,NX).GT.0.0.AND.SSIN(NY,NX).LT.TWILGT)SSIN(NY,NX)=TWILGT
      IF(RADN(NY,NX).LE.0.0)SSIN(NY,NX)=0.0
      IF(SSIN(NY,NX).LE.-TWILGT)RADN(NY,NX)=0.0
      RADX=4.896*AMAX1(0.0,SSIN(NY,NX))
      RADN(NY,NX)=AMIN1(RADX,RADN(NY,NX))
C
C     DIRECT VS DIFFUSE RADIATION IN SOLAR OR SKY BEAMS
C
C     RADZ=diffuse radiation at horizontal surface
C     RADS,RADY,RAPS,RAPY=direct,diffuse SW,PAR in solar beam 
C
      RADZ=AMIN1(RADN(NY,NX),0.5*(RADX-RADN(NY,NX)))
      RADS(NY,NX)=(RADN(NY,NX)-RADZ)/SSIN(NY,NX)
      RADS(NY,NX)=AMIN1(4.167,RADS(NY,NX))
      RADY(NY,NX)=RADZ/TYSIN
      RAPS(NY,NX)=RADS(NY,NX)*CDIR*PDIR
      RAPY(NY,NX)=RADY(NY,NX)*CDIF*PDIF
C
C     ATMOSPHERIC RADIATIVE PROPERTIES AFM 139:171
C
C     CLD=cloudiness factor for EMM
C     EMM=sky emissivity
C     VPK,TKA=vapor pressure,temperature     
C
      IF(RADX.GT.ZERO)THEN
      CLD=AMIN1(1.0,AMAX1(0.2,2.33-3.33*RADN(NY,NX)/RADX))
      ELSE
      CLD=0.2
      ENDIF
      EMM=0.625*AMAX1(1.0,(1.0E+03*VPK(NY,NX)/TKA(NY,NX))**0.131)
      EMM=EMM*(1.0+0.242*CLD**0.583) 
C
C     IF PHYTOTRON
C
      ELSE
      IF(RADN(NY,NX).LE.0.0)THEN
      SSIN(NY,NX)=0.0
      ELSE
      SSIN(NY,NX)=1.0
      ENDIF
      SSINN(NY,NX)=1.0
      CLD=0.0
      EMM=0.96
      ENDIF
C
C     LONGWAVE RADIATION
C
C     XRADH=longwave radiation
C     THSX=longwave radiation from weather file or calculated from
C     atmospheric properties  
C
      IF(XRADH(J,I).GT.0.0)THEN
C     THSX(NY,NX)=EMM*(2.04E-10*TKA(NY,NX)**4)
C     THSX(NY,NX)=THSX(NY,NX)+XRADH(J,I)
      THSX(NY,NX)=XRADH(J,I)
      ELSE
      THSX(NY,NX)=EMM*(2.04E-10*TKA(NY,NX)**4)
      ENDIF
C
C     INSERT CESM WEATHER HERE
C
C     ELSE
C     RADS=DIRECT SW RADIATION (MJ M-2 H-1)
C     RADY=INDIRECT SW RADIATION (MJ M-2 H-1)
C     RAPS=DIRECT PAR (UMOL M-2 S-1)
C     RAPY=INDIRECT PAR (UMOL M-2 S-1)
C     THSX=LW RADIATION (MJ M-2 H-1)      
C     TCA=AIR TEMPERATURE (C)
C     TKA=AIR TEMPERATURE (K)
C     VPK=VAPOR PRESSURE (KPA)
C     UA=WINDSPEED (M H-1)
C     PRECRI(NY,NX)=RAIN (M H-1)
C     PRECWI(NY,NX)=SNOW (M H-1)
C     SSIN=SOLAR ANGLE CURRENT HOUR (SINE)
C     SSINN=SOLAR ANGLE NEXT HOUR (SINE)
C     ENDIF
C
C     ADD IRRIGATION
C
C     PRECII,PRECUI=surface,subsurface irrigation
C     RRIG=irrigation from soil management file in reads.f
C
      WDPTHD=WDPTH(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
C     IF(WDPTHD.LE.CDPTH(NU(NY,NX),NY,NX))THEN
      PRECII(NY,NX)=RRIG(J,I,NY,NX)
      PRECUI(NY,NX)=0.0
C     ELSE
C     PRECII(NY,NX)=0.0
C     PRECUI(NY,NX)=RRIG(J,I,NY,NX)
C     ENDIF
9960  CONTINUE
9965  CONTINUE
C
C     IMPLEMENT CLIMATE CHANGES READ IN 'READS' TO HOURLY TEMPERATURE,
C     RADIATION, WINDSPEED,VAPOR PRESSURE, PRECIPITATION, IRRIGATION
C     AND CO2
C
C     ICLM=changes to weather data (0=none,1=step,2=transient)
C     N=season or month
C
      IF(ICLM.EQ.1.OR.ICLM.EQ.2)THEN
C
C     SEASONAL CHANGES
C
      IF(I.GT.334.OR.I.LE.59)THEN
      N=1
      ELSEIF(I.GT.59.AND.I.LE.151)THEN
      N=2
      ELSEIF(I.GT.151.AND.I.LE.243)THEN
      N=3
      ELSE
      N=4
      ENDIF
C
C     MONTHLY CHANGES
C
C     IF(I.GT.0.AND.I.LE.31)THEN
C     N=1
C     ELSEIF(I.GT.31.AND.I.LE.59)THEN
C     N=2
C     ELSEIF(I.GT.59.AND.I.LE.90)THEN
C     N=3
C     ELSEIF(I.GT.90.AND.I.LE.120)THEN
C     N=4
C     ELSEIF(I.GT.120.AND.I.LE.151)THEN
C     N=5
C     ELSEIF(I.GT.151.AND.I.LE.181)THEN
C     N=6
C     ELSEIF(I.GT.181.AND.I.LE.212)THEN
C     N=7
C     ELSEIF(I.GT.212.AND.I.LE.243)THEN
C     N=8
C     ELSEIF(I.GT.243.AND.I.LE.273)THEN
C     N=9
C     ELSEIF(I.GT.273.AND.I.LE.304)THEN
C     N=10
C     ELSEIF(I.GT.304.AND.I.LE.334)THEN
C     N=11
C     ELSE
C     N=12
C     ENDIF
      DO 9925 NX=NHW,NHE
      DO 9920 NY=NVN,NVS
C
C     TEMPERATURE CHANGE
C
C     TDTPX,TDTPN=change in max,min temperature
C     DTA,AMP,DHR=change in daily average,amplitude of air temperature
C     DHR=diurnal effect on AMP         
C
      IF(TDTPX(NY,NX,N).NE.0.0.OR.TDTPN(NY,NX,N).NE.0.0)THEN
      DTA=0.5*(TDTPX(NY,NX,N)+TDTPN(NY,NX,N))
      AMP=0.5*(TDTPX(NY,NX,N)-TDTPN(NY,NX,N))
      DHR=SIN(0.2618*(J-(ZNOON(NY,NX)+3.0))+1.5708)
      TCA(NY,NX)=TCA(NY,NX)+DTA+AMP*DHR
      TKA(NY,NX)=TCA(NY,NX)+273.15
C
C     ACCLIMATION TO GRADUAL CLIMATE CHANGE
C
C     DTS=change in daily average soil temperature
C     ATCA,ATCS=mean annual air,soil temperature
C     OFFSET=shift in Arrhenius curve for MFT activity in nitro.f 
C     OFFST=shift in Arrhenius curve for PFT activity in uptake.f
C     ZTYP=PFT thermal adaptation zone
C     HTC=high temperature threshold for grain number loss (oC)
C     GROUPI,XTLI=node number at floral initiation,planting (maturity group) 
C
      IF(ICLM.EQ.2.AND.J.EQ.1)THEN
      DTS=0.5*DTA
      ATCA(NY,NX)=ATCAI(NY,NX)+DTA
      ATCS(NY,NX)=ATCAI(NY,NX)+DTS
      OFFSET(NY,NX)=0.33*(12.5-AMAX1(0.0,AMIN1(25.0,ATCS(NY,NX))))
      DO 9900 NZ=1,NP(NY,NX)
      ZTYP(NZ,NY,NX)=ZTYPI(NZ,NY,NX)+0.30/2.667*DTA
      OFFST(NZ,NY,NX)=2.667*(2.5-ZTYP(NZ,NY,NX))
C     TCZ(NZ,NY,NX)=TCZD-OFFST(NZ,NY,NX)
C     TCX(NZ,NY,NX)=AMIN1(15.0,TCZ(NZ,NY,NX)+TCXD)
      IF(ICTYP(NZ,NY,NX).EQ.3)THEN
      HTC(NZ,NY,NX)=27.0+3.0*ZTYP(NZ,NY,NX)
      ELSE
      HTC(NZ,NY,NX)=30.0+3.0*ZTYP(NZ,NY,NX)
      ENDIF
      GROUPI(NZ,NY,NX)=GROUPX(NZ,NY,NX)+0.30*DTA
      IF(IBTYP(NZ,NY,NX).NE.0)THEN
      GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)/25.0
      ENDIF
      GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)-XTLI(NZ,NY,NX)
C     IF(I.EQ.180)THEN
C     WRITE(*,1111)'OFFSET',IYRC,I,J,NZ,N,OFFSET(NY,NX),OFFST(NZ,NY,NX)
C    2,DTA,DTS,ATCA(NY,NX),ATCS(NY,NX),ZTYP(NZ,NY,NX)
C    3,GROUPI(NZ,NY,NX),TDTPX(NY,NX,N),TDTPN(NY,NX,N) 
1111  FORMAT(A8,5I4,12E12.4)
C     ENDIF
9900  CONTINUE
      ENDIF
C
C     ADJUST VAPOR PRESSURE FOR TEMPERATURE CHANGE
C
      IF(DHUM(N).EQ.1.0)THEN
      VPX=VPS(NY,NX)
      VPS(NY,NX)=0.61*EXP(5360.0*(3.661E-03-1.0/TKA(NY,NX)))
     2*EXP(-ALTI(NY,NX)/7272.0)
      VPK(NY,NX)=VPK(NY,NX)*VPS(NY,NX)/VPX
      ENDIF
      ENDIF
C
C     CHANGES IN VAPOR PRESSURE,RADIATION,WIND SPEED,PRECIPITATION
C     IRRIGATION ,CO2,NH4,NO3
C
C     TDRAD,TDWND,TDHUM=change in radiation,windspeed,vapor pressure
C     TDPRC,TDIRRI=change in precipitation,irrigation
C     TDCO2,TDCN4,TDCNO=change in atm CO2,NH4,NO3 concn in precipitation
C
      RADS(NY,NX)=RADS(NY,NX)*TDRAD(NY,NX,N)
      RADY(NY,NX)=RADY(NY,NX)*TDRAD(NY,NX,N)
      RAPS(NY,NX)=RAPS(NY,NX)*TDRAD(NY,NX,N)
      RAPY(NY,NX)=RAPY(NY,NX)*TDRAD(NY,NX,N)
      UA(NY,NX)=UA(NY,NX)*TDWND(NY,NX,N)
      VPK(NY,NX)=AMIN1(VPS(NY,NX),VPK(NY,NX)*TDHUM(NY,NX,N))
      PRECRI(NY,NX)=PRECRI(NY,NX)*TDPRC(NY,NX,N)
      PRECWI(NY,NX)=PRECWI(NY,NX)*TDPRC(NY,NX,N)
      PRECII(NY,NX)=PRECII(NY,NX)*TDIRI(NY,NX,N)
      PRECUI(NY,NX)=PRECUI(NY,NX)*TDIRI(NY,NX,N)
      CO2E(NY,NX)=CO2EI(NY,NX)*TDCO2(NY,NX,N)
      CN4R(NY,NX)=CN4RI(NY,NX)*TDCN4(NY,NX,N)
      CNOR(NY,NX)=CNORI(NY,NX)*TDCNO(NY,NX,N)
9920  CONTINUE
9925  CONTINUE
      ENDIF
C
C     DAILY WEATHER TOTALS, MAXIMA AND MINIMA FOR DAILY OUTPUT
C     CHECK AGAINST INPUTS FROM WEATHER FILE
C
      DO 9945 NX=NHW,NHE
      DO 9940 NY=NVN,NVS
      IF(SSIN(NY,NX).GT.0.0)TRAD(NY,NX)=TRAD(NY,NX)+RADS(NY,NX)
     2*SSIN(NY,NX)+RADY(NY,NX)*TYSIN
      TAMX(NY,NX)=AMAX1(TAMX(NY,NX),TCA(NY,NX))
      TAMN(NY,NX)=AMIN1(TAMN(NY,NX),TCA(NY,NX))
      HUDX(NY,NX)=AMAX1(HUDX(NY,NX),VPK(NY,NX))
      HUDN(NY,NX)=AMIN1(HUDN(NY,NX),VPK(NY,NX))
      TWIND(NY,NX)=TWIND(NY,NX)+UA(NY,NX)
      VPA(NY,NX)=VPK(NY,NX)*2.173E-03/TKA(NY,NX)
      TRAI(NY,NX)=TRAI(NY,NX)+(PRECRI(NY,NX)+PRECWI(NY,NX)
     2+PRECII(NY,NX)+PRECUI(NY,NX))*1000.0
C
C     WATER AND HEAT INPUTS TO GRID CELLS
C
C     AREA=area of grid cell (m2)
C
C     PRECR,PRECW=rainfall,snowfall
C     PRECI,PRECU=surface,subsurface irrigation
C     PRECA,PRECQ=rain+irrigation,rain+snow
C     THS=sky LW radiation
C
      PRECR(NY,NX)=PRECRI(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PRECW(NY,NX)=PRECWI(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PRECI(NY,NX)=PRECII(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PRECU(NY,NX)=PRECUI(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PRECA(NY,NX)=PRECR(NY,NX)+PRECI(NY,NX)
      PRECQ(NY,NX)=PRECR(NY,NX)+PRECW(NY,NX)
      THS(NY,NX)=THSX(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
9940  CONTINUE
9945  CONTINUE
      RETURN
      END
