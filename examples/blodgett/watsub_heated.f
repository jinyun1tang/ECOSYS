      SUBROUTINE watsub(I,J,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE CACULATES ENERGY BALANCES OF SNOW, RESIDUE
C     AND SOIL SURFACES, FREEZING, THAWING, AND HEAT AND WATER
C     TRANSFER THROUGH SOIL PROFILES
C
      include "parameters.h"
      include "blkc.h"
      include "blk2a.h"
      include "blk2b.h"
      include "blk2c.h"
      include "blk5.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk10.h"
      include "blk11a.h"
      include "blk11b.h"
      include "blk13a.h"
      include "blk13b.h"
      include "blk13c.h"
      include "blk15a.h"
      include "blk15b.h"
      include "blk22a.h"
      include "blk22b.h"
      include "blk22c.h"
      include "blktest.h"
      DIMENSION VOLWX1(JZ,JY,JX),VPQ(JY,JX),TKQ(JY,JX)
     2,XVOLT(JY,JX),XVOLW(JY,JX),XVOLI(JY,JX),FMAC(JZ,JY,JX)
     3,FGRD(JZ,JY,JX),VOLW1(0:JZ,JY,JX),VOLI1(0:JZ,JY,JX)
     4,VHCP1(0:JZ,JY,JX),VHCP1A(JZ,JY,JX),VHCP1B(JZ,JY,JX)
     4,TK1(0:JZ,JY,JX),TWFLXL(JZ,JY,JX),VOLW2(JZ,JY,JX)
     5,VOLP1(0:JZ,JY,JX),TWFLXH(JZ,JY,JX),PRECM(JY,JX)
     6,VOLS0(JS,JY,JX),VOLI0(JS,JY,JX),VOLW0(JS,JY,JX)
     7,VOLS1(JS,JY,JX),DLYRS0(JS,JY,JX),VOLP1Z(JZ,JY,JX)
     8,TK0(JS,JY,JX),AREAU(JZ,JY,JX),AREAUD(JZ,JY,JX),FLQ0S(JY,JX)
     9,FLQ0I(JY,JX),FLQ0W(JY,JX),FLQ1(JY,JX),FLH1(JY,JX)
     9,FLY1(JY,JX),HWFLQ0(JY,JX),HWFLQ1(JY,JX),HWFLY1(JY,JX)
     1,RAR(JY,JX),RAGS(JY,JX),BAREW(JY,JX),CVRDW(JY,JX)
     2,RARG(JY,JX),RAGR(JY,JX),RAGW(JY,JX),PAREW(JY,JX),PAREG(JY,JX)
     3,RAG(JY,JX),PARSW(JY,JX),PARSG(JY,JX),PARER(JY,JX),PARSR(JY,JX)
     4,QR1(2,2,JV,JH),HQR1(2,2,JV,JH),VOLPH1Z(JZ,JY,JX)
     5,QS1(2,JV,JH),QW1(2,JV,JH),QI1(2,JV,JH),HQS1(2,JV,JH)
     6,TQR1(JY,JX),THQR1(JY,JX),EVAPG(JY,JX),TTFLXL(JZ,JY,JX)
     7,EVAPW(JY,JX),EVAPS(JY,JX),EVAPR(JY,JX),FLWRL(JY,JX)
     8,HFLWRL(JY,JX),FINHL(JZ,JY,JX),FLWL(3,JD,JV,JH)
      DIMENSION FLWHL(3,JD,JV,JH),HFLWL(3,JD,JV,JH)
     2,TFLWL(JZ,JY,JX),TFLWHL(JZ,JY,JX),THFLWL(JZ,JY,JX)
     3,WFLXL(JZ,JY,JX),TFLXL(JZ,JY,JX),AVCNHL(3,JD,JV,JH)
     5,THRYW(JY,JX),THRMW(JY,JX),THRMS(JY,JX),THRMR(JY,JX)
     6,THRYG(JY,JX),THRYR(JY,JX),RADXW(JY,JX),RADXG(JY,JX)
     7,RADXR(JY,JX),FLWLX(3,JD,JV,JH),TFLWLX(JZ,JY,JX)
     8,FLU1(JZ,JY,JX),HWFLU1(JZ,JY,JX),PSISM1(0:JZ,JY,JX)
     9,ALTG(JY,JX),WFLXLH(JZ,JY,JX),DLYRR(JY,JX),WFLXR(JY,JX)
     1,TFLXR(JY,JX),HCNDR(JY,JX),CNDH1(JZ,JY,JX),VOLA1(0:JZ,JY,JX)
     2,THETWX(0:JZ,JY,JX),THETIX(0:JZ,JY,JX),THETPX(0:JZ,JY,JX)
     4,VOLAH1(JZ,JY,JX),VOLWH1(JZ,JY,JX),VOLPH1(JZ,JY,JX)
     5,VOLIH1(JZ,JY,JX),THETPY(0:JZ,JY,JX),FLWNX(JY,JX)
     6,FLWXNX(JY,JX),FLWHNX(JY,JX),HFLWNX(JY,JX),N6X(JY,JX)
     7,PSISA1(JZ,JY,JX)
      DIMENSION TQS1(JY,JX),TQW1(JY,JX),TQI1(JY,JX),THQS1(JY,JX)
     2,TFLWS(JS,JY,JX),TFLWW(JS,JY,JX),TFLWI(JS,JY,JX)
     3,THFLWW(JS,JY,JX),TFLX0(JS,JY,JX),WFLXS(JS,JY,JX)
     4,WFLXI(JS,JY,JX),FLW0W(JS,JY,JX),FLW0S(JS,JY,JX)
     5,FLW0I(JS,JY,JX),HFLW0W(JS,JY,JX)
     6,VOLS0M(JS,JY,JX),VOLW0M(JS,JY,JX)
     7,VOLI0M(JS,JY,JX),VHCPWMM(JS,JY,JX),TK0M(JS,JY,JX)
C
C     EMMS,EMMW,EMMR=emissivities of surface soil, snow and litter
C     RACX,RARX=minimum boundary layer resistances of canopy,litter (h m-1)
C     RZ=minimum resistance to evaporation of surface water (h m-1)
C     RAM=minimum boundary layer resistance (h m-1)
C     DPTHSX=minimum snowpack depth for full cover (m)
C     Z1S,Z2SW,Z2SD,Z3SX=parameters for air-water gas transfers in soil
C     Z1R,Z2RW,Z2RD,Z3RX=parameters for air-water gas transfers in litter
C
      PARAMETER (EMMS=0.97,EMMW=0.97,EMMR=0.97
     2,RACX=0.0139,RARX=0.0139,RZ=0.0139,RAM=1.39-03,DPTHSX=0.075)
      PARAMETER (Z1S=0.010,Z2SW=12.0,Z2SD=12.0,Z3SX=0.50
     2,Z1R=0.01,Z2RW=12.0,Z2RD=12.0,Z3R=0.50)
C
C     Parameters for calculating convective effects on heat transfer
C     in porous media (air and water)
C     VISCW,VISCA=water,air viscosity (Mg m-1 s)
C
      PARAMETER (VISCW=1.0E-06,VISCA=2.0E-08,DIFFW=1.45E-07
     2,DIFFA=2.01E-05,EXPNW=2.07E-04,EXPNA=3.66E-03,GRAV=9.8
     3,RYLXW=GRAV*EXPNW/(VISCW*DIFFW),RYLXA=GRAV*EXPNA/(VISCA*DIFFA)
     4,PRNTW=VISCW/DIFFW,PRNTA=VISCA/DIFFA
     5,DNUSW=(1.0+(0.492/PRNTW)**0.5625)**0.4444
     6,DNUSA=(1.0+(0.492/PRNTA)**0.5625)**0.4444
     7,TRBW=0.375,TRBA=0.000)
C
C     FVOLAH=parameter for clay effect on macropore volume
C     DTHETW=difference between saturation and effective saturation
C     HCNDRR=saturated hydraulic conductivity of surface litter
C     FENGYP=rate constant for restoring surface Ksat
C
      PARAMETER (FVOLAH=0.0,DTHETW=1.0E-06,HCNDRR=25.0
     2,FENGYP=1.0E-03)
      REAL*4 RI,THETWR,THETW1,THETA1,THETAL,THETWL
     2,TKR1,TKS1,TKY,TKW1,TK11,TK12,TK0X,TKXR,TK1X,TKX1,TFND1

      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      NUM(NY,NX)=NU(NY,NX)
C
C     ADJUST SURFACE ELEVATION USED IN RUNOFF FOR FREEZE-THAW, EROSION
C     AND SOC
C
C     ALTG,ALT=current,initial elevation of ground surface
C     CDPTH(NUM(NY,NX)-1,=depth of ground surface
C     ENGYP=cumulative rainfall energy impact on soil surface
C
      ALTG(NY,NX)=ALT(NY,NX)-CDPTH(NUM(NY,NX)-1,NY,NX)
      ENGYP(NY,NX)=ENGYP(NY,NX)*(1.0-FENGYP)
C
C     ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
C     FOR USE AT INTERNAL TIME STEP
C
C     SET INITIAL SNOWPACK VALUES
C
C     VOLS0,VOLSSL=snowpack snow content (water equivalent)
C     VOLI0,VOLISSL=snowpack ice content
C     VOLW0,VOLWSL=snowpack water content
C     VOLS1,VOLSL=snowpack volume
C     DLYRS0,DLYRS=snowpack depth
C     VHCPWM,VHCPW=snowpack heat capacity
C     TK0,TKW=snowpack temperature
C
      DO 60 L=1,JS
      VOLS0(L,NY,NX)=VOLSSL(L,NY,NX)
      VOLI0(L,NY,NX)=VOLISL(L,NY,NX)
      VOLW0(L,NY,NX)=VOLWSL(L,NY,NX)
      VOLS1(L,NY,NX)=VOLSL(L,NY,NX)
      DLYRS0(L,NY,NX)=DLYRS(L,NY,NX)
      VHCPWM(1,L,NY,NX)=VHCPW(L,NY,NX)
      TK0(L,NY,NX)=TKW(L,NY,NX)
60    CONTINUE
C
C     SET INITIAL SOIL VALUES
C
C     WFLXR,TFLXR=initialize surface litter freeze,thaw,latent heat
C     CDPTH=depth to bottom of soil layer
C     WDPTH,LWDPTH=depth,layer of subsurface irrigation
C
      WFLXR(NY,NX)=0.0
      TFLXR(NY,NX)=0.0
      DO 65 L=NUM(NY,NX),NL(NY,NX)
      IF(CDPTH(L,NY,NX).GE.WDPTH(I,NY,NX))THEN
      LWDPTH=L
      GO TO 55
      ENDIF
65    CONTINUE
55    CONTINUE
      DO 30 L=NUM(NY,NX),NL(NY,NX)
C
C     ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
C     FOR USE AT INTERNAL TIME STEP IN30    con SOIL LAYERS
C
C     PSISM1,PSISM=matric water potential
C     VOLA*,VOLW*,VOLI*,VOLP*=pore,water,ice,air volumes of micropores
C     VOLWX1=VOLW1 accounting for wetting front
C     VOLAH*,VOLWH*,VOLIH*,VOLPH*=pore,water,ice,air macropores
C     BKDS=bulk density
C     CCLAY=clay concentration
C     FVOLAH=parameter for clay effect on macropore volume
C     VOLX,VOLT=soil,total volumes
C     WP=wilting point
C     THETW*,THETI*,THETP*=water,ice,air-filled porosity
C     VHCP1,VHCM=volumetric heat capacities of total volume, solid
C     VHCP1A,VHCP1B=volumetric heat capacities of micropore,macropore
C
      PSISM1(L,NY,NX)=PSISM(L,NY,NX)
      VOLA1(L,NY,NX)=VOLA(L,NY,NX)
      VOLW1(L,NY,NX)=VOLW(L,NY,NX)
      VOLWX1(L,NY,NX)=VOLWX(L,NY,NX)
      VOLI1(L,NY,NX)=VOLI(L,NY,NX)
      VOLWH1(L,NY,NX)=VOLWH(L,NY,NX)
      VOLIH1(L,NY,NX)=VOLIH(L,NY,NX)
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
C VOLA1: total volume in micropores
C VOLP1: air-filled volume
      VOLP1Z(L,NY,NX)=VOLA1(L,NY,NX)-VOLW1(L,NY,NX)-VOLI1(L,NY,NX)
      VOLP1(L,NY,NX)=AMAX1(0.0,VOLP1Z(L,NY,NX))
      ELSE
      VOLP1Z(L,NY,NX)=0.0
      VOLP1(L,NY,NX)=0.0
      ENDIF
      VOLAH1(L,NY,NX)=AMAX1(0.0,VOLAH(L,NY,NX)-FVOLAH*CCLAY(L,NY,NX)
     2*(VOLW1(L,NY,NX)/VOLY(L,NY,NX)-WP(L,NY,NX))*VOLT(L,NY,NX))
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLPH1Z(L,NY,NX)=VOLAH1(L,NY,NX)-VOLWH1(L,NY,NX)-VOLIH1(L,NY,NX)
      VOLPH1(L,NY,NX)=AMAX1(0.0,VOLPH1Z(L,NY,NX))
      ELSE
      VOLPH1Z(L,NY,NX)=0.0
      VOLPH1(L,NY,NX)=0.0
      ENDIF
      VOLWM(1,L,NY,NX)=VOLW1(L,NY,NX)
      VOLWHM(1,L,NY,NX)=VOLWH1(L,NY,NX)
      VOLPM(1,L,NY,NX)=VOLP1(L,NY,NX)+VOLPH1(L,NY,NX)
     2+THETPI*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
      VOLTX=VOLY(L,NY,NX)+VOLAH1(L,NY,NX)
      IF(VOLTX.GT.ZEROS2(NY,NX))THEN
      THETWX(L,NY,NX)=AMAX1(0.0,(VOLW1(L,NY,NX)+VOLWH1(L,NY,NX))
     2/VOLTX)
      THETIX(L,NY,NX)=AMAX1(0.0,(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
     2/VOLTX)
      THETPX(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX))
     2/VOLTX)
      ELSE
      THETWX(L,NY,NX)=POROS(L,NY,NX)
      THETIX(L,NY,NX)=0.0
      THETPX(L,NY,NX)=0.0
      ENDIF
      THETPM(1,L,NY,NX)=THETPX(L,NY,NX)
      IF(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETPY(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX))
     2/(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX)))
      ELSE
      THETPY(L,NY,NX)=0.0
      ENDIF
      VHCP1(L,NY,NX)=VHCM(L,NY,NX)+4.19*(VOLW1(L,NY,NX)
     2+VOLWH1(L,NY,NX))+1.9274*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
      VHCP1A(L,NY,NX)=VHCM(L,NY,NX)+4.19*VOLW1(L,NY,NX)
     2+1.9274*VOLI1(L,NY,NX)
      VHCP1B(L,NY,NX)=4.19*VOLWH1(L,NY,NX)+1.9274*VOLIH1(L,NY,NX)
C     IF(I.GT.331)THEN
C     WRITE(*,3376)'VOLWI',I,J,NX,NY,L,VOLW1(L,NY,NX)
C    2,VOLI1(L,NY,NX),VOLP1(L,NY,NX),VOLA1(L,NY,NX)
C    3,VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),VOLPH1(L,NY,NX)
C    3,VOLAH1(L,NY,NX),VOLT(L,NY,NX),VOLY(L,NY,NX)
C    4,THETWX(L,NY,NX),THETIX(L,NY,NX),THETPX(L,NY,NX)
3376  FORMAT(A8,5I4,40E14.6)
C     ENDIF
C
C     MACROPOROSITY
C
C     FMAC,FGRD=macropore,micropore volume fractions
C     CNDH*=macropore hydraulic conductivity
C     TKS,TK1=soil temperature
C     FLU,HWFLU=subsurface water,convective heat fluxes
C     AREAU,AREAD=fractions of layer below natural,artifl water table
C
      IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      FMAC(L,NY,NX)=FHOL(L,NY,NX)*VOLAH1(L,NY,NX)/VOLAH(L,NY,NX)
      CNDH1(L,NY,NX)=CNDH(L,NY,NX)
     2*(VOLAH1(L,NY,NX)/VOLAH(L,NY,NX))**2
      ELSE
      FMAC(L,NY,NX)=0.0
      CNDH1(L,NY,NX)=0.0
      ENDIF
      FGRD(L,NY,NX)=1.0-FMAC(L,NY,NX)
      TK1(L,NY,NX)=TKS(L,NY,NX)
      IF(L.EQ.LWDPTH)THEN
      FLU(L,NY,NX)=PRECU(NY,NX)
      HWFLU(L,NY,NX)=4.19*TKA(NY,NX)*PRECU(NY,NX)
      FLU1(L,NY,NX)=FLU(L,NY,NX)*XNPH
      HWFLU1(L,NY,NX)=HWFLU(L,NY,NX)*XNPH
      ELSE
      FLU(L,NY,NX)=0.0
      HWFLU(L,NY,NX)=0.0
      FLU1(L,NY,NX)=0.0
      HWFLU1(L,NY,NX)=0.0
      ENDIF
      IF(CDPTH(L,NY,NX).GE.DTBLX(NY,NX))THEN
      AREAU(L,NY,NX)=AMIN1(1.0,AMAX1(0.0
     2,(CDPTH(L,NY,NX)-DTBLX(NY,NX))
     2/DLYR(3,L,NY,NX)))
      ELSE
      AREAU(L,NY,NX)=0.0
      ENDIF
      IF(CDPTH(L,NY,NX).GE.DTBLY(NY,NX))THEN
      AREAUD(L,NY,NX)=AMIN1(1.0,AMAX1(0.0
     2,(CDPTH(L,NY,NX)-DTBLY(NY,NX))
     2/DLYR(3,L,NY,NX)))
      ELSE
      AREAUD(L,NY,NX)=0.0
      ENDIF
30    CONTINUE
C
C     ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
C     FOR USE AT INTERNAL TIME STEP IN SURFACE LITTER
C
C     THRMG=longwave emission from litter surface
C     VHCP1=volumetric heat capacity of litter
C     VOLA*,VOLW*,VOLI*,VOLP*=pore,water,ice,air volumes of litter
C     VOLWRX=maximum water retention by litter
C     XVOLT,XVOLW=free surface water+ice,water
C     VHCPRX=min heat capacity for litter water,heat fluxes
C     VOLR=litter volume
C     THETW*,THETI*,THETP*=water,ice,air concentrations
C     PSISM*=litter matric water potential
C
      THRMG(NY,NX)=0.0
      VHCP1(0,NY,NX)=2.496E-06*ORGC(0,NY,NX)+4.19*VOLW(0,NY,NX)
     2+1.9274*VOLI(0,NY,NX)
      VOLA1(0,NY,NX)=VOLA(0,NY,NX)
      VOLW1(0,NY,NX)=AMAX1(0.0,VOLW(0,NY,NX))
      VOLI1(0,NY,NX)=AMAX1(0.0,VOLI(0,NY,NX))
      VOLP1(0,NY,NX)=AMAX1(0.0,VOLA1(0,NY,NX)-VOLW1(0,NY,NX)
     2-VOLI1(0,NY,NX))
      VOLWM(1,0,NY,NX)=VOLW1(0,NY,NX)
      VOLPM(1,0,NY,NX)=VOLP1(0,NY,NX)
      TVOLWI=VOLW1(0,NY,NX)+VOLI1(0,NY,NX)
      XVOLT(NY,NX)=AMAX1(0.0,TVOLWI-VOLWRX(NY,NX))
      IF(TVOLWI.GT.ZEROS(NY,NX))THEN
      VOLWRZ=VOLW1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      VOLIRZ=VOLI1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      XVOLW(NY,NX)=AMAX1(0.0,VOLW1(0,NY,NX)-VOLWRZ)
      XVOLI(NY,NX)=AMAX1(0.0,VOLI1(0,NY,NX)-VOLIRZ)
      ELSE
      XVOLW(NY,NX)=0.0
      XVOLI(NY,NX)=0.0
      ENDIF
      XVOLTM(1,NY,NX)=XVOLT(NY,NX)
      XVOLWM(1,NY,NX)=XVOLW(NY,NX)
      XVOLIM(1,NY,NX)=XVOLI(NY,NX)
      IF(VOLR(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWX(0,NY,NX)=AMAX1(0.0,VOLW1(0,NY,NX)/VOLR(NY,NX))
      THETIX(0,NY,NX)=AMAX1(0.0,VOLI1(0,NY,NX)/VOLR(NY,NX))
      THETPX(0,NY,NX)=AMAX1(0.0,VOLP1(0,NY,NX)/VOLR(NY,NX))
     2*AMAX1(0.0,(1.0-XVOLT(NY,NX)/VOLWD(NY,NX)))
      ELSE
      THETWX(0,NY,NX)=0.0
      THETIX(0,NY,NX)=0.0
      THETPX(0,NY,NX)=1.0
      ENDIF
      THETPM(1,0,NY,NX)=THETPX(0,NY,NX)
      PSISM1(0,NY,NX)=PSISM(0,NY,NX)
      TK1(0,NY,NX)=TKS(0,NY,NX)
C     WRITE(*,7751)'THETPX',I,J,NX,NY
C    2,VOLW1(0,NY,NX),VOLI1(0,NY,NX),VOLP1(0,NY,NX),VOLA1(0,NY,NX)
C    3,THETWX(0,NY,NX),THETIX(0,NY,NX),THETPX(0,NY,NX)
C    3,XVOLW(NY,NX),XVOLT(NY,NX),VOLWD(NY,NX),VOLWG(NY,NX),VOLR(NY,NX)
C    4,VOLW1(1,NY,NX),VOLI1(1,NY,NX),VOLP1(1,NY,NX)
C    5,VOLWRX(NY,NX)
7751  FORMAT(A8,4I4,20E12.4)
C
C     SNOW AND RESIDUE COVERAGE OF SOIL SURFACE
C
C     FSNW,FSNX=fractions of snow,snow-free cover
C     DPTHS=snowpack depth
C     DPTHSX=minimum snowpack depth for full cover
C     BARE,CVRD=fractions of soil,litter cover
C     PRECA=precipitation+irrigation
C     PRECD,PRECB=direct,indirect precipn+irrign at soil surface
C
      FSNW(NY,NX)=AMIN1(1.0,SQRT((DPTHS(NY,NX)/DPTHSX)))
      FSNX(NY,NX)=1.0-FSNW(NY,NX)
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      BARE(NY,NX)=AMIN1(1.0,AMAX1(0.0,EXP(-0.8E-02
     2*(ORGC(0,NY,NX)/AREA(3,0,NY,NX)))))
      ELSE
      BARE(NY,NX)=1.0
      ENDIF
      CVRD(NY,NX)=1.0-BARE(NY,NX)
      PRECM(NY,NX)=1.0E+03*PRECA(NY,NX)
     2/AREA(3,NU(NY,NX),NY,NX)
      PRECD(NY,NX)=1.0E+03*(PRECA(NY,NX)-TFLWCI(NY,NX))
     2/AREA(3,NU(NY,NX),NY,NX)
      PRECB(NY,NX)=1.0E+03*(TFLWCI(NY,NX)-TFLWC(NY,NX))
     2/AREA(3,NU(NY,NX),NY,NX)
C     IF(PRECA(NY,NX).GT.0.0)THEN
C     WRITE(*,3112)'BARE',I,J,NX,NY,BARE(NY,NX)
C    2,FSNX(NY,NX),ORGC(0,NY,NX)/AREA(3,0,NY,NX),VOLWRX(NY,NX)
C    3,XVOLW(NY,NX),VOLWD(NY,NX)
C    4,PRECA(NY,NX),TFLWCI(NY,NX),TFLWC(NY,NX)
C    5,PRECA(NY,NX)*XNPH*1000*BARE(NY,NX),PRECD(NY,NX),PRECB(NY,NX)
3112  FORMAT(A8,4I4,20E12.4)
C     ENDIF
C
C     RESIDUE WATER ABSORPTION CAPACITY
C
C     HCNDR=litter saturated hydraulic conductivity
C     DLYRR=litter depth
C
      HCNDR(NY,NX)=HCNDRR
      DLYRR(NY,NX)=AMAX1(2.5E-03,DLYR(3,0,NY,NX))
C
C     DISTRIBUTION OF PRECIPITATION AND ITS HEAT AMONG SURFACE
C     RESIDUE, SOIL SURFACE, AND MACROPORES
C
C     PRECA,PRECW=rainfall+irrigation,snowfall (water equiv)
C     FLWQW=rainfall to snowpack
C     FLWSW=snowfall to snowpack
C     HFLWSW=convective heat flux to snowpack
C     FLWQB=precip to litter+soil surfaces
C     FLWQAX,FLWQBX=precip to soil,litter surfaces
C     HFLWQA,HFLWQB=convective heat flux to soil,litter surfaces
C     FLWQAS,FLWQAH=precip to soil micropores,macropores
C
      IF(PRECA(NY,NX).GT.0.0.OR.PRECW(NY,NX).GT.0.0)THEN
      FLWQW=(PRECA(NY,NX)-TFLWC(NY,NX))*FSNW(NY,NX)
      FLWSW=PRECW(NY,NX)
      HFLWSW=2.095*TKA(NY,NX)*FLWSW+4.19*TKA(NY,NX)*FLWQW
      FLWQB=(PRECA(NY,NX)-TFLWC(NY,NX))*FSNX(NY,NX)
      FLWQBX=FLWQB*CVRD(NY,NX)
      HFLWQB=4.19*TKA(NY,NX)*FLWQBX
      FLWQAX=FLWQB*BARE(NY,NX)
      HFLWQA=4.19*TKA(NY,NX)*FLWQAX
      FLWQAS=FLWQAX*FGRD(NUM(NY,NX),NY,NX)
      FLWQAH=FLWQAX*FMAC(NUM(NY,NX),NY,NX)
      ELSE
      FLWQW=-TFLWC(NY,NX)*FSNW(NY,NX)
      FLWSW=0.0
      HFLWSW=4.19*TKA(NY,NX)*FLWQW
      FLWQB=-TFLWC(NY,NX)*FSNX(NY,NX)
      FLWQBX=FLWQB*CVRD(NY,NX)
      HFLWQB=4.19*TKA(NY,NX)*FLWQBX
      FLWQAX=FLWQB*BARE(NY,NX)
      HFLWQA=4.19*TKA(NY,NX)*FLWQAX
      FLWQAS=FLWQAX*FGRD(NUM(NY,NX),NY,NX)
      FLWQAH=FLWQAX*FMAC(NUM(NY,NX),NY,NX)
      ENDIF
C
C     PRECIP ON SNOW ARRAYS EXPORTED TO TRNSFR.F, TRNSFRS.F
C     FOR SOLUTE FLUX CALCULATIONS
C
C     PRECW,PRECR,PRECQ,PRECI=snow,rain,snow+rain,irrigation
C     VHCPW,VHCPWX=current, minimum snowpack heat capacities
C     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
C     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
C
      IF(PRECW(NY,NX).GT.0.0.OR.(PRECR(NY,NX).GT.0.0
     2.AND.VHCPW(1,NY,NX).GT.VHCPWX(NY,NX)))THEN
      FLQRQ(NY,NX)=0.0
      FLQRI(NY,NX)=0.0
      FLQGQ(NY,NX)=PRECQ(NY,NX)
      FLQGI(NY,NX)=PRECI(NY,NX)
      ELSEIF((PRECQ(NY,NX).GT.0.0.OR.PRECI(NY,NX).GT.0.0)
     2.AND.VHCPW(1,NY,NX).LE.VHCPWX(NY,NX))THEN
      FLQRQ(NY,NX)=FLWQBX*PRECQ(NY,NX)/(PRECQ(NY,NX)+PRECI(NY,NX))
      FLQRI(NY,NX)=FLWQBX*PRECI(NY,NX)/(PRECQ(NY,NX)+PRECI(NY,NX))
      FLQGQ(NY,NX)=PRECQ(NY,NX)-FLQRQ(NY,NX)
      FLQGI(NY,NX)=PRECI(NY,NX)-FLQRI(NY,NX)
      ELSE
      FLQRQ(NY,NX)=0.0
      FLQRI(NY,NX)=0.0
      FLQGQ(NY,NX)=0.0
      FLQGI(NY,NX)=0.0
      ENDIF
C
C     GATHER PRECIPITATION AND MELTWATER FLUXES AND THEIR HEATS
C     AMONG ATMOSPHERE, SNOWPACK, RESIDUE AND SOIL SURFACES
C     INTO LOCAL ARRAYS FOR USE IN MASS AND ENERGY EXCHANGE
C     ALGORITHMS
C
C     XNPH=internal time step for fluxes through soil profile
C
C     FLW0S,FLQ0I,FLQ0W=snow,ice,water input to snowpack
C     HWFLQ0=convective heat flux to snowpack
C     FLQ1,FLH1,FLY1=rain+irrigation to micropores,macropores,litter
C     HWFLQ1,HWFLY1=convective heat flux to soil,litter surfaces
C
      FLQ0S(NY,NX)=FLWSW*XNPH
      FLQ0I(NY,NX)=0.0
      FLQ0W(NY,NX)=FLWQW*XNPH
      HWFLQ0(NY,NX)=HFLWSW*XNPH
      FLQ1(NY,NX)=FLWQAS*XNPH
      FLH1(NY,NX)=FLWQAH*XNPH
      FLY1(NY,NX)=FLWQBX*XNPH
      HWFLQ1(NY,NX)=HFLWQA*XNPH
      HWFLY1(NY,NX)=HFLWQB*XNPH
C     IF(I.EQ.118.AND.NX.EQ.3.AND.NY.EQ.4)THEN
C     WRITE(*,4422)'FLQ0W',I,J,FLQ0W(NY,NX),FLWQW,XNPH
C     WRITE(*,4422)'FLY',I,J,PRECA(NY,NX),TFLWC(NY,NX),FLY1(NY,NX)
C    2,PSISM1(0,NY,NX),PSISM(0,NY,NX)
C    2,FLQ1(NY,NX),FLH1(NY,NX),FLWQBX
C    2,FLWQAS,FLWQAH
C    3,FGRD(NUM(NY,NX),NY,NX),FMAC(NUM(NY,NX),NY,NX)
C    4,FHOL(L,NY,NX),VOLAH1(L,NY,NX),VOLAH1(L,NY,NX)
C    5,FLWQAX,PRECA(NY,NX),TFLWC(NY,NX),FLWQBX
C    6,BARE(NY,NX),ORGC(0,NY,NX),XVOLW(NY,NX),VOLWG(NY,NX)
C    7,VOLW1(0,NY,NX),VOLWRX(NY,NX)
4422  FORMAT(A8,2I4,60F18.6)
C     ENDIF
C
C     INITIALIZE PARAMETERS, FLUXES FOR ENERGY EXCHANGE
C     AT SNOW, RESIDUE AND SOIL SURFACES
C
C     RADGX=shortwave radiation at ground surface
C     RADXW,RADXG,RADXR= shortwave radn at snowpack,soil,litter
C     FRADG=fraction of shortwave radiation at ground surface
C     FSNW,FSNX=fractions of snow,snow-free cover
C     BARE,CVRD=fractions of soil,litter cover
C     XNPS=internal time step for fluxes through snowpack
C     THRYX=longwave radiation at ground surface
C     THRYW,THRYG,THRYR=longwave radn incident at snowpack,soil,litter
C     THRMW,THRMS,THRMR=longwave radn emitted by snowpack,soil,litter
C     EMMW,EMMS,EMMR=emissivity of snowpack,soil,litter surfaces
C
      RADGX=RADG(NY,NX)*XNPH
      RADXW(NY,NX)=RADGX*FSNW(NY,NX)*XNPS
      RADXG(NY,NX)=RADGX*FSNX(NY,NX)*BARE(NY,NX)
      RADXR(NY,NX)=RADGX*FSNX(NY,NX)*CVRD(NY,NX)*XNPR
      THRYX=(THS(NY,NX)*FRADG(NY,NX)+THRMCX(NY,NX))*XNPH
      THRYW(NY,NX)=THRYX*FSNW(NY,NX)*XNPS
      THRYG(NY,NX)=THRYX*FSNX(NY,NX)*BARE(NY,NX)
      THRYR(NY,NX)=THRYX*FSNX(NY,NX)*CVRD(NY,NX)*XNPR
      THRMW(NY,NX)=EMMW*2.04E-10*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNW(NY,NX)*XNPY
      THRMS(NY,NX)=EMMS*2.04E-10*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNX(NY,NX)*BARE(NY,NX)*XNPH
      THRMR(NY,NX)=EMMR*2.04E-10*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNX(NY,NX)*CVRD(NY,NX)*XNPZ
C
C     AERODYNAMIC RESISTANCE OF CANOPY TO SNOW/RESIDUE/SOIL
C     SURFACE ENERGY EXCHANGE WITH ATMOSPHERE
C
C     ALFZ=parameter for canopy effect on windspeed
C     FRADG=fraction of shortwave radiation at ground surface
C     RAB,RAC=isothermal blr above canopy, canopy blr
C     ZT,ZS=canopy, surface roughness heights
C     UA,UAG=windspeeds above,below canopy
C     VPQ,VPA=vapor pressure within,above canopy
C     TKQ,TKA=temperature within,above canopy
C     TLEX,TSHX=net latent,sensible heat fluxes x blrs from prev hour
C     VAP=latent heat of evaporation
C     1.25E-03=heat capacity of air
C     AREA=surface area of grid cell
C
      ALFZ=2.0*(1.0-FRADG(NY,NX))
      IF(RAB(NY,NX).GT.ZERO.AND.ZT(NY,NX).GT.ZS(NY,NX)
     2.AND.ALFZ.GT.ZERO)THEN
      RAC(NY,NX)=AMIN1(RACX,AMAX1(0.0,ZT(NY,NX)*EXP(ALFZ)
     2/(ALFZ/RAB(NY,NX))*AMAX1(0.0,EXP(-ALFZ*ZS(NY,NX)/ZT(NY,NX))
     3-EXP(-ALFZ*(ZD(NY,NX)+ZR(NY,NX))/ZT(NY,NX)))))
      UAG=UA(NY,NX)*EXP(-ALFZ)
      ELSE
      RAC(NY,NX)=0.0
      UAG=UA(NY,NX)
      ENDIF
      VPQ(NY,NX)=VPA(NY,NX)-TLEX(NY,NX)
     2/(VAP*AREA(3,NUM(NY,NX),NY,NX))
      TKQ(NY,NX)=TKA(NY,NX)-TSHX(NY,NX)
     2/(1.25E-03*AREA(3,NUM(NY,NX),NY,NX))
C     WRITE(*,3114)'RAC',I,J,NX,NY,RAC(NY,NX),FRADG(NY,NX)
C    2,RAB(NY,NX),ZT(NY,NX),ZS(NY,NX),ALFZ
C    3,VPQ(NY,NX),TKQ(NY,NX)
3114  FORMAT(A8,4I4,20E12.4)
C
C     AERODYNAMIC RESISTANCE OF SNOWPACK, RESIDUE AND SOIL
C     SURFACES TO ENERGY EXCHANGE WITH ATMOSPHERE
C     Soil Sci. Soc. Am. J. 48:25-32
C
C     RAR=porosity-unlimited litter blr
C     DLYRR=litter depth
C     WGSGR=vapor diffusivity in litter
C     RAG,RAGW,RAGR=isothermal blrs at ground,snowpack,litter surfaces
C     RARX=blr of litter surface
C     THETPX*=air-filled porosity of litter
C     DFVR=porosity limitation to diffusion through litter
C     POROQ=litter tortuosity
C     RAR1=porosity-limited litter blr
C     PAREX,PARSX=conductances for latent,sensible heat fluxes
C     PAREW,PARSW=conductances for snowpack latent,sensible heatfluxes
C     PAREG,PARSG=conductances for soil latent,sensible heat fluxes
C     PARER,PARSR=conductances for litter latent,sensible heat fluxes
C     XNPR=internal time step for fluxes through litter
C
      RAR(NY,NX)=DLYRR(NY,NX)/WGSGR(NY,NX)
      RAG(NY,NX)=RAC(NY,NX)+RAB(NY,NX)
      RAGW(NY,NX)=RAG(NY,NX)
      RAGR(NY,NX)=RAG(NY,NX)+RARX
      RARG(NY,NX)=RAGR(NY,NX)
      THETPX0=AMAX1(ZERO2,THETPX(0,NY,NX))
      DFVR=THETPX0*POROQ*THETPX0/POROS(0,NY,NX)
      RAR1=RAG(NY,NX)+RAR(NY,NX)/DFVR
      PAREX=AREA(3,NUM(NY,NX),NY,NX)*XNPH
      PARSX=1.25E-03*AREA(3,NUM(NY,NX),NY,NX)*XNPH
      PAREW(NY,NX)=PAREX*FSNW(NY,NX)*XNPS
      PARSW(NY,NX)=PARSX*FSNW(NY,NX)*XNPS
      PAREG(NY,NX)=PAREX*FSNX(NY,NX)
      PARER(NY,NX)=PAREX*FSNX(NY,NX)*XNPR*CVRD(NY,NX)
      PARSG(NY,NX)=PARSX*FSNX(NY,NX)
      PARSR(NY,NX)=PARSX*FSNX(NY,NX)*XNPR*CVRD(NY,NX)
C     WRITE(*,3115)'RAR',I,J,NX,NY,RAR1,RAR(NY,NX),DFVR
C    2,THETPX(0,NY,NX),THETPX0,POROS(0,NY,NX),DLYRR(NY,NX)
C    3,RAG(NY,NX),PARER(NY,NX),PAREX,FSNX(NY,NX),XNPR,CVRD(NY,NX)
C    4,BARE(NY,NX)
3115  FORMAT(A8,4I4,30E12.4)
C
C     BOUNDARY LAYER CONDUCTANCES FOR EXPORT TO TRNSFR.F
C
C     RAS,RASL=blrs of snowpack,snowpack layer
C     VOLS,VOLS1=volume of snowpack,snowpack layer
C     DLYRS=snowpack later depth
C     WGSGW=vapor diffusivity in snowpack
C     THETPL=snowpack air-filled porosity
C     VOLS0,VOLI0,VOLW0,VOLS1=snow,ice,water,total volumes of snowpack
C     PARR=boundary layer conductance above litter,soil surfaces
C
      RAS=0.0
      IF(VOLS(NY,NX).GT.ZEROS2(NY,NX))THEN
      DO 9775 L=1,JS
      IF(VOLS1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      RASX=DLYRS(L,NY,NX)/WGSGW(L,NY,NX)
      THETPL=AMAX1(THETPI,1.0-(VOLS0(L,NY,NX)+VOLI0(L,NY,NX)
     2+VOLW0(L,NY,NX))/VOLS1(L,NY,NX))
      RASL=RASX/AMAX1(ZERO,THETPL)**2.0
      RAS=RAS+RASL
C     WRITE(*,3113)'RAS',I,J,NX,NY,L,ALFZ,RAS,RASL,RASX
C    2,DLYRS(L,NY,NX),WGSGW(L,NY,NX),THETPL,THETPI,VOLS0(L,NY,NX)
C    3,VOLI0(L,NY,NX),VOLW0(L,NY,NX),VOLS1(L,NY,NX),TKW(L,NY,NX)
3113  FORMAT(A8,5I4,40E12.4)
      ENDIF
9775  CONTINUE
      ENDIF
      PARR(NY,NX)=AREA(3,NUM(NY,NX),NY,NX)*XNPH/(RAGR(NY,NX)+RAS)
C     IF(NX.EQ.1)THEN
C     WRITE(*,3111)'RAC',I,J,NX,NY,ALFZ,RAC(NY,NX)
C    2,PARR(NY,NX),ZT(NY,NX),RAB(NY,NX)
C    2,RAGR(NY,NX),RAS,VOLS(NY,NX)
C    3,DLYRR(NY,NX),RAG(NY,NX),RAGR(NY,NX)
C    4,THETX,THETPX(0,NY,NX),VHCP1(0,NY,NX)
C    5,WGSGR(NY,NX),VOLW1(0,NY,NX)
C    5,VOLI1(0,NY,NX),VOLP1(0,NY,NX),VOLR(NY,NX),VOLA1(0,NY,NX)
C    4,TLEX(NY,NX),TSHX(NY,NX),RADG(NY,NX),THS(NY,NX)
C    5,FRADG(NY,NX),THRMCX(NY,NX),ZS(NY,NX)
C    6,XVOLW(NY,NX),VHCPRX(NY,NX)/4.19,VOLWD(NY,NX),ORGC(0,NY,NX)
3111  FORMAT(A8,4I4,40E12.4)
C     ENDIF
9990  CONTINUE
9995  CONTINUE
C
C     INITIALIZE SOIL HYDRAULIC PARAMETERS IN LOCAL ARRAYS
C     FOR LATER USE IN WATER TRANSFER ALGORITHMS
C
C     N3,N2,N1=L,NY,NX of source grid cell
C     N6,N5,N4=L,NY,NX of destination grid cell
C
      DO 9985 NX=NHW,NHE
      DO 9980 NY=NVN,NVS
      DO 35 L=NUM(NY,NX),NL(NY,NX)
      DO 40 N=NCN(NY,NX),3
      N1=NX
      N2=NY
      N3=L
      IF(N.EQ.1)THEN
      IF(NX.EQ.NHE)THEN
      GO TO 50
      ELSE
      N4=NX+1
      N5=NY
      N6=L
      ENDIF
      ELSEIF(N.EQ.2)THEN
      IF(NY.EQ.NVS)THEN
      GO TO 50
      ELSE
      N4=NX
      N5=NY+1
      N6=L
      ENDIF
      ELSEIF(N.EQ.3)THEN
      IF(L.EQ.NL(NY,NX))THEN
      GO TO 50
      ELSE
      N4=NX
      N5=NY
      N6=L+1
      ENDIF
      ENDIF
C
C     MACROPORE CONDUCTIVITY FROM 'HOUR1' AND GRAVITATIONAL
C     GRADIENT USED TO CALCULATE MACROPORE FLOW FOR USE BELOW
C
C     CNDH1=macropore hydraulic conductivity
C     AVCNHL=macropore hydraulic conductance
C     DLYR=layer depth
C
      IF(CNDH1(N3,N2,N1).GT.ZERO.AND.CNDH1(N6,N5,N4)
     2.GT.ZERO)THEN
      AVCNHL(N,N6,N5,N4)=2.0*CNDH1(N3,N2,N1)*CNDH1(N6,N5,N4)
     2/(CNDH1(N3,N2,N1)*DLYR(N,N6,N5,N4)+CNDH1(N6,N5,N4)
     3*DLYR(N,N3,N2,N1))
      ELSE
      AVCNHL(N,N6,N5,N4)=0.0
      ENDIF
50    CONTINUE
40    CONTINUE
35    CONTINUE
9980  CONTINUE
9985  CONTINUE
C
C     DYNAMIC LOOP FOR FLUX CALCULATIONS
C
      DO 3320 M=1,NPH
      DO 9895 NX=NHW,NHE
      DO 9890 NY=NVN,NVS
C
C     INITIALIZE NET SURFACE FLUX ACCUMULATORS
C
C     TQR1,TQS1,TQW1,TQI1=net water and snowpack snow,water,ice runoff
C     THQR1,THQS1=net convective heat from surface water and snow runoff
C     BAREW,CVRDW=fractions of soil,litter cover including free water+ice
C     RAGS= boundary layer resistance at soil surface
C     PARG=boundary layer conductance above soil surface
C
      TQR1(NY,NX)=0.0
      THQR1(NY,NX)=0.0
      TQS1(NY,NX)=0.0
      TQW1(NY,NX)=0.0
      TQI1(NY,NX)=0.0
      THQS1(NY,NX)=0.0
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      BAREW(NY,NX)=AMAX1(0.0,BARE(NY,NX)
     2-AMIN1(1.0,AMAX1(0.0,XVOLT(NY,NX)/VOLWD(NY,NX))))
      ELSE
      BAREW(NY,NX)=1.0
      ENDIF
      CVRDW(NY,NX)=1.0-BAREW(NY,NX)
      RAGS(NY,NX)=1.0/(BAREW(NY,NX)/RAGR(NY,NX)+CVRDW(NY,NX)/RAR1)
      PARG(M,NY,NX)=AREA(3,NUM(NY,NX),NY,NX)*XNPH/(RAGS(NY,NX)+RAS)
C
C     REDISTRIBUTE INCOMING PRECIPITATION
C     BETWEEN RESIDUE AND SOIL SURFACE
C
C     BKDS=bulk density
C     FLQRS,FLQRH=water flux from soil micropores,macropores to litter
C     FLQ1,FLH1,FLY1=rain+irrigation to micropores,macropores,litter
C     VOLP1,VOLPH1=air-filled microporosity,macroporosity
C     HFLQR1=convective heat flux from soil to litter
C     FLYM,HWFLYM=total water flux, convective heat flux to litter
C     FLQM,FLHM=total water flux to soil micropores, macropores
C     HWFLQM=total convective heat flux to soil micropores, macropores
C     XNPR=time step for litter water,heat flux calculations
C
      IF(BKDS(NUM(NY,NX),NY,NX).GT.ZERO)THEN
      FLQRS=AMAX1(0.0,FLQ1(NY,NX)-VOLP1(NUM(NY,NX),NY,NX))
      FLQRH=AMAX1(0.0,FLH1(NY,NX)-VOLPH1(NUM(NY,NX),NY,NX))
      HFLQR1=4.19*TKA(NY,NX)*(FLQRS+FLQRH)
      FLYM=FLY1(NY,NX)+FLQRS+FLQRH
      HWFLYM=HWFLY1(NY,NX)+HFLQR1
      FLQM=FLQ1(NY,NX)-FLQRS
      FLHM=FLH1(NY,NX)-FLQRH
      HWFLQM=HWFLQ1(NY,NX)-HFLQR1
      ELSE
      FLYM=FLY1(NY,NX)
      HWFLYM=HWFLY1(NY,NX)
      FLQM=FLQ1(NY,NX)
      FLHM=FLH1(NY,NX)
      HWFLQM=HWFLQ1(NY,NX)
      ENDIF
      FLYM2=FLYM*XNPR
      HWFLM2=HWFLYM*XNPR
C
C     WATER GAS EXCHANGE COEFFICIENTS IN SURFACE LITTER
C
C     VOLA1,VOLI1,VOLW1,VOLPM=total,ice-,water-,air-filled porosity
C     TFND1=temperature effect on gas diffusivity
C     DFGS=rate constant for air-water gas exchange
C     Z1R,Z2RW,Z2RD,Z3RX=parameters for litter air-water gas transfers
C     XNPD=time step for gas transfer calculations
C     TORT=tortuosity for aqueous diffusivity
C
      VOLAT0=VOLA1(0,NY,NX)-VOLI1(0,NY,NX)
      IF(VOLAT0.GT.ZEROS2(NY,NX)
     2.AND.VOLPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWA=AMAX1(0.0,AMIN1(1.0
     2,VOLW1(0,NY,NX)/VOLAT0))
      TFND1=(TK1(0,NY,NX)/298.15)**6
      IF(THETWA.GT.Z3R)THEN
      DFGS(M,0,NY,NX)=AMAX1(0.0
     2,TFND1*XNPD/((Z1R**-1)*EXP(Z2RW*(THETWA-Z3R))))
      ELSE
      DFGS(M,0,NY,NX)=AMIN1(1.0
     2,TFND1*XNPD/((Z1R**-1)*EXP(Z2RD*(THETWA-Z3R))))
      ENDIF
      ELSE
      DFGS(M,0,NY,NX)=0.0
      ENDIF
      IF(VOLWRX(NY,NX).GT.ZEROS(NY,NX))THEN
      THETWT=AMIN1(1.0,VOLW(0,NY,NX)/VOLWRX(NY,NX))
      ELSE
      THETWT=1.0
      ENDIF
      TORT(M,0,NY,NX)=0.7*THETWT**2
C
C     KINETIC ENERGY OF DIRECT RAINFALL AND THROUGHFALL
C
C     PRECD,PRECB=direct,indirect precipn+irrign at soil surface
C     ENGYD,ENGYB=energy impact of direct,indirect precipn+irrign at soil surface
C     VOLWG=ground surface water retention capacity
C     XVOLW=free surface water
C     ZT=canopy height
C     ENGYPM=total energy impact for use in erosion.f
C     ENGYP=cumulative rainfall energy impact on soil surface
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
C
      IF(PRECD(NY,NX).GT.ZERO)THEN
      ENGYD=AMAX1(0.0,8.95+8.44*LOG(PRECM(NY,NX)))
      ELSE
      ENGYD=0.0
      ENDIF
      IF(PRECB(NY,NX).GT.ZERO)THEN
      ENGYB=AMAX1(0.0,15.8*SQRT(AMIN1(2.5,ZT(NY,NX)))-5.87)
      ELSE
      ENGYB=0.0
      ENDIF
      IF(ENGYD+ENGYB.GT.ZERO)THEN
      HV=1.0E+03*AMAX1(0.0,XVOLT(NY,NX)-VOLWG(NY,NX))
     2/AREA(3,NU(NY,NX),NY,NX)
      ENGYPM(M,NY,NX)=(ENGYD*PRECD(NY,NX)+ENGYB*PRECB(NY,NX))
     2*EXP(-2.0*HV)*BARE(NY,NX)*XNPH
      ENGYP(NY,NX)=ENGYP(NY,NX)+ENGYPM(M,NY,NX)
      ELSE
      ENGYPM(M,NY,NX)=0.0
      ENDIF
      FKSAT=EXP(-2.0E-03*(CSILT(NU(NY,NX),NY,NX)
     2+CCLAY(NU(NY,NX),NY,NX))*ENGYP(NY,NX))
C     IF(ENGYD+ENGYB.GT.ZERO)THEN
C     WRITE(*,1117)'FKSAT',I,J,M,NX,NY,FKSAT,ENGYP(NY,NX)
C    2,ENGYPM(M,NY,NX),ENGYD,PRECD(NY,NX),ENGYB,PRECB(NY,NX)
C    3,PRECM(NY,NX),HV,XVOLWM(M,NY,NX),XVOLIM(M,NY,NX)
C    4,XVOLT(NY,NX),VOLWG(NY,NX),ORGC(0,NY,NX),BARE(NY,NX)
C    5,CCLAY(NU(NY,NX),NY,NX),CSILT(NU(NY,NX),NY,NX),ZT(NY,NX)
1117  FORMAT(A8,5I4,20E12.4)
C     ENDIF
C
C     SNOWPACK FLUX ACCUMULATORS
C
C     TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
C     THFLWW=convective heat fluxes of snow,water,ice in snowpack
C
      DO 9875 L=1,JS
      TFLWS(L,NY,NX)=0.0
      TFLWW(L,NY,NX)=0.0
      TFLWI(L,NY,NX)=0.0
      THFLWW(L,NY,NX)=0.0
9875  CONTINUE
C
C     SURFACE FLUX ACCUMULATORS
C
C     TWFLXL,TWFLXH=total freeze-thaw in micropores,macropores
C     TTFLXL=total latent heat from freeze-thaw
C     TFLWL,TFLWHL=net water flux in micropores,macropores
C     THFLWL=net heat flux
C
      DO 9885 L=NUM(NY,NX),NL(NY,NX)
      TWFLXL(L,NY,NX)=0.0
      TWFLXH(L,NY,NX)=0.0
      TTFLXL(L,NY,NX)=0.0
      TFLWL(L,NY,NX)=0.0
      TFLWLX(L,NY,NX)=0.0
      TFLWHL(L,NY,NX)=0.0
      THFLWL(L,NY,NX)=0.0
      VOLW2(L,NY,NX)=VOLW1(L,NY,NX)
C
C     GAS EXCHANGE COEFFICIENTS SOIL LAYERS
C
C     VOLA1,VOLI1,VOLW1=total,ice-,water-filled microporosity
C     VOLAH1,VOLIH1,VOLWH1=total,ice-,water-filled macroporosity
C     VOLPM=air-filled porosity
C     TFND1=temperature effect on gas diffusivity
C     DFGS=rate constant for air-water gas exchange
C     Z1S,Z2SW,Z2SD,Z3SX=parameters for soil air-water gas transfers
C     XNPD=time step for gas transfer calculations
C     TORT,TORTH=tortuosity for aqueous diffn in micropores,macropres
C
      VOLWT=VOLW1(L,NY,NX)+VOLWH1(L,NY,NX)
      VOLAT=VOLA1(L,NY,NX)+VOLAH1(L,NY,NX)
     2-VOLI1(L,NY,NX)-VOLIH1(L,NY,NX)
      IF(VOLAT.GT.ZEROS2(NY,NX)
     2.AND.VOLPM(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWA=AMAX1(0.0,AMIN1(1.0,VOLWT/VOLAT))
      TFND1=(TK1(L,NY,NX)/298.15)**6
      Z3S=AMAX1(Z3SX,FC(L,NY,NX)/POROS(L,NY,NX))
      IF(THETWA.GT.Z3S)THEN
      DFGS(M,L,NY,NX)=AMAX1(0.0
     2,TFND1*XNPD/((Z1S**-1)*EXP(Z2SW*(THETWA-Z3S))))
      ELSE
      DFGS(M,L,NY,NX)=AMIN1(1.0
     2,TFND1*XNPD/((Z1S**-1)*EXP(Z2SD*(THETWA-Z3S))))
      ENDIF
      ELSE
      DFGS(M,L,NY,NX)=0.0
      ENDIF
C     IF(I.EQ.121.AND.L.EQ.2)THEN
C     WRITE(*,3371)'DFGS',I,J,M,NX,NY,L,DFGS(M,L,NY,NX)
C    2,THETWA,VOLWT,VOLAT,VOLW1(L,NY,NX),VOLA1(L,NY,NX)
C    3,VOLWH1(L,NY,NX),VOLAH1(L,NY,NX),BKDS(L,NY,NX),FHOL(L,NY,NX)
C    3,Z3S,Z2S*(THETWA-Z3S),EXP(Z2S*(THETWA-Z3S)),Z1S**-1
C    4,(Z1S**-1)*EXP(Z2S*(THETWA-Z3S))
3371  FORMAT(A8,6I4,20E14.6)
C     ENDIF
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      THETWT=VOLWM(M,L,NY,NX)/VOLY(L,NY,NX)
      TORT(M,L,NY,NX)=0.7*THETWT**2*(1.0-FHOL(L,NY,NX))
      ELSE
      TORT(M,L,NY,NX)=0.7
      ENDIF
      IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWH=VOLWHM(M,L,NY,NX)/VOLAH1(L,NY,NX)
      TORTH(M,L,NY,NX)=AMIN1(1.0,2.8*THETWH**3)*FHOL(L,NY,NX)
      ELSE
      TORTH(M,L,NY,NX)=0.0
      ENDIF
9885  CONTINUE
C     IF(NX.EQ.4.AND.NY.EQ.5)THEN
C     WRITE(*,3132)'FLQR1',I,J,M,NX,NY,FLY1(NY,NX),FLQ1(NY,NX)
C    2,VHCPWM(M,1,NY,NX),VHCPWX(NY,NX)
C    2,FLH1(NY,NX),FLYM,FLQM,FLHM,FLQR1
C    3,FMAC(NUM(NY,NX),NY,NX),FGRD(NUM(NY,NX),NY,NX)
C    5,VOLAH1(NUM(NY,NX),NY,NX),FVOLAH,CCLAY(NUM(NY,NX),NY,NX)
C    4,VOLW1(NUM(NY,NX),NY,NX),VOLX(NUM(NY,NX),NY,NX),WP(L,NY,NX)
C    2,VOLT(NUM(NY,NX),NY,NX),VOLAH1(NUM(NY,NX),NY,NX)
C    5,VOLWRX(NY,NX),VOLW1(0,NY,NX),VOLI1(0,NY,NX)
C    6,PSISM1(0,NY,NX),PSISM1(NUM(NY,NX),NY,NX)
3132  FORMAT(A8,5I4,40E12.4)
C     ENDIF
C
C     ENERGY EXCHANGE VARIABLES AT SNOW SURFACE IF PRESENT
C
C     RFLXW,EFLXW,VFLXW,SFLXW,HFLXW=netradn,latent,convective,sensible
C     and storage heat fluxes
C     FLWLW,FLWHLW=water from snowpack to soil micropores,macropores
C     HFLWLW=conv heat from snowpack to soil micropores,macropores
C     FLWRLW=water flux from snowpack to litter
C     HFLWRLW=convective heat flux from snowpack to litter
C     FLWVLW=snowpack-litter water flux accounting for wetting front
C
      RFLXW=0.0
      EFLXW=0.0
      VFLXW=0.0
      SFLXW=0.0
      HFLXW=0.0
      FLWLW=0.0
      FLWLXW=0.0
      FLWHLW=0.0
      HFLWLW=0.0
      FLWRLW=0.0
      HFLWRLW=0.0
      FLWVLW=0.0
C
C     EVAPS,EVAPW=evaporation from soil,snowpack surfaces
C     FLQRM,FLQSM,FLQHM=water into litter,soil micropores,micropores for use in trnsfr.f
C
      EVAPS(NY,NX)=0.0
      EVAPW(NY,NX)=0.0
      FLQRM(M,NY,NX)=0.0
      FLQSM(M,NY,NX)=0.0
      FLQHM(M,NY,NX)=0.0
C
C     FLUX VARIABLES IN SNOWPACK
C
C     TFLX0=latent heat from freeze-thaw
C     WFLXS,WFLXI=freeze-thaw between snow,ice and water
C     FLW0S,FLW0W,FLW0I=snow,water,ice fluxes
C     HFLW0W=convective heat flux from snow,water,ice fluxes
C     FLQWM=snowpack water flux
C     VOLS0M,VOLW0M,VOLI0M=snow,water,ice contents
C     VHCPWMM=volumetric heat capacity
C     TK0M=temperature
C
      DO 9765 L=1,JS
      TFLX0(L,NY,NX)=0.0
      WFLXS(L,NY,NX)=0.0
      WFLXI(L,NY,NX)=0.0
      FLW0S(L,NY,NX)=0.0
      FLW0W(L,NY,NX)=0.0
      FLW0I(L,NY,NX)=0.0
      HFLW0W(L,NY,NX)=0.0
      FLQWM(M,L,NY,NX)=0.0
      VOLS0M(L,NY,NX)=VOLS0(L,NY,NX)
      VOLW0M(L,NY,NX)=VOLW0(L,NY,NX)
      VOLI0M(L,NY,NX)=VOLI0(L,NY,NX)
      VHCPWMM(L,NY,NX)=VHCPWM(M,L,NY,NX)
      TK0M(L,NY,NX)=TK0(L,NY,NX)
9765  CONTINUE
C
C     HEAT AND VAPOR FLUXES BETWEEN SNOWPACK AND ATMOSPHERE
C
C     VHCPWM=volumetric heat capacity of snowpack
C     NPS=number of cycles for solving snowpack heat and water fluxes
C     ALBW=snowpack albedo
C     VOLS0M,VOLI0M,VOLW0M=snow,ice,water volumes
C     RFLX0=net radiation input
C     RADXW=shortwave radiation at snowpack surface
C     THRYW=longwave radn incident at snowpack surface
C     THRMX=longwave radn emitted by snowpack surface
C     TK0M=snowpack surface temperature
C     RFLXW2=net radiation
C
      IF(VHCPWM(M,1,NY,NX).GT.VHCPWX(NY,NX))THEN
      DO 3000 MM=1,NPS
      ALBW=(0.85*VOLS0M(1,NY,NX)+0.30*VOLI0M(1,NY,NX)
     2+0.06*VOLW0M(1,NY,NX))
     2/(VOLS0M(1,NY,NX)+VOLI0M(1,NY,NX)+VOLW0M(1,NY,NX))
      RFLX0=(1.0-ALBW)*RADXW(NY,NX)+THRYW(NY,NX)
      THRMX=THRMW(NY,NX)*TK0M(1,NY,NX)**4
      RFLXW2=RFLX0-THRMX
C
C     AERODYNAMIC RESISTANCE ABOVE SNOWPACK INCLUDING
C     RESISTANCE IMPOSED BY PLANT CANOPY
C
C     RI=Richardsons number
C     RIB=isothermal RI
C     TKQ=canopy air temperature
C     RAGX,RA=snowpack blr
C     RAG,RAGW=isothermal blrs at ground,snowpack surfaces
C
      RI=AMAX1(-0.3,AMIN1(0.075
     2,RIB(NY,NX)*(TKQ(NY,NX)-TK0M(1,NY,NX))))
      RAGX=AMAX1(RAM,0.8*RAGW(NY,NX),AMIN1(1.2*RAGW(NY,NX)
     2,RAG(NY,NX)/(1.0-10.0*RI)))
      RAGW(NY,NX)=RAGX
      RA=RAGX
C
C     PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
C
C     PARE,PARS=blcs for snowpack latent,sensible heat fluxes
C     PAREW,PARSW=conductances for latent,sensible heat fluxes
C     RZ=surface resistance
C     VP0,VPQ=vapor pressure at snowpack surface, canopy air
C     EVAPT2,EVAPW2,EVAPS2=evaporation total, water,snow
C     XNPS=1/NPS
C     EFLXW2=latent heat flux
C     VAP,VAPS=latent heat of evaporation,sublimation
C     VFLXW2=convective heat of evaporation flux
C
      PARE=PAREW(NY,NX)/(RA+RZ)
      PARS=PARSW(NY,NX)/RA
      VP0=2.173E-03/TK0M(1,NY,NX)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK0M(1,NY,NX)))
      EVAPT2=PARE*(VPQ(NY,NX)-VP0)
      EVAPW2=AMAX1(EVAPT2,-AMAX1(0.0,VOLW0M(1,NY,NX)*XNPA))
      EVAPX2=AMIN1(0.0,EVAPT2-EVAPW2)
      EVAPS2=AMAX1(EVAPX2,-AMAX1(0.0,VOLS0M(1,NY,NX)*XNPA))
      EFLXW2=EVAPW2*VAP+EVAPS2*VAPS
      IF(EVAPT2.LT.0.0)THEN
      VFLXW2=(EVAPW2*4.19+EVAPS2*2.095)*TK0M(1,NY,NX)
      ELSE
      VFLXW2=(EVAPW2*4.19+EVAPS2*2.095)*TKQ(NY,NX)
      ENDIF
C
C     SOLVE FOR SNOWPACK SURFACE TEMPERATURE AT WHICH ENERGY
C     BALANCE OCCURS, SOLVE AND ACCUMULATE LATENT, SENSIBLE
C     STORAGE HEAT FLUXES AND EVAPORATION
C
C     SFLXW2,EFLXW2,RFLXW2=sensible,latent heat fluxes, net radiation
C     VFLXW2=convective heat flux from EFLXW2
C     HFLX02=storage heat flux
C     FLQ0S2,FLQ0W2,FLQ0I2=snow,water,ice input to snowpack
C     HWFLQ02=convective heat from snow,water,ice input to snowpack
C
      SFLXW2=PARS*(TKQ(NY,NX)-TK0M(1,NY,NX))
      HFLX02=RFLXW2+EFLXW2+SFLXW2
      HFLXW2=HFLX02+VFLXW2
      RFLXW=RFLXW+RFLXW2
      EFLXW=EFLXW+EFLXW2
      VFLXW=VFLXW+VFLXW2
      SFLXW=SFLXW+SFLXW2
      HFLXW=HFLXW+HFLXW2
      EVAPS(NY,NX)=EVAPS(NY,NX)+EVAPS2
      EVAPW(NY,NX)=EVAPW(NY,NX)+EVAPW2
      FLQ0S2=FLQ0S(NY,NX)*XNPS
      FLQ0W2=FLQ0W(NY,NX)*XNPS
      FLQ0I2=FLQ0I(NY,NX)*XNPS
      HWFLQ02=HWFLQ0(NY,NX)*XNPS
      FLW0S2=FLQ0S2+EVAPS2
      FLW0W2=FLQ0W2+EVAPW2
      FLW0I2=FLQ0I2
      HFLW0W2=HWFLQ02+HFLXW2
      FLW0S(1,NY,NX)=FLW0S2
      FLW0W(1,NY,NX)=FLW0W2
      FLW0I(1,NY,NX)=FLW0I2
      HFLW0W(1,NY,NX)=HFLW0W2
      FLQWM(M,1,NY,NX)=FLQWM(M,1,NY,NX)+FLQ0S2+FLQ0I2+FLQ0W2
C     IF(NX.EQ.3.AND.NY.EQ.3)THEN
C     WRITE(*,7759)'EVAP',I,J,M,MM,FLW0S2
C    2,FLQ0S2,EVAPS2,FLW0W2,FLQ0W2
C    3,FSNW(NY,NX),FLW0I2,FLQ0I2,RFLXW2,EFLXW2
C    4,SFLXW2,VFLXW2,RA,EVAPT2,EVAPX2,VPQ(NY,NX),VP0
C    5,VOLW0M(1,NY,NX),VOLS0M(1,NY,NX),VOLI0M(1,NY,NX)
C    6,HFLW0W(1,NY,NX),HFLW0W2,HWFLQ02,HFLXW2,RFLXW2,EFLXW2
C    7,SFLXW2,VFLXW2,TK0M(1,NY,NX),TKQ(NY,NX)
C    8,PARE,RA,RZ,EVAPS2,EVAPW2,EVAPT2
7759  FORMAT(A8,4I4,40E14.6)
C     ENDIF
C
C     PHYSICAL AND HYDRAULIC PROPERTIES OF SNOWPACK INCLUDING
C     AIR AND WATER-FILLED POROSITY, WATER POTENTIAL OF UNDERLYING
C     SOIL SURFACE USED IN FLUX CALCULATIONS
C
C     VHCPW,VHCPWX=current, minimum snowpack heat capacities
C     VOLS0M,VOLI0M,VOLW0M,VOLS1=snow,ice,water,total snowpack volume
C     DENSS,DENSI,DENS0=snow,ice,minimum snow density
C     AREA=area of grid cell
C     VOLP01=snowpack air volume
C     THETP1=snowpack air concentration
C     CNV1=snowpack vapor conductivity
C     VP1=snowpack vapor concentration
C     TK0M=snowpack temperature
C     WGSGW=vapor diffusivity
C     DENSW1=snowpack density
C
      ICHKL=0
      DO 9880 L=1,JS
      IF(VHCPWMM(L,NY,NX).GT.VHCPWX(NY,NX))THEN
      VOLS1(L,NY,NX)=VOLS0M(L,NY,NX)/DENSS(L,NY,NX)
     2+VOLW0M(L,NY,NX)+VOLI0M(L,NY,NX)
      DLYRS0(L,NY,NX)=VOLS1(L,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
      VOLP01=AMAX1(0.0,VOLS1(L,NY,NX)-VOLS0M(L,NY,NX)-VOLI0M(L,NY,NX)
     2-VOLW0M(L,NY,NX))
      THETP1=AMAX1(THETPI,VOLP01/VOLS1(L,NY,NX))
      CNV1=THETP1**2.0*WGSGW(L,NY,NX)
      VP1=2.173E-03/TK0M(L,NY,NX)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK0M(L,NY,NX)))
      IF(VOLS1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DENSW1=AMIN1(0.6,(VOLS0M(L,NY,NX)+VOLW0M(L,NY,NX)
     2+VOLI0M(L,NY,NX)*DENSI)/VOLS1(L,NY,NX))
      ELSE
      DENSW1=DENS0(NY,NX)
      ENDIF
C
C     SNOW THERMAL CONDUCTIVITY FROM J GLACIOL 43:26-41
C
C     TCND1W=snow thermal conductivity
C     FLWQX=porosity-unconstrained snow water flux
C
      TCND1W=0.0036*10**(2.650*DENSW1-1.652)
C
C     DISCHARGE OF MELTWATER AND ITS HEAT FROM SNOWPACK LAYER
C     TO LOWER SNOWPACK LAYER
C
      FLWQX=AMAX1(0.0,AMAX1(0.0,VOLW0M(L,NY,NX))
     2-0.05*AMAX1(0.0,VOLS0M(L,NY,NX)))*XNPA
C
C     WATER AND HEAT FLUXES IN SNOWPACK
C
C     DLYRS0=snow layer thickness
C     FLWQM=porosity-constrained snow water flux
C     HFLWQM=convective heat flux from water flux
C
      L2=MIN(JS,L+1)
      IF(L.LT.JS.AND.VHCPWMM(L2,NY,NX).GT.VHCPWX(NY,NX))THEN
      VOLS1(L2,NY,NX)=VOLS0M(L2,NY,NX)/DENSS(L2,NY,NX)
     2+VOLW0M(L2,NY,NX)+VOLI0M(L2,NY,NX)
      DLYRS0(L2,NY,NX)=VOLS1(L2,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
      VOLP02=VOLS1(L2,NY,NX)-VOLS0M(L2,NY,NX)-VOLI0M(L2,NY,NX)
     2-VOLW0M(L2,NY,NX)
      THETP2=AMAX1(THETPI,VOLP02/VOLS1(L2,NY,NX))
      FLWQM=AMIN1(THETP2,FLWQX)
      HFLWQM=4.19*TK0M(L,NY,NX)*FLWQM
C
C     VAPOR FLUX IN SNOWPACK
C
C     VOLP01,VOLP02=air-filled volumes of source, destination layers
C     L2=destination layer
C     CNV1,CNV2=vapor conductivities of source, destination layers
C     VP1,VP2=vapor concentrations of source, destination layers
C     TK0M=soil temperature
C     ATCNVW=snow vapor conductance
C     DLYRS0=snow layer thickness
C     FLVC,FLVX=vapor-unconstrained,vapor-constrained vapor flux
C     FLVSS,HFLVSS=vapor flux and its convective heat flux
C
      IF(VOLP01.GT.ZEROS2(NY,NX).AND.VOLP02.GT.ZEROS2(NY,NX))THEN
      CNV2=THETP2**2.0*WGSGW(L2,NY,NX)
      VP2=2.173E-03/TK0M(L2,NY,NX)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK0M(L2,NY,NX)))
      ATCNVW=2.0*CNV1*CNV2/(CNV1*DLYRS0(L2,NY,NX)
     2+CNV2*DLYRS0(L,NY,NX))
      FLVC=ATCNVW*(VP1-VP2)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)
     2*XNPY
      VPY=(VP1*VOLP01+VP2*VOLP02)/(VOLP01+VOLP02)
      FLVX=(VP1-VPY)*VOLP01*XNPA
      IF(FLVC.GE.0.0)THEN
      FLVSS=AMAX1(0.0,AMIN1(FLVC,FLVX,VOLW0M(L,NY,NX)*XNPX))
      HFLVSS=(4.19*TK0M(L,NY,NX)+VAP)*FLVSS
      ELSE
      FLVSS=AMIN1(0.0,AMAX1(FLVC,FLVX,-VOLW0M(L2,NY,NX)*XNPX))
      HFLVSS=(4.19*TK0M(L2,NY,NX)+VAP)*FLVSS
      ENDIF
      ELSE
      FLVSS=0.0
      HFLVSS=0.0
      ENDIF
C
C     HEAT FLUX IN SNOWPACK
C
C     DENSW2,TCNDW2=density,thermal conductivity in destination layer
C     ATCNDW=thermal conductance
C     DLYRS0=layer thickness
C     TKY=equilibrium temperature
C     HFLWX,HFLWC=heat-constrained,heat-unconstrained heat fluxes
C     VHCPWMM,TK0M=volumetric heat capacity,temperature
C     XNPX=time step for flux calculations
C     FSNW=snow cover fraction
C     XNPY=time step for snowpack flux calculations
C     HFLWSS=snowpack heat flux
C     FLW0S,FLQ0I,FLQ0W=snow,ice,water fluxes through snowpack
C     HFLW0W=convective heat flux snow,water,ice fluxes
C
      IF(VOLS1(L2,NY,NX).GT.ZEROS2(NY,NX))THEN
      DENSW2=AMIN1(0.6,(VOLS0M(L2,NY,NX)+VOLW0M(L2,NY,NX)
     2+VOLI0M(L2,NY,NX)*DENSI)/VOLS1(L2,NY,NX))
      ELSE
      DENSW2=DENS0(NY,NX)
      ENDIF
      TCND2W=0.0036*10**(2.650*DENSW2-1.652)
      ATCNDW=2.0*TCND1W*TCND2W/(TCND1W*DLYRS0(L2,NY,NX)
     2+TCND2W*DLYRS0(L,NY,NX))
      TKY=(TK0M(L,NY,NX)*VHCPWMM(L,NY,NX)+TK0M(L2,NY,NX)
     2*VHCPWMM(L2,NY,NX))/(VHCPWMM(L,NY,NX)+VHCPWMM(L2,NY,NX))
      HFLWX=(TK0M(L,NY,NX)-TKY)*VHCPWMM(L,NY,NX)*XNPA
      HFLWC=ATCNDW*(TK0M(L,NY,NX)-TK0M(L2,NY,NX))
     2*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*XNPY
      IF(HFLWC.GE.0.0)THEN
      HFLWSS=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWSS=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
      FLW0T=FLWQM+FLVSS
      HFLW0T=HFLWQM+HFLVSS+HFLWSS
      FLW0S(L2,NY,NX)=0.0
      FLW0W(L2,NY,NX)=FLW0T
      FLW0I(L2,NY,NX)=0.0
      HFLW0W(L2,NY,NX)=HFLW0T
      FLQWM(M,L2,NY,NX)=FLQWM(M,L2,NY,NX)+FLWQM
C     IF(NX.EQ.3.AND.NY.EQ.3.AND.L.EQ.1)THEN
C     WRITE(*,7757)'FLW0',I,J,M,MM,L2,FLW0W(L2,NY,NX),FLW0T,FLWQM,FLVSS
C    2,HFLW0W(L2,NY,NX),HFLW0T,HFLWQM,HFLVSS,HFLWSS,VP1,VP2,FLVX
C    3,TCND1W,TCND2W,DENSW1,DENSW2,HFLXW2,VHCPWM2,THETP2,FLWQX
C    2,VOLS0M(L,NY,NX),VOLW0M(L,NY,NX),TK0M(L2,NY,NX),TK0M(L,NY,NX)
7757  FORMAT(A8,5I4,30E14.6)
C     ENDIF
C
C     DISCHARGE OF MELTWATER AND ITS HEAT FROM LOWEST SNOWPACK LAYER
C     TO RESIDUE, SURFACE SOIL MICROPORES AND MACROPORES
C
C     FLWQX,FLWQR=porosity-unconstrained water flux to soil,litter
C     FLWQGX,FLWQGS,FLWQGH=water flux to soil surface,
C     micropores,macropores
C     VOLP1,VOLPH1=air volumes of soil micropores,macropores
C     FMAC,FGRD=macropore,micropore volume fractions
C     HFLWQG,HFLWQR=convective heat fluxes to soil,litter
C     THETWR,THETW1=litter, soil water concentration
C     VOLWRX=litter water retention capacity
C     PSISM1(0,PSISM1(NUM=litter,soil water potentials
C     THETY=hygroscopic water concentration
C     POROS=soil porosity
C     FC,WP,FCL,WPL=field capacity,wilting point, log(FC),log(WP)
C     FCI,WPI=FC,WP of ice
C     THETIX=ice concentration
C     BKVL=bulk density x volume of soil layer
C
      ELSE
      IF(ICHKL.EQ.0)THEN
      FLWQGX=FLWQX*BARE(NY,NX)
      FLWQGS=AMIN1(VOLP1(NUM(NY,NX),NY,NX)*XNPX
     2,FLWQGX*FGRD(NUM(NY,NX),NY,NX))
      FLWQGH=AMIN1(VOLPH1(NUM(NY,NX),NY,NX)*XNPX
     2,FLWQGX*FMAC(NUM(NY,NX),NY,NX))
      FLWQG=FLWQGS+FLWQGH
      HFLWQG=4.19*TK0M(L,NY,NX)*FLWQG
      FLWQR=FLWQX-FLWQG
      HFLWQR=4.19*TK0M(L,NY,NX)*FLWQR
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX)
     2.AND.VOLW1(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
      IF(THETWR.LT.FC(0,NY,NX))THEN
      PSISM1(0,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCL(0,NY,NX)-LOG(THETWR))
     3/FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETWR.LT.POROS0(NY,NX))THEN
      PSISM1(0,NY,NX)=-EXP(PSIMS(NY,NX)
     2+(((PSL(0,NY,NX)-LOG(THETWR))
     3/PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
      ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
      ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
      THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX)
     2,AMIN1(POROS(NUM(NY,NX),NY,NX)
     2,VOLW1(NUM(NY,NX),NY,NX)/VOLY(NUM(NY,NX),NY,NX)))
      IF(BKVL(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      IF(THETW1.LT.FC(NUM(NY,NX),NY,NX))THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCL(NUM(NY,NX),NY,NX)-LOG(THETW1))
     3/FCD(NUM(NY,NX),NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETW1.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX)
     2+(((PSL(NUM(NY,NX),NY,NX)-LOG(THETW1))
     3/PSD(NUM(NY,NX),NY,NX))**SRP(NUM(NY,NX),NY,NX)*PSISD(NY,NX)))
      ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
      ELSEIF(VOLX(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      FCX=FCI*THETIX(NUM(NY,NX),NY,NX)
      WPX=WPI*THETIX(NUM(NY,NX),NY,NX)
      FCLX=LOG(FCX)
      WPLX=LOG(WPX)
      PSDX=PSL(L,NY,NX)-FCLX
      FCDX=FCLX-WPLX
      IF(THETWX(NUM(NY,NX),NY,NX).LT.FCX)THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCLX-LOG(THETWX(NUM(NY,NX),NY,NX)))
     3/FCDX*PSIMD(NY,NX))))
      ELSEIF(THETWX(NUM(NY,NX),NY,NX)
     2.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX)
     2+(((PSL(NUM(NY,NX),NY,NX)-LOG(THETWX(NUM(NY,NX),NY,NX)))
     3/PSDX)*PSISD(NY,NX)))
      ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
C     WRITE(*,1119)'PSISMS',I,J,M,N,NX,NY,NUM(NY,NX)
C    2,PSISM(NUM(NY,NX),NY,NX),THETW1,VOLW1(NUM(NY,NX),NY,NX)
C    3,VOLX(NUM(NY,NX),NY,NX)
C    2,THETW(NUM(NY,NX),NY,NX),THETI(NUM(NY,NX),NY,NX)
C    3,FCX,WPX,POROS(NUM(NY,NX),NY,NX)
1119  FORMAT(A8,7I4,20E12.4)
      ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
      PSISV1=PSISM1(NUM(NY,NX),NY,NX)+PSISO(NUM(NY,NX),NY,NX)
C
C     VAPOR FLUX BETWEEN SNOWPACK AND SOIL SURFACE
C
C     VOLP01,THETPM=air volume,concentration
C     CNV1,CNV2=vapor conductances of source, destination layers
C     VP1,VP2=vapor concentrations of source, destination layers
C     POROS,POROQ=porosity, tortuosity
C     WGSGL=vapor diffusivity
C     TK0M,TK1=snow,soil surface temperature
C     PSISV1=soil matric+osmotic potential
C     ATCNVS=snow-soil vapor conductance
C     DLYR=soil surface layer depth
C     FLVC,FLVX=vapor flux unlimited,limited by vapor
C     VPY=equilibrium vapor concentration
C     XNPX=time step for flux calculations
C     FLVS1,HFLVS1=vapor flux and its convective heat flux
C
      IF(VOLP01.GT.ZEROS2(NY,NX)
     2.AND.THETPM(M,NUM(NY,NX),NY,NX).GT.THETX)THEN
      CNV2=WGSGL(NUM(NY,NX),NY,NX)*THETPM(M,NUM(NY,NX),NY,NX)*POROQ
     2*THETPM(M,NUM(NY,NX),NY,NX)/POROS(NUM(NY,NX),NY,NX)
      VP2=2.173E-03/TK1(NUM(NY,NX),NY,NX)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK1(NUM(NY,NX),NY,NX)))
     3*EXP(18.0*PSISV1/(8.3143*TK1(NUM(NY,NX),NY,NX)))
      ATCNVS=2.0*CNV1*CNV2
     2/(CNV1*DLYR(3,NUM(NY,NX),NY,NX)+CNV2*DLYRS0(L,NY,NX))
      FLVC=ATCNVS*(VP1-VP2)*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNW(NY,NX)*BARE(NY,NX)*XNPY
      VPY=(VP1*VOLP01+VP2*VOLPM(M,NUM(NY,NX),NY,NX))
     2/(VOLP01+VOLPM(M,NUM(NY,NX),NY,NX))
      FLVX=(VP1-VPY)*VOLP01*XNPA
      IF(FLVC.GE.0.0)THEN
      FLVS1=AMAX1(0.0,AMIN1(FLVC,FLVX))
      HFLVS1=(4.19*TK0M(L,NY,NX)+VAP)*FLVS1
      ELSE
      FLVS1=AMIN1(0.0,AMAX1(FLVC,FLVX))
      HFLVS1=(4.19*TK1(NUM(NY,NX),NY,NX)+VAP)*FLVS1
      ENDIF
      ELSE
      CNV2=0.0
      FLVS1=0.0
      HFLVS1=0.0
      ENDIF
C
C     HEAT FLUX BETWEEN SNOWPACK AND SURFACE SOIL
C
C     WTHET2=multiplier for air concentration in thermal conductivity
C     TCND1W,TCNDS=thermal conductivity of snowpack, soil surface
C     STC,DTC=mineral component of thermal conductivity
C     THETWX,THETIX,THETPX=soil surface water,ice,air concentrations
C     BAREW=soil surface fraction
C     ATCNDS=snowpack-soil thermal conductance
C     TKWX1=interim snowpack temperature
C     TKY=equilibrium temperature
C     HFLWX,HFLWC=heat-constrained,heat-unconstrained heat fluxes
C     XNPY=time step for snowpack flux calculations
C     HFLWS1=snowpack-soil heat flux
C
      WTHET2=1.467-0.467*THETPY(NUM(NY,NX),NY,NX)
      TCNDS=(STC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX)
     2*2.067E-03+0.611*THETIX(NUM(NY,NX),NY,NX)*7.844E-03
     3+WTHET2*THETPX(NUM(NY,NX),NY,NX)*9.050E-05)
     4/(DTC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX)
     5+0.611*THETIX(NUM(NY,NX),NY,NX)
     6+WTHET2*THETPX(NUM(NY,NX),NY,NX))
      IF(BARE(NY,NX).GT.ZERO)THEN
      ATCNDS=2.0*TCND1W*TCNDS/(TCND1W*DLYR(3,NUM(NY,NX),NY,NX)
     2+TCNDS*DLYRS0(L,NY,NX))
      ELSE
      ATCNDS=0.0
      ENDIF
      TKWX1=TK1(NUM(NY,NX),NY,NX)+HFLVS1/VHCP1(NUM(NY,NX),NY,NX)
      TKY=(TK0M(L,NY,NX)*VHCPWMM(L,NY,NX)+TKWX1*VHCP1(NUM(NY,NX),NY,NX))
     2/(VHCPWMM(L,NY,NX)+VHCP1(NUM(NY,NX),NY,NX))
      HFLWX=(TK0M(L,NY,NX)-TKY)*VHCPWMM(L,NY,NX)*XNPA
      HFLWC=ATCNDS*(TK0M(L,NY,NX)-TKWX1)*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNW(NY,NX)*BARE(NY,NX)*XNPY
      IF(HFLWC.GE.0.0)THEN
      HFLWS1=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWS1=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
C     IF(J.EQ.15.AND.M.EQ.NPH)THEN
C     WRITE(*,1113)'HFLWS1',I,J,M,MM,L,FLVS1,FLVX,HFLVS1
C    2,HFLWS1,ATCNVS,VP1,VP2,CNV1,CNV2,PSISV1
C    3,HFLWX,HFLWC,ATCNDS,TKW(L,NY,NX),TK1(NUM(NY,NX),NY,NX)
C    4,THETPX(NUM(NY,NX),NY,NX),WGSGL(NUM(NY,NX),NY,NX)
C    5,VHCPWMM(L,NY,NX),TCND1W,TCNDS,PSISV1,TKY,TK0M(L,NY,NX),TKWX1
C    6,VOLP1(NUM(NY,NX),NY,NX),VOLPH1(NUM(NY,NX),NY,NX)
C    6,VOLT(NUM(NY,NX),NY,NX),VOLA1(NUM(NY,NX),NY,NX)
C    7,VOLW1(NUM(NY,NX),NY,NX),VOLI1(NUM(NY,NX),NY,NX)
C    8,POROS(NUM(NY,NX),NY,NX)
1113  FORMAT(A8,5I4,60E14.6)
C     ENDIF
C
C     HEAT FLUX AMONG SNOWPACK, SURFACE RESIDUE AND SURFACE SOIL
C
C     FLVSR=snowpack-litter vapor flux
C     HFLVSR,HFLWSR=snowpack-litter convective,conductive heat fluxes
C     FLVS1=snowpack-soil vapor flux
C     HFLVS1,HFLWS1=snowpack-soil convective,conductive heat fluxes
C     VHCP1,VHCPRX=current,minimum litter heat capacities
C     TK0X,TKXR,TK1X=snowpack,litter,soil temperatures
C     CNVR,CNV1,CNV2=litter,snowpack,soil vapor conductivity
C     THETP*,THETWX,THETIX=litter air,water,ice concentration
C     POROS,POROQ=litter porosity, tortuosity
C     CVRD=litter cover fraction
C     WGSGR=litter vapor diffusivity
C     ATCNVR,ATCNVS=snowpack-litter,litter-soil vapor conductance
C     DLYRR,DLYRS0,DLYR=litter,snowpack,soil depths
C     THETRR=dry litter concentration
C     TCNDR,TCND1W,TCNDS=litter,snowpack,soil thermal conductivity
C     ATCNDR,ATCNDS=snow-litter,litter-soil thermal conductance
C
      FLVSR=0.0
      HFLVSR=0.0
      HFLWSR=0.0
      FLVR1=0.0
      HFLVR1=0.0
      HFLWR1=0.0
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      TK0X=TK0M(L,NY,NX)
      TKXR=TK1(0,NY,NX)
      TK1X=TK1(NUM(NY,NX),NY,NX)
      CNVR=WGSGR(NY,NX)*THETPM(M,0,NY,NX)*POROQ
     2*THETPM(M,0,NY,NX)/POROS(0,NY,NX)
      IF(CVRD(NY,NX).GT.ZERO)THEN
      IF(CNV1.GT.ZERO.AND.CNVR.GT.ZERO)THEN
      ATCNVR=2.0*CNVR*CNV1/(CNV1*DLYRR(NY,NX)+CNVR*DLYRS0(L,NY,NX))
      ELSE
      ATCNVR=2.0*CNV1/(DLYRR(NY,NX)+DLYRS0(L,NY,NX))
      ENDIF
      IF(CNVR.GT.ZERO.AND.CNV2.GT.ZERO)THEN
      ATCNVS=2.0*CNVR*CNV2
     2/(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+CNV2*DLYRR(NY,NX))
      ELSE
      ATCNVS=2.0*CNV2/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR(NY,NX))
      ENDIF
      THETRR=AMAX1(0.0,1.0-THETPX(0,NY,NX)-THETWX(0,NY,NX)
     2-THETIX(0,NY,NX))
      TCNDR=(0.779*THETRR*9.050E-04+0.622*THETWX(0,NY,NX)
     2*2.067E-03+0.380*THETIX(0,NY,NX)*7.844E-03+THETPX(0,NY,NX)
     3*9.050E-05)/(0.779*THETRR+0.622*THETWX(0,NY,NX)
     4+0.380*THETIX(0,NY,NX)+THETPX(0,NY,NX))
      IF(TCND1W.GT.ZERO.AND.TCNDR.GT.ZERO)THEN
      ATCNDR=2.0*TCND1W*TCNDR
     2/(TCND1W*DLYRR(NY,NX)+TCNDR*DLYRS0(L,NY,NX))
      ELSE
      ATCNDR=0.0
      ENDIF
      IF(TCNDR.GT.ZERO.AND.TCNDS.GT.ZERO)THEN
      ATCNDS=2.0*TCNDR*TCNDS
     2/(TCNDR*DLYR(3,NUM(NY,NX),NY,NX)+TCNDS*DLYRR(NY,NX))
      ELSE
      ATCNDS=0.0
      ENDIF
      ELSE
      ATCNVR=0.0
      ATCNVS=0.0
      ATCNDR=0.0
      ATCNDS=0.0
      ENDIF
C
C     SHORTER TIME STEP FOR SURFACE RESIDUE FLUX CALCULATIONS
C
      DO 4000 NN=1,NPR
C
C     VAPOR FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
C
C     VP0,VPR,VPY=snowpack,litter, equilibrium vapor concentration
C     TK0X,TKXR=snowpack,litter temperature
C     PSISM1=litter matric water potential
C     FLVC,FLVX=vapor-unconstrained,vapor-constrained vapor flux
C     AREA=area of grid cell
C     FSNW,CVRD=snow,litter cover fraction
C     XNPQ=time step for flux calculation
C     FLVRSX=snow-litter vapor flux
C     HFLVRSX=convective heat flux from snow-litter vapor flux
C
      VPR=2.173E-03/TKXR
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKXR))
     3*EXP(18.0*PSISM1(0,NY,NX)/(8.3143*TKXR))
      IF(VOLP01.GT.ZEROS2(NY,NX)
     2.AND.THETPM(M,0,NY,NX).GT.THETX)THEN
      VP0=2.173E-03/TK0X
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK0X))
      FLVC=ATCNVR*(VP0-VPR)*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ
      VPY=(VP0*VOLP01+VPR*VOLPM(M,0,NY,NX))
     2/(VOLP01+VOLPM(M,0,NY,NX))
      FLVX=(VP0-VPY)*VOLP01*XNPC
      IF(FLVC.GE.0.0)THEN
      FLVSRX=AMAX1(0.0,AMIN1(FLVC,FLVX))
      HFLVSRX=(4.19*TK0X+VAP)*FLVSRX
      ELSE
      FLVSRX=AMIN1(0.0,AMAX1(FLVC,FLVX))
      HFLVSRX=(4.19*TKXR+VAP)*FLVSRX
      ENDIF
      ELSE
      FLVSRX=0.0
      HFLVSRX=0.0
      ENDIF
C
C     HEAT FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
C
C     TKY=snow-litter equilibrium temperature
C     HFLWC,HFLWX=snow-litter heat flux unltd,ltd by heat
C     HFLWSRX=snow-litter heat flux
C
      TKY=(TK0X*VHCPWMM(L,NY,NX)+TKXR*VHCP1(0,NY,NX))
     2/(VHCPWMM(L,NY,NX)+VHCP1(0,NY,NX))
      HFLWX=(TK0X-TKY)*VHCPWMM(L,NY,NX)*XNPC
      HFLWC=ATCNDR*(TK0X-TKXR)*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ
      IF(HFLWC.GE.0.0)THEN
      HFLWSRX=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWSRX=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
C     WRITE(*,7752)'TKXR',I,J,M,MM,NX,NY,L
C    2,FLVC,FLVX,VP1,VPR,VPY,ATCNVS,FSNW(NY,NX),BARE(NY,NX)
C    3,VOLP01,VOLPM(M,NUM(NY,NX),NY,NX),TK1X,TKXR
C    4,HFLVR1X,VHCP1(NUM(NY,NX),NY,NX),VHCP1(0,NY,NX)
C    3,VOLPM(M,0,NY,NX),VOLPM(M,NUM(NY,NX),NY,NX),VOLP01
C    2,THETPM(M,0,NY,NX),TKXR
C
C     VAPOR FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
C
C     THETPM,VOLPM=air-filled porosity,volume
C     VP1,VPY=soil,litter-soil equilibrium vapor concentration
C     TK1X=soil temperature
C     PSISV1=soil matric+osmotic water potentials
C     FLVC,FLVX=vapor-unconstrained,vapor-constrained vapor flux
C     FLVR1X=litter-soil vapor flux
C     HFLVR1X=convective heat of litter-soil vapor flux
C     TKXR,TK1X=interim calculation of litter,soil temperatures
C
      IF(VOLPM(M,0,NY,NX).GT.ZEROS(NY,NX)
     2.AND.VOLPM(M,NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      VP1=2.173E-03/TK1X
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK1X))
     3*EXP(18.0*PSISV1/(8.3143*TK1X))
      FLVC=ATCNVS*(VPR-VP1)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)
     2*CVRD(NY,NX)*XNPQ
      VPY=(VPR*VOLPM(M,0,NY,NX)+VP1*VOLPM(M,NUM(NY,NX),NY,NX))
     2/(VOLPM(M,0,NY,NX)+VOLPM(M,NUM(NY,NX),NY,NX))
      FLVX=(VPR-VPY)*VOLPM(M,0,NY,NX)*XNPC
      IF(FLVC.GE.0.0)THEN
      FLVR1X=AMAX1(0.0,AMIN1(FLVC,FLVX,VOLW0M(L,NY,NX)*XNPB))
      HFLVR1X=(4.19*TKXR+VAP)*FLVR1X
      ELSE
      FLVR1X=AMIN1(0.0,AMAX1(FLVC,FLVX))
      HFLVR1X=(4.19*TK1X+VAP)*FLVR1X
      ENDIF
      ELSE
      FLVR1X=0.0
      HFLVR1X=0.0
      ENDIF
      TKXR=TKXR-HFLVR1X/VHCP1(0,NY,NX)
      TK1X=TK1X+HFLVR1X/VHCP1(NUM(NY,NX),NY,NX)
C     IF(NY.EQ.6)THEN
C     WRITE(*,7752)'TK1X',I,J,M,MM,NX,NY,NN
C    2,FLVC,FLVX,VP1,VPR,VPY,ATCNVS,FSNW(NY,NX),BARE(NY,NX)
C    3,VOLP01,VOLPM(M,NUM(NY,NX),NY,NX),TK1X,TKXR,TK1(0,NY,NX)
C    4,HFLVR1X,VHCP1(NUM(NY,NX),NY,NX),VHCP1(0,NY,NX)
C    3,VOLPM(M,0,NY,NX),VOLPM(M,NUM(NY,NX),NY,NX),VOLP01
C    2,THETPM(M,0,NY,NX)
C     ENDIF
C
C     HEAT FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
C
C     TKY=litter-soil equilibrium temperature
C     HFLWC,HFLWX=litter-soil heat flux unltd,ltg by heat
C     HFLWR1X=litter-soil heat flux
C
      TKY=(TKXR*VHCP1(0,NY,NX)+TK1X*VHCP1(NUM(NY,NX),NY,NX))
     2/(VHCP1(0,NY,NX)+VHCP1(NUM(NY,NX),NY,NX))
      HFLWX=(TKXR-TKY)*VHCP1(0,NY,NX)*XNPC
      HFLWC=ATCNDS*(TKXR-TK1X)*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNW(NY,NX)*CVRD(NY,NX)*XNPQ
      IF(HFLWC.GE.0.0)THEN
      HFLWR1X=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWR1X=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
C
C     ACCUMULATE SNOW-LITTER, LITTER-SOIL HEAT FLUXES
C     WITHIN LONGER TIME STEP FOR SNOWPACK FLUX CALCULATIONS
C
      FLVSR=FLVSR+FLVSRX
      HFLVSR=HFLVSR+HFLVSRX
      HFLWSR=HFLWSR+HFLWSRX
      FLVR1=FLVR1+FLVR1X
      HFLVR1=HFLVR1+HFLVR1X
      HFLWR1=HFLWR1+HFLWR1X
      TK0X=TK0X-HFLVSRX/VHCPWMM(L,NY,NX)
      TKXR=TKXR+(HFLVSRX-HFLWR1X)/VHCP1(0,NY,NX)
      TK1X=TK1X+HFLWR1X/VHCP1(NUM(NY,NX),NY,NX)
C     IF(I.EQ.53)THEN
C     WRITE(*,1114)'FLVR0',I,J,M,N,NX,NY
C    2,FLVSR,FLVSRX,FLVC,FLVX,ATCNVR,VP0,VPR,VPY,FSNW(NY,NX),CVRD(NY,NX)
C    3,XNPQ,THETPM(M,0,NY,NX),VOLP01,XNPX,XNPV
C    2,TK0M(1,NY,NX),TK1(0,NY,NX)
C    2,TK1(NUM(NY,NX),NY,NX),TK0X,TKXR,TK1X,FLVR1,HWFLVR1,FLVS1
C    4,HWFLVS1,HFLC0R1,HFLCR11,FLVR,HWFLVR,FLVS,HWFLVS
C    3,HFLC0R,HFLCR1,VPQ(NY,NX),VP0,VPR,VP1,PSISM1(0,NY,NX),PSISV1
C    5,AVCNVR,ATCNDR,AVCNVS,ATCNDS,VHCPWMM(L,NY,NX),VHCP1(0,NY,NX)
C    6,VHCP1(NUM(NY,NX),NY,NX),DLYRR(NY,NX),DLYRS0(NY,NX),CNV01,CNVR1
C    7,CNV11,CNV1,THETPX(NUM(NY,NX),NY,NX),POROQ
C    2,WGSGL(NUM(NY,NX),NY,NX),CVRD(NY,NX),HFLXR,HFLVS1,HFLCR1
1114  FORMAT(A8,6I4,60E12.4)
C     ENDIF
4000  CONTINUE
      ENDIF
C
C     GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
C     FOR LATER UPDATES TO STATE VARIABLES
C
C     FLWLT,FLWLW=total,accumulated water flux to soil micropores
C     FLWLXW,FLWHLW=total,accumd snow-soil micropore,macropore water
C     HFLWLT,HFLWLW=total,accumulated snow+litter heat flux to soil
C     FLWRT,FLWRLW=total,accumulated snow+soil water flux to litter
C     HFLWRT,HFLWRLW=total,accumulated snow+soil heat flux to litter
C     FLQRM,FLQSM,FLQHM=total water flux to litter,soil micropore,macropore
C     FLSW,FLSWH,FLSWR=water flux from lowest snow layer to soil macropore,micropore,litter
C     HFLSW,HFLSWR=heat flux from lowest snow layer to soil,litter
C
      FLWLT=FLWQGS+FLVS1+FLVR1
      FLWLW=FLWLW+FLWLT
      FLWLXW=FLWLXW+FLWQGS
      FLWHLW=FLWHLW+FLWQGH
      HFLWLT=HFLWQG+HFLVS1+HFLWS1+HFLVR1+HFLWR1
      HFLWLW=HFLWLW+HFLWLT
      FLWRT=FLWQR+FLVSR-FLVR1
      FLWRLW=FLWRLW+FLWRT
      HFLWRT=HFLWQR+HFLVSR+HFLWSR-HFLVR1-HFLWR1
      HFLWRLW=HFLWRLW+HFLWRT
      FLWVLW=0.0
      FLQRM(M,NY,NX)=FLQRM(M,NY,NX)+FLWQR
      FLQSM(M,NY,NX)=FLQSM(M,NY,NX)+FLWQGS
      FLQHM(M,NY,NX)=FLQHM(M,NY,NX)+FLWQGH
      FLSW(L,NY,NX)=FLSW(L,NY,NX)+FLWLT
      FLSWH(L,NY,NX)=FLSWH(L,NY,NX)+FLWQGH
      HFLSW(L,NY,NX)=HFLSW(L,NY,NX)+HFLWLT
      FLSWR(L,NY,NX)=FLSWR(L,NY,NX)+FLWRT
      HFLSWR(L,NY,NX)=HFLSWR(L,NY,NX)+HFLWRT
C     IF(I.EQ.53)THEN
C     WRITE(*,7752)'FLWLW',I,J,M,MM,NX,NY,L
C    2,FLWLW,FLWLT,FLWQGS,FLVS1,FLVR1,FLSW(L,NY,NX)
C    2,FLWRLW,FLWRT,FLWQR,FLVSR,FLVR1,FLWQX,FLWQG,FLSWR(L,NY,NX)
C    2,FLWQX,FLWQGX,BARE(NY,NX),VOLW0M(L,NY,NX),VOLS0M(L,NY,NX)
C    2,HFLWLW,HFLWLT,HFLWQG,HFLVS1,HFLWS1,HFLVR1,HFLWR1
C    3,HFLSW(L,NY,NX)
C    3,HFLWRLW,HFLWRT,HFLWQR,HFLVSR,HFLWSR,HFLVR1,HFLWR1
C    3,HFLSWR(L,NY,NX)
C    2,VOLP01,THETPM(M,NUM(NY,NX),NY,NX),THETX
C    3,CNV2,VP2,TK1(NUM(NY,NX),NY,NX),ATCNVS,FLVC,VPY,FLVX
C    4,VP1,VP2,TK1X,PSISV1,HFLVR1X,HFLWR1X
C    5,VHCP1(NUM(NY,NX),NY,NX)
C    3,HFLWRLW,HFLWLW,VP0,VPR,VPY
C    3,THETPX(NUM(NY,NX),NY,NX),FLWQX,BARE(NY,NX)
C    4,VOLW0M(L,NY,NX),VOLS0M(L,NY,NX)
C    2,HFLWX,HFLWC,ATCNDS,TK0M(L,NY,NX),TKWX1
C    2,TCND1W,TCNDS,DLYR(3,NUM(NY,NX),NY,NX),DLYRS0(L,NY,NX)
C    2,THETWX(NUM(NY,NX),NY,NX),THETIX(NUM(NY,NX),NY,NX)
C    3,WTHET2,THETPX(NUM(NY,NX),NY,NX),VOLP1(NUM(NY,NX),NY,NX)
C    2,VOLPH1(NUM(NY,NX),NY,NX),VOLA1(NUM(NY,NX),NY,NX)
C    2,VOLW1(NUM(NY,NX),NY,NX),VOLI1(NUM(NY,NX),NY,NX),BARE(NY,NX)
C    3,FLQRM(M,NY,NX),FLWQR,FLQSM(M,NY,NX),FLWQGS
C    4,FLQHM(M,NY,NX),FLWQG,VOLS0M(L,NY,NX),VOLW0M(L,NY,NX)
C    5,VOLI0M(L,NY,NX),DLYRS0(L,NY,NX),FLWQX,TK0M(L,NY,NX)
7752  FORMAT(A8,7I4,40E12.4)
C     ENDIF
      ICHKL=1
      ENDIF
      ENDIF
      ENDIF
9880  CONTINUE
C
C     ACCUMULATE SNOWPACK FLUXES TO LONGER TIME STEP FOR
C     LITTER, SOIL FLUX CALCULATIONS
C
C     THRMG=total longwave emission
C     XFLWS,XFLWW,XFLWI=hourly accumulated snow,water,ice transfer
C     XHFLWW=hourly convective heat flux from snow,water,ice transfer
C     TFLWSX,TFLWWX,TFLWIX=net snow,water,ice transfer
C     THFLWWX=convective heat flux from net snow,water,ice transfer
C     TFLWS,TFLWW,TFLWI=accumulated net snow,water,ice transfer
C     THFLWW=convective heat flux from accumd snow,water,ice transfer
C
      THRMG(NY,NX)=THRMG(NY,NX)+THRMX
      DO 9860 L=1,JS
      XFLWS(L,NY,NX)=XFLWS(L,NY,NX)+FLW0S(L,NY,NX)
      XFLWW(L,NY,NX)=XFLWW(L,NY,NX)+FLW0W(L,NY,NX)
      XFLWI(L,NY,NX)=XFLWI(L,NY,NX)+FLW0I(L,NY,NX)
      XHFLWW(L,NY,NX)=XHFLWW(L,NY,NX)+HFLW0W(L,NY,NX)
      L2=MIN(JS,L+1)
C
C     IF WITHIN SNOWPACK
C
      IF(L.LT.JS.AND.VHCPWMM(L2,NY,NX).GT.VHCPWX(NY,NX))THEN
      TFLWSX=FLW0S(L,NY,NX)-FLW0S(L2,NY,NX)
      TFLWWX=FLW0W(L,NY,NX)-FLW0W(L2,NY,NX)
      TFLWIX=FLW0I(L,NY,NX)-FLW0I(L2,NY,NX)
      THFLWWX=HFLW0W(L,NY,NX)-HFLW0W(L2,NY,NX)
      TFLWS(L,NY,NX)=TFLWS(L,NY,NX)+TFLWSX
      TFLWW(L,NY,NX)=TFLWW(L,NY,NX)+TFLWWX
      TFLWI(L,NY,NX)=TFLWI(L,NY,NX)+TFLWIX
      THFLWW(L,NY,NX)=THFLWW(L,NY,NX)+THFLWWX
C
C     IF AT BOTTOM OF SNOWPACK
C
      ELSEIF(VHCPWMM(L,NY,NX).GT.VHCPWX(NY,NX))THEN
      TFLWSX=FLW0S(L,NY,NX)
      TFLWWX=FLW0W(L,NY,NX)-FLWRT-FLWLT-FLWQGH
      TFLWIX=FLW0I(L,NY,NX)
      THFLWWX=HFLW0W(L,NY,NX)-HFLWRT-HFLWLT
      TFLWS(L,NY,NX)=TFLWS(L,NY,NX)+TFLWSX
      TFLWW(L,NY,NX)=TFLWW(L,NY,NX)+TFLWWX
      TFLWI(L,NY,NX)=TFLWI(L,NY,NX)+TFLWIX
      THFLWW(L,NY,NX)=THFLWW(L,NY,NX)+THFLWWX
      ELSE
      TFLWSX=0.0
      TFLWWX=0.0
      TFLWIX=0.0
      THFLWWX=0.0
      ENDIF
C     IF(L.EQ.1)THEN
C     WRITE(*,7763)'TFLWW',I,J,M,MM,NX,NY,L,L2,TFLWW(L,NY,NX)
C    2,TFLWWX,FLW0W(L,NY,NX),FLW0W(L2,NY,NX),FLWRLW,FLWLW,FLWHLW
C    3,VOLW0(L,NY,NX),FLWRT,FLWLT,FLWQGH
C    2,THFLWW(L,NY,NX),THFLWWX
C    3,HFLW0W(L,NY,NX),HFLW0W(L2,NY,NX),HFLWRLW,HFLWLW
C    4,VHCPWMM(L,NY,NX)
7763  FORMAT(A8,8I4,30E14.6)
C     ENDIF
C
C     FREEZE-THAW IN SNOWPACK FROM NET CHANGE IN SNOWPACK
C     HEAT STORAGE
C
C     VOLS0M,VOLW0M,VOLI0M=snow,water,ice volume
C     VHCPWMM,VHCPWMX,VHCPWX=previous,current,minimum heat capacity
C     TK0M=snowpack temperature
C     THFLWWX=net conductive+convective heat flux
C     TFLX1=unconstrained latent heat flux from freeze-thaw
C     FVOLS0,FVOLI0=fractions of total water in water,ice
C     TFLX0X=source-limited latent heat flux from freeze-thaw
C     WFLXSX,WFLXIX=freeze-thaw changes in water,ice
C     WFLXS,WFLXI=accumulated freeze-thaw
C     TFLX0=accumulated latent heat flux from freeze-thaw
C     XWFLXS,XWFLXI=hourly accumulated freeze-thaw
C     XTHAWW=hourly accumulated latent heat flux from freeze-thaw
C
      VOLS0X=AMAX1(0.0,VOLS0M(L,NY,NX))
      VOLW0X=AMAX1(0.0,VOLW0M(L,NY,NX))
      VOLI0X=AMAX1(0.0,VOLI0M(L,NY,NX))
      ENGY0=VHCPWMM(L,NY,NX)*TK0M(L,NY,NX)
      VHCPWMX=2.095*VOLS0X+4.19*VOLW0X+1.9274*VOLI0X
      IF(VHCPWMX.GT.VHCPWX(NY,NX))THEN
      TK0X=(ENGY0+THFLWWX)/VHCPWMX
      IF((TK0X.LT.273.15
     2.AND.VOLW0X.GT.ZERO*VOLS(NY,NX))
     3.OR.(TK0X.GT.273.15
     4.AND.VOLI0X+VOLS0X.GT.ZERO*VOLS(NY,NX)))THEN
      TFLX1=VHCPWMX*(273.15-TK0X)/2.7185*XNPX
      IF(TFLX1.LT.0.0)THEN
      TVOLWS=VOLS0X+VOLI0X*DENSI
      IF(TVOLWS.GT.ZEROS2(NY,NX))THEN
      FVOLS0=VOLS0X/TVOLWS
      FVOLI0=VOLI0X*DENSI/TVOLWS
      ELSE
      FVOLS0=0.0
      FVOLI0=0.0
      ENDIF
      TFLX0X=AMAX1(-333.0*TVOLWS*XNPX,TFLX1)
      WFLXSX=-TFLX0X*FVOLS0/333.0
      WFLXIX=-TFLX0X*FVOLI0/333.0
      ELSE
      FVOLS0=0.0
      FVOLI0=0.0
      TFLX0X=AMIN1(333.0*VOLW0X*XNPX,TFLX1)
      WFLXSX=0.0
      WFLXIX=-TFLX0X/333.0
      ENDIF
      ELSE
      TFLX1=0.0
      FVOLS0=0.0
      FVOLI0=0.0
      TFLX0X=0.0
      WFLXSX=0.0
      WFLXIX=0.0
      ENDIF
      WFLXS(L,NY,NX)=WFLXS(L,NY,NX)+WFLXSX
      WFLXI(L,NY,NX)=WFLXI(L,NY,NX)+WFLXIX
      TFLX0(L,NY,NX)=TFLX0(L,NY,NX)+TFLX0X
      XWFLXS(L,NY,NX)=XWFLXS(L,NY,NX)+WFLXSX
      XWFLXI(L,NY,NX)=XWFLXI(L,NY,NX)+WFLXIX
      XTHAWW(L,NY,NX)=XTHAWW(L,NY,NX)+TFLX0X
      ELSE
      TFLX0X=0.0
      WFLXSX=0.0
      WFLXIX=0.0
      ENDIF
C     IF(L.EQ.5)THEN
C     WRITE(*,7758)'TFLX0',I,J,M,MM,NX,NY,L,TK0M(L,NY,NX)
C    4,TFLX0(L,NY,NX),WFLXS(L,NY,NX),WFLXI(L,NY,NX)
C    4,XTHAWW(L,NY,NX),XWFLXS(L,NY,NX),XWFLXI(L,NY,NX)
C    2,TK0X,TKW(L,NY,NX),VHCPWMX,TFLX1,VOLS0X
C    3,VOLW0X,VOLI0X,TFLWS(L,NY,NX),TFLWW(L,NY,NX)
C    4,TFLWI(L,NY,NX),FVOLS0,FVOLI0
C    5,TFLWW(L,NY,NX),THFLWW(L,NY,NX),FLW0W(L,NY,NX)
7758  FORMAT(A8,7I4,30E14.6)
C     ENDIF
C
C     INTERNAL SNOWPACK SNOW, WATER, ICE, TEMPERATURE
C
C     VOLS0M,VOLW0M,VOLI0M=snow,water,ice volume
C     TFLWSX,TFLWWX,TFLWIX=net snow,water,ice transfer
C     THFLWWX=conductive+convective heat from snow,water,ice transfer
C     WFLXSX,WFLXIX=freeze-thaw changes in water,ice
C     TFLX0X=source-limited latent heat flux from freeze-thaw
C     DENSI=ice density
C     TK0M,TKA=snowpack,air temperature
C     VHCPWMM,VHCPWX=snowpack, minimum heat capacity
C
      VOLS0M(L,NY,NX)=VOLS0M(L,NY,NX)+TFLWSX-WFLXSX
      VOLW0M(L,NY,NX)=VOLW0M(L,NY,NX)+TFLWWX+WFLXSX+WFLXIX
      VOLI0M(L,NY,NX)=VOLI0M(L,NY,NX)-WFLXIX/DENSI
      ENGY0=VHCPWMM(L,NY,NX)*TK0M(L,NY,NX)
      VHCPWMM(L,NY,NX)=2.095*VOLS0M(L,NY,NX)+4.19*VOLW0M(L,NY,NX)
     2+1.9274*VOLI0M(L,NY,NX)
      IF(VHCPWMM(L,NY,NX).GT.VHCPWX(NY,NX))THEN
      TK0M(L,NY,NX)=(ENGY0+THFLWWX+TFLX0X)/VHCPWMM(L,NY,NX)
      ELSEIF(L.EQ.1)THEN
      TK0M(L,NY,NX)=TKA(NY,NX)
      ELSE
      TK0M(L,NY,NX)=TK0M(L-1,NY,NX)
      ENDIF
C     IF(L.EQ.5)THEN
C     WRITE(*,7758)'TK0M',I,J,M,MM,NX,NY,L,TK0M(L,NY,NX)
C    2,THFLWWX,TFLX0X
C    3,VOLS0M(L,NY,NX),VOLW0M(L,NY,NX),VOLI0M(L,NY,NX),VOLS1(L,NY,NX)
C    2,WFLXSX,WFLXIX
C    3,TFLWSX,TFLWWX,TFLWIX
C    3,XFLWS(L,NY,NX),XFLWW(L,NY,NX),XFLWI(L,NY,NX)
C    4,XHFLWW(L,NY,NX),VHCPWMM(L,NY,NX)
C     ENDIF
9860  CONTINUE
3000  CONTINUE
      ENDIF
C
C     ENERGY EXCHANGE AT SOIL SURFACE IF EXPOSED UNDER SNOWPACK
C
      IF(FSNX(NY,NX).GT.0.0.AND.(BKDS(NUM(NY,NX),NY,NX).GT.ZERO
     2.OR.VHCP1(NUM(NY,NX),NY,NX).GT.VHCPNX(NY,NX)))THEN
C
C     PHYSICAL AND HYDRAULIC PROPERTIES OF SOIL SURFACE INCLUDING
C     AIR AND WATER-FILLED POROSITY, AND WATER POTENTIAL USED IN
C     FLUX CALCULATIONS
C
C     THETW1,THETY=current,hygroscopic water concentration
C     POROS=porosity
C     VOLW1,VOLXI=volume of micropore water
C     VOLXI=soil volume less macropore,rock
C     FC,WP=water contents at field capacity,wilting point
C     FCL,WPL=log FC,WP
C     FCD,PSD=FCL-WPL,log(POROS)-FCL
C     PSISM1,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
C     PSIMX,PSIMN,PSIMS=log water potential at FC,WP,POROS
C     PSISD,PSIMD=PSIMX-PSIMS,PSIMN-PSIMX
C     SRP=parameter for deviation from linear log-log water retention
C     function from hour1.f
C     PSISO=osmotic potential
C
      IF(BKVL(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX)
     2,AMIN1(POROS(NUM(NY,NX),NY,NX)
     2,VOLW1(NUM(NY,NX),NY,NX)/VOLY(NUM(NY,NX),NY,NX)))
      IF(THETW1.LT.FC(NUM(NY,NX),NY,NX))THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCL(NUM(NY,NX),NY,NX)-LOG(THETW1))
     3/FCD(NUM(NY,NX),NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETW1.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX)
     2+(((PSL(NUM(NY,NX),NY,NX)-LOG(THETW1))
     3/PSD(NUM(NY,NX),NY,NX))**SRP(NUM(NY,NX),NY,NX)*PSISD(NY,NX)))
      ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
C
C     SURFACE WATER LAYER
C
C     THETIX,THETWX=ice,water concentration
C     FCI,WPI=ice field capacity,wilting point
C     PSISM1=matric water potential
C
      ELSEIF(VOLX(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      FCX=FCI*THETIX(NUM(NY,NX),NY,NX)
      WPX=WPI*THETIX(NUM(NY,NX),NY,NX)
      FCLX=LOG(FCX)
      WPLX=LOG(WPX)
      PSDX=PSL(NU(NY,NX),NY,NX)-FCLX
      FCDX=FCLX-WPLX
      IF(THETWX(NUM(NY,NX),NY,NX).LT.FCX)THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCLX-LOG(THETWX(NUM(NY,NX),NY,NX)))
     3/FCDX*PSIMD(NY,NX))))
      ELSEIF(THETWX(NUM(NY,NX),NY,NX)
     2.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX)
     2+(((PSL(NUM(NY,NX),NY,NX)-LOG(THETWX(NUM(NY,NX),NY,NX)))
     3/PSDX)*PSISD(NY,NX)))
      ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
C     WRITE(*,1119)'PSISMG',I,J,M,N,NX,NY,NUM(NY,NX)
C    2,PSISM(NUM(NY,NX),NY,NX)
C    2,THETW(NUM(NY,NX),NY,NX),THETI(NUM(NY,NX),NY,NX)
C    3,FCX,WPX,POROS(NUM(NY,NX),NY,NX)
      ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
      PSISV1=PSISM1(NUM(NY,NX),NY,NX)+PSISO(NUM(NY,NX),NY,NX)
C     IF(NX.EQ.4.AND.NY.EQ.5)THEN
C     WRITE(*,3232)'PSISV1',I,J,M,NX,NY,NUM(NY,NX),PSISV1
C    2,PSISM1(NUM(NY,NX),NY,NX),PSISO(NUM(NY,NX),NY,NX)
C    3,THETWX(NUM(NY,NX),NY,NX),THETW1,POROS(NUM(NY,NX),NY,NX)
C    4,PSL(NUM(NY,NX),NY,NX),LOG(THETW1),PSD(NUM(NY,NX),NY,NX)
C    5,VOLW1(NUM(NY,NX),NY,NX),VOLY(NUM(NY,NX),NY,NX)
C    5,VOLX(NUM(NY,NX),NY,NX)
C    5,SRP(NUM(NY,NX),NY,NX)
3232  FORMAT(A8,6I4,20E14.6)
C     ENDIF
C
C     SOIL SURFACE ALBEDO, NET RADIATION
C
C     VOLW1,VOLI1=water,ice volume in micopores
C     VOLWH1,VOLIH1=water,ice volume in macopores
C     ALBG,ALBS=albedo of ground surface,soil
C     BKVL=soil mass
C     RADXG,THRYG,RFLXG=incoming shortwave,longwave,net radiation
C     THRMA,THRMS=emitted longwave radiation, emissivity
C     TK1=soil temperature
C
      VOLWXG=VOLW1(NUM(NY,NX),NY,NX)+VOLWH1(NUM(NY,NX),NY,NX)
      VOLIXG=VOLI1(NUM(NY,NX),NY,NX)+VOLIH1(NUM(NY,NX),NY,NX)
      IF(VOLWXG+VOLIXG.GT.ZEROS2(NY,NX))THEN
      ALBG=(ALBS(NY,NX)*BKVL(NUM(NY,NX),NY,NX)+0.06*VOLWXG
     2+0.30*VOLIXG)/(BKVL(NUM(NY,NX),NY,NX)+VOLWXG+VOLIXG)
      ELSE
      ALBG=ALBS(NY,NX)
      ENDIF
      RFLX0=(1.0-ALBG)*RADXG(NY,NX)+THRYG(NY,NX)
      THRMA=THRMS(NY,NX)*TK1(NUM(NY,NX),NY,NX)**4
      RFLXG=RFLX0-THRMA
C
C     AERODYNAMIC RESISTANCE ABOVE SOIL SURFACE INCLUDING
C     RESISTANCE IMPOSED BY PLANT CANOPY
C
C     THETPX*=air-filled porosity of soil
C     DFVR=porosity limitation to diffusion through soil
C     POROQ=soil tortuosity
C     RAR1=porosity-limited litter blr
C     RAGZ=combined soil+litter blr
C     RI=Richardsons number
C     RIB=isothermal RI
C     TKQ,TK1=canopy air,soil temperature
C     RAGZ,RA=soil+litter blr
C     RAGS=isothermal blr at ground surface
C
      THETPX0=AMAX1(ZERO,THETPX(0,NY,NX))
      DFVR=THETPX0*POROQ*THETPX0/POROS(0,NY,NX)
      RAR1=RAG(NY,NX)+RAR(NY,NX)/DFVR
      RI=AMAX1(-0.3,AMIN1(0.075,RIB(NY,NX)*(TKQ(NY,NX)
     2-TK1(NUM(NY,NX),NY,NX))))
      RAGX=AMAX1(RAM,0.8*RAGS(NY,NX),AMIN1(1.2*RAGS(NY,NX)
     2,RAR1/(1.0-10.0*RI)))
      RAGS(NY,NX)=RAGX
      RA=RAGR(NY,NX)+RAGS(NY,NX)
C     IF(I.EQ.63.AND.NX.EQ.1)THEN
C     WRITE(*,7776)'RAGX',I,J,M,NX,NY,RAGZ,BARE(NY,NX),RAG(NY,NX)
C    2,CVRDW(NY,NX),RAR1,RI,RIB(NY,NX),TKQ(NY,NX),TK1(NUM(NY,NX),NY,NX)
C    3,TK1(0,NY,NX),RAGX,RAM,RAGS(NY,NX),RA
C    4,RAR(NY,NX),DFVR,THETPX0,POROQ,THETPX(0,NY,NX)
C    5,DLYRR(NY,NX),WGSGR(NY,NX)
7776  FORMAT(A8,5I4,30E12.4)
C     ENDIF
C
C     PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
C
C     PARE,PARS=blcs for latent,sensible heat fluxes over soil
C     PAREG,PARSG=conductances for latent,sensible heat fluxes
C     RZ=minimum surface resistance
C     VP1,VPQ=vapor pressure at soil surface, canopy air
C     EVAPG=evaporation
C     EFLXG=latent heat flux
C     XH=rate constant
C     VOLW2=soil water volume
C     VAP=latent heat of evaporation
C     VFLXG=convective heat of evaporation flux
C
      PARE=PAREG(NY,NX)/(RA+RZ)
      PARS=PARSG(NY,NX)/RA
      TKX1=TK1(NUM(NY,NX),NY,NX)
      VP1=2.173E-03/TKX1
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKX1))
     3*EXP(18.0*PSISV1/(8.3143*TKX1))
      EVAPG(NY,NX)=AMAX1(PARE*(VPQ(NY,NX)-VP1)
     2,-AMAX1(0.0,VOLW2(NUM(NY,NX),NY,NX)*XNPX))
      EFLXG=EVAPG(NY,NX)*VAP
      IF(EVAPG(NY,NX).LT.0.0)THEN
      VFLXG=EVAPG(NY,NX)*4.19*TK1(NUM(NY,NX),NY,NX)
      ELSE
      VFLXG=EVAPG(NY,NX)*4.19*TKQ(NY,NX)
      ENDIF
      VOLW2(NUM(NY,NX),NY,NX)=VOLW2(NUM(NY,NX),NY,NX)+EVAPG(NY,NX)
C
C     SOLVE FOR SOIL SURFACE TEMPERATURE AT WHICH ENERGY
C     BALANCE OCCURS, SOLVE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
C
C     SFLXG,EFLXG,RFLXG=sensible,latent heat fluxes, net radiation
C     VFLXG=convective heat flux from EFLXG
C     HFLXG=storage heat flux
C
      SFLXG=PARS*(TKQ(NY,NX)-TK1(NUM(NY,NX),NY,NX))
      HFLX0=RFLXG+EFLXG+SFLXG
      HFLXG=HFLX0+VFLXG
C     IF(I.GT.331)THEN
C     WRITE(*,1112)'EFLXG',I,J,M,NX,NY,NUM(NY,NX),EVAPG(NY,NX)
C    2,TK1(NUM(NY,NX),NY,NX),TK1(0,NY,NX),DLYRR(NY,NX)
C    2,RFLXG,EFLXG,SFLXG,VFLXG,HFLX0,HFLXG,RA,RAC(NY,NX),RAG(NY,NX)
C    3,RAGZ,RAR1,RAGX,RI,RAGS(NY,NX),BKVL(NUM(NY,NX),NY,NX)
C    3,VOLW2(NUM(NY,NX),NY,NX),VOLI1(NUM(NY,NX),NY,NX),RFLX0,ALBG
C    4,RADXG(NY,NX),THRYG(NY,NX),THRMA,THRYW(NY,NX),THS(NY,NX)
C    5,VPQ(NY,NX),VP1,FRADG(NY,NX),THRMCX(NY,NX)
C    6,PSISO(NUM(NY,NX),NY,NX)
C    6,FLQM,PARE,PARS,PARSG(NY,NX),HWFLQM
C    7,THETPY(NUM(NY,NX),NY,NX),RAR(NY,NX),THETPY(0,NY,NX)
C    8,VHCP1(0,NY,NX),VHCPRX(NY,NX),VHCP1(NUM(NY,NX),NY,NX)
C    3,TKQ(NY,NX),BARE(NY,NX),THETWX(NUM(NY,NX),NY,NX)
C    5,PSISM1(NUM(NY,NX),NY,NX),THETW1,VOLW1(NUM(NY,NX),NY,NX)
1112  FORMAT(A8,6I4,60E12.4)
C     ENDIF
C
C     ENERGY BALANCE AT RESIDUE SURFACE
C
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
C
C     PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
C
      EVAPR(NY,NX)=0.0
      RFLXR=0.0
      EFLXR=0.0
      VFLXR=0.0
      SFLXR=0.0
      HFLXR=0.0
      FLV1=0.0
      HWFLV1=0.0
      HFLCR1=0.0
      THRMZ=0.0
C
C     NET RADIATION AT RESIDUE SURFACE
C
C     ALBR=litter albedo
C     BKVL=litter mass
C     VOLW1,VOLI1=water,ice volume in litter
C     RADXR,THRYR=incoming shortwave,longwave radiation
C     TKR1,TKS1=litter,soil temperature
C     VOLWR2,VOLW12=litter,soil water volume
C     VHCPR2,VHCP12=litter,soil heat capacity
C
C
      ALBR=(0.20*BKVL(0,NY,NX)+0.06*VOLW1(0,NY,NX)+0.30
     2*VOLI1(0,NY,NX))/(BKVL(0,NY,NX)+VOLW1(0,NY,NX)+VOLI1(0,NY,NX))
      RFLX0=(1.0-ALBR)*RADXR(NY,NX)+THRYR(NY,NX)
      TKR1=TK1(0,NY,NX)
      VOLWR2=VOLW1(0,NY,NX)
      VHCPR2=VHCP1(0,NY,NX)
      TKS1=TK1(NUM(NY,NX),NY,NX)
      VOLW12=VOLW1(NUM(NY,NX),NY,NX)
      VHCP12=VHCP1(NUM(NY,NX),NY,NX)
C
C     THERMAL CONDUCTIVITY BETWEEN SURFACE RESIDUE AND SOIL SURFACE
C
C     CNVR,CNV1=litter,soil vapor conductivity
C     THETPM=litter air concentration
C     POROS,POROQ=litter porosity, tortuosity
C     WGSGR,WGSGL=litter,soil vapor diffusivity
C     CVRD=litter cover fraction
C     ATCNVR=litter-soil vapor conductance
C     DLYRR,DLYR=litter,soil depths
C     THETRR=dry litter concentration
C     DTH*,RYL*,DNU*,TRB*=turbulence effects on thermal conductivity
C     WTHET0,WTHET1=multiplier for air concn in thermal conductivity
C     TCNDW*,TCNDA*=thermal conductivity of water,air
C     TCNDR,TCND1=litter,soil thermal conductivity
C     ATCNDR=litter-soil thermal conductance
C
      CNVR=WGSGR(NY,NX)*THETPM(M,0,NY,NX)*POROQ
     2*THETPM(M,0,NY,NX)/POROS(0,NY,NX)
      CNV1=WGSGL(NUM(NY,NX),NY,NX)*THETPM(M,NUM(NY,NX),NY,NX)*POROQ
     2*THETPM(M,NUM(NY,NX),NY,NX)/POROS(NUM(NY,NX),NY,NX)
      IF(CVRD(NY,NX).GT.ZERO)THEN
      IF(CNVR.GT.ZERO.AND.CNV1.GT.ZERO)THEN
      ATCNVR=2.0*CNVR*CNV1
     2/(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+CNV1*DLYRR(NY,NX))
      ELSE
      ATCNVR=2.0*CNVR
     2/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR(NY,NX))*CVRD(NY,NX)
      ENDIF
      ELSE
      ATCNVR=0.0
      ENDIF
      THETRR=AMAX1(0.0,1.0-THETPX(0,NY,NX)-THETWX(0,NY,NX)
     2-THETIX(0,NY,NX))
      DTKX=ABS(TK1(0,NY,NX)-TK1(NUM(NY,NX),NY,NX))*1.0E-06
      DTHW0=AMAX1(0.0,THETWX(0,NY,NX)-TRBW)**3
      DTHA0=AMAX1(0.0,THETPX(0,NY,NX)-TRBA)**3
      DTHW1=AMAX1(0.0,THETWX(NUM(NY,NX),NY,NX)-TRBW)**3
      DTHA1=AMAX1(0.0,THETPX(NUM(NY,NX),NY,NX)-TRBA)**3
      RYLXW0=DTKX*DTHW0
      RYLXA0=DTKX*DTHA0
      RYLXW1=DTKX*DTHW1
      RYLXA1=DTKX*DTHA1
      RYLNW0=AMIN1(1.0E+04,RYLXW*RYLXW0)
      RYLNA0=AMIN1(1.0E+04,RYLXA*RYLXA0)
      RYLNW1=AMIN1(1.0E+04,RYLXW*RYLXW1)
      RYLNA1=AMIN1(1.0E+04,RYLXA*RYLXA1)
      XNUSW0=AMAX1(1.0,0.68+0.67*RYLNW0**0.25/DNUSW)
      XNUSA0=AMAX1(1.0,0.68+0.67*RYLNA0**0.25/DNUSA)
      XNUSW1=AMAX1(1.0,0.68+0.67*RYLNW1**0.25/DNUSW)
      XNUSA1=AMAX1(1.0,0.68+0.67*RYLNA1**0.25/DNUSA)
      TCNDW0=2.067E-03*XNUSW0
      TCNDA0=9.050E-05*XNUSA0
      TCNDW1=2.067E-03*XNUSW1
      TCNDA1=9.050E-05*XNUSA1
      WTHET0=1.467-0.467*THETPY(0,NY,NX)
      WTHET1=1.467-0.467*THETPY(NUM(NY,NX),NY,NX)
      TCNDR=(0.779*THETRR*9.050E-04+0.622*THETWX(0,NY,NX)*TCNDW0
     2+0.380*THETIX(0,NY,NX)*7.844E-03
     3+WTHET0*THETPX(0,NY,NX)*TCNDA0)
     4/(0.779*THETRR+0.622*THETWX(0,NY,NX)
     5+0.380*THETIX(0,NY,NX)+WTHET0*THETPX(0,NY,NX))
      TCND1=(STC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX)*TCNDW1
     2+0.611*THETIX(NUM(NY,NX),NY,NX)*7.844E-03
     3+WTHET1*THETPX(NUM(NY,NX),NY,NX)*TCNDA1)
     4/(DTC(NUM(NY,NX),NY,NX)+THETWX(NUM(NY,NX),NY,NX)
     5+0.611*THETIX(NUM(NY,NX),NY,NX)+WTHET1*THETPX(NUM(NY,NX),NY,NX))
      ATCNDR=2.0*TCNDR*TCND1/(TCNDR*DLYR(3,NUM(NY,NX),NY,NX)
     2+TCND1*DLYRR(NY,NX))
C
C     SMALLER TIME STEP FOR SOLVING SURFACE RESIDUE ENERGY EXCHANGE
C
      DO 5000 NN=1,NPR
      IF(VHCPR2.GT.VHCPRX(NY,NX))THEN
C
C     AERODYNAMIC RESISTANCE ABOVE RESIDUE INCLUDING
C     RESISTANCE IMPOSED BY PLANT CANOPY
C
C     RI=Richardsons number
C     RIB=isothermal RI
C     TKQ,TKR1=canopy air,litter temperature
C     RZ=surface resistance to evaporation
C     RAGX,RA=litter blr
C     RAG,RAGR=isothermal blr at ground surface
C     PARE,PARS=blcs for litter latent,sensible heat fluxes
C
      RI=AMAX1(-0.3,AMIN1(0.075,RIB(NY,NX)*(TKQ(NY,NX)-TKR1)))
      RAGX=AMAX1(RAM,0.8*RAGR(NY,NX),AMIN1(1.2*RAGR(NY,NX)
     2,RARG(NY,NX)/(1.0-10.0*RI)))
      RAGR(NY,NX)=RAGX
      RA=RAGX
      PARE=PARER(NY,NX)/(RA+RZ)
      PARS=PARSR(NY,NX)/RA
C
C     NET RADIATION AT RESIDUE SURFACE
C
C     THRMZ2=longwave radiation emitted by litter
C     RFLXR2=litter net radiation
C     THETWR=litter water content
C     VOLWRX=maximum water retention by litter
C     PSISM1=litter matric water potential
C
      THRMZ2=THRMR(NY,NX)*TKR1**4
      RFLXR2=RFLX0-THRMZ2
      IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLWR2)/VOLR(NY,NX)
      ELSE
      THETWR=POROS0(NY,NX)
      ENDIF
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX)
     2.AND.VOLW1(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
      IF(THETWR.LT.FC(0,NY,NX))THEN
      PSISM1(0,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCL(0,NY,NX)-LOG(THETWR))
     3/FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETWR.LT.POROS0(NY,NX))THEN
      PSISM1(0,NY,NX)=-EXP(PSIMS(NY,NX)
     2+(((PSL(0,NY,NX)-LOG(THETWR))
     3/PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
      ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
      ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
C
C     VAPOR FLUX AT RESIDUE SURFACE
C
C     VPR,VP1,VPQ=vapor pressure in litter,soil,canopy air
C     TKS1=soil temperature
C     EVAPR2=litter evaporation
C     EFLXR2=litter latent heat flux
C     VAP=latent heat of evaporation
C     VFLXR2=convective heat of evaporation flux
C
      VPR=2.173E-03/TKR1
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKR1))
     3*EXP(18.0*PSISM1(0,NY,NX)/(8.3143*TKR1))
      VP1=2.173E-03/TKS1
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKS1))
     3*EXP(18.0*PSISV1/(8.3143*TKS1))
      EVAPR2=AMAX1(-AMAX1(0.0,VOLWR2*XNPX)
     2,PARE*(VPQ(NY,NX)-VPR))
      EFLXR2=EVAPR2*VAP
      VFLXR2=EVAPR2*4.19*TKR1
C
C     SOLVE FOR RESIDUE TO SOIL SURFACE HEAT FLUXES
C
C     FLVC,FLVX=vapor unconstrained,vapor constrained vapor flux
C     XNPZ=time step for litter flux calculations
C     VPY=equilibrium vapor concentration
C     VOLPM=litter,soil air filled porosity
C     FLV2=litter soil vapor flux
C     HWFLV2=convective heat of litter soil vapor flux
C     TKXR,TK1X=interim calculation of litter,soil temperatures
C     TKY=equilibrium litter-soil temperature
C     HFLWC,HFLWX=litter-soil heat flux unlimited,limited by heat
C     HFLCR2=litter-soil heat flux
C
      IF(THETPM(M,0,NY,NX).GT.THETX
     2.AND.THETPM(M,NUM(NY,NX),NY,NX).GT.THETX)THEN
      FLVC=ATCNVR*(VPR-VP1)*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)
     2*CVRD(NY,NX)*XNPZ
      VPY=(VPR*VOLPM(M,0,NY,NX)+VP1*VOLPM(M,NUM(NY,NX),NY,NX))
     2/(VOLPM(M,0,NY,NX)+VOLPM(M,NUM(NY,NX),NY,NX))
      FLVX=(VPR-VPY)*VOLPM(M,0,NY,NX)*XNPB
      IF(FLVC.GE.0.0)THEN
      FLV2=AMAX1(0.0,AMIN1(FLVC,FLVX))
      HWFLV2=(4.19*TKR1+VAP)*FLV2
      ELSE
      FLV2=AMIN1(0.0,AMAX1(FLVC,FLVX))
      HWFLV2=(4.19*TKS1+VAP)*FLV2
      ENDIF
      ELSE
      FLV2=0.0
      HWFLV2=0.0
      ENDIF
      TKXR=TKR1-HWFLV2/VHCPR2
      TK1X=TKS1+HWFLV2/VHCP12
      TKY=(TKXR*VHCPR2+TK1X*VHCP12)/(VHCPR2+VHCP12)
      HFLWX=(TKXR-TKY)*VHCPR2*XNPB
      HFLWC=ATCNDR*(TKXR-TK1X)*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)
     2*CVRD(NY,NX)*XNPZ
      IF(HFLWC.GE.0.0)THEN
      HFLCR2=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLCR2=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
C
C     SOLVE FOR RESIDUE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
C
C     SFLXR2,RFLXR2,EFLXR2=litter sensible,net radn,latent heat fluxes
C     HFLX02,HFLXR2=storage,total litter heat flux
C
      SFLXR2=PARS*(TKQ(NY,NX)-TKR1)
      HFLX02=RFLXR2+EFLXR2+SFLXR2
      HFLXR2=HFLX02+VFLXR2
C
C     AGGREGATE WATER AND ENERGY FLUXES FROM TIME STEP FOR LITTER
C     CALCULATIONS TO THAT FOR SOIL PROFILE
C
      EVAPR(NY,NX)=EVAPR(NY,NX)+EVAPR2
      RFLXR=RFLXR+RFLXR2
      EFLXR=EFLXR+EFLXR2
      VFLXR=VFLXR+VFLXR2
      SFLXR=SFLXR+SFLXR2
      HFLXR=HFLXR+HFLXR2
      FLV1=FLV1+FLV2
      HWFLV1=HWFLV1+HWFLV2
      HFLCR1=HFLCR1+HFLCR2
      THRMZ=THRMZ+THRMZ2
      ELSE
      EVAPR2=0.0
      RFLXR2=0.0
      EFLXR2=0.0
      VFLXR2=0.0
      SFLXR2=0.0
      HFLXR2=0.0
      FLV2=0.0
      HWFLV2=0.0
      HFLCR2=0.0
      THRMZ2=0.0
      ENDIF
      VOLWR2=VOLWR2+FLYM2+EVAPR2-FLV2
      VOLW12=VOLW12+FLV2
      ENGYR=VHCPR2*TKR1
      TKR10=TKR1
      VHCPRXX=VHCPR2
      VHCPR2=2.496E-06*ORGC(0,NY,NX)+4.19*VOLWR2
     2+1.9274*VOLI1(0,NY,NX)
      VHCP12=VHCP12+4.19*FLV2
      TKR1=(ENGYR+HFLXR2+HWFLM2-HWFLV2-HFLCR2)/VHCPR2
C      if(I==162)then
C      print*,'tkrcp',VHCPRXX,VHCPR2,VHCPRXX/VHCPR2
C      print*,'tkr1',j,TKR1,TKR10,TK1(0,NY,NX)
C     2,TK1(NUM(NY,NX),NY,NX),TKS(0,NY,NX)
C     2,HFLXR2+HWFLM2-HWFLV2-HFLCR2,VHCPR2
C      print*,'tkr11',j,HFLXR2,HWFLM2,HWFLV2,HFLCR2
C      endif
!     IF(ABS(VHCPRXX/VHCPR2-1.)>0.025.or.
!    2abs(TKR1/TKR10-1.)>0.025)then
!     TKR1=TK1(0,NY,NX)
C      print*,'vhcpx',VHCPRXX,VHCPR2,TKR1,TKR10
C      print*,'tkr1',j,TKA(NY,NX),TK1(0,NY,NX)
C     2,TK1(NUM(NY,NX),NY,NX),TKS(0,NY,NX)
C      pause
!     endif
      TKS1=TKS1+(HWFLV2+HFLCR2)/VHCP12
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,1111)'EFLXR2',I,J,M,NX,NY,NUM(NY,NX),NN
C    2,TKR1,TKS1,TKQ(NY,NX),RFLXR2,EFLXR2,SFLXR2,VFLXR2
C    3,HFLX02,HFLXR2,HWFLM2,HWFLV2,HFLCR2,VHCPR2,HFLXR
C    3,RA,RAGX,RAG(NY,NX),RAB(NY,NX),RAC(NY,NX)
C    4,RAR1,PARE,HFLWX,HFLWC,TKXR,TK1X,TKY
C    2,VPQ(NY,NX),VPR,EVAPR(NY,NX),EVAPR2,VOLWR2,XNPX
C    3,HFLCR2,HWFLV2,VHCPR2,VHCP12,RFLX0,THRMZ2,VHCPRX(NY,NX)
C    4,ALBR,RADXR(NY,NX),THRYR(NY,NX),THRMR(NY,NX),XNPZ
C    3,FLV1,FLV2,VPR,VP1,CNVR,CNV1,FLVC,FLVX,XNPZ,XNPR
C    3,PSISM1(0,NY,NX),PSISV1,THETWR,VOLWRX(NY,NX),ORGC(0,NY,NX)
C    4,VHCPRX(NY,NX),PARS,PARE,RA,RZ,RI,TKQ(NY,NX),VOLW1(0,NY,NX)
C    5,VOLW1(NUM(NY,NX),NY,NX),VOLT(NUM(NY,NX),NY,NX),FLV1
C    5,CNVR,CNV1,VOLX(0,NY,NX),POROQ,WGSGR(NY,NX)
C    5,ATCNDR,TCNDR
C    6,TCND1,STC(NUM(NY,NX),NY,NX),THETWX(NUM(NY,NX),NY,NX)
C    2,THETIX(NUM(NY,NX),NY,NX),WTHET1,THETPX(NUM(NY,NX),NY,NX),TCNDA1
C    4,DTC(NUM(NY,NX),NY,NX),VOLP1(0,NY,NX),VOLR(NY,NX)
C    7,THETWX(0,NY,NX),THETIX(0,NY,NX),THETPY(0,NY,NX),ORGC(0,NY,NX)
C    6,DLYR(3,0,NY,NX),DLYR(3,NUM(NY,NX),NY,NX)
C    8,CVRD(NY,NX),XVOLW(NY,NX),VOLWG(NY,NX)
1111  FORMAT(A8,7I4,100E12.4)
C     ENDIF
5000  CONTINUE
C
C     IF NO SURFACE RESIDUE
C
      ELSE
      TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
      EVAPR(NY,NX)=0.0
      RFLXR=0.0
      EFLXR=0.0
      VFLXR=0.0
      SFLXR=0.0
      HFLXR=0.0
      FLV1=0.0
      HWFLV1=0.0
      HFLCR1=0.0
      THRMZ=0.0
      ENDIF
C
C     GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
C     FOR LATER UPDATES TO STATE VARIABLES
C
C     FLWLG,FLWHLG=water flux from atm to soil micropores,macropores
C     HFLWLG=convective heat flux from atm to soil
C     FLWRLG=water flux from atm to litter
C     HFLWRLG=convective heat flux from atm to litter
C     FLWVLS=water flux within soil accounting for wetting front
C
      FLWLG=FLQM+EVAPG(NY,NX)+FLV1
C      if(I==162)then
C      print*,'FLWLG=',j,FLQM,EVAPG(NY,NX),FLV1
C      endif
      FLWLXG=FLQM+EVAPG(NY,NX)+FLV1
      FLWHLG=FLHM
      HFLWLG=HWFLQM+HFLXG+HWFLV1+HFLCR1
      FLWRLG=FLYM+EVAPR(NY,NX)-FLV1
      HFLWRLG=HWFLYM+HFLXR-HWFLV1-HFLCR1
      FLWVLS=(VOLW1(NUM(NY,NX),NY,NX)-VOLWX1(NUM(NY,NX),NY,NX))*XNPH
C     IF(IYRC.EQ.1813.AND.I.EQ.231)THEN
C     WRITE(*,7754)'FLWLG',I,J,M,NX,NY
C    2,FLWLG,FLQM,EVAPG(NY,NX),FLV1
C    3,FLWRLG,FLYM,EVAPR(NY,NX),FLV1
C    4,HFLWLG,HWFLQM,HFLXG,HWFLV1,HFLCR1
C    5,HFLWRLG,HWFLYM,HFLXR
C     ENDIF
C
C     GENERATE NEW SNOWPACK
C
C     XFLWS,XFLWW,XFLWI=hourly snow,water,ice transfer
C     FLQ0S,FLQ0W,FLQ0I=snow,water,ice input to snowpack
C     XHFLWW=hourly convective heat flux from snow,water,ice transfer
C     HWFLQ0=convective heat flux from snow,water,ice to snowpack
C
      IF(VHCPW(1,NY,NX).LE.VHCPWX(NY,NX)
     2.AND.FLQ0S(NY,NX).GT.ZEROS(NY,NX))THEN
      XFLWS(1,NY,NX)=XFLWS(1,NY,NX)+FLQ0S(NY,NX)
      XFLWW(1,NY,NX)=XFLWW(1,NY,NX)+FLQ0W(NY,NX)
      XFLWI(1,NY,NX)=XFLWI(1,NY,NX)+FLQ0I(NY,NX)
      XHFLWW(1,NY,NX)=XHFLWW(1,NY,NX)+HWFLQ0(NY,NX)
C     WRITE(*,4422)'INIT',I,J,FLQ0S(NY,NX),FLQ0W(NY,NX)
C    3,FLQ0I(NY,NX),HWFLQ0(NY,NX),XFLWS(1,NY,NX),XFLWW(1,NY,NX)
C    2,XFLWI(1,NY,NX),XHFLWW(1,NY,NX),HFLWL(3,NUM(NY,NX),NY,NX)
C    3,HFLW(3,NUM(NY,NX),NY,NX),FSNX(NY,NX),VHCP1(NUM(NY,NX),NY,NX)
C    4*TK1(NUM(NY,NX),NY,NX),HFLWRL(NY,NX),HFLWR(NY,NX)
C    5,VHCP1(0,NY,NX)*TK1(0,NY,NX),HEATH(NY,NX),RFLXG,RFLXR,RFLXW
C    2,SFLXG,SFLXR,SFLXW,EFLXG,EFLXR,EFLXW,VFLXG,VFLXR,VFLXW
      ENDIF
      THRMG(NY,NX)=THRMG(NY,NX)+THRMA+THRMZ
      ELSE
      RFLXG=0.0
      EFLXG=0.0
      VFLXG=0.0
      SFLXG=0.0
      HFLXG=0.0
      RFLXR=0.0
      EFLXR=0.0
      VFLXR=0.0
      SFLXR=0.0
      HFLXR=0.0
      FLWLG=0.0
      FLWLXG=0.0
      FLWHLG=0.0
      HFLWLG=0.0
      FLWRLG=0.0
      HFLWRLG=0.0
      FLWVLS=0.0
      EVAPG(NY,NX)=0.0
      EVAPR(NY,NX)=0.0
      ENDIF
C
C     AGGREGATE RESIDUE AND SOIL SURFACE FLUXES BENEATH SNOW
C     AND ATMOSPHERE
C
C     FLWL,FLWLX=total water flux into soil micropores
C     FLWHL=total water flux into soil macropores
C     HFLWL=total heat flux into soil
C     FLWRL,FLWLX=total water flux into litter
C     HFLWRL=total heat flux into litter
C     FLWV*=total internal vapor flux in soil
C
      FLWL(3,NUM(NY,NX),NY,NX)=FLWLW+FLWLG
C      if(I==162)then
C      print*,'flwljj',j,M,FLWLW,FLWLG,FLWLW+FLWLG
C      endif
      FLWLX(3,NUM(NY,NX),NY,NX)=FLWLXW+FLWLXG
      FLWHL(3,NUM(NY,NX),NY,NX)=FLWHLW+FLWHLG
      HFLWL(3,NUM(NY,NX),NY,NX)=HFLWLW+HFLWLG
      FLWRL(NY,NX)=FLWRLW+FLWRLG
      HFLWRL(NY,NX)=HFLWRLW+HFLWRLG
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,7756)'FLWT',I,J,M,NX,NY,NUM(NY,NX)
C    2,FLWL(3,NUM(NY,NX),NY,NX),FLWLW,FLWLG
C    4,FLWHL(3,NUM(NY,NX),NY,NX),FLWHLW,FLWHLG
C    5,HFLWL(3,NUM(NY,NX),NY,NX),HFLWLW,HFLWLG
C    6,FLWRL(NY,NX),FLWRLW,FLWRLG
C    7,HFLWRL(NY,NX),HFLWRLW,HFLWRLG
C    9,HFLWQG,HFLVS1,HFLWS1,HFLVR1,HFLWR1
C    1,HFLWQR,HFLVSR,HFLWSR,HFLVR1,HFLWR1
C    8,FLQ0S(NY,NX),FLQ0W(NY,NX),FLQ0I(NY,NX),FLQM,FLHM,FLYM
C    9,PRECA(NY,NX)*XNPH,PRECW(NY,NX)*XNPH
C    7,RFLXW,EFLXW,SFLXW,-HFLXW+VFLXW,RFLXG,EFLXG,SFLXG,-HFLXG+VFLXG
C    7,RFLXR,EFLXR,SFLXR,-HFLXR+VFLXR
C    8,FLW0S(1,NY,NX),FLQ0S(NY,NX),EVAPS(NY,NX)
C    9,FLW0W(1,NY,NX),FLQ0W(NY,NX),EVAPW(NY,NX)
C    1,FLW0I(1,NY,NX),FLQ0I(NY,NX)
C    2,HFLW0W(1,NY,NX),HWFLQ0(NY,NX),HFLXW
C    3,XFLWS(1,NY,NX),XFLWW(1,NY,NX),XFLWI(1,NY,NX),XHFLWW(1,NY,NX)
C    8,FSNW(NY,NX),FSNX(NY,NX),CVRD(NY,NX),BARE(NY,NX)
C    9,THETPM(M,NUM(NY,NX),NY,NX)
7756  FORMAT(A8,6I4,60E12.4)
C     ENDIF
C
C     CAPILLARY EXCHANGE OF WATER BETWEEN SOIL SURFACE AND RESIDUE
C
C     CNDR,HCNDR=current,saturated litter hydraulic conductivity
C     PSISE,PSISM1(0,=air entry,current litter water potential
C     VOLW1(0,VOLWRX=current,maximum litter water volume
C     CND1,HCND=soil hydraulic conductivity
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
C     K1=soil relative water-filled porosity
C     THETWX,POROS=soil water content,porosity
C     AVCNDR=litter-soil hydraulic conductance
C     DLYRR,DLYR=litter,soil thicknesses
C     FLQR,FLQL=litter-soil water flux unltd,ltd by water
C     XNPH=time step of flux calculations
C     CVRD=fraction of litter cover
C     PSISM1(NUM=soil matric water potential
C     VOLW1(NUM=soil water volume
C     HFLQL=convective heat from litter-soil water flux
C     FLWL,HFLWL=micropore water,heat flux
C     FLWRL,HFLWRL=total litter water,heat flux
C     FLWRM=litter-soil water flux for solute transfer in trnsfr.f
C     CND1,CNDL=hydraulic conductivity of source,destination layer
C     HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity
C
      IF(BKDS(NUM(NY,NX),NY,NX).GT.ZERO)THEN
C  VOLWRX: surface layer water holding capacity
C  THETWR: surface layer saturation
      IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
      ELSE
      THETWR=POROS0(NY,NX)
      ENDIF

      THETW1=AMAX1(THETY(NUM(NY,NX),NY,NX)
     2,AMIN1(POROS(NUM(NY,NX),NY,NX)
     2,VOLW1(NUM(NY,NX),NY,NX)/VOLY(NUM(NY,NX),NY,NX)))
      K0=MAX(1,MIN(100,INT(100.0*(AMAX1(0.0,POROS0(NY,NX)
     2-THETWR))/POROS0(NY,NX))+1))
      K1=MAX(1,MIN(100,INT(100.0*(AMAX1(0.0,POROS(NUM(NY,NX),NY,NX)
     2-THETW1))/POROS(NUM(NY,NX),NY,NX))+1))
      CNDR=HCND(3,K0,0,NY,NX)
      CND1=HCND(3,K1,NUM(NY,NX),NY,NX)*FKSAT
      AVCNDR=2.0*CNDR*CND1/(CNDR*DLYR(3,NUM(NY,NX),NY,NX)
     2+CND1*DLYRR(NY,NX))
      PSIST0=PSISM1(0,NY,NX)+PSISH(0,NY,NX)
     2+PSISO(0,NY,NX)
      PSIST1=PSISM1(NUM(NY,NX),NY,NX)+PSISH(NUM(NY,NX),NY,NX)
     2+PSISO(NUM(NY,NX),NY,NX)
C    psist0 : pressure in residual layer layer
C    psist1 : pressure in topsoil layer
C    flqx : water flux from layer 0 into layer 1 (>0)
      FLQX=AVCNDR*(PSIST0-PSIST1)
     2*AREA(3,NUM(NY,NX),NY,NX)*CVRDW(NY,NX)*XNPH

      IF(FLQX.GE.0.0)THEN
C layer 0 into layer 1
C THETW1: saturation of layer 1
      IF(THETWR.GT.THETS(0,NY,NX))THEN
      FLQZ=FLQX+AMIN1((THETWR-THETS(0,NY,NX))
     2*VOLR(NY,NX),AMAX1(0.0,(THETS(NUM(NY,NX),NY,NX)-THETW1)
     3*VOLY(NUM(NY,NX),NY,NX)))*XNPX
      ELSE
      FLQZ=FLQX
      ENDIF
C positve
      FLQR=AMAX1(0.0,AMIN1(FLQZ,VOLW1(0,NY,NX)*XNPX
     2,VOLP1(NUM(NY,NX),NY,NX)))
      FLQ2=AMAX1(0.0,AMIN1(FLQX,VOLW1(0,NY,NX)*XNPX
     2,VOLP1(NUM(NY,NX),NY,NX)))

      ELSE
C layer 1 into layer 0, FLQX<0
      IF(THETW1.GT.THETS(NUM(NY,NX),NY,NX))THEN
C more water than saturation in layer 1,remove excessive water
      FLQZ=FLQX+AMAX1((THETS(NUM(NY,NX),NY,NX)-THETW1)
     2*VOLY(NUM(NY,NX),NY,NX),AMIN1(0.0,(THETWR-THETS(0,NY,NX))
     3*VOLR(NY,NX)))*XNPX
      ELSE
      FLQZ=FLQX
      ENDIF
      FLQR=AMIN1(0.0,AMAX1(FLQZ,-VOLW1(NUM(NY,NX),NY,NX)*XNPX
     2,-VOLP1(0,NY,NX)))
      FLQ2=AMIN1(0.0,AMAX1(FLQX,-VOLW1(NUM(NY,NX),NY,NX)*XNPX
     2,-VOLP1(0,NY,NX)))
      ENDIF

      IF(VOLP1Z(NUM(NY,NX),NY,NX).LT.0.0)THEN
      FLQR=FLQR+AMIN1(0.0,AMAX1(-VOLW1(NUM(NY,NX),NY,NX)*XNPX
     2,VOLP1Z(NUM(NY,NX),NY,NX)))
      FLQ2=FLQ2+AMIN1(0.0,AMAX1(-VOLW1(NUM(NY,NX),NY,NX)*XNPX
     2,VOLP1Z(NUM(NY,NX),NY,NX)))
      ENDIF

      IF(FLQR.GT.0.0)THEN
C from layer 0 to layer 1
      HFLQR=4.19*TK1(0,NY,NX)*FLQR
      ELSE
C from layer 1 to layer 0
      HFLQR=4.19*TK1(NUM(NY,NX),NY,NX)*FLQR
      ENDIF

      FLWL(3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)+FLQR
      HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQR
      FLWRL(NY,NX)=FLWRL(NY,NX)-FLQR
      HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQR
      FLWRM(M,NY,NX)=FLQR
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,4322)'FLQR',I,J,M,NX,NY,NUM(NY,NX),K0,K1
C    2,FLQR,FLWRL(NY,NX),FLWL(3,NUM(NY,NX),NY,NX),FLQX,FLQZ,FLQ2
C    2,PSISM1(0,NY,NX),PSISM1(NUM(NY,NX),NY,NX),PSIST0,PSIST1
C    3,CVRDW(NY,NX),CNDR,CND1,AVCNDR,FKSAT
C    3,POROS0(NY,NX),VOLW1(0,NY,NX)
C    2,VOLP1(NUM(NY,NX),NY,NX),FLWL(3,NUM(NY,NX),NY,NX)
C    3,VOLP1Z(NUM(NY,NX),NY,NX),THETWX(0,NY,NX)
C    4,THETWX(NUM(NY,NX),NY,NX)
C    4,THETS(0,NY,NX),THETS(NUM(NY,NX),NY,NX)
C    4,THETWR,THETW1,XVOLT(NY,NX)
C    3,HFLQR,HFLWRL(NY,NX),HFLWL(3,NUM(NY,NX),NY,NX)
C    6,THETWR,VHCP1(0,NY,NX),VHCPRX(NY,NX)
C    2,FLWLYH,VOLX(NUM(NY,NX),NY,NX),VOLA1(NUM(NY,NX),NY,NX)
C    3,VOLP1(NUM(NY,NX),NY,NX),VOLW1(NUM(NY,NX),NY,NX)
C    3,VOLI1(NUM(NY,NX),NY,NX),VOLP1(0,NY,NX),VOLW1(0,NY,NX)
C    3,VOLI1(0,NY,NX),PSISM1(0,NY,NX),VOLP1Z(NUM(NY,NX),NY,NX)
C    4,VOLX(NUM(NY,NX),NY,NX),VOLAH(NUM(NY,NX),NY,NX)
C    4,PSISM1(NUM(NY,NX),NY,NX),AVCNDR
C    2,VOLAH1(NUM(NY,NX),NY,NX),VOLPH1(NUM(NY,NX),NY,NX)
C    2,VOLWH1(NUM(NY,NX),NY,NX),VOLIH1(NUM(NY,NX),NY,NX)
4322  FORMAT(A8,8I4,40E12.4)
C     ENDIF
      ELSE
      FLQR=XVOLW(NY,NX)*XNPX
      HFLQR=4.19*TK1(0,NY,NX)*FLQR
      FLWL(3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)+FLQR
      HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQR
      FLWRL(NY,NX)=FLWRL(NY,NX)-FLQR
      HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQR
      FLWRM(M,NY,NX)=FLQR
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,4323)'FLQ0',I,J,M,NX,NY
C    2,FLQR,XVOLW(NY,NX),XNPX,FLWRL(NY,NX)
4323  FORMAT(A8,5I4,12E12.4)
C     ENDIF
      ENDIF
C     IF((I/10)*10.EQ.I)THEN
C     WRITE(*,4321)'HCNDR',I,J,M,NX,NY,K1,AVCND1,CNDR,CND1,DLYRR(NY,NX)
C    2,PSISM1(0,NY,NX),PSISM1(NUM(NY,NX),NY,NX),FLQL,HFLQL
C    3,VOLWR2,EVAPR(NY,NX),VOLWRX(NY,NX)-VOLW1(0,NY,NX)
C    2-VOLI1(0,NY,NX),VOLW1(NUM(NY,NX),NY,NX),VOLW1(0,NY,NX)
C    4,VOLP1(NUM(NY,NX),NY,NX),POROS(NUM(NY,NX),NY,NX)
C    5,VOLWG(NY,NX),FLYM,HCNDR(NY,NX),PSISE(0,NY,NX),PSISM1(0,NY,NX)
C    6,THETWR,VHCP1(0,NY,NX),VHCPRX(NY,NX)
C    7,VOLPH1(NUM(NY,NX),NY,NX),VOLW1(0,NY,NX),VOLWRX(NY,NX)
4321  FORMAT(A8,6I4,30E12.4)
C     ENDIF
C
C     OVERLAND FLOW INTO SOIL MACROPORES WHEN WATER STORAGE CAPACITY
C     OF THE LITTER IS EXCEEDED
C
C     VOLPH1=air-filled macroporosity
C     FINHR,HFINHR=water,convective heat from litter to macropores
C     VOLW1(0,VOLWRX=current,maximum litter water volume
C     FLWL,HFLWL=micropore water,heat flux
C     FLWRL,HFLWRL=total litter water,heat flux
C
      IF(VOLPH1(NUM(NY,NX),NY,NX).GT.0.0
     2.AND.XVOLW(NY,NX).GT.0.0)THEN
      FLQHR=AMIN1(XVOLW(NY,NX)*XNPX,VOLPH1(NUM(NY,NX),NY,NX))
      HFLQHR=FLQHR*4.19*TK1(0,NY,NX)
      FLWHL(3,NUM(NY,NX),NY,NX)=FLWHL(3,NUM(NY,NX),NY,NX)+FLQHR
      HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQHR
      FLWRL(NY,NX)=FLWRL(NY,NX)-FLQHR
      HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQHR
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,4357)'FLQHR',I,J,M,NX,NY,NUM(NY,NX),FLQHR,FLWRL(NY,NX)
C    2,FLWHL(3,NUM(NY,NX),NY,NX),VOLPH1(NUM(NY,NX),NY,NX)
C    3,XVOLW(NY,NX),VOLW1(0,NY,NX),VOLWRX(NY,NX)
C    4,HFLQHR,HFLWRL(NY,NX),HFLWL(3,NUM(NY,NX),NY,NX),TK1(0,NY,NX)
4357  FORMAT(A8,6I4,40E12.4)
C     ENDIF
      ENDIF
C
C     FREEZE-THAW IN RESIDUE SURFACE FROM NET CHANGE IN RESIDUE
C     SURFACE HEAT STORAGE
C
C     TFREEZ=litter freezing temperature
C     PSISM1=litter water potential
C     VOLW1*,VOLI1=litter water,ice volume
C     VHCP1*=litter volumetric heat capacity
C     TK1*=litter temperature
C     ORGC=litter organic C
C     HFLWRL=total litter conductive, convective heat flux
C     TFLX1,TFLX=unltd,ltd latent heat from freeze-thaw
C     TFLXR,WFLXR=litter water,latent heat flux from freeze-thaw
C
      TFREEZ=-9.0959E+04/(PSISM1(0,NY,NX)-333.0)
      VOLW1X=AMAX1(0.0,VOLW1(0,NY,NX)+FLWRL(NY,NX))
      ENGYR=VHCP1(0,NY,NX)*TK1(0,NY,NX)
      VHCP1X=2.496E-06*ORGC(0,NY,NX)+4.19*VOLW1X
     2+1.9274*VOLI1(0,NY,NX)
      IF(VHCP1X.GT.ZEROS(NY,NX))THEN
      TK1X=(ENGYR+HFLWRL(NY,NX))/VHCP1X
      ELSE
      TK1X=TK1(0,NY,NX)
      ENDIF
      IF((TK1X.LT.TFREEZ
     2.AND.VOLW1(0,NY,NX).GT.ZERO*VOLT(0,NY,NX))
     3.OR.(TK1X.GT.TFREEZ
     4.AND.VOLI1(0,NY,NX).GT.ZERO*VOLT(0,NY,NX)))THEN
      TFLX1=VHCP1(0,NY,NX)*(TFREEZ-TK1X)
     2/((1.0+TFREEZ*6.2913E-03)*(1.0-0.10*PSISM1(0,NY,NX)))*XNPX
      IF(TFLX1.LT.0.0)THEN
      TFLX=AMAX1(-333.0*DENSI*VOLI1(0,NY,NX)*XNPX,TFLX1)
      ELSE
      TFLX=AMIN1(333.0*VOLW1X*XNPX,TFLX1)
      ENDIF
      TFLXR(NY,NX)=TFLX
      WFLXR(NY,NX)=-TFLX/333.0
C     IF(NX.EQ.1.AND.NY.EQ.1.AND.TFLX1.NE.0.0)THEN
C     WRITE(*,5352)'TFLXR',I,J,M,NX,NY,WFLXR(NY,NX),TFLXR(NY,NX)
C    2,TFLX1,TFREEZ,TK1(0,NY,NX),TK1X
C    2,TFREEZ-TK1X,FLWRL(NY,NX),HFLWRL(NY,NX)
C    2,THETWR,TFLX,VOLI1(0,NY,NX),VOLW1(0,NY,NX),VHCP1(0,NY,NX)
C    3,PSISM1(0,NY,NX),THAWR(NY,NX),HTHAWR(NY,NX)
5352  FORMAT(A8,5I4,30E12.4)
C     ENDIF
      ELSE
      WFLXR(NY,NX)=0.0
      TFLXR(NY,NX)=0.0
      ENDIF
C
C     THICKNESS OF WATER FILMS IN LITTER AND SOIL SURFACE
C     FROM WATER POTENTIALS FOR GAS EXCHANGE IN TRNSFR.F
C
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      FILM(M,0,NY,NX)=AMAX1(1.0E-06
     2,EXP(-13.650-0.857*LOG(-PSISM1(0,NY,NX))))
      ELSE
      FILM(M,0,NY,NX)=1.0E-03
      ENDIF
      FILM(M,NUM(NY,NX),NY,NX)=AMAX1(1.0E-06
     2,EXP(-13.650-0.857*LOG(-PSISM1(NUM(NY,NX),NY,NX))))
C
C     OVERLAND FLOW WHEN WATER STORAGE CAPACITY
C     OF THE SOIL SURFACE PLUS MACROPORES IS EXCEEDED
C
      N1=NX
      N2=NY
C
C     SURFACE WATER FLUX
C
C     N2,N1=NY,NX of source grid cell
C     XVOLT,XVOLW=excess water+ice,water in source grid cell
C     VOLWG=ground surface water retention capacity
C     VX=ponded water volume above surface retention capacity
C     D,R,S,V=depth,perimeter,slope,velocity of runoff
C     DIST=distance between source,destination
C     ZM=surface roughness height for runoff
C     Q=runoff from Mannings equation
C     QRM,QRV=runoff,velocity for erosion, solute transfer
C
      IF(XVOLT(N2,N1).GT.VOLWG(N2,N1))THEN
      VX=XVOLT(N2,N1)-VOLWG(N2,N1)
      D=VX/AREA(3,0,N2,N1)
      R=D/2.828
      V=R**0.67*SQRT(SLOPE(0,N2,N1))/ZM(N2,N1)
      Q=V*D*AREA(3,NUM(N2,N1),N2,N1)*3.6E+03*XNPH
      VOLW1X=AMAX1(0.0,VOLW1(0,N2,N1)+WFLXR(N2,N1))
      QRM(M,N2,N1)=AMIN1(Q,VX*XNPX,VOLW1X*XNPX)
     2*XVOLW(N2,N1)/XVOLT(N2,N1)
      QRV(M,N2,N1)=V
C     IF(I.EQ.232)THEN
C     WRITE(*,5554)'QRINT',I,J,M,N1,N2,QRM(M,N2,N1),QRV(M,N2,N1)
C    2,Q,V,D,VX*XNPX,SLOPE(0,N2,N1),ZM(N2,N1)
C    3,XVOLW(N2,N1),XVOLT(N2,N1),VOLWG(N2,N1),VOLW1(0,N2,N1)
C    4,WFLXR(N2,N1)
5554  FORMAT(A8,5I4,20E12.4)
C     ENDIF
      ELSE
      QRM(M,N2,N1)=0.0
      QRV(M,N2,N1)=0.0
      ENDIF
C
C     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
C
      DO 4310 N=1,2
      DO 4305 NN=1,2
      IF(N.EQ.1)THEN
      IF(NX.EQ.NHE.AND.NN.EQ.1
     2.OR.NX.EQ.NHW.AND.NN.EQ.2)THEN
      GO TO 4305
      ELSE
      N4=NX+1
      N5=NY
      N4B=NX-1
      N5B=NY
      ENDIF
      ELSEIF(N.EQ.2)THEN
      IF(NY.EQ.NVS.AND.NN.EQ.1
     2.OR.NY.EQ.NVN.AND.NN.EQ.2)THEN
      GO TO 4305
      ELSE
      N4=NX
      N5=NY+1
      N4B=NX
      N5B=NY-1
      ENDIF
      ENDIF
C
C     ELEVATION OF EACH PAIR OF ADJACENT GRID CELLS
C
C     XVOLT,XVOLW=excess water+ice,water in destination grid cell
C     ALT1,ALT2=elevation of source,destination
C     QRQ1=equilibrium runoff
C     QR1,HQR1=runoff, convective heat from runoff
C     QR,HQR=hourly-accumulated runoff, convective heat from runoff
C
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
      ALT1=ALTG(N2,N1)+XVOLT(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
C
C     EAST OR SOUTH RUNOFF
C
      IF(NN.EQ.1)THEN
      ALT2=ALTG(N5,N4)+XVOLT(N5,N4)/AREA(3,NU(N5,N4),N5,N4)
      IF(ALT1.GT.ALT2)THEN
      QRQ1=AMAX1(0.0,((ALT1-ALT2)*AREA(3,NUM(N2,N1),N2,N1)
     2*AREA(3,NU(N5,N4),N5,N4)-XVOLT(N5,N4)*AREA(3,NUM(N2,N1),N2,N1)
     3+XVOLT(N2,N1)*AREA(3,NU(N5,N4),N5,N4))
     4/(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5,N4),N5,N4)))
      QR1(N,2,N5,N4)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
      HQR1(N,2,N5,N4)=4.19*TK1(0,N2,N1)*QR1(N,2,N5,N4)
      QR(N,2,N5,N4)=QR(N,2,N5,N4)+QR1(N,2,N5,N4)
      HQR(N,2,N5,N4)=HQR(N,2,N5,N4)+HQR1(N,2,N5,N4)
      QRMN(M,N,2,N5,N4)=QR1(N,2,N5,N4)
      IFLBM(M,N,2,N5,N4)=0
C     IF(NX.EQ.1.AND.NY.EQ.1)THEN
C     WRITE(*,5555)'QRFOR',I,J,M,N1,N2,N4,N5,N,NN
C    2,QRM(M,N2,N1),QR1(N,2,N5,N4),QR(N,2,N5,N4)
C    2,ALT1,ALT2,ALTG(N2,N1),ALTG(N5,N4),ALT(N2,N1),ALT(N5,N4)
C    3,QRQ1,FSLOPE(N,N2,N1),QR1(2,2,4,1)
5555  FORMAT(A8,9I4,30E12.4)
C     ENDIF
      ELSE
      QR1(N,2,N5,N4)=0.0
      HQR1(N,2,N5,N4)=0.0
      QRMN(M,N,2,N5,N4)=0.0
      IFLBM(M,N,2,N5,N4)=1
      ENDIF
      ENDIF
C
C     WEST OR NORTH RUNOFF
C
      IF(NN.EQ.2)THEN
      IF(N4B.GT.0.AND.N5B.GT.0)THEN
      ALTB=ALTG(N5B,N4B)+XVOLT(N5B,N4B)/AREA(3,NU(N5,N4B),N5B,N4B)
      IF(ALT1.GT.ALTB)THEN
      QRQ1=AMAX1(0.0,((ALT1-ALTB)*AREA(3,NUM(N2,N1),N2,N1)
     2*AREA(3,NU(N5B,N4B),N5B,N4B)-XVOLT(N5B,N4B)
     2*AREA(3,NUM(N2,N1),N2,N1)
     3+XVOLT(N2,N1)*AREA(3,NU(N5B,N4B),N5B,N4B))
     4/(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5B,N4B),N5B,N4B)))
      QR1(N,1,N5B,N4B)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
      HQR1(N,1,N5B,N4B)=4.19*TK1(0,N2,N1)*QR1(N,1,N5B,N4B)
      QR(N,1,N5B,N4B)=QR(N,1,N5B,N4B)+QR1(N,1,N5B,N4B)
      HQR(N,1,N5B,N4B)=HQR(N,1,N5B,N4B)+HQR1(N,1,N5B,N4B)
      QRMN(M,N,1,N5B,N4B)=QR1(N,1,N5B,N4B)
      IFLBM(M,N,1,N5B,N4B)=1
C     WRITE(*,5555)'QRBAK',I,J,M,N1,N2,N4B,N5B,N,NN
C    2,QRM(M,N2,N1),QR1(N,1,N5B,N4B),QR(N,1,N5B,N4B)
C    2,ALT1,ALTB,ALTG(N2,N1),ALTG(N5B,N4B),QRQ1,FSLOPE(N,N2,N1)
      ELSE
      QR1(N,1,N5B,N4B)=0.0
      HQR1(N,1,N5B,N4B)=0.0
      QRMN(M,N,1,N5B,N4B)=0.0
      IFLBM(M,N,1,N5B,N4B)=0
      ENDIF
      ENDIF
      ENDIF
      ELSE
      QR1(N,2,N5,N4)=0.0
      HQR1(N,2,N5,N4)=0.0
      QRMN(M,N,2,N5,N4)=0.0
      IFLBM(M,N,2,N5,N4)=0
      IF(N4B.GT.0.AND.N5B.GT.0)THEN
      QR1(N,1,N5B,N4B)=0.0
      HQR1(N,1,N5B,N4B)=0.0
      QRMN(M,N,1,N5B,N4B)=0.0
      IFLBM(M,N,1,N5B,N4B)=0
      ENDIF
      ENDIF
C     WRITE(*,5557)'QRFORA',I,J,M,N1,N2,N4,N5,N,NN,IFLBM(M,N,2,N5,N4)
C    2,QRM(M,N2,N1),QRMN(M,N,2,N5,N4),QR1(N,2,N5,N4),QR(N,2,N5,N4)
5557  FORMAT(A8,10I4,30E12.4)
C     IF(N4B.GT.0.AND.N5B.GT.0)THEN
C     WRITE(*,5557)'QRBAKA',I,J,M,N1,N2,N4B,N5B,N,NN
C    2,IFLBM(M,N,1,N5B,N4B)
C    2,QRM(M,N2,N1),QRMN(M,N,1,N5B,N4B),QR1(N,1,N5B,N4B)
C    3,QR(N,1,N5B,N4B)
C     ENDIF
C
C     SNOW REDISTRIBUTION FROM SNOWPACK
C
C     N2,N1=NY,NX of source grid cell
C     N5,N4=NY,NX of destination grid cell
C     ALTS1,ALTS2=elevation of source,destination snowpack surfaces
C     SS,DIST=slope,distance between source,destination
C     QSX=transfer fraction
C     QS1,QW1,QI1=snow,water,ice transfer
C     HQS1=convective heat transfer from snow,water,ice transfer
C     VOLS0,VOLW0,VOLI0=snow,water,ice volume
C     DPTHSX=minimum snowpack depth for full cover
C     QS,QW,QI=hourly-accumulated snow,water,ice transfer
C     HQS=hourly-accumd convective heat from snow,water,ice transfer
C     QSM=snow transfer for solute flux calculation
C
      IF(NN.EQ.1)THEN
      ALTS1=ALTG(N2,N1)+DPTHS(N2,N1)
      ALTS2=ALTG(N5,N4)+DPTHS(N5,N4)
      SS=(ALTS1-ALTS2)/DIST(N,NU(N5,N4),N5,N4)
      QSX=SS/AMAX1(1.0,DIST(N,NU(N5,N4),N5,N4))*XNPH
      IF(SS.GT.0.0.AND.DPTHS(N2,N1).GT.DPTHSX)THEN
      QS1(N,N5,N4)=QSX*VOLS0(1,N2,N1)
      QW1(N,N5,N4)=QSX*VOLW0(1,N2,N1)
      QI1(N,N5,N4)=QSX*VOLI0(1,N2,N1)
      HQS1(N,N5,N4)=TK0(1,N2,N1)*(2.095*QS1(N,N5,N4)
     2+4.19*QW1(N,N5,N4)+1.9274*QI1(N,N5,N4))
      ELSEIF(SS.LT.0.0.AND.DPTHS(N5,N4).GT.DPTHSX)THEN
      QS1(N,N5,N4)=QSX*VOLS0(1,N5,N4)
      QW1(N,N5,N4)=QSX*VOLW0(1,N5,N4)
      QI1(N,N5,N4)=QSX*VOLI0(1,N5,N4)
      HQS1(N,N5,N4)=TK0(1,N5,N4)*(2.095*QS1(N,N5,N4)
     2+4.19*QW1(N,N5,N4)+1.9274*QI1(N,N5,N4))
      ELSE
      QS1(N,N5,N4)=0.0
      QW1(N,N5,N4)=0.0
      QI1(N,N5,N4)=0.0
      HQS1(N,N5,N4)=0.0
      ENDIF
      QS(N,N5,N4)=QS(N,N5,N4)+QS1(N,N5,N4)
      QW(N,N5,N4)=QW(N,N5,N4)+QW1(N,N5,N4)
      QI(N,N5,N4)=QI(N,N5,N4)+QI1(N,N5,N4)
      HQS(N,N5,N4)=HQS(N,N5,N4)+HQS1(N,N5,N4)
      QSM(M,N,N5,N4)=QS1(N,N5,N4)
C     IF(NX.EQ.2.AND.NY.EQ.5)THEN
C     WRITE(*,5556)'QS1',I,J,M,N1,N2,N4,N5,N,QSX,QS1(N,N5,N4)
C    2,QW1(N,N5,N4),QI1(N,N5,N4),VOLS0(N3,N2,N1),VOLW0(N3,N2,N1)
C    3,VOLI0(N3,N2,N1),ALTS1,ALTS2,ALTG(N2,N1),ALTG(N5,N4)
C    4,DIST(N,NU(N5,N4),N5,N4),SS,DLYRS0(N2,N1),DLYRS0(N5,N4)
C    5,VOLS1(N2,N1),VOLS1(N5,N4),VOLWG(N2,N1),VOLWG(N5,N4)
5556  FORMAT(A8,8I4,30E12.4)
C     ENDIF
      ENDIF
4305  CONTINUE
4310  CONTINUE
C
C     HOURLY-ACCUMULATED WATER, VAPOR AND HEAT FLUXES THROUGH
C     SURFACE RESIDUE AND SOIL SURFACE
C
C     THAWR,HTHAWR=litter water,heat fluxes from freeze-thaw
C     FLW,FLWH,HFLW=soil micropore,macropore,heat fluxes
C     FLWR,HFLWR=litter water,heat fluxes
C     FLSW,FLSWH=water from snowpack to soil micropores,macropores
C     HFLSW=convective heat from snowpack to soil
C     FLSWR=water flux from snowpack to litter
C     HFLSWR=convective heat flux from snowpack to litter
C     HEATI,HEATE,HEATS,HEATG=net radiation,latent,sensible,storage heat
C     TEVAPG=total evaporation
C     FLWM,FLWHM=water flux into soil micropore,macropore for use in trnsfr.f
C     VOLWX1=VOLW1 accounting for wetting front
C
      THAWR(NY,NX)=THAWR(NY,NX)+WFLXR(NY,NX)
      HTHAWR(NY,NX)=HTHAWR(NY,NX)+TFLXR(NY,NX)
      FLW(3,NUM(NY,NX),NY,NX)=FLW(3,NUM(NY,NX),NY,NX)
     2+FLWL(3,NUM(NY,NX),NY,NX)
      FLWX(3,NUM(NY,NX),NY,NX)=FLWX(3,NUM(NY,NX),NY,NX)
     2+FLWLX(3,NUM(NY,NX),NY,NX)
      FLWH(3,NUM(NY,NX),NY,NX)=FLWH(3,NUM(NY,NX),NY,NX)
     2+FLWHL(3,NUM(NY,NX),NY,NX)
      HFLW(3,NUM(NY,NX),NY,NX)=HFLW(3,NUM(NY,NX),NY,NX)
     2+HFLWL(3,NUM(NY,NX),NY,NX)
      FLWR(NY,NX)=FLWR(NY,NX)+FLWRL(NY,NX)
      HFLWR(NY,NX)=HFLWR(NY,NX)+HFLWRL(NY,NX)
      HEATI(NY,NX)=HEATI(NY,NX)+RFLXG+RFLXR+RFLXW
      HEATS(NY,NX)=HEATS(NY,NX)+SFLXG+SFLXR+SFLXW
      HEATE(NY,NX)=HEATE(NY,NX)+EFLXG+EFLXR+EFLXW
      HEATV(NY,NX)=HEATV(NY,NX)+VFLXG+VFLXR+VFLXW
      HEATH(NY,NX)=HEATH(NY,NX)+RFLXG+RFLXR+RFLXW
     2+SFLXG+SFLXR+SFLXW+EFLXG+EFLXR+EFLXW+VFLXG+VFLXR+VFLXW
      TEVAPG(NY,NX)=TEVAPG(NY,NX)+EVAPG(NY,NX)+EVAPR(NY,NX)
     2+EVAPS(NY,NX)+EVAPW(NY,NX)
      FLWM(M,3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)
      FLWHM(M,3,NUM(NY,NX),NY,NX)=FLWHL(3,NUM(NY,NX),NY,NX)
C
C     WATER AND ENERGY TRANSFER THROUGH SOIL PROFILE
C
C     N3,N2,N1=L,NY,NX of source grid cell
C     N6,N5,N4=L,NY,NX of destination grid cell
C
      IFLGH=0
      DO 4400 L=1,NL(NY,NX)
      N1=NX
      N2=NY
      N3=L
C
C     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
C
      DO 4320 N=NCN(N2,N1),3
      IF(N.EQ.1)THEN
      IF(NX.EQ.NHE)THEN
      GO TO 4320
      ELSE
      N4=NX+1
      N5=NY
      N6=L
C
C     ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
C
C     IF(N2.EQ.2.AND.(N1.EQ.2.OR.N1.EQ.3).AND.L.LE.15)THEN
C     GO TO 4320
C     ENDIF
      ENDIF
      ELSEIF(N.EQ.2)THEN
      IF(NY.EQ.NVS)THEN
      GO TO 4320
      ELSE
      N4=NX
      N5=NY+1
      N6=L
C
C     ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
C
C     IF(N1.EQ.3.AND.(N2.EQ.1.OR.N2.EQ.2).AND.L.LE.15)THEN
C     GO TO 4320
C     ENDIF
C
C     END ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
C
      ENDIF
      ELSEIF(N.EQ.3)THEN
      IF(L.EQ.NL(NY,NX))THEN
      GO TO 4320
      ELSE
      N4=NX
      N5=NY
      N6=L+1
      ENDIF
      ENDIF
C
C     SKIP NON-EXISTENT DESTINATION SOIL LAYERS
C
      DO 1100 LL=N6,NL(NY,NX)
      IF(VOLX(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
      N6=LL
      GO TO 1101
      ENDIF
1100  CONTINUE
1101  CONTINUE
      IF(N3.EQ.NU(N2,N1))N6X(N2,N1)=N6
C
C     POROSITIES 'THETP*', WATER CONTENTS 'THETA*', AND POTENTIALS
C     'PSIS*' FOR EACH GRID CELL
C
C     THETA1,THETAL=micropore water concn in source,destination cells
C     THETY=hygroscopic water concentration
C     POROS=soil porosity
C     VOLXI=soil volume excluding rock, macropore
C
      IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4)
     2.AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN
      THETA1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1)
     2,VOLW1(N3,N2,N1)/VOLY(N3,N2,N1)))
      THETAL=AMAX1(THETY(N6,N5,N4),AMIN1(POROS(N6,N5,N4)
     2,VOLW1(N6,N5,N4)/VOLY(N6,N5,N4)))
C
C     WATER POTENTIAL OF UPPER LAYER
C
C     BKVL=soil mass
C     FC,WP=water contents at field capacity,wilting point
C     FCL,WPL=log FC,WP
C     FCD,PSD=FCL-WPL,log(POROS)-FCL
C     PSISA1,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
C     PSIMX,PSIMD,PSIMS=log water potential at FC,WP,saturation
C     PSISD=PSIMX-PSIMS
C     SRP=parameter for deviation from linear log-log water retention
C     PSISO=osmotic potential
C
      IF(BKVL(N3,N2,N1).GT.ZEROS(NY,NX))THEN
      IF(THETA1.LT.FC(N3,N2,N1))THEN
      PSISA1(N3,N2,N1)=AMAX1(PSIHY,-EXP(PSIMX(N2,N1)
     2+((FCL(N3,N2,N1)-LOG(THETA1))
     3/FCD(N3,N2,N1)*PSIMD(N2,N1))))
      ELSEIF(THETA1.LT.POROS(N3,N2,N1)-DTHETW)THEN
      PSISA1(N3,N2,N1)=-EXP(PSIMS(N2,N1)
     2+(((PSL(N3,N2,N1)-LOG(THETA1))
     3/PSD(N3,N2,N1))**SRP(N3,N2,N1)*PSISD(N2,N1)))
      ELSE
      PSISA1(N3,N2,N1)=PSISE(N3,N2,N1)
      ENDIF
C
C     SUBSURFCE UPPER WATER LAYER
C
C     THETIX,THETWX=ice,water concentration
C     FCI,WPI=ice field capacity,wilting point
C     PSISA1=matric water potential
C
      ELSEIF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      FCX=FCI*THETIX(N3,N2,N1)
      WPX=WPI*THETIX(N3,N2,N1)
      FCLX=LOG(FCX)
      WPLX=LOG(WPX)
      PSDX=PSL(N3,N2,N1)-FCLX
      FCDX=FCLX-WPLX
      IF(THETWX(N3,N2,N1).LT.FCX)THEN
      PSISA1(N3,N2,N1)=AMAX1(PSIHY,-EXP(PSIMX(N2,N1)
     2+((FCLX-LOG(THETWX(N3,N2,N1)))
     3/FCDX*PSIMD(NY,NX))))
      ELSEIF(THETWX(N3,N2,N1).LT.POROS(N3,N2,N1)-DTHETW)THEN
      PSISA1(N3,N2,N1)=-EXP(PSIMS(N2,N1)
     2+(((PSL(N3,N2,N1)-LOG(THETWX(N3,N2,N1)))
     3/PSDX)*PSISD(N2,N1)))
      ELSE
      PSISA1(N3,N2,N1)=PSISE(N3,N2,N1)
      ENDIF
      ELSE
      PSISA1(N3,N2,N1)=PSISE(N3,N2,N1)
      ENDIF
C     IF(N1.EQ.4.AND.N2.EQ.1.AND.N3.EQ.11)THEN
C     WRITE(*,1119)'PSISA1',I,J,M,N,N1,N2,N3,PSISA1(N3,N2,N1)
C    2,THETWX(N3,N2,N1),THETIX(N3,N2,N1),FCX,WPX
C    3,BKVL(N3,N2,N1),VOLX(N3,N2,N1)
C     ENDIF
C
C     WATER POTENTIAL OF LOWER LAYER
C
      IF(BKVL(N6,N5,N4).GT.ZEROS(NY,NX))THEN
      IF(THETAL.LT.FC(N6,N5,N4))THEN
      PSISA1(N6,N5,N4)=AMAX1(PSIHY,-EXP(PSIMX(N5,N4)
     2+((FCL(N6,N5,N4)-LOG(THETAL))
     3/FCD(N6,N5,N4)*PSIMD(N5,N4))))
      ELSEIF(THETAL.LT.POROS(N6,N5,N4)-DTHETW)THEN
      PSISA1(N6,N5,N4)=-EXP(PSIMS(N5,N4)
     2+(((PSL(N6,N5,N4)-LOG(THETAL))
     3/PSD(N6,N5,N4))**SRP(N6,N5,N4)*PSISD(N5,N4)))
      ELSE
      PSISA1(N6,N5,N4)=PSISE(N6,N5,N4)
      ENDIF
C
C     SUBSURFCE LOWER WATER LAYER
C
      ELSEIF(VOLX(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      FCX=FCI*THETIX(N6,N5,N4)
      WPX=WPI*THETIX(N6,N5,N4)
      FCLX=LOG(FCX)
      WPLX=LOG(WPX)
      PSDX=PSL(N6,N5,N4)-FCLX
      FCDX=FCLX-WPLX
      IF(THETWX(N6,N5,N4).LT.FCX)THEN
      PSISA1(N6,N5,N4)=AMAX1(PSIHY,-EXP(PSIMX(N5,N4)
     2+((FCLX-LOG(THETWX(N6,N5,N4)))
     3/FCDX*PSIMD(NY,NX))))
      ELSEIF(THETWX(N6,N5,N4).LT.POROS(N6,N5,N4)-DTHETW)THEN
      PSISA1(N6,N5,N4)=-EXP(PSIMS(NY,NX)
     2+(((PSL(N6,N5,N4)-LOG(THETWX(N6,N5,N4)))
     3/PSDX)*PSISD(NY,NX)))
      ELSE
      PSISA1(N6,N5,N4)=PSISE(N6,N5,N4)
      ENDIF
      ELSE
      PSISA1(N6,N5,N4)=PSISE(N6,N5,N4)
      ENDIF
C     IF(N1.EQ.4.AND.N2.EQ.1.AND.N3.EQ.11)THEN
C     WRITE(*,1119)'PSISAL',I,J,M,N,N4,N5,N6,PSISA1(N6,N5,N4)
C    2,THETWX(N6,N5,N4),THETIX(N6,N5,N4),FCX,WPX
C    3,BKVL(N6,N5,N4),VOLX(N6,N5,N4)
C     ENDIF
C
C     ACCOUNT FOR WETTING FRONTS WHEN CALCULATING WATER CONTENTS,
C     MATRIC WATER POTENTIALS AND HYDRAULIC CONDUCTIVITIES USED
C     IN WATER FLUX CALCULATIONS
C
C     THETW1,THETWL=water concentrations in source,destination cells
C     CND1,CNDL=hydraulic conductivities in source,destination cells
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
C     PSISM1=soil matric potential
C     VOLWX1=VOLW1 accounting for wetting front
C
C     DARCY FLOW IF BOTH CELLS ARE SATURATED
C     (CURRENT WATER POTENTIAL > AIR ENTRY WATER POTENTIAL)
C
      IF(PSISA1(N3,N2,N1).GT.PSISA(N3,N2,N1)
     2.AND.PSISA1(N6,N5,N4).GT.PSISA(N6,N5,N4))THEN
      THETW1=THETA1
      THETWL=THETAL
      K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)
     2-THETW1)/POROS(N3,N2,N1))+1))
      KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4)
     2-THETWL)/POROS(N6,N5,N4))+1))
      PSISM1(N3,N2,N1)=PSISA1(N3,N2,N1)
      PSISM1(N6,N5,N4)=PSISA1(N6,N5,N4)
C
C     GREEN-AMPT FLOW IF ONE LAYER IS SATURATED
C     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POENTIAL)
C
C     GREEN-AMPT FLOW IF SOURCE CELL SATURATED
C
      ELSEIF(PSISA1(N3,N2,N1).GT.PSISA(N3,N2,N1))THEN
      THETW1=THETA1
      THETWL=AMAX1(THETY(N6,N5,N4),AMIN1(POROS(N6,N5,N4)
     2,VOLWX1(N6,N5,N4)/VOLY(N6,N5,N4)))
      K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)
     2-THETW1)/POROS(N3,N2,N1))+1))
      KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4)
     2-AMIN1(THETS(N6,N5,N4),THETWL))/POROS(N6,N5,N4))+1))
      PSISM1(N3,N2,N1)=PSISA1(N3,N2,N1)
      IF(BKVL(N6,N5,N4).GT.ZEROS(NY,NX))THEN
      IF(THETWL.LT.FC(N6,N5,N4))THEN
      PSISM1(N6,N5,N4)=AMAX1(PSIHY,-EXP(PSIMX(N5,N4)
     2+((FCL(N6,N5,N4)-LOG(THETWL))
     3/FCD(N6,N5,N4)*PSIMD(N5,N4))))
      ELSEIF(THETWL.LT.POROS(N6,N5,N4)-DTHETW)THEN
      PSISM1(N6,N5,N4)=-EXP(PSIMS(N5,N4)
     2+(((PSL(N6,N5,N4)-LOG(THETWL))
     3/PSD(N6,N5,N4))**SRP(N6,N5,N4)*PSISD(N5,N4)))
      ELSE
      THETWL=POROS(N6,N5,N4)
      PSISM1(N6,N5,N4)=PSISE(N6,N5,N4)
      ENDIF
      ELSE
      THETWL=POROS(N6,N5,N4)
      PSISM1(N6,N5,N4)=PSISE(N6,N5,N4)
      ENDIF
C     IF(N3.EQ.NUM(NY,NX))THEN
C     WRITE(*,1116)'GA',I,J,M,N1,N2,N3,N4,N5,N6,N
C    3,PSISM1(N3,N2,N1),PSISM1(N6,N5,N4)
C    3,VOLW1(N3,N2,N1),VOLW1(N6,N5,N4)
C    3,VOLWX1(N3,N2,N1),VOLWX1(N6,N5,N4)
C    6,THETW1,THETWL
1116  FORMAT(A8,10I4,100E12.4)
C     ENDIF
C
C     GREEN-AMPT FLOW IF ADJACENT CELL SATURATED
C
      ELSEIF(PSISA1(N6,N5,N4).GT.PSISA(N6,N5,N4))THEN
      THETW1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1)
     2,VOLWX1(N3,N2,N1)/VOLY(N3,N2,N1)))
      THETWL=THETAL
      K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)
     2-AMIN1(THETS(N3,N2,N1),THETW1))/POROS(N3,N2,N1))+1))
      KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4)
     2-THETWL)/POROS(N6,N5,N4))+1))
      IF(BKVL(N3,N2,N1).GT.ZEROS(NY,NX))THEN
      IF(THETW1.LT.FC(N3,N2,N1))THEN
      PSISM1(N3,N2,N1)=AMAX1(PSIHY,-EXP(PSIMX(N2,N1)
     2+((FCL(N3,N2,N1)-LOG(THETW1))
     3/FCD(N3,N2,N1)*PSIMD(N2,N1))))
      ELSEIF(THETW1.LT.POROS(N3,N2,N1)-DTHETW)THEN
      PSISM1(N3,N2,N1)=-EXP(PSIMS(N2,N1)
     2+(((PSL(N3,N2,N1)-LOG(THETW1))
     3/PSD(N3,N2,N1))**SRP(N3,N2,N1)*PSISD(N2,N1)))
      ELSE
      THETW1=POROS(N3,N2,N1)
      PSISM1(N3,N2,N1)=PSISE(N3,N2,N1)
      ENDIF
      ELSE
      THETW1=POROS(N3,N2,N1)
      PSISM1(N3,N2,N1)=PSISE(N3,N2,N1)
      ENDIF
C
C     RICHARDS FLOW IF NEITHER CELL IS SATURATED
C     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POTENTIAL)
C
      ELSE
      THETW1=THETA1
      THETWL=THETAL
      K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)
     2-THETW1)/POROS(N3,N2,N1))+1))
      KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4)
     2-THETWL)/POROS(N6,N5,N4))+1))
      PSISM1(N3,N2,N1)=PSISA1(N3,N2,N1)
      PSISM1(N6,N5,N4)=PSISA1(N6,N5,N4)
      ENDIF
C
C     HYDRAULIC CONUCTIVITY
C
C     CND1,CNDL=hydraulic conductivity of source,destination layer
C     HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity
C
      IF(N3.EQ.NUM(NY,NX))THEN
      CND1=HCND(N,K1,N3,N2,N1)*FKSAT
      ELSE
      CND1=HCND(N,K1,N3,N2,N1)
      ENDIF
      CNDL=HCND(N,KL,N6,N5,N4)
C
C     TOTAL SOIL WATER POTENTIAL = MATRIC, GRAVIMETRIC + OSMOTIC
C
C     PSISM1,PSISH,PSISO=soil matric,gravitational,osmotic potentials
C
      PSIST1=PSISM1(N3,N2,N1)+PSISH(N3,N2,N1)+PSISO(N3,N2,N1)
      PSISTL=PSISM1(N6,N5,N4)+PSISH(N6,N5,N4)+PSISO(N6,N5,N4)
      PSISV1=PSISM1(N3,N2,N1)+PSISO(N3,N2,N1)
      PSISVL=PSISM1(N6,N5,N4)+PSISO(N6,N5,N4)
C     IF(N6.EQ.12)THEN
C     WRITE(*,7272)'PSIM',I,J,M,N1,N2,N3,N4,N5,N6
C    2,PSISM1(N3,N2,N1),PSISM1(N6,N5,N4),PSISA1(N3,N2,N1)
C    3,PSISA1(N6,N5,N4),THETWL,THETAL
C    2,PSIMX(N5,N4),FCL(N6,N5,N4),FCD(N6,N5,N4),PSIMD(N5,N4)
C    3,POROS(N6,N5,N4),PSIMS(N5,N4),PSL(N6,N5,N4),PSD(N6,N5,N4)
C    4,SRP(N6,N5,N4),PSISD(N5,N4),PSISE(N6,N5,N4)
C    5,THETY(N6,N5,N4),POROS(N6,N5,N4),VOLW1(N6,N5,N4),VOLX(N6,N5,N4)
7272  FORMAT(A8,9I4,30E12.4)
C     ENDIF
C
C     HYDRAULIC CONDUCTIVITY FROM CURRENT WATER CONTENT
C     AND LOOKUP ARRAY GENERATED IN 'HOUR1'
C
C     CND1,CNDL=hydraulic conductivities in source,destination cells
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
C     AVCNDL=source-destination hydraulic conductance
C     DLYR=layer thickness
C
      IF(CND1.GT.ZERO.AND.CNDL.GT.ZERO)THEN
      AVCNDL=2.0*CND1*CNDL/(CND1*DLYR(N,N6,N5,N4)
     2+CNDL*DLYR(N,N3,N2,N1))
      ELSE
      AVCNDL=0.0
      ENDIF
C
C     WATER FLUX FROM WATER POTENTIALS, HYDRAULIC CONDUCTIVITY
C     CONSTRAINED BY WATER POTENTIAL GRADIENT, COUPLED WITH
C     CONVECTIVE HEAT FLUX FROM WATER FLUX
C
C     FLQX,FLQL=micropore water flux unlimited,limited by source water
C     XNPH=time step of flux calculations
C     VOLW2,VOLP1=water,air contents of source,destination micropores
C     HWFLWL=convective heat flux from micropore water flux
C     VOLP1Z=excess water+ice relative to porosity
C
      FLQX=AVCNDL*(PSIST1-PSISTL)*AREA(N,N3,N2,N1)*XNPH
      IF(FLQX.GE.0.0)THEN
      IF(THETW1.GT.THETS(N3,N2,N1))THEN
      FLQZ=FLQX+AMIN1((THETW1-THETS(N3,N2,N1))
     2*VOLY(N3,N2,N1),AMAX1(0.0,(THETS(N6,N5,N4)-THETWL)
     3*VOLY(N6,N5,N4)))*XNPX
      ELSE
      FLQZ=FLQX
      ENDIF
      FLQL=AMAX1(0.0,AMIN1(FLQZ,VOLW2(N3,N2,N1)*XNPX
     2,VOLP1(N6,N5,N4)*XNPX))
      FLQ2=AMAX1(0.0,AMIN1(FLQX,VOLW2(N3,N2,N1)*XNPX
     2,VOLP1(N6,N5,N4)*XNPX))
C     FLQL1=(THETW1-THETS(N3,N2,N1))*VOLY(N3,N2,N1)
C     FLQL2=(THETS(N6,N5,N4)-THETWL)*VOLY(N6,N5,N4)
C     FLQL3=FLQX+AMIN1(FLQL1,AMAX1(0.0,FLQL2))*XNPX
C     FLQL4=AMAX1(0.0,AMIN1(FLQL3,VOLP1(N6,N5,N4)*XNPX))
      ELSE
      IF(THETWL.GT.THETS(N6,N5,N4))THEN
      FLQZ=FLQX+AMAX1((THETS(N6,N5,N4)-THETWL)
     2*VOLY(N6,N5,N4),AMIN1(0.0,(THETW1-THETS(N3,N2,N1))
     3*VOLY(N3,N2,N1)))*XNPX
      ELSE
      FLQZ=FLQX
      ENDIF
      FLQL=AMIN1(0.0,AMAX1(FLQZ,-VOLW2(N6,N5,N4)*XNPX
     2,-VOLP1(N3,N2,N1)*XNPX))
      FLQ2=AMIN1(0.0,AMAX1(FLQX,-VOLW2(N6,N5,N4)*XNPX
     2,-VOLP1(N3,N2,N1)*XNPX))
C     FLQL1=(THETS(N6,N5,N4)-THETWL)*VOLY(N6,N5,N4)
C     FLQL2=(THETW1-THETS(N3,N2,N1))*VOLY(N3,N2,N1)
C     FLQL3=FLQX+AMAX1(FLQL1,AMIN1(0.0,FLQL2))*XNPX
C     FLQL4=AMIN1(0.0,AMAX1(FLQL3,-VOLP1(N3,N2,N1)*XNPX))
      ENDIF
      IF(N.EQ.3.AND.VOLP1Z(N6,N5,N4).LT.0.0)THEN
      FLQL=FLQL+AMIN1(0.0,AMAX1(-VOLW2(N6,N5,N4)*XNPX
     2,VOLP1Z(N6,N5,N4)))
      FLQ2=FLQ2+AMIN1(0.0,AMAX1(-VOLW2(N6,N5,N4)*XNPX
     2,VOLP1Z(N6,N5,N4)))
      ENDIF
      IF(FLQL.GT.0.0)THEN
      HWFLQL=4.19*TK1(N3,N2,N1)*FLQL
      ELSE
      HWFLQL=4.19*TK1(N6,N5,N4)*FLQL
      ENDIF
      VOLW2(N3,N2,N1)=VOLW2(N3,N2,N1)-FLQL
      VOLW2(N6,N5,N4)=VOLW2(N6,N5,N4)+FLQL
C
C     MACROPORE FLOW FROM POISEUILLE FLOW IF MACROPORES PRESENT
C
C     PSISH1,PSISHL=macropore total water potl in source,destination
C     DLYR=layer thickness
C     VOLWH1,VOLPH1=macropore water,air content
C
      IF(VOLAH1(N3,N2,N1).GT.ZEROS2(N2,N1)
     2.AND.VOLAH1(N6,N5,N4).GT.ZEROS2(N5,N4).AND.IFLGH.EQ.0)THEN
      PSISH1=PSISH(N3,N2,N1)+0.0098*DLYR(3,N3,N2,N1)
     2*(AMIN1(1.0,AMAX1(0.0,VOLWH1(N3,N2,N1)/VOLAH1(N3,N2,N1)))-0.5)
      PSISHL=PSISH(N6,N5,N4)+0.0098*DLYR(3,N6,N5,N4)
     2*(AMIN1(1.0,AMAX1(0.0,VOLWH1(N6,N5,N4)/VOLAH1(N6,N5,N4)))-0.5)
C
C     MACROPORE FLOW IF GRAVITATIONAL GRADIENT IS POSITIVE
C     AND MACROPORE POROSITY EXISTS IN ADJACENT CELL
C
C     FLWHX,FLWHL=macropore water flux unltd,ltd by source water
C     XNPH=time step of flux calculations
C     VOLW2,VOLP1=water,air contents of source,destination micropores
C     HWFLHL=convective heat flux from micropore water flux
C
      FLWHX=AVCNHL(N,N6,N5,N4)*(PSISH1-PSISHL)*AREA(N,N3,N2,N1)*XNPH
      IF(N.NE.3)THEN
      IF(PSISH1.GT.PSISHL)THEN
      FLWHL(N,N6,N5,N4)=AMAX1(0.0,AMIN1(AMIN1(VOLWH1(N3,N2,N1)
     2,VOLPH1(N6,N5,N4))*XNPX,FLWHX))
      ELSEIF(PSISH1.LT.PSISHL)THEN
      FLWHL(N,N6,N5,N4)=AMIN1(0.0,AMAX1(AMAX1(-VOLWH1(N6,N5,N4)
     2,-VOLPH1(N3,N2,N1))*XNPX,FLWHX))
      ELSE
      FLWHL(N,N6,N5,N4)=0.0
      ENDIF
      ELSE
      FLWHL(N,N6,N5,N4)=AMAX1(0.0,AMIN1(AMIN1(VOLWH1(N3,N2,N1)*XNPX
     2+FLWHL(N,N3,N2,N1),VOLPH1(N6,N5,N4)*XNPX),FLWHX))
      ENDIF
      IF(N.EQ.3)THEN
      FLWHL(N,N6,N5,N4)=FLWHL(N,N6,N5,N4)+AMIN1(0.0,VOLPH1Z(N6,N5,N4))
      ENDIF
      FLWHM(M,N,N6,N5,N4)=FLWHL(N,N6,N5,N4)
C     IF(N4.EQ.1)THEN
C     WRITE(*,5478)'FLWH',I,J,M,N1,N2,N3,IFLGH
C    2,FLHM,FLWHX,FLWHL(N,N3,N2,N1),FLWHL(N,N6,N5,N4)
C    2,AVCNHL(N,N6,N5,N4),PSISH(N3,N2,N1),PSISH(N6,N5,N4)
C    3,VOLPH1(N3,N2,N1),VOLPH1(N6,N5,N4),VOLWH1(N3,N2,N1)
C    4,VOLWH1(N6,N5,N4),VOLAH1(N3,N2,N1),VOLAH1(N6,N5,N4)
C    5,DLYR(N,N6,N5,N4),DLYR(N,N3,N2,N1),AREA(N,N3,N2,N1)
C    7,CNDH1(N3,N2,N1),CNDH1(N6,N5,N4),XNPH,HWFLHL
5478  FORMAT(A8,7I4,30E12.4)
C     ENDIF
      ELSE
      FLWHL(N,N6,N5,N4)=0.0
      FLWHM(M,N,N6,N5,N4)=0.0
      IF(VOLPH1(N6,N5,N4).LE.0.0)IFLGH=1
      ENDIF
      IF(FLWHL(N,N6,N5,N4).GT.0.0)THEN
      HWFLHL=4.19*TK1(N3,N2,N1)*FLWHL(N,N6,N5,N4)
      ELSE
      HWFLHL=4.19*TK1(N6,N5,N4)*FLWHL(N,N6,N5,N4)
      ENDIF
C
C     VAPOR PRESSURE AND DIFFUSIVITY IN EACH GRID CELL
C
C     THETPM,THETX=current, minimum air-filled porosity
C     TK11,TK12=interim soil temperature in source,destination
C     VP1,VPL=vapor concentration in source,destination
C     PSISV1,PSISVL=matric+osmotic water potl in source,destination
C     CNV1,CNV2=vapor conductivities of source, destination
C     POROS,POROQ=porosity, tortuosity
C     WGSGL=vapor diffusivity
C     ATCNVL=source,destination vapor conductance
C     DLYR=soil layer depth
C     FLVC,FLVX=vapor flux unlimited,limited by vapor
C     VPY=equilibrium vapor concentration
C     XNPX=time step for flux calculations
C     FLVL,HWFLVL=vapor flux and its convective heat flux
C
      IF(THETPM(M,N3,N2,N1).GT.THETX
     2.AND.THETPM(M,N6,N5,N4).GT.THETX)THEN
      TK11=TK1(N3,N2,N1)
      if(abs(tks(n3,n2,n1))>400.)pause
      TK12=TK1(N6,N5,N4)
      VP1=2.173E-03/TK11
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK11))
     3*EXP(18.0*PSISV1/(8.3143*TK11))
      VPL=2.173E-03/TK12
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK12))
     3*EXP(18.0*PSISVL/(8.3143*TK12))
      CNV1=WGSGL(N3,N2,N1)*THETPM(M,N3,N2,N1)*POROQ
     2*THETPM(M,N3,N2,N1)/POROS(N3,N2,N1)
      CNVL=WGSGL(N6,N5,N4)*THETPM(M,N6,N5,N4)*POROQ
     2*THETPM(M,N6,N5,N4)/POROS(N6,N5,N4)
      ATCNVL=2.0*CNV1*CNVL
     2/(CNV1*DLYR(N,N6,N5,N4)+CNVL*DLYR(N,N3,N2,N1))
C
C     VAPOR FLUX FROM VAPOR PRESSURE AND DIFFUSIVITY,
C     AND CONVECTIVE HEAT FLUX FROM VAPOR FLUX
C
      FLVC=ATCNVL*(VP1-VPL)*AREA(N,N3,N2,N1)*XNPH
      VPY=(VP1*VOLPM(M,N3,N2,N1)+VPL*VOLPM(M,N6,N5,N4))
     2/(VOLPM(M,N3,N2,N1)+VOLPM(M,N6,N5,N4))
      FLVX=(VP1-VPY)*VOLPM(M,N3,N2,N1)*XNPX

      IF(FLVC.GE.0.0)THEN
      FLVL=AMAX1(0.0,AMIN1(FLVC,FLVX))
      HWFLVL=(4.19*TK1(N3,N2,N1)+VAP)*FLVL
      ELSE
      FLVL=AMIN1(0.0,AMAX1(FLVC,FLVX))
      HWFLVL=(4.19*TK1(N6,N5,N4)+VAP)*FLVL
      ENDIF
      ELSE
      FLVL=0.0
      HWFLVL=0.0
      ENDIF
C
C     FLWL=total water+vapor flux to destination
C     FLWLX=total unsaturated water+vapor flux to destination
C     HWFLWL=total convective heat flux from water+vapor flux
C
      FLWL(N,N6,N5,N4)=FLQL+FLVL
      FLWLX(N,N6,N5,N4)=FLQ2+FLVL
      HWFLWL=HWFLQL+HWFLVL
C     IF(I.EQ.232.AND.N3.LE.NUM(N2,N1))THEN
C     WRITE(*,1115)'FLWL',I,J,M,N1,N2,N3,N4,N5,N6,N,K1,KL
C    2,FLWL(N,N3,N2,N1),FLWL(N,N6,N5,N4),FLQX,FLQZ,FLQL,FLQ2,FLVL
C    3,PSIST1,PSISTL
C    4,AVCNDL,CND1,CNDL,FKSAT
C    2,VOLP1(N3,N2,N1),VOLP1(N6,N5,N4),VOLW1(N3,N2,N1)
C    3,VOLWX1(N3,N2,N1),VOLW1(N6,N5,N4),VOLWX1(N6,N5,N4)
C    6,THETW1,THETS(N3,N2,N1),THETWL,THETS(N6,N5,N4)
C    7,PSISA1(N3,N2,N1),PSISA1(N6,N5,N4),PSISM1(N3,N2,N1)
C    7,PSISM1(N6,N5,N4),PSISH(N3,N2,N1),PSISH(N6,N5,N4)
C    3,VOLY(N3,N2,N1),VOLY(N6,N5,N4)
C    3,FLVL,FLVX,VP1,VPL,VPY,CNV1,CNVL,ATCNVL
C    4,VOLA1(N6,N5,N4),VOLI1(N6,N5,N4),SCNV(N3,N2,N1)
C    5,SCNV(N6,N5,N4),VOLP1(N3,N2,N1),VOLP1(N6,N5,N4)
C    7,VOLP1Z(N6,N5,N4),VOLT(N3,N2,N1),VOLT(N6,N5,N4)
C    8,DLYR(N,N3,N2,N1),DLYR(N,N6,N5,N4),AREA(N,N3,N2,N1)
C    3,HWFLWL,HWFLQL,HWFLVL
C    9,TK1(N3,N2,N1),TK1(N6,N5,N4),TKY
C    9,VHCP1(N3,N2,N1),VHCP1(N6,N5,N4),POROS(N6,N5,N4)
C    9,VOLP1(N3,N2,N1),VOLX(N3,N2,N1),VOLA1(N6,N5,N4),VOLA1(N3,N2,N1)
C    8,VOLP1(N6,N5,N4),VOLX(N6,N5,N4),FLW(N,N3,N2,N1),FLW(N,N6,N5,N4)
C    1,THETPM(M,N3,N2,N1),THETPM(M,N6,N5,N4),THETX
C    4,FLQL1,FLQL2,FLQL3,FLQL4
1115  FORMAT(A8,12I4,100E12.4)
C     ENDIF
C
C     THERMAL CONDUCTIVITY IN EACH GRID CELL
C
C     DTH*,RYL*,DNU*,TRB*,XNUS*=turbulence effects on thermal conductivity
C     THETWX,THETPX=water,air concentration
C     TCNDW*,TCNDA*=thermal conductivity of water,air
C     TCND1,TCNDL=soil thermal conductivity in source,destination
C     WTHET*=multiplier for air concn in thermal conductivity
C     ATCNDL=source-destination thermal conductance
C
      DTKX=ABS(TK1(N3,N2,N1)-TK1(N6,N5,N4))*1.0E-06
      IF(BKDS(N3,N2,N1).GT.ZERO.OR.THETWX(N3,N2,N1)
     2+THETIX(N3,N2,N1).GT.ZERO)THEN
      DTHW1=AMAX1(0.0,THETWX(N3,N2,N1)-TRBW)**3
      DTHA1=AMAX1(0.0,THETPX(N3,N2,N1)-TRBA)**3
      RYLXW1=DTKX*DTHW1
      RYLXA1=DTKX*DTHA1
      RYLNW1=AMIN1(1.0E+04,RYLXW*RYLXW1)
      RYLNA1=AMIN1(1.0E+04,RYLXA*RYLXA1)
      XNUSW1=AMAX1(1.0,0.68+0.67*RYLNW1**0.25/DNUSW)
      XNUSA1=AMAX1(1.0,0.68+0.67*RYLNA1**0.25/DNUSA)
      TCNDW1=2.067E-03*XNUSW1
      TCNDA1=9.050E-05*XNUSA1
      WTHET1=1.467-0.467*THETPY(N3,N2,N1)
      TCND1=(STC(N3,N2,N1)+THETWX(N3,N2,N1)*TCNDW1
     2+0.611*THETIX(N3,N2,N1)*7.844E-03
     3+WTHET1*THETPX(N3,N2,N1)*TCNDA1)
     4/(DTC(N3,N2,N1)+THETWX(N3,N2,N1)+0.611*THETIX(N3,N2,N1)
     5+WTHET1*THETPX(N3,N2,N1))
      ELSE
      TCND1=0.0
      ENDIF
      IF(BKDS(N6,N5,N4).GT.ZERO.OR.THETWX(N6,N5,N4)
     2+THETIX(N6,N5,N4).GT.ZERO)THEN
      DTHW2=AMAX1(0.0,THETWX(N6,N5,N4)-TRBW)**3
      DTHA2=AMAX1(0.0,THETPX(N6,N5,N4)-TRBA)**3
      RYLXW2=DTKX*DTHW2
      RYLXA2=DTKX*DTHA2
      RYLNW2=AMIN1(1.0E+04,RYLXW*RYLXW2)
      RYLNA2=AMIN1(1.0E+04,RYLXA*RYLXA2)
      XNUSW2=AMAX1(1.0,0.68+0.67*RYLNW2**0.25/DNUSW)
      XNUSA2=AMAX1(1.0,0.68+0.67*RYLNA2**0.25/DNUSA)
      TCNDW2=2.067E-03*XNUSW2
      TCNDA2=9.050E-05*XNUSA2
      WTHET2=1.467-0.467*THETPY(N6,N5,N4)
      TCND2=(STC(N6,N5,N4)+THETWX(N6,N5,N4)*TCNDW2
     2+0.611*THETIX(N6,N5,N4)*7.844E-03
     3+WTHET2*THETPX(N6,N5,N4)*TCNDA2)
     4/(DTC(N6,N5,N4)+THETWX(N6,N5,N4)+0.611*THETIX(N6,N5,N4)
     5+WTHET2*THETPX(N6,N5,N4))
      ELSE
      TCND2=0.0
      ENDIF
      ATCNDL=(2.0*TCND1*TCND2)/(TCND1*DLYR(N,N6,N5,N4)
     3+TCND2*DLYR(N,N3,N2,N1))
C
C     HEAT FLOW FROM THERMAL CONDUCTIVITY AND TEMPERATURE GRADIENT
C
C     VHCP1,VHCPW=volumetric heat capacity of soil,snowpack
C     TK1X,TKLX=interim temperatures of source,destination
C     HWFLVL,HFLXG=convective heat from soil vapor flux
C     HFLXG=storage heat flux from snowpack
C     TKY=equilibrium source-destination temperature
C     HFLWC,HFLWX=source-destination heat flux unltd,ltd by heat
C     ATCNDL=source-destination thermal conductance
C     HFLWSX=source-destination conductive heat flux
C     HFLWL=total conductive+convective source-destination heat flux
C
      IF(VHCP1(N3,N2,N1).GT.VHCPNX(NY,NX))THEN
      IF(N3.EQ.NUM(NY,NX).AND.VHCPW(1,N2,N1).LE.VHCPWX(N2,N1))THEN
      TK1X=TK1(N3,N2,N1)-(HWFLVL-HFLXG)/VHCP1(N3,N2,N1)
      ELSE
      TK1X=TK1(N3,N2,N1)-HWFLVL/VHCP1(N3,N2,N1)
      ENDIF
      ELSE
      TK1X=TK1(N3,N2,N1)
      ENDIF
      IF(VHCP1(N6,N5,N4).GT.ZEROS(NY,NX))THEN
      TKLX=TK1(N6,N5,N4)+HWFLVL/VHCP1(N6,N5,N4)
      ELSE
      TKLX=TK1(N6,N5,N4)
      ENDIF
      TKY=(VHCP1(N3,N2,N1)*TK1X+VHCP1(N6,N5,N4)*TKLX)
     2/(VHCP1(N3,N2,N1)+VHCP1(N6,N5,N4))
      HFLWX=(TK1X-TKY)*VHCP1(N3,N2,N1)*XNPX
      HFLWC=ATCNDL*(TK1X-TKLX)*AREA(N,N3,N2,N1)*XNPH
      IF(HFLWC.GE.0.0)THEN
      HFLWSX=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWSX=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
      HFLWL(N,N6,N5,N4)=HWFLWL+HWFLHL+HFLWSX
C     IF(N3.EQ.1)THEN
C     WRITE(*,8765)'HFLWL',I,J,M,N1,N2,N3,N4,N5,N6,N
C    2,HFLWL(N,N3,N2,N1),HFLWL(N,N6,N5,N4)
C    2,HWFLWL,HWFLQL,HWFLVL,HWFLHL,HFLWC,HFLWX,HFLWSX
C    3,ATCNDL,TK1X,TKLX,TKY,HWFLVL,TK1(N3,N2,N1),TK1(N6,N5,N4)
C    2,TCND1,TCND2,DLYR(N,N3,N2,N1),DLYR(N,N6,N5,N4)
C    4,VHCP1(N3,N2,N1),VHCP1(N6,N5,N4),VOLY(N3,N2,N1)
C    5,VOLY(N6,N5,N4),VOLW1(N3,N2,N1),VOLW1(N6,N5,N4)
C    3,THETPX(N3,N2,N1),THETIX(N3,N2,N1),THETWX(N3,N2,N1)
C    3,THETPX(N6,N5,N4),THETIX(N6,N5,N4),THETWX(N6,N5,N4)
C    3,THETPY(N3,N2,N1),THETPY(N6,N5,N4)
C    4,RYLNA2,XNUSA2,XNUSW2
C    4,STC(N6,N5,N4),DTC(N6,N5,N4)
C    2,WTHET2,TCNDA2,TCNDW2
8765  FORMAT(A8,10I4,60E14.6)
C     ENDIF
C
C     TOTAL WATER, VAPOR AND HEAT FLUXES
C
C     FLW,FLWX,FLWH=total water flux through micropores,macropores
C     HFLW=total heat flux
C     FLWM=water flux used for solute flux calculations in trnsfr.f
C
      FLW(N,N6,N5,N4)=FLW(N,N6,N5,N4)+FLWL(N,N6,N5,N4)
      FLWX(N,N6,N5,N4)=FLWX(N,N6,N5,N4)+FLWLX(N,N6,N5,N4)
      FLWH(N,N6,N5,N4)=FLWH(N,N6,N5,N4)+FLWHL(N,N6,N5,N4)
      HFLW(N,N6,N5,N4)=HFLW(N,N6,N5,N4)+HFLWL(N,N6,N5,N4)
      FLWM(M,N,N6,N5,N4)=FLWL(N,N6,N5,N4)
      IF(N.EQ.3)THEN
C     IF(I.EQ.55)THEN
C     WRITE(*,1115)'FLWL2',I,J,M,N1,N2,N3,N4,N5,N6,N,FLWL(N,N3,N2,N1)
C    2,FLWL(N,N6,N5,N4),FLW(N,N3,N2,N1),FLW(N,N6,N5,N4)
C    3,FLQL,FLVL,FLQX,FLVX,HFLWX
C    3,CND1,CNDL,AVCNDL,AVCNVL,VP1,VPL,PSIST1,PSISTL
C    4,UAG,VOLA1(N6,N5,N4),VOLI1(N6,N5,N4),SCNV(N6,N5,N4)
C    5,VOLP1(N3,N2,N1),VOLP1(N6,N5,N4),TKY
C    7,TK1(N3,N2,N1),TK1(N6,N5,N4),VOLT(N3,N2,N1),VOLT(N6,N5,N4)
C    8,VOLW1(N6,N5,N4),VOLP1(N6,N5,N4),VOLX(N6,N5,N4),VOLW1(N3,N2,N1)
C    9,VOLP1(N3,N2,N1),VOLX(N3,N2,N1),VOLA1(N6,N5,N4),VOLA1(N3,N2,N1)
C    6,THETW1,THETWL,PSISA1(N3,N2,N1)
C    7,PSISA1(N6,N5,N4),PSISM1(N3,N2,N1)
C    7,PSISM1(N6,N5,N4),PSISH(N3,N2,N1),PSISH(N6,N5,N4)
C    8,DLYR(N,N3,N2,N1),DLYR(N,N6,N5,N4),AREA(N,N3,N2,N1)
C    9,VHCP1(N3,N2,N1),VHCP1(N6,N5,N4),POROS(N6,N5,N4)
C     ENDIF
C
C     WATER FILM THICKNESS FOR CALCULATING GAS EXCHANGE IN TRNSFR.F
C
      FILM(M,N6,N5,N4)=AMAX1(1.0E-06
     2,EXP(-13.833-0.857*LOG(-PSISA1(N6,N5,N4))))
      ENDIF
      ELSEIF(N.NE.3)THEN
      FLWL(N,N6,N5,N4)=0.0
      FLWLX(N,N6,N5,N4)=0.0
      FLWHL(N,N6,N5,N4)=0.0
      HFLWL(N,N6,N5,N4)=0.0
      FLWM(M,N,N6,N5,N4)=0.0
      FLWHM(M,N,N6,N5,N4)=0.0
      ENDIF
      ELSE
      IF(N.EQ.3)THEN
      FLWL(N,N3,N2,N1)=0.0
      FLWLX(N,N3,N2,N1)=0.0
      FLWHL(N,N3,N2,N1)=0.0
      HFLWL(N,N3,N2,N1)=0.0
      FLWHM(M,N,N3,N2,N1)=0.0
      FLWHM(M,N,N3,N2,N1)=0.0
      ELSE
      FLWL(N,N6,N5,N4)=0.0
      FLWLX(N,N6,N5,N4)=0.0
      FLWHL(N,N6,N5,N4)=0.0
      HFLWL(N,N6,N5,N4)=0.0
      FLWM(M,N,N6,N5,N4)=0.0
      FLWHM(M,N,N6,N5,N4)=0.0
      ENDIF
C     IF(I.EQ.336)THEN
C     WRITE(*,1115)'FLWLX',I,J,M,N1,N2,N3,N4,N5,N6,N
C    2,FLWL(N,N3,N2,N1),FLW(N,N3,N2,N1)
C    2,FLWL(N,N6,N5,N4),FLW(N,N6,N5,N4)
C    3,VOLX(N3,N2,N1),VOLX(N6,N5,N4)
C     ENDIF
      ENDIF
4320  CONTINUE
4400  CONTINUE
9890  CONTINUE
9895  CONTINUE
C
C     BOUNDARY WATER AND HEAT FLUXES
C
C     XVOLT,XVOLW=excess water+ice,water in source grid cell
C     VOLP2,VOLPH2=air-filled porosity in micropores,macropores
C
      DO 9595 NX=NHW,NHE
      DO 9590 NY=NVN,NVS
      DO 9585 L=NUM(NY,NX),NL(NY,NX)
      VOLP2=VOLA1(L,NY,NX)-VOLW1(L,NY,NX)-VOLI1(L,NY,NX)
      VOLPX2=VOLP2
      VOLPH2=VOLAH1(L,NY,NX)-VOLWH1(L,NY,NX)-VOLIH1(L,NY,NX)
C
C     IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO WATER TABLE
C
C     IDTBL=water table flag
C     DPTH,DTBLX=depth to layer midpoint,natural water table
C     PSISM1,PSISE=matric,air entry water potential
C     DTBLXX=equilibrium water potential with natural water table
C     DPTHA=active layer depth
C     IFLGU=micropore discharge flag to natural water table
C
      IF(IDTBL(NY,NX).NE.0.AND.DPTH(L,NY,NX).LT.DTBLX(NY,NX))THEN
      IF(PSISM1(L,NY,NX).GT.0.0098*(DPTH(L,NY,NX)-DTBLX(NY,NX)))THEN
      IFLGU=0
      DO 9565 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
      DTBLXX=DTBLX(NY,NX)+PSISE(LL,NY,NX)/0.0098
      IF(DPTH(LL,NY,NX).LT.DTBLXX)THEN
      IF((PSISM1(LL,NY,NX).LE.0.0098*(DPTH(LL,NY,NX)-DTBLXX)
     2.AND.L.NE.NL(NY,NX)).OR.DPTH(LL,NY,NX).GT.DPTHA(NY,NX))THEN
      IFLGU=1
      ENDIF
      ENDIF
9565  CONTINUE
      ELSE
      IFLGU=1
      ENDIF
      ELSE
      IFLGU=1
      ENDIF
C
C     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO WATER TABLE
C
C     VOLAH1,VOLWH1,VOLIH1=macropore volume,water,ice content
C     DPTHH depth to layer macropore water
C     CDPTH=depth to layer bottom
C     DLYR=layer thickness
C     IFLGUH=macropore discharge flag to natural water table
C
      IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DPTHH=CDPTH(L,NY,NX)-(VOLWH1(L,NY,NX)+VOLIH1(L,NY,NX))
     2/VOLAH1(L,NY,NX)*DLYR(3,L,NY,NX)
      ELSE
      DPTHH=CDPTH(L,NY,NX)
      ENDIF
      IF(IDTBL(NY,NX).NE.0.AND.DPTHH.LT.DTBLX(NY,NX)
     2.AND.VOLWH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IFLGUH=0
C     DO 9566 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
C     IF(DPTH(LL,NY,NX).LT.DTBLX(NY,NX))THEN
C     IF(VOLAH1(LL,NY,NX).LE.ZEROS(NY,NX))THEN
C     IFLGUH=1
C     ENDIF
C     ENDIF
9566  CONTINUE
      ELSE
      IFLGUH=1
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,9567)'IFLGU',I,J,M,NX,NY,L,IFLGU,IFLGUH,PSISM1(L,NY,NX)
C    2,0.0098*(DPTH(L,NY,NX)-DTBLX(NY,NX)),DTBLXX
C    2,DPTH(L,NY,NX),DTBLX(NY,NX),DTBLZ(NY,NX),DTBLI(NY,NX)
C    3,VOLAH1(L,NY,NX),VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),CDPTH(L,NY,NX)
C    4,DLYR(3,L,NY,NX),DTBLZ(NY,NX),DPTHH,THETX,DPTHA(NY,NX)
9567  FORMAT(A8,8I4,30E12.4)
C     ENDIF
C
C     IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO TILE DRAIN
C
C     IDTBL=water table flag
C     DPTH,DTBLY=depth to layer midpoint, artificial water table
C     PSISM1,PSISE=soil,air entry matric potential
C     DTBLYX=equilibrium water potential with artificial water table
C     IFLGD=micropore discharge flag to artificial water table
C
      IF(IDTBL(NY,NX).GE.3.AND.DPTH(L,NY,NX).LT.DTBLY(NY,NX))THEN
      IF(PSISM1(L,NY,NX).GT.0.0098*(DPTH(L,NY,NX)-DTBLY(NY,NX)))THEN
      IFLGD=0
      IF(L.LT.NL(NY,NX))THEN
      DO 9568 LL=L+1,NL(NY,NX)
      DTBLYX=DTBLY(NY,NX)+PSISE(LL,NY,NX)/0.0098
      IF(DPTH(LL,NY,NX).LT.DTBLYX)THEN
      IF((PSISM1(LL,NY,NX).LE.0.0098*(DPTH(LL,NY,NX)-DTBLYX)
     2.AND.L.NE.NL(NY,NX)).OR.DPTH(LL,NY,NX).GT.DPTHA(NY,NX))THEN
      IFLGD=1
      ENDIF
      ENDIF
9568  CONTINUE
      ENDIF
      ELSE
      IFLGD=1
      ENDIF
      ELSE
      IFLGD=1
      ENDIF
C
C     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO TILE DRAIN
C
C     VOLAH1,VOLWH1,VOLIH1=macropore volume,water,ice content
C     CDPTH=depth to layer bottom
C     DLYR=layer thickness
C     IFLGDH=macropore discharge flag to artificial water table
C
      IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DPTHH=CDPTH(L,NY,NX)-(VOLWH1(L,NY,NX)+VOLIH1(L,NY,NX))
     2/VOLAH1(L,NY,NX)*DLYR(3,L,NY,NX)
      ELSE
      DPTHH=CDPTH(L,NY,NX)
      ENDIF
      IF(IDTBL(NY,NX).GE.3.AND.DPTHH.LT.DTBLY(NY,NX)
     2.AND.VOLWH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IFLGDH=0
      IF(L.LT.NL(NY,NX))THEN
      DO 9569 LL=L+1,NL(NY,NX)
      IF(DPTH(LL,NY,NX).LT.DTBLY(NY,NX))THEN
      IF(VOLAH1(LL,NY,NX).LE.ZEROS(NY,NX))THEN
      IFLGDH=1
      ENDIF
      ENDIF
9569  CONTINUE
      ENDIF
      ELSE
      IFLGDH=1
      ENDIF
C     IF(L.EQ.12)THEN
C     WRITE(*,9567)'IFLGD',I,J,M,NX,NY,L,IFLGD,IFLGDH,PSISM1(L,NY,NX)
C    2,PSISE(L,NY,NX)+0.0098*(DPTH(L,NY,NX)-DTBLY(NY,NX))
C    2,DPTH(L,NY,NX),DTBLY(NY,NX),DTBLD(NY,NX),DTBLDI(NY,NX)
C    3,VOLAH1(L,NY,NX),VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),CDPTH(L,NY,NX)
C    4,DLYR(3,L,NY,NX),DPTHH,THETX,DPTHA(NY,NX),RCHGFT
C     ENDIF
C
C     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
C     ENTERED IN 'READS'
C
C     N3,N2,N1=L,NY,NX of source grid cell
C     M6,M5,M4=L,NY,NX of destination grid cell
C
      N1=NX
      N2=NY
      N3=L
C
C     LOCATE EXTERNAL BOUNDARIES
C
      DO 9580 N=1,3
      DO 9575 NN=1,2
      IF(N.EQ.1)THEN
      N4=NX+1
      N5=NY
      N4B=NX-1
      N5B=NY
      N6=L
      IF(NN.EQ.1)THEN
      IF(NX.EQ.NHE)THEN
      M1=NX
      M2=NY
      M3=L
      M4=NX+1
      M5=NY
      M6=L
      XN=-1.0
      RCHQF=RCHQE(M2,M1)
      RCHGFU=RCHGEU(M2,M1)
      RCHGFT=RCHGET(M2,M1)
      ELSE
      GO TO 9575
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NX.EQ.NHW)THEN
      M1=NX+1
      M2=NY
      M3=L
      M4=NX
      M5=NY
      M6=L
      XN=1.0
      RCHQF=RCHQW(M5,M4)
      RCHGFU=RCHGWU(M5,M4)
      RCHGFT=RCHGWT(M5,M4)
      ELSE
      GO TO 9575
      ENDIF
      ENDIF
      ELSEIF(N.EQ.2)THEN
      N4=NX
      N5=NY+1
      N4B=NX
      N5B=NY-1
      N6=L
      IF(NN.EQ.1)THEN
      IF(NY.EQ.NVS)THEN
      M1=NX
      M2=NY
      M3=L
      M4=NX
      M5=NY+1
      M6=L
      XN=-1.0
      RCHQF=RCHQS(M2,M1)
      RCHGFU=RCHGSU(M2,M1)
      RCHGFT=RCHGST(M2,M1)
      ELSE
      GO TO 9575
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NY.EQ.NVN)THEN
      M1=NX
      M2=NY+1
      M3=L
      M4=NX
      M5=NY
      M6=L
      XN=1.0
      RCHQF=RCHQN(M5,M4)
      RCHGFU=RCHGNU(M5,M4)
      RCHGFT=RCHGNT(M5,M4)
      ELSE
      GO TO 9575
      ENDIF
      ENDIF
      ELSEIF(N.EQ.3)THEN
      N4=NX
      N5=NY
      N6=L+1
      IF(NN.EQ.1)THEN
      IF(L.EQ.NL(NY,NX))THEN
      M1=NX
      M2=NY
      M3=L
      M4=NX
      M5=NY
      M6=L+1
      XN=-1.0
      RCHGFU=RCHGD(M2,M1)
      RCHGFT=1.0
      ELSE
      GO TO 9575
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      GO TO 9575
      ENDIF
      ENDIF
C
C     BOUNDARY SURFACE RUNOFF DEPENDING ON ASPECT, SLOPE
C     VELOCITY, HYDRAULIC RADIUS AND SURFACE WATER STORAGE
C
C     CDPTH,CDPTHI=current,initial surface elevation
C     BKDS=bulk density
C     IRCHG,RCHQ*=runoff boundary flags
C
      IF(L.EQ.NUM(N2,N1).AND.N.NE.3
     2.AND.(CDPTH(NU(N2,N1)-1,N2,N1).LE.CDPTHI(N2,N1)
     3.OR.BKDS(NUI(N2,N1),N2,N1).GT.ZERO))THEN
C     WRITE(*,7744)'QRCHK',I,J,M,N1,N2,N4,N5,M4,M5,N,NN
C    1,IRCHG(NN,N,N2,N1),RCHQF,QRM(M,N2,N1)
      IF(IRCHG(NN,N,N2,N1).EQ.0.OR.RCHQF.EQ.0.0
     2.OR.ABS(QRM(M,N2,N1)).LT.ZEROS(N2,N1))THEN
      QR1(N,NN,M5,M4)=0.0
      HQR1(N,NN,M5,M4)=0.0
      ELSE
C
C     SURFACE BOUNDARY WATER FLUX
C
C     DPTHW1,DPTHW2=surface water depth of source,destination
C     ALT1,ALT2=elevation of source,destination
C     XVOLT=excess surface water+ice
C     VOLWG=ground surface water retention capacity
C     DTBLX=natural water table depth
C     QR1,HQR1=runoff, convective heat from runoff
C     QR,HQR=hourly-accumulated runoff, convective heat from runoff
C     QRM,QRV=runoff,velocity for erosion, solute transfer
C     XN=direction
C
C     RUNOFF
C
      DPTHW1=XVOLT(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
      DPTHW2=VOLWG(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
      ALT1=ALTG(N2,N1)+DPTHW1
      ALT2=ALTG(N2,N1)+DPTHW2-XN*SLOPE(N,N2,N1)*DLYR(N,NUM(N2,N1),N2,N1)
      IF(ALT1.GT.ALT2
     2.AND.CDPTH(NU(N2,N1)-1,N2,N1)-DPTHW1.LT.DTBLX(N2,N1))THEN
      QR1(N,NN,M5,M4)=-XN*QRM(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
      HQR1(N,NN,M5,M4)=4.19*TK1(0,N2,N1)*QR1(N,NN,M5,M4)
      QR(N,NN,M5,M4)=QR(N,NN,M5,M4)+QR1(N,NN,M5,M4)
      HQR(N,NN,M5,M4)=HQR(N,NN,M5,M4)+HQR1(N,NN,M5,M4)
C     IF(N1.EQ.1.AND.N2.EQ.1)THEN
C     WRITE(*,7744)'QRBND',I,J,M,N1,N2,N4,N5,M4,M5,N,NN
C    1,IRCHG(NN,N,N2,N1),QRM(M,N2,N1)
C    2,QR1(N,NN,M5,M4),QR(N,NN,M5,M4)
C    2,ALT1,ALT2,ALTG(N2,N1),FSLOPE(N,N2,N1)
C    3,VOLWG(N2,N1),VOLWRX(N2,N1),ZM(N2,N1),ZS(N2,N1)
C    4,VOLW1(0,N2,N1),VOLI1(0,N2,N1),DLYR(N,NUM(N2,N1),N2,N1)
C    5,XVOLT(N2,N1)-VOLWG(N2,N1),DTBLX(N2,N1)
C    6,DTBLX(N2,N1)-CDPTH(NU(N2,N1)-1,N2,N1)+DPTHW1
C    7,XVOLTM(M,N2,N1),XVOLWM(M,N2,N1)
7744  FORMAT(A8,12I4,40E12.4)
C     ENDIF
C
C     RUNON
C
      ELSEIF(CDPTH(NU(N2,N1)-1,N2,N1)-DPTHW1.GT.DTBLX(N2,N1))THEN
      VX=AMIN1(0.0,(DTBLX(N2,N1)-CDPTH(NU(N2,N1)-1,N2,N1)+DPTHW1)
     2*AREA(3,NUM(N2,N1),N2,N1))
      QRM(M,N2,N1)=VX*XNPX
      QRV(M,N2,N1)=0.0
      QR1(N,NN,M5,M4)=-XN*QRM(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
      HQR1(N,NN,M5,M4)=4.19*TK1(0,N2,N1)*QR1(N,NN,M5,M4)
      QR(N,NN,M5,M4)=QR(N,NN,M5,M4)+QR1(N,NN,M5,M4)
      HQR(N,NN,M5,M4)=HQR(N,NN,M5,M4)+HQR1(N,NN,M5,M4)
C     WRITE(*,7744)'QRBNB',I,J,M,N1,N2,N4,N5,M4,M5,N,NN
C    1,IRCHG(NN,N,N2,N1),QRM(M,N2,N1)
C    2,QR1(N,NN,M5,M4),QR(N,NN,M5,M4),QRMN(M,N,NN,M5,M4)
C    2,ALTG(N2,N1),FSLOPE(N,N2,N1)
C    3,VOLWG(N2,N1),VOLWRX(N2,N1),ZM(N2,N1),ZS(N2,N1)
C    4,VOLW1(0,N2,N1),VOLI1(0,N2,N1),DLYR(N,NUM(N2,N1),N2,N1)
C    5,XVOLT(N2,N1)-VOLWG(N2,N1),DTBLX(N2,N1)
C    6,CDPTH(NU(N2,N1)-1,N2,N1),DPTHW1
C    7,XVOLTM(M,N2,N1),XVOLWM(M,N2,N1)
      ELSE
      QR1(N,NN,M5,M4)=0.0
      HQR1(N,NN,M5,M4)=0.0
      ENDIF
      QRMN(M,N,NN,M5,M4)=QR1(N,NN,M5,M4)
      IFLBM(M,N,NN,M5,M4)=0
C     WRITE(*,7745)'QRB',I,J,M,N1,N2,N4,N5,M4,M5,N,NN
C    2,IFLBM(M,N,NN,M5,M4)
C    2,QRM(M,N2,N1),QRMN(M,N,NN,M5,M4),QR(N,NN,M5,M4)
7745  FORMAT(A8,12I4,40E12.4)
C
C     BOUNDARY SNOW FLUX
C
C     QS1,QW1,QI1=snow,water,ice transfer
C     HQS1=convective heat transfer from snow,water,ice transfer
C     QS,QW,QI=cumulative hourly snow,water,ice transfer
C     HQS=cumulative hourly convective heat transfer from snow,water,ice transfer
C
      IF(NN.EQ.1)THEN
      QS1(N,M5,M4)=0.0
      QW1(N,M5,M4)=0.0
      QI1(N,M5,M4)=0.0
      HQS1(N,M5,M4)=0.0
C     QS(N,M5,M4)=QS(N,M5,M4)+QS1(N,M5,M4)
C     QW(N,M5,M4)=QW(N,M5,M4)+QW1(N,M5,M4)
C     QI(N,M5,M4)=QI(N,M5,M4)+QI1(N,M5,M4)
C     HQS(N,M5,M4)=HQS(N,M5,M4)+HQS1(N,M5,M4)
      QSM(M,N,M5,M4)=QS1(N,M5,M4)
      ENDIF
      ENDIF
      ELSE
      IF(N.NE.3)THEN
      QR1(N,NN,M5,M4)=0.0
      HQR1(N,NN,M5,M4)=0.0
      ENDIF
      ENDIF
C
C     BOUNDARY SUBSURFACE WATER AND HEAT TRANSFER DEPENDING
C     ON LEVEL OF WATER TABLE
C
      IF(VOLX(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
      IF(NCN(N2,N1).NE.3.OR.N.EQ.3)THEN
C
C     IF NO WATER TABLE
C
C     IDTBL=water table flag
C     THETA1,THETAX=water content ahead,behind wetting front
C     K1,KL=pore water class ahead,behind wetting front
C     CND1,CNDL=hydraulic conductivity ahead,behind wetting front
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
C     FLWL,FLWLX=lower boundary micropore water flux
C     FLWHL=lower boundary macropore water flux
C     HFLWL=convective heat from lower boundary water flux
C     XH,XN,XNPH=rate constant,direction indicator,time step
C     SLOPE=sin(vertical slope)=1
C     RCHG*=boundary flags
C
      IF(IDTBL(N2,N1).EQ.0.OR.N.EQ.3)THEN
      THETA1=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1)
     2,VOLW1(N3,N2,N1)/VOLY(N3,N2,N1)))
      THETAX=AMAX1(THETY(N3,N2,N1),AMIN1(POROS(N3,N2,N1)
     2,VOLWX1(N3,N2,N1)/VOLY(N3,N2,N1)))
      K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)
     2-THETA1)/POROS(N3,N2,N1))+1))
      KL=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)
     2-THETAX)/POROS(N3,N2,N1))+1))
      IF(N3.EQ.NUM(NY,NX))THEN
      CND1=HCND(N,K1,N3,N2,N1)*FKSAT
      ELSE
      CND1=HCND(N,K1,N3,N2,N1)
      ENDIF
      CNDL=HCND(N,KL,N3,N2,N1)
      FLWL(N,M6,M5,M4)=AMIN1(VOLW1(N3,N2,N1)*XNPX
     2,XN*0.0098*-ABS(SLOPE(N,N2,N1))*CND1*AREA(3,N3,N2,N1))
     3*RCHGFU*RCHGFT*XNPH
      FLWLX(N,M6,M5,M4)=FLWL(N,M6,M5,M4)
      FLWHL(N,M6,M5,M4)=AMIN1(VOLWH1(L,NY,NX)*XNPX
     2,XN*0.0098*-ABS(SLOPE(N,N2,N1))*CNDH1(L,NY,NX)*AREA(3,N3,N2,N1))
     3*RCHGFU*RCHGFT*XNPH
      HFLWL(N,M6,M5,M4)=4.19*TK1(N3,N2,N1)
     2*(FLWL(N,M6,M5,M4)+FLWHL(N,M6,M5,M4))
C     IF(I.EQ.336)THEN
C     WRITE(*,4443)'ABV',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,XN
C    2,FLWL(N,M6,M5,M4)
C    2,VOLP2,RCHGFU,VOLX(N3,N2,N1),VOLW1(N3,N2,N1)
C    3,VOLWH1(N3,N2,N1),VOLPH1(N3,N2,N1),VOLPH2,VOLI1(N3,N2,N1)
C    4,VOLIH1(N3,N2,N1),VOLP1(N3,N2,N1),HFLWL(N,M6,M5,M4)
C    5,PSISM1(N3,N2,N1),PSISE(N3,N2,N1),FLWHL(N,M6,M5,M4),DTBLD(N2,N1)
C    6,SLOPE(N,N2,N1)
4443  FORMAT(A8,11I4,30E12.4)
C     ENDIF
      ELSE
C
C     MICROPORE DISCHARGE ABOVE WATER TABLE
C
C     IFLGU=micropore discharge flag to natural water table
C     PSISWD=water potential from water table slope
C     XN,RCHG*=direction indicator,boundary flag
C     SLOPE=sin(lateral slope)
C     DLYR=layer width
C     DTBLG=water table slope
C     PSISWT=water potential driving micropore discharge
C     PSISA1,PSISO=matric,osmotic water potential
C     DPTH,DTBLX=depth to layer midpoint,natural water table
C     DPTHT=depth to internal water table
C     FLWL=micropore discharge to natural water table
C     HFLWL=convective heat from discharge to natural water table
C     HCND=saturated hydraulic conductivity
C     AREAU=fraction of layer below natural water table
C
      IF(IFLGU.EQ.0.AND.RCHGFT.NE.0.0)THEN
      PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISWT=AMIN1(0.0,-PSISA1(N3,N2,N1)-0.03*PSISO(N3,N2,N1)
     2+0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1))
     3-0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1)))
      IF(PSISWT.LT.0.0)PSISWT=PSISWT-PSISWD
      FLWT=PSISWT*HCND(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)
     2*(1.0-AREAU(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*XNPH
      FLWL(N,M6,M5,M4)=XN*FLWT
      FLWLX(N,M6,M5,M4)=XN*FLWT
      HFLWL(N,M6,M5,M4)=4.19*TK1(N3,N2,N1)*XN*FLWT
C     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,4445)'DISCHMI',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGU
C    2,XN,FLWL(N,M6,M5,M4),FLWT,PSISWT,PSISWD,HCND(N,1,N3,N2,N1)
C    3,AREA(N,N3,N2,N1),AREAU(N3,N2,N1),RCHGFU,RCHGFT
C    4,PSISE(N3,N2,N1),PSISA1(N3,N2,N1),DPTH(N3,N2,N1),DTBLX(N2,N1)
C    2,PSISO(N3,N2,N1),0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1))
C    3,0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1))
4445  FORMAT(A8,12I4,30E14.6)
C     ENDIF
      ELSE
      FLWL(N,M6,M5,M4)=0.0
      FLWLX(N,M6,M5,M4)=0.0
      HFLWL(N,M6,M5,M4)=0.0
      ENDIF
C
C     MACROPORE DISCHARGE ABOVE WATER TABLE
C
C     IFLGUH=macropore discharge flag to natural water table
C     PSISWD=water potential from water table slope
C     XN,RCHG*=direction indicator,boundary flag
C     SLOPE=sin(lateral slope)
C     DLYR=layer width
C     DTBLG=water table slope
C     PSISWTH=water potential driving macropore discharge
C     PSISO=osmotic water potential
C     DPTHH,DTBLX=depth to layer macropore water,natural water table
C     DPTHT=depth to internal water table
C     FLWTH,FLWTHL=macropore discharge unltd,ltd by macropore water
C     CNDH1=macropore hydraulic conductivity
C     FLWHL=macropore discharge to natural water table
C     HFLWL=convective heat from discharge to natural water table
C     HCND=saturated hydraulic conductivity
C     AREAU=fraction of layer below natural water table
C
      IF(IFLGUH.EQ.0.AND.RCHGFT.NE.0.0
     2.AND.VOLAH1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISWTH=-0.03*PSISO(N3,N2,N1)
     2+0.0098*(DPTHH-DTBLX(N2,N1))
     3-0.0098*AMAX1(0.0,DPTHH-DPTHT(N2,N1))
      IF(PSISWTH.LT.0.0)PSISWTH=PSISWTH-PSISWD
      FLWTH=PSISWTH*CNDH1(N3,N2,N1)*AREA(N,N3,N2,N1)
     2*(1.0-AREAU(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*XNPH
      FLWTHL=AMAX1(FLWTH,AMIN1(0.0,-(VOLWH1(N3,N2,N1)*XNPX
     2+FLWHL(3,N3,N2,N1)-FLWHL(3,N3+1,N2,N1))))
      FLWHL(N,M6,M5,M4)=XN*FLWTHL
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)
     2+4.19*TK1(N3,N2,N1)*XN*FLWTHL
C     WRITE(*,4446)'DISCHMA',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGUH
C    2,XN,FLWHL(N,M6,M5,M4),FLWTHL,FLWTH,PSISWTH,CNDH1(N3,N2,N1)
C    3,DPTH(N3,N2,N1),DLYR(3,N3,N2,N1),DPTHH,VOLWH1(N3,N2,N1)
C    4,VOLIH1(L,NY,NX),VOLAH1(N3,N2,N1),DTBLX(N2,N1),PSISWD
4446  FORMAT(A8,12I4,30E14.6)
      ELSE
      FLWHL(N,M6,M5,M4)=0.0
      ENDIF
C
C     MICROPORE DISCHARGE ABOVE TILE DRAIN
C
C     IFLGD=micropore discharge flag to artificial water table
C     PSISWD=water potential from water table slope
C     XN,RCHG*=direction indicator,boundary flag
C     SLOPE=sin(lateral slope)
C     DLYR=layer width
C     DTBLG=water table slope
C     PSISWT=water potential driving micropore discharge
C     PSISA1,PSISO=matric,osmotic water potential
C     DPTH,DTBLY=depth to layer midpoint,artificial water table
C     DPTHT=depth to internal water table
C     FLWL=micropore discharge to natural+artificial water table
C     HFLWL=convective heat from dischg to natural+artifl water table
C     HCND=saturated hydraulic conductivity
C     AREAUD=fraction of layer below artificial water table
C
      IF(IFLGD.EQ.0.AND.RCHGFT.NE.0.0)THEN
      PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISWT=AMIN1(0.0,-PSISA1(N3,N2,N1)-0.03*PSISO(N3,N2,N1)
     2+0.0098*(DPTH(N3,N2,N1)-DTBLY(N2,N1))
     3-0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1)))
      IF(PSISWT.LT.0.0)PSISWT=PSISWT-PSISWD
      FLWT=PSISWT*HCND(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)
     2*(1.0-AREAUD(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*XNPH
      FLWL(N,M6,M5,M4)=FLWL(N,M6,M5,M4)+XN*FLWT
      FLWLX(N,M6,M5,M4)=FLWLX(N,M6,M5,M4)+XN*FLWT
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+4.19*TK1(N3,N2,N1)*XN*FLWT
C     IF(L.EQ.12)THEN
C     WRITE(*,4445)'DISCHMD',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGD
C    2,XN,FLWL(N,M6,M5,M4),FLWT,PSISWT,PSISWD
C    3,0.0098*(DPTH(N3,N2,N1)-DTBLY(N2,N1)),RCHGFU,RCHGFT
C    3,HCND(N,1,N3,N2,N1),AREA(N,N3,N2,N1),AREAUD(N3,N2,N1)
C    4,PSISE(N3,N2,N1),PSISA1(N3,N2,N1),PSISO(N3,N2,N1)
C    5,DPTH(N3,N2,N1),DTBLY(N2,N1),DPTHT(N2,N1)
C    3,0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1))
C     ENDIF
      ENDIF
C
C     MACROPORE DISCHARGE ABOVE TILE DRAIN
C
C     IFLGDH=macropore discharge flag to artificial water table
C     PSISWD=water potential from water table slope
C     XN,RCHG*=direction indicator,boundary flag
C     SLOPE=sin(lateral slope)
C     DLYR=layer width
C     DTBLG=water table slope
C     PSISWTH=water potential driving macropore discharge
C     PSISO=osmotic water potential
C     DPTHH,DTBLY=depth to layer macropore water,artificl water table
C     DPTHT=depth to internal water table
C     FLWTH,FLWTHL=macropore discharge unltd,ltd by macropore water
C     CNDH1=macropore hydraulic conductivity
C     FLWHL=macropore discharge to artificial water table
C     HFLWL=convective heat from discharge to artificial water table
C     HCND=saturated hydraulic conductivity
C     AREAUD=fraction of layer below artificial water table
C
      IF(IFLGDH.EQ.0.AND.RCHGFT.NE.0.0
     2.AND.VOLAH1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISWTH=-0.03*PSISO(N3,N2,N1)
     2+0.0098*(DPTHH-DTBLY(N2,N1))
     3-0.0098*AMAX1(0.0,DPTHH-DPTHT(N2,N1))
      IF(PSISWTH.LT.0.0)PSISWTH=PSISWTH-PSISWD
      FLWTH=PSISWTH*CNDH1(N3,N2,N1)*AREA(N,N3,N2,N1)
     2*(1.0-AREAUD(N3,N2,N1))/(RCHGFU+1.0)*RCHGFT*XNPH
      FLWTHL=AMAX1(FLWTH,AMIN1(0.0
     2,-(VOLWH1(N3,N2,N1)*XNPX+FLWHL(3,N3,N2,N1)-FLWHL(3,N3+1,N2,N1))))
      FLWHL(N,M6,M5,M4)=FLWHL(N,M6,M5,M4)+XN*FLWTHL
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)
     2+4.19*TK1(N3,N2,N1)*XN*FLWTHL
C     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,4446)'DISCHDH',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGDH
C    2,XN,FLWHL(N,M6,M5,M4),FLWTHL,FLWTH,PSISWTH,CNDH1(N3,N2,N1)
C    3,DPTH(N3,N2,N1),DLYR(3,N3,N2,N1),DPTHH,VOLWH1(N3,N2,N1)
C    4,VOLIH1(L,NY,NX),VOLAH1(N3,N2,N1),DTBLY(N2,N1),PSISWD
C    2,0.0098*(DPTHH-DTBLY(N2,N1))
C    3,0.0098*AMAX1(0.0,DPTHH-DPTHT(N2,N1))
C     ENDIF
      ENDIF
C
C     MICROPORE RECHARGE BELOW WATER TABLE
C
C     DPTH,DTBLX=depth to layer midpoint,natural water table
C     DPTHA=active layer depth
C     VOLP2=air volume
C     PSISWD=water potential from water table slope
C     XN,RCHG*=direction indicator,boundary flag
C     SLOPE=sin(lateral slope)
C     DLYR=layer width
C     DTBLG=water table slope
C     PSISUT=water potential driving micropore recharge
C     PSISA1,PSISO=matric,osmotic water potential
C     DPTH,DTBLX=depth to layer midpoint,natural water table
C     FLWU,FLWUL=micropore recharge unltd,ltd by micropore air volume
C     FLWL=micropore recharge from natural water table
C     HFLWL=convective heat from recharge from natural water table
C     HCND=saturated hydraulic conductivity
C     AREAU=fraction of layer below natural water table
C
      IF(DPTH(N3,N2,N1).GE.DTBLX(N2,N1)
     2.AND.DPTHA(N2,N1).GT.DTBLX(N2,N1)
     3.AND.DPTH(N3,N2,N1).LT.DPTHA(N2,N1)
     4.AND.(VOLP2.GT.ZEROS2(N2,N1).OR.BKDS(N3,N2,N1).LE.ZERO)
     4.AND.VOLP1Z(N3,N2,N1).GT.0.0
     5.AND.RCHGFT.NE.0.0)THEN
      PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISUT=AMAX1(0.0,-PSISA1(N3,N2,N1)-0.03*PSISO(N3,N2,N1)
     2+0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1)))
      IF(PSISUT.GT.0.0)PSISUT=PSISUT+PSISWD
      FLWU=PSISUT*HCND(N,1,N3,N2,N1)*AREA(N,N3,N2,N1)
     2*AREAU(N3,N2,N1)/(RCHGFU+1.0)*RCHGFT*XNPH
      IF(BKDS(N3,N2,N1).GT.ZERO)THEN
      FLWUL=AMIN1(FLWU,VOLP2)
      FLWUX=AMIN1(FLWU,VOLPX2)
      ELSE
      FLWUL=FLWU
      FLWUX=FLWU
      ENDIF
      FLWL(N,M6,M5,M4)=FLWL(N,M6,M5,M4)+XN*FLWUL
      FLWLX(N,M6,M5,M4)=FLWLX(N,M6,M5,M4)+XN*FLWUX
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+4.19*TK1(N3,N2,N1)
     2*XN*FLWUL
C     IF(J.EQ.15)THEN
C     WRITE(*,4444)'RECHGMI',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGU
C    2,XN,FLWL(N,M6,M5,M4),FLWU,FLWUL
C    2,PSISUT,PSISA1(N3,N2,N1),PSISO(N3,N2,N1)
C    3,0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1)),RCHGFU,RCHGFT
C    4,HCND(N,1,N3,N2,N1),AREA(N,N3,N2,N1),AREAU(N3,N2,N1)
C    5,VOLP2,VOLA1(L,NY,NX),VOLW1(L,NY,NX),VOLI1(L,NY,NX),DTBLX(N2,N1)
C    6,DPTH(N3,N2,N1),DPTHA(N2,N1),DPTHT(N2,N1)
C    6,DTBLD(N2,N1),VOLW1(N3,N2,N1),VOLI1(N3,N2,N1)
C    7,VOLX(N3,N2,N1),VOLP1(N3,N2,N1),VOLP1Z(N3,N2,N1)
C    8,FLWL(3,N3,N2,N1),FLWL(3,N6,N5,N4),QR1(N,NN,M5,M4)
C    9,DLYR(N,N3,N2,N1),DLYR(3,N3,N2,N1),PSISWD
C    1,CDPTH(N3,N2,N1),AREA(N,N3,N2,N1),SLOPE(N,N2,N1)
4444  FORMAT(A8,12I4,40E14.6)
C     ENDIF
      VOLP2=VOLP2-XN*FLWL(N,M6,M5,M4)
      VOLPX2=VOLPX2-XN*FLWLX(N,M6,M5,M4)
      ENDIF
C
C     MACROPORE RECHARGE BELOW WATER TABLE
C
C     PSISWD=water potential from water table slope
C     XN,RCHG*=direction indicator,boundary flag
C     SLOPE=sin(lateral slope)
C     DLYR=layer width
C     DTBLG=water table slope
C     PSISUTH=water potential driving macropore recharge
C     PSISO=osmotic water potential
C     DPTHH,DTBLX=depth to layer macropore water,natural water table
C     DPTHT=depth to internal water table
C     CNDH1=macropore hydraulic conductivity
C     FLWUH,FLWUHL=macropore recharge unltd,ltd by macropore air volume
C     FLWHL=macropore discharge to natural water table
C     HFLWL=convective heat from discharge to natural water table
C     HCND=saturated hydraulic conductivity
C     AREAU=fraction of layer below natural water table
C
      IF(DPTHH.GT.DTBLX(N2,N1)
     2.AND.DPTHA(N2,N1).GT.DTBLX(N2,N1)
     2.AND.DPTH(N3,N2,N1).LT.DPTHA(N2,N1)
     2.AND.VOLPH2.GT.ZEROS2(NY,NX)
     2.AND.RCHGFT.NE.0.0)THEN
      PSISWD=XN*0.005*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISUTH=-0.03*PSISO(N3,N2,N1)
     2+0.0098*(DPTHH-DTBLX(N2,N1))
      IF(PSISUTH.GT.0.0)PSISUTH=PSISUTH+PSISWD
      FLWUH=PSISUTH*CNDH1(N3,N2,N1)*AREA(N,N3,N2,N1)
     2*AREAU(N3,N2,N1)/(RCHGFU+1.0)*RCHGFT*XNPH
      FLWUHL=AMIN1(FLWUH,VOLPH2*XNPX)
      FLWHL(N,M6,M5,M4)=FLWHL(N,M6,M5,M4)+XN*FLWUHL
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+4.19*TK1(N3,N2,N1)
     2*XN*FLWUHL
C     IF(M6.EQ.7)THEN
C     WRITE(*,4447)'RECHGMA',I,J,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGU
C    2,XN,FLWHL(N,M6,M5,M4),FLWUH,FLWUHL,DPTHH,PSISUTH,VOLPH2
C    3,CNDH1(N3,N2,N1),DTBLX(N2,N1),CDPTH(N3,N2,N1),DPTHT(N2,N1)
C    6,DTBLD(N2,N1),DPTH(N3,N2,N1),VOLWH1(N3,N2,N1),VOLPH1(N3,N2,N1)
C    8,FLWHL(3,N3,N2,N1),FLWHL(3,N3+1,N2,N1),RCHGFU,AREA(N,N3,N2,N1)
C    9,DLYR(N,N3,N2,N1),DLYR(3,N3,N2,N1),PSISWD
C    1,SLOPE(N,N2,N1),AREAU(N3,N2,N1)
4447  FORMAT(A8,12I4,40E14.6)
C     ENDIF
      VOLPH2=VOLPH2-XN*FLWHL(N,M6,M5,M4)
      ENDIF
      ENDIF
C
C     SUBSURFACE HEAT SOURCE/SINK
C
C     HFLWL=heat flux across lower boundary
C     TK1=lower boundary soil temperature
C     TKSD=deep source/sink temperature from geothermal flux
C     TCNDG=thermal conductivity below lower boundary
C     DPTHSK,CDPTH=depth of thermal sink/source, lower boundary
C
      IF(N.EQ.3.AND.IETYP(N2,N1).NE.-2)THEN
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+(TK1(N3,N2,N1)
     2-TKSD(N2,N1))*TCNDG/(DPTHSK(N2,N1)-CDPTH(N3,N2,N1))
     3*AREA(N,N3,N2,N1)*XNPH
      ENDIF
      FLW(N,M6,M5,M4)=FLW(N,M6,M5,M4)+FLWL(N,M6,M5,M4)
      FLWX(N,M6,M5,M4)=FLWX(N,M6,M5,M4)+FLWLX(N,M6,M5,M4)
      FLWH(N,M6,M5,M4)=FLWH(N,M6,M5,M4)+FLWHL(N,M6,M5,M4)
      HFLW(N,M6,M5,M4)=HFLW(N,M6,M5,M4)+HFLWL(N,M6,M5,M4)
      FLWM(M,N,M6,M5,M4)=FLWL(N,M6,M5,M4)
      FLWHM(M,N,M6,M5,M4)=FLWHL(N,M6,M5,M4)
      ENDIF
      ELSE
      FLWL(N,M6,M5,M4)=0.0
      FLWLX(N,M6,M5,M4)=0.0
      FLWHL(N,M6,M5,M4)=0.0
      HFLWL(N,M6,M5,M4)=0.0
      FLWM(M,N,M6,M5,M4)=0.0
      FLWHM(M,N,M6,M5,M4)=0.0
      ENDIF
9575  CONTINUE
C
C     NET WATER AND HEAT FLUXES IN RUNOFF AND SNOW DRIFT
C
C     TQR1,THQR1=net runoff,convective heat from runoff
C     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from snowpack runoff
C     QR1,HQR1=runoff, convective heat from runoff
C     QS1,QW1,QI1=snow,water,ice transfer
C     HQS1=convective heat transfer from snow,water,ice transfer
C
      IF(L.EQ.NUM(N2,N1).AND.N.NE.3)THEN
      DO 1202 NN=1,2
      TQR1(N2,N1)=TQR1(N2,N1)+QR1(N,NN,N2,N1)
      THQR1(N2,N1)=THQR1(N2,N1)+HQR1(N,NN,N2,N1)
      IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      TQR1(N2,N1)=TQR1(N2,N1)-QR1(N,NN,N5,N4)
      THQR1(N2,N1)=THQR1(N2,N1)-HQR1(N,NN,N5,N4)
      ENDIF
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,6631)'TQR1',I,J,M,N1,N2,N4,N5,N,NN
C    2,IFLBM(M,N,NN,N5,N4),TQR1(N2,N1),THQR1(N2,N1)
C    2,QR1(N,NN,N2,N1),QR1(N,NN,N5,N4)
C    3,QR(N,NN,N2,N1),QR(N,NN,N5,N4)
C    2,HQR1(N,NN,N2,N1),HQR1(N,NN,N5,N4)
C    3,HQR(N,NN,N2,N1),HQR(N,NN,N5,N4)
6631  FORMAT(A8,10I4,12E12.4)
C     ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      TQR1(N2,N1)=TQR1(N2,N1)-QR1(N,NN,N5B,N4B)
      THQR1(N2,N1)=THQR1(N2,N1)-HQR1(N,NN,N5B,N4B)
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,6631)'TQRB1',I,J,M,N1,N2,N4B,N5B,N,NN
C    2,IFLBM(M,N,NN,N5B,N4B),TQR1(N2,N1),THQR1(N2,N1)
C    2,QR1(N,NN,N5B,N4B),HQR1(N,NN,N5B,N4B)
C    2,QR(N,NN,N5B,N4B),HQR(N,NN,N5B,N4B)
C     ENDIF
      ENDIF
      IF(M.EQ.NPH)THEN
      IFLBH(N,NN,N5,N4)=IFLBM(M,N,NN,N5,N4)
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      IFLBH(N,NN,N5B,N4B)=IFLBM(M,N,NN,N5B,N4B)
      ENDIF
      ENDIF
1202  CONTINUE
      TQS1(N2,N1)=TQS1(N2,N1)+QS1(N,N2,N1)-QS1(N,N5,N4)
      TQW1(N2,N1)=TQW1(N2,N1)+QW1(N,N2,N1)-QW1(N,N5,N4)
      TQI1(N2,N1)=TQI1(N2,N1)+QI1(N,N2,N1)-QI1(N,N5,N4)
      THQS1(N2,N1)=THQS1(N2,N1)+HQS1(N,N2,N1)-HQS1(N,N5,N4)
      ENDIF
C
C     NET WATER AND HEAT FLUXES THROUGH SOIL AND SNOWPACK
C
C     TFLWL,THFLWL=net water micropore,macropore flux
C     THFLWL=net convective+conductive heat flux
C     FLWL =micropore water,heat flux
C     FLWHL=macropore water,heat flux
C     HFLWL=soil heat flux
C
      IF(NCN(N2,N1).NE.3.OR.N.EQ.3)THEN
      DO 1200 LL=N6,NL(N5,N4)
      IF(VOLX(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
      N6=LL
      GO TO 1201
      ENDIF
1200  CONTINUE
1201  CONTINUE
      IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      TFLWL(N3,N2,N1)=TFLWL(N3,N2,N1)+FLWL(N,N3,N2,N1)
     2-FLWL(N,N6,N5,N4)
      TFLWLX(N3,N2,N1)=TFLWLX(N3,N2,N1)+FLWLX(N,N3,N2,N1)
     2-FLWLX(N,N6,N5,N4)
      TFLWHL(N3,N2,N1)=TFLWHL(N3,N2,N1)+FLWHL(N,N3,N2,N1)
     2-FLWHL(N,N6,N5,N4)
      THFLWL(N3,N2,N1)=THFLWL(N3,N2,N1)+HFLWL(N,N3,N2,N1)
     2-HFLWL(N,N6,N5,N4)
C     IF(N3.EQ.1)THEN
C     WRITE(*,3378)'THFLW',I,J,M,N1,N2,N3,N4,N5,N6,N
C    2,TFLWL(N3,N2,N1),FLWL(N,N3,N2,N1),FLWL(N,N6,N5,N4)
C    2,THFLWL(N3,N2,N1),HFLWL(N,N3,N2,N1),HFLWL(N,N6,N5,N4)
C    2,FLW(N,N3,N2,N1),FLW(N,N6,N5,N4)
C    2,HFLW(N,N3,N2,N1),HFLW(N,N6,N5,N4)
3378  FORMAT(A8,10I4,20E14.6)
C     ENDIF
      ELSE
      TFLWL(N3,N2,N1)=0.0
      TFLWLX(N3,N2,N1)=0.0
      TFLWHL(N3,N2,N1)=0.0
      THFLWL(N3,N2,N1)=0.0
      ENDIF
      ENDIF
9580  CONTINUE
C
C     INFILTRATION OF WATER FROM MACROPORES INTO MICROPORES
C
C     VOLWH1=macropore volume
C     FINHX,FINHL=macro-micropore transfer unltd,ltd by water,air volume
C     FINHM=macro-micropore transfer for use in trnsfr.f
C     HCND=saturated hydraulic conductivity
C     PSISE,PSISA1=air entry,matric water potentials
C     PHOL,HRAD=path length between,radius of macropores from hour1.f
C     XNPH=time step
C     VOLW1X,VOLP1X=current micropore water,air volume
C     VOLWH1X,VOLPH1X=current macropore water,air volume
C
      IF(VOLWH1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      FINHX=6.283*HCND(2,1,N3,N2,N1)*AREA(3,N3,N2,N1)
     2*(PSISE(N3,N2,N1)-PSISA1(N3,N2,N1))
     3/LOG(PHOL(N3,N2,N1)/HRAD(N3,N2,N1))*XNPH
      VOLW1X=VOLW1(N3,N2,N1)+TFLWL(N3,N2,N1)+FLU1(N3,N2,N1)
      VOLP1X=AMAX1(0.0,VOLA1(N3,N2,N1)-VOLW1X-VOLI1(N3,N2,N1))
      VOLWH1X=VOLWH1(N3,N2,N1)+TFLWHL(N3,N2,N1)
      VOLPH1X=AMAX1(0.0,VOLAH1(N3,N2,N1)-VOLWH1X-VOLIH1(N3,N2,N1))
      IF(FINHX.GT.0.0)THEN
      FINHL(N3,N2,N1)=AMAX1(0.0,AMIN1(FINHX,VOLWH1X,VOLP1X))
      ELSE
      FINHL(N3,N2,N1)=AMIN1(0.0,AMAX1(FINHX,-VOLPH1X,-VOLW1X))
      ENDIF
      FINHM(M,N3,N2,N1)=FINHL(N3,N2,N1)
      FINH(N3,N2,N1)=FINH(N3,N2,N1)+FINHL(N3,N2,N1)
C     IF(NX.EQ.1.AND.NY.EQ.1)THEN
C     WRITE(*,3366)'FINHL',I,J,M,N4,N5,N6,IFLGH,FINHL(N3,N2,N1)
C    3,FINHX,VOLWH1(N3,N2,N1),VOLPH1(N3,N2,N1),VOLP1(N3,N2,N1)
C    4,PSISA1(N3,N2,N1),HCND(2,1,N3,N2,N1),PHOL(N3,N2,N1)
C    5,HRAD(N3,N2,N1)
3366  FORMAT(A8,7I4,20E12.4)
C     ENDIF
      ELSE
      FINHL(N3,N2,N1)=0.0
      FINHM(M,N3,N2,N1)=0.0
      ENDIF
C
C     FREEZE-THAW IN SOIL LAYER MICROPORE FROM NET CHANGE IN SOIL
C     LAYER HEAT STORAGE
C
C     TFREEZ=micropore freezing temperature
C     PSISA1,PSISO=micropore matric,osmotic potential
C     VOLW1*,VOLI1=micropore water,ice volume
C     VOLWH1*,VOLIH1=macropore water,ice volume
C     VHCP1X,VHCP1AX,VHCP1BX=total soil,micropore,macropore heat capacity
C     VHCM=soil solid volumetric heat capacity
C     TK1*=soil temperature
C     THFLWL=total soil conductive, convective heat flux
C     HWFLU1=subsurface convective heat flux
C     TFLX1,TFLX=latent heat from micro freeze-thaw unltd,ltd by water,ice
C     WFLXL=soil water flux from micropore freeze-thaw
C     TFLXH1,TFLXL=latent heat from macro freeze-thaw unltd,ltd by water,ice
C     WFLXLH=soil water flux from macropore freeze-thaw
C
      PSISMX=PSISA1(N3,N2,N1)+PSISO(N3,N2,N1)
      TFREEZ=-9.0959E+04/(PSISMX-333.0)
      VOLW1X=VOLW1(N3,N2,N1)+TFLWL(N3,N2,N1)+FINHL(N3,N2,N1)
     2+FLU1(N3,N2,N1)
      VOLWH1X=VOLWH1(N3,N2,N1)+TFLWHL(N3,N2,N1)-FINHL(N3,N2,N1)
      ENGY1=VHCP1(N3,N2,N1)*TK1(N3,N2,N1)
      VHCP1X=VHCM(N3,N2,N1)+4.19*(VOLW1X+VOLWH1X)
     2+1.9274*(VOLI1(N3,N2,N1)+VOLIH1(N3,N2,N1))
      IF(VHCP1X.GT.ZEROS(NY,NX))THEN
      TK1X=(ENGY1+THFLWL(N3,N2,N1)+HWFLU1(N3,N2,N1))/VHCP1X
      IF((TK1X.LT.TFREEZ
     2.AND.VOLW1(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1))
     4.OR.(TK1X.GT.TFREEZ
     5.AND.VOLI1(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1)))THEN
      VHCP1AX=VHCM(N3,N2,N1)+4.19*VOLW1X+1.9274*VOLI1(N3,N2,N1)
      TFLX1=VHCP1AX*(TFREEZ-TK1X)/((1.0+6.2913E-03*TFREEZ)
     2*(1.0-0.10*PSISMX))*XNPX
      IF(TFLX1.LT.0.0)THEN
      TFLX=AMAX1(-333.0*DENSI*VOLI1(N3,N2,N1)*XNPX,TFLX1)
      ELSE
      TFLX=AMIN1(333.0*VOLW1X*XNPX,TFLX1)
      ENDIF
      WFLXL(N3,N2,N1)=-TFLX/333.0
      ELSE
      TFLX=0.0
      WFLXL(N3,N2,N1)=0.0
      ENDIF
      ELSE
      TFLX=0.0
      WFLXL(N3,N2,N1)=0.0
      ENDIF
C
C     FREEZE-THAW IN SOIL LAYER MACROPORE FROM NET CHANGE IN SOIL
C     LAYER HEAT STORAGE
C
      IF((TK1X.LT.273.15.AND.VOLWH1(N3,N2,N1)
     2.GT.ZERO*VOLT(N3,N2,N1)).OR.(TK1X.GT.273.15
     3.AND.VOLIH1(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1)))THEN
      VHCP1BX=4.19*VOLWH1X+1.9274*VOLIH1(L,NY,NX)
      TFLXH1=VHCP1BX*(TFREEZ-TK1X)/((1.0+6.2913E-03*TFREEZ)
     2*(1.0-0.10*PSISMX))*XNPX
      IF(TFLXH1.LT.0.0)THEN
      TFLXH=AMAX1(-333.0*DENSI*VOLIH1(N3,N2,N1)*XNPX,TFLXH1)
      ELSE
      TFLXH=AMIN1(333.0*VOLWH1X*XNPX,TFLXH1)
      ENDIF
      WFLXLH(N3,N2,N1)=-TFLXH/333.0
      ELSE
      TFLXH=0.0
      WFLXLH(N3,N2,N1)=0.0
      ENDIF
      TFLXL(N3,N2,N1)=TFLX+TFLXH
C
C     TOTAL AND HOURLY ACCUMULATED FREEZE-THAW FLUXES
C
C     THAW,THAWH=hourly accumulated freeze-thaw flux in micropores,macropores
C     HTHAW=hourly accumulated freeze-thaw latent heat flux
C     TWFLXL,TWFLXH=total accumulated freeze-thaw in micropores,macropores
C     TTFLXL=total latent heat flux
C
      THAW(N3,N2,N1)=THAW(N3,N2,N1)+WFLXL(N3,N2,N1)
      THAWH(N3,N2,N1)=THAWH(N3,N2,N1)+WFLXLH(N3,N2,N1)
      HTHAW(N3,N2,N1)=HTHAW(N3,N2,N1)+TFLXL(N3,N2,N1)
      TWFLXL(N3,N2,N1)=TWFLXL(N3,N2,N1)+WFLXL(N3,N2,N1)
      TWFLXH(N3,N2,N1)=TWFLXH(N3,N2,N1)+WFLXLH(N3,N2,N1)
      TTFLXL(N3,N2,N1)=TTFLXL(N3,N2,N1)+TFLXL(N3,N2,N1)
C     IF(N3.EQ.1)THEN
C     WRITE(*,4359)'TFLX',I,J,M,N1,N2,N3,TWFLXL(N3,N2,N1)
C    2,WFLXL(N3,N2,N1),THAW(N3,N2,N1),VHCP1X,VHCP1AX
C    2,TFLX1,TFLX,TFREEZ,TK1(N3,N2,N1),TK1X,VHCP1(N3,N2,N1)
C    3,TFLXL(N3,N2,N1),TTFLXL(N3,N2,N1)
C    2,TFLWL(N3,N2,N1),FINHL(N3,N2,N1),FLU1(N3,N2,N1),THFLWL(N3,N2,N1)
C    4,HWFLU1(N3,N2,N1),VOLW1X,VOLW1(N3,N2,N1),VOLI1(N3,N2,N1)
C    5,PSISA1(N3,N2,N1),PSISM(N3,N2,N1),PSISO(N3,N2,N1),PSISMX
C    5,VOLWH1(N3,N2,N1),VOLIH1(N3,N2,N1)
C    5,FGRD(N3,N2,N1),FMAC(N3,N2,N1)
4359  FORMAT(A8,6I4,40E14.6)
C     ENDIF
C
C     DISSIPATE WETTING FRONT
C
C     VOLW1=soil micropore water content
C     VOLWX1=soil micropore water content behind wetting front
C     FLWVL=water flux from wetted to drier soil
C
9585  CONTINUE
9590  CONTINUE
9595  CONTINUE
C
C     UPDATE STATE VARIABLES FROM FLUXES CALCULATED ABOVE
C
      IF(M.NE.NPH)THEN
      DO 9795 NX=NHW,NHE
      DO 9790 NY=NVN,NVS
C
C     SNOWPACK WATER, ICE, SNOW AND TEMPERATURE
C
C     VOLS0,VOLW0,VOLI0=snow,water,ice volumes in snowpack
C     TFLWS,TFLWW=net snow,water flux
C     WFLXS,WFLXI=snow-water,ice-water freeze-thaw flux
C     DENSI=ice density
C     VHCPWM=snowpack volumetric heat capacity
C     TK0=snowpack temperature
C     THFLWW=total snowpack conductive+convective heat flux
C     TFLX0=snowpack latent heat flux from freeze-thaw
C
      DO 9780 L=1,JS
      VOLS0(L,NY,NX)=VOLS0(L,NY,NX)+TFLWS(L,NY,NX)-WFLXS(L,NY,NX)
      VOLW0(L,NY,NX)=VOLW0(L,NY,NX)+TFLWW(L,NY,NX)+WFLXS(L,NY,NX)
     2+WFLXI(L,NY,NX)
      VOLI0(L,NY,NX)=VOLI0(L,NY,NX)-WFLXI(L,NY,NX)/DENSI
      ENGY0=VHCPWM(M,L,NY,NX)*TK0(L,NY,NX)
      VHCPWM(M+1,L,NY,NX)=2.095*VOLS0(L,NY,NX)+4.19*VOLW0(L,NY,NX)
     2+1.9274*VOLI0(L,NY,NX)
      IF(VHCPWM(M+1,L,NY,NX).GT.VHCPWX(NY,NX))THEN
      TK0(L,NY,NX)=(ENGY0+THFLWW(L,NY,NX)+TFLX0(L,NY,NX))
     2/VHCPWM(M+1,L,NY,NX)
      ELSEIF(L.EQ.1)THEN
      TK0(L,NY,NX)=TKA(NY,NX)
      ELSE
      TK0(L,NY,NX)=TK0(L-1,NY,NX)
      ENDIF
C     IF(L.EQ.1)THEN
C     WRITE(*,7753)'TKM',I,J,M,NX,NY,L,TK0(L,NY,NX)
C    2,THFLWW(L,NY,NX),TFLX0(L,NY,NX),HFLW0W(L,NY,NX)
C    3,HFLWRLW,HFLWLW
C    3,VOLS0(L,NY,NX),VOLW0(L,NY,NX),VOLI0(L,NY,NX),VOLS1(L,NY,NX)
C    2,WFLXS(L,NY,NX),WFLXI(L,NY,NX)
C    3,TFLWS(L,NY,NX),TFLWW(L,NY,NX),TFLWI(L,NY,NX)
C    3,XFLWS(L,NY,NX),XFLWW(L,NY,NX),XFLWI(L,NY,NX)
C    4,XHFLWW(L,NY,NX),HFLSW(L,NY,NX)
C    5,VHCPWM(M+1,L,NY,NX)
7753  FORMAT(A8,6I4,80E14.6)
C     ENDIF
9780  CONTINUE
C
C     SNOW RUNOFF
C
C     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from snowpack runoff
C
      VOLS0(1,NY,NX)=VOLS0(1,NY,NX)+TQS1(NY,NX)
      VOLW0(1,NY,NX)=VOLW0(1,NY,NX)+TQW1(NY,NX)
      VOLI0(1,NY,NX)=VOLI0(1,NY,NX)+TQI1(NY,NX)
      ENGY0=VHCPWM(M+1,1,NY,NX)*TK0(1,NY,NX)
      VHCPWM(M+1,1,NY,NX)=2.095*VOLS0(1,NY,NX)+4.19*VOLW0(1,NY,NX)
     2+1.9274*VOLI0(1,NY,NX)
      IF(VHCPWM(M+1,1,NY,NX).GT.VHCPWX(NY,NX))THEN
      TK0(1,NY,NX)=(ENGY0+THQS1(NY,NX))
     2/VHCPWM(M+1,1,NY,NX)
      ELSE
      TK0(1,NY,NX)=TKA(NY,NX)
      ENDIF
C
C     IF SNOWPACK DISAPPEARS ALL MATERIAL,HEAT TRANSFERRED TO SOIL SURFACE
C
C     VHCP1,VHCP1A,VHCP1P=total soil,soil+micropore,macropore heat capacity
C     TK1=soil surface temperature
C
      IF(VHCPW(1,NY,NX).LE.VHCPWX(NY,NX)
     2.AND.TKA(NY,NX).GT.273.15)THEN
      FLWS=VOLS0(1,NY,NX)
      FLWW=VOLW0(1,NY,NX)
      FLWI=VOLI0(1,NY,NX)
      HFLWS=(4.19*FLWW+2.095*FLWS+1.9274*FLWI)*TK0(1,NY,NX)
      VOLS0(1,NY,NX)=VOLS0(1,NY,NX)-FLWS
      VOLW0(1,NY,NX)=VOLW0(1,NY,NX)-FLWW
      VOLI0(1,NY,NX)=VOLI0(1,NY,NX)-FLWI
      VOLW1(NUM(NY,NX),NY,NX)=VOLW1(NUM(NY,NX),NY,NX)+FLWW
      VOLI1(NUM(NY,NX),NY,NX)=VOLI1(NUM(NY,NX),NY,NX)+FLWI
     2+FLWS/DENSI
      ENGY1=VHCP1(NUM(NY,NX),NY,NX)*TK1(NUM(NY,NX),NY,NX)
      VHCP1x0=VHCP1(NUM(NY,NX),NY,NX)
      VHCP1(NUM(NY,NX),NY,NX)=VHCM(NUM(NY,NX),NY,NX)
     2+4.19*(VOLW1(NUM(NY,NX),NY,NX)+VOLWH1(NUM(NY,NX),NY,NX))
     3+1.9274*(VOLI1(NUM(NY,NX),NY,NX)+VOLIH1(NUM(NY,NX),NY,NX))
      VHCP1A(NUM(NY,NX),NY,NX)=VHCM(NUM(NY,NX),NY,NX)
     2+4.19*VOLW1(NUM(NY,NX),NY,NX)+1.9274*VOLI1(NUM(NY,NX),NY,NX)
      VHCP1B(NUM(NY,NX),NY,NX)=4.19*VOLWH1(NUM(NY,NX),NY,NX)
     3+1.9274*VOLIH1(NUM(NY,NX),NY,NX)
      IF(VHCP1(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      TK10=TK1(NUM(NY,NX),NY,NX)
      TK1(NUM(NY,NX),NY,NX)=(ENGY1+HFLWS)
     2/VHCP1(NUM(NY,NX),NY,NX)
!     if(abs(tK10/TK1(NUM(NY,NX),NY,NX)-1.)>0.025)then
!     TK1(NUM(NY,NX),NY,NX)=TKA(NY,NX)
!     endif
      ELSE
      TK1(NUM(NY,NX),NY,NX)=TKA(NY,NX)
      ENDIF
      ENDIF
C
C     SURFACE RESIDUE WATER AND TEMPERATURE
C
C     XVOLT,XVOLW=free water+ice,water in litter layer
C     VOLWM,VOLPM=surface water,air content for use in trnsfr.f
C     VOLWRX=maximum water retention by litter
C     VHCP1=volumetric heat capacity of litter
C     VOLA1,VOLW1,VOLI1,VOLP1=pore,water,ice,air volumes of litter
C     VOLWRX=maximum water retention by litter
C     TFLXR,WFLXR=litter water,latent heat flux from freeze-thaw
C     VOLR=dry litter volume
C     THETWX,THETIX,THETPX=water,ice,air concentrations
C     VHCP1=volumetric heat capacity of litter
C     TK1=litter temperature
C     HFLWRL,TFLXR,THQR1=litter total cond+conv,latent,runoff heat flux
C
      VOLW1xx=VOLW1(0,NY,NX)
      VOLW1x=VOLW1(0,NY,NX)+FLWRL(NY,NX)
     2+WFLXR(NY,NX)+TQR1(NY,NX)
      VOLW1(0,NY,NX)=AMAX1(0.0,VOLW1(0,NY,NX)+FLWRL(NY,NX)
     2+WFLXR(NY,NX)+TQR1(NY,NX))
      VOLI1(0,NY,NX)=AMAX1(0.0,VOLI1(0,NY,NX)-WFLXR(NY,NX)/DENSI)
      VOLP1(0,NY,NX)=AMAX1(0.0,VOLA1(0,NY,NX)-VOLW1(0,NY,NX)
     2-VOLI1(0,NY,NX))
      VOLWM(M+1,0,NY,NX)=VOLW1(0,NY,NX)
      VOLPM(M+1,0,NY,NX)=VOLP1(0,NY,NX)
      TVOLWI=VOLW1(0,NY,NX)+VOLI1(0,NY,NX)
      XVOLT(NY,NX)=AMAX1(0.0,TVOLWI-VOLWRX(NY,NX))
      IF(TVOLWI.GT.ZEROS(NY,NX))THEN
      VOLWRZ=VOLW1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      VOLIRZ=VOLI1(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      XVOLW(NY,NX)=AMAX1(0.0,VOLW1(0,NY,NX)-VOLWRZ)
      XVOLI(NY,NX)=AMAX1(0.0,VOLI1(0,NY,NX)-VOLIRZ)
      ELSE
      XVOLW(NY,NX)=0.0
      XVOLI(NY,NX)=0.0
      ENDIF
      XVOLTM(M+1,NY,NX)=XVOLT(NY,NX)
      XVOLWM(M+1,NY,NX)=XVOLW(NY,NX)
      XVOLIM(M+1,NY,NX)=XVOLI(NY,NX)
      IF(VOLR(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWX(0,NY,NX)=AMAX1(0.0,VOLW1(0,NY,NX)/VOLR(NY,NX))
      THETIX(0,NY,NX)=AMAX1(0.0,VOLI1(0,NY,NX)/VOLR(NY,NX))
      THETPX(0,NY,NX)=AMAX1(0.0,VOLP1(0,NY,NX)/VOLR(NY,NX))
     2*AMAX1(0.0,(1.0-XVOLT(NY,NX)/VOLWD(NY,NX)))
      ELSE
      THETWX(0,NY,NX)=0.0
      THETIX(0,NY,NX)=0.0
      THETPX(0,NY,NX)=1.0
      ENDIF
      THETPM(M+1,0,NY,NX)=THETPX(0,NY,NX)
      VHCPXX=VHCP1(0,NY,NX)
      TK0XX=TK1(0,NY,NX)
      ENGYR=VHCP1(0,NY,NX)*TK1(0,NY,NX)
      VHCP1(0,NY,NX)=2.496E-06*ORGC(0,NY,NX)+4.19*VOLW1(0,NY,NX)
     2+1.9274*VOLI1(0,NY,NX)
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      TK1(0,NY,NX)=(ENGYR+HFLWRL(NY,NX)+TFLXR(NY,NX)
     2+THQR1(NY,NX))/VHCP1(0,NY,NX)

!     IF(ABS(VHCP1(0,NY,NX)/VHCPXX-1.)>0.025.or.
!    2abs(TK1(0,NY,NX)/tk0xx-1)>0.025)THEN
!     TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
!     ENDIF

      ELSE
      TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
      ENDIF
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,7754)'VOLW0',I,J,M,NX,NY,NUM(NY,NX)
C    2,VOLW1(0,NY,NX),VOLI1(0,NY,NX),VOLP1(0,NY,NX)
C    3,THETWX(0,NY,NX),THETIX(0,NY,NX),THETPX(0,NY,NX)
C    2,FLWRL(NY,NX),WFLXR(NY,NX),TQR1(NY,NX),EVAPR(NY,NX)
C    4,XVOLW(NY,NX),XVOLI(NY,NX),XVOLT(NY,NX)
C    5,TVOLWI,VOLWRX(NY,NX)
C    5,FLWRLG,FLWRLW,FLYM,FLQR
C    3,TK1(0,NY,NX),HFLWRL(NY,NX),TFLXR(NY,NX)
C    3,THQR1(NY,NX),ENGYR,VHCP1(0,NY,NX),VHCPRX(NY,NX)
C    2,HFLWRLG,HFLWRLW,HWFLYM,HFLXR
C    3,VOLPM(M,0,NY,NX),VOLPM(M,NUM(NY,NX),NY,NX)
C    4,VOLI1(NUM(NY,NX),NY,NX),TK1(NUM(NY,NX),NY,NX)
7754  FORMAT(A8,6I4,60E12.4)
C     ENDIF
C
C     SOIL LAYER WATER, ICE AND TEMPERATURE
C
C     VOLW1,VOLI1=micropore water,ice volume
C     VOLWX1=micropore water volume behind wetting front
C     VOLWH1,VOLIH1=macropore water,ice volume
C     TFLWL=net water flux
C     FINHL=micropore-macropore flux
C     TWFLXL,TWFLXH=total accumulated freeze-thaw in micropores,macropores
C     FLU1=subsurface water input
C     DENSI=ice density
C     VOLA1,VOLAH1=micropore,macropore volume
C     VOLP1,VOLPH1=micropore,macropore air volume
C     VOLWM,VOLWHM,VOLPM,FLPM=micropore,macropore water volume, air volume
C     and change in air volume for use in trnsfr.f
C     THETWX,THETIX,THETPX,THETPY=bulk water,ice,air concn,air-filled porosity
C     THETPM=air concentration for use in trnsfr.f
C     FMAC,FGRD=macropore,micropore fraction
C     CNDH1=maropore hydraulic conductivity
C     VHCP1,VHCM=volumetric heat capacities of total volume, solid
C     VHCP1A,VHCP1B=volumetric heat capacities of soil+micropore,macropore
C     TK1=soil temperature
C
      DO 9785 L=NUM(NY,NX),NL(NY,NX)
      IF(VOLT(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      VOLw10=VOLW1(L,NY,NX)
      VOLW1(L,NY,NX)=VOLW1(L,NY,NX)+TFLWL(L,NY,NX)
     2+FINHL(L,NY,NX)+TWFLXL(L,NY,NX)+FLU1(L,NY,NX)
C      if(VOLW1(L,NY,NX)/=VOLW1(L,NY,NX))then
C      print*,'volw',L,VOLw10
C      print*,TFLWL(L,NY,NX),
C    2FINHL(L,NY,NX),TWFLXL(L,NY,NX),FLU1(L,NY,NX)
C      pause
C      endif
      VOLWX1(L,NY,NX)=VOLWX1(L,NY,NX)+TFLWLX(L,NY,NX)
     2+FINHL(L,NY,NX)+TWFLXL(L,NY,NX)+FLU1(L,NY,NX)
      VOLWX1(L,NY,NX)=AMIN1(VOLW1(L,NY,NX),VOLWX1(L,NY,NX))
      VOLI1(L,NY,NX)=VOLI1(L,NY,NX)-TWFLXL(L,NY,NX)/DENSI
      VOLWH1(L,NY,NX)=VOLWH1(L,NY,NX)+TFLWHL(L,NY,NX)
     2-FINHL(L,NY,NX)+TWFLXH(L,NY,NX)
      VOLIH1(L,NY,NX)=VOLIH1(L,NY,NX)-TWFLXH(L,NY,NX)/DENSI
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLP1Z(L,NY,NX)=VOLA1(L,NY,NX)-VOLW1(L,NY,NX)-VOLI1(L,NY,NX)
      VOLP1(L,NY,NX)=AMAX1(0.0,VOLP1Z(L,NY,NX))
      VOLPH1Z(L,NY,NX)=VOLAH1(L,NY,NX)-VOLWH1(L,NY,NX)-VOLIH1(L,NY,NX)
      VOLPH1(L,NY,NX)=AMAX1(0.0,VOLPH1Z(L,NY,NX))
      VOLAH1(L,NY,NX)=AMAX1(0.0,VOLAH(L,NY,NX)-FVOLAH*CCLAY(L,NY,NX)
     2*(VOLW1(L,NY,NX)/VOLY(L,NY,NX)-WP(L,NY,NX))*VOLT(L,NY,NX))
      ELSE
      VOLP1Z(L,NY,NX)=0.0
      VOLP1(L,NY,NX)=0.0
      VOLPH1Z(L,NY,NX)=0.0
      VOLPH1(L,NY,NX)=0.0
      VOLA1(L,NY,NX)=VOLW1(L,NY,NX)+VOLI1(L,NY,NX)
      VOLAH1(L,NY,NX)=0.0
      ENDIF
      VOLWM(M+1,L,NY,NX)=VOLW1(L,NY,NX)
      VOLWHM(M+1,L,NY,NX)=VOLWH1(L,NY,NX)
      VOLPM(M+1,L,NY,NX)=VOLP1(L,NY,NX)+VOLPH1(L,NY,NX)
     2+THETPI*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
      FLPM(M,L,NY,NX)=VOLPM(M,L,NY,NX)-VOLPM(M+1,L,NY,NX)
      VOLTX=VOLY(L,NY,NX)+VOLAH1(L,NY,NX)
      THETWX(L,NY,NX)=AMAX1(0.0,(VOLW1(L,NY,NX)+VOLWH1(L,NY,NX))
     2/VOLTX)
      THETIX(L,NY,NX)=AMAX1(0.0,(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
     2/VOLTX)
      THETPX(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX))
     2/VOLTX)
      THETPM(M+1,L,NY,NX)=THETPX(L,NY,NX)
      IF(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETPY(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX))
     2/(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX)))
      ELSE
      THETPY(L,NY,NX)=0.0
      ENDIF
      IF(VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      FMAC(L,NY,NX)=FHOL(L,NY,NX)*VOLAH1(L,NY,NX)/VOLAH(L,NY,NX)
      CNDH1(L,NY,NX)=CNDH(L,NY,NX)
     2*(VOLAH1(L,NY,NX)/VOLAH(L,NY,NX))**2
      ELSE
      FMAC(L,NY,NX)=0.0
      CNDH1(L,NY,NX)=0.0
      ENDIF
      FGRD(L,NY,NX)=1.0-FMAC(L,NY,NX)
      TKXX=TK1(L,NY,NX)
      VHXX=VHCP1(L,NY,NX)
      ENGY1=VHCP1(L,NY,NX)*TK1(L,NY,NX)
      VHCP1(L,NY,NX)=VHCM(L,NY,NX)+4.19*(VOLW1(L,NY,NX)
     2+VOLWH1(L,NY,NX))+1.9274*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
      VHCP1A(L,NY,NX)=VHCM(L,NY,NX)+4.19*VOLW1(L,NY,NX)
     2+1.9274*VOLI1(L,NY,NX)
      VHCP1B(L,NY,NX)=4.19*VOLWH1(L,NY,NX)+1.9274*VOLIH1(L,NY,NX)
C
! WJR modified
C     BEGIN ARTIFICIAL SOIL WARMING
C
C     THFLWL=THFLWL incremented for soil warming
C     TKSZ=temperature used to calculate additional heat flux for warming
C
      IF (IYRC .ge. 2014 .and. L .le. 8) then
         THFLWL(L,NY,NX)=THFLWL(L,NY,NX)
     2       +(TKSZ(I,J,L)-TK1(L,NY,NX))*VHCP1(L,NY,NX)*XNPH
!        if (L .eq. 3) then
!           write (777, *) IYRC, I, J, L, TKSZ(I,J,L),TK1(L,NY,NX),
!    1       TKSZ(I,J,L)-TK1(L,NY,NX),
!    2       (TKSZ(I,J,L)-TK1(L,NY,NX))*VHCP1(L,NY,NX)*XNPH,
!    3       THFLWL(L,NY,NX)
!        endif
      ENDIF
C
C     END ARTIFICIAL SOIL WARMING
C
! END WJR modified

      IF(VHCP1(L,NY,NX).GT.ZEROS(NY,NX))THEN
!     tk1l=TK1(L,NY,NX)
      TK1(L,NY,NX)=(ENGY1+THFLWL(L,NY,NX)+TTFLXL(L,NY,NX)
     2+HWFLU1(L,NY,NX))/VHCP1(L,NY,NX)
!     if(abs(TK1(L,NY,NX)/tk1l-1.)>0.025)then
!     TK1(L,NY,NX)=tk1l
!     endif
      ELSEIF(L.EQ.1)THEN
      TK1(L,NY,NX)=TKA(NY,NX)
      ELSE
      TK1(L,NY,NX)=TK1(L-1,NY,NX)
      ENDIF
      ELSE
      VOLWM(M+1,L,NY,NX)=0.0
      VOLWHM(M+1,L,NY,NX)=0.0
      VOLPM(M+1,L,NY,NX)=0.0
      FLPM(M,L,NY,NX)=VOLPM(M,L,NY,NX)
      THETPM(M+1,L,NY,NX)=0.0
      ENDIF
C     IF(L.LE.NUM(NY,NX)+1)THEN
C     WRITE(*,3377)'VOLW1',I,J,M,NX,NY,L,N6X(NY,NX)
C    2,VOLW1(L,NY,NX),VOLWX1(L,NY,NX)
C    3,VOLI1(L,NY,NX),VOLA1(L,NY,NX),VOLP1(L,NY,NX)
C    2,TFLWL(L,NY,NX),TFLWLX(L,NY,NX)
C    3,TWFLXL(L,NY,NX),FINHL(L,NY,NX),FLU1(L,NY,NX)
C    3,VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),VOLAH1(L,NY,NX),VOLPH1(L,NY,NX)
C    4,VOLA1(L,NY,NX)-VOLW1(L,NY,NX)-VOLI1(L,NY,NX)
C    5,VOLPM(M,L,NY,NX),VOLPM(M+1,L,NY,NX)
C    5,PSISM1(L,NY,NX),THETPX(L,NY,NX)
C    6,FLWL(3,L,NY,NX),FLWL(3,L+1,NY,NX)
C    7,FLWL(2,L,NY,NX),FLWL(2,L,NY+1,NX)
C    8,FLWL(1,L,NY,NX),FLWL(1,L,NY,NX+1)
C    6,FLWLX(3,L,NY,NX),FLWLX(3,L+1,NY,NX)
C    7,FLWLX(2,L,NY,NX),FLWLX(2,L,NY+1,NX)
C    8,FLWLX(1,L,NY,NX),FLWLX(1,L,NY,NX+1)
C    6,FLW(3,L,NY,NX),FLW(3,L+1,NY,NX)
C    7,FLW(2,L,NY,NX),FLW(2,L,NY+1,NX)
C    8,FLW(1,L,NY,NX),FLW(1,L,NY,NX+1)
C    9,WFLXL(L,NY,NX),THAW(L,NY,NX)
C    9,FLPM(M,L,NY,NX),FLSW(L,NY,NX)
C     WRITE(*,3377)'VOLWH1',I,J,M,NX,NY,L,N6X(NY,NX),VOLWH1(L,NY,NX)
C    2,TFLWHL(L,NY,NX),FINHL(L,NY,NX),VOLIH1(L,NY,NX)
C    4,TWFLXH(L,NY,NX),TQR1(NY,NX),VOLPH1(L,NY,NX)
C    6,FLWHL(3,L,NY,NX),FLWHL(3,L+1,NY,NX)
C    7,FLWHL(2,L,NY,NX),FLWHL(2,L,NY+1,NX)
C    8,FLWHL(1,L,NY,NX),FLWHL(1,L,NY,NX+1)
C    3,VHCM(L,NY,NX),VOLW1(L,NY,NX),VOLWH1(L,NY,NX),VOLI1(L,NY,NX)
C    4,THETW(L,NY,NX),THETI(L,NY,NX),FINHL(L,NY,NX),THQR1(NY,NX)
C    5,HFLWL(3,L,NY,NX),HFLWL(3,N6X(NY,NX),NY,NX)
C    6,HFLWL(1,L,NY,NX),HFLWL(1,L,NY,NX+1)
C     WRITE(*,3377)'TK1',I,J,M,NX,NY,L,N6X(NY,NX),TK1(L,NY,NX)
C    2,THFLWL(L,NY,NX),TTFLXL(L,NY,NX),HWFLU1(L,NY,NX)
C    3,VHCP1(L,NY,NX),HFLWL(3,L,NY,NX),HFLWL(3,L+1,NY,NX)
C    4,HFLWLW,HFLWLG,TKXX,VHXX,ENGY1,TK1(0,NY,NX)
3377  FORMAT(A8,7I4,40E12.4)
C     ENDIF
9785  CONTINUE
C
C     RESET SURFACE LAYER NUMBER AND TRANSFER ALL WATER TO SOIL SURFACE LAYER
C     IF LAKE SURFACE LAYER IS LOST TO EVAPORATION
C
C     NUM=new surface layer number after complete lake evaporation
C     FLWNU,FLWHNU,HFLWNU=lake surface water flux, heat flux if lake surface disappears
C
      IF(BKDS(NUM(NY,NX),NY,NX).LE.ZERO
     2.AND.VHCP1(NUM(NY,NX),NY,NX).LE.VHCPNX(NY,NX))THEN
      NUX=NUM(NY,NX)
      DO 9970 LL=NUX+1,NL(NY,NX)
      IF(VOLX(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
      NUM(NY,NX)=LL
      FLWNX(NY,NX)=FLW(3,NUM(NY,NX),NY,NX)
      FLWXNX(NY,NX)=FLWX(3,NUM(NY,NX),NY,NX)
      FLWHNX(NY,NX)=FLWH(3,NUM(NY,NX),NY,NX)
      HFLWNX(NY,NX)=HFLW(3,NUM(NY,NX),NY,NX)
C     WRITE(*,5598)'SURFM',I,J,M,NX,NY,LL,NUX,NUM(NY,NX)
C    2,FLWNX(NY,NX)
C    2,VOLW1(NUX,NY,NX),VOLW1(NUM(NY,NX),NY,NX)
C    2,VHCP1(NUX,NY,NX),VHCP1(NUM(NY,NX),NY,NX)
C    2,TK1(NUX,NY,NX),TK1(NUM(NY,NX),NY,NX)
C    3,FLW(3,NUX,NY,NX),FLW(3,NUM(NY,NX),NY,NX)
C    3,HFLW(3,NUX,NY,NX),HFLW(3,NUM(NY,NX),NY,NX)
C    4,VHCPNX(NY,NX),VHCP1(0,NY,NX)
5598  FORMAT(A8,8I4,20E12.4)
      GO TO 9971
      ENDIF
9970  CONTINUE
      ENDIF
9971  CONTINUE
9790  CONTINUE
9795  CONTINUE
      ELSE
      DO 9695 NX=NHW,NHE
      DO 9690 NY=NVN,NVS
      IF(NUM(NY,NX).EQ.NU(NY,NX))THEN
      FLWNU(NY,NX)=FLW(3,N6X(NY,NX),NY,NX)
      FLWXNU(NY,NX)=FLWX(3,N6X(NY,NX),NY,NX)
      FLWHNU(NY,NX)=FLWH(3,N6X(NY,NX),NY,NX)
      HFLWNU(NY,NX)=HFLW(3,N6X(NY,NX),NY,NX)
      ELSE
      FLWNU(NY,NX)=FLWNX(NY,NX)
      FLWXNU(NY,NX)=FLWXNX(NY,NX)
      FLWHNU(NY,NX)=FLWHNX(NY,NX)
      HFLWNU(NY,NX)=HFLWNX(NY,NX)
      ENDIF
9690  CONTINUE
9695  CONTINUE
      ENDIF
3320  CONTINUE
      RETURN
      END
