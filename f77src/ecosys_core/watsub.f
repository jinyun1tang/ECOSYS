      SUBROUTINE watsub(I,J,NFZ,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE CACULATES ENERGY BALANCES OF SNOW, RESIDUE
C     AND SOIL SURFACES, FREEZING, THAWING, AND HEAT AND WATER
C     TRANSFER THROUGH SOIL PROFILES
C
      include "parameters.h"
      include "blkc.h"
      include "blk1g.h"
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
      include "blk18a.h"
      include "blk22a.h"
      include "blk22b.h"
      include "blk22c.h"
      include "blktest.h"
      DIMENSION VOLWX1(JZ,JY,JX),ALBZ(JY,JX) 
     2,XVOLT(JY,JX),XVOLW(JY,JX),XVOLI(JY,JX),FMAC(JZ,JY,JX)
     3,FGRD(JZ,JY,JX),VOLW1(0:JZ,JY,JX),VOLI1(0:JZ,JY,JX) 
     4,VHCP1(0:JZ,JY,JX),VHCP1A(JZ,JY,JX),VHCP1B(JZ,JY,JX)
     4,TK1(0:JZ,JY,JX),TWFLFL(JZ,JY,JX),VOLW2(0:JZ,JY,JX) 
     5,VOLP1(0:JZ,JY,JX),TWFLFH(JZ,JY,JX),PRECM(JY,JX) 
     6,VOLS0(JS,JY,JX),VOLI0(JS,JY,JX),VOLW0(JS,JY,JX)
     7,VOLS1(JS,JY,JX),DLYRS0(JS,JY,JX),VOLP1Z(JZ,JY,JX) 
     8,TK0(JS,JY,JX),AREAU(JZ,JY,JX),AREAUD(JZ,JY,JX),FLQ0S(JY,JX) 
     9,FLQ0I(JY,JX),FLQ0W(JY,JX),FLQ1(JY,JX),FLH1(JY,JX) 
     9,FLY1(JY,JX),HWFLQ0(JY,JX),HWFLQ1(JY,JX),HWFLY1(JY,JX)
     1,RZR(JY,JX),BAREW(JY,JX),CVRDW(JY,JX),DPTHS0(JY,JX)
     4,QR1(2,2,JV,JH),HQR1(2,2,JV,JH),VOLPH1Z(JZ,JY,JX)
     5,QS1(2,2,JV,JH),QW1(2,2,JV,JH),QI1(2,2,JV,JH),HQS1(2,2,JV,JH)
     6,TQR1(JY,JX),THQR1(JY,JX),EVAPG(JY,JX),THFLFL(JZ,JY,JX)
     7,VOLTX(JZ,JY,JX),EVAP0(JY,JX),EVAPR(JY,JX),FLWRL(JY,JX) 
     8,FLVRL(JY,JX),HFLWRL(JY,JX),FINHL(JZ,JY,JX),HFLWL(3,JD,JV,JH)
     9,FLWL(3,JD,JV,JH),FLVL(3,JD,JV,JH),TWFLVL(JZ,JY,JX)
     1,THFLVL(JZ,JY,JX),RAS(JY,JX),FLWLY(3,JD,JV,JH)
     2,VOLV2(0:JZ,JY,JX),VOLI2(0:JZ,JY,JX),VOLP2(0:JZ,JY,JX)
     3,VOLWH2(JZ,JY,JX),VOLIH2(JZ,JY,JX),VOLPH2(JZ,JY,JX) 
      DIMENSION FLWHL(3,JD,JV,JH),FLWHLY(3,JD,JV,JH)
     2,TFLWL(JZ,JY,JX),TFLWHL(JZ,JY,JX),THFLWL(JZ,JY,JX)
     3,WFLFL(JZ,JY,JX),HFLFL(JZ,JY,JX),AVCNHL(3,JD,JV,JH)
     5,THRYW(JY,JX),EMMGW(JY,JX),EMMGS(JY,JX),EMMGR(JY,JX)
     5,EMMCW(JY,JX),EMMCS(JY,JX),EMMCR(JY,JX)
     6,THRYG(JY,JX),THRYR(JY,JX),RADXW(JY,JX),RADXG(JY,JX)
     7,RADXR(JY,JX),FLWLX(3,JD,JV,JH),TFLWLX(JZ,JY,JX) 
     8,FLU1(JZ,JY,JX),HWFLU1(JZ,JY,JX),PSISM1(0:JZ,JY,JX)
     9,ALTG(JY,JX),WFLFLH(JZ,JY,JX),DLYRR(JY,JX),WFLFR(JY,JX)
     1,HFLFR(JY,JX),CNDH1(JZ,JY,JX),VOLA1(0:JZ,JY,JX) 
     2,THETWX(0:JZ,JY,JX),THETIX(0:JZ,JY,JX),THETPX(0:JZ,JY,JX)
     4,VOLAH1(JZ,JY,JX),VOLWH1(JZ,JY,JX),VOLPH1(JZ,JY,JX)
     5,VOLIH1(JZ,JY,JX),THETPY(0:JZ,JY,JX),FLWNX(JY,JX)
     6,FLWXNX(JY,JX),FLWHNX(JY,JX),HFLWNX(JY,JX),N6X(JY,JX)
     7,PSISA1(JZ,JY,JX),VOLV0(JS,JY,JX),VOLV1(0:JZ,JY,JX)
     8,WFLVL(JZ,JY,JX),HFLVL(JZ,JY,JX),TFLVL(JZ,JY,JX)
     9,FLVNX(JY,JX),HFLXF(0:JZ,JY,JX),HFLXH(JY,JX)
      DIMENSION TQS1(JY,JX),TQW1(JY,JX),TQI1(JY,JX),THQS1(JY,JX)
     2,TFLWS(JS,JY,JX),TFLWW(JS,JY,JX),TFLWI(JS,JY,JX)
     3,THFLWW(JS,JY,JX),HFLF0(JS,JY,JX),WFLFS(JS,JY,JX)
     4,WFLFI(JS,JY,JX),FLW0W(JS,JY,JX),FLW0S(JS,JY,JX)
     5,FLW0I(JS,JY,JX),HFLW0W(JS,JY,JX),HFLV0(JS,JY,JX) 
     6,WFLVW(JS,JY,JX),WFLVS(JS,JY,JX),VOLS02(JS,JY,JX)
     7,VOLV02(JS,JY,JX),VOLW02(JS,JY,JX),VOLI02(JS,JY,JX)
     8,VOLP02(JS,JY,JX),VHCPWM2(JS,JY,JX),TK02(JS)
     9,TFLWV(JS,JY,JX),WFLVR(JY,JX),HFLVR(JY,JX),FLW0V(JS,JY,JX)
     1,QSM(JY,JX),QWM(JY,JX),QIM(JY,JX),VHCPRXS(JY,JX)
     2,VHCPRXG(JY,JX),DTKQC(JZ),DVPQC(JZ),DTKQD(JZ),DVPQD(JZ)
C
C     EMMS,EMMW,EMMR=emissivities of surface soil, snow and litter
C     DPTHSX=minimum snowpack depth for full cover (m)
C     Z1S,Z2SW,Z2SD,Z3SX=parameters for air-water gas transfers 
C        in soil
C     Z1R,Z2RW,Z2RD,Z3RX=parameters for air-water gas transfers 
C        in surface litter
C
      PARAMETER (EMMS=0.97,EMMW=0.97,EMMR=0.97,DPTHSX=0.05)
      PARAMETER (Z1S=0.4,Z2SW=12.0,Z2SD=12.0,Z3SX=0.50
     2,Z1R=0.4,Z2RW=12.0,Z2RD=12.0,Z3R=0.50)
C
C     Parameters for calculating convective effects on heat transfer
C        in porous media (air and water):
C     VISCW,VISCA=water,air viscosity (Mg m-1 s)
C     RAY*,*NUS*=Rayleigh,Nusslelt numbers
C     TRBW,TRBA=minimum water,air concentrations for convective
C        effects (m3 m-3)
C
      PARAMETER (VISCW=1.0E-06,VISCA=2.0E-08,DIFFW=1.45E-07
     2,DIFFA=2.01E-05,EXPNW=2.07E-04,EXPNA=3.66E-03,GRAV=9.8
     3,RYLXW=GRAV*EXPNW/(VISCW*DIFFW),RYLXA=GRAV*EXPNA/(VISCA*DIFFA)
     4,PRNTW=VISCW/DIFFW,PRNTA=VISCA/DIFFA
     5,DNUSW=(1.0+(0.492/PRNTW)**0.5625)**0.4444
     6,DNUSA=(1.0+(0.492/PRNTA)**0.5625)**0.4444
     7,TRBW=0.625,TRBA=0.625)
C
C     FVOLAH=parameter for clay effect on macropore volume
C     DTHETW=difference between saturation and effective saturation
C        for calculating water potentials (m3 m-3)
C     FENGYP=rate constant for restoring surface Ksat after compaction
C        (h-1)
C     ORFLN=reflection coefficient for osmotic potential-driven water
C        flux 
C     RAH,RAHW,RAE=isothermal resistance to heat,evaporation at ground
C        and snowpack surfaces (h m-1)
C     RACGM,RACGZ=minimum,maximum canopy aerodynamic resistance
C        for canopy,ground-atmosphere heat,vapor exchange (h m-1)  
C
      PARAMETER (FVOLAH=0.0,DTHETW=1.0E-06,FENGYP=1.0E-03,ORFLN=0.03)
      PARAMETER (RAH=0.278E-02,RAHW=0.278E-01,RAE=0.139E-01
     2,RACGM=0.278E-02,RACGZ=0.278E-01)
      REAL*4 THETWR,THETW1,THETA1,THETAL,THETWL
     2,TKR2,TKS2,TKY 
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      NUM(NY,NX)=NU(NY,NX)
C
C     ADJUST SURFACE ELEVATION FOR FREEZE-THAW, EROSION
C     AND SOC   
C
C     ALTG,ALT=current,initial elevation of ground surface (m)
C     CDPTH(NUM(NY,NX)-1,=current ground surface elevation (m)
C     ENGYP=cumulative rainfall energy impact on soil surface Ksat 
C        (J m-2)
C     FENGYP=rate constant for restoring surface Ksat (h-1)
C     XNFH=time step from ‘wthr.f’ (h t-1)
C 
      ALTG(NY,NX)=ALT(NY,NX)-CDPTH(NUM(NY,NX)-1,NY,NX)
      ENGYP(NY,NX)=ENGYP(NY,NX)*(1.0-FENGYP*XNFH)
C
C     ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
C     FOR USE AT INTERNAL TIME STEP
C
C     SET INITIAL SNOWPACK VALUES
C
C     VOLS0,VOLSSL=snowpack snow content (water equivalent)
C     VOLI0,VOLISSL=snowpack ice content (m3)
C     VOLW0,VOLWSL=snowpack water content (m3)
C     VOLV0,VOLVSL=snowpack vapor content (m3)
C     VOLS1,VOLSL=snowpack volume (m3)
C     DLYRS0,DLYRS=snowpack depth (m)
C     VHCPWM,VHCPW=snowpack heat capacity (MJ K-1)
C     TK0,TKW=snowpack temperature (K)
C
      DPTHS0(NY,NX)=0.0
      DO 60 L=1,JS
      VOLS0(L,NY,NX)=VOLSSL(L,NY,NX)
      VOLI0(L,NY,NX)=VOLISL(L,NY,NX)
      VOLW0(L,NY,NX)=VOLWSL(L,NY,NX)
      VOLV0(L,NY,NX)=VOLVSL(L,NY,NX)
      VOLS1(L,NY,NX)=VOLSL(L,NY,NX)
      DLYRS0(L,NY,NX)=DLYRS(L,NY,NX)
      VHCPWM(1,L,NY,NX)=VHCPW(L,NY,NX)
      TK0(L,NY,NX)=TKW(L,NY,NX)
      DPTHS0(NY,NX)=DPTHS0(NY,NX)+DLYRS0(L,NY,NX)
60    CONTINUE
C
C     SET SOIL LAYER IN WHICH IRRIGATION IS ADDED
C
C     CDPTH=depth to bottom of soil layer (m)
C     WDPTH,LWDPTH=depth,layer of subsurface irrigation (m)
C
      LWDPTH=NUM(NY,NX)
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
C     FOR USE AT INTERNAL TIME STEP (M=1,NPH)
C
C     PSISM1,PSISM=matric water potential (MPa)
C     VOLA*,VOLV*,VOLW*,VOLI*,VOLP*=pore,vapor,water,ice,air 
C        micropore volumes (m3)
C     VOLWX1=VOLW1 accounting for wetting front
C     VOLP1Z,VOLPH1Z=excess water+ice in micropores,macropores
C        during freezing (if –ve) (m3)
C     VOLAH*,VOLWH*,VOLIH*,VOLPH*=pore,water,ice,air macropores (m3)
C     BKDS=bulk density (Mg m-3) (0=pond)
C     CCLAY=clay concentration (Mg Mg-1)
C     FVOLAH=parameter for clay effect on macropore volume
C     VOLX,VOLT=soil,total volumes (m3)
C     THETW*,THETI*,THETP*=water,ice,air-filled porosity (m3 m-3)
C     VHCP1,VHCM=volumetric heat capacities of total volume, solid 
C        (MJ m-3 K-1) 
C     VHCP1A,VHCP1B=volumetric heat capacities of micropore,macropore
C        (MJ m-3 K-1) 
C
      PSISM1(L,NY,NX)=PSISM(L,NY,NX)
      VOLA1(L,NY,NX)=VOLA(L,NY,NX)
      VOLV1(L,NY,NX)=VOLV(L,NY,NX)
      VOLW1(L,NY,NX)=VOLW(L,NY,NX)
      VOLWX1(L,NY,NX)=VOLWX(L,NY,NX)
      VOLI1(L,NY,NX)=VOLI(L,NY,NX) 
      VOLWH1(L,NY,NX)=VOLWH(L,NY,NX)
      VOLIH1(L,NY,NX)=VOLIH(L,NY,NX)
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
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
      VOLTX(L,NY,NX)=VOLY(L,NY,NX)+VOLAH1(L,NY,NX)
      IF(VOLTX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWX(L,NY,NX)=AMAX1(0.0,(VOLW1(L,NY,NX)+VOLWH1(L,NY,NX))
     2/VOLTX(L,NY,NX))
      THETIX(L,NY,NX)=AMAX1(0.0,(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
     2/VOLTX(L,NY,NX))
      THETPX(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX))
     2/VOLTX(L,NY,NX))
      THETPM(1,L,NY,NX)=AMAX1(0.0,VOLPM(1,L,NY,NX)/VOLTX(L,NY,NX))
      ELSE
      THETWX(L,NY,NX)=POROS(L,NY,NX)
      THETIX(L,NY,NX)=0.0
      THETPX(L,NY,NX)=1.0
      THETPM(1,L,NY,NX)=1.0
      ENDIF
      IF(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX).GT.ZEROS2(NY,NX))THEN 
      THETPY(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX))
     2/(VOLA1(L,NY,NX)+VOLAH1(L,NY,NX)))
      ELSE
      THETPY(L,NY,NX)=0.0
      ENDIF
      VHCP1(L,NY,NX)=VHCM(L,NY,NX)
     2+4.19*(VOLW1(L,NY,NX)+VOLV1(L,NY,NX)+VOLWH1(L,NY,NX))
     2+1.9274*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
      VHCP1A(L,NY,NX)=VHCM(L,NY,NX)
     2+4.19*(VOLW1(L,NY,NX)+VOLW1(L,NY,NX))
     2+1.9274*VOLI1(L,NY,NX)
      VHCP1B(L,NY,NX)=4.19*VOLWH1(L,NY,NX)+1.9274*VOLIH1(L,NY,NX)
C     IF(I.EQ.196.AND.L.EQ.1)THEN
C     WRITE(*,3376)'VOLWI',I,J,NFZ,NX,NY,L,VOLW1(L,NY,NX)
C    2,VOLI1(L,NY,NX),VOLP1(L,NY,NX),VOLA1(L,NY,NX)
C    3,VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),VOLPH1(L,NY,NX) 
C    3,VOLAH1(L,NY,NX),VOLT(L,NY,NX),VOLY(L,NY,NX) 
C    4,THETWX(L,NY,NX),THETIX(L,NY,NX),THETPX(L,NY,NX)
C    5,VHCM(L,NY,NX),VHCP1(L,NY,NX)
3376  FORMAT(A8,6I4,40E14.6)
C     ENDIF
C
C     INITIALIZE MACROPOROSITY, TEMPERATURES, INCOMING WATER 
C     AND HEAT FLUXES
C
C     VOLAH1=total macropore volume (m3)
C     FMAC,FGRD=macropore,micropore volume fractions
C     CNDH*=macropore hydraulic conductivity (m2 h-1 MPa-1)
C     TKAM=air temperature (K)
C     TKS,TK1=soil temperature (K)
C     FLU,HWFLU=subsurface water,convective heat fluxes from
C        irrigation file (m3 t-1,MJ t-1)
C     AREAU,AREAD=area of layer below natural,artificial 
C        water table
C     XNPH=time step from ‘wthr.f’(h t-1)
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
      TKSM(1,L,NY,NX)=TKS(L,NY,NX)
      IF(L.EQ.LWDPTH)THEN
      FLU(L,NY,NX)=PRECU(NY,NX)*XNFH
      HWFLU(L,NY,NX)=4.19*TKAM(NY,NX)*FLU(L,NY,NX)
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
     2,(CDPTH(L,NY,NX)-DTBLY(NY,NX))/DLYR(3,L,NY,NX)))
      ELSE
      AREAUD(L,NY,NX)=0.0
      ENDIF
30    CONTINUE
C
C     ENTER STATE VARIABLES AND DRIVERS INTO LOCAL ARRAYS
C     FOR USE AT INTERNAL TIME STEP IN SURFACE LITTER
C
C     VHCP1=volumetric heat capacity of litter (MJ K-1)
C     ORGC,ORGCC=litter organic C, charcoal (g)
C     ALBZ=dry surface litter albedo
C     VOLA*,VOLV*,VOLW*,VOLI*,VOLP*=pore,vapor,water,ice,air 
C        volumes in litter (m3)
C     THETW*,THETI*,THETP*=litter water,ice,air concentrations (m3 m-3)
C     VOLR=litter volume (m3)
C     VOLWRX=litter water retention capacity (m3)
C     XVOLT,XVOLW=surface water+ice,surface water in excess of VOLWRX
C        (m3) 
C     VHCPRX=minimum heat capacity for solving litter water,heat fluxes 
C        (MJ K-1)
C     VOLR=litter volume (m3)
C     PSISM*=litter matric water potential (MPa)
C     TK*=litter temperature (K)
C
      VHCP1(0,NY,NX)=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))
     2+4.19*(VOLW(0,NY,NX)+VOLV(0,NY,NX))
     2+1.9274*VOLI(0,NY,NX)
      IF(ORGC(0,NY,NX).GT.ZEROS(NY,NX))THEN
      ALBZ(NY,NX)=(0.30*ORGC(0,NY,NX)+0.00*ORGCC(0,NY,NX))
     2/(ORGC(0,NY,NX)+ORGCC(0,NY,NX))
      ELSE
      ALBZ(NY,NX)=0.30
      ENDIF
      VOLA1(0,NY,NX)=VOLA(0,NY,NX)
      VOLW1(0,NY,NX)=VOLW(0,NY,NX)
      VOLV1(0,NY,NX)=VOLV(0,NY,NX)
      VOLI1(0,NY,NX)=VOLI(0,NY,NX)
      VOLP1(0,NY,NX)=AMAX1(0.0,VOLA1(0,NY,NX)-VOLW1(0,NY,NX)
     2-VOLI1(0,NY,NX))
      VOLWM(1,0,NY,NX)=VOLW1(0,NY,NX)
      VOLPM(1,0,NY,NX)=VOLP1(0,NY,NX)
      TVOLWI=VOLW(0,NY,NX)+VOLI(0,NY,NX)
      XVOLT(NY,NX)=AMAX1(0.0,TVOLWI-VOLWRX(NY,NX))
      IF(VOLR(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWX(0,NY,NX)=AMAX1(0.0,VOLW1(0,NY,NX)/VOLR(NY,NX))
      THETIX(0,NY,NX)=AMAX1(0.0,VOLI1(0,NY,NX)/VOLR(NY,NX))
      THETPX(0,NY,NX)=AMAX1(0.0,VOLP1(0,NY,NX)/VOLR(NY,NX))
     2*AMAX1(0.0,(1.0-XVOLT(NY,NX)/VOLWD(NY,NX)))
      THETPM(1,0,NY,NX)=AMAX1(0.0,VOLPM(1,0,NY,NX)/VOLR(NY,NX))
     2*AMAX1(0.0,(1.0-XVOLT(NY,NX)/VOLWD(NY,NX)))
      ELSE
      THETWX(0,NY,NX)=0.0
      THETIX(0,NY,NX)=0.0
      THETPX(0,NY,NX)=1.0
      THETPM(1,0,NY,NX)=1.0
      ENDIF
      PSISM1(0,NY,NX)=PSISM(0,NY,NX)
      TK1(0,NY,NX)=TKS(0,NY,NX)
      TKSM(1,0,NY,NX)=TKS(0,NY,NX)
      TKQG(1,NY,NX)=TKQGX(NY,NX)
      VPQG(1,NY,NX)=VPQGX(NY,NX)
C
C     SURFACE COVER BY SNOWPACK, LITTER
C
C     FSNW,FSNX=snow,snow-free cover fractions 
C     DPTHSX=minimum snowpack depth for full cover (m)
C     BARE,CVRD=soil,litter cover fractions 
C     BAREW,CVRDW=soil,litter cover fractions accounting for excess
C        surface water
C     TKGS=composite surface temperature (K) 
C
      FSNW(NY,NX)=AMIN1(1.0,SQRT(DPTHS(NY,NX)/DPTHSX))
      IF(FSNW(NY,NX).LT.1.0)THEN
      FSNX(NY,NX)=AMAX1(1.0E-03,1.0-FSNW(NY,NX))
      FSNW(NY,NX)=1.0-FSNX(NY,NX)
      ELSE
      FSNX(NY,NX)=0.0
      ENDIF
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      BARE(NY,NX)=AMIN1(1.0,AMAX1(0.0,EXP(-0.5E-02
     2*((ORGC(0,NY,NX)+ORGCC(0,NY,NX))/AREA(3,0,NY,NX)))))
      ELSE
      BARE(NY,NX)=1.0
      ENDIF
      CVRD(NY,NX)=1.0-BARE(NY,NX)
      BAREW(NY,NX)=AMAX1(0.0,BARE(NY,NX)
     2-AMIN1(1.0,AMAX1(0.0,XVOLT(NY,NX)/VOLWD(NY,NX))))
      CVRDW(NY,NX)=1.0-BAREW(NY,NX)
      TKGS(1,NY,NX)=FSNW(NY,NX)*TKW(1,NY,NX)+FSNX(NY,NX)
     2*(TKS(0,NY,NX)*CVRDW(NY,NX)+TKS(NUM(NY,NX),NY,NX)*BAREW(NY,NX))
C     IF(I.GT.270)THEN
C     WRITE(*,7751)'THETPX',I,J,NFZ,NX,NY
C    2,VOLW1(0,NY,NX),VOLI1(0,NY,NX),VOLP1(0,NY,NX),VOLA1(0,NY,NX)
C    3,THETWX(0,NY,NX),THETIX(0,NY,NX),THETPX(0,NY,NX) 
C    3,XVOLW(NY,NX),XVOLT(NY,NX),VOLWD(NY,NX),VOLWG(NY,NX),VOLR(NY,NX) 
C    4,VOLW1(1,NY,NX),VOLI1(1,NY,NX),VOLP1(1,NY,NX)
C    5,VOLWRX(NY,NX),TKS(0,NY,NX)
C    6,TKQG(1,NY,NX),VPQG(1,NY,NX)
C    2,TKGS(1,NY,NX),FSNW(NY,NX),TKW(1,NY,NX),FSNX(NY,NX)
C    2,TKS(0,NY,NX),CVRDW(NY,NX),TKS(NUM(NY,NX),NY,NX),BAREW(NY,NX)
C    1,BARE(NY,NX),XVOLT(NY,NX),VOLWD(NY,NX),VOLW(0,NY,NX)
C    2,VOLI(0,NY,NX),VOLWRX(NY,NX)
7751  FORMAT(A8,5I4,20E12.4)
C     ENDIF
C
C     PRECIPITATION AT SOIL SURFACE USED TO CALCULATE WATER EROSION
C     IN ‘EROSION.F’
C
C     PRECM,PRECA=precipitation+irrigation at top,bottom of canopy 
C        (mm h-1)
C     TFLWC,TFLWCI=total water retention,interception by canopy 
C        from ‘hour1.f’(m h-1)
C     PRECD,PRECB=direct,indirect precipn+irrign at soil surface 
C        (mm h-1) 
C
      PRECM(NY,NX)=1.0E+03*PRECA(NY,NX)
     2/AREA(3,NU(NY,NX),NY,NX)
      PRECD(NY,NX)=1.0E+03*(PRECA(NY,NX)-TFLWCI(NY,NX))
     2/AREA(3,NU(NY,NX),NY,NX)
      PRECB(NY,NX)=1.0E+03*(TFLWCI(NY,NX)-TFLWC(NY,NX))
     2/AREA(3,NU(NY,NX),NY,NX)
C     IF(ICHKF.EQ.1)THEN
C     WRITE(*,3112)'BARE',I,J,NFZ,NX,NY,BARE(NY,NX),BAREY
C    2,FSNX(NY,NX),ORGC(0,NY,NX),ORGCC(0,NY,NX),VOLWRX(NY,NX)
C    3,XVOLW(NY,NX),VOLWD(NY,NX)
C    4,PRECA(NY,NX),TFLWCI(NY,NX),TFLWC(NY,NX)
C    5,PRECA(NY,NX)*XNPHX*1000*BARE(NY,NX),PRECD(NY,NX),PRECB(NY,NX)
3112  FORMAT(A8,5I4,20E12.4)
C     ENDIF 
C
C     LITTER DEPTH
C
C     DLYRR=surface litter depth (m)
C
      DLYRR(NY,NX)=AMAX1(1.0E-06,DLYR(3,0,NY,NX))
C
C     DISTRIBUTION OF PRECIPITATION AND ITS HEAT AMONG SURFACE
C     RESIDUE, SOIL SURFACE, AND MACROPORES
C
C     PRECA,PRECW=rainfall+irrigation,snowfall (water equiv) from
C        ‘wthr.f’(m3)
C     FLWQW=rainfall to snowpack (m3 h-1)
C     FLWSW=snowfall to snowpack (m3 h-1)
C     HFLWSW=convective heat flux to snowpack (MJ h-1)
C     TKAM=air temperature (K)
C     FLWQB=precip to litter+soil surfaces (m3 h-1) 
C     FLWQAX,FLWQBX=precip to soil,litter surfaces (m3 h-1)
C     HFLWQA,HFLWQB=convective heat flux to soil,litter surfaces 
C        (MJ h-1)
C     FLWQAS,FLWQAH=precip to soil micropores,macropores (m3 h-1)
C
      IF(PRECA(NY,NX).GT.0.0.OR.PRECW(NY,NX).GT.0.0)THEN
      FLWQW=(PRECA(NY,NX)-TFLWC(NY,NX))*FSNW(NY,NX) 
      FLWSW=PRECW(NY,NX)
      HFLWSW=2.095*TKAM(NY,NX)*FLWSW+4.19*TKAM(NY,NX)*FLWQW
      FLWQB=(PRECA(NY,NX)-TFLWC(NY,NX))*FSNX(NY,NX)
      IF(BKDS(NUI(NY,NX),NY,NX).LE.ZERO
     2.AND.CDPTH(NU(NY,NX)-1,NY,NX).LE.CDPTHI(NY,NX))THEN
      FLWQBX=FLWQB
      FLWQAX=0.0
      ELSE     
      FLWQBX=AMIN1((VOLWRX(NY,NX)-VOLW1(0,NY,NX))*FSNX(NY,NX)*XNFH
     2,FLWQB*CVRD(NY,NX)) 
      FLWQAX=FLWQB-FLWQBX 
      ENDIF
      HFLWQB=4.19*TKAM(NY,NX)*FLWQBX
      HFLWQA=4.19*TKAM(NY,NX)*FLWQAX
      FLWQAS=FLWQAX*FGRD(NUM(NY,NX),NY,NX)
      FLWQAH=FLWQAX*FMAC(NUM(NY,NX),NY,NX)
      ELSE
      FLWQW=-TFLWC(NY,NX)*FSNW(NY,NX) 
      FLWSW=0.0
      HFLWSW=4.19*TKAM(NY,NX)*FLWQW
      FLWQB=-TFLWC(NY,NX)*FSNX(NY,NX)
      IF(BKDS(NUI(NY,NX),NY,NX).LE.ZERO
     2.AND.CDPTH(NU(NY,NX)-1,NY,NX).LE.CDPTHI(NY,NX))THEN
      FLWQBX=FLWQB
      FLWQAX=0.0
      ELSE     
      FLWQBX=AMIN1((VOLWRX(NY,NX)-VOLW1(0,NY,NX))*FSNX(NY,NX)*XNFH
     2,FLWQB*CVRD(NY,NX)) 
      FLWQAX=FLWQB-FLWQBX 
      ENDIF
      HFLWQB=4.19*TKAM(NY,NX)*FLWQBX
      HFLWQA=4.19*TKAM(NY,NX)*FLWQAX
      FLWQAS=FLWQAX*FGRD(NUM(NY,NX),NY,NX)
      FLWQAH=FLWQAX*FMAC(NUM(NY,NX),NY,NX)
      ENDIF
C
C     PRECIP ON SNOW ARRAYS EXPORTED TO ‘TRNSFR.F’, ‘TRNSFRS.F’ 
C     FOR SOLUTE FLUX CALCULATIONS AT NFH TIME STEP
C
C     PRECW,PRECR,PRECQ,PRECI=snow,rain,snow+rain,irrigation (m3 t-1)
C     VHCPW,VHCPWX=current, minimum snowpack heat capacities 
C        (MJ m-3 K-1)
C     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
C        (m3 t-1)
C     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation (m3 t-1)
C     XNFH=time step from ‘wthr.f’(h t-1)
C
      IF(PRECW(NY,NX).GT.0.0.OR.(PRECR(NY,NX).GT.0.0
     2.AND.VHCPW(1,NY,NX).GT.VHCPWX(NY,NX)))THEN
      FLQRQ(NY,NX)=0.0
      FLQRI(NY,NX)=0.0
      FLQGQ(NY,NX)=PRECQ(NY,NX)*XNFH
      FLQGI(NY,NX)=PRECI(NY,NX)*XNFH
      ELSEIF((PRECQ(NY,NX).GT.0.0.OR.PRECI(NY,NX).GT.0.0)
     2.AND.VHCPW(1,NY,NX).LE.VHCPWX(NY,NX))THEN
      FLQRQ(NY,NX)=FLWQBX*PRECQ(NY,NX)
     2/(PRECQ(NY,NX)+PRECI(NY,NX))*XNFH
      FLQRI(NY,NX)=FLWQBX*PRECI(NY,NX)
     2/(PRECQ(NY,NX)+PRECI(NY,NX))*XNFH
      FLQGQ(NY,NX)=(PRECQ(NY,NX)*XNFH-FLQRQ(NY,NX)) 
      FLQGI(NY,NX)=(PRECI(NY,NX)*XNFH-FLQRI(NY,NX))
C     WRITE(*,7763)'FLQRI',I,J,NFZ,FLQRI(NY,NX),FLQRQ(NY,NX)
C    2,FLQGI(NY,NX),FLQGQ(NY,NX),FLWQBX,PRECI(NY,NX),PRECQ(NY,NX)
7763  FORMAT(A8,3I4,20E12.4) 
      ELSE
      FLQRQ(NY,NX)=0.0
      FLQRI(NY,NX)=0.0
      FLQGQ(NY,NX)=0.0
      FLQGI(NY,NX)=0.0
      ENDIF
C
C     GATHER PRECIPITATION AND MELTWATER FLUXES AND THEIR HEATS
C     AMONG ATMOSPHERE, SNOWPACK, RESIDUE AND SOIL SURFACES
C     INTO LOCAL ARRAYS AT CURRENT TIME STEP FOR USE IN MASS AND ENERGY
C     EXCHANGE ALGORITHMS
C     
C     FLW0S,FLQ0I,FLQ0W=snow,ice,water input to snowpack (m3 t-1)
C     HWFLQ0=convective heat flux to snowpack (MJ t-1)
C     FLQ1,FLH1,FLY1=rain+irrigation to micropores,macropores,litter
C        (m3 t-1)
C     HWFLQ1,HWFLY1=convective heat flux to soil,litter surfaces 
C        (MJ t-1)
C     XNPHX=time step for water fluxes from ‘wthr.f’(h t-1)
C
      FLQ0S(NY,NX)=FLWSW*XNPHX 
      FLQ0I(NY,NX)=0.0
      FLQ0W(NY,NX)=FLWQW*XNPHX 
      HWFLQ0(NY,NX)=HFLWSW*XNPHX 
      FLQ1(NY,NX)=FLWQAS*XNPHX
      FLH1(NY,NX)=FLWQAH*XNPHX
      FLY1(NY,NX)=FLWQBX*XNPHX
      HWFLQ1(NY,NX)=HFLWQA*XNPHX
      HWFLY1(NY,NX)=HFLWQB*XNPHX
C
C     INITIALIZE PARAMETERS, FLUXES FOR ENERGY EXCHANGE
C     AT SNOW, RESIDUE AND SOIL SURFACES
C
C     RADG=shortwave radiation at ground surface from ‘hour1.f’(MJ h-1)
C     RADGX=shortwave radiation at ground surface (MJ t-1)
C     RADXW,RADXG,RADXR= shortwave radn at snowpack,soil,litter 
C        (MJ t-1) 
C     FRADG=fraction of shortwave radiation at ground surface
C     FSNW,FSNX=fractions of snow,snow-free cover
C     BARE,CVRD=fractions of soil,litter cover
C     XNPS=internal time step for fluxes through snowpack 
C        from ‘wthr.f’(h t-1)
C     THS=sky LW radiation from ‘wthr.f’(MJ h-1)
C     THRYX=longwave radiation at ground surface (MJ t-1)
C
      RADGX=RADG(NY,NX)*XNPHX 
      RADXW(NY,NX)=RADGX*FSNW(NY,NX)*XNPS
      RADXG(NY,NX)=RADGX*FSNX(NY,NX)*BAREW(NY,NX)
      RADXR(NY,NX)=RADGX*FSNX(NY,NX)*CVRDW(NY,NX)*XNPR
      THRYX=THS(NY,NX)*XNPHX*FRADG(NY,NX)
C
C     FIRE IGNITION
C
C     ITILL=22:fire event from disturbance file
C     ZNOON=hour of solar noon
C     THRYX=sky + fire longwave radiation at ground surface (MJ t-1)
C     DCORP=fire ignition intensity from disturbance file (kW m-2)
C     XNFZ,XNFH,XNPH,XNPHX=time steps for water,heat fluxes from
C        ‘wthr.f’(h t-1)
C     HCBFH=heat released in canopy air by combustion in previous time
C        step (MJ t-1)
C     HCBFX=heat released in soil layer by combustion in previous time
C        step (MJ t-1)
C     HFLXH,HFLXF=heat from HCBFH,HCBFX added in current time step 
C        (MJ t-1)
C 
      IF(ITILL(I,NY,NX).EQ.22)THEN
      IF(J.EQ.INT(ZNOON(NY,NX)))THEN
      THRYX=THRYX+3.6*DCORP(I,NY,NX)*AMIN1(1.0,4.0*XNFZ*XNFH)
     2*AREA(3,NU(NY,NX),NY,NX)*XNPHX
      ELSEIF(J.EQ.INT(ZNOON(NY,NX)+1))THEN
      THRYX=THRYX+3.6*DCORP(I,NY,NX)*AMAX1(0.0,1.0-4.0*XNFZ*XNFH)
     2*AREA(3,NU(NY,NX),NY,NX)*XNPHX
      ENDIF
C     WRITE(*,1109)'THRYX',I,J,NFZ,NX,NY
C    2,THRYX,DCORP(I,NY,NX),XNPHX
C    3,AMIN1(1.0,4.0*XNFZ*XNFH),AMAX1(0.0,1.0-4.0*XNFZ*XNFH)
C    3,THS(NY,NX),FRADG(NY,NX) 
1109  FORMAT(A8,5I4,10E12.4)
      ENDIF
      HFLXH(NY,NX)=HCBFH(NY,NX)*XNPH
      HFLXF(0,NY,NX)=HCBFX(0,NY,NX)*XNPH
      DO 975 L=NUM(NY,NX),NL(NY,NX)
      HFLXF(L,NY,NX)=HCBFX(L,NY,NX)*XNPH
975   CONTINUE
C
C     END FIRE IGNITION
C
C     THRYW,THRYG,THRYR=longwave radiationn incident at
C        snowpack,soil,litter (MJ t-1)
C     EMMW,EMMS,EMMR=emissivity of snowpack,soil surface,litter
C     EMMGW,EMMCW=emission of snowpack,canopy over snowpack 
C     EMMGS,EMMCG=emission of soil surface,canopy over soil surface 
C     EMMGR,EMMCR=emission of litter,canopy over litter
C     BAREW,CVRDW=soil,litter cover fractions accounting for excess
C        surface water
C     XNP*=time step from ‘wthr.f’ (h t-1)
C
      THRYW(NY,NX)=THRYX*FSNW(NY,NX)*XNPS
      THRYG(NY,NX)=THRYX*FSNX(NY,NX)*BAREW(NY,NX)
      THRYR(NY,NX)=THRYX*FSNX(NY,NX)*CVRDW(NY,NX)*XNPR
      EMMGW(NY,NX)=EMMW*2.04E-10*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNW(NY,NX)*XNPYX*FRADG(NY,NX)
      EMMCW(NY,NX)=EMMW*2.04E-10*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNW(NY,NX)*XNPYX
      EMMGS(NY,NX)=EMMS*2.04E-10*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNX(NY,NX)*BAREW(NY,NX)*XNPHX*FRADG(NY,NX) 
      EMMCS(NY,NX)=EMMS*2.04E-10*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNX(NY,NX)*BAREW(NY,NX)*XNPHX 
      EMMGR(NY,NX)=EMMR*2.04E-10*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNX(NY,NX)*CVRDW(NY,NX)*XNPZX*FRADG(NY,NX) 
      EMMCR(NY,NX)=EMMR*2.04E-10*AREA(3,NUM(NY,NX),NY,NX)
     2*FSNX(NY,NX)*CVRDW(NY,NX)*XNPZX 
C
C     AERODYNAMIC RESISTANCE OF SNOWPACK, RESIDUE AND SOIL
C     SURFACES TO ENERGY EXCHANGE WITH ATMOSPHERE
C     Soil Sci. Soc. Am. J. 48:25-32
C
C     RZR=porosity-unlimited litter boundary layer resistance (h m-1)
C     DLYRR=litter depth (m)
C     WGSGR=vapor diffusivity in litter from ‘hour1.f’ (m2 h-1)
C
      RZR(NY,NX)=DLYRR(NY,NX)/WGSGR(NY,NX)
C
C     SNOWPACK DIFFUSIVITY
C
C     RAS=snowpack aerodynamic resistance (h m-1)
C     VOLS,VOLS1=snowpack total,layer volume (m3)
C     THETPL=snowpack layer air-filled porosity (m3 m-3)
C     VOLS0,VOLI0,VOLW0=snowpack snow,ice,water volume (m3)
C
      RAS(NY,NX)=0.0
      IF(VOLS(NY,NX).GT.ZEROS2(NY,NX))THEN
      DO 9775 L=1,JS
      IF(VOLS1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      RASX=DLYRS(L,NY,NX)/WGSGW(L,NY,NX)
      THETPL=AMAX1(THETPI,1.0-(VOLS0(L,NY,NX)+VOLI0(L,NY,NX)
     2+VOLW0(L,NY,NX))/VOLS1(L,NY,NX))
      RASL=RASX/AMAX1(ZERO,THETPL)**2.0
      RAS(NY,NX)=RAS(NY,NX)+RASL
C     WRITE(*,3113)'RAS',I,J,NFZ,NX,NY,L,RAS(NY,NX),RASL,RASX
C    2,DLYRS(L,NY,NX),WGSGW(L,NY,NX),THETPL,THETPI,VOLS0(L,NY,NX)
C    3,VOLI0(L,NY,NX),VOLW0(L,NY,NX),VOLS1(L,NY,NX),TKW(L,NY,NX)
3113  FORMAT(A8,6I4,40E12.4)
      ENDIF
9775  CONTINUE
      ENDIF
C     IF(NX.EQ.1)THEN
C     WRITE(*,3111)'RAC',I,J,NX,NY,RAC(NY,NX)
C    2,ZT(NY,NX),RAS(NY,NX),VOLS(NY,NX) 
C    3,DLYRR(NY,NX),THETX,THETPX(0,NY,NX),VHCP1(0,NY,NX)
C    5,WGSGR(NY,NX),VOLW1(0,NY,NX)
C    5,VOLI1(0,NY,NX),VOLP1(0,NY,NX),VOLR(NY,NX),VOLA1(0,NY,NX)
C    4,TEVX(NY,NX),TSHX(NY,NX),RADG(NY,NX),THS(NY,NX)
C    5,FRADG(NY,NX),ZS(NY,NX)
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
      N1=NX
      N2=NY
      DO 35 L=NUM(NY,NX),NL(NY,NX)
      N3=L
      DO 40 N=NCN(NY,NX),3
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
C     CNDH1=macropore hydraulic conductivity (m2 h MPa-1)
C     AVCNHL=macropore hydraulic conductance (m h MPa-1)
C     DLYR=layer depth (m)
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
      EVAPG(NY,NX)=0.0
      EVAPR(NY,NX)=0.0
      EVAP0(NY,NX)=0.0
      WFLVR(NY,NX)=0.0
      HFLVR(NY,NX)=0.0
      WFLFR(NY,NX)=0.0
      HFLFR(NY,NX)=0.0
      WFLVL(NUM(NY,NX),NY,NX)=0.0
      HFLVL(NUM(NY,NX),NY,NX)=0.0
      WFLFL(NUM(NY,NX),NY,NX)=0.0
      HFLFL(NUM(NY,NX),NY,NX)=0.0
      FLQRM(M,NY,NX)=0.0
      FLQSM(M,NY,NX)=0.0
      FLQHM(M,NY,NX)=0.0
      TQR1(NY,NX)=0.0
      THQR1(NY,NX)=0.0
      TQS1(NY,NX)=0.0
      TQW1(NY,NX)=0.0
      TQI1(NY,NX)=0.0
      THQS1(NY,NX)=0.0
      FLWRM(M,NY,NX)=0.0
C
C     CALCULATE CANOPY AIR TEMPERATURE, VAPOR CONCENTRATION
C
C     DTKQ=atmosphere-ground surface air temperature difference (K)
C     TKAM=air temperature (K)
C     TKQG=air temperature at ground surface (K)
C     RI=Richardson number
C     RAB,RABX=aerodynamic,isothermal (from ‘hour1.f’) boundary 
C        layer resistance (h m-1)
C     ZT,ZS=canopy,soil surface roughness height
C        from ‘hour1.f’ (m)
C     RACG,RACGX=canopy aerodynamic resistance below 
C        maximum canopy height (h m-1)
C     RATG=aerodynamic+canopy boundary layer resistance (h m-1)
C     DTKG=ground surface air-ground surface temperature difference (K)
C     ARLSS=total leaf,stalk,standing dead of all canopies (m2)
C     ZT=canopy height (m)
C     AREA=area of grid cell (m2)
C     RARH,RAH=current,isothermal ground surface resistance (h m-1)
C     RARW,RAHW=current,isothermal snowpack surface resistance (h m-1)
C     THETP0,THETPG=litter air-filled porosity without,with free
C        surface water (m3 m-3)
C     RZR=porosity-unlimited litter boundary layer resistance (h m-1)
C     DFVR,DFVG=porosity effects on litter boundary layer resistance
C        without,with free surface water
C     POROS,POROQ=litter porosity (m3 m-3), tortuosity from ‘starts.f’ 
C     RZRX,RZRG=litter vapor resistance without,with free surface water
C        (h m-1)  
C     RAGS=soil surface boundary layer resistance (h m-1)
C
      DTKQ=TKAM(NY,NX)-TKQG(M,NY,NX) 
      RI=AMAX1(RIX,AMIN1(RIY
     2,RIBX(NY,NX)/TKAM(NY,NX)*DTKQ))
      RIC=1.0-10.0*RI
      RAB(NY,NX)=AMIN1(RABZ,AMAX1(RABM,RABX(NY,NX)/RIC))
      IF(ZT(NY,NX).GT.ZS(NY,NX))THEN
      ARDNS=ARLSS(NY,NX)/(ZT(NY,NX)*AREA(3,0,NY,NX))
      RACGX=ZT(NY,NX)*ARDNS*2.0E-04/WGSGA(NY,NX) 
      ELSE
      RACGX=0.0
      ENDIF
      RACG(M,NY,NX)=AMIN1(RACGZ,AMAX1(RACGM,RACGX/RIC))
C     WRITE(*,101)'RACG',I,J,NFZ,M,NY,NX,RACGX
C    2,RACG(M,NY,NX),ZT(NY,NX),ARDNS,ARLSS(NY,NX),WGSGA(NY,NX),RIC
101   FORMAT(A8,6I4,20E12.4) 
      RATG=RAB(NY,NX)+RACG(M,NY,NX)
      DTKG=TKQG(M,NY,NX)-TKGS(M,NY,NX)
      RI=AMAX1(RIX,AMIN1(RIY
     2,RIBX(NY,NX)/TKQG(M,NY,NX)*DTKG))
      RIG=1.0-3.2*RI
      RARH=AMIN1(RABZ,AMAX1(RABM,RAH/RIG))
      RARW=AMIN1(RABZ,AMAX1(RABM,RAHW/RIG))
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX))THEN
      THETP0=AMAX1(0.0,VOLP1(0,NY,NX)/VOLR(NY,NX))
      THETPG=THETP0*AMAX1(0.0,(1.0-XVOLT(NY,NX)/VOLWD(NY,NX)))
      DFVR=AMAX1(ZERO2,POROQ*THETP0**2/POROS(0,NY,NX)) 
      DFVG=AMAX1(ZERO2,POROQ*THETPG**2/POROS(0,NY,NX)) 
      RZRX=RZR(NY,NX)/DFVR
      RZRE=RARH+0.25*RZRX*THETP0
      RZRG=RZR(NY,NX)/DFVG*CVRD(NY,NX)
      ELSE
      RZRE=RARH
      RZRG=RARH
      ENDIF
      RAGS=RARH+RZRG
C     WRITE(*,4424)'RAGS',I,J,NFZ,M,NY,NX
C    2,RAGS,RARH,RZRG,RZR(NY,NX),DFVG,VOLA1(0,NY,NX),VOLW(0,NY,NX)
C    2,VOLI(0,NY,NX),CVRD(NY,NX)
C    3,VOLR(NY,NX),DLYRR(NY,NX),THETPX(0,NY,NX),THETP0,VOLP1(0,NY,NX)
C    2,POROS(0,NY,NX),VOLX(0,NY,NX),XVOLT(NY,NX),VOLWD(NY,NX)
4424  FORMAT(A8,6I4,20E12.4)
C
C     SURFACE CONDUCTANCES FOR LATENT,SENSIBLE HEAT FLUXES
C
C     FSNW,FSNX=snow,snow-free cover fraction
C     BAREW,CVRDW=soil,litter cover fractions accounting for excess
C        surface water
C     VHCPRX=minimum heat capacity for solving litter water,heat fluxes 
C        (MJ K-1)
C     PAREX,PARSX=terms used to calculate boundary layer
C        conductance for latent,sensible heat from ‘hour1.f’ 
C        (m2 h t-1,MJ h m-1 K-1 t-1)
C     PAREWM,PARSWM=conductances for snowpack latent,sensible
C        heat fluxes (m3 t-1,MJ K-1 t-1) 
C     PAREGM,PARSGM=conductances for soil latent,sensible 
C        heat fluxes (m3 t-1,MJ K-1 t-1) 
C     PARERM,PARSRM=conductances for litter latent,sensible 
C        heat fluxes (m3 t-1,MJ K-1 t-1) 
C     RARH=current surface layer resistance (h m-1)
C     RAGS=soil surface boundary layer resistance (h m-1)
C     RAE=resistance to evaporation at soil and litter surfaces (h m-1)
C     XNPS,XNPR=time step  for fluxes from ‘wthr.f’(h t-1)
C
      FSNW(NY,NX)=AMIN1(1.0,SQRT(DPTHS0(NY,NX)/DPTHSX))
      IF(FSNW(NY,NX).LT.1.0)THEN
      FSNX(NY,NX)=AMAX1(1.0E-03,1.0-FSNW(NY,NX))
      FSNW(NY,NX)=1.0-FSNX(NY,NX)
      ELSE
      FSNX(NY,NX)=0.0
      ENDIF
      VHCPRXS(NY,NX)=VHCPRX(NY,NX)*FSNW(NY,NX)
      VHCPRXG(NY,NX)=VHCPRX(NY,NX)*FSNX(NY,NX)
      BAREW(NY,NX)=AMAX1(0.0,BARE(NY,NX)
     2-AMIN1(1.0,AMAX1(0.0,XVOLT(NY,NX)/VOLWD(NY,NX))))
      CVRDW(NY,NX)=1.0-BAREW(NY,NX)
      PAREW=PAREX(NY,NX)*FSNW(NY,NX)*XNPS
      PARSW=PARSX(NY,NX)*FSNW(NY,NX)*XNPS 
      PAREG=PAREX(NY,NX)*FSNX(NY,NX)*BAREW(NY,NX) 
      PARSG=PARSX(NY,NX)*FSNX(NY,NX)*BAREW(NY,NX) 
      PARER=PAREX(NY,NX)*FSNX(NY,NX)*CVRDW(NY,NX)*XNPR
      PARSR=PARSX(NY,NX)*FSNX(NY,NX)*CVRDW(NY,NX)*XNPR
C     WRITE(*,3115)'RZR',I,J,NFZ,M,NX,NY
C    2,RZR(NY,NX),RZRX,DFVR
C    2,THETPX(0,NY,NX),VOLP1(0,NY,NX),VOLR(NY,NX),DLYRR(NY,NX)
C    3,PARER,PAREG,FSNW(NY,NX),FSNX(NY,NX),XNPR
C    4,BARE(NY,NX),BAREW(NY,NX),CVRD(NY,NX),CVRDW(NY,NX)
C    5,VHCP1(0,NY,NX),ORGC(0,NY,NX),ORGCC(0,NY,NX)
C    6,XVOLT(NY,NX),VOLWD(NY,NX),VOLW(0,NY,NX),VOLWRX(NY,NX)
C    7,DPTHS0(NY,NX),DPTHSX 
3115  FORMAT(A8,6I4,40E12.4)
      PAREWM=PAREW/(RARW+RAE)
      PARSWM=PARSW/RARW 
      PAREGM=PAREG/(RAGS+RAE)             
      PARSGM=PARSG/RAGS
      PARERM=PARER/(RZRE+RAE)
      PARSRM=PARSR/RARH
C     IF(J.EQ.14)THEN
C     WRITE(*,4423)'RAB',I,J,NFZ,M,NY,NX,PAREGM,PARSGM,PARERM,PARSRM 
C    2,PAREG,PARSG,PARER,PARSR
C    2,RAH,RAE,RZRE,RARH,RZRX,RZR(NY,NX),RAGS,CVRDW(NY,NX),THETP0
C    2,RAB(NY,NX),RABM,RABX(NY,NX),RI
C    2,RI,RIX,RIY,RIBX(NY,NX),TKGS(M,NY,NX),TKQG(M,NY,NX) 
4423  FORMAT(A8,6I4,30E12.4)
C     ENDIF 
C
C     REDISTRIBUTE INCOMING PRECIPITATION
C     BETWEEN RESIDUE AND SOIL SURFACE
C
C     BKDS=bulk density (0=pond) (Mg m-3)
C     FLQRS,FLQRH=water flux from soil micropores,macropores to litter
C        (m3 t-1)
C     FLQ1,FLH1,FLY1=rain+irrigation to micropores,macropores,litter
C        (m3 t-1)
C     VOLP1,VOLPH1=air-filled microporosity,macroporosity (m3)
C     HFLQR1=convective heat flux from soil to litter (MJ t-1)
C     FLYM,HWFLYM=total water flux, convective heat flux to litter
C        (m3 t-1,MJ t-1))
C     FLQM,FLHM=total water flux to soil micropores, macropores
C        (m3 t-1)
C     HWFLQM=total convective heat flux to soil micropores, macropores
C        (MJ t-1)
C     XNPR=time step for litter water,heat fluxes from ‘wthr.f’(h t-1)
C
      IF(BKDS(NUM(NY,NX),NY,NX).GT.ZERO)THEN
      FLQRS=AMAX1(0.0,FLQ1(NY,NX)-VOLP1(NUM(NY,NX),NY,NX)*FSNX(NY,NX))
      FLQRH=AMAX1(0.0,FLH1(NY,NX)-VOLPH1(NUM(NY,NX),NY,NX)*FSNX(NY,NX))
      HFLQR1=4.19*TKAM(NY,NX)*(FLQRS+FLQRH)
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
C     IF(I.EQ.178)THEN
C     WRITE(*,4422)'FLQ0W',I,J,NFZ,M,FLQ0W(NY,NX),FLWQW,XNPH
C     WRITE(*,4422)'FLY',I,J,NFZ,M,PRECA(NY,NX),TFLWC(NY,NX) 
C    2,FLY1(NY,NX),FLQ1(NY,NX),FLH1(NY,NX),FLYM,FLQM,FLHM 
C    3,FGRD(NUM(NY,NX),NY,NX),CVRD(NY,NX)
C    4,FHOL(L,NY,NX),VOLAH1(L,NY,NX) 
C    5,FLWQAX,PRECA(NY,NX),TFLWC(NY,NX),FLWQBX
C    6,BARE(NY,NX),ORGC(0,NY,NX),XVOLW(NY,NX),VOLWG(NY,NX)
C    7,VOLW1(0,NY,NX),VOLWRX(NY,NX)
4422  FORMAT(A8,4I4,60E14.6)
C     ENDIF
C
C     WATER GAS EXCHANGE COEFFICIENTS IN SURFACE LITTER
C
C     VOLA1,VOLI1,VOLW1,VOLPM=total,ice-,water-,air-filled porosity 
C        (m3)
C     THETWA,THETWT=litter water concentration relative to air-filled
C        porosity,water holding capacity (m3 m-3)
C     DFGS=rate constant for air-water gas exchange in ‘nitro.f’ and
C        ‘trnsfr.f’ (t-1)
C     Z1R,Z2RW,Z2RD,Z3RX=parameters for litter air-water gas transfers 
C     VOLWRX=litter water retention capacity (m3)
C     XNPGX =time step for gas transfer from ‘wthr.f’(h t-1)
C     TORT=tortuosity for aqueous diffusivity
C
      VOLAT0=VOLA1(0,NY,NX)-VOLI1(0,NY,NX)
      IF(VOLAT0.GT.ZEROS2(NY,NX)
     2.AND.VOLPM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWA=AMAX1(0.0,AMIN1(1.0
     2,VOLW1(0,NY,NX)/VOLAT0))
      IF(THETWA.GT.Z3R)THEN
      DFGS(M,0,NY,NX)=AMAX1(0.0
     2,1.0/((Z1R**-1)*EXP(Z2RW*(THETWA-Z3R)))*XNPT)
      ELSE
      DFGS(M,0,NY,NX)=AMIN1(1.0
     2,1.0/((Z1R**-1)*EXP(Z2RD*(THETWA-Z3R)))*XNPT)
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
C        (mm h-1)
C     ENGYD,ENGYB=energy impact of direct,indirect 
C        precipitation+irrigation at soil surface (J t-1)
C     HV=free water depth on soil surface (m)
C     VOLWG=ground surface water retention capacity (m3)
C     XVOLW,XVOLWM=surface water in excess of litter water retention
C        capacity (m3) 
C     ZT=canopy height (m)
C     ENGYPM=total rainfall energy impact for use in ‘erosion.f’ (J t-1)
C     ENGYP=cumulative rainfall energy impact on soil surface (J)
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
C     CSILT,CCLAY=soil surface silt,clay concentration (Mg Mg-1)
C     XNPHX=time step for water fluxes from ‘wthr.f’ (h t-1)
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
     2*EXP(-2.0*HV)*BARE(NY,NX)*XNPHX
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
      DO 9875 L=1,JS
      TFLWS(L,NY,NX)=0.0
      TFLWW(L,NY,NX)=0.0
      TFLWV(L,NY,NX)=0.0
      TFLWI(L,NY,NX)=0.0
      THFLWW(L,NY,NX)=0.0  
      TK02(L)=TK0(L,NY,NX) 
9875  CONTINUE
C
C     SURFACE FLUX ACCUMULATORS
C
      DO 9885 L=NUM(NY,NX),NL(NY,NX)
      TWFLVL(L,NY,NX)=0.0
      THFLVL(L,NY,NX)=0.0
      TWFLFL(L,NY,NX)=0.0
      TWFLFH(L,NY,NX)=0.0
      THFLFL(L,NY,NX)=0.0
      TFLWL(L,NY,NX)=0.0
      TFLVL(L,NY,NX)=0.0
      TFLWLX(L,NY,NX)=0.0
      TFLWHL(L,NY,NX)=0.0
      THFLWL(L,NY,NX)=0.0
      VOLW2(L,NY,NX)=VOLW1(L,NY,NX)
      VOLV2(L,NY,NX)=VOLV1(L,NY,NX)
      VOLI2(L,NY,NX)=VOLI1(L,NY,NX)
      VOLP2(L,NY,NX)=AMAX1(0.0,VOLA1(L,NY,NX)
     2-VOLW2(L,NY,NX)-VOLI2(L,NY,NX))
C
C     GAS EXCHANGE COEFFICIENTS SOIL LAYERS
C
C     VOLA1,VOLI1,VOLW1=total,ice-,water-filled micropore volume (m3)
C     VOLAH1,VOLIH1,VOLWH1=total,ice-,water-filled macropore volume
C        (m3) 
C     VOLPM=air-filled porosity (m3)
C     THETWA=soil water concentration (m3 m-3)
C     DFGS=rate constant for air-water gas exchange in ‘trnsfr.f’ (t-1)
C     Z1S,Z2SW,Z2SD,Z3SX=parameters for soil air-water gas transfers 
C     TORT,TORTH=tortuosity for aqueous diffusivity in micropores,
C        macropores
C     XNPGX=time step for gas transfer from ‘wthr.f’(h t-1)
C
      VOLWT=VOLW1(L,NY,NX)+VOLWH1(L,NY,NX)
      VOLAT=VOLA1(L,NY,NX)+VOLAH1(L,NY,NX)
     2-VOLI1(L,NY,NX)-VOLIH1(L,NY,NX)
      IF(VOLAT.GT.ZEROS2(NY,NX)
     2.AND.VOLPM(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWA=AMAX1(0.0,AMIN1(1.0,VOLWT/VOLAT))
      Z3S=AMAX1(Z3SX,FC(L,NY,NX)/POROS(L,NY,NX))
      IF(THETWA.GT.Z3S)THEN 
      DFGS(M,L,NY,NX)=AMAX1(0.0
     2,1.0/((Z1S**-1)*EXP(Z2SW*(THETWA-Z3S)))*XNPT) 
      ELSE
      DFGS(M,L,NY,NX)=AMIN1(1.0
     2,1.0/((Z1S**-1)*EXP(Z2SD*(THETWA-Z3S)))*XNPT) 
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
      IF(BKDS(L,NY,NX).GT.ZERO
     2.AND.VOLY(L,NY,NX).GT.ZEROS(NY,NX))THEN
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
C     WRITE(*,3133)'TORT',I,J,NFZ,M,NX,NY,L
C    7,THETWT,VOLWM(M,L,NY,NX),VOLY(L,NY,NX)
C    8,TORT(M,L,NY,NX),FHOL(L,NY,NX)
3133  FORMAT(A8,7I4,40E12.4)
9885  CONTINUE
C     IF(NX.EQ.4.AND.NY.EQ.5)THEN
C     WRITE(*,3132)'FLQR1',I,J,NFZ,M,NX,NY
C    2,FLY1(NY,NX),FLQ1(NY,NX)
C    2,VHCPWM(M,1,NY,NX),VHCPWX(NY,NX)
C    2,FLH1(NY,NX),FLYM,FLQM,FLHM,FLQR1
C    3,FMAC(NUM(NY,NX),NY,NX),FGRD(NUM(NY,NX),NY,NX)
C    5,VOLAH1(NUM(NY,NX),NY,NX),FVOLAH,CCLAY(NUM(NY,NX),NY,NX)
C    4,VOLW1(NUM(NY,NX),NY,NX),VOLX(NUM(NY,NX),NY,NX),WP(L,NY,NX)
C    2,VOLT(NUM(NY,NX),NY,NX),VOLAH1(NUM(NY,NX),NY,NX)
C    5,VOLWRX(NY,NX),VOLW1(0,NY,NX),VOLI1(0,NY,NX)
C    6,PSISM1(0,NY,NX),PSISM1(NUM(NY,NX),NY,NX)
3132  FORMAT(A8,6I4,40E12.4)
C     ENDIF
C
C     ENERGY EXCHANGE VARIABLES AT SNOW SURFACE IF PRESENT
C
      RFLXW=0.0
      EFLXW=0.0
      VFLXW=0.0
      SFLXW=0.0
      HFLXW=0.0
      FLWLW=0.0 
      FLWLV=0.0 
      FLWLXW=0.0 
      FLWHLW=0.0 
      HFLWLW=0.0
      FLWRLW=0.0 
      FLVRLW=0.0 
      HFLWRLW=0.0
      FLWVLW=0.0
C
C     FLUX VARIABLES USED FOR TIME STEP IN SNOWPACK
C
C     VOLS02,VOLW02,VOLI02,VOLV02=snow,water,ice,vapor contents (m3)
C     VOLP0M=snowpack air-filled content (m3)
C     VHCPWM2=snowpack volumetric heat capacity (MJ K-1)
C     TK0M=snowpack temperature (K)
C
      DO 9765 L=1,JS
      VOLS02(L,NY,NX)=VOLS0(L,NY,NX) 
      VOLW02(L,NY,NX)=VOLW0(L,NY,NX) 
      VOLI02(L,NY,NX)=VOLI0(L,NY,NX) 
      VOLV02(L,NY,NX)=VOLV0(L,NY,NX) 
      VOLP02(L,NY,NX)=AMAX1(0.0,VOLS1(L,NY,NX)-VOLS02(L,NY,NX) 
     2-VOLI02(L,NY,NX)-VOLW02(L,NY,NX))
      HFLV0(L,NY,NX)=0.0
      WFLVW(L,NY,NX)=0.0
      WFLVS(L,NY,NX)=0.0
      HFLF0(L,NY,NX)=0.0 
      WFLFS(L,NY,NX)=0.0
      WFLFI(L,NY,NX)=0.0
      FLW0S(L,NY,NX)=0.0 
      FLW0W(L,NY,NX)=0.0 
      FLW0V(L,NY,NX)=0.0 
      FLW0I(L,NY,NX)=0.0 
      HFLW0W(L,NY,NX)=0.0
      FLQWM(M,L,NY,NX)=0.0
      VHCPWM2(L,NY,NX)=2.095*VOLS02(L,NY,NX)
     2+4.19*(VOLW02(L,NY,NX)+VOLV02(L,NY,NX))
     3+1.9274*VOLI02(L,NY,NX)
9765  CONTINUE
      VOLW2(0,NY,NX)=VOLW1(0,NY,NX)*FSNW(NY,NX)
      VOLV2(0,NY,NX)=VOLV1(0,NY,NX)*FSNW(NY,NX)
      VOLI2(0,NY,NX)=VOLI1(0,NY,NX)*FSNW(NY,NX)
      VOLA2=VOLA1(0,NY,NX)*FSNW(NY,NX)
      VOLP2(0,NY,NX)=AMAX1(0.0,VOLA2-VOLW2(0,NY,NX)
     2-VOLI2(0,NY,NX))
C
C     HEAT AND VAPOR FLUXES BETWEEN SNOWPACK AND GROUND AIR
C
C     VHCPWM2,VHCPR2,VHCPG2=volumetric heat capacity of
C        snowpack,surface litter,soil surface (MJ K-1)
C     VHCPWX=minimum heat capacity for solving snowpack water and
C        heat fluxes (MJ K-1)
C     TK02,TKR2,TKS2=snowpack, surface litter,soil surface temperatures
C        (K)
C     ALBW=snowpack albedo
C     VOLS02,VOLI02,VOLW02=snow,ice,water volumes (m3)
C     RFLX0=net radiation (MJ t-1)
C     RADXW=shortwave radiation at snowpack surface (MJ t-1)
C     THRYW=longwave radiation incident at snowpack surface (MJ t-1)
C     THRMXW=longwave radiation emitted by snowpack surface (MJ t-1)
C     THRMCW,THRMDW=net LW radiation exchange between 
C        snowpack and canopy,standing dead surfaces (MJ t-1)
C     RFLXW2=net radiation at snowpack surface (MJ t-1)
C
      IF(VHCPWM2(1,NY,NX).GT.VHCPWX(NY,NX))THEN
      TKR2=TK1(0,NY,NX)
      TKS2=TK1(NUM(NY,NX),NY,NX)
      VHCPR2=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))*FSNW(NY,NX)
     2+4.19*(VOLW2(0,NY,NX)+VOLV2(0,NY,NX))
     2+1.9274*VOLI2(0,NY,NX)
      VHCPG2=VHCP1(NUM(NY,NX),NY,NX)
      DO 3000 MM=1,NPS
      ALBW=(0.90*VOLS02(1,NY,NX)+0.30*VOLI02(1,NY,NX)
     2+0.06*VOLW02(1,NY,NX))
     2/(VOLS02(1,NY,NX)+VOLI02(1,NY,NX)+VOLW02(1,NY,NX))
      RFLX0=(1.0-ALBW)*RADXW(NY,NX)+THRYW(NY,NX)
      THRMXW=EMMGW(NY,NX)*TK02(1)**4
      RFLXW2=RFLX0-THRMXW
      DO 905 NZ=1,NP(NY,NX)
      THRMCW=EMMCW(NY,NX)*(TKC(NZ,NY,NX)**4-TK02(1)**4)
     2*FRADP(NZ,NY,NX)
      THRMDW=EMMCW(NY,NX)*(TKD(NZ,NY,NX)**4-TK02(1)**4)
     2*FRADQ(NZ,NY,NX)
      RFLXW2=RFLXW2+THRMCW+THRMDW
C     IF(IYRC.EQ.1983.AND.I.GT.300)THEN
C     WRITE(*,7760)'RFLXW2',I,J,NFZ,M,MM,NX,NY,NZ
C    2,RFLXW2,RFLX0,THRMXW,THRMCW,THRMDW,RADXW(NY,NX),THRYW(NY,NX)
C    3,ALBW,VOLS02(1,NY,NX),VOLI02(1,NY,NX),VOLW02(1,NY,NX)
C    4,VHCPWM2(1,NY,NX),VHCPWX(NY,NX)
C    5,TKC(NZ,NY,NX),TKD(NZ,NY,NX),TK02(1)
7760  FORMAT(A8,8I4,40E14.6)
C     ENDIF 
905   CONTINUE 
C
C     EVAPORATION-CONDENSATION IN SNOWPACK SURFACE
C
C     VPSV=saturated vapor concentration at snowpack surface (m3 m-3)
C     VOLV02,VOLW02,VOLS02,VOLP02=snowpack surface vapor,water,snow,air
C        content (m3) 
C     WFLVW2,WFLVS2=condensation(+ve) or evaporation(-ve)
C        in snowpack surface water,snow (m3 t-1)
C     HFLVX=latent heat of condensation,evaporation from
C        snowpack water+snow (MJ t-1)
C     XWFLVW,XWFLVS=aggregated condensation(+ve), evaporation(-ve) from
C        water,snow used in ‘redist.f’(m3 t-1)
C     XHFLV0=aggregated latent heat of condensation,evaporation from
C        water+snow used in ‘redist.f’ (MJ t-1)
C     XNPSX=time step for snowpack fluxes from ‘wthr.f’ (h t-1)
C
      VPSV=2.173E-03/TK02(1)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK02(1)))
      IF(VOLP02(1,NY,NX).GT.ZEROS(NY,NX))THEN
      WFLVT=VOLV02(1,NY,NX)-VPSV*VOLP02(1,NY,NX)
      WFLVW2=AMAX1(WFLVT,-AMAX1(0.0,VOLW02(1,NY,NX))*XNPSX)
      WFLVX=AMIN1(0.0,WFLVT-WFLVW2)
      WFLVS2=AMAX1(WFLVX,-AMAX1(0.0,VOLS02(1,NY,NX))*XNPSX)
      HFLVX=VAP*WFLVW2+VAPS*WFLVS2
      ELSE
      WFLVT=0.0
      WFLVX=0.0
      WFLVW2=0.0
      WFLVS2=0.0
      HFLVX=0.0
      ENDIF
      WFLVW(1,NY,NX)=WFLVW(1,NY,NX)+WFLVW2
      WFLVS(1,NY,NX)=WFLVS(1,NY,NX)+WFLVS2 
      HFLV0(1,NY,NX)=HFLV0(1,NY,NX)+HFLVX      
      XWFLVW(1,NY,NX)=XWFLVW(1,NY,NX)+WFLVW2      
      XWFLVS(1,NY,NX)=XWFLVS(1,NY,NX)+WFLVS2      
      XHFLV0(1,NY,NX)=XHFLV0(1,NY,NX)+HFLVX
      VOLV02(1,NY,NX)=VOLV02(1,NY,NX)-WFLVS2-WFLVW2
      VOLS02(1,NY,NX)=VOLS02(1,NY,NX)+WFLVS2
      VOLW02(1,NY,NX)=VOLW02(1,NY,NX)+WFLVW2 
      VOLI02(1,NY,NX)=VOLI02(1,NY,NX)
      VOLP02(1,NY,NX)=AMAX1(0.0,VOLS1(1,NY,NX)-VOLS02(1,NY,NX) 
     2-VOLI02(1,NY,NX)-VOLW02(1,NY,NX))
C     WRITE(*,7759)'WFLVWS',I,J,NFZ,M,MM,NX,NY
C    2,WFLVT,VOLV02(1,NY,NX),VPSV*VOLP02(1,NY,NX),VPSV,VOLP02(1,NY,NX)
C    3,WFLVW2,WFLVT,VOLW02(1,NY,NX),TK02(1)  
C    4,WFLVX,VOLS02(1,NY,NX),HFLVX,VOLV0(1,NY,NX)
C 
C     VAPOR FLUX AT SNOWPACK SURFACE
C
C     VPSV=saturated vapor concentration at snowpack surface (m3 m-3)
C     TK02=snowpack surface temperature (K)
C     VOLV02,VOLW02,VOLS02,VOLP02=snowpack surface vapor,water,snow,air
C        content (m3) 
C     VP0,VPQG=vapor concentration at snowpack surface, in air above
C        ground surface (m3 m-3)
C     EVAP02X=unlimited vapor flux at snowpack surface (m3 t-1)
C     PAREWM=conductance for latent heat flux (m t-1) 
C     EVAP02V,EVAP02W,EVAP02S=vapor flux from vapor,water,snow at
C        snowpack surface (m3 t-1)   
C     VFLXW2=convective heat of vapor flux (MJ t-1)    
C     FLW0S,FLQ0I,FLQ0W=snow,ice,water input to snowpack (m3 h-1)
C     XNPS=time step from ‘wthr.f’(h t-1)
C
      IF(VOLP02(1,NY,NX).GT.ZEROS(NY,NX))THEN
      VP0=AMAX1(0.0,VOLV02(1,NY,NX)/VOLP02(1,NY,NX))
      ELSE
      VP0=VPSV
      ENDIF
      EVAP02X=PAREWM*(VPQG(M,NY,NX)-VP0)
      EVAP02V=AMAX1(EVAP02X
C    2,-AMAX1(0.0,VOLV02(1,NY,NX)*XNPSX))
     2,-AMAX1(0.0,VOLV02(1,NY,NX)))
      EVAP02W=AMAX1(EVAP02X-EVAP02V
     2,-AMAX1(0.0,VOLW02(1,NY,NX)*XNPSX))
      EVAP02S=AMAX1(EVAP02X-EVAP02V-EVAP02W
     2,-AMAX1(0.0,VOLS02(1,NY,NX)*XNPSX))
      EVAP02=EVAP02V+EVAP02W+EVAP02S 
      IF(EVAP02.LT.0.0)THEN
      VFLXW2=((EVAP02V+EVAP02W)*4.19+EVAP02S*2.095)*TK02(1)
      ELSE
      VFLXW2=((EVAP02V+EVAP02W)*4.19+EVAP02S*2.095)*TKQG(M,NY,NX)
      ENDIF
      FLQ0S2=FLQ0S(NY,NX)*XNPS 
      FLQ0W2=FLQ0W(NY,NX)*XNPS 
      FLQ0I2=FLQ0I(NY,NX)*XNPS
      HWFLQ02=HWFLQ0(NY,NX)*XNPS 
      VOLV02(1,NY,NX)=VOLV02(1,NY,NX)+EVAP02V
      VOLS02(1,NY,NX)=VOLS02(1,NY,NX)+FLQ0S2+EVAP02S 
      VOLW02(1,NY,NX)=VOLW02(1,NY,NX)+FLQ0W2+EVAP02W 
      VOLI02(1,NY,NX)=VOLI02(1,NY,NX)+FLQ0I2
C
C     SOLVE FOR SNOWPACK SURFACE TEMPERATURE AT WHICH ENERGY
C     BALANCE OCCURS, SOLVE AND ACCUMULATE LATENT, SENSIBLE 
C     STORAGE HEAT FLUXES AND EVAPORATION
C
C     SFLXW2,EFLXW2,RFLXW2=sensible,latent heat fluxes, 
C        net radiation at snowpack surface (MJ t-1)
C     HFLX02,HFLXW2=storage,total heat flux (MJ t-1)
C     VAP,VAPS=latent heat of evaporation,sublimation from ‘starts.f’
C        (MJ m-3)
C     EVAP02W,EVAP02S=vapor flux from water,snow at
C        snowpack surface (m3 t-1)   
C     WFLVW2,WFLVS2=condensation(+ve) or evaporation(-ve)
C        in snowpack surface water,snow (m3 t-1)
C     VFLXW2=convective heat flux from EFLXW2 (MJ t-1)
C     PARSWM=conductance for sensible heat flux (MJ K-1 t-1)
C     TKQG,TK02=air temperature at ground surface,
C        snowpack surface temperature (K)
C     FLQ0S2,FLQ0W2,FLQ0I2=snow,water,ice input to snowpack (m3 t-1)
C     HWFLQ02=convective heat from snow,water,ice input to snowpack 
C        (MJ t-1)
C     VOLV02,VOLW02,VOLS02,VOLP02=snowpack surface vapor,water,snow,air
C        content (m3) 
C
      EFLXW2=(WFLVW2+EVAP02W)*VAP+(WFLVS2+EVAP02S)*VAPS
      SFLXW2=PARSWM*(TKQG(M,NY,NX)-TK02(1))
      HFLX02=RFLXW2+EFLXW2+SFLXW2 
      HFLXW2=HFLX02+VFLXW2
      RFLXW=RFLXW+RFLXW2
      EFLXW=EFLXW+EFLXW2
      VFLXW=VFLXW+VFLXW2
      SFLXW=SFLXW+SFLXW2
      HFLXW=HFLXW+HFLXW2
      EVAP0(NY,NX)=EVAP0(NY,NX)+EVAP02
      FLW0V(1,NY,NX)=EVAP02V
      FLW0S(1,NY,NX)=FLQ0S2+EVAP02S 
      FLW0W(1,NY,NX)=FLQ0W2+EVAP02W 
      FLW0I(1,NY,NX)=FLQ0I2
      HFLW0W(1,NY,NX)=HWFLQ02+HFLXW2
      FLQWM(M,1,NY,NX)=FLQWM(M,1,NY,NX)+FLQ0S2+FLQ0I2+FLQ0W2
C     IF(I.GT.270)THEN
C     WRITE(*,7759)'EFLXW',I,J,NFZ,M,MM,NX,NY
C    2,EVAP02X,EVAP02V,EVAP02W,EVAP02S,EVAP02,WFLVT,WFLVW2,WFLVS2 
C    3,FSNW(NY,NX),DPTHS0(NY,NX),DPTHSX
C    4,RFLXW2,EFLXW2,SFLXW2,VFLXW2,HFLXW2,HWFLQ02 
C    5,VP0,VPQG(M,NY,NX),VPQGX(NY,NX),VPSV,PAREWM,VHCPWM2(1,NY,NX) 
C    5,VOLV02(1,NY,NX),VOLW02(1,NY,NX),VOLP02(1,NY,NX)
C    7,VOLI02(1,NY,NX),TK02(1),TKR2,TKS2,TKQG(M,NY,NX),TKQGX(NY,NX)
C    8,RAH,PARSWM,PAREWM,RARH,RAH
C    9,TKGS(M,NY,NX),TK0(1,NY,NX),TKW(1,NY,NX) 
C    9,(VOLSSL(L,NY,NX),L=1,5),(VOLWSL(L,NY,NX),L=1,5)
C    9,(VOLISL(L,NY,NX),L=1,5)
7759  FORMAT(A8,7I4,120E14.6) 
C     ENDIF 
C
C     PHYSICAL AND HYDRAULIC PROPERTIES OF SNOWPACK LAYERS INCLUDING
C     AIR AND WATER-FILLED POROSITY, WATER POTENTIAL OF UNDERLYING
C     SOIL SURFACE USED IN FLUX CALCULATIONS
C
C     VHCPWMM,VHCPWX=current, minimum snowpack layer heat capacities 
C        (MJ K-1)
C     VOLS02,VOLI02,VOLW02,VOLS1=snow,ice,water,total snowpack 
C        layer volume (m3)
C     DLYRS0=snowpack layer depth (m)
C     DENSS,DENSI,DENS0=snow,ice,minimum snow density (Mg m-3)
C     AREA=area of grid cell (m2)
C     VOLP0M=snowpack layer air volume (m3)
C     THETP1=snowpack layer air concentration (m3 m-3)   
C     CNV1=snowpack layer vapor conductivity (m2 h-1)
C     WGSGW=snowpack layer vapor diffusivity from ‘hour1.f’ (m2 h-1)
C     VP1=snowpack layer vapor concentration (m3 m-3)
C     TK0M=snowpack layer temperature (K)
C     DENSW1=snowpack layer density (Mg m-3)
C     
      ICHKL=0
      DO 9880 L=1,JS
      IF(VHCPWM2(L,NY,NX).GT.VHCPWX(NY,NX))THEN
      VOLS1(L,NY,NX)=VOLS02(L,NY,NX)/DENSS(L,NY,NX)
     2+VOLW02(L,NY,NX)+VOLI02(L,NY,NX)
      DLYRS0(L,NY,NX)=VOLS1(L,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
      THETP1=AMAX1(THETPI,VOLP02(L,NY,NX)/VOLS1(L,NY,NX))
      CNV1=THETP1**2.0*WGSGW(L,NY,NX)
      IF(VOLP02(L,NY,NX).GT.ZEROS(NY,NX))THEN
      VP1=AMAX1(0.0,VOLV02(L,NY,NX)/VOLP02(L,NY,NX))
      ELSE
      VP1=0.0
      ENDIF
      IF(VOLS1(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DENSW1=AMIN1(0.6,(VOLS02(L,NY,NX)+VOLW02(L,NY,NX)
     2+VOLI02(L,NY,NX)*DENSI)/VOLS1(L,NY,NX))
      ELSE
      DENSW1=DENS0(NY,NX)
      ENDIF
C
C     SNOW THERMAL CONDUCTIVITY FROM J GLACIOL 43:26-41
C
C     TCND1W=snow thermal conductivity (m MJ h-1 K-1)
C     DENSW1=snowpack layer density (Mg m-3)
C
      TCND1W=0.0036*10**(2.650*DENSW1-1.652)
C
C     DISCHARGE OF MELTWATER AND ITS HEAT FROM SNOWPACK LAYER
C     TO LOWER SNOWPACK LAYER 
C
C     FLWQX=porosity-unlimited snow water flux (m3 t-1)
C     VOLS02,VOLW02 =snow,water volume of snowpack layer (m3)
C     XNPSX=time step for snowpack fluxes from ‘wthr.f’(h t-1) 
C
      FLWQX=AMAX1(0.0,AMAX1(0.0,VOLW02(L,NY,NX))
     2-0.05*AMAX1(0.0,VOLS02(L,NY,NX)))*XNPSX
C
C     WATER AND HEAT FLUXES IN SNOWPACK
C
C     VHCPWMM,VHCPWX=current, minimum snowpack layer heat capacities
C        (MK K-1)
C     VOLS02,VOLI02,VOLW02,VOLS1=snow,ice,water,total snowpack 
C        layer volume (m3)
C     DENSS=snow density (Mg m-3)
C     DLYRS0=snow layer thickness (m)
C     VOLP0M=snowpack layer air volume (m3)
C     THETP2=snowpack layer air concentration (m3 m-3)   
C     FLWQM=porosity-limited snowpack water flux (m3 t-1)
C     HFLWQM=convective heat flux from water flux (MJ t-1)
C     TK02=snowpack temperature (K)
C
      L2=MIN(JS,L+1)
      IF(L.LT.JS.AND.VHCPWM2(L2,NY,NX).GT.VHCPWX(NY,NX))THEN
      VOLS1(L2,NY,NX)=VOLS02(L2,NY,NX)/DENSS(L2,NY,NX)
     2+VOLW02(L2,NY,NX)+VOLI02(L2,NY,NX)
      DLYRS0(L2,NY,NX)=VOLS1(L2,NY,NX)/AREA(3,NUM(NY,NX),NY,NX)
      THETP2=AMAX1(THETPI,VOLP02(L2,NY,NX)/VOLS1(L2,NY,NX))
      FLWQM=AMIN1(THETP2,FLWQX)
      HFLWQM=4.19*TK02(L)*FLWQM
C
C     VAPOR FLUX IN SNOWPACK 
C
C     VOLP0M=air-filled volumes of snowpack layers (m3)
C     L2=destination layer
C     CNV1,CNV2=vapor conductivities of source, destination layers 
C        (m2 h-1)
C     VP1,VP2=vapor concentrations of source, destination layers 
C        (m3 m-3)
C     TK0M=soil layer temperature (K)
C     AVCNVW=snow vapor conductance (m h-1)
C     DLYRS0=snow layer thickness (m)
C     FLVC=vapor-unlimited vapor flux (m3 t-1)
C     AREA=area of grid cell (m2)
C     FLVSS,HFLVSS=vapor flux and its convective heat flux (m3,MJ t-1)
C     XNPYX=time step for snowpack fluxes from ‘wthr.f’(h t-1)
C
      IF(VOLP02(L,NY,NX).GT.ZEROS2(NY,NX)
     2.AND.VOLP02(L2,NY,NX).GT.ZEROS2(NY,NX))THEN
      CNV2=THETP2**2.0*WGSGW(L2,NY,NX)
      VP2=AMAX1(0.0,VOLV02(L2,NY,NX)/VOLP02(L2,NY,NX))
      AVCNVW=2.0*CNV1*CNV2/(CNV1*DLYRS0(L2,NY,NX)
     2+CNV2*DLYRS0(L,NY,NX)) 
      FLVC=AVCNVW*(VP1-VP2)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)
     2*XNPYX 
      IF(FLVC.GE.0.0)THEN
      FLVSS=AMAX1(0.0,AMIN1(FLVC,VOLV02(L,NY,NX)*XNPXX))
      HFLVSS=4.19*TK02(L)*FLVSS
      ELSE
      FLVSS=AMIN1(0.0,AMAX1(FLVC,-VOLV02(L2,NY,NX)*XNPXX))
      HFLVSS=4.19*TK02(L2)*FLVSS
      ENDIF
      ELSE
      FLVSS=0.0
      HFLVSS=0.0
      ENDIF
C
C     HEAT FLUX IN SNOWPACK
C
C     VOLS02,VOLI02,VOLW02,VOLS1=snow,ice,water,total snowpack 
C        layer volume (m3)
C     DENSW2=density in destination layer (Mg m-3)
C     TCNDW2=thermal conductivity in destination layer (m MJ h-1 K-1)
C     ATCNDW=thermal conductance (MJ h-1 K-1)
C     DLYRS0=layer thickness (m)
C     TKY=equilibrium temperature (K)
C     HFLWX,HFLWC=heat-limited,heat-unlimited heat fluxes 
C        (MJ t-1)
C     AREA=area of grid cell (m2)
C     VHCPWMM=snowpack heat capacity (MJ K-1)
C     TK0M=snowpack temperature (K)
C     FSNW=snow cover fraction
C     HFLWSS=snowpack heat flux
C     FLW0S,FLQ0I,FLQ0W,FLQ0V=snow,ice,water,vapor fluxes (m3 t-1)
C     HFLW0W=convective heat flux (MJ t-1)
C     TFLWS,TFLWI,TFLWW,TFLWV=net snow,ice,water,vapor fluxes (m3 t-1)
C     THFLWW=net convective heat flux (MJ t-1)
C     XNPYX=time step for snowpack flux from ‘wthr.f’(h t-1)
C
      IF(VOLS1(L2,NY,NX).GT.ZEROS2(NY,NX))THEN
      DENSW2=AMIN1(0.6,(VOLS02(L2,NY,NX)+VOLW02(L2,NY,NX)
     2+VOLI02(L2,NY,NX)*DENSI)/VOLS1(L2,NY,NX))
      ELSE
      DENSW2=DENS0(NY,NX)
      ENDIF
      TCND2W=0.0036*10**(2.650*DENSW2-1.652)
      ATCNDW=2.0*TCND1W*TCND2W/(TCND1W*DLYRS0(L2,NY,NX)
     2+TCND2W*DLYRS0(L,NY,NX)) 
      HFLWC=ATCNDW*(TK02(L)-TK02(L2)) 
     2*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX)*XNPYX 
      TKY=(TK02(L)*VHCPWM2(L,NY,NX)+TK02(L2) 
     2*VHCPWM2(L2,NY,NX))/(VHCPWM2(L,NY,NX)+VHCPWM2(L2,NY,NX))
      HFLWX=(TK02(L)-TKY)*VHCPWM2(L,NY,NX) 
      IF(HFLWC.GE.0.0)THEN
      HFLWSS=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWSS=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
      HFLW0T=HFLWQM+HFLVSS+HFLWSS 
      FLW0S(L2,NY,NX)=0.0 
      FLW0W(L2,NY,NX)=FLWQM
      FLW0V(L2,NY,NX)=FLVSS
      FLW0I(L2,NY,NX)=0.0 
      HFLW0W(L2,NY,NX)=HFLW0T 
      FLQWM(M,L2,NY,NX)=FLQWM(M,L2,NY,NX)+FLWQM
      TFLWSX=FLW0S(L,NY,NX)-FLW0S(L2,NY,NX)
      TFLWWX=FLW0W(L,NY,NX)-FLW0W(L2,NY,NX)
      TFLWVX=FLW0V(L,NY,NX)-FLW0V(L2,NY,NX)
      TFLWIX=FLW0I(L,NY,NX)-FLW0I(L2,NY,NX)
      THFLWWX=HFLW0W(L,NY,NX)-HFLW0W(L2,NY,NX)
      TFLWS(L,NY,NX)=TFLWS(L,NY,NX)+TFLWSX 
      TFLWW(L,NY,NX)=TFLWW(L,NY,NX)+TFLWWX 
      TFLWV(L,NY,NX)=TFLWV(L,NY,NX)+TFLWVX 
      TFLWI(L,NY,NX)=TFLWI(L,NY,NX)+TFLWIX
      THFLWW(L,NY,NX)=THFLWW(L,NY,NX)+THFLWWX
C
C     UPDATE INTERMEDIATE STATE VARIABLES FOR SNOWPACK WATER,VAPOR,AIR
C
C     VOLW02,VOLV02,	voli02,VOLP02,VOLS1= snowpack
C        water,vapor,ice,air,layer volume (m3)
C     FLQ0W,FLQ0V=water,vapor fluxes (m3 t-1)
C
      VOLW02(L,NY,NX)=VOLW02(L,NY,NX)-FLW0W(L2,NY,NX)
      VOLV02(L,NY,NX)=VOLV02(L,NY,NX)-FLW0V(L2,NY,NX)
      VOLW02(L2,NY,NX)=VOLW02(L2,NY,NX)+FLW0W(L2,NY,NX)
      VOLV02(L2,NY,NX)=VOLV02(L2,NY,NX)+FLW0V(L2,NY,NX)
      VOLP02(L,NY,NX)=AMAX1(0.0,VOLS1(L,NY,NX)-VOLS02(L,NY,NX) 
     2-VOLI02(L,NY,NX)-VOLW02(L,NY,NX))
      VOLP02(L2,NY,NX)=AMAX1(0.0,VOLS1(L2,NY,NX)-VOLS02(L2,NY,NX) 
     2-VOLI02(L2,NY,NX)-VOLW02(L2,NY,NX))
C     IF(NX.EQ.3.AND.NY.EQ.3.AND.L.EQ.1)THEN
C     WRITE(*,7757)'FLW0',I,J,M,MM,L2,FLW0W(L2,NY,NX),FLW0V(L2,NY,NX)
C    2,HFLW0W(L2,NY,NX),HFLW0T,HFLWQM,HFLVSS,HFLWSS,VP1,VP2
C    3,TCND1W,TCND2W,DENSW1,DENSW2,HFLXW2,THETP2,FLWQX
C    2,VOLS02(L,NY,NX),VOLW02(L,NY,NX),TK02(L2),TK02(L)
7757  FORMAT(A8,5I4,30E14.6)
C     ENDIF
C
C     DISCHARGE OF MELTWATER AND ITS HEAT FROM LOWEST SNOWPACK LAYER
C     TO RESIDUE, SURFACE SOIL MICROPORES AND MACROPORES
C
C     FLWQX,FLWQR=porosity-unlimited water flux to soil,litter 
C        (m3 t-1)
C     FLWQGX,FLWQGS,FLWQGH=water flux to soil surface,
C        micropores,macropores (m3 t-1)
C     VOLP1,VOLPH1=air volumes of soil micropores,macropores (m3)
C     FMAC,FGRD=macropore,micropore volume fractions 
C     HFLWQG,HFLWQR=convective heat fluxes to soil,litter (MJ t-1)
C     THETWR,THETW1=litter, soil water concentration (m3 m-3)
C     VOLR,VOLWRX=litter volume,water retention capacity (m3)
C     PSISM1(0,PSISM1(NUM=litter,soil surface water potentials (MPa)
C     POROS=soil porosity (m3 m-3)
C     FC,WP=field capacity,wilting point from soil file (m3 m-3)
C     FCL,WPL=log(FC),log(WP)
C     FCI,WPI=FC,WP of ice
C     PSL,FCL,WPL=log POROS0,FC,WP
C     FCD,PSD=FCL-WPL,PSL-FCL
C     SRP=parameter for deviation from linear log-log water retention 
C     THETIX=ice concentration (m3 m-3)
C     BKVL=bulk density x volume of soil layer (Mg)
C
      ELSE
      IF(ICHKL.EQ.0)THEN
      FLWQGX=FLWQX*BARE(NY,NX)
      FLWQGS=AMIN1(VOLP1(NUM(NY,NX),NY,NX)*XNPSX 
     2,FLWQGX*FGRD(NUM(NY,NX),NY,NX))
      FLWQGH=AMIN1(VOLPH1(NUM(NY,NX),NY,NX)*XNPSX 
     2,FLWQGX*FMAC(NUM(NY,NX),NY,NX))
      FLWQG=FLWQGS+FLWQGH
      HFLWQG=4.19*TK02(L)*FLWQG
      FLWQR=FLWQX-FLWQG 
      HFLWQR=4.19*TK02(L)*FLWQR
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX)
     2.AND.VOLW1(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW1(0,NY,NX))/VOLR(NY,NX)
      IF(THETWR.LT.FC(0,NY,NX))THEN
      PSISM1(0,NY,NX)=AMAX1(PSISX,-EXP(PSIMX(NY,NX)
     2+((FCL(0,NY,NX)-LOG(THETWR))
     3/FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETWR.LT.POROS0(NY,NX))THEN 
      PSISM1(0,NY,NX)=-EXP(PSIMS(NY,NX)
     2+((AMAX1(0.0,(PSL(0,NY,NX)-LOG(THETWR)))
     3/PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
      ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
      ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
      THETW1=AMAX1(THETZ(NUM(NY,NX),NY,NX)
     2,AMIN1(POROS(NUM(NY,NX),NY,NX)
     2,VOLW1(NUM(NY,NX),NY,NX)/VOLY(NUM(NY,NX),NY,NX)))
      IF(BKVL(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      IF(THETW1.LT.FC(NUM(NY,NX),NY,NX))THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSISX,-EXP(PSIMX(NY,NX)
     2+((FCL(NUM(NY,NX),NY,NX)-LOG(THETW1))
     3/FCD(NUM(NY,NX),NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETW1.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN 
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX)
     2+((AMAX1(0.0,(PSL(NUM(NY,NX),NY,NX)-LOG(THETW1)))
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
      PSDX=PSL(NUM(NY,NX),NY,NX)-FCLX
      FCDX=FCLX-WPLX
      IF(THETWX(NUM(NY,NX),NY,NX).LT.FCX)THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSISX,-EXP(PSIMX(NY,NX)
     2+((FCLX-LOG(THETWX(NUM(NY,NX),NY,NX)))
     3/FCDX*PSIMD(NY,NX))))
      ELSEIF(THETWX(NUM(NY,NX),NY,NX)
     2.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN 
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX)
     2+((AMAX1(0.0,(PSL(NUM(NY,NX),NY,NX)
     2-LOG(THETWX(NUM(NY,NX),NY,NX))))
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
C
C     VAPOR FLUX BETWEEN SNOWPACK AND SOIL SURFACE
C
C     VOLP02,VOLP2=snowpack,soil air volume (m3)
C     CNV1,CNV2=vapor conductances of source, destination layers
C        (m2 h-1)
C     VP1,VP2=vapor concentrations of source, destination layers 
C        (m3 m-3)
C     POROS,POROQ=porosity (m3 m-3), tortuosity from ‘starts.f’  
C     WGSGL=vapor diffusivity (m2 h-1)
C     TK0M,TK1=snow,soil surface temperature (K) 
C     VOLV2=soil vapor content (m3)
C     AVCNVS=snow-soil vapor conductance (m h-1)
C     DLYR,DLYRS0=soil surface,snowpack layer depth (m)
C     FLVC=vapor flux unlimited by vapor (m3 t-1)
C     AREA=area of grid cell (m2)
C     FSNW,CVRD=snow,litter cover fraction 
C     FLVS1,HFLVS1=vapor flux and its convective heat flux 
C        (m3 t-1,MJ t-1)
C     XNPYX=time step for flux calculations from ‘wthr.f’ (h t-1)    
C
      IF(VOLP02(L,NY,NX).GT.ZEROS2(NY,NX)
     2.AND.VOLP2(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      CNV2=WGSGL(NUM(NY,NX),NY,NX)*POROQ
     2*THETPM(M,NUM(NY,NX),NY,NX)**2/POROS(NUM(NY,NX),NY,NX)
      VP2=AMAX1(0.0,VOLV2(NUM(NY,NX),NY,NX)/VOLP2(NUM(NY,NX),NY,NX))
      AVCNVS=2.0*CNV1*CNV2
     2/(CNV1*DLYR(3,NUM(NY,NX),NY,NX)+CNV2*DLYRS0(L,NY,NX)) 
      FLVC=AVCNVS*(VP1-VP2)*AREA(3,NUM(NY,NX),NY,NX) 
     2*FSNW(NY,NX)*BARE(NY,NX)*XNPYX 
      IF(FLVC.GE.0.0)THEN
      FLVS1=AMAX1(0.0,AMIN1(FLVC,VOLV02(L,NY,NX)*XNPXX))
      HFLVS1=4.19*TK02(L)*FLVS1
      ELSE
      FLVS1=AMIN1(0.0,AMAX1(FLVC,-VOLV2(NUM(NY,NX),NY,NX)*XNPXX))
      HFLVS1=4.19*TK1(NUM(NY,NX),NY,NX)*FLVS1
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
C        (m MJ h-1 K-1)
C     STC,DTC=mineral component of thermal conductivity from ‘hour1.f’
C        (m MJ h-1 K-1)
C     THETWX,THETIX,THETPX=soil surface water,ice,air concentrations
C        (m3 m-3)
C     BARE=soil surface fraction
C     ATCNDS=snowpack-soil thermal conductance (MJ h-1 K-1)
C     TKY=equilibrium temperature (K)
C     HFLWX,HFLWC=heat-limited,heat-unlimited heat fluxes 
C        (MJ t-1)
C     HFLWS1=snowpack-soil heat flux (MJ t-1)
C     XNPYX=time step for snowpack flux from ‘wthr.f’(h t-1) (h t-1)     
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
      HFLWC=ATCNDS*(TK02(L)-TKS2)*AREA(3,NUM(NY,NX),NY,NX) 
     2*FSNW(NY,NX)*BARE(NY,NX)*XNPYX 
      TKY=(TK02(L)*VHCPWM2(L,NY,NX)
     2+TKS2*VHCP1(NUM(NY,NX),NY,NX))
     3/(VHCPWM2(L,NY,NX)+VHCP1(NUM(NY,NX),NY,NX))
      HFLWX=(TK02(L)-TKY)*VHCPWM2(L,NY,NX) 
      IF(HFLWC.GE.0.0)THEN
      HFLWS1=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWS1=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
C     IF(I.GT.200)THEN
C     WRITE(*,1113)'HFLWS1',I,J,NFZ,M,MM,L 
C    2,HFLWS1,HFLWX,HFLWC,TK02(L),TKS2,TKY,ATCNDS,FSNW(NY,NX)
C    2,BARE(NY,NX),VHCPWM2(L,NY,NX),VHCP1(NUM(NY,NX),NY,NX)
C    3,TCND1W,TCNDS,DLYR(3,NUM(NY,NX),NY,NX),DLYRS0(L,NY,NX)
C    3,VOLW1(L,NY,NX),VOLV1(L,NY,NX),VOLWH1(L,NY,NX)
C    3,VOLI1(L,NY,NX),VOLIH1(L,NY,NX),VHCM(L,NY,NX) 
1113  FORMAT(A8,6I4,30E12.4)
C     ENDIF
C
C     HEAT FLUX AMONG SNOWPACK, SURFACE RESIDUE AND SURFACE SOIL
C
C     FLVSR=snowpack-litter vapor flux (m3 t-1)
C     HFLVSR,HFLWSR=snowpack-litter convective,conductive heat fluxes
C        (MJ t-1)
C     FLVS1=snowpack-soil vapor flux (m3 t-1)
C     HFLVS1,HFLWS1=snowpack-soil convective,conductive heat fluxes
C        (MJ t-1)
C     VHCP1,VHCPRX=current,minimum litter heat capacities (MJ K-1)
C     TK0X,TKXR,TK1X=snowpack,litter,soil temperatures (K)
C     CNVR,CNV1,CNV2=litter,snowpack,soil vapor conductivity (m2 h-1)
C     THETP*,THETWX,THETIX=litter air,water,ice concentration (m3 m-3)
C     POROS,POROQ=litter porosity (m3 m-3), tortuosity from ‘starts.f’ 
C     CVRD=litter cover fraction
C     WGSGR=litter vapor diffusivity (m2 h-1)
C     AVCNVR,AVCNVS=snowpack-litter,litter-soil vapor conductance 
C        (m h-1)
C     DLYRR,DLYRS0,DLYR=litter,snowpack,soil depths (m)
C     THETRR=dry litter concentration (m3 m-3)
C     TCNDR,TCND1W,TCNDS=litter,snowpack,soil thermal conductivity
C        (m MJ h-1 K-1)
C     ATCNDR,ATCNDS=snow-litter,litter-soil thermal conductance
C        (MJ h-1 K-1)
C
      FLVSR=0.0
      HFLVSR=0.0
      HFLWSR=0.0
      FLVR1=0.0
      HFLVR1=0.0
      HFLWR1=0.0
      FLVR2T=0.0
      HFLVR2T=0.0
      FLFR2T=0.0
      HFLFR2T=0.0
      FLVG2T=0.0
      HFLVG2T=0.0
      FLFG2T=0.0
      HFLFG2T=0.0
      IF(VHCP1(0,NY,NX).GT.VHCPRXS(NY,NX))THEN
      CNVR=WGSGR(NY,NX)*POROQ
     2*THETPM(M,0,NY,NX)**2/POROS(0,NY,NX)
      IF(CVRD(NY,NX).GT.ZERO)THEN
      IF(CNV1.GT.ZERO.AND.CNVR.GT.ZERO)THEN
      AVCNVR=2.0*CNVR*CNV1/(CNV1*DLYRR(NY,NX)+CNVR*DLYRS0(L,NY,NX)) 
      ELSE
      AVCNVR=2.0*CNV1/(DLYRR(NY,NX)+DLYRS0(L,NY,NX)) 
      ENDIF
      IF(CNVR.GT.ZERO.AND.CNV2.GT.ZERO)THEN
      AVCNVS=2.0*CNVR*CNV2
     2/(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+CNV2*DLYRR(NY,NX))
      ELSE
      AVCNVS=2.0*CNV2/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR(NY,NX)) 
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
      AVCNVR=0.0
      AVCNVS=0.0
      ATCNDR=0.0
      ATCNDS=0.0
      ENDIF
C
C     SHORTER TIME STEP FOR SURFACE RESIDUE FLUX CALCULATIONS
C
C     INITIALIZE INTERMEDIATE WATER, VAPOR, ICE CONTENTS AND 
C     TEMPERATURES OF LOWEST SNOW LAYER, SURFACE LITTER AND SOIL
C     SURFACE LAYER
C 
      VOLS02S=VOLS02(L,NY,NX)
      VOLW02S=VOLW02(L,NY,NX) 
      VOLV02S=VOLV02(L,NY,NX) 
      VOLI02S=VOLI02(L,NY,NX)
      VOLP02S=AMAX1(0.0,VOLS1(L,NY,NX)-VOLS02S-VOLI02S-VOLW02S)
      VOLW2R=VOLW2(0,NY,NX)
      VOLV2R=VOLV2(0,NY,NX)
      VOLI2R=VOLI2(0,NY,NX)
      VOLP2R=VOLP2(0,NY,NX)
      VOLW2G=VOLW2(NUM(NY,NX),NY,NX) 
      VOLV2G=VOLV2(NUM(NY,NX),NY,NX) 
      VOLI2G=VOLI2(NUM(NY,NX),NY,NX) 
      VOLP2G=AMAX1(0.0,VOLA1(NUM(NY,NX),NY,NX)-VOLW2G-VOLI2G)
      TK022=TK02(L)
      TKR22=TKR2
      TKS22=TKS2
      PSISVR=PSISM1(0,NY,NX)+PSISO(0,NY,NX)
      PSISVG=PSISM1(NUM(NY,NX),NY,NX)+PSISO(NUM(NY,NX),NY,NX)
      DO 4000 NN=1,NPRS
C
C     EVAPORATION-CONDENSATION IN SURFACE LITTER UNDER SNOW
C
C     VPRV=litter saturated vapor concentration (m3 m-3)
C     TKR22,PSISVR=litter temperature,water potential (K,MPa)
C     FLVRX,FLVR2=litter evaporation-condensation unlimited,
C        limited by vapor (m3 t-1)
C     HFLVR2=litter latent heat flux (MJ t-1)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C
      VPRV=2.173E-03/TKR22 
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKR22))
     3*EXP(18.0*PSISVR/(8.3143*TKR22))
      FLVRX=(VOLV2R-VPRV*VOLP2R)
      FLVR2=AMAX1(FLVRX,-AMAX1(0.0,VOLW2R)*XNPSRX)
      HFLVR2=VAP*FLVR2
      FLVR2T=FLVR2T+FLVR2
      HFLVR2T=HFLVR2T+HFLVR2
C     WRITE(*,7769)'FLVR2',I,J,NFZ,M,MM,NN,NX,NY
C    2,FLVR2,FLVRX,FLVR2T,VOLW2R,XNPSRX,VOLV2R,VPRV*VOLP2R
C    3,VPRV,VOLP2R,FSNW(NY,NX),TKR22,VOLV2(0,NY,NX)
C
C     FREEZE-THAW IN SURFACE LITTER UNDER SNOW FROM NET CHANGE IN 
C     LITTER SURFACE HEAT STORAGE
C
C     TFREEZ=litter freezing temperature (K)
C     333.0=latent heat of freezing (MJ m-3)
C     TKR22=litter temperature (K)
C     PSISVR=litter natric+osmotic potential (MPa)
C     VOLW2R,VOLI2R=litter water,ice volume (m3)
C     VOLT=litter volume (m3)
C     VHCPR2=litter heat capacity (MJ K-1)
C     HFLFRX,HFLFR2=litter freeze-thaw latent heat flux
C        unlimited,limited by water,ice (MJ t-1)
C     DENSI=ice density (Mg m-3)
C     FLFR2=litter freeze-thaw (m3 t-1)
C     XNPRS=time step for litter flux from ‘wthr.f’(h t-1)
C     
      TFREEZ=-9.0959E+04/(PSISVR-333.0)
      IF((TKR22.LT.TFREEZ
     2.AND.VOLW2R.GT.ZERO*VOLT(0,NY,NX))
     3.OR.(TKR22.GT.TFREEZ 
     4.AND.VOLI2R.GT.ZERO*VOLT(0,NY,NX)))THEN
      HFLFRX=VHCPR2*(TFREEZ-TKR22)*XNPRS
     2/(1.0+6.2913E-03*TFREEZ)
      IF(HFLFRX.LT.0.0)THEN
      HFLFR2=AMAX1(-333.0*DENSI*VOLI2R*XNPSRX,HFLFRX)
      ELSE
      HFLFR2=AMIN1(333.0*VOLW2R*XNPSRX,HFLFRX)
      ENDIF
      FLFR2=-HFLFR2/333.0
      FLFR2T=FLFR2T+FLFR2
      HFLFR2T=HFLFR2T+HFLFR2
      ELSE
      HFLFR2=0.0
      FLFR2=0.0
      ENDIF
C
C     VAPOR FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
C
C     VOLP02S,VOLP2R=snowpack,litter air content (m3)
C     VP0,VPR=snowpack,litter vapor concentration (m3 m-3)
C     VOLV2R,VOLP2R=litter vapor,air-filled content (m3)
C     FLVC=vapor-unlimited vapor flux (m3 t-1)
C     AVCNVR=snowpack-litter vapor conductance (m h-1) 
C     TK0X,TKXR=snowpack,litter interim temperature (K)
C     HFLXF=heat released by combustion in previous time step (MJ t-1)
C     PSISM1=litter matric water potential (MPa)
C     AREA=area of grid cell (m2)
C     FSNW,CVRD=snow,litter cover fraction 
C     FLVSRX=snow-litter vapor flux (m3 t-1) 
C     HFLVSRX=convective heat flux from snow-litter vapor flux (MJ t-1)
C     XNPQX=time step for snow-litter fluxes from ‘wthr.f’     
C 
      IF(VOLP02S.GT.ZEROS2(NY,NX)
     2.AND.VOLP2R.GT.ZEROS2(NY,NX))THEN
      VPR=AMAX1(0.0,VOLV2R/VOLP2R)
      FLVC=AVCNVR*(VP0-VPR)*AREA(3,NUM(NY,NX),NY,NX) 
     2*FSNW(NY,NX)*CVRD(NY,NX)*XNPQX
      IF(FLVC.GE.0.0)THEN
      FLVSRX=AMAX1(0.0,AMIN1(FLVC,VOLV02S*XNPXX))
      HFLVSRX=4.19*TK022*FLVSRX
      ELSE
      FLVSRX=AMIN1(0.0,AMAX1(FLVC,-VOLV2R*XNPXX))
      HFLVSRX=4.19*TKR22*FLVSRX
      ENDIF
      ELSE
      FLVSRX=0.0
      HFLVSRX=0.0
      ENDIF
C     WRITE(*,7768)'FLVSRX',I,J,NFZ,M,MM,NN,NX,NY,L
C    2,FLVSRX,FLVC,VOLV2R,AVCNVR,VP0,VPR,AREA(3,NUM(NY,NX),NY,NX) 
C    2,FSNW(NY,NX),CVRD(NY,NX),XNPQX 
C
C     HEAT FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
C
C     TKY=snow-litter equilibrium temperature (K)
C     TK022,TKR22=snow,litter interim temperature (K)
C     VHCPWM2,VHCPR2=current snowpack,litter layer heat capacity 
C        (MJ K-1)
C     HFLWX,HFLWC=snow-litter heat flux unlimited,limited by
C        temperature (MJ K-1)
C     ATCNDR=snow-litter thermal conductance (MJ h-1 K-1)
C     FSNW,CVRD=snow,litter cover fraction 
C     AREA=area of grid cell (m2)
C     HFLWSRX=snow-litter heat flux(MJ t-1)
C     XNPQX=time step for snows-litter fluxes from ‘wthr.f’(h t-1)
C
      TKY=(TK022*VHCPWM2(L,NY,NX)+TKR22*VHCPR2)
     2/(VHCPWM2(L,NY,NX)+VHCPR2) 
      HFLWX=(TK022-TKY)*VHCPWM2(L,NY,NX) 
      HFLWC=ATCNDR*(TK022-TKR22)*AREA(3,NUM(NY,NX),NY,NX) 
     2*FSNW(NY,NX)*CVRD(NY,NX)*XNPQX
      IF(HFLWC.GE.0.0)THEN
      HFLWSRX=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWSRX=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
C     WRITE(*,7769)'HFLWSRX',I,J,NFZ,M,MM,NN,NX,NY
C    2,HFLWSRX,HFLWX,HFLWC,ATCNDR,TK022,TKR22 
C    2,FSNW(NY,NX),CVRD(NY,NX),XNPQX,TKY,VHCPWM2(L,NY,NX)
C    3,HFLWSR,TK02(L) 
C
C     EVAPORATION-CONDENSATION IN SOIL SURFACE UNDER SNOW 
C
C     VPGV=soil saturation vapor concentration (m3 m-3)
C     TKS22,PSISVG=soil temperature,matric+osmotic potential (K,MPa)
C     FLVGX,FLVG2=soil evaporation-condensation unlimited,
C       limited by vapor (m3 t-1)
C     VOLW2G=soil water content (m3)
C     HFLVG2=soil latent heat flux (MJ t-1)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C
      VPGV=2.173E-03/TKS22 
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKS22))
     3*EXP(18.0*PSISVG/(8.3143*TKS22))
      FLVGX=(VOLV2G-VPGV*VOLP2G)
      FLVG2=AMAX1(FLVGX,-AMAX1(0.0,VOLW2G)*XNPSRX)
      HFLVG2=VAP*FLVG2
      FLVG2T=FLVG2T+FLVG2
      HFLVG2T=HFLVG2T+HFLVG2
C     WRITE(*,7769)'FLVG2',I,J,NFZ,M,MM,NN,NX,NY
C    2,FLVG2,FLVGX,FLVG2T,VOLW2G,XNPSRX,VOLV2G,VPGV*VOLP2G
C    3,VPGV,VOLP2G,FSNW(NY,NX),TKS22,TKS2,VOLV2(NUM(NY,NX),NY,NX)
C    4,HFLVG2,HFLVG2T
C    5,VHCPG2,(TKQG(M,NY,NX)-TKS2)*VHCPG2 
C
C     FREEZE-THAW IN SOIL SURFACE FROM NET CHANGE IN SOIL 
C     SURFACE HEAT STORAGE
C
C     TFREEZ=soil freezing temperature (K)
C     333.0=latent heat of freezing (MJ m-3)
C     PSISVG=soil water potential (MPa)
C     TKS22=soil temperature (K)
C     VOLW2G,VOLI2G=soil water,ice volume (m3)
C     VOLT=soil surface layer volume (m3)
C     HFLFGX,HFLFG2=soil freeze-thaw latent heat flux
C        unlimited,limited by water,ice (MJ t-1)
C     FLFG2=soil freeze-thaw flux (m3 t-1)
C     HFLFR,WFLFR=total litter freeze-thaw,latent heat flux
C     
      TFREEZ=-9.0959E+04/(PSISVG-333.0)
      IF((TKS22.LT.TFREEZ
     2.AND.VOLW2G.GT.ZERO*VOLT(NUM(NY,NX),NY,NX))
     3.OR.(TKS22.GT.TFREEZ 
     4.AND.VOLI2G.GT.ZERO*VOLT(NUM(NY,NX),NY,NX)))THEN
      HFLFGX=VHCPR2*(TFREEZ-TKS22)*XNPRS
     2/(1.0+6.2913E-03*TFREEZ)
      IF(HFLFGX.LT.0.0)THEN
      HFLFG2=AMAX1(-333.0*DENSI*VOLI2G*XNPSRX,HFLFGX)
      ELSE
      HFLFG2=AMIN1(333.0*VOLW2G*XNPSRX,HFLFGX)
      ENDIF
      FLFG2=-HFLFG2/333.0
      FLFG2T=FLFG2T+FLFG2
      HFLFG2T=HFLFG2T+HFLFG2
      ELSE
      HFLFG2=0.0
      FLFG2=0.0
      ENDIF
C
C     VAPOR FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE 
C
C     VOLP2R,VOLP2G=litter,soil air-filled volume (m3)
C     VPR,VP1=litter,soil vapor concentration (m3 m-3)
C     VOLV2R,VOLP2R=soil vapor,air-filled content (m3)
C     FLVC=vapor-unlimited vapor flux (m3 t-1)
C     AVCNVS=litter-soil vapor conductance (m h-1) 
C     FLVR1X=litter-soil vapor flux (m3 t-1)
C     HFLVR1X=convective heat of litter-soil vapor flux (m3 t-1) 
C     TKXR,TK1X=interim calculation of litter,soil temperatures
C     XNPRX=time step for litter fluxes from ‘wthr.f’(h t-1)     
C
      IF(VOLP2R.GT.ZEROS(NY,NX).AND.VOLP2G.GT.ZEROS(NY,NX))THEN
      VP1=AMAX1(0.0,VOLV2G/VOLP2G)
      FLVC=AVCNVS*(VPR-VP1)*AREA(3,NUM(NY,NX),NY,NX)*FSNW(NY,NX) 
     2*CVRD(NY,NX)*XNPQX 
      IF(FLVC.GE.0.0)THEN
      FLVR1X=AMAX1(0.0,AMIN1(FLVC,VOLV2R*XNPXX))
      HFLVR1X=4.19*TKR22*FLVR1X
      ELSE
      FLVR1X=AMIN1(0.0,AMAX1(FLVC,-VOLV2G*XNPXX))
      HFLVR1X=4.19*TKS22*FLVR1X
      ENDIF
      ELSE
      FLVR1X=0.0
      HFLVR1X=0.0 
      ENDIF
C     TKXR=TKXR-HFLVR1X/VHCP1(0,NY,NX)
C     TK1X=TK1X+HFLVR1X/VHCP1(NUM(NY,NX),NY,NX) 
C
C     HEAT FLUX BETWEEN SURFACE LITTER AND SOIL SURFACE 
C
C     TKY=litter-soil equilibrium temperature
C     VHCPR2,VHCPG2=litter,soil surface layer heat capacity (MJ K-1)
C     TKR22,TKS22=litter,soil surface interim temperature (K)
C     HFLWC,HFLWX=litter-soil heat flux unlimited,limited by heat 
C        (MJ t-1)
C     ATCNDS=litter-soil thermal conductance (MJ h-1 K-1)
C     AREA=area of grid cell (m2)
C     FSNW,CVRD=snow,litter cover fraction
C     HFLWR1X=litter-soil heat flux (MJ t-1)
C     XNPQX=time step for litter-soil fluxes from ‘wthr.f’ (h t-1)
C
      TKY=(TKR22*VHCPR2+TKS22*VHCPG2)/(VHCPR2+VHCPG2)
      HFLWX=(TKR22-TKY)*VHCPR2 
      HFLWC=ATCNDS*(TKR22-TKS22)*AREA(3,NUM(NY,NX),NY,NX) 
     2*FSNW(NY,NX)*CVRD(NY,NX)*XNPQX 
      IF(HFLWC.GE.0.0)THEN
      HFLWR1X=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWR1X=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
C
C     ACCUMULATE SNOW-LITTER, LITTER-SOIL HEAT FLUXES
C     WITHIN LONGER TIME STEP FOR SNOWPACK FLUX CALCULATIONS
C
C     FLVSR=snow-litter vapor flux (m3 t-1)
C     HFLVSR=convective heat flux from snow-litter vapor flux (MJ t-1)
C     HFLWSR=snow-litter heat flux (MJ t-1)
C     FLVR1=litter-soil vapor flux (m3 t-1)
C     HFLVR1=convective heat of litter-soil vapor flux (MJ t-1)
C     HFLWR1=litter-soil heat flux (MJ t-1)
C
      FLVSR=FLVSR+FLVSRX
      HFLVSR=HFLVSR+HFLVSRX
      HFLWSR=HFLWSR+HFLWSRX
      FLVR1=FLVR1+FLVR1X
      HFLVR1=HFLVR1+HFLVR1X
      HFLWR1=HFLWR1+HFLWR1X
C
C     UPDATE INTERMEDIATE WATER, VAPOR, ICE CONTENTS OF LOWEST SNOW
C     LAYER, SURFACE LITTER AND SOIL SURFACE LAYER
C
C     VOLW02S,VOLW02S,VOLP02S=snowpack water,vapor,air content (m3)
C     FLWQG,FLWQR=snowpack-soil,litter water flux (m3 t-1) 
C     FLVS1,FLVSRX=snowpack-soil,litter vapor flux (m3 t-1)
C     VOLS1,VOLA2,VOLA1=snowpack,litter,soil layer volume (m3)
C     VOLW2R,VOLV2R,VOLI2R,VOLP2R=litter water,vapor,ice,air content
C        (m3) 
C     FLVR2,FLFR2=litter evaporation-condensation,freeze-thaw (m3 t-1)
C     FLVSRX,FLVR1X=snowpack-litter,litter-soil vapor flux (m3 t-1)
C     VOLW2G,VOLV2G,VOLI2G,VOLP2G=soil water,vapor,ice,air content
C        (m3) 
C     FLWQGS=water flux to soil surface micropores (m3 t-1)
C     FLVG2,FLFG2=soil evaporation-condensation,freeze-thaw (m3 t-1)
C 
      VOLW02S=VOLW02S-(FLWQG+FLWQR)*XNPRS
      VOLV02S=VOLV02S-FLVS1*XNPRS-FLVSRX 
      VOLP02S=AMAX1(0.0,VOLS1(L,NY,NX)-VOLS02S-VOLI02S-VOLW02S)
      VOLW2R=VOLW2R+FLWQR*XNPRS+FLVR2+FLFR2
      VOLV2R=VOLV2R+FLVSRX-FLVR1X-FLVR2 
      VOLI2R=VOLI2R-FLFR2/DENSI
      VOLP2R=AMAX1(0.0,VOLA2-VOLW2R-VOLI2R)
      VOLW2G=VOLW2G+FLWQGS*XNPRS+FLVG2+FLFG2 
      VOLV2G=VOLV2G+FLVR1X+FLVS1*XNPRS-FLVG2 
      VOLI2G=VOLI2G-FLFG2/DENSI
      VOLP2G=AMAX1(0.0,VOLA1(NUM(NY,NX),NY,NX)-VOLW2G-VOLI2G)
7768  FORMAT(A8,9I4,40E14.6) 
7769  FORMAT(A8,8I4,40E14.6) 
4000  CONTINUE
      ENDIF
C
C     SNOW WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
C     FOR LATER UPDATES TO STATE VARIABLES
C        units:FL*=m3 t-1,HF*=MJ t-1
C
C     FLWLT,FLWLW=total,accumulated water flux to soil micropores 
C     FLVLT,FLWLV=total,accumulated vapor flux to soil 
C     FLWLXW,FLWHLW=total,accumd snow-soil micropore,macropore water 
C     HFLWLT,HFLWLW=total,accumulated snow+litter heat flux to soil
C     FLWRT,FLWRLW=total,accumulated snow+soil water flux to litter
C     FLVRT,FLVRLW=total,accumulated snow+soil vapor flux to litter
C     HFLWRT,HFLWRLW=total,accumulated snow+soil heat flux to litter
C     FLQRM,FLQSM,FLQHM=total water flux to litter,soil
C        micropore,macropore
C     FLSW,FLSWH,FLSWR=water flux from lowest snow layer to soil
C        macropore,micropore,litter
C     FLSV,FLSVR=vapor flux from lowest snow layer to soil,litter
C     HFLSW,HFLSWR=heat flux from lowest snow layer to soil,litter
C     FLWQGS,FLWQGH=snowpack-soil micro-,macro-pore water flux
C     FLWQR=snowpack-litter water flux  
C     FLVS1,FLVR1=snowpack-soil,litter vapor flux 
C     HFLWQG,HFLWQR=convective heat fluxes to soil,litter
C     HFLVS1,HFLWS1=snowpack-soil convective,conductive heat fluxes
C     HFLVSR,HFLWSR=snowpack-litter convective,conductive heat fluxes
C     HFLVR1=convective heat of litter-soil vapor flux
C     HFLWR1=litter-soil heat flux
C
      FLWLT=FLWQGS 
      FLWLW=FLWLW+FLWLT 
      FLVLT=FLVS1+FLVR1 
      FLWLV=FLWLV+FLVLT 
      FLWLXW=FLWLXW+FLWQGS 
      FLWHLW=FLWHLW+FLWQGH
      HFLWLT=HFLWQG+HFLVS1+HFLWS1+HFLVR1+HFLWR1 
      HFLWLW=HFLWLW+HFLWLT
      FLWRT=FLWQR 
      FLWRLW=FLWRLW+FLWRT
      FLVRT=FLVSR-FLVR1
      FLVRLW=FLVRLW+FLVRT
      HFLWRT=HFLWQR+HFLVSR+HFLWSR-HFLVR1-HFLWR1 
      HFLWRLW=HFLWRLW+HFLWRT
      FLWVLW=0.0
      FLQRM(M,NY,NX)=FLQRM(M,NY,NX)+FLWQR
      FLQSM(M,NY,NX)=FLQSM(M,NY,NX)+FLWQGS
      FLQHM(M,NY,NX)=FLQHM(M,NY,NX)+FLWQGH
      FLSW(L,NY,NX)=FLSW(L,NY,NX)+FLWLT 
      FLSV(L,NY,NX)=FLSV(L,NY,NX)+FLVLT 
      FLSWH(L,NY,NX)=FLSWH(L,NY,NX)+FLWQGH 
      HFLSW(L,NY,NX)=HFLSW(L,NY,NX)+HFLWLT 
      FLSWR(L,NY,NX)=FLSWR(L,NY,NX)+FLWRT 
      FLSVR(L,NY,NX)=FLSVR(L,NY,NX)+FLVRT 
      HFLSWR(L,NY,NX)=HFLSWR(L,NY,NX)+HFLWRT
      TFLWSX=FLW0S(L,NY,NX)
      TFLWWX=FLW0W(L,NY,NX)-FLWRT-FLWLT-FLWQG
      TFLWVX=FLW0V(L,NY,NX)-FLVRT-FLVLT 
      TFLWIX=FLW0I(L,NY,NX)
C     THFLWWX=HFLW0W(L,NY,NX)-HFLWRT-HFLWLT
      THFLWWX=HFLW0W(L,NY,NX)-HFLVSR-HFLWSR-HFLVS1-HFLWS1
     2-HFLWQG-HFLWQR
      TFLWS(L,NY,NX)=TFLWS(L,NY,NX)+TFLWSX
      TFLWW(L,NY,NX)=TFLWW(L,NY,NX)+TFLWWX 
      TFLWV(L,NY,NX)=TFLWV(L,NY,NX)+TFLWVX 
      TFLWI(L,NY,NX)=TFLWI(L,NY,NX)+TFLWIX
      THFLWW(L,NY,NX)=THFLWW(L,NY,NX)+THFLWWX
C
C     UPDATE INTERMEDIATE STATE VARIABLES FOR SNOWPACK WATER,VAPOR
C
C     VOLW02,VOLV02,VOLI02,VOLP02,VOLS1=snowpack
C        water,vapor,ice,air,layer volume (m3)
C     FLWLT,FLWRT=total water flux to soil micropores,litter (m3)
C     FLWQGH= total water flux to soil macropore (m3)
C     FLVLT,FLVRT=total vapor flux to soil,litter (m3)
C
      VOLW02(L,NY,NX)=VOLW02(L,NY,NX)-FLWLT-FLWRT-FLWQGH 
      VOLV02(L,NY,NX)=VOLV02(L,NY,NX)-FLVLT-FLVRT
      VOLP02(L,NY,NX)=AMAX1(0.0,VOLS1(L,NY,NX)-VOLS02(L,NY,NX) 
     2-VOLI02(L,NY,NX)-VOLW02(L,NY,NX))
C     IF(I.GT.270)THEN
C     WRITE(*,7752)'FLWLW',I,J,NFZ,M,MM,NX,NY,L 
C    2,FLWLW,FLWLT,FLWQGS,FLVS1,FLVR1,FLSW(L,NY,NX)
C    2,FLWRLW,FLWRT,FLWQR,FLVSR,FLVR1,FLWQX,FLWQG,FLSWR(L,NY,NX)
C    2,FLVRLW,FLVRT,FLVSR,FLVR1,FLVR2T,VOLV2(0,NY,NX)
C    2,FLWQX,FLWQGX,BARE(NY,NX),VOLW02(L,NY,NX),VOLS02(L,NY,NX)
C    2,HFLWLW,HFLWLT,HFLWQG,HFLVS1,HFLWS1,HFLVR1,HFLWR1
C    3,HFLSW(L,NY,NX)
C    3,HFLWRLW,HFLWRT,HFLWQR,HFLVSR,HFLWSR,HFLVR1,HFLWR1
C    3,THFLWW(L,NY,NX),XHFLF0(L,NY,NX),XHFLV0(L,NY,NX),HFLSWR(L,NY,NX)
C    4,THFLWW(L,NY,NX),THFLWWX,HFLW0W(L,NY,NX),HFLWRT,HFLWLT
C    5,VOLV02(L,NY,NX),FLVLT,FLVRT,FLVSR,FLVR1 
7752  FORMAT(A8,8I4,40E14.6)
C     ENDIF
      ICHKL=1
      ENDIF
      ENDIF
      ENDIF
C
C     EVAPORATION-CONDENSATION IN SNOWPACK
C
C     VOLS02,VOLW02,VOLI02,VOLV02=snow,water,ice,vapor volume (m3)
C     VPSV=equilibrium vapor concentration (m3 m-3)
C     WFLVW2,WFLVS2=condensation(+ve) or evaporation(-ve)
C        from snowpack water,snow (m3 t-1)
C     HFLVX=latent heat of condensation,evaporation from
C        snowpack water+snow (MJ t-1)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C     XWFLVW,XWFLVS=aggregated condensation(+ve), evaporation(-ve) from
C        water,snow used in ‘redist.f’ (MJ t-1)
C     XHFLV0=aggregated latent heat of condensation,evaporation from
C        water+snow used in ‘redist.f’ (MJ t-1)
C     XNPSX=time step for snowpack fluxes from ‘wthr.f’ (h t-1)
C
      IF(L.GT.1)THEN
      VPSV=2.173E-03/TK02(L)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK02(L)))
      IF(VOLP02(L,NY,NX).GT.ZEROS(NY,NX))THEN
      WFLVT=VOLV02(L,NY,NX)-VPSV*VOLP02(L,NY,NX)
      WFLVW2=AMAX1(WFLVT,-AMAX1(0.0,VOLW02(L,NY,NX))*XNPSX)
      WFLVX=AMIN1(0.0,WFLVT-WFLVW2)
      WFLVS2=AMAX1(WFLVX,-AMAX1(0.0,VOLS02(L,NY,NX))*XNPSX)
      HFLVX=VAP*WFLVW2+VAPS*WFLVS2
      ELSE
      WFLVT=0.0
      WFLVX=0.0
      WFLVW2=0.0
      WFLVS2=0.0
      HFLVX=0.0
      ENDIF
      WFLVW(L,NY,NX)=WFLVW(L,NY,NX)+WFLVW2
      WFLVS(L,NY,NX)=WFLVS(L,NY,NX)+WFLVS2 
      HFLV0(L,NY,NX)=HFLV0(L,NY,NX)+HFLVX      
      XWFLVW(L,NY,NX)=XWFLVW(L,NY,NX)+WFLVW2      
      XWFLVS(L,NY,NX)=XWFLVS(L,NY,NX)+WFLVS2      
      XHFLV0(L,NY,NX)=XHFLV0(L,NY,NX)+HFLVX
C
C     UPDATE INTERMEDIATE STATE VARIABLES FOR SNOWPACK WATER,VAPOR
C
C     VOLS02,VOLW02,VOLI02,VOLV02,VOLP02=snow,water,ice,vapor,air
C        volume (m3)
C     WFLVW2,WFLVS2=condensation(+ve) or evaporation(-ve)
C        from snowpack water,snow (m3 t-1)
C
      VOLS02(L,NY,NX)=VOLS02(L,NY,NX)+WFLVS2 
      VOLW02(L,NY,NX)=VOLW02(L,NY,NX)+WFLVW2 
      VOLV02(L,NY,NX)=VOLV02(L,NY,NX)-WFLVS2-WFLVW2 
      VOLP02(L,NY,NX)=AMAX1(0.0,VOLS1(L,NY,NX)-VOLS02(L,NY,NX) 
     2-VOLI02(L,NY,NX)-VOLW02(L,NY,NX))
C     WRITE(*,7761)'WFLVWL',I,J,NFZ,M,MM,NX,NY,L
C    2,WFLVT,VOLV0X,VP1
C    2,VOLP02(L,NY,NX),VP1*VOLP02(L,NY,NX),WFLVW2,VOLW02(L,NY,NX)
C    3,WFLVS2,WFLVS2,VOLS02(L,NY,NX),VOLV02(L,NY,NX),TK02(L)     
7761  FORMAT(A8,8I4,80E14.6)
      ENDIF
C
C     FREEZE-THAW IN SNOWPACK LAYER FROM NET CHANGE IN SNOWPACK
C     HEAT STORAGE
C
C     VHCPWMM,VHCPWMX,VHCPWX=previous,current,minimum heat capacity 
C        (MJ K-1)
C     VOLS02,VOLW02,VOLV02,VOLI02,VOLS=snow,water,vapor,ice,
C        total snowpack volume (m3)
C     DENSI=ice density (Mg m-3)
C     TK02=snowpack interim temperature (K)
C     HFLF1=unlimited latent heat flux from freeze-thaw (MJ t-1)
C     FVOLS0,FVOLI0=fractions of total water in water,ice 
C     HFLF0X=source-limited latent heat flux from freeze-thaw (MJ t-1)
C     WFLFSX,WFLFIX=freeze-thaw changes in water,ice (m3 t-1)
C     WFLFS,WFLFI=accumulated freeze-thaw (m3 t-1)
C     HFLF0=accumulated latent heat flux from freeze-thaw (MJ t-1)
C     XWFLFS,XWFLFI=aggregated freeze-thaw for ‘redist.f’ (m3 t-1) 
C     XTHAWW=aggregated latent heat flux from freeze-thaw for
C        ‘redist.f’ (MJ t-1) 
C
      VHCPWMX=2.095*VOLS02(L,NY,NX)+4.19*(VOLW02(L,NY,NX)
     2+VOLV02(L,NY,NX))+1.9274*VOLI02(L,NY,NX)
      IF(VHCPWMX.GT.VHCPWX(NY,NX))THEN
      IF((TK02(L).LT.273.15
     2.AND.VOLW02(L,NY,NX).GT.ZERO*VOLS(NY,NX))
     3.OR.(TK02(L).GT.273.15 
     4.AND.VOLI02(L,NY,NX)+VOLS02(L,NY,NX).GT.ZERO*VOLS(NY,NX)))THEN
      HFLF1=VHCPWMX*(273.15-TK02(L))/2.7185*XNPS
      IF(HFLF1.LT.0.0)THEN
      TVOLWS=VOLS02(L,NY,NX)+VOLI02(L,NY,NX)*DENSI 
      IF(TVOLWS.GT.ZEROS2(NY,NX))THEN
      FVOLS0=VOLS02(L,NY,NX)/TVOLWS
      FVOLI0=VOLI02(L,NY,NX)*DENSI/TVOLWS
      ELSE
      FVOLS0=0.0
      FVOLI0=0.0
      ENDIF
      HFLF0X=AMAX1(-333.0*TVOLWS*XNPSX,HFLF1)
      WFLFSX=-HFLF0X*FVOLS0/333.0
      WFLFIX=-HFLF0X*FVOLI0/333.0
      ELSE
      FVOLS0=0.0
      FVOLI0=0.0
      HFLF0X=AMIN1(333.0*VOLW02(L,NY,NX)*XNPSX,HFLF1)
      WFLFSX=0.0
      WFLFIX=-HFLF0X/333.0
      ENDIF
      ELSE
      HFLF0X=0.0
      WFLFSX=0.0
      WFLFIX=0.0
      ENDIF
      WFLFS(L,NY,NX)=WFLFS(L,NY,NX)+WFLFSX
      WFLFI(L,NY,NX)=WFLFI(L,NY,NX)+WFLFIX
      HFLF0(L,NY,NX)=HFLF0(L,NY,NX)+HFLF0X
      XWFLFS(L,NY,NX)=XWFLFS(L,NY,NX)+WFLFSX 
      XWFLFI(L,NY,NX)=XWFLFI(L,NY,NX)+WFLFIX 
      XHFLF0(L,NY,NX)=XHFLF0(L,NY,NX)+HFLF0X
C
C     UPDATE INTERMEDIATE STATE VARIABLES FOR SNOWPACK WATER,VAPOR
C
C     VOLS02,VOLW02,VOLI02,VOLV02,VOLP02=snow,water,ice,vapor,air
C        volume (m3)
C     WFLFSX,WFLFIX=freeze-thaw changes in water,ice (m3 t-1)
C
      VOLS02(L,NY,NX)=VOLS02(L,NY,NX)-WFLFSX 
      VOLW02(L,NY,NX)=VOLW02(L,NY,NX)+WFLFSX+WFLFIX 
      VOLI02(L,NY,NX)=VOLI02(L,NY,NX)-WFLFIX/DENSI
      VOLP02(L,NY,NX)=AMAX1(0.0,VOLS1(L,NY,NX)-VOLS02(L,NY,NX) 
     2-VOLI02(L,NY,NX)-VOLW02(L,NY,NX))
      ELSE
      HFLF0X=0.0
      WFLFSX=0.0
      WFLFIX=0.0
      ENDIF
C
C     INTERNAL SNOWPACK TEMPERATURE
C
C     VHCPWM2,VHCPWX=snowpack heat capacity,minimum (MJ K-1)
C     TK02=snowpack interim temperature (K)
C     VOLS02,VOLW02,VOLI02,VOLV02=snow,water,ice,vapor volume (m3)
C     THFLWWX=convective heat flux from net snow,water,ice transfer 
C        (MJ t-1)
C     HFLF0X=source-limited latent heat flux from freeze-thaw (MJ t-1)
C     HFLVX=latent heat of condensation,evaporation from
C        snowpack water+snow (MJ t-1)
C
      ENGY0=VHCPWM2(L,NY,NX)*TK02(L)
      VHCPWM2(L,NY,NX)=2.095*VOLS02(L,NY,NX)
     2+4.19*(VOLW02(L,NY,NX)+VOLV02(L,NY,NX))
     3+1.9274*VOLI02(L,NY,NX)
      IF(VHCPWM2(L,NY,NX).GT.VHCPWX(NY,NX))THEN
      TK02(L)=(ENGY0+THFLWWX+HFLF0X+HFLVX)/VHCPWM2(L,NY,NX)
C     IF(IYRC.EQ.1979)THEN
C     WRITE(*,7758)'TK02',I,J,NFZ,M,MM,NX,NY,L
C    2,TK02(L),ENGY0,THFLWWX,HFLW0W(L,NY,NX) 
C    3,HFLVSR,HFLWSR,HFLVS1,HFLWS1
C    2,HFLWQG,HFLWQR,VHCPWM2(L,NY,NX)
C    4,TKAM(NY,NX),TKR2,TKS2,TK0(L,NY,NX),TKW(L,NY,NX) 
C    2,TKR2,TKS2,TFLWVX,WFLVS2,WFLVW2 
C    2,WFLFSX,WFLFIX,WFLFS(L,NY,NX),WFLFI(L,NY,NX),XWFLFS(L,NY,NX) 
C    3,TFLWSX,TFLWWX,TFLWIX
C    4,HFLF0X,HFLVX,VHCPWM2(L,NY,NX)
C    4,HFLW0W(L,NY,NX)
C    3,XFLWS(L,NY,NX),XFLWW(L,NY,NX),XFLWI(L,NY,NX)
C    3,VOLS02(L,NY,NX),VOLW02(L,NY,NX),VOLV02(L,NY,NX),VOLI02(L,NY,NX)
C    4,VOLW0(L,NY,NX),VOLWSL(L,NY,NX)
C     ENDIF
      ELSEIF(L.EQ.1)THEN
      TK02(L)=TKQG(M,NY,NX)
      ELSE
      TK02(L)=TK02(L-1)
      ENDIF
C     IF(L.EQ.5)THEN
C     WRITE(*,7758)'HFLF0',I,J,NFZ,M,MM,NX,NY,L,TK02(L) 
C    4,HFLF0(L,NY,NX),WFLFS(L,NY,NX),WFLFI(L,NY,NX)
C    4,TFLXV(L,NY,NX),WFLXV(L,NY,NX) 
C    4,XHFLXW(L,NY,NX),XWFLFS(L,NY,NX),XWFLFI(L,NY,NX)
C    2,TK0X,TKW(L,NY,NX),VHCPWMX,HFLF1,VOLS02(L,NY,NX)
C    3,VOLW02(L,NY,NX),VOLW02(L,NY,NX),TFLWS(L,NY,NX),TFLWW(L,NY,NX)
C    4,TFLWI(L,NY,NX),FVOLS0,FVOLI0
C    5,TFLWW(L,NY,NX),THFLWW(L,NY,NX),FLW0W(L,NY,NX)
7758  FORMAT(A8,8I4,60E14.6) 
C     ENDIF
C
C     ACCUMULATE SNOWPACK FLUXES TO LONGER TIME STEP FOR
C     LITTER, SOIL FLUX CALCULATIONS
C
C     XFLWS,XFLWW,XFLWV,XFLWI=aggregated snow,water,vapor,ice transfer
C        used in ‘redist.f’ (m3 t-1)
C     XHFLWW=aggregated convective heat flux from snow,water,ice
C        transfer used in ‘redist.f’ (MJ t-1) 
C
      XFLWS(L,NY,NX)=XFLWS(L,NY,NX)+FLW0S(L,NY,NX)
      XFLWW(L,NY,NX)=XFLWW(L,NY,NX)+FLW0W(L,NY,NX)
      XFLWV(L,NY,NX)=XFLWV(L,NY,NX)+FLW0V(L,NY,NX)
      XFLWI(L,NY,NX)=XFLWI(L,NY,NX)+FLW0I(L,NY,NX)
      XHFLWW(L,NY,NX)=XHFLWW(L,NY,NX)+HFLW0W(L,NY,NX) 
9880  CONTINUE
C
C     WATER, VAPOR, ICE OF SURFACE LITTER AND SOIL
C     SURFACE BENEATH SNPOWPACK
C
C     WFLVR,HFLVR=litter evaporation-condensation,
C        latent heat flux (m3 t-1,MJ t-1)
C     WFLFR,HFLFR=litter freeze-thaw,latent 
C        heat flux (m3 t-1,MJ t-1)
C     WFLVL,HFLVL=soil evaporation-condensation,
C        latent heat flux (m3 t-1,MJ t-1)
C     WFLFL,HFLFL=soil freeze-thaw,latent 
C        heat flux (m3 t-1,MJ t-1)
C     FLVR2T,HFLVR2T=litter evaporation-condensation,
C        latent heat flux (m3 t-1,MJ t-1)
C     FLFR2T,HFLFR2T=litter freeze-thaw,latent 
C        heat flux (m3 t-1,MJ t-1)
C     FLVG2T,HFLVG2T=soil evaporation-condensation,
C        latent heat flux (m3 t-1,MJ t-1)
C     FLFG2T,HFLFG2T=soil freeze-thaw,latent 
C        heat flux (m3 t-1,MJ t-1)
C     VOLW2,VOLV2,VOLI2,VOLP2=litter(0,),soil(NUM,)
C        water,vapor,ice,air content (m3)
C
      WFLVR(NY,NX)=WFLVR(NY,NX)+FLVR2T
      HFLVR(NY,NX)=HFLVR(NY,NX)+HFLVR2T
      WFLFR(NY,NX)=WFLFR(NY,NX)+FLFR2T
      HFLFR(NY,NX)=HFLFR(NY,NX)+HFLFR2T
      WFLVL(NUM(NY,NX),NY,NX)=WFLVL(NUM(NY,NX),NY,NX)+FLVG2T
      HFLVL(NUM(NY,NX),NY,NX)=HFLVL(NUM(NY,NX),NY,NX)+HFLVG2T
      WFLFL(NUM(NY,NX),NY,NX)=WFLFL(NUM(NY,NX),NY,NX)+FLFG2T
      HFLFL(NUM(NY,NX),NY,NX)=HFLFL(NUM(NY,NX),NY,NX)+HFLFG2T
C     WRITE(*,5352)'WFLVRS',I,J,NFZ,M,MM,NX,NY
C    2,WFLVR(NY,NX),HFLVR(NY,NX),WFLFR(NY,NX),HFLFR(NY,NX),FLVR2T 
C     WRITE(*,5352)'WFLVGS',I,J,NFZ,M,MM,NX,NY
C    2,WFLVL(NUM(NY,NX),NY,NX),HFLVL(NUM(NY,NX),NY,NX)
C    3,WFLFL(NUM(NY,NX),NY,NX),HFLFL(NUM(NY,NX),NY,NX)
C    4,FLVG2T,FLFG2T 
5352  FORMAT(A8,7I4,30E14.6)
      VOLW2(0,NY,NX)=VOLW2(0,NY,NX)+FLWRT+FLVR2T+FLFR2T
      VOLV2(0,NY,NX)=VOLV2(0,NY,NX)+FLVRT-FLVR2T 
      VOLI2(0,NY,NX)=VOLI2(0,NY,NX)-FLFR2T/DENSI
      VOLP2(0,NY,NX)=AMAX1(0.0,VOLA2-VOLW2(0,NY,NX)
     2-VOLI2(0,NY,NX))
      VOLW2(NUM(NY,NX),NY,NX)=VOLW2(NUM(NY,NX),NY,NX)+FLWLT
     2+FLVG2T+FLFG2T 
      VOLV2(NUM(NY,NX),NY,NX)=VOLV2(NUM(NY,NX),NY,NX)+FLVLT-FLVG2T 
      VOLI2(NUM(NY,NX),NY,NX)=VOLI2(NUM(NY,NX),NY,NX)-FLFG2T/DENSI
      VOLP2(NUM(NY,NX),NY,NX)=AMAX1(0.0,VOLA1(NUM(NY,NX),NY,NX)
     2-VOLW2(NUM(NY,NX),NY,NX)-VOLI2(NUM(NY,NX),NY,NX))
C
C     LITTER AND SOIL SURFACE TEMPERATURES
C
C     VHCPR2,VHCPG2=litter,soil heat capacity (MJ K-1)
C     TKR2,TKS2=litter,soil temperature (K)
C     ORGC,ORGCC=litter C, charcoal content (g)
C     VOLW2,VOLV2,VOLI2,VOLP2=litter(0,),soil(NUM,)
C        water,vapor,ice,air content (m3)
C     VOLWH1,VOLIH1=soil macropore water,ice content (m3)
C     FSNW=snow cover fraction
C     HFLWRT,HFLWLT=total heat flux to litter,soil (MJ t-1)
C     HFLVR2T,HFLVG2T =litter,soil evaporation-condensation
C        latent heat flux (MJ t-1)
C     HFLFR2T,HFLFG2T=litter,soil freeze-thaw latent 
C        heat flux (MJ t-1)
C
      ENGYR=VHCPR2*TKR2
      VHCPR2=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))*FSNW(NY,NX)
     2+4.19*(VOLW2(0,NY,NX)+VOLV2(0,NY,NX))
     2+1.9274*VOLI2(0,NY,NX)
      IF(VHCPR2.GT.VHCPRXS(NY,NX))THEN
      TKR2=(ENGYR+HFLWRT+HFLVR2T+HFLFR2T)/VHCPR2
      ELSE
      TKR2=TK1(0,NY,NX)
      ENDIF
      ENGY2=VHCPG2*TKS2
      VHCPG2=VHCM(NUM(NY,NX),NY,NX)
     2+4.19*(VOLW2(NUM(NY,NX),NY,NX)+VOLV2(NUM(NY,NX),NY,NX)
     2+VOLWH1(NUM(NY,NX),NY,NX))
     2+1.9274*(VOLI2(NUM(NY,NX),NY,NX)+VOLIH1(NUM(NY,NX),NY,NX))
      IF(VHCPG2.GT.VHCPRX(NY,NX)
     2.AND.VOLX(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      TKS2=(ENGY2+HFLWLT+HFLVG2T+HFLFG2T)/VHCPG2
      ELSE
      TKS2=TK1(NUM(NY,NX),NY,NX)
      ENDIF
C     IF(I.GT.270)THEN
C     WRITE(*,7759)'VOLV2S',I,J,NFZ,M,MM,NX,NY
C    2,VOLV2(0,NY,NX),FLVRT,FLVR2T,VOLV1(0,NY,NX)
C     WRITE(*,7759)'VOLW2S',I,J,NFZ,M,MM,NX,NY
C    2,VOLW2(0,NY,NX),FLWRT,FLVR2T,FLFR2T,FLWRLW,FLWRLG
C     WRITE(*,7759)'TKR2S',I,J,NFZ,M,MM,NX,NY
C    2,TKR2,ENGYR,HFLWRT,HFLVR2T,HFLFR2T,VHCPR2,FSNW(NY,NX)
C    3,VOLW2(0,NY,NX),VOLV2(0,NY,NX),VOLI2(0,NY,NX)
C    4,FLWRT,FLVR2T,FLFR2T,VOLW1(0,NY,NX)
C    5,TK1(0,NY,NX),TKS(0,NY,NX),VHCPRXS(NY,NX) 
C     WRITE(*,7759)'TKS2S',I,J,NFZ,M,MM,NX,NY
C    2,TKS2,ENGY2,HFLWLT,HFLVG2T,HFLFG2T,VHCPG2,TK1(NUM(NY,NX),NY,NX)
C     3,TKS(NUM(NY,NX),NY,NX)
C     ENDIF
9860  CONTINUE
3000  CONTINUE
      ENDIF
C
C     ENERGY EXCHANGE AT EXPOSED LITTER AND SOIL SURFACE 
C
C     FSNX=snow-free surface fraction
C     BKDS=bulk density (0=pond) (Mg m-3)
C     VHCP1,VHCPRX=current,minimum soil heat capacities (MJ K-1)
C
      IF(FSNX(NY,NX).GT.0.0.AND.(BKDS(NUM(NY,NX),NY,NX).GT.ZERO 
     2.OR.VHCP1(NUM(NY,NX),NY,NX).GT.VHCPRX(NY,NX)))THEN
C
C     PHYSICAL AND HYDRAULIC PROPERTIES OF EXPOSESOIL SURFACE INCLUDING
C     AIR AND WATER-FILLED POROSITY, AND WATER POTENTIAL USED IN
C     FLUX CALCULATIONS
C
C     VOLW2,VOLV2,VOLI2,VOLA2,VOLP2=exposed litter
C        water,vapor,ice,porosity,air content (m3)
C     TKR2,TKS2=exposed litter,soil temperature (K)
C     VHCPR2,VHCPG2=exposed litter,soil heat capacity (MJ K-1)
C     BKVL=soil mass (0=pond)  (Mg m-3)
C     ORGC,ORGCC=litter C, charcoal g)
C     THETW1,THETZ=current,minimum water concentration (m3 m-3)
C     POROS=porosity (m3 m-3)
C     FC,WP=water contents at field capacity,wilting point (m3 m-3)
C     FCL,WPL=log FC,WP
C     FCD,PSD=FCL-WPL,log(POROS)-FCL
C     PSISM1,PSISX,PSISE=soil matric,minimum (set in ‘starts.f’),
C        saturation potential(MPa)
C     PSIMX,PSIMN,PSIMS=log water potential at FC,WP,POROS
C     PSISD,PSIMD=PSIMX-PSIMS,PSIMN-PSIMX
C     SRP=parameter for deviation from linear log-log water retention 
C        function from ‘hour1.f’
C
      VOLW2(0,NY,NX)=VOLW1(0,NY,NX)*FSNX(NY,NX)
      VOLV2(0,NY,NX)=VOLV1(0,NY,NX)*FSNX(NY,NX)
      VOLI2(0,NY,NX)=VOLI1(0,NY,NX)*FSNX(NY,NX)
      VOLA2=VOLA1(0,NY,NX)*FSNX(NY,NX)
      VOLP2(0,NY,NX)=AMAX1(0.0,VOLA2-VOLW2(0,NY,NX)
     2-VOLI2(0,NY,NX))
      TKR2=TK1(0,NY,NX)
      TKS2=TK1(NUM(NY,NX),NY,NX)
      VHCPR2=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))*FSNX(NY,NX)
     2+4.19*(VOLW2(0,NY,NX)+VOLV2(0,NY,NX))
     2+1.9274*VOLI2(0,NY,NX)
      VHCPG2=VHCP1(NUM(NY,NX),NY,NX)
      IF(BKVL(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      THETW1=AMAX1(THETZ(NUM(NY,NX),NY,NX)
     2,AMIN1(POROS(NUM(NY,NX),NY,NX)
     2,VOLW1(NUM(NY,NX),NY,NX)/VOLY(NUM(NY,NX),NY,NX)))
      IF(THETW1.LT.FC(NUM(NY,NX),NY,NX))THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSISX,-EXP(PSIMX(NY,NX)
     2+((FCL(NUM(NY,NX),NY,NX)-LOG(THETW1))
     3/FCD(NUM(NY,NX),NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETW1.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN 
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX)
     2+((AMAX1(0.0,(PSL(NUM(NY,NX),NY,NX)-LOG(THETW1)))
     3/PSD(NUM(NY,NX),NY,NX))**SRP(NUM(NY,NX),NY,NX)*PSISD(NY,NX)))
      ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
C
C     ENERGY EXCHANGE AT POND SURFACE IF EXPOSED UNDER SNOWPACK
C
C     THETIX,THETWX=soil ice,water concentration (m3 m-3)
C     FCI,WPI=ice field capacity,wilting point (m3 m-3)
C     PSL,FCL,WPL=log POROS0,FC,WP
C     FCD,PSD=FCL-WPL,PSL-FCL
C     PSISM1=matric water potential (MPa)
C     POROS=porosity (m3 m-3)
C     PSISM1,PSISX,PSISE=soil matric,minimum (set in ‘starts.f’),
C        saturation potential(MPa)
C
      ELSEIF(VOLX(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      FCX=FCI*THETIX(NUM(NY,NX),NY,NX)
      WPX=WPI*THETIX(NUM(NY,NX),NY,NX)
      FCLX=LOG(FCX)
      WPLX=LOG(WPX)
      PSDX=PSL(NUM(NY,NX),NY,NX)-FCLX
      FCDX=FCLX-WPLX
      IF(THETWX(NUM(NY,NX),NY,NX).LT.FCX)THEN
      PSISM1(NUM(NY,NX),NY,NX)=AMAX1(PSISX,-EXP(PSIMX(NY,NX)
     2+((FCLX-LOG(THETWX(NUM(NY,NX),NY,NX)))
     3/FCDX*PSIMD(NY,NX))))
      ELSEIF(THETWX(NUM(NY,NX),NY,NX)
     2.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN 
      PSISM1(NUM(NY,NX),NY,NX)=-EXP(PSIMS(NY,NX)
     2+((AMAX1(0.0,(PSL(NUM(NY,NX),NY,NX)
     2-LOG(THETWX(NUM(NY,NX),NY,NX))))
     3/PSDX)*PSISD(NY,NX)))
      ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
      ELSE
      THETW1=POROS(NUM(NY,NX),NY,NX)
      PSISM1(NUM(NY,NX),NY,NX)=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
C     IF(NX.EQ.4.AND.NY.EQ.5)THEN
C     WRITE(*,3232)'PSISM1',I,J,M,NX,NY,NUM(NY,NX)
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
C     VOLW1,VOLI1=water,ice volume in micopores (m3)
C     VOLWH1,VOLIH1=water,ice volume in macopores (m3)
C     ALBG,ALBS=albedo of ground,soil surface
C     BKVL=soil mass (Mg)
C     RFLX0,THRYG=net SW,sky LW radiation at soil surface (MJ t-1)
C     THRMXS=longwave radiation emitted from soil surface (MJ t-1)
C     EMMGS=soil surface emissivity
C     THRMCS,THRMDS=net LW exchange between soil and canopy, standing
C        dead surfaces (MJ t-1)
C     TKC,TKD,TK1=canopy,standing dead,soil surface temperature (K)
C     RFLXG=net radiation (MJ t-1)
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
      THRMXS=EMMGS(NY,NX)*TK1(NUM(NY,NX),NY,NX)**4
      RFLXG=RFLX0-THRMXS
      DO 910 NZ=1,NP(NY,NX)
      THRMCS=EMMCS(NY,NX)*(TKC(NZ,NY,NX)**4-TK1(NUM(NY,NX),NY,NX)**4)
     2*FRADP(NZ,NY,NX)
      THRMDS=EMMCS(NY,NX)*(TKD(NZ,NY,NX)**4-TK1(NUM(NY,NX),NY,NX)**4)
     2*FRADQ(NZ,NY,NX)
      RFLXG=RFLXG+THRMCS+THRMDS
C     IF(NX.EQ.5)THEN
C     WRITE(*,1108)'RFLXG',I,J,NFZ,M,NX,NY,NZ 
C    3,RFLXG,RFLX0,THRMXS,THRMCS,THRMDS
C    4,TKC(NZ,NY,NX),TKD(NZ,NY,NX),TK1(NUM(NY,NX),NY,NX)
C    5,FRADP(NZ,NY,NX),FRADQ(NZ,NY,NX)
1108  FORMAT(A8,7I4,60E12.4)
C     ENDIF
910   CONTINUE
C
C     EVAPORATION-CONDENSATION IN SOIL SURFACE
C
C     VOLV2,VOLPM=soil vapor,air volume (m3)
C     FLVGX,FLVGS=soil evaporation-condensation unlimited,
C        limited by vapor (m3 t-1)
C     VPGV=soil saturated vapor concentration (m3 m-3)
C     HFLVGS=soil latent heat flux (MJ t-1)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C     WFLVL,HFLVL=total soil evaporation-condensation,
C        latent heat flux (m3 t-1,MJ t-1)
C
      PSISVG=PSISM1(NUM(NY,NX),NY,NX)+PSISO(NUM(NY,NX),NY,NX)
      VPGV=2.173E-03/TKS2
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKS2))
     3*EXP(18.0*PSISVG/(8.3143*TKS2))
      IF(VOLPM(M,NUM(NY,NX),NY,NX).GT.ZEROS(N2,N1))THEN
      FLVGX=VOLV2(NUM(NY,NX),NY,NX)-VPGV*VOLPM(M,NUM(NY,NX),NY,NX)
      FLVGS=AMAX1(FLVGX,-AMAX1(0.0,VOLW2(NUM(NY,NX),NY,NX))*XNPXX)
      HFLVGS=VAP*FLVGS
      WFLVL(NUM(NY,NX),NY,NX)=WFLVL(NUM(NY,NX),NY,NX)+FLVGS
      HFLVL(NUM(NY,NX),NY,NX)=HFLVL(NUM(NY,NX),NY,NX)+HFLVGS
      ELSE
      FLVGS=0.0
      HFLVGS=0.0
      ENDIF 
C
C     FREEZE-THAW IN SOIL SURFACE FROM NET CHANGE IN SOIL 
C     SURFACE HEAT STORAGE
C
C     TFREEZ=soil freezing temperature (K)
C     333.0=latent heat of freezing (MJ m-3)
C     TKS2=soil temperature (K) 
C     PSISVG=soil matric+osmotic potential (MPa)
C     VOLW2,VOLIw=soil water,ice volume (m3)
C     HFLFGX,HFLFGS=soil freeze-thaw latent heat flux
C        unlimited,limited by water,ice (MJ t-1)
C     FLFGS=soil freeze-thaw flux (m3 t-1)
C     WFLFL,HFLFL=total litter freeze-thaw,latent heat flux 
C        (m3 t-1,MJ t-1)
C     
      TFREEZ=-9.0959E+04/(PSISVG-333.0)
      IF((TKS2.LT.TFREEZ.AND.VOLW2(NUM(NY,NX),NY,NX)
     2.GT.ZERO*VOLT(NUM(NY,NX),NY,NX))
     3.OR.(TKS2.GT.TFREEZ.AND.VOLI2(NUM(NY,NX),NY,NX) 
     4.GT.ZERO*VOLT(NUM(NY,NX),NY,NX)))THEN
      HFLFGX=VHCP1(NUM(NY,NX),NY,NX)*(TFREEZ-TKS2)*XNPR
     2/(1.0+6.2913E-03*TFREEZ)
      IF(HFLFGX.LT.0.0)THEN
      HFLFGS=AMAX1(-333.0*DENSI*VOLI2(NUM(NY,NX),NY,NX)*XNPXX,HFLFGX)
      ELSE
      HFLFGS=AMIN1(333.0*VOLW2(NUM(NY,NX),NY,NX)*XNPSRX,HFLFGX)
      ENDIF
      FLFGS=-HFLFGS/333.0
      WFLFL(NUM(NY,NX),NY,NX)=WFLFL(NUM(NY,NX),NY,NX)+FLFGS
      HFLFL(NUM(NY,NX),NY,NX)=HFLFL(NUM(NY,NX),NY,NX)+HFLFGS
      ELSE
      HFLFGS=0.0
      FLFGS=0.0
      ENDIF
      VOLW2(NUM(NY,NX),NY,NX)=VOLW2(NUM(NY,NX),NY,NX)+FLVGS+FLFGS 
      VOLV2(NUM(NY,NX),NY,NX)=VOLV2(NUM(NY,NX),NY,NX)-FLVGS 
      VOLI2(NUM(NY,NX),NY,NX)=VOLI2(NUM(NY,NX),NY,NX)-FLFGS/DENSI 
      VOLP2(NUM(NY,NX),NY,NX)=AMAX1(0.0,VOLA1(NUM(NY,NX),NY,NX)
     2-VOLW2(NUM(NY,NX),NY,NX)-VOLI2(NUM(NY,NX),NY,NX))
C     IF(I.GT.300)THEN
C     WRITE(*,7764)'WFLVGG',I,J,NFZ,M,NX,NY
C    2,WFLVL(NUM(NY,NX),NY,NX),HFLVL(NUM(NY,NX),NY,NX)
C    3,WFLFL(NUM(NY,NX),NY,NX),HFLFL(NUM(NY,NX),NY,NX) 
C    2,FLVGS,FLFGS 
C    2,VOLW2(NUM(NY,NX),NY,NX),VOLV2(NUM(NY,NX),NY,NX)
C    2,VOLPM(M,NUM(NY,NX),NY,NX),VOLP2(NUM(NY,NX),NY,NX)
C    2,TK1(NUM(NY,NX),NY,NX),VPLV,VP1
C    3,VPLV*VOLPM(M,NUM(NY,NX),NY,NX)
C    2,PSISM1(NUM(NY,NX),NY,NX)
C    2,VPLV*VOLP2(NUM(NY,NX),NY,NX),WFLVL(NUM(NY,NX),NY,NX)
C    3,VOLV1(NUM(NY,NX),NY,NX),TFLVL(NUM(NY,NX),NY,NX)
C    2,VHCP1(NUM(NY,NX),NY,NX),HFLVL(NUM(NY,NX),NY,NX)
7764  FORMAT(A8,6I4,20E14.6)
C     ENDIF
C
C     HEAT AND VAPOR FLUXES BETWEEN SOIL SURFACE AND GROUND AIR 
C
C     TKS2,PSISVG=soil temperature,matric+osmotic potential (K,MPa)
C     VPGV=soil saturation vapor concentration (m3 m-3)
C     VP1,VPQG=vapor concentration at soil surface, canopy air 
C        above ground surface (m3 m-3)
C     PAREGM=conductance for soil vapor fluxes (m h-1) 
C     EVAPG,EVAPGV,EVAPGW=evaporation from soil surface:total,
C        from vapor,from water (m3 t-1)
C     VOLV2,VOLW2=soil vapor,water content (m3)
C     VFLXG=convective heat of evaporation flux (MJ t-1)     
C
      IF(VOLPM(M,NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      VP1=AMAX1(0.0,VOLV2(NUM(NY,NX),NY,NX)/VOLPM(M,NUM(NY,NX),NY,NX))
      ELSE
      VP1=VPGV
      ENDIF
      EVAPGX=PAREGM*(VPQG(M,NY,NX)-VP1)
      EVAPGV=AMAX1(EVAPGX,
C    2-AMAX1(0.0,VOLV2(NUM(NY,NX),NY,NX))*XNPXX)
     2-AMAX1(0.0,VOLV2(NUM(NY,NX),NY,NX)))
      EVAPGW=AMAX1(EVAPGX-EVAPGV
     2,-AMAX1(0.0,VOLW2(NUM(NY,NX),NY,NX))*XNPXX)
      EVAPG(NY,NX)=EVAPGV+EVAPGW
      IF(EVAPG(NY,NX).LT.0.0)THEN
      VFLXG=EVAPG(NY,NX)*4.19*TK1(NUM(NY,NX),NY,NX)
      ELSE
      VFLXG=EVAPG(NY,NX)*4.19*TKQG(M,NY,NX)
      ENDIF
C     IF(J.EQ.14)THEN
C     WRITE(*,7764)'EVAPG',I,J,NFZ,M,NX,NY
C    3,EVAPGX,EVAPGV,EVAPGW,PAREGM,VPQG(M,NY,NX),VP1
C    2,VOLV2(NUM(NY,NX),NY,NX),VOLW2(NUM(NY,NX),NY,NX)
C    3,VOLV1(NUM(NY,NX),NY,NX),XNPXX,EVAPG(NY,NX)
C     ENDIF
      VOLW2(NUM(NY,NX),NY,NX)=VOLW2(NUM(NY,NX),NY,NX)+EVAPGW
      VOLV2(NUM(NY,NX),NY,NX)=VOLV2(NUM(NY,NX),NY,NX)+EVAPGV 
C
C     SOIL SURFACE LATENT, SENSIBLE AND STORAGE HEAT FLUXES
C
C     EFLXG,SFLXG=soil latent,sensible heat flux (MJ t-1)
C     RFLXG=net radiation (MJ t-1)
C     VFLXG=convective heat of evaporation flux (MJ t-1)     
C     EVAPGW=evaporation directly from soil water (m3 t-1)
C     FLVGS=soil evaporation-condensation (m3 t-1)
C     FLFGS=soil freeze-thaw flux (m3 t-1)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C     PARSGM=conductance for soil sensible heat flux (m h-1) 
C     TKQG,TK1=temperature of near-surface canopy air,soil (K) 
C     HFLXG=total soil storage heat flux (MJ t-1)
C     VOLW2,VOLV2,VOLI2,VOLP2=soil water,vapor,ice,air content (m3)     
C
      EFLXG=(EVAPGW+FLVGS)*VAP
      SFLXG=AMIN1(PARSGM,VHCPG2)
     2*(TKQG(M,NY,NX)-TK1(NUM(NY,NX),NY,NX))
      HFLX0=RFLXG+EFLXG+SFLXG
      HFLXG=HFLX0+VFLXG
C     IF(J.EQ.14)THEN
C     WRITE(*,1112)'EFLXG',I,J,NFZ,M,NX,NY
C    2,EVAPG(NY,NX),EVAPGX,EVAPGW,EVAPGV,FLVGX,FLVGS,VPGV,PSISVG
C    2,PSISM1(NUM(NY,NX),NY,NX),PSISO(NUM(NY,NX),NY,NX)
C    2,TK1(NUM(NY,NX),NY,NX),VPQG(M,NY,NX),VP1,TKQG(M,NY,NX),TKR2,TKS2 
C    2,DLYRR(NY,NX),RFLXR,RFLXG,EFLXG,SFLXG,VFLXG,HFLX0,HFLXG 
C    3,RARH,RAGS,RAE,RI,RZR(NY,NX),PAREGM,PARSGM,VHCPG2
C    4,FSNX(NY,NX),BAREW(NY,NX)
C    3,RFLX0,THRMXS,THRMCS,THRMDS,ALBG
C    3,VOLV2(NUM(NY,NX),NY,NX),VOLV1(NUM(NY,NX),NY,NX)
C    4,VOLV(NUM(NY,NX),NY,NX),VOLPM(M,NUM(NY,NX),NY,NX)
C    3,FSNX(NY,NX),VOLI1(NUM(NY,NX),NY,NX),RFLX0,ALBG 
C    4,RADG(NY,NX),RADXG(NY,NX),THRYG(NY,NX),THRYW(NY,NX),THS(NY,NX)
C    5,FRADG(NY,NX),VPQG(M,NY,NX),VP1,FLQM 
C    6,PARSGM,PAREGM,HWFLQM,BAREW(NY,NX),RARH,CVRDW(NY,NX),RZR(NY,NX)
C    7,THETPY(NUM(NY,NX),NY,NX),THETPY(0,NY,NX)
C    8,VHCP1(0,NY,NX),VHCPRX(NY,NX),VHCP1(NUM(NY,NX),NY,NX)
C    3,TKQG(M,NY,NX),BARE(NY,NX),THETWX(NUM(NY,NX),NY,NX)
C    5,PSISM1(NUM(NY,NX),NY,NX),THETW1,VOLW1(NUM(NY,NX),NY,NX)
C    6,DFVR,THETPX0,POROQ,THETPX(0,NY,NX),POROS(0,NY,NX) 
1112  FORMAT(A8,6I4,60E12.4)
C     ENDIF
C
C     ENERGY BALANCE AT RESIDUE SURFACE
C
C     VHCP1,VHCPRX=current,minimum litter heat capacities
C
      IF(VHCP1(0,NY,NX).GT.VHCPRXG(NY,NX))THEN
      EVAPRW=0.0
      EVAPRV=0.0
      RFLXR=0.0
      EFLXR=0.0
      VFLXR=0.0
      SFLXR=0.0
      HFLXR=0.0
      FLV1=0.0
      HWFLV1=0.0
      HFLCR1=0.0
      HWFLQ2=HWFLQM*XNPR
      HFLXG2=HFLXG*XNPR
      HFLVG2=HFLVGS*XNPR
      HFLFG2=HFLFGS*XNPR
C
C     NET RADIATION AT RESIDUE SURFACE
C
C     ALBR=litter albedo
C     ALBZ=dry surface litter albedo
C     BKVL=litter mass (Mg)
C     VOLW1,VOLI1,VOLP1=water,ice,air volume in litter (m3)
C     RADXR,THRYR=incoming shortwave,longwave radiation (MJ t-1)
C     RFLX0=net radiation (MJ t-1)
C
      ALBR=(ALBZ(NY,NX)*BKVL(0,NY,NX)+0.06*VOLW1(0,NY,NX)+0.30
     2*VOLI1(0,NY,NX))/(BKVL(0,NY,NX)+VOLW1(0,NY,NX)+VOLI1(0,NY,NX))
      RFLX0=(1.0-ALBR)*RADXR(NY,NX)+THRYR(NY,NX)
C
C     VAPOR CONDUCTIVITY BETWEEN SURFACE RESIDUE AND SOIL SURFACE
C
C     CNVR,CNV1=litter,soil vapor conductivity (m2 h-1)
C     WGSGR,WGSGL=litter,soil vapor diffusivity (m2 h-1)
C     THETPM=litter air concentration (m3 m-3)
C     POROS,POROQ=litter porosity m3 m-3), tortuosity from ‘starts.f’ 
C     CVRD=litter cover fraction
C     AVCNVR=litter-soil vapor conductance (m h-1)
C     DLYRR,DLYR=litter,soil depths (m)
C
      CNVR=WGSGR(NY,NX)*POROQ
     2*THETPM(M,0,NY,NX)**2/POROS(0,NY,NX) 
      CNV1=WGSGL(NUM(NY,NX),NY,NX)*POROQ
     2*THETPM(M,NUM(NY,NX),NY,NX)**2/POROS(NUM(NY,NX),NY,NX)
      IF(CVRD(NY,NX).GT.ZERO)THEN
      IF(CNVR.GT.ZERO.AND.CNV1.GT.ZERO)THEN
      AVCNVR=2.0*CNVR*CNV1
     2/(CNVR*DLYR(3,NUM(NY,NX),NY,NX)+CNV1*DLYRR(NY,NX)) 
      ELSE
      AVCNVR=2.0*CNVR
     2/(DLYR(3,NUM(NY,NX),NY,NX)+DLYRR(NY,NX))*CVRD(NY,NX)
      ENDIF
      ELSE
      AVCNVR=0.0
      ENDIF
C
C     THERMAL CONDUCTIVITY BETWEEN SURFACE RESIDUE AND SOIL SURFACE
C
C     THETRR=dry litter concentration (m3 m-3)
C     DTH*,RYL*,XNU*,TRB*=turbulence effects on thermal conductivity
C        of air (*A) and water (*W)
C     WTHET0,WTHET1=multiplier for air concentration 
C        on thermal conductivity
C     TCNDW*,TCNDA*=thermal conductivity of water,air (m MJ h-1 K-1)
C     TCNDR,TCND1=litter,soil thermal conductivity (m MJ h-1 K-1)
C     ATCNDR=litter-soil thermal conductance (MJ h-1 K-1)
C     
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
C     VHCPR2,VHCPRXG=current,minimum litter heat capacities (MJ K-1)
C     VHCPG2,VHCPRX=current,minimum soil heat capacities(MJ K-1)
C     TKR2,TKS2=litter,soil surface interim temperature (k)
C     HFLXF=heat released by combustion in previous time step (MJ t-1)
C
      DO 5000 NN=1,NPR
      IF(VHCPR2.GT.VHCPRXG(NY,NX).AND.VHCPG2.GT.VHCPRX(NY,NX))THEN
      TKR2=TKR2+HFLXF(0,NY,NX)*XNPR/VHCPR2
      TKS2=TKS2+HFLXF(NUM(NY,NX),NY,NX)*XNPR/VHCPG2
C
C     NET RADIATION AT RESIDUE SURFACE
C
C     THRMXR2=longwave radiation emitted by litter (MJ t-1)
C     RFLXR2=litter net radiation (MJ t-1)
C     THRMCR,THRMDR=net LW exchange between litter and canopy,
C        standing dead surfaces (MJ t-1)
C     TKC,TKD=surface temperatures of canopy,
C        standing dead surfaces (K) 
C     FRADP,FRADQ=fraction of incoming radiation 
C        received by each PFT canopy, standing dead from ‘hour1.f’ 
C     THETWR=litter water content (m3 m-3)
C     VOLWRX=litter water retention capacity (m3)
C     PSISM1=litter matric water potential (MPa)
C     PSL,FCL,WPL=log POROS0,FC,WP
C     FCD,PSD=FCL-WPL,PSL-FCL
C     SRP=parameter for deviation from linear log-log water retention 
C
      THRMXR2=EMMGR(NY,NX)*TKR2**4
      RFLXR2=RFLX0-THRMXR2
      DO 915 NZ=1,NP(NY,NX)
      THRMCR2=EMMCR(NY,NX)*(TKC(NZ,NY,NX)**4-TKR2**4)
     2*FRADP(NZ,NY,NX)
      THRMDR2=EMMCR(NY,NX)*(TKD(NZ,NY,NX)**4-TKR2**4)
     2*FRADQ(NZ,NY,NX)
      RFLXR2=RFLXR2+THRMCR2+THRMDR2
915   CONTINUE 
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX)
     2.AND.VOLW2(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW2(0,NY,NX))/VOLR(NY,NX)
      IF(THETWR.LT.FC(0,NY,NX))THEN
      PSISM1(0,NY,NX)=AMAX1(PSISX,-EXP(PSIMX(NY,NX)
     2+((FCL(0,NY,NX)-LOG(THETWR))
     3/FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETWR.LT.POROS0(NY,NX))THEN 
      PSISM1(0,NY,NX)=-EXP(PSIMS(NY,NX)
     2+((AMAX1(0.0,(PSL(0,NY,NX)-LOG(THETWR)))
     3/PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
C     WRITE(*,2111)'PSISM0',I,J,NFZ,M,NX,NY 
C    4,PSISM(0,NY,NX),PSISM1(0,NY,NX),PSIMS(NY,NX)
C    2,PSL(0,NY,NX),LOG(THETWR),PSL(0,NY,NX)-LOG(THETWR)
C    3,PSD(0,NY,NX),SRP(0,NY,NX),PSISD(NY,NX)
C    4,THETWR,FC(0,NY,NX),POROS0(NY,NX),PSISE(0,NY,NX)
2111  FORMAT(A8,6I4,40E16.8)
      ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
      ELSE
      THETWR=POROS0(NY,NX)
      PSISM1(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
C
C     EVAPORATION-CONDENSATION IN SURFACE LITTER
C
C     VOLW2,VOLV2,VOLP2=litter water,vapor,air content (m3)
C     VPRV=litter saturated vapor concentration (m3 m-3)
C     FLVRX,FLVR2=litter evaporation-condensation unlimited,
C        limited by vapor (m3 t-1)
C     HFLVR2=litter latent heat flux (MJ t-1)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C     WFLVR,HFLVR=litter evaporation-condensation,
C        latent heat flux (m3 t-1,MJ t-1)
C
      PSISVR=PSISM1(0,NY,NX)+PSISO(0,NY,NX)
      VPRV=2.173E-03/TKR2 
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKR2))
     3*EXP(18.0*PSISVR/(8.3143*TKR2))
      FLVRX=VOLV2(0,NY,NX)-VPRV*VOLP2(0,NY,NX)
      FLVR2=AMAX1(FLVRX,-AMAX1(0.0,VOLW2(0,NY,NX))*XNPRX)
      HFLVR2=VAP*FLVR2
      WFLVR(NY,NX)=WFLVR(NY,NX)+FLVR2
      HFLVR(NY,NX)=HFLVR(NY,NX)+HFLVR2
C
C     FREEZE-THAW IN SURFACE LITTER FROM NET CHANGE IN RESIDUE 
C     SURFACE HEAT STORAGE
C
C     TFREEZ=litter freezing temperature (K)
C     PSISVR=litter natric+osmotic potential (MPa)
C     333.0=latent heat of freezing (MJ m-3)
C     TKR2=litter temperature (K)
C     VOLW2,VOLI2=litter water,ice volume (m3)
C     HFLFRX,HFLFR2=litter freeze-thaw latent heat flux
C        unlimited,limited by water,ice (MJ t-1)
C     FLFR2=litter freeze-thaw flux (m3 t-1)
C     HFLFR,WFLFR=total litter freeze-thaw,latent heat flux
C        (m3 t-1,MJ t-1)
C     
      TFREEZ=-9.0959E+04/(PSISVR-333.0)
      IF((TKR2.LT.TFREEZ
     2.AND.VOLW2(0,NY,NX).GT.ZERO*VOLT(0,NY,NX))
     3.OR.(TKR2.GT.TFREEZ 
     4.AND.VOLI2(0,NY,NX).GT.ZERO*VOLT(0,NY,NX)))THEN
      HFLFRX=VHCPR2*(TFREEZ-TKR2)*XNPR
     2/(1.0+6.2913E-03*TFREEZ)
      IF(HFLFRX.LT.0.0)THEN
      HFLFR2=AMAX1(-333.0*DENSI*VOLI2(0,NY,NX)*XNPRX,HFLFRX)
      ELSE
      HFLFR2=AMIN1(333.0*VOLW2(0,NY,NX)*XNPRX,HFLFRX)
      ENDIF
      FLFR2=-HFLFR2/333.0
      HFLFR(NY,NX)=HFLFR(NY,NX)+HFLFR2
      WFLFR(NY,NX)=WFLFR(NY,NX)+FLFR2
      ELSE
      HFLFR2=0.0
      FLFR2=0.0
      ENDIF
      VOLW2(0,NY,NX)=VOLW2(0,NY,NX)+FLVR2+FLFR2 
      VOLI2(0,NY,NX)=VOLI2(0,NY,NX)-FLFR2/DENSI 
      VOLV2(0,NY,NX)=VOLV2(0,NY,NX)-FLVR2
      VOLP2(0,NY,NX)=AMAX1(0.0,VOLA2-VOLW2(0,NY,NX)
     2-VOLI2(0,NY,NX))
C     IF(NX.EQ.1.AND.NY.EQ.1.AND.HFLF1.NE.0.0)THEN
C     WRITE(*,5353)'WFLVRG',I,J,NFZ,M,NN,NX,NY
C    2,WFLVR(NY,NX),HFLVR(NY,NX),WFLFR(NY,NX),HFLFR(NY,NX) 
C    3,FLVR2,FLVRX,VOLV2(0,NY,NX)
C    4,VPRV,VOLP2(0,NY,NX),VPRV*VOLP2(0,NY,NX)
C    5,TKR2,PSISVR,EVAPR2V 
5353  FORMAT(A8,7I4,30E12.4)
C     ENDIF
C
C     VAPOR FLUX AT LITTER SURFACE
C
C     VPR,VP1,VPQG=vapor concentration in litter,soil,canopy air 
C        (m3 m-3)
C     TKR2=litter interim temperature (K)
C     PSISVR=litter matric+osmotic potential (MPa)
C     VPRV=litter saturated vapor concentration (m3 m-3)
C     VOLW2,VOLV2,VOLP2=litter water,vapor,air content (m3)
C     EVAPR2X=unlimited vapor flux at litter surface (m3 t-1)
C     PARERM=conductance for latent heat flux (m t-1) 
C     EVAPR2V,EVAPR2W=vapor flux from vapor,water at
C        litter surface (m3 t-1)   
C     VFLXR2=convective heat of vapor flux (MJ t-1)
C     FLYM2=precipitation flux to litter (m3 t-1)     
C
      IF(VOLP2(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      VPR=AMAX1(0.0,VOLV2(0,NY,NX)/VOLP2(0,NY,NX))
      ELSE
      VPR=VPRV
      ENDIF
      EVAPR2X=PARERM*(VPQG(M,NY,NX)-VPR)
      EVAPR2V=AMAX1(EVAPR2X
C    2,-AMAX1(0.0,VOLV2(0,NY,NX))*XNPXX)
     2,-AMAX1(0.0,VOLV2(0,NY,NX)))
      EVAPR2W=AMAX1(EVAPR2X-EVAPR2V
     2,-AMAX1(0.0,VOLW2(0,NY,NX))*XNPXX)
      EVAPR2=EVAPR2V+EVAPR2W
      IF(EVAPR2X.LT.0.0)THEN
      VFLXR2=EVAPR2*4.19*TKR2
      ELSE
      VFLXR2=EVAPR2*4.19*TKQG(M,NY,NX)
      ENDIF
C     IF(J.EQ.14)THEN
C     WRITE(*,5353)'EVAPR2',I,J,NFZ,M,NN,NX,NY
C    2,EVAPR2X,EVAPR2V,EVAPR2W,PARERM,VPQG(M,NY,NX),VPR 
C    3,VOLV2(0,NY,NX),VOLW2(0,NY,NX),FLVR2,PSISVR
C    4,XNPXX,EVAPR2,EVAPR(NY,NX)
C     ENDIF
      VOLW2(0,NY,NX)=VOLW2(0,NY,NX)+EVAPR2W+FLYM2 
      VOLV2(0,NY,NX)=VOLV2(0,NY,NX)+EVAPR2V 
C
C     HEAT AND VAPOR FLUXES BETWEEN LITTER AND GROUND AIR 
C
C     SFLXR2,EFLXR2,RFLXR2=sensible,latent heat fluxes, 
C        net radiation at litter surface (MJ t-1)
C     EVAPR2W=vapor flux from water at litter surface (m3 t-1)   
C     FLVR2=litter evaporation-condensation (m3 t-1)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C     HFLX02,HFLXR2=storage,total litter heat flux (MJ t-1) 
C     PARSRM=conductance for litter sensible heat flux (m h-1)
C     TKQG,TKR2=temperature of near-surface canopy air,litter 
C     VFLXR2=convective heat of vapor flux (MJ t-1)
C
      EFLXR2=(EVAPR2W+FLVR2)*VAP
      SFLXR2=PARSRM*(TKQG(M,NY,NX)-TKR2)
      HFLX02=RFLXR2+EFLXR2+SFLXR2
      HFLXR2=HFLX02+VFLXR2
C
C     VAPOR AND HEAT FLUXES BETWEEN RESIDUE AND SOIL SURFACE 
C
C     VOLW2,VOLV2,VOLI2,VOLP2=litter(0,),soil(NUM,)
C        water,vapor,ice,air content (m3)
C     VPR,VP1=litter,soil vapor concentration (m3 m-3)
C     FLVC=vapor unlimited vapor flux (m3 t-1)
C     AVCNVR=litter-soil vapor conductance (m h-1)
C     AREA=area of grid cell (m2)
C     FSNX=snow-free surface fraction
C     CVRD=fraction of litter cover
C     FLV2=vapor unlimited litter-soil vapor flux (m3 t-1)
C     HWFLV2=convective heat of litter-soil vapor flux (MJ t-1)
C     TKY=equilibrium litter-soil temperature (K)     
C     VHCPR2,VHCPG2=litter,soil heat capacity (MJ K-1)
C     TKR2,TKS2=litter,soil temperature (K)
C     HFLWC,HFLWX=litter-soil heat flux unlimited,limited by heat 
C        (MJ t-1)
C     ATCNDR=litter-soil thermal conductance (MJ h-1 K-1)
C     HFLCR2=litter-soil conductive heat flux (MJ t-1) 
C     XNPZX=time step for litter flux from ‘wthr.f’(h t-1)
C
      IF(VOLP2(0,NY,NX).GT.ZEROS(NY,NX)
     2.AND.VOLP2(NUM(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      VP1=AMAX1(0.0,VOLV2(NUM(NY,NX),NY,NX)/VOLP2(NUM(NY,NX),NY,NX))
      FLVC=AVCNVR*(VPR-VP1)*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)
     2*CVRD(NY,NX)*XNPZX
      IF(FLVC.GE.0.0)THEN
      FLV2=AMAX1(0.0,AMIN1(FLVC,VOLV1(0,NY,NX)*XNPXX))
      HWFLV2=4.19*TKR2*FLV2
      ELSE
      FLV2=AMIN1(0.0,AMAX1(FLVC,-VOLV2(NUM(NY,NX),NY,NX)*XNPXX))
      HWFLV2=4.19*TKS2*FLV2
      ENDIF
      ELSE
      FLV2=0.0
      HWFLV2=0.0 
      ENDIF
      HFLWC=ATCNDR*(TKR2-TKS2)*AREA(3,NUM(NY,NX),NY,NX)*FSNX(NY,NX)
     2*CVRD(NY,NX)*XNPZX 
      TKY=(TKR2*VHCPR2+TKS2*VHCPG2)/(VHCPR2+VHCPG2)
      HFLWX=(TKR2-TKY)*VHCPR2 
      IF(HFLWC.GE.0.0)THEN
      HFLCR2=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLCR2=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
      VOLV2(0,NY,NX)=VOLV2(0,NY,NX)-FLV2
      VOLV2(NUM(NY,NX),NY,NX)=VOLV2(NUM(NY,NX),NY,NX)+FLV2 
C
C     AGGREGATE WATER AND ENERGY FLUXES FROM TIME STEP FOR LITTER
C     CALCULATIONS TO THAT FOR SOIL PROFILE
C
C     EVAPR=litter evaporation (m3 t-1)
C     SFLXR,EFLXR,RFLXR=sensible,latent heat fluxes, net radiation 
C        (MJ t-1)
C     VFLXR=convective heat flux from evaporation (MJ t-1)
C     HFLXR=storage heat flux (MJ t-1)     
C     HWFLV1=convective heat of litter-soil vapor flux (MJ m-1)
C     HFLCR1=litter-soil conductive heat flux (MJ t-1)
C
      EVAPR(NY,NX)=EVAPR(NY,NX)+EVAPR2
      EVAPRW=EVAPRW+EVAPR2W
      EVAPRV=EVAPRV+EVAPR2V
      RFLXR=RFLXR+RFLXR2
      EFLXR=EFLXR+EFLXR2
      VFLXR=VFLXR+VFLXR2
      SFLXR=SFLXR+SFLXR2
      HFLXR=HFLXR+HFLXR2
      FLV1=FLV1+FLV2
      HWFLV1=HWFLV1+HWFLV2
      HFLCR1=HFLCR1+HFLCR2 
      ELSE
      EVAPR2=0.0
      EVAPR2V=0.0
      EVAPR2W=0.0
      FLVR2=0.0
      HFLVR2=0.0
      RFLXR2=0.0
      EFLXR2=0.0
      VFLXR2=0.0
      SFLXR2=0.0
      HFLXR2=0.0
      FLV2=0.0
      HWFLV2=0.0
      HFLCR2=0.0
      ENDIF
C
C     UPDATE LITTER WATER,VAPOR,ICE,AIR,TEMPERATURE
C
C     VOLW2,VOLV2,VOLI2,VOLP2=litter(0,),soil(NUM,)
C        water,vapor,ice,air content (m3)
C     FLVR2=litter evaporation-condensation (m3 t-1)
C     FLFR2=litter freeze-thaw (m3 t-1)
C     DENSI=ice density (Mg m-3)
C     FLV2=litter-soil vapor flux (m3 t-1)
C     TKR2,TKS2=litter,soil surface temperature (K)
C     VHCPR2,VHCPG2=litter,soil heat capacity (MJ K-1)
C     VHCM=soil solid fraction heat capacity from ‘hour1.f’ (MJ K-1)
C     ORGC,ORGCC=litter organic C, charcoal
C     FSNX=snow-free surface fraction
C     HFLXR2=total litter heat flux (MJ t-1) 
C     HWFLM2=convective heat flux from precipitation to litter (MJ t-1)    
C     HWFLV2=convective heat of litter-soil vapor flux (MJ t-1)
C     HFLCR2=litter-soil conductive heat flux (MJ t-1) 
C     HFLVR2=litter latent heat flux (MJ t-1)
C     HFLFR2=litter freeze-thaw latent heat flux (MJ t-1)
C     HWFLQ2=total convective heat flux to soil micropores, macropores
C        (MJ t-1)
C     HFLXR2=total soil heat flux (MJ t-1) 
C     HFLVG2=soil latent heat flux (MJ t-1)
C     HFLFG2=soil freeze-thaw latent heat flux (MJ t-1)
C
      ENGYR=VHCPR2*TKR2
      VHCPR2=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))*FSNX(NY,NX)
     2+4.19*(VOLW2(0,NY,NX)+VOLV2(0,NY,NX))
     2+1.9274*VOLI2(0,NY,NX)
      IF(VHCPR2.GT.VHCPRXG(NY,NX))THEN
      TKR2=(ENGYR+HFLXR2+HWFLM2-HWFLV2-HFLCR2+HFLVR2+HFLFR2)/VHCPR2
      ELSE
      TKR2=TK1(0,NY,NX)
      ENDIF
      ENGY2=VHCPG2*TKS2
      VHCPG2=VHCM(NUM(NY,NX),NY,NX)
     2+4.19*(VOLW2(NUM(NY,NX),NY,NX)+VOLV2(NUM(NY,NX),NY,NX)
     2+VOLWH1(NUM(NY,NX),NY,NX))
     2+1.9274*(VOLI2(NUM(NY,NX),NY,NX)+VOLIH1(NUM(NY,NX),NY,NX))
      IF(VHCPG2.GT.VHCPRX(NY,NX)
     2.AND.VOLX(NUM(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      TKS2=(ENGY2+HFLXG2+HWFLQ2+HFLVG2+HFLFG2+HWFLV2+HFLCR2)/VHCPG2
      ELSE
      TKS2=TK1(NUM(NY,NX),NY,NX)
      ENDIF
C     IF(I.GT.300)THEN
C     WRITE(*,1111)'VOLW2G',I,J,NFZ,M,NN,NX,NY
C    2,VOLW2(0,NY,NX),FLVR2,FLFR2,EVAPR2W,FLYM2
C    3,VOLW1(0,NY,NX),FSNX(NY,NX)
C     WRITE(*,1111)'TKR2G',I,J,NFZ,M,NN,NX,NY
C    2,TKR2,ENGYR,HFLXR2,HWFLM2,HWFLV2,HFLCR2,HFLVR2,HFLFR2 
C    3,VHCPR2,FSNX(NY,NX),VOLW2(0,NY,NX),VOLV2(0,NY,NX)
C    4,VOLI2(0,NY,NX),TK1(0,NY,NX),TKS(0,NY,NX),VHCPRXG(NY,NX)
C     WRITE(*,1111)'TKS2G',I,J,NFZ,M,NN,NX,NY
C    2,TKS2,ENGY2,HWFLQ2,HFLXG2,HFLVG2,HFLFG2,HWFLV2,HFLCR2,VHCPG2
C     WRITE(*,1111)'EFLXR2',I,J,NFZ,M,NN,NX,NY
C    2,EVAPR2,EVAPR2X,EVAPR2W,EVAPR2V,FLVRX,FLVR2,FLVC,FLV2,FLV1
C    3,EVAPR(NY,NX),EVAPRW,EVAPRV,WFLVR(NY,NX)
C    4,VPRV,VPR,VP1,VPQG(M,NY,NX),TKR2,TKS2,PSISVR,FSNX(NY,NX)
C    4,PARERM,VPQG(M,NY,NX),TKQG(M,NY,NX)
C    5,PSISM1(0,NY,NX),TK1(0,NY,NX),TK1(NUM(NY,NX),NY,NX)
C    6,HFLXF(0,NY,NX),ENGYR,HWFLM2,HWFLV2,HFLCR2,VHCPR2
C    3,VOLV2(NUM(NY,NX),NY,NX),VOLV1(NUM(NY,NX),NY,NX)
C    4,VOLV(NUM(NY,NX),NY,NX),VOLPM(M,NUM(NY,NX),NY,NX) 
C    2,RFLXR2,EFLXR2,SFLXR2,VFLXR2,HFLXR2,SFLXR 
C    3,VOLV2(0,NY,NX),VOLP2(NUM(NY,NX),NY,NX),VPRV*VOLP2(0,NY,NX)
C    3,HFLXR2,HWFLM2,HWFLV2,HFLCR2,HFLXR2+HWFLM2-HWFLV2-HFLCR2
C    3,RFLX0,THRMXR2,THRMCR2,EMMGR(NY,NX),RADXR(NY,NX),THRYR(NY,NX)
C    4,TKCT(NY,NX),CVRD(NY,NX),CVRDW(NY,NX),VHCPR2,PARERM,PARSRM 
C    3,HFLX02,HFLXR2,HWFLM2,HWFLV2,HFLCR2,HFLWX,HFLWC 
1111  FORMAT(A8,7I4,100E12.4)
C     ENDIF
5000  CONTINUE
C
C     IF NO SURFACE LITTER
C
      ELSE
      TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
      TKSM(1,0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
      EVAPR(NY,NX)=0.0
      EVAPRW=0.0
      EVAPRV=0.0
      RFLXR=0.0
      EFLXR=0.0
      VFLXR=0.0
      SFLXR=0.0
      HFLXR=0.0
      FLV1=0.0
      HWFLV1=0.0
      HFLCR1=0.0
      ENDIF
C
C     GATHER LITTER AND SOIL SURFACE WATER, VAPOR AND HEAT FLUXES WITH
C     ATMOSPHERE INTO FLUX ARRAYS FOR LATER UPDATES TO STATE VARIABLES
C
C     FLWLG,FLWHLG=net water flux from canopy air to soil
C        micropores,macropores (m3 t-1)
C     FLVLG=vapor flux from canopy air to soil (m3 t-1)
C     HFLWLG=net convective heat flux from canopy air to soil (MJ t-1) 
C     FLWRLG,FLVRLG=net water,vapor flux from canopy air to litter 
C        (m3 t-1)
C     HFLWRLG=net convective heat flux from canopy air to litter 
C        (MJ t-1)
C     FLWVLS=water flux within soil accounting for wetting front 
C        (m3 t-1)
C     FLQM,FLHM=total water flux to soil micropores, macropores
C        (m3 t-1)
C     EVAPGW,EVAPGV=evaporation from soil water,vapor (m3 t-1)
C     HWFLQM=total convective heat flux to soil micropores, macropores
C        (MJ t-1)
C     HFLXG=total soil storage heat flux (MJ t-1)
C     FLV1,HWFLV1=vapor,convective heat of litter-soil vapor flux 
C        (MJ m-1)
C     HFLCR1=litter-soil conductive heat flux (MJ t-1)
C     EVAPRW,EVAPRV=evaporation from litter water,vapor (m3 t-1)
C     FLYM,HWFLYM=total precipitation flux, convective heat flux to
C        litter (m3 t-1,MJ t-1))
C
      FLWLG=FLQM+EVAPGW 
      FLVLG=FLV1+EVAPGV 
      FLWLXG=FLQM+EVAPGW 
      FLWHLG=FLHM
      HFLWLG=HWFLQM+HFLXG+HWFLV1+HFLCR1
      FLWRLG=FLYM+EVAPRW
      FLVRLG=-FLV1+EVAPRV
      HFLWRLG=HFLXR+HWFLYM-HWFLV1-HFLCR1
      FLWVLS=(VOLW1(NUM(NY,NX),NY,NX)-VOLWX1(NUM(NY,NX),NY,NX))*XNPHX
C     IF(I.GT.270)THEN
C     WRITE(*,7749)'FLWLG',I,J,NFZ,M,NX,NY
C    2,FLWLG,FLQM 
C    3,FLWRLG,FLYM
C    4,FLVLG,FLV1,EVAPG(NY,NX),EVAPGW,EVAPGV 
C    5,FLVRLG,FLV1,EVAPR(NY,NX),EVAPRW,EVAPRV
C    4,HFLWLG,HWFLQM,HFLXG,HWFLV1,HFLCR1
C    5,HFLWRLG,HWFLYM,HFLXR,HWFLV1,HFLCR1
7749  FORMAT(A8,6I4,60E12.4) 
C     ENDIF
C
C     GENERATE NEW SNOWPACK
C
C     VHCPW,VHCPWX=current,minimum snowpack heat capacity (MJ K-1)
C     XFLWS,XFLWW,XFLWI=net snow,water,ice transfer (m3 t-1)
C     FLQ0S,FLQ0W,FLQ0I=snow,water,ice input to snowpack (m3 t-1)
C     XHFLWW=net convective heat flux from snow,water,ice transfer 
C        (MJ t-1)
C     HWFLQ0=convective heat flux from snow,water,ice to snowpack
C        (MJ t-1)
C
      IF(VHCPW(1,NY,NX).LE.VHCPWX(NY,NX)
     2.AND.FLQ0S(NY,NX).GT.ZEROS(NY,NX))THEN
      XFLWS(1,NY,NX)=XFLWS(1,NY,NX)+FLQ0S(NY,NX)
      XFLWW(1,NY,NX)=XFLWW(1,NY,NX)+FLQ0W(NY,NX)
      XFLWI(1,NY,NX)=XFLWI(1,NY,NX)+FLQ0I(NY,NX)
      XHFLWW(1,NY,NX)=XHFLWW(1,NY,NX)+HWFLQ0(NY,NX)
C     WRITE(*,4422)'INIT',I,J,NFZ,M,FLQ0S(NY,NX),FLQ0W(NY,NX)
C    3,FLQ0I(NY,NX),HWFLQ0(NY,NX),XFLWS(1,NY,NX),XFLWW(1,NY,NX)
C    2,XFLWI(1,NY,NX),XHFLWW(1,NY,NX),HFLWL(3,NUM(NY,NX),NY,NX)
C    3,HFLW(3,NUM(NY,NX),NY,NX),FSNX(NY,NX),VHCP1(NUM(NY,NX),NY,NX)
C    4*TK1(NUM(NY,NX),NY,NX),HFLWRL(NY,NX),HFLWR(NY,NX)
C    5,VHCP1(0,NY,NX)*TK1(0,NY,NX),HEATH(NY,NX),RFLXG,RFLXR,RFLXW 
C    2,SFLXG,SFLXR,SFLXW,EFLXG,EFLXR,EFLXW,VFLXG,VFLXR,VFLXW
      ENDIF 
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
      FLVLG=0.0 
      FLWLXG=0.0 
      FLWHLG=0.0 
      HFLWLG=0.0
      FLWRLG=0.0 
      FLVRLG=0.0 
      HFLWRLG=0.0
      FLWVLS=0.0
      EVAPRW=0.0
      EVAPRV=0.0
      ENDIF
C
C     AGGREGATE LITTER AND SOIL SURFACE FLUXES BENEATH SNOW 
C     AND ATMOSPHERE
C
C     FLWL,FLWLX=total water flux into soil micropores (m3 t-1)
C     FLVL=total vapor flux into soil (m3 t-1) 
C     FLWHL=total water flux into soil macropores (m3 t-1) 
C     HFLWL=total heat flux into soil (MJ t-1)
C     FLWRL,FLWLX=total water flux into litter (m3 t-1) 
C     FLVRL=total vapor flux into litter (m3 t-1) 
C     HFLWRL=total heat flux into litter (MJ t-1)
C     F*W,F*G=water or vapor fluxes under snowpack,
C        canopy air (m3 t-1) 
C     H*W,H*G=heat fluxes under snowpack,
C        canopy air (MJ t-1)
C
      FLWL(3,NUM(NY,NX),NY,NX)=FLWLW+FLWLG 
      FLVL(3,NUM(NY,NX),NY,NX)=FLWLV+FLVLG 
      FLWLX(3,NUM(NY,NX),NY,NX)=FLWLXW+FLWLXG 
      FLWHL(3,NUM(NY,NX),NY,NX)=FLWHLW+FLWHLG
      HFLWL(3,NUM(NY,NX),NY,NX)=HFLWLW+HFLWLG 
      FLWRL(NY,NX)=FLWRLW+FLWRLG 
      FLVRL(NY,NX)=FLVRLW+FLVRLG 
      HFLWRL(NY,NX)=HFLWRLW+HFLWRLG 
C     IF(I.GT.270)THEN
C     WRITE(*,7756)'FLWT',I,J,NFZ,M,NX,NY,NUM(NY,NX)
C    2,FLWL(3,NUM(NY,NX),NY,NX),FLWLW,FLWLG
C    4,FLWHL(3,NUM(NY,NX),NY,NX),FLWHLW,FLWHLG 
C    5,HFLWL(3,NUM(NY,NX),NY,NX),HFLWLW,HFLWLG 
C    6,FLWRL(NY,NX),FLWRLW,FLWRLG
C    6,FLVL(3,NUM(NY,NX),NY,NX),FLWLV,FLVLG 
C    6,FLVRL(NY,NX),FLVRLW,FLVRLG
C    7,HFLWRL(NY,NX),HFLWRLW,HFLWRLG
C    8,HWFLQM,HFLXG,HWFLV1,HFLCR1
C    8,HFLXR,HWFLYM,HWFLV1,HFLCR1
C    9,HFLWQG,HFLVS1,HFLWS1,HFLVR1,HFLWR1
C    1,HFLWQR,HFLVSR,HFLWSR,HFLVR1,HFLWR1
C    8,FLQ0S(NY,NX),FLQ0W(NY,NX),FLQ0I(NY,NX),FLQM,FLHM,FLYM
7756  FORMAT(A8,7I4,60E12.4) 
C     ENDIF
C
C     CAPILLARY EXCHANGE OF WATER BETWEEN SOIL SURFACE AND RESIDUE
C
C     BKDS=soil bulk density (0=pond) (Mg m-3)
C     VOLW10,VOLP10=litter water,air volume (m3)
C     VOLW1N,VOLP1N=soil surface water,air volume (m3)
C     VOLP1Z=excess water+ice in micropores during freezing (if –ve) 
C        (m3)
C     VOLR=litter volume (m3)
C     THETWR,THETW1=litter,soil surface water concentration (m3 m-3)
C     PSISE,PSISM10,PSISM1N=saturation,current litter,soil surface
C        water potential (MPa)
C     POROS,FC,WP=saturation,field capacity,wilting point 
C        water concentration (m3 m-3)
C     CNDR,CND1=current litter,soil surface hydraulic conductivity
C        (m2 h MPa-1) 
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
C     AVCNDR=hydraulic conductuctance between litter and soil surface
C        (m h MPa-1) 
C     DLYRR,DLYR=litter,soil thicknesses (m)
C     PSIST,PSISH,PSISO=total,gravimetric,osmotic potentials (MPa)
C     ORFLN=reflection coefficient for osmotic potential-driven water
C        flux 
C     THETWX,POROS=soil water content,porosity (m3 m-3)
C     FLQX=litter-soil water flux unlimited by water content (m3 t-1)
C     FLQZ=FLQX + saturated litter-soil water flux (m3 t-1)
C     THETS=water concentration at air entry potential from ‘hour1.f’
C        (m3 m-3) 
C     FLQR,FLQ2=soil water flux limited by water content (m3 t-1)
C     XNPRX=time step of litter flux from ‘wthr.f’(h t-1)
C     CVRD=fraction of litter cover
C     HFLQR=convective heat from litter-soil water flux (MJ t-1)
C 
      IF(BKDS(NUM(NY,NX),NY,NX).GT.ZERO)THEN
      VOLW10=VOLW2(0,NY,NX)
      VOLW1N=VOLW2(NUM(NY,NX),NY,NX) 
      VOLP10=AMAX1(0.0,VOLWRX(NY,NX)-VOLW10)
      VOLP1ZN=VOLP1Z(NUM(NY,NX),NY,NX)
      VOLP1N=AMAX1(0.0,VOLP1(NUM(NY,NX),NY,NX)
     2-FLWL(3,NUM(NY,NX),NY,NX)) 
      DO 6000 NN=1,NPR
      IF(VOLR(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW10)/VOLR(NY,NX)
      ELSE
      THETWR=POROS0(NY,NX)
      ENDIF
      IF(THETWR.LT.FC(0,NY,NX))THEN
      PSISM10=AMAX1(PSISX,-EXP(PSIMX(NY,NX)
     2+((FCL(0,NY,NX)-LOG(THETWR))
     3/FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETWR.LT.POROS0(NY,NX))THEN 
      PSISM10=-EXP(PSIMS(NY,NX)
     2+((AMAX1(0.0,(PSL(0,NY,NX)-LOG(THETWR)))
     3/PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
      ELSE
      PSISM10=PSISE(0,NY,NX)
      ENDIF
      THETW1=AMAX1(THETZ(NUM(NY,NX),NY,NX)
     2,AMIN1(POROS(NUM(NY,NX),NY,NX),VOLW1N/VOLY(NUM(NY,NX),NY,NX)))
      IF(THETW1.LT.FC(NUM(NY,NX),NY,NX))THEN
      PSISM1N=AMAX1(PSISX,-EXP(PSIMX(NY,NX)
     2+((FCL(NUM(NY,NX),NY,NX)-LOG(THETW1))
     3/FCD(NUM(NY,NX),NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETW1.LT.POROS(NUM(NY,NX),NY,NX)-DTHETW)THEN 
      PSISM1N=-EXP(PSIMS(NY,NX)
     2+((AMAX1(0.0,(PSL(NUM(NY,NX),NY,NX)-LOG(THETW1)))
     3/PSD(NUM(NY,NX),NY,NX))**SRP(NUM(NY,NX),NY,NX)*PSISD(NY,NX)))
      ELSE
      PSISM1N=PSISE(NUM(NY,NX),NY,NX)
      ENDIF
      K0=MAX(1,MIN(100,INT(100.0*(AMAX1(0.0,POROS0(NY,NX)
     2-THETWR))/POROS0(NY,NX))+1))
      K1=MAX(1,MIN(100,INT(100.0*(AMAX1(0.0,POROS(NUM(NY,NX),NY,NX)
     2-THETW1))/POROS(NUM(NY,NX),NY,NX))+1))
      CNDR=HCND(3,K0,0,NY,NX) 
      CND1=HCND(3,K1,NUM(NY,NX),NY,NX)*FKSAT 
      AVCNDR=2.0*CNDR*CND1/(CNDR*DLYR(3,NUM(NY,NX),NY,NX)
     2+CND1*DLYRR(NY,NX)) 
      PSIST0=PSISM10+PSISH(0,NY,NX)
     2+ORFLN*PSISO(0,NY,NX)
      PSIST1=PSISM1N+PSISH(NUM(NY,NX),NY,NX)
     2+ORFLN*PSISO(NUM(NY,NX),NY,NX)
      FLQX=AVCNDR*(PSIST0-PSIST1)
     2*AREA(3,NUM(NY,NX),NY,NX)*CVRDW(NY,NX)*XNPZX
      IF(FLQX.GE.0.0)THEN
      IF(THETWR.GT.THETS(0,NY,NX))THEN
      FLQZ=FLQX+AMIN1((THETWR-THETS(0,NY,NX))
     2*VOLR(NY,NX),AMAX1(0.0,(THETS(NUM(NY,NX),NY,NX)-THETW1)
     3*VOLY(NUM(NY,NX),NY,NX)))*XNPRX
      ELSE
      FLQZ=FLQX
      ENDIF
      FLQR=AMAX1(0.0,AMIN1(FLQZ,VOLW10*XNPRX,VOLP1N*XNPRX))
      FLQ2=AMAX1(0.0,AMIN1(FLQX,VOLW10*XNPRX,VOLP1N*XNPRX))
C     WRITE(*,4322)'FLQ1',I,J,NFZ,M,NN,NX,NY,NUM(NY,NX),K0,K1
C    2,FLWRM(M,NY,NX),FLQR,FLQ2,FLQX,FLQZ,THETQ,THETWR
C    3,VOLR(NY,NX),VOLP1N,PSIST0,PSIST1
      ELSE
      FLQZ=FLQX
      FLQR=AMIN1(0.0,AMAX1(FLQZ,-VOLW1N*XNPRX,-VOLP10*XNPRX))
      FLQ2=AMIN1(0.0,AMAX1(FLQX,-VOLW1N*XNPRX,-VOLP10*XNPRX))
C     WRITE(*,4322)'FLQ2',I,J,NFZ,M,NN,NX,NY,NUM(NY,NX),K0,K1
C    2,FLWRM(M,NY,NX),FLQR,FLQ2,FLQX,FLQZ,THETQ,THETW1
C    3,VOLY(NUM(NY,NX),NY,NX),-VOLP10
C    4,PSIST0,PSIST1 
      ENDIF
C
C     ACCOUNT FOR EXCESS WATER+ICE VOLUME DURING FREEZING
C
C     VOLWP1ZN=excess water+ice (-ve) (m3)
C     FLQR,FLQ2=soil water flux limited by water content (m3 t-1)
C     HFLQR=convective heat from litter-soil water flux (MJ t-1)
C
      IF(VOLP1ZN.LT.0.0)THEN
      FLQR=FLQR+AMIN1(0.0,AMAX1(-VOLW1N*XNPXX,VOLP1ZN))
      FLQ2=FLQ2+AMIN1(0.0,AMAX1(-VOLW1N*XNPXX,VOLP1ZN))
      ENDIF
      IF(FLQR.GT.0.0)THEN
      HFLQR=4.19*TK1(0,NY,NX)*FLQR
      ELSE
      HFLQR=4.19*TK1(NUM(NY,NX),NY,NX)*FLQR
      ENDIF
C
C     AGGREGATE LITTER-SOIL WATER AND HEAT FLUXES
C
C     FLWL,HFLWL=micropore water,heat flux (m3 t-1,MJ t-1)
C     FLWRL,HFLWRL=total litter water,heat flux (m3 t-1,MJ t-1)
C     FLWRM=litter-soil water flux for solute transfer in ‘trnsfr.f’
C        (m3 t-1)     
C
      FLWRL(NY,NX)=FLWRL(NY,NX)-FLQR
      HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQR
      FLWL(3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)+FLQR
      HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQR
      FLWRM(M,NY,NX)=FLWRM(M,NY,NX)+FLQR
      VOLW10=VOLW10-FLQR
      VOLW1N=VOLW1N+FLQR
      VOLP10=AMAX1(0.0,VOLWRX(NY,NX)-VOLW10)
      VOLP1ZN=VOLA1(NUM(NY,NX),NY,NX)-VOLW1N-VOLI1(NUM(NY,NX),NY,NX)
      VOLP1N=AMAX1(0.0,VOLP1ZN)
C     IF(NY.EQ.5)THEN
C     WRITE(*,4322)'FLQR',I,J,NFZ,M,NN,NX,NY,NUM(NY,NX),K0,K1
C    2,FLWRM(M,NY,NX),FLQR,FLQX,FLQZ,FLQ2,VOLW1(0,NY,NX)
C    3,VOLW1(NUM(NY,NX),NY,NX),FLWRL(NY,NX),FLWL(3,NUM(NY,NX),NY,NX) 
C    2,THETWR,VOLWRX(NY,NX),VOLR(NY,NX)
C    3,THETW1,PSISM10,PSISM1N,PSIST0,PSIST1
C    4,VOLP1(0,NY,NX),VOLP1(NUM(NY,NX),NY,NX)
C    3,PSISO(0,NY,NX),PSISO(NUM(NY,NX),NY,NX)
C    4,CVRDW(NY,NX),CNDR,CND1,AVCNDR,FKSAT
C    3,POROS0(NY,NX),VOLW1(0,NY,NX),VOLI1(0,NY,NX)  
C    2,FLWL(3,NUM(NY,NX),NY,NX)
C    3,VOLW10,VOLW1N,VOLP10,VOLP1N,VOLP1ZN 
C    4,THETWX(0,NY,NX),THETWX(NUM(NY,NX),NY,NX)
C    4,THETS(0,NY,NX),THETS(NUM(NY,NX),NY,NX)
C    4,THETWR,THETW1,XVOLT(NY,NX),CDPTH(NU(NY,NX)-1,NY,NX)
C    3,HFLQR,HFLWRL(NY,NX),HFLWL(3,NUM(NY,NX),NY,NX) 
C    6,THETWR,VHCP1(0,NY,NX),VHCPRXG(NY,NX)
C    2,FLWLYH,VOLX(NUM(NY,NX),NY,NX),VOLA1(NUM(NY,NX),NY,NX)
4322  FORMAT(A8,10I4,40E12.4) 
C     ENDIF
6000  CONTINUE
      ELSE
      FLQR=AMAX1(0.0,(VOLW1(0,NY,NX)-0.01/AREA(3,0,NY,NX))*XNPXX)
      HFLQR=4.19*TK1(0,NY,NX)*FLQR
      FLWRL(NY,NX)=FLWRL(NY,NX)-FLQR
      HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQR
      FLWL(3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)+FLQR
      HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQR
      FLWRM(M,NY,NX)=FLWRM(M,NY,NX)+FLQR
C     IF(NY.EQ.5)THEN
C     WRITE(*,4323)'FLQR0',I,J,NFZ,M,NX,NY,NUM(NY,NX),FLWRM(M,NY,NX)
C    2,FLQR,VOLW1(0,NY,NX),VOLW1(NUM(NY,NX),NY,NX)
C    3,FLWRL(NY,NX),FLWL(3,NUM(NY,NX),NY,NX),CVRDW(NY,NX),BAREW(NY,NX)
C    4,BARE(NY,NX),XVOLT(NY,NX),VOLWD(NY,NX),ORGC(0,NY,NX)
4323  FORMAT(A8,7I4,30E12.4) 
C     ENDIF
      ENDIF
C
C     OVERLAND FLOW INTO SOIL MACROPORES WHEN WATER STORAGE CAPACITY
C     OF THE LITTER IS EXCEEDED
C
C     VOLPH1=air-filled macroporosity (m3)
C     XVOLW=surface water in excess of litter water retention
C        capacity (m3) 
C     FLQHR,HFLQHR=water,convective heat from litter to macropores
C        (m3 t-1,MJ t-1)
C     FLWHL,HFLWL=total macropore water,heat flux (m3 t-1,MJ t-1)
C     FLWRL,HFLWRL=total litter water,heat flux (m3 t-1,MJ t-1)
C
      IF(VOLPH1(NUM(NY,NX),NY,NX).GT.0.0
     2.AND.XVOLW(NY,NX).GT.0.0)THEN
      FLQHR=AMIN1(XVOLW(NY,NX)*XNPXX,VOLPH1(NUM(NY,NX),NY,NX))
      HFLQHR=FLQHR*4.19*TK1(0,NY,NX)
      FLWHL(3,NUM(NY,NX),NY,NX)=FLWHL(3,NUM(NY,NX),NY,NX)+FLQHR
      HFLWL(3,NUM(NY,NX),NY,NX)=HFLWL(3,NUM(NY,NX),NY,NX)+HFLQHR
      FLWRL(NY,NX)=FLWRL(NY,NX)-FLQHR
      HFLWRL(NY,NX)=HFLWRL(NY,NX)-HFLQHR
C     IF(NY.EQ.5)THEN
C     WRITE(*,4357)'FLQHR',I,J,M,NX,NY,NUM(NY,NX),FLQHR,FLWRL(NY,NX)
C    2,FLWHL(3,NUM(NY,NX),NY,NX),VOLPH1(NUM(NY,NX),NY,NX)
C    3,VOLWH1(NUM(NY,NX),NY,NX),XVOLW(NY,NX)
C    4,VOLW1(0,NY,NX),VOLWRX(NY,NX)
C    4,HFLQHR,HFLWRL(NY,NX),HFLWL(3,NUM(NY,NX),NY,NX),TK1(0,NY,NX)
4357  FORMAT(A8,6I4,40E12.4)
C     ENDIF
      ENDIF
C
C     THICKNESS OF WATER FILMS IN LITTER AND SOIL SURFACE 
C     FROM WATER POTENTIALS FOR GAS EXCHANGE IN TRNSFR.F
C
C     VHCP1,VHCPRX=litter current,minimum heat capacity (MJ K-1)
C     FILM=litter (0) and soil (NUM) water film thickness for nutrient
C        and gas uptake in ‘nitro.f’ and uptake.f’ (m)
C     PSISM1=litter matric potential (MPa)
C
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      FILM(M,0,NY,NX)=AMAX1(1.0E-06
     2,0.5*EXP(-13.650-0.857*LOG(-PSISM1(0,NY,NX))))
      ELSE
      FILM(M,0,NY,NX)=1.0E-06
      ENDIF
      FILM(M,NUM(NY,NX),NY,NX)=AMAX1(1.0E-06
     2,0.5*EXP(-13.650-0.857*LOG(-PSISM1(NUM(NY,NX),NY,NX))))
C
C     OVERLAND FLOW WHEN WATER STORAGE CAPACITY
C     OF THE SOIL SURFACE PLUS MACROPORES IS EXCEEDED
C
C     SURFACE WATER FLUX
C
C     VOLW10,VOLI10=litter water,ice volume (m3)
C     XVOLW,XVOLI=surface water,ice in excess of litter water retention
C        capacity (m3) 
C     VOLWRX=litter water retention capacity (m3)
C     VOLWG=ground surface water retention capacity (m3)
C     VX=ponded water volume above surface retention capacity (m3)
C     D,R=depth,perimeter of runoff (m)
C     SLOPE=ground surface slope from ‘starts.f’ (m m-1) 
C     V=runoff velocity (m s-1)
C     DIST=distance between source,destination (m)
C     ZM=surface roughness height for runoff (m)
C     Q=runoff from Mannings equation (m3 t-1)
C     QRM=downslope runoff used in ‘trnsfr.f’ and ‘trnsfrs.f’ (m3 t-1)
C     QRV=downslope runoff velocity used in ‘erosion.f’ (m s-1)
C
      VOLW10=VOLW1(0,NY,NX)+FLWRL(NY,NX)+WFLVR(NY,NX)
     2+WFLFR(NY,NX) 
      VOLI10=VOLI1(0,NY,NX)-WFLFR(NY,NX)/DENSI
      TVOLWI=VOLW10+VOLI10
      XVOLT(NY,NX)=AMAX1(0.0,TVOLWI-VOLWRX(NY,NX))
      IF(TVOLWI.GT.ZEROS(NY,NX))THEN
      VOLWRZ=VOLW10/TVOLWI*VOLWRX(NY,NX)
      VOLIRZ=VOLI10/TVOLWI*VOLWRX(NY,NX)
      XVOLW(NY,NX)=AMAX1(0.0,VOLW10-VOLWRZ)
      XVOLI(NY,NX)=AMAX1(0.0,VOLI10-VOLIRZ)
      ELSE
      XVOLW(NY,NX)=0.0
      XVOLI(NY,NX)=0.0
      ENDIF
      XVOLTM(M,NY,NX)=XVOLT(NY,NX)
      XVOLWM(M,NY,NX)=XVOLW(NY,NX)
      XVOLIM(M,NY,NX)=XVOLI(NY,NX)
      IF(XVOLW(NY,NX).GT.VOLWG(NY,NX))THEN
      VX=XVOLW(NY,NX)-VOLWG(NY,NX)
      D=VX/AREA(3,0,NY,NX)
      R=D/2.828
      V=R**0.67*SQRT(SLOPE(0,NY,NX))/ZM(NY,NX) 
      Q=V*D*AREA(3,NUM(NY,NX),NY,NX)*3.6E+03*XNPHX 
      QRM(M,NY,NX)=AMIN1(Q,VX*XNPXX,VOLW2(N3,N2,N1)*XNPXX)
C    2*XVOLW(NY,NX)/XVOLT(NY,NX)
C     IF(NY.EQ.6.OR.NY.EQ.7)THEN
C     WRITE(*,5554)'QRM',I,J,NFZ,M,NY,NX
C    2,QRM(M,NY,NX),QRV(M,NY,NX) 
C    2,Q,V,D,VX,XNPXX,SLOPE(0,NY,NX),ZM(NY,NX),FLYM
C    3,XVOLW(NY,NX),XVOLI(NY,NX),XVOLT(NY,NX),VOLWG(NY,NX)
C    4,VOLW1(0,NY,NX),VOLI1(0,NY,NX)
C    4,VOLW10,FLWRL(NY,NX),WFLVR(NY,NX),WFLFR(NY,NX)
5554  FORMAT(A8,6I4,20E12.4)
C     ENDIF
      QRV(M,NY,NX)=V
      ELSE
      QRM(M,NY,NX)=0.0 
      QRV(M,NY,NX)=0.0
      ENDIF
C
C     DOWNSLOPE SNOW REDISTRIBUTION FROM SNOWPACK
C
C     QSX=snow transfer fraction (t-1)
C     UA=wind speed from weather file (m h-1)
C     QSM,QWM,QIM,QST=snow,water,ice,total transfer (m3 t-1)
C     VOLS0,VOLW0,VOLI0=snow,water,ice volume (m3)    
C     XNPHX=time step from ‘wthr.f’(h t-1)
C
      IF(DPTHS(NY,NX).GT.ZERO)THEN
      QSX=1.0E-07*UA(NY,NX)*XNPHX
      QSM(NY,NX)=QSX*VOLS0(1,NY,NX)
      QWM(NY,NX)=QSX*VOLW0(1,NY,NX)
      QIM(NY,NX)=QSX*VOLI0(1,NY,NX)
      QST(M,NY,NX)=QSM(NY,NX)+QWM(NY,NX)+QIM(NY,NX)
C     WRITE(*,5554)'QST',I,J,NFZ,M,NX,NY,QSX,QST(M,NY,NX)
C    2,UA(NY,NX),SLOPE(0,NY,NX),XNPHX,QSM(NY,NX),QWM(NY,NX),QIM(NY,NX)
C    3,VOLS0(1,NY,NX),VOLW0(1,NY,NX),VOLI0(1,NY,NX)
      ELSE
      QSM(NY,NX)=0.0
      QWM(NY,NX)=0.0
      QIM(NY,NX)=0.0
      QST(M,NY,NX)=0.0
      ENDIF
C
C     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
C
      N1=NX
      N2=NY
      DO 4310 N=1,2
      DO 4305 NN=1,2
      IF(N.EQ.1)THEN
C
C     N2,N1=NY,NX of source grid cell
C     N5,N4=NY,NX of destination grid cell E or S
C     N5B,N4B=NY,NX of destination grid cell W or N
C     NN=boundary:N=1:NN=1 east,NN=2 west, N=2:NN=1 south,NN=2 north
C
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
C     XVOLT,XVOLW=surface water+ice,water in excess of litter
C        retention capacity in destination grid cell (m3)
C     ALT1,ALT2,ALTB=elevation of source,destination grid cell 
C        E or S,N or W (m)
C     QRQ1=equilibrium runoff between grid cells (m3 t-1)
C     FSLOPE=partitions surface water flow in (N=1)EW,(N=2)NS
C       directions from ‘starts.f’
C     QR1,HQR1=runoff, convective heat from runoff (m3 t-1,MJ t-1)
C     QR,HQR=aggregated runoff, convective heat from runoff
C        used in ‘redist.f’ (m3 t-1,MJ t-1)
C     QRM=downslope runoff used in ‘erosion.f’, ‘trnsfr.f’ 
C        and ‘trnsfrs.f’ (m3 t-1)
C     QRMN=EW(N=1),NS(N=2) runoff used in ‘erosion.f’, ‘trnsfr.f’ 
C        and ‘trnsfrs.f’ (m3 t-1)
C     IFLBM=runoff direction flag used in ‘erosion.f’, ‘trnsfr.f’ 
C        and ‘trnsfrs.f’ (0 = E or S, 1 = W or N)
C
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
      ALT1=ALTG(N2,N1)+XVOLT(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)
C
C     EAST OR SOUTH RUNOFF
C
      IF(NN.EQ.1)THEN
      ALT2=ALTG(N5,N4)+XVOLT(N5,N4)/AREA(3,NU(N5,N4),N5,N4)
      IF(ALT1.GT.ALT2)THEN
      QRQ1=AMAX1(0.0,(0.5*(ALT1-ALT2)*AREA(3,NUM(N2,N1),N2,N1)
     2*AREA(3,NU(N5,N4),N5,N4))
     4/(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5,N4),N5,N4)))
      QR1(N,2,N5,N4)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
      HQR1(N,2,N5,N4)=4.19*TK1(0,N2,N1)*QR1(N,2,N5,N4)
      QR(N,2,N5,N4)=QR(N,2,N5,N4)+QR1(N,2,N5,N4)
      HQR(N,2,N5,N4)=HQR(N,2,N5,N4)+HQR1(N,2,N5,N4)
      QRMN(M,N,2,N5,N4)=QR1(N,2,N5,N4)
      IFLBM(M,N,2,N5,N4)=0
      ELSE
      QR1(N,2,N5,N4)=0.0
      HQR1(N,2,N5,N4)=0.0
      QRMN(M,N,2,N5,N4)=0.0
      IFLBM(M,N,2,N5,N4)=0
      ENDIF
C     IF(N2.EQ.6.OR.N2.EQ.7)THEN
C     WRITE(*,5555)'QRFOR',I,J,NFZ,M,N1,N2,N4,N5,N,NN
C    2,QRM(M,N2,N1),QR1(N,2,N5,N4),QR(N,2,N5,N4) 
C    2,ALT1,ALT2,ALTG(N2,N1),ALTG(N5,N4),ALT(N2,N1),ALT(N5,N4)
C    3,QRQ1,FSLOPE(N,N2,N1),VOLW1(0,N2,N1),VOLW1(0,N5,N4)
C    4,XVOLT(N2,N1),XVOLT(N5,N4)
5555  FORMAT(A8,10I4,30E14.6)
C     ENDIF
      ENDIF
C
C     WEST OR NORTH RUNOFF
C
      IF(NN.EQ.2)THEN
      IF(N4B.GT.0.AND.N5B.GT.0)THEN
      ALTB=ALTG(N5B,N4B)+XVOLT(N5B,N4B)/AREA(3,NU(N5,N4B),N5B,N4B)
      IF(ALT1.GT.ALTB)THEN
      QRQ1=AMAX1(0.0,(0.5*(ALT1-ALTB)*AREA(3,NUM(N2,N1),N2,N1)
     2*AREA(3,NU(N5B,N4B),N5B,N4B))
     4/(AREA(3,NUM(N2,N1),N2,N1)+AREA(3,NU(N5B,N4B),N5B,N4B)))
      QR1(N,1,N5B,N4B)=AMIN1(QRQ1,QRM(M,N2,N1))*FSLOPE(N,N2,N1)
      HQR1(N,1,N5B,N4B)=4.19*TK1(0,N2,N1)*QR1(N,1,N5B,N4B)
      QR(N,1,N5B,N4B)=QR(N,1,N5B,N4B)+QR1(N,1,N5B,N4B)
      HQR(N,1,N5B,N4B)=HQR(N,1,N5B,N4B)+HQR1(N,1,N5B,N4B)
      QRMN(M,N,1,N5B,N4B)=QR1(N,1,N5B,N4B)
      IFLBM(M,N,1,N5B,N4B)=1
      ELSE
      QR1(N,1,N5B,N4B)=0.0
      HQR1(N,1,N5B,N4B)=0.0
      QRMN(M,N,1,N5B,N4B)=0.0
      IFLBM(M,N,1,N5B,N4B)=1
      ENDIF
C     IF(N2.EQ.6.OR.N2.EQ.7)THEN
C     WRITE(*,5556)'QRBAK',I,J,NFZ,M,N1,N2,N4B,N5B,N,NN
C    2,QRM(M,N2,N1),QR1(N,1,N5B,N4B),QR(N,1,N5B,N4B) 
C    2,ALT1,ALTB,ALTG(N2,N1),ALTG(N5B,N4B),QRQ1,FSLOPE(N,N2,N1)
C    4,VOLW1(0,N2,N1),VOLW1(0,N5B,N4B),XVOLT(N2,N1),XVOLT(N5B,N4B)
5556  FORMAT(A8,10I4,30E14.6)
C     ENDIF
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
C     IF(N4B.GT.0.AND.N5B.GT.0.AND.(N2.EQ.6.OR.N2.EQ.7))THEN
C     WRITE(*,5557)'QREND',I,J,NFZ,M,N1,N2,N4,N5,N4B,N5B,N,NN
C    2,IFLBM(M,N,NN,N5,N4)
C    2,QRM(M,N2,N1),QR1(N,NN,N5,N4),QR(N,NN,N5,N4)
C    3,QR1(N,NN,N5B,N4B),QR(N,NN,N5B,N4B),ALT1,ALT2,ALTB 
5557  FORMAT(A8,13I4,30E14.6)
C     ENDIF
C
C     SNOW REDISTRIBUTION FROM SNOWPACK AMONG GRID CELLS
C
C     N2,N1=NY,NX of source grid cell
C     N5,N4=NY,NX of destination grid cell
C     ALTS1,ALTS2=elevation of source,destination snowpack surfaces (m)
C     QS1,QW1,QI1=E,S,W,N snow,water,ice transfer (m3 t-1)
C     QSM,QWM,QIM=downslope snow,water,ice transfer (m3 t-1)
C     FSLOPE=partitions surface water flow in (N=1)EW,(N=2)NS
C       directions from ‘starts.f’
C     HQS1=convective heat transfer from snow,water,ice transfer 
C        (MJ t-1)
C     QS,QW,QI=aggregated snow,water,ice transfer used in ‘redist.f’ 
C        (m3 t-1)
C     HQS=aggregated convective heat from snow,water,ice transfer
C        used in ‘redist.f’ (MJ t-1)
C     QSTN=EW(N=1),NS(N=2) snow transfer for solute flux calculation in
C        ‘trnsfr.f’ and ‘trnsfrs.f’ (m3 t-1)
C     IFLBMS=snow drift direction flag used in ‘trnsfr.f’ 
C        and ‘trnsfrs.f’ (0 = E or S, 1 = W or N)
C
      IF(QST(M,N2,N1).GT.ZEROS(N2,N1))THEN
      ALTS1=ALTG(N2,N1)+DPTHS(N2,N1)
C
C     EAST OR SOUTH DRIFT
C
      IF(NN.EQ.1)THEN
      ALTS2=ALTG(N5,N4)+DPTHS(N5,N4)
      IF(ALTS1.GT.ALTS2)THEN
      QS1(N,2,N5,N4)=QSM(N2,N1)*FSLOPE(N,N2,N1)
      QW1(N,2,N5,N4)=QWM(N2,N1)*FSLOPE(N,N2,N1)
      QI1(N,2,N5,N4)=QIM(N2,N1)*FSLOPE(N,N2,N1)
      HQS1(N,2,N5,N4)=TK0(1,N2,N1)*(2.095*QS1(N,2,N5,N4)
     2+4.19*QW1(N,2,N5,N4)+1.9274*QI1(N,2,N5,N4))
      QS(N,2,N5,N4)=QS(N,2,N5,N4)+QS1(N,2,N5,N4)
      QW(N,2,N5,N4)=QW(N,2,N5,N4)+QW1(N,2,N5,N4)
      QI(N,2,N5,N4)=QI(N,2,N5,N4)+QI1(N,2,N5,N4)
      HQS(N,2,N5,N4)=HQS(N,2,N5,N4)+HQS1(N,2,N5,N4)
      QSTN(M,N,2,N5,N4)=QS1(N,2,N5,N4)
     2+QW1(N,2,N5,N4)+QI1(N,2,N5,N4)
      IFLBMS(M,N,2,N5,N4)=0
C     WRITE(*,5555)'QSFOR',I,J,NFZ,M,N1,N2,N4,N5,N,NN
C    2,QSM(NY,NX),QWM(NY,NX),QIM(NY,NX) 
C    3,QS1(N,2,N5,N4),QW1(N,2,N5,N4),QI1(N,2,N5,N4)
C    4,FSLOPE(N,N2,N1)
      ELSE
      QS1(N,2,N5,N4)=0.0
      QW1(N,2,N5,N4)=0.0
      QI1(N,2,N5,N4)=0.0
      HQS1(N,2,N5,N4)=0.0
      QSTN(M,N,2,N5,N4)=0.0
      IFLBMS(M,N,2,N5,N4)=0
      ENDIF
      ENDIF
C
C     WEST OR NORTH DRIFT
C
      IF(NN.EQ.2)THEN
      IF(N4B.GT.0.AND.N5B.GT.0)THEN
      ALTSB=ALTG(N5B,N4B)+DPTHS(N5B,N4B) 
      IF(ALTS1.GT.ALTSB)THEN
      QS1(N,1,N5B,N4B)=QSM(N2,N1)*FSLOPE(N,N2,N1)
      QW1(N,1,N5B,N4B)=QWM(N2,N1)*FSLOPE(N,N2,N1)
      QI1(N,1,N5B,N4B)=QIM(N2,N1)*FSLOPE(N,N2,N1)
      HQS1(N,1,N5B,N4B)=TK0(1,N2,N1)*(2.095*QS1(N,1,N5B,N4B)
     2+4.19*QW1(N,1,N5B,N4B)+1.9274*QI1(N,1,N5B,N4B))
      QS(N,1,N5B,N4B)=QS(N,1,N5B,N4B)+QS1(N,1,N5B,N4B)
      QW(N,1,N5B,N4B)=QW(N,1,N5B,N4B)+QW1(N,1,N5B,N4B)
      QI(N,1,N5B,N4B)=QI(N,1,N5B,N4B)+QI1(N,1,N5B,N4B)
      HQS(N,1,N5B,N4B)=HQS(N,1,N5B,N4B)+HQS1(N,1,N5B,N4B)
      QSTN(M,N,1,N5B,N4B)=QS1(N,1,N5B,N4B)
     2+QW1(N,1,N5B,N4B)+QI1(N,1,N5B,N4B)
      IFLBMS(M,N,1,N5B,N4B)=1
      ELSE
      QS1(N,1,N5B,N4B)=0.0
      QW1(N,1,N5B,N4B)=0.0
      QI1(N,1,N5B,N4B)=0.0
      HQS1(N,1,N5B,N4B)=0.0
      QSTN(M,N,1,N5B,N4B)=0.0
      IFLBMS(M,N,1,N5B,N4B)=1
      ENDIF
      ENDIF
      ENDIF
      ELSE
      QS1(N,2,N5,N4)=0.0
      QW1(N,2,N5,N4)=0.0
      QI1(N,2,N5,N4)=0.0
      HQS1(N,2,N5,N4)=0.0
      QSTN(M,N,2,N5,N4)=0.0
      IFLBMS(M,N,2,N5,N4)=0
      IF(N4B.GT.0.AND.N5B.GT.0)THEN
      QS1(N,1,N5B,N4B)=0.0
      QW1(N,1,N5B,N4B)=0.0
      QI1(N,1,N5B,N4B)=0.0
      HQS1(N,1,N5B,N4B)=0.0
      QSTN(M,N,1,N5B,N4B)=0.0
      IFLBMS(M,N,1,N5B,N4B)=0
      ENDIF
      ENDIF
4305  CONTINUE
4310  CONTINUE
C
C     ACCUMULATED WATER, VAPOR AND HEAT FLUXES THROUGH 
C     SURFACE RESIDUE AND SOIL SURFACE
C
C     WFLVR,HFLVR=surface litter evaporation-condensation,latent 
C        heat flux (m3 t-1,MJ t-1)
C     WFLFR,HFLFR=surface litter freeze-thaw,latent 
C        heat flux (m3 t-1,MJ t-1)
C     XW*,XH*=aggregated water,heat fluxes used in ‘redist.f’
C        (m3 t-1,MJ t-1)
C     FLW,FLWH,FLV,HFLW=soil surface micropore and macropore water,
C        vapor,heat fluxes (m3 t-1)
C     FLWX=FLW accounting for wetting front during infiltration
C        (m3)
C     FLWR,FLVR,HFLWR=litter water,vapor,heat fluxes (m3 t-1,MJ t-1) 
C     FLSW,FLSV,FLSWH=water,vapor from snowpack to soil
C        micropores,macropores (m3 t-1)
C     HEATI,HEATE,HEATV,HEATS,HEATG=total net radiation,latent,
C        convective,sensible,storage heat at all ground surfaces 
C        (MJ t-1)
C     TEVAPG=total evaporation at all ground surfaces (m3 t-1)
C     FLWM,FLWHM=water flux into soil micropore,macropore 
C        for use in ‘trnsfr.f’ (m3 h-1)
C
      XWFLVR(NY,NX)=XWFLVR(NY,NX)+WFLVR(NY,NX)
      XHFLVR(NY,NX)=XHFLVR(NY,NX)+HFLVR(NY,NX)
      XWFLFR(NY,NX)=XWFLFR(NY,NX)+WFLFR(NY,NX)
      XHFLFR(NY,NX)=XHFLFR(NY,NX)+HFLFR(NY,NX)
      FLW(3,NUM(NY,NX),NY,NX)=FLW(3,NUM(NY,NX),NY,NX)
     2+FLWL(3,NUM(NY,NX),NY,NX)
      FLV(3,NUM(NY,NX),NY,NX)=FLV(3,NUM(NY,NX),NY,NX)
     2+FLVL(3,NUM(NY,NX),NY,NX)
      FLWX(3,NUM(NY,NX),NY,NX)=FLWX(3,NUM(NY,NX),NY,NX)
     2+FLWLX(3,NUM(NY,NX),NY,NX)
      FLWH(3,NUM(NY,NX),NY,NX)=FLWH(3,NUM(NY,NX),NY,NX)
     2+FLWHL(3,NUM(NY,NX),NY,NX)
      HFLW(3,NUM(NY,NX),NY,NX)=HFLW(3,NUM(NY,NX),NY,NX)
     2+HFLWL(3,NUM(NY,NX),NY,NX)
      FLWR(NY,NX)=FLWR(NY,NX)+FLWRL(NY,NX)
      FLVR(NY,NX)=FLVR(NY,NX)+FLVRL(NY,NX)
      HFLWR(NY,NX)=HFLWR(NY,NX)+HFLWRL(NY,NX)
      HEATI(NY,NX)=HEATI(NY,NX)+RFLXG+RFLXR+RFLXW
      HEATS(NY,NX)=HEATS(NY,NX)+SFLXG+SFLXR+SFLXW
      HEATE(NY,NX)=HEATE(NY,NX)+EFLXG+EFLXR+EFLXW
      HEATV(NY,NX)=HEATV(NY,NX)+VFLXG+VFLXR+VFLXW
      HEATH(NY,NX)=HEATH(NY,NX)+RFLXG+RFLXR+RFLXW 
     2+SFLXG+SFLXR+SFLXW+EFLXG+EFLXR+EFLXW+VFLXG+VFLXR+VFLXW
      TEVAPG(NY,NX)=TEVAPG(NY,NX)+EVAPG(NY,NX)+EVAPR(NY,NX)
     2+EVAP0(NY,NX) 
      FLWM(M,3,NUM(NY,NX),NY,NX)=FLWL(3,NUM(NY,NX),NY,NX)
      FLWHM(M,3,NUM(NY,NX),NY,NX)=FLWHL(3,NUM(NY,NX),NY,NX)
C
C     AERODYNAMIC ENERGY EXCHANGE BETWEEN GROUND SURFACE AIR
C     AND ATMOSPHERE
C
C     TKQC,VPQC=live canopy air temperature,vapor
C        concentration (K,m3 m-3)
C     TKQD,VPQD=dead canopy air temperature,vapor
C        concentration (K,m3 m-3)
C     TKQG,VPQG=ground surface air temperature,vapor
C        concentration (K,m3 m-3)
C     DTKQC,DVPQC=live canopy-ground surface air temperature,vapor
C        concentration difference (K,m3 m-3)
C     DTKQD,DVPQD=dead canopy-ground surface air temperature,vapor
C        concentration difference (K,m3 m-3)
C     RI=Richardson number
C     RAGC,RAGD=canopy aerodynamic resistance below 
C        PFT canopy height (h m-1)
C     ZC,ZG=canopy live, standing dead height (m)
C        resistance from ‘hour1.f’ 
C     PAREX,PARSX=terms used to calculate boundary layer
C        conductance for latent,sensible heat from ‘hour1.f’
C        (m2 h t-1,MJ h m-1 K-1 t-1)
C     FRADT=fraction of radiation received by all PFT canopies 
C        from ‘hour1.f’
C     PARSCM,PARECM=live canopy-ground conductance for vapor,heat
C        exchange(m t-1,m MJ K-1 t-1)
C     PARSDM,PAREDM=dead canopy-ground conductance for vapor,heat
C        exchange(m t-1,m MJ K-1 t-1)
C     SHGCM,VPGCM=total live canopy sensible heat,vapor exchange 
C        (MJ t-1,m3 t-1)
C     SHGDM,VPGDM=total dead canopy sensible heat,vapor exchange 
C        (MJ t-1,m3 t-1)
C     FLAIP,FLAIQ=fraction of total canopy + standing dead radiation
C        received by each PFT canopy, standing dead from ‘hour1.f’ 
C     SFLXG,SFLXR,SFLXW=sensible heat flux at soil,litter,snowpack
C        surfaces (MJ t-1)
C     HFLXH=heat from combustion in canopy air (MJ t-1)
C     SFLXA,EVAPA=sensible heat,vapor exchange between surface air 
C        and atmosphere (MJ t-1,m3 t-1)
C     RATG=aerodynamic+canopy boundary layer resistance (h m-1)
C     VPAM=atmospheric vapor concentration from ‘hour1.f’(m3 m-3)
C     SHGCM,SHGDM=total ground surface-live canopy,standing dead
C        sensible heat flux (MJ t-1)
C     VPGCM,VPGDM=total ground surface-live canopy,standing dead 
C        vapor flux (m3 t-1)
C     TKQG,TKQC,TKQD=air temperature at ground surface,
C        live canopy,standing dead (K)
C     VPQG,VPQC,VPQD=vapor concentration at ground surface,
C        canopy,standing dead (m3 m-3)
C     VHCPQ,EVAPQ=canopy heat capacity,volume used to calculate 
C        canopy air temperature,vapor concentration from ‘hour1.f’ 
C        (MJ K-1,m3)
C     TKQGX,VPQGX=TKQC,VPQC used in ‘redist.f’ (K,m3 m-3)
C
      DO 920 NZ=1,NP(NY,NX)
      DTKQC(NZ)=TKQC(NZ,NY,NX)-TKQG(M,NY,NX)
      DVPQC(NZ)=VPQC(NZ,NY,NX)-VPQG(M,NY,NX)
      RI=AMAX1(RIX,AMIN1(RIY
     2,RIBX(NY,NX)/TKQC(NZ,NY,NX)*DTKQC(NZ)))
      RIC=1.0-10.0*RI 
      DTKQD(NZ)=TKQD(NZ,NY,NX)-TKQG(M,NY,NX) 
      DVPQD(NZ)=VPQD(NZ,NY,NX)-VPQG(M,NY,NX) 
      RI=AMAX1(RIX,AMIN1(RIY
     2,RIBX(NY,NX)/TKQD(NZ,NY,NX)*DTKQD(NZ)))
      RID=1.0-10.0*RI 
      IF(ZC(NZ,NY,NX).EQ.ZT(NY,NX))THEN
      RAGCX=RACGX
      ELSEIF(ZC(NZ,NY,NX).GT.ZERO)THEN
      RAGCX=ZC(NZ,NY,NX)/ZT(NY,NX)*RACGX
      ELSE
      RAGCX=0.0
      ENDIF
      RAGC(M,NZ,NY,NX)=AMIN1(RACGZ,AMAX1(RACGM,RAGCX/RIC))
      IF(ZG(NZ,NY,NX).EQ.ZT(NY,NX))THEN
      RAGDX=RACGX
      ELSEIF(ZG(NZ,NY,NX).GT.ZERO)THEN
      RAGDX=ZG(NZ,NY,NX)/ZT(NY,NX)*RACGX
      ELSE
      RAGDX=0.0
      ENDIF
      RAGD(M,NZ,NY,NX)=AMIN1(RACGZ,AMAX1(RACGM,RAGDX/RID))
C     WRITE(*,102)'RAGC',I,J,NFZ,M,NY,NX,NZ,RACGX,RAGCX,RAGDX
C    2,RAGC(M,NZ,NY,NX),RAGD(M,NZ,NY,NX),ZC(NZ,NY,NX),ZG(NZ,NY,NX)
C    3,RIC,RID
102   FORMAT(A8,7I4,20E12.4) 
920   CONTINUE
      IF(M.LT.NPH)THEN
      SHGCM=0.0
      SHGDM=0.0
      VPGCM=0.0
      VPGDM=0.0
      DO 921 NZ=1,NP(NY,NX)
      PARSCM=PARSX(NY,NX)/RAGC(M,NZ,NY,NX)*FRADT(NY,NX)
      PARECM=PAREX(NY,NX)/RAGC(M,NZ,NY,NX)*FRADT(NY,NX)
      PARSDM=PARSX(NY,NX)/RAGD(M,NZ,NY,NX)*FRADT(NY,NX)
      PAREDM=PAREX(NY,NX)/RAGD(M,NZ,NY,NX)*FRADT(NY,NX)
      SHGCM=SHGCM+PARSCM*DTKQC(NZ)*FLAIP(NZ,NY,NX)
      VPGCM=VPGCM+PARECM*DVPQC(NZ)*FLAIP(NZ,NY,NX)
      SHGDM=SHGDM+PARSDM*DTKQD(NZ)*FLAIQ(NZ,NY,NX)
      VPGDM=VPGDM+PAREDM*DVPQD(NZ)*FLAIQ(NZ,NY,NX)
921   CONTINUE
      TSHGM2=SHGCM+SHGDM-SFLXG-SFLXR-SFLXW+HFLXH(NY,NX) 
      TKQG2=TKQG(M,NY,NX)+TSHGM2/VHCPQ(NY,NX)
      DTKQ2=TKAM(NY,NX)-TKQG2 
      SFLXA=DTKQ2*PARSX(NY,NX)/RATG
      TSHGM=SFLXA+TSHGM2 
      TKQG(M+1,NY,NX)=TKQG(M,NY,NX)+TSHGM/VHCPQ(NY,NX)
      TEVGM2=VPGCM+VPGDM 
     2-EVAPG(NY,NX)-EVAPR(NY,NX)-EVAP0(NY,NX)
      VPQG2=AMAX1(0.0,VPQG(M,NY,NX)+TEVGM2/EVAPQ(NY,NX))
      DVPQ2=VPAM(NY,NX)-VPQG2
      EVAPA=DVPQ2*PAREX(NY,NX)/RATG 
      TEVGM=EVAPA+TEVGM2
      VPSG=2.173E-03/TKQG(M+1,NY,NX)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKQG(M+1,NY,NX)))
      VPQG(M+1,NY,NX)=AMAX1(0.0,AMIN1(VPSG
     2,VPQG(M,NY,NX)+TEVGM/EVAPQ(NY,NX)))
      IF(M+1.EQ.NPH)THEN
      TKQGX(NY,NX)=TKQG(M+1,NY,NX)
      VPQGX(NY,NX)=VPQG(M+1,NY,NX)
      ENDIF
C     IF(J.EQ.14)THEN
C     WRITE(*,3114)'TKQG',I,J,NFZ,M,NX,NY,NUM(NY,NX)
C    2,TKQG(M+1,NY,NX),TKQG(M,NY,NX),TKGS(M,NY,NX),TKAM(NY,NX) 
C    3,TSHGM,TSHGM2,SFLXA,SHGCM,SHGDM,SFLXG,SFLXR,SFLXW,HFLXH(NY,NX) 
C    3,VPSG,VPQG(M,NY,NX),VPAM(NY,NX),EVAPA
C    3,TK1(0,NY,NX),TK1(NUM(NY,NX),NY,NX),RAB(NY,NX),RACG(M,NY,NX) 
C    4,BAREW(NY,NX),VHCP1(0,NY,NX),VHCP1(NUM(NY,NX),NY,NX) 
C    3,VHCPRXG(NY,NX),HFLXH(NY,NX),PSISO(0,NY,NX)
C    4,DVPQ2,VPQG2,TSHGM2,TSHGM,SFLXA,SHGCM,SHGDM,SFLXG,SFLXR,SFLXW 
C    5,PAREX(NY,NX),RATG,PAREX(NY,NX)/RATG,VHCPQ(NY,NX)
C    5,VOLW10,VOLI10,FLWRL(NY,NX),WFLVR(NY,NX),WFLFR(NY,NX) 
C    3,FRADT(NY,NX),VHCPQ(NY,NX),EVAPQ(NY,NX)
C    2,RAB(NY,NX),RACGX,RACG(M,NY,NX),RAH,RARH,RZR(NY,NX)
C    4,(RAGC(M,NZ,NY,NX),NZ=1,NP(NY,NX))
C    5,(RAGD(M,NZ,NY,NX),NZ=1,NP(NY,NX)) 
C     WRITE(*,3114)'TKQC',I,J,NFZ,M,NX,NY,NUM(NY,NX)
C    5,(TKC(NZ,NY,NX),NZ=1,NP(NY,NX))
C    5,(TKD(NZ,NY,NX),NZ=1,NP(NY,NX))
C    5,(TKQC(NZ,NY,NX),NZ=1,NP(NY,NX))
C    5,(TKQD(NZ,NY,NX),NZ=1,NP(NY,NX))
C    6,(FLAIP(NZ,NY,NX),NZ=1,NP(NY,NX))
C    6,(FLAIQ(NZ,NY,NX),NZ=1,NP(NY,NX))
C    7,(RAGC(M,NZ,NY,NX),NZ=1,NP(NY,NX))
C    7,(RAGD(M,NZ,NY,NX),NZ=1,NP(NY,NX))
3114  FORMAT(A8,7I4,80E12.4)
C     ENDIF
      ENDIF
C
C     BOUNDARY LAYER CONDUCTANCES FOR EXPORT TO ‘TRNSFR.F’
C
C     PARR,PARG=boundary layer conductances above litter,soil 
C        surfaces used in ‘trnsfr.f’ (m3 t-1)
C     PAREX=term used to calculate boundary layer
C        conductance for latent heat from ‘hour1.f’ (m2 h t-1)
C     RATG=aerodynamic+canopy boundary layer resistance (h m-1)
C     RARH=surface litter boundary layer resistance (h m-1)
C     RAGS=soil surface boundary layer resistance (h m-1)
C     RAS=snowpack aerodynamic resistance (h m-1)
C
      PARR(M,NY,NX)=PAREX(NY,NX)/(RATG+RARH+RAS(NY,NX))
      PARG(M,NY,NX)=PAREX(NY,NX)/(RATG+RAGS+RAS(NY,NX)) 
C     WRITE(*,7748)'PARR',I,J,NFZ,M,NY,NX
C    2,PARR(M,NY,NX),PAREX(NY,NX),RATG,RARH,RAGS,RAS(NY,NX)
7748  FORMAT(A8,6I4,12E12.4)
C
C     WATER AND ENERGY TRANSFER THROUGH SOIL PROFILE
C
C     N3,N2,N1=L,NY,NX of source grid cell
C     N6,N5,N4=L,NY,NX of destination grid cell
C     
      IFLGH=0
      N1=NX
      N2=NY
      DO 4400 L=1,NL(NY,NX)
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
C
C     END ARTIFICIAL SOIL WARMING PREVENT LATERAL FLOW
C
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
C     FIND NEXT EXISTING DESTINATION SOIL LAYER
C
C     VOLX=soil volume (m3)
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
C     VHCP1=soil layer heat capacity (MJ K-1)
C     TK1N3,TK1N6=soil temperature in source,destination cells (K)
C     HFLXF=heat released by combustion in previous time step (MJ t-1)
C     THETA1,THETAL=micropore water concentration 
C        in source,destination cells (m3 m-3)
C     THETZ,POROS=minimum,saturated soil water concentration 
C        from ‘hour1.f’ (m3 m-3)
C
C     ADD HEAT FROM COMBUSTION
C
      IF(VHCP1(N3,N2,N1).GT.VHCPRX(N2,N1))THEN
      TK1N3=TK1(N3,N2,N1)+HFLXF(N3,N2,N1)/VHCP1(N3,N2,N1)
      ELSE
      TK1N3=TK1(N3,N2,N1)
      ENDIF
      IF(VHCP1(N6,N5,N4).GT.VHCPRX(N5,N4))THEN
      TK1N6=TK1(N6,N5,N4)+HFLXF(N6,N5,N4)/VHCP1(N6,N5,N4)
      ELSE
      TK1N6=TK1(N6,N5,N4)
      ENDIF
      IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4)
     2.AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN
      THETA1=AMAX1(THETZ(N3,N2,N1),AMIN1(POROS(N3,N2,N1)
     2,VOLW1(N3,N2,N1)/VOLY(N3,N2,N1)))
      THETAL=AMAX1(THETZ(N6,N5,N4),AMIN1(POROS(N6,N5,N4)
     2,VOLW1(N6,N5,N4)/VOLY(N6,N5,N4)))
C
C     WATER POTENTIAL OF UPPER SOIL LAYER
C
C     BKVL=soil mass (0=pond)  (Mg m-3)
C     FC,WP=water contents at field capacity,wilting point from soil
C        file (m3 m-3)
C     FCL,WPL=log FC,WP
C     FCD,PSD=FCL-WPL,log(POROS)-FCL
C     PSISA1,PSISX,PSISE=soil matric,minimum,saturation potential (MPa)
C     PSIMX,PSIMD,PSIMS=log water potential at FC,WP,saturation (MPa)
C     PSISD=PSIMX-PSIMS
C     SRP=parameter for deviation from linear log-log water retention
C
      IF(BKVL(N3,N2,N1).GT.ZEROS(NY,NX))THEN
      IF(THETA1.LT.FC(N3,N2,N1))THEN
      PSISA1(N3,N2,N1)=AMAX1(PSISX,-EXP(PSIMX(N2,N1)
     2+((FCL(N3,N2,N1)-LOG(THETA1))
     3/FCD(N3,N2,N1)*PSIMD(N2,N1))))
      ELSEIF(THETA1.LT.POROS(N3,N2,N1)-DTHETW)THEN 
      PSISA1(N3,N2,N1)=-EXP(PSIMS(N2,N1)
     2+((AMAX1(0.0,(PSL(N3,N2,N1)-LOG(THETA1)))
     3/PSD(N3,N2,N1))**SRP(N3,N2,N1)*PSISD(N2,N1)))
      ELSE
      PSISA1(N3,N2,N1)=PSISE(N3,N2,N1) 
      ENDIF
C
C     WATER POTENTIAL OF UPPER POND LAYER
C
C     THETIX,THETWX=pond ice,water concentration (m3 m-3)
C     FCI,WPI=ice field capacity,wilting point (m3 m-3)
C     PSISA1,PSISX,PSISE=soil matric,minimum,saturation potential (MPa)
C
      ELSEIF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      FCX=FCI*THETIX(N3,N2,N1)
      WPX=WPI*THETIX(N3,N2,N1)
      FCLX=LOG(FCX)
      WPLX=LOG(WPX)
      PSDX=PSL(N3,N2,N1)-FCLX
      FCDX=FCLX-WPLX
      IF(THETWX(N3,N2,N1).LT.FCX)THEN
      PSISA1(N3,N2,N1)=AMAX1(PSISX,-EXP(PSIMX(N2,N1)
     2+((FCLX-LOG(THETWX(N3,N2,N1)))
     3/FCDX*PSIMD(NY,NX))))
      ELSEIF(THETWX(N3,N2,N1).LT.POROS(N3,N2,N1)-DTHETW)THEN 
      PSISA1(N3,N2,N1)=-EXP(PSIMS(N2,N1)
     2+((AMAX1(0.0,(PSL(N3,N2,N1)-LOG(THETWR)))
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
C     WATER POTENTIAL OF LOWER SOIL LAYER
C
C     BKVL=soil mass (0=pond) (Mg m-3)
C     FC,WP=water contents at field capacity,wilting point from soil
C        file (m3 m-3)
C     FCL,WPL=log FC,WP
C     FCD,PSD=FCL-WPL,log(POROS)-FCL
C     PSISA1,PSISX,PSISE=soil matric,minimum,saturation potential (MPa)
C     PSIMX,PSIMD,PSIMS=log water potential at FC,WP,saturation (MPa)
C     PSISD=PSIMX-PSIMS
C     SRP=parameter for deviation from linear log-log water retention
C
      IF(BKVL(N6,N5,N4).GT.ZEROS(NY,NX))THEN
      IF(THETAL.LT.FC(N6,N5,N4))THEN
      PSISA1(N6,N5,N4)=AMAX1(PSISX,-EXP(PSIMX(N5,N4)
     2+((FCL(N6,N5,N4)-LOG(THETAL))
     3/FCD(N6,N5,N4)*PSIMD(N5,N4))))
      ELSEIF(THETAL.LT.POROS(N6,N5,N4)-DTHETW)THEN 
      PSISA1(N6,N5,N4)=-EXP(PSIMS(N5,N4)
     2+((AMAX1(0.0,(PSL(N6,N5,N4)-LOG(THETAL)))
     3/PSD(N6,N5,N4))**SRP(N6,N5,N4)*PSISD(N5,N4)))
      ELSE
      PSISA1(N6,N5,N4)=PSISE(N6,N5,N4) 
      ENDIF
C
C     WATER POTENTIAL OF LOWER POND LAYER
C
C     THETIX,THETWX=pond ice,water concentration (m3 m-3)
C     FCI,WPI=ice field capacity,wilting point (m3 m-3)
C     PSISA1,PSISX,PSISE=soil matric,minimum,saturation potential (MPa)
C
      ELSEIF(VOLX(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      FCX=FCI*THETIX(N6,N5,N4)
      WPX=WPI*THETIX(N6,N5,N4)
      FCLX=LOG(FCX)
      WPLX=LOG(WPX)
      PSDX=PSL(N6,N5,N4)-FCLX
      FCDX=FCLX-WPLX
      IF(THETWX(N6,N5,N4).LT.FCX)THEN
      PSISA1(N6,N5,N4)=AMAX1(PSISX,-EXP(PSIMX(N5,N4)
     2+((FCLX-LOG(THETWX(N6,N5,N4)))
     3/FCDX*PSIMD(NY,NX))))
      ELSEIF(THETWX(N6,N5,N4).LT.POROS(N6,N5,N4)-DTHETW)THEN 
      PSISA1(N6,N5,N4)=-EXP(PSIMS(NY,NX)
     2+((AMAX1(0.0,(PSL(N6,N5,N4)-LOG(THETWX(N6,N5,N4))))
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
C        (m3 m-3)
C     CND1,CNDL=hydraulic conductivities in source,destination cells
C        (m2 h MPa-1) 
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
C     PSISM1=soil matric potential (MPa)
C     VOLWX1=VOLW1 accounting for wetting front (m3)
C
C     DARCY FLOW IF BOTH CELLS ARE SATURATED
C     (CURRENT WATER POTENTIAL > AIR ENTRY WATER POTENTIAL)
C
C     PSISA1,PSISA=soil matric,air entry potential from ‘hour1.f’(MPa)
C     THETW1,THETWL=soil water concentration in 
C        source,destination layers (m3 m-3)
C     K1,KL=pore water class for calculating hydraulic conductivity
C     PSISM1=PSISA1 accounting for PSISA (MPa)
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
C     PSISA1,PSISA=soil matric,air entry potential from ‘hour1.f’(MPa)
C     THETW1,THETWL=soil water concentration in 
C        source,destination layers (m3 m-3)
C     K1,KL=pore water class for calculating hydraulic conductivity
C     PSISM1=PSISA1 accounting for PSISA (MPa)
C
      ELSEIF(PSISA1(N3,N2,N1).GT.PSISA(N3,N2,N1))THEN
      THETW1=THETA1
      THETWL=AMAX1(THETZ(N6,N5,N4),AMIN1(POROS(N6,N5,N4)
     2,VOLWX1(N6,N5,N4)/VOLY(N6,N5,N4)))
      K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)
     2-THETW1)/POROS(N3,N2,N1))+1))
      KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4)
     2-AMIN1(THETS(N6,N5,N4),THETWL))/POROS(N6,N5,N4))+1))
      PSISM1(N3,N2,N1)=PSISA1(N3,N2,N1)
C
C     IF DESTINATION CELL IS SOIL
C
      IF(BKVL(N6,N5,N4).GT.ZEROS(NY,NX))THEN
      IF(THETWL.LT.FC(N6,N5,N4))THEN
      PSISM1(N6,N5,N4)=AMAX1(PSISX,-EXP(PSIMX(N5,N4)
     2+((FCL(N6,N5,N4)-LOG(THETWL))
     3/FCD(N6,N5,N4)*PSIMD(N5,N4))))
      ELSEIF(THETWL.LT.POROS(N6,N5,N4)-DTHETW)THEN 
      PSISM1(N6,N5,N4)=-EXP(PSIMS(N5,N4)
     2+((AMAX1(0.0,(PSL(N6,N5,N4)-LOG(THETWL)))
     3/PSD(N6,N5,N4))**SRP(N6,N5,N4)*PSISD(N5,N4)))
      ELSE
      THETWL=POROS(N6,N5,N4)
      PSISM1(N6,N5,N4)=PSISE(N6,N5,N4)
      ENDIF
C
C     IF DESTINATION CELL IS POND
C
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
C     PSISA1,PSISA=soil matric,air entry potential from ‘hour1.f’(MPa)
C     THETW1,THETWL=soil water concentration in 
C        source,destination layers (m3 m-3)
C     K1,KL=pore water class for calculating hydraulic conductivity
C     PSISM1=PSISA1 accounting for PSISA (MPa)
C
      ELSEIF(PSISA1(N6,N5,N4).GT.PSISA(N6,N5,N4))THEN
      THETW1=AMAX1(THETZ(N3,N2,N1),AMIN1(POROS(N3,N2,N1)
     2,VOLWX1(N3,N2,N1)/VOLY(N3,N2,N1)))
      THETWL=THETAL
      K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)
     2-AMIN1(THETS(N3,N2,N1),THETW1))/POROS(N3,N2,N1))+1))
      KL=MAX(1,MIN(100,INT(100.0*(POROS(N6,N5,N4)
     2-THETWL)/POROS(N6,N5,N4))+1))
C
C     IF DESTINATION CELL IS SOIL
C
      IF(BKVL(N3,N2,N1).GT.ZEROS(NY,NX))THEN
      IF(THETW1.LT.FC(N3,N2,N1))THEN
      PSISM1(N3,N2,N1)=AMAX1(PSISX,-EXP(PSIMX(N2,N1)
     2+((FCL(N3,N2,N1)-LOG(THETW1))
     3/FCD(N3,N2,N1)*PSIMD(N2,N1))))
      ELSEIF(THETW1.LT.POROS(N3,N2,N1)-DTHETW)THEN 
      PSISM1(N3,N2,N1)=-EXP(PSIMS(N2,N1)
     2+((AMAX1(0.0,(PSL(N3,N2,N1)-LOG(THETW1)))
     3/PSD(N3,N2,N1))**SRP(N3,N2,N1)*PSISD(N2,N1)))
      ELSE
      THETW1=POROS(N3,N2,N1)
      PSISM1(N3,N2,N1)=PSISE(N3,N2,N1) 
      ENDIF
C
C     IF DESTINATION CELL IS POND
C
      ELSE
      THETW1=POROS(N3,N2,N1)
      PSISM1(N3,N2,N1)=PSISE(N3,N2,N1) 
      ENDIF
C
C     RICHARDS FLOW IF NEITHER CELL IS SATURATED
C     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POTENTIAL)
C
C     THETW1,THETWL=soil water concentration in 
C        source,destination layers (m3 m-3)
C     K1,KL=pore water class for calculating hydraulic conductivity
C     PSISM1=PSISA1 accounting for PSISA (MPa)
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
C        (m2 h MPa-1) 
C     HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity 
C        (m2 h MPa-1) 
C     K1,KL=pore water class for calculating hydraulic conductivity
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
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
C        (MPa)
C
      PSIST1=PSISM1(N3,N2,N1)+PSISH(N3,N2,N1)+ORFLN*PSISO(N3,N2,N1)
      PSISTL=PSISM1(N6,N5,N4)+PSISH(N6,N5,N4)+ORFLN*PSISO(N6,N5,N4)
C     IF(N6.EQ.12)THEN
C     WRITE(*,7272)'PSIM',I,J,M,N1,N2,N3,N4,N5,N6
C    2,PSISM1(N3,N2,N1),PSISM1(N6,N5,N4),PSISA1(N3,N2,N1)
C    3,PSISA1(N6,N5,N4),PSISO(N3,N2,N1),PSISO(N6,N5,N4)
C    2,PSIMX(N5,N4),FCL(N6,N5,N4),FCD(N6,N5,N4),PSIMD(N5,N4)
C    3,POROS(N6,N5,N4),PSIMS(N5,N4),PSL(N6,N5,N4),PSD(N6,N5,N4)
C    4,SRP(N6,N5,N4),PSISD(N5,N4),PSISE(N6,N5,N4)
C    5,THETZ(N6,N5,N4),POROS(N6,N5,N4),VOLW1(N6,N5,N4),VOLX(N6,N5,N4)
7272  FORMAT(A8,9I4,30E12.4)
C     ENDIF
C
C     HYDRAULIC CONDUCTIVITY FROM CURRENT WATER CONTENT
C
C     CND1,CNDL=hydraulic conductivities in source,destination cells
C        (m2 h MPa-1) 
C     AVCNDL=source-destination hydraulic conductance
C        (m h MPa-1) 
C     DLYR=layer thickness (m)
C
C     IF(CND1.GT.0.0.AND.CNDL.GT.0.0)THEN
      AVCNDL=2.0*CND1*CNDL/(CND1*DLYR(N,N6,N5,N4)
     2+CNDL*DLYR(N,N3,N2,N1)) 
C     ELSE
C     AVCNDL=0.0
C     ENDIF
C
C     WATER FLUX FROM WATER POTENTIALS, HYDRAULIC CONDUCTIVITY
C     LIMITED BY WATER POTENTIAL GRADIENT, COUPLED WITH
C     CONVECTIVE HEAT FLUX FROM WATER FLUX
C
C     FLQX=micropore water flux unlimited by source water (m3 t-1)
C     PSIST1,PSISTL=total water potentials of source,destination 
C        micropores (MPa) 
C     FLQZ=FLQX + saturated soil water flux (m3 t-1)
C     THETS=water concentration at air entry potential from ‘hour1.f’ 
C        (m3 m-3)
C     VOLW2,VOLP1=water,air contents of source,destination micropores 
C        (m3)
C     VOLY=volume of soil micropore fraction (m3)
C     HWFLWL=total convective heat from micropore water+vapor flux 
C        (MJ t-1) 
C     VOLP1Z=excess water+ice relative to porosity (m3)
C     FLQL,FLQ2=soil water flux limited by soil water content (m3 t-1)
C     XNPHX=time step of flux calculations from ‘wthr.f’(h t-1)
C
      FLQX=AVCNDL*(PSIST1-PSISTL)*AREA(N,N3,N2,N1)*XNPHX
      IF(FLQX.GE.0.0)THEN
      IF(THETW1.GT.THETS(N3,N2,N1))THEN
      FLQZ=FLQX+AMIN1((THETW1-THETS(N3,N2,N1))
     2*VOLY(N3,N2,N1),AMAX1(0.0,(THETS(N6,N5,N4)-THETWL)
     3*VOLY(N6,N5,N4)))*XNPXX
      ELSE
      FLQZ=FLQX
      ENDIF
      VOLW2N=VOLW2(N3,N2,N1)+FLWL(3,N3,N2,N1)
      FLQL=AMAX1(0.0,AMIN1(FLQZ,VOLW2N*XNPXX
     2,VOLP1(N6,N5,N4)*XNPXX))
      FLQ2=AMAX1(0.0,AMIN1(FLQX,VOLW2N*XNPXX
     2,VOLP1(N6,N5,N4)*XNPXX))
      ELSE
      IF(THETWL.GT.THETS(N6,N5,N4))THEN
      FLQZ=FLQX+AMAX1((THETS(N6,N5,N4)-THETWL) 
     2*VOLY(N6,N5,N4),AMIN1(0.0,(THETW1-THETS(N3,N2,N1))
     3*VOLY(N3,N2,N1)))*XNPXX
      ELSE
      FLQZ=FLQX
      ENDIF
      VOLP2N=VOLP1(N3,N2,N1)-FLWL(3,N3,N2,N1)
      FLQL=AMIN1(0.0,AMAX1(FLQZ,-VOLW2(N6,N5,N4)*XNPXX
     2,-VOLP2N*XNPXX))
      FLQ2=AMIN1(0.0,AMAX1(FLQX,-VOLW2(N6,N5,N4)*XNPXX
     2,-VOLP2N*XNPXX))
      ENDIF
C
C     ACCOUNT FOR EXCESS WATER+ICE VOLUME DURING FREEZING
C
C     VOLWP1Z=excess water+ice (-ve) (m3)
C     FLQL,FLQ2=soil water flux limited by soil water content (m3 t-1)
C     HWFLQL=convective heat from soil water flux (MJ t-1)
C
      IF(N.EQ.3.AND.VOLP1Z(N6,N5,N4).LT.0.0)THEN
      FLQL=FLQL+AMIN1(0.0,AMAX1(-VOLW2(N6,N5,N4)*XNPXX
     2,VOLP1Z(N6,N5,N4)))
      FLQ2=FLQ2+AMIN1(0.0,AMAX1(-VOLW2(N6,N5,N4)*XNPXX
     2,VOLP1Z(N6,N5,N4)))
      ENDIF
      IF(FLQL.GT.0.0)THEN
      HWFLQL=4.19*TK1N3*FLQL
      ELSE
      HWFLQL=4.19*TK1N6*FLQL
      ENDIF
C
C     UPDATE SOIL WATER,AIR CONTENTS
C
C     VOLW2,VOLI2,VOLP2=soil water,ice,air contents (m3)
C     FLQL=soil water flux limited by soil water content (m3 t-1)
C
      VOLW2(N3,N2,N1)=VOLW2(N3,N2,N1)-FLQL
      VOLW2(N6,N5,N4)=VOLW2(N6,N5,N4)+FLQL
      VOLP2(N3,N2,N1)=AMAX1(0.0,VOLA1(N3,N2,N1)
     2-VOLW2(N3,N2,N1)-VOLI2(N3,N2,N1))
      VOLP2(N6,N5,N4)=AMAX1(0.0,VOLA1(N6,N5,N4)
     2-VOLW2(N6,N5,N4)-VOLI2(N6,N5,N4))
C
C     MACROPORE FLOW FROM POISEUILLE FLOW IF MACROPORES PRESENT
C
C     VOLAH1,VOLIH1,VOLWH1,VOLPH1=total,ice-,water-,air-filled
C        macropore volume (m3) 
C     PSISH1,PSISHL=macropore total water potential in
C        source,destination (MPa)
C     DLYR=layer thickness (m)
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
C     FLWHX,FLWHL=macropore water flux unlimited,limited 
C        by macropore source water (m3 t-1)
C     AVCNHL=macropore hydraulic conductance (m h MPa-1)
C     PSISH1,PSISHL=macropore total water potential in
C        source,destination (MPa)
C     VOLWH2,VOLPH1=water,air contents of source,destination macropores
C        (m3)
C     XNPHX=time step of flux calculations from ’wthr.f’(h t-1)
C
      FLWHX=AVCNHL(N,N6,N5,N4)*(PSISH1-PSISHL)*AREA(N,N3,N2,N1)*XNPHX
      IF(N.NE.3)THEN
      IF(PSISH1.GT.PSISHL)THEN
      FLWHL(N,N6,N5,N4)=AMAX1(0.0,AMIN1(AMIN1(VOLWH1(N3,N2,N1)
     2,VOLPH1(N6,N5,N4))*XNPXX,FLWHX))
      ELSEIF(PSISH1.LT.PSISHL)THEN
      FLWHL(N,N6,N5,N4)=AMIN1(0.0,AMAX1(AMAX1(-VOLWH1(N6,N5,N4)
     2,-VOLPH1(N3,N2,N1))*XNPXX,FLWHX))
      ELSE
      FLWHL(N,N6,N5,N4)=0.0
      ENDIF
      ELSE
      FLWHL(N,N6,N5,N4)=AMAX1(0.0,AMIN1(AMIN1(VOLWH1(N3,N2,N1)*XNPXX 
     2+FLWHL(N,N3,N2,N1),VOLPH1(N6,N5,N4)*XNPXX),FLWHX))
      ENDIF
C
C     ACCOUNT FOR EXCESS WATER+ICE DURING FREEZING
C
C     VOLPH1Z=excess water+ice in macropores during freezing (–ve)(m3)
C     FLWHL=macropore water flux limited 
C        by macropore source water (m3 t-1)
C     VOLWH2,VOLPH1=water,air contents of source,destination macropores
C        (m3)
C     HWFLHL=convective heat flux from macropore water flux (MJ t-1)
C
      IF(N.EQ.3)THEN
      FLWHL(N,N6,N5,N4)=FLWHL(N,N6,N5,N4)
     2+AMIN1(0.0,AMAX1(-VOLWH1(N6,N5,N4)*XNPXX,VOLPH1Z(N6,N5,N4)))
      ENDIF
      FLWHM(M,N,N6,N5,N4)=FLWHL(N,N6,N5,N4)
C     IF(I.EQ.70)THEN
C     WRITE(*,5478)'FLWH',I,J,M,N1,N2,N3,IFLGH
C    2,FLHM,FLWHX,FLWHL(N,N3,N2,N1),FLWHL(N,N6,N5,N4) 
C    2,AVCNHL(N,N6,N5,N4),PSISH(N3,N2,N1),PSISH(N6,N5,N4) 
C    3,VOLPH1(N3,N2,N1),VOLPH1(N6,N5,N4),VOLWH1(N3,N2,N1)
C    4,VOLWH1(N6,N5,N4),VOLAH1(N3,N2,N1),VOLAH1(N6,N5,N4)
C    5,DLYR(N,N6,N5,N4),DLYR(N,N3,N2,N1),AREA(N,N3,N2,N1)
C    7,CNDH1(N3,N2,N1),CNDH1(N6,N5,N4),XNPH,HWFLHL 
C    8,VOLPH1Z(N6,N5,N4)
5478  FORMAT(A8,7I4,30E12.4)
C     ENDIF
      ELSE
      FLWHL(N,N6,N5,N4)=0.0
      FLWHM(M,N,N6,N5,N4)=0.0
      IF(VOLPH1(N6,N5,N4).LE.0.0)IFLGH=1
      ENDIF
      IF(FLWHL(N,N6,N5,N4).GT.0.0)THEN
      HWFLHL=4.19*TK1N3*FLWHL(N,N6,N5,N4)
      ELSE
      HWFLHL=4.19*TK1N6*FLWHL(N,N6,N5,N4)
      ENDIF
C
C     VAPOR CONCENTRATION AND DIFFUSIVITY IN EACH GRID CELL
C
C     VOLPM,VOLV2=soil air,vapor volume (m3)
C     VP1,VPL=vapor concentration in source,destination (m3 m-3)
C     CNV1,CNV2=vapor conductivities of source, destination (m2 h-1) 
C     WGSGL=vapor diffusivity from ‘hour1.f’ (m2 t-1)
C     POROS,POROQ=porosity (m3 m-3), tortuosity from ‘starts.f’ 
C     AVCNVL=source,destination vapor conductance (m t-1)
C     DLYR=soil layer depth (m)
C
      IF(VOLPM(M,N3,N2,N1).GT.ZEROS(N2,N1)
     2.AND.VOLPM(M,N6,N5,N4).GT.ZEROS(N2,N1))THEN
      VP1=AMAX1(0.0,VOLV2(N3,N2,N1)/VOLPM(M,N3,N2,N1))
      VPL=AMAX1(0.0,VOLV2(N6,N5,N4)/VOLPM(M,N6,N5,N4))
      CNV1=WGSGL(N3,N2,N1)*POROQ
     2*THETPM(M,N3,N2,N1)**2/POROS(N3,N2,N1) 
      CNVL=WGSGL(N6,N5,N4)*POROQ
     2*THETPM(M,N6,N5,N4)**2/POROS(N6,N5,N4) 
      AVCNVL=2.0*CNV1*CNVL
     2/(CNV1*DLYR(N,N6,N5,N4)+CNVL*DLYR(N,N3,N2,N1))
C
C     VAPOR FLUX FROM VAPOR CONCENTRATION AND DIFFUSIVITY,
C     AND CONVECTIVE HEAT FLUX FROM VAPOR FLUX
C
C     FLVC,FLVQ=vapor flux unlimited,limited by vapor content (m3 t-1)
C     VOLV2=soil vapor content (m3)
C     HWFLVQ=convective heat of vapor flux (MJ t-1)
C     TK1N3,TK1N6=soil temperature (K)
C     XNPHX=time step for flux calculations from ‘wthr.f’(h t-1)
C
      FLVC=AVCNVL*(VP1-VPL)*AREA(N,N3,N2,N1)*XNPHX
      IF(FLVC.GE.0.0)THEN
      FLVQ=AMAX1(0.0,AMIN1(FLVC,VOLV2(N3,N2,N1)*XNPXX)) 
      HWFLVQ=4.19*TK1N3*FLVQ
      ELSE
      FLVQ=AMIN1(0.0,AMAX1(FLVC,-VOLV2(N6,N5,N4)*XNPXX))
      HWFLVQ=4.19*TK1N6*FLVQ
      ENDIF
      VOLV2(N3,N2,N1)=VOLV2(N3,N2,N1)-FLVQ
      VOLV2(N6,N5,N4)=VOLV2(N6,N5,N4)+FLVQ
      ELSE
      FLVC=0.0
      FLVQ=0.0
      HWFLVQ=0.0
      VP1=0.0
      VPL=0.0
      AVCNVL=0.0
      ENDIF
C
C     FLWL=total water+vapor flux between source,destination (m3 t-1)
C     FLWLX=total unsaturated water+vapor flux between source,
C        destination (m3 t-1)
C     HWFLWL=total convective heat from micropore water+vapor flux 
C        (MJ t-1) 
C
      FLWL(N,N6,N5,N4)=FLQL 
      FLVL(N,N6,N5,N4)=FLVQ
      FLWLX(N,N6,N5,N4)=FLQ2 
      HWFLWL=HWFLQL+HWFLVQ
C     IF(N3.EQ.1)THEN
C     WRITE(*,1115)'FLWL',I,J,NFZ,N1,N2,N3,N4,N5,N6,M,N,K1,KL
C    2,FLWL(N,N3,N2,N1),FLWL(N,N6,N5,N4)
C    3,FLQX,FLQZ,FLQL,FLQ2,AVCNDL,PSIST1,PSISTL,CND1,CNDL,FKSAT
C    2,FLVL(N,N3,N2,N1),FLVL(N,N6,N5,N4)
C    3,AVCNVL,THETPM(M,N3,N2,N1),THETPM(M,N6,N5,N4)
C    4,VP1,VPL,VPQG(M,NY,NX),FLVC,FLVQ
C    3,VOLV2(N3,N2,N1),VOLV2(N6,N5,N4),XNPHX
C    4,VOLV1(N3,N2,N1),VOLV1(N6,N5,N4)
C    6,THETW1,THETS(N3,N2,N1),THETWL,THETS(N6,N5,N4)
C    6,THETZ(N3,N2,N1),THETZ(N6,N5,N4),THETA1,THETAL  
C    6,POROS(N3,N2,N1),POROS(N6,N5,N4)  
C    7,PSISA1(N3,N2,N1),PSISA1(N6,N5,N4),PSISM1(N3,N2,N1)
C    7,PSISM1(N6,N5,N4),PSISH(N3,N2,N1),PSISH(N6,N5,N4)
C    3,VOLY(N3,N2,N1),VOLY(N6,N5,N4),XNPXX 
C    3,FLVQ,VP1,VPL,CNV1,CNVL,AVCNVL 
C    4,VOLA1(N6,N5,N4),VOLW1(N6,N5,N4),VOLI1(N6,N5,N4) 
C    5,SCNV(N3,N2,N1),SCNV(N6,N5,N4),VOLP1(N3,N2,N1) 
C    7,VOLP1(N6,N5,N4),VOLT(N3,N2,N1),VOLT(N6,N5,N4)
C    8,DLYR(N,N3,N2,N1),DLYR(N,N6,N5,N4),AREA(N,N3,N2,N1)
C    3,HWFLWL,HWFLQL,HWFLVQ
C    9,TK1N3,TK1N6,TKY 
C    9,VHCP1(N3,N2,N1),VHCP1(N6,N5,N4),POROS(N6,N5,N4)
C    9,VOLP1(N3,N2,N1),VOLX(N3,N2,N1),VOLA1(N6,N5,N4),VOLA1(N3,N2,N1)
C    8,VOLP1(N6,N5,N4),VOLX(N6,N5,N4),FLW(N,N3,N2,N1),FLW(N,N6,N5,N4)
C    1,THETPM(M,N3,N2,N1),THETPM(M,N6,N5,N4),THETX 
C    4,FLQL1,FLQL2,FLQL3,FLQL4
1115  FORMAT(A8,13I4,100F16.8)
C     ENDIF
C
C     THERMAL CONDUCTIVITY IN EACH GRID CELL
C
C     DTH*,RYL*,XNU*,TRB*=turbulence effects on thermal conductivity
C        of air (*A) and water (*W)
C     THETWX,THETPX=water,air concentration (m3 m-3)
C     TCNDW*,TCNDA*=thermal conductivity of water,air (m MJ h-1 K-1)
C     TCND1,TCNDL=soil thermal conductivity in source,destination 
C        (m MJ h-1 K-1)
C     WTHET*=multiplier for air concentration in thermal conductivity
C     ATCNDL=source-destination thermal conductance (MJ h-1 K-1)
C     DLYR=layer thickness (m)
C
      DTKX=ABS(TK1N3-TK1N6)*1.0E-06
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
C     VHCP1,VHCPRX=current,minimum soil volumetric heat capacity 
C        (MJ K-1)
C     TK1X,TKLX=interim temperatures of source,destination (K)
C     HWFLVQ=convective heat from soil vapor flux (MJ t-1)
C     HFLXG=storage heat flux from snowpack (MJ t-1)   
C     TKY=equilibrium source-destination temperature (K)     
C     HFLWX,HFLWC=source-destination conductive heat flux
C        limited,unlimited by temperature (MJ t-1)
C     ATCNDL=source-destination thermal conductance (MJ h-1 K-1)
C     HFLWSX=source-destination conductive heat flux (MJ t-1)
C     HWFLWL,HWFLHL=convective heat flux through micro,-macro-pores 
C        (MJ t-1) 
C     HFLWL=total conductive+convective source-destination heat flux
C       (MJ t-1)
C
      IF(VHCP1(N3,N2,N1).GT.VHCPRX(N2,N1))THEN
      IF(N3.EQ.NUM(NY,NX).AND.VHCPW(1,N2,N1).LE.VHCPWX(N2,N1))THEN
      TK1X=TK1N3-(HWFLVQ-HFLXG)/VHCP1(N3,N2,N1)
      ELSE
      TK1X=TK1N3-HWFLVQ/VHCP1(N3,N2,N1)
      ENDIF
      ELSE
      TK1X=TK1N3
      ENDIF
      IF(VHCP1(N6,N5,N4).GT.VHCPRX(N5,N4))THEN 
      TKLX=TK1N6+HWFLVQ/VHCP1(N6,N5,N4) 
      ELSE
      TKLX=TK1N6
      ENDIF
      TKY=(VHCP1(N3,N2,N1)*TK1X+VHCP1(N6,N5,N4)*TKLX)
     2/(VHCP1(N3,N2,N1)+VHCP1(N6,N5,N4))
      HFLWX=(TK1X-TKY)*VHCP1(N3,N2,N1) 
      HFLWC=ATCNDL*(TK1X-TKLX)*AREA(N,N3,N2,N1)*XNPHX
      IF(HFLWC.GE.0.0)THEN
      HFLWSX=AMAX1(0.0,AMIN1(HFLWX,HFLWC))
      ELSE
      HFLWSX=AMIN1(0.0,AMAX1(HFLWX,HFLWC))
      ENDIF
      HFLWL(N,N6,N5,N4)=HWFLWL+HWFLHL+HFLWSX
C     IF(I.GT.270.AND.N3.EQ.1)THEN
C     WRITE(*,8765)'HFLWL',I,J,NFZ,M,N1,N2,N3,N4,N5,N6,N
C    2,HFLWL(N,N3,N2,N1),HFLWL(N,N6,N5,N4)
C    2,HWFLWL,HWFLQL,HWFLVQ,HWFLHL,HFLWC,HFLWX,HFLWSX,HFLXG 
C    3,ATCNDL,TK1X,TKLX,TKY,HWFLVQ,TK1N3,TK1N6
C    2,TCND1,TCND2,DLYR(N,N3,N2,N1),DLYR(N,N6,N5,N4)
C    4,VHCP1(N3,N2,N1),VHCP1(N6,N5,N4),VOLY(N3,N2,N1)
C    5,STC(N3,N2,N1),THETWX(N3,N2,N1),TCNDW1
C    2,THETIX(N3,N2,N1) 
C    3,WTHET1,THETPX(N3,N2,N1),TCNDA1
C    4,DTC(N3,N2,N1) 
C    5,VOLY(N6,N5,N4),VOLW1(N3,N2,N1),VOLW1(N6,N5,N4) 
C    3,THETPX(N3,N2,N1),THETIX(N3,N2,N1),THETWX(N3,N2,N1)
C    3,THETPX(N6,N5,N4),THETIX(N6,N5,N4),THETWX(N6,N5,N4)
C    3,THETPY(N3,N2,N1),THETPY(N6,N5,N4)
C    4,RYLNA2,XNUSA2,XNUSW2 
C    4,STC(N6,N5,N4),DTC(N6,N5,N4) 
C    2,WTHET2,TCNDA2,TCNDW2
8765  FORMAT(A8,11I4,60E12.4)
C     ENDIF
C
C     TOTAL WATER, VAPOR AND HEAT FLUXES
C
C     FLW,FLV,FLWX,FLWH=total water,vapor flux through
C        micropores,macropores (m3 t-1) 
C     HFLW=total heat flux (MJ t-1)
C     FLWM=water flux used for solute flux calculations in ‘trnsfr.f’
C        (m3 t-1)
C
      FLW(N,N6,N5,N4)=FLW(N,N6,N5,N4)+FLWL(N,N6,N5,N4)
      FLV(N,N6,N5,N4)=FLV(N,N6,N5,N4)+FLVL(N,N6,N5,N4)
      FLWX(N,N6,N5,N4)=FLWX(N,N6,N5,N4)+FLWLX(N,N6,N5,N4)
      FLWH(N,N6,N5,N4)=FLWH(N,N6,N5,N4)+FLWHL(N,N6,N5,N4)
      HFLW(N,N6,N5,N4)=HFLW(N,N6,N5,N4)+HFLWL(N,N6,N5,N4)
      FLWM(M,N,N6,N5,N4)=FLWL(N,N6,N5,N4)
C     IF(I.EQ.55)THEN
C     WRITE(*,1114)'FLWL2',I,J,NFZ,N4,N5,N6,M,N,FLWL(N,N3,N2,N1)
C    2,FLWL(N,N6,N5,N4),FLW(N,N3,N2,N1),FLW(N,N6,N5,N4)
C    3,FLQL,FLVQ,FLQX,HFLWX 
C    3,CND1,CNDL,AVCNDL,AVCNVL,VP1,VPL,PSIST1,PSISTL 
C    4,UAG,VOLA1(N6,N5,N4),VOLI1(N6,N5,N4),SCNV(N6,N5,N4) 
C    5,VOLP1(N3,N2,N1),VOLP1(N6,N5,N4),TKY
C    7,TK1N3,TK1N6,VOLT(N3,N2,N1),VOLT(N6,N5,N4)
C    8,VOLW1(N6,N5,N4),VOLP1(N6,N5,N4),VOLX(N6,N5,N4),VOLW1(N3,N2,N1) 
C    9,VOLP1(N3,N2,N1),VOLX(N3,N2,N1),VOLA1(N6,N5,N4),VOLA1(N3,N2,N1)
C    6,THETW1,THETWL,PSISA1(N3,N2,N1)
C    7,PSISA1(N6,N5,N4),PSISM1(N3,N2,N1)
C    7,PSISM1(N6,N5,N4),PSISH(N3,N2,N1),PSISH(N6,N5,N4)
C    8,DLYR(N,N3,N2,N1),DLYR(N,N6,N5,N4),AREA(N,N3,N2,N1)
C    9,VHCP1(N3,N2,N1),VHCP1(N6,N5,N4),POROS(N6,N5,N4) 
1114  FORMAT(A8,8I4,50F16.8)
C     ENDIF
C
C     WATER FILM THICKNESS FOR CALCULATING GAS EXCHANGE IN TRNSFR.F
C
C     FILM=soil water film thickness for nutrient
C        and gas uptake in ‘nitro.f’ and uptake.f’ (m)
C     PSISA1=soil matric potential (MPa)
C
      IF(N.EQ.3)THEN
      FILM(M,N6,N5,N4)=AMAX1(1.0E-06
     2,0.5*EXP(-13.833-0.857*LOG(-PSISA1(N6,N5,N4))))
      ENDIF
      ELSEIF(N.NE.3)THEN
      FLWL(N,N6,N5,N4)=0.0
      FLVL(N,N6,N5,N4)=0.0
      FLWLX(N,N6,N5,N4)=0.0
      FLWHL(N,N6,N5,N4)=0.0
      HFLWL(N,N6,N5,N4)=0.0
      FLWM(M,N,N6,N5,N4)=0.0
      FLWHM(M,N,N6,N5,N4)=0.0
      ENDIF
      ELSE
      IF(N.EQ.3)THEN
      FLWL(N,N3,N2,N1)=0.0
      FLVL(N,N3,N2,N1)=0.0
      FLWLX(N,N3,N2,N1)=0.0
      FLWHL(N,N3,N2,N1)=0.0
      HFLWL(N,N3,N2,N1)=0.0
      FLWHM(M,N,N3,N2,N1)=0.0
      FLWHM(M,N,N3,N2,N1)=0.0
      ELSE
      FLWL(N,N6,N5,N4)=0.0
      FLVL(N,N6,N5,N4)=0.0
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
C     VOLPD,VOLPHD=air-filled porosity in micropores,macropores
C        used to constrain boundary water flux (m3)
C
      DO 9595 NX=NHW,NHE
      DO 9590 NY=NVN,NVS
      DO 9585 L=NUM(NY,NX),NL(NY,NX)
      VOLPD=VOLP1(L,NY,NX)-FLWL(3,L,NY,NX)+FLWL(3,L+1,NY,NX)
      VOLPXD=VOLPD
      VOLPHD=VOLPH1(L,NY,NX)
C
C     IDENTIFY CONDITIONS FOR MICROPORE DISCHARGE TO WATER TABLE
C
C     IDTBL=water table flag from site file
C     DPTH,DTBLX=depth to layer midpoint,natural water table (m)
C     PSISM1,PSISE=matric,saturation water potential (MPa)
C     DPTHA=active layer depth (m)
C     IFLGU=micropore discharge flag to natural water table
C        :0=discharge enabled,1=discharge disabled
C
      IF(IDTBL(NY,NX).NE.0.AND.DPTH(L,NY,NX).LT.DTBLX(NY,NX))THEN
      IF(PSISM1(L,NY,NX).GE.0.0098*(DPTH(L,NY,NX)-DTBLX(NY,NX)))THEN
      IFLGU=0
      LU=L
      DO 9565 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
C     WRITE(*,9567)'IFLGCHK',I,J,NFZ,M,NX,NY,LL,IFLGU,IFLGUH
C    2,PSISM1(LL,NY,NX),0.0098*(DPTH(LL,NY,NX)-DTBLX(NY,NX))
C    3,DPTH(LL,NY,NX),DTBLX(NY,NX) 
      IF(DPTH(LL,NY,NX).LT.DTBLX(NY,NX))THEN
      LU=LL
      IF((PSISM1(LL,NY,NX).LT.0.0098*(DPTH(LL,NY,NX)-DTBLX(NY,NX)) 
     2.AND.L.NE.NL(NY,NX)).OR.DPTH(LL,NY,NX).GT.DPTHA(NY,NX))THEN
      IFLGU=1
      GO TO 9566
      ENDIF
      ENDIF
9565  CONTINUE
9566  CONTINUE
      ELSE
      IFLGU=1
      ENDIF
      ELSE
      IFLGU=1
      ENDIF
C
C     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO WATER TABLE
C
C     VOLAH1,VOLWH1,VOLIH1=macropore volume,water,ice content (m3)
C     DPTHH depth to layer macropore water (m) 
C     CDPTH=depth to layer bottom (m)
C     DLYR=layer thickness (m)
C     IFLGUH=macropore discharge flag to natural water table
C        :0=discharge enabled,1=discharge disabled
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
      ELSE
      IFLGUH=1
      ENDIF
C     IF(NFZ.EQ.1)THEN
C     WRITE(*,9567)'IFLGU',I,J,NFZ,M,NX,NY,L,IFLGU,IFLGUH
C    2,PSISM1(L,NY,NX),0.0098*(DPTH(L,NY,NX)-DTBLX(NY,NX)) 
C    2,DPTH(L,NY,NX),DTBLX(NY,NX),DTBLZ(NY,NX),DTBLI(NY,NX) 
C    3,VOLAH1(L,NY,NX),VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),CDPTH(L,NY,NX)
C    4,DLYR(3,L,NY,NX),DTBLZ(NY,NX),DPTHH,THETX,DPTHA(NY,NX)
9567  FORMAT(A8,9I4,40E12.4)
C     ENDIF 
C
C     IDENTIFY CONDITIONS FOR MICROPORE DISCHARGE TO TILE DRAIN
C
C     IDTBL=water table flag from site file
C     DPTH,DTBLY=depth to layer midpoint, artificial water table (m)
C     PSISM1,PSISE=soil,saturation matric potential (MPa)
C     DPTHA=active layer depth (m)
C     IFLGD=micropore discharge flag to artificial water table
C        :0=discharge enabled,1=discharge disabled
C
      IF(IDTBL(NY,NX).GE.3.AND.DPTH(L,NY,NX).LT.DTBLY(NY,NX))THEN
      IF(PSISM1(L,NY,NX).GE.0.0098*(DPTH(L,NY,NX)-DTBLY(NY,NX)))THEN
      IFLGD=0
      LD=0
      IF(L.LT.NL(NY,NX))THEN
      DO 9568 LL=L+1,NL(NY,NX)
      IF(DPTH(LL,NY,NX).LT.DTBLY(NY,NX))THEN
      LD=LL
      IF((PSISM1(LL,NY,NX).LT.0.0098*(DPTH(LL,NY,NX)-DTBLY(NY,NX)) 
     2.AND.L.NE.NL(NY,NX)).OR.DPTH(LL,NY,NX).GT.DPTHA(NY,NX))THEN
      IFLGD=1
      GO TO 9569
      ENDIF
      ENDIF
9568  CONTINUE
9569  CONTINUE
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
C     VOLAH1,VOLWH1,VOLIH1=macropore volume,water,ice content (m3)
C     DPTHH depth to layer macropore water (m) 
C     CDPTH=depth to layer bottom (m)
C     DLYR=layer thickness (m)
C     IFLGDH=macropore discharge flag to artificial water table
C        :0=discharge enabled,1=discharge disabled
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
      ELSE
      IFLGDH=1
      ENDIF
C     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,9567)'IFLGD',I,J,NFZ,M,NX,NY,L,IFLGD,IFLGDH
C    2,PSISM1(L,NY,NX),PSISE(L,NY,NX)
C    3+0.0098*(DPTH(L,NY,NX)-DTBLY(NY,NX)) 
C    2,DPTH(L,NY,NX),DTBLY(NY,NX),DTBLD(NY,NX),DTBLDI(NY,NX) 
C    3,VOLAH1(L,NY,NX),VOLWH1(L,NY,NX),VOLIH1(L,NY,NX),CDPTH(L,NY,NX)
C    4,DLYR(3,L,NY,NX),DPTHH,THETX,DPTHA(NY,NX) 
C    5,RCHGNA(NY,NX),RCHGEA(NY,NX),RCHGSA(NY,NX),RCHGWA(NY,NX)
C    6,RCHGNB(NY,NX),RCHGEB(NY,NX),RCHGSB(NY,NX),RCHGWB(NY,NX)
C     ENDIF 
C
C     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
C     ENTERED IN 'READS'
C
C     N3,N2,N1=L,NY,NX of source grid cell
C     M6,M5,M4=L,NY,NX of destination grid cell
C     N6,N5,N4=L,NY,NX of destination grid cell E OR S
C     N6,N5B,N4B=L,NY,NX of destination grid cell W or N
C     XN=direction parameter
C     RCHQE*,RCHQW*,RCHQN*,RCHQS*=flux conditions at E,W,N,S
C        lateral boundaries from site or soil management file
C        :U=distance to natural natural water table
C        :T=modifier for discharge/recharge to/from 
C            artificial water table 
C        :A=distance to natural natural water table
C        :B=modifier for discharge/recharge to/from 
C            artificial water table
C     RCHGD=flux conditions at bottom boundary from site file 
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
      RCHGFA=RCHGEA(M2,M1)
      RCHGFB=RCHGEB(M2,M1)
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
      RCHGFA=RCHGWA(M5,M4)
      RCHGFB=RCHGWB(M5,M4)
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
      RCHGFA=RCHGSA(M2,M1)
      RCHGFB=RCHGSB(M2,M1)
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
      RCHGFA=RCHGNA(M5,M4)
      RCHGFB=RCHGNB(M5,M4)
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
      RCHGFU=1.0 
      RCHGFT=RCHGD(M2,M1)
      RCHGFA=0.0
      RCHGFB=0.0
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
C     CDPTH,CDPTHI=current,initial surface elevation (m)
C     DLYR=layer thickness (m)
C     BKDS=bulk density (0=pond) (Mg m-3)
C     IRCHG,RCHQ*=topographic, user-set runoff boundary constraints
C     NN=boundary:N=1:NN=1 east,NN=2 west, 
C                :N=2:NN=1 south,NN=2 north
C     QR1,HQR1=runoff, convective heat from runoff (m3 t-1,MJ t-1)
C
      IF(L.EQ.NUM(N2,N1).AND.N.NE.3)THEN
C     WRITE(*,7744)'QRCHK1',I,J,NFZ,M,N1,N2,N4,N5,M4,M5,N,NN
C    2,IRCHG(NN,N,N2,N1),RCHQF,QRM(M,N2,N1) 
C    3,CDPTH(NU(N2,N1)-1,N2,N1),CDPTHI(N2,N1)
C    3,BKDS(NUI(N2,N1),N2,N1) 
C    3,DPTHW1,DTBLX(N2,N1),QRM(M,N2,N1)
      IF(IRCHG(NN,N,N2,N1).EQ.0.OR.RCHQF.EQ.0.0
     2.OR.(CDPTH(NU(N2,N1)-1,N2,N1).GT.CDPTHI(N2,N1)
     3.AND.BKDS(NUI(N2,N1),N2,N1).LE.ZERO))THEN
      QR1(N,NN,M5,M4)=0.0
      HQR1(N,NN,M5,M4)=0.0 
      ELSE
C
C     SURFACE BOUNDARY WATER FLUX
C
C     DPTHW1,DPTHW2=surface water depth of source,destination (m) 
C     ALT1,ALT2=elevation of source,destination (m)
C     XVOLT=surface water+ice in excess of litter
C        retention capacity (m3)
C     VOLWG=ground surface water retention capacity (m)
C     CDPTH(NU(N2,N1)-1,N2,N1)=soil surface elevation (m)
C     DTBLX=natural water table depth (m)
C     FSLOPE=partitions surface water flow in (N=1)EW,(N=2)NS
C       directions from ‘starts.f’
C     QR1,HQR1=runoff, convective heat from runoff in (1)EW,(2)NS
C        directions (m3 t-1,MJ t-1)
C     QR,HQR=aggregated runoff, convective heat from runoff 
C        used in ‘redist.f’ (m3 t-1,MJ t-1)
C     QRM=downslope runoff used in ‘trnsfr.f’ and ‘trnsfrs.f’ (m3 t-1)
C     QRV=downslope runoff velocity used in ‘erosion.f’ (m s-1)
C     QRMN=EW(N=1),NS(N=2) runoff used in ‘erosion.f’, ‘trnsfr.f’ 
C        and ‘trnsfrs.f’ (m3 t-1)
C     XN=direction parameter
C
C     RUNOFF
C
      DPTHW1=XVOLT(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)  
      DPTHW2=VOLWG(N2,N1)/AREA(3,NUM(N2,N1),N2,N1)  
      ALT1=ALTG(N2,N1)+DPTHW1 
      ALT2=ALTG(N2,N1)+DPTHW2
     2-XN*SLOPE(N,N2,N1)*DLYR(N,NUM(N2,N1),N2,N1)
C     WRITE(*,7744)'QRCHK2',I,J,NFZ,M,N1,N2,N4,N5,M4,M5,N,NN
C    2,IRCHG(NN,N,N2,N1),ALT1,ALT2,CDPTH(NU(N2,N1)-1,N2,N1)
C    3,DPTHW1,DTBLX(N2,N1),QRM(M,N2,N1)
      IF(ALT1.GT.ALT2.AND.XVOLT(N2,N1).GT.VOLWG(N2,N1)
     2.AND.CDPTH(NU(N2,N1)-1,N2,N1)-DPTHW1.LT.DTBLX(N2,N1))THEN
      QRQ1=(XVOLT(N2,N1)-VOLWG(N2,N1))
     2*VOLW(0,NY,NX)/XVOLT(N2,N1)*XNPXX
      QR1(N,NN,M5,M4)=-XN*(QRM(M,N2,N1)*FSLOPE(N,N2,N1)+QRQ1)*RCHQF
      HQR1(N,NN,M5,M4)=4.19*TK1(0,N2,N1)*QR1(N,NN,M5,M4)
      QR(N,NN,M5,M4)=QR(N,NN,M5,M4)+QR1(N,NN,M5,M4)
      HQR(N,NN,M5,M4)=HQR(N,NN,M5,M4)+HQR1(N,NN,M5,M4)
C     IF(QRM(M,N2,N1).GT.0.0)THEN
C     WRITE(*,7744)'QRBND',I,J,NFZ,M,N1,N2,N4,N5,M4,M5,N,NN
C    1,IRCHG(NN,N,N2,N1),RCHQF
C    2,QRM(M,N2,N1),QR1(N,NN,M5,M4),QRQ1,QR(N,NN,M5,M4) 
C    2,ALT1,ALT2,ALTG(N2,N1),XVOLT(N2,N1)
C    3,VOLWG(N2,N1),VOLWRX(N2,N1),ZM(N2,N1),ZS(N2,N1) 
C    4,VOLW1(0,N2,N1),VOLI1(0,N2,N1),DLYR(N,NUM(N2,N1),N2,N1)
C    5,XVOLT(N2,N1)-VOLWG(N2,N1),DTBLX(N2,N1)
C    6,DTBLX(N2,N1)-CDPTH(NU(N2,N1)-1,N2,N1)+DPTHW1
C    7,XVOLTM(M,N2,N1),XVOLWM(M,N2,N1),FSLOPE(N,N2,N1)
7744  FORMAT(A8,13I4,40E12.4)
C     ENDIF
C
C     RUNON
C
C     CDPTH(NU(N2,N1)-1,N2,N1)=soil surface elevation (m)
C     DPTHW1=surface water depth of source (m)
C     VX=water volume available for runon (m3)
C     XNPHX=time step from ‘wthr.f’(h t-1) 
C     DTBLX=natural water table depth (m)
C     FSLOPE=partitions surface water flow in (N=1)EW,(N=2)NS
C       directions from ‘starts.f’
C     QR1,HQR1=runon, convective heat from runon in (1)EW,(2)NS
C        directions (m3 t-1,MJ t-1)
C     QR,HQR=aggregated runon, convective heat from runon 
C        used in ‘redist.f’ (m3 t-1,MJ t-1)
C     QRM=downslope runoff used in ‘trnsfr.f’ and ‘trnsfrs.f’ (m3 t-1)
C     QRV=downslope runoff velocity used in ‘erosion.f’ (m s-1)
C     QRMN=EW(N=1),NS(N=2) runoff used in ‘erosion.f’, ‘trnsfr.f’ 
C        and ‘trnsfrs.f’ (m3 t-1)
C
      ELSEIF(CDPTH(NU(N2,N1)-1,N2,N1)-DPTHW1.GT.DTBLX(N2,N1))THEN
      VX=AMIN1(0.0,(DTBLX(N2,N1)-CDPTH(NU(N2,N1)-1,N2,N1)+DPTHW1)
     2*AREA(3,NUM(N2,N1),N2,N1))
      QRM(M,N2,N1)=VX*XNPHX 
      QRV(M,N2,N1)=0.0
      QR1(N,NN,M5,M4)=-XN*QRM(M,N2,N1)*FSLOPE(N,N2,N1)*RCHQF
      HQR1(N,NN,M5,M4)=4.19*TK1(0,N2,N1)*QR1(N,NN,M5,M4) 
      QR(N,NN,M5,M4)=QR(N,NN,M5,M4)+QR1(N,NN,M5,M4)
      HQR(N,NN,M5,M4)=HQR(N,NN,M5,M4)+HQR1(N,NN,M5,M4)
C     WRITE(*,7744)'QRBNB',I,J,NFZ,M,N1,N2,N4,N5,M4,M5,N,NN
C    1,IRCHG(NN,N,N2,N1),QRM(M,N2,N1)
C    2,QR1(N,NN,M5,M4),QR(N,NN,M5,M4),QRMN(M,N,NN,M5,M4)
C    2,ALTG(N2,N1),FSLOPE(N,N2,N1)
C    3,VOLWG(N2,N1),VOLWRX(N2,N1),ZM(N2,N1),ZS(N2,N1)
C    4,VOLW1(0,N2,N1),VOLI1(0,N2,N1),DLYR(N,NUM(N2,N1),N2,N1)
C    5,XVOLT(N2,N1)-VOLWG(N2,N1),DTBLX(N2,N1)
C    6,CDPTH(NU(N2,N1)-1,N2,N1),DPTHW1
C    7,XVOLTM(M,N2,N1),XVOLWM(M,N2,N1)
      ELSE
      QRM(M,N2,N1)=0.0 
      QRV(M,N2,N1)=0.0
      QR1(N,NN,M5,M4)=0.0
      HQR1(N,NN,M5,M4)=0.0
      ENDIF
      QRMN(M,N,NN,M5,M4)=QR1(N,NN,M5,M4)
      IFLBM(M,N,NN,M5,M4)=0
      ENDIF
C     IF(N5.EQ.6)THEN
C     WRITE(*,7745)'QRB',I,J,NFZ,M,N5,N4,M5,M4,N,NN
C    2,IFLBM(M,2,1,N5,N4) 
C    2,QRM(M,N2,N1),QRMN(M,N,NN,M5,M4),QR(N,NN,M5,M4) 
7745  FORMAT(A8,11I4,40E12.4)
C     ENDIF
C
C     BOUNDARY SNOW REDISTRIBUTION FROM SNOWPACK (DISABLED)
C
C     QS1,QW1,QI1=E,S,W,N snow,water,ice transfer (m3 t-1)
C     HQS1=convective heat transfer from snow,water,ice transfer 
C        (MJ t-1)
C     QS,QW,QI=aggregated snow,water,ice transfer used in ‘redist.f’
C        (m3 t-1)
C     HQS=aggregated convective heat transfer from
C        snow,water,ice transfer used in ‘redist.f’ (MJ t-1)
C     QSTN=EW(N=1),NS(N=2) snow transfer for solute flux calculation in
C        ‘trnsfr.f’ and ‘trnsfrs.f’ (m3 t-1)
C
C     IF(IRCHG(NN,N,N2,N1).EQ.0.OR.RCHQF.EQ.0.0
C    2.OR.ABS(QST(M,N2,N1)).LE.ZEROS(N2,N1)
C    3.OR.DPTHS(N2,N1).LE.ZERO)THEN
      QS1(N,NN,M5,M4)=0.0
      QW1(N,NN,M5,M4)=0.0
      QI1(N,NN,M5,M4)=0.0
      HQS1(N,NN,M5,M4)=0.0
C     ELSE
C     ALTS1=ALTG(N2,N1)+DPTHS(N2,N1) 
C     ALTS2=ALTG(N2,N1)+DPTHS(N2,N1)
C    2-XN*SLOPE(N,N2,N1)*DLYR(N,NUM(N2,N1),N2,N1)
C     IF(ALTS1.GT.ALTS2)THEN
C     QS1(N,NN,M5,M4)=-XN*QSM(N2,N1)*FSLOPE(N,N2,N1)*RCHQF 
C     QW1(N,NN,M5,M4)=-XN*QWM(N2,N1)*FSLOPE(N,N2,N1)*RCHQF 
C     QI1(N,NN,M5,M4)=-XN*QIM(N2,N1)*FSLOPE(N,N2,N1)*RCHQF 
C     HQS1(N,NN,M5,M4)=TK0(1,N2,N1)*(2.095*QS1(N,NN,M5,M4)
C    2+4.19*QW1(N,NN,M5,M4)+1.9274*QI1(N,NN,M5,M4))
C     QS(N,NN,M5,M4)=QS(N,NN,M5,M4)+QS1(N,NN,M5,M4)
C     QW(N,NN,M5,M4)=QW(N,NN,M5,M4)+QW1(N,NN,M5,M4)
C     QI(N,NN,M5,M4)=QI(N,NN,M5,M4)+QI1(N,NN,M5,M4)
C     HQS(N,NN,M5,M4)=HQS(N,NN,M5,M4)+HQS1(N,NN,M5,M4)
C     WRITE(*,7744)'QSBND',I,J,NFZ,M,N1,N2,N4,N5,M4,M5,N,NN
C    2,IRCHG(NN,N,N2,N1),QSM(N2,N1),QWM(N2,N1),QIM(N2,N1)
C    3,QS1(N,NN,M5,M4),QW1(N,NN,M5,M4),QI1(N,NN,M5,M4)
C    3,FSLOPE(N,N2,N1),RCHQF
C     ELSE
C     QS1(N,NN,M5,M4)=0.0
C     QW1(N,NN,M5,M4)=0.0
C     QI1(N,NN,M5,M4)=0.0
C     HQS1(N,NN,M5,M4)=0.0
C     ENDIF
      QSTN(M,N,NN,M5,M4)=0.0
      IFLBMS(M,N,NN,M5,M4)=0
C     ENDIF
      ELSE
      IF(N.NE.3)THEN
      QR1(N,NN,M5,M4)=0.0
      HQR1(N,NN,M5,M4)=0.0
      QS1(N,NN,M5,M4)=0.0
      QW1(N,NN,M5,M4)=0.0
      QI1(N,NN,M5,M4)=0.0
      HQS1(N,NN,M5,M4)=0.0
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
C     IDTBL=water table flag from site file
C     THETA1,THETAX=water content ahead,behind wetting front (m3 m-3)
C     VOLW1,VOLWH1=soil water micro-,macro-pore volume (m3)
C     VOLY=volume of soil micropore fraction (m3)
C     K1,KL=pore water class ahead,behind wetting front
C     CND1,CNDL=hydraulic conductivity ahead,behind wetting front
C        (m2 h MPa-1) 
C     POROS=soil porosity (m3 m-3)
C     FKSAT=reduction in soil surface Ksat from rainfall energy impact
C     FLWL,FLWLX=lower boundary micropore water flux (m3 t-1)
C     FLWHL=lower boundary macropore water flux (m3 t-1)
C     FLWLY=micropore discharge to artificial water table (m3 t-1)
C     FLWHLY=macropore discharge to artificial water table (m3 t-1)
C     HFLWL=convective heat from lower boundary water flux (MJ t-1)
C     XN=direction parameter
C     SLOPE for N=3:sin(vertical slope)=1 from ‘starts.f’
C     RCHG*=boundary flags from site file
C     XNPHX=time step water for fluxes from ‘wthr.f’(h t-1)
C
      IF(IDTBL(N2,N1).EQ.0.OR.N.EQ.3)THEN
      THETA1=AMAX1(THETZ(N3,N2,N1),AMIN1(POROS(N3,N2,N1)
     2,VOLW1(N3,N2,N1)/VOLY(N3,N2,N1)))
      THETAX=AMAX1(THETZ(N3,N2,N1),AMIN1(POROS(N3,N2,N1)
     2,VOLWX1(N3,N2,N1)/VOLY(N3,N2,N1)))
      K1=MAX(1,MIN(100,INT(100.0*(POROS(N3,N2,N1)
     2-THETA1)/POROS(N3,N2,N1))+1))
      IF(N3.EQ.NUM(NY,NX))THEN
      CND1=HCND(N,K1,N3,N2,N1)*FKSAT
      ELSE
      CND1=HCND(N,K1,N3,N2,N1)
      ENDIF
      VOLWZ=AMAX1(0.0,VOLW1(N3,N2,N1)*XNPXX) 
      FLWL(N,M6,M5,M4)=AMIN1(VOLWZ,AMAX1(-VOLWZ
     2,XN*0.0098*-ABS(SLOPE(N,N2,N1))*CND1
     3*AREA(3,N3,N2,N1)*XNPHX))/RCHGFU*RCHGFT 
      FLWLY(N,M6,M5,M4)=0.0 
      FLWLX(N,M6,M5,M4)=FLWL(N,M6,M5,M4) 
      VOLWHZ=AMAX1(0.0,VOLWH1(N3,N2,N1)*XNPXX) 
      FLWHL(N,M6,M5,M4)=AMIN1(VOLWHZ,AMAX1(-VOLWHZ 
     2,XN*0.0098*-ABS(SLOPE(N,N2,N1))*CNDH1(N,N2,N1) 
     3*AREA(3,N3,N2,N1)))/RCHGFU*RCHGFT*XNPHX 
      FLWHLY(N,M6,M5,M4)=0.0 
      HFLWL(N,M6,M5,M4)=4.19*TK1(N3,N2,N1)
     2*(FLWL(N,M6,M5,M4)+FLWHL(N,M6,M5,M4))
C     IF(L.EQ.2)THEN
C     WRITE(*,4443)'ABV',I,J,NFZ,M,N1,N2,N3,M4,M5,M6,N,NN,K1,XN
C    2,FLWL(N,M6,M5,M4),XNPXX,XNPHX,CND1,THETA1,THETZ(N3,N2,N1)
C    2,POROS(N3,N2,N1),VOLY(N3,N2,N1)
C    2,VOLPD,RCHGFU,RCHGFT,VOLX(N3,N2,N1),VOLW1(N3,N2,N1)
C    3,VOLWH1(N3,N2,N1),VOLPH1(N3,N2,N1),VOLPHD,VOLI1(N3,N2,N1)
C    4,VOLIH1(N3,N2,N1),VOLP1(N3,N2,N1),HFLWL(N,M6,M5,M4)
C    5,PSISM1(N3,N2,N1),PSISE(N3,N2,N1),FLWHL(N,M6,M5,M4),DTBLD(N2,N1)
C    6,SLOPE(N,N2,N1)
4443  FORMAT(A8,13I4,30E12.4)
C     ENDIF
      ELSE
C
C     MICROPORE DISCHARGE ABOVE WATER TABLE
C
C     IFLGU=micropore discharge flag to natural water table
C        :0=discharge enabled,1=discharge disabled 
C     PSISWD=water potential from water table slope (MPa)
C     XN,RCHG*=direction parameter,boundary flag
C     SLOPE=sin(lateral slope) from ‘starts.f’
C     DLYR=layer width (m)
C     DTBLG=water table slope from site file (m m-1) 
C     PSISWT=total water potential driving micropore discharge (MPa)
C     PSISA1,PSISO=matric,osmotic water potential (MPa)
C     DPTH,DTBLX=depth to layer midpoint,natural water table (MPa)
C     DPTHT=depth to internal water table (m)
C     FLWL=micropore discharge to natural water table (m3 t-1)
C     FLWLY=micropore discharge to natural water table (m3 t-1)
C     HFLWL=convective heat from discharge to natural water table 
C        (MJ t-1)
C     HCND=saturated hydraulic conductivity (m2 h MPa-1) 
C     AREA=lateral area of grid cell soil layer (m2)       
C     AREAU=fraction of AREA below natural water table
C
      IF(IFLGU.EQ.0.AND.RCHGFT.NE.0.0)THEN
      PSISWD=XN*0.0049*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISWT=AMIN1(0.0,-PSISA1(N3,N2,N1) 
     2+0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1))
     3-0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1)))
      IF(PSISWT.LT.0.0)PSISWT=PSISWT-PSISWD 
      THETAU=AMIN1(POROS(LU,N2,N1)
     2,VOLW1(LU,N2,N1)/VOLY(LU,N2,N1))
      KU=MAX(1,MIN(100,INT(100.0*(POROS(LU,N2,N1)
     2-THETAU)/POROS(LU,N2,N1))+1))
      FLWT=PSISWT*HCND(N,KU,LU,N2,N1)*AREA(N,N3,N2,N1)
     2*(1.0-AREAU(N3,N2,N1))/AMAX1(RCHGFU,1.0)*RCHGFT*XNPHX 
      FLWL(N,M6,M5,M4)=XN*FLWT 
      FLWLY(N,M6,M5,M4)=0.0 
      FLWLX(N,M6,M5,M4)=XN*FLWT 
      HFLWL(N,M6,M5,M4)=4.19*TK1(N3,N2,N1)*XN*FLWT
C     IF(NFZ.EQ.1)THEN
C     WRITE(*,4445)'DISCHMI',I,J,NFZ,M,N1,N2,N3,M4,M5,M6,N,NN,KU,LU
C    2,XN,FLWL(N,M6,M5,M4),FLWT,PSISWT,PSISWD,HCND(N,KU,LU,N2,N1)
C    3,AREA(N,N3,N2,N1),AREAU(N3,N2,N1),RCHGFU,RCHGFT
C    4,PSISE(N3,N2,N1),PSISA1(N3,N2,N1),DPTH(N3,N2,N1),DTBLX(N2,N1) 
C    2,PSISO(N3,N2,N1),0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1))
C    3,0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1))
4445  FORMAT(A8,14I4,30E14.6)
C     ENDIF
      ELSE
      FLWL(N,M6,M5,M4)=0.0
      FLWLY(N,M6,M5,M4)=0.0
      FLWLX(N,M6,M5,M4)=0.0
      HFLWL(N,M6,M5,M4)=0.0
      ENDIF
C
C     MACROPORE DISCHARGE ABOVE WATER TABLE
C
C     IFLGUH=macropore discharge flag to natural water table
C        :0=discharge enabled,1=discharge disabled 
C     PSISWD=water potential from water table slope (MPa)
C     XN,RCHG*=direction indicator,boundary flag from site file
C     SLOPE=sin(lateral slope) from ‘starts.f’
C     DLYR=layer width (m)
C     DTBLG=water table slope (m m-1) 
C     PSISWTH=water potential driving macropore discharge (MPa)
C     DPTHH,DTBLX=depth to layer macropore water,natural water table
C        (m)
C     DPTHT=depth to internal water table (m)
C     FLWTH,FLWTHL=macropore discharge unlimited,limited by macropore
C         water (m3 t-1)
C     CNDH1=macropore hydraulic conductivity (m2 h MPa-1)
C     FLWHL=macropore discharge to natural water table (m3 t-1)
C     FLWHLY=macropore discharge to artificial water table (m3 t-1)
C     HFLWL=convective heat from discharge to natural water table 
C        (MJ t-1)
C     HCND=saturated hydraulic conductivity (m2 h MPa-1)      
C     AREA=lateral area of grid cell soil layer (m2)       
C     AREAU=fraction of AREA below natural water table
C
      IF(IFLGUH.EQ.0.AND.RCHGFT.NE.0.0
     2.AND.VOLAH1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      PSISWD=XN*0.0049*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISWTH=0.0
     2+0.0098*(DPTHH-DTBLX(N2,N1))
     3-0.0098*AMAX1(0.0,DPTHH-DPTHT(N2,N1))
      IF(PSISWTH.LT.0.0)PSISWTH=PSISWTH-PSISWD
      FLWTH=PSISWTH*CNDH1(N3,N2,N1)*AREA(N,N3,N2,N1)  
     2*(1.0-AREAU(N3,N2,N1))/AMAX1(RCHGFU,1.0)*RCHGFT*XNPHX 
      FLWTHL=AMAX1(FLWTH,AMIN1(0.0,-(VOLWH1(N3,N2,N1)*XNPXX
     2+FLWHL(3,N3,N2,N1)-FLWHL(3,N3+1,N2,N1))))
      FLWHL(N,M6,M5,M4)=XN*FLWTHL 
      FLWHLY(N,M6,M5,M4)=0.0 
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)
     2+4.19*TK1(N3,N2,N1)*XN*FLWTHL
C     WRITE(*,4446)'DISCHMA',I,J,NFZ,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGUH
C    2,XN,FLWHL(N,M6,M5,M4),FLWTHL,FLWTH,PSISWTH,CNDH1(N3,N2,N1)
C    3,DPTH(N3,N2,N1),DLYR(3,N3,N2,N1),DPTHH,VOLWH1(N3,N2,N1)
C    4,VOLIH1(L,NY,NX),VOLAH1(N3,N2,N1),DTBLX(N2,N1),PSISWD
4446  FORMAT(A8,13I4,30E14.6)
      ELSE
      FLWHL(N,M6,M5,M4)=0.0
      FLWHLY(N,M6,M5,M4)=0.0 
      ENDIF
C
C     MICROPORE DISCHARGE ABOVE TILE DRAIN
C
C     IFLGD=micropore discharge flag to artificial water table
C        :0=discharge enabled,1=discharge disabled
C     PSISWD=water potential from water table slope (MPa)
C     XN,RCHG*=direction indicator,boundary flag from site 
C        or management file
C     SLOPE=sin(lateral slope) from ‘starts.f’
C     DLYR=layer width (m)
C     DTBLG=water table slope from site file (m m-1)
C     PSISWT=water potential driving micropore discharge (MPa)
C     PSISA1,PSISO=matric,osmotic water potential (MPa)
C     DPTH,DTBLY=depth to layer midpoint,artificial water table (m)
C     DPTHT=depth to internal water table (m)
C     FLWL=micropore discharge to artificial water table (m3 t-1)
C     FLWLY=micropore discharge to artificial water table (m3 t-1)
C     HFLWL=convective heat from discharge to artificial water table
C        (MJ t-1)
C     HCND=saturated hydraulic conductivity (m2 h MPa-1) 
C     AREA=lateral area of grid cell soil layer (m2)      
C     AREAUD=fraction of AREA below artificial water table
C
      IF(IFLGD.EQ.0.AND.RCHGFB.NE.0.0)THEN
      PSISWD=XN*0.0049*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISWT=AMIN1(0.0,-PSISA1(N3,N2,N1) 
     2+0.0098*(DPTH(N3,N2,N1)-DTBLY(N2,N1))
     3-0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1)))
      IF(PSISWT.LT.0.0)PSISWT=PSISWT-PSISWD 
      THETAD=AMIN1(POROS(LD,N2,N1)
     2,VOLW1(LD,N2,N1)/VOLY(LD,N2,N1))
      KD=MAX(1,MIN(100,INT(100.0*(POROS(LD,N2,N1)
     2-THETAD)/POROS(LD,N2,N1))+1))
      FLWT=PSISWT*HCND(N,KD,LD,N2,N1)*AREA(N,N3,N2,N1)
     2*(1.0-AREAUD(N3,N2,N1))/AMAX1(RCHGFA,1.0)*RCHGFB*XNPHX 
      FLWL(N,M6,M5,M4)=FLWL(N,M6,M5,M4)+XN*FLWT 
      FLWLY(N,M6,M5,M4)=FLWLY(N,M6,M5,M4)+XN*FLWT 
      FLWLX(N,M6,M5,M4)=FLWLX(N,M6,M5,M4)+XN*FLWT 
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+4.19*TK1(N3,N2,N1)*XN*FLWT
C     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,4445)'DISCHMD',I,J,NFZ,M,N1,N2,N3,M4,M5,M6,N,NN,KD,LD
C    2,XN,FLWL(N,M6,M5,M4),FLWLY(N,M6,M5,M4),FLWT,PSISWT,PSISWD
C    3,0.0098*(DPTH(N3,N2,N1)-DTBLY(N2,N1)),RCHGFU,RCHGFB
C    3,HCND(N,KD,LD,N2,N1),AREA(N,N3,N2,N1),AREAUD(N3,N2,N1)
C    4,PSISE(N3,N2,N1),PSISA1(N3,N2,N1),PSISO(N3,N2,N1)
C    5,DPTH(N3,N2,N1),DTBLY(N2,N1),DPTHT(N2,N1) 
C    3,0.0098*AMAX1(0.0,DPTH(N3,N2,N1)-DPTHT(N2,N1))
C    4,RCHGFA,RCHGFB
C     ENDIF
      ENDIF
C
C     MACROPORE DISCHARGE ABOVE TILE DRAIN
C
C     IFLGDH=macropore discharge flag to artificial water table
C        :0=discharge enabled,1=discharge disabled 
C     PSISWD=water potential from water table slope (MPa)
C     XN,RCHG*=direction indicator,boundary flag from site 
C        or management files
C     SLOPE=sin(lateral slope) from ‘starts.f’
C     DLYR=layer width (m)
C     DTBLG=water table slope from site file (m m-1)
C     PSISWTH=water potential driving macropore discharge (MPa)
C     DPTHH,DTBLY=depth to layer macropore water,artificial water table
C        (m) 
C     DPTHT=depth to internal water table (m)
C     FLWTH,FLWTHL=macropore discharge unlimited,limited by 
C        macropore water (m3 t-1)
C     CNDH1=macropore hydraulic conductivity (m2 h MPa-1) 
C     FLWHL=macropore discharge to artificial water table (m3 t-1)
C     FLWHLY=macropore discharge to artificial water table (m3 t-1)
C     HFLWL=convective heat from discharge to artificial water table
C        (MJ t-1)
C     HCND=saturated hydraulic conductivity (m2 h MPa-1)       
C     AREA=lateral area of grid cell soil layer (m2)       
C     AREAUD=fraction of AREA below artificial water table
C
      IF(IFLGDH.EQ.0.AND.RCHGFB.NE.0.0
     2.AND.VOLAH1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      PSISWD=XN*0.0049*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISWTH=0.0
     2+0.0098*(DPTHH-DTBLY(N2,N1))
     3-0.0098*AMAX1(0.0,DPTHH-DPTHT(N2,N1))
      IF(PSISWTH.LT.0.0)PSISWTH=PSISWTH-PSISWD
      FLWTH=PSISWTH*CNDH1(N3,N2,N1)*AREA(N,N3,N2,N1)  
     2*(1.0-AREAUD(N3,N2,N1))/AMAX1(RCHGFA,1.0)*RCHGFB*XNPHX 
      FLWTHL=AMAX1(FLWTH,AMIN1(0.0
     2,-(VOLWH1(N3,N2,N1)*XNPXX+FLWHL(3,N3,N2,N1)
     3-FLWHL(3,N3+1,N2,N1))))
      FLWHL(N,M6,M5,M4)=FLWHL(N,M6,M5,M4)+XN*FLWTHL 
      FLWHLY(N,M6,M5,M4)=FLWHLY(N,M6,M5,M4)+XN*FLWTHL 
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
C     DPTH,DTBLX=depth to layer midpoint,natural water table (m)
C     DPTHA=active layer depth (m)
C     VOLPD=air volume (m3)
C     PSISWD=water potential from water table slope (MPa)
C     XN,RCHG*=direction indicator,boundary flag from site file
C     SLOPE=sin(lateral slope) from ‘starts.f’
C     DLYR=layer width (m)
C     DTBLG=water table slope (m m-1) 
C     PSISUT=water potential driving micropore recharge (MPa)
C     PSISA1=matric water potential (MPa)
C     DPTH,DTBLX=depth to layer midpoint,natural water table (m) 
C     FLWU,FLWUL=micropore recharge unlimited,limited by micropore 
C        air volume (m3 t-1)
C     FLWL=micropore recharge from natural water table (m3 t-1)
C     HFLWL=convective heat from recharge from natural water table 
C        (MJ t-1)
C     HCND=saturated hydraulic conductivity (m2 h MPa-1)     
C     AREA=lateral area of grid cell soil layer (m2)      
C     AREAU=fraction of AREA below natural water table
C
C     IF(NFZ.EQ.1.AND.L.EQ.14)THEN
C     WRITE(*,4444)'RECHECK',I,J,NFZ,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGU 
C    2,DPTH(N3,N2,N1),DTBLX(N2,N1),DPTHA(N2,N1),VOLPD,VOLA1(N3,N2,N1)
C    4,BKDS(N3,N2,N1),VOLP1Z(N3,N2,N1),RCHGFT 
C     ENDIF
      IF(IFLGU.EQ.1
     2.AND.DPTH(N3,N2,N1).GE.DTBLX(N2,N1)
     2.AND.DPTHA(N2,N1).GT.DTBLX(N2,N1)
     3.AND.DPTH(N3,N2,N1).LT.DPTHA(N2,N1)
     4.AND.(VOLPD.GT.1.0E-03*VOLA1(N3,N2,N1)
     4.OR.BKDS(N3,N2,N1).LE.ZERO)
     4.AND.VOLP1Z(N3,N2,N1).GT.0.0
     5.AND.RCHGFT.NE.0.0)THEN
      PSISWD=XN*0.0049*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISUT=AMAX1(0.0,-PSISA1(N3,N2,N1) 
     2+0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1)))
      IF(PSISUT.GT.0.0)PSISUT=PSISUT+PSISWD 
      FLWU=PSISUT*HCND(N,1,N3,N2,N1)*AREA(N,N3,N2,N1) 
     2*AREAU(N3,N2,N1)/AMAX1(RCHGFU,1.0)*RCHGFT*XNPHX
      IF(BKDS(N3,N2,N1).GT.ZERO)THEN 
      FLWUL=AMIN1(FLWU,VOLPD*XNPXX)
      FLWUX=AMIN1(FLWU,VOLPXD*XNPXX)
      ELSE
      FLWUL=FLWU
      FLWUX=FLWU
      ENDIF
      FLWL(N,M6,M5,M4)=FLWL(N,M6,M5,M4)+XN*FLWUL 
      FLWLY(N,M6,M5,M4)=0.0 
      FLWLX(N,M6,M5,M4)=FLWLX(N,M6,M5,M4)+XN*FLWUX 
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+4.19*TK1(N3,N2,N1)
     2*XN*FLWUL
C     IF(N3.EQ.2)THEN
C     WRITE(*,4444)'RECHGMI',I,J,NFZ,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGU 
C    2,XN,FLWL(N,M6,M5,M4),FLWU,FLWUL
C    2,PSISUT,PSISA1(N3,N2,N1),THETW1,THETWL,PSISO(N3,N2,N1) 
C    3,0.0098*(DPTH(N3,N2,N1)-DTBLX(N2,N1)),RCHGFU,RCHGFT 
C    5,VOLPD,VOLA1(N3,N2,N1),VOLW1(N3,N2,N1),VOLW2(N3,N2,N1)
C    6,VOLI1(N3,N2,N1),VOLP1(N3,N2,N1) 
C    6,FLWL(3,N3,N2,N1),FLWL(3,N3+1,N2,N1)
C    4,HCND(N,1,N3,N2,N1),AREA(N,N3,N2,N1),AREAU(N3,N2,N1) 
C    6,DPTH(N3,N2,N1),DPTHA(N2,N1),DPTHT(N2,N1),DTBLX(N2,N1) 
C    6,DTBLD(N2,N1),VOLW1(N3,N2,N1),VOLI1(N3,N2,N1) 
C    7,VOLX(N3,N2,N1),VOLP1(N3,N2,N1),VOLP1Z(N3,N2,N1) 
C    8,FLWL(3,N3,N2,N1),FLWL(3,N6,N5,N4),QR1(N,NN,M5,M4)
C    9,DLYR(N,N3,N2,N1),DLYR(3,N3,N2,N1),PSISWD
C    1,CDPTH(N3,N2,N1),AREA(N,N3,N2,N1),SLOPE(N,N2,N1)
4444  FORMAT(A8,13I4,40E14.6)
C     ENDIF
      VOLPD=VOLPD-XN*FLWL(N,M6,M5,M4)
      VOLPXD=VOLPXD-XN*FLWLX(N,M6,M5,M4)
      ENDIF
C
C     MACROPORE RECHARGE BELOW WATER TABLE
C
C     PSISWD=water potential from water table slope (MPa)
C     XN,RCHG*=direction indicator,boundary flag from site file
C     SLOPE=sin(lateral slope) from ‘starts.f’
C     DLYR=layer width (m)
C     DTBLG=water table slope (m m-1)
C     PSISUTH=water potential driving macropore recharge (MPa)
C     DPTHH,DTBLX=depth to layer macropore water,natural water table
C        (m)
C     DPTHT=depth to internal water table (m)
C     CNDH1=macropore hydraulic conductivity
C     FLWUH,FLWUHL=macropore recharge unlimited,limited by 
C        macropore air volume
C     FLWHL=macropore discharge to natural water table
C     HFLWL=convective heat from discharge to natural water table
C     HCND=saturated hydraulic conductivity (m2 h MPa-1)      
C     AREA=lateral area of grid cell soil layer (m2)       
C     AREAU=fraction of AREA below natural water table
C
      IF(DPTHH.GE.DTBLX(N2,N1)
     2.AND.DPTHA(N2,N1).GT.DTBLX(N2,N1)
     2.AND.DPTH(N3,N2,N1).LT.DPTHA(N2,N1)
     2.AND.VOLPHD.GT.1.0E-03*VOLAH1(N3,N2,N1)
     2.AND.RCHGFT.NE.0.0)THEN
      PSISWD=XN*0.0049*SLOPE(N,N2,N1)*DLYR(N,N3,N2,N1)
     2*(1.0-DTBLG(N2,N1))
      PSISUTH=0.0
     2+0.0098*(DPTHH-DTBLX(N2,N1))
      IF(PSISUTH.GT.0.0)PSISUTH=PSISUTH+PSISWD
      FLWUH=PSISUTH*CNDH1(N3,N2,N1)*AREA(N,N3,N2,N1)  
     2*AREAU(N3,N2,N1)/AMAX1(RCHGFU,1.0)*RCHGFT*XNPHX 
      FLWUHL=AMIN1(FLWUH,VOLPHD*XNPXX)
      FLWHL(N,M6,M5,M4)=FLWHL(N,M6,M5,M4)+XN*FLWUHL
      FLWHLY(N,M6,M5,M4)=0.0 
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+4.19*TK1(N3,N2,N1)
     2*XN*FLWUHL
C     IF(N3.EQ.2)THEN
C     WRITE(*,4447)'RECHGMA',I,J,NFZ,M,N1,N2,N3,M4,M5,M6,N,NN,IFLGU 
C    2,XN,FLWHL(N,M6,M5,M4),FLWUH,FLWUHL,DPTHH,PSISUTH,VOLPHD
C    3,CNDH1(N3,N2,N1),DTBLX(N2,N1),CDPTH(N3,N2,N1),DPTHT(N2,N1) 
C    6,DTBLD(N2,N1),DPTH(N3,N2,N1),VOLAH1(L,NY,NX),VOLWH1(L,NY,NX)
C    2,VOLIH1(L,NY,NX) 
C    8,FLWHL(3,N3,N2,N1),FLWHL(3,N3+1,N2,N1),RCHGFU,AREA(N,N3,N2,N1)
C    9,DLYR(N,N3,N2,N1),DLYR(3,N3,N2,N1),PSISWD
C    1,SLOPE(N,N2,N1),AREAU(N3,N2,N1)
4447  FORMAT(A8,13I4,40E14.6)
C     ENDIF
      VOLPHD=VOLPHD-XN*FLWHL(N,M6,M5,M4)
      ENDIF
      ENDIF
C
C     SUBSURFACE HEAT SOURCE/SINK
C
C     HFLWL=heat flux across lower boundary (MJ t-1)
C     TK1=lower boundary soil temperature (K)
C     TKSD=deep source/sink temperature from geothermal flux
C         =mean annual temperature from site file (K)
C     TCNDG=thermal conductivity below lower boundary (m MJ h-1 K-1)
C     DPTHSK,CDPTH=depth of thermal sink/source, lower boundary (m)
C     XNPHX=time step for fluxes from ‘wthr.f’(h t-1)
C     
      IF(N.EQ.3.AND.IETYP(N2,N1).NE.-2)THEN
      HFLWL(N,M6,M5,M4)=HFLWL(N,M6,M5,M4)+(TK1(N3,N2,N1)
     2-TKSD(N2,N1))*TCNDG/(DPTHSK(N2,N1)-CDPTH(N3,N2,N1))
     3*AREA(N,N3,N2,N1)*XNPHX
      ENDIF
C
C     ADD BOUNDARY WATER, HEAT FLUXES TO INTERNAL WATER, HEAT FLUXES
C
      FLW(N,M6,M5,M4)=FLW(N,M6,M5,M4)+FLWL(N,M6,M5,M4)
      FLWX(N,M6,M5,M4)=FLWX(N,M6,M5,M4)+FLWLX(N,M6,M5,M4)
      FLWH(N,M6,M5,M4)=FLWH(N,M6,M5,M4)+FLWHL(N,M6,M5,M4)
      HFLW(N,M6,M5,M4)=HFLW(N,M6,M5,M4)+HFLWL(N,M6,M5,M4)
      FLWM(M,N,M6,M5,M4)=FLWL(N,M6,M5,M4)
      FLWHM(M,N,M6,M5,M4)=FLWHL(N,M6,M5,M4)
      FLWY(N,M6,M5,M4)=FLWY(N,M6,M5,M4)+FLWLY(N,M6,M5,M4)
      FLWHY(N,M6,M5,M4)=FLWHY(N,M6,M5,M4)+FLWHLY(N,M6,M5,M4)
      ENDIF
      ELSE
      FLWL(N,M6,M5,M4)=0.0
      FLWLY(N,M6,M5,M4)=0.0
      FLWLX(N,M6,M5,M4)=0.0
      FLWHL(N,M6,M5,M4)=0.0
      FLWHLY(N,M6,M5,M4)=0.0
      HFLWL(N,M6,M5,M4)=0.0
      FLWM(M,N,M6,M5,M4)=0.0
      FLWHM(M,N,M6,M5,M4)=0.0
      ENDIF
9575  CONTINUE
C
C     NET WATER AND HEAT FLUXES IN RUNOFF AND SNOW DRIFT 
C
C     TQR1,THQR1=net runoff,convective heat from runoff (m3 h-1,MJ h-1)
C     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, convective heat from 
C        snowpack drift (m3 h-1,MJ h-1)
C     QR1,HQR1=runoff, convective heat from runoff (m3 h-1,MJ h-1)
C     QS1,QW1,QI1,HQS1=snow,water,ice, convective heat from snowpack 
C        drift (m3 h-1,MJ h-1)
C
C     E and S fluxes
C
      IF(L.EQ.NUM(N2,N1).AND.N.NE.3)THEN
      DO 1202 NN=1,2
      TQR1(N2,N1)=TQR1(N2,N1)+QR1(N,NN,N2,N1) 
      THQR1(N2,N1)=THQR1(N2,N1)+HQR1(N,NN,N2,N1)
      TQS1(N2,N1)=TQS1(N2,N1)+QS1(N,NN,N2,N1) 
      TQW1(N2,N1)=TQW1(N2,N1)+QW1(N,NN,N2,N1) 
      TQI1(N2,N1)=TQI1(N2,N1)+QI1(N,NN,N2,N1)
      THQS1(N2,N1)=THQS1(N2,N1)+HQS1(N,NN,N2,N1) 
      IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      TQR1(N2,N1)=TQR1(N2,N1)-QR1(N,NN,N5,N4)
      THQR1(N2,N1)=THQR1(N2,N1)-HQR1(N,NN,N5,N4)
      ENDIF
      IF(IFLBMS(M,N,NN,N5,N4).EQ.0)THEN
      TQS1(N2,N1)=TQS1(N2,N1)-QS1(N,NN,N5,N4)
      TQW1(N2,N1)=TQW1(N2,N1)-QW1(N,NN,N5,N4)
      TQI1(N2,N1)=TQI1(N2,N1)-QI1(N,NN,N5,N4)
      THQS1(N2,N1)=THQS1(N2,N1)-HQS1(N,NN,N5,N4)
      ENDIF
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,6631)'TQR1',I,J,NFZ,M,N1,N2,N4,N5,N,NN
C    2,IFLBM(M,N,NN,N5,N4),TQR1(N2,N1),THQR1(N2,N1)
C    2,QR1(N,NN,N2,N1),QR1(N,NN,N5,N4)
C    3,QR(N,NN,N2,N1),QR(N,NN,N5,N4)
C    2,HQR1(N,NN,N2,N1),HQR1(N,NN,N5,N4)
C    3,HQR(N,NN,N2,N1),HQR(N,NN,N5,N4)
6631  FORMAT(A8,10I4,12E12.4) 
C     ENDIF
C
C     W and N fluxes
C
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      TQR1(N2,N1)=TQR1(N2,N1)-QR1(N,NN,N5B,N4B)
      THQR1(N2,N1)=THQR1(N2,N1)-HQR1(N,NN,N5B,N4B)
      TQS1(N2,N1)=TQS1(N2,N1)-QS1(N,NN,N5B,N4B)
      TQW1(N2,N1)=TQW1(N2,N1)-QW1(N,NN,N5B,N4B)
      TQI1(N2,N1)=TQI1(N2,N1)-QI1(N,NN,N5B,N4B)
      THQS1(N2,N1)=THQS1(N2,N1)-HQS1(N,NN,N5B,N4B)
C     IF(I.GT.350.AND.NX.EQ.1)THEN
C     WRITE(*,6631)'TQRB1',I,J,NFZ,M,N1,N2,N4B,N5B,N,NN
C    2,IFLBM(M,N,NN,N5B,N4B),TQR1(N2,N1),THQR1(N2,N1)
C    2,QR1(N,NN,N5B,N4B),HQR1(N,NN,N5B,N4B)  
C    2,QR(N,NN,N5B,N4B),HQR(N,NN,N5B,N4B)  
C     ENDIF
      ENDIF
C
C     IFLBH=runoff direction flag used in ‘erosion.f’, ‘trnsfr.f’ 
C        and ‘trnsfrs.f’ (0 = E or S, 1 = W or N)
C     IFLBHS=snow drift direction flag used in ‘trnsfr.f’ 
C        and ‘trnsfrs.f’ (0 = E or S, 1 = W or N)
C
      IF(M.EQ.NPH)THEN
      IFLBH(N,NN,N5,N4)=IFLBM(M,N,NN,N5,N4)
      IFLBHS(N,NN,N5,N4)=IFLBMS(M,N,NN,N5,N4)
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      IFLBH(N,NN,N5B,N4B)=IFLBM(M,N,NN,N5B,N4B)
      IFLBHS(N,NN,N5B,N4B)=IFLBMS(M,N,NN,N5B,N4B)
      ENDIF
      ENDIF
1202  CONTINUE
      ENDIF
C
C     NET WATER,VAPOR AND HEAT FLUXES THROUGH SOIL AND SNOWPACK     
C
C     TFLWL,THFLWL=net soil water micropore,macropore flux (m3 t-1)
C     TFLVL=net soil vapor flux (m3 t-1)
C     THFLWL=net soil convective+conductive heat flux (MJ t-1)
C     FLWL,FLWHL=soil water micropore,macropore flux (m3 t-1) 
C     HFLWL=soil convective+conductive heat flux (MJ t-1) 
C
C     FIND NEXT EXISTING DESTINATION SOIL LAYER
C
C     VOLX=soil volume (m3)
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
      TFLVL(N3,N2,N1)=TFLVL(N3,N2,N1)+FLVL(N,N3,N2,N1)
     2-FLVL(N,N6,N5,N4)
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
      TFLVL(N3,N2,N1)=0.0
      TFLWLX(N3,N2,N1)=0.0
      TFLWHL(N3,N2,N1)=0.0
      THFLWL(N3,N2,N1)=0.0
      ENDIF
      ENDIF
9580  CONTINUE
      VOLWH2(N3,N2,N1)=AMAX1(0.0,VOLWH1(N3,N2,N1)+TFLWHL(N3,N2,N1)) 
      VOLIH2(N3,N2,N1)=AMAX1(0.0,VOLIH1(N3,N2,N1))
C
C     EVAPORATION-CONDENSATION IN SOIL BELOW SURFACE LAYER
C
C     VOLW2,VOLI2=soil micropore water,ice volume (m3)
C     VOLWH1,VOLIH1=macropore water,ice volume (m3)
C     VOLV2,VOLPM=soil vapor,air volume (m3)
C     VPLV=soil saturated vapor concentration (m3 m-3) 
C     TK1,PSISV1=soil temperature,matric+osmotic potential (K,MPa)
C     FLVLX,FLVLW=soil evaporation-condensation unlimited,
C        limited by vapor (m3 t-1)
C     HFLVLW=soil latent heat flux (MJ t-1)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C     WFLVL,HFLVL=soil evaporation-condensation,
C        latent heat flux (m3 t-1,MJ t-1)
C
      PSISV1=PSISM1(N3,N2,N1)+PSISO(N3,N2,N1)
      IF(N3.GT.NUM(NY,NX))THEN
      IF(VOLPM(M,N3,N2,N1).GT.ZEROS(N2,N1))THEN
      VPLV=2.173E-03/TK1(N3,N2,N1)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TK1(N3,N2,N1)))
     3*EXP(18.0*PSISV1/(8.3143*TK1(N3,N2,N1)))
      FLVLX=VOLV2(N3,N2,N1)-VPLV*VOLPM(M,N3,N2,N1)
      FLVLW=AMAX1(FLVLX,-AMAX1(0.0,VOLW2(N3,N2,N1))*XNPXX)
      HFLVLW=VAP*FLVLW
      WFLVL(N3,N2,N1)=FLVLW 
      HFLVL(N3,N2,N1)=HFLVLW 
      VOLW2(N3,N2,N1)=VOLW2(N3,N2,N1)+FLVLW
      VOLV2(N3,N2,N1)=VOLV2(N3,N2,N1)-FLVLW
      ELSE
      WFLVL(N3,N2,N1)=0.0
      HFLVL(N3,N2,N1)=0.0
      ENDIF
C     IF(N3.EQ.2)THEN     
C     WRITE(*,7762)'WFLVL',I,J,NFZ,M,N1,N2,N3,WFLVL(N3,N2,N1) 
C    2,FLVLX,FLVLW,HFLVLW,VPLV,TK1(N3,N2,N1),PSISV1 
C     2,PSISM1(N3,N2,N1),PSISO(N3,N2,N1),VPLV*VOLP2(N3,N2,N1)
C    3,VOLV2(N3,N2,N1),VOLV1(N3,N2,N1),VOLPM(M,N3,N2,N1)
C     4,VOLP1(N3,N2,N1),HFLVL(N3,N2,N1)
7762  FORMAT(A8,7I4,20E12.4)
C     ENDIF
C
C     FREEZE-THAW IN SOIL LAYER MICROPORE FROM NET CHANGE IN SOIL 
C     LAYER HEAT STORAGE
C
C     TFREEZ=micropore freezing temperature (K)
C     333.0=latent heat of freezing (MJ m-3)
C     TK1,PSISV1=soil temperature,matric+osmotic potential (K,MPa)
C     VOLW2,VOLI2=micropore water,ice volume (m3)
C     VOLT=soil volume (m3)
C     VHCP1=soil micropore heat capacity (MJ K-1)
C     HFLFM1,HFLFM=latent heat from micropore freeze-thaw
C        unlimited,limited by water,ice (MJ t-1)
C     WFLFL=soil water flux from micropore freeze-thaw (m3 t-1)
C
      TFREEZ=-9.0959E+04/(PSISV1-333.0)
      IF((TK1(N3,N2,N1).LT.TFREEZ
     2.AND.VOLW2(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1))
     4.OR.(TK1(N3,N2,N1).GT.TFREEZ
     5.AND.VOLI2(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1)))THEN
      HFLFM1=VHCP1(N3,N2,N1)*(TFREEZ-TK1(N3,N2,N1))
     2/(1.0+6.2913E-03*TFREEZ)
      IF(HFLFM1.LT.0.0)THEN
      HFLFM=AMAX1(-333.0*DENSI*VOLI2(N3,N2,N1)*XNPXX,HFLFM1)
      ELSE
      HFLFM=AMIN1(333.0*VOLW2(N3,N2,N1)*XNPXX,HFLFM1)
      ENDIF
      WFLFL(N3,N2,N1)=-HFLFM/333.0
      ELSE
      HFLFM=0.0
      WFLFL(N3,N2,N1)=0.0
      ENDIF
      ELSE
      HFLFM=0.0
      ENDIF 
C
C     FREEZE-THAW IN SOIL LAYER MACROPORE FROM NET CHANGE IN SOIL 
C     LAYER HEAT STORAGE
C
C     VOLWH12,VOLIH2=macropore water,ice volume (m3)
C     TK1=soil temperature (K)
C     VHCP1HX=soil macropore heat capacity (MJ K-1)
C     HFLFH1,HFLFH=latent heat from macropore freeze-thaw
C        unlimited,limited by water,ice (MJ t-1)
C     WFLFLH=soil water flux from macropore freeze-thaw (m3 t-1)
C
      IF((TK1(N3,N2,N1).LT.273.15
     2.AND.VOLWH2(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1))
     3.OR.(TK1(N3,N2,N1).GT.273.15
     4.AND.VOLIH2(N3,N2,N1).GT.ZERO*VOLT(N3,N2,N1)))THEN
      VHCP1HX=4.19*VOLWH2(N3,N2,N1)+1.9274*VOLIH2(N3,N2,N1)
      HFLFH1=VHCP1HX*(TFREEZ-TK1(N3,N2,N1))
     2/(1.0+6.2913E-03*TFREEZ)
      IF(HFLFH1.LT.0.0)THEN
      HFLFH=AMAX1(-333.0*DENSI*VOLIH2(N3,N2,N1)*XNPXX,HFLFH1)
      ELSE
      HFLFH=AMIN1(333.0*VOLWH2(N3,N2,N1)*XNPXX,HFLFH1)
      ENDIF
      WFLFLH(N3,N2,N1)=-HFLFH/333.0
      ELSE
      HFLFH=0.0
      WFLFLH(N3,N2,N1)=0.0
      ENDIF
      IF(N3.EQ.NUM(NY,NX))THEN
      HFLFL(N3,N2,N1)=HFLFL(N3,N2,N1)+HFLFH
      ELSE
      HFLFL(N3,N2,N1)=HFLFM+HFLFH
      ENDIF
C     IF(N3.EQ.9)THEN
C     WRITE(*,4359)'WFLFL',I,J,NFZ,M,N1,N2,N3
C    2,HFLFL(N3,N2,N1),WFLFL(N3,N2,N1),WFLFLH(N3,N2,N1),HFLFM,HFLFH
C    3,HFLFM1,TFREEZ,PSISV1,TK1(N3,N2,N1),DENSI 
C    4,VOLW2(N3,N2,N1),VOLI2(N3,N2,N1),VHCP1(N3,N2,N1)
C     ENDIF
C
C     TOTAL AND ACCUMULATED FREEZE-THAW FLUXES
C
C     XWFLVL=total evaporation-condensation 
C        used in ‘redist.f’ (m3 t-1)
C     XHFLVL=total evaporation-condensation latent heat flux 
C        used in ‘redist.f’ (MJ t-1)
C     XWFLFL,XWFLFH=total freeze-thaw flux in micropores,macropores
C        used in ‘redist.f’(m3 t-1)
C     XHFLFL=total freeze-thaw latent heat flux
C        used in ‘redist.f’ (MJ t-1)
C     TWFLVL=net evaporation-condensation (m3 t-1) 
C     THFLVL=net evaporation-condensation latent heat flux (MJ t-1) 
C     TWFLFL,TWFLFH=net freeze-thaw in micropores,macropores (m3 t-1)
C     THFLFL=net freeze-thaw latent heat flux (MJ t-1) 
C
      XWFLVL(N3,N2,N1)=XWFLVL(N3,N2,N1)+WFLVL(N3,N2,N1)
      XHFLVL(N3,N2,N1)=XHFLVL(N3,N2,N1)+HFLVL(N3,N2,N1)
      TWFLVL(N3,N2,N1)=TWFLVL(N3,N2,N1)+WFLVL(N3,N2,N1)
      THFLVL(N3,N2,N1)=THFLVL(N3,N2,N1)+HFLVL(N3,N2,N1)
      XWFLFL(N3,N2,N1)=XWFLFL(N3,N2,N1)+WFLFL(N3,N2,N1)
      XWFLFH(N3,N2,N1)=XWFLFH(N3,N2,N1)+WFLFLH(N3,N2,N1)
      XHFLFL(N3,N2,N1)=XHFLFL(N3,N2,N1)+HFLFL(N3,N2,N1)
      TWFLFL(N3,N2,N1)=TWFLFL(N3,N2,N1)+WFLFL(N3,N2,N1)
      TWFLFH(N3,N2,N1)=TWFLFH(N3,N2,N1)+WFLFLH(N3,N2,N1)
      THFLFL(N3,N2,N1)=THFLFL(N3,N2,N1)+HFLFL(N3,N2,N1)
C     IF(N3.EQ.1)THEN
C     WRITE(*,4359)'HFLF',I,J,NFZ,M,N1,N2,N3,TWFLFL(N3,N2,N1)
C    2,WFLFL(N3,N2,N1),XWFLFL(N3,N2,N1),VHCP1X,VHCP1AX 
C    3,TWFLVL(N3,N2,N1),WFLVL(N3,N2,N1)
C    2,HFLFM1,HFLFM,TFREEZ,TK1(N3,N2,N1),TK1X,VHCP1(N3,N2,N1)
C    3,HFLFL(N3,N2,N1),THFLFL(N3,N2,N1),PSISML 
C    2,TFLWL(N3,N2,N1),FINHL(N3,N2,N1),FLU1(N3,N2,N1),THFLWL(N3,N2,N1)
C    4,VOLW2(N3,N2,N1),VOLI2(N3,N2,N1),VOLW1(N3,N2,N1),VOLI1(N3,N2,N1)
C    5,PSISA1(N3,N2,N1),PSISM(N3,N2,N1),PSISO(N3,N2,N1),PSISML 
C    5,VOLWH1(N3,N2,N1),VOLIH1(N3,N2,N1)
C    5,FGRD(N3,N2,N1),FMAC(N3,N2,N1)
4359  FORMAT(A8,7I4,40E14.6)
C     ENDIF
C
C     INFILTRATION OF WATER FROM MACROPORES INTO MICROPORES
C
C     VOLWH1=macropore volume (m3)
C     FINHX,FINHL=macro-micropore transfer unlimited,limited by
C        water,air volume (m3 t-1)
C     FINHM=macro-micropore transfer for use in ‘trnsfr.f’ (m3 t-1) 
C     HCND=saturated hydraulic conductivity (m2 h MPa-1)
C     PSISE,PSISA1=saturation,matric water potentials (MPa)
C     PHOL,HRAD=path length between,radius of macropores from
C        ‘hour1.f’(m)
C     VOLW2,VOLP2=current micropore water,air volume (m3)
C     VOLWH2,VOLPH2=current macropore water,air volume (m3)
C     XNPHX=time step for water fluxes from ‘wthr.f’(h t-1)
C
      IF(VOLWH1(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      FINHX=6.283*HCND(2,1,N3,N2,N1)*AREA(3,N3,N2,N1)
     2*(PSISE(N3,N2,N1)-PSISA1(N3,N2,N1))
     3/LOG(PHOL(N3,N2,N1)/HRAD(N3,N2,N1))*XNPHX 
      VOLW2(N3,N2,N1)=VOLW2(N3,N2,N1)+TFLWL(N3,N2,N1)+TWFLVL(N3,N2,N1) 
     2+TWFLFL(N3,N2,N1)+FLU1(N3,N2,N1)
      VOLP2(N3,N2,N1)=AMAX1(0.0,VOLA1(N3,N2,N1)-VOLW2(N3,N2,N1)
     2-VOLI2(N3,N2,N1))
      VOLWH2(N3,N2,N1)=VOLWH2(N3,N2,N1)+TFLWHL(N3,N2,N1)
     2+TWFLFH(N3,N2,N1) 
      VOLIH2(N3,N2,N1)=VOLIH2(N3,N2,N1)-TWFLFH(N3,N2,N1)/DENSI 
      VOLPH2(N3,N2,N1)=AMAX1(0.0,VOLAH1(N3,N2,N1)-VOLWH2(N3,N2,N1)
     2-VOLIH2(N3,N2,N1))
      IF(FINHX.GT.0.0)THEN
      FINHL(N3,N2,N1)=AMAX1(0.0,AMIN1(FINHX,VOLWH2(N3,N2,N1)
     2,VOLP2(N3,N2,N1)))
      ELSE
      FINHL(N3,N2,N1)=AMIN1(0.0,AMAX1(FINHX,-VOLPH2(N3,N2,N1)
     2,-VOLW2(N3,N2,N1)))
      ENDIF
      FINHM(M,N3,N2,N1)=FINHL(N3,N2,N1)
      FINH(N3,N2,N1)=FINH(N3,N2,N1)+FINHL(N3,N2,N1)
C     IF(N3.LE.5)THEN
C     WRITE(*,3366)'FINHL',I,J,NFZ,M,N1,N2,N3,IFLGH,FINHL(N3,N2,N1)
C    3,FINHX,VOLWH2(N3,N2,N1),VOLP2(N3,N2,N1),VOLPH2(N3,N2,N1)
C    2,VOLW2(N3,N2,N1)
C    4,VOLWH1(N3,N2,N1),VOLPH1(N3,N2,N1),VOLP1(N3,N2,N1)
C    4,PSISA1(N3,N2,N1),HCND(2,1,N3,N2,N1),PHOL(N3,N2,N1)
C    5,HRAD(N3,N2,N1),TFLWHL(N3,N2,N1),TWFLFH(N3,N2,N1)
C    6,VOLA1(N3,N2,N1),VOLW2(N3,N2,N1),VOLI1(N3,N2,N1)
C    7,VOLW1(N3,N2,N1),TFLWL(N3,N2,N1),TWFLVL(N3,N2,N1) 
C    2,TWFLFL(N3,N2,N1),FLU1(N3,N2,N1),VOLA(N3,N2,N1) 
3366  FORMAT(A8,8I4,30E12.4)
C     ENDIF
      ELSE
      FINHL(N3,N2,N1)=0.0
      FINHM(M,N3,N2,N1)=0.0
      ENDIF
9585  CONTINUE
9590  CONTINUE
9595  CONTINUE
C
C     UPDATE STATE VARIABLES FROM NET FLUXES CALCULATED ABOVE
C
      IF(M.LT.NPH)THEN
      DO 9795 NX=NHW,NHE
      DO 9790 NY=NVN,NVS
C
C     SNOWPACK WATER, ICE, SNOW AND TEMPERATURE
C
C     VOLS0,VOLW0,VOLW0,VOLI0=snow,water,vapor,ice volumes in snowpack
C        (m3)
C     TFLWS,TFLWW,TFLWV,TFLWI=net snow,water,vapor,ice flux (m3 t-1)
C     WFLFS,WFLFI=snow-water,ice-water freeze-thaw flux (m3 t-1)
C     DENSI=ice density (Mg m-3)
C     VHCPWM=snowpack volumetric heat capacity (MJ K-1)
C     TK0=snowpack temperature (K)
C     THFLWW=total snowpack conductive+convective heat flux (MJ K-1)
C     HFLF0,HFLV0=snowpack latent heat flux from freeze-thaw,
C        evaporation-condensation (MJ t-1)
C
      DPTHS0(NY,NX)=0.0
      DO 9780 L=1,JS
      VOLS0(L,NY,NX)=VOLS0(L,NY,NX)+TFLWS(L,NY,NX)-WFLFS(L,NY,NX)
     2+WFLVS(L,NY,NX) 
      VOLW0(L,NY,NX)=VOLW0(L,NY,NX)+TFLWW(L,NY,NX)+WFLFS(L,NY,NX)
     2+WFLFI(L,NY,NX)+WFLVW(L,NY,NX) 
      VOLV0(L,NY,NX)=VOLV0(L,NY,NX)+TFLWV(L,NY,NX)-WFLVS(L,NY,NX)
     2-WFLVW(L,NY,NX) 
      VOLI0(L,NY,NX)=VOLI0(L,NY,NX)+TFLWI(L,NY,NX)
     2-WFLFI(L,NY,NX)/DENSI
      DPTHS0(NY,NX)=DPTHS0(NY,NX)+DLYRS0(L,NY,NX)
      ENGY0=VHCPWM(M,L,NY,NX)*TK0(L,NY,NX)
      VHCPWM(M+1,L,NY,NX)=2.095*VOLS0(L,NY,NX)
     2+4.19*(VOLW0(L,NY,NX)+VOLV0(L,NY,NX))
     2+1.9274*VOLI0(L,NY,NX)
      IF(VHCPWM(M+1,L,NY,NX).GT.VHCPWX(NY,NX))THEN
      TK0(L,NY,NX)=(ENGY0+THFLWW(L,NY,NX)+HFLF0(L,NY,NX) 
     2+HFLV0(L,NY,NX))/VHCPWM(M+1,L,NY,NX)
      ELSEIF(L.EQ.1)THEN
      TK0(L,NY,NX)=TKQG(M,NY,NX)
      ELSE
      TK0(L,NY,NX)=TK0(L-1,NY,NX)
      ENDIF
C     IF(L.EQ.1)THEN
C     WRITE(*,7753)'TK0',I,J,NFZ,M,NX,NY,L
C    2,TK0(L,NY,NX),ENGY0,THFLWW(L,NY,NX),HFLF0(L,NY,NX) 
C    2,HFLV0(L,NY,NX),VHCPWM(M+1,L,NY,NX),TKW(L,NY,NX)
C    3,VOLS1(L,NY,NX)
C    3,VOLS0(L,NY,NX),VOLW0(L,NY,NX),VOLV0(L,NY,NX),VOLI0(L,NY,NX)
C    4,TFLWV(L,NY,NX),WFLVS(L,NY,NX),WFLVW(L,NY,NX) 
C    2,THFLWW(L,NY,NX),HFLF0(L,NY,NX),HFLV0(L,NY,NX) 
C    3,VHCPWM(M+1,L,NY,NX),HFLW0W(L,NY,NX),HFLWRLW,HFLWLW
C    4,HFLW0W(L,NY,NX) 
C    2,WFLFS(L,NY,NX),WFLFI(L,NY,NX),WFLVS(L,NY,NX),WFLVW(L,NY,NX)
C    3,TFLWS(L,NY,NX),TFLWW(L,NY,NX),TFLWV(L,NY,NX),TFLWI(L,NY,NX)
C    4,TQS1(NY,NX),TQW1(NY,NX)
C    3,XFLWS(L,NY,NX),XFLWW(L,NY,NX),XFLWI(L,NY,NX)
C    4,XWFLVS(L,NY,NX),XWFLVW(L,NY,NX) 
C    4,XHFLWW(L,NY,NX),HFLSW(L,NY,NX)
C    5,VHCPWM(M+1,L,NY,NX) 
7753  FORMAT(A8,7I4,80E14.6)
C     ENDIF
9780  CONTINUE
C
C     SNOW DRIFT
C
C     VOLS0,VOLI0,VOLW0,VOLV0=snowpack snow,ice,water,vapor volume (m3)
C     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from 
C        snowpack runoff (m3 t-1)
C     VHCPWM=snowpack volumetric heat capacity (MJ K-1)
C     TK0=snowpack temperature (K)
C     TKQG=air temperature at ground surface (K)
C
      VOLS0(1,NY,NX)=VOLS0(1,NY,NX)+TQS1(NY,NX)
      VOLW0(1,NY,NX)=VOLW0(1,NY,NX)+TQW1(NY,NX)
      VOLI0(1,NY,NX)=VOLI0(1,NY,NX)+TQI1(NY,NX)
      ENGY0=VHCPWM(M+1,1,NY,NX)*TK0(1,NY,NX)
      VHCPWM(M+1,1,NY,NX)=2.095*VOLS0(1,NY,NX)
     2+4.19*(VOLW0(1,NY,NX)+VOLV0(1,NY,NX))
     2+1.9274*VOLI0(1,NY,NX)
      IF(VHCPWM(M+1,1,NY,NX).GT.VHCPWX(NY,NX))THEN
      TK0(1,NY,NX)=(ENGY0+THQS1(NY,NX))
     2/VHCPWM(M+1,1,NY,NX)
      ELSE
      TK0(1,NY,NX)=TKQG(M,NY,NX)
      ENDIF
C
C     IF SNOWPACK DISAPPEARS ALL MATERIAL,HEAT IS TRANSFERRED 
C     TO LITTER
C
C     VHCPW,VHCPWX=current, minimum snowpack layer heat capacities 
C        (MJ K-1)
C     TKQG=air temperature at ground surface (K)
C     XFLWSX,XFLWWX,XFLWVX,XFLWIX=snow,water,vapor,ice flux 
C        to litter (m3 t-1)
C     XHFLWX=heat flux to litter (MJ t-1)
C     VOLS0,VOLI0,VOLW0,VOLV0=snowpack snow,ice,water,vapor volume (m3) 
C     VOLI1,VOLW1,VOLV1=litter ice,water,vapor volume (m3) 
C     TK1=litter temperature (K)
C     VHCP1=litter heat capacity (MJ K-1)
C     ORGC,ORGCC=litter organic C, charcoal (g)
C
      IF(VHCPW(1,NY,NX).LE.VHCPWX(NY,NX)
     2.AND.TKQG(M,NY,NX).GT.273.15)THEN
      FLWS=VOLS0(1,NY,NX) 
      FLWW=VOLW0(1,NY,NX)
      FLWV=VOLV0(1,NY,NX)
      FLWI=VOLI0(1,NY,NX) 
      HFLWS=(4.19*(FLWW+FLWV)+2.095*FLWS+1.9274*FLWI)*TK0(1,NY,NX)
      XFLWSX(NY,NX)=XFLWSX(NY,NX)+FLWS
      XFLWWX(NY,NX)=XFLWWX(NY,NX)+FLWW
      XFLWVX(NY,NX)=XFLWVX(NY,NX)+FLWV
      XFLWIX(NY,NX)=XFLWIX(NY,NX)+FLWI
      XHFLWX(NY,NX)=XHFLWX(NY,NX)+HFLWS
      VOLS0(1,NY,NX)=VOLS0(1,NY,NX)-FLWS
      VOLW0(1,NY,NX)=VOLW0(1,NY,NX)-FLWW
      VOLV0(1,NY,NX)=VOLV0(1,NY,NX)-FLWV
      VOLI0(1,NY,NX)=VOLI0(1,NY,NX)-FLWI
      VOLW1(0,NY,NX)=VOLW1(0,NY,NX)+FLWW
      VOLV1(0,NY,NX)=VOLV1(0,NY,NX)+FLWV
      VOLI1(0,NY,NX)=VOLI1(0,NY,NX)+FLWI+FLWS/DENSI
      ENGY1=VHCP1(0,NY,NX)*TK1(0,NY,NX)
      VHCP1(0,NY,NX)=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))
     2+4.19*(VOLW1(0,NY,NX)+VOLV1(0,NY,NX))
     3+1.9274*VOLI1(0,NY,NX)
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      TK1(0,NY,NX)=(ENGY1+HFLWS)/VHCP1(0,NY,NX)
      ELSE
      TK1(0,NY,NX)=TKQG(M,NY,NX)
      ENDIF
      ENDIF
C
C     SURFACE LITTER WATER AND TEMPERATURE
C
C     VOLW1,VOLV1,VOLI1,VOLP1,VOLA1=litter ice,water,vapor,air,porous
C        volume (m3)
C     FLWRL=total water flux into litter (m3 t-1)
C     WFLVR,HFLVR=litter evaporation-condensation,
C        latent heat flux (m3 t-1,MJ t-1)
C     WFLFR,HFLFR=surface litter freeze-thaw,latent 
C        heat flux (m3 t-1,MJ t-1)
C     TQR1,THQR1=net runoff,convective heat from runoff (m3 t-1,MJ t-1)
C     DENSI=ice density (Mg m-3)
C     VOLWM,VOLPM=surface water,air content for use in ‘trnsfr.f’ (m3)
C     VHCP1=volumetric heat capacity of litter (MJ K-1)
C     VOLR=dry litter volume (m3)
C     THETWX,THETIX,THETPX=water,ice,air concentrations (m3 m-3)
C     VHCP1=volumetric heat capacity of litter (MJ K-1)
C     TK1=litter temperature (K)
C     HFLWRL=litter total conductive+convective heat flux (MJ K-1)
C     HFLXF=heat added by litter combustion (MJ t-1)
C     TKSM=soil temperature used in ‘trnsfr.f’(K)
C     TKGS=ground surface temperature (K)
C     FSNW,FSNX=snow,snow-free cover fractions 
C     BAREW,CVRDW=soil,litter cover fractions accounting for excess
C        surface water
C
      VOLW1(0,NY,NX)=VOLW1(0,NY,NX)+FLWRL(NY,NX)+WFLVR(NY,NX)
     2+WFLFR(NY,NX)+TQR1(NY,NX)
      VOLV1(0,NY,NX)=VOLV1(0,NY,NX)+FLVRL(NY,NX)-WFLVR(NY,NX)
      VOLI1(0,NY,NX)=VOLI1(0,NY,NX)-WFLFR(NY,NX)/DENSI
      VOLP1(0,NY,NX)=AMAX1(0.0,VOLA1(0,NY,NX)-VOLW1(0,NY,NX)
     2-VOLI1(0,NY,NX))
      VOLWM(M+1,0,NY,NX)=VOLW1(0,NY,NX)
      VOLPM(M+1,0,NY,NX)=VOLP1(0,NY,NX)
      IF(VOLR(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWX(0,NY,NX)=AMAX1(0.0,VOLW1(0,NY,NX)/VOLR(NY,NX))
      THETIX(0,NY,NX)=AMAX1(0.0,VOLI1(0,NY,NX)/VOLR(NY,NX))
      THETPX(0,NY,NX)=AMAX1(0.0,VOLP1(0,NY,NX)/VOLR(NY,NX))
     2*AMAX1(0.0,(1.0-XVOLT(NY,NX)/VOLWD(NY,NX)))
      THETPM(M+1,0,NY,NX)=AMAX1(0.0,VOLPM(M+1,0,NY,NX)/VOLR(NY,NX))
     2*AMAX1(0.0,(1.0-XVOLT(NY,NX)/VOLWD(NY,NX)))
      ELSE
      THETWX(0,NY,NX)=0.0
      THETIX(0,NY,NX)=0.0
      THETPX(0,NY,NX)=1.0
      THETPM(M+1,0,NY,NX)=1.0
      ENDIF
      VHCPXX=VHCP1(0,NY,NX)
      TK0XX=TK1(0,NY,NX)
      ENGYR=VHCP1(0,NY,NX)*TK1(0,NY,NX)
      VHCP1(0,NY,NX)=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))
     2+4.19*(VOLW1(0,NY,NX)+VOLV1(0,NY,NX))
     2+1.9274*VOLI1(0,NY,NX)
      IF(VHCP1(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      TK1(0,NY,NX)=(ENGYR+HFLWRL(NY,NX)+HFLVR(NY,NX)+HFLFR(NY,NX)
     2+THQR1(NY,NX)+HFLXF(0,NY,NX))/VHCP1(0,NY,NX)
      ELSE
      TK1(0,NY,NX)=TK1(NUM(NY,NX),NY,NX)
      ENDIF
      TKSM(M+1,0,NY,NX)=TK1(0,NY,NX)
      TKGS(M+1,NY,NX)=FSNW(NY,NX)*TK0(1,NY,NX)+FSNX(NY,NX)
     2*(TK1(0,NY,NX)*CVRDW(NY,NX)+TK1(NUM(NY,NX),NY,NX)*BAREW(NY,NX))
C     IF(J.EQ.14.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,7754)'TKGS',I,J,NFZ,M,NX,NY,NUM(NY,NX),ICHKF
C    2,TKGS(M+1,NY,NX),FSNW(NY,NX),TK0(1,NY,NX),FSNX(NY,NX)
C    2,TK1(0,NY,NX),CVRDW(NY,NX),TK1(NUM(NY,NX),NY,NX),BAREW(NY,NX)
C     WRITE(*,7754)'VOLV1R',I,J,NFZ,M,NX,NY,NUM(NY,NX),ICHKF
C    2,VOLW1(0,NY,NX),VOLV1(0,NY,NX),VOLI1(0,NY,NX),VOLP1(0,NY,NX)
C    2,FLVRL(NY,NX),WFLVR(NY,NX),WFLFR(NY,NX),TQR1(NY,NX)
C    3,EVAPR(NY,NX),FLVRLG,FLV1,EVAPRV
C    4,XVOLW(NY,NX),XVOLI(NY,NX),XVOLT(NY,NX)
C    5,TVOLWI,VOLWRX(NY,NX),FLWRLG,FLWRLW,FLQR,FLYM,FLY1(NY,NX)
C    6,FLQ1(NY,NX),VOLP1(NUM(NY,NX),NY,NX),FLQRS,BARE(NY,NX) 
C    3,THETWX(0,NY,NX),THETIX(0,NY,NX),THETPX(0,NY,NX),VHCP1(0,NY,NX)
C    2,VOLW2(0,NY,NX),VOLV2(0,NY,NX),VOLI2(0,NY,NX),VOLP2(0,NY,NX)
C     WRITE(*,7754)'TKS0',I,J,NFZ,M,NX,NY,NUM(NY,NX),ICHKF
C    3,TK1(0,NY,NX),ENGYR,HFLWRL(NY,NX),HFLVR(NY,NX),HFLFR(NY,NX)
C    2,THQR1(NY,NX),HFLXF(0,NY,NX),VHCP1(0,NY,NX),VHCPXX
C    3,VHCPRX(NY,NX),ENGYR,CVRD(NY,NX)
C    2,HFLWRLW,HFLWRLG
C    3,HWFLYM,HFLXR,HFLWR(NY,NX)
C    4,HFLWRL(NY,NX)+HFLVR(NY,NX)+HFLFR(NY,NX)
C    2+THQR1(NY,NX)+HFLXF(0,NY,NX)
C    3,BARE(NY,NX),ORGC(0,NY,NX),ORGCC(0,NY,NX)
C    4,TK1(NUM(NY,NX),NY,NX),VHCP1(NUM(NY,NX),NY,NX) 
C    5,VHCM(NUM(NY,NX),NY,NX),VOLW1(NUM(NY,NX),NY,NX)
C    6,VOLV1(NUM(NY,NX),NY,NX),VOLWH1(NUM(NY,NX),NY,NX)
C    3,ALBZ(NY,NX),ORGC(0,NY,NX),ORGCC(0,NY,NX),VOLPM(M,0,NY,NX)
C    4,VOLPM(M,NUM(NY,NX),NY,NX),VOLI1(NUM(NY,NX),NY,NX)
C    4,TK1(NUM(NY,NX),NY,NX),HFLXF(NUM(NY,NX),NY,NX)
C    5,TKGS(M+1,NY,NX),FSNW(NY,NX),TK0(1,NY,NX),FSNX(NY,NX)
C    2,TK1(0,NY,NX),CVRDW(NY,NX),TK1(NUM(NY,NX),NY,NX),BAREW(NY,NX)
7754  FORMAT(A8,8I4,60E12.4)
C     ENDIF      
C
C     SOIL LAYER WATER, ICE AND TEMPERATURE
C
C     VOLW1,VOLV1,VOLI1=micropore water,vapor,ice volume (m3)
C     VOLWX1=micropore water volume behind wetting front (m3)
C     VOLWH1,VOLIH1=macropore water,ice volume (m3)
C     TFLWL,TFLWHL=net micropore,macropore water flux (m3 t-1)
C     TFLVL=net soil vapor flux (m3 t-1)
C     TWFLVL=total evaporation-condensation in micropores+macropores
C        (m3 t-1)
C     TWFLFL,TWFLFH=total freeze-thaw in micropores,macropores (m3 t-1)
C     FINHL=micropore-macropore flux (m3 t-1)
C     FLU1=subsurface water input (m3 t-1)
C     DENSI=ice density (Mg m-3)
C     VOLA1,VOLAH1=micropore,macropore volume (m3)
C     VOLP1,VOLPH1=micropore,macropore air volume (m3)
C     VOLWM,VOLWHM,VOLPM=micropore,macropore water volume,
C        air volume (m3)
C     FLPM=change in air volume for use in ‘trnsfr.f’ (m3 t-1)
C     THETWX,THETIX,THETPX,THETPY=bulk water,ice,air concentration,
C        air-filled porosity (m3 m-3)
C     THETPM=air concentration for use in ‘nitro.f’ and ‘trnsfr.f’ 
C        (m3 m-3)
C     FMAC,FGRD=macropore,micropore fraction
C     CNDH1=macropore hydraulic conductivity (m2 h MPa-1)
C     VHCP1,VHCM=volumetric heat capacities of total,solid volume 
C        (MJ K-1)
C     VHCP1A,VHCP1B=volumetric heat capacities of soil+micropore,
C        macropore (MJ K-1)
C
      DO 9785 L=NUM(NY,NX),NL(NY,NX)
      IF(VOLT(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      VOLW1(L,NY,NX)=VOLW1(L,NY,NX)+TFLWL(L,NY,NX)+TWFLVL(L,NY,NX) 
     2+TWFLFL(L,NY,NX)+FINHL(L,NY,NX)+FLU1(L,NY,NX)
      VOLV1(L,NY,NX)=VOLV1(L,NY,NX)+TFLVL(L,NY,NX)-TWFLVL(L,NY,NX) 
      VOLWX1(L,NY,NX)=VOLWX1(L,NY,NX)+TFLWLX(L,NY,NX) 
     2+FINHL(L,NY,NX)+TWFLFL(L,NY,NX)+FLU1(L,NY,NX) 
      VOLWX1(L,NY,NX)=AMIN1(VOLW1(L,NY,NX),VOLWX1(L,NY,NX))
      VOLI1(L,NY,NX)=VOLI1(L,NY,NX)-TWFLFL(L,NY,NX)/DENSI
      VOLWH1(L,NY,NX)=VOLWH1(L,NY,NX)+TFLWHL(L,NY,NX) 
     2-FINHL(L,NY,NX)+TWFLFH(L,NY,NX)
      VOLIH1(L,NY,NX)=VOLIH1(L,NY,NX)-TWFLFH(L,NY,NX)/DENSI
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
      VOLTX(L,NY,NX)=VOLY(L,NY,NX)+VOLAH1(L,NY,NX)
      IF(VOLTX(L,NY,NX).GT.ZEROS(NY,NX))THEN
      THETWX(L,NY,NX)=AMAX1(0.0,(VOLW1(L,NY,NX)+VOLWH1(L,NY,NX))
     2/VOLTX(L,NY,NX))
      THETIX(L,NY,NX)=AMAX1(0.0,(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
     2/VOLTX(L,NY,NX))
      THETPX(L,NY,NX)=AMAX1(0.0,(VOLP1(L,NY,NX)+VOLPH1(L,NY,NX))
     2/VOLTX(L,NY,NX))
      THETPM(M+1,L,NY,NX)=AMAX1(0.0,VOLPM(M+1,L,NY,NX)/VOLTX(L,NY,NX))
      ELSE
      THETWX(L,NY,NX)=0.0
      THETIX(L,NY,NX)=0.0
      THETPX(L,NY,NX)=1.0
      THETPM(M+1,L,NY,NX)=1.0
      ENDIF
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
      ENGY1=VHCP1(L,NY,NX)*TK1(L,NY,NX)
      VHCP1(L,NY,NX)=VHCM(L,NY,NX)
     2+4.19*(VOLW1(L,NY,NX)+VOLV1(L,NY,NX)+VOLWH1(L,NY,NX))
     2+1.9274*(VOLI1(L,NY,NX)+VOLIH1(L,NY,NX))
      VHCP1A(L,NY,NX)=VHCM(L,NY,NX)
     2+4.19*(VOLW1(L,NY,NX)+VOLV1(L,NY,NX))
     2+1.9274*VOLI1(L,NY,NX)
      VHCP1B(L,NY,NX)=4.19*VOLWH1(L,NY,NX)+1.9274*VOLIH1(L,NY,NX)
C
C     BEGIN ARTIFICIAL SOIL WARMING
C
C     THFLWL=THFLWL incremented for soil warming (MJ t-1)
C     TKSZ=temperature used to calculate additional heat flux 
C        for warming (K)
C     XNPHX=time step for heat fluxes from ‘wthr.f’ (h t-1)
C
C     IF(NX.EQ.3.AND.NY.EQ.2.AND.L.GT.NUM(NY,NX)
C    3.AND.L.LE.17.AND.I.GE.152.AND.I.LE.304)THEN
C     THFLWL(L,NY,NX)=THFLWL(L,NY,NX)
C    2+(TKSZ(I,J,L)-TK1(L,NY,NX))*VHCP1(L,NY,NX)*XNPHX
C     WRITE(*,3379)'WATSUB',I,J,NFZ,M,NX,NY,L,TKSZ(I,J,L)
C    2,TK1(L,NY,NX),VHCP1(L,NY,NX),THFLWL(L,NY,NX)
3379  FORMAT(A8,7I4,12E12.4)
C     ENDIF
C
C     END ARTIFICIAL SOIL WARMING
C
C     TK1=soil temperature (K)
C     THFLWL=soil conductive+convective heat flux (MJ t-1)
C     THFLFL=soil latent heat flux from freeze-thaw (MJ t-1)
C     THFLVL=soil latent heat flux from evaporation-condensation 
C        (MJ t-1)
C     HWFLU1=convective heat from subsurface water input (MJ t-1)
C     HFLXF=heat released by soil combustion (MJ t-1)
C     TKSM=soil temperature used in ‘trnsfr.f’(K)
C
      IF(VHCP1(L,NY,NX).GT.VHCPRX(NY,NX)
     2.AND.VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      TK1(L,NY,NX)=(ENGY1+THFLWL(L,NY,NX)+THFLFL(L,NY,NX)
     2+THFLVL(L,NY,NX)+HWFLU1(L,NY,NX)+HFLXF(L,NY,NX))/VHCP1(L,NY,NX)
      ELSE
      TK1(L,NY,NX)=TKS(L,NY,NX)
      ENDIF
      TKSM(M+1,L,NY,NX)=TK1(L,NY,NX)
      ELSE
      VOLWM(M+1,L,NY,NX)=0.0
      VOLWHM(M+1,L,NY,NX)=0.0
      VOLPM(M+1,L,NY,NX)=0.0
      FLPM(M,L,NY,NX)=VOLPM(M,L,NY,NX) 
      THETPM(M+1,L,NY,NX)=1.0
      IF(L.EQ.1)THEN
      TKSM(M+1,L,NY,NX)=TKQG(M,NY,NX)
      ELSE
      TKSM(M+1,L,NY,NX)=TK1(L-1,NY,NX)
      ENDIF
      ENDIF
C     IF(L.EQ.NU(NY,NX))THEN
C     WRITE(*,3377)'VOLV1',I,J,NFZ,M,NX,NY,L,N6X(NY,NX) 
C    2,VOLW1(L,NY,NX),VOLV1(L,NY,NX),VOLI1(L,NY,NX)
C    3,VOLA1(L,NY,NX),VOLP1(L,NY,NX),VOLPH1(L,NY,NX) 
C    4,VOLA1(L,NY,NX)-VOLW1(L,NY,NX)-VOLI1(L,NY,NX)
C    2,TFLWL(L,NY,NX),TFLVL(L,NY,NX),TWFLVL(L,NY,NX)
C    3,FLVL(3,L,NY,NX),FLVL(3,L+1,NY,NX),FLWLV,FLVLG 
C    2,TWFLFL(L,NY,NX),FINHL(L,NY,NX),FLU1(L,NY,NX)
C    5,VOLPM(M,L,NY,NX),VOLPM(M,L+1,NY,NX) 
C    5,PSISM1(L,NY,NX),THETPX(L,NY,NX)  
C    6,FLVL(3,L,NY,NX),FLVL(3,L+1,NY,NX)
C    7,FLWL(2,L,NY,NX),FLWL(2,L,NY+1,NX)
C    8,FLWL(1,L,NY,NX),FLWL(1,L,NY,NX+1)
C    6,FLWLX(3,L,NY,NX),FLWLX(3,L+1,NY,NX)
C    7,FLWLX(2,L,NY,NX),FLWLX(2,L,NY+1,NX)
C    8,FLWLX(1,L,NY,NX),FLWLX(1,L,NY,NX+1)
C    6,FLW(3,L,NY,NX),FLW(3,L+1,NY,NX) 
C    7,FLW(2,L,NY,NX),FLW(2,L,NY+1,NX)
C    8,FLW(1,L,NY,NX),FLW(1,L,NY,NX+1)
C    9,WFLFL(L,NY,NX),XWFLFL(L,NY,NX)
C    9,FLPM(M,L,NY,NX),FLSW(L,NY,NX)
C     WRITE(*,3377)'VOLWH1',I,J,NFZ,M,NX,NY,L,N6X(NY,NX)
C    2,VOLWH1(L,NY,NX)
C    2,TFLWHL(L,NY,NX),FINHL(L,NY,NX),VOLIH1(L,NY,NX) 
C    4,TWFLFH(L,NY,NX),TQR1(NY,NX),VOLPH1(L,NY,NX)
C    6,FLWHL(3,L,NY,NX),FLWHL(3,L+1,NY,NX)
C    7,FLWHL(2,L,NY,NX),FLWHL(2,L,NY+1,NX)
C    8,FLWHL(1,L,NY,NX),FLWHL(1,L,NY,NX+1)
C     WRITE(*,3377)'TKS1',I,J,NFZ,M,NX,NY,L,N6X(NY,NX),TK1(L,NY,NX)
C    2,THFLWL(L,NY,NX),THFLFL(L,NY,NX),THFLVL(L,NY,NX)
C    3,HWFLU1(L,NY,NX),HFLXF(L,NY,NX)
C    3,VHCP1(L,NY,NX)/AREA(3,L,NY,NX),VHCPRX(NY,NX)
C    3,HFLWL(3,L,NY,NX),HFLWL(3,N6X(NY,NX),NY,NX),HFLWLW,HFLWLG 
C    4,ENGY1,TKS(L,NY,NX),HFLW(3,L,NY,NX),HFLW(3,N6X(NY,NX),NY,NX)
C    5,VHCM(L,NY,NX),VOLW1(L,NY,NX),VOLV1(L,NY,NX),VOLWH1(L,NY,NX)
C    2,VOLI1(L,NY,NX),VOLIH1(L,NY,NX) 
3377  FORMAT(A8,8I4,40E12.4)
C     ENDIF
9785  CONTINUE
C
C     RESET SURFACE LAYER NUMBER AND TRANSFER ALL WATER TO SOIL 
C     SURFACE LAYER IF POND SURFACE LAYER IS LOST TO EVAPORATION 
C
C     BKDS=bulk density (0=pond)(Mg m-3)
C     NUM=new soil surface layer number after complete lake
C        evaporation
C     FLWNU,FLVNU,FLWHNU,HFLWNU=soil surface water,vapor flux, 
C        heat flux if pond surface disappears used in ‘redist.f’
C        (m3 t-1,MJ t-1)
C     
      IF(BKDS(NUM(NY,NX),NY,NX).LE.ZERO
     2.AND.VHCP1(NUM(NY,NX),NY,NX).LE.VHCPNX(NY,NX))THEN
      NUX=NUM(NY,NX)
      DO 9970 LL=NUX+1,NL(NY,NX)
      IF(VOLX(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
      NUM(NY,NX)=LL
      FLWNX(NY,NX)=FLW(3,NUM(NY,NX),NY,NX) 
      FLVNX(NY,NX)=FLV(3,NUM(NY,NX),NY,NX) 
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
      FLVNU(NY,NX)=FLV(3,N6X(NY,NX),NY,NX) 
      FLWXNU(NY,NX)=FLWX(3,N6X(NY,NX),NY,NX) 
      FLWHNU(NY,NX)=FLWH(3,N6X(NY,NX),NY,NX) 
      HFLWNU(NY,NX)=HFLW(3,N6X(NY,NX),NY,NX) 
      ELSE
      FLWNU(NY,NX)=FLWNX(NY,NX)
      FLVNU(NY,NX)=FLVNX(NY,NX)
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


