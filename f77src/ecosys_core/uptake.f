
      SUBROUTINE uptake(I,J,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE CALCULATES EXCHANGES OF ENERGY, C, N AND P
C     BETWEEN THE CANOPY AND THE ATMOSPHERE AND BETWEEN ROOTS AND THE SOIL
C
      include "parameters.h"
      include "blkc.h"
      include "blk1cp.h"
      include "blk1cr.h"
      include "blk1g.h"
      include "blk1n.h"
      include "blk1p.h"
      include "blk2a.h"
      include "blk2b.h"
      include "blk2c.h"
      include "blk3.h"
      include "blk5.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk9a.h"
      include "blk9b.h"
      include "blk9c.h"
      include "blk10.h"
      include "blk11a.h"
      include "blk11b.h"
      include "blk12a.h"
      include "blk12b.h"
      include "blk13a.h"
      include "blk13b.h"
      include "blk13c.h"
      include "blk13d.h"
      include "blk14.h"
      include "blk16.h"
      include "blk18a.h"
      include "blk18b.h"
      include "blk1u.h"
      DIMENSION PSIST1(JZ),PATH(2,JZ),RRADL(2,JZ)
     2,RSRT(2,JZ),ILYR(2,JZ),RSRG(2,JZ),RSR1(2,JZ),RSR2(2,JZ)
     3,RSSX(2,JZ),RSRS(2,JZ),WTRTG(JZ),FPQ(2,JZ,05),FPP(2,JZ,05)
     4,RACZ(JP,JY,JX),FRTDPX(JZ,05),RTARR(2,JZ)
     5,VOLPU(JZ),VOLWU(JZ)
C
C     MXN=max number of cycles in convergence soln for water uptake
C     DIFFX,DIFFY=acceptance criteria in convergence soln
C     FMN=min PFT:total population ratio
C     RACM,RACX=min,max canopy boundary layer resistance (h m-1)
C     RZ=surface resistance to evaporation (h m-1)
C     EMODW=wood modulus of elasticity (MPa)
C     DSTK,VSTK=stalk density (Mg m-3),specific volume (m3 g-1)
C     SNH3X=NH3 solubility at 25 oC (g m-3 water/(g m-3 air))
C     EMMC=canopy emissivity
C     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition on root,myco N,P uptake(g g-1)
C     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation (h-1)
C
      PARAMETER(MXN=200,DIFFX=1.0E-09,DIFFY=0.5E-02)
      PARAMETER (FMN=1.0E-06,RACM=0.00139,RACX=0.0278,RZ=0.0139)
      PARAMETER(DSTK=0.225,VSTK=1.0E-06/DSTK)
      PARAMETER(SNH3X=2.852E+02,EMMC=0.97,EMODW=50.0)
      PARAMETER(ZCKI=0.5E-01,PCKI=0.5E-02,ZPKI=ZCKI/PCKI
     2,PZKI=PCKI/ZCKI)
      PARAMETER(FEXUC=0.5E-03,FEXUN=1.0E-02,FEXUP=1.0E-02)
      REAL*4 RI,TKGO,TKSO
C     REAL*16 B,C

      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
C
C     RESET TOTAL UPTAKE ARRAYS
C
C     ARLFP,ARSTP=leaf,stalk areas
C
      ARLSC=0.0
      DO 9984 NZ=1,NP0(NY,NX)
C     TKC(NZ,NY,NX)=TKA(NY,NX)+DTKC(NZ,NY,NX)
C     TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-273.15
      ARLSC=ARLSC+ARLFP(NZ,NY,NX)+ARSTP(NZ,NY,NX)
      RAD1(NZ,NY,NX)=0.0
      EFLXC(NZ,NY,NX)=0.0
      SFLXC(NZ,NY,NX)=0.0
      HFLXC(NZ,NY,NX)=0.0
      THRM1(NZ,NY,NX)=0.0
      EP(NZ,NY,NX)=0.0
      EVAPC(NZ,NY,NX)=0.0
      UPOMC(NZ,NY,NX)=0.0
      UPOMN(NZ,NY,NX)=0.0
      UPOMP(NZ,NY,NX)=0.0
      UPNH4(NZ,NY,NX)=0.0
      UPNO3(NZ,NY,NX)=0.0
      UPH2P(NZ,NY,NX)=0.0
      UPH1P(NZ,NY,NX)=0.0
      UPNF(NZ,NY,NX)=0.0
C
C     RESET UPTAKE ARRAYS
C
      DO 9984 L=NU(NY,NX),NJ(NY,NX)
      DO 9984 N=1,MY(NZ,NY,NX)
      UPWTR(N,L,NZ,NY,NX)=0.0
      RCO2P(N,L,NZ,NY,NX)=0.0
      RUPOXP(N,L,NZ,NY,NX)=0.0
      RCO2S(N,L,NZ,NY,NX)=0.0
      RUPOXS(N,L,NZ,NY,NX)=0.0
      RUPCHS(N,L,NZ,NY,NX)=0.0
      RUPN2S(N,L,NZ,NY,NX)=0.0
      RUPN3S(N,L,NZ,NY,NX)=0.0
      RUPN3B(N,L,NZ,NY,NX)=0.0
      RUPHGS(N,L,NZ,NY,NX)=0.0
      RCOFLA(N,L,NZ,NY,NX)=0.0
      ROXFLA(N,L,NZ,NY,NX)=0.0
      RCHFLA(N,L,NZ,NY,NX)=0.0
      RN2FLA(N,L,NZ,NY,NX)=0.0
      RNHFLA(N,L,NZ,NY,NX)=0.0
      RHGFLA(N,L,NZ,NY,NX)=0.0
      RCODFA(N,L,NZ,NY,NX)=0.0
      ROXDFA(N,L,NZ,NY,NX)=0.0
      RCHDFA(N,L,NZ,NY,NX)=0.0
      RN2DFA(N,L,NZ,NY,NX)=0.0
      RNHDFA(N,L,NZ,NY,NX)=0.0
      RHGDFA(N,L,NZ,NY,NX)=0.0
9984  CONTINUE
C
C     PSIST1=total soil water potential PSIST adjusted for surf elevn
C     ALT=surface elevation
C     VOLWU,VOLWM=water volume available for uptake,total water volume
C     THETY,VOLX=hygroscopic SWC,soil volume
C     VOLPU=air volume
C     WTRTG=total biome root mass
C
      DO 9000 L=NU(NY,NX),NJ(NY,NX)
      PSIST1(L)=PSIST(L,NY,NX)-0.0098*ALT(NY,NX)
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLWU(L)=VOLWM(NPH,L,NY,NX)-THETY(L,NY,NX)*VOLY(L,NY,NX)
      VOLPU(L)=AMAX1(0.0,VOLA(L,NY,NX)-VOLW(L,NY,NX)-VOLI(L,NY,NX))
      ELSE
      VOLWU(L)=VOLWM(NPH,L,NY,NX)
      VOLPU(L)=0.0
      ENDIF
      WTRTG(L)=0.0
      DO 9005 NZ=1,NP(NY,NX)
      DO 9005 N=1,MY(NZ,NY,NX)
C     IF(IFLGC(NZ,NY,NX).EQ.1.AND.PP(NZ,NY,NX).GT.0.0)THEN
      WTRTG(L)=WTRTG(L)+AMAX1(0.0,WTRTD(N,L,NZ,NY,NX))
C     ENDIF
9005  CONTINUE
9000  CONTINUE
C
C     IF PLANT SPECIES EXISTS
C
      DO 9985 NZ=1,NP(NY,NX)
      OSTRN=0.0
      OSTRD=0.0
      IF(IFLGC(NZ,NY,NX).EQ.1.AND.PP(NZ,NY,NX).GT.0.0)THEN
C
C     APPLY CLUMPING FACTOR TO LEAF SURFACE AREA DEFINED BY
C     INCLINATION N, LAYER L, NODE K, BRANCH NB, SPECIES NZ,
C     N-S POSITION NY, E-W POSITION NX(AZIMUTH M ASSUMED UNIFORM)
C
      DO 500 NB=1,NBR(NZ,NY,NX)
      DO 550 K=1,25
C
C     NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
C
C     ARLF=leaf area
C     WSLF=leaf protein content
C     SURFX,SURF=unself-shaded,total leaf surface area
C     CFX=clumping factor from PFT file
C
      IF(ARLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WSLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      KLEAFX(NB,NZ,NY,NX)=K
      ENDIF
      DO 600 L=JC,1,-1
      DO 650 N=1,4
      SURFX(N,L,K,NB,NZ,NY,NX)=SURF(N,L,K,NB,NZ,NY,NX)*CFX(NZ,NY,NX)
650   CONTINUE
600   CONTINUE
550   CONTINUE
500   CONTINUE
C
C     CANOPY HEIGHT FROM HEIGHT OF MAXIMUM LEAF LAYER
C
C     DIFFUSIVE RESISTANCE OF OTHER TALLER CANOPIES TO HEAT AND VAPOR
C     TRANSFER OF CURRENT CANOPY ADDED TO BOUNDARY LAYER RESISTANCE
C     OF TALLEST CANOPY CALCULATED IN 'HOUR1'
C
C     IETYP=Koppen climate zone
C     ZC,ZT,ZR=PFT canopy,biome,surface roughness height
C     FRADP=fraction of radiation received by each PFT canopy
C     ALFZ=shape parameter for windspeed attenuation in canopy
C     RAB=biome canopy isothermal boundary layer resistance
C     RACZ,RAZ=additional,total PFT canopy isothermal blr
C
      IF(ARLFS(NZ,NY,NX).GT.0.0)THEN
      IF(IETYP(NY,NX).GE.0)THEN
      TFRADP=0.0
      DO 700 NZZ=1,NP(NY,NX)
      IF(ZC(NZZ,NY,NX).GT.ZC(NZ,NY,NX)+ZR(NY,NX))THEN
      TFRADP=TFRADP+FRADP(NZZ,NY,NX)
      ENDIF
700   CONTINUE
      ALFZ=2.0*TFRADP
      IF(RAB(NY,NX).GT.ZERO.AND.ZT(NY,NX).GT.ZERO
     2.AND.ALFZ.GT.ZERO)THEN
      RACZ(NZ,NY,NX)=AMIN1(RACX,AMAX1(0.0,ZT(NY,NX)*EXP(ALFZ)
     2/(ALFZ/RAB(NY,NX))*(EXP(-ALFZ*ZC(NZ,NY,NX)/ZT(NY,NX))
     3-EXP(-ALFZ*(ZD(NY,NX)+ZR(NY,NX))/ZT(NY,NX)))))
      ELSE
      RACZ(NZ,NY,NX)=0.0
      ENDIF
      ELSE
      RACZ(NZ,NY,NX)=0.0
      ENDIF
      ELSE
      RACZ(NZ,NY,NX)=RACX
      ENDIF
      RAZ(NZ,NY,NX)=RAB(NY,NX)+RACZ(NZ,NY,NX)
C     IF(NX.EQ.1.AND.NY.EQ.5.AND.NZ.EQ.3)THEN
C     WRITE(*,4411)'RAC',I,J,NZ,RACZ(NZ,NY,NX),RAB(NY,NX),RAZ(NZ,NY,NX)
C    2,ZT(NY,NX),ZD(NY,NX),ZR(NY,NX),ZC(NZ,NY,NX),ALFZ
C    3,TFRADP,RACX,EXP(-ALFZ*ZC(NZ,NY,NX)/ZT(NY,NX))
C    4,EXP(-ALFZ*(ZD(NY,NX)+ZR(NY,NX))/ZT(NY,NX)),RADP(NZ,NY,NX)
C    5,FRADP(NZ,NY,NX),DPTHS(NY,NX)
4411  FORMAT(A8,3I4,20E12.4)
C     ENDIF
C
C     INITIALIZE CANOPY TEMPERATURE WITH CURRENT AIR TEMPERATURE AND
C     LAST HOUR'S CANOPY-AIR TEMPERATURE DIFFERENCE, AND CALL A
C     SUBROUTINE TO CALCULATE MINIMUM CANOPY STOMATAL RESISTANCE
C     FOR SUBSEQUENT USE IN ENERGY EXCHANGE CALCULATIONS
C
C     TKCZ=initial estimate of canopy temperature TKC
C     TKA=current air temperature
C     DTKC=TKC-TKA from previous hour
C     STOMATE=solve for minimum canopy stomatal resistance
C
      TKCZ(NZ,NY,NX)=TKA(NY,NX)+DTKC(NZ,NY,NX)
      CALL STOMATE(I,J,NZ,NY,NX)
C
C     CALCULATE VARIABLES USED IN ROOT UPTAKE OF WATER AND NUTRIENTS
C
C     RTDPZ,RTDP1=primary root depth
C     FRTDPX=fraction of each soil layer with primary root
C     DLYR=layer thickness
C     CDPTHZ=depth from soil surface to layer bottom
C     SDPTH=seeding depth
C     HTCTL=hypocotyledon height
C     FPQ=PFT fraction of biome root mass
C     RTDNP,RTLGP=root length density,root length per plant
C     RRADL=root radius
C     PORT=root porosity
C     FMPR=micropore fraction excluding macropore,rock
C     RTVLW=root aqueous volume
C     PP=plant population
C     PATH=path length of water and nutrient uptake
C     RTARR=root surface area/radius for uptake,diffusivity
C
      DO 2000 N=1,MY(NZ,NY,NX)
      DO 2000 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(N.EQ.1)THEN
      RTDPZ=0.0
      DO 2005 NR=1,NRT(NZ,NY,NX)
      RTDPZ=AMAX1(RTDPZ,RTDP1(1,NR,NZ,NY,NX))
2005  CONTINUE
      IF(L.EQ.NU(NY,NX))THEN
      FRTDPX(L,NZ)=1.0
      ELSE
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
      RTDPX=AMAX1(0.0,RTDPZ-CDPTHZ(L-1,NY,NX))
      RTDPX=AMAX1(0.0,AMIN1(DLYR(3,L,NY,NX),RTDPX)
     2-AMAX1(0.0,SDPTH(NZ,NY,NX)-CDPTHZ(L-1,NY,NX)-HTCTL(NZ,NY,NX)))
      FRTDPX(L,NZ)=RTDPX/DLYR(3,L,NY,NX)
      ELSE
      FRTDPX(L,NZ)=0.0
      ENDIF
      ENDIF
C     IF(NZ.EQ.1.OR.NZ.EQ.2)THEN
C     WRITE(*,4413)'FRTDPX',I,J,NZ,L,N,FRTDPX(L,NZ)
C    2,RTDPX,DLYR(3,L,NY,NX),RTDPZ,RTDP1(1,1,NZ,NY,NX)
C    3,SDPTH(NZ,NY,NX),CDPTHZ(L-1,NY,NX),HTCTL(NZ,NY,NX)
C    4,WTRTD(1,L,NZ,NY,NX),WTRTG(L),FPQ(N,L,NZ)
C     ENDIF
      ENDIF
      IF(WTRTG(L).GT.ZEROS(NY,NX))THEN
      FPQ(N,L,NZ)=AMAX1(0.0,WTRTD(N,L,NZ,NY,NX))/WTRTG(L)
      ELSE
      FPQ(N,L,NZ)=1.0
      ENDIF
      FPP(N,L,NZ)=FMN*FPQ(N,L,NZ)
      IF(RTDNP(N,L,NZ,NY,NX).GT.ZERO
     2.AND.FRTDPX(L,NZ).GT.ZERO)THEN
      RRADL(N,L)=AMAX1(RRAD2X(N,NZ,NY,NX),SQRT((RTVLW(N,L,NZ,NY,NX)
     2/(1.0-PORT(N,NZ,NY,NX)))/(3.1416*PP(NZ,NY,NX)
     3*RTLGP(N,L,NZ,NY,NX))))
      PATH(N,L)=AMAX1(1.001*RRADL(N,L)
     2,1.0/(SQRT(3.1416*(RTDNP(N,L,NZ,NY,NX)/FRTDPX(L,NZ))
     3/FMPR(L,NY,NX))))
      RTARR(N,L)=6.283*RTLGP(N,L,NZ,NY,NX)/FRTDPX(L,NZ)
      ELSE
      RRADL(N,L)=RRAD2M(N,NZ,NY,NX)
      PATH(N,L)=1.001*RRADL(N,L)
      RTARR(N,L)=6.283*RTLGP(N,L,NZ,NY,NX)
      ENDIF
C     IF(NZ.EQ.1.OR.NZ.EQ.2)THEN
C     WRITE(*,4413)'RTAR',I,J,NX,NY,NZ,L,N
C    2,RTARR(N,L),RRADL(N,L),PATH(N,L),WTRTD(N,L,NZ,NY,NX)
C    2,RTLGP(N,L,NZ,NY,NX),RTDNP(N,L,NZ,NY,NX)
C    2,PP(NZ,NY,NX),RTDNP(N,L,NZ,NY,NX)*PP(NZ,NY,NX)
C    2/AREA(3,L,NY,NX),FRTDPX(L,NZ),RTLG2X(N,NZ,NY,NX)
C    2,FPQ(N,L,NZ),WTRTG(L)
C    3,RRADL(N,L),RRAD2X(N,NZ,NY,NX),RTVLW(N,L,NZ,NY,NX)
C    2,PORT(N,NZ,NY,NX),PP(NZ,NY,NX),RTLGP(N,L,NZ,NY,NX)
4413  FORMAT(A8,7I4,30E12.4)
C     ENDIF
2000  CONTINUE
C
C     CALCULATE CANOPY WATER STATUS FROM CONVERGENCE SOLUTION FOR
C     TRANSPIRATION - ROOT WATER UPTAKE = CHANGE IN CANOPY WATER CONTENT
C
C     IF(NX.EQ.4.AND.NY.EQ.5)THEN
C     WRITE(*,2123)'START',I,J,NX,NY,NZ
C    2,IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX),ARLFS(NZ,NY,NX),PP(NZ,NY,NX)
C    3,FRADP(NZ,NY,NX),RTDP1(1,1,NZ,NY,NX),RADC(NZ,NY,NX)
C    4,RADP(NZ,NY,NX),FRTDPX(L,NZ),RTDPX,RAZ(NZ,NY,NX)
2123  FORMAT(A8,6I4,20E12.4)
C     ENDIF
      IF((IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0)
     3.AND.(ARLFS(NZ,NY,NX).GT.ZEROL(NZ,NY,NX)
     4.AND.FRADP(NZ,NY,NX).GT.0.0)
     5.AND.(RTDP1(1,1,NZ,NY,NX).GT.SDPTH(NZ,NY,NX)+CDPTHZ(0,NY,NX)))THEN
C
C     GRAVIMETRIC WATER POTENTIAL FROM CANOPY HEIGHT
C
C     HTSTZ=canopy height for water uptake
C     PSILH=gravimetric water potential at HTSTZ
C     FRADW=conducting elements of stalk relative to those of primary root
C     PSILT=canopy total water potential
C     EMODW=wood modulus of elasticity (MPa)
C
      CNDT=0.0
      HTSTZ(NZ,NY,NX)=0.80*ZC(NZ,NY,NX)
      PSILH=-0.0098*HTSTZ(NZ,NY,NX)
      FRADW=1.0E+04*(AMAX1(0.5,1.0+PSILT(NZ,NY,NX)/EMODW))**4
C
C     SOIL AND ROOT HYDRAULIC RESISTANCES TO ROOT WATER UPTAKE
C
C     VOLX,VOLWM,THETW=soil,water volume,content
C     RTDNP,RTLGP=root length density,root length per plant
C     CNDU=soil hydraulic conductivity for root uptake
C     RTN1,RTNL=number of root,myco primary,secondary axes
C     ILYR:1=rooted,0=not rooted
C     N:1=root,2=mycorrhizae
C
      DO 3880 N=1,MY(NZ,NY,NX)
      DO 3880 L=NU(NY,NX),NI(NZ,NY,NX)
C     IF(NZ.EQ.2)THEN
C     WRITE(*,2124)'ILYR',I,J,NX,NY,NZ,L,N,RTDNP(N,L,NZ,NY,NX)
C    2,CNDU(L,NY,NX),RTN1(1,L,NZ,NY,NX),RTNL(N,L,NZ,NY,NX)
C    3,THETW(L,NY,NX),ZEROP(NZ,NY,NX)
2124  FORMAT(A8,7I4,20E12.4)
C      ENDIF
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX)
     2.AND.VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX)
     2.AND.RTDNP(N,L,NZ,NY,NX).GT.ZERO
     2.AND.CNDU(L,NY,NX).GT.ZERO
     3.AND.RTN1(1,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     4.AND.RTNL(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     5.AND.THETW(L,NY,NX).GT.ZERO)THEN
      ILYR(N,L)=1
C
C     SOIL HYDRAULIC RESISTANCE FROM RADIAL UPTAKE GEOMETRY
C     AND SOIL HYDRAULIC CONDUCTIVITY
C
C     RSSX=soil hydraulic resistance
C     PP=plant population
C     PATH=path length of water and nutrient uptake
C     RRADL,RTARR=root radius,surface/radius area
C
      RSSL=(LOG(PATH(N,L)/RRADL(N,L))/RTARR(N,L))/PP(NZ,NY,NX)
      RSSX(N,L)=RSSL/CNDU(L,NY,NX)
C
C     RADIAL ROOT RESISTANCE FROM ROOT AREA AND RADIAL RESISTIVITY
C     ENTERED IN 'READQ'
C
C     RRAD2=secondary root radius
C     RTLGP=root length per plant
C     RSRG=radial resistance
C     RSRR=radial resistivity from PFT file
C     VOLA,VOLWM=soil micropore,water volume
C
      RTAR2=6.283*RRAD2(N,L,NZ,NY,NX)*RTLGP(N,L,NZ,NY,NX)*PP(NZ,NY,NX)
      RSRG(N,L)=RSRR(N,NZ,NY,NX)/RTAR2
     2*VOLA(L,NY,NX)/VOLWM(NPH,L,NY,NX)
C
C     ROOT AXIAL RESISTANCE FROM RADII AND LENGTHS OF PRIMARY AND
C     SECONDARY ROOTS AND FROM AXIAL RESISTIVITY ENTERED IN 'READQ'
C
C     FRAD1,FRAD2=primary,secondary root radius relative to maximum
C     secondary radius from PFT file RRAD2M at which RSRA is defined
C     RRAD1,RRAD2=primary,secondary root radius
C     RSRA=axial resistivity from PFT file
C     DPTHZ=depth of primary root from surface
C     RSR1,RSR2=axial resistance of primary,secondary roots
C     RTLGA=average secondary root length
C     RTN1,RTNL=number of primary,secondary axes
C
      FRAD1=(RRAD1(N,L,NZ,NY,NX)/RRAD2M(N,NZ,NY,NX))**4
      RSR1(N,L)=RSRA(N,NZ,NY,NX)*DPTHZ(L,NY,NX)
     2/(FRAD1*RTN1(1,L,NZ,NY,NX))
     3+RSRA(1,NZ,NY,NX)*HTSTZ(NZ,NY,NX)
     4/(FRADW*RTN1(1,L,NZ,NY,NX))
      FRAD2=(RRAD2(N,L,NZ,NY,NX)/RRAD2M(N,NZ,NY,NX))**4
      RSR2(N,L)=RSRA(N,NZ,NY,NX)*RTLGA(N,L,NZ,NY,NX)
     2/(FRAD2*RTNL(N,L,NZ,NY,NX))
      ELSE
      ILYR(N,L)=0
      ENDIF
3880  CONTINUE
      DO 3890 N=1,MY(NZ,NY,NX)
      DO 3890 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(ILYR(N,L).EQ.1)THEN
C
C     TOTAL ROOT RESISTANCE = SOIL + RADIAL + AXIAL
C
C     RSRT=root radial+axial resistance
C     RSRS=total soil+root resistance
C     CNDT=total soil+root conductance for all layers
C
      RSRT(N,L)=RSRG(N,L)+RSR1(N,L)+RSR2(N,L)
      RSRS(N,L)=RSSX(N,L)+RSRT(N,L)
      CNDT=CNDT+1.0/RSRS(N,L)
C     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,8855)'RSRT',I,J,NX,NY,NZ,L,N,RSRT(N,L),RSRG(N,L)
C    2,RSR1(N,L),RSR2(N,L),RSSX(N,L),RSRS(N,L),DPTHZ(L,NY,NX)
C    3,HTSTZ(NZ,NY,NX),RSRA(1,NZ,NY,NX)*HTSTZ(NZ,NY,NX)
C    4/(FRADW*AMAX1(PP(NZ,NY,NX),RTN1(1,L,NZ,NY,NX)))
C    4,RTNL(N,L,NZ,NY,NX),RTLGP(N,L,NZ,NY,NX)
C    7,RTLGA(N,L,NZ,NY,NX),FRAD1,PP(NZ,NY,NX)
C    8,RTN1(1,L,NZ,NY,NX),FRADW,RTNL(N,L,NZ,NY,NX),CNDT
8855  FORMAT(A8,7I4,30E14.6)
C     ENDIF
      ENDIF
3890  CONTINUE
      ICHK=0
      PSIL2=0.0
      EPX=0.0
      UPRTX=0.0
      VOLWPX=0.0
C
C     INITIALIZE CANOPY WATER POTENTIAL, OTHER VARIABLES USED IN ENERGY
C     BALANCE THAT DON'T NEED TO BE RECALCULATED DURING CONVERGENCE
C
C     PSILT=initial estimate of total canopy water potential
C     FLWC,HFLWC1=convective water,heat flux from precip to canopy
C     FTHRM,FDTHS=LW emitted,absorbed by canopy
C     FPC=fraction of ground area AREA occupied by PFT
C     CCPOLT=total nonstructural canopy C,N,P concentration
C     CCPOLP,CZPOLP,CPPOLP=nonstructural C,N,P concn in canopy
C     OSWT=molar mass of CCPOLT
C     TKCX=intermediate estimate of TKC used in convergence solution
C     WVPLT=leaf+petiole+stalk mass
C     VSTK=specific stalk volume
C     VHCPX=canopy heat capacity
C     VOLWP,VOLWC=water volume in canopy,on canopy surfaces
C     RAZ=canopy isothermal boundary later resistance
C
      PSILT(NZ,NY,NX)=AMIN1(-1.0E-06,0.667*PSILT(NZ,NY,NX))
      EP(NZ,NY,NX)=0.0
      EVAPC(NZ,NY,NX)=0.0
      UPRT=0.0
      HFLWC1=FLWC(NZ,NY,NX)*4.19*TKA(NY,NX)
      FTHRM=EMMC*2.04E-10*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      FDTHS=(THS(NY,NX)+THRMGX(NY,NX))*FRADP(NZ,NY,NX)
      IF(ARLSS(NY,NX).GT.ZEROS(NY,NX))THEN
      FPC=ARLFS(NZ,NY,NX)/ARLSS(NY,NX)*AMIN1(1.0,0.5*ARLFC(NY,NX)
     2/AREA(3,NU(NY,NX),NY,NX))
      ELSEIF(PPT(NY,NX).GT.ZEROS(NY,NX))THEN
      FPC=PP(NZ,NY,NX)/PPT(NY,NX)
      ELSE
      FPC=1.0/NP(NY,NX)
      ENDIF
      PAREX=FPC*AREA(3,NU(NY,NX),NY,NX)
      PARHX=FPC*AREA(3,NU(NY,NX),NY,NX)*1.25E-03
      CCPOLT=CCPOLP(NZ,NY,NX)+CZPOLP(NZ,NY,NX)+CPPOLP(NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      TKCX=TKC(NZ,NY,NX)
      WVPLT=AMAX1(0.0,WTLS(NZ,NY,NX)+WVSTK(NZ,NY,NX))
      VHCPX=4.19*(WVPLT*VSTK+VOLWC(NZ,NY,NX)+VOLWP(NZ,NY,NX))
      RA1=RAZ(NZ,NY,NX)
      IC=0
      XC=0.5
C
C     CONVERGENCE SOLUTION
C
      DO 4000 NN=1,MXN
C
C     NET RADIATION FROM ABSORBED SW AND NET LW
C
C     THRM1=LW emitted by canopy
C     DTHS1=net LW absorbed by canopy
C     RADC=total SW absorbed by canopy
C     RAD1=net SW+LW absorbed by canopy
C
      TKC1=TKCZ(NZ,NY,NX)
      THRM1(NZ,NY,NX)=FTHRM*TKC1**4
      DTHS1=FDTHS-THRM1(NZ,NY,NX)*2.0
      RAD1(NZ,NY,NX)=RADC(NZ,NY,NX)+DTHS1
C
C     BOUNDARY LAYER RESISTANCE FROM RICHARDSON NUMBER
C
C     RI=Ricardson�s number
C     RA=canopy boundary layer resistance
C     PAREC,PARHC=canopy latent,sensible heat conductance
C
      RI=AMAX1(-0.3,AMIN1(0.075,RIB(NY,NX)*(TKA(NY,NX)-TKC1)))
C     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,4443)'RI',I,J,NX,NY,NZ,NN,RI,RIB(NY,NX)
C    2,TKA(NY,NX),TKC1
4443  FORMAT(A8,6I4,12E24.16)
C     ENDIF
      RA(NZ,NY,NX)=AMAX1(RACM,0.9*RA1,AMIN1(1.1*RA1
     2,RAZ(NZ,NY,NX)/(1.0-10.0*RI)))
      RA1=RA(NZ,NY,NX)
      PAREC=PAREX/RA(NZ,NY,NX)
      PARHC=PARHX/RA(NZ,NY,NX)
C
C     CANOPY WATER AND OSMOTIC POTENTIALS
C
C     PSILT=canopy total water potential
C     FDMP=dry matter content
C     OSMO=osmotic potential at PSILT=0 from PFT file
C     PSILO,PSILG=canopy osmotic,turgor water potential
C
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
      PSILO(NZ,NY,NX)=FDMP/0.16*OSMO(NZ,NY,NX)
     2-8.3143*TKC1*FDMP*CCPOLT/OSWT
      PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
C
C     CANOPY STOMATAL RESISTANCE
C
C     RCS=shape parameter for RC vs PSILG from PFT file
C     RC=canopy stomatal resistance
C     RSMN=minimum RC at PSILT=0 from stomate.f
C     RSMX=cuticular resistance from PFT file
C
      WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      RC(NZ,NY,NX)=RSMN(NZ,NY,NX)+(RSMH(NZ,NY,NX)-RSMN(NZ,NY,NX))*WFNC
C
C     CANOPY VAPOR PRESSURE AND EVAPORATION OF INTERCEPTED WATER
C     OR TRANSPIRATION OF UPTAKEN WATER
C
C     VPC,VPA=vapor pressure inside canopy, in atmosphere
C     TKC1=canopy temperature
C     PSILT=canopy total water potential
C     EX=canopy-atmosphere water flux
C     RA,RZ=canopy boundary layer,surface resistance
C     EVAPC=water flux to,from canopy surfaces
C     VOLWC=water volume on canopy surfaces
C     EP=water flux to,from inside canopy
C     EFLXC=canopy latent heat flux
C     VFLXC=convective heat flux from EFLXC
C     VAP=latent heat of evaporation
C
      VPC=2.173E-03/TKC1
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKC1))
     3*EXP(18.0*PSILT(NZ,NY,NX)/(8.3143*TKC1))
      EX=PAREC*(VPA(NY,NX)-VPC)
      IF(EX.GT.0.0)THEN
      EVAPC(NZ,NY,NX)=EX*RA(NZ,NY,NX)/(RA(NZ,NY,NX)+RZ)
      EX=0.0
      ELSEIF(EX.LE.0.0.AND.VOLWC(NZ,NY,NX).GT.0.0)THEN
      EVAPC(NZ,NY,NX)=AMAX1(EX*RA(NZ,NY,NX)/(RA(NZ,NY,NX)+RZ)
     2,-VOLWC(NZ,NY,NX))
      EX=EX-EVAPC(NZ,NY,NX)
      ENDIF
      EP(NZ,NY,NX)=EX*RA(NZ,NY,NX)/(RA(NZ,NY,NX)+RC(NZ,NY,NX))
      EFLXC(NZ,NY,NX)=(EP(NZ,NY,NX)+EVAPC(NZ,NY,NX))*VAP
      VFLXC=EVAPC(NZ,NY,NX)*4.19*TKC1
C
C     SENSIBLE + STORAGE HEAT FROM RN, LE AND CONVECTIVE HEAT FLUXES
C
C     HFLXS=initial estimate of sensible+storage heat flux
C     HFLWC1=convective heat flux from precip to canopy
C
      HFLXS=RAD1(NZ,NY,NX)+EFLXC(NZ,NY,NX)+VFLXC+HFLWC1
C
C     SOLVE FOR CANOPY TEMPERATURE CAUSED BY SENSIBLE + STORAGE HEAT
C
C     VHCPC=canopy heat capacity
C     TKCY=equilibrium canopy temperature for HFLXS
C
      VHCPC(NZ,NY,NX)=VHCPX+4.19*(EVAPC(NZ,NY,NX)+FLWC(NZ,NY,NX))
      TKCY=(TKCX*VHCPX+TKA(NY,NX)*PARHC+HFLXS)/(VHCPC(NZ,NY,NX)+PARHC)
      TKCY=AMIN1(TKA(NY,NX)+10.0,AMAX1(TKA(NY,NX)-10.0,TKCY))
C
C     RESET CANOPY TEMPERATURE FOR NEXT ITERATION
C
C     XC,IC=magnitude,direction of change in canopy temp for next cycle
C
      IF((IC.EQ.0.AND.TKCY.GT.TKC1).OR.(IC.EQ.1.AND.TKCY.LT.TKC1))THEN
      XC=0.5*XC
      ENDIF
      TKCZ(NZ,NY,NX)=TKC1+0.1*(TKCY-TKC1)
      IF(TKCY.GT.TKC1)THEN
      IC=1
      ELSE
      IC=0
      ENDIF
C     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,4444)'TKZ',I,J,NX,NY,NZ,NN,XC,TKC1,TKCY,TKCZ(NZ,NY,NX)
C    2,TKA(NY,NX),TKCX,VHCPX,PARHC,HFLXS,VHCPC(NZ,NY,NX),WVPLT,EX
C    2,FLWC(NZ,NY,NX),VOLWC(NZ,NY,NX),VOLWP(NZ,NY,NX),EVAPC(NZ,NY,NX)
C    2,RAD1(NZ,NY,NX),EFLXC(NZ,NY,NX),RA(NZ,NY,NX),RC(NZ,NY,NX)
C    2,EP(NZ,NY,NX),HFLXS,VFLXC,HFLWC1,RADC(NZ,NY,NX),FRADP(NZ,NY,NX)
C    3,THS(NY,NX),THRMGX(NY,NX)
C    2,RSMN(NZ,NY,NX),CCPOLT,OSWT,CCPOLP(NZ,NY,NX),CPOOLP(NZ,NY,NX)
C    4,DCO2(NZ,NY,NX),AREA(3,NU(NY,NX),NY,NX),WTLS(NZ,NY,NX)
C    2,PSILT(NZ,NY,NX),PSILG(NZ,NY,NX),RACZ(NZ,NY,NX),RAZ(NZ,NY,NX),RI
C    3,RIB(NY,NX),RA1,ARLFV(1,NZ,NY,NX),ARSTV(1,NZ,NY,NX)
C    4,WVPLT,VSTK
4444  FORMAT(A8,6I4,60E16.8)
C     ENDIF
C
C     IF CONVERGENCE CRITERION IS MET OR ON EVERY TENTH ITERATION,
C     PROCEED TO WATER BALANCE
C
C     PSILC=canopy water potential adjusted for canopy height
C
      IF(ABS(TKCY-TKC1).LT.0.05.OR.(NN/10)*10.EQ.NN)THEN
      UPRT=0.0
      PSILC=PSILT(NZ,NY,NX)-PSILH
C
C     ROOT WATER UPTAKE FROM SOIL-CANOPY WATER POTENTIALS,
C     SOIL + ROOT HYDRAULIC RESISTANCES
C
C     ILYR=rooted layer flag
C     UPWTR=root water uptake from soil layer
C     VOLWU,VOLPU=water volume available for uptake,air volume
C     FPQ=PFT fraction of biome root mass
C     PSILC=canopy water potential adjusted for canopy height
C     PSIST1=total soil water potential PSIST adjusted for surf elevn
C     RSRS=total soil+root resistance
C     UPRT=total water uptake from soil profile
C
      DO 4200 N=1,MY(NZ,NY,NX)
      DO 4200 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(ILYR(N,L).EQ.1)THEN
      UPWTR(N,L,NZ,NY,NX)=AMAX1(AMIN1(0.0,-VOLWU(L)*FPQ(N,L,NZ))
     2,AMIN1((PSILC-PSIST1(L))/RSRS(N,L),VOLPU(L)*FPQ(N,L,NZ)))
      IF(UPWTR(N,L,NZ,NY,NX).GT.0.0)THEN
      UPWTR(N,L,NZ,NY,NX)=0.1*UPWTR(N,L,NZ,NY,NX)
      ENDIF
      UPRT=UPRT+UPWTR(N,L,NZ,NY,NX)
      ELSE
      UPWTR(N,L,NZ,NY,NX)=0.0
      ENDIF
C     IF(NZ.EQ.2)THEN
C     WRITE(*,6565)'UPRT',I,J,NX,NY,NZ,L,N,NN,ILYR(N,L)
C    2,UPRT,UPWTR(N,L,NZ,NY,NX)
C    2,PSILC,PSIST1(L),PSISM(L,NY,NX),RSRS(N,L),RSSX(N,L)
C    2,RSRT(N,L),RSRG(N,L),RSR1(N,L),RSR2(N,L),PSILH,RTAR2
C    3,RSRR(N,NZ,NY,NX),VOLA(L,NY,NX),VOLWM(NPH,L,NY,NX)
C    4,VOLWU(L),THETY(L,NY,NX),VOLY(L,NY,NX),FPQ(N,L,NZ)
6565  FORMAT(A8,9I4,30E12.4)
C     ENDIF
4200  CONTINUE
C
C     TEST TRANSPIRATION - ROOT WATER UPTAKE VS. CHANGE IN CANOPY
C     WATER STORAGE
C
C     VOLWPZ,VOLWP=canopy water content
C     DIFFZ,DIFFU=change in canopy water content,transpiration-uptake
C     DIFF=normalized difference between DIFFZ and DIFFU
C     5.0E-03=acceptance criterion for DIFF
C     RSSZ=change in canopy water potl vs change in canopy water cnt
C     RSSU=change in canopy water potl vs change in transpiration
C
      VOLWPZ=1.0E-06*WVPLT/FDMP
      DIFFZ=VOLWPZ-VOLWP(NZ,NY,NX)
      DIFFU=EP(NZ,NY,NX)-UPRT
      IF(UPRT.NE.0.0)THEN
      DIFF=ABS((DIFFU-DIFFZ)/UPRT)
      ELSE
      DIFF=ABS((DIFFU-DIFFZ)/VOLWPZ)
      ENDIF
      IF(DIFF.LT.5.0E-03)GO TO 4250
      IF(ABS(VOLWPZ-VOLWPX).GT.ZEROP(NZ,NY,NX))THEN
      RSSZ=ABS((PSILT(NZ,NY,NX)-PSIL2)/(VOLWPZ-VOLWPX))
      ELSEIF(CNDT.GT.ZEROP(NZ,NY,NX))THEN
      RSSZ=1.0/CNDT
      ELSE
      RSSZ=ZEROL(NZ,NY,NX)
      ENDIF
      IF(ABS(EP(NZ,NY,NX)-EPX).GT.ZEROP(NZ,NY,NX))THEN
      RSSUX=ABS((PSILT(NZ,NY,NX)-PSIL2)/(EP(NZ,NY,NX)-EPX))
      IF(CNDT.GT.ZEROP(NZ,NY,NX))THEN
      RSSU=AMIN1(1.0/CNDT,RSSUX)
      ELSE
      RSSU=RSSUX
      ENDIF
      ELSEIF(ABS(UPRT-UPRTX).GT.ZEROP(NZ,NY,NX))THEN
      RSSUX=ABS((PSILT(NZ,NY,NX)-PSIL2)/(UPRT-UPRTX))
      IF(CNDT.GT.ZEROP(NZ,NY,NX))THEN
      RSSU=AMIN1(1.0/CNDT,RSSUX)
      ELSE
      RSSU=RSSUX
      ENDIF
      ELSEIF(CNDT.GT.ZEROP(NZ,NY,NX))THEN
      RSSU=1.0/CNDT
      ELSE
      RSSU=ZEROL(NZ,NY,NX)
      ENDIF
C
C     CHANGE IN CANOPY WATER POTENTIAL REQUIRED TO BRING AGREEMENT
C     BETWEEN TRANSPIRATION - ROOT WATER UPTAKE AND CHANGE IN CANOPY
C     WATER STORAGE
C
C     DPSI=change in PSILT for next convergence cycle
C     1.0E-03=acceptance criterion for DPSI
C
      DPSI=AMIN1(AMIN1(RSSZ,RSSU)*(DIFFU-DIFFZ),ABS(PSILT(NZ,NY,NX)))
C     IF(NX.EQ.1.AND.NY.EQ.5.AND.NZ.EQ.3)THEN
C     WRITE(*,2222)'PSI',I,J,NX,NY,NZ,NN,PSILT(NZ,NY,NX),PSIL2,DPSI
C    2,RSSUX,RSSU,RSSZ,1.0/CNDT,UPRT,UPRTX,EP(NZ,NY,NX),EPX,EX
C    3,EVAPC(NZ,NY,NX),RC(NZ,NY,NX),RA(NZ,NY,NX),FRADP(NZ,NY,NX)
C    3,PAREC,VPA(NY,NX),VPC,TKA(NY,NX),TKC1,VOLWP(NZ,NY,NX)
C    4,VOLWPZ,VOLWPX,WVPLT,DIFF,WFNC,PSILG(NZ,NY,NX)
C    5,FDMP,CCPOLT,OSWT,RAZ(NZ,NY,NX),RI,RIB(NY,NX)
C    5,((UPWTR(N,L,NZ,NY,NX),L=1,8),N=1,1)
C    6,((RSRS(N,L),L=1,8),N=1,1)
C    7,(PSIST1(L),L=1,8)
C    4,DIFFZ,DIFFU,DIFF
2222  FORMAT(A8,6I4,80E12.4)
C     ENDIF
C
C     IF CONVERGENCE CRITERION IS MET THEN FINISH,
C     OTHERWISE START NEXT ITERATION WITH CANOPY WATER POTENTIAL
C     TRANSPIRATION, UPTAKE AND WATER CONTENT FROM CURRENT ITERATION
C
      IF((NN.GE.30.AND.ABS(DPSI).LT.1.0E-03).OR.NN.GE.MXN)GO TO 4250
      PSIL2=PSILT(NZ,NY,NX)
      EPX=EP(NZ,NY,NX)
      UPRTX=UPRT
      VOLWPX=VOLWPZ
      PSILT(NZ,NY,NX)=AMIN1(0.0,PSILT(NZ,NY,NX)+0.5*DPSI)
      XC=0.50
      GO TO 4000
C
C     RESET MIN STOMATAL RESISTANCE IN STOMATE.F BEFORE FINAL ITERATION
C
4250  IF(ICHK.EQ.1)THEN
      GO TO 4500
      ELSE
      ICHK=1
      CALL STOMATE(I,J,NZ,NY,NX)
      ENDIF
      ENDIF
4000  CONTINUE
4500  CONTINUE
C
C     FINAL CANOPY TEMPERATURE, DIFFERENCE WITH AIR TEMPERATURE
C
C     TKC=final estimate of canopy temperature TKCZ
C     TKA=current air temperature
C     DTKC=TKC-TKA for next hour
C
      TKC(NZ,NY,NX)=TKCZ(NZ,NY,NX)
      TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-273.15
      DTKC(NZ,NY,NX)=TKC(NZ,NY,NX)-TKA(NY,NX)
C
C     IF CONVERGENCE NOT ACHIEVED (RARE), SET DEFAULT
C     TEMPERATURES, ENERGY FLUXES, WATER POTENTIALS, RESISTANCES
C
      IF(NN.GE.MXN)THEN
      WRITE(*,9999)IYRC,I,J,NX,NY,NZ
9999  FORMAT('CONVERGENCE FOR WATER UPTAKE NOT ACHIEVED ON   ',6I4)
      IF(DIFF.GT.0.5)THEN
      RAD1(NZ,NY,NX)=0.0
      EFLXC(NZ,NY,NX)=0.0
      SFLXC(NZ,NY,NX)=0.0
      HFLXC(NZ,NY,NX)=0.0
      EVAPC(NZ,NY,NX)=0.0
      EP(NZ,NY,NX)=0.0
      TKC(NZ,NY,NX)=TKA(NY,NX)+DTKC(NZ,NY,NX)
      TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-273.15
      FTHRM=EMMC*2.04E-10*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      THRM1(NZ,NY,NX)=FTHRM*TKC(NZ,NY,NX)**4
      PSILT(NZ,NY,NX)=PSIST1(NG(NZ,NY,NX))
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
      CCPOLT=CCPOLP(NZ,NY,NX)+CZPOLP(NZ,NY,NX)+CPPOLP(NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSILO(NZ,NY,NX)=FDMP/0.16*OSMO(NZ,NY,NX)
     2-8.3143*TKC(NZ,NY,NX)*FDMP*CCPOLT/OSWT
      PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
      WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      RC(NZ,NY,NX)=RSMN(NZ,NY,NX)+(RSMH(NZ,NY,NX)
     2-RSMN(NZ,NY,NX))*WFNC
      RA(NZ,NY,NX)=RAZ(NZ,NY,NX)
      VHCPC(NZ,NY,NX)=4.19*(WTSHT(NZ,NY,NX)*10.0E-06)
      DTKC(NZ,NY,NX)=0.0
      DO 4290 N=1,MY(NZ,NY,NX)
      DO 4290 L=NU(NY,NX),NI(NZ,NY,NX)
      PSIRT(N,L,NZ,NY,NX)=PSIST1(L)
      APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
      FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
      CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSIRO(N,L,NZ,NY,NX)=FDMR/0.16*OSMO(NZ,NY,NX)
     2-8.3143*TKS(L,NY,NX)*FDMR*CCPOLT/OSWT
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)
     2-PSIRO(N,L,NZ,NY,NX))
      UPWTR(N,L,NZ,NY,NX)=0.0
4290  CONTINUE
      ENDIF
      ENDIF
C
C     CANOPY SURFACE WATER STORAGE, SENSIBLE AND STORAGE HEAT FLUXES
C     (NOT EXPLICITLY CALCULATED IN CONVERGENCE SOLUTION)
C
C     VOLWP,VOLWC=water volume in canopy,on canopy surfaces
C     SFLXC,HFLXC=canopy sensible,storage heat fluxes
C     VHCPX,VHCPC=previous,current canopy heat capacity
C     PARHC=canopy sensible heat conductance
C     VFLXC=convective heat flux from latent heat flux
C     HFLWC1=convective heat flux from precip to canopy
C
      VOLWP(NZ,NY,NX)=VOLWP(NZ,NY,NX)+EP(NZ,NY,NX)-UPRT
      VOLWC(NZ,NY,NX)=VOLWC(NZ,NY,NX)+FLWC(NZ,NY,NX)+EVAPC(NZ,NY,NX)
      SFLXC(NZ,NY,NX)=PARHC*(TKA(NY,NX)-TKCZ(NZ,NY,NX))
      HFLXC(NZ,NY,NX)=TKCX*VHCPX-TKCZ(NZ,NY,NX)*VHCPC(NZ,NY,NX)
     2+VFLXC+HFLWC1
C
C     ROOT TOTAL, OSMOTIC AND TURGOR WATER POTENTIALS
C
C     PSIRT,PSILT=root,canopy total water potential
C     PSIST1=total soil water potential PSIST adjusted for surf elevn
C     RSSX,RSRS,RSRT=soil,soil+root,root radial+axial resistance
C     PSIRO,PSIRG=root osmotic,turgor water potential
C     FDMR=dry matter content
C     OSMO=osmotic potential at PSIRT=0 from PFT file
C
C
      DO 4505 N=1,MY(NZ,NY,NX)
      DO 4510 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(ILYR(N,L).EQ.1)THEN
      PSIRT(N,L,NZ,NY,NX)=AMIN1(0.0,(PSIST1(L)*RSRT(N,L)
     2+PSILT(NZ,NY,NX)*RSSX(N,L))/RSRS(N,L))
      APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
      FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
      CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSIRO(N,L,NZ,NY,NX)=FDMR/0.16*OSMO(NZ,NY,NX)
     2-8.3143*TKS(L,NY,NX)*FDMR*CCPOLT/OSWT
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)
     2-PSIRO(N,L,NZ,NY,NX))
      ELSE
      PSIRT(N,L,NZ,NY,NX)=PSIST1(L)
      APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
      FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
      CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSIRO(N,L,NZ,NY,NX)=FDMR/0.16*OSMO(NZ,NY,NX)
     2-8.3143*TKS(L,NY,NX)*FDMR*CCPOLT/OSWT
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)
     2-PSIRO(N,L,NZ,NY,NX))
      ENDIF
C     IF(I.EQ.284)THEN
C     WRITE(*,1256)'PSIRT',I,J,NX,NY,NZ,NN,PSIRT(N,L,NZ,NY,NX)
C    2,PSIST1(L),RSRT(N,L),PSILT(NZ,NY,NX),RSSX(N,L),RSRS(N,L)
C    3,RSRG(N,L),RSR1(N,L),RSR2(N,L),RTAR2,VOLWM(NPH,L,NY,NX)
1256  FORMAT(A8,6I4,20E12.4)
C     ENDIF
4510  CONTINUE
4505  CONTINUE
C
C     DEFAULT VALUES IF PLANT SPECIES DOES NOT EXIST
C
      ELSE
      RAD1(NZ,NY,NX)=0.0
      EFLXC(NZ,NY,NX)=0.0
      SFLXC(NZ,NY,NX)=0.0
      HFLXC(NZ,NY,NX)=0.0
      EVAPC(NZ,NY,NX)=0.0
      EP(NZ,NY,NX)=0.0
      IF(ZC(NZ,NY,NX).GE.DPTHS(NY,NX)-ZERO)THEN
      TKC(NZ,NY,NX)=TKA(NY,NX)
      ELSE
      TKC(NZ,NY,NX)=TKW(1,NY,NX)
      ENDIF
      TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-273.15
      FTHRM=EMMC*2.04E-10*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      THRM1(NZ,NY,NX)=FTHRM*TKC(NZ,NY,NX)**4
      PSILT(NZ,NY,NX)=PSIST1(NG(NZ,NY,NX))
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDMP=0.16+0.10*APSILT/(0.05*APSILT+2.0)
      CCPOLT=CCPOLP(NZ,NY,NX)+CZPOLP(NZ,NY,NX)+CPPOLP(NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSILO(NZ,NY,NX)=FDMP/0.16*OSMO(NZ,NY,NX)
     2-8.3143*TKC(NZ,NY,NX)*FDMP*CCPOLT/OSWT
      PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
      WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      RC(NZ,NY,NX)=RSMN(NZ,NY,NX)+(RSMH(NZ,NY,NX)
     2-RSMN(NZ,NY,NX))*WFNC
      RA(NZ,NY,NX)=RAZ(NZ,NY,NX)
      VHCPC(NZ,NY,NX)=4.19*(WTSHT(NZ,NY,NX)*10.0E-06)
      DTKC(NZ,NY,NX)=0.0
      DO 4300 N=1,MY(NZ,NY,NX)
      DO 4300 L=NU(NY,NX),NI(NZ,NY,NX)
      PSIRT(N,L,NZ,NY,NX)=PSIST1(L)
      APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
      FDMR=0.16+0.10*APSIRT/(0.05*APSIRT+2.0)
      CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSIRO(N,L,NZ,NY,NX)=FDMR/0.16*OSMO(NZ,NY,NX)
     2-8.3143*TKS(L,NY,NX)*FDMR*CCPOLT/OSWT
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)
     2-PSIRO(N,L,NZ,NY,NX))
      UPWTR(N,L,NZ,NY,NX)=0.0
4300  CONTINUE
      ENDIF
C
C     SET CANOPY GROWTH TEMPERATURE FROM SOIL SURFACE
C     OR CANOPY TEMPERATURE DEPENDING ON GROWTH STAGE
C
      IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0)THEN
      TKG(NZ,NY,NX)=TKS(NU(NY,NX),NY,NX)
C     ELSEIF((IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)
C    2.AND.IDAY(2,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0)THEN
C     TKG(NZ,NY,NX)=TKS(NU(NY,NX),NY,NX)
      ELSE
      TKG(NZ,NY,NX)=TKC(NZ,NY,NX)
      ENDIF
      TCG(NZ,NY,NX)=TKG(NZ,NY,NX)-273.15
C
C     ARRHENIUS FUNCTION FOR CANOPY AND ROOT GROWTH WITH OFFSET
C     FOR ZONE OF THERMAL ADAPTATION ENTERED IN 'READQ'
C
C     TKG,TKGO=canopy temperature,canopy temp used in Arrhenius eqn
C     TKS,TKSO=soil temperature,soil temp used in Arrhenius eqn
C     OFFST=shift in Arrhenius curve for thermal adaptation
C     TFN3,TFN4=temperature function for canopy,root growth (25 oC =1)
C     8.3143,710.0=gas constant,enthalpy
C     62500,197500,222500=energy of activn,high,low temp inactivn(KJ mol-1)
C     PSILZ=minimum daily canopy water potential
C
      TKGO=TKG(NZ,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKGO
      STK=710.0*TKGO
      ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFN3(NZ,NY,NX)=EXP(25.229-62500/RTK)/ACTV
      DO 100 L=NU(NY,NX),NI(NZ,NY,NX)
      TKSO=TKS(L,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKSO
      STK=710.0*TKSO
      ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFN4(L,NZ,NY,NX)=EXP(25.229-62500/RTK)/ACTV
100   CONTINUE
      PSILZ(NZ,NY,NX)=AMIN1(PSILZ(NZ,NY,NX),PSILT(NZ,NY,NX))
C
C     DIURNAL CHILLING
C
C     CTC=chilling temperature from PFT file
C     CHILL=accumulated chilling hours used to limit CO2 fixn in stomate.f
C
      IF(TCC(NZ,NY,NX).LT.CTC(NZ,NY,NX))THEN
      CHILL(NZ,NY,NX)=AMIN1(24.0,CHILL(NZ,NY,NX)+1.0)
      ELSE
      CHILL(NZ,NY,NX)=AMAX1(0.0,CHILL(NZ,NY,NX)-1.0)
      ENDIF
C
C     NH3 EXCHANGE BETWEEN CANOPY AND ATMOSPHERE FROM NH3
C     CONCENTRATION DIFFERENCES 'CNH3E' (ATMOSPHERE FROM 'READS') AND
C     'CNH3P' (CANOPY), AND FROM STOMATAL + BOUNDARY LAYER RESISTANCE
C
C     SNH3P,SNH3X=NH3 solubility at TCC, 25 oC
C     TCC=canopy temperature (oC)
C     FDMP,FNH3P=canopy dry matter content,NH3 concentration
C     ARLFB,ARLFP=branch,canopy leaf area
C     CNH3P,CNH3E=gaseous NH3 concentration in branch,atmosphere
C     CZPOLB,ZPOOLB=nonstructural N concentration,content in branch
C     RNH3B=NH3 flux between atmosphere and branch
C     RA,RC=canopy boundary layer,stomatal resistance
C     FRADP=fraction of radiation received by each PFT canopy
C
      SNH3P=SNH3X*EXP(0.513-0.0171*TCC(NZ,NY,NX))
      FNH3P=1.0E-04*FDMP
      DO 105 NB=1,NBR(NZ,NY,NX)
      IF(WTLSB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.ARLFB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     3.AND.ARLFP(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CNH3P=AMAX1(0.0,FNH3P*CZPOLB(NB,NZ,NY,NX)/SNH3P)
      ZPOOLB=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
      RNH3B(NB,NZ,NY,NX)=AMIN1(0.1*ZPOOLB
     2,AMAX1((CNH3E(NY,NX)-CNH3P)/(RA(NZ,NY,NX)+RC(NZ,NY,NX))
     3*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
     3*ARLFB(NB,NZ,NY,NX)/ARLFP(NZ,NY,NX),-0.1*ZPOOLB))
      ELSE
      RNH3B(NB,NZ,NY,NX)=0.0
      ENDIF
C     WRITE(*,7777)'RNH3',I,J,NZ,NB,RNH3B(NB,NZ,NY,NX)
C    2,RNH3C(NZ,NY,NX),CNH3E(NY,NX),CNH3P,RA(NZ,NY,NX),RC(NZ,NY,NX)
C    2,ARLFB(NB,NZ,NY,NX),ARLFP(NZ,NY,NX),SNH3P
C    4,ZPOOL(NB,NZ,NY,NX),WTLSB(NB,NZ,NY,NX),FRADP(NZ,NY,NX)
7777  FORMAT(A8,4I4,40E24.16)
105   CONTINUE
C
C     ROOT(N=1) AD MYCORRHIZAL(N=2) O2 AND NUTRIENT UPTAKE
C
      DO 955 N=1,MY(NZ,NY,NX)
      DO 950 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX)
     2.AND.RTDNP(N,L,NZ,NY,NX).GT.ZERO
     3.AND.RTVLW(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     4.AND.THETW(L,NY,NX).GT.ZERO)THEN
      TFOXYX=0.0
      TFNH4X=0.0
      TFNHBX=0.0
      TFNO3X=0.0
      TFNOBX=0.0
      TFPO4X=0.0
      TFPOBX=0.0
      TFP14X=0.0
      TFP1BX=0.0
C
C     UPTAKE CAPACITY 'FWSRT' DEPENDS ON ROOT,MYCORRHIZAL
C     PROTEIN CONTENT RELATIVE TO 5% FOR WHICH ACTIVE UPTAKE
C     PARAMETERS ARE DEFINED
C
C     CWSRTL,CWSRT=current,maximum protein concentration
C     WSRTL,WTRTL=protein content,mass
C     FWSRT=protein concentration relative to 5%
C
      IF(WTRTL(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CWSRTL(N,L,NZ,NY,NX)=AMIN1(CWSRT(NZ,NY,NX)
     2,WSRTL(N,L,NZ,NY,NX)/WTRTL(N,L,NZ,NY,NX))
      FWSRT=CWSRTL(N,L,NZ,NY,NX)/0.05
      ELSE
      CWSRTL(N,L,NZ,NY,NX)=CWSRT(NZ,NY,NX)
      FWSRT=1.0
      ENDIF
C
C     RESPIRATION CONSTRAINT ON UPTAKE FROM NON-STRUCTURAL C
C
C     RCO2N=total respiration from CPOOLR
C     FCUP=limitation to active uptake respiration from CPOOLR
C     CPOOLR=nonstructural C content
C
      IF(RCO2N(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FCUP=AMAX1(0.0,AMIN1(1.0,0.25*CPOOLR(N,L,NZ,NY,NX)
     2/RCO2N(N,L,NZ,NY,NX)))
      ELSE
      FCUP=0.0
      ENDIF
C
C     FEEDBACK CONSTRAINT ON N UPTAKE FROM NON-STRUCTURAL N AND P
C
C     FZUP,FPUP=limitn to active uptake respiration from CZPOLR,CPPOLR
C     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
C     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition effect on N,P uptake
C     UPWTRH=water uptake at time step for gas flux calculations
C
      IF(CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      FZUP=AMIN1(CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX)
     2+CZPOLR(N,L,NZ,NY,NX)/ZCKI)
     3,CPPOLR(N,L,NZ,NY,NX)/(CPPOLR(N,L,NZ,NY,NX)
     4+CZPOLR(N,L,NZ,NY,NX)/ZPKI))
      FPUP=AMIN1(CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)/PCKI)
     3,CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX)
     4+CPPOLR(N,L,NZ,NY,NX)/PZKI))
      ELSE
      FZUP=0.0
      FPUP=0.0
      ENDIF
C     NN=0
      UPWTRP=AMAX1(0.0,-UPWTR(N,L,NZ,NY,NX)/PP(NZ,NY,NX))
      UPWTRH=UPWTRP*XNPG
C
C     FACTORS CONSTRAINING O2 AND NUTRIENT UPTAKE AMONG
C     COMPETING ROOT,MYCORRHIZAL AND MICROBIAL POPULATIONS
C     IN BAND AND NON-BAND SOIL ZONES FROM DEMAND CALCULATED
C     IN PREVIOUS HOUR
C
C     ROXYY=O2 demand by all microbial,root,myco populations
C     ROXYP=O2 demand by each root,myco population
C     FOXYX=fraction of ROXYY by each root,myco population
C     RNH4Y=NH4 demand in non-band by all microbial,root,myco populations
C     RUNNHP=NH4 demand in non-band by each root,myco population
C     FNH4X=fraction of RNH4Y by each root,myco populn
C     RNHBY=NH4 demand in band by all microbial,root,myco populations
C     RUNNBP=NH4 demand in band by each root,myco population
C     FNHBX=fraction of RNHBY by each root,myco populn
C     RNO3Y=NO3 demand in non-band by all microbial,root,myco populations
C     RUNNOP=NO3 demand in non-band by each root,myco population
C     FNO3X=fraction of RNO3Y by each root,myco populn
C     RN3BY=NO3 demand in band by all microbial,root,myco populations
C     RUNNXB=NO3 demand in band by each root,myco population
C     FNOBX=fraction of RN3BY by each root,myco populn
C     RPO4Y=H2PO4 demand in non-band by all microbial,root,myco populations
C     RUPP2P=H2PO4 demand in non-band by each root,myco population
C     FPO4X=fraction of RPO4Y by each root,myco populn
C     RPOBY=H2PO4 demand in band by all microbial,root,myco populations
C     RUPP2B=H2PO4 demand in band by each root,myco population
C     FPOBX=fraction of RPOBY by each root,myco populn
C     RP14Y=HPO4 demand in non-band by all microbial,root,myco populations
C     RUPP1P=HPO4 demand in non-band by each root,myco population
C     FP14X=fraction of RP14Y by each root,myco populn
C     RP1BY=HPO4 demand in band by all microbial,root,myco populations
C     RUPP1B=HPO4 demand in band by each root,myco population
C     FP1BX=fraction of RP1BY by each root,myco populn
C     FPP=minimum uptake fraction
C     FPQ=PFT fraction of biome root mass
C
      IF(ROXYY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOXYX=AMAX1(FPP(N,L,NZ),ROXYP(N,L,NZ,NY,NX)/ROXYY(L,NY,NX))
      ELSE
      FOXYX=FPQ(N,L,NZ)
      ENDIF
      IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4X=AMAX1(FPP(N,L,NZ),RUNNHP(N,L,NZ,NY,NX)/RNH4Y(L,NY,NX))
      ELSE
      FNH4X=FPQ(N,L,NZ)
      ENDIF
      IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNHBX=AMAX1(FPP(N,L,NZ),RUNNBP(N,L,NZ,NY,NX)/RNHBY(L,NY,NX))
      ELSE
      FNHBX=FPQ(N,L,NZ)
      ENDIF
      IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3X=AMAX1(FPP(N,L,NZ),RUNNOP(N,L,NZ,NY,NX)/RNO3Y(L,NY,NX))
      ELSE
      FNO3X=FPQ(N,L,NZ)
      ENDIF
      IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNOBX=AMAX1(FPP(N,L,NZ),RUNNXP(N,L,NZ,NY,NX)/RN3BY(L,NY,NX))
      ELSE
      FNOBX=FPQ(N,L,NZ)
      ENDIF
      IF(RPO4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FPO4X=AMAX1(FPP(N,L,NZ),RUPP2P(N,L,NZ,NY,NX)/RPO4Y(L,NY,NX))
      ELSE
      FPO4X=FPQ(N,L,NZ)
      ENDIF
      IF(RPOBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FPOBX=AMAX1(FPP(N,L,NZ),RUPP2B(N,L,NZ,NY,NX)/RPOBY(L,NY,NX))
      ELSE
      FPOBX=FPQ(N,L,NZ)
      ENDIF
      IF(RP14Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FP14X=AMAX1(FPP(N,L,NZ),RUPP1P(N,L,NZ,NY,NX)/RP14Y(L,NY,NX))
      ELSE
      FP14X=FPQ(N,L,NZ)
      ENDIF
      IF(RP1BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FP1BX=AMAX1(FPP(N,L,NZ),RUPP1B(N,L,NZ,NY,NX)/RP1BY(L,NY,NX))
      ELSE
      FP1BX=FPQ(N,L,NZ)
      ENDIF
      TFOXYX=TFOXYX+FOXYX
      TFNH4X=TFNH4X+FNH4X
      TFNO3X=TFNO3X+FNO3X
      TFPO4X=TFPO4X+FPO4X
      TFP14X=TFP14X+FP14X
      TFNHBX=TFNHBX+FNHBX
      TFNOBX=TFNOBX+FNOBX
      TFPOBX=TFPOBX+FPOBX
      TFP1BX=TFP1BX+FP1BX
C
C     ROOT O2 DEMAND CALCULATED FROM O2 NON-LIMITED RESPIRATION RATE
C
C     ROXYP=O2 demand
C     RCO2M=respiration unlimited by O2
C     RTVLW=root or myco aqueous volume
C     FOXYX=fraction of total O2 demand from previous hour
C
      ROXYP(N,L,NZ,NY,NX)=2.667*RCO2M(N,L,NZ,NY,NX)
      IF(RCO2M(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.RTVLW(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.FOXYX.GT.ZEROQ(NZ,NY,NX))THEN
C
C     INITIALIZE VARIABLES USED IN ROOT GAS EXCHANGE
C     (CO2, O2, CH4, N2, N2O, NH3, H2)
C
C     CO2A1,CO2P1,CO2G1,CO2S1=gaseous,aqueous CO2 in root,soil
C     OXYA1,OXYP1,OXYG1,OXYS1=gaseous,aqueous O2 in root,soil
C     CH4A1,CH4P1,CH4G1,CH4S1=gaseous,aqueous CH4 in root,soil
C     Z2OA1,Z2OP1,Z2OG1,Z2OS1=gaseous,aqueous N2O in root,soil
C     ZH3A1,ZH3P1,ZH3G1,ZH3S1=gaseous,aqueous NH3 in root,soil
C     H2GA1,H2GP1,H2GG1,H2GS1=gaseous,aqueous H2 in root,soil
C     CCH4S1,CCH4P1=aqueous CH4 concentration in soil,root
C     CN2OS1,CN2OP1=aqueous N2O concentration in soil,root
C     CNH3S1,CNH3B1,CNH3P1=aqueous NH3 concn in soil non-band,band,root
C     CH2GS1,CH2GP1=aqueous H2 concentration in soil,root
C     RTVLWA,RTVLWB=root aqueous volume in non-band,band
C     XNPG=time step of flux calculation
C     UPMXP=O2 demand per plant at time step of flux calculation
C     ROXYFX=net O2 gas flux at time step of flux calculation
C     RCO2FX=net CO2 gas flux at time step of flux calculation
C     ROXYLX=net O2 aqueous flux at time step of flux calculation
C
      CO2A1=AMAX1(ZEROP(NZ,NY,NX),CO2A(N,L,NZ,NY,NX))
      CO2P1=AMAX1(ZEROP(NZ,NY,NX),CO2P(N,L,NZ,NY,NX))
      CO2G1=AMAX1(ZEROP(NZ,NY,NX),CO2G(L,NY,NX)*FPQ(N,L,NZ))
      CO2S1=AMAX1(ZEROP(NZ,NY,NX),CO2S(L,NY,NX)*FPQ(N,L,NZ))
      OXYA1=AMAX1(ZEROP(NZ,NY,NX),OXYA(N,L,NZ,NY,NX))
      OXYP1=AMAX1(ZEROP(NZ,NY,NX),OXYP(N,L,NZ,NY,NX))
      OXYG1=AMAX1(ZEROP(NZ,NY,NX),OXYG(L,NY,NX)*FOXYX)
      OXYS1=OXYS(L,NY,NX)*FOXYX
      CH4A1=CH4A(N,L,NZ,NY,NX)
      CH4P1=CH4P(N,L,NZ,NY,NX)
      CH4S1=CH4S(L,NY,NX)*FPQ(N,L,NZ)
      CCH4S1=CCH4S(L,NY,NX)
      CCH4P1=AMAX1(0.0,CH4P1/RTVLW(N,L,NZ,NY,NX))
      Z2OA1=Z2OA(N,L,NZ,NY,NX)
      Z2OP1=Z2OP(N,L,NZ,NY,NX)
      Z2OS1=Z2OS(L,NY,NX)*FPQ(N,L,NZ)
      CN2OS1=CZ2OS(L,NY,NX)
      CN2OP1=AMAX1(0.0,Z2OP1/RTVLW(N,L,NZ,NY,NX))
      ZH3A1=ZH3A(N,L,NZ,NY,NX)
      ZH3P1=ZH3P(N,L,NZ,NY,NX)
      ZH3S1=ZNH3S(L,NY,NX)*FPQ(N,L,NZ)
      ZH3B1=ZNH3B(L,NY,NX)*FPQ(N,L,NZ)
      CNH3S1=CNH3S(L,NY,NX)
      CNH3B1=CNH3B(L,NY,NX)
      CNH3P1=AMAX1(0.0,ZH3P1/RTVLW(N,L,NZ,NY,NX))
      H2GA1=H2GA(N,L,NZ,NY,NX)
      H2GP1=H2GP(N,L,NZ,NY,NX)
      H2GS1=H2GS(L,NY,NX)*FPQ(N,L,NZ)
      CH2GS1=CH2GS(L,NY,NX)
      CH2GP1=AMAX1(0.0,H2GP1/RTVLW(N,L,NZ,NY,NX))
      RTVLWA=RTVLW(N,L,NZ,NY,NX)*VLNH4(L,NY,NX)
      RTVLWB=RTVLW(N,L,NZ,NY,NX)*VLNHB(L,NY,NX)
      UPMXP=ROXYP(N,L,NZ,NY,NX)*XNPG/PP(NZ,NY,NX)
      ROXYFX=ROXYF(L,NY,NX)*FOXYX*XNPG
      RCO2FX=RCO2F(L,NY,NX)*FOXYX*XNPG
      ROXYLX=ROXYL(L,NY,NX)*FOXYX*XNPG
C
C     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
C
C     *SGL1=diffusivity
C     PORTX=tortuosity effect of root porosity on diffusivity
C     CG=CO2g,OG=O2g,CH=CH4g,Z2=N2Og,ZH=NH3g,HG=H2g
C     CL=CO2s,OL=O2s,CQ=CH4s,ZV=N2Os,ZN=NH3s,HL=H2s
C
      CGSGL1=CGSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
      OGSGL1=OGSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
      CHSGL1=CHSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
      Z2SGL1=Z2SGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
      ZHSGL1=ZHSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
      HGSGL1=HGSGL(L,NY,NX)*XNPG*PORTX(N,NZ,NY,NX)
      CLSGL1=CLSGL(L,NY,NX)*XNPG*FOXYX
      OLSGL1=OLSGL(L,NY,NX)*XNPG*FOXYX
      CQSGL1=CQSGL(L,NY,NX)*XNPG*FOXYX
      ZVSGL1=ZVSGL(L,NY,NX)*XNPG*FOXYX
      ZNSGL1=ZNSGL(L,NY,NX)*XNPG*FOXYX
      HLSGL1=HLSGL(L,NY,NX)*XNPG*FOXYX
      OLSGLP=OLSGL(L,NY,NX)*XNPG
      ROXDFQ=0.0
      RCHDFQ=0.0
      RN2DFQ=0.0
      RNHDFQ=0.0
      RHGDFQ=0.0
      ROXDF1=0.0
      RCHDF1=0.0
      RN2DF1=0.0
      RNHDF1=0.0
      RHGDF1=0.0
C
C     ROOT CONDUCTANCE TO GAS TRANSFER
C
C     WTRTS=total root,myco mass
C     FRTDPX=fraction of each soil layer with primary root
C     RTCR1,RTCR2,RTCRA=cross-sectional area/length of
C     primary,secondary,total root,myco system
C     RTN1,RTNL=number of root,myco primary,secondary axes
C     RRAD1,RRAD2=primary,secondary root radius
C     DPTHZ=depth of primary root from surface
C     RTLGA=average secondary root length
C
      IF(WTRTS(NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.FRTDPX(L,NZ).GT.ZERO)THEN
      RTCR1=AMAX1(PP(NZ,NY,NX),RTN1(N,L,NZ,NY,NX))
     2*3.1416*RRAD1(N,L,NZ,NY,NX)**2
     2/DPTHZ(L,NY,NX)
      RTCR2=(RTNL(N,L,NZ,NY,NX)*3.1416*RRAD2(N,L,NZ,NY,NX)**2
     2/RTLGA(N,L,NZ,NY,NX))/FRTDPX(L,NZ)
      IF(RTCR2.GT.RTCR1)THEN
      RTCRA=RTCR1*RTCR2/(RTCR1+RTCR2)
      ELSE
      RTCRA=RTCR1
      ENDIF
      ELSE
      RTCRA=0.0
      ENDIF
C
C     VARIABLES USED TO CALCULATE ROOT GAS TRANSFER
C     BETWEEN AQUEOUS AND GASEOUS PHASES
C
C     RTLGP=root,myco length per plant
C     IDAY(1,=emergence date
C     RTARR=root surface area/radius for uptake
C     RRADP=path length for radial diffusion within root
C     DIFOP=aqueous diffusivity of O2 within root
C     RTVLW=root,myco aqueous volume
C     S*L=solubility of gas in water from hour1.f:
C     CO2=CO2,OXY=O2,CH4=CH4,N2O=N2O,NH3=NH3,H2G=H2
C     DF*A=root-atmosphere gas conductance
C     DFGP=rate const for equilibrn of gas concn in gaseous-aqueous phases
C     RCO2PX=root CO2 gas flux at time step for gas flux calculations
C     RCO2A=root CO2 flux from grosub.f
C
      IF(N.EQ.1.AND.IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).GT.0
     2.AND.RTLGP(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RTARRX=RTARR(N,L)/RRADP(N,NZ,NY,NX)
      DIFOP=OLSGLP*RTARRX
      VOLWCA=RTVLW(N,L,NZ,NY,NX)*SCO2L(L,NY,NX)
      VOLWOA=RTVLW(N,L,NZ,NY,NX)*SOXYL(L,NY,NX)
      VOLWC4=RTVLW(N,L,NZ,NY,NX)*SCH4L(L,NY,NX)
      VOLWZA=RTVLW(N,L,NZ,NY,NX)*SN2OL(L,NY,NX)
      VOLWNA=RTVLW(N,L,NZ,NY,NX)*SNH3L(L,NY,NX)
      VOLWH2=RTVLW(N,L,NZ,NY,NX)*SH2GL(L,NY,NX)
      DFCOA=CGSGL1*RTCRA
      DFOXA=OGSGL1*RTCRA
      DFCHA=CHSGL1*RTCRA
      DFN2A=Z2SGL1*RTCRA
      DFNHA=ZHSGL1*RTCRA
      DFHGA=HGSGL1*RTCRA
      ELSE
      RTARRX=0.0
      DIFOP=0.0
      VOLWCA=0.0
      VOLWOA=0.0
      VOLWC4=0.0
      VOLWZA=0.0
      VOLWNA=0.0
      VOLWH2=0.0
      DFCOA=0.0
      DFOXA=0.0
      DFCHA=0.0
      DFN2A=0.0
      DFNHA=0.0
      DFHGA=0.0
      ENDIF
      DFGP=AMIN1(1.0,XNPD*SQRT(PORT(N,NZ,NY,NX))*TFND(L,NY,NX))
      RCO2PX=-RCO2A(N,L,NZ,NY,NX)*XNPG
C
C     SOLVE FOR GAS EXCHANGE IN SOIL AND ROOTS DURING ROOT UPTAKE
C     AT SMALLER TIME STEP NPH
C
      DO 99 M=1,NPH
C
C     AQUEOUS GAS DIFFUSIVITY THROUGH SOIL WATER TO ROOT
C
C     gas code:CO2=CO2,OXY=O2,CH4=CH4,Z2O=N2O,NH3=NH3 non-band,
C     NHB=NH3 band,H2G=H2
C     VOLWMM,VOLPMM=soil micropore water,air volume
C     FOXYX=root fraction of total O2 demand from previous hour
C     FPQ=PFT fraction of biome root mass
C     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
C     VOLX=soil volume excluding rock,macropores
C     THETW1=soil water concentration
C     TORT=soil tortuosity
C     FILM=soil water film thickness
C     RRADL=root radius
C     RRADS=path length for radial diffusion from soil to root
C     DIF*=aqueous diffusivity from soil to root:OL=O2,CL=CH4
C     ZL=N2O,NL=NH3 non-band,NB=NH4 band,HL=H2
C     C*G=soil gaseous concentration
C     VOLW*,VOLP*=VOLWMM,VOLPMM*gas solubility
C
      VOLWMO=VOLWM(M,L,NY,NX)*FOXYX
      VOLWMM=VOLWM(M,L,NY,NX)*FPQ(N,L,NZ)
      VOLPMM=VOLPM(M,L,NY,NX)*FPQ(N,L,NZ)
      VOLWSP=RTVLW(N,L,NZ,NY,NX)+VOLWMM
      VOLWMA=VOLWMM*VLNH4(L,NY,NX)
      VOLWMB=VOLWMM*VLNHB(L,NY,NX)
      VOLWSA=RTVLWA+VOLWMA
      VOLWSB=RTVLWB+VOLWMB
      THETW1=AMAX1(0.0,VOLWM(M,L,NY,NX)/VOLY(L,NY,NX))
      IF(THETW1.GT.THETY(L,NY,NX)
     2.AND.FPQ(N,L,NZ).GT.ZEROQ(NZ,NY,NX))THEN
      THETM=TORT(M,L,NY,NX)*THETW1
      RRADS=LOG((FILM(M,L,NY,NX)+RRADL(N,L))/RRADL(N,L))
      RTARRX=RTARR(N,L)/RRADS
      DIFOL=THETM*OLSGL1*RTARRX
      DIFCL=THETM*CQSGL1*RTARRX
      DIFZL=THETM*ZVSGL1*RTARRX
      DIFNL=THETM*ZNSGL1*RTARRX*VLNH4(L,NY,NX)
      DIFNB=THETM*ZNSGL1*RTARRX*VLNHB(L,NY,NX)
      DIFHL=THETM*HLSGL1*RTARRX
      CH4G1=CCH4G(L,NY,NX)*VOLPMM
      Z2OG1=CZ2OG(L,NY,NX)*VOLPMM
      ZH3G1=CNH3G(L,NY,NX)*VOLPMM
      H2GG1=CH2GG(L,NY,NX)*VOLPMM
      VOLWCO=VOLWMM*SCO2L(L,NY,NX)
      VOLWOX=VOLWMM*SOXYL(L,NY,NX)
      VOLWCH=VOLWMM*SCH4L(L,NY,NX)
      VOLWN2=VOLWMM*SN2OL(L,NY,NX)
      VOLWNH=VOLWMM*SNH3L(L,NY,NX)*VLNH4(L,NY,NX)
      VOLWNB=VOLWMM*SNH3L(L,NY,NX)*VLNHB(L,NY,NX)
      VOLWHG=VOLWMM*SH2GL(L,NY,NX)
      VOLPNH=VOLPMM*VLNH4(L,NY,NX)
      VOLPNB=VOLPMM*VLNHB(L,NY,NX)
C
C     MASS FLOW OF GAS FROM SOIL TO ROOT AT SHORTER TIME STEP NPT
C
C     C*S1=soil aqueous concentration non-band
C     C*B1=soil aqueous concentration band
C     C*A1=root gaseous concentration
C     C*P1=root aqueous concentration
C     ROXYLX=soil net O2 aqueous flux
C     VOLWMM=micropore water volume
C     RTVLW,RTVLP=root aqueous,gaseous volume
C     RMF*=soil convective solute flux:COS=CO2,OXS=O2,CHS=CH4,
C     N2S=N2O,NHS=NH3 non-band,NHB=NH3 band,HGS=H2
C     UPWTRH=water uptake
C
      DO 90 MX=1,NPT
      OXYS1=OXYS1+ROXYLX
      CCO2S1=AMAX1(0.0,CO2S1/VOLWMM)
      COXYS1=AMIN1(COXYE(NY,NX)*SOXYL(L,NY,NX)
     2,AMAX1(0.0,OXYS1/VOLWMO))
      CCH4S1=AMAX1(0.0,CH4S1/VOLWMM)
      CN2OS1=AMAX1(0.0,Z2OS1/VOLWMM)
      CNH3S1=AMAX1(0.0,ZH3S1/VOLWMM)
      CNH3B1=AMAX1(0.0,ZH3B1/VOLWMM)
      CH2GS1=AMAX1(0.0,H2GS1/VOLWMM)
      IF(RTVLP(N,L,NZ,NY,NX).GT.ZERO)THEN
      CCO2A1=AMAX1(0.0,CO2A1/RTVLP(N,L,NZ,NY,NX))
      COXYA1=AMAX1(0.0,OXYA1/RTVLP(N,L,NZ,NY,NX))
      CCH4A1=AMAX1(0.0,CH4A1/RTVLP(N,L,NZ,NY,NX))
      CZ2OA1=AMAX1(0.0,Z2OA1/RTVLP(N,L,NZ,NY,NX))
      CNH3A1=AMAX1(0.0,ZH3A1/RTVLP(N,L,NZ,NY,NX))
      CH2GA1=AMAX1(0.0,H2GA1/RTVLP(N,L,NZ,NY,NX))
      ELSE
      CCO2A1=0.0
      COXYA1=0.0
      CCH4A1=0.0
      CZ2OA1=0.0
      CNH3A1=0.0
      CH2GA1=0.0
      ENDIF
      CCO2P1=AMAX1(0.0,CO2P1/RTVLW(N,L,NZ,NY,NX))
      COXYP1=AMIN1(COXYE(NY,NX)*SOXYL(L,NY,NX)
     2,AMAX1(0.0,OXYP1/RTVLW(N,L,NZ,NY,NX)))
      CCH4P1=AMAX1(0.0,CH4P1/RTVLW(N,L,NZ,NY,NX))
      CN2OP1=AMAX1(0.0,Z2OP1/RTVLW(N,L,NZ,NY,NX))
      CNH3P1=AMAX1(0.0,ZH3P1/RTVLW(N,L,NZ,NY,NX))
      CH2GP1=AMAX1(0.0,H2GP1/RTVLW(N,L,NZ,NY,NX))
      DIFOX=DIFOL+DIFOP
      RMFCOS=UPWTRH*CCO2S1
      RMFOXS=UPWTRH*COXYS1
      RMFCHS=UPWTRH*CCH4S1
      RMFN2S=UPWTRH*CN2OS1
      RMFN3S=UPWTRH*CNH3S1*VLNH4(L,NY,NX)
      RMFN3B=UPWTRH*CNH3B1*VLNHB(L,NY,NX)
      RMFHGS=UPWTRH*CH2GS1
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF O2 IN AQUEOUS PHASES OF
C     SOIL AND ROOT = ACTIVE UPTAKE OF O2 BY ROOT
C
C     DIFOL=O2 aqueous diffusivity from soil to root
C     UPWTRH=water uptake
C     DIFOP=aqueous diffusivity of O2 within root
C     COXYS1,COXYP1=soil,root aqueous O2 concentration
C     UPMXP=O2 demand per plant
C     RUPOXR=root O2 uptake per plant
C     COXYR=aqueous O2 concentration at root surface
C     RDFOXS,RDFOXP=aqueous O2 diffusion per plant:soil-root,within root
C
      X=(DIFOL+UPWTRH)*COXYS1+DIFOP*COXYP1
      IF(X.GT.ZERO.AND.OXYS1.GT.ZEROP(NZ,NY,NX))THEN
      B=-UPMXP-DIFOX*OXKM-X
      C=X*UPMXP
      delta=B*B-4.0*C
      if(abs(delta)<ZERO)delta=0.0
      RUPOXR=(-B-SQRT(delta))/2.0
      COXYR=(X-RUPOXR)/DIFOX
      RDFOXS=RMFOXS+DIFOL*(COXYS1-COXYR)
      RDFOXP=DIFOP*(COXYP1-COXYR)
      ELSE
      X=DIFOP*COXYP1
      IF(X.GT.ZERO.AND.OXYP1.GT.ZEROP(NZ,NY,NX))THEN
      B=-UPMXP-DIFOP*OXKM-X
      C=X*UPMXP
      DELTA=B*B-4.0*C
      IF(ABS(DELTA)<ZERO)DELTA=ZERO
      RUPOXR=(-B-SQRT(DELTA))/2.0
      COXYR=(X-RUPOXR)/DIFOP
      RDFOXS=0.0
      RDFOXP=DIFOP*(COXYP1-COXYR)
      ELSE
      RUPOXR=0.0
      COXYR=0.0
      RDFOXS=0.0
      RDFOXP=0.0
      ENDIF
      ENDIF
C
C     MASS FLOW + DIFFUSIVE EXCHANGE OF OTHER GASES
C     BETWEEN ROOT AND SOIL, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RUPOSX,RUPOPX=aqueous O2 uptake from soil,root
C     PP=PFT population
C     RDFCOS,RCO2SX=aqueous CO2 soil-root diffusion,root uptake
C     RDFCHS,RUPCSX=aqueous CH4 soil-root diffusion,root uptake
C     RDFN2S,RUPZSX=aqueous N2O soil-root diffusion,root uptake
C     RDFN3S,RUPNSX=aqueous NH3 soil-root diffusion,root uptake:non-band
C     RDFN3B,RUPNBX=aqueous NH3 soil-root diffusion,root uptake:band
C     RDFHGS,RUPHGX=aqueous H2 soil-root diffusion,root uptake
C     RMF*=soil convective solute flux
C     DIF*=aqueous diffusivity from soil to root
C     C*S1=soil aqueous concentration non-band
C     C*B1=soil aqueous concentration band
C     C*P1=root aqueous concentration
C
C     IF(IYRC.EQ.2000.AND.I.LE.180)THEN
C     WRITE(*,5555)'COXYR',I,J,NX,NY,NZ,L,N,M,MX,COXYR,RUPOXR
C    2,RMFOXS,RDFOXS,RDFOXP,COXYS1,COXYS1-COXYR,COXYP1,FOXYX
C    3,WTRTG(L),DIFOL,DIFOP,THETM,OLSGL1,UPWTRH,RTARR(N,L)
C    5,RTARRX,UPMXP,THETW(L,NY,NX),OXYS1,OXYS(L,NY,NX),OXYP1
C    3,OXYP(N,L,NZ,NY,NX),ROXYY(L,NY,NX),RTLGP(N,L,NZ,NY,NX)
C    2,UPMXP,DIFOX,THETW1,THETM,RRADS,FPQ(N,L,NZ)
C    4,RUPOXS(N,L,NZ,NY,NX),RUPOXP(N,L,NZ,NY,NX)
C    5,COXYE(NY,NX),SOXYL(L,NY,NX),FRTDPX(L,NZ)
5555  FORMAT(A8,9I4,40E12.4)
C     ENDIF
      RUPOSX=RDFOXS*PP(NZ,NY,NX)
      RUPOPX=RDFOXP*PP(NZ,NY,NX)
      RDFCOS=RMFCOS+DIFCL*(CCO2S1-CCO2P1)
      RDXCOS=(RTVLW(N,L,NZ,NY,NX)*AMAX1(ZEROP(NZ,NY,NX),CO2S1)
     2-VOLWMM*AMAX1(ZEROP(NZ,NY,NX),CO2P1))/VOLWSP
      IF(RDFCOS.GT.0.0)THEN
      RCO2SX=AMIN1(AMAX1(0.0,RDXCOS),RDFCOS*PP(NZ,NY,NX))
      ELSE
      RCO2SX=AMAX1(AMIN1(0.0,RDXCOS),RDFCOS*PP(NZ,NY,NX))
      ENDIF
      IF(N.EQ.1)THEN
      RDFCHS=RMFCHS+DIFCL*(CCH4S1-CCH4P1)
      RDXCHS=(RTVLW(N,L,NZ,NY,NX)*AMAX1(ZEROP(NZ,NY,NX),CH4S1)
     2-VOLWMM*AMAX1(ZEROP(NZ,NY,NX),CH4P1))/VOLWSP
      IF(RDFCHS.GT.0.0)THEN
      RUPCSX=AMIN1(AMAX1(0.0,RDXCHS),RDFCHS*PP(NZ,NY,NX))
      ELSE
      RUPCSX=AMAX1(AMIN1(0.0,RDXCHS),RDFCHS*PP(NZ,NY,NX))
      ENDIF
      RDFN2S=RMFN2S+DIFZL*(CN2OS1-CN2OP1)
      RDXN2S=(RTVLW(N,L,NZ,NY,NX)*AMAX1(ZEROP(NZ,NY,NX),Z2OS1)
     2-VOLWMM*AMAX1(ZEROP(NZ,NY,NX),Z2OP1))/VOLWSP
      IF(RDFN2S.GT.0.0)THEN
      RUPZSX=AMIN1(AMAX1(0.0,RDXN2S),RDFN2S*PP(NZ,NY,NX))
      ELSE
      RUPZSX=AMAX1(AMIN1(0.0,RDXN2S),RDFN2S*PP(NZ,NY,NX))
      ENDIF
      RDFN3S=RMFN3S+DIFNL*(CNH3S1-CNH3P1)
      IF(VOLWSA.GT.ZEROP(NZ,NY,NX))THEN
      ZH3PA=ZH3P1*VLNH4(L,NY,NX)
      RDXNHS=(RTVLWA*AMAX1(ZEROP(NZ,NY,NX),ZH3S1)
     2-VOLWMA*AMAX1(ZEROP(NZ,NY,NX),ZH3PA))/VOLWSA
      ELSE
      RDXNHS=0.0
      ENDIF
      IF(RDFN3S.GT.0.0)THEN
      RUPNSX=AMIN1(AMAX1(0.0,RDXNHS),RDFN3S*PP(NZ,NY,NX))
      ELSE
      RUPNSX=AMAX1(AMIN1(0.0,RDXNHS),RDFN3S*PP(NZ,NY,NX))
      ENDIF
      RDFN3B=RMFN3B+DIFNB*(CNH3B1-CNH3P1)
      IF(VOLWSB.GT.ZEROP(NZ,NY,NX))THEN
      ZH3PB=ZH3P1*VLNHB(L,NY,NX)
      RDXNHB=(RTVLWB*AMAX1(ZEROP(NZ,NY,NX),ZH3B1)
     2-VOLWMB*AMAX1(ZEROP(NZ,NY,NX),ZH3PB))/VOLWSB
      ELSE
      RDXNHB=0.0
      ENDIF
      IF(RDFN3B.GT.0.0)THEN
      RUPNBX=AMIN1(AMAX1(0.0,RDXNHB),RDFN3B*PP(NZ,NY,NX))
      ELSE
      RUPNBX=AMAX1(AMIN1(0.0,RDXNHB),RDFN3B*PP(NZ,NY,NX))
      ENDIF
      RDFHGS=RMFHGS+DIFHL*(CH2GS1-CH2GP1)
      RDXHGS=(RTVLW(N,L,NZ,NY,NX)*AMAX1(ZEROP(NZ,NY,NX),H2GS1)
     2-VOLWMM*AMAX1(ZEROP(NZ,NY,NX),H2GP1))/VOLWSP
      IF(RDFHGS.GT.0.0)THEN
      RUPHGX=AMIN1(AMAX1(0.0,RDXHGS),RDFHGS*PP(NZ,NY,NX))
      ELSE
      RUPHGX=AMAX1(AMIN1(0.0,RDXHGS),RDFHGS*PP(NZ,NY,NX))
      ENDIF
      ELSE
      RUPCSX=0.0
      RUPZSX=0.0
      RUPNSX=0.0
      RUPNBX=0.0
      RUPHGX=0.0
      ENDIF
C
C     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN SOIL
C     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
C     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENTS
C     FROM 'WATSUB'
C
C     THETPM,THETX=air-filled porosity,minimum THETPM
C     R*DFQ=soil gas exchange between gaseous-aqueous phases
C     DFGS=rate constant for soil gas exchange from watsub.f
C     CO2G1,CO2S1=gaseous,aqueous CO2 in soil
C     OXYG1,OXYS1=gaseous,aqueous O2 in soil
C     CH4G1,CH4S1=gaseous,aqueous CH4 in soil
C     Z2OG1,Z2OS1=gaseous,aqueous N2O in soil
C     ZH3G1,ZH3S1,ZH3B1=gaseous,aqueous NH3 in soil non-band,band
C     H2GG1,H2GS1=gaseous,aqueous H2 in soil
C     RUPOSX=root aqueous O2 uptake
C     ROXYLX=root net O2 aqueous flux
C     RCO2SX=root aqueous CO2 uptake
C     RUPCSX=root aqueous CH4 uptake
C     RUPZSX=root aqueous N2O uptake
C     RUPNSX=root aqueous NH3 uptake non-band
C     RUPNBX=root aqueous NH3 uptake band
C     RUPHGX=root aqueous H2 uptake
C     VOLWMM,VOLPMM=soil micropore water,air volume
C     VOLW*=VOLWMM*gas solubility
C
      IF(THETPM(M,L,NY,NX).GT.THETX)THEN
      DFGSP=FPQ(N,L,NZ)*DFGS(M,L,NY,NX)
      RCODFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),CO2G1)*VOLWCO
     2-(AMAX1(ZEROS(NY,NX),CO2S1)-RCO2SX)*VOLPMM)/(VOLWCO+VOLPMM)
      RUPOST=RUPOSX-ROXYLX
      ROXDFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),OXYG1)*VOLWOX
     2-(AMAX1(ZEROS(NY,NX),OXYS1)-RUPOST)*VOLPMM)/(VOLWOX+VOLPMM)
      IF(N.EQ.1)THEN
      RCHDFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),CH4G1)*VOLWCH
     2-(AMAX1(ZEROS(NY,NX),CH4S1)-RUPCSX)*VOLPMM)/(VOLWCH+VOLPMM)
      RN2DFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),Z2OG1)*VOLWN2
     2-(AMAX1(ZEROS(NY,NX),Z2OS1)-RUPZSX)*VOLPMM)/(VOLWN2+VOLPMM)
      IF(VOLWNH+VOLPNH.GT.ZEROP(NZ,NY,NX))THEN
      ZH3GA=ZH3G1*VLNH4(L,NY,NX)
      RNHDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX
     2,DFGSP*(AMAX1(ZEROP(NZ,NY,NX),ZH3GA)*VOLWNH
     3-(AMAX1(ZEROS(NY,NX),ZH3S1)-RUPNSX)*VOLPNH)/(VOLWNH+VOLPNH)))
      ELSE
      RNHDFQ=0.0
      ENDIF
      IF(VOLWNB+VOLPNB.GT.ZEROP(NZ,NY,NX))THEN
      ZH3GB=ZH3G1*VLNHB(L,NY,NX)
      RNBDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX
     2,DFGSP*(AMAX1(ZEROP(NZ,NY,NX),ZH3GB)*VOLWNB
     3-(AMAX1(ZEROS(NY,NX),ZH3B1)-RUPNBX)*VOLPNB)/(VOLWNB+VOLPNB)))
      ELSE
      RNBDFQ=0.0
      ENDIF
      RHGDFQ=DFGSP*(AMAX1(ZEROP(NZ,NY,NX),H2GG1)*VOLWHG
     2-(AMAX1(ZEROS(NY,NX),H2GS1)-RUPHGX)*VOLPMM)/(VOLWHG+VOLPMM)
      ELSE
      RCHDFQ=0.0
      RN2DFQ=0.0
      RNHDFQ=0.0
      RNBDFQ=0.0
      RHGDFQ=0.0
      ENDIF
      ELSE
      RCODFQ=0.0
      ROXDFQ=0.0
      RCHDFQ=0.0
      RN2DFQ=0.0
      RNHDFQ=0.0
      RNBDFQ=0.0
      RHGDFQ=0.0
      ENDIF
C
C     UPDATE GASEOUS, AQUEOUS GAS CONTENTS AND CONCENTRATIONS
C     FROM GASEOUS-AQUEOUS EXCHANGE, SOIL GAS TRANSFERS
C
      OXYG1=OXYG1-ROXDFQ+ROXYFX
      OXYS1=OXYS1+ROXDFQ-RUPOSX
      CO2G1=CO2G1-RCODFQ+RCO2FX
      CO2S1=CO2S1+RCODFQ-RCO2SX
      CH4S1=CH4S1+RCHDFQ-RUPCSX
      Z2OS1=Z2OS1+RN2DFQ-RUPZSX
      ZH3S1=ZH3S1+RNHDFQ-RUPNSX
      ZH3B1=ZH3B1+RNBDFQ-RUPNBX
      H2GS1=H2GS1+RHGDFQ-RUPHGX
C
C     GAS TRANSFER THROUGH ROOTS
C
      IF(N.EQ.1.AND.RTVLP(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RUPNTX=RUPNSX+RUPNBX
C
C     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN ROOTS
C     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
C     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENT
C
C     R*DF1=root gas exchange between gaseous-aqueous phases
C     R*FL1=root gas exchange with atmosphere
C     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
C     CO2A1,CO2P1=gaseous,aqueous CO2 in root
C     OXYA1,OXYP1=gaseous,aqueous O2 in root
C     CH4A1,CH4P1=gaseous,aqueous CH4 in root
C     Z2OA1,Z2OP1=gaseous,aqueous N2O in root
C     ZH3A1,ZH3P1=gaseous,aqueous NH3 in root
C     H2GA1,H2GP1=gaseous,aqueous H2 in root
C     RTVLW,RTVLP=root aqueous,gaseous volume
C     VOLW*=RTVLW*gas solubility
C     C*E,C*A1=atmosphere,root gas concentration
C     DF*A=root-atmosphere gas conductance
C
      CO2PX=CO2P1+RCO2PX
      RCODF1=AMAX1(-CO2PX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),CO2A1)*VOLWCA
     2-CO2PX*RTVLP(N,L,NZ,NY,NX))/(VOLWCA+RTVLP(N,L,NZ,NY,NX)))
      OXYPX=OXYP1-RUPOPX
      ROXDF1=AMAX1(-OXYPX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),OXYA1)*VOLWOA
     2-OXYPX*RTVLP(N,L,NZ,NY,NX))/(VOLWOA+RTVLP(N,L,NZ,NY,NX)))
      CH4PX=CH4P1+RUPCSX
      RCHDF1=AMAX1(-CH4PX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),CH4A1)*VOLWC4
     2-CH4PX*RTVLP(N,L,NZ,NY,NX))/(VOLWC4+RTVLP(N,L,NZ,NY,NX)))
      Z2OPX=Z2OP1+RUPZSX
      RN2DF1=AMAX1(-Z2OPX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),Z2OA1)*VOLWZA
     2-Z2OPX*RTVLP(N,L,NZ,NY,NX))/(VOLWZA+RTVLP(N,L,NZ,NY,NX)))
      ZH3PX=ZH3P1+RUPNTX
      RNHDF1=AMAX1(-ZH3PX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),ZH3A1)*VOLWNA
     2-ZH3PX*RTVLP(N,L,NZ,NY,NX))/(VOLWNA+RTVLP(N,L,NZ,NY,NX)))
      H2GPX=H2GP1+RUPHGX
      RHGDF1=AMAX1(-H2GPX,DFGP*(AMAX1(ZEROP(NZ,NY,NX),H2GA1)*VOLWH2
     2-H2GPX*RTVLP(N,L,NZ,NY,NX))/(VOLWH2+RTVLP(N,L,NZ,NY,NX)))
      RCOFL1=AMIN1(DFCOA,RTVLP(N,L,NZ,NY,NX))*(CCO2E(NY,NX)-CCO2A1)
      ROXFL1=AMIN1(DFOXA,RTVLP(N,L,NZ,NY,NX))*(COXYE(NY,NX)-COXYA1)
      RCHFL1=AMIN1(DFCHA,RTVLP(N,L,NZ,NY,NX))*(CCH4E(NY,NX)-CCH4A1)
      RN2FL1=AMIN1(DFN2A,RTVLP(N,L,NZ,NY,NX))*(CZ2OE(NY,NX)-CZ2OA1)
      RNHFL1=AMIN1(DFNHA,RTVLP(N,L,NZ,NY,NX))*(CNH3E(NY,NX)-CNH3A1)
      RHGFL1=AMIN1(DFHGA,RTVLP(N,L,NZ,NY,NX))*(CH2GE(NY,NX)-CH2GA1)
      ELSE
      RCODF1=0.0
      ROXDF1=0.0
      RCHDF1=0.0
      RN2DF1=0.0
      RNHDF1=0.0
      RHGDF1=0.0
      RCOFL1=0.0
      ROXFL1=0.0
      RCHFL1=0.0
      RN2FL1=0.0
      RNHFL1=0.0
      RHGFL1=0.0
      ENDIF
C
C     UPDATE ROOT AQUEOUS, GASEOUS GAS CONTENTS AND CONCENTRATIONS
C     FOR ROOT AQUEOUS-GASEOUS, GASEOUS-ATMOSPHERE EXCHANGES
C
      CO2A1=CO2A1-RCODF1+RCOFL1
      OXYA1=OXYA1-ROXDF1+ROXFL1
      CH4A1=CH4A1-RCHDF1+RCHFL1
      Z2OA1=Z2OA1-RN2DF1+RN2FL1
      ZH3A1=ZH3A1-RNHDF1+RNHFL1
      H2GA1=H2GA1-RHGDF1+RHGFL1
      CO2P1=CO2P1+RCODF1+RCO2SX+RCO2PX
      OXYP1=OXYP1+ROXDF1-RUPOPX
      CH4P1=CH4P1+RCHDF1+RUPCSX
      Z2OP1=Z2OP1+RN2DF1+RUPZSX
      ZH3P1=ZH3P1+RNHDF1+RUPNSX+RUPNBX
      H2GP1=H2GP1+RHGDF1+RUPHGX
C
C     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
C
C     RCO2S=soil-root CO2 exchange
C     RUPOXS=soil-root O2 exchange
C     RUPCHS=soil-root CH4 exchange
C     RUPN2S=soil-root N2O exchange
C     RUPN3S=soil-root NH3 exchange non-band
C     RUPN3B=soil-root NH3 exchange band
C     RUPHGS=soil-root H2 exchange
C
      RCO2S(N,L,NZ,NY,NX)=RCO2S(N,L,NZ,NY,NX)+RCO2SX
      RUPOXS(N,L,NZ,NY,NX)=RUPOXS(N,L,NZ,NY,NX)+RUPOSX
      RUPCHS(N,L,NZ,NY,NX)=RUPCHS(N,L,NZ,NY,NX)+RUPCSX
      RUPN2S(N,L,NZ,NY,NX)=RUPN2S(N,L,NZ,NY,NX)+RUPZSX
      RUPN3S(N,L,NZ,NY,NX)=RUPN3S(N,L,NZ,NY,NX)+RUPNSX
      RUPN3B(N,L,NZ,NY,NX)=RUPN3B(N,L,NZ,NY,NX)+RUPNBX
      RUPHGS(N,L,NZ,NY,NX)=RUPHGS(N,L,NZ,NY,NX)+RUPHGX
C     IF(NZ.EQ.1.AND.L.EQ.7.AND.N.EQ.1)THEN
C     WRITE(*,5547)'RCO2S',I,J,NX,NY,NZ,L,N,RCO2S(N,L,NZ,NY,NX)
C    2,RCO2SX,RDXCOS,RDFCOS,RMFCOS,DIFCL,CCO2S1,CCO2P1,VOLWSP,VOLWMM
C     WRITE(*,5547)'RCH4S',I,J,NX,NY,NZ,L,N,RUPCHS(N,L,NZ,NY,NX)
C    2,RUPCSX,RDXCHS,RDFCHS,RMFCHS,DIFCL,CCH4S1,CCH4P1,VOLWSP,VOLWMM
C     WRITE(*,5547)'RUPNBX',I,J,NX,NY,NZ,L,N,RUPNBX,RDXNHB,RDFN3B
C    2,RTVLWB,ZH3B1,VOLWMB,ZH3PB,VOLWSB,RUPN3B(N,L,NZ,NY,NX)
C    3,RNBDFQ,RNHDF1,RUPNSX
C     WRITE(*,5547)'RUPN2S',I,J,NX,NY,NZ,L,N,RUPN2S(N,L,NZ,NY,NX)
C    2,RUPZSX,RDXN2S,RDFN2S,RTVLW(N,L,NZ,NY,NX),Z2OS1,Z2OP1,VOLWMM
C    2,RMFN2S,DIFZL,CN2OS1,CN2OP1,PP(NZ,NY,NX)
5547  FORMAT(A8,7I4,30E12.4)
C     ENDIF
C
C     ACCUMULATE ROOT-ATMOSPHERE GAS EXCHANGE TO HOURLY TIME SCALE
C
C     R*DFA=root aqueous-gaseous CO2 exchange
C     R*FLA=root gaseous-atmosphere CO2 exchange
C     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
C
      RCODFA(N,L,NZ,NY,NX)=RCODFA(N,L,NZ,NY,NX)+RCODF1
      ROXDFA(N,L,NZ,NY,NX)=ROXDFA(N,L,NZ,NY,NX)+ROXDF1
      RCHDFA(N,L,NZ,NY,NX)=RCHDFA(N,L,NZ,NY,NX)+RCHDF1
      RN2DFA(N,L,NZ,NY,NX)=RN2DFA(N,L,NZ,NY,NX)+RN2DF1
      RNHDFA(N,L,NZ,NY,NX)=RNHDFA(N,L,NZ,NY,NX)+RNHDF1
      RHGDFA(N,L,NZ,NY,NX)=RHGDFA(N,L,NZ,NY,NX)+RHGDF1
      RCOFLA(N,L,NZ,NY,NX)=RCOFLA(N,L,NZ,NY,NX)+RCOFL1
      ROXFLA(N,L,NZ,NY,NX)=ROXFLA(N,L,NZ,NY,NX)+ROXFL1
      RCHFLA(N,L,NZ,NY,NX)=RCHFLA(N,L,NZ,NY,NX)+RCHFL1
      RN2FLA(N,L,NZ,NY,NX)=RN2FLA(N,L,NZ,NY,NX)+RN2FL1
      RNHFLA(N,L,NZ,NY,NX)=RNHFLA(N,L,NZ,NY,NX)+RNHFL1
      RHGFLA(N,L,NZ,NY,NX)=RHGFLA(N,L,NZ,NY,NX)+RHGFL1
C
C     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
C
C     RCO2P=root CO2 emission into root
C     RUPOXP=root O2 uptake from root
C     ROXSK=total O2 uptake from soil by all microbial,root popns
C
      RCO2P(N,L,NZ,NY,NX)=RCO2P(N,L,NZ,NY,NX)+RCO2PX+RCO2SX
      RUPOXP(N,L,NZ,NY,NX)=RUPOXP(N,L,NZ,NY,NX)+RUPOPX
      ROXSK(M,L,NY,NX)=ROXSK(M,L,NY,NX)+RUPOSX
C     IF(NZ.EQ.1.AND.L.EQ.7.AND.N.EQ.1)THEN
C     WRITE(*,5566)'CO2P1',I,J,NX,NY,NZ,L,N,M,MX,RCO2S(N,L,NZ,NY,NX)
C    2,RCO2SX,RDFCOS,RDXCOS,RMFCOS,DIFCL,DFGS(M,L,NY,NX)
C    3,CCO2S1,CCO2P1,RTVLW(N,L,NZ,NY,NX),CO2S1,VOLWMM
C    4,CO2P1,VOLWSP,PP(NZ,NY,NX),FPQ(N,L,NZ),RCODF1,RCO2PX
C    5,CO2PX,RTVLP(N,L,NZ,NY,NX),DFGP,VOLWCA,CO2A1
C    6,XNPD,PORT(N,NZ,NY,NX),TFND(L,NY,NX),RCO2FX
C    7,RCODFQ,DFGSP,FOXYX,CO2G1,VOLWCO,VOLPMM,THETPM(M,L,NY,NX),THETX
C    8,ROXYP(N,L,NZ,NY,NX),ROXYY(L,NY,NX)
C     WRITE(*,5566)'OXYP1',I,J,NX,NY,NZ,L,N,M,MX,UPMXP*PP(NZ,NY,NX)
C    2,RUPOSX,ROXDFQ,OXYS1,RUPOPX,ROXDF1,OXYP1
C    3,FOXYX,DFGS(M,L,NY,NX),DFGP,ROXYFX,ROXYLX,ROXFL1
C    3,OXYG1,OXYA1,COXYS1,COXYP1,COXYR,ROXSK(M,L,NY,NX),XS,XR
C    4,OXYPY,VOLWOA,RTVLP(N,L,NZ,NY,NX),RTVLW(N,L,NZ,NY,NX)
C    5,DFOXA,COXYE(NY,NX),COXYA1,RUPOXP(N,L,NZ,NY,NX)
C    6,RUPOXS(N,L,NZ,NY,NX),ROXYP(N,L,NZ,NY,NX),THETPM(M,L,NY,NX)
C     WRITE(*,5566)'CH4S1',I,J,NX,NY,NZ,L,N,M,MX,CH4S1,RCHDFQ,RUPCSX
C    2,RDFCHS,RMFCHS,DIFCL,CCH4S1,CCH4P1,CH4P1,CH4PX,RCHDF1,RCHFL1
C    3,DFCHA,RTVLP(N,L,NZ,NY,NX),CCH4E(NY,NX),CCH4A1,THETM,CQSGL1,RTARRX
C    4,RRADS,FILM(M,L,NY,NX),RRADL(N,L),RRAD2M(N,NZ,NY,NX)
C    5,RTDNP(N,L,NZ,NY,NX)
C    4,DFGS(M,L,NY,NX),CH4G1,VOLWCH,CH4S1,VOLPMM,THETPM(M,L,NY,NX)
5566  FORMAT(A8,9I4,50E12.4)
C     ENDIF
90    CONTINUE
      ENDIF
99    CONTINUE
C
C     O2 CONSTRAINTS TO ROOT RESPIRATION DEPENDS UPON RATIO
C     OF ROOT O2 UPTAKE 'RUPOXT' TO ROOT O2 DEMAND 'ROXYP'
C
C     RUPOXT=O2 uptake from soil+root by each root,myco population
C     ROXYP=O2 demand by each root,myco population
C     WFR=constraint by O2 consumption on all root processes
C     imposed by O2 uptake
C
      RUPOXT=RUPOXP(N,L,NZ,NY,NX)+RUPOXS(N,L,NZ,NY,NX)
      WFR(N,L,NZ,NY,NX)=AMIN1(1.0,AMAX1(0.0
     2,RUPOXT/ROXYP(N,L,NZ,NY,NX)))
C     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,3368)'WFR',I,J,NX,NY,NZ,L,N,WFR(N,L,NZ,NY,NX)
C    2,RUPOXP(N,L,NZ,NY,NX),RUPOXS(N,L,NZ,NY,NX)
C    3,ROXYP(N,L,NZ,NY,NX)
3368  FORMAT(A8,7I4,12E12.4)
C     ENDIF
      ELSE
      RUPOXT=0.0
      IF(L.GT.NG(NZ,NY,NX))THEN
      WFR(N,L,NZ,NY,NX)=WFR(N,L-1,NZ,NY,NX)
      ELSE
      WFR(N,L,NZ,NY,NX)=1.0
      ENDIF
      ENDIF
      OSTRD=OSTRD+ROXYP(N,L,NZ,NY,NX)
      OSTRN=OSTRN+RUPOXT
C
C     ROOT EXUDATION OF C, N AND P DEPENDS ON CONCN DIFFERENCES
C     BETWEEN ROOT NON-STRUCTURAL POOLS AND SOIL DISSOLVED POOLS
C
C     VOLWMM=soil micropore water volume
C     FOSRH=fraction of total SOC in each substrate K from nitro.f
C     RTVLW=root aqueous volume
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P in root,myco
C     XFRC,XFRN,XFRP=nonstructural C,N,P exchg at root-soil DOC equilibrium
C     OQC=soil DOC
C     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
C     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation
C     TLEC,TSHC=total fluxes x blr for calculating canopy air temperature,
C     vapor pressure in watsub.f
C     EFLXC,SFLXC=canopylatent,sensible heat fluxes
C     RA=canopy boundary layer resistance
C     OSTR=O2 stress indicator
C
      DO 195 K=0,4
      VOLWK=VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)
      IF(VOLWK.GT.ZEROS2(NY,NX)
     2.AND.RTVLW(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      VOLWT=VOLWK+RTVLW(N,L,NZ,NY,NX)
      CPOOLX=AMIN1(1.25E+03*RTVLW(N,L,NZ,NY,NX),CPOOLR(N,L,NZ,NY,NX))
      XFRC=(OQC(K,L,NY,NX)*RTVLW(N,L,NZ,NY,NX)
     2-CPOOLX*VOLWK)/VOLWT
      RDFOMC(N,K,L,NZ,NY,NX)=FEXUC*XFRC
      IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX)
     2.AND.CPOOLR(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CPOOLT=OQC(K,L,NY,NX)+CPOOLR(N,L,NZ,NY,NX)
      ZPOOLX=0.1*ZPOOLR(N,L,NZ,NY,NX)
      PPOOLX=0.1*PPOOLR(N,L,NZ,NY,NX)
      XFRN=(OQN(K,L,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
     2-ZPOOLX*OQC(K,L,NY,NX))/CPOOLT
      XFRP=(OQP(K,L,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
     2-PPOOLX*OQC(K,L,NY,NX))/CPOOLT
      RDFOMN(N,K,L,NZ,NY,NX)=FEXUN*XFRN
      RDFOMP(N,K,L,NZ,NY,NX)=FEXUP*XFRP
      ELSE
      RDFOMN(N,K,L,NZ,NY,NX)=0.0
      RDFOMP(N,K,L,NZ,NY,NX)=0.0
      ENDIF
      ELSE
      RDFOMC(N,K,L,NZ,NY,NX)=0.0
      RDFOMN(N,K,L,NZ,NY,NX)=0.0
      RDFOMP(N,K,L,NZ,NY,NX)=0.0
      ENDIF
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.1)THEN
C     WRITE(*,2224)'RDFOMC',I,J,NX,NY,L,NZ,K,N,RDFOMC(N,K,L,NZ,NY,NX)
C    2,RDFOMN(N,K,L,NZ,NY,NX),RDFOMP(N,K,L,NZ,NY,NX)
C    3,OQC(K,L,NY,NX),OQN(K,L,NY,NX),OQP(K,L,NY,NX)
C    2,CPOOLR(N,L,NZ,NY,NX),ZPOOLR(N,L,NZ,NY,NX),PPOOLR(N,L,NZ,NY,NX)
C    3,VOLWM(NPH,L,NY,NX),RTVLW(N,L,NZ,NY,NX),RTAR1X(N,NZ,NY,NX)
C    4,RTAR2X(N,NZ,NY,NX),RTLGP(N,L,NZ,NY,NX)*PP(NZ,NY,NX)
C    4,WTRTD(N,L,NZ,NY,NX)
C    5,VOLWK,VOLWM(NPH,L,NY,NX),FOSRH(K,L,NY,NX)
C    5,OQC(K,L,NY,NX)/VOLWK
C    5,OQN(K,L,NY,NX)/OQC(K,L,NY,NX)
C    5,OQP(K,L,NY,NX)/OQC(K,L,NY,NX)
C    6,CPOOLR(N,L,NZ,NY,NX)/RTVLW(N,L,NZ,NY,NX)
C    6,ZPOOLX/CPOOLR(N,L,NZ,NY,NX)
C    6,PPOOLX/CPOOLR(N,L,NZ,NY,NX)
2224  FORMAT(A8,8I4,30E12.4)
C     ENDIF
195   CONTINUE
C
C     NUTRIENT UPTAKE
C
C     WFR=constraint by O2 consumption on all biological processes
C     FCUP=limitation to active uptake respiration from CPOOLR
C     FWSRT=protein concentration relative to 5%
C     RTLGP=root,myco length per plant
C
      IF(WFR(N,L,NZ,NY,NX).GT.ZERO
     2.AND.FCUP.GT.ZERO.AND.FWSRT.GT.ZERO
     3.AND.RTLGP(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C
C     FZUP=limitn to active uptake respiration from CZPOLR
C
      IF(FZUP.GT.ZERO2)THEN
C
C     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NH4,NO3
C     FROM SOIL TO ROOT
C
C     ZNSGL=NH4 diffusivity
C     TORT=soil tortuosity
C     PATH=path length of water and nutrient uptake
C     RRADL=root radius
C     DIFFL=NH4 diffusion per plant
C
      ZNSGX=ZNSGL(L,NY,NX)*TORT(NPH,L,NY,NX)
      PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*ZNSGX))
      DIFFL=ZNSGX*RTARR(N,L)/LOG(PATHL/RRADL(N,L))
C
C     NH4 UPTAKE IN NON-BAND SOIL ZONE
C
C     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
C     CNH4S=NH4 concentration in non-band
C     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
C     UPWTRP=root water uptake per plant
C     RMFNH4=soil-root convective NH4 flux per plant in non-band
C     DIFNH4=soil-root NH4 diffusion per plant in non-band
C
      IF(VLNH4(L,NY,NX).GT.ZERO.AND.CNH4S(L,NY,NX)
     2.GT.UPMNZH(N,NZ,NY,NX))THEN
      RMFNH4=UPWTRP*CNH4S(L,NY,NX)*VLNH4(L,NY,NX)
      DIFNH4=DIFFL*VLNH4(L,NY,NX)
C
C     NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
C     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
C
C     UPMXP,UPMX=max NH4 uptake in non-band unlimited,limited by O2
C     RTARP=root surface area per plant from grosub.f
C     FWSRT=protein concentration relative to 5%
C     TFN4=temperature function for root growth
C     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
C     WFR=constraint by O2 consumption on all biological processes
C
      UPMXP=UPMXZH(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLNH4(L,NY,NX)*AMIN1(FCUP,FZUP)
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFNH4=soil-root convective NH4 flux per plant in non-band
C     DIFNH4=soil-root NH4 diffusion per plant in non-band
C     CNH4S=NH4 concentration in non-band
C     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
C     RTKNH4,RTKNHP=NH4 uptake per plant in non-band lmtd,unlmtd by O2
C     ZNH4M,ZNH4X=minimum,maximum NH4 available for uptake in non-band
C     FNH4X=fraction of total NH4 uptake in non-band by root,myco populn
C     RUNNHP,RUPNH4=NH4 uptake in non-band unlimited,limited by NH4
C     RUONH4=NH4 uptake in non-band unlimited by O2
C     RUCNH4=NH4 uptake in non-band unlimited by nonstructural C
C
      X=(DIFNH4+RMFNH4)*CNH4S(L,NY,NX)
      Y=DIFNH4*UPMNZH(N,NZ,NY,NX)
      B=-UPMX-DIFNH4*UPKMZH(N,NZ,NY,NX)-X+Y
      C=(X-Y)*UPMX
      RTKNH4=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFNH4*UPKMZH(N,NZ,NY,NX)-X+Y
      CP=(X-Y)*UPMXP
      RTKNHP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      ZNH4M=UPMNZH(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLNH4(L,NY,NX)
      ZNH4X=AMAX1(0.0,FNH4X*(ZNH4S(L,NY,NX)-ZNH4M))
      RUNNHP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNH4*PP(NZ,NY,NX))
      RUPNH4(N,L,NZ,NY,NX)=AMIN1(ZNH4X,RUNNHP(N,L,NZ,NY,NX))
      RUONH4(N,L,NZ,NY,NX)=AMIN1(ZNH4X,AMAX1(0.0
     2,RTKNHP*PP(NZ,NY,NX)))
      RUCNH4(N,L,NZ,NY,NX)=RUPNH4(N,L,NZ,NY,NX)/FCUP
C     IF(NX.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,1110)'UPNH4',I,J,NZ,L,N,RUNNHP(N,L,NZ,NY,NX)
C    2,RUPNH4(N,L,NZ,NY,NX),RTKNH4,RMFNH4,X,Y,B,C,UPMX,UPMXP
C    2,WFR(N,L,NZ,NY,NX),CNH4S(L,NY,NX),DIFNH,RTDNP(N,L,NZ,NY,NX)
C    2,WTRTD(N,L,NZ,NY,NX),CNH4S(L,NY,NX),RDFOMN(N,L,NZ,NY,NX)
C    3,CCPOLR(N,L,NZ,NY,NX),CZPOLR(N,L,NZ,NY,NX),CPPOLR(N,L,NZ,NY,NX)
C    4,THETW(L,NY,NX),TKS(L,NY,NX),RSCS(L,NY,NX),UPMXP,FWSRT
C    5,FZUP,FCUP,COXYS(L,NY,NX),COXYG(L,NY,NX)
C    6,CCPOLP(NZ,NY,NX)
C    7,CZPOLP(NZ,NY,NX),CPPOLP(NZ,NY,NX),FDBK(1,NZ,NY,NX),PSIST1(L)
C    2,PSIRT(N,L,NZ,NY,NX),ZPOOLR(N,L,NZ,NY,NX),WTRTL(N,L,NZ,NY,NX)
C    3,RTARP(N,L,NZ,NY,NX),RRADL(N,L),PATH(N,L)
C    4,DIFFL,ZNSGX,RTARR(N,L),PATHL,RRADL(N,L),VLNH4(L,NY,NX)
1110  FORMAT(A8,5I4,100E24.16)
C     ENDIF
      ELSE
      RUNNHP(N,L,NZ,NY,NX)=0.0
      RUPNH4(N,L,NZ,NY,NX)=0.0
      RUONH4(N,L,NZ,NY,NX)=0.0
      RUCNH4(N,L,NZ,NY,NX)=0.0
      ENDIF
C
C     NH4 UPTAKE IN BAND SOIL ZONE
C
C     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
C     CNH4B=NH4 concentration in band
C     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
C     UPWTRP=root water uptake per plant
C     RMFNHB=soil-root convective NH4 flux per plant in band
C     DIFNHB=soil-root NH4 diffusion per plant in band
C

      IF(VLNHB(L,NY,NX).GT.ZERO.AND.CNH4B(L,NY,NX)
     2.GT.UPMNZH(N,NZ,NY,NX))THEN
      RMFNHB=UPWTRP*CNH4B(L,NY,NX)*VLNHB(L,NY,NX)
      DIFNHB=DIFFL*VLNHB(L,NY,NX)
C
C     NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
C     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
C
C     UPMXP,UPMX=maximum NH4 uptake in band unlimited,limited by O2
C     RTARP=root surface area per plant from grosub.f
C     FWSRT=protein concentration relative to 5%
C     TFN4=temperature function for root growth
C     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
C     WFR=constraint by O2 consumption on all biological processes
C
      UPMXP=UPMXZH(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLNHB(L,NY,NX)*AMIN1(FCUP,FZUP)
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFNHB=soil-root convective NH4 flux per plant in band
C     DIFNHB=soil-root NH4 diffusion per plant in band
C     CNH4B=NH4 concentration in band
C     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
C     RTKNHB,RTKNBP=NH4 uptake per plant in band lmtd,unlmtd by O2
C     ZNHBM,ZNHBX=minimum,maximum NH4 available for uptake in band
C     FNHBX=fraction of total NH4 uptake in band by root,myco populn
C     RUNNBP,RUPNHB=NH4 uptake in band unlimited,limited by NH4
C     RUONHB=NH4 uptake in band unlimited by O2
C     RUCNHB=NH4 uptake in band unlimited by nonstructural C
C
      X=(DIFNHB+RMFNHB)*CNH4B(L,NY,NX)
      Y=DIFNHB*UPMNZH(N,NZ,NY,NX)
      B=-UPMX-DIFNHB*UPKMZH(N,NZ,NY,NX)-X+Y
      C=(X-Y)*UPMX
      RTKNHB=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFNHB*UPKMZH(N,NZ,NY,NX)-X+Y
      CP=(X-Y)*UPMXP
      RTKNBP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      ZNHBM=UPMNZH(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLNHB(L,NY,NX)
      ZNHBX=AMAX1(0.0,FNHBX*(ZNH4B(L,NY,NX)-ZNHBM))
      RUNNBP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNHB*PP(NZ,NY,NX))
      RUPNHB(N,L,NZ,NY,NX)=AMIN1(ZNHBX,RUNNBP(N,L,NZ,NY,NX))
      RUONHB(N,L,NZ,NY,NX)=AMIN1(ZNHBX,AMAX1(0.0
     2,RTKNBP*PP(NZ,NY,NX)))
      RUCNHB(N,L,NZ,NY,NX)=RUPNHB(N,L,NZ,NY,NX)/FCUP
      ELSE
      RUNNBP(N,L,NZ,NY,NX)=0.0
      RUPNHB(N,L,NZ,NY,NX)=0.0
      RUONHB(N,L,NZ,NY,NX)=0.0
      RUCNHB(N,L,NZ,NY,NX)=0.0
      ENDIF
C
C     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NO3
C     FROM SOIL TO ROOT
C
C     ZOSGL=NO3 diffusivity
C     TORT=soil tortuosity
C     RRADL=root radius
C     PATH=path length of water and nutrient uptake
C     DIFFL=NO3 diffusion per plant
C
      ZOSGX=ZOSGL(L,NY,NX)*TORT(NPH,L,NY,NX)
      PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*ZOSGX))
      DIFFL=ZOSGX*RTARR(N,L)/LOG(PATHL/RRADL(N,L))
C
C     NO3 UPTAKE IN NON-BAND SOIL ZONE
C
C     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
C     CNO3S=NO3 concentration in non-band
C     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
C     UPWTRP=root water uptake per plant
C     RMFNO3=soil-root convective NO3 flux per plant in non-band
C     DIFNO3=soil-root NO3 diffusion per plant in non-band
C
      IF(VLNO3(L,NY,NX).GT.ZERO.AND.CNO3S(L,NY,NX)
     2.GT.UPMNZO(N,NZ,NY,NX))THEN
      RMFNO3=UPWTRP*CNO3S(L,NY,NX)*VLNO3(L,NY,NX)
      DIFNO3=DIFFL*VLNO3(L,NY,NX)
C
C     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
C     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
C
C     UPMXP,UPMX=max NO3 uptake in non-band unlimited,limited by O2
C     RTARP=root surface area per plant from grosub.f
C     FWSRT=protein concentration relative to 5%
C     TFN4=temperature function for root growth
C     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
C     WFR=constraint by O2 consumption on all biological processes
C
      UPMXP=UPMXZO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLNO3(L,NY,NX)*AMIN1(FCUP,FZUP)
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFNO3=soil-root convective N03 flux per plant in non-band
C     DIFNO3=soil-root N03 diffusion per plant in non-band
C     CNO3S=NO3 concentration in non-band
C     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
C     RTKNO3,RTKNOP=NO3 uptake per plant in non-band lmtd,unlmtd by O2
C     ZNO3M,ZNO3X=minimum,maximum NO3 available for uptake in non-band
C     FNO3X=fraction of total NH4 uptake in non-band by root,myco populn
C     RUNNOP,RUPNO3=NO3 uptake in non-band unlimited,limited by NO3
C     RUONO3=NO3 uptake in non-band unlimited by O2
C     RUCNO3=NO3 uptake in non-band unlimited by nonstructural C
C
      X=(DIFNO3+RMFNO3)*CNO3S(L,NY,NX)
      Y=DIFNO3*UPMNZO(N,NZ,NY,NX)
      B=-UPMX-DIFNO3*UPKMZO(N,NZ,NY,NX)-X+Y
      C=(X-Y)*UPMX
      RTKNO3=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFNO3*UPKMZO(N,NZ,NY,NX)-X+Y
      CP=(X-Y)*UPMXP
      RTKNOP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      ZNO3M=UPMNZO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLNO3(L,NY,NX)
      ZNO3X=AMAX1(0.0,FNO3X*(ZNO3S(L,NY,NX)-ZNO3M))
      RUNNOP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNO3*PP(NZ,NY,NX))
      RUPNO3(N,L,NZ,NY,NX)=AMIN1(ZNO3X,RUNNOP(N,L,NZ,NY,NX))
      RUONO3(N,L,NZ,NY,NX)=AMIN1(ZNO3X,AMAX1(0.0
     2,RTKNOP*PP(NZ,NY,NX)))
      RUCNO3(N,L,NZ,NY,NX)=RUPNO3(N,L,NZ,NY,NX)/FCUP
C     IF(NX.EQ.4.AND.NY.EQ.2)THEN
C     WRITE(*,1111)'UPNO3',I,J,NZ,L,N,RUPNO3(N,L,NZ,NY,NX),FNO3X
C    2,ZNO3S(L,NY,NX),ZNO3M,RTDNP(N,L,NZ,NY,NX),RTKNO3,RMFNO3,X,Y,B,C
C    2,UPMX,CNO3S(L,NY,NX),DIFNO,RUONO3(N,L,NZ,NY,NX)
C    3,CCPOLR(N,L,NZ,NY,NX),CZPOLR(N,L,NZ,NY,NX),CPPOLR(N,L,NZ,NY,NX)
C    4,THETW(L,NY,NX),TKS(L,NY,NX),RSCS(L,NY,NX),UPMXP,FWSRT
C    5,FZUP,FCUP,COXYS(L,NY,NX),COXYG(L,NY,NX),WFR(N,L,NZ,NY,NX)
C    6,CCPOLP(NZ,NY,NX),CZPOLP(NZ,NY,NX),CPPOLP(NZ,NY,NX)
C    7,FDBK(1,NZ,NY,NX),PSIST1(L),PSIRT(N,L,NZ,NY,NX)
C    2,ZPOOLR(N,L,NZ,NY,NX),WTRTL(N,L,NZ,NY,NX)
C    3,RUNNOP(N,L,NZ,NY,NX),RNO3Y(L,NY,NX)
1111  FORMAT(A8,5I4,40E12.4)
C     ENDIF
      ELSE
      RUNNOP(N,L,NZ,NY,NX)=0.0
      RUPNO3(N,L,NZ,NY,NX)=0.0
      RUONO3(N,L,NZ,NY,NX)=0.0
      RUCNO3(N,L,NZ,NY,NX)=0.0
      ENDIF
C
C     NO3 UPTAKE IN BAND SOIL ZONE
C
C     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
C     CNO3B=NO3 concentration in band
C     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
C     UPWTRP=root water uptake per plant
C     RMFNOB=soil-root convective NO3 flux per plant in band
C     DIFNOB=soil-root NO3 diffusion per plant in band
C

      IF(VLNOB(L,NY,NX).GT.ZERO.AND.CNO3B(L,NY,NX)
     2.GT.UPMNZO(N,NZ,NY,NX))THEN
      RMFNOB=UPWTRP*CNO3B(L,NY,NX)*VLNOB(L,NY,NX)
      DIFNOB=DIFFL*VLNOB(L,NY,NX)
C
C     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
C     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
C
C     UPMXP,UPMX=maximum NO3 uptake in band unlimited,limited by O2
C     RTARP=root surface area per plant from grosub.f
C     FWSRT=protein concentration relative to 5%
C     TFN4=temperature function for root growth
C     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
C     WFR=constraint by O2 consumption on all biological processes
C
      UPMXP=UPMXZO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLNOB(L,NY,NX)*AMIN1(FCUP,FZUP)
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFNOB=soil-root convective NO3 flux per plant in band
C     DIFNOB=soil-root NO3 diffusion per plant in band
C     CNO3B=NH4 concentration in band
C     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
C     RTKNOB,RTKNPB=NO3 uptake per plant in band lmtd,unlmtd by O2
C     ZNOBM,ZNOBX=minimum,maximum NO3 available for uptake in band
C     FNOBX=fraction of total NO3 uptake in band by root,myco populn
C     RUNNXP,RUPNOB=NO3 uptake in band unlimited,limited by NH4
C     RUONOB=NO3 uptake in band unlimited by O2
C     RUCNOB=NO3 uptake in band unlimited by nonstructural C
C
      X=(DIFNOB+RMFNOB)*CNO3B(L,NY,NX)
      Y=DIFNOB*UPMNZO(N,NZ,NY,NX)
      B=-UPMX-DIFNOB*UPKMZO(N,NZ,NY,NX)-X+Y
      C=(X-Y)*UPMX
      RTKNOB=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFNOB*UPKMZO(N,NZ,NY,NX)-X+Y
      CP=(X-Y)*UPMXP
      RTKNPB=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      ZNOBM=UPMNZO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLNOB(L,NY,NX)
      ZNOBX=AMAX1(0.0,FNOBX*(ZNO3B(L,NY,NX)-ZNOBM))
      RUNNXP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNOB*PP(NZ,NY,NX))
      RUPNOB(N,L,NZ,NY,NX)=AMIN1(ZNOBX,RUNNXP(N,L,NZ,NY,NX))
      RUONOB(N,L,NZ,NY,NX)=AMIN1(ZNOBX
     2,AMAX1(0.0,RTKNPB*PP(NZ,NY,NX)))
      RUCNOB(N,L,NZ,NY,NX)=RUPNOB(N,L,NZ,NY,NX)/FCUP
      ELSE
      RUNNXP(N,L,NZ,NY,NX)=0.0
      RUPNOB(N,L,NZ,NY,NX)=0.0
      RUONOB(N,L,NZ,NY,NX)=0.0
      RUCNOB(N,L,NZ,NY,NX)=0.0
      ENDIF
      ELSE
      RUNNHP(N,L,NZ,NY,NX)=0.0
      RUPNH4(N,L,NZ,NY,NX)=0.0
      RUONH4(N,L,NZ,NY,NX)=0.0
      RUCNH4(N,L,NZ,NY,NX)=0.0
      RUNNBP(N,L,NZ,NY,NX)=0.0
      RUPNHB(N,L,NZ,NY,NX)=0.0
      RUONHB(N,L,NZ,NY,NX)=0.0
      RUCNHB(N,L,NZ,NY,NX)=0.0
      RUNNOP(N,L,NZ,NY,NX)=0.0
      RUPNO3(N,L,NZ,NY,NX)=0.0
      RUONO3(N,L,NZ,NY,NX)=0.0
      RUCNO3(N,L,NZ,NY,NX)=0.0
      RUNNXP(N,L,NZ,NY,NX)=0.0
      RUPNOB(N,L,NZ,NY,NX)=0.0
      RUONOB(N,L,NZ,NY,NX)=0.0
      RUCNOB(N,L,NZ,NY,NX)=0.0
      ENDIF
C
C     FPUP=limitn to active uptake respiration from CPPOLR
C
      IF(FPUP.GT.ZERO2)THEN
C
C     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF H2PO4,HPO4
C     FROM SOIL TO ROOT
C
C     POSGL=PO4 diffusivity
C     TORT=soil tortuosity
C     PATH=path length of water and nutrient uptake
C     RRADL=root radius
C     DIFFL=PO4 diffusion per plant
C
      POSGX=POSGL(L,NY,NX)*TORT(NPH,L,NY,NX)
      PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*POSGX))
      DIFFL=POSGX*RTARR(N,L)/LOG(PATHL/RRADL(N,L))
C
C     H2PO4 UPTAKE IN NON-BAND SOIL ZONE
C
C     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
C     CH2P4=H2PO4 concentration in non-band
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
C     UPWTRP=root water uptake per plant
C     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
C     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
C
      IF(VLPO4(L,NY,NX).GT.ZERO.AND.CH2P4(L,NY,NX)
     2.GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH2P=UPWTRP*CH2P4(L,NY,NX)*VLPO4(L,NY,NX)
      DIFH2P=DIFFL*VLPO4(L,NY,NX)
C
C     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
C     AND FROM ROOT SURFACE AREA, C AND P CONSTRAINTS CALCULATED ABOVE
C
C     UPMXP,UPMX=max H2PO4 uptake in non-band unlimited,limited by O2
C     RTARP=root surface area per plant from grosub.f
C     FWSRT=protein concentration relative to 5%
C     TFN4=temperature function for root growth
C     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
C     WFR=constraint by O2 consumption on all biological processes
C
      UPMXP=UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLPO4(L,NY,NX)*AMIN1(FCUP,FPUP)
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF H2PO4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY
C     COMPETITION WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
C     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
C     CH2P4=H2PO4 concentration in non-band
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
C     RTKH2P,RTKHPP=H2PO4 uptake per plant in non-band lmtd,unlmtd by O2
C     H2POM,H2POX=minimum,maximum H2PO4 available for uptake in non-band
C     FPO4X=fraction of total H2PO4 uptake in non-band by root,myco populn
C     RUPP2P,RUPH2P=H2PO4 uptake in non-band unlimited,limited by H2PO4
C     RUOH2P=H2PO4 uptake in non-band unlimited by O2
C     RUCH2P=H2PO4 uptake in non-band unlimited by nonstructural C
C
      X=(DIFH2P+RMFH2P)*CH2P4(L,NY,NX)
      Y=DIFH2P*UPMNPO(N,NZ,NY,NX)
      B=-UPMX-DIFH2P*UPKMPO(N,NZ,NY,NX)-X+Y
      C=(X-Y)*UPMX
      RTKH2P=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFH2P*UPKMPO(N,NZ,NY,NX)-X+Y
      CP=(X-Y)*UPMXP
      RTKHPP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      H2POM=UPMNPO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
      H2POX=AMAX1(0.0,FPO4X*(H2PO4(L,NY,NX)-H2POM))
      RUPP2P(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH2P*PP(NZ,NY,NX))
      RUPH2P(N,L,NZ,NY,NX)=AMIN1(H2POX,RUPP2P(N,L,NZ,NY,NX))
      RUOH2P(N,L,NZ,NY,NX)=AMIN1(H2POX,AMAX1(0.0
     2,RTKHPP*PP(NZ,NY,NX)))
      RUCH2P(N,L,NZ,NY,NX)=RUPH2P(N,L,NZ,NY,NX)/FCUP
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.3)THEN
C     WRITE(*,2223)'UPPO4',I,J,NZ,L,N,RUPH2P(N,L,NZ,NY,NX),FPO4X
C    2,H2PO4(L,NY,NX),RUPP2P(N,L,NZ,NY,NX),UPMX,DIFPO,UPKMPO(N,NZ,NY,NX)
C    3,UPMNPO(N,NZ,NY,NX),RMFH2P,CH2P4(L,NY,NX),UPMXP,WFR(N,L,NZ,NY,NX)
C    4,FCUP,FZUP,FPUP,UPMXPO(N,NZ,NY,NX),RTARP(N,L,NZ,NY,NX),FWSRT
C    5,TFN4(L,NZ,NY,NX),DIFFL,CPO4S(L,NY,NX),CPOOLR(N,L,NZ,NY,NX)
C    6,PPOOLR(N,L,NZ,NY,NX),RTKH2P,PP(NZ,NY,NX)
C    2,RTLGP(N,L,NZ,NY,NX)
2223  FORMAT(A8,5I4,40E12.4)
C     ENDIF
      ELSE
      RUPP2P(N,L,NZ,NY,NX)=0.0
      RUPH2P(N,L,NZ,NY,NX)=0.0
      RUOH2P(N,L,NZ,NY,NX)=0.0
      RUCH2P(N,L,NZ,NY,NX)=0.0
      ENDIF
C
C     H2PO4 UPTAKE IN BAND SOIL ZONE
C
C     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
C     CH2P4B=H2PO4 concentration in band
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
C     UPWTRP=root water uptake per plant
C     RMFH2B=soil-root convective H2PO4 flux per plant in band
C     DIFH2B=soil-root H2PO4 diffusion per plant in band
C

      IF(VLPOB(L,NY,NX).GT.ZERO.AND.CH2P4B(L,NY,NX)
     2.GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH2B=UPWTRP*CH2P4B(L,NY,NX)*VLPOB(L,NY,NX)
      DIFH2B=DIFFL*VLPOB(L,NY,NX)
C
C     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
C     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
C
C     UPMXP,UPMX=maximum H2PO4 uptake in band unlimited,limited by O2
C     RTARP=root surface area per plant from grosub.f
C     FWSRT=protein concentration relative to 5%
C     TFN4=temperature function for root growth
C     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
C     WFR=constraint by O2 consumption on all biological processes
C
      UPMXP=UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLPOB(L,NY,NX)*AMIN1(FCUP,FPUP)
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF PO4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFH2B=soil-root convective H2PO4 flux per plant in band
C     DIFH2B=soil-root H2PO4 diffusion per plant in band
C     CH2P4B=H2PO4 concentration in band
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
C     RTKH2B,RTKHPB=H2PO4 uptake per plant in band lmtd,unlmtd by O2
C     H2PXM,H2PXB=minimum,maximum H2PO4 available for uptake in band
C     FPOBX=fraction of total H2PO4 uptake in band by root,myco populn
C     RUPP2B,RUPH2B=H2PO4 uptake in band unlimited,limited by H2PO4
C     RUOH2B=H2PO4 uptake in band unlimited by O2
C     RUCH2B=H2PO4 uptake in band unlimited by nonstructural C
C
      X=(DIFH2B+RMFH2B)*CH2P4B(L,NY,NX)
      Y=DIFH2B*UPMNPO(N,NZ,NY,NX)
      B=-UPMX-DIFH2B*UPKMPO(N,NZ,NY,NX)-X+Y
      C=(X-Y)*UPMX
      RTKH2B=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFH2B*UPKMPO(N,NZ,NY,NX)-X+Y
      CP=(X-Y)*UPMXP
      RTKHPB=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      H2PXM=UPMNPO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
      H2PXB=AMAX1(0.0,FPOBX*(H2POB(L,NY,NX)-H2PXM))
      RUPP2B(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH2B*PP(NZ,NY,NX))
      RUPH2B(N,L,NZ,NY,NX)=AMIN1(H2PXB,RUPP2B(N,L,NZ,NY,NX))
      RUOH2B(N,L,NZ,NY,NX)=AMIN1(H2PXB
     2,AMAX1(0.0,RTKHPB*PP(NZ,NY,NX)))
      RUCH2B(N,L,NZ,NY,NX)=RUPH2B(N,L,NZ,NY,NX)/FCUP
      ELSE
      RUPP2B(N,L,NZ,NY,NX)=0.0
      RUPH2B(N,L,NZ,NY,NX)=0.0
      RUOH2B(N,L,NZ,NY,NX)=0.0
      RUCH2B(N,L,NZ,NY,NX)=0.0
      ENDIF
C
C     HPO4 UPTAKE IN NON-BAND SOIL ZONE
C
C     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
C     CH1P4=HPO4 concentration in non-band
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
C     UPWTRP=root water uptake per plant
C     RMFH1P=soil-root convective HPO4 flux per plant in non-band
C     DIFH1P=soil-root HPO4 diffusion per plant in non-band
C
      IF(VLPO4(L,NY,NX).GT.ZERO.AND.CH1P4(L,NY,NX)
     2.GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH1P=UPWTRP*CH1P4(L,NY,NX)*VLPO4(L,NY,NX)
      DIFH1P=DIFFL*VLPO4(L,NY,NX)
C
C     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
C     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
C
C     UPMXP,UPMX=max HPO4 uptake in non-band unlimited,limited by O2
C     RTARP=root surface area per plant from grosub.f
C     FWSRT=protein concentration relative to 5%
C     TFN4=temperature function for root growth
C     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
C     WFR=constraint by O2 consumption on all biological processes
C
      UPMXP=0.1*UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLPO4(L,NY,NX)*AMIN1(FCUP,FPUP)
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFH1P=soil-root convective HPO4 flux per plant in non-band
C     DIFH1P=soil-root HPO4 diffusion per plant in non-band
C     CH1P4=HPO4 concentration in non-band
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
C     RTKH1P,RTKHP1=HPO4 uptake per plant in non-band lmtd,unlmtd by O2
C     H1POM,H1POX=minimum,maximum HPO4 available for uptake in non-band
C     FP14X=fraction of total HPO4 uptake in non-band by root,myco populn
C     RUPP1P,RUPH1P=HPO4 uptake in non-band unlimited,limited by HPO4
C     RUOH1P=HPO4 uptake in non-band unlimited by O2
C     RUCH1P=HPO4 uptake in non-band unlimited by nonstructural C
C
      X=(DIFH1P+RMFH1P)*CH1P4(L,NY,NX)
      Y=DIFH1P*UPMNPO(N,NZ,NY,NX)
      B=-UPMX-DIFH1P*UPKMPO(N,NZ,NY,NX)-X+Y
      C=(X-Y)*UPMX
      RTKH1P=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFH1P*UPKMPO(N,NZ,NY,NX)-X+Y
      CP=(X-Y)*UPMXP
      RTKHP1=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      H1POM=UPMNPO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLPO4(L,NY,NX)
      H1POX=AMAX1(0.0,FP14X*(H1PO4(L,NY,NX)-H1POM))
      RUPP1P(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH1P*PP(NZ,NY,NX))
      RUPH1P(N,L,NZ,NY,NX)=AMIN1(H1POX,RUPP1P(N,L,NZ,NY,NX))
      RUOH1P(N,L,NZ,NY,NX)=AMIN1(H1POX,AMAX1(0.0
     2,RTKHP1*PP(NZ,NY,NX)))
      RUCH1P(N,L,NZ,NY,NX)=RUPH1P(N,L,NZ,NY,NX)/FCUP
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.3)THEN
C     WRITE(*,2226)'UPPO4',I,J,NZ,L,N,RUPH2P(N,L,NZ,NY,NX),FPO4X
C    2,H2PO4(L,NY,NX),RUPP2P(N,L,NZ,NY,NX),UPMX,DIFPO,UPKMPO(N,NZ,NY,NX)
C    3,UPMNPO(N,NZ,NY,NX),RMFH2P,CH2P4(L,NY,NX),UPMXP,WFR(N,L,NZ,NY,NX)
C    4,FCUP,FZUP,FPUP,UPMXPO(N,NZ,NY,NX),RTARP(N,L,NZ,NY,NX),FWSRT
C    5,TFN4(L,NZ,NY,NX),DIFFL,FH2P,CPO4S(L,NY,NX),CPOOLR(N,L,NZ,NY,NX)
C    6,PPOOLR(N,L,NZ,NY,NX),RTKH2P,PP(NZ,NY,NX)
C    2,RTLGP(N,L,NZ,NY,NX)
2226  FORMAT(A8,5I4,40E12.4)
C     ENDIF
      ELSE
      RUPP1P(N,L,NZ,NY,NX)=0.0
      RUPH1P(N,L,NZ,NY,NX)=0.0
      RUOH1P(N,L,NZ,NY,NX)=0.0
      RUCH1P(N,L,NZ,NY,NX)=0.0
      ENDIF
C
C     HPO4 UPTAKE IN BAND SOIL ZONE
C
C     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
C     CH1P4B=HPO4 concentration in band
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
C     UPWTRP=root water uptake per plant
C     RMFH1B=soil-root convective HPO4 flux per plant in band
C     DIFH1B=soil-root HPO4 diffusion per plant in band
C
      IF(VLPOB(L,NY,NX).GT.ZERO.AND.CH1P4B(L,NY,NX)
     2.GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH2B=UPWTRP*CH1P4B(L,NY,NX)*VLPOB(L,NY,NX)
      DIFH1B=DIFFL*VLPOB(L,NY,NX)
C
C     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
C     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
C
C     UPMXP,UPMX=maximum HPO4 uptake in band unlimited,limited by O2
C     RTARP=root surface area per plant from grosub.f
C     FWSRT=protein concentration relative to 5%
C     TFN4=temperature function for root growth
C     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
C     WFR=constraint by O2 consumption on all biological processes
C
      UPMXP=0.1*UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLPOB(L,NY,NX)*AMIN1(FCUP,FPUP)
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX)
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFH1B=soil-root convective HPO4 flux per plant in band
C     DIFH1B=soil-root HPO4 diffusion per plant in band
C     CH1P4B=HPO4 concentration in band
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
C     RTKH1B,RTKHB1=HPO4 uptake per plant in band lmtd,unlmtd by O2
C     H1PXM,H1PXB=minimum,maximum HPO4 available for uptake in band
C     FP1BX=fraction of total HPO4 uptake in band by root,myco populn
C     RUPP1B,RUPH1B=HPO4 uptake in band unlimited,limited by H2PO4
C     RUOH1B=HPO4 uptake in band unlimited by O2
C     RUCH1B=HPO4 uptake in band unlimited by nonstructural C
C
      X=(DIFH1B+RMFH2B)*CH1P4B(L,NY,NX)
      Y=DIFH1B*UPMNPO(N,NZ,NY,NX)
      B=-UPMX-DIFH1B*UPKMPO(N,NZ,NY,NX)-X+Y
      C=(X-Y)*UPMX
      RTKH1B=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFH1B*UPKMPO(N,NZ,NY,NX)-X+Y
      CP=(X-Y)*UPMXP
      RTKHB1=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      H1PXM=UPMNPO(N,NZ,NY,NX)*VOLW(L,NY,NX)*VLPOB(L,NY,NX)
      H1PXB=AMAX1(0.0,FP1BX*(H1POB(L,NY,NX)-H1PXM))
      RUPP1B(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH1B*PP(NZ,NY,NX))
      RUPH1B(N,L,NZ,NY,NX)=AMIN1(H1PXB,RUPP1B(N,L,NZ,NY,NX))
      RUOH1B(N,L,NZ,NY,NX)=AMIN1(H1PXB
     2,AMAX1(0.0,RTKHB1*PP(NZ,NY,NX)))
      RUCH1B(N,L,NZ,NY,NX)=RUPH1B(N,L,NZ,NY,NX)/FCUP
      ELSE
      RUPP1B(N,L,NZ,NY,NX)=0.0
      RUPH1B(N,L,NZ,NY,NX)=0.0
      RUOH1B(N,L,NZ,NY,NX)=0.0
      RUCH1B(N,L,NZ,NY,NX)=0.0
      ENDIF
      ELSE
      RUPP2P(N,L,NZ,NY,NX)=0.0
      RUPH2P(N,L,NZ,NY,NX)=0.0
      RUOH2P(N,L,NZ,NY,NX)=0.0
      RUCH2P(N,L,NZ,NY,NX)=0.0
      RUPP2B(N,L,NZ,NY,NX)=0.0
      RUPH2B(N,L,NZ,NY,NX)=0.0
      RUOH2B(N,L,NZ,NY,NX)=0.0
      RUCH2B(N,L,NZ,NY,NX)=0.0
      RUPP1P(N,L,NZ,NY,NX)=0.0
      RUPH1P(N,L,NZ,NY,NX)=0.0
      RUOH1P(N,L,NZ,NY,NX)=0.0
      RUCH1P(N,L,NZ,NY,NX)=0.0
      RUPP1B(N,L,NZ,NY,NX)=0.0
      RUPH1B(N,L,NZ,NY,NX)=0.0
      RUOH1B(N,L,NZ,NY,NX)=0.0
      RUCH1B(N,L,NZ,NY,NX)=0.0
      ENDIF
      ELSE
      RUNNHP(N,L,NZ,NY,NX)=0.0
      RUPNH4(N,L,NZ,NY,NX)=0.0
      RUONH4(N,L,NZ,NY,NX)=0.0
      RUCNH4(N,L,NZ,NY,NX)=0.0
      RUNNBP(N,L,NZ,NY,NX)=0.0
      RUPNHB(N,L,NZ,NY,NX)=0.0
      RUONHB(N,L,NZ,NY,NX)=0.0
      RUCNHB(N,L,NZ,NY,NX)=0.0
      RUNNOP(N,L,NZ,NY,NX)=0.0
      RUPNO3(N,L,NZ,NY,NX)=0.0
      RUONO3(N,L,NZ,NY,NX)=0.0
      RUCNO3(N,L,NZ,NY,NX)=0.0
      RUNNXP(N,L,NZ,NY,NX)=0.0
      RUPNOB(N,L,NZ,NY,NX)=0.0
      RUONOB(N,L,NZ,NY,NX)=0.0
      RUCNOB(N,L,NZ,NY,NX)=0.0
      RUPP2P(N,L,NZ,NY,NX)=0.0
      RUPH2P(N,L,NZ,NY,NX)=0.0
      RUOH2P(N,L,NZ,NY,NX)=0.0
      RUCH2P(N,L,NZ,NY,NX)=0.0
      RUPP2B(N,L,NZ,NY,NX)=0.0
      RUPH2B(N,L,NZ,NY,NX)=0.0
      RUOH2B(N,L,NZ,NY,NX)=0.0
      RUCH2B(N,L,NZ,NY,NX)=0.0
      RUPP1P(N,L,NZ,NY,NX)=0.0
      RUPH1P(N,L,NZ,NY,NX)=0.0
      RUOH1P(N,L,NZ,NY,NX)=0.0
      RUCH1P(N,L,NZ,NY,NX)=0.0
      RUPP1B(N,L,NZ,NY,NX)=0.0
      RUPH1B(N,L,NZ,NY,NX)=0.0
      RUOH1B(N,L,NZ,NY,NX)=0.0
      RUCH1B(N,L,NZ,NY,NX)=0.0
      ENDIF
      ELSE
      RCOFLA(N,L,NZ,NY,NX)=0.0
      ROXFLA(N,L,NZ,NY,NX)=0.0
      RCHFLA(N,L,NZ,NY,NX)=0.0
      RN2FLA(N,L,NZ,NY,NX)=0.0
      RNHFLA(N,L,NZ,NY,NX)=0.0
      RCODFA(N,L,NZ,NY,NX)=0.0
      ROXDFA(N,L,NZ,NY,NX)=0.0
      RCHDFA(N,L,NZ,NY,NX)=0.0
      RN2DFA(N,L,NZ,NY,NX)=0.0
      RNHDFA(N,L,NZ,NY,NX)=0.0
      RCO2S(N,L,NZ,NY,NX)=0.0
      RUPOXS(N,L,NZ,NY,NX)=0.0
      RUPCHS(N,L,NZ,NY,NX)=0.0
      RUPN2S(N,L,NZ,NY,NX)=0.0
      RUPN3S(N,L,NZ,NY,NX)=0.0
      RCO2P(N,L,NZ,NY,NX)=0.0
      RUPOXP(N,L,NZ,NY,NX)=0.0
      DO 395 K=0,4
      RDFOMC(N,K,L,NZ,NY,NX)=0.0
      RDFOMN(N,K,L,NZ,NY,NX)=0.0
      RDFOMP(N,K,L,NZ,NY,NX)=0.0
395   CONTINUE
      WFR(N,L,NZ,NY,NX)=1.0
      RUNNHP(N,L,NZ,NY,NX)=0.0
      RUPNH4(N,L,NZ,NY,NX)=0.0
      RUONH4(N,L,NZ,NY,NX)=0.0
      RUCNH4(N,L,NZ,NY,NX)=0.0
      RUNNBP(N,L,NZ,NY,NX)=0.0
      RUPNHB(N,L,NZ,NY,NX)=0.0
      RUONHB(N,L,NZ,NY,NX)=0.0
      RUCNHB(N,L,NZ,NY,NX)=0.0
      RUNNOP(N,L,NZ,NY,NX)=0.0
      RUPNO3(N,L,NZ,NY,NX)=0.0
      RUONO3(N,L,NZ,NY,NX)=0.0
      RUCNO3(N,L,NZ,NY,NX)=0.0
      RUNNXP(N,L,NZ,NY,NX)=0.0
      RUPNOB(N,L,NZ,NY,NX)=0.0
      RUONOB(N,L,NZ,NY,NX)=0.0
      RUCNOB(N,L,NZ,NY,NX)=0.0
      RUPP2P(N,L,NZ,NY,NX)=0.0
      RUPH2P(N,L,NZ,NY,NX)=0.0
      RUOH2P(N,L,NZ,NY,NX)=0.0
      RUCH2P(N,L,NZ,NY,NX)=0.0
      RUPP2B(N,L,NZ,NY,NX)=0.0
      RUPH2B(N,L,NZ,NY,NX)=0.0
      RUOH2B(N,L,NZ,NY,NX)=0.0
      RUCH2B(N,L,NZ,NY,NX)=0.0
      RUPP1P(N,L,NZ,NY,NX)=0.0
      RUPH1P(N,L,NZ,NY,NX)=0.0
      RUOH1P(N,L,NZ,NY,NX)=0.0
      RUCH1P(N,L,NZ,NY,NX)=0.0
      RUPP1B(N,L,NZ,NY,NX)=0.0
      RUPH1B(N,L,NZ,NY,NX)=0.0
      RUOH1B(N,L,NZ,NY,NX)=0.0
      RUCH1B(N,L,NZ,NY,NX)=0.0
      IF(N.EQ.1)RUPNF(L,NZ,NY,NX)=0.0
      ENDIF
C
C     TOTAL C,N,P EXCHANGE BETWEEN ROOTS AND SOIL
C
C     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
C     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
C     XOQCS,XOQNZ,XOQPS=accumulated change in DOC,DON,DOP from nitro.f
C     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
C     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
C     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
C
      DO 295 K=0,4
      UPOMC(NZ,NY,NX)=UPOMC(NZ,NY,NX)+RDFOMC(N,K,L,NZ,NY,NX)
      UPOMN(NZ,NY,NX)=UPOMN(NZ,NY,NX)+RDFOMN(N,K,L,NZ,NY,NX)
      UPOMP(NZ,NY,NX)=UPOMP(NZ,NY,NX)+RDFOMP(N,K,L,NZ,NY,NX)
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-RDFOMC(N,K,L,NZ,NY,NX)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-RDFOMN(N,K,L,NZ,NY,NX)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-RDFOMP(N,K,L,NZ,NY,NX)
C     WRITE(*,8766)'XOQCSU',I,J,NX,NY,L,NZ,K,N
C    2,XOQCS(K,L,NY,NX),RDFOMC(N,K,L,NZ,NY,NX)
8766  FORMAT(A8,8I4,12E12.4)
295   CONTINUE
      UPNH4(NZ,NY,NX)=UPNH4(NZ,NY,NX)+RUPNH4(N,L,NZ,NY,NX)
     2+RUPNHB(N,L,NZ,NY,NX)
      UPNO3(NZ,NY,NX)=UPNO3(NZ,NY,NX)+RUPNO3(N,L,NZ,NY,NX)
     2+RUPNOB(N,L,NZ,NY,NX)
      UPH2P(NZ,NY,NX)=UPH2P(NZ,NY,NX)+RUPH2P(N,L,NZ,NY,NX)
     2+RUPH2B(N,L,NZ,NY,NX)
      UPH1P(NZ,NY,NX)=UPH1P(NZ,NY,NX)+RUPH1P(N,L,NZ,NY,NX)
     2+RUPH1B(N,L,NZ,NY,NX)
C     IF(J.EQ.12)THEN
C     WRITE(*,8765)'PLANT',I,J,NX,NY,L,NZ,N,TFOXYX,TFNH4X
C    2,TFNO3X,TFPO4X,TFNHBX,TFNOBX,TFPOBX
8765  FORMAT(A8,7I4,7F15.6)
C     ENDIF
950   CONTINUE
955   CONTINUE
      TLEC(NY,NX)=TLEC(NY,NX)+EFLXC(NZ,NY,NX)*RA(NZ,NY,NX)
      TSHC(NY,NX)=TSHC(NY,NX)+SFLXC(NZ,NY,NX)*RA(NZ,NY,NX)
      IF(OSTRD.GT.ZEROP(NZ,NY,NX))THEN
      OSTR(NZ,NY,NX)=OSTRN/OSTRD
      ELSE
      OSTR(NZ,NY,NX)=0.0
      ENDIF
      ENDIF
9985  CONTINUE
9990  CONTINUE
9995  CONTINUE
      RETURN
      END
