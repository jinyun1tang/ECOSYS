
      SUBROUTINE uptake(I,J,NFZ,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE CALCULATES EXCHANGES OF ENERGY, C, N AND P
C     BETWEEN THE CANOPY AND THE ATMOSPHERE AND BETWEEN ROOTS AND 
C     THE SOIL
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
      include "blk19a.h"
      include "blk19b.h"
      include "blk1u.h"
      DIMENSION PSIST1(JZ),PATH(2,JZ),RRADL(2,JZ)
     2,RSRT(2,JZ),ILYR(2,JZ),RSRG(2,JZ),RSR1(2,JZ),RSR2(2,JZ) 
     3,RSSX(2,JZ),RSRS(2,JZ),WTRTG(JZ),FPQ(2,JZ,05),FPP(2,JZ,05) 
     4,FRTDPX(JZ,05),RTARR(2,JZ),VOLPU(JZ),VOLWU(JZ),UPWTRM(2,JZ)
     5,DTHRMT(JY,JX),DSHQT(JY,JX),DVPQT(JY,JX)
     5,DCO2QT(JY,JX),DCH4QT(JY,JX),DOXYQT(JY,JX)
     6,DTHRM(2,JV,JH),DSHQ(2,JV,JH),DVPQ(2,JV,JH)  
     6,DCO2Q(2,JV,JH),DCH4Q(2,JV,JH),DOXYQ(2,JV,JH)  
C
C     MXN=max number of cycles in convergence solution 
C        for water uptake
C     DIFFX=acceptance criteria for DIFFU-DIFFZ in convergence
C        solution for water uptake 
C     DIFFT=acceptance criteria for temperature change in convergence
C        solution for water uptake 
C     VHCPYM,VHCPXM=minimum heat capacity for solving canopy 
C        water potential,aerodymamic energy exchange (MJ m-2 K-1)
C     DTKX1,DTKX2,DTKD1,DTKD2=step change in canopy(X),standing
C        dead(D) surface temperature during convergence for energy
C        balance when canopy heat capacity > (1) or < (2) VHCPXM (K)         
C     SNH3X=NH3 solubility at 25 oC (g m-3 water/(g m-3 air))
C     EMMC=canopy emissivity
C     EMODW=wood modulus of elasticity (MPa)
C     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition on root,mycorrhizal 
C        N,P uptake (g g-1)
C     FMN=min PFT:total population ratio
C     FEXU=rate constant for root C,N,P exudation (h-1)
C     RAH,RAE=canopy isothermal surface resistance to sensible,latent
C        heat (h m-1)
C
      PARAMETER(MXN=200,DIFFX=1.0E-06,DIFFT=1.0E-03)
      PARAMETER(VHCPYM=0.419E-04,VHCPXM=0.838E-03
     2,DTKX1=0.125,DTKX2=0.025,DTKD1=0.125,DTKD2=0.025)
      PARAMETER(SNH3X=2.852E+02,EMMC=0.97,EMODW=50.0)
      PARAMETER(ZCKI=0.1,PCKI=0.1,ZPKI=1.0,PZKI=1.0)
      PARAMETER(FMN=1.0E-04,FEXU=0.75E-03)
      PARAMETER(RAH=0.278E-02,RAE=0.139E-01)
      PARAMETER(CZALKI=1.0E-04,CZFEKI=1.0E-04,CZCAKI=1.0E-04
     2,CZMGKI=1.0E-04,CZNAKI=1.0E-04,CZKAKI=1.0E-04
     2,CZSOKI=1.0E-04,CZCLKI=1.0E-04) 
      REAL*4 TKGO,TKSO
C
C     GAS, ENERGY TRANSFER BETWEEN GRID CELLS (EG CANOPY FIRE SPREAD)
C
      DO 9895 NX=NHW,NHE
      DO 9890 NY=NVN,NVS
      DTHRMT(NY,NX)=0.0
      DSHQT(NY,NX)=0.0
      DVPQT(NY,NX)=0.0
      DCO2QT(NY,NX)=0.0
      DCH4QT(NY,NX)=0.0
      DOXYQT(NY,NX)=0.0
C
C     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
C
C     N2,N1=NY,NX of source grid cell
C     N5,N4=NY,NX of destination grid cell
C
      N1=NX
      N2=NY
      DO 9885 N=1,2
      IF(N.EQ.1)THEN
      IF(NX.EQ.NHE)THEN
      GO TO 9885
      ELSE
      N4=NX+1
      N5=NY
      ENDIF
      ELSEIF(N.EQ.2)THEN
      IF(NY.EQ.NVS)THEN
      GO TO 9885
      ELSE
      N4=NX
      N5=NY+1
      ENDIF
      ENDIF
C
C     EXCHANGE OF GAS, THERMAL AND SENSIBLE HEAT FLUXES 
C     BETWEEN ADJACENT CANOPIES
C
C     DTHRM=net lateral LW flux (MJ t-1)
C     TKCT=bulk canopy surface temperature (K)
C     EMMC=canopy emissivity
C     FRADT=fraction of radiation received by all PFT canopies
C        from ‘hour1.f’
C     TKY=canopy air temperature (K)
C     VHCPQ=canopy heat capacity used to calculate 
C        canopy air temperature from ‘hour1.f’ 
C        (MJ K-1,m3)
C     DSHQ=net lateral sensible heat flux (MJ t-1)
C     TKQT=bulk canopy air temperature (K)
C     RACG=canopy aerodynamic resistance (h m-1)
C     DVPQ=net lateral vapor flux (m3 t-1)
C     VPQT=bulk canopy vapor concentration (m3 m-3) 
C     DCO2Q=net lateral CO2 flux (g C t-1)
C     CO2Q=canopy CO2 concentration (g C m-3)
C     ZT=bulk canopy height (m)
C     DCH4Q=net lateral CH4 flux (g C t-1)
C     CH4Q=canopy CH4 concentration (g C m-3)
C     DOXYQ=net lateral O2 flux (g O t-1)
C     COXYQ=canopy O2 concentration (g O m-3)
C
      IF(TKQT(N2,N1).GT.ZERO.AND.TKQT(N5,N4).GT.ZERO)THEN
      DTHRM(N,N5,N4)=(TKCT(N2,N1)**4*AREA(3,NU(N2,N1),N2,N1)
     2-TKCT(N5,N4)**4*AREA(3,NU(N5,N4),N5,N4))
     3*EMMC*2.04E-10*AMIN1(FRADT(N2,N1),FRADT(N5,N4))
      RACGQ=5.0*(AMAX1(RABM,AMIN1(RABZ
     2,0.5*(RACG(NPH,N2,N1)+RACG(NPH,N5,N4)))))
      TKY=(TKQT(N2,N1)*VHCPQ(N2,N1)+TKQT(N5,N4) 
     2*VHCPQ(N5,N4))/(VHCPQ(N2,N1)+VHCPQ(N5,N4))
      IF(TKQT(N2,N1).GT.TKQT(N5,N4))THEN
      DSHQ(N,N5,N4)=1.25E-03*(TKQT(N2,N1)*AREA(3,NU(N2,N1),N2,N1)
     2-AMAX1(TKQT(N5,N4),TKY)*AREA(3,NU(N5,N4),N5,N4))
     3/RACGQ*AMIN1(FRADT(N2,N1),FRADT(N5,N4))
      ELSE
      DSHQ(N,N5,N4)=1.25E-03*(TKQT(N2,N1)*AREA(3,NU(N2,N1),N2,N1)
     2-AMIN1(TKQT(N5,N4),TKY)*AREA(3,NU(N5,N4),N5,N4))
     3/RACGQ*AMIN1(FRADT(N2,N1),FRADT(N5,N4))
      ENDIF
      DVPQ(N,N5,N4)=(VPQT(N2,N1)*AREA(3,NU(N2,N1),N2,N1)
     2-VPQT(N5,N4)*AREA(3,NU(N5,N4),N5,N4))
     3/RACGQ*AMIN1(FRADT(N2,N1),FRADT(N5,N4))
      DCO2Q(N,N5,N4)=(CO2Q(N2,N1)-CO2Q(N5,N4))/RACGQ
     3*5.357E-04*273.15/AMAX1(TKQT(N2,N1),TKQT(N5,N4))
     3*AMIN1(ZT(N2,N1)*AREA(3,NU(N2,N1),N2,N1)
     3,ZT(N5,N4)*AREA(3,NU(N5,N4),N5,N4))
      DCH4Q(N,N5,N4)=(CH4Q(N2,N1)-CH4Q(N5,N4))/RACGQ 
     3*5.357E-04*273.15/AMAX1(TKQT(N2,N1),TKQT(N5,N4))
     3*AMIN1(ZT(N2,N1)*AREA(3,NU(N2,N1),N2,N1)
     3,ZT(N5,N4)*AREA(3,NU(N5,N4),N5,N4))
      DOXYQ(N,N5,N4)=(OXYQ(N2,N1)-OXYQ(N5,N4))/RACGQ 
     3*1.429E-03*273.15/AMAX1(TKQT(N2,N1),TKQT(N5,N4))
     3*AMIN1(ZT(N2,N1)*AREA(3,NU(N2,N1),N2,N1)
     3,ZT(N5,N4)*AREA(3,NU(N5,N4),N5,N4))
      ELSE
      DTHRM(N,N5,N4)=0.0
      DSHQ(N,N5,N4)=0.0
      DVPQ(N,N5,N4)=0.0
      DCO2Q(N,N5,N4)=0.0
      DCH4Q(N,N5,N4)=0.0
      DOXYQ(N,N5,N4)=0.0
      ENDIF
C     IF(ICHKF.EQ.1)THEN
C     WRITE(*,4410)'DTHRM',I,J,NFZ,N1,N2,N4,N5,N 
C    2,DTHRM(N,N5,N4),TKCT(N2,N1),TKCT(N5,N4)
C    3,DSHQ(N,N5,N4),TKQT(N2,N1),TKQT(N5,N4),TKY
C    4,DVPQ(N,N5,N4),VPQT(N2,N1),VPQT(N5,N4)
C    5,DCO2Q(N,N5,N4),CO2Q(N2,N1),CO2Q(N5,N4)
C    6,DCH4Q(N,N5,N4),CH4Q(N2,N1),CH4Q(N5,N4)
C    7,DOXYQ(N,N5,N4),OXYQ(N2,N1),OXYQ(N5,N4)
C    8,RACG(NPH,N2,N1),RACG(NPH,N5,N4),RACGQ,RABM,RABZ
C    9,FRADT(N2,N1),FRADT(N5,N4),VHCPQ(N2,N1),VHCPQ(N5,N4)
C    3,AMIN1(ZT(N2,N1)*AREA(3,NU(N2,N1),N2,N1)
C    3,ZT(N5,N4)*AREA(3,NU(N5,N4),N5,N4))
4410  FORMAT(A8,8I4,30E12.4)
C     ENDIF 
9885  CONTINUE
9890  CONTINUE
9895  CONTINUE
C
C     LOCATE EXTERNAL BOUNDARIES
C
      DO 9795 NX=NHW,NHE
      DO 9790 NY=NVN,NVS
C
C     N2,N1=NY,NX of source grid cell
C     M5,M4=NY,NX of destination grid cell
C     N5,N4=NY,NX of destination grid cell 
C
      N1=NX
      N2=NY
      DO 9785 N=1,2
      DO 9780 NN=1,2
      IF(N.EQ.1)THEN
      N4=NX+1
      N5=NY
      IF(NN.EQ.1)THEN
      IF(NX.EQ.NHE)THEN
      M1=NX
      M2=NY
      M4=NX+1
      M5=NY
      XN=-1.0
      ELSE
      GO TO 9780
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NX.EQ.NHW)THEN
      M1=NX+1
      M2=NY
      M4=NX
      M5=NY
      XN=1.0
      ELSE
      GO TO 9780
      ENDIF
      ENDIF
      ELSEIF(N.EQ.2)THEN
      N4=NX
      N5=NY+1
      IF(NN.EQ.1)THEN
      IF(NY.EQ.NVS)THEN
      M1=NX
      M2=NY
      M4=NX
      M5=NY+1
      XN=-1.0
      ELSE
      GO TO 9780
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NY.EQ.NVN)THEN
      M1=NX
      M2=NY+1
      M4=NX
      M5=NY
      XN=1.0
      ELSE
      GO TO 9780
      ENDIF
      ENDIF
      ENDIF
C
C     BOUNDARY EXCHANGE OF MASS, THERMAL AND SENSIBLE HEAT FLUXES 
C     IN EACH CANOPY (SET TO 0)
C
      DTHRM(N,M5,M4)=0.0
      DSHQ(N,M5,M4)=0.0
      DVPQ(N,M5,M4)=0.0
      DCO2Q(N,M5,M4)=0.0
      DCH4Q(N,M5,M4)=0.0
      DOXYQ(N,M5,M4)=0.0
C
C     NET EXCHANGE OF MASS, THERMAL AND SENSIBLE HEAT FLUXES 
C     IN EACH CANOPY 
C
C     DTHRMT=total lateral LW radiation (MJ t-1)
C     DSHQT=total lateral sensible heat flux (MJ t-1)
C     DVPQT=total lateral vapor flux (m3 t-1)
C     DCO2QT=total lateral CO2 flux g C m-3)
C     DCH4QT=total lateral CH4 flux (g C m-3)
C     DOXYQT=total lateral O2 flux (g O t-1)
C
      DTHRMT(N2,N1)=DTHRMT(N2,N1)+DTHRM(N,N2,N1)
     2-DTHRM(N,N5,N4)
      DSHQT(N2,N1)=DSHQT(N2,N1)+DSHQ(N,N2,N1)
     2-DSHQ(N,N5,N4)
      DVPQT(N2,N1)=DVPQT(N2,N1)+DVPQ(N,N2,N1)
     2-DVPQ(N,N5,N4)
      DCO2QT(N2,N1)=DCO2QT(N2,N1)+DCO2Q(N,N2,N1)
     2-DCO2Q(N,N5,N4)
      DCH4QT(N2,N1)=DCH4QT(N2,N1)+DCH4Q(N,N2,N1)
     2-DCH4Q(N,N5,N4)
      DOXYQT(N2,N1)=DOXYQT(N2,N1)+DOXYQ(N,N2,N1)
     2-DOXYQ(N,N5,N4)
C     WRITE(*,4412)'DTHRMT',I,J,NFZ,N1,N2,N4,N5,M4,M5,N
C    2,DTHRMT(N2,N1),DTHRM(N,N2,N1),DTHRM(N,N5,N4),DTHRM(N,M5,M4)
C    2,DSHQT(N2,N1),DSHQ(N,N2,N1),DSHQ(N,N5,N4),DSHQ(N,M5,M4)
C    2,DVPQT(N2,N1),DVPQ(N,N2,N1),DVPQ(N,N5,N4),DVPQ(N,M5,M4)
C    2,DCO2QT(N2,N1),DCO2Q(N,N2,N1),DCO2Q(N,N5,N4),DCO2Q(N,M5,M4)
C    2,DCH4QT(N2,N1),DCH4Q(N,N2,N1),DCH4Q(N,N5,N4),DCH4Q(N,M5,M4)
C    2,DOXYQT(N2,N1),DOXYQ(N,N2,N1),DOXYQ(N,N5,N4),DOXYQ(N,M5,M4)
4412  FORMAT(A8,10I4,12E12.4) 
9780  CONTINUE
9785  CONTINUE
9790  CONTINUE
9795  CONTINUE
C
C     NET CANOPY GAS EXCHANGE WITHIN GRID CELL 
C
C     XCNET,XHNET,XONET=total CO2,CH4,O2 exchange by all PFT canopies
C        (g C,C,O t-1)
C     DCO2QT=total lateral CO2 flux (g C t-1)
C     DCH4QT=total lateral CH4 flux (g C t-1)
C     DOXYQT=total lateral O2 flux (g O t-1)
C     XNFH=time step for biological activity from ‘wthr.f’ (h t-1)
C
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      XCNET(NY,NX)=XCNET(NY,NX)-DCO2QT(NY,NX)*XNFH
      XHNET(NY,NX)=XHNET(NY,NX)-DCH4QT(NY,NX)*XNFH
      XONET(NY,NX)=XONET(NY,NX)-DOXYQT(NY,NX)*XNFH
C     WRITE(*,4409)'XONET',I,J,NFZ,NX,NY
C    2,XCNET(NY,NX),DCO2QT(NY,NX)*XNFH
C    2,XHNET(NY,NX),DCH4QT(NY,NX)*XNFH
C    2,XONET(NY,NX),DOXYQT(NY,NX)*XNFH
C    3,OXYC(NY,NX),TKQ(NY,NX)
4409  FORMAT(A8,5I4,20E12.4)
C
C     RESET TOTAL UPTAKE ARRAYS
C
      DO 9984 NZ=1,NP0(NY,NX)
      PSILT(NZ,NY,NX)=0.95*PSILT(NZ,NY,NX)
      RFLXC(NZ,NY,NX)=0.0
      EFLXC(NZ,NY,NX)=0.0
      SFLXC(NZ,NY,NX)=0.0
      HFLXC(NZ,NY,NX)=0.0
      VFLXC(NZ,NY,NX)=0.0
      THRMP(NZ,NY,NX)=0.0
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
      IF(ISALTG.NE.0)THEN
      RUPZAL(N,L,NZ,NY,NX)=0.0
      RUPZFE(N,L,NZ,NY,NX)=0.0
      RUPZCA(N,L,NZ,NY,NX)=0.0
      RUPZMG(N,L,NZ,NY,NX)=0.0
      RUPZNA(N,L,NZ,NY,NX)=0.0
      RUPZKA(N,L,NZ,NY,NX)=0.0
      RUPZSO(N,L,NZ,NY,NX)=0.0
      RUPZCL(N,L,NZ,NY,NX)=0.0
      ENDIF
9984  CONTINUE
C
C     PSIST1=total soil water potential PSIST adjusted for 
C        surface elevation (MPa)
C     ALT=surface elevation (m)
C     VOLWU,VOLWM=water volume available for root uptake,
C        total water volume (m3)
C     THETY,VOLWY=water concentration,content at hygroscopic (PSIHY)
C        water potential set in ‘starts.f’ (m3 m-3)
C     VOLPU=air volume (m3)
C     WTRTG=total biome root mass (g C) 
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
C     EQUILIBRATE CANOPY ET WITH ROOT WATER UPTAKE
C
      DO 9985 NZ=1,NP(NY,NX)
      OSTRN=0.0
      OSTRD=0.0
C
C     IF PLANT SPECIES EXISTS
C
C     IFLGC=PFT flag:0=not active,1=active
C     PP=PFT population (n) 
C
      IF(IFLGC(NZ,NY,NX).EQ.1.AND.PP(NZ,NY,NX).GT.0.0)THEN
C
C     APPLY CLUMPING FACTOR TO LEAF SURFACE AREA DEFINED BY
C     INCLINATION N, LAYER L, NODE K, BRANCH NB, SPECIES NZ, 
C     N-S POSITION NY, E-W POSITION NX(AZIMUTH M ASSUMED UNIFORM)
C
      DO 500 NB=1,NBR(NZ,NY,NX)
      DO 550 K=1,25
C
C     ARLF=node leaf area (m2 m-2)
C     WSLF=leaf protein content (g m-2)
C     KLEAFX=lowest leafed node number
C     SURFX,SURF=unself-shaded,total leaf surface area used in
C        ‘hour1.f’,stomaTe.f’,grosub.f’ (m2 m-2)
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
C     INITIALIZE LIVING CANOPY TEMPERATURE AND VAPOR CONCENTRATION, 
C     AND CALL ‘STOMATE.F’ TO CALCULATE MINIMUM CANOPY STOMATAL 
C     RESISTANCE FOR SUBSEQUENT USE IN ENERGY EXCHANGE CALCULATIONS
C
C     TKQC,VPQC=canopy air temperature, vapor concentration (K,m3 m-3)
C     STOMATE.F:solve for minimum canopy stomatal resistance
C
      TKQY=TKQC(NZ,NY,NX)
      VPQY=VPQC(NZ,NY,NX)
      CALL STOMATE(I,J,NFZ,NZ,NY,NX)
C
C     CALCULATE VARIABLES USED IN ROOT UPTAKE OF WATER AND NUTRIENTS
C
C     RTDPZ,RTDP1=primary root depth (m)
C     FRTDPX=fraction of each soil layer with primary root
C     DLYR=layer thickness (m)
C     CDPTHZ=depth from soil surface to layer bottom (m)
C     SDPTH=seeding depth (m)
C     HTCTL=hypocotyledon height (m)
C     FPQ=PFT fraction of biome root mass
C     RTDNP,RTLGP=root length density,root length per plant (m m-3,m)
C     RRADL=root radius (m)
C     PORT=root porosity (m3 m-3)
C     FMPR=micropore fraction excluding macropore,rock
C     RTVLW=root aqueous volume (m3)
C     PP=plant population (n)
C     PATH=path length of water and nutrient uptake (m)
C     RTARR=root surface area/radius for uptake,diffusivity (m)
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
C     IF(N.EQ.1)THEN
C     WRITE(*,4413)'FRTDPX',I,J,NFZ,NX,NY,NZ,L,N,FRTDPX(L,NZ)
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
      PATH(N,L)=1.0/(SQRT(3.1416*(RTDNP(N,L,NZ,NY,NX)/FRTDPX(L,NZ))
     3/FMPR(L,NY,NX)))
      RTARR(N,L)=6.283*RTLGP(N,L,NZ,NY,NX)/FRTDPX(L,NZ)
      ELSE
      RRADL(N,L)=RRAD2M(N,NZ,NY,NX)
      PATH(N,L)=DLYR(3,L,NY,NX)
      RTARR(N,L)=6.283*RTLGP(N,L,NZ,NY,NX)
      ENDIF
      UPWTRM(N,L)=0.0
C     IF(NZ.EQ.1.OR.NZ.EQ.2)THEN
C     WRITE(*,4413)'RTAR',I,J,NFZ,NX,NY,NZ,L,N
C    2,RTARR(N,L),RRADL(N,L),PATH(N,L),WTRTD(N,L,NZ,NY,NX) 
C    2,RTLGP(N,L,NZ,NY,NX),RTDNP(N,L,NZ,NY,NX) 
C    2,PP(NZ,NY,NX),RTDNP(N,L,NZ,NY,NX)*PP(NZ,NY,NX)
C    2/AREA(3,L,NY,NX),FRTDPX(L,NZ),RTLG2X(N,NZ,NY,NX) 
C    2,FPQ(N,L,NZ),WTRTG(L)
C    3,RRADL(N,L),RRAD2X(N,NZ,NY,NX),RTVLW(N,L,NZ,NY,NX)
C    2,PORT(N,NZ,NY,NX),PP(NZ,NY,NX),RTLGP(N,L,NZ,NY,NX)
4413  FORMAT(A8,8I4,30E12.4)
C     ENDIF 
2000  CONTINUE
C
C     CALCULATE CANOPY WATER STATUS FROM CONVERGENCE SOLUTION FOR
C     TRANSPIRATION - ROOT WATER UPTAKE = CHANGE IN CANOPY WATER
C     CONTENT
C
C     FLWC=foliar water retention of precipitation from
C        ‘hour1.f’(m3 h-1)
C     FLWCCM=foliar water retention of precipitation (m3 t-1)
C     WVPLT=total hydrologically active canopy C (m3)
C     VOLWP,VOLPY=water volume in canopy (m3) 
C     WTLS,WTSTK=leaf+petiole,stalk mass (g C)
C     FARS=stalk sapwood thickness from ‘startq.f’(m)
C     ARSTP=canopy stalk surface area (m2)
C     VSTK=stalk volume:mass from ‘startq.f’ (m3 g-1)
C     PSILT=canopy total water potential (MPa)
C     FDMPM=minimum canopy dry matter concentration (g C g C-1) 
C        from ‘startq.f’
C     VOLWPZ=canopy water capacity at current PSILT (m3)
C     XNPHX,XNPH=time steps for water fluxes from ‘wthr.f’ (h t-1) 
C
      FLWCCM=FLWC(NZ,NY,NX)*XNPHX
      VOLWPY=VOLWP(NZ,NY,NX)
      WVPLTY=WVPLT(NZ,NY,NX)
      WVPLT(NZ,NY,NX)=AMAX1(0.0,WTLS(NZ,NY,NX)
     2+AMIN1(WTSTK(NZ,NY,NX),FARS*ARSTP(NZ,NY,NX)/VSTK))
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDMP=FDMPM+0.10*APSILT/(0.05*APSILT+2.0)
      VOLWPZ=1.0E-06*WVPLT(NZ,NY,NX)/FDMP
      VOLWPX=VOLWPZ
C     IF(NZ.EQ.1)THEN
C     WRITE(*,2221)'VOLWP',I,J,NFZ,NX,NY,NZ
C    2,WVPLTY,WVPLT(NZ,NY,NX),VOLWPY,VOLWP(NZ,NY,NX),VOLWPZ
C    3,VOLWPD,AMIN1(0.0,VOLWPD),APSILT,FDMP,FLWR(NY,NX)
2221  FORMAT(A8,6I4,20E14.6)
C     ENDIF
C
C     CANOPY HEAT CAPACITY
C
C     VHCPCP=dry canopy heat capacity (MJ K-1)
C     WVPLT=total hydrologically active canopy C (m3)
C     VSTK=stalk volume:mass from ‘startq.f’ (m3 g-1)
C     VHCPCX=wet canopy heat capacity (MJ K-1)
C     VOLWP,VOLWC=water volume in canopy,on canopy surfaces (m3)
C     VHCPYM,VHCPXM=minimum heat capacity for solving canopy 
C        water potential,aerodymamic energy exchange (MJ m-2 K-1)
C     AREA=grid cell area (m2)
C     VHCPXZ,VHCPYZ=minimum heat capacity for solving canopy
C        water potential,aerodymamic energy exchange (MJ K-1) 
C     TKCY=intermediate estimate of TKC used in convergence (K) 
C     DTKX=step change in canopy surface temperature used in
C        convergence (K) 
C
      VHCPCP=2.496*WVPLT(NZ,NY,NX)*VSTK
      VHCPCX=VHCPCP+4.19*(AMAX1(0.0,VOLWC(NZ,NY,NX))
     2+AMAX1(0.0,VOLWP(NZ,NY,NX)))
      VHCPXZ=VHCPXM*AREA(3,NU(NY,NX),NY,NX)
      VHCPYZ=VHCPYM*AREA(3,NU(NY,NX),NY,NX)
      TKCY=TKC(NZ,NY,NX)
      TKCX=TKCY
      IF(VHCPCP.GT.VHCPXZ)THEN
      DTKX=DTKX1
      ELSE
      DTKX=DTKX2
      ENDIF
C     IF(ICHKF.EQ.1.AND.NZ.EQ.1)THEN
C     WRITE(*,2123)'STARTC',I,J,NFZ,NX,NY,NZ
C    2,IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX)
C    2,TKC(NZ,NY,NX)
C    2,ARLFS(NZ,NY,NX),ARLSS(NY,NX) 
C    3,PP(NZ,NY,NX),FRADP(NZ,NY,NX),ZEROP(NZ,NY,NX)
C    4,RTDP1(1,1,NZ,NY,NX),VHCPCX,VHCPYZ
C    4,WVPLT(NZ,NY,NX),WTLF(NZ,NY,NX),WTSHE(NZ,NY,NX),WTLS(NZ,NY,NX)
C    5,WTSTK(NZ,NY,NX),ARSTP(NZ,NY,NX),VSTK
C    3,VOLWC(NZ,NY,NX),VOLWP(NZ,NY,NX)
2123  FORMAT(A8,7I4,20E12.4)
C     ENDIF
C
C     INITIALIZE CONVERGENCE SOLUTION FOR EQUILIBRATING ROOT WATER 
C     UPTAKE WITH CANOPY TRANSPIRATION
C
C     IDAY(1,=emergence date
C     VHCPCX=wet canopy heat capacity (MJ K-1)
C     VHCPYZ=minimum canopy heat capacity for solving 
C        canopy energy exchange (MJ K-1)
C     FRADP=fraction of incoming radiation received by each PFT canopy 
C        from ‘hour1.f’ 
C     RTDP1=primary root depth from soil surface (m)
C     SDPTHI=seeding depth (m)
C     CDPTHZ(0=soil surface elevation (m)
C
      IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0
     2.AND.VHCPCX.GT.VHCPYZ
     3.AND.FRADP(NZ,NY,NX).GT.ZERO2
     4.AND.RTDP1(1,1,NZ,NY,NX).GT.SDPTH(NZ,NY,NX)
     5+CDPTHZ(0,NY,NX))THEN
C
C     GRAVIMETRIC WATER POTENTIAL FROM CANOPY HEIGHT
C
C     HTSTZ=canopy height for calculating water uptake (m)
C     ZC=canopy height (m) 
C     PSILH=gravimetric water potential at HTSTZ (MPa)
C     FRADW=conducting elements of stalk relative to those of 
C        primary root
C     PSILT=canopy total water potential (MPa)
C     EMODW=wood modulus of elasticity (MPa)
C
      CNDT=0.0
      HTSTZ(NZ,NY,NX)=0.80*ZC(NZ,NY,NX)
      PSILH=-0.0098*HTSTZ(NZ,NY,NX)
      FRADW=3.75E+03*(AMAX1(0.5,1.0+PSILT(NZ,NY,NX)/EMODW))**4
C
C     SOIL AND ROOT HYDRAULIC RESISTANCES TO ROOT WATER UPTAKE
C
      DO 3880 N=1,MY(NZ,NY,NX)
      DO 3880 L=NU(NY,NX),NI(NZ,NY,NX)
C     IF(NZ.EQ.2)THEN
C     WRITE(*,2124)'ILYR',I,J,NX,NY,NZ,L,N,RTDNP(N,L,NZ,NY,NX)
C    2,CNDU(L,NY,NX),RTN1(1,L,NZ,NY,NX),RTNL(N,L,NZ,NY,NX)
C    3,THETW(L,NY,NX),ZEROP(NZ,NY,NX)
2124  FORMAT(A8,7I4,20E12.4)
C     ENDIF
C
C     VOLX,VOLWM,THETW=soil,water volume,content (m3,m3,m3 m-3)
C     RTDNP,RTLGP=root length density,root length per plant (m m-3,m)
C     CNDU=soil hydraulic conductivity for root uptake from ‘hour1.f’
C        (m2 h-1 MPa-1)
C     RTN1,RTNL=number of root,mycorrhizal primary,secondary axes
C     ILYR:1=rooted,0=not rooted
C     N:1=root,2=mycorrhizae
C
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
C     RSSL=root density effect on soil hydraulic resistance (m-1)
C     PATH=path length of water and nutrient uptake (m)
C     RRADL,RTARR=root radius,surface area/radius (m,m)
C     RSSX=soil hydraulic resistance (MPa h m-1)
C     CNDU=soil hydraulic conductivity for root uptake from ‘hour1.f’
C        (m2 h-1 MPa-1)
C
      RSSL=(LOG((PATH(N,L)+RRADL(N,L))/RRADL(N,L))
     2/RTARR(N,L)) 
      RSSX(N,L)=(RSSL/CNDU(L,NY,NX))
C
C     RADIAL ROOT RESISTANCE FROM ROOT AREA AND RADIAL RESISTIVITY
C     ENTERED IN 'READQ'
C
C     RRAD2=secondary root radius (m)
C     RTLGP=root length per plant (m)
C     PP=plant population (n)
C     RSRG=root radial resistance (MPa h m-1)
C     RSRR=root radial resistivity from PFT file (MPa h m-3)
C     VOLA,VOLWM=soil micropore,water volume (m3,m3)
C
      RTAR2=6.283*RRAD2(N,L,NZ,NY,NX)*RTLGP(N,L,NZ,NY,NX) 
      RSRG(N,L)=(RSRR(N,NZ,NY,NX)/RTAR2
     2*VOLA(L,NY,NX)/VOLWM(NPH,L,NY,NX))
C
C     ROOT AXIAL RESISTANCE FROM RADII AND LENGTHS OF PRIMARY AND
C     SECONDARY ROOTS AND FROM AXIAL RESISTIVITY ENTERED IN 'READQ'
C
C     FRAD1,FRAD2=primary,secondary root radius relative to maximum 
C        secondary radius from PFT file RRAD2M at which RSRA 
C        is defined
C     RRAD1,RRAD2=primary,secondary root radius (m)
C     RSRA=axial resistivity from PFT file (MPa h m-2)
C     DPTHZ=depth of primary root from surface (m)
C     RSR1,RSR2=axial resistance of primary,secondary roots (MPa h m-1)
C     RTLGA=average secondary root length (m)
C     RTN1,RTNL=number of primary,secondary axes
C     HTSTZ=canopy height for calculating water uptake (m)
C     FRADW=conducting elements of stalk relative to those of 
C        primary root
C     ILYR:1=rooted,0=not rooted
C
      FRAD1=(RRAD1(N,L,NZ,NY,NX)/RRAD2M(N,NZ,NY,NX))**4
      RSR1(N,L)=(RSRA(N,NZ,NY,NX)*DPTHZ(L,NY,NX) 
     2/(FRAD1*RTN1(1,L,NZ,NY,NX)/PP(NZ,NY,NX))
     3+RSRA(1,NZ,NY,NX)*HTSTZ(NZ,NY,NX)
     4/(FRADW*RTN1(1,L,NZ,NY,NX)/PP(NZ,NY,NX)))
      FRAD2=(RRAD2(N,L,NZ,NY,NX)/RRAD2M(N,NZ,NY,NX))**4
      RSR2(N,L)=(RSRA(N,NZ,NY,NX)*RTLGA(N,L,NZ,NY,NX)
     2/(FRAD2*RTNL(N,L,NZ,NY,NX)/PP(NZ,NY,NX)))
      ELSE
      ILYR(N,L)=0
      ENDIF
3880  CONTINUE
C
C     TOTAL ROOT RESISTANCE = SOIL + RADIAL + AXIAL
C
      DO 3890 N=1,MY(NZ,NY,NX)
      DO 3890 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(ILYR(N,L).EQ.1)THEN
C
C     RSRT=root radial+axial resistance (MPa h m-1)
C     RSRG=root radial resistance (MPa h m-1)
C     RSR1,RSR2=axial resistance of primary,secondary roots (MPa h m-1)
C     RSRS=total soil+root resistance (MPa h m-1)
C     CNDT=total soil+root conductance for all layers (m h-1 MPa-1)
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
      PSILX=PSILT(NZ,NY,NX)
      EPCX=0.0
      UPRTX=0.0
C
C     INITIALIZE CANOPY WATER POTENTIAL, OTHER VARIABLES USED IN 
C     ENERGY BALANCE TERMS THAT DON'T NEED TO BE RECALCULATED DURING
C     CONVERGENCE
C
C     RADCC=SW radiation absorbed by living canopy (MJ t-1)
C     RADC=total SW absorbed by living canopy from ‘hour1.f’
C        (MJ h-1)
C     FRADP=fraction of incoming radiation received by each PFT canopy 
C        from ‘hour1.f’ 
C     EMMC=canopy emissivity
C     EMMCX=EMMC x SB constant (MJ h-1 m-2 K-4) x FRADP x AREA (m2) 
C        x time step (MJ t-1 K-4)
C     AREA=grid cell area (m2)
C     THRMCXM,THRMCYM=LW radiation absorbed by canopy 
C        from sky, adjacent grid cell (MJ t-1)
C     THS=sky LW radiation from ‘wthr.f’ (MJ h-1)
C     DTHRMT=total lateral LW radiation (MJ t-1)
C     FLAIP=fraction of total canopy+standing dead radiation
C        received by each PFT canopy from ‘hour1.f’
C     CCPOLT=total nonstructural canopy C,N,P concentration (g C g C-1)  
C     CCPOLP,CZPOLP,CPPOLP=nonstructural canopy C,N,P concentration 
C        (g C g C-1)  
C     OSWT=molar mass of CCPOLT (g mol-1)
C     PAREX,PARSX=terms used to calculate canopy boundary 
C        layer conductance for latent,sensible heat from ‘hour1.f’ 
C        (m2 t-1,MJ m-1 K-1 t-1)
C     PAREZ,PARSZ=canopy boundary layer conductance for latent,sensible
C        heat (m2 t-1,MJ m-1 K-1 t-1)
C     XNPHX=time step for heat fluxes from ‘wthr.f’ (h t-1)
C
      RADCC=RADC(NZ,NY,NX)*XNPHX
      EMMCX=EMMC*2.04E-10*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
     2*XNPHX
      THRMCXM=THS(NY,NX)*FRADP(NZ,NY,NX)*XNPHX
      THRMCYM=DTHRMT(NY,NX)*FLAIP(NZ,NY,NX)*XNPHX
      CCPOLT=CCPOLP(NZ,NY,NX)+CZPOLP(NZ,NY,NX)+CPPOLP(NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      NNM=0
      RFLXCC=0.0
      EFLXCC=0.0
      SFLXCC=0.0
      HFLXCC=0.0
      VFLXCC=0.0
      THRMCC=0.0
      EVAPCC=0.0
      EPCC=0.0
      TEVCG=0.0
      TSHCG=0.0
      UPRTM=0.0
      PAREZ=FRADP(NZ,NY,NX)*PAREX(NY,NX)
      PARSZ=FRADP(NZ,NY,NX)*PARSX(NY,NX)
C
C     CONVERGENCE SOLUTION FOR CANOPY ENERGY BALANCE AND WATER UPTAKE
C
C     VHCPQ=canopy heat capacity used to calculate 
C        canopy air temperature from ‘hour1.f’ 
C        (MJ K-1)
C     FLAIP=fraction of total canopy+standing dead radiation
C        received by each PFT canopy from ‘hour1.f’
C     TKQY,TKCY=canopy air,surface temperature (K)
C     HCBFCY=heat released by canopy+standing dead combustion 
C        in previous time step (MJ t-1)
C     VOLWC=water volume on canopy surfaces (m3)
C     FLWCCM,HFLWCCM=foliar water,convective heat retention 
C        of precipitation (m3 t-1,MJ t-1) 
C
      DO 4005 M=1,NPH
      IC=0
      XC=1.0
      IF(VHCPQ(NY,NX).GT.ZEROS(NY,NX)
     2.AND.FLAIP(NZ,NY,NX).GT.ZERO2)THEN
      TKQY=TKQY+HCBFCY(NZ,NY,NX)
     2/(VHCPQ(NY,NX)*FLAIP(NZ,NY,NX))*XNPH
      ENDIF
      VOLWC(NZ,NY,NX)=VOLWC(NZ,NY,NX)+FLWCCM 
      HFLWCCM=FLWCCM*4.19*TKCY 
      DO 4000 NN=1,MXN
C
C     NET RADIATION FROM ABSORBED SW AND NET LW
C
C     VOLWPDM=change in difference between canopy water content 
C        and canopy water holding capacity (m3 t-1)
C     THRMCZM=LW emitted by canopy (MJ t-1)
C     EMMCX=EMMC x SB constant (MJ h-1 m-2 K-4) x FRADP x AREA (m2) 
C        x time step (MJ t-1 K-4)
C     TKCY=canopy surface temperature (K)
C     THRMCGM=LW radiation emitted by canopy to ground surface 
C        from ‘watsub.f’(MJ t-1)
C     DTHRMC=net LW radiation absorbed by canopy (MJ t-1) 
C     THRMCXM,THRMCYM=LW radiation absorbed by canopy 
C        from sky,to adjacent grid cell (MJ t-1)
C     RADCC=SW radiation absorbed by canopy (MJ t-1)
C     RFLXCCM=net SW+LW radiation absorbed by canopy (MJ t-1)
C
      VOLWPDM=(VOLWPZ-VOLWP(NZ,NY,NX))*XNPH
      THRMCZM=EMMCX*TKCY**4
      THRMCGM=EMMCX*(TKCY**4-TKGS(M,NY,NX)**4)
     2*FRADP(NZ,NY,NX)
      DTHRMC=THRMCXM+THRMCYM-THRMCZM-THRMCGM
      RFLXCCM=RADCC+DTHRMC 
C     IF(NZ.EQ.3.AND.ICHKF.EQ.1)THEN
C     WRITE(*,4447)'THRMCG',I,J,NFZ,M,NN,NX,NY,NZ
C    2,RFLXCCM,RADCC,DTHRMC,THRMCXM,THRMCYM,THRMCZM,THRMCGM
C    3,TKCY,FRADP(NZ,NY,NX),FLAIP(NZ,NY,NX)
C     ENDIF 
C
C     CANOPY BOUNDARY LAYER RESISTANCE 
C
C     DTKQC=air-canopy temperature difference (K)
C     TKAM,TKQC,TKCY=air,canopy air,surface temperature (K)
C     RI=Richardson number
C     RIX,RIY=minimum, maximum values used to calculate Richardson’s
C        number from ‘starts.f’
C     RABC=canopy boundary layer resistance (h m-1) 
C     RABX=biome canopy isothermal boundary layer resistance
C        from ‘hour1.f’ (h m-1)
C     RABM,RABX=minimum, maximum values used to calculate  boundary
C        layer resistances from ‘starts.f’ (h m-1)
C     RAGC=canopy aerodynamic resistance below 
C        PFT canopy height from ‘watsub.f’ (h m-1)
C     ZC,ZT,ZR=PFT canopy,biome,surface roughness height (m) 
C     RACC=canopy aerodynamic resistance between PFT canopy height
C        and maximum canopy height (h m-1)
C     RACG,RACGX=canopy aerodynamic resistance below 
C        maximum canopy height height from ‘watsub.f’ (h m-1)
C     RA=canopy boundary layer+aerodynamic resistance (h m-1)
C     PAREX,PARSX=terms used to calculate total canopy boundary 
C        layer conductance from ‘hour1.f’ (m2 t-1,MJ m-1 K-1 t-1)
C     FRADT=fraction of radiation received by all PFT canopies
C        from ‘hour1.f’
C     PAREZ,PARSZ=terms used to calculate PFT canopy boundary 
C        layer conductance (m2 t-1,MJ m-1 K-1 t-1)
C     RAHC=canopy surface resistance adjusted for canopy 
C        aerodynamic resistance (h m-1)
C     RAH,RAE=surface isothermal resistance to sensible,latent heat 
C       (h m-1) 
C     DTKC=difference between canopy,surface air temperature (K)
C     PAREC,PARSC=canopy surface conductance to sensible,latent heat 
C        (m3 t-1,MJ K-1 t-1)
C 
      DTKQC=TKAM(NY,NX)-TKQY 
      RI=AMAX1(RIX,AMIN1(RIY
     2,RIBX(NY,NX)/TKAM(NY,NX)*DTKQC))
      RIC=1.0-10.0*RI
      RABC=AMIN1(RABZ,AMAX1(RABM,RABX(NY,NX)/RIC))
      RACC=AMAX1(0.0,RACG(M,NY,NX)-RAGC(M,NZ,NY,NX))
      RA(NZ,NY,NX)=RABC+RACC
      PARSGC=PARSX(NY,NX)/RAGC(M,NZ,NY,NX)*FRADT(NY,NX)
      PAREGC=PAREX(NY,NX)/RAGC(M,NZ,NY,NX)*FRADT(NY,NX)
      DTKC=TKQY-TKCY
      RI=AMAX1(RIX,AMIN1(RIY
     2,RIBX(NY,NX)/TKQY*DTKC))
      RIS=1.0-3.2*RI
      RAHC=AMIN1(RABZ,AMAX1(RABM,RAH/RIS)) 
      PARSC=PARSZ/RAHC
      PAREC=PAREZ/(RAHC+RAE) 
C
C     CANOPY WATER AND OSMOTIC POTENTIALS
C
C     PSILT=canopy total water potential (K)
C     FDMP=canopy dry matter content (g C g C-1) 
C     FDMPM=minimum canopy dry matter concentration (g C g C-1)
C        from ‘startq.f’
C     OSMO=osmotic potential at PSILT=0 from PFT file (MPa)
C     OSWT=molar mass of CCPOLT (g mol-1)
C     TKCY=canopy surface temperature (K)
C     CCPOLT=total nonstructural canopy C,N,P concentration (g g C-1)  
C     PSILO,PSILG=canopy osmotic,turgor water potential (MPa)
C     VOLWPZ=canopy water content at current PSILT (m3)
C
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDMP=FDMPM+0.10*APSILT/(0.05*APSILT+2.0)
      VOLWPZ=1.0E-06*WVPLT(NZ,NY,NX)/FDMP
      PSILO(NZ,NY,NX)=FDMP/FDMPM*OSMO(NZ,NY,NX)
     2-8.3143*TKCY*FDMP*(CCPOLT/OSWT+CSALTP(NZ,NY,NX))
      PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
C
C     CANOPY STOMATAL RESISTANCE
C
C     RCS=shape parameter for RC vs PSILG from PFT file
C     PSILG=canopy turgor potential (MPa)
C     RC=canopy stomatal resistance (h m-1)
C     RSMN=minimum RC at PSILT=0 from ‘stomate.f’ (h m-1)
C     RSMH=cuticular resistance from PFT file and ‘startq.f’ (h m-1)
C
      WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      RC(NZ,NY,NX)=RSMN(NZ,NY,NX)+(RSMH(NZ,NY,NX)-RSMN(NZ,NY,NX))*WFNC
C
C     CANOPY VAPOR CONCENTRATION AND EVAPORATION OF INTERCEPTED WATER
C     OR TRANSPIRATION OF UPTAKEN WATER
C
C     TKCY=canopy surface temperature (K)
C     VPCY=canopy surface vapor concentration (m3 m-3)
C     PSILT=canopy total water potential (MPa)
C     EX=canopy-atmosphere vapor flux (m3 t-1)
C     PAREC,PARSC=canopy surface conductance to sensible,latent heat 
C        (m3 t-1,MJ K-1 t-1)
C     VPQY=canopy vapor concentration (m3 m-3)
C     EVAPCCM=canopy surface evaporation (m3 t-1)
C     VOLWC=water volume on canopy surfaces (m3)
C     EPCCMX,EPCCM=canopy actual,net transpiration (m3 t-1) 
C     RAHC=canopy surface resistance adjusted for canopy 
C        aerodynamic resistance (h m-1)
C     RAE=canopy surface resistance to latent heat (h m-1)
C     EFLXCCM=canopy surface latent heat flux (MJ t-1)
C     VFLXCCM=convective heat flux from EFLXCCM (MJ t-1)
C     DTKC=difference between canopy,surface air temperature (K)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C     VOLWPDM=change in difference between canopy water content 
C        and canopy water holding capacity (m3 t-1)
C     SFLXCCM=canopy surface sensible heat flux (MJ t-1)
C
      VPCY=2.173E-03/TKCY
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKCY))
     3*EXP(18.0*PSILT(NZ,NY,NX)/(8.3143*TKCY))
      EX=PAREC*(VPQY-VPCY)
      IF(EX.GT.0.0)THEN
      EVAPCCM=EX
      EX=0.0
      ELSE 
      EVAPCCM=AMAX1(EX,-AMAX1(0.0,VOLWC(NZ,NY,NX)*XNPG))
      EX=EX-EVAPCCM
      ENDIF
      EPCCMX=EX*(RAHC+RAE)/(RAHC+RC(NZ,NY,NX))
      EPCCM=EPCCMX+VOLWPDM 
      EFLXCCM=(EPCCMX+EVAPCCM)*VAP
      VFLXCCM=EVAPCCM*4.19*TKCY
      SFLXCCM=PARSC*DTKC
C
C     CANOPY STORAGE HEAT FROM RN,LE AND CONVECTIVE HEAT
C     FLUXES
C
C     HFLXCCM=canopy storage heat flux (MJ t-1) 
C     RFLXCCM=net SW+LW radiation absorbed by canopy (MJ t-1)
C     EFLXCCM=canopy surface latent heat flux (MJ t-1)
C     SFLXCCM=canopy surface sensible heat flux (MJ t-1)
C     VFLXCCM=convective heat flux from EFLXCCM (MJ t-1)
C     HFLWCCM=convective heat flux from precipitation to canopy(MJ t-1)
C
      HFLXCCM=RFLXCCM+EFLXCCM+SFLXCCM+VFLXCCM+HFLWCCM
C
C     SOLVE FOR CANOPY TEMPERATURE CAUSED BY SENSIBLE+STORAGE HEAT
C
C     VHCPCX=wet canopy heat capacity (MJ K-1)
C     VHCPCC=canopy heat capacity (MJ K-1)
C     EVAPCCM=canopy surface evaporation (m3 t-1)
C     FLWCCM=foliar water retention of precipitation (m3 t-1)
C     TKCY=canopy surface temperature (K)
C     HFLXCCM=canopy storage heat flux (MJ t-1) 
C     DTKX=step change in canopy surface temperature used in
C        convergence (K)
C     XC,IC=magnitude,direction of change in canopy surface 
C        temperature for next cycle 
C
      VHCPCC=VHCPCX+4.19*(EVAPCCM+FLWCCM)
      IF(VHCPCC.GT.VHCPYZ)THEN
      TKCX=TKCY
      TKCZ=(TKCX*VHCPCX+HFLXCCM)/VHCPCC
      TKCY=TKCX+XC*DTKX*(TKCZ-TKCX)
      ELSE
      TKCZ=TKCX
      TKCY=TKCX
      ENDIF
      IF((IC.EQ.0.AND.TKCY.GT.TKCX).OR.(IC.EQ.1.AND.TKCY.LT.TKCX))THEN
      XC=0.5*XC
      ENDIF
      IF(TKCY.GT.TKCX)THEN
      IC=1
      ELSE
      IC=0
      ENDIF
C     IF(NZ.EQ.1)THEN
C     WRITE(*,4444)'TKCY',I,J,NFZ,M,NN,NX,NY,NZ 
C    2,XC,TKCX,TKCY,TKCZ,TKCY-TKCX,DTKX,TKC(NZ,NY,NX),TKQY
C    2,TKQC(NZ,NY,NX),TKQG(M,NY,NX),VHCPCX,VHCPCC,VHCPXZ
C    3,FRADP(NZ,NY,NX),WVPLT(NZ,NY,NX),EX,VPCY,VPQY  
C    3,VPQC(NZ,NY,NX),PAREC,PARSC 
C    3,EPCCM,UPRTM,VOLWPDM,VOLWPZ,VOLWP(NZ,NY,NX),PSILT(NZ,NY,NX)
C    3,FLWCCM,RC(NZ,NY,NX),PSILG(NZ,NY,NX),PSILO(NZ,NY,NX)
C    3,FDMP,FDMPM,CCPOLT,OSWT,CCPOLP(NZ,NY,NX),CPOOLP(NZ,NY,NX)
C    3,RADCC,DTHRMC,THRMCXM,THRMCYM,THRMCZM,THRMCGM 
C    2,TKAM(NY,NX),TKS(0,NY,NX),TKS(NU(NY,NX),NY,NX)
C    2,RAHC,RAE,RAH,RA(NZ,NY,NX),RABC,RACC,RABX(NY,NX),RI,RIC
C    3,RACG(M,NY,NX),RAGC(M,NZ,NY,NX)
C    3,RIBX(NY,NX),TKAM(NY,NX),DTKQC,DTKC,RIS 
C    4,HCBFCY(NZ,NY,NX),VHCPQ(NY,NX),FRADP(NZ,NY,NX)
C    3,HFLXCCM,RFLXCCM,EFLXCCM,SFLXCCM,VFLXCCM,HFLWCCM 
C    3,HFLXCC,RFLXCC,EFLXCC,SFLXCC,VFLXCC,HFLWCC 
C    2,RADC(NZ,NY,NX),FRADP(NZ,NY,NX) 
C    3,THS(NY,NX),THRMCCM,WTLS(NZ,NY,NX) 
C    2,CCPOLT,OSWT,CCPOLP(NZ,NY,NX),CPOOLP(NZ,NY,NX)
C    4,DCO2(NZ,NY,NX),AREA(3,NU(NY,NX),NY,NX),WTLS(NZ,NY,NX)
C    2,PSILT(NZ,NY,NX),PSILG(NZ,NY,NX),RAZ(NZ,NY,NX),RI
C    3,ARLFV(1,NZ,NY,NX),ARSTV(1,NZ,NY,NX)
C    4,WVPLT(NZ,NY,NX),VSTK,HCBFCY(NZ,NY,NX) 
4444  FORMAT(A8,8I4,80E14.6)
C     ENDIF
C
C     IF CONVERGENCE FOR ENERGY EXCHANGE IS ACHIEVED
C
C     TKCY,TKCX=intermediate,previous estimate of TKC used in
C        convergence (K) 
C     DIFFT=acceptance criteria for temperature change in convergence
C        solution for water uptake 
C     PSILC=canopy water potential adjusted for canopy height (MPa)
C     PSILT=canopy total water potential (MPa)
C     PSILH=gravimetric water potential at HTSTZ (MPa)
C
      IF(ABS(TKCY-TKCX).LT.DIFFT)THEN 
      UPRTM=0.0
      PSILC=PSILT(NZ,NY,NX)-PSILH
C
C     ROOT WATER UPTAKE FROM SOIL-CANOPY WATER POTENTIALS,
C     SOIL + ROOT HYDRAULIC RESISTANCES
C
C     ILYR:1=rooted,0=not rooted
C     UPWTRM=root water uptake from soil layer (m3 t-1)
C     VOLWU,VOLPU=water volume available for uptake,air volume (m3)
C     FPQ=PFT fraction of biome root mass
C     PSILC=canopy water potential adjusted for canopy height (MPa)
C     PSIST1=total soil water potential PSIST adjusted for 
C        surface elevation (MPa)
C     RSRS=total soil+root resistance (MPa h m-1)
C     PP=PFT population (n)
C     UPRTM=total water uptake from soil profile (m3 t-1)     
C     XNPHX=time step for water fluxes from ‘wthr.f’ (h t-1)
C
      DO 4200 N=1,MY(NZ,NY,NX)
      DO 4200 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(ILYR(N,L).EQ.1)THEN
      UPWTRX=VOLWU(L)*FPQ(N,L,NZ)*XNPHX
      UPWTRM(N,L)=AMAX1(AMIN1(0.0,-UPWTRX)
     2,AMIN1((PSILC-PSIST1(L))/RSRS(N,L)*PP(NZ,NY,NX)*XNPHX,UPWTRX))
      IF(UPWTRM(N,L).GT.0.0)THEN
      UPWTRM(N,L)=0.1*UPWTRM(N,L)
      ENDIF
      UPRTM=UPRTM+UPWTRM(N,L)
      ELSE
      UPWTRM(N,L)=0.0
      ENDIF
C     IF(NZ.EQ.1)THEN
C     WRITE(*,6565)'UPRTM',I,J,NFZ,M,NX,NY,NZ,NN,N,L,ILYR(N,L)
C    2,UPRTM,UPWTRX,UPWTRM(N,L),PSILC,PSIST1(L),PSISM(L,NY,NX)
C    2,(PSILC-PSIST1(L))/RSRS(N,L)*PP(NZ,NY,NX)*XNPHX,VOLWU(L)
C    3,RSRS(N,L),RSSX(N,L)
C    2,RSRT(N,L),RSRG(N,L),RSR1(N,L),RSR2(N,L),PSILH,RTAR2 
C    3,RSRR(N,NZ,NY,NX),VOLA(L,NY,NX),VOLWM(NPH,L,NY,NX)
C    4,VOLWU(L),THETY(L,NY,NX),VOLY(L,NY,NX),FPQ(N,L,NZ)
6565  FORMAT(A8,11I4,30E14.6)
C     ENDIF
4200  CONTINUE
C
C     TEST TRANSPIRATION - ROOT WATER UPTAKE VS. CHANGE IN CANOPY
C     WATER STORAGE
C
C     PSILT,PSILX=current,previous canopy water potential (MPa)
C     VOLWPZ,VOLWPX=current,previous canopy water potential (m3)
C     AREA=grid cell area (m2) 
C     RSSZ=change in PSILT vs VOLWPZ (MPa m-1)
C     EPCCM,EPCX=current,previous canopy net transpiration (m3 t-1) 
C     UPRTM,UPRTX=current,previous water uptake from soil profile 
C        (m3 t-1)
C     RSSU=change in PSILT vs transpiration or uptake (MPa m-1)    
C     CNDT=total soil+root conductance for all layers (m h-1 MPa-1)
C
      DPSILT=ABS(PSILT(NZ,NY,NX)-PSILX)
      DVOLWP=ABS(VOLWPZ-VOLWPX)/AREA(3,NUM(NY,NX),NY,NX)
      IF(DVOLWP.GT.ZERO2.AND.DPSILT.GT.ZERO)THEN
      RSSZ=DPSILT/DVOLWP 
      ELSEIF(CNDT.GT.ZERO)THEN
      RSSZ=1.0/CNDT
      ELSE
      RSSZ=1.0E+03
      ENDIF
      DEPCCM=AMAX1(ABS(EPCCM-EPCX),ABS(UPRTM-UPRTX))
     2/AREA(3,NUM(NY,NX),NY,NX) 
      IF(DEPCCM.GT.ZERO2.AND.DPSILT.GT.ZERO)THEN
      RSSU=DPSILT/DEPCCM
      ELSEIF(CNDT.GT.ZERO)THEN
      RSSU=1.0/CNDT
      ELSE
      RSSU=1.0E+03
      ENDIF
C
C     CHANGE IN CANOPY WATER POTENTIAL REQUIRED TO BRING AGREEMENT
C     BETWEEN TRANSPIRATION - ROOT WATER UPTAKE AND CHANGE IN CANOPY
C     WATER STORAGE
C
C     DIFFZ,DIFFU=change in canopy water content,transpiration-uptake
C        (m3 t-1)
C     VOLWPZ,VOLWPX=current,previous canopy water content (m3)
C     EPCCM,EPCX=current,previous canopy net transpiration (m3 t-1) 
C     DIFFD=difference in canopy water content and transpiration-uptake
C        (m3 m-2 t-1)
C     AREA=grid cell area (m2) 
C     DPSI=change in canopy water potential for next convergence cycle
C        (MPa t-1)
C     RSSZ=change in canopy water potential vs change in canopy water
C        content (MPa m-1)
C     RSSU=change in canopy water potential vs change in transpiration 
C        (MPa m-1)
C
      DIFFZ=VOLWPZ-VOLWP(NZ,NY,NX)
      DIFFU=EPCCM-UPRTM
      DIFFD=(DIFFU-DIFFZ)/AREA(3,NUM(NY,NX),NY,NX)
      DPSI=AMIN1(RSSZ,RSSU)*DIFFD 
C     IF(I.EQ.181.AND.NZ.EQ.4)THEN
C     WRITE(*,2222)'PSI',I,J,NFZ,NX,NY,NZ,M,NN
C    2,DPSI,PSILT(NZ,NY,NX),PSILX,UPRTM,UPRTX,EPCCM,EPCX,VOLWPDM
C    3,VOLWPZ,VOLWPX,VOLWP(NZ,NY,NX),VOLWC(NZ,NY,NX)
C    3,EVAPCCM,DIFFU,DIFFZ,DIFFD
C    3,ABS(PSILT(NZ,NY,NX)-PSILX),DPSILT
C    3,DVOLWP,DEPCCM,ABS(EPCCM-EPCX),ABS(UPRTM-UPRTX) 
C    2,RSSZ,RSSU,1.0/CNDT 
C    3,RC(NZ,NY,NX),FRADP(NZ,NY,NX)
C    3,PAREC,VPQY,VPCY,TKCY
C    5,FDMP,CCPOLT,OSWT,RC(NZ,NY,NX),RAHC 
C    5,((UPWTRM(N,L),L=1,8),N=1,1)
C    6,((RSRS(N,L),L=1,8),N=1,1)
C    7,(PSIST1(L),L=1,8)
2222  FORMAT(A8,8I4,80E12.4)
C     ENDIF
      EPCX=EPCCM
      UPRTX=UPRTM
      VOLWPX=VOLWPZ
      PSILX=PSILT(NZ,NY,NX)
      PSILT(NZ,NY,NX)=AMIN1(0.5*PSILT(NZ,NY,NX)
     2,PSILT(NZ,NY,NX)+0.5*DPSI)
      IF(NN.GE.MXN.OR.ABS(DIFFD).LT.DIFFX)GO TO 4250
      GO TO 4000
4250  IF(ICHK.EQ.1)THEN
      NNM=MAX(NN,NNM)
      GO TO 4500
      ELSE
      ICHK=1
      IF(M.EQ.NPH)THEN
C
C     RESET MINIMUM STOMATAL RESISTANCE IN ‘STOMATE.F’ 
C     BEFORE FINAL ITERATION
C
      CALL STOMATE(I,J,NFZ,NZ,NY,NX)
      ENDIF
      ENDIF
      ENDIF
4000  CONTINUE
4500  CONTINUE
C
C     CALCULATE CANOPY AIR TEMPERATURE, VAPOR CONCENTRATION
C
C     VHCPCP=dry canopy heat capacity (MJ K-1)
C     SFLXA,EVAPA=canopy-atmosphere sensible heat flux,vapor flux
C        (MJ t-1,m3 t-1)
C     DTKQC=air-canopy temperature difference (K)
C     PAREZ,PARSZ=canopy boundary layer conductance 
C        (m2 t-1,MJ m-1 K-1 t-1)
C     RA=canopy boundary layer+aerodynamic resistance (h m-1)
C     VHCPQ,EVAPQ=canopy heat capacity,volume used to calculate 
C        canopy air temperature,vapor concentration from ‘hour1.f’ 
C        (MJ K-1,m3)
C     FLAIP=fraction of total canopy+standing dead radiation
C        received by each PFT canopy from ‘hour1.f’
C     VPAM,VPQY=atmosphere,canopy vapor concentration (m3 m-3)
C     TSHCNM,TEVCNM=canopy net sensible heat flux,vapor flux 
C        (MJ t-1,m3 t-1)
C     SFLXCCM=canopy surface sensible heat flux (MJ t-1)
C     EPCCM,EVAPCCM=canopy net transpiration,evaporation (m3 t-1)
C     TSHCGM,TEVCGM=ground surface-canopy sensible heat flux, 
C        vapor flux from ‘watsub.f’ (MJ t-1,m3 t-1)
C     TSHCYM,TEVCYM=total lateral canopy sensible heat flux,vapor flux
C        (MJ t-1,m3 t-1)
C     TKQY=canopy air temperature (K)
C     VPQY,VPSY=canopy,saturated vapor concentration (m3 m-3)
C     VOLWP,VOLWC=water volume in canopy,on canopy surfaces (m3)
C     XNPHX=time step for energy fluxes from ‘wthr.f’ (h t-1)
C
      IF(VHCPCP.GT.VHCPXZ.AND.FLAIP(NZ,NY,NX).GT.1.0E-03)THEN
      SFLXA=DTKQC*PARSZ/RA(NZ,NY,NX)*FLAIP(NZ,NY,NX)
      DVPQA=VPAM(NY,NX)-VPQY
      EVAPA=DVPQA*AMIN1(PAREZ/RA(NZ,NY,NX)
     2,EVAPQ(NY,NX)*FLAIP(NZ,NY,NX))
      TSHCNM=SFLXA-SFLXCCM 
      TEVCNM=EVAPA-EPCCM-EVAPCCM 
      TSHCGM=PARSGC*(TKQY-TKQG(M,NY,NX))*FLAIP(NZ,NY,NX)
      TEVCGM=PAREGC*(VPQY-VPQG(M,NY,NX))*FLAIP(NZ,NY,NX)
      TSHCYM=DSHQT(NY,NX)*FLAIP(NZ,NY,NX)*XNPHX 
      TEVCYM=DVPQT(NY,NX)*FLAIP(NZ,NY,NX)*XNPHX
      TKQY=TKQY+(TSHCNM-TSHCGM-TSHCYM)
     2/(VHCPQ(NY,NX)*FLAIP(NZ,NY,NX)) 
      VPSY=2.173E-03/TKQY
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKQY))
      VPQY=AMAX1(0.0,AMIN1(VPSY,VPQY+(TEVCNM-TEVCGM-TEVCYM)
     2/(EVAPQ(NY,NX)*FLAIP(NZ,NY,NX))))
      ELSE
      TSHCNM=0.0
      TEVCNM=0.0 
      TSHCGM=0.0 
      TEVCGM=0.0
      TSHCYM=0.0 
      TEVCYM=0.0
      TKQY=TKAM(NY,NX)
      VPQY=VPAM(NY,NX)
      ENDIF   
      VOLWC(NZ,NY,NX)=VOLWC(NZ,NY,NX)+EVAPCCM
      VOLWP(NZ,NY,NX)=VOLWP(NZ,NY,NX)+EPCCM-UPRTM
C
C     AGGREGATE CANOPY ENERGY BALANCE FOR ALL M
C
C     RFLXCC=net SW+LW radiation absorbed by canopy (MJ t-1)
C     EFLXCC=canopy surface latent heat flux (MJ t-1)
C     SFLXCC=canopy surface sensible heat flux (MJ t-1)
C     HFLXCC=canopy storage heat flux (MJ t-1) 
C     VFLXCC=convective heat flux from EFLXCC (MJ t-1)
C     EPCC,EVAPCC=canopy surface transpiration,evaporation (m3 t-1)
C     TSHCG,TEVCG=ground surface-canopy sensible heat flux, 
C        vapor flux from ‘watsub.f’ (MJ t-1,m3 t-1)
C     VHCPCX,VHCPCP=wet,dry canopy heat capacity (MJ K-1)
C
      RFLXCC=RFLXCC+RFLXCCM
      EFLXCC=EFLXCC+EFLXCCM
      SFLXCC=SFLXCC+SFLXCCM
      HFLXCC=HFLXCC+HFLXCCM
      VFLXCC=VFLXCC+VFLXCCM
      THRMCC=THRMCC+THRMCGM
      EVAPCC=EVAPCC+EVAPCCM
      EPCC=EPCC+EPCCM
      TEVCG=TEVCG+TEVCGM
      TSHCG=TSHCG+TSHCGM
      VHCPCX=VHCPCP+4.19*(AMAX1(0.0,VOLWC(NZ,NY,NX))
     2+AMAX1(0.0,VOLWP(NZ,NY,NX)))
C     IF(NZ.EQ.1)THEN
C     WRITE(*,3114)'TKQC',I,J,NFZ,M,NN,NX,NY,NZ 
C    2,TKQY,TKAM(NY,NX),VPQY,VPAM(NY,NX),TKCY,TKC(NZ,NY,NX)
C    3,TKQC(NZ,NY,NX),TKQG(M,NY,NX),TSHCNM,TSHCGM,TSHCYM,SFLXCCM,SFLXA
C    4,EPCCM,EPCC
C    2,RABX(NY,NX),RIBX(NY,NX),RA(NZ,NY,NX),RABC,RACC,RI
C    3,RACG(M,NY,NX),RAGC(M,NZ,NY,NX),ZC(NZ,NY,NX)
C    4,ZG(NZ,NY,NX),ZT(NY,NX),DTKQC
C    3,VPSY,VPQC(NZ,NY,NX),TEVCNM,TEVCN,EPCCM,UPRTM,EVAPCCM 
C    2,EVAPA,FLAIP(NZ,NY,NX),VHCPCP,VHCPXZ
C    3,VHCPQ(NY,NX),FRADP(NZ,NY,NX),HCBFCY(NZ,NY,NX)
C    4,FRADP(NZ,NY,NX),FRADQ(NZ,NY,NX),VOLWC(NZ,NY,NX),TSHCN,TEVCN 
C    3,RADC(NZ,NY,NX),RA(NZ,NY,NX),RADP(NZ,NY,NX),RADQ(NZ,NY,NX)
3114  FORMAT(A8,8I4,60E12.4)
C     ENDIF
C
C     AGGREGATE ROOT WATER UPTAKE FOR ALL M
C
C     UPWTR=root water uptake from soil layer (m3 t-1)
C
      DO 4201 N=1,MY(NZ,NY,NX)
      DO 4201 L=NU(NY,NX),NI(NZ,NY,NX)
      UPWTR(N,L,NZ,NY,NX)=UPWTR(N,L,NZ,NY,NX)+UPWTRM(N,L)
4201  CONTINUE
C     WRITE(*,4443)'VOLWP',I,J,NFZ,M,NX,NY,NZ,NN,NI(NZ,NY,NX)
C    2,VOLWC(NZ,NY,NX),EVAPCCM,FLWCCM,FLWC(NZ,NY,NX)
C    3,VOLWP(NZ,NY,NX),EPCCM,UPRTM
C    4,EPCC,EVAPCC
C    5,((UPWTRM(N,L),L=NU(NY,NX),NI(NZ,NY,NX)),N=1,MY(NZ,NY,NX))
C    6,(UPWTR(N,L,NZ,NY,NX),L=NU(NY,NX),NI(NZ,NY,NX))    
4443  FORMAT(A8,9I4,60F16.8)
4005  CONTINUE
C
C     FINAL CANOPY TEMPERATURE 
C
C     TKQC,VPQC=final estimate of canopy air temperature,
C        vapor concentration (K,m3 m-3)
C     TKC=final estimate of canopy surface temperature (K)
C
      IF(NNM.LT.MXN)THEN
      TKQC(NZ,NY,NX)=TKQY
      VPQC(NZ,NY,NX)=VPQY
      TKC(NZ,NY,NX)=TKCY
      TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-273.15
      ELSE
C
C     IF CONVERGENCE NOT ACHIEVED (RARE), SET DEFAULT 
C     TEMPERATURES, ENERGY FLUXES, WATER POTENTIALS, RESISTANCES
C
C     TKQC,VPQC=default estimate of canopy air temperature,
C        vapor concentration (K,m3 m-3)
C     TKC=default estimate of canopy surface temperature (K)
C     THRMCC=default estimate of canopy LW emission (MJ t-1) 
C     PSILT,PSILO,PSILG=default estimate of canopy total,
C        osmotic,turgor potential (MPa)
C     RC,RA=default estimate of canopy stomatal,aerodynamic
C        +boundary layer resistance (h m-1)
C     PSIRT,PSIRO,PSIRG=default estimate of root total,
C        osmotic,turgor potential (MPa)
C     TKS=soil temperature (K)
C
      IF(ABS(DIFFD).GT.5.0*DIFFX)THEN 
      WRITE(*,9999)IYRC,I,J,NFZ,NX,NY,NZ
9999  FORMAT('CONVERGENCE FOR WATER UPTAKE NOT ACHIEVED ON   ',7I6)
      ENDIF
C     TKQC(NZ,NY,NX)=TKQ(NY,NX)
C     VPQC(NZ,NY,NX)=VPQ(NY,NX)
C     TKC(NZ,NY,NX)=TKQ(NY,NX)
C     TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-273.15
      EMMCX=EMMC*2.04E-10*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
     2*XNPHX 
      THRMCC=EMMCX*TKC(NZ,NY,NX)**4 
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDMP=FDMPM+0.10*APSILT/(0.05*APSILT+2.0)
      CCPOLT=CCPOLP(NZ,NY,NX)+CZPOLP(NZ,NY,NX)+CPPOLP(NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSILO(NZ,NY,NX)=FDMP/FDMPM*OSMO(NZ,NY,NX)
     2-8.3143*TKC(NZ,NY,NX)*FDMP*(CCPOLT/OSWT+CSALTP(NZ,NY,NX))
      PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
      WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      RC(NZ,NY,NX)=RSMN(NZ,NY,NX)+(RSMH(NZ,NY,NX)
     2-RSMN(NZ,NY,NX))*WFNC
      RA(NZ,NY,NX)=RAB(NY,NX)
      VHCPCC=2.496*WVPLT(NZ,NY,NX)*VSTK
      DO 4290 N=1,MY(NZ,NY,NX)
      DO 4290 L=NU(NY,NX),NI(NZ,NY,NX)
      PSIRT(N,L,NZ,NY,NX)=PSIST1(L)
      APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
      FDMR=FDMPM+0.10*APSIRT/(0.05*APSIRT+2.0)
      CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSIRO(N,L,NZ,NY,NX)=FDMR/FDMPM*OSMO(NZ,NY,NX)
     2-8.3143*TKS(L,NY,NX)*FDMR*(CCPOLT/OSWT+CSALTR(N,L,NZ,NY,NX))
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)
     2-PSIRO(N,L,NZ,NY,NX))
4290  CONTINUE
      ENDIF
C     IF(I.GT.175)THEN
C     WRITE(*,4444)'EBC',I,J,NFZ,NX,NY,NZ,NN
C    2,RFLXCC,EFLXCC,SFLXCC,HFLXCC
C    2,RFLXCC+EFLXCC+SFLXCC+HFLXCC
C    3,PARSC,TKCY
C    3,TKCX*VHCPCX-TKCY*VHCPCC,VFLXCC,HFLWCC
C    4,VOLWP(NZ,NY,NX),EPCC,UPRTM,VHCPCX,ZEROP2(NZ,NY,NX) 
C     ENDIF
C
C     ROOT TOTAL, OSMOTIC AND TURGOR WATER POTENTIALS
C     CALCULATED AFTER CONVERGENCE
C
C     ILYR:1=rooted,0=not rooted
C     PSIRT,PSIRO,PSIRG=root total,osmotic,turgor potential (MPa)
C     PSIST1=total soil water potential PSIST adjusted for 
C        surface elevation (MPa)
C     CCPOLT=total nonstructural root C,N,P concentration 
C        (g C,N,P g C-1)  
C     FDMR=dry matter content (g C g C-1)
C     FDMPM=minimum canopy dry matter concentration (g C g C-1) 
C        from ‘startq.f’
C     OSMO=osmotic potential at PSIRT=0 from PFT file (MPa)
C     TKS=soil temperature (K)
C     
      DO 4505 N=1,MY(NZ,NY,NX)
      DO 4510 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(ILYR(N,L).EQ.1)THEN
      PSIRT(N,L,NZ,NY,NX)=AMIN1(0.0,(PSIST1(L)*RSRT(N,L)
     2+PSILT(NZ,NY,NX)*RSSX(N,L))/RSRS(N,L))
      APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
      FDMR=FDMPM+0.10*APSIRT/(0.05*APSIRT+2.0)
      CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSIRO(N,L,NZ,NY,NX)=FDMR/FDMPM*OSMO(NZ,NY,NX)
     2-8.3143*TKS(L,NY,NX)*FDMR*(CCPOLT/OSWT+CSALTR(N,L,NZ,NY,NX))
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)
     2-PSIRO(N,L,NZ,NY,NX))
      ELSE
      PSIRT(N,L,NZ,NY,NX)=PSIST1(L)
      APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
      FDMR=FDMPM+0.10*APSIRT/(0.05*APSIRT+2.0)
      CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSIRO(N,L,NZ,NY,NX)=FDMR/FDMPM*OSMO(NZ,NY,NX)
     2-8.3143*TKS(L,NY,NX)*FDMR*(CCPOLT/OSWT+CSALTR(N,L,NZ,NY,NX))
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)
     2-PSIRO(N,L,NZ,NY,NX))
      ENDIF
C     IF(NZ.EQ.1)THEN
C     WRITE(*,1256)'PSIRT',I,J,NFZ,NX,NY,NZ,NN,PSIRT(N,L,NZ,NY,NX)
C    2,PSIST1(L),RSRT(N,L),PSILT(NZ,NY,NX),RSSX(N,L),RSRS(N,L)
C    3,RSRG(N,L),RSR1(N,L),RSR2(N,L),RTAR2,VOLWM(NPH,L,NY,NX)
C    4,PSIRG(N,L,NZ,NY,NX),PSIRO(N,L,NZ,NY,NX),PSIST(L,NY,NX)
C    5,ALT(NY,NX)
1256  FORMAT(A8,7I4,20E12.4)
C     ENDIF
4510  CONTINUE
4505  CONTINUE
C
C     DEFAULT VALUES IF PLANT CANOPY IS TOO SMALL FOR ENERGY FLUX 
C     CALCULATIONS AS DEFINED UNDER ‘INITIALIZE CONVERGENCE SOLUTION’
C     ABOVE
C
      ELSE
      VOLWC(NZ,NY,NX)=VOLWC(NZ,NY,NX)+FLWC(NZ,NY,NX)*XNFH
      RFLXCC=0.0
      EFLXCC=0.0
      SFLXCC=0.0
      HFLXCC=0.0
      VFLXCC=0.0
      THRMCC=0.0
      EVAPCC=0.0
      EPCC=0.0
      TEVCG=0.0
      TSHCG=0.0
      IF(ZC(NZ,NY,NX).GE.DPTHS(NY,NX)-ZERO)THEN
      TKQC(NZ,NY,NX)=TKAM(NY,NX)
      TKC(NZ,NY,NX)=TKAM(NY,NX)
      ELSE
      TKQC(NZ,NY,NX)=TKW(1,NY,NX)
      TKC(NZ,NY,NX)=TKW(1,NY,NX)
      ENDIF
      TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-273.15
      VPQC(NZ,NY,NX)=VPAM(NY,NX)
      EMMCX=EMMC*2.04E-10*FRADP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
     2*XNPHX
      THRMCC=EMMCX*TKC(NZ,NY,NX)**4 
      PSILT(NZ,NY,NX)=PSIST1(NG(NZ,NY,NX))
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDMP=FDMPM+0.10*APSILT/(0.05*APSILT+2.0)
      CCPOLT=CCPOLP(NZ,NY,NX)+CZPOLP(NZ,NY,NX)+CPPOLP(NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSILO(NZ,NY,NX)=FDMP/FDMPM*OSMO(NZ,NY,NX)
     2-8.3143*TKCY*FDMP*(CCPOLT/OSWT+CSALTP(NZ,NY,NX))
      PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
      WFNC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      RC(NZ,NY,NX)=RSMN(NZ,NY,NX)+(RSMH(NZ,NY,NX)
     2-RSMN(NZ,NY,NX))*WFNC
      RA(NZ,NY,NX)=RAB(NY,NX)
      VHCPCC=2.496*WVPLT(NZ,NY,NX)*VSTK
      DO 4300 N=1,MY(NZ,NY,NX)
      DO 4300 L=NU(NY,NX),NI(NZ,NY,NX)
      PSIRT(N,L,NZ,NY,NX)=PSIST1(L)
      APSIRT=ABS(PSIRT(N,L,NZ,NY,NX))
      FDMR=FDMPM+0.10*APSIRT/(0.05*APSIRT+2.0)
      CCPOLT=CCPOLR(N,L,NZ,NY,NX)+CZPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)
      OSWT=36.0+840.0*AMAX1(0.0,CCPOLT)
      PSIRO(N,L,NZ,NY,NX)=FDMR/FDMPM*OSMO(NZ,NY,NX)
     2-8.3143*TKS(L,NY,NX)*FDMR*(CCPOLT/OSWT+CSALTR(N,L,NZ,NY,NX))
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)
     2-PSIRO(N,L,NZ,NY,NX))
      UPWTR(N,L,NZ,NY,NX)=0.0
4300  CONTINUE
      ENDIF
C
C     SET CANOPY GROWTH TEMPERATURE FROM SOIL SURFACE
C     OR CANOPY TEMPERATURE DEPENDING ON GROWTH STAGE
C
C     TKG=canopy growth temperature (K)
C     TKS,TKC=soil,canopy surface temperature (K)
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
C     TKG,TKGO=canopy temperature,canopy temperature used in Arrhenius
C        equation (K)
C     TKS,TKSO=soil temperature,soil temperature used in Arrhenius
C        equation (K)
C     OFFST=shift in Arrhenius curve for thermal adaptation (K)
C     TFN3,TFN4=temperature function for canopy,root growth (25 oC =1)
C     8.3143,710.0=gas constant,enthalpy (J mol-1)
C     62500,197500,222500=energy of activn,high,low temp inactivation
C        (J mol-1)
C     PSILZ=minimum daily canopy water potential (MPa)
C
      TKGO=TKG(NZ,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKGO
      STK=710.0*TKGO
      ACTV=1.0+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFN3(NZ,NY,NX)=EXP(25.229-62500/RTK)/ACTV
C     WRITE(*,4438)'TFN3',I,J,NFZ,NX,NY,NZ
C    2,TFN3(NZ,NY,NX),TKG(NZ,NY,NX),OFFST(NZ,NY,NX)
C    3,TKS(NU(NY,NX),NY,NX),TKC(NZ,NY,NX),TKCY
4438  FORMAT(A8,6I4,12E12.4)
      DO 100 L=NU(NY,NX),NI(NZ,NY,NX)
      TKSO=TKS(L,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKSO
      STK=710.0*TKSO
      ACTV=1.0+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFN4(L,NZ,NY,NX)=EXP(25.229-62500/RTK)/ACTV
100   CONTINUE
      PSILZ(NZ,NY,NX)=AMIN1(PSILZ(NZ,NY,NX),PSILT(NZ,NY,NX))
C
C     CHILLING AND HEATING EFFECTS
C
C     TCC=canopy surface temperature (oC)
C     CTC=chilling temperature from PFT file (oC)
C     CHILL=accumulated chilling hours during low temperature events
C        used to limit CO2 fixation in ‘stomate.f’(h)
C     HEAT=accumulated heating hours during high temperature events
C        (eg fire) used to limit CO2 fixation in ‘stomate.f’(h)
C     XNFH=time step for biological activity from ‘wthr.f’ (h t-1)
C
      IF(TCC(NZ,NY,NX).LT.CTC(NZ,NY,NX))THEN
      CHILL(NZ,NY,NX)=AMIN1(24.0,CHILL(NZ,NY,NX)+1.0*XNFH)
      ELSE
      CHILL(NZ,NY,NX)=AMAX1(0.0,CHILL(NZ,NY,NX)-1.0*XNFH)
      ENDIF
      IF(TCC(NZ,NY,NX).GT.60.0)THEN
      HEAT(NZ,NY,NX)=HEAT(NZ,NY,NX)+(TCC(NZ,NY,NX)-60.0)*XNFH
      ELSE 
      HEAT(NZ,NY,NX)=AMAX1(0.0,HEAT(NZ,NY,NX)-0.02*XNFH)
      ENDIF
C
C     NH3 EXCHANGE BETWEEN CANOPY AND ATMOSPHERE FROM NH3
C     CONCENTRATION DIFFERENCES 'CNH3E' (ATMOSPHERE FROM 'READS') AND
C     'CNH3P' (CANOPY), AND FROM STOMATAL + BOUNDARY LAYER RESISTANCE
C
C     SNH3P,SNH3X=NH3 solubility at TCC, 25 oC from ‘hour1.f’ 
C        (g N m-3/(g N m-3)) 
C     TCC=canopy temperature (oC)
C     FDMP=canopy dry matter content (g C g C-1) 
C     FNH3P=canopy structural C:air volume ratio (g NH3-N 
C        g non-structural N g C m-3)
C     ARLFB,ARLFP=branch,canopy leaf area (m2)
C     CNH3P,CNH3E=gaseous NH3 concentration in branch,atmosphere 
C        (g N m-3)
C     CZPOLB,ZPOOLB=nonstructural N concentration,content 
C        (g N g C-1,g N) 
C     RNH3B=NH3 flux between atmosphere and branch (g N t-1)
C     RA,RC=canopy aerodynamic+boundary layer,stomatal resistance 
C        (h m-1)
C     FLAIP=fraction of total canopy+standing dead radiation
C        received by each PFT canopy from ‘hour1.f’
C     XNFH=time step for biological activity from ‘wthr.f’ (h t-1)
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
     3*FLAIP(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)*XNFH
     3*ARLFB(NB,NZ,NY,NX)/ARLFP(NZ,NY,NX),-0.1*ZPOOLB))
      ELSE
      RNH3B(NB,NZ,NY,NX)=0.0
      ENDIF
C     WRITE(*,7777)'RNH3',I,J,NFZ,NZ,NB,RNH3B(NB,NZ,NY,NX)
C    2,RNH3C(NZ,NY,NX),CNH3E(NY,NX),CNH3P,RABC(NY,NX),RC(NZ,NY,NX)
C    2,ARLFB(NB,NZ,NY,NX),ARLFP(NZ,NY,NX),SNH3P
C    4,ZPOOL(NB,NZ,NY,NX),WTLSB(NB,NZ,NY,NX),FRADP(NZ,NY,NX)
7777  FORMAT(A8,5I4,40E24.16)
105   CONTINUE
C
C     ROOT(N=1) AND MYCORRHIZAL(N=2) O2 AND NUTRIENT UPTAKE
C
C     RTDNP=root length density (m m-3)
C     RTVLW=root aqueous volume (m3)
C     THETW=soil water concentration (m3 m-3)
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
C     RELATIVE UPTAKE CAPACITY 'FWSRT' DEPENDS ON ROOT,MYCORRHIZAL 
C     PROTEIN CONTENT RELATIVE TO 5% FOR WHICH ACTIVE UPTAKE
C     PARAMETERS IN THE PFT FILE ARE DEFINED
C
C     CWSRTL,CWSRT=current,maximum protein concentration from
C        ‘startq.f’ (g C g C-1)
C     WSRTL,WTRTL=protein content,mass (g,g C)
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
C     RCO2N=total root respiration unlimited by nonstructural C
C        from ‘grosub.f’ (g C t-1) 
C     FCUP=limitation to root active uptake respiration from CPOOLR
C     CPOOLR=root nonstructural C (g C) 
C
      IF(RCO2N(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FCUP=AMAX1(0.0,AMIN1(1.0,CPOOLR(N,L,NZ,NY,NX)*XNFH
     2/RCO2N(N,L,NZ,NY,NX)))
      ELSE
      FCUP=0.0
      ENDIF
C
C     FEEDBACK CONSTRAINT ON N UPTAKE FROM NON-STRUCTURAL N AND P
C
C     FZUP,FPUP=limitation to root active uptake respiration 
C        from CZPOLR,CPPOLR
C     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration 
C        (g C,N,P g C-1)
C     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition effect on N,P uptake 
C        (g N,P g C-1)
C     UPWTRH=root water uptake (m3 t-1)
C     XNPG=time step for gas fluxes from ‘wthr.f’ (h t-1)
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
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.NFZ.EQ.1)THEN
C     WRITE(*,4414)'FZUP',I,J,NFZ,NX,NY,NZ,L,N,FZUP,FCUP,FPUP
C    2,CCPOLR(N,L,NZ,NY,NX),CZPOLR(N,L,NZ,NY,NX),CPPOLR(N,L,NZ,NY,NX) 
C    3,ZCKI,ZPKI,PCKI,PZKI
C    4,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX)
C    2+CZPOLR(N,L,NZ,NY,NX)/ZCKI)
C    5,CPPOLR(N,L,NZ,NY,NX)/(CPPOLR(N,L,NZ,NY,NX)
C    4+CZPOLR(N,L,NZ,NY,NX)/ZPKI)
C    6,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX)
C    2+CPPOLR(N,L,NZ,NY,NX)/PCKI)
C    7,CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX)
C    4+CPPOLR(N,L,NZ,NY,NX)/PZKI)
4414  FORMAT(A8,8I4,20E12.4)
C     ENDIF
C
C     FACTORS CONSTRAINING O2 AND NUTRIENT UPTAKE AMONG
C     COMPETING ROOT,MYCORRHIZAL AND MICROBIAL POPULATIONS 
C     IN BAND AND NON-BAND SOIL ZONES FROM DEMAND CALCULATED 
C     IN PREVIOUS HOUR
C
C     ROXYY=O2 demand by all microbial,root,mycorrhizal populations 
C        (g O t-1) 
C     ROXYP=O2 demand by each root,mycorrhizal population
C        (g O t-1) 
C     FOXYX=fraction of ROXYY by each root,mycorrhizal population
C     RNH4Y=NH4 demand in non-band by all microbial,root,mycorrhizal
C        populations (g N t-1) 
C     RUNNHP=NH4 demand in non-band by each root,mycorrhizal
C        population (g N t-1)
C     FNH4X=fraction of RNH4Y by each root,mycorrhizal population 
C     RNHBY=NH4 demand in band by all microbial,root, mycorrhizal
C        populations (g N t-1)
C     RUNNBP=NH4 demand in band by each root,mycorrhizal population
C          (g N t-1)
C     FNHBX=fraction of RNHBY by each root,mycorrhizal population 
C     RNO3Y=NO3 demand in non-band by all microbial,root,mycorrhizal
C        populations (g N t-1)
C     RUNNOP=NO3 demand in non-band by each root,mycorrhizal
C        population (g N t-1)
C     FNO3X=fraction of RNO3Y by each root,mycorrhizal populn 
C     RN3BY=NO3 demand in band by all microbial,root,mycorrhizal
C        populations (g N t-1)
C     RUNNXB=NO3 demand in band by each root,mycorrhizal population
C        (g N t-1)
C     FNOBX=fraction of RN3BY by each root,mycorrhizal population 
C     RPO4Y=H2PO4 demand in non-band by all microbial,root,mycorrhizal
C        populations (g P t-1)
C     RUPP2P=H2PO4 demand in non-band by each root,mycorrhizal
C        population (g P t-1)
C     FPO4X=fraction of RPO4Y by each root,mycorrhizal population 
C     RPOBY=H2PO4 demand in band by all microbial,root,mycorrhizal
C        populations (g P t-1)
C     RUPP2B=H2PO4 demand in band by each root,mycorrhizal population
C        (g P t-1)
C     FPOBX=fraction of RPOBY by each root,mycorrhizal population 
C     RP14Y=HPO4 demand in non-band by all microbial,root,mycorrhizal
C        populations (g P t-1)
C     RUPP1P=HPO4 demand in non-band by each root,mycorrhizal 
C        population (g P t-1)
C     FP14X=fraction of RP14Y by each root,mycorrhizal population 
C     RP1BY=HPO4 demand in band by all microbial,root,mycorrhizal 
C        populations (g P t-1)
C     RUPP1B=HPO4 demand in band by each root,mycorrhizal population
C        (g P t-1)
C     FP1BX=fraction of RP1BY by each root,mycorrhizal population 
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
C     ROXYP=O2 demand g O t-1) 
C     RCO2M=respiration unlimited by O2 from ‘grosub.f’ (g C t-1)
C     RTVLW=root or mycorrhizal aqueous volume (m3)
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
C     CO2A1,CO2P1,CO2G1,CO2S1=gaseous,aqueous CO2 in root,soil (g C)
C     OXYA1,OXYP1,OXYG1,OXYS1=gaseous,aqueous O2 in root,soil (g O) 
C     CH4A1,CH4P1,CH4G1,CH4S1=gaseous,aqueous CH4 in root,soil (g C) 
C     Z2OA1,Z2OP1,Z2OG1,Z2OS1=gaseous,aqueous N2O in root,soil (g N) 
C     ZH3A1,ZH3P1,ZH3G1,ZH3S1=gaseous,aqueous NH3 in root,soil (g N) 
C     H2GA1,H2GP1,H2GG1,H2GS1=gaseous,aqueous H2 in root,soil (g H)
C     CCH4S1,CCH4P1=aqueous CH4 concentration in soil,root (g Cm-3) 
C     CN2OS1,CN2OP1=aqueous N2O concentration in soil,root (g N m-3)  
C     CNH3S1,CNH3B1,CNH3P1=aqueous NH3 concentration in soil 
C        non-band,band,root (g N m-3) 
C     CH2GS1,CH2GP1=aqueous H2 concentration in soil,root (g H m-3)
C     RTVLWA,RTVLWB=root aqueous volume in non-band,band (m3)
C     UPMXP=O2 demand per plant (g O t-1) 
C     ROXYF=net diffusive-convective O2 gas flux 
C        from previous time step (g O t-1) 
C     RCO2F=net diffusive-convective CO2 gas flux 
C        from previous time step (g C t-1) 
C     ROXYL=net diffusive-convective O2 aqueous flux 
C        from previous time step (g O t-1) 
C     FOXYX=fraction of total O2 demand by each root,mycorrhizal
C        population
C     FPQ=PFT fraction of biome root mass
C     XNPG=time step of gas flux calculation (h t-1)
C
      CO2A1=AMAX1(0.0,CO2A(N,L,NZ,NY,NX))
      CO2P1=AMAX1(0.0,CO2P(N,L,NZ,NY,NX))
      CO2G1=AMAX1(0.0,CO2G(L,NY,NX)*FPQ(N,L,NZ))
      CO2S1=AMAX1(0.0,CO2S(L,NY,NX)*FPQ(N,L,NZ))
      OXYA1=AMAX1(0.0,OXYA(N,L,NZ,NY,NX))
      OXYP1=AMAX1(0.0,OXYP(N,L,NZ,NY,NX))
      OXYG1=AMAX1(0.0,OXYG(L,NY,NX)*FOXYX)
      OXYS1=AMAX1(0.0,OXYS(L,NY,NX)*FOXYX)
      CH4A1=CH4A(N,L,NZ,NY,NX)
      CH4P1=CH4P(N,L,NZ,NY,NX)
      CH4S1=CH4S(L,NY,NX)*FPQ(N,L,NZ)
      Z2OA1=Z2OA(N,L,NZ,NY,NX)
      Z2OP1=Z2OP(N,L,NZ,NY,NX)
      Z2OS1=Z2OS(L,NY,NX)*FPQ(N,L,NZ)
      ZH3A1=ZH3A(N,L,NZ,NY,NX)
      ZH3P1=ZH3P(N,L,NZ,NY,NX)
      ZH3S1=ZNH3S(L,NY,NX)*FPQ(N,L,NZ)
      ZH3B1=ZNH3B(L,NY,NX)*FPQ(N,L,NZ)
      H2GA1=H2GA(N,L,NZ,NY,NX)
      H2GP1=H2GP(N,L,NZ,NY,NX)
      H2GS1=H2GS(L,NY,NX)*FPQ(N,L,NZ)
      RTVLWA=RTVLW(N,L,NZ,NY,NX)*VLNH4(L,NY,NX)
      RTVLWB=RTVLW(N,L,NZ,NY,NX)*VLNHB(L,NY,NX)
      UPMXP=ROXYP(N,L,NZ,NY,NX)*XNPG/PP(NZ,NY,NX)
      ROXYFX=ROXYF(L,NY,NX)*FOXYX*XNPG
      RCO2FX=RCO2F(L,NY,NX)*FOXYX*XNPG
      ROXYLX=ROXYL(L,NY,NX)*FOXYX*XNPG
C
C     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
C
C     *SGL(=diffusivity from ‘hour1.f’ (m2 h-1)
C     *SGL1=diffusivity (m2 t-1)
C     PORTX=tortuosity effect of root porosity on diffusivity
C     PORT=root porosity from PFT file (m3 m-3)
C     XNPGX=time step for gas transfer calculations from ‘wthr.f’ 
C     soil gas code:CG=CO2g,OG=O2g,CH=CH4g,Z2=N2Og,ZH=NH3g,HG=H2g
C                   CL=CO2s,OL=O2s,CQ=CH4s,ZV=N2Os,ZN=NH3s,HL=H2s
C                   *g=gaseous,*s=aqueous     
C
      PORTX=PORT(N,NZ,NY,NX)**1.33
      CGSGL1=CGSGL(L,NY,NX)*XNPGX*PORTX
      OGSGL1=OGSGL(L,NY,NX)*XNPGX*PORTX
      CHSGL1=CHSGL(L,NY,NX)*XNPGX*PORTX
      Z2SGL1=Z2SGL(L,NY,NX)*XNPGX*PORTX
      ZHSGL1=ZHSGL(L,NY,NX)*XNPGX*PORTX
      HGSGL1=HGSGL(L,NY,NX)*XNPGX*PORTX
      CLSGL1=CLSGL(L,NY,NX)*XNPGX*FOXYX
      OLSGL1=OLSGL(L,NY,NX)*XNPGX*FOXYX
      CQSGL1=CQSGL(L,NY,NX)*XNPGX*FOXYX
      ZVSGL1=ZVSGL(L,NY,NX)*XNPGX*FOXYX
      ZNSGL1=ZNSGL(L,NY,NX)*XNPGX*FOXYX
      HLSGL1=HLSGL(L,NY,NX)*XNPGX*FOXYX
      OLSGLP=OLSGL(L,NY,NX)*XNPGX
      RCODF1=0.0
      ROXDF1=0.0
      RCHDF1=0.0
      RN2DF1=0.0
      RNHDF1=0.0
      RHGDF1=0.0
C
C     SALTS
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C     Z*=soil salt content (mol)
C     Z*R=root salt content (mol)
C     *SGL=salt diffusivity (m2 t-1)
C     salt code:AL=Al,FE=Fe,Ca=Ca,MG=Mg,NA=Na,KA=K,SO=SO4,CL=Cl
C
      IF(ISALTG.NE.0)THEN
      ZALS1=ZAL(L,NY,NX)*FPQ(N,L,NZ)
      ZALP1=ZALR(N,L,NZ,NY,NX)
      ALSGL1=ALSGL(L,NY,NX)*XNPHX
      ZFES1=ZFE(L,NY,NX)*FPQ(N,L,NZ)
      ZFEP1=ZFER(N,L,NZ,NY,NX)
      FESGL1=FESGL(L,NY,NX)*XNPHX
      ZCAS1=ZCA(L,NY,NX)*FPQ(N,L,NZ)
      ZCAP1=ZCAR(N,L,NZ,NY,NX)
      CASGL1=CASGL(L,NY,NX)*XNPHX
      ZMGS1=ZMG(L,NY,NX)*FPQ(N,L,NZ)
      ZMGP1=ZMGR(N,L,NZ,NY,NX)
      GMSGL1=GMSGL(L,NY,NX)*XNPHX
      ZNAS1=ZNA(L,NY,NX)*FPQ(N,L,NZ)
      ZNAP1=ZNAR(N,L,NZ,NY,NX)
      ANSGL1=ANSGL(L,NY,NX)*XNPHX
      ZKAS1=ZKA(L,NY,NX)*FPQ(N,L,NZ)
      ZKAP1=ZKAR(N,L,NZ,NY,NX)
      AKSGL1=AKSGL(L,NY,NX)*XNPHX
      ZSOS1=ZSO4(L,NY,NX)*FPQ(N,L,NZ)
      ZSOP1=ZSOR(N,L,NZ,NY,NX)
      SOSGL1=SOSGL(L,NY,NX)*XNPHX
      ZCLS1=ZCL(L,NY,NX)*FPQ(N,L,NZ)
      ZCLP1=ZCLR(N,L,NZ,NY,NX)
      CLSGL1=CLSGL(L,NY,NX)*XNPHX
      ENDIF
C
C     ROOT CONDUCTANCE TO GAS TRANSFER DETERMINED BY ROOT DIMENSIONS
C
C     WTRTS=total root,mycorrhizal mass (g C)
C     FRTDPX=fraction of each soil layer with primary root
C     RTCR1,RTCR2,RTCRA=cross-sectional area/length of
C        primary,secondary,total root,mycorrhizal system (m)
C     RTN1,RTNL=number of root,mycorrhizal primary,secondary axes
C     RRAD1,RRAD2=primary,secondary root radius (m)
C     DPTHZ=depth of primary root from surface (m)
C     RTLGA=average secondary root length (m)
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
C     RTLGP=root,mycorrhizal length per plant (m)
C     IDAY(1,=emergence date
C     RTARR=root surface area/radius for uptake (m)
C     RRADP=path length for radial diffusion within root
C     OLSGLP=aqueous diffusivity of O2 (m2 t-1)  
C     DIFOP=aqueous diffusivity of O2 within root (m3 t-1)
C     RTVLW=root,mycorrhizal aqueous volume (m3)
C     S*L=solubility of gas in water from ‘hour1.f’ (g m-3/(g m-3)):
C         CO2=CO2,OXY=O2,CH4=CH4,N2O=N2O,NH3=NH3,H2G=H2
C     DF*A=root-atmosphere gas conductance (m3 t-1)
C     root gas code:CO=CO2a,OX=O2a,CH=CH4a,N2=N2Oa,NH=NH3a,HG=H2a
C                   CL=CO2s,OL=O2s,CQ=CH4s,ZV=N2Os,ZN=NH3s,HL=H2s
C                   *a=gaseous,*s=aqueous     
C     DFGP=rate constant for equilibration of gas concentration 
C        in gaseous-aqueous phases (t-1)
C     RCO2PX=root CO2 gas flux (g C t-1)
C     RCO2A=root CO2 flux from grosub.f (g C t-1)
C     XNPGX,XNPG=time step for gas transfer calculations from ‘wthr.f’ 
C        (h t-1)
C
      IF(N.EQ.1.AND.IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).GT.0
     2.AND.RTLGP(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RRADS=LOG((FILM(NPH,L,NY,NX)+RRADL(N,L))/RRADL(N,L))
      RTARRX=RTARR(N,L)/RRADS 
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
      DFGP=AMIN1(1.0,PORT(N,NZ,NY,NX)*XNPT)
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
C              NHB=NH3 band,H2G=H2
C     VOLWMM,VOLPMM=soil micropore water,air volume (m3)
C     FOXYX=root fraction of total O2 demand from previous hour
C     FPQ=PFT fraction of biome root mass
C     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
C     VOLX=soil volume excluding rock,macropores (m3)
C     THETW1=soil water concentration (m3 m-3)
C     TORT=soil tortuosity from ‘watsub.f’
C     FILM=soil water film thickness from ‘watsub.f’ (m)
C     RRADL=root radius (m)
C     RRADS=path length effect on radial diffusion from soil to root
C     RTARR=root surface area/radius for uptake,diffusivity (m)
C     DIF*=aqueous diffusivity from soil to root (m3 t-1)
C        :OL=O2,CL=CH4,ZL=N2O,NL=NH3 non-band,NB=NH4 band,HL=H2
C     *G1=soil gaseous content (g) 
C     C*G1=soil gaseous concentration (g m-3) 
C     VOLW*,VOLP*=VOLWMM,VOLPMM*gas solubility (m3)
C
      VOLWMO=VOLWM(M,L,NY,NX)*FOXYX
      VOLWMM=VOLWM(M,L,NY,NX)*FPQ(N,L,NZ)
      VOLPMM=VOLPM(M,L,NY,NX)*FPQ(N,L,NZ)
      VOLWSP=RTVLW(N,L,NZ,NY,NX)+VOLWMM
      VOLWMA=VOLWMM*VLNH4(L,NY,NX)
      VOLWMB=VOLWMM*VLNHB(L,NY,NX)
      VOLWSA=RTVLWA+VOLWMA
      VOLWSB=RTVLWB+VOLWMB
      IF(VOLY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      THETW1=AMAX1(0.0,VOLWM(M,L,NY,NX)/VOLY(L,NY,NX))
      ELSE
      THETW1=0.0
      ENDIF
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
C     C*S1=soil aqueous concentration non-band (g m-3) 
C     C*B1=soil aqueous concentration band (g m-3)
C     C*A1=root gaseous concentration (g m-3) 
C     C*P1=root aqueous concentration (g m-3)
C     ROXYLX=soil net O2 aqueous flux (g O t-1) 
C     VOLWMM=micropore water volume (m3)
C     RTVLW,RTVLP=root aqueous,gaseous volume (m3)
C     RMF*=soil convective solute flux (g t-1)
C        :COS=CO2,OXS=O2,CHS=CH4,N2S=N2O,NHS=NH3 non-band
C        :NHB=NH3 band,HGS=H2
C     UPWTRH=root water uptake (m3 t-1) 
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
C     DIFOL=O2 aqueous diffusivity from soil to root (m3 t-1)
C     UPWTRH=water uptake (m3 t-1)
C     DIFOP=aqueous diffusivity of O2 within root (m3 t-1)
C     COXYS1,COXYP1=soil,root aqueous O2 concentration (g O m-3)
C     UPMXP=O2 demand per plant (g O t-1)
C     RUPOXR=root O2 uptake per plant (g O t-1)
C     COXYR=aqueous O2 concentration at root surface (g O m-3)
C     RDFOXS,RDFOXP=aqueous O2 diffusion per plant in aqueous phase
C         of soil,root (g O t-1)
C
      X=(DIFOL+UPWTRH)*COXYS1+DIFOP*COXYP1
      IF(X.GT.ZERO.AND.OXYS1.GT.ZEROP(NZ,NY,NX))THEN
      B=-UPMXP-DIFOX*OXKM-X
      C=X*UPMXP
      RUPOXR=(-B-SQRT(B*B-4.0*C))/2.0
      COXYR=(X-RUPOXR)/DIFOX
      RDFOXS=RMFOXS+DIFOL*(COXYS1-COXYR)
      RDFOXP=DIFOP*(COXYP1-COXYR)
      ELSE
      X=DIFOP*COXYP1
      IF(X.GT.ZERO.AND.OXYP1.GT.ZEROP(NZ,NY,NX))THEN
      B=-UPMXP-DIFOP*OXKM-X
      C=X*UPMXP
      RUPOXR=(-B-SQRT(B*B-4.0*C))/2.0
      COXYR=(X-RUPOXR)/DIFOP
      RDFOXS=0.0
      RDFOXP=AMIN1(AMAX1(0.0,FOXYX*OXYP1/PP(NZ,NY,NX))
     2,DIFOP*(COXYP1-COXYR))
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
C     RUPOSX,RUPOPX=root aqueous O2 uptake from soil,root (g O t-1)
C     PP=PFT population
C     RDFCOS,RCO2SX=aqueous CO2 soil-root diffusion,uptake (g C t-1)
C     RDFCHS,RUPCSX=aqueous CH4 soil-root diffusion,uptake (g C t-1)
C     RDFN2S,RUPZSX=aqueous N2O soil-root diffusion,uptake (g N t-1)
C     RDFN3S,RUPNSX=aqueous NH3 soil-root diffusion,uptake:non-band 
C        (g N t-1)
C     RDFN3B,RUPNBX=aqueous NH3 soil-root diffusion,uptake:band
C        (g N t-1)
C     RDFHGS,RUPHGX=aqueous H2 soil-root diffusion,uptake (g H t-1)
C     RMF*=soil convective solute flux (g t-1)
C     DIF*=aqueous diffusivity from soil to root (m3 t-1)
C     C*S1=soil aqueous gas concentration non-band (g m-3) 
C     C*B1=soil aqueous gas concentration band (g m-3)
C     C*P1=root aqueous gas concentration (g m-3)
C
C     IF(NZ.EQ.4)THEN
C     WRITE(*,5555)'COXYR',I,J,NFZ,M,MX,NX,NY,NZ,L,N,COXYR,RUPOXR
C    2,RMFOXS,RDFOXS,RDFOXP,COXYS1,COXYP1,FOXYX,PATH(N,L)
C    3,DIFOL,DIFOP,THETM,OLSGL1,UPWTRH,RTARR(N,L)
C    5,RTARRX,UPMXP,THETW(L,NY,NX),OXYS1,OXYS(L,NY,NX),OXYP1
C    3,OXYP(N,L,NZ,NY,NX),ROXYY(L,NY,NX) 
C    2,UPMXP,DIFOX,THETW1,THETM,RRADS,FPQ(N,L,NZ) 
C    4,RUPOXS(N,L,NZ,NY,NX),RUPOXP(N,L,NZ,NY,NX)
C    5,COXYE(NY,NX),SOXYL(L,NY,NX),RTARRX,RTARR(N,L),RRADS
C    6,RRADP(N,NZ,NY,NX),OLSGL1,OLSGP,THETM,FRTDPX(L,NZ)
C    7,RTLGP(N,L,NZ,NY,NX) 
5555  FORMAT(A8,10I4,50E12.4)
C     ENDIF
      RUPOSX=RDFOXS*PP(NZ,NY,NX)
      RUPOPX=RDFOXP*PP(NZ,NY,NX)
      RDFCOS=RMFCOS+DIFCL*(CCO2S1-CCO2P1)
      RDXCOS=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,CO2S1)
     2-VOLWMM*AMAX1(0.0,CO2P1))/VOLWSP
      IF(RDFCOS.GT.0.0)THEN
      RCO2SX=AMIN1(AMAX1(0.0,RDXCOS),RDFCOS*PP(NZ,NY,NX))
      ELSE
      RCO2SX=AMAX1(AMIN1(0.0,RDXCOS),RDFCOS*PP(NZ,NY,NX))
      ENDIF
      IF(N.EQ.1)THEN
      RDFCHS=RMFCHS+DIFCL*(CCH4S1-CCH4P1)
      RDXCHS=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,CH4S1)
     2-VOLWMM*AMAX1(0.0,CH4P1))/VOLWSP
      IF(RDFCHS.GT.0.0)THEN
      RUPCSX=AMIN1(AMAX1(0.0,RDXCHS),RDFCHS*PP(NZ,NY,NX))
      ELSE
      RUPCSX=AMAX1(AMIN1(0.0,RDXCHS),RDFCHS*PP(NZ,NY,NX))
      ENDIF
      RDFN2S=RMFN2S+DIFZL*(CN2OS1-CN2OP1)
      RDXN2S=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,Z2OS1)
     2-VOLWMM*AMAX1(0.0,Z2OP1))/VOLWSP
      IF(RDFN2S.GT.0.0)THEN
      RUPZSX=AMIN1(AMAX1(0.0,RDXN2S),RDFN2S*PP(NZ,NY,NX))
      ELSE
      RUPZSX=AMAX1(AMIN1(0.0,RDXN2S),RDFN2S*PP(NZ,NY,NX))
      ENDIF
      RDFN3S=RMFN3S+DIFNL*(CNH3S1-CNH3P1)
      IF(VOLWSA.GT.ZEROP(NZ,NY,NX))THEN
      ZH3PA=ZH3P1*VLNH4(L,NY,NX)
      RDXNHS=(RTVLWA*AMAX1(0.0,ZH3S1)
     2-VOLWMA*AMAX1(0.0,ZH3PA))/VOLWSA
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
      RDXNHB=(RTVLWB*AMAX1(0.0,ZH3B1)
     2-VOLWMB*AMAX1(0.0,ZH3PB))/VOLWSB
      ELSE
      RDXNHB=0.0
      ENDIF
      IF(RDFN3B.GT.0.0)THEN
      RUPNBX=AMIN1(AMAX1(0.0,RDXNHB),RDFN3B*PP(NZ,NY,NX))
      ELSE
      RUPNBX=AMAX1(AMIN1(0.0,RDXNHB),RDFN3B*PP(NZ,NY,NX))
      ENDIF
      RDFHGS=RMFHGS+DIFHL*(CH2GS1-CH2GP1)
      RDXHGS=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,H2GS1)
     2-VOLWMM*AMAX1(0.0,H2GP1))/VOLWSP
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
C     THETPM,THETX=air-filled porosity,minimum THETPM (m3 m-3)
C     R*DFQ=soil gas exchange between gaseous-aqueous phases (g t-1)
C        gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,HG=H2
C     DFGS=rate constant for gaseous-aqueous exchange from ‘watsub.f’
C        (t-1)
C     FPQ=PFT fraction of biome root mass
C     CO2G1,CO2S1=gaseous,aqueous CO2 in soil (g C)
C     OXYG1,OXYS1=gaseous,aqueous O2 in soil (g O) 
C     CH4G1,CH4S1=gaseous,aqueous CH4 in soil (g C) 
C     Z2OG1,Z2OS1=gaseous,aqueous N2O in soil (g N) 
C     ZH3G1,ZH3S1,ZH3B1=gaseous,aqueous NH3 in soil non-band,band (g N)
C     H2GG1,H2GS1=gaseous,aqueous H2 in soil (g H)
C     RUPOSX=root aqueous O2 uptake (g O t-1)
C     ROXYLX=root net O2 aqueous flux (g O t-1) 
C     RCO2SX=soil aqueous CO2 exchange (g C t-1) 
C     RUPCSX=root aqueous CH4 uptake (g C t-1)
C     RUPZSX=root aqueous N2O uptake (g N t-1)
C     RUPNSX=root aqueous NH3 uptake non-band (g N t-1) 
C     RUPNBX=root aqueous NH3 uptake band (g N t-1)
C     RUPHGX=root aqueous H2 uptake (g H t-1)
C     VOLWMM,VOLPMM=soil micropore water,air volume (m3)
C     VOLW*=VOLWMM*gas solubility (m3)
C
      IF(THETPM(M,L,NY,NX).GT.THETX)THEN
      DFGSP=FPQ(N,L,NZ)*DFGS(M,L,NY,NX)
      RCODFQ=DFGSP*(AMAX1(0.0,CO2G1)*VOLWCO
     2-(AMAX1(0.0,CO2S1)-RCO2SX)*VOLPMM)/(VOLWCO+VOLPMM)
      RUPOST=RUPOSX-ROXYLX
      ROXDFQ=DFGSP*(AMAX1(0.0,OXYG1)*VOLWOX
     2-(AMAX1(0.0,OXYS1)-RUPOST)*VOLPMM)/(VOLWOX+VOLPMM)
      IF(N.EQ.1)THEN
      RCHDFQ=DFGSP*(AMAX1(0.0,CH4G1)*VOLWCH
     2-(AMAX1(0.0,CH4S1)-RUPCSX)*VOLPMM)/(VOLWCH+VOLPMM)
      RN2DFQ=DFGSP*(AMAX1(0.0,Z2OG1)*VOLWN2
     2-(AMAX1(0.0,Z2OS1)-RUPZSX)*VOLPMM)/(VOLWN2+VOLPMM)
      IF(VOLWNH+VOLPNH.GT.ZEROP(NZ,NY,NZ))THEN
      ZH3GA=ZH3G1*VLNH4(L,NY,NX)
      RNHDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX
     2,DFGSP*(AMAX1(0.0,ZH3GA)*VOLWNH
     3-(AMAX1(0.0,ZH3S1)-RUPNSX)*VOLPNH)/(VOLWNH+VOLPNH)))
      ELSE
      RNHDFQ=0.0
      ENDIF
      IF(VOLWNB+VOLPNB.GT.ZEROP(NZ,NY,NZ))THEN
      ZH3GB=ZH3G1*VLNHB(L,NY,NX)
      RNBDFQ=AMIN1(RUPNSX,AMAX1(-RUPNSX
     2,DFGSP*(AMAX1(0.0,ZH3GB)*VOLWNB
     3-(AMAX1(0.0,ZH3B1)-RUPNBX)*VOLPNB)/(VOLWNB+VOLPNB)))
      ELSE
      RNBDFQ=0.0
      ENDIF
      RHGDFQ=DFGSP*(AMAX1(0.0,H2GG1)*VOLWHG
     2-(AMAX1(0.0,H2GS1)-RUPHGX)*VOLPMM)/(VOLWHG+VOLPMM)
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
C     UPDATE SOIL GASEOUS, AQUEOUS GAS CONTENTS AND CONCENTRATIONS
C     FROM GASEOUS-AQUEOUS EXCHANGE, SOIL GAS TRANSFERS
C
C     *G1,*S1=soil gaseous,aqueous gas content (g)
C     R*DFQ=soil gas exchange between gaseous-aqueous phases (g t-1)
C     R*FX=soil net diffusive-convective gaseous flux (g t-1) 
C     R*SX=root aqueous gas uptake from soil (g t-1)
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
C     R*DF1=root gas exchange between gaseous-aqueous phases (g t-1)
C     R*FL1=root diffusive-convective gas exchange with atmosphere 
C        (g t-1)
C        gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
C     CO2A1,CO2P1=gaseous,aqueous CO2 in root (g C)
C     OXYA1,OXYP1=gaseous,aqueous O2 in root (g O) 
C     CH4A1,CH4P1=gaseous,aqueous CH4 in root (g C) 
C     Z2OA1,Z2OP1=gaseous,aqueous N2O in root (g N) 
C     ZH3A1,ZH3P1=gaseous,aqueous NH3 in root (g N)
C     H2GA1,H2GP1=gaseous,aqueous H2 in root (g H)
C     RTVLW,RTVLP=root aqueous,gaseous volume (m3)
C     VOLW*=RTVLW*gas solubility (m3)
C     C*E,C*A1=atmosphere,root gas concentration (g m-3)
C     DF*A=root-atmosphere gas conductance (m3 t-1)
C
      IF(RTVLP(N,L,NZ,NY,NX).GT.ZERO)THEN
      CCO2A1=AMAX1(0.0,(CO2A1-RCODF1)/RTVLP(N,L,NZ,NY,NX))
      COXYA1=AMAX1(0.0,(OXYA1-ROXDF1)/RTVLP(N,L,NZ,NY,NX))
      CCH4A1=AMAX1(0.0,(CH4A1-RCHDF1)/RTVLP(N,L,NZ,NY,NX))
      CZ2OA1=AMAX1(0.0,(Z2OA1-RN2DF1)/RTVLP(N,L,NZ,NY,NX))
      CNH3A1=AMAX1(0.0,(ZH3A1-RNHDF1)/RTVLP(N,L,NZ,NY,NX))
      CH2GA1=AMAX1(0.0,(H2GA1-RHGDF1)/RTVLP(N,L,NZ,NY,NX))
      ELSE
      CCO2A1=0.0
      COXYA1=0.0
      CCH4A1=0.0
      CZ2OA1=0.0
      CNH3A1=0.0
      CH2GA1=0.0
      ENDIF
      RCOFL1=AMIN1(DFCOA,RTVLP(N,L,NZ,NY,NX))*(CCO2E(NY,NX)-CCO2A1)
      ROXFL1=AMIN1(DFOXA,RTVLP(N,L,NZ,NY,NX))*(COXYE(NY,NX)-COXYA1)
      RCHFL1=AMIN1(DFCHA,RTVLP(N,L,NZ,NY,NX))*(CCH4E(NY,NX)-CCH4A1)
      RN2FL1=AMIN1(DFN2A,RTVLP(N,L,NZ,NY,NX))*(CZ2OE(NY,NX)-CZ2OA1)
      RNHFL1=AMIN1(DFNHA,RTVLP(N,L,NZ,NY,NX))*(CNH3E(NY,NX)-CNH3A1)
      RHGFL1=AMIN1(DFHGA,RTVLP(N,L,NZ,NY,NX))*(CH2GE(NY,NX)-CH2GA1)
      CO2PX=CO2P1+RCO2PX
      RCODF1=AMAX1(-CO2PX,DFGP*(AMAX1(0.0,CO2A1)*VOLWCA 
     2-CO2PX*RTVLP(N,L,NZ,NY,NX))/(VOLWCA+RTVLP(N,L,NZ,NY,NX)))
      OXYPX=OXYP1-RUPOPX
      ROXDF1=AMAX1(-OXYPX,DFGP*(AMAX1(0.0,OXYA1)*VOLWOA
     2-OXYPX*RTVLP(N,L,NZ,NY,NX))/(VOLWOA+RTVLP(N,L,NZ,NY,NX)))
      ROXDF1=AMIN1(ROXDF1,ROXFL1+OXYA1)
      CH4PX=CH4P1+RUPCSX
      RCHDF1=AMAX1(-CH4PX,DFGP*(AMAX1(0.0,CH4A1)*VOLWC4 
     2-CH4PX*RTVLP(N,L,NZ,NY,NX))/(VOLWC4+RTVLP(N,L,NZ,NY,NX)))
      Z2OPX=Z2OP1+RUPZSX
      RN2DF1=AMAX1(-Z2OPX,DFGP*(AMAX1(0.0,Z2OA1)*VOLWZA 
     2-Z2OPX*RTVLP(N,L,NZ,NY,NX))/(VOLWZA+RTVLP(N,L,NZ,NY,NX)))
      ZH3PX=ZH3P1+RUPNTX
      RNHDF1=AMAX1(-ZH3PX,DFGP*(AMAX1(0.0,ZH3A1)*VOLWNA 
     2-ZH3PX*RTVLP(N,L,NZ,NY,NX))/(VOLWNA+RTVLP(N,L,NZ,NY,NX)))
      H2GPX=H2GP1+RUPHGX 
      RHGDF1=AMAX1(-H2GPX,DFGP*(AMAX1(0.0,H2GA1)*VOLWH2 
     2-H2GPX*RTVLP(N,L,NZ,NY,NX))/(VOLWH2+RTVLP(N,L,NZ,NY,NX)))
      ELSE
      RCOFL1=0.0
      ROXFL1=0.0
      RCHFL1=0.0
      RN2FL1=0.0
      RNHFL1=0.0
      RHGFL1=0.0
      RCODF1=0.0
      ROXDF1=0.0
      RCHDF1=0.0
      RN2DF1=0.0
      RNHDF1=0.0
      RHGDF1=0.0
      ENDIF
C
C     UPDATE ROOT AQUEOUS, GASEOUS GAS CONTENTS AND CONCENTRATIONS
C     FOR ROOT AQUEOUS-GASEOUS, GASEOUS-ATMOSPHERE EXCHANGES
C
C     *A1,*P1=root gaseous,aqueous gas content (g) 
C     R*DF1=root gas exchange between gaseous-aqueous phases (g t-1)
C     R*FL1=root net diffusive-convective gaseous flux (g t-1)
C     R*X=root aqueous gas uptake from root (g t-1)
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
C     ACCUMULATE SOIL-ROOT GAS EXCHANGE FOR ALL M
C
C     RCO2S=soil-root CO2 exchange (g C t-1)
C     RUPOXS=soil-root O2 exchange (g O t-1)
C     RUPCHS=soil-root CH4 exchange (g C t-1)
C     RUPN2S=soil-root N2O exchange (g N t-1)
C     RUPN3S=soil-root NH3 exchange non-band (g N t-1)
C     RUPN3B=soil-root NH3 exchange band (g N t-1)
C     RUPHGS=soil-root H2 exchange (g H t-1)
C     R*X=root gas uptake from soil (g t-1)
C
      RCO2S(N,L,NZ,NY,NX)=RCO2S(N,L,NZ,NY,NX)+RCO2SX
      RUPOXS(N,L,NZ,NY,NX)=RUPOXS(N,L,NZ,NY,NX)+RUPOSX
      RUPCHS(N,L,NZ,NY,NX)=RUPCHS(N,L,NZ,NY,NX)+RUPCSX
      RUPN2S(N,L,NZ,NY,NX)=RUPN2S(N,L,NZ,NY,NX)+RUPZSX
      RUPN3S(N,L,NZ,NY,NX)=RUPN3S(N,L,NZ,NY,NX)+RUPNSX
      RUPN3B(N,L,NZ,NY,NX)=RUPN3B(N,L,NZ,NY,NX)+RUPNBX
      RUPHGS(N,L,NZ,NY,NX)=RUPHGS(N,L,NZ,NY,NX)+RUPHGX
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.1)THEN
C     WRITE(*,5547)'RUPOXS',I,J,NFZ,NX,NY,NZ,L,N,RUPOXS(N,L,NZ,NY,NX)
C    2,RUPOSX 
C     WRITE(*,5547)'RCO2S',I,J,NFZ,NX,NY,NZ,L,N,RCO2S(N,L,NZ,NY,NX)
C    2,RCO2SX,RDXCOS,RDFCOS,RMFCOS,DIFCL,CCO2S1,CCO2P1,VOLWSP,VOLWMM
C     WRITE(*,5547)'RCH4S',I,J,NFZ,NX,NY,NZ,L,N,RUPCHS(N,L,NZ,NY,NX)
C    2,RUPCSX,RDXCHS,RDFCHS,RMFCHS,DIFCL,CCH4S1,CCH4P1,VOLWSP,VOLWMM
C     WRITE(*,5547)'RUPNBX',I,J,NFZ,NX,NY,NZ,L,N,RUPNBX,RDXNHB,RDFN3B
C    2,RTVLWB,ZH3B1,VOLWMB,ZH3PB,VOLWSB,RUPN3B(N,L,NZ,NY,NX)
C    3,RNBDFQ,RNHDF1,RUPNSX
C     WRITE(*,5547)'RUPN2S',I,J,NFZ,NX,NY,NZ,L,N,RUPN2S(N,L,NZ,NY,NX)
C    2,RUPZSX,RDXN2S,RDFN2S,RTVLW(N,L,NZ,NY,NX),Z2OS1,Z2OP1,VOLWMM
C    2,RMFN2S,DIFZL,CN2OS1,CN2OP1,PP(NZ,NY,NX)
C     WRITE(*,5547)'RUPN3S',I,J,NFZ,NX,NY,NZ,L,N,RUPN3S(N,L,NZ,NY,NX)
C    2,RUPNSX,RDXN3S,RDFN3S,RTVLW(N,L,NZ,NY,NX),ZH3S1,ZH3P1,VOLWMM
C    2,RMFN3S,DIFZL,CNH3S1,CNH3P1,PP(NZ,NY,NX),THETM,ZVSGL1,RTARRX
C    3,VOLWM(M,L,NY,NX),VOLY(L,NY,NX) 
5547  FORMAT(A8,8I4,30F16.8)
C     ENDIF
C
C     ACCUMULATE ROOT-ATMOSPHERE GAS EXCHANGE FOR ALL M
C
C     R*DFA,R*DF1=root aqueous-gaseous CO2 exchange (g t-1)
C     R*FLA,R*FL1=root gaseous-atmosphere CO2 exchange (g t-1)
C        gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,HG=H2
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
C     ACCUMULATE SOIL-ROOT GAS EXCHANGE FOR ALL M
C
C     RCO2P=root CO2 emission into root (g C t-1)
C     RCO2PX=root aqueous CO2 exchange with root (g C t-1) 
C     RCO2SX=root aqueous CO2 exchange with soil (g C t-1) 
C     RUPOSX,RUPOPX=root aqueous O2 uptake from soil,root (g O t-1)
C     ROXSK=total soil O2 uptake from soil by all microbial,root
C        populations used in ‘trnsfr.f’ (g O t-1)
C
      RCO2P(N,L,NZ,NY,NX)=RCO2P(N,L,NZ,NY,NX)+RCO2PX+RCO2SX
      RUPOXP(N,L,NZ,NY,NX)=RUPOXP(N,L,NZ,NY,NX)+RUPOPX
      ROXSK(M,L,NY,NX)=ROXSK(M,L,NY,NX)+RUPOSX
C     IF(J.EQ.15.AND.NFZ.EQ.NFH.AND.NZ.EQ.2.AND.L.EQ.7)THEN
C     WRITE(*,5566)'CO2P1',I,J,NFZ,M,MX,NZ,L,N,RCO2S(N,L,NZ,NY,NX) 
C    2,RCO2SX,RDFCOS,RDXCOS,RMFCOS,DIFCL,DFGS(M,L,NY,NX) 
C    3,CCO2S1,CCO2P1,RTVLW(N,L,NZ,NY,NX),CO2S1,VOLWMM
C    4,CO2P1,VOLWSP,PP(NZ,NY,NX),FPQ(N,L,NZ),RCODF1,RCO2PX
C    5,CO2PX,RTVLP(N,L,NZ,NY,NX),DFGP,VOLWCA,CO2A1
C    6,PORT(N,NZ,NY,NX),RCO2FX
C    7,RCODFQ,DFGSP,FOXYX,CO2G1,VOLWCO,VOLPMM,THETPM(M,L,NY,NX),THETX
C    8,ROXYP(N,L,NZ,NY,NX),ROXYY(L,NY,NX) 
C     WRITE(*,5566)'OXYP1',I,J,NFZ,M,MX,NZ,L,N 
C    2,OXYA1,OXYP1,OXYG1,OXYS1,RUPOSX,RUPOPX,ROXDF1,ROXFL1,OXYPX
C    3,DFGP,VOLWOA,RTVLP(N,L,NZ,NY,NX),ROXDFQ 
C    3,FOXYX,DFGS(M,L,NY,NX),DFGP,ROXYFX,ROXYLX 
C    3,COXYS1,COXYP1,COXYR,ROXSK(M,L,NY,NX),XS,XR
C    4,OXYPY,VOLWOA,RTVLP(N,L,NZ,NY,NX),RTVLW(N,L,NZ,NY,NX)
C    5,DFOXA,COXYE(NY,NX),COXYA1,RUPOXP(N,L,NZ,NY,NX)
C    6,RUPOXS(N,L,NZ,NY,NX),ROXYP(N,L,NZ,NY,NX),THETPM(M,L,NY,NX)
C    7,ROXSK(M,L,NY,NX),ROXDFA(N,L,NZ,NY,NX),ROXFLA(N,L,NZ,NY,NX)
C    8,FRTDPX(L,NZ)
C     WRITE(*,5566)'CH4S1',I,J,NFZ,M,MX,NZ,L,N,CH4S1,RCHDFQ,RUPCSX
C    2,RDFCHS,RMFCHS,DIFCL,CCH4S1,CCH4P1,CH4P1,CH4PX,RCHDF1,RCHFL1 
C    3,DFCHA,RTVLP(N,L,NZ,NY,NX),CCH4E(NY,NX),CCH4A1
C    3,THETM,CQSGL1,RTARRX 
C    4,RRADS,FILM(M,L,NY,NX),RRADL(N,L),RRAD2M(N,NZ,NY,NX) 
C    5,RTDNP(N,L,NZ,NY,NX)
C    4,DFGS(M,L,NY,NX),CH4G1,VOLWCH,CH4S1,VOLPMM,THETPM(M,L,NY,NX)
5566  FORMAT(A8,8I4,50E12.4)
C      
90    CONTINUE
C
C     ROOT SALT UPTAKE
C
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C     C*S1,C*P1=salt concentration in soil, root (mol m-3)
C     DIF*=aqueous diffusivity from soil to root (m3 t-1)
C     RMF*=soil convective solute flux (g t-1)
C     RDF*=aqueous soil-root diffusion (mol t-1)
C     UPWTRH=root water uptake (m3 t-1)
C     RUP*=root uptake (mol t-1) 
C     salt code:AL=Al,FE=Fe,Ca=Ca,MG=Mg,NA=Na,KA=K,SO=SO4,CL=Cl
C
      IF(ISALTG.NE.0)THEN
      CZALS1=AMAX1(0.0,ZALS1/VOLWMM)
      CZALP1=AMAX1(0.0,ZALP1/RTVLW(N,L,NZ,NY,NX))
      DIFAL=THETM*ALSGL1*RTARRX
      RMFZAL=UPWTRH*CZALS1
      RDFZAL=RMFZAL+DIFAL*(CZALS1-CZALP1) 
      RDXZAL=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,ZALS1)
     2-VOLWMM*AMAX1(0.0,ZALP1))/VOLWSP
      IF(RDFZAL.GT.0.0)THEN
      RUPZAX=AMIN1(AMAX1(0.0,RDXZAL),RDFZAL*PP(NZ,NY,NX))
     2/(1.0+CZALP1/CZALKI)
      ELSE
      RUPZAX=AMAX1(AMIN1(0.0,RDXZAL),RDFZAL*PP(NZ,NY,NX))
      ENDIF
      ZALS1=ZALS1-RUPZAX
      ZALP1=ZALP1+RUPZAX
      RUPZAL(N,L,NZ,NY,NX)=RUPZAL(N,L,NZ,NY,NX)+RUPZAX
      CZFES1=AMAX1(0.0,ZFES1/VOLWMM)
      CZFEP1=AMAX1(0.0,ZFEP1/RTVLW(N,L,NZ,NY,NX))
      DIFFE=THETM*FESGL1*RTARRX
      RMFZFE=UPWTRH*CZFES1
      RDFZFE=RMFZFE+DIFFE*(CZFES1-CZFEP1) 
      RDXZFE=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,ZFES1)
     2-VOLWMM*AMAX1(0.0,ZFEP1))/VOLWSP
      IF(RDFZFE.GT.0.0)THEN
      RUPZFX=AMIN1(AMAX1(0.0,RDXZFE),RDFZFE*PP(NZ,NY,NX))
     2/(1.0+CZFEP1/CZFEKI)
      ELSE
      RUPZFX=AMAX1(AMIN1(0.0,RDXZFE),RDFZFE*PP(NZ,NY,NX))
      ENDIF
      ZFES1=ZFES1-RUPZFX
      ZFEP1=ZFEP1+RUPZFX
      RUPZFE(N,L,NZ,NY,NX)=RUPZFE(N,L,NZ,NY,NX)+RUPZFX
      CZCAS1=AMAX1(0.0,ZCAS1/VOLWMM)
      CZCAP1=AMAX1(0.0,ZCAP1/RTVLW(N,L,NZ,NY,NX))
      DIFCA=THETM*CASGL1*RTARRX
      RMFZCA=UPWTRH*CZCAS1
      RDFZCA=RMFZCA+DIFCA*(CZCAS1-CZCAP1) 
      RDXZCA=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,ZCAS1)
     2-VOLWMM*AMAX1(0.0,ZCAP1))/VOLWSP
      IF(RDFZCA.GT.0.0)THEN
      RUPZCX=AMIN1(AMAX1(0.0,RDXZCA),RDFZCA*PP(NZ,NY,NX))
     2/(1.0+CZCAP1/CZCAKI)
      ELSE
      RUPZCX=AMAX1(AMIN1(0.0,RDXZCA),RDFZCA*PP(NZ,NY,NX))
      ENDIF
      ZCAS1=ZCAS1-RUPZCX
      ZCAP1=ZCAP1+RUPZCX
      RUPZCA(N,L,NZ,NY,NX)=RUPZCA(N,L,NZ,NY,NX)+RUPZCX
      CZMGS1=AMAX1(0.0,ZMGS1/VOLWMM)
      CZMGP1=AMAX1(0.0,ZMGP1/RTVLW(N,L,NZ,NY,NX))
      DIFMG=THETM*GMSGL1*RTARRX
      RMFZMG=UPWTRH*CZMGS1
      RDFZMG=RMFZMG+DIFMG*(CZMGS1-CZMGP1) 
      RDXZMG=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,ZMGS1)
     2-VOLWMM*AMAX1(0.0,ZMGP1))/VOLWSP
      IF(RDFZMG.GT.0.0)THEN
      RUPZMX=AMIN1(AMAX1(0.0,RDXZMG),RDFZMG*PP(NZ,NY,NX))
     2/(1.0+CZMGP1/CZMGKI)
      ELSE
      RUPZMX=AMAX1(AMIN1(0.0,RDXZMG),RDFZMG*PP(NZ,NY,NX))
      ENDIF
      ZMGS1=ZMGS1-RUPZMX
      ZMGP1=ZMGP1+RUPZMX
      RUPZMG(N,L,NZ,NY,NX)=RUPZMG(N,L,NZ,NY,NX)+RUPZMX
      CZNAS1=AMAX1(0.0,ZNAS1/VOLWMM)
      CZNAP1=AMAX1(0.0,ZNAP1/RTVLW(N,L,NZ,NY,NX))
      DIFNA=THETM*ANSGL1*RTARRX
      RMFZNA=UPWTRH*CZNAS1
      RDFZNA=RMFZNA+DIFNA*(CZNAS1-CZNAP1) 
      RDXZNA=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,ZNAS1)
     2-VOLWMM*AMAX1(0.0,ZNAP1))/VOLWSP
      IF(RDFZNA.GT.0.0)THEN
      RUPZNX=AMIN1(AMAX1(0.0,RDXZNA),RDFZNA*PP(NZ,NY,NX))
     2/(1.0+CZNAP1/CZNAKI)
      ELSE
      RUPZNX=AMAX1(AMIN1(0.0,RDXZNA),RDFZNA*PP(NZ,NY,NX))
      ENDIF
      ZNAS1=ZNAS1-RUPZNX
      ZNAP1=ZNAP1+RUPZNX
      RUPZNA(N,L,NZ,NY,NX)=RUPZNA(N,L,NZ,NY,NX)+RUPZNX
      CZKAS1=AMAX1(0.0,ZKAS1/VOLWMM)
      CZKAP1=AMAX1(0.0,ZKAP1/RTVLW(N,L,NZ,NY,NX))
      DIFKA=THETM*AKSGL1*RTARRX
      RMFZKA=UPWTRH*CZKAS1
      RDFZKA=RMFZKA+DIFKA*(CZKAS1-CZKAP1) 
      RDXZKA=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,ZKAS1)
     2-VOLWMM*AMAX1(0.0,ZKAP1))/VOLWSP
      IF(RDFZKA.GT.0.0)THEN
      RUPZKX=AMIN1(AMAX1(0.0,RDXZKA),RDFZKA*PP(NZ,NY,NX))
     2/(1.0+CZKAP1/CZKAKI)
      ELSE
      RUPZKX=AMAX1(AMIN1(0.0,RDXZKA),RDFZKA*PP(NZ,NY,NX))
      ENDIF
      ZKAS1=ZKAS1-RUPZKX
      ZKAP1=ZKAP1+RUPZKX
      RUPZKA(N,L,NZ,NY,NX)=RUPZKA(N,L,NZ,NY,NX)+RUPZKX
      CZSOS1=AMAX1(0.0,ZSOS1/VOLWMM)
      CZSOP1=AMAX1(0.0,ZSOP1/RTVLW(N,L,NZ,NY,NX))
      DIFSO=THETM*SOSGL1*RTARRX
      RMFZSO=UPWTRH*CZSOS1
      RDFZSO=RMFZSO+DIFSO*(CZSOS1-CZSOP1) 
      RDXZSO=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,ZSOS1)
     2-VOLWMM*AMAX1(0.0,ZSOP1))/VOLWSP
      IF(RDFZSO.GT.0.0)THEN
      RUPZSX=AMIN1(AMAX1(0.0,RDXZSO),RDFZSO*PP(NZ,NY,NX))
     2/(1.0+CZSOP1/CZSOKI)
      ELSE
      RUPZSX=AMAX1(AMIN1(0.0,RDXZSO),RDFZSO*PP(NZ,NY,NX))
      ENDIF
      ZSOS1=ZSOS1-RUPZSX
      ZSOP1=ZSOP1+RUPZSX
      RUPZSO(N,L,NZ,NY,NX)=RUPZSO(N,L,NZ,NY,NX)+RUPZSX
      CZCLS1=AMAX1(0.0,ZCLS1/VOLWMM)
      CZCLP1=AMAX1(0.0,ZCLP1/RTVLW(N,L,NZ,NY,NX))
      DIFCL=THETM*CLSGL1*RTARRX
      RMFZCL=UPWTRH*CZCLS1
      RDFZCL=RMFZCL+DIFCL*(CZCLS1-CZCLP1) 
      RDXZCL=(RTVLW(N,L,NZ,NY,NX)*AMAX1(0.0,ZCLS1)
     2-VOLWMM*AMAX1(0.0,ZCLP1))/VOLWSP
      IF(RDFZCL.GT.0.0)THEN
      RUPZCX=AMIN1(AMAX1(0.0,RDXZCL),RDFZCL*PP(NZ,NY,NX))
     2/(1.0+CZCLP1/CZCLKI)
      ELSE
      RUPZCX=AMAX1(AMIN1(0.0,RDXZCL),RDFZCL*PP(NZ,NY,NX))
      ENDIF
      ZCLS1=ZCLS1-RUPZCX
      ZCLP1=ZCLP1+RUPZCX
      RUPZCL(N,L,NZ,NY,NX)=RUPZCL(N,L,NZ,NY,NX)+RUPZCX
C     IF(NZ.EQ.2.AND.L.EQ.4)THEN
C     WRITE(*,5546)'RUPZAL',I,J,NFZ,NX,NY,NZ,L,N,M
C    2,RUPZAL(N,L,NZ,NY,NX),RUPZAX,RMFZAL,DIFAL*(CZALS1-CZALP1)
C    3,RDFZAL*PP(NZ,NY,NX)
C    3,RDXZAL,DIFAL,RTVLW(N,L,NZ,NY,NX),VOLWMM,VOLWSP
C    3,ZALS1,ZALP1,CZALS1,CZALP1,CZALKI,ZAL(L,NY,NX)
C     WRITE(*,5546)'RUPZSO',I,J,NFZ,NX,NY,NZ,L,N,M
C    2,RUPZSO(N,L,NZ,NY,NX),RUPZSX,RMFZSO,DIFSO*(CZSOS1-CZSOP1)
C    3,RDFZSO,RDFZSO*PP(NZ,NY,NX)
C    3,RDXZSO,DIFSO,RTVLW(N,L,NZ,NY,NX),VOLWMM,VOLWSP
C    3,ZSOS1,ZSOP1,CZSOS1,CZSOP1,CZSOKI,ZSO4(L,NY,NX),FPQ(N,L,NZ)
C    4,ZSOR(N,L,NZ,NY,NX)
5546  FORMAT(A8,9I4,30E12.4)
C     ENDIF
      ENDIF
      ENDIF
99    CONTINUE
C
C     O2 CONSTRAINTS TO ROOT RESPIRATION DEPENDS UPON RATIO
C     OF ROOT O2 UPTAKE 'RUPOXT' TO ROOT O2 DEMAND 'ROXYP'
C
C     RUPOXT=O2-limited uptake from soil+root by each root,mycorrhizal
C        population (g O t-1)
C     ROXYP=O2-unlimited uptake from soil+root by each root,mycorrhizal
C        population (g O t-1)
C     WFR=constraint by O2 consumption on all root processes
C        imposed by O2 uptake
C     OSTRD,OSTRN=O2 demand,uptake for entire soil profile (g O t-1)
C
      RUPOXT=RUPOXP(N,L,NZ,NY,NX)+RUPOXS(N,L,NZ,NY,NX)
      WFR(N,L,NZ,NY,NX)=AMIN1(1.0,AMAX1(0.0
     2,RUPOXT/ROXYP(N,L,NZ,NY,NX)))
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
C     IF(NZ.EQ.4.AND.(OSTRN.LT.0.0.OR.OSTRD.LT.0.0))THEN
C     WRITE(*,3368)'WFR',I,J,NFZ,NX,NY,NZ,L,N,WFR(N,L,NZ,NY,NX)
C    2,RUPOXP(N,L,NZ,NY,NX),RUPOXS(N,L,NZ,NY,NX),RUPOXT
C    3,ROXYP(N,L,NZ,NY,NX),OSTRN,OSTRD
3368  FORMAT(A8,8I4,12E12.4)
C     ENDIF
C
C     ROOT EXUDATION OF C, N AND P DEPENDS ON CONCN DIFFERENCES
C     BETWEEN ROOT NON-STRUCTURAL POOLS AND SOIL DISSOLVED POOLS
C
C     VOLWMM=soil micropore water volume (m3)
C     FOSAH=fraction of total SOC in each substrate K from ‘nitro.f’
C     RTVLW=root aqueous volume (m3)
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P in root,mycorrhizae 
C        (g C,N,P)
C     XFRC,XFRN,XFRP=nonstructural C,N,P exchange at root-soil 
C        DOC equilibrium (g C,N,P t-1)
C     OQC=soil DOC (g C)
C     RDFOMC,RDFOMN,RDFOMP=nonstructural C,N,P exchange (g C,N,P t-1)
C        :-ve=exudn,+ve=uptake
C     FEXU=rate constant for root C,N,P exudation (h-1)
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h t-1)
C     
      DO 195 K=0,4
      VOLWK=VOLWM(NPH,L,NY,NX)*FOSAH(K,L,NY,NX)
      IF(VOLWK.GT.ZEROS2(NY,NX)
     2.AND.RTVLW(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      VOLWT=VOLWK+RTVLW(N,L,NZ,NY,NX)
      CPOOLX=AMIN1(1.0E+03*RTVLW(N,L,NZ,NY,NX),CPOOLR(N,L,NZ,NY,NX))
      XFRC=(OQC(K,L,NY,NX)*RTVLW(N,L,NZ,NY,NX)
     2-CPOOLX*VOLWK)/VOLWT
      RDFOMC(N,K,L,NZ,NY,NX)=FEXU*XFRC*XNFH
C     IF(I.EQ.146.AND.L.EQ.1)THEN
C     WRITE(*,2224)'RDFOMC',I,J,NX,NY,L,NZ,K,N,RDFOMC(N,K,L,NZ,NY,NX)
C    2,RDFOMN(N,K,L,NZ,NY,NX),RDFOMP(N,K,L,NZ,NY,NX) 
C    3,OQC(K,L,NY,NX),OQN(K,L,NY,NX),OQP(K,L,NY,NX)
C    2,CPOOLR(N,L,NZ,NY,NX),ZPOOLR(N,L,NZ,NY,NX),PPOOLR(N,L,NZ,NY,NX)
C    3,VOLWM(NPH,L,NY,NX),RTVLW(N,L,NZ,NY,NX),RTAR1X(N,NZ,NY,NX)
C    4,RTAR2X(N,NZ,NY,NX),RTLGP(N,L,NZ,NY,NX)*PP(NZ,NY,NX)
C    4,WTRTD(N,L,NZ,NY,NX)
C    5,VOLWK,VOLWM(NPH,L,NY,NX),FOSAH(K,L,NY,NX)
C    5,OQC(K,L,NY,NX)/VOLWK
C    5,OQN(K,L,NY,NX)/OQC(K,L,NY,NX)
C    5,OQP(K,L,NY,NX)/OQC(K,L,NY,NX)
C    6,CPOOLR(N,L,NZ,NY,NX)/RTVLW(N,L,NZ,NY,NX)
C    6,ZPOOLX/CPOOLR(N,L,NZ,NY,NX)
C    6,PPOOLX/CPOOLR(N,L,NZ,NY,NX)
2224  FORMAT(A8,8I4,30E12.4)
C     ENDIF
      IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX)
     2.AND.CPOOLR(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CPOOLT=OQC(K,L,NY,NX)+CPOOLR(N,L,NZ,NY,NX)
      ZPOOLX=0.1*ZPOOLR(N,L,NZ,NY,NX)
      PPOOLX=0.1*PPOOLR(N,L,NZ,NY,NX)
      XFRN=(OQN(K,L,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
     2-ZPOOLX*OQC(K,L,NY,NX))/CPOOLT
      XFRP=(OQP(K,L,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
     2-PPOOLX*OQC(K,L,NY,NX))/CPOOLT
      RDFOMN(N,K,L,NZ,NY,NX)=FEXU*XFRN*XNFH 
      RDFOMP(N,K,L,NZ,NY,NX)=FEXU*XFRP*XNFH 
      ELSE
      RDFOMN(N,K,L,NZ,NY,NX)=0.0
      RDFOMP(N,K,L,NZ,NY,NX)=0.0
      ENDIF
      ELSE
      RDFOMC(N,K,L,NZ,NY,NX)=0.0
      RDFOMN(N,K,L,NZ,NY,NX)=0.0
      RDFOMP(N,K,L,NZ,NY,NX)=0.0
      ENDIF
195   CONTINUE
C
C     NUTRIENT UPTAKE
C
C     WFR=constraint by O2 consumption on all biological processes
C     FCUP=limitation to active uptake respiration from CPOOLR
C     FWSRT=protein concentration relative to 5% 
C     RTLGP=root,mycorrhizal length per plant (m)
C
      IF(WFR(N,L,NZ,NY,NX).GT.ZERO
     2.AND.FCUP.GT.ZERO.AND.FWSRT.GT.ZERO
     3.AND.RTLGP(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C
C     FZUP=limitation to active uptake respiration from CZPOLR
C
      IF(FZUP.GT.ZERO2)THEN
C
C     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NH4,NO3
C     FROM SOIL TO ROOT
C
C     ZNSGL=NH4 diffusivity from ‘hour1.f’ (m2 t-1)
C     TORT=soil tortuosity from ‘watsub.f’
C     PATH=path length of water and nutrient uptake (m)
C     RRADL=root radius (m)
C     DIFFL=NH4 diffusivity per plant (m3 t-1)
C     RTARR=root surface area/radius for uptake,diffusivity (m)
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      ZNSGX=ZNSGL(L,NY,NX)*TORT(NPH,L,NY,NX)*XNFH 
      PATHL=AMIN1(PATH(N,L),SQRT(2.0*ZNSGX))
      DIFFL=ZNSGX*RTARR(N,L)/LOG((PATHL+RRADL(N,L))/RRADL(N,L)) 
C
C     NH4 UPTAKE IN NON-BAND SOIL ZONE
C
C     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
C     CNH4S=NH4 concentration in non-band (g N m-3)
C     UPMXZH,UPKMZH,UPMNZH=NH4 maximum uptake,Km,minimum concentration
C        from PFT file (g N m-2 h-1,g N m-3,g N m-3)
C     UPWTRP=root water uptake per plant (m3 t-1)
C     RMFNH4=soil-root convective NH4 flux per plant in non-band 
C        (g N t-1)
C     DIFNH4=soil-root NH4 diffusivity per plant in non-band (m3 t-1)       
C     DIFFL=NH4 diffusivity per plant (m3 t-1)
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
C        (g N t-1)
C     RTARP=root surface area per plant from ‘grosub.f’ (m2)
C     FWSRT=protein concentration relative to 5% 
C     TFN4=temperature function for root growth 
C     FCUP,FZUP=limitation to active uptake respiration from
C        CCPOLR,CZPOLR
C     WFR=constraint by O2 consumption on all biological processes
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      UPMXP=UPMXZH(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) 
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLNH4(L,NY,NX)*AMIN1(FCUP,FZUP)*XNFH
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX) 
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFNH4=soil-root convective NH4 flux per plant in non-band 
C        (g N t-1)
C     DIFNH4=soil-root NH4 diffusivity per plant in non-band (m3 t-1)      
C     CNH4S=NH4 concentration in non-band (g N m-3)
C     UPMXZH,UPKMZH,UPMNZH=NH4 maximum uptake,Km,minimum concentration
C        from PFT file (g N m-2 h-1,g N m-3,g N m-3)
C     RTKNH4,RTKNHP=NH4 uptake per plant in non-band limited,unlimited
C        by O2 (g N t-1)
C     ZNH4M,ZNH4X=minimum,maximum NH4 available for uptake in non-band
C        (g N)
C     FNH4X=fraction of total NH4 uptake in non-band by
C        root,mycorrhizal population
C     RUNNHP,RUPNH4=NH4 uptake in non-band unlimited,limited by NH4 
C        (g N t-1)
C     RUONH4=NH4 uptake in non-band unlimited by O2 (g N t-1)
C     RUCNH4=NH4 uptake in non-band unlimited by nonstructural C 
C        (g N t-1)  
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
      ZNH4X=AMAX1(0.0,FNH4X*(ZNH4S(L,NY,NX)-ZNH4M)*XNFH)
      RUNNHP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNH4*PP(NZ,NY,NX))
      RUPNH4(N,L,NZ,NY,NX)=AMIN1(ZNH4X,RUNNHP(N,L,NZ,NY,NX))
      RUONH4(N,L,NZ,NY,NX)=AMIN1(ZNH4X,AMAX1(0.0
     2,RTKNHP*PP(NZ,NY,NX)))
      RUCNH4(N,L,NZ,NY,NX)=RUPNH4(N,L,NZ,NY,NX)/FCUP
C     IF(NX.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,1110)'UPNH4',I,J,NZ,L,N,RUNNHP(N,L,NZ,NY,NX)
C    2,RUPNH4(N,L,NZ,NY,NX),RTKNH4,DIFNH4,RMFNH4,X,Y,B,C,UPMX,UPMXP
C    2,DIFFL,ZNSGX,RTARR(N,L),PATHL,RRADL(N,L)
C    3,ZNSGL(L,NY,NX),TORT(NPH,L,NY,NX),XNFH,PATH(N,L) 
C    2,WFR(N,L,NZ,NY,NX),CNH4S(L,NY,NX),DIFNH,RTDNP(N,L,NZ,NY,NX) 
C    2,WTRTD(N,L,NZ,NY,NX),CNH4S(L,NY,NX) 
C    3,CCPOLR(N,L,NZ,NY,NX),CZPOLR(N,L,NZ,NY,NX),CPPOLR(N,L,NZ,NY,NX)
C    4,THETW(L,NY,NX),TKS(L,NY,NX),RSCS(L,NY,NX),UPMXP,FWSRT
C    5,FZUP,FCUP,COXYS(L,NY,NX),COXYG(L,NY,NX) 
C    6,CCPOLP(NZ,NY,NX)
C    7,CZPOLP(NZ,NY,NX),CPPOLP(NZ,NY,NX),FDBK(1,NZ,NY,NX),PSIST1(L)
C    2,PSIRT(N,L,NZ,NY,NX),ZPOOLR(N,L,NZ,NY,NX),WTRTL(N,L,NZ,NY,NX)
C    3,RTARP(N,L,NZ,NY,NX),RRADL(N,L),PATH(N,L)
C    4,DIFFL,ZNSGX,RTARR(N,L),PATHL,RRADL(N,L),VLNH4(L,NY,NX)
1110  FORMAT(A8,5I4,100E16.8)
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
C     CNH4B=NH4 concentration in band (g N m-3)
C     UPMXZH,UPKMZH,UPMNZH=NH4 maximum uptake,Km,minimum concentration
C        from PFT file (g N m-2 h-1,g N m-3,g N m-3)
C     UPWTRP=root water uptake per plant (m3 t-1)
C     RMFNHB=soil-root convective NH4 flux per plant in band
C        (g N t-1)
C     DIFNHB=soil-root NH4 diffusivity per plant in band (m3 t-1)      
C     DIFFL=NH4 diffusivity per plant (m3 t-1)
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
C        (g N t-1)
C     RTARP=root surface area per plant from ‘grosub.f’ (m2)
C     FWSRT=protein concentration relative to 5% 
C     TFN4=temperature function for root growth 
C     FCUP,FZUP=limitation to active uptake respiration from
C        CCPOLR,CZPOLR
C     WFR=constraint by O2 consumption on all biological processes
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      UPMXP=UPMXZH(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX) 
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLNHB(L,NY,NX)*AMIN1(FCUP,FZUP)*XNFH
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX) 
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFNHB=soil-root convective NH4 flux per plant in band
C        (g N t-1)
C     DIFNHB=soil-root NH4 diffusivity per plant in band (m3 t-1)      
C     CNH4B=NH4 concentration in band (g N m-3)
C     UPMXZH,UPKMZH,UPMNZH=NH4 maximum uptake,Km,minimum concentration
C        from PFT file (g N m-2 h-1,g N m-3,g N m-3)
C     RTKNHB,RTKNBP=NH4 uptake per plant in band limited,unlimited
C        by O2 (g N t-1)
C     ZNHBM,ZNHBX=minimum,maximum NH4 available for uptake in band 
C        (g N) 
C     FNHBX=fraction of total NH4 uptake in band by root,mycorrhizal
C        population
C     RUNNBP,RUPNHB=NH4 uptake in band unlimited,limited by NH4
C        (g N t-1)
C     RUONHB=NH4 uptake in band unlimited by O2 (g N t-1)
C     RUCNHB=NH4 uptake in band unlimited by nonstructural C  
C        (g N t-1)
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
      ZNHBX=AMAX1(0.0,FNHBX*(ZNH4B(L,NY,NX)-ZNHBM)*XNFH)
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
C     ZOSGL=NO3 diffusivity from ‘hour1.f’ (m2 t-1)
C     TORT=soil tortuosity from ‘watsub.f’
C     RRADL=root radius (m)
C     PATH=path length of water and nutrient uptake (m)
C     DIFFL=NO3 diffusivity per plant (m3 t-1)
C     RTARR=root surface area/radius for uptake,diffusivity (m)
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      ZOSGX=ZOSGL(L,NY,NX)*TORT(NPH,L,NY,NX)*XNFH
      PATHL=AMIN1(PATH(N,L),SQRT(2.0*ZOSGX))
      DIFFL=ZOSGX*RTARR(N,L)/LOG((PATHL+RRADL(N,L))/RRADL(N,L))
C
C     NO3 UPTAKE IN NON-BAND SOIL ZONE
C
C     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
C     CNO3S=NO3 concentration in non-band (g N m-3)
C     UPMXZO,UPKMZO,UPMNZO=NO3 maximum uptake,Km,minimum concentration 
C        from PFT file (g N m-2 h-1,g N m-3,g N m-3)
C     UPWTRP=root water uptake per plant (m3 t-1)
C     RMFNO3=soil-root convective NO3 flux per plant in non-band
C        (g N t-1)
C     DIFNO3=soil-root NO3 diffusivity per plant in non-band (m3 t-1)       
C     DIFFL=NO3 diffusivity per plant (m3 t-1)
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
C        (g N t-1)
C     RTARP=root surface area per plant from ‘grosub.f’ (m2)
C     FWSRT=protein concentration relative to 5% 
C     TFN4=temperature function for root growth 
C     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
C     WFR=constraint by O2 consumption on all biological processes
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      UPMXP=UPMXZO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLNO3(L,NY,NX)*AMIN1(FCUP,FZUP)*XNFH 
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX) 
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFNO3=soil-root convective NO3 flux per plant in non-band
C        (g N t-1)
C     DIFNO3=soil-root N03 diffusivity per plant in non-band (m3 t-1)      
C     CNO3S=NO3 concentration in non-band (g N m-3)
C     UPMXZO,UPKMZO,UPMNZO=NO3 maximum uptake,Km,minimum concentration
C        from PFT file (g N m-2 h-1,g N m-3,g N m-3)
C     RTKNO3,RTKNOP=NO3 uptake per plant in non-band limited,unlimited
C        by O2 (g N t-1)
C     ZNO3M,ZNO3X=minimum,maximum NO3 available for uptake in non-band 
C        (g N)
C     FNO3X=fraction of total NH4 uptake in non-band by root,
C        mycorrhizal population
C     RUNNOP,RUPNO3=NO3 uptake in non-band unlimited,limited by NO3
C        (g N t-1)
C     RUONO3=NO3 uptake in non-band unlimited by O2 (g N t-1)
C     RUCNO3=NO3 uptake in non-band unlimited by nonstructural C  
C        (g N t-1)
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
      ZNO3X=AMAX1(0.0,FNO3X*(ZNO3S(L,NY,NX)-ZNO3M)*XNFH)
      RUNNOP(N,L,NZ,NY,NX)=AMAX1(0.0,RTKNO3*PP(NZ,NY,NX))
      RUPNO3(N,L,NZ,NY,NX)=AMIN1(ZNO3X,RUNNOP(N,L,NZ,NY,NX))
      RUONO3(N,L,NZ,NY,NX)=AMIN1(ZNO3X,AMAX1(0.0
     2,RTKNOP*PP(NZ,NY,NX)))
      RUCNO3(N,L,NZ,NY,NX)=RUPNO3(N,L,NZ,NY,NX)/FCUP
C     IF(L.EQ.12)THEN
C     WRITE(*,1111)'UPNO3',I,J,NFZ,NZ,L,N,RUPNO3(N,L,NZ,NY,NX),FNO3X
C    2,ZNO3S(L,NY,NX),ZNO3M,RTDNP(N,L,NZ,NY,NX),RTKNO3,RMFNO3,X,Y,B,C 
C    2,UPMX,CNO3S(L,NY,NX),DIFNO3,DIFFL,RTARR(N,L),PATHL,RRADL(N,L) 
C    3,CCPOLR(N,L,NZ,NY,NX),CZPOLR(N,L,NZ,NY,NX),CPPOLR(N,L,NZ,NY,NX)
C    4,THETW(L,NY,NX),TKS(L,NY,NX),RSCS(L,NY,NX),UPMXP,FWSRT
C    5,FZUP,FCUP,COXYS(L,NY,NX),COXYG(L,NY,NX),WFR(N,L,NZ,NY,NX)
C    6,CCPOLP(NZ,NY,NX),CZPOLP(NZ,NY,NX),CPPOLP(NZ,NY,NX)
C    7,FDBK(1,NZ,NY,NX),PSIST1(L),PSIRT(N,L,NZ,NY,NX),ZOSGX
C    2,ZPOOLR(N,L,NZ,NY,NX),WTRTL(N,L,NZ,NY,NX),VLNO3(L,NY,NX)
C    3,RUONO3(N,L,NZ,NY,NX),RUNNOP(N,L,NZ,NY,NX),RNO3Y(L,NY,NX)
1111  FORMAT(A8,6I4,40E12.4)
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
C     CNO3B=NO3 concentration in band (g N m-3)
C     UPMXZO,UPKMZO,UPMNZO=NO3 maximum uptake,Km,minimum concentration
C        from PFT file (g N m-2 h-1,g N m-3,g N m-3)
C     UPWTRP=root water uptake per plant (m3 t-1)
C     RMFNOB=soil-root convective NO3 flux per plant in band
C        (g N t-1)
C     DIFNOB=soil-root NO3 diffusivity per plant in band (m3 t-1)       
C     DIFFL=NO3 diffusivity per plant (m3 t-1)
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
C        (g N t-1)
C     RTARP=root surface area per plant from ‘grosub.f’ (m2)
C     FWSRT=protein concentration relative to 5% 
C     TFN4=temperature function for root growth 
C     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
C     WFR=constraint by O2 consumption on all biological processes
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      UPMXP=UPMXZO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLNOB(L,NY,NX)*AMIN1(FCUP,FZUP)*XNFH
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX) 
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFNOB=soil-root convective NO3 flux per plant in band
C        (g N t-1)
C     DIFNOB=soil-root NO3 diffusivity per plant in band (m3 t-1)      
C     CNO3B=NH4 concentration in band (g N m-3)
C     UPMXZO,UPKMZO,UPMNZO=NO3 maximum uptake,Km,minimum concentration
C        from PFT file (g N m-2 h-1,g N m-3,g N m-3)
C     RTKNOB,RTKNPB=NO3 uptake per plant in band limited,unlimited
C        by O2 (g N t-1)
C     ZNOBM,ZNOBX=minimum,maximum NO3 available for uptake in band 
C        (g N) 
C     FNOBX=fraction of total NO3 uptake in band by root,mycorrhizal
C        population
C     RUNNXP,RUPNOB=NO3 uptake in band unlimited,limited by NH4
C        (g N t-1)
C     RUONOB=NO3 uptake in band unlimited by O2 (g N t-1)
C     RUCNOB=NO3 uptake in band unlimited by nonstructural C  
C        (g N t-1)
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
      ZNOBX=AMAX1(0.0,FNOBX*(ZNO3B(L,NY,NX)-ZNOBM)*XNFH)
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
C     FPUP=limitation to active uptake respiration from CPPOLR
C
      IF(FPUP.GT.ZERO2)THEN
C
C     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF H2PO4,HPO4
C     FROM SOIL TO ROOT
C
C     POSGL=PO4 diffusivity from ‘hour1.f’ (m2 t-1)
C     TORT=soil tortuosity from ‘watsub.f’
C     PATH=path length of water and nutrient uptake (m)
C     RRADL=root radius (m)
C     DIFFL=PO4 diffusivity per plant (m3 t-1)
C     RTARR=root surface area/radius for uptake,diffusivity (m)
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      POSGX=POSGL(L,NY,NX)*TORT(NPH,L,NY,NX)*XNFH
      PATHL=AMIN1(PATH(N,L),SQRT(2.0*POSGX))
      DIFFL=POSGX*RTARR(N,L)/LOG((PATHL+RRADL(N,L))/RRADL(N,L))
C
C     H2PO4 UPTAKE IN NON-BAND SOIL ZONE
C
C     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
C     CH2P4=H2PO4 concentration in non-band (g P m-3)
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 maximum uptake,Km,minimum
C        concentration from PFT file (g P m-2 h-1,g P m-3,g P m-3)
C     UPWTRP=root water uptake per plant (m3 t-1)
C     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
C        (g P t-1)
C     DIFH2P=soil-root H2PO4 diffusivity per plant in non-band (m3 t-1)       
C     DIFFL=PO4 diffusivity per plant (m3 t-1)
C
      IF(VLPO4(L,NY,NX).GT.ZERO.AND.CH2P4(L,NY,NX)
     2.GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH2P=UPWTRP*CH2P4(L,NY,NX)*VLPO4(L,NY,NX)
      DIFH2P=DIFFL*VLPO4(L,NY,NX)
C
C     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 
C     'READQ'AND FROM ROOT SURFACE AREA, C AND P CONSTRAINTS
C     CALCULATED ABOVE
C
C     UPMXP,UPMX=maximum H2PO4 uptake in non-band unlimited,limited 
C        by O2 (g P t-1)
C     RTARP=root surface area per plant from ‘grosub.f’(m2)
C     FWSRT=protein concentration relative to 5% 
C     TFN4=temperature function for root growth 
C     FCUP,FPUP=limitation to active uptake respiration from
C        CCPOLR,CPPOLR
C     WFR=constraint by O2 consumption on all biological processes
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      UPMXP=UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLPO4(L,NY,NX)*AMIN1(FCUP,FPUP)*XNFH
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX) 
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF H2PO4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY 
C     COMPETITION WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
C        (g P t-1)
C     DIFH2P=soil-root H2PO4 diffusivity per plant in non-band (m3 t-1)      
C     CH2P4=H2PO4 concentration in non-band (g P m-3)
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 maximum uptake,Km,minimum
C        concentration from PFT file (g P m-2 h-1,g P m-3,g P m-3)
C     RTKH2P,RTKHPP=H2PO4 uptake per plant in non-band
C        limited,unlimited by O2 (g P t-1)
C     H2POM,H2POX=minimum,maximum H2PO4 available for uptake 
C        in non-band (g P)
C     FPO4X=fraction of total H2PO4 uptake in non-band by
C        root,mycorrhizal popululation
C     RUPP2P,RUPH2P=H2PO4 uptake in non-band unlimited,limited 
C        by H2PO4 (g P t-1)
C     RUOH2P=H2PO4 uptake in non-band unlimited by O2 (g P t-1)
C     RUCH2P=H2PO4 uptake in non-band unlimited by nonstructural C 
C        (g P t-1)  
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
      H2POX=AMAX1(0.0,FPO4X*(H2PO4(L,NY,NX)-H2POM)*XNFH)
      RUPP2P(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH2P*PP(NZ,NY,NX))
      RUPH2P(N,L,NZ,NY,NX)=AMIN1(H2POX,RUPP2P(N,L,NZ,NY,NX))
      RUOH2P(N,L,NZ,NY,NX)=AMIN1(H2POX,AMAX1(0.0
     2,RTKHPP*PP(NZ,NY,NX)))
      RUCH2P(N,L,NZ,NY,NX)=RUPH2P(N,L,NZ,NY,NX)/FCUP
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.1)THEN
C     WRITE(*,2223)'UPPO4',I,J,NZ,L,N,RUPH2P(N,L,NZ,NY,NX) 
C    2,FPO4X,H2PO4(L,NY,NX),RUPP2P(N,L,NZ,NY,NX)
C    3,UPMX,DIFPO,UPKMPO(N,NZ,NY,NX)
C    3,UPMNPO(N,NZ,NY,NX),RMFH2P,CH2P4(L,NY,NX)
C    4,UPMXP,WFR(N,L,NZ,NY,NX)
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
C     CH2P4B=H2PO4 concentration in band (g P m-3)
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 maximum uptake,Km,minimum
C        concentration from PFT file (g P m-2 h-1,g P m-3,g P m-3)
C     UPWTRP=root water uptake per plant (m3 t-1)
C     RMFH2B=soil-root convective H2PO4 flux per plant in band 
C        (g P t-1)
C     DIFH2B=soil-root H2PO4 diffusivity per plant in band (m3 t-1)      
C     DIFFL=PO4 diffusivity per plant (m3 t-1)
C
      IF(VLPOB(L,NY,NX).GT.ZERO.AND.CH2P4B(L,NY,NX)
     2.GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH2B=UPWTRP*CH2P4B(L,NY,NX)*VLPOB(L,NY,NX)
      DIFH2B=DIFFL*VLPOB(L,NY,NX)
C
C     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 
C     'READQ' AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS
C     CALCULATED ABOVE
C
C     UPMXP,UPMX=maximum H2PO4 uptake in band unlimited,limited by O2
C        (g P t-1)
C     RTARP=root surface area per plant from ‘grosub.f’
C     FWSRT=protein concentration relative to 5% 
C     TFN4=temperature function for root growth 
C     FCUP,FPUP=limitation to active uptake respiration from
C        CCPOLR,CPPOLR
C     WFR=constraint by O2 consumption on all biological processes
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      UPMXP=UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLPOB(L,NY,NX)*AMIN1(FCUP,FPUP)*XNFH
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX) 
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF PO4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY 
C     COMPETITION WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFH2B=soil-root convective H2PO4 flux per plant in band
C        (g P t-1)
C     DIFH2B=soil-root H2PO4 diffusivity per plant in band (m3 t-1)      
C     CH2P4B=H2PO4 concentration in band (g P m-3)
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 maximum uptake,Km,minimum
C        concentration from PFT file (g P m-2 h-1,g P m-3,g P m-3)
C     RTKH2B,RTKHPB=H2PO4 uptake per plant in band limited,unlimited
C        by O2 (g P t-1)
C     H2PXM,H2PXB=minimum,maximum H2PO4 available for uptake in band 
C        (g P) 
C     FPOBX=fraction of total H2PO4 uptake in band by root,mycorrhizal
C        population
C     RUPP2B,RUPH2B=H2PO4 uptake in band unlimited,limited by H2PO4
C        (g P t-1)
C     RUOH2B=H2PO4 uptake in band unlimited by O2 (g P t-1)
C     RUCH2B=H2PO4 uptake in band unlimited by nonstructural C  
C        (g P t-1)
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
      H2PXB=AMAX1(0.0,FPOBX*(H2POB(L,NY,NX)-H2PXM)*XNFH)
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
C     CH1P4=HPO4 concentration in non-band (g P m-3)
C     UPMXPO,UPKMPO,UPMNPO=HPO4 maximum uptake,Km,minimum
C        concentration from PFT file (g P m-2 h-1,g P m-3,g P m-3)
C        assumed same as H2PO4
C     UPWTRP=root water uptake per plant (m3 t-1)
C     RMFH1P=soil-root convective HPO4 flux per plant in non-band
C        (g P t-1)
C     DIFH1P=soil-root HPO4 diffusivity per plant in non-band (m3 t-1)      
C     DIFFL=PO4 diffusivity per plant (m3 t-1)
C
      IF(VLPO4(L,NY,NX).GT.ZERO.AND.CH1P4(L,NY,NX)
     2.GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH1P=UPWTRP*CH1P4(L,NY,NX)*VLPO4(L,NY,NX)
      DIFH1P=DIFFL*VLPO4(L,NY,NX)
C
C     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 
C     'READQ' AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS
C     CALCULATED ABOVE
C
C     UPMXP,UPMX= maximum HPO4 uptake in non-band unlimited,limited 
C        by O2 (g P t-1)
C     RTARP=root surface area per plant from ‘grosub.f’ (m2)
C     FWSRT=protein concentration relative to 5% 
C     TFN4=temperature function for root growth 
C     FCUP,FPUP=limitation to active uptake respiration from
C        CCPOLR,CPPOLR
C     WFR=constraint by O2 consumption on all biological processes
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      UPMXP=0.1*UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLPO4(L,NY,NX)*AMIN1(FCUP,FPUP)*XNFH
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX) 
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFH1P=soil-root convective HPO4 flux per plant in non-band
C        (g P t-1)
C     DIFH1P=soil-root HPO4 diffusivity per plant in non-band (m3 t-1)      
C     CH1P4=HPO4 concentration in non-band (g P m-3)
C     UPMXPO,UPKMPO,UPMNPO=HPO4 maximum uptake,Km,minimum
C        concentration from PFT file (g P m-2 h-1,g P m-3,g P m-3)
C        assumed same as H2PO4
C     RTKH1P,RTKHP1=HPO4 uptake per plant in non-band
C        limited,unlimited by O2 (g P t-1)
C     H1POM,H1POX=minimum,maximum HPO4 available for uptake 
C        in non-band (g P)
C     FP14X=fraction of total HPO4 uptake in non-band by
C        root,mycorrhizal population
C     RUPP1P,RUPH1P=HPO4 uptake in non-band unlimited,limited by HPO4
C        (g P t-1)
C     RUOH1P=HPO4 uptake in non-band unlimited by O2 (g P t-1)
C     RUCH1P=HPO4 uptake in non-band unlimited by nonstructural C  
C        (g P t-1)
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
      H1POX=AMAX1(0.0,FP14X*(H1PO4(L,NY,NX)-H1POM)*XNFH)
      RUPP1P(N,L,NZ,NY,NX)=AMAX1(0.0,RTKH1P*PP(NZ,NY,NX))
      RUPH1P(N,L,NZ,NY,NX)=AMIN1(H1POX,RUPP1P(N,L,NZ,NY,NX))
      RUOH1P(N,L,NZ,NY,NX)=AMIN1(H1POX,AMAX1(0.0
     2,RTKHP1*PP(NZ,NY,NX)))
      RUCH1P(N,L,NZ,NY,NX)=RUPH1P(N,L,NZ,NY,NX)/FCUP
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.3)THEN
C     WRITE(*,2226)'UPPO4',I,J,NZ,L,N,RUPH2P(N,L,NZ,NY,NX),FPO4X
C    2,H2PO4(L,NY,NX),RUPP2P(N,L,NZ,NY,NX)
C    3,UPMX,DIFPO,UPKMPO(N,NZ,NY,NX)
C    3,UPMNPO(N,NZ,NY,NX),RMFH2P,CH2P4(L,NY,NX)
C    4,UPMXP,WFR(N,L,NZ,NY,NX)
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
C     CH1P4B=HPO4 concentration in band (g P m-3)
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 maximum uptake,Km,minimum
C        concentration from PFT file (g P m-2 h-1,g P m-3,g P m-3)
C        assumed same as H2PO4
C     UPWTRP=root water uptake per plant (m3 t-1)
C     RMFH1B=soil-root convective HPO4 flux per plant in band
C        (g P t-1)
C     DIFH1B=soil-root HPO4 diffusivity per plant in band (m3 t-1)      
C     DIFFL=PO4 diffusivity per plant (m3 t-1)
C
      IF(VLPOB(L,NY,NX).GT.ZERO.AND.CH1P4B(L,NY,NX)
     2.GT.UPMNPO(N,NZ,NY,NX))THEN
      RMFH2B=UPWTRP*CH1P4B(L,NY,NX)*VLPOB(L,NY,NX)
      DIFH1B=DIFFL*VLPOB(L,NY,NX)
C
C     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 
C     'READQ' AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS
C     CALCULATED ABOVE
C
C     UPMXP,UPMX=maximum HPO4 uptake in band unlimited,limited by O2
C        (g P t-1)
C     RTARP=root surface area per plant from ‘grosub.f’ (m2)
C     FWSRT=protein concentration relative to 5% 
C     TFN4=temperature function for root growth 
C     FCUP,FPUP=limitation to active uptake respiration from
C        CCPOLR,CPPOLR
C     WFR=constraint by O2 consumption on all biological processes
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h -1)
C
      UPMXP=0.1*UPMXPO(N,NZ,NY,NX)*RTARP(N,L,NZ,NY,NX)
     2*FWSRT*TFN4(L,NZ,NY,NX)*VLPOB(L,NY,NX)*AMIN1(FCUP,FPUP)*XNFH
      UPMX=UPMXP*WFR(N,L,NZ,NY,NX) 
C
C     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
C     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
C     WITH OTHER ROOT AND MICROBIAL POPULATIONS
C
C     RMFH1B=soil-root convective HPO4 flux per plant in band
C        (g P t-1)
C     DIFH1B=soil-root HPO4 diffusivity per plant in band (m3 t-1)     
C     CH1P4B=HPO4 concentration in band (g P m-3)
C     UPMXPO,UPKMPO,UPMNPO=HPO4 maximum uptake,Km,minimum
C        concentration from PFT file (g P m-2 h-1,g P m-3,g P m-3)
C        assumed same as H2PO4
C     RTKH1B,RTKHB1=HPO4 uptake per plant in band limited,unlimited 
C        by O2 (g P t-1)
C     H1PXM,H1PXB=minimum,maximum HPO4 available for uptake in band 
C        (g P)
C     FP1BX=fraction of total HPO4 uptake in band by root,mycorrhizal
C        population
C     RUPP1B,RUPH1B=HPO4 uptake in band unlimited,limited by H2PO4
C        (g P t-1)
C     RUOH1B=HPO4 uptake in band unlimited by O2 (g P t-1)
C     RUCH1B=HPO4 uptake in band unlimited by nonstructural C 
C        (g P t-1)
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
      H1PXB=AMAX1(0.0,FP1BX*(H1POB(L,NY,NX)-H1PXM)*XNFH)
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
C     TOTAL C,N,P EXCHANGE BETWEEN ROOTS AND SOIL FOR ALL N,K
C
C     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange
C        :-ve=exudn,+ve=uptake (g C,N,P t-1)
C     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructural C,N,P exchange
C        (g C,N,P t-1)
C     XOQCS,XOQNZ,XOQPS=accumulated change in DOC,DON,DOP 
C        from ‘nitro.f’(g C,N,P t-1)
C     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band 
C        of NH4,NO3 (g N t-1)
C     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band 
C        of H2PO4,HPO4 (g P t-1)
C     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
C        (g C,N,P t-1) 
C
      DO 295 K=0,4
      UPOMC(NZ,NY,NX)=UPOMC(NZ,NY,NX)+RDFOMC(N,K,L,NZ,NY,NX)
      UPOMN(NZ,NY,NX)=UPOMN(NZ,NY,NX)+RDFOMN(N,K,L,NZ,NY,NX)
      UPOMP(NZ,NY,NX)=UPOMP(NZ,NY,NX)+RDFOMP(N,K,L,NZ,NY,NX)
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-RDFOMC(N,K,L,NZ,NY,NX)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-RDFOMN(N,K,L,NZ,NY,NX)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-RDFOMP(N,K,L,NZ,NY,NX)
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NZ.EQ.1)THEN
C     WRITE(*,8766)'UPOMC',I,J,NX,NY,L,NZ,K,N
C    2,UPOMC(NZ,NY,NX),RDFOMC(N,K,L,NZ,NY,NX)
8766  FORMAT(A8,8I4,12E12.4)
C     ENDIF
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
C
C     ROOT O2 LIMITATION USED TO REPRESENT O2 STRESS IN MODEL OUTPUT
C
C     OSTR=O2 root O2-limited uptake:root O2-unlimited uptake 
C
      IF(OSTRD.GT.ZEROP(NZ,NY,NX))THEN
      OSTR(NZ,NY,NX)=OSTRN/OSTRD
      ELSE
      OSTR(NZ,NY,NX)=0.0
      ENDIF
      ELSE
C
C     DEFAULT VALUES IF PLANT SPECIES DOES NOT EXIST
C
      VOLWC(NZ,NY,NX)=VOLWC(NZ,NY,NX)+FLWC(NZ,NY,NX)*XNFH
      RFLXCC=0.0
      EFLXCC=0.0
      SFLXCC=0.0
      HFLXCC=0.0
      VFLXCC=0.0
      THRMCC=0.0
      EVAPCC=0.0
      EPCC=0.0
      TEVCG=0.0
      TSHCG=0.0
      TKQC(NZ,NY,NX)=TKQ(NY,NX)
      VPQC(NZ,NY,NX)=VPQ(NY,NX)
      IF(ZC(NZ,NY,NX).GE.DPTHS(NY,NX)-ZERO)THEN
      TKC(NZ,NY,NX)=TKQ(NY,NX)
      ELSE
      TKC(NZ,NY,NX)=TKW(1,NY,NX)
      ENDIF
      TCC(NZ,NY,NX)=TKC(NZ,NY,NX)-273.15
      ENDIF
C
C     INITIALIZE STANDING DEAD CANOPY TEMPERATURE AND VAPOR 
C     CONCENTRATION
C
C     ARSTG=standing dead surface area (m2)
C     TKQD,VPQD=standing dead air temperature, vapor concentration
C        (K,m3 m-3)
C
      IF(ARSTG(NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
C
      TKQY=TKQD(NZ,NY,NX)
      VPQY=VPQD(NZ,NY,NX)
C
C     CALCULATE CANOPY WATER STATUS FROM CONVERGENCE SOLUTION FOR
C     TRANSPIRATION - ROOT WATER UPTAKE = CHANGE IN CANOPY WATER
C     CONTENT
C
C     FLWD=standing dead water retention of precipitation from
C        ‘hour1.f’(m3 h-1) 
C     FLWCDM=standing dead water retention of precipitation (m3 t-1) 
C     VHCPDP,VHCPDX=dry,wet standing dead heat capacity (MJ K-1)
C     ARSTG,WTSTG=standing dead surface area, mass (m2,g C)
C     VSTK=stalk volume:mass from ‘startq.f’ (m3 g-1)
C     FARS=stalk sapwood thickness from ‘startq.f’(m)
C     VOLWQ=water on standing dead surfaces (m3)
C     VHCPXZ,VHCPYZ=minimum heat capacity for solving standing dead
C        water potential,aerodymamic energy exchange(MJ K-1)
C     AREA=grid cell area (m2) 
C     TKDY=standing dead surface temperature (K) 
C     DTKX=step change in standing dead surface temperature used in
C        convergence (K)
C     XNPHX,XNPH=time steps for water fluxes from ‘wthr.f’ (h t-1) 
C
      FLWCDM=FLWD(NZ,NY,NX)*XNPHX
      VHCPDP=2.496*AMIN1(FARS*ARSTG(NZ,NY,NX),WTSTG(NZ,NY,NX)*VSTK) 
      VHCPDX=VHCPDP+4.19*AMAX1(0.0,VOLWQ(NZ,NY,NX))
      VHCPXZ=VHCPXM*AREA(3,NU(NY,NX),NY,NX)
      VHCPYZ=VHCPYM*AREA(3,NU(NY,NX),NY,NX)
      TKDY=TKD(NZ,NY,NX)
      TKDX=TKDY 
      IF(VHCPDP.GT.VHCPXZ)THEN
      DTKX=DTKD1
      ELSE
      DTKX=DTKD2
      ENDIF
C
C     INITIALIZE STANDING DEAD ENERGY BALANCE
C
C     FRADQ=fraction of incoming radiation received by each PFT
C        standing dead from ‘hour1.f’ 
C     RADCD=SW radiation absorbed by standing dead (MJ t-1)
C     EMMC=standing dead emissivity
C     EMMDX=EMMC x SB constant (MJ h-1 m-2 K-4) x FRADQ x AREA (m2) 
C        x time step (MJ t-1 K-4)
C     AREA=grid cell area (m2)
C     THRMDXM,THRMDYM=LW radiation absorbed by standing dead
C        from sky,to adjacent grid cell (MJ t-1)
C     THS=sky LW radiation from ‘wthr.f’ (MJ h-1)
C     DTHRMT=total lateral LW radiation (MJ t-1)
C     FLAIQ=fraction of total canopy+standing dead radiation
C        received by each PFT standing dead from ‘hour1.f’
C     PAREX,PARSX=terms used to calculate canopy boundary 
C        layer conductance for latent,sensible heat from ‘hour1.f’ 
C        (m2 t-1,MJ m-1 K-1 t-1)
C     PAREY,PARSY=standing dead boundary layer conductance for
C        latent,sensible heat (m2 t-1,MJ m-1 K-1 t-1)
C     XNPHX=time step for heat fluxes from ‘wthr.f’ (h t-1)
C
      IF(VHCPDX.GT.VHCPYZ
     2.AND.FRADQ(NZ,NY,NX).GT.ZERO2)THEN 
      EVAPD=0.0
      RADCD=RADD(NZ,NY,NX)*XNPHX 
      EMMDX=EMMC*2.04E-10*FRADQ(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
     2*XNPHX 
      THRMDXM=THS(NY,NX)*FRADQ(NZ,NY,NX)*XNPHX
      THRMDYM=DTHRMT(NY,NX)*FLAIQ(NZ,NY,NX)*XNPHX
      NNM=0
      RFLXCD=0.0
      EFLXCD=0.0
      SFLXCD=0.0
      HFLXCD=0.0
      VFLXCD=0.0
      THRMCD=0.0
      EVAPCD=0.0
      EPCD=0.0
      TEVDG=0.0
      TSHDG=0.0
      PAREY=FRADQ(NZ,NY,NX)*PAREX(NY,NX)
      PARSY=FRADQ(NZ,NY,NX)*PARSX(NY,NX)
C
C     CONVERGENCE SOLUTION FOR STANDING DEAD ENERGY BALANCE
C
C     VHCPQ=canopy heat capacity used to calculate canopy air
C        temperature from ‘hour1.f’ (MJ K-1)
C     FLAIQ=fraction of total canopy+standing dead radiation
C        received by each PFT standing dead from ‘hour1.f’
C     TKQY,TKDY=standing dead air,surface temperature (K) 
C     HCBFCY=heat released by canopy+standing dead combustion 
C        in previous time step (MJ t-1)
C     VOLWQ=water volume on standing dead surfaces (m3)
C     FLWCDM,HFLWCDM=standing dead water,heat retention 
C        of precipitation (m3 t-1,MJ t-1) 
C
      DO 5005 M=1,NPH
      IC=0
      XC=1.0
      IF(VHCPQ(NY,NX).GT.ZEROS(NY,NX)
     2.AND.FLAIQ(NZ,NY,NX).GT.ZERO2)THEN
      TKQY=TKQY+HCBFDY(NZ,NY,NX)*XNPH
     2/(VHCPQ(NY,NX)*FLAIQ(NZ,NY,NX)) 
      ENDIF
      VOLWQ(NZ,NY,NX)=VOLWQ(NZ,NY,NX)+FLWCDM 
      HFLWCDM=FLWCDM*4.19*TKDY 
      DO 5000 NN=1,MXN
C
C     NET RADIATION FROM ABSORBED SW AND NET LW
C
C     THRMDZM=LW radiation emitted by standing dead (MJ t-1)
C     EMMDX=EMMC x SB constant (MJ h-1 m-2 K-4) x FRADQ x AREA (m2) 
C        x time step (MJ t-1 K-4)
C     TKDY=standing dead surface temperature (K)
C     THRMDGM=LW radiation emitted from standing dead 
C        to ground surface from ‘watsub.f’ (MJ t-1)
C     DTHRMD=net LW radiation absorbed by standing dead (MJ t-1) 
C     RADCD=total SW radiation absorbed by standing dead (MJ t-1) 
C     RFLXCDM=net SW+LW radiation absorbed by standing dead (MJ t-1)
C
      THRMDZM=EMMDX*TKDY**4 
      THRMDGM=EMMDX*(TKDY**4-TKGS(M,NY,NX)**4)
     2*FRADQ(NZ,NY,NX)
      DTHRMD=THRMDXM+THRMDYM-THRMDZM-THRMDGM 
      RFLXCDM=RADCD+DTHRMD
C     IF(NZ.EQ.2.AND.ICHKF.EQ.1)THEN
C     WRITE(*,4447)'THRMDG',I,J,NFZ,M,NN,NX,NY,NZ
C    2,RFLXCDM,RADCD,DTHRMD,THRMDXM,THRMDYM,THRMDZM,THRMDGM
C    3,TKDY,FRADQ(NZ,NY,NX),FLAIQ(NZ,NY,NX),THRMGD(M,NZ,NY,NX) 
4447  FORMAT(A8,8I4,40E14.6)
C     ENDIF
C
C     STANDING DEAD BOUNDARY LAYER RESISTANCE 
C
C     DTKQD=air-standing dead temperature difference (K)
C     TKAM,TKQY,TKDY=air,standing dead air,surface temperature (K)
C     RI=Richardson number
C     RIX,RIY=minimum, maximum values used to calculate Richardson’s
C        number from ‘starts.f’
C     RABD=standing dead boundary layer resistance (h m-1) 
C     RABX=biome standing dead isothermal boundary layer resistance
C        from ‘hour1.f’ (h m-1)
C     RABM,RABX=minimum, maximum values used to calculate  boundary
C        layer resistances from ‘starts.f’ (h m-1)
C     RAGD=standing dead aerodynamic resistance below canopy height
C     ZG,ZT,ZR=PFT standing dead,biome,surface roughness height 
C     RACD=standing dead aerodynamic resistance between canopy height
C        and maximum standing dead height (h m-1) 
C     RACG,RACGX=standing dead aerodynamic resistance below 
C        maximum standing dead height (h m-1)
C     RATD=standing dead boundary layer+aerodynamic resistance (h m-1)
C     RAHD=standing dead surface resistance adjusted for standing dead
C        aerodynamic resistance (h m-1)
C     RAH,RAE=surface isothermal resistance to sensible,latent heat 
C        (h m-1)
C     DTKD=difference between standing dead,surface air temperature (K)
C     PARSD,PARED=standing dead surface conductance to sensible,
C        latent heat (m3 t-1,MJ K-1 t-1)
C
      DTKQD=TKAM(NY,NX)-TKQY 
      RI=AMAX1(RIX,AMIN1(RIY
     2,RIBX(NY,NX)/TKAM(NY,NX)*DTKQD))
      RID=1.0-10.0*RI 
      RABD=AMIN1(RABZ,AMAX1(RABM,RABX(NY,NX)/RID))
      RACD=AMAX1(0.0,RACG(M,NY,NX)-RAGD(M,NZ,NY,NX))
      RATD=RABD+RACD 
      PARSGD=PARSX(NY,NX)/RAGD(M,NZ,NY,NX)*FRADT(NY,NX)
      PAREGD=PAREX(NY,NX)/RAGD(M,NZ,NY,NX)*FRADT(NY,NX)
      DTKD=TKQY-TKDY
      RI=AMAX1(RIX,AMIN1(RIY
     2,RIBX(NY,NX)/TKQY*DTKD))
      RIS=1.0-3.2*RI 
      RAHD=AMIN1(RABZ,AMAX1(RABM,RAH/RIS)) 
      PARSD=PARSY/RAHD
      PARED=PAREY/(RAHD+RAE)  
C
C     STANDING DEAD VAPOR CONCENTRATION AND EVAPORATION OF INTERCEPTED 
C     WATER OR TRANSPIRATION OF UPTAKEN WATER
C
C     TKDY,VPDY=standing dead surface temperature,vapor concentration
C        (K,m3 m-3)
C     EX=standing dead-atmosphere water flux (m3 t-1)
C     PARSD,PARED=standing dead surface conductance to sensible,
C        latent heat (m3 t-1,MJ K-1 t-1)
C     VPQY=standing dead vapor concentration (m3 m-3)
C     EVAPCDM=standing dead surface evaporation (m3 t-1)
C     VOLWQ=water volume on standing dead surfaces (m3)
C     EPCDM=standing dead surface transpiration (m3 t-1) 
C     EFLXCDM=standing dead surface latent heat flux (MJ t-1)
C     VFLXCDM=convective heat flux from EFLXCDM (MJ t-1)
C     VAP=latent heat of evaporation from ‘starts.f’ (MJ m-3)
C     DTKD=difference between standing dead,surface air temperature (K)
C     SFLXCDM=standing dead surface sensible heat flux (MJ t-1)
C     XNPG=time step for heat fluxes from ‘wthr.f’ (h t-1)
C
      VPDY=2.173E-03/TKDY
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKDY))
      EX=PARED*(VPQY-VPDY)
      IF(EX.GT.0.0)THEN
      EVAPCDM=EX
      EX=0.0
      ELSE 
      EVAPCDM=AMAX1(EX,-AMAX1(0.0,VOLWQ(NZ,NY,NX)*XNPG))
      EX=EX-EVAPCDM
      ENDIF
      EPCDM=0.0
      EFLXCDM=(EPCDM+EVAPCDM)*VAP
      VFLXCDM=EVAPCDM*4.19*TKDY
      SFLXCDM=PARSD*DTKD
C
C     STANDING DEAD SENSIBLE+STORAGE HEAT FROM RN,LE AND CONVECTIVE
C     HEAT FLUXES
C
C     HFLXCDM=standing dead storage heat flux (MJ t-1)  
C     RFLXCDM=net SW+LW radiation absorbed by standing dead (MJ t-1) 
C     EFLXCDM=standing dead surface latent heat flux (MJ t-1) 
C     SFLXCDM=standing dead surface sensible heat flux (MJ t-1) 
C     VFLXCDM=convective heat flux from EFLXCDM (MJ t-1) 
C     HFLWCDM=convective heat flux from precipitation to standing dead
C        (MJ t-1) 
C 
      HFLXCDM=RFLXCDM+EFLXCDM+SFLXCDM+VFLXCDM+HFLWCDM
C
C     SOLVE FOR STANDING DEAD TEMPERATURE CAUSED BY SENSIBLE+STORAGE
C     HEAT
C
C     VHCPDX=wet standing dead heat capacity (MJ K-1)
C     VHCPDC=standing dead heat capacity (MJ K-1)
C     EVAPCDM=standing dead surface evaporation (m3 t-1)
C     FLWCDM=standing dead water retention of precipitation (m3 t-1) 
C     TKDY=standing dead surface temperature (K)
C     HFLXCDM=standing dead storage heat flux (MJ K-1) 
C     DTKX=step change in standing dead surface temperature used in
C        convergence (K) 
C     XC,IC=magnitude,direction of change in standing dead surface 
C        temperature for next cycle 
C
      VHCPDC=VHCPDX+4.19*(EVAPCDM+FLWCDM)
      TKDX=TKDY
      IF(VHCPDC.GT.VHCPYZ)THEN
      TKDZ=(TKDX*VHCPDX+HFLXCDM)/VHCPDC
      TKDY=TKDX+XC*DTKX*(TKDZ-TKDX)
      ELSE
      TKDZ=TKDX
      TKDY=TKDX
      ENDIF
      IF((IC.EQ.0.AND.TKDY.GT.TKDX).OR.(IC.EQ.1.AND.TKDY.LT.TKDX))THEN
      XC=0.5*XC
      ENDIF
      IF(TKDY.GT.TKDX)THEN
      IC=1
      ELSE
      IC=0
      ENDIF
C     IF(NZ.EQ.1.AND.ICHKF.EQ.1)THEN
C     WRITE(*,4446)'TKDY',I,J,NFZ,M,NN,NX,NY,NZ 
C    2,XC,TKDX,TKDY,TKDY-TKDX,DTKX,TKD(NZ,NY,NX),TKQG(M,NY,NX)
C    2,TKQD(NZ,NY,NX),DTKD,VHCPDC,VHCPYZ,FRADQ(NZ,NY,NX),EX
C    3,EPCDM,EVAPCDM,FLWCDM,PARED,PARSD,RAHD,RABD,RACD,RAH,RAE 
C    2,VOLWQ(NZ,NY,NX),ARSTG(NZ,NY,NX),RADD(NZ,NY,NX)
C    3,RADCD,DTHRMD,THRMDXM,THRMDYM,THRMDZM,THRMDGM,TKGS(M,NY,NX) 
C    3,HFLXCDM,RFLXCDM,EFLXCDM,SFLXCDM,VFLXCDM,HFLWCDM
C    5,VPQD(NZ,NY,NX) 
C    3,HFLXCD,RFLXCD,EFLXCD,SFLXCD,VFLXCD,HFLWCD 
C    4,HCBFDY(NZ,NY,NX)
4446  FORMAT(A8,8I4,60E14.6)
C     ENDIF
      NNM=MAX(NN,NNM)
C
C     IF CONVERGENCE FOR ENERGY EXCHANGE IS ACHIEVED
C
      IF(ABS(TKDY-TKDX).LT.DIFFT)GO TO 5500
5000  CONTINUE
5500  CONTINUE
C
C     CALCULATE CANOPY AIR TEMPERATURE, VAPOR CONCENTRATION
C
C     VHCPDP=dry standing dead heat capacity (MJ K-1)
C     SFLXA,EVAPA=standing dead-atmosphere sensible heat flux,
C        vapor flux (MJ t-1,m3 t-1)
C     DTKQD=atmosphere-standing dead air temperature difference (K)
C     PAREY,PARSY=standing dead boundary layer conductance for
C        latent,sensible heat (m2 t-1,MJ m-1 K-1 t-1)
C     RATD=standing dead boundary layer+aerodynamic resistance (h t-1)
C     VHCPQ,EVAPQ=canopy heat capacity,volume used to calculate 
C        canopy air temperature,vapor concentration from ‘hour1.f’ 
C        (MJ K-1,m3)
C     FLAIQ=fraction of total canopy+standing dead radiation
C        received by each PFT standing dead from ‘hour1.f’
C     VPAM,VPQY=atmosphere,standing dead vapor concentration (m3 m-3)
C     TSHDNM,TEVDNM=standing dead net sensible heat flux,vapor flux
C        (MJ t-1,m3 t-1)
C     SFLXCDM=standing dead surface sensible heat flux (MJ t-1)
C     EPCDM,EVAPCDM=standing dead surface transpiration,evaporation 
C        (m3 t-1)
C     TSHDGM,TEVDGM=ground surface-standing dead sensible heat flux, 
C        vapor flux from ‘watsub.f’ (MJ t-1,m3 t-1)
C     TSHDYM,TEVDYM=total lateral standing dead sensible heat flux,
C        vapor flux (MJ t-1,m3 t-1)
C     TKQY=standing dead air temperature (K)
C     VPQY,VPSY=standing dead,saturated vapor concentration (m3 m-3)
C     VOLWQ=water volume on canopy surfaces (m3)
C     XNPHX=time step for heat fluxes from ‘wthr.f’ (h t-1)
C
      IF(VHCPDP.GT.VHCPXZ.AND.FLAIQ(NZ,NY,NX).GT.1.0E-03)THEN
      SFLXA=DTKQD*PARSY/RATD*FLAIQ(NZ,NY,NX)
      DVPQA=VPAM(NY,NX)-VPQY
      EVAPA=DVPQA*AMIN1(PAREY/RATD
     2,EVAPQ(NY,NX)*FLAIQ(NZ,NY,NX))
      TSHDNM=SFLXA-SFLXCDM
      TEVDNM=EVAPA-EPCDM-EVAPCDM 
      TSHDGM=PARSGD*(TKQY-TKQG(M,NY,NX))*FLAIQ(NZ,NY,NX)
      TEVDGM=PAREGD*(VPQY-VPQG(M,NY,NX))*FLAIQ(NZ,NY,NX)
      TSHDYM=DSHQT(NY,NX)*FLAIQ(NZ,NY,NX)*XNPHX 
      TEVDYM=DVPQT(NY,NX)*FLAIQ(NZ,NY,NX)*XNPHX
      TKQY=TKQY+(TSHDNM-TSHDGM-TSHDYM)
     2/(VHCPQ(NY,NX)*FLAIQ(NZ,NY,NX)) 
      VPSY=2.173E-03/TKQY
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKQY))
      VPQY=AMAX1(0.0,AMIN1(VPSY,VPQY+(TEVDNM-TEVDGM-TEVDYM)
     2/(EVAPQ(NY,NX)*FLAIQ(NZ,NY,NX))))
      ELSE
      TSHDNM=0.0
      TEVDNM=0.0 
      TSHDGM=0.0
      TEVDGM=0.0 
      TSHDYM=0.0 
      TEVDYM=0.0
      TKQY=TKAM(NY,NX)
      VPQY=VPAM(NY,NX)
      ENDIF
      VOLWQ(NZ,NY,NX)=VOLWQ(NZ,NY,NX)+EVAPCDM
C
C     AGGREGATE STANDING DEAD ENERGY BALANCE FOR ALL M
C
C     RFLXCD=net SW+LW radiation absorbed by standing dead (MJ t-1)
C     EFLXCD=standing dead surface latent heat flux (MJ t-1)
C     SFLXCD=standing dead surface sensible heat flux (MJ t-1)
C     HFLXCD=standing dead storage heat flux (MJ t-1)
C     VFLXCD=convective heat flux from EFLXCD (MJ t-1)
C     EPCD,EVAPCD=standing dead surface transpiration,evaporation 
C        (m3 t-1)
C     TSHDG,TEVDG=ground surface-standing dead sensible heat flux, 
C        vapor flux from ‘watsub.f’ (MJ t-1,m3 t-1)
C     VHCPDX,VHCPDP=wet,dry standing dead heat capacity (MJ K-1)
C
      RFLXCD=RFLXCD+RFLXCDM
      EFLXCD=EFLXCD+EFLXCDM
      SFLXCD=SFLXCD+SFLXCDM
      HFLXCD=HFLXCD+HFLXCDM
      VFLXCD=VFLXCD+VFLXCDM
      THRMCD=THRMCD+THRMDGM
      EVAPCD=EVAPCD+EVAPCDM
      EPCD=EPCD+EPCDM
      TEVDG=TEVDG+TEVDGM
      TSHDG=TSHDG+TSHDGM
      VHCPDX=VHCPDP+4.19*AMAX1(0.0,VOLWQ(NZ,NY,NX))
C     IF(NZ.EQ.3)THEN
C     WRITE(*,4448)'TKQD',I,J,NFZ,M,NN,NX,NY,NZ,TKQD(NZ,NY,NX)
C    2,TKQY,TKAM(NY,NX),TKDY,TKD(NZ,NY,NX),TKQG(M,NY,NX),VHCPDP,VHCPXZ 
C    3,TSHDNM,TSHDGM,TSHDYM,SFLXCDM,SFLXA,PARSGD,RACG(M,NY,NX),PARSY
C    3,RATD,RABD,RACD,RAGD(M,NZ,NY,NX),RID,RIS,ZG(NZ,NY,NX),ZT(NY,NX) 
C    3,VHCPQ(NY,NX),FLAIQ(NZ,NY,NX),VHCPQ(NY,NX)*FLAIQ(NZ,NY,NX) 
C    3,HCBFDY(NZ,NY,NX)*XNPH,WTSTG(NZ,NY,NX)
C    4,ARSTG(NZ,NY,NX),DTKQD,PARSY/RATD 
C    3,VPQY,VPAM(NY,NX),VPSY,VPQD(NZ,NY,NX),TEVDNM,EPCDM,EVAPCDM 
C    4,EVAPA,EVAPQ(NY,NX)*FLAIQ(NZ,NY,NX),PAREY/RATD,TEVDN,EVAPCD
C    5,DVPQA,PAREY,FRADQ(NZ,NY,NX) 
C    2,WTSTG(NZ,NY,NX),ZT(NY,NX),ARLSS(NY,NX)
C    4,FRADP(NZ,NY,NX),FRADQ(NZ,NY,NX) 
C    4,VOLWQ(NZ,NY,NX),FLWCDM,TSHDN 
C    3,RADD(NZ,NY,NX),RADQ(NZ,NY,NX)
4448  FORMAT(A8,8I4,60E14.6)
C     ENDIF
5005  CONTINUE
C
C     FINAL STANDING DEAD TEMPERATURE 
C
C     TKQD,VPQD=final estimate of standing dead air temperature,
C        vapor concentration (K,m3 m-3)
C     TKD=final estimate of standing dead surface temperature (K)
C
      TKQD(NZ,NY,NX)=TKQY
      VPQD(NZ,NY,NX)=VPQY
      TKD(NZ,NY,NX)=TKDY
      TCD(NZ,NY,NX)=TKD(NZ,NY,NX)-273.15
      ELSE
C
C     VOLWQ=water volume on standing dead surfaces (m3)
C
      VOLWQ(NZ,NY,NX)=VOLWQ(NZ,NY,NX)+FLWD(NZ,NY,NX)*XNFH
      RFLXCD=0.0
      EFLXCD=0.0
      SFLXCD=0.0
      HFLXCD=0.0
      VFLXCD=0.0
      THRMCD=0.0 
      EVAPCD=0.0
      EPCD=0.0
      TEVDG=0.0
      TSHDG=0.0
      ENDIF
      ELSE
      VOLWQ(NZ,NY,NX)=VOLWQ(NZ,NY,NX)+FLWD(NZ,NY,NX)*XNFH
      RFLXCD=0.0
      EFLXCD=0.0
      SFLXCD=0.0
      HFLXCD=0.0
      VFLXCD=0.0
      THRMCD=0.0 
      EVAPCD=0.0
      EPCD=0.0
      TEVDG=0.0
      TSHDG=0.0
      TKQD(NZ,NY,NX)=TKAM(NY,NX)
      VPQD(NZ,NY,NX)=VPAM(NY,NX)
      IF(ZG(NZ,NY,NX).GE.DPTHS(NY,NX)-ZERO)THEN
      TKD(NZ,NY,NX)=TKAM(NY,NX)
      ELSE
      TKD(NZ,NY,NX)=TKW(1,NY,NX)
      ENDIF
      TCD(NZ,NY,NX)=TKD(NZ,NY,NX)-273.15
      ENDIF
C
C     TOTAL CANOPY + STANDING DEAD ENERGY BALANCE
C
C     RFLXC=net SW+LW radiation absorbed by canopy+standing dead 
C        (MJ t-1)
C     EFLXC=canopy+standing dead surface latent heat flux (MJ t-1)
C     SFLXC=canopy+standing dead surface sensible heat flux (MJ t-1)
C     HFLXC=canopy+standing dead storage heat flux (MJ t-1)
C     VFLXC=convective heat flux from EFLXC (MJ t-1)
C     EP,EVAPC=canopy+standing dead surface transpiration,evaporation
C        (m3 t-1)
C
      RFLXC(NZ,NY,NX)=RFLXCC+RFLXCD
      EFLXC(NZ,NY,NX)=EFLXCC+EFLXCD 
      SFLXC(NZ,NY,NX)=SFLXCC+SFLXCD 
      HFLXC(NZ,NY,NX)=HFLXCC+HFLXCD 
      VFLXC(NZ,NY,NX)=VFLXCC+VFLXCD 
      THRMP(NZ,NY,NX)=THRMCC+THRMCD 
      EVAPC(NZ,NY,NX)=EVAPCC+EVAPCD
      EP(NZ,NY,NX)=EPCC+EPCD
C     IF(NZ.EQ.1)THEN
C     WRITE(*,4445)'TQCT',I,J,NFZ,NX,NY,NZ
C    2,RFLXC(NZ,NY,NX),RFLXCC,RFLXCD
C    3,EFLXC(NZ,NY,NX),EFLXCC,EFLXCD 
C    4,SFLXC(NZ,NY,NX),SFLXCC,SFLXCD 
C    5,HFLXC(NZ,NY,NX),HFLXCC,HFLXCD 
C    6,THRMP(NZ,NY,NX),THRMCC,THRMCD 
C    7,EP(NZ,NY,NX),EPCC,EPCD
C    8,EVAPC(NZ,NY,NX),EVAPCC,EVAPCD
C    2,TKC(NZ,NY,NX),TKD(NZ,NY,NX) 
C    2,TKQC(NZ,NY,NX),TKQD(NZ,NY,NX)
C    2,VPQC(NZ,NY,NX),VPQD(NZ,NY,NX) 
C    3,VOLWP(NZ,NY,NX),FLWC(NZ,NY,NX),FLWCCM 
C    3,VOLWC(NZ,NY,NX),FLWC(NZ,NY,NX),FLWCCM 
C    3,VOLWQ(NZ,NY,NX),FLWD(NZ,NY,NX),FLWCDM 
C    4,TKQT(NY,NX),VPQT(NY,NX),FLAIP(NZ,NY,NX),FLAIQ(NZ,NY,NX)
C    4,HCBFCY(NZ,NY,NX),HCBFDY(NZ,NY,NX)
4445  FORMAT(A8,6I4,40E14.6)
C     ENDIF
9985  CONTINUE
9990  CONTINUE
9995  CONTINUE
      RETURN
      END

