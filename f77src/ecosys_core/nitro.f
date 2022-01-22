      SUBROUTINE nitro(I,J,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE CALCULATES ALL SOIL BIOLOGICAL TRANSFORMATIONS
C
      include "parameters.h"
      include "blkc.h"
      include "blk2a.h"
      include "blk2b.h"
      include "blk2c.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk10.h"
      include "blk11a.h"
      include "blk11b.h"
      include "blk13a.h"
      include "blk13b.h"
      include "blk13c.h"
      include "blk13d.h"
      include "blk15a.h"
      include "blk15b.h"
      include "blk16.h"
      include "blk18a.h"
      include "blk18b.h"
      include "blk19a.h"
      include "blk19b.h"
      include "blk21b.h"
      DIMENSION CNOMA(JG,7,0:5),CPOMA(JG,7,0:5),OMA(JG,7,0:5)
     2,FOMA(JG,7,0:5),FOMN(JG,7,0:5),RDOSC(4,0:4),RDOSN(4,0:4)
     3,RDOSP(4,0:4),RHOSC(4,0:4),RHOSN(4,0:4),RHOSP(4,0:4)
     4,RCOSC(4,0:4),RCOSN(4,0:4),RCOSP(4,0:4),SPOSC(4,0:4)
     5,RDORC(2,0:4),RDORN(2,0:4),RDORP(2,0:4),SPORC(2)
     6,RDOHC(0:4),RDOHN(0:4),RDOHP(0:4),RDOHA(0:4),CSORP(0:4),ZSORP(0:4)
     8,PSORP(0:4),CSORPA(0:4),OSRH(0:4),RUPOX(JG,7,0:5),RGN2F(JG,7,0:5)
     9,RGOMO(JG,7,0:5),ROXYM(JG,7,0:5),ROXYP(JG,7,0:5),ROXYO(JG,7,0:5)
     1,RDNO3(JG,7,0:5),RDNOB(JG,7,0:5),RDNO2(JG,7,0:5),RDN2B(JG,7,0:5)
     2,RDN2O(JG,7,0:5),RGOMD(JG,7,0:5),RMOMC(2,JG,7,0:5),RINH4(JG,7,0:5)
     3,RINO3(JG,7,0:5),RIPO4(JG,7,0:5),RINB4(JG,7,0:5),RINB3(JG,7,0:5)
     4,RIPOB(JG,7,0:5),FOMK(JG,7,0:5),RDOMC(2,JG,7,0:5)
     5,RDOMN(2,JG,7,0:5),RDOMP(2,JG,7,0:5),RHOMC(2,JG,7,0:5)
     6,RHOMN(2,JG,7,0:5),RHOMP(2,JG,7,0:5),RCOMC(2,JG,7,0:5)
     7,RCOMN(2,JG,7,0:5),RCOMP(2,JG,7,0:5),CGOMC(JG,7,0:5)
     8,CGOMN(JG,7,0:5),RH2GX(JG,7,0:5),CGOMP(JG,7,0:5)
     1,RDMMC(2,JG,7,0:5),RHMMC(2,JG,7,0:5),RCMMC(2,JG,7,0:5)
     2,RDMMN(2,JG,7,0:5),RHMMN(2,JG,7,0:5),RCMMN(2,JG,7,0:5)
     3,RDMMP(2,JG,7,0:5),RHMMP(2,JG,7,0:5),RCMMP(2,JG,7,0:5)
     4,RCCMC(2,JG,7,0:4),RCCMN(2,JG,7,0:4),RCCMP(2,JG,7,0:4)
     5,RN2FX(JG,7,0:5),TOMK(0:5)
     6,TONK(0:5),TOPK(0:5),SPOMC(2),OMC2(JG,7,0:5)
     7,TFNG(JG,7,0:5),TFNR(JG,7,0:5)
     2,OMN2(JG,7,0:5),FOM2(JG,7,0:5),FOCA(0:4),FOAA(0:4)
     3,RXOMC(2,JG,7,0:5),RXOMN(2,JG,7,0:5),RXOMP(2,JG,7,0:5)
     4,R3OMC(2,JG,7,0:5),R3OMN(2,JG,7,0:5)
     5,R3OMP(2,JG,7,0:5),RXMMC(2,JG,7,0:5),RXMMN(2,JG,7,0:5)
     6,RXMMP(2,JG,7,0:5),R3MMC(2,JG,7,0:5),R3MMN(2,JG,7,0:5)
     7,R3MMP(2,JG,7,0:5),WFN(JG,7,0:5)
      DIMENSION CGOQC(JG,7,0:5),CGOAC(JG,7,0:5),ROQCK(0:4),XOQCK(0:4)
     2,EN2F(7),ORCT(0:4),OSCT(0:4),OSAT(0:4),ZNH4T(0:JZ),ZNO3T(0:JZ)
     3,ZNO2T(0:JZ),H2P4T(0:JZ),RINH4R(JG,7,0:5),RINO3R(JG,7,0:5)
     4,RIPO4R(JG,7,0:5),FNH4XR(JG,7,0:5),FNO3XR(JG,7,0:5)
     5,FPO4XR(JG,7,0:5),RGOMY(JG,7,0:5),CNQ(0:4),CPQ(0:4)
     6,CNH(0:4),CPH(0:4),CNS(4,0:4),CPS(4,0:4),ROQCD(JG,7,0:4)
     7,FORC(0:5),SPOMK(2),RMOMK(2)
     8,CGOMS(2,JG,7,0:5),CGONS(2,JG,7,0:5),CGOPS(2,JG,7,0:5)
     9,H1P4T(0:JZ)
     1,TONX(0:5),TOPX(0:5),FCNK(0:4),FCPK(0:4),FP14XR(JG,7,0:5)
     2,RCO2X(JG,7,0:5),RCH3X(JG,7,0:5),RCH4X(JG,7,0:5)
     3,RVOXA(JG,7),RVOXB(JG,7)
     2,XOQCZ(0:4),XOQNZ(0:4),XOQPZ(0:4),XOQAZ(0:4)
     3,XOMCZ(3,JG,7,0:4),XOMNZ(3,JG,7,0:4),XOMPZ(3,JG,7,0:4)
     4,FCN(JG,7,0:5),FCP(JG,7,0:5),FCNP(JG,7,0:5),RIP14(JG,7,0:5)
     5,RIP1B(JG,7,0:5)
     5,TCGOQC(0:5),TCGOAC(0:5),TCGOMN(0:5),TCGOMP(0:5)
     6,TRN2ON(JY,JX),TRN2OD(JY,JX),TRN2GD(JY,JX),RIP14R(JG,7,0:5)

      DIMENSION ONL(4,0:4),OPL(4,0:4),EFIRE(2,21:22),DOSA(0:4)
      DIMENSION TOMCNK(2)
C
C     SUBSTRATE DECOMPOSITION BY MICROBIAL POPULATIONS
C
C     ORAD=microbial radius (m), BIOS=microbial density (n m-3)
C     BIOA=microbial surface area (m2 m-3), DCKI=inhibition of
C     decomposition by microbial concentration (g C m-3)
C     RCCX=maximum remobilization of microbial N (-)
C     RCCY=maximum remobilization of microbial P (-)
C     RCCZ, RCCY = minimum, maximum remobilization of microbial C (-)
C     FPRIM, FPRIMM=fraction of nonstructural, microbial C,N,P
C     transferred with priming (-), OMGR=rate constant for
C     transferring nonstructural to structural microbial C (h-1)
C     OQKI=DOC product inhibition constant for decomposition (g C m-3)
C     H2KI=H2 product inhibition for methanogenesis (g H m-3)
C     COMKI, COMKM= Km to slow microbial decomposition, maintenance
C     respiration with low microbial C (g micr C g-1 subs C)
C     CKC=controls C remobilization of microbial C (g C g-1 C)
C     FOSCZ0, FOSCZL=rate constants for mixing surface (0) and
C     subsurface (L) litter (h-1),FMN=minimum ratio of total
C     biological demand for any substrate by any microbial population
C     DCKM0, DCKML=Km for SOC decomposition (g C g-1 soil)
C
      PARAMETER (ORAD=1.0E-06,BIOS=1.0E-06/(4.19*ORAD**3)
     2,BIOA=BIOS*12.57*ORAD**2,DCKI=2.5,RCCX=0.833
     3,RCCQ=0.833,RCCZ=0.167,RCCY=0.833,FPRIM=5.0E-02,FPRIMM=1.0E-06
     4,OMGR=2.5E-01,OQKI=1.2E+03,H2KI=1.0,OAKI=12.0,COMKI=1.0E-03
     5,COMKM=1.0E-04,CKC=1.0E-03,FOSCZ0=2.0E-02,FOSCZL=2.0E-06
     6,FMN=1.0E-03,DCKM0=1.0E+03,DCKML=1.0E+03)
C
C     SPECIFIC RESPIRATION RATES, M-M UPTAKE CONSTANTS,
C     STOICHIOMETRIC CONSTANTS FOR MICROBIAL REDOX REACTIONS
C
C     VMX*=specific oxidation rates (g C g-1C h-1_
C        O=all bacteria, F=fungi, M=acetotrophic methanogens
C        H=ammonia oxidizers, N=nitrite oxidizers, 4=methanotrophs
C        C=hydrogenotrophic methanogens
C     OQK*=Km for DOC uptake by heterotrophs (g C m-3)
C        M=all bacteria and fungi, A=acetate by fermenters
C        AM=acetate by acetotrophic methanogens
C     CCKM=Km for CO2 uptake, CCK4=Km for CH4 uptake (g C m-3)
C     Z*KM=Km for N uptake (g N m-3)
C        H=NH4 by nitrifiers, N=NO2 by nitrifiers
C        3=NO3 by denitrifiers, 2=NO2 by denitrifiers
C        1=N2O uptake by denitrifiers
C     Z4*=NH4 uptake kinetics by all MFTs(g N m-2 h-1, g N m-3)
C       MX=maximum uptake rate, KU=Km, MN= minimum concentration
C     ZO*=NO3 uptake kinetics by all MFTs(g N m-2 h-1, g N m-3)
C       MX=maximum uptake rate, KU=Km, MN= minimum concentration
C     HP*=H2PO4 uptake kinetics by all MFTs(g P m-2 h-1, g P m-3)
C       MX=maximum uptake rate, KU=Km, MN= minimum concentration
C     ZFKM=Km for N2 uptake by diazotrophs (g N m-3)
C     H2KM=Km for H2 uptake by hydrogenotrophic methanogens (g H m-3)
C     ECNH=efficiency CO2 conversion to biomass by ammonia oxidizers
C     ECNO=efficiency CO2 conversion to biomass by nitrite oxidizers
C     ECHO=efficiency CO2 conversion to biomass by methane oxidizers
C     ECN3,ECN2,ECN1=N2:O2 ratios for e- transfers to NO3, NO2 and N2O
C     by denitrifiers, RNFNI=parameter for nitrification inhibition
C     ZHKI=inhibition of nitrification inhibition by NH3 (g N m-3)
C     VMKI=product inhibn for NOx reduction by denitrifiers(g N m-3)
C     VHKI=product inhibn for NH3 oxidation by nitrifiers (g N m-3)
C     OXKA=Km for O2 uptake by nitrifiers(g O m-3)
C
      PARAMETER (VMXO=0.125,VMXF=0.125,VMXM=0.125,VMXH=0.375
     2,VMXN=0.25,VMX4=0.375,VMXC=0.125,OQKM=1.2E+01,OQKA=1.2E+01
     3,OQKAM=1.2E+01,CCKM=0.15,CCK4=1.2E-03,ZHKM=1.4
     4,ZNKM=1.4,Z3KM=1.4,Z2KM=1.4,Z1KM=0.014,Z4MX=5.0E-03
     5,Z4KU=0.40,Z4MN=0.0125,ZOMX=5.0E-03,ZOKU=0.35,ZOMN=0.03
     7,HPMX=1.0E-03,HPKU=0.075,HPMN=0.002,ZFKM=0.14,H2KM=0.01
     8,ECNH=0.30,ECNO=0.10,ECHO=0.75,ECN3=0.857,ECN2=0.857,ECN1=0.429
     9,RNFNI=2.0E-04,ZHKI=7.0E+03,VMKI=0.25,VHKI=15.0,OXKA=0.16)
C
C     ENERGY REQUIREMENTS FOR MICROBIAL GROWTH AND
C     ENERGY YIELDS FROM REDUCTION OF O2, OC, CH4, NO3, N2
C
C     EOM*=energy requirements for microbial growth (kJ g-1 C)
C        C=aerobic bacteria, D=denitrifiers, G=fungi, F=fermenters
C        H=methanogens, N=diazotrophs
C     G*=free energy yields of redox reactions (kJ g-1 C or N)
C        O2X=DOC-CO2, H4X=CO2-CH4, CHX=DOC-acetate, O2A=acetate-CO2
C        C4X=acetate-CH4, COX=CO2-CH4, NOX=NO3-NO2,NO2-N2O,N2O-N2
C        N2X=N2-NH3
C     E*=growth respiration efficiency (-)(growth yield=1.0-E*)
C        N2X=aerobic N2 fixation, N2Y=anaerobic N2 fixation
C        O2X=aerobic bacteria (DOC), H4X=fermenters, O2G=fungi
C        O2D=denitrifiers (aerobic), NFX=diazotrophs
C        NOX= denitrifiers (anaerobic),O2A=aerobic bacteria (acetate)
C
      PARAMETER (EOMC=25.0,EOMD=37.5,EOMG=37.5,EOMF=75.0,EOMH=25.0
     2,EOMN=75.0,GO2X=37.5,GH4X=66.5,GCHX=4.50
     3,GO2A=GO2X-GCHX,GC4X=3.00,GCOX=11.00,GNOX=10.0
     3,GN2X=187.5,EN2X=GO2X/GN2X,EN2Y=GCHX/GN2X
     4,EO2X=1.0/(1.0+GO2X/EOMC),EH4X=1.0/(1.0+GH4X/EOMC)
     5,EO2G=1.0/(1.0+GO2X/EOMG),EO2D=1.0/(1.0+GO2X/EOMD)
     6,ENFX=1.0/(1.0+GO2X/EOMN),ENOX=1.0/(1.0+GNOX/EOMC)
     7,EO2A=1.0/(1.0+GO2A/EOMC))
C
C     SORPTION COEFFICIENTS
C
      PARAMETER (TSORP=0.5,HSORP=1.0)
C
C     SPECIFIC DECOMPOSITION RATES
C
C     DOSA=rate constant for litter colonization by heterotrophs (g C g-1 C)
C     SP*= specific decomposition rate constant (g subs. C g-1 micr. C)
C       OHC=adsorbed SOC, OHA=adsorbed acetate, OSC=SOC
C       (K=0,M=1,4 woody litter, K=1,M=1,4 non-woody litter,
C       K=2,M=1,4 manure, K=3,M=1,1 POC, K=4,M=1,2 humus)
C       ORC (M=1,2) microbial residue, OMC (M=1,2) microbial biomass
C     RMOM=specific maintenance respiration (g C g-1 N h-1)
C
      PARAMETER (SPOHC=0.25,SPOHA=0.25,RMOM=0.010)
      DATA DOSA/0.25E-03,0.25,0.25,0.25,0.25/
      DATA SPOSC/7.5,7.5,1.5,0.5,7.5,7.5,1.5,0.5
     2,7.5,7.5,1.5,0.5,0.05,0.00,0.00,0.00
     3,0.05,0.0167,0.00,0.00/
      DATA SPORC/7.5,1.5/
      DATA SPOMC/1.0E-02,0.1E-02/
      DATA EN2F/0.0,0.0,0.0,0.0,0.0,EN2X,EN2Y/
      DATA EFIRE/1.0,1.0,0.917,0.167/
      REAL*4 WFNG,TFNX,TFNY,TFNG,TFNR,CNSHZ,CPSHZ,FRM


      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
C     IF(I.EQ.1.AND.J.EQ.1)THEN
C     TRN2ON(NY,NX)=0.0
C     TRN2OD(NY,NX)=0.0
C     TRN2GD(NY,NX)=0.0
C     ENDIF
C
C     VOLWZ=water volume used to calculate aqueous microbial
C     concentrations that drive microbial density effects on
C     decomposition
C
      DO 998 L=0,NL(NY,NX)
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IF(L.EQ.0.OR.L.GE.NU(NY,NX))THEN
      IF(L.EQ.0)THEN
      KL=2
      IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETR=VOLW(0,NY,NX)/VOLR(NY,NX)
      THETZ=AMAX1(0.0,THETR-THETY(L,NY,NX))
      VOLWZ=THETZ*VOLR(NY,NX)
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
C     WRITE(*,8825)'THETZ',I,J,L,THETR,THETZ,VOLWZ,VOLWRX(NY,NX)
C    2,VOLW(0,NY,NX),POROS(L,NY,NX),FC(0,NY,NX),WP(0,NY,NX)
C    3,THETY(L,NY,NX),PSISM(0,NY,NX),ORGC(0,NY,NX),VOLR(NY,NX)
8825  FORMAT(A8,3I4,20E12.4)
C     ENDIF
      ELSE
      VOLWZ=0.0
      ENDIF
      ELSE
      KL=4
      THETZ=AMAX1(0.0,(AMIN1(AMAX1(0.5*POROS(L,NY,NX),FC(L,NY,NX))
     2,THETW(L,NY,NX))-THETY(L,NY,NX)))
      VOLWZ=THETZ*VOLY(L,NY,NX)
C     IF((I/120)*120.EQ.I.AND.J.EQ.24.AND.L.LE.6)THEN
C     WRITE(*,8824)'THETZ',I,J,NX,NY,L,THETZ,THETW(L,NY,NX),VOLWZ
C    2,POROS(L,NY,NX),FC(L,NY,NX),WP(L,NY,NX),THETY(L,NY,NX)
C    3,VOLW(L,NY,NX),VOLWH(L,NY,NX),VOLY(L,NY,NX),VOLT(L,NY,NX)
C    4,DTBLX(NY,NX)
8824  FORMAT(A8,5I4,20E12.4)
C     ENDIF
      ENDIF
C
C     TEMPERATURE FUNCTIONS FOR GROWTH AND MAINTENANCE
C     WITH OFFSET FOR THERMAL ADAPTATION
C
C     TKS=soil temperature
C     OFFSET=adjustment for acclimation based on MAT in starts.f
C     8.313,710.0=gas constant,enthalpy
C     62500=activation energy
C     197500,195000 low temp inactivation for growth,maintenance
C     222500,232500 high temp inactivation for growth,maintenance
C     TFNX,TFNY=temperature function for growth,maintenance respiration
C
      TKSO=TKS(L,NY,NX)+OFFSET(NY,NX)
      RTK=8.3143*TKSO
      STK=710.0*TKSO
      ACTV=1+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFNX=EXP(25.229-62500/RTK)/ACTV
      ACTVM=1+EXP((195000-STK)/RTK)+EXP((STK-232500)/RTK)
      TFNY=EXP(25.214-62500/RTK)/ACTVM
C
C     OXYI=inhibition of fermenters by O2
C     ORGCL=SOC used to calculate microbial concentration
C
      OXYI=1.0-1.0/(1.0+EXP(1.0*(-COXYS(L,NY,NX)+2.5)))
      ORGCL=AMIN1(1.0E+05*BKVL(L,NY,NX),ORGC(L,NY,NX))
C
C     TOTAL MINERAL NH4, NO3 AND PO4
C
C     allocate NH4, NO3, HPO4, H2PO4 to non-band and band fractions
C
      ZNH4T(L)=AMAX1(0.0,ZNH4S(L,NY,NX))+AMAX1(0.0,ZNH4B(L,NY,NX))
C     IF(ZNH4T(L).GT.ZEROS(NY,NX))THEN
C     FNH4S=AMAX1(0.0,ZNH4S(L,NY,NX))/ZNH4T(L)
C     FNHBS=AMAX1(0.0,ZNH4B(L,NY,NX))/ZNH4T(L)
C     ELSE
      FNH4S=VLNH4(L,NY,NX)
      FNHBS=VLNHB(L,NY,NX)
C     ENDIF
      ZNO3T(L)=AMAX1(0.0,ZNO3S(L,NY,NX))+AMAX1(0.0,ZNO3B(L,NY,NX))
C     IF(ZNO3T(L).GT.ZEROS(NY,NX))THEN
C     FNO3S=AMAX1(0.0,ZNO3S(L,NY,NX))/ZNO3T(L)
C     FNO3B=AMAX1(0.0,ZNO3B(L,NY,NX))/ZNO3T(L)
C     ELSE
      FNO3S=VLNO3(L,NY,NX)
      FNO3B=VLNOB(L,NY,NX)
C     ENDIF
      ZNO2T(L)=AMAX1(0.0,ZNO2S(L,NY,NX))+AMAX1(0.0,ZNO2B(L,NY,NX))
C     IF(ZNO2T(L).GT.ZEROS(NY,NX))THEN
C     FNO2S=AMAX1(0.0,ZNO2S(L,NY,NX))/ZNO2T(L)
C     FNO2B=AMAX1(0.0,ZNO2B(L,NY,NX))/ZNO2T(L)
C     ELSE
      FNO2S=VLNO3(L,NY,NX)
      FNO2B=VLNOB(L,NY,NX)
C     ENDIF
      H1P4T(L)=AMAX1(0.0,H1PO4(L,NY,NX))+AMAX1(0.0,H1POB(L,NY,NX))
C     IF(H1P4T(L).GT.ZEROS(NY,NX))THEN
C     FH1PS=AMAX1(0.0,H1PO4(L,NY,NX))/H1P4T(L)
C     FH1PB=AMAX1(0.0,H1POB(L,NY,NX))/H1P4T(L)
C     ELSE
      FH1PS=VLPO4(L,NY,NX)
      FH1PB=VLPOB(L,NY,NX)
C     ENDIF
      H2P4T(L)=AMAX1(0.0,H2PO4(L,NY,NX))+AMAX1(0.0,H2POB(L,NY,NX))
C     IF(H2P4T(L).GT.ZEROS(NY,NX))THEN
C     FH2PS=AMAX1(0.0,H2PO4(L,NY,NX))/H2P4T(L)
C     FH2PB=AMAX1(0.0,H2POB(L,NY,NX))/H2P4T(L)
C     ELSE
      FH2PS=VLPO4(L,NY,NX)
      FH2PB=VLPOB(L,NY,NX)
C     ENDIF
C
C     CCO2S=aqueous CO2 concentration
C
      XCO2=CCO2S(L,NY,NX)/(CCO2S(L,NY,NX)+CCKM)
C
C     TOTAL SUBSTRATE
C
C     TOSC=total SOC, TOSA=total colonized SOC
C     TORC=total microbial residue, TOHC=total adsorbed C
C     in each K:
C     OSCT=total SOC n each K, OSAT=total colonized SOC
C     ORCT=total microbial residue, OHCT=total adsorbed C
C
      TOSC=0.0
      TOSA=0.0
      TORC=0.0
      TOHC=0.0
C
C     TOTAL SOLID SUBSTRATE
C
      DO 870 K=0,KL
      OSCT(K)=0.0
      OSAT(K)=0.0
      DO 865 M=1,4
      OSCT(K)=OSCT(K)+OSC(M,K,L,NY,NX)
      OSAT(K)=OSAT(K)+OSA(M,K,L,NY,NX)
865   CONTINUE
      TOSC=TOSC+OSCT(K)
      TOSA=TOSA+OSAT(K)
870   CONTINUE
C
C     TOTAL BIORESIDUE
C
      DO 880 K=0,KL
      ORCT(K)=0.0
      DO 875 M=1,2
      ORCT(K)=ORCT(K)+ORC(M,K,L,NY,NX)
C     IF(L.EQ.4.AND.K.EQ.2)THEN
C     WRITE(*,876)'ORCT',I,J,NX,NY,L,K,M,ORCT(K)
C    2,ORC(M,K,L,NY,NX)
876   FORMAT(A8,7I4,60E12.4)
C     ENDIF
875   CONTINUE
      TORC=TORC+ORCT(K)
C
C     TOTAL ADSORBED AND DISSOLVED SUBSTRATE
C
C     OSRH=total SOC
C
      TOHC=TOHC+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
880   CONTINUE
      DO 860 K=0,KL
      OSRH(K)=OSAT(K)+ORCT(K)+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
C     WRITE(*,861)'OSRH',I,J,NX,NY,L,K,OSRH(K),OSCT(K)
C    2,OSAT(K),ORCT(K),OHC(K,L,NY,NX),OHA(K,L,NY,NX)
861   FORMAT(A8,6I4,20E12.4)
C     ENDIF
860   CONTINUE
      TSRH=TOSA+TORC+TOHC
C
C     C:N AND C:P RATIOS OF TOTAL BIOMASS
C     CNOMA,CPOMA=N,P contents of active biomass OMA
C     FCN,FCP=effects of N,P limitations on biomass activity
C
      TOMA=0.0
      TOMN=0.0
      DO 890 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 895 N=1,7
      DO 895 NGL=1,JG
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      IF(OMC(1,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNOMA(NGL,N,K)=AMAX1(0.0,OMN(1,NGL,N,K,L,NY,NX)
     2/OMC(1,NGL,N,K,L,NY,NX))
      CPOMA(NGL,N,K)=AMAX1(0.0,OMP(1,NGL,N,K,L,NY,NX)
     2/OMC(1,NGL,N,K,L,NY,NX))
      ELSE
      CNOMA(NGL,N,K)=CNOMC(1,NGL,N,K)
      CPOMA(NGL,N,K)=CPOMC(1,NGL,N,K)
      ENDIF
      OMA(NGL,N,K)=AMAX1(0.0,OMC(1,NGL,N,K,L,NY,NX)/FL(1))
      FCN(NGL,N,K)=AMIN1(1.0,AMAX1(0.50,SQRT(CNOMA(NGL,N,K)
     2/CNOMC(1,NGL,N,K))))
      FCP(NGL,N,K)=AMIN1(1.0,AMAX1(0.50,SQRT(CPOMA(NGL,N,K)
     2/CPOMC(1,NGL,N,K))))
      FCNP(NGL,N,K)=AMIN1(FCN(NGL,N,K),FCP(NGL,N,K))
C
C     TOTAL BIOMASS
C     OMC2=active biomass in recalcitrant fraction
C
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      TOMA=TOMA+OMA(NGL,N,K)
      ENDIF
      IF((K.LE.4.AND.N.EQ.2).OR.(K.EQ.5.AND.N.EQ.1))THEN
      TOMN=TOMN+OMA(NGL,N,K)
      ENDIF
      OMC2(NGL,N,K)=AMAX1(0.0,AMIN1(OMA(NGL,N,K)*FL(2)
     2,OMC(2,NGL,N,K,L,NY,NX)))
      IF(OMC(2,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOM2(NGL,N,K)=AMAX1(0.0,OMC2(NGL,N,K)
     2/OMC(2,NGL,N,K,L,NY,NX))
      OMN2(NGL,N,K)=AMAX1(0.0,FOM2(NGL,N,K)
     2*OMN(2,NGL,N,K,L,NY,NX))
      ELSE
      FOM2(NGL,N,K)=0.0
      OMN2(NGL,N,K)=0.0
      ENDIF
      ENDIF
895   CONTINUE
      ENDIF
890   CONTINUE
      DO 690 K=0,KL
      TOMK(K)=0.0
      TONK(K)=0.0
      TOPK(K)=0.0
      TONX(K)=0.0
      TOPX(K)=0.0
      DO 685 N=1,7
      DO 685 NGL=1,JG
      TOMK(K)=TOMK(K)+OMA(NGL,N,K)
      TONK(K)=TONK(K)+OMA(NGL,N,K)*CNOMA(NGL,N,K)
      TOPK(K)=TOPK(K)+OMA(NGL,N,K)*CPOMA(NGL,N,K)
      TONX(K)=TONX(K)+OMA(NGL,N,K)*CNOMC(1,NGL,N,K)
      TOPX(K)=TOPX(K)+OMA(NGL,N,K)*CPOMC(1,NGL,N,K)
685   CONTINUE
690   CONTINUE
C
C     FOSRH=fraction of total SOC in each substrate complex K
C
      DO 790 K=0,KL
      IF(TSRH.GT.ZEROS(NY,NX))THEN
      FOSRH(K,L,NY,NX)=OSRH(K)/TSRH
      ELSE
      FOSRH(K,L,NY,NX)=1.0
      ENDIF
C
C     DOC CONCENTRATIONS
C
C     COQC,COQA=aqueous DOC,acetate concentrations
C     VOLWM=soil water content, FOSRH=fraction of total SOC
C     occupied by each substrate complex K
C
      IF(VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IF(FOSRH(K,L,NY,NX).GT.ZERO)THEN
      COQC(K,L,NY,NX)=AMAX1(0.0,OQC(K,L,NY,NX)
     2/(VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)))
      COQA(K,L,NY,NX)=AMAX1(0.0,OQA(K,L,NY,NX)
     2/(VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)))
      ELSE
      COQC(K,L,NY,NX)=AMAX1(0.0,OQC(K,L,NY,NX)/VOLWM(NPH,L,NY,NX))
      COQA(K,L,NY,NX)=AMAX1(0.0,OQA(K,L,NY,NX)/VOLWM(NPH,L,NY,NX))
      ENDIF
      ELSE
      COQC(K,L,NY,NX)=0.0
      COQA(K,L,NY,NX)=0.0
      OHCQ=0.0
      ENDIF
C
C     CNQ,CPQ=DON:DOC,DOP:DOC,FOCA,FOAA=DOC,DOA:(DOC+DOA)
C
      IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNQ(K)=AMAX1(0.0,OQN(K,L,NY,NX)/OQC(K,L,NY,NX))
      CPQ(K)=AMAX1(0.0,OQP(K,L,NY,NX)/OQC(K,L,NY,NX))
      ELSE
      CNQ(K)=0.0
      CPQ(K)=0.0
      ENDIF
      IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX).AND.OQA(K,L,NY,NX)
     2.GT.ZEROS(NY,NX))THEN
      FOCA(K)=OQC(K,L,NY,NX)/(OQC(K,L,NY,NX)+OQA(K,L,NY,NX))
      FOAA(K)=1.0-FOCA(K)
      ELSEIF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOCA(K)=1.0
      FOAA(K)=0.0
      ELSE
      FOCA(K)=0.0
      FOAA(K)=1.0
      ENDIF
790   CONTINUE
C
C     nitrous acid concn CHNO2 and energy yield of hydrogenotrophic
C     methanogenesis GH2X at ambient H2 concentration CH2GS
C
      CHY1=AMAX1(ZERO,10.0**(-(PH(L,NY,NX)-3.0)))
      CHNO2=CNO2S(L,NY,NX)*CHY1/0.5
      CHNOB=CNO2B(L,NY,NX)*CHY1/0.5
      GH2X=8.3143E-03*TKS(L,NY,NX)
     2*LOG((AMAX1(1.0E-03,CH2GS(L,NY,NX))/H2KI)**4)
C
C     RESPIRATION BY MICROBIAL POPULATIONS
C
      TFOXYX=0.0
      TFNH4X=0.0
      TFNO3X=0.0
      TFNO2X=0.0
      TFN2OX=0.0
      TFP14X=0.0
      TFPO4X=0.0
      TFNH4B=0.0
      TFNO3B=0.0
      TFNO2B=0.0
      TFP14B=0.0
      TFPO4B=0.0
      TCH4H=0.0
      TCH4A=0.0
      TFOQC=0.0
      TFOQA=0.0
      TRH2G=0.0
      IF(L.NE.0)THEN
      LL=L
      ELSE
      LL=NU(NY,NX)
      ENDIF
      DO 760 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      TCGOQC(K)=0.0
      TCGOAC(K)=0.0
      TCGOMN(K)=0.0
      TCGOMP(K)=0.0
      DO 750 N=1,7
      TOMCNK(:)=0.0
      DO NGL=1,JG
      DO M=1,2
      TOMCNK(M)=TOMCNK(M)+OMC(M,NGL,N,K,L,NY,NX)
      ENDDO
      ENDDO
      DO 750 NGL=1,JG
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      IF(K.LE.4)THEN
      IF(N.EQ.3)THEN
C
C     WFNG=water potential (PSISM) effect on microbial respiration
C     OXKX=Km for O2 uptake
C     OXKM=Km for heterotrophic O2 uptake set in starts.f
C     TFNG=combined temp and water stress effect on growth respiration
C     TFNR=temperature effect on maintenance respiration
C
      WFNG=EXP(0.1*PSISM(L,NY,NX))
      ELSE
      WFNG=EXP(0.2*PSISM(L,NY,NX))
      ENDIF
      OXKX=OXKM
      ELSE
      WFNG=EXP(0.2*PSISM(L,NY,NX))
      OXKX=OXKA
      ENDIF
      TFNG(NGL,N,K)=TFNX*WFNG
      TFNR(NGL,N,K)=TFNY
C
C     FOMA,FOMN=fraction of total active biomass C,N in each N and K
C
      IF(OMA(NGL,N,K).GT.0.0)THEN
      IF(TOMA.GT.ZEROS(NY,NX))THEN
      FOMA(NGL,N,K)=OMA(NGL,N,K)/TOMA
      ELSE
      FOMA(NGL,N,K)=1.0
      ENDIF
      IF(TOMN.GT.ZEROS(NY,NX))THEN
      FOMN(NGL,N,K)=OMA(NGL,N,K)/TOMN
      ELSE
      FOMN(NGL,N,K)=1.0
      ENDIF
      IF(TOMK(K).GT.ZEROS(NY,NX))THEN
      FOMK(NGL,N,K)=OMA(NGL,N,K)/TOMK(K)
      ELSE
      FOMK(NGL,N,K)=1.0
      ENDIF
C
C     ADJUST MCROBIAL GROWTH AND DECOMPOSITION RATES FOR BIOMASS
C
C     COMC=microbial C concentration relative to substrate
C     SPOMK=effect of microbial C concentration on microbial decay
C     RMOMK=effect of microbial C concentration on maintenance respn
C
      IF(ORGCL.GT.ZEROS(NY,NX))THEN
      DO 765 M=1,2
C      COMC=OMC(M,NGL,N,K,L,NY,NX)/ORGCL
      COMC=TOMCNK(M)/ORGCL
      SPOMK(M)=COMC/(COMC+COMKI)
      RMOMK(M)=COMC/(COMC+COMKM)
765   CONTINUE
      ELSE
      DO 770 M=1,2
      SPOMK(M)=1.0
      RMOMK(M)=1.0
770   CONTINUE
      ENDIF
C
C     FACTORS CONSTRAINING DOC, ACETATE, O2, NH4, NO3, PO4 UPTAKE
C     AMONG COMPETING MICROBIAL AND ROOT POPULATIONS IN SOIL LAYERS
C
C     F*=fraction of substrate uptake relative to total uptake from
C     previous hour. OXYX=O2, NH4X=NH4 non-band, NB4X=NH4 band
C     NO3X=NO3 non-band, NB3X=NO3 band, PO4X=H2PO4 non-band
C     POBX=H2PO4 band,P14X=HPO4 non-band, P1BX=HPO4 band, OQC=DOC
C     oxidation, OQA=acetate oxidation
C
      IF(ROXYY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOXYX=AMAX1(FMN,ROXYS(NGL,N,K,L,NY,NX)/ROXYY(L,NY,NX))
      ELSE
      FOXYX=AMAX1(FMN,FOMA(NGL,N,K))
      ENDIF
      IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4X=AMAX1(FMN,RINHO(NGL,N,K,L,NY,NX)/RNH4Y(L,NY,NX))
      ELSE
      FNH4X=AMAX1(FMN,FOMA(NGL,N,K)*VLNH4(L,NY,NX))
      ENDIF
      IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB4X=AMAX1(FMN,RINHB(NGL,N,K,L,NY,NX)/RNHBY(L,NY,NX))
      ELSE
      FNB4X=AMAX1(FMN,FOMA(NGL,N,K)*VLNHB(L,NY,NX))
      ENDIF
      IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3X=AMAX1(FMN,RINOO(NGL,N,K,L,NY,NX)/RNO3Y(L,NY,NX))
      ELSE
      FNO3X=AMAX1(FMN,FOMA(NGL,N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB3X=AMAX1(FMN,RINOB(NGL,N,K,L,NY,NX)/RN3BY(L,NY,NX))
      ELSE
      FNB3X=AMAX1(FMN,FOMA(NGL,N,K)*VLNOB(L,NY,NX))
      ENDIF
      IF(RPO4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FPO4X=AMAX1(FMN,RIPOO(NGL,N,K,L,NY,NX)/RPO4Y(L,NY,NX))
      ELSE
      FPO4X=AMAX1(FMN,FOMA(NGL,N,K)*VLPO4(L,NY,NX))
      ENDIF
      IF(RPOBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FPOBX=AMAX1(FMN,RIPBO(NGL,N,K,L,NY,NX)/RPOBY(L,NY,NX))
      ELSE
      FPOBX=AMAX1(FMN,FOMA(NGL,N,K)*VLPOB(L,NY,NX))
      ENDIF
      IF(RP14Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FP14X=AMAX1(FMN,RIPO1(NGL,N,K,L,NY,NX)/RP14Y(L,NY,NX))
      ELSE
      FP14X=AMAX1(FMN,FOMA(NGL,N,K)*VLPO4(L,NY,NX))
      ENDIF
      IF(RP1BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FP1BX=AMAX1(FMN,RIPB1(NGL,N,K,L,NY,NX)/RP1BY(L,NY,NX))
      ELSE
      FP1BX=AMAX1(FMN,FOMA(NGL,N,K)*VLPOB(L,NY,NX))
      ENDIF
      IF(K.LE.4)THEN
      IF(ROQCY(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOQC=AMAX1(FMN,ROQCS(NGL,N,K,L,NY,NX)/ROQCY(K,L,NY,NX))
      ELSE
      FOQC=AMAX1(FMN,FOMK(NGL,N,K))
      ENDIF
      TFOQC=TFOQC+FOQC
      IF(ROQAY(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOQA=AMAX1(FMN,ROQAS(NGL,N,K,L,NY,NX)/ROQAY(K,L,NY,NX))
      ELSE
      FOQA=AMAX1(FMN,FOMK(NGL,N,K))
      ENDIF
      TFOQA=TFOQA+FOQA
      ENDIF
      TFOXYX=TFOXYX+FOXYX
      TFNH4X=TFNH4X+FNH4X
      TFNO3X=TFNO3X+FNO3X
      TFPO4X=TFPO4X+FPO4X
      TFP14X=TFP14X+FP14X
      TFNH4B=TFNH4B+FNB4X
      TFNO3B=TFNO3B+FNB3X
      TFPO4B=TFPO4B+FPOBX
      TFP14B=TFP14B+FP1BX
C
C     FACTORS CONSTRAINING NH4, NO3, PO4 UPTAKE AMONG COMPETING
C     MICROBIAL POPULATIONS IN SURFACE RESIDUE
C     F*=fraction of substrate uptake relative to total uptake from
C     previous hour in surface litter, labels as for soil layers above
C
      IF(L.EQ.0)THEN
      IF(RNH4Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4XR(NGL,N,K)=AMAX1(FMN,RINHOR(NGL,N,K,NY,NX)
     2/RNH4Y(NU(NY,NX),NY,NX))
      ELSE
      FNH4XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
      ENDIF
      IF(RNO3Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3XR(NGL,N,K)=AMAX1(FMN,RINOOR(NGL,N,K,NY,NX)
     2/RNO3Y(NU(NY,NX),NY,NX))
      ELSE
      FNO3XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
      ENDIF
      IF(RPO4Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FPO4XR(NGL,N,K)=AMAX1(FMN,RIPOOR(NGL,N,K,NY,NX)
     2/RPO4Y(NU(NY,NX),NY,NX))
      ELSE
      FPO4XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
      ENDIF
      IF(RP14Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FP14XR(NGL,N,K)=AMAX1(FMN,RIPO1R(NGL,N,K,NY,NX)
     2/RP14Y(NU(NY,NX),NY,NX))
      ELSE
      FP14XR(NGL,N,K)=AMAX1(FMN,FOMK(NGL,N,K))
      ENDIF
      ENDIF
      IF(L.EQ.NU(NY,NX).AND.K.NE.3.AND.K.NE.4
     2.AND.BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
      TFNH4X=TFNH4X+FNH4XR(NGL,N,K)
      TFNO3X=TFNO3X+FNO3XR(NGL,N,K)
      TFPO4X=TFPO4X+FPO4XR(NGL,N,K)
      TFP14X=TFP14X+FP14XR(NGL,N,K)
      ENDIF
C
C     HETEROTROPHIC BIOMASS RESPIRATION
C
      IF(K.LE.4)THEN
C
C     RESPIRATION BY HETEROTROPHIC AEROBES:
C     N=(1)OBLIGATE AEROBES,(2)FACULTATIVE ANAEROBES,(3)FUNGI
C    (6)N2 FIXERS
C
      IF(N.LE.3.OR.N.EQ.6)THEN
C
C     ENERGY YIELDS OF O2 REDOX REACTIONS
C     E* = growth respiration efficiency calculated in PARAMETERS
C
      IF(N.EQ.1)THEN
      EO2Q=EO2X
      ELSEIF(N.EQ.2)THEN
      EO2Q=EO2D
      ELSEIF(N.EQ.3)THEN
      EO2Q=EO2G
      ELSEIF(N.EQ.6)THEN
      EO2Q=ENFX
      ENDIF
C
C     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
C     'RGO*Z'FROM SPECIFIC RESPIRATION RATE, ACTIVE BIOMASS, DOC OR
C     ACETATE CONCENTRATION,MICROBIAL C:N:P FACTOR, AND TEMPERATURE
C     FOLLOWED BY POTENTIAL RESPIRATION RATES 'RGO*P' WITH UNLIMITED
C     SUBSTRATE USED FOR MICROBIAL COMPETITION FACTOR
C
C     COQC,COQA=DOC,DOA concentration, FOCA,FOAA=DOC,DOA vs DOC+DOA
C     FCNP=N,P limitation,VMXO=specific respiration rate
C     WFNG=water stress effect, OMA=active biomass
C     TFNX=temp stress effect,FOQC,FOQA=OQC,OQA limitation
C     RGOMP=O2-unlimited respiration of DOC+DOA
C     RGOCP,RGOAP,RGOMP=O2-unlimited respiration of DOC, DOA, DOC+DOA
C
      FSBSTC=COQC(K,L,NY,NX)/(COQC(K,L,NY,NX)+OQKM)
      FSBSTA=COQA(K,L,NY,NX)/(COQA(K,L,NY,NX)+OQKA)
      FSBST=FOCA(K)*FSBSTC+FOAA(K)*FSBSTA
      RGOCY=AMAX1(0.0,FCNP(NGL,N,K)*VMXO*WFNG*OMA(NGL,N,K))
      RGOCZ=RGOCY*FSBSTC*FOCA(K)*TFNX
      RGOAZ=RGOCY*FSBSTA*FOAA(K)*TFNX
      RGOCX=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC*EO2Q)
      RGOAX=AMAX1(0.0,OQA(K,L,NY,NX)*FOQA*EO2A)
      RGOCP=AMIN1(RGOCX,RGOCZ)
      RGOAP=AMIN1(RGOAX,RGOAZ)
      RGOMP=RGOCP+RGOAP
      IF(RGOMP.GT.ZEROS(NY,NX))THEN
      FGOCP=RGOCP/RGOMP
      FGOAP=RGOAP/RGOMP
      ELSE
      FGOCP=1.0
      FGOAP=0.0
      ENDIF
C
C     ENERGY YIELD AND O2 DEMAND FROM DOC AND ACETATE OXIDATION
C     BY HETEROTROPHIC AEROBES
C
C     ECHZ=growth respiration yield
C     ROXYM,ROXYP,ROXYS=O2 demand from DOC,DOA oxidation
C     ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation
C     ROQCD=microbial respiration used to represent microbial activity
C
      ECHZ=EO2Q*FGOCP+EO2A*FGOAP
      ROXYM(NGL,N,K)=2.667*RGOMP
      ROXYP(NGL,N,K)=ROXYM(NGL,N,K)
      ROXYSX=ROXYS(NGL,N,K,L,NY,NX)
      ROQCSX=ROQCS(NGL,N,K,L,NY,NX)
      ROQASX=ROQAS(NGL,N,K,L,NY,NX)
      ROXYS(NGL,N,K,L,NY,NX)=ROXYP(NGL,N,K)
      ROQCS(NGL,N,K,L,NY,NX)=RGOCZ
      ROQAS(NGL,N,K,L,NY,NX)=RGOAZ
      ROQCD(NGL,N,K)=RGOCY
C     IF((I/1)*1.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
C     WRITE(*,5555)'RGOMP',I,J,NX,NY,L,K,N,RGOMP,RGOCX,RGOAX,RGOCZ
C    2,RGOAZ,RGOCX,RGOAX,FCNP(NGL,N,K),TFNG(NGL,N,K),VMXO,OMA(NGL,N,K),OSRH(K)
C    2,FOQC,FOQA,COQC(K,L,NY,NX),OQC(K,L,NY,NX),EO2Q,TKS(L,NY,NX)
C    3,COXYS(L,NY,NX),OQKM,OMC(1,NGL,N,K,L,NY,NX),OMC(2,NGL,N,K,L,NY,NX)
C    4,OMC(3,NGL,N,K,L,NY,NX),VOLWM(NPH,L,NY,NX),FOSRH(K,L,NY,NX)
C    5,FSBST,SPOMK(1),RMOMK(1),ROQCD(NGL,N,K),ROXYSX,ROXYS(NGL,N,K,L,NY,NX)
C    6,ROQCSX,ROQCS(NGL,N,K,L,NY,NX),ROQASX,ROQAS(NGL,N,K,L,NY,NX)
C    7,TFNX,WFNG,PSISM(L,NY,NX),TKS(L,NY,NX),TKSO
5555  FORMAT(A8,7I4,60E12.4)
C     ENDIF
C
C     RESPIRATION BY HETEROTROPHIC ANAEROBES:
C     N=(4)ACETOGENIC FERMENTERS (7) ACETOGENIC N2 FIXERS
C
C     ENERGY YIELD FROM FERMENTATION DEPENDS ON H2 AND
C     ACETATE CONCENTRATION
C
C     GH2F=energy yield of acetotrophic methanogenesis per g C
C     GHAX=H2 effect on energy yield of fermentation
C     GOAX=acetate effect on energy yield of fermentation
C     ECHZ=growth respiration efficiency of fermentation
C
      ELSEIF(N.EQ.4.OR.N.EQ.7)THEN
      GH2F=GH2X/72.0
      GOAX=8.3143E-03*TKS(L,NY,NX)
     2*LOG((AMAX1(ZERO,COQA(K,L,NY,NX))/OAKI)**2)
      GOAF=GOAX/72.0
      GHAX=GH2F+GOAF
      IF(N.EQ.4)THEN
      ECHZ=AMAX1(EO2X,AMIN1(1.0,1.0
     2/(1.0+AMAX1(0.0,(GCHX-GHAX))/EOMF)))
      ELSE
      ECHZ=AMAX1(ENFX,AMIN1(1.0,1.0
     2/(1.0+AMAX1(0.0,(GCHX-GHAX))/EOMN)))
      ENDIF
C
C     RESPIRATION RATES BY HETEROTROPHIC ANAEROBES 'RGOMP' FROM
C     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
C     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
C     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
C     MICROBIAL COMPETITION FACTOR
C
C     OXYI=O2 inhibition of fermentation
C     FCNP=N,P limitation on respiration
C     VMXF=maximum respiration rate by fermenters
C     WFNG=water stress effect on respiration
C     OMA=active fermenter biomass
C     TFNX=temp stress effect, FOQC=OQC limitation
C     RFOMP=O2-unlimited respiration of DOC
C     ROQCD=microbial respiration used to represent microbial activity
C
      FSBST=COQC(K,L,NY,NX)/(COQC(K,L,NY,NX)+OQKM)*OXYI
      RGOFY=AMAX1(0.0,FCNP(NGL,N,K)*VMXF*WFNG*OMA(NGL,N,K))
      RGOFZ=RGOFY*FSBST*TFNX
      RGOFX=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC*ECHZ)
      RGOMP=AMIN1(RGOFX,RGOFZ)
      FGOCP=1.0
      FGOAP=0.0
      ROXYM(NGL,N,K)=0.0
      ROXYP(NGL,N,K)=0.0
      ROXYS(NGL,N,K,L,NY,NX)=0.0
      ROQCS(NGL,N,K,L,NY,NX)=RGOFZ
      ROQAS(NGL,N,K,L,NY,NX)=0.0
      ROQCD(NGL,N,K)=RGOFY
      TRH2G=TRH2G+RGOMP
C     IF((I/120)*120.EQ.I.AND.J.EQ.24.AND.L.LE.6)THEN
C     WRITE(*,5554)'FERM',I,J,NX,NY,L,K,N,RGOMP,RGOFX,RGOFZ,GHAX,GOAF
C    2,ECHZ,FCNP(NGL,N,K),TFNG(NGL,N,K),OMA(NGL,N,K),OSRH(K),FOQC,COQC(K,L,NY,NX)
C    3,OQKM,OMC(1,NGL,N,K,L,NY,NX),OMC(2,NGL,N,K,L,NY,NX),OMC(3,NGL,N,K,L,NY,NX)
C    3,OMN(1,NGL,N,K,L,NY,NX),OMN(2,NGL,N,K,L,NY,NX),OMN(3,NGL,N,K,L,NY,NX)
C    5,VOLWM(NPH,L,NY,NX),PSISM(L,NY,NX),WFNG,COXYS(L,NY,NX),OXYI
C    6,FSBST,FOSRH(K,L,NY,NX),SPOMK(1),RMOMK(1),ROQCD(NGL,N,K)
5554  FORMAT(A8,7I4,60E12.4)
C     ENDIF
C
C     ENERGY YIELD FROM ACETOTROPHIC METHANOGENESIS
C
C     GOMX=acetate effect on energy yield
C     ECHZ=growth respiration efficiency of aceto. methanogenesis
C
      ELSEIF(N.EQ.5)THEN
      GOMX=8.3143E-03*TKS(L,NY,NX)
     2*LOG((AMAX1(ZERO,COQA(K,L,NY,NX))/OAKI))
      GOMM=GOMX/24.0
      ECHZ=AMAX1(EO2X,AMIN1(1.0
     2,1.0/(1.0+AMAX1(0.0,(GC4X+GOMM))/EOMH)))
C
C     RESPIRATION RATES BY ACETOTROPHIC METHANOGENS 'RGOMP' FROM
C     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
C     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL C
C     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
C     MICROBIAL COMPETITION FACTOR
C
C     COQA=DOA concentration
C     OQKAM=Km for acetate uptake,FCNP=N,P limitation
C     VMXM=specific respiration rate
C     WFNG=water stress effect, OMA=active biomass
C     TFNX=temp stress effect, FOQA= acetate limitation
C     RGOGX=substrate-limited respiration of acetate
C     RGOGX=competition-limited respiration of acetate
C     OQA=acetate, FOQA=fraction of biological demand for acetate
C     RGOMP=O2-unlimited respiration of acetate
C     ROXY*=O2 demand, ROQCS,ROQCA=DOC, acetate demand
C     ROQCD=microbial respiration used to represent microbial activity
C
      FSBST=COQA(K,L,NY,NX)/(COQA(K,L,NY,NX)+OQKAM)
      RGOGY=AMAX1(0.0,FCNP(NGL,N,K)*VMXM*WFNG*OMA(NGL,N,K))
      RGOGZ=RGOGY*FSBST*TFNX
      RGOGX=AMAX1(0.0,OQA(K,L,NY,NX)*FOQA*ECHZ)
      RGOMP=AMIN1(RGOGX,RGOGZ)
      FGOCP=0.0
      FGOAP=1.0
      ROXYM(NGL,N,K)=0.0
      ROXYP(NGL,N,K)=0.0
      ROXYS(NGL,N,K,L,NY,NX)=0.0
      ROQCS(NGL,N,K,L,NY,NX)=0.0
      ROQAS(NGL,N,K,L,NY,NX)=RGOGZ
      ROQCD(NGL,N,K)=0.0
      TCH4H=TCH4H+0.5*RGOMP
C     IF((I/30)*30.EQ.I.AND.NX.EQ.3.AND.NY.EQ.1.AND.J.EQ.24)THEN
C     WRITE(*,5552)'ACMETH',I,J,NX,NY,L,K,N,RGOMP,RGOGZ,RGOGX,GOMM
C    2,ECHZ,FCNP(NGL,N,K),TFNG(NGL,N,K),OMA(NGL,N,K),FOQA,COQA(K,L,NY,NX)
C    2,OQA(K,L,NY,NX)
C    3,OMC(1,NGL,N,K,L,NY,NX),OMC(2,NGL,N,K,L,NY,NX),OMC(3,NGL,N,K,L,NY,NX)
C    3,OMN(1,NGL,N,K,L,NY,NX),OMN(2,NGL,N,K,L,NY,NX),OMN(3,NGL,N,K,L,NY,NX)
C    5,VOLWM(NPH,L,NY,NX),PSISM(L,NY,NX),WFNG,COXYS(L,NY,NX)
C    6,OHA(K,L,NY,NX),FSBST,SPOMK(1),RMOMK(1)
5552  FORMAT(A8,7I4,40E12.4)
C     ENDIF
      ENDIF
C
C     RESPIRATION RATES BY AUTOTROPHS 'RGOMP' FROM SPECIFIC
C     OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
C     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL
C     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR
C     MICROBIAL COMPETITION FACTOR. N=(1) NH4 OXIDIZERS (2) NO2
C     OXIDIZERS,(3) CH4 OXIDIZERS, (5) H2TROPHIC METHANOGENS
C
      ELSEIF(K.EQ.5)THEN
C
C     NH3 OXIDIZERS
C
      IF(N.EQ.1)THEN
C
C     FACTOR TO REGULATE COMPETITION FOR NH4 AMONG DIFFERENT
C     MICROBIAL AND ROOT POPULATIONS FNH4
C
C     FNH4,FNB4=frac of total biol demand for NH4 in non-band, band
C
      IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4=AMAX1(FMN,RVMX4(NGL,N,K,L,NY,NX)/RNH4Y(L,NY,NX))
      ELSE
      FNH4=AMAX1(FMN,VLNH4(L,NY,NX)*FOMA(NGL,N,K))
      ENDIF
      IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB4=AMAX1(FMN,RVMB4(NGL,N,K,L,NY,NX)/RNHBY(L,NY,NX))
      ELSE
      FNB4=AMAX1(FMN,VLNHB(L,NY,NX)*FOMA(NGL,N,K))
      ENDIF
      TFNH4X=TFNH4X+FNH4
      TFNH4B=TFNH4B+FNB4
C
C     NITRIFICATION INHIBITION
C
C     ZNFN0=inhibition when fertilizer added
C     ZNFNI=reduction in inhibition since fertilizer added
C     CNH4S,CNH4B=NH4 concentrations in non-band, band
C     TFNX=temperature effect
C     RNFNI=rate constant for inhibition decline
C     ZHKI=inhibition from high CNH4
C     ZNFN4S,ZNFN4B=inhibition in non-band, band
C
      IF(ZNFN0(L,NY,NX).GT.ZEROS(NY,NX))THEN
      ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)*(1.0-RNFNI*TFNX)
      ZNFN4S=ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX)/(1.0+CNH4S(L,NY,NX)/ZHKI)
      ZNFN4B=ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX)/(1.0+CNH4B(L,NY,NX)/ZHKI)
      ELSE
      ZNFN4S=1.0
      ZNFN4B=1.0
      ENDIF
C
C     NH3 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
C     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
C     NH3 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
C
C     ECHZ=growth respiration efficiency
C     VMXX=potential NH3 oxidation, VMXH=specific oxidation
C     TFNG=temperature+water limitation, FCNP=N,P limitation
C     XCO2=aqueous CO2 limitation, OMA=active biomass
C     VMXA= non-substrate limited NH3 oxidation
C     VHKI=nonlinear increase in VMXA with VMXH
C     FNH4S,FNHBS=fractions of NH4 in non-band, band
C     CNH4S,CNH4B=NH4 concentration in non-band, band
C     ZHKM=Km for NH4 uptake
C     FNH4,FNB4=fractions of total NH4 demand in non-band, band
C     ZNH4S,ZNH4B=NH4 amount in non-band, band
C     RNNH4,RNNHB=NH3 oxidation in non-band, band
C     RGOMP=O2-unlimited respiration
C     ECNH=efficiency CO2 conversion to biomass
C     RVMX4,RVMXB=nitrifier demand for NH4 in non-band, band
C
      ECHZ=EO2X
      VMXX=VMXH*TFNG(NGL,N,K)*FCNP(NGL,N,K)*XCO2*OMA(NGL,N,K)
      IF(VOLWZ.GT.ZEROS2(NY,NX))THEN
      VMXA=VMXX/(1.0+VMXX/(VHKI*VOLWZ))
      ELSE
      VMXA=0.0
      ENDIF
      FCN4S=FNH4S*CNH4S(L,NY,NX)/(CNH4S(L,NY,NX)+ZHKM)
      FCN4B=FNHBS*CNH4B(L,NY,NX)/(CNH4B(L,NY,NX)+ZHKM)
      FSBST=FCN4S+FCN4B
      VMX4S=VMXA*FCN4S
      VMX4B=VMXA*FCN4B
      RNNH4=AMAX1(0.0,AMIN1(VMX4S,FNH4*ZNH4S(L,NY,NX)))*ZNFN4S
      RNNHB=AMAX1(0.0,AMIN1(VMX4B,FNB4*ZNH4B(L,NY,NX)))*ZNFN4B
      RVOXP=RNNH4+RNNHB
      RVOXPA=RNNH4
      RVOXPB=RNNHB
      RGOMP=AMAX1(0.0,RVOXP*ECNH*ECHZ)
      RVMX4(NGL,N,K,L,NY,NX)=VMX4S
      RVMB4(NGL,N,K,L,NY,NX)=VMX4B
C
C     O2 DEMAND FROM NH3 OXIDATION
C
C     ROXYM=O2 demand from respiration by nitrifiers
C     ROXYP,ROXYM=O2 demand from respiration + NH3 oxidation
C
      ROXYM(NGL,N,K)=2.667*RGOMP
      ROXYP(NGL,N,K)=ROXYM(NGL,N,K)+3.429*RVOXP
      ROXYS(NGL,N,K,L,NY,NX)=ROXYP(NGL,N,K)
C     IF(IYRC.EQ.2012.AND.I.EQ.151.AND.NX.EQ.1)THEN
C     WRITE(*,6666)'NITRI',I,J,L,K,N,RNNH4,RNNHB,VMXX,VMXA,VOLWZ
C    2,CNH4S(L,NY,NX),CNH4B(L,NY,NX)
C    2,14.0*XN4(L,NY,NX),14.0*XNB(L,NY,NX)
C    3,ZNH4S(L,NY,NX),ZNH4B(L,NY,NX),COXYS(L,NY,NX),RGOMP
C    4,PH(L,NY,NX),TFNX,FCNP(NGL,N,K),XCO2,ROXYM(NGL,N,K)
C    5,VMX4S,VMX4B,FCN4S,FCN4B,FNH4S,FNHBS,OMA(NGL,N,K)
C    6,FNH4,FNB4,ZNFN4S,ZNFN4B,ZNFNI(L,NY,NX),ZNFN0(L,NY,NX)
6666  FORMAT(A8,5I4,40E12.4)
C     ENDIF
C
C     NO2 OXIDIZERS
C
      ELSEIF(N.EQ.2)THEN
C
C     FACTOR TO REGULATE COMPETITION FOR NO2 AMONG DIFFERENT
C     MICROBIAL POPULATIONS
C
C     FNO2=fraction of total biological demand for NO2 in non-band, band
C
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMX2(NGL,N,K,L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=AMAX1(FMN,FOMN(NGL,N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMB2(NGL,N,K,L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=AMAX1(FMN,FOMN(NGL,N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
C
C     NO2 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
C     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
C     NO2 CONCENTRATIONS
C
C     ECHZ=growth respiration efficiency
C     VMXA= non-substrate limited NH3 oxidation
C     VMXN=specific oxidation
C     TFNG=temperature+water limitation, FCNP=N,P limitation
C     XCO2=aqueous CO2 limitation, OMA=active biomass
C     OMA=active biomass
C     FNH4S,FNHBS=fractions of NH4 in non-band, band
C     CNO2S,CNO2B=NO2 concentration in non-band, band
C     ZNKM=Km for NO2 uptake
C     FNO2,FNB2=fractions of total NO2 demand in non-band, band
C     ZNO2S,ZNO2B=NO2 amount in non-band, band
C     RNNO2,RNNOB=NO2 oxidation in non-band, band
C     RGOMP=O2-unlimited respiration
C     ECNO=efficiency CO2 conversion to biomass
C     RVMX2,RVMB2=nitrifier demand for NO2 in non-band, band
C
      ECHZ=EO2X
      VMXA=TFNG(NGL,N,K)*FCNP(NGL,N,K)*XCO2*OMA(NGL,N,K)*VMXN
      FCN2S=FNH4S*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+ZNKM)
      FCN2B=FNHBS*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+ZNKM)
      FSBST=FCN2S+FCN2B
      VMX2S=VMXA*FCN2S
      VMX2B=VMXA*FCN2B
      RNNO2=AMAX1(0.0,AMIN1(VMX2S,FNO2*ZNO2S(L,NY,NX)))
      RNNOB=AMAX1(0.0,AMIN1(VMX2B,FNB2*ZNO2B(L,NY,NX)))
      RVOXP=RNNO2+RNNOB
      RVOXPA=RNNO2
      RVOXPB=RNNOB
      RGOMP=AMAX1(0.0,RVOXP*ECNO*ECHZ)
      RVMX2(NGL,N,K,L,NY,NX)=VMX2S
      RVMB2(NGL,N,K,L,NY,NX)=VMX2B
C
C     O2 DEMAND FROM NO2 OXIDATION
C
C     ROXYM=O2 demand from respiration by nitrifiers
C     ROXYP,ROXYM=O2 demand from respiration + NO2 oxidation
C
      ROXYM(NGL,N,K)=2.667*RGOMP
      ROXYP(NGL,N,K)=ROXYM(NGL,N,K)+1.143*RVOXP
      ROXYS(NGL,N,K,L,NY,NX)=ROXYP(NGL,N,K)
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.LE.6)THEN
C     WRITE(*,6667)'NO2OX',I,J,L,K,N,RNNO2,RNNOB,ZNO2S(L,NY,NX)
C    2,ZNO2B(L,NY,NX),CNO2S(L,NY,NX),CNO2B(L,NY,NX),CNH3S(L,NY,NX)
C    3,CNH3B(L,NY,NX),CNH4S(L,NY,NX),CNH4B(L,NY,NX),CNO3S(L,NY,NX)
C    3,CNO3B(L,NY,NX),CHNO2,CHNOB,VMXA,TFNG(NGL,N,K),FCNP(NGL,N,K),VMXN,ZNKM
C    4,FCN2S,FCN2B,OMA(NGL,N,K),FOMN(NGL,N,K),TOMN,RVMX2(NGL,N,K,L,NY,NX)
C    5,RNO2Y(L,NY,NX),FNO2,FNB2,ROXYM(NGL,N,K),ROXYP(NGL,N,K)
C    6,ROXYS(NGL,N,K,L,NY,NX),VLNHB(L,NY,NX),VLNOB(L,NY,NX)
C    7,SPOMK(1),RMOMK(1)
6667  FORMAT(A8,5I4,50E12.4)
C     ENDIF
C
C     H2TROPHIC METHANOGENS
C
      ELSEIF(N.EQ.5)THEN
C
C     CO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND H2
C
C     GH2H=energy yield of hydrogenotrophic methanogenesis per g C
C     ECHZ=growth respiration efficiency of hydrogen. methanogenesis
C     VMXA=substrate-unlimited H2 oxidation rate
C     H2GSX=aqueous H2 (H2GS) + total H2 from fermentation (TRH2G)
C     CH2GS=H2 concentration, H2KM=Km for H2 uptake
C     RGOMP=H2 oxidation, ROXY*=O2 demand
C
      GH2H=GH2X/12.0
      ECHZ=AMAX1(EO2X,AMIN1(1.0
     2,1.0/(1.0+AMAX1(0.0,(GCOX+GH2H))/EOMH)))
      VMXA=TFNG(NGL,N,K)*FCNP(NGL,N,K)*XCO2*OMA(NGL,N,K)*VMXC
      H2GSX=H2GS(L,NY,NX)+0.111*TRH2G
      FSBST=CH2GS(L,NY,NX)/(CH2GS(L,NY,NX)+H2KM)
      RGOMP=AMAX1(0.0,AMIN1(1.5*H2GSX,VMXA*FSBST))
      ROXYM(NGL,N,K)=0.0
      ROXYP(NGL,N,K)=0.0
      ROXYS(NGL,N,K,L,NY,NX)=0.0
      TCH4A=TCH4A+RGOMP
C     IF((I/30)*30.EQ.I.AND.NX.EQ.3.AND.NY.EQ.1.AND.J.EQ.24)THEN
C     WRITE(*,5553)'H2METH',I,J,NX,NY,L,K,N,RGOMP,H2GS(L,NY,NX)
C    2,H2GSX,CH2GS(L,NY,NX),VMXA,TFNG(NGL,N,K),FCNP(NGL,N,K),XCO2
C    3,OMA(NGL,N,K),VMXC,ECHZ,GCOX,GH2H,TKS(L,NY,NX),FSBST
C    4,SPOMK(1),RMOMK(1)
5553  FORMAT(A8,7I4,20E12.4)
C     ENDIF
C
C     METHANOTROPHS
C
      ELSEIF(N.EQ.3)THEN
C
C     CH4 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
C     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
C     CH4 CONCENTRATIONS IN BAND AND NON-BAND SOIL ZONES
C
C     ECHZ=growth respiration efficiency
C     VMXA=potential oxidation
C     TFNG=temperature+water effect,FCNP=N,P limitation
C     OMA=active biomass,VMX4=specific respiration rate
C     RCH4L=total aqueous CH4 exchange from previous hour
C     RCH4F=total gaseous CH4 exchange from previous hour
C     TCH4H+TCH4A=total CH4 generated from methanogenesis
C     XNPG=1.0/(NPH*NPT)
C     CH4G1,CH4S1=CH4 gaseous, aqueous amounts
C     CCH4E,CCH4G=CH4 gas concentration in atmosphere, soil
C     VOLPM,VOLWM=air,water-filled porosity
C     SCH4L=CH4 aqueous solubility
C     CCK4=Km for CH4 uptake
C     ECHO=efficiency CO2 conversion to biomass
C     RGOMP1=substrate-limited CH4 oxidation
C     RCHDF=gaseous-aqueous CH4 exchange
C     DFGS=rate constant for gaseous-aqueous exchange
C
      ECHZ=EH4X
      VMXA=TFNG(NGL,N,K)*FCNP(NGL,N,K)*OMA(NGL,N,K)*VMX4
      RCH4L1=RCH4L(L,NY,NX)*XNPG
      RCH4F1=RCH4F(L,NY,NX)*XNPG
      RCH4S1=(TCH4H+TCH4A)*XNPG
      IF(L.EQ.0)THEN
      CH4G1=CCH4E(NY,NX)*VOLPM(1,L,NY,NX)
      ELSE
      CH4G1=CCH4G(L,NY,NX)*VOLPM(1,L,NY,NX)
      ENDIF
      CH4S1=CH4S(L,NY,NX)
      VMXA1=VMXA*XNPG
      RVOXP=0.0
      RGOMP=0.0
C
C     CH4 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
C     TO MAINTAIN AQUEOUS CH4 CONCENTRATION DURING OXIDATION
C
      DO 320 M=1,NPH
      IF(VOLWM(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      VOLWCH=VOLWM(M,L,NY,NX)*SCH4L(L,NY,NX)
      VOLWPM=VOLWCH+VOLPM(M,L,NY,NX)
      DO 325 MM=1,NPT
      CH4G1=CH4G1+RCH4F1
      CH4S1=CH4S1+RCH4L1+RCH4S1
      CCH4S1=AMAX1(0.0,CH4S1/VOLWM(M,L,NY,NX))
      FSBST=CCH4S1/(CCH4S1+CCK4)
      RVOXP1=AMIN1(AMAX1(0.0,CH4S1)/(1.0+ECHO*ECHZ)
     2,VMXA1*FSBST)
      RGOMP1=RVOXP1*ECHO*ECHZ
      CH4S1=CH4S1-RVOXP1-RGOMP1
      IF(THETPM(M,L,NY,NX).GT.THETX)THEN
      RCHDF=DFGS(M,L,NY,NX)*(AMAX1(ZEROS(NY,NX),CH4G1)*VOLWCH
     2-CH4S1*VOLPM(M,L,NY,NX))/VOLWPM
      ELSE
      RCHDF=0.0
      ENDIF
      CH4G1=CH4G1-RCHDF
      CH4S1=CH4S1+RCHDF
      RVOXP=RVOXP+RVOXP1
      RGOMP=RGOMP+RGOMP1
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.0
C    2.AND.MM.EQ.NPT)THEN
C     WRITE(*,5547)'CH4OX',I,J,NX,NY,L,K,N,M,MM,RVOXP1,RGOMP1,CH4G1
C     2,CH4S1,VMXA1,RVOXP,RGOMP,RCHDF,RCH4L1,RCH4F1,RCH4S1,CCH4S1
C    3,ECHO,ECHZ,OMA(NGL,N,K),VOLWM(M,L,NY,NX),VOLPM(M,L,NY,NX),VOLWCH
C    4,THETPM(M,L,NY,NX),SCH4L(L,NY,NX),DFGS(M,L,NY,NX)
C    5,COXYS(L,NY,NX),CCH4E(NY,NX),FSBST,SPOMK(1),RMOMK(1)
C    6,CH4G1/VOLPM(M,L,NY,NX)
5547  FORMAT(A8,9I4,30E12.4)
C     ENDIF
325   CONTINUE
      ENDIF
320   CONTINUE
      RVOXPA=RVOXP
      RVOXPB=0.0
C
C     O2 DEMAND FROM CH4 OXIDATION
C
C     ROXYM=O2 demand from respiration
C     ROXYP=O2 demand from respiration + CH4 oxidation
C
      ROXYM(NGL,N,K)=2.667*RGOMP
      ROXYP(NGL,N,K)=ROXYM(NGL,N,K)+4.00*RVOXP
      ROXYS(NGL,N,K,L,NY,NX)=ROXYP(NGL,N,K)
      ELSE
      RGOMP=0.0
      ROXYM(NGL,N,K)=0.0
      ROXYP(NGL,N,K)=0.0
      ROXYS(NGL,N,K,L,NY,NX)=0.0
      ENDIF
      ELSE
      RGOMP=0.0
      ROXYM(NGL,N,K)=0.0
      ROXYP(NGL,N,K)=0.0
      ROXYS(NGL,N,K,L,NY,NX)=0.0
      ENDIF
C
C     O2 UPTAKE BY AEROBES
C
C     RUPOX, ROXYP=O2-limited, O2-unlimited rates of O2 uptake
C     RUPMX=O2-unlimited rate of O2 uptake
C     FOXYX=fraction of O2 uptake by N,K relative to total
C     XNPG=1/(NPH*NPT)
C     ROXYF,ROXYL=net O2 gaseous, aqueous fluxes from previous hour
C     OLSGL=aqueous O2 diffusivity
C     OXYG,OXYS=gaseous, aqueous O2 amounts
C     FLQRQ,FLQRI=surface water flux from precipitation, irrigation
C     COXR,COXQ=O2 concentration in FLQRQ,FLQRI
C
      RUPOX(NGL,N,K)=0.0
      IF(N.LE.3.OR.N.EQ.6)THEN
      IF(ROXYP(NGL,N,K).GT.ZEROS(NY,NX).AND.FOXYX.GT.ZERO)THEN
      IF(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX))THEN
C
C     MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC
C     POPULATION
C
      RUPMX=ROXYP(NGL,N,K)*XNPG
      ROXYFX=ROXYF(L,NY,NX)*XNPG*FOXYX
      OLSGL1=OLSGL(L,NY,NX)*XNPG
      IF(L.NE.0)THEN
      OXYG1=OXYG(L,NY,NX)*FOXYX
      ROXYLX=ROXYL(L,NY,NX)*XNPG*FOXYX
      ELSE
      OXYG1=COXYG(L,NY,NX)*VOLPM(1,L,NY,NX)*FOXYX
      ROXYLX=(ROXYL(L,NY,NX)+FLQRQ(NY,NX)*COXR(NY,NX)
     2+FLQRI(NY,NX)*COXQ(NY,NX))*XNPG*FOXYX
      ENDIF
      OXYS1=OXYS(L,NY,NX)*FOXYX
C
C     O2 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
C     TO MAINTAIN AQUEOUS O2 CONCENTRATION DURING REDUCTION
C
      DO 420 M=1,NPH
C
C     ACTUAL REDUCTION OF AQUEOUS BY AEROBES CALCULATED
C     FROM MASS FLOW PLUS DIFFUSION = ACTIVE UPTAKE
C     COUPLED WITH DISSOLUTION OF GASEOUS O2 DURING REDUCTION
C     OF AQUEOUS O2 FROM DISSOLUTION RATE CONSTANT 'DFGS'
C     CALCULATED IN 'WATSUB'
C
C     VOLWM,VOLPM,VOLX=water, air and total volumes
C     ORAD=microbial radius,FILM=water film thickness
C     DIFOX=aqueous O2 diffusion, TORT=tortuosity
C     BIOS=microbial number, OMA=active biomass
C     SOXYL=O2 solubility, OXKX=Km for O2 uptake
C     OXYS,COXYS=aqueous O2 amount, concentration
C     OXYG,COXYG=gaseous O2 amount, concentration
C     RMPOX,ROXSK=O2 uptake
C
      THETW1=AMAX1(0.0,VOLWM(M,L,NY,NX)/VOLY(L,NY,NX))
      RRADO=ORAD*(FILM(M,L,NY,NX)+ORAD)/FILM(M,L,NY,NX)
      DIFOX=TORT(M,L,NY,NX)*OLSGL1*12.57*BIOS*OMA(NGL,N,K)*RRADO
      VOLWOX=VOLWM(M,L,NY,NX)*SOXYL(L,NY,NX)
      VOLPOX=VOLPM(M,L,NY,NX)
      VOLWPM=VOLWOX+VOLPOX
      DO 425 MX=1,NPT
      OXYG1=OXYG1+ROXYFX
      OXYS1=OXYS1+ROXYLX
      COXYS1=AMIN1(COXYE(NY,NX)*SOXYL(L,NY,NX)
     2,AMAX1(0.0,OXYS1/(VOLWM(M,L,NY,NX)*FOXYX)))
      X=DIFOX*COXYS1
      IF(X.GT.ZEROS(NY,NX).AND.OXYS1.GT.ZEROS(NY,NX))THEN
      B=-RUPMX-DIFOX*OXKX-X
      C=X*RUPMX
      RMPOX=(-B-SQRT(B*B-4.0*C))/2.0
      ELSE
      RMPOX=0.0
      ENDIF
      OXYS1=OXYS1-RMPOX
      IF(THETPM(M,L,NY,NX).GT.THETX.AND.VOLPOX.GT.ZEROS(NY,NX))THEN
      ROXDFQ=DFGS(M,L,NY,NX)*(AMAX1(ZEROS(NY,NX),OXYG1)*VOLWOX
     2-OXYS1*VOLPOX)/VOLWPM
      ELSE
      ROXDFQ=0.0
      ENDIF
      OXYG1=OXYG1-ROXDFQ
      OXYS1=OXYS1+ROXDFQ
      RUPOX(NGL,N,K)=RUPOX(NGL,N,K)+RMPOX
      ROXSK(M,L,NY,NX)=ROXSK(M,L,NY,NX)+RMPOX
C     IF(I.EQ.151.AND.J.EQ.24.AND.L.LE.5.AND.M.EQ.NPH.AND.MX.EQ.NPT)THEN
C     WRITE(*,5545)'RMPOX',I,J,L,K,N,M,MX,OXYS1,ROXDFQ,ROXYLX,RMPOX
C    2,DFGS(M,L,NY,NX),OXYG1,VOLWOX,VOLPOX,VOLWPM,X,B,C
C    3,RUPMX,DIFOX,OXKX,COXYS1,FOXYX,ROXYL(L,NY,NX)
C    4,ROXSK(M,L,NY,NX),VOLWM(M,L,NY,NX)/VOLY(L,NY,NX)
C    5,OXYS(L,NY,NX)
5545  FORMAT(A8,7I4,30E16.6)
C     ENDIF
C     IF((I/120)*120.EQ.I.AND.J.EQ.24.AND.L.LE.3
C    2.AND.K.GE.3.AND.N.EQ.3)THEN
C     WRITE(*,5544)'OXY',I,J,L,K,N,M,MX,RUPOX(NGL,N,K),ROXYP(NGL,N,K)
C    2,ROXSK(M,L,NY,NX),RUPMX,RMPOX,DIFOX,OLSGL1,BIOS,OMA(NGL,N,K),X
C    2,ROXDFQ,ROXYLX,ROXYFX,FOXYX,COXYS1,OXYS1,OXYG1,OXYS1
C    4/(VOLWM(M,L,NY,NX)*FOXYX),OXYG1/(VOLPM(M,L,NY,NX)*FOXYX)
C    5,THETW1,THETPM(M,L,NY,NX),DFGS(M,L,NY,NX),ROXSK(M,L,NY,NX)
C    6,VOLPM(M,L,NY,NX),VOLWM(M,L,NY,NX),VOLA(L,NY,NX)
C    7,COXYS(L,NY,NX),COXYG(L,NY,NX),ROXYY(L,NY,NX)
5544  FORMAT(A8,7I4,50E12.4)
C     ENDIF
425   CONTINUE
420   CONTINUE

C
C     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (WFN)
C
C     WFN=ratio of O2-limited to O2-unlimited uptake
C     RVMX4,RVNHB,RVMX2,RVMB2=NH3,NO2 oxidation in non-band, band
C
      WFN(NGL,N,K)=AMIN1(1.0,AMAX1(0.0,RUPOX(NGL,N,K)/ROXYP(NGL,N,K)))
C     IF(K.LE.4)THEN
C     ROQCS(NGL,N,K,L,NY,NX)=ROQCS(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
C     ROQAS(NGL,N,K,L,NY,NX)=ROQAS(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
C     ROQCD(NGL,N,K)=ROQCD(NGL,N,K)*WFN(NGL,N,K)
C     ENDIF
      IF(K.EQ.5)THEN
      IF(N.EQ.1)THEN
      RVMX4(NGL,N,K,L,NY,NX)=RVMX4(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
      RVMB4(NGL,N,K,L,NY,NX)=RVMB4(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
      ELSEIF(N.EQ.2)THEN
      RVMX2(NGL,N,K,L,NY,NX)=RVMX2(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
      RVMB2(NGL,N,K,L,NY,NX)=RVMB2(NGL,N,K,L,NY,NX)*WFN(NGL,N,K)
      ENDIF
      ENDIF
      ELSE
      RUPOX(NGL,N,K)=ROXYP(NGL,N,K)
      WFN(NGL,N,K)=1.0
      ENDIF
      ELSE
      RUPOX(NGL,N,K)=0.0
      WFN(NGL,N,K)=1.0
      ENDIF

C
C     RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2
C
C     RGOMO,RGOMP=O2-limited, O2-unlimited respiration
C     RCO2X,RCH3X,RCH4X,RH2GX=CO2,acetate,CH4,H2 production from RGOMO
C     ROXYO=O2-limited O2 uptake
C     RVOXA,RVOXB=total O2-lmited (1)NH4,(2)NO2,(3)CH4 oxidation
C
      RGOMO(NGL,N,K)=RGOMP*WFN(NGL,N,K)
      RCO2X(NGL,N,K)=RGOMO(NGL,N,K)
      RCH3X(NGL,N,K)=0.0
      RCH4X(NGL,N,K)=0.0
      ROXYO(NGL,N,K)=ROXYM(NGL,N,K)*WFN(NGL,N,K)
      RH2GX(NGL,N,K)=0.0
      IF(K.EQ.5)THEN
      RVOXA(NGL,N)=RVOXPA*WFN(NGL,N,K)
      RVOXB(NGL,N)=RVOXPB*WFN(NGL,N,K)
      ENDIF
      ELSEIF(N.EQ.4.OR.N.EQ.7)THEN
      RGOMO(NGL,N,K)=RGOMP
      RCO2X(NGL,N,K)=0.333*RGOMO(NGL,N,K)
      RCH3X(NGL,N,K)=0.667*RGOMO(NGL,N,K)
      RCH4X(NGL,N,K)=0.0
      ROXYO(NGL,N,K)=ROXYM(NGL,N,K)
      IF(K.LE.4)THEN
      RH2GX(NGL,N,K)=0.111*RGOMO(NGL,N,K)
      ELSE
      RH2GX(NGL,N,K)=0.0
      ENDIF
      ELSEIF(N.EQ.5)THEN
      RGOMO(NGL,N,K)=RGOMP
      IF(K.LE.4)THEN
      RCO2X(NGL,N,K)=0.50*RGOMO(NGL,N,K)
      RCH3X(NGL,N,K)=0.00
      RCH4X(NGL,N,K)=0.50*RGOMO(NGL,N,K)
      ROXYO(NGL,N,K)=ROXYM(NGL,N,K)
      RH2GX(NGL,N,K)=0.0
      ELSEIF(K.EQ.5)THEN
      RCO2X(NGL,N,K)=0.00
      RCH3X(NGL,N,K)=0.00
      RCH4X(NGL,N,K)=RGOMO(NGL,N,K)
      ROXYO(NGL,N,K)=ROXYM(NGL,N,K)
      RH2GX(NGL,N,K)=0.0
      RH2GZ=0.667*RGOMO(NGL,N,K)
      ENDIF
      ENDIF
C
C     HETEROTROPHIC DENITRIFICATION
C
      IF(K.LE.4.AND.N.EQ.2.AND.ROXYM(NGL,N,K).GT.0.0
     2.AND.(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX)))THEN
C
C     FACTOR TO CONSTRAIN NO3 UPAKE AMONG COMPETING MICROBIAL
C     AND ROOT POPULATIONS
C
C     FNO3,FNB3=fraction of total biological demand for NO3
C
      IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3=AMAX1(FMN,RVMX3(NGL,N,K,L,NY,NX)/RNO3Y(L,NY,NX))
      ELSE
      FNO3=AMAX1(FMN,FOMA(NGL,N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB3=AMAX1(FMN,RVMB3(NGL,N,K,L,NY,NX)/RN3BY(L,NY,NX))
      ELSE
      FNB3=AMAX1(FMN,FOMA(NGL,N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO3X=TFNO3X+FNO3
      TFNO3B=TFNO3B+FNB3
C
C     NO3 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO3
C     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
C     NOT ACCEPTED BY O2 IN BAND AND NON-BAND SOIL ZONES
C
C     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO
C     VMXD3=demand for NO3-N reduction
C     VMXDXS,VMXDXB=maximum NO3 reduction in non-band, band
C     FNO3S,FNO3B=fractions of total NO3 in non-band, band
C     CNO3S,CNO3B=NO3 concentrations in non-band, band
C     Z3KM,Z2KM=Km for NO3, NO2 uptake
C     FVMXDX=nonlinear effect of product inhibition for NOx reduction
C     VMKI=product inhibition for NOx reduction
C     VMXD3S,VMXD3B=substrate-unlimited NO3 reduction in non-band,band
C     OQCD3S,OQCD3B=DOC limitation to NO3 reduction in non-band, band
C     RDNO3,RDNOB=substrate-limited NO3 reduction in non-band,band
C     RGOM3X,RGOMD3=substrate-unltd,-ltd respn from NO3 reduction
C     RVMX3,RVMB3=demand for NO3 reduction in non-band,band
C
      ROXYD=AMAX1(0.0,ROXYM(NGL,N,K)-ROXYO(NGL,N,K))
      VMXD3=0.875*ROXYD
      IF(CNO3S(L,NY,NX).GT.ZERO)THEN
      VMXDXS=FNO3S*VMXD3*CNO3S(L,NY,NX)/(CNO3S(L,NY,NX)+Z3KM)
     2/(1.0+(CNO2S(L,NY,NX)*Z3KM)/(CNO3S(L,NY,NX)*Z2KM))
      ELSE
      VMXDXS=0.0
      ENDIF
      IF(CNO3B(L,NY,NX).GT.ZERO)THEN
      VMXDXB=FNO3B*VMXD3*CNO3B(L,NY,NX)/(CNO3B(L,NY,NX)+Z3KM)
     2/(1.0+(CNO2B(L,NY,NX)*Z3KM)/(CNO3B(L,NY,NX)*Z2KM))
      ELSE
      VMXDXB=0.0
      ENDIF
      VMXDXT=VMXDXS+VMXDXB
      IF(VOLWZ.GT.ZEROS2(NY,NX).AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
      FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ*FOSRH(K,L,NY,NX)))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD3S=VMXDXS*FVMXDX
      VMXD3B=VMXDXB*FVMXDX
      OQCZ3=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC-RGOCP*WFN(NGL,N,K))
      OQCD3=OQCZ3/ECN3
      OQCD3S=OQCD3*FNO3S
      OQCD3B=OQCD3*FNO3B
      ZNO3SX=ZNO3S(L,NY,NX)*FNO3
      ZNO3BX=ZNO3B(L,NY,NX)*FNB3
      RDNO3X=AMAX1(0.0,AMIN1(ZNO3SX,VMXD3S))
      RDNOBX=AMAX1(0.0,AMIN1(ZNO3BX,VMXD3B))
      RDNO3(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD3S,OQCD3S,ZNO3SX))
      RDNOB(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD3B,OQCD3B,ZNO3BX))
      RDNOX=RDNO3X+RDNOBX
      RDNOT=RDNO3(NGL,N,K)+RDNOB(NGL,N,K)
      RGOM3X=ECN3*RDNOX
      RGOMD3=ECN3*RDNOT
      RVMX3(NGL,N,K,L,NY,NX)=VMXD3S
      RVMB3(NGL,N,K,L,NY,NX)=VMXD3B
C
C     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
C     POPULATIONS
C
C     FNO2,FNB2=fraction of total biological demand for NO2
C
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMX2(NGL,N,K,L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=AMAX1(FMN,FOMA(NGL,N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMB2(NGL,N,K,L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=AMAX1(FMN,FOMA(NGL,N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
C
C     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2
C     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
C     NOT ACCEPTED BY O2 AND NO3 IN BAND AND NON-BAND SOIL ZONES
C
C     VMXD2=demand for NO2-N reduction
C     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
C     FNO2S,FNO2B=fractions of total NO2 in non-band, band
C     CNO2S,CNO2B=NO2 concentrations in non-band, band
C     Z2KM,Z1KM=Km for NO2, N2O uptake
C     FVMXDX=nonlinear effect of product inhibition for NOx reduction
C     VMKI=product inhibition for NOx reduction
C     VMXD2S,VMXD2B=substrate-unlimited NO2 reduction in non-band,band
C     OQCD2S,OQCD2B=DOC limitation to NO2 reduction in non-band, band
C     RDNO2,RDN2B=substrate-limited NO2 reduction in non-band,band
C     RGOM2X,RGOMD2=substrate-unltd,-ltd respn from NO2 reduction
C
      VMXD2=VMXD3-RDNOT
      IF(CNO2S(L,NY,NX).GT.ZERO)THEN
      VMXDXS=FNO2S*VMXD2*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+Z2KM)
     2/(1.0+(CZ2OS(L,NY,NX)*Z2KM)/(CNO2S(L,NY,NX)*Z1KM))
      ELSE
      VMXDXS=0.0
      ENDIF
      IF(CNO2B(L,NY,NX).GT.ZERO)THEN
      VMXDXB=FNO2B*VMXD2*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+Z2KM)
     2/(1.0+(CZ2OS(L,NY,NX)*Z2KM)/(CNO2B(L,NY,NX)*Z1KM))
      ELSE
      VMXDXB=0.0
      ENDIF
      VMXDXT=VMXDXS+VMXDXB
      IF(VOLWZ.GT.ZEROS2(NY,NX).AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
      FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ*FOSRH(K,L,NY,NX)))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD2S=VMXDXS*FVMXDX
      VMXD2B=VMXDXB*FVMXDX
      OQCZ2=AMAX1(0.0,OQCZ3-RGOMD3)
      OQCD2=OQCZ2/ECN2
      OQCD2S=OQCD2*FNO3S
      OQCD2B=OQCD2*FNO3B
      ZNO2SX=(ZNO2S(L,NY,NX)+RDNO3(NGL,N,K))*FNO2
      ZNO2BX=(ZNO2B(L,NY,NX)+RDNOB(NGL,N,K))*FNB2
      RDNO2X=AMAX1(0.0,AMIN1(ZNO2SX,VMXD2S))
      RDNOBX=AMAX1(0.0,AMIN1(ZNO2BX,VMXD2B))
      RDNO2(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD2S,OQCD2S,ZNO2SX))
      RDN2B(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD2B,OQCD2B,ZNO2BX))
      RDN2X=RDNO2X+RDNOBX
      RDN2T=RDNO2(NGL,N,K)+RDN2B(NGL,N,K)
      RGOM2X=ECN2*RDN2X
      RGOMD2=ECN2*RDN2T
      RVMX2(NGL,N,K,L,NY,NX)=VMXD2S
      RVMB2(NGL,N,K,L,NY,NX)=VMXD2B
C
C     FACTOR TO CONSTRAIN N2O UPAKE AMONG COMPETING MICROBIAL
C     AND ROOT POPULATIONS
C
C     FN2O=fraction of total biological demand for N2O
C
      IF(RN2OY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FN2O=AMAX1(FMN,RVMX1(NGL,N,K,L,NY,NX)/RN2OY(L,NY,NX))
      ELSE
      FN2O=AMAX1(FMN,FOMA(NGL,N,K))
      ENDIF
      TFN2OX=TFN2OX+FN2O
C
C     N2O REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS N2O
C     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
C     NOT ACCEPTED BY O2, NO3 AND NO2 IN BAND AND NON-BAND SOIL ZONES
C
C     VMXD1=demand for N2O-N reduction
C     VMXDXS=maximum N2O reduction
C     CZ2OS=N2O concentrations
C     Z1KM=Km for N2O uptake
C     FVMXDX=nonlinear effect of product inhibition for NOx reduction
C     VMKI=product inhibition for NOx reduction
C     VMXD1S=substrate-unlimited N2O reduction
C     OQCD1=DOC limitation to N2O reduction
C     RDN2O=substrate-limited N2O reduction
C     RGOM1X,RGOMD1=substrate-unltd,-ltd  respn from N2O reduction
C     RGOMY,RGOMD=total substrate-unltd,-ltd respn from NOx reduction
C     RVMX1=demand for N2O reduction
C
      VMXD1=(VMXD2-RDN2T)*2.0
      VMXDXS=VMXD1*CZ2OS(L,NY,NX)/(CZ2OS(L,NY,NX)+Z1KM)
      IF(VOLWZ.GT.ZEROS2(NY,NX).AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
      FVMXDX=1.0/(1.0+VMXDXS/(VMKI*VOLWZ*FOSRH(K,L,NY,NX)))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD1S=VMXDXS*FVMXDX
      OQCZ1=AMAX1(0.0,OQCZ2-RGOMD2)
      OQCD1=OQCZ1/ECN1
      Z2OSX=(Z2OS(L,NY,NX)+RDN2T)*FN2O
      RDN2OX=AMAX1(0.0,AMIN1(Z2OSX,VMXD1S))
      RDN2O(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD1S,OQCD1,Z2OSX))
      RGOM1X=ECN1*RDN2OX
      RGOMD1=ECN1*RDN2O(NGL,N,K)
      RGOMY(NGL,N,K)=RGOM3X+RGOM2X+RGOM1X
      RGOMD(NGL,N,K)=RGOMD3+RGOMD2+RGOMD1
      RVMX1(NGL,N,K,L,NY,NX)=VMXD1S
C     TRN2OD(NY,NX)=TRN2OD(NY,NX)+RDNO2(NGL,N,K)+RDN2B(NGL,N,K)
C     TRN2GD(NY,NX)=TRN2GD(NY,NX)+RDN2O(NGL,N,K)
C     IF((I/1)*1.EQ.I.AND.L.LE.5)THEN
C     WRITE(*,2222)'DENIT',I,J,L,K,N,RDNO3(NGL,N,K),RDNOB(NGL,N,K),RDNO2(NGL,N,K)
C    2,RDN2B(NGL,N,K),RDN2O(NGL,N,K),TRN2OD(NY,NX),TRN2GD(NY,NX)
C    3,COXYS(L,NY,NX),COXYG(L,NY,NX),ROXYM(NGL,N,K)
C    3,ROXYO(NGL,N,K),OMA(NGL,N,K),VMXD,CNO3S(L,NY,NX),CNO3B(L,NY,NX)
C    4,CNO2S(L,NY,NX),CNO2B(L,NY,NX),CZ2OS(L,NY,NX),VLNO3(L,NY,NX)
C    5,VLNOB(L,NY,NX),THETW(L,NY,NX),THETI(L,NY,NX),FOMA(NGL,N,K)
C    5,ZNO3S(L,NY,NX),ZNO3B(L,NY,NX),ZNO2S(L,NY,NX),ZNO2B(L,NY,NX)
C    6,Z2OS(L,NY,NX),RGOMY(NGL,N,K),RGOMD(NGL,N,K),TOMA,FOXYX,FNO23S,FNO23B
C    7,OQC(K,L,NY,NX),FOQC,RGOCP,WFN(NGL,N,K),VOLWZ,FOSRH(K,L,NY,NX),ZERO
C    9,RGOM3X,RGOM2X,RGOM1X,FNO3,FNO2,FN2O,ZNO3SX,ZNO2SX,Z2OSX
C    3,OQCD3S,OQCD2S,OQCD1,VMXD3S,VMXD2S,VMXD1S,VMXD3,VMXD2,VMXD1
C    4,ROXYD,VMXDX,TFNX,WFNG,TFNG(NGL,N,K),PSISM(L,NY,NX)
C    2,(1.0+(CNO2S(L,NY,NX)*Z3KM)/(CNO3S(L,NY,NX)*Z2KM))
C    2,(1.0+(CZ2OS(L,NY,NX)*Z2KM)/(CNO2S(L,NY,NX)*Z1KM))
2222  FORMAT(A8,5I4,70E12.4)
C     ENDIF
C
C     AUTOTROPHIC DENITRIFICATION
C
      ELSEIF(K.EQ.5.AND.N.EQ.1.AND.ROXYM(NGL,N,K).GT.0.0
     2.AND.(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX)))THEN
C
C     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
C     POPULATIONS
C
C     FNO2,FNB2=fraction of total biological demand for NO2
C
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMX2(NGL,N,K,L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=AMAX1(FMN,FOMN(NGL,N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMB2(NGL,N,K,L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=AMAX1(FMN,FOMN(NGL,N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
C
C     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE NITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2 AND CO2
C     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
C     NOT ACCEPTED BY O2
C
C     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO
C     VMXD4=demand for NO2-N reduction
C     VMXDXS,VMXDXB=maximum NO2 reduction in non-band, band
C     FNO2S,FNO2B=fractions of total NO2 in non-band, band
C     CNO2S,CNO2B=NO2 concentrations in non-band, band
C     Z2KM=Km for NO2 uptake
C     FVMXDX=nonlinear effect of product inhibition for NOx reduction
C     VMKI=product inhibition for NOx reduction
C     VMXD4S,VMXD4B=substrate-unlimited NO2 reduction in non-band,band
C     RDNO2,RDN2B=substrate-limited NO2 reduction in non-band,band
C     RGOMY,RGOMD=total substrate-unltd,-ltd respn from NO2 reduction
C     ECNO=efficiency CO2 conversion to biomass
C     ECHZ=growth respiration efficiency
C     RVOXA,RVOXB=total O2-limited (1)NH4,(2)NO2,(3)CH4 oxidation
C
      ROXYD=AMAX1(0.0,ROXYM(NGL,N,K)-ROXYO(NGL,N,K))
      VMXD4=0.875*ROXYD*XCO2
      VMXDXS=FNO2S*VMXD4*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+Z2KM)
      VMXDXB=FNO2B*VMXD4*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+Z2KM)
      VMXDXT=VMXDXS+VMXDXB
      IF(VOLWZ.GT.ZEROS2(NY,NX))THEN
      FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWZ))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD4S=VMXDXS*FVMXDX
      VMXD4B=VMXDXB*FVMXDX
      ZNO2SX=ZNO2S(L,NY,NX)+RVOXA(NGL,1)
      ZNO2BX=ZNO2B(L,NY,NX)+RVOXB(NGL,1)
      RDNO2(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD4S,ZNO2SX))
      RDN2B(NGL,N,K)=AMAX1(0.0,AMIN1(VMXD4B,ZNO2BX))
      RDNOT=RDNO2(NGL,N,K)+RDN2B(NGL,N,K)
      RGOMY(NGL,N,K)=0.0
      RGOMD(NGL,N,K)=RDNOT*ECNO*ENOX
      RDNO3(NGL,N,K)=0.0
      RDNOB(NGL,N,K)=0.0
      RDN2O(NGL,N,K)=0.0
      RVMX2(NGL,N,K,L,NY,NX)=VMXD4S
      RVMB2(NGL,N,K,L,NY,NX)=VMXD4B
      RVOXA(NGL,N)=RVOXA(NGL,N)+0.333*RDNO2(NGL,N,K)
      RVOXB(NGL,N)=RVOXB(NGL,N)+0.333*RDN2B(NGL,N,K)
C     TRN2ON(NY,NX)=TRN2ON(NY,NX)+RDNO2(NGL,N,K)+RDN2B(NGL,N,K)
C     IF((I/1)*1.EQ.I.AND.J.EQ.19.AND.L.LE.5)THEN
C     WRITE(*,7777)'AUTO',I,J,L,K,N,RDNO2(NGL,N,K)
C    2,RDN2B(NGL,N,K),TRN2ON(NY,NX)
C    2,CNO2S(L,NY,NX),CNO2B(L,NY,NX),CNO3S(L,NY,NX),CNO3B(L,NY,NX)
C    3,Z2OS(L,NY,NX),VLNOB(L,NY,NX),ZNO2S(L,NY,NX),ZNO2B(L,NY,NX)
C    3,XCO2,FNO2,FNB2,TFNG(NGL,N,K),OMA(NGL,N,K),ROXYP(NGL,N,K)
C    2,ROXYM(NGL,N,K),ROXYO(NGL,N,K),WFN(NGL,N,K),FOXYX
C    3,THETW(L,NY,NX),COXYS(L,NY,NX),COXYG(L,NY,NX)
C    4,ROXYD,VMXD4,VMXDXS,VMXDXB,VMXD4S,VMXD4B,FNO2S,FNO2B
C    5,ZNFN4S,ZNFN4B
7777  FORMAT(A8,5I4,50E12.4)
C     ENDIF
      ELSE
      RDNO3(NGL,N,K)=0.0
      RDNOB(NGL,N,K)=0.0
      RDNO2(NGL,N,K)=0.0
      RDN2B(NGL,N,K)=0.0
      RDN2O(NGL,N,K)=0.0
      RGOMY(NGL,N,K)=0.0
      RGOMD(NGL,N,K)=0.0
      ENDIF
C
C     BIOMASS DECOMPOSITION AND MINERALIZATION
C
C     MINERALIZATION-IMMOBILIZATION OF NH4 IN SOIL FROM MICROBIAL
C     C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
C
C     RINHP=NH4 mineralization (-ve) or immobilization (+ve) demand
C     OMC,OMN=microbial nonstructural C,N
C     CNOMC=maximum microbial N:C ratio
C     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
C     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
C     minimum NH4 concentration and Km for NH4 uptake
C     RINHX=microbially limited NH4 demand
C     BIOA=microbial surface area, OMA=active biomass
C     TFNG=temp+water stress
C     FNH4S,FNHBS=fractions of NH4 in non-band, band
C     RINHO,RINHB=substrate-unlimited NH4 mineraln-immobiln
C     VOLW=water content
C     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band
C     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
C     RINH4,RINB4=substrate-limited NH4 mineraln-immobiln in non-band, band
C     TRINH4=total NH4 net mineraln (-ve) or immobiln (+ve)
C

      RINHP=(OMC(3,NGL,N,K,L,NY,NX)*CNOMC(3,NGL,N,K)
     2-OMN(3,NGL,N,K,L,NY,NX))
      IF(RINHP.GT.0.0)THEN
      CNH4X=AMAX1(0.0,CNH4S(L,NY,NX)-Z4MN)
      CNH4Y=AMAX1(0.0,CNH4B(L,NY,NX)-Z4MN)
      RINHX=AMIN1(RINHP,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*Z4MX)
      RINHO(NGL,N,K,L,NY,NX)=FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
      RINHB(NGL,N,K,L,NY,NX)=FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
      ZNH4M=Z4MN*VOLW(L,NY,NX)*FNH4S
      ZNHBM=Z4MN*VOLW(L,NY,NX)*FNHBS
      RINH4(NGL,N,K)=AMIN1(FNH4X*AMAX1(0.0,(ZNH4S(L,NY,NX)-ZNH4M))
     2,RINHO(NGL,N,K,L,NY,NX))
      RINB4(NGL,N,K)=AMIN1(FNB4X*AMAX1(0.0,(ZNH4B(L,NY,NX)-ZNHBM))
     2,RINHB(NGL,N,K,L,NY,NX))
      ELSE
      RINHO(NGL,N,K,L,NY,NX)=0.0
      RINHB(NGL,N,K,L,NY,NX)=0.0
      RINH4(NGL,N,K)=RINHP*FNH4S
      RINB4(NGL,N,K)=RINHP*FNHBS
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+(RINH4(NGL,N,K)+RINB4(NGL,N,K))
C    2/AREA(3,L,NY,NX)
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,7776)'RINH4',I,J,NX,NY,L,K,N,RINH4(NGL,N,K),RINHP
C    1,BIOA*OMA(NGL,N,K)*Z4MX*TFNG(NGL,N,K),BIOA,OMA(NGL,N,K),Z4MX,TFNG(NGL,N,K)
C    2,OMC(3,NGL,N,K,L,NY,NX),CNOMC(3,NGL,N,K),OMN(3,NGL,N,K,L,NY,NX)
C    3,RINHO(NGL,N,K,L,NY,NX),CNH4S(L,NY,NX),FNH4X,ZNH4S(L,NY,NX)
C    4,ZNH4B(L,NY,NX),ZNH4T(L),OQN(K,L,NY,NX),TRINH4(NY,NX)
7776  FORMAT(A8,7I4,30E12.4)
C     ENDIF
C
C     MINERALIZATION-IMMOBILIZATION OF NO3 IN SOIL FROM MICROBIAL
C     C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
C
C     RINOP=NO3 immobilization (+ve) demand
C     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
C     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
C     min NO3 concentration and Km for NO3 uptake
C     RINOX=microbially limited NO3 demand
C     BIOA=microbial surface area, OMA=active biomass
C     TFNG=temp+water stress
C     FNO3S,FNO3B=fractions of NO3 in non-band, band
C     RINOO,RINOB=substrate-unlimited NO3 immobiln
C     VOLW=water content
C     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band
C     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
C     RINO3,RINB3=substrate-limited NO3 immobiln in non-band, band
C     TRINH4=total net NH4+NO3 mineraln (-ve) or immobiln (+ve)
C
      RINOP=AMAX1(0.0,RINHP-RINH4(NGL,N,K)-RINB4(NGL,N,K))
      IF(RINOP.GT.0.0)THEN
      CNO3X=AMAX1(0.0,CNO3S(L,NY,NX)-ZOMN)
      CNO3Y=AMAX1(0.0,CNO3B(L,NY,NX)-ZOMN)
      RINOX=AMIN1(RINOP,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*ZOMX)
      RINOO(NGL,N,K,L,NY,NX)=FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
      RINOB(NGL,N,K,L,NY,NX)=FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
      ZNO3M=ZOMN*VOLW(L,NY,NX)*FNO3S
      ZNOBM=ZOMN*VOLW(L,NY,NX)*FNO3B
      RINO3(NGL,N,K)=AMIN1(FNO3X*AMAX1(0.0,(ZNO3S(L,NY,NX)-ZNO3M))
     2,RINOO(NGL,N,K,L,NY,NX))
      RINB3(NGL,N,K)=AMIN1(FNB3X*AMAX1(0.0,(ZNO3B(L,NY,NX)-ZNOBM))
     2,RINOB(NGL,N,K,L,NY,NX))
      ELSE
      RINOO(NGL,N,K,L,NY,NX)=0.0
      RINOB(NGL,N,K,L,NY,NX)=0.0
      RINO3(NGL,N,K)=RINOP*FNO3S
      RINB3(NGL,N,K)=RINOP*FNO3B
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+(RINO3(NGL,N,K)+RINB3(NGL,N,K))
C     IF(RINO3(NGL,N,K).LT.0.0.OR.RINB3(NGL,N,K).LT.0.0)THEN
C     WRITE(*,4321)'RINO3',I,J,NX,NY,L,K,N
C    2,RINOO(NGL,N,K,L,NY,NX),RINO3(NGL,N,K)
C    2,RINOP,BIOA,OMA(NGL,N,K),TFNG(NGL,N,K),ZOMX,WFN(NGL,N,K),FNO3X,FNO3B
C    2,VLNO3(L,NY,NX),VLNOB(L,NY,NX),CNO3S(L,NY,NX),CNO3B(L,NY,NX)
C    3,RINOB(NGL,N,K,L,NY,NX),RINB3,ZNO3S(L,NY,NX),ZNO3B(L,NY,NX)
C    3,OMC(3,NGL,N,K,L,NY,NX),CPOMC(3,NGL,N,K),OMP(3,NGL,N,K,L,NY,NX),WFNG
4321  FORMAT(A8,7I4,60E12.4)
C     ENDIF
C
C     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SOIL FROM MICROBIAL
C     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
C
C     RIPOP=H2PO4 mineralization (-ve) or immobilization (+ve) demand
C     OMC,OMP=microbial nonstructural C,P
C     CPOMC=maximum microbial P:C ratio
C     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
C     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
C     min H2PO4 concentration and Km for H2PO4 uptake
C     RIPOX=microbially limited H2PO4 demand
C     BIOA=microbial surface area, OMA=active biomass
C     TFNG=temp+water stress
C     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
C     RIPOO,RIPBO=substrate-unlimited H2PO4 mineraln-immobiln
C     H2POM,H2PBM=H2PO4 not available for uptake in non-band, band
C     VOLW=water content
C     FPO4X,FPOBX=fractions of biol H2PO4 demand in non-band, band
C     RIPO4,RIPOB=substrate-limited H2PO4 mineraln-immobn in non-band, band
C     TRIPO4=total H2PO4 net mineraln (-ve) or immobiln (+ve)
C
      RIPOP=(OMC(3,NGL,N,K,L,NY,NX)*CPOMC(3,NGL,N,K)
     2-OMP(3,NGL,N,K,L,NY,NX))
      IF(RIPOP.GT.0.0)THEN
      CH2PX=AMAX1(0.0,CH2P4(L,NY,NX)-HPMN)
      CH2PY=AMAX1(0.0,CH2P4B(L,NY,NX)-HPMN)
      RIPOX=AMIN1(RIPOP,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX)
      RIPOO(NGL,N,K,L,NY,NX)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU)
      RIPBO(NGL,N,K,L,NY,NX)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU)
      H2POM=HPMN*VOLW(L,NY,NX)*FH2PS
      H2PBM=HPMN*VOLW(L,NY,NX)*FH2PB
      RIPO4(NGL,N,K)=AMIN1(FPO4X*AMAX1(0.0,(H2PO4(L,NY,NX)-H2POM))
     2,RIPOO(NGL,N,K,L,NY,NX))
      RIPOB(NGL,N,K)=AMIN1(FPOBX*AMAX1(0.0,(H2POB(L,NY,NX)-H2PBM))
     2,RIPBO(NGL,N,K,L,NY,NX))
      ELSE
      RIPOO(NGL,N,K,L,NY,NX)=0.0
      RIPBO(NGL,N,K,L,NY,NX)=0.0
      RIPO4(NGL,N,K)=RIPOP*FH2PS
      RIPOB(NGL,N,K)=RIPOP*FH2PB
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+(RIPO4(NGL,N,K)+RIPOB(NGL,N,K))
C     IF(NY.EQ.5.AND.L.EQ.10.AND.K.EQ.3.AND.N.EQ.2)THEN
C     WRITE(*,4322)'RIPO4',I,J,NX,NY,L,K,N,RIPO4(NGL,N,K),FPO4X,H2P4T(L)
C    2,RIPOO(NGL,N,K,L,NY,NX),RIPOP,BIOA,OMA(NGL,N,K),TFNG(NGL,N,K),HPMX,WFN(NGL,N,K)
C    2,VLPO4(L,NY,NX),VLPOB(L,NY,NX),CH2P4(L,NY,NX),CH2P4B(L,NY,NX)
C    3,OMC(3,NGL,N,K,L,NY,NX),CPOMC(3,NGL,N,K),OMP(3,NGL,N,K,L,NY,NX),WFNG
4322  FORMAT(A8,7I4,30E12.4)
C     ENDIF
C
C     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SOIL FROM MICROBIAL
C     C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
C
C     RIP1P=HPO4 mineralization (-ve) or immobilization (+ve) demand
C     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
C     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
C     min HPO4 concentration and Km for HPO4 uptake
C     RIP1X=microbially limited HPO4 demand
C     BIOA=microbial surface area, OMA=active biomass
C     TFNG=temp+water stress
C     FH1PS,FH1PB=fractions of HPO4 in non-band, band
C     RIPO1,RIPB1=substrate-unlimited HPO4 mineraln-immobiln
C     H1POM,H1PBM=HPO4 not available for uptake in non-band, band
C     VOLW=water content
C     FP14X,FP1BX=fractions of biol HPO4 demand in non-band, band
C     RIP14,RIP1B=substrate-limited HPO4 mineraln-immobn in non-band, band
C     TRIPO4=total H2PO4+HPO4 net mineraln (-ve) or immobiln (+ve)
C
      RIP1P=0.1*AMAX1(0.0,RIPOP-RIPO4(NGL,N,K)-RIPOB(NGL,N,K))
      IF(RIP1P.GT.0.0)THEN
      CH1PX=AMAX1(0.0,CH1P4(L,NY,NX)-HPMN)
      CH1PY=AMAX1(0.0,CH1P4B(L,NY,NX)-HPMN)
      RIP1X=AMIN1(RIP1P,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)*HPMX)
      RIPO1(NGL,N,K,L,NY,NX)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU)
      RIPB1(NGL,N,K,L,NY,NX)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU)
      H1POM=HPMN*VOLW(L,NY,NX)*FH1PS
      H1PBM=HPMN*VOLW(L,NY,NX)*FH1PB
      RIP14(NGL,N,K)=AMIN1(FP14X*AMAX1(0.0,(H1PO4(L,NY,NX)-H1POM))
     2,RIPO1(NGL,N,K,L,NY,NX))
      RIP1B(NGL,N,K)=AMIN1(FP1BX*AMAX1(0.0,(H1POB(L,NY,NX)-H1PBM))
     2,RIPB1(NGL,N,K,L,NY,NX))
      ELSE
      RIPO1(NGL,N,K,L,NY,NX)=0.0
      RIPB1(NGL,N,K,L,NY,NX)=0.0
      RIP14(NGL,N,K)=RIP1P*FH1PS
      RIP1B(NGL,N,K)=RIP1P*FH1PB
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+(RIP14(NGL,N,K)+RIP1B(NGL,N,K))
C     IF(NY.EQ.5.AND.L.EQ.10.AND.K.EQ.3.AND.N.EQ.2)THEN
C     WRITE(*,4323)'RIP14',I,J,NX,NY,L,K,N,RIP14(NGL,N,K),FP14X,H1P4T(L)
C    2,RIPO1(NGL,N,K,L,NY,NX),RIP1P,BIOA,OMA(NGL,N,K),TFNG(NGL,N,K),HPMX,WFN(NGL,N,K)
C    2,VLPO4(L,NY,NX),VLPOB(L,NY,NX),CH1P4(L,NY,NX),CH1P4B(L,NY,NX)
C    3,OMC(3,NGL,N,K,L,NY,NX),CPOMC(3,NGL,N,K),OMP(3,NGL,N,K,L,NY,NX),WFNG
4323  FORMAT(A8,7I4,30E12.4)
C     ENDIF
C
C     MINERALIZATION-IMMOBILIZATION OF NH4 IN SURFACE RESIDUE FROM
C     MICROBIAL C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL
C     ZONES OF SOIL SURFACE
C
C     RINHPR=NH4 mineralization (-ve) or immobilization (+ve) demand
C     NU=surface layer number
C     CNH4S,CNH4B=aqueous NH4 concentrations in non-band, band
C     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate,
C     minimum NH4 concentration and Km for NH4 uptake
C     BIOA=microbial surface area, OMA=active biomass
C     TFNG=temp+water stress
C     FNH4S,FNHBS=fractions of NH4 in non-band, band
C     RINHOR=substrate-unlimited NH4 mineraln-immobiln
C     VOLW=water content
C     ZNH4M=NH4 not available for uptake
C     FNH4XR=fractions of biological NH4 demand
C     RINH4R=substrate-limited NH4 mineraln-immobiln
C     TRINH4=total NH4 net mineraln (-ve) or immobiln (+ve)
C
      IF(L.EQ.0)THEN
      RINHPR=RINHP-RINH4(NGL,N,K)-RINO3(NGL,N,K)
      IF(RINHPR.GT.0.0)THEN
      CNH4X=AMAX1(0.0,CNH4S(NU(NY,NX),NY,NX)-Z4MN)
      CNH4Y=AMAX1(0.0,CNH4B(NU(NY,NX),NY,NX)-Z4MN)
      RINHOR(NGL,N,K,NY,NX)=AMIN1(RINHPR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)
     2*Z4MX)*(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y/(CNH4Y+Z4KU))
      ZNH4M=Z4MN*VOLW(NU(NY,NX),NY,NX)
      RINH4R(NGL,N,K)=AMIN1(FNH4XR(NGL,N,K)*AMAX1(0.0
     2,(ZNH4T(NU(NY,NX))-ZNH4M)),RINHOR(NGL,N,K,NY,NX))
      ELSE
      RINHOR(NGL,N,K,NY,NX)=0.0
      RINH4R(NGL,N,K)=RINHPR
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+RINH4R(NGL,N,K)
C    2/AREA(3,L,NY,NX)
C     IF(K.EQ.2.AND.N.EQ.1)THEN
C     WRITE(*,7778)'RINH4R',I,J,NX,NY,L,K,N,RINH4R(NGL,N,K),RINHPR
C    2,BIOA*OMA(NGL,N,K)*Z4MX,RINHP,RINH4(NGL,N,K),RINO3(NGL,N,K)
C    3,RINHOR(NGL,N,K,NY,NX),CNH4S(NU(NY,NX),NY,NX),FNH4XR(NGL,N,K)
C    4,ZNH4T(NU(NY,NX))
7778  FORMAT(A8,7I4,20E12.4)
C     ENDIF
C
C     MINERALIZATION-IMMOBILIZATION OF NO3 IN SURFACE RESIDUE FROM
C     MICROBIAL C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL
C     ZONES OF SOIL SURFACE
C
C     RINOPR=NH4 mineralization (-ve) or immobilization (+ve) demand
C     NU=surface layer number
C     CNO3S,CNO3B=aqueous NO3 concentrations in non-band, band
C     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate,
C     minimum NO3 concentration and Km for NO3 uptake
C     RINOOR=microbially limited NO3 demand
C     BIOA=microbial surface area, OMA=active biomass
C     TFNG=temp+water stress
C     FNO3S,FNO3B=fractions of NO3 in non-band, band
C     RINO3R=substrate-unlimited NO3 immobiln
C     VOLW=water content
C     ZNO3M=NO3 not available for uptake
C     FNO3XR=fraction of biological NO3 demand
C     RINO3R=substrate-limited NO3 immobiln
C     TRINH4=total NH4+NO3 net mineraln (-ve) or immobiln (+ve)
C
      RINOPR=AMAX1(0.0,RINHPR-RINH4R(NGL,N,K))
      IF(RINOPR.GT.0.0)THEN
      CNO3X=AMAX1(0.0,CNO3S(NU(NY,NX),NY,NX)-ZOMN)
      CNO3Y=AMAX1(0.0,CNO3B(NU(NY,NX),NY,NX)-ZOMN)
      RINOOR(NGL,N,K,NY,NX)=AMAX1(RINOPR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)
     2*ZOMX)*(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y/(CNO3Y+ZOKU))
      ZNO3M=ZOMN*VOLW(NU(NY,NX),NY,NX)
      RINO3R(NGL,N,K)=AMIN1(FNO3XR(NGL,N,K)*AMAX1(0.0
     2,(ZNO3T(NU(NY,NX))-ZNO3M)),RINOOR(NGL,N,K,NY,NX))
      ELSE
      RINOOR(NGL,N,K,NY,NX)=0.0
      RINO3R(NGL,N,K)=RINOPR
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+RINO3R(NGL,N,K)
C
C     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SURFACE RESIDUE FROM
C     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
C     ZONES OF SOIL SURFACE
C
C     RIPOPR=H2PO4 mineralization (-ve) or immobilization (+ve) demand
C     NU=surface layer number
C     CH2P4,CH2P4B=aqueous H2PO4 concentrations in non-band, band
C     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate,
C     minimum H2PO4 concentration and Km for H2PO4 uptake
C     RIPOOR=microbially limited H2PO4 demand
C     BIOA=microbial surface area, OMA=active biomass
C     TFNG=temp+water stress
C     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
C     RIPOOR=substrate-unlimited H2PO4 mineraln-immobiln
C     VOLW=water content
C     H2P4M=H2PO4 not available for uptake
C     FPO4XR=fractions of biological H2PO4 demand
C     RIPO4R=substrate-limited H2PO4 mineraln-immobiln
C     TRIPO4=total H2PO4 net mineraln (-ve) or immobiln (+ve)
C
      RIPOPR=RIPOP-RIPO4(NGL,N,K)
      IF(RIPOPR.GT.0.0)THEN
      CH2PX=AMAX1(0.0,CH2P4(NU(NY,NX),NY,NX)-HPMN)
      CH2PY=AMAX1(0.0,CH2P4B(NU(NY,NX),NY,NX)-HPMN)
      RIPOOR(NGL,N,K,NY,NX)=AMIN1(RIPOPR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)
     2*HPMX)*(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY/(CH2PY+HPKU))
      H2P4M=HPMN*VOLW(NU(NY,NX),NY,NX)
      RIPO4R(NGL,N,K)=AMIN1(FPO4XR(NGL,N,K)*AMAX1(0.0
     2,(H2P4T(NU(NY,NX))-H2P4M)),RIPOOR(NGL,N,K,NY,NX))
      ELSE
      RIPOOR(NGL,N,K,NY,NX)=0.0
      RIPO4R(NGL,N,K)=RIPOPR
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+RIPO4R(NGL,N,K)
C     WRITE(*,7778)'RIPO4R',I,J,NX,NY,L,K,N,RIPO4R(NGL,N,K),FPO4XR(NGL,N,K)
C    2,H2P4T(NU(NY,NX)),RIPOOR(NGL,N,K,NY,NX),RIPOPR
C
C     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SURFACE RESIDUE FROM
C     MICROBIAL C:P AND PO4 CONCENTRATION IN BAND AND NON-BAND SOIL
C     ZONES OF SOIL SURFACE
C
C     RIP1PR=HPO4 mineralization (-ve) or immobilization (+ve) demand
C     NU=surface layer number
C     CH1P4,CH1P4B=aqueous HPO4 concentrations in non-band, band
C     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate,
C     minimum HPO4 concentration and Km for HPO4 uptake
C     RIPO1R=microbially limited HPO4 demand
C     BIOA=microbial surface area, OMA=active biomass
C     TFNG=temp+water stress
C     FH1PS,FH1PB=fractions of HPO4 in non-band, band
C     RIPO1R=substrate-unlimited HPO4 mineraln-immobiln
C     VOLW=water content
C     H1P4M=HPO4 not available for uptake
C     FP14XR=fraction of biological HPO4 demand
C     RIP14R=substrate-limited HPO4 minereraln-immobiln
C     TRIPO4=total HPO4 net mineraln (-ve) or immobiln (+ve)
C
      RIP1PR=0.1*AMAX1(0.0,RIPOPR-RIPO4R(NGL,N,K))
      IF(RIP1PR.GT.0.0)THEN
      CH1PX=AMAX1(0.0,CH1P4(NU(NY,NX),NY,NX)-HPMN)
      CH1PY=AMAX1(0.0,CH1P4B(NU(NY,NX),NY,NX)-HPMN)
      RIPO1R(NGL,N,K,NY,NX)=AMIN1(RIP1PR,BIOA*OMA(NGL,N,K)*TFNG(NGL,N,K)
     2*HPMX)*(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY/(CH1PY+HPKU))
      H1P4M=HPMN*VOLW(NU(NY,NX),NY,NX)
      RIP14R(NGL,N,K)=AMIN1(FP14XR(NGL,N,K)*AMAX1(0.0
     2,(H1P4T(NU(NY,NX))-H1P4M)),RIPO1R(NGL,N,K,NY,NX))
      ELSE
      RIPO1R(NGL,N,K,NY,NX)=0.0
      RIP14R(NGL,N,K)=RIP1PR
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+RIP14R(NGL,N,K)
C     WRITE(*,7778)'RIP14R',I,J,NX,NY,L,K,N,RIP14R(NGL,N,K),FP14XR(NGL,N,K)
C    2,H1P4T(NU(NY,NX)),RIPO1R(NGL,N,K,NY,NX),RIP1PR
      ENDIF
C
C     pH EFFECT ON MAINTENANCE RESPIRATION
C
C     FPH=pH effect on maintenance respiration
C     RMOM=specific maintenance respiration rate
C     TFNR=temperature effect on maintenance respiration
C     OMN=microbial N biomass
C     RMOMK=effect of low microbial C concentration on mntc respn
C
      FPH=1.0+AMAX1(0.0,0.25*(6.5-PH(L,NY,NX)))
      RMOMX=RMOM*TFNR(NGL,N,K)*FPH
      RMOMC(1,NGL,N,K)=OMN(1,NGL,N,K,L,NY,NX)*RMOMX*RMOMK(1)
      RMOMC(2,NGL,N,K)=OMN2(NGL,N,K)*RMOMX*RMOMK(2)

C
C     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
C
C     RMOMT=total maintenance respiration
C     RGOMT=growth respiration
C     RXOMT=senescence respiration
C
      RMOMT=RMOMC(1,NGL,N,K)+RMOMC(2,NGL,N,K)
      RGOMT=AMAX1(0.0,RGOMO(NGL,N,K)-RMOMT)
      RXOMT=AMAX1(0.0,RMOMT-RGOMO(NGL,N,K))
C
C     N2 FIXATION: N=(6) AEROBIC, (7) ANAEROBIC
C     FROM GROWTH RESPIRATION, FIXATION ENERGY REQUIREMENT,
C     MICROBIAL N REQUIREMENT IN LABILE (1) AND
C     RESISTANT (2) FRACTIONS
C
C     RGN2P=respiration to meet N2 fixation demand
C     OMC,OMN=microbial nonstructural C,N
C     CNOMC=maximum microbial N:C ratio
C     EN2F=N2 fixation yield per unit nonstructural C
C     RGOMT=growth respiration
C     RGN2F=respiration for N2 fixation
C     CZ2GS=aqueous N2 concentration
C     ZFKM=Km for N2 uptake
C     OMGR*OMC(3,NGL,N,K,L,NY,NX)=nonstructural C limitation to RGN2F
C     RN2FX=N2 fixation rate
C
      IF(K.LE.4.AND.(N.EQ.6.OR.N.EQ.7))THEN
      RGN2P=AMAX1(0.0,OMC(3,NGL,N,K,L,NY,NX)*CNOMC(3,NGL,N,K)
     2-OMN(3,NGL,N,K,L,NY,NX))/EN2F(N)
      IF(RGOMT.GT.ZEROS(NY,NX))THEN
      RGN2F(NGL,N,K)=AMIN1(RGOMT*RGN2P/(RGOMT+RGN2P)
     2*CZ2GS(L,NY,NX)/(CZ2GS(L,NY,NX)+ZFKM)
     2,OMGR*OMC(3,NGL,N,K,L,NY,NX))
      ELSE
      RGN2F(NGL,N,K)=0.0
      ENDIF
      RN2FX(NGL,N,K)=RGN2F(NGL,N,K)*EN2F(N)
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,5566)'N2 FIX',I,J,NX,NY,L,K,N,RN2FX(NGL,N,K),EN2F(N)
C    2,OMC(3,NGL,N,K,L,NY,NX)*CNOMC(3,NGL,N,K),OMN(3,NGL,N,K,L,NY,NX)
C    3,RINH4(NGL,N,K),RINO3(NGL,N,K),RGN2P,RGN2F(NGL,N,K),RGOMT
C    4,CZ2GS(L,NY,NX)
5566  FORMAT(A8,7I4,30E12.4)
C     ENDIF
      ELSE
      RN2FX(NGL,N,K)=0.0
      RGN2F(NGL,N,K)=0.0
      ENDIF
C
C     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
C     FROM O2, NOX AND C REDUCTION
C
C     CGOMX=DOC+acetate uptake from aerobic growth respiration
C     CGOMD=DOC+acetate uptake from denitrifier growth respiration
C     RMOMT=maintenance respiration
C     RGOMO=total respiration
C     RGOMD=respiration for denitrifcation
C     RGN2F=respiration for N2 fixation
C     ECHZ,ENOX=growth respiration efficiencies for O2, NOx reduction
C     CGOMC,CGOQC,CGOAC=total DOC+acetate, DOC, acetate uptake(heterotrophs
C     CGOMC=total CO2,CH4 uptake (autotrophs)
C     CGOMN,CGOMP=DON, DOP uptake
C     FGOCP,FGOAP=DOC,acetate/(DOC+acetate)
C     OQN,OPQ=DON,DOP
C     FOMK=faction of OMA in total OMA
C     CNQ,CPQ=DON/DOC, DOP/DOC
C     FCN,FCP=limitation from N,P
C
      CGOMX=AMIN1(RMOMT,RGOMO(NGL,N,K))+RGN2F(NGL,N,K)
     2+(RGOMT-RGN2F(NGL,N,K))/ECHZ
      CGOMD=RGOMD(NGL,N,K)/ENOX
      CGOMC(NGL,N,K)=CGOMX+CGOMD
      IF(K.LE.4)THEN
      CGOQC(NGL,N,K)=CGOMX*FGOCP+CGOMD
      CGOAC(NGL,N,K)=CGOMX*FGOAP
      CGOXC=CGOQC(NGL,N,K)+CGOAC(NGL,N,K)
      CGOMN(NGL,N,K)=AMAX1(0.0,AMIN1(OQN(K,L,NY,NX)*FOMK(NGL,N,K)
     2,CGOXC*CNQ(K)/FCN(NGL,N,K)))
      CGOMP(NGL,N,K)=AMAX1(0.0,AMIN1(OQP(K,L,NY,NX)*FOMK(NGL,N,K)
     2,CGOXC*CPQ(K)/FCP(NGL,N,K)))
      ELSE
      CGOQC(NGL,N,K)=CGOMX+CGOMD
      CGOAC(NGL,N,K)=0.0
      CGOMN(NGL,N,K)=0.0
      CGOMP(NGL,N,K)=0.0
      ENDIF
      TCGOQC(K)=TCGOQC(K)+CGOQC(NGL,N,K)
      TCGOAC(K)=TCGOAC(K)+CGOAC(NGL,N,K)
      TCGOMN(K)=TCGOMN(K)+CGOMN(NGL,N,K)
      TCGOMP(K)=TCGOMP(K)+CGOMP(NGL,N,K)
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.3)THEN
C     WRITE(*,5557)'CGOQC',I,J,NX,NY,L,K,N,CGOQC(NGL,N,K),CGOMX
C    2,FGOCP,FGOAP,CGOMD,RMOMT,RGN2F(NGL,N,K),ECHZ
C    3,RGOMD(NGL,N,K),ENOX,RGOMO(NGL,N,K),WFN(NGL,N,K),FOXYX
C     WRITE(*,5557)'CGOMP',I,J,NX,NY,L,K,N,CGOMP(NGL,N,K),OQP(K,L,NY,NX)
C    2,FOMK(NGL,N,K),CGOXC,CPQ(K),FCP(NGL,N,K),CGOQC(NGL,N,K),CGOAC(NGL,N,K)
5557  FORMAT(A8,7I4,30E12.4)
C     ENDIF
C
C     TRANSFER UPTAKEN C,N,P FROM STORAGE TO ACTIVE BIOMASS
C
C     OMC,OMN,OMP=nonstructural C,N,P
C     CCC,CNC,CPC=C:N:P ratios used to calculate C,N,P recycling
C     CNOMC,CPOMC=maximum microbial N:C, P:C ratios
C     RCCC,RCCN,RCCP=C,N,P recycling fractions
C     RCCZ,RCCY=min, max C recycling fractions
C     RCCX,RCCQ=max N,P recycling fractions
C
      IF(OMC(3,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX)
     2.AND.OMC(1,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CCC=AMAX1(0.0,AMIN1(1.0
     2,OMN(3,NGL,N,K,L,NY,NX)/(OMN(3,NGL,N,K,L,NY,NX)
     2+OMC(3,NGL,N,K,L,NY,NX)*CNOMC(3,NGL,N,K))
     3,OMP(3,NGL,N,K,L,NY,NX)/(OMP(3,NGL,N,K,L,NY,NX)
     4+OMC(3,NGL,N,K,L,NY,NX)*CPOMC(3,NGL,N,K))))
      CXC=OMC(3,NGL,N,K,L,NY,NX)/OMC(1,NGL,N,K,L,NY,NX)
      C3C=1.0/(1.0+CXC/CKC)
      CNC=AMAX1(0.0,AMIN1(1.0
     2,OMC(3,NGL,N,K,L,NY,NX)/(OMC(3,NGL,N,K,L,NY,NX)
     2+OMN(3,NGL,N,K,L,NY,NX)/CNOMC(3,NGL,N,K))))
      CPC=AMAX1(0.0,AMIN1(1.0
     2,OMC(3,NGL,N,K,L,NY,NX)/(OMC(3,NGL,N,K,L,NY,NX)
     3+OMP(3,NGL,N,K,L,NY,NX)/CPOMC(3,NGL,N,K))))
      RCCC=RCCZ+AMAX1(CCC,C3C)*RCCY

      RCCN=CNC*RCCX
      RCCP=CPC*RCCQ
      ELSE
      RCCC=RCCZ
      RCCN=0.0
      RCCP=0.0
      ENDIF
C     IF((I/120)*120.EQ.I.AND.J.EQ.24)THEN
C     WRITE(*,5555)'RCCC',I,J,NX,NY,L,K,N,RCCC,RCCN,RCCP
C    2,OMC(3,NGL,N,K,L,NY,NX),OMN(3,NGL,N,K,L,NY,NX),OMP(3,NGL,N,K,L,NY,NX)
C    3,CCC,C3C,CNC,CPC
C     ENDIF
C
C     MICROBIAL ASSIMILATION OF NONSTRUCTURAL C,N,P
C
C     CGOMZ=transfer from nonstructural to structural microbial C
C     TFNG=temperature+water stress function
C     OMGR=rate constant for transferring nonstructural to structural C
C     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
C     FL=partitioning between labile and resistant microbial components
C     OMC,OMN,OMP=nonstructural microbial C,N,P
C
      CGOMZ=TFNG(NGL,N,K)*OMGR*AMAX1(0.0,OMC(3,NGL,N,K,L,NY,NX))

      DO 745 M=1,2
      CGOMS(M,NGL,N,K)=FL(M)*CGOMZ
      IF(OMC(3,NGL,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CGONS(M,NGL,N,K)=AMIN1(FL(M)*AMAX1(0.0,OMN(3,NGL,N,K,L,NY,NX))
     2,CGOMS(M,NGL,N,K)*OMN(3,NGL,N,K,L,NY,NX)/OMC(3,NGL,N,K,L,NY,NX))
      CGOPS(M,NGL,N,K)=AMIN1(FL(M)*AMAX1(0.0,OMP(3,NGL,N,K,L,NY,NX))
     2,CGOMS(M,NGL,N,K)*OMP(3,NGL,N,K,L,NY,NX)/OMC(3,NGL,N,K,L,NY,NX))
      ELSE
      CGONS(M,NGL,N,K)=0.0
      CGOPS(M,NGL,N,K)=0.0
      ENDIF
C
C     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
C     RATE, TEMPERATURE
C
C     SPOMX=rate constant for microbial decomposition
C     SPOMC=basal decomposition rate
C     SPOMK=effect of low microbial C concentration on microbial decay
C     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
C     RDOMC,RDOMN,RDOMP=microbial C,N,P litterfall
C     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
C
      SPOMX=SQRT(TFNG(NGL,N,K))*SPOMC(M)*SPOMK(M)
      RXOMC(M,NGL,N,K)=AMAX1(0.0,OMC(M,NGL,N,K,L,NY,NX)*SPOMX)
      RXOMN(M,NGL,N,K)=AMAX1(0.0,OMN(M,NGL,N,K,L,NY,NX)*SPOMX)
      RXOMP(M,NGL,N,K)=AMAX1(0.0,OMP(M,NGL,N,K,L,NY,NX)*SPOMX)
      RDOMC(M,NGL,N,K)=RXOMC(M,NGL,N,K)*(1.0-RCCC)
      RDOMN(M,NGL,N,K)=RXOMN(M,NGL,N,K)*(1.0-RCCC)*(1.0-RCCN)
      RDOMP(M,NGL,N,K)=RXOMP(M,NGL,N,K)*(1.0-RCCC)*(1.0-RCCP)
      R3OMC(M,NGL,N,K)=RXOMC(M,NGL,N,K)-RDOMC(M,NGL,N,K)
      R3OMN(M,NGL,N,K)=RXOMN(M,NGL,N,K)-RDOMN(M,NGL,N,K)
      R3OMP(M,NGL,N,K)=RXOMP(M,NGL,N,K)-RDOMP(M,NGL,N,K)
C
C     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
C     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
C
C     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
C     EHUM=humus transfer fraction from hour1.f
C     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
C
      RHOMC(M,NGL,N,K)=AMAX1(0.0,RDOMC(M,NGL,N,K)*EHUM(L,NY,NX))
      RHOMN(M,NGL,N,K)=AMAX1(0.0,RDOMN(M,NGL,N,K)*EHUM(L,NY,NX))
      RHOMP(M,NGL,N,K)=AMAX1(0.0,RDOMP(M,NGL,N,K)*EHUM(L,NY,NX))
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,8821)'RHOMC',I,J,L,K,N,M
C    2,RXOMC(M,NGL,N,K),RXOMN(M,NGL,N,K),RXOMP(M,NGL,N,K)
C    2,RDOMC(M,NGL,N,K),RDOMN(M,NGL,N,K),RDOMP(M,NGL,N,K)
C    2,R3OMC(M,NGL,N,K),R3OMN(M,NGL,N,K),R3OMP(M,NGL,N,K)
C    2,RHOMC(M,NGL,N,K),RHOMN(M,NGL,N,K),RHOMP(M,NGL,N,K)
C    4,OMC(M,NGL,N,K,L,NY,NX),OMN(M,NGL,N,K,L,NY,NX)
C    5,OMP(M,NGL,N,K,L,NY,NX)
C    4,OMC(3,NGL,N,K,L,NY,NX),OMN(3,NGL,N,K,L,NY,NX)
C    5,OMP(3,NGL,N,K,L,NY,NX)
C    6,OQC(K,L,NY,NX),OQN(K,L,NY,NX),OQP(K,L,NY,NX)
C    2,SPOMX,RCCC,RCCN,RCCP
C     ENDIF
C
C     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
C
      RCOMC(M,NGL,N,K)=RDOMC(M,NGL,N,K)-RHOMC(M,NGL,N,K)
      RCOMN(M,NGL,N,K)=RDOMN(M,NGL,N,K)-RHOMN(M,NGL,N,K)
      RCOMP(M,NGL,N,K)=RDOMP(M,NGL,N,K)-RHOMP(M,NGL,N,K)
745   CONTINUE
C
C     MICROBIAL DECOMPOSITION WHEN MAINTENANCE RESPIRATION
C     EXCEEDS UPTAKE
C
C     OMC,OMN,OMP=microbial C,N,P
C     RMOMT=total maintenance respiration
C     RXOMT=senescence respiration
C     RCCC=C recycling fraction
C     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
C     RMOMC=maintenance respiration
C     CNOMA,CPOMA=N:C,P:C ratios of active biomass
C     RDMMC,RDMMN,RDMMP=microbial C,N,P litterfall from senescence
C     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
C
      IF(RXOMT.GT.ZEROS(NY,NX).AND.RMOMT.GT.ZEROS(NY,NX)
     2.AND.RCCC.GT.ZERO)THEN
      FRM=RXOMT/RMOMT
      DO 730 M=1,2
      RXMMC(M,NGL,N,K)=AMIN1(OMC(M,NGL,N,K,L,NY,NX)
     2,AMAX1(0.0,FRM*RMOMC(M,NGL,N,K)/RCCC))
      RXMMN(M,NGL,N,K)=AMIN1(OMN(M,NGL,N,K,L,NY,NX)
     2,AMAX1(0.0,RXMMC(M,NGL,N,K)*CNOMA(NGL,N,K)))
      RXMMP(M,NGL,N,K)=AMIN1(OMP(M,NGL,N,K,L,NY,NX)
     2,AMAX1(0.0,RXMMC(M,NGL,N,K)*CPOMA(NGL,N,K)))
      RDMMC(M,NGL,N,K)=RXMMC(M,NGL,N,K)*(1.0-RCCC)
      RDMMN(M,NGL,N,K)=RXMMN(M,NGL,N,K)*(1.0-RCCN)*(1.0-RCCC)
      RDMMP(M,NGL,N,K)=RXMMP(M,NGL,N,K)*(1.0-RCCP)*(1.0-RCCC)
      R3MMC(M,NGL,N,K)=RXMMC(M,NGL,N,K)-RDMMC(M,NGL,N,K)
      R3MMN(M,NGL,N,K)=RXMMN(M,NGL,N,K)-RDMMN(M,NGL,N,K)
      R3MMP(M,NGL,N,K)=RXMMP(M,NGL,N,K)-RDMMP(M,NGL,N,K)
C
C     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
C     PRODUCTS
C
C     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
C     EHUM=humus transfer fraction
C     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
C
      RHMMC(M,NGL,N,K)=AMAX1(0.0,RDMMC(M,NGL,N,K)*EHUM(L,NY,NX))
      RHMMN(M,NGL,N,K)=AMAX1(0.0,RDMMN(M,NGL,N,K)*EHUM(L,NY,NX))
      RHMMP(M,NGL,N,K)=AMAX1(0.0,RDMMP(M,NGL,N,K)*EHUM(L,NY,NX))
      RCMMC(M,NGL,N,K)=RDMMC(M,NGL,N,K)-RHMMC(M,NGL,N,K)
      RCMMN(M,NGL,N,K)=RDMMN(M,NGL,N,K)-RHMMN(M,NGL,N,K)
      RCMMP(M,NGL,N,K)=RDMMP(M,NGL,N,K)-RHMMP(M,NGL,N,K)
C     IF(L.EQ.11.AND.K.EQ.1)THEN
C     WRITE(*,8821)'RCMMC',I,J,L,K,N,M,RCMMC(M,NGL,N,K)
C    2,RDMMC(M,NGL,N,K),RHMMC(M,NGL,N,K),OMC(M,N,K,L,NY,NX)
C    3,FRM,RMOMC(M,NGL,N,K),OMN(1,N,K,L,NY,NX),OMN2(NGL,N,K)
C    4,RMOM,TFNR(NGL,N,K),FPH,RDMMN(M,NGL,N,K),CNSHZ,RDMMP(M,NGL,N,K)
C    5,CPSHZ,EHUM(L,NY,NX),RXOMT,RMOMT,RMOMT,RGOMO(NGL,N,K)
C    6,RGOMP,WFN(NGL,N,K)
C     WRITE(*,8821)'RCMMP',I,J,L,K,N,M,RCMMP(M,NGL,N,K)
C    2,RDMMP(M,NGL,N,K),RHMMP(M,NGL,N,K),EHUM(L,NY,NX)
C    3,RCCC,RCCN,RCCP,RXMMP(M,NGL,N,K)
C     ENDIF
730   CONTINUE
      ELSE
      DO 720 M=1,2
      RXMMC(M,NGL,N,K)=0.0
      RXMMN(M,NGL,N,K)=0.0
      RXMMP(M,NGL,N,K)=0.0
      RDMMC(M,NGL,N,K)=0.0
      RDMMN(M,NGL,N,K)=0.0
      RDMMP(M,NGL,N,K)=0.0
      R3MMC(M,NGL,N,K)=0.0
      R3MMN(M,NGL,N,K)=0.0
      R3MMP(M,NGL,N,K)=0.0
      RHMMC(M,NGL,N,K)=0.0
      RHMMN(M,NGL,N,K)=0.0
      RHMMP(M,NGL,N,K)=0.0
      RCMMC(M,NGL,N,K)=0.0
      RCMMN(M,NGL,N,K)=0.0
      RCMMP(M,NGL,N,K)=0.0
720   CONTINUE
      ENDIF
      ELSE
      RUPOX(NGL,N,K)=0.0
      RGOMO(NGL,N,K)=0.0
      RCO2X(NGL,N,K)=0.0
      RCH3X(NGL,N,K)=0.0
      RCH4X(NGL,N,K)=0.0
      RGOMY(NGL,N,K)=0.0
      RGOMD(NGL,N,K)=0.0
      CGOMC(NGL,N,K)=0.0
      CGOMN(NGL,N,K)=0.0
      CGOMP(NGL,N,K)=0.0
      CGOQC(NGL,N,K)=0.0
      CGOAC(NGL,N,K)=0.0
      RDNO3(NGL,N,K)=0.0
      RDNOB(NGL,N,K)=0.0
      RDNO2(NGL,N,K)=0.0
      RDN2B(NGL,N,K)=0.0
      RDN2O(NGL,N,K)=0.0
      RN2FX(NGL,N,K)=0.0
      RINH4(NGL,N,K)=0.0
      RINO3(NGL,N,K)=0.0
      RIPO4(NGL,N,K)=0.0
      RIP14(NGL,N,K)=0.0
      RINB4(NGL,N,K)=0.0
      RINB3(NGL,N,K)=0.0
      RIPOB(NGL,N,K)=0.0
      RIP1B(NGL,N,K)=0.0
      IF(L.EQ.0)THEN
      RINH4R(NGL,N,K)=0.0
      RINO3R(NGL,N,K)=0.0
      RIPO4R(NGL,N,K)=0.0
      RIP14R(NGL,N,K)=0.0
      FNH4XR(NGL,N,K)=0.0
      FNO3XR(NGL,N,K)=0.0
      FPO4XR(NGL,N,K)=0.0
      FP14XR(NGL,N,K)=0.0
      ENDIF
      DO 725 M=1,2
      CGOMS(M,NGL,N,K)=0.0
      CGONS(M,NGL,N,K)=0.0
      CGOPS(M,NGL,N,K)=0.0
      RMOMC(M,NGL,N,K)=0.0
      RXMMC(M,NGL,N,K)=0.0
      RXMMN(M,NGL,N,K)=0.0
      RXMMP(M,NGL,N,K)=0.0
      RDMMC(M,NGL,N,K)=0.0
      RDMMN(M,NGL,N,K)=0.0
      RDMMP(M,NGL,N,K)=0.0
      R3MMC(M,NGL,N,K)=0.0
      R3MMN(M,NGL,N,K)=0.0
      R3MMP(M,NGL,N,K)=0.0
      RHMMC(M,NGL,N,K)=0.0
      RHMMN(M,NGL,N,K)=0.0
      RHMMP(M,NGL,N,K)=0.0
      RCMMC(M,NGL,N,K)=0.0
      RCMMN(M,NGL,N,K)=0.0
      RCMMP(M,NGL,N,K)=0.0
      RXOMC(M,NGL,N,K)=0.0
      RXOMN(M,NGL,N,K)=0.0
      RXOMP(M,NGL,N,K)=0.0
      RDOMC(M,NGL,N,K)=0.0
      RDOMN(M,NGL,N,K)=0.0
      RDOMP(M,NGL,N,K)=0.0
      R3OMC(M,NGL,N,K)=0.0
      R3OMN(M,NGL,N,K)=0.0
      R3OMP(M,NGL,N,K)=0.0
      RHOMC(M,NGL,N,K)=0.0
      RHOMN(M,NGL,N,K)=0.0
      RHOMP(M,NGL,N,K)=0.0
      RCOMC(M,NGL,N,K)=0.0
      RCOMN(M,NGL,N,K)=0.0
      RCOMP(M,NGL,N,K)=0.0
725   CONTINUE
      RH2GX(NGL,N,K)=0.0
      IF(K.EQ.5)THEN
      RVOXA(NGL,N)=0.0
      RVOXB(NGL,N)=0.0
      IF(N.EQ.5)THEN
      RH2GZ=0.0
      ENDIF
      ENDIF
      ENDIF
      ENDIF
750   CONTINUE
      ENDIF
760   CONTINUE
C
C     CHEMODENITRIFICATION
C
C     FNO2,FNB2=fraction of total NO2 demand in non-band,band
C     VMXC4S,VMXC4B=substrate-unlimited NO2 reduction in non-band,band
C     CHNO2,CHNOB=nitrous acid concentration in non-band,band
C     VOLWM=soil water content
C     FNO3S,FNO3B=fractions of NO2 in non-band,band
C     TFNX=temperature stress function
C     RCNO2,RCNOB=substrate-limited nitrous acid reduction in non-band,band
C     RCN2O,RCN2B=N2O production from nitrous acid reduction in non-band,band
C     RCNO3,RCN3B=NO3 production from nitrous acid reduction in non-band,band
C     RCOQN=DON production from nitrous acid reduction
C     RVMXC,RVMBC=demand for NO2 reduction in non-band,band
C
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMXC(L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=FMN*VLNO3(L,NY,NX)
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMBC(L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=FMN*VLNOB(L,NY,NX)
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
      VMXC4S=7.5E-02*CHNO2*VOLWM(NPH,L,NY,NX)*FNO3S*TFNX
      VMXC4B=7.5E-02*CHNOB*VOLWM(NPH,L,NY,NX)*FNO3B*TFNX
      RCNO2=AMAX1(0.0,AMIN1(ZNO2S(L,NY,NX)*FNO2,VMXC4S))
      RCNOB=AMAX1(0.0,AMIN1(ZNO2B(L,NY,NX)*FNB2,VMXC4B))
      RCN2O=0.10*RCNO2
      RCN2B=0.10*RCNOB
      RCNO3=0.80*RCNO2
      RCN3B=0.80*RCNOB
      RCOQN=0.10*(RCNO2+RCNOB)
      RVMXC(L,NY,NX)=VMXC4S
      RVMBC(L,NY,NX)=VMXC4B
C     IF((I/1)*1.EQ.I.AND.L.LE.5)THEN
C     WRITE(*,7779)'CHEMO',I,J,L,RCNO2,RCNOB,CHY1,CHNO2,CHNOB
C    2,CNO2S(L,NY,NX),CNO2B(L,NY,NX),VOLWM(NPH,L,NY,NX),FNO2
C    3,VMXC4S,VMXC4B,RVMXC(L,NY,NX),RNO2Y(L,NY,NX),RCN2O,RCN2B
C    4,RCNO3,RCNOB,RCOQN,VLNO3(L,NY,NX),VLNOB(L,NY,NX)
7779  FORMAT(A8,3I4,30E12.4)
C     ENDIF
C
C     DECOMPOSITION
C
C     ROQCK=total respiration of DOC+DOA used to represent microbial activity
C
      DO 1870 K=0,KL
      ROQCK(K)=0.0
      DO 1875 N=1,7
      DO 1875 NGL=1,JG
      ROQCK(K)=ROQCK(K)+ROQCD(NGL,N,K)
1875  CONTINUE
      XOQCK(K)=0.0
      XOQCZ(K)=0.0
      XOQNZ(K)=0.0
      XOQPZ(K)=0.0
      XOQAZ(K)=0.0
      DO 845 N=1,7
      DO 845 NGL=1,JG
      DO 845 M=1,3
      XOMCZ(M,NGL,N,K)=0.0
      XOMNZ(M,NGL,N,K)=0.0
      XOMPZ(M,NGL,N,K)=0.0
845   CONTINUE
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.LE.1)THEN
C     WRITE(*,4443)'PRIM1',I,J,NX,NY,L,K,ROQCK(K)
C    2,XOQCK(K),OQC(K,L,NY,NX),XOQCZ(K),OQN(K,L,NY,NX),XOQNZ(K)
C    3,OQP(K,L,NY,NX),XOQPZ(K),OQA(K,L,NY,NX),XOQAZ(K)
C     ENDIF
1870  CONTINUE
C
C     PRIMING of DOC,DON,DOP BETWEEN LITTER AND NON-LITTER C
C
C     OSRH=total SOC in each K
C     XFRK,XFRC,XFRN,XFRP,XFRA=transfer of respiration,DOC,DON,DOP,acetate
C     between each K and KK, FPRIM=priming transfer rate constant
C     TFND=temperature effect on priming transfers
C     ROQCK,OQC,OQN,OQP=respiration,DOC,DON,DOP
C     XOQCK,XOQCZ,XOQNZ,XOQPZ,XOQAZ=total XFRK,XFRC,XFRN,XFRP,XFRA for all K
C
      DO 795 K=0,KL
      IF(K.LE.KL-1)THEN
      DO 800 KK=K+1,KL
      OSRT=OSRH(K)+OSRH(KK)
      IF(OSRH(K).GT.ZEROS(NY,NX).AND.OSRH(KK).GT.ZEROS(NY,NX))THEN
      XFRK=FPRIM*TFND(L,NY,NX)*(ROQCK(K)*OSRH(KK)
     2-ROQCK(KK)*OSRH(K))/OSRT
      XFRC=FPRIM*TFND(L,NY,NX)*(OQC(K,L,NY,NX)*OSRH(KK)
     2-OQC(KK,L,NY,NX)*OSRH(K))/OSRT
      XFRN=FPRIM*TFND(L,NY,NX)*(OQN(K,L,NY,NX)*OSRH(KK)
     2-OQN(KK,L,NY,NX)*OSRH(K))/OSRT
      XFRP=FPRIM*TFND(L,NY,NX)*(OQP(K,L,NY,NX)*OSRH(KK)
     2-OQP(KK,L,NY,NX)*OSRH(K))/OSRT
      XFRA=FPRIM*TFND(L,NY,NX)*(OQA(K,L,NY,NX)*OSRH(KK)
     2-OQA(KK,L,NY,NX)*OSRH(K))/OSRT
      IF(ROQCK(K)+XOQCK(K)-XFRK.GT.0.0
     2.AND.ROQCK(KK)+XOQCK(KK)+XFRK.GT.0.0)THEN
      XOQCK(K)=XOQCK(K)-XFRK
      XOQCK(KK)=XOQCK(KK)+XFRK
C     IF(I.EQ.116)THEN
C     WRITE(*,4442)'XOQCK',I,J,NX,NY,L,K,KK,XFRC,ROQCK(K)
C    2,OSRH(K),ROQCK(KK),OSRH(KK),XOQCK(K),XOQCK(KK)
C    3,OQC(K,L,NY,NX),OQC(KK,L,NY,NX)
4442  FORMAT(A8,7I4,12E12.4)
C     ENDIF
      ENDIF
      IF(OQC(K,L,NY,NX)+XOQCZ(K)-XFRC.GT.0.0
     2.AND.OQC(KK,L,NY,NX)+XOQCZ(KK)+XFRC.GT.0.0)THEN
      XOQCZ(K)=XOQCZ(K)-XFRC
      XOQCZ(KK)=XOQCZ(KK)+XFRC
C     IF(I.EQ.116)THEN
C     WRITE(*,4442)'XOQCZ',I,J,NX,NY,L,K,KK,XFRC,OQC(K,L,NY,NX)
C    2,OSRH(K),OQC(KK,L,NY,NX),OSRH(KK),XOQCZ(K),XOQCZ(KK)
C    3,OQC(K,L,NY,NX),OQC(KK,8,NY,NX)
C     ENDIF
      ENDIF
      IF(OQN(K,L,NY,NX)+XOQNZ(K)-XFRN.GT.0.0
     2.AND.OQN(KK,L,NY,NX)+XOQNZ(KK)+XFRN.GT.0.0)THEN
      XOQNZ(K)=XOQNZ(K)-XFRN
      XOQNZ(KK)=XOQNZ(KK)+XFRN
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
C     WRITE(*,4442)'XOQNZ',I,J,NX,NY,L,K,KK,XFRN,OQN(K,L,NY,NX)
C    2,OSRH(K),OQN(KK,L,NY,NX),OSRH(KK),XOQNZ(K),XOQNZ(KK)
C     ENDIF
      ENDIF
      IF(OQP(K,L,NY,NX)+XOQPZ(K)-XFRP.GT.0.0
     2.AND.OQP(KK,L,NY,NX)+XOQPZ(KK)+XFRP.GT.0.0)THEN
      XOQPZ(K)=XOQPZ(K)-XFRP
      XOQPZ(KK)=XOQPZ(KK)+XFRP
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
C     WRITE(*,4442)'XOQPZ',I,J,NX,NY,L,K,KK,XFRP,OQP(K,L,NY,NX)
C    2,OSRH(K),OQP(KK,L,NY,NX),OSRH(KK),XOQPZ(K),XOQPZ(KK)
C     ENDIF
      ENDIF
      IF(OQA(K,L,NY,NX)+XOQAZ(K)-XFRA.GT.0.0
     2.AND.OQA(KK,L,NY,NX)+XOQAZ(KK)+XFRA.GT.0.0)THEN
      XOQAZ(K)=XOQAZ(K)-XFRA
      XOQAZ(KK)=XOQAZ(KK)+XFRA
C     IF((I/1)*1.EQ.I.AND.L.EQ.3.AND.K.EQ.1)THEN
C     WRITE(*,4442)'XOQAZ',I,J,NX,NY,L,K,KK,XFRA,OQA(K,L,NY,NX)
C    2,OSRH(K),OQA(KK,L,NY,NX),OSRH(KK),XOQAZ(K),XOQAZ(KK)
C     ENDIF
      ENDIF
C
C     PRIMING of MICROBIAL C,N,P BETWEEN LITTER AND NON-LITTER C
C
C     XFMC,XFMN,XFMP=transfer of microbial C,N,P
C     between each K and KK, FPRIMM=priming transfer rate constant
C     TFNG=temperature+water effect
C     OMC,OMN,OMP=microbial C,N,P
C     OSRH=total SOC in each K
C     XOMCZ,XOMNZ,XOMPZ=total microbial C,N,P transfer for all K
C
      DO 850 N=1,7
      DO 850 NGL=1,JG
      DO 850 M=1,3
      XFMC=FPRIMM*TFNG(NGL,N,K)*(OMC(M,NGL,N,K,L,NY,NX)*OSRH(KK)
     2-OMC(M,NGL,N,KK,L,NY,NX)*OSRH(K))/OSRT
      XFMN=FPRIMM*TFNG(NGL,N,K)*(OMN(M,NGL,N,K,L,NY,NX)*OSRH(KK)
     2-OMN(M,NGL,N,KK,L,NY,NX)*OSRH(K))/OSRT
      XFMP=FPRIMM*TFNG(NGL,N,K)*(OMP(M,NGL,N,K,L,NY,NX)*OSRH(KK)
     2-OMP(M,NGL,N,KK,L,NY,NX)*OSRH(K))/OSRT
      IF(OMC(M,NGL,N,K,L,NY,NX)+XOMCZ(M,NGL,N,K)-XFMC.GT.0.0
     2.AND.OMC(M,NGL,N,KK,L,NY,NX)+XOMCZ(M,NGL,N,KK)+XFMC.GT.0.0)THEN
      XOMCZ(M,NGL,N,K)=XOMCZ(M,NGL,N,K)-XFMC
      XOMCZ(M,NGL,N,KK)=XOMCZ(M,NGL,N,KK)+XFMC
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,4447)'XOMCZ',I,J,NX,NY,L,K,KK,N,M,XFMC
C    2,OMC(M,NGL,N,K,L,NY,NX)
C    2,OQC(K,L,NY,NX),OMC(M,NGL,N,KK,L,NY,NX),OQC(KK,8,NY,NX)
C    3,XOMCZ(M,NGL,N,K),XOMCZ(M,NGL,N,KK)
4447  FORMAT(A8,9I4,20E12.4)
C     ENDIF
      ENDIF
      IF(OMN(M,NGL,N,K,L,NY,NX)+XOMNZ(M,NGL,N,K)-XFMN.GT.0.0
     2.AND.OMN(M,NGL,N,KK,L,NY,NX)+XOMNZ(M,NGL,N,KK)+XFMN.GT.0.0)THEN
      XOMNZ(M,NGL,N,K)=XOMNZ(M,NGL,N,K)-XFMN
      XOMNZ(M,NGL,N,KK)=XOMNZ(M,NGL,N,KK)+XFMN
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,4447)'XOMNZ',I,J,NX,NY,L,K,KK,N,M,XFMN
C    2,OMN(M,NGL,N,K,L,NY,NX)
C    2,OSRH(K),OMN(M,NGL,N,KK,L,NY,NX),OSRH(KK),XOMNZ(M,NGL,N,K)
C    3,XOMNZ(M,NGL,N,KK)
C     ENDIF
      ENDIF
      IF(OMP(M,NGL,N,K,L,NY,NX)+XOMPZ(M,NGL,N,K)-XFMP.GT.0.0
     2.AND.OMP(M,NGL,N,KK,L,NY,NX)+XOMPZ(M,NGL,N,KK)+XFMP.GT.0.0)THEN
      XOMPZ(M,NGL,N,K)=XOMPZ(M,NGL,N,K)-XFMP
      XOMPZ(M,NGL,N,KK)=XOMPZ(M,NGL,N,KK)+XFMP
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,4447)'XOMPZ',I,J,NX,NY,L,K,KK,N,M,XFMP
C    2,OMP(M,NGL,N,K,L,NY,NX),OSRH(K),OMP(M,NGL,N,KK,L,NY,NX),OSRH(KK)
C    2,XOMPZ(M,NGL,N,K),XOMPZ(M,NGL,N,KK)
C     ENDIF
      ENDIF
850   CONTINUE
      ENDIF
800   CONTINUE
      ENDIF
795   CONTINUE
C
C     TRANSFER ALL PRIMING AMONG ALL K
C
C     TOQCK=total respiration of DOC+DOA in soil layer
C     ROQCK=total respiration of DOC+DOA in substrate complex
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
C     OMC,OMN,OMP=microbial C,N,P
C
      TOQCK(L,NY,NX)=0.0
      DO 1790 K=0,KL
      ROQCK(K)=ROQCK(K)+XOQCK(K)
      TOQCK(L,NY,NX)=TOQCK(L,NY,NX)+ROQCK(K)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+XOQCZ(K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+XOQNZ(K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+XOQPZ(K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+XOQAZ(K)
      DO 840 N=1,7
      DO 840 NGL=1,JG
      DO 840 M=1,3
      OMC(M,NGL,N,K,L,NY,NX)=OMC(M,NGL,N,K,L,NY,NX)+XOMCZ(M,NGL,N,K)
      OMN(M,NGL,N,K,L,NY,NX)=OMN(M,NGL,N,K,L,NY,NX)+XOMNZ(M,NGL,N,K)
      OMP(M,NGL,N,K,L,NY,NX)=OMP(M,NGL,N,K,L,NY,NX)+XOMPZ(M,NGL,N,K)
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,5559)'XOM',I,J,NX,NY,L,K,N,M,OMC(M,NGL,N,K,L,NY,NX)
C    2,OMN(M,NGL,N,K,L,NY,NX),OMP(M,NGL,N,K,L,NY,NX)
C    3,XOMCZ(M,NGL,N,K),XOMNZ(M,NGL,N,K),XOMPZ(M,NGL,N,K)
5559  FORMAT(A8,8I4,12E12.4)
C     ENDIF
840   CONTINUE

C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
C     WRITE(*,4443)'PRIM2',I,J,NX,NY,L,K,ROQCK(K)
C    2,XOQCK(K),OQC(K,L,NY,NX),XOQCZ(K),OQN(K,L,NY,NX),XOQNZ(K)
C    3,OQP(K,L,NY,NX),XOQPZ(K),OQA(K,L,NY,NX),XOQAZ(K),TOMK(K)
C    3,TONK(K),TOPK(K),TONX(K),TOPX(K),CNOMX,CPOMX,FCNK(K),FCPK(K)
C    4,TOQCK(L,NY,NX)
4443  FORMAT(A8,6I4,20E12.4)
C     ENDIF
C
C     DECOMPOSITION OF ORGANIC SUBSTRATES
C
C     FCPK=N,P limitation to microbial activity in each K
C     CNOMX,CPOMX=N:C,P:C ratios relative to set maximum values
C     COQCK=aqueous concentration of microbial activity
C     DCKD=Km for decomposition of SOC at current COQCK
C     DCKM0,DCKML=Km for decomposition of SOC at zero COQCK
C     DCKI=inhibition of decomposition by microbial concentration
C     OSRH=total SOC
C     COSC=concentration of total SOC
C     BKVL,VOLX=mass, volume of soil layer
C     DFNS=effect of microbial concentration on decomposition
C     OQCI=DOC product inhibition for decomposition
C     OQKI=DOC product inhibition constant for decomposition
C
      IF(TOMK(K).GT.ZEROS(NY,NX))THEN
      CNOMX=TONK(K)/TONX(K)
      CPOMX=TOPK(K)/TOPX(K)
      FCNK(K)=AMIN1(1.0,AMAX1(0.50,CNOMX))
      FCPK(K)=AMIN1(1.0,AMAX1(0.50,CPOMX))
      ELSE
      FCNK(K)=1.0
      FCPK(K)=1.0
      ENDIF
C
C     AQUEOUS CONCENTRATION OF BIOMASS TO CACULATE INHIBITION
C     CONSTANT FOR DECOMPOSITION
C
      IF(VOLWZ.GT.ZEROS2(NY,NX))THEN
      COQCK=AMIN1(0.1E+06,ROQCK(K)/VOLWZ)
      ELSE
      COQCK=0.1E+06
      ENDIF
      IF(L.EQ.0)THEN
      DCKD=DCKM0*(1.0+COQCK/DCKI)
      ELSE
      DCKD=DCKML*(1.0+COQCK/DCKI)
      ENDIF
      IF(OSRH(K).GT.ZEROS(NY,NX))THEN
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      COSC=OSRH(K)/BKVL(L,NY,NX)
      ELSE
      COSC=OSRH(K)/VOLY(L,NY,NX)
      ENDIF
      DFNS=COSC/(COSC+DCKD)
      OQCI=1.0/(1.0+COQC(K,L,NY,NX)/OQKI)
C     IF(L.EQ.0.AND.J.EQ.15)THEN
C     WRITE(*,4242)'COSC',I,J,L,K,DFNS,COSC,COQCK,DCKD,OSRH(K)
C    2,OSAT(K),OSCT(K),ORCT(K),OHC(K,L,NY,NX),BKVL(L,NY,NX),ROQCK(K)
C    3,VOLWZ,VOLWRX(NY,NX),VOLW(0,NY,NX),FCR(NY,NX)
C    4,THETY(L,NY,NX)
4242  FORMAT(A8,4I4,30E12.4)
C     ENDIF
C
C     C, N, P DECOMPOSITION RATE OF SOLID SUBSTRATES 'RDOS*' FROM
C     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
C     TEMPERATURE, SUBSTRATE C:N, C:P
C
C     CNS,CPS=N:C,P:C ratios of SOC
C     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
C     OSA,OSN,OSP=active biomass C,N,P
C     SPOSC=specific decomposition rate constant
C     ROQCK=total respiration of DOC+DOA used to represent microbial activity
C     DFNS=effect of microbial concentration on decomposition
C     OQCI=DOC product inhibition for decomposition
C     TFNX=temperature stress effect
C     OSRH=total SOC
C     FCNK,FCPK=N,P limitation to microbial activity in each K
C
      DO 785 M=1,4
      IF(OSC(M,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNS(M,K)=AMAX1(0.0,OSN(M,K,L,NY,NX)/OSC(M,K,L,NY,NX))
      CPS(M,K)=AMAX1(0.0,OSP(M,K,L,NY,NX)/OSC(M,K,L,NY,NX))
      RDOSC(M,K)=AMAX1(0.0,AMIN1(0.5*OSA(M,K,L,NY,NX)
     2,SPOSC(M,K)*ROQCK(K)*DFNS*OQCI*TFNX*OSA(M,K,L,NY,NX)/OSRH(K)))
C    3*AMIN1(FCNK(K),FCPK(K))
      RDOSN(M,K)=AMAX1(0.0,AMIN1(OSN(M,K,L,NY,NX)
     2,CNS(M,K)*RDOSC(M,K)))/FCNK(K)
      RDOSP(M,K)=AMAX1(0.0,AMIN1(OSP(M,K,L,NY,NX)
     2,CPS(M,K)*RDOSC(M,K)))/FCPK(K)
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
C     WRITE(*,4444)'RDOSC',I,J,NX,NY,L,K,M,RDOSC(M,K),RDOSN(M,K)
C    2,RDOSP(M,K),CNS(M,K),CPS(M,K),SPOSC(M,K),ROQCK(K),DFNS,TFNX
C    3,OQCI,OSA(M,K,L,NY,NX),OSRH(K),COSC,COQCK,DCKD,VOLWZ
C    4,COXYS(L,NY,NX),TKS(L,NY,NX),PSISM(L,NY,NX),THETW(L,NY,NX)
C    5,WFN(1,K),WFN(3,K),OXYI,VOLW(0,NY,NX),VOLWRX(NY,NX)
C    4,FOSRH(K,L,NY,NX),VOLY(L,NY,NX),ORGC(L,NY,NX),OSC(M,K,L,NY,NX)
C    2,OSN(M,K,L,NY,NX),OSP(M,K,L,NY,NX),TONK(K),TONX(K),FCNK(K)
C    6,FCPK(K),WFN(1,K),WFN(3,K),WFN(4,K),COQC(K,L,NY,NX)
C    7,OQC(K,L,NY,NX),VOLWM(NPH,L,NY,NX)
4444  FORMAT(A8,7I4,50E12.4)
C     ENDIF
      ELSE
      CNS(M,K)=CNOSC(M,K,L,NY,NX)
      CPS(M,K)=CPOSC(M,K,L,NY,NX)
      RDOSC(M,K)=0.0
      RDOSN(M,K)=0.0
      RDOSP(M,K)=0.0
      ENDIF
785   CONTINUE
C
C     HUMIFICATION OF DECOMPOSED RESIDUE LIGNIN WITH PROTEIN,
C     CH2O AND CELLULOSE 'RHOS*' WITH REMAINDER 'RCOS*' TO DOC,DON,DOP
C
C     RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
C     RDOSC,RDOSN,RDOSP=decomposition of SOC,SON,SOP
C     CNRH,CPRH=N:C,P:C in POC
C     EPOC=fraction of RDOSC allocated to POC from hour1.f
C     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
C
      IF(K.LE.2)THEN
      RHOSC(4,K)=AMAX1(0.0,AMIN1(RDOSN(4,K)/CNRH(3)
     2,RDOSP(4,K)/CPRH(3),EPOC(L,NY,NX)*RDOSC(4,K)))
      RHOSCM=0.10*RHOSC(4,K)
      RHOSC(1,K)=AMAX1(0.0,AMIN1(RDOSC(1,K),RDOSN(1,K)/CNRH(3)
     2,RDOSP(1,K)/CPRH(3),RHOSCM))
      RHOSC(2,K)=AMAX1(0.0,AMIN1(RDOSC(2,K),RDOSN(2,K)/CNRH(3)
     2,RDOSP(2,K)/CPRH(3),RHOSCM))
      RHOSC(3,K)=AMAX1(0.0,AMIN1(RDOSC(3,K),RDOSN(3,K)/CNRH(3)
     2,RDOSP(3,K)/CPRH(3),RHOSCM-RHOSC(2,K)))
      DO 805 M=1,4
      RHOSN(M,K)=AMIN1(RDOSN(M,K),RHOSC(M,K)*CNRH(3))
      RHOSP(M,K)=AMIN1(RDOSP(M,K),RHOSC(M,K)*CPRH(3))
      RCOSC(M,K)=RDOSC(M,K)-RHOSC(M,K)
      RCOSN(M,K)=RDOSN(M,K)-RHOSN(M,K)
      RCOSP(M,K)=RDOSP(M,K)-RHOSP(M,K)
805   CONTINUE
      ELSE
      DO 810 M=1,4
      RHOSC(M,K)=0.0
      RHOSN(M,K)=0.0
      RHOSP(M,K)=0.0
      RCOSC(M,K)=RDOSC(M,K)
      RCOSN(M,K)=RDOSN(M,K)
      RCOSP(M,K)=RDOSP(M,K)
810   CONTINUE
      ENDIF
      ELSE
      DO 780 M=1,4
      RDOSC(M,K)=0.0
      RDOSN(M,K)=0.0
      RDOSP(M,K)=0.0
      RHOSC(M,K)=0.0
      RHOSN(M,K)=0.0
      RHOSP(M,K)=0.0
      RCOSC(M,K)=0.0
      RCOSN(M,K)=0.0
      RCOSP(M,K)=0.0
780   CONTINUE
      ENDIF
C
C     C, N, P DECOMPOSITION RATE OF BIORESIDUE 'RDOR*' FROM
C     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
C     TEMPERATURE, SUBSTRATE C:N, C:P
C
C     ORC,ORN,ORP=microbial residue C,N,P
C     CNR,CPR=N:C,P:C ratios of microbial residue
C     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
C     SPORC=specific decomposition rate constant for microbial residue
C     ROQCK=total respiration of DOC+DOA used to represent microbial activity
C     DFNS=effect of microbial concentration on decomposition
C     OQCI=DOC product inhibition for decomposition
C     TFNX=temperature stress effect
C     OSRH=total SOC
C     FCNK,FCPK=N,P limitation to microbial activity in each K
C
      IF(OSRH(K).GT.ZEROS(NY,NX))THEN
      DO 775 M=1,2
      IF(ORC(M,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNR=AMAX1(0.0,ORN(M,K,L,NY,NX)/ORC(M,K,L,NY,NX))
      CPR=AMAX1(0.0,ORP(M,K,L,NY,NX)/ORC(M,K,L,NY,NX))
      RDORC(M,K)=AMAX1(0.0,AMIN1(ORC(M,K,L,NY,NX)
     2,SPORC(M)*ROQCK(K)*DFNS*OQCI*TFNX*ORC(M,K,L,NY,NX)/OSRH(K)))
C    3*AMIN1(FCNK(K),FCPK(K))
      RDORN(M,K)=AMAX1(0.0,AMIN1(ORN(M,K,L,NY,NX),CNR*RDORC(M,K)))
     2/FCNK(K)
      RDORP(M,K)=AMAX1(0.0,AMIN1(ORP(M,K,L,NY,NX),CPR*RDORC(M,K)))
     2/FCPK(K)
      ELSE
      RDORC(M,K)=0.0
      RDORN(M,K)=0.0
      RDORP(M,K)=0.0
      ENDIF
775   CONTINUE
      ELSE
      DO 776 M=1,2
      RDORC(M,K)=0.0
      RDORN(M,K)=0.0
      RDORP(M,K)=0.0
776   CONTINUE
      ENDIF
C
C     C, N, P DECOMPOSITION RATE OF SORBED SUBSTRATES 'RDOH*' FROM
C     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
C     TEMPERATURE, SUBSTRATE C:N, C:P
C
C     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
C     CNH,CPH=N:C,P:C ratios of adsorbed C,N,P
C     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
C     SPOHC=specific decomposition rate constant for adsorbed C
C     ROQCK=total respiration of DOC+DOA used to represent microbial activity
C     DFNS=effect of microbial concentration on decomposition
C     OQCI=DOC product inhibition for decomposition
C     TFNX=temperature stress effect
C     OSRH=total SOC
C     FCNK,FCPK=N,P limitation to microbial activity in each K
C
      IF(OSRH(K).GT.ZEROS(NY,NX))THEN
      IF(OHC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNH(K)=AMAX1(0.0,OHN(K,L,NY,NX)/OHC(K,L,NY,NX))
      CPH(K)=AMAX1(0.0,OHP(K,L,NY,NX)/OHC(K,L,NY,NX))
      RDOHC(K)=AMAX1(0.0,AMIN1(OHC(K,L,NY,NX)
     2,SPOHC*ROQCK(K)*DFNS*OQCI*TFNX*OHC(K,L,NY,NX)/OSRH(K)))
C    3*AMIN1(FCNK(K),FCPK(K))
      RDOHN(K)=AMAX1(0.0,AMIN1(OHN(K,L,NY,NX),CNH(K)*RDOHC(K)))
     2/FCNK(K)
      RDOHP(K)=AMAX1(0.0,AMIN1(OHP(K,L,NY,NX),CPH(K)*RDOHC(K)))
     2/FCPK(K)
      RDOHA(K)=AMAX1(0.0,AMIN1(OHA(K,L,NY,NX)
     2,SPOHA*ROQCK(K)*DFNS*TFNX*OHA(K,L,NY,NX)/OSRH(K)))
C    3*AMIN1(FCNK(K),FCPK(K))
      ELSE
      CNH(K)=0.0
      CPH(K)=0.0
      RDOHC(K)=0.0
      RDOHN(K)=0.0
      RDOHP(K)=0.0
      RDOHA(K)=0.0
      ENDIF
      ELSE
      CNH(K)=0.0
      CPH(K)=0.0
      RDOHC(K)=0.0
      RDOHN(K)=0.0
      RDOHP(K)=0.0
      RDOHA(K)=0.0
      ENDIF
C
C     DOC ADSORPTION - DESORPTION
C
C     VOLWM=soil water content, FOSRH=fraction of total SOC
C     AEC,AECX=anion exchange capacity
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
C     TCGOQC,TCGOMN,TCGOMP,TCGOAC=total uptake of DOC,DON,DOP,acetate
C     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
C     TSORP,HSORP=sorption rate constant and coefficient for OHC
C     FOCA,FOAA=fractions of DOC and acetate vs. DOC+acetate
C     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
C
      IF(VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX)
     2.AND.FOSRH(K,L,NY,NX).GT.ZERO)THEN
      IF(L.EQ.0)THEN
      AECX=0.5E+03
      ELSE
      AECX=AEC(L,NY,NX)
      ENDIF
      OQCX=AMAX1(ZEROS(NY,NX),OQC(K,L,NY,NX)-TCGOQC(K))
      OQNX=AMAX1(ZEROS(NY,NX),OQN(K,L,NY,NX)-TCGOMN(K))
      OQPX=AMAX1(ZEROS(NY,NX),OQP(K,L,NY,NX)-TCGOMP(K))
      OQAX=AMAX1(ZEROS(NY,NX),OQA(K,L,NY,NX)-TCGOAC(K))
      OHCX=AMAX1(ZEROS(NY,NX),OHC(K,L,NY,NX))
      OHNX=AMAX1(ZEROS(NY,NX),OHN(K,L,NY,NX))
      OHPX=AMAX1(ZEROS(NY,NX),OHP(K,L,NY,NX))
      OHAX=AMAX1(ZEROS(NY,NX),OHA(K,L,NY,NX))
      VOLXX=BKVL(L,NY,NX)*AECX*HSORP*FOSRH(K,L,NY,NX)
      VOLXW=VOLWM(NPH,L,NY,NX)*FOSRH(K,L,NY,NX)
      IF(FOCA(K).GT.ZERO)THEN
      VOLCX=FOCA(K)*VOLXX
      VOLCW=FOCA(K)*VOLXW
      CSORP(K)=TSORP*(OQCX*VOLCX-OHCX*VOLCW)/(VOLCX+VOLCW)
      ELSE
      CSORP(K)=TSORP*(OQCX*VOLXX-OHCX*VOLXW)/(VOLXX+VOLXW)
      ENDIF
      IF(FOAA(K).GT.ZERO)THEN
      VOLAX=FOAA(K)*VOLXX
      VOLAW=FOAA(K)*VOLXW
      CSORPA(K)=TSORP*(OQAX*VOLAX-OHAX*VOLAW)/(VOLAX+VOLAW)
      ELSE
      CSORPA(K)=TSORP*(OQAX*VOLXX-OHAX*VOLXW)/(VOLXX+VOLXW)
      ENDIF
      ZSORP(K)=TSORP*(OQNX*VOLXX-OHNX*VOLXW)/(VOLXX+VOLXW)
      PSORP(K)=TSORP*(OQPX*VOLXX-OHPX*VOLXW)/(VOLXX+VOLXW)
      ELSE
      CSORP(K)=0.0
      CSORPA(K)=0.0
      ZSORP(K)=0.0
      PSORP(K)=0.0
      ENDIF
C     IF(I.EQ.116)THEN
C     WRITE(*,591)'CSORP',I,J,NX,NY,L,K,CSORP(K),CSORPA(K)
C    1,OQC(K,L,NY,NX),OHC(K,L,NY,NX),OQA(K,L,NY,NX),OHA(K,L,NY,NX)
C    2,OQC(K,L,NY,NX)/VOLWM(NPH,L,NY,NX)
C    2,OQA(K,L,NY,NX)/VOLWM(NPH,L,NY,NX)
C    3,OHC(K,L,NY,NX)/BKVL(L,NY,NX),OHA(K,L,NY,NX)/BKVL(L,NY,NX)
C    4,BKVL(L,NY,NX),VOLWM(NPH,L,NY,NX),FOCA(K),FOAA(K)
C    5,FOSRH(K,L,NY,NX),TCGOQC(K),OQCX
591   FORMAT(A8,6I4,40E12.4)
C     ENDIF
1790  CONTINUE
C
C     REDISTRIBUTE AUTOTROPHIC DECOMPOSITION PRODUCTS AMONG
C     HETEROTROPHIC SUBSTRATE-MICROBE COMPLEXES
C
C     FORC=fraction of total microbial residue
C     ORCT=microbial residue
C     RCCMC,RCCMN,RCCMP=transfer of auto litterfall C,N,P to each hetero K
C     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
C     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
C
      DO 1690 K=0,KL
      IF(TORC.GT.ZEROS(NY,NX))THEN
      FORC(K)=ORCT(K)/TORC
      ELSE
      IF(K.EQ.3)THEN
      FORC(K)=1.0
      ELSE
      FORC(K)=0.0
      ENDIF
      ENDIF
      DO 1685 N=1,7
      DO 1685 NGL=1,JG
      DO 1680 M=1,2
      RCCMC(M,NGL,N,K)=(RCOMC(M,NGL,N,5)+RCMMC(M,NGL,N,5))*FORC(K)
      RCCMN(M,NGL,N,K)=(RCOMN(M,NGL,N,5)+RCMMN(M,NGL,N,5))*FORC(K)
      RCCMP(M,NGL,N,K)=(RCOMP(M,NGL,N,5)+RCMMP(M,NGL,N,5))*FORC(K)
C     IF(L.EQ.0)THEN
C     WRITE(*,8821)'RCCMC',I,J,L,K,N,M,RCCMC(M,NGL,N,K)
C    2,RCOMC(M,NGL,N,5),RCMMC(M,NGL,N,5),FORC(K)
C     ENDIF
1680  CONTINUE
1685  CONTINUE
1690  CONTINUE
C
C     REDISTRIBUTE C,N AND P TRANSFORMATIONS AMONG STATE
C     VARIABLES IN SUBSTRATE-MICROBE COMPLEXES
C
      DO 590 K=0,KL
      DO 580 M=1,4
C
C     SUBSTRATE DECOMPOSITION PRODUCTS
C
C     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
C     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP
C     OQC,OQN,OQP,OQA=DOC,DON,DOP
C     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
C
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-RDOSC(M,K)
C     OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-RDOSC(M,K)
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-RDOSN(M,K)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-RDOSP(M,K)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RCOSC(M,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RCOSN(M,K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RCOSP(M,K)
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.K.EQ.4)THEN
C     WRITE(*,4444)'RDOSC',I,J,NX,NY,L,K,M,OSC(M,K,L,NY,NX)
C    2,RDOSC(M,K)
C     ENDIF
C
C     LIGNIFICATION PRODUCTS
C
C     RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
C
      IF(L.NE.0)THEN
      OSC(1,3,L,NY,NX)=OSC(1,3,L,NY,NX)+RHOSC(M,K)
C     OSA(1,3,L,NY,NX)=OSA(1,3,L,NY,NX)+RHOSC(M,K)
      OSN(1,3,L,NY,NX)=OSN(1,3,L,NY,NX)+RHOSN(M,K)
      OSP(1,3,L,NY,NX)=OSP(1,3,L,NY,NX)+RHOSP(M,K)
      ELSE
      OSC(1,3,NU(NY,NX),NY,NX)=OSC(1,3,NU(NY,NX),NY,NX)+RHOSC(M,K)
C     OSA(1,3,NU(NY,NX),NY,NX)=OSA(1,3,NU(NY,NX),NY,NX)+RHOSC(M,K)
      OSN(1,3,NU(NY,NX),NY,NX)=OSN(1,3,NU(NY,NX),NY,NX)+RHOSN(M,K)
      OSP(1,3,NU(NY,NX),NY,NX)=OSP(1,3,NU(NY,NX),NY,NX)+RHOSP(M,K)
      ENDIF
580   CONTINUE
C
C     MICROBIAL RESIDUE DECOMPOSITION PRODUCTS
C
C     ORC,ORN,ORP=microbial residue C,N,P
C     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
C     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
C     RCOQN=DON production from nitrous acid reduction
C
      DO 575 M=1,2
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-RDORC(M,K)
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-RDORN(M,K)
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-RDORP(M,K)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RDORC(M,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RDORN(M,K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RDORP(M,K)
575   CONTINUE
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RDOHC(K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RDOHN(K)+RCOQN*FORC(K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RDOHP(K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+RDOHA(K)
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-RDOHC(K)
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-RDOHN(K)
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-RDOHP(K)
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-RDOHA(K)
C
C     MICROBIAL UPTAKE OF DISSOLVED C, N, P
C
C     CGOQC,CGOAC,CGOMN,CGOMP=DOC,acetate,DON,DOP uptake
C     RCH3X=acetate production from fermentation
C
      DO 570 N=1,7
      DO 570 NGL=1,JG
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-CGOQC(NGL,N,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-CGOMN(NGL,N,K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-CGOMP(NGL,N,K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-CGOAC(NGL,N,K)+RCH3X(NGL,N,K)
C
C     MICROBIAL DECOMPOSITION PRODUCTS
C
C     ORC,ORN,ORP=microbial residue C,N,P
C     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to residue
C     RCCMC,RCCMN,RCCMP=transfer of auto litterfall C,N,P to each hetero K
C     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to residue
C
      DO 565 M=1,2
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)+RCOMC(M,NGL,N,K)
     2+RCCMC(M,NGL,N,K)+RCMMC(M,NGL,N,K)
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)+RCOMN(M,NGL,N,K)
     2+RCCMN(M,NGL,N,K)+RCMMN(M,NGL,N,K)
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)+RCOMP(M,NGL,N,K)
     2+RCCMP(M,NGL,N,K)+RCMMP(M,NGL,N,K)
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4.AND.K.EQ.2)THEN
C     WRITE(*,8821)'ORC',I,J,L,K,N,M,ORC(M,K,L,NY,NX)
C    2,RCOMC(M,NGL,N,K),RCCMC(M,NGL,N,K),RCMMC(M,NGL,N,K),RDORC(M,K)
C     WRITE(*,8821)'ORP',I,J,L,K,N,M,ORP(M,K,L,NY,NX)
C    2,RCOMP(M,NGL,N,K),RCCMP(M,NGL,N,K),RCMMP(M,NGL,N,K),RDORP(M,K)
8821  FORMAT(A8,6I4,40E12.4)
C     ENDIF
565   CONTINUE
570   CONTINUE
C
C     SORPTION PRODUCTS
C
C     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of OQC,acetate,DON,DOP
C
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-CSORP(K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-ZSORP(K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-PSORP(K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-CSORPA(K)
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)+CSORP(K)
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)+ZSORP(K)
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)+PSORP(K)
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)+CSORPA(K)
C     IF(L.EQ.1)THEN
C     WRITE(*,592)'OQC',I,J,NX,NY,L,K,OQC(K,L,NY,NX)
C    2,OQA(K,L,NY,NX),(RCOSC(M,K),M=1,4),(RDORC(M,K),M=1,2)
C    3,RDOHC(K),(CGOQC(NGL,N,K),NGL=1,JG,N=1,7),CSORP(K),OHC(K,L,NY,NX)
C    4,(WFN(NGL,N,K),NGL=1,JG,N=1,7),RDOHA(K),(RCH3X(NGL,N,K),NGL=1,JG,N=1,7)
C    3,(CGOAC(NGL,N,K),NGL=1,JG,N=1,7),CSORPA(K),OHA(K,L,NY,NX)
C     WRITE(*,592)'OQN',I,J,NX,NY,L,K,OQN(K,L,NY,NX)
C    2,(RCOSN(M,K),M=1,4),(RDORN(M,K),M=1,2),RDOHN(K)
C    2,RCOQN*FORC(K),(CGOMN(NGL,N,K),N=1,NGL,N=1,7),ZSORP(K),OHN(K,L,NY,NX)
592   FORMAT(A8,6I4,80E12.4)
C     ENDIF
590   CONTINUE
C
C     MICROBIAL GROWTH FROM RESPIRATION, MINERALIZATION
C
C     OMC,OMN,OMP=microbial C,N,P
C     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
C     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition
C     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence
C
      DO 550 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 545 N=1,7
      DO 545 NGL=1,JG
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      DO 540 M=1,2
      OMC(M,NGL,N,K,L,NY,NX)=OMC(M,NGL,N,K,L,NY,NX)+CGOMS(M,NGL,N,K)
     2-RXOMC(M,NGL,N,K)-RXMMC(M,NGL,N,K)
      OMN(M,NGL,N,K,L,NY,NX)=OMN(M,NGL,N,K,L,NY,NX)+CGONS(M,NGL,N,K)
     2-RXOMN(M,NGL,N,K)-RXMMN(M,NGL,N,K)
      OMP(M,NGL,N,K,L,NY,NX)=OMP(M,NGL,N,K,L,NY,NX)+CGOPS(M,NGL,N,K)
     2-RXOMP(M,NGL,N,K)-RXMMP(M,NGL,N,K)


C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,4488)'OMN2',I,J,NX,NY,L,K,N,M,OMC(M,NGL,N,K,L,NY,NX)
C    2,OMN(M,NGL,N,K,L,NY,NX),OMP(M,NGL,N,K,L,NY,NX)
C     WRITE(*,4488)'RDOMC',I,J,NX,NY,L,K,N,M,CGOMS(M,NGL,N,K)
C    2,CGOQC(NGL,N,K)
C    4,CGOAC(NGL,N,K),RGOMO(NGL,N,K),RGOMD(NGL,N,K),RXOMC(M,NGL,N,K)
C    5,RXMMC(M,NGL,N,K)
C    3,RMOMC(M,NGL,N,K),TFNX,OMGR,OMC(3,NGL,N,K,L,NY,NX),WFN(NGL,N,K)
C    3,OMC(M,NGL,N,K,L,NY,NX),OMA(NGL,N,K),TSRH
C    4,RCH3X(NGL,N,K),RH2GZ,RH2GX(4,K),FOCA(K),FOAA(K)
C    6,OQA(K,L,NY,NX),OHA(K,L,NY,NX),OQC(K,L,NY,NX),OHC(K,L,NY,NX)
C    7,OMP(M,NGL,N,K,L,NY,NX),CGOPS(M,NGL,N,K),RDOMP(M,NGL,N,K)
C    8,RDMMP(M,NGL,N,K)
C    8,OMP(3,NGL,N,K,L,NY,NX),CGOMP(NGL,N,K),RIPO4(NGL,N,K)
4488  FORMAT(A8,8I4,40E12.4)
C     ENDIF
C
C     HUMIFICATION PRODUCTS
C
C     CFOMC=fractions allocated to humic vs fulvic humus
C     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall to humus
C     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P to humus
C
      IF(L.NE.0)THEN
      OSC(1,4,L,NY,NX)=OSC(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
     2*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
C     OSA(1,4,L,NY,NX)=OSA(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
C    2*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
      OSN(1,4,L,NY,NX)=OSN(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
     2*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
      OSP(1,4,L,NY,NX)=OSP(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
     2*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
      OSC(2,4,L,NY,NX)=OSC(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
     2*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
C     OSA(2,4,L,NY,NX)=OSA(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
C    2*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
      OSN(2,4,L,NY,NX)=OSN(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
     2*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
      OSP(2,4,L,NY,NX)=OSP(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
     2*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
C     IF((I/10)*10.EQ.I.AND.J.EQ.24)THEN
C     WRITE(*,4445)'RHOMC',I,J,NX,NY,L,K,M,N,OSC(1,4,L,NY,NX)
C    2,OSC(2,4,L,NY,NX),CFOMC(1,L,NY,NX),CFOMC(2,L,NY,NX)
C    3,RHOMC(M,NGL,N,K),RHMMC(M,NGL,N,K)
4445  FORMAT(A8,8I4,40E12.4)
C     ENDIF
      ELSE
      OSC(1,4,NU(NY,NX),NY,NX)=OSC(1,4,NU(NY,NX),NY,NX)
     2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
C     OSA(1,4,NU(NY,NX),NY,NX)=OSA(1,4,NU(NY,NX),NY,NX)
C    2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
      OSN(1,4,NU(NY,NX),NY,NX)=OSN(1,4,NU(NY,NX),NY,NX)
     2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
      OSP(1,4,NU(NY,NX),NY,NX)=OSP(1,4,NU(NY,NX),NY,NX)
     2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
      OSC(2,4,NU(NY,NX),NY,NX)=OSC(2,4,NU(NY,NX),NY,NX)
     2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
C     OSA(2,4,NU(NY,NX),NY,NX)=OSA(2,4,NU(NY,NX),NY,NX)
C    2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMC(M,NGL,N,K)+RHMMC(M,NGL,N,K))
      OSN(2,4,NU(NY,NX),NY,NX)=OSN(2,4,NU(NY,NX),NY,NX)
     2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMN(M,NGL,N,K)+RHMMN(M,NGL,N,K))
      OSP(2,4,NU(NY,NX),NY,NX)=OSP(2,4,NU(NY,NX),NY,NX)
     2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMP(M,NGL,N,K)+RHMMP(M,NGL,N,K))
      ENDIF
540   CONTINUE
C
C     INPUTS TO NONSTRUCTURAL POOLS
C
C     CGOMC=total DOC+acetate uptake
C     RGOMO=total respiration
C     RGOMD=respiration for denitrifcation
C     RGN2F=respiration for N2 fixation
C     RCO2X=total CO2 emission
C     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural C,N,P
C     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling
C     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
C     CGOMN,CGOMP=DON, DOP uptake
C     RINH4,RINB4=substrate-limited NH4 mineraln-immobiln in non-band, band
C     RINO3,RINB3=substrate-limited NO3 immobiln in non-band, band
C     RIPO4,RIPOB=substrate-limited H2PO4 mineraln-immobn in non-band, band
C     RIP14,RIP1B=substrate-limited HPO4 mineraln-immobn in non-band, band
C     RINH4R,RINO3R =substrate-limited NH4,NO3 mineraln-immobiln
C     RIPO4R,RIP14R=substrate-limited H2PO4,HPO4 mineraln-immobiln
C
      CGROMC=CGOMC(NGL,N,K)-RGOMO(NGL,N,K)-RGOMD(NGL,N,K)-RGN2F(NGL,N,K)
      RCO2X(NGL,N,K)=RCO2X(NGL,N,K)+RGN2F(NGL,N,K)
      DO 555 M=1,2
      OMC(3,NGL,N,K,L,NY,NX)=OMC(3,NGL,N,K,L,NY,NX)-CGOMS(M,NGL,N,K)
     2+R3OMC(M,NGL,N,K)
      OMN(3,NGL,N,K,L,NY,NX)=OMN(3,NGL,N,K,L,NY,NX)-CGONS(M,NGL,N,K)
     2+R3OMN(M,NGL,N,K)+R3MMN(M,NGL,N,K)
      OMP(3,NGL,N,K,L,NY,NX)=OMP(3,NGL,N,K,L,NY,NX)-CGOPS(M,NGL,N,K)
     2+R3OMP(M,NGL,N,K)+R3MMP(M,NGL,N,K)
      RCO2X(NGL,N,K)=RCO2X(NGL,N,K)+R3MMC(M,NGL,N,K)
555   CONTINUE
      OMC(3,NGL,N,K,L,NY,NX)=OMC(3,NGL,N,K,L,NY,NX)+CGROMC
      OMN(3,NGL,N,K,L,NY,NX)=OMN(3,NGL,N,K,L,NY,NX)+CGOMN(NGL,N,K)
     2+RINH4(NGL,N,K)+RINB4(NGL,N,K)+RINO3(NGL,N,K)+RINB3(NGL,N,K)
     2+RN2FX(NGL,N,K)
      OMP(3,NGL,N,K,L,NY,NX)=OMP(3,NGL,N,K,L,NY,NX)+CGOMP(NGL,N,K)
     2+RIPO4(NGL,N,K)+RIPOB(NGL,N,K)+RIP14(NGL,N,K)+RIP1B(NGL,N,K)
      IF(L.EQ.0)THEN
      OMN(3,NGL,N,K,L,NY,NX)=OMN(3,NGL,N,K,L,NY,NX)+RINH4R(NGL,N,K)
     2+RINO3R(NGL,N,K)
      OMP(3,NGL,N,K,L,NY,NX)=OMP(3,NGL,N,K,L,NY,NX)+RIPO4R(NGL,N,K)
     2+RIP14R(NGL,N,K)
      ENDIF
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,5556)'OMC3',I,J,NX,NY,L,K,N,OMC(3,NGL,N,K,L,NY,NX)
C    2,OMN(3,NGL,N,K,L,NY,NX),OMP(3,NGL,N,K,L,NY,NX),OMC(1,NGL,N,K,L,NY,NX)
C    2,OMN(1,NGL,N,K,L,NY,NX),OMP(1,NGL,N,K,L,NY,NX),WFN(NGL,N,K),OXYI
C    2,COXYS(L,NY,NX)
C    2,CGOMS(1,NGL,N,K),CGOMS(2,NGL,N,K),CGROMC
C    3,CGOPS(1,NGL,N,K),CGOPS(2,NGL,N,K),CGOMP(NGL,N,K),RIPO4(NGL,N,K)
C    4,CGOMC(NGL,N,K),RGOMO(NGL,N,K),RGOMD(NGL,N,K),RMOMT,WFN(NGL,N,K)
C    5,(CGONS(M,NGL,N,K),M=1,2),(R3OMN(M,NGL,N,K),M=1,2),(R3MMN(M,NGL,N,K)
C    6,M=1,2),(XOMCZ(M,NGL,N,K),M=1,2)
C    6,CGOMN(NGL,N,K),RINH4(NGL,N,K),RINB4(NGL,N,K),RINO3(NGL,N,K),RINB3(NGL,N,K)
C    7,RN2FX(NGL,N,K)
5556  FORMAT(A8,7I4,60E12.4)
C     ENDIF
      ENDIF
545   CONTINUE
      ENDIF
550   CONTINUE

C
C     MICROBIAL COLONIZATION OF NEW LITTER
C
C     OSCT,OSAT,OSCX=total,colonized,uncolonized SOC
C     OSA,OSC=colonized,total litter
C     DOSA=rate constant for litter colonization
C     ROQCK=total respiration of DOC+DOA used to represent microbial activity
C
      DO 475 K=0,KL
      OSCT(K)=0.0
      OSAT(K)=0.0
      DO 475 M=1,4
      OSCT(K)=OSCT(K)+OSC(M,K,L,NY,NX)
      OSAT(K)=OSAT(K)+OSA(M,K,L,NY,NX)
475   CONTINUE
      DO 480 K=0,KL
      IF(OSCT(K).GT.ZEROS(NY,NX))THEN
      DOSAK=DOSA(K)*AMAX1(0.0,ROQCK(K))
      DO 485 M=1,4
      OSA(M,K,L,NY,NX)=AMIN1(OSC(M,K,L,NY,NX)
     2,OSA(M,K,L,NY,NX)+DOSAK*OSC(M,K,L,NY,NX)/OSCT(K))
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
C     WRITE(*,8822)'OSA',I,J,NX,NY,L,K,M,OSA(M,K,L,NY,NX)
C    2,OSC(M,K,L,NY,NX),DOSA(K),ROQCK(K),DOSAK,OSAT(K),OSCT(K)
8822  FORMAT(A8,7I4,30E12.4)
C     ENDIF
485   CONTINUE
      ELSE
      DO 490 M=1,4
      OSA(M,K,L,NY,NX)=AMIN1(OSC(M,K,L,NY,NX),OSA(M,K,L,NY,NX))
490   CONTINUE
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
C     WRITE(*,8823)'OSC',I,J,L,K,((OMC(M,NGL,N,K,L,NY,NX),NGL=1,JG,N=1,7),M=1,3)
C    2,(ORC(M,K,L,NY,NX),M=1,2),OQC(K,L,NY,NX),OQCH(K,L,NY,NX)
C    3,OHC(K,L,NY,NX),OQA(K,L,NY,NX),OQAH(K,L,NY,NX),OHA(K,L,NY,NX)
C    4,(OSC(M,K,L,NY,NX),M=1,4)
8823  FORMAT(A8,4I4,100E12.4)
C     ENDIF
480   CONTINUE
C
C     AGGREGATE ALL TRANSFORMATIONS CALCULATED ABOVE FOR EACH N,K
C
      TRINH=0.0
      TRINO=0.0
      TRIPO=0.0
      TRIP1=0.0
      TRINB=0.0
      TRIOB=0.0
      TRIPB=0.0
      TRIB1=0.0
      TRGOM=0.0
      TRGOC=0.0
      TRGOD=0.0
      TRGOA=0.0
      TRGOH=0.0
      TUPOX=0.0
      TRDN3=0.0
      TRDNB=0.0
      TRDN2=0.0
      TRD2B=0.0
      TRDNO=0.0
      TRN2F=0.0
      DO 650 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 640 N=1,7
      DO 640 NGL=1,JG
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      TRINH=TRINH+RINH4(NGL,N,K)
      TRINO=TRINO+RINO3(NGL,N,K)
      TRIPO=TRIPO+RIPO4(NGL,N,K)
      TRIP1=TRIP1+RIP14(NGL,N,K)
      TRINB=TRINB+RINB4(NGL,N,K)
      TRIOB=TRIOB+RINB3(NGL,N,K)
      TRIPB=TRIPB+RIPOB(NGL,N,K)
      TRIB1=TRIB1+RIP1B(NGL,N,K)
      TRN2F=TRN2F+RN2FX(NGL,N,K)
      IF(L.EQ.NU(NY,NX))THEN
      TRINH=TRINH+RINH4R(NGL,N,K)
      TRINO=TRINO+RINO3R(NGL,N,K)
      TRIPO=TRIPO+RIPO4R(NGL,N,K)
      TRIP1=TRIP1+RIP14R(NGL,N,K)
      ENDIF
C     IF(NY.EQ.5.AND.L.EQ.10.AND.K.EQ.3.AND.N.EQ.2)THEN
C     WRITE(*,4469)'TRINH',I,J,NX,NY,L,K,NGL,N,TRINH,RINH4(NGL,N,K),RINH4R(NGL,N,K)
C     WRITE(*,4469)'TRIPO',I,J,NX,NY,L,K,NGL,N,TRIPO,RIPO4(NGL,N,K),RIPO4R(NGL,N,K)
C    2,CGOMP(NGL,N,K)
4469  FORMAT(A8,7I4,20E12.4)
C     ENDIF
      TRGOM=TRGOM+RCO2X(NGL,N,K)
      TRGOC=TRGOC+RCH4X(NGL,N,K)
      TRGOD=TRGOD+RGOMD(NGL,N,K)
      TUPOX=TUPOX+RUPOX(NGL,N,K)
      TRDN3=TRDN3+RDNO3(NGL,N,K)
      TRDNB=TRDNB+RDNOB(NGL,N,K)
      TRDN2=TRDN2+RDNO2(NGL,N,K)
      TRD2B=TRD2B+RDN2B(NGL,N,K)
      TRDNO=TRDNO+RDN2O(NGL,N,K)
      TRGOH=TRGOH+RH2GX(NGL,N,K)
C     IF(IYRC.EQ.2012.AND.I.EQ.151.AND.NX.EQ.1)THEN
C     WRITE(*,3333)'TRGOM',I,J,NX,NY,L,K,N,TRGOM
C    2,RCO2X(NGL,N,K),TRGOA,RGOMO(NGL,N,K),WFN(NGL,N,K),RGOMP
C     WRITE(*,3333)'TUPOX',I,J,NX,NY,L,K,N,TUPOX,RUPOX(NGL,N,K)
C     ENDIF
C     IF(J.EQ.12.AND.L.LE.4)THEN
C     WRITE(*,3333)'N2O',I,J,NX,NY,L,K,N,TRDN2,TRD2B,TRDNO
C    2,RDNO2(NGL,N,K),RDN2B(NGL,N,K),RDN2O(NGL,N,K),COXYS(L,NY,NX)
C    3,COXYG(L,NY,NX)
C     WRITE(*,3333)'TRGOH',I,J,NX,NY,L,K,N,TRGOH,RH2GX(NGL,N,K)
C    2,RGOMO(NGL,N,K)
3333  FORMAT(A8,7I4,20E12.4)
C     ENDIF
      ENDIF
640   CONTINUE
      ENDIF
650   CONTINUE
      DO 645 N=1,7
      DO 645 NGL=1,JG
      IF(N.LE.3.OR.N.EQ.5)THEN
      IF(N.NE.3)THEN
      TRGOA=TRGOA+CGOMC(NGL,N,5)
      ENDIF
      ENDIF
645   CONTINUE
C
C     ALLOCATE AGGREGATED TRANSFORMATIONS INTO ARRAYS TO UPDATE
C     STATE VARIABLES IN 'REDIST'
C
C     RCO2O=net CO2 uptake
C     TRGOA=total CO2 uptake by autotrophs
C     TRGOM total CO2 emission by heterotrophs reducing O2
C     TRGOD=total CO2 emission by denitrifiers reducing NOx
C     RVOXA(3)=CH4 oxidation
C     RCH4O=net CH4 uptake
C     CGOMC=total CH4 uptake by autotrophs
C     TRGOC=total CH4 emission
C     RH2GO=net H2 uptake
C     RH2GZ,TRGOH=total H2 uptake, emission
C     RUPOXO,TUPOX=total O2 uptake
C     RN2G=total N2 production
C     TRDNO=total N2O reduction
C     RN2O=total N2O uptake
C     TRDN2,TRD2B=total NO2 reduction in non-band,band
C     RCN2O,RCN2B=nitrous acid reduction in non-band,band
C
      RCO2O(L,NY,NX)=TRGOA-TRGOM-TRGOD
      RCH4O(L,NY,NX)=-TRGOC
      DO NGL=1,JG
      RCO2O(L,NY,NX)=RCO2O(L,NY,NX)-RVOXA(NGL,3)
      RCH4O(L,NY,NX)=RCH4O(L,NY,NX)+RVOXA(NGL,3)+CGOMC(NGL,3,5)
      ENDDO
      RH2GO(L,NY,NX)=RH2GZ-TRGOH
      RUPOXO(L,NY,NX)=TUPOX
      RN2G(L,NY,NX)=-TRDNO
      RN2O(L,NY,NX)=-TRDN2-TRD2B-RCN2O-RCN2B+TRDNO
C     IF(L.EQ.10)THEN
C     WRITE(*,2468)'RCO2O',I,J,NX,NY,L,RCO2O(L,NY,NX)
C    2,TRGOA,TRGOM,TRGOD,RVOXA(3),RCH4O(L,NY,NX)
C    3,CGOMC(3,5),TRGOC
C     WRITE(*,2468)'RN2O',I,J,NX,NY,L
C    2,RN2O(L,NY,NX),TRDN2,TRD2B,RCN2O,RCN2B,TRDNO
C    2,RCH4O(L,NY,NX),RVOXA(3)
C    2,CGOMC(3,5),TRGOC,(OMA(NGL,N,1),NGL=1,JG,N=1,7)
2468  FORMAT(A8,5I4,20E12.4)
C     ENDIF
C
C     XOQCS,XOQNZ,XOQPS,XOQAS=net change in DOC,DON,DOP,acetate
C
      DO 655 K=0,4
      DO 660 M=1,4
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)+RCOSC(M,K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)+RCOSN(M,K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)+RCOSP(M,K)
660   CONTINUE
      DO 665 M=1,2
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)+RDORC(M,K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)+RDORN(M,K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)+RDORP(M,K)
665   CONTINUE
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)+RDOHC(K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)+RDOHN(K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)+RDOHP(K)
      XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)+RDOHA(K)
      DO 670 N=1,7
      DO 670 NGL=1,JG
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-CGOQC(NGL,N,K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-CGOMN(NGL,N,K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-CGOMP(NGL,N,K)
      XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)-CGOAC(NGL,N,K)
     2+RCH3X(NGL,N,K)
670   CONTINUE
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-CSORP(K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-ZSORP(K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-PSORP(K)
      XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)-CSORPA(K)
655   CONTINUE
C
C     XNH4S,XNH4B=net change in NH4 in band,non-band
C     TRINH,TRINB=total NH4 mineraln-immobn in non-band,band
C     RVOXA(1),RVOXB(1)=total NH4 oxidation in non-band,band
C     XNO3S,XNO3B=net change in NO3 in band,non-band
C     TRINO,TRIOB=total NO3 immobn in non-band,band
C     RVOXA(2),RVOXB(2)=total NO2 oxidation in non-band,band
C     TRDN3,TRDNB=total NO3 reduction in non-band,band
C     RCNO3,RCN3B=NO3 production from nitrous acid reduction in non-band,band
C     XNO2S,XNO2B=net change in NO3 in band,non-band
C     TRDN2,TRD2B=total NO2 reduction in non-band,band
C     RCNO2,RCNOB=substrate-limited nitrous acid reduction in non-band,band
C     XH2PS,XH2BS=net change in H2PO4 in band,non-band
C     TRIPO,TRIPB=total H2PO4 mineraln-immobn in non-band,band
C     XH1PS,XH1BS=net change in HPO4 in band,non-band
C     TRIP1,TRIB1=total HPO4 mineraln-immobn in non-band,band
C     XN2GS=total N2 fixation
C     XZHYS=total H+ production
C     TRN2F=total N2 fixation
C
      XNH4S(L,NY,NX)=-TRINH
      XNO3S(L,NY,NX)=-TRINO-TRDN3+RCNO3
      XNO2S(L,NY,NX)=+TRDN3-TRDN2-RCNO2
      XH2PS(L,NY,NX)=-TRIPO
      XH1PS(L,NY,NX)=-TRIP1
      XNH4B(L,NY,NX)=-TRINB
      XNO3B(L,NY,NX)=-TRIOB-TRDNB+RCN3B
      XNO2B(L,NY,NX)=TRDNB-TRD2B-RCNOB
      DO NGL=1,JG
      XNH4S(L,NY,NX)=XNH4S(L,NY,NX)-RVOXA(NGL,1)
      XNO3S(L,NY,NX)=XNO3S(L,NY,NX)+RVOXA(NGL,2)
      XNO2S(L,NY,NX)=XNO2S(L,NY,NX)+RVOXA(NGL,1)-RVOXA(NGL,2)
      XNH4B(L,NY,NX)=XNH4B(L,NY,NX)-RVOXB(NGL,1)
      XNO3B(L,NY,NX)=XNO3B(L,NY,NX)+RVOXB(NGL,2)
      XNO2B(L,NY,NX)=XNO2B(L,NY,NX)+RVOXB(NGL,1)-RVOXB(NGL,2)
      ENDDO
      XH2BS(L,NY,NX)=-TRIPB
      XH1BS(L,NY,NX)=-TRIB1
      XN2GS(L,NY,NX)=TRN2F
      TFNQ(L,NY,NX)=TFNX
      VOLQ(L,NY,NX)=VOLWZ
C     IF(ISALTG.NE.0)THEN
C     XZHYS(L,NY,NX)=XZHYS(L,NY,NX)+0.1429*(RVOXA(1)+RVOXB(1)
C    2-TRDN3-TRDNB)-0.0714*(TRDN2+TRD2B+TRDNO)
C     ENDIF
C     IF(L.EQ.0)THEN
C     WRITE(*,2323)'XNH4S',I,J,L,XNH4S(L,NY,NX)
C    2,TRINH,RVOXA(1),VLNH4(L,NY,NX),TRDN2
C     WRITE(*,2323)'XNO3S',I,J,L,XNO3S(L,NY,NX)
C    2,TRINO,RVOXA(2),VLNO3(L,NY,NX),TRDN3,RCNO3
C     WRITE(*,2323)'XH2PS',I,J,L,XH2PS(L,NY,NX)
C    2,RIPOT,TRIPO,VLPO4(L,NY,NX)
C     WRITE(*,2323)'XNO2B',I,J,L,XNO2B(L,NY,NX),RVOXB(1)
C    2,VLNHB(L,NY,NX),RVOXB(2),VLNOB(L,NY,NX),TRDNB,TRD2B,RCNOB
C     ENDIF
C     WRITE(*,2324)'XOQCS',I,J,NX,NY,L,(XOQCS(K,L,NY,NX),K=0,4)
2324  FORMAT(A8,5I4,12E12.4)
      ELSE
      RCO2O(L,NY,NX)=0.0
      RCH4O(L,NY,NX)=0.0
      RH2GO(L,NY,NX)=0.0
      RUPOXO(L,NY,NX)=0.0
      RN2G(L,NY,NX)=0.0
      RN2O(L,NY,NX)=0.0
      XNH4S(L,NY,NX)=0.0
      XNO3S(L,NY,NX)=0.0
      XNO2S(L,NY,NX)=0.0
      XH2PS(L,NY,NX)=0.0
      XH1PS(L,NY,NX)=0.0
      XNH4B(L,NY,NX)=0.0
      XNO3B(L,NY,NX)=0.0
      XNO2B(L,NY,NX)=0.0
      XH2BS(L,NY,NX)=0.0
      XH1BS(L,NY,NX)=0.0
      XN2GS(L,NY,NX)=0.0
      ENDIF
C
C     MIX LITTER C BETWEEN ADJACENT SOIL LAYERS L AND LL
C
      IF(FOSCZ0.GT.ZERO)THEN
C     ORGR=total litter C
C     FOSCZ0=rate constant for mixing surface litter
C     FOSCXS=mixing fraction for surface litter
C     TOQCK=total active biomass respiration activity
C     TFNX=temperature function
C     VOLX=soil layer volume
C     OSCXD=mixing required for equilibrating litter concentration
C     FOSCXD=mixing fraction for equilibrating subsurface litter
C     FOSCXS=mixing fraction for subsurface litter
C
C     IF(I.EQ.116)THEN
C     WRITE(*,336)'LAYER',I,J,L,TOQCK(L,NY,NX),TOMA,TFNX,TOMA*TFNX
336   FORMAT(A8,3I4,20E12.4)
C     ENDIF
      IF(L.LT.NL(NY,NX))THEN
      IF(L.EQ.0)THEN
      LL=NU(NY,NX)
      IF(ORGR(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOSCXS=AMIN1(1.0,FOSCZ0/ORGR(L,NY,NX)*TOQCK(L,NY,NX))
      ELSE
      FOSCXS=0.0
      ENDIF
      ELSE
      DO 1100 LN=L+1,NL(NY,NX)
      IF(VOLX(LN,NY,NX).GT.ZEROS2(NY,NX))THEN
      LL=LN
      GO TO 1101
      ENDIF
1100  CONTINUE
1101  CONTINUE
      ORGRL=AMAX1(0.0,ORGR(L,NY,NX))
      ORGRLL=AMAX1(0.0,ORGR(LL,NY,NX))
      OSCXD=(ORGRL*VOLT(LL,NY,NX)-ORGRLL*VOLT(L,NY,NX))
     2/(VOLT(L,NY,NX)+VOLT(LL,NY,NX))
      IF(OSCXD.GT.0.0.AND.ORGR(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOSCXD=OSCXD/ORGR(L,NY,NX)
      ELSEIF(OSCXD.LT.0.0.AND.ORGR(LL,NY,NX).GT.ZEROS(NY,NX))THEN
      FOSCXD=OSCXD/ORGR(LL,NY,NX)
      ELSE
      FOSCXD=0.0
      ENDIF
      IF(VOLT(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      FOSCXS=FOSCZL*FOSCXD*TOQCK(L,NY,NX)/VOLT(L,NY,NX)
      ELSE
      FOSCXS=0.0
      ENDIF
      ENDIF
C     IF(L.EQ.3)THEN
C     WRITE(*,1115)'MIX',I,J,NX,NY,L,LL,FOSCXS,FOSCZ0
C    2,FOSCZL,FOSCXD,TOQCK(L,NY,NX),VOLT(L,NY,NX)
C    2,OSCXD,ORGR(L,NY,NX),ORGR(LL,NY,NX)
C    3,TKS(L,NY,NX)
1115  FORMAT(A8,6I4,30E12.4)
C     ENDIF
      IF(FOSCXS.GT.ZERO)THEN
      DO 7971 K=1,2
      DO 7961 N=1,7
      DO 7961 NGL=1,JG
      DO 7962 M=1,3
      IF(FOSCXS.GT.0.0)THEN
      OMCXS=FOSCXS*AMAX1(0.0,OMC(M,NGL,N,K,L,NY,NX))
      OMNXS=FOSCXS*AMAX1(0.0,OMN(M,NGL,N,K,L,NY,NX))
      OMPXS=FOSCXS*AMAX1(0.0,OMP(M,NGL,N,K,L,NY,NX))
      ELSE
      OMCXS=FOSCXS*AMAX1(0.0,OMC(M,NGL,N,K,LL,NY,NX))
      OMNXS=FOSCXS*AMAX1(0.0,OMN(M,NGL,N,K,LL,NY,NX))
      OMPXS=FOSCXS*AMAX1(0.0,OMP(M,NGL,N,K,LL,NY,NX))
      ENDIF
      OMC(M,NGL,N,K,L,NY,NX)=OMC(M,NGL,N,K,L,NY,NX)-OMCXS
      OMN(M,NGL,N,K,L,NY,NX)=OMN(M,NGL,N,K,L,NY,NX)-OMNXS
      OMP(M,NGL,N,K,L,NY,NX)=OMP(M,NGL,N,K,L,NY,NX)-OMPXS
      OMC(M,NGL,N,K,LL,NY,NX)=OMC(M,NGL,N,K,LL,NY,NX)+OMCXS
      OMN(M,NGL,N,K,LL,NY,NX)=OMN(M,NGL,N,K,LL,NY,NX)+OMNXS
      OMP(M,NGL,N,K,LL,NY,NX)=OMP(M,NGL,N,K,LL,NY,NX)+OMPXS
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,5558)'OMX',I,J,NX,NY,L,LL,K,N,M,OMC(M,NGL,N,K,L,NY,NX)
C    2,OMN(M,NGL,N,K,L,NY,NX),OMP(M,NGL,N,K,L,NY,NX),
C    2,OMC(M,NGL,N,K,LL,NY,NX)
C    2,OMN(M,NGL,N,K,LL,NY,NX),OMP(M,NGL,N,K,LL,NY,NX)
C    3,OMCXS,OMNXS,OMPXS,FOSCXS
5558  FORMAT(A8,9I4,12E12.4)
C     ENDIF
7962  CONTINUE
7961  CONTINUE
7971  CONTINUE

      DO 7901 K=1,2
      DO 7941 M=1,2
      IF(FOSCXS.GT.0.0)THEN
      ORCXS=FOSCXS*AMAX1(0.0,ORC(M,K,L,NY,NX))
      ORNXS=FOSCXS*AMAX1(0.0,ORN(M,K,L,NY,NX))
      ORPXS=FOSCXS*AMAX1(0.0,ORP(M,K,L,NY,NX))
      ELSE
      ORCXS=FOSCXS*AMAX1(0.0,ORC(M,K,LL,NY,NX))
      ORNXS=FOSCXS*AMAX1(0.0,ORN(M,K,LL,NY,NX))
      ORPXS=FOSCXS*AMAX1(0.0,ORP(M,K,LL,NY,NX))
      ENDIF
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-ORCXS
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-ORNXS
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-ORPXS
      ORC(M,K,LL,NY,NX)=ORC(M,K,LL,NY,NX)+ORCXS
      ORN(M,K,LL,NY,NX)=ORN(M,K,LL,NY,NX)+ORNXS
      ORP(M,K,LL,NY,NX)=ORP(M,K,LL,NY,NX)+ORPXS
C     IF(L.EQ.3.AND.K.EQ.2)THEN
C     WRITE(*,7942)'ORC',I,J,L,LL,K,M,ORC(M,K,L,NY,NX)
C    2,ORC(M,K,LL,NY,NX),ORCXS,FOSCXS
7942  FORMAT(A8,6I4,20E12.4)
C     ENDIF
7941  CONTINUE
      IF(FOSCXS.GT.0.0)THEN
      OQCXS=FOSCXS*AMAX1(0.0,OQC(K,L,NY,NX))
      OQCHXS=FOSCXS*AMAX1(0.0,OQCH(K,L,NY,NX))
      OHCXS=FOSCXS*AMAX1(0.0,OHC(K,L,NY,NX))
      OQAXS=FOSCXS*AMAX1(0.0,OQA(K,L,NY,NX))
      OQAHXS=FOSCXS*AMAX1(0.0,OQAH(K,L,NY,NX))
      OHAXS=FOSCXS*AMAX1(0.0,OHA(K,L,NY,NX))
      OQNXS=FOSCXS*AMAX1(0.0,OQN(K,L,NY,NX))
      OQNHXS=FOSCXS*AMAX1(0.0,OQNH(K,L,NY,NX))
      OHNXS=FOSCXS*AMAX1(0.0,OHN(K,L,NY,NX))
      OQPXS=FOSCXS*AMAX1(0.0,OQP(K,L,NY,NX))
      OQPHXS=FOSCXS*AMAX1(0.0,OQPH(K,L,NY,NX))
      OHPXS=FOSCXS*AMAX1(0.0,OHP(K,L,NY,NX))
      ELSE
      OQCXS=FOSCXS*AMAX1(0.0,OQC(K,LL,NY,NX))
      OQCHXS=FOSCXS*AMAX1(0.0,OQCH(K,LL,NY,NX))
      OHCXS=FOSCXS*AMAX1(0.0,OHC(K,LL,NY,NX))
      OQAXS=FOSCXS*AMAX1(0.0,OQA(K,LL,NY,NX))
      OQAHXS=FOSCXS*AMAX1(0.0,OQAH(K,LL,NY,NX))
      OHAXS=FOSCXS*AMAX1(0.0,OHA(K,LL,NY,NX))
      OQNXS=FOSCXS*AMAX1(0.0,OQN(K,LL,NY,NX))
      OQNHXS=FOSCXS*AMAX1(0.0,OQNH(K,LL,NY,NX))
      OHNXS=FOSCXS*AMAX1(0.0,OHN(K,LL,NY,NX))
      OQPXS=FOSCXS*AMAX1(0.0,OQP(K,LL,NY,NX))
      OQPHXS=FOSCXS*AMAX1(0.0,OQPH(K,LL,NY,NX))
      OHPXS=FOSCXS*AMAX1(0.0,OHP(K,LL,NY,NX))
      ENDIF
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-OQCXS
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)-OQCHXS
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-OHCXS
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-OQAXS
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)-OQAHXS
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-OHAXS
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-OQNXS
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)-OQNHXS
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-OHNXS
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-OQPXS
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)-OQPHXS
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-OHPXS
      OQC(K,LL,NY,NX)=OQC(K,LL,NY,NX)+OQCXS
      OQCH(K,LL,NY,NX)=OQCH(K,LL,NY,NX)+OQCHXS
      OHC(K,LL,NY,NX)=OHC(K,LL,NY,NX)+OHCXS
      OQA(K,LL,NY,NX)=OQA(K,LL,NY,NX)+OQAXS
      OQAH(K,LL,NY,NX)=OQAH(K,LL,NY,NX)+OQAHXS
      OHA(K,LL,NY,NX)=OHA(K,LL,NY,NX)+OHAXS
      OQN(K,LL,NY,NX)=OQN(K,LL,NY,NX)+OQNXS
      OQNH(K,LL,NY,NX)=OQNH(K,LL,NY,NX)+OQNHXS
      OHN(K,LL,NY,NX)=OHN(K,LL,NY,NX)+OHNXS
      OQP(K,LL,NY,NX)=OQP(K,LL,NY,NX)+OQPXS
      OQPH(K,LL,NY,NX)=OQPH(K,LL,NY,NX)+OQPHXS
      OHP(K,LL,NY,NX)=OHP(K,LL,NY,NX)+OHPXS
      DO 7931 M=1,4
      IF(FOSCXS.GT.0.0)THEN
      OSCXS=FOSCXS*AMAX1(0.0,OSC(M,K,L,NY,NX))
      OSAXS=FOSCXS*AMAX1(0.0,OSA(M,K,L,NY,NX))
      OSNXS=FOSCXS*AMAX1(0.0,OSN(M,K,L,NY,NX))
      OSPXS=FOSCXS*AMAX1(0.0,OSP(M,K,L,NY,NX))
      ELSE
      OSCXS=FOSCXS*AMAX1(0.0,OSC(M,K,LL,NY,NX))
      OSAXS=FOSCXS*AMAX1(0.0,OSA(M,K,LL,NY,NX))
      OSNXS=FOSCXS*AMAX1(0.0,OSN(M,K,LL,NY,NX))
      OSPXS=FOSCXS*AMAX1(0.0,OSP(M,K,LL,NY,NX))
      ENDIF
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-OSCXS
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-OSAXS
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-OSNXS
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-OSPXS
      OSC(M,K,LL,NY,NX)=OSC(M,K,LL,NY,NX)+OSCXS
      OSA(M,K,LL,NY,NX)=OSA(M,K,LL,NY,NX)+OSAXS
      OSN(M,K,LL,NY,NX)=OSN(M,K,LL,NY,NX)+OSNXS
      OSP(M,K,LL,NY,NX)=OSP(M,K,LL,NY,NX)+OSPXS
7931  CONTINUE
7901  CONTINUE
      ENDIF
      ENDIF
C     IF((I/1)*1.EQ.I.AND.J.EQ.19.AND.L.LE.5)THEN
C     WRITE(*,2123)'TOTALL',I,J,NX,NY,L,TFOXYX,TFNH4X
C    2,TFNO3X,TFPO4X,TFNH4B,TFNO3B,TFPO4B,TFNO2X,TFNO2B
C    3,TFOQC,TFOQA
2123  FORMAT(A8,5I4,12E15.4)
C     ENDIF
      ENDIF
      ELSE
      RCO2O(L,NY,NX)=0.0
      RCH4O(L,NY,NX)=0.0
      RH2GO(L,NY,NX)=0.0
      RUPOXO(L,NY,NX)=0.0
      RN2G(L,NY,NX)=0.0
      RN2O(L,NY,NX)=0.0
      XNH4S(L,NY,NX)=0.0
      XNO3S(L,NY,NX)=0.0
      XNO2S(L,NY,NX)=0.0
      XH2PS(L,NY,NX)=0.0
      XH1PS(L,NY,NX)=0.0
      XNH4B(L,NY,NX)=0.0
      XNO3B(L,NY,NX)=0.0
      XNO2B(L,NY,NX)=0.0
      XH2BS(L,NY,NX)=0.0
      XH1BS(L,NY,NX)=0.0
      XN2GS(L,NY,NX)=0.0
      ENDIF
998   CONTINUE
C     WRITE(20,3434)'RN2O',IYRC,I,J,(RN2O(L,NY,NX),L=0,NL(NY,NX))
3434  FORMAT(A8,3I4,20E12.4)
C
C     SOC LOSS IF FIRE OR REMOVAL EVENT IS ENTERED IN DISTURBANCE FILE
C
      IF(J.EQ.INT(ZNOON(NY,NX)).AND.(ITILL(I,NY,NX).EQ.21
     2.OR.ITILL(I,NY,NX).EQ.22))THEN
      IF(ITILL(I,NY,NX).EQ.22)THEN
      IFLGS(NY,NX)=1
      IFLGJ=0
      NLL=-1
      DO 2945 L=0,NL(NY,NX)
C     WRITE(*,9494)'FIRE',I,J,L,IFLGJ,NLL,ORGC(L,NY,NX),THETW(L,NY,NX)
C    2,FVLWB,CORGC(L,NY,NX),FORGC,DPTH(L,NY,NX),BKDS(L,NY,NX)
C    3,VOLY(L,NY,NX),DTBLX(NY,NX),DCORP(I,NY,NX)
9494  FORMAT(A8,5I6,12E12.4)
      IF(L.EQ.0.OR.L.GE.NUM(NY,NX))THEN
      IF(IFLGJ.EQ.1)THEN
      GO TO 2946
      ELSEIF(THETW(L,NY,NX).GT.FVLWB.OR.CORGC(L,NY,NX).LE.FORGC)THEN
      IFLGJ=1
      ELSE
      NLL=L
      ENDIF
      ENDIF
2945  CONTINUE
      ELSE
      NLL=0
      ENDIF
2946  CONTINUE
      DO 2950 L=0,NLL
      IF(NLL.GE.0)THEN
      IF(ITILL(I,NY,NX).EQ.22)THEN
      IF(L.EQ.0)THEN
      FORGCX=0.0
      ELSE
      FORGCX=FORGC
      ENDIF
      DCORPC=AMIN1(0.999,DCORP(I,NY,NX))*(CORGC(L,NY,NX)-FORGCX)
     2/(AMAX1(CORGC(L,NY,NX),0.55E+06)-FORGCX)
      ELSE
      DCORPC=AMIN1(0.999,DCORP(I,NY,NX))
      ENDIF
C     VOLWOU=VOLWOU+DCORPC*VOLW(L,NY,NX)
C     HEATOU=HEATOU+DCORPC*4.19*TKS(L,NY,NX)*VOLW(L,NY,NX)
C     VOLW(L,NY,NX)=VOLW(L,NY,NX)-DCORPC*VOLW(L,NY,NX)
C     WRITE(*,9696)'BURN',I,J,L,NLL,ITILL(I,NY,NX)
C    2,CORGC(L,NY,NX),ORGC(L,NY,NX)
C    2,FORGCX,DCORPC,DCORP(I,NY,NX),VOLW(L,NY,NX),BKDS(L,NY,NX)
9696  FORMAT(A8,5I6,12E12.4)
      OC=0.0
      ON=0.0
      OP=0.0
      DC=0.0
      DN=0.0
      DP=0.0
      DO 2955 K=0,4
      DO 2955 M=1,4
      ONL(M,K)=0.0
      OPL(M,K)=0.0
2955  CONTINUE
      DO 2970 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
C
C     REMOVE MICROBIAL BIOMASS
C
      DO 2960 N=1,7
      DO 2960 NGL=1,JG
      DO 2960 M=1,3
      OCH=DCORPC*OMC(M,NGL,N,K,L,NY,NX)
      ONH=DCORPC*OMN(M,NGL,N,K,L,NY,NX)
      OPH=DCORPC*OMP(M,NGL,N,K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSEIF(K.LE.4)THEN
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ELSEIF(K.EQ.5)THEN
      ONL(4,1)=ONL(4,1)+ONH-ONX
      OPL(4,1)=OPL(4,1)+OPH-OPX
      ENDIF
      OMC(M,NGL,N,K,L,NY,NX)=OMC(M,NGL,N,K,L,NY,NX)-OCH
      OMN(M,NGL,N,K,L,NY,NX)=OMN(M,NGL,N,K,L,NY,NX)-ONH
      OMP(M,NGL,N,K,L,NY,NX)=OMP(M,NGL,N,K,L,NY,NX)-OPH
      DC=DC+OMC(M,NGL,N,K,L,NY,NX)
      DN=DN+OMN(M,NGL,N,K,L,NY,NX)
      DP=DP+OMP(M,NGL,N,K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONX
      OP=OP+OPX
2960  CONTINUE
      ENDIF
2970  CONTINUE

C
C     REMOVE MICROBIAL RESIDUE
C
      DO 2900 K=0,4
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 2940 M=1,2
      OCH=DCORPC*ORC(M,K,L,NY,NX)
      ONH=DCORPC*ORN(M,K,L,NY,NX)
      OPH=DCORPC*ORP(M,K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSE
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ENDIF
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-OCH
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-ONH
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-OPH
      DC=DC+ORC(M,K,L,NY,NX)
      DN=DN+ORN(M,K,L,NY,NX)
      DP=DP+ORP(M,K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONX
      OP=OP+OPX
2940  CONTINUE
C
C     REMOVE DOC, DON, DOP
C
      OCH=DCORPC*OQC(K,L,NY,NX)
      OCA=DCORPC*OQA(K,L,NY,NX)
      ONH=DCORPC*OQN(K,L,NY,NX)
      OPH=DCORPC*OQP(K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSE
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ENDIF
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-OCH
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-OCA
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-ONH
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-OPH
      OC=OC+OCH+OCA
      ON=ON+ONX
      OP=OP+OPX
      OCH=DCORPC*OQCH(K,L,NY,NX)
      ONH=DCORPC*OQNH(K,L,NY,NX)
      OPH=DCORPC*OQPH(K,L,NY,NX)
      OAH=DCORPC*OQAH(K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSE
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ENDIF
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)-OCH
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)-ONH
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)-OPH
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)-OAH
      OC=OC+OCH+OAH
      ON=ON+ONX
      OP=OP+OPX
C
C     REMOVE ADSORBED OM
C
      OCH=DCORPC*OHC(K,L,NY,NX)
      ONH=DCORPC*OHN(K,L,NY,NX)
      OPH=DCORPC*OHP(K,L,NY,NX)
      OAH=DCORPC*OHA(K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      IF(K.LE.2)THEN
      ONL(4,K)=ONL(4,K)+ONH-ONX
      OPL(4,K)=OPL(4,K)+OPH-OPX
      ELSE
      ONL(1,K)=ONL(1,K)+ONH-ONX
      OPL(1,K)=OPL(1,K)+OPH-OPX
      ENDIF
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-OCH
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-ONH
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-OPH
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-OAH
      DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
     2+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      DP=DP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONX
      OP=OP+OPX
C
C     REMOVE RESIDUE
C
      DO 2930 M=1,4
      OCH=DCORPC*OSC(M,K,L,NY,NX)
      OCA=DCORPC*OSA(M,K,L,NY,NX)
      ONH=DCORPC*OSN(M,K,L,NY,NX)
      OPH=DCORPC*OSP(M,K,L,NY,NX)
      ONX=EFIRE(1,ITILL(I,NY,NX))*ONH
      OPX=EFIRE(2,ITILL(I,NY,NX))*OPH
      ONL(M,K)=ONL(M,K)+ONH-ONX
      OPL(M,K)=OPL(M,K)+OPH-OPX
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-OCH
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-OCA
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-ONH
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-OPH
      DC=DC+OSC(M,K,L,NY,NX)
      DN=DN+OSN(M,K,L,NY,NX)
      DP=DP+OSP(M,K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONX
      OP=OP+OPX
2930  CONTINUE
      ENDIF
2900  CONTINUE
C
C     ADD UNBURNED N,P TO ORG N, ORG P
C
      DO 2905 K=0,4
      DO 2905 M=1,4
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)+ONL(M,K)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)+OPL(M,K)
      DN=DN+ONL(M,K)
      DP=DP+OPL(M,K)
2905  CONTINUE
C
C     REMOVE FERTILIZER IN RESIDUE
C
      IF(ITILL(I,NY,NX).EQ.21)THEN
      ON=ON+DCORPC*(ZNH4S(L,NY,NX)+ZNH3S(L,NY,NX)
     2+ZNO3S(L,NY,NX)+ZNO2S(L,NY,NX))
      OP=OP+DCORPC*(H1PO4(L,NY,NX)+H2PO4(L,NY,NX))
      ZNH4S(L,NY,NX)=(1.0-DCORPC)*ZNH4S(L,NY,NX)
      ZNH3S(L,NY,NX)=(1.0-DCORPC)*ZNH3S(L,NY,NX)
      ZNO3S(L,NY,NX)=(1.0-DCORPC)*ZNO3S(L,NY,NX)
      ZNO2S(L,NY,NX)=(1.0-DCORPC)*ZNO2S(L,NY,NX)
      H1PO4(L,NY,NX)=(1.0-DCORPC)*H1PO4(L,NY,NX)
      H2PO4(L,NY,NX)=(1.0-DCORPC)*H2PO4(L,NY,NX)
      XN4(L,NY,NX)=(1.0-DCORPC)*XN4(L,NY,NX)
      PALPO(L,NY,NX)=(1.0-DCORPC)*PALPO(L,NY,NX)
      PFEPO(L,NY,NX)=(1.0-DCORPC)*PFEPO(L,NY,NX)
      PCAPD(L,NY,NX)=(1.0-DCORPC)*PCAPD(L,NY,NX)
      PCAPH(L,NY,NX)=(1.0-DCORPC)*PCAPH(L,NY,NX)
      PCAPM(L,NY,NX)=(1.0-DCORPC)*PCAPM(L,NY,NX)
      ZNH4FA(L,NY,NX)=(1.0-DCORPC)*ZNH4FA(L,NY,NX)
      ZNH3FA(L,NY,NX)=(1.0-DCORPC)*ZNH3FA(L,NY,NX)
      ZNHUFA(L,NY,NX)=(1.0-DCORPC)*ZNHUFA(L,NY,NX)
      ZNO3FA(L,NY,NX)=(1.0-DCORPC)*ZNO3FA(L,NY,NX)
      ENDIF
      ORGC(L,NY,NX)=DC
      ORGN(L,NY,NX)=DN
      IF(L.EQ.0)THEN
      HFLXD=4.19E-06*(ORGCX(L,NY,NX)-ORGC(L,NY,NX))*TKS(L,NY,NX)
      HEATOU=HEATOU+HFLXD
      ENDIF
C     IF(L.EQ.0)THEN
C     VHCP(0,NY,NX)=2.496E-06*ORGC(0,NY,NX)+4.19*VOLW(0,NY,NX)
C    2+1.9274*VOLI(0,NY,NX)
C     ELSE
C     VHCP(L,NY,NX)=VHCM(L,NY,NX)+4.19*(VOLW(L,NY,NX)+VOLWH(L,NY,NX))
C    2+1.9274*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
C     ENDIF
      IF(ITILL(I,NY,NX).EQ.21)THEN
      TCOU=TCOU+OC
      TZOU=TZOU+ON
      TPOU=TPOU+OP
      UDOCQ(NY,NX)=UDOCQ(NY,NX)+OC
      UDONQ(NY,NX)=UDONQ(NY,NX)+ON
      UDOPQ(NY,NX)=UDOPQ(NY,NX)+OP
      TNBP(NY,NX)=TNBP(NY,NX)-OC
      ELSEIF(ITILL(I,NY,NX).EQ.22)THEN
      CO2GIN=CO2GIN-OC
      OXYGIN=OXYGIN+2.667*OC
      OXYGOU=OXYGOU+2.667*OC
      TZOU=TZOU+ON
      TPOU=TPOU+OP
      UCO2F(NY,NX)=UCO2F(NY,NX)-(1.0-FCH4F)*OC
      UCH4F(NY,NX)=UCH4F(NY,NX)-FCH4F*OC
      UOXYF(NY,NX)=UOXYF(NY,NX)+(1.0-FCH4F)*2.667*OC
      UNH3F(NY,NX)=UNH3F(NY,NX)-ON
      UN2OF(NY,NX)=UN2OF(NY,NX)-0.0
      UPO4F(NY,NX)=UPO4F(NY,NX)-OP
      TNBP(NY,NX)=TNBP(NY,NX)-OC
      ENDIF
      ENDIF
2950  CONTINUE
      ENDIF
9990  CONTINUE
9995  CONTINUE
      RETURN
      END
