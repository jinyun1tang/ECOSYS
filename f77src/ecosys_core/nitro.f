      SUBROUTINE nitro(I,J,NFZ,NHW,NHE,NVN,NVS)
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
      DIMENSION CNOMA(7,0:5),CPOMA(7,0:5),OMA(7,0:5),FOMA(7,0:5)
     2,FOMN(7,0:5),RDOSC(4,0:4),RDOSN(4,0:4),RDOSP(4,0:4),RHOSC(4,0:4)
     3,RHOSN(4,0:4),RHOSP(4,0:4),RCOSC(4,0:4),RCOSN(4,0:4),RCOSP(4,0:4)
     4,SPOSC(4,0:4),RDORC(2,0:4),RDORN(2,0:4),RDORP(2,0:4),SPORC(2)
     5,RDOHC(0:4),RDOHN(0:4),RDOHP(0:4),RDOHA(0:4),CSORP(0:4),ZSORP(0:4)
     6,PSORP(0:4),CSORPA(0:4),RUPOX(7,0:5),RGN2F(7,0:5)
     8,RGOMO(7,0:5),ROXYM(7,0:5),ROXYP(7,0:5),ROXYO(7,0:5) 
     9,RDNO3(7,0:5),RDNOB(7,0:5),RDNO2(7,0:5),RDN2B(7,0:5),RDN2O(7,0:5) 
     1,RGOMD(7,0:5),RMOMC(2,7,0:5),RINH4(7,0:5),RINO3(7,0:5) 
     2,RIPO4(7,0:5),RINB4(7,0:5),RINB3(7,0:5),RIPOB(7,0:5),FOMK(7,0:5)
     3,RDOMC(2,7,0:5),RDOMN(2,7,0:5),RDOMP(2,7,0:5),RHOMC(2,7,0:5)
     4,RHOMN(2,7,0:5),RHOMP(2,7,0:5),RCOMC(2,7,0:5),RCOMN(2,7,0:5)
     5,RCOMP(2,7,0:5),CGOMC(7,0:5),CGOMN(7,0:5),RH2GX(7,0:5)
     6,CGOMP(7,0:5),RDMMC(2,7,0:5),RHMMC(2,7,0:5),RCMMC(2,7,0:5)
     7,RDMMN(2,7,0:5),RHMMN(2,7,0:5),RCMMN(2,7,0:5),RDMMP(2,7,0:5)
     8,RHMMP(2,7,0:5),RCMMP(2,7,0:5),RCCMC(2,7,0:4)
     9,RCCMN(2,7,0:4),RCCMP(2,7,0:4),RN2FX(7,0:5),TOMK(0:5) 
     1,TONK(0:5),TOPK(0:5),SPOMC(2),OMC2(7,0:5) 
     2,OMN2(7,0:5),FOM2(7,0:5),FOCA(0:4),FOAA(0:4),RXOMC(2,7,0:5)
     3,RXOMN(2,7,0:5),RXOMP(2,7,0:5),R3OMC(2,7,0:5),R3OMN(2,7,0:5)
     4,R3OMP(2,7,0:5),RXMMC(2,7,0:5),RXMMN(2,7,0:5),RXMMP(2,7,0:5) 
     5,R3MMC(2,7,0:5),R3MMN(2,7,0:5),R3MMP(2,7,0:5),WFN(7,0:5)
     6,TFNG(7,0:5),TFNR(7,0:5),OSCH(0:4),OSAH(0:4)  
      DIMENSION CGOQC(7,0:5),CGOAC(7,0:5),ROQCK(0:4),XOQCK(0:4)  
     2,EN2F(7),ORCT(0:4),OSCT(0:4),OSAT(0:4),ZNH4T(0:JZ),ZNO3T(0:JZ)
     3,ZNO2T(0:JZ),H2P4T(0:JZ),RINH4R(7,0:5),RINO3R(7,0:5)
     4,RIPO4R(7,0:5),FNH4XR(7,0:5),FNO3XR(7,0:5),FPO4XR(7,0:5)
     5,RGOMY(7,0:5),CNQ(0:4),CPQ(0:4),CNH(0:4),CPH(0:4)
     6,CNS(4,0:4),CPS(4,0:4),ROQCD(7,0:4),FORC(0:5),SPOMK(2),RMOMK(2) 
     8,CGOMS(2,7,0:5),CGONS(2,7,0:5),CGOPS(2,7,0:5),H1P4T(0:JZ) 
     1,TONX(0:5),TOPX(0:5),FCNK(0:4),FCPK(0:4),FP14XR(7,0:5) 
     2,RCO2X(7,0:5),RCH3X(7,0:5),RCH4X(7,0:5),RVOXA(7),RVOXB(7)
     2,XOQCZ(0:4),XOQNZ(0:4),XOQPZ(0:4),XOQAZ(0:4)
     3,XOMCZ(3,7,0:4),XOMNZ(3,7,0:4),XOMPZ(3,7,0:4)    
     4,FCN(7,0:5),FCP(7,0:5),FCNP(7,0:5),RIP14(7,0:5),RIP1B(7,0:5)
     5,RIP14R(7,0:5),SPCMB(0:5),FRCBCO(0:5) 
      DIMENSION DOSA(0:4) 
C
C     SUBSTRATE DECOMPOSITION BY MICROBIAL POPULATIONS
C
C     ORAD=microbial radius (m)
C     BIOS=microbial density (n m-3)
C     BIOA=microbial surface area (m2 m-3)
C     DCKI=inhibition of decomposition by microbial activity
C        vs substrate (g C m-3 t-1)
C     RCCX=maximum remobilization of microbial N (-)
C     RCCY=maximum remobilization of microbial P (-)
C     RCCZ,RCCY=minimum, maximum remobilization of microbial C (-)
C     FPRIM, FPRIMM=fraction of nonstructural, microbial C,N,P
C        transferred with priming (-)
C     OMGR=rate constant for transferring nonstructural to 
C       structural microbial C (h-1)
C     OQKI=DOC product inhibition constant for decomposition (g C m-3)
C     H2KI=H2 product inhibition for methanogenesis (g H m-3)
C     OAKI=acetate product inhibition for methanogenesis (g C m-3) 
C     COMKI,COMKM=inhibition constant for microbial decomposition,
C        maintenance (g C m-3)
C     respiration with low microbial C (g micr C g-1 subs C)
C     FPRIMB=rate constant for mixing soil microbial biomass 
C        among soil layers (h-1)
C     FMN=minimum ratio of total biological demand for any substrate
C        by any microbial population
C     DCKM=Km for SOC decomposition (g C g-1 soil)
C     TCMCX=maximum value of Arrhenius function for combustion
C     CNKI,CPKI=nonstructural N,P inhibition constant on microbial
C        nutrient recycling (g N,P g-1 C)
C     
      PARAMETER (ORAD=1.0E-06,BIOS=1.0E-06/(4.19*ORAD**3) 
     2,BIOA=BIOS*12.57*ORAD**2,DCKI=2.5,FPRIM=5.0E-02
     3,FPRIMM=1.0E-06,RCCZ=0.167,RCCY=0.333,RCCX=0.833,RCCQ=0.833
     4,OMGR=2.5E-01,OQKI=2.4E+03,H2KI=1.0,OAKI=12.0,COMKI=2.5E-03
     5,COMKM=1.0E-04,FPRIMB=2.5E-08,FMN=1.0E-03,DCKM=1.0E+03
     7,TFNCX=2.0)
      PARAMETER(CNKI=1.0E-01,CPKI=1.0E-02)
C
C     SPECIFIC RESPIRATION RATES, M-M UPTAKE CONSTANTS,
C     STOICHIOMETRIC CONSTANTS FOR MICROBIAL REDOX REACTIONS
C
C     VMX*=specific oxidation rates (g C g-1 C h-1)
C        O=all bacteria, 
C        F=fungi
C        M=acetotrophic methanogens
C        H=ammonia oxidizers
C        N=nitrite oxidizers
C        4=methanotrophs
C        C=hydrogenotrophic methanogens
C     OQK*=Km for DOC,acetate uptake by heterotrophs (g C m-3)
C        M=DOC by all bacteria and fungi 
C        A=acetate by fermenters
C        AM=acetate by acetotrophic methanogens
C     CCKM=Km for CO2 uptake(g C m-3)
C     CCK4=Km for CH4 uptake (g C m-3)
C     Z*KM=Km for N uptake (g N m-3)
C        H=NH4 by nitrifiers, N=NO2 by nitrifiers
C        3=NO3 by denitrifiers, 2=NO2 by denitrifiers
C        1=N2O uptake by denitrifiers
C     Z4*=NH4 uptake kinetics by all MFTs(g N m-2 h-1, g N m-3)
C        MX=maximum uptake rate, 
C        KU=Km
C        MN= minimum concentration  
C     ZO*=NO3 uptake kinetics by all MFTs(g N m-2 h-1, g N m-3)
C        MX=maximum uptake rate
C        KU=Km
C        MN= minimum concentration  
C     HP*=H2PO4 uptake kinetics by all MFTs(g P m-2 h-1, g P m-3)
C        MX=maximum uptake rate
C        KU=Km
C        MN= minimum concentration
C     ZFKM=Km for N2 uptake by diazotrophs (g N m-3)
C     H2KM=Km for H2 uptake by hydrogenotrophic methanogens (g H m-3)  
C     ECNH=efficiency CO2 conversion to biomass by ammonia oxidizers 
C        (g C g N-1) 
C     ECNO=efficiency CO2 conversion to biomass by nitrite oxidizers
C        (g C g N-1) 
C     ECHO=efficiency CO2 conversion to biomass by methane oxidizers
C        (g C g C-1) 
C     ECN3,ECN2,ECN1=C:N ratios for e- transfers to NO3, NO2 and N2O
C        by denitrifiers (g C g N-1) 
C     RNFNI=rate constant for decline in nitrification inhibition (h-1)
C     ZHKI=inhibition of nitrification inhibition by NH3 (g N m-3) 
C     VMKI=product inhibition for NOx reduction by denitrifiers
C        (g N m-3 t-1)
C     VHKI=product inhibition for NH3 oxidation by nitrifiers 
C        (g N m-3)  
C
      PARAMETER (VMXO=0.125,VMXF=0.125,VMXM=0.125,VMXH=0.375
     2,VMXN=0.25,VMX4=0.375,VMXC=0.125,OQKM=1.2E+01,OQKA=1.2E+01
     3,OQKAM=1.2E+01,CCKM=0.15,CCK4=1.2E-03,ZHKM=1.4
     4,ZNKM=1.4,Z3KM=1.4,Z2KM=1.4,Z1KM=0.014,Z4MX=2.5E-03
     5,Z4KU=0.40,Z4MN=0.0125,ZOMX=2.5E-03,ZOKU=0.35,ZOMN=0.03
     7,HPMX=2.5E-03,HPKU=0.075,HPMN=0.002,ZFKM=0.14,H2KM=0.01
     8,ECNH=0.30,ECNO=0.10,ECHO=0.75,ECN3=0.857,ECN2=0.857,ECN1=0.429
     9,RNFNI=2.0E-04,ZHKI=7.0E+03,VMKI=0.25,VHKI=15.0)
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
C     E*=growth respiration requirement (g C g C-1)(growth yield=1.0-E*)
C        N2X=aerobic N2 fixation, N2Y=anaerobic N2 fixation
C        O2X=aerobic bacteria (DOC), H4X=fermenters, O2G=fungi
C        O2D=denitrifiers (aerobic), NFX=diazotrophs
C        NOX=denitrifiers (anaerobic),O2A=aerobic bacteria (acetate)
C
      PARAMETER (EOMC=25.0,EOMD=37.5,EOMG=37.5,EOMF=75.0,EOMH=25.0
     2,EOMN=112.5,GO2X=37.5,GH4X=66.5,GCHX=4.50
     3,GO2A=GO2X-GCHX,GC4X=3.00,GCOX=11.00,GNOX=10.0
     3,GN2X=150.0,EN2X=GO2X/GN2X,EN2Y=GCHX/GN2X
     4,EO2X=1.0/(1.0+GO2X/EOMC),EH4X=1.0/(1.0+GH4X/EOMC)
     5,EO2G=1.0/(1.0+GO2X/EOMG),EO2D=1.0/(1.0+GO2X/EOMD)
     6,ENFX=1.0/(1.0+GO2X/EOMN),ENOX=1.0/(1.0+GNOX/EOMC)
     7,EO2A=1.0/(1.0+GO2A/EOMC))
C
C     SORPTION COEFFICIENTS 
C
C     TSORP,HSORP=sorption rate constant (h-1),isotherm shape
C        parameter for DOC,DON,DOP,DOA adsorption
C     
      PARAMETER (TSORP=0.1,HSORP=1.0E+00)
C
C     SPECIFIC DECOMPOSITION RATES
C
C     DOSA=rate constant for litter colonization by heterotrophs 
C        (g C g-1 C)
C     SP*= specific decomposition rate constant (g subs. C g-1 micr. C)
C        :OHC=adsorbed SOC
C        :OHA=adsorbed acetate (g)
C        :OSC=bulk SOC (g) 
C           K=0,M=1,4 woody litter
C           K=1,M=1,4 non-woody litter,
C           K=2,M=1,4 manure
C           K=3,M=1,1 POC
C           K=4,M=1,2 humus
C        :ORC (M=1,2) microbial residue
C        :OMC (M=1,2) microbial biomass
C     RMOM=specific maintenance respiration (g C g-1 N h-1)
C     EN2X,EN2Y=energy efficiency of N2 fixation by aerobic,anaerobic
C        diazotrophs (g N g C-1)  
C
      PARAMETER (SPOHC=0.25,SPOHA=0.25,RMOM=0.010)
      DATA DOSA/0.10,0.25,0.25,0.25,0.25/
      DATA SPOSC/6.75,6.75,1.35,0.45
     1,6.75,6.75,1.35,0.45
     2,6.75,6.75,1.35,0.45
     3,0.045,0.0,0.0,0.0
     4,0.045,0.015,0.0,0.0/
      DATA SPORC/6.75,1.35/
      DATA SPOMC/1.0E-02,1.0E-03/
      DATA EN2F/0.0,0.0,0.0,0.0,0.0,EN2X,EN2Y/
C
C     SPECIFIC COMBUSTION RATES
C
C     SPCMB=specific combustion rate of SOC(K=0,5)at 600K (g C m-2 h-1)
C     SPCMBH=specific combustion rate of charcoal at 700K (g C m-2 h-1)
C
      DATA SPCMB/5.0E+02,1.0E+03,5.0E+02,5.0E+02,5.0E+02,5.0E+02/
      DATA SPCMBH/5.0E+02/
      REAL*4 WFNG,TFNX,TFNY,TFNG,TFNR,CNSHZ,CPSHZ,FRM
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
C
C     TOTAL MINERAL N
C
C     ZNH4T,ZNH4S,ZNH4B=total,NH4 in non-band,band zones (g N)
C     ZNO3T,ZNO3S,ZNO3B=total,NO3 in non-band,band zones (g N)
C     ZNO2T,ZNO2S,ZNO2B=total,NO2 in non-band,band zones (g N)
C     H1P4T,H1PO4,H1POB=total,HPO4 in non-band,band zones (g P)
C     H2P4T,H2PO4,H2POB=total,H2PO4 in non-band,band zones (g P)
C     THETWY,VOLWY=water concentration,volume used to calculate
C        aqueous microbial concentrations that drive microbial density
C        effects on decomposition (m3 m-3,m3)
C     PSISM,PSISO=matric,osmotic water potential (MPa)
C
      DO 998 L=0,NL(NY,NX)
      PSISG=PSISM(L,NY,NX)+PSISO(L,NY,NX)
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IF(L.EQ.0.OR.L.GE.NU(NY,NX))THEN
      IF(L.EQ.0)THEN
      KL=2
      ZNH4T(NU(NY,NX))=AMAX1(0.0,ZNH4S(NU(NY,NX),NY,NX)
     2+ZNH4B(NU(NY,NX),NY,NX))
      ZNO3T(NU(NY,NX))=AMAX1(0.0,ZNO3S(NU(NY,NX),NY,NX)
     2+ZNO3B(NU(NY,NX),NY,NX))
      ZNO2T(NU(NY,NX))=AMAX1(0.0,ZNO2S(NU(NY,NX),NY,NX)
     2+ZNO2B(NU(NY,NX),NY,NX))
      H1P4T(NU(NY,NX))=AMAX1(0.0,H1PO4(NU(NY,NX),NY,NX)
     2+H1POB(NU(NY,NX),NY,NX))
      H2P4T(NU(NY,NX))=AMAX1(0.0,H2PO4(NU(NY,NX),NY,NX)
     2+H2POB(NU(NY,NX),NY,NX))
      IF(VOLWRX(NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=VOLW(0,NY,NX)/VOLWRX(NY,NX)
      THETWY=AMAX1(0.0,(AMIN1(FC(L,NY,NX),THETWR)-THETY(L,NY,NX)))
      VOLWY=THETWY/(1.0+THETWY)*VOLWRX(NY,NX)
C     IF((I/30)*30.EQ.I.AND.NFZ.EQ.1.AND.J.EQ.15.AND.L.EQ.0)THEN
C     WRITE(*,8825)'THETWY0',I,J,NFZ,NX,NY,L,THETWR,THETWY
C    2,VOLWY,VOLWRX(NY,NX),VOLW(0,NY,NX)/VOLR(NY,NX)
C    2,VOLW(0,NY,NX),POROS(L,NY,NX),FC(0,NY,NX),WP(0,NY,NX)
C    3,THETY(L,NY,NX),PSISM(0,NY,NX),ORGC(0,NY,NX),VOLR(NY,NX)
8825  FORMAT(A8,6I4,20E12.4)
C     ENDIF
      ELSE
      VOLWY=0.0
      ENDIF
      ELSE
      KL=4
      THETWY=AMAX1(0.0,(AMIN1(AMAX1(0.5*POROS(L,NY,NX),FC(L,NY,NX))
     2,THETW(L,NY,NX))-THETY(L,NY,NX)))
      VOLWY=THETWY*VOLY(L,NY,NX)
C     IF((I/120)*120.EQ.I.AND.J.EQ.24.AND.L.LE.6)THEN
C     WRITE(*,8824)'THETWYL',I,J,NFZ,NX,NY,L,THETWY,THETW(L,NY,NX)
C    2,VOLWY,POROS(L,NY,NX),FC(L,NY,NX),WP(L,NY,NX),THETY(L,NY,NX)
C    3,VOLW(L,NY,NX),VOLWH(L,NY,NX),VOLY(L,NY,NX),VOLT(L,NY,NX)
C    4,DTBLX(NY,NX)
8824  FORMAT(A8,6I4,20E12.4)
C     ENDIF
      ENDIF
C
C     TEMPERATURE FUNCTIONS FOR GROWTH AND MAINTENANCE
C     WITH OFFSET FOR THERMAL ADAPTATION
C
C     TKS=soil temperature (K)
C     OFFSET=adjustment for acclimation based on MAT in ‘starts.f’
C     8.313,710.0=gas constant,enthalpy (J mol-1 K-1)
C     62500=activation energy (J mol-1)
C     197500=low temp inactivation for growth and maintenance
C        respiration (J mol-1)
C     222500=high temp inactivation for growth respiration (J mol-1)
C     TFNX,TFNY=temperature function for growth,maintenance
C        respiration 
C
      TKSO=TKS(L,NY,NX)+OFFSET(NY,NX)
      RTK=8.3143*TKSO
      STK=710.0*TKSO
      ACTV=1.0+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFNX=EXP(25.229-62500/RTK)/ACTV
      ACTVM=1.0+EXP((197500-STK)/RTK)
      TFNY=AMIN1(1.0E+03,EXP(25.216-62500/RTK)/ACTVM)
C
C     OXYI=inhibition of fermenters by O2
C     COXYS=aqueous O2 concentration (g m-3)
C     ORGCL=SOC used to calculate microbial concentration (g C Mg-1)
C     BKVL=soil layer mass (Mg)
C
      OXYI=1.0-1.0/(1.0+EXP(1.0*(-COXYS(L,NY,NX)+2.5)))
      ORGCL=AMIN1(1.0E+05*BKVL(L,NY,NX),ORGC(L,NY,NX))
C
C     TOTAL MINERAL NH4, NO3 AND PO4
C
C     ALLOCATE NH4, NO3, HPO4, H2PO4 TO NON-BAND AND BAND FRACTIONS
C
C     VLNH4,VLNO3,VLPO4=fraction of soil volume in NH4,NO3,PO4 
C        non-band
C     VLNHB,VLNOB,VLPOB=fraction of soil volume in NH4,NO3,PO4 
C        band
C
      FNH4S=VLNH4(L,NY,NX)
      FNHBS=VLNHB(L,NY,NX)
      FNO3S=VLNO3(L,NY,NX)
      FNO3B=VLNOB(L,NY,NX)
      FNO2S=VLNO3(L,NY,NX)
      FNO2B=VLNOB(L,NY,NX)
      FH1PS=VLPO4(L,NY,NX)
      FH1PB=VLPOB(L,NY,NX)
      FH2PS=VLPO4(L,NY,NX)
      FH2PB=VLPOB(L,NY,NX)
C
C     CO2 CONSTRAINT ON CO2 UPTAKE BY AUTOTROPHS
C
C     CCO2S=aqueous CO2 concentration (g C m-3)
C     CCKM=Km for CO2 uptake (g C m-3)
C     XCO2=aqueous CO2 limitation to CO2 reduction
C
      XCO2=CCO2S(L,NY,NX)/(CCO2S(L,NY,NX)+CCKM)
C
C     TOTAL SUBSTRATE
C
C     TOSC=total SOC (g C)
C     TOSA=total colonized SOC (g C)
C     TORC=total microbial residue in each K (g C)
C     TOHC=total adsorbed C in each K (g C) 
C     OSCT=total SOC in each K (g C)
C     OSAT=total colonized SOC in each K (g C)  
C     ORCT=total microbial residue in each K (g C)
C     OHCT=total adsorbed C in each K (g C) 
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
C     WRITE(*,876)'ORCT',I,J,NFZ,NX,NY,L,K,M,ORCT(K)
C    2,ORC(M,K,L,NY,NX)
876   FORMAT(A8,8I4,60E12.4)
C     ENDIF
875   CONTINUE
      TORC=TORC+ORCT(K)
C
C     TOTAL ADSORBED AND DISSOLVED SUBSTRATE
C
C     TOHC=total adsorbed C (g C)
C     OSCH=total SOC in each K (g C) 
C     OSAH=total colonized SOC in each K (g C)
C     TSRH=total colonized SOC+microbial litter+adsorbed C (g C)
C
      TOHC=TOHC+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
880   CONTINUE
      DO 860 K=0,KL
      OSCH(K)=OSCT(K)+ORCT(K)+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
      OSAH(K)=OSAT(K)+ORCT(K)+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
C     IF((I/10)*10.EQ.I.AND.J.EQ.15.AND.NFZ.EQ.1)THEN
C     WRITE(*,861)'OSAH',I,J,NFZ,NX,NY,L,K,OSAH(K),OSCH(K)
C    2,OSAT(K),ORCT(K),OHC(K,L,NY,NX),OHA(K,L,NY,NX)
861   FORMAT(A8,7I4,20E12.4)
C     ENDIF      
860   CONTINUE
      TSRH=TOSA+TORC+TOHC
C
C     C:N AND C:P RATIOS OF TOTAL BIOMASS
C
C     OMC,OMN,OMP=microbial C,N,P biomass (g C)
C     CNOMA,CPOMA=N,P concentrations of active biomass (g N,P g C-1)
C     FCN,FCP,FCNP=effects of N,P limitations on biomass activity 
C
      TOMA=0.0
      TOMN=0.0
      DO 890 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 895 N=1,7
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      IF(OMC(1,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNOMA(N,K)=AMAX1(0.0,OMN(1,N,K,L,NY,NX)/OMC(1,N,K,L,NY,NX))
      CPOMA(N,K)=AMAX1(0.0,OMP(1,N,K,L,NY,NX)/OMC(1,N,K,L,NY,NX))
      ELSE
      CNOMA(N,K)=CNOMC(1,N,K)
      CPOMA(N,K)=CPOMC(1,N,K)
      ENDIF
      OMA(N,K)=AMAX1(0.0,OMC(1,N,K,L,NY,NX)/FL(1))
      FCN(N,K)=AMIN1(1.0,AMAX1(0.50,SQRT(CNOMA(N,K)/CNOMC(1,N,K))))
      FCP(N,K)=AMIN1(1.0,AMAX1(0.50,SQRT(CPOMA(N,K)/CPOMC(1,N,K))))
      FCNP(N,K)=AMIN1(FCN(N,K),FCP(N,K))
C
C     TOTAL BIOMASS
C
C     TOMA=total active microbial biomass (g C) 
C     OMA=active microbial biomass (g C)
C     OMC2,OMN2=active biomass C,N in recalcitrant fraction (g)
C     FL=allocation to labile(1),recalcitrant(2) fractions
C     TOMK,TONK,TOPK=total active biomass C,N,P (g)
C     TONX,TOPX=maximum total active biomass N,P (g)  
C     CNOMC,CPOMC=maximum N:C and P:C ratios in microbial biomass
C        from ‘starts.f’ (g N,P g C-1) 
C
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      TOMA=TOMA+OMA(N,K)
      ENDIF
      IF((K.LE.4.AND.N.EQ.2).OR.(K.EQ.5.AND.N.EQ.1))THEN
      TOMN=TOMN+OMA(N,K)
      ENDIF
      OMC2(N,K)=AMAX1(0.0,AMIN1(OMA(N,K)*FL(2),OMC(2,N,K,L,NY,NX)))
      IF(OMC(2,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOM2(N,K)=AMAX1(0.0,OMC2(N,K)/OMC(2,N,K,L,NY,NX))
      OMN2(N,K)=AMAX1(0.0,FOM2(N,K)*OMN(2,N,K,L,NY,NX))
      ELSE
      FOM2(N,K)=0.0
      OMN2(N,K)=0.0
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
      TOMK(K)=TOMK(K)+OMA(N,K)
      TONK(K)=TONK(K)+OMA(N,K)*CNOMA(N,K)
      TOPK(K)=TOPK(K)+OMA(N,K)*CPOMA(N,K)
      TONX(K)=TONX(K)+OMA(N,K)*CNOMC(1,N,K)
      TOPX(K)=TOPX(K)+OMA(N,K)*CPOMC(1,N,K)
685   CONTINUE
690   CONTINUE
C
C     TSRH=total colonized SOC+microbial litter+adsorbed C (g C)
C     FOSAH=fraction of TSRH in each substrate complex K
C     OSAH=total colonized SOC in each K (g C)
C
      DO 790 K=0,KL
      IF(TSRH.GT.ZEROS(NY,NX))THEN
      FOSAH(K,L,NY,NX)=OSAH(K)/TSRH
      ELSE
      FOSAH(K,L,NY,NX)=1.0
      ENDIF
C
C     DOC CONCENTRATIONS
C
C     VOLWM=soil water volume (m3)
C     FOSAH=fraction of TSRH in each substrate complex K
C     TSRH=total colonized SOC+microbial litter+adsorbed C (g C)
C     COQC,COQA=aqueous DOC,acetate concentrations (g C m-3)
C     OQC,OQA=DOC,acetate mass (g C)
C
      IF(VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IF(FOSAH(K,L,NY,NX).GT.ZERO)THEN
      COQC(K,L,NY,NX)=AMAX1(0.0,OQC(K,L,NY,NX)
     2/(VOLWM(NPH,L,NY,NX)*FOSAH(K,L,NY,NX)))
      COQA(K,L,NY,NX)=AMAX1(0.0,OQA(K,L,NY,NX)
     2/(VOLWM(NPH,L,NY,NX)*FOSAH(K,L,NY,NX)))
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
C     DOC RATIOS
C
C     CNQ,CPQ=DON:DOC,DOP:DOC
C     FOCA,FOAA=DOC,DOA:(DOC+DOA)
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
C     INPUTS FOR FERMENTATION, METHANOGENESIS
C
C     CHY1=aqueous H+ concentration (mol m-3)
C     CHNO2=nitrous acid concentration (g N m-3) 
C     GH2X=CH2GS feedback effect on energy yield of fermentation,
C        hydrogenotrophic methanogenesis (kJ mol-1)
C     TKS=soil temperature (K)
C     CH2GS=aqueous H2 concentration (g H m-3)
C     H2KI=H2 product inhibition for methanogenesis (g H m-3)
C 
      CHY1=AMAX1(ZERO,10.0**(-(PH(L,NY,NX)-3.0)))
      CHNO2=CNO2S(L,NY,NX)*CHY1/0.5
      CHNOB=CNO2B(L,NY,NX)*CHY1/0.5
      GH2X=8.3143E-03*TKS(L,NY,NX)
     2*LOG((AMAX1(1.0E-03,CH2GS(L,NY,NX))/H2KI)**4)
C
C     INITIALIZE TOTAL SUBSTRATE UPTAKE BY MICROBIAL POPULATIONS
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
C
C     TEMPERATURE AND WATER CONSTRINTS ON MICROBIAL ACTIVITY
C
      IF(L.NE.0)THEN
      LL=L
      ELSE
      LL=NU(NY,NX)
      ENDIF
      DO 760 K=0,5
C
C     WFNG=water potential effect on microbial respiration
C     PSISG=soil matric + osmotic potential (MPa)
C     TFNG=temperature+water limitation to biological activity
C     TFNR=temperature effect on maintenance respiration
C
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 750 N=1,7
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      IF(K.LE.4)THEN
      IF(N.EQ.3)THEN
      WFNG=EXP(0.10*PSISG)
      ELSE
      WFNG=EXP(0.20*PSISG)
      ENDIF
      ELSE
      WFNG=EXP(0.20*PSISG)
      ENDIF
      TFNG(N,K)=TFNX*WFNG
      TFNR(N,K)=TFNY
C
C     FOMA,FOMN=fraction of total active biomass C in each N and K
C     FOMK=fraction of total active biomass C in each N 
C
      IF(OMA(N,K).GT.0.0)THEN
      IF(TOMA.GT.ZEROS(NY,NX))THEN
      FOMA(N,K)=OMA(N,K)/TOMA
      ELSE
      FOMA(N,K)=1.0
      ENDIF
      IF(TOMN.GT.ZEROS(NY,NX))THEN
      FOMN(N,K)=OMA(N,K)/TOMN
      ELSE
      FOMN(N,K)=1.0
      ENDIF
      IF(TOMK(K).GT.ZEROS(NY,NX))THEN
      FOMK(N,K)=OMA(N,K)/TOMK(K)
      ELSE
      FOMK(N,K)=1.0
      ENDIF
C
C     ADJUST MICROBIAL RESPIRATION AND DECOMPOSITION RATES FOR BIOMASS
C
C     ORGCL=SOC used to calculate microbial concentration (g C)
C     COMC=microbial C concentration relative to substrate (g Mg-1) 
C     SPOMK=effect of COMC on microbial decay
C     RMOMK=effect of COMC on maintenance respiration
C     COMKI,COMKM=inhibition constant for microbial decomposition,
C        maintenance (g C Mg-1)
C
      IF(ORGCL.GT.ZEROS(NY,NX))THEN
      DO 765 M=1,2
      COMC=OMC(M,N,K,L,NY,NX)/ORGCL
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
C     R*Y=total substrate demand from R*X by all microbial, root and
C        mycorrhizal populations in ‘redist.f’ (g O,N,P,C t-1)
C     F*X=fraction of substrate uptake by each microbial population
C        relative to total demand 
C     substrate code: OXY=O2, NH4=NH4 non-band, NB4=NH4 band
C        NO3=NO3 non-band, NB3=NO3 band, PO4=H2PO4 non-band
C        POB=H2PO4 band,P14=HPO4 non-band, P1B=HPO4 band, OQC=DOC
C        oxidation, OQA=acetate oxidation  
C
      IF(ROXYY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOXYX=AMAX1(FMN,ROXYS(N,K,L,NY,NX)/ROXYY(L,NY,NX))
      ELSE
      FOXYX=AMAX1(FMN,FOMA(N,K))
      ENDIF
      IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4X=AMAX1(FMN,RINHO(N,K,L,NY,NX)/RNH4Y(L,NY,NX))
      ELSE
      FNH4X=AMAX1(FMN,FOMA(N,K)*VLNH4(L,NY,NX))
      ENDIF
      IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB4X=AMAX1(FMN,RINHB(N,K,L,NY,NX)/RNHBY(L,NY,NX))
      ELSE
      FNB4X=AMAX1(FMN,FOMA(N,K)*VLNHB(L,NY,NX))
      ENDIF
      IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3X=AMAX1(FMN,RINOO(N,K,L,NY,NX)/RNO3Y(L,NY,NX))
      ELSE
      FNO3X=AMAX1(FMN,FOMA(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB3X=AMAX1(FMN,RINOB(N,K,L,NY,NX)/RN3BY(L,NY,NX))
      ELSE
      FNB3X=AMAX1(FMN,FOMA(N,K)*VLNOB(L,NY,NX))
      ENDIF
      IF(RPO4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FPO4X=AMAX1(FMN,RIPOO(N,K,L,NY,NX)/RPO4Y(L,NY,NX))
      ELSE
      FPO4X=AMAX1(FMN,FOMA(N,K)*VLPO4(L,NY,NX))
      ENDIF
      IF(RPOBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FPOBX=AMAX1(FMN,RIPBO(N,K,L,NY,NX)/RPOBY(L,NY,NX))
      ELSE
      FPOBX=AMAX1(FMN,FOMA(N,K)*VLPOB(L,NY,NX))
      ENDIF
      IF(RP14Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FP14X=AMAX1(FMN,RIPO1(N,K,L,NY,NX)/RP14Y(L,NY,NX))
      ELSE
      FP14X=AMAX1(FMN,FOMA(N,K)*VLPO4(L,NY,NX))
      ENDIF
      IF(RP1BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FP1BX=AMAX1(FMN,RIPB1(N,K,L,NY,NX)/RP1BY(L,NY,NX))
      ELSE
      FP1BX=AMAX1(FMN,FOMA(N,K)*VLPOB(L,NY,NX))
      ENDIF
      IF(K.LE.4)THEN
      IF(ROQCY(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOQC=AMAX1(FMN,ROQCS(N,K,L,NY,NX)/ROQCY(K,L,NY,NX))
      ELSE
      FOQC=AMAX1(FMN,FOMK(N,K))
      ENDIF
      TFOQC=TFOQC+FOQC
      IF(ROQAY(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      FOQA=AMAX1(FMN,ROQAS(N,K,L,NY,NX)/ROQAY(K,L,NY,NX))
      ELSE
      FOQA=AMAX1(FMN,FOMK(N,K))
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
C     MICROBIAL POPULATIONS IN SURFACE LITTER
C
C     R*Y=total substrate demand from R*X in ‘redist.f’ (g N,P t-1)
C     F*R=fraction of substrate uptake relative to total uptake from
C        previous hour in surface litter 
C        :substrate code: NH4=NH4,NO3=NO3,PO4=H2PO4,P14=HPO4  
C
      IF(L.EQ.0)THEN
      IF(RNH4Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4XR(N,K)=AMAX1(FMN,RINHOR(N,K,NY,NX)
     2/RNH4Y(NU(NY,NX),NY,NX))
      ELSE
      FNH4XR(N,K)=AMAX1(FMN,FOMK(N,K))
      ENDIF
      IF(RNO3Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3XR(N,K)=AMAX1(FMN,RINOOR(N,K,NY,NX)
     2/RNO3Y(NU(NY,NX),NY,NX))
      ELSE
      FNO3XR(N,K)=AMAX1(FMN,FOMK(N,K))
      ENDIF
      IF(RPO4Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FPO4XR(N,K)=AMAX1(FMN,RIPOOR(N,K,NY,NX)
     2/RPO4Y(NU(NY,NX),NY,NX))
      ELSE
      FPO4XR(N,K)=AMAX1(FMN,FOMK(N,K))
      ENDIF
      IF(RP14Y(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      FP14XR(N,K)=AMAX1(FMN,RIPO1R(N,K,NY,NX)
     2/RP14Y(NU(NY,NX),NY,NX))
      ELSE
      FP14XR(N,K)=AMAX1(FMN,FOMK(N,K))
      ENDIF
      ENDIF
      IF(L.EQ.NU(NY,NX).AND.K.NE.3.AND.K.NE.4
     2.AND.BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
      TFNH4X=TFNH4X+FNH4XR(N,K)
      TFNO3X=TFNO3X+FNO3XR(N,K)
      TFPO4X=TFPO4X+FPO4XR(N,K)
      TFP14X=TFP14X+FP14XR(N,K)
      ENDIF
C
C     HETEROTROPHIC BIOMASS RESPIRATION
C
      IF(K.LE.4)THEN
C
C     RESPIRATION BY HETEROTROPHIC AEROBES
C
C        N=1;obligate aerobes
C         =2:facultative anaerobes
C         =3:fungi
C         =6:N2 fixers
C
      IF(N.LE.3.OR.N.EQ.6)THEN
C
C     ENERGY YIELDS OF O2 REDOX REACTIONS
C
C     EO2X,EO2D,EO2G,ENFX=growth respiration requirement from O2
C        reduction calculated in PARAMETERS 
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
C     COQC,COQA=DOC,DOA concentration (g C m-3)
C     OQKM,OQKA=Km for DOC,DOA uptake by heterotrophs(g C m-3)
C     FOCA,FOAA=DOC,DOA fractions of DOC+DOA
C     RGOCY=DOC,DOA-unlimited DOC+DOA oxidation (g C t-1)
C     RGOCZ,RGOAZ=microbial limitation to DOC,DOA oxidation (g C t-1)
C     RGOCX,RGOAX=DOC,DOA supply limitation to DOC,DOA oxidation 
C        (g C t-1)
C     RGOCP,RGOAP,RGOMP=O2-unlimited DOC,DOA,DOC+DOA oxidation
C        (g C t-1)
C     FCNP=N,P limitation to biological activity
C     VMXO=microbial specific oxidation rate (g C g C-1 h-1) 
C     WFNG=water stress effect
C     OMA=active bacterial biomass (g C)
C     TFNX=temperature effect on biological activity
C     FOQC,FOQA=fraction of OQC,OQA biological demand by aerobic
C        heterotrophic population  
C     EO2Q,EO2A=microbial respiration requirement from DOC,DOA
C        oxidation (g C g C-1)
C     OQC,OQA=DOC,DOA mass (g C)
C     FGOCP,FGOAP=fraction of RGOMP that oxidizes DOC,DOA
C     ECHZ=growth respiration requirement from O2+DOC reduction 
C        (g C g C-1) 
C     EO2Q,EO2A=growth respiration requirement from O2,DOC reduction
C        (g C g C-1) 
C     ROXYM,ROXYP,ROXYS=O2 demand from DOC oxidation
C        unconstrained,constained by DOC (g O t-1)
C     ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation (g C t-1)
C     ROQCD=DOC,DOA-unlimited microbial respiration used to represent 
C        microbial activity in decomposition (g C t-1) 
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      FSBSTC=COQC(K,L,NY,NX)/(COQC(K,L,NY,NX)+OQKM)
      FSBSTA=COQA(K,L,NY,NX)/(COQA(K,L,NY,NX)+OQKA)
      FSBST=FOCA(K)*FSBSTC+FOAA(K)*FSBSTA
      RGOCY=AMAX1(0.0,FCNP(N,K)*VMXO*WFNG*OMA(N,K))*XNFH 
      RGOCZ=RGOCY*FSBSTC*FOCA(K)*TFNX 
      RGOAZ=RGOCY*FSBSTA*FOAA(K)*TFNX 
      RGOCX=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC*EO2Q*XNFH)
      RGOAX=AMAX1(0.0,OQA(K,L,NY,NX)*FOQA*EO2A*XNFH)
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
      ECHZ=EO2Q*FGOCP+EO2A*FGOAP
      ROXYM(N,K)=2.667*RGOMP
      ROXYP(N,K)=ROXYM(N,K)
      ROXYSX=ROXYS(N,K,L,NY,NX)
      ROQCSX=ROQCS(N,K,L,NY,NX)
      ROQASX=ROQAS(N,K,L,NY,NX)
      ROXYS(N,K,L,NY,NX)=ROXYP(N,K)
      ROQCS(N,K,L,NY,NX)=RGOCZ 
      ROQAS(N,K,L,NY,NX)=RGOAZ 
      ROQCD(N,K)=RGOCY 
C     IF(I.EQ.146.AND.L.EQ.1)THEN
C     WRITE(*,5555)'RGOMP',I,J,NFZ,NX,NY,L,K,N,RGOMP,RGOCX,RGOAX,RGOCZ 
C    2,RGOAZ,RGOCX,RGOAX,FCNP(N,K),TFNG(N,K),VMXO,OMA(N,K),OSAH(K)
C    2,FOQC,FOQA,COQC(K,L,NY,NX),OQC(K,L,NY,NX),EO2Q,TKS(L,NY,NX)
C    3,COXYS(L,NY,NX),OQKM,OMC(1,N,K,L,NY,NX),OMC(2,N,K,L,NY,NX)
C    4,OMC(3,N,K,L,NY,NX),VOLWM(NPH,L,NY,NX),FOSAH(K,L,NY,NX)
C    5,FSBST,SPOMK(1),RMOMK(1),ROQCD(N,K),ROXYSX,ROXYS(N,K,L,NY,NX)
C    6,ROQCSX,ROQCS(N,K,L,NY,NX),ROQASX,ROQAS(N,K,L,NY,NX)
C    7,TFNX,WFNG,PSISM(L,NY,NX),TKS(L,NY,NX),COQC(K,L,NY,NX),OXYI
C    8,OQC(K,L,NY,NX)
5555  FORMAT(A8,8I4,60E12.4)
C     ENDIF
C
C     RESPIRATION BY HETEROTROPHIC ANAEROBES:
C     N=(4)ACETOGENIC FERMENTERS (7) ACETOGENIC N2 FIXERS
C
C     ENERGY YIELD FROM FERMENTATION DEPENDS ON H2 AND 
C     ACETATE CONCENTRATION
C
C     GH2X,GH2F=CH2GS feedback effect on energy yield of fermentation 
C        (kJ mol-1,kJ g C-1) 
C     GOAX,GOAF=COQA effect on energy yield of fermentation 
C        (kJ mol-1,kJ g C-1) 
C     TKS=soil temperature (K)
C     GHAX=CH2GS+COQA effect on energy yield of fermentation (kJ g C-1)
C     GCHX=energy yield of fermentation (kJ g C-1) 
C     ECHZ=growth respiration requirement of fermentation (g C g C-1)
C     EOMF,EOMN=energy requirement for fermenter,anaerobic diazotroph
C        growth (kJ g-1 C)
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
C     RESPIRATION RATES BY HETEROTROPHIC FERMENTERS 'RGOMP' FROM 
C     SPECIFIC OXIDATION RATE, ACTIVE BIOMASS, DOC CONCENTRATION,
C     MICROBIAL C:N:P FACTOR, AND TEMPERATURE FOLLOWED BY POTENTIAL 
C     RESPIRATION RATES 'RGOMP' WITH UNLIMITED SUBSTRATE USED FOR 
C     MICROBIAL COMPETITION FACTOR
C
C     COQC=aqueous DOC concentration (g C m-3)
C     OQKM=Km for DOC uptake by heterotrophs (g C m-3)
C     OXYI=O2 inhibition of fermentation (g O m-3)
C     RGOFY=DOC-unlimited DOC oxidation (g C t-1)
C     RGOFZ=microbial limitation to DOC oxidation (g C t-1)
C     RGOFX=DOC supply limitation to DOC oxidation (g C t-1)
C     RGOMP=O2-unlimited DOC oxidation (g C t-1)
C     FCNP=N,P limitation to oxidation 
C     VMXF=fermenter specific oxidation rate (g C g C-1 h-1)
C     WFNG=water stress effect on oxidation 
C     OMA=active fermenter biomass (g C) 
C     TFNX=temperature effect on biological activity
C     FOQC=fraction of OQC biological demand by fermenters
C     ECHZ=fermenter respiration requirement from DOC
C        oxidation (g C g C-1)
C     ROXYM,ROXYP,ROXYS=O2 demand from DOC oxidation
C        unconstrained,constained by DOC (g O t-1)
C     ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation (g C t-1)
C     ROQCD=substrate-unlimited microbial respiration used to represent 
C        microbial activity in decomposition (g C t-1) 
C     TRH2G=total H2 uptake (g H t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      FSBST=COQC(K,L,NY,NX)/(COQC(K,L,NY,NX)+OQKM)*OXYI
      RGOFY=AMAX1(0.0,FCNP(N,K)*VMXF*WFNG*OMA(N,K)*XNFH) 
      RGOFZ=RGOFY*FSBST*TFNX 
      RGOFX=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC*ECHZ*XNFH)
      RGOMP=AMIN1(RGOFX,RGOFZ)
      FGOCP=1.0
      FGOAP=0.0
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      ROQCS(N,K,L,NY,NX)=RGOFZ 
      ROQAS(N,K,L,NY,NX)=0.0 
      ROQCD(N,K)=RGOFY 
      TRH2G=TRH2G+RGOMP
C     IF(I.EQ.322)THEN
C     WRITE(*,5554)'FERM',I,J,NFZ,NX,NY,L,K,N
C    2,RGOMP,RGOFX,RGOFZ,GHAX,GOAF 
C    2,ECHZ,FCNP(N,K),TFNG(N,K),OMA(N,K),OSAH(K),FOQC,COQC(K,L,NY,NX)
C    3,OQKM,OMC(1,N,K,L,NY,NX),OMC(2,N,K,L,NY,NX),OMC(3,N,K,L,NY,NX)
C    3,OMN(1,N,K,L,NY,NX),OMN(2,N,K,L,NY,NX),OMN(3,N,K,L,NY,NX)
C    5,VOLWM(NPH,L,NY,NX),PSISM(L,NY,NX),WFNG,COXYS(L,NY,NX),OXYI
C    6,FSBST,FOSAH(K,L,NY,NX),SPOMK(1),RMOMK(1),ROQCD(N,K)
C    7,OXYS(L,NY,NX),VOLW(L,NY,NX)
5554  FORMAT(A8,8I4,60E12.4)
C     ENDIF
C
C     ENERGY YIELD FROM ACETOTROPHIC METHANOGENESIS
C
C     GOMX,GOMM=acetate feedback effect on methanogenesis energy yield 
C        (kJ mol,kJ g C-1)
C     TKS=soil temperature (K)
C     COQA=DOA concentration (g C m-3)
C     OAKI=acetate product inhibition for acetotrophic 
C        methanogenesis (g C m-3) 
C     ECHZ=growth respiration requirement of acetotrophic
C        methanogenesis (g C g C-1)
C     GC4X=energy yield of acetotrophic methanogenesis (kJ g C-1)
C     EOMH=growth respiration requirement (g C g C-1)
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
C     COQA=DOA concentration (g C m-3)
C     OQKAM=Km for acetate uptake (g C m-3)
C     FCNP=N,P limitation to biological activity
C     VMXM=specific respiration rate of acetotrophic methanogens 
C        (g C g C-1 h-1)
C     WFNG=water stress effect
C     OMA=active methanogenic biomass (g C)
C     TFNX=temperature effect on biological activity
C     FOQA=DOA supply limitation
C     RGOGY=DOA-unlimited DOA oxidation (g C t-1)
C     RGOGZ=microbial limitation to DOA oxidation (g C t-1)
C     RGOGX=DOA supply limitation to DOA oxidation (g C t-1)
C     RGOMP=O2-unlimited DOA oxidation (g C t-1) 
C     OQA=acetate DOA (g)
C     FOQA=fraction of biological demand for acetate
C     ROXYM,ROXYP,ROXYS=O2 demand from DOA oxidation
C        unconstrained,constained by DOC (g O t-1)
C     ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation (g C t-1)
C     ROQCD=DOA-unlimited microbial respiration used to represent 
C        microbial activity in decomposition (g C t-1) 
C     TCH4H=total heterotrophic CH4 production (g C t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      FSBST=COQA(K,L,NY,NX)/(COQA(K,L,NY,NX)+OQKAM)
      RGOGY=AMAX1(0.0,FCNP(N,K)*VMXM*WFNG*OMA(N,K)*XNFH)
      RGOGZ=RGOGY*FSBST*TFNX 
      RGOGX=AMAX1(0.0,OQA(K,L,NY,NX)*FOQA*ECHZ*XNFH) 
      RGOMP=AMIN1(RGOGX,RGOGZ)
      FGOCP=0.0
      FGOAP=1.0
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      ROQCS(N,K,L,NY,NX)=0.0 
      ROQAS(N,K,L,NY,NX)=RGOGZ 
      ROQCD(N,K)=0.0 
      TCH4H=TCH4H+0.5*RGOMP
C     IF((I/30)*30.EQ.I.AND.NX.EQ.3.AND.NY.EQ.1.AND.J.EQ.24)THEN
C     WRITE(*,5552)'ACMETH',I,J,NX,NY,L,K,N,RGOMP,RGOGZ,RGOGX,GOMM 
C    2,ECHZ,FCNP(N,K),TFNG(N,K),OMA(N,K),FOQA,COQA(K,L,NY,NX)
C    2,OQA(K,L,NY,NX)
C    3,OMC(1,N,K,L,NY,NX),OMC(2,N,K,L,NY,NX),OMC(3,N,K,L,NY,NX)
C    3,OMN(1,N,K,L,NY,NX),OMN(2,N,K,L,NY,NX),OMN(3,N,K,L,NY,NX)
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
C     MICROBIAL COMPETITION FACTOR. 
C        N=1:NH3 oxidizers 
C         =2:NO2 oxidizers
C         =3:CH4 oxidizers
C         =5:H2trophic methanogens
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
C     FNH4,FNB4=fraction of total biological demand for NH4 
C        in non-band, band by NH3 oxidizers 
C     RNH4Y,RVMX4=total,NH3 oxidizer NH4 demand in non-band from 
C        previous time step (g N t-1)
C     RNHBY,RVMB4=total,NH3 oxidizer NH4 demand in band from 
C        previous time step (g N t-1)
C
      IF(RNH4Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNH4=AMAX1(FMN,RVMX4(N,K,L,NY,NX)/RNH4Y(L,NY,NX))
      ELSE
      FNH4=AMAX1(FMN,VLNH4(L,NY,NX)*FOMA(N,K))
      ENDIF
      IF(RNHBY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB4=AMAX1(FMN,RVMB4(N,K,L,NY,NX)/RNHBY(L,NY,NX))
      ELSE
      FNB4=AMAX1(FMN,VLNHB(L,NY,NX)*FOMA(N,K))
      ENDIF
      TFNH4X=TFNH4X+FNH4
      TFNH4B=TFNH4B+FNB4
C
C     NITRIFICATION INHIBITION
C
C     ZNFN0,ZNFNI=initial (from ‘hour1.f’),current 
C        nitrification inhibition activity 
C     TFNX=temperature effect on biological activity
C     ZNFN4S,ZNFN4B=inhibition in non-band, band
C     CNH4S,CNH4B=NH4 concentrations in non-band, band (g N m-3)
C     RNFNI=rate constant for decline in nitrification inhibition (h-1)
C     ZHKI=inhibition of decline in ZNFNI from high CNH4 (g N m-3) 
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      IF(ZNFN0(L,NY,NX).GT.ZEROS(NY,NX))THEN
      ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)*(1.0-RNFNI*TFNX*XNFH)
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
C     ECHZ=growth respiration requirement of NH3 oxidation
C        (g C g C-1)
C     VMXX=NH4 unlimited NH3 oxidation (g N t-1)
C     VMXH=specific oxidation (g N g C-1 h-1)
C     TFNG=temperature+water limitation to biological activity
C     FCNP=N,P limitation to biological activity
C     XCO2=aqueous CO2 limitation to CO2 reduction
C     OMA=active biomass (g C)
C     VMXA=NH4 unlimited NH3 oxidation (g N t-1)
C     VHKI=nonlinear increase in VMXA with VMXX (g N t-1)
C     FNH4S,FNHBS=fractions of NH4 in non-band, band
C     CNH4S,CNH4B=NH4 concentration in non-band, band (g N m-3)
C     ZHKM=Km for NH4 uptake (g N m-3)
C     VMX4S,VMX4B=NH3 oxidation in non-band, band unlimited by NH4
C        supply (g N t-1)
C     RNNH4,RNNHB=NH3 oxidation in non-band, band limited by NH4
C        supply (g N t-1)
C     FNH4,FNB4=fractions of total NH4 demand in non-band, band
C        by NH3 oxidizers
C     ZNH4S,ZNH4B=NH4 amount in non-band, band (g N)
C     ZNFN4S,ZNFN4B=inhibition in non-band, band
C     RVOXP=total NH3 oxidation in non-band+band (g N t-1)
C     RGOMP=O2-unlimited respiration (g C t-1)
C     ECNH=efficiency of CO2 conversion to biomass (g C g N-1)
C     RVMX4,RVMXB=NH3 oxidation in non-band, band unlimited by NH4
C        supply (g N t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      ECHZ=EO2X
      VMXX=VMXH*TFNG(N,K)*FCNP(N,K)*XCO2*OMA(N,K)*XNFH 
      IF(VOLWY.GT.ZEROS2(NY,NX))THEN
      VMXA=VMXX/(1.0+VMXX/(VHKI*VOLWY)) 
      ELSE
      VMXA=0.0
      ENDIF
      FCN4S=FNH4S*CNH4S(L,NY,NX)/(CNH4S(L,NY,NX)+ZHKM) 
      FCN4B=FNHBS*CNH4B(L,NY,NX)/(CNH4B(L,NY,NX)+ZHKM) 
      FSBST=FCN4S+FCN4B
      VMX4S=VMXA*FCN4S 
      VMX4B=VMXA*FCN4B 
      RNNH4=AMAX1(0.0,AMIN1(VMX4S,FNH4*ZNH4S(L,NY,NX)*XNFH))*ZNFN4S 
      RNNHB=AMAX1(0.0,AMIN1(VMX4B,FNB4*ZNH4B(L,NY,NX)*XNFH))*ZNFN4B
      RVOXP=RNNH4+RNNHB
      RVOXPA=RNNH4
      RVOXPB=RNNHB
      RGOMP=AMAX1(0.0,RVOXP*ECNH*ECHZ) 
      RVMX4(N,K,L,NY,NX)=VMX4S 
      RVMB4(N,K,L,NY,NX)=VMX4B 
C
C     O2 DEMAND FROM NH3 OXIDATION
C
C     ROXYM=O2 demand from respiration (g O t-1) 
C     ROXYP=O2 demand from respiration + NH3 oxidation (g O t-1) 
C
      ROXYM(N,K)=2.667*RGOMP
      ROXYP(N,K)=ROXYM(N,K)+3.429*RVOXP 
      ROXYS(N,K,L,NY,NX)=ROXYP(N,K)
C     IF(IYRC.EQ.2012.AND.I.EQ.151.AND.NX.EQ.1)THEN
C     WRITE(*,6666)'NITRI',I,J,NFZ,L,K,N
C    2,RNNH4,RNNHB,VMXX,VMXA,VOLWY 
C    2,CNH4S(L,NY,NX),CNH4B(L,NY,NX),CNH3S(L,NY,NX),CNH3B(L,NY,NX)
C    3,ZNH4S(L,NY,NX),ZNH4B(L,NY,NX),ZNH3S(L,NY,NX),ZNH3B(L,NY,NX)
C    2,14.0*XN4(L,NY,NX),14.0*XNB(L,NY,NX)
C    3,ZNHUFA(L,NY,NX),ZNHUFB(L,NY,NX),VLNH4(L,NY,NX),VLNHB(L,NY,NX)
C    4,COXYS(L,NY,NX),RGOMP,PH(L,NY,NX),TFNX,FCNP(N,K),XCO2,ROXYM(N,K)
C    5,VMX4S,VMX4B,FCN4S,FCN4B,FNH4S,FNHBS,OMA(N,K)
C    6,FNH4,FNB4,ZNFN4S,ZNFN4B,ZNFNI(L,NY,NX),ZNFN0(L,NY,NX)
6666  FORMAT(A8,6I4,40E12.4)
C     ENDIF
C
C     NO2 OXIDIZERS
C
      ELSEIF(N.EQ.2)THEN
C
C     FACTOR TO REGULATE COMPETITION FOR NO2 AMONG DIFFERENT
C     MICROBIAL POPULATIONS
C
C     FNO2,FNB2=fraction of total biological demand for NO2 
C        in non-band, band by NO2 oxidizers 
C     RNO2Y,RVMX2=total,NO2 oxidizer NO2 demand in non-band from 
C        previous time step
C     RN2BY,RVMB2=total,NO2 oxidizer NO2 demand in band from 
C        previous time step
C
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMX2(N,K,L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=AMAX1(FMN,FOMN(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMB2(N,K,L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=AMAX1(FMN,FOMN(N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
C
C     NO2 OXIDATION FROM SPECIFIC OXIDATION RATE, ENERGY YIELD,
C     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND
C     NO2 CONCENTRATIONS
C
C     ECHZ=growth respiration requirement of NO2 oxidation
C        (g C g C-1)
C     VMXA=NO2 unlimited NO2 oxidation (g N t-1)
C     VMXN=specific nitrifier oxidation (g N g C-1 h-1)
C     TFNG=temperature+water limitation to biological activity
C     FCNP=N,P limitation to biological activity
C     XCO2=aqueous CO2 limitation to CO2 reduction
C     OMA=active biomass (g)
C     FNH4S,FNHBS=fractions of NH4 in non-band, band
C     CNO2S,CNO2B=NO2 concentration in non-band, band (g N m-3)
C     ZNKM=Km for NO2 uptake (g N m-3)
C     VMX2S,VMX2B=NO2 oxidation in non-band, band unlimited by NO2
C        supply (g N t-1)
C     FNO2,FNB2=fractions of total NO2 demand in non-band, band
C        by NO2 oxidizers
C     ZNO2S,ZNO2B=NO2 amount in non-band, band (g N)
C     RVOXP=total NO2 oxidation in non-band+band (g N t-1)
C     RNNO2,RNNOB=NO2 oxidation in non-band, band limited by NO2
C        supply (g N t-1)
C     RGOMP=O2-unlimited respiration (g C t-1)
C     ECNO=efficiency CO2 conversion to biomass (g C g N-1)
C     RVMX2,RVMB2=NO2 oxidation in non-band, band unlimited by NO2
C        supply (g N t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      ECHZ=EO2X
      VMXA=TFNG(N,K)*FCNP(N,K)*XCO2*OMA(N,K)*VMXN*XNFH 
      FCN2S=FNH4S*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+ZNKM) 
      FCN2B=FNHBS*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+ZNKM) 
      FSBST=FCN2S+FCN2B
      VMX2S=VMXA*FCN2S
      VMX2B=VMXA*FCN2B
      RNNO2=AMAX1(0.0,AMIN1(VMX2S,FNO2*ZNO2S(L,NY,NX)*XNFH))
      RNNOB=AMAX1(0.0,AMIN1(VMX2B,FNB2*ZNO2B(L,NY,NX)*XNFH))
      RVOXP=RNNO2+RNNOB
      RVOXPA=RNNO2
      RVOXPB=RNNOB
      RGOMP=AMAX1(0.0,RVOXP*ECNO*ECHZ) 
      RVMX2(N,K,L,NY,NX)=VMX2S 
      RVMB2(N,K,L,NY,NX)=VMX2B
C
C     O2 DEMAND FROM NO2 OXIDATION
C
C     ROXYM=O2 demand from respiration by nitrifiers (g O t-1)
C     ROXYP=O2 demand from respiration + NO2 oxidation (g O t-1) 
C
      ROXYM(N,K)=2.667*RGOMP
      ROXYP(N,K)=ROXYM(N,K)+1.143*RVOXP 
      ROXYS(N,K,L,NY,NX)=ROXYP(N,K)
C     IF(L.EQ.NU(NY,NX))THEN
C     WRITE(*,6667)'NO2OX',I,J,NFZ,L,K,N,RNNO2,RNNOB,ZNO2S(L,NY,NX)
C    2,ZNO2B(L,NY,NX),CNO2S(L,NY,NX),CNO2B(L,NY,NX),CNH3S(L,NY,NX)
C    3,CNH3B(L,NY,NX),CNH4S(L,NY,NX),CNH4B(L,NY,NX),CNO3S(L,NY,NX) 
C    3,CNO3B(L,NY,NX),CHNO2,CHNOB,VMXA,TFNG(N,K),FCNP(N,K),VMXN,ZNKM
C    4,FCN2S,FCN2B,OMA(N,K),FOMN(N,K),TOMN,RVMX2(N,K,L,NY,NX)
C    5,RNO2Y(L,NY,NX),FNO2,FNB2,ROXYM(N,K),ROXYP(N,K) 
C    6,ROXYS(N,K,L,NY,NX),VLNHB(L,NY,NX),VLNOB(L,NY,NX)
C    7,SPOMK(1),RMOMK(1),RGOMP 
6667  FORMAT(A8,6I4,50E12.4)
C     ENDIF
C
C     H2TROPHIC METHANOGENS
C
      ELSEIF(N.EQ.5)THEN
C
C     CO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE OXIDIZER BIOMASS, TEMPERATURE, AQUEOUS CO2 AND H2
C
C     GH2X,GH2H=CH2GS feedback effect on energy yield of
C        hydrogenotrophic methanogenesis (kJ mol-1,kJ g C-1)
C     CH2GS=aqueous H2 concentration (g H m-3)
C     ECHZ=growth respiration requirement of hydrogenotrophic
C        methanogenesis (g C g C-1)
C     GCOX=free energy yield of hydrogenotrophic
C        methanogenesis (kJ g-1 C)
C     EOMH=energy requirements for hydrogenotrophic
C        methanogen growth (kJ g-1 C)
C     VMXA=H2-unlimited CO2 reduction rate (g t-1)
C     TFNG=temperature+water limitation to biological activity
C     FCNP=N,P limitation to biological activity
C     XCO2=aqueous CO2 limitation to CO2 reduction
C     OMA=active biomass (g)
C     VMXC=specific hydrogenotrophic methanogenic CO2 reduction rate 
C        (g C g C-1 h-1) 
C     H2GSX=aqueous H2 (H2GS) + H2 from fermentation (TRH2G) (g H)
C     H2KM=Km for H2 uptake (g H m-3)
C     RGOMP=H2-limited CO2 reduction rate (g C t-1)
C     ROXYM,ROXYS=O2 demand (=0)
C     TCH4A=total autotrophic CH4 production (g C t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      GH2H=GH2X/12.0
      ECHZ=AMAX1(EO2X,AMIN1(1.0
     2,1.0/(1.0+AMAX1(0.0,(GCOX+GH2H))/EOMH)))
      VMXA=TFNG(N,K)*FCNP(N,K)*XCO2*OMA(N,K)*VMXC*XNFH 
      H2GSX=H2GS(L,NY,NX)+0.111*TRH2G
      FSBST=CH2GS(L,NY,NX)/(CH2GS(L,NY,NX)+H2KM)
      RGOMP=AMAX1(0.0,AMIN1(VMXA*FSBST,1.5*H2GSX*XNFH))
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      TCH4A=TCH4A+RGOMP 
C     IF((I/30)*30.EQ.I.AND.NX.EQ.3.AND.NY.EQ.1.AND.J.EQ.24)THEN
C     WRITE(*,5553)'H2METH',I,J,NFZ,NX,NY,L,K,N,RGOMP,H2GS(L,NY,NX)
C    2,H2GSX,CH2GS(L,NY,NX),VMXA,TFNG(N,K),FCNP(N,K),XCO2 
C    3,OMA(N,K),VMXC,ECHZ,GCOX,GH2H,TKS(L,NY,NX),FSBST
C    4,SPOMK(1),RMOMK(1) 
5553  FORMAT(A8,8I4,20E12.4)
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
C     ECHZ=growth respiration requirement from CH4 oxidation 
C        (g C g C-1)  
C     VMXA=CH4-unlimited CH4 oxidation
C     TFNG=temperature+water limitation to biological activity
C     FCNP=N,P limitation to biological activity
C     OMA=active biomass (g)
C     VMX4=specific methanotroph respiration rate (g C g C-1 h-1)
C     RCH4L,RCH4F=total aqueous,gaseous CH4 exchange from previous 
C        time step (g C t-1)
C     TCH4H+TCH4A=total CH4 production from methanogenesis (g C t-1)
C     CH4G1,CH4S1=gaseous, aqueous CH4 (g C) 
C     CCH4E,CCH4G=CH4 gas concentration in atmosphere, soil (g C m-3)
C     XNPG,XNFH=time step for gas fluxes, biological activity 
C        from ‘wthr.f’ (h t-1)     
C
      ECHZ=EH4X
      VMXA=TFNG(N,K)*FCNP(N,K)*OMA(N,K)*VMX4*XNFH 
      VMXA1=VMXA*XNPG
      RCH4L1=RCH4L(L,NY,NX)*XNPG
      RCH4F1=RCH4F(L,NY,NX)*XNPG
      RCH4S1=(TCH4H+TCH4A)*XNPG
      IF(L.EQ.0)THEN
      CH4G1=CCH4E(NY,NX)*VOLPM(1,L,NY,NX)
      ELSE
      CH4G1=CCH4G(L,NY,NX)*VOLPM(1,L,NY,NX)
      ENDIF
      CH4S1=CH4S(L,NY,NX) 
      RVOXP=0.0
      RGOMP=0.0
C
C     CH4 DISSOLUTION FROM GASEOUS PHASE SOLVED IN SHORTER TIME STEP
C     TO MAINTAIN AQUEOUS CH4 CONCENTRATION DURING OXIDATION
C
C     VOLPM,VOLWM=air,water-filled porosity (m3)
C     SCH4L=CH4 aqueous solubility (g C m-3/(g C m-3))
C     CH4G1,CH4S1=gaseous, aqueous CH4 (g C) 
C     RCH4L1,RCH4F1=total aqueous,gaseous CH4 exchange from previous 
C        time step (g C t-1)
C     RCH4S1=total CH4 production from methanogenesis (g C t-1)
C     CCH4S1=aqueous CH4 concentration (g C m-3)
C     CCK4=Km for CH4 uptake (g C m-3)
C     VMXA1=CH4-unlimited CH4 oxidation
C     RVOXP1=CH4 oxidation rate limited by CH4 supply (g C t-1)
C     ECHO=efficiency of CO2 conversion to biomass (g C g C-1)
C     ECHZ=growth respiration requirement from CH4 oxidation 
C        (g C g C-1)  
C     RGOMP1=CH4-limited CH4 respiration (g C t-1)
C     RCHDF=gaseous-aqueous CH4 exchange (g C t-1)
C     DFGS=rate constant for gaseous-aqueous exchange from ‘watsub.f’
C        (t-1)
C     XNPG,XNFH=time step for gas fluxes, biological activity 
C        from ‘wthr.f’ (h t-1) 
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
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.1)THEN
C    2.AND.MM.EQ.NPT)THEN
C     WRITE(*,5547)'CH4OX',I,J,NFZ,M,MM,NX,NY,L,K,N
C    2,RVOXP1,RGOMP1,CH4G1 
C    2,CH4S1,VMXA1,RVOXP,RGOMP,RCHDF,RCH4L1,RCH4F1,RCH4S1,CCH4S1 
C    3,ECHO,ECHZ,OMA(N,K),VOLWM(M,L,NY,NX),VOLPM(M,L,NY,NX),VOLWCH
C    4,THETPM(M,L,NY,NX),SCH4L(L,NY,NX),DFGS(M,L,NY,NX)
C    5,COXYS(L,NY,NX),CCH4E(NY,NX),FSBST,SPOMK(1),RMOMK(1)
C    6,CH4S(L,NY,NX) 
5547  FORMAT(A8,10I4,30E12.4)
C     ENDIF
325   CONTINUE
      ENDIF
320   CONTINUE
      RVOXPA=RVOXP
      RVOXPB=0.0
C
C     O2 DEMAND FROM CH4 OXIDATION
C
C     ROXYM=O2 demand from CH4 respiration (g O t-1)
C     ROXYP=O2 demand from CH4 respiration + oxidation (g O t-1)
C
      ROXYM(N,K)=2.667*RGOMP
      ROXYP(N,K)=ROXYM(N,K)+5.333*RVOXP 
      ROXYS(N,K,L,NY,NX)=ROXYP(N,K) 
      ELSE
      RGOMP=0.0
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      ENDIF
      ELSE
      RGOMP=0.0
      ROXYM(N,K)=0.0
      ROXYP(N,K)=0.0
      ROXYS(N,K,L,NY,NX)=0.0
      ENDIF
C
C     O2 UPTAKE BY AEROBES
C
C     RUPOX,ROXYP=O2-limited,O2-unlimited O2 uptake rates (g O t-1)
C     RUPMX=O2-unlimited maximum O2 uptake rate (g O t-1)
C     FOXYX=fraction of O2 uptake by microbial population
C        relative to total
C     ROXYP=O2 demand from all oxidation reactions (g O t-1) 
C     OXYG1,OXYS1=gaseous,aqueous O2 mass (g O)
C     COXYG,COXYS=gaseous,aqueous O2 concentration (g O m-3)
C     ROXYF,ROXYL=net O2 gaseous,aqueous flux through soil layer
C        from ‘redist.f’ (g O t-1)
C     ROXYFX,ROXYLX=net O2 gaseous,aqueous flux through soil layer 
C        during gas flux time step (g O t-1)  
C     OLSGL=aqueous O2 diffusivity (m2 h-1)
C     FLQRQ,FLQRI=surface water flux from precipitation, irrigation 
C        (m3 t-1)
C     COXR,COXQ=O2 concentration in FLQRQ,FLQRI (g O m-3)
C     XNPG,XNFH=time step for gas fluxes, biological activity 
C        from ‘wthr.f’ (h t-1)     
C
      RUPOX(N,K)=0.0
      IF(N.LE.3.OR.N.EQ.6)THEN
      IF(ROXYP(N,K).GT.ZEROS(NY,NX).AND.FOXYX.GT.ZERO)THEN
      IF(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX))THEN
C
C     MAXIMUM O2 UPAKE FROM POTENTIAL RESPIRATION OF EACH AEROBIC
C     POPULATION
C
      RUPMX=ROXYP(N,K)*XNPG 
      ROXYFX=ROXYF(L,NY,NX)*XNPG*FOXYX
      OLSGL1=OLSGL(L,NY,NX)*XNPGX 
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
C     ACTIVE O2 UPTAKE
C
C     VOLWM,VOLPM,VOLY=water, air and total volumes (m3)
C     ORAD=microbial radius from ‘starts.f’(m)
C     FILM=water film thickness from ‘watsub.f’(m)
C     DIFOX=aqueous O2 diffusivity (m2 t-1)
C     TORT=diffusion tortuosity
C     OLSGL1=aqueous O2 diffusivity (m2 t-1)
C     BIOS=microbial number *g-1)
C     OMA=active biomass (g)
C     SOXYL=O2 solubility from ‘hour1.f’ (g O m-3/(g O m-3))
C     OXYG1,OXYS1=gaseous,aqueous O2 mass (g O)
C     COXYG1,COXYS1=gaseous,aqueous O2 concentration (g O m-3)
C     ROXYFX,ROXYLX=net O2 gaseous,aqueous flux through soil layer 
C        during time step (g O t-1) 
C     OXKM=Km for O2 uptake from ‘starts.f’(g O m-3)
C     RUPMX=O2-unlimited maximum O2 uptake rate (g O t-1)
C     RMPOX=O2-limited O2 uptake rate (g O t-1)
C
      RRADO=ORAD*(FILM(M,L,NY,NX)+ORAD)/FILM(M,L,NY,NX)
      DIFOX=TORT(M,L,NY,NX)*OLSGL1*12.57*BIOS*OMA(N,K)*RRADO
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
      B=-RUPMX-DIFOX*OXKM-X
      C=X*RUPMX
      RMPOX=(-B-SQRT(B*B-4.0*C))/2.0
      ELSE
      RMPOX=0.0
      ENDIF
      OXYS1=OXYS1-RMPOX
C
C     GASEOUS-AQUEOUS O2 EXCHANGE
C
C     THETPM,VOLWPM=air concentration (m3 m-3),volume (m3)
C        from ‘watsub.f’
C     ROXDFQ=gaseous-aqueous O2 exchange (g O t-1)
C     OXYG1,OXYS1=gaseous,aqueous O2 (g O) 
C     DFGS=rate constant for air-water gas exchange from ‘watsub.f’ 
C        (t-1)
C     RUPOX=accumulated O2 uptake (g O t-1)
C     RMPOX=O2-limited O2 uptake rate (g O t-1)
C     ROXSK=total soil O2 uptake from soil by all microbial,root
C        populations used in ‘trnsfr.f’ (g O t-1)
C    
      IF(THETPM(M,L,NY,NX).GT.THETX.AND.VOLPOX.GT.ZEROS(NY,NX))THEN
      ROXDFQ=DFGS(M,L,NY,NX)*(AMAX1(ZEROS(NY,NX),OXYG1)*VOLWOX
     2-OXYS1*VOLPOX)/VOLWPM
      ELSE
      ROXDFQ=0.0
      ENDIF
      OXYG1=OXYG1-ROXDFQ 
      OXYS1=OXYS1+ROXDFQ 
      RUPOX(N,K)=RUPOX(N,K)+RMPOX
      ROXSK(M,L,NY,NX)=ROXSK(M,L,NY,NX)+RMPOX
C     IF(I.EQ.150.AND.L.EQ.1
C    2.AND.K.EQ.1.AND.N.EQ.1)THEN
C     WRITE(*,5544)'OXY',I,J,NFZ,M,MX,L,K,N,COXYS1
C    2,OXYS1,RMPOX,ROXDFQ,ROXYLX,ROXYFX,FOXYX 
C    2,RUPOX(N,K),ROXYP(N,K),ROXSK(M,L,NY,NX)
C    2,RUPMX,DIFOX,OLSGL1,BIOS,OMA(N,K),X 
C    2,OXYG1/(VOLPM(M,L,NY,NX)*FOXYX)
C    5,THETPM(M,L,NY,NX),DFGS(M,L,NY,NX) 
C    6,VOLPM(M,L,NY,NX),VOLWM(M,L,NY,NX),VOLA(L,NY,NX)
C    7,COXYS(L,NY,NX),COXYG(L,NY,NX),ROXYY(L,NY,NX) 
5544  FORMAT(A8,8I4,50E12.4)
C     ENDIF
425   CONTINUE
420   CONTINUE
C
C     RATIO OF ACTUAL O2 UPAKE TO BIOLOGICAL DEMAND (WFN)
C
C     WFN=ratio of O2-limited to O2-unlimited uptake
C     RUPOX,ROXYP=O2-limited,O2-unlimited uptake (g O t-1) 
C     RVMX4,RVNHB,RVMX2,RVMB2=O2-limited NH3,NO2 oxidation 
C        in non-band,band (g N t-1)
C
      WFN(N,K)=AMIN1(1.0,AMAX1(0.0,RUPOX(N,K)/ROXYP(N,K)))
      IF(K.EQ.5)THEN
      IF(N.EQ.1)THEN
      RVMX4(N,K,L,NY,NX)=RVMX4(N,K,L,NY,NX)*WFN(N,K)
      RVMB4(N,K,L,NY,NX)=RVMB4(N,K,L,NY,NX)*WFN(N,K)
      ELSEIF(N.EQ.2)THEN
      RVMX2(N,K,L,NY,NX)=RVMX2(N,K,L,NY,NX)*WFN(N,K)
      RVMB2(N,K,L,NY,NX)=RVMB2(N,K,L,NY,NX)*WFN(N,K)
      ENDIF
      ENDIF
      ELSE
      RUPOX(N,K)=ROXYP(N,K)
      WFN(N,K)=1.0
      ENDIF
      ELSE
      RUPOX(N,K)=0.0
      WFN(N,K)=1.0
      ENDIF
C
C     RESPIRATION PRODUCTS ALLOCATED TO O2, CO2, ACETATE, CH4, H2
C
C     RGOMO,RGOMP=O2-limited, O2-unlimited respiration (g C t-1)
C     WFN=ratio of O2-limited to O2-unlimited uptake
C     RCO2X,RCH3X,RCH4X,RH2GX=CO2,acetate,CH4,H2 production from RGOMO
C        (g C t-1)
C     ROXYO=O2-limited O2 uptake (g O t-1)
C     RVOXA,RVOXB=O2-limited oxidation in non-band,band (g O t-1) of:
C        N=1:NH4
C        N=2:NO2
C        N=3:CH4 
C
      RGOMO(N,K)=RGOMP*WFN(N,K)
      RCO2X(N,K)=RGOMO(N,K)
      RCH3X(N,K)=0.0
      RCH4X(N,K)=0.0
      ROXYO(N,K)=ROXYM(N,K)*WFN(N,K)
      RH2GX(N,K)=0.0
      IF(K.EQ.5)THEN
      RVOXA(N)=RVOXPA*WFN(N,K)
      RVOXB(N)=RVOXPB*WFN(N,K)
      ENDIF
      ELSEIF(N.EQ.4.OR.N.EQ.7)THEN
      RGOMO(N,K)=RGOMP
      RCO2X(N,K)=0.333*RGOMO(N,K)
      RCH3X(N,K)=0.667*RGOMO(N,K)
      RCH4X(N,K)=0.0
      ROXYO(N,K)=ROXYM(N,K)
      IF(K.LE.4)THEN
      RH2GX(N,K)=0.111*RGOMO(N,K)
      ELSE
      RH2GX(N,K)=0.0
      ENDIF
      ELSEIF(N.EQ.5)THEN
      RGOMO(N,K)=RGOMP
      IF(K.LE.4)THEN
      RCO2X(N,K)=0.50*RGOMO(N,K)
      RCH3X(N,K)=0.00
      RCH4X(N,K)=0.50*RGOMO(N,K)
      ROXYO(N,K)=ROXYM(N,K)
      RH2GX(N,K)=0.0
      ELSEIF(K.EQ.5)THEN
      RCO2X(N,K)=0.00
      RCH3X(N,K)=0.00
      RCH4X(N,K)=RGOMO(N,K)
      ROXYO(N,K)=ROXYM(N,K)
      RH2GX(N,K)=0.0
      RH2GZ=0.667*RGOMO(N,K)
      ENDIF
      ENDIF
C
C     HETEROTROPHIC DENITRIFICATION
C
      IF(K.LE.4.AND.N.EQ.2.AND.ROXYM(N,K).GT.0.0
     2.AND.(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX)))THEN
C
C     FACTOR TO CONSTRAIN NO3 UPAKE AMONG COMPETING MICROBIAL
C     AND ROOT POPULATIONS
C
C     FNO3,FNB3=fraction of total biological demand for NO3
C        by denitrifiers
C     RVMX3,RVMB3=demand for NO3 reduction in non-band,band (g N t-1)
C     RNO3Y,RN3BY=total demand for NO3 reduction by all microbial, root
C        and mycorrhizal populations in non-band,band 
C        from RNO3X,RN3BX in ‘redist.f’ (g N t-1)
C
      IF(RNO3Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO3=AMAX1(FMN,RVMX3(N,K,L,NY,NX)/RNO3Y(L,NY,NX))
      ELSE
      FNO3=AMAX1(FMN,FOMA(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN3BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB3=AMAX1(FMN,RVMB3(N,K,L,NY,NX)/RN3BY(L,NY,NX))
      ELSE
      FNB3=AMAX1(FMN,FOMA(N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO3X=TFNO3X+FNO3
      TFNO3B=TFNO3B+FNB3
C
C     NO3 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO3
C     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
C     NOT ACCEPTED BY O2 IN BAND AND NON-BAND SOIL ZONES
C
C     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO (g O t-1)
C     VMXD3=demand for NO3-N reduction (g N t-1)
C     VMXDXS,VMXDXB=NO3 concentration-limited NO3 reduction in 
C        non-band, band with NO2 product feedback (g N t-1)
C     FNO3S,FNO3B=fractions of total NO3 in non-band, band
C     CNO3S,CNO3B=NO3 concentrations in non-band, band (g N m-3)
C     CNO2S,CNO2B=NO2 concentrations in non-band, band (g N m-3)
C     Z3KM,Z2KM=Km for NO3, NO2 uptake (g N m-3)
C     FVMXDX=nonlinear effect of product inhibition for NOx reduction
C     VMKI=product inhibition for NOx reduction (g N m-3 t-1)
C     VOLWY=biologically active water volume (m3)
C     FOSAH=fraction of TSRH in each substrate complex K
C     TSRH=total colonized SOC+microbial litter+adsorbed C (g C)
C     VMXD3S,VMXD3B=NO3 concentration-limited reduction with non-linear
C        effect in non-band,band (g N t-1)
C     OQCZ3=DOC supply limitation to NO3 reduction (g C t-1) 
C     FOQC=fraction of OQC biological demand by denitrifiers
C     OQCD3S,OQCD3B=DOC supply limitation to NO3 reduction in 
C        non-band, band (g N t-1)
C     FNO3S,FNO3B=fractions of NO3 in non-band, band
C     ZNO3SX,ZNO3BX=NO3 supply limitation to NO3 reduction (g N t-1)
C     FNO3,FNB3=fraction of total biological demand for NO3
C        by denitrifiers
C     ECN3=C:N ratio for e- transfers to NO3 (g C g N-1) 
C     RDNO3,RDNOB=NO3+DOC supply-limited NO3 reduction in non-band,band
C        (g N t-1)
C     RDNOT=total NO3 reduction in non-band+band (g N t-1)
C     RGOM3X,RGOMD3=NO3 and DOC-unlimited,-limited respiration 
C        from NO3 reduction (g C t-1)
C     RVMX3,RVMB3=demand for NO3 reduction in non-band,band (g N t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      ROXYD=AMAX1(0.0,ROXYM(N,K)-ROXYO(N,K))
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
      IF(VOLWY.GT.ZEROS2(NY,NX).AND.FOSAH(K,L,NY,NX).GT.ZERO)THEN
      FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWY*FOSAH(K,L,NY,NX)))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD3S=VMXDXS*FVMXDX 
      VMXD3B=VMXDXB*FVMXDX 
      OQCZ3=AMAX1(0.0,OQC(K,L,NY,NX)*FOQC-RGOCP*WFN(N,K))
      OQCD3=OQCZ3/ECN3*XNFH
      OQCD3S=OQCD3*FNO3S
      OQCD3B=OQCD3*FNO3B
      ZNO3SX=ZNO3S(L,NY,NX)*FNO3*XNFH 
      ZNO3BX=ZNO3B(L,NY,NX)*FNB3*XNFH 
      RDNO3X=AMAX1(0.0,AMIN1(ZNO3SX,VMXD3S))
      RDNOBX=AMAX1(0.0,AMIN1(ZNO3BX,VMXD3B))
      RDNO3(N,K)=AMAX1(0.0,AMIN1(VMXD3S,OQCD3S,ZNO3SX))
      RDNOB(N,K)=AMAX1(0.0,AMIN1(VMXD3B,OQCD3B,ZNO3BX))
      RDNOX=RDNO3X+RDNOBX
      RDNOT=RDNO3(N,K)+RDNOB(N,K)
      RGOM3X=ECN3*RDNOX
      RGOMD3=ECN3*RDNOT
      RVMX3(N,K,L,NY,NX)=VMXD3S 
      RVMB3(N,K,L,NY,NX)=VMXD3B
C
C     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
C     POPULATIONS
C
C     FNO2,FNB2=fraction of total biological demand for NO2
C        by denitrifiers in non-band,band 
C     RVMX2,RVMB2=demand for NO2 reduction in non-band,band (g N t-1)
C     RNO2Y,RN2BY=total demand for NO2 reduction by all microbial, root
C        and mycorrhizal populations in non-band,band 
C        from RNO2X,RN2BX in ‘redist.f’ (g N t-1)
C
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMX2(N,K,L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=AMAX1(FMN,FOMA(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMB2(N,K,L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=AMAX1(FMN,FOMA(N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
C
C     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2
C     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
C     NOT ACCEPTED BY O2 AND NO3 IN BAND AND NON-BAND SOIL ZONES
C
C     VMXD2=demand for NO2-N reduction (g N t-1)
C     VMXD3=demand for NO3-N reduction (g N t-1)
C     RDNOT=total NO3 reduction in non-band+band (g N t-1)
C     VMXDXS,VMXDXB=NO2 concentration-limited NO2 reduction in 
C        non-band, band with N2O product feedback (g N t-1)
C     FNO2S,FNO2B=fractions of total NO2 in non-band, band
C     CNO2S,CNO2B=NO2 concentrations in non-band, band (g N m-3)
C     CZ2OS=N2O concentration (g N m-3)
C     Z2KM,Z1KM=Km for NO2, N2O uptake (g N m-3)
C     FVMXDX=nonlinear effect of product inhibition for NOx reduction
C     VMKI=product inhibition for NOx reduction (g N m-3 t-1)
C     VOLWY=biologically active water volume (m3)
C     FOSAH=fraction of TSRH in each substrate complex K
C     TSRH=total colonized SOC+microbial litter+adsorbed C (g C)
C     VMXD2S,VMXD2B=NO2 concentration-limited reduction with non-linear
C        effect in non-band,band (g N t-1)
C     OQCZ2=DOC supply limitation to NO2 reduction (g C t-1) 
C     FOQC=fraction of OQC biological demand by denitrifiers
C     OQCD2S,OQCD2B=DOC supply limitation to NO2 reduction in 
C        non-band, band (g N t-1)
C     FNO3S,FNO3B=fractions of NO3 in non-band, band
C     ZNO2SX,ZNO2BX=NO2 supply limitation to NO2 reduction (g N t-1)
C     FNO2,FNB2=fraction of total biological demand for NO2
C        by denitrifiers
C     ECN2=C:N ratio for e- transfers to NO2 (g C g N-1) 
C     RDNO2,RDN2B=NO2+DOC supply-limited NO2 reduction in non-band,band
C        (g N t-1)
C     RDN2T=total NO2 reduction in non-band+band (g N t-1)
C     RGOM2X,RGOMD2=NO2 and DOC-unlimited,-limited respiration 
C        from NO2 reduction (g C t-1)
C     RVMX2,RVMB2=demand for NO2 reduction in non-band,band (g N t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
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
      IF(VOLWY.GT.ZEROS2(NY,NX).AND.FOSAH(K,L,NY,NX).GT.ZERO)THEN
      FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWY*FOSAH(K,L,NY,NX)))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD2S=VMXDXS*FVMXDX 
      VMXD2B=VMXDXB*FVMXDX 
      OQCZ2=AMAX1(0.0,OQCZ3-RGOMD3)
      OQCD2=OQCZ2/ECN2*XNFH 
      OQCD2S=OQCD2*FNO3S
      OQCD2B=OQCD2*FNO3B
      ZNO2SX=(ZNO2S(L,NY,NX)*XNFH+RDNO3(N,K))*FNO2
      ZNO2BX=(ZNO2B(L,NY,NX)*XNFH+RDNOB(N,K))*FNB2 
      RDNO2X=AMAX1(0.0,AMIN1(ZNO2SX,VMXD2S))
      RDNOBX=AMAX1(0.0,AMIN1(ZNO2BX,VMXD2B))
      RDNO2(N,K)=AMAX1(0.0,AMIN1(VMXD2S,OQCD2S,ZNO2SX))
      RDN2B(N,K)=AMAX1(0.0,AMIN1(VMXD2B,OQCD2B,ZNO2BX))
      RDN2X=RDNO2X+RDNOBX
      RDN2T=RDNO2(N,K)+RDN2B(N,K)
      RGOM2X=ECN2*RDN2X
      RGOMD2=ECN2*RDN2T
      RVMX2(N,K,L,NY,NX)=VMXD2S 
      RVMB2(N,K,L,NY,NX)=VMXD2B
C
C     FACTOR TO CONSTRAIN N2O UPAKE AMONG COMPETING MICROBIAL
C     AND ROOT POPULATIONS
C 
C     FN2O=fraction of total biological demand for N2O
C        by denitrifiers
C     RVMX1=demand for N2O2 reduction (g N t-1)
C     RN2OY=total demand for N2O reduction by all microbial, root
C        and mycorrhizal populations in non-band,band 
C        from RN2OX in ‘redist.f’ (g N t-1)
C
      IF(RN2OY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FN2O=AMAX1(FMN,RVMX1(N,K,L,NY,NX)/RN2OY(L,NY,NX))
      ELSE
      FN2O=AMAX1(FMN,FOMA(N,K))
      ENDIF
      TFN2OX=TFN2OX+FN2O
C
C     N2O REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE DENITRIFIER BIOMASS, TEMPERATURE, AQUEOUS N2O
C     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
C     NOT ACCEPTED BY O2, NO3 AND NO2 IN BAND AND NON-BAND SOIL ZONES
C
C     VMXD1=demand for N2O-N reduction
C     VMXD2=demand for NO2-N reduction (g N t-1)
C     RDN2T=total NO2 reduction in non-band+band (g N t-1)
C     VMXDXS=N2O concentration-limited N2O reduction (g N t-1)
C     CZ2OS=N2O concentration (g N m-3)
C     Z1KM=Km for N2O uptake (g N m-3)
C     FVMXDX=nonlinear effect of product inhibition for NOx reduction
C     VMKI=product inhibition for NOx reduction (g N m-3 t-1)
C     VOLWY=biologically active water volume (m3)
C     FOSAH=fraction of TSRH in each substrate complex K
C     TSRH=total colonized SOC+microbial litter+adsorbed C (g C)
C     VMXD1S=NO2 concentration-limited N2O reduction (g N t-1)
C     OQCZ1=DOC supply limitation to N2O reduction (g C t-1) 
C     OQCZ2=DOC supply limitation to NO2 reduction (g C t-1) 
C     RGOMD2=NO2 and DOC-limited respiration 
C        from NO2 reduction (g C t-1)
C     OQCD1=DOC limitation to N2O reduction (g N t-1)
C     ECN1=C:N ratio for e- transfers to N2O (g C g N-1) 
C     ZN2OSX=N2O supply limitation to N2O reduction (g N t-1)
C     RDN2O=N2O+DOC supply-limited N2O reduction (g N g t-1)
C     RGOM1X,RGOMD1=NO2 and DOC-unlimited,-limited respiration 
C        from N2O reduction (g C t-1)
C     RGOMY,RGOMD=total substrate-unlimited,-limited respiration 
C        from NOx reduction (g C t-1) 
C     RVMX1=demand for N2O reduction (g N t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      VMXD1=(VMXD2-RDN2T)*2.0
      VMXDXS=VMXD1*CZ2OS(L,NY,NX)/(CZ2OS(L,NY,NX)+Z1KM)
      IF(VOLWY.GT.ZEROS2(NY,NX).AND.FOSAH(K,L,NY,NX).GT.ZERO)THEN
      FVMXDX=1.0/(1.0+VMXDXS/(VMKI*VOLWY*FOSAH(K,L,NY,NX)))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD1S=VMXDXS*FVMXDX 
      OQCZ1=AMAX1(0.0,OQCZ2-RGOMD2)
      OQCD1=OQCZ1/ECN1*XNFH
      Z2OSX=(Z2OS(L,NY,NX)*XNFH+RDN2T)*FN2O 
      RDN2OX=AMAX1(0.0,AMIN1(Z2OSX,VMXD1S))
      RDN2O(N,K)=AMAX1(0.0,AMIN1(VMXD1S,OQCD1,Z2OSX))
      RGOM1X=ECN1*RDN2OX
      RGOMD1=ECN1*RDN2O(N,K)
      RGOMY(N,K)=RGOM3X+RGOM2X+RGOM1X
      RGOMD(N,K)=RGOMD3+RGOMD2+RGOMD1
      RVMX1(N,K,L,NY,NX)=VMXD1S 
C     IF((I/1)*1.EQ.I.AND.L.LE.5)THEN
C     WRITE(*,2222)'DENIT',I,J,L,K,N,RDNO3(N,K),RDNOB(N,K),RDNO2(N,K)
C    2,RDN2B(N,K),RDN2O(N,K) 
C    3,COXYS(L,NY,NX),COXYG(L,NY,NX),ROXYM(N,K)
C    3,ROXYO(N,K),OMA(N,K),VMXD,CNO3S(L,NY,NX),CNO3B(L,NY,NX) 
C    4,CNO2S(L,NY,NX),CNO2B(L,NY,NX),CZ2OS(L,NY,NX),VLNO3(L,NY,NX)
C    5,VLNOB(L,NY,NX),THETW(L,NY,NX),THETI(L,NY,NX),FOMA(N,K)
C    5,ZNO3S(L,NY,NX),ZNO3B(L,NY,NX),ZNO2S(L,NY,NX),ZNO2B(L,NY,NX)
C    6,Z2OS(L,NY,NX),RGOMY(N,K),RGOMD(N,K),TOMA,FOXYX,FNO23S,FNO23B 
C    7,OQC(K,L,NY,NX),FOQC,RGOCP,WFN(N,K),VOLWY,FOSAH(K,L,NY,NX),ZERO 
C    9,RGOM3X,RGOM2X,RGOM1X,FNO3,FNO2,FN2O,ZNO3SX,ZNO2SX,Z2OSX 
C    3,OQCD3S,OQCD2S,OQCD1,VMXD3S,VMXD2S,VMXD1S,VMXD3,VMXD2,VMXD1
C    4,ROXYD,VMXDX,TFNX,WFNG,TFNG(N,K),PSISM(L,NY,NX) 
C    2,(1.0+(CNO2S(L,NY,NX)*Z3KM)/(CNO3S(L,NY,NX)*Z2KM))
C    2,(1.0+(CZ2OS(L,NY,NX)*Z2KM)/(CNO2S(L,NY,NX)*Z1KM))
2222  FORMAT(A8,5I4,70E12.4)
C     ENDIF
C
C     AUTOTROPHIC DENITRIFICATION
C
      ELSEIF(K.EQ.5.AND.N.EQ.1.AND.ROXYM(N,K).GT.0.0
     2.AND.(L.NE.0.OR.VOLX(L,NY,NX).GT.ZEROS(NY,NX)))THEN
C
C     FACTOR TO CONSTRAIN NO2 UPAKE AMONG COMPETING MICROBIAL
C     POPULATIONS
C
C     FNO2,FNB2=fraction of total biological demand for NO2
C        by nitrifiers in non-band,band 
C     RVMX2,RVMB2=demand for NO2 reduction in non-band,band (g N t-1)
C     RNO2Y,RN2BY=total demand for NO2 reduction by all microbial, root
C        and mycorrhizal populations in non-band,band 
C        from RNO2X,RN2BX in ‘redist.f’ (g N t-1)
C
      IF(RNO2Y(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNO2=AMAX1(FMN,RVMX2(N,K,L,NY,NX)/RNO2Y(L,NY,NX))
      ELSE
      FNO2=AMAX1(FMN,FOMN(N,K)*VLNO3(L,NY,NX))
      ENDIF
      IF(RN2BY(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FNB2=AMAX1(FMN,RVMB2(N,K,L,NY,NX)/RN2BY(L,NY,NX))
      ELSE
      FNB2=AMAX1(FMN,FOMN(N,K)*VLNOB(L,NY,NX))
      ENDIF
      TFNO2X=TFNO2X+FNO2
      TFNO2B=TFNO2B+FNB2
C
C     NO2 REDUCTION FROM SPECIFIC REDUCTION RATE, ENERGY YIELD,
C     ACTIVE NITRIFIER BIOMASS, TEMPERATURE, AQUEOUS NO2 AND CO2
C     CONCENTRATIONS AND STOICHIOMETRY OF REDOX ELECTRON TRANSFER
C     NOT ACCEPTED BY O2
C
C     ROXYD=O2 demand ROXYM not met by O2 uptake ROXYO (g O t-1)
C     VMXD4=demand for NO2-N reduction (g N t-1)
C     XCO2=aqueous CO2 limitation to CO2 reduction
C     VMXDXS,VMXDXB=NO2 concentration-limited NO2 reduction in 
C        non-band, band with N2O product feedback (g N t-1)
C     FNO2S,FNO2B=fractions of total NO2 in non-band, band
C     CNO2S,CNO2B=NO2 concentrations in non-band, band (g N m-3)
C     Z2KM=Km for NO2 uptake (g N m-3)
C     FVMXDX=nonlinear effect of product inhibition for NOx reduction
C     VMKI=product inhibition for NOx reduction (g N m-3 t-1)  
C     VOLWY=biologically active water volume (m3)
C     VMXD4S,VMXD4B=NO2 concentration-limited reduction with non-linear
C        effect in non-band,band (g N t-1)
C     ZNO2SX,ZNO2BX=NO2 supply limitation to NO2 reduction (g N t-1)
C     RVOXA,RVOXB=O2-limited oxidation in non-band,band (g O t-1) of:
C        N=1:NH4
C        N=2:NO2
C     RDNO2,RDN2B=NO2 supply-limited NO2 reduction in non-band,band 
C        (g N t-1)
C     RDNOT=total NO2 reduction in non-band+band (g N t-1)
C     RGOMD=total respiration from NO2 reduction (g C t-1)
C     ECNO=efficiency CO2 conversion to biomass by nitrite oxidizers
C        (g C g N-1) 
C     RVMX2,RVMB2=demand for NO2 reduction in non-band,band (g N t-1)
C     ENOX=growth respiration efficiency for NO2 reduction (g C g C-1)
C
      ROXYD=AMAX1(0.0,ROXYM(N,K)-ROXYO(N,K))
      VMXD4=0.875*ROXYD*XCO2
      VMXDXS=FNO2S*VMXD4*CNO2S(L,NY,NX)/(CNO2S(L,NY,NX)+Z2KM) 
      VMXDXB=FNO2B*VMXD4*CNO2B(L,NY,NX)/(CNO2B(L,NY,NX)+Z2KM)
      VMXDXT=VMXDXS+VMXDXB
      IF(VOLWY.GT.ZEROS2(NY,NX))THEN 
      FVMXDX=1.0/(1.0+VMXDXT/(VMKI*VOLWY))
      ELSE
      FVMXDX=0.0
      ENDIF
      VMXD4S=VMXDXS*FVMXDX 
      VMXD4B=VMXDXB*FVMXDX 
      ZNO2SX=ZNO2S(L,NY,NX)+RVOXA(1)
      ZNO2BX=ZNO2B(L,NY,NX)+RVOXB(1)
      RDNO2(N,K)=AMAX1(0.0,AMIN1(VMXD4S,ZNO2SX))
      RDN2B(N,K)=AMAX1(0.0,AMIN1(VMXD4B,ZNO2BX))
      RDNOT=RDNO2(N,K)+RDN2B(N,K)
      RGOMY(N,K)=0.0
      RGOMD(N,K)=RDNOT*ECNO*ENOX
      RDNO3(N,K)=0.0
      RDNOB(N,K)=0.0
      RDN2O(N,K)=0.0
      RVMX2(N,K,L,NY,NX)=VMXD4S 
      RVMB2(N,K,L,NY,NX)=VMXD4B
      RVOXA(N)=RVOXA(N)+0.333*RDNO2(N,K)
      RVOXB(N)=RVOXB(N)+0.333*RDN2B(N,K)
C     IF((I/1)*1.EQ.I.AND.J.EQ.19.AND.L.LE.5)THEN
C     WRITE(*,7777)'AUTO',I,J,L,K,N,RDNO2(N,K)
C    2,RDN2B(N,K)
C    2,CNO2S(L,NY,NX),CNO2B(L,NY,NX),CNO3S(L,NY,NX),CNO3B(L,NY,NX)
C    3,Z2OS(L,NY,NX),VLNOB(L,NY,NX),ZNO2S(L,NY,NX),ZNO2B(L,NY,NX)
C    3,XCO2,FNO2,FNB2,TFNG(N,K),OMA(N,K),ROXYP(N,K)
C    2,ROXYM(N,K),ROXYO(N,K),WFN(N,K),FOXYX
C    3,THETW(L,NY,NX),COXYS(L,NY,NX),COXYG(L,NY,NX)
C    4,ROXYD,VMXD4,VMXDXS,VMXDXB,VMXD4S,VMXD4B,FNO2S,FNO2B
C    5,ZNFN4S,ZNFN4B 
7777  FORMAT(A8,5I4,50E12.4)
C     ENDIF
      ELSE
      RDNO3(N,K)=0.0
      RDNOB(N,K)=0.0
      RDNO2(N,K)=0.0
      RDN2B(N,K)=0.0
      RDN2O(N,K)=0.0
      RGOMY(N,K)=0.0
      RGOMD(N,K)=0.0
      ENDIF
C
C     N AND P MINERALIZATION-IMMOBILIZATION
C
C     MINERALIZATION-IMMOBILIZATION OF NH4 IN SOIL FROM MICROBIAL
C     C:N AND NH4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
C
C     RINHP=NH4 mineralization (-ve) or immobilization (+ve) demand 
C        (g N t-1) 
C     OMC,OMN=microbial nonstructural C,N (g C,N)
C     CNOMC=maximum microbial N:C ratio from ’starts.f’ (g N g C-1)
C     CNH4S,CNH4B=NH4 concentrations in non-band, band (g N m-3)
C     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate (g N m-2 h-1,
C        minimum NH4 concentration (g N m-3) and Km for NH4 uptake 
C        (g N m-3) 
C     RINHX=microbial limitation to NH4 demand (g N t-1)
C     BIOA=microbial surface area (m2 g C-1)
C     OMA=active biomass (g C)
C     TFNG=temperature+water limitation to biological activity
C     FNH4S,FNHBS=fractions of NH4 in non-band, band
C     RINHO,RINHB=NH4 supply-unlimited NH4 mineralization-
C        immobilization in non-band, band (g N t-1) 
C     ZNH4M,ZNHBM=NH4 not available for uptake in non-band, band (g N)
C     VOLW=soil water content (m3)
C     RINH4,RINB4=NH4 supply-limited NH4 mineralization-immobilization
C        in non-band,band (g N t-1)
C     FNH4X,FNB4X=fractions of biological NH4 demand in non-band, band
C        by microbial population
C     TRINH4=total NH4 net mineralization(-ve) or immobilization(+ve)
C        (g N t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C     
      RINHP=(OMC(3,N,K,L,NY,NX)*CNOMC(3,N,K)-OMN(3,N,K,L,NY,NX))
      IF(RINHP.GT.0.0)THEN
      CNH4X=AMAX1(0.0,CNH4S(L,NY,NX)-Z4MN)
      CNH4Y=AMAX1(0.0,CNH4B(L,NY,NX)-Z4MN)
      RINHX=AMIN1(RINHP,BIOA*OMA(N,K)*TFNG(N,K)*Z4MX*XNFH) 
      RINHO(N,K,L,NY,NX)=FNH4S*RINHX*CNH4X/(CNH4X+Z4KU)
      RINHB(N,K,L,NY,NX)=FNHBS*RINHX*CNH4Y/(CNH4Y+Z4KU)
      ZNH4M=Z4MN*VOLW(L,NY,NX)*FNH4S
      ZNHBM=Z4MN*VOLW(L,NY,NX)*FNHBS
      RINH4(N,K)=AMIN1(FNH4X*AMAX1(0.0,(ZNH4S(L,NY,NX)-ZNH4M)*XNFH)
     2,RINHO(N,K,L,NY,NX)) 
      RINB4(N,K)=AMIN1(FNB4X*AMAX1(0.0,(ZNH4B(L,NY,NX)-ZNHBM)*XNFH)
     2,RINHB(N,K,L,NY,NX))
      ELSE
      RINHO(N,K,L,NY,NX)=0.0
      RINHB(N,K,L,NY,NX)=0.0
      RINH4(N,K)=RINHP*FNH4S
      RINB4(N,K)=RINHP*FNHBS
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+(RINH4(N,K)+RINB4(N,K))
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,7776)'RINH4',I,J,NFZ,NX,NY,L,K,N,RINH4(N,K),RINHP
C    1,BIOA*OMA(N,K)*Z4MX*TFNG(N,K),BIOA,OMA(N,K),Z4MX,TFNG(N,K) 
C    2,OMC(3,N,K,L,NY,NX),CNOMC(3,N,K),OMN(3,N,K,L,NY,NX)
C    3,RINHO(N,K,L,NY,NX),CNH4S(L,NY,NX),FNH4X,ZNH4S(L,NY,NX)
C    4,ZNH4B(L,NY,NX),ZNH4T(L),OQN(K,L,NY,NX),TRINH4(NY,NX)
7776  FORMAT(A8,8I4,30E12.4)
C     ENDIF
C
C     MINERALIZATION-IMMOBILIZATION OF NO3 IN SOIL FROM MICROBIAL
C     C:N AND NO3 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
C
C     RINOP=NO3 immobilization (+ve) demand (g N t-1)
C     RINH4,RINB4=NH4 supply-limited NH4 mineralization-immobilization
C        in non-band,band (g N t-1)
C     CNO3S,CNO3B=NO3 concentrations in non-band, band (g N m-3)
C     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate (g N m-2 h-1,
C        minimum NO3 concentration (g N m-3) and Km for NO3 uptake 
C        (g N m-3) 
C     RINOX=microbial limitation to NO3 demand (g N t-1)
C     BIOA=microbial surface area (m2 g C-1)
C     OMA=active biomass (g C)
C     TFNG=temperature+water limitation to biological activity
C     FNO3S,FNO3B=fractions of NO3 in non-band, band
C     RINOO,RINOB=NO3 supply-unlimited NO3 immobilization (g N t-1) 
C     VOLW=soil water content (m3)
C     ZNO3M,ZNOBM=NO3 not available for uptake in non-band, band (g N)
C     RINO3,RINB3=NO3 supply-limited NO3 immobilization in non-band,
C        band 
C     FNO3X,FNB3X=fractions of biological NO3 demand in non-band, band
C        by microbial population (g N t-1)
C     TRINH4=total net NH4+NO3 mineralization(-ve) or immobilization
C        (+ve) (g N t-1) 
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      RINOP=AMAX1(0.0,RINHP-RINH4(N,K)-RINB4(N,K))
      IF(RINOP.GT.0.0)THEN
      CNO3X=AMAX1(0.0,CNO3S(L,NY,NX)-ZOMN)
      CNO3Y=AMAX1(0.0,CNO3B(L,NY,NX)-ZOMN)
      RINOX=AMIN1(RINOP,BIOA*OMA(N,K)*TFNG(N,K)*ZOMX*XNFH)
      RINOO(N,K,L,NY,NX)=FNO3S*RINOX*CNO3X/(CNO3X+ZOKU)
      RINOB(N,K,L,NY,NX)=FNO3B*RINOX*CNO3Y/(CNO3Y+ZOKU)
      ZNO3M=ZOMN*VOLW(L,NY,NX)*FNO3S
      ZNOBM=ZOMN*VOLW(L,NY,NX)*FNO3B
      RINO3(N,K)=AMIN1(FNO3X*AMAX1(0.0,(ZNO3S(L,NY,NX)-ZNO3M)*XNFH)
     2,RINOO(N,K,L,NY,NX)) 
      RINB3(N,K)=AMIN1(FNB3X*AMAX1(0.0,(ZNO3B(L,NY,NX)-ZNOBM)*XNFH)
     2,RINOB(N,K,L,NY,NX)) 
      ELSE
      RINOO(N,K,L,NY,NX)=0.0
      RINOB(N,K,L,NY,NX)=0.0
      RINO3(N,K)=RINOP*FNO3S
      RINB3(N,K)=RINOP*FNO3B
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+(RINO3(N,K)+RINB3(N,K))
C     IF(RINO3(N,K).LT.0.0.OR.RINB3(N,K).LT.0.0)THEN
C     WRITE(*,4321)'RINO3',I,J,NX,NY,L,K,N
C    2,RINOO(N,K,L,NY,NX),RINO3(N,K)
C    2,RINOP,BIOA,OMA(N,K),TFNG(N,K),ZOMX,WFN(N,K),FNO3X,FNO3B 
C    2,VLNO3(L,NY,NX),VLNOB(L,NY,NX),CNO3S(L,NY,NX),CNO3B(L,NY,NX)
C    3,RINOB(N,K,L,NY,NX),RINB3,ZNO3S(L,NY,NX),ZNO3B(L,NY,NX)
C    3,OMC(3,N,K,L,NY,NX),CPOMC(3,N,K),OMP(3,N,K,L,NY,NX),WFNG 
4321  FORMAT(A8,7I4,60E12.4)
C     ENDIF
C
C     MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SOIL FROM MICROBIAL
C     C:P AND H2PO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
C
C     RIPOP=H2PO4 mineralization (-ve) or immobilization (+ve) demand
C        (g P t-1)
C     OMC,OMP=microbial nonstructural C,P (g C,P)
C     CPOMC=maximum microbial P:C ratio from ‘starts.f’ (g P g C-1)
C     CH2P4,CH2P4B=H2PO4 concentrations in non-band, band (g P m-3)
C     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate (g P m-2 h-1,
C        minimum H2PO4 concentration (g P m-3) and Km for H2PO4 uptake 
C        (g P m-3) 
C     RIPOX=microbial limitation to H2PO4 demand (g P t-1)
C     BIOA=microbial surface area (m2 g C-1)
C     OMA=active biomass (g C)
C     TFNG=temperature+water limitation to biological activity
C     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
C     RIPOO,RIPBO=H2PO4 supply-unlimited H2PO4 mineralization-
C        immobilization (g P t-1) 
C     H2POM,H2PBM=H2PO4 not available for uptake in non-band, band 
C        (g P)
C     VOLW=soil water content (m3)
C     FPO4X,FPOBX=fractions of biological H2PO4 demand in non-band,
C        band by microbial population
C     RIPO4,RIPOB=H2PO4 supply-limited H2PO4 mineralization-
C        immobilization in non-band,band (g P t-1)
C     TRIPO4=total H2PO4 net mineralization (-ve) or immobilization
C        (+ve) (g P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      RIPOP=(OMC(3,N,K,L,NY,NX)*CPOMC(3,N,K)-OMP(3,N,K,L,NY,NX))
      IF(RIPOP.GT.0.0)THEN
      CH2PX=AMAX1(0.0,CH2P4(L,NY,NX)-HPMN)
      CH2PY=AMAX1(0.0,CH2P4B(L,NY,NX)-HPMN)
      RIPOX=AMIN1(RIPOP,BIOA*OMA(N,K)*TFNG(N,K)*HPMX*XNFH)
      RIPOO(N,K,L,NY,NX)=FH2PS*RIPOX*CH2PX/(CH2PX+HPKU) 
      RIPBO(N,K,L,NY,NX)=FH2PB*RIPOX*CH2PY/(CH2PY+HPKU) 
      H2POM=HPMN*VOLW(L,NY,NX)*FH2PS 
      H2PBM=HPMN*VOLW(L,NY,NX)*FH2PB
      RIPO4(N,K)=AMIN1(FPO4X*AMAX1(0.0,(H2PO4(L,NY,NX)-H2POM)*XNFH)
     2,RIPOO(N,K,L,NY,NX)) 
      RIPOB(N,K)=AMIN1(FPOBX*AMAX1(0.0,(H2POB(L,NY,NX)-H2PBM)*XNFH)
     2,RIPBO(N,K,L,NY,NX)) 
      ELSE
      RIPOO(N,K,L,NY,NX)=0.0
      RIPBO(N,K,L,NY,NX)=0.0
      RIPO4(N,K)=RIPOP*FH2PS 
      RIPOB(N,K)=RIPOP*FH2PB 
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+(RIPO4(N,K)+RIPOB(N,K))
C     IF(NY.EQ.5.AND.L.EQ.10.AND.K.EQ.3.AND.N.EQ.2)THEN
C     WRITE(*,4322)'RIPO4',I,J,NX,NY,L,K,N,RIPO4(N,K),FPO4X,H2P4T(L)
C    2,RIPOO(N,K,L,NY,NX),RIPOP,BIOA,OMA(N,K),TFNG(N,K),HPMX,WFN(N,K) 
C    2,VLPO4(L,NY,NX),VLPOB(L,NY,NX),CH2P4(L,NY,NX),CH2P4B(L,NY,NX)
C    3,OMC(3,N,K,L,NY,NX),CPOMC(3,N,K),OMP(3,N,K,L,NY,NX),WFNG 
4322  FORMAT(A8,7I4,30E12.4)
C     ENDIF
C
C     MINERALIZATION-IMMOBILIZATION OF HPO4 IN SOIL FROM MICROBIAL
C     C:P AND HPO4 CONCENTRATION IN BAND AND NON-BAND SOIL ZONES
C
C     RIP1P=HPO4 mineralization (-ve) or immobilization (+ve) demand
C        (g P t-1)
C     RIPO4,RIPOB=H2PO4 supply-limited H2PO4 mineralization-
C        immobilization in non-band,band (g P t-1)
C     CH1P4,CH1P4B=HPO4 concentrations in non-band, band (g P m-3)
C     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate (g P m-2 h-1,
C        minimum HPO4 concentration (g P m-3) and Km for HPO4 uptake 
C        (g P m-3), assumed same as H2PO4
C     RIP1X=microbial limitation to HPO4 demand (g P m-3)
C     BIOA=microbial surface area (m2 g C-1)
C     OMA=active biomass (g C)
C     TFNG=temperature+water limitation to biological activity
C     FH1PS,FH1PB=fractions of HPO4 in non-band, band
C     RIPO1,RIPB1=HPO4 supply-unlimited HPO4 mineralization-
C        immobilization in non-band,band (g P t-1) 
C     H1POM,H1PBM=HPO4 not available for uptake in non-band, band
C        (g P)
C     VOLW=soil water content (m3)
C     FP14X,FP1BX=fractions of biological HPO4 demand in non-band, band
C     RIP14,RIP1B=HPO4 supply-limited HPO4 mineralization-
C        immobilizationin non-band,band (g P t-1) 
C     TRIPO4=total H2PO4+HPO4 net mineralization (-ve) or
C        immobilization (+ve) (g P t-1) 
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      RIP1P=0.1*AMAX1(0.0,RIPOP-RIPO4(N,K)-RIPOB(N,K))
      IF(RIP1P.GT.0.0)THEN
      CH1PX=AMAX1(0.0,CH1P4(L,NY,NX)-HPMN)
      CH1PY=AMAX1(0.0,CH1P4B(L,NY,NX)-HPMN)
      RIP1X=AMIN1(RIP1P,BIOA*OMA(N,K)*TFNG(N,K)*HPMX*XNFH)
      RIPO1(N,K,L,NY,NX)=FH1PS*RIP1X*CH1PX/(CH1PX+HPKU) 
      RIPB1(N,K,L,NY,NX)=FH1PB*RIP1X*CH1PY/(CH1PY+HPKU) 
      H1POM=HPMN*VOLW(L,NY,NX)*FH1PS 
      H1PBM=HPMN*VOLW(L,NY,NX)*FH1PB
      RIP14(N,K)=AMIN1(FP14X*AMAX1(0.0,(H1PO4(L,NY,NX)-H1POM)*XNFH)
     2,RIPO1(N,K,L,NY,NX)) 
      RIP1B(N,K)=AMIN1(FP1BX*AMAX1(0.0,(H1POB(L,NY,NX)-H1PBM)*XNFH)
     2,RIPB1(N,K,L,NY,NX)) 
      ELSE
      RIPO1(N,K,L,NY,NX)=0.0
      RIPB1(N,K,L,NY,NX)=0.0
      RIP14(N,K)=RIP1P*FH1PS 
      RIP1B(N,K)=RIP1P*FH1PB 
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+(RIP14(N,K)+RIP1B(N,K))
C     IF(NY.EQ.5.AND.L.EQ.10.AND.K.EQ.3.AND.N.EQ.2)THEN
C     WRITE(*,4323)'RIP14',I,J,NX,NY,L,K,N,RIP14(N,K),FP14X,H1P4T(L)
C    2,RIPO1(N,K,L,NY,NX),RIP1P,BIOA,OMA(N,K),TFNG(N,K),HPMX,WFN(N,K) 
C    2,VLPO4(L,NY,NX),VLPOB(L,NY,NX),CH1P4(L,NY,NX),CH1P4B(L,NY,NX)
C    3,OMC(3,N,K,L,NY,NX),CPOMC(3,N,K),OMP(3,N,K,L,NY,NX),WFNG 
4323  FORMAT(A8,7I4,30E12.4)
C     ENDIF
C
C     ADDITIONAL MINERALIZATION-IMMOBILIZATION OF NH4 IN SURFACE
C     LITTER FROM MICROBIAL C:N AND NH4 CONCENTRATION IN BAND AND 
C     NON-BAND SOIL ZONES OF SOIL SURFACE
C
C     RINHPR=NH4 mineralization (-ve) or immobilization (+ve) demand 
C        in litter (g N t-1) 
C     RINH4,RINO3=NH4,NO3 supply-limited NH4,NO3 mineralization-
C        immobilization (g N t-1)
C     NU=surface layer number
C     CNH4S,CNH4B=NH4 concentrations in non-band, band (g N m-3)
C     Z4MX,Z4MN,Z4KU=parameters for max NH4 uptake rate (g N m-2 h-1,
C        minimum NH4 concentration (g N m-3) and Km for NH4 uptake 
C        (g N m-3) 
C     BIOA=microbial surface area (m2 g C-1)
C     OMA=active biomass (g C)
C     TFNG=temperature+water limitation to biological activity
C     FNH4S,FNHBS=fractions of NH4 in non-band, band
C     RINHOR=NH4 supply-unlimited NH4 mineralization-immobilization 
C        (g N t-1) 
C     ZNH4M=NH4 not available for uptake (g N)
C     VOLW=soil water content (m3)
C     RINH4R=NH4 supply-limited NH4 mineralization-immobilization 
C        (g N t-1)
C     FNH4XR=fraction of biological NH4 demand by microbial population
C     TRINH4=total NH4 net mineralization (-ve) or immobilization
C        (+ve) (g t-1)
C
      IF(L.EQ.0)THEN
      RINHPR=RINHP-RINH4(N,K)-RINO3(N,K)
      IF(RINHPR.GT.0.0)THEN
      CNH4X=AMAX1(0.0,CNH4S(NU(NY,NX),NY,NX)-Z4MN)
      CNH4Y=AMAX1(0.0,CNH4B(NU(NY,NX),NY,NX)-Z4MN)
      RINHOR(N,K,NY,NX)=AMIN1(RINHPR
     2,BIOA*OMA(N,K)*TFNG(N,K)*Z4MX*XNFH) 
     3*(FNH4S*CNH4X/(CNH4X+Z4KU)+FNHBS*CNH4Y
     4/(CNH4Y+Z4KU)) 
      ZNH4M=Z4MN*VOLW(NU(NY,NX),NY,NX) 
      RINH4R(N,K)=AMIN1(FNH4XR(N,K)*AMAX1(0.0
     2,(ZNH4T(NU(NY,NX))-ZNH4M)),RINHOR(N,K,NY,NX))
      ELSE
      RINHOR(N,K,NY,NX)=0.0
      RINH4R(N,K)=RINHPR
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+RINH4R(N,K)
C     IF(K.EQ.2.AND.N.EQ.1)THEN
C     WRITE(*,7778)'RINH4R',I,J,NFZ,NX,NY,L,K,N,RINH4R(N,K),RINHPR
C    2,BIOA*OMA(N,K)*Z4MX,RINHP,RINH4(N,K),RINO3(N,K) 
C    3,RINHOR(N,K,NY,NX),CNH4S(NU(NY,NX),NY,NX),FNH4XR(N,K)
C    4,ZNH4T(NU(NY,NX))
7778  FORMAT(A8,8I4,20E12.4)
C     ENDIF
C
C     ADDITIONAL MINERALIZATION-IMMOBILIZATION OF NO3 IN SURFACE
C     LITTER FROM MICROBIAL C:N AND NO3 CONCENTRATION IN BAND AND 
C     NON-BAND SOIL ZONES OF SOIL SURFACE
C
C     RINOPR=NO3 immobilization (+ve) demand in litter (g N t-1) 
C     NU=surface layer number
C     CNO3S,CNO3B=NO3 concentrations in non-band, band (g N m-3)
C     ZOMX,ZOMN,ZOKU=parameters for max NO3 uptake rate (g N m-2 h-1,
C        minimum NO3 concentration (g N m-3) and Km for NO3 uptake 
C        (g N m-3) 
C     RINOOR=microbial limitation to NO3 demand
C     BIOA=microbial surface area (m2 g C-1)
C     OMA=active biomass (g C)
C     TFNG=temperature+water limitation to biological activity
C     FNO3S,FNO3B=fractions of NO3 in non-band, band
C     RINO3R=NO3 supply-unlimited NO3 immobilization 
C     VOLW=soil water content (m3)
C     ZNO3M=NO3 not available for uptake (g N) 
C     FNO3XR=fraction of biological NO3 demand 
C        by microbial population (g N t-1)
C     RINO3R=NO3 supply-limited NO3 immobilization (g N t-1)
C     TRINH4=total NH4+NO3 net mineralization (-ve) or immobilization
C        (+ve) (g N t-1)
C
      RINOPR=AMAX1(0.0,RINHPR-RINH4R(N,K))
      IF(RINOPR.GT.0.0)THEN
      CNO3X=AMAX1(0.0,CNO3S(NU(NY,NX),NY,NX)-ZOMN)
      CNO3Y=AMAX1(0.0,CNO3B(NU(NY,NX),NY,NX)-ZOMN)
      RINOOR(N,K,NY,NX)=AMAX1(RINOPR
     2,BIOA*OMA(N,K)*TFNG(N,K)*ZOMX*XNFH) 
     3*(FNO3S*CNO3X/(CNO3X+ZOKU)+FNO3B*CNO3Y
     4/(CNO3Y+ZOKU)) 
      ZNO3M=ZOMN*VOLW(NU(NY,NX),NY,NX) 
      RINO3R(N,K)=AMIN1(FNO3XR(N,K)*AMAX1(0.0
     2,(ZNO3T(NU(NY,NX))-ZNO3M)),RINOOR(N,K,NY,NX))
      ELSE
      RINOOR(N,K,NY,NX)=0.0
      RINO3R(N,K)=RINOPR
      ENDIF
      TRINH4(NY,NX)=TRINH4(NY,NX)+RINO3R(N,K)
C
C     ADDITIONAL MINERALIZATION-IMMOBILIZATION OF H2PO4 IN SURFACE
C     LITTER FROM MICROBIAL C:P AND H2PO4 CONCENTRATION IN BAND AND 
C     NON-BAND SOIL ZONES OF SOIL SURFACE
C
C     RIPOPR=H2PO4 mineralization (-ve) or immobilization (+ve) demand
C        in litter (g P t-1) 
C     NU=surface layer number
C     CH2P4,CH2P4B=H2PO4 concentrations in non-band, band (g P m-3)
C     HPMX,HPMN,HPKU=parameters for max H2PO4 uptake rate (g P m-2 h-1,
C        minimum H2PO4 concentration (g P m-3) and Km for H2PO4 uptake 
C        (g P m-3) 
C     RIPOOR=microbial limitation to H2PO4 demand
C     BIOA=microbial surface area (m2 g C-1)
C     OMA=active biomass (g C)
C     TFNG=temperature+water limitation to biological activity
C     FH2PS,FH2PB=fractions of H2PO4 in non-band, band
C     VOLW=soil water content (m3)
C     H2P4M=H2PO4 not available for uptake (g P) 
C     FPO4XR=fraction of biological H2PO4 demand by microbial
C        population 
C     RIPO4R=H2PO4 supply-limited H2PO4 mineralization-immobilization
C        (g P t-1)
C     TRIPO4=total H2PO4 net mineralization (-ve) or immobilization
C        (+ve) (g P t-1)
C
      RIPOPR=RIPOP-RIPO4(N,K)
      IF(RIPOPR.GT.0.0)THEN
      CH2PX=AMAX1(0.0,CH2P4(NU(NY,NX),NY,NX)-HPMN)
      CH2PY=AMAX1(0.0,CH2P4B(NU(NY,NX),NY,NX)-HPMN)
      RIPOOR(N,K,NY,NX)=AMIN1(RIPOPR
     2,BIOA*OMA(N,K)*TFNG(N,K)*HPMX*XNFH) 
     3*(FH2PS*CH2PX/(CH2PX+HPKU)+FH2PB*CH2PY
     4/(CH2PY+HPKU)) 
      H2P4M=HPMN*VOLW(NU(NY,NX),NY,NX) 
      RIPO4R(N,K)=AMIN1(FPO4XR(N,K)*AMAX1(0.0
     2,(H2P4T(NU(NY,NX))-H2P4M)),RIPOOR(N,K,NY,NX))
      ELSE
      RIPOOR(N,K,NY,NX)=0.0
      RIPO4R(N,K)=RIPOPR
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+RIPO4R(N,K) 
C     WRITE(*,7778)'RIPO4R',I,J,NX,NY,L,K,N,RIPO4R(N,K),FPO4XR(N,K)
C    2,H2P4T(NU(NY,NX)),RIPOOR(N,K,NY,NX),RIPOPR 
C
C     ADDITIONAL MINERALIZATION-IMMOBILIZATION OF HPO4 IN SURFACE
C     LITTER FROM MICROBIAL C:P AND HPO4 CONCENTRATION IN BAND AND 
C     NON-BAND SOIL ZONES OF SOIL SURFACE
C
C     RIP1PR=HPO4 mineralization (-ve) or immobilization (+ve) demand
C        in litter (g P t-1) 
C     RIPO4R=H2PO4 supply-limited H2PO4 mineralization-immobilization
C        (g P t-1)
C     NU=surface layer number
C     CH1P4,CH1P4B=HPO4 concentrations in non-band, band (g P m-3)
C     HPMX,HPMN,HPKU=parameters for max HPO4 uptake rate (g P m-2 h-1,
C        minimum HPO4 concentration (g P m-3) and Km for HPO4 uptake 
C        (g P m-3), assumed same as H2PO4
C     RIPO1R=microbially limited HPO4 demand
C     BIOA=microbial surface area (m2 g C-1)
C     OMA=active biomass (g C)
C     TFNG=temperature+water limitation to biological activity
C     FH1PS,FH1PB=fractions of HPO4 in non-band, band
C     RIPO1R=HPO4 supply-unlimited HPO4 mineralization-immobilization
C        (g P t-1) 
C     VOLW=soil water content (m3)
C     H1P4M=HPO4 not available for uptake (g P)
C     FP14XR=fraction of biological HPO4 demand
C     RIP14R=HPO4 supply-limited HPO4 minereralation-immobilization
C        (g P t-1)
C     TRIPO4=total HPO4 net mineralization (-ve) or immobilization
C        (+ve) (g P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      RIP1PR=0.1*AMAX1(0.0,RIPOPR-RIPO4R(N,K))
      IF(RIP1PR.GT.0.0)THEN
      CH1PX=AMAX1(0.0,CH1P4(NU(NY,NX),NY,NX)-HPMN)
      CH1PY=AMAX1(0.0,CH1P4B(NU(NY,NX),NY,NX)-HPMN)
      RIPO1R(N,K,NY,NX)=AMIN1(RIP1PR
     2,BIOA*OMA(N,K)*TFNG(N,K)*HPMX*XNFH) 
     3*(FH1PS*CH1PX/(CH1PX+HPKU)+FH1PB*CH1PY
     4/(CH1PY+HPKU)) 
      H1P4M=HPMN*VOLW(NU(NY,NX),NY,NX) 
      RIP14R(N,K)=AMIN1(FP14XR(N,K)*AMAX1(0.0
     2,(H1P4T(NU(NY,NX))-H1P4M)),RIPO1R(N,K,NY,NX))
      ELSE
      RIPO1R(N,K,NY,NX)=0.0
      RIP14R(N,K)=RIP1PR
      ENDIF
      TRIPO4(NY,NX)=TRIPO4(NY,NX)+RIP14R(N,K) 
C     WRITE(*,7778)'RIP14R',I,J,NX,NY,L,K,N,RIP14R(N,K),FP14XR(N,K)
C    2,H1P4T(NU(NY,NX)),RIPO1R(N,K,NY,NX),RIP1PR 
      ENDIF
C
C     pH EFFECT ON MAINTENANCE RESPIRATION
C
C     FPH=pH effect on maintenance respiration
C     RMOM=specific maintenance respiration rate (g C g N-1 t-1)
C     TFNR=temperature effect on maintenance respiration
C     RMOMC=maintenance respiration rate of labile(1) and resistant(2)
C        fractions (g C t-1)
C     OMN=microbial N biomass in labile(1) and resistant(2)
C        fractions (g N)
C     RMOMK=effect of low microbial C concentration on maintenance
C        respiration
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      FPH=1.0+AMAX1(0.0,0.25*(6.5-PH(L,NY,NX)))
      RMOMX=RMOM*TFNR(N,K)*FPH*XNFH
      RMOMC(1,N,K)=OMN(1,N,K,L,NY,NX)*RMOMX*RMOMK(1) 
      RMOMC(2,N,K)=OMN2(N,K)*RMOMX*RMOMK(2)
C
C     MICROBIAL MAINTENANCE AND GROWTH RESPIRATION
C
C     RMOMT=total maintenance respiration (g C t-1)
C     RGOMT=growth respiration (g C t-1)
C     RGOMO=O2-limited respiration (g C t-1)
C     RXOMT=senescence respiration (g C t-1)
C
      RMOMT=RMOMC(1,N,K)+RMOMC(2,N,K)
      RGOMT=AMAX1(0.0,RGOMO(N,K)-RMOMT)
      RXOMT=AMAX1(0.0,RMOMT-RGOMO(N,K))
C
C     NON-SYMBIOTIC N2 FIXATION FROM GROWTH RESPIRATION,  
C     FIXATION ENERGY REQUIREMENT, MICROBIAL N REQUIREMENT 
C        N=6:aerobic
C         =7:anaerobic
C
C     RGN2P=respiration required to meet N2 fixation demand (g C t-1)
C     OMC(3,OMN(3=microbial nonstructural C,N (g C,N)
C     CNOMC=maximum microbial N:C ratio from ‘starts.f’(g N g C-1)
C     EN2F=N2 fixation yield per unit nonstructural C (g N g C-1)
C     RGOMT=growth respiration (g C t-1)
C     RGN2F=respiration for N2 fixation (g C t-1)
C     CZ2GS=aqueous N2 concentration (g N m-3)
C     ZFKM=Km for N2 uptake (g N m-3)
C     OMGR=rate constant for transferring nonstructural to 
C       structural microbial C (h-1)
C     OMGR*OMC(3,N,K,L,NY,NX)*XNFH=nonstructural C limitation to RGN2F
C        (g C t-1)
C     RN2FX=N2 fixation rate (g N t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      IF(K.LE.4.AND.(N.EQ.6.OR.N.EQ.7))THEN
      RGN2P=AMAX1(0.0,OMC(3,N,K,L,NY,NX)*CNOMC(3,N,K)
     2-OMN(3,N,K,L,NY,NX))/EN2F(N)
      IF(RGOMT.GT.ZEROS(NY,NX))THEN
      RGN2F(N,K)=AMIN1(RGOMT*RGN2P/(RGOMT+RGN2P)
     2*CZ2GS(L,NY,NX)/(CZ2GS(L,NY,NX)+ZFKM)
     2,OMGR*OMC(3,N,K,L,NY,NX)*XNFH)
      ELSE
      RGN2F(N,K)=0.0
      ENDIF
      RN2FX(N,K)=RGN2F(N,K)*EN2F(N) 
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,5566)'N2 FIX',I,J,NX,NY,L,K,N,RN2FX(N,K),EN2F(N) 
C    2,OMC(3,N,K,L,NY,NX)*CNOMC(3,N,K),OMN(3,N,K,L,NY,NX)
C    3,RINH4(N,K),RINO3(N,K),RGN2P,RGN2F(N,K),RGOMT
C    4,CZ2GS(L,NY,NX) 
5566  FORMAT(A8,7I4,30E12.4)
C     ENDIF
      ELSE
      RN2FX(N,K)=0.0
      RGN2F(N,K)=0.0
      ENDIF
C
C     DOC, DON, DOP AND ACETATE UPTAKE DRIVEN BY GROWTH RESPIRATION
C     FROM O2, NOX AND C REDUCTION
C
C     CGOMX=DOC+acetate uptake from aerobic growth respiration 
C        (g C t-1) 
C     CGOMD=DOC+acetate uptake from denitrifier growth respiration
C        (g C t-1) 
C     CGOQC,CGOAC=DOC,acetate uptake from aerobic growth respiration 
C        (g C t-1) 
C     RMOMT=maintenance respiration (g C t-1) 
C     RGOMO=total aerobic respiration (g C t-1) 
C     RGOMD=respiration for denitrification (g C t-1) 
C     RGN2F=respiration for N2 fixation (g C t-1) 
C     ECHZ,ENOX=growth respiration efficiencies for O2, NOx reduction
C        (g C g C-1)
C     CGOMC,CGOQC,CGOAC=total DOC+acetate, DOC, acetate uptake
C        (heterotrophs) (g C t-1)
C     CGOMC=total CO2,CH4 uptake (autotrophs) (g C t-1)
C     CGOMN,CGOMP=DON, DOP uptake (g N,P t-1)
C     FGOCP,FGOAP=DOC,acetate/(DOC+acetate)
C     OQN,OPQ=DON,DOP (g N,P)
C     FOMK=faction of OMA in total OMA
C     CNQ,CPQ=DON/DOC, DOP/DOC
C     FCN,FCP=limitation to respiration from N,P
C
      CGOMX=AMIN1(RMOMT,RGOMO(N,K))+RGN2F(N,K)
     2+(RGOMT-RGN2F(N,K))/ECHZ
      CGOMD=RGOMD(N,K)/ENOX
      CGOMC(N,K)=CGOMX+CGOMD
      IF(K.LE.4)THEN
      CGOQC(N,K)=CGOMX*FGOCP+CGOMD
      CGOAC(N,K)=CGOMX*FGOAP
      CGOXC=CGOQC(N,K)+CGOAC(N,K)
      CGOMN(N,K)=AMAX1(0.0,AMIN1(OQN(K,L,NY,NX)*FOMK(N,K)
     2,CGOXC*CNQ(K)/FCN(N,K)))
      CGOMP(N,K)=AMAX1(0.0,AMIN1(OQP(K,L,NY,NX)*FOMK(N,K)
     2,CGOXC*CPQ(K)/FCP(N,K)))
C     IF(I.GT.140.AND.L.LE.3.)THEN
C     WRITE(*,5557)'CGOQC',I,J,NFZ,NX,NY,L,K,N,CGOQC(N,K),CGOMX 
C    2,FGOCP,FGOAP,CGOMD,RMOMT,RGN2F(N,K),ECHZ 
C    3,RGOMD(N,K),ENOX,RGOMO(N,K),WFN(N,K),FOXYX
C     WRITE(*,5557)'CGOMP',I,J,NFZ,NX,NY,L,K,N,CGOMP(N,K)
C    2,OQP(K,L,NY,NX),FOMK(N,K),CGOXC,CPQ(K),FCP(N,K)
C    2,CGOQC(N,K),CGOAC(N,K),CGOMX,RMOMT,RGOMO(N,K),RGN2F(N,K)
C    2,RGOMT,RGN2F(N,K),ECHZ,FGOCP,FGOAP,CGOMD,RGOMP,WFN(N,K) 
5557  FORMAT(A8,8I4,30E12.4)
C     ENDIF
      ELSE
      CGOQC(N,K)=CGOMX+CGOMD
      CGOAC(N,K)=0.0
      CGOMN(N,K)=0.0
      CGOMP(N,K)=0.0
      ENDIF
C
C     RECYCLE C,N,P FROM SENESCING BIOMASS TO NONSTRUCTURAL STORAGE
C
C     OMC(3,OMN(3,OMP(3=nonstructural C,N,P (g C,N,P)
C     CCC,CNC,CPC=C:N:P ratios used to calculate C,N,P recycling
C     CNKI,CPKI=nonstructural N,P inhibition constant on microbial
C        nutrient recycling (g N,P g-1 C)
C     RCCC,RCCN,RCCP=C,N,P recycling fractions
C     RCCZ,RCCY=min, max C recycling fractions
C     RCCX,RCCQ=max N,P recycling fractions
C
      IF(OMC(3,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CCC=AMAX1(0.0,AMIN1(1.0
     2,OMN(3,N,K,L,NY,NX)/(OMN(3,N,K,L,NY,NX)
     2+OMC(3,N,K,L,NY,NX)*CNKI)
     3,OMP(3,N,K,L,NY,NX)/(OMP(3,N,K,L,NY,NX)
     4+OMC(3,N,K,L,NY,NX)*CNKI)))
      CNC=AMAX1(0.0,AMIN1(1.0
     2,OMC(3,N,K,L,NY,NX)/(OMC(3,N,K,L,NY,NX)
     2+OMN(3,N,K,L,NY,NX)/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0
     2,OMC(3,N,K,L,NY,NX)/(OMC(3,N,K,L,NY,NX)
     3+OMP(3,N,K,L,NY,NX)/CPKI)))
      RCCC=RCCZ+CCC*RCCY
      RCCN=CNC*RCCX
      RCCP=CPC*RCCQ
      ELSE
      RCCC=0.0
      RCCN=0.0
      RCCP=0.0
      ENDIF
C     IF((I/120)*120.EQ.I.AND.J.EQ.24)THEN
C     WRITE(*,5555)'RCCC',I,J,NX,NY,L,K,N,RCCC,RCCN,RCCP
C    2,OMC(3,N,K,L,NY,NX),OMN(3,N,K,L,NY,NX),OMP(3,N,K,L,NY,NX)
C    3,CCC,CNC,CPC
C     ENDIF
C
C     MICROBIAL ASSIMILATION OF NONSTRUCTURAL C,N,P
C
C     CGOMZ=transfer from nonstructural to structural microbial C 
C        (g C t-1)
C     TFNG=temperature+water limitation to biological activity
C     OMGR=rate constant for transferring nonstructural to 
C       structural microbial C (h-1)
C     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural
C        C,N,P (g C,N,P t-1)
C     FL=partitioning between labile and resistant microbial
C        components
C     OMC(3,OMN(3,OMP(3=nonstructural C,N,P (g C,N,P)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      CGOMZ=TFNG(N,K)*OMGR*AMAX1(0.0,OMC(3,N,K,L,NY,NX))*XNFH
      DO 745 M=1,2
      CGOMS(M,N,K)=FL(M)*CGOMZ 
      IF(OMC(3,N,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CGONS(M,N,K)=AMIN1(FL(M)*AMAX1(0.0,OMN(3,N,K,L,NY,NX))
     2,CGOMS(M,N,K)*OMN(3,N,K,L,NY,NX)/OMC(3,N,K,L,NY,NX))
      CGOPS(M,N,K)=AMIN1(FL(M)*AMAX1(0.0,OMP(3,N,K,L,NY,NX))
     2,CGOMS(M,N,K)*OMP(3,N,K,L,NY,NX)/OMC(3,N,K,L,NY,NX))
      ELSE
      CGONS(M,N,K)=0.0
      CGOPS(M,N,K)=0.0
      ENDIF 
C
C     MICROBIAL DECOMPOSITION FROM BIOMASS, SPECIFIC DECOMPOSITION
C     RATE, TEMPERATURE
C
C     SPOMX=rate constant for microbial decomposition (t-1)
C     SPOMC=basal decomposition rate (h-1)
C     SPOMK=effect of microbial C concentration on microbial decay
C     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition (g C,N,P t-1)
C     RDOMC,RDOMN,RDOMP=microbial C,N,P litterfall to microbial
C        residue,humus (g C,N,P t-1)
C     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling to nonstructural
C        C,N,P (g C,N,P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      SPOMX=SQRT(TFNG(N,K))*SPOMC(M)*SPOMK(M)*XNFH
      RXOMC(M,N,K)=AMAX1(0.0,OMC(M,N,K,L,NY,NX)*SPOMX)
      RXOMN(M,N,K)=AMAX1(0.0,OMN(M,N,K,L,NY,NX)*SPOMX)
      RXOMP(M,N,K)=AMAX1(0.0,OMP(M,N,K,L,NY,NX)*SPOMX)
      RDOMC(M,N,K)=RXOMC(M,N,K)*(1.0-RCCC)
      RDOMN(M,N,K)=RXOMN(M,N,K)*(1.0-RCCC)*(1.0-RCCN)
      RDOMP(M,N,K)=RXOMP(M,N,K)*(1.0-RCCC)*(1.0-RCCP)
      R3OMC(M,N,K)=RXOMC(M,N,K)-RDOMC(M,N,K)
      R3OMN(M,N,K)=RXOMN(M,N,K)-RDOMN(M,N,K)
      R3OMP(M,N,K)=RXOMP(M,N,K)-RDOMP(M,N,K)
C
C     HUMIFICATION OF MICROBIAL DECOMPOSITION PRODUCTS FROM
C     DECOMPOSITION RATE, SOIL CLAY AND OC 'EHUM' FROM 'HOUR1'
C
C     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall 
C        to humus (g C,N,P t-1)
C     RDOMC,RDOMN,RDOMP=microbial C,N,P litterfall to microbial
C        residue,humus (g C,N,P t-1)
C     EHUM=fraction of microbial decomposition product allocated to
C        humus from ‘hour1.f’
C     CNRH,CPRH=default N:C,P:C ratios in SOC complexes
C        woody(0),fine(1),manure(2),POC(3),humus (4) from ‘starts.f’
C
      RHOMC(M,N,K)=AMAX1(0.0,AMIN1(RDOMC(M,N,K)*EHUM(L,NY,NX)
     2,RDOMN(M,N,K)/CNRH(4),RDOMP(M,N,K)/CPRH(4)))
      RHOMN(M,N,K)=AMAX1(0.0,AMIN1(RDOMN(M,N,K)*EHUM(L,NY,NX)
     2,RHOMC(M,N,K)*CNRH(4)))
      RHOMP(M,N,K)=AMAX1(0.0,AMIN1(RDOMP(M,N,K)*EHUM(L,NY,NX)
     2,RHOMC(M,N,K)*CPRH(4))) 
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,8821)'RHOMC',I,J,L,K,N,M
C    2,RXOMC(M,N,K),RXOMN(M,N,K),RXOMP(M,N,K)
C    2,RDOMC(M,N,K),RDOMN(M,N,K),RDOMP(M,N,K)
C    2,R3OMC(M,N,K),R3OMN(M,N,K),R3OMP(M,N,K)
C    2,RHOMC(M,N,K),RHOMN(M,N,K),RHOMP(M,N,K) 
C    4,OMC(M,N,K,L,NY,NX),OMN(M,N,K,L,NY,NX)
C    5,OMP(M,N,K,L,NY,NX)
C    4,OMC(3,N,K,L,NY,NX),OMN(3,N,K,L,NY,NX)
C    5,OMP(3,N,K,L,NY,NX)
C    6,OQC(K,L,NY,NX),OQN(K,L,NY,NX),OQP(K,L,NY,NX)
C    2,SPOMX,RCCC,RCCN,RCCP
C     ENDIF
C
C     NON-HUMIFIED PRODUCTS TO MICROBIAL RESIDUE
C
C     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall 
C        to microbial residue (g C,N,P t-1)
C     RHOMC,RHOMN,RHOMP=transfer of microbial C,N,P litterfall 
C        to humus (g C,N,P t-1)
C     RDOMC,RDOMN,RDOMP=microbial C,N,P litterfall to microbial
C        residue,humus (g C,N,P t-1)
C
      RCOMC(M,N,K)=RDOMC(M,N,K)-RHOMC(M,N,K)
      RCOMN(M,N,K)=RDOMN(M,N,K)-RHOMN(M,N,K)
      RCOMP(M,N,K)=RDOMP(M,N,K)-RHOMP(M,N,K)
745   CONTINUE
C
C     MICROBIAL DECOMPOSITION WHEN C MAINTENANCE RESPIRATION
C     EXCEEDS C UPTAKE
C
C     OMC,OMN,OMP=microbial C,N,P (g C,N,P)
C     RMOMT=total maintenance respiration (g C t-1)
C     RXOMT=senescence respiration (g C t-1)
C     RCCC=C recycling fraction
C     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence 
C        (g C,N,P t-1)
C     RMOMC=maintenance respiration (g C t-1)
C     CNOMA,CPOMA=N,P concentrations of active biomass (g N,P g C-1)
C     RDMMC,RDMMN,RDMMP=microbial C,N,P litterfall from senescence
C        to humus, microbial residue (g C,N,P t-1)
C     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence
C        to nonstructural C,N,P (g C,N,P t-1)
C
      IF(RXOMT.GT.ZEROS(NY,NX).AND.RMOMT.GT.ZEROS(NY,NX)
     2.AND.RCCC.GT.ZERO)THEN
      FRM=RXOMT/RMOMT
      DO 730 M=1,2
      RXMMC(M,N,K)=AMIN1(OMC(M,N,K,L,NY,NX)
     2,AMAX1(0.0,FRM*RMOMC(M,N,K)/RCCC))
      RXMMN(M,N,K)=AMIN1(OMN(M,N,K,L,NY,NX)
     2,AMAX1(0.0,RXMMC(M,N,K)*CNOMA(N,K))) 
      RXMMP(M,N,K)=AMIN1(OMP(M,N,K,L,NY,NX)
     2,AMAX1(0.0,RXMMC(M,N,K)*CPOMA(N,K)))
      RDMMC(M,N,K)=RXMMC(M,N,K)*(1.0-RCCC)
      RDMMN(M,N,K)=RXMMN(M,N,K)*(1.0-RCCN)*(1.0-RCCC)
      RDMMP(M,N,K)=RXMMP(M,N,K)*(1.0-RCCP)*(1.0-RCCC)
      R3MMC(M,N,K)=RXMMC(M,N,K)-RDMMC(M,N,K)
      R3MMN(M,N,K)=RXMMN(M,N,K)-RDMMN(M,N,K)
      R3MMP(M,N,K)=RXMMP(M,N,K)-RDMMP(M,N,K)
C
C     HUMIFICATION AND RECYCLING OF RESPIRATION DECOMPOSITION
C     PRODUCTS
C
C     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P 
C        to humus (g C,N,P t-1)
C     EHUM=fraction of microbial decomposition product allocated to
C        humus from ‘hour1.f’
C     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P 
C        to microbial residue (g C,N,P t-1)
C
      RHMMC(M,N,K)=AMAX1(0.0,RDMMC(M,N,K)*EHUM(L,NY,NX))
      RHMMN(M,N,K)=AMAX1(0.0,AMIN1(RDMMN(M,N,K)*EHUM(L,NY,NX)
     2,RDMMC(M,N,K)*CNRH(4)))
      RHMMP(M,N,K)=AMAX1(0.0,AMIN1(RDMMP(M,N,K)*EHUM(L,NY,NX)
     2,RDMMC(M,N,K)*CPRH(4))) 
      RCMMC(M,N,K)=RDMMC(M,N,K)-RHMMC(M,N,K)
      RCMMN(M,N,K)=RDMMN(M,N,K)-RHMMN(M,N,K)
      RCMMP(M,N,K)=RDMMP(M,N,K)-RHMMP(M,N,K)
C     IF(L.EQ.11.AND.K.EQ.1)THEN
C     WRITE(*,8821)'RCMMC',I,J,L,K,N,M,RCMMC(M,N,K)
C    2,RDMMC(M,N,K),RHMMC(M,N,K),OMC(M,N,K,L,NY,NX)
C    3,FRM,RMOMC(M,N,K),OMN(1,N,K,L,NY,NX),OMN2(N,K)
C    4,RMOM,TFNR(N,K),FPH,RDMMN(M,N,K),CNSHZ,RDMMP(M,N,K)
C    5,CPSHZ,EHUM(L,NY,NX),RXOMT,RMOMT,RMOMT,RGOMO(N,K)
C    6,RGOMP,WFN(N,K)
C     WRITE(*,8821)'RCMMP',I,J,L,K,N,M,RCMMP(M,N,K)
C    2,RDMMP(M,N,K),RHMMP(M,N,K),EHUM(L,NY,NX) 
C    3,RCCC,RCCN,RCCP,RXMMP(M,N,K)
C     ENDIF
730   CONTINUE
      ELSE
      DO 720 M=1,2
      RXMMC(M,N,K)=0.0
      RXMMN(M,N,K)=0.0
      RXMMP(M,N,K)=0.0
      RDMMC(M,N,K)=0.0
      RDMMN(M,N,K)=0.0
      RDMMP(M,N,K)=0.0
      R3MMC(M,N,K)=0.0
      R3MMN(M,N,K)=0.0
      R3MMP(M,N,K)=0.0
      RHMMC(M,N,K)=0.0
      RHMMN(M,N,K)=0.0
      RHMMP(M,N,K)=0.0
      RCMMC(M,N,K)=0.0
      RCMMN(M,N,K)=0.0
      RCMMP(M,N,K)=0.0
720   CONTINUE
      ENDIF
      ELSE
      RUPOX(N,K)=0.0
      RGOMO(N,K)=0.0
      RCO2X(N,K)=0.0
      RCH3X(N,K)=0.0
      RCH4X(N,K)=0.0
      RGOMY(N,K)=0.0
      RGOMD(N,K)=0.0
      CGOMC(N,K)=0.0
      CGOMN(N,K)=0.0
      CGOMP(N,K)=0.0
      CGOQC(N,K)=0.0
      CGOAC(N,K)=0.0
      RDNO3(N,K)=0.0
      RDNOB(N,K)=0.0
      RDNO2(N,K)=0.0
      RDN2B(N,K)=0.0
      RDN2O(N,K)=0.0
      RN2FX(N,K)=0.0
      RINH4(N,K)=0.0
      RINO3(N,K)=0.0
      RIPO4(N,K)=0.0
      RIP14(N,K)=0.0
      RINB4(N,K)=0.0
      RINB3(N,K)=0.0
      RIPOB(N,K)=0.0
      RIP1B(N,K)=0.0
      IF(L.EQ.0)THEN
      RINH4R(N,K)=0.0
      RINO3R(N,K)=0.0
      RIPO4R(N,K)=0.0
      RIP14R(N,K)=0.0
      FNH4XR(N,K)=0.0
      FNO3XR(N,K)=0.0
      FPO4XR(N,K)=0.0
      FP14XR(N,K)=0.0
      ENDIF
      DO 725 M=1,2
      CGOMS(M,N,K)=0.0
      CGONS(M,N,K)=0.0
      CGOPS(M,N,K)=0.0
      RMOMC(M,N,K)=0.0
      RXMMC(M,N,K)=0.0
      RXMMN(M,N,K)=0.0
      RXMMP(M,N,K)=0.0
      RDMMC(M,N,K)=0.0
      RDMMN(M,N,K)=0.0
      RDMMP(M,N,K)=0.0
      R3MMC(M,N,K)=0.0
      R3MMN(M,N,K)=0.0
      R3MMP(M,N,K)=0.0
      RHMMC(M,N,K)=0.0
      RHMMN(M,N,K)=0.0
      RHMMP(M,N,K)=0.0
      RCMMC(M,N,K)=0.0
      RCMMN(M,N,K)=0.0
      RCMMP(M,N,K)=0.0
      RXOMC(M,N,K)=0.0
      RXOMN(M,N,K)=0.0
      RXOMP(M,N,K)=0.0
      RDOMC(M,N,K)=0.0
      RDOMN(M,N,K)=0.0
      RDOMP(M,N,K)=0.0
      R3OMC(M,N,K)=0.0
      R3OMN(M,N,K)=0.0
      R3OMP(M,N,K)=0.0
      RHOMC(M,N,K)=0.0
      RHOMN(M,N,K)=0.0
      RHOMP(M,N,K)=0.0
      RCOMC(M,N,K)=0.0
      RCOMN(M,N,K)=0.0
      RCOMP(M,N,K)=0.0
725   CONTINUE
      RH2GX(N,K)=0.0
      IF(K.EQ.5)THEN
      RVOXA(N)=0.0
      RVOXB(N)=0.0
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
C     FNO2,FNB2=fraction of NO2 demand in non-band,band
C     RNO2Y,RVMXC=total,NO2 chemodenitrification NO2 demand in non-band
C        from previous time step (g N t-1)
C     RN2BY,RVMBC=total,NO2 chemodenitrification NO2 demand in band
C        from previous time step (g N t-1)
C     VMXC4S,VMXC4B=NO2 supply-unlimited NO2 reduction in non-band,band 
C        (g N t-1)
C     CHNO2,CHNOB=nitrous acid concentration in non-band,band (g N m-3)
C     VOLWM=soil water content (m3)
C     FNO3S,FNO3B=fractions of NO2 in non-band,band
C     TFNX=temperature effect on biological activity
C     RCNO2,RCNOB=NO2 supply-limited nitrous acid reduction in 
C        non-band,band (g N t-1)
C     RCN2O,RCN2B=N2O production from nitrous acid reduction in 
C        non-band,band (g N t-1)  
C     RCNO3,RCN3B=NO3 production from nitrous acid reduction in 
C        non-band,band (g N t-1)
C     RCOQN=DON production from nitrous acid reduction (g N t-1) 
C     RVMXC,RVMBC=NO2 supply-unlimited NO2 reduction in non-band,band
C        (g N t-1) 
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
      VMXC4S=7.5E-02*CHNO2*VOLWM(NPH,L,NY,NX)*FNO3S*TFNX*XNFH 
      VMXC4B=7.5E-02*CHNOB*VOLWM(NPH,L,NY,NX)*FNO3B*TFNX*XNFH 
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
C     WRITE(*,7779)'CHEMO',I,J,NFZ,L,NPH,RCNO2,RCNOB,CHY1,CHNO2,CHNOB
C    2,CNO2S(L,NY,NX),CNO2B(L,NY,NX),VOLWM(NPH,L,NY,NX),VOLW(0,NY,NX) 
C    3,FNO2,VMXC4S,VMXC4B,RVMXC(L,NY,NX),RNO2Y(L,NY,NX),RCN2O,RCN2B 
C    4,RCNO3,RCNOB,RCOQN,VLNO3(L,NY,NX),VLNOB(L,NY,NX) 
7779  FORMAT(A8,5I4,30E12.4)
C     ENDIF
C
C     DECOMPOSITION    
C
C     ROQCK=total respiration of DOC+DOA used to represent microbial
C        activity (g C t-1)
C     ROQCD=DOC,DOA-unlimited microbial respiration (g C t-1) 
C
      DO 1870 K=0,KL
      ROQCK(K)=0.0
      DO 1875 N=1,7
      ROQCK(K)=ROQCK(K)+ROQCD(N,K) 
1875  CONTINUE
      XOQCK(K)=0.0
      XOQCZ(K)=0.0
      XOQNZ(K)=0.0
      XOQPZ(K)=0.0
      XOQAZ(K)=0.0
      DO 845 N=1,7
      DO 845 M=1,3
      XOMCZ(M,N,K)=0.0
      XOMNZ(M,N,K)=0.0
      XOMPZ(M,N,K)=0.0
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
C     OSCH=total SOC in each K (g C)
C     XFRK,XFRC,XFRN,XFRP,XFRA=transfer of activity,DOC,DON,DOP,
C        acetate between each K and KK (g C,C,N,P,C t-1) 
C     FPRIM=priming transfer rate constant (h-1)
C     TFND=temperature effect on priming transfers
C     ROQCK=total respiration of DOC+DOA used to represent microbial
C        activity (g C t-1)
C     OQC,OQN,OQP=DOC,DON,DOP (g C,N,P)
C     XOQCK,XOQCZ,XOQNZ,XOQPZ,XOQAZ=total XFRK,XFRC,XFRN,XFRP,XFRA 
C        for all K (g C,C,N,P,C t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      DO 795 K=0,KL
      IF(K.LE.KL-1)THEN
      DO 800 KK=K+1,KL
      OSRT=OSCH(K)+OSCH(KK)      
      IF(OSCH(K).GT.ZEROS(NY,NX).AND.OSCH(KK).GT.ZEROS(NY,NX))THEN
      XFRK=FPRIM*TFND(L,NY,NX)*(ROQCK(K)*OSCH(KK)
     2-ROQCK(KK)*OSCH(K))/OSRT*XNFH 
      XFRC=FPRIM*TFND(L,NY,NX)*(OQC(K,L,NY,NX)*OSCH(KK)
     2-OQC(KK,L,NY,NX)*OSCH(K))/OSRT*XNFH 
      XFRN=FPRIM*TFND(L,NY,NX)*(OQN(K,L,NY,NX)*OSCH(KK)
     2-OQN(KK,L,NY,NX)*OSCH(K))/OSRT*XNFH 
      XFRP=FPRIM*TFND(L,NY,NX)*(OQP(K,L,NY,NX)*OSCH(KK)
     2-OQP(KK,L,NY,NX)*OSCH(K))/OSRT*XNFH 
      XFRA=FPRIM*TFND(L,NY,NX)*(OQA(K,L,NY,NX)*OSCH(KK)
     2-OQA(KK,L,NY,NX)*OSCH(K))/OSRT*XNFH 
      IF(ROQCK(K)+XOQCK(K)-XFRK.GT.0.0
     2.AND.ROQCK(KK)+XOQCK(KK)+XFRK.GT.0.0)THEN 
      XOQCK(K)=XOQCK(K)-XFRK
      XOQCK(KK)=XOQCK(KK)+XFRK
C     IF(I.EQ.116)THEN
C     WRITE(*,4442)'XOQCK',I,J,NX,NY,L,K,KK,XFRC,ROQCK(K)
C    2,OSCH(K),ROQCK(KK),OSCH(KK),XOQCK(K),XOQCK(KK)
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
C    2,OSCH(K),OQC(KK,L,NY,NX),OSCH(KK),XOQCZ(K),XOQCZ(KK)
C    3,OQC(K,L,NY,NX),OQC(KK,8,NY,NX)
C     ENDIF
      ENDIF
      IF(OQN(K,L,NY,NX)+XOQNZ(K)-XFRN.GT.0.0
     2.AND.OQN(KK,L,NY,NX)+XOQNZ(KK)+XFRN.GT.0.0)THEN 
      XOQNZ(K)=XOQNZ(K)-XFRN
      XOQNZ(KK)=XOQNZ(KK)+XFRN
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
C     WRITE(*,4442)'XOQNZ',I,J,NX,NY,L,K,KK,XFRN,OQN(K,L,NY,NX)
C    2,OSCH(K),OQN(KK,L,NY,NX),OSCH(KK),XOQNZ(K),XOQNZ(KK)
C     ENDIF
      ENDIF
      IF(OQP(K,L,NY,NX)+XOQPZ(K)-XFRP.GT.0.0
     2.AND.OQP(KK,L,NY,NX)+XOQPZ(KK)+XFRP.GT.0.0)THEN 
      XOQPZ(K)=XOQPZ(K)-XFRP
      XOQPZ(KK)=XOQPZ(KK)+XFRP
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.L.EQ.4)THEN
C     WRITE(*,4442)'XOQPZ',I,J,NX,NY,L,K,KK,XFRP,OQP(K,L,NY,NX)
C    2,OSCH(K),OQP(KK,L,NY,NX),OSCH(KK),XOQPZ(K),XOQPZ(KK)
C     ENDIF
      ENDIF
      IF(OQA(K,L,NY,NX)+XOQAZ(K)-XFRA.GT.0.0
     2.AND.OQA(KK,L,NY,NX)+XOQAZ(KK)+XFRA.GT.0.0)THEN 
      XOQAZ(K)=XOQAZ(K)-XFRA
      XOQAZ(KK)=XOQAZ(KK)+XFRA
C     IF((I/1)*1.EQ.I.AND.L.EQ.3.AND.K.EQ.1)THEN
C     WRITE(*,4442)'XOQAZ',I,J,NX,NY,L,K,KK,XFRA,OQA(K,L,NY,NX)
C    2,OSCH(K),OQA(KK,L,NY,NX),OSCH(KK),XOQAZ(K),XOQAZ(KK)
C     ENDIF
      ENDIF
C
C     PRIMING of MICROBIAL C,N,P BETWEEN LITTER AND NON-LITTER C
C
C     XFMC,XFMN,XFMP=transfer of microbial C,N,P between each K and KK
C        (g C t-1) 
C     FPRIMM=priming transfer rate constant (h-1)
C     TFNG=temperature+water limitation to biological activity
C     OMC,OMN,OMP=microbial C,N,P (g C,N,P)
C     OSCH=total SOC in each K (g C)
C     XOMCZ,XOMNZ,XOMPZ=total microbial C,N,P transfer for all K 
C        (g C t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      DO 850 N=1,7
      DO 850 M=1,3
      XFMC=FPRIMM*TFNG(N,K)*(OMC(M,N,K,L,NY,NX)*OSCH(KK)
     2-OMC(M,N,KK,L,NY,NX)*OSCH(K))/OSRT*XNFH 
      XFMN=FPRIMM*TFNG(N,K)*(OMN(M,N,K,L,NY,NX)*OSCH(KK)
     2-OMN(M,N,KK,L,NY,NX)*OSCH(K))/OSRT*XNFH 
      XFMP=FPRIMM*TFNG(N,K)*(OMP(M,N,K,L,NY,NX)*OSCH(KK)
     2-OMP(M,N,KK,L,NY,NX)*OSCH(K))/OSRT*XNFH 
      IF(OMC(M,N,K,L,NY,NX)+XOMCZ(M,N,K)-XFMC.GT.0.0
     2.AND.OMC(M,N,KK,L,NY,NX)+XOMCZ(M,N,KK)+XFMC.GT.0.0)THEN 
      XOMCZ(M,N,K)=XOMCZ(M,N,K)-XFMC
      XOMCZ(M,N,KK)=XOMCZ(M,N,KK)+XFMC
C     IF((I/10)*10.EQ.I.AND.J.EQ.15.AND.NFZ.EQ.1)THEN
C     WRITE(*,4447)'XOMCZ',I,J,NFZ,NX,NY,L,K,KK,N,M,XFMC
C    2,OMC(M,N,K,L,NY,NX),OMC(M,N,KK,L,NY,NX)
C    2,OQC(K,L,NY,NX),OQC(KK,L,NY,NX) 
C    3,XOMCZ(M,N,K),XOMCZ(M,N,KK) 
4447  FORMAT(A8,10I4,20E12.4)
C     ENDIF
      ENDIF
      IF(OMN(M,N,K,L,NY,NX)+XOMNZ(M,N,K)-XFMN.GT.0.0
     2.AND.OMN(M,N,KK,L,NY,NX)+XOMNZ(M,N,KK)+XFMN.GT.0.0)THEN 
      XOMNZ(M,N,K)=XOMNZ(M,N,K)-XFMN
      XOMNZ(M,N,KK)=XOMNZ(M,N,KK)+XFMN
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,4447)'XOMNZ',I,J,NX,NY,L,K,KK,N,M,XFMN
C    2,OMN(M,N,K,L,NY,NX)
C    2,OSCH(K),OMN(M,N,KK,L,NY,NX),OSCH(KK),XOMNZ(M,N,K),XOMNZ(M,N,KK)
C     ENDIF
      ENDIF
      IF(OMP(M,N,K,L,NY,NX)+XOMPZ(M,N,K)-XFMP.GT.0.0
     2.AND.OMP(M,N,KK,L,NY,NX)+XOMPZ(M,N,KK)+XFMP.GT.0.0)THEN 
      XOMPZ(M,N,K)=XOMPZ(M,N,K)-XFMP
      XOMPZ(M,N,KK)=XOMPZ(M,N,KK)+XFMP
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,4447)'XOMPZ',I,J,NX,NY,L,K,KK,N,M,XFMP
C    2,OMP(M,N,K,L,NY,NX),OSCH(K),OMP(M,N,KK,L,NY,NX),OSCH(KK)
C    2,XOMPZ(M,N,K),XOMPZ(M,N,KK)
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
C     ROQCK=total respiration of DOC+DOA used to represent microbial
C        activity (g C t-1)
C     TOQCK=total respiration of DOC+DOA in soil layer (g C t-1)
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores (g C,N,P,C)
C     OMC,OMN,OMP=microbial C,N,P (g C,N,P)
C     XOQCK,XOQCZ,XOQNZ,XOQPZ,XOQAZ=total priming transfer of
C        activity,DOC,DON,DOP,acetate (g C,C,N,P,C t-1)
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
      DO 840 M=1,3
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)+XOMCZ(M,N,K)
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)+XOMNZ(M,N,K)
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)+XOMPZ(M,N,K)
C     IF((I/10)*10.EQ.I.AND.J.EQ.15.AND.NFZ.EQ.1)THEN
C     WRITE(*,5559)'XOM',I,J,NFZ,NX,NY,L,K,M,N,OMC(M,N,K,L,NY,NX)
C    2,OMN(M,N,K,L,NY,NX),OMP(M,N,K,L,NY,NX)
C    3,XOMCZ(M,N,K),XOMNZ(M,N,K),XOMPZ(M,N,K)
5559  FORMAT(A8,9I4,12E12.4)
C     ENDIF 
840   CONTINUE
C     IF((I/10)*10.EQ.I.AND.J.EQ.15.AND.NFZ.EQ.1)THEN
C     WRITE(*,4443)'PRIM2',I,J,NFZ,NX,NY,L,K,ROQCK(K)
C    2,XOQCK(K),OQC(K,L,NY,NX),XOQCZ(K),OQN(K,L,NY,NX),XOQNZ(K)
C    3,OQP(K,L,NY,NX),XOQPZ(K),OQA(K,L,NY,NX),XOQAZ(K),TOMK(K)
C    3,TONK(K),TOPK(K),TONX(K),TOPX(K),FCNK(K),FCPK(K)
C    4,TOQCK(L,NY,NX)
4443  FORMAT(A8,7I4,20E12.4)
C     ENDIF
C
C     DECOMPOSITION OF ORGANIC SUBSTRATES
C
C     TONK,TOPK=total active biomass N,P (g N,P)
C     TONX,TOPX=maximum total active biomass N,P (g N,P)
C     FCNK,FCPK=N,P limitation to microbial activity
C     CNOMK,CPOMK=N:C,P:C ratios relative to set maximum values 
C
      IF(TOMK(K).GT.ZEROS(NY,NX))THEN
      CNOMK=TONK(K)/TONX(K)
      CPOMK=TOPK(K)/TOPX(K)
      FCNK(K)=AMIN1(1.0,AMAX1(0.50,CNOMK)) 
      FCPK(K)=AMIN1(1.0,AMAX1(0.50,CPOMK))
      ELSE
      FCNK(K)=1.0
      FCPK(K)=1.0
      ENDIF
C
C     AQUEOUS CONCENTRATION OF BIOMASS TO CACULATE INHIBITION
C     CONSTANT FOR DECOMPOSITION
C
C     OSAH=total colonized SOC (g C)
C     COQCK=aqueous concentration of microbial activity (g C m-3 t-1)
C     ROQCK=total respiration of DOC+DOA used to represent microbial
C        activity (g C t-1)
C     VOLWY=biologically active water volume (m3)
C     DCKD=Km for decomposition of SOC at current COQCK (g C m-3)
C     DCKM=Km for decomposition of SOC at zero COQCK (g C m-3)
C     DCKI=inhibition of decomposition by COQCK (g C m-3 t-1)
C     COSC=concentration of total SOC (g C Mg-1)
C     BKVL=mass of soil layer (Mg)
C     DFNS=effect of microbial concentration on decomposition
C     OQCI=DOC product inhibition for decomposition 
C     OQKI=DOC product inhibition constant for decomposition (g C m-3) 
C
      IF(OSAH(K).GT.ZEROS(NY,NX))THEN
      IF(VOLWY.GT.ZEROS2(NY,NX))THEN
      COQCK=AMIN1(0.1E+06,ROQCK(K)/VOLWY)
      ELSE
      COQCK=0.1E+06
      ENDIF
      DCKD=DCKM*(1.0+COQCK/DCKI)
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      COSC=OSAH(K)/BKVL(L,NY,NX)
      ELSE
      COSC=OSAH(K)/VOLY(L,NY,NX)
      ENDIF
      DFNS=COSC/(COSC+DCKD)
      OQCI=1.0/(1.0+COQC(K,L,NY,NX)/OQKI)
C     IF((I/10)*10.EQ.I.AND.J.EQ.15.AND.NFZ.EQ.1)THEN
C     WRITE(*,4242)'COSC',I,J,NFZ,NX,NY,L,K
C    2,DFNS,COSC,COQCK,DCKD,OSAH(K)
C    2,OSAT(K),OSCT(K),ORCT(K),OHC(K,L,NY,NX),BKVL(L,NY,NX),ROQCK(K)
C    3,VOLWY,VOLWRX(NY,NX),VOLW(0,NY,NX),THETY(L,NY,NX) 
4242  FORMAT(A8,7I4,30E12.4)
C     ENDIF
C
C     C, N, P DECOMPOSITION RATE OF SOLID SUBSTRATES FROM
C     RATE CONSTANT, TOTAL ACTIVE BIOMASS, DENSITY FACTOR,
C     TEMPERATURE, SUBSTRATE C:N, C:P
C
C     OSA,OSN,OSP=active biomass C,N,P (g C,N,P)
C     CNS,CPS=N:C,P:C ratios of SOC (g N,P g C-1)
C     RDOSC,RDOSN,RDOSP=decomposition rates of SOC,SON,SOP 
C        (g C,N,P t-1)
C     SPOSC=specific decomposition rate constant for SOC 
C        (g C g C-1 h-1)
C     ROQCK=total respiration of DOC+DOA used to represent microbial
C        activity (g C t-1)
C     DFNS=effect of microbial concentration on decomposition
C     OQCI=DOC product inhibition for decomposition
C     TFNX=temperature effect on biological activity
C     OSAH=total colonized SOC in each K (g C)
C     FCNK,FCPK=N,P limitation to microbial activity in each K
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      DO 785 M=1,4
      IF(OSC(M,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNS(M,K)=AMAX1(0.0,OSN(M,K,L,NY,NX)/OSC(M,K,L,NY,NX))
      CPS(M,K)=AMAX1(0.0,OSP(M,K,L,NY,NX)/OSC(M,K,L,NY,NX))
      RDOSC(M,K)=AMAX1(0.0,AMIN1(OSA(M,K,L,NY,NX)*XNFH
     2,SPOSC(M,K)*ROQCK(K)*DFNS*OQCI*TFNX
     3*OSA(M,K,L,NY,NX)/OSAH(K)))
C    4*AMIN1(FCNK(K),FCPK(K)) 
      RDOSN(M,K)=AMAX1(0.0,AMIN1(OSN(M,K,L,NY,NX)*XNFH
     2,CNS(M,K)*RDOSC(M,K)))
     3/FCNK(K)
      RDOSP(M,K)=AMAX1(0.0,AMIN1(OSP(M,K,L,NY,NX)*XNFH
     2,CPS(M,K)*RDOSC(M,K)))
     3/FCPK(K)
C     IF((I/10)*10.EQ.I.AND.J.EQ.15.AND.NFZ.EQ.1.AND.L.EQ.0)THEN
C     WRITE(*,4444)'RDOSC',I,J,NFZ,NX,NY,L,K,M,RDOSC(M,K),RDOSN(M,K) 
C    2,RDOSP(M,K),CNS(M,K),CPS(M,K),SPOSC(M,K),ROQCK(K),DFNS,TFNX 
C    3,OQCI,OSA(M,K,L,NY,NX),OSAH(K),COSC,COQCK,DCKD,VOLWY
C    4,COXYS(L,NY,NX),TKS(L,NY,NX),PSISM(L,NY,NX),THETW(L,NY,NX)
C    5,WFN(1,K),WFN(3,K),OXYI,VOLW(0,NY,NX),VOLWRX(NY,NX) 
C    4,FOSAH(K,L,NY,NX),VOLY(L,NY,NX),ORGC(L,NY,NX),OSC(M,K,L,NY,NX)
C    2,OSN(M,K,L,NY,NX),OSP(M,K,L,NY,NX),TONK(K),TONX(K),FCNK(K)
C    6,FCPK(K),WFN(1,K),WFN(3,K),WFN(6,K),COQC(K,L,NY,NX)
C    7,OQC(K,L,NY,NX),BKVL(L,NY,NX),OSC(M,0,L,NY,NX)
4444  FORMAT(A8,8I4,50E14.6)
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
C     CH2O AND CELLULOSE WITH REMAINDER TO DOC,DON,DOP
C
C     RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
C        (g C,N,P h-1)
C     RDOSC,RDOSN,RDOSP=decomposition of SOC,SON,SOP (g C,N,P h-1)
C     CNRH,CPRH=N:C,P:C in POC (g N,P g C-1)
C     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
C        (g C,N,P h-1)
C
      IF(K.LE.2)THEN
      RHOSC(4,K)=AMAX1(0.0,AMIN1(RDOSN(4,K)/CNRH(3)
     2,RDOSP(4,K)/CPRH(3),RDOSC(4,K)))
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
C     ORC,ORN,ORP=microbial residue C,N,P (g C,N,P)
C     CNR,CPR=N:C,P:C ratios of microbial residue (g N,P g C-1)
C     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P 
C        (g C,N,P t-1)
C     SPORC=specific decomposition rate constant for microbial residue
C        (g C g C-1 h-1)
C     ROQCK=total respiration of DOC+DOA used to represent microbial
C        activity (g C t-1)
C     DFNS=effect of microbial concentration on decomposition
C     OQCI=DOC product inhibition for decomposition
C     TFNX=temperature effect on biological activity
C     OSAH=total colonized SOC in each K (g C)
C     FCNK,FCPK=N,P limitation to microbial activity in each K
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      IF(OSAH(K).GT.ZEROS(NY,NX))THEN
      DO 775 M=1,2
      IF(ORC(M,K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNR=AMAX1(0.0,ORN(M,K,L,NY,NX)/ORC(M,K,L,NY,NX))
      CPR=AMAX1(0.0,ORP(M,K,L,NY,NX)/ORC(M,K,L,NY,NX))
      RDORC(M,K)=AMAX1(0.0,AMIN1(ORC(M,K,L,NY,NX)*XNFH
     2,SPORC(M)*ROQCK(K)*DFNS*OQCI*TFNX
     3*ORC(M,K,L,NY,NX)/OSAH(K)))
C    4*AMIN1(FCNK(K),FCPK(K)) 
      RDORN(M,K)=AMAX1(0.0,AMIN1(ORN(M,K,L,NY,NX)*XNFH
     2,CNR*RDORC(M,K)))
     3/FCNK(K)
      RDORP(M,K)=AMAX1(0.0,AMIN1(ORP(M,K,L,NY,NX)*XNFH
     2,CPR*RDORC(M,K)))
     3/FCPK(K) 
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
C     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate (g C,N,P,C)
C     CNH,CPH=N:C,P:C ratios of adsorbed C,N,P (g N,P g C-1)
C     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
C        (g C,N,P,C t-1)
C     SPOHC=specific decomposition rate constant for adsorbed C
C        (g C g C-1 h-1)
C     ROQCK=total respiration of DOC+DOA used to represent microbial
C        activity (g C t-1)
C     DFNS=effect of microbial concentration on decomposition
C     OQCI=DOC product inhibition for decomposition
C     TFNX=temperature effect on biological activity
C     OSAH=total colonized SOC in each K (g C)
C     FCNK,FCPK=N,P limitation to microbial activity in each K
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      IF(OSAH(K).GT.ZEROS(NY,NX))THEN
      IF(OHC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN 
      CNH(K)=AMAX1(0.0,OHN(K,L,NY,NX)/OHC(K,L,NY,NX))
      CPH(K)=AMAX1(0.0,OHP(K,L,NY,NX)/OHC(K,L,NY,NX))
      RDOHC(K)=AMAX1(0.0,AMIN1(OHC(K,L,NY,NX)*XNFH
     2,SPOHC*ROQCK(K)*DFNS*OQCI*TFNX
     3*OHC(K,L,NY,NX)/OSAH(K)))
C    4*AMIN1(FCNK(K),FCPK(K)) 
      RDOHN(K)=AMAX1(0.0,AMIN1(OHN(K,L,NY,NX)*XNFH
     2,CNH(K)*RDOHC(K)))
     3/FCNK(K)
      RDOHP(K)=AMAX1(0.0,AMIN1(OHP(K,L,NY,NX)*XNFH
     2,CPH(K)*RDOHC(K)))
     3/FCPK(K)
      RDOHA(K)=AMAX1(0.0,AMIN1(OHA(K,L,NY,NX)*XNFH
     2,SPOHA*ROQCK(K)*DFNS*TFNX
     3*OHA(K,L,NY,NX)/OSAH(K)))
C    4*AMIN1(FCNK(K),FCPK(K))
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
C     WRITE(*,4446)'RDOHC',I,J,NFZ,NX,NY,L,K,RDOHC(K),RDOHA(K)
C    2,OHC(K,L,NY,NX),XNFH,SPOHC,ROQCK(K),DFNS,OQCI,TFNX
C    3,OSAH(K),OHA(K,L,NY,NX),SPOHA
4446  FORMAT(A8,7I4,50E12.4)
C     ENDIF
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
C     VOLWM=soil water content (m3)
C     BKVL=soil mass (0= pond) (Mg)
C     FOSAH=fraction of TSRH in each substrate complex K
C     TSRH=total colonized SOC+microbial litter+adsorbed C (g C)
C     AEC,AECX=anion exchange capacity from soil file (mol Mg-1)
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores (g C,N,P,C)
C     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate (g C,N,P,C)
C     TSORP,HSORP=sorption rate constant (h-1) and coefficient 
C        for adsorption
C     FOCA,FOAA=fractions of DOC and acetate vs. DOC+acetate
C     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) of
C        OQC,acetate,DON,DOP (g C,C,N,P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      IF(VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX)
     2.AND.FOSAH(K,L,NY,NX).GT.ZERO)THEN
      IF(L.EQ.0)THEN
      AECX=0.5E+03
      ELSE
      AECX=AEC(L,NY,NX) 
      ENDIF 
      OQCX=AMAX1(ZEROS(NY,NX),OQC(K,L,NY,NX))
      OQNX=AMAX1(ZEROS(NY,NX),OQN(K,L,NY,NX))
      OQPX=AMAX1(ZEROS(NY,NX),OQP(K,L,NY,NX))
      OQAX=AMAX1(ZEROS(NY,NX),OQA(K,L,NY,NX))
      OHCX=AMAX1(ZEROS(NY,NX),OHC(K,L,NY,NX))
      OHNX=AMAX1(ZEROS(NY,NX),OHN(K,L,NY,NX))
      OHPX=AMAX1(ZEROS(NY,NX),OHP(K,L,NY,NX))
      OHAX=AMAX1(ZEROS(NY,NX),OHA(K,L,NY,NX))
      VOLXX=BKVL(L,NY,NX)*AECX*HSORP*FOSAH(K,L,NY,NX)
      VOLXW=VOLWM(NPH,L,NY,NX)*FOSAH(K,L,NY,NX)
      IF(FOCA(K).GT.ZERO)THEN
      VOLCX=FOCA(K)*VOLXX
      VOLCW=FOCA(K)*VOLXW
      CSORP(K)=TSORP*(OQCX*VOLCX-OHCX*VOLCW)/(VOLCX+VOLCW)*XNFH
      ELSE
      CSORP(K)=TSORP*(OQCX*VOLXX-OHCX*VOLXW)/(VOLXX+VOLXW)*XNFH
      ENDIF
      IF(FOAA(K).GT.ZERO)THEN
      VOLAX=FOAA(K)*VOLXX
      VOLAW=FOAA(K)*VOLXW
      CSORPA(K)=TSORP*(OQAX*VOLAX-OHAX*VOLAW)/(VOLAX+VOLAW)*XNFH
      ELSE
      CSORPA(K)=TSORP*(OQAX*VOLXX-OHAX*VOLXW)/(VOLXX+VOLXW)*XNFH
      ENDIF
      ZSORP(K)=TSORP*(OQNX*VOLXX-OHNX*VOLXW)/(VOLXX+VOLXW)*XNFH
      PSORP(K)=TSORP*(OQPX*VOLXX-OHPX*VOLXW)/(VOLXX+VOLXW)*XNFH
      ELSE
      CSORP(K)=0.0 
      CSORPA(K)=0.0
      ZSORP(K)=0.0
      PSORP(K)=0.0
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.0)THEN
C     WRITE(*,591)'CSORP',I,J,NFZ,NX,NY,L,K,CSORP(K),CSORPA(K) 
C    1,OQC(K,L,NY,NX),OHC(K,L,NY,NX),OQA(K,L,NY,NX),OHA(K,L,NY,NX)
C    2,OQC(K,L,NY,NX)/VOLWM(NPH,L,NY,NX)
C    2,OQA(K,L,NY,NX)/VOLWM(NPH,L,NY,NX)
C    3,OHC(K,L,NY,NX)/BKVL(L,NY,NX),OHA(K,L,NY,NX)/BKVL(L,NY,NX)
C    4,BKVL(L,NY,NX),VOLWM(NPH,L,NY,NX),FOCA(K),FOAA(K) 
C    5,FOSAH(K,L,NY,NX),OQCX
591   FORMAT(A8,7I4,40E12.4)
C     ENDIF
1790  CONTINUE
C
C     REDISTRIBUTE AUTOTROPHIC DECOMPOSITION PRODUCTS AMONG
C     HETEROTROPHIC SUBSTRATE-MICROBE COMPLEXES
C
C     FORC=fraction of total microbial residue in complex K
C     ORCT,TORC=population,total microbial residue (g C)
C     RCCMC,RCCMN,RCCMP=transfer of autotrophic litterfall C,N,P 
C        (K=5) to heterotrophic residue (g C,N,P t-1)
C     RCOMC,RCOMN,RCOMP=transfer of microbial litterfall C,N,P 
C        (K=5) to heterotrophic residue (g C,N,P t-1)
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
      DO 1680 M=1,2
      RCCMC(M,N,K)=(RCOMC(M,N,5)+RCMMC(M,N,5))*FORC(K)
      RCCMN(M,N,K)=(RCOMN(M,N,5)+RCMMN(M,N,5))*FORC(K)
      RCCMP(M,N,K)=(RCOMP(M,N,5)+RCMMP(M,N,5))*FORC(K)
C     IF(L.EQ.0)THEN
C     WRITE(*,8821)'RCCMC',I,J,L,K,N,M,RCCMC(M,N,K)
C    2,RCOMC(M,N,5),RCMMC(M,N,5),FORC(K) 
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
C     OSC,OSA,OSN,OSP=SOC,colonized SOC,SON,SOP (g C,C,N,P)
C     RDOSC,RDOSN,RDOSP=decomposition of SOC,SON,SOP (g C,N,P h-1)
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores (g C,C,N,P)
C     RCOSC,RCOSN,RCOSP=transfer of decomposition C,N,P to DOC,DON,DOP
C        (g C,N,P h-1)
C
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-RDOSC(M,K)
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-RDOSC(M,K)
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-RDOSN(M,K)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-RDOSP(M,K)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+RCOSC(M,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+RCOSN(M,K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+RCOSP(M,K)
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.L.EQ.1.AND.K.EQ.0)THEN
C     WRITE(*,4444)'RDOSC',I,J,NFZ,NX,NY,L,K,M,OSC(M,K,L,NY,NX)
C    2,RDOSC(M,K),OQN(K,L,NY,NX),RCOSN(M,K)
C     ENDIF
C
C     LIGNIFICATION PRODUCTS
C
C     RHOSC,RHOSN,RHOSP=transfer of decomposition C,N,P to POC,PON,POP
C        (g C,N,P h-1)
C     OSC,OSA,OSN,OSP=POC,colonized POC,PON,POP (g C,C,N,P) 
C
      IF(L.NE.0)THEN
      OSC(1,3,L,NY,NX)=OSC(1,3,L,NY,NX)+RHOSC(M,K)
      OSA(1,3,L,NY,NX)=OSA(1,3,L,NY,NX)+RHOSC(M,K)
      OSN(1,3,L,NY,NX)=OSN(1,3,L,NY,NX)+RHOSN(M,K)
      OSP(1,3,L,NY,NX)=OSP(1,3,L,NY,NX)+RHOSP(M,K)
      ELSE
      OSC(1,3,NU(NY,NX),NY,NX)=OSC(1,3,NU(NY,NX),NY,NX)+RHOSC(M,K)
      OSA(1,3,NU(NY,NX),NY,NX)=OSA(1,3,NU(NY,NX),NY,NX)+RHOSC(M,K)
      OSN(1,3,NU(NY,NX),NY,NX)=OSN(1,3,NU(NY,NX),NY,NX)+RHOSN(M,K)
      OSP(1,3,NU(NY,NX),NY,NX)=OSP(1,3,NU(NY,NX),NY,NX)+RHOSP(M,K)
      ENDIF
580   CONTINUE
C
C     MICROBIAL RESIDUE DECOMPOSITION PRODUCTS
C
C     ORC,ORN,ORP=microbial residue C,N,P (g C,N,P)
C     RDORC,RDORN,RDORP=decomposition of microbial residue C,N,P
C        (g C,N,P t-1)
C     RDOHC,RDOHN,RDOHP,RDOHA=decomposition of adsorbed C,N,P,acetate
C        (g C,N,P,C t-1)
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
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate (g C,N,P,C)
C     CGOQC,CGOAC,CGOMN,CGOMP=DOC,acetate,DON,DOP uptake from growth
C        respiration (g C,C,N,P t-1) 
C     RCH3X=acetate production from fermentation (g C t-1)
C
      DO 570 N=1,7
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-CGOQC(N,K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-CGOMN(N,K) 
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-CGOMP(N,K) 
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-CGOAC(N,K)+RCH3X(N,K)
C
C     MICROBIAL DECOMPOSITION PRODUCTS
C
C     ORC,ORN,ORP=microbial residue C,N,P (g C,N,P) 
C     RCOMC,RCOMN,RCOMP=transfer of microbial C,N,P litterfall to
C        microbial residue (g C,N,P t-1)
C     RCCMC,RCCMN,RCCMP=transfer of autotrophic litterfall C,N,P to
C        each heterotrophic population (g C,N,P t-1)
C     RCMMC,RCMMN,RCMMC=transfer of senesence litterfall C,N,P to
C        microbial residue (g C,N,P t-1)
C
      DO 565 M=1,2
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)+RCOMC(M,N,K)+RCCMC(M,N,K)
     2+RCMMC(M,N,K)
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)+RCOMN(M,N,K)+RCCMN(M,N,K)
     2+RCMMN(M,N,K)
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)+RCOMP(M,N,K)+RCCMP(M,N,K)
     2+RCMMP(M,N,K)
C     IF(K.EQ.2)THEN
C     WRITE(*,8821)'ORC',I,J,NFZ,NX,NY,L,K,M,N,ORC(M,K,L,NY,NX)
C    2,RCOMC(M,N,K),RCCMC(M,N,K),RCMMC(M,N,K),RDORC(M,K)
C     WRITE(*,8821)'ORP',I,J,L,K,N,M,ORP(M,K,L,NY,NX)
C    2,RCOMP(M,N,K),RCCMP(M,N,K),RCMMP(M,N,K),RDORP(M,K)
8821  FORMAT(A8,9I4,40E12.4)
C     ENDIF
565   CONTINUE
570   CONTINUE
C
C     SORPTION PRODUCTS
C
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores (g C,N,P,C)
C     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate (g C,N,P,C)
C     CSORP,CSORPA,ZSORP,PSORP=sorption(ad=+ve,de=-ve) 
C        of OQC,acetate,DON,DOP (g C,C,N,P t-1) 
C
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-CSORP(K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-ZSORP(K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-PSORP(K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-CSORPA(K)
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)+CSORP(K)
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)+ZSORP(K)
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)+PSORP(K)
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)+CSORPA(K)
C     IF(L.EQ.14)THEN
C     WRITE(*,592)'OQC',I,J,NFZ,NX,NY,L,K
C    2,OQC(K,L,NY,NX),(RCOSC(M,K),M=1,4),(RDORC(M,K),M=1,2)
C    3,RDOHC(K),(CGOQC(N,K),N=1,7),CSORP(K),OHC(K,L,NY,NX) 
C    4,OQA(K,L,NY,NX),RDOHA(K),(RCH3X(N,K),N=1,7)
C    3,(CGOAC(N,K),N=1,7),CSORPA(K),OHA(K,L,NY,NX) 
C    5,(WFN(N,K),N=1,7)
C     WRITE(*,592)'OQN',I,J,NFZ,NX,NY,L,K,OQN(K,L,NY,NX)
C    2,(RCOSN(M,K),M=1,4),(RDORN(M,K),M=1,2),RDOHN(K)
C    2,RCOQN,FORC(K),(CGOMN(N,K),N=1,7),ZSORP(K),OHN(K,L,NY,NX)  
592   FORMAT(A8,7I4,80E12.4)
C     ENDIF 
590   CONTINUE
C
C     MICROBIAL GROWTH FROM RESPIRATION, MINERALIZATION
C
C     OMC,OMN,OMP=microbial C,N,P (g C,N,P)
C     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural
C        C,N,P (g C,N,P t-1)
C     RXOMC,RXOMN,RXOMP=microbial C,N,P decomposition (g C,N,P t-1)
C     RXMMC,RXMMN,RXMMP=microbial C,N,P loss from senescence 
C        (g C,N,P t-1)
C
      DO 550 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 545 N=1,7
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      DO 540 M=1,2
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)+CGOMS(M,N,K) 
     2-RXOMC(M,N,K)-RXMMC(M,N,K) 
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)+CGONS(M,N,K) 
     2-RXOMN(M,N,K)-RXMMN(M,N,K)
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)+CGOPS(M,N,K) 
     2-RXOMP(M,N,K)-RXMMP(M,N,K) 
C     IF(L.EQ.3.AND.K.EQ.1.AND.N.EQ.6)THEN
C     WRITE(*,4488)'OMN2',I,J,NFZ,NX,NY,L,K,N,M,OMC(M,N,K,L,NY,NX)
C    2,OMN(M,N,K,L,NY,NX),OMP(M,N,K,L,NY,NX),CGONS(M,N,K) 
C    2,RXOMN(M,N,K),RXMMN(M,N,K)
C     WRITE(*,4488)'RDOMC',I,J,NFZ,NX,NY,L,K,N,M,CGOMS(M,N,K),CGOQC(N,K)
C    4,CGOAC(N,K),RGOMO(N,K),RGOMD(N,K),RXOMC(M,N,K),RXMMC(M,N,K)
C    3,RMOMC(M,N,K),TFNX,OMGR,OMC(3,N,K,L,NY,NX),WFN(N,K)
C    3,OMC(M,N,K,L,NY,NX),OMA(N,K),TSRH 
C    4,RCH3X(N,K),RH2GZ,RH2GX(4,K),FOCA(K),FOAA(K) 
C    6,OQA(K,L,NY,NX),OHA(K,L,NY,NX),OQC(K,L,NY,NX),OHC(K,L,NY,NX)
C    7,OMP(M,N,K,L,NY,NX),CGOPS(M,N,K),RDOMP(M,N,K),RDMMP(M,N,K)
C    8,OMP(3,N,K,L,NY,NX),CGOMP(N,K),RIPO4(N,K)
4488  FORMAT(A8,9I4,40E12.4)
C     ENDIF
C
C     HUMIFICATION PRODUCTS
C
C     OSC,OSA,OSN,OSP=SOC,colonized SOC,SON,SOP (g C,C,N,P)
C     CFOMC=fractions allocated to humic(1) vs fulvic(2) humus
C     RHOMC,RHOMN,RHOMP=transfer of microbial litterfall C,N,P 
C        to humus (g C,N,P t-1)
C     RHMMC,RHMMN,RHMMC=transfer of senesence litterfall C,N,P 
C        to humus (g C,N,P t-1)
C     NU=surface layer number for litter humification products
C
      IF(L.NE.0)THEN
      OSC(1,4,L,NY,NX)=OSC(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
     2*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSA(1,4,L,NY,NX)=OSA(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
     2*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSN(1,4,L,NY,NX)=OSN(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
     2*(RHOMN(M,N,K)+RHMMN(M,N,K))
      OSP(1,4,L,NY,NX)=OSP(1,4,L,NY,NX)+CFOMC(1,L,NY,NX)
     2*(RHOMP(M,N,K)+RHMMP(M,N,K))
      OSC(2,4,L,NY,NX)=OSC(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
     2*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSA(2,4,L,NY,NX)=OSA(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
     2*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSN(2,4,L,NY,NX)=OSN(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
     2*(RHOMN(M,N,K)+RHMMN(M,N,K))
      OSP(2,4,L,NY,NX)=OSP(2,4,L,NY,NX)+CFOMC(2,L,NY,NX)
     2*(RHOMP(M,N,K)+RHMMP(M,N,K))
C     IF((I/10)*10.EQ.I.AND.J.EQ.24)THEN
C     WRITE(*,4445)'RHOMC',I,J,NFZ,NX,NY,L,K,M,N,OSC(1,4,L,NY,NX)
C    2,OSC(2,4,L,NY,NX),CFOMC(1,L,NY,NX),CFOMC(2,L,NY,NX)
C    3,RHOMC(M,N,K),RHMMC(M,N,K)
4445  FORMAT(A8,9I4,40E12.4)
C     ENDIF 
      ELSE
      OSC(1,4,NU(NY,NX),NY,NX)=OSC(1,4,NU(NY,NX),NY,NX)
     2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSA(1,4,NU(NY,NX),NY,NX)=OSA(1,4,NU(NY,NX),NY,NX)
     2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSN(1,4,NU(NY,NX),NY,NX)=OSN(1,4,NU(NY,NX),NY,NX)
     2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMN(M,N,K)+RHMMN(M,N,K))
      OSP(1,4,NU(NY,NX),NY,NX)=OSP(1,4,NU(NY,NX),NY,NX)
     2+CFOMC(1,NU(NY,NX),NY,NX)*(RHOMP(M,N,K)+RHMMP(M,N,K))
      OSC(2,4,NU(NY,NX),NY,NX)=OSC(2,4,NU(NY,NX),NY,NX)
     2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSA(2,4,NU(NY,NX),NY,NX)=OSA(2,4,NU(NY,NX),NY,NX)
     2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMC(M,N,K)+RHMMC(M,N,K))
      OSN(2,4,NU(NY,NX),NY,NX)=OSN(2,4,NU(NY,NX),NY,NX)
     2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMN(M,N,K)+RHMMN(M,N,K))
      OSP(2,4,NU(NY,NX),NY,NX)=OSP(2,4,NU(NY,NX),NY,NX)
     2+CFOMC(2,NU(NY,NX),NY,NX)*(RHOMP(M,N,K)+RHMMP(M,N,K))
      ENDIF
540   CONTINUE
C
C     INPUTS TO NONSTRUCTURAL POOLS
C
C     CGOMC=total DOC+acetate uptake (g C t-1)
C     RGOMO=total aerobic respiration (g C t-1)
C     RGOMD=respiration for denitrifcation (g C t-1)
C     RGN2F=respiration for N2 fixation (g C t-1)
C     OMC,OMN,OMP=nonstructural C,N,P (g C,N,P)
C     RCO2X=total CO2 emission (g C t-1)
C     CGOMS,CGONS,CGOPS=transfer from nonstructural to structural
C        C,N,P (g C,N,P t-1)
C     R3OMC,R3OMN,R3OMP=microbial C,N,P recycling (g C,N,P t-1)
C     R3MMC,R3MMN,R3MMP=microbial C,N,P recycling from senescence 
C        (g C,N,P t-1)
C     CGOMN,CGOMP=DON, DOP uptake (g N,P t-1)
C     RINH4,RINB4=NH4 mineralization-immobilization in non-band, band
C        (g N t-1)
C     RINO3,RINB3=NO3 immobilization in non-band, band (g N t-1)
C     RIPO4,RIPOB=H2PO4 mineralization-immobilization in non-band, band 
C        (g P t-1)
C     RIP14,RIP1B=HPO4 mineralization-immobilization in non-band, band 
C        (g P t-1)
C     RINH4R,RINO3R=NH4,NO3 mineralization-immobilization in surface
C         litter (g P t-1)
C     RIPO4R,RIP14R=H2PO4,HPO4 mineralization-immobilization in surface
C         litter (g P t-1)
C     RN2FX=N2 fixation rate (g N t-1)
C
      CGROMC=CGOMC(N,K)-RGOMO(N,K)-RGOMD(N,K)-RGN2F(N,K)
      RCO2X(N,K)=RCO2X(N,K)+RGN2F(N,K) 
      DO 555 M=1,2
      OMC(3,N,K,L,NY,NX)=OMC(3,N,K,L,NY,NX)-CGOMS(M,N,K)
     2+R3OMC(M,N,K) 
      OMN(3,N,K,L,NY,NX)=OMN(3,N,K,L,NY,NX)-CGONS(M,N,K)
     2+R3OMN(M,N,K)+R3MMN(M,N,K)
      OMP(3,N,K,L,NY,NX)=OMP(3,N,K,L,NY,NX)-CGOPS(M,N,K)
     2+R3OMP(M,N,K)+R3MMP(M,N,K)
      RCO2X(N,K)=RCO2X(N,K)+R3MMC(M,N,K) 
555   CONTINUE
      OMC(3,N,K,L,NY,NX)=OMC(3,N,K,L,NY,NX)+CGROMC
      OMN(3,N,K,L,NY,NX)=OMN(3,N,K,L,NY,NX)+CGOMN(N,K) 
     2+RINH4(N,K)+RINB4(N,K)+RINO3(N,K)+RINB3(N,K)+RN2FX(N,K) 
      OMP(3,N,K,L,NY,NX)=OMP(3,N,K,L,NY,NX)+CGOMP(N,K) 
     2+RIPO4(N,K)+RIPOB(N,K)+RIP14(N,K)+RIP1B(N,K) 
      IF(L.EQ.0)THEN
      OMN(3,N,K,L,NY,NX)=OMN(3,N,K,L,NY,NX)+RINH4R(N,K)+RINO3R(N,K)
      OMP(3,N,K,L,NY,NX)=OMP(3,N,K,L,NY,NX)+RIPO4R(N,K)+RIP14R(N,K)
      ENDIF
C     IF(L.EQ.1.AND.K.EQ.1)THEN
C     WRITE(*,5556)'OMC3',I,J,NFZ,NX,NY,L,K,N,OMC(3,N,K,L,NY,NX)
C    2,OMN(3,N,K,L,NY,NX),OMP(3,N,K,L,NY,NX),OMC(1,N,K,L,NY,NX)
C    2,OMN(1,N,K,L,NY,NX),OMP(1,N,K,L,NY,NX),WFN(N,K),OXYI
C    2,COXYS(L,NY,NX) 
C    2,CGOMS(1,N,K),CGOMS(2,N,K),CGROMC
C    3,CGOPS(1,N,K),CGOPS(2,N,K),CGOMP(N,K),RIPO4(N,K)
C    4,CGOMC(N,K),RGOMO(N,K),RGOMD(N,K),RMOMT,WFN(N,K)
C    5,(CGONS(M,N,K),M=1,2),(R3OMN(M,N,K),M=1,2),(R3MMN(M,N,K)
C    6,M=1,2),(XOMCZ(M,N,K),M=1,2)
C    6,CGOMN(N,K),RINH4(N,K),RINB4(N,K),RINO3(N,K),RINB3(N,K)
C    7,RN2FX(N,K),RINH4R(N,K),RINO3R(N,K) 
5556  FORMAT(A8,8I4,60E12.4)
C     ENDIF
      ENDIF
545   CONTINUE
      ENDIF
550   CONTINUE
C
C     MICROBIAL COLONIZATION OF NEW LITTER
C
C     OSCT,OSAT,OSCX=total,colonized,uncolonized SOC (g C)
C     OSA,OSC=colonized,total litter (g C)
C     DOSA=rate constant for litter colonization (h-1) 
C     ROQCK=total respiration of DOC+DOA used to represent microbial
C        activity (g C t-1)
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
C     IF((I/10)*10.EQ.I.AND.J.EQ.15.AND.NFZ.EQ.1)THEN
C     WRITE(*,8822)'OSAT',I,J,NFZ,NX,NY,L,K,M,OSA(M,K,L,NY,NX)
C    2,OSC(M,K,L,NY,NX),DOSA(K),ROQCK(K),DOSAK,OSAT(K),OSCT(K)
8822  FORMAT(A8,8I4,30E12.4)
C     ENDIF
485   CONTINUE
      ELSE
      DO 490 M=1,4
      OSA(M,K,L,NY,NX)=AMIN1(OSC(M,K,L,NY,NX),OSA(M,K,L,NY,NX)) 
490   CONTINUE
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.L.EQ.1.AND.K.EQ.0)THEN
C     WRITE(*,8823)'OSC',I,J,NFZ,NX,NY,L,K
C    2,((OMC(M,N,K,L,NY,NX),N=1,7),M=1,3)
C    2,(ORC(M,K,L,NY,NX),M=1,2),OQC(K,L,NY,NX),OQCH(K,L,NY,NX)
C    3,OHC(K,L,NY,NX),OQA(K,L,NY,NX),OQAH(K,L,NY,NX),OHA(K,L,NY,NX)
C    4,(OSC(M,K,L,NY,NX),M=1,4)
8823  FORMAT(A8,7I4,100E14.6)
C     ENDIF
480   CONTINUE
C
C     AGGREGATE ALL TRANSFORMATIONS CALCULATED ABOVE FOR EACH 
C     POPULATION N IN EACH COMPLEX K
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
      IF(K.NE.5.OR.(N.LE.3.OR.N.EQ.5))THEN
      TRINH=TRINH+RINH4(N,K)
      TRINO=TRINO+RINO3(N,K)
      TRIPO=TRIPO+RIPO4(N,K)
      TRIP1=TRIP1+RIP14(N,K)
      TRINB=TRINB+RINB4(N,K)
      TRIOB=TRIOB+RINB3(N,K)
      TRIPB=TRIPB+RIPOB(N,K)
      TRIB1=TRIB1+RIP1B(N,K)
      TRN2F=TRN2F+RN2FX(N,K)
      IF(L.EQ.NU(NY,NX))THEN
      TRINH=TRINH+RINH4R(N,K)
      TRINO=TRINO+RINO3R(N,K)
      TRIPO=TRIPO+RIPO4R(N,K)
      TRIP1=TRIP1+RIP14R(N,K)
      ENDIF
C     IF(I.EQ.322)THEN
C     WRITE(*,4469)'TRINH',I,J,NFZ,NX,NY,L,K,N
C    2,TRINH,RINH4(N,K),RINH4R(N,K)
C     WRITE(*,4469)'TRIPO',I,J,NFZ,NX,NY,L,K,N
C    2,TRIPO,RIPO4(N,K),RIPO4R(N,K)
C    2,CGOMP(N,K)
4469  FORMAT(A8,8I4,20E12.4)
C     ENDIF
      TRGOM=TRGOM+RCO2X(N,K)
      TRGOC=TRGOC+RCH4X(N,K)
      TRGOD=TRGOD+RGOMD(N,K)
      TUPOX=TUPOX+RUPOX(N,K)
      TRDN3=TRDN3+RDNO3(N,K)
      TRDNB=TRDNB+RDNOB(N,K)
      TRDN2=TRDN2+RDNO2(N,K)
      TRD2B=TRD2B+RDN2B(N,K)
      TRDNO=TRDNO+RDN2O(N,K)
      TRGOH=TRGOH+RH2GX(N,K)
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.1)THEN
C     WRITE(*,3333)'TRGOM',I,J,NFZ,NX,NY,L,K,N,TRGOM
C    2,RCO2X(N,K),TRGOA,RGOMO(N,K),WFN(N,K),TRGOC,RCH4X(N,K) 
C     WRITE(*,3333)'TUPOX',I,J,NFZ,NX,NY,L,K,N,TUPOX,RUPOX(N,K)
C    2,ROXYM(N,K),WFN(N,K) 
C     WRITE(*,3333)'N2O',I,J,NFZ,NX,NY,L,K,N,TRDN2,TRD2B,TRDNO
C    2,RDNO2(N,K),RDN2B(N,K),RDN2O(N,K),COXYS(L,NY,NX)
C    3,COXYG(L,NY,NX),ROXYM(N,K)-ROXYO(N,K)
C     WRITE(*,3333)'TRGOH',I,J,NFZ,NX,NY,L,K,N,TRGOH,RH2GX(N,K)
C    2,RGOMO(N,K)
3333  FORMAT(A8,8I4,20E12.4)
C     ENDIF
      ENDIF
640   CONTINUE
      ENDIF
650   CONTINUE
      DO 645 N=1,7
      IF(N.LE.3.OR.N.EQ.5)THEN
      IF(N.NE.3)THEN
      TRGOA=TRGOA+CGOMC(N,5) 
      ENDIF
      ENDIF
645   CONTINUE
C
C     AGGREGATE TRANSFORMATIONS INTO ARRAYS FOR TRANSFER TO 'REDIST.F'
C
C     RCO2O=net CO2 uptake (g C t-1)
C     TRGOA=total CO2 uptake by autotrophs (g C t-1) 
C     TRGOM total CO2 emission by heterotrophs reducing O2 (g C t-1)
C     TRGOD=total CO2 emission by denitrifiers reducing NOx (g C t-1)
C     RVOXA(3)=CH4 oxidation (g C t-1)
C     RCH4O=net CH4 uptake (g C t-1)
C     CGOMC=total CH4 uptake by autotrophs (g C t-1)
C     TRGOC=total CH4 emission (g C t-1)
C     RH2GO=net H2 uptake (g H t-1)
C     RH2GZ,TRGOH=total H2 uptake, emission (g H t-1)
C     RUPOXO,TUPOX=total O2 uptake (g O t-1)
C     RN2G=total N2 production (g N t-1)
C     TRDNO=total N2O reduction (g N t-1)
C     RN2O=total N2O uptake (g N t-1)
C     TRDN2,TRD2B=total NO2 reduction in non-band,band (g N t-1)
C     RCN2O,RCN2B=nitrous acid reduction in non-band,band (g N t-1)
C
      RCO2O(L,NY,NX)=TRGOA-TRGOM-TRGOD-RVOXA(3)
      RCH4O(L,NY,NX)=RVOXA(3)+CGOMC(3,5)-TRGOC
      RH2GO(L,NY,NX)=RH2GZ-TRGOH
      RUPOXO(L,NY,NX)=TUPOX
      RN2G(L,NY,NX)=-TRDNO
      RN2O(L,NY,NX)=-TRDN2-TRD2B-RCN2O-RCN2B+TRDNO
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.1)THEN
C     WRITE(*,2468)'RCO2O',I,J,NFZ,NX,NY,L,RCO2O(L,NY,NX)
C    2,TRGOA,TRGOM,TRGOD,RVOXA(3),RCH4O(L,NY,NX)
C    3,CGOMC(3,5),TRGOC
C     WRITE(*,2468)'RCH4O',I,J,NFZ,NX,NY,L,RCH4O(L,NY,NX)
C    2,RVOXA(3),CGOMC(3,5),TRGOC,WFN(3,5)
C     WRITE(*,2468)'RN2O',I,J,NFZ,NX,NY,L
C    2,RN2O(L,NY,NX),TRDN2,TRD2B,RCN2O,RCN2B,TRDNO
C    2,RCH4O(L,NY,NX),RVOXA(3)
C    2,CGOMC(3,5),TRGOC,(OMA(N,1),N=1,7)
C     WRITE(*,2468)'RN2G',I,J,NFZ,NX,NY,L
C    2,RN2G(L,NY,NX),TRDNO,TRDN2,TRD2B,RCN2O,RCN2B 
2468  FORMAT(A8,6I4,20E12.4)
c     ENDIF
C
C     AGGREGATE FLUXES FOR TRANSFER TO ‘TRNSFR.F’ AND ‘REDIST.F’
C
C     XOQCS,XOQNZ,XOQPS,XOQAS=net change in DOC,DON,DOP,acetate 
C        (g C,N,P,C t-1)
C     XNH4S,XNH4B=net change in NH4 in band,non-band (g N t-1)
C     TRINH,TRINB=total NH4 mineralization-immobilization 
C        in non-band,band (g N t-1)
C     RVOXA(1),RVOXB(1)=total NH4 oxidation in non-band,band (g O t-1) 
C     XNO3S,XNO3B=net change in NO3 in band,non-band (g N t-1)
C     TRINO,TRIOB=total NO3 immobilization in non-band,band (g N t-1)
C     RVOXA(2),RVOXB(2)=total NO2 oxidation in non-band,band (g O t-1)
C     TRDN3,TRDNB=total NO3 reduction in non-band,band (g N t-1)
C     RCNO3,RCN3B=NO3 production from nitrous acid reduction 
C        in non-band,band (g N t-1)
C     XNO2S,XNO2B=net change in NO3 in band,non-band (g N t-1) 
C     TRDN2,TRD2B=total NO2 reduction in non-band,band (g N t-1)
C     RCNO2,RCNOB=nitrous acid reduction in non-band,band (g N t-1)
C     XH2PS,XH2BS=net change in H2PO4 in band,non-band (g P t-1) 
C     TRIPO,TRIPB=total H2PO4 mineralization-immobilization 
C        in non-band,band (g P t-1)
C     XH1PS,XH1BS=net change in HPO4 in band,non-band (g P t-1)
C     TRIP1,TRIB1=total HPO4 mineralization-immobilization   
C        in non-band,band (g P t-1)
C     XN2GS=total N2 fixation (g N t-1)
C     XZHYS=net H+ production (g H t-1)
C
      DO 655 K=0,KL
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
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-CGOQC(N,K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-CGOMN(N,K) 
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-CGOMP(N,K) 
      XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)-CGOAC(N,K)+RCH3X(N,K)
670   CONTINUE
      XOQCS(K,L,NY,NX)=XOQCS(K,L,NY,NX)-CSORP(K)
      XOQNS(K,L,NY,NX)=XOQNS(K,L,NY,NX)-ZSORP(K)
      XOQPS(K,L,NY,NX)=XOQPS(K,L,NY,NX)-PSORP(K)
      XOQAS(K,L,NY,NX)=XOQAS(K,L,NY,NX)-CSORPA(K)
655   CONTINUE
      XNH4S(L,NY,NX)=-TRINH-RVOXA(1) 
      XNO3S(L,NY,NX)=-TRINO+RVOXA(2)-TRDN3+RCNO3
      XNO2S(L,NY,NX)=RVOXA(1)-RVOXA(2)+TRDN3-TRDN2-RCNO2
      XH2PS(L,NY,NX)=-TRIPO
      XH1PS(L,NY,NX)=-TRIP1
      XNH4B(L,NY,NX)=-TRINB-RVOXB(1) 
      XNO3B(L,NY,NX)=-TRIOB+RVOXB(2)-TRDNB+RCN3B
      XNO2B(L,NY,NX)=RVOXB(1)-RVOXB(2)+TRDNB-TRD2B-RCNOB
      XH2BS(L,NY,NX)=-TRIPB
      XH1BS(L,NY,NX)=-TRIB1
      XN2GS(L,NY,NX)=TRN2F
      TFNQ(L,NY,NX)=TFNX
      VOLQ(L,NY,NX)=VOLWY
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C     XZHYS=net H+ production-consumption from oxidation-reduction 
C        reactions driving pH changes in ‘solute.f’
C
      IF(ISALTG.NE.0)THEN 
      XZHYS(L,NY,NX)=XZHYS(L,NY,NX)
     2+AMAX1(-ZHY(L,NY,NX)*XNFH,0.1429*(RVOXA(1)+RVOXB(1) 
     2-TRDN3-TRDNB)-0.0714*(TRDN2+TRD2B+TRDNO))
C     WRITE(*,2323)'XZHYSN',I,J,NFZ,NX,NY,L,XZHYS(L,NY,NX)
C    2,-ZHY(L,NY,NX),RVOXA(1),RVOXB(1) 
C    2,-TRDN3,-TRDNB,-(TRDN2+TRD2B+TRDNO)
2323  FORMAT(A8,6I4,20F16.8)
      ENDIF
C     IF(I.EQ.322)THEN
C     WRITE(*,2324)'XNH4S',I,J,NFZ,NX,NY,L,XNH4S(L,NY,NX)
C    2,TRINH,RVOXA(1),VLNH4(L,NY,NX),TRDN2
C     WRITE(*,2324)'XNO3S',I,J,NFZ,NX,NY,L,NPH,XNO3S(L,NY,NX)
C    2,TRINO,RVOXA(2),VLNO3(L,NY,NX),TRDN3,RCNO3
C     WRITE(*,2324)'XNO2S',I,J,NFZ,NX,NY,L,NPH,XNO2S(L,NY,NX)
C    2,RVOXA(1),RVOXA(2),TRDN3,TRDN2,RCNO2
C    3,ZNO2S(L,NY,NX),FNO2,VMXC4S
C    4,CHNO2,VOLWM(NPH,L,NY,NX),FNO3S,TFNX,XNFH 
C     WRITE(*,2324)'XH2PS',I,J,NFZ,NX,NY,L,XH2PS(L,NY,NX)
C    2,TRIPO,VLPO4(L,NY,NX)
C     WRITE(*,2324)'XNO2B',I,J,NFZ,NX,NY,L,XNO2B(L,NY,NX),RVOXB(1)
C    2,VLNHB(L,NY,NX),RVOXB(2),VLNOB(L,NY,NX),TRDNB,TRD2B,RCNOB
C     WRITE(*,2324)'XOQCS',I,J,NFZ,NX,NY,L,(XOQCS(K,L,NY,NX),K=0,KL)
2324  FORMAT(A8,6I4,30E12.4)
C     ENDIF
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
C     MIX MICROBIAL C,N,P BETWEEN ADJACENT SOIL LAYERS L AND LL
C
C     IDENTIFY NEXT LOWER SOIL LAYER
C
C     BKVL=soil mass (Mg)
C
      IF(L.EQ.0.OR.(L.GE.NU(NY,NX).AND.L.LT.NL(NY,NX)))THEN
      IF(BKDS(L,NY,NX).GT.ZERO.AND.BKDS(L+1,NY,NX).GT.ZERO)THEN
      IF(L.EQ.0)THEN
      LL=NU(NY,NX)
      ELSE
C     DO 1100 LN=NU(NY,NX)+1,NL(NY,NX)
      DO 1100 LN=L+1,NL(NY,NX)
      IF(VOLX(LN,NY,NX).GT.ZEROS2(NY,NX))THEN
      LL=LN
      GO TO 1101
      ENDIF
1100  CONTINUE
1101  CONTINUE
      ENDIF
C
C     MIX MICROBIAL BIOMASS BETWEEN ADJACENT SOIL LAYERS L AND LL
C
C     TOQCK=total active biomass respiration activity (g C t-1)
C     TOQCKL,TOQCKLL=total active biomass respiration activity 
C        per unit soil volume (g C m-3 t-1)
C     CORGCI=total SOC concentration from soil file (g Mg-1)
C     DLYR=thickness of soil layer (m)
C     OMC,OMN,OMP=microbial C,N,P (g C,N,P)
C     VOLX=soil layer volume (m3)
C     FOMCX=mixing fraction for OMC,OMN,OMP
C     FPRIMB=rate constant for mixing soil microbial biomass 
C        among soil layers (h-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      TOQCKL=TOQCK(L,NY,NX)/VOLT(L,NY,NX)
      TOQCKLL=TOQCK(LL,NY,NX)/VOLT(LL,NY,NX)
      CORGCLL=AMIN1(CORGCI(L,NY,NX),CORGCI(LL,NY,NX))*1.82E-06
      IF(BKDS(LL,NY,NX).GT.ZERO)THEN 
      FOMCX=(2.0*FPRIMB*(TOQCKL-TOQCKLL)*XNFH
     2/(DLYR(3,L,NY,NX)+DLYR(3,LL,NY,NX)))*CORGCLL
      ELSE
      FOMCX=0.0
      ENDIF
      DO 7960 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 7961 N=1,7
      DO 7962 M=1,3
      IF(ABS(FOMCX).GT.ZERO)THEN
      IF(FOMCX.GT.0.0)THEN 
      OMCXS=FOMCX*AMAX1(0.0,OMC(M,N,K,L,NY,NX))
      OMNXS=FOMCX*AMAX1(0.0,OMN(M,N,K,L,NY,NX))
      OMPXS=FOMCX*AMAX1(0.0,OMP(M,N,K,L,NY,NX))
      ELSE
      OMCXS=FOMCX*AMAX1(0.0,OMC(M,N,K,LL,NY,NX))
      OMNXS=FOMCX*AMAX1(0.0,OMN(M,N,K,LL,NY,NX))
      OMPXS=FOMCX*AMAX1(0.0,OMP(M,N,K,LL,NY,NX))
      ENDIF
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)-OMCXS
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)-OMNXS
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)-OMPXS
      OMC(M,N,K,LL,NY,NX)=OMC(M,N,K,LL,NY,NX)+OMCXS
      OMN(M,N,K,LL,NY,NX)=OMN(M,N,K,LL,NY,NX)+OMNXS
      OMP(M,N,K,LL,NY,NX)=OMP(M,N,K,LL,NY,NX)+OMPXS
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.NFZ.EQ.1.AND.L.EQ.0)THEN
C     WRITE(*,5561)'OMX',I,J,NFZ,NX,NY,L,LL,K,N,M,OMC(M,N,K,L,NY,NX)
C    2,OMN(M,N,K,L,NY,NX),OMP(M,N,K,L,NY,NX),OMC(M,N,K,LL,NY,NX)
C    2,OMN(M,N,K,LL,NY,NX),OMP(M,N,K,LL,NY,NX)
C    3,OMCXS,OMNXS,OMPXS,FOMCX,OMCXD 
5561  FORMAT(A8,10I4,20E12.4)
C     ENDIF 
      ENDIF
7962  CONTINUE
7961  CONTINUE
      ENDIF
7960  CONTINUE
      ENDIF
      ENDIF
C     IF((I/1)*1.EQ.I.AND.J.EQ.19.AND.L.LE.5)THEN
C     WRITE(*,2123)'TOTALL',I,J,NX,NY,L,TFOXYX,TFNH4X
C    2,TFNO3X,TFPO4X,TFNH4B,TFNO3B,TFPO4B,TFNO2X,TFNO2B
C    3,TFOQC,TFOQA 
2123  FORMAT(A8,5I4,12E15.4)
C     ENDIF
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
C     IF(NX.EQ.1.AND.J.EQ.24)THEN
C     DO 4344 K=0,5
C     WRITE(*,4343)'ORGNIT',I,J,NFZ,NX,NY,L,K
C    2,ORGC(L,NY,NX),ORGR(L,NY,NX),DC,OC 
C    2,((OMC(M,N,K,L,NY,NX),M=1,3),N=1,7)
C    3,(ORC(M,K,L,NY,NX),M=1,2),(OSC(M,K,L,NY,NX),M=1,4)
C    4,OQC(K,L,NY,NX),OQCH(K,L,NY,NX),OHC(K,L,NY,NX)
C    2,OQA(K,L,NY,NX),OQAH(K,L,NY,NX),OHA(K,L,NY,NX)
C    5,(OSC(M,K,L,NY,NX),M=1,4)
C     WRITE(*,4343)'ORGN',I,J,NX,NY,L,K,ORGN(L,NY,NX),DN,ON
C    2,((OMN(M,N,K,L,NY,NX),M=1,3),N=1,7)
C    3,(ORN(M,K,L,NY,NX),M=1,2),(OSN(M,K,L,NY,NX),M=1,4)
C    4,OQN(K,L,NY,NX),OQNH(K,L,NY,NX),OHN(K,L,NY,NX)
4343  FORMAT(A8,7I4,120E12.4)
4344  CONTINUE
C     ENDIF
998   CONTINUE
C     WRITE(20,3434)'RN2O',IYRC,I,J,(RN2O(L,NY,NX),L=0,NL(NY,NX))
3434  FORMAT(A8,3I4,20E12.4)
C
C     IF FIRE EVENT IS IN PROGRESS
C
C     ICHKF=fire flag :0=no fire,1=fire
C     TKS=soil temperature (K)
C     TFNCO=temperature function for combusting SOC (600K = 1)
C     TFNCOS=TFNCO used in ‘trnsfr.f’ and ‘redist.f’
C     SPCMB=specific combustion rate of SOC(K=0,4) at 600K,
C        charcoal(K=5) at 700K (g C m-2 h-1)
C     ORGC=total SOC (g C)
C     FRCBCO=combustion fraction of total SOC (g C g C-1 t-1)
C     AREA=grid cell area (m2)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      IF(ICHKF.EQ.1)THEN
      DO 2850 L=0,NL(NY,NX)
      RTK=8.3143*TKS(L,NY,NX)
      TFNCO=AMIN1(TFNCX,EXP(12.028-60000/RTK))
      TFNCOS(L,NY,NX)=AMIN1(1.0,TFNCO)
      DO 2865 K=0,5
      IF(ORGC(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FRCBCO(K)=AMIN1(1.0,SPCMB(K)*TFNCO/ORGC(L,NY,NX)
     2*AREA(3,L,NY,NX)*XNFH)
      ELSE
      FRCBCO(K)=0.0
      ENDIF
2865  CONTINUE
C
C     COMBUST MICROBIAL BIOMASS
C
C     FRCBCO=combustion fraction of total SOC
C     SPCMB=specific combustion rate of SOC(K=0,4) at 600K,
C        charcoal(K=5) at 700K (g C g C-1 h-1)
C     TFNCO,TFNCH=temperature function for combusting SOC (600K = 1)
C        charcoal (700K = 1)
C     VOLX=soil layer volume (m3)
C     OMC,OMN,OMP=microbial C,N,P biomass (g C,N,P)
C     RCBOMC,RCBOMN,RCBOMP=combustion rate of OMC,OMN,OMP (g C,N,P t-1)
C     RCGSK=total soil combustion for use in ‘trnsfr.f’ and 
C        ‘redist.f’(g C t-1)
C     VCOXFS,VNOXFS,VPOXFS=total C,N,P combustion loss for output 
C        (g C,N,P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      DO 2870 K=0,5
      IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
      DO 2860 N=1,7
      DO 2860 M=1,3
      IF(OMC(M,N,K,L,NY,NX).GT.ZEROS(NY,NX)
     2.AND.TKS(L,NY,NX).GT.TCMBX)THEN
      RCBOMC(M,N,K,L,NY,NX)=FRCBCO(K)*OMC(M,N,K,L,NY,NX)
      RCBOMN(M,N,K,L,NY,NX)=RCBOMC(M,N,K,L,NY,NX)
     2*OMN(M,N,K,L,NY,NX)/OMC(M,N,K,L,NY,NX)
      RCBOMP(M,N,K,L,NY,NX)=RCBOMC(M,N,K,L,NY,NX)
     2*OMP(M,N,K,L,NY,NX)/OMC(M,N,K,L,NY,NX)
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)-RCBOMC(M,N,K,L,NY,NX) 
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)-RCBOMN(M,N,K,L,NY,NX)
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)-RCBOMP(M,N,K,L,NY,NX)
      RCGSK(L,NY,NX)=RCGSK(L,NY,NX)+RCBOMC(M,N,K,L,NY,NX)
      VCOXFS(NY,NX)=VCOXFS(NY,NX)-RCBOMC(M,N,K,L,NY,NX)
      VNOXFS(NY,NX)=VNOXFS(NY,NX)-RCBOMN(M,N,K,L,NY,NX)
      VPOXFS(NY,NX)=VPOXFS(NY,NX)-RCBOMP(M,N,K,L,NY,NX)
      ELSE
      RCBOMC(M,N,K,L,NY,NX)=0.0
      RCBOMN(M,N,K,L,NY,NX)=0.0
      RCBOMP(M,N,K,L,NY,NX)=0.0
      ENDIF
C     WRITE(*,4448)'RCBOMN',I,J,NFZ,NX,NY,L,K,M,N
C    2,RCBOMN(M,N,K,L,NY,NX),OMN(M,N,K,L,NY,NX),FRCBCO(K)
C    4,RCGSK(L,NY,NX),VCOXFS(NY,NX),VNOXFS(NY,NX),VPOXFS(NY,NX)
C    5,SPCMB(K),TFNCO,XNFH,ORGC(L,NY,NX)
4448  FORMAT(A8,9I4,30E12.4)
2860  CONTINUE
      ENDIF
2870  CONTINUE
C
C     COMBUST MICROBIAL RESIDUE
C
C     FRCBCO=combustion fraction of total SOC
C     SPCMB=specific combustion rate of SOC(K=0,4) at 600K,
C        charcoal(K=5) at 700K (g C g C-1 h-1)
C     TFNCO,TFNCH=temperature function for combusting SOC (600K = 1)
C        charcoal (700K = 1)
C     VOLX=soil layer volume (m3)
C     ORC,ORN,ORP=microbial residue C,N,P (g C,N,P)
C     RCBORC,RCBORN,RCBORP=combustion rate of ORC,ORN,ORP (g C,N,P t-1)
C     RCGSK=total combustion for use in ‘trnsfr.f’ (g C t-1)
C     VCOXFS,VNOXFS,VPOXFS=total C,N,P combustion loss for output 
C        (g C,N,P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      DO 2800 K=0,4
      DO 2840 M=1,2
      IF(ORC(M,K,L,NY,NX).GT.ZEROS(NY,NX) 
     2.AND.TKS(L,NY,NX).GT.TCMBX)THEN
      RCBORC(M,K,L,NY,NX)=FRCBCO(K)*ORC(M,K,L,NY,NX)
      RCBORN(M,K,L,NY,NX)=RCBORC(M,K,L,NY,NX)
     2*ORN(M,K,L,NY,NX)/ORC(M,K,L,NY,NX)
      RCBORP(M,K,L,NY,NX)=RCBORC(M,K,L,NY,NX)
     2*ORP(M,K,L,NY,NX)/ORC(M,K,L,NY,NX)
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-RCBORC(M,K,L,NY,NX)
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-RCBORN(M,K,L,NY,NX)
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-RCBORP(M,K,L,NY,NX)
      RCGSK(L,NY,NX)=RCGSK(L,NY,NX)+RCBORC(M,K,L,NY,NX)
      VCOXFS(NY,NX)=VCOXFS(NY,NX)-RCBORC(M,K,L,NY,NX)
      VNOXFS(NY,NX)=VNOXFS(NY,NX)-RCBORN(M,K,L,NY,NX)
      VPOXFS(NY,NX)=VPOXFS(NY,NX)-RCBORP(M,K,L,NY,NX)
      ELSE
      RCBORC(M,K,L,NY,NX)=0.0
      RCBORN(M,K,L,NY,NX)=0.0
      RCBORP(M,K,L,NY,NX)=0.0
      ENDIF
C     WRITE(*,4449)'RCBORN',I,J,NFZ,NX,NY,L,K,M
C    2,RCBORN(M,K,L,NY,NX),ORN(M,K,L,NY,NX)
C    4,RCGSK(L,NY,NX),VCOXFS(NY,NX),VNOXFS(NY,NX),VPOXFS(NY,NX)
2840  CONTINUE
C
C     COMBUST DOC, DOA, DON, DOP
C
C     FRCBCO=combustion fraction of total SOC
C     SPCMB=specific combustion rate of SOC(K=0,4) at 600K,
C        charcoal(K=5) at 700K (g C g C-1 h-1)
C     TFNCO,TFNCH=temperature function for combusting SOC (600K = 1)
C        charcoal (700K = 1)
C     VOLX=soil layer volume (m3)
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate (g C,N,P,C) 
C     RCBOQC,RCBOQA,RCBOQN,RCBOQP=combustion rate of OQC,OQA,OQN,OQP
C        (g C,N,P,C t-1) 
C     RCGSK=total combustion for use in ‘trnsfr.f’ (g C -1)
C     VCOXFS,VNOXFS,VPOXFS=total C,N,P combustion loss for output
C        (g C,N,P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C     
      IF(OQC(K,L,NY,NX).GT.ZEROS(NY,NX) 
     2.AND.TKS(L,NY,NX).GT.TCMBX)THEN
      RCBOQC(K,L,NY,NX)=FRCBCO(K)*OQC(K,L,NY,NX)
      RCBOQA(K,L,NY,NX)=RCBOQC(K,L,NY,NX)
     2*OQA(K,L,NY,NX)/OQC(K,L,NY,NX)
      RCBOQN(K,L,NY,NX)=RCBOQC(K,L,NY,NX)
     2*OQN(K,L,NY,NX)/OQC(K,L,NY,NX)
      RCBOQP(K,L,NY,NX)=RCBOQC(K,L,NY,NX)
     2*OQP(K,L,NY,NX)/OQC(K,L,NY,NX)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-RCBOQC(K,L,NY,NX)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-RCBOQA(K,L,NY,NX)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-RCBOQN(K,L,NY,NX)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-RCBOQP(K,L,NY,NX)
      RCGSK(L,NY,NX)=RCGSK(L,NY,NX)+RCBOQC(K,L,NY,NX)
      VCOXFS(NY,NX)=VCOXFS(NY,NX)-RCBOQC(K,L,NY,NX)
      VNOXFS(NY,NX)=VNOXFS(NY,NX)-RCBOQN(K,L,NY,NX)
      VPOXFS(NY,NX)=VPOXFS(NY,NX)-RCBOQP(K,L,NY,NX)
      ELSE
      RCBOQC(K,L,NY,NX)=0.0
      RCBOQA(K,L,NY,NX)=0.0
      RCBOQN(K,L,NY,NX)=0.0
      RCBOQP(K,L,NY,NX)=0.0
      ENDIF
C     WRITE(*,4449)'RCBOQN',I,J,NFZ,NX,NY,L,K,M
C    2,RCBOQN(K,L,NY,NX),OQN(K,L,NY,NX)
C    4,RCGSK(L,NY,NX),VCOXFS(NY,NX),VNOXFS(NY,NX),VPOXFS(NY,NX)
C
C     COMBUST ADSORBED OM 
C
C     FRCBCO=combustion fraction of total SOC
C     SPCMB=specific combustion rate of SOC(K=0,4) at 600K,
C        charcoal(K=5) at 700K (g C g C-1 h-1)
C     TFNCO,TFNCH=temperature function for combusting SOC (600K = 1)
C        charcoal (700K = 1)
C     VOLX=soil layer volume (m3)
C     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate (g C,N,P,C)
C     RCBOHC,RCBOHA,RCBOHN,RCBOHP=combustion rate of OHC,OHA,OHN,OHP
C        (g C,N,P,C t-1)
C     RCGSK=total combustion for use in ‘trnsfr.f’ (g C t-1)
C     VCOXFS,VNOXFS,VPOXFS=total C,N,P combustion loss for output
C        (g C,N,P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      IF(OHC(K,L,NY,NX).GT.ZEROS(NY,NX) 
     2.AND.TKS(L,NY,NX).GT.TCMBX)THEN
      RCBOHC(K,L,NY,NX)=FRCBCO(K)*OHC(K,L,NY,NX)
      RCBOHA(K,L,NY,NX)=RCBOHC(K,L,NY,NX)
     2*OHA(K,L,NY,NX)/OHC(K,L,NY,NX)
      RCBOHN(K,L,NY,NX)=RCBOHC(K,L,NY,NX)
     2*OHN(K,L,NY,NX)/OHC(K,L,NY,NX)
      RCBOHP(K,L,NY,NX)=RCBOHC(K,L,NY,NX)
     2*OHP(K,L,NY,NX)/OHC(K,L,NY,NX)
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-RCBOHC(K,L,NY,NX)
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-RCBOHA(K,L,NY,NX)
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-RCBOHN(K,L,NY,NX)
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-RCBOHP(K,L,NY,NX)
      RCGSK(L,NY,NX)=RCGSK(L,NY,NX)+RCBOHC(K,L,NY,NX)
      VCOXFS(NY,NX)=VCOXFS(NY,NX)-RCBOHC(K,L,NY,NX)
      VNOXFS(NY,NX)=VNOXFS(NY,NX)-RCBOHN(K,L,NY,NX)
      VPOXFS(NY,NX)=VPOXFS(NY,NX)-RCBOHP(K,L,NY,NX)
      ELSE
      RCBOHC(K,L,NY,NX)=0.0
      RCBOHA(K,L,NY,NX)=0.0
      RCBOHN(K,L,NY,NX)=0.0
      RCBOHP(K,L,NY,NX)=0.0
      ENDIF
C     WRITE(*,4449)'RCBOHN',I,J,NFZ,NX,NY,L,K,M
C    2,RCBOHN(K,L,NY,NX),OHN(K,L,NY,NX)
C    4,RCGSK(L,NY,NX),VCOXFS(NY,NX),VNOXFS(NY,NX),VPOXFS(NY,NX)
C
C     COMBUST SOM
C
C     FRCBCO=combustion fraction of total SOC
C     SPCMB=specific combustion rate of SOC(K=0,4) at 600K,
C        charcoal(K=5) at 700K (g C g C-1 h-1)
C     TFNCO,TFNCH=temperature function for combusting SOC (600K = 1)
C        charcoal (700K = 1)
C     VOLX=soil layer volume (m3)
C     OSC,OSA,OSN,OSP=SOC,colonized SOC,SON,SOP (g C,C,N,P)
C     RCBOSC,RCBOSA,RCBOSN,RCBOSP=combustion rate of OSC,OSA,OSN,OSP
C        (g C,C,N,P t-1)
C     RCGSK=total combustion for use in ‘trnsfr.f’(g C t-1)
C     VCOXFS,VNOXFS,VPOXFS=total C,N,P combustion loss for output
C        (g C,N,P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      DO 2830 M=1,4
      IF(OSC(M,K,L,NY,NX).GT.ZEROS(NY,NX) 
     2.AND.TKS(L,NY,NX).GT.TCMBX)THEN
      RCBOSC(M,K,L,NY,NX)=FRCBCO(K)*OSC(M,K,L,NY,NX) 
      RCBOSA(M,K,L,NY,NX)=RCBOSC(M,K,L,NY,NX)
     2*OSA(M,K,L,NY,NX)/OSC(M,K,L,NY,NX)
      RCBOSN(M,K,L,NY,NX)=RCBOSC(M,K,L,NY,NX)
     2*OSN(M,K,L,NY,NX)/OSC(M,K,L,NY,NX)
      RCBOSP(M,K,L,NY,NX)=RCBOSC(M,K,L,NY,NX)
     2*OSP(M,K,L,NY,NX)/OSC(M,K,L,NY,NX)
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-RCBOSC(M,K,L,NY,NX)
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-RCBOSA(M,K,L,NY,NX)
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-RCBOSN(M,K,L,NY,NX)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-RCBOSP(M,K,L,NY,NX)
      RCGSK(L,NY,NX)=RCGSK(L,NY,NX)+RCBOSC(M,K,L,NY,NX)
      VCOXFS(NY,NX)=VCOXFS(NY,NX)-RCBOSC(M,K,L,NY,NX)
      VNOXFS(NY,NX)=VNOXFS(NY,NX)-RCBOSN(M,K,L,NY,NX)
      VPOXFS(NY,NX)=VPOXFS(NY,NX)-RCBOSP(M,K,L,NY,NX)
C     IF(L.EQ.0)THEN
C     WRITE(*,4449)'RCBOSN',I,J,NFZ,NX,NY,L,K,M
C    2,RCBOSN(M,K,L,NY,NX),OSN(M,K,L,NY,NX)
C    2,RCBOSC(M,K,L,NY,NX),OSC(M,K,L,NY,NX)
C    3,RCGSK(L,NY,NX),VCOXFS(NY,NX),VNOXFS(NY,NX),VPOXFS(NY,NX)
C    4,RCBCO,SPCMB(M),TFNCO,XNFH,TKS(L,NY,NX)
4449  FORMAT(A8,8I4,30E12.4)
C     ENDIF 
      ELSE
      RCBOSC(M,K,L,NY,NX)=0.0 
      RCBOSA(M,K,L,NY,NX)=0.0
      RCBOSN(M,K,L,NY,NX)=0.0
      RCBOSP(M,K,L,NY,NX)=0.0
      ENDIF
2830  CONTINUE
2800  CONTINUE
C
C     COMBUST CHARCOAL
C
C     FRCBCO=combustion fraction of total charcoal
C     SPCMBH=specific combustion rate of charcoal at 700K
C     TFNCH=temperature function for combusting charcoal (700K = 1)
C     VOLX=soil layer volume (m3)
C     OSC,OSA,OSN,OSP=charcoal SOC,colonized SOC,SON,SOP (g C,C,N,P)
C     RCBOSC,RCBOSA,RCBOSN,RCBOSP=combustion rate of OSC,OSA,OSN,OSP
C        (g C,C,N,P t-1)
C     RCGSK=total combustion for use in ‘trnsfr.f’ (g C t-1)
C     VCOXFS,VNOXFS,VPOXFS=total C,N,P combustion loss for output
C        (g C,N,P t-1)
C     XNFH=time step for biological fluxes from ‘wthr.f’ (h t-1)    
C
      TFNCH=AMIN1(TFNCX,EXP(20.620-120000/RTK))
      IF(ORGCC(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FRCBCH=AMIN1(1.0,SPCMBH*TFNCH*XNFH/ORGCC(L,NY,NX))
      ELSE
      FRCBCH=0.0
      ENDIF
      DO 2810 K=0,4
      IF(OSC(5,K,L,NY,NX).GT.ZEROS(NY,NX) 
     2.AND.TKS(L,NY,NX).GT.TCMBX)THEN
      RCBOSC(5,K,L,NY,NX)=FRCBCH*OSC(5,K,L,NY,NX) 
      RCBOSA(5,K,L,NY,NX)=RCBOSC(5,K,L,NY,NX)
     2*OSA(5,K,L,NY,NX)/OSC(5,K,L,NY,NX) 
      RCBOSN(5,K,L,NY,NX)=RCBOSC(5,K,L,NY,NX)
     2*OSN(5,K,L,NY,NX)/OSC(5,K,L,NY,NX)
      RCBOSP(5,K,L,NY,NX)=RCBOSC(5,K,L,NY,NX)
     2*OSP(5,K,L,NY,NX)/OSC(5,K,L,NY,NX)
      OSC(5,K,L,NY,NX)=OSC(5,K,L,NY,NX)-RCBOSC(5,K,L,NY,NX)
      OSA(5,K,L,NY,NX)=OSA(5,K,L,NY,NX)-RCBOSA(5,K,L,NY,NX)
      OSN(5,K,L,NY,NX)=OSN(5,K,L,NY,NX)-RCBOSN(5,K,L,NY,NX)
      OSP(5,K,L,NY,NX)=OSP(5,K,L,NY,NX)-RCBOSP(5,K,L,NY,NX)
      RCGSK(L,NY,NX)=RCGSK(L,NY,NX)+RCBOSC(5,K,L,NY,NX)
      VCOXFS(NY,NX)=VCOXFS(NY,NX)-RCBOSC(5,K,L,NY,NX)
      VNOXFS(NY,NX)=VNOXFS(NY,NX)-RCBOSN(5,K,L,NY,NX)
      VPOXFS(NY,NX)=VPOXFS(NY,NX)-RCBOSP(5,K,L,NY,NX)
      ELSE
      RCBOSC(5,K,L,NY,NX)=0.0 
      RCBOSN(5,K,L,NY,NX)=0.0
      RCBOSP(5,K,L,NY,NX)=0.0
      ENDIF
C     WRITE(*,4449)'RCBOSC',I,J,NFZ,NX,NY,L,K,M
C    2,RCBOSC(5,K,L,NY,NX),OSC(5,K,L,NY,NX),FRCBCH
C    3,RCGSK(L,NY,NX),VCOXFS(NY,NX),VNOXFS(NY,NX),VPOXFS(NY,NX)
C    4,SPCMBH,TFNCH,XNFH,ORGCC(L,NY,NX)
2810  CONTINUE
2850  CONTINUE
      ENDIF
9990  CONTINUE
9995  CONTINUE
      RETURN
      END



