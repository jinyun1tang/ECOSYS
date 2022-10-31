
      SUBROUTINE grosub(I,J,NFZ,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE CALCULATES ALL PLANT BIOLOGICAL TRANSFORMATIONS
C
      include "parameters.h"
      include "files.h"
      include "blkc.h"
      include "blk1cp.h"
      include "blk1cr.h"
      include "blk1g.h"
      include "blk1n.h"
      include "blk1p.h"
      include "blk1s.h"
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
      include "blk11a.h"
      include "blk11b.h"
      include "blk12a.h"
      include "blk12b.h"
      include "blk13a.h"
      include "blk13b.h"
      include "blk13c.h"
      include "blk14.h"
      include "blk16.h"
      include "blk18a.h"
      include "blk18b.h"
      DIMENSION PART(7),TFN6(JZ),ARSTKB(JB),NRX(2,JZ),ICHK1(2,JZ)
     2,SLA2(0:5),FXFB(0:5),FXFR(0:5),FXRT(0:1),FXSH(0:1),FXRN(6)
     3,FPART1(0:5),FPART2(0:5),CNDLI(6),CMOSC(4),WTLSBZ(JB),CPOOLZ(JB)
     4,ZPOOLZ(JB),PPOOLZ(JB),ZCX(JP,JY,JX),UPNFC(JP,JY,JX)
     5,FRSV(0:5),FXFY(0:1),FXFZ(0:1),NBZ(JB)
      DIMENSION RTNT(2),RLNT(2,JZ),RTSK1(2,JZ,10),RTSK2(2,JZ,10)
     2,RTDPL(10,JZ),FWTR(JZ),FWTB(JP),FRTDP(0:3)
     3,RCCZ(0:5),RCCY(0:5),RCCX(0:5),RCCQ(0:5)
     4,RCCZR(0:5),RCCYR(0:5),RCCXR(0:5),RCCQR(0:5)
     5,WTSHTA(JZ,JY,JX),FLG4Y(0:3),ATRPX(0:1),VMXG(0:1),RTSK(0:3)
     6,WGLFBL(JL,JB,JP,JY,JX)
      DIMENSION CH2O3(25),CH2O4(25),CPOOLK(JB,JP,JY,JX),FHVSTK(0:25)
     2,FHVSHK(0:25),WFNRG(2,JZ),WFNRR(2,JZ),PSILY(0:3),FVRN(0:5)
      DIMENSION FWOOD(0:1,JZ),FWOODN(0:1,JZ),FWOODP(0:1,JZ)
     2,FWODB(0:1,JZ),FWODLN(0:1,JZ),FWODLP(0:1,JZ),FWODSN(0:1,JZ)
     3,FWODSP(0:1,JZ),FWODR(0:1,JZ),FWODRN(0:1,JZ),FWODRP(0:1,JZ)
      DIMENSION TWPOOL(JY,JX),TWTLF(JY,JX),TWTSHE(JY,JX)
     2,TWTSTK(JY,JX),TWTRSV(JY,JX),TWTHSK(JY,JX),TWTEAR(JY,JX)
     3,TWTGR(JY,JX),TWTNDC(JY,JX),TWPOLN(JY,JX),TWTRVC(JY,JX)
     4,TWTSTG(JY,JX),TWPOLR(JZ,JY,JX),TWTRTL(JZ,JY,JX)
     5,TWPODL(JZ,JY,JX),TWTNDL(JZ,JY,JX),SPCMB(5)
C
C     PART1,PART2=leaf,sheath partitioning coefficients at emergence
C     PART1X,PART2X=minimum leaf,sheath partitioning coefficients
C     VMXC=rate constant for nonstructural C oxidation in respiration
C       (h-1)
C     FSNR=rate constant for litterfall at end of growing season (h-1)
C     FLG4X=number of hours with no grain filling required for
C        physiological maturity
C     FLGZX=number of hours until full senescence after
C        physiological maturity
C     XFRX=maximum storage C concentration for remobilization from
C        stalk,root reserves
C     XFRY=rate constant for physiological to storage from
C        stalk,root reserves (h-1)
C     IFLGQX=required hours after physiological maturity
C        until start of litterfall
C     FSNK=minimum ratio of branch or mycorrhizae to root for
C        calculating shoot-root C transfer
C     FMYC=rate constant for root-mycorrhizal C,N,P exchange (h-1)
C
      PARAMETER(PART1=0.75,PART2=0.25,PART1X=0.015,PART2X=0.005
     2,VMXC=0.015,FSNR=2.884E-03,FLG4X=168.0
     3,FLGZX=480.0,XFRX=2.5E-02,XFRY=2.5E-03,IFLGQX=960
     4,FSNK=0.05,FMYC=0.1,TFNCX=1.0)
C
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling (g N,P g-1 C)
C     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
C     PSILM=minimum water potential for organ expansion,extension
C        (MPa)
C     RCMN=minimum stomatal resistance to CO2 (s m-1)
C     RTDPX=distance behind growing point for secondary roots (m)
C     RTLGAX=minimum average secondary root length (m)
C     EMODR=root modulus of elasticity (MPa)
C
      PARAMETER(CNKI=1.0E-01,CPKI=1.0E-02)
      PARAMETER(RMPLT=0.010,PSILM=0.1,RCMN=1.560E+01,RTDPX=0.00
     2,RTLGAX=1.0E-03,EMODR=5.0)
C
C     QNTM=quantum efficiency (umol e- umol-1 PAR)
C     CURV=shape parameter for e- transport response to PAR
C     ELEC3,ELEC4=e- requirement for CO2 fixn by rubisco,PEP
C        carboxylase (umol e- umol CO2)
C     CO2KI=Ki for C3 leakage from bundle sheath to mesophyll
C        in C4 (uM)
C     FCO2B,FHCOB=partition decarboxylation,leakage to CO2,HCO3 in C4
C     COMP4=C4 CO2 compensation point (uM)
C     FDML=leaf water content (g H2O g-1 C)
C     FBS,FMP=leaf water content in bundle sheath, mesophyll
C        in C4 CO2 fixn
C
      PARAMETER(QNTM=0.45,CURV=0.70,CURV2=2.0*CURV,CURV4=4.0*CURV
     2,ELEC3=4.5,ELEC4=3.0,CO2KI=1.0E+03,FCO2B=0.02,FHCOB=1.0-FCO2B)
      PARAMETER(COMP4=0.5,FDML=6.0,FBS=0.2*FDML,FMP=0.8*FDML)
C
C     ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
C     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
C     GY,GZ=partitioning of grazed material to removal,respiration
C
      PARAMETER(ZPLFM=0.333,ZPLFD=1.0-ZPLFM,ZPGRM=0.750
     2,ZPGRD=1.0-ZPGRM,GY=1.0,GZ=1.0-GY)
      PARAMETER(YSTK=2.5E-03)
C
C     SETC,SETN,SETP=Km for nonstructural C,N,P concn on seed set
C        (g g-1)
C     SSL2,SNL2=parameter for calculating petiole and internode
C        extension vs petiole, internode growth
C     CNMX,CPMX,CNMN,CPMN=max,min N:C,P:C for nonstructural C,N,P
C        transfers
C
      PARAMETER(SETC=1.0E-02,SETN=1.0E-03,SETP=1.0E-04,TCMBX=505.0)
      PARAMETER(SSL2=-0.50,SNL2=-0.625,PPM=0.375)
      PARAMETER(CNMX=0.20,CPMX=0.020,CNMN=0.050,CPMN=0.005)
C
C     EN2F=N fixation yield from C oxidation (g N g-1 C)
C     VMXO=specific respiration rate by bacterial N2 fixers
C        (g g-1 h-1)
C     SPNDL=specific decomposition rate by canopy,root bacterial
C        N2 fixers (g g-1 h-1)
C     WTNDI=initial bacterial mass at infection (g C m-2)
C     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
C        (g N,P g-1 C)
C     RCCZN,RCCYN=min,max fractions for bacteria C recycling
C     RCCXN,RCCQN=max fractions for bacteria N,P recycling
C     RCCZ,RCCY=min,max fractions for shoot,bacteria C recycling
C     RCCX,RCCQ=max fractions for shoot,bacteria N,P recycling
C     RCCZR,RCCYR=min,max fractions for root C recycling
C     RCCXR,RCCQR=max fractions for root N,P recycling
C     FXFB=rate constant for leaf-storage nonstructural C,N,P exchange
C          from IBTYP in PFT file (h-1)
C     FXFR=rate constant for root-storage nonstructural C,N,P exchange
C          from IBTYP in PFT file (h-1)
C     FPART1,FPART2=phenology effect on reducing PART1,PART2
C     SLA2=parameter for calculating leaf area expansion
C        vs leaf growth
C
      PARAMETER(GN2X=150.0,EOMC=37.5,EN2F=EOMC/GN2X
     2,VMXO=0.125,WTNDI=1.0E-04,CZKM=1.0E-04
     3,CPKM=1.0E-05,SPNDL=1.0E-02,ZCKI=1.0E+01,ZPKI=1.0E+03
     4,RCCZN=0.167,RCCYN=0.167,RCCXN=0.833,RCCQN=0.833)
      DATA RCCZ/0.167,0.167,0.167,0.167,0.167,0.167/
      DATA RCCY/0.167,0.167,0.167,0.167,0.167,0.167/
      DATA RCCX/0.833,0.833,0.833,0.833,0.833,0.833/
      DATA RCCQ/0.833,0.833,0.833,0.833,0.833,0.833/
      DATA RCCZR/0.167,0.167,0.167,0.167,0.167,0.167/
      DATA RCCYR/0.167,0.167,0.167,0.167,0.167,0.167/
      DATA RCCXR/0.833,0.833,0.833,0.833,0.833,0.833/
      DATA RCCQR/0.833,0.833,0.833,0.833,0.833,0.833/
      DATA FXFB/5.0E-03,5.0E-03,5.0E-06,5.0E-06,5.0E-05,2.5E-04/
      DATA FXFR/2.0E-05,2.0E-05,2.0E-05,2.0E-05,2.0E-05,2.0E-05/
      DATA FPART1/0.75,1.00,1.50,1.50,1.25,1.00/
      DATA FPART2/0.25,0.33,0.50,0.50,0.42,0.33/
      DATA SLA2/-0.333,-0.333,-0.083,-0.333,-0.333,-0.333/
C
C     RTSK=relative primary root sink strength 0.25=shallow,
C        4.0=deep root profile
C     FXRN=rate constant for plant-bacteria nonstructural C,N,P
C        exchange (h-1)
C     CNDLI=bacteria:plant ratio controlling bacterial decomposition
C     FXSH,FXRT=shoot-root partitioning of storage C during leafout
C     FRSV=rate constant for remobilization of storage C,N,P during
C        leafout (h-1)
C     FXFY,FXFZ=rate constant for leaf-reserve nonstructural C,N,P
C        exchange (h-1)
C     PSILY=canopy water potential below which leafoff is induced
C       (MPa)
C     FLG4Y=number of hours after physiolical maturity required for
C        senescence
C     ATRPX=number of hours required to initiate remobilization of
C        storage C for leafout
C     VMXG=specific oxidation rate of nonstructural C during leafout
C        at 25 C
C     FVRN=fraction of hours required for leafoff to initiate
C        remobilization
C     SPCMB=specific combustion rate of living plant material (0-4)
C        at 600K, charcoal(5) at 700K (g C m-2 h-1)
C
      DATA RTSK/0.25,1.0,1.0,4.0/
      DATA FXRN/5.0E-02,2.5E-02,1.25E-02,5.0E-02,2.5E-02,1.25E-02/
      DATA CNDLI/10.0E-01,2.0E-01,2.0E-01,1.0E-01,0.2E-01,0.2E-01/
      DATA FXSH/0.25,0.75/,FXRT/0.75,0.25/
      DATA FRSV/0.025,0.025,0.0025,0.025,0.025,0.025/
      DATA FXFY/0.125E-02,0.125E-02/,FXFZ/0.10E-01,0.10E-01/
      DATA PSILY/-200.0,-2.0,-2.0,-2.0/
      DATA FLG4Y/360.0,1440.0,720.0,720.0/
      DATA ATRPX/45.8,138.4/,VMXG/0.015,0.005/
      DATA FVRN/0.75,0.5,0.5,0.5,0.5,0.5/
      DATA SPCMB/2.5E+02,2.5E+02,2.5E+02,2.5E+02,2.5E+02/
C     DATA TC4,TLK/0.0,0.0/
      REAL*4 TFN5,WFNSG,WFNSC,WFNST,WFNSR,WFN4,WFNB
     2,WFNRT,WFNRG,WFNRR,FSNC2
      XNFHQ=XNFH**0.667
      DO 2995 NX=NHW,NHE
      DO 2990 NY=NVN,NVS
      DO 2985 NZ=1,NP(NY,NX)
C
C     TOTAL AGB FOR GRAZING IN LANDSCAPE SECTION
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     LSG=landscape grazing section number
C     WTSHTZ,WTSHTA=total,average biomass in landscape grazing section
C
      IF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      WTSHTZ=0
      NN=0
      DO 1995 NX1=NHW,NHE
      DO 1990 NY1=NVN,NVS
      IF(LSG(NZ,NY1,NX1).EQ.LSG(NZ,NY,NX))THEN
      IF(IFLGC(NZ,NY1,NX1).EQ.1)THEN
      WTSHTZ=WTSHTZ+WTSHT(NZ,NY1,NX1)
      NN=NN+1
      ENDIF
      ENDIF
1990  CONTINUE
1995  CONTINUE
      IF(NN.GT.0)THEN
      WTSHTA(NZ,NY,NX)=WTSHTZ/NN
      ELSE
      WTSHTA(NZ,NY,NX)=WTSHT(NZ,NY,NX)
      ENDIF
      ENDIF
2985  CONTINUE
2990  CONTINUE
2995  CONTINUE
C
C     INITIALIZE SENESCENCE ARRAYS
C
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      DO 9980 NZ=1,NP0(NY,NX)
      DO 1 L=0,NJ(NY,NX)
      DO 1 K=0,1
      DO 1 M=1,5
      CSNC(M,K,L,NZ,NY,NX)=0.0
      ZSNC(M,K,L,NZ,NY,NX)=0.0
      PSNC(M,K,L,NZ,NY,NX)=0.0
1     CONTINUE
      UPNFC(NZ,NY,NX)=0.0
      ZCX(NZ,NY,NX)=ZC(NZ,NY,NX)
      ZC(NZ,NY,NX)=0.0
9980  CONTINUE
C
C     TRANSFORMATIONS IN LIVING PLANT POPULATIONS
C
      DO 9985 NZ=1,NP(NY,NX)
C     IF(J.EQ.INT(ZNOON(NY,NX)))THEN
      XHVST=1.0
      WHVSBL=0.0
      WTHTH0=0.0
      WTHNH0=0.0
      WTHPH0=0.0
      WTHTH1=0.0
      WTHNH1=0.0
      WTHPH1=0.0
      WTHTH2=0.0
      WTHNH2=0.0
      WTHPH2=0.0
      WTHTH3=0.0
      WTHNH3=0.0
      WTHPH3=0.0
      WTHTH4=0.0
      WTHNH4=0.0
      WTHPH4=0.0
      WTHTR1=0.0
      WTHNR1=0.0
      WTHPR1=0.0
      WTHTR2=0.0
      WTHNR2=0.0
      WTHPR2=0.0
      WTHTR3=0.0
      WTHNR3=0.0
      WTHPR3=0.0
      WTHTR4=0.0
      WTHNR4=0.0
      WTHPR4=0.0
      WTHTX0=0.0
      WTHNX0=0.0
      WTHPX0=0.0
      WTHTX1=0.0
      WTHNX1=0.0
      WTHPX1=0.0
      WTHTX2=0.0
      WTHNX2=0.0
      WTHPX2=0.0
      WTHTX3=0.0
      WTHNX3=0.0
      WTHPX3=0.0
      WTHTX4=0.0
      WTHNX4=0.0
      WTHPX4=0.0
      WTHTG=0.0
      WTHNG=0.0
      WTHPG=0.0
C     ENDIF
      IF(IFLGC(NZ,NY,NX).EQ.1)THEN
      IF(IDTHP(NZ,NY,NX).EQ.0.OR.IDTHR(NZ,NY,NX).EQ.0)THEN
      DO 2 L=1,JC
      ARLFV(L,NZ,NY,NX)=0.0
      WGLFV(L,NZ,NY,NX)=0.0
      ARSTV(L,NZ,NY,NX)=0.0
2     CONTINUE
      DO 5 NR=1,NRT(NZ,NY,NX)
      DO 5 N=1,MY(NZ,NY,NX)
      NRX(N,NR)=0
      ICHK1(N,NR)=0
5     CONTINUE
      DO 9 N=1,MY(NZ,NY,NX)
      RTNT(N)=0.0
      DO 6 L=NU(NY,NX),NJ(NY,NX)
      WSRTL(N,L,NZ,NY,NX)=0.0
      RTN1(N,L,NZ,NY,NX)=0.0
      RTNL(N,L,NZ,NY,NX)=0.0
      RCO2M(N,L,NZ,NY,NX)=0.0
      RCO2N(N,L,NZ,NY,NX)=0.0
      RCO2A(N,L,NZ,NY,NX)=0.0
      RLNT(N,L)=0.0
      DO 6 NR=1,NRT(NZ,NY,NX)
      RTSK1(N,L,NR)=0.0
      RTSK2(N,L,NR)=0.0
6     CONTINUE
9     CONTINUE
C
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     WTSTK,WVSTK=stalk,sapwood mass
C     FWOOD,FWODB=C woody fraction in stalk,other organs
C        :0=woody,1=non-woody
C     CN*,CP*=N:C,P:C ratios in plant organs from PFT files
C     CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood
C        content
C        plant organ code:*LF=leaf,*SHE=petiole,*STK=stalk,*RT=root
C     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
C     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
C     FWOODN,FWOODP=N,P woody fraction in stalk:0=woody,1=non-woody
C
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1
     2.OR.WTSTK(NZ,NY,NX).LE.ZEROP(NZ,NY,NX))THEN
      FWODB(1,NZ)=1.0
      FWOOD(1,NZ)=1.0
      FWODR(1,NZ)=1.0
      ELSE
      FWODB(1,NZ)=1.0
      FWOOD(1,NZ)=WVSTK(NZ,NY,NX)/WTSTK(NZ,NY,NX)
      FWODR(1,NZ)=FWOOD(1,NZ)**0.167
      ENDIF
      FWODB(0,NZ)=1.0-FWODB(1,NZ)
      FWOOD(0,NZ)=1.0-FWOOD(1,NZ)
      FWODR(0,NZ)=1.0-FWODR(1,NZ)
C     WRITE(*,8822)'FWOOD',I,J,NFZ,NX,NY,NZ
C    2,IBTYP(NZ,NY,NX),IGTYP(NZ,NY,NX)
C    2,FWOOD(0,NZ),FWOOD(1,NZ),FWODB(0,NZ),FWODB(1,NZ)
C    3,FWODR(0,NZ),FWODR(1,NZ)
C    3,WVSTK(NZ,NY,NX),WTSTK(NZ,NY,NX),WTRVC(NZ,NY,NX)
8822  FORMAT(A8,8I4,12E12.4)
      CNLFW=FWODB(0,NZ)*CNSTK(NZ,NY,NX)+FWODB(1,NZ)*CNLF(NZ,NY,NX)
      CPLFW=FWODB(0,NZ)*CPSTK(NZ,NY,NX)+FWODB(1,NZ)*CPLF(NZ,NY,NX)
      CNSHW=FWODB(0,NZ)*CNSTK(NZ,NY,NX)+FWODB(1,NZ)*CNSHE(NZ,NY,NX)
      CPSHW=FWODB(0,NZ)*CPSTK(NZ,NY,NX)+FWODB(1,NZ)*CPSHE(NZ,NY,NX)
      CNRTW=FWODR(0,NZ)*CNSTK(NZ,NY,NX)+FWODR(1,NZ)*CNRT(NZ,NY,NX)
      CPRTW=FWODR(0,NZ)*CPSTK(NZ,NY,NX)+FWODR(1,NZ)*CPRT(NZ,NY,NX)
      FWODLN(1,NZ)=FWODB(1,NZ)
      FWODLP(1,NZ)=FWODB(1,NZ)
      FWODSN(1,NZ)=FWODB(1,NZ)
      FWODSP(1,NZ)=FWODB(1,NZ)
      FWOODN(1,NZ)=FWOOD(1,NZ)
      FWOODP(1,NZ)=FWOOD(1,NZ)
      FWODRN(1,NZ)=FWODR(1,NZ)
      FWODRP(1,NZ)=FWODR(1,NZ)
      FWODLN(0,NZ)=1.0-FWODLN(1,NZ)
      FWODLP(0,NZ)=1.0-FWODLP(1,NZ)
      FWODSN(0,NZ)=1.0-FWODSN(1,NZ)
      FWODSP(0,NZ)=1.0-FWODSP(1,NZ)
      FWOODN(0,NZ)=1.0-FWOODN(1,NZ)
      FWOODP(0,NZ)=1.0-FWOODP(1,NZ)
      FWODRN(0,NZ)=1.0-FWODRN(1,NZ)
      FWODRP(0,NZ)=1.0-FWODRP(1,NZ)
C
C     SHOOT AND ROOT TEMPERATURE FUNCTIONS FOR MAINTENANCE
C     RESPIRATION FROM TEMPERATURES WITH OFFSETS FOR THERMAL ADAPTATION
C
C     TKC,TKCM=canopy temperature,canopy temperature used in
C        temperature function
C     TKS,TKSO=soil temperature,soil temp used in temperature function
C     OFFST=shift in Arrhenius curve for thermal adaptation
C     TFN5,TFN6=temperature function for canopy,root maintenance
C        respiration (25oC =1)
C     8.3143,710.0=gas constant,enthalpy
C     62500,195000=energy of activation,low temperature inactivation
C        (KJ mol-1)
C
      TKCM=TKC(NZ,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKCM
      STK=710.0*TKCM
      ACTVM=1.0+EXP((197500-STK)/RTK)
      TFN5=AMIN1(1.0E+03,EXP(25.216-62500/RTK)/ACTVM)
      DO 7 L=NU(NY,NX),NJ(NY,NX)
      TKSO=TKS(L,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKSO
      STK=710.0*TKSO
      ACTVM=1.0+EXP((197500-STK)/RTK)
      TFN6(L)=AMIN1(1.0E+03,EXP(25.216-62500/RTK)/ACTVM)
7     CONTINUE
      GROGR=0.0
C
C     PRIMARY ROOT NUMBER
C
C     WTRTA=root mass per plant used to calculate primary root number
C     WTRT,PP=root mass,PFT population
C     XRTN1=multiplier for number of primary root axes
C
      WTRTA(NZ,NY,NX)=AMAX1(0.999992087*WTRTA(NZ,NY,NX)
     2,WTRT(NZ,NY,NX)/PP(NZ,NY,NX))
      XRTN1=AMAX1(1.0,WTRTA(NZ,NY,NX)**0.667)*PP(NZ,NY,NX)
C
C     WATER STRESS FUNCTIONS FOR EXPANSION AND GROWTH RESPIRATION
C     FROM CANOPY TURGOR
C
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     WFNST=turgor expansion,extension function
C     PSILG,PSILM=current,minimum canopy turgor potential for
C        expansion,extension
C     WFNSC=stomatal resistance function of canopy turgor
C     PSILT=canopy water potential
C     WFNSG=growth function of canopy water potential
C     WFNSR=expansion,extension function of canopy water potential
C
      WFNST=AMIN1(1.0,AMAX1(0.0,PSILG(NZ,NY,NX)-PSILM))
      IF(IGTYP(NZ,NY,NX).EQ.0)THEN
      WFNSC=1.0
      WFNSG=EXP(0.05*PSILT(NZ,NY,NX))
      WFNSR=WFNSG**0.25
      ELSE
      WFNSC=EXP(RCS(NZ,NY,NX)*PSILG(NZ,NY,NX))
      WFNSG=EXP(0.10*PSILT(NZ,NY,NX))
      WFNSR=WFNSG**0.25
      ENDIF
C
C     CALCULATE GROWTH OF EACH BRANCH
C
C     WTLFB,WTSHEB,WTLSB=leaf,petiole,leaf+petiole mass
C     IDTHB=branch living flag: 0=alive,1=dead
C
      DO 105 NB=1,NBR(NZ,NY,NX)
      WTLSB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX)
     2+WTSHEB(NB,NZ,NY,NX))
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
C
C     PARTITION GROWTH WITHIN EACH BRANCH FROM GROWTH STAGE
C
C     PART=organ partitioning fraction
C        :1=leaf
C        :2=sheath or petiole
C        :3=stalk
C        :4=reserve,
C        :5,6=reproductive organs
C        :7=grain
C
      ARSTKB(NB)=0.0
      TOTAL=0.0
      DO 10 N=1,7
      PART(N)=0.0
10    CONTINUE
C
C     IF BEFORE FLORAL INDUCTION
C
C     IDAY(2,=floral initiation date
C
      IF(IDAY(2,NB,NZ,NY,NX).EQ.0)THEN
      PART(1)=PART1
      PART(2)=PART2
C
C     IF BEFORE ANTHESIS
C
C     IDAY(6,=start of anthesis and setting final seed number
C     PART1,PART2=leaf,sheath partitioning coefficients at emergence
C     PART1X,PART2X=minimum leaf,sheath partitioning coefficients
C     FPART1,FPART2=phenology effect on reducing PART1,PART2
C     TGSTGI=total change in vegetative node number normalized for
C        maturity group
C
      ELSEIF(IDAY(6,NB,NZ,NY,NX).EQ.0)THEN
      PART(1)=AMAX1(PART1X
     2,PART1-FPART1(IBTYP(NZ,NY,NX))*TGSTGI(NB,NZ,NY,NX))
      PART(2)=AMAX1(PART2X
     2,PART2-FPART2(IBTYP(NZ,NY,NX))*TGSTGI(NB,NZ,NY,NX))
      PARTS=1.0-PART(1)-PART(2)
      IF(SNL1(NZ,NY,NX).GT.0.0)THEN
      PART(3)=0.60*PARTS
      PART(4)=0.30*PARTS
      ELSE
      PART(3)=0.00
      PART(4)=0.90*PARTS
      ENDIF
      PARTX=PARTS-PART(3)-PART(4)
      PART(5)=0.5*PARTX
      PART(6)=0.5*PARTX
C
C     IF BEFORE GRAIN FILLING, DETERMINATE OR INDETERMINATE
C
C     IDAY(7,=start of grain filling and setting max seed size
C     IDTYP=growth habit:0=determinate,1=indetermimate from PFT file
C     PART1,PART2=leaf,sheath partitioning coefficients at emergence
C     PART1X,PART2X=minimum leaf,sheath partitioning coefficients
C     FPART1,FPART2=phenology effect on reducing PART1,PART2
C     TGSTGF=total change in reproductive node number normalized for
C        maturity group
C
      ELSEIF(IDAY(7,NB,NZ,NY,NX).EQ.0)THEN
      IF(IDTYP(NZ,NY,NX).EQ.0)THEN
      PART(1)=0.0
      PART(2)=0.0
      ELSE
      PART(1)=AMAX1(PART1X,(PART1-FPART1(IBTYP(NZ,NY,NX)))
     2*(1.0-TGSTGF(NB,NZ,NY,NX)))
      PART(2)=AMAX1(PART2X,(PART2-FPART2(IBTYP(NZ,NY,NX)))
     2*(1.0-TGSTGF(NB,NZ,NY,NX)))
      ENDIF
      PARTS=1.0-PART(1)-PART(2)
      IF(SNL1(NZ,NY,NX).GT.0.0)THEN
      PART(3)=AMAX1(0.0,0.60*PARTS*(1.0-TGSTGF(NB,NZ,NY,NX)))
      PART(4)=AMAX1(0.0,0.30*PARTS*(1.0-TGSTGF(NB,NZ,NY,NX)))
      ELSE
      PART(3)=0.00
      PART(4)=AMAX1(0.0,0.90*PARTS*(1.0-TGSTGF(NB,NZ,NY,NX)))
      ENDIF
      PARTX=PARTS-PART(3)-PART(4)
      PART(5)=0.5*PARTX
      PART(6)=0.5*PARTX
C
C     DURING GRAIN FILLING, DETERMINATE OR INDETERMINATE
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IDTYP=growth habit:0=determinate,1=indetermimate from PFT file
C     PART1X,PART2X=minimum leaf,sheath partitioning coefficients
C
      ELSE
      IF(IDTYP(NZ,NY,NX).EQ.0)THEN
      PART(7)=1.0
      ELSE
      PART(1)=PART1X
      PART(2)=PART2X
      PARTS=1.0-PART(1)-PART(2)
      IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      IF(SNL1(NZ,NY,NX).GT.0.0)THEN
      PART(3)=0.125*PARTS
      PART(5)=0.125*PARTS
      PART(6)=0.125*PARTS
      PART(7)=0.625*PARTS
      ELSE
      PART(3)=0.000
      PART(5)=0.125*PARTS
      PART(6)=0.125*PARTS
      PART(7)=0.750*PARTS
      ENDIF
      ELSE
      IF(SNL1(NZ,NY,NX).GT.0.0)THEN
      PART(3)=0.75*PARTS
      PART(7)=0.25*PARTS
      ELSE
      PART(3)=0.00
      PART(7)=1.00*PARTS
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
C     IF AFTER GRAIN FILLING
C
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IDAY(10,=physiological maturity date
C
      IF(IBTYP(NZ,NY,NX).EQ.0.AND.IDAY(10,NB,NZ,NY,NX).NE.0)THEN
      IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      PART(4)=0.0
      PART(3)=0.0
      PART(7)=0.0
      ELSE
      PART(4)=PART(4)+PART(3)
      PART(3)=0.0
      PART(7)=0.0
      ENDIF
      ENDIF
C
C     REDIRECT FROM STALK TO STALK RESERVES IF RESERVES BECOME LOW
C
C     IDAY(2,=floral initiation date
C     WTRSVB,WVSTKB=stalk reserve,sapwood mass
C     XFRX=maximum storage C concentration for remobilization from
C        stalk,root reserves
C
      IF(IDAY(2,NB,NZ,NY,NX).NE.0)THEN
      IF(WTRSVB(NB,NZ,NY,NX).LT.XFRX*WVSTKB(NB,NZ,NY,NX))THEN
      DO 1020 N=1,7
      IF(N.NE.4)THEN
      PART(4)=PART(4)+0.10*PART(N)
      PART(N)=PART(N)-0.10*PART(N)
      ENDIF
1020  CONTINUE
C
C     REDIRECT FROM STALK RESERVES TO STALK IF RESERVES BECOME TOO
C     LARGE
C
      ELSEIF(WTRSVB(NB,NZ,NY,NX).GT.WVSTKB(NB,NZ,NY,NX))THEN
      IF(SNL1(NZ,NY,NX).GT.0.0)THEN
      PART(3)=PART(3)+PART(4)+PART(7)
      PART(4)=0.0
      PART(7)=0.0
      ELSE
      PART(1)=PART(1)+PART1*(PART(4)+PART(7))
      PART(2)=PART(2)+PART2*(PART(4)+PART(7))
      PART(4)=0.0
      PART(7)=0.0
      ENDIF
      ENDIF
      ENDIF
C
C     REDIRECT FROM LEAVES TO STALK IF LAI BECOMES TOO LARGE
C
C     ARLFP=PFT leaf area (m2 m-2)
C     PTRT=allocation to leaf+petiole used to simulate phenology
C        effect on shoot-root C transfer
C
      ARLFI=ARLFP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      IF(ARLFI.GT.5.0)THEN
      FPARTL=AMIN1(1.0,(ARLFI-5.0)/5.0)
      IF(SNL1(NZ,NY,NX).GT.0.0)THEN
      PART(3)=PART(3)+FPARTL*(PART(1)+PART(2))
      ELSE
      PART(4)=PART(4)+FPARTL*(PART(1)+PART(2))
      ENDIF
      PART(1)=(1.0-FPARTL)*PART(1)
      PART(2)=(1.0-FPARTL)*PART(2)
      ENDIF
      IF(NB.EQ.NB1(NZ,NY,NX))THEN
      PTRT=AMAX1(0.0,PART(1)+PART(2))
      ENDIF
C
C     DECIDUOUS LEAF FALL AFTER GRAIN FILL IN DETERMINATES,
C     AFTER AUTUMNIZATION IN INDETERMINATES, OR AFTER SUSTAINED
C     WATER STRESS
C
C     ISTYP=growth habit:0=annual,1=perennial
C     VRNF,VRNX=leafoff hours,hours required for leafoff
C     FVRN=fraction of hours required for leafoff to initiate
C        remobilization
C     IDAY(8,=end date for setting final seed number
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     IFLGY,IFLGZ=remobilization flags
C     FLGZ=control rate of remobilization
C     XNFH=time step from �wthr.f�
C
      IF((ISTYP(NZ,NY,NX).NE.0.AND.VRNF(NB,NZ,NY,NX)
     2.GE.FVRN(IWTYP(NZ,NY,NX))*VRNX(NB,NZ,NY,NX))
     3.OR.(ISTYP(NZ,NY,NX).EQ.0
     4.AND.IDAY(8,NB,NZ,NY,NX).NE.0))THEN
      IFLGZ=1
      IF(ISTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.0)THEN
      IFLGY=1
      FLGZ(NB,NZ,NY,NX)=FLGZ(NB,NZ,NY,NX)+1.0*XNFH
      ELSEIF((IWTYP(NZ,NY,NX).EQ.1.OR.IWTYP(NZ,NY,NX).EQ.3)
     2.AND.TCC(NZ,NY,NX).LT.TCX(NZ,NY,NX))THEN
      IFLGY=1
      FLGZ(NB,NZ,NY,NX)=FLGZ(NB,NZ,NY,NX)+1.0*XNFH
      ELSEIF(IWTYP(NZ,NY,NX).GE.2
     2.AND.PSILT(NZ,NY,NX).LT.PSILY(IGTYP(NZ,NY,NX)))THEN
      IFLGY=1
      FLGZ(NB,NZ,NY,NX)=FLGZ(NB,NZ,NY,NX)+1.0*XNFH
      ENDIF
      IF(ISTYP(NZ,NY,NX).NE.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      IF(SNL1(NZ,NY,NX).GT.0.0)THEN
      PART(3)=PART(3)+0.5*(PART(1)+PART(2))
      PART(4)=PART(4)+0.5*(PART(1)+PART(2))
      ELSE
      PART(4)=PART(4)+1.0*(PART(1)+PART(2))
      ENDIF
      PART(1)=0.0
      PART(2)=0.0
      ENDIF
      ELSE
      IFLGZ=0
      IFLGY=0
      FLGZ(NB,NZ,NY,NX)=0.0
      ENDIF
      IF(NB.EQ.NB1(NZ,NY,NX))THEN
      IFLGZR=IFLGZ
      IFLGYR=IFLGY
      FLGZR=FLGZ(NB,NZ,NY,NX)
      ENDIF
C
C     CHECK PARTITIONING COEFFICIENTS
C
      DO 1000 N=1,7
      PART(N)=AMAX1(0.0,PART(N))
      TOTAL=TOTAL+PART(N)
1000  CONTINUE
      IF(TOTAL.GT.ZERO)THEN
      DO 1010 N=1,7
      PART(N)=PART(N)/TOTAL
1010  CONTINUE
      ELSE
      DO 1015 N=1,7
      PART(N)=0.0
1015  CONTINUE
      ENDIF
C
C     SHOOT COEFFICIENTS FOR GROWTH RESPIRATION AND N,P CONTENTS
C     FROM GROWTH YIELDS ENTERED IN 'READQ', AND FROM PARTITIONING
C     COEFFICIENTS ABOVE
C
C     DM*B=C production vs nonstructural C consumption
C     CN*W,CP*W=N:C,P:C ratios in plant organs weighted for wood
C        content
C     plant organ code:*LF=leaf,*SHE=petiole,*STK=stalk,*RSV=stalk
C        reserve,*HSK=husk,*EAR=ear,*GR=grain from PFT file,*SH=shoot
C     DMSHT=branch C production vs nonstructural C consumption
C     DMSHD=branch C respiration vs nonstructural C consumption
C     CN*M,CP*M=minimum N,P production vs nonstructural C consumption
C     CNLFX,CPLFX=N,P production vs nonstructural C consumption
C        in leaf
C     CNSHX,CPSHX=N,P production vs nonstructural C consumption
C        in rest of shoot
C     ZPLFM=minimum N:C,P:C in leaves relative to maximum values
C        from PFT file
C     ZPLFD=1.0-ZPLFM
C
      IF(IDAY(1,NB,NZ,NY,NX).NE.0)THEN
      DMLFB=DMLF(NZ,NY,NX)
      DMSHB=DMSHE(NZ,NY,NX)
      CNLFB=CNLFW
      CNSHB=CNSHW
      CPLFB=CPLFW
      CPSHB=CPSHW
      ELSE
      DMLFB=DMRT(NZ,NY,NX)
      DMSHB=DMRT(NZ,NY,NX)
      CNLFB=CNRTW
      CNSHB=CNRTW
      CPLFB=CPRTW
      CPSHB=CPRTW
      ENDIF
      DMSHT=PART(1)*DMLFB+PART(2)*DMSHB+PART(3)*DMSTK(NZ,NY,NX)
     2+PART(4)*DMRSV(NZ,NY,NX)+PART(5)*DMHSK(NZ,NY,NX)
     3+PART(6)*DMEAR(NZ,NY,NX)+PART(7)*DMGR(NZ,NY,NX)
      DMSHD=1.0-DMSHT
      CNLFM=PART(1)*DMLFB*ZPLFM*CNLFB
      CPLFM=PART(1)*DMLFB*ZPLFM*CPLFB
      CNLFX=PART(1)*DMLFB*ZPLFD*CNLFB
      CPLFX=PART(1)*DMLFB*ZPLFD*CPLFB
      CNSHX=PART(2)*DMSHB*CNSHB
     2+PART(3)*DMSTK(NZ,NY,NX)*CNSTK(NZ,NY,NX)
     3+PART(4)*DMRSV(NZ,NY,NX)*CNRSV(NZ,NY,NX)
     4+PART(5)*DMHSK(NZ,NY,NX)*CNHSK(NZ,NY,NX)
     5+PART(6)*DMEAR(NZ,NY,NX)*CNEAR(NZ,NY,NX)
     6+PART(7)*DMGR(NZ,NY,NX)*CNRSV(NZ,NY,NX)
      CPSHX=PART(2)*DMSHB*CPSHB
     2+PART(3)*DMSTK(NZ,NY,NX)*CPSTK(NZ,NY,NX)
     3+PART(4)*DMRSV(NZ,NY,NX)*CPRSV(NZ,NY,NX)
     4+PART(5)*DMHSK(NZ,NY,NX)*CPHSK(NZ,NY,NX)
     5+PART(6)*DMEAR(NZ,NY,NX)*CPEAR(NZ,NY,NX)
     6+PART(7)*DMGR(NZ,NY,NX)*CPRSV(NZ,NY,NX)
C
C     TOTAL SHOOT STRUCTURAL N MASS FOR MAINTENANCE RESPIRATION
C
C     WTSHXN=shoot structural N mass
C     WTLFBN,WTSHBN,WTHSBN,WTEARN,WTFRBN=leaf,petiole,husk,ear,grain
C        N mass
C     CNSTK,WVSTKB=stalk N:C,sapwood mass
C     IDAY(10=date of physiological maturity
C
      WTSHXN=AMAX1(0.0,WTLFBN(NB,NZ,NY,NX)+WTSHBN(NB,NZ,NY,NX)
     2+CNSTK(NZ,NY,NX)*WVSTKB(NB,NZ,NY,NX))
      IF(IDAY(10,NB,NZ,NY,NX).EQ.0)THEN
      WTSHXN=WTSHXN+AMAX1(0.0,WTHSBN(NB,NZ,NY,NX)
     2+WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX))
      ENDIF
C
C     GROSS PRIMARY PRODUCTIVITY
C
C     IDAY(1,=emergence date
C     FDBK=N,P feedback inhibition on C3 CO2 fixation from �stomate.f�
C     SSIN=sine of solar angle from �hour1.f�
C     RADP=total PAR absorbed by canopy from �hour1.f�
C     CO2Q=canopy air CO2 concentration
C     WFNSC=stomatal resistance function of canopy turgor
C
      IF(IDAY(1,NB,NZ,NY,NX).NE.0)THEN
      IF(FDBK(NB,NZ,NY,NX).NE.0.0)THEN
      IF(SSIN(NY,NX).GT.0.0.AND.RADP(NZ,NY,NX).GT.0.0
     2.AND.CO2Q(NY,NX).GT.0.0)THEN
      CO2F=0.0
      CH2O=0.0
      IF(IGTYP(NZ,NY,NX).NE.0.OR.WFNSC.GT.0.0)THEN
C
C     FOR EACH NODE
C
      DO 100 K=1,25
      CH2O3(K)=0.0
      CH2O4(K)=0.0
      IF(ARLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C
C     C4 PHOTOSYNTHESIS
C
C     ARLF,ARLFL=leaf area
C     ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
C     VCGR4=PEP carboxylation rate unlimited by CO2
C
      IF(ICTYP(NZ,NY,NX).EQ.4.AND.VCGR4(K,NB,NZ,NY,NX).GT.0.0)THEN
C
C     FOR EACH CANOPY LAYER
C
      DO 110 L=JC,1,-1
      IF(ARLFL(L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C
C
C     FOR EACH LEAF AZIMUTH AND INCLINATION
C
      DO 115 N = 1,4
      DO 120 M = 1,4
C
C     CO2 FIXATION IN MESOPHYLL BY SUNLIT LEAVES
C
C     SURFX=unself-shaded leaf surface area
C     PAR=direct PAR flux from �hour1.f�
C
      IF(SURFX(N,L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      IF(PAR(N,M,L,NZ,NY,NX).GT.0.0)THEN
C
C     C4 CARBOXYLATION REACTIONS IN MESOPHYLL
C
C     QNTM=quantum efficiency
C     PAR=direct PAR flux from �hour1.f�
C     ETGR4=light saturated e- transport rate from stomate.f
C     ETLF4=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO4=light-limited PEP carboxylation rate
C     CBXN4=PEP caboxylation efficiency
C     VL=PEP carboxylation rate limited by light,CO2,N,P
C     VGRO4=PEP carboxylation rate limited by CO2 from stomate.f
C     FDBK4=N,P feedback inhibition on C4 CO2 fixation
C
      PARX=QNTM*PAR(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGR4(K,NB,NZ,NY,NX)
      ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ,NY,NX)))/CURV2
      EGRO4=ETLF4*CBXN4(K,NB,NZ,NY,NX)
      VL=AMIN1(VGRO4(K,NB,NZ,NY,NX),EGRO4)*FDBK4(K,NB,NZ,NY,NX)
C
C     STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
C
C     RS,RSL=leaf stomatal resistance to CO2 at zero,current water
C        potential
C     RCMN=minimum stomatal resistance to CO2 (s m-1)
C     RCMX=cuticular resistance to CO2 from �startq.f� (s m-1)
C     DCO2=difference between atmospheric and intercellular CO2
C        concentration (umol m-3)
C     GSL=leaf stomatal conductance (mol m-2 s-1)
C     WFNSC=stomatal resistance function of canopy turgor
C     FMOL=number of moles of air per m3
C
      IF(VL.GT.ZERO)THEN
      RS=AMIN1(RCMX(NZ,NY,NX),AMAX1(RCMN,DCO2(NZ,NY,NX)/VL))
      RSL=RS+(RCMX(NZ,NY,NX)-RS)*WFNSC
      GSL=1.0/RSL*FMOL(NZ,NY,NX)
C
C     EFFECT OF WATER DEFICIT IN MESOPHYLL
C
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     WFNB=non-stomatal effects of water stress on C4,C3 CO2 fixation
C
      IF(IGTYP(NZ,NY,NX).NE.0)THEN
      WFN4=RS/RSL
      WFNB=SQRT(RS/RSL)
      ELSE
      WFN4=WFNSG
      WFNB=WFNSG
      ENDIF
C
C     CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
C     EQUALS DIFFUSION IN MESOPHYLL
C
C     CO2I=intercellular,mesophyll CO2 concentration at zero water
C        potential
C     CO2X,CO2C=intercellular,mesophyll CO2 concentration during
C        convergence for CO2I
C     SCO2=solubility of CO2 (uM/(umol mol-1))
C     COMP4=C4 CO2 compensation point (uM)
C     CBXNX=PEP carboxylation efficiency
C     ELEC4=e- requirement for CO2 fixn by PEP carboxylase
C     VCGR4,VGROX=PEP carboxylation rate unlimited,limited by CO2
C     XKCO24=Km for VCMX4 from PFT file (uM)
C     EGROX=light-limited PEP carboxylation rate
C     ETLF4=light-limited e- transport rate
C     VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
C     VG=CO2 diffusion rate limited by water stress
C     GSL=leaf stomatal conductance (mol m-2 s-1)
C
      CO2X=CO2I(NZ,NY,NX)
      DO 125 NN=1,100
      CO2C=CO2X*SCO2(NZ,NY,NX)
      CO2Y=AMAX1(0.0,CO2C-COMP4)
      CBXNX=CO2Y/(ELEC4*CO2C+10.5*COMP4)
      VGROX=VCGR4(K,NB,NZ,NY,NX)*CO2Y/(CO2C+XKCO24(NZ,NY,NX))
      EGROX=ETLF4*CBXNX
      VL=AMIN1(VGROX,EGROX)*WFN4*FDBK4(K,NB,NZ,NY,NX)
      VG=(CO2Q(NY,NX)-CO2X)*GSL
      IF(VL+VG.GT.ZERO)THEN
      DIFF=(VL-VG)/(VL+VG)
      IF(ABS(DIFF).LT.0.005)GO TO 130
      VA=0.95*VG+0.05*VL
      CO2X=CO2Q(NY,NX)-VA/GSL
      ELSE
      VL=0.0
      GO TO 130
      ENDIF
125   CONTINUE
C
C     ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
C
C     CH2O4=total C4 CO2 fixation
C     SURFX=unself-shaded leaf surface area from �hour1.f�
C     TAUS=fraction of direct radiation transmitted from layer above
C        from �hour1.f�
C
130   CH2O4(K)=CH2O4(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX)
     2*TAUS(L+1,NY,NX)
C     IF(NB.EQ.1.AND.M.EQ.1.AND.N.EQ.3.AND.K.EQ.KLEAF(NB,NZ,NY,NX)
C    2.AND.(I/10)*10.EQ.I.AND.J.EQ.12)THEN
C     WRITE(*,4444)'VLD4',IYRC,I,J,NZ,L,M,N,K,VL,PAR(N,M,L,NZ,NY,NX)
C    2,PAR(N,M,L,NZ,NY,NX)*TAUS(L+1,NY,NX)+PARDIF(N,M,L,NZ,NY,NX)
C    3*TAU0(L+1,NY,NX)
C    2,TKC(NZ,NY,NX),CO2Q(NY,NX),ETGR4(K,NB,NZ,NY,NX)
C    3,CBXN4(K,NB,NZ,NY,NX),VGRO4(K,NB,NZ,NY,NX),EGRO
C    3,FDBK4(K,NB,NZ,NY,NX),CH2O4(K),WFN4,VGROX,EGROX
C    4,VCGR4(K,NB,NZ,NY,NX),CO2X,CO2C,CBXNX
C    5,RS,RSL,SURFX(N,L,K,NB,NZ,NY,NX)
4444  FORMAT(A8,8I4,40E12.4)
C     ENDIF
C
C     C3 CARBOXYLATION REACTIONS IN BUNDLE SHEATH OF C4 PLANTS
C
C     ETGRO=light saturated e- transport rate from �stomate.f�
C     ETLF=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO=light-limited rubisco carboxylation rate
C     CBXN=rubisco caboxylation efficiency
C     VL=rubisco carboxylation rate limited by light,CO2,N,P
C     VGRO=rubisco carboxylation rate limited by CO2 from �stomate.f�
C     FDBK=N,P feedback inhibition on C3 CO2 fixation from �stomate.f�
C
      PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
      ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
      EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
      VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*WFNB*FDBK(NB,NZ,NY,NX)
C
C     ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
C
C     CH2O3=total C3 CO2 fixation
C     SURFX=unself-shaded leaf surface area
C     TAUS=fraction of direct radiation transmitted from layer above
C
      CH2O3(K)=CH2O3(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX)
     2*TAUS(L+1,NY,NX)
C     IF(L.EQ.NC-1.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
C     WRITE(*,4445)'VLD3',IYRC,I,J,NZ,L,M,N,K,VL,PAR(N,M,L,NZ,NY,NX)
C    2,TKC(NZ,NY,NX),CO2Q(NY,NX),ETGRO(K,NB,NZ,NY,NX)
C    3,CBXN(K,NB,NZ,NY,NX),VGRO(K,NB,NZ,NY,NX),EGRO
C    3,FDBK(NB,NZ,NY,NX),WFNB
4445  FORMAT(A8,8I4,20E12.4)
C     ENDIF
      ENDIF
      ENDIF
C
C     CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
C
      IF(PARDIF(N,M,L,NZ,NY,NX).GT.0.0)THEN
C
C     C4 CARBOXYLATION REACTIONS IN MESOPHYLL
C
C     QNTM=quantum efficiency
C     PARDIF=diffuse PAR flux from �hour1.f�
C     ETGR4=light saturated e- transport rate from �stomate.f�
C     ETLF4=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO4=light-limited PEP carboxylation rate
C     CBXN4=PEP caboxylation efficiency
C     VL=PEP carboxylation rate limited by light,CO2,N,P
C     VGRO4=PEP carboxylation rate limited by CO2 from �stomate.f�
C     FDBK4=N,P feedback inhibition on C4 CO2 fixation
C
      PARX=QNTM*PARDIF(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGR4(K,NB,NZ,NY,NX)
      ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ,NY,NX)))/CURV2
      EGRO4=ETLF4*CBXN4(K,NB,NZ,NY,NX)
      VL=AMIN1(VGRO4(K,NB,NZ,NY,NX),EGRO4)*FDBK4(K,NB,NZ,NY,NX)
C
C     STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
C
C     RS,RSL=leaf stomatal resistance to CO2 at zero,current water
C        potential
C     RCMN=minimum stomatal resistance to CO2 (s m-1)
C     RCMX=cuticular resistance to CO2 from startq.f (s m-1)
C     DCO2=difference between atmospheric and intercellular CO2
C        concentration (umol m-3)
C     GSL=leaf stomatal conductance (mol m-2 s-1)
C     WFNSC=stomatal resistance function of canopy turgor
C     FMOL=number of moles of air per m3
C
      IF(VL.GT.ZERO)THEN
      RS=AMIN1(RCMX(NZ,NY,NX),AMAX1(RCMN,DCO2(NZ,NY,NX)/VL))
      RSL=RS+(RCMX(NZ,NY,NX)-RS)*WFNSC
      GSL=1.0/RSL*FMOL(NZ,NY,NX)
C
C     EFFECT OF WATER DEFICIT IN MESOPHYLL
C
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     WFN4,WFNB=non-stomatal effects of water stress on C4,C3 CO2
C        fixation
C     RS,RSL=leaf stomatal resistance to CO2 at zero,current water
C        potential
C     WFNSG=growth function of canopy water potential
C
      IF(IGTYP(NZ,NY,NX).NE.0)THEN
      WFN4=(RS/RSL)**1.00
      WFNB=SQRT(RS/RSL)
      ELSE
      WFN4=WFNSG
      WFNB=WFNSG
      ENDIF
C
C     CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
C     EQUALS DIFFUSION IN MESOPHYLL
C
C     CO2I=intercellular,mesophyll CO2 concentration at zero water
C        potential
C     CO2X,CO2C=intercellular,mesophyll CO2 concentration during
C        convergence for CO2I
C     SCO2=solubility of CO2 (uM/(umol mol-1))
C     COMP4=C4 CO2 compensation point (uM)
C     CBXNX=PEP caboxylation efficiency
C     ELEC4=e- requirement for CO2 fixn by PEP carboxylase
C     VCGR4,VGRO4=PEP carboxylation rate unlimited,limited by CO2
C     XKCO24=Km for VCMX4 from PFT file (uM)
C     EGROX=light-limited PEP carboxylation rate
C     ETLF4=light-limited e- transport rate
C     VL=PEP carboxylation rate limited by light,CO2,N,P,water stress
C     VG=CO2 diffusion rate limited by water stress
C     GSL=leaf stomatal conductance (mol m-2 s-1)
C
      CO2X=CO2I(NZ,NY,NX)
      DO 135 NN=1,100
      CO2C=CO2X*SCO2(NZ,NY,NX)
      CO2Y=AMAX1(0.0,CO2C-COMP4)
      CBXNX=CO2Y/(ELEC4*CO2C+10.5*COMP4)
      VGROX=VCGR4(K,NB,NZ,NY,NX)*CO2Y/(CO2C+XKCO24(NZ,NY,NX))
      EGROX=ETLF4*CBXNX
      VL=AMIN1(VGROX,EGROX)*WFN4*FDBK4(K,NB,NZ,NY,NX)
      VG=(CO2Q(NY,NX)-CO2X)*GSL
      IF(VL+VG.GT.ZERO)THEN
      DIFF=(VL-VG)/(VL+VG)
      IF(ABS(DIFF).LT.0.005)GO TO 140
      VA=0.95*VG+0.05*VL
      CO2X=CO2Q(NY,NX)-VA/GSL
      ELSE
      VL=0.0
      GO TO 140
      ENDIF
135   CONTINUE
C
C     ACCUMULATE C4 FIXATION PRODUCT IN MESOPHYLL
C
C     CH2O4=total C4 CO2 fixation
C     SURFX=unself-shaded leaf surface area from �hour1.f�
C     TAU0=fraction of diffuse radiation transmitted from layer above
C        from �hour1.f�
C
140   CH2O4(K)=CH2O4(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX)
     2*TAU0(L+1,NY,NX)
C     WRITE(*,4455)'VLB4',IYRC,I,J,NZ,L,M,N,K,VL,PAR(N,M,L,NZ,NY,NX)
C    2,TKC(NZ,NY,NX),CO2Q(NY,NX),ETGR4(K,NB,NZ,NY,NX)
C    3,CBXN4(K,NB,NZ,NY,NX),VGRO4(K,NB,NZ,NY,NX),EGRO
C    3,FDBK4(K,NB,NZ,NY,NX),CH2O4(K),WFN4,VGROX,EGROX
C    4,VCGR4(K,NB,NZ,NY,NX),CO2X,CO2C,CBXNX
C    5,RS,RSL,SURFX(N,L,K,NB,NZ,NY,NX)
4455  FORMAT(A8,8I4,40E12.4)
C
C     C3 CARBOXYLATION REACTIONS IN IN BUNDLE SHEATH OF C4 PLANTS
C
C     ETGRO=light saturated e- transport rate from stomate.f
C     ETLF=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO=light-limited rubisco carboxylation rate
C     CBXN=rubisco caboxylation efficiency
C     VL=rubisco carboxylation rate limited by light,CO2,N,P
C     VGRO=rubisco carboxylation rate limited by CO2 from stomate.f
C     FDBK=N,P feedback inhibition on C3 CO2 fixation
C
      PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
      ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
      EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
      VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*WFNB*FDBK(NB,NZ,NY,NX)
C
C     ACCUMULATE C3 FIXATION PRODUCT IN BUNDLE SHEATH
C
C     CH2O3=total C3 CO2 fixation
C     SURFX=unself-shaded leaf surface area from �hour1.f�
C     TAU0=fraction of diffuse radiation transmitted from layer above
C        from �hour1.f�
C
      CH2O3(K)=CH2O3(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX)
     2*TAU0(L+1,NY,NX)
C     IF(J.EQ.13.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
C     WRITE(*,4444)'VLB4',IYRC,I,J,NZ,L,K,VL,PARDIF(N,M,L,NZ,NY,NX)
C    2,RAPY,TKC(NZ,NY,NX),CO2Q(NY,NX),CO2X,FMOL(NZ,NY,NX)/GSL
C    3,VCGRO(K,NB,NZ,NY,NX),ETLF,FDBK(NB,NZ,NY,NX),WFNB
C     ENDIF
      ENDIF
      ENDIF
      ENDIF
120   CONTINUE
115   CONTINUE
      ENDIF
110   CONTINUE
      CO2F=CO2F+CH2O4(K)
      CH2O=CH2O+CH2O3(K)
C
C     C3 PHOTOSYNTHESIS
C
C     ICTYP=photosynthesis type:3=C3,4=C4
C     VCGRO=rubisco carboxylation rate unlimited by CO2
C     ARLFL=node leaf area in canopy layer
C
      ELSEIF(ICTYP(NZ,NY,NX).NE.4.AND.VCGRO(K,NB,NZ,NY,NX).GT.0.0)THEN
C
C     FOR EACH CANOPY LAYER
C
      DO 210 L=JC,1,-1
      IF(ARLFL(L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C
C     FOR EACH LEAF AZIMUTH AND INCLINATION
C
      DO 215 N=1,4
      DO 220 M=1,4
C
C     CO2 FIXATION BY SUNLIT LEAVES
C
C     SURFX=unself-shaded leaf surface area from �hour1.f�
C     PAR=direct PAR flux from �hour1.f�
C
      IF(SURFX(N,L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      IF(PAR(N,M,L,NZ,NY,NX).GT.0.0)THEN
C
C     C3 CARBOXYLATION REACTIONS IN MESOPHYLL
C
C     QNTM=quantum efficiency
C     PAR=direct PAR flux from �hour1.f�
C     ETGRO=light saturated e- transport rate from stomate.f
C     ETLF=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO=light-limited rubisco carboxylation rate
C     CBXN=rubisco caboxylation efficiency
C     VL=rubisco carboxylation rate limited by light,CO2,N,P
C     VGRO=rubisco carboxylation rate limited by CO2 from �stomate.f�
C     FDBK=N,P feedback inhibition on C4 CO2 fixation from �stomate.f�
C
      PARX=QNTM*PAR(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
      ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
      EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
      VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*FDBK(NB,NZ,NY,NX)
C
C     STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
C
C     RS,RSL=leaf stomatal resistance to CO2 at zero,current water
C        potential (s m-1)
C     RCMN=minimum stomatal resistance to CO2 (s m-1)
C     RCMX=cuticular resistance to CO2 from startq.f (s m-1)
C     DCO2=difference between atmospheric and intercellular CO2
C        concentration (umol m-3)
C     GSL=leaf stomatal conductance (mol m-2 s-1)
C     WFNSC=stomatal resistance function of canopy turgor
C     FMOL=number of moles of air per m3
C
      IF(VL.GT.ZERO)THEN
      RS=AMIN1(RCMX(NZ,NY,NX),AMAX1(RCMN,DCO2(NZ,NY,NX)/VL))
      RSL=RS+(RCMX(NZ,NY,NX)-RS)*WFNSC
      GSL=1.0/RSL*FMOL(NZ,NY,NX)
C
C     EFFECT OF WATER DEFICIT IN MESOPHYLL
C
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     WFNB=non-stomatal effects of water stress on CO2 fixation
C
      IF(IGTYP(NZ,NY,NX).NE.0)THEN
      WFNB=SQRT(RS/RSL)
      ELSE
      WFNB=WFNSG
      ENDIF
C
C     CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
C     EQUALS DIFFUSION IN MESOPHYLL
C
C     CO2I=intercellular,mesophyll CO2 concentration at zero water
C        potential
C     CO2X,CO2C=intercellular,mesophyll CO2 concentration during
C        convergence for CO2I
C     SCO2=solubility of CO2 (uM/(umol mol-1))
C     COMPL=C3 CO2 compensation point (uM)
C     CBXNX=rubisco carboxylation efficiency
C     ELEC3=e- requirement for CO2 fixn by rubisco
C     VCGRO,VGROX=rubisco carboxylation rate unlimited,limited by CO2
C     XKCO2O=Km for rubisco carboxylation
C     EGROX=light-limited rubisco carboxylation rate
C     ETLF=light-limited e- transport rate
C     VL=rubisco carboxylation rate limited by light,CO2,N,P,water
C        stress
C     VG=CO2 diffusion rate limited by water stress
C     GSL=leaf stomatal conductance (mol m-2 s-1)
C
      CO2X=CO2I(NZ,NY,NX)
      DO 225 NN=1,100
      CO2C=CO2X*SCO2(NZ,NY,NX)
      CO2Y=AMAX1(0.0,CO2C-COMPL(K,NB,NZ,NY,NX))
      CBXNX=CO2Y/(ELEC3*CO2C+10.5*COMPL(K,NB,NZ,NY,NX))
      VGROX=VCGRO(K,NB,NZ,NY,NX)*CO2Y/(CO2C+XKCO2O(NZ,NY,NX))
      EGROX=ETLF*CBXNX
      VL=AMIN1(VGROX,EGROX)*WFNB*FDBK(NB,NZ,NY,NX)
      VG=(CO2Q(NY,NX)-CO2X)*GSL
      IF(VL+VG.GT.ZERO)THEN
      DIFF=(VL-VG)/(VL+VG)
      IF(ABS(DIFF).LT.0.005)GO TO 230
      VA=0.95*VG+0.05*VL
      CO2X=CO2Q(NY,NX)-VA/GSL
      ELSE
      VL=0.0
      GO TO 230
      ENDIF
225   CONTINUE
C
C     ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
C
C     CH2O3=total C4 CO2 fixation
C     SURFX=unself-shaded leaf surface area from �hour1.f�
C     TAUS=fraction of direct radiation transmitted from layer above
C        from �hour1.f�
C
230   CH2O3(K)=CH2O3(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX)
     2*TAUS(L+1,NY,NX)
C     IF(J.EQ.20.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
C     WRITE(20,3335)'VLD',I,J,NFZ,NZ,L,M,N,K
C    2,VL,PAR(N,M,L,NZ,NY,NX),PARDIF(N,M,L,NZ,NY,NX)
C    3,SURFX(N,L,K,NB,NZ,NY,NX)
C    2,TKC(NZ,NY,NX),TKAM(NY,NX),CO2Q(NY,NX),CO2X,CO2C
C    3,FMOL(NZ,NY,NX)/GSL,VGROX,EGROX,ETLF,CBXNX
C    4,FDBK(NB,NZ,NY,NX),WFNB,PSILG(NZ,NY,NX)
C    4,VCGRO(K,NB,NZ,NY,NX),ETGRO(K,NB,NZ,NY,NX),COMPL(K,NB,NZ,NY,NX)
C    5,SURFX(N,L,K,NB,NZ,NY,NX),TAUS(L+1,NY,NX),CH2O3(K)
3335  FORMAT(A8,8I4,240E12.4)
C     ENDIF
      ENDIF
      ENDIF
C
C     CO2 FIXATION IN MESOPHYLL BY SHADED LEAVES
C
      IF(PARDIF(N,M,L,NZ,NY,NX).GT.0.0)THEN
C
C     C3 CARBOXYLATION REACTIONS USING VARIABLES FROM 'STOMATE'
C
C     QNTM=quantum efficiency
C     PARDIF=diffuse PAR flux from �hour1.f�
C     ETGRO=light saturated e- transport rate from �stomate.f�
C     ETLF=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO=light-limited rubisco carboxylation rate
C     CBXN=rubisco caboxylation efficiency
C     VL=rubisco carboxylation rate limited by light,CO2,N,P
C     VGRO=rubisco carboxylation rate limited by CO2 from �stomate.f�
C     FDBK=N,P feedback inhibition on C3 CO2 fixation from �stomate.f�
C
      PARX=QNTM*PARDIF(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
      ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
      EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
      VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*FDBK(NB,NZ,NY,NX)
C
C     STOMATAL EFFECT OF WATER DEFICIT IN MESOPHYLL
C
C     RS,RSL=leaf stomatal resistance to CO2 at zero,current water
C        potential (s m-1)
C     RCMN=minimum stomatal resistance to CO2 (s m-1)
C     RCMX=cuticular resistance to CO2 from startq.f (s m-1)
C     DCO2=difference between atmospheric and intercellular CO2
C        concentration (umol m-3)
C     GSL=leaf stomatal conductance (mol m-2 s-1)
C     WFNSC=stomatal resistance function of canopy turgor
C     FMOL=number of moles of air per m3
C
      IF(VL.GT.ZERO)THEN
      RS=AMIN1(RCMX(NZ,NY,NX),AMAX1(RCMN,DCO2(NZ,NY,NX)/VL))
      RSL=RS+(RCMX(NZ,NY,NX)-RS)*WFNSC
      GSL=1.0/RSL*FMOL(NZ,NY,NX)
C
C     EFFECT OF WATER DEFICIT IN MESOPHYLL
C
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     WFNB=non-stomatal effects of water stress on C3 CO2 fixation
C
      IF(IGTYP(NZ,NY,NX).NE.0)THEN
      WFNB=SQRT(RS/RSL)
      ELSE
      WFNB=WFNSG
      ENDIF
C
C     CONVERGENCE SOLUTION FOR CO2I AT WHICH CARBOXYLATION
C     EQUALS DIFFUSION IN MESOPHYLL
C
C     CO2I=intercellular,mesophyll CO2 concentration at zero water
C        potential
C     CO2X,CO2C=intercellular,mesophyll CO2 concentration during
C        convergence for CO2I
C     SCO2=solubility of CO2 (uM/(umol mol-1))
C     COMPL=C3 CO2 compensation point (uM)
C     CBXNX=rubisco caboxylation efficiency
C     ELEC3=e- requirement for CO2 fixation by rubisco carboxylase
C     VCGRO,VGROX=rubisco carboxylation rate unlimited,limited by CO2
C     XKCO2O=Km for rubisco carboxylation from stomate.f (uM)
C     EGROX=light-limited rubisco carboxylation rate
C     ETLF=light-limited e- transport rate
C     VL=rubisco carboxylation rate limited by light,CO2,N,P,water
C        stress
C     VG=CO2 diffusion rate limited by water stress
C     GSL=leaf stomatal conductance (mol m-2 s-1)
C
      CO2X=CO2I(NZ,NY,NX)
      DO 235 NN=1,100
      CO2C=CO2X*SCO2(NZ,NY,NX)
      CO2Y=AMAX1(0.0,CO2C-COMPL(K,NB,NZ,NY,NX))
      CBXNX=CO2Y/(ELEC3*CO2C+10.5*COMPL(K,NB,NZ,NY,NX))
      VGROX=VCGRO(K,NB,NZ,NY,NX)*CO2Y/(CO2C+XKCO2O(NZ,NY,NX))
      EGROX=ETLF*CBXNX
      VL=AMIN1(VGROX,EGROX)*WFNB*FDBK(NB,NZ,NY,NX)
      VG=(CO2Q(NY,NX)-CO2X)*GSL
      IF(VL+VG.GT.ZERO)THEN
      DIFF=(VL-VG)/(VL+VG)
      IF(ABS(DIFF).LT.0.005)GO TO 240
      VA=0.95*VG+0.05*VL
      CO2X=CO2Q(NY,NX)-VA/GSL
      ELSE
      VL=0.0
      GO TO 240
      ENDIF
235   CONTINUE
C
C     ACCUMULATE C3 FIXATION PRODUCT IN MESOPHYLL
C
C     CH2O3=total C3 CO2 fixation
C     SURFX=unself-shaded leaf surface area from �hour1.f�
C     TAU0=fraction of diffuse radiation transmitted from layer above
C        from �hour1.f�
C
240   CH2O3(K)=CH2O3(K)+VL*SURFX(N,L,K,NB,NZ,NY,NX)
     2*TAU0(L+1,NY,NX)
C     IF(J.EQ.20.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1)THEN
C     WRITE(*,3335)'VLB',I,J,NFZ,NZ,L,M,N,K
C    2,VL,PARDIF(N,M,L,NZ,NY,NX),SURFX(N,L,K,NB,NZ,NY,NX)
C    2,RAPY,TKC(NZ,NY,NX),TKAM(NY,NX)
C    3,CO2Q(NY,NX),CO2X,FMOL(NZ,NY,NX)/GSL
C    3,VCGRO(K,NB,NZ,NY,NX),ETLF,FDBK(NB,NZ,NY,NX),WFNB
C     ENDIF
      ENDIF
      ENDIF
      ENDIF
220   CONTINUE
215   CONTINUE
      ENDIF
210   CONTINUE
      CO2F=CO2F+CH2O3(K)
      CH2O=CH2O+CH2O3(K)
      ENDIF
      ENDIF
100   CONTINUE
C
C     CO2F,CH2O=total CO2 fixation,CH2O production
C
      CO2F=CO2F*0.0432*XNFH
      CH2O=CH2O*0.0432*XNFH
C
C     CONVERT UMOL M-2 S-1 TO G C M-2 H-1
C
      DO 150 K=1,25
      CH2O3(K)=CH2O3(K)*0.0432*XNFH
      CH2O4(K)=CH2O4(K)*0.0432*XNFH
150   CONTINUE
      ELSE
      CO2F=0.0
      CH2O=0.0
      IF(ICTYP(NZ,NY,NX).EQ.4)THEN
      DO 155 K=1,25
      CH2O3(K)=0.0
      CH2O4(K)=0.0
155   CONTINUE
      ENDIF
      ENDIF
      ELSE
      CO2F=0.0
      CH2O=0.0
      IF(ICTYP(NZ,NY,NX).EQ.4)THEN
      DO 160 K=1,25
      CH2O3(K)=0.0
      CH2O4(K)=0.0
160   CONTINUE
      ENDIF
      ENDIF
      ELSE
      CO2F=0.0
      CH2O=0.0
      IF(ICTYP(NZ,NY,NX).EQ.4)THEN
      DO 165 K=1,25
      CH2O3(K)=0.0
      CH2O4(K)=0.0
165   CONTINUE
      ENDIF
      ENDIF
C
C     SHOOT AUTOTROPHIC RESPIRATION AFTER EMERGENCE
C
C     N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
C
C     CNPG=N,P constraint on growth respiration
C     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concentration
C        in branch(g g-1)
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C
      IF(CCPOLB(NB,NZ,NY,NX).GT.ZERO)THEN
      CNPG=AMIN1(CZPOLB(NB,NZ,NY,NX)/(CZPOLB(NB,NZ,NY,NX)
     2+CCPOLB(NB,NZ,NY,NX)*CNKI)
     3,CPPOLB(NB,NZ,NY,NX)/(CPPOLB(NB,NZ,NY,NX)
     2+CCPOLB(NB,NZ,NY,NX)*CPKI))
      ELSE
      CNPG=1.0
      ENDIF
C
C     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
C     NON-STRUCTURAL C:N:P
C     MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
C
C     RCO2C=respiration from non-structural C
C     VMXC=rate constant for nonstructural C oxidation in respiration
C        (h-1)
C     CPOOL=non-structural C mass
C     TFN3=temperature function for canopy growth
C     WFNSG=growth function of canopy water potential
C     CNPG=N,P constraint on respiration
C     FDBKX=termination feedback inhibition on C3 CO2
C     RMNCS=maintenance respiration
C     RMPLT=specific maintenance respiration rate
C     TFN5=temperature function for canopy maintenance respiration
C     WTSHXN=shoot structural N mass
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     XNFH=time step from �wthr.f�
C
      IF(ISTYP(NZ,NY,NX).NE.0.OR.IDAY(10,NB,NZ,NY,NX).EQ.0)THEN
      RCO2C=AMAX1(0.0,VMXC*CPOOL(NB,NZ,NY,NX)
     2*TFN3(NZ,NY,NX)*CNPG*FDBKX(NB,NZ,NY,NX)*WFNSG*XNFH)
      RMNCS=AMAX1(0.0,RMPLT*TFN5*WTSHXN*XNFH)
      IF(IGTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.2)THEN
      RMNCS=RMNCS*WFNSG
      ELSE
      RMNCS=RMNCS*WFNSR
      ENDIF
      ELSE
      RCO2C=0.0
      RMNCS=0.0
      ENDIF
C
C     GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
C     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
C
C     RCO2X=difference between non-structural C respiration and
C        maintenance respiration
C     RCO2Y=growth respiration unlimited by N,P
C     WFNSR=expansion,extension function of canopy water potential
C     SNCR=excess maintenance respiration
C
      RCO2X=RCO2C-RMNCS
      RCO2Y=AMAX1(0.0,RCO2X)
      SNCR=AMAX1(0.0,-RCO2X)
C
C     GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
C     AVAILABLE FOR GROWTH
C
C     RCO2Y,RCO2G=growth respiration unlimited,limited by N,P
C     CNLFX,CPLFX=N,P production vs nonstructural C consumption
C        in leaf
C     CNSHX,CPSHX=N,P production vs nonstructural C consumption
C        in rest of shoot
C     ZPOOL,PPOOL=nonstructural N,P mass
C     DMSHD=branch C respiration vs nonstructural C consumption
C     CNLFM,CPLFM=minimum leaf N,P production vs nonstructural C
C        consumption
C     CNPG=N,P constraint on growth respiration
C
      IF(RCO2Y.GT.0.0.AND.(CNSHX.GT.0.0.OR.CNLFX.GT.0.0))THEN
      ZPOOLB=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
      PPOOLB=AMAX1(0.0,PPOOL(NB,NZ,NY,NX))
      RCO2G=AMIN1(RCO2Y,ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG)
     2,PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
      ELSE
      RCO2G=0.0
      ENDIF
C
C     TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
C     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
C     ENTERED IN 'READQ'
C
C     CGROS=total non-structural C used in growth and growth
C        respiration
C     RCO2G=growth respiration limited by N,P
C     DMSHD=branch C respiration vs nonstructural C consumption
C     ZADDB,PADDB=nonstructural N,P used in growth
C     ZPOOL,PPOOL=nonstructural N,P mass
C     CNSHX,CPSHX=N,P production vs nonstructural C consumption
C        in rest of shoot
C     CNLFM,CPLFM=minimum leaf N,P production vs nonstructural C
C        consumption
C     CNLFX,CPLFX=N,P production vs nonstructural C consumption
C        in leaf
C     CNPG=N,P constraint on growth respiration
C     CNRDA=respiration for N assimilation
C     CH2O=total CH2O production
C
      CGROS=RCO2G/DMSHD
      ZADDB=AMAX1(0.0,AMIN1(ZPOOL(NB,NZ,NY,NX)
     2,CGROS*(CNSHX+CNLFM+CNLFX*CNPG)))
      PADDB=AMAX1(0.0,AMIN1(PPOOL(NB,NZ,NY,NX)
     2,CGROS*(CPSHX+CPLFM+CPLFX*CNPG)))
      CNRDA=AMAX1(0.0,1.70*ZADDB-0.025*CH2O)
C
C     TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
C     ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
C
C     RCO2T=total C respiration
C     RMNCS=maintenance respiration
C     RCO2C=respiration from non-structural C
C     RCO2G=growth respiration limited by N,P
C     SNCR=excess maintenance respiration
C     CNRDA=respiration for N assimilation
C     CARBN=total PFT CO2 fixation
C     CO2F=total CO2 fixation
C     TCO2T,TCO2A=total,above-ground PFT respiration
C     CNET=PFT net CO2 fixation
C     TGPP=ecosystem GPP
C     RECO=ecosystem respiration
C     TRAU=total autotrophic respiration
C
      RCO2T=AMIN1(RMNCS,RCO2C)+RCO2G+SNCR+CNRDA
      CARBN(NZ,NY,NX)=CARBN(NZ,NY,NX)+CO2F
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)-RCO2T
      TCO2A(NZ,NY,NX)=TCO2A(NZ,NY,NX)-RCO2T
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)+CO2F-RCO2T
      HCNET(NZ,NY,NX)=HCNET(NZ,NY,NX)+CO2F-RCO2T
      TGPP(NY,NX)=TGPP(NY,NX)+CO2F
      RECO(NY,NX)=RECO(NY,NX)-RCO2T
      TRAU(NY,NX)=TRAU(NY,NX)-RCO2T
C     IF(NZ.EQ.1)THEN
C     WRITE(*,4477)'RCO2',I,J,NFZ,NX,NY,NZ,NB,IFLGZ
C    2,RCO2C,CNPG,CPOOL(NB,NZ,NY,NX),ZPOOL(NB,NZ,NY,NX)
C    3,PPOOL(NB,NZ,NY,NX),CCPOLB(NB,NZ,NY,NX)
C    3,CZPOLB(NB,NZ,NY,NX),CPPOLB(NB,NZ,NY,NX)
C    2,CO2F,CH2O,RMNCS,CGROS,CNRDA,RCO2T,RCO2X,SNCR
C    3,RCO2G,DMSHD,ZADDB,DMLFB,DMSHB
C    4,WTRSVB(NB,NZ,NY,NX),WVSTKB(NB,NZ,NY,NX),WTSHXN
C    5,PSILT(NZ,NY,NX)
C    6,ZADDB,RNH3B(NB,NZ,NY,NX),WFR(1,NG(NZ,NY,NX),NZ,NY,NX)
C    7,WFNSG,TFN3(NZ,NY,NX),TFN5,TKC(NZ,NY,NX),FDBK(NB,NZ,NY,NX)
C    8,FDBKX(NB,NZ,NY,NX),VMXC
4477  FORMAT(A8,8I4,40E12.4)
C     ENDIF
C
C     SHOOT AUTOTROPHIC RESPIRATION BEFORE EMERGENCE
C
      ELSE
C
C     N,P CONSTRAINT ON RESPIRATION FROM NON-STRUCTURAL C:N:P
C
C     CNPG=N,P constraint on growth respiration
C     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concentration in branch
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C
      IF(CCPOLB(NB,NZ,NY,NX).GT.ZERO)THEN
      CNPG=AMIN1(CZPOLB(NB,NZ,NY,NX)/(CZPOLB(NB,NZ,NY,NX)
     2+CCPOLB(NB,NZ,NY,NX)*CNKI)
     3,CPPOLB(NB,NZ,NY,NX)/(CPPOLB(NB,NZ,NY,NX)
     4+CCPOLB(NB,NZ,NY,NX)*CPKI))
      ELSE
      CNPG=1.0
      ENDIF
C
C     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
C     NON-STRUCTURAL C:N:P, O2 UPTAKE
C
C     RCO2CM,RCO2C=respiration from non-structural C unlimited,limited
C        by O2
C     VMXC=rate constant for nonstructural C oxidation in respiration
C       (h-1)
C     CPOOL=non-structural C mass
C     TFN4=temperature function for root growth
C     WFNSG=growth function of canopy water potential
C     CNPG=N,P constraint on respiration
C     FDBKX=termination feedback inhibition on C3 CO2
C     WFR=constraint by O2 consumption on all root processes
C
      RCO2CM=AMAX1(0.0,VMXC*CPOOL(NB,NZ,NY,NX)
     2*TFN4(NG(NZ,NY,NX),NZ,NY,NX))*CNPG*WFNSG*FDBKX(NB,NZ,NY,NX)*XNFH
      RCO2C=RCO2CM*WFR(1,NG(NZ,NY,NX),NZ,NY,NX)
C
C     MAINTENANCE RESPIRATION FROM TEMPERATURE, PLANT STRUCTURAL N
C
C     RMNCS=maintenance respiration
C     TFN6=temperature function for root maintenance respiration
C     WTSHXN=shoot structural N mass
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     WFNSG=growth function of canopy water potential
C
      RMNCS=AMAX1(0.0,RMPLT*TFN6(NG(NZ,NY,NX))*WTSHXN*XNFH)
      IF(IGTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.2)THEN
      RMNCS=RMNCS*WFNSG
      ELSE
      RMNCS=RMNCS*WFNSR
      ENDIF
C
C     GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
C     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
C
C     RCO2XM,RCO2X=diff between C respn unltd,ltd by O2 and mntc respn
C     RCO2YM,RCO2Y=growth respiration unltd,ltd by O2 and unlimited by
C        N,P
C     WFNSR=expansion,extension function of canopy water potential
C     SNCRM,SNCR=excess maintenance respiration unltd,ltd by O2
C
      RCO2XM=RCO2CM-RMNCS
      RCO2X=RCO2C-RMNCS
      RCO2YM=AMAX1(0.0,RCO2XM)
      RCO2Y=AMAX1(0.0,RCO2X)
      SNCRM=AMAX1(0.0,-RCO2XM)
      SNCR=AMAX1(0.0,-RCO2X)
C
C     GROWTH RESPIRATION MAY BE LIMITED BY NON-STRUCTURAL N,P
C     AVAILABLE FOR GROWTH
C
C     RCO2YM,RCO2Y=growth respiration unlimited,limited by O2
C        and unlimited by N,P
C     CNLFX,CPLFX=N,P production vs nonstructural C consumption
C        in leaf
C     CNSHX,CPSHX=N,P production vs nonstructural C consumption
C        in rest of shoot
C     ZPOOL,PPOOL=nonstructural N,P mass
C     FNP=growth respiration limited by O2 and N,P
C     DMSHD=branch C respiration vs nonstructural C consumption
C     CNLFM,CPLFM=minimum leaf N,P production vs nonstructural C
C        consumption
C     CNPG=N,P constraint on growth respiration
C     RCO2GM,RCO2G=growth respiration unlimited,limited by O2
C        and limited by N,P
C     WFR=constraint by O2 consumption on all root processes
C
      IF(CNSHX.GT.0.0.OR.CNLFX.GT.0.0)THEN
      ZPOOLB=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
      PPOOLB=AMAX1(0.0,PPOOL(NB,NZ,NY,NX))
      FNP=AMIN1(ZPOOLB*DMSHD/(CNSHX+CNLFM+CNLFX*CNPG)
     2,PPOOLB*DMSHD/(CPSHX+CPLFM+CPLFX*CNPG))
      IF(RCO2YM.GT.0.0)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
      ELSE
      RCO2GM=0.0
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
      RCO2G=AMIN1(RCO2Y,FNP*WFR(1,NG(NZ,NY,NX),NZ,NY,NX))
      ELSE
      RCO2G=0.0
      ENDIF
      ELSE
      RCO2GM=0.0
      RCO2G=0.0
      ENDIF
C
C     TOTAL NON-STRUCTURAL C,N,P USED IN GROWTH
C     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELDS
C     ENTERED IN 'READQ'
C
C     CGROSM,CGROS=total non-structural C used in growth and
C        respiration unlimited,limited by O2
C     RCO2GM,RCO2G=growth respiration unlimited,limited by O2 and
C        limited by N,P
C     DMSHD=branch C respiration vs nonstructural C consumption
C     ZADDBM,ZADDB,PADDB=nonstructural N,P unlimited,limited by O2
C        used in growth
C     ZPOOL,PPOOL=nonstructural N,P mass
C     CNSHX,CPSHX=N,P production vs nonstructural C consumption
C        in rest of shoot
C     CNLFM,CPLFM=minimum leaf N,P production vs nonstructural C
C        consumption
C     CNLFX,CPLFX=N,P production vs nonstructural C consumption
C        in leaf
C     CNPG=N,P constraint on growth respiration
C     CNRDM,CNRDA=respiration for N assimilation unlimited,limited
C        by O2
C
      CGROSM=RCO2GM/DMSHD
      CGROS=RCO2G/DMSHD
      ZADDBM=AMAX1(0.0,CGROSM*(CNSHX+CNLFM+CNLFX*CNPG))
      ZADDB=AMAX1(0.0,CGROS*(CNSHX+CNLFM+CNLFX*CNPG))
      PADDB=AMAX1(0.0,CGROS*(CPSHX+CPLFM+CPLFX*CNPG))
      CNRDM=AMAX1(0.0,1.70*ZADDBM)
      CNRDA=AMAX1(0.0,1.70*ZADDB)
C
C     TOTAL ABOVE-GROUND AUTOTROPHIC RESPIRATION BY BRANCH
C     ACCUMULATE GPP, SHOOT AUTOTROPHIC RESPIRATION, NET C EXCHANGE
C
C     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
C     RMNCS=maintenance respiration
C     RCO2GM,RCO2G=growth respiration limited by N,P unlimited,limited
C        by O2
C     SNCRM,SNCR=excess maintenance respiration unlimited,limited
C        by O2
C     CNRDM,CNRDA=respiration for N assimilation unlimited,limited
C        by O2
C     RCO2A=total branch respiration
C     RCO2M,RCO2N=RCO2A unlimited by O2,nonstructural C
C
      RCO2TM=RMNCS+RCO2GM+SNCRM+CNRDM
      RCO2T=RMNCS+RCO2G+SNCR+CNRDA
      RCO2M(1,NG(NZ,NY,NX),NZ,NY,NX)=RCO2M(1,NG(NZ,NY,NX),NZ,NY,NX)
     2+RCO2TM
      RCO2N(1,NG(NZ,NY,NX),NZ,NY,NX)=RCO2N(1,NG(NZ,NY,NX),NZ,NY,NX)
     2+RCO2T
      RCO2A(1,NG(NZ,NY,NX),NZ,NY,NX)=RCO2A(1,NG(NZ,NY,NX),NZ,NY,NX)
     2-RCO2T
      CH2O=0.0
C     IF(NZ.EQ.2)THEN
C     WRITE(*,4478)'RCO2E',I,J,NFZ,NX,NY,NZ,NB,IFLGC(NZ,NY,NX)
C    2,RCO2CM,VMXC,CPOOL(NB,NZ,NY,NX)
C    2,TFN4(NG(NZ,NY,NX),NZ,NY,NX),CNPG,WFNSG,FDBKX(NB,NZ,NY,NX)
C    3,RCO2C,WFR(1,NG(NZ,NY,NX),NZ,NY,NX)
4478  FORMAT(A8,8I4,20E12.4)
C     ENDIF
      ENDIF
C
C     REMOVE C,N,P USED IN MAINTENANCE + GROWTH REPIRATION AND GROWTH
C     FROM NON-STRUCTURAL POOLS
C
C     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
C     CH2O=total CH2O production
C     RMNCS=maintenance respiration
C     RCO2C=respiration from non-structural C
C     CGROS=total non-structural C used in growth and respiration
C     CNRDA=respiration for N assimilation
C     ZADDB,PADDB=nonstructural N,P used in growth
C     RNH3B=NH3 flux between atmosphere and branch from �uptake.f�
C     XFRC,XFRN,XFRP=branch-root layer C,N,P transfer
C
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+CH2O-AMIN1(RMNCS,RCO2C)
     2-CGROS-CNRDA
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-ZADDB+RNH3B(NB,NZ,NY,NX)
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-PADDB
C
C     TRANSFER OF C4 FIXATION PRODUCTS FROM NON-STRUCTURAL POOLS
C     IN MESOPHYLL TO THOSE IN BUNDLE SHEATH, DECARBOXYLATION
C     OF C4 FIXATION PRODUCTS IN BUNDLE SHEATH, LEAKAGE OF
C     DECARBOXYLATION PRODUCTS BACK TO MESOPHYLL IN C4 PLANTS
C
C     ICTYP=photosynthesis type:3=C3,4=C4
C
      IF(ICTYP(NZ,NY,NX).EQ.4)THEN
      DO 170 K=1,25
      IF(WGLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C
C     MESOPHYLL TO BUNDLE SHEATH TRANSFER
C
C     WGLF=node leaf C mass
C     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
C     CH2O3,CH2O4=total CO2 fixation in bundle sheath,mesophyll
C     CPL4M=mesophyll to bundle sheath transfer of nonstructural C4
C     FBS,FMP=leaf water content in bundle sheath, mesophyll
C
      CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)-CH2O3(K)
      CPOOL4(K,NB,NZ,NY,NX)=CPOOL4(K,NB,NZ,NY,NX)+CH2O4(K)
      CPL4M=1.0*(CPOOL4(K,NB,NZ,NY,NX)*WGLF(K,NB,NZ,NY,NX)*FBS
     2-CPOOL3(K,NB,NZ,NY,NX)*WGLF(K,NB,NZ,NY,NX)*FMP)
     2/(WGLF(K,NB,NZ,NY,NX)*(FBS+FMP))*XNFH
      CPOOL4(K,NB,NZ,NY,NX)=CPOOL4(K,NB,NZ,NY,NX)-CPL4M
      CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)+CPL4M
C
C     BUNDLE SHEATH CO2 DECARBOXYLATION
C
C     CCBS=CO2 concn in bundle sheath (uM)
C     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
C     WGLF=node leaf C mass
C     FBS,FMP=leaf water content in bundle sheath, mesophyll
C     CPL3K=bundle sheath CO2 decarboxylation
C     CO2KI=Ki for C3 leakage from bundle sheath to mesophyll in C4
C        (uM)
C     FCO2B,FHCOB=partition decarboxylation to CO2,HCO3
C     CPOOL3=C4 nonstructural C mass in bundle sheath
C
      CCBS=AMAX1(0.0,0.083E+09*CO2B(K,NB,NZ,NY,NX)
     2/(WGLF(K,NB,NZ,NY,NX)*FBS))
      CPL3K=2.5E-02*CPOOL3(K,NB,NZ,NY,NX)/(1.0+CCBS/CO2KI)*XNFH
      CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)-CPL3K
      CO2B(K,NB,NZ,NY,NX)=CO2B(K,NB,NZ,NY,NX)+FCO2B*CPL3K
      HCOB(K,NB,NZ,NY,NX)=HCOB(K,NB,NZ,NY,NX)+FHCOB*CPL3K
C
C     BUNDLE SHEATH LEAKAGE
C
C     CO2LK=bundle sheath CO2 leakage
C     CCBS=CO2 concn in bundle sheath (uM)
C     CO2L=intercellular CO2 concentration (uM)
C     WGLF=node leaf C mass
C     FBS=leaf water content in bundle sheath
C     FCO2B,FHCOB=partition decarboxylation to CO2,HCO3
C     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
C
      CO2LK=5.0E-07*(CCBS-CO2L(NZ,NY,NX))*WGLF(K,NB,NZ,NY,NX)*FBS*XNFH
      CO2B(K,NB,NZ,NY,NX)=CO2B(K,NB,NZ,NY,NX)-FCO2B*CO2LK
      HCOB(K,NB,NZ,NY,NX)=HCOB(K,NB,NZ,NY,NX)-FHCOB*CO2LK
C     IF(NB.EQ.1.AND.K.EQ.14)THEN
C     WRITE(*,6667)'CO2K',I,J,NB,K,CO2LK,CO2B(K,NB,NZ,NY,NX)
C    2,HCOB(K,NB,NZ,NY,NX),CPOOL3(K,NB,NZ,NY,NX),CH2O3(K),CH2O4(K)
C    3,CCBS,CO2L(NZ,NY,NX),WGLF(K,NB,NZ,NY,NX),CPL3Z
C    4,CPL3K,CPL4M,CPOOL4(K,NB,NZ,NY,NX),FBS,FMP
6667  FORMAT(A8,4I4,30E14.6)
C     ENDIF
C
C     TOTAL C EXCHANGE
C
C     TCO2T,TCO2A=total,above-ground PFT respiration
C     CNET=PFT net CO2 fixation
C     RECO=ecosystem respiration
C     TRAU=total autotrophic respiration
C     CO2LK=bundle sheath CO2 leakage
C
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)-CO2LK
      TCO2A(NZ,NY,NX)=TCO2A(NZ,NY,NX)-CO2LK
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-CO2LK
      HCNET(NZ,NY,NX)=HCNET(NZ,NY,NX)-CO2LK
      RECO(NY,NX)=RECO(NY,NX)-CO2LK
      TRAU(NY,NX)=TRAU(NY,NX)-CO2LK
      ENDIF
170   CONTINUE
      ENDIF
C
C     C,N,P GROWTH OF LEAF, SHEATH OR PETIOLE, STALK,
C     STALK RESERVES, REPRODUCTIVE ORGANS, GRAIN
C
C     GRO*,GRO*N,GRO*P=organ C,N,P growth rate
C     DM*=C production vs nonstructural C consumption
C     organ key:LF=leaf,SHE=petiole,STK=stalk,RSV=reserve
C        HSK=husk,EAR=ear,GR=grain,SHT=shoot
C     PART=organ partitioning fraction
C     CGROS=total non-structural C used in growth and growth
C        respiration
C     CN*,CP*=N:C,P:C ratios in plant organs
C     ZPLFM=minimum N:C,P:C in leaves relative to maximum values
C        from PFT file
C     ZPLFD=1.0-ZPLFM
C     CNPG=N,P constraint on growth respiration
C     WT*,WT*N,WT*P=organ C,N,P mass
C
      GROLF=PART(1)*CGROS*DMLFB
      GROSHE=PART(2)*CGROS*DMSHB
      GROSTK=PART(3)*CGROS*DMSTK(NZ,NY,NX)
      GRORSV=PART(4)*CGROS*DMRSV(NZ,NY,NX)
      GROHSK=PART(5)*CGROS*DMHSK(NZ,NY,NX)
      GROEAR=PART(6)*CGROS*DMEAR(NZ,NY,NX)
      GROGR=PART(7)*CGROS*DMGR(NZ,NY,NX)
      GROSHT=CGROS*DMSHT
      GROLFN=GROLF*CNLFB*(ZPLFM+ZPLFD*CNPG)
      GROSHN=GROSHE*CNSHB
      GROSTN=GROSTK*CNSTK(NZ,NY,NX)
      GRORSN=GRORSV*CNRSV(NZ,NY,NX)
      GROHSN=GROHSK*CNHSK(NZ,NY,NX)
      GROEAN=GROEAR*CNEAR(NZ,NY,NX)
      GROGRN=GROGR*CNRSV(NZ,NY,NX)
      GROLFP=GROLF*CPLFB*(ZPLFM+ZPLFD*CNPG)
      GROSHP=GROSHE*CPSHB
      GROSTP=GROSTK*CPSTK(NZ,NY,NX)
      GRORSP=GRORSV*CPRSV(NZ,NY,NX)
      GROHSP=GROHSK*CPHSK(NZ,NY,NX)
      GROEAP=GROEAR*CPEAR(NZ,NY,NX)
      GROGRP=GROGR*CPRSV(NZ,NY,NX)
      WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)+GROLF
      WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX)+GROSHE
      WTSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)+GROSTK
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+GRORSV
      WTHSKB(NB,NZ,NY,NX)=WTHSKB(NB,NZ,NY,NX)+GROHSK
      WTEARB(NB,NZ,NY,NX)=WTEARB(NB,NZ,NY,NX)+GROEAR
      WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)+GROLFN
      WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX)+GROSHN
      WTSTBN(NB,NZ,NY,NX)=WTSTBN(NB,NZ,NY,NX)+GROSTN
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+GRORSN
      WTHSBN(NB,NZ,NY,NX)=WTHSBN(NB,NZ,NY,NX)+GROHSN
      WTEABN(NB,NZ,NY,NX)=WTEABN(NB,NZ,NY,NX)+GROEAN
      WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)+GROLFP
      WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX)+GROSHP
      WTSTBP(NB,NZ,NY,NX)=WTSTBP(NB,NZ,NY,NX)+GROSTP
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+GRORSP
      WTHSBP(NB,NZ,NY,NX)=WTHSBP(NB,NZ,NY,NX)+GROHSP
      WTEABP(NB,NZ,NY,NX)=WTEABP(NB,NZ,NY,NX)+GROEAP
C
C     ETOLIATION
C
C     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concentration
C        in branch(g g-1)
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C     ETOL=coefficient for etoliation effects on expansion,extension
C
      CCE=AMIN1(CZPOLB(NB,NZ,NY,NX)/(CZPOLB(NB,NZ,NY,NX)
     2+CCPOLB(NB,NZ,NY,NX)*CNKI)
     3,CPPOLB(NB,NZ,NY,NX)/(CPPOLB(NB,NZ,NY,NX)
     4+CCPOLB(NB,NZ,NY,NX)*CPKI))
      ETOL=1.0+CCE
C
C     DISTRIBUTE LEAF GROWTH AMONG CURRENTLY GROWING NODES
C
C     MXNOD,MNNOD=maximum,minimum node number currently growing
C     KVSTG=integer of most recent leaf number
C     KNOD,GNOD=number of currently growing nodes
C     ALLOCL=fraction of leaf growth allocated to each node
C     GRO,GRON,GROP=leaf C,N,P growth at each node
C     NNOD=number of concurrently growing nodes
C
      IF(NB.EQ.NB1(NZ,NY,NX)
     2.AND.HTCTL(NZ,NY,NX).LE.SDPTH(NZ,NY,NX))THEN
      NNOD1=0
      ELSE
      NNOD1=1
      ENDIF
      IF(GROLF.GT.0.0)THEN
      MXNOD=KVSTG(NB,NZ,NY,NX)
      MNNOD=MAX(NNOD1,MXNOD-NNOD(NZ,NY,NX)+1)
      MXNOD=MAX(MXNOD,MNNOD)
      KNOD=MXNOD-MNNOD+1
      GNOD=KNOD
      ALLOCL=1.0/GNOD
      GRO=ALLOCL*GROLF
      GRON=ALLOCL*GROLFN
      GROP=ALLOCL*GROLFP
C
C     GROWTH AT EACH CURRENT NODE
C
C     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
C     GRO,GRON,GROP=leaf C,N,P growth at each node
C     CNWS,CPWS=protein:N,protein:P ratios from startq.f
C
      DO 490 KK=MNNOD,MXNOD
      K=MOD(KK,25)
      IF(K.EQ.0.AND.KK.NE.0)K=25
      WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX)+GRO
      WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)+GRON
      WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)+GROP
      WSLF(K,NB,NZ,NY,NX)=WSLF(K,NB,NZ,NY,NX)
     2+AMIN1(GRON*CNWS(NZ,NY,NX),GROP*CPWS(NZ,NY,NX))
C
C     SPECIFIC LEAF AREA FUNCTION OF CURRENT LEAF MASS
C     AT EACH NODE
C
C     SLA=specific area of leaf growth
C     ETOL=coefficient for etoliation effects on expansion,extension
C     SLA1=growth in leaf area vs mass from PFT file
C     SLA2=parameter for calculating leaf area expansion
C     WGLF=leaf C mass
C     PP=PFT population
C     WFNST=turgor expansion,extension function
C     GROA,GRO=leaf area,mass growth
C     ARLFB,ARLF=branch,node leaf area
C
      SLA=ETOL*SLA1(NZ,NY,NX)*(AMAX1(ZEROP2(NZ,NY,NX)
     2,WGLF(K,NB,NZ,NY,NX)*NBR(NZ,NY,NX))
     3/AMAX1(PPM*AREA(3,NU(NY,NX),NY,NX),PP(NZ,NY,NX)))
     3**SLA2(IBTYP(NZ,NY,NX))*WFNST
      GROA=GRO*SLA
      ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)+GROA
      ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX)+GROA
490   CONTINUE
      ENDIF
C
C     DISTRIBUTE SHEATH OR PETIOLE GROWTH AMONG
C     CURRENTLY GROWING NODES
C
C     MXNOD,MNNOD=max,min node number currently growing
C     KVSTG=integer of most recent leaf number
C     GNOD=number of currently growing nodes
C     ALLOCS=fraction of petiole growth allocated to each node
C     GRO,GRON,GROP=petiole C,N,P growth at each node
C     NNOD=number of concurrently growing nodes
C
      IF(GROSHE.GT.0.0)THEN
      MXNOD=KVSTG(NB,NZ,NY,NX)
      MNNOD=MAX(NNOD1,MXNOD-NNOD(NZ,NY,NX)+1)
      MXNOD=MAX(MXNOD,MNNOD)
      GNOD=MXNOD-MNNOD+1
      ALLOCS=1.0/GNOD
      GRO=ALLOCS*GROSHE
      GRON=ALLOCS*GROSHN
      GROP=ALLOCS*GROSHP
C
C     GROWTH AT EACH CURRENT NODE
C
C     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
C     GRO,GRON,GROP=petiole C,N,P growth at each node
C     CNWS,CPWS=protein:N,protein:P ratios from startq.f
C
      DO 505 KK=MNNOD,MXNOD
      K=MOD(KK,25)
      IF(K.EQ.0.AND.KK.NE.0)K=25
      WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX)+GRO
      WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX)+GRON
      WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX)+GROP
      WSSHE(K,NB,NZ,NY,NX)=WSSHE(K,NB,NZ,NY,NX)
     2+AMIN1(GRON*CNWS(NZ,NY,NX),GROP*CPWS(NZ,NY,NX))
C
C     SPECIFIC SHEATH OR PETIOLE LENGTH FUNCTION OF CURRENT MASS
C     AT EACH NODE
C
C     SSL=specific length of petiole growth
C     ETOL=coefficient for etoliation effects on expansion,extension
C     SSL1=growth in petiole length vs mass from PFT file
C     SSL2=parameter for calculating petiole extension
C     WGSHE=petiole C mass
C     PP=PFT population
C     WFNST=turgor expansion,extension function
C     GROS,GRO=petiole length,mass growth
C     HTSHE=petiole length
C
      IF(WGLF(K,NB,NZ,NY,NX).GT.0.0)THEN
      SSL=ETOL*SSL1(NZ,NY,NX)*(AMAX1(ZEROP2(NZ,NY,NX)
     2,WGSHE(K,NB,NZ,NY,NX)*NBR(NZ,NY,NX))
     3/AMAX1(PPM*AREA(3,NU(NY,NX),NY,NX),PP(NZ,NY,NX)))**SSL2*WFNST
      GROS=GRO/PP(NZ,NY,NX)*SSL
      HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX)+GROS*ANGSH(NZ,NY,NX)
C     IF(I.EQ.120.AND.J.EQ.24)THEN
C     WRITE(*,2526)'HTSHE',I,J,NZ,NB,K,SSL,WGSHE(K,NB,NZ,NY,NX)
C    2,HTSHE(K,NB,NZ,NY,NX),PP(NZ,NY,NX),SSL1(NZ,NY,NX)
C    3,SSL3,WFNST,GROS,GRO,ANGSH(NZ,NY,NX),ZEROP2(NZ,NY,NX)
C    4,CCPOLB(NB,NZ,NY,NX),ETOL
2526  FORMAT(A8,5I4,20E12.4)
C     ENDIF
      ENDIF
505   CONTINUE
      ENDIF
C
C     DISTRIBUTE STALK GROWTH AMONG CURRENTLY GROWING NODES
C
C     MXNOD,MNNOD=maximum,minimum node number currently growing
C     KVSTG=integer of most recent leaf number
C     GNOD=number of currently growing nodes
C     ALLOCN=fraction of stalk growth allocated to each node
C     GRO,GRON,GROP=stalk C,N,P growth at each node
C
      IF(IDAY(1,NB,NZ,NY,NX).EQ.0)THEN
      NN=0
      ELSE
      NN=1
      ENDIF
      MXNOD=KVSTG(NB,NZ,NY,NX)
      MNNOD=MAX(MIN(NN,MAX(NN,MXNOD-NNOD(NZ,NY,NX)))
     2,KVSTG(NB,NZ,NY,NX)-23)
      MXNOD=MAX(MXNOD,MNNOD)
      IF(GROSTK.GT.0.0)THEN
      GNOD=MXNOD-MNNOD+1
      ALLOCN=1.0/GNOD
      GRO=ALLOCN*GROSTK
      GRON=ALLOCN*GROSTN
      GROP=ALLOCN*GROSTP
C
C     SPECIFIC INTERNODE LENGTH FUNCTION OF CURRENT STALK MASS
C     AT EACH NODE
C
C     SNL=specific length of stalk growth
C     ETOL=coefficient for etoliation effects on expansion,extension
C     SNL1=growth in stalk length vs mass from PFT file
C     SNL2=parameter for calculating stalk extension
C     WTSKB=stalk C mass
C     PP=PFT population
C     GROH,GRO=stalk length,mass growth
C
      SNL=ETOL*SNL1(NZ,NY,NX)*(AMAX1(ZEROP2(NZ,NY,NX)
     2,WTSTKB(NB,NZ,NY,NX)*NBR(NZ,NY,NX))/PP(NZ,NY,NX))**SNL2
      GROH=GRO/PP(NZ,NY,NX)*SNL
      KX=MOD(MNNOD-1,25)
      IF(KX.EQ.0.AND.MNNOD-1.NE.0)KX=25
C
C     GROWTH AT EACH CURRENT NODE
C
C     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
C     GRO,GRON,GROP=stalk C,N,P growth at each node
C     HTNODX,HTNODE=stalk height,stalk internode length
C     ANGBR=sine of stalk angle from horizontal from PFT file
C
      DO 510 KK=MNNOD,MXNOD
      K1=MOD(KK,25)
      IF(K1.EQ.0.AND.KK.NE.0)K1=25
      K2=MOD(KK-1,25)
      IF(K2.EQ.0.AND.KK-1.NE.0)K2=25
      WGNODE(K1,NB,NZ,NY,NX)=WGNODE(K1,NB,NZ,NY,NX)+GRO
      WGNODN(K1,NB,NZ,NY,NX)=WGNODN(K1,NB,NZ,NY,NX)+GRON
      WGNODP(K1,NB,NZ,NY,NX)=WGNODP(K1,NB,NZ,NY,NX)+GROP
      HTNODX(K1,NB,NZ,NY,NX)=HTNODX(K1,NB,NZ,NY,NX)
     2+GROH*ANGBR(NZ,NY,NX)
      IF(K1.NE.0)THEN
      HTNODE(K1,NB,NZ,NY,NX)=HTNODX(K1,NB,NZ,NY,NX)
     2+HTNODE(K2,NB,NZ,NY,NX)
      ELSE
      HTNODE(K1,NB,NZ,NY,NX)=HTNODX(K1,NB,NZ,NY,NX)
      ENDIF
C     IF(NZ.EQ.1)THEN
C     WRITE(*,515)'HTNODE',I,J,NX,NY,NZ,NB,KK,K1,K2,MNNOD,MXNOD
C    1,NNOD(NZ,NY,NX),ARLF(K1,NB,NZ,NY,NX)
C    2,HTNODE(K1,NB,NZ,NY,NX),HTNODE(K2,NB,NZ,NY,NX),SNL,GRO
C    3,ALLOCN,WTSTKB(NB,NZ,NY,NX),WGNODE(K1,NB,NZ,NY,NX)
C    4,HTNODX(K1,NB,NZ,NY,NX),PP(NZ,NY,NX),GROSTK
515   FORMAT(A8,12I4,20E12.4)
C     ENDIF
510   CONTINUE
      ENDIF
C
C     RECOVERY OF REMOBILIZABLE N,P DURING REMOBILIZATION DEPENDS
C     ON SHOOT NON-STRUCTURAL C:N:P
C
C     IDAY(1,=emergence date
C     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concentration
C        in branch(g g-1)
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     RCCZ,RCCY=minimum,maximum fractions for shoot C recycling
C     RCCX,RCCQ=maximum fractions for shoot N,P recycling
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C
      IF(IDAY(1,NB,NZ,NY,NX).NE.0
     2.AND.CCPOLB(NB,NZ,NY,NX).GT.ZERO)THEN
      CCC=AMAX1(0.0,AMIN1(1.0
     1,CZPOLB(NB,NZ,NY,NX)/(CZPOLB(NB,NZ,NY,NX)
     2+CCPOLB(NB,NZ,NY,NX)*CNKI)
     3,CPPOLB(NB,NZ,NY,NX)/(CPPOLB(NB,NZ,NY,NX)
     4+CCPOLB(NB,NZ,NY,NX)*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLB(NB,NZ,NY,NX)/(CCPOLB(NB,NZ,NY,NX)
     2+CZPOLB(NB,NZ,NY,NX)/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLB(NB,NZ,NY,NX)/(CCPOLB(NB,NZ,NY,NX)
     2+CPPOLB(NB,NZ,NY,NX)/CPKI)))
      ELSE
      CCC=0.0
      CNC=0.0
      CPC=0.0
      ENDIF
      RCCC=RCCZ(IBTYP(NZ,NY,NX))+CCC*RCCY(IBTYP(NZ,NY,NX))
      RCCN=CNC*RCCX(IBTYP(NZ,NY,NX))
      RCCP=CPC*RCCQ(IBTYP(NZ,NY,NX))
C
C     WITHDRAW REMOBILIZABLE C,N,P FROM LOWEST NODE AFTER
C     MAXIMUM NODE NUMBER OF 25 IS REACHED
C
C     IFLGG=PFT senescence flag
C     KVSTG=integer of most recent leaf number
C     TFN3=temperature function for canopy growth
C     XRLA=rate of leaf appearance at 25 oC (h-1)
C     FSNC=fraction of lowest leaf to be remobilized
C
      IF(IFLGG(NB,NZ,NY,NX).EQ.1)THEN
      KVSTGX=MAX(0,KVSTG(NB,NZ,NY,NX)-24)
      IF(KVSTGX.GT.0)THEN
      K=MOD(KVSTGX,25)
      IF(K.EQ.0.AND.KVSTGX.GT.0)K=25
      KX=MOD(KVSTG(NB,NZ,NY,NX),25)
      IF(KX.EQ.0.AND.KVSTG(NB,NZ,NY,NX).NE.0)KX=25
      FSNC=TFN3(NZ,NY,NX)*XRLA(NZ,NY,NX)
C
C     REMOBILIZATION OF LEAF C,N,P ALSO DEPENDS ON STRUCTURAL C:N:P
C
C     IFLGP=flag for remobilization
C     WGLF,WGLFN,WGLFP=node leaf C,N,P mass
C     ARLF=node leaf area
C     RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
C
      IF(IFLGP(NB,NZ,NY,NX).EQ.1)THEN
      WGLFX(NB,NZ,NY,NX)=AMAX1(0.0,WGLF(K,NB,NZ,NY,NX))
      WGLFNX(NB,NZ,NY,NX)=AMAX1(0.0,WGLFN(K,NB,NZ,NY,NX))
      WGLFPX(NB,NZ,NY,NX)=AMAX1(0.0,WGLFP(K,NB,NZ,NY,NX))
      ARLFZ(NB,NZ,NY,NX)=AMAX1(0.0,ARLF(K,NB,NZ,NY,NX))
      IF(WGLFX(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RCCLX(NB,NZ,NY,NX)=RCCC*WGLFX(NB,NZ,NY,NX)
      RCZLX(NB,NZ,NY,NX)=WGLFNX(NB,NZ,NY,NX)*(RCCN+(1.0-RCCN)*RCCC)
      RCPLX(NB,NZ,NY,NX)=WGLFPX(NB,NZ,NY,NX)*(RCCP+(1.0-RCCP)*RCCC)
      ELSE
      RCCLX(NB,NZ,NY,NX)=0.0
      RCZLX(NB,NZ,NY,NX)=0.0
      RCPLX(NB,NZ,NY,NX)=0.0
      ENDIF
      ENDIF
C
C     FRACTION OF CURRENT LEAF TO BE REMOBILIZED
C
C     FSNC,FSNCL=fraction of lowest leaf to be remobilized
C     WGLF,WGLFX=node living,senescing leaf C mass
C
      IF(FSNC*WGLFX(NB,NZ,NY,NX).GT.WGLF(K,NB,NZ,NY,NX)
     2.AND.WGLFX(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FSNCL=AMAX1(0.0,WGLF(K,NB,NZ,NY,NX)/WGLFX(NB,NZ,NY,NX))
      ELSE
      FSNCL=FSNC
      ENDIF
      FSNAL=0.0
C
C     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
C     TO FRACTIONS SET IN 'STARTQ'
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FSNCL=fraction of lowest leaf to be remobilized
C     RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
C     WGLFX,WGLFNX,WGLFPX=senescing leaf C,N,P mass
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
C
      DO 6300 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*FSNCL*WGLFX(NB,NZ,NY,NX)*FWODB(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*FSNCL*WGLFNX(NB,NZ,NY,NX)*FWODLN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*FSNCL*WGLFPX(NB,NZ,NY,NX)*FWODLP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(1,M,NZ,NY,NX)
     2*FSNCL*(WGLFX(NB,NZ,NY,NX)-RCCLX(NB,NZ,NY,NX))*FWODB(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(1,M,NZ,NY,NX)
     2*FSNCL*(WGLFNX(NB,NZ,NY,NX)-RCZLX(NB,NZ,NY,NX))*FWODLN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(1,M,NZ,NY,NX)
     2*FSNCL*(WGLFPX(NB,NZ,NY,NX)-RCPLX(NB,NZ,NY,NX))*FWODLP(1,NZ)
6300  CONTINUE
C
C     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
C
C     FSNCL=fraction of lowest leaf to be remobilized
C     ARLFB,ARLFZ=branch living,senescing leaf area
C     WTLFB,WTLFBN,WTLFBP,WGLFX,WGLFNX,WGLFPX=C,N,P mass in
C        living,senescing leaf
C     WSLF=leaf protein mass
C     CNWS,CPWS=protein:N,protein:P ratios from startq.f
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
C     RCCLX,RCZLX,RCPLX=remobilization of C,N,P from senescing leaf
C
      ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)
     2-FSNAL*ARLFZ(NB,NZ,NY,NX)
      WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)
     2-FSNCL*WGLFX(NB,NZ,NY,NX)
      WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)
     2-FSNCL*WGLFNX(NB,NZ,NY,NX)
      WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)
     2-FSNCL*WGLFPX(NB,NZ,NY,NX)
      ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX)
     2-FSNAL*ARLFZ(NB,NZ,NY,NX)
      WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX)
     2-FSNCL*WGLFX(NB,NZ,NY,NX)
      WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)
     2-FSNCL*WGLFNX(NB,NZ,NY,NX)
      WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)
     2-FSNCL*WGLFPX(NB,NZ,NY,NX)
      WSLF(K,NB,NZ,NY,NX)=AMAX1(0.0,WSLF(K,NB,NZ,NY,NX)
     2-FSNCL*AMAX1(WGLFNX(NB,NZ,NY,NX)*CNWS(NZ,NY,NX)
     3,WGLFPX(NB,NZ,NY,NX)*CPWS(NZ,NY,NX)))
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)
     2+FSNCL*RCCLX(NB,NZ,NY,NX)*FWODB(1,NZ)
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)
     2+FSNCL*RCZLX(NB,NZ,NY,NX)*FWODLN(1,NZ)
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)
     2+FSNCL*RCPLX(NB,NZ,NY,NX)*FWODLP(1,NZ)
C
C     REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P ALSO DEPENDS ON
C     STRUCTURAL C:N:P
C
C     IFLGP=flag for remobilization
C     WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
C     HTSHEX=petiole length
C     RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
C
      IF(IFLGP(NB,NZ,NY,NX).EQ.1)THEN
      WGSHEX(NB,NZ,NY,NX)=AMAX1(0.0,WGSHE(K,NB,NZ,NY,NX))
      WGSHNX(NB,NZ,NY,NX)=AMAX1(0.0,WGSHN(K,NB,NZ,NY,NX))
      WGSHPX(NB,NZ,NY,NX)=AMAX1(0.0,WGSHP(K,NB,NZ,NY,NX))
      HTSHEX(NB,NZ,NY,NX)=AMAX1(0.0,HTSHE(K,NB,NZ,NY,NX))
      IF(WGSHEX(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RCCSX(NB,NZ,NY,NX)=RCCC*WGSHEX(NB,NZ,NY,NX)
      RCZSX(NB,NZ,NY,NX)=WGSHNX(NB,NZ,NY,NX)
     2*(RCCN+(1.0-RCCN)*RCCSX(NB,NZ,NY,NX)/WGSHEX(NB,NZ,NY,NX))
      RCPSX(NB,NZ,NY,NX)=WGSHPX(NB,NZ,NY,NX)
     2*(RCCP+(1.0-RCCP)*RCCSX(NB,NZ,NY,NX)/WGSHEX(NB,NZ,NY,NX))
      ELSE
      RCCSX(NB,NZ,NY,NX)=0.0
      RCZSX(NB,NZ,NY,NX)=0.0
      RCPSX(NB,NZ,NY,NX)=0.0
      ENDIF
      WTSTXB(NB,NZ,NY,NX)=WTSTXB(NB,NZ,NY,NX)+WGNODE(K,NB,NZ,NY,NX)
      WTSTXN(NB,NZ,NY,NX)=WTSTXN(NB,NZ,NY,NX)+WGNODN(K,NB,NZ,NY,NX)
      WTSTXP(NB,NZ,NY,NX)=WTSTXP(NB,NZ,NY,NX)+WGNODP(K,NB,NZ,NY,NX)
C     IF(NZ.EQ.2)THEN
C     WRITE(*,2358)'WTSTXB',I,J,NZ,NB,K,WTSTXB(NB,NZ,NY,NX)
C    2,WTSTKB(NB,NZ,NY,NX),WGNODE(K,NB,NZ,NY,NX)
2358  FORMAT(A8,5I4,12E12.4)
C     ENDIF
      WGNODE(K,NB,NZ,NY,NX)=0.0
      WGNODN(K,NB,NZ,NY,NX)=0.0
      WGNODP(K,NB,NZ,NY,NX)=0.0
      HTNODX(K,NB,NZ,NY,NX)=0.0
      ENDIF
C
C     FRACTION OF CURRENT SHEATH TO BE REMOBILIZED
C
C     FSNCS=fraction of lowest petiole to be remobilized
C     WGSHE,WGSHEX=node living,senescing sheath C mass
C
      IF(FSNC*WGSHEX(NB,NZ,NY,NX).GT.WGSHE(K,NB,NZ,NY,NX)
     2.AND.WGSHEX(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FSNCS=AMAX1(0.0,WGSHE(K,NB,NZ,NY,NX)/WGSHEX(NB,NZ,NY,NX))
      ELSE
      FSNCS=FSNC
      ENDIF
C
C     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
C     TO FRACTIONS SET IN 'STARTQ'
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FSNCS=fraction of lowest petiole to be remobilized
C     RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
C     WGSHX,WGSHNX,WGSHPX=senescing petiole C,N,P mass
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
C
      DO 6305 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*FSNCS*WGSHEX(NB,NZ,NY,NX)*FWODB(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*FSNCS*WGSHNX(NB,NZ,NY,NX)*FWODSN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*FSNCS*WGSHPX(NB,NZ,NY,NX)*FWODSP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(2,M,NZ,NY,NX)
     2*FSNCS*(WGSHEX(NB,NZ,NY,NX)-RCCSX(NB,NZ,NY,NX))*FWODB(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(2,M,NZ,NY,NX)
     2*FSNCS*(WGSHNX(NB,NZ,NY,NX)-RCZSX(NB,NZ,NY,NX))*FWODSN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(2,M,NZ,NY,NX)
     2*FSNCS*(WGSHPX(NB,NZ,NY,NX)-RCPSX(NB,NZ,NY,NX))*FWODSP(1,NZ)
6305  CONTINUE
C
C     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
C
C     FSNCS=fraction of lowest petiole to be remobilized
C     HTSHE,HTSHEX=living,senescing petiole length
C     WTSHB,WTSHBN,WTSHBP,WGSHEX,WGSHNX,WGSHPX=C,N,P mass in
C        living,senescing petiole
C     WSSHE=petiole protein mass
C     CNWS,CPWS=protein:N,protein:P ratios from startq.f
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
C     RCCSX,RCZSX,RCPSX=remobilization of C,N,P from senescing petiole
C
      WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX)
     2-FSNCS*WGSHEX(NB,NZ,NY,NX)
      WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX)
     2-FSNCS*WGSHNX(NB,NZ,NY,NX)
      WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX)
     2-FSNCS*WGSHPX(NB,NZ,NY,NX)
      HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX)
     2-FSNCS*HTSHEX(NB,NZ,NY,NX)
      WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX)
     2-FSNCS*WGSHEX(NB,NZ,NY,NX)
      WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX)
     2-FSNCS*WGSHNX(NB,NZ,NY,NX)
      WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX)
     2-FSNCS*WGSHPX(NB,NZ,NY,NX)
      WSSHE(K,NB,NZ,NY,NX)=AMAX1(0.0,WSSHE(K,NB,NZ,NY,NX)
     2-FSNCS*AMAX1(WGSHNX(NB,NZ,NY,NX)*CNWS(NZ,NY,NX)
     3,WGSHPX(NB,NZ,NY,NX)*CPWS(NZ,NY,NX)))
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)
     2+FSNCS*RCCSX(NB,NZ,NY,NX)*FWODB(1,NZ)
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)
     2+FSNCS*RCZSX(NB,NZ,NY,NX)*FWODSN(1,NZ)
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)
     2+FSNCS*RCPSX(NB,NZ,NY,NX)*FWODSP(1,NZ)
      ENDIF
      ENDIF
C
C     REMOBILIZATION OF STALK RESERVE C,N,P IF GROWTH RESPIRATION < 0
C
C     IFLGZ=remobilization flag
C     SNCR=excess maintenance respiration
C     WTRSVB=stalk reserve C mass
C     RCO2V=remobilization of stalk reserve C
C     VMXC=rate constant for nonstructural C oxidation in respiration
C     TFN3=temperature function for canopy growth
C
      IF(IFLGZ.EQ.0)THEN
      IF(SNCR.GT.0.0.AND.WTRSVB(NB,NZ,NY,NX).GT.0.0)THEN
      RCO2V=AMIN1(SNCR,VMXC*WTRSVB(NB,NZ,NY,NX)*TFN3(NZ,NY,NX))*XNFH
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)-RCO2V
      SNCR=SNCR-RCO2V
      ENDIF
      ENDIF
C
C     TOTAL REMOBILIZATION = GROWTH RESPIRATION < 0 + DECIDUOUS LEAF
C     FALL DURING AUTUMN + REMOBILZATION DURING GRAIN FILL IN ANNUALS
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IFLGY,IFLGZ=remobilization flags
C     SNCZ=phenologically-driven respiration senescence during
C        late-season
C     FXFB=rate constant for plant-storage nonstructural C,N,P
C        exchange
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     WTLSB=leaf+petiole mass
C     FLGZ=control rate of remobilization
C     FLGZX=number of hours until full senescence after physl maturity
C     SNCX=total senescence respiration
C     KVSTG,KVSTGN=integer of highest,lowest leaf number currently
C        growing
C     KSNC=number of nodes undergoing remobilization
C     SNCF=ratio of phenologically-driven vs total senescence
C        respiration
C
      IF(IFLGZ.EQ.1.AND.IFLGY.EQ.1.AND.ISTYP(NZ,NY,NX).NE.0)THEN
      SNCZ=FXFB(IBTYP(NZ,NY,NX))*WTLSB(NB,NZ,NY,NX)
     2*AMIN1(1.0,FLGZ(NB,NZ,NY,NX)/FLGZX)*XNFH
      ELSE
      SNCZ=0.0
      ENDIF
      SNCX=SNCR+SNCZ
      IF(SNCX.GT.ZEROP(NZ,NY,NX))THEN
      SNCF=SNCZ/SNCX
      KSNC=INT(0.5*(KVSTG(NB,NZ,NY,NX)-KVSTGN(NB,NZ,NY,NX)))+1
      XKSNC=KSNC
      KN=MAX(0,KVSTGN(NB,NZ,NY,NX)-1)
C     IF(NZ.EQ.2.OR.NZ.EQ.3)THEN
C     WRITE(*,1266)'SNCX0',I,J,NX,NY,NZ,NB,SNCY,SNCR,SNCX,SNCF
C    2,CPOOL(NB,NZ,NY,NX),WTLSB(NB,NZ,NY,NX),RCCC
1266  FORMAT(A8,6I4,12E16.8)
C     ENDIF
C
C     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH STALK RESERVES
C     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
C
C     WVSTKB=stalk sapwood mass
C     WTSTKB=branch stalk mass
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     XFRC,XFRN,XFRC=C,N,P transfer among branches
C     FXFY,FXFZ=rate constant for leaf-reserve nonstructural C,N,P
C        exchange (h-1)

      IF(NB.EQ.NB1(NZ,NY,NX))THEN
      DO 580 NBL=1,NBR(NZ,NY,NX)
      IF(NBL.NE.NB1(NZ,NY,NX))THEN
      IF(WVSTKB(NBL,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WVSTKT=WVSTKB(NB1(NZ,NY,NX),NZ,NY,NX)+WVSTKB(NBL,NZ,NY,NX)
      WTRSVT=WTRSVB(NB1(NZ,NY,NX),NZ,NY,NX)+WTRSVB(NBL,NZ,NY,NX)
      WTRSVD=(WTRSVB(NB1(NZ,NY,NX),NZ,NY,NX)*WVSTKB(NBL,NZ,NY,NX)
     2-WTRSVB(NBL,NZ,NY,NX)*WVSTKB(NB1(NZ,NY,NX),NZ,NY,NX))/WVSTKT
      XFRC=FXFY(ISTYP(NZ,NY,NX))*WTRSVD
      WTRSVB(NB1(NZ,NY,NX),NZ,NY,NX)=WTRSVB(NB1(NZ,NY,NX),NZ,NY,NX)
     2-XFRC
      WTRSVB(NBL,NZ,NY,NX)=WTRSVB(NBL,NZ,NY,NX)+XFRC
      IF(WTRSVT.GT.ZEROP(NZ,NY,NX))THEN
      WTRSND=(WTRSBN(NB1(NZ,NY,NX),NZ,NY,NX)*WTRSVB(NBL,NZ,NY,NX)
     2-WTRSBN(NBL,NZ,NY,NX)*WTRSVB(NB1(NZ,NY,NX),NZ,NY,NX))/WTRSVT
      WTRSPD=(WTRSBP(NB1(NZ,NY,NX),NZ,NY,NX)*WTRSVB(NBL,NZ,NY,NX)
     2-WTRSBP(NBL,NZ,NY,NX)*WTRSVB(NB1(NZ,NY,NX),NZ,NY,NX))/WTRSVT
      XFRN=FXFZ(ISTYP(NZ,NY,NX))*WTRSND
      XFRP=FXFZ(ISTYP(NZ,NY,NX))*WTRSPD
      WTRSBN(NB1(NZ,NY,NX),NZ,NY,NX)=WTRSBN(NB1(NZ,NY,NX),NZ,NY,NX)
     2-XFRN
      WTRSBN(NBL,NZ,NY,NX)=WTRSBN(NBL,NZ,NY,NX)+XFRN
      WTRSBP(NB1(NZ,NY,NX),NZ,NY,NX)=WTRSBP(NB1(NZ,NY,NX),NZ,NY,NX)
     2-XFRP
      WTRSBP(NBL,NZ,NY,NX)=WTRSBP(NBL,NZ,NY,NX)+XFRP
C     WRITE(*,4487)'BRANCH',I,J,NFZ,NX,NY,NZ,NB1(NZ,NY,NX),NBL
C    2,XFRC,XFRN,XFRP
C    2,WVSTKB(NB1(NZ,NY,NX),NZ,NY,NX),WTRSVB(NB1(NZ,NY,NX),NZ,NY,NX)
C    3,WTRSBN(NB1(NZ,NY,NX),NZ,NY,NX),WTRSBP(NB1(NZ,NY,NX),NZ,NY,NX)
C    2,WVSTKB(NBL,NZ,NY,NX),WTRSVB(NBL,NZ,NY,NX)
C    3,WTRSBN(NBL,NZ,NY,NX),WTRSBP(NBL,NZ,NY,NX)
4487  FORMAT(A8,8I4,20E12.4)
      ENDIF
      ENDIF
      ENDIF
580   CONTINUE
      ENDIF
C
C     REMOBILIZATION AND LITTERFALL WHEN GROWTH RESPIRATION < 0
C     STARTING FROM LOWEST LEAFED NODE AND PROCEEDING UPWARDS
C
C     SNCX,SNCT=branch,node senescence respiration
C     KSNC=number of nodes undergoing remobilization
C
C     IF(NZ.EQ.2.OR.NZ.EQ.3)THEN
C     WRITE(*,1266)'SNCX1',I,J,NX,NY,NZ,NB,SNCY,SNCR,SNCX,SNCF
C    2,CPOOL(NB,NZ,NY,NX),WTLSB(NB,NZ,NY,NX),RCCC
C     ENDIF
      DO 575 N=1,KSNC
      SNCT=SNCX/XKSNC
      DO 650 KK=KN,KVSTG(NB,NZ,NY,NX)
      SNCLF=0.0
      SNCSH=0.0
      K=MOD(KK,25)
      IF(K.EQ.0.AND.KK.NE.0)K=25
C
C     REMOBILIZATION OF LEAF C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
C
C     WGLF,WGSHE=node leaf,petiole C mass
C     SCNF,SCNSH=leaf,petiole senescence respiration
C     RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     RCCZ,RCCY=minimum,maximum fractions for shoot C recycling
C
      IF(WGLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FNCLF=WGLF(K,NB,NZ,NY,NX)/(WGLF(K,NB,NZ,NY,NX)
     2+WGSHE(K,NB,NZ,NY,NX))
      SNCLF=FNCLF*SNCT
      SNCSH=SNCT-SNCLF
      RCCL=RCCC*WGLF(K,NB,NZ,NY,NX)
      RCZL=WGLFN(K,NB,NZ,NY,NX)*(RCCN+(1.0-RCCN)*RCCC)
      RCPL=WGLFP(K,NB,NZ,NY,NX)*(RCCP+(1.0-RCCP)*RCCC)
C
C     FRACTION OF CURRENT LEAF TO BE REMOBILIZED
C
C     FSNCL,FSNAL=fraction of current leaf C,area to be remobilized
C     SCNLF=leaf senescence respiration
C     RCCL=remobilization of C from senescing leaf
C
      IF(RCCL.GT.ZEROP(NZ,NY,NX))THEN
      FSNCL=AMAX1(0.0,AMIN1(1.0,SNCLF/RCCL))
      ELSE
      FSNCL=1.0
      ENDIF
      FSNAL=0.0
C
C     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
C     TO FRACTIONS SET IN 'STARTQ'
C
C     CSNC,ZSNC,PSNC=literfall C,N,P
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FSNCL=fraction of current leaf to be remobilized
C     WGLF,WGLFN,WGLFP=node leaf C,N,P mass
C     RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
C     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
C
      DO 6310 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*FSNCL*WGLF(K,NB,NZ,NY,NX)*FWODB(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*FSNCL*WGLFN(K,NB,NZ,NY,NX)*FWODLN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*FSNCL*WGLFP(K,NB,NZ,NY,NX)*FWODLP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(1,M,NZ,NY,NX)
     2*FSNCL*(WGLF(K,NB,NZ,NY,NX)-RCCL)*FWODB(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(1,M,NZ,NY,NX)
     2*FSNCL*(WGLFN(K,NB,NZ,NY,NX)-RCZL)*FWODLN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(1,M,NZ,NY,NX)
     2*FSNCL*(WGLFP(K,NB,NZ,NY,NX)-RCPL)*FWODLP(1,NZ)
6310  CONTINUE
      IF(K.NE.0)THEN
      CSNC(2,1,0,NZ,NY,NX)=CSNC(2,1,0,NZ,NY,NX)
     2+FSNCL*(CPOOL3(K,NB,NZ,NY,NX)+CPOOL4(K,NB,NZ,NY,NX))
      CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)
     2-FSNCL*CPOOL3(K,NB,NZ,NY,NX)
      CPOOL4(K,NB,NZ,NY,NX)=CPOOL4(K,NB,NZ,NY,NX)
     2-FSNCL*CPOOL4(K,NB,NZ,NY,NX)
      ENDIF
C
C     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
C
C     ARLFB=leaf area
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
C     FSNCL=fraction of current leaf to be remobilized
C     CNWS,CPWS=protein:N,protein:P ratios from startq.f
C
      ARLFB(NB,NZ,NY,NX)=AMAX1(0.0,ARLFB(NB,NZ,NY,NX)
     2-FSNAL*ARLF(K,NB,NZ,NY,NX))
      WTLFB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX)
     2-FSNCL*WGLF(K,NB,NZ,NY,NX))
      WTLFBN(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBN(NB,NZ,NY,NX)
     2-FSNCL*WGLFN(K,NB,NZ,NY,NX))
      WTLFBP(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBP(NB,NZ,NY,NX)
     2-FSNCL*WGLFP(K,NB,NZ,NY,NX))
      ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX)
     2-FSNAL*ARLF(K,NB,NZ,NY,NX)
      WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX)
     2-FSNCL*WGLF(K,NB,NZ,NY,NX)
      WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)
     2-FSNCL*WGLFN(K,NB,NZ,NY,NX)
      WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)
     2-FSNCL*WGLFP(K,NB,NZ,NY,NX)
      WSLF(K,NB,NZ,NY,NX)=AMAX1(0.0,WSLF(K,NB,NZ,NY,NX)
     2-FSNCL*AMAX1(WGLFN(K,NB,NZ,NY,NX)*CNWS(NZ,NY,NX)
     3,WGLFP(K,NB,NZ,NY,NX)*CPWS(NZ,NY,NX)))
C
C     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
C     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
C
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
C     FSNCL=fraction of current leaf C to be remobilized
C     RCCL,RCZL,RCPL=remobilization of C,N,P from senescing leaf
C     SNCLF,SNCT=remaining senescence respiration carried to next node
C     WTLFB=leaf C mass
C
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+FSNCL*RCCL*SNCF
     2*FWODB(1,NZ)
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+FSNCL*RCZL*FWODLN(1,NZ)
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+FSNCL*RCPL*FWODLP(1,NZ)
      SNCLF=SNCLF-FSNCL*RCCL*FWODB(1,NZ)
      SNCT=SNCT-FSNCL*RCCL*FWODB(1,NZ)
      IF(WTLFB(NB,NZ,NY,NX).LE.ZEROP2(NZ,NY,NX))THEN
      WTLFB(NB,NZ,NY,NX)=0.0
      ARLFB(NB,NZ,NY,NX)=0.0
      ENDIF
C
C     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
C
      IF(SNCLF.LE.ZEROP(NZ,NY,NX))GO TO 564
C
C     OTHERWISE REMAINING C,N,P IN LEAF GOES TO LITTERFALL
C
C     CSNC,ZSNC,PSNC=literfall C,N,P
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
C     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
C     ARLFB=leaf area
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
C
      ELSE
      DO 6315 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*WGLF(K,NB,NZ,NY,NX)*FWODB(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*WGLFN(K,NB,NZ,NY,NX)*FWODLN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*WGLFP(K,NB,NZ,NY,NX)*FWODLP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(1,M,NZ,NY,NX)
     2*WGLF(K,NB,NZ,NY,NX)*FWODB(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(1,M,NZ,NY,NX)
     2*WGLFN(K,NB,NZ,NY,NX)*FWODLN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(1,M,NZ,NY,NX)
     2*WGLFP(K,NB,NZ,NY,NX)*FWODLP(1,NZ)
6315  CONTINUE
      IF(K.NE.0)THEN
      CSNC(2,1,0,NZ,NY,NX)=CSNC(2,1,0,NZ,NY,NX)
     2+CPOOL3(K,NB,NZ,NY,NX)+CPOOL4(K,NB,NZ,NY,NX)
      CPOOL3(K,NB,NZ,NY,NX)=0.0
      CPOOL4(K,NB,NZ,NY,NX)=0.0
      ENDIF
      ARLFB(NB,NZ,NY,NX)=AMAX1(0.0,ARLFB(NB,NZ,NY,NX)
     2-ARLF(K,NB,NZ,NY,NX))
      WTLFB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX)
     2-WGLF(K,NB,NZ,NY,NX))
      WTLFBN(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBN(NB,NZ,NY,NX)
     2-WGLFN(K,NB,NZ,NY,NX))
      WTLFBP(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBP(NB,NZ,NY,NX)
     2-WGLFP(K,NB,NZ,NY,NX))
      ARLF(K,NB,NZ,NY,NX)=0.0
      WGLF(K,NB,NZ,NY,NX)=0.0
      WGLFN(K,NB,NZ,NY,NX)=0.0
      WGLFP(K,NB,NZ,NY,NX)=0.0
      WSLF(K,NB,NZ,NY,NX)=0.0
      IF(WTLFB(NB,NZ,NY,NX).LE.ZEROP2(NZ,NY,NX))THEN
      WTLFB(NB,NZ,NY,NX)=0.0
      ARLFB(NB,NZ,NY,NX)=0.0
      ENDIF
      ENDIF
564   CONTINUE
C
C     REMOBILIZATION OF SHEATHS OR PETIOLE C,N,P DEPENDS ON
C     NON-STRUCTURAL C:N:P
C
C     WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
C     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C
      IF(WGSHE(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RCCS=RCCC*WGSHE(K,NB,NZ,NY,NX)
      RCZS=WGSHN(K,NB,NZ,NY,NX)*(RCCN+(1.0-RCCN)*RCCC)
      RCPS=WGSHP(K,NB,NZ,NY,NX)*(RCCP+(1.0-RCCP)*RCCC)
C
C     FRACTION OF REMOBILIZATION THAT CAN BE MET FROM CURRENT SHEATH
C     OR PETIOLE
C
C     FSNCS,FSNAS=fraction of current petiole C,length to be
C        remobilized
C     SCNSH=petiole senescence respiration
C     RCCS=remobilization of C from senescing petiole
C
      IF(RCCS.GT.ZEROP(NZ,NY,NX))THEN
      FSNCS=AMAX1(0.0,AMIN1(1.0,SNCSH/RCCS))
      ELSE
      FSNCS=1.0
      ENDIF
      FSNAS=FSNCS
C
C     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
C     TO FRACTIONS SET IN 'STARTQ'
C
C     CSNC,ZSNC,PSNC=literfall C,N,P
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FSNCS=fraction of current petiole to be remobilized
C     WGSHE,WGSHN,WGSHP=node petiole C,N,P mass
C     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
C
      DO 6320 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*FSNCS*WGSHE(K,NB,NZ,NY,NX)*FWODB(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*FSNCS*WGSHN(K,NB,NZ,NY,NX)*FWODSN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*FSNCS*WGSHP(K,NB,NZ,NY,NX)*FWODSP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(2,M,NZ,NY,NX)
     2*FSNCS*(WGSHE(K,NB,NZ,NY,NX)-RCCS)*FWODB(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(2,M,NZ,NY,NX)
     2*FSNCS*(WGSHN(K,NB,NZ,NY,NX)-RCZS)*FWODSN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(2,M,NZ,NY,NX)
     2*FSNCS*(WGSHP(K,NB,NZ,NY,NX)-RCPS)*FWODSP(1,NZ)
6320  CONTINUE
C
C     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
C
C     HTSHE=petiole length
C     WTSHEB,WTLFBN,WTSHBP=branch petiole C,N,P mass
C     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
C     FSNCS=fraction of current petiole to be remobilized
C     CNWS,CPWS=protein:N,protein:P ratios from startq.f
C
      WTSHEB(NB,NZ,NY,NX)=AMAX1(0.0,WTSHEB(NB,NZ,NY,NX)
     2-FSNCS*WGSHE(K,NB,NZ,NY,NX))
      WTSHBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSHBN(NB,NZ,NY,NX)
     2-FSNCS*WGSHN(K,NB,NZ,NY,NX))
      WTSHBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSHBP(NB,NZ,NY,NX)
     2-FSNCS*WGSHP(K,NB,NZ,NY,NX))
      HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX)
     2-FSNAS*HTSHE(K,NB,NZ,NY,NX)
      WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX)
     2-FSNCS*WGSHE(K,NB,NZ,NY,NX)
      WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX)
     2-FSNCS*WGSHN(K,NB,NZ,NY,NX)
      WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX)
     2-FSNCS*WGSHP(K,NB,NZ,NY,NX)
      WSSHE(K,NB,NZ,NY,NX)=AMAX1(0.0,WSSHE(K,NB,NZ,NY,NX)
     2-FSNCS*AMAX1(WGSHN(K,NB,NZ,NY,NX)*CNWS(NZ,NY,NX)
     3,WGSHP(K,NB,NZ,NY,NX)*CPWS(NZ,NY,NX)))
C
C     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
C     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
C
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
C     FSNCS=fraction of current petiole C to be remobilized
C     RCCS,RCZS,RCPS=remobilization of C,N,P from senescing petiole
C     SNCSH,SNCT=remaining senescence respiration carried to next node
C
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+FSNCS*RCCS*SNCF
     2*FWODB(1,NZ)
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+FSNCS*RCZS*FWODSN(1,NZ)
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+FSNCS*RCPS*FWODSN(1,NZ)
      SNCSH=SNCSH-FSNCS*RCCS*FWODB(1,NZ)
      SNCT=SNCT-FSNCS*RCCS*FWODB(1,NZ)
      IF(WTSHEB(NB,NZ,NY,NX).LE.ZEROP2(NZ,NY,NX))THEN
      WTSHEB(NB,NZ,NY,NX)=0.0
      ENDIF
C
C     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
C
C     SNCSH=remaining senescence respiration carried to next node
C
      IF(SNCSH.LE.ZEROP(NZ,NY,NX))GO TO 565
C
C     OTHERWISE REMAINING C,N,P IN SHEATH OR PETIOLE GOES TO
C     LITTERFALL
C
C     CSNC,ZSNC,PSNC=literfall C,N,P
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FWODB=C woody fraction in branch:0=woody,1=non-woody
C     FWODSN,FWODSP=N,P woody fraction in petiole:0=woody,1=non-woody
C     HTSHE=petiole length
C     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
C     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
C
      ELSE
      DO 6325 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*WGSHE(K,NB,NZ,NY,NX)*FWODB(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*WGSHN(K,NB,NZ,NY,NX)*FWODSN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*WGSHP(K,NB,NZ,NY,NX)*FWODSP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(2,M,NZ,NY,NX)
     2*WGSHE(K,NB,NZ,NY,NX)*FWODB(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(2,M,NZ,NY,NX)
     2*WGSHN(K,NB,NZ,NY,NX)*FWODSN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(2,M,NZ,NY,NX)
     2*WGSHP(K,NB,NZ,NY,NX)*FWODSP(1,NZ)
6325  CONTINUE
      WTSHEB(NB,NZ,NY,NX)=AMAX1(0.0,WTSHEB(NB,NZ,NY,NX)
     2-WGSHE(K,NB,NZ,NY,NX))
      WTSHBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSHBN(NB,NZ,NY,NX)
     2-WGSHN(K,NB,NZ,NY,NX))
      WTSHBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSHBP(NB,NZ,NY,NX)
     2-WGSHP(K,NB,NZ,NY,NX))
      HTSHE(K,NB,NZ,NY,NX)=0.0
      WGSHE(K,NB,NZ,NY,NX)=0.0
      WGSHN(K,NB,NZ,NY,NX)=0.0
      WGSHP(K,NB,NZ,NY,NX)=0.0
      WSSHE(K,NB,NZ,NY,NX)=0.0
      IF(WTSHEB(NB,NZ,NY,NX).LE.ZEROP2(NZ,NY,NX))THEN
      WTSHEB(NB,NZ,NY,NX)=0.0
      ENDIF
      ENDIF
650   CONTINUE
      KN=KN+1
      SNCR=SNCT*(1.0-SNCF)
C
C     REMOBILIZATION OF RESERVE C
C
C     WTRSVB=stalk reserve C mass
C     SNCR=excess maintenance respiration
C
      IF(WTRSVB(NB,NZ,NY,NX).GT.SNCR)THEN
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)-SNCR
      SNCR=0.0
      GO TO 565
      ENDIF
C
C     REMOBILIZATION OF STALK C,N,P
C
C     SNCZ=phenologically-driven respiration senescence during
C        late-season
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     WTSTKB,WVSTKB=stalk,sapwood C mass
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     MXNOD,MNNOD=maximum,minimum node number currently growing
C     NNOD=number of concurrently growing nodes
C     KVSTG=integer of most recent leaf number
C
      SNCZ=0.0
      SNCT=SNCR+SNCZ
      IF(ISTYP(NZ,NY,NX).NE.0.AND.SNCT.GT.ZEROP(NZ,NY,NX)
     2.AND.WTSTKB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      SNCF=SNCZ/SNCT
      FRCC=WVSTKB(NB,NZ,NY,NX)/WTSTKB(NB,NZ,NY,NX)
      RCSC=RCCC*FRCC
      RCSN=RCCN*FRCC
      RCSP=RCCP*FRCC
      MXNOD=KVSTG(NB,NZ,NY,NX)
      MNNOD=MAX(MIN(0,MAX(0,MXNOD-NNOD(NZ,NY,NX)))
     2,KVSTG(NB,NZ,NY,NX)-23)
      MXNOD=MAX(MXNOD,MNNOD)
      DO 1650 KK=MXNOD,MNNOD,-1
      K=MAX(0,MOD(KK,25))
      IF(K.EQ.0.AND.KK.NE.0)K=25
C
C     REMOBILIZATION OF STALK C,N,P DEPENDS ON NON-STRUCTURAL C:N:P
C
C     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
C     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
C
      IF(WGNODE(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RCCK=RCSC*WGNODE(K,NB,NZ,NY,NX)
      RCZK=WGNODN(K,NB,NZ,NY,NX)*(RCSN+(1.0-RCSN)*RCSC)
      RCPK=WGNODP(K,NB,NZ,NY,NX)*(RCSP+(1.0-RCSP)*RCSC)
C
C     FRACTION OF CURRENT NODE TO BE REMOBILIZED
C
C     FSNCS=fraction of lowest internode to be remobilized
C     SNCT=remaining node senescence respiration
C     RCCK=remobilization of C from senescing internode
C
      IF(RCCK.GT.ZEROP(NZ,NY,NX))THEN
      FSNCK=AMAX1(0.0,AMIN1(1.0,SNCT/RCCK))
      ELSE
      FSNCK=1.0
      ENDIF
C
C     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
C     TO FRACTIONS SET IN 'STARTQ'
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FSNCK=fraction of lowest internode to be remobilized
C     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
C     WGNODE,WGNODN,WGNODP=senescing internode C,N,P mass
C
      DO 7310 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*FSNCK*WGNODE(K,NB,NZ,NY,NX)*FWOOD(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*FSNCK*WGNODN(K,NB,NZ,NY,NX)*FWOODN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*FSNCK*WGNODP(K,NB,NZ,NY,NX)*FWOODP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)
     2*FSNCK*(WGNODE(K,NB,NZ,NY,NX)-RCCK)*FWOOD(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)
     2*FSNCK*(WGNODN(K,NB,NZ,NY,NX)-RCZK)*FWOODN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)
     2*FSNCK*(WGNODP(K,NB,NZ,NY,NX)-RCPK)*FWOODP(1,NZ)
7310  CONTINUE
C
C     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
C
C     FSNCK=fraction of lowest internode to be remobilized
C     HTNODE,HTNODX=living,senescing internode length
C     WTSTKB,WTSTBN,WTSTBP,WGNODE,WGNODN,WGNODP=C,N,P mass in
C        senescing internode
C
      WTSTKB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTKB(NB,NZ,NY,NX)
     2-FSNCK*WGNODE(K,NB,NZ,NY,NX))
      WTSTBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBN(NB,NZ,NY,NX)
     2-FSNCK*WGNODN(K,NB,NZ,NY,NX))
      WTSTBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBP(NB,NZ,NY,NX)
     2-FSNCK*WGNODP(K,NB,NZ,NY,NX))
      HTNODE(K,NB,NZ,NY,NX)=HTNODE(K,NB,NZ,NY,NX)
     2-FSNCK*HTNODX(K,NB,NZ,NY,NX)
      WGNODE(K,NB,NZ,NY,NX)=WGNODE(K,NB,NZ,NY,NX)
     2-FSNCK*WGNODE(K,NB,NZ,NY,NX)
      WGNODN(K,NB,NZ,NY,NX)=WGNODN(K,NB,NZ,NY,NX)
     2-FSNCK*WGNODN(K,NB,NZ,NY,NX)
      WGNODP(K,NB,NZ,NY,NX)=WGNODP(K,NB,NZ,NY,NX)
     2-FSNCK*WGNODP(K,NB,NZ,NY,NX)
      HTNODX(K,NB,NZ,NY,NX)=HTNODX(K,NB,NZ,NY,NX)
     2-FSNCK*HTNODX(K,NB,NZ,NY,NX)
C     IF(NZ.EQ.1.AND.ICHKF.EQ.1)THEN
C     WRITE(*,2356)'WGNODE1',I,J,NFZ,NX,NY,NZ,NB,K,KK,MXNOD,MNNOD
C    2,KSNC,RCCC,FRCC,RCSC,SNCT,WGNODE(K,NB,NZ,NY,NX)
C    3,HTNODX(K,NB,NZ,NY,NX),WTSTKB(NB,NZ,NY,NX)
C    4,CPOOL(NB,NZ,NY,NX)
C     ENDIF
C
C     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
C     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
C
C     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     FSNCK=fraction of lowest internode to be remobilized
C     SNCT=remaining node senescence respiration
C
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+FSNCK*RCCK*SNCF
     2*FWOOD(1,NZ)
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+FSNCK*RCZK*FWOODN(1,NZ)
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+FSNCK*RCPK*FWOODP(1,NZ)
      SNCT=SNCT-FSNCK*RCCK*FWOOD(1,NZ)
C
C     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
C
      IF(SNCT.LE.ZEROP(NZ,NY,NX))GO TO 565
C
C     OTHERWISE REMAINING C,N,P IN NODE GOES TO LITTERFALL
C
C     CSNC,ZSNC,PSNC=literfall C,N,P
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTSTKB,WTSTBN,WTSTBP,WGNODE,WGNODN,WGNODP=C,N,P mass in
C        senescing internode
C     HTNODE,HTNODX=living,senescing internode length
C
      ELSE
      DO 7315 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*WGNODE(K,NB,NZ,NY,NX)*FWOOD(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*WGNODN(K,NB,NZ,NY,NX)*FWOODN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*WGNODP(K,NB,NZ,NY,NX)*FWOODP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)
     2*WGNODE(K,NB,NZ,NY,NX)*FWOOD(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)
     2*WGNODN(K,NB,NZ,NY,NX)*FWOODN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)
     2*WGNODP(K,NB,NZ,NY,NX)*FWOODP(1,NZ)
7315  CONTINUE
      WTSTKB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTKB(NB,NZ,NY,NX)
     2-WGNODE(K,NB,NZ,NY,NX))
      WTSTBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBN(NB,NZ,NY,NX)
     2-WGNODN(K,NB,NZ,NY,NX))
      WTSTBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBP(NB,NZ,NY,NX)
     2-WGNODP(K,NB,NZ,NY,NX))
      HTNODE(K,NB,NZ,NY,NX)=HTNODE(K,NB,NZ,NY,NX)
     2-HTNODX(K,NB,NZ,NY,NX)
      WGNODE(K,NB,NZ,NY,NX)=0.0
      WGNODN(K,NB,NZ,NY,NX)=0.0
      WGNODP(K,NB,NZ,NY,NX)=0.0
      HTNODX(K,NB,NZ,NY,NX)=0.0
C     IF(NZ.EQ.1.AND.ICHKF.EQ.1)THEN
C     WRITE(*,2356)'WGNODE2',I,J,NFZ,NX,NY,NZ,NB,K,KK,MXNOD,MNNOD
C    2,KSNC,RCCC,FRCC,RCSC,SNCT,WGNODE(K,NB,NZ,NY,NX)
C    3,HTNODX(K,NB,NZ,NY,NX),WTSTKB(NB,NZ,NY,NX)
C    4,CPOOL(NB,NZ,NY,NX)
2356  FORMAT(A8,12I4,12E16.8)
C     ENDIF
      ENDIF
1650  CONTINUE
C
C     RESIDUAL STALK
C
C     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
C     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
C
      IF(WTSTXB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RCCK=RCSC*WTSTXB(NB,NZ,NY,NX)
      RCZK=WTSTXN(NB,NZ,NY,NX)*(RCSN+(1.0-RCSN)*RCSC)
      RCPK=WTSTXP(NB,NZ,NY,NX)*(RCSP+(1.0-RCSP)*RCSC)
C
C     FRACTION OF RESIDUAL STALK TO BE REMOBILIZED
C
C     FSNCR=fraction of residual stalk to be remobilized
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
C
      IF(RCCK.GT.ZEROP(NZ,NY,NX))THEN
      FSNCR=AMAX1(0.0,AMIN1(1.0,SNCT/RCCK))
      ELSE
      FSNCR=1.0
      ENDIF
C
C     NON-REMOBILIZABLE C,N,P BECOMES LITTERFALL ALLOCATED
C     TO FRACTIONS SET IN 'STARTQ'
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FSNCR=fraction of residual stalk to be remobilized
C     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
C     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,
C        1=non-woody
C
      DO 8310 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*FSNCR*WTSTXB(NB,NZ,NY,NX)*FWOOD(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*FSNCR*WTSTXN(NB,NZ,NY,NX)*FWOODN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*FSNCR*WTSTXP(NB,NZ,NY,NX)*FWOODP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)
     2*FSNCR*(WTSTXB(NB,NZ,NY,NX)-RCCK)*FWOOD(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)
     2*FSNCR*(WTSTXN(NB,NZ,NY,NX)-RCZK)*FWOODN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)
     2*FSNCR*(WTSTXP(NB,NZ,NY,NX)-RCPK)*FWOODP(1,NZ)
8310  CONTINUE
C
C     UPDATE STATE VARIABLES FOR REMOBILIZATION AND LITTERFALL
C
C     FSNCR=fraction of residual stalk to be remobilized
C     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in senescing stalk
C     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
C     HTNODE,HTNODX=living,senescing internode length
C
      WTSTKB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTKB(NB,NZ,NY,NX)
     2-FSNCR*WTSTXB(NB,NZ,NY,NX))
      WTSTBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBN(NB,NZ,NY,NX)
     2-FSNCR*WTSTXN(NB,NZ,NY,NX))
      WTSTBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBP(NB,NZ,NY,NX)
     2-FSNCR*WTSTXP(NB,NZ,NY,NX))
      WTSTXB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTXB(NB,NZ,NY,NX)
     2-FSNCR*WTSTXB(NB,NZ,NY,NX))
      WTSTXN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTXN(NB,NZ,NY,NX)
     2-FSNCR*WTSTXN(NB,NZ,NY,NX))
      WTSTXP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTXP(NB,NZ,NY,NX)
     2-FSNCR*WTSTXP(NB,NZ,NY,NX))
      HTNODZ=0.0
      DO 8320 K=0,25
      HTNODZ=AMAX1(HTNODZ,HTNODE(K,NB,NZ,NY,NX))
8320  CONTINUE
      HTNODZ=AMAX1(0.0,HTNODZ-FSNCR*HTNODZ)
      DO 8325 K=0,25
      HTNODE(K,NB,NZ,NY,NX)=AMIN1(HTNODZ,HTNODE(K,NB,NZ,NY,NX))
8325  CONTINUE
C
C     FRACTION OF C REMOBILIZED FOR GROWTH RESPIRATION < 0 IS
C     RESPIRED AND NOT TRANSFERRED TO NON-STRUCTURAL POOLS
C
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     RCCK,RCZK,RCPK=remobilization of C,N,P from senescing internode
C     FSNCR=fraction of residual stalk to be remobilized
C     SNCT=remaining node senescence respiration
C
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+FSNCR*RCCK*SNCF
     2*FWOOD(1,NZ)
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+FSNCR*RCZK*FWOODN(1,NZ)
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+FSNCR*RCPK*FWOODP(1,NZ)
      SNCT=SNCT-FSNCR*RCCK*FWOOD(1,NZ)
C     IF(NZ.EQ.1.AND.ICHKF.EQ.1)THEN
C     WRITE(*,2357)'WTSTXB1',I,J,NFZ,NX,NY,NZ,NB,K,FSNCR,SNCT
C    3,WTSTKB(NB,NZ,NY,NX),WTSTXB(NB,NZ,NY,NX)
C    4,(HTNODE(K,NB,NZ,NY,NX),K=0,25)
2357  FORMAT(A8,8I4,40E12.4)
C     ENDIF
      ENDIF
C
C     EXIT LOOP IF REMOBILIZATION REQUIREMENT HAS BEEN MET
C
C     SNCT=remaining node senescence respiration
C
      IF(SNCT.LE.ZEROP(NZ,NY,NX))GO TO 565
C
C     OTHERWISE REMAINING C,N,P IN NODE GOES TO LITTERFALL
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in senescing stalk
C
      ELSE
      IF(ICHKF.EQ.0)THEN
      DO 8315 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*WTSTXB(NB,NZ,NY,NX)*FWOOD(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*WTSTXN(NB,NZ,NY,NX)*FWOODN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*WTSTXP(NB,NZ,NY,NX)*FWOODP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)
     2*WTSTXB(NB,NZ,NY,NX)*FWOOD(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)
     2*WTSTXN(NB,NZ,NY,NX)*FWOODN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)
     2*WTSTXP(NB,NZ,NY,NX)*FWOODP(1,NZ)
8315  CONTINUE
C     IF(NZ.EQ.1.AND.ICHKF.EQ.1)THEN
C     WRITE(*,2357)'WTSTXB2',I,J,NFZ,NX,NY,NZ,NB,K,FSNCR,SNCT
C    3,WTSTKB(NB,NZ,NY,NX),WTSTXB(NB,NZ,NY,NX)
C    4,(HTNODE(K,NB,NZ,NY,NX),K=0,25)
C     ENDIF
      WTSTKB(NB,NZ,NY,NX)=AMAX1(0.0,WTSTKB(NB,NZ,NY,NX)
     2-WTSTXB(NB,NZ,NY,NX))
      WTSTBN(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBN(NB,NZ,NY,NX)
     2-WTSTXN(NB,NZ,NY,NX))
      WTSTBP(NB,NZ,NY,NX)=AMAX1(0.0,WTSTBP(NB,NZ,NY,NX)
     2-WTSTXP(NB,NZ,NY,NX))
      WTSTXB(NB,NZ,NY,NX)=0.0
      WTSTXN(NB,NZ,NY,NX)=0.0
      WTSTXP(NB,NZ,NY,NX)=0.0
      ENDIF
      ENDIF
C
C     REMOBILIZATION OF STORAGE C,N,P
C
C     WTRVC=storage C
C     IDTHB=branch living flag: 0=alive,1=dead
C     SNCR=remaining excess maintenance respiration
C
      SNCR=SNCT
      IF(WTRVC(NZ,NY,NX).GT.SNCR)THEN
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)-SNCR
      SNCR=0.0
      GO TO 565
      ELSE
      SNCR=SNCR-WTRVC(NZ,NY,NX)
      WTRVC(NZ,NY,NX)=0.0
      WTSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)-SNCR
      IF(ISTYP(NZ,NY,NX).NE.0.OR.IWTYP(NZ,NY,NX).EQ.0)THEN
      IDTHB(NB,NZ,NY,NX)=1
      ELSE
      VRNF(NB,NZ,NY,NX)=VRNX(NB,NZ,NY,NX)+0.5
      ENDIF
      WRITE(*,1113)'DEATHC',IYRC,I,J,NFZ,NX,NY,NZ,NB
     2,IFLGF(NB,NZ,NY,NX),IDTHB(NB,NZ,NY,NX),IDTHP(NZ,NY,NX)
     3,IDTHR(NZ,NY,NX),SNCT,SNCR,WTRVC(NZ,NY,NX),VRNF(NB,NZ,NY,NX)
1113  FORMAT(A8,12I4,12E12.4)
      ENDIF
565   CONTINUE
575   CONTINUE
      ENDIF
595   CONTINUE
C
C     DEATH IF MAIN STALK OF TREE DIES
C
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     IDTHB=branch living flag: 0=alive,1=dead
C     KVSTGX,KVSTG=integer of lowest,highest leaf number currently
C        growing
C     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
C     ZPLFM=min N:C,P:C in leaves relative to max values from PFT file
C     CNLFB,CPLFB=N:C,P:C ratios in leaf
C
      IF(IBTYP(NZ,NY,NX).NE.0.AND.IGTYP(NZ,NY,NX).GT.1
     2.AND.IDTHB(NB1(NZ,NY,NX),NZ,NY,NX).EQ.1)IDTHB(NB,NZ,NY,NX)=1
C
C     REMOBILIZE EXCESS LEAF STRUCTURAL N,P
C
      KVSTGX=MAX(0,KVSTG(NB,NZ,NY,NX)-24)
      DO 495 KK=KVSTGX,KVSTG(NB,NZ,NY,NX)
      K=MAX(0,MOD(KK,25))
      IF(K.EQ.0.AND.KK.NE.0)K=25
      IF(WGLF(K,NB,NZ,NY,NX).GT.0.0)THEN
      CPOOLT=WGLF(K,NB,NZ,NY,NX)+CPOOL(NB,NZ,NY,NX)
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLD=WGLFN(K,NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX)
     2-ZPOOL(NB,NZ,NY,NX)*WGLF(K,NB,NZ,NY,NX)
      XFRN1=AMAX1(0.0,AMIN1(1.0E-03*ZPOOLD/CPOOLT,WGLFN(K,NB,NZ,NY,NX)
     2-ZPLFM*CNLFB*WGLF(K,NB,NZ,NY,NX)))
      PPOOLD=WGLFP(K,NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX)
     2-PPOOL(NB,NZ,NY,NX)*WGLF(K,NB,NZ,NY,NX)
      XFRP1=AMAX1(0.0,AMIN1(1.0E-03*PPOOLD/CPOOLT,WGLFP(K,NB,NZ,NY,NX)
     2-ZPLFM*CPLFB*WGLF(K,NB,NZ,NY,NX)))
      XFRN=AMAX1(XFRN1,10.0*XFRP1)
      XFRP=AMAX1(XFRP1,0.10*XFRN1)
      WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)-XFRN
      WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)-XFRN
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+XFRN
      WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)-XFRP
      WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)-XFRP
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+XFRP
      WSLF(K,NB,NZ,NY,NX)=AMAX1(0.0,WSLF(K,NB,NZ,NY,NX)
     2-AMAX1(XFRN*CNWS(NZ,NY,NX),XFRP*CPWS(NZ,NY,NX)))
      ENDIF
      ENDIF
495   CONTINUE
C
C     ALLOCATION OF LEAF AREA TO CANOPY LAYERS
C
C     HTCTL=hypocotyledon height
C     SDPTH=seeding depth
C     ARLF=node leaf area
C     HTSHE=petiole length
C
      KVSTGN(NB,NZ,NY,NX)=0
      IF(HTCTL(NZ,NY,NX).LE.SDPTH(NZ,NY,NX)
     2.AND.ARLF(0,NB1(NZ,NY,NX),NZ,NY,NX).GT.0.0)THEN
      XLGLF=SQRT(1.0E+02*ARLF(0,NB1(NZ,NY,NX),NZ,NY,NX)
     2/PP(NZ,NY,NX))
      HTCTL(NZ,NY,NX)=XLGLF+HTSHE(0,NB1(NZ,NY,NX),NZ,NY,NX)
     2+HTNODE(0,NB1(NZ,NY,NX),NZ,NY,NX)
      ENDIF
C
C     IF CANOPY HAS EMERGED INITIALIZE LEAF VARIABLES
C
      IF(HTCTL(NZ,NY,NX).GT.SDPTH(NZ,NY,NX))THEN
      DO 540 K=0,25
      DO 540 L=1,JC
      ARLFL(L,K,NB,NZ,NY,NX)=0.0
      WGLFL(L,K,NB,NZ,NY,NX)=0.0
      WGLFLN(L,K,NB,NZ,NY,NX)=0.0
      WGLFLP(L,K,NB,NZ,NY,NX)=0.0
540   CONTINUE
      DO 535 L=1,JC
      ARSTK(L,NB,NZ,NY,NX)=0.0
535   CONTINUE
C
C     BRANCH HEIGHT
C
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     KVSTG1,KVSTGN=integer of highest,lowest leaf number currently
C        growing
C     HTNODE=internode length
C     HTBR=branch base height
C
      IF(IBTYP(NZ,NY,NX).NE.0.AND.IGTYP(NZ,NY,NX).GT.1)THEN
      IF(NB.NE.NB1(NZ,NY,NX))THEN
      KVSTG1=MAX(1,KVSTG(NB1(NZ,NY,NX),NZ,NY,NX)-24)
      IF(NBTB(NB,NZ,NY,NX).GE.KVSTG1)THEN
      K=MAX(0,MOD(NBTB(NB,NZ,NY,NX),25))
      IF(K.EQ.0.AND.KK.NE.0)K=25
      HTBR=HTNODE(K,NB1(NZ,NY,NX),NZ,NY,NX)
      ELSE
      HTBR=0.0
      ENDIF
      ELSE
      HTBR=0.0
      ENDIF
      ELSE
      HTBR=0.0
      ENDIF
      KVSTGX=MAX(0,KVSTG(NB,NZ,NY,NX)-24)
C
C     FOR ALL LEAFED NODES
C
      DO 560 KK=KVSTGX,KVSTG(NB,NZ,NY,NX)
      K=MAX(0,MOD(KK,25))
      IF(K.EQ.0.AND.KK.NE.0)K=25
C
C     HEIGHT OF STALK INTERNODE + SHEATH OR PETIOLE
C     AND LENGTH OF LEAF
C
C     HTSTK=stalk height
C     HTNODE=internode length
C     HTLF=leaf node height
C     ARLF=leaf node area
C     PP=plant population
C     XLGLF=leaf length
C
      HTSTK=HTBR+HTNODE(K,NB,NZ,NY,NX)
      HTLF=HTSTK+HTSHE(K,NB,NZ,NY,NX)
      XLGLF=AMAX1(0.0,SQRT(WDLF(NZ,NY,NX)*AMAX1(0.0
     2,ARLF(K,NB,NZ,NY,NX))/PP(NZ,NY,NX)))
      TLGLF=0.0
C
C     ALLOCATE FRACTIONS OF LEAF IN EACH INCLINATION CLASS
C     FROM HIGHEST TO LOWEST TO CANOPY LAYER
C
C     YLGLF=leaf elevation
C     CLASS=leaf inclination class
C     XLGLF=leaf length
C     ZC,ZCX=canopy height
C     HTLFL,HTLFU=height of leaf base,tip
C     ZL=height to bottom of each canopy layer
C     LHTLFL,LHTLFU=layer number of leaf base,tip
C     FRACL=leaf fraction in each layer
C
      DO 555 N=4,1,-1
      YLGLF=ZSIN(N)*CLASS(N,NZ,NY,NX)*XLGLF
      HTLFL=AMIN1(ZCX(NZ,NY,NX)+0.01-YLGLF,HTLF+TLGLF)
      HTLFU=AMIN1(ZCX(NZ,NY,NX)+0.01,HTLFL+YLGLF)
      LU=0
      LL=0
      DO 550 L=JC,1,-1
      IF(LU.EQ.1.AND.LL.EQ.1)GO TO 551
      IF((HTLFU.GT.ZL(L-1,NY,NX).OR.ZL(L-1,NY,NX).LE.ZERO)
     2.AND.LU.EQ.0)THEN
      LHTLFU=MAX(1,L)
      LU=1
      ENDIF
      IF((HTLFL.GT.ZL(L-1,NY,NX).OR.ZL(L-1,NY,NX).LE.ZERO)
     2.AND.LL.EQ.0)THEN
      LHTLFL=MAX(1,L)
      LL=1
      ENDIF
550   CONTINUE
551   CONTINUE
      DO 570 L=LHTLFL,LHTLFU
      IF(LHTLFU.EQ.LHTLFL)THEN
      FRACL=CLASS(N,NZ,NY,NX)
      ELSEIF(HTLFU.GT.HTLFL.AND.ZL(L,NY,NX).GT.HTLFL)THEN
      FRACL=CLASS(N,NZ,NY,NX)*(AMIN1(HTLFU,ZL(L,NY,NX))
     2-AMAX1(HTLFL,ZL(L-1,NY,NX)))/(HTLFU-HTLFL)
      ELSE
      FRACL=CLASS(N,NZ,NY,NX)
      ENDIF
      YARLF=FRACL*ARLF(K,NB,NZ,NY,NX)
      YWGLF=FRACL*WGLF(K,NB,NZ,NY,NX)
      YWGLFN=FRACL*WGLFN(K,NB,NZ,NY,NX)
      YWGLFP=FRACL*WGLFP(K,NB,NZ,NY,NX)
C
C     ACCUMULATE LAYER LEAF AREAS, C, N AND P CONTENTS
C
C     ARLFL=node leaf area in canopy layer
C     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
C     ARLFV,WGLFV=total leaf area,C in canopy layer
C     HTNODE=internode length
C
      ARLFL(L,K,NB,NZ,NY,NX)=ARLFL(L,K,NB,NZ,NY,NX)+YARLF
      WGLFL(L,K,NB,NZ,NY,NX)=WGLFL(L,K,NB,NZ,NY,NX)+YWGLF
      WGLFLN(L,K,NB,NZ,NY,NX)=WGLFLN(L,K,NB,NZ,NY,NX)+YWGLFN
      WGLFLP(L,K,NB,NZ,NY,NX)=WGLFLP(L,K,NB,NZ,NY,NX)+YWGLFP
      ARLFV(L,NZ,NY,NX)=ARLFV(L,NZ,NY,NX)+YARLF
      WGLFV(L,NZ,NY,NX)=WGLFV(L,NZ,NY,NX)+YWGLF
C     IF(NZ.EQ.2)THEN
C     WRITE(*,4813)'GRO',I,J,NX,NY,NZ,NB,K,KK,L,LHTLFL,LHTLFU
C    2,FRACL,HTLFU,HTLFL,ZL(L-1,NY,NX),ARLFB(NB,NZ,NY,NX)
C    3,ARLF(K,NB,NZ,NY,NX),WTLFB(NB,NZ,NY,NX),WGLF(K,NB,NZ,NY,NX)
C    4,ARLFP(NZ,NY,NX),ZL(L,NY,NX),HTLF,TLGLF,HTSTK,HTBR
C    4,HTNODE(K,NB,NZ,NY,NX),HTSHE(K,NB,NZ,NY,NX),YLGLF
C    5,ZSIN(N),CLASS(N,NZ,NY,NX),XLGLF,ZC(NZ,NY,NX)
C    6,ZCX(NZ,NY,NX)
4813  FORMAT(A8,11I4,30E12.4)
C     ENDIF
570   CONTINUE
      TLGLF=TLGLF+YLGLF
      ZC(NZ,NY,NX)=AMAX1(ZC(NZ,NY,NX),HTLFU)
555   CONTINUE
      IF(WSSHE(K,NB,NZ,NY,NX).GT.0.0)THEN
      IF(KVSTGN(NB,NZ,NY,NX).EQ.0)KVSTGN(NB,NZ,NY,NX)
     2=MIN(KK,KVSTG(NB,NZ,NY,NX))
      ENDIF
560   CONTINUE
      IF(KVSTGN(NB,NZ,NY,NX).EQ.0)KVSTGN(NB,NZ,NY,NX)
     2=KVSTG(NB,NZ,NY,NX)
      K1=MAX(0,MOD(KVSTG(NB,NZ,NY,NX),25))
      IF(K1.EQ.0.AND.KVSTG(NB,NZ,NY,NX).NE.0)K1=25
      K2=MAX(0,MOD(KVSTG(NB,NZ,NY,NX)-1,25))
      IF(K2.EQ.0.AND.KVSTG(NB,NZ,NY,NX)-1.NE.0)K2=25
      IF(HTNODE(K1,NB,NZ,NY,NX).EQ.0.0)THEN
      HTNODE(K1,NB,NZ,NY,NX)=HTNODE(K2,NB,NZ,NY,NX)
      ENDIF
      HTLFB=HTBR
     2+AMAX1(0.0,HTNODE(K1,NB,NZ,NY,NX))
C
C     ALLOCATE STALK SURFACE AREA TO CANOPY LAYERS
C
C     HTNODE=internode length
C     HTLFB=leaf base height
C     ZL=height to bottom of each canopy layer
C     LHTBRL,LHTBRU=layer number of branch base,tip
C     WTSTKB,ARSTKB=branch stalk mass,surface area
C     WVSTKB=stalk sapwood mass
C     FRACL=stalk fraction in each layer
C     ARSTK=total branch stalk surface area in each layer
C
C     IF(NZ.EQ.1)THEN
C     WRITE(*,6679)'K1',I,J,NZ,NB,K1,KVSTG(NB,NZ,NY,NX)
C    2,HTNODE(K1,NB,NZ,NY,NX)
6679  FORMAT(A8,6I4,12E12.4)
C     ENDIF
      IF(HTNODE(K1,NB,NZ,NY,NX).GT.0.0)THEN
      LU=0
      LL=0
      DO 545 L=JC,1,-1
      IF(LU.EQ.1.AND.LL.EQ.1)GO TO 546
      IF((HTLFB.GT.ZL(L-1,NY,NX).OR.ZL(L-1,NY,NX).LE.ZERO)
     2.AND.LU.EQ.0)THEN
      LHTBRU=MAX(1,L)
      LU=1
      ENDIF
      IF((HTBR.GT.ZL(L-1,NY,NX).OR.ZL(L-1,NY,NX)
     2.LE.ZERO).AND.LL.EQ.0)THEN
      LHTBRL=MAX(1,L)
      LL=1
      ENDIF
545   CONTINUE
546   CONTINUE
      RSTK=SQRT(VSTK*(AMAX1(0.0,WTSTKB(NB,NZ,NY,NX))/PP(NZ,NY,NX))
     3/(3.1416*HTNODE(K1,NB,NZ,NY,NX)))
      ARSTKB(NB)=6.2832*RSTK*HTNODE(K1,NB,NZ,NY,NX)*PP(NZ,NY,NX)
      IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      WVSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)
      ELSE
      ZSTK=AMIN1(YSTK,0.1*RSTK)
      ASTV=3.1416*(2.0*RSTK*ZSTK-ZSTK**2)
      WVSTKB(NB,NZ,NY,NX)=ASTV/VSTK
     2*HTNODE(K1,NB,NZ,NY,NX)*PP(NZ,NY,NX)
      ENDIF
C     IF((I/210)*210.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.1.AND.NZ.EQ.1)THEN
C     WRITE(*,6677)'WVSTK',I,J,NX,NY,NZ,NB,K1,WVSTKB(NB,NZ,NY,NX)
C    2,ASTV,VSTK,ZSTK,RSTK,HTNODE(K1,NB,NZ,NY,NX),PP(NZ,NY,NX)
C    3,WTSTKB(NB,NZ,NY,NX),WTLSB(NB,NZ,NY,NX)
6677  FORMAT(A8,7I4,12E12.4)
C     ENDIF
      DO 445 L=LHTBRL,LHTBRU
      IF(HTLFB.GT.HTBR)THEN
      IF(HTLFB.GT.ZL(L-1,NY,NX))THEN
      FRACL=(AMIN1(HTLFB,ZL(L,NY,NX))-AMAX1(HTBR
     2,ZL(L-1,NY,NX)))/(HTLFB-HTBR)
      ELSE
      FRACL=0.0
      ENDIF
      ELSE
      FRACL=1.0
      ENDIF
      ARSTK(L,NB,NZ,NY,NX)=FRACL*ARSTKB(NB)
445   CONTINUE
      ELSE
      IF(SNL1(NZ,NY,NX).GT.0.0)THEN
      WVSTKB(NB,NZ,NY,NX)=0.0
      ELSE
      WVSTKB(NB,NZ,NY,NX)=WTLSB(NB,NZ,NY,NX)
      ENDIF
      DO 450 L=1,JC
      ARSTK(L,NB,NZ,NY,NX)=0.0
450   CONTINUE
      ENDIF
      ELSE
      IF(SNL1(NZ,NY,NX).GT.0.0)THEN
      WVSTKB(NB,NZ,NY,NX)=0.0
      ELSE
      WVSTKB(NB,NZ,NY,NX)=WTLSB(NB,NZ,NY,NX)
      ENDIF
      DO 455 L=1,JC
      ARSTK(L,NB,NZ,NY,NX)=0.0
455   CONTINUE
      ENDIF
C
C     ALLOCATE LEAF AREA TO INCLINATION CLASSES ACCORDING TO
C     DISTRIBUTION ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
C
C     SSIN=sine of solar angle
C     SURF=leaf node surface area in canopy layer
C     ARLF,ARLFL=leaf node surface area in canopy layer
C     ZC,DPTHS=canopy,snowpack height
C     CLASS=leaf inclination class
C
      IF(SSINN(NY,NX).GT.0.0)THEN
      DO 900 K=1,25
      DO 900 L=1,JC
      DO 900 N=1,4
      SURF(N,L,K,NB,NZ,NY,NX)=0.0
900   CONTINUE
C     ARLFXB=0.0
C     ARLFXL=0.0
C     SURFXX=0.0
      DO 500 K=1,25
C     ARLFXB=ARLFXB+ARLF(K,NB,NZ,NY,NX)
      IF(ARLF(K,NB,NZ,NY,NX).GT.0.0)THEN
      DO 700 L=JC,1,-1
C     ARLFXL=ARLFXL+ARLFL(L,K,NB,NZ,NY,NX)
      DO 800 N=1,4
      SURF(N,L,K,NB,NZ,NY,NX)=AMAX1(0.0,CLASS(N,NZ,NY,NX)
     2*0.25*ARLFL(L,K,NB,NZ,NY,NX))
C     SURFXX=SURFXX+SURF(N,L,K,NB,NZ,NY,NX)
C     IF(NZ.EQ.2)THEN
C     WRITE(*,6363)'SURF',I,J,NX,NY,NZ,NB,K,L,N
C    2,ARLFL(L,K,NB,NZ,NY,NX)
C    2,SURF(N,L,K,NB,NZ,NY,NX),CLASS(N,NZ,NY,NX),ARLF(K,NB,NZ,NY,NX)
C    3,ARLFXB,ARLFXL,SURFXX,ARLF(0,NB,NZ,NY,NX)
C    4,ARLFB(NB,NZ,NY,NX),ZC(NZ,NY,NX)
6363  FORMAT(A8,9I4,20E12.4)
C     ENDIF
800   CONTINUE
700   CONTINUE
      ENDIF
500   CONTINUE
C
C     ALLOCATE STALK AREA TO INCLINATION CLASSES ACCORDING TO
C     BRANCH ANGLE ENTERED IN 'READQ' ASSUMING AZIMUTH IS UNIFORM
C
C     SURFB=stalk surface area in canopy layer used in �hour1.f�
C     ANGBR=stem angle from horizontal from PFT file
C     ARSTK=total branch stalk surface area in each layer
C
      DO 910 L=1,JC
      DO 910 N=1,4
      SURFB(N,L,NB,NZ,NY,NX)=0.0
910   CONTINUE
      IF(NB.EQ.NB1(NZ,NY,NX))THEN
      N=4
      ELSE
      N=MIN(4,INT(ASIN(ANGBR(NZ,NY,NX))/0.3927)+1)
      ENDIF
      DO 710 L=1,JC
      SURFB(N,L,NB,NZ,NY,NX)=0.25*ARSTK(L,NB,NZ,NY,NX)
710   CONTINUE
      ENDIF
C
C     SET MAXIMUM GRAIN NUMBER FROM SHOOT MASS BEFORE ANTHESIS
C
C     IDAY(3,=start of stem elongation and setting max seed number
C     IDAY(6,=start of anthesis and setting final seed number
C     GRNXB=potential number of seed set sites
C     STMX=maximum potential seed number from PFT file
C     GROSTK=stalk growth rate
C
      IF(IDAY(3,NB,NZ,NY,NX).NE.0.AND.IDAY(6,NB,NZ,NY,NX).EQ.0
     2.AND.WTSHT(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      SETM=AMAX1(0.0,DWTSHT(NZ,NY,NX)/XNFH)**0.5*XNFH
     2*WTSHTB(NB,NZ,NY,NX)/WTSHT(NZ,NY,NX)
      GRNXB(NB,NZ,NY,NX)=GRNXB(NB,NZ,NY,NX)+STMX(NZ,NY,NX)*SETM
C     WRITE(*,4246)'GRNX',I,J,NFZ,NZ,NB,IDAY(3,NB,NZ,NY,NX)
C    2,GRNXB(NB,NZ,NY,NX),STMX(NZ,NY,NX),SETM,DWTSHT(NZ,NY,NX)
C    3,WTSHTB(NB,NZ,NY,NX),WTSHT(NZ,NY,NX),XNFH
      ENDIF
C
C     SET FINAL GRAIN NUMBER FROM C,N,P NON-STRUCTURAL POOLS
C     AFTER ANTHESIS
C
C     IDAY(6,=start of anthesis and setting final seed number
C     IDAY(7,=start of grain filling and setting max seed size
C     IDAY(8,=end date for setting final seed number
C     IDAY(9,=end of setting max seed size
C     SET=seed set limited by nonstructural C,N,P
C     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concn in branch(g g-1)
C     TCC=canopy temperature
C     CTC=chilling temperature for CO2 fixation, seed loss (oC)
C        from PFT file
C     HTC=high temperature threshold for grain number loss
C     FGRNX=loss of seed set from high or low temperature
C     SSTX=sensitivity to TCC > HTC,TCC < CTC from �startq.f�
C        (seeds oC-1)
C     GRNOB=seed set number
C     GRNXB=potential number of seed set sites
C     SDMX=maximum seed number per STMX from PFT file
C     DGSTGF=change in reproductive node number normalized for maturity group
C
      IF(IDAY(6,NB,NZ,NY,NX).NE.0.AND.IDAY(9,NB,NZ,NY,NX).EQ.0)THEN
      SET=AMIN1(CCPOLB(NB,NZ,NY,NX)/(CCPOLB(NB,NZ,NY,NX)+SETC)
     2,CZPOLB(NB,NZ,NY,NX)/(CZPOLB(NB,NZ,NY,NX)+SETN)
     3,CPPOLB(NB,NZ,NY,NX)/(CPPOLB(NB,NZ,NY,NX)+SETP))
      IF(TCC(NZ,NY,NX).LT.CTC(NZ,NY,NX))THEN
      IF(IDAY(7,NB,NZ,NY,NX).EQ.0)THEN
      FGRNX=SSTX(NZ,NY,NX)*(CTC(NZ,NY,NX)-TCC(NZ,NY,NX))*XNFH
      ELSEIF(IDAY(8,NB,NZ,NY,NX).EQ.0)THEN
      FGRNX=SSTX(NZ,NY,NX)*(CTC(NZ,NY,NX)-TCC(NZ,NY,NX))*XNFH
      ELSE
      FGRNX=0.0
      ENDIF
      ELSEIF(TCC(NZ,NY,NX).GT.HTC(NZ,NY,NX))THEN
      IF(IDAY(7,NB,NZ,NY,NX).EQ.0)THEN
      FGRNX=SSTX(NZ,NY,NX)*(TCC(NZ,NY,NX)-HTC(NZ,NY,NX))*XNFH
      ELSEIF(IDAY(8,NB,NZ,NY,NX).EQ.0)THEN
      FGRNX=SSTX(NZ,NY,NX)*(TCC(NZ,NY,NX)-HTC(NZ,NY,NX))*XNFH
      ELSE
      FGRNX=0.0
      ENDIF
      ELSE
      FGRNX=0.0
      ENDIF
      IF(IDAY(6,NB,NZ,NY,NX).NE.0.AND.IDAY(8,NB,NZ,NY,NX).EQ.0)THEN
      SETG=AMIN1(SET,WFNSR)
      GRNOB(NB,NZ,NY,NX)=AMIN1(SDMX(NZ,NY,NX)*GRNXB(NB,NZ,NY,NX)
     2,GRNOB(NB,NZ,NY,NX)+(SDMX(NZ,NY,NX)*GRNXB(NB,NZ,NY,NX)
     3*SETG*DGSTGF(NB,NZ,NY,NX)-FGRNX*GRNOB(NB,NZ,NY,NX)))
C     IF(FGRNX.LT.1.0)THEN
C     WRITE(*,4246)'GRNO',I,J,NFZ,NZ,NB,IDAY(6,NB,NZ,NY,NX)
C    2,GRNOB(NB,NZ,NY,NX),GRNXB(NB,NZ,NY,NX)
C    2,TCC(NZ,NY,NX),HTC(NZ,NY,NX),FGRNX
C    3,SET,CCPOLB(NB,NZ,NY,NX),CZPOLB(NB,NZ,NY,NX)
C    4,CPPOLB(NB,NZ,NY,NX)
4246  FORMAT(A8,6I4,20E12.4)
C     ENDIF
      ENDIF
C
C     SET MAXIMUM GRAIN SIZE FROM C,N,P NON-STRUCTURAL POOLS
C     AFTER ANTHESIS
C
C     GRMX=maximum individual seed size from PFT file (g)
C     DGSTGF=change in reproductive node number normalized for
C        maturity group
C     GRWTB=individual seed size
C     SET=seed set limited by nonstructural C,N,P
C
      IF(IDAY(7,NB,NZ,NY,NX).NE.0.AND.IDAY(9,NB,NZ,NY,NX).EQ.0)THEN
      SETW=AMIN1(SET,WFNSR)**0.25
      GRWTB(NB,NZ,NY,NX)=AMIN1(GRMX(NZ,NY,NX),GRWTB(NB,NZ,NY,NX)
     2+GRMX(NZ,NY,NX)*AMAX1(0.50,SETW)*DGSTGF(NB,NZ,NY,NX))
C     IF(FGRNX.LT.1.0)THEN
C     WRITE(*,4246)'GRWT',I,J,NFZ,NZ,NB,IDAY(8,NB,NZ,NY,NX)
C    2,TCC(NZ,NY,NX)
C    2,HTC(NZ,NY,NX),FGRNX,GRMX(NZ,NY,NX),GRWTB(NB,NZ,NY,NX)
C     ENDIF
      ENDIF
      ENDIF
C
C     RESET PHENOLOGY AT EMERGENCE ('VRNS' > 'VRNL')
C     AND END OF SEASON ('VRNF' > 'VRNX')
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     IFLGE=flag for enabling leafout:0=enable,1=disable
C     VRNS,VRNL=leafout hours,hours required for leafout
C     IFLGF=flag for enabling leafoff:0=enable,1=disable
C     VRNF,VRNX=leafoff hours,hours required for leafoff
C
      IF(NFZ.EQ.1)THEN
      IF(ISTYP(NZ,NY,NX).NE.0
     2.OR.(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).GT.1))THEN
      IF((IFLGE(NB,NZ,NY,NX).EQ.0
     2.AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))
     3.OR.(IFLGF(NB,NZ,NY,NX).EQ.0
     4.AND.VRNF(NB,NZ,NY,NX).GE.VRNX(NB,NZ,NY,NX)))THEN
C
C     SPRING PHENOLOGY RESET
C
C     GROUP,GROUPI=node number required for floral initiation
C     PSTGI=node number at floral initiation
C     PSTGF=node number at flowering
C     VSTGX=leaf number on date of floral initiation
C     TGSTGI=total change in vegetative node number normalized
C        for maturity group
C     TGSTGF=total change in reproductive node number normalized
C        for maturity group
C     IDAY(1,=emergence date
C
      IF((IFLGE(NB,NZ,NY,NX).EQ.0.AND.ISTYP(NZ,NY,NX).NE.0)
     2.AND.(VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX)))THEN
      IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      GROUP(NB,NZ,NY,NX)=AMAX1(0.0,GROUPI(NZ,NY,NX)
     2-NBTB(NB,NZ,NY,NX))
      ELSE
      GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
      ENDIF
      PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
      PSTGF(NB,NZ,NY,NX)=0.0
      VSTGX(NB,NZ,NY,NX)=0.0
      TGSTGI(NB,NZ,NY,NX)=0.0
      TGSTGF(NB,NZ,NY,NX)=0.0
      IDAY(1,NB,NZ,NY,NX)=I
      DO 2005 M=2,10
      IDAY(M,NB,NZ,NY,NX)=0
2005  CONTINUE
      IF(NB.EQ.NB1(NZ,NY,NX))THEN
      WSTR(NZ,NY,NX)=0.0
      ENDIF
C
C     SPRING LEAF AND SHEATH RESET
C
C     IFLGA,IFLGE=flag for initializing,enabling leafout
C     VRNS,VRNL=leafout hours,hours required for leafout
C     PSTG=node number
C     VSTG=number of leaves appeared
C     KVSTG=integer of most recent leaf number currently growing
C     FLG4=number of hours with no grain fill
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WT*,WT*N,WT*P=branch organ C,N,P mass
C     WG,WG*N,WG*P=node organ C,N,P mass
C     organ key:LF=leaf,SHE=petiole,STK=stalk,RSV=reserve
C        HSK=husk,EAR=ear,GR=grain,SHT=shoot
C     ARLFB,ARLF=branch,node leaf area
C     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
C     HTNODX,HTNODE=stalk height,stalk internode length
C     GRNOB=seed set number
C     GRNXB=potential number of seed set sites
C     GRWTB=individual seed size
C
      IF(IFLGE(NB,NZ,NY,NX).EQ.0.AND.ISTYP(NZ,NY,NX).NE.0
     2.AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
      IF(IBTYP(NZ,NY,NX).EQ.0)THEN
      PSTG(NB,NZ,NY,NX)=XTLI(NZ,NY,NX)
      VSTG(NB,NZ,NY,NX)=0.0
      KLEAF(NB,NZ,NY,NX)=1
      KVSTG(NB,NZ,NY,NX)=1
      FLG4(NB,NZ,NY,NX)=0.0
      DO 5330 M=1,4
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)
     2+CFOPC(5,M,NZ,NY,NX)*WTLFB(NB,NZ,NY,NX)*FWODB(0,NZ)
     3+CFOPC(5,M,NZ,NY,NX)*WTSHEB(NB,NZ,NY,NX)*FWODB(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)
     2+CFOPN(5,M,NZ,NY,NX)*WTLFBN(NB,NZ,NY,NX)*FWODLN(0,NZ)
     3+CFOPN(5,M,NZ,NY,NX)*WTSHBN(NB,NZ,NY,NX)*FWODSN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)
     2+CFOPP(5,M,NZ,NY,NX)*WTLFBP(NB,NZ,NY,NX)*FWODLP(0,NZ)
     3+CFOPP(5,M,NZ,NY,NX)*WTSHBP(NB,NZ,NY,NX)*FWODSP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+CFOPC(1,M,NZ,NY,NX)*WTLFB(NB,NZ,NY,NX)*FWODB(1,NZ)
     3+CFOPC(2,M,NZ,NY,NX)*WTSHEB(NB,NZ,NY,NX)*FWODB(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+CFOPN(1,M,NZ,NY,NX)*WTLFBN(NB,NZ,NY,NX)*FWODLN(1,NZ)
     3+CFOPN(2,M,NZ,NY,NX)*WTSHBN(NB,NZ,NY,NX)*FWODSN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+CFOPP(1,M,NZ,NY,NX)*WTLFBP(NB,NZ,NY,NX)*FWODLP(1,NZ)
     3+CFOPP(2,M,NZ,NY,NX)*WTSHBP(NB,NZ,NY,NX)*FWODSP(1,NZ)
5330  CONTINUE
      ARLFB(NB,NZ,NY,NX)=0.0
      WTLFB(NB,NZ,NY,NX)=0.0
      WTLFBN(NB,NZ,NY,NX)=0.0
      WTLFBP(NB,NZ,NY,NX)=0.0
      WTSHEB(NB,NZ,NY,NX)=0.0
      WTSHBN(NB,NZ,NY,NX)=0.0
      WTSHBP(NB,NZ,NY,NX)=0.0
      DO 5335 K=0,25
      ARLF(K,NB,NZ,NY,NX)=0.0
      HTSHE(K,NB,NZ,NY,NX)=0.0
      WGLF(K,NB,NZ,NY,NX)=0.0
      WSLF(K,NB,NZ,NY,NX)=0.0
      WGLFN(K,NB,NZ,NY,NX)=0.0
      WGLFP(K,NB,NZ,NY,NX)=0.0
      WGSHE(K,NB,NZ,NY,NX)=0.0
      WSSHE(K,NB,NZ,NY,NX)=0.0
      WGSHN(K,NB,NZ,NY,NX)=0.0
      WGSHP(K,NB,NZ,NY,NX)=0.0
5335  CONTINUE
      ENDIF
      ENDIF
C
C     REPRODUCTIVE ORGANS BECOME LITTERFALL IN GRASSES, SHRUBS
C     AT START OF SEASON
C
      IF(IFLGE(NB,NZ,NY,NX).EQ.0.AND.ISTYP(NZ,NY,NX).NE.0
     2.AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
      DO 6245 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(2,M,NZ,NY,NX)
     2*(WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)+WTGRB(NB,NZ,NY,NX))
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(2,M,NZ,NY,NX)
     2*(WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX))
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(2,M,NZ,NY,NX)
     2*(WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX))
6245  CONTINUE
      WTHSKB(NB,NZ,NY,NX)=0.0
      WTEARB(NB,NZ,NY,NX)=0.0
      WTGRB(NB,NZ,NY,NX)=0.0
      WTHSBN(NB,NZ,NY,NX)=0.0
      WTEABN(NB,NZ,NY,NX)=0.0
      WTGRBN(NB,NZ,NY,NX)=0.0
      WTHSBP(NB,NZ,NY,NX)=0.0
      WTEABP(NB,NZ,NY,NX)=0.0
      WTGRBP(NB,NZ,NY,NX)=0.0
      GRNXB(NB,NZ,NY,NX)=0.0
      GRNOB(NB,NZ,NY,NX)=0.0
      GRWTB(NB,NZ,NY,NX)=0.0
C
C     REMAINING STALKS BECOME STANDING DEAD IN GRASSES, SHRUBS
C     AT START OF SEASON
C
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     WTSTDG,WTSTDN,WTSTDP=standing dead C,N,P mass
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     WGNODE,WGNODN,WGNODP=internode C,N,P mass
C
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      DO 6345 M=1,4
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)
     5*(WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX))
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)
     5*(WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX))
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)
     5*(WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX))
6345  CONTINUE
      WTSTKB(NB,NZ,NY,NX)=0.0
      WTSTBN(NB,NZ,NY,NX)=0.0
      WTSTBP(NB,NZ,NY,NX)=0.0
      WTSTXB(NB,NZ,NY,NX)=0.0
      WTSTXN(NB,NZ,NY,NX)=0.0
      WTSTXP(NB,NZ,NY,NX)=0.0
      WTRSVB(NB,NZ,NY,NX)=0.0
      WTRSBN(NB,NZ,NY,NX)=0.0
      WTRSBP(NB,NZ,NY,NX)=0.0
      DO 6340 K=0,25
      HTNODE(K,NB,NZ,NY,NX)=0.0
      HTNODX(K,NB,NZ,NY,NX)=0.0
      WGNODE(K,NB,NZ,NY,NX)=0.0
      WGNODN(K,NB,NZ,NY,NX)=0.0
      WGNODP(K,NB,NZ,NY,NX)=0.0
6340  CONTINUE
      ENDIF
      ENDIF
      ENDIF
C
C     SPRING OR FALL FLAG RESET
C
C     IFLGE=flag for enabling leafout:0=enable,1=disable
C     VRNS,VRNL=leafout hours,hours required for leafout
C     IFLGF=flag for enabling leafoff:0=enable,1=disable
C     IFLGR=flag for reproductive growth
C     IFLGQ=current hours after physlogical maturity
C        until start of litterfall
C     IFLGA=flag for initializing leafout
C
      IF(IFLGE(NB,NZ,NY,NX).EQ.0
     2.AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
      IFLGE(NB,NZ,NY,NX)=1
      IFLGF(NB,NZ,NY,NX)=0
      IFLGR(NB,NZ,NY,NX)=0
      IFLGQ(NB,NZ,NY,NX)=0
      ELSE
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).EQ.1)THEN
      IFLGE(NB,NZ,NY,NX)=1
      ELSE
      IFLGE(NB,NZ,NY,NX)=0
      ENDIF
      IFLGF(NB,NZ,NY,NX)=1
      IFLGR(NB,NZ,NY,NX)=1
      IFLGQ(NB,NZ,NY,NX)=0
      IFLGA(NB,NZ,NY,NX)=0
      ENDIF
      ENDIF
      ENDIF
C
C     REPRODUCTIVE MATERIAL BECOMES LITTERFALL AT END OF SEASON
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     FSNR=rate constant for litterfall at end of growing season (h-1)
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTHSKB,WTEARB=branch husk,ear C mass
C     WTHSBN,WTEABN=branch husk,ear N mass
C     WTHSBP,WTEABP=branch husk,ear P mass
C
      IF(IFLGR(NB,NZ,NY,NX).EQ.1)THEN
      IFLGQ(NB,NZ,NY,NX)=IFLGQ(NB,NZ,NY,NX)+1
      IF(IFLGQ(NB,NZ,NY,NX).EQ.IFLGQX)THEN
      IFLGR(NB,NZ,NY,NX)=0
      IFLGQ(NB,NZ,NY,NX)=0
      ENDIF
      DO 6330 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+FSNR*CFOPC(2,M,NZ,NY,NX)
     2*(WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX))
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+FSNR*CFOPN(2,M,NZ,NY,NX)
     2*(WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX))
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+FSNR*CFOPP(2,M,NZ,NY,NX)
     2*(WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX))
C
C     IF WINTER ANNUAL
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     FSNR=rate constant for litterfall (h-1)
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTGRB,WTGRBN,WTGRBP=total seed C,N,P mass
C
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)
     2+FSNR*CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)
     2+FSNR*CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)
     2+FSNR*CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ELSE
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+FSNR*CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+FSNR*CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+FSNR*CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ENDIF
6330  CONTINUE
      WTHSKB(NB,NZ,NY,NX)=(1.0-FSNR)*WTHSKB(NB,NZ,NY,NX)
      WTEARB(NB,NZ,NY,NX)=(1.0-FSNR)*WTEARB(NB,NZ,NY,NX)
      WTGRB(NB,NZ,NY,NX)=(1.0-FSNR)*WTGRB(NB,NZ,NY,NX)
      WTHSBN(NB,NZ,NY,NX)=(1.0-FSNR)*WTHSBN(NB,NZ,NY,NX)
      WTEABN(NB,NZ,NY,NX)=(1.0-FSNR)*WTEABN(NB,NZ,NY,NX)
      WTGRBN(NB,NZ,NY,NX)=(1.0-FSNR)*WTGRBN(NB,NZ,NY,NX)
      WTHSBP(NB,NZ,NY,NX)=(1.0-FSNR)*WTHSBP(NB,NZ,NY,NX)
      WTEABP(NB,NZ,NY,NX)=(1.0-FSNR)*WTEABP(NB,NZ,NY,NX)
      WTGRBP(NB,NZ,NY,NX)=(1.0-FSNR)*WTGRBP(NB,NZ,NY,NX)
      GRNXB(NB,NZ,NY,NX)=(1.0-FSNR)*GRNXB(NB,NZ,NY,NX)
      GRNOB(NB,NZ,NY,NX)=(1.0-FSNR)*GRNOB(NB,NZ,NY,NX)
      GRWTB(NB,NZ,NY,NX)=(1.0-FSNR)*GRWTB(NB,NZ,NY,NX)
C
C     REMAINING STALKS BECOME STANDING DEAD IN GRASSES, SHRUBS
C     AT END OF SEASON
C
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FSNR=rate constant for litterfall (h-1)
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
C
      IF((IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)
     2.AND.ISTYP(NZ,NY,NX).NE.0)THEN
      DO 6335 M=1,4
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)
     5*FSNR*(WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX))
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)
     5*FSNR*(WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX))
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)
     5*FSNR*(WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX))
6335  CONTINUE
      WTSTKB(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTKB(NB,NZ,NY,NX)
      WTSTBN(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTBN(NB,NZ,NY,NX)
      WTSTBP(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTBP(NB,NZ,NY,NX)
      WTSTXB(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTXB(NB,NZ,NY,NX)
      WTSTXN(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTXN(NB,NZ,NY,NX)
      WTSTXP(NB,NZ,NY,NX)=(1.0-FSNR)*WTSTXP(NB,NZ,NY,NX)
      WTRSVB(NB,NZ,NY,NX)=(1.0-FSNR)*WTRSVB(NB,NZ,NY,NX)
      WTRSBN(NB,NZ,NY,NX)=(1.0-FSNR)*WTRSBN(NB,NZ,NY,NX)
      WTRSBP(NB,NZ,NY,NX)=(1.0-FSNR)*WTRSBP(NB,NZ,NY,NX)
      DO 2010 K=0,25
C     HTNODE(K,NB,NZ,NY,NX)=(1.0-FSNR)*HTNODE(K,NB,NZ,NY,NX)
      HTNODX(K,NB,NZ,NY,NX)=(1.0-FSNR)*HTNODX(K,NB,NZ,NY,NX)
      WGNODE(K,NB,NZ,NY,NX)=(1.0-FSNR)*WGNODE(K,NB,NZ,NY,NX)
      WGNODN(K,NB,NZ,NY,NX)=(1.0-FSNR)*WGNODN(K,NB,NZ,NY,NX)
      WGNODP(K,NB,NZ,NY,NX)=(1.0-FSNR)*WGNODP(K,NB,NZ,NY,NX)
2010  CONTINUE
      ENDIF
C
C     SELF-SEEDING ANNUALS IF COLD OR DROUGHT DECIDUOUS
C
C     ISTYP=growth habit:0=annual,1=perennial
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     IDAYH,IYRH=day,year of harvesting
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
C     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
C          IHVST=3:reduction of clumping factor
C          IHVST=4 or 6:animal or insect biomass(g LM m-2)
C     THIN=IHVST=0-3,5: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C           leaf,non-foliar,woody, standing dead removed from PFT
C     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
C           leaf,non-foliar,woody, standing dead removed from
C         ecosystem
C     IDAY0,IYR0=day,year of planting
C     IFLGI=PFT initialization flag:0=no,1=yes
C
      IF(NB.EQ.NB1(NZ,NY,NX))THEN
C
C     IF WINTER ANNUAL
C
C     ISTYP=growth habit:0=annual,1=perennial
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      IDAYH(NZ,NY,NX)=I
      IYRH(NZ,NY,NX)=IYRC
      IHVST(NZ,I,NY,NX)=1
      JHVST(NZ,I,NY,NX)=2
      HVST(NZ,I,NY,NX)=0.0
      THIN(NZ,I,NY,NX)=0.0
      EHVST(1,1,NZ,I,NY,NX)=1.0
      EHVST(1,2,NZ,I,NY,NX)=1.0
      EHVST(1,3,NZ,I,NY,NX)=1.0
      EHVST(1,4,NZ,I,NY,NX)=1.0
      EHVST(2,1,NZ,I,NY,NX)=0.0
      EHVST(2,2,NZ,I,NY,NX)=1.0
      EHVST(2,3,NZ,I,NY,NX)=0.0
      EHVST(2,4,NZ,I,NY,NX)=0.0
      IDAY0(NZ,NY,NX)=-1E+06
      IYR0(NZ,NY,NX)=-1E+06
      IFLGI(NZ,NY,NX)=1
C     WRITE(*,3366)'HVST',I,J,NFZ,NZ,IYRC
C    2,IDAYH(NZ,NY,NX),IYRH(NZ,NY,NX)
C    2,IHVST(NZ,I,NY,NX),JHVST(NZ,I,NY,NX),IFLGI(NZ,NY,NX)
3366  FORMAT(A8,10I8)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
C     TRANSFER C,N,P FROM SEASONAL STORAGE TO SHOOT AND ROOT
C     NON-STRUCTURAL C DURING SEED GERMINATION OR LEAFOUT
C
C     IF(NZ.EQ.2)THEN
C     WRITE(*,2322)'VRNS',I,J,NX,NY,NZ,NB,NB1(NZ,NY,NX),IFLGZ
C    2,ISTYP(NZ,NY,NX),IFLGI(NZ,NY,NX),IFLGE(NB,NZ,NY,NX)
C    3,IFLGF(NB,NZ,NY,NX),IDAY0(NZ,NY,NX),IYR0(NZ,NY,NX)
C    3,VRNS(NB1(NZ,NY,NX),NZ,NY,NX),VRNL(NB,NZ,NY,NX)
C    3,VRNF(NB,NZ,NY,NX),VRNX(NB,NZ,NY,NX),WTRVC(NZ,NY,NX)
2322  FORMAT(A8,14I4,20E12.4)
C     ENDIF
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IDAY0,IYR0=day,year of planting
C     IYRC=current year
C     VRNF,VRNX=leafoff hours,hours required for leafoff
C     FVRN=fraction of hours required for leafoff to initiate
C        remobilization
C     VRNS,VRNL=leafout hours,hours required for leafout
C     WTRTD,CPOOLR=root structural,nonstructural C
C
      IF((ISTYP(NZ,NY,NX).EQ.0.AND.IFLGI(NZ,NY,NX).EQ.0)
     2.OR.(I.GE.IDAY0(NZ,NY,NX).AND.IYRC.EQ.IYR0(NZ,NY,NX)
     3.AND.VRNF(NB,NZ,NY,NX)
     3.LT.FVRN(IWTYP(NZ,NY,NX))*VRNX(NB,NZ,NY,NX))
     4.OR.(VRNS(NB1(NZ,NY,NX),NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX)
     4.AND.VRNF(NB,NZ,NY,NX)
     4.LT.FVRN(IWTYP(NZ,NY,NX))*VRNX(NB,NZ,NY,NX)))THEN
      WTRTM=0.0
      CPOOLM=0.0
      DO 4 L=NU(NY,NX),NI(NZ,NY,NX)
      WTRTM=WTRTM+AMAX1(0.0,WTRTD(1,L,NZ,NY,NX))
      CPOOLM=CPOOLM+AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX))
4     CONTINUE
C
C     RESET TIME COUNTER
C
C     ATRP=hourly leafout counter
C     IFLGA=flag for initializing leafout
C
      IF(IFLGA(NB,NZ,NY,NX).EQ.0)THEN
      ATRP(NB,NZ,NY,NX)=0.0
      IFLGA(NB,NZ,NY,NX)=1
      ENDIF
C
C     INCREMENT TIME COUNTER FOR GERMINATION, LEAFOUT
C
C     IPTYP=photoperiod type:0=day neutral,1=short day,2=long day
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     XDL=critical photoperiod (h):<0=maximum daylength from site file
C     XPPD=photoperiod sensitivity (node h-1)
C     DYLN=daylength
C     WFNSR=expansion,extension function of canopy water potential
C     TFN3=temperature function for canopy growth
C     ATRPX=number of hours required to initiate remobilization of
C        storage C for leafout
C
      IF(NB.EQ.NB1(NZ,NY,NX))THEN
      IF(IPTYP(NZ,NY,NX).EQ.2
     2.AND.(IWTYP(NZ,NY,NX).EQ.1.OR.IWTYP(NZ,NY,NX).EQ.3))THEN
      PPDX=AMAX1(0.0,XDL(NZ,NY,NX)-XPPD(NZ,NY,NX)-DYLN(NY,NX))
      ATRPPD=EXP(-0.0*PPDX)
      ELSE
      ATRPPD=1.0
      ENDIF
      DATRP=ATRPPD*TFN3(NZ,NY,NX)*WFNSG*XNFH
      ATRP(NB,NZ,NY,NX)=ATRP(NB,NZ,NY,NX)+DATRP
C     IF(NZ.EQ.2)THEN
C     WRITE(*,2323)'ATRP',I,J,NFZ,NX,NY,NZ,NB,NG(NZ,NY,NX)
C    2,NI(NZ,NY,NX)
C    2,ATRP(NB,NZ,NY,NX),DATRP,PSILT(NZ,NY,NX),SDPTHI(NZ,NY,NX)
C    2,ATRPPD,TFN3(NZ,NY,NX),WFNSG,PPDX,XDL(NZ,NY,NX),XPPD(NZ,NY,NX)
C    3,DYLN(NY,NX),WTLFB(NB,NZ,NY,NX),ARLFB(NB,NZ,NY,NX)
C    4,WTRVC(NZ,NY,NX),HTCTL(NZ,NY,NX),PSILT(NZ,NY,NX)
2323  FORMAT(A8,9I4,20E12.4)
C     ENDIF
C
C     IF GERMINATION OR LEAFOUT IS ENABLED,
C     REMOBILIZE C FROM SEASONAL STORAGE AT FIRST-ORDER RATE
C     MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
C
      IF(ATRP(NB,NZ,NY,NX).LE.ATRPX(ISTYP(NZ,NY,NX))
     2.OR.(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).EQ.0))THEN
C
C     VMXG=specific oxidation rate of storage C during leafout at 25 C
C     WTRVC=storage C
C     CH2OH=storage C oxidation rate during leafout
C     CPOOL,CPOOLR=non-structural C mass in branch,root
C     FXFC=fraction of total root mass in layer
C     FXSH,FXRT=shoot-root partitioning of storage C during leafout
C     WTRTD=root C mass
C
      IF(WTRVC(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CPOOLT=CPOOLM+CPOOL(NB,NZ,NY,NX)
      GFNX=VMXG(ISTYP(NZ,NY,NX))*DATRP
      CH2OH=AMAX1(0.0,GFNX*WTRVC(NZ,NY,NX))
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)-CH2OH
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)
     2+CH2OH*FXSH(ISTYP(NZ,NY,NX))
      IF(WTRTM.GT.ZEROP(NZ,NY,NX).AND.CPOOLM.GT.ZEROP(NZ,NY,NX))THEN
      DO 50 L=NU(NY,NX),NI(NZ,NY,NX)
      FXFC=AMAX1(0.0,WTRTD(1,L,NZ,NY,NX))/WTRTM
      CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)
     2+FXFC*CH2OH*FXRT(ISTYP(NZ,NY,NX))
50    CONTINUE
      ELSE
      CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)
     2+CH2OH*FXRT(ISTYP(NZ,NY,NX))
      ENDIF
C     IF(NZ.EQ.2)THEN
C     WRITE(*,2123)'GERM',I,J,NFZ,NX,NY,NZ,NB,NG(NZ,NY,NX)
C    2,NIX(NZ,NY,NX),GFNX,VMXG(ISTYP(NZ,NY,NX)),DATRP
C    2,CH2OH,WTRVC(NZ,NY,NX),WTRVN(NZ,NY,NX),WTRVP(NZ,NY,NX)
C    2,CPOOL(NB,NZ,NY,NX),CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)
C    3,FXSH(ISTYP(NZ,NY,NX)),FXRT(ISTYP(NZ,NY,NX))
2123  FORMAT(A8,9I4,20E12.4)
C     ENDIF
      ELSE
      CH2OH=0.0
      ENDIF
      ELSE
      CH2OH=0.0
      ENDIF
C
C     REMOBILIZE N,P FROM SEASONAL STORAGE AT FIRST-ORDER RATE
C     MODIFIED BY SOIL TEMPERATURE AT SEED DEPTH
C
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
C     UPNH4B,UPPO4B=N,P transfer from storage to shoot
C     CH2OH=storage C oxidation rate during leafout
C     FRSV=rate constant for remobilization of storage C,N,P
C        during leafout
C     FXSH=shoot partitioning of storage C during leafout
C
      IF(WTRVC(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      IF(ISTYP(NZ,NY,NX).NE.0)THEN
      CPOOLT=AMAX1(0.0,WTRVC(NZ,NY,NX)+CPOOL(NB,NZ,NY,NX))
      ZPOOLD=(WTRVN(NZ,NY,NX)*CPOOL(NB,NZ,NY,NX)
     2-ZPOOL(NB,NZ,NY,NX)*WTRVC(NZ,NY,NX))/CPOOLT
      PPOOLD=(WTRVP(NZ,NY,NX)*CPOOL(NB,NZ,NY,NX)
     2-PPOOL(NB,NZ,NY,NX)*WTRVC(NZ,NY,NX))/CPOOLT
      UPNH4B=AMAX1(0.0,FRSV(IBTYP(NZ,NY,NX))*ZPOOLD)
      UPPO4B=AMAX1(0.0,FRSV(IBTYP(NZ,NY,NX))*PPOOLD)
      ELSE
      UPNH4B=AMAX1(0.0,FXSH(ISTYP(NZ,NY,NX))
     2*CH2OH*WTRVN(NZ,NY,NX)/WTRVC(NZ,NY,NX))
      UPPO4B=AMAX1(0.0,FXSH(ISTYP(NZ,NY,NX))
     2*CH2OH*WTRVP(NZ,NY,NX)/WTRVC(NZ,NY,NX))
      ENDIF
      ELSE
      UPNH4B=AMAX1(0.0,FXSH(ISTYP(NZ,NY,NX))*WTRVN(NZ,NY,NX))
      UPPO4B=AMAX1(0.0,FXSH(ISTYP(NZ,NY,NX))*WTRVP(NZ,NY,NX))
      ENDIF
C
C     ADD TO NON-STRUCTURAL POOLS IN ROOT
C
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     UPNH4R,UPPO4R=N,P transfer from storage to root
C     FRSV=rate constant for remobilization of storage C,N,P
C        during leafout
C     FXRT=root partitioning of storage C during leafout
C
      CPOOLM=0.0
      ZPOOLM=0.0
      PPOOLM=0.0
      DO 3 L=NU(NY,NX),NI(NZ,NY,NX)
      CPOOLM=CPOOLM+AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX))
      ZPOOLM=ZPOOLM+AMAX1(0.0,ZPOOLR(1,L,NZ,NY,NX))
      PPOOLM=PPOOLM+AMAX1(0.0,PPOOLR(1,L,NZ,NY,NX))
3     CONTINUE
      IF(WTRVC(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      IF(ISTYP(NZ,NY,NX).NE.0)THEN
      CPOOLT=AMAX1(ZEROP(NZ,NY,NX),WTRVC(NZ,NY,NX)+CPOOLM)
      ZPOOLD=(WTRVN(NZ,NY,NX)*CPOOLM
     2-ZPOOLM*WTRVC(NZ,NY,NX))/CPOOLT
      PPOOLD=(WTRVP(NZ,NY,NX)*CPOOLM
     2-PPOOLM*WTRVC(NZ,NY,NX))/CPOOLT
      UPNH4R=AMAX1(0.0,FRSV(IBTYP(NZ,NY,NX))*ZPOOLD)
      UPPO4R=AMAX1(0.0,FRSV(IBTYP(NZ,NY,NX))*PPOOLD)
C     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,9878)'GERM1',I,J,NZ,UPNH4R,FRSV(IBTYP(NZ,NY,NX))
C    2,ZPOOLD,WTRVN(NZ,NY,NX),CPOOLM,ZPOOLM,WTRVC(NZ,NY,NX)
C    3,CPOOLT
9878  FORMAT(A8,3I4,12E24.16)
C     ENDIF
      ELSE
      UPNH4R=AMAX1(0.0,FXRT(ISTYP(NZ,NY,NX))
     2*CH2OH*WTRVN(NZ,NY,NX)/WTRVC(NZ,NY,NX))
      UPPO4R=AMAX1(0.0,FXRT(ISTYP(NZ,NY,NX))
     2*CH2OH*WTRVP(NZ,NY,NX)/WTRVC(NZ,NY,NX))
      ENDIF
      ELSE
      UPNH4R=AMAX1(0.0,FXRT(ISTYP(NZ,NY,NX))*WTRVN(NZ,NY,NX))
      UPPO4R=AMAX1(0.0,FXRT(ISTYP(NZ,NY,NX))*WTRVP(NZ,NY,NX))
      ENDIF
C
C     TRANSFER STORAGE FLUXES
C
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
C     UPNH4B,UPPO4B=N,P transfer from storage to shoot
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     UPNH4R,UPPO4R=N,P transfer from storage to root
C     FXFN=root layer allocation
C
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)-UPNH4B-UPNH4R
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)-UPPO4B-UPPO4R
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+UPNH4B
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+UPPO4B
      IF(WTRTM.GT.ZEROP(NZ,NY,NX)
     2.AND.CPOOLM.GT.ZEROP(NZ,NY,NX))THEN
      DO 51 L=NU(NY,NX),NI(NZ,NY,NX)
      FXFN=AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX))/CPOOLM
C     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,9879)'GERM2',I,J,NZ,L,UPNH4R,FXFN
C    2,ZPOOLR(1,L,NZ,NY,NX),CPOOLR(1,L,NZ,NY,NX),CPOOLM
9879  FORMAT(A8,4I4,12E24.16)
C     ENDIF
      ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)+FXFN*UPNH4R
      PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)+FXFN*UPPO4R
51    CONTINUE
      ELSE
C     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,9879)'GERM3',I,J,NZ,L,UPNH4R,FXFN
C    2,ZPOOLR(1,L,NZ,NY,NX),CPOOLR(1,L,NZ,NY,NX),CPOOLM
C     ENDIF
      ZPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=ZPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)
     2+UPNH4R
      PPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=PPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)
     2+UPPO4R
      ENDIF
      ENDIF
C
C     REDISTRIBUTE TRANFERRED C FROM MAIN STEM TO OTHER BRANCHES
C
C     ATRP=hourly leafout counter
C     TFN3=temperature function for canopy growth
C     ATRPX=number of hours required for remobilization of storage C
C        during leafout
C     WFNSG=growth function of canopy water potential
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
C     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
C
      IF(NB.NE.NB1(NZ,NY,NX).AND.ATRP(NB,NZ,NY,NX)
     2.LE.ATRPX(ISTYP(NZ,NY,NX)))THEN
      ATRP(NB,NZ,NY,NX)=ATRP(NB,NZ,NY,NX)+TFN3(NZ,NY,NX)*WFNSG*XNFH
      XFRC=AMAX1(0.0,0.05*TFN3(NZ,NY,NX)
     2*(0.5*(CPOOL(NB1(NZ,NY,NX),NZ,NY,NX)+CPOOL(NB,NZ,NY,NX))
     3-CPOOL(NB,NZ,NY,NX)))
      XFRN=AMAX1(0.0,0.05*TFN3(NZ,NY,NX)
     2*(0.5*(ZPOOL(NB1(NZ,NY,NX),NZ,NY,NX)+ZPOOL(NB,NZ,NY,NX))
     2-ZPOOL(NB,NZ,NY,NX)))
      XFRP=AMAX1(0.0,0.05*TFN3(NZ,NY,NX)
     2*(0.5*(PPOOL(NB1(NZ,NY,NX),NZ,NY,NX)+PPOOL(NB,NZ,NY,NX))
     3-PPOOL(NB,NZ,NY,NX)))
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+XFRC
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+XFRP
      CPOOL(NB1(NZ,NY,NX),NZ,NY,NX)=CPOOL(NB1(NZ,NY,NX),NZ,NY,NX)-XFRC
      ZPOOL(NB1(NZ,NY,NX),NZ,NY,NX)=ZPOOL(NB1(NZ,NY,NX),NZ,NY,NX)-XFRN
      PPOOL(NB1(NZ,NY,NX),NZ,NY,NX)=PPOOL(NB1(NZ,NY,NX),NZ,NY,NX)-XFRP
      ENDIF
      ENDIF
C
C     TRANSFER LEAF AND STALK NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
C     IN PERENNIALS AFTER GRAIN FILL IN DETERMINATES, AFTER AUTUMNIZ'N
C     IN INDETERMINATES, OR AFTER SUSTAINED WATER STRESS
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IFLGZ=remobilization flag
C     WVSTKB=stalk sapwood mass
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C     FXFB=rate constant for plant-storage nonstructural C,N,P
C        exchange
C     XFRCX,XFRNX,XFRPX=nonstructural C,N,P transfer unconstrained
C        by C,N,P
C     XFRC,XFRN,XFRP=nonstructural C,N,P transfer constrained
C        by C,N,P
C     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concentration in branch
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concentration in branch
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
C     XNFH=time step from �wthr.f�
C
      IF(IFLGZ.EQ.1.AND.ISTYP(NZ,NY,NX).NE.0)THEN
      XFRCX=FXFB(IBTYP(NZ,NY,NX))
     2*AMAX1(0.0,WTRSVB(NB,NZ,NY,NX))*XNFH
      XFRNX=FXFB(IBTYP(NZ,NY,NX))
     2*AMAX1(0.0,WTRSBN(NB,NZ,NY,NX))*XNFH
      XFRPX=FXFB(IBTYP(NZ,NY,NX))
     2*AMAX1(0.0,WTRSBP(NB,NZ,NY,NX))*XNFH
      XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
      XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN)
      XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN)
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)-XFRC
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+XFRC
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)-XFRN
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+XFRN
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)-XFRP
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+XFRP
      XFRCX=FXFB(IBTYP(NZ,NY,NX))
     2*AMAX1(0.0,CPOOL(NB,NZ,NY,NX))*XNFH
      XFRNX=FXFB(IBTYP(NZ,NY,NX))
     2*AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))*XNFH
      XFRPX=FXFB(IBTYP(NZ,NY,NX))
     2*AMAX1(0.0,PPOOL(NB,NZ,NY,NX))*XNFH
      XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
      XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN)
      XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN)
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)-XFRC
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+XFRC
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-XFRN
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-XFRP
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+XFRP
C     IF(NZ.EQ.1)THEN
C     WRITE(*,4490)'RSV',I,J,NZ,NB,XFRC,XFRN,WTRSVB(NB,NZ,NY,NX)
C    2,WTRSBN(NB,NZ,NY,NX),WTRVC(NZ,NY,NX),WTRVN(NZ,NY,NX)
C    3,CNR,CNL,CPOOL(NB,NZ,NY,NX),ZPOOL(NB,NZ,NY,NX)
C    4,FXFB(IBTYP(NZ,NY,NX))
4490  FORMAT(A8,4I4,20E12.4)
C     ENDIF
      ENDIF
C
C     TRANSFER NON-STRUCTURAL C,N,P FROM LEAVES TO RESERVES
C     IN STALKS DURING GRAIN FILL IN ANNUALS OR BETWEEN LEAVES AND
C     STALK RESERVES IN PERENNIALS ACCORDING TO CONCENTRATION
C     DIFFERENCES
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IDAY(3,=start of stem elongation and setting max seed number
C     IDAY(8,=end date for setting final seed number
C     WTGRB=total seed C mass
C     GRWTB=individual seed size
C     GRNOB=seed set number
C     WTLSB=leaf+petiole mass
C     WVSTKB=stalk sapwood mass
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     FXFY,FXFZ=rate constant for plant-reserve nonstructural C,N,P
C        exchange
C     XFRC,XFRN,XFRC=nonstructural C,N,P transfer
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C
      IF((ISTYP(NZ,NY,NX).EQ.0.AND.IDAY(8,NB,NZ,NY,NX).NE.0)
     4.OR.(ISTYP(NZ,NY,NX).NE.0.AND.IDAY(3,NB,NZ,NY,NX).NE.0))THEN
      IF(WTLSB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WTPLTT=WTLSB(NB,NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
      CPOOLT=CPOOL(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
      CPOOLD=(CPOOL(NB,NZ,NY,NX)*WVSTKB(NB,NZ,NY,NX)
     2-WTRSVB(NB,NZ,NY,NX)*WTLSB(NB,NZ,NY,NX))/WTPLTT
      XFRC=FXFY(ISTYP(NZ,NY,NX))*CPOOLD
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)-XFRC
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+XFRC
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLD=(ZPOOL(NB,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX)
     2-WTRSBN(NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX))/CPOOLT
      PPOOLD=(PPOOL(NB,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX)
     2-WTRSBP(NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX))/CPOOLT
      XFRN=FXFZ(ISTYP(NZ,NY,NX))*ZPOOLD
      XFRP=FXFZ(ISTYP(NZ,NY,NX))*PPOOLD
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-XFRN
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-XFRP
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+XFRP
      ENDIF
C     IF(NZ.EQ.2)THEN
C     WRITE(*,4488)'EXCHC',I,J,NFZ,NX,NY,NZ,NB,XFRC,XFRN
C    2,WTRSVB(NB,NZ,NY,NX),WVSTKB(NB,NZ,NY,NX),CPOOL(NB,NZ,NY,NX)
C    3,WTLSB(NB,NZ,NY,NX)
C    4,CPOOLT,CPOOLD,ZPOOLD,ZPOOL(NB,NZ,NY,NX),WTRSBN(NB,NZ,NY,NX)
4488  FORMAT(A8,7I4,20E12.4)
C     ENDIF
      ENDIF
      ENDIF
C
C     REPLENISH BRANCH NON-STRUCTURAL POOL FROM
C     SEASONAL STORAGE POOL
C
C     WVSTKB,WVSTK=stalk,total stalk sapwood mass
C     WTRT=total root mass
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     XFRX=maximum storage C content for remobilization from
C        stalk,root reserves
C     XFRC=C transfer
C
      IF(WVSTKB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WVSTK(NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     3.AND.WTRT(NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     4.AND.WTRSVB(NB,NZ,NY,NX).LE.XFRX*WVSTKB(NB,NZ,NY,NX))THEN
      FWTBR=WVSTKB(NB,NZ,NY,NX)/WVSTK(NZ,NY,NX)
      WVSTBX=WVSTKB(NB,NZ,NY,NX)
      WTRTTX=WTRT(NZ,NY,NX)*FWTBR
      WTPLTT=WVSTBX+WTRTTX
      WTRSBX=AMAX1(0.0,WTRSVB(NB,NZ,NY,NX))
      WTRVCX=AMAX1(0.0,WTRVC(NZ,NY,NX)*FWTBR)
      CPOOLD=(WTRVCX*WVSTBX-WTRSBX*WTRTTX)/WTPLTT
      XFRC=AMAX1(0.0,XFRY*CPOOLD)
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+XFRC
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)-XFRC
      ENDIF
C
C     TRANSFER NON-STRUCTURAL C,N,P FROM ROOTS TO RESERVES
C     IN STALKS DURING GRAIN FILL IN ANNUALS OR BETWEEN LEAVES AND
C     STALK RESERVES IN PERENNIALS ACCORDING TO CONCENTRATION
C     DIFFERENCES
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IDAY(8,=end date for setting final seed number
C     WTRTL=active root C mass
C     FWODR=root woody fraction
C     WVSTKB,WTRSVB=stalk sapwood C
C     WTRSVB,WTRSBN,WTRSBP=stalk nonstructural C,N,P mass
C     CPOOLR,ZPOOLR,PPOOLR=root nonstructural C,N,P
C     XFRC,XFRN,XFRP=root-reserve nonstructural C,N,P transfer
C
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IDAY(8,NB,NZ,NY,NX).NE.0)THEN
      DO 2050 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IF(WTRTL(1,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WTRTRX=AMAX1(ZEROP(NZ,NY,NX),WTRTL(1,L,NZ,NY,NX)*FWODR(1,NZ))
      WTPLTX=WTRTRX+WVSTKB(NB,NZ,NY,NX)
      CPOOLD=(CPOOLR(1,L,NZ,NY,NX)*WVSTKB(NB,NZ,NY,NX)
     2-WTRSVB(NB,NZ,NY,NX)*WTRTRX)/WTPLTX
      XFRC=AMAX1(0.0,FXFY(ISTYP(NZ,NY,NX))*CPOOLD)
      CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)-XFRC
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+XFRC
      CPOOLT=CPOOLR(1,L,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLD=(ZPOOLR(1,L,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX)
     2-WTRSBN(NB,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
      PPOOLD=(PPOOLR(1,L,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX)
     2-WTRSBP(NB,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
      XFRN=AMAX1(0.0,FXFZ(ISTYP(NZ,NY,NX))*ZPOOLD)
      XFRP=AMAX1(0.0,FXFZ(ISTYP(NZ,NY,NX))*PPOOLD)
      ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)-XFRN
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+XFRN
      PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)-XFRP
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+XFRP
C     IF(L.EQ.12)THEN
C     WRITE(*,4489)'EXCHR',I,J,NFZ,NZ,NB,L,XFRC,XFRN
C    2,WTRSVB(NB,NZ,NY,NX),WVSTKB(NB,NZ,NY,NX),CPOOLR(1,L,NZ,NY,NX)
C    3,WTRTL(1,L,NZ,NY,NX),FWOOD(1,NZ),WTRTRX,WTPLTX
C    4,CPOOLT,CPOOLD,XFRC,FXFZ(ISTYP(NZ,NY,NX))
4489  FORMAT(A8,6I4,29E12.4)
C     ENDIF
C     IF(NZ.EQ.1.OR.NZ.EQ.4)THEN
C     WRITE(*,4489)'EXCHN',I,J,NZ,NB,L,WTRSBN(NB,NZ,NY,NX)
C    2,WTRSVB(NB,NZ,NY,NX),ZPOOLR(1,L,NZ,NY,NX)
C    3,CPOOLR(1,L,NZ,NY,NX),FWOOD(1,NZ),ZPOOLD,XFRN
C     ENDIF
      ENDIF
      ENDIF
      ENDIF
2050  CONTINUE
      ENDIF
C
C     GRAIN FILL BY TRANSLOCATION FROM STALK RESERVES
C     UNTIL GRAIN SINK (=FINAL GRAIN NUMBER X MAXIMUM
C     GRAIN SIZE) IS FILLED OR RESERVES ARE EXHAUSTED
C
C     IDAY(7,=start of grain filling and setting max seed size
C     WTGRB=total seed C mass
C     GRWTB=individual seed size
C     GRNOB=seed set number
C     GROLM=maximum grain fill rate
C     GFILL=grain filling rate at 25 oC from PFT file
C     TFN3=temperature function for canopy growth
C     TFN4=temperature function for root growth
C
      IF(IDAY(7,NB,NZ,NY,NX).NE.0)THEN
      IF(WTGRB(NB,NZ,NY,NX).GE.GRWTB(NB,NZ,NY,NX)
     2*GRNOB(NB,NZ,NY,NX))THEN
      GROLM=0.0
      ELSEIF(IRTYP(NZ,NY,NX).EQ.0)THEN
      GROLM=AMAX1(0.0,GFILL(NZ,NY,NX)*GRNOB(NB,NZ,NY,NX)
     2*SQRT(TFN3(NZ,NY,NX))*XNFH)
      ELSE
      GROLM=AMAX1(0.0,GFILL(NZ,NY,NX)*GRNOB(NB,NZ,NY,NX)
     2*SQRT(TFN4(NG(NZ,NY,NX),NZ,NY,NX))*XNFH)
      ENDIF
C
C     GRAIN FILL RATE MAY BE CONSTRAINED BY HIGH GRAIN C:N OR C:P
C
C     WTGRB,WTGRBN,WTGRBP=total seed C,N,P mass
C     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
C     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
C     GROLM,GROLC=maximum,actual grain fill rate
C     XLOCM,XLOCC=maximum,actual C translocation rate from reserve
C        to grain
C
      IF(WTGRBN(NB,NZ,NY,NX).LT.ZPGRM*CNGR(NZ,NY,NX)
     2*WTGRB(NB,NZ,NY,NX).OR.WTGRBP(NB,NZ,NY,NX).LT.ZPGRM
     3*CPGR(NZ,NY,NX)*WTGRB(NB,NZ,NY,NX))THEN
      GROLC=0.0
      ELSE
      GROLC=GROLM
      ENDIF
      XLOCM=AMIN1(GROLM,WTRSVB(NB,NZ,NY,NX))
      XLOCC=AMIN1(GROLC,WTRSVB(NB,NZ,NY,NX))
C
C     GRAIN N OR P FILL RATE MAY BE LIMITED BY LOW C:N OR C:P RATIOS
C     OF STALK RESERVES
C
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     ZNPGN,ZNPGP=effect of reserve N:C,P:C on grain fill N:C,P:C
C     SETN,SETP=Km for nonstructural N,P concn on seed set (g g-1)
C     ZPGRM=min N:C,P:C in grain relative to max values from PFT file
C     ZPGRD=1.0-ZPGRM
C     ZPGRN,ZPGRP=N:C,P:C ratios during grain fill
C     XLOCM,XLOCC=maximum,actual C translocation rate from reserve
C        to grain
C     CNGR,CPGR=maximum N:C,P:C ratios in grain from PFT file
C     XLOCN,XLOCP=N,P translocation rate from reserve to grain
C
      IF(WTRSBN(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      ZNPGN=WTRSBN(NB,NZ,NY,NX)/(WTRSBN(NB,NZ,NY,NX)
     2+SETN*WTRSVB(NB,NZ,NY,NX))
      ZPGRN=ZPGRM+ZPGRD*AMAX1(0.0,AMIN1(1.0,ZNPGN))
      XLOCN=AMIN1(XLOCM*CNGR(NZ,NY,NX)
     2,AMAX1(0.0,WTRSBN(NB,NZ,NY,NX)*ZPGRN)
     3,(WTGRB(NB,NZ,NY,NX)+XLOCC)*CNGR(NZ,NY,NX)-WTGRBN(NB,NZ,NY,NX))
      ELSE
      XLOCN=0.0
      ENDIF
      IF(WTRSBP(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      ZNPGP=WTRSBP(NB,NZ,NY,NX)/(WTRSBP(NB,NZ,NY,NX)
     3+SETP*WTRSVB(NB,NZ,NY,NX))
      ZPGRP=ZPGRM+ZPGRD*AMAX1(0.0,AMIN1(1.0,ZNPGP))
      XLOCP=AMIN1(XLOCM*CPGR(NZ,NY,NX)
     2,AMAX1(0.0,WTRSBP(NB,NZ,NY,NX)*ZPGRP)
     3,(WTGRB(NB,NZ,NY,NX)+XLOCC)*CPGR(NZ,NY,NX)-WTGRBP(NB,NZ,NY,NX))
      ELSE
      XLOCP=0.0
      ENDIF
C     IF(NX.EQ.1.AND.NY.EQ.6.AND.NZ.EQ.3)THEN
C     WRITE(*,85)'XLOC',I,J,NFZ,NZ,NB
C    2,WTGRB(NB,NZ,NY,NX),WTGRBN(NB,NZ,NY,NX)
C    2,WTRSVB(NB,NZ,NY,NX),WTRSBN(NB,NZ,NY,NX),XLOCC,XLOCN,XLOCP,XLOCM
C    3,CNGR(NZ,NY,NX),ZPGRX,ZNPG,GROLC,GROLM,GROGR,GROGRN
C    3,XLOCM*CNGR(NZ,NY,NX),AMAX1(0.0,WTRSBN(NB,NZ,NY,NX)*ZPGRX)
C    4,(WTGRB(NB,NZ,NY,NX)+XLOCC)*CNGR(NZ,NY,NX)-WTGRBN(NB,NZ,NY,NX)
C    4,GRNOB(NB,NZ,NY,NX),GRWTB(NB,NZ,NY,NX),GFILL(NZ,NY,NX)
C    5,SQRT(TFN3(NZ,NY,NX)),FLG4(NB,NZ,NY,NX)
85    FORMAT(A8,5I4,30E12.4)
C     ENDIF
C
C     TRANSLOCATE C,N,P FROM STALK RESERVES TO GRAIN
C
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     GROGR=grain growth rate
C     XLOCC,XLOCN,XLOCP=C,N,P translocation rate from reserve to grain
C
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)+GROGR-XLOCC
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)+GROGRN-XLOCN
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)+GROGRP-XLOCP
      WTGRB(NB,NZ,NY,NX)=WTGRB(NB,NZ,NY,NX)+XLOCC
      WTGRBN(NB,NZ,NY,NX)=WTGRBN(NB,NZ,NY,NX)+XLOCN
      WTGRBP(NB,NZ,NY,NX)=WTGRBP(NB,NZ,NY,NX)+XLOCP
      ELSE
      XLOCC=0.0
      XLOCN=0.0
      XLOCP=0.0
      ENDIF
C
C     SET DATE OF PHYSIOLOGICAL MATURITY WHEN GRAIN FILL
C     HAS STOPPED FOR SET PERIOD OF TIME
C
C     IDAY(8,=end date for setting final seed number
C     XLOCC=C translocation rate from reserve to grain
C     PP=PFT population
C     FLG4=number of hours with no grain fill
C     FLG4X=number of hours with no grain filling until
C        physiological maturity
C     IDAY(10,=date of physiological maturity
C
      IF(IDAY(8,NB,NZ,NY,NX).NE.0)THEN
      IF(XLOCC.LE.1.0E-09*PP(NZ,NY,NX))THEN
      FLG4(NB,NZ,NY,NX)=FLG4(NB,NZ,NY,NX)+1.0*XNFH
      ELSE
      FLG4(NB,NZ,NY,NX)=0.0
      ENDIF
      IF(FLG4(NB,NZ,NY,NX).GE.FLG4X)THEN
      IF(IDAY(10,NB,NZ,NY,NX).EQ.0)THEN
      IDAY(10,NB,NZ,NY,NX)=I
      ENDIF
      ENDIF
C
C     TERMINATE ANNUALS AFTER GRAIN FILL
C
C     ISTYP=growth habit:0=annual,1=perennial
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     FLG4=number of hours with no grain fill
C     FLG4X=number of hours with no grain filling until
C        physiological maturity
C     FLG4Y=number of hours after physiological maturity required
C        for senescence
C     VRNF,VRNX=leafoff hours,hours required for leafoff
C
C     IF WINTER ANNUAL
C
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      IF(FLG4(NB,NZ,NY,NX).GT.FLG4X+FLG4Y(IWTYP(NZ,NY,NX)))THEN
      VRNF(NB,NZ,NY,NX)=VRNX(NB,NZ,NY,NX)+0.5
      ENDIF
      ENDIF
      ENDIF
C
C     CANOPY N2 FIXATION (CYANOBACTERIA)
C
C     INTYP=N2 fixation: 4,5,6=rapid to slow canopy symbiosis
C
      IF(INTYP(NZ,NY,NX).GE.4)THEN
C
C     INITIAL INFECTION
C
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     WTNDI=initial bacterial mass at infection
C     AREA=grid cell area
C     CNND,CPND=bacterial N:C,P:C ratio from PFT file
C
      IF(NFZ.EQ.1.AND.WTNDB(NB,NZ,NY,NX).LE.0.0
     2.AND.ICHKF.EQ.0)THEN
      WTNDB(NB,NZ,NY,NX)=WTNDB(NB,NZ,NY,NX)
     2+WTNDI*AREA(3,NU(NY,NX),NY,NX)
      WTNDBN(NB,NZ,NY,NX)=WTNDBN(NB,NZ,NY,NX)
     2+WTNDI*AREA(3,NU(NY,NX),NY,NX)*CNND(NZ,NY,NX)
      WTNDBP(NB,NZ,NY,NX)=WTNDBP(NB,NZ,NY,NX)
     2+WTNDI*AREA(3,NU(NY,NX),NY,NX)*CPND(NZ,NY,NX)
      ENDIF
C
C     CONSTRAINTS ON BACTERIAL ACTIVITY IMPOSED BY
C     C,P AND BIOMASS
C
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P mass in bacteria
C     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concentration
C        in bacteria
C     ZCPOLI,ZPPOLI=nonstructural N:C,N:P
C     ZCKI,ZPKI=nonstructural C,P inhibition on bacterial growth,
C        nutrient recycling
C     FCNPF=N,P constraint to bacterial growth
C     CNDLB=bacteria:leaf+petiole ratio
C
      IF(WTNDB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CCPOLN=AMAX1(0.0,CPOLNB(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX))
      CZPOLN=AMAX1(0.0,ZPOLNB(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX))
      CPPOLN=AMAX1(0.0,PPOLNB(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX))
      ELSE
      CCPOLN=1.0
      CZPOLN=1.0
      CPPOLN=1.0
      ENDIF
      IF(CCPOLN.GT.ZERO)THEN
      ZCPOLI=CZPOLN/CCPOLN
      ELSE
      ZCPOLI=0.0
      ENDIF
      IF(CPPOLN.GT.ZERO)THEN
      ZPPOLI=CZPOLN/CPPOLN
      ELSE
      ZPPOLI=0.0
      ENDIF
      IF(CCPOLN.GT.ZERO)THEN
      CCC=AMAX1(0.0,AMIN1(1.0
     1,CZPOLN/(CZPOLN+CCPOLN*CNKI)
     2,CPPOLN/(CPPOLN+CCPOLN*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLN/(CCPOLN+CZPOLN/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLN/(CCPOLN+CPPOLN/CPKI)))
      ELSE
      CCC=0.0
      CNC=0.0
      CPC=0.0
      ENDIF
      IF(WTNDB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FCNPF=AMIN1(1.0
     2,SQRT(WTNDBN(NB,NZ,NY,NX)/(WTNDB(NB,NZ,NY,NX)*CNND(NZ,NY,NX)))
     3,SQRT(WTNDBP(NB,NZ,NY,NX)/(WTNDB(NB,NZ,NY,NX)*CPND(NZ,NY,NX))))
      ELSE
      FCNPF=1.0
      ENDIF
      IF(ARLFB(NB,NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
      CNDLB=AMAX1(0.0,WTNDB(NB,NZ,NY,NX)/ARLFB(NB,NZ,NY,NX))
      ELSE
      CNDLB=0.0
      ENDIF
C
C     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
C     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
C     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
C     AND TEMPERATURE WITH EXCESS N FEEDBACK
C
C     RCNDL=respiration from non-structural C
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     VMXO=specific respiration rate by bacterial N2 fixers
C     WTNDB=bacterial C mass
C     TFN3=temperature function for canopy growth
C     FCNPF=N,P constraint to bacterial activity
C     WFNSG=growth function of canopy water potential
C     ZCPOLI,ZPPOLI=nonstructural N:C,N:P
C
      RCNDL=AMAX1(0.0,AMIN1(CPOLNB(NB,NZ,NY,NX)
     2,VMXO*WTNDB(NB,NZ,NY,NX))*FCNPF*TFN3(NZ,NY,NX)*WFNSG*XNFH)
     2/(1.0+AMAX1(ZCPOLI/ZCKI,ZPPOLI/ZPKI))
C
C     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
C     NODULE STRUCTURAL N
C
C     RMNDL=bacterial maintenance respiration
C     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
C     TFN5=temperature function for canopy maintenance respiration
C     WTNDBN=bacterial N mass
C     WFNSR=expansion,extension function of canopy water potential
C
      RMNDL=AMAX1(0.0,RMPLT*TFN5*WTNDBN(NB,NZ,NY,NX)*WFNSR*XNFH)
C
C     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
C     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
C
C     RXNDL=difference between non-structural C respiration and
C        maintenance respiration
C     RGNDL=growth respiration unlimited by N,P
C     RSNDL=excess maintenance respiration
C
      RXNDL=RCNDL-RMNDL
      RGNDL=AMAX1(0.0,RXNDL)
      RSNDL=AMAX1(0.0,-RXNDL)
C
C     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
C     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
C     CONSTRAINED BY MICROBIAL N REQUIREMENT
C
C     RGN2P=respiration requirement to maintain bacterial N:C ratio
C     WTNDB,WTNDBN=bacterial C,N mass
C     CNND=bacterial N:C ratio from PFT file
C     EN2F=N fixation yield from C oxidation (g N g-1 C)
C     RGNDL=growth respiration unlimited by N,P
C     RGN2F=respiration for N2 fixation
C     RUPNFB,UPNFC=branch,total N2 fixation
C
      RGN2P=AMAX1(0.0,WTNDB(NB,NZ,NY,NX)*CNND(NZ,NY,NX)
     2-WTNDBN(NB,NZ,NY,NX))/EN2F
      IF(RGNDL.GT.ZEROP(NZ,NY,NX))THEN
      RGN2F=RGNDL*RGN2P/(RGNDL+RGN2P)
      ELSE
      RGN2F=0.0
      ENDIF
      RUPNFB=RGN2F*EN2F
      UPNFC(NZ,NY,NX)=UPNFC(NZ,NY,NX)+RUPNFB
C
C     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION AND LEAKAGE
C
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     RCCZN,RCCYN=min,max fractions for bacteria C recycling
C     RCCXN,RCCQN=max fractions for bacteria N,P recycling
C     SPNDX=specific bacterial decomposition rate at current CNDLB
C     SPNDL=specific decomposition rate by bacterial N2 fixers
C     CNDLB=bacteria:leaf+petiole ratio
C     CNDLI=bacteria: leaf+petiole ratio controlling bacterial
C        decomposition
C     TFN3=temperature function for canopy growth
C     WFNSG=growth function of canopy water potential
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
C     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
C     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
C
      RCCC=RCCZN+CCC*RCCYN
      RCCN=CNC*RCCXN
      RCCP=CPC*RCCQN
      SPNDX=AMIN1(1.0,SPNDL*CNDLB/CNDLI(INTYP(NZ,NY,NX)))
     2*SQRT(TFN3(NZ,NY,NX)*WFNSG)*XNFH
      RXNDLC=SPNDX*WTNDB(NB,NZ,NY,NX)
      RXNDLN=SPNDX*WTNDBN(NB,NZ,NY,NX)
      RXNDLP=SPNDX*WTNDBP(NB,NZ,NY,NX)
      RDNDLC=RXNDLC*(1.0-RCCC)
      RDNDLN=RXNDLN*(1.0-RCCC)*(1.0-RCCN)
      RDNDLP=RXNDLP*(1.0-RCCC)*(1.0-RCCP)
      RCNDLC=RXNDLC-RDNDLC
      RCNDLN=RXNDLN-RDNDLN
      RCNDLP=RXNDLP-RDNDLP
C
C     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
C     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
C     ENTERED IN 'READQ'
C
C     CGNDL=total non-structural C used in bacterial growth
C        and growth respiration
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     RMNDL=bacterial maintenance respiration
C     RCNDL=respiration from non-structural C
C     RCNDLC=bacterial C decomposition to recycling
C     RGNDL=growth respiration ltd by O2
C     RGN2F=respiration for N2 fixation
C     GRNDG=bacterial growth
C     DMND=bacterial growth yield
C     RGNDG=bacterial respiration for growth and N2 fixation
C     ZADDN,PADDN=nonstructural N,P used in growth
C     CNND,CPND=bacterial N:C,P:C ratio from PFT file
C     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concentration
C        in bacteria
C     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
C
      CGNDL=AMIN1(CPOLNB(NB,NZ,NY,NX)-AMIN1(RMNDL,RCNDL)
     2-RGN2F+RCNDLC,(RGNDL-RGN2F)/(1.0-DMND(NZ,NY,NX)))
      GRNDG=CGNDL*DMND(NZ,NY,NX)
      RGNDG=RGN2F+CGNDL*(1.0-DMND(NZ,NY,NX))
      ZADDN=AMAX1(0.0,AMIN1(ZPOLNB(NB,NZ,NY,NX)
     2,GRNDG*CNND(NZ,NY,NX)))*CZPOLN/(CZPOLN+CZKM)
      PADDN=AMAX1(0.0,AMIN1(PPOLNB(NB,NZ,NY,NX)
     2,GRNDG*CPND(NZ,NY,NX)))*CPPOLN/(CPPOLN+CPKM)
C
C     NODULE SENESCENCE
C
C     RSNDL=excess maintenance respiration
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
C     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
C     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
C
      IF(RSNDL.GT.0.0.AND.WTNDB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RXNSNC=RSNDL
      RXNSNN=RXNSNC*WTNDBN(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX)
      RXNSNP=RXNSNC*WTNDBP(NB,NZ,NY,NX)/WTNDB(NB,NZ,NY,NX)
      RDNSNC=RXNSNC*(1.0-RCCC)
      RDNSNN=RXNSNN*(1.0-RCCC)*(1.0-RCCN)
      RDNSNP=RXNSNP*(1.0-RCCC)*(1.0-RCCP)
      RCNSNC=RXNSNC-RDNSNC
      RCNSNN=RXNSNN-RDNSNN
      RCNSNP=RXNSNP-RDNSNP
      ELSE
      RXNSNC=0.0
      RXNSNN=0.0
      RXNSNP=0.0
      RDNSNC=0.0
      RDNSNN=0.0
      RDNSNP=0.0
      RCNSNC=0.0
      RCNSNN=0.0
      RCNSNP=0.0
      ENDIF
C
C     TOTAL NODULE RESPIRATION
C
C     RCO2T=total C respiration
C     RMNDL=bacterial maintenance respiration
C     RCNDL=respiration from non-structural C
C     RGNDG=bacterial respiration for growth and N2 fixation
C     RCNSNC=bacterial C senescence to recycling
C     TCO2T,TCO2A=total,above-ground PFT respiration
C     CNET=PFT net CO2 fixation
C     RECO=ecosystem respiration
C     TRAU=total autotrophic respiration
C
      RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNC
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)-RCO2T
      TCO2A(NZ,NY,NX)=TCO2A(NZ,NY,NX)-RCO2T
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-RCO2T
      HCNET(NZ,NY,NX)=HCNET(NZ,NY,NX)-RCO2T
      RECO(NY,NX)=RECO(NY,NX)-RCO2T
      TRAU(NY,NX)=TRAU(NY,NX)-RCO2T
C
C     NODULE LITTERFALL CAUSED BY REMOBILIZATION
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and
C        senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
C     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
C
      DO 6470 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+CFOPC(1,M,NZ,NY,NX)
     2*(RDNDLC+RDNSNC)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+CFOPN(1,M,NZ,NY,NX)
     2*(RDNDLN+RDNSNN)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+CFOPP(1,M,NZ,NY,NX)
     2*(RDNDLP+RDNSNP)
C     WRITE(*,2119)'NODLIT',I,J,NFZ,NZ,NB,M,INTYP(NZ,NY,NX)
C    2,CSNC(M,1,0,NZ,NY,NX),CFOPC(1,M,NZ,NY,NX)
C    2,RDNDLC,RDNSNC,RXNDLC,(1.0-RCCC),SPNDX,WTNDB(NB,NZ,NY,NX)
2119  FORMAT(A8,7I4,12E12.4)
6470  CONTINUE
C
C     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
C
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     RMNDL=bacterial maintenance respiration
C     RCNDL=respiration from non-structural C
C     RGN2F=respiration for N2 fixation
C     CGNDL=total non-structural C used in bacterial growth
C        and growth respiration
C     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
C     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
C     ZADDN,PADDN=nonstructural N,P used in growth
C     RUPNFB=branch N2 fixation
C
      CPOLNB(NB,NZ,NY,NX)=CPOLNB(NB,NZ,NY,NX)-AMIN1(RMNDL,RCNDL)
     2-RGN2F-CGNDL+RCNDLC
      ZPOLNB(NB,NZ,NY,NX)=ZPOLNB(NB,NZ,NY,NX)-ZADDN+RCNDLN+RCNSNN
     2+RUPNFB
      PPOLNB(NB,NZ,NY,NX)=PPOLNB(NB,NZ,NY,NX)-PADDN+RCNDLP+RCNSNP
C
C     UPDATE STATE VARIABLES FOR NODULE C, N, P
C
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     GRNDG=bacterial growth
C     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
C     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
C     ZADDN,PADDN=nonstructural N,P used in growth
C
      WTNDB(NB,NZ,NY,NX)=WTNDB(NB,NZ,NY,NX)+GRNDG-RXNDLC-RXNSNC
      WTNDBN(NB,NZ,NY,NX)=WTNDBN(NB,NZ,NY,NX)+ZADDN-RXNDLN-RXNSNN
      WTNDBP(NB,NZ,NY,NX)=WTNDBP(NB,NZ,NY,NX)+PADDN-RXNDLP-RXNSNP
C     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,2121)'NODGB',I,J,NFZ,NZ,NB,INTYP(NZ,NY,NX)
C    2,CNDLB,WTNDB(NB,NZ,NY,NX),ARLFB(NB,NZ,NY,NX)
C    2,GRNDG,RXNDLC,RXNSNC,CNDLB/CNDLI(INTYP(NZ,NY,NX))
C    2,RCNDL,RMNDL,RGNDL,RGN2P,RCO2T,RXNDLC
C    2,RGN2P,RGN2F,CGNDL,RSNDL,RGNDG,RCNSNC
C    3,ZADDN,PADDN,RCCC,RCCN,SPNDX,SPNDL
C    2,TFN3(NZ,NY,NX),WFNSG,CNDLB,CNDLI(INTYP(NZ,NY,NX))
C    8,RCCP,RDNDLC,RDNDLN,RDNDLP,RDNDLX,WTLSB(NB,NZ,NY,NX)
C    3,WTNDB(NB,NZ,NY,NX),WTNDBN(NB,NZ,NY,NX),WTNDBP(NB,NZ,NY,NX)
C    4,CPOLNB(NB,NZ,NY,NX),ZPOLNB(NB,NZ,NY,NX),PPOLNB(NB,NZ,NY,NX)
C    5,CCPOLN,CZPOLN,CPPOLN,ZCPOLI,ZPPOLI
C    6,(1.0+AMAX1(ZCPOLI/ZCKI,ZPPOLI/ZPKI))
C    6,TFN3(NZ,NY,NX),FCNPF,WFNSG,CPOOLNX,VMXOX
2121  FORMAT(A8,6I4,60F16.8)
C     ENDIF
C
C     TRANSFER NON-STRUCTURAL C,N,P BETWEEN BRANCH AND NODULES
C     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
C
C     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     WTLSB=leaf+petiole C mass
C     IDAY(10,=physiological maturity date
C     WTNDB=bacterial C mass
C     WTNDI=initial bacterial mass at infection
C     FXRN=rate constant for plant-bacteria nonstructural C,N,P
C        exchange
C     XFRC,XFRN,XFRC=nonstructural canopy-bacteria C,N,P transfer
C
      IF(CPOOL(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WTLSB(NB,NZ,NY,NX).GT.ZEROP2(NZ,NY,NX)
     3.AND.(ISTYP(NZ,NY,NX).NE.0
     4.OR.IDAY(10,NB,NZ,NY,NX).EQ.0))THEN
      WTLSB1=WTLSB(NB,NZ,NY,NX)
      WTNDB1=AMIN1(WTLSB(NB,NZ,NY,NX)
     2,AMAX1(WTNDI*AREA(3,NU(NY,NX),NY,NX),WTNDB(NB,NZ,NY,NX)))
      WTLSBT=WTLSB1+WTNDB1
      IF(WTLSBT.GT.ZEROP(NZ,NY,NX))THEN
      CPOOLD=(CPOOL(NB,NZ,NY,NX)*WTNDB1
     2-CPOLNB(NB,NZ,NY,NX)*WTLSB1)/WTLSBT
      XFRC=FXRN(INTYP(NZ,NY,NX))*CPOOLD
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)-XFRC
      CPOLNB(NB,NZ,NY,NX)=CPOLNB(NB,NZ,NY,NX)+XFRC
      CPOOLT=CPOOL(NB,NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX)
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLD=(ZPOOL(NB,NZ,NY,NX)*CPOLNB(NB,NZ,NY,NX)
     2-ZPOLNB(NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX))/CPOOLT
      XFRN=FXRN(INTYP(NZ,NY,NX))*ZPOOLD
      PPOOLD=(PPOOL(NB,NZ,NY,NX)*CPOLNB(NB,NZ,NY,NX)
     2-PPOLNB(NB,NZ,NY,NX)*CPOOL(NB,NZ,NY,NX))/CPOOLT
      XFRP=FXRN(INTYP(NZ,NY,NX))*PPOOLD
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-XFRP
      ZPOLNB(NB,NZ,NY,NX)=ZPOLNB(NB,NZ,NY,NX)+XFRN
      PPOLNB(NB,NZ,NY,NX)=PPOLNB(NB,NZ,NY,NX)+XFRP
C     IF((I/30)*30.EQ.I.AND.J.EQ.12)THEN
C     WRITE(*,2120)'NODEX',I,J,NZ,NB,IFLGA(NB,NZ,NY,NX)
C    2,XFRC,XFRN,XFRP
C    3,WTLSB(NB,NZ,NY,NX),WTNDB(NB,NZ,NY,NX),CPOOLT
C    4,CPOLNB(NB,NZ,NY,NX),ZPOLNB(NB,NZ,NY,NX),PPOLNB(NB,NZ,NY,NX)
C    4,CPOOL(NB,NZ,NY,NX),ZPOOL(NB,NZ,NY,NX),PPOOL(NB,NZ,NY,NX)
C    5,WTLSB1,WTNDB1,WTLSBT
2120  FORMAT(A8,5I4,40E12.4)
C     ENDIF
C     WRITE(*,2121)'NODBAL',I,J,NZ,NB,CPOLNB(NB,NZ,NY,NX)
C    2,WTNDB(NB,NZ,NY,NX),CPOLNB(NB,NZ,NY,NX)+WTNDB(NB,NZ,NY,NX)
C    3,RMNDL,RCNDL,RGN2F,CGNDL,RCNDLC,GRNDG,RXNDLC,RXNSNC,RCO2T
C    4,RGNDG,RGNDL,RCNSNC
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
105   CONTINUE
C
C     ROOT GROWTH
C
      NIX(NZ,NY,NX)=NG(NZ,NY,NX)
C
C     FOR ROOTS (N=1) AND MYCORRHIZAE (N=2) IN EACH SOIL LAYER
C
      DO 4995 N=1,MY(NZ,NY,NX)
      DO 4990 L=NU(NY,NX),NI(NZ,NY,NX)
C
C     RESPIRATION FROM NUTRIENT UPTAKE CALCULATED IN 'UPTAKE':
C     ACTUAL, O2-UNLIMITED AND C-UNLIMITED
C
C     VOLX=soil layer volume excluding macropore, rocks
C     CUPRL=C respiration for nutrient uptake
C     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
C     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of
C        NH4,NO3
C     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of
C        H2PO4,HPO4
C     RUONH4,RUONHB,RUON03,RUONOB=uptake from non-band,band of
C        NH4,NO3 unlimited by O2
C     RUOH2P,RUOH2B,RUOH1P,RUOH1B=uptake from non-band,band of
C        H2PO4,HPO4 unlimited by O2
C     RUCNH4,RUCNHB,RUCN03,RUCNOB=uptake from non-band,band of
C        NH4,NO3 unlimited by nonstructural C
C     RUCH2P,RUCH2B,RUCH1P,RUCH1B=uptake from non-band,band of
C        H2PO4,HPO4 unlimited by nonstructural C
C
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      CUPRL=0.86*(RUPNH4(N,L,NZ,NY,NX)+RUPNHB(N,L,NZ,NY,NX)
     2+RUPNO3(N,L,NZ,NY,NX)+RUPNOB(N,L,NZ,NY,NX)+RUPH2P(N,L,NZ,NY,NX)
     3+RUPH2B(N,L,NZ,NY,NX)+RUPH1P(N,L,NZ,NY,NX)+RUPH1B(N,L,NZ,NY,NX))
      CUPRO=0.86*(RUONH4(N,L,NZ,NY,NX)+RUONHB(N,L,NZ,NY,NX)
     2+RUONO3(N,L,NZ,NY,NX)+RUONOB(N,L,NZ,NY,NX)+RUOH2P(N,L,NZ,NY,NX)
     3+RUOH2B(N,L,NZ,NY,NX)+RUOH1P(N,L,NZ,NY,NX)+RUOH1B(N,L,NZ,NY,NX))
      CUPRC=0.86*(RUCNH4(N,L,NZ,NY,NX)+RUCNHB(N,L,NZ,NY,NX)
     2+RUCNO3(N,L,NZ,NY,NX)+RUCNOB(N,L,NZ,NY,NX)+RUCH2P(N,L,NZ,NY,NX)
     3+RUCH2B(N,L,NZ,NY,NX)+RUCH1P(N,L,NZ,NY,NX)+RUCH1B(N,L,NZ,NY,NX))
C
C     ACCUMULATE RESPIRATION IN FLUX ARRAYS
C
C     RCO2A=total root respiration
C     RCO2M,RCO2N=RCO2A unlimited by O2,nonstructural C
C     CUPRL=C respiration for nutrient uptake
C     CUPRO,CUPRC=CUPRL unlimited by O2,root nonstructural C
C     CPOOLR=non-structural C mass in root
C
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)+CUPRO
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)+CUPRC
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)-CUPRL
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-CUPRL
C     WRITE(*,8876)'CUPRL',I,J,NFZ,NZ,L,N
C    2,CPOOLR(N,L,NZ,NY,NX),CUPRL,RUPNH4(N,L,NZ,NY,NX)
C    3,RUPNHB(N,L,NZ,NY,NX)
C    2,RUPNO3(N,L,NZ,NY,NX),RUPNOB(N,L,NZ,NY,NX),RUPH2P(N,L,NZ,NY,NX)
C    3,RUPH2B(N,L,NZ,NY,NX),RUPH1P(N,L,NZ,NY,NX),RUPH1B(N,L,NZ,NY,NX)
8876  FORMAT(A8,6I4,20E12.4)
C
C     EXUDATION AND UPTAKE OF C, N AND P TO/FROM SOIL AND ROOT
C     OR MYCORRHIZAL NON-STRUCTURAL C,N,P POOLS
C
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange
C        -ve=exudn,+ve=uptake
C     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of
C        NH4,NO3
C     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of
C        H2PO4,HPO4
C
      DO 195 K=0,4
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)+RDFOMC(N,K,L,NZ,NY,NX)
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)+RDFOMN(N,K,L,NZ,NY,NX)
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)+RDFOMP(N,K,L,NZ,NY,NX)
C     WRITE(*,8877)'RDFOMC',I,J,NFZ,NZ,L,N,K
C    2,CPOOLR(N,L,NZ,NY,NX),RDFOMC(N,K,L,NZ,NY,NX)
8877  FORMAT(A8,7I4,12E12.4)
195   CONTINUE
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)
     2+RUPNH4(N,L,NZ,NY,NX)+RUPNHB(N,L,NZ,NY,NX)
     3+RUPNO3(N,L,NZ,NY,NX)+RUPNOB(N,L,NZ,NY,NX)
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)
     2+RUPH2P(N,L,NZ,NY,NX)+RUPH2B(N,L,NZ,NY,NX)
     3+RUPH1P(N,L,NZ,NY,NX)+RUPH1B(N,L,NZ,NY,NX)
C     IF(NZ.EQ.2.AND.L.EQ.6)THEN
C     WRITE(*,9881)'RUPNH4',I,J,NFZ,NZ,L,N,CPOOLR(N,L,NZ,NY,NX)
C    2,ZPOOLR(N,L,NZ,NY,NX),PPOOLR(N,L,NZ,NY,NX),CUPRL
C    2,RDFOMC(N,K,L,NZ,NY,NX),RDFOMN(N,K,L,NZ,NY,NX)
C    3,RDFOMP(N,K,L,NZ,NY,NX)
C    2,RUPNH4(N,L,NZ,NY,NX),RUPNHB(N,L,NZ,NY,NX),RUPNO3(N,L,NZ,NY,NX)
C    2,RUPNOB(N,L,NZ,NY,NX),RUPH2P(N,L,NZ,NY,NX),RUPH2B(N,L,NZ,NY,NX)
C    3,RUPH1P(N,L,NZ,NY,NX),RUPH1B(N,L,NZ,NY,NX),WFR(N,L,NZ,NY,NX)
9881  FORMAT(A8,6I4,30E12.4)
C     ENDIF
C
C     GROWTH OF EACH ROOT AXIS
C
      DO 4985 NR=1,NRT(NZ,NY,NX)
C
C     PRIMARY ROOT SINK STRENGTH FROM ROOT RADIUS AND ROOT DEPTH
C
C     RTDP1=primary root depth from soil surface
C     RTDPP=primary root depth from canopy
C     CDPTHZ=depth from soil surface to layer bottom
C     HTSTZ=canopy height for water uptake
C     RTSK=relative primary root sink strength
C     RTSK1=primary root sink strength
C     XRTN1=number of primary root axes
C     RRAD1,RRAD2=primary,secondary root radius
C     RTNT,RLNT=total root sink strength
C
      IF(N.EQ.1)THEN
      IF(RTDP1(N,NR,NZ,NY,NX).GT.CDPTHZ(L-1,NY,NX))THEN
      IF(RTDP1(N,NR,NZ,NY,NX).LE.CDPTHZ(L,NY,NX))THEN
      RTDPP=RTDP1(N,NR,NZ,NY,NX)+HTSTZ(NZ,NY,NX)
      RTSK1(N,L,NR)=RTSK(IGTYP(NZ,NY,NX))*XRTN1
     2*RRAD1(N,L,NZ,NY,NX)**2/RTDPP
      RTNT(N)=RTNT(N)+RTSK1(N,L,NR)
      RLNT(N,L)=RLNT(N,L)+RTSK1(N,L,NR)
      ENDIF
      ENDIF
      ENDIF
C
C     SECONDARY ROOT SINK STRENGTH FROM ROOT RADIUS, ROOT AXIS NUMBER,
C     AND ROOT LENGTH IN SERIES WITH PRIMARY ROOT SINK STRENGTH
C
C     RTDPL=depth of primary root axis in layer
C     RTDP1=primary root depth from soil surface
C     CDPTHZ=depth from soil surface to layer bottom
C     RTDPX=distance behind growing point for secondary roots
C     DLYR=layer thickness
C     SDPTH=seeding depth
C     HTCTL=hypocotyledon height
C     HTSTZ=canopy height for water uptake
C     RTDPS=secondary root depth from canopy
C     RTSKP,RTSKS=primary,secondary root sink strength
C     RTN2=number of secondary root axes
C     RTSK2=total secondary root sink strength
C     RTLGA=average secondary root length
C     RTNT,RLNT=total root sink strength
C
      IF(N.EQ.1)THEN
      RTDPL(NR,L)=AMAX1(0.0,RTDP1(1,NR,NZ,NY,NX)-CDPTHZ(L-1,NY,NX)
     2-RTDPX)
      RTDPL(NR,L)=AMAX1(0.0,AMIN1(DLYR(3,L,NY,NX),RTDPL(NR,L))
     2-AMAX1(0.0,SDPTH(NZ,NY,NX)-CDPTHZ(L-1,NY,NX)-HTCTL(NZ,NY,NX)))
      RTDPS=AMAX1(SDPTH(NZ,NY,NX),CDPTHZ(L-1,NY,NX))
     2+0.5*RTDPL(NR,L)+HTSTZ(NZ,NY,NX)
      IF(RTDPS.GT.ZERO)THEN
      RTSKP=XRTN1*RRAD1(N,L,NZ,NY,NX)**2/RTDPS
      RTSKS=RTN2(N,L,NR,NZ,NY,NX)*RRAD2(N,L,NZ,NY,NX)**2
     2/RTLGA(N,L,NZ,NY,NX)
      IF(RTSKP+RTSKS.GT.ZEROP(NZ,NY,NX))THEN
      RTSK2(N,L,NR)=RTSKP*RTSKS/(RTSKP+RTSKS)
      ELSE
      RTSK2(N,L,NR)=0.0
      ENDIF
      ELSE
      RTSK2(N,L,NR)=0.0
      ENDIF
      ELSE
      RTSK2(N,L,NR)=RTN2(N,L,NR,NZ,NY,NX)*RRAD2(N,L,NZ,NY,NX)**2
     2/RTLGA(N,L,NZ,NY,NX)
      ENDIF
      RTNT(N)=RTNT(N)+RTSK2(N,L,NR)
      RLNT(N,L)=RLNT(N,L)+RTSK2(N,L,NR)
C     IF(IYRC.EQ.2000.AND.I.LE.160)THEN
C     WRITE(*,3341)'SINK',I,J,NX,NY,NZ,L,NR,N,RTDP1(N,NR,NZ,NY,NX)
C    2,HTCTL(NZ,NY,NX),RTSK1(N,L,NR),RTSK2(N,L,NR),RLNT(N,L),RTNT(N)
C    3,XRTN1,PP(NZ,NY,NX),RRAD1(N,L,NZ,NY,NX),RTDPS,RTDPP
C    4,RTDPL(NR,L),RTN2(N,L,NR,NZ,NY,NX),RRAD2(N,L,NZ,NY,NX)
C    2,RTLGA(N,L,NZ,NY,NX),CDPTHZ(L-1,NY,NX),CDPTHZ(L,NY,NX)
3341  FORMAT(A8,8I4,30E12.4)
C     ENDIF
4985  CONTINUE
      ENDIF
4990  CONTINUE
4995  CONTINUE
C
C     ADAPT ROOT POROSITY TO O2 STRESS
C
C     PORT,PORTI=current,initial root porosity from PFT file
C     OSTR=O2 root uptake:demand
C
      DO 5010 N=1,MY(NZ,NY,NX)
      PORT(N,NZ,NY,NX)=AMIN1(0.75,PORT(N,NZ,NY,NX)
     2+(1.0E-02*PORTI(N,NZ,NY,NX)*(1.0-OSTR(NZ,NY,NX))
     3-1.0E-02*(PORT(N,NZ,NY,NX)-PORTI(N,NZ,NY,NX)))*XNFH)
C     IF((I/10)*10.EQ.I.AND.J.EQ.12)THEN
C     WRITE(*,2127)'PORT',I,J,NZ,N,PORT(N,NZ,NY,NX)
C    2,PORTI(N,NZ,NY,NX),OSTR(NZ,NY,NX)
2127  FORMAT(A8,4I4,12E12.4)
C     ENDIF
C
C     RESPIRATION AND GROWTH OF ROOT, MYCORRHIZAE IN EACH LAYER
C
      DO 5000 L=NU(NY,NX),NI(NZ,NY,NX)
C
C     IDENTIFY NEXT LOWER ROOT LAYER L1
C
C     VOLX=soil layer volume excluding macropore, rocks (0=pond)
C
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
C     WRITE(*,4994)'5004',I,J,NZ,N,L,NI(NZ,NY,NX)
C    2,NL(NY,NX),VOLX(L,NY,NX),CDPTHZ(L-1,NY,NX)
      DO 5003 LZ=L+1,NL(NY,NX)
C     WRITE(*,4994)'5003',I,J,NZ,N,L,LZ
C    2,LZ,VOLX(L,NY,NX),CDPTHZ(LZ,NY,NX)
      IF(VOLX(LZ,NY,NX).GT.ZEROS2(NY,NX)
     2.OR.LZ.EQ.NL(NY,NX))THEN
      L1=LZ
      GO TO 5004
      ENDIF
5003  CONTINUE
5004  CONTINUE
C     WRITE(*,4994)'5005',I,J,NZ,N,L,LZ
C    2,L1,VOLX(L,NY,NX),CDPTHZ(L1,NY,NX)
4994  FORMAT(A8,7I4,12E12.4)
C
C     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
C     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
C
C     RSCS,RSCS2=soil resistance to secondary root penetration (MPa)
C     RRAD2=secondary root radius
C     WFNRT=water function for root extension
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     WFNRG,WFNRR=growth,respiration function of root water potential
C     PSIRT,PSIRG=root total,turgor water potential
C     DMRT=root growth yield
C
      RSCS2=RSCS(L,NY,NX)*RRAD2(N,L,NZ,NY,NX)/1.0E-03
      WFNRT=AMIN1(1.0,AMAX1(0.0,PSIRG(N,L,NZ,NY,NX)-PSILM-RSCS2))
      IF(IGTYP(NZ,NY,NX).EQ.0)THEN
      WFNRG(N,L)=EXP(0.05*PSIRT(N,L,NZ,NY,NX))
      WFNRR(N,L)=WFNRG(N,L)**0.25
      ELSE
      WFNRG(N,L)=EXP(0.10*PSIRT(N,L,NZ,NY,NX))
      WFNRR(N,L)=WFNRG(N,L)**0.25
      ENDIF
      DMRTD=1.0-DMRT(NZ,NY,NX)
      RTLGL=0.0
      RTLGZ=0.0
      WTRTX=0.0
      WTRTZ=0.0
C
C     FOR EACH ROOT AXIS
C
      DO 5050 NR=1,NRT(NZ,NY,NX)
C
C     SECONDARY ROOT EXTENSION
C
      IF(L.LE.NINR(NR,NZ,NY,NX).AND.NRX(N,NR).EQ.0)THEN
C
C     FRACTION OF SECONDARY ROOT SINK IN SOIL LAYER ATTRIBUTED
C     TO CURRENT AXIS
C
C     RTSK2=total secondary root sink strength
C     RLNT=total root sink strength
C     FRTN=fraction of secondary root sink strength in axis
C
      IF(RLNT(N,L).GT.ZEROP(NZ,NY,NX))THEN
      FRTN=RTSK2(N,L,NR)/RLNT(N,L)
      ELSE
      FRTN=1.0
      ENDIF
C
C     N,P CONSTRAINT ON SECONDARY ROOT RESPIRATION FROM
C     NON-STRUCTURAL C:N:P
C
C     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
C     CNPG=N,P constraint on growth respiration
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C
      IF(CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      CNPG=AMIN1(CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX)
     2+CCPOLR(N,L,NZ,NY,NX)*CNKI),CPPOLR(N,L,NZ,NY,NX)
     3/(CPPOLR(N,L,NZ,NY,NX)+CCPOLR(N,L,NZ,NY,NX)*CPKI))
      ELSE
      CNPG=1.0
      ENDIF
C
C     O2-UNLIMITED SECONDARY ROOT RESPIRATION FROM NON-STRUCTURAL C
C     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
C
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     WFNRG=growth function of root water potential
C     RCO2RM=respiration from non-structural C unlimited by O2
C     VMXC=rate constant for nonstructural C oxidation in respiration
C     FRTN=fraction of secondary root sink strength in axis
C     CPOOL=non-structural C mass
C     TFN4=temperature function for root growth
C     CNPG=N,P constraint on respiration
C     FDBKX=termination feedback inhibition on C3 CO2
C     WFNRG=growth function of root water potential
C     RMNCR=root maintenance respiration
C     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
C     WTRT2N=secondary root N mass
C     TFN6=temperature function for root maintenance respiration
C     XNFH=time step from �wthr.f�
C
      IF(ISTYP(NZ,NY,NX).NE.0
     2.OR.IDAY(10,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0)THEN
      CPOOLRX=CPOOLR(N,L,NZ,NY,NX)
      RCO2RM=AMAX1(0.0,VMXC*FRTN*CPOOLR(N,L,NZ,NY,NX)
     2*TFN4(L,NZ,NY,NX))*CNPG*FDBKX(NB1(NZ,NY,NX),NZ,NY,NX)
     3*WFNRG(N,L)*XNFH
      RMNCR=AMAX1(0.0,RMPLT*WTRT2N(N,L,NR,NZ,NY,NX)*TFN6(L)*XNFH)
      IF(IGTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.2)THEN
      RMNCR=RMNCR*WFNRG(N,L)
      ELSE
      RMNCR=RMNCR*WFNRR(N,L)
      ENDIF
      ELSE
      RCO2RM=0.0
      RMNCR=0.0
      ENDIF
C
C     O2-LIMITED SECONDARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
C
C     RCO2R=respiration from non-structural C limited by O2
C     WFR=constraint by O2 consumption on all root processes
C     RCO2XM,RCO2X=diff between C respiration unlimited,limited by O2
C        and maintenance respiration
C     RCO2YM,RCO2Y=growth respiration unlimited,limited by O2 and
C        unlimited by N,P
C
      RCO2R=RCO2RM*WFR(N,L,NZ,NY,NX)
      RCO2XM=RCO2RM-RMNCR
      RCO2X=RCO2R-RMNCR
      RCO2YM=AMAX1(0.0,RCO2XM)
      RCO2Y=AMAX1(0.0,RCO2X)
C
C     SECONDARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
C     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
C
C     FRTN=fraction of secondary root sink strength in axis
C     ZPOOLR,PPOOLR=non-structural N,P mass in root
C     CNRTS,CPRTS=N,P root growth yield
C     FNP=growth respiration limited by non-structural N,P
C     RCO2GM,RCO2G=growth respiration limited by N,P unlimited,limited
C        by O2
C
      DMRTR=DMRTD*FRTN
      ZPOOLB=AMAX1(0.0,ZPOOLR(N,L,NZ,NY,NX))
      PPOOLB=AMAX1(0.0,PPOOLR(N,L,NZ,NY,NX))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ,NY,NX)
     2,PPOOLB*DMRTR/CPRTS(NZ,NY,NX))
      IF(RCO2YM.GT.0.0)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
      ELSE
      RCO2GM=0.0
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
      RCO2G=AMIN1(RCO2Y,FNP*WFR(N,L,NZ,NY,NX))
      ELSE
      RCO2G=0.0
      ENDIF
C
C     TOTAL NON-STRUCTURAL C,N,P USED IN SECONDARY ROOT GROWTH
C     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD ENTERED IN
C     'READQ'
C
C     CGRORM,CGROR=total non-structural C used in growth and
C        respiration unlimited,limited by O2
C     RCO2GM,RCO2G=growth respiration limited by N,P unlimited,limited
C        by O2
C     DMRTD=root C respiration vs nonstructural C consumption
C     GRTWGM,GRTWTG=root C growth unlimited,limited by O2
C     DMRT=root growth yield
C     ZADD2M,ZADD2,PADD2=nonstructural N,P unlimited,limited
C        by O2 used in growth
C     CNRDM,CNRDA=respiration for N assimilation unlimited,limited
C        by O2
C
      CGRORM=RCO2GM/DMRTD
      CGROR=RCO2G/DMRTD
      GRTWGM=CGRORM*DMRT(NZ,NY,NX)
      GRTWTG=CGROR*DMRT(NZ,NY,NX)
      ZADD2M=AMAX1(0.0,GRTWGM*CNRTW)
      ZADD2=AMAX1(0.0,AMIN1(FRTN*ZPOOLR(N,L,NZ,NY,NX),GRTWTG*CNRTW))
      PADD2=AMAX1(0.0,AMIN1(FRTN*PPOOLR(N,L,NZ,NY,NX),GRTWTG*CPRTW))
      CNRDM=AMAX1(0.0,1.70*ZADD2M)
      CNRDA=AMAX1(0.0,1.70*ZADD2)
C
C     SECONDARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
C     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
C     SECONDARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT
C     LITTERFALL
C
C     IDAY(1,=emergence date
C     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     RCCZR,RCCYR=minimum,maximum fractions for root C recycling
C     RCCXR,RCCQR=maximum fractions for root N,P recycling
C
      IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0
     2.AND.CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      CCC=AMAX1(0.0,AMIN1(1.0
     1,CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX)
     2+CCPOLR(N,L,NZ,NY,NX)*CNKI)
     3,CPPOLR(N,L,NZ,NY,NX)/(CPPOLR(N,L,NZ,NY,NX)
     4+CCPOLR(N,L,NZ,NY,NX)*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX)
     2+CZPOLR(N,L,NZ,NY,NX)/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)/CPKI)))
      ELSE
      CCC=0.0
      CNC=0.0
      CPC=0.0
      ENDIF
      RCCC=RCCZR(IBTYP(NZ,NY,NX))+CCC*RCCYR(IBTYP(NZ,NY,NX))
      RCCN=CNC*RCCXR(IBTYP(NZ,NY,NX))
      RCCP=CPC*RCCQR(IBTYP(NZ,NY,NX))
C
C     RECOVERY OF REMOBILIZABLE N,P FROM SECONDARY ROOT DURING
C     REMOBILIZATION DEPENDS ON ROOT NON-STRUCTURAL C:N:P
C
C     RCO2XM,RCO2X=diff between C respiration unlimited,limited
C        by O2 and maintenance respiration
C     SNCRM,SNCR=excess maintenance respiration unlimited,limited
C        by O2
C     SCNZ=phenologically driven root senescence
C     FXFR=rate constant for root-storage nonstructural C,N,P exchange
C          from IBTYP in PFT file
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
C     WFR=constraint by O2 consumption on all root processes
C     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     FSNC2=fraction of secondary root C to be remobilized
C
      IF(-RCO2XM.GT.0.0)THEN
      IF(-RCO2XM.LT.WTRT2(N,L,NR,NZ,NY,NX)*RCCC)THEN
      SNCRM=-RCO2XM
      ELSE
      SNCRM=AMAX1(0.0,WTRT2(N,L,NR,NZ,NY,NX)*RCCC)
      ENDIF
      ELSE
      SNCRM=0.0
      ENDIF
      IF(-RCO2X.GT.0.0)THEN
      IF(-RCO2X.LT.WTRT2(N,L,NR,NZ,NY,NX)*RCCC)THEN
      SNCR=-RCO2X
      ELSE
      SNCR=AMAX1(0.0,WTRT2(N,L,NR,NZ,NY,NX)*RCCC)
     2*WFR(N,L,NZ,NY,NX)
      ENDIF
      ELSE
      SNCR=0.0
      ENDIF
C     IF(IFLGZR.EQ.1.AND.IFLGYR.EQ.1.AND.ISTYP(NZ,NY,NX).NE.0)THEN
      IF(IFLGZR.EQ.1.AND.IFLGYR.EQ.1)THEN
      SNCZ=FXFR(IBTYP(NZ,NY,NX))*WTRT2(N,L,NR,NZ,NY,NX)
     2*AMIN1(1.0,FLGZR/FLGZX)*XNFH
      SNCR=SNCR+SNCZ
      IF(J.EQ.24)THEN
C     WRITE(*,9776)'IFLGR',I,J,NFZ,NZ,NR,L,N
C    2,IFLGZR,IFLGYR,ISTYP(NZ,NY,NX),FLGZR,FLGZX
C    3,WTRT2(N,L,NR,NZ,NY,NX),SNCZ,SNCR,XNFH
9776  FORMAT(A8,10I4,20E12.4)
      ENDIF
      ENDIF
      IF(SNCR.GT.0.0.AND.WTRT2(N,L,NR,NZ,NY,NX)
     2.GT.ZEROP(NZ,NY,NX))THEN
      RCCR=RCCC*WTRT2(N,L,NR,NZ,NY,NX)
      RCZR=WTRT2N(N,L,NR,NZ,NY,NX)*(RCCN+(1.0-RCCN)*RCCC)
      RCPR=WTRT2P(N,L,NR,NZ,NY,NX)*(RCCP+(1.0-RCCP)*RCCC)
      IF(RCCR.GT.ZEROP(NZ,NY,NX))THEN
      FSNC2=AMAX1(0.0,AMIN1(1.0,SNCR/RCCR))
      ELSE
      FSNC2=1.0
      ENDIF
      ELSE
      RCCR=0.0
      RCZR=0.0
      RCPR=0.0
      FSNC2=0.0
      ENDIF
C
C     SECONDARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
C
C     CSNC,ZSNC,PSNC=literfall C,N,P
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FSNC2=fraction of secondary root C to be remobilized
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
C     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
C     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,
C        1=non-woody
C
      DO 6350 M=1,4
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*FSNC2*WTRT2(N,L,NR,NZ,NY,NX)*FWODR(0,NZ)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*FSNC2*WTRT2N(N,L,NR,NZ,NY,NX)*FWODRN(0,NZ)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*FSNC2*WTRT2P(N,L,NR,NZ,NY,NX)*FWODRP(0,NZ)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX)
     2*FSNC2*(WTRT2(N,L,NR,NZ,NY,NX)-RCCR)*FWODR(1,NZ)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX)
     2*FSNC2*(WTRT2N(N,L,NR,NZ,NY,NX)-RCZR)*FWODRN(1,NZ)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX)
     2*FSNC2*(WTRT2P(N,L,NR,NZ,NY,NX)-RCPR)*FWODRP(1,NZ)
C     IF(ZSNC(M,1,L,NZ,NY,NX).LT.0.0)THEN
C     WRITE(*,9777)'ZSNC',I,J,NFZ,NZ,NR,L,N,M
C    2,CSNC(M,1,L,NZ,NY,NX),CFOPC(4,M,NZ,NY,NX)
C    2,FSNC2,WTRT2(N,L,NR,NZ,NY,NX),RCCR,FWODR(1,NZ)
C    2,ZSNC(M,1,L,NZ,NY,NX),CFOPN(4,M,NZ,NY,NX)
C    2,FSNC2,WTRT2N(N,L,NR,NZ,NY,NX),RCZR,FWODRN(1,NZ)
C    3,FWODRN(0,NZ),FWODR(0,NZ),FWODR(1,NZ),CNRT(NZ,NY,NX),CNRTW
9777  FORMAT(A8,8I4,30E14.6)
C     ENDIF
6350  CONTINUE
C
C     CONSUMPTION OF NON-STRUCTURAL C,N,P BY SECONDARY ROOT
C
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     RMNCR=root maintenance respiration
C     RCO2R=respiration from non-structural C limited by O2
C     CGROR=total non-structural C used in growth and respiration
C        limited by O2
C     CNRDA=respiration for N assimilation unlimited,limited by O2
C     SNCR=excess maintenance respiration limited by O2
C     FSNC2=fraction of secondary root C to be remobilized
C     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
C     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
C
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-AMIN1(RMNCR,RCO2R)
     2-CGROR-CNRDA-SNCR+FSNC2*RCCR*FWODR(1,NZ)
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)
     2-ZADD2+FSNC2*RCZR*FWODRN(1,NZ)
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)
     2-PADD2+FSNC2*RCPR*FWODRP(1,NZ)
C
C     TOTAL SECONDARY ROOT RESPIRATION
C
C     RCO2TM,RCO2T=total C respiration unlimited,limited
C        by O2
C     RMNCR=root maintenance respiration
C     RCO2RM,RCO2R=respiration from non-structural C unlimited,limited
C        by O2
C     RCO2GM,RCO2G=growth respiration limited by N,P unlimited,limited
C        by O2
C     SNCRM,SNCR=excess maintenance respiration unlimited,limited
C        by O2
C     CNRDM,CNRDA=respiration for N assimilation unlimited,limited
C        by O2
C     RCO2A=total root respiration
C     RCO2M,RCO2N=RCO2A unlimited by O2,nonstructural C
C
      RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
      RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)+RCO2TM
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)+RCO2T
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)-RCO2T
C
C     SECONDARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
C
C     GRTLGL=secondary root length extension
C     GRTWTG=secondary root C growth limited by O2
C     RTLG2X=specific secondary root length from startq.f
C     WFNRT=water function for root extension
C     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
C     FSNC2=fraction of secondary root C to be remobilized
C     RTLG2=secondary root length
C     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
C     ZADD2,PADD2=nonstructural N,P ltd by O2 used in growth
C
      GRTLGL=GRTWTG*RTLG2X(N,NZ,NY,NX)*WFNRT
     2-FSNC2*RTLG2(N,L,NR,NZ,NY,NX)
      GRTWTL=GRTWTG-FSNC2*WTRT2(N,L,NR,NZ,NY,NX)
      GRTWTN=ZADD2-FSNC2*WTRT2N(N,L,NR,NZ,NY,NX)
      GRTWTP=PADD2-FSNC2*WTRT2P(N,L,NR,NZ,NY,NX)
C
C     UPDATE STATE VARIABLES FOR SECONDARY ROOT LENGTH, C, N, P
C     AND AXIS NUMBER
C
C     RTLG2=secondary root length
C     GRTLGL=secondary root length extension
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
C     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
C     WSRTL=total root protein C mass
C     CNWS,CPWS=protein:N,protein:P ratios from startq.f
C     RTFQ=root branching frequency from PFT file
C     RTN2,RTNL=number of secondary root axes
C
      RTLG2(N,L,NR,NZ,NY,NX)=RTLG2(N,L,NR,NZ,NY,NX)+GRTLGL
      WTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)+GRTWTL
      WTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)+GRTWTN
      WTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)+GRTWTP
      WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX)
     2+AMIN1(CNWS(NZ,NY,NX)*WTRT2N(N,L,NR,NZ,NY,NX)
     2,CPWS(NZ,NY,NX)*WTRT2P(N,L,NR,NZ,NY,NX))
      RTLGL=RTLGL+RTLG2(N,L,NR,NZ,NY,NX)
      WTRTX=WTRTX+WTRT2(N,L,NR,NZ,NY,NX)
      RTN2X=RTFQ(NZ,NY,NX)*XRTN1
      RTN2Y=RTFQ(NZ,NY,NX)*RTN2X
      RTN2(N,L,NR,NZ,NY,NX)=(RTN2X+RTN2Y)*DLYR(3,L,NY,NX)
      RTNL(N,L,NZ,NY,NX)=RTNL(N,L,NZ,NY,NX)+RTN2(N,L,NR,NZ,NY,NX)
C     IF(NZ.EQ.2.AND.L.EQ.6.AND.N.EQ.1)THEN
C     WRITE(*,9876)'RCO22',I,J,NFZ,NZ,NR,L,N,NINR(NR,NZ,NY,NX)
C    2,RCO2TM,RCO2T,RMNCR,RCO2RM,RCO2R,RCO2GM,RCO2G,CPOOLRX
C    3,RCO2XM,RCO2X,CGROR,SNCRM,SNCR,CNRDA,CPOOLR(N,L,NZ,NY,NX),FRTN
C    4,TFN4(L,NZ,NY,NX),CNPG,FDBKX(NB1(NZ,NY,NX),NZ,NY,NX),WFNRG(N,L)
C    5,TFN6(L),GRTWTG,GRTWTL,GRTLGL,RTLGL,RTLG2(N,L,NR,NZ,NY,NX)
C    5,WTRT2(N,L,NR,NZ,NY,NX),RTLG2X(N,NZ,NY,NX),WFNRT,FWODR(1,NZ)
C    4,FSNC2,RCO2M(N,L,NZ,NY,NX),RCO2A(N,L,NZ,NY,NX),WFR(N,L,NZ,NY,NX)
C    5,PSIRG(N,L,NZ,NY,NX),PSILM,RSCS1
C    8,ZPOOLR(N,L,NZ,NY,NX),PPOOLR(N,L,NZ,NY,NX)
C    9,ZADD2,FSNC2,RCZR,FWODRN(1,NZ)
C    9,FSNC2,RLNT(N,L),RTSK1(N,L,NR),RTSK2(N,L,NR)
C    4,RTN2X,RTN2Y,XRTN1
C    5,RTDPL(NR,L),RTDNP(N,L,NZ,NY,NX)
C    5,RTDP1(1,NR,NZ,NY,NX),CDPTHZ(L-1,NY,NX),DLYR(3,L,NY,NX)
C    6,SDPTH(NZ,NY,NX),HTCTL(NZ,NY,NX)
C    5,WFNRR(N,L),RSCS2,PSILM,PSIRG(N,L,NZ,NY,NX),PSIRT(N,L,NZ,NY,NX)
C    6,FNP,RTLGP(N,L,NZ,NY,NX),ZADD2,PADD2,CUPRO,CUPRL
C    7,RUPNH4(N,L,NZ,NY,NX),RUPNHB(N,L,NZ,NY,NX)
C    8,RUPNO3(N,L,NZ,NY,NX),RUPNOB(N,L,NZ,NY,NX)
C    9,RUPH2P(N,L,NZ,NY,NX),RUPH2B(N,L,NZ,NY,NX)
C    9,RUPH1P(N,L,NZ,NY,NX),RUPH1B(N,L,NZ,NY,NX)
C    6,RDFOMN(N,L,NZ,NY,NX),RDFOMP(N,L,NZ,NY,NX)
C    2,RTN1(N,L,NZ,NY,NX),RTN2(N,L,NR,NZ,NY,NX)
C    3,RTNL(N,L,NZ,NY,NX),DLYR(3,L,NY,NX)
9876  FORMAT(A8,8I4,100F16.8)
C     ENDIF
C
C     PRIMARY ROOT EXTENSION
C
C     BKDS=soil bulk density (0=pond)
C     RTDP1,RTDP1X=primary root depth from soil surface
C     CDPTHZ=depth from soil surface to layer bottom
C     ICHKL=flag for identifying layer with primary root tip
C     RTN1=number of primary root axes
C     XRTN1=multiplier for number of primary root axes
C
      IF(N.EQ.1)THEN
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      RTDP1X=RTDP1(N,NR,NZ,NY,NX)-CDPTHZ(0,NY,NX)
      ELSE
      RTDP1X=RTDP1(N,NR,NZ,NY,NX)
      ENDIF
      IF(RTDP1X.GT.CDPTHZ(L-1,NY,NX).AND.ICHK1(N,NR).EQ.0)THEN
      RTN1(N,L,NZ,NY,NX)=RTN1(N,L,NZ,NY,NX)+XRTN1
      IF(RTDP1X.LE.CDPTHZ(L,NY,NX).OR.L.EQ.NJ(NY,NX))THEN
      ICHK1(N,NR)=1
C     IF(J.EQ.24.AND.NZ.EQ.2)THEN
C     WRITE(*,9874)'RTDP1',I,J,NZ,NR,L,L-1,L1,N,NINR(NR,NZ,NY,NX)
C    2,ICHK1(N,NR),RTDP1(N,NR,NZ,NY,NX),RTDP1X,RTN1(N,L,NZ,NY,NX)
C    3,CDPTHZ(L-1,NY,NX),CDPTHZ(L,NY,NX)
9874  FORMAT(A8,10I4,12E12.4)
C     ENDIF
C
C     FRACTION OF PRIMARY ROOT SINK IN SOIL LAYER
C     ATTRIBUTED TO CURRENT AXIS
C
C     RTSK1=primary root sink strength
C     RLNT=total root sink strength
C     FRTN=fraction of primary root sink strength in axis
C
      IF(RLNT(N,L).GT.ZEROP(NZ,NY,NX))THEN
      FRTN=RTSK1(N,L,NR)/RLNT(N,L)
      ELSE
      FRTN=1.0
      ENDIF
C
C     WATER STRESS CONSTRAINT ON SECONDARY ROOT EXTENSION IMPOSED
C     BY ROOT TURGOR AND SOIL PENETRATION RESISTANCE
C
C     RSCS,RSCS1=soil resistance to primary root penetration (MPa)
C     RRAD1=primary root radius
C     PSIRG=root turgor potential
C     PSILM=minimum water potential for organ expansion,extension
C     WFNRT=water stress function for root extension
C     WFNRG,WFNRR=growth,respiration function of root water potential
C
      RSCS1=RSCS(L,NY,NX)*RRAD1(N,L,NZ,NY,NX)/1.0E-03
      WFNRT=AMIN1(1.0,AMAX1(0.0,PSIRG(N,L,NZ,NY,NX)-PSILM-RSCS1))
      IF(IGTYP(NZ,NY,NX).EQ.0)THEN
      WFNRG(N,L)=EXP(0.05*PSIRT(N,L,NZ,NY,NX))
      WFNRR(N,L)=WFNRG(N,L)**0.25
      ELSE
      WFNRG(N,L)=EXP(0.10*PSIRT(N,L,NZ,NY,NX))
      WFNRR(N,L)=WFNRG(N,L)**0.25
      ENDIF
C
C     N,P CONSTRAINT ON PRIMARY ROOT RESPIRATION FROM
C     NON-STRUCTURAL C:N:P
C
C     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
C     CNPG=N,P constraint on growth respiration
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C
      IF(CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      CNPG=AMIN1(CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX)
     2+CCPOLR(N,L,NZ,NY,NX)*CNKI),CPPOLR(N,L,NZ,NY,NX)
     3/(CPPOLR(N,L,NZ,NY,NX)+CCPOLR(N,L,NZ,NY,NX)*CPKI))
      ELSE
      CNPG=1.0
      ENDIF
C
C     O2-UNLIMITED PRIMARY ROOT RESPIRATION FROM ROOT NON-STRUCTURAL C
C     CONSTRAINED BY TEMPERATURE AND NON-STRUCTURAL C:N:P
C
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     WFNRG=growth function of root water potential
C     RCO2RM=respiration from non-structural C unlimited by O2
C     VMXC=rate constant for nonstructural C oxidation in respiration C     FRTN=fraction of primary root sink strength in axis
C     CPOOL=non-structural C mass
C     TFN4=temperature function for root growth
C     CNPG=N,P constraint on respiration
C     FDBKX=termination feedback inhibition on C3 CO2
C     WFNRG=growth function of root water potential
C     RMNCR=root maintenance respiration
C     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
C     WTRT1N=primary root N mass
C     TFN6=temperature function for root maintenance respiration
C     XNFH=time step from �wthr.f�
C
      IF(ISTYP(NZ,NY,NX).NE.0
     2.OR.IDAY(10,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0)THEN
      RCO2RM=AMAX1(0.0,VMXC*FRTN*CPOOLR(N,L,NZ,NY,NX)
     2*TFN4(L,NZ,NY,NX))*CNPG*FDBKX(NB1(NZ,NY,NX),NZ,NY,NX)
     3*WFNRG(N,L)*XNFH
      IF(RTDP1X.GE.CDPTHZ(NJ(NY,NX),NY,NX))THEN
      RCO2RM=AMIN1(RMNCR,RCO2RM)
      ENDIF
      RMNCR=AMAX1(0.0,RMPLT*RTWT1N(N,NR,NZ,NY,NX)*TFN6(L)*XNFH)
      IF(IGTYP(NZ,NY,NX).EQ.0.OR.IWTYP(NZ,NY,NX).EQ.2)THEN
      RMNCR=RMNCR*WFNRG(N,L)
      ELSE
      RMNCR=RMNCR*WFNRR(N,L)
      ENDIF
      ELSE
      RCO2RM=0.0
      RMNCR=0.0
      ENDIF
C
C     O2-LIMITED PRIMARY ROOT RESPIRATION FROM 'WFR' IN 'UPTAKE'
C
C     RCO2R=respiration from non-structural C limited by O2
C     WFR=constraint by O2 consumption on all root processes
C     RCO2XM,RCO2X=diff between C respiration unlimited,limited by O2
C        and maintenance respiration
C     RCO2YM,RCO2Y=growth respiration unlimited,limited by O2 and
C        unlimited by N,P
C
      RCO2R=RCO2RM*WFR(N,L,NZ,NY,NX)
      RCO2XM=RCO2RM-RMNCR
      RCO2X=RCO2R-RMNCR
      RCO2YM=AMAX1(0.0,RCO2XM)
      RCO2Y=AMAX1(0.0,RCO2X)
C
C     PRIMARY ROOT GROWTH RESPIRATION MAY BE LIMITED BY
C     NON-STRUCTURAL N,P AVAILABLE FOR GROWTH
C
C     FRTN=fraction of secondary root sink strength in axis
C     ZPOOLR,PPOOLR=non-structural N,P mass in root
C     CNRTS,CPRTS=N,P root growth yield
C     FNP=growth respiration limited by non-structural N,P
C     RCO2GM,RCO2G=growth respiration limited by N,P unltd,ltd by O2
C
      DMRTR=DMRTD*FRTN
      ZPOOLB=AMAX1(0.0,ZPOOLR(N,L,NZ,NY,NX))
      PPOOLB=AMAX1(0.0,PPOOLR(N,L,NZ,NY,NX))
      FNP=AMIN1(ZPOOLB*DMRTR/CNRTS(NZ,NY,NX)
     2,PPOOLB*DMRTR/CPRTS(NZ,NY,NX))
      IF(RCO2YM.GT.0.0)THEN
      RCO2GM=AMIN1(RCO2YM,FNP)
      ELSE
      RCO2GM=0.0
      ENDIF
      IF(RCO2Y.GT.0.0)THEN
      RCO2G=AMIN1(RCO2Y,FNP*WFR(N,L,NZ,NY,NX))
      ELSE
      RCO2G=0.0
      ENDIF
C
C     TOTAL NON-STRUCTURAL C,N,P USED IN PRIMARY ROOT GROWTH
C     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
C     ENTERED IN 'READQ'
C
C     CGRORM,CGROR=total non-structural C used in growth and
C        respiration unlimited,limited by O2
C     RCO2GM,RCO2G=growth respiration limited by N,P unlimited,limited
C        by O2
C     DMRTD=root C respiration vs nonstructural C consumption
C     GRTWGM,GRTWTG=root C growth unlimited,limited by O2
C     DMRT=root growth yield
C     ZADD1M,ZADD1,PADD1=nonstructural N,P unlimited,limited
C        by O2 used in growth
C     CNRDM,CNRDA=respiration for N assimilation unlimited,limited
C        by O2
C
      CGRORM=RCO2GM/DMRTD
      CGROR=RCO2G/DMRTD
      GRTWGM=CGRORM*DMRT(NZ,NY,NX)
      GRTWTG=CGROR*DMRT(NZ,NY,NX)
      ZADD1M=AMAX1(0.0,GRTWGM*CNRTW)
      ZADD1=AMAX1(0.0,AMIN1(FRTN*ZPOOLR(N,L,NZ,NY,NX),GRTWTG*CNRTW))
      PADD1=AMAX1(0.0,AMIN1(FRTN*PPOOLR(N,L,NZ,NY,NX),GRTWTG*CPRTW))
      CNRDM=AMAX1(0.0,1.70*ZADD1M)
      CNRDA=AMAX1(0.0,1.70*ZADD1)
C
C     PRIMARY ROOT GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
C     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION, ALSO
C     PRIMARY ROOT C LOSS FROM REMOBILIZATION AND CONSEQUENT
C     LITTERFALL
C
C     IDAY(1,=emergence date
C     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     RCCZR,RCCYR=min,max fractions for root C recycling
C     RCCXR,RCCQR=max fractions for root N,P recycling
C
      IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0
     2.AND.CCPOLR(N,L,NZ,NY,NX).GT.ZERO)THEN
      CCC=AMAX1(0.0,AMIN1(1.0
     1,CZPOLR(N,L,NZ,NY,NX)/(CZPOLR(N,L,NZ,NY,NX)
     2+CCPOLR(N,L,NZ,NY,NX)*CNKI)
     3,CPPOLR(N,L,NZ,NY,NX)/(CPPOLR(N,L,NZ,NY,NX)
     4+CCPOLR(N,L,NZ,NY,NX)*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX)
     2+CZPOLR(N,L,NZ,NY,NX)/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLR(N,L,NZ,NY,NX)/(CCPOLR(N,L,NZ,NY,NX)
     2+CPPOLR(N,L,NZ,NY,NX)/CPKI)))
      ELSE
      CCC=0.0
      CNC=0.0
      CPC=0.0
      ENDIF
      RCCC=RCCZR(IBTYP(NZ,NY,NX))+CCC*RCCYR(IBTYP(NZ,NY,NX))
      RCCN=CNC*RCCXR(IBTYP(NZ,NY,NX))
      RCCP=CPC*RCCQR(IBTYP(NZ,NY,NX))
C
C     RECOVERY OF REMOBILIZABLE N,P DURING PRIMARY ROOT REMOBILIZATION
C     DEPENDS ON ROOT NON-STRUCTURAL C:N:P
C
C     RCO2XM,RCO2X=diff between C respiration unlimited,limited
C        by O2 and maintenance respiration
C     SNCRM,SNCR=excess maintenance respiration unlimited,limited
C        by O2
C     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
C     WFR=constraint by O2 consumption on all root processes
C     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     FSNC1=fraction of primary root C to be remobilized
C
      IF(-RCO2XM.GT.0.0)THEN
      IF(-RCO2XM.LT.RTWT1(N,NR,NZ,NY,NX)*RCCC)THEN
      SNCRM=-RCO2XM
      ELSE
      SNCRM=AMAX1(0.0,RTWT1(N,NR,NZ,NY,NX)*RCCC)
      ENDIF
      ELSE
      SNCRM=0.0
      ENDIF
      IF(-RCO2X.GT.0.0)THEN
      IF(-RCO2X.LT.RTWT1(N,NR,NZ,NY,NX)*RCCC)THEN
      SNCR=-RCO2X
      ELSE
      SNCR=AMAX1(0.0,RTWT1(N,NR,NZ,NY,NX)*RCCC)
     2*WFR(N,L,NZ,NY,NX)
      ENDIF
      ELSE
      SNCR=0.0
      ENDIF
      IF(SNCR.GT.0.0.AND.RTWT1(N,NR,NZ,NY,NX)
     2.GT.ZEROP(NZ,NY,NX))THEN
      RCCR=RCCC*RTWT1(N,NR,NZ,NY,NX)
      RCZR=RTWT1N(N,NR,NZ,NY,NX)*(RCCN+(1.0-RCCN)*RCCC)
      RCPR=RTWT1P(N,NR,NZ,NY,NX)*(RCCP+(1.0-RCCP)*RCCC)
      IF(RCCR.GT.ZEROP(NZ,NY,NX))THEN
      FSNC1=AMAX1(0.0,AMIN1(1.0,SNCR/RCCR))
      ELSE
      FSNC1=1.0
      ENDIF
      ELSE
      RCCR=0.0
      RCZR=0.0
      RCPR=0.0
      FSNC1=0.0
      ENDIF
C
C     PRIMARY ROOT LITTERFALL CAUSED BY REMOBILIZATION
C
C     CSNC,ZSNC,PSNC=literfall C,N,P
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     FSNC1=fraction of primary root C to be remobilized
C     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
C     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
C     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,
C        1=non-woody
C
      DO 6355 M=1,4
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*FSNC1*RTWT1(N,NR,NZ,NY,NX)*FWODR(0,NZ)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*FSNC1*RTWT1N(N,NR,NZ,NY,NX)*FWODRN(0,NZ)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*FSNC1*RTWT1P(N,NR,NZ,NY,NX)*FWODRP(0,NZ)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX)
     2*FSNC1*(RTWT1(N,NR,NZ,NY,NX)-RCCR)*FWODR(1,NZ)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX)
     2*FSNC1*(RTWT1N(N,NR,NZ,NY,NX)-RCZR)*FWODRN(1,NZ)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX)
     2*FSNC1*(RTWT1P(N,NR,NZ,NY,NX)-RCPR)*FWODRP(1,NZ)
6355  CONTINUE
C
C     CONSUMPTION OF NON-STRUCTURAL C,N,P BY PRIMARY ROOTS
C
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     RMNCR=root maintenance respiration
C     RCO2R=respiration from non-structural C limited by O2
C     CGROR=total non-structural C used in growth and respiration
C        limited by O2
C     CNRDA=respiration for N assimilation unlimited,limited by O2
C     SNCR=excess maintenance respiration limited by O2
C     FSNC1=fraction of primary root C to be remobilized
C     RCCR,RCZR,RCPR=remobilization of C,N,P from senescing root
C     ZADD1,PADD1=nonstructural N,P ltd by O2 used in growth
C
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-AMIN1(RMNCR,RCO2R)
     2-CGROR-CNRDA-SNCR+FSNC1*RCCR*FWODR(1,NZ)
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)
     2-ZADD1+FSNC1*RCZR*FWODRN(1,NZ)
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)
     2-PADD1+FSNC1*RCPR
C
C     TOTAL PRIMARY ROOT RESPIRATION
C
C     RCO2TM,RCO2T=total C respiration unlimited,limited
C        by O2
C     RMNCR=root maintenance respiration
C     RCO2RM,RCO2R=respiration from non-structural C unlimited,limited
C        by O2
C     RCO2GM,RCO2G=growth respiration limited by N,P unlimited,limited
C        by O2
C     SNCRM,SNCR=excess maintenance respiration unlimited,limited
C        by O2
C     CNRDM,CNRDA=respiration for N assimilation unlimited,limited
C        by O2
C     RCO2A=total root respiration
C     RCO2M,RCO2N=RCO2A unlimited by O2,nonstructural C
C
      RCO2TM=AMIN1(RMNCR,RCO2RM)+RCO2GM+SNCRM+CNRDM
      RCO2T=AMIN1(RMNCR,RCO2R)+RCO2G+SNCR+CNRDA
C
C     ALLOCATE PRIMARY ROOT TOTAL RESPIRATION TO ALL SOIL LAYERS
C     THROUGH WHICH PRIMARY ROOTS GROW
C
C     RTDP1=primary root depth from soil surface
C     CDPTHZ=depth from soil surface to layer bottom
C     RTLG1=primary root length
C     SDPTH=seeding depth
C     FRCO2=fraction of primary root respiration attributed to layer
C     RCO2A=total root respiration
C     RCO2M,RCO2N=RCO2A unlimited by O2,nonstructural C
C     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
C
      IF(RTDP1(N,NR,NZ,NY,NX).GT.CDPTHZ(NG(NZ,NY,NX),NY,NX))THEN
      TFRCO2=0.0
      DO 5100 LL=NG(NZ,NY,NX),NINR(NR,NZ,NY,NX)
      IF(LL.LT.NINR(NR,NZ,NY,NX))THEN
      FRCO2=AMIN1(1.0,RTLG1(N,LL,NR,NZ,NY,NX)
     2/(RTDP1(N,NR,NZ,NY,NX)-SDPTH(NZ,NY,NX)))
      ELSE
      FRCO2=1.0-TFRCO2
      ENDIF
      TFRCO2=TFRCO2+FRCO2
      RCO2M(N,LL,NZ,NY,NX)=RCO2M(N,LL,NZ,NY,NX)+RCO2TM*FRCO2
      RCO2N(N,LL,NZ,NY,NX)=RCO2N(N,LL,NZ,NY,NX)+RCO2T*FRCO2
      RCO2A(N,LL,NZ,NY,NX)=RCO2A(N,LL,NZ,NY,NX)-RCO2T*FRCO2
C     IF(NZ.EQ.2)THEN
C     WRITE(*,9877)'RCO2A',I,J,NZ,NR,L,LL,N,NG(NZ,NY,NX)
C    2,NINR(NR,NZ,NY,NX),RCO2T,FRCO2,TFRCO2,RCO2A(N,LL,NZ,NY,NX)
C    3,RTLG1(N,LL,NR,NZ,NY,NX),RTDP1(N,NR,NZ,NY,NX)
C    4,SDPTH(NZ,NY,NX)
C     ENDIF
5100  CONTINUE
      ELSE
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)+RCO2TM
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)+RCO2T
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)-RCO2T
      ENDIF
C
C     ALLOCATE ANY NEGATIVE PRIMARY ROOT C,N,P GROWTH TO SECONDARY
C     ROOTS ON THE SAME AXIS IN THE SAME LAYER UNTIL SECONDARY ROOTS
C     HAVE DISAPPEARED
C
C     GRTWTG=primary root C growth limited by O2
C     GRTWTL,GRTWTN,GRTWTP=net primary root C,N,P growth
C     FSNC1=fraction of primary root C to be remobilized
C     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
C     ZADD1,PADD1=nonstructural N,P limited by O2 used in growth
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
C     RTLG2=secondary root length
C
      GRTWTL=GRTWTG-FSNC1*RTWT1(N,NR,NZ,NY,NX)
      GRTWTN=ZADD1-FSNC1*RTWT1N(N,NR,NZ,NY,NX)
      GRTWTP=PADD1-FSNC1*RTWT1P(N,NR,NZ,NY,NX)
      IF(GRTWTL.LT.0.0)THEN
      LX=MAX(1,L-1)
      DO 5105 LL=L,LX,-1
      GRTWTM=GRTWTL
      IF(GRTWTL.LT.0.0)THEN
      IF(GRTWTL.GT.-WTRT2(N,LL,NR,NZ,NY,NX))THEN
      RTLG2(N,LL,NR,NZ,NY,NX)=RTLG2(N,LL,NR,NZ,NY,NX)+GRTWTL
     2*RTLG2(N,LL,NR,NZ,NY,NX)/WTRT2(N,LL,NR,NZ,NY,NX)
      WTRT2(N,LL,NR,NZ,NY,NX)=WTRT2(N,LL,NR,NZ,NY,NX)+GRTWTL
      GRTWTL=0.0
      ELSE
      GRTWTL=GRTWTL+WTRT2(N,LL,NR,NZ,NY,NX)
      RTLG2(N,LL,NR,NZ,NY,NX)=0.0
      WTRT2(N,LL,NR,NZ,NY,NX)=0.0
      ENDIF
      ENDIF
      IF(GRTWTN.LT.0.0)THEN
      IF(GRTWTN.GT.-WTRT2N(N,LL,NR,NZ,NY,NX))THEN
      WTRT2N(N,LL,NR,NZ,NY,NX)=WTRT2N(N,LL,NR,NZ,NY,NX)+GRTWTN
      GRTWTN=0.0
      ELSE
      GRTWTN=GRTWTN+WTRT2N(N,LL,NR,NZ,NY,NX)
      WTRT2N(N,LL,NR,NZ,NY,NX)=0.0
      ENDIF
      ENDIF
      IF(GRTWTP.LT.0.0)THEN
      IF(GRTWTP.GT.-WTRT2P(N,LL,NR,NZ,NY,NX))THEN
      WTRT2P(N,LL,NR,NZ,NY,NX)=WTRT2P(N,LL,NR,NZ,NY,NX)+GRTWTP
      GRTWTP=0.0
      ELSE
      GRTWTP=GRTWTP+WTRT2P(N,LL,NR,NZ,NY,NX)
      WTRT2P(N,LL,NR,NZ,NY,NX)=0.0
      ENDIF
      ENDIF
C     WRITE(*,9876)'WTRT2',I,J,NZ,NR,LL,N
C    2,GRTWTL,GRTWTM,GRTWTG,FSNC1,SNCR,RCCR,RTWT1(N,NR,NZ,NY,NX)
C    3,WTRT2(1,LL,NR,NZ,NY,NX),WTRTL(1,LL,NZ,NY,NX)
C    3,WTRT2(2,LL,NR,NZ,NY,NX),WTRTL(2,LL,NZ,NY,NX)
C    4,RTLG2(1,LL,NR,NZ,NY,NX),RTLG1(1,LL,NR,NZ,NY,NX)
C    4,RTLG2(2,LL,NR,NZ,NY,NX),RTLG1(2,LL,NR,NZ,NY,NX)
C
C     CONCURRENT LOSS OF MYCORRHIZAE AND NODULES WITH LOSS
C     OF SECONDARY ROOTS
C
C     GRTWTM=negative primary root C growth
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass
C     FSNCM,FSNCP=fraction of mycorrhizal structural,nonstructural C
C        to be remobilized
C     WTRTL=active root C mass
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTRT2,WTRT2N,WTRT2P=mycorrhizal C,N,P mass
C     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,
C        1=non-woody
C     RTLG2=mycorrhizal length
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in mycorrhizae
C
      IF(GRTWTM.LT.0.0)THEN
      IF(WTRT2(1,LL,NR,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FSNCM=AMIN1(1.0,ABS(GRTWTM)/WTRT2(1,LL,NR,NZ,NY,NX))
      ELSE
      FSNCM=1.0
      ENDIF
      IF(WTRTL(1,LL,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FSNCP=AMIN1(1.0,ABS(GRTWTM)/WTRTL(1,LL,NZ,NY,NX))
      ELSE
      FSNCP=1.0
      ENDIF
      DO 6450 M=1,4
      CSNC(M,0,LL,NZ,NY,NX)=CSNC(M,0,LL,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*FSNCM*AMAX1(0.0,WTRT2(2,LL,NR,NZ,NY,NX))*FWODR(0,NZ)
      ZSNC(M,0,LL,NZ,NY,NX)=ZSNC(M,0,LL,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*FSNCM*AMAX1(0.0,WTRT2N(2,LL,NR,NZ,NY,NX))*FWODRN(0,NZ)
      PSNC(M,0,LL,NZ,NY,NX)=PSNC(M,0,LL,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*FSNCM*AMAX1(0.0,WTRT2P(2,LL,NR,NZ,NY,NX))*FWODRP(0,NZ)
      CSNC(M,1,LL,NZ,NY,NX)=CSNC(M,1,LL,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX)
     2*FSNCM*AMAX1(0.0,WTRT2(2,LL,NR,NZ,NY,NX))*FWODR(1,NZ)
      ZSNC(M,1,LL,NZ,NY,NX)=ZSNC(M,1,LL,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX)
     2*FSNCM*AMAX1(0.0,WTRT2N(2,LL,NR,NZ,NY,NX))*FWODRN(1,NZ)
      PSNC(M,1,LL,NZ,NY,NX)=PSNC(M,1,LL,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX)
     2*FSNCM*AMAX1(0.0,WTRT2P(2,LL,NR,NZ,NY,NX))*FWODRP(1,NZ)
      CSNC(M,1,LL,NZ,NY,NX)=CSNC(M,1,LL,NZ,NY,NX)+CFOPC(0,M,NZ,NY,NX)
     2*FSNCP*AMAX1(0.0,CPOOLR(2,LL,NZ,NY,NX))
      ZSNC(M,1,LL,NZ,NY,NX)=ZSNC(M,1,LL,NZ,NY,NX)+CFOPN(0,M,NZ,NY,NX)
     2*FSNCP*AMAX1(0.0,ZPOOLR(2,LL,NZ,NY,NX))
      PSNC(M,1,LL,NZ,NY,NX)=PSNC(M,1,LL,NZ,NY,NX)+CFOPP(0,M,NZ,NY,NX)
     2*FSNCP*AMAX1(0.0,PPOOLR(2,LL,NZ,NY,NX))
6450  CONTINUE
      RTLG2(2,LL,NR,NZ,NY,NX)=AMAX1(0.0,RTLG2(2,LL,NR,NZ,NY,NX))
     2*(1.0-FSNCM)
      WTRT2(2,LL,NR,NZ,NY,NX)=AMAX1(0.0,WTRT2(2,LL,NR,NZ,NY,NX))
     2*(1.0-FSNCM)
      WTRT2N(2,LL,NR,NZ,NY,NX)=AMAX1(0.0,WTRT2N(2,LL,NR,NZ,NY,NX))
     2*(1.0-FSNCM)
      WTRT2P(2,LL,NR,NZ,NY,NX)=AMAX1(0.0,WTRT2P(2,LL,NR,NZ,NY,NX))
     2*(1.0-FSNCM)
      CPOOLR(2,LL,NZ,NY,NX)=AMAX1(0.0,CPOOLR(2,LL,NZ,NY,NX))
     2*(1.0-FSNCP)
      ZPOOLR(2,LL,NZ,NY,NX)=AMAX1(0.0,ZPOOLR(2,LL,NZ,NY,NX))
     2*(1.0-FSNCP)
      PPOOLR(2,LL,NZ,NY,NX)=AMAX1(0.0,PPOOLR(2,LL,NZ,NY,NX))
     2*(1.0-FSNCP)
      ENDIF
C     IF(LL.EQ.12)THEN
C     WRITE(*,2126)'WITHMY',I,J,NZ,NR,LL,NN,CPOOLR(2,LL,NZ,NY,NX)
C    2,FSNCP
2126  FORMAT(A8,6I4,12E12.4)
C     ENDIF
5105  CONTINUE
      ENDIF
C
C     PRIMARY ROOT EXTENSION FROM ROOT GROWTH AND ROOT TURGOR
C
C     GRTLGL=primary root length extension
C     GRTWTG=primary root C growth limited by O2
C     RTLG1X=specific primary root length from �startq.f�
C     PP=PFT population
C     WFNRT=water function for root extension
C     FWOOD=C,N,P woody fraction in root:0=woody,1=non-woody
C     GRTWTL,GRTWTN,GRTWTP=net primary root C,N,P growth
C     RTDP1=primary root depth from soil surface
C     SDPTH=seeding depth
C     FSNC1=fraction of primary root C to be remobilized
C     RTLG1=primary root length
C     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
C     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
C     DLYR=soil layer thickness
C
      IF(GRTWTL.LT.0.0.AND.RTWT1(N,NR,NZ,NY,NX)
     2.GT.ZEROP(NZ,NY,NX))THEN
      GRTLGL=GRTWTG*RTLG1X(N,NZ,NY,NX)/PP(NZ,NY,NX)*WFNRT
     2+GRTWTL*(RTDP1(N,NR,NZ,NY,NX)-SDPTH(NZ,NY,NX))
     3/RTWT1(N,NR,NZ,NY,NX)
      ELSE
      GRTLGL=GRTWTG*RTLG1X(N,NZ,NY,NX)/PP(NZ,NY,NX)*WFNRT
      ENDIF
      IF(L.LT.NJ(NY,NX))THEN
      GRTLGL=AMIN1(DLYR(3,L1,NY,NX),GRTLGL)
      ENDIF
C
C     ALLOCATE PRIMARY ROOT GROWTH TO CURRENT
C     AND NEXT SOIL LAYER WHEN PRIMARY ROOTS EXTEND ACROSS LOWER
C     BOUNDARY OF CURRENT LAYER
C
C     GRTLGL=primary root length extension
C     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
C     CDPTHZ=depth from soil surface to layer bottom
C     RTDP1,RTDP1X=primary root depth from soil surface
C
      IF(GRTLGL.GT.ZEROP(NZ,NY,NX).AND.L.LT.NJ(NY,NX))THEN
      FGROL=AMAX1(0.0,AMIN1(1.0,(CDPTHZ(L,NY,NX)
     2-RTDP1(N,NR,NZ,NY,NX))/GRTLGL))
      IF(FGROL.LT.1.0)FGROL=0.0
      FGROZ=AMAX1(0.0,1.0-FGROL)
      ELSE
      FGROL=1.0
      FGROZ=0.0
      ENDIF
C
C     UPDATE STATE VARIABLES FOR PRIMARY ROOT LENGTH, GROWTH
C     AND AXIS NUMBER
C
C     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
C     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
C     GRTLGL=primary root length extension
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     FGROL,FGROZ=fraction of GRTLGL in current,next lower soil layer
C     WSRTL=total root protein C mass
C     CNWS,CPWS=protein:N,protein:P ratios from �startq.f�
C     RTLG1=primary root length
C
      RTWT1(N,NR,NZ,NY,NX)=RTWT1(N,NR,NZ,NY,NX)+GRTWTL
      RTWT1N(N,NR,NZ,NY,NX)=RTWT1N(N,NR,NZ,NY,NX)+GRTWTN
      RTWT1P(N,NR,NZ,NY,NX)=RTWT1P(N,NR,NZ,NY,NX)+GRTWTP
      RTDP1(N,NR,NZ,NY,NX)=RTDP1(N,NR,NZ,NY,NX)+GRTLGL
      WTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)+GRTWTL*FGROL
      WTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)+GRTWTN*FGROL
      WTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)+GRTWTP*FGROL
      WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX)
     2+AMIN1(CNWS(NZ,NY,NX)*WTRT1N(N,L,NR,NZ,NY,NX)
     2,CPWS(NZ,NY,NX)*WTRT1P(N,L,NR,NZ,NY,NX))
      RTLG1(N,L,NR,NZ,NY,NX)=RTLG1(N,L,NR,NZ,NY,NX)+GRTLGL*FGROL
C
C     TRANSFER STRUCTURAL, NONSTRUCTURAL C,N,P INTO NEXT SOIL LAYER
C     WHEN PRIMARY ROOT EXTENDS ACROSS LOWER BOUNDARY
C     OF CURRENT SOIL LAYER
C
C     FGROZ=fraction of GRTLGL in next lower soil layer
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     GRTWTL,GRTWTN,GRTWTP=net root C,N,P growth
C     WSRTL=total root protein C mass
C     CNWS,CPWS=protein:N,protein:P ratios from �startq.f�
C     WTRTD=root C mass
C     RTLG1=primary root length
C     GRTLGL=primary root length extension
C     FRTN=fraction of primary root sink strength in axis
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     PSIRT,PSIRG,PSIRO=root total,turgor,osmotic water potential
C     NINR=deepest root layer
C
      IF(FGROZ.GT.0.0)THEN
      WTRT1(N,L1,NR,NZ,NY,NX)=WTRT1(N,L1,NR,NZ,NY,NX)
     2+GRTWTL*FGROZ
      WTRT1N(N,L1,NR,NZ,NY,NX)=WTRT1N(N,L1,NR,NZ,NY,NX)
     2+GRTWTN*FGROZ
      WTRT1P(N,L1,NR,NZ,NY,NX)=WTRT1P(N,L1,NR,NZ,NY,NX)
     2+GRTWTP*FGROZ
      WSRTL(N,L1,NZ,NY,NX)=WSRTL(N,L1,NZ,NY,NX)
     2+AMIN1(CNWS(NZ,NY,NX)*WTRT1N(N,L1,NR,NZ,NY,NX)
     2,CPWS(NZ,NY,NX)*WTRT1P(N,L1,NR,NZ,NY,NX))
      WTRTD(N,L1,NZ,NY,NX)=WTRTD(N,L1,NZ,NY,NX)
     2+WTRT1(N,L1,NR,NZ,NY,NX)
      RTLG1(N,L1,NR,NZ,NY,NX)=RTLG1(N,L1,NR,NZ,NY,NX)+GRTLGL*FGROZ
      RRAD1(N,L1,NZ,NY,NX)=RRAD1(N,L,NZ,NY,NX)
      RTLGZ=RTLGZ+RTLG1(N,L1,NR,NZ,NY,NX)
      WTRTZ=WTRTZ+WTRT1(N,L1,NR,NZ,NY,NX)
      XFRC=FRTN*CPOOLR(N,L,NZ,NY,NX)
      XFRN=FRTN*ZPOOLR(N,L,NZ,NY,NX)
      XFRP=FRTN*PPOOLR(N,L,NZ,NY,NX)
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-XFRC
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)-XFRN
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)-XFRP
      CPOOLR(N,L1,NZ,NY,NX)=CPOOLR(N,L1,NZ,NY,NX)+XFRC
      ZPOOLR(N,L1,NZ,NY,NX)=ZPOOLR(N,L1,NZ,NY,NX)+XFRN
      PPOOLR(N,L1,NZ,NY,NX)=PPOOLR(N,L1,NZ,NY,NX)+XFRP
      PSIRT(N,L1,NZ,NY,NX)=PSIRT(N,L,NZ,NY,NX)
      PSIRO(N,L1,NZ,NY,NX)=PSIRO(N,L,NZ,NY,NX)
      PSIRG(N,L1,NZ,NY,NX)=PSIRG(N,L,NZ,NY,NX)
      NINR(NR,NZ,NY,NX)=MAX(NG(NZ,NY,NX),L+1)
C     IF(NZ.EQ.2)THEN
C     WRITE(*,9877)'INFIL',I,J,NZ,NR,L,N,NINR(NR,NZ,NY,NX)
C    2,FRTN,WTRTD(N,L1,NZ,NY,NX),CPOOLR(N,L1,NZ,NY,NX)
C    2,FGROZ,RTDP1(N,NR,NZ,NY,NX),GRTLGL,CDPTHZ(L,NY,NX)
C     ENDIF
      ENDIF
C     IF(L.EQ.12)THEN
C     WRITE(*,9877)'RCO21',I,J,NFZ,NX,NY,NZ,NR,L,L-1,L1,N
C    2,NINR(NR,NZ,NY,NX)
C    2,RCO2TM,RCO2T,RMNCR,RCO2RM,RCO2R,RCO2GM,RCO2G
C    3,RCO2XM,RCO2X,CGROR,SNCRM,SNCR,CNRDA,CPOOLR(N,L,NZ,NY,NX),FRTN
C    4,TFN4(L,NZ,NY,NX),CNPG,FDBKX(NB1(NZ,NY,NX),NZ,NY,NX),WFNRG(N,L)
C    5,TFN6(L),GRTWTG,GRTWTL,GRTLGL,RTLGL,FGROL,RTLG1(N,L,NR,NZ,NY,NX)
C    6,WTRT1(N,L,NR,NZ,NY,NX),RTDP1(N,NR,NZ,NY,NX),RTDP1X
C    3,RCO2M(N,L,NZ,NY,NX),RCO2A(N,L,NZ,NY,NX),WFR(N,L,NZ,NY,NX)
C    4,RTSK1(N,L,NR),RRAD1(N,L,NZ,NY,NX),RTDPP
C    5,PSIRG(N,L,NZ,NY,NX),WFNRT,WFNRR(N,L),FWOOD(1,NZ)
C    6,FGROZ,RTWT1(N,NR,NZ,NY,NX),FSNC1
C    9,ZADD1,PADD1,ZPOOLR(N,L,NZ,NY,NX),PPOOLR(N,L,NZ,NY,NX)
C    1,RUPNH4(N,L,NZ,NY,NX),RUPNO3(N,L,NZ,NY,NX)
C    2,CDPTHZ(L,NY,NX),CDPTHZ(L-1,NY,NX),CDPTHZ(L1,NY,NX)
9877  FORMAT(A8,12I4,100E12.4)
C     ENDIF
      ENDIF
C
C     TRANSFER PRIMARY ROOT C,N,P TO NEXT SOIL LAYER ABOVE THE
C     CURRENT SOIL LAYER WHEN NEGATIVE PRIMARY ROOT GROWTH FORCES
C     WITHDRAWAL FROM THE CURRENT SOIL LAYER AND ALL SECONDARY ROOTS
C     IN THE CURRENT SOIL LAYER HAVE BEEN LOST
C
C     NINR=deepest root layer
C     VOLX=soil layer volume excluding macropore, rocks
C     RTDP1X=primary root depth from soil surface
C     CDPTHZ=depth from soil surface to layer bottom
C     SDPTH=seeding depth
C     FRTN=fraction of primary root sink strength in axis
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C     RTLG1=primary root length
C     WSRTL=root protein C mass
C     WTRTD=root C mass
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C
      IF(L.EQ.NINR(NR,NZ,NY,NX))THEN
      DO 5115 LL=L,NG(NZ,NY,NX)+1,-1
      IF(VOLX(LL-1,NY,NX).GT.ZEROS2(NY,NX)
     2.AND.(RTDP1X.LT.CDPTHZ(LL-1,NY,NX)
     2.OR.RTDP1X.LT.SDPTH(NZ,NY,NX)))THEN
      IF(RLNT(N,LL).GT.ZEROP(NZ,NY,NX))THEN
      FRTN=(RTSK1(N,LL,NR)+RTSK2(N,LL,NR))/RLNT(N,LL)
      ELSE
      FRTN=1.0
      ENDIF
      DO 5110 NN=1,MY(NZ,NY,NX)
      WTRT1(NN,LL-1,NR,NZ,NY,NX)=WTRT1(NN,LL-1,NR,NZ,NY,NX)
     2+WTRT1(NN,LL,NR,NZ,NY,NX)
      WTRT1N(NN,LL-1,NR,NZ,NY,NX)=WTRT1N(NN,LL-1,NR,NZ,NY,NX)
     2+WTRT1N(NN,LL,NR,NZ,NY,NX)
      WTRT1P(NN,LL-1,NR,NZ,NY,NX)=WTRT1P(NN,LL-1,NR,NZ,NY,NX)
     2+WTRT1P(NN,LL,NR,NZ,NY,NX)
      WTRT2(NN,LL-1,NR,NZ,NY,NX)=WTRT2(NN,LL-1,NR,NZ,NY,NX)
     2+WTRT2(NN,LL,NR,NZ,NY,NX)
      WTRT2N(NN,LL-1,NR,NZ,NY,NX)=WTRT2N(NN,LL-1,NR,NZ,NY,NX)
     2+WTRT2N(NN,LL,NR,NZ,NY,NX)
      WTRT2P(NN,LL-1,NR,NZ,NY,NX)=WTRT2P(NN,LL-1,NR,NZ,NY,NX)
     2+WTRT2P(NN,LL,NR,NZ,NY,NX)
      RTLG1(NN,LL-1,NR,NZ,NY,NX)=RTLG1(NN,LL-1,NR,NZ,NY,NX)
     2+RTLG1(NN,LL,NR,NZ,NY,NX)
      WTRT1(NN,LL,NR,NZ,NY,NX)=0.0
      WTRT1N(NN,LL,NR,NZ,NY,NX)=0.0
      WTRT1P(NN,LL,NR,NZ,NY,NX)=0.0
      WTRT2(NN,LL,NR,NZ,NY,NX)=0.0
      WTRT2N(NN,LL,NR,NZ,NY,NX)=0.0
      WTRT2P(NN,LL,NR,NZ,NY,NX)=0.0
      RTLG1(NN,LL,NR,NZ,NY,NX)=0.0
      XFRC=FRTN*CPOOLR(NN,LL,NZ,NY,NX)
      XFRN=FRTN*ZPOOLR(NN,LL,NZ,NY,NX)
      XFRP=FRTN*PPOOLR(NN,LL,NZ,NY,NX)
      XFRW=FRTN*WSRTL(NN,L,NZ,NY,NX)
      XFRD=FRTN*WTRTD(NN,LL,NZ,NY,NX)
      CPOOLR(NN,LL,NZ,NY,NX)=CPOOLR(NN,LL,NZ,NY,NX)-XFRC
      ZPOOLR(NN,LL,NZ,NY,NX)=ZPOOLR(NN,LL,NZ,NY,NX)-XFRN
      PPOOLR(NN,LL,NZ,NY,NX)=PPOOLR(NN,LL,NZ,NY,NX)-XFRP
      WSRTL(NN,LL,NZ,NY,NX)=WSRTL(NN,LL,NZ,NY,NX)-XFRW
      WTRTD(NN,LL,NZ,NY,NX)=WTRTD(NN,LL,NZ,NY,NX)-XFRD
      CPOOLR(NN,LL-1,NZ,NY,NX)=CPOOLR(NN,LL-1,NZ,NY,NX)+XFRC
      ZPOOLR(NN,LL-1,NZ,NY,NX)=ZPOOLR(NN,LL-1,NZ,NY,NX)+XFRN
      PPOOLR(NN,LL-1,NZ,NY,NX)=PPOOLR(NN,LL-1,NZ,NY,NX)+XFRP
      WSRTL(NN,LL-1,NZ,NY,NX)=WSRTL(NN,LL-1,NZ,NY,NX)+XFRW
      WTRTD(NN,LL-1,NZ,NY,NX)=WTRTD(NN,LL-1,NZ,NY,NX)+XFRD
C     IF(LL.EQ.12)THEN
C     WRITE(*,2128)'WITHPR',I,J,NZ,NR,LL,NN,CPOOLR(N,LL,NZ,NY,NX)
C    2,FRTN,XFRC
2128  FORMAT(A8,6I4,12E12.4)
C     ENDIF
C
C     WITHDRAW GASES IN PRIMARY ROOTS
C
C     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2,O2,CH4,
C        N2O,NH3,H2
C     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,
C        N2O,NH3,H2
C     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,
C        N2O,NH3,H2
C     FRTN=fraction of primary root sink strength in axis
C
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-FRTN*(CO2A(NN,LL,NZ,NY,NX)
     2+CO2P(NN,LL,NZ,NY,NX))
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-FRTN*(OXYA(NN,LL,NZ,NY,NX)
     2+OXYP(NN,LL,NZ,NY,NX))
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-FRTN*(CH4A(NN,LL,NZ,NY,NX)
     2+CH4P(NN,LL,NZ,NY,NX))
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-FRTN*(Z2OA(NN,LL,NZ,NY,NX)
     2+Z2OP(NN,LL,NZ,NY,NX))
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-FRTN*(ZH3A(NN,LL,NZ,NY,NX)
     2+ZH3P(NN,LL,NZ,NY,NX))
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-FRTN*(H2GA(NN,LL,NZ,NY,NX)
     2+H2GP(NN,LL,NZ,NY,NX))
      CO2A(NN,LL,NZ,NY,NX)=(1.0-FRTN)*CO2A(NN,LL,NZ,NY,NX)
      OXYA(NN,LL,NZ,NY,NX)=(1.0-FRTN)*OXYA(NN,LL,NZ,NY,NX)
      CH4A(NN,LL,NZ,NY,NX)=(1.0-FRTN)*CH4A(NN,LL,NZ,NY,NX)
      Z2OA(NN,LL,NZ,NY,NX)=(1.0-FRTN)*Z2OA(NN,LL,NZ,NY,NX)
      ZH3A(NN,LL,NZ,NY,NX)=(1.0-FRTN)*ZH3A(NN,LL,NZ,NY,NX)
      H2GA(NN,LL,NZ,NY,NX)=(1.0-FRTN)*H2GA(NN,LL,NZ,NY,NX)
      CO2P(NN,LL,NZ,NY,NX)=(1.0-FRTN)*CO2P(NN,LL,NZ,NY,NX)
      OXYP(NN,LL,NZ,NY,NX)=(1.0-FRTN)*OXYP(NN,LL,NZ,NY,NX)
      CH4P(NN,LL,NZ,NY,NX)=(1.0-FRTN)*CH4P(NN,LL,NZ,NY,NX)
      Z2OP(NN,LL,NZ,NY,NX)=(1.0-FRTN)*Z2OP(NN,LL,NZ,NY,NX)
      ZH3P(NN,LL,NZ,NY,NX)=(1.0-FRTN)*ZH3P(NN,LL,NZ,NY,NX)
      H2GP(NN,LL,NZ,NY,NX)=(1.0-FRTN)*H2GP(NN,LL,NZ,NY,NX)
C     IF(L.EQ.12)THEN
C     WRITE(*,9868)'WITHDR',I,J,NZ,NR,LL,NN,NINR(NR,NZ,NY,NX)
C    2,FRTN,RTSK1(N,LL,NR),RTSK2(N,LL,NR),RLNT(N,LL)
C    2,WTRTD(NN,LL-1,NZ,NY,NX),WTRTD(NN,LL,NZ,NY,NX)
C    2,RTLG1(NN,LL-1,NR,NZ,NY,NX),RTLG1(NN,LL,NR,NZ,NY,NX)
C    2,RTLG2(NN,LL-1,NR,NZ,NY,NX),RTLG2(NN,LL,NR,NZ,NY,NX)
C    3,RTDP1(N,NR,NZ,NY,NX),RTDP1(NN,NR,NZ,NY,NX)
C    4,CPOOLR(NN,LL-1,NZ,NY,NX),CPOOLR(NN,LL,NZ,NY,NX)
C    4,WTRT1(NN,LL-1,NR,NZ,NY,NX),WTRT1(NN,LL,NR,NZ,NY,NX)
C    4,WTRT2(NN,LL-1,NR,NZ,NY,NX),WTRT2(NN,LL,NR,NZ,NY,NX)
9868  FORMAT(A8,7I4,100E24.16)
C      ENDIF
5110  CONTINUE
C
C     RESET ROOT NUMBER AND PRIMARY ROOT LENGTH
C
C     RTN2,RTNL=number of secondary root axes
C     RTN1=number of primary root axes
C     RTLG1=primary root length
C     CDPTHZ=depth from soil surface to layer bottom
C     RTDP1=primary root depth from soil surface
C     SDPTH=seeding depth
C
      RTNL(N,LL,NZ,NY,NX)=RTNL(N,LL,NZ,NY,NX)
     2-RTN2(N,LL,NR,NZ,NY,NX)
      RTNL(N,LL-1,NZ,NY,NX)=RTNL(N,LL-1,NZ,NY,NX)
     2+RTN2(N,LL,NR,NZ,NY,NX)
      RTN2(N,LL,NR,NZ,NY,NX)=0.0
      RTN1(N,LL,NZ,NY,NX)=RTN1(N,LL,NZ,NY,NX)-XRTN1
      IF(LL-1.GT.NG(NZ,NY,NX))THEN
      RTLG1(N,LL-1,NR,NZ,NY,NX)=DLYR(3,LL-1,NY,NX)
     2-(CDPTHZ(LL-1,NY,NX)-RTDP1(N,NR,NZ,NY,NX))
      ELSE
      RTLG1(N,LL-1,NR,NZ,NY,NX)=DLYR(3,LL-1,NY,NX)
     2-(CDPTHZ(LL-1,NY,NX)-RTDP1(N,NR,NZ,NY,NX))
     3-(SDPTH(NZ,NY,NX)-CDPTHZ(LL-2,NY,NX))
      ENDIF
C
C     WITHDRAW C,N,P FROM ROOT NODULES IN LEGUMES
C
C     INTYP=N2 fixation:1,2,3=rapid to slow root symbiosis
C     FRTN=fraction of primary root sink strength in axis
C     WTNDL,WTNDLN,WTNDLP=root bacterial C,N,P mass
C     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in root bacteria
C
      IF(INTYP(NZ,NY,NX).GE.1.AND.INTYP(NZ,NY,NX).LE.3)THEN
      XFRC=FRTN*WTNDL(LL,NZ,NY,NX)
      XFRN=FRTN*WTNDLN(LL,NZ,NY,NX)
      XFRP=FRTN*WTNDLP(LL,NZ,NY,NX)
      WTNDL(LL,NZ,NY,NX)=WTNDL(LL,NZ,NY,NX)-XFRC
      WTNDLN(LL,NZ,NY,NX)=WTNDLN(LL,NZ,NY,NX)-XFRN
      WTNDLP(LL,NZ,NY,NX)=WTNDLP(LL,NZ,NY,NX)-XFRP
      WTNDL(LL-1,NZ,NY,NX)=WTNDL(LL-1,NZ,NY,NX)+XFRC
      WTNDLN(LL-1,NZ,NY,NX)=WTNDLN(LL-1,NZ,NY,NX)+XFRN
      WTNDLP(LL-1,NZ,NY,NX)=WTNDLP(LL-1,NZ,NY,NX)+XFRP
      XFRC=FRTN*CPOOLN(LL,NZ,NY,NX)
      XFRN=FRTN*ZPOOLN(LL,NZ,NY,NX)
      XFRP=FRTN*PPOOLN(LL,NZ,NY,NX)
      CPOOLN(LL,NZ,NY,NX)=CPOOLN(LL,NZ,NY,NX)-XFRC
      ZPOOLN(LL,NZ,NY,NX)=ZPOOLN(LL,NZ,NY,NX)-XFRN
      PPOOLN(LL,NZ,NY,NX)=PPOOLN(LL,NZ,NY,NX)-XFRP
      CPOOLN(LL-1,NZ,NY,NX)=CPOOLN(LL-1,NZ,NY,NX)+XFRC
      ZPOOLN(LL-1,NZ,NY,NX)=ZPOOLN(LL-1,NZ,NY,NX)+XFRN
      PPOOLN(LL-1,NZ,NY,NX)=PPOOLN(LL-1,NZ,NY,NX)+XFRP
C     IF(NZ.EQ.2.AND.LL-1.EQ.6)THEN
C     WRITE(*,9869)'WITHDRN',I,J,NFZ,NZ,NR,LL,NN,NINR(NR,NZ,NY,NX)
C    2,XFRN,WTNDLN(LL,NZ,NY,NX),ZPOOLN(LL,NZ,NY,NX)
C    3,RTDP1(N,NR,NZ,NY,NX)
9869  FORMAT(A8,8I4,12E12.4)
C     ENDIF
      ENDIF
      NINR(NR,NZ,NY,NX)=MAX(NG(NZ,NY,NX),LL-1)
      ELSE
      GO TO 5120
      ENDIF
5115  CONTINUE
      ENDIF
5120  CONTINUE
C
C     REMOVE ANY NEGATIVE ROOT MASS FROM NONSTRUCTURAL C
C
C     WTRT1,WTRT2=primary,secondary root C mass in soil layer
C     CPOOLR=non-structural C mass in root
C
      IF(WTRT1(N,L,NR,NZ,NY,NX).LT.0.0)THEN
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)+WTRT1(N,L,NR,NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=0.0
      ENDIF
      IF(WTRT2(N,L,NR,NZ,NY,NX).LT.0.0)THEN
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX)
      WTRT2(N,L,NR,NZ,NY,NX)=0.0
      ENDIF
C     IF(L.EQ.12)THEN
C     WRITE(*,2129)'REMOVE',I,J,NZ,NR,L,N,CPOOLR(N,L,NZ,NY,NX)
2129  FORMAT(A8,6I4,12E12.4)
C     ENDIF
C
C     TOTAL PRIMARY ROOT LENGTH AND MASS
C
C     RTLGZ=total primary root length
C     WTRTZ=total primary root C mass
C     RTLG1=primary root length in soil layer
C     WTRT1=primary root C mass in soil layer
C     NINR=deepest root layer
C
      RTLGZ=RTLGZ+RTLG1(N,L,NR,NZ,NY,NX)
      WTRTZ=WTRTZ+WTRT1(N,L,NR,NZ,NY,NX)
      NINR(NR,NZ,NY,NX)=MIN(NINR(NR,NZ,NY,NX),NJ(NY,NX))
      IF(L.EQ.NINR(NR,NZ,NY,NX))NRX(N,NR)=1
      ENDIF
      ENDIF
      RTLGZ=RTLGZ+RTLG1(N,L,NR,NZ,NY,NX)
      WTRTZ=WTRTZ+WTRT1(N,L,NR,NZ,NY,NX)
C     ENDIF
      ENDIF
      NIX(NZ,NY,NX)=MAX(NIX(NZ,NY,NX),NINR(NR,NZ,NY,NX))
5050  CONTINUE
C
C     DRAW FROM ROOT NON-STRUCTURAL TO SEASONAL STORAGE POOL WHEN
C     SEASONAL STORAGE POOL IS DEPLETED
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     WTRTL,WTRT=total root C mass
C     WTRVC=storage C
C     XFRX=maximum storage C content for remobilization
C        from stalk,root reserves
C     CPOOLR=non-structural C mass in root
C     XFRC=root-storage nonstructural C transfer
C
      IF(ISTYP(NZ,NY,NX).NE.0
     2.AND.L.LE.NIX(NZ,NY,NX))THEN
      IF(WTRTL(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WTRT(NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WTRVC(NZ,NY,NX).LT.XFRX*WTRT(NZ,NY,NX))THEN
      FWTRT=WTRTL(N,L,NZ,NY,NX)/WTRT(NZ,NY,NX)
      WTRTLX=WTRTL(N,L,NZ,NY,NX)
      WTRTTX=WTRT(NZ,NY,NX)*FWTRT
      WTRTTT=WTRTLX+WTRTTX
      CPOOLX=AMAX1(0.0,CPOOLR(N,L,NZ,NY,NX))
      WTRVCX=AMAX1(0.0,WTRVC(NZ,NY,NX)*FWTRT)
      CPOOLD=(WTRVCX*WTRTLX-CPOOLX*WTRTTX)/WTRTTT
      XFRC=AMIN1(0.0,XFRY*CPOOLD)
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)+XFRC
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)-XFRC
C     WRITE(*,3471)'RVC',I,J,NX,NY,NZ,L,N
C    2,XFRC,CPOOLR(N,L,NZ,NY,NX),WTRTD(N,L,NZ,NY,NX)
C    3,WTRVC(NZ,NY,NX),WTRT(NZ,NY,NX),FWTRT,XFRY,CPOOLD
C    4,WTRVCX,WTRTLX,CPOOLX,WTRTTX,WTRTTT
3471  FORMAT(A8,7I4,30E12.4)
      ENDIF
      ENDIF
C
C     ROOT AND MYCORRHIZAL LENGTH, DENSITY, VOLUME, RADIUS, AREA
C     TO CALCULATE WATER AND NUTRIENT UPTAKE IN 'UPTAKE'
C
C     RTLGZ,RTLGX=total primary root length per plant,population
C     PP=PFT population
C     WTRTZ=total primary root C mass
C     RTLGL=total secondary root length
C     WTRTX=total secondary root C mass
C     RTLGT=total root length
C     WTRTT=total root C mass
C     FWOOD=C woody fraction in root:0=woody,1=non-woody
C     RTDNP,RTLGP=root length density,root length per plant
C     RTVL,RTVLW,RTVLP=root or myco total,aqueous,gaseous volume
C     RRAD1,RRAD2=primary,secondary root radius
C     RTARP=root surface area per plant
C     RTLGA=average secondary root length
C     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=loss of root CO2,O2,CH4,
C        N2O,NH3,H2
C     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
C     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
C
      RTLGX=RTLGZ*PP(NZ,NY,NX)
      RTLGT=RTLGL+RTLGX
      WTRTT=WTRTX+WTRTZ
      IF(RTLGT.GT.ZEROP(NZ,NY,NX).AND.WTRTT.GT.ZEROP(NZ,NY,NX)
     2.AND.PP(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RTLGP(N,L,NZ,NY,NX)=RTLGT/PP(NZ,NY,NX)
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
      RTDNP(N,L,NZ,NY,NX)=RTLGP(N,L,NZ,NY,NX)/DLYR(3,L,NY,NX)
      ELSE
      RTDNP(N,L,NZ,NY,NX)=0.0
      ENDIF
      RTVL=AMAX1(RTAR1X(N,NZ,NY,NX)*RTLGX+RTAR2X(N,NZ,NY,NX)*RTLGL
     2,WTRTT*DMVL(N,NZ,NY,NX)*PSIRG(N,L,NZ,NY,NX))
      RTVLP(N,L,NZ,NY,NX)=PORT(N,NZ,NY,NX)*RTVL
      RTVLW(N,L,NZ,NY,NX)=(1.0-PORT(N,NZ,NY,NX))*RTVL
      RRAD1(N,L,NZ,NY,NX)=AMAX1(RRAD1X(N,NZ,NY,NX)
     2,(1.0+PSIRT(N,L,NZ,NY,NX)/EMODR)*RRAD1M(N,NZ,NY,NX))
      RRAD2(N,L,NZ,NY,NX)=AMAX1(RRAD2X(N,NZ,NY,NX)
     2,(1.0+PSIRT(N,L,NZ,NY,NX)/EMODR)*RRAD2M(N,NZ,NY,NX))
      RTAR=6.283*RRAD1(N,L,NZ,NY,NX)*RTLGX
     2+6.283*RRAD2(N,L,NZ,NY,NX)*RTLGL
      IF(RTNL(N,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RTLGA(N,L,NZ,NY,NX)=AMAX1(RTLGAX,RTLGL/RTNL(N,L,NZ,NY,NX))
      ELSE
      RTLGA(N,L,NZ,NY,NX)=RTLGAX
      ENDIF
      RTARP(N,L,NZ,NY,NX)=RTAR/PP(NZ,NY,NX)
      ELSE
      RTLGP(N,L,NZ,NY,NX)=0.0
      RTDNP(N,L,NZ,NY,NX)=0.0
      RTVLP(N,L,NZ,NY,NX)=0.0
      RTVLW(N,L,NZ,NY,NX)=0.0
      RRAD1(N,L,NZ,NY,NX)=RRAD1M(N,NZ,NY,NX)
      RRAD2(N,L,NZ,NY,NX)=RRAD2M(N,NZ,NY,NX)
      RTARP(N,L,NZ,NY,NX)=0.0
      RTLGA(N,L,NZ,NY,NX)=RTLGAX
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-(CO2A(N,L,NZ,NY,NX)
     2+CO2P(N,L,NZ,NY,NX))
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-(OXYA(N,L,NZ,NY,NX)
     2+OXYP(N,L,NZ,NY,NX))
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-(CH4A(N,L,NZ,NY,NX)
     2+CH4P(N,L,NZ,NY,NX))
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-(Z2OA(N,L,NZ,NY,NX)
     2+Z2OP(N,L,NZ,NY,NX))
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-(ZH3A(N,L,NZ,NY,NX)
     2+ZH3P(N,L,NZ,NY,NX))
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-(H2GA(N,L,NZ,NY,NX)
     2+H2GP(N,L,NZ,NY,NX))
      CO2A(N,L,NZ,NY,NX)=0.0
      OXYA(N,L,NZ,NY,NX)=0.0
      CH4A(N,L,NZ,NY,NX)=0.0
      Z2OA(N,L,NZ,NY,NX)=0.0
      ZH3A(N,L,NZ,NY,NX)=0.0
      H2GA(N,L,NZ,NY,NX)=0.0
      CO2P(N,L,NZ,NY,NX)=0.0
      OXYP(N,L,NZ,NY,NX)=0.0
      CH4P(N,L,NZ,NY,NX)=0.0
      Z2OP(N,L,NZ,NY,NX)=0.0
      ZH3P(N,L,NZ,NY,NX)=0.0
      H2GP(N,L,NZ,NY,NX)=0.0
      ENDIF
      ENDIF
5000  CONTINUE
5010  CONTINUE
C
C     ADD SEED DIMENSIONS TO ROOT DIMENSIONS (ONLY IMPORTANT DURING
C     GERMINATION)
C
C     RTLGP=root length per plant
C     SDVL,SDLG,SDAR=seed volume,length,area from �startq.f�
C     DLYR=layer thickness
C     RTDNP=root length density per plant
C     RTVLW,RTVLP=root water,air volume
C     RTARP=root surface area per plant
C
      RTLGP(1,NG(NZ,NY,NX),NZ,NY,NX)=RTLGP(1,NG(NZ,NY,NX),NZ,NY,NX)
     2+SDLG(NZ,NY,NX)
      IF(DLYR(3,NG(NZ,NY,NX),NY,NX).GT.ZERO)THEN
      RTDNP(1,NG(NZ,NY,NX),NZ,NY,NX)=RTLGP(1,NG(NZ,NY,NX),NZ,NY,NX)
     2/DLYR(3,NG(NZ,NY,NX),NY,NX)
      ELSE
      RTDNP(1,NG(NZ,NY,NX),NZ,NY,NX)=0.0
      ENDIF
      RTVL=RTVLP(1,NG(NZ,NY,NX),NZ,NY,NX)
     2+RTVLW(1,NG(NZ,NY,NX),NZ,NY,NX)+SDVL(NZ,NY,NX)*PP(NZ,NY,NX)
      RTVLP(1,NG(NZ,NY,NX),NZ,NY,NX)=PORT(1,NZ,NY,NX)*RTVL
      RTVLW(1,NG(NZ,NY,NX),NZ,NY,NX)=(1.0-PORT(1,NZ,NY,NX))*RTVL
      RTARP(1,NG(NZ,NY,NX),NZ,NY,NX)=RTARP(1,NG(NZ,NY,NX),NZ,NY,NX)
     2+SDAR(NZ,NY,NX)
C
C     PLANT DEATH IF STORAGE EXHAUSTED
C
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IDTHB=branch living flag: 0=alive,1=dead
C     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
C     PPX,PP=PFT population per m2,grid cell
C     VOLWP,VOLWC,VOLWQ=water volume in canopy,on canopy,standing dead
C        surfaces
C
      IF((WTRVC(NZ,NY,NX).LE.ZEROP2(NZ,NY,NX)
     2.OR.WTRVN(NZ,NY,NX).LE.ZEROP2(NZ,NY,NX)
     2.OR.WTRVP(NZ,NY,NX).LE.ZEROP2(NZ,NY,NX))
     2.AND.ISTYP(NZ,NY,NX).NE.0)THEN
      DO 301 NB=1,NBR(NZ,NY,NX)
      IDTHB(NB,NZ,NY,NX)=1
301   CONTINUE
      IDTHP(NZ,NY,NX)=1
      IDTHR(NZ,NY,NX)=1
      PPX(NZ,NY,NX)=0.0
      PP(NZ,NY,NX)=0.0
      VOLWQ(NZ,NY,NX)=VOLWQ(NZ,NY,NX)+VOLWC(NZ,NY,NX)+VOLWP(NZ,NY,NX)
      VOLWC(NZ,NY,NX)=0.0
      VOLWP(NZ,NY,NX)=0.0
      WRITE(*,1115)'DEATHR',IYRC,I,J,NFZ,NX,NY,NZ,VOLWP(NZ,NY,NX)
     2,WTRVC(NZ,NY,NX),WTRVN(NZ,NY,NX),WTRVP(NZ,NY,NX)
     3,ZEROP2(NZ,NY,NX)
1115  FORMAT(A8,7I4,12E12.4)
      ENDIF
C
C     ROOT N2 FIXATION (RHIZOBIA)
C
C     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
C
      IF(INTYP(NZ,NY,NX).GE.1.AND.INTYP(NZ,NY,NX).LE.3)THEN
      DO 5400 L=NU(NY,NX),NIX(NZ,NY,NX)
      IF(WTRTD(1,L,NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
C
C     INITIAL INFECTION
C
C     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C     WTNDI=initial bacterial mass at infection
C     AREA=grid cell area
C     CNND,CPND=bacterial N:C,P:C ratio from PFT file
C
      IF(NFZ.EQ.1.AND.WTNDL(L,NZ,NY,NX).LE.0.0
     2.AND.ICHKF.EQ.0)THEN
      WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)
     2+WTNDI*AREA(3,NU(NY,NX),NY,NX)
      WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)
     2+WTNDI*AREA(3,NU(NY,NX),NY,NX)*CNND(NZ,NY,NX)
      WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)
     2+WTNDI*AREA(3,NU(NY,NX),NY,NX)*CPND(NZ,NY,NX)
      ENDIF
C
C     O2-UNCONSTRAINED RESPIRATION RATES BY HETEROTROPHIC AEROBES
C     IN NODULE FROM SPECIFIC OXIDATION RATE, ACTIVE BIOMASS,
C     NON-STRUCTURAL C CONCENTRATION, MICROBIAL C:N:P FACTOR,
C     AND TEMPERATURE
C
C     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concentration
C        in bacteria
C     ZCPOLI,ZPPOLI=nonstructural N:C,N:P
C     CNKI,CPKI=nonstructural N,P inhibition constant on growth,
C        nutrient recycling
C     FCNPF=N,P constraint to bacterial activity
C     CCNDLF=nodule:root mass ratio
C
      IF(WTNDL(L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      CCPOLN=AMAX1(0.0,CPOOLN(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX))
      CZPOLN=AMAX1(0.0,ZPOOLN(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX))
      CPPOLN=AMAX1(0.0,PPOOLN(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX))
      ELSE
      CCPOLN=1.0
      CZPOLN=1.0
      CPPOLN=1.0
      ENDIF
      IF(CCPOLN.GT.ZERO)THEN
      ZCPOLI=CZPOLN/CCPOLN
      ELSE
      ZCPOLI=0.0
      ENDIF
      IF(CPPOLN.GT.ZERO)THEN
      ZPPOLI=CZPOLN/CPPOLN
      ELSE
      ZPPOLI=0.0
      ENDIF
      IF(CCPOLN.GT.ZERO)THEN
      CCC=AMAX1(0.0,AMIN1(1.0
     1,CZPOLN/(CZPOLN+CCPOLN*CNKI)
     2,CPPOLN/(CPPOLN+CCPOLN*CPKI)))
      CNC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLN/(CCPOLN+CZPOLN/CNKI)))
      CPC=AMAX1(0.0,AMIN1(1.0
     1,CCPOLN/(CCPOLN+CPPOLN/CPKI)))
      ELSE
      CCC=0.0
      CNC=0.0
      CPC=0.0
      ENDIF
      IF(WTNDL(L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FCNPF=AMIN1(1.0
     2,SQRT(WTNDLN(L,NZ,NY,NX)/(WTNDL(L,NZ,NY,NX)*CNND(NZ,NY,NX)))
     3,SQRT(WTNDLP(L,NZ,NY,NX)/(WTNDL(L,NZ,NY,NX)*CPND(NZ,NY,NX))))
      ELSE
      FCNPF=1.0
      ENDIF
      IF(WTRTD(1,L,NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
      CNDLR=AMAX1(0.0,WTNDL(L,NZ,NY,NX)/WTRTD(1,L,NZ,NY,NX))
      ELSE
      CNDLR=0.0
      ENDIF
C
C     RESPIRATION FROM NON-STRUCTURAL C DETERMINED BY TEMPERATURE,
C     NON-STRUCTURAL C:N:P WITH EXCESS N FEEDBACK
C
C     RCNDLM=respiration from non-structural C unltd by O2
C     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C     VMXO=specific respiration rate by bacterial N2 fixers
C     WTNDL=bacterial C mass
C     TFN4=temperature function for root growth
C     FCNPF=N,P constraint to bacterial activity
C     WFNRG=growth function of root water potential
C     ZCPOLI,ZPPOLI=nonstructural N:C,N:P
C
      RCNDLM=AMAX1(0.0,AMIN1(CPOOLN(L,NZ,NY,NX)
     2,VMXO*WTNDL(L,NZ,NY,NX))*FCNPF*TFN4(L,NZ,NY,NX)
     2/(1.0+AMAX1(ZCPOLI/ZCKI,ZPPOLI/ZPKI))
     3*WFNRG(1,L)*XNFH)
C     CPOOLNX=CPOOLN(L,NZ,NY,NX)
C
C     O2-LIMITED NODULE RESPIRATION FROM 'WFR' IN 'UPTAKE'
C
C     RCNDL=respiration from non-structural C ltd by O2
C     WFR=constraint by O2 consumption on all root processes
C
      RCNDL=RCNDLM*WFR(1,L,NZ,NY,NX)
C
C     NODULE MAINTENANCE RESPIRATION FROM SOIL TEMPERATURE,
C     NODULE STRUCTURAL N
C
C     RMNDL=bacterial maintenance respiration
C     RMPLT=specific maintenance respiration rate (g C g-1 N h-1)
C     TFN6=temperature function for root maintenance respiration
C     WTNDLN=bacterial N mass
C     WFNRR=respiration function of root water potential
C
      RMNDL=AMAX1(0.0,RMPLT*TFN6(L)*WTNDLN(L,NZ,NY,NX)
     2*WFNRR(1,L)*XNFH)
C
C     NODULE GROWTH RESPIRATION FROM TOTAL - MAINTENANCE
C     IF > 0 DRIVES GROWTH, IF < 0 DRIVES REMOBILIZATION
C
C     RXNDLM,RXNDL=difference between non-structural C respiration and
C        maintenance respiration unlimited,limited by O2
C     RGNDLM,RGNDL=growth respiration unlimited by N,P and
C        unlimited,limited by O2
C     RSNDLM,RSNDL=excess maintenance respiration unlimited,limited
C        by O2
C
      RXNDLM=RCNDLM-RMNDL
      RXNDL=RCNDL-RMNDL
      RGNDLM=AMAX1(0.0,RXNDLM)
      RGNDL=AMAX1(0.0,RXNDL)
      RSNDLM=AMAX1(0.0,-RXNDLM)
      RSNDL=AMAX1(0.0,-RXNDL)
C
C     NODULE N2 FIXATION FROM GROWTH RESPIRATION, FIXATION ENERGY
C     REQUIREMENT AND NON-STRUCTURAL C:N:P PRODUCT INHIBITION,
C     CONSTRAINED BY MICROBIAL N REQUIREMENT
C
C     RGN2P=respiration requirement to maintain bacterial N:C ratio
C     WTNDL,WTNDLN=bacterial C,N mass
C     CNND=bacterial N:C ratio from PFT file
C     EN2F=N fixation yield from C oxidation (g N g-1 C)
C     RGNDL=growth respiration unlimited by N,P
C     RGN2F=respiration for N2 fixation
C     RUPNF,UPNF=layer,total root N2 fixation
C
      RGN2P=AMAX1(0.0,WTNDL(L,NZ,NY,NX)*CNND(NZ,NY,NX)
     2-WTNDLN(L,NZ,NY,NX))/EN2F
      IF(RGNDL.GT.ZEROP(NZ,NY,NX))THEN
      RGN2F=(RGNDL*RGN2P/(RGNDL+RGN2P))
      ELSE
      RGN2F=0.0
      ENDIF
      RUPNF(L,NZ,NY,NX)=RGN2F*EN2F
      UPNF(NZ,NY,NX)=UPNF(NZ,NY,NX)+RUPNF(L,NZ,NY,NX)
C
C     NODULE C,N,P REMOBILIZATION AND DECOMPOSITION
C
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     RCCZN,RCCYN=min,max fractions for bacteria C recycling
C     RCCXN,RCCQN=max fractions for bacteria N,P recycling
C     WTRTD=root C mass
C     SPNDX=specific bacterial decomposition rate at current CNDLR
C     SPNDL=specific decomposition rate by bacterial N2 fixers
C     CNDLR=bacteria:root ratio
C     CNDLI=bacteria:root ratio controlling bacterial decomposition
C     TFN4=temperature function for root growth
C     WFNRG=growth function of root water potential
C     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
C     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
C     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
C
      RCCC=RCCZN+CCC*RCCYN
      RCCN=CNC*RCCXN
      RCCP=CPC*RCCQN
      SPNDX=AMIN1(1.0,SPNDL*CNDLR/CNDLI(INTYP(NZ,NY,NX)))
     2*SQRT(TFN4(L,NZ,NY,NX)*WFNRG(1,L))*XNFH
      RXNDLC=SPNDX*WTNDL(L,NZ,NY,NX)
      RXNDLN=SPNDX*WTNDLN(L,NZ,NY,NX)
      RXNDLP=SPNDX*WTNDLP(L,NZ,NY,NX)
      RDNDLC=RXNDLC*(1.0-RCCC)
      RDNDLN=RXNDLN*(1.0-RCCC)*(1.0-RCCN)
      RDNDLP=RXNDLP*(1.0-RCCC)*(1.0-RCCP)
      RCNDLC=RXNDLC-RDNDLC
      RCNDLN=RXNDLN-RDNDLN
      RCNDLP=RXNDLP-RDNDLP
C
C     TOTAL NON-STRUCTURAL C,N,P USED IN NODULE GROWTH
C     AND GROWTH RESPIRATION DEPENDS ON GROWTH YIELD
C     ENTERED IN 'READQ'
C
C     CGNDL=total non-structural C used in bacterial growth
C        and growth respiration
C     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C     RMNDL=bacterial maintenance respiration
C     RCNDL=respiration from non-structural C
C     RCNDLC=bacterial C decomposition to recycling
C     RGNDL=growth respiration ltd by O2
C     RGN2F=respiration for N2 fixation
C     GRNDG=bacterial growth
C     DMND=bacterial growth yield
C     RGNDG=bacterial respiration for growth and N2 fixation
C     ZADDN,PADDN=nonstructural N,P used in growth
C     CNND,CPND=bacterial N:C,P:C ratio from PFT file
C     CCPOLN,CZPOLN,CPPOLN=nonstructural C,N,P concentration
C        in bacteria
C     CZKM,CPKM=Km for nonstructural N,P uptake by bacteria
C
      CGNDL=AMIN1(CPOOLN(L,NZ,NY,NX)-AMIN1(RMNDL,RCNDL)
     2-RGN2F+RCNDLC,(RGNDL-RGN2F)/(1.0-DMND(NZ,NY,NX)))
      GRNDG=CGNDL*DMND(NZ,NY,NX)
      RGNDG=RGN2F+CGNDL*(1.0-DMND(NZ,NY,NX))
      ZADDN=AMAX1(0.0,AMIN1(ZPOOLN(L,NZ,NY,NX)
     2,GRNDG*CNND(NZ,NY,NX)))*CZPOLN/(CZPOLN+CZKM)
      PADDN=AMAX1(0.0,AMIN1(PPOOLN(L,NZ,NY,NX)
     2,GRNDG*CPND(NZ,NY,NX)))*CPPOLN/(CPPOLN+CPKM)
C
C     NODULE SENESCENCE
C
C     RSNDL=excess maintenance respiration
C     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C     RCCC,RCCN,RCCP=remobilization coefficient for C,N,P
C     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
C     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
C     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
C
      IF(RSNDL.GT.0.0.AND.WTNDL(L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      RXNSNC=RSNDL
      RXNSNN=RXNSNC*WTNDLN(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX)
      RXNSNP=RXNSNC*WTNDLP(L,NZ,NY,NX)/WTNDL(L,NZ,NY,NX)
      RDNSNC=RXNSNC*(1.0-RCCC)
      RDNSNN=RXNSNN*(1.0-RCCC)*(1.0-RCCN)
      RDNSNP=RXNSNP*(1.0-RCCC)*(1.0-RCCP)
      RCNSNC=RXNSNC-RDNSNC
      RCNSNN=RXNSNN-RDNSNN
      RCNSNP=RXNSNP-RDNSNP
      ELSE
      RXNSNC=0.0
      RXNSNN=0.0
      RXNSNP=0.0
      RDNSNC=0.0
      RDNSNN=0.0
      RDNSNP=0.0
      RCNSNC=0.0
      RCNSNN=0.0
      RCNSNP=0.0
      ENDIF
C
C     TOTAL NODULE RESPIRATION
C
C     RCO2TM,RCO2T=total C respiration unlimited,limited by O2
C     TCO2T,TCO2A=total,above-ground PFT respiration
C     RMNDL=bacterial maintenance respiration
C     RCNDL=respiration from non-structural C
C     RGNDG=bacterial respiration for growth and N2 fixation
C     RCNSNC=bacterial C senescence to recycling
C     RCO2A=total root respiration
C     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
C
      RCO2TM=AMIN1(RMNDL,RCNDLM)+RGNDLM+RCNSNC
      RCO2T=AMIN1(RMNDL,RCNDL)+RGNDG+RCNSNC
      RCO2M(1,L,NZ,NY,NX)=RCO2M(1,L,NZ,NY,NX)+RCO2TM
      RCO2N(1,L,NZ,NY,NX)=RCO2N(1,L,NZ,NY,NX)+RCO2T
      RCO2A(1,L,NZ,NY,NX)=RCO2A(1,L,NZ,NY,NX)-RCO2T
C
C     NODULE LITTERFALL CAUSED BY REMOBILIZATION
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition and
C        senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     RDNDLC,RDNDLN,RDNDLP=bacterial C,N,P decomposition to litterfall
C     RDNSNC,RDNSNC,RDNSNP=bacterial C,N,P senescence to litterfall
C
      DO 6370 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX)
     2*(RDNDLC+RDNSNC)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX)
     2*(RDNDLN+RDNSNN)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX)
     2*(RDNDLP+RDNSNP)
6370  CONTINUE
C
C     CONSUMPTION OF NON-STRUCTURAL C,N,P BY NODULE
C
C     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C     RMNDL=bacterial maintenance respiration
C     RCNDL=respiration from non-structural C
C     RGN2F=respiration for N2 fixation
C     CGNDL=total non-structural C used in bacterial growth
C        and growth respiration
C     RCNDLC,RCNDLN,RCNDLP=bacterial C,N,P decomposition to recycling
C     RCNSNC,RCNSNC,RCNSNP=bacterial C,N,P senescence to recycling
C     ZADDN,PADDN=nonstructural N,P used in growth
C     RUPNF=root N2 fixation
C
      CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)-AMIN1(RMNDL,RCNDL)
     2-RGN2F-CGNDL+RCNDLC
      ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)-ZADDN+RCNDLN+RCNSNN
     2+RUPNF(L,NZ,NY,NX)
      PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)-PADDN+RCNDLP+RCNSNP
C
C     UPDATE STATE VARIABLES FOR NODULE C, N, P
C
C     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C     GRNDG=bacterial growth
C     RXNDLC,RXNDLN,RXNDLP=bacterial C,N,P loss from decomposition
C     RXNSNC,RXNSNC,RXNSNP=bacterial C,N,P loss from senescence
C     ZADDN,PADDN=nonstructural N,P used in growth
C
      WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)+GRNDG-RXNDLC-RXNSNC
      WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)+ZADDN-RXNDLN-RXNSNN
      WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)+PADDN-RXNDLP-RXNSNP
C     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,2122)'NODGR',I,J,NFZ,NZ,L
C    2,CNDLR,WTNDL(L,NZ,NY,NX),WTRTD(1,L,NZ,NY,NX)
C    3,CNDLR/CNDLI(INTYP(NZ,NY,NX))
C    2,RCNDLM,RCNDL,RMNDL,RGNDL,RGN2P
C    2,RGN2F,CGNDL,RSNDL,GRNDG,ZADDN,PADDN,RCCC,RCCN,RCCP
C    8,RDNDLC,RDNDLN,RDNDLP,RCNDLC,RDNDLX,WFR(1,L,NZ,NY,NX)
C    5,RCNSNN,RUPNF(L,NZ,NY,NX)
C    3,WTNDL(L,NZ,NY,NX),WTNDLN(L,NZ,NY,NX),WTNDLP(L,NZ,NY,NX)
C    3,ZADDN,RCNDLN,RCNSNN,RUPNF(L,NZ,NY,NX)
C    2,CPOOLN(L,NZ,NY,NX),ZPOOLN(L,NZ,NY,NX),PPOOLN(L,NZ,NY,NX)
C    5,FCNPF,TFN4(L,NZ,NY,NX),WFNRG(1,L),PSIRT(1,L,NZ,NY,NX)
C    5,CCPOLN,CZPOLN,CPPOLN,ZCPOLI,ZPPOLI
C    6,(1.0+AMAX1(ZCPOLI/ZCKI,ZPPOLI/ZPKI))
C    6,VMXO*WTNDL(L,NZ,NY,NX)*TFN4(L,NZ,NY,NX)*FCNPF*WFNRG(1,L)
2122  FORMAT(A8,5I4,60E14.6)
C     ENDIF
C
C     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND NODULES
C     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
C
C     CPOOLR,ZPOOLR,PPOOLR=root non-structural C,N,P mass
C     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C     WTRTD=root C mass
C     WTNDL=bacterial C mass
C     WTNDI=initial bacterial mass at infection
C     FXRN=rate constant for plant-bacteria nonstructural C,N,P
C        exchange
C     XFRC,XFRN,XFRC=nonstructural root-bacteria C,N,P transfer
C
      IF(CPOOLR(1,L,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WTRTD(1,L,NZ,NY,NX).GT.ZEROP2(NZ,NY,NX)
     3.AND.(ISTYP(NZ,NY,NX).NE.0
     4.OR.IDAY(10,NB1(NZ,NY,NX),NZ,NY,NX).EQ.0))THEN
      WTRTD1=WTRTD(1,L,NZ,NY,NX)
      WTNDL1=AMIN1(WTRTD(1,L,NZ,NY,NX)
     2,AMAX1(WTNDI*AREA(3,NU(NY,NX),NY,NX),WTNDL(L,NZ,NY,NX)))
      WTRTDT=WTRTD1+WTNDL1
      IF(WTRTDT.GT.ZEROP(NZ,NY,NX))THEN
      CPOOLD=(CPOOLR(1,L,NZ,NY,NX)*WTNDL1
     2-CPOOLN(L,NZ,NY,NX)*WTRTD1)/WTRTDT
      XFRC=FXRN(INTYP(NZ,NY,NX))*CPOOLD
      CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)-XFRC
      CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)+XFRC
      CPOOLT=CPOOLR(1,L,NZ,NY,NX)+CPOOLN(L,NZ,NY,NX)
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLD=(ZPOOLR(1,L,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX)
     2-ZPOOLN(L,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
      XFRN=FXRN(INTYP(NZ,NY,NX))*ZPOOLD
      PPOOLD=(PPOOLR(1,L,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX)
     2-PPOOLN(L,NZ,NY,NX)*CPOOLR(1,L,NZ,NY,NX))/CPOOLT
      XFRP=FXRN(INTYP(NZ,NY,NX))*PPOOLD
      ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)-XFRN
      PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)-XFRP
      ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)+XFRN
      PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)+XFRP
C     IF(NZ.EQ.2.AND.L.EQ.6)THEN
C     WRITE(*,2122)'NODEX',I,J,NFZ,NZ,L
C    3,XFRC,XFRN,XFRP
C    3,WTRTD(1,L,NZ,NY,NX),WTNDL1,CPOOLT
C    4,WTNDL(L,NZ,NY,NX),WTNDLN(L,NZ,NY,NX),WTNDLP(L,NZ,NY,NX)
C    2,CPOOLN(L,NZ,NY,NX),ZPOOLN(L,NZ,NY,NX),PPOOLN(L,NZ,NY,NX)
C    3,CPOOLR(1,L,NZ,NY,NX),ZPOOLR(1,L,NZ,NY,NX),PPOOLR(1,L,NZ,NY,NX)
C     ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
5400  CONTINUE
      ENDIF
C
C     TRANSFER NON-STRUCTURAL C,N,P AMONG BRANCH LEAVES
C     FROM NON-STRUCTURAL C,N,P CONCENTRATION DIFFERENCES
C     WHEN SEASONAL STORAGE C IS NOT BEING MOBILIZED
C
C     IDTHB=branch living flag: 0=alive,1=dead
C     ATRP=hourly leafout counter
C     ATRPX=number of hours required to initiate remobilization of
C        storage C for leafout
C     WTLSB=leaf+petiole mass
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
C     XFRC,XFRN,XFRP=interbranch nonstructural C,N,P transfer
C
      IF(NBR(NZ,NY,NX).GT.1)THEN
      WTPLTT=0.0
      CPOOLT=0.0
      ZPOOLT=0.0
      PPOOLT=0.0
      DO 300 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(ATRP(NB,NZ,NY,NX).GT.ATRPX(ISTYP(NZ,NY,NX)))THEN
      WTLSBZ(NB)=AMAX1(0.0,WTLSB(NB,NZ,NY,NX))
      CPOOLZ(NB)=AMAX1(0.0,CPOOL(NB,NZ,NY,NX))
      ZPOOLZ(NB)=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
      PPOOLZ(NB)=AMAX1(0.0,PPOOL(NB,NZ,NY,NX))
      WTPLTT=WTPLTT+WTLSBZ(NB)
      CPOOLT=CPOOLT+CPOOLZ(NB)
      ZPOOLT=ZPOOLT+ZPOOLZ(NB)
      PPOOLT=PPOOLT+PPOOLZ(NB)
      ENDIF
      ENDIF
300   CONTINUE
      DO 305 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(ATRP(NB,NZ,NY,NX).GT.ATRPX(ISTYP(NZ,NY,NX)))THEN
      IF(WTPLTT.GT.ZEROP(NZ,NY,NX)
     2.AND.CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      CPOOLD=CPOOLT*WTLSBZ(NB)-CPOOLZ(NB)*WTPLTT
      ZPOOLD=ZPOOLT*CPOOLZ(NB)-ZPOOLZ(NB)*CPOOLT
      PPOOLD=PPOOLT*CPOOLZ(NB)-PPOOLZ(NB)*CPOOLT
      XFRC=0.01*CPOOLD/WTPLTT
      XFRN=0.01*ZPOOLD/CPOOLT
      XFRP=0.01*PPOOLD/CPOOLT
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)+XFRC
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)+XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)+XFRP
      ENDIF
      ENDIF
      ENDIF
305   CONTINUE
      ENDIF
C
C     TRANSFER NON-STRUCTURAL C,N,P BWTWEEN ROOT AND MYCORRHIZAE
C     IN EACH ROOTED SOIL LAYER FROM NON-STRUCTURAL C,N,P
C     CONCENTRATION DIFFERENCES
C
C     MY=mycorrhizal:1=no,2=yes
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in
C        1:root,2:mycorrhizae
C     WTRTD=1:root,2:mycorrhizal C mass
C     FSNK=min ratio of branch or mycorrhizae to root for calculating
C        C transfer
C     FMYC=rate constant for root-mycorrhizal C,N,P exchange
C     XFRC,XFRN,XFRP=root-mycorrhizal nonstructural C,N,P transfer
C
      IF(MY(NZ,NY,NX).EQ.2)THEN
      DO 425 L=NU(NY,NX),NIX(NZ,NY,NX)
      CPOOLR1=AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX))
      CPOOLR2=AMAX1(0.0,CPOOLR(2,L,NZ,NY,NX))
      WTRTD1=AMAX1(0.0,WTRTD(1,L,NZ,NY,NX))
      WTRTD2=AMAX1(0.0,AMIN1(WTRTD(1,L,NZ,NY,NX),AMAX1(FSNK
     2*WTRTD(1,L,NZ,NY,NX),WTRTD(2,L,NZ,NY,NX))))
      WTPLTT=WTRTD1+WTRTD2
      IF(WTPLTT.GT.ZEROP(NZ,NY,NX))THEN
      CPOOLD=(CPOOLR1*WTRTD2-CPOOLR2*WTRTD1)/WTPLTT
      XFRC=FMYC*CPOOLD
      CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)-XFRC
      CPOOLR(2,L,NZ,NY,NX)=CPOOLR(2,L,NZ,NY,NX)+XFRC
      CPOOLT=CPOOLR1+CPOOLR2
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLR1=AMAX1(0.0,ZPOOLR(1,L,NZ,NY,NX))
      ZPOOLR2=AMAX1(0.0,ZPOOLR(2,L,NZ,NY,NX))
      PPOOLR1=AMAX1(0.0,PPOOLR(1,L,NZ,NY,NX))
      PPOOLR2=AMAX1(0.0,PPOOLR(2,L,NZ,NY,NX))
      ZPOOLD=(ZPOOLR1*CPOOLR2-ZPOOLR2*CPOOLR1)/CPOOLT
      XFRN=FMYC*ZPOOLD
      PPOOLD=(PPOOLR1*CPOOLR2-PPOOLR2*CPOOLR1)/CPOOLT
      XFRP=FMYC*PPOOLD
      ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)-XFRN
      ZPOOLR(2,L,NZ,NY,NX)=ZPOOLR(2,L,NZ,NY,NX)+XFRN
      PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)-XFRP
      PPOOLR(2,L,NZ,NY,NX)=PPOOLR(2,L,NZ,NY,NX)+XFRP
C     IF(NZ.EQ.2.AND.L.EQ.6)THEN
C     WRITE(*,9873)'MYCO',I,J,NFZ,NZ,L,XFRC,XFRN,XFRP
C    2,CPOOLR(1,L,NZ,NY,NX),CPOOLR(2,L,NZ,NY,NX)
C    3,WTRTD(1,L,NZ,NY,NX),WTRTD2,WTPLTT
C    3,ZPOOLR(1,L,NZ,NY,NX),ZPOOLR(2,L,NZ,NY,NX)
C    4,ZPOOLD,CPOOLT
C    4,PPOOLR(1,L,NZ,NY,NX),PPOOLR(2,L,NZ,NY,NX)
9873  FORMAT(A8,5I4,20E18.6)
C     ENDIF
      ENDIF
      ENDIF
425   CONTINUE
      ENDIF
C
C     TRANSFER ROOT NON-STRUCTURAL C,N,P TO SEASONAL STORAGE
C     IN PERENNIALS
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IFLGZR=remobilization flag
C     FXFR=rate constant for root-storage nonstructural C,N,P
C        exchange
C     XFRCX,XFRNX,XFRPX=nonstructural C,N,P transfer unconstrained
C        by C,N,P
C     XFRC,XFRN,XFRP=nonstructural C,N,P transfer constrained
C        by C,N,P
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     XNFH=time step from �wthr.f�
C
      IF(IFLGZR.EQ.1.AND.ISTYP(NZ,NY,NX).NE.0)THEN
      DO 5550 L=NU(NY,NX),NI(NZ,NY,NX)
      XFRCX=FXFR(IBTYP(NZ,NY,NX))
     2*AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX))*XNFH
      XFRNX=FXFR(IBTYP(NZ,NY,NX))
     2*AMAX1(0.0,ZPOOLR(1,L,NZ,NY,NX))*XNFH
      XFRPX=FXFR(IBTYP(NZ,NY,NX))
     2*AMAX1(0.0,PPOOLR(1,L,NZ,NY,NX))*XNFH
      XFRC=AMIN1(XFRCX,XFRNX/CNMN,XFRPX/CPMN)
      XFRN=AMIN1(XFRNX,XFRC*CNMX,XFRPX*CNMX/CPMN)
      XFRP=AMIN1(XFRPX,XFRC*CPMX,XFRNX*CPMX/CNMN)
      CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)-XFRC
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+XFRC
      ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)-XFRN
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+XFRN
      PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)-XFRP
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+XFRP
C     IF(NZ.EQ.1)THEN
C     WRITE(*,9871)'STORE',I,J,NFZ,NZ,L,FXFR(IBTYP(NZ,NY,NX))
C    2,XFRC,XFRN,XFRP,XFRCX,XFRNX,XFRPX
C    2,CPOOLR(1,L,NZ,NY,NX),ZPOOLR(1,L,NZ,NY,NX),PPOOLR(1,L,NZ,NY,NX)
C    3,WTRVC(NZ,NY,NX),WTRVN(NZ,NY,NX),WTRVP(NZ,NY,NX)
9871  FORMAT(A8,5I4,20E12.4)
C     ENDIF
5550  CONTINUE
      ENDIF
C
C     ROOT AND NODULE TOTALS
C
C     WTRTL,WTRTD=active,actual root C mass
C     WTRT1,WTRT2=primary,secondary root C mass in soil layer
C     TCO2T=total PFT respiration
C     RCO2A=total root respiration
C     RECO=ecosystem respiration
C     TRAU=total autotrophic respiration
C
      DO 5445 N=1,MY(NZ,NY,NX)
      DO 5450 L=NU(NY,NX),NI(NZ,NY,NX)
      WTRTL(N,L,NZ,NY,NX)=0.0
      WTRTD(N,L,NZ,NY,NX)=0.0
      DO 5460 NR=1,NRT(NZ,NY,NX)
      WTRTL(N,L,NZ,NY,NX)=WTRTL(N,L,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX)
      WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX)
     2+WTRT1(N,L,NR,NZ,NY,NX)
5460  CONTINUE
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)+RCO2A(N,L,NZ,NY,NX)
      RECO(NY,NX)=RECO(NY,NX)+RCO2A(N,L,NZ,NY,NX)
      TRAU(NY,NX)=TRAU(NY,NX)+RCO2A(N,L,NZ,NY,NX)
5450  CONTINUE
      DO 5470 NR=1,NRT(NZ,NY,NX)
      WTRTL(N,NINR(NR,NZ,NY,NX),NZ,NY,NX)
     2=WTRTL(N,NINR(NR,NZ,NY,NX),NZ,NY,NX)
     3+RTWT1(N,NR,NZ,NY,NX)
5470  CONTINUE
5445  CONTINUE
C
C     TRANSFER NON-STRUCTURAL C,N,P BETWEEN ROOT AND SHOOT
C
C     SINK STRENGTH OF ROOTS IN EACH SOIL LAYER AS A FRACTION
C     OF TOTAL SINK STRENGTH OF ROOTS IN ALL SOIL LAYERS
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     WTLS,WTRT=total PFT leaf+petiole,root C mass
C     FWTC,FWTS,FWTR=canopy,root system,root layer sink
C        weighting factor
C     RLNT,RTNT=root layer,root system sink strength
C
C     IF(ISTYP(NZ,NY,NX).EQ.0
C    2.OR.(VRNS(NB1(NZ,NY,NX),NZ,NY,NX)
C    3.GE.VRNL(NB1(NZ,NY,NX),NZ,NY,NX)
C    4.AND.VRNF(NB1(NZ,NY,NX),NZ,NY,NX)
C    5.LT.VRNX(NB1(NZ,NY,NX),NZ,NY,NX)))THEN
      IF(WTLS(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FWTC=AMIN1(1.0,0.667*WTRT(NZ,NY,NX)/WTLS(NZ,NY,NX))
      ELSE
      FWTC=1.0
      ENDIF
      IF(WTRT(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FWTS=AMIN1(1.0,WTLS(NZ,NY,NX)/(0.667*WTRT(NZ,NY,NX)))
      ELSE
      FWTS=1.0
      ENDIF
      DO 290 L=NU(NY,NX),NI(NZ,NY,NX)
      IF(RTNT(1).GT.ZEROP(NZ,NY,NX))THEN
      FWTR(L)=AMAX1(0.0,RLNT(1,L)/RTNT(1))
      ELSE
      FWTR(L)=1.0
      ENDIF
290   CONTINUE
C
C     RATE CONSTANT FOR TRANSFER IS SET FROM INPUT IN 'READQ'
C     BUT IS NOT USED FOR ANNUALS DURING GRAIN FILL
C
C     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
C
      WTLS(NZ,NY,NX)=0.0
      DO 309 NB=1,NBR(NZ,NY,NX)
      WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
309   CONTINUE
C
C     SINK STRENGTH OF BRANCHES IN EACH CANOPY AS A FRACTION
C     OF TOTAL SINK STRENGTH OF THE CANOPY
C
C     IDTHB=branch living flag: 0=alive,1=dead
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IDAY(8,=end date for setting final seed number
C     FWTB=branch sink weighting factor
C     PTSHT=rate constant for equilibrating shoot-root
C        nonstructural C concentration from PFT file
C     PTRT=allocation to leaf+petiole used to simulate phenology
C        effect on shoot-root C transfer
C     FWTC,FWTS,FWTR=canopy,root system,root layer sink
C        weighting factor
C     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,
C        1=non-woody
C     FWODB=C woody fraction in branch:0=woody,1=non-woody
C     FSNK=min ratio of branch or mycorrhizae to root for
C        calculating C transfer
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     XFRC,XFRN,XFRC=nonstructural shoot-root C,N,P transfer
C
      DO 310 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(WTLS(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FWTB(NB)=AMAX1(0.0,WTLSB(NB,NZ,NY,NX)/WTLS(NZ,NY,NX))
      ELSE
      FWTB(NB)=1.0
      ENDIF
      IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      PTSHTR=PTSHT(NZ,NY,NX)*PTRT**0.25
      ELSE
      PTSHTR=PTSHT(NZ,NY,NX)
      ENDIF
      DO 415 L=NU(NY,NX),NI(NZ,NY,NX)
      WTLSBX=WTLSB(NB,NZ,NY,NX)*FWODB(1,NZ)*FWTR(L)*FWTC
      WTRTLX=WTRTL(1,L,NZ,NY,NX)*FWODR(1,NZ)*FWTB(NB)*FWTS
      WTLSBB=AMAX1(0.0,WTLSBX,FSNK*WTRTLX)
      WTRTLR=AMAX1(0.0,WTRTLX,FSNK*WTLSBX)
      WTPLTT=WTLSBB+WTRTLR
      IF(WTPLTT.GT.ZEROP(NZ,NY,NX))THEN
      CPOOLB=AMAX1(0.0,CPOOL(NB,NZ,NY,NX)*FWTR(L))
      CPOOLS=AMAX1(0.0,CPOOLR(1,L,NZ,NY,NX)*FWTB(NB))
      CPOOLD=(CPOOLB*WTRTLR-CPOOLS*WTLSBB)/WTPLTT
      XFRC=PTSHTR*CPOOLD*XNFH
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)-XFRC
      CPOOLR(1,L,NZ,NY,NX)=CPOOLR(1,L,NZ,NY,NX)+XFRC
      CPOOLT=CPOOLS+CPOOLB
      IF(CPOOLT.GT.ZEROP(NZ,NY,NX))THEN
      ZPOOLB=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX)*FWTR(L))
      ZPOOLS=AMAX1(0.0,ZPOOLR(1,L,NZ,NY,NX)*FWTB(NB))
      ZPOOLD=(ZPOOLB*CPOOLS-ZPOOLS*CPOOLB)/CPOOLT
      XFRN=PTSHTR*ZPOOLD*XNFH
      PPOOLB=AMAX1(0.0,PPOOL(NB,NZ,NY,NX)*FWTR(L))
      PPOOLS=AMAX1(0.0,PPOOLR(1,L,NZ,NY,NX)*FWTB(NB))
      PPOOLD=(PPOOLB*CPOOLS-PPOOLS*CPOOLB)/CPOOLT
      XFRP=PTSHTR*PPOOLD*XNFH
      ELSE
      XFRN=0.0
      XFRP=0.0
      ENDIF
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-XFRN
      ZPOOLR(1,L,NZ,NY,NX)=ZPOOLR(1,L,NZ,NY,NX)+XFRN
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-XFRP
      PPOOLR(1,L,NZ,NY,NX)=PPOOLR(1,L,NZ,NY,NX)+XFRP
C     IF(NZ.EQ.2.AND.L.EQ.6)THEN
C     WRITE(*,3344)'ROOT',I,J,NFZ,NX,NY,NZ,NB,L
C    2,XFRC,XFRN,XFRP,CPOOL(NB,NZ,NY,NX)
C    3,CPOOLR(1,L,NZ,NY,NX),ZPOOL(NB,NZ,NY,NX)
C    3,ZPOOLR(1,L,NZ,NY,NX),FWTB(NB),FWTR(L)
C    3,FWTC,FWTS,WTLSBX,WTRTLX,FSNK,FDBK(NB,NZ,NY,NX)
C    4,CPOOLD,CPOOLB,WTLSBB,CPOOLS,WTRTLR
C    5,FWOOD(1,NZ),FWODB(1,NZ),WTRTL(1,L,NZ,NY,NX)
C    6,WTLSB(NB,NZ,NY,NX),RLNT(1,L),RTNT(1)
3344  FORMAT(A8,8I4,30E12.4)
C     ENDIF
      ENDIF
415   CONTINUE
      ENDIF
310   CONTINUE
C     ENDIF
C
C     TOTAL C,N,P IN EACH BRANCH
C
C     CPOOLK=total C4 nonstructural C in branch
C     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
C     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
C     CPOOL,ZPOOL,PPOOL=C3 non-structural C,N,P mass
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
C     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
C     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
C     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C
C
      DO 320 NB=1,NBR(NZ,NY,NX)
      CPOOLK(NB,NZ,NY,NX)=0.0
      DO 325 K=1,25
      CPOOLK(NB,NZ,NY,NX)=CPOOLK(NB,NZ,NY,NX)
     2+CPOOL3(K,NB,NZ,NY,NX)+CPOOL4(K,NB,NZ,NY,NX)
     3+CO2B(K,NB,NZ,NY,NX)+HCOB(K,NB,NZ,NY,NX)
325   CONTINUE
      WTSHTB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)
     2+WTSHEB(NB,NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
     3+WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)+WTGRB(NB,NZ,NY,NX)
     4+CPOOL(NB,NZ,NY,NX)+CPOOLK(NB,NZ,NY,NX)
      WTSHTN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)
     2+WTSHBN(NB,NZ,NY,NX)+WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX)
     3+WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX)
     4+ZPOOL(NB,NZ,NY,NX)
      WTSHTP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)
     2+WTSHBP(NB,NZ,NY,NX)+WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX)
     3+WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX)
     4+PPOOL(NB,NZ,NY,NX)
C     IF(NB.EQ.1)THEN
C     WRITE(*,1009)'WTSHB',I,J,NFZ,NX,NY,NZ,NB
C    2,WTSHTB(NB,NZ,NY,NX),WTLFB(NB,NZ,NY,NX)
C    2,WTSHEB(NB,NZ,NY,NX),WTSTKB(NB,NZ,NY,NX),WTRSVB(NB,NZ,NY,NX)
C    3,WTHSKB(NB,NZ,NY,NX),WTEARB(NB,NZ,NY,NX),WTGRB(NB,NZ,NY,NX)
C    4,CPOOL(NB,NZ,NY,NX),CPOOLK(NB,NZ,NY,NX)
C     ENDIF
1009  FORMAT(A8,7I4,30F16.8)
320   CONTINUE
C
C     TOTAL C,N,P IN ROOTS AND MYCORRHIZAE IN EACH SOIL LAYER
C
C     WTRTD=root C mass
C     CPOOLR=non-structural C mass in root
C
      DO 345 N=1,MY(NZ,NY,NX)
      DO 345 L=NU(NY,NX),NI(NZ,NY,NX)
      WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)+CPOOLR(N,L,NZ,NY,NX)
345   CONTINUE
      ENDIF
C
C     HARVEST, DISTURBANCE, DEATH
C
      IF(NFZ.EQ.1)THEN
C
C     TRANSFER ABOVE-GROUND C,N,P
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C
      DO 11 M=1,4
      CSNM(M,NZ,NY,NX)=0.0
      ZSNM(M,NZ,NY,NX)=0.0
      PSNM(M,NZ,NY,NX)=0.0
11    CONTINUE
      ZSNI(NZ,NY,NX)=0.0
      PSNI(NZ,NY,NX)=0.0
      IF((IHVST(NZ,I,NY,NX).GE.0.AND.J.EQ.INT(ZNOON(NY,NX))
     2.AND.IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)
     3.OR.(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6))THEN
C
C     ACCUMULATE ALL HARVESTED MATERIAL ABOVE CUTTING HEIGHT
C     ACCOUNTING FOR HARVEST EFFICIENCY ENTERED IN 'READQ'
C
C     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
C     PPX,PP=PFT population per m2,grid cell
C     THIN=thinning:fraction of population removed
C     CF=clumping factor
C     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
C          IHVST=3:reduction of clumping factor
C          IHVST=4 or 6:animal or insect biomass(g LM m-2)
C     THIN=IHVST=0-3,5: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C     ARLFC,ARLFT=leaf area of combined canopy, canopy layer
C     ARLFR,ARLFY=leaf area harvested,remaining
C     ZL=height to bottom of each canopy layer
C
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(JHVST(NZ,I,NY,NX).NE.2)THEN
      PPX(NZ,NY,NX)=PPX(NZ,NY,NX)*(1.0-THIN(NZ,I,NY,NX))
      PP(NZ,NY,NX)=PP(NZ,NY,NX)*(1.0-THIN(NZ,I,NY,NX))
      ELSE
      PPX(NZ,NY,NX)=PPZ(NZ,NY,NX)
      PP(NZ,NY,NX)=PPX(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      ENDIF
      IF(IHVST(NZ,I,NY,NX).EQ.3)THEN
      CF(NZ,NY,NX)=CF(NZ,NY,NX)*HVST(NZ,I,NY,NX)
      ENDIF
      IF(IHVST(NZ,I,NY,NX).LE.2.AND.HVST(NZ,I,NY,NX).LT.0.0)THEN
      ARLFY=(1.0-ABS(HVST(NZ,I,NY,NX)))*ARLFC(NY,NX)
      ARLFR=0.0
      DO 9875 L=1,JC
      IF(ZL(L,NY,NX).GT.ZL(L-1,NY,NX)
     2.AND.ARLFT(L,NY,NX).GT.ZEROS(NY,NX)
     3.AND.ARLFR.LT.ARLFY)THEN
      IF(ARLFR+ARLFT(L,NY,NX).GT.ARLFY)THEN
      HVST(NZ,I,NY,NX)=ZL(L-1,NY,NX)+((ARLFY-ARLFR)
     2/ARLFT(L,NY,NX))*(ZL(L,NY,NX)-ZL(L-1,NY,NX))
      ELSE
      HVST(NZ,I,NY,NX)=0.0
      ENDIF
      ARLFR=ARLFR+ARLFT(L,NY,NX)
      ENDIF
C     IF(NZ.EQ.2)THEN
C     WRITE(*,6544)'HVST',I,J,NFZ,L,NZ,IHVST(NZ,I,NY,NX),ARLFC(NY,NX)
C    2,ARLFT(L,NY,NX),ARLFY,ARLFR,ZL(L,NY,NX),ZL(L-1,NY,NX)
C    3,ARLFV(L,NZ,NY,NX)
6544  FORMAT(A8,6I4,20E12.4)
C     ENDIF
9875  CONTINUE
      ENDIF
      WHVSTT=0.0
      WHVSLF=0.0
      WHVHSH=0.0
      WHVEAH=0.0
      WHVGRH=0.0
      WHVSCP=0.0
      WHVSTH=0.0
      WHVRVH=0.0
      ELSE
C
C     GRAZING REMOVAL
C
C     WTSHTA=average biomass in landscape grazing section
C     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
C          IHVST=3:reduction of clumping factor
C          IHVST=4 or 6:animal or insect biomass(g LM m-2)
C     THIN=IHVST=0-3: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C     WHVSTT=total phytomass grazed or removed
C     TFN3=temperature function for canopy growth
C     CCPOLP=nonstructural C concentration in canopy
C     CCPLNP=nonstructural C concentration in canopy nodules
C
      IF(WTSHTA(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSTT=HVST(NZ,I,NY,NX)*THIN(NZ,I,NY,NX)*0.45/24.0
     2*AREA(3,NU(NY,NX),NY,NX)*WTSHT(NZ,NY,NX)/WTSHTA(NZ,NY,NX)
      ELSE
      WHVSTT=0.0
      ENDIF
      IF(IHVST(NZ,I,NY,NX).EQ.6)THEN
      WHVSTT=WHVSTT*TFN3(NZ,NY,NX)
      ENDIF
      CCPOLX=CCPOLP(NZ,NY,NX)/(1.0+CCPOLP(NZ,NY,NX))
      CCPLNX=CCPLNP(NZ,NY,NX)/(1.0+CCPLNP(NZ,NY,NX))
C
C     LEAF,BACTERIA GRAZED,REMOVED
C
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C           leaf,non-foliar,woody, standing dead removed from PFT
C     EHVST(2,1,EHVST(2,2,EHVST(2,3,EHVST(2,4=fraction of
C           leaf,non-foliar,woody, standing dead removed from
C        ecosystem
C     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
C     WTLF=PFT leaf C mass
C     WHVXXX=grazing requirement unmet by leaf
C
      WHVSLX=WHVSTT*EHVST(1,1,NZ,I,NY,NX)
      WHVSLY=AMIN1(WTLF(NZ,NY,NX),WHVSLX)
      WHVSLF=WHVSLY*(1.0-CCPOLX)
      WHVSCL=WHVSLY*CCPOLX
      WHVSNL=WHVSLY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVSLX-WHVSLY)
      WHVSSX=WHVSTT*EHVST(1,2,NZ,I,NY,NX)
C
C     OTHER NON-FOLIAR GRAZED,REMOVED
C
C     WTSHE,WTHSK,WTEAR,WTGR=PFT petiole,husk,ear,grain C mass
C     WHVSH*,WHVHS*,WHVEA*,WHVGR*,WHVSC*
C        =petiole,husk,ear,grain,nonstructural C removed
C     WHVXXX=grazing requirement unmet by non-foliar removal
C
      WTSHTT=WTSHE(NZ,NY,NX)+WTHSK(NZ,NY,NX)+WTEAR(NZ,NY,NX)
     2+WTGR(NZ,NY,NX)
      IF(WTSHTT.GT.ZEROP(NZ,NY,NX))THEN
      WHVSHX=WHVSSX*WTSHE(NZ,NY,NX)/WTSHTT+WHVXXX
      WHVSHY=AMIN1(WTSHE(NZ,NY,NX),WHVSHX)
      WHVSHH=WHVSHY*(1.0-CCPOLX)
      WHVSCS=WHVSHY*CCPOLX
      WHVSNS=WHVSHY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVSHX-WHVSHY)
      WHVHSX=WHVSSX*WTHSK(NZ,NY,NX)/WTSHTT+WHVXXX
      WHVHSY=AMIN1(WTHSK(NZ,NY,NX),WHVHSX)
      WHVHSH=WHVHSY
      WHVXXX=AMAX1(0.0,WHVHSX-WHVHSY)
      WHVEAX=WHVSSX*WTEAR(NZ,NY,NX)/WTSHTT+WHVXXX
      WHVEAY=AMIN1(WTEAR(NZ,NY,NX),WHVEAX)
      WHVEAH=WHVEAY
      WHVXXX=AMAX1(0.0,WHVEAX-WHVEAY)
      WHVGRX=WHVSSX*WTGR(NZ,NY,NX)/WTSHTT+WHVXXX
      WHVGRY=AMIN1(WTGR(NZ,NY,NX),WHVGRX)
      WHVGRH=WHVGRY
      WHVXXX=AMAX1(0.0,WHVGRX-WHVGRY)
      ELSE
      WHVSHH=0.0
      WHVSCS=0.0
      WHVSNS=0.0
      WHVHSH=0.0
      WHVEAH=0.0
      WHVGRH=0.0
      WHVXXX=WHVXXX+WHVSSX
      ENDIF
      WHVSCP=WHVSCL+WHVSCS
      WHVSNP=WHVSNL+WHVSNS
      WHVSKX=WHVSTT*EHVST(1,3,NZ,I,NY,NX)
C
C     STALK GRAZED, REMOVED
C
C     WTSTK,WTRSV=stalk,reserve C mass
C     WHVST*,WHVRV*=stalk,reserve C removed
C     WHVXXX=grazing requirement unmet by stalk,reserve
C
      WTSTKT=WTSTK(NZ,NY,NX)+WTRSV(NZ,NY,NX)
      IF(WTSTKT.GT.WHVSKX+WHVXXX)THEN
      WHVSTX=WHVSKX*WTSTK(NZ,NY,NX)/WTSTKT+WHVXXX
      WHVSTY=AMIN1(WTSTK(NZ,NY,NX),WHVSTX)
      WHVSTH=WHVSTY
      WHVXXX=AMAX1(0.0,WHVSTX-WHVSTY)
      WHVRVX=WHVSKX*WTRSV(NZ,NY,NX)/WTSTKT+WHVXXX
      WHVRVY=AMIN1(WTRSV(NZ,NY,NX),WHVRVX)
      WHVRVH=WHVRVY
      WHVXXX=AMAX1(0.0,WHVRVX-WHVRVY)
      ELSE
      WHVSTH=0.0
      WHVRVH=0.0
      WHVXXX=AMAX1(0.0,WHVSKX)
C
C     ALLOCATE UNMET DEMAND FOR GRAZING TO LEAF,PETIOLE,HUSK
C     EAR,GRAIN
C
C     WHVSL*,WHVSC*,WHVSN=leaf,nonstructural,bacteria removed
C     WHVSH*,WHVHS,WHVEA,WHVGR,WHVSC
C        =petiole,husk,ear,grain,nonstructural C removed
C
      IF(WHVXXX.GT.0.0)THEN
      WHVSLY=AMIN1(WTLF(NZ,NY,NX)-WHVSLF-WHVSCL,WHVXXX)
      WHVSLF=WHVSLF+WHVSLY*(1.0-CCPOLX)
      WHVSCL=WHVSCL+WHVSLY*CCPOLX
      WHVSNL=WHVSNL+WHVSLY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVXXX-WHVSLY)
      IF(WTSHTT.GT.ZEROP(NZ,NY,NX))THEN
      WHVSHX=WHVXXX*WTSHE(NZ,NY,NX)/WTSHTT
      WHVSHY=AMIN1(WTSHE(NZ,NY,NX),WHVSHX)
      WHVSHH=WHVSHH+WHVSHY*(1.0-CCPOLX)
      WHVSCS=WHVSCS+WHVSHY*CCPOLX
      WHVSNS=WHVSNS+WHVSHY*CCPLNX
      WHVXXX=AMAX1(0.0,WHVXXX-WHVSHY)
      WHVHSX=WHVXXX*WTHSK(NZ,NY,NX)/WTSHTT
      WHVHSY=AMIN1(WTHSK(NZ,NY,NX),WHVHSX)
      WHVHSH=WHVHSH+WHVHSY
      WHVXXX=AMAX1(0.0,WHVXXX-WHVHSY)
      WHVEAX=WHVXXX*WTEAR(NZ,NY,NX)/WTSHTT
      WHVEAY=AMIN1(WTEAR(NZ,NY,NX),WHVEAX)
      WHVEAH=WHVEAH+WHVEAY
      WHVXXX=AMAX1(0.0,WHVEAX-WHVEAY)
      WHVGRX=WHVXXX*WTGR(NZ,NY,NX)/WTSHTT
      WHVGRY=AMIN1(WTGR(NZ,NY,NX),WHVGRX)
      WHVGRH=WHVGRH+WHVGRY
      WHVXXX=AMAX1(0.0,WHVGRX-WHVGRY)
      ENDIF
      ENDIF
      ENDIF
C
C     ALL HARVEST REMOVALS
C
C     WGLFBL=branch leaf C mass in canopy layer
C
      DO 9860 NB=1,NBR(NZ,NY,NX)
      DO 9860 L=1,JC
      DO 9860 K=0,25
      WGLFBL(L,NB,NZ,NY,NX)=0.0
9860  CONTINUE
      DO 9870 NB=1,NBR(NZ,NY,NX)
      DO 9870 L=1,JC
      DO 9870 K=0,25
      WGLFBL(L,NB,NZ,NY,NX)=WGLFBL(L,NB,NZ,NY,NX)
     2+WGLFL(L,K,NB,NZ,NY,NX)
9870  CONTINUE
      ENDIF
C
C     HARVEST REMOVAL FROM TOP TO BOTTOM OF CANOPY
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     ZL=height to bottom of each canopy layer
C     FHGT=fraction of canopy layer height not harvested
C     FHVST=fraction of canopy layer mass not harvested
C     THIN=IHVST=0-3,5: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C        leaf,non-foliar,woody, standing dead removed from PFT
C
      DO 9865 L=JC,1,-1
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(IHVST(NZ,I,NY,NX).NE.3)THEN
      IF(ZL(L,NY,NX).GT.ZL(L-1,NY,NX))THEN
      FHGT=AMAX1(0.0,AMIN1(1.0,1.0-((ZL(L,NY,NX))
     2-HVST(NZ,I,NY,NX))/(ZL(L,NY,NX)-ZL(L-1,NY,NX))))
      ELSE
      FHGT=1.0
      ENDIF
      ELSE
      FHGT=0.0
      ENDIF
      IF(THIN(NZ,I,NY,NX).EQ.0.0)THEN
      FHVST=AMAX1(0.0,1.0-(1.0-FHGT)*EHVST(1,1,NZ,I,NY,NX))
      FHVSH=FHVST
      ELSE
      FHVST=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      FHVSH=1.0-(1.0-FHGT)*EHVST(1,1,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
      ELSE
      FHVSH=FHVST
      ENDIF
      ENDIF
      ELSE
      FHVST=0.0
      FHVSH=0.0
      ENDIF
C
C     CUT LEAVES AT HARVESTED NODES AND LAYERS
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     WTLF=PFT leaf C mass
C     WGLFBL=branch leaf C mass in canopy layer
C     WHVBSL,WHVSLF=layer,total leaf C mass removed
C     WGLFL=leaf node C in canopy layer
C     FHVST=fraction of leaf node mass not harvested
C
      DO 9855 NB=1,NBR(NZ,NY,NX)
      IF((IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)
     2.AND.WTLF(NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
      WHVSBL=WHVSLF*AMAX1(0.0,WGLFBL(L,NB,NZ,NY,NX))/WTLF(NZ,NY,NX)
      ELSE
      WHVSBL=0.0
      ENDIF
      DO 9845 K=25,0,-1
      IF((IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)
     2.OR.WHVSBL.GT.0.0)THEN
      IF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      IF(WGLFL(L,K,NB,NZ,NY,NX).GT.WHVSBL)THEN
      FHVST=AMAX1(0.0,AMIN1(1.0,(WGLFL(L,K,NB,NZ,NY,NX)-WHVSBL)
     2/WGLFL(L,K,NB,NZ,NY,NX)))
      FHVSH=FHVST
      ELSE
      FHVST=1.0
      FHVSH=1.0
      ENDIF
      ENDIF
C
C     HARVESTED LEAF AREA, C, N, P
C
C     FHVST=fraction of leaf node mass not harvested
C     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
C     ARLFL,ARSTK=leaf,stalk node area in canopy layer
C     WTHTH1,WTHNH1,WTHPH1=harvested leaf C,N,P
C     WTHTX1,WTHNX1,WTHPX1=harvested leaf C,N,P to litter
C     WTHTH3,WTHNH3,WTHPH3=harvested woody C,N,P
C     WTHTX3,WTHNX3,WTHPX3=harvested woody C,N,P to litter
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
C
      WHVSBL=WHVSBL-(1.0-FHVST)*WGLFL(L,K,NB,NZ,NY,NX)
      WTHTH1=WTHTH1+(1.0-FHVSH)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(1,NZ)
      WTHNH1=WTHNH1+(1.0-FHVSH)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(1,NZ)
      WTHPH1=WTHPH1+(1.0-FHVSH)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(1,NZ)
      WTHTX1=WTHTX1+(FHVSH-FHVST)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(1,NZ)
      WTHNX1=WTHNX1+(FHVSH-FHVST)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(1,NZ)
      WTHPX1=WTHPX1+(FHVSH-FHVST)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(1,NZ)
      WTHTH3=WTHTH3+(1.0-FHVSH)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(0,NZ)
      WTHNH3=WTHNH3+(1.0-FHVSH)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(0,NZ)
      WTHPH3=WTHPH3+(1.0-FHVSH)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(0,NZ)
      WTHTX3=WTHTX3+(FHVSH-FHVST)*WGLFL(L,K,NB,NZ,NY,NX)*FWODB(0,NZ)
      WTHNX3=WTHNX3+(FHVSH-FHVST)*WGLFLN(L,K,NB,NZ,NY,NX)*FWODLN(0,NZ)
      WTHPX3=WTHPX3+(FHVSH-FHVST)*WGLFLP(L,K,NB,NZ,NY,NX)*FWODLP(0,NZ)
C
C     REMAINING LEAF C,N,P AND AREA
C
      WGLFL(L,K,NB,NZ,NY,NX)=FHVST*WGLFL(L,K,NB,NZ,NY,NX)
      WGLFLN(L,K,NB,NZ,NY,NX)=FHVST*WGLFLN(L,K,NB,NZ,NY,NX)
      WGLFLP(L,K,NB,NZ,NY,NX)=FHVST*WGLFLP(L,K,NB,NZ,NY,NX)
      ARLFL(L,K,NB,NZ,NY,NX)=FHVST*ARLFL(L,K,NB,NZ,NY,NX)
      IF(K.EQ.1)THEN
      ARSTK(L,NB,NZ,NY,NX)=FHVST*ARSTK(L,NB,NZ,NY,NX)
      ENDIF
      ENDIF
C     IF(I.EQ.262.AND.K.EQ.5)THEN
C     WRITE(*,6543)'GRAZ',I,J,NZ,NB,K,L,IHVST(NZ,I,NY,NX)
C    2,ZL(L,NY,NX),ZL(L-1,NY,NX),HVST(NZ,I,NY,NX),FHVST,FHVSH
C    5,WGLFBL(L,NB,NZ,NY,NX),WTLF(NZ,NY,NX),CPOOLP(NZ,NY,NX)
C    6,ARLFL(L,K,NB,NZ,NY,NX),WGLF(K,NB,NZ,NY,NX),ARLF(K,NB,NZ,NY,NX)
C    7,HTNODE(K,NB,NZ,NY,NX)
C    7,WTSHTA(NZ,NY,NX),WHVSBL,WHVSTT,WHVSLF,WHVSHH
C    3,WHVHSH,WHVEAH,WHVGRH,WHVSCP,WHVSTH,WHVRVH,WHVXXX
C    4,WTSHTT,WHVSSX,CCPOLX
6543  FORMAT(A8,7I4,30E12.4)
C     ENDIF
9845  CONTINUE
9855  CONTINUE
      ARLFV(L,NZ,NY,NX)=0.0
      WGLFV(L,NZ,NY,NX)=0.0
      ARSTV(L,NZ,NY,NX)=ARSTV(L,NZ,NY,NX)*FHVST
9865  CONTINUE
      DO 9835 NB=1,NBR(NZ,NY,NX)
      IF(JHVST(NZ,I,NY,NX).NE.0)THEN
      IDTHB(NB,NZ,NY,NX)=1
      ENDIF
      CPOOLG=0.0
      ZPOOLG=0.0
      PPOOLG=0.0
      CPOLNG=0.0
      ZPOLNG=0.0
      PPOLNG=0.0
      WTNDG=0.0
      WTNDNG=0.0
      WTNDPG=0.0
      WGLFGX=0.0
      WGSHGX=0.0
      WGLFGY=0.0
      WGSHGY=0.0
      DO 9825 K=0,25
      ARLFG=0.0
      WGLFG=0.0
      WGLFNG=0.0
      WGLFPG=0.0
C
C     ACCUMULATE REMAINING LEAF AREA, C, N, P
C
C     WGLFL,WGLFLN,WGLFLP=leaf node C,N,P in canopy layer
C     ARLFL,ARLFV=leaf node,total area in canopy layer
C
      DO 9815 L=1,JC
      ARLFG=ARLFG+ARLFL(L,K,NB,NZ,NY,NX)
      WGLFG=WGLFG+WGLFL(L,K,NB,NZ,NY,NX)
      WGLFNG=WGLFNG+WGLFLN(L,K,NB,NZ,NY,NX)
      WGLFPG=WGLFPG+WGLFLP(L,K,NB,NZ,NY,NX)
      ARLFV(L,NZ,NY,NX)=ARLFV(L,NZ,NY,NX)+ARLFL(L,K,NB,NZ,NY,NX)
      WGLFV(L,NZ,NY,NX)=WGLFV(L,NZ,NY,NX)+WGLFL(L,K,NB,NZ,NY,NX)
9815  CONTINUE
C
C     CUT STALK AT HARVESTED NODES AND LAYERS
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     WGLF=leaf node C mass
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C        leaf,non-foliar,woody, standing dead removed from PFT
C     FHVSTK=fraction of internode layer mass not harvested
C     THIN=IHVST=0-3,5: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(WGLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.EHVST(1,1,NZ,I,NY,NX).GT.0.0)THEN
      FHVSTK(K)=AMAX1(0.0,AMIN1(1.0,(1.0-(1.0-AMAX1(0.0,WGLFG)
     2/WGLF(K,NB,NZ,NY,NX))*EHVST(1,2,NZ,I,NY,NX)
     3/EHVST(1,1,NZ,I,NY,NX))))
      FHVSHK(K)=FHVSTK(K)
      ELSE
      IF(THIN(NZ,I,NY,NX).EQ.0.0)THEN
      FHVSTK(K)=1.0-EHVST(1,2,NZ,I,NY,NX)
      FHVSHK(K)=FHVSTK(K)
      ELSE
      FHVSTK(K)=1.0-THIN(NZ,I,NY,NX)
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      FHVSHK(K)=1.0-EHVST(1,2,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
      ELSE
      FHVSHK(K)=FHVSTK(K)
      ENDIF
      ENDIF
      ENDIF
      ELSE
      FHVSTK(K)=0.0
      FHVSHK(K)=0.0
      ENDIF
C
C     ACCUMULATE REMAINING BRANCH LEAF AREA, C, N, P
C
C     WGLF=leaf node C mass
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     ARLFB,ARLF=branch,node leaf area
C     WSLF=leaf protein mass
C
      WGLFGY=WGLFGY+WGLF(K,NB,NZ,NY,NX)
      WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)
     2-WGLF(K,NB,NZ,NY,NX)+WGLFG
      WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)
     2-WGLFN(K,NB,NZ,NY,NX)+WGLFNG
      WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)
     2-WGLFP(K,NB,NZ,NY,NX)+WGLFPG
      ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)-ARLF(K,NB,NZ,NY,NX)+ARLFG
      IF(ARLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WSLF(K,NB,NZ,NY,NX)=WSLF(K,NB,NZ,NY,NX)
     2*ARLFG/ARLF(K,NB,NZ,NY,NX)
      ELSE
      WSLF(K,NB,NZ,NY,NX)=0.0
      ENDIF
      ARLF(K,NB,NZ,NY,NX)=ARLFG
      WGLF(K,NB,NZ,NY,NX)=WGLFG
      WGLFN(K,NB,NZ,NY,NX)=WGLFNG
      WGLFP(K,NB,NZ,NY,NX)=WGLFPG
      WGLFGX=WGLFGX+WGLF(K,NB,NZ,NY,NX)
9825  CONTINUE
C
C     CUT SHEATHS OR PETIOLES AND STALKS HARVESTED NODES AND LAYERS
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     WTSHE,WTSHEB=PFT,branch petiole C mass
C     WHVSBS,WHVSHH=branch, PFT petiole C mass removed
C     HTNODE=internode length
C     HTSTKX=internode length removed
C
      HTSTKX=0.0
      IF((IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)
     2.AND.WTSHE(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSBS=WHVSHH*WTSHEB(NB,NZ,NY,NX)/WTSHE(NZ,NY,NX)
      ELSE
      WHVSBS=0.0
      ENDIF
      DO 9805 K=25,0,-1
112   FORMAT(A8,8I4,12E12.4)
      IF(HTNODE(K,NB,NZ,NY,NX).GT.0.0)
     2HTSTKX=AMAX1(HTSTKX,HTNODE(K,NB,NZ,NY,NX))
C     WRITE(*,112)'VSTG',I,J,NX,NY,NZ,NB,K,IDTHB(NB,NZ,NY,NX)
C    2,VSTG(NB,NZ,NY,NX),FHVSTK(K),HTSTKX,HTNODE(K,NB,NZ,NY,NX)
C    3,HVST(NZ,I,NY,NX)
C
C     HARVESTED SHEATH OR PETIOLE C,N,P
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     WHVSBS=branch petiole C mass removed
C     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
C     FHVSTK=fraction of internode layer mass not harvested
C     WTHTH2,WTHNH2,WTHPH2=harvested petiole C,N,P
C     WTHTX2,WTHNX2,WTHPX2=harvested petiole C,N,P to litter
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     FWODLN,FWODLP=N,P woody fraction in leaf:0=woody,1=non-woody
C     HTSHE,HTNODE=petiole,internode length
C
      IF((IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)
     2.OR.WHVSBS.GT.0.0)THEN
      IF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      IF(WGSHE(K,NB,NZ,NY,NX).GT.WHVSBS)THEN
      FHVSTK(K)=AMAX1(0.0,AMIN1(1.0,(WGSHE(K,NB,NZ,NY,NX)-WHVSBS)
     2/WGSHE(K,NB,NZ,NY,NX)))
      FHVSHK(K)=FHVSTK(K)
      ELSE
      FHVSTK(K)=0.0
      FHVSHK(K)=0.0
      ENDIF
      ENDIF
      WHVSBS=WHVSBS-(1.0-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)
      WTHTH2=WTHTH2+(1.0-FHVSHK(K))*WGSHE(K,NB,NZ,NY,NX)*FWODB(1,NZ)
      WTHNH2=WTHNH2+(1.0-FHVSHK(K))*WGSHN(K,NB,NZ,NY,NX)*FWODSN(1,NZ)
      WTHPH2=WTHPH2+(1.0-FHVSHK(K))*WGSHP(K,NB,NZ,NY,NX)*FWODSP(1,NZ)
      WTHTX2=WTHTX2+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)
     2*FWODB(1,NZ)
      WTHNX2=WTHNX2+(FHVSHK(K)-FHVSTK(K))*WGSHN(K,NB,NZ,NY,NX)
     2*FWODSN(1,NZ)
      WTHPX2=WTHPX2+(FHVSHK(K)-FHVSTK(K))*WGSHP(K,NB,NZ,NY,NX)
     2*FWODSP(1,NZ)
      WTHTH3=WTHTH3+(1.0-FHVSHK(K))*WGSHE(K,NB,NZ,NY,NX)*FWODB(0,NZ)
      WTHNH3=WTHNH3+(1.0-FHVSHK(K))*WGSHN(K,NB,NZ,NY,NX)*FWODSN(0,NZ)
      WTHPH3=WTHPH3+(1.0-FHVSHK(K))*WGSHP(K,NB,NZ,NY,NX)*FWODSP(0,NZ)
      WTHTX3=WTHTX3+(FHVSHK(K)-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)
     2*FWODB(0,NZ)
      WTHNX3=WTHNX3+(FHVSHK(K)-FHVSTK(K))*WGSHN(K,NB,NZ,NY,NX)
     2*FWODSN(0,NZ)
      WTHPX3=WTHPX3+(FHVSHK(K)-FHVSTK(K))*WGSHP(K,NB,NZ,NY,NX)
     2*FWODSP(0,NZ)
C
C     ACCUMULATE REMAINING SHEATH OR PETIOLE C,N,P AND LENGTH
C
C     WGSHE=petiole node C mass
C     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
C     HTSHE=node petiole height
C     WSSHE=petiole protein mass
C
      WGSHGY=WGSHGY+WGSHE(K,NB,NZ,NY,NX)
      WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX)
     2-(1.0-FHVSTK(K))*WGSHE(K,NB,NZ,NY,NX)
      WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX)
     2-(1.0-FHVSTK(K))*WGSHN(K,NB,NZ,NY,NX)
      WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX)
     2-(1.0-FHVSTK(K))*WGSHP(K,NB,NZ,NY,NX)
      WGSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*WGSHE(K,NB,NZ,NY,NX)
      WSSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*WSSHE(K,NB,NZ,NY,NX)
      WGSHN(K,NB,NZ,NY,NX)=FHVSTK(K)*WGSHN(K,NB,NZ,NY,NX)
      WGSHP(K,NB,NZ,NY,NX)=FHVSTK(K)*WGSHP(K,NB,NZ,NY,NX)
      IF(IHVST(NZ,I,NY,NX).LE.2
     2.AND.HTSHE(K,NB,NZ,NY,NX).GT.0.0)THEN
      FHGT=AMAX1(0.0,AMIN1(1.0,(HTNODE(K,NB,NZ,NY,NX)
     2+HTSHE(K,NB,NZ,NY,NX)-HVST(NZ,I,NY,NX))/HTSHE(K,NB,NZ,NY,NX)))
      HTSHE(K,NB,NZ,NY,NX)=(1.0-FHGT)*HTSHE(K,NB,NZ,NY,NX)
      ELSE
      HTSHE(K,NB,NZ,NY,NX)=FHVSTK(K)*HTSHE(K,NB,NZ,NY,NX)
      ENDIF
      WGSHGX=WGSHGX+WGSHE(K,NB,NZ,NY,NX)
C     IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
C     IF(HTNODE(K,NB,NZ,NY,NX).GT.HVST(NZ,I,NY,NX)
C    2.OR.IHVST(NZ,I,NY,NX).EQ.3)THEN
C     IF(FHVSTK(K).EQ.0.0.AND.K.GT.0)THEN
C     IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
C     VSTG(NB,NZ,NY,NX)=AMAX1(0.0,VSTG(NB,NZ,NY,NX)-1.0)
C     ELSE
C     VSTG(NB,NZ,NY,NX)=AMAX1(0.0,VSTG(NB,NZ,NY,NX)-0.04)
C     ENDIF
C     ENDIF
C     ENDIF
C     ENDIF
      ENDIF
9805  CONTINUE
C
C     CUT NON-STRUCTURAL C,N,P IN HARVESTED BRANCHES
C
C     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     FHVST=fraction of leaf+petiole node mass not harvested
C     CPOOLG,ZPOOLG,PPOOLG=branch non-structural C,N,P mass
C        after harvest
C     CPOLNG,ZPOLNG,PPOLNG=nonstructural C,N,P in bacteria
C        after harvest
C     WTNDG,WTNDNG,WTNDPG=bacterial C,N,P mass after harvest
C     WTLS,WTLSB=total,branch PFT leaf+petiole C mass
C     WHVSC*=nonstructural C removed
C
      CPOOLX=AMAX1(0.0,CPOOL(NB,NZ,NY,NX))
      ZPOOLX=AMAX1(0.0,ZPOOL(NB,NZ,NY,NX))
      PPOOLX=AMAX1(0.0,PPOOL(NB,NZ,NY,NX))
      CPOLNX=AMAX1(0.0,CPOLNB(NB,NZ,NY,NX))
      ZPOLNX=AMAX1(0.0,ZPOLNB(NB,NZ,NY,NX))
      PPOLNX=AMAX1(0.0,PPOLNB(NB,NZ,NY,NX))
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(WGLFGY+WGSHGY.GT.ZEROP(NZ,NY,NX))THEN
      FHVST=AMAX1(0.0,AMIN1(1.0,(WGLFGX+WGSHGX)
     2/(WGLFGY+WGSHGY)))
      CPOOLG=CPOOLX*FHVST
      ZPOOLG=ZPOOLX*FHVST
      PPOOLG=PPOOLX*FHVST
      CPOLNG=CPOLNX*FHVST
      ZPOLNG=ZPOLNX*FHVST
      PPOLNG=PPOLNX*FHVST
      WTNDG=WTNDB(NB,NZ,NY,NX)*FHVST
      WTNDNG=WTNDBN(NB,NZ,NY,NX)*FHVST
      WTNDPG=WTNDBP(NB,NZ,NY,NX)*FHVST
      ELSE
      CPOOLG=0.0
      ZPOOLG=0.0
      PPOOLG=0.0
      CPOLNG=0.0
      ZPOLNG=0.0
      PPOLNG=0.0
      WTNDG=0.0
      WTNDNG=0.0
      WTNDPG=0.0
      ENDIF
      ELSE
      IF(WTLS(NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
      WTLSBX=AMAX1(0.0,WTLSB(NB,NZ,NY,NX))
      IF(CPOOL(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSCX=AMAX1(0.0,WHVSCP)*WTLSBX/WTLS(NZ,NY,NX)
      CPOOLG=AMAX1(0.0,CPOOLX-WHVSCX)
      ZPOOLG=AMAX1(0.0,ZPOOLX-WHVSCX*ZPOOLX/CPOOL(NB,NZ,NY,NX))
      PPOOLG=AMAX1(0.0,PPOOLX-WHVSCX*PPOOLX/CPOOL(NB,NZ,NY,NX))
      ELSE
      CPOOLG=0.0
      ZPOOLG=0.0
      PPOOLG=0.0
      ENDIF
      IF(CPOLNB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSNX=AMAX1(0.0,WHVSNP)*WTLSBX/WTLS(NZ,NY,NX)
      CPOLNG=AMAX1(0.0,CPOLNX-WHVSNX)
      ZPOLNG=AMAX1(0.0,ZPOLNX-WHVSNX*ZPOLNX/CPOLNB(NB,NZ,NY,NX))
      PPOLNG=AMAX1(0.0,PPOLNX-WHVSNX*PPOLNX/CPOLNB(NB,NZ,NY,NX))
      WTNDG=WTNDB(NB,NZ,NY,NX)*(1.0-WHVSNX/CPOLNX)
      WTNDNG=WTNDBN(NB,NZ,NY,NX)*(1.0-WHVSNX/CPOLNX)
      WTNDPG=WTNDBP(NB,NZ,NY,NX)*(1.0-WHVSNX/CPOLNX)
      ELSE
      CPOLNG=0.0
      ZPOLNG=0.0
      PPOLNG=0.0
      WTNDG=0.0
      WTNDNG=0.0
      WTNDPG=0.0
      ENDIF
      ELSE
      CPOOLG=0.0
      ZPOOLG=0.0
      PPOOLG=0.0
      CPOLNG=0.0
      ZPOLNG=0.0
      PPOLNG=0.0
      WTNDG=0.0
      WTNDNG=0.0
      WTNDPG=0.0
      ENDIF
      ENDIF
C
C     HARVESTED NON-STRUCTURAL C, N, P
C
C     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
C
      WTHTH0=WTHTH0+CPOOLX-CPOOLG+CPOLNX-CPOLNG
      WTHNH0=WTHNH0+ZPOOLX-ZPOOLG+ZPOLNX-ZPOLNG
      WTHPH0=WTHPH0+PPOOLX-PPOOLG+PPOLNX-PPOLNG
      WTHTH0=WTHTH0+WTNDB(NB,NZ,NY,NX)-WTNDG
      WTHNH0=WTHNH0+WTNDBN(NB,NZ,NY,NX)-WTNDNG
      WTHPH0=WTHPH0+WTNDBP(NB,NZ,NY,NX)-WTNDPG
C
C     REMAINING NON-STRUCTURAL C, N, P
C
C     CPOOL,ZPOOL,PPOOL=branch non-structural C,N,P mass
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C
      CPOOL(NB,NZ,NY,NX)=CPOOLG
      ZPOOL(NB,NZ,NY,NX)=ZPOOLG
      PPOOL(NB,NZ,NY,NX)=PPOOLG
      CPOLNB(NB,NZ,NY,NX)=CPOLNG
      ZPOLNB(NB,NZ,NY,NX)=ZPOLNG
      PPOLNB(NB,NZ,NY,NX)=PPOLNG
      WTNDB(NB,NZ,NY,NX)=WTNDG
      WTNDBN(NB,NZ,NY,NX)=WTNDNG
      WTNDBP(NB,NZ,NY,NX)=WTNDPG
C
C     REMOVE C4 NON-STRUCTURAL C
C
C     ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
C     FHVST4=fraction of nonstructural mass not harvested
C     CPOOLG=branch non-structural C mass after harvest
C     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
C     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
C     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
C
      IF(ICTYP(NZ,NY,NX).EQ.4.AND.CPOOLX.GT.ZEROP(NZ,NY,NX))THEN
      FHVST4=CPOOLG/CPOOLX
      DO 9810 K=1,25
      WTHTH0=WTHTH0+(1.0-FHVST4)*CPOOL3(K,NB,NZ,NY,NX)
      WTHTH0=WTHTH0+(1.0-FHVST4)*CPOOL4(K,NB,NZ,NY,NX)
      WTHTH0=WTHTH0+(1.0-FHVST4)*CO2B(K,NB,NZ,NY,NX)
      WTHTH0=WTHTH0+(1.0-FHVST4)*HCOB(K,NB,NZ,NY,NX)
      CPOOL3(K,NB,NZ,NY,NX)=FHVST4*CPOOL3(K,NB,NZ,NY,NX)
      CPOOL4(K,NB,NZ,NY,NX)=FHVST4*CPOOL4(K,NB,NZ,NY,NX)
      CO2B(K,NB,NZ,NY,NX)=FHVST4*CO2B(K,NB,NZ,NY,NX)
      HCOB(K,NB,NZ,NY,NX)=FHVST4*HCOB(K,NB,NZ,NY,NX)
9810  CONTINUE
      ENDIF
C
C     CUT STALKS
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     HTSTKX=internode length removed
C     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
C          IHVST=3:reduction of clumping factor
C          IHVST=4 or 6:animal or insect biomass(g LM m-2)
C     FHGT=fraction of canopy layer height not harvested
C     FHVST=fraction of canopy layer mass not harvested
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C        leaf,non-foliar,woody, standing dead removed from PFT
C     THIN=IHVST=0-3,5: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C     WTSTK=stalk C mass
C
C
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(HTSTKX.GT.ZERO)THEN
      IF(IHVST(NZ,I,NY,NX).NE.3)THEN
      FHGT=AMAX1(0.0,AMIN1(1.0,HVST(NZ,I,NY,NX)/HTSTKX))
      ELSE
      FHGT=0.0
      ENDIF
      IF(THIN(NZ,I,NY,NX).EQ.0.0)THEN
      FHVST=AMAX1(0.0,1.0-(1.0-FHGT)*EHVST(1,3,NZ,I,NY,NX))
      FHVSH=FHVST
      ELSE
      FHVST=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      FHVSH=1.0-(1.0-FHGT)*EHVST(1,3,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
      ELSE
      FHVSH=FHVST
      ENDIF
      ENDIF
      ELSE
      FHVST=1.0
      FHVSH=1.0
      ENDIF
      ELSE
      IF(WTSTK(NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
      FHVST=AMAX1(0.0,AMIN1(1.0,1.0-WHVSTH/WTSTK(NZ,NY,NX)))
      FHVSH=FHVST
      ELSE
      FHVST=1.0
      FHVSH=1.0
      ENDIF
      ENDIF
C
C     HARVESTED STALK C,N,P
C
C     WTHTH3,WTHNH3,WTHPH3=harvested stalk C,N,P
C     WTHTX3,WTHNX3,WTHPX3=harvested stalk C,N,P to litter
C     WTSTKB,WTSTBN,WTSTBP=C,N,P mass remaining in harvested stalk
C
      WTHTH3=WTHTH3+(1.0-FHVSH)*WTSTKB(NB,NZ,NY,NX)
      WTHNH3=WTHNH3+(1.0-FHVSH)*WTSTBN(NB,NZ,NY,NX)
      WTHPH3=WTHPH3+(1.0-FHVSH)*WTSTBP(NB,NZ,NY,NX)
      WTHTX3=WTHTX3+(FHVSH-FHVST)*WTSTKB(NB,NZ,NY,NX)
      WTHNX3=WTHNX3+(FHVSH-FHVST)*WTSTBN(NB,NZ,NY,NX)
      WTHPX3=WTHPX3+(FHVSH-FHVST)*WTSTBP(NB,NZ,NY,NX)
C
C     REMAINING STALK C,N,P
C
C     WTSTXB,WTSTXN,WTSTXP=residual C,N,P mass in harvested stalk
C
      WTSTKB(NB,NZ,NY,NX)=FHVST*WTSTKB(NB,NZ,NY,NX)
      WTSTBN(NB,NZ,NY,NX)=FHVST*WTSTBN(NB,NZ,NY,NX)
      WTSTBP(NB,NZ,NY,NX)=FHVST*WTSTBP(NB,NZ,NY,NX)
      WVSTKB(NB,NZ,NY,NX)=FHVST*WVSTKB(NB,NZ,NY,NX)
      WTSTXB(NB,NZ,NY,NX)=FHVST*WTSTXB(NB,NZ,NY,NX)
      WTSTXN(NB,NZ,NY,NX)=FHVST*WTSTXN(NB,NZ,NY,NX)
      WTSTXP(NB,NZ,NY,NX)=FHVST*WTSTXP(NB,NZ,NY,NX)
C
C     CUT STALK NODES
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     HTNODX,HTNODE=stalk height,stalk internode length
C     FHGTK=fraction of internode length not harvested
C     THIN=IHVST=0-3,5: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C        leaf,non-foliar,woody, standing dead removed from PFT
C     WTSTK=stalk C mass
C     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
C
      DO 9820 K=25,0,-1
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(HTNODX(K,NB,NZ,NY,NX).GT.ZERO)THEN
      IF(IHVST(NZ,I,NY,NX).NE.3)THEN
      FHGTK=AMAX1(0.0,AMIN1(1.0,(HTNODE(K,NB,NZ,NY,NX)
     2-HVST(NZ,I,NY,NX))/HTNODX(K,NB,NZ,NY,NX)))
      ELSE
      FHGTK=0.0
      ENDIF
      IF(THIN(NZ,I,NY,NX).EQ.0.0)THEN
      FHVSTS=AMAX1(0.0,1.0-FHGTK*EHVST(1,3,NZ,I,NY,NX))
      ELSE
      FHVSTS=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
      ENDIF
      ELSE
      FHVSTS=1.0
      ENDIF
      ELSE
      IF(WTSTK(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVSTS=AMAX1(0.0,AMIN1(1.0,1.0-WHVSTH/WTSTK(NZ,NY,NX)))
      ELSE
      FHVSTS=1.0
      ENDIF
      ENDIF
      WGNODE(K,NB,NZ,NY,NX)=FHVSTS*WGNODE(K,NB,NZ,NY,NX)
      WGNODN(K,NB,NZ,NY,NX)=FHVSTS*WGNODN(K,NB,NZ,NY,NX)
      WGNODP(K,NB,NZ,NY,NX)=FHVSTS*WGNODP(K,NB,NZ,NY,NX)
      IF(IHVST(NZ,I,NY,NX).LE.2.AND.THIN(NZ,I,NY,NX).EQ.0.0)THEN
      HTNODX(K,NB,NZ,NY,NX)=FHVSTS*HTNODX(K,NB,NZ,NY,NX)
      HTNODE(K,NB,NZ,NY,NX)=AMIN1(HTNODE(K,NB,NZ,NY,NX)
     2,HVST(NZ,I,NY,NX))
      ENDIF
C     IF(NZ.EQ.2)THEN
C     WRITE(*,4811)'STK2',I,J,NX,NY,NZ,NB,K,IHVST(NZ,I,NY,NX)
C    2,HTNODX(K,NB,NZ,NY,NX),HTNODE(K,NB,NZ,NY,NX)
C    3,HVST(NZ,I,NY,NX),FHGTK,FHVSTS,ARLF(K,NB,NZ,NY,NX)
C    4,EHVST(1,3,NZ,I,NY,NX),THIN(NZ,I,NY,NX)
4811  FORMAT(A8,8I4,12E12.4)
C     ENDIF
9820  CONTINUE
C
C     CUT STALK RESERVES
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     WTSTKB=C mass remaining in harvested stalk
C     WTRSV=stalk reserve C mass
C     WHVRVH=remaining stalk reserve C mass
C     FHVST=fraction of reserve mass not harvested
C
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(WTSTKB(NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVST=FHVST
      FHVSH=FHVSH
      ELSE
      FHVST=0.0
      FHVSH=0.0
      ENDIF
      ELSE
      IF(WTRSV(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVST=AMAX1(0.0,AMIN1(1.0,1.0-WHVRVH/WTRSV(NZ,NY,NX)))
      FHVSH=FHVST
      ELSE
      FHVST=0.0
      FHVSH=0.0
      ENDIF
      ENDIF
C
C     HARVESTED STALK RESERVE C,N,P
C
C     WTHTH3,WTHNH3,WTHPH3=harvested stalk C,N,P
C     WTHTX3,WTHNX3,WTHPX3=harvested stalk C,N,P to litter
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C
      WTHTH3=WTHTH3+(1.0-FHVSH)*WTRSVB(NB,NZ,NY,NX)
      WTHNH3=WTHNH3+(1.0-FHVSH)*WTRSBN(NB,NZ,NY,NX)
      WTHPH3=WTHPH3+(1.0-FHVSH)*WTRSBP(NB,NZ,NY,NX)
      WTHTX3=WTHTX3+(FHVSH-FHVST)*WTRSVB(NB,NZ,NY,NX)
      WTHNX3=WTHNX3+(FHVSH-FHVST)*WTRSBN(NB,NZ,NY,NX)
      WTHPX3=WTHPX3+(FHVSH-FHVST)*WTRSBP(NB,NZ,NY,NX)
C
C     REMAINING STALK RESERVE C,N,P IF STALK REMAINING
C
      WTRSVB(NB,NZ,NY,NX)=FHVST*WTRSVB(NB,NZ,NY,NX)
      WTRSBN(NB,NZ,NY,NX)=FHVST*WTRSBN(NB,NZ,NY,NX)
      WTRSBP(NB,NZ,NY,NX)=FHVST*WTRSBP(NB,NZ,NY,NX)
C
C     CUT REPRODUCTIVE ORGANS
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
C          IHVST=3:reduction of clumping factor
C          IHVST=4 or 6:animal or insect biomass(g LM m-2)
C     THIN=IHVST=0-3,5: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C     FHVSTG,FHVSTH,FHVSTE=fraction of grain,husk,ear mass
C        not harvested
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C        leaf,non-foliar,woody, standing dead removed from PFT
C     WTHSK,WTEAR,WTGR=PFT husk,ear,grain C mass
C
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(HVST(NZ,I,NY,NX).LT.HTSTKX
     2.OR.IHVST(NZ,I,NY,NX).EQ.1
     3.OR.IHVST(NZ,I,NY,NX).EQ.3)THEN
      IF(THIN(NZ,I,NY,NX).EQ.0.0)THEN
      FHVSTG=1.0-EHVST(1,2,NZ,I,NY,NX)
      FHVSHG=FHVSTG
      ELSE
      FHVSTG=1.0-THIN(NZ,I,NY,NX)
      FHVSHG=1.0-EHVST(1,2,NZ,I,NY,NX)*THIN(NZ,I,NY,NX)
      ENDIF
      ELSE
      FHVSTG=1.0-THIN(NZ,I,NY,NX)
      FHVSHG=FHVSTG
      ENDIF
      FHVSTH=FHVSTG
      FHVSTE=FHVSTG
      FHVSHH=FHVSHG
      FHVSHE=FHVSHG
      ELSE
      IF(WTHSK(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVSTH=AMAX1(0.0,AMIN1(1.0,1.0-WHVHSH/WTHSK(NZ,NY,NX)))
      FHVSHH=FHVSTH
      ELSE
      FHVSTH=1.0
      FHVSHH=1.0
      ENDIF
      IF(WTEAR(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVSTE=AMAX1(0.0,AMIN1(1.0,1.0-WHVEAH/WTEAR(NZ,NY,NX)))
      FHVSHE=FHVSTE
      ELSE
      FHVSTE=1.0
      FHVSHE=1.0
      ENDIF
      IF(WTGR(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      FHVSTG=AMAX1(0.0,AMIN1(1.0,1.0-WHVGRH/WTGR(NZ,NY,NX)))
      FHVSHG=FHVSTG
      ELSE
      FHVSTG=1.0
      FHVSHG=1.0
      ENDIF
      ENDIF
C
C     HARVESTED REPRODUCTIVE C,N,P
C
C     WTHTH2,WTHNH2,WTHPH2=reproductive C,N,P removed
C     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
C     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
C     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
C     WTHTG,WTHNG,WTHPG=grain harvested
C
      WTHTH2=WTHTH2+(1.0-FHVSHH)*WTHSKB(NB,NZ,NY,NX)+(1.0-FHVSHE)
     2*WTEARB(NB,NZ,NY,NX)+(1.0-FHVSHG)*WTGRB(NB,NZ,NY,NX)
      WTHNH2=WTHNH2+(1.0-FHVSHH)*WTHSBN(NB,NZ,NY,NX)+(1.0-FHVSHE)
     2*WTEABN(NB,NZ,NY,NX)+(1.0-FHVSHG)*WTGRBN(NB,NZ,NY,NX)
      WTHPH2=WTHPH2+(1.0-FHVSHH)*WTHSBP(NB,NZ,NY,NX)+(1.0-FHVSHE)
     2*WTEABP(NB,NZ,NY,NX)+(1.0-FHVSHG)*WTGRBP(NB,NZ,NY,NX)
      WTHTX2=WTHTX2+(FHVSHH-FHVSTH)*WTHSKB(NB,NZ,NY,NX)
     2+(FHVSHE-FHVSTE)
     2*WTEARB(NB,NZ,NY,NX)+(FHVSHG-FHVSTG)*WTGRB(NB,NZ,NY,NX)
      WTHNX2=WTHNX2+(FHVSHH-FHVSTH)*WTHSBN(NB,NZ,NY,NX)
     2+(FHVSHE-FHVSTE)
     2*WTEABN(NB,NZ,NY,NX)+(FHVSHG-FHVSTG)*WTGRBN(NB,NZ,NY,NX)
      WTHPX2=WTHPX2+(FHVSHH-FHVSTH)*WTHSBP(NB,NZ,NY,NX)
     2+(FHVSHE-FHVSTE)
     2*WTEABP(NB,NZ,NY,NX)+(FHVSHG-FHVSTG)*WTGRBP(NB,NZ,NY,NX)
      WTHTG=WTHTG+(1.0-FHVSTG)*WTGRB(NB,NZ,NY,NX)
      WTHNG=WTHNG+(1.0-FHVSTG)*WTGRBN(NB,NZ,NY,NX)
      WTHPG=WTHPG+(1.0-FHVSTG)*WTGRBP(NB,NZ,NY,NX)
C
C     REMAINING REPRODUCTIVE C,N,P
C
C     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
C     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
C     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
C
      WTHSKB(NB,NZ,NY,NX)=FHVSTH*WTHSKB(NB,NZ,NY,NX)
      WTEARB(NB,NZ,NY,NX)=FHVSTE*WTEARB(NB,NZ,NY,NX)
      WTGRB(NB,NZ,NY,NX)=FHVSTG*WTGRB(NB,NZ,NY,NX)
      WTHSBN(NB,NZ,NY,NX)=FHVSTH*WTHSBN(NB,NZ,NY,NX)
      WTEABN(NB,NZ,NY,NX)=FHVSTE*WTEABN(NB,NZ,NY,NX)
      WTGRBN(NB,NZ,NY,NX)=FHVSTG*WTGRBN(NB,NZ,NY,NX)
      WTHSBP(NB,NZ,NY,NX)=FHVSTH*WTHSBP(NB,NZ,NY,NX)
      WTEABP(NB,NZ,NY,NX)=FHVSTE*WTEABP(NB,NZ,NY,NX)
      WTGRBP(NB,NZ,NY,NX)=FHVSTG*WTGRBP(NB,NZ,NY,NX)
      GRNXB(NB,NZ,NY,NX)=FHVSTG*GRNXB(NB,NZ,NY,NX)
      GRNOB(NB,NZ,NY,NX)=FHVSTG*GRNOB(NB,NZ,NY,NX)
      GRWTB(NB,NZ,NY,NX)=FHVSTG*GRWTB(NB,NZ,NY,NX)
C
C     REMAINING TOTAL BRANCH C,N,P AND LEAF, STALK AREA
C
C     CPOOLK=total C4 nonstructural C in branch
C     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
C     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
C     WTLSB=leaf+petiole mass
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
C     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
C     WTEARB,WTEABN,WTEABP=ear C,N,P mass
C     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
C     WVSTKB=stalk sapwood mass
C     PSILT=canopy water potential
C
      CPOOLK(NB,NZ,NY,NX)=0.0
      DO 1325 K=1,25
      CPOOLK(NB,NZ,NY,NX)=CPOOLK(NB,NZ,NY,NX)
     2+CPOOL3(K,NB,NZ,NY,NX)+CPOOL4(K,NB,NZ,NY,NX)
     2+CO2B(K,NB,NZ,NY,NX)+HCOB(K,NB,NZ,NY,NX)
1325  CONTINUE
      WTLSB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX)
     2+WTSHEB(NB,NZ,NY,NX))
      WTSHTB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX)
     2+WTSHEB(NB,NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
     3+WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)+WTGRB(NB,NZ,NY,NX)
     4+CPOOL(NB,NZ,NY,NX)+CPOOLK(NB,NZ,NY,NX))
      WTSHTN(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBN(NB,NZ,NY,NX)
     2+WTSHBN(NB,NZ,NY,NX)+WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX)
     3+WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX)
     4+ZPOOL(NB,NZ,NY,NX))
      WTSHTP(NB,NZ,NY,NX)=AMAX1(0.0,WTLFBP(NB,NZ,NY,NX)
     2+WTSHBP(NB,NZ,NY,NX)+WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX)
     3+WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX)
     4+PPOOL(NB,NZ,NY,NX))
C
C     RESET PHENOLOGY, GROWTH STAGE IF STALKS ARE CUT
C
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
C          IHVST=3:reduction of clumping factor
C          IHVST=4 or 6:animal or insect biomass(g LM m-2)
C     ZC=canopy height
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     VRNF,VRNX=leafoff hours,hours required for leafoff
C     IDAY(1,=emergence date
C     GROUP=node number required for floral initiation
C     PSTGI=node number at floral initiation
C     PSTGF=node number at flowering
C     VSTGX=leaf number on date of floral initiation
C     TGSTGI=total change in vegetative node number normalized
C        for maturity group
C     TGSTGF=total change in reproductive node number normalized
C        for maturity group
C     FLG4=number of hours with no grain fill
C     IFLGA=flag for initializing leafout
C
      IF((IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)
     2.AND.(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)
     3.AND.ZC(NZ,NY,NX).GT.HVST(NZ,I,NY,NX))THEN
      IF((IWTYP(NZ,NY,NX).NE.0.AND.VRNF(NB,NZ,NY,NX)
     2.LE.FVRN(IWTYP(NZ,NY,NX))*VRNX(NB,NZ,NY,NX))
     3.OR.(IWTYP(NZ,NY,NX).EQ.0
     4.AND.IDAY(1,NB,NZ,NY,NX).NE.0))THEN
      GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
      PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
      PSTGF(NB,NZ,NY,NX)=0.0
      VSTGX(NB,NZ,NY,NX)=0.0
      TGSTGI(NB,NZ,NY,NX)=0.0
      TGSTGF(NB,NZ,NY,NX)=0.0
      FLG4(NB,NZ,NY,NX)=0.0
      IDAY(1,NB,NZ,NY,NX)=I
      DO 3005 M=2,10
      IDAY(M,NB,NZ,NY,NX)=0
3005  CONTINUE
      IFLGA(NB,NZ,NY,NX)=0
      IF(NB.EQ.NB1(NZ,NY,NX))THEN
      DO 3010 NBX=1,NBR(NZ,NY,NX)
      IF(NBX.NE.NB1(NZ,NY,NX))THEN
      GROUP(NBX,NZ,NY,NX)=GROUPI(NZ,NY,NX)
      PSTGI(NBX,NZ,NY,NX)=PSTG(NBX,NZ,NY,NX)
      PSTGF(NBX,NZ,NY,NX)=0.0
      VSTGX(NBX,NZ,NY,NX)=0.0
      TGSTGI(NBX,NZ,NY,NX)=0.0
      TGSTGF(NBX,NZ,NY,NX)=0.0
      FLG4(NBX,NZ,NY,NX)=0.0
      IDAY(1,NBX,NZ,NY,NX)=I
      DO 3015 M=2,10
      IDAY(M,NBX,NZ,NY,NX)=0
3015  CONTINUE
      IFLGA(NBX,NZ,NY,NX)=0
      ENDIF
3010  CONTINUE
      ENDIF
      ENDIF
      ENDIF
C
C     DEATH OF BRANCH IF KILLING HARVEST ENTERED IN 'READQ'
C
C     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
C     IDTHB=branch living flag: 0=alive,1=dead
C     PP=PFT population
C     WTLS=total PFT leaf+petiole C mass
C     WTSTK=total PFT stalk C mass
C     WVSTK=total PFT sapwood C mass
C     ARSTK=total PFT stalk surface area
C
9835  CONTINUE
      WTLS(NZ,NY,NX)=0.0
      WTSTK(NZ,NY,NX)=0.0
      WVSTK(NZ,NY,NX)=0.0
      ARSTP(NZ,NY,NX)=0.0
      DO 9840 NB=1,NBR(NZ,NY,NX)
      WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
      WTSTK(NZ,NY,NX)=WTSTK(NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)
      WVSTK(NZ,NY,NX)=WVSTK(NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
      DO 9830 L=1,JC
      ARSTP(NZ,NY,NX)=ARSTP(NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
9830  CONTINUE
9840  CONTINUE
C
C     ROOT LITTERFALL FROM HARVESTING
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining
C        after disturbance
C     THIN=IHVST=0-3,5: fraction of population removed,
C          IHVST=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
C     THETW=soil water concentration
C     CORGC=SOC concentration
C     ITILL=soil disturbance type 1-20:tillage,21=litter removal
C                                 22=fire,23-24=drainage
C     FHVST,FHVSN,FHVSP=fraction of root layer C,N,P not removed
C        by disturbance
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
C     VCOXF,VNOXF,VPOXF=C,N,P gaseous emissions from disturbance
C     CNET=PFT net CO2 fixation
C     TNBP=total net biome productivity
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,
C        1=non-woody
C
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      XHVST=1.0-THIN(NZ,I,NY,NX)
      XHVSN=XHVST
      XHVSP=XHVST
      DO 3985 N=1,MY(NZ,NY,NX)
      DO 3980 L=NU(NY,NX),NJ(NY,NX)
      DO 3385 M=1,4
      FHVST=(1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
      FHVSN=(1.0-XHVSN)*CFOPN(0,M,NZ,NY,NX)*ZPOOLR(N,L,NZ,NY,NX)
      FHVSP=(1.0-XHVSP)*CFOPP(0,M,NZ,NY,NX)*PPOOLR(N,L,NZ,NY,NX)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+FHVST
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+FHVSN
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+FHVSP
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-FHVST
      HCNET(NZ,NY,NX)=HCNET(NZ,NY,NX)-FHVST
      TNBP(NY,NX)=TNBP(NY,NX)-0.0
      DO 3385 NR=1,NRT(NZ,NY,NX)
      FHVST=(1.0-XHVST)*CFOPC(5,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX)
     3+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0,NZ)
      FHVSN=(1.0-XHVSN)*CFOPN(5,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX)
     3+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0,NZ)
      FHVSP=(1.0-XHVSP)*CFOPP(5,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX)
     3+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0,NZ)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+FHVST
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+FHVSN
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+FHVSP
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-FHVST
      HCNET(NZ,NY,NX)=HCNET(NZ,NY,NX)-FHVST
      TNBP(NY,NX)=TNBP(NY,NX)-0.0
      FHVST=(1.0-XHVST)*CFOPC(4,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX)
     3+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1,NZ)
      FHVSN=(1.0-XHVSN)*CFOPN(4,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX)
     3+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1,NZ)
      FHVSP=(1.0-XHVSP)*CFOPP(4,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX)
     3+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1,NZ)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+FHVST
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+FHVSN
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+FHVSP
      CNET(NZ,NY,NX)=CNET(NZ,NY,NX)-FHVST
      HCNET(NZ,NY,NX)=HCNET(NZ,NY,NX)-FHVST
      TNBP(NY,NX)=TNBP(NY,NX)-0.0
3385  CONTINUE
C
C     RELEASE ROOT GAS CONTENTS DURING HARVESTING
C
C     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
C     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
C     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous
C        CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
C
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-(1.0-XHVST)
     2*(CO2A(N,L,NZ,NY,NX)+CO2P(N,L,NZ,NY,NX))
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-(1.0-XHVST)
     2*(OXYA(N,L,NZ,NY,NX)+OXYP(N,L,NZ,NY,NX))
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-(1.0-XHVST)
     2*(CH4A(N,L,NZ,NY,NX)+CH4P(N,L,NZ,NY,NX))
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-(1.0-XHVST)
     2*(Z2OA(N,L,NZ,NY,NX)+Z2OP(N,L,NZ,NY,NX))
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-(1.0-XHVST)
     2*(ZH3A(N,L,NZ,NY,NX)+ZH3P(N,L,NZ,NY,NX))
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-(1.0-XHVST)
     2*(H2GA(N,L,NZ,NY,NX)+H2GP(N,L,NZ,NY,NX))
      CO2A(N,L,NZ,NY,NX)=XHVST*CO2A(N,L,NZ,NY,NX)
      OXYA(N,L,NZ,NY,NX)=XHVST*OXYA(N,L,NZ,NY,NX)
      CH4A(N,L,NZ,NY,NX)=XHVST*CH4A(N,L,NZ,NY,NX)
      Z2OA(N,L,NZ,NY,NX)=XHVST*Z2OA(N,L,NZ,NY,NX)
      ZH3A(N,L,NZ,NY,NX)=XHVST*ZH3A(N,L,NZ,NY,NX)
      H2GA(N,L,NZ,NY,NX)=XHVST*H2GA(N,L,NZ,NY,NX)
      CO2P(N,L,NZ,NY,NX)=XHVST*CO2P(N,L,NZ,NY,NX)
      OXYP(N,L,NZ,NY,NX)=XHVST*OXYP(N,L,NZ,NY,NX)
      CH4P(N,L,NZ,NY,NX)=XHVST*CH4P(N,L,NZ,NY,NX)
      Z2OP(N,L,NZ,NY,NX)=XHVST*Z2OP(N,L,NZ,NY,NX)
      ZH3P(N,L,NZ,NY,NX)=XHVST*ZH3P(N,L,NZ,NY,NX)
      H2GP(N,L,NZ,NY,NX)=XHVST*H2GP(N,L,NZ,NY,NX)
C
C     REDUCE ROOT STATE VARIABLES DURING HARVESTING
C
C     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining
C        after disturbance
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
C     RTLG1,RTLG2=primary,secondary root length
C     RTN2=number of secondary root axes
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     WTRTL,WTRTD=active,actual root C mass
C     WSRTL=root protein C mass
C     RTN1,RTNL=number of primary,secondary root axes
C     RTDNP,RTLGP=root length density,root length per plant
C     RTVLW,RTVLP=root or myco aqueous,gaseous volume
C     RTARP=root surface area per plant
C     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
C
      DO 3960 NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)*XHVST
      WTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)*XHVST
      WTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)*XHVSN
      WTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)*XHVSN
      WTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)*XHVSP
      WTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)*XHVSP
      RTWT1(N,NR,NZ,NY,NX)=RTWT1(N,NR,NZ,NY,NX)*XHVST
      RTWT1N(N,NR,NZ,NY,NX)=RTWT1N(N,NR,NZ,NY,NX)*XHVST
      RTWT1P(N,NR,NZ,NY,NX)=RTWT1P(N,NR,NZ,NY,NX)*XHVST
      RTLG1(N,L,NR,NZ,NY,NX)=RTLG1(N,L,NR,NZ,NY,NX)*XHVST
      RTLG2(N,L,NR,NZ,NY,NX)=RTLG2(N,L,NR,NZ,NY,NX)*XHVST
      RTN2(N,L,NR,NZ,NY,NX)=RTN2(N,L,NR,NZ,NY,NX)*XHVST
3960  CONTINUE
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)*XHVST
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)*XHVSN
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)*XHVSP
      WTRTL(N,L,NZ,NY,NX)=WTRTL(N,L,NZ,NY,NX)*XHVST
      WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)*XHVST
      WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX)*XHVST
      RTN1(N,L,NZ,NY,NX)=RTN1(N,L,NZ,NY,NX)*XHVST
      RTNL(N,L,NZ,NY,NX)=RTNL(N,L,NZ,NY,NX)*XHVST
      RTLGP(N,L,NZ,NY,NX)=RTLGP(N,L,NZ,NY,NX)*XHVST
      RTDNP(N,L,NZ,NY,NX)=RTDNP(N,L,NZ,NY,NX)*XHVST
      RTVLP(N,L,NZ,NY,NX)=RTVLP(N,L,NZ,NY,NX)*XHVST
      RTVLW(N,L,NZ,NY,NX)=RTVLW(N,L,NZ,NY,NX)*XHVST
      RTARP(N,L,NZ,NY,NX)=RTARP(N,L,NZ,NY,NX)*XHVST
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)*XHVST
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)*XHVST
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)*XHVST
C
C     NODULE LITTERFALL AND STATE VARIABLES DURING HARVESTING
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
C     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining
C        after disturbance
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C
      IF(INTYP(NZ,NY,NX).NE.0.AND.N.EQ.1)THEN
      DO 3395 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*(CFOPC(4,M,NZ,NY,NX)*WTNDL(L,NZ,NY,NX)
     3+CFOPC(0,M,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX))
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVSN)
     2*(CFOPN(4,M,NZ,NY,NX)*WTNDLN(L,NZ,NY,NX)
     3+CFOPN(0,M,NZ,NY,NX)*ZPOOLN(L,NZ,NY,NX))
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVSP)
     2*(CFOPP(4,M,NZ,NY,NX)*WTNDLP(L,NZ,NY,NX)
     3+CFOPP(0,M,NZ,NY,NX)*PPOOLN(L,NZ,NY,NX))
3395  CONTINUE
      WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)*XHVST
      WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)*XHVSN
      WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)*XHVSP
      CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)*XHVST
      ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)*XHVSN
      PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)*XHVSP
      ENDIF
3980  CONTINUE
3985  CONTINUE
C
C     STORAGE LITTERFALL AND STATE VARIABLES DURING HARVESTING
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining
C        after disturbance
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C
      IF(ISTYP(NZ,NY,NX).NE.0)THEN
      DO 3400 M=1,4
      CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(0,NZ)
      ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVSN)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(0,NZ)
      PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVSP)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(0,NZ)
      CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(1,NZ)
      ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVSN)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(1,NZ)
      PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVSP)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(1,NZ)
3400  CONTINUE
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)*XHVST
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)*XHVSN
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)*XHVSP
      ENDIF
      ENDIF
      ENDIF
C
C     REDUCE OR REMOVE PLANT POPULATIONS DURING TILLAGE
C
C     ZNOON=hour of solar noon
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     IDAY0,IYR0=day,year of planting
C     IYRC=current year
C     ITILL=soil disturbance type 1-20:tillage,21=litter removal
C                                 22=fire,23-24=drainage
C     XHVST=fraction of PFT remaining after disturbance
C     PPX,PP=PFT population per m2,grid cell
C     FRADP=fraction of radiation received by each PFT canopy
C
      IF(J.EQ.INT(ZNOON(NY,NX)).AND.(IBTYP(NZ,NY,NX).EQ.0
     2.OR.IGTYP(NZ,NY,NX).LE.1).AND.(I.NE.IDAY0(NZ,NY,NX)
     3.OR.IYRC.NE.IYR0(NZ,NY,NX)))THEN
      IF(ITILL(I,NY,NX).LE.10.OR.NZ.NE.1)THEN
      IF(I.GT.IDAY0(NZ,NY,NX).OR.IYRC.GT.IYR0(NZ,NY,NX))THEN
      XHVST=XCORP(NY,NX)
      PPX(NZ,NY,NX)=PPX(NZ,NY,NX)*XHVST
      PP(NZ,NY,NX)=PP(NZ,NY,NX)*XHVST
      FRADP(NZ,NY,NX)=FRADP(NZ,NY,NX)*XHVST
      WTLS(NZ,NY,NX)=0.0
      WVSTK(NZ,NY,NX)=0.0
C
C     TERMINATE BRANCHES IF TILLAGE IMPLEMENT 10 IS SELECTED
C
C     IDTHB=branch living flag: 0=alive,1=dead
C     PP=PFT population
C
      DO 8975 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      IF(PP(NZ,NY,NX).LE.0.0)IDTHB(NB,NZ,NY,NX)=1
C
C     LITTERFALL FROM BRANCHES DURING TILLAGE
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
C     XHVST=fraction of PFT remaining after disturbance
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     CPOOLK=total C4 nonstructural C in branch
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
C     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
C     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
C     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C
      DO 6380 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST)
     2*(CFOPC(0,M,NZ,NY,NX)*(CPOOL(NB,NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX)
     3+CPOOLK(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX))
     4+CFOPC(1,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(1,NZ)
     5+WTNDB(NB,NZ,NY,NX))
     6+CFOPC(2,M,NZ,NY,NX)*(WTSHEB(NB,NZ,NY,NX)*FWODB(1,NZ)
     7+WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)))
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPC(5,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(0,NZ)
     3+WTSHEB(NB,NZ,NY,NX)*FWODB(0,NZ))
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST)
     2*(CFOPN(0,M,NZ,NY,NX)*(ZPOOL(NB,NZ,NY,NX)+ZPOLNB(NB,NZ,NY,NX)
     3+WTRSBN(NB,NZ,NY,NX))
     4+CFOPN(1,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(1,NZ)
     5+WTNDBN(NB,NZ,NY,NX))
     6+CFOPN(2,M,NZ,NY,NX)*(WTSHBN(NB,NZ,NY,NX)*FWODSN(1,NZ)
     7+WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)))
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPN(5,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(0,NZ)
     3+WTSHBN(NB,NZ,NY,NX)*FWODSN(0,NZ))
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST)
     2*(CFOPP(0,M,NZ,NY,NX)*(PPOOL(NB,NZ,NY,NX)+PPOLNB(NB,NZ,NY,NX)
     3+WTRSBP(NB,NZ,NY,NX))
     4+CFOPP(1,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(1,NZ)
     5+WTNDBP(NB,NZ,NY,NX))
     6+CFOPP(2,M,NZ,NY,NX)*(WTSHBP(NB,NZ,NY,NX)*FWODSP(1,NZ)
     7+WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)))
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPP(5,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(0,NZ)
     3+WTSHBP(NB,NZ,NY,NX)*FWODSP(0,NZ))
C
C     IF WINTER ANNUAL
C
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+(1.0-XHVST)
     2*CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+(1.0-XHVST)
     2*CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+(1.0-XHVST)
     2*CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ELSE
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ENDIF
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPC(5,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)*FWOOD(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPN(5,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)*FWOODN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPP(5,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)*FWOODP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPC(3,M,NZ,NY,NX)*WTSTKB(NB,NZ,NY,NX)*FWOOD(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPN(3,M,NZ,NY,NX)*WTSTBN(NB,NZ,NY,NX)*FWOODN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPP(3,M,NZ,NY,NX)*WTSTBP(NB,NZ,NY,NX)*FWOODP(1,NZ)
6380  CONTINUE
C
C     PLANT STATE VARIABLES REMAINING AFTER TILLAGE
C
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)*XHVST
      CPOOLK(NB,NZ,NY,NX)=CPOOLK(NB,NZ,NY,NX)*XHVST
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)*XHVST
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)*XHVST
      CPOLNB(NB,NZ,NY,NX)=CPOLNB(NB,NZ,NY,NX)*XHVST
      ZPOLNB(NB,NZ,NY,NX)=ZPOLNB(NB,NZ,NY,NX)*XHVST
      PPOLNB(NB,NZ,NY,NX)=PPOLNB(NB,NZ,NY,NX)*XHVST
      WTSHTB(NB,NZ,NY,NX)=WTSHTB(NB,NZ,NY,NX)*XHVST
      WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)*XHVST
      WTNDB(NB,NZ,NY,NX)=WTNDB(NB,NZ,NY,NX)*XHVST
      WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX)*XHVST
      WTSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)*XHVST
      WVSTKB(NB,NZ,NY,NX)=WVSTKB(NB,NZ,NY,NX)*XHVST
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)*XHVST
      WTHSKB(NB,NZ,NY,NX)=WTHSKB(NB,NZ,NY,NX)*XHVST
      WTEARB(NB,NZ,NY,NX)=WTEARB(NB,NZ,NY,NX)*XHVST
      WTGRB(NB,NZ,NY,NX)=WTGRB(NB,NZ,NY,NX)*XHVST
      WTSHTN(NB,NZ,NY,NX)=WTSHTN(NB,NZ,NY,NX)*XHVST
      WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)*XHVST
      WTNDBN(NB,NZ,NY,NX)=WTNDBN(NB,NZ,NY,NX)*XHVST
      WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX)*XHVST
      WTSTBN(NB,NZ,NY,NX)=WTSTBN(NB,NZ,NY,NX)*XHVST
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)*XHVST
      WTHSBN(NB,NZ,NY,NX)=WTHSBN(NB,NZ,NY,NX)*XHVST
      WTEABN(NB,NZ,NY,NX)=WTEABN(NB,NZ,NY,NX)*XHVST
      WTGRBN(NB,NZ,NY,NX)=WTGRBN(NB,NZ,NY,NX)*XHVST
      WTSHTP(NB,NZ,NY,NX)=WTSHTP(NB,NZ,NY,NX)*XHVST
      WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)*XHVST
      WTNDBP(NB,NZ,NY,NX)=WTNDBP(NB,NZ,NY,NX)*XHVST
      WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX)*XHVST
      WTSTBP(NB,NZ,NY,NX)=WTSTBP(NB,NZ,NY,NX)*XHVST
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)*XHVST
      WTHSBP(NB,NZ,NY,NX)=WTHSBP(NB,NZ,NY,NX)*XHVST
      WTEABP(NB,NZ,NY,NX)=WTEABP(NB,NZ,NY,NX)*XHVST
      WTGRBP(NB,NZ,NY,NX)=WTGRBP(NB,NZ,NY,NX)*XHVST
      GRNXB(NB,NZ,NY,NX)=GRNXB(NB,NZ,NY,NX)*XHVST
      GRNOB(NB,NZ,NY,NX)=GRNOB(NB,NZ,NY,NX)*XHVST
      GRWTB(NB,NZ,NY,NX)=GRWTB(NB,NZ,NY,NX)*XHVST
      ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)*XHVST
      WTLSB(NB,NZ,NY,NX)=AMAX1(0.0,WTLFB(NB,NZ,NY,NX)
     2+WTSHEB(NB,NZ,NY,NX))
      WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
      WTSTXB(NB,NZ,NY,NX)=WTSTXB(NB,NZ,NY,NX)*XHVST
      WTSTXN(NB,NZ,NY,NX)=WTSTXN(NB,NZ,NY,NX)*XHVST
      WTSTXP(NB,NZ,NY,NX)=WTSTXP(NB,NZ,NY,NX)*XHVST
      WVSTK(NZ,NY,NX)=WVSTK(NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
      DO 8970 K=0,25
      IF(K.NE.0)THEN
      CPOOL3(K,NB,NZ,NY,NX)=CPOOL3(K,NB,NZ,NY,NX)*XHVST
      CPOOL4(K,NB,NZ,NY,NX)=CPOOL4(K,NB,NZ,NY,NX)*XHVST
      CO2B(K,NB,NZ,NY,NX)=CO2B(K,NB,NZ,NY,NX)*XHVST
      HCOB(K,NB,NZ,NY,NX)=HCOB(K,NB,NZ,NY,NX)*XHVST
      ENDIF
      ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX)*XHVST
      WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX)*XHVST
      WSLF(K,NB,NZ,NY,NX)=WSLF(K,NB,NZ,NY,NX)*XHVST
C     HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX)*XHVST
      WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX)*XHVST
      WSSHE(K,NB,NZ,NY,NX)=WSSHE(K,NB,NZ,NY,NX)*XHVST
C     HTNODE(K,NB,NZ,NY,NX)=HTNODE(K,NB,NZ,NY,NX)*XHVST
C     HTNODX(K,NB,NZ,NY,NX)=HTNODX(K,NB,NZ,NY,NX)*XHVST
      WGNODE(K,NB,NZ,NY,NX)=WGNODE(K,NB,NZ,NY,NX)*XHVST
      WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)*XHVST
      WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX)*XHVST
      WGNODN(K,NB,NZ,NY,NX)=WGNODN(K,NB,NZ,NY,NX)*XHVST
      WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)*XHVST
      WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX)*XHVST
      WGNODP(K,NB,NZ,NY,NX)=WGNODP(K,NB,NZ,NY,NX)*XHVST
      DO 8965 L=1,JC
      ARLFL(L,K,NB,NZ,NY,NX)=ARLFL(L,K,NB,NZ,NY,NX)*XHVST
      WGLFL(L,K,NB,NZ,NY,NX)=WGLFL(L,K,NB,NZ,NY,NX)*XHVST
      WGLFLN(L,K,NB,NZ,NY,NX)=WGLFLN(L,K,NB,NZ,NY,NX)*XHVST
      WGLFLP(L,K,NB,NZ,NY,NX)=WGLFLP(L,K,NB,NZ,NY,NX)*XHVST
8965  CONTINUE
8970  CONTINUE
      ENDIF
8975  CONTINUE
C
C     TERMINATE ROOTS IF TILLAGE IMPLEMENT 10 IS SELECTED
C
C     PP=PFT population
C     IDTHR,IDTHP=PFT root,shoot living flag: 0=alive,1=dead
C     IDTH=PFT living flag: 0=alive,1=dead
C     JHVST=terminate PFT:0=no,1=yes,2=yes,and reseed
C     IDAYH,IYRH=day,year of harvesting
C     IYRC=current year
C
      IF(PP(NZ,NY,NX).LE.ZERO)THEN
      IDTHR(NZ,NY,NX)=1
      IDTHP(NZ,NY,NX)=1
      IDTH(NZ,NY,NX)=1
      JHVST(NZ,I,NY,NX)=1
      IDAYH(NZ,NY,NX)=I
      IYRH(NZ,NY,NX)=IYRC
      ENDIF
C
C     LITTERFALL FROM ROOTS DURING TILLAGE
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
C     XHVST=fraction of PFT remaining after disturbance
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C
      DO 8985 N=1,MY(NZ,NY,NX)
      DO 8980 L=NU(NY,NX),NJ(NY,NX)
      DO 6385 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPC(0,M,NZ,NY,NX)*CPOOLR(N,L,NZ,NY,NX)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPN(0,M,NZ,NY,NX)*ZPOOLR(N,L,NZ,NY,NX)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPP(0,M,NZ,NY,NX)*PPOOLR(N,L,NZ,NY,NX)
      DO 6385 NR=1,NRT(NZ,NY,NX)
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPC(5,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX)
     3+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0,NZ)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPN(5,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX)
     3+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0,NZ)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPP(5,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX)
     3+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0,NZ)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPC(4,M,NZ,NY,NX)*(WTRT1(N,L,NR,NZ,NY,NX)
     3+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1,NZ)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPN(4,M,NZ,NY,NX)*(WTRT1N(N,L,NR,NZ,NY,NX)
     3+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1,NZ)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*CFOPP(4,M,NZ,NY,NX)*(WTRT1P(N,L,NR,NZ,NY,NX)
     3+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1,NZ)
6385  CONTINUE
C
C     RELEASE ROOT GAS CONTENTS DURING TILLAGE
C
C     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
C     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
C     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous
C        CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
C
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-(1.0-XHVST)
     2*(CO2A(N,L,NZ,NY,NX)+CO2P(N,L,NZ,NY,NX))
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-(1.0-XHVST)
     2*(OXYA(N,L,NZ,NY,NX)+OXYP(N,L,NZ,NY,NX))
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-(1.0-XHVST)
     2*(CH4A(N,L,NZ,NY,NX)+CH4P(N,L,NZ,NY,NX))
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-(1.0-XHVST)
     2*(Z2OA(N,L,NZ,NY,NX)+Z2OP(N,L,NZ,NY,NX))
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-(1.0-XHVST)
     2*(ZH3A(N,L,NZ,NY,NX)+ZH3P(N,L,NZ,NY,NX))
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-(1.0-XHVST)
     2*(H2GA(N,L,NZ,NY,NX)+H2GP(N,L,NZ,NY,NX))
      CO2A(N,L,NZ,NY,NX)=XHVST*CO2A(N,L,NZ,NY,NX)
      OXYA(N,L,NZ,NY,NX)=XHVST*OXYA(N,L,NZ,NY,NX)
      CH4A(N,L,NZ,NY,NX)=XHVST*CH4A(N,L,NZ,NY,NX)
      Z2OA(N,L,NZ,NY,NX)=XHVST*Z2OA(N,L,NZ,NY,NX)
      ZH3A(N,L,NZ,NY,NX)=XHVST*ZH3A(N,L,NZ,NY,NX)
      H2GA(N,L,NZ,NY,NX)=XHVST*H2GA(N,L,NZ,NY,NX)
      CO2P(N,L,NZ,NY,NX)=XHVST*CO2P(N,L,NZ,NY,NX)
      OXYP(N,L,NZ,NY,NX)=XHVST*OXYP(N,L,NZ,NY,NX)
      CH4P(N,L,NZ,NY,NX)=XHVST*CH4P(N,L,NZ,NY,NX)
      Z2OP(N,L,NZ,NY,NX)=XHVST*Z2OP(N,L,NZ,NY,NX)
      ZH3P(N,L,NZ,NY,NX)=XHVST*ZH3P(N,L,NZ,NY,NX)
      H2GP(N,L,NZ,NY,NX)=XHVST*H2GP(N,L,NZ,NY,NX)
C
C     ROOT STATE VARIABLES REMAINING AFTER TILLAGE
C
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
C     RTLG1,RTLG2=primary,secondary root length
C     RTN2=number of secondary root axes
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     WTRTL,WTRTD=active,actual root C mass
C     WSRTL=root protein C mass
C     RTN1,RTNL=number of primary,secondary root axes
C     RTDNP,RTLGP=root length density,root length per plant
C     RTVLW,RTVLP=root or myco aqueous,gaseous volume
C     RTARP=root surface area per plant
C     RCO2M,RCO2N,RCO2A unlimited by O2,nonstructural C
C
      DO 8960 NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)*XHVST
      WTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)*XHVST
      WTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)*XHVST
      WTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)*XHVST
      WTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)*XHVST
      WTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)*XHVST
      RTWT1(N,NR,NZ,NY,NX)=RTWT1(N,NR,NZ,NY,NX)*XHVST
      RTWT1N(N,NR,NZ,NY,NX)=RTWT1N(N,NR,NZ,NY,NX)*XHVST
      RTWT1P(N,NR,NZ,NY,NX)=RTWT1P(N,NR,NZ,NY,NX)*XHVST
      RTLG1(N,L,NR,NZ,NY,NX)=RTLG1(N,L,NR,NZ,NY,NX)*XHVST
      RTLG2(N,L,NR,NZ,NY,NX)=RTLG2(N,L,NR,NZ,NY,NX)*XHVST
      RTN2(N,L,NR,NZ,NY,NX)=RTN2(N,L,NR,NZ,NY,NX)*XHVST
8960  CONTINUE
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)*XHVST
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)*XHVST
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)*XHVST
      WTRTL(N,L,NZ,NY,NX)=WTRTL(N,L,NZ,NY,NX)*XHVST
      WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)*XHVST
      WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX)*XHVST
      RTN1(N,L,NZ,NY,NX)=RTN1(N,L,NZ,NY,NX)*XHVST
      RTNL(N,L,NZ,NY,NX)=RTNL(N,L,NZ,NY,NX)*XHVST
      RTLGP(N,L,NZ,NY,NX)=RTLGP(N,L,NZ,NY,NX)*XHVST
      RTDNP(N,L,NZ,NY,NX)=RTDNP(N,L,NZ,NY,NX)*XHVST
      RTVLP(N,L,NZ,NY,NX)=RTVLP(N,L,NZ,NY,NX)*XHVST
      RTVLW(N,L,NZ,NY,NX)=RTVLW(N,L,NZ,NY,NX)*XHVST
      RTARP(N,L,NZ,NY,NX)=RTARP(N,L,NZ,NY,NX)*XHVST
      RCO2M(N,L,NZ,NY,NX)=RCO2M(N,L,NZ,NY,NX)*XHVST
      RCO2N(N,L,NZ,NY,NX)=RCO2N(N,L,NZ,NY,NX)*XHVST
      RCO2A(N,L,NZ,NY,NX)=RCO2A(N,L,NZ,NY,NX)*XHVST
C
C     LITTERFALL AND STATE VARIABLES FOR NODULES DURING TILLAGE
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
C     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining
C        after disturbance
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C
      IF(INTYP(NZ,NY,NX).NE.0.AND.N.EQ.1)THEN
      DO 6395 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*(CFOPC(4,M,NZ,NY,NX)*WTNDL(L,NZ,NY,NX)
     3+CFOPC(0,M,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX))
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*(CFOPN(4,M,NZ,NY,NX)*WTNDLN(L,NZ,NY,NX)
     3+CFOPN(0,M,NZ,NY,NX)*ZPOOLN(L,NZ,NY,NX))
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+(1.0-XHVST)
     2*(CFOPP(4,M,NZ,NY,NX)*WTNDLP(L,NZ,NY,NX)
     3+CFOPP(0,M,NZ,NY,NX)*PPOOLN(L,NZ,NY,NX))
6395  CONTINUE
      WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)*XHVST
      WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)*XHVST
      WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)*XHVST
      CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)*XHVST
      ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)*XHVST
      PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)*XHVST
      ENDIF
8980  CONTINUE
8985  CONTINUE
C
C     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE RESERVES
C     DURING TILLAGE
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     XHVST,XHVSN,XHVSP=fraction of root C,N,P remaining
C        after disturbance
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C
      DO 6400 M=1,4
      CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(0,NZ)
      ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVST)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(0,NZ)
      PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVST)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(0,NZ)
      CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVST)*CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX))*FWOOD(1,NZ)
      ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVST)*CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX))*FWOODN(1,NZ)
      PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)
     2+((1.0-XHVST)*CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX))*FWOODP(1,NZ)
6400  CONTINUE
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)*XHVST
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)*XHVST
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)*XHVST
      ENDIF
      ENDIF
      ENDIF
C
C     HARVEST STANDING DEAD
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     THIN=thinning:fraction of population removed
C     FHVST=fraction of standing dead mass not harvested
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C           leaf,non-foliar,woody, standing dead removed from PFT
C     HVST=IHVST=0-2:>0=cutting height,<0=fraction of LAI removed
C          IHVST=3:reduction of clumping factor
C          IHVST=4 or 6:animal or insect biomass(g LM m-2)
C     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
C     WTHTH4,WTHNH4,WTHPH4=harvested standing dead C,N,P
C     WTHTX4,WTHNX4,WTHPX4=harvested standing dead C,N,P to litter
C
      IF(IHVST(NZ,I,NY,NX).GE.0)THEN
      IF(J.EQ.INT(ZNOON(NY,NX)).AND.IHVST(NZ,I,NY,NX).NE.4
     2.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
      IF(THIN(NZ,I,NY,NX).EQ.0.0)THEN
      FHVST=AMAX1(0.0,1.0-EHVST(1,4,NZ,I,NY,NX))
      FHVSH=FHVST
      ELSE
      FHVST=AMAX1(0.0,1.0-THIN(NZ,I,NY,NX))
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      FHVSH=AMAX1(0.0,1.0-EHVST(1,4,NZ,I,NY,NX)*THIN(NZ,I,NY,NX))
      ELSE
      FHVSH=FHVST
      ENDIF
      ENDIF
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      IF(WTSTG(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WHVSTD=HVST(NZ,I,NY,NX)*THIN(NZ,I,NY,NX)*0.45/24.0
     2*AREA(3,NU(NY,NX),NY,NX)*EHVST(1,4,NZ,I,NY,NX)
      FHVST=AMAX1(0.0,1.0-WHVSTD/WTSTG(NZ,NY,NX))
      FHVSH=FHVST
      ELSE
      FHVST=1.0
      FHVSH=1.0
      ENDIF
      ELSE
      FHVST=1.0
      FHVSH=1.0
      ENDIF
      DO 6475 M=1,5
      WTHTH4=WTHTH4+(1.0-FHVSH)*WTSTDG(M,NZ,NY,NX)
      WTHNH4=WTHNH4+(1.0-FHVSH)*WTSTDN(M,NZ,NY,NX)
      WTHPH4=WTHPH4+(1.0-FHVSH)*WTSTDP(M,NZ,NY,NX)
      WTHTX4=WTHTX4+(FHVSH-FHVST)*WTSTDG(M,NZ,NY,NX)
      WTHNX4=WTHNX4+(FHVSH-FHVST)*WTSTDN(M,NZ,NY,NX)
      WTHPX4=WTHPX4+(FHVSH-FHVST)*WTSTDP(M,NZ,NY,NX)
      WTSTDG(M,NZ,NY,NX)=FHVST*WTSTDG(M,NZ,NY,NX)
      WTSTDN(M,NZ,NY,NX)=FHVST*WTSTDN(M,NZ,NY,NX)
      WTSTDP(M,NZ,NY,NX)=FHVST*WTSTDP(M,NZ,NY,NX)
6475  CONTINUE
C
C     IF NO PLANT C,N,P REMOVED AT HARVEST (ALL RESIDUE RETURNED)
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     WTHTH0,WTHNH0,WTHPH0=nonstructural C,N,P removed
C     WTHTH1,WTHNH1,WTHPH1=leaf C,N,P removed
C     WTHTH2,WTHNH2,WTHPH2=fine,non-leaf C,N,P removed
C     WTHTH3,WTHNH3,WTHPH3=woody C,N,P removed
C     WTHTH4,WTHNH4,WTHPH4=standing dead C,N,P removed
C     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
C     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
C     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
C     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
C     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
C     EHVST(1,1,EHVST(1,2,EHVST(1,3,EHVST(1,4=fraction of
C        leaf,non-foliar,woody, standing dead removed from PFT
C
      IF(IHVST(NZ,I,NY,NX).EQ.0)THEN
      WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
C
C     IF ONLY GRAIN C,N,P REMOVED AT HARVEST
C
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.1)THEN
      WTHTR0=WTHTH0
      WTHNR0=WTHNH0
      WTHPR0=WTHPH0
      WTHTR1=WTHTH1
      WTHNR1=WTHNH1
      WTHPR1=WTHPH1
      WTHTR2=WTHTH2-WTHTG*EHVST(2,2,NZ,I,NY,NX)
      WTHNR2=WTHNH2-WTHNG*EHVST(2,2,NZ,I,NY,NX)
      WTHPR2=WTHPH2-WTHPG*EHVST(2,2,NZ,I,NY,NX)
      WTHTR3=WTHTH3
      WTHNR3=WTHNH3
      WTHPR3=WTHPH3
      WTHTR4=WTHTH4
      WTHNR4=WTHNH4
      WTHPR4=WTHPH4
C
C     IF ONLY WOOD C,N,P REMOVED AT HARVEST
C
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.2)THEN
      WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
C
C     IF ALL PLANT C,N,P REMOVED AT HARVEST (NO RESIDUE RETURNED)
C
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.3)THEN
      WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
C
C     IF PLANT C,N,P REMOVED BY GRAZING
C
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.4.OR.IHVST(NZ,I,NY,NX).EQ.6)THEN
      WTHTR0=WTHTH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR0=WTHNH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR0=WTHPH0*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR1=WTHTH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHNR1=WTHNH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHPR1=WTHPH1*(1.0-EHVST(2,1,NZ,I,NY,NX))
      WTHTR2=WTHTH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHNR2=WTHNH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHPR2=WTHPH2*(1.0-EHVST(2,2,NZ,I,NY,NX))
      WTHTR3=WTHTH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHNR3=WTHNH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHPR3=WTHPH3*(1.0-EHVST(2,3,NZ,I,NY,NX))
      WTHTR4=WTHTH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHNR4=WTHNH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      WTHPR4=WTHPH4*(1.0-EHVST(2,4,NZ,I,NY,NX))
      ENDIF
C
C     TOTAL C,N,P REMOVAL FROM DISTURBANCE
C
C     WTHTHT,WTHNHT,WTHPHT=total C,N,P removed
C     WTHTRT,WTHNRT,WTHPRT=total removed C,N,P to litter
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
C     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
C     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem
C        from all PFT
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C
      WTHTHT=WTHTH0+WTHTH1+WTHTH2+WTHTH3+WTHTH4
      WTHTRT=WTHTR0+WTHTR1+WTHTR2+WTHTR3+WTHTR4
      WTHNHT=WTHNH0+WTHNH1+WTHNH2+WTHNH3+WTHNH4
      WTHNRT=WTHNR0+WTHNR1+WTHNR2+WTHNR3+WTHNR4
      WTHPHT=WTHPH0+WTHPH1+WTHPH2+WTHPH3+WTHPH4
      WTHPRT=WTHPR0+WTHPR1+WTHPR2+WTHPR3+WTHPR4
      WTHTXT=WTHTX0+WTHTX1+WTHTX2+WTHTX3+WTHTX4
      WTHNXT=WTHNX0+WTHNX1+WTHNX2+WTHNX3+WTHNX4
      WTHPXT=WTHPX0+WTHPX1+WTHPX2+WTHPX3+WTHPX4
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
C     IF(IHVST(NZ,I,NY,NX).NE.5)THEN
      IF(JHVST(NZ,I,NY,NX).NE.2)THEN
      HVSTC(NZ,NY,NX)=HVSTC(NZ,NY,NX)+WTHTHT-WTHTRT
      HVSTN(NZ,NY,NX)=HVSTN(NZ,NY,NX)+WTHNHT-WTHNRT
      HVSTP(NZ,NY,NX)=HVSTP(NZ,NY,NX)+WTHPHT-WTHPRT
      TNBP(NY,NX)=TNBP(NY,NX)+WTHTRT-WTHTHT
      XHVSTC(NY,NX)=XHVSTC(NY,NX)+WTHTHT-WTHTRT
      XHVSTN(NY,NX)=XHVSTN(NY,NX)+WTHNHT-WTHNRT
      XHVSTP(NY,NX)=XHVSTP(NY,NX)+WTHPHT-WTHPRT
      ELSE
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+WTHTHT-WTHTRT
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+WTHNHT-WTHNRT
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+WTHPHT-WTHPRT
      ENDIF
C     ENDIF
C
C     C,N,P REMOVED FROM GRAZING
C
C     HVSTC,HVSTN,HVSTP=total C,N,P removed from ecosystem from PFT
C     XHVSTC,XHVSTN,XHVSTP=total C,N,P removed from ecosystem
C        from all PFT
C     GY=growth yield of grazers
C     WTHTHT,WTHNHT,WTHPHT=total C,N,P removed
C     WTHTRT,WTHNRT,WTHPRT=total removed C,N,P to litter
C     RECO=ecosystem respiration
C     TRAU=total autotrophic respiration
C
      ELSE
      HVSTC(NZ,NY,NX)=HVSTC(NZ,NY,NX)+GY*(WTHTHT-WTHTRT)
      HVSTN(NZ,NY,NX)=HVSTN(NZ,NY,NX)+WTHNHT-WTHNRT
      HVSTP(NZ,NY,NX)=HVSTP(NZ,NY,NX)+WTHPHT-WTHPRT
      TCO2T(NZ,NY,NX)=TCO2T(NZ,NY,NX)-GZ*(WTHTHT-WTHTRT)
      TCO2A(NZ,NY,NX)=TCO2A(NZ,NY,NX)-GZ*(WTHTHT-WTHTRT)
C     TNBP(NY,NX)=TNBP(NY,NX)+GY*(WTHTRT-WTHTHT)
C     CNET(NZ,NY,NX)=CNET(NZ,NY,NX)+GZ*(WTHTRT-WTHTHT)
C     HCNET(NZ,NY,NX)=HCNET(NZ,NY,NX)+GZ*(WTHTRT-WTHTHT)
      XHVSTC(NY,NX)=XHVSTC(NY,NX)+GY*(WTHTHT-WTHTRT)
      XHVSTN(NY,NX)=XHVSTN(NY,NX)+WTHNHT-WTHNRT
      XHVSTP(NY,NX)=XHVSTP(NY,NX)+WTHPHT-WTHPRT
      RECO(NY,NX)=RECO(NY,NX)-GZ*(WTHTHT-WTHTRT)
      TRAU(NY,NX)=TRAU(NY,NX)-GZ*(WTHTHT-WTHTRT)
C     WRITE(*,6542)'GRAZ',I,J,NFZ,NX,NY,NZ,HVSTC(NZ,NY,NX)
C    2,GY,GZ,WTHTHT,WTHTRT
      ENDIF
C
C     ABOVE-GROUND LITTERFALL FROM HARVESTING
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     CSNC,ZSNC,PSNC=C,N,P litterfall from disturbance
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTHTR0,WTHNR0,WTHPR0=nonstructural C,N,P to litter
C     WTHTR1,WTHNR1,WTHPR1=leaf C,N,P to litter
C     WTHTR2,WTHNR2,WTHPR2=fine,non-leaf C,N,P to litter
C     WTHTR3,WTHNR3,WTHPR3=woody C,N,P to litter
C     WTHTR4,WTHNR4,WTHPR4=standing dead C,N,P to litter
C     WTHTX1,WTHNX1,WTHPX1=harvested leaf C,N,P to litter
C     WTHTX2,WTHNX2,WTHPX2=harvested petiole C,N,P to litter
C     WTHTX3,WTHNX3,WTHPX3=harvested woody C,N,P to litter
C     WTHTX4,WTHNX4,WTHPX4=harvested standing dead C,N,P to litter
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C
      IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
C     IF(IHVST(NZ,I,NY,NX).NE.5)THEN
      DO 6375 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+CFOPC(0,M,NZ,NY,NX)*(WTHTR0+WTHTX0)
     3+CFOPC(1,M,NZ,NY,NX)*(WTHTR1+WTHTX1)
     4+CFOPC(2,M,NZ,NY,NX)*(WTHTR2+WTHTX2)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+CFOPN(0,M,NZ,NY,NX)*(WTHNR0+WTHNX0)
     3+CFOPN(1,M,NZ,NY,NX)*(WTHNR1+WTHNX1)
     4+CFOPN(2,M,NZ,NY,NX)*(WTHNR2+WTHNX2)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+CFOPP(0,M,NZ,NY,NX)*(WTHPR0+WTHPX0)
     3+CFOPP(1,M,NZ,NY,NX)*(WTHPR1+WTHPX1)
     4+CFOPP(2,M,NZ,NY,NX)*(WTHPR2+WTHPX2)
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+CFOPC(3,M,NZ,NY,NX)*(WTHTR3+WTHTX3+WTHTR4+WTHTX4)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+CFOPN(3,M,NZ,NY,NX)*(WTHNR3+WTHNX3+WTHNR4+WTHNX4)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+CFOPP(3,M,NZ,NY,NX)*(WTHPR3+WTHPX3+WTHPR4+WTHPX4)
      ELSE
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)
     2+CFOPC(5,M,NZ,NY,NX)*(WTHTX3+WTHTX4)
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)
     2+CFOPN(5,M,NZ,NY,NX)*(WTHNX3+WTHNX4)
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)
     2+CFOPP(5,M,NZ,NY,NX)*(WTHPX3+WTHPX4)
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)
     2+CFOPC(5,M,NZ,NY,NX)*(WTHTR3+WTHTR4)*FWOOD(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)
     2+CFOPN(5,M,NZ,NY,NX)*(WTHNR3+WTHNR4)*FWOODN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)
     2+CFOPP(5,M,NZ,NY,NX)*(WTHPR3+WTHPR4)*FWOODP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+CFOPC(2,M,NZ,NY,NX)*(WTHTR3+WTHTR4)*FWOOD(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+CFOPN(2,M,NZ,NY,NX)*(WTHNR3+WTHNR4)*FWOODN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+CFOPP(2,M,NZ,NY,NX)*(WTHPR3+WTHPR4)*FWOODP(1,NZ)
      ENDIF
6375  CONTINUE
C     ENDIF
      ELSE
C
C     ABOVE-GROUND LITTERFALL FROM GRAZING
C
C     HCSNC,HZSNC,HPSNC=hourly C,N,P litterfall
C     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
C     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
C
      HCSNC(NZ,NY,NX)=HCSNC(NZ,NY,NX)+WTHTRT+WTHTXT
      HZSNC(NZ,NY,NX)=HZSNC(NZ,NY,NX)+WTHNRT+WTHNXT
      HPSNC(NZ,NY,NX)=HPSNC(NZ,NY,NX)+WTHPRT+WTHPXT
      TCSNC(NZ,NY,NX)=TCSNC(NZ,NY,NX)+WTHTRT+WTHTXT
      TZSNC(NZ,NY,NX)=TZSNC(NZ,NY,NX)+WTHNRT+WTHNXT
      TPSNC(NZ,NY,NX)=TPSNC(NZ,NY,NX)+WTHPRT+WTHPXT
      TCSN0(NZ,NY,NX)=TCSN0(NZ,NY,NX)+WTHTRT+WTHTXT
      TZSN0(NZ,NY,NX)=TZSNC(NZ,NY,NX)+WTHNRT+WTHNXT
      TPSN0(NZ,NY,NX)=TPSNC(NZ,NY,NX)+WTHPRT+WTHPXT
C     WRITE(*,4812)'GRAZ',I,J,NFZ,NX,NY,NZ,NB,IHVST(NZ,I,NY,NX)
C    2,TCSN0(NZ,NY,NX),WTHTRT,WTHTXT
4812  FORMAT(A8,8I4,12E12.4)
C
C     ADD MANURE FROM GRAZING TO SURFACE LITTERFALL
C
C     RUMINANT
C
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C     CMOSC=manure protein(1),carbohydrate(2),cellulose(3),
C        lignin(4) fractions
C     CSNM,ZSNM,PSNM=manure organic C,N,P deposition
C     ZSNI,PSNI=manure inorganic N,P deposition
C     WTHTRT,WTHNRT,WTHPRT=total removed C,N,P to litter
C     UORGF,UFERTN,UFERTP=accumulated litter C,N,P application
C
      IF(IHVST(NZ,I,NY,NX).EQ.4)THEN
      CMOSC(1)=0.036
      CMOSC(2)=0.044
      CMOSC(3)=0.630
      CMOSC(4)=0.290
C
C     NON-RUMINANT
C
      ELSEIF(IHVST(NZ,I,NY,NX).EQ.6)THEN
      CMOSC(1)=0.138
      CMOSC(2)=0.401
      CMOSC(3)=0.316
      CMOSC(4)=0.145
      ENDIF
      DO 6500 M=1,4
      CSNM(M,NZ,NY,NX)=CMOSC(M)*(WTHTRT+WTHTXT)
      ZSNM(M,NZ,NY,NX)=CMOSC(M)*(WTHNRT+WTHNXT)*0.5
      PSNM(M,NZ,NY,NX)=CMOSC(M)*(WTHPRT+WTHPXT)*0.5
C     IF(NX.EQ.2)THEN
C     WRITE(*,6542)'MANURE',I,J,NFZ,NX,NY,NZ,M,CSNM(M,NZ,NY,NX)
C    2,ZSNM(M,NZ,NY,NX),PSNM(M,NZ,NY,NX),WTHTRT,WTHNRT,WTHPRT
C    3,WTHTXT,WTHNXT,WTHPXT
6542  FORMAT(A8,7I4,20E12.4)
C     ENDIF
6500  CONTINUE
      ZSNI(NZ,NY,NX)=(WTHNRT+WTHNXT)*0.5
      PSNI(NZ,NY,NX)=(WTHPRT+WTHPXT)*0.5
      UORGF(NY,NX)=UORGF(NY,NX)+WTHTRT+WTHTXT
      UFERTN(NY,NX)=UFERTN(NY,NX)+WTHNRT+WTHNXT
      UFERTP(NY,NX)=UFERTP(NY,NX)+WTHPRT+WTHPXT
      ENDIF
      ZEROP(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)
      ZEROQ(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      ZEROP2(NZ,NY,NX)=ZERO2*PP(NZ,NY,NX)
      ENDIF
      ENDIF
C
C     END OF HARVEST
C
C     RESET DEAD BRANCHES
C
C     ZNOON=hour of solar noon
C     IDAY(1,=emergence date
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IDAYH,IYRH=day,year of harvesting
C     IYRC=current year
C     IDTHB=branch living flag: 0=alive,1=dead
C     GROUP=node number required for floral initiation
C     PSTGI=node number at floral initiation
C     PSTGF=node number at flowering
C     VSTG=number of leaves appeared
C     KVSTG=integer of most recent leaf number currently growing
C     VSTGX=leaf number on date of floral initiation
C     TGSTGI=total change in vegetative node number normalized
C        for maturity group
C     TGSTGF=total change in reproductive node number normalized
C        for maturity group
C     FLG4=number of hours with no grain fill
C     IFLGA=flag for initializing leafout
C     VRNS,VRNL=leafout hours,hours required for leafout
C     VRNF,VRNX=leafoff hours,hours required for leafoff
C     ATRP=hourly leafout counter
C     FDBK,FDBKX=N,P feedback inhibition on C3 CO2 fixation
C     IFLGA,IFLGE=flag for initializing,enabling leafout
C     IFLGF=flag for enabling leafoff:0=enable,1=disable
C     IFLGQ=current hours after physlogical maturity
C        until start of litterfall
C
      IF(IDAY(1,NB1(NZ,NY,NX),NZ,NY,NX).NE.0
     3.AND.(ISTYP(NZ,NY,NX).NE.0.OR.(I.GT.IDAYH(NZ,NY,NX)
     4.AND.IYRC.GE.IYRH(NZ,NY,NX))))THEN
      IDTHY=0
C
C     RESET PHENOLOGY AND GROWTH STAGE OF DEAD BRANCHES
C
      DO 8845 NB=1,NBR(NZ,NY,NX)
      IF(IDTHB(NB,NZ,NY,NX).EQ.1)THEN
      GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
      PSTG(NB,NZ,NY,NX)=XTLI(NZ,NY,NX)
      PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
      PSTGF(NB,NZ,NY,NX)=0.0
      VSTG(NB,NZ,NY,NX)=0.0
      VSTGX(NB,NZ,NY,NX)=0.0
      KLEAF(NB,NZ,NY,NX)=1
      KVSTG(NB,NZ,NY,NX)=1
      TGSTGI(NB,NZ,NY,NX)=0.0
      TGSTGF(NB,NZ,NY,NX)=0.0
      VRNS(NB,NZ,NY,NX)=0.0
      VRNF(NB,NZ,NY,NX)=0.0
      VRNY(NB,NZ,NY,NX)=0.0
      VRNZ(NB,NZ,NY,NX)=0.0
      ATRP(NB,NZ,NY,NX)=0.0
      FLG4(NB,NZ,NY,NX)=0.0
      FDBK(NB,NZ,NY,NX)=1.0
      FDBKX(NB,NZ,NY,NX)=1.0
      IFLGA(NB,NZ,NY,NX)=0
      IFLGE(NB,NZ,NY,NX)=1
      IFLGF(NB,NZ,NY,NX)=0
      IFLGR(NB,NZ,NY,NX)=0
      IFLGQ(NB,NZ,NY,NX)=0
      NBTB(NB,NZ,NY,NX)=0
      DO 8850 M=1,10
      IDAY(M,NB,NZ,NY,NX)=0
8850  CONTINUE
C
C     LITTERFALL FROM DEAD BRANCHES
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
C     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
C     WTEARB,WTEABN,WTEABP=ear C,N,P mass
C     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C
      DO 6405 M=1,4
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+CFOPC(0,M,NZ,NY,NX)*CPOLNB(NB,NZ,NY,NX)
     4+CFOPC(1,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(1,NZ)
     5+WTNDB(NB,NZ,NY,NX))
     6+CFOPC(2,M,NZ,NY,NX)*(WTSHEB(NB,NZ,NY,NX)*FWODB(1,NZ)
     7+WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX))
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)
     2+CFOPC(5,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(0,NZ)
     3+WTSHEB(NB,NZ,NY,NX)*FWODB(0,NZ))
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+CFOPN(0,M,NZ,NY,NX)*ZPOLNB(NB,NZ,NY,NX)
     4+CFOPN(1,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(1,NZ)
     5+WTNDBN(NB,NZ,NY,NX))
     6+CFOPN(2,M,NZ,NY,NX)*(WTSHBN(NB,NZ,NY,NX)*FWODSN(1,NZ)
     7+WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX))
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)
     2+CFOPN(5,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(0,NZ)
     3+WTSHBN(NB,NZ,NY,NX)*FWODSN(0,NZ))
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+CFOPP(0,M,NZ,NY,NX)*PPOLNB(NB,NZ,NY,NX)
     4+CFOPP(1,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(1,NZ)
     5+WTNDBP(NB,NZ,NY,NX))
     6+CFOPP(2,M,NZ,NY,NX)*(WTSHBP(NB,NZ,NY,NX)*FWODSP(1,NZ)
     7+WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX))
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)
     2+CFOPP(5,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(0,NZ)
     3+WTSHBP(NB,NZ,NY,NX)*FWODSP(0,NZ))
C
C     IF WINTER ANNUAL
C
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)
     2+CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)
     2+CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)
     2+CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ELSE
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ENDIF
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)
     5*(WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX))
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)
     5*(WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX))
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)
     5*(WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX))
6405  CONTINUE
C
C     RECOVER NON-STRUCTURAL C,N,P FROM BRANCH TO
C     SEASONAL STORAGE RESERVES
C
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
C     CPOOLK=total C4 nonstructural C in branch
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     IHVST=harvest type:0=none,1=grain,2=all above-ground
C                       ,3=pruning,4=grazing,6=herbivory
C
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+CPOOL(NB,NZ,NY,NX)
     2+CPOOLK(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+ZPOOL(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+PPOOL(NB,NZ,NY,NX)
C     IF(IHVST(NZ,I,NY,NX).NE.4.AND.IHVST(NZ,I,NY,NX).NE.6)THEN
C     DO 6406 M=1,4
C     CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
C    2+CFOPC(0,M,NZ,NY,NX)*WTRSVB(NB,NZ,NY,NX)
C     ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
C    2+CFOPN(0,M,NZ,NY,NX)*WTRSBN(NB,NZ,NY,NX)
C     PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
C    2+CFOPP(0,M,NZ,NY,NX)*WTRSBP(NB,NZ,NY,NX)
6406  CONTINUE
C     ELSE
C     WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
C     WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX)
C     WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX)
C     ENDIF
C
C     RESET STATE VARIABLES FROM DEAD BRANCHES
C
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
C     WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
C     WTEARB,WTEABN,WTEABP=ear C,N,P mass
C     WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     GRNOB=seed set number
C     GRNXB=potential number of seed set sites
C     GRWTB=individual seed size
C     CPOOL3,CPOOL4=C4 nonstructural C mass in bundle sheath,mesophyll
C     CO2B,HCOB=aqueous CO2,HCO3-C mass in bundle sheath
C     WSLF=leaf protein mass
C     ARLFB=branch leaf area
C     WGLF,WGLFN,WGLFP,WSLF=node leaf C,N,P,protein mass
C     WGSHE,WGSHN,WGSHP,WSSHE=node petiole C,N,P,protein mass
C     WGNODE,WGNODN,WGNODP=node stalk C,N,P mass
C
      CPOOL(NB,NZ,NY,NX)=0.0
      CPOOLK(NB,NZ,NY,NX)=0.0
      ZPOOL(NB,NZ,NY,NX)=0.0
      PPOOL(NB,NZ,NY,NX)=0.0
      CPOLNB(NB,NZ,NY,NX)=0.0
      ZPOLNB(NB,NZ,NY,NX)=0.0
      PPOLNB(NB,NZ,NY,NX)=0.0
      WTSHTB(NB,NZ,NY,NX)=0.0
      WTLFB(NB,NZ,NY,NX)=0.0
      WTNDB(NB,NZ,NY,NX)=0.0
      WTSHEB(NB,NZ,NY,NX)=0.0
      WTSTKB(NB,NZ,NY,NX)=0.0
      WVSTKB(NB,NZ,NY,NX)=0.0
      WTRSVB(NB,NZ,NY,NX)=0.0
      WTHSKB(NB,NZ,NY,NX)=0.0
      WTEARB(NB,NZ,NY,NX)=0.0
      WTGRB(NB,NZ,NY,NX)=0.0
      WTLSB(NB,NZ,NY,NX)=0.0
      WTSHTN(NB,NZ,NY,NX)=0.0
      WTLFBN(NB,NZ,NY,NX)=0.0
      WTNDBN(NB,NZ,NY,NX)=0.0
      WTSHBN(NB,NZ,NY,NX)=0.0
      WTSTBN(NB,NZ,NY,NX)=0.0
      WTRSBN(NB,NZ,NY,NX)=0.0
      WTHSBN(NB,NZ,NY,NX)=0.0
      WTEABN(NB,NZ,NY,NX)=0.0
      WTGRBN(NB,NZ,NY,NX)=0.0
      WTSHTP(NB,NZ,NY,NX)=0.0
      WTLFBP(NB,NZ,NY,NX)=0.0
      WTNDBP(NB,NZ,NY,NX)=0.0
      WTSHBP(NB,NZ,NY,NX)=0.0
      WTSTBP(NB,NZ,NY,NX)=0.0
      WTRSBP(NB,NZ,NY,NX)=0.0
      WTHSBP(NB,NZ,NY,NX)=0.0
      WTEABP(NB,NZ,NY,NX)=0.0
      WTGRBP(NB,NZ,NY,NX)=0.0
      GRNXB(NB,NZ,NY,NX)=0.0
      GRNOB(NB,NZ,NY,NX)=0.0
      GRWTB(NB,NZ,NY,NX)=0.0
      ARLFB(NB,NZ,NY,NX)=0.0
      WTSTXB(NB,NZ,NY,NX)=0.0
      WTSTXN(NB,NZ,NY,NX)=0.0
      WTSTXP(NB,NZ,NY,NX)=0.0
      DO 8855 K=0,25
      IF(K.NE.0)THEN
      CPOOL3(K,NB,NZ,NY,NX)=0.0
      CPOOL4(K,NB,NZ,NY,NX)=0.0
      CO2B(K,NB,NZ,NY,NX)=0.0
      HCOB(K,NB,NZ,NY,NX)=0.0
      ENDIF
      ARLF(K,NB,NZ,NY,NX)=0.0
      HTNODE(K,NB,NZ,NY,NX)=0.0
      HTNODX(K,NB,NZ,NY,NX)=0.0
      HTSHE(K,NB,NZ,NY,NX)=0.0
      WGLF(K,NB,NZ,NY,NX)=0.0
      WSLF(K,NB,NZ,NY,NX)=0.0
      WGLFN(K,NB,NZ,NY,NX)=0.0
      WGLFP(K,NB,NZ,NY,NX)=0.0
      WGSHE(K,NB,NZ,NY,NX)=0.0
      WSSHE(K,NB,NZ,NY,NX)=0.0
      WGSHN(K,NB,NZ,NY,NX)=0.0
      WGSHP(K,NB,NZ,NY,NX)=0.0
      WGNODE(K,NB,NZ,NY,NX)=0.0
      WGNODN(K,NB,NZ,NY,NX)=0.0
      WGNODP(K,NB,NZ,NY,NX)=0.0
      DO 8865 L=1,JC
      ARLFV(L,NZ,NY,NX)=ARLFV(L,NZ,NY,NX)-ARLFL(L,K,NB,NZ,NY,NX)
      WGLFV(L,NZ,NY,NX)=WGLFV(L,NZ,NY,NX)-WGLFL(L,K,NB,NZ,NY,NX)
      ARLFL(L,K,NB,NZ,NY,NX)=0.0
      WGLFL(L,K,NB,NZ,NY,NX)=0.0
      WGLFLN(L,K,NB,NZ,NY,NX)=0.0
      WGLFLP(L,K,NB,NZ,NY,NX)=0.0
      IF(K.NE.0)THEN
      DO 8860 N=1,4
      SURF(N,L,K,NB,NZ,NY,NX)=0.0
8860  CONTINUE
      ENDIF
8865  CONTINUE
8855  CONTINUE
      DO 8875 L=1,JC
      ARSTK(L,NB,NZ,NY,NX)=0.0
      DO 8875 N=1,4
      SURFB(N,L,NB,NZ,NY,NX)=0.0
8875  CONTINUE
      IDTHY=IDTHY+1
      ENDIF
8845  CONTINUE
      IF(IDTHY.EQ.NBR(NZ,NY,NX))THEN
      IDTHP(NZ,NY,NX)=1
      IDTHR(NZ,NY,NX)=1
      NBT(NZ,NY,NX)=0
      WSTR(NZ,NY,NX)=0.0
C
C     IF WINTER ANNUAL
C
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      NBR(NZ,NY,NX)=1
      ELSE
      NBR(NZ,NY,NX)=0
      PPX(NZ,NY,NX)=0.0
      PP(NZ,NY,NX)=0.0
      ENDIF
      HTCTL(NZ,NY,NX)=0.0
C     VOLWOU=VOLWOU+VOLWP(NZ,NY,NX)
C     UVOLO(NY,NX)=UVOLO(NY,NX)+VOLWP(NZ,NY,NX)
C     VOLWP(NZ,NY,NX)=0.0
C
C     RESET LIVING FLAGS
C
C     WTRVC,WTRT=PFT storage,root C
C     ISTYP=growth habit:0=annual,1=perennial
C     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
C     PP=PFT population
C     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
C
C     IF(WTRVC(NZ,NY,NX).LT.1.0E-04*WTRT(NZ,NY,NX)
C    2.AND.ISTYP(NZ,NY,NX).NE.0)IDTHR(NZ,NY,NX)=1
C     IF(ISTYP(NZ,NY,NX).EQ.0)IDTHR(NZ,NY,NX)=1
C     IF(JHVST(NZ,I,NY,NX).NE.0)IDTHR(NZ,NY,NX)=1
C     IF(PP(NZ,NY,NX).LE.0.0)IDTHR(NZ,NY,NX)=1
C     IF(IDTHP(NZ,NY,NX).EQ.1)IDTHR(NZ,NY,NX)=1
      ENDIF
C
C     DEAD ROOTS
C
C     LITTERFALL FROM DEAD ROOTS
C
C     IDTHR=PFT root living flag: 0=alive,1=dead
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,
C        1=non-woody
C     CO2A,OXYA,CH4A,Z2OA,ZH3A,H2GA=root gaseous CO2,O2,CH4,N2O,NH3,H2
C     CO2P,OXYP,CH4P,Z2OP,ZH3P,H2GP=root aqueous CO2,O2,CH4,N2O,NH3,H2
C     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=root gaseous
C        CO2,O2,CH4,N2O,NH3,H2 loss from disturbance
C
      IF(IDTHR(NZ,NY,NX).EQ.1)THEN
      DO 8900 N=1,MY(NZ,NY,NX)
      DO 8895 L=NU(NY,NX),NJ(NY,NX)
      DO 6410 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(0,M,NZ,NY,NX)
     2*CPOOLR(N,L,NZ,NY,NX)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(0,M,NZ,NY,NX)
     2*ZPOOLR(N,L,NZ,NY,NX)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(0,M,NZ,NY,NX)
     2*PPOOLR(N,L,NZ,NY,NX)
      DO 6410 NR=1,NRT(NZ,NY,NX)
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*(WTRT1(N,L,NR,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0,NZ)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*(WTRT1N(N,L,NR,NZ,NY,NX)+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0,NZ)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*(WTRT1P(N,L,NR,NZ,NY,NX)+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0,NZ)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX)
     2*(WTRT1(N,L,NR,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1,NZ)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX)
     2*(WTRT1N(N,L,NR,NZ,NY,NX)+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1,NZ)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX)
     2*(WTRT1P(N,L,NR,NZ,NY,NX)+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1,NZ)
6410  CONTINUE
C
C     RELEASE GAS CONTENTS OF DEAD ROOTS
C
      RCO2Z(NZ,NY,NX)=RCO2Z(NZ,NY,NX)-CO2A(N,L,NZ,NY,NX)
     2-CO2P(N,L,NZ,NY,NX)
      ROXYZ(NZ,NY,NX)=ROXYZ(NZ,NY,NX)-OXYA(N,L,NZ,NY,NX)
     2-OXYP(N,L,NZ,NY,NX)
      RCH4Z(NZ,NY,NX)=RCH4Z(NZ,NY,NX)-CH4A(N,L,NZ,NY,NX)
     2-CH4P(N,L,NZ,NY,NX)
      RN2OZ(NZ,NY,NX)=RN2OZ(NZ,NY,NX)-Z2OA(N,L,NZ,NY,NX)
     2-Z2OP(N,L,NZ,NY,NX)
      RNH3Z(NZ,NY,NX)=RNH3Z(NZ,NY,NX)-ZH3A(N,L,NZ,NY,NX)
     2-ZH3P(N,L,NZ,NY,NX)
      RH2GZ(NZ,NY,NX)=RH2GZ(NZ,NY,NX)-H2GA(N,L,NZ,NY,NX)
     2-H2GP(N,L,NZ,NY,NX)
      CO2A(N,L,NZ,NY,NX)=0.0
      OXYA(N,L,NZ,NY,NX)=0.0
      CH4A(N,L,NZ,NY,NX)=0.0
      Z2OA(N,L,NZ,NY,NX)=0.0
      ZH3A(N,L,NZ,NY,NX)=0.0
      H2GA(N,L,NZ,NY,NX)=0.0
      CO2P(N,L,NZ,NY,NX)=0.0
      OXYP(N,L,NZ,NY,NX)=0.0
      CH4P(N,L,NZ,NY,NX)=0.0
      Z2OP(N,L,NZ,NY,NX)=0.0
      ZH3P(N,L,NZ,NY,NX)=0.0
      H2GP(N,L,NZ,NY,NX)=0.0
C
C     RESET STATE VARIABLES OF DEAD ROOTS
C
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
C     RTLG1,RTLG2=primary,secondary root length
C     RTN2=number of secondary root axes
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     WTRTL,WTRTD=active,actual root C mass
C     WSRTL=root protein C mass
C     RTN1,RTNL=number of primary,secondary root axes
C     RTDNP,RTLGP=root length density,root length per plant
C     RTVLW,RTVLP=root or myco aqueous,gaseous volume
C     RTARP=root surface area per plant
C     RTLGA=average secondary root length
C
      DO 8870 NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=0.0
      WTRT1N(N,L,NR,NZ,NY,NX)=0.0
      WTRT1P(N,L,NR,NZ,NY,NX)=0.0
      WTRT2(N,L,NR,NZ,NY,NX)=0.0
      WTRT2N(N,L,NR,NZ,NY,NX)=0.0
      WTRT2P(N,L,NR,NZ,NY,NX)=0.0
      RTWT1(N,NR,NZ,NY,NX)=0.0
      RTWT1N(N,NR,NZ,NY,NX)=0.0
      RTWT1P(N,NR,NZ,NY,NX)=0.0
      RTLG1(N,L,NR,NZ,NY,NX)=0.0
      RTLG2(N,L,NR,NZ,NY,NX)=0.0
      RTN2(N,L,NR,NZ,NY,NX)=0.0
8870  CONTINUE
      CPOOLR(N,L,NZ,NY,NX)=0.0
      ZPOOLR(N,L,NZ,NY,NX)=0.0
      PPOOLR(N,L,NZ,NY,NX)=0.0
      WTRTL(N,L,NZ,NY,NX)=0.0
      WTRTD(N,L,NZ,NY,NX)=0.0
      WSRTL(N,L,NZ,NY,NX)=0.0
      RTN1(N,L,NZ,NY,NX)=0.0
      RTNL(N,L,NZ,NY,NX)=0.0
      RTLGP(N,L,NZ,NY,NX)=0.0
      RTDNP(N,L,NZ,NY,NX)=0.0
      RTVLP(N,L,NZ,NY,NX)=0.0
      RTVLW(N,L,NZ,NY,NX)=0.0
      RRAD1(N,L,NZ,NY,NX)=RRAD1M(N,NZ,NY,NX)
      RRAD2(N,L,NZ,NY,NX)=RRAD2M(N,NZ,NY,NX)
      RTARP(N,L,NZ,NY,NX)=0.0
      RTLGA(N,L,NZ,NY,NX)=RTLGAX
C
C     LITTERFALL AND STATE VARIABLES FROM DEAD NODULES
C
C     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
C     CSNC,ZSNC,PSNC=C,N,P litterfall from decomposition
C        and senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C     CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C
      IF(INTYP(NZ,NY,NX).NE.0.AND.N.EQ.1)THEN
      DO 6420 M=1,4
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX)
     2*WTNDL(L,NZ,NY,NX)+CFOPC(0,M,NZ,NY,NX)*CPOOLN(L,NZ,NY,NX)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX)
     2*WTNDLN(L,NZ,NY,NX)+CFOPN(0,M,NZ,NY,NX)*ZPOOLN(L,NZ,NY,NX)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX)
     2*WTNDLP(L,NZ,NY,NX)+CFOPP(0,M,NZ,NY,NX)*PPOOLN(L,NZ,NY,NX)
6420  CONTINUE
      WTNDL(L,NZ,NY,NX)=0.0
      WTNDLN(L,NZ,NY,NX)=0.0
      WTNDLP(L,NZ,NY,NX)=0.0
      CPOOLN(L,NZ,NY,NX)=0.0
      ZPOOLN(L,NZ,NY,NX)=0.0
      PPOOLN(L,NZ,NY,NX)=0.0
      ENDIF
8895  CONTINUE
8900  CONTINUE
C
C     RESET DEPTH VARIABLES OF DEAD ROOTS
C
C     NINR=deepest root layer
C     RTDP1=primary root depth from soil surface
C     RTWT1,RTWT1N,RTWT1P=primary root C,N,P mass
C
      DO 8795 NR=1,NRT(NZ,NY,NX)
      NINR(NR,NZ,NY,NX)=NG(NZ,NY,NX)
      DO 8790 N=1,MY(NZ,NY,NX)
      RTDP1(N,NR,NZ,NY,NX)=SDPTH(NZ,NY,NX)
      RTWT1(N,NR,NZ,NY,NX)=0.0
      RTWT1N(N,NR,NZ,NY,NX)=0.0
      RTWT1P(N,NR,NZ,NY,NX)=0.0
8790  CONTINUE
8795  CONTINUE
      NIX(NZ,NY,NX)=NG(NZ,NY,NX)
      NRT(NZ,NY,NX)=0
      ENDIF
C
C     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE
C     RESERVES FROM SHOOT AT DEATH
C
C     IDTHP,IDTHR=PFT shoot,root living flag: 0=alive,1=dead
C     IFLGI=PFT initialization flag:0=no,1=yes
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P in branch
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     CPOOLK=total C4 nonstructural C in branch
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     FWODB=C woody fraction in other organs:0=woody,1=non-woody
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
C     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
C     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
C     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
C
      IF(IDTHP(NZ,NY,NX).EQ.1.AND.IDTHR(NZ,NY,NX).EQ.1)THEN
      IF(ISTYP(NZ,NY,NX).NE.0.OR.IWTYP(NZ,NY,NX).EQ.0)THEN
      IDTH(NZ,NY,NX)=1
      DO 6426 M=1,4
      CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)
     2+CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX)*FWOOD(0,NZ)
      ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)
     2+CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX)*FWOODN(0,NZ)
      PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,0,NG(NZ,NY,NX),NZ,NY,NX)
     2+CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX)*FWOODP(0,NZ)
      CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=CSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)
     2+CFOPC(0,M,NZ,NY,NX)*WTRVC(NZ,NY,NX)*FWOOD(1,NZ)
      ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=ZSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)
     2+CFOPN(0,M,NZ,NY,NX)*WTRVN(NZ,NY,NX)*FWOODN(1,NZ)
      PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)=PSNC(M,1,NG(NZ,NY,NX),NZ,NY,NX)
     2+CFOPP(0,M,NZ,NY,NX)*WTRVP(NZ,NY,NX)*FWOODP(1,NZ)
6426  CONTINUE
      WTRVC(NZ,NY,NX)=0.0
      WTRVN(NZ,NY,NX)=0.0
      WTRVP(NZ,NY,NX)=0.0
      ENDIF
      DO 6425 M=1,4
      DO 8825 NB=1,NBR(NZ,NY,NX)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+CFOPC(0,M,NZ,NY,NX)*(CPOOL(NB,NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX)
     3+CPOOLK(NB,NZ,NY,NX))
     4+CFOPC(1,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(1,NZ)
     5+WTNDB(NB,NZ,NY,NX))
     6+CFOPC(2,M,NZ,NY,NX)*(WTSHEB(NB,NZ,NY,NX)*FWODB(1,NZ)
     7+WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX))
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)
     2+CFOPC(5,M,NZ,NY,NX)*(WTLFB(NB,NZ,NY,NX)*FWODB(0,NZ)
     3+WTSHEB(NB,NZ,NY,NX)*FWODB(0,NZ))
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+CFOPN(0,M,NZ,NY,NX)*(ZPOOL(NB,NZ,NY,NX)+ZPOLNB(NB,NZ,NY,NX))
     4+CFOPN(1,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(1,NZ)
     5+WTNDBN(NB,NZ,NY,NX))
     6+CFOPN(2,M,NZ,NY,NX)*(WTSHBN(NB,NZ,NY,NX)*FWODSN(1,NZ)
     7+WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX))
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)
     2+CFOPN(5,M,NZ,NY,NX)*(WTLFBN(NB,NZ,NY,NX)*FWODLN(0,NZ)
     3+WTSHBN(NB,NZ,NY,NX)*FWODSN(0,NZ))
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+CFOPP(0,M,NZ,NY,NX)*(PPOOL(NB,NZ,NY,NX)+PPOLNB(NB,NZ,NY,NX))
     4+CFOPP(1,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(1,NZ)
     5+WTNDBP(NB,NZ,NY,NX))
     6+CFOPP(2,M,NZ,NY,NX)*(WTSHBP(NB,NZ,NY,NX)*FWODSP(1,NZ)
     7+WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX))
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)
     2+CFOPP(5,M,NZ,NY,NX)*(WTLFBP(NB,NZ,NY,NX)*FWODLP(0,NZ)
     3+WTSHBP(NB,NZ,NY,NX)*FWODSP(0,NZ))
C
C     IF WINTER ANNUAL
C
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).NE.0)THEN
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)
     2+CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)
     2+CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)
     2+CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ELSE
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)
     2+CFOPC(2,M,NZ,NY,NX)*WTGRB(NB,NZ,NY,NX)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)
     2+CFOPN(2,M,NZ,NY,NX)*WTGRBN(NB,NZ,NY,NX)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)
     2+CFOPP(2,M,NZ,NY,NX)*WTGRBP(NB,NZ,NY,NX)
      ENDIF
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)+CFOPC(3,M,NZ,NY,NX)
     5*(WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX))
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)+CFOPN(3,M,NZ,NY,NX)
     5*(WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX))
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)+CFOPP(3,M,NZ,NY,NX)
     5*(WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX))
8825  CONTINUE
C
C     LITTERFALL AND STATE VARIABLES FOR SEASONAL STORAGE
C     RESERVES FROM ROOT AND STORGE AT DEATH
C
C     CSNC,ZSNC,PSNC=C,N,P litterfall from senescence
C     CFOPC,CFOPN,CFOPC=fraction of plant litter C,N,P allocated
C        to nonstructural(0,*),foliar(1,*),non-foliar(2,*),
C        stalk(3,*),root(4,*),coarse woody (5,*)
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,
C        1=non-woody
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C
      DO 6415 L=NU(NY,NX),NJ(NY,NX)
      DO 6415 N=1,MY(NZ,NY,NX)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(0,M,NZ,NY,NX)
     2*CPOOLR(N,L,NZ,NY,NX)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(0,M,NZ,NY,NX)
     2*ZPOOLR(N,L,NZ,NY,NX)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(0,M,NZ,NY,NX)
     2*PPOOLR(N,L,NZ,NY,NX)
      DO 6415 NR=1,NRT(NZ,NY,NX)
      CSNC(M,0,L,NZ,NY,NX)=CSNC(M,0,L,NZ,NY,NX)+CFOPC(5,M,NZ,NY,NX)
     2*(WTRT1(N,L,NR,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(0,NZ)
      ZSNC(M,0,L,NZ,NY,NX)=ZSNC(M,0,L,NZ,NY,NX)+CFOPN(5,M,NZ,NY,NX)
     2*(WTRT1N(N,L,NR,NZ,NY,NX)+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(0,NZ)
      PSNC(M,0,L,NZ,NY,NX)=PSNC(M,0,L,NZ,NY,NX)+CFOPP(5,M,NZ,NY,NX)
     2*(WTRT1P(N,L,NR,NZ,NY,NX)+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(0,NZ)
      CSNC(M,1,L,NZ,NY,NX)=CSNC(M,1,L,NZ,NY,NX)+CFOPC(4,M,NZ,NY,NX)
     2*(WTRT1(N,L,NR,NZ,NY,NX)+WTRT2(N,L,NR,NZ,NY,NX))*FWODR(1,NZ)
      ZSNC(M,1,L,NZ,NY,NX)=ZSNC(M,1,L,NZ,NY,NX)+CFOPN(4,M,NZ,NY,NX)
     2*(WTRT1N(N,L,NR,NZ,NY,NX)+WTRT2N(N,L,NR,NZ,NY,NX))*FWODRN(1,NZ)
      PSNC(M,1,L,NZ,NY,NX)=PSNC(M,1,L,NZ,NY,NX)+CFOPP(4,M,NZ,NY,NX)
     2*(WTRT1P(N,L,NR,NZ,NY,NX)+WTRT2P(N,L,NR,NZ,NY,NX))*FWODRP(1,NZ)
6415  CONTINUE
6425  CONTINUE
C
C     RESET BRANCH STATE VARIABLES
C
      DO 8835 NB=1,NBR(NZ,NY,NX)
      CPOOL(NB,NZ,NY,NX)=0.0
      CPOOLK(NB,NZ,NY,NX)=0.0
      ZPOOL(NB,NZ,NY,NX)=0.0
      PPOOL(NB,NZ,NY,NX)=0.0
      CPOLNB(NB,NZ,NY,NX)=0.0
      ZPOLNB(NB,NZ,NY,NX)=0.0
      PPOLNB(NB,NZ,NY,NX)=0.0
      WTSHTB(NB,NZ,NY,NX)=0.0
      WTLFB(NB,NZ,NY,NX)=0.0
      WTNDB(NB,NZ,NY,NX)=0.0
      WTSHEB(NB,NZ,NY,NX)=0.0
      WTSTKB(NB,NZ,NY,NX)=0.0
      WVSTKB(NB,NZ,NY,NX)=0.0
      WTRSVB(NB,NZ,NY,NX)=0.0
      WTHSKB(NB,NZ,NY,NX)=0.0
      WTEARB(NB,NZ,NY,NX)=0.0
      WTGRB(NB,NZ,NY,NX)=0.0
      WTLSB(NB,NZ,NY,NX)=0.0
      WTSHTN(NB,NZ,NY,NX)=0.0
      WTLFBN(NB,NZ,NY,NX)=0.0
      WTNDBN(NB,NZ,NY,NX)=0.0
      WTSHBN(NB,NZ,NY,NX)=0.0
      WTSTBN(NB,NZ,NY,NX)=0.0
      WTRSBN(NB,NZ,NY,NX)=0.0
      WTHSBN(NB,NZ,NY,NX)=0.0
      WTEABN(NB,NZ,NY,NX)=0.0
      WTGRBN(NB,NZ,NY,NX)=0.0
      WTSHTP(NB,NZ,NY,NX)=0.0
      WTLFBP(NB,NZ,NY,NX)=0.0
      WTNDBP(NB,NZ,NY,NX)=0.0
      WTSHBP(NB,NZ,NY,NX)=0.0
      WTSTBP(NB,NZ,NY,NX)=0.0
      WTRSBP(NB,NZ,NY,NX)=0.0
      WTHSBP(NB,NZ,NY,NX)=0.0
      WTEABP(NB,NZ,NY,NX)=0.0
      WTGRBP(NB,NZ,NY,NX)=0.0
      WTSTXB(NB,NZ,NY,NX)=0.0
      WTSTXN(NB,NZ,NY,NX)=0.0
      WTSTXP(NB,NZ,NY,NX)=0.0
8835  CONTINUE
C
C     RESET ROOT STATE VARIABLES
C
      DO 6416 L=NU(NY,NX),NJ(NY,NX)
      DO 6416 N=1,MY(NZ,NY,NX)
      CPOOLR(N,L,NZ,NY,NX)=0.0
      ZPOOLR(N,L,NZ,NY,NX)=0.0
      PPOOLR(N,L,NZ,NY,NX)=0.0
      DO 6416 NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L,NR,NZ,NY,NX)=0.0
      WTRT1N(N,L,NR,NZ,NY,NX)=0.0
      WTRT1P(N,L,NR,NZ,NY,NX)=0.0
      WTRT2(N,L,NR,NZ,NY,NX)=0.0
      WTRT2N(N,L,NR,NZ,NY,NX)=0.0
      WTRT2P(N,L,NR,NZ,NY,NX)=0.0
      RTWT1(N,NR,NZ,NY,NX)=0.0
      RTWT1N(N,NR,NZ,NY,NX)=0.0
      RTWT1P(N,NR,NZ,NY,NX)=0.0
      RTLG1(N,L,NR,NZ,NY,NX)=0.0
      RTLG2(N,L,NR,NZ,NY,NX)=0.0
      RTN2(N,L,NR,NZ,NY,NX)=0.0
6416  CONTINUE
C
C     RESEED DEAD PERENNIALS
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     JHVST=terminate PFT:0=no,1=yes,2=yes,but reseed
C     LYRC=number of days in current year
C     IDAY0,IYR0=day,year of planting
C
      IF(ISTYP(NZ,NY,NX).NE.0.AND.JHVST(NZ,I,NY,NX).EQ.0)THEN
      IF(I.LT.LYRC)THEN
      IDAY0(NZ,NY,NX)=I+1
      IYR0(NZ,NY,NX)=IDATA(3)
      ELSE
      IDAY0(NZ,NY,NX)=1
      IYR0(NZ,NY,NX)=IDATA(3)+1
      ENDIF
      IDTH(NZ,NY,NX)=1
      WRITE(*,1114)'RESEED',IYRC,I,J,NFZ,NX,NY,NZ
     2,IDAY0(NZ,NY,NX),IYR0(NZ,NY,NX)
1114  FORMAT(A8,9I4)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
C     ACCUMULATE TOTAL SOIL-PLANT C,N,P EXCHANGE
C
C     HCUPTK,HZUPTK,HPUPTK=net PFT root-soil C,N,P exchange
C     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
C     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
C     UPNF=PFT N2 fixation
C     TCUPTK,TZUPTK,TPUPTK=cumulative PFT root-soil C,N,P exchange
C     TZUPFX=cumulative PFT N2 fixation
C     ZNPP=net primary productivity
C
      HCUPTK(NZ,NY,NX)=UPOMC(NZ,NY,NX)
      HZUPTK(NZ,NY,NX)=UPOMN(NZ,NY,NX)+UPNH4(NZ,NY,NX)
     2+UPNO3(NZ,NY,NX)+UPNF(NZ,NY,NX)
      HPUPTK(NZ,NY,NX)=UPOMP(NZ,NY,NX)+UPH2P(NZ,NY,NX)
     2+UPH1P(NZ,NY,NX)
      TCUPTK(NZ,NY,NX)=TCUPTK(NZ,NY,NX)+UPOMC(NZ,NY,NX)
      TZUPTK(NZ,NY,NX)=TZUPTK(NZ,NY,NX)+UPOMN(NZ,NY,NX)
     2+UPNH4(NZ,NY,NX)+UPNO3(NZ,NY,NX)
      TPUPTK(NZ,NY,NX)=TPUPTK(NZ,NY,NX)+UPOMP(NZ,NY,NX)
     2+UPH2P(NZ,NY,NX)+UPH1P(NZ,NY,NX)
      TZUPFX(NZ,NY,NX)=TZUPFX(NZ,NY,NX)+UPNF(NZ,NY,NX)
     2+UPNFC(NZ,NY,NX)
      ZNPP(NZ,NY,NX)=CARBN(NZ,NY,NX)+TCO2T(NZ,NY,NX)
C
C     TRANSFORMATIONS IN LIVING OR DEAD PLANT POPULATIONS
C
C     TOTAL STANDING DEAD
C
C     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
C     ZG=standing dead height
C     ARSTG,ARSTD=standing dead surface area:total,layer
C     SURFD=standing dead layer surface area accounting for
C        inclination,azimuth used in �hour1.f�
C
      WTSTG(NZ,NY,NX)=0.0
      WTSTGN(NZ,NY,NX)=0.0
      WTSTGP(NZ,NY,NX)=0.0
      DO 6435 M=1,5
      WTSTG(NZ,NY,NX)=WTSTG(NZ,NY,NX)+WTSTDG(M,NZ,NY,NX)
      WTSTGN(NZ,NY,NX)=WTSTGN(NZ,NY,NX)+WTSTDN(M,NZ,NY,NX)
      WTSTGP(NZ,NY,NX)=WTSTGP(NZ,NY,NX)+WTSTDP(M,NZ,NY,NX)
6435  CONTINUE
      IF(WTSTG(NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      VSTD=VSTK*WTSTG(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      RSTD=(0.01*VSTD/3.1416)**0.4
      ZG(NZ,NY,NX)=100.0*RSTD
      ARSTG(NZ,NY,NX)=6.2832*RSTD*ZG(NZ,NY,NX)
      DO 6440 L=1,JC
      IF(ZG(NZ,NY,NX).GT.ZERO.AND.ZL(L-1,NY,NX).LT.ZG(NZ,NY,NX)
     2.AND.ZL(L,NY,NX).GT.ZL(L-1,NY,NX))THEN
      FARSTD=AMIN1(1.0,(ZG(NZ,NY,NX)-ZL(L-1,NY,NX))
     2/(ZL(L,NY,NX)-ZL(L-1,NY,NX)))
      ARSTD(L,NZ,NY,NX)=FARSTD*ARSTG(NZ,NY,NX)
     2*(ZL(L,NY,NX)-ZL(L-1,NY,NX))/AMIN1(ZG(NZ,NY,NX),ZL(JC,NY,NX))
      ELSE
      FARSTD=0.0
      ARSTD(L,NZ,NY,NX)=0.0
      ENDIF
      SURFD(4,L,NZ,NY,NX)=0.25*ARSTD(L,NZ,NY,NX)
C     WRITE(*,4248)'ARSTG',I,J,NFZ,NZ,L
C    2,ARSTG(NZ,NY,NX),ARSTD(L,NZ,NY,NX),FARSTD,ZL(L,NY,NX)
C    3,ZL(L-1,NY,NX),ZT(NY,NX),ZG(NZ,NY,NX),WTSTG(NZ,NY,NX)
C    4/AREA(3,NU(NY,NX),NY,NX),SURFD(4,L,NZ,NY,NX),RSTD
4248  FORMAT(A8,5I4,20E12.4)
6440  CONTINUE
      ELSE
      ZG(NZ,NY,NX)=0.0
      ARSTG(NZ,NY,NX)=0.0
      DO 6441 L=1,JC
      ARSTD(L,NZ,NY,NX)=0.0
      SURFD(4,L,NZ,NY,NX)=0.0
6441  CONTINUE
      ENDIF
9985  CONTINUE
C
C     IF FIRE EVENT IS IN PROGRESS
C
C     ICHKF=fire flag
C
      IF(ICHKF.EQ.1)THEN
C
C     COMBUST CANOPY CH4
C
C     TKQ=canopy air temperature
C     TCMBX=minimum temperature for combustion (K)
C     TFNCO=temperature function for combusting C (600K = 1)
C     RC4CK=CH4 combustion fraction for use in �trnsfr.f�
C
      IF(TKQ(NY,NX).GT.TCMBX)THEN
      RTK=8.3143*TKQ(NY,NX)
      TFNCO=AMIN1(TFNCX,EXP(12.028-60000/RTK))
      RC4CK(NY,NX)=AMIN1(1.0,TFNCO)
      ELSE
      RC4CK(NY,NX)=0.0
      ENDIF
C
C     TOTAL SHOOT C POOLS FOR ALL PFTS USED TO CALCULATE
C     COMBUSTION RATES
C
C     CPOOLP,WTLF,WTSHE,WTSTK,WTRSV,WTHSK,WTEAR,WTGR,CPOLNP
C        WTND,WTSTG=total canopy nonstructural,leaf,petiole,
C        stalk,reserve,husk,ear,grain,bacterial nonstructural,
C        biomass,standing dead C
C
      TWPOOL(NY,NX)=0.0
      TWTLF(NY,NX)=0.0
      TWTSHE(NY,NX)=0.0
      TWTSTK(NY,NX)=0.0
      TWTRSV(NY,NX)=0.0
      TWTHSK(NY,NX)=0.0
      TWTEAR(NY,NX)=0.0
      TWTGR(NY,NX)=0.0
      TWPOLN(NY,NX)=0.0
      TWTNDC(NY,NX)=0.0
      TWTRVC(NY,NX)=0.0
      TWTSTG(NY,NX)=0.0
      DO 2980 L=NU(NY,NX),NJ(NY,NX)
      TWPOLR(L,NY,NX)=0.0
      TWTRTL(L,NY,NX)=0.0
      TWPODL(L,NY,NX)=0.0
      TWTNDL(L,NY,NX)=0.0
2980  CONTINUE
      DO 9965 NZ=1,NP0(NY,NX)
      TWPOOL(NY,NX)=TWPOOL(NY,NX)+CPOOLP(NZ,NY,NX)
      TWTLF(NY,NX)=TWTLF(NY,NX)+WTLF(NZ,NY,NX)
      TWTSHE(NY,NX)=TWTSHE(NY,NX)+WTSHE(NZ,NY,NX)
      TWTSTK(NY,NX)=TWTSTK(NY,NX)+WTSTK(NZ,NY,NX)
      TWTRSV(NY,NX)=TWTRSV(NY,NX)+WTRSV(NZ,NY,NX)
      TWTHSK(NY,NX)=TWTHSK(NY,NX)+WTHSK(NZ,NY,NX)
      TWTEAR(NY,NX)=TWTEAR(NY,NX)+WTEAR(NZ,NY,NX)
      TWTGR(NY,NX)=TWTGR(NY,NX)+WTGR(NZ,NY,NX)
      IF(INTYP(NZ,NY,NX).GE.4)THEN
      TWPOLN(NY,NX)=TWPOLN(NY,NX)+CPOLNP(NZ,NY,NX)
      TWTNDC(NY,NX)=TWTNDC(NY,NX)+WTND(NZ,NY,NX)
      ENDIF
      TWTSTG(NY,NX)=TWTSTG(NY,NX)+WTSTG(NZ,NY,NX)
C
C     TOTAL ROOT C POOLS FOR ALL PFTS
C
C     WTRVC,CPOOLN,WTNDL,CPOOLR,WTRT1,WTRT2=storage,bacterial
C        nonstructural,biomass,root nonstructural,primary,secondary C
C
      TWTRVC(NY,NX)=TWTRVC(NY,NX)+WTRVC(NZ,NY,NX)
      DO 2975 L=NU(NY,NX),NJ(NY,NX)
      TWPODL(L,NY,NX)=TWPODL(L,NY,NX)+CPOOLN(L,NZ,NY,NX)
      TWTNDL(L,NY,NX)=TWTNDL(L,NY,NX)+WTNDL(L,NZ,NY,NX)
      DO 2975 N=1,MY(NZ,NY,NX)
      TWPOLR(L,NY,NX)=TWPOLR(L,NY,NX)+CPOOLR(N,L,NZ,NY,NX)
      DO 2975 NR=1,NRT(NZ,NY,NX)
      TWTRTL(L,NY,NX)=TWTRTL(L,NY,NX)+WTRT1(N,L,NR,NZ,NY,NX)
     2+WTRT2(N,L,NR,NZ,NY,NX)
2975  CONTINUE
9965  CONTINUE
C
C     COMBUST SHOOT
C
C     TKC,TKD=canopy,standing dead surface temperature
C     TCMBX=minimum temperature for combustion (K)
C     TFNCO,TFNDO=temperature function for combusting living,
C        standing dead C (600K = 1)
C     SPCMB=specific combustion rate of different combustion types
C        living plant material (0-4) at 600K, charcoal(5) at 700K
C     RCBC*=combustion rate for each combustion type
C
      DO 9970 NZ=1,NP0(NY,NX)
      IF(TKC(NZ,NY,NX).GT.TCMBX.OR.TKD(NZ,NY,NX).GT.TCMBX)THEN
      RTK=8.3143*TKC(NZ,NY,NX)
      RTKD=8.3143*TKD(NZ,NY,NX)
      TFNCO=AMIN1(TFNCX,EXP(12.028-60000/RTK))
      TFNDO=AMIN1(TFNCX,EXP(12.028-60000/RTKD))
      RCBCX=TFNCO*AREA(3,NU(NY,NX),NY,NX)*XNFH
      RCBDX=TFNDO*AREA(3,NU(NY,NX),NY,NX)*XNFH
      RCBC1=SPCMB(1)*RCBCX
      RCBC2=SPCMB(2)*RCBCX
      RCBC3=SPCMB(3)*RCBCX
      RCBC4=SPCMB(4)*RCBCX
      RCBC5=SPCMB(5)*RCBDX
C
C     SHOOT C POOL COMBUSTION FRACTIONS
C
C     FW*=fraction of C combusted for each total C pool
C     R*=C,N,P combustion rate for each PFT C pool
C     PFT pool code:
C        CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
C        WTLFB,WTLFBN,WTLFBP=leaf C,N,P mass
C        WTSHEB,WTSHBN,WTSHBP=petiole C,N,P mass
C        WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
C        WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C        WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
C        WTEARB,WTEABN,WTEABP=ear C,N,P mass
C        WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
C        CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C        WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C
      IF(TWPOOL(NY,NX).GT.ZEROS(NY,NX))THEN
      FWPOOL=AMIN1(1.0,RCBC1/TWPOOL(NY,NX))
      ELSE
      FWPOOL=0.0
      ENDIF
      IF(TWTLF(NY,NX).GT.ZEROS(NY,NX))THEN
      FWTLF=AMIN1(1.0,RCBC1/TWTLF(NY,NX))
      ELSE
      FWTLF=0.0
      ENDIF
      IF(TWTSHE(NY,NX).GT.ZEROS(NY,NX))THEN
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      FWTSHE=AMIN1(1.0,RCBC2/TWTSHE(NY,NX))
      ELSE
      FWTSHE=AMIN1(1.0,RCBC3/TWTSHE(NY,NX))
      ENDIF
      ELSE
      FWTSHE=0.0
      ENDIF
      IF(TWTSTK(NY,NX).GT.ZEROS(NY,NX))THEN
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      FWTSTK=AMIN1(1.0,RCBC2/TWTSTK(NY,NX))
      ELSE
      FWTSTK=AMIN1(1.0,RCBC4/TWTSTK(NY,NX))
      ENDIF
      ELSE
      FWTSTK=0.0
      ENDIF
      FWTRSV=FWTSTK
      IF(TWTHSK(NY,NX).GT.ZEROS(NY,NX))THEN
      FWTHSK=AMIN1(1.0,RCBC1/TWTHSK(NY,NX))
      ELSE
      FWTHSK=0.0
      ENDIF
      IF(TWTEAR(NY,NX).GT.ZEROS(NY,NX))THEN
      FWTEAR=AMIN1(1.0,RCBC2/TWTEAR(NY,NX))
      ELSE
      FWTEAR=0.0
      ENDIF
      IF(TWTGR(NY,NX).GT.ZEROS(NY,NX))THEN
      FWTGR=AMIN1(1.0,RCBC2/TWTGR(NY,NX))
      ELSE
      FWTGR=0.0
      ENDIF
      IF(TWPOLN(NY,NX).GT.ZEROS(NY,NX))THEN
      FWPOLN=AMIN1(1.0,RCBC1/TWPOLN(NY,NX))
      ELSE
      FWPOLN=0.0
      ENDIF
      IF(TWTNDC(NY,NX).GT.ZEROS(NY,NX))THEN
      FWTNDC=AMIN1(1.0,RCBC2/TWTNDC(NY,NX))
      ELSE
      FWTNDC=0.0
      ENDIF
      IF(TWTSTG(NY,NX).GT.ZEROS(NY,NX))THEN
      FWTSTG=AMIN1(1.0,RCBC5/TWTSTG(NY,NX))
      ELSE
      FWTSTG=0.0
      ENDIF
C
C     SHOOT COMBUSTION RATES
C
C     R*=C,N,P combustion rate for each PFT C pool
C     PFT pool code:
C        CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass
C        WTLFB,WTLFBN,WTLFBP=leaf C,N,P mass
C        WTSHEB,WTSHBN,WTSHBP=petiole C,N,P mass
C        WTSTKB,WTSTBN,WTSTBP=stalk C,N,P mass
C        WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C        WTHSKB,WTHSBN,WTHSBP=husk C,N,P mass
C        WTEARB,WTEABN,WTEABP=ear C,N,P mass
C        WTGRB,WTGRBN,WTGRBP=grain C,N,P mass
C        CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C        WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     ARLFB=leaf area
C
      DO 8750 NB=1,NBR(NZ,NY,NX)
      RCPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)*FWPOOL
      RZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)*FWPOOL
      RPPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)*FWPOOL
      RWTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)*FWTLF
      RWTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)*FWTLF
      RWTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)*FWTLF
      RWTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX)*FWTSHE
      RWTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX)*FWTSHE
      RWTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX)*FWTSHE
      RWTSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)*FWTSTK
      RWTSTBN(NB,NZ,NY,NX)=WTSTBN(NB,NZ,NY,NX)*FWTSTK
      RWTSTBP(NB,NZ,NY,NX)=WTSTBP(NB,NZ,NY,NX)*FWTSTK
      RWTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)*FWTRSV
      RWTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)*FWTRSV
      RWTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)*FWTRSV
      RWTHSKB(NB,NZ,NY,NX)=WTHSKB(NB,NZ,NY,NX)*FWTHSK
      RWTHSBN(NB,NZ,NY,NX)=WTHSBN(NB,NZ,NY,NX)*FWTHSK
      RWTHSBP(NB,NZ,NY,NX)=WTHSBP(NB,NZ,NY,NX)*FWTHSK
      RWTEARB(NB,NZ,NY,NX)=WTEARB(NB,NZ,NY,NX)*FWTEAR
      RWTEABN(NB,NZ,NY,NX)=WTEABN(NB,NZ,NY,NX)*FWTEAR
      RWTEABP(NB,NZ,NY,NX)=WTEABP(NB,NZ,NY,NX)*FWTEAR
      RWTGRB(NB,NZ,NY,NX)=WTGRB(NB,NZ,NY,NX)*FWTGR
      RWTGRBN(NB,NZ,NY,NX)=WTGRBN(NB,NZ,NY,NX)*FWTGR
      RWTGRBP(NB,NZ,NY,NX)=WTGRBP(NB,NZ,NY,NX)*FWTGR
      RCPOLNB(NB,NZ,NY,NX)=CPOLNB(NB,NZ,NY,NX)*FWPOLN
      RZPOLNB(NB,NZ,NY,NX)=ZPOLNB(NB,NZ,NY,NX)*FWPOLN
      RPPOLNB(NB,NZ,NY,NX)=PPOLNB(NB,NZ,NY,NX)*FWPOLN
      RWTNDB(NB,NZ,NY,NX)=WTNDB(NB,NZ,NY,NX)*FWTNDC
      RWTNDBN(NB,NZ,NY,NX)=WTNDBN(NB,NZ,NY,NX)*FWTNDC
      RWTNDBP(NB,NZ,NY,NX)=WTNDBP(NB,NZ,NY,NX)*FWTNDC
C
C     REMAINING UNCOMBUSTED SHOOT C POOLS AND ATTRIBUTES
C
      CPOOL(NB,NZ,NY,NX)=CPOOL(NB,NZ,NY,NX)-RCPOOL(NB,NZ,NY,NX)
      ZPOOL(NB,NZ,NY,NX)=ZPOOL(NB,NZ,NY,NX)-RZPOOL(NB,NZ,NY,NX)
      PPOOL(NB,NZ,NY,NX)=PPOOL(NB,NZ,NY,NX)-RPPOOL(NB,NZ,NY,NX)
      WTLFB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)-RWTLFB(NB,NZ,NY,NX)
      WTLFBN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)-RWTLFBN(NB,NZ,NY,NX)
      WTLFBP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)-RWTLFBP(NB,NZ,NY,NX)
      WTSHEB(NB,NZ,NY,NX)=WTSHEB(NB,NZ,NY,NX)-RWTSHEB(NB,NZ,NY,NX)
      WTSHBN(NB,NZ,NY,NX)=WTSHBN(NB,NZ,NY,NX)-RWTSHBN(NB,NZ,NY,NX)
      WTSHBP(NB,NZ,NY,NX)=WTSHBP(NB,NZ,NY,NX)-RWTSHBP(NB,NZ,NY,NX)
      WTSTKB(NB,NZ,NY,NX)=WTSTKB(NB,NZ,NY,NX)-RWTSTKB(NB,NZ,NY,NX)
      WTSTBN(NB,NZ,NY,NX)=WTSTBN(NB,NZ,NY,NX)-RWTSTBN(NB,NZ,NY,NX)
      WTSTBP(NB,NZ,NY,NX)=WTSTBP(NB,NZ,NY,NX)-RWTSTBP(NB,NZ,NY,NX)
      WTRSVB(NB,NZ,NY,NX)=WTRSVB(NB,NZ,NY,NX)-RWTRSVB(NB,NZ,NY,NX)
      WTRSBN(NB,NZ,NY,NX)=WTRSBN(NB,NZ,NY,NX)-RWTRSBN(NB,NZ,NY,NX)
      WTRSBP(NB,NZ,NY,NX)=WTRSBP(NB,NZ,NY,NX)-RWTRSBP(NB,NZ,NY,NX)
      WTHSKB(NB,NZ,NY,NX)=WTHSKB(NB,NZ,NY,NX)-RWTHSKB(NB,NZ,NY,NX)
      WTHSBN(NB,NZ,NY,NX)=WTHSBN(NB,NZ,NY,NX)-RWTHSBN(NB,NZ,NY,NX)
      WTHSBP(NB,NZ,NY,NX)=WTHSBP(NB,NZ,NY,NX)-RWTHSBP(NB,NZ,NY,NX)
      WTEARB(NB,NZ,NY,NX)=WTEARB(NB,NZ,NY,NX)-RWTEARB(NB,NZ,NY,NX)
      WTEABN(NB,NZ,NY,NX)=WTEABN(NB,NZ,NY,NX)-RWTEABN(NB,NZ,NY,NX)
      WTEABP(NB,NZ,NY,NX)=WTEABP(NB,NZ,NY,NX)-RWTEABP(NB,NZ,NY,NX)
      WTGRB(NB,NZ,NY,NX)=WTGRB(NB,NZ,NY,NX)-RWTGRB(NB,NZ,NY,NX)
      WTGRBN(NB,NZ,NY,NX)=WTGRBN(NB,NZ,NY,NX)-RWTGRBN(NB,NZ,NY,NX)
      WTGRBP(NB,NZ,NY,NX)=WTGRBP(NB,NZ,NY,NX)-RWTGRBP(NB,NZ,NY,NX)
      CPOLNB(NB,NZ,NY,NX)=CPOLNB(NB,NZ,NY,NX)-RCPOLNB(NB,NZ,NY,NX)
      ZPOLNB(NB,NZ,NY,NX)=ZPOLNB(NB,NZ,NY,NX)-RZPOLNB(NB,NZ,NY,NX)
      PPOLNB(NB,NZ,NY,NX)=PPOLNB(NB,NZ,NY,NX)-RPPOLNB(NB,NZ,NY,NX)
      WTNDB(NB,NZ,NY,NX)=WTNDB(NB,NZ,NY,NX)-RWTNDB(NB,NZ,NY,NX)
      WTNDBN(NB,NZ,NY,NX)=WTNDBN(NB,NZ,NY,NX)-RWTNDBN(NB,NZ,NY,NX)
      WTNDBP(NB,NZ,NY,NX)=WTNDBP(NB,NZ,NY,NX)-RWTNDBP(NB,NZ,NY,NX)
      WTSHTB(NB,NZ,NY,NX)=WTLFB(NB,NZ,NY,NX)
     2+WTSHEB(NB,NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
     3+WTHSKB(NB,NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)+WTGRB(NB,NZ,NY,NX)
     4+CPOOL(NB,NZ,NY,NX)+CPOOLK(NB,NZ,NY,NX)
      WTSHTN(NB,NZ,NY,NX)=WTLFBN(NB,NZ,NY,NX)
     2+WTSHBN(NB,NZ,NY,NX)+WTSTBN(NB,NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX)
     3+WTHSBN(NB,NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX)
     4+ZPOOL(NB,NZ,NY,NX)
      WTSHTP(NB,NZ,NY,NX)=WTLFBP(NB,NZ,NY,NX)
     2+WTSHBP(NB,NZ,NY,NX)+WTSTBP(NB,NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX)
     3+WTHSBP(NB,NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX)
     4+PPOOL(NB,NZ,NY,NX)
      ARLFB(NB,NZ,NY,NX)=ARLFB(NB,NZ,NY,NX)*(1.0-FWTLF)
C
C     TOTAL SHOOT COMBUSTION LOSSES
C
C     RCGCKZ=total canopy C combustion of PFT
C     VCOXF,VNOXF,VPOXF=total C,N,P combustion loss for output
C
      RCGCKB=RCPOOL(NB,NZ,NY,NX)
     2+RWTLFB(NB,NZ,NY,NX)+RWTSHEB(NB,NZ,NY,NX)+RWTSTKB(NB,NZ,NY,NX)
     3+RWTRSVB(NB,NZ,NY,NX)+RWTHSKB(NB,NZ,NY,NX)+RWTEARB(NB,NZ,NY,NX)
     4+RWTGRB(NB,NZ,NY,NX)+RCPOLNB(NB,NZ,NY,NX)+RWTNDB(NB,NZ,NY,NX)
      RCGNKB=RZPOOL(NB,NZ,NY,NX)
     2+RWTLFBN(NB,NZ,NY,NX)+RWTSHBN(NB,NZ,NY,NX)+RWTSTBN(NB,NZ,NY,NX)
     3+RWTRSBN(NB,NZ,NY,NX)+RWTHSBN(NB,NZ,NY,NX)+RWTEABN(NB,NZ,NY,NX)
     4+RWTGRBN(NB,NZ,NY,NX)+RZPOLNB(NB,NZ,NY,NX)+RWTNDBN(NB,NZ,NY,NX)
      RCGPKB=RPPOOL(NB,NZ,NY,NX)
     2+RWTLFBP(NB,NZ,NY,NX)+RWTSHBP(NB,NZ,NY,NX)+RWTSTBP(NB,NZ,NY,NX)
     3+RWTRSBP(NB,NZ,NY,NX)+RWTHSBP(NB,NZ,NY,NX)+RWTEABP(NB,NZ,NY,NX)
     4+RWTGRBP(NB,NZ,NY,NX)+RPPOLNB(NB,NZ,NY,NX)+RWTNDBP(NB,NZ,NY,NX)
      RCGCKZ(NZ,NY,NX)=RCGCKZ(NZ,NY,NX)+RCGCKB
      VCOXF(NZ,NY,NX)=VCOXF(NZ,NY,NX)-RCGCKB
      VNOXF(NZ,NY,NX)=VNOXF(NZ,NY,NX)-RCGNKB
      VPOXF(NZ,NY,NX)=VPOXF(NZ,NY,NX)-RCGPKB
C     WRITE(*,1910)'BURNSH',I,J,NFZ,NX,NY,NZ,NB,ICHKF
C    2,IDTHB(NB,NZ,NY,NX),RCGCKZ(NZ,NY,NX)
C    2,TKC(NZ,NY,NX),TFNCO,XNFH
C    2,VCOXF(NZ,NY,NX),VNOXF(NZ,NY,NX),VPOXF(NZ,NY,NX)
C    2,RCBCX,RCBC1,RCBC2,RCBC3,RCBC4,RCBC5
C    3,FWPOOL,FWTLF,FWTSHE,FWTSTK,FWTRSV,FWTHSK,FWTEAR,FWTGR
C    4,FWPOLN,FWTNDC,FWTSTG,WTLFB(NB,NZ,NY,NX),WTLF(NZ,NY,NX)
C    4,TWTLF(NY,NX),TKC(NZ,NY,NX),TKD(NZ,NY,NX)
1910  FORMAT(A8,9I4,40E12.4)
C
C     REMAINING UNCOMBUSTED SHOOT C POOLS AT EACH NODE
C
      DO 8745 K=0,25
      ARLF(K,NB,NZ,NY,NX)=ARLF(K,NB,NZ,NY,NX)*(1.0-FWTLF)
      HTSHE(K,NB,NZ,NY,NX)=HTSHE(K,NB,NZ,NY,NX)*(1.0-FWTSHE)
      WGLF(K,NB,NZ,NY,NX)=WGLF(K,NB,NZ,NY,NX)*(1.0-FWTLF)
      WSLF(K,NB,NZ,NY,NX)=WSLF(K,NB,NZ,NY,NX)*(1.0-FWTLF)
      WGLFN(K,NB,NZ,NY,NX)=WGLFN(K,NB,NZ,NY,NX)*(1.0-FWTLF)
      WGLFP(K,NB,NZ,NY,NX)=WGLFP(K,NB,NZ,NY,NX)*(1.0-FWTLF)
      WGSHE(K,NB,NZ,NY,NX)=WGSHE(K,NB,NZ,NY,NX)*(1.0-FWTSHE)
      WSSHE(K,NB,NZ,NY,NX)=WSSHE(K,NB,NZ,NY,NX)*(1.0-FWTSHE)
      WGSHN(K,NB,NZ,NY,NX)=WGSHN(K,NB,NZ,NY,NX)*(1.0-FWTSHE)
      WGSHP(K,NB,NZ,NY,NX)=WGSHP(K,NB,NZ,NY,NX)*(1.0-FWTSHE)
      WGNODE(K,NB,NZ,NY,NX)=WGNODE(K,NB,NZ,NY,NX)*(1.0-FWTSTK)
      WGNODN(K,NB,NZ,NY,NX)=WGNODN(K,NB,NZ,NY,NX)*(1.0-FWTSTK)
      WGNODP(K,NB,NZ,NY,NX)=WGNODP(K,NB,NZ,NY,NX)*(1.0-FWTSTK)
      DO 8745 L=1,JC
      ARLFL(L,K,NB,NZ,NY,NX)=ARLFL(L,K,NB,NZ,NY,NX)*(1.0-FWTLF)
      WGLFL(L,K,NB,NZ,NY,NX)=WGLFL(L,K,NB,NZ,NY,NX)*(1.0-FWTLF)
      WGLFLN(L,K,NB,NZ,NY,NX)=WGLFLN(L,K,NB,NZ,NY,NX)*(1.0-FWTLF)
      WGLFLP(L,K,NB,NZ,NY,NX)=WGLFLP(L,K,NB,NZ,NY,NX)*(1.0-FWTLF)
8745  CONTINUE
8750  CONTINUE
C
C     COMBUST STANDING DEAD
C
C     FW*=fraction of C combusted for each total C pool
C     R*=C,N,P combustion rate for each PFT C pool
C     PFT pool code:
C        WTSTDG,WTSTDN,WTSTDP=standing dead C,N,P mass
C
      DO 8740 M=1,4
      RWTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)*FWTSTG
      RWTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)*FWTSTG
      RWTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)*FWTSTG
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)-RWTSTDG(M,NZ,NY,NX)
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)-RWTSTDN(M,NZ,NY,NX)
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)-RWTSTDP(M,NZ,NY,NX)
C
C     TOTAL STANDING DEAD COMBUSTION LOSSES
C
C     RCGDKZ=total standing dead C combustion of PFT
C     VCOXF,VNOXF,VPOXF=total C,N,P combustion loss for output
C
      RCGDKZ(NZ,NY,NX)=RCGDKZ(NZ,NY,NX)+RWTSTDG(M,NZ,NY,NX)
      VCOXF(NZ,NY,NX)=VCOXF(NZ,NY,NX)-RWTSTDG(M,NZ,NY,NX)
      VNOXF(NZ,NY,NX)=VNOXF(NZ,NY,NX)-RWTSTDN(M,NZ,NY,NX)
      VPOXF(NZ,NY,NX)=VPOXF(NZ,NY,NX)-RWTSTDP(M,NZ,NY,NX)
C     WRITE(*,1912)'BURNST',I,J,NFZ,NX,NY,NZ,M,ICHKF
C    2,RCGDKZ(NZ,NY,NX),RWTSTDG(M,NZ,NY,NX),WTSTDG(M,NZ,NY,NX),FWTSTG
C    3,TFNDO,RCBC5,TWTSTG(NY,NX),WTSTG(NZ,NY,NX),RCGCK(NY,NX)
1912  FORMAT(A8,8I4,30E12.4)
8740  CONTINUE
C
C     COMBUST CHARCOAL
C
C     TFNCH=temperature function for combusting charcoal C
C        (700K = 1)
C     SPCMB=specific combustion rate of different combustion types
C        living plant material (0-4) at 600K, charcoal(5) at 700K
C     RCBC*=combustion rate for each combustion type
C     FWSTG=fraction of charcoal combusted
C     WTSTDG(5,WTSTDN(5,WTSTDP(5=charcoal C,N,P in standing dead
C
      TFNCH=AMIN1(TFNCX,EXP(20.620-120000/RTKD))
      RCBCX=TFNCH*AREA(3,NU(NY,NX),NY,NX)*XNFH
      RCBC5=SPCMB(5)*RCBCX
      IF(TWTSTG(NY,NX).GT.ZEROS(NY,NX))THEN
      FWTSTG=AMIN1(1.0,RCBC5/TWTSTG(NY,NX))
      ELSE
      FWTSTG=0.0
      ENDIF
      RWTSTDG(5,NZ,NY,NX)=WTSTDG(5,NZ,NY,NX)*FWTSTG
      RWTSTDN(5,NZ,NY,NX)=WTSTDN(5,NZ,NY,NX)*FWTSTG
      RWTSTDP(5,NZ,NY,NX)=WTSTDP(5,NZ,NY,NX)*FWTSTG
      WTSTDG(5,NZ,NY,NX)=WTSTDG(5,NZ,NY,NX)-RWTSTDG(5,NZ,NY,NX)
      WTSTDN(5,NZ,NY,NX)=WTSTDN(5,NZ,NY,NX)-RWTSTDN(5,NZ,NY,NX)
      WTSTDP(5,NZ,NY,NX)=WTSTDP(5,NZ,NY,NX)-RWTSTDP(5,NZ,NY,NX)
C
C     TOTAL CHARCOAL COMBUSTION LOSSES
C
C     RCGCKZ=total C combustion of each PFT
C     RCGCK=total canopy+standing dead C combustion of all PFT
C     VCOXF,VNOXF,VPOXF=total C,N,P combustion loss for output
C
      RCGDKZ(NZ,NY,NX)=RCGDKZ(NZ,NY,NX)+RWTSTDG(5,NZ,NY,NX)
      VCOXF(NZ,NY,NX)=VCOXF(NZ,NY,NX)-RWTSTDG(5,NZ,NY,NX)
      VNOXF(NZ,NY,NX)=VNOXF(NZ,NY,NX)-RWTSTDN(5,NZ,NY,NX)
      VPOXF(NZ,NY,NX)=VPOXF(NZ,NY,NX)-RWTSTDP(5,NZ,NY,NX)
      RCGCK(NY,NX)=RCGCK(NY,NX)+RCGCKZ(NZ,NY,NX)+RCGDKZ(NZ,NY,NX)
C
C     IF NO FIRE
C
      ELSE
      RCGCKZ(NZ,NY,NX)=0.0
      RCGDKZ(NZ,NY,NX)=0.0
      DO 8751 NB=1,NBR(NZ,NY,NX)
      RCPOOL(NB,NZ,NY,NX)=0.0
      RZPOOL(NB,NZ,NY,NX)=0.0
      RPPOOL(NB,NZ,NY,NX)=0.0
      RWTLFB(NB,NZ,NY,NX)=0.0
      RWTLFBN(NB,NZ,NY,NX)=0.0
      RWTLFBP(NB,NZ,NY,NX)=0.0
      RWTSHEB(NB,NZ,NY,NX)=0.0
      RWTSHBN(NB,NZ,NY,NX)=0.0
      RWTSHBP(NB,NZ,NY,NX)=0.0
      RWTSTKB(NB,NZ,NY,NX)=0.0
      RWTSTBN(NB,NZ,NY,NX)=0.0
      RWTSTBP(NB,NZ,NY,NX)=0.0
      RWTRSVB(NB,NZ,NY,NX)=0.0
      RWTRSBN(NB,NZ,NY,NX)=0.0
      RWTRSBP(NB,NZ,NY,NX)=0.0
      RWTHSKB(NB,NZ,NY,NX)=0.0
      RWTHSBN(NB,NZ,NY,NX)=0.0
      RWTHSBP(NB,NZ,NY,NX)=0.0
      RWTEARB(NB,NZ,NY,NX)=0.0
      RWTEABN(NB,NZ,NY,NX)=0.0
      RWTEABP(NB,NZ,NY,NX)=0.0
      RWTGRB(NB,NZ,NY,NX)=0.0
      RWTGRBN(NB,NZ,NY,NX)=0.0
      RWTGRBP(NB,NZ,NY,NX)=0.0
      RCPOLNB(NB,NZ,NY,NX)=0.0
      RZPOLNB(NB,NZ,NY,NX)=0.0
      RPPOLNB(NB,NZ,NY,NX)=0.0
      RWTNDB(NB,NZ,NY,NX)=0.0
      RWTNDBN(NB,NZ,NY,NX)=0.0
      RWTNDBP(NB,NZ,NY,NX)=0.0
8751  CONTINUE
      DO 8741 M=1,5
      RWTSTDG(M,NZ,NY,NX)=0.0
      RWTSTDN(M,NZ,NY,NX)=0.0
      RWTSTDP(M,NZ,NY,NX)=0.0
8741  CONTINUE
C     WRITE(*,1913)'BURNCH',I,J,NFZ,NX,NY,NZ,ICHKF
C    2,RCGDKZ(NZ,NY,NX),RWTSTDG(5,NZ,NY,NX),WTSTDG(5,NZ,NY,NX),FWTSTG
C    3,TFNCH,RCBC5,TWTSTG(NY,NX),WTSTG(NZ,NY,NX),RCGCK(NY,NX)
1913  FORMAT(A8,7I4,30E12.4)
      ENDIF
C
C     COMBUST ROOTS AND STORAGE
C
C     TKS=soil temperature
C     TCMBX=minimum temperature for combustion (K)
C     TFNCO,TFNDO=temperature function for combusting living C
C        (600K = 1)
C     SPCMB=specific combustion rate of different combustion types
C        living plant material (0-4) at 600K, charcoal(5) at 700K
C     RCBC*=combustion rate for each combustion type
C
      DO 8880 L=NU(NY,NX),NJ(NY,NX)
      IF(TKS(L,NY,NX).GT.TCMBX)THEN
      RTK=8.3143*TKS(L,NY,NX)
      TFNCO=AMIN1(TFNCX,EXP(12.028-60000/RTK))
      RCBCX=TFNCO*AREA(3,NU(NY,NX),NY,NX)*XNFH
      RCBC1=SPCMB(1)*RCBCX
      RCBC2=SPCMB(2)*RCBCX
      RCBC3=SPCMB(3)*RCBCX
      RCBC4=SPCMB(4)*RCBCX
      RCBC5=SPCMB(5)*RCBCX
C
C     ROOT AND STORAGE C POOL COMBUSTION FRACTIONS
C
C     FW*=fraction of C combusted for each total C pool
C     R*=C,N,P combustion rate for each PFT C pool
C     PFT pool code:
C        CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C        WTRTL=total active root C
C        WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C        WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C        CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C        WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C        WTRVC=storage C,N,P
C
      IF(TWPOLR(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FWPOLR=AMIN1(1.0,RCBC1/TWPOLR(L,NY,NX))
      ELSE
      FWPOLR=0.0
      ENDIF
      IF(TWTRTL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FWTRTL=AMIN1(1.0,RCBC4/TWTRTL(L,NY,NX))
      ELSE
      FWTRTL=0.0
      ENDIF
      IF(TWPODL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FWPODL=AMIN1(1.0,RCBC1/TWPODL(L,NY,NX))
      ELSE
      FWPODL=0.0
      ENDIF
      IF(TWTNDL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FWTNDL=AMIN1(1.0,RCBC2/TWTNDL(L,NY,NX))
      ELSE
      FWTNDL=0.0
      ENDIF
      IF(L.EQ.NU(NY,NX))THEN
C     IF(TWTRVC(NY,NX).GT.ZEROS(NY,NX))THEN
C     FWTRVC=AMIN1(1.0,RCBC1/TWTRVC(NY,NX))
C     ELSE
      FWTRVC=FWTRTL
C     ENDIF
C
C     STORAGE COMBUSTION RATES
C
C     R*=C,N,P combustion rate for each PFT C pool
C
      RWTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)*FWTRVC
      RWTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)*FWTRVC
      RWTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)*FWTRVC
C
C     REMAINING UNCOMBUSTED STORAGE C POOLS
C
      WTRVC(NZ,NY,NX)=WTRVC(NZ,NY,NX)-RWTRVC(NZ,NY,NX)
      WTRVN(NZ,NY,NX)=WTRVN(NZ,NY,NX)-RWTRVN(NZ,NY,NX)
      WTRVP(NZ,NY,NX)=WTRVP(NZ,NY,NX)-RWTRVP(NZ,NY,NX)
C
C     TOTAL STORAGE COMBUSTION LOSSES
C
C     RCGSKR=total root C combustion of PFT
C     VCOXF,VNOXF,VPOXF=total C,N,P combustion loss for output
C
      RCGSKR(L,NZ,NY,NX)=RCGSKR(L,NZ,NY,NX)+RWTRVC(NZ,NY,NX)
      VCOXF(NZ,NY,NX)=VCOXF(NZ,NY,NX)-RWTRVC(NZ,NY,NX)
      VNOXF(NZ,NY,NX)=VNOXF(NZ,NY,NX)-RWTRVN(NZ,NY,NX)
      VPOXF(NZ,NY,NX)=VPOXF(NZ,NY,NX)-RWTRVP(NZ,NY,NX)
C     WRITE(*,1909)'BURNRV',I,J,NFZ,NX,NY,NZ,L,NR,N,ICHKF
C    2,WTRVC(NZ,NY,NX),RWTRVC(NZ,NY,NX),FWTRVC,FWTRTL
C    3,RCBC4,TWTRTL(L,NY,NX)
      ENDIF
C
C     ROOT COMBUSTION RATES
C
C     R*=C,N,P combustion rate for each PFT C pool
C
      DO 8885 N=1,MY(NZ,NY,NX)
      RCPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)*FWPOLR
      RZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)*FWPOLR
      RPPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)*FWPOLR
      CPOOLR(N,L,NZ,NY,NX)=CPOOLR(N,L,NZ,NY,NX)-RCPOOLR(N,L,NZ,NY,NX)
      ZPOOLR(N,L,NZ,NY,NX)=ZPOOLR(N,L,NZ,NY,NX)-RZPOOLR(N,L,NZ,NY,NX)
      PPOOLR(N,L,NZ,NY,NX)=PPOOLR(N,L,NZ,NY,NX)-RPPOOLR(N,L,NZ,NY,NX)
C
C     TOTAL ROOT COMBUSTION LOSSES
C
C     RCGSKR=total root C combustion of PFT
C     VCOXF,VNOXF,VPOXF=total C,N,P combustion loss for output
C
      RCGSKR(L,NZ,NY,NX)=RCGSKR(L,NZ,NY,NX)+RCPOOLR(N,L,NZ,NY,NX)
      VCOXF(NZ,NY,NX)=VCOXF(NZ,NY,NX)-RCPOOLR(N,L,NZ,NY,NX)
      VNOXF(NZ,NY,NX)=VNOXF(NZ,NY,NX)-RZPOOLR(N,L,NZ,NY,NX)
      VPOXF(NZ,NY,NX)=VPOXF(NZ,NY,NX)-RPPOOLR(N,L,NZ,NY,NX)
C
C     REMAINING UNCOMBUSTED ROOT C POOLS AND ATTRIBUTES
C
      WTRTL(N,L,NZ,NY,NX)=WTRTL(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      WTRTD(N,L,NZ,NY,NX)=WTRTD(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      WSRTL(N,L,NZ,NY,NX)=WSRTL(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      RTN1(N,L,NZ,NY,NX)=RTN1(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      RTNL(N,L,NZ,NY,NX)=RTNL(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      RTLGP(N,L,NZ,NY,NX)=RTLGP(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      RTDNP(N,L,NZ,NY,NX)=RTDNP(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      RTVLP(N,L,NZ,NY,NX)=RTVLP(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      RTVLW(N,L,NZ,NY,NX)=RTVLW(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      RTARP(N,L,NZ,NY,NX)=RTARP(N,L,NZ,NY,NX)*(1.0-FWTRTL)
      DO 8785 NR=1,NRT(NZ,NY,NX)
      RWTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)*FWTRTL
      RWTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)*FWTRTL
      RWTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)*FWTRTL
      RWTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)*FWTRTL
      RWTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)*FWTRTL
      RWTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)*FWTRTL
      WTRT1(N,L,NR,NZ,NY,NX)=WTRT1(N,L,NR,NZ,NY,NX)
     2-RWTRT1(N,L,NR,NZ,NY,NX)
      WTRT2(N,L,NR,NZ,NY,NX)=WTRT2(N,L,NR,NZ,NY,NX)
     2-RWTRT2(N,L,NR,NZ,NY,NX)
      WTRT1N(N,L,NR,NZ,NY,NX)=WTRT1N(N,L,NR,NZ,NY,NX)
     2-RWTRT1N(N,L,NR,NZ,NY,NX)
      WTRT2N(N,L,NR,NZ,NY,NX)=WTRT2N(N,L,NR,NZ,NY,NX)
     2-RWTRT2N(N,L,NR,NZ,NY,NX)
      WTRT1P(N,L,NR,NZ,NY,NX)=WTRT1P(N,L,NR,NZ,NY,NX)
     2-RWTRT1P(N,L,NR,NZ,NY,NX)
      WTRT2P(N,L,NR,NZ,NY,NX)=WTRT2P(N,L,NR,NZ,NY,NX)
     2-RWTRT2P(N,L,NR,NZ,NY,NX)
C
C     TOTAL ROOT COMBUSTION LOSSES
C
C     RCGSKR=total root C combustion of PFT
C     VCOXF,VNOXF,VPOXF=total C,N,P combustion loss for output
C
      RCGSKR(L,NZ,NY,NX)=RCGSKR(L,NZ,NY,NX)+RWTRT1(N,L,NR,NZ,NY,NX)
     2+RWTRT2(N,L,NR,NZ,NY,NX)
      VCOXF(NZ,NY,NX)=VCOXF(NZ,NY,NX)-(RWTRT1(N,L,NR,NZ,NY,NX)
     2+RWTRT2(N,L,NR,NZ,NY,NX))
      VNOXF(NZ,NY,NX)=VNOXF(NZ,NY,NX)-(RWTRT1N(N,L,NR,NZ,NY,NX)
     2+RWTRT2N(N,L,NR,NZ,NY,NX))
      VPOXF(NZ,NY,NX)=VPOXF(NZ,NY,NX)-(RWTRT1P(N,L,NR,NZ,NY,NX)
     2+RWTRT2P(N,L,NR,NZ,NY,NX))
      RTWT1(N,NR,NZ,NY,NX)=RTWT1(N,NR,NZ,NY,NX)*(1.0-FWTRTL)
      RTWT1N(N,NR,NZ,NY,NX)=RTWT1N(N,NR,NZ,NY,NX)*(1.0-FWTRTL)
      RTWT1P(N,NR,NZ,NY,NX)=RTWT1P(N,NR,NZ,NY,NX)*(1.0-FWTRTL)
      RTLG1(N,L,NR,NZ,NY,NX)=RTLG1(N,L,NR,NZ,NY,NX)*(1.0-FWTRTL)
      RTLG2(N,L,NR,NZ,NY,NX)=RTLG2(N,L,NR,NZ,NY,NX)*(1.0-FWTRTL)
      RTN2(N,L,NR,NZ,NY,NX)=RTN2(N,L,NR,NZ,NY,NX)*(1.0-FWTRTL)
C     WRITE(*,1909)'BURNRT',I,J,NFZ,NX,NY,NZ,L,NR,N,ICHKF
C    2,RCGSKR(L,NZ,NY,NX),RWTRT1(N,L,NR,NZ,NY,NX)
C    3,RWTRT2(N,L,NR,NZ,NY,NX)
1909  FORMAT(A8,10I4,30E12.4)
8785  CONTINUE
8885  CONTINUE
C
C     ROOT NODULE COMBUSTION RATES
C
C     R*=C,N,P combustion rate for each PFT C pool
C
      RCPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)*FWPODL
      RZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)*FWPODL
      RPPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)*FWPODL
      RWTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)*FWTNDL
      RWTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)*FWTNDL
      RWTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)*FWTNDL
C
C     REMAINING UNCOMBUSTED ROOT NODULE C,N,P
C
      CPOOLN(L,NZ,NY,NX)=CPOOLN(L,NZ,NY,NX)-RCPOOLN(L,NZ,NY,NX)
      ZPOOLN(L,NZ,NY,NX)=ZPOOLN(L,NZ,NY,NX)-RZPOOLN(L,NZ,NY,NX)
      PPOOLN(L,NZ,NY,NX)=PPOOLN(L,NZ,NY,NX)-RPPOOLN(L,NZ,NY,NX)
      WTNDL(L,NZ,NY,NX)=WTNDL(L,NZ,NY,NX)-RWTNDL(L,NZ,NY,NX)
      WTNDLN(L,NZ,NY,NX)=WTNDLN(L,NZ,NY,NX)-RWTNDLN(L,NZ,NY,NX)
      WTNDLP(L,NZ,NY,NX)=WTNDLP(L,NZ,NY,NX)-RWTNDLP(L,NZ,NY,NX)
C
C     TOTAL ROOT NODULE COMBUSTION LOSSES
C
C     RCGSKR=total root C combustion of PFT
C     VCOXF,VNOXF,VPOXF=total C,N,P combustion loss for output
C     RCGSK=total root+soil C combustion
C
      RCGSKR(L,NZ,NY,NX)=RCGSKR(L,NZ,NY,NX)+RCPOOLN(L,NZ,NY,NX)
     2+RWTNDL(L,NZ,NY,NX)
      VCOXF(NZ,NY,NX)=VCOXF(NZ,NY,NX)-(RCPOOLN(L,NZ,NY,NX)
     2+RWTNDL(L,NZ,NY,NX))
      VNOXF(NZ,NY,NX)=VNOXF(NZ,NY,NX)-(RZPOOLN(L,NZ,NY,NX)
     2+RWTNDLN(L,NZ,NY,NX))
      VPOXF(NZ,NY,NX)=VPOXF(NZ,NY,NX)-(RPPOOLN(L,NZ,NY,NX)
     2+RWTNDLP(L,NZ,NY,NX))
      RCGSK(L,NY,NX)=RCGSK(L,NY,NX)+RCGSKR(L,NZ,NY,NX)
C     IF(L.EQ.NU(NY,NX)+1)THEN
C     WRITE(*,1911)'BURNRT',I,J,NFZ,NX,NY,NZ,L,ICHKF,NG(NZ,NY,NX)
C    2,RCGSKR(L,NZ,NY,NX),RCGSK(L,NY,NX)
C    3,RCPOOLN(L,NZ,NY,NX),RWTNDL(L,NZ,NY,NX)
C    3,RWTRVC(NZ,NY,NX),(RCPOOLR(N,L,NZ,NY,NX),N=1,2)
C    4,WTRVC(NZ,NY,NX),TFNCO,TKS(L,NY,NX)
1911  FORMAT(A8,9I4,30E12.4)
C     ENDIF
C
C     IF NO FIRE
C
      ELSE
      RCGSKR(L,NZ,NY,NX)=0.0
      IF(L.EQ.NU(NY,NX)+1)THEN
      RWTRVC(NZ,NY,NX)=0.0
      RWTRVN(NZ,NY,NX)=0.0
      RWTRVP(NZ,NY,NX)=0.0
      ENDIF
      DO 8886 N=1,MY(NZ,NY,NX)
      RCPOOLR(N,L,NZ,NY,NX)=0.0
      RZPOOLR(N,L,NZ,NY,NX)=0.0
      RPPOOLR(N,L,NZ,NY,NX)=0.0
      DO 8786 NR=1,NRT(NZ,NY,NX)
      RWTRT1(N,L,NR,NZ,NY,NX)=0.0
      RWTRT2(N,L,NR,NZ,NY,NX)=0.0
      RWTRT1N(N,L,NR,NZ,NY,NX)=0.0
      RWTRT2N(N,L,NR,NZ,NY,NX)=0.0
      RWTRT1P(N,L,NR,NZ,NY,NX)=0.0
      RWTRT2P(N,L,NR,NZ,NY,NX)=0.0
8786  CONTINUE
8886  CONTINUE
      RCPOOLN(L,NZ,NY,NX)=0.0
      RZPOOLN(L,NZ,NY,NX)=0.0
      RPPOOLN(L,NZ,NY,NX)=0.0
      RWTNDL(L,NZ,NY,NX)=0.0
      RWTNDLN(L,NZ,NY,NX)=0.0
      RWTNDLP(L,NZ,NY,NX)=0.0
      ENDIF
8880  CONTINUE
9970  CONTINUE
      ENDIF
C
C     END FIRE
C
C
C     ACTIVATE DORMANT SEEDS
C
C     IFLGI=PFT initialization flag:0=no,1=yes
C     IFLGE=flag for enabling leafout:0=enable,1=disable
C     VRNS,VRNL=leafout hours,hours required for leafout
C     IDAY0,IYR0=day,year of planting
C     IYRC=current year
C     SDPTHI=seeding depth
C
      DO 9975 NZ=1,NP0(NY,NX)
      DO 205 NB=1,NBR(NZ,NY,NX)
      IF(NFZ.EQ.1)THEN
      IF(IFLGI(NZ,NY,NX).EQ.1)THEN
      IF(IFLGE(NB,NZ,NY,NX).EQ.0
     2.AND.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX))THEN
      IDAY0(NZ,NY,NX)=I
      IYR0(NZ,NY,NX)=IYRC
      SDPTHI(NZ,NY,NX)=0.005+CDPTHZ(0,NY,NX)
      IFLGI(NZ,NY,NX)=0
      ENDIF
      ENDIF
      ENDIF
205   CONTINUE
C
C     LITTERFALL FROM STANDING DEAD
C
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IGTYP=growth type:0=bryophyte,1=grass,herb,2=woody shrub,tree
C     XFRC,XFRN,XFRP=litterfall from standing dead
C     TFN3=temperature function for canopy growth
C     WTSTG,WTSTDN,WTSTDP=standing dead C,N,P mass
C     CSNC,ZSNC,PSNC=C,N,P litterfall
C     FWOOD,FWOODN,FWOODP=C,N,P woody fraction in root:0=woody,
C        1=non-woody
C
      DO 6235 M=1,5
      TFNS=SQRT(TFN3(NZ,NY,NX))
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      XFRC=1.5814E-04*TFNS*WTSTDG(M,NZ,NY,NX)*XNFH
      XFRN=1.5814E-04*TFNS*WTSTDN(M,NZ,NY,NX)*XNFH
      XFRP=1.5814E-04*TFNS*WTSTDP(M,NZ,NY,NX)*XNFH
      ELSE
      XFRC=1.5814E-05*TFNS*WTSTDG(M,NZ,NY,NX)*XNFH
      XFRN=1.5814E-05*TFNS*WTSTDN(M,NZ,NY,NX)*XNFH
      XFRP=1.5814E-05*TFNS*WTSTDP(M,NZ,NY,NX)*XNFH
      ENDIF
      CSNC(M,0,0,NZ,NY,NX)=CSNC(M,0,0,NZ,NY,NX)+XFRC*FWOOD(0,NZ)
      ZSNC(M,0,0,NZ,NY,NX)=ZSNC(M,0,0,NZ,NY,NX)+XFRN*FWOODN(0,NZ)
      PSNC(M,0,0,NZ,NY,NX)=PSNC(M,0,0,NZ,NY,NX)+XFRP*FWOODP(0,NZ)
      CSNC(M,1,0,NZ,NY,NX)=CSNC(M,1,0,NZ,NY,NX)+XFRC*FWOOD(1,NZ)
      ZSNC(M,1,0,NZ,NY,NX)=ZSNC(M,1,0,NZ,NY,NX)+XFRN*FWOODN(1,NZ)
      PSNC(M,1,0,NZ,NY,NX)=PSNC(M,1,0,NZ,NY,NX)+XFRP*FWOODP(1,NZ)
      WTSTDG(M,NZ,NY,NX)=WTSTDG(M,NZ,NY,NX)-XFRC
      WTSTDN(M,NZ,NY,NX)=WTSTDN(M,NZ,NY,NX)-XFRN
      WTSTDP(M,NZ,NY,NX)=WTSTDP(M,NZ,NY,NX)-XFRP
C     IF(NZ.EQ.2)THEN
C     WRITE (*,1919)'STGLIT',I,J,NFZ,NZ,M,XFRC,XFRN,XFRP
C    2,FWOOD(1,NZ),FWOOD(0,NZ),TFN3(NZ,NY,NX)
C    2,WTSTDG(M,NZ,NY,NX),WTSTDN(M,NZ,NY,NX),WTSTDP(M,NZ,NY,NX)
C    3,CSNC(M,1,0,NZ,NY,NX),ZSNC(M,1,0,NZ,NY,NX),PSNC(M,1,0,NZ,NY,NX)
C    3,CSNC(M,0,0,NZ,NY,NX),ZSNC(M,0,0,NZ,NY,NX),PSNC(M,0,0,NZ,NY,NX)
1919  FORMAT(A8,5I4,20E12.4)
C     ENDIF
6235  CONTINUE
C
C     ACCUMULATE TOTAL SURFACE, SUBSURFACE LITTERFALL
C
C     TCSN0,TZSN0,TPSN0=cumulative above-ground C,N,P litterfall
C     TCSNC,TZSNC,TPSNC=cumulative C,N,P litterfall
C     HCSNC,HZSNC,HPSNC=hourly C,N,P litterfall
C     CSNC,ZSNC,PSNC=C,N,P litterfall
C
      DO 6430 K=0,1
      DO 6430 M=1,5
      TCSN0(NZ,NY,NX)=TCSN0(NZ,NY,NX)+CSNC(M,K,0,NZ,NY,NX)
      TZSN0(NZ,NY,NX)=TZSN0(NZ,NY,NX)+ZSNC(M,K,0,NZ,NY,NX)
      TPSN0(NZ,NY,NX)=TPSN0(NZ,NY,NX)+PSNC(M,K,0,NZ,NY,NX)
C     WRITE(*,4814)'CSNC0',I,J,NFZ,NX,NY,NZ,NB,K,M
C    2,TCSN0(NZ,NY,NX),CSNC(M,K,0,NZ,NY,NX)
4814  FORMAT(A8,9I4,12E12.4)
      DO 8955 L=0,NJ(NY,NX)
      HCSNC(NZ,NY,NX)=HCSNC(NZ,NY,NX)+CSNC(M,K,L,NZ,NY,NX)
      HZSNC(NZ,NY,NX)=HZSNC(NZ,NY,NX)+ZSNC(M,K,L,NZ,NY,NX)
      HPSNC(NZ,NY,NX)=HPSNC(NZ,NY,NX)+PSNC(M,K,L,NZ,NY,NX)
      TCSNC(NZ,NY,NX)=TCSNC(NZ,NY,NX)+CSNC(M,K,L,NZ,NY,NX)
      TZSNC(NZ,NY,NX)=TZSNC(NZ,NY,NX)+ZSNC(M,K,L,NZ,NY,NX)
      TPSNC(NZ,NY,NX)=TPSNC(NZ,NY,NX)+PSNC(M,K,L,NZ,NY,NX)
8955  CONTINUE
6430  CONTINUE
C
C     CALCULATE PFT STATE VARIABLES
C
      CPOOLP(NZ,NY,NX)=0.0
      ZPOOLP(NZ,NY,NX)=0.0
      PPOOLP(NZ,NY,NX)=0.0
      WTSHTX=WTSHT(NZ,NY,NX)
      WTSHT(NZ,NY,NX)=0.0
      WTSHN(NZ,NY,NX)=0.0
      WTSHP(NZ,NY,NX)=0.0
      WTLF(NZ,NY,NX)=0.0
      WTSHE(NZ,NY,NX)=0.0
      WTSTK(NZ,NY,NX)=0.0
      WVSTK(NZ,NY,NX)=0.0
      WTRSV(NZ,NY,NX)=0.0
      WTHSK(NZ,NY,NX)=0.0
      WTEAR(NZ,NY,NX)=0.0
      WTGR(NZ,NY,NX)=0.0
      WTLS(NZ,NY,NX)=0.0
      WTRT(NZ,NY,NX)=0.0
      WTRTS(NZ,NY,NX)=0.0
      WTRTN(NZ,NY,NX)=0.0
      WTRTP(NZ,NY,NX)=0.0
      WTLFN(NZ,NY,NX)=0.0
      WTSHEN(NZ,NY,NX)=0.0
      WTSTKN(NZ,NY,NX)=0.0
      WTRSVN(NZ,NY,NX)=0.0
      WTHSKN(NZ,NY,NX)=0.0
      WTEARN(NZ,NY,NX)=0.0
      WTGRNN(NZ,NY,NX)=0.0
      WTLFP(NZ,NY,NX)=0.0
      WTSHEP(NZ,NY,NX)=0.0
      WTSTKP(NZ,NY,NX)=0.0
      WTRSVP(NZ,NY,NX)=0.0
      WTHSKP(NZ,NY,NX)=0.0
      WTEARP(NZ,NY,NX)=0.0
      WTGRNP(NZ,NY,NX)=0.0
      GRNO(NZ,NY,NX)=0.0
      ARLFP(NZ,NY,NX)=0.0
      ARSTP(NZ,NY,NX)=0.0
      DO 8940 L=1,JC
      ARSTV(L,NZ,NY,NX)=0.0
8940  CONTINUE
C
C     ACCUMULATE PFT STATE VARIABLES FROM BRANCH STATE VARIABLES
C
C     CPOOL,ZPOOL,PPOOL=non-structural C,N,P mass in branch
C     WTSHTB,WTSHTN,WTSHTP=branch total C,N,P mass
C     WTRSVB,WTRSBN,WTRSBP=stalk reserve C,N,P mass
C     WTLFB,WTLFBN,WTLFBP=branch leaf C,N,P mass
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C     WTSHEB,WTSHBN,WTSHBP=branch petiole C,N,P mass
C     WTHSKB,WTEARB,WTGRB=branch husk,ear,grain C mass
C     WTHSBN,WTEABN,WTGRBN=branch husk,ear,grain N mass
C     WTHSBP,WTEABP,WTGRBP=branch husk,ear,grain P mass
C     WTRVC,WTRVN,WTRVP=storage C,N,P
C     ARLFB=branch leaf area
C     ARSTK=total branch stalk surface area in each layer
C     GRNOB=seed set number
C
      DO 8950 NB=1,NBR(NZ,NY,NX)
      CPOOLP(NZ,NY,NX)=CPOOLP(NZ,NY,NX)+CPOOL(NB,NZ,NY,NX)
      ZPOOLP(NZ,NY,NX)=ZPOOLP(NZ,NY,NX)+ZPOOL(NB,NZ,NY,NX)
      PPOOLP(NZ,NY,NX)=PPOOLP(NZ,NY,NX)+PPOOL(NB,NZ,NY,NX)
      WTSHT(NZ,NY,NX)=WTSHT(NZ,NY,NX)+WTSHTB(NB,NZ,NY,NX)
      WTLF(NZ,NY,NX)=WTLF(NZ,NY,NX)+WTLFB(NB,NZ,NY,NX)
      WTSHE(NZ,NY,NX)=WTSHE(NZ,NY,NX)+WTSHEB(NB,NZ,NY,NX)
      WTSTK(NZ,NY,NX)=WTSTK(NZ,NY,NX)+WTSTKB(NB,NZ,NY,NX)
      WVSTK(NZ,NY,NX)=WVSTK(NZ,NY,NX)+WVSTKB(NB,NZ,NY,NX)
      WTRSV(NZ,NY,NX)=WTRSV(NZ,NY,NX)+WTRSVB(NB,NZ,NY,NX)
      WTHSK(NZ,NY,NX)=WTHSK(NZ,NY,NX)+WTHSKB(NB,NZ,NY,NX)
      WTEAR(NZ,NY,NX)=WTEAR(NZ,NY,NX)+WTEARB(NB,NZ,NY,NX)
      WTGR(NZ,NY,NX)=WTGR(NZ,NY,NX)+WTGRB(NB,NZ,NY,NX)
      WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(NB,NZ,NY,NX)
      WTSHN(NZ,NY,NX)=WTSHN(NZ,NY,NX)+WTSHTN(NB,NZ,NY,NX)
      WTLFN(NZ,NY,NX)=WTLFN(NZ,NY,NX)+WTLFBN(NB,NZ,NY,NX)
      WTSHEN(NZ,NY,NX)=WTSHEN(NZ,NY,NX)+WTSHBN(NB,NZ,NY,NX)
      WTSTKN(NZ,NY,NX)=WTSTKN(NZ,NY,NX)+WTSTBN(NB,NZ,NY,NX)
      WTRSVN(NZ,NY,NX)=WTRSVN(NZ,NY,NX)+WTRSBN(NB,NZ,NY,NX)
      WTHSKN(NZ,NY,NX)=WTHSKN(NZ,NY,NX)+WTHSBN(NB,NZ,NY,NX)
      WTEARN(NZ,NY,NX)=WTEARN(NZ,NY,NX)+WTEABN(NB,NZ,NY,NX)
      WTGRNN(NZ,NY,NX)=WTGRNN(NZ,NY,NX)+WTGRBN(NB,NZ,NY,NX)
      WTSHP(NZ,NY,NX)=WTSHP(NZ,NY,NX)+WTSHTP(NB,NZ,NY,NX)
      WTLFP(NZ,NY,NX)=WTLFP(NZ,NY,NX)+WTLFBP(NB,NZ,NY,NX)
      WTSHEP(NZ,NY,NX)=WTSHEP(NZ,NY,NX)+WTSHBP(NB,NZ,NY,NX)
      WTSTKP(NZ,NY,NX)=WTSTKP(NZ,NY,NX)+WTSTBP(NB,NZ,NY,NX)
      WTRSVP(NZ,NY,NX)=WTRSVP(NZ,NY,NX)+WTRSBP(NB,NZ,NY,NX)
      WTHSKP(NZ,NY,NX)=WTHSKP(NZ,NY,NX)+WTHSBP(NB,NZ,NY,NX)
      WTEARP(NZ,NY,NX)=WTEARP(NZ,NY,NX)+WTEABP(NB,NZ,NY,NX)
      WTGRNP(NZ,NY,NX)=WTGRNP(NZ,NY,NX)+WTGRBP(NB,NZ,NY,NX)
      ARLFP(NZ,NY,NX)=ARLFP(NZ,NY,NX)+ARLFB(NB,NZ,NY,NX)
      GRNO(NZ,NY,NX)=GRNO(NZ,NY,NX)+GRNOB(NB,NZ,NY,NX)
      DO 8945 L=1,JC
      ARSTP(NZ,NY,NX)=ARSTP(NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
      ARSTV(L,NZ,NY,NX)=ARSTV(L,NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
8945  CONTINUE
8950  CONTINUE
      DWTSHT(NZ,NY,NX)=WTSHT(NZ,NY,NX)-WTSHTX
C
C     ACCUMULATE ROOT STATE VARIABLES FROM ROOT LAYER STATE VARIABLES
C
C     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C     WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C     WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C
      DO 8925 N=1,MY(NZ,NY,NX)
      DO 8930 L=NU(NY,NX),NJ(NY,NX)
      WTRT(NZ,NY,NX)=WTRT(NZ,NY,NX)+CPOOLR(N,L,NZ,NY,NX)
      WTRTN(NZ,NY,NX)=WTRTN(NZ,NY,NX)+ZPOOLR(N,L,NZ,NY,NX)
      WTRTP(NZ,NY,NX)=WTRTP(NZ,NY,NX)+PPOOLR(N,L,NZ,NY,NX)
      DO 8935 NR=1,NRT(NZ,NY,NX)
      WTRT(NZ,NY,NX)=WTRT(NZ,NY,NX)+WTRT1(N,L,NR,NZ,NY,NX)
     2+WTRT2(N,L,NR,NZ,NY,NX)
      WTRTS(NZ,NY,NX)=WTRTS(NZ,NY,NX)+WTRT1(N,L,NR,NZ,NY,NX)
     2+WTRT2(N,L,NR,NZ,NY,NX)
      WTRTN(NZ,NY,NX)=WTRTN(NZ,NY,NX)+WTRT1N(N,L,NR,NZ,NY,NX)
     2+WTRT2N(N,L,NR,NZ,NY,NX)
      WTRTP(NZ,NY,NX)=WTRTP(NZ,NY,NX)+WTRT1P(N,L,NR,NZ,NY,NX)
     2+WTRT2P(N,L,NR,NZ,NY,NX)
8935  CONTINUE
8930  CONTINUE
8925  CONTINUE
C
C     ACCUMULATE NODULE STATE VATIABLES FROM NODULE LAYER VARIABLES
C
C     INTYP=N2 fixation: 1,2,3=rapid to slow root symbiosis
C     CPOLNB,ZPOLNB,PPOLNB=nonstructural C,N,P in bacteria
C     WTNDB,WTNDBN,WTNDBP=bacterial C,N,P mass
C
      IF(INTYP(NZ,NY,NX).NE.0)THEN
      WTND(NZ,NY,NX)=0.0
      WTNDN(NZ,NY,NX)=0.0
      WTNDP(NZ,NY,NX)=0.0
      IF(INTYP(NZ,NY,NX).GE.4)THEN
      DO 7950 NB=1,NBR(NZ,NY,NX)
      CPOLNP(NZ,NY,NX)=CPOLNP(NZ,NY,NX)+CPOLNB(NB,NZ,NY,NX)
      ZPOLNP(NZ,NY,NX)=ZPOLNP(NZ,NY,NX)+ZPOLNB(NB,NZ,NY,NX)
      PPOLNP(NZ,NY,NX)=PPOLNP(NZ,NY,NX)+PPOLNB(NB,NZ,NY,NX)
      WTND(NZ,NY,NX)=WTND(NZ,NY,NX)+WTNDB(NB,NZ,NY,NX)
     2+CPOLNB(NB,NZ,NY,NX)
      WTNDN(NZ,NY,NX)=WTNDN(NZ,NY,NX)+WTNDBN(NB,NZ,NY,NX)
     2+ZPOLNB(NB,NZ,NY,NX)
      WTNDP(NZ,NY,NX)=WTNDP(NZ,NY,NX)+WTNDBP(NB,NZ,NY,NX)
     2+PPOLNB(NB,NZ,NY,NX)
7950  CONTINUE
      ELSEIF(INTYP(NZ,NY,NX).GE.1.AND.INTYP(NZ,NY,NX).LE.3)THEN
      DO 8920 L=NU(NY,NX),NI(NZ,NY,NX)
      WTND(NZ,NY,NX)=WTND(NZ,NY,NX)+WTNDL(L,NZ,NY,NX)
     2+CPOOLN(L,NZ,NY,NX)
      WTNDN(NZ,NY,NX)=WTNDN(NZ,NY,NX)+WTNDLN(L,NZ,NY,NX)
     2+ZPOOLN(L,NZ,NY,NX)
      WTNDP(NZ,NY,NX)=WTNDP(NZ,NY,NX)+WTNDLP(L,NZ,NY,NX)
     2+PPOOLN(L,NZ,NY,NX)
8920  CONTINUE
      ENDIF
      ENDIF
C
C     PLANT C BALANCE = TOTAL C STATE VARIABLES + TOTAL
C     AUTOTROPHIC RESPIRATION + TOTAL LITTERFALL - TOTAL EXUDATION
C     - TOTAL CO2 FIXATION
C
C     IFLGC=PFT flag:0=not active,1=active
C     BALC=PFT C balance
C     WTSHT,WTRT,WTND,WTRVC,WTSTG=PFT
C        shoot,root,bacteria,storage,standing dead C
C     ZNPP=cumulative PFT NPP
C     TCSNC=cumulative PFT C litterfall
C     TCUPTK=cumulative PFT root-soil C exchange
C     RSETC=cumulative C balance from previous year
C     THVSTC=cumulative PFT C removed from ecosystem from
C        previous year
C     HVSTC=total PFT C removed from ecosystem in current year
C     VCOXF=CO2 emission from disturbance
C
      IF(IFLGC(NZ,NY,NX).EQ.1)THEN
      BALC(NZ,NY,NX)=WTSHT(NZ,NY,NX)+WTRT(NZ,NY,NX)+WTND(NZ,NY,NX)
     2+WTRVC(NZ,NY,NX)-ZNPP(NZ,NY,NX)+TCSNC(NZ,NY,NX)-TCUPTK(NZ,NY,NX)
     3-RSETC(NZ,NY,NX)+WTSTG(NZ,NY,NX)+THVSTC(NZ,NY,NX)
     4+HVSTC(NZ,NY,NX)-VCOXF(NZ,NY,NX)
C     IF(NZ.EQ.1.AND.ICHKF.EQ.1)THEN
C     WRITE(*,1111)'BALC',I,J,NFZ,NX,NY,NZ
C    2,BALC(NZ,NY,NX),WTSHT(NZ,NY,NX),WTRT(NZ,NY,NX)
C    2,WTND(NZ,NY,NX),WTRVC(NZ,NY,NX),ZNPP(NZ,NY,NX)
C    3,TCSNC(NZ,NY,NX),TCUPTK(NZ,NY,NX),RSETC(NZ,NY,NX)
C    3,WTSTG(NZ,NY,NX),UPOMC(NZ,NY,NX),THVSTC(NZ,NY,NX)
C    2,HVSTC(NZ,NY,NX),CARBN(NZ,NY,NX),TCO2T(NZ,NY,NX)
C    3,TCO2A(NZ,NY,NX),RCGCKZ(NZ,NY,NX),VCOXF(NZ,NY,NX)
C    4,(WTSTDG(M,NZ,NY,NX),M=1,5)
C    3,(((CSNC(M,K,L,NZ,NY,NX),M=1,4),K=0,1),L=0,0)
C    3,WTLF(NZ,NY,NX),WTSHE(NZ,NY,NX),WTSTK(NZ,NY,NX),WTRSV(NZ,NY,NX)
C    3,WTHSK(NZ,NY,NX),WTEAR(NZ,NY,NX),WTGR(NZ,NY,NX),CPOOLP(NZ,NY,NX)
C    4,(WTSTKB(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
C    4,(CPOOL(NB,NZ,NY,NX),NB=1,NBR(NZ,NY,NX))
C    4,((CPOOLR(N,L,NZ,NY,NX),L=NU(NY,NX),NJ(NY,NX)),N=1,2)
C    5,(((WTRT1(N,L,NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
C    2,L=NU(NY,NX),NJ(NY,NX)),N=1,2)
C    5,(((WTRT2(N,L,NR,NZ,NY,NX),NR=1,NRT(NZ,NY,NX))
C    2,L=NU(NY,NX),NJ(NY,NX)),N=1,2)
1111  FORMAT(A8,6I4,360F14.6)
C     ENDIF
C
C     PLANT N BALANCE = TOTAL N STATE VARIABLES + TOTAL N LITTERFALL
C     - TOTAL N UPTAKE FROM SOIL - TOTAL N ABSORPTION FROM ATMOSPHERE
C
C     BALN=PFT N balance
C     WTSHN,WTRTN,WTNDN,WTRVN,WTSTGN=PFT
C        shoot,root,bacteria,storage,standing dead N
C     TZSNC=cumulative PFT N litterfall
C     TZUPTK=cumulative PFT root-soil N exchange
C     TNH3C=cumulative NH3 exchange
C     RSETN=cumulative N balance from previous year
C     THVSTN=cumulative PFT N removed from ecosystem from
C        previous year
C     HVSTN=total PFT N removed from ecosystem in current year
C     VNH3F=NH3 emission from disturbance
C     TZUPFX=cumulative PFT N2 fixation
C
      BALN(NZ,NY,NX)=WTSHN(NZ,NY,NX)+WTRTN(NZ,NY,NX)+WTNDN(NZ,NY,NX)
     2+WTRVN(NZ,NY,NX)+TZSNC(NZ,NY,NX)-TZUPTK(NZ,NY,NX)
     3-TNH3C(NZ,NY,NX)-RSETN(NZ,NY,NX)+WTSTGN(NZ,NY,NX)
     4+HVSTN(NZ,NY,NX)+THVSTN(NZ,NY,NX)
     4-VNOXF(NZ,NY,NX)-TZUPFX(NZ,NY,NX)
C     IF(NZ.EQ.2)THEN
C     WRITE(*,1112)'BALN',I,J,NFZ,NX,NY,NZ
C    2,BALN(NZ,NY,NX),WTSHN(NZ,NY,NX),WTRTN(NZ,NY,NX)
C    2,WTNDN(NZ,NY,NX),WTRVN(NZ,NY,NX),TZSNC(NZ,NY,NX)
C    3,TZUPTK(NZ,NY,NX),TNH3C(NZ,NY,NX),RSETN(NZ,NY,NX)
C    4,WTSTGN(NZ,NY,NX),VNOXF(NZ,NY,NX),TZUPFX(NZ,NY,NX)
C    5,HVSTN(NZ,NY,NX),THVSTN(NZ,NY,NX)
C    4,(WTSTDN(M,NZ,NY,NX),M=1,5)
C    4,WTLFN(NZ,NY,NX),WTSHEN(NZ,NY,NX)
C    5,WTSTKN(NZ,NY,NX),WTRSVN(NZ,NY,NX),WTHSKN(NZ,NY,NX)
C    3,WTEARN(NZ,NY,NX),WTGRNN(NZ,NY,NX)
C    4,UPOMN(NZ,NY,NX),UPNH4(NZ,NY,NX),UPNO3(NZ,NY,NX)
C    5,(WTNDLN(L,NZ,NY,NX),L=NU(NY,NX),NI(NZ,NY,NX))
C    5,(ZPOOLR(1,L,NZ,NY,NX),L=NU(NY,NX),NI(NZ,NY,NX))
C    4,(((RDFOMN(N,K,L,NZ,NY,NX),N=1,2),K=0,4)
C    5,L=NU(NY,NX),NI(NZ,NY,NX))
C    4,((ZPOOLR(N,L,NZ,NY,NX),N=1,2),L=NU(NY,NX),NI(NZ,NY,NX))
1112  FORMAT(A8,6I4,200F18.6)
C     ENDIF
C
C     PLANT P BALANCE = TOTAL P STATE VARIABLES + TOTAL P LITTERFALL
C     - TOTAL P UPTAKE FROM SOIL
C
C     BALP=PFT N balance
C     WTSHP,WTRTP,WTNDP,WTRVP,WTSTGP=PFT
C        shoot,root,bacteria,storage,standing dead P
C     TPSNC=cumulative PFT P litterfall
C     TPUPTK=cumulative PFT root-soil P exchange
C     RSETP=cumulative P balance from previous year
C     THVSTP=cumulative PFT P removed from ecosystem from
C        previous year
C     HVSTP=total PFT P removed from ecosystem in current year
C     VPO4F=PO4 emission from disturbance
C
      BALP(NZ,NY,NX)=WTSHP(NZ,NY,NX)+WTRTP(NZ,NY,NX)+WTNDP(NZ,NY,NX)
     2+WTRVP(NZ,NY,NX)+TPSNC(NZ,NY,NX)-TPUPTK(NZ,NY,NX)
     3-RSETP(NZ,NY,NX)+WTSTDP(1,NZ,NY,NX)+WTSTGP(NZ,NY,NX)
     4+HVSTP(NZ,NY,NX)+THVSTP(NZ,NY,NX)-VPOXF(NZ,NY,NX)
C     IF(NZ.EQ.4)THEN
C     WRITE(*,1112)'BALP',I,J,NFZ,NX,NY,NZ
C    2,BALP(NZ,NY,NX),WTSHP(NZ,NY,NX)
C    2,WTRTP(NZ,NY,NX),WTNDP(NZ,NY,NX),WTRVP(NZ,NY,NX),TPSNC(NZ,NY,NX)
C    3,TPUPTK(NZ,NY,NX),RSETP(NZ,NY,NX)
C    4,WTSTDP(1,NZ,NY,NX),WTSTGP(NZ,NY,NX),HVSTP(NZ,NY,NX)
C    5,THVSTP(NZ,NY,NX),VPOXF(NZ,NY,NX)
C     ENDIF
      ENDIF
9975  CONTINUE
9990  CONTINUE
9995  CONTINUE
      RETURN
      END
