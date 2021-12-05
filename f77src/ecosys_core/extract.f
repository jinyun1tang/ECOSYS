      SUBROUTINE extract(I,J,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE AGGREGATES ALL SOIL-PLANT C,N,P EXCHANGES
C     FROM 'UPTAKE' AMD 'GROSUB' AND SENDS RESULTS TO 'REDIST'
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
      include "blk8a.h"
      include "blk8b.h"
      include "blk9a.h"
      include "blk9b.h"
      include "blk9c.h"
      include "blk11a.h"
      include "blk11b.h"
      include "blk12a.h"
      include "blk12b.h"
      include "blk13c.h"
      include "blk14.h"
      include "blk16.h"
      include "blk18a.h"
      include "blk18b.h"
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      DO 9985 NZ=1,NP0(NY,NX)
C
C     TOTAL LITTERFALL OF ALL PLANT SPECIES
C
C     ZCSNC,ZZSNC,ZPSNC=total C,N,P litterfall 
C     HCSNC,HZSNC,HPSNC=hourly PFT C,N,P litterfall from grosub.f 
C     WTSTGT=total standing dead C,N,P mass
C     WTSTG=PFT standing dead C,N,P mass
C     CSNC,ZSNC,PSNC=cumulative PFT C,N,P litterfall from grosub.f
C     CSNT,ZSNT,PSNT=cumulative total C,N,P litterfall
C
      ZCSNC(NY,NX)=ZCSNC(NY,NX)+HCSNC(NZ,NY,NX)
      ZZSNC(NY,NX)=ZZSNC(NY,NX)+HZSNC(NZ,NY,NX)
      ZPSNC(NY,NX)=ZPSNC(NY,NX)+HPSNC(NZ,NY,NX)
      WTSTGT(NY,NX)=WTSTGT(NY,NX)+WTSTG(NZ,NY,NX)
      DO 90 L=0,NI(NZ,NY,NX)
      DO 95 K=0,1
      DO 95 M=1,4
      CSNT(M,K,L,NY,NX)=CSNT(M,K,L,NY,NX)+CSNC(M,K,L,NZ,NY,NX)
      ZSNT(M,K,L,NY,NX)=ZSNT(M,K,L,NY,NX)+ZSNC(M,K,L,NZ,NY,NX)
      PSNT(M,K,L,NY,NX)=PSNT(M,K,L,NY,NX)+PSNC(M,K,L,NZ,NY,NX)
95    CONTINUE
90    CONTINUE
9985  CONTINUE
      ARLFC(NY,NX)=0.0
      ARSTC(NY,NX)=0.0
      DO 915 L=1,JC
      ARLFT(L,NY,NX)=0.0
      WGLFT(L,NY,NX)=0.0
      ARSTT(L,NY,NX)=0.0
915   CONTINUE
      DO 9975 NZ=1,NP(NY,NX)
      IF(IFLGC(NZ,NY,NX).EQ.1)THEN
C
C     TOTAL LEAF AREA OF ALL PLANT SPECIES
C
C     ARLFT,ARSTT=total leaf,stalk area of combined canopy layer
C     ARLFV,ARSTV=PFT leaf,stalk area in canopy layer
C     WGLFT=total leaf C of combined canopy layer
C     WGLFV=PFT leaf C in canopy layer
C
      DO 910 L=1,JC
      ARLFT(L,NY,NX)=ARLFT(L,NY,NX)+ARLFV(L,NZ,NY,NX)
      WGLFT(L,NY,NX)=WGLFT(L,NY,NX)+WGLFV(L,NZ,NY,NX)
      ARSTT(L,NY,NX)=ARSTT(L,NY,NX)+ARSTV(L,NZ,NY,NX)
910   CONTINUE
C
C     TOTAL GAS AND SOLUTE UPTAKE BY ALL PLANT SPECIES
C
      DO 100 N=1,MY(NZ,NY,NX)
      DO 100 L=NU(NY,NX),NI(NZ,NY,NX)
C
C     TOTAL ROOT DENSITY
C
C     RTDNT=total root length density 
C     RTDNP=PFT root length density per plant
C     UPWTR=total water uptake
C     UPWTR=PFT root water uptake
C     TUPHT=total convective heat in root water uptake
C     TKS=soil temperature
C     PP=PFT population
C
      IF(N.EQ.1)THEN
      RTDNT(L,NY,NX)=RTDNT(L,NY,NX)+RTDNP(N,L,NZ,NY,NX)
     2*PP(NZ,NY,NX)/AREA(3,L,NY,NX)
      ENDIF
C
C     TOTAL WATER UPTAKE
C
      TUPWTR(L,NY,NX)=TUPWTR(L,NY,NX)+UPWTR(N,L,NZ,NY,NX)
      TUPHT(L,NY,NX)=TUPHT(L,NY,NX)+UPWTR(N,L,NZ,NY,NX)
     2*4.19*TKS(L,NY,NX)
C
C     ROOT GAS CONTENTS FROM FLUXES IN 'UPTAKE'
C
C     *A,*P=PFT root gaseous, aqueous gas content
C     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
C     R*FLA=root gaseous-atmosphere CO2 exchange
C     R*DFA=root aqueous-gaseous CO2 exchange
C
      CO2A(N,L,NZ,NY,NX)=CO2A(N,L,NZ,NY,NX)+RCOFLA(N,L,NZ,NY,NX)
     2-RCODFA(N,L,NZ,NY,NX)
      OXYA(N,L,NZ,NY,NX)=OXYA(N,L,NZ,NY,NX)+ROXFLA(N,L,NZ,NY,NX)
     2-ROXDFA(N,L,NZ,NY,NX)
      CH4A(N,L,NZ,NY,NX)=CH4A(N,L,NZ,NY,NX)+RCHFLA(N,L,NZ,NY,NX)
     2-RCHDFA(N,L,NZ,NY,NX)
      Z2OA(N,L,NZ,NY,NX)=Z2OA(N,L,NZ,NY,NX)+RN2FLA(N,L,NZ,NY,NX)
     2-RN2DFA(N,L,NZ,NY,NX)
      ZH3A(N,L,NZ,NY,NX)=ZH3A(N,L,NZ,NY,NX)+RNHFLA(N,L,NZ,NY,NX)
     2-RNHDFA(N,L,NZ,NY,NX)
      H2GA(N,L,NZ,NY,NX)=H2GA(N,L,NZ,NY,NX)+RHGFLA(N,L,NZ,NY,NX)
     2-RHGDFA(N,L,NZ,NY,NX)
      CO2P(N,L,NZ,NY,NX)=CO2P(N,L,NZ,NY,NX)+RCODFA(N,L,NZ,NY,NX)
     2+RCO2P(N,L,NZ,NY,NX)
      OXYP(N,L,NZ,NY,NX)=OXYP(N,L,NZ,NY,NX)+ROXDFA(N,L,NZ,NY,NX)
     2-RUPOXP(N,L,NZ,NY,NX)
      CH4P(N,L,NZ,NY,NX)=CH4P(N,L,NZ,NY,NX)+RCHDFA(N,L,NZ,NY,NX)
     2+RUPCHS(N,L,NZ,NY,NX)
      Z2OP(N,L,NZ,NY,NX)=Z2OP(N,L,NZ,NY,NX)+RN2DFA(N,L,NZ,NY,NX)
     2+RUPN2S(N,L,NZ,NY,NX)
      ZH3P(N,L,NZ,NY,NX)=ZH3P(N,L,NZ,NY,NX)+RNHDFA(N,L,NZ,NY,NX)
     2+RUPN3S(N,L,NZ,NY,NX)+RUPN3B(N,L,NZ,NY,NX)
      H2GP(N,L,NZ,NY,NX)=H2GP(N,L,NZ,NY,NX)+RHGDFA(N,L,NZ,NY,NX)
     2+RUPHGS(N,L,NZ,NY,NX)
C
C     TOTAL ROOT GAS CONTENTS
C
C     TL*P=total root gas content
C     *A,*P=PFT root gaseous, aqueous gas content
C
      TLCO2P(L,NY,NX)=TLCO2P(L,NY,NX)+CO2P(N,L,NZ,NY,NX)
     2+CO2A(N,L,NZ,NY,NX)
      TLOXYP(L,NY,NX)=TLOXYP(L,NY,NX)+OXYP(N,L,NZ,NY,NX)
     2+OXYA(N,L,NZ,NY,NX)
      TLCH4P(L,NY,NX)=TLCH4P(L,NY,NX)+CH4P(N,L,NZ,NY,NX)
     2+CH4A(N,L,NZ,NY,NX)
      TLN2OP(L,NY,NX)=TLN2OP(L,NY,NX)+Z2OP(N,L,NZ,NY,NX)
     2+Z2OA(N,L,NZ,NY,NX)
      TLNH3P(L,NY,NX)=TLNH3P(L,NY,NX)+ZH3P(N,L,NZ,NY,NX)
     2+ZH3A(N,L,NZ,NY,NX)
      TLH2GP(L,NY,NX)=TLH2GP(L,NY,NX)+H2GP(N,L,NZ,NY,NX)
     2+H2GA(N,L,NZ,NY,NX)
C     IF(NX.EQ.1.AND.NY.EQ.1)THEN
C     WRITE(*,5566)'TLOXYP',I,J,NX,NY,NZ,L,N
C    2,TLOXYP(L,NY,NX),OXYP(N,L,NZ,NY,NX),OXYA(N,L,NZ,NY,NX) 
C    3,ROXDFA(N,L,NZ,NY,NX),RUPOXP(N,L,NZ,NY,NX),ROXFLA(N,L,NZ,NY,NX)
C    4,TLCO2P(L,NY,NX),CO2P(N,L,NZ,NY,NX),CO2A(N,L,NZ,NY,NX)
C    5,RCODFA(N,L,NZ,NY,NX),RCO2P(N,L,NZ,NY,NX)
C    6,RCO2S(N,L,NZ,NY,NX)
C    4,TLCH4P(L,NY,NX),CH4P(N,L,NZ,NY,NX),CH4A(N,L,NZ,NY,NX)
C    5,RCHDFA(N,L,NZ,NY,NX),RUPCHS(N,L,NZ,NY,NX),RCHFLA(N,L,NZ,NY,NX) 
5566  FORMAT(A8,7I4,20E12.4)
C     ENDIF
C
C     TOTAL ROOT BOUNDARY GAS FLUXES
C
C     T*FLA=total root gaseous-atmosphere CO2 exchange
C     R*FLA=PFT root gaseous-atmosphere CO2 exchange
C     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
C     TUP*S,TUP*B=total root-soil gas, solute exchange in non-band,band
C     RUP*S,RUP*B*=PFT root-soil gas, solute exchange in non-band,band
C     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
C     solute code:NH4=NH4,NO3=NO3,H2P=H2PO4,H1P=H1PO4 in non-band
C                :NHB=NH4,NOB=NO3,H2B=H2PO4,H1B=H1PO4 in band
C
      TCOFLA(L,NY,NX)=TCOFLA(L,NY,NX)+RCOFLA(N,L,NZ,NY,NX)
      TOXFLA(L,NY,NX)=TOXFLA(L,NY,NX)+ROXFLA(N,L,NZ,NY,NX)
      TCHFLA(L,NY,NX)=TCHFLA(L,NY,NX)+RCHFLA(N,L,NZ,NY,NX)
      TN2FLA(L,NY,NX)=TN2FLA(L,NY,NX)+RN2FLA(N,L,NZ,NY,NX)
      TNHFLA(L,NY,NX)=TNHFLA(L,NY,NX)+RNHFLA(N,L,NZ,NY,NX)
      THGFLA(L,NY,NX)=THGFLA(L,NY,NX)+RHGFLA(N,L,NZ,NY,NX)
      TCO2P(L,NY,NX)=TCO2P(L,NY,NX)-RCO2P(N,L,NZ,NY,NX)
      TUPOXP(L,NY,NX)=TUPOXP(L,NY,NX)+RUPOXP(N,L,NZ,NY,NX)
      TCO2S(L,NY,NX)=TCO2S(L,NY,NX)+RCO2S(N,L,NZ,NY,NX)
      TUPOXS(L,NY,NX)=TUPOXS(L,NY,NX)+RUPOXS(N,L,NZ,NY,NX)
      TUPCHS(L,NY,NX)=TUPCHS(L,NY,NX)+RUPCHS(N,L,NZ,NY,NX)
      TUPN2S(L,NY,NX)=TUPN2S(L,NY,NX)+RUPN2S(N,L,NZ,NY,NX)
      TUPN3S(L,NY,NX)=TUPN3S(L,NY,NX)+RUPN3S(N,L,NZ,NY,NX)
      TUPN3B(L,NY,NX)=TUPN3B(L,NY,NX)+RUPN3B(N,L,NZ,NY,NX)
      TUPHGS(L,NY,NX)=TUPHGS(L,NY,NX)+RUPHGS(N,L,NZ,NY,NX)
      TUPNH4(L,NY,NX)=TUPNH4(L,NY,NX)+RUPNH4(N,L,NZ,NY,NX)
      TUPNO3(L,NY,NX)=TUPNO3(L,NY,NX)+RUPNO3(N,L,NZ,NY,NX)
      TUPH2P(L,NY,NX)=TUPH2P(L,NY,NX)+RUPH2P(N,L,NZ,NY,NX)
      TUPH1P(L,NY,NX)=TUPH1P(L,NY,NX)+RUPH1P(N,L,NZ,NY,NX)
      TUPNHB(L,NY,NX)=TUPNHB(L,NY,NX)+RUPNHB(N,L,NZ,NY,NX)
      TUPNOB(L,NY,NX)=TUPNOB(L,NY,NX)+RUPNOB(N,L,NZ,NY,NX)
      TUPH2B(L,NY,NX)=TUPH2B(L,NY,NX)+RUPH2B(N,L,NZ,NY,NX)
      TUPH1B(L,NY,NX)=TUPH1B(L,NY,NX)+RUPH1B(N,L,NZ,NY,NX)
C     IF(NY.EQ.5.AND.L.EQ.1)THEN
C     WRITE(*,4141)'TCO2S',I,J,NX,NY,L,NZ,N,TCO2S(L,NY,NX)
C    2,RCO2S(N,L,NZ,NY,NX),TCO2P(L,NY,NX),RCO2P(N,L,NZ,NY,NX)
C     WRITE(*,4141)'TUPOX',I,J,NX,NY,L,NZ,N,TUPOXS(L,NY,NX)
C    2,RUPOXS(N,L,NZ,NY,NX),TUPOXP(L,NY,NX),RUPOXP(N,L,NZ,NY,NX)
C     WRITE(*,4141)'TCHFLA',I,J,NX,NY,L,NZ,N,TCHFLA(L,NY,NX)
C    2,RCHFLA(N,L,NZ,NY,NX),RUPCHS(N,L,NZ,NY,NX)
C     WRITE(*,4141)'TUPN2S',I,J,NX,NY,L,NZ,N,TUPN2S(L,NY,NX)
C    2,RUPN2S(N,L,NZ,NY,NX) 
4141  FORMAT(A8,7I4,12E12.4)
C     ENDIF
C
C     TOTAL ROOT C,N,P EXUDATION
C
C     TDFOMC,TDFOMN,TDFOMP=total nonstructl C,N,P exchange
C     RDFOMC,RDFOMN,RDFOMP=PFT nonstructl C,N,P exchange 
C
      DO 195 K=0,4
      TDFOMC(K,L,NY,NX)=TDFOMC(K,L,NY,NX)-RDFOMC(N,K,L,NZ,NY,NX)
      TDFOMN(K,L,NY,NX)=TDFOMN(K,L,NY,NX)-RDFOMN(N,K,L,NZ,NY,NX)
      TDFOMP(K,L,NY,NX)=TDFOMP(K,L,NY,NX)-RDFOMP(N,K,L,NZ,NY,NX)
195   CONTINUE
C
C     TOTAL ROOT O2, NH4, NO3, PO4 UPTAKE CONTRIBUTES TO
C     TOTAL ROOT + MICROBIAL UPTAKE USED TO CALCULATE
C     COMPETITION CONSTRAINTS
C
C     ROXYX=O2 demand by all microbial,root,myco populations 
C     RNH4X=NH4 demand in non-band by all microbial,root,myco populations
C     RNO3X=NO3 demand in non-band by all microbial,root,myco populations
C     RPO4X=H2PO4 demand in non-band by all microbial,root,myco populations
C     RP14X=HPO4 demand in non-band by all microbial,root,myco populations
C     RNHBX=NH4 demand in band by all microbial,root,myco populations
C     RN3BX=NO3 demand in band by all microbial,root,myco populations
C     RPOBX=H2PO4 demand in band by all microbial,root,myco populations
C     RP1BX=HPO4 demand in band by all microbial,root,myco populations
C     ROXYP=O2 demand by each root,myco population
C     RUNNHP=NH4 demand in non-band by each root population
C     RUNNOP=NO3 demand in non-band by each root population
C     RUPP2P=H2PO4 demand in non-band by each root population
C     RUPP1P=HPO4 demand in non-band by each root population
C     RUNNBP=NH4 demand in band by each root population
C     RUNNXB=NO3 demand in band by each root population
C     RUPP2B=H2PO4 demand in band by each root population
C     RUPP1B=HPO4 demand in band by each root population
C
      ROXYX(L,NY,NX)=ROXYX(L,NY,NX)+ROXYP(N,L,NZ,NY,NX)
      RNH4X(L,NY,NX)=RNH4X(L,NY,NX)+RUNNHP(N,L,NZ,NY,NX)
      RNO3X(L,NY,NX)=RNO3X(L,NY,NX)+RUNNOP(N,L,NZ,NY,NX)
      RPO4X(L,NY,NX)=RPO4X(L,NY,NX)+RUPP2P(N,L,NZ,NY,NX)
      RP14X(L,NY,NX)=RP14X(L,NY,NX)+RUPP1P(N,L,NZ,NY,NX)
      RNHBX(L,NY,NX)=RNHBX(L,NY,NX)+RUNNBP(N,L,NZ,NY,NX)
      RN3BX(L,NY,NX)=RN3BX(L,NY,NX)+RUNNXP(N,L,NZ,NY,NX)
      RPOBX(L,NY,NX)=RPOBX(L,NY,NX)+RUPP2B(N,L,NZ,NY,NX)
      RP1BX(L,NY,NX)=RP1BX(L,NY,NX)+RUPP1B(N,L,NZ,NY,NX)
100   CONTINUE
C     IF(ISALTG.NE.0)THEN
C     XZHYS(L,NY,NX)=XZHYS(L,NY,NX)
C    2+0.0714*(TUPNH4(L,NY,NX)+TUPNHB(L,NY,NX))
C    3-0.0714*(TUPNO3(L,NY,NX)+TUPNOB(L,NY,NX))
C     ENDIF
C
C     TOTAL ROOT N2 FIXATION BY ALL PLANT SPECIES
C
C     TUPNF=total root N2 fixation
C     RUPNF=PFT root N2 fixation
C
      DO 85 L=NU(NY,NX),NI(NZ,NY,NX)
      TUPNF(L,NY,NX)=TUPNF(L,NY,NX)+RUPNF(L,NZ,NY,NX)
85    CONTINUE
C
C     TOTAL ENERGY, WATER, CO2 FLUXES
C
C     TRN=total net SW+LW absorbed by canopy
C     RAD1=PFT net SW+LW absorbed by canopy
C     TLE=total canopy latent heat flux
C     EFLXC=PFT canopy latent heat flux
C     TSH=total canopy sensible heat flux
C     SFLXC=PFT canopy sensible heat flux
C     TGH=total canopy storage heat flux
C     HFLXC=PFT canopy storage heat flux
C     TCCAN=total net CO2 fixation
C     CNET=PFT net CO2 fixation
C     TVOLWP,TVOLWC=total water volume in canopy,on canopy surfaces
C     VOLWP,VOLWC=PFT water volume in canopy,on canopy surfaces
C     TEVAPP,TEVAPC=total water flux to,from canopy,canopy surfaces 
C     EVAPC,EP=water flux to,from canopy surfaces, inside canopy
C     TENGYC=total canopy water heat content
C     ENGYC=PFT canopy water heat content
C     ARLFC,ARSTC=total leaf,stalk area
C     ARLFP,ARSTP=PFT leaf,stalk area
C     ZCSNC,ZZSNC,ZPSNC=total net root-soil C,N,P exchange 
C     HCUPTK,HZUPTK,HPUPTK=PFT net root-soil C,N,P exchange 
C     TBALC,TBALN,TBALP=total C,N,P balance
C     BALC,BALN,BALP=PFT C,N,P balance
C     TCO2Z,TOXYZ,TCH4Z,TN2OZ,TNH3Z,TH2GZ=total loss of root CO2, O2, CH4, N2O, NH3, H2
C     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=PFT loss of root CO2, O2, CH4, N2O, NH3, H2
C
      TRN(NY,NX)=TRN(NY,NX)+RAD1(NZ,NY,NX)
      TLE(NY,NX)=TLE(NY,NX)+EFLXC(NZ,NY,NX)
      TSH(NY,NX)=TSH(NY,NX)+SFLXC(NZ,NY,NX)
      TGH(NY,NX)=TGH(NY,NX)+HFLXC(NZ,NY,NX)
      TCCAN(NY,NX)=TCCAN(NY,NX)+CNET(NZ,NY,NX)
      CTRAN(NZ,NY,NX)=CTRAN(NZ,NY,NX)+EP(NZ,NY,NX)+EVAPC(NZ,NY,NX)
      TVOLWP(NY,NX)=TVOLWP(NY,NX)+VOLWP(NZ,NY,NX)
      TVOLWC(NY,NX)=TVOLWC(NY,NX)+VOLWC(NZ,NY,NX)
      TEVAPP(NY,NX)=TEVAPP(NY,NX)+EP(NZ,NY,NX)+EVAPC(NZ,NY,NX)
      TEVAPC(NY,NX)=TEVAPC(NY,NX)+EVAPC(NZ,NY,NX)
      ENGYC=4.19*(VOLWC(NZ,NY,NX)+FLWC(NZ,NY,NX)+EVAPC(NZ,NY,NX))
     2*TKC(NZ,NY,NX)
      TENGYC(NY,NX)=TENGYC(NY,NX)+ENGYC
      THFLXC(NY,NX)=THFLXC(NY,NX)+ENGYC-ENGYX(NZ,NY,NX)
     2-(FLWC(NZ,NY,NX)*4.19*TKA(NY,NX))
      ENGYX(NZ,NY,NX)=ENGYC
      THRMC(NY,NX)=THRMC(NY,NX)+THRM1(NZ,NY,NX)
      ARLFC(NY,NX)=ARLFC(NY,NX)+ARLFP(NZ,NY,NX)
      ARSTC(NY,NX)=ARSTC(NY,NX)+ARSTP(NZ,NY,NX)
      ZCSNC(NY,NX)=ZCSNC(NY,NX)-HCUPTK(NZ,NY,NX)
      ZZSNC(NY,NX)=ZZSNC(NY,NX)-HZUPTK(NZ,NY,NX)
      ZPSNC(NY,NX)=ZPSNC(NY,NX)-HPUPTK(NZ,NY,NX)
      TBALC=TBALC+BALC(NZ,NY,NX)
      TBALN=TBALN+BALN(NZ,NY,NX)
      TBALP=TBALP+BALP(NZ,NY,NX)
      TCO2Z(NY,NX)=TCO2Z(NY,NX)+RCO2Z(NZ,NY,NX)
      TOXYZ(NY,NX)=TOXYZ(NY,NX)+ROXYZ(NZ,NY,NX)
      TCH4Z(NY,NX)=TCH4Z(NY,NX)+RCH4Z(NZ,NY,NX)
      TN2OZ(NY,NX)=TN2OZ(NY,NX)+RN2OZ(NZ,NY,NX)
      TNH3Z(NY,NX)=TNH3Z(NY,NX)+RNH3Z(NZ,NY,NX)
      TH2GZ(NY,NX)=TH2GZ(NY,NX)+RH2GZ(NZ,NY,NX)
C
C     TOTAL CANOPY NH3 EXCHANGE AND EXUDATION
C
C     RNH3B,RNH3C=PFT NH3 flux between atmosphere and branch,canopy
C     TNH3C=total NH3 flux between atmosphere and canopy
C
      RNH3C(NZ,NY,NX)=0.0
      DO 80 NB=1,NBR(NZ,NY,NX)
      RNH3C(NZ,NY,NX)=RNH3C(NZ,NY,NX)+RNH3B(NB,NZ,NY,NX)
      TNH3C(NZ,NY,NX)=TNH3C(NZ,NY,NX)+RNH3B(NB,NZ,NY,NX)
80    CONTINUE
      ENDIF
9975  CONTINUE
9990  CONTINUE
9995  CONTINUE
      RETURN
      END
