      SUBROUTINE extract(I,J,NFZ,NHW,NHE,NVN,NVS)
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
      include "blk15a.h"
      include "blk16.h"
      include "blk18a.h"
      include "blk18b.h"
      include "blk19a.h"
C
C     FCHMCX=maximum fraction of anaerobic C combustion to charcoal
C     FCOMCX=maximum fraction of aerobic C combustion to charcoal
C     FCOMNX=maximum fraction of N combustion to NH4
C     FCOMPX=maximum fraction of P combustion to H2PO4
C     FCOMNY=minimum fraction of N combustion to NH4
C     FCOMPY=minimum fraction of P combustion to H2PO4
C
      PARAMETER (FCOMCX=0.0,FCHMCX=0.5)
      PARAMETER (FCOMNX=0.4,FCOMPX=0.9,FCOMNY=0.1,FCOMPY=0.7) 
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      DO 9985 NZ=1,NP0(NY,NX)
C
C     TOTAL LITTERFALL OF ALL PLANT SPECIES
C
C     ZCSNC,ZZSNC,ZPSNC=C,N,P litterfall from all PFT (g C,N,P t-1)
C     HCSNC,HZSNC,HPSNC=PFT C,N,P litterfall from ‘grosub.f’ 
C        (g C,N,P t-1)
C     WTSTGT=total standing dead C,N,P mass (g C,N,P)
C     WTSTG=PFT standing dead C,N,P mass (g C,N,P)
C     CSNC,ZSNC,PSNC=PFT C,N,P litterfall from ‘grosub.f’ (g C,N,P t-1)
C     CSNT,ZSNT,PSNT=C,N,P litterfall from all PFT (g C,N,P t-1)
C
      ZCSNC(NY,NX)=ZCSNC(NY,NX)+HCSNC(NZ,NY,NX)
      ZZSNC(NY,NX)=ZZSNC(NY,NX)+HZSNC(NZ,NY,NX)
      ZPSNC(NY,NX)=ZPSNC(NY,NX)+HPSNC(NZ,NY,NX)
      WTSTGT(NY,NX)=WTSTGT(NY,NX)+WTSTG(NZ,NY,NX)
      DO 90 L=0,NJ(NY,NX)
      DO 95 K=0,1
      DO 95 M=1,5
      CSNT(M,K,L,NY,NX)=CSNT(M,K,L,NY,NX)+CSNC(M,K,L,NZ,NY,NX)
      ZSNT(M,K,L,NY,NX)=ZSNT(M,K,L,NY,NX)+ZSNC(M,K,L,NZ,NY,NX)
      PSNT(M,K,L,NY,NX)=PSNT(M,K,L,NY,NX)+PSNC(M,K,L,NZ,NY,NX)
C     IF(NZ.EQ.2)THEN
C     WRITE(*,8484)'CSNC',I,J,NFZ,NX,NY,NZ,L,K,M
C    2,CSNT(M,K,L,NY,NX),CSNC(M,K,L,NZ,NY,NX)
C    2,ZSNT(M,K,L,NY,NX),ZSNC(M,K,L,NZ,NY,NX)
C    2,PSNT(M,K,L,NY,NX),PSNC(M,K,L,NZ,NY,NX)
8484  FORMAT(A8,9I4,12E14.6)
C     ENDIF
95    CONTINUE
C
C     SALTS IN LITTERFALL
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C     *SNC=PFT salt litterfall from ‘grosub.f’ (mol t-1)
C     *SNT=salt litterfall from all PFT (mol t-1)
C     salt code:AL=Al,FE=Fe,CA=Ca,GM=Mg,AN=Na,AK=K,SO=SO4,CL=Cl
C
      IF(ISALTG.NE.0)THEN
      ALSNT(L,NY,NX)=ALSNT(L,NY,NX)+ALSNC(L,NZ,NY,NX)
      FESNT(L,NY,NX)=FESNT(L,NY,NX)+FESNC(L,NZ,NY,NX)
      CASNT(L,NY,NX)=CASNT(L,NY,NX)+CASNC(L,NZ,NY,NX)
      GMSNT(L,NY,NX)=GMSNT(L,NY,NX)+GMSNC(L,NZ,NY,NX)
      ANSNT(L,NY,NX)=ANSNT(L,NY,NX)+ANSNC(L,NZ,NY,NX)
      AKSNT(L,NY,NX)=AKSNT(L,NY,NX)+AKSNC(L,NZ,NY,NX)
      SOSNT(L,NY,NX)=SOSNT(L,NY,NX)+SOSNC(L,NZ,NY,NX)
      CLSNT(L,NY,NX)=CLSNT(L,NY,NX)+CLSNC(L,NZ,NY,NX)
C     WRITE(*,1194)'ALSNT',I,J,NFZ,NX,NY,NZ,L
C    2,ALSNT(L,NY,NX),ALSNC(L,NZ,NY,NX) 
C     WRITE(*,1194)'CASNT',I,J,NFZ,NX,NY,NZ,L
C    2,CASNT(L,NY,NX),CASNC(L,NZ,NY,NX) 
C     WRITE(*,1194)'SOSNT',I,J,NFZ,NX,NY,NZ,L
C    2,SOSNT(L,NY,NX),SOSNC(L,NZ,NY,NX) 
1194  FORMAT(A8,7I4,20F16.8)
      ENDIF
90    CONTINUE
C
C     TOTAL MANURE DEPOSITION OF ALL PLANT SPECIES
C
C     CSNMT,ZSNMT,PSNMT=total manure organic C,N,P deposition 
C        (g C,N,P h-1)
C     CSNM,ZSNM,PSNM=manure organic C,N,P deposition from ‘grosub.f’
C        (g C,N,P h-1)
C     ZSNIT,PSNIT=total manure inorganic N,P deposition 
C        (g C,N,P h-1)  
C     ZSNI,PSNI=manure inorganic N,P deposition from ‘grosub.f’
C        (g C,N,P h-1)  
C
      IF(NFZ.EQ.1)THEN
      DO 96 M=1,4
      CSNMT(M,NY,NX)=CSNMT(M,NY,NX)+CSNM(M,NZ,NY,NX)
      ZSNMT(M,NY,NX)=ZSNMT(M,NY,NX)+ZSNM(M,NZ,NY,NX)
      PSNMT(M,NY,NX)=PSNMT(M,NY,NX)+PSNM(M,NZ,NY,NX)
96    CONTINUE
      ZSNIT(NY,NX)=ZSNIT(NY,NX)+ZSNI(NZ,NY,NX)
      PSNIT(NY,NX)=PSNIT(NY,NX)+PSNI(NZ,NY,NX)
      ENDIF
9985  CONTINUE
      ARLFC(NY,NX)=0.0
      ARSTC(NY,NX)=0.0
      ARSDC(NY,NX)=0.0
      DO 915 L=1,JC
      ARLFT(L,NY,NX)=0.0
      WGLFT(L,NY,NX)=0.0
      ARSTT(L,NY,NX)=0.0
      ARSDT(L,NY,NX)=0.0
915   CONTINUE
C
C     FIRE
C
C     ICHKF=fire flag:0=no fire,1=fire)
C
      IF(ICHKF.EQ.1)THEN
C
C     TOTAL CANOPY C COMBUSTION AND O2,CO2,CH4 FLUXES
C
C     RCGCK=total canopy+standing dead C combustion of all PFT 
C        from ‘grosub.f’(g C t-1)
C     TFNCO=temperature function for combusting living C (600K = 1)
C     FCOMC=fraction of aerobic C combustion to charcoal
C     FCHMC=fraction of anaerobic C combustion to charcoal
C     ROGCK=O2-unlimited O2 consumption from RCGCK (g O t-1)
C     RCGOK=O2-limited O2 consumption (g O t-1)
C     OXYQ,OXYC=canopy air O2 concentration,content (umol mol-1,g O)
C     COXYGK,CCH4GK=M-M constant for O2,CH4 combustion (umol mol-1)
C     RCGOK=O2-limited total C combustion of all PFT (g C t-1)
C     RCHOK=CH4,O2-unlimited CH4 combustion (g C t-1)
C     FCOGC,FCHGC=fraction of aerobic,anaerobic C combustion to CO2
C     FCC=fraction of aerobic C combustion to charcoal
C     RC4OK=CH4,O2-limited CH4 combustion (g C t-1)
C     CH4Q=canopy air CH4 concentration (umol mol-1)
C
      IF(RCGCK(NY,NX).GT.ZEROS2(NY,NX))THEN
      RTK=8.3143*TKQ(NY,NX)
      TFNCO=AMIN1(1.0,EXP(12.028-60000/RTK))
      FCOMC=FCOMCX*(1.0-TFNCO)
      FCOGC=1.0-FCOMC
      FCHMC=FCHMCX*(1.0-TFNCO)
      FCHGC=1.0-FCHMC
      ROGCK=RCGCK(NY,NX)*2.667
      ROGOK=AMIN1(ROGCK*OXYQ(NY,NX)/(OXYQ(NY,NX)+COXYGK)
     2,OXYC(NY,NX))*FCOGC
      RCGOK=ROGOK/2.667
      RCHOK=(RCGCK(NY,NX)*FCOGC-RCGOK)*FCHGC
      FCC=(RCGCK(NY,NX)-RCGOK-RCHOK)/RCGCK(NY,NX) 
      RC4OK=AMIN1(RCHOK*CH4Q(NY,NX)/(CH4Q(NY,NX)+CCH4GK)
     2*OXYQ(NY,NX)/(OXYQ(NY,NX)+COXYGK),(OXYC(NY,NX)-ROGOK)/2.667)
C
C     NET O2,CO2,CH4 FLUXES
C
C     XCNET,XHNET,XONET=total CO2,CH4,O2 exchange by all PFT canopies
C        (g C,C,O t-1)
C     FCBOX,FCBCH=fraction of C combusted aerobically,anaerobically
C     HCBFC=heat released by canopy+standing dead combustion (MJ t-1) 
C     GCBCO,GCBCH,GCBC4=combustion energy of organic C (aerobic,
C        anaerobic), and CH4 from ‘starts.f’ (MJ g C-1)
C
      XONET(NY,NX)=XONET(NY,NX)+ROGOK+RC4OK*2.667
      XCNET(NY,NX)=XCNET(NY,NX)-RCGOK-RC4OK
      XHNET(NY,NX)=XHNET(NY,NX)-RCHOK+RC4OK
      FCBOX=RCGOK/(RCGOK+RCHOK)
      FCBCH=RCHOK/(RCGOK+RCHOK)
      HCBFC=RCGOK*GCBCO+RCHOK*GCBCH+RC4OK*GCBC4 
C     WRITE(*,1193)'XHNET',I,J,NFZ,NX,NY,XHNET(NY,NX),RCHOK,RC4OK
C    2,FCOGC,RCGOK,FCHGC,CH,CH4Q(NY,NX),CCH4GK,OXYQ(NY,NX),COXYGK
C    3,TKQ(NY,NX),TFNCO,OXYC(NY,NX)
1193  FORMAT(A8,5I4,20E12.4)
      ELSE
      ROGOK=0.0
      RCGOK=0.0
      RCHOK=0.0
      RC4OK=0.0
      FCC=0.0
      FCBOX=0.0
      FCBCH=0.0
      HCBFC=0.0
      ENDIF
      DO 9970 NZ=1,NP(NY,NX)
      COX=0.0
      CHX=0.0
      ZOX=0.0
      POX=0.0
      COR=0.0
      Z4M=0.0
      P4M=0.0
C
C     COMBUST SHOOT
C
C     FCOMN,FGOMN=fraction of N combustion to NH4, gaseous NOX
C     FCOMP,FCOGP=fraction of P combustion to H2PO4, gaseous POX
C     TFNCOC,TFNCOD=temperature function used to partition combusted
C        living,standing dead N,P from ‘grosub.f’ 
C     HCBFCZ=heat released by living canopy combustion (MJ t-1)
C     RCGCKZ,RCGDKZ=total living,standing dead C combustion 
C        from ‘grosub.f’ (g C t-1)
C     R*=C,N,P combustion from each PFT C pool from ‘grosub.f’(g C t-1)
C        PFT pool code (g C,N,P):
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
C     WTSTDG(5,=standing dead charcoal (g C)
C     FCC=fraction of aerobic C combustion to charcoal
C     ZNH4S(0=litter NH4 (g N)
C     H2PO4(0=litter H2PO4 (g P)
C     FCBOX,FCBCH=fraction of C combusted aerobically,anaerobically
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C        (g C,C,N,P t-1)
C     COR=charcoal addition to litter (g C t-1)
C     Z4M,P4M=N,P mineral additions to litter (g N,P t-1)
C
      FCOMN=FCOMNY+(FCOMNX-FCOMNY)*(1.0-TFNCOC(NZ,NY,NX))
      FCOGN=1.0-FCOMN
      FCOMP=FCOMPY+(FCOMPX-FCOMPY)*(1.0-TFNCOC(NZ,NY,NX))
      FCOGP=1.0-FCOMP
      IF(RCGCKZ(NZ,NY,NX)+RCGDKZ(NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
      HCBFCZ(NZ,NY,NX)=HCBFC*RCGCKZ(NZ,NY,NX)/RCGCK(NY,NX) 
      DO 8750 NB=1,NBR(NZ,NY,NX)
      RCMBC=RCPOOL(NB,NZ,NY,NX)
     2+RWTLFB(NB,NZ,NY,NX)+RWTSHEB(NB,NZ,NY,NX)+RWTSTKB(NB,NZ,NY,NX)
     3+RWTRSVB(NB,NZ,NY,NX)+RWTHSKB(NB,NZ,NY,NX)+RWTEARB(NB,NZ,NY,NX)
     4+RWTGRB(NB,NZ,NY,NX)+RCPOLNB(NB,NZ,NY,NX)+RWTNDB(NB,NZ,NY,NX)
      RCMBN=RZPOOL(NB,NZ,NY,NX)
     2+RWTLFBN(NB,NZ,NY,NX)+RWTSHBN(NB,NZ,NY,NX)+RWTSTBN(NB,NZ,NY,NX)
     3+RWTRSBN(NB,NZ,NY,NX)+RWTHSBN(NB,NZ,NY,NX)+RWTEABN(NB,NZ,NY,NX)
     4+RWTGRBN(NB,NZ,NY,NX)+RZPOLNB(NB,NZ,NY,NX)+RWTNDBN(NB,NZ,NY,NX)
      RCMBP=RPPOOL(NB,NZ,NY,NX)
     2+RWTLFBP(NB,NZ,NY,NX)+RWTSHBP(NB,NZ,NY,NX)+RWTSTBP(NB,NZ,NY,NX)
     3+RWTRSBP(NB,NZ,NY,NX)+RWTHSBP(NB,NZ,NY,NX)+RWTEABP(NB,NZ,NY,NX)
     4+RWTGRBP(NB,NZ,NY,NX)+RPPOLNB(NB,NZ,NY,NX)+RWTNDBP(NB,NZ,NY,NX)
      WTSTDG(5,NZ,NY,NX)=WTSTDG(5,NZ,NY,NX)+RCMBC*FCC
      ZNH4S(0,NY,NX)=ZNH4S(0,NY,NX)+RCMBN*FCOMN     
      H2PO4(0,NY,NX)=H2PO4(0,NY,NX)+RCMBP*FCOMP
      COX=COX+RCMBC*(1.0-FCC)*FCBOX
      CHX=CHX+RCMBC*(1.0-FCC)*FCBCH 
      COR=COR+RCMBC*FCC 
      ZOX=ZOX+RCMBN*FCOGN 
      POX=POX+RCMBP*FCOGP 
      Z4M=Z4M+RCMBN*FCOMN 
      P4M=P4M+RCMBP*FCOMP
C
C     SALT FROM SHOOT
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C     R*Q=salt combustion rate for each PFT shoot pool (mol t-1)
C     Z*=salt in soil (mol)
C     salt code:AL=Al,FE=Fe,CA=Ca,MG=Mg,NA=Na,KA=K,SO=SO4,CL=Cl
C
      IF(ISALTG.NE.0)THEN
      ZAL(0,NY,NX)=ZAL(0,NY,NX)+RZALQ(NB,NZ,NY,NX)
      ZFE(0,NY,NX)=ZFE(0,NY,NX)+RZFEQ(NB,NZ,NY,NX)
      ZCA(0,NY,NX)=ZCA(0,NY,NX)+RZCAQ(NB,NZ,NY,NX)
      ZMG(0,NY,NX)=ZMG(0,NY,NX)+RZMGQ(NB,NZ,NY,NX)
      ZNA(0,NY,NX)=ZNA(0,NY,NX)+RZNAQ(NB,NZ,NY,NX)
      ZKA(0,NY,NX)=ZKA(0,NY,NX)+RZKAQ(NB,NZ,NY,NX)
      ZSO4(0,NY,NX)=ZSO4(0,NY,NX)+RZSOQ(NB,NZ,NY,NX)
      ZCL(0,NY,NX)=ZCL(0,NY,NX)+RZCLQ(NB,NZ,NY,NX)
      TIONIN=TIONIN+RZALQ(NB,NZ,NY,NX)+RZFEQ(NB,NZ,NY,NX)
     2+RZCAQ(NB,NZ,NY,NX)+RZMGQ(NB,NZ,NY,NX)+RZNAQ(NB,NZ,NY,NX)
     5+RZKAQ(NB,NZ,NY,NX)+RZSOQ(NB,NZ,NY,NX)+RZCLQ(NB,NZ,NY,NX)
C     WRITE(*,5571)'RZQ',I,J,NFZ,NX,NY,NZ,NB
C    2,TIONIN,RZALQ(NB,NZ,NY,NX),RZFEQ(NB,NZ,NY,NX)
C    2,RZCAQ(NB,NZ,NY,NX),RZMGQ(NB,NZ,NY,NX),RZNAQ(NB,NZ,NY,NX)
C    5,RZKAQ(NB,NZ,NY,NX),RZSOQ(NB,NZ,NY,NX),RZCLQ(NB,NZ,NY,NX)
5571  FORMAT(A8,7I4,50F16.8)
      ENDIF
8750  CONTINUE
      ENDIF
C
C     COMBUST STANDING DEAD
C
C     HCBFDZ=heat released by standing dead combustion (MJ t-1)
C     HCBFC=heat released by living+standing dead combustion (MJ t-1) 
C     RCGCKZ,RCGDKZ=total canopy,standing dead C combustion (g C t-1)
C     R*=C,N,P combustion rate for each PFT C pool (g C,N,P t-1)
C        PFT pool code (g C,N,P):
C        WTSTDG,WTSTDN,WTSTDP=standing dead C,N,P mass
C     WTSTDG(5,=standing dead charcoal (g C)
C     FCC=fraction of aerobic C combustion to charcoal
C     ZNH4S(0=litter NH4 (g N)
C     H2PO4(0=litter H2PO4 (g P)
C     FCBOX,FCBCH=fraction of C combusted aerobically,anaerobically
C     FCOMN,FGOMN=fraction of N combustion to NH4, gaseous NOX
C     FCOMP,FCOGP=fraction of P combustion to H2PO4, gaseous POX
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C        (g C,C,N,P t-1)
C     Z4M,P4M=N,P mineral additions to litter (g N,P t-1)
C
      IF(RCGCKZ(NZ,NY,NX)+RCGDKZ(NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
      HCBFDZ(NZ,NY,NX)=HCBFC*RCGDKZ(NZ,NY,NX)/RCGCK(NY,NX) 
      FCOMN=FCOMNY+(FCOMNX-FCOMNY)*(1.0-TFNCOD(NZ,NY,NX))
      FCOGN=1.0-FCOMN
      FCOMP=FCOMPY+(FCOMPX-FCOMPY)*(1.0-TFNCOD(NZ,NY,NX))
      FCOGP=1.0-FCOMP
      DO 8740 M=1,4
      WTSTDG(5,NZ,NY,NX)=WTSTDG(5,NZ,NY,NX)+RWTSTDG(M,NZ,NY,NX)*FCC
      ZNH4S(0,NY,NX)=ZNH4S(0,NY,NX)+RWTSTDN(M,NZ,NY,NX)*FCOMN     
      H2PO4(0,NY,NX)=H2PO4(0,NY,NX)+RWTSTDP(M,NZ,NY,NX)*FCOMP
      COX=COX+RWTSTDG(M,NZ,NY,NX)*(1.0-FCC)*FCBOX
      CHX=CHX+RWTSTDG(M,NZ,NY,NX)*(1.0-FCC)*FCBCH
      COR=COR+RWTSTDG(M,NZ,NY,NX)*FCC 
      ZOX=ZOX+RWTSTDN(M,NZ,NY,NX)*FCOGN
      POX=POX+RWTSTDP(M,NZ,NY,NX)*FCOGP 
      Z4M=Z4M+RWTSTDN(M,NZ,NY,NX)*FCOMN 
      P4M=P4M+RWTSTDP(M,NZ,NY,NX)*FCOMP
C     IF(NZ.EQ.1)THEN
C     WRITE(*,5568)'DEAD',I,J,NFZ,NX,NY,NZ,M
C    2,RWTSTDG(M,NZ,NY,NX),RWTSTDN(M,NZ,NY,NX),RWTSTDP(M,NZ,NY,NX)
C    3,FCOMN,FCOMP,COX,CHX,COR,ZOX,POX,Z4M,P4M
C    4,ZNH4S(0,NY,NX),H2PO4(0,NY,NX)
5568  FORMAT(A8,7I4,50F16.8)
C     ENDIF
8740  CONTINUE
C
C     COMBUST STANDING DEAD CHARCOAL
C
C     R*=C,N,P combustion rate for each PFT C pool from ‘grosub.f’ 
C        (g C,N,P t-1)
C        PFT pool code (g C,N,P):
C        WTSTDG(5,WTSTDN(5,WTSTDP(5=charcoal C,N,P mass
C     FCC=fraction of aerobic C combustion to charcoal
C     ZNH4S(0=litter NH4 (g N)
C     H2PO4(0=litter H2PO4 (g P)
C     FCBOX,FCBCH=fraction of C combusted aerobically,anaerobically
C     FCOMN,FGOMN=fraction of N combustion to NH4, gaseous NOX
C     FCOMP,FCOGP=fraction of P combustion to H2PO4, gaseous POX
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C        (g C,C,N,P t-1)
C     Z4M,P4M=N,P mineral additions to litter (g N,P t-1)
C     VCO2R,VCH4R,VNOXR,VPOXR=total CO2,CH4,NOX,POX emissions from fire
C        (g C,C,N,P t-1)
C
      WTSTDG(5,NZ,NY,NX)=WTSTDG(5,NZ,NY,NX)+RWTSTDG(5,NZ,NY,NX)*FCC
      ZNH4S(0,NY,NX)=ZNH4S(0,NY,NX)+RWTSTDN(5,NZ,NY,NX)*FCOMN     
      H2PO4(0,NY,NX)=H2PO4(0,NY,NX)+RWTSTDP(5,NZ,NY,NX)*FCOMP
      COX=COX+RWTSTDG(5,NZ,NY,NX)*(1.0-FCC)*FCBOX
      CHX=CHX+RWTSTDG(5,NZ,NY,NX)*(1.0-FCC)*FCBCH
      COR=COR+RWTSTDG(5,NZ,NY,NX)*FCC
      ZOX=ZOX+RWTSTDN(5,NZ,NY,NX)*FCOGN     
      POX=POX+RWTSTDP(5,NZ,NY,NX)*FCOGP
      Z4M=Z4M+RWTSTDN(5,NZ,NY,NX)*FCOMN 
      P4M=P4M+RWTSTDP(5,NZ,NY,NX)*FCOMP
      TZIN=TZIN+Z4M
      TPIN=TPIN+P4M
      VCOXF(NZ,NY,NX)=VCOXF(NZ,NY,NX)+COR
      VCO2R(NY,NX)=VCO2R(NY,NX)-COX
      VCH4R(NY,NX)=VCH4R(NY,NX)-CHX
      VNOXR(NY,NX)=VNOXR(NY,NX)-ZOX
      VPOXR(NY,NX)=VPOXR(NY,NX)-POX
C     IF(NZ.EQ.1)THEN
C     WRITE(*,5567)'SHOOT',I,J,NFZ,NX,NY,NZ
C    2,RCGCK(NY,NX),ROGOK,RCGOK,RCHOK,RC4OK,HCBFC
C    2,HCBFCZ(NZ,NY,NX),XONET(NY,NX),XCNET(NY,NX),XHNET(NY,NX) 
C    2,FCC,FCBOX,FCBCH,FCBOX+FCBCH,RCGCKZ(NZ,NY,NX),RCGDKZ(NZ,NY,NX)
C    3,COR,COX,CHX,COR+COX+CHX,ZOX,Z4M,POX,P4M
C    3,VCO2R(NY,NX),VCH4R(NY,NX),VNOXR(NY,NX),VPOXR(NY,NX)
C    4,VCOXF(NZ,NY,NX),VNOXF(NZ,NY,NX),VPOXF(NZ,NY,NX)
C    5,RCGCK(NY,NX),RCGCKZ(NZ,NY,NX) 
C    6,OXYC(NY,NX),OXYQ(NY,NX),CH4Q(NY,NX)
C    7,FCOMN,FCOMP,TFNCOC(NZ,NY,NX)
C    8,TKQ(NY,NX),TFNCO
5567  FORMAT(A8,6I4,50F14.6)
C     ENDIF
      ENDIF
C
C     COMBUST ROOT
C
C     RCGSKR=total root C combustion of PFT from ‘grosub.f’ (g C t-1)
C     R*=C,N,P combustion rate for each PFT C pool from ‘grosub.f’ 
C        (g C,N,P t-1)
C     PFT pool code (g C,N,P):
C        CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P mass in root
C        WTRTL=total active root C
C        WTRT1,WTRT1N,WTRT1P=primary root C,N,P mass in soil layer
C        WTRT2,WTRT2N,WTRT2P=secondary root C,N,P mass in soil layer
C        CPOOLN,ZPOOLN,PPOOLN=nonstructural C,N,P in bacteria
C        WTNDL,WTNDLN,WTNDLP=bacterial C,N,P mass
C        WTRVC=storage C,N,P
C     OSC(5=soil charcoal (g C)
C     FCC=fraction of aerobic C combustion to charcoal
C     FCBOX,FCBCH=fraction of soil+root C combusted aerobically,
C        anaerobically
C     RCGOX,RCHOX=O2-limited soil+root aerobic,anaerobic combustion 
C        from ‘trnsfr.f’ (g C t-1)
C     RCGSK=O2-unlimited total soil+root combustion (g C t-1)
C     ZNH4S=soil NH4 (g N)
C     H2PO4=soil H2PO4 (g P)
C     FCOMN,FGOMN=fraction of N combustion to NH4, gaseous NOX
C     FCOMP,FCOGP=fraction of P combustion to H2PO4, gaseous POX
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C        (g C,C,N,P t-1)
C     Z4M,P4M=N,P mineral additions to litter (g N t-1)
C
      COX=0.0
      CHX=0.0
      ZOX=0.0
      POX=0.0
      COR=0.0
      Z4M=0.0
      P4M=0.0
      DO 8880 L=NU(NY,NX),NJ(NY,NX)
      IF(RCGSKR(L,NZ,NY,NX).GT.ZEROP2(NZ,NY,NX))THEN
      FCOMN=FCOMNY+(FCOMNX-FCOMNY)*(1.0-TFNCOS(L,NY,NX))
      FCOGN=1.0-FCOMN
      FCOMP=FCOMPY+(FCOMPX-FCOMPY)*(1.0-TFNCOS(L,NY,NX))
      FCOGP=1.0-FCOMP
      FCC=(RCGSK(L,NY,NX)-RCGOX(L,NY,NX)-RCHOX(L,NY,NX))/RCGSK(L,NY,NX)
      FCBOXR=RCGOX(L,NY,NX)/(RCGOX(L,NY,NX)+RCHOX(L,NY,NX))
      FCBCHR=RCHOX(L,NY,NX)/(RCGOX(L,NY,NX)+RCHOX(L,NY,NX))
      IF(L.EQ.NU(NY,NX))THEN
      OSC(5,1,L,NY,NX)=OSC(5,1,L,NY,NX)+RWTRVC(NZ,NY,NX)*FCC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RWTRVN(NZ,NY,NX)*FCOMN     
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RWTRVP(NZ,NY,NX)*FCOMP
      COX=COX+RWTRVC(NZ,NY,NX)*(1.0-FCC)*FCBOXR
      CHX=CHX+RWTRVC(NZ,NY,NX)*(1.0-FCC)*FCBCHR
      COR=COR+RWTRVC(NZ,NY,NX)*FCC 
      ZOX=ZOX+RWTRVN(NZ,NY,NX)*FCOGN 
      POX=POX+RWTRVP(NZ,NY,NX)*FCOGP
      Z4M=Z4M+RWTRVN(NZ,NY,NX)*FCOMN 
      P4M=P4M+RWTRVP(NZ,NY,NX)*FCOMP
      ENDIF
      DO 8885 N=1,MY(NZ,NY,NX)
      OSC(5,1,L,NY,NX)=OSC(5,1,L,NY,NX)+RCPOOLR(N,L,NZ,NY,NX)*FCC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RZPOOLR(N,L,NZ,NY,NX)*FCOMN     
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RPPOOLR(N,L,NZ,NY,NX)*FCOMP
      COX=COX+RCPOOLR(N,L,NZ,NY,NX)*(1.0-FCC)*FCBOXR
      CHX=CHX+RCPOOLR(N,L,NZ,NY,NX)*(1.0-FCC)*FCBCHR
      COR=COR+RCPOOLR(N,L,NZ,NY,NX)*FCC 
      ZOX=ZOX+RZPOOLR(N,L,NZ,NY,NX)*FCOGN     
      POX=POX+RPPOOLR(N,L,NZ,NY,NX)*FCOGP
      Z4M=Z4M+RZPOOLR(N,L,NZ,NY,NX)*FCOMN 
      P4M=P4M+RPPOOLR(N,L,NZ,NY,NX)*FCOMP
      DO 8785 NR=1,NRT(NZ,NY,NX)
      RCMBC=RWTRT1(N,L,NR,NZ,NY,NX)+RWTRT2(N,L,NR,NZ,NY,NX)
      RCMBN=RWTRT1N(N,L,NR,NZ,NY,NX)+RWTRT2N(N,L,NR,NZ,NY,NX)
      RCMBP=RWTRT1P(N,L,NR,NZ,NY,NX)+RWTRT2P(N,L,NR,NZ,NY,NX)
      OSC(5,1,L,NY,NX)=OSC(5,1,L,NY,NX)+RCMBC*FCC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCMBN*FCOMN     
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCMBP*FCOMP
      COX=COX+RCMBC*(1.0-FCC)*FCBOXR
      CHX=CHX+RCMBC*(1.0-FCC)*FCBCHR
      COR=COR+RCMBC*FCC 
      ZOX=ZOX+RCMBN*FCOGN     
      POX=POX+RCMBP*FCOGP
      Z4M=Z4M+RCMBN*FCOMN 
      P4M=P4M+RCMBP*FCOMP
8785  CONTINUE
8885  CONTINUE
      RCMBC=RCPOOLN(L,NZ,NY,NX)+RWTNDL(L,NZ,NY,NX)
      RCMBN=RZPOOLN(L,NZ,NY,NX)+RWTNDLN(L,NZ,NY,NX)
      RCMBP=RPPOOLN(L,NZ,NY,NX)+RWTNDLP(L,NZ,NY,NX)
      OSC(5,1,L,NY,NX)=OSC(5,1,L,NY,NX)+RCMBC*FCC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCMBN*FCOMN     
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCMBP*FCOMP
      COX=COX+RCMBC*(1.0-FCC)*FCBOXR
      CHX=CHX+RCMBC*(1.0-FCC)*FCBCHR
      COR=COR+RCMBC*FCC
      ZOX=ZOX+RCMBN*FCOGN     
      POX=POX+RCMBP*FCOGP
      Z4M=Z4M+RCMBN*FCOMN 
      P4M=P4M+RCMBP*FCOMP
C     IF(NZ.EQ.3)THEN
C     WRITE(*,5569)'ROOT',I,J,NFZ,NX,NY,NZ,L
C    3,RCGSKR(L,NZ,NY,NX),RCGSK(L,NY,NX),RCGOX(L,NY,NX)
C    4,RCHOX(L,NY,NX),RCGOX(L,NY,NX)+RCHOX(L,NY,NX)
C    5,FCC,FCBOXR,FCBCHR,FCBOXR+FCBCHR
C    2,COR,COX,CHX,COR+COX+CHX,ZOX,Z4M,POX,P4M
C    3,VCO2R(NY,NX),VCH4R(NY,NX),VNOXR(NY,NX),VPOXR(NY,NX)
C    4,VCOXF(NZ,NY,NX),VNOXF(NZ,NY,NX),VPOXF(NZ,NY,NX)
C    3,FCBOXR*FCOGN,FCBOXR*FCOGP 
C    7,FCOMN,FCOMP,TFNCOS(L,NY,NX)
C    8,FCBOXR*FCOGP
5569  FORMAT(A8,7I4,40F14.6)
C     ENDIF
C
C     SALT FROM ROOTS
C
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C     R*R=salt combustion rate for each PFT root pool (mol t-1)
C     Z*=salt in soil (mol)
C     salt code:AL=Al,FE=Fe,CA=Ca,MG=Mg,NA=Na,KA=K,SO=SO4,CL=Cl
C
      IF(ISALTG.NE.0)THEN
      DO 8886 N=1,MY(NZ,NY,NX)
      ZAL(L,NY,NX)=ZAL(L,NY,NX)+RZALR(N,L,NZ,NY,NX)
      ZFE(L,NY,NX)=ZFE(L,NY,NX)+RZFER(N,L,NZ,NY,NX)
      ZCA(L,NY,NX)=ZCA(L,NY,NX)+RZCAR(N,L,NZ,NY,NX)
      ZMG(L,NY,NX)=ZMG(L,NY,NX)+RZMGR(N,L,NZ,NY,NX)
      ZNA(L,NY,NX)=ZNA(L,NY,NX)+RZNAR(N,L,NZ,NY,NX)
      ZKA(L,NY,NX)=ZKA(L,NY,NX)+RZKAR(N,L,NZ,NY,NX)
      ZSO4(L,NY,NX)=ZSO4(L,NY,NX)+RZSOR(N,L,NZ,NY,NX)
      ZCL(L,NY,NX)=ZCL(L,NY,NX)+RZCLR(N,L,NZ,NY,NX)
      TIONIN=TIONIN+RZALR(N,L,NZ,NY,NX)+RZFER(N,L,NZ,NY,NX)
     2+RZCAR(N,L,NZ,NY,NX)+RZMGR(N,L,NZ,NY,NX)+RZNAR(N,L,NZ,NY,NX)
     3+RZKAR(N,L,NZ,NY,NX)+RZSOR(N,L,NZ,NY,NX)+RZCLR(N,L,NZ,NY,NX)
C     WRITE(*,5572)'RZR',I,J,NFZ,NX,NY,NZ,N,L
C    2,TIONIN,RZALR(N,L,NZ,NY,NX),RZFER(N,L,NZ,NY,NX)
C    2,RZCAR(N,L,NZ,NY,NX),RZMGR(N,L,NZ,NY,NX),RZNAR(N,L,NZ,NY,NX)
C    3,RZKAR(N,L,NZ,NY,NX),RZSOR(N,L,NZ,NY,NX),RZCLR(N,L,NZ,NY,NX)
5572  FORMAT(A8,8I4,50F16.8)
8886  CONTINUE
      ENDIF
      ENDIF
8880  CONTINUE
C
C     BOUNDARY FLUXES FROM ROOTS TO SOIL
C
      CO2GIN=CO2GIN+COX+CHX
      ZCSNC(NY,NX)=ZCSNC(NY,NX)+COR
      TZIN=TZIN+Z4M
      TPIN=TPIN+P4M
      VCO2R(NY,NX)=VCO2R(NY,NX)-COX
      VCH4R(NY,NX)=VCH4R(NY,NX)-CHX
      VNOXR(NY,NX)=VNOXR(NY,NX)-ZOX
      VPOXR(NY,NX)=VPOXR(NY,NX)-POX
9970  CONTINUE
C     WRITE(*,1192)'VCO2R',I,J,NFZ,NX,NY,VCO2R(NY,NX),VCH4R(NY,NX)
C    2,COX,CHX,COR,RCGSK(L,NY,NX),RCGOX(L,NY,NX),RCHOX(L,NY,NX)
1192  FORMAT(A8,5I4,12F16.8)
      ENDIF
C
C     END FIRE
C
      DO 9975 NZ=1,NP(NY,NX)
C
C     DEAD AND LIVING TRANSFERS
C
C     ARSTT,ARSTD=total,PFT standing dead surface area (m2)
C     TRN=all canopy+standing dead net SW+LW radiation (MJ h-1)
C     TLE=all canopy+standing dead latent heat flux (MJ h-1)
C     TSH=all canopy+standing dead sensible heat flux (MJ h-1)
C     TGH=all canopy+standing dead storage heat flux (MJ h-1) 
C     TRNP=PFT canopy+standing dead net SW+LW radiation (MJ h-1)
C     TLEP=PFT canopy+standing dead latent heat flux (MJ h-1)
C     TSHP=PFT canopy+standing dead sensible heat flux (MJ h-1)
C     TGHP=PFT canopy+standing dead storage heat flux (MJ h-1) 
C     TTRAN=all canopy+standing dead surface transpiration,evaporation
C        (m3 H-1)
C     CTRAN=cumulative canopy+standing dead surface transpiration,
C        evaporation (m3)
C     TVOLWP=all water volume in canopy+standing dead (m3)
C     TVOLWC=all water volume on canopy+standing dead surfaces (m3)
C     TEVAPP=all canopy+standing dead transpiration+evaporation(m3 t-1)
C     TEVAPC=all canopy+standing dead evaporation (m3 t-1)
C     RFLXC=PFT canopy+standing dead net SW+LW radiation (MJ t-1)
C     EFLXC=PFT canopy+standing dead latent heat flux (MJ t-1)
C     SFLXC=PFT canopy+standing dead sensible heat flux(MJ t-1)
C     HFLXC=PFT canopy+standing dead storage heat flux (MJ t-1) 
C     VFLXC=convective heat flux from EFLXC (MJ t-1)
C     EP,EVAPC=PFT canopy+standing dead surface transpiration,
C        evaporation (m3 t-1)
C     VOLWP,VOLWC=PFT water volume in canopy,on canopy surfaces (m3)
C     VOLWQ=PFT water volume on standing dead surfaces (m3)
C     FLWC,FLWD=PFT foliar,standing dead water retention of
C        precipitation (m3 t-1)
C     TKAM,TKC=air,canopy surface temperature (K)
C
      DO 900 L=1,JC
      ARSDT(L,NY,NX)=ARSDT(L,NY,NX)+ARSTD(L,NZ,NY,NX)
900   CONTINUE
      ARSDC(NY,NX)=ARSDC(NY,NX)+ARSTG(NZ,NY,NX)
      TRN(NY,NX)=TRN(NY,NX)+RFLXC(NZ,NY,NX)
      TLE(NY,NX)=TLE(NY,NX)+EFLXC(NZ,NY,NX)
      TSH(NY,NX)=TSH(NY,NX)+SFLXC(NZ,NY,NX)
      TGH(NY,NX)=TGH(NY,NX)-(HFLXC(NZ,NY,NX)-VFLXC(NZ,NY,NX))
      TRNP(NZ,NY,NX)=TRNP(NZ,NY,NX)+RFLXC(NZ,NY,NX)
      TLEP(NZ,NY,NX)=TLEP(NZ,NY,NX)+EFLXC(NZ,NY,NX)
      TSHP(NZ,NY,NX)=TSHP(NZ,NY,NX)+SFLXC(NZ,NY,NX)
      TGHP(NZ,NY,NX)=TGHP(NZ,NY,NX)-(HFLXC(NZ,NY,NX)-VFLXC(NZ,NY,NX))
      TTRAN(NZ,NY,NX)=TTRAN(NZ,NY,NX)+EP(NZ,NY,NX)+EVAPC(NZ,NY,NX)
      CTRAN(NZ,NY,NX)=CTRAN(NZ,NY,NX)+EP(NZ,NY,NX)+EVAPC(NZ,NY,NX)
      TVOLWP(NY,NX)=TVOLWP(NY,NX)+VOLWP(NZ,NY,NX)
      TVOLWC(NY,NX)=TVOLWC(NY,NX)+VOLWC(NZ,NY,NX)+VOLWQ(NZ,NY,NX)
      TEVAPP(NY,NX)=TEVAPP(NY,NX)+EP(NZ,NY,NX)+EVAPC(NZ,NY,NX)
      TEVAPC(NY,NX)=TEVAPC(NY,NX)+EVAPC(NZ,NY,NX)
      ENGYC=4.19*(VOLWC(NZ,NY,NX)+VOLWQ(NZ,NY,NX)
     2+FLWC(NZ,NY,NX)+FLWD(NZ,NY,NX)+EVAPC(NZ,NY,NX))
     2*TKC(NZ,NY,NX)
      TENGYC(NY,NX)=TENGYC(NY,NX)+ENGYC
      THFLXC(NY,NX)=THFLXC(NY,NX)+ENGYC-ENGYX(NZ,NY,NX)
     2-(FLWC(NZ,NY,NX)+FLWD(NZ,NY,NX))*4.19*TKAM(NY,NX)*XNFH
      ENGYX(NZ,NY,NX)=ENGYC
C     IF(I.EQ.227.AND.NX.EQ.1.AND.NZ.EQ.1)THEN
C     WRITE(*,5570)'TEVAPC',I,J,NFZ,NX,NY,NZ,IFLGC(NZ,NY,NX)
C    2,TRNP(NZ,NY,NX),RFLXC(NZ,NY,NX)
C    3,TVOLWP(NY,NX),VOLWP(NZ,NY,NX)
C    2,TEVAPC(NY,NX),TEVAPP(NY,NX),EP(NZ,NY,NX),EVAPC(NZ,NY,NX)
C    3,(UPWTR(1,L,NZ,NY,NX),L=NU(NY,NX),NI(NZ,NY,NX))
5570  FORMAT(A8,7I4,20F16.8) 
C     ENDIF
C
C     LIVING TRANSFERS
C
C     IFLGC=PFT flag:0=not active,1=active
C
      IF(IFLGC(NZ,NY,NX).EQ.1)THEN
C
C     TOTAL LEAF AREA OF ALL PLANT SPECIES
C
C     ARLFT,ARSTT=total leaf,stalk area of combined canopy layer (m2)
C     ARLFV,ARSTV=PFT leaf,stalk area in canopy layer (m2)
C     WGLFT=total leaf C of combined canopy layer (g C)
C     WGLFV=PFT leaf C in canopy layer (g C)
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
C     RTDNT=total root length density (m m-3) 
C     RTDNP=PFT root length density per plant (m m-3)
C     TUPWTR=total root water uptake (m3 t-1)
C     UPWTR=PFT root water uptake (m3 t-1)
C     TUPHT=total convective heat in root water uptake (MJ t-1)
C     TKS=soil temperature (K)
C     PP=PFT population (n)
C     AREA=grid cell atea (m2)
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
C        gas code (g C,O,C,N,N,H t-1):
C        CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
C     R*FLA=root gaseous-atmosphere CO2 exchange (g C t-1)
C     R*DFA=root aqueous-gaseous CO2 exchange (g C t-1)
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
C        gas code (g C,O,C,N,N,H):
C        CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
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
C     IF(J.EQ.15.AND.NFZ.EQ.NFH.AND.NZ.EQ.2)THEN
C     WRITE(*,5566)'TLOXYP',I,J,NFZ,NX,NY,NZ,L,N
C    2,TLOXYP(L,NY,NX),OXYP(N,L,NZ,NY,NX),OXYA(N,L,NZ,NY,NX) 
C    3,ROXDFA(N,L,NZ,NY,NX),RUPOXP(N,L,NZ,NY,NX),ROXFLA(N,L,NZ,NY,NX)
C    4,TLCO2P(L,NY,NX),CO2P(N,L,NZ,NY,NX),CO2A(N,L,NZ,NY,NX)
C    5,RCODFA(N,L,NZ,NY,NX),RCO2P(N,L,NZ,NY,NX)
C    6,RCO2S(N,L,NZ,NY,NX)
C    4,TLCH4P(L,NY,NX),CH4P(N,L,NZ,NY,NX),CH4A(N,L,NZ,NY,NX)
C    5,RCHDFA(N,L,NZ,NY,NX),RUPCHS(N,L,NZ,NY,NX),RCHFLA(N,L,NZ,NY,NX) 
5566  FORMAT(A8,8I4,20E12.4)
C     ENDIF
C
C     TOTAL ROOT BOUNDARY GAS AND SOLUTE FLUXES
C
C     T*FLA=total root gaseous-atmosphere gas exchange
C     R*FLA=PFT root gaseous-atmosphere gas exchange
C        gas code (g C,O,C,N,N,H t-1):
C        CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
C     TUP*S,TUP*B=total root-soil gas, solute exchange in 
C        non-band,band
C     RUP*S,RUP*B*=PFT root-soil gas, solute exchange in non-band,band
C        solute code (g N,P t-1):
C        NH4=NH4,NO3=NO3,H2P=H2PO4,H1P=H1PO4 in non-band
C        NHB=NH4,NOB=NO3,H2B=H2PO4,H1B=H1PO4 in band
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
C     WRITE(*,4141)'TCO2S',I,J,NFZ,NX,NY,L,NZ,N,TCO2S(L,NY,NX)
C    2,RCO2S(N,L,NZ,NY,NX),TCO2P(L,NY,NX),RCO2P(N,L,NZ,NY,NX)
C     WRITE(*,4141)'TUPOX',I,J,NFZ,NX,NY,L,NZ,N,TUPOXS(L,NY,NX)
C    2,RUPOXS(N,L,NZ,NY,NX),TUPOXP(L,NY,NX),RUPOXP(N,L,NZ,NY,NX)
C     WRITE(*,4141)'TCHFLA',I,J,NFZ,NX,NY,L,NZ,N,TCHFLA(L,NY,NX)
C    2,RCHFLA(N,L,NZ,NY,NX),RUPCHS(N,L,NZ,NY,NX)
C     WRITE(*,4141)'TUPN2S',I,J,NFZ,NX,NY,L,NZ,N,TUPN2S(L,NY,NX)
C    2,RUPN2S(N,L,NZ,NY,NX) 
4141  FORMAT(A8,8I4,12F16.8)
C     ENDIF
C
C     SALT UPTAKE
C
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C     RUP*=root salt uptake by each PFT from ‘uptake.f’ (mol t-1)
C     TUP*=total root salt uptake (mol t-1)
C
      IF(ISALTG.NE.0)THEN
      TUPZAL(L,NY,NX)=TUPZAL(L,NY,NX)+RUPZAL(N,L,NZ,NY,NX)
      TUPZFE(L,NY,NX)=TUPZFE(L,NY,NX)+RUPZFE(N,L,NZ,NY,NX)
      TUPZCA(L,NY,NX)=TUPZCA(L,NY,NX)+RUPZCA(N,L,NZ,NY,NX)
      TUPZMG(L,NY,NX)=TUPZMG(L,NY,NX)+RUPZMG(N,L,NZ,NY,NX)
      TUPZNA(L,NY,NX)=TUPZNA(L,NY,NX)+RUPZNA(N,L,NZ,NY,NX)
      TUPZKA(L,NY,NX)=TUPZKA(L,NY,NX)+RUPZKA(N,L,NZ,NY,NX)
      TUPZSO(L,NY,NX)=TUPZSO(L,NY,NX)+RUPZSO(N,L,NZ,NY,NX)
      TUPZCL(L,NY,NX)=TUPZCL(L,NY,NX)+RUPZCL(N,L,NZ,NY,NX)
C     IF(L.LE.3)THEN
C     WRITE(*,4141)'TUPZAL',I,J,NFZ,NX,NY,NZ,L,N,TUPZAL(L,NY,NX)
C    2,RUPZAL(N,L,NZ,NY,NX)
C     WRITE(*,4141)'TUPZSO',I,J,NFZ,NX,NY,NZ,L,N,TUPZSO(L,NY,NX)
C    2,RUPZSO(N,L,NZ,NY,NX)
C     ENDIF
      ENDIF
C
C     TOTAL ROOT C,N,P EXUDATION
C
C     TDFOMC,TDFOMN,TDFOMP=total root-soil nonstructural C,N,P
C        exchange (g C,N,P t-1)
C     RDFOMC,RDFOMN,RDFOMP=PFT root-soil nonstructural C,N,P 
C        exchange (g C,N,P t-1) 
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
C     ROXYX=O2 demand by all microbial,root,mycorrhizal populations
C        (g O t-1) 
C     RNH4X=NH4 demand in non-band by all microbial,root,mycorrhizal
C        populations (g N t-1)
C     RNO3X=NO3 demand in non-band by all microbial,root,mycorrhizal
C        populations (g N t-1)
C     RPO4X=H2PO4 demand in non-band by all microbial,root,mycorrhizal
C        populations (g P t-1)
C     RP14X=HPO4 demand in non-band by all microbial,root,mycorrhizal
C        populations (g P t-1)
C     RNHBX=NH4 demand in band by all microbial,root,mycorrhizal
C        populations (g N t-1)
C     RN3BX=NO3 demand in band by all microbial,root,mycorrhizal
C        populations (g N t-1)
C     RPOBX=H2PO4 demand in band by all microbial,root,mycorrhizal 
C        populations (g P t-1)
C     RP1BX=HPO4 demand in band by all microbial,root,mycorrhizal
C        populations (g P t-1)
C     ROXYP=O2 demand by each root,mycorrhizal population (g O t-1)
C     RUNNHP=NH4 demand in non-band by each root,mycorrhizal population
C        (g N t-1)
C     RUNNOP=NO3 demand in non-band by each root,mycorrhizal population
C        (g N t-1)
C     RUPP2P=H2PO4 demand in non-band by each root,mycorrhizal
C        population (g P t-1)
C     RUPP1P=HPO4 demand in non-band by each root,mycorrhizal
C        population (g P t-1)
C     RUNNBP=NH4 demand in band by each root,mycorrhizal population
C        (g N t-1)
C     RUNNXB=NO3 demand in band by each root,mycorrhizal population
C        (g N t-1)
C     RUPP2B=H2PO4 demand in band by each root,mycorrhizal population
C        (g P t-1)
C     RUPP1B=HPO4 demand in band by each root,mycorrhizal population
C        (g P t-1)
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
C
C     TOTAL ROOT N2 FIXATION BY ALL PLANT SPECIES
C
C     TUPNF=total root N2 fixation (g N t-1)
C     RUPNF=PFT root N2 fixation (g N t-1)
C
      DO 85 L=NU(NY,NX),NI(NZ,NY,NX)
      TUPNF(L,NY,NX)=TUPNF(L,NY,NX)+RUPNF(L,NZ,NY,NX)
85    CONTINUE
C
C     TOTAL WATER, CO2 FLUXES
C
C     TCCAN=total net CO2 fixation (g C h-1)
C     HCNET=PFT net CO2 fixation (g C h-1) 	
C     XCNET,XONET=total CO2,O2 exchange by all PFT canopies
C        (g C,C,O t-1)
C     CNET=PFT net CO2 fixation (g C t-1)
C     ARLFC,ARSTC=total leaf,stalk area (m2)
C     ARLFP,ARSTP=PFT leaf,stalk area (m2)
C     ZCSNC,ZZSNC,ZPSNC=total net root-soil C,N,P exchange 
C        (g C,N,P t-1) 
C     HCUPTK,HZUPTK,HPUPTK=PFT net root-soil C,N,P exchange 
C        (g C,N,P t-1) 
C     TBALC,TBALN,TBALP=total C,N,P balance (g C,N,P)
C     BALC,BALN,BALP=PFT C,N,P balance (g C,N,P)
C     TCO2Z,TOXYZ,TCH4Z,TN2OZ,TNH3Z,TH2GZ=total loss of root CO2, O2,
C        CH4, N2O, NH3, H2 (g C,O,C,N,N,H t-1)
C     RCO2Z,ROXYZ,RCH4Z,RN2OZ,RNH3Z,RH2GZ=PFT loss of root CO2, O2,
C        CH4, N2O, NH3, H2 (g C,O,C,N,N,H t-1)
C
      IF(NFZ.EQ.NFH)THEN
      TCCAN(NY,NX)=TCCAN(NY,NX)+HCNET(NZ,NY,NX)
      ENDIF
      XCNET(NY,NX)=XCNET(NY,NX)+CNET(NZ,NY,NX)
      XONET(NY,NX)=XONET(NY,NX)-CNET(NZ,NY,NX)*2.667
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
C        (g N t-1)
C     TNH3C=total NH3 flux between atmosphere and canopy (g N t-1)
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
