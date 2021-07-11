      SUBROUTINE hour1(I,J,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE REINITIALIZES HOURLY VARIABLES USED IN OTHER
C     SUBROUTINES
C
      include "parameters.h"
      include "files.h"
      include "filec.h"
      include "blkc.h"
      include "blk1cp.h"
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
      include "blk22b.h"
      include "blk22c.h"
      CHARACTER*16 DATA(30)
      DIMENSION OFC(2),OFN(2),OFP(2),CNOF(4),CPOF(4)
      DIMENSION PSISK(0:100),THETK(100),TAUY(0:JC+1),RABSL(0:JC+1)
     2,RABPL(0:JC+1),RAFSL(0:JC+1),RAFPL(0:JC+1),RADSL(JP,JY,JX)
     3,RADPL(JP,JY,JX),RADS1(JP,JY,JX),RADP1(JP,JY,JX),RADS2(JP,JY,JX)
     4,RADP2(JP,JY,JX),RAYSL(JP,JY,JX),RAYPL(JP,JY,JX),RAYS1(JP,JY,JX)
     5,RAYP1(JP,JY,JX),RAYS2(JP,JY,JX),RAYP2(JP,JY,JX),RADSW(JP,JY,JX)
     6,RADPW(JP,JY,JX),RADW1(JP,JY,JX),RADQ1(JP,JY,JX),RADW2(JP,JY,JX)
     7,RADQ2(JP,JY,JX),RAYSW(JP,JY,JX),RAYPW(JP,JY,JX),RAYW1(JP,JY,JX)
     8,RAYQ1(JP,JY,JX),RAYW2(JP,JY,JX),RAYQ2(JP,JY,JX),RADSA(JP,JY,JX)
     9,RAPSA(JP,JY,JX),RDNDIR(4,4,JP,JY,JX),PARDIR(4,4,JP,JY,JX)
     1,RADWA(JP,JY,JX),RAPWA(JP,JY,JX),RDNDIW(4,4,JP,JY,JX)
     2,PARDIW(4,4,JP,JY,JX),TSURF(4,JZ,JP,JY,JX),TSURFB(4,JZ,JP,JY,JX)
     3,XVOLWC(0:3),ZL1(0:JZ,JY,JX),THETPZ(JZ,JY,JX),TRADC(JY,JX)
     4,TRAPC(JY,JX),TRADG(JY,JX),TRAPG(JY,JX),THETRX(0:2),DPTH0(JY,JX)
     5,IALBS(4,4)
C
C     *SG=diffusivity (m2 h-1):CG=CO2g,CL=CO2s,CH=CH4g,CQ=CH4s,OG=O2g
C     OL=O2s,ZG=N2g,ZL=N2s,Z2=N2Og,ZV=N2Os,ZH=NH3g,ZN=NH3s,ZO=NO3
C     PO=H2PO4,OC=DOC,ON=DON,OP=DOP,OA=acetate,WG=H2Og,AL=Al,FE=Fe
C     HY=H,CA=Ca,GM=Mg,AN=Na,AK=K,OH=OH,C3=CO3,HC=HCO3,SO=SO4,CL=Cl
C     HG=H2g,HL=H2s
C
      PARAMETER (CGSG=4.68E-02,CLSG=4.25E-06,CHSG=7.80E-02
     2,CQSG=7.08E-06,OGSG=6.43E-02,OLSG=8.57E-06,ZGSG=5.57E-02
     3,ZLSG=7.34E-06,Z2SG=5.57E-02,ZVSG=5.72E-06,ZHSG=6.67E-02
     4,ZNSG=4.00E-06,ZOSG=6.00E-06,POSG=3.00E-06,OCSG=1.0E-08
     5,ONSG=1.0E-08,OPSG=1.0E-08,OASG=3.64E-06,WGSG=7.70E-02
     6,ALSG=5.0E-06,FESG=5.0E-06,HYSG=5.0E-06,CASG=5.0E-06
     7,GMSG=5.0E-06,ANSG=5.0E-06,AKSG=5.0E-06,OHSG=5.0E-06
     8,C3SG=5.0E-06,HCSG=5.0E-06,SOSG=5.0E-06,CLSX=5.0E-06
     9,HGSG=5.57E-02,HLSG=7.34E-06)
C
C     SC*X=solubility (g m-3/g m-3):CO2=CO2,CH4=CH4,OXY=O2,N2G=N2g
C     N2O=N2O,NH3=NH3,H2G=H2
C     AC*X=activity (g m-3):CO2=CO2,CH4=CH4,OXY=O2,N2G=N2g
C     N2O=N2O,NH3=NH3,H2G=H2
C
      PARAMETER (SCO2X=7.391E-01,SCH4X=3.156E-02,SOXYX=2.925E-02
     2,SN2GX=1.510E-02,SN2OX=5.241E-01,SNH3X=2.852E+02,SH2GX=3.156E-02
     3,ACO2X=0.14,ACH4X=0.14,AOXYX=0.31,AN2GX=0.23,AN2OX=0.23
     4,ANH3X=0.07,AH2GX=0.14)
C
C     ALBRW,ALBPW=stalk albedo for shortwave,PAR
C     VISCW=water viscosity (Mg m-1 s)
C     BKDSX=maximm soil bulk density, ZW=snowpack surface roughness (m)
C     CFW=stalk clumping factor,FORGW=minimum SOC or organic soil (g Mg-1)
C     THETPW=minimum air-filled porosity for saturation (m3 m-3)
C     RAM=minimum boundary layer resistance (h m-1)
C
      PARAMETER (ALBRW=0.1,ALBPW=0.1,ABSRW=1.0-ALBRW,ABSPW=1.0-ALBPW)
      PARAMETER (VISCW=1.0E-06,BKDSX=1.89,ZW=0.01,CFW=0.5
     2,FORGW=0.25E+06,DTHETW=1.0E-06,THETPW=0.01,THETWP=1.0-THETPW 
     3,RAM=2.78E-03)
C
C     XVOLWC=foliar water retention capacity (m3 m-2)
C     THETRX=litter water retention capacity (m3 g C-1) 
C
      DATA XVOLWC/5.0E-04,2.5E-04,2.5E-04,2.5E-04/
      DATA THETRX/4.0E-06,8.0E-06,8.0E-06/
      REAL*4 TFACL,TFACG,TFACW,TFACR,TFACA
      XJ=J
      DOY=I-1+XJ/24
C
C     RESET HOURLY SOIL ACCUMULATORS FOR WATER, HEAT, GASES, SOLUTES
C
      VOLWSO=0.0
      HEATSO=0.0
      OXYGSO=0.0
      TLH2G=0.0
      TSEDSO=0.0
      TLRSDC=0.0
      TLORGC=0.0
      TLCO2G=0.0
      TLRSDN=0.0
      TLORGN=0.0
      TLN2G=0.0
      TLRSDP=0.0
      TLORGP=0.0
      TLNH4=0.0
      TLNO3=0.0
      TLPO4=0.0
      TION=0.0
      TBALC=0.0
      TBALN=0.0
      TBALP=0.0
C
C     CONCENTRATIONS OF CO2, CH4, O2, N2, N2O, NH3, H2 IN ATMOSPHERE,
C     PRECIPITATION AND IRRIGATION FROM MIXING RATIOS READ IN 'READS'
C
C     C*E,C*R,C*Q=atmospheric,precipitation,irrigation solute concentrations
C     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
C             :*ZN3*=NH3,*H2G*=H2
C
      DO 9145 NX=NHW,NHE
      DO 9140 NY=NVN,NVS
      CCO2E(NY,NX)=CO2E(NY,NX)*5.36E-04*273.15/TKA(NY,NX)
      CCH4E(NY,NX)=CH4E(NY,NX)*5.36E-04*273.15/TKA(NY,NX)
      COXYE(NY,NX)=OXYE(NY,NX)*1.43E-03*273.15/TKA(NY,NX)
      CZ2GE(NY,NX)=Z2GE(NY,NX)*1.25E-03*273.15/TKA(NY,NX)
      CZ2OE(NY,NX)=Z2OE(NY,NX)*1.25E-03*273.15/TKA(NY,NX)
      CNH3E(NY,NX)=ZNH3E(NY,NX)*6.25E-04*273.15/TKA(NY,NX)
      CH2GE(NY,NX)=H2GE(NY,NX)*8.92E-05*273.15/TKA(NY,NX)
      CCOR(NY,NX)=CCO2E(NY,NX)*SCO2X/(EXP(ACO2X*CSTRR(NY,NX)))
     2*EXP(0.843-0.0281*TCA(NY,NX))
      CCHR(NY,NX)=CCH4E(NY,NX)*SCH4X/(EXP(ACH4X*CSTRR(NY,NX)))
     2*EXP(0.597-0.0199*TCA(NY,NX))
      COXR(NY,NX)=COXYE(NY,NX)*SOXYX/(EXP(AOXYX*CSTRR(NY,NX)))
     2*EXP(0.516-0.0172*TCA(NY,NX))
      CNNR(NY,NX)=CZ2GE(NY,NX)*SN2GX/(EXP(AN2GX*CSTRR(NY,NX)))
     2*EXP(0.456-0.0152*TCA(NY,NX))
      CN2R(NY,NX)=CZ2OE(NY,NX)*SN2OX/(EXP(AN2OX*CSTRR(NY,NX)))
     2*EXP(0.897-0.0299*TCA(NY,NX))
      CCOQ(NY,NX)=CCO2E(NY,NX)*SCO2X/(EXP(ACO2X*CSTRQ(I,NY,NX)))
     2*EXP(0.843-0.0281*TCA(NY,NX))
      CCHQ(NY,NX)=CCH4E(NY,NX)*SCH4X/(EXP(ACH4X*CSTRQ(I,NY,NX)))
     2*EXP(0.597-0.0199*TCA(NY,NX))
      COXQ(NY,NX)=COXYE(NY,NX)*SOXYX/(EXP(AOXYX*CSTRQ(I,NY,NX)))
     2*EXP(0.516-0.0172*TCA(NY,NX))
      CNNQ(NY,NX)=CZ2GE(NY,NX)*SN2GX/(EXP(AN2GX*CSTRQ(I,NY,NX)))
     2*EXP(0.456-0.0152*TCA(NY,NX))
      CN2Q(NY,NX)=CZ2OE(NY,NX)*SN2OX/(EXP(AN2OX*CSTRQ(I,NY,NX)))
     2*EXP(0.897-0.0299*TCA(NY,NX))
9140  CONTINUE
9145  CONTINUE
C
C     RESET FLUX ARRAYS USED IN OTHER SUBROUTINES
C
      DO 9895 NX=NHW,NHE+1
      DO 9890 NY=NVN,NVS+1
      DO 9885 L=0,NL(NY,NX)+1
      DO 9880 N=1,3
C
C     WATER,SNOW,SOLUTE RUNOFF
C
      IF(L.EQ.0.AND.N.NE.3)THEN
      DO 9835 NN=1,2
      QR(N,NN,NY,NX)=0.0
      HQR(N,NN,NY,NX)=0.0
      DO 9870 K=0,4
      XOCQRS(K,N,NN,NY,NX)=0.0
      XONQRS(K,N,NN,NY,NX)=0.0
      XOPQRS(K,N,NN,NY,NX)=0.0
      XOAQRS(K,N,NN,NY,NX)=0.0
9870  CONTINUE
      XCOQRS(N,NN,NY,NX)=0.0
      XCHQRS(N,NN,NY,NX)=0.0
      XOXQRS(N,NN,NY,NX)=0.0
      XNGQRS(N,NN,NY,NX)=0.0
      XN2QRS(N,NN,NY,NX)=0.0
      XHGQRS(N,NN,NY,NX)=0.0
      XN4QRW(N,NN,NY,NX)=0.0
      XN3QRW(N,NN,NY,NX)=0.0
      XNOQRW(N,NN,NY,NX)=0.0
      XNXQRS(N,NN,NY,NX)=0.0
      XP1QRW(N,NN,NY,NX)=0.0
      XP4QRW(N,NN,NY,NX)=0.0
9835  CONTINUE
      QS(N,NY,NX)=0.0
      QW(N,NY,NX)=0.0
      QI(N,NY,NX)=0.0
      HQS(N,NY,NX)=0.0
      XCOQSS(N,NY,NX)=0.0
      XCHQSS(N,NY,NX)=0.0
      XOXQSS(N,NY,NX)=0.0
      XNGQSS(N,NY,NX)=0.0
      XN2QSS(N,NY,NX)=0.0
      XN4QSS(N,NY,NX)=0.0
      XN3QSS(N,NY,NX)=0.0
      XNOQSS(N,NY,NX)=0.0
      XP1QSS(N,NY,NX)=0.0
      XP4QSS(N,NY,NX)=0.0
C
C     IF EROSION FLAG SET
C
      IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
      DO 9855 NN=1,2
      XSEDER(N,NN,NY,NX)=0.0
      XSANER(N,NN,NY,NX)=0.0
      XSILER(N,NN,NY,NX)=0.0
      XCLAER(N,NN,NY,NX)=0.0
      XCECER(N,NN,NY,NX)=0.0
      XAECER(N,NN,NY,NX)=0.0
      XNH4ER(N,NN,NY,NX)=0.0
      XNH3ER(N,NN,NY,NX)=0.0
      XNHUER(N,NN,NY,NX)=0.0
      XNO3ER(N,NN,NY,NX)=0.0
      XNH4EB(N,NN,NY,NX)=0.0
      XNH3EB(N,NN,NY,NX)=0.0
      XNHUEB(N,NN,NY,NX)=0.0
      XNO3EB(N,NN,NY,NX)=0.0
      XN4ER(N,NN,NY,NX)=0.0
      XNBER(N,NN,NY,NX)=0.0
      XHYER(N,NN,NY,NX)=0.0
      XALER(N,NN,NY,NX)=0.0
      XCAER(N,NN,NY,NX)=0.0
      XMGER(N,NN,NY,NX)=0.0
      XNAER(N,NN,NY,NX)=0.0
      XKAER(N,NN,NY,NX)=0.0
      XHCER(N,NN,NY,NX)=0.0
      XAL2ER(N,NN,NY,NX)=0.0
      XOH0ER(N,NN,NY,NX)=0.0
      XOH1ER(N,NN,NY,NX)=0.0
      XOH2ER(N,NN,NY,NX)=0.0
      XH1PER(N,NN,NY,NX)=0.0
      XH2PER(N,NN,NY,NX)=0.0
      XOH0EB(N,NN,NY,NX)=0.0
      XOH1EB(N,NN,NY,NX)=0.0
      XOH2EB(N,NN,NY,NX)=0.0
      XH1PEB(N,NN,NY,NX)=0.0
      XH2PEB(N,NN,NY,NX)=0.0
      PALOER(N,NN,NY,NX)=0.0
      PFEOER(N,NN,NY,NX)=0.0
      PCACER(N,NN,NY,NX)=0.0
      PCASER(N,NN,NY,NX)=0.0
      PALPER(N,NN,NY,NX)=0.0
      PFEPER(N,NN,NY,NX)=0.0
      PCPDER(N,NN,NY,NX)=0.0
      PCPHER(N,NN,NY,NX)=0.0
      PCPMER(N,NN,NY,NX)=0.0
      PALPEB(N,NN,NY,NX)=0.0
      PFEPEB(N,NN,NY,NX)=0.0
      PCPDEB(N,NN,NY,NX)=0.0
      PCPHEB(N,NN,NY,NX)=0.0
      PCPMEB(N,NN,NY,NX)=0.0
      DO 9480 K=0,5
      DO 9480 NO=1,7
      OMCER(3,NO,K,N,NN,NY,NX)=0.0
      DO 9480 M=1,2
      OMCER(M,NO,K,N,NN,NY,NX)=0.0
      OMNER(M,NO,K,N,NN,NY,NX)=0.0
      OMPER(M,NO,K,N,NN,NY,NX)=0.0
9480  CONTINUE
      DO 9475 K=0,4
      DO 9470 M=1,2
      ORCER(M,K,N,NN,NY,NX)=0.0
      ORNER(M,K,N,NN,NY,NX)=0.0
      ORPER(M,K,N,NN,NY,NX)=0.0
9470  CONTINUE
      OHCER(K,N,NN,NY,NX)=0.0
      OHNER(K,N,NN,NY,NX)=0.0
      OHPER(K,N,NN,NY,NX)=0.0
      DO 9465 M=1,4
      OSCER(M,K,N,NN,NY,NX)=0.0
      OSAER(M,K,N,NN,NY,NX)=0.0
      OSNER(M,K,N,NN,NY,NX)=0.0
      OSPER(M,K,N,NN,NY,NX)=0.0
9465  CONTINUE
9475  CONTINUE
9855  CONTINUE
      ENDIF
      ENDIF
C
C     GAS AND SOLUTE FLUXES
C
      XCOFLS(N,L,NY,NX)=0.0
      XCHFLS(N,L,NY,NX)=0.0
      XOXFLS(N,L,NY,NX)=0.0
      XNGFLS(N,L,NY,NX)=0.0
      XN2FLS(N,L,NY,NX)=0.0
      XHGFLS(N,L,NY,NX)=0.0
      XN4FLW(N,L,NY,NX)=0.0
      XN3FLW(N,L,NY,NX)=0.0
      XNOFLW(N,L,NY,NX)=0.0
      XNXFLS(N,L,NY,NX)=0.0
      XH1PFS(N,L,NY,NX)=0.0
      XH2PFS(N,L,NY,NX)=0.0
      DO 9860 K=0,4
      XOCFLS(K,N,L,NY,NX)=0.0
      XONFLS(K,N,L,NY,NX)=0.0
      XOPFLS(K,N,L,NY,NX)=0.0
      XOAFLS(K,N,L,NY,NX)=0.0
9860  CONTINUE
9880  CONTINUE
C
C     BAND AND MACROPORE FLUXES
C
      IF(L.NE.0)THEN
      DO 9840 N=1,3
      FLW(N,L,NY,NX)=0.0
      FLWX(N,L,NY,NX)=0.0
      FLWH(N,L,NY,NX)=0.0
      HFLW(N,L,NY,NX)=0.0
      XN4FLB(N,L,NY,NX)=0.0
      XN3FLB(N,L,NY,NX)=0.0
      XNOFLB(N,L,NY,NX)=0.0
      XNXFLB(N,L,NY,NX)=0.0
      XH1BFB(N,L,NY,NX)=0.0
      XH2BFB(N,L,NY,NX)=0.0
      XCOFHS(N,L,NY,NX)=0.0
      XCHFHS(N,L,NY,NX)=0.0
      XOXFHS(N,L,NY,NX)=0.0
      XNGFHS(N,L,NY,NX)=0.0
      XN2FHS(N,L,NY,NX)=0.0
      XHGFHS(N,L,NY,NX)=0.0
      XN4FHW(N,L,NY,NX)=0.0
      XN3FHW(N,L,NY,NX)=0.0
      XNOFHW(N,L,NY,NX)=0.0
      XNXFHS(N,L,NY,NX)=0.0
      XH1PHS(N,L,NY,NX)=0.0
      XH2PHS(N,L,NY,NX)=0.0
      XN4FHB(N,L,NY,NX)=0.0
      XN3FHB(N,L,NY,NX)=0.0
      XNOFHB(N,L,NY,NX)=0.0
      XNXFHB(N,L,NY,NX)=0.0
      XH1BHB(N,L,NY,NX)=0.0
      XH2BHB(N,L,NY,NX)=0.0
      XCOFLG(N,L,NY,NX)=0.0
      XCHFLG(N,L,NY,NX)=0.0
      XOXFLG(N,L,NY,NX)=0.0
      XNGFLG(N,L,NY,NX)=0.0
      XN2FLG(N,L,NY,NX)=0.0
      XN3FLG(N,L,NY,NX)=0.0
      XHGFLG(N,L,NY,NX)=0.0
      DO 9820 K=0,4
      XOCFHS(K,N,L,NY,NX)=0.0
      XONFHS(K,N,L,NY,NX)=0.0
      XOPFHS(K,N,L,NY,NX)=0.0
      XOAFHS(K,N,L,NY,NX)=0.0
9820  CONTINUE
9840  CONTINUE
      ENDIF
9885  CONTINUE
9890  CONTINUE
9895  CONTINUE
C
C     IF SALT FLAG SET
C
      DO 8895 NX=NHW,NHE+1
      DO 8890 NY=NVN,NVS+1
      IF(ISALTG.NE.0)THEN
      DO 8885 L=1,NL(NY,NX)+1
      DO 8880 N=1,3
      IF(L.EQ.1.AND.N.NE.3)THEN
      DO 9836 NN=1,2
      XQRAL(N,NN,NY,NX)=0.0
      XQRFE(N,NN,NY,NX)=0.0
      XQRHY(N,NN,NY,NX)=0.0
      XQRCA(N,NN,NY,NX)=0.0
      XQRMG(N,NN,NY,NX)=0.0
      XQRNA(N,NN,NY,NX)=0.0
      XQRKA(N,NN,NY,NX)=0.0
      XQROH(N,NN,NY,NX)=0.0
      XQRSO(N,NN,NY,NX)=0.0
      XQRCL(N,NN,NY,NX)=0.0
      XQRC3(N,NN,NY,NX)=0.0
      XQRHC(N,NN,NY,NX)=0.0
      XQRAL1(N,NN,NY,NX)=0.0
      XQRAL2(N,NN,NY,NX)=0.0
      XQRAL3(N,NN,NY,NX)=0.0
      XQRAL4(N,NN,NY,NX)=0.0
      XQRALS(N,NN,NY,NX)=0.0
      XQRFE1(N,NN,NY,NX)=0.0
      XQRFE2(N,NN,NY,NX)=0.0
      XQRFE3(N,NN,NY,NX)=0.0
      XQRFE4(N,NN,NY,NX)=0.0
      XQRFES(N,NN,NY,NX)=0.0
      XQRCAO(N,NN,NY,NX)=0.0
      XQRCAC(N,NN,NY,NX)=0.0
      XQRCAH(N,NN,NY,NX)=0.0
      XQRCAS(N,NN,NY,NX)=0.0
      XQRMGO(N,NN,NY,NX)=0.0
      XQRMGC(N,NN,NY,NX)=0.0
      XQRMGH(N,NN,NY,NX)=0.0
      XQRMGS(N,NN,NY,NX)=0.0
      XQRNAC(N,NN,NY,NX)=0.0
      XQRNAS(N,NN,NY,NX)=0.0
      XQRKAS(N,NN,NY,NX)=0.0
      XQRH0P(N,NN,NY,NX)=0.0
      XQRH3P(N,NN,NY,NX)=0.0
      XQRF1P(N,NN,NY,NX)=0.0
      XQRF2P(N,NN,NY,NX)=0.0
      XQRC0P(N,NN,NY,NX)=0.0
      XQRC1P(N,NN,NY,NX)=0.0
      XQRC2P(N,NN,NY,NX)=0.0
      XQRM1P(N,NN,NY,NX)=0.0
9836  CONTINUE
      XQSAL(N,NY,NX)=0.0
      XQSFE(N,NY,NX)=0.0
      XQSHY(N,NY,NX)=0.0
      XQSCA(N,NY,NX)=0.0
      XQSMG(N,NY,NX)=0.0
      XQSNA(N,NY,NX)=0.0
      XQSKA(N,NY,NX)=0.0
      XQSOH(N,NY,NX)=0.0
      XQSSO(N,NY,NX)=0.0
      XQSCL(N,NY,NX)=0.0
      XQSC3(N,NY,NX)=0.0
      XQSHC(N,NY,NX)=0.0
      XQSAL1(N,NY,NX)=0.0
      XQSAL2(N,NY,NX)=0.0
      XQSAL3(N,NY,NX)=0.0
      XQSAL4(N,NY,NX)=0.0
      XQSALS(N,NY,NX)=0.0
      XQSFE1(N,NY,NX)=0.0
      XQSFE2(N,NY,NX)=0.0
      XQSFE3(N,NY,NX)=0.0
      XQSFE4(N,NY,NX)=0.0
      XQSFES(N,NY,NX)=0.0
      XQSCAO(N,NY,NX)=0.0
      XQSCAC(N,NY,NX)=0.0
      XQSCAH(N,NY,NX)=0.0
      XQSCAS(N,NY,NX)=0.0
      XQSMGO(N,NY,NX)=0.0
      XQSMGC(N,NY,NX)=0.0
      XQSMGH(N,NY,NX)=0.0
      XQSMGS(N,NY,NX)=0.0
      XQSNAC(N,NY,NX)=0.0
      XQSNAS(N,NY,NX)=0.0
      XQSKAS(N,NY,NX)=0.0
      XQSH0P(N,NY,NX)=0.0
      XQSH1P(N,NY,NX)=0.0
      XQSH3P(N,NY,NX)=0.0
      XQSF1P(N,NY,NX)=0.0
      XQSF2P(N,NY,NX)=0.0
      XQSC0P(N,NY,NX)=0.0
      XQSC1P(N,NY,NX)=0.0
      XQSC2P(N,NY,NX)=0.0
      XQSM1P(N,NY,NX)=0.0
      ENDIF
      XALFLS(N,L,NY,NX)=0.0
      XFEFLS(N,L,NY,NX)=0.0
      XHYFLS(N,L,NY,NX)=0.0
      XCAFLS(N,L,NY,NX)=0.0
      XMGFLS(N,L,NY,NX)=0.0
      XNAFLS(N,L,NY,NX)=0.0
      XKAFLS(N,L,NY,NX)=0.0
      XOHFLS(N,L,NY,NX)=0.0
      XSOFLS(N,L,NY,NX)=0.0
      XCLFLS(N,L,NY,NX)=0.0
      XC3FLS(N,L,NY,NX)=0.0
      XHCFLS(N,L,NY,NX)=0.0
      XAL1FS(N,L,NY,NX)=0.0
      XAL2FS(N,L,NY,NX)=0.0
      XAL3FS(N,L,NY,NX)=0.0
      XAL4FS(N,L,NY,NX)=0.0
      XALSFS(N,L,NY,NX)=0.0
      XFE1FS(N,L,NY,NX)=0.0
      XFE2FS(N,L,NY,NX)=0.0
      XFE3FS(N,L,NY,NX)=0.0
      XFE4FS(N,L,NY,NX)=0.0
      XFESFS(N,L,NY,NX)=0.0
      XCAOFS(N,L,NY,NX)=0.0
      XCACFS(N,L,NY,NX)=0.0
      XCAHFS(N,L,NY,NX)=0.0
      XCASFS(N,L,NY,NX)=0.0
      XMGOFS(N,L,NY,NX)=0.0
      XMGCFS(N,L,NY,NX)=0.0
      XMGHFS(N,L,NY,NX)=0.0
      XMGSFS(N,L,NY,NX)=0.0
      XNACFS(N,L,NY,NX)=0.0
      XNASFS(N,L,NY,NX)=0.0
      XKASFS(N,L,NY,NX)=0.0
      XH0PFS(N,L,NY,NX)=0.0
      XH3PFS(N,L,NY,NX)=0.0
      XF1PFS(N,L,NY,NX)=0.0
      XF2PFS(N,L,NY,NX)=0.0
      XC0PFS(N,L,NY,NX)=0.0
      XC1PFS(N,L,NY,NX)=0.0
      XC2PFS(N,L,NY,NX)=0.0
      XM1PFS(N,L,NY,NX)=0.0
      XH0BFB(N,L,NY,NX)=0.0
      XH3BFB(N,L,NY,NX)=0.0
      XF1BFB(N,L,NY,NX)=0.0
      XF2BFB(N,L,NY,NX)=0.0
      XC0BFB(N,L,NY,NX)=0.0
      XC1BFB(N,L,NY,NX)=0.0
      XC2BFB(N,L,NY,NX)=0.0
      XM1BFB(N,L,NY,NX)=0.0
      XALFHS(N,L,NY,NX)=0.0
      XFEFHS(N,L,NY,NX)=0.0
      XHYFHS(N,L,NY,NX)=0.0
      XCAFHS(N,L,NY,NX)=0.0
      XMGFHS(N,L,NY,NX)=0.0
      XNAFHS(N,L,NY,NX)=0.0
      XKAFHS(N,L,NY,NX)=0.0
      XOHFHS(N,L,NY,NX)=0.0
      XSOFHS(N,L,NY,NX)=0.0
      XCLFHS(N,L,NY,NX)=0.0
      XC3FHS(N,L,NY,NX)=0.0
      XHCFHS(N,L,NY,NX)=0.0
      XAL1HS(N,L,NY,NX)=0.0
      XAL2HS(N,L,NY,NX)=0.0
      XAL3HS(N,L,NY,NX)=0.0
      XAL4HS(N,L,NY,NX)=0.0
      XALSHS(N,L,NY,NX)=0.0
      XFE1HS(N,L,NY,NX)=0.0
      XFE2HS(N,L,NY,NX)=0.0
      XFE3HS(N,L,NY,NX)=0.0
      XFE4HS(N,L,NY,NX)=0.0
      XFESHS(N,L,NY,NX)=0.0
      XCAOHS(N,L,NY,NX)=0.0
      XCACHS(N,L,NY,NX)=0.0
      XCAHHS(N,L,NY,NX)=0.0
      XCASHS(N,L,NY,NX)=0.0
      XMGOHS(N,L,NY,NX)=0.0
      XMGCHS(N,L,NY,NX)=0.0
      XMGHHS(N,L,NY,NX)=0.0
      XMGSHS(N,L,NY,NX)=0.0
      XNACHS(N,L,NY,NX)=0.0
      XNASHS(N,L,NY,NX)=0.0
      XKASHS(N,L,NY,NX)=0.0
      XH0PHS(N,L,NY,NX)=0.0
      XH3PHS(N,L,NY,NX)=0.0
      XF1PHS(N,L,NY,NX)=0.0
      XF2PHS(N,L,NY,NX)=0.0
      XC0PHS(N,L,NY,NX)=0.0
      XC1PHS(N,L,NY,NX)=0.0
      XC2PHS(N,L,NY,NX)=0.0
      XM1PHS(N,L,NY,NX)=0.0
      XH0BHB(N,L,NY,NX)=0.0
      XH3BHB(N,L,NY,NX)=0.0
      XF1BHB(N,L,NY,NX)=0.0
      XF2BHB(N,L,NY,NX)=0.0
      XC0BHB(N,L,NY,NX)=0.0
      XC1BHB(N,L,NY,NX)=0.0
      XC2BHB(N,L,NY,NX)=0.0
      XM1BHB(N,L,NY,NX)=0.0
8880  CONTINUE
8885  CONTINUE
      ENDIF
8890  CONTINUE
8895  CONTINUE
C
C     RESET SOIL PROPERTIES AND PEDOTRANSFER FUNCTIONS
C     FOLLOWING ANY SOIL DISTURBANCE
C
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      IF(J.EQ.1)THEN
      IFLGT(NY,NX)=0
      DO 9905 NZ=1,NP(NY,NX)
      PSILZ(NZ,NY,NX)=0.0
9905  CONTINUE
      ENDIF
C
C     HYDROLOGICAL PRPOERTIES OR SURFACE LITTER
C
C     VOLWRX=liter water holding capacity
C     VOLR=dry litter volume
C     POROS0,FC,WP=litter porosity,field capacity,wilting point
C     
      VOLWRX(NY,NX)=AMAX1(0.0,THETRX(0)*RC0(0,NY,NX)
     2+THETRX(1)*RC0(1,NY,NX)+THETRX(2)*RC0(2,NY,NX))
      VOLR(NY,NX)=AMAX1(0.0,RC0(0,NY,NX)*1.0E-06/BKRS(0)
     2+RC0(1,NY,NX)*1.0E-06/BKRS(1)+RC0(2,NY,NX)*1.0E-06/BKRS(2))
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX))THEN
      FVOLR=VOLWRX(NY,NX)/VOLR(NY,NX)
      ELSE
      FVOLR=THETRX(1)/BKRS(1)
      ENDIF
      POROS0(NY,NX)=FVOLR
      FC(0,NY,NX)=0.500*FVOLR
      WP(0,NY,NX)=0.125*FVOLR
      PSL(0,NY,NX)=LOG(POROS0(NY,NX))
      FCL(0,NY,NX)=LOG(FC(0,NY,NX))
      WPL(0,NY,NX)=LOG(WP(0,NY,NX))
      PSD(0,NY,NX)=PSL(0,NY,NX)-FCL(0,NY,NX)
      FCD(0,NY,NX)=FCL(0,NY,NX)-WPL(0,NY,NX)
      SRP(0,NY,NX)=1.00
C
C     RESET SURFACE LITTER PHYSICAL PROPERTIES (DENSITY, TEXTURE)
C     AFTER DISTURBANCES (E.G. TILLAGE, EROSION)
C
C     IFLGS=disturbance flag
C     BKDS,BKDSI=current,initial bulk density
C
C     WRITE(*,1116)'IFLGS',IYRC,I,J,IFLGS(NY,NX)
1116  FORMAT(A8,4I6) 
      IF(IFLGS(NY,NX).NE.0)THEN
      IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      BKDS(0,NY,NX)=BKVL(0,NY,NX)/VOLT(0,NY,NX)
      ELSE
      BKDS(0,NY,NX)=BKRS(1)
      ENDIF      
      THETY(0,NY,NX)=EXP((PSIMX(NY,NX)-LOG(-PSIHY))
     2*FCD(0,NY,NX)/PSIMD(NY,NX)+FCL(0,NY,NX))
      SUM2=0.0
      DO 1220 K=1,100
      XK=K-1
      THETK(K)=POROS0(NY,NX)-(XK/100.0*POROS0(NY,NX))
      IF(THETK(K).LT.FC(0,NY,NX))THEN
      PSISK(K)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCL(0,NY,NX)-LOG(THETK(K)))
     3/FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETK(K).LT.POROS0(NY,NX))THEN 
      PSISK(K)=-EXP(PSIMS(NY,NX)
     2+(((PSL(0,NY,NX)-LOG(THETK(K)))
     3/PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
      ELSE
      PSISK(K)=PSISE(0,NY,NX)
      ENDIF
      SUM2=SUM2+(2*K-1)/(PSISK(K)**2)
1220  CONTINUE
      DO 1235 K=1,100
      SUM1=0.0
      XK=K-1
      YK=((100.0-XK)/100.0)**1.33
      DO 1230 M=K,100
      SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2)
1230  CONTINUE
      HCND(3,K,0,NY,NX)=SCNV(0,NY,NX)*YK*SUM1/SUM2
      HCND(1,K,0,NY,NX)=0.0
      HCND(2,K,0,NY,NX)=0.0
      IF(K.GT.1.AND.PSISK(K).LT.PSISA(0,NY,NX)
     2.AND.PSISK(K-1).GE.PSISA(0,NY,NX))THEN
      THETS(0,NY,NX)=THETK(K)
      ENDIF
C     WRITE(*,3534)'PSI0',K,THETK(K),PSISK(K),HCND(3,K,0,NY,NX) 
C    2,PSL(0,NY,NX),LOG(THETK(K)),POROS0(NY,NX),FC(0,NY,NX),WP(0,NY,NX)
C    3,SRP(0,NY,NX),VOLWRX(NY,NX),VOLR(NY,NX),THETS(0,NY,NX)
C    4,PSISA(0,NY,NX)
3534  FORMAT(A8,1I4,20E12.4)
1235  CONTINUE
C
C     RESET SOIL PHYSICAL PROPERTIES (DENSITY, TEXTURE)
C     AFTER DISTURBANCES (E.G. TILLAGE, EROSION)
C
      DO 9975 L=NUI(NY,NX),NLI(NY,NX)
C
C     AREA,DLYR=lateral(1,2), vertical(3) area,thickness of soil layer 
C     VOLT,VOLX,VOLY=layer volume including,excluding rock,macropores
C
      IF(BKDS(L,NY,NX).LE.ZERO.AND.DLYR(3,L,NY,NX).LE.ZERO2)THEN
      VOLW(L,NY,NX)=0.0
      VOLI(L,NY,NX)=0.0
      ENDIF
      AREA(1,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(2,L,NY,NX)
      AREA(2,L,NY,NX)=DLYR(3,L,NY,NX)*DLYR(1,L,NY,NX)
      VOLT(L,NY,NX)=AREA(3,L,NY,NX)*DLYR(3,L,NY,NX)
      VOLX(L,NY,NX)=VOLT(L,NY,NX)*FMPR(L,NY,NX)
      IF(BKDS(L,NY,NX).LE.ZERO)THEN
      VOLY(L,NY,NX)=VOLX(L,NY,NX)
      ENDIF
C     IF(NX.EQ.1)THEN
C     WRITE(*,443)'VOLT',I,J,NX,NY,L,VOLT(L,NY,NX)
C    2,VOLX(L,NY,NX),DLYR(3,L,NY,NX),AREA(3,L,NY,NX)
443   FORMAT(A8,5I4,20E14.6)
C     ENDIF
C
C     BKVL=soil mass
C     C*=concentration,ORGC=SOC,SAND=sand,SILT=silt,CLAY=clay
C     PTDS=particle density
C     PTDSNU=particle density of surface layer for use in erosion.f 
C     POROS=porosity used in diffusivity
C     VOLA,VOLW,VOLI,VOLP=total,water-,ice-,air-filled micropore volume
C     VOLAH,VOLWH,VOLIH,VOLPH=total,water-,ice-,air-filled macropore volume
C     EHUM=fraction of microbial decomposition product allocated to humus 
C     EPOC=fraction of SOC decomposition product allocated to POC
C     SRP=parameter for deviation from linear log-log water retention 
C     PSIMX,PSIMN,PSIMS=log water potential at FC,WP,POROS
C     PSISD,PSIMD=PSIMX-PSIMS,PSIMN-PSIMX
C     FC,WP=water contents at field capacity,wilting point 
C     FCL,WPL=log FC,WP
C     FCD,PSD=FCL-WPL,log(POROS)-FCL
C
      BKVL(L,NY,NX)=BKDS(L,NY,NX)*VOLX(L,NY,NX)
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CORGC(L,NY,NX)=AMIN1(0.55E+06,ORGC(L,NY,NX)/BKVL(L,NY,NX))
      CSAND(L,NY,NX)=SAND(L,NY,NX)/BKVL(L,NY,NX)
      CSILT(L,NY,NX)=SILT(L,NY,NX)/BKVL(L,NY,NX)
      CCLAY(L,NY,NX)=CLAY(L,NY,NX)/BKVL(L,NY,NX)
      ELSE
      CORGC(L,NY,NX)=0.0
      CSAND(L,NY,NX)=0.0
      CSILT(L,NY,NX)=0.0
      CCLAY(L,NY,NX)=0.0
      ENDIF
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      CORGCM=AMAX1(0.0,AMIN1(1.0,1.82E-06*CORGC(L,NY,NX)))
      PTDS=1.30*CORGCM+2.66*(1.0-CORGCM)
      IF(L.EQ.NU(NY,NX))THEN
      POROS(L,NY,NX)=AMAX1(POROS(L,NY,NX),1.0-(BKDS(L,NY,NX)/PTDS))
      ELSE
      POROS(L,NY,NX)=1.0-(BKDS(L,NY,NX)/PTDS)
      ENDIF
      ELSE
      PTDS=0.0
      POROS(L,NY,NX)=1.0
      ENDIF
C     VOLA(L,NY,NX)=AMAX1(POROS(L,NY,NX)*VOLY(L,NY,NX)
C    2,VOLW(L,NY,NX)+VOLI(L,NY,NX))
C     VOLAH(L,NY,NX)=AMAX1(FHOL(L,NY,NX)*VOLT(L,NY,NX)
C    2,VOLWH(L,NY,NX)+VOLIH(L,NY,NX))
      VOLA(L,NY,NX)=POROS(L,NY,NX)*VOLY(L,NY,NX)
      VOLAH(L,NY,NX)=FHOL(L,NY,NX)*VOLT(L,NY,NX)
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLP(L,NY,NX)=AMAX1(0.0,VOLA(L,NY,NX)-VOLW(L,NY,NX)
     2-VOLI(L,NY,NX))+AMAX1(0.0,VOLAH(L,NY,NX)-VOLWH(L,NY,NX)
     3-VOLIH(L,NY,NX))
      ELSE
      VOLP(L,NY,NX)=0.0
      ENDIF
      EHUM(L,NY,NX)=0.200+0.333*AMIN1(0.5,CCLAY(L,NY,NX))
C    2+0.167E-06*CORGC(L,NY,NX)
C     WRITE(*,3331)'EHUM',I,J,L,EHUM(L,NY,NX),CCLAY(L,NY,NX)
C    3,CSILT(L,NY,NX),CSAND(L,NY,NX),CORGC(L,NY,NX),ORGC(L,NY,NX)
C    4,VOLA(L,NY,NX),POROS(L,NY,NX),VOLX(L,NY,NX),BKVL(L,NY,NX)
3331  FORMAT(A8,3I4,12E12.4)
      EPOC(L,NY,NX)=1.0
      IF(CORGC(L,NY,NX).GT.FORGC)THEN
      SRP(L,NY,NX)=0.25
      ELSEIF(CORGC(L,NY,NX).GT.0.5*FORGC)THEN
      SRP(L,NY,NX)=0.33
      ELSE
      SRP(L,NY,NX)=1.00
      ENDIF
      PSL(L,NY,NX)=LOG(POROS(L,NY,NX))
      IF((ISOIL(1,L,NY,NX).EQ.0.AND.ISOIL(2,L,NY,NX).EQ.0)
     2.OR.DATA(20).EQ.'YES')THEN
      FCL(L,NY,NX)=LOG(FC(L,NY,NX))
      WPL(L,NY,NX)=LOG(WP(L,NY,NX))
      PSD(L,NY,NX)=PSL(L,NY,NX)-FCL(L,NY,NX)
      FCD(L,NY,NX)=FCL(L,NY,NX)-WPL(L,NY,NX)
      ELSE 
C
C     DEFAULT SOIL HYDROLOGIC PPTYS (FIELD CAPACITY, WILTING POINT)
C     IF ACTUAL VALUES WERE NOT INPUT TO THE SOIL FILE
C
C     THW,THI=initial soil water,ice content from soil file
C
      IF(DATA(20).EQ.'NO')THEN
      IF(ISOIL(1,L,NY,NX).EQ.1.OR.ISOIL(2,L,NY,NX).EQ.1)THEN
      IF(CORGC(L,NY,NX).LT.FORGW)THEN
      FC(L,NY,NX)=0.2576-0.20*CSAND(L,NY,NX)
     2+0.36*CCLAY(L,NY,NX)+0.60E-06*CORGC(L,NY,NX)
      ELSE
      IF(BKDS(L,NY,NX).LT.0.075)THEN
      FC(L,NY,NX)=0.27
      ELSEIF(BKDS(L,NY,NX).LT.0.195)THEN
      FC(L,NY,NX)=0.62
      ELSE
      FC(L,NY,NX)=0.71
      ENDIF
      ENDIF
      FC(L,NY,NX)=FC(L,NY,NX)/(1.0-FHOL(L,NY,NX))
      FC(L,NY,NX)=AMIN1(0.75*POROS(L,NY,NX),FC(L,NY,NX))
C     WRITE(*,3332)'FC',IYRC,I,J,L,FC(L,NY,NX),CCLAY(L,NY,NX)
C    2,CORGC(L,NY,NX),CSAND(L,NY,NX)
3332  FORMAT(A8,4I6,20E12.4)
      IF(CORGC(L,NY,NX).LT.FORGW)THEN
      WP(L,NY,NX)=0.0260+0.50*CCLAY(L,NY,NX)
     2+0.32E-06*CORGC(L,NY,NX)
      ELSE
      IF(BKDS(L,NY,NX).LT.0.075)THEN
      WP(L,NY,NX)=0.04
      ELSEIF(BKDS(L,NY,NX).LT.0.195)THEN
      WP(L,NY,NX)=0.15
      ELSE
      WP(L,NY,NX)=0.22
      ENDIF
      ENDIF
      WP(L,NY,NX)=WP(L,NY,NX)/(1.0-FHOL(L,NY,NX))
      WP(L,NY,NX)=AMIN1(0.75*FC(L,NY,NX),WP(L,NY,NX))
C     WRITE(*,3332)'WP',IYRC,I,J,L,WP(L,NY,NX),CCLAY(L,NY,NX)
C    2,CORGC(L,NY,NX),FC(L,NY,NX)
      ENDIF
      FCL(L,NY,NX)=LOG(FC(L,NY,NX))
      WPL(L,NY,NX)=LOG(WP(L,NY,NX))
      PSD(L,NY,NX)=PSL(L,NY,NX)-FCL(L,NY,NX)
      FCD(L,NY,NX)=FCL(L,NY,NX)-WPL(L,NY,NX)
      ENDIF
      IF(I.EQ.IBEGIN.AND.J.EQ.1.AND.IYRC.EQ.IDATA(9))THEN
      IF(THW(L,NY,NX).GT.1.0.OR.DPTH(L,NY,NX).GE.DTBLZ(NY,NX))THEN
      THETW(L,NY,NX)=POROS(L,NY,NX)
      ELSEIF(THW(L,NY,NX).EQ.1.0)THEN 
      THETW(L,NY,NX)=FC(L,NY,NX)
      ELSEIF(THW(L,NY,NX).EQ.0.0)THEN 
      THETW(L,NY,NX)=WP(L,NY,NX)
      ELSEIF(THW(L,NY,NX).LT.0.0)THEN 
      THETW(L,NY,NX)=0.0
      ENDIF
      IF(THI(L,NY,NX).GT.1.0.OR.DPTH(L,NY,NX).GE.DTBLZ(NY,NX))THEN
      THETI(L,NY,NX)=AMAX1(0.0,AMIN1(POROS(L,NY,NX)
     2,POROS(L,NY,NX)-THW(L,NY,NX)))
      ELSEIF(THI(L,NY,NX).EQ.1.0)THEN 
      THETI(L,NY,NX)=AMAX1(0.0,AMIN1(FC(L,NY,NX)
     2,POROS(L,NY,NX)-THW(L,NY,NX)))
      ELSEIF(THI(L,NY,NX).EQ.0.0)THEN 
      THETI(L,NY,NX)=AMAX1(0.0,AMIN1(WP(L,NY,NX)
     2,POROS(L,NY,NX)-THW(L,NY,NX)))
      ELSEIF(THI(L,NY,NX).LT.0.0)THEN 
      THETI(L,NY,NX)=0.0
      ENDIF
      IF(DATA(20).EQ.'NO')THEN
      VOLW(L,NY,NX)=THETW(L,NY,NX)*VOLX(L,NY,NX)
      VOLWX(L,NY,NX)=VOLW(L,NY,NX)
      VOLWH(L,NY,NX)=THETW(L,NY,NX)*VOLAH(L,NY,NX)
      VOLI(L,NY,NX)=THETI(L,NY,NX)*VOLX(L,NY,NX)
      VOLIH(L,NY,NX)=THETI(L,NY,NX)*VOLAH(L,NY,NX)
      VHCP(L,NY,NX)=VHCM(L,NY,NX)+4.19*(VOLW(L,NY,NX)
     2+VOLWH(L,NY,NX))+1.9274*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
      THETWZ(L,NY,NX)=THETW(L,NY,NX) 
      THETIZ(L,NY,NX)=THETI(L,NY,NX)
      ENDIF
      ENDIF
      ENDIF
      VOLP(L,NY,NX)=AMAX1(0.0,VOLA(L,NY,NX)-VOLW(L,NY,NX)
     2-VOLI(L,NY,NX))+AMAX1(0.0,VOLAH(L,NY,NX)-VOLWH(L,NY,NX)
     3-VOLIH(L,NY,NX))
      IF(VOLT(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETP(L,NY,NX)=VOLP(L,NY,NX)/VOLY(L,NY,NX)
      ELSE
      THETP(L,NY,NX)=0.0
      ENDIF
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      THETY(L,NY,NX)=EXP((PSIMX(NY,NX)-LOG(-PSIHY))
     2*FCD(L,NY,NX)/PSIMD(NY,NX)+FCL(L,NY,NX))
      ELSE
      THETY(L,NY,NX)=ZERO2
      ENDIF
C
C     SATURATED HYDRAULIC CONDUCTIVITY FROM SWC AT SATURATION VS.
C     -0.033 MPA (MINERAL SOILS) IF NOT ENTERED IN SOIL FILE IN 'READS'
C
C     SCNV,SCNH=vertical,lateral saturated hydraulic conductivity
C
      IF(ISOIL(3,L,NY,NX).EQ.1)THEN
      IF(CORGC(L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS(L,NY,NX),EXP((PSIMS(NY,NX)-LOG(0.033))
     2*(PSL(L,NY,NX)-FCL(L,NY,NX))/PSISD(NY,NX)+PSL(L,NY,NX)))
      SCNV(L,NY,NX)=1.54*((POROS(L,NY,NX)-THETF)/THETF)**2
      ELSE
      SCNV(L,NY,NX)=0.10+75.0*1.0E-15**BKDS(L,NY,NX)
      SCNV(L,NY,NX)=SCNV(L,NY,NX)*FMPR(L,NY,NX)
      ENDIF
C     WRITE(*,3332)'SCNV',IYRC,I,J,L,SCNV(L,NY,NX),POROS(L,NY,NX) 
C    2,THETF,FMPR(L,NY,NX),PSIMS(NY,NX),LOG(0.033)
C    3,PSL(L,NY,NX),FCL(L,NY,NX),PSISD(NY,NX)
      ENDIF
      IF(ISOIL(4,L,NY,NX).EQ.1)THEN
      IF(CORGC(L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS(L,NY,NX),EXP((PSIMS(NY,NX)-LOG(0.033))
     2*(PSL(L,NY,NX)-FCL(L,NY,NX))/PSISD(NY,NX)+PSL(L,NY,NX)))
      SCNH(L,NY,NX)=1.54*((POROS(L,NY,NX)-THETF)/THETF)**2
      ELSE
      SCNH(L,NY,NX)=0.10+75.0*1.0E-15**BKDS(L,NY,NX)
      SCNH(L,NY,NX)=SCNH(L,NY,NX)*FMPR(L,NY,NX)
      ENDIF
C     WRITE(*,3332)'SCNH',IYRC,I,J,L,SCNH(L,NY,NX),POROS(L,NY,NX) 
C    2,THETF,FMPR(L,NY,NX)
      ENDIF
C     WRITE(*,3333)'PPTYS',I,J,NX,NY,L,IBEGIN,IYRC,IDATA(9)
C    2,ISOIL(1,L,NY,NX),ISOIL(2,L,NY,NX)
C    3,ISOIL(3,L,NY,NX),ISOIL(4,L,NY,NX)
C    3,SCNV(L,NY,NX),SCNH(L,NY,NX),POROS(L,NY,NX) 
C    2,FC(L,NY,NX),WP(L,NY,NX),BKDS(L,NY,NX),THW(L,NY,NX)
C    3,THETW(L,NY,NX),THI(L,NY,NX),THETI(L,NY,NX),VOLI(L,NY,NX)
3333  FORMAT(A8,12I6,20E12.4)
C
C     HYDRAULIC CONDUCTIVITY FUNCTION FROM KSAT AND SOIL WATER RELEASE CURVE
C
C     THETK,PSISK=micropore class water content,potential
C     HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity 
C     
C     IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      SUM2=0.0
      DO 1320 K=1,100
      XK=K-1
      THETK(K)=POROS(L,NY,NX)-(XK/100.0*POROS(L,NY,NX))
      IF(THETK(K).LT.FC(L,NY,NX))THEN
      PSISK(K)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCL(L,NY,NX)-LOG(THETK(K)))
     3/FCD(L,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETK(K).LT.POROS(L,NY,NX)-DTHETW)THEN 
      PSISK(K)=-EXP(PSIMS(NY,NX)
     2+(((PSL(L,NY,NX)-LOG(THETK(K)))
     3/PSD(L,NY,NX))**SRP(L,NY,NX)*PSISD(NY,NX)))
      ELSE
      PSISK(K)=PSISE(L,NY,NX)
      ENDIF
      SUM2=SUM2+(2*K-1)/(PSISK(K)**2)
1320  CONTINUE
      DO 1335 K=1,100
      SUM1=0.0
      XK=K-1
      YK=((100.0-XK)/100.0)**1.33
      DO 1330 M=K,100
      SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2)
1330  CONTINUE
      DO 1340 N=1,3
      IF(N.EQ.3)THEN
      HCND(N,K,L,NY,NX)=SCNV(L,NY,NX)*YK*SUM1/SUM2
      IF(K.GT.1.AND.PSISK(K).LT.PSISA(L,NY,NX)
     2.AND.PSISK(K-1).GE.PSISA(L,NY,NX))THEN
      THETS(L,NY,NX)=THETK(K)
      ENDIF
C     WRITE(*,3536)'PSI',L,K,THETK(K),PSISK(K),HCND(N,K,L,NY,NX) 
C    2,PSL(L,NY,NX),LOG(THETK(K)),POROS(L,NY,NX),FC(L,NY,NX),WP(L,NY,NX)
C    3,SRP(L,NY,NX),THETS(L,NY,NX),PSISA(L,NY,NX)
3536  FORMAT(A8,2I4,12E12.4)
      ELSE
      HCND(N,K,L,NY,NX)=SCNH(L,NY,NX)*YK*SUM1/SUM2
      ENDIF
1340  CONTINUE
1335  CONTINUE
2340  CONTINUE
2335  CONTINUE
C     ENDIF
C
C     SOIL MACROPORE DIMENSIONS AND CONDUCTIVITY FROM MACROPORE FRACTION
C     ENTERED IN 'READS'
C
C     PHOL,NHOL,HRAD=path length between, number,radius of macropores
C     CNDH=macropore hydraulic conductivity
C
      HRAD(L,NY,NX)=0.5E-03
      NHOL(L,NY,NX)=INT(VOLAH(L,NY,NX)/(3.1416*HRAD(L,NY,NX)**2
     2*VOLTI(L,NY,NX)))
      IF(NHOL(L,NY,NX).GT.0.0)THEN
      PHOL(L,NY,NX)=1.0/(SQRT(3.1416*NHOL(L,NY,NX)))
      ELSE
      PHOL(L,NY,NX)=1.0
      ENDIF
      VISCWL=VISCW*EXP(0.533-0.0267*TCS(L,NY,NX))
      CNDH(L,NY,NX)=3.6E+03*3.1416*NHOL(L,NY,NX)*HRAD(L,NY,NX)**4
     2/(8.0*VISCWL) 
C
C     SOIL HEAT CAPACITY AND THERMAL CONDUCTIVITY OF SOLID PHASE
C     FROM SOC AND TEXTURE
C
C     VORGC,VMINL,VSAND=volume fractions of SOC,mineral,sand
C     STC,DTC=weighted thermal conductivity of soil solid component
C
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VORGC=CORGCM*BKDS(L,NY,NX)/PTDS
      VMINL=(CSILT(L,NY,NX)+CCLAY(L,NY,NX))*BKDS(L,NY,NX)/PTDS
      VSAND=CSAND(L,NY,NX)*BKDS(L,NY,NX)/PTDS
      STC(L,NY,NX)=(1.253*VORGC*9.050E-04+0.514*VMINL*1.056E-02
     2+0.386*VSAND*2.112E-02)*FMPR(L,NY,NX)
     3+0.514*ROCK(L,NY,NX)*1.056E-02
      DTC(L,NY,NX)=(1.253*VORGC+0.514*VMINL+0.386*VSAND)
     2*FMPR(L,NY,NX)+0.514*ROCK(L,NY,NX)
      ELSE
      STC(L,NY,NX)=0.0
      DTC(L,NY,NX)=0.0
      ENDIF
9975  CONTINUE
C
C     SURFACE RESIDUE PROPERTIES
C
C     FCR=litter water content at -0.01 MPa
C     THETY=litter hygroscopic water content
C
      CORGC(0,NY,NX)=0.55E+06
C
C     SOIL SURFACE WATER STORAGE CAPACITY
C
C     IDTBL=water table flag from site file
C     DTBLX,DTBLZ=current,initial natural water table depth
C     DTBLY,DTBLD=current,initial artificial water table depth
C     ZS,ZW=soil,water surface roughness
C     VOLWD=soil surface water retention capacity
C     VOLWG=VOLWD accounting for above-ground water table
C     EHUM=fraction of microbial decompn product allocated to surface humus 
C     EPOC=fraction of SOC decomposition product allocated to surface POC
C     
      IF(IDTBL(NY,NX).LE.1.OR.IDTBL(NY,NX).EQ.3)THEN
      DTBLX(NY,NX)=DTBLZ(NY,NX)
      ELSEIF(IDTBL(NY,NX).EQ.2.OR.IDTBL(NY,NX).EQ.4)THEN
      DTBLX(NY,NX)=DTBLZ(NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      ENDIF
      IF(IDTBL(NY,NX).EQ.3.OR.IDTBL(NY,NX).EQ.4)THEN
      DTBLY(NY,NX)=DTBLD(NY,NX)
      ENDIF
C     IF(J.EQ.24)THEN
C     WRITE(*,1114)'DTBLX',I,J,IYRC,NX,NY,IDTBL(NY,NX),DTBLX(NY,NX)
C    2,DTBLZ(NY,NX),CDPTH(NU(NY,NX)-1,NY,NX)
C    3,(DPTH(L,NY,NX),L=1,NL(NY,NX))
C    3,BKDS(NU(NY,NX),NY,NX),BKDSI(NU(NY,NX),NY,NX)
C    4,DLYR(3,NU(NY,NX),NY,NX),DLYRI(3,NU(NY,NX),NY,NX)   
1114  FORMAT(A8,6I4,30E12.4)
C     ENDIF
      IF(BKDS(NU(NY,NX),NY,NX).GT.ZERO)THEN
      ZS(NY,NX)=0.020
      ELSE
      ZS(NY,NX)=ZW
      ENDIF
      VOLWD(NY,NX)=AMAX1(0.001,0.112*ZS(NY,NX)+3.10*ZS(NY,NX)**2
     2-0.012*ZS(NY,NX)*SLOPE(0,NY,NX))*AREA(3,NU(NY,NX),NY,NX)
      VOLWG(NY,NX)=AMAX1(VOLWD(NY,NX)
     2,-(DTBLX(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX))
     3*AREA(3,NU(NY,NX),NY,NX))
C     WRITE(*,6631)'VOLWG',I,J,NY,NX,NU(NY,NX),VOLWG(NY,NX)
C    2,DTBLX(NY,NX),CDPTH(NU(NY,NX)-1,NY,NX)
C    3,AREA(3,NU(NY,NX),NY,NX),VOLWD(NY,NX),ZS(NY,NX),SLOPE(0,NY,NX)
6631  FORMAT(A8,5I4,12E12.4)
      DPTH(NU(NY,NX),NY,NX)=CDPTH(NU(NY,NX),NY,NX)
     2-0.5*DLYR(3,NU(NY,NX),NY,NX)
      IF(BKVL(NU(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      CCLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)
     2/BKVL(NU(NY,NX),NY,NX)
      CSILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)
     2/BKVL(NU(NY,NX),NY,NX)
      CSAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)
     2/BKVL(NU(NY,NX),NY,NX)
      ELSE
      CCLAY(NU(NY,NX),NY,NX)=0.0
      CSILT(NU(NY,NX),NY,NX)=0.0
      CSAND(NU(NY,NX),NY,NX)=0.0
      ENDIF
      EHUM(0,NY,NX)=0.200+0.333*AMIN1(0.5,CCLAY(NU(NY,NX),NY,NX))
C    2+0.167E-06*CORGC(NU(NY,NX),NY,NX)
      EPOC(0,NY,NX)=0.150
C
C     IFLGS=reset disturbance flag
C
      IFLGS(NY,NX)=0
      ENDIF
C
C     END OF RESET AFTER DISTURBANCE
C
C     PARAMETERS FOR COHESION, EROSIVITY, AND ROUGHNESS OF SURFACE SOIL USED
C     FOR SURFACE WATER AND SEDIMENT TRANSPORT IN 'EROSION'
C
C     DETS=soil detachability from rainfall impact 
C     D50=average particle size  
C     CER,XER=parameters for runoff transport capacity
C     ZD50=particle size effect on surface roughness 
C     VLS=hourly sinking rate 
C     COHS=soil cohesion
C     DETE=soil detachability    
C     ZM=surface roughness used in runoff velocity calculation in watsub.f
C
      BKVLNU(NY,NX)=AMAX1(0.0,BKVLNM(NY,NX)
     3+1.82E-06*ORGC(NU(NY,NX),NY,NX)) 
C     WRITE(*,2423)'BKVLH',I,J,NX,NY,BKVL(NU(N2,N1),N2,N1)
C    2,BKVLNM(N2,N1),BKVLNU(N2,N1),ORGR(NU(N2,N1),N2,N1)
C    3,ORGC(NU(N2,N1),N2,N1),SAND(NU(N2,N1),N2,N1)
C    3,SILT(NU(N2,N1),N2,N1),CLAY(NU(N2,N1),N2,N1)
2423  FORMAT(A8,4I4,20E12.4)
      BKVLNX=SAND(NU(NY,NX),NY,NX)+SILT(NU(NY,NX),NY,NX)
     2+CLAY(NU(NY,NX),NY,NX)+1.82E-06*ORGC(NU(NY,NX),NY,NX) 
      IF(BKVLNX.GT.ZEROS(NY,NX))THEN
      CORGM=1.82E-06*ORGC(NU(NY,NX),NY,NX)/BKVLNX
      CORGC(NU(NY,NX),NY,NX)=0.55E+06*CORGM
      CSAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)/BKVLNX
      CSILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)/BKVLNX
      CCLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)/BKVLNX
      ELSE
      CORGM=0.0
      CORGC(NU(NY,NX),NY,NX)=0.0
      CSAND(NU(NY,NX),NY,NX)=0.0
      CSILT(NU(NY,NX),NY,NX)=1.0
      CCLAY(NU(NY,NX),NY,NX)=0.0
      ENDIF
      IF(IERSNG.EQ.2.OR.IERSNG.EQ.3)THEN
      D50=1.0*CCLAY(NU(NY,NX),NY,NX)+10.0*CSILT(NU(NY,NX),NY,NX)
     2+100.0*CSAND(NU(NY,NX),NY,NX)+100.0*CORGM
      ZD50=0.041*(1.0E-06*D50)**0.167
      ZM(NY,NX)=ZS(NY,NX)+ZD50+1.0*VOLR(NY,NX)/AREA(3,0,NY,NX)
      CER(NY,NX)=((D50+5.0)/0.32)**(-0.6)
      XER(NY,NX)=((D50+5.0)/300.0)**0.25
      DETS(NY,NX)=1.0E-06*(1.0+2.0*(1.0-CSILT(NU(NY,NX),NY,NX)-CORGM))
      COHS=2.0+10.0*(CCLAY(NU(NY,NX),NY,NX)+CORGM)
     2+5.0*(1.0-EXP(-2.0E-06*RTDNT(NU(NY,NX),NY,NX)))
      DETE(NY,NX)=0.79*EXP(-0.85*AMAX1(1.0,COHS))
      PTDSNU(NY,NX)=1.30*CORGM+2.66*(1.0-CORGM)
      VISCWL=VISCW*EXP(0.533-0.0267*TCS(0,NY,NX))
      VLS(NY,NX)=3.6E+03*9.8*(PTDSNU(NY,NX)-1.0)
     2*(1.0E-06*D50)**2/(18.0*VISCWL)
C     WRITE(*,1118)'COHS',I,J,NX,NY,NU(NY,NX),COHS,DETE(NY,NX)
C    2,ZM(NY,NX),VLS(NY,NX),D50,ZD50,PTDSNU(NY,NX)
C    3,RTDNT(NU(NY,NX),NY,NX),VOLR(NY,NX)/AREA(3,0,NY,NX)
C    3,ORGC(0,NY,NX)
C    3,CCLAY(NU(NY,NX),NY,NX),CSILT(NU(NY,NX),NY,NX)
C    4,CSAND(NU(NY,NX),NY,NX),CORGM,CORGC(NU(NY,NX),NY,NX)
C    5,BKVL(NU(NY,NX),NY,NX),BKVLNX
C    3,VISCWL,TCS(0,NY,NX)
1118  FORMAT(A8,5I4,20E12.4)
      ENDIF
C
C     RESET HOURLY ACCUMULATORS
C
      UCO2S(NY,NX)=0.0
      TOMT(NY,NX)=0.0
      TONT(NY,NX)=0.0
      TOPT(NY,NX)=0.0
      UVOLW(NY,NX)=0.0
      URSDC(NY,NX)=0.0
      UORGC(NY,NX)=0.0
      URSDN(NY,NX)=0.0
      UORGN(NY,NX)=0.0
      URSDP(NY,NX)=0.0
      UORGP(NY,NX)=0.0
      UNH4(NY,NX)=0.0
      UNO3(NY,NX)=0.0
      UPO4(NY,NX)=0.0
      UPP4(NY,NX)=0.0
      UION(NY,NX)=0.0
      HVOLO(NY,NX)=0.0
      HCO2G(NY,NX)=0.0
      HCH4G(NY,NX)=0.0
      HOXYG(NY,NX)=0.0
      HN2GG(NY,NX)=0.0
      HN2OG(NY,NX)=0.0
      HNH3G(NY,NX)=0.0
      FLWR(NY,NX)=0.0
      HFLWR(NY,NX)=0.0
      THAWR(NY,NX)=0.0
      HTHAWR(NY,NX)=0.0
      HEATI(NY,NX)=0.0
      HEATS(NY,NX)=0.0
      HEATE(NY,NX)=0.0
      HEATV(NY,NX)=0.0
      HEATH(NY,NX)=0.0
      TEVAPG(NY,NX)=0.0
      XCODFS(NY,NX)=0.0
      XCHDFS(NY,NX)=0.0
      XOXDFS(NY,NX)=0.0
      XNGDFS(NY,NX)=0.0
      XN2DFS(NY,NX)=0.0
      XN3DFS(NY,NX)=0.0
      XNBDFS(NY,NX)=0.0
      XHGDFS(NY,NX)=0.0
      XCODFR(NY,NX)=0.0
      XCHDFR(NY,NX)=0.0
      XOXDFR(NY,NX)=0.0
      XNGDFR(NY,NX)=0.0
      XN2DFR(NY,NX)=0.0
      XN3DFR(NY,NX)=0.0
      XHGDFR(NY,NX)=0.0
      TVOLWP(NY,NX)=0.0
      TVOLWC(NY,NX)=0.0
      TFLWCI(NY,NX)=0.0
      TFLWC(NY,NX)=0.0
      TEVAPP(NY,NX)=0.0
      TEVAPC(NY,NX)=0.0
      THFLXC(NY,NX)=0.0
      TENGYC(NY,NX)=0.0
      TCO2Z(NY,NX)=0.0
      TOXYZ(NY,NX)=0.0
      TCH4Z(NY,NX)=0.0
      TN2OZ(NY,NX)=0.0
      TNH3Z(NY,NX)=0.0
      TH2GZ(NY,NX)=0.0
      ZCSNC(NY,NX)=0.0
      ZZSNC(NY,NX)=0.0
      ZPSNC(NY,NX)=0.0
      WTSTGT(NY,NX)=0.0
      PPT(NY,NX)=0.0
      DO 9865 L=1,JS
      FLSW(L,NY,NX)=0.0
      FLSWH(L,NY,NX)=0.0
      HFLSW(L,NY,NX)=0.0
      FLSWR(L,NY,NX)=0.0 
      HFLSWR(L,NY,NX)=0.0
      XFLWS(L,NY,NX)=0.0
      XFLWW(L,NY,NX)=0.0
      XFLWI(L,NY,NX)=0.0
      XHFLWW(L,NY,NX)=0.0  
      XWFLXS(L,NY,NX)=0.0
      XWFLXI(L,NY,NX)=0.0
      XTHAWW(L,NY,NX)=0.0
      XCOBLS(L,NY,NX)=0.0
      XCHBLS(L,NY,NX)=0.0
      XOXBLS(L,NY,NX)=0.0
      XNGBLS(L,NY,NX)=0.0
      XN2BLS(L,NY,NX)=0.0
      XN4BLW(L,NY,NX)=0.0
      XN3BLW(L,NY,NX)=0.0
      XNOBLW(L,NY,NX)=0.0
      XH1PBS(L,NY,NX)=0.0
      XH2PBS(L,NY,NX)=0.0
      IF(ISALTG.NE.0)THEN
      XALBLS(L,NY,NX)=0.0
      XFEBLS(L,NY,NX)=0.0
      XHYBLS(L,NY,NX)=0.0
      XCABLS(L,NY,NX)=0.0
      XMGBLS(L,NY,NX)=0.0
      XNABLS(L,NY,NX)=0.0
      XKABLS(L,NY,NX)=0.0
      XOHBLS(L,NY,NX)=0.0
      XSOBLS(L,NY,NX)=0.0
      XCLBLS(L,NY,NX)=0.0
      XC3BLS(L,NY,NX)=0.0
      XHCBLS(L,NY,NX)=0.0
      XAL1BS(L,NY,NX)=0.0
      XAL2BS(L,NY,NX)=0.0
      XAL3BS(L,NY,NX)=0.0
      XAL4BS(L,NY,NX)=0.0
      XALSBS(L,NY,NX)=0.0
      XFE1BS(L,NY,NX)=0.0
      XFE2BS(L,NY,NX)=0.0
      XFE3BS(L,NY,NX)=0.0
      XFE4BS(L,NY,NX)=0.0
      XFESBS(L,NY,NX)=0.0
      XCAOBS(L,NY,NX)=0.0
      XCACBS(L,NY,NX)=0.0
      XCAHBS(L,NY,NX)=0.0
      XCASBS(L,NY,NX)=0.0
      XMGOBS(L,NY,NX)=0.0
      XMGCBS(L,NY,NX)=0.0
      XMGHBS(L,NY,NX)=0.0
      XMGSBS(L,NY,NX)=0.0
      XNACBS(L,NY,NX)=0.0
      XNASBS(L,NY,NX)=0.0
      XKASBS(L,NY,NX)=0.0
      XH0PBS(L,NY,NX)=0.0
      XH3PBS(L,NY,NX)=0.0
      XF1PBS(L,NY,NX)=0.0
      XF2PBS(L,NY,NX)=0.0
      XC0PBS(L,NY,NX)=0.0
      XC1PBS(L,NY,NX)=0.0
      XC2PBS(L,NY,NX)=0.0
      XM1PBS(L,NY,NX)=0.0
      ENDIF
9865  CONTINUE
C
C     RESET ARRAYS TO TRANSFER MATERIALS WITHIN SOILS
C     AND BETWEEN SOILS AND PLANTS
C
      DO 9875 L=0,NL(NY,NX)
      DO 9950 K=0,1
      DO 9950 M=1,4
      CSNT(M,K,L,NY,NX)=0.0
      ZSNT(M,K,L,NY,NX)=0.0
      PSNT(M,K,L,NY,NX)=0.0
9950  CONTINUE
      DO 7775 K=0,4
      XOQCS(K,L,NY,NX)=0.0
      XOQNS(K,L,NY,NX)=0.0
      XOQPS(K,L,NY,NX)=0.0
      XOQAS(K,L,NY,NX)=0.0
7775  CONTINUE
      XZHYS(L,NY,NX)=0.0
      TRN4S(L,NY,NX)=0.0
      TRN3S(L,NY,NX)=0.0
      TRN3G(L,NY,NX)=0.0
      TRNO3(L,NY,NX)=0.0
      TRNO2(L,NY,NX)=0.0
      TRH1P(L,NY,NX)=0.0
      TRH2P(L,NY,NX)=0.0
      TRXN4(L,NY,NX)=0.0
      TRXH0(L,NY,NX)=0.0
      TRXH1(L,NY,NX)=0.0
      TRXH2(L,NY,NX)=0.0
      TRX1P(L,NY,NX)=0.0
      TRX2P(L,NY,NX)=0.0
      TRALPO(L,NY,NX)=0.0
      TRFEPO(L,NY,NX)=0.0
      TRCAPD(L,NY,NX)=0.0
      TRCAPH(L,NY,NX)=0.0
      TRCAPM(L,NY,NX)=0.0
      TUPWTR(L,NY,NX)=0.0
      TUPHT(L,NY,NX)=0.0
      XCODFG(L,NY,NX)=0.0
      XCHDFG(L,NY,NX)=0.0
      XOXDFG(L,NY,NX)=0.0
      XNGDFG(L,NY,NX)=0.0
      XN2DFG(L,NY,NX)=0.0
      XN3DFG(L,NY,NX)=0.0
      XNBDFG(L,NY,NX)=0.0
      XHGDFG(L,NY,NX)=0.0
      IF(L.GE.NU(NY,NX))THEN
      DO 195 K=0,4
      TDFOMC(K,L,NY,NX)=0.0
      TDFOMN(K,L,NY,NX)=0.0
      TDFOMP(K,L,NY,NX)=0.0
195   CONTINUE
      ENDIF
      DO 9795 M=1,NPH
      ROXSK(M,L,NY,NX)=0.0
9795  CONTINUE
C
C     IF SOC FLAG IS SET
C
      IF(IERSNG.EQ.2.OR.IERSNG.EQ.3)THEN
C
C     TOTAL SOC FOR CALCULATING CHANGES IN SOC CALCULATED IN NITRO.F
C    
C     OMC=microbial biomass, ORC=microbial residue
C     OQC,OQCH=DOC in micropores,macropores
C     OQA,OQAH=acetate in micropores,macropores
C     OHC,OHA=adsorbed SOC,acetate 
C     OSC=SOC(K=0:woody litter, K=1:non-woody litter,
C     K=2:manure, K=3:POC, K=4:humus)
C
      DC=0.0
      OC=0.0
      DO 7970 K=0,5
      DO 7950 N=1,7
      DO 7950 M=1,3
      OC=OC+OMC(M,N,K,L,NY,NX)
7950  CONTINUE
7970  CONTINUE
      DO 7900 K=0,4
      DO 7920 M=1,2
      OC=OC+ORC(M,K,L,NY,NX)
7920  CONTINUE
      OC=OC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
     2+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      DO 7910 M=1,4
      OC=OC+OSC(M,K,L,NY,NX)
7910  CONTINUE
7900  CONTINUE
      ORGCX(L,NY,NX)=OC
      ENDIF
9875  CONTINUE
      IFLGL=0
      IFLGY=0
      ICHKA=0
      DO 9985 L=NUI(NY,NX),NLI(NY,NX)
      FINH(L,NY,NX)=0.0
      TCO2S(L,NY,NX)=0.0
      TCO2P(L,NY,NX)=0.0
      TCOFLA(L,NY,NX)=0.0
      TCHFLA(L,NY,NX)=0.0
      TLCO2P(L,NY,NX)=0.0
      TUPOXP(L,NY,NX)=0.0
      TUPOXS(L,NY,NX)=0.0
      TUPCHS(L,NY,NX)=0.0
      TUPN2S(L,NY,NX)=0.0
      TUPN3S(L,NY,NX)=0.0
      TUPN3B(L,NY,NX)=0.0
      TUPHGS(L,NY,NX)=0.0
      TOXFLA(L,NY,NX)=0.0
      TCHFLA(L,NY,NX)=0.0
      TN2FLA(L,NY,NX)=0.0
      TNHFLA(L,NY,NX)=0.0
      THGFLA(L,NY,NX)=0.0
      TLOXYP(L,NY,NX)=0.0
      TLCH4P(L,NY,NX)=0.0
      TLN2OP(L,NY,NX)=0.0
      TLNH3P(L,NY,NX)=0.0
      TLH2GP(L,NY,NX)=0.0
      TUPNH4(L,NY,NX)=0.0
      TUPNO3(L,NY,NX)=0.0
      TUPH2P(L,NY,NX)=0.0
      TUPH1P(L,NY,NX)=0.0
      TUPNHB(L,NY,NX)=0.0
      TUPNOB(L,NY,NX)=0.0
      TUPH2B(L,NY,NX)=0.0
      TUPH1B(L,NY,NX)=0.0
      TUPNF(L,NY,NX)=0.0
      TRN4B(L,NY,NX)=0.0
      TRN3B(L,NY,NX)=0.0
      TRNOB(L,NY,NX)=0.0
      TRN2B(L,NY,NX)=0.0
      TRH1B(L,NY,NX)=0.0
      TRH2B(L,NY,NX)=0.0
      TRAL(L,NY,NX)=0.0
      TRFE(L,NY,NX)=0.0
      TRHY(L,NY,NX)=0.0
      TRCA(L,NY,NX)=0.0
      TRMG(L,NY,NX)=0.0
      TRNA(L,NY,NX)=0.0
      TRKA(L,NY,NX)=0.0
      TROH(L,NY,NX)=0.0
      TRSO4(L,NY,NX)=0.0
      TRCO3(L,NY,NX)=0.0
      TRHCO(L,NY,NX)=0.0
      TRCO2(L,NY,NX)=0.0
      TBCO2(L,NY,NX)=0.0
      TRAL1(L,NY,NX)=0.0
      TRAL2(L,NY,NX)=0.0
      TRAL3(L,NY,NX)=0.0
      TRAL4(L,NY,NX)=0.0
      TRALS(L,NY,NX)=0.0
      TRFE1(L,NY,NX)=0.0
      TRFE2(L,NY,NX)=0.0
      TRFE3(L,NY,NX)=0.0
      TRFE4(L,NY,NX)=0.0
      TRFES(L,NY,NX)=0.0
      TRCAO(L,NY,NX)=0.0
      TRCAC(L,NY,NX)=0.0
      TRCAH(L,NY,NX)=0.0
      TRCAS(L,NY,NX)=0.0
      TRMGO(L,NY,NX)=0.0
      TRMGC(L,NY,NX)=0.0
      TRMGH(L,NY,NX)=0.0
      TRMGS(L,NY,NX)=0.0
      TRNAC(L,NY,NX)=0.0
      TRNAS(L,NY,NX)=0.0
      TRKAS(L,NY,NX)=0.0
      TRH0P(L,NY,NX)=0.0
      TRH3P(L,NY,NX)=0.0
      TRC0P(L,NY,NX)=0.0
      TRF1P(L,NY,NX)=0.0
      TRF2P(L,NY,NX)=0.0
      TRC1P(L,NY,NX)=0.0
      TRC2P(L,NY,NX)=0.0
      TRM1P(L,NY,NX)=0.0
      TRH0B(L,NY,NX)=0.0
      TRH3B(L,NY,NX)=0.0
      TRF1B(L,NY,NX)=0.0
      TRF2B(L,NY,NX)=0.0
      TRC0B(L,NY,NX)=0.0
      TRC1B(L,NY,NX)=0.0
      TRC2B(L,NY,NX)=0.0
      TRM1B(L,NY,NX)=0.0
      TRXNB(L,NY,NX)=0.0
      TRXHY(L,NY,NX)=0.0
      TRXAL(L,NY,NX)=0.0
      TRXFE(L,NY,NX)=0.0
      TRXCA(L,NY,NX)=0.0
      TRXMG(L,NY,NX)=0.0
      TRXNA(L,NY,NX)=0.0
      TRXKA(L,NY,NX)=0.0
      TRXHC(L,NY,NX)=0.0
      TRXAL2(L,NY,NX)=0.0
      TRXFE2(L,NY,NX)=0.0
      TRBH0(L,NY,NX)=0.0
      TRBH1(L,NY,NX)=0.0
      TRBH2(L,NY,NX)=0.0
      TRB1P(L,NY,NX)=0.0
      TRB2P(L,NY,NX)=0.0
      TRALOH(L,NY,NX)=0.0
      TRFEOH(L,NY,NX)=0.0
      TRCACO(L,NY,NX)=0.0
      TRCASO(L,NY,NX)=0.0
      TRALPB(L,NY,NX)=0.0
      TRFEPB(L,NY,NX)=0.0
      TRCPDB(L,NY,NX)=0.0
      TRCPHB(L,NY,NX)=0.0
      TRCPMB(L,NY,NX)=0.0
      XCOFXS(L,NY,NX)=0.0
      XCHFXS(L,NY,NX)=0.0
      XOXFXS(L,NY,NX)=0.0
      XNGFXS(L,NY,NX)=0.0
      XN2FXS(L,NY,NX)=0.0
      XHGFXS(L,NY,NX)=0.0
      XN4FXW(L,NY,NX)=0.0
      XN3FXW(L,NY,NX)=0.0
      XNOFXW(L,NY,NX)=0.0
      XNXFXS(L,NY,NX)=0.0
      XH1PXS(L,NY,NX)=0.0
      XH2PXS(L,NY,NX)=0.0
      XN4FXB(L,NY,NX)=0.0
      XN3FXB(L,NY,NX)=0.0
      XNOFXB(L,NY,NX)=0.0
      XNXFXB(L,NY,NX)=0.0
      XH1BXB(L,NY,NX)=0.0
      XH2BXB(L,NY,NX)=0.0
      XALFXS(L,NY,NX)=0.0
      XFEFXS(L,NY,NX)=0.0
      XHYFXS(L,NY,NX)=0.0
      XCAFXS(L,NY,NX)=0.0
      XMGFXS(L,NY,NX)=0.0
      XNAFXS(L,NY,NX)=0.0
      XKAFXS(L,NY,NX)=0.0
      XOHFXS(L,NY,NX)=0.0
      XSOFXS(L,NY,NX)=0.0
      XCLFXS(L,NY,NX)=0.0
      XC3FXS(L,NY,NX)=0.0
      XHCFXS(L,NY,NX)=0.0
      XAL1XS(L,NY,NX)=0.0
      XAL2XS(L,NY,NX)=0.0
      XAL3XS(L,NY,NX)=0.0
      XAL4XS(L,NY,NX)=0.0
      XALSXS(L,NY,NX)=0.0
      XFE1XS(L,NY,NX)=0.0
      XFE2XS(L,NY,NX)=0.0
      XFE3XS(L,NY,NX)=0.0
      XFE4XS(L,NY,NX)=0.0
      XFESXS(L,NY,NX)=0.0
      XCAOXS(L,NY,NX)=0.0
      XCACXS(L,NY,NX)=0.0
      XCAHXS(L,NY,NX)=0.0
      XCASXS(L,NY,NX)=0.0
      XMGOXS(L,NY,NX)=0.0
      XMGCXS(L,NY,NX)=0.0
      XMGHXS(L,NY,NX)=0.0
      XMGSXS(L,NY,NX)=0.0
      XNACXS(L,NY,NX)=0.0
      XNASXS(L,NY,NX)=0.0
      XKASXS(L,NY,NX)=0.0
      XH0PXS(L,NY,NX)=0.0
      XH3PXS(L,NY,NX)=0.0
      XF1PXS(L,NY,NX)=0.0
      XF2PXS(L,NY,NX)=0.0
      XC0PXS(L,NY,NX)=0.0
      XC1PXS(L,NY,NX)=0.0
      XC2PXS(L,NY,NX)=0.0
      XM1PXS(L,NY,NX)=0.0
      XH0BXB(L,NY,NX)=0.0
      XH3BXB(L,NY,NX)=0.0
      XF1BXB(L,NY,NX)=0.0
      XF2BXB(L,NY,NX)=0.0
      XC0BXB(L,NY,NX)=0.0
      XC1BXB(L,NY,NX)=0.0
      XC2BXB(L,NY,NX)=0.0
      XM1BXB(L,NY,NX)=0.0
      DO 9955 K=0,4
      XOCFXS(K,L,NY,NX)=0.0
      XONFXS(K,L,NY,NX)=0.0
      XOPFXS(K,L,NY,NX)=0.0
      XOAFXS(K,L,NY,NX)=0.0
9955  CONTINUE
      THAW(L,NY,NX)=0.0
      THAWH(L,NY,NX)=0.0
      HTHAW(L,NY,NX)=0.0
      XCOBBL(L,NY,NX)=0.0
      XCHBBL(L,NY,NX)=0.0
      XOXBBL(L,NY,NX)=0.0
      XNGBBL(L,NY,NX)=0.0
      XN2BBL(L,NY,NX)=0.0
      XN3BBL(L,NY,NX)=0.0
      XNBBBL(L,NY,NX)=0.0
      XHGBBL(L,NY,NX)=0.0
      RTDNT(L,NY,NX)=0.0
C
C     CALCULATE SOIL CONCENTRATIONS OF SOLUTES, GASES
C
C     THETW,THETI,THETP=soil micropore water,ice,air concentration
C     THETPZ=soil micropore+macropore air concn for output
C
      IF(VOLX(L,NY,NX).LE.ZEROS(NY,NX))THEN
      THETW(L,NY,NX)=POROS(L,NY,NX) 
      THETI(L,NY,NX)=0.0
      THETP(L,NY,NX)=0.0
      ELSE
      THETW(L,NY,NX)=AMAX1(0.0,AMIN1(POROS(L,NY,NX) 
     2,VOLW(L,NY,NX)/VOLY(L,NY,NX)))
      THETI(L,NY,NX)=AMAX1(0.0,AMIN1(POROS(L,NY,NX)
     2,VOLI(L,NY,NX)/VOLY(L,NY,NX)))
      THETP(L,NY,NX)=AMAX1(0.0,VOLP(L,NY,NX)/VOLY(L,NY,NX))
      ENDIF
      THETPZ(L,NY,NX)=AMAX1(0.0,POROS(L,NY,NX)-THETW(L,NY,NX)
     2-THETI(L,NY,NX)) 
C     IF(L.EQ.7)THEN
C     WRITE(*,1117)'BKDS',I,J,L 
C    3,BKDS(L,NY,NX),BKDSI(L,NY,NX),BKVL(L,NY,NX)
C    4,DLYR(3,L,NY,NX),DLYRI(3,L,NY,NX)
C    2,CORGC(L,NY,NX),ORGC(L,NY,NX)
C    6,VOLT(L,NY,NX),VOLX(L,NY,NX),VOLA(L,NY,NX)
C    7,VOLY(L,NY,NX),VOLW(L,NY,NX),VOLI(L,NY,NX)
C    7,VOLWH(L,NY,NX),VOLIH(L,NY,NX)
C    8,THETWZ(L,NY,NX),THETIZ(L,NY,NX)
C    9,THETWZ(L,NY,NX)+THETIZ(L,NY,NX)
1117  FORMAT(A8,3I4,30E14.6)
C     ENDIF
C
C     GAS CONCENTRATIONS
C
C     C*G=soil gas gaseous concentration
C     C*S=soil gas aqueous concentration
C
      IF(THETP(L,NY,NX).GT.THETX)THEN
      CCO2G(L,NY,NX)=AMAX1(0.0,CO2G(L,NY,NX)/VOLP(L,NY,NX))
      CCH4G(L,NY,NX)=AMAX1(0.0,CH4G(L,NY,NX)/VOLP(L,NY,NX))
      COXYG(L,NY,NX)=AMAX1(0.0,OXYG(L,NY,NX)/VOLP(L,NY,NX))
      CZ2GG(L,NY,NX)=AMAX1(0.0,Z2GG(L,NY,NX)/VOLP(L,NY,NX))
      CZ2OG(L,NY,NX)=AMAX1(0.0,Z2OG(L,NY,NX)/VOLP(L,NY,NX))
      CNH3G(L,NY,NX)=AMAX1(0.0,ZNH3G(L,NY,NX)/VOLP(L,NY,NX))
      CH2GG(L,NY,NX)=AMAX1(0.0,H2GG(L,NY,NX)/VOLP(L,NY,NX))
      ELSE
      CCO2G(L,NY,NX)=0.0
      CCH4G(L,NY,NX)=0.0
      COXYG(L,NY,NX)=0.0
      CZ2GG(L,NY,NX)=0.0
      CZ2OG(L,NY,NX)=0.0
      CNH3G(L,NY,NX)=0.0
      CH2GG(L,NY,NX)=0.0
      ENDIF
      IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      CCO2S(L,NY,NX)=AMAX1(0.0,CO2S(L,NY,NX)/VOLW(L,NY,NX))
      CCH4S(L,NY,NX)=AMAX1(0.0,CH4S(L,NY,NX)/VOLW(L,NY,NX))
      COXYS(L,NY,NX)=AMAX1(0.0,OXYS(L,NY,NX)/VOLW(L,NY,NX))
      CZ2GS(L,NY,NX)=AMAX1(0.0,Z2GS(L,NY,NX)/VOLW(L,NY,NX))
      CZ2OS(L,NY,NX)=AMAX1(0.0,Z2OS(L,NY,NX)/VOLW(L,NY,NX))
      CH2GS(L,NY,NX)=AMAX1(0.0,H2GS(L,NY,NX)/VOLW(L,NY,NX))
      ELSE
      CCO2S(L,NY,NX)=0.0
      CCH4S(L,NY,NX)=0.0
      COXYS(L,NY,NX)=0.0
      CZ2GS(L,NY,NX)=0.0
      CZ2OS(L,NY,NX)=0.0
      CH2GS(L,NY,NX)=0.0
      ENDIF
C
C     CORGC=SOC concentration
C
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CORGC(L,NY,NX)=AMIN1(0.55E+06,ORGC(L,NY,NX)/BKVL(L,NY,NX))
      ELSE
      CORGC(L,NY,NX)=0.0
      ENDIF
C     WRITE(*,1113)'CORGC',I,J,NX,NY,L,CORGC(L,NY,NX)
C    2,ORGC(L,NY,NX),BKVL(L,NY,NX),BKDS(L,NY,NX),VOLX(L,NY,NX)
C
C     CALCULATE SOIL CONCENTRATIONS OF NH4, NH3, NO3, PO4
C     IN BAND AND NON-BAND ZONES
C
C     C*S=solute concentration in non-band
C     CH1P4,CH2P4=HPO4,H2PO4 concentration in non-band
C     Z*P=P ion pair amounts in non-band (see solute.f)
C     VLNH4,VLNO3,VLPO4=fraction of soil volume in NH4,NO3,PO4 non-band
C
      IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      IF(VLNH4(L,NY,NX).GT.ZERO)THEN
      CNH4S(L,NY,NX)=AMAX1(0.0,ZNH4S(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLNH4(L,NY,NX)))
      CNH3S(L,NY,NX)=AMAX1(0.0,ZNH3S(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLNH4(L,NY,NX)))
      ELSE
      CNH4S(L,NY,NX)=0.0
      CNH3S(L,NY,NX)=0.0
      ENDIF
      IF(VLNO3(L,NY,NX).GT.ZERO)THEN
      CNO3S(L,NY,NX)=AMAX1(0.0,ZNO3S(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLNO3(L,NY,NX)))
      CNO2S(L,NY,NX)=AMAX1(0.0,ZNO2S(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLNO3(L,NY,NX)))
      ELSE
      CNO3S(L,NY,NX)=0.0
      CNO2S(L,NY,NX)=0.0
      ENDIF
      IF(VLPO4(L,NY,NX).GT.ZERO)THEN
      CH1P4(L,NY,NX)=AMAX1(0.0,H1PO4(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLPO4(L,NY,NX)))
      CH2P4(L,NY,NX)=AMAX1(0.0,H2PO4(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLPO4(L,NY,NX)))
      CPO4S(L,NY,NX)=AMAX1(0.0,((H0PO4(L,NY,NX)+H3PO4(L,NY,NX)
     2+ZFE1P(L,NY,NX)+ZFE2P(L,NY,NX)+ZCA0P(L,NY,NX)
     3+ZCA1P(L,NY,NX)+ZCA2P(L,NY,NX)+ZMG1P(L,NY,NX))*31.0
     4+H1PO4(L,NY,NX)+H2PO4(L,NY,NX))/(VOLW(L,NY,NX)*VLPO4(L,NY,NX)))
      ELSE
      CH1P4(L,NY,NX)=0.0
      CH2P4(L,NY,NX)=0.0
      CPO4S(L,NY,NX)=0.0
      ENDIF
C
C     C*B=solute concentration in band
C     CH1PB,CH2PB=HPO4,H2PO4 concentration in band
C     Z*B=P ion pair amounts in band (see solute.f)
C     VLNHB,VLNOB,VLPOB=fraction of soil volume in NH4,NO3,PO4 band
C
      IF(VLNHB(L,NY,NX).GT.ZERO)THEN
      CNH4B(L,NY,NX)=AMAX1(0.0,ZNH4B(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLNHB(L,NY,NX)))
      CNH3B(L,NY,NX)=AMAX1(0.0,ZNH3B(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLNHB(L,NY,NX)))
      ELSE
      CNH4B(L,NY,NX)=0.0
      CNH3B(L,NY,NX)=0.0
      ENDIF
      IF(VLNOB(L,NY,NX).GT.ZERO)THEN
      CNO3B(L,NY,NX)=AMAX1(0.0,ZNO3B(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLNOB(L,NY,NX)))
      CNO2B(L,NY,NX)=AMAX1(0.0,ZNO2B(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLNOB(L,NY,NX)))
      ELSE
      CNO3B(L,NY,NX)=0.0
      CNO2B(L,NY,NX)=0.0
      ENDIF
      IF(VLPOB(L,NY,NX).GT.ZERO)THEN
      CH1P4B(L,NY,NX)=AMAX1(0.0,H1POB(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLPOB(L,NY,NX)))
      CH2P4B(L,NY,NX)=AMAX1(0.0,H2POB(L,NY,NX)
     2/(VOLW(L,NY,NX)*VLPOB(L,NY,NX)))
      CPO4B(L,NY,NX)=AMAX1(0.0,((H0POB(L,NY,NX)+H3POB(L,NY,NX) 
     2+ZFE1PB(L,NY,NX)+ZFE2PB(L,NY,NX)+ZCA0PB(L,NY,NX)
     3+ZCA1PB(L,NY,NX)+ZCA2PB(L,NY,NX)+ZMG1PB(L,NY,NX))*31.0
     4+H1POB(L,NY,NX)+H2POB(L,NY,NX))/(VOLW(L,NY,NX)*VLPOB(L,NY,NX)))
      ELSE
      CH1P4B(L,NY,NX)=0.0
      CH2P4B(L,NY,NX)=0.0
      CPO4B(L,NY,NX)=0.0
      ENDIF
      ELSE
      CNH4S(L,NY,NX)=0.0
      CNH3S(L,NY,NX)=0.0
      CNO3S(L,NY,NX)=0.0
      CNO2S(L,NY,NX)=0.0
      CH1P4(L,NY,NX)=0.0
      CH2P4(L,NY,NX)=0.0
      CPO4S(L,NY,NX)=0.0
      CNH4B(L,NY,NX)=0.0
      CNH3B(L,NY,NX)=0.0
      CNO3B(L,NY,NX)=0.0
      CNO2B(L,NY,NX)=0.0
      CH1P4B(L,NY,NX)=0.0
      CH2P4B(L,NY,NX)=0.0
      CPO4B(L,NY,NX)=0.0
      ENDIF
C
C     PREPARE ARRAYS FOR TOTAL O2 UPTAKE AND NH4,NO3.NO2,N2O,HPO4,H2PO4
C     UPTAKE IN NON-BAND,BAND AND DOC,DON,DOP,ACETATE UPTAKE
C
C     R*Y,R*X=total substrate uptake from previous,current hour 
C     used in nitro.f, uptake.f   
C
      ROXYY(L,NY,NX)=ROXYX(L,NY,NX)
      RNH4Y(L,NY,NX)=RNH4X(L,NY,NX)
      RNO3Y(L,NY,NX)=RNO3X(L,NY,NX)
      RNO2Y(L,NY,NX)=RNO2X(L,NY,NX)
      RN2OY(L,NY,NX)=RN2OX(L,NY,NX)
      RP14Y(L,NY,NX)=RP14X(L,NY,NX)
      RPO4Y(L,NY,NX)=RPO4X(L,NY,NX)
      RNHBY(L,NY,NX)=RNHBX(L,NY,NX)
      RN3BY(L,NY,NX)=RN3BX(L,NY,NX)
      RN2BY(L,NY,NX)=RN2BX(L,NY,NX)
      RP1BY(L,NY,NX)=RP1BX(L,NY,NX)
      RPOBY(L,NY,NX)=RPOBX(L,NY,NX)
      ROXYX(L,NY,NX)=0.0
      RNH4X(L,NY,NX)=0.0
      RNO3X(L,NY,NX)=0.0
      RNO2X(L,NY,NX)=0.0
      RN2OX(L,NY,NX)=0.0
      RP14X(L,NY,NX)=0.0
      RPO4X(L,NY,NX)=0.0
      RNHBX(L,NY,NX)=0.0
      RN3BX(L,NY,NX)=0.0
      RN2BX(L,NY,NX)=0.0
      RP1BX(L,NY,NX)=0.0
      RPOBX(L,NY,NX)=0.0
      DO 5050 K=0,4
      ROQCY(K,L,NY,NX)=ROQCX(K,L,NY,NX) 
      ROQAY(K,L,NY,NX)=ROQAX(K,L,NY,NX)
      ROQCX(K,L,NY,NX)=0.0 
      ROQAX(K,L,NY,NX)=0.0
5050  CONTINUE
C
C     DIFFUSIVITY
C
C     TFACG,TFACL=temperature effects on gaseous,aqueous diffusivity
C
C     *SGL= gaseous,aqueous diffusivity for gases,solutes listed in
C     *SG PARAMETER statement above
C
      TFACG=(TKS(L,NY,NX)/298.15)**1.75
      TFACL=(TKS(L,NY,NX)/298.15)**6
      TFND(L,NY,NX)=TFACL
      CGSGL(L,NY,NX)=CGSG*TFACG
      CHSGL(L,NY,NX)=CHSG*TFACG
      OGSGL(L,NY,NX)=OGSG*TFACG
      ZGSGL(L,NY,NX)=ZGSG*TFACG
      Z2SGL(L,NY,NX)=Z2SG*TFACG
      ZHSGL(L,NY,NX)=ZHSG*TFACG
      HGSGL(L,NY,NX)=HGSG*TFACG
      CLSGL(L,NY,NX)=CLSG*TFACL
      CQSGL(L,NY,NX)=CQSG*TFACL
      OLSGL(L,NY,NX)=OLSG*TFACL
      ZLSGL(L,NY,NX)=ZLSG*TFACL
      ZNSGL(L,NY,NX)=ZNSG*TFACL
      HLSGL(L,NY,NX)=HLSG*TFACL
      ZVSGL(L,NY,NX)=ZVSG*TFACL
      ZOSGL(L,NY,NX)=ZOSG*TFACL
      POSGL(L,NY,NX)=POSG*TFACL
      OCSGL(L,NY,NX)=OCSG*TFACL
      ONSGL(L,NY,NX)=ONSG*TFACL
      OPSGL(L,NY,NX)=OPSG*TFACL
      OASGL(L,NY,NX)=OASG*TFACL
      WGSGL(L,NY,NX)=WGSG*TFACG
      IF(ISALTG.NE.0)THEN
      ALSGL(L,NY,NX)=ALSG*TFACL
      FESGL(L,NY,NX)=FESG*TFACL
      HYSGL(L,NY,NX)=HYSG*TFACL
      CASGL(L,NY,NX)=CASG*TFACL
      GMSGL(L,NY,NX)=GMSG*TFACL
      ANSGL(L,NY,NX)=ANSG*TFACL
      AKSGL(L,NY,NX)=AKSG*TFACL
      OHSGL(L,NY,NX)=OHSG*TFACL
      C3SGL(L,NY,NX)=C3SG*TFACL
      HCSGL(L,NY,NX)=HCSG*TFACL
      SOSGL(L,NY,NX)=SOSG*TFACL
      CLSXL(L,NY,NX)=CLSX*TFACL
C
C     TOTAL ION CONCENTRATION
C
C     ZC3,ZA3,ZC2,ZA2,ZC1,ZA1=total tri-,di-,univalent cations C,anions A
C     CSTR,CION=ion strength, total ion concentration
C
      ZC3=ZAL(L,NY,NX)+ZFE(L,NY,NX)
      ZA3=H0PO4(L,NY,NX)+H0POB(L,NY,NX)
      ZC2=ZCA(L,NY,NX)+ZMG(L,NY,NX)+ZALOH1(L,NY,NX)+ZFEOH1(L,NY,NX)
     2+ZFE2P(L,NY,NX)+ZFE2PB(L,NY,NX)
      ZA2=ZSO4(L,NY,NX)+ZCO3(L,NY,NX)+H1PO4(L,NY,NX)+H1POB(L,NY,NX)
      ZC1=(ZNH4S(L,NY,NX)+ZNH4B(L,NY,NX))/14.0+ZHY(L,NY,NX)
     2+ZNA(L,NY,NX)+ZKA(L,NY,NX)+ZALOH2(L,NY,NX)+ZFEOH2(L,NY,NX)
     3+ZALS(L,NY,NX)+ZFES(L,NY,NX)+ZCAO(L,NY,NX)+ZCAH(L,NY,NX)
     4+ZMGO(L,NY,NX)+ZMGH(L,NY,NX)+ZFE1P(L,NY,NX)+ZFE1PB(L,NY,NX)
     5+ZCA2P(L,NY,NX)+ZCA2PB(L,NY,NX)
      ZA1=(ZNO3S(L,NY,NX)+ZNO3B(L,NY,NX))/14.0+ZOH(L,NY,NX)
     2+ZHCO3(L,NY,NX)+ZCL(L,NY,NX)+ZALOH4(L,NY,NX)+ZFEOH4(L,NY,NX)
     3+ZNAC(L,NY,NX)+ZNAS(L,NY,NX)+ZKAS(L,NY,NX)+(H2PO4(L,NY,NX)
     4+H2POB(L,NY,NX))/31.0+ZCA0P(L,NY,NX)+ZCA0PB(L,NY,NX)
      ZN=CO2S(L,NY,NX)/12.0+CH4S(L,NY,NX)/12.0+OXYS(L,NY,NX)/32.0
     2+(Z2GS(L,NY,NX)+Z2OS(L,NY,NX)+ZNH3S(L,NY,NX)+ZNH3B(L,NY,NX))/14.0
     3+ZALOH3(L,NY,NX)+ZFEOH3(L,NY,NX)+ZCAC(L,NY,NX)+ZCAS(L,NY,NX)
     4+ZMGC(L,NY,NX)+ZMGS(L,NY,NX)+H3PO4(L,NY,NX)+ZCA1P(L,NY,NX)
     5+ZMG1P(L,NY,NX)+H3POB(L,NY,NX)+ZCA1PB(L,NY,NX)+ZMG1PB(L,NY,NX)
      ZION1=ABS(3.0*(ZC3-ZA3)+2.0*(ZC2-ZA2)+ZC1-ZA1)
      IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      CSTR(L,NY,NX)=AMAX1(0.0,0.5E-03*(9.0*(ZC3+ZA3)+4.0*(ZC2+ZA2)
     2+ZC1+ZA1+ZION1)/VOLW(L,NY,NX))
      CION(L,NY,NX)=AMAX1(0.0,(ZC3+ZA3+ZC2+ZA2+ZC1+ZA1+ZN)
     2/VOLW(L,NY,NX))
      ELSE
      CSTR(L,NY,NX)=0.0
      CION(L,NY,NX)=0.0
      ENDIF
      ENDIF
C     IF(L.EQ.1)THEN
C     WRITE(*,1113)'CION',I,J,NX,NY,L,CION(L,NY,NX)
C    2,CSTR(L,NY,NX),ZC3,ZA3,ZC2,ZA2,ZC1,ZA1,ZN,VOLW(L,NY,NX)
C    3,ZAL(L,NY,NX),ZFE(L,NY,NX),ZNO3S(L,NY,NX),ZNO3B(L,NY,NX)
C    4,ZOH(L,NY,NX),ZHCO3(L,NY,NX),ZCL(L,NY,NX),ZALOH4(L,NY,NX)
C    5,ZFEOH4(L,NY,NX),ZNAC(L,NY,NX),ZNAS(L,NY,NX)
C    6,ZKAS(L,NY,NX),H2PO4(L,NY,NX),H2POB(L,NY,NX)
C    7,ZCA0P(L,NY,NX),ZCA0PB(L,NY,NX),ZCO3(L,NY,NX)
C     ENDIF
C
C     OSTWALD COEFFICIENTS FOR CO2, CH4, O2, N2, N2O, NH3 AND H2
C     SOLUBILITY IN WATER
C
C     S*L=solubility of gas in water
C     TCS=soil temperature (oC)
C
      FH2O=5.56E+04/(5.56E+04+CION(L,NY,NX))
      SCO2L(L,NY,NX)=SCO2X/(EXP(ACO2X*CSTR(L,NY,NX)))
     2*EXP(0.843-0.0281*TCS(L,NY,NX))*FH2O
      SCH4L(L,NY,NX)=SCH4X/(EXP(ACH4X*CSTR(L,NY,NX)))
     2*EXP(0.597-0.0199*TCS(L,NY,NX))*FH2O
      SOXYL(L,NY,NX)=SOXYX/(EXP(AOXYX*CSTR(L,NY,NX)))
     2*EXP(0.516-0.0172*TCS(L,NY,NX))*FH2O
      SN2GL(L,NY,NX)=SN2GX/(EXP(AN2GX*CSTR(L,NY,NX)))
     2*EXP(0.456-0.0152*TCS(L,NY,NX))*FH2O
      SN2OL(L,NY,NX)=SN2OX/(EXP(AN2OX*CSTR(L,NY,NX)))
     2*EXP(0.897-0.0299*TCS(L,NY,NX))*FH2O
      SNH3L(L,NY,NX)=SNH3X/(EXP(ANH3X*CSTR(L,NY,NX)))
     2*EXP(0.513-0.0171*TCS(L,NY,NX))*FH2O
      SH2GL(L,NY,NX)=SH2GX/(EXP(AH2GX*CSTR(L,NY,NX)))
     2*EXP(0.597-0.0199*TCS(L,NY,NX))*FH2O
C
C     WATER POTENTIALS
C
C     FC,WP=water contents at field capacity,wilting point,saturation
C     PSISM,PSISE=matric,saturation water potential
C     SRP=parameter for deviation from linear log-log water retention 
C     FC,WP=water contents at field capacity,wilting point from soil file
C     FCL,WPL=log FC,WP
C     FCD,PSD=FCL-WPL,log(POROS)-FCL
C     FCI,WPI=FC,WP of ice
C     THETIX=ice concentration
C
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX)
     2.AND.VOLX(L,NY,NX).GT.ZEROS(NY,NX))THEN
      THETW1=AMAX1(0.0,AMIN1(POROS(L,NY,NX) 
     2,VOLW(L,NY,NX)/VOLY(L,NY,NX)))
      IF(THETW1.LT.FC(L,NY,NX))THEN
      PSISM(L,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCL(L,NY,NX)-LOG(THETW1))
     3/FCD(L,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETW1.LT.POROS(L,NY,NX)-DTHETW)THEN 
      PSISM(L,NY,NX)=-EXP(PSIMS(NY,NX)
     2+(((PSL(L,NY,NX)-LOG(THETW1))
     3/PSD(L,NY,NX))**SRP(L,NY,NX)*PSISD(NY,NX)))
      ELSE
      PSISM(L,NY,NX)=PSISE(L,NY,NX)
      ENDIF
      ELSEIF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      FCX=FCI*THETI(L,NY,NX)
      WPX=WPI*THETI(L,NY,NX)
      FCLX=LOG(FCX)
      WPLX=LOG(WPX)
      PSDX=PSL(L,NY,NX)-FCLX
      FCDX=FCLX-WPLX
      IF(THETW(L,NY,NX).LT.FCX)THEN
      PSISM(L,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCLX-LOG(THETW(L,NY,NX)))
     3/FCDX*PSIMD(NY,NX))))
      ELSEIF(THETW(L,NY,NX).LT.POROS(L,NY,NX)-DTHETW)THEN 
      PSISM(L,NY,NX)=-EXP(PSIMS(NY,NX)
     2+(((PSL(L,NY,NX)-LOG(THETW(L,NY,NX)))
     3/PSDX)*PSISD(NY,NX)))
      ELSE
      PSISM(L,NY,NX)=PSISE(L,NY,NX)
      ENDIF
      ELSE
      PSISM(L,NY,NX)=PSISE(L,NY,NX)
      ENDIF
C     WRITE(*,443)'PSISM',I,J,NX,NY,L,PSISM(L,NY,NX)
C    2,THETW(L,NY,NX),THETI(L,NY,NX),FCX,WPX,POROS(L,NY,NX)
C
C     SOIL OSMOTIC, GRAVIMETRIC AND MATRIC WATER POTENTIALS
C
C     PSISM,PSISO,PSISH,PSIST=matric,osmotic,gravimetric,total water potential
C
      PSISO(L,NY,NX)=-8.3143E-06*TKS(L,NY,NX)*CION(L,NY,NX)
      PSISH(L,NY,NX)=0.0098*(ALT(NY,NX)-DPTH(L,NY,NX))
      PSIST(L,NY,NX)=AMIN1(0.0,PSISM(L,NY,NX)+PSISO(L,NY,NX)
     2+PSISH(L,NY,NX))
C     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,1113)'PSISM1',I,J,NX,NY,L,PSISM(L,NY,NX)
C    2,THETW(L,NY,NX),THETW1,VOLX(L,NY,NX) 
C    2,FC(L,NY,NX),WP(L,NY,NX),POROS(L,NY,NX),VOLP(L,NY,NX)
C    3,VOLW(L,NY,NX),VOLI(L,NY,NX),VOLA(L,NY,NX),VOLT(L,NY,NX)
C    4,CDPTH(L,NY,NX),DPTH(L,NY,NX),CDPTHZ(L,NY,NX),DPTHZ(L,NY,NX)
C    5,DLYR(3,L,NY,NX),PSIST(L,NY,NX),PSISO(L,NY,NX)
C    2,PSISH(L,NY,NX),TKS(L,NY,NX),CION(L,NY,NX)
1113  FORMAT(A8,5I4,50E12.4)
C     ENDIF
C
C     SOIL RESISTANCE TO ROOT PENETRATION
C
C     RSCS=soil resistance to root penetration (MPa)
C
C     IF(BKDS(L,NY,NX).GT.ZERO)THEN
C     CCLAYT=CCLAY(L,NY,NX)*1.0E+02
C     CORGCT=CORGC(L,NY,NX)*1.0E-04
C     CC=EXP(-3.6733-0.1447*CCLAYT+0.7653*CORGCT)
C     DD=-0.4805-0.1239*CCLAYT+0.2080*CORGCT
C     EE=3.8521+0.0963*CCLAYT
C     RSCS(L,NY,NX)=CC*THETW(L,NY,NX)**DD*BKDS(L,NY,NX)**EE
C     ELSE
      RSCS(L,NY,NX)=0.0
C     ENDIF
C     WRITE(*,2442)'RSCS',I,J,NX,NY,L,RSCS(L,NY,NX),THETW(L,NY,NX)
C    2,BKDS(L,NY,NX),CCLAY(L,NY,NX),CORGC(L,NY,NX)
2442  FORMAT(A8,5I4,12E12.4)
C
C     SOIL HYDRAULIC CONDUCTIVITIES FROM AMBIENT SOIL WATER CONTENTS
C
C     CNDU=soil hydraulic conductivity for root uptake
C
      K=MAX(1,MIN(100,INT(100.0*(POROS(L,NY,NX)-THETW(L,NY,NX))
     2/POROS(L,NY,NX))+1))
      CNDU(L,NY,NX)=0.5*(HCND(1,K,L,NY,NX)+HCND(3,K,L,NY,NX))
C
C     CALCULATE ACTIVE LAYER DEPTH
C
C     VOLI,VOLIH=ice volume in micropores,macropores 
C     VOLW,VOLWH=water volume in micropores,macropores 
C     VOLA,VOLAH=total volume in micropores,macropores
C     DPTHA=active layer depth
C     CDPTH,DLYR=depth to bottom,thickness of soil layer 
C
      IF(ICHKA.EQ.0)THEN
      VOLIT=VOLI(L,NY,NX)+VOLIH(L,NY,NX)
      VOLAT=VOLA(L,NY,NX)+VOLAH(L,NY,NX)
      IF(VOLAT.GT.ZEROS2(NY,NX).AND.VOLIT.GT.0.01*VOLAT)THEN
      DO 5700 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
      VOLITL=VOLI(LL,NY,NX)+VOLIH(LL,NY,NX)
      VOLWTL=VOLW(LL,NY,NX)+VOLWH(LL,NY,NX)
      VOLATL=VOLA(LL,NY,NX)+VOLAH(LL,NY,NX)
      IF(VOLATL.GT.ZEROS2(NY,NX).AND.VOLITL.LT.0.01*VOLATL)THEN
      GO TO 5701
      ENDIF
5700  CONTINUE
      IF(VOLAT.GT.ZEROS2(NY,NX))THEN
      DPTHA(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)
     2*AMIN1(1.0,VOLIT/VOLAT)
      ELSE
      DPTHA(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)
      ENDIF
      ICHKA=1
      GO TO 5702
5701  DPTHA(NY,NX)=9999.0
5702  CONTINUE
      ENDIF
      ENDIF
C
C     OUTPUT FOR WATER TABLE DEPTH
C
C     IDTBL=water table flag from site file
C     THETPZ,THETPW=current,minimum air-filled, porosity for water table
C     DPTH,DTBLX=depth of soil layer midpoint, water table
C     PSIS1=water potential in hydraulic equilibrium with layer below
C     THETW1,THETWP=water content at PSIS1,minimum SWC for water table
C     DPTHT=water table depth 
C
      IF(IDTBL(NY,NX).NE.0)THEN
      IF(IFLGY.EQ.0)THEN
      IF(THETPZ(L,NY,NX).LT.THETPW.OR.L.EQ.NL(NY,NX))THEN
      IFLGY=1
      IF(DPTH(L,NY,NX).LT.DTBLX(NY,NX))THEN
      DO 5705 LL=MIN(L+1,NL(NY,NX)),NL(NY,NX)
      IF(THETPZ(LL,NY,NX).GE.THETPW.AND.LL.NE.NL(NY,NX))THEN
      IFLGY=0
      GO TO 5706
      ELSEIF(DPTH(LL,NY,NX).GE.DTBLX(NY,NX))THEN
      GO TO 5706
      ENDIF
5705  CONTINUE
      ENDIF
5706  CONTINUE
      IF(IFLGY.EQ.1)THEN
      IF(THETPZ(L,NY,NX).GE.THETPW.AND.L.NE.NL(NY,NX))THEN
      PSIS1=PSISM(L+1,NY,NX)-0.0098*(DPTH(L+1,NY,NX)-DPTH(L,NY,NX))
      THETWM=THETWP*POROS(L,NY,NX)
      THETW1=AMIN1(THETWM,EXP((PSIMS(NY,NX)-LOG(-PSIS1))
     2*PSD(L,NY,NX)/PSISD(NY,NX)+PSL(L,NY,NX)))
      IF(THETWM.GT.THETW1)THEN
      THETPX=AMIN1(1.0,AMAX1(0.0,(THETWM-THETW(L,NY,NX))
     2/(THETWM-THETW1)))
      DPTHT(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)*(1.0-THETPX)
      ELSE
      DPTHT(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX)
      ENDIF
      ELSEIF(L.GT.NU(NY,NX))THEN
      PSIS1=PSISM(L,NY,NX)-0.0098*(DPTH(L,NY,NX)-DPTH(L-1,NY,NX))
      THETWM=THETWP*POROS(L-1,NY,NX)
      THETW1=AMIN1(THETWM,EXP((PSIMS(NY,NX)-LOG(-PSIS1))
     2*PSD(L-1,NY,NX)/PSISD(NY,NX)+PSL(L-1,NY,NX)))
      IF(THETWM.GT.THETW1)THEN
      THETPX=AMIN1(1.0,AMAX1(0.0,(THETWM-THETW(L-1,NY,NX))
     2/(THETWM-THETW1)))
      DPTHT(NY,NX)=CDPTH(L-1,NY,NX)-DLYR(3,L-1,NY,NX)*(1.0-THETPX)
      ELSE
      DPTHT(NY,NX)=CDPTH(L-1,NY,NX)-DLYR(3,L-1,NY,NX) 
      ENDIF
      ELSE
      DPTHT(NY,NX)=CDPTH(L,NY,NX)-DLYR(3,L,NY,NX) 
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C     IF(NX.EQ.3.AND.NY.EQ.3)THEN
C     WRITE(*,5353)'DPTHT',I,J,NX,NY,L,LL,IFLGY,DPTHT(NY,NX)
C    2,CDPTH(L,NY,NX),DLYR(3,L,NY,NX),DPTH(L,NY,NX),DPTH(L-1,NY,NX)
C    3,PSIS1,PSISM(L,NY,NX),THETWM,THETW1,THETW(L-1,NY,NX)
C    4,THETPZ(L,NY,NX),THETPW 
5353  FORMAT(A8,7I4,30E12.4)
C     ENDIF
      ENDIF
9985  CONTINUE
C
C     PHYSICAL PROPERTIES, AND WATER, GAS, AND MINERAL CONTENTS
C     OF SURFACE RESIDUE
C
C     VOLWRX=liter water holding capacity
C     THETRX,RC0=specific WHC,mass of woody(0),fine(1),manure(2) litter 
C     VOLR=dry litter volume
C     BKRS=dry bulk density of woody(0),fine(1),manure(2) litter
C     TVOLG0=excess litter water+ice
C     VOLT,VOLX=wet litter volume
C     BKVL=litter mass
C     VOLW,VOLI,VOLA,VOLP=litter water,ice,porosity,air volume
C     THETW,THETI,THETA,THETP=litter water,ice,porosity,air concentration
C     POROS=litter porosity
C     THETW0,THETI0,DPTH0=litter excess water,ice,water+ice depth
C     DLYR=litter thickness
C     PSISM,PSISE=litter matric,saturation water potential
C
      TVOLG0=AMAX1(0.0,VOLW(0,NY,NX)+VOLI(0,NY,NX)-VOLWRX(NY,NX))
      VOLT(0,NY,NX)=TVOLG0+VOLR(NY,NX)
      IF(VOLT(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      VOLX(0,NY,NX)=VOLT(0,NY,NX)
      BKVL(0,NY,NX)=1.82E-06*ORGC(0,NY,NX)
      VOLA(0,NY,NX)=AMAX1(0.0,VOLR(NY,NX)-BKVL(0,NY,NX)/1.30)
      VOLP(0,NY,NX)=AMAX1(0.0,VOLA(0,NY,NX)-VOLW(0,NY,NX)
     2-VOLI(0,NY,NX))
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX))THEN
      POROS(0,NY,NX)=VOLA(0,NY,NX)/VOLR(NY,NX)
      THETW(0,NY,NX)=AMAX1(0.0,AMIN1(1.0
     2,VOLW(0,NY,NX)/VOLR(NY,NX)))
      THETI(0,NY,NX)=AMAX1(0.0,AMIN1(1.0
     2,VOLI(0,NY,NX)/VOLR(NY,NX)))
      THETP(0,NY,NX)=AMAX1(0.0,AMIN1(1.0
     2,VOLP(0,NY,NX)/VOLR(NY,NX)))
      ELSE
      POROS(0,NY,NX)=1.0
      THETW(0,NY,NX)=0.0
      THETI(0,NY,NX)=0.0
      THETP(0,NY,NX)=0.0
      ENDIF
      TVOLWI=VOLW(0,NY,NX)+VOLI(0,NY,NX)
      IF(TVOLWI.GT.ZEROS(NY,NX))THEN
      VOLWRZ=VOLW(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      VOLIRZ=VOLI(0,NY,NX)/TVOLWI*VOLWRX(NY,NX)
      XVOLW0=AMAX1(0.0,VOLW(0,NY,NX)-VOLWRZ)/AREA(3,NU(NY,NX),NY,NX)
      XVOLI0=AMAX1(0.0,VOLI(0,NY,NX)-VOLIRZ)/AREA(3,NU(NY,NX),NY,NX)
      ELSE
      XVOLW0=0.0
      XVOLI0=0.0
      ENDIF
      DPTH0(NY,NX)=XVOLW0+XVOLI0
      DLYR(3,0,NY,NX)=VOLX(0,NY,NX)/AREA(3,0,NY,NX)
      IF(VOLR(NY,NX).GT.ZEROS(NY,NX)
     2.AND.VOLW(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      THETWR=AMIN1(VOLWRX(NY,NX),VOLW(0,NY,NX))/VOLR(NY,NX)
      IF(THETWR.LT.FC(0,NY,NX))THEN
      PSISM(0,NY,NX)=AMAX1(PSIHY,-EXP(PSIMX(NY,NX)
     2+((FCL(0,NY,NX)-LOG(THETWR))
     3/FCD(0,NY,NX)*PSIMD(NY,NX))))
      ELSEIF(THETWR.LT.POROS(0,NY,NX))THEN 
      PSISM(0,NY,NX)=-EXP(PSIMS(NY,NX)
     2+(((PSL(0,NY,NX)-LOG(THETWR))
     3/PSD(0,NY,NX))**SRP(0,NY,NX)*PSISD(NY,NX)))
      ELSE
      PSISM(0,NY,NX)=PSISE(0,NY,NX)
      ENDIF
      PSISO(0,NY,NX)=0.0
      PSISH(0,NY,NX)=0.0098*(ALT(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX)
     2+0.5*DLYR(3,0,NY,NX))
      PSIST(0,NY,NX)=AMIN1(0.0,PSISM(0,NY,NX)+PSISO(0,NY,NX)
     2+PSISH(0,NY,NX))
C
C     LITTER NH4,NH3,NO3,NO2,HPO4,H2PO4 CONCENTRATIONS
C
C     C*=litter solute concentrations
C 
      CNH4S(0,NY,NX)=AMAX1(0.0,ZNH4S(0,NY,NX)/VOLW(0,NY,NX))
      CNH3S(0,NY,NX)=AMAX1(0.0,ZNH3S(0,NY,NX)/VOLW(0,NY,NX))
      CNO3S(0,NY,NX)=AMAX1(0.0,ZNO3S(0,NY,NX)/VOLW(0,NY,NX))
      CNO2S(0,NY,NX)=AMAX1(0.0,ZNO2S(0,NY,NX)/VOLW(0,NY,NX))
      CH1P4(0,NY,NX)=AMAX1(0.0,H1PO4(0,NY,NX)/VOLW(0,NY,NX))
      CH2P4(0,NY,NX)=AMAX1(0.0,H2PO4(0,NY,NX)/VOLW(0,NY,NX))
      ELSE
      PSISM(0,NY,NX)=PSISM(NU(NY,NX),NY,NX)
      CNH4S(0,NY,NX)=0.0
      CNH3S(0,NY,NX)=0.0
      CNO3S(0,NY,NX)=0.0
      CNO2S(0,NY,NX)=0.0
      CH1P4(0,NY,NX)=0.0
      CH2P4(0,NY,NX)=0.0
      ENDIF
      ELSE
      VOLX(0,NY,NX)=0.0
      BKVL(0,NY,NX)=0.0
      VOLA(0,NY,NX)=0.0
      VOLP(0,NY,NX)=0.0
      POROS(0,NY,NX)=1.0
      DLYR(3,0,NY,NX)=0.0
      THETW(0,NY,NX)=0.0
      THETI(0,NY,NX)=0.0
      THETP(0,NY,NX)=1.0
      VOLWRX(NY,NX)=0.0
      PSISM(0,NY,NX)=PSISM(NU(NY,NX),NY,NX)
      CNH4S(0,NY,NX)=0.0
      CNH3S(0,NY,NX)=0.0
      CNO3S(0,NY,NX)=0.0
      CNO2S(0,NY,NX)=0.0
      CH1P4(0,NY,NX)=0.0
      CH2P4(0,NY,NX)=0.0
      CCO2S(0,NY,NX)=0.0
      CCH4S(0,NY,NX)=0.0
      COXYS(0,NY,NX)=0.0
      CZ2GS(0,NY,NX)=0.0
      CZ2OS(0,NY,NX)=0.0
      CH2GS(0,NY,NX)=0.0
      ENDIF
C
C     LITTER GAS CONCENTRATIOS 
C    
C     C*G=soil gas gaseous concentration
C     *E=atmospheric concentration
C     TKS,TCS=litter temperature (K,C)
C     S*L=gas solubility
C     C*S=soil gas aqueous concentration
C     
      CCO2G(0,NY,NX)=CO2E(NY,NX)*5.36E-04*273.15/TKS(0,NY,NX)
      CCH4G(0,NY,NX)=CH4E(NY,NX)*5.36E-04*273.15/TKS(0,NY,NX)
      COXYG(0,NY,NX)=OXYE(NY,NX)*1.43E-03*273.15/TKS(0,NY,NX)
      CZ2GG(0,NY,NX)=Z2GE(NY,NX)*1.25E-03*273.15/TKS(0,NY,NX)
      CZ2OG(0,NY,NX)=Z2OE(NY,NX)*1.25E-03*273.15/TKS(0,NY,NX)
      CNH3G(0,NY,NX)=ZNH3E(NY,NX)*6.25E-04*273.15/TKS(0,NY,NX)
      CH2GG(0,NY,NX)=H2GE(NY,NX)*8.92E-05*273.15/TKS(0,NY,NX)
      CNH4B(0,NY,NX)=0.0
      CNH3B(0,NY,NX)=0.0
      CNO3B(0,NY,NX)=0.0
      CNO2B(0,NY,NX)=0.0
      CH2P4B(0,NY,NX)=0.0
      SCO2L(0,NY,NX)=SCO2X*EXP(0.843-0.0281*TCS(0,NY,NX)) 
      SCH4L(0,NY,NX)=SCH4X*EXP(0.597-0.0199*TCS(0,NY,NX))
      SOXYL(0,NY,NX)=SOXYX*EXP(0.516-0.0172*TCS(0,NY,NX)) 
      SN2GL(0,NY,NX)=SN2GX*EXP(0.456-0.0152*TCS(0,NY,NX)) 
      SN2OL(0,NY,NX)=SN2OX*EXP(0.897-0.0299*TCS(0,NY,NX)) 
      SNH3L(0,NY,NX)=SNH3X*EXP(0.513-0.0171*TCS(0,NY,NX)) 
      SH2GL(0,NY,NX)=SH2GX*EXP(0.597-0.0199*TCS(0,NY,NX)) 
      IF(VOLW(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      CCO2S(0,NY,NX)=AMAX1(0.0,CO2S(0,NY,NX)/VOLW(0,NY,NX))
      CCH4S(0,NY,NX)=AMAX1(0.0,CH4S(0,NY,NX)/VOLW(0,NY,NX))
      COXYS(0,NY,NX)=AMAX1(0.0,OXYS(0,NY,NX)/VOLW(0,NY,NX))
      CZ2GS(0,NY,NX)=AMAX1(0.0,Z2GS(0,NY,NX)/VOLW(0,NY,NX))
      CZ2OS(0,NY,NX)=AMAX1(0.0,Z2OS(0,NY,NX)/VOLW(0,NY,NX))
      CH2GS(0,NY,NX)=AMAX1(0.0,H2GS(0,NY,NX)/VOLW(0,NY,NX))
      ELSE
      CCO2S(0,NY,NX)=0.0
      CCH4S(0,NY,NX)=0.0 
      COXYS(0,NY,NX)=0.0
      CZ2GS(0,NY,NX)=0.0
      CZ2OS(0,NY,NX)=0.0
      CH2GS(0,NY,NX)=0.0
      ENDIF
C
C     TFACL=temperature effect on diffusivity
C     *SGL= gaseous,aqueous diffusivity for gases,solutes listed in
C     *SG PARAMETER statement above
C
      TFACL=(TKS(0,NY,NX)/298.15)**6
      TFND(0,NY,NX)=TFACL
      CQSGL(0,NY,NX)=CQSG*TFACL
      OLSGL(0,NY,NX)=OLSG*TFACL
      ZLSGL(0,NY,NX)=ZLSG*TFACL
      ZNSGL(0,NY,NX)=ZNSG*TFACL
      HLSGL(0,NY,NX)=HLSG*TFACL
      ZVSGL(0,NY,NX)=ZVSG*TFACL
      ZOSGL(0,NY,NX)=ZOSG*TFACL
      POSGL(0,NY,NX)=POSG*TFACL
      OCSGL(0,NY,NX)=OCSG*TFACL
      ONSGL(0,NY,NX)=ONSG*TFACL
      OPSGL(0,NY,NX)=OPSG*TFACL
      OASGL(0,NY,NX)=OASG*TFACL
C
C     R*Y,R*X=total substrate uptake from previous,current hour    
C     used in nitro.f, uptake.f   
C
      ROXYY(0,NY,NX)=ROXYX(0,NY,NX)
      RNH4Y(0,NY,NX)=RNH4X(0,NY,NX)
      RNO3Y(0,NY,NX)=RNO3X(0,NY,NX)
      RNO2Y(0,NY,NX)=RNO2X(0,NY,NX)
      RN2OY(0,NY,NX)=RN2OX(0,NY,NX)
      RP14Y(0,NY,NX)=RP14X(0,NY,NX)
      RPO4Y(0,NY,NX)=RPO4X(0,NY,NX)
      ROXYX(0,NY,NX)=0.0
      RNH4X(0,NY,NX)=0.0
      RNO3X(0,NY,NX)=0.0
      RNO2X(0,NY,NX)=0.0
      RN2OX(0,NY,NX)=0.0
      RP14X(0,NY,NX)=0.0
      RPO4X(0,NY,NX)=0.0
      DO 5055 K=0,4
      ROQCY(K,0,NY,NX)=ROQCX(K,0,NY,NX) 
      ROQAY(K,0,NY,NX)=ROQAX(K,0,NY,NX)
      ROQCX(K,0,NY,NX)=0.0 
      ROQAX(K,0,NY,NX)=0.0
5055  CONTINUE
C
C     WGSGA,WGSGR,WGSGW=vapor diffusivity in air,litter,snowpack
C
      TFACA=(TKA(NY,NX)/298.15)**1.75
      WGSGA(NY,NX)=WGSG*TFACA
      TFACR=(TKS(0,NY,NX)/298.15)**1.75
      WGSGR(NY,NX)=WGSG*TFACR
      DO 5060 L=1,JS
      TFACW=(TKW(L,NY,NX)/298.15)**1.75
      WGSGW(L,NY,NX)=WGSG*TFACW
5060  CONTINUE
C
C     MULTILAYER CANOPY INTERECEPTION OF DIRECT AND DIFFUSE RADIATION
C     IN SW AND VISIBLE BANDS BY INCLINATION N, AZIMUTH M, LAYER L,
C     NODE K, BRANCH NB, PFT NZ
C
C     ARLFS,ARLSS=leaf+stalk area of combined,each PFT canopy 
C     ZL=height to bottom of canopy layer
C     DPTHS,DPTH0=snowpack,surface water depths
C     ARLFL,ARSTK=leaf,stalk areas of PFT
C     RAD,RAP=vertical direct+diffuse SW,PAR
C     RADS,RADY,RAPS,RAPY=solar beam direct,diffuse SW,PAR 
C     SSIN,TYSIN=sine of solar,sky angles
C     RADC,RADP=total SW,PAR absorbed by canopy
C     CFX=clumping factor for self-shading
C
      ARLSS(NY,NX)=0.0
      DO 1135 NZ=1,NP(NY,NX)
      ARLFS(NZ,NY,NX)=0.0
      DO 1135 NB=1,NBR(NZ,NY,NX)
      DO 1135 L=1,JC
      IF(ZL(L-1,NY,NX).GE.DPTHS(NY,NX)-ZERO
     2.AND.ZL(L-1,NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
      DO 1130 K=1,25
      ARLFS(NZ,NY,NX)=ARLFS(NZ,NY,NX)+ARLFL(L,K,NB,NZ,NY,NX)
      ARLSS(NY,NX)=ARLSS(NY,NX)+ARLFL(L,K,NB,NZ,NY,NX)
1130  CONTINUE
      ARLFS(NZ,NY,NX)=ARLFS(NZ,NY,NX)+ARSTK(L,NB,NZ,NY,NX)
      ARLSS(NY,NX)=ARLSS(NY,NX)+ARSTK(L,NB,NZ,NY,NX)
      ENDIF
1135  CONTINUE
      IF(SSIN(NY,NX).GT.ZERO)THEN
      RAD(NY,NX)=RADS(NY,NX)*SSIN(NY,NX)+RADY(NY,NX)*TYSIN
      RAP(NY,NX)=RAPS(NY,NX)*SSIN(NY,NX)+RAPY(NY,NX)*TYSIN
      ELSE
      RADS(NY,NX)=0.0
      RADY(NY,NX)=0.0
      RAPS(NY,NX)=0.0
      RAPY(NY,NX)=0.0
      RAD(NY,NX)=0.0
      RAP(NY,NX)=0.0
      ENDIF
      TRADC(NY,NX)=0.0
      TRAPC(NY,NX)=0.0
      DO 1025 NZ=1,NP(NY,NX)
      RADC(NZ,NY,NX)=0.0
      RADP(NZ,NY,NX)=0.0
      CFX(NZ,NY,NX)=CF(NZ,NY,NX)*(1.0-0.025
     2*ARLFP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX))
1025  CONTINUE
C
C     ANGLE BETWEEN SUN AND GROUND SURFACE
C
C     SAZI,SCOS=solar azimuth,cosine of solar angle
C     BETAG=incident solar angle at ground surface
C     GCOS,GSIN=cos,sin of ground surface
C     ZNOON=hour of solar noon from weather file
C
      IF(SSIN(NY,NX).GT.ZERO)THEN
      SAZI=0.2618*(ZNOON(NY,NX)-J)+4.7124
      SCOS=SQRT(1.0-SSIN(NY,NX)**2)
      DGAZI=COS(GAZI(NY,NX)-SAZI)
      BETAG=AMAX1(0.0,AMIN1(1.0,GCOS(NY,NX)*SSIN(NY,NX)
     2+GSIN(NY,NX)*SCOS*DGAZI))
      IF(ARLSS(NY,NX).GT.0.0)THEN
      SAGL=ASIN(SSIN(NY,NX))
C
C     ABSORBED RADIATION FROM OPTICAL PROPERTIES ENTERED IN 'READS'
C
C     RADSA,RADWA,RAPSA,RAPWA=SW,PAR absorbed at leaf,stalk surface    
C     perpendicular to incoming radiation
C
      DO 1050 NZ=1,NP(NY,NX)
      RADSA(NZ,NY,NX)=RADS(NY,NX)*ABSR(NZ,NY,NX)
      RADWA(NZ,NY,NX)=RADS(NY,NX)*ABSRW
      RAPSA(NZ,NY,NX)=RAPS(NY,NX)*ABSP(NZ,NY,NX)
      RAPWA(NZ,NY,NX)=RAPS(NY,NX)*ABSPW
1050  CONTINUE
C
C     ANGLES BETWEEN SUN OR SKY ZONES AND FOLIAR SURFACES
C
C     ZAZI=leaf azimuth
C     BETA,BETX=incident angle of direct radiation at leaf,horizontal surface
C     ZAGL=determines forward vs backscattering
C     IALBS=flag for forward vs backscattering
C
      DO 1100 M=1,4
      ZAZI=SAZI+(M-0.5)*3.1416/4
      DAZI=COS(ZAZI-SAZI)
      DO 1100 N=1,4
      BETY=ZCOS(N)*SSIN(NY,NX)+ZSIN(N)*SCOS*DAZI
      BETA(N,M)=ABS(BETY)
      BETX(N,M)=BETA(N,M)/SSIN(NY,NX)
      IF(ZCOS(N).GT.SSIN(NY,NX))THEN
      BETZ=ACOS(BETY)
      ELSE
      BETZ=-ACOS(BETY)
      ENDIF
      IF(BETZ.GT.-1.5708)THEN
      ZAGL=SAGL+2.0*BETZ
      ELSE
      ZAGL=SAGL-2.0*(3.1416+BETZ)
      ENDIF
      IF(ZAGL.GT.0.0.AND.ZAGL.LT.3.1416)THEN
      IALBS(N,M)=1
      ELSE
      IALBS(N,M)=2
      ENDIF
C
C     INTENSITY OF ABSORBED DIRECT RADIATION AT LEAF SURFACES
C
C     RDNDIR,RDNDIW,PARDIR,PARDIW=SW,PAR flux absorbed by leaf,stalk surfaces
C     PAR,PARDIF=direct,diffuse PAR flux
C     RADYL,RAPYL=solar beam diffuse SW,PAR flux
C     RAFYL,RAFPL=forward scattered diffuse SW,PAR flux
C     TAUS,TAUY=fraction of direct,diffuse radiation transmitted 
C
      DO 1100 NZ=1,NP(NY,NX)
      RDNDIR(N,M,NZ,NY,NX)=RADSA(NZ,NY,NX)*ABS(BETA(N,M))
      RDNDIW(N,M,NZ,NY,NX)=RADWA(NZ,NY,NX)*ABS(BETA(N,M))
      PARDIR(N,M,NZ,NY,NX)=RAPSA(NZ,NY,NX)*ABS(BETA(N,M))
      PARDIW(N,M,NZ,NY,NX)=RAPWA(NZ,NY,NX)*ABS(BETA(N,M))
      DO 1100 L=1,JC
      PARDIF(N,M,L,NZ,NY,NX)=0.0
      PAR(N,M,L,NZ,NY,NX)=PARDIR(N,M,NZ,NY,NX)
1100  CONTINUE
      XAREA=1.00/AREA(3,NU(NY,NX),NY,NX)
      YAREA=0.25/AREA(3,NU(NY,NX),NY,NX)
      RADYL=RADY(NY,NX)
      RAPYL=RAPY(NY,NX)
      TAUS(JC+1,NY,NX)=1.0
      TAUY(JC+1)=1.0
      RAFSL(JC+1)=0.0
      RAFPL(JC+1)=0.0
      STOPS=0.0
C
C     RESET ARRAYS OF SUNLIT AND SHADED LEAF AREAS IN DIFFERENT
C     LAYERS AND ANGLE CLASSES
C
C     TSURF,TSURFB,SURF,SURFB=leaf,stalk total,PFT surface area
C
      DO 1150 NZ=1,NP(NY,NX)
      DO 1150 L=1,JC
      DO 1150 N=1,4
      TSURF(N,L,NZ,NY,NX)=0.0
      TSURFB(N,L,NZ,NY,NX)=0.0
1150  CONTINUE
      DO 1200 NZ=1,NP(NY,NX)
      DO 1200 NB=1,NBR(NZ,NY,NX)
      DO 1200 L=1,JC
      IF(ZL(L-1,NY,NX).GT.DPTHS(NY,NX)-ZERO
     2.AND.ZL(L-1,NY,NX).GT.DPTH0(NY,NX)-ZERO)THEN
      DO 1205 N=1,4
      DO 1210 K=1,25
      TSURF(N,L,NZ,NY,NX)=TSURF(N,L,NZ,NY,NX)+SURF(N,L,K,NB,NZ,NY,NX)
1210  CONTINUE
      TSURFB(N,L,NZ,NY,NX)=TSURFB(N,L,NZ,NY,NX)
     2+SURFB(N,L,NB,NZ,NY,NX)
1205  CONTINUE
      ENDIF
1200  CONTINUE
C
C     CALCULATE ABSORPTION, REFLECTION AND TRANSMISSION OF DIRECT AND 
C     DIFFUSE DOWNWARD TOTAL AND VISIBLE RADIATION BY EACH SPECIES 
C     NZ IN EACH LAYER L
C
C     RAFYL,RAFPL=forward scattered diffuse SW,PAR
C     RABYL,RABPL=backscattered diffuse SW,PAR
C     RADYL,RAPYL=solar beam diffuse SW,PAR
C     STOPY,STOPSZ,STOPYZ=fraction of direct,diffuse radiation intercepted
C
      DO 1800 L=JC,1,-1
      IF(ZL(L-1,NY,NX).GE.DPTHS(NY,NX)-ZERO
     2.AND.ZL(L-1,NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
      RADYL=RADYL*TAUY(L+1)+RAFSL(L+1)
      RAPYL=RAPYL*TAUY(L+1)+RAFPL(L+1)
      RAFSL(L)=0.0
      RAFPL(L)=0.0
      RABSL(L)=0.0
      RABPL(L)=0.0
      STOPY=0.0
      STOPSZ=0.0
      STOPYZ=0.0
C
C     RESET ACCUMULATORS OB ABSORBED, REFLECTED AND TRANSMITTED RADIATION
C
C     RADSL,RADSW,RADPL,RADPW=direct atmosph SW,PAR absbd by leaf,stalk surf
C     RAYSL,RAYSW,RAYPL,RAYPW=diffuse atmosph SW,PAR absbd by leaf,stalk surf
C     RADS1,RADW1,RADP1,RADQ1=backscattered direct SW,PAR absbd by leaf,stalk surf
C     RAYS1,RAYW1,RAYP1,RAYQ1=backscattered diffuse SW,PAR absbd by leaf,stalk surf
C     RADS2,RADW2,RADP2,RADQ2=fwd scattered direct SW,PAR absbd by leaf,stalk surf
C     RAYS2,RAYW2,RAYP2,RAYQ2=fwd scattered diffuse SW,PAR absbd by leaf,stalk surf
C
      DO 1500 NZ=1,NP(NY,NX)
      RADSL(NZ,NY,NX)=0.0
      RADSW(NZ,NY,NX)=0.0
      RADPL(NZ,NY,NX)=0.0
      RADPW(NZ,NY,NX)=0.0
      RAYSL(NZ,NY,NX)=0.0
      RAYSW(NZ,NY,NX)=0.0
      RAYPL(NZ,NY,NX)=0.0
      RAYPW(NZ,NY,NX)=0.0
      RADS1(NZ,NY,NX)=0.0
      RADW1(NZ,NY,NX)=0.0
      RADP1(NZ,NY,NX)=0.0
      RADQ1(NZ,NY,NX)=0.0
      RAYS1(NZ,NY,NX)=0.0
      RAYW1(NZ,NY,NX)=0.0
      RAYP1(NZ,NY,NX)=0.0
      RAYQ1(NZ,NY,NX)=0.0
      RADS2(NZ,NY,NX)=0.0
      RADW2(NZ,NY,NX)=0.0
      RADP2(NZ,NY,NX)=0.0
      RADQ2(NZ,NY,NX)=0.0
      RAYS2(NZ,NY,NX)=0.0
      RAYW2(NZ,NY,NX)=0.0
      RAYP2(NZ,NY,NX)=0.0
      RAYQ2(NZ,NY,NX)=0.0
C
C     LEAF SURFACE AREA IN EACH INCLINATION CLASS N, AZIMUTH CLASS M,
C     LAYER L AND SPECIES NZ
C
C     TSURFY=unself-shaded leaf area 
C     TSURFZ=unself-shaded leaf area m-2 in each azimuth class
C     TSURFS=TSURFY with shading from canopy layers above
C     TSURFX=TSURFS m-2 
C     TSURWY=unself-shaded stalk area
C     TSURWZ=unself-shaded stalk area m-2 in each azimuth class
C     TSURWS=TSURWY with shading from canopy layers above
C     TSURWX=TSURWS m-2 
C
      DO 1600 N=1,4
      TSURFY=TSURF(N,L,NZ,NY,NX)*CFX(NZ,NY,NX)
      TSURFZ=TSURFY*YAREA
      TSURFS=TSURFY*TAUS(L+1,NY,NX)
      TSURFX=TSURFS*XAREA
      TSURWY=TSURFB(N,L,NZ,NY,NX)*CFW
      TSURWZ=TSURWY*YAREA
      TSURWS=TSURWY*TAUS(L+1,NY,NX)
      TSURWX=TSURWS*XAREA
C
C     ABSORPTION OF DIRECT RADIATION BY SUNLIT LEAF SURFACES
C
C     STOPZ=accumulated horizontal area of intercepted direct radiation
C
      DO 1700 M=1,4
      RADSL(NZ,NY,NX)=RADSL(NZ,NY,NX)+TSURFS*RDNDIR(N,M,NZ,NY,NX)
      RADSW(NZ,NY,NX)=RADSW(NZ,NY,NX)+TSURWS*RDNDIW(N,M,NZ,NY,NX)
      RADPL(NZ,NY,NX)=RADPL(NZ,NY,NX)+TSURFS*PARDIR(N,M,NZ,NY,NX)
      RADPW(NZ,NY,NX)=RADPW(NZ,NY,NX)+TSURWS*PARDIW(N,M,NZ,NY,NX)
      STOPSZ=STOPSZ+(TSURFX+TSURWX)*BETX(N,M)
C
C     BACKSCATTERING OF REFLECTED DIRECT RADIATION
C
      IF(IALBS(N,M).EQ.1)THEN
      RADS1(NZ,NY,NX)=RADS1(NZ,NY,NX)+TSURFS*RDNDIR(N,M,NZ,NY,NX)
      RADW1(NZ,NY,NX)=RADW1(NZ,NY,NX)+TSURWS*RDNDIW(N,M,NZ,NY,NX)
      RADP1(NZ,NY,NX)=RADP1(NZ,NY,NX)+TSURFS*PARDIR(N,M,NZ,NY,NX)
      RADQ1(NZ,NY,NX)=RADQ1(NZ,NY,NX)+TSURWS*PARDIW(N,M,NZ,NY,NX)
C
C     FORWARD SCATTERING OF REFLECTED DIRECT RADIATION
C
      ELSE
      RADS2(NZ,NY,NX)=RADS2(NZ,NY,NX)+TSURFS*RDNDIR(N,M,NZ,NY,NX)
      RADW2(NZ,NY,NX)=RADW2(NZ,NY,NX)+TSURWS*RDNDIW(N,M,NZ,NY,NX)
      RADP2(NZ,NY,NX)=RADP2(NZ,NY,NX)+TSURFS*PARDIR(N,M,NZ,NY,NX)
      RADQ2(NZ,NY,NX)=RADQ2(NZ,NY,NX)+TSURWS*PARDIW(N,M,NZ,NY,NX)
      ENDIF
C
C     INTENSITY OF ABSORBED DIFFUSE RADIATION AT LEAF SURFACES
C
C     RADYN,RADYW,RAPYN,RAPYW=diffuse SW,PAR flux absorbed by leaf,stalk surf
C     OMEGA,OMEGX=incident angle of diffuse radn at leaf,horizontal surface
C
      DO 1750 NN=1,4
      RADYN=RADYL*OMEGA(M,N,NN)*ABSR(NZ,NY,NX)
      RADYW=RADYL*OMEGA(M,N,NN)*ABSRW
      RAPYN=RAPYL*OMEGA(M,N,NN)*ABSP(NZ,NY,NX)
      RAPYW=RAPYL*OMEGA(M,N,NN)*ABSPW
      PARDIF(N,M,L,NZ,NY,NX)=PARDIF(N,M,L,NZ,NY,NX)+RAPYN
      PAR(N,M,L,NZ,NY,NX)=PAR(N,M,L,NZ,NY,NX)+RAPYN
C
C     ABSORPTION OF DIFFUSE RADIATION BY SHADED LEAF SURFACES
C
C     STOPYZ=accumulated horizontal area of intercepted diffuse radiation
C
      RAYSL(NZ,NY,NX)=RAYSL(NZ,NY,NX)+TSURFY*RADYN
      RAYSW(NZ,NY,NX)=RAYSW(NZ,NY,NX)+TSURWY*RADYW
      RAYPL(NZ,NY,NX)=RAYPL(NZ,NY,NX)+TSURFY*RAPYN
      RAYPW(NZ,NY,NX)=RAYPW(NZ,NY,NX)+TSURWY*RAPYW
      STOPYZ=STOPYZ+(TSURFZ+TSURWZ)*OMEGX(M,N,NN)
C
C     BACKSCATTERING OF REFLECTED DIFFUSE RADIATION
C
      IF(IALBY(M,N,NN).EQ.1)THEN
      RAYS1(NZ,NY,NX)=RAYS1(NZ,NY,NX)+TSURFY*RADYN
      RAYW1(NZ,NY,NX)=RAYW1(NZ,NY,NX)+TSURWY*RADYW
      RAYP1(NZ,NY,NX)=RAYP1(NZ,NY,NX)+TSURFY*RAPYN
      RAYQ1(NZ,NY,NX)=RAYQ1(NZ,NY,NX)+TSURWY*RAPYW
C
C     FORWARD SCATTERING OF REFLECTED DIFFUSE RADIATION
C
      ELSE
      RAYS2(NZ,NY,NX)=RAYS2(NZ,NY,NX)+TSURFY*RADYN
      RAYW2(NZ,NY,NX)=RAYW2(NZ,NY,NX)+TSURWY*RADYW
      RAYP2(NZ,NY,NX)=RAYP2(NZ,NY,NX)+TSURFY*RAPYN
      RAYQ2(NZ,NY,NX)=RAYQ2(NZ,NY,NX)+TSURWY*RAPYW
      ENDIF
1750  CONTINUE
1700  CONTINUE
1600  CONTINUE
1500  CONTINUE
C
C     ACCUMULATED INTERCEPTION BY CANOPY LAYER
C
C     XTAUS=interception of direct radiation in current layer
C     STOPZ=accumulated interception of direct radiation from topmost layer
C     TAUS=transmission of direct radiation to next lower layer   
C     
      IF(STOPS+STOPSZ.GT.1.0)THEN
      IF(STOPSZ.GT.ZERO)THEN
      XTAUS=(1.0-STOPS)/((1.0-STOPS)-(1.0-STOPS-STOPSZ))
      ELSE
      XTAUS=0.0
      ENDIF
      TAUS(L+1,NY,NX)=TAUS(L+1,NY,NX)*XTAUS
      STOPSZ=STOPSZ*XTAUS
      DO 1510 NZ=1,NP(NY,NX)
      RADSL(NZ,NY,NX)=RADSL(NZ,NY,NX)*XTAUS
      RADSW(NZ,NY,NX)=RADSW(NZ,NY,NX)*XTAUS
      RADPL(NZ,NY,NX)=RADPL(NZ,NY,NX)*XTAUS
      RADPW(NZ,NY,NX)=RADPW(NZ,NY,NX)*XTAUS
      RADS1(NZ,NY,NX)=RADS1(NZ,NY,NX)*XTAUS
      RADW1(NZ,NY,NX)=RADW1(NZ,NY,NX)*XTAUS
      RADP1(NZ,NY,NX)=RADP1(NZ,NY,NX)*XTAUS
      RADQ1(NZ,NY,NX)=RADQ1(NZ,NY,NX)*XTAUS
      RADS2(NZ,NY,NX)=RADS2(NZ,NY,NX)*XTAUS
      RADW2(NZ,NY,NX)=RADW2(NZ,NY,NX)*XTAUS
      RADP2(NZ,NY,NX)=RADP2(NZ,NY,NX)*XTAUS
      RADQ2(NZ,NY,NX)=RADQ2(NZ,NY,NX)*XTAUS
1510  CONTINUE
      ENDIF
C
C     XTAUY=interception of diffuse radiation in current layer
C     STOPYZ=accumulated interception of diffuse radiation from topmost layer
C     TAUY=transmission of diffuse radiation to next lower layer   
C
      IF(STOPY+STOPYZ.GT.1.0)THEN
      XTAUY=(1.0-STOPY)/((1.0-STOPY)-(1.0-STOPY-STOPYZ))
      TAUY(L+1)=TAUY(L+1)*XTAUY
      STOPYZ=STOPYZ*XTAUY
      DO 1520 NZ=1,NP(NY,NX)
      RAYSL(NZ,NY,NX)=RAYSL(NZ,NY,NX)*XTAUY
      RAYSW(NZ,NY,NX)=RAYSW(NZ,NY,NX)*XTAUY
      RAYPL(NZ,NY,NX)=RAYPL(NZ,NY,NX)*XTAUY
      RAYPW(NZ,NY,NX)=RAYPW(NZ,NY,NX)*XTAUY
      RAYS1(NZ,NY,NX)=RAYS1(NZ,NY,NX)*XTAUY
      RAYW1(NZ,NY,NX)=RAYW1(NZ,NY,NX)*XTAUY
      RAYP1(NZ,NY,NX)=RAYP1(NZ,NY,NX)*XTAUY
      RAYQ1(NZ,NY,NX)=RAYQ1(NZ,NY,NX)*XTAUY
      RAYS2(NZ,NY,NX)=RAYS2(NZ,NY,NX)*XTAUY
      RAYW2(NZ,NY,NX)=RAYW2(NZ,NY,NX)*XTAUY
      RAYP2(NZ,NY,NX)=RAYP2(NZ,NY,NX)*XTAUY
      RAYQ2(NZ,NY,NX)=RAYQ2(NZ,NY,NX)*XTAUY
      DO 1730 N=1,4
      DO 1730 M=1,4
      PARDIF(N,M,L,NZ,NY,NX)=PARDIF(N,M,L,NZ,NY,NX)*XTAUY
      PAR(N,M,L,NZ,NY,NX)=PARDIR(N,M,NZ,NY,NX)+PARDIF(N,M,L,NZ,NY,NX)
1730  CONTINUE
1520  CONTINUE
      ENDIF
C
C     TOTAL RADIATION ABSORBED, REFLECTED AND TRANSMITTED BY ALL PFTs
C
C     RADST,RADWT,RADPT,RADQT=total atmospheric SW,PAR absorbed by leaf,stalk 
C     RA1ST,RA1WT,RA1PT,RA1QT=total backscattered SW,PAR absd by leaf,stalk 
C     RA2ST,RA2WT,RA2PT,RA2QT=total fwd scattered SW,PAR absd by leaf,stalk
C     RAFSL,RAFPL=total fwd scattered SW,PAR to next layer 
C     RABSL,RABPL=total back scattered SW,PAR to next layer
C     RADC,TRADC,RADP,TRADP=total atmospheric SW,PAR absbd by each,all PFT 
C     STOPS,STOPY=accumulated interception of direct,diffuse radiation
C     TAUS,TAUY=transmission of direct,diffuse radiation to next lower layer  
C
      DO 1530 NZ=1,NP(NY,NX)
      RADST=RADSL(NZ,NY,NX)+RAYSL(NZ,NY,NX)
      RADWT=RADSW(NZ,NY,NX)+RAYSW(NZ,NY,NX)
      RADPT=RADPL(NZ,NY,NX)+RAYPL(NZ,NY,NX)
      RADQT=RADPW(NZ,NY,NX)+RAYPW(NZ,NY,NX)
      RA1ST=RADS1(NZ,NY,NX)+RAYS1(NZ,NY,NX)
      RA1WT=RADW1(NZ,NY,NX)+RAYW1(NZ,NY,NX)
      RA1PT=RADP1(NZ,NY,NX)+RAYP1(NZ,NY,NX)
      RA1QT=RADP1(NZ,NY,NX)+RAYQ1(NZ,NY,NX)
      RA2ST=RADS2(NZ,NY,NX)+RAYS2(NZ,NY,NX)
      RA2WT=RADW2(NZ,NY,NX)+RAYW2(NZ,NY,NX)
      RA2PT=RADP2(NZ,NY,NX)+RAYP2(NZ,NY,NX)
      RA2QT=RADQ2(NZ,NY,NX)+RAYQ2(NZ,NY,NX)
      RAFSL(L)=RAFSL(L)+(RADST*TAUR(NZ,NY,NX)
     2+RA2ST*ALBR(NZ,NY,NX)+RA2WT*ALBRW)*YAREA
      RAFPL(L)=RAFPL(L)+(RADPT*TAUP(NZ,NY,NX)
     2+RA2PT*ALBP(NZ,NY,NX)+RA2QT*ALBPW)*YAREA
      RABSL(L)=RABSL(L)+(RA1ST*ALBR(NZ,NY,NX)+RA1WT*ALBRW)*YAREA
      RABPL(L)=RABPL(L)+(RA1PT*ALBP(NZ,NY,NX)+RA1QT*ALBPW)*YAREA
      RADC(NZ,NY,NX)=RADC(NZ,NY,NX)+RADST+RADWT 
      RADP(NZ,NY,NX)=RADP(NZ,NY,NX)+RADPT+RADQT 
      TRADC(NY,NX)=TRADC(NY,NX)+RADST+RADWT 
      TRAPC(NY,NX)=TRAPC(NY,NX)+RADPT+RADQT 
1530  CONTINUE
      STOPS=STOPS+STOPSZ
      STOPY=STOPY+STOPYZ
      TAUS(L,NY,NX)=1.0-STOPS
      TAU0(L,NY,NX)=1.0-TAUS(L,NY,NX)
      TAUY(L)=1.0-STOPY
      ELSE
      RAFSL(L)=RAFSL(L+1)
      RAFPL(L)=RAFPL(L+1)
      TAUS(L,NY,NX)=TAUS(L+1,NY,NX)
      TAU0(L,NY,NX)=1.0-TAUS(L,NY,NX)
      TAUY(L)=TAUY(L+1)
      ENDIF
1800  CONTINUE
C
C     DIRECT AND DIFFUSE RADIATION ABSORBED AT GROUND SURFACE
C
C     RADSG,RADYG,RAPSG,RAPYG=direct,diffuse SW,PAR at horizl ground surface
C     RADS,RAPS =solar beam direct SW,PAR flux
C     TAUS,TAUY=transmission of direct,diffuse radiation below canopy   
C     RADYL,RAPYL=solar beam diffuse SW,PAR flux
C     RASG,RAPG=total SW,PAR at ground surface
C     BETAG,OMEGAG=incident solar,sky angle at ground surface
C
      RADSG=RADS(NY,NX)*TAUS(1,NY,NX)
      RADYG=RADYL*TAUY(1)+RAFSL(1)
      RAPSG=RAPS(NY,NX)*TAUS(1,NY,NX)
      RAPYG=RAPYL*TAUY(1)+RAFPL(1)
      RASG=ABS(BETAG)*RADSG
      RAPG=ABS(BETAG)*RAPSG
      DO 20 N=1,4
      RASG=RASG+ABS(OMEGAG(N,NY,NX))*RADYG
      RAPG=RAPG+ABS(OMEGAG(N,NY,NX))*RAPYG
20    CONTINUE
      RADG(NY,NX)=RASG*AREA(3,NU(NY,NX),NY,NX)
C
C     RADIATION REFLECTED FROM GROUND SURFACE
C
C     VHCPW,VHCPWX=current,minimum snowpack heat capacity
C     ALBW,VOLSS,VOLWS,VOLIS=snowpack albedo,snow,water,ice volume
C     ALBG,ALBS,FSNOW=ground,soil albedo,snow cover fraction
C     THETW1=soil surface water content
C     RABSL,RADPL=SW,PAR backscatter from ground surface
C     TRADG,TRAPG=SW,PAR absorbed by ground surface 
C
      IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
      ALBW=(0.80*VOLSS(NY,NX)+0.30*VOLIS(NY,NX)+0.06*VOLWS(NY,NX))
     2/(VOLSS(NY,NX)+VOLIS(NY,NX)+VOLWS(NY,NX))
      FSNOW=AMIN1((DPTHS(NY,NX)/0.07)**2,1.0)
      ALBG=FSNOW*ALBW+(1.0-FSNOW)*ALBS(NY,NX)
      ELSE
      IF(VOLX(NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      THETW1=AMIN1(POROS(NU(NY,NX),NY,NX)
     2,VOLW(NU(NY,NX),NY,NX)/VOLY(NU(NY,NX),NY,NX))
      ELSE
      THETW1=0.0
      ENDIF
      ALBG=AMIN1(ALBX(NY,NX),ALBS(NY,NX)
     2+AMAX1(0.0,ALBX(NY,NX)-THETW1))
      ENDIF
      RABSL(0)=RASG*ALBG*0.25
      RABPL(0)=RAPG*ALBG*0.25
      TRADG(NY,NX)=(1.0-ALBG)*RASG*AREA(3,NU(NY,NX),NY,NX)
      TRAPG(NY,NX)=(1.0-ALBG)*RAPG*AREA(3,NU(NY,NX),NY,NX)
C
C     ADD RADIATION FROM SCATTERING THROUGH CANOPY LAYERS
C
C     RABSL,RABPL=total backscattered SW,PAR to next layer
C     RAFSL,RAFPL=total fwd scattered SW,PAR to next layer
C     RADYN,RADYW,RAPYN,RAPYW=leaf,stalk SW,PAR absbd fwd+back flux
C     RAYSL,RAYSW,RAYPL,RAYPW=total leaf,stalk SW,PAR absbd fwd+back
C     RADC,TRADC,RADP,TRADP=total SW,PAR absbd by each,all PFT 
C
      RADYL=0.0
      RAPYL=0.0
      TAUY(0)=1.0
      RAFSL(0)=0.0
      RAFPL(0)=0.0
      DO 2800 L=1,JC
      IF(ZL(L-1,NY,NX).GE.DPTHS(NY,NX)-ZERO
     2.AND.ZL(L-1,NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
      RADYL=RADYL*TAUY(L-1)+RAFSL(L-1)+RABSL(L-1)
      RAPYL=RAPYL*TAUY(L-1)+RAFPL(L-1)+RABPL(L-1)
      RAFSL(L)=0.0
      RAFPL(L)=0.0
      DO 2500 NZ=1,NP(NY,NX)
      RAYSL(NZ,NY,NX)=0.0
      RAYSW(NZ,NY,NX)=0.0
      RAYPL(NZ,NY,NX)=0.0
      RAYPW(NZ,NY,NX)=0.0
      DO 2600 N=1,4
      TSURFY=TSURF(N,L,NZ,NY,NX)*CFX(NZ,NY,NX)
      TSURWY=TSURFB(N,L,NZ,NY,NX)*CFW
      DO 2700 M=1,4
      DO 2750 NN=1,4
      RADYN=RADYL*OMEGA(M,N,NN)*ABSR(NZ,NY,NX)
      RADYW=RADYL*OMEGA(M,N,NN)*ABSRW
      RAPYN=RAPYL*OMEGA(M,N,NN)*ABSP(NZ,NY,NX)
      RAPYW=RAPYL*OMEGA(M,N,NN)*ABSPW
      PARDIF(N,M,L,NZ,NY,NX)=PARDIF(N,M,L,NZ,NY,NX)+RAPYN
      PAR(N,M,L,NZ,NY,NX)=PAR(N,M,L,NZ,NY,NX)+RAPYN
      RAYSL(NZ,NY,NX)=RAYSL(NZ,NY,NX)+TSURFY*RADYN
      RAYSW(NZ,NY,NX)=RAYSW(NZ,NY,NX)+TSURWY*RADYW
      RAYPL(NZ,NY,NX)=RAYPL(NZ,NY,NX)+TSURFY*RAPYN
      RAYPW(NZ,NY,NX)=RAYPW(NZ,NY,NX)+TSURWY*RAPYW
2750  CONTINUE
2700  CONTINUE
2600  CONTINUE
      RAFSL(L)=RAFSL(L)+RAYSL(NZ,NY,NX)*TAUR(NZ,NY,NX)*YAREA
      RAFPL(L)=RAFPL(L)+RAYPL(NZ,NY,NX)*TAUP(NZ,NY,NX)*YAREA
      RADC(NZ,NY,NX)=RADC(NZ,NY,NX)+RAYSL(NZ,NY,NX)+RAYSW(NZ,NY,NX)
      RADP(NZ,NY,NX)=RADP(NZ,NY,NX)+RAYPL(NZ,NY,NX)+RAYPW(NZ,NY,NX)
      TRADC(NY,NX)=TRADC(NY,NX)+RAYSL(NZ,NY,NX)+RAYSW(NZ,NY,NX)
      TRAPC(NY,NX)=TRAPC(NY,NX)+RAYPL(NZ,NY,NX)+RAYPW(NZ,NY,NX)
2500  CONTINUE
      ELSE
      RAFSL(L)=RAFSL(L-1)
      RAFPL(L)=RAFPL(L-1)
      RABSL(L)=RABSL(L-1)
      RABPL(L)=RABPL(L-1)
      ENDIF
2800  CONTINUE
C
C     RADIATION AT GROUND SURFACE IF NO CANOPY
C
      ELSE
      RASG=ABS(BETAG)*RADS(NY,NX)
      DO 120 N=1,4
      RASG=RASG+ABS(OMEGAG(N,NY,NX))*RADY(NY,NX)
120   CONTINUE
      RADG(NY,NX)=RASG*AREA(3,NU(NY,NX),NY,NX)
      DO 135 NZ=1,NP(NY,NX)
      RADC(NZ,NY,NX)=0.0
      RADP(NZ,NY,NX)=0.0
135   CONTINUE
      ENDIF
C
C     IF NO RADIATION 
C
      ELSE
      RADG(NY,NX)=0.0
      DO 125 NZ=1,NP(NY,NX)
      RADC(NZ,NY,NX)=0.0
      RADP(NZ,NY,NX)=0.0
125   CONTINUE
      ENDIF
C
C     CANOPY AND GROUND SKY FRACTIONS USED FOR BOUNDARY LAYER CALCULNS
C
C     FRADG=fraction of radiation received by ground surface
C     FRADP=fraction of radiation received by each PFT canopy
C     ARLSS,ARLFS=leaf+stalk area of all PFTs,each PFT
C
      FRADG(NY,NX)=1.0
      IF(ARLSS(NY,NX).GT.ZEROS(NY,NX))THEN
      FRADPT=1.0-EXP(-0.65*ARLSS(NY,NX)/AREA(3,NU(NY,NX),NY,NX))
      DO 145 NZ=1,NP(NY,NX)
      FRADP(NZ,NY,NX)=FRADPT*ARLFS(NZ,NY,NX)/ARLSS(NY,NX)
      FRADG(NY,NX)=FRADG(NY,NX)-FRADP(NZ,NY,NX)
145   CONTINUE
      ELSE
      FRADG(NY,NX)=1.0
      DO 146 NZ=1,NP(NY,NX)
      FRADP(NZ,NY,NX)=0.0
146   CONTINUE
      ENDIF
C     IF(NX.EQ.4.AND.NY.EQ.5)THEN
C     DO 140 NZ=1,NP(NY,NX)
C     WRITE(*,1926)'CANOPY',IYRC,I,J,NX,NY,NZ,FRADP(NZ,NY,NX)
C    2,RADP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX),RAP(NY,NX)
C    2,FRADG(NY,NX),ARLFS(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
C    3,ARLFP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
C    4,ARSTP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
C    4,ARLSS(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
C    5,SSIN(NY,NX),DPTHS(NY,NX),FRADPT 
1926  FORMAT(A10,6I6,30E12.4)
140   CONTINUE
C     ENDIF
C
C     DIVISION OF CANOPY INTO LAYERS WITH EQUAL LAI 
C
C     ZT,ZC=heights of combined canopy,PFT canopy
C     ZL=height to bottom of each canopy layer
C     ARLFC,ARSTC=leaf,stalk area of combined canopy
C     ARLFT,ARSTT=leaf,stalk area of combined canopy layer
C
      ZT(NY,NX)=0.0
      DO 9685 NZ=1,NP(NY,NX)
      ZT(NY,NX)=AMAX1(ZT(NY,NX),ZC(NZ,NY,NX))
9685  CONTINUE
      ZL(JC,NY,NX)=ZT(NY,NX)+0.01
      ZL1(JC,NY,NX)=ZL(JC,NY,NX)
      ZL1(0,NY,NX)=0.0
      ART=(ARLFC(NY,NX)+ARSTC(NY,NX))/JC
      IF(ART.GT.ZEROS(NY,NX))THEN
      DO 2765 L=JC,2,-1
      ARL=ARLFT(L,NY,NX)+ARSTT(L,NY,NX)
      IF(ARL.GT.1.01*ART)THEN
      DZL=ZL(L,NY,NX)-ZL(L-1,NY,NX)
      ZL1(L-1,NY,NX)=ZL(L-1,NY,NX)+0.5*AMIN1(1.0,(ARL-ART)/ARL)*DZL
      ELSEIF(ARL.LT.0.99*ART)THEN
      ARX=ARLFT(L-1,NY,NX)+ARSTT(L-1,NY,NX)
      DZL=ZL(L-1,NY,NX)-ZL(L-2,NY,NX)
      IF(ARX.GT.ZEROS(NY,NX))THEN
      ZL1(L-1,NY,NX)=ZL(L-1,NY,NX)-0.5*AMIN1(1.0,(ART-ARL)/ARX)*DZL
      ENDIF
      ELSE
      ZL1(L-1,NY,NX)=ZL(L-1,NY,NX)
      ENDIF
C     IF(J.EQ.12)THEN
C     WRITE(*,3233)'ZL',I,J,NX,NY,L,ZL1(L,NY,NX),ZL1(L-1,NY,NX)
C    3,ZL(L,NY,NX),ZL(L-1,NY,NX),ART,ARL,ARX
C    2,DZL,ARLFC(NY,NX),ARSTC(NY,NX),ARLFT(L,NY,NX),ARSTT(L,NY,NX)
C    3,ARLFT(L-1,NY,NX),ARSTT(L-1,NY,NX)
C    3,ZL(JC,NY,NX),ZT(NY,NX),ZC(1,NY,NX)
3233  FORMAT(A8,5I4,30E12.4)
C     ENDIF
2765  CONTINUE
      DO 2770 L=JC,2,-1
      ZL(L-1,NY,NX)=ZL1(L-1,NY,NX)
C     ZL(L-1,NY,NX)=AMAX1(0.0,AMIN1(ZL(L,NY,NX)-1.0E-06
C    2,ZL(L-1,NY,NX)))
2770  CONTINUE
      ENDIF
C
C     CANOPY ZERO PLANE AND ROUGHNESS HEIGHTS
C
C     ARLFC,ARSTC=leaf,stalk area of combined canopy
C     DPTHS,DPTH0=snowpack,surface water depths
C     ZT,ZD,ZR=canopy,zero plane displacement,roughness height
C     ZZ=reference height for wind speed
C
      ARLSC=ARLFC(NY,NX)+ARSTC(NY,NX)
      IF(ARLSC.GT.ZEROS(NY,NX)
     2.AND.ZT(NY,NX).GE.DPTHS(NY,NX)-ZERO
     3.AND.ZT(NY,NX).GE.DPTH0(NY,NX)-ZERO)THEN
      ARLSG=ARLSC/AREA(3,NU(NY,NX),NY,NX)
      ZX=EXP(-0.5*ARLSG)
      ZY=1.0-ZX
      ZD(NY,NX)=ZT(NY,NX)*AMAX1(0.0,1.0-2.0/ARLSG*ZY)
      ZE=ZT(NY,NX)*AMAX1(0.05,ZX*ZY)
      ELSE
      ZD(NY,NX)=0.0
      ZE=0.0
      ENDIF
      IF(IFLGW.EQ.1)THEN
      ZZ=Z0(NY,NX)+ZT(NY,NX)
      ELSE
      ZZ=AMAX1(Z0(NY,NX),ZD(NY,NX)+2.0)
      ENDIF
      IF(IETYP(NY,NX).GE.0)THEN
      IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
      ZR(NY,NX)=AMAX1(0.001,ZE,ZW)
      ELSE
      ZR(NY,NX)=AMAX1(0.001,ZE,ZS(NY,NX))
      ENDIF
C
C     CANOPY ISOTHERMAL BOUNDARY LAYER RESISTANCE 
C
C     RAB,RAM=biome canopy,minimum isothermal boundary layer resistance
C     UA=wind speed
C     RIB=canopy isothermal Richardson number
C
      RAB(NY,NX)=AMAX1(RAM,
     2(LOG((ZZ-ZD(NY,NX))/ZR(NY,NX)))**2/(0.168*UA(NY,NX)))
      RIB(NY,NX)=1.27E+08*(ZZ-ZR(NY,NX))/(UA(NY,NX)**2*TKA(NY,NX))
      ELSE
      RAB(NY,NX)=RAM
      RIB(NY,NX)=0.0
      ENDIF
C
C     RESET HOURLY INDICATORS
C
      THRMCX(NY,NX)=THRMC(NY,NX)
      THRMGX(NY,NX)=THRMG(NY,NX)
      CNETX(NY,NX)=TCNET(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      THRMC(NY,NX)=0.0
      THRMG(NY,NX)=0.0
      TLEX(NY,NX)=TLEC(NY,NX)
      TSHX(NY,NX)=TSHC(NY,NX)
      TLEC(NY,NX)=0.0
      TSHC(NY,NX)=0.0
      TRN(NY,NX)=0.0
      TLE(NY,NX)=0.0
      TSH(NY,NX)=0.0
      TGH(NY,NX)=0.0
      TCCAN(NY,NX)=0.0
      TCNET(NY,NX)=0.0
      RECO(NY,NX)=0.0
C
C     CANOPY RETENTION OF PRECIPITATION
C
C     XVOLWC=foliar surface water retention capacity 
C     ARLFP,ARSTP=leaf,stalk area of PFT
C     FLWC,TFLWC=water retention of PFT,combined canopy
C     PRECA=precipitation+irrigation
C     FRADP=fraction of radiation received by each PFT canopy
C     VOLWC=canopy surface water retention
C
      DO 1930 NZ=1,NP(NY,NX)
      VOLWCX=XVOLWC(IGTYP(NZ,NY,NX))
     2*(ARLFP(NZ,NY,NX)+ARSTP(NZ,NY,NX))
      FLWC(NZ,NY,NX)=AMAX1(0.0,AMIN1(PRECA(NY,NX)*FRADP(NZ,NY,NX)
     2,VOLWCX-VOLWC(NZ,NY,NX)))
      TFLWCI(NY,NX)=TFLWCI(NY,NX)+PRECA(NY,NX)*FRADP(NZ,NY,NX)
      TFLWC(NY,NX)=TFLWC(NY,NX)+FLWC(NZ,NY,NX)
C
C     NUMBERS OF TOP AND BOTTOM ROOTED SOIL LAYERS 
C
C     NG=number of uppermost rooted layer
C     NINR=number of lowest rooted layer
C
      NG(NZ,NY,NX)=MAX(NG(NZ,NY,NX),NU(NY,NX))
      NIX(NZ,NY,NX)=MAX(NIX(NZ,NY,NX),NU(NY,NX))
      DO 9790 NR=1,10
      NINR(NR,NZ,NY,NX)=MAX(NINR(NR,NZ,NY,NX),NU(NY,NX))
9790  CONTINUE
1930  CONTINUE
C
C     WRITE SW AND PAR ALBEDO
C
C     IF(ABS(J-ZNOON(NY,NX)).LT.1)THEN
C     IF(RAD(NY,NX).GT.0.0.AND.RAP(NY,NX).GT.0.0)THEN
C     WRITE(19,1927)'ALBEDO',IYRC,I,J,NX,NY
C    2,RAD(NY,NX),RAP(NY,NX),TRADC(NY,NX),TRAPC(NY,NX)
C    3,TRADG(NY,NX),TRAPG(NY,NX)
C    4,(RAD(NY,NX)-TRADC(NY,NX)-TRADG(NY,NX))/RAD(NY,NX) 
C    5,(RAP(NY,NX)-TRAPC(NY,NX)-TRAPG(NY,NX))/RAP(NY,NX) 
1927  FORMAT(A10,5I6,30E12.4)
C     ENDIF
C     ENDIF
9990  CONTINUE
9995  CONTINUE
C
C     FERTILIZER APPLICATIONS OCCUR AT SOLAR NOON
C
      DO 8990 NX=NHW,NHE
      DO 8995 NY=NVN,NVS
      IF(J.EQ.INT(ZNOON(NY,NX)))THEN
C
C     NH4,NH3,UREA,NO3 FERTILIZER APPLICATION
C
C     *A,*B=broadcast,banded
C     Z4,Z3,ZU,ZO=NH4,NH3,urea,NO3
C
      Z4A=FERT(1,I,NY,NX)
      Z3A=FERT(2,I,NY,NX)
      ZUA=FERT(3,I,NY,NX)
      ZOA=FERT(4,I,NY,NX)
      Z4B=FERT(5,I,NY,NX)
      Z3B=FERT(6,I,NY,NX)
      ZUB=FERT(7,I,NY,NX)
      ZOB=FERT(8,I,NY,NX)
C
C     MONOCALCIUM PHOSPHATE OR HYDROXYAPATITE 
C
C     PM*,PH*=Ca(H2PO4)2,apatite
C
      PMA=FERT(9,I,NY,NX)
      PMB=FERT(10,I,NY,NX)
      PHA=FERT(11,I,NY,NX)
C
C     LIME AND GYPSUM
C
C     CAC,CAS=CaCO3,CaSO4
C
      CAC=FERT(12,I,NY,NX)
      CAS=FERT(13,I,NY,NX)
C
C     PLANT(1) AND ANIMAL(2) RESIDUE C, N AND P
C
      OFC(1)=FERT(14,I,NY,NX)
      OFN(1)=FERT(15,I,NY,NX)
      OFP(1)=FERT(16,I,NY,NX)
      OFC(2)=FERT(17,I,NY,NX)
      OFN(2)=FERT(18,I,NY,NX)
      OFP(2)=FERT(19,I,NY,NX)
C’
C     SOIL LAYER NUMBER AT DEPTH OF FERTILIZER APPLICATION
C
C     LFDPTH=layer number
C     CVRDF=fraction of fertilizer applied to surface litter
C
      IF(Z4A+Z3A+ZUA+ZOA+Z4B+Z3B+ZUB+ZOB
     2+PMA+PMB+PHA+CAC+CAS.GT.0.0)THEN
      FDPTHF=FDPTH(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      IF(FDPTHF.LE.0.0.AND.Z4B+Z3B+ZUB+ZOB+PMB.EQ.0.0)THEN
      LFDPTH=0
      CVRDF=1.0-EXP(-0.8E-02*(ORGC(0,NY,NX)/AREA(3,0,NY,NX)))
      ELSE
      DO 65 L=NUI(NY,NX),JZ
      IF(CDPTH(L,NY,NX).GE.FDPTHF)THEN
      LFDPTH=L
      CVRDF=1.0
      GO TO 55
      ENDIF
65    CONTINUE
55    CONTINUE
      ENDIF
      BAREF=1.0-CVRDF
C
C     RESET WIDTH AND DEPTH OF NH4 FERTILIZER BAND IF NEW BAND 
C     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
C
C     ROWN=width of NH4 band row
C     DPNHB,WDNHB=depth,width of NH4 band
C     VLNHB,VLNH4=soil volume in NH4 band,non-band
C
      IF((Z4B+Z3B+ZUB.GT.0.0).OR.((ZNH4B(LFDPTH,NY,NX).GT.0.0
     2.OR.ZNH3B(LFDPTH,NY,NX).GT.0.0).AND.IFNHB(NY,NX).EQ.0))THEN
      IFNHB(NY,NX)=1
      ROWN(NY,NX)=ROWI(I,NY,NX)
      DO 50 L=NUI(NY,NX),JZ
      IF(L.LT.LFDPTH)THEN
      DPNHB(L,NY,NX)=DLYR(3,L,NY,NX)
      WDNHB(L,NY,NX)=0.0
      ELSEIF(L.EQ.LFDPTH)THEN
      DPNHB(L,NY,NX)=AMAX1(0.025,FDPTHF-CDPTH(L-1,NY,NX))
      WDNHB(L,NY,NX)=AMIN1(0.025,ROWN(NY,NX))
      ELSE
      DPNHB(L,NY,NX)=0.0
      WDNHB(L,NY,NX)=0.0
      ENDIF
      IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
      VLNHB(L,NY,NX)=AMIN1(0.999,WDNHB(L,NY,NX)/ROWN(NY,NX)
     2*DPNHB(L,NY,NX)/DLYR(3,L,NY,NX))
      ELSE
      VLNHB(L,NY,NX)=0.0
      ENDIF
      VLNH4(L,NY,NX)=1.0-VLNHB(L,NY,NX)
      ZNH4T=ZNH4S(L,NY,NX)+ZNH4B(L,NY,NX)
      ZNH3T=ZNH3S(L,NY,NX)+ZNH3B(L,NY,NX)
      XN4T=XN4(L,NY,NX)+XNB(L,NY,NX)
      ZNH4S(L,NY,NX)=ZNH4T*VLNH4(L,NY,NX)
      ZNH3S(L,NY,NX)=ZNH3T*VLNH4(L,NY,NX)
      ZNH4B(L,NY,NX)=ZNH4T*VLNHB(L,NY,NX)
      ZNH3B(L,NY,NX)=ZNH3T*VLNHB(L,NY,NX)
      XN4(L,NY,NX)=XN4T*VLNH4(L,NY,NX)
      XNB(L,NY,NX)=XN4T*VLNHB(L,NY,NX)
50    CONTINUE
      DPNH4(NY,NX)=DPNHB(LFDPTH,NY,NX)+CDPTH(LFDPTH-1,NY,NX)
      ENDIF
C
C     RESET WIDTH AND DEPTH OF NO3 FERTILIZER BAND IF NEW BAND 
C     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
C
C     ROWO=width of NO3 band row
C     DPNOB,WDNOB=depth,width of NO3 band
C     VLNOB,VLNO3=soil volume in NO3 band,non-band
C
      IF((Z4B+Z3B+ZUB+ZOB.GT.0.0).OR.((ZNO3B(LFDPTH,NY,NX).GT.0.0
     2.OR.ZNO2B(LFDPTH,NY,NX).GT.0.0).AND.IFNOB(NY,NX).EQ.0))THEN
      IFNOB(NY,NX)=1
      ROWO(NY,NX)=ROWI(I,NY,NX)
      DO 45 L=NUI(NY,NX),JZ
      IF(L.LT.LFDPTH)THEN
      DPNOB(L,NY,NX)=DLYR(3,L,NY,NX)
      WDNOB(L,NY,NX)=0.0
      ELSEIF(L.EQ.LFDPTH)THEN
      DPNOB(L,NY,NX)=AMAX1(0.01,FDPTHF-CDPTH(L-1,NY,NX))
      WDNOB(L,NY,NX)=AMIN1(0.01,ROWO(NY,NX))
      ELSE
      DPNOB(L,NY,NX)=0.0
      WDNOB(L,NY,NX)=0.0
      ENDIF
      IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
      VLNOB(L,NY,NX)=AMIN1(0.999,WDNOB(L,NY,NX)/ROWO(NY,NX)
     2*DPNOB(L,NY,NX)/DLYR(3,L,NY,NX))
      ELSE
      VLNOB(L,NY,NX)=0.0
      ENDIF
      VLNO3(L,NY,NX)=1.0-VLNOB(L,NY,NX)
      ZNO3T=ZNO3S(L,NY,NX)+ZNO3B(L,NY,NX)
      ZNO2T=ZNO2S(L,NY,NX)+ZNO2B(L,NY,NX)
      ZNO3S(L,NY,NX)=ZNO3T*VLNO3(L,NY,NX)
      ZNO2S(L,NY,NX)=ZNO2T*VLNO3(L,NY,NX)
      ZNO3B(L,NY,NX)=ZNO3T*VLNOB(L,NY,NX)
      ZNO2B(L,NY,NX)=ZNO2T*VLNOB(L,NY,NX)
45    CONTINUE
      DPNO3(NY,NX)=DPNOB(LFDPTH,NY,NX)+CDPTH(LFDPTH-1,NY,NX)
      ENDIF
C
C     RESET WIDTH AND DEPTH OF PO4 FERTILIZER BAND IF NEW BAND 
C     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
C
C     ROWP=width of H2PO4 band row
C     DPPOB,WDPOB=depth,width of H2PO4 band
C     VLPOB,VLPO4=soil volume in H2PO4 band,non-band
C
      IF((PMB.GT.0.0).OR.(H2POB(LFDPTH,NY,NX).GT.0.0
     2.AND.IFPOB(NY,NX).EQ.0))THEN
      IFPOB(NY,NX)=1
      ROWP(NY,NX)=ROWI(I,NY,NX)
      DO 40 L=NUI(NY,NX),JZ
      IF(L.LT.LFDPTH)THEN
      DPPOB(L,NY,NX)=DLYR(3,L,NY,NX)
      WDPOB(L,NY,NX)=AMIN1(0.01,ROWP(NY,NX))
      ELSEIF(L.EQ.LFDPTH)THEN
      DPPOB(L,NY,NX)=AMAX1(0.01,FDPTHF-CDPTH(L-1,NY,NX))
      WDPOB(L,NY,NX)=AMIN1(0.01,ROWP(NY,NX))
      ELSE
      DPPOB(L,NY,NX)=0.0
      WDPOB(L,NY,NX)=0.0
      ENDIF
      IF(DLYR(3,L,NY,NX).GT.ZERO2)THEN
      VLPOB(L,NY,NX)=AMIN1(0.999,WDPOB(L,NY,NX)/ROWP(NY,NX)
     2*DPPOB(L,NY,NX)/DLYR(3,L,NY,NX))
      ELSE
      VLPOB(L,NY,NX)=0.0
      ENDIF
      VLPO4(L,NY,NX)=1.0-VLPOB(L,NY,NX)
      H0PO4T=H0PO4(L,NY,NX)+H0POB(L,NY,NX)
      H1PO4T=H1PO4(L,NY,NX)+H1POB(L,NY,NX)
      H2PO4T=H2PO4(L,NY,NX)+H2POB(L,NY,NX)
      H3PO4T=H3PO4(L,NY,NX)+H3POB(L,NY,NX)
      ZFE1PT=ZFE1P(L,NY,NX)+ZFE1PB(L,NY,NX)
      ZFE2PT=ZFE2P(L,NY,NX)+ZFE2PB(L,NY,NX)
      ZCA0PT=ZCA0P(L,NY,NX)+ZCA0PB(L,NY,NX)
      ZCA1PT=ZCA1P(L,NY,NX)+ZCA1PB(L,NY,NX)
      ZCA2PT=ZCA2P(L,NY,NX)+ZCA2PB(L,NY,NX)
      ZMG1PT=ZMG1P(L,NY,NX)+ZMG1PB(L,NY,NX)
      XOH0T=XOH0(L,NY,NX)+XOH0B(L,NY,NX)
      XOH1T=XOH1(L,NY,NX)+XOH1B(L,NY,NX)
      XOH2T=XOH2(L,NY,NX)+XOH2B(L,NY,NX)
      XH1PT=XH1P(L,NY,NX)+XH1PB(L,NY,NX)
      XH2PT=XH2P(L,NY,NX)+XH2PB(L,NY,NX)
      PALPOT=PALPO(L,NY,NX)+PALPB(L,NY,NX)
      PFEPOT=PFEPO(L,NY,NX)+PFEPB(L,NY,NX)
      PCAPDT=PCAPD(L,NY,NX)+PCPDB(L,NY,NX)
      PCAPHT=PCAPH(L,NY,NX)+PCPHB(L,NY,NX)
      PCAPMT=PCAPM(L,NY,NX)+PCPMB(L,NY,NX)
      H0PO4(L,NY,NX)=H0PO4T*VLPO4(L,NY,NX)
      H1PO4(L,NY,NX)=H1PO4T*VLPO4(L,NY,NX)
      H2PO4(L,NY,NX)=H2PO4T*VLPO4(L,NY,NX)
      H3PO4(L,NY,NX)=H3PO4T*VLPO4(L,NY,NX)
      ZFE1P(L,NY,NX)=ZFE1PT*VLPO4(L,NY,NX)
      ZFE2P(L,NY,NX)=ZFE2PT*VLPO4(L,NY,NX)
      ZCA0P(L,NY,NX)=ZCA0PT*VLPO4(L,NY,NX)
      ZCA1P(L,NY,NX)=ZCA1PT*VLPO4(L,NY,NX)
      ZCA2P(L,NY,NX)=ZCA2PT*VLPO4(L,NY,NX)
      ZMG1P(L,NY,NX)=ZMG1PT*VLPO4(L,NY,NX)
      H0POB(L,NY,NX)=H0PO4T*VLPOB(L,NY,NX)
      H1POB(L,NY,NX)=H1PO4T*VLPOB(L,NY,NX)
      H2POB(L,NY,NX)=H2PO4T*VLPOB(L,NY,NX)
      H3POB(L,NY,NX)=H3PO4T*VLPOB(L,NY,NX)
      ZFE1PB(L,NY,NX)=ZFE1PT*VLPOB(L,NY,NX)
      ZFE2PB(L,NY,NX)=ZFE2PT*VLPOB(L,NY,NX)
      ZCA0PB(L,NY,NX)=ZCA0PT*VLPOB(L,NY,NX)
      ZCA1PB(L,NY,NX)=ZCA1PT*VLPOB(L,NY,NX)
      ZCA2PB(L,NY,NX)=ZCA2PT*VLPOB(L,NY,NX)
      ZMG1PB(L,NY,NX)=ZMG1PT*VLPOB(L,NY,NX)
      XOH0(L,NY,NX)=XOH0T*VLPO4(L,NY,NX)
      XOH1(L,NY,NX)=XOH1T*VLPO4(L,NY,NX)
      XOH2(L,NY,NX)=XOH2T*VLPO4(L,NY,NX)
      XH1P(L,NY,NX)=XH1PT*VLPO4(L,NY,NX)
      XH2P(L,NY,NX)=XH2PT*VLPO4(L,NY,NX)
      XOH0B(L,NY,NX)=XOH0T*VLPOB(L,NY,NX)
      XOH1B(L,NY,NX)=XOH1T*VLPOB(L,NY,NX)
      XOH2B(L,NY,NX)=XOH2T*VLPOB(L,NY,NX)
      XH1PB(L,NY,NX)=XH1PT*VLPOB(L,NY,NX)
      XH2PB(L,NY,NX)=XH2PT*VLPOB(L,NY,NX)
      PALPO(L,NY,NX)=PALPOT*VLPO4(L,NY,NX)
      PFEPO(L,NY,NX)=PFEPOT*VLPO4(L,NY,NX)
      PCAPD(L,NY,NX)=PCAPDT*VLPO4(L,NY,NX)
      PCAPH(L,NY,NX)=PCAPHT*VLPO4(L,NY,NX)
      PCAPM(L,NY,NX)=PCAPMT*VLPO4(L,NY,NX)
      PALPB(L,NY,NX)=PALPOT*VLPOB(L,NY,NX)
      PFEPB(L,NY,NX)=PFEPOT*VLPOB(L,NY,NX)
      PCPDB(L,NY,NX)=PCAPDT*VLPOB(L,NY,NX)
      PCPHB(L,NY,NX)=PCAPHT*VLPOB(L,NY,NX)
      PCPMB(L,NY,NX)=PCAPMT*VLPOB(L,NY,NX)
40    CONTINUE
      DPPO4(NY,NX)=DPPOB(LFDPTH,NY,NX)+CDPTH(LFDPTH-1,NY,NX)
      ENDIF
C
C     UPDATE STATE VARIABLES FOR BROADCAST AND BANDED FERTILIZER
C     NH4, NH3, UREA, NO3, PO4, LIME AND GYPSUM IN SOIL
C     AND CONVERT FROM G TO MOLE
C
C     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=bdcast NH4,NH3,urea,NO3 fertilizer
C     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
C     PCAPM1,PCAPD1,PCAPH1=concn of precip CaH2PO4,CaHPO4,apatite in non-band
C     PCAPMB,PCAPDB,PCAPHB=concn of precip CaH2PO4,CaHPO4,apatite in band
C     PCACO,PCASO=precipitated CaCO3,CaSO4
C
      Z4AX=Z4A*AREA(3,LFDPTH,NY,NX)/14.0
      Z3AX=Z3A*AREA(3,LFDPTH,NY,NX)/14.0
      ZUAX=ZUA*AREA(3,LFDPTH,NY,NX)/14.0
      ZOAX=ZOA*AREA(3,LFDPTH,NY,NX)/14.0
      Z4BX=Z4B*AREA(3,LFDPTH,NY,NX)/14.0
      Z3BX=Z3B*AREA(3,LFDPTH,NY,NX)/14.0
      ZUBX=ZUB*AREA(3,LFDPTH,NY,NX)/14.0
      ZOBX=ZOB*AREA(3,LFDPTH,NY,NX)/14.0
      PMAX=PMA*AREA(3,LFDPTH,NY,NX)/62.0
      PMBX=PMB*AREA(3,LFDPTH,NY,NX)/62.0
      PHAX=PHA*AREA(3,LFDPTH,NY,NX)/93.0
      CACX=CAC*AREA(3,LFDPTH,NY,NX)/40.0
      CASX=CAS*AREA(3,LFDPTH,NY,NX)/40.0
      ZNH4FA(LFDPTH,NY,NX)=ZNH4FA(LFDPTH,NY,NX)+Z4AX*CVRDF
      ZNHUFA(LFDPTH,NY,NX)=ZNHUFA(LFDPTH,NY,NX)+ZUAX*CVRDF
      ZNO3FA(LFDPTH,NY,NX)=ZNO3FA(LFDPTH,NY,NX)+ZOAX*CVRDF
      ZNH4FB(LFDPTH,NY,NX)=ZNH4FB(LFDPTH,NY,NX)+Z4BX*CVRDF
      ZNHUFB(LFDPTH,NY,NX)=ZNHUFB(LFDPTH,NY,NX)+ZUBX*CVRDF
      ZNO3FB(LFDPTH,NY,NX)=ZNO3FB(LFDPTH,NY,NX)+ZOBX*CVRDF
      PCAPM(LFDPTH,NY,NX)=PCAPM(LFDPTH,NY,NX)
     2+PMAX*VLPO4(LFDPTH,NY,NX)*CVRDF
      PCPMB(LFDPTH,NY,NX)=PCPMB(LFDPTH,NY,NX)
     2+PMAX*VLPOB(LFDPTH,NY,NX)*CVRDF+PMBX*CVRDF
      PCAPH(LFDPTH,NY,NX)=PCAPH(LFDPTH,NY,NX)
     2+PHAX*VLPO4(LFDPTH,NY,NX)*CVRDF
      PCPHB(LFDPTH,NY,NX)=PCPHB(LFDPTH,NY,NX)
     2+PHAX*VLPOB(LFDPTH,NY,NX)*CVRDF
      IF(LFDPTH.EQ.0)THEN
      ZNH4FA(NU(NY,NX),NY,NX)=ZNH4FA(NU(NY,NX),NY,NX)+Z4AX*BAREF
      ZNH3FA(NU(NY,NX),NY,NX)=ZNH3FA(NU(NY,NX),NY,NX)+Z3AX 
      ZNHUFA(NU(NY,NX),NY,NX)=ZNHUFA(NU(NY,NX),NY,NX)+ZUAX*BAREF
      ZNO3FA(NU(NY,NX),NY,NX)=ZNO3FA(NU(NY,NX),NY,NX)+ZOAX*BAREF
      ZNH4FB(NU(NY,NX),NY,NX)=ZNH4FB(NU(NY,NX),NY,NX)+Z4BX*BAREF
      ZNH3FB(NU(NY,NX),NY,NX)=ZNH3FB(NU(NY,NX),NY,NX)+Z3BX 
      ZNHUFB(NU(NY,NX),NY,NX)=ZNHUFB(NU(NY,NX),NY,NX)+ZUBX*BAREF
      ZNO3FB(NU(NY,NX),NY,NX)=ZNO3FB(NU(NY,NX),NY,NX)+ZOBX*BAREF
      PCAPM(NU(NY,NX),NY,NX)=PCAPM(NU(NY,NX),NY,NX)
     2+PMAX*VLPO4(NU(NY,NX),NY,NX)*BAREF
      PCPMB(NU(NY,NX),NY,NX)=PCPMB(NU(NY,NX),NY,NX)
     2+PMAX*VLPOB(NU(NY,NX),NY,NX)*BAREF+PMBX*BAREF
      PCAPH(NU(NY,NX),NY,NX)=PCAPH(NU(NY,NX),NY,NX)
     2+PHAX*VLPO4(NU(NY,NX),NY,NX)*BAREF
      PCPHB(NU(NY,NX),NY,NX)=PCPHB(NU(NY,NX),NY,NX)
     2+PHAX*VLPOB(NU(NY,NX),NY,NX)*BAREF
      ELSE
      ZNH3FA(LFDPTH,NY,NX)=ZNH3FA(LFDPTH,NY,NX)+Z3AX*CVRDF
      ZNH3FB(LFDPTH,NY,NX)=ZNH3FB(LFDPTH,NY,NX)+Z3BX*CVRDF
      ENDIF
      PCACO(NU(NY,NX),NY,NX)=PCACO(NU(NY,NX),NY,NX)+CACX
      PCASO(NU(NY,NX),NY,NX)=PCASO(NU(NY,NX),NY,NX)+CASX
      TZIN=TZIN+14.0*(Z4AX+Z3AX+ZUAX+ZOAX+Z4BX+Z3BX+ZUBX+ZOBX)
      TPIN=TPIN+62.0*(PMAX+PMBX)+93.0*PHAX
      TIONIN=TIONIN+2.0*(CACX+CASX)
      UFERTN(NY,NX)=UFERTN(NY,NX)+14.0*(Z4AX+Z4BX+Z3AX+Z3BX
     2+ZUAX+ZUBX+ZOAX+ZOBX)
      UFERTP(NY,NX)=UFERTP(NY,NX)+62.0*(PMAX+PMBX)+93.0*PHAX
      ENDIF
C
C     SOIL LAYER NUMBER IN WHICH PLANT OR ANIMAL RESIDUES ARE APPLIED
C
C     LFDPTH=layer number
C
      IF(OFC(1)+OFC(2).GT.0.0)THEN
      DO 2985 L=0,JZ
      FDPTHM=FDPTH(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      IF(FDPTHM.LE.0.0)THEN
      LFDPTH=0
      GO TO 2980
      ELSEIF(CDPTH(L,NY,NX).GE.FDPTHM)THEN
      LFDPTH=L
      GO TO 2980
      ENDIF
2985  CONTINUE
2980  CONTINUE
C
C     ALLOCATION OF PLANT RESIDUE APPLICATION TO
C     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
C
C     CFOSC=fraction of litter allocated to protein(1)
C     soluble CH2O(2), cellulose(3) and lignin(4) 
C     ITYPE=litter type entered in fertilizer input file
C
C     MAIZE
C
      IF(IYTYP(1,I,NY,NX).EQ.1)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.080
      CFOSC(2,1,LFDPTH,NY,NX)=0.245
      CFOSC(3,1,LFDPTH,NY,NX)=0.613
      CFOSC(4,1,LFDPTH,NY,NX)=0.062
C
C     WHEAT
C
      ELSEIF(IYTYP(1,I,NY,NX).EQ.2)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.125
      CFOSC(2,1,LFDPTH,NY,NX)=0.171
      CFOSC(3,1,LFDPTH,NY,NX)=0.560
      CFOSC(4,1,LFDPTH,NY,NX)=0.144
C
C     SOYBEAN
C
      ELSEIF(IYTYP(1,I,NY,NX).EQ.3)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.138
      CFOSC(2,1,LFDPTH,NY,NX)=0.426
      CFOSC(3,1,LFDPTH,NY,NX)=0.316
      CFOSC(4,1,LFDPTH,NY,NX)=0.120
C
C     OLD STRAW
C
      ELSEIF(IYTYP(1,I,NY,NX).EQ.4)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.075
      CFOSC(2,1,LFDPTH,NY,NX)=0.125
      CFOSC(3,1,LFDPTH,NY,NX)=0.550
      CFOSC(4,1,LFDPTH,NY,NX)=0.250
C
C     STRAW
C
      ELSEIF(IYTYP(1,I,NY,NX).EQ.5)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.036
      CFOSC(2,1,LFDPTH,NY,NX)=0.044
      CFOSC(3,1,LFDPTH,NY,NX)=0.767
      CFOSC(4,1,LFDPTH,NY,NX)=0.153
C
C     COMPOST
C
      ELSEIF(IYTYP(1,I,NY,NX).EQ.6)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.143
      CFOSC(2,1,LFDPTH,NY,NX)=0.015
      CFOSC(3,1,LFDPTH,NY,NX)=0.640
      CFOSC(4,1,LFDPTH,NY,NX)=0.202
C
C     GREEN MANURE
C
      ELSEIF(IYTYP(1,I,NY,NX).EQ.7)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.202
      CFOSC(2,1,LFDPTH,NY,NX)=0.013
      CFOSC(3,1,LFDPTH,NY,NX)=0.560
      CFOSC(4,1,LFDPTH,NY,NX)=0.225
C
C     SIMPLE SUBSTRATE
C
      ELSEIF(IYTYP(1,I,NY,NX).EQ.10)THEN
      CFOSC(1,1,LFDPTH,NY,NX)=0.000
      CFOSC(2,1,LFDPTH,NY,NX)=1.000
      CFOSC(3,1,LFDPTH,NY,NX)=0.000
      CFOSC(4,1,LFDPTH,NY,NX)=0.000
      ELSE
      CFOSC(1,1,LFDPTH,NY,NX)=0.075
      CFOSC(2,1,LFDPTH,NY,NX)=0.125
      CFOSC(3,1,LFDPTH,NY,NX)=0.550
      CFOSC(4,1,LFDPTH,NY,NX)=0.250
      ENDIF
C
C     ALLOCATION OF ANIMAL MANURE APPLICATION TO
C     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
C
C     RUMINANT
C
      IF(IYTYP(2,I,NY,NX).EQ.1)THEN
      CFOSC(1,2,LFDPTH,NY,NX)=0.036
      CFOSC(2,2,LFDPTH,NY,NX)=0.044
      CFOSC(3,2,LFDPTH,NY,NX)=0.630
      CFOSC(4,2,LFDPTH,NY,NX)=0.290
C
C     NON-RUMINANT
C
      ELSEIF(IYTYP(2,I,NY,NX).EQ.2)THEN
      CFOSC(1,2,LFDPTH,NY,NX)=0.138
      CFOSC(2,2,LFDPTH,NY,NX)=0.401
      CFOSC(3,2,LFDPTH,NY,NX)=0.316
      CFOSC(4,2,LFDPTH,NY,NX)=0.145
C
C     GRAZING
C
      ELSEIF(IYTYP(2,I,NY,NX).EQ.3)THEN
      CFOSC(1,2,LFDPTH,NY,NX)=0.036
      CFOSC(2,2,LFDPTH,NY,NX)=0.044
      CFOSC(3,2,LFDPTH,NY,NX)=0.630
      CFOSC(4,2,LFDPTH,NY,NX)=0.290
C
C     OTHER
C
      ELSE
      CFOSC(1,2,LFDPTH,NY,NX)=0.138
      CFOSC(2,2,LFDPTH,NY,NX)=0.401
      CFOSC(3,2,LFDPTH,NY,NX)=0.316
      CFOSC(4,2,LFDPTH,NY,NX)=0.145
      ENDIF
C
C     DISTRIBUTE RESIDUE APPLICATION AMONG COMPONENTS OF RESIDUE COMPLEX
C
C     OFC,OFN,OFP=litter C,N,P application from fertilizer file
C
      DO 2965 K=1,2
      OSCI=OFC(K)*AREA(3,LFDPTH,NY,NX)
      OSNI=OFN(K)*AREA(3,LFDPTH,NY,NX)
      OSPI=OFP(K)*AREA(3,LFDPTH,NY,NX)
      IF(BKVL(LFDPTH,NY,NX).GT.ZEROS(NY,NX))THEN
      CORGCX=OSCI/BKVL(LFDPTH,NY,NX)
      ELSE
      CORGCX=0.55E+06
      ENDIF
      OSCX=0.0
      OSNX=0.0
      OSPX=0.0
C
C     BIOMASSES OF MICROBIAL POPULATIONS IN RESIDUE
C
C     OMC,OMN,OMP=microbial biomass in litter application
C     OMCI=microbial biomass content in litter 
C     OMCF,OMCA=hetero,autotrophic biomass composition in litter 
C
      DO 2960 N=1,7
      DO 2961 M=1,3
      OMC1=AMAX1(0.0,AMIN1(OSCI*OMCI(M,K)*OMCF(N),OSCI-OSCX))
      OMN1=AMAX1(0.0,AMIN1(OMC1*CNOMC(M,N,K),OSNI-OSNX))
      OMP1=AMAX1(0.0,AMIN1(OMC1*CPOMC(M,N,K),OSPI-OSPX))
      OMC(M,N,K,LFDPTH,NY,NX)=OMC(M,N,K,LFDPTH,NY,NX)+OMC1
      OMN(M,N,K,LFDPTH,NY,NX)=OMN(M,N,K,LFDPTH,NY,NX)+OMN1
      OMP(M,N,K,LFDPTH,NY,NX)=OMP(M,N,K,LFDPTH,NY,NX)+OMP1
C     WRITE(*,2345)'OMCI',I,J,LFDPTH,K,N,M
C    2,OMC1,OMN1,OMP1,OSCI,OMCI(M,K)
C    2,OMCF(N),OSCX,CNOMC(M,N,K),CPOMC(M,N,K),OSNI,OSPI
C    2,OMC(M,N,K,LFDPTH,NY,NX),OMN(M,N,K,LFDPTH,NY,NX)
2345  FORMAT(A8,6I4,20E12.4)
      OSCX=OSCX+OMC1
      OSNX=OSNX+OMN1
      OSPX=OSPX+OMP1
      DO 2962 NN=1,7
      OMC(M,NN,5,LFDPTH,NY,NX)=OMC(M,NN,5,LFDPTH,NY,NX)+OMC1*OMCA(NN)
      OMN(M,NN,5,LFDPTH,NY,NX)=OMN(M,NN,5,LFDPTH,NY,NX)+OMN1*OMCA(NN)
      OMP(M,NN,5,LFDPTH,NY,NX)=OMP(M,NN,5,LFDPTH,NY,NX)+OMP1*OMCA(NN)
C     WRITE(*,2346)'OMCA',I,J,LFDPTH,K,N,NN,M
C    2,OMC1,OMCA(NN),OMC(M,NN,5,LFDPTH,NY,NX)
2346  FORMAT(A8,7I4,20E12.4)
      OSCX=OSCX+OMC1*OMCA(NN)
      OSNX=OSNX+OMN1*OMCA(NN)
      OSPX=OSPX+OMP1*OMCA(NN)
2962  CONTINUE
2961  CONTINUE
2960  CONTINUE
C
C     DOC, DON AND DOP IN RESIDUE
C
C     OQC,OQN,OQP=DOC,DON,DOP in litter
C
      OQC1=AMIN1(0.1*OSCX,OSCI-OSCX)
      OQN1=AMIN1(0.1*OSNX,OSNI-OSNX)
      OQP1=AMIN1(0.1*OSPX,OSPI-OSPX)
      OQC(K,LFDPTH,NY,NX)=OQC(K,LFDPTH,NY,NX)+OQC1
      OQN(K,LFDPTH,NY,NX)=OQN(K,LFDPTH,NY,NX)+OQN1
      OQP(K,LFDPTH,NY,NX)=OQP(K,LFDPTH,NY,NX)+OQP1
C
C     REMAINDER DISTRIBUTED TO RESIDUE FRACTIONS
C
C     OSC,OSN,OSP,OSA=SOC,SON,SOP,colonized SOC in litter
C     VOLT=litter volume
C     UORGF,UFERTN,UFERTP=accumulated litter C,N,P application
C     TNBP=accumulated net biome productivity
C
      OSCX=OSCX+OQC1
      OSNX=OSNX+OQN1
      OSPX=OSPX+OQP1
      CNOFT=0.0
      CPOFT=0.0
      IF(OSCI-OSCX.GT.ZEROS(NY,NX))THEN
      RNT=0.0
      RPT=0.0
      DO 965 M=1,4
      RNT=RNT+(OSCI-OSCX)*CFOSC(M,K,LFDPTH,NY,NX)*CNOFC(M,K)
      RPT=RPT+(OSCI-OSCX)*CFOSC(M,K,LFDPTH,NY,NX)*CPOFC(M,K)
965   CONTINUE
      FRNT=(OSNI-OSNX)/RNT
      FRPT=(OSPI-OSPX)/RPT
      DO 970 M=1,4
      CNOF(M)=CNOFC(M,K)*FRNT 
      CPOF(M)=CPOFC(M,K)*FRPT 
      CNOFT=CNOFT+CFOSC(M,K,LFDPTH,NY,NX)*CNOF(M)
      CPOFT=CPOFT+CFOSC(M,K,LFDPTH,NY,NX)*CPOF(M)
970   CONTINUE
      ELSE
      DO 975 M=1,4
      CNOF(M)=0.0
      CPOF(M)=0.0
975   CONTINUE
      ENDIF
      DO 2970 M=1,4
      OSC1=CFOSC(M,K,LFDPTH,NY,NX)*(OSCI-OSCX)
      IF(CNOFT.GT.ZERO)THEN
      OSN1=CFOSC(M,K,LFDPTH,NY,NX)*CNOF(M)/CNOFT*(OSNI-OSNX)
      ELSE
      OSN1=0.0
      ENDIF
      IF(CPOFT.GT.ZERO)THEN
      OSP1=CFOSC(M,K,LFDPTH,NY,NX)*CPOF(M)/CPOFT*(OSPI-OSPX)
      ELSE
      OSP1=0.0
      ENDIF
      OSC(M,K,LFDPTH,NY,NX)=OSC(M,K,LFDPTH,NY,NX)+OSC1
      OSA(M,K,LFDPTH,NY,NX)=OSA(M,K,LFDPTH,NY,NX)+OSC1*OMCI(1,K)
      OSN(M,K,LFDPTH,NY,NX)=OSN(M,K,LFDPTH,NY,NX)+OSN1
      OSP(M,K,LFDPTH,NY,NX)=OSP(M,K,LFDPTH,NY,NX)+OSP1
      IF(LFDPTH.EQ.0)THEN
      VOLT(LFDPTH,NY,NX)=VOLT(LFDPTH,NY,NX)+OSC1*1.0E-06/BKRS(1)
      ENDIF
2970  CONTINUE
      TORGF=TORGF+OSCI
      TORGN=TORGN+OSNI
      TORGP=TORGP+OSPI
      UORGF(NY,NX)=UORGF(NY,NX)+OSCI
      UFERTN(NY,NX)=UFERTN(NY,NX)+OSNI
      UFERTP(NY,NX)=UFERTP(NY,NX)+OSPI
      IF(IYTYP(2,I,NY,NX).LT.3)THEN
      TNBP(NY,NX)=TNBP(NY,NX)+OSCI
      ENDIF
2965  CONTINUE
      ENDIF
C
C     FERTILIZER UREA, NITRIFICATION INHIBITORS
C
C     IYTYP=fertilizer release type from fertilizer input file
C     FERT=fertilizer type from fertilizer input file
C     IUTYP=urea hydrolysis inhibitor type (1=no,2=yes)
C     ZNHU0,ZNHUI=initial,current urea hydrolysis inhibition activity
C     ZNFN0,ZNFNI=initial,current nitrification inhibition activity 
C
      IF(FERT(3,I,NY,NX).GT.0.0.OR.FERT(7,I,NY,NX).GT.0.0)THEN
      IF(IYTYP(0,I,NY,NX).EQ.0)THEN
      IUTYP(NY,NX)=0
      ELSEIF(IYTYP(0,I,NY,NX).EQ.1.OR.IYTYP(0,I,NY,NX).EQ.3)THEN
      IUTYP(NY,NX)=1
      ELSE
      IUTYP(NY,NX)=2
      ENDIF 
      DO 9964 L=0,NL(NY,NX)
      IF(L.EQ.LFDPTH)THEN 
      ZNHU0(L,NY,NX)=1.0
      ZNHUI(L,NY,NX)=1.0
      ELSE
      ZNHU0(L,NY,NX)=0.0
      ZNHUI(L,NY,NX)=0.0
      ENDIF
9964  CONTINUE
      ENDIF
      IF(IYTYP(0,I,NY,NX).EQ.3.OR.IYTYP(0,I,NY,NX).EQ.4)THEN
      DO 9965 L=0,NL(NY,NX)
      IF(L.EQ.LFDPTH)THEN 
      ZNFN0(L,NY,NX)=1.0
      ZNFNI(L,NY,NX)=1.0
      ELSE
      ZNFN0(L,NY,NX)=0.0
      ZNFNI(L,NY,NX)=0.0
      ENDIF
9965  CONTINUE
      ENDIF
      ENDIF
8995  CONTINUE
8990  CONTINUE
      RETURN
      END




