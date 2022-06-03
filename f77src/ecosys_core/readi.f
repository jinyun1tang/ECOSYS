
      SUBROUTINE readi(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NTZ
     2,NTZX,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE READS ALL SOIL AND TOPOGRAPHIC INPUT FILES
C
      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"
      include "blk2a.h"
      include "blk2b.h"
      include "blk2c.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk17.h"
      include "blk11a.h"
      DIMENSION NA(10),ND(10),NM(JY,JX),DHI(JX),DVI(JY)
      CHARACTER*16 DATA(30),DATAC(30,250,250),DATAP(JP,JY,JX)
     2,DATAM(JP,JY,JX),DATAX(JP),DATAY(JP),DATAZ(JP,JY,JX)
     3,OUTS(10),OUTP(10),OUTFILS(10,JY,JX),OUTFILP(10,JP,JY,JX)
      CHARACTER*3 CHOICE(102,20)
      CHARACTER*8 CDATE
      CHARACTER*16 OUTW,OUTI,OUTT,OUTN,OUTF
      CHARACTER*4 CHARY
      CHARACTER*1 TTYPE,CTYPE,IVAR(20),VAR(50),TYP(50)
      CHARACTER*80 PREFIX
      DIMENSION IDAT(20),DAT(50),DATK(50)
	  dimension datav(40)
      PARAMETER (TWILGT=0.06976)
C
C     OPEN OUTPUT LOGFILES,AND SITE,TOPOGRAPHY FILES FROM
C     FILE NAMES IN DATA ARRAYS LOADED IN 'MAIN'
C
      OPEN(18,FILE='logfile1',STATUS='UNKNOWN')
      OPEN(19,FILE='logfile2',STATUS='UNKNOWN')
      OPEN(20,FILE='logfile3',STATUS='UNKNOWN')
      OPEN(1,FILE=TRIM(PREFIX)//DATA(1),STATUS='OLD')
      OPEN(7,FILE=TRIM(PREFIX)//DATA(2),STATUS='OLD')
      WRITE(18,5000)' 08 JUL 2021'
5000  FORMAT(A16)
      NF=1
      NFX=1
      NTZ=0
C
C     READ SITE DATA
C
C     ALATG,ALTIG,ATCAG=latitude,altitude,MAT(oC)
C     IDTBLG=water table flag
C        :0=none
C        :1,2=natural stationary,mobile
C        :3,4=artificial stationary,mobile
C     OXYEG,Z2GEG,CO2EIG,CH4EG,Z2OEG,ZNH3EG=atm O2,N2,CO2,CH4,N2O,NH3 (ppm)
C     IETYPG,ISALTG,IERSNG=Koppen climate zone,salt,erosion options
C     NCNG=1:lateral connections between grid cells,3:no connections
C     DTBLIG,DTBLDIG=depth of natural,artificial water table (IDTBLG)
C     DTBLGG=slope of natural water table relative to landscape surface
C     RCHQNG,RCHQEG,RCHQSG,RCHQWG=boundary condns for N,E,S,W surface runoff
C     RCHGNUG,RCHGEUG,RCHGSUG,RCHGWUG=bound condns for N,E,S,W subsurf flow
C     RCHGNTG,RCHGETG,RCHGSTG,RCHGWTG=N,E,S,W distance to water table (m)
C     RCHGDG=lower boundary conditions for water flow
C     DHI=width of each W-E landscape column
C     DVI=width of each N-S landscape row
C  IERSNG=erosion options
C       :=0 means allowing  freeze- thaw to change elevation
C       :=1 means allowing freeze-thaw plus erosion to change elevation
C       :=2 means allowing freeze-thaw plus SOC accumulation to change elevation
C       :=3 means allowing freeze-thaw plus SOC accumulation, plus erosion to change elevation
C       :=-1 means no change in elevation.
      READ(1,*)(datav(jj),jj=1,4)
	  ALATG=datav(1)
	  ALTIG=datav(2)
	  ATCAG=datav(3)
	  IDTBLG=int(datav(4))
C      READ(1,*)ALATG,ALTIG,ATCAG,IDTBLG
      READ(1,*)OXYEG,Z2GEG,CO2EIG,CH4EG,Z2OEG,ZNH3EG
      READ(1,*)IETYPG,ISALTG,IERSNG,NCNG,DTBLIG,DTBLDIG,DTBLGG
      READ(1,*)RCHQNG,RCHQEG,RCHQSG,RCHQWG,RCHGNUG,RCHGEUG,RCHGSUG
     2,RCHGWUG,RCHGNTG,RCHGETG,RCHGSTG,RCHGWTG,RCHGDG
      READ(1,*)(DHI(NX),NX=1,NHE)
      READ(1,*)(DVI(NY),NY=1,NVS)
      CLOSE(1)
      DO 9895 NX=NHW,NHE
      DO 9890 NY=NVN,NVS
      ALAT(NY,NX)=ALATG
      ALTI(NY,NX)=ALTIG
      ATCAI(NY,NX)=ATCAG
      IDTBL(NY,NX)=IDTBLG
      OXYE(NY,NX)=OXYEG
      Z2GE(NY,NX)=Z2GEG
      CO2EI(NY,NX)=CO2EIG
      CH4E(NY,NX)=CH4EG
      Z2OE(NY,NX)=Z2OEG
      ZNH3E(NY,NX)=ZNH3EG
      IETYP(NY,NX)=IETYPG
      NCN(NY,NX)=NCNG
      DTBLI(NY,NX)=DTBLIG
      DTBLDI(NY,NX)=DTBLDIG
      DTBLG(NY,NX)=DTBLGG
      RCHQN(NY,NX)=RCHQNG
      RCHQE(NY,NX)=RCHQEG
      RCHQS(NY,NX)=RCHQSG
      RCHQW(NY,NX)=RCHQWG
      RCHGNU(NY,NX)=RCHGNUG
      RCHGEU(NY,NX)=RCHGEUG
      RCHGSU(NY,NX)=RCHGSUG
      RCHGWU(NY,NX)=RCHGWUG
      RCHGNT(NY,NX)=RCHGNTG
      RCHGET(NY,NX)=RCHGETG
      RCHGST(NY,NX)=RCHGSTG
      RCHGWT(NY,NX)=RCHGWTG
      RCHGD(NY,NX)=RCHGDG
      DH(NY,NX)=DHI(NX)
      DV(NY,NX)=DVI(NY)
      CO2E(NY,NX)=CO2EI(NY,NX)
      H2GE(NY,NX)=1.0E-03
C
C     CALCULATE MAXIMUM DAYLENTH FOR PLANT PHENOLOGY
C
C     DYLM=maximum daylength (h)
C
      IF(ALAT(NY,NX).GT.0.0)THEN
      XI=173
      ELSE
      XI=356
      ENDIF
      DECDAY=XI+100
      DECLIN=SIN((DECDAY*0.9863)*1.7453E-02)*(-23.47)
      AZI=SIN(ALAT(NY,NX)*1.7453E-02)*SIN(DECLIN*1.7453E-02)
      DEC=COS(ALAT(NY,NX)*1.7453E-02)*COS(DECLIN*1.7453E-02)
      IF(AZI/DEC.GE.1.0-TWILGT)THEN
      DYLM(NY,NX)=24.0
      ELSEIF(AZI/DEC.LE.-1.0+TWILGT)THEN
      DYLM(NY,NX)=0.0
      ELSE
      DYLM(NY,NX)=12.0*(1.0+2.0/3.1416*ASIN(TWILGT+AZI/DEC))
      ENDIF
9890  CONTINUE
9895  CONTINUE
C
C     READ TOPOGRAPHY DATA AND SOIL FILE NAME FOR EACH GRID CELL
C
C     for each unit within the landscape:
C     NH1,NV1,NH2,NV2=NW,SE column,row
C     ASPX=N,E,S,W aspect (o)
C     SL1,SL2=EW,NS slope (o)
C     DPTHSX=initial snowpack depth
C
50    READ(7,*,END=20)NH1,NV1,NH2,NV2,ASPX,SL0,SLX,DPTHSX
      READ(7,52)DATA(7)
52    FORMAT(A16)
C
C     OPEN AND READ SOIL FILE
C
      OPEN(9,FILE=TRIM(PREFIX)//DATA(7),STATUS='OLD')
      DO 9995 NX=NH1,NH2
      DO 9990 NY=NV1,NV2
C
C     SURFACE SLOPES AND ASPECTS
C
      ASP(NY,NX)=ASPX
      SL(NY,NX)=SL0
      DPTHS(NY,NX)=DPTHSX
C
C     CONVERT ASPECT TO GEOMETRIC FORMAT
C
      ASP(NY,NX)=450.0-ASP(NY,NX)
      IF(ASP(NY,NX).GE.360.0)ASP(NY,NX)=ASP(NY,NX)-360.0
C
C     SURFACE PROPERTIES
C
C     PSIFC,PSIWP=water potentials at field capacity,wilting point (MPa)
C     ALBS=wet soil albedo
C     PH=litter pH
C     RSC,RSC,RSP=C,N,P in fine(1,0),woody(0,0),manure(2,0) surface litter (g m-2)
C     IXTYP=surface litter type:1=plant,2=manure
C     NUI,NJ=number of soil surface layer,maximum rooting layer
C     NL1,NL2=number of additional layers below NJ with,without data in file
C     ISOILR=natural(0),reconstructed(1) soil profile
C
C      READ(9,*)PSIFC(NY,NX),PSIWP(NY,NX),ALBS(NY,NX),PH(0,NY,NX)
C     2,RSC(1,0,NY,NX),RSN(1,0,NY,NX),RSP(1,0,NY,NX)
C     3,RSC(0,0,NY,NX),RSN(0,0,NY,NX),RSP(0,0,NY,NX)
C     4,RSC(2,0,NY,NX),RSN(2,0,NY,NX),RSP(2,0,NY,NX)
C     5,IXTYP(1,NY,NX),IXTYP(2,NY,NX)
C     6,NUI(NY,NX),NJ(NY,NX),NL1,NL2,ISOILR(NY,NX)
	  READ(9,*)(datav(jj),jj=1,20)
      PSIFC(NY,NX)=datav(1)
	  PSIWP(NY,NX)=datav(2)
	  ALBS(NY,NX) =datav(3)
	  PH(0,NY,NX) =datav(4)
      RSC(1,0,NY,NX) =datav(5)
	  RSN(1,0,NY,NX) =datav(6)
	  RSP(1,0,NY,NX) =datav(7)
      RSC(0,0,NY,NX) =datav(8)
	  RSN(0,0,NY,NX) =datav(9)
	  RSP(0,0,NY,NX) =datav(10)
      RSC(2,0,NY,NX) =datav(11)
	  RSN(2,0,NY,NX) =datav(12)
	  RSP(2,0,NY,NX) =datav(13)
      IXTYP(1,NY,NX) =int(datav(14))
	  IXTYP(2,NY,NX) =int(datav(15))
      NUI(NY,NX) = int(datav(16))
	  NJ(NY,NX)  =int(datav(17))
	  NL1=int(datav(18))
	  NL2=int(datav(19))
	  ISOILR(NY,NX)=int(datav(20))

      NU(NY,NX)=NUI(NY,NX)
      NK(NY,NX)=NJ(NY,NX)+1
      NM(NY,NX)=NJ(NY,NX)+NL1
      NLI(NY,NX)=NM(NY,NX)+NL2
      NL(NY,NX)=NLI(NY,NX)
C
C     PHYSICAL PROPERTIES
C
C     CDPTH=depth to bottom (m)
C     BKDSI=initial bulk density (Mg m-3,0=water)
C
      READ(9,*)(CDPTH(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      DO L=NU(NY,NX),NM(NY,NX)
        if(CDPTH(L,NY,NX)<ZERO)then
          write(*,*)'not sufficient input data in line 2'//
     2'of the soil file'
          stop
        endif
      ENDDO
      READ(9,*)(BKDSI(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     HYDROLOGIC PROPERTIES
C
C     FC,WP=field capacity,wilting point:<0=unknown (m3 m-3)
C     SCNV,SCNH=vertical,lateral Ksat:<0=unknown (mm h-1)
C
      READ(9,*)(FC(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(WP(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(SCNV(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(SCNH(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     PHYSICAL PROPERTIES
C
C     CSAND,CSILT=sand,silt contents (kg Mg-1)
C     FHOL,ROCK=macropore,rock fraction
C
      READ(9,*)(CSAND(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CSILT(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(FHOL(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(ROCK(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     CHEMICAL PROPERTIES
C
C     PH=pH
C     CEC,AEC=cation,anion exchange capacity:CEC<0=unknown (cmol Kg-1)
C
      READ(9,*)(PH(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CEC(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(AEC(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     ORGANIC C, N AND P CONCENTRATIONS
C
C     CORGC,CORGR=total SOC,POC(part of SOC) (kg Mg-1)
C     CORGN,CORGP=SON,SOP:<0=unknown (g Mg-1)
C
      READ(9,*)(CORGC(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CORGR(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CORGN(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CORGP(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     INORGANIC N AND P CONCENTRATIONS
C
C     CNH4,CNO3,CPO4=soluble+exchangeable NH4,NO3,H2PO4 (g Mg-1)
C
      READ(9,*)(CNH4(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CNO3(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CPO4(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     CATION AND ANION CONCENTRATIONS
C
C     C*=soluble concentration from sat. paste extract (g Mg-1)
C     AL,FE,CA,MG,NA,KA,SO4,CL=Al,Fe,Ca,Mg,Na,K,SO4-S,Cl
C
      READ(9,*)(CAL(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CFE(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CCA(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CMG(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CNA(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CKA(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CSO4(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CCL(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     PRECIPITATED MINERAL CONCENTRATIONS
C
C     CALPO,CFEPO,CCAPD,CCAPH=AlPO4,FePO4,CaHPO4,apatite (g Mg-1)
C     CALOH,CFEOH,CCACO,CCASO=AlOH3,FeOH3,CaSO4,CaCO3 (g Mg-1)
C
      READ(9,*)(CALPO(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CFEPO(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CCAPD(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CCAPH(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CALOH(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CFEOH(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CCACO(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(CCASO(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     GAPON SELECTIVITY CO-EFFICIENTS
C
C     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
C     Ca-NH4,Ca-H,Ca-Al,Ca-Mg,Ca-Na,Ca-K
C
      READ(9,*)(GKC4(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(GKCH(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(GKCA(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(GKCM(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(GKCN(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(GKCK(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     INITIAL WATER, ICE CONTENTS
C
      READ(9,*)(THW(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(THI(L,NY,NX),L=NU(NY,NX),NM(NY,NX))
C
C     THW,THI=initial water,ice:>1=satd,1=FC,0=WP,<0=0,0-1=m3 m-3
C
C     INITIAL PLANT AND ANIMAL RESIDUE C, N AND P
C
C     RSC,RSC,RSP=C,N,P in fine(1),woody(0),manure(2) litter (g m-2)
C
      READ(9,*)(RSC(1,L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(RSN(1,L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(RSP(1,L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(RSC(0,L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(RSN(0,L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(RSP(0,L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(RSC(2,L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(RSN(2,L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      READ(9,*)(RSP(2,L,NY,NX),L=NU(NY,NX),NM(NY,NX))
      REWIND(9)
      RSC(1,0,NY,NX)=AMAX1(1.0E-06,RSC(1,0,NY,NX))
      RSN(1,0,NY,NX)=AMAX1(0.04E-06,RSN(1,0,NY,NX))
      RSP(1,0,NY,NX)=AMAX1(0.004E-06,RSP(1,0,NY,NX))
      SCNV(0,NY,NX)=10.0*0.098
C
C     SET FLAGS FOR ESTIMATING FC,WP,SCNV,SCNH IF UNKNOWN
C
C     ISOIL=flag for calculating FC(1),WP(2),SCNV(3),SCNH(4)
C
      DO 25 L=NU(NY,NX),NM(NY,NX)
      IF(FC(L,NY,NX).LT.0.0)THEN
      ISOIL(1,L,NY,NX)=1
      PSIFC(NY,NX)=-0.033
      ELSE
      ISOIL(1,L,NY,NX)=0
      ENDIF
      IF(WP(L,NY,NX).LT.0.0)THEN
      ISOIL(2,L,NY,NX)=1
      PSIWP(NY,NX)=-1.5
      ELSE
      ISOIL(2,L,NY,NX)=0
      ENDIF
      IF(SCNV(L,NY,NX).LT.0.0)THEN
      ISOIL(3,L,NY,NX)=1
      ELSE
      ISOIL(3,L,NY,NX)=0
      ENDIF
      IF(SCNH(L,NY,NX).LT.0.0)THEN
      ISOIL(4,L,NY,NX)=1
      ELSE
      ISOIL(4,L,NY,NX)=0
      ENDIF
25    CONTINUE
C
C     FILL OUT SOIL BOUNDARY LAYERS ABOVE ROOTING ZONE (NOT USED)
C
      IF(NU(NY,NX).GT.1)THEN
      DO 31 L=NU(NY,NX)-1,0,-1
      IF(BKDSI(L+1,NY,NX).GT.0.025)THEN
      CDPTH(L,NY,NX)=CDPTH(L+1,NY,NX)-0.01
      ELSE
      CDPTH(L,NY,NX)=CDPTH(L+1,NY,NX)-0.02
      ENDIF
      IF(L.GT.0)THEN
      BKDSI(L,NY,NX)=BKDSI(L+1,NY,NX)
      FC(L,NY,NX)=FC(L+1,NY,NX)
      WP(L,NY,NX)=WP(L+1,NY,NX)
      SCNV(L,NY,NX)=SCNV(L+1,NY,NX)
      SCNH(L,NY,NX)=SCNH(L+1,NY,NX)
      CSAND(L,NY,NX)=CSAND(L+1,NY,NX)
      CSILT(L,NY,NX)=CSILT(L+1,NY,NX)
      CCLAY(L,NY,NX)=CCLAY(L+1,NY,NX)
      FHOL(L,NY,NX)=FHOL(L+1,NY,NX)
      ROCK(L,NY,NX)=ROCK(L+1,NY,NX)
      PH(L,NY,NX)=PH(L+1,NY,NX)
      CEC(L,NY,NX)=CEC(L+1,NY,NX)
      AEC(L,NY,NX)=AEC(L+1,NY,NX)
      CORGC(L,NY,NX)=1.00*CORGC(L+1,NY,NX)
      CORGR(L,NY,NX)=1.00*CORGR(L+1,NY,NX)
      CORGN(L,NY,NX)=1.00*CORGN(L+1,NY,NX)
      CORGP(L,NY,NX)=1.00*CORGP(L+1,NY,NX)
      CNH4(L,NY,NX)=CNH4(L+1,NY,NX)
      CNO3(L,NY,NX)=CNO3(L+1,NY,NX)
      CPO4(L,NY,NX)=CPO4(L+1,NY,NX)
      CAL(L,NY,NX)=CAL(L+1,NY,NX)
      CFE(L,NY,NX)=CFE(L+1,NY,NX)
      CCA(L,NY,NX)=CCA(L+1,NY,NX)
      CMG(L,NY,NX)=CMG(L+1,NY,NX)
      CNA(L,NY,NX)=CNA(L+1,NY,NX)
      CKA(L,NY,NX)=CKA(L+1,NY,NX)
      CSO4(L,NY,NX)=CSO4(L+1,NY,NX)
      CCL(L,NY,NX)=CCL(L+1,NY,NX)
      CALOH(L,NY,NX)=CALOH(L+1,NY,NX)
      CFEOH(L,NY,NX)=CFEOH(L+1,NY,NX)
      CCACO(L,NY,NX)=CCACO(L+1,NY,NX)
      CCASO(L,NY,NX)=CCASO(L+1,NY,NX)
      CALPO(L,NY,NX)=CALPO(L+1,NY,NX)
      CFEPO(L,NY,NX)=CFEPO(L+1,NY,NX)
      CCAPD(L,NY,NX)=CCAPD(L+1,NY,NX)
      CCAPH(L,NY,NX)=CCAPH(L+1,NY,NX)
      GKC4(L,NY,NX)=GKC4(L+1,NY,NX)
      GKCH(L,NY,NX)=GKCH(L+1,NY,NX)
      GKCA(L,NY,NX)=GKCA(L+1,NY,NX)
      GKCM(L,NY,NX)=GKCM(L+1,NY,NX)
      GKCN(L,NY,NX)=GKCN(L+1,NY,NX)
      GKCK(L,NY,NX)=GKCK(L+1,NY,NX)
      THW(L,NY,NX)=THW(L+1,NY,NX)
      THI(L,NY,NX)=THI(L+1,NY,NX)
      ISOIL(1,L,NY,NX)=ISOIL(1,L+1,NY,NX)
      ISOIL(2,L,NY,NX)=ISOIL(2,L+1,NY,NX)
      ISOIL(3,L,NY,NX)=ISOIL(3,L+1,NY,NX)
      ISOIL(4,L,NY,NX)=ISOIL(4,L+1,NY,NX)
      RSC(1,L,NY,NX)=0.0
      RSN(1,L,NY,NX)=0.0
      RSP(1,L,NY,NX)=0.0
      RSC(0,L,NY,NX)=0.0
      RSN(0,L,NY,NX)=0.0
      RSP(0,L,NY,NX)=0.0
      RSC(2,L,NY,NX)=0.0
      RSN(2,L,NY,NX)=0.0
      RSP(2,L,NY,NX)=0.0
      ENDIF
31    CONTINUE
      ENDIF
C
C     ADD SOIL BOUNDARY LAYERS BELOW SOIL ZONE
C
      DO 32 L=NM(NY,NX)+1,JZ
      CDPTH(L,NY,NX)=2.0*CDPTH(L-1,NY,NX)-1.0*CDPTH(L-2,NY,NX)
      BKDSI(L,NY,NX)=BKDSI(L-1,NY,NX)
      FC(L,NY,NX)=FC(L-1,NY,NX)
      WP(L,NY,NX)=WP(L-1,NY,NX)
      SCNV(L,NY,NX)=SCNV(L-1,NY,NX)
      SCNH(L,NY,NX)=SCNH(L-1,NY,NX)
      CSAND(L,NY,NX)=CSAND(L-1,NY,NX)
      CSILT(L,NY,NX)=CSILT(L-1,NY,NX)
      CCLAY(L,NY,NX)=CCLAY(L-1,NY,NX)
      FHOL(L,NY,NX)=FHOL(L-1,NY,NX)
      ROCK(L,NY,NX)=ROCK(L-1,NY,NX)
      PH(L,NY,NX)=PH(L-1,NY,NX)
      CEC(L,NY,NX)=CEC(L-1,NY,NX)
      AEC(L,NY,NX)=AEC(L-1,NY,NX)
C     IF(IDTBL(NY,NX).EQ.0)THEN
      CORGC(L,NY,NX)=0.25*CORGC(L-1,NY,NX)
      CORGR(L,NY,NX)=0.25*CORGR(L-1,NY,NX)
      CORGN(L,NY,NX)=0.25*CORGN(L-1,NY,NX)
      CORGP(L,NY,NX)=0.25*CORGP(L-1,NY,NX)
C     ELSE
C     CORGC(L,NY,NX)=CORGC(L-1,NY,NX)
C     CORGR(L,NY,NX)=CORGR(L-1,NY,NX)
C     CORGN(L,NY,NX)=CORGN(L-1,NY,NX)
C     CORGP(L,NY,NX)=CORGP(L-1,NY,NX)
C     ENDIF
      CNH4(L,NY,NX)=CNH4(L-1,NY,NX)
      CNO3(L,NY,NX)=CNO3(L-1,NY,NX)
      CPO4(L,NY,NX)=CPO4(L-1,NY,NX)
      CAL(L,NY,NX)=CAL(L-1,NY,NX)
      CFE(L,NY,NX)=CFE(L-1,NY,NX)
      CCA(L,NY,NX)=CCA(L-1,NY,NX)
      CMG(L,NY,NX)=CMG(L-1,NY,NX)
      CNA(L,NY,NX)=CNA(L-1,NY,NX)
      CKA(L,NY,NX)=CKA(L-1,NY,NX)
      CSO4(L,NY,NX)=CSO4(L-1,NY,NX)
      CCL(L,NY,NX)=CCL(L-1,NY,NX)
      CALOH(L,NY,NX)=CALOH(L-1,NY,NX)
      CFEOH(L,NY,NX)=CFEOH(L-1,NY,NX)
      CCACO(L,NY,NX)=CCACO(L-1,NY,NX)
      CCASO(L,NY,NX)=CCASO(L-1,NY,NX)
      CALPO(L,NY,NX)=CALPO(L-1,NY,NX)
      CFEPO(L,NY,NX)=CFEPO(L-1,NY,NX)
      CCAPD(L,NY,NX)=CCAPD(L-1,NY,NX)
      CCAPH(L,NY,NX)=CCAPH(L-1,NY,NX)
      GKC4(L,NY,NX)=GKC4(L-1,NY,NX)
      GKCH(L,NY,NX)=GKCH(L-1,NY,NX)
      GKCA(L,NY,NX)=GKCA(L-1,NY,NX)
      GKCM(L,NY,NX)=GKCM(L-1,NY,NX)
      GKCN(L,NY,NX)=GKCN(L-1,NY,NX)
      GKCK(L,NY,NX)=GKCK(L-1,NY,NX)
      THW(L,NY,NX)=THW(L-1,NY,NX)
      THI(L,NY,NX)=THI(L-1,NY,NX)
      ISOIL(1,L,NY,NX)=ISOIL(1,L-1,NY,NX)
      ISOIL(2,L,NY,NX)=ISOIL(2,L-1,NY,NX)
      ISOIL(3,L,NY,NX)=ISOIL(3,L-1,NY,NX)
      ISOIL(4,L,NY,NX)=ISOIL(4,L-1,NY,NX)
      RSC(1,L,NY,NX)=0.0
      RSN(1,L,NY,NX)=0.0
      RSP(1,L,NY,NX)=0.0
      RSC(0,L,NY,NX)=0.0
      RSN(0,L,NY,NX)=0.0
      RSP(0,L,NY,NX)=0.0
      RSC(2,L,NY,NX)=0.0
      RSN(2,L,NY,NX)=0.0
      RSP(2,L,NY,NX)=0.0
32    CONTINUE
C
C     CALCULATE DERIVED SOIL PROPERTIES FROM INPUT SOIL PROPERTIES
C
C     FMPR=micropore fraction excluding macropore,rock
C     SCNV,SCNH=vertical,lateral Ksat converted to m2 MPa-1 h-1
C     CSAND,CSILT,CCLAY=sand,silt,clay content converted to g Mg-1
C     CORGC,CORGR=SOC,POC converted to g Mg-1
C     CEC,AEC=cation,anion exchange capacity converted to mol Mg-1
C     CNH4...=solute concentrations converted to mol Mg-1
C
      DO 28 L=1,NL(NY,NX)
C     BKDSI(L,NY,NX)=BKDSI(L,NY,NX)/(1.0-FHOL(L,NY,NX))
      BKDS(L,NY,NX)=BKDSI(L,NY,NX)
      IF(BKDS(L,NY,NX).EQ.0.0)FHOL(L,NY,NX)=0.0
      FMPR(L,NY,NX)=(1.0-ROCK(L,NY,NX))*(1.0-FHOL(L,NY,NX))
C     FC(L,NY,NX)=FC(L,NY,NX)/(1.0-FHOL(L,NY,NX))
C     WP(L,NY,NX)=WP(L,NY,NX)/(1.0-FHOL(L,NY,NX))
      SCNV(L,NY,NX)=0.098*SCNV(L,NY,NX)*FMPR(L,NY,NX)
      SCNH(L,NY,NX)=0.098*SCNH(L,NY,NX)*FMPR(L,NY,NX)
      CCLAY(L,NY,NX)=AMAX1(0.0,1.0E+03-(CSAND(L,NY,NX)
     2+CSILT(L,NY,NX)))
      CORGC(L,NY,NX)=CORGC(L,NY,NX)*1.0E+03
      CORGR(L,NY,NX)=CORGR(L,NY,NX)*1.0E+03
      CORGCI(L,NY,NX)=CORGC(L,NY,NX)
      FHOLI(L,NY,NX)=FHOL(L,NY,NX)
      CSAND(L,NY,NX)=CSAND(L,NY,NX)
     2*1.0E-03*AMAX1(0.0,(1.0-CORGC(L,NY,NX)/0.55E+06))
      CSILT(L,NY,NX)=CSILT(L,NY,NX)
     2*1.0E-03*AMAX1(0.0,(1.0-CORGC(L,NY,NX)/0.55E+06))
      CCLAY(L,NY,NX)=CCLAY(L,NY,NX)
     2*1.0E-03*AMAX1(0.0,(1.0-CORGC(L,NY,NX)/0.55E+06))
      CEC(L,NY,NX)=CEC(L,NY,NX)*10.0
      AEC(L,NY,NX)=AEC(L,NY,NX)*10.0
      CNH4(L,NY,NX)=CNH4(L,NY,NX)/14.0
      CNO3(L,NY,NX)=CNO3(L,NY,NX)/14.0
      CPO4(L,NY,NX)=CPO4(L,NY,NX)/31.0
      CAL(L,NY,NX)=CAL(L,NY,NX)/27.0
      CFE(L,NY,NX)=CFE(L,NY,NX)/56.0
      CCA(L,NY,NX)=CCA(L,NY,NX)/40.0
      CMG(L,NY,NX)=CMG(L,NY,NX)/24.3
      CNA(L,NY,NX)=CNA(L,NY,NX)/23.0
      CKA(L,NY,NX)=CKA(L,NY,NX)/39.1
      CSO4(L,NY,NX)=CSO4(L,NY,NX)/32.0
      CCL(L,NY,NX)=CCL(L,NY,NX)/35.5
      CALPO(L,NY,NX)=CALPO(L,NY,NX)/31.0
      CFEPO(L,NY,NX)=CFEPO(L,NY,NX)/31.0
      CCAPD(L,NY,NX)=CCAPD(L,NY,NX)/31.0
      CCAPH(L,NY,NX)=CCAPH(L,NY,NX)/(31.0*3.0)
      CALOH(L,NY,NX)=CALOH(L,NY,NX)/27.0
      CFEOH(L,NY,NX)=CFEOH(L,NY,NX)/56.0
      CCACO(L,NY,NX)=CCACO(L,NY,NX)/40.0
      CCASO(L,NY,NX)=CCASO(L,NY,NX)/40.0
C
C     ESTIMATE SON,SOP,CEC IF UNKNOWN
C     BIOCHEMISTRY 130:117-131
C
      IF(CORGN(L,NY,NX).LT.0.0)THEN
      CORGN(L,NY,NX)=AMIN1(0.125*CORGC(L,NY,NX)
     2,8.9E+02*(CORGC(L,NY,NX)/1.0E+04)**0.80)
C     WRITE(*,1111)'CORGN',L,CORGN(L,NY,NX),CORGC(L,NY,NX)
      ENDIF
      IF(CORGP(L,NY,NX).LT.0.0)THEN
      CORGP(L,NY,NX)=AMIN1(0.0125*CORGC(L,NY,NX)
     2,1.2E+02*(CORGC(L,NY,NX)/1.0E+04)**0.52)
C     WRITE(*,1111)'CORGP',L,CORGP(L,NY,NX),CORGC(L,NY,NX)
      ENDIF
      IF(CEC(L,NY,NX).LT.0.0)THEN
      CEC(L,NY,NX)=10.0*(200.0*2.0*CORGC(L,NY,NX)/1.0E+06
     2+80.0*CCLAY(L,NY,NX)+20.0*CSILT(L,NY,NX)
     3+5.0*CSAND(L,NY,NX))
      ENDIF
28    CONTINUE
      CORGC(0,NY,NX)=0.55E+06
      FMPR(0,NY,NX)=1.0
9990  CONTINUE
9995  CONTINUE
      CLOSE(9)
      GO TO 50
20    CONTINUE
      CLOSE(7)
      DO 9975 NX=NHW,NHE
      NL(NVS+1,NX)=NL(NVS,NX)
C     WRITE(*,2223)'NHE',NX,NHW,NHE,NVS,NL(NVS,NX)
9975  CONTINUE
      DO 9970 NY=NVN,NVS
      NL(NY,NHE+1)=NL(NY,NHE)
C     WRITE(*,2223)'NVS',NY,NVN,NVS,NHE,NL(NY,NHE)
2223  FORMAT(A8,6I4)
9970  CONTINUE
      NL(NVS+1,NHE+1)=NL(NVS,NHE)
      IOLD=0
      RETURN
      END
