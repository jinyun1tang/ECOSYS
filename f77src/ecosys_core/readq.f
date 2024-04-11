
      SUBROUTINE readq(NA,ND,NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NTZ
     2,NTZX,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE READS INPUT DATA FROM PLANT SPECIES
C     AND MANAGEMENT FILES IDENTIFIED IN 'ROUTQ'
C
      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"
      include "blk9a.h"
      include "blk9b.h"
      include "blk9c.h"
      include "blk17.h"
      DIMENSION NA(10),ND(10)
      CHARACTER*16 DATA(30),DATAC(30,250,250),DATAP(JP,JY,JX)
     2,DATAM(JP,JY,JX),DATAX(JP),DATAY(JP),DATAZ(JP,JY,JX)
     3,OUTS(10),OUTP(10),OUTFILS(10,JY,JX),OUTFILP(10,JP,JY,JX)
      CHARACTER*3 CHOICE(102,20)
      CHARACTER*8 CDATE
      CHARACTER*80 PREFIX
      PARAMETER (TWILGT=0.06976)
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      DO 9985 NZ=1,NP(NY,NX)	
C
C     OPEN PLANT SPECIES (11), PLANT MANAGEMENT (12) FILES FROM
C     FILE NAMES IN DATA ARRAYS LOADED IN MAIN.F
C
C     PREFIX=path for files in current or higher level directory
C     DATAP=plant species file name
C
      IF(DATAP(NZ,NY,NX).NE.'NO')THEN
      OPEN(11,FILE=TRIM(PREFIX)//DATAP(NZ,NY,NX),STATUS='OLD')
      ENDIF
      IF(DATAM(NZ,NY,NX).NE.'NO')THEN
      OPEN(12,FILE=TRIM(PREFIX)//DATAM(NZ,NY,NX),STATUS='OLD')
      ENDIF
C
C     READ INPUTS FOR EACH PLANT SPECIES
C
      IF(DATAP(NZ,NY,NX).NE.'NO')THEN
C
C     PLANT FUNCTIONAL TYPE (SPECIES)
C
C     ICTYP=photosynthesis type:3=C3,4=C4
C     IGTYP=root profile
C              :0=shallow (eg bryophytes)
C              :1=intermediate(herbaceous)
C              :2=deep (trees)
C              :3=deeper 
C          =turnover of above-ground biomass
C              :0 or 1:all above-ground (herbaceous),
C              :>1: foliage and twigs (trees) 
C     ISTYP=growth habit:0=annual,1=perennial
C     IDTYP=growth habit:0=determinate,1=indetermimate
C     INTYP=N2 fixation:1,2,3=rapid to slow root symbiosis 
C                             (e.g.legumes),
C                       4,5,6=rapid to slow canopy symbiosis 
C                             (e.g. cyanobacteria)
C     IWTYP=phenology type:0=evergreen
C                         :1=cold deciduous
C                         :2=drought deciduous
C                         :3=cold and drought deciduous
C     IPTYP=photoperiod type:0=day neutral
C                           :1=short day
C                           :2=long day
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     IRTYP=storage organ:0=above ground,1=below ground
C     MY=mycorrhizal:1=no,2=yes
C     ZTYPI=thermal adaptation zone
C              :1=arctic,boreal
C              :2=cool temperate,
C              :3=warm temperate
C              :4=subtropical
C              :5=tropical
C
      READ(11,*)ICTYP(NZ,NY,NX),IGTYP(NZ,NY,NX),ISTYP(NZ,NY,NX)
     2,IDTYP(NZ,NY,NX),INTYP(NZ,NY,NX),IWTYP(NZ,NY,NX)
     3,IPTYP(NZ,NY,NX),IBTYP(NZ,NY,NX),IRTYP(NZ,NY,NX),MY(NZ,NY,NX)
     4,ZTYPI(NZ,NY,NX)
C
C     PHOTOSYNTHETIC PROPERTIES
C
C     VCMX,VOMX=specific rubisco carboxylase,oxygenase activity 
C        at 25 oC (umol C,O g-1 s-1)
C     VCMX4=specific PEP carboxylase activity at 25 oC (umol g-1 s-1)
C     XKCO2,XKO2,XKCO24=Km for VCMX,VOMX,VCMX4 at 25 oC (uM)
C     RUBP,PEPC=fraction of leaf protein in rubisco, PEP carboxylase
C     ETMX=specific chlorophyll activity at 25 oC (umol e- g-1 s-1)
C     CHL=fraction of leaf protein in mesophyll(C3), bundle sheath(C4)
C        chlorophyll
C     CHL4=fraction of leaf protein in mesophyll chlorophyll(C4)
C     FCO2=intercellular:atmospheric CO2 concentration ratio
C          usual values C3=0.70, C4=0.45 
C
      READ(11,*)VCMX(NZ,NY,NX),VOMX(NZ,NY,NX),VCMX4(NZ,NY,NX)
     2,XKCO2(NZ,NY,NX),XKO2(NZ,NY,NX),XKCO24(NZ,NY,NX)
     3,RUBP(NZ,NY,NX),PEPC(NZ,NY,NX),ETMX(NZ,NY,NX),CHL(NZ,NY,NX)
     4,CHL4(NZ,NY,NX),FCO2(NZ,NY,NX)
C
C     OPTICAL PROPERTIES
C
C     ALBP,ALBP=leaf SW,PAR albedo
C     TAUR,TAUP=leaf SW,PAR transmission
C
      READ(11,*)ALBR(NZ,NY,NX),ALBP(NZ,NY,NX),TAUR(NZ,NY,NX)
     2,TAUP(NZ,NY,NX)
C
C     PHENOLOGICAL PROPERTIES
C
C     XRNI,XRLA=rate of node initiation,leaf appearance at 25oC (h-1)
C     CTC=chilling temperature for CO2 fixation, seed loss (oC)
C     VRNLI,VRNXI=hour requirement for spring leafout,autumn leafoff
C     WDLF=leaf length:width ratio
C     PB=nonstructural C concentration needed for branching
C     GROUPX=node number required to start floral initiation
C     XTLI=node number in seed at planting  
C     XDL=critical photoperiod for floral initiation (h) if IPTYP>0
C        :<0=maximum daylength at summer solstice for latitude from site file
C     XPPD=difference between XDL and photoperiod for floral induction (h)
C
      READ(11,*)XRNI(NZ,NY,NX),XRLA(NZ,NY,NX),CTC(NZ,NY,NX)
     2,VRNLI,VRNXI,WDLF(NZ,NY,NX),PB(NZ,NY,NX)
      READ(11,*)GROUPX(NZ,NY,NX),XTLI(NZ,NY,NX),XDL(NZ,NY,NX)
     2,XPPD(NZ,NY,NX)
C
C     MORPHOLOGICAL PROPERTIES
C
C     SLA1,SSL1,SNL1=growth in leaf area,petiole length,internode
C        length vs growth in mass
C     CLASS=fraction of leaf area in 0-22.5,22.5-45,45-67.5,67.5-90 deg. 
C        inclination classes
C     CFI=initial clumping factor used for self-shading in radiation interception
C     ANGBR,ANGSH=stem,petiole angle from horizontal (deg.)
C     STMX=maximum potential seed sites from pre-anthesis 
C        stalk growth
C     SDMX=maximum seed number per STMX during post-anthesis seed set
C     GRMX=maximum seed size per SDMX (g)
C     GRDM=seed size at planting (g)
C     GFILL=grain filling rate at 25 oC (g seed-1 h-1)
C     WTSTDI=mass of dead standing biomass at planting (g m-2) 
C
      READ(11,*)SLA1(NZ,NY,NX),SSL1(NZ,NY,NX),SNL1(NZ,NY,NX)
      READ(11,*)(CLASS(N,NZ,NY,NX),N=1,4),CFI(NZ,NY,NX) 
     2,ANGBR(NZ,NY,NX),ANGSH(NZ,NY,NX)
      READ(11,*)STMX(NZ,NY,NX),SDMX(NZ,NY,NX),GRMX(NZ,NY,NX)
     2,GRDM(NZ,NY,NX),GFILL(NZ,NY,NX),WTSTDI(NZ,NY,NX)
C
C     ROOT CHARACTERISTICS
C
C     RRAD1M,RRAD2M=radius of primary,secondary roots (m)
C     PORT=root porosity
C     PR=nonstructural C concentration needed for root branching
C     RSRR,RSRA=radial,axial root resistivity (MPa h m-3,MPa h m-2)
C     PTSHT=rate constant for equilibrating shoot-root nonstructural 
C        C concentration
C     RTFQ=secondary root branching frequency (m-1)
C
      READ(11,*)RRAD1M(1,NZ,NY,NX),RRAD2M(1,NZ,NY,NX) 
     2,PORTI(1,NZ,NY,NX),PR(NZ,NY,NX),RSRR(1,NZ,NY,NX) 
     3,RSRA(1,NZ,NY,NX),PTSHT(NZ,NY,NX),RTFQ(NZ,NY,NX)
C
C     ROOT UPTAKE PARAMETERS
C
C     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake (g m-2 h-1),
C        Km (uM), min concentration (uM)      
C     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake (g m-2 h-1),
C        Km (uM), min concentration (uM)      
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake (g m-2 h-1),
C        Km (uM), min concentration (uM)      
C
      READ(11,*)UPMXZH(1,NZ,NY,NX),UPKMZH(1,NZ,NY,NX)
     2,UPMNZH(1,NZ,NY,NX)
      READ(11,*)UPMXZO(1,NZ,NY,NX),UPKMZO(1,NZ,NY,NX)
     2,UPMNZO(1,NZ,NY,NX)
      READ(11,*)UPMXPO(1,NZ,NY,NX),UPKMPO(1,NZ,NY,NX)
     2,UPMNPO(1,NZ,NY,NX)
C
C     WATER RELATIONS
C
C     OSMO=leaf osmotic potential at zero leaf water potential (MPa)
C     RCS=shape parameter for stomatal resistance vs leaf turgor
C        potential
C     RSMX=cuticular resistance (s m-1)
C
      READ(11,*)OSMO(NZ,NY,NX),RCS(NZ,NY,NX),RSMX(NZ,NY,NX)
C
C     ORGAN GROWTH YIELDS
C
C     DM*=DM-C production vs nonstructural C consumption (g g-1)
C     *LF=leaf,*SHE=petiole,*STK=stalk,*RSV=stalk reserve,*HSK=husk
C     *EAR=ear,*GR=grain,*RT=root,*ND=bacteria in root nodule,canopy
C
      READ(11,*)DMLF(NZ,NY,NX),DMSHE(NZ,NY,NX),DMSTK(NZ,NY,NX)
     2,DMRSV(NZ,NY,NX),DMHSK(NZ,NY,NX),DMEAR(NZ,NY,NX)
     3,DMGR(NZ,NY,NX),DMRT(NZ,NY,NX),DMND(NZ,NY,NX)
C
C     ORGAN N AND P CONCENTRATIONS
C
C     CN*,CP*=N:C,P:C ratios in plant organs
C     *LF=leaf,*SHE=petiole,*STK=stalk,*RSV=stalk reserve,*HSK=husk
C     *EAR=ear,*GR=grain,*RT=root,*ND=bacteria in root nodule,canopy
C
      READ(11,*)CNLF(NZ,NY,NX),CNSHE(NZ,NY,NX),CNSTK(NZ,NY,NX)
     2,CNRSV(NZ,NY,NX),CNHSK(NZ,NY,NX),CNEAR(NZ,NY,NX)
     3,CNGR(NZ,NY,NX),CNRT(NZ,NY,NX),CNND(NZ,NY,NX)
      READ(11,*)CPLF(NZ,NY,NX),CPSHE(NZ,NY,NX),CPSTK(NZ,NY,NX)
     2,CPRSV(NZ,NY,NX),CPHSK(NZ,NY,NX),CPEAR(NZ,NY,NX)
     3,CPGR(NZ,NY,NX),CPRT(NZ,NY,NX),CPND(NZ,NY,NX)
C
C     RE-CALCULATE PLANT INPUTS IN MODEL UNITS
C
C     ABSR,ABSP=leaf SW,PAR absorbtivity
C
C     2.5 converts g enzyme mass into g enzyme C
C
      VCMX(NZ,NY,NX)=2.5*VCMX(NZ,NY,NX)
      VOMX(NZ,NY,NX)=2.5*VOMX(NZ,NY,NX)
      VCMX4(NZ,NY,NX)=2.5*VCMX4(NZ,NY,NX)
      ETMX(NZ,NY,NX)=2.5*ETMX(NZ,NY,NX)
      ABSR(NZ,NY,NX)=1.0-ALBR(NZ,NY,NX)-TAUR(NZ,NY,NX)
      ABSP(NZ,NY,NX)=1.0-ALBP(NZ,NY,NX)-TAUP(NZ,NY,NX)
      ALBR(NZ,NY,NX)=ALBR(NZ,NY,NX)/ABSR(NZ,NY,NX)
      ALBP(NZ,NY,NX)=ALBP(NZ,NY,NX)/ABSP(NZ,NY,NX)
      TAUR(NZ,NY,NX)=TAUR(NZ,NY,NX)/ABSR(NZ,NY,NX)
      TAUP(NZ,NY,NX)=TAUP(NZ,NY,NX)/ABSP(NZ,NY,NX)
      ANGBR(NZ,NY,NX)=SIN(ANGBR(NZ,NY,NX)/57.29578)
      ANGSH(NZ,NY,NX)=SIN(ANGSH(NZ,NY,NX)/57.29578)
      GROUPI(NZ,NY,NX)=GROUPX(NZ,NY,NX)
      IF(IBTYP(NZ,NY,NX).NE.0)THEN
      XRNI(NZ,NY,NX)=XRNI(NZ,NY,NX)/25.0
      XRLA(NZ,NY,NX)=XRLA(NZ,NY,NX)/25.0
      GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)/25.0
      XTLI(NZ,NY,NX)=XTLI(NZ,NY,NX)/25.0
      ENDIF
      GROUPI(NZ,NY,NX)=GROUPI(NZ,NY,NX)-XTLI(NZ,NY,NX)
      IF(XDL(NZ,NY,NX).LT.0.0)THEN
      XDL(NZ,NY,NX)=DYLM(NY,NX)
      ENDIF
      DO 5 NB=1,JB
      IF(IWTYP(NZ,NY,NX).EQ.0.AND.ISTYP(NZ,NY,NX).NE.0)THEN
      VRNL(NB,NZ,NY,NX)=AMIN1(4380.0,VRNLI
     2+144.0*ZTYPI(NZ,NY,NX)*(NB-1))
      VRNX(NB,NZ,NY,NX)=AMIN1(4380.0,VRNXI
     2+144.0*ZTYPI(NZ,NY,NX)*(NB-1))
      ELSE
      VRNL(NB,NZ,NY,NX)=VRNLI
      VRNX(NB,NZ,NY,NX)=VRNXI
      ENDIF
5     CONTINUE         
      CLOSE(11)
      ENDIF
C     WRITE(*,1111)'CRITICAL DAYLENGTH',IGO,NZ,XDL(NZ,NY,NX)
1111  FORMAT(A20,2I8,E12.4)
C
C     INITIALIZE INPUTS FOR PLANT MANAGEMENT 
C
C     IHVST=harvest type
C     JHVST=termination type
C     HVST=harvest management
C     THIN=plant removal during harvest
C     EHVST=harvest partitioning fractions
C     IDAY0,IYR0=date,year of planting
C     IDAYH,IYRH=date,year of harvesting 
C
      DO 15 M=1,366
      IHVST(NZ,M,NY,NX)=-1
      JHVST(NZ,M,NY,NX)=0
      HVST(NZ,M,NY,NX)=1.0E+06
      THIN(NZ,M,NY,NX)=-1.0
      EHVST(1,1,NZ,M,NY,NX)=1.0
      EHVST(1,2,NZ,M,NY,NX)=1.0
      EHVST(1,3,NZ,M,NY,NX)=1.0
      EHVST(1,4,NZ,M,NY,NX)=1.0
      EHVST(2,1,NZ,M,NY,NX)=1.0
      EHVST(2,2,NZ,M,NY,NX)=1.0
      EHVST(2,3,NZ,M,NY,NX)=1.0
      EHVST(2,4,NZ,M,NY,NX)=1.0
15    CONTINUE
      IF(DATAM(NZ,NY,NX).NE.'NO')THEN
      N=0
      NN=0
      IDAY0(NZ,NY,NX)=-1E+06
      IYR0(NZ,NY,NX)=-1E+06
      IDAYH(NZ,NY,NX)=1E+06
      IYRH(NZ,NY,NX)=1E+06
10    IF(N.EQ.0)THEN
C
C     READ INPUTS FOR PLANT MANAGEMENT IN DATAM FROM �ROUTQ.F�
C
C     PLANTING
C
C     DY=planting date DDMMYYYY
C     PPI=initial planting density (m-2)
C     SDPTHI=seeding depth (m)
C
      READ(12,*,END=540)DY,PPI(NZ,NY,NX),SDPTHI(NZ,NY,NX)
      SDPTHI(NZ,NY,NX)=AMAX1(1.0E-06,SDPTHI(NZ,NY,NX))
      ELSE
C
C     HARVESTING
C
C     DY=if harvesting: harvesting date DDMMYYYY
C       =if grazing or insect: first grazing date DDMMYYYY,
C        followed by last grazing date in next line
C     ICUT=harvest type:0=none
C                      :1=grain
C                      :2=all above-ground
C                      :3=pruning
C                      :4=animal grazing
C                      :5=not used
C                      :6=insect grazing
C     JCUT=terminate PFT:0=no
C                       :1=yes
C                       :2=yes,but reseed for the next growing season
C     HCUT=if harvesting:>0=cutting height (m)
C                       :<0=fraction of LAI removed
C         =if grazing   :grazer biomass (g FW m-2)
C     PCUT=if harvesting:thinning,fraction of population removed
C         =if grazing   :grazer consumption rate (g DW g FW-1 d-1)
C     ECUT11,ECUT12,ECUT13,ECUT14=fraction of harvest removed from 
C       plant as 1:leaf, 2:non-foliar, 3:woody, 4:standing dead (sum=1) 
C     ECUT21,ECUT22,ECUT23,ECUT24=fraction of ECUT11,ECUT12,ECUT13,
C        ECUT14 removed from ecosystem as biomass, CO2
C        ECUT1 � ECUT2 = fraction of harvest returned as litter
C
      READ(12,*,END=540)DY,ICUT,JCUT,HCUT,PCUT
     2,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23,ECUT24
      IF(JCUT.EQ.1)PCUT=1.0
      ENDIF
C
C     DERIVE PLANTING,HARVESTING DATES,YEARS
C
C     IDAY0,IYR0,IDAYH,IYRH=day,year of planting,harvesting
C
      IDX=INT(DY/1.0E+06)
      IMO=INT(DY/1.0E+04-IDX*1.0E+02)
      IYR=INT(DY-(IDX*1.0E+06+IMO*1.0E+04))
      LPY=0
      IF(MOD(IYR,4))520,510,520
510   IF(IMO.GT.2)LPY=1
520   IF(IMO.EQ.1)GO TO 525
      IDY=30*(IMO-1)+ICOR(IMO-1)+IDX+LPY
      GO TO 527
525   IDY=IDX
527   IF(N.EQ.0)THEN
      IF(IDY.GT.0.AND.IYR.GT.0)THEN
C     IDY=IDY-0.5*(NTX-1)
      IDAY0(NZ,NY,NX)=IDY
C     IF(ISTYP(NZ,NY,NX).EQ.0)THEN
      IYR=IYR+(NT-1)*NF+(NTX-1)*NFX-NTZX 
C     ENDIF
      IYR0(NZ,NY,NX)=MIN(IYR,IYRC)
      IDAYX(NZ,NY,NX)=IDAY0(NZ,NY,NX)
      IYRX(NZ,NY,NX)=IYR0(NZ,NY,NX)
      PPZ(NZ,NY,NX)=PPI(NZ,NY,NX)
      ENDIF
      ELSE
      IF(IDY.GT.0.AND.JCUT.EQ.1)THEN
C     IDY=IDY+0.5*(NTX-1)
      IDAYH(NZ,NY,NX)=IDY
      IYR=IYR+(NT-1)*NF+(NTX-1)*NFX-NTZX 
      IYRH(NZ,NY,NX)=MIN(IYR,IYRC)
      ENDIF
C
C     LOAD MODEL ARRAYS WITH PFT MANAGEMENT FOR SPECIFIED DATES
C     FOR USE IN �GROSUB.F�
C
C     IHVST=harvest type
C     JHVST=termination type
C     HVST=harvest management
C     THIN=plant removal during harvest
C     EHVST=(1,:fraction of harvest removed from plant
C           (2,:fraction of EHVST(1, removed from ecosystem
C     EHVST(1,- EHVST(2,: fraction of harvest returned as litter

      IHVST(NZ,IDY,NY,NX)=ICUT
      JHVST(NZ,IDY,NY,NX)=JCUT
      HVST(NZ,IDY,NY,NX)=HCUT
      THIN(NZ,IDY,NY,NX)=PCUT
      EHVST(1,1,NZ,IDY,NY,NX)=ECUT11
      EHVST(1,2,NZ,IDY,NY,NX)=ECUT12
      EHVST(1,3,NZ,IDY,NY,NX)=ECUT13
      EHVST(1,4,NZ,IDY,NY,NX)=ECUT14
      EHVST(2,1,NZ,IDY,NY,NX)=ECUT21
      EHVST(2,2,NZ,IDY,NY,NX)=ECUT22
      EHVST(2,3,NZ,IDY,NY,NX)=ECUT23
      EHVST(2,4,NZ,IDY,NY,NX)=ECUT24
C
C     IF ANIMAL OR INSECT GRAZING LOAD GRAZING PARAMETERS FOR ALL DATES 
C     BETWEEN FIRST AND LAST
C
      IF(IHVST(NZ,IDY,NY,NX).EQ.4.OR.IHVST(NZ,IDY,NY,NX).EQ.6)THEN
      NN=NN+1
      IF(MOD(NN,2))570,560,570
560   IDYE=IDY
      DO 580 IDYG=IDYS+1,IDYE-1
      IHVST(NZ,IDYG,NY,NX)=ICUT
      JHVST(NZ,IDYG,NY,NX)=JCUT
      HVST(NZ,IDYG,NY,NX)=HCUT
      THIN(NZ,IDYG,NY,NX)=PCUT
      EHVST(1,1,NZ,IDYG,NY,NX)=ECUT11
      EHVST(1,2,NZ,IDYG,NY,NX)=ECUT12
      EHVST(1,3,NZ,IDYG,NY,NX)=ECUT13
      EHVST(1,4,NZ,IDYG,NY,NX)=ECUT14
      EHVST(2,1,NZ,IDYG,NY,NX)=ECUT21
      EHVST(2,2,NZ,IDYG,NY,NX)=ECUT22
      EHVST(2,3,NZ,IDYG,NY,NX)=ECUT23
      EHVST(2,4,NZ,IDYG,NY,NX)=ECUT24      
580   CONTINUE
570   IDYS=IDY
      ENDIF
      ENDIF
      N=N+1
      GO TO 10
540   CLOSE(12)
      ELSE
      SDPTHI(NZ,NY,NX)=1.0E-06
      ENDIF
      IDAYY(NZ,NY,NX)=IDAYH(NZ,NY,NX)
      IYRY(NZ,NY,NX)=IYRH(NZ,NY,NX)
      IFLGC(NZ,NY,NX)=0
      IDTH(NZ,NY,NX)=0
      IFLGN(NZ,NY,NX)=0
9985  CONTINUE
9990  CONTINUE
9995  CONTINUE
      RETURN
      END

