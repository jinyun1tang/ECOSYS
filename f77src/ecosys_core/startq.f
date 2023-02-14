      SUBROUTINE startq(NHWQ,NHEQ,NVNQ,NVSQ,NZ1Q,NZ2Q)
C
C     THIS SUBROUTINE INITIALIZES ALL PLANT VARIABLES
C
      include "parameters.h"
      include "filec.h"
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
      include "blk14.h"
      include "blk16.h"
      include "blk18a.h"
      include "blk18b.h"
      CHARACTER*16 DATA(30),DATAP(JP,JY,JX)
      CHARACTER*3 CHOICE(102,20)
      CHARACTER*8 CDATE
      DIMENSION CNOPC(4),CPOPC(4)
C
C     VSTK=stalk volume:mass (g m-3)
C     FDMPM=minimum canopy dry matter concentration (g g-1)
C     FARS=stalk sapwood thickness (m)
C
      DATA VSTK,FDMPM,FARS/4.0E-06,0.16,0.01/
C
C     INITIALIZE SHOOT GROWTH VARIABLES
C
C     IFLGN=0:take values from readq.f
C          =1:retain existing values
C     IFLGC=PFT flag:0=not alive,1=alive
C     IYR0,IDAY0,IYRH,IDAYH=year,day of planting,harvesting
C     PPI,PPX=initial,current population (m-2)
C     WTRVC,WTRVN,WTRVP=initial seed stotage C,N,P stocks (g)
C     CF,CFI=current,initial clumping factor
C     RSMH=cuticular resistance to water (h m-1)
C     RCMX=cuticular resistance to CO2 (s m-1)
C     CNWS,CPWS=protein:N,protein:P ratios
C     CWSRT=maximum root protein concentration (g g-1)
C     O2I=intercellular O2 concentration in C3,C4 PFT (umol mol-1)
C
      DO 9995 NX=NHWQ,NHEQ
      DO 9990 NY=NVNQ,NVSQ
      NZ2X=MIN(NZ2Q,NP(NY,NX)) 
      DO 9985 NZ=NZ1Q,NZ2X
      IF(IFLGN(NZ,NY,NX).EQ.0)THEN
      IYR0(NZ,NY,NX)=IYRX(NZ,NY,NX)
      IDAY0(NZ,NY,NX)=IDAYX(NZ,NY,NX)
      IYRH(NZ,NY,NX)=IYRY(NZ,NY,NX)
      IDAYH(NZ,NY,NX)=IDAYY(NZ,NY,NX)
      ENDIF
      IF(IYR0(NZ,NY,NX).EQ.IYRC)THEN
      IF(PPI(NZ,NY,NX).GT.0.0)THEN
      PP(NZ,NY,NX)=PPI(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      PPZ(NZ,NY,NX)=PPI(NZ,NY,NX)
      ELSE
      PP(NZ,NY,NX)=PPZ(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      ENDIF
      ELSE
      PP(NZ,NY,NX)=PPX(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      ENDIF
      DPP(NZ,NY,NX)=PP(NZ,NY,NX)
C     WRITE(*,3232)'STARTQ',IYRC,NX,NY,NZ
C    2,IDAY0(NZ,NY,NX),IYR0(NZ,NY,NX) 
C    3,IDAYH(NZ,NY,NX),IYRH(NZ,NY,NX) 
C    4,IYRC,IDAYX(NZ,NY,NX),IDAYY(NZ,NY,NX)
C    5,IYRX(NZ,NY,NX),IYRY(NZ,NY,NX),IFLGC(NZ,NY,NX)
C    5,PP(NZ,NY,NX),PPZ(NZ,NY,NX),PPI(NZ,NY,NX),PPX(NZ,NY,NX)
3232  FORMAT(A8,14I8,20E12.4)
      IF(IFLGC(NZ,NY,NX).EQ.0)THEN
      WTRVX(NZ,NY,NX)=GRDM(NZ,NY,NX)*PP(NZ,NY,NX)
      WTRVC(NZ,NY,NX)=WTRVX(NZ,NY,NX)
      WTRVN(NZ,NY,NX)=CNGR(NZ,NY,NX)*WTRVC(NZ,NY,NX)
      WTRVP(NZ,NY,NX)=CPGR(NZ,NY,NX)*WTRVC(NZ,NY,NX)
      PPX(NZ,NY,NX)=PPI(NZ,NY,NX)
      CF(NZ,NY,NX)=CFI(NZ,NY,NX)
C     IF(DATAP(NZ,NY,NX).NE.'NO')THEN
      RSMH(NZ,NY,NX)=RSMX(NZ,NY,NX)/3600.0
      RCMX(NZ,NY,NX)=RSMX(NZ,NY,NX)*1.56
      CNWS(NZ,NY,NX)=2.5
      CPWS(NZ,NY,NX)=25.0
      CWSRT(NZ,NY,NX)=AMIN1(CNRT(NZ,NY,NX)*CNWS(NZ,NY,NX)
     2,CPRT(NZ,NY,NX)*CPWS(NZ,NY,NX))
      IF(ICTYP(NZ,NY,NX).EQ.3)THEN
      O2I(NZ,NY,NX)=2.10E+05
      ELSE
      O2I(NZ,NY,NX)=3.96E+05
      ENDIF
C
C     FRACTIONS OF PLANT LITTER ALLOCATED TO KINETIC COMPONENTS
C     PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
C
C     CFOPC=fraction of plant litter allocated in nonstructural(0,*),
C        foliar(1,*),non-foliar(2,*),stalk(3,*),root(4,*), 
C        coarse woody (5,*)
C
C     NONSTRUCTURAL
C
      CFOPC(0,1,NZ,NY,NX)=0.00
      CFOPC(0,2,NZ,NY,NX)=1.00
      CFOPC(0,3,NZ,NY,NX)=0.00
      CFOPC(0,4,NZ,NY,NX)=0.00
C
C     NON-VASCULAR (E.G. MOSSES)
C
      IF(IGTYP(NZ,NY,NX).EQ.0)THEN
      CFOPC(1,1,NZ,NY,NX)=0.07
      CFOPC(1,2,NZ,NY,NX)=0.25
      CFOPC(1,3,NZ,NY,NX)=0.30
      CFOPC(1,4,NZ,NY,NX)=0.38
      CFOPC(2,1,NZ,NY,NX)=0.07
      CFOPC(2,2,NZ,NY,NX)=0.25
      CFOPC(2,3,NZ,NY,NX)=0.30
      CFOPC(2,4,NZ,NY,NX)=0.38
C
C     LEGUMES
C
      ELSEIF(INTYP(NZ,NY,NX).NE.0)THEN
      CFOPC(1,1,NZ,NY,NX)=0.16
      CFOPC(1,2,NZ,NY,NX)=0.38
      CFOPC(1,3,NZ,NY,NX)=0.34
      CFOPC(1,4,NZ,NY,NX)=0.12
      CFOPC(2,1,NZ,NY,NX)=0.07
      CFOPC(2,2,NZ,NY,NX)=0.41
      CFOPC(2,3,NZ,NY,NX)=0.37
      CFOPC(2,4,NZ,NY,NX)=0.15
C
C     ANNUALS, GRASSES, SHRUBS
C
      ELSEIF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      CFOPC(1,1,NZ,NY,NX)=0.08
      CFOPC(1,2,NZ,NY,NX)=0.41
      CFOPC(1,3,NZ,NY,NX)=0.36
      CFOPC(1,4,NZ,NY,NX)=0.15
      CFOPC(2,1,NZ,NY,NX)=0.07
      CFOPC(2,2,NZ,NY,NX)=0.41
      CFOPC(2,3,NZ,NY,NX)=0.36
      CFOPC(2,4,NZ,NY,NX)=0.16
C
C     DECIDUOUS TREES
C
      ELSEIF(IBTYP(NZ,NY,NX).EQ.1.OR.IBTYP(NZ,NY,NX).GE.3)THEN
      CFOPC(1,1,NZ,NY,NX)=0.07
      CFOPC(1,2,NZ,NY,NX)=0.34
      CFOPC(1,3,NZ,NY,NX)=0.36
      CFOPC(1,4,NZ,NY,NX)=0.23
      CFOPC(2,1,NZ,NY,NX)=0.000
      CFOPC(2,2,NZ,NY,NX)=0.045
      CFOPC(2,3,NZ,NY,NX)=0.660
      CFOPC(2,4,NZ,NY,NX)=0.295
C
C     CONIFEROUS TREES
C
      ELSE
      CFOPC(1,1,NZ,NY,NX)=0.07
      CFOPC(1,2,NZ,NY,NX)=0.25
      CFOPC(1,3,NZ,NY,NX)=0.38
      CFOPC(1,4,NZ,NY,NX)=0.30
      CFOPC(2,1,NZ,NY,NX)=0.000
      CFOPC(2,2,NZ,NY,NX)=0.045
      CFOPC(2,3,NZ,NY,NX)=0.660
      CFOPC(2,4,NZ,NY,NX)=0.295
      ENDIF
C
C     FRACTIONS OF STALK LITTER ALLOCATED TO
C     PROTEIN, CH2O, CELLULOSE, LIGNIN
C
C     NON-VASCULAR
C
      IF(IGTYP(NZ,NY,NX).EQ.0)THEN
      CFOPC(3,1,NZ,NY,NX)=0.07
      CFOPC(3,2,NZ,NY,NX)=0.25
      CFOPC(3,3,NZ,NY,NX)=0.30
      CFOPC(3,4,NZ,NY,NX)=0.38
C
C     ANNUALS, GRASSES, SHRUBS
C
      ELSEIF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      CFOPC(3,1,NZ,NY,NX)=0.03
      CFOPC(3,2,NZ,NY,NX)=0.25
      CFOPC(3,3,NZ,NY,NX)=0.57
      CFOPC(3,4,NZ,NY,NX)=0.15
C
C     DECIDUOUS AND CONIFEROUS TREES
C
      ELSE
      CFOPC(3,1,NZ,NY,NX)=0.00
      CFOPC(3,2,NZ,NY,NX)=0.045
      CFOPC(3,3,NZ,NY,NX)=0.660
      CFOPC(3,4,NZ,NY,NX)=0.295
      ENDIF
C
C     FRACTIONS OF FINE ROOT LITTER ALLOCATED TO
C     PROTEIN, CH2O, CELLULOSE, LIGNIN PC&E 25:601-608
C
C     NON-VASCULAR
C
      IF(IGTYP(NZ,NY,NX).EQ.0)THEN
      CFOPC(4,1,NZ,NY,NX)=0.07
      CFOPC(4,2,NZ,NY,NX)=0.25
      CFOPC(4,3,NZ,NY,NX)=0.30
      CFOPC(4,4,NZ,NY,NX)=0.38
C
C     ANNUALS, GRASSES, SHRUBS
C
      ELSEIF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      CFOPC(4,1,NZ,NY,NX)=0.057
      CFOPC(4,2,NZ,NY,NX)=0.263
      CFOPC(4,3,NZ,NY,NX)=0.542
      CFOPC(4,4,NZ,NY,NX)=0.138
C
C     DECIDUOUS TREES
C
      ELSEIF(IBTYP(NZ,NY,NX).EQ.1.OR.IBTYP(NZ,NY,NX).GE.3)THEN
      CFOPC(4,1,NZ,NY,NX)=0.059
      CFOPC(4,2,NZ,NY,NX)=0.308
      CFOPC(4,3,NZ,NY,NX)=0.464
      CFOPC(4,4,NZ,NY,NX)=0.169
C
C     CONIFEROUS TREES
C
      ELSE
      CFOPC(4,1,NZ,NY,NX)=0.059
      CFOPC(4,2,NZ,NY,NX)=0.308
      CFOPC(4,3,NZ,NY,NX)=0.464
      CFOPC(4,4,NZ,NY,NX)=0.169
      ENDIF
C
C     COARSE WOODY LITTER FROM BOLES AND ROOTS
C
      CFOPC(5,1,NZ,NY,NX)=0.00
      CFOPC(5,2,NZ,NY,NX)=0.045
      CFOPC(5,3,NZ,NY,NX)=0.660
      CFOPC(5,4,NZ,NY,NX)=0.295
C
C     INITIALIZE C-N AND C-P RATIOS IN PLANT LITTER
C
C     CNOPC,CPOPC=fractions to allocate N,P to kinetic components
C     CFOPN,CFOPP=distribution of litter N,P to kinetic components
C
      CNOPC(1)=0.020
      CNOPC(2)=0.010
      CNOPC(3)=0.010
      CNOPC(4)=0.020
      CPOPC(1)=0.0020
      CPOPC(2)=0.0010
      CPOPC(3)=0.0010
      CPOPC(4)=0.0020
      DO 110 N=0,5
      CNOPCT=0.0
      CPOPCT=0.0
      DO 100 M=1,4
      CNOPCT=CNOPCT+CFOPC(N,M,NZ,NY,NX)*CNOPC(M)
      CPOPCT=CPOPCT+CFOPC(N,M,NZ,NY,NX)*CPOPC(M)
100   CONTINUE
      DO 105 M=1,4
      CFOPN(N,M,NZ,NY,NX)=CFOPC(N,M,NZ,NY,NX)*CNOPC(M)/CNOPCT
      CFOPP(N,M,NZ,NY,NX)=CFOPC(N,M,NZ,NY,NX)*CPOPC(M)/CPOPCT
105   CONTINUE
110   CONTINUE
C
C     CONCURRENT NODE GROWTH
C
C     FNOD=scales node number for perennial vegetation (e.g. trees)
C     NNOD=number of concurrently growing nodes
C
      IF(IBTYP(NZ,NY,NX).EQ.0.OR.IGTYP(NZ,NY,NX).LE.1)THEN
      FNOD(NZ,NY,NX)=1.0
      IF(GROUPI(NZ,NY,NX).LE.10)THEN
      NNOD(NZ,NY,NX)=3
      ELSEIF(GROUPI(NZ,NY,NX).LE.15)THEN
      NNOD(NZ,NY,NX)=4
      ELSE
      NNOD(NZ,NY,NX)=5
      ENDIF
      ELSE
      FNOD(NZ,NY,NX)=AMAX1(1.0,0.04/XRLA(NZ,NY,NX))
      NNOD(NZ,NY,NX)=24
      ENDIF
C
C     PFT THERMAL ACCLIMATION
C
C     TCZD,TCXC=base threshold temperature for leafout,leafoff (oC)
C     ZTYP,ZTYPI=dynamic,initial thermal adaptation zone from PFT file
C     OFFST=shift in Arrhenius curve for thermal adaptation (oC)
C     TCZ,TCX=threshold temperature for leafout,leafoff (oC)
C     HTC=high temperature threshold for grain number loss (oC)
C     SSTX=sensitivity to HTC (seeds oC-1 h-1 above HTC)
C     ICTYP=3:C3,=4:C4 from plant species file
C
      TCZD=5.00
      TCXD=7.50
      ZTYP(NZ,NY,NX)=ZTYPI(NZ,NY,NX)
      OFFST(NZ,NY,NX)=2.50*(3.0-ZTYP(NZ,NY,NX))
      TCZ(NZ,NY,NX)=TCZD-OFFST(NZ,NY,NX)
      TCX(NZ,NY,NX)=AMIN1(15.0,TCXD-OFFST(NZ,NY,NX))
      IF(ICTYP(NZ,NY,NX).EQ.3)THEN
      IF(DATAP(NZ,NY,NX)(1:4).EQ.'soyb')THEN
      HTC(NZ,NY,NX)=30.0+3.0*ZTYP(NZ,NY,NX)
      SSTX(NZ,NY,NX)=0.002
      ELSE
      HTC(NZ,NY,NX)=27.0+3.0*ZTYP(NZ,NY,NX)
      SSTX(NZ,NY,NX)=0.002
      ENDIF
      ELSE
      HTC(NZ,NY,NX)=27.0+3.0*ZTYP(NZ,NY,NX)
      SSTX(NZ,NY,NX)=0.005
      ENDIF
C
C     SEED CHARACTERISTICS
C
C     SDVL,SDLG,SDAR=seed volume(m3),length(m),area(m2)
C     GRDM=seed C mass (g) from PFT file
C
      SDVL(NZ,NY,NX)=GRDM(NZ,NY,NX)*5.0E-06
      SDLG(NZ,NY,NX)=2.0*(0.75*SDVL(NZ,NY,NX)/3.1416)**0.33
      SDAR(NZ,NY,NX)=4.0*3.1416*(SDLG(NZ,NY,NX)/2.0)**2
C
C     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) DIMENSIONS, 
C     UPTAKE PARAMETERS
C
C     SDPTH=seeding depth(m) from PFT management file
C     CDPTHZ=depth to soil layer bottom from surface(m)
C     NG,NIX,NINR=seeding,upper,lower rooting layer
C     CNRTS,CPRTS=N,P root growth yield
C     RRAD1M,RRAD2M=maximum primary,secondary mycorrhizal radius (m)
C     PORT=mycorrhizal porosity
C     UPMXZH,UPKMZH,UPMNZH=NH4 max uptake(g m-2 h-1),
C        Km(uM),min concn (uM)      
C     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake(g m-2 h-1),
C        Km(uM), min concn (uM)      
C     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake(g m-2 h-1),
C        Km(uM),min concn (uM)      
C     RSRR,RSRA=radial,axial root resistivity (m2 MPa-1 h-1)
C
      SDPTH(NZ,NY,NX)=SDPTHI(NZ,NY,NX)
      DO 9795 L=NU(NY,NX),NL(NY,NX)
      IF(SDPTH(NZ,NY,NX).GE.CDPTHZ(L-1,NY,NX)
     2.AND.SDPTH(NZ,NY,NX).LT.CDPTHZ(L,NY,NX))THEN
      NG(NZ,NY,NX)=L
      NIX(NZ,NY,NX)=L
      DO 9790 NR=1,10
      NINR(NR,NZ,NY,NX)=L
9790  CONTINUE
      ENDIF
9795  CONTINUE
      CNRTS(NZ,NY,NX)=CNRT(NZ,NY,NX)*DMRT(NZ,NY,NX)
      CPRTS(NZ,NY,NX)=CPRT(NZ,NY,NX)*DMRT(NZ,NY,NX)
      RRAD1M(2,NZ,NY,NX)=2.5E-06
      RRAD2M(2,NZ,NY,NX)=2.5E-06
      PORTI(2,NZ,NY,NX)=PORTI(1,NZ,NY,NX)
      UPMXZH(2,NZ,NY,NX)=UPMXZH(1,NZ,NY,NX)
      UPKMZH(2,NZ,NY,NX)=UPKMZH(1,NZ,NY,NX)
      UPMNZH(2,NZ,NY,NX)=UPMNZH(1,NZ,NY,NX)
      UPMXZO(2,NZ,NY,NX)=UPMXZO(1,NZ,NY,NX)
      UPKMZO(2,NZ,NY,NX)=UPKMZO(1,NZ,NY,NX)
      UPMNZO(2,NZ,NY,NX)=UPMNZO(1,NZ,NY,NX)
      UPMXPO(2,NZ,NY,NX)=UPMXPO(1,NZ,NY,NX)
      UPKMPO(2,NZ,NY,NX)=UPKMPO(1,NZ,NY,NX)
      UPMNPO(2,NZ,NY,NX)=UPMNPO(1,NZ,NY,NX)
      RSRR(2,NZ,NY,NX)=1.0E+04
      RSRA(2,NZ,NY,NX)=1.0E+12
C
C     PORTX=tortuosity for gas transport
C     RRADP=path length for radial diffusion within root (m)
C     DMVL=volume:C ratio (m3 g-1)
C     RTLG1X,RTLG2X=specific primary,secondary root length (m g-1)
C     RTAR1X,RTAR2X=specific primary,secondary root area (m2 g-1)
C
      DO 500 N=1,2
      PORT(N,NZ,NY,NX)=PORTI(N,NZ,NY,NX)
      RRADP(N,NZ,NY,NX)=LOG(1.0/SQRT(AMAX1(0.01,PORTI(N,NZ,NY,NX))))
      DMVL(N,NZ,NY,NX)=1.0E-06/(0.05*(1.0-PORTI(N,NZ,NY,NX)))
      RTLG1X(N,NZ,NY,NX)=DMVL(N,NZ,NY,NX)
     2/(3.142*RRAD1M(N,NZ,NY,NX)**2)
      RTLG2X(N,NZ,NY,NX)=DMVL(N,NZ,NY,NX)
     2/(3.142*RRAD2M(N,NZ,NY,NX)**2)
      RRAD1X(N,NZ,NY,NX)=RRAD1M(N,NZ,NY,NX)
C    2*SQRT(0.25*(1.0-PORTI(N,NZ,NY,NX)))
      RRAD2X(N,NZ,NY,NX)=RRAD2M(N,NZ,NY,NX)
C    2*SQRT(0.25*(1.0-PORTI(N,NZ,NY,NX)))
      RTAR1X(N,NZ,NY,NX)=3.142*RRAD1X(N,NZ,NY,NX)**2
      RTAR2X(N,NZ,NY,NX)=3.142*RRAD2X(N,NZ,NY,NX)**2
500   CONTINUE
C
C     INITIALIZE PLANT PHENOLOGY 
C
      IFLGI(NZ,NY,NX)=0
      IDTHP(NZ,NY,NX)=0
      IDTHR(NZ,NY,NX)=0
      NBT(NZ,NY,NX)=0
      NBR(NZ,NY,NX)=0
      HTCTL(NZ,NY,NX)=0.0
      ZC(NZ,NY,NX)=0.0
      DO 10 NB=1,JB
      IFLGA(NB,NZ,NY,NX)=0
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.IWTYP(NZ,NY,NX).EQ.1)THEN
      IFLGE(NB,NZ,NY,NX)=1
      ELSE
      IFLGE(NB,NZ,NY,NX)=0
      ENDIF
      IFLGF(NB,NZ,NY,NX)=0
      IFLGR(NB,NZ,NY,NX)=0
      IFLGQ(NB,NZ,NY,NX)=0
      GROUP(NB,NZ,NY,NX)=GROUPI(NZ,NY,NX)
      PSTG(NB,NZ,NY,NX)=XTLI(NZ,NY,NX)
      PSTGI(NB,NZ,NY,NX)=PSTG(NB,NZ,NY,NX)
      PSTGF(NB,NZ,NY,NX)=0.0
      VSTG(NB,NZ,NY,NX)=0.0
      VSTGX(NB,NZ,NY,NX)=0.0
      KLEAF(NB,NZ,NY,NX)=1
      KLEAFX(NB,NZ,NY,NX)=1
      KVSTG(NB,NZ,NY,NX)=1
      KVSTGN(NB,NZ,NY,NX)=0
      GSTGI(NB,NZ,NY,NX)=0.0
      GSTGF(NB,NZ,NY,NX)=0.0
      TGSTGI(NB,NZ,NY,NX)=0.0
      TGSTGF(NB,NZ,NY,NX)=0.0
      VRNY(NB,NZ,NY,NX)=0.0
      VRNZ(NB,NZ,NY,NX)=0.0
      VRNS(NB,NZ,NY,NX)=VRNY(NB,NZ,NY,NX)
      VRNF(NB,NZ,NY,NX)=VRNZ(NB,NZ,NY,NX)
      ATRP(NB,NZ,NY,NX)=0.0
      FDBK(NB,NZ,NY,NX)=1.0
      FDBKX(NB,NZ,NY,NX)=1.0
      FLG4(NB,NZ,NY,NX)=0.0
      FLGZ(NB,NZ,NY,NX)=0.0
      NBTB(NB,NZ,NY,NX)=0
      IDTHB(NB,NZ,NY,NX)=1
      DO 15 M=1,10
      IDAY(M,NB,NZ,NY,NX)=0
15    CONTINUE
10    CONTINUE
C
C     INITIALIZE PLANT MORPHOLOGY AND BIOMASS
C
      WSTR(NZ,NY,NX)=0.0
      CHILL(NZ,NY,NX)=0.0
      HEAT(NZ,NY,NX)=0.0
      DO 25 NB=1,JB
      CPOOL(NB,NZ,NY,NX)=0.0
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
      RNH3B(NB,NZ,NY,NX)=0.0
      RCZLX(NB,NZ,NY,NX)=0.0
      RCPLX(NB,NZ,NY,NX)=0.0
      RCCLX(NB,NZ,NY,NX)=0.0
      WGLFX(NB,NZ,NY,NX)=0.0
      WGLFNX(NB,NZ,NY,NX)=0.0
      WGLFPX(NB,NZ,NY,NX)=0.0
      ARLFZ(NB,NZ,NY,NX)=0.0
      RCZSX(NB,NZ,NY,NX)=0.0
      RCPSX(NB,NZ,NY,NX)=0.0
      RCCSX(NB,NZ,NY,NX)=0.0
      WTSTXB(NB,NZ,NY,NX)=0.0
      WTSTXN(NB,NZ,NY,NX)=0.0
      WTSTXP(NB,NZ,NY,NX)=0.0
      WGSHEX(NB,NZ,NY,NX)=0.0
      WGSHNX(NB,NZ,NY,NX)=0.0
      WGSHPX(NB,NZ,NY,NX)=0.0
      HTSHEX(NB,NZ,NY,NX)=0.0
      DO 5 L=1,JC
      ARSTK(L,NB,NZ,NY,NX)=0.0
      DO 5 N=1,4
      SURFB(N,L,NB,NZ,NY,NX)=0.0
      IF(NB.EQ.1)THEN
      SURFD(N,L,NZ,NY,NX)=0.0
      ENDIF
5     CONTINUE
      DO 25 K=0,25
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
      DO 55 L=1,JC
      ARLFL(L,K,NB,NZ,NY,NX)=0.0
      WGLFL(L,K,NB,NZ,NY,NX)=0.0
      WGLFLN(L,K,NB,NZ,NY,NX)=0.0
      WGLFLP(L,K,NB,NZ,NY,NX)=0.0
55    CONTINUE
      IF(K.NE.0)THEN
      CPOOL3(K,NB,NZ,NY,NX)=0.0
      CO2B(K,NB,NZ,NY,NX)=0.0
      HCOB(K,NB,NZ,NY,NX)=0.0
      CPOOL4(K,NB,NZ,NY,NX)=0.0
      DO 45 L=1,JC
      DO 45 N=1,4
      SURF(N,L,K,NB,NZ,NY,NX)=0.0
45    CONTINUE
      ENDIF
25    CONTINUE
      DO 35 L=1,JC
      ARLFV(L,NZ,NY,NX)=0.0
      WGLFV(L,NZ,NY,NX)=0.0
      ARSTV(L,NZ,NY,NX)=0.0
35    CONTINUE
      CPOOLP(NZ,NY,NX)=0.0
      ZPOOLP(NZ,NY,NX)=0.0
      PPOOLP(NZ,NY,NX)=0.0
      CCPOLP(NZ,NY,NX)=0.0
      CCPLNP(NZ,NY,NX)=0.0
      CZPOLP(NZ,NY,NX)=0.0
      CPPOLP(NZ,NY,NX)=0.0
      WTSHT(NZ,NY,NX)=0.0
      DWTSHT(NZ,NY,NX)=0.0
      WTLF(NZ,NY,NX)=0.0
      WTSHE(NZ,NY,NX)=0.0
      WTSTK(NZ,NY,NX)=0.0
      WVSTK(NZ,NY,NX)=0.0
      WTRSV(NZ,NY,NX)=0.0
      WTHSK(NZ,NY,NX)=0.0
      WTEAR(NZ,NY,NX)=0.0
      WTGR(NZ,NY,NX)=0.0
      WTRT(NZ,NY,NX)=0.0
      WTRTS(NZ,NY,NX)=0.0
      WTND(NZ,NY,NX)=0.0
      WTLS(NZ,NY,NX)=0.0
      WTSHN(NZ,NY,NX)=0.0
      WTLFN(NZ,NY,NX)=0.0
      WTSHEN(NZ,NY,NX)=0.0
      WTSTKN(NZ,NY,NX)=0.0
      WTRSVN(NZ,NY,NX)=0.0
      WTHSKN(NZ,NY,NX)=0.0
      WTEARN(NZ,NY,NX)=0.0
      WTGRNN(NZ,NY,NX)=0.0
      WTNDN(NZ,NY,NX)=0.0
      WTSHP(NZ,NY,NX)=0.0
      WTLFP(NZ,NY,NX)=0.0
      WTSHEP(NZ,NY,NX)=0.0
      WTSTKP(NZ,NY,NX)=0.0
      WTRSVP(NZ,NY,NX)=0.0
      WTHSKP(NZ,NY,NX)=0.0
      WTEARP(NZ,NY,NX)=0.0
      WTGRNP(NZ,NY,NX)=0.0
      WTNDP(NZ,NY,NX)=0.0
      ARLFP(NZ,NY,NX)=0.0
      WTRTA(NZ,NY,NX)=0.0
      ARSTP(NZ,NY,NX)=0.0
C
C     INITIALIZE MASS BALANCE CHECKS
C
      IF(DATA(20).EQ.'NO'.AND.IGO.EQ.0)THEN
      CARBN(NZ,NY,NX)=0.0
      TCSN0(NZ,NY,NX)=0.0
      TZSN0(NZ,NY,NX)=0.0
      TPSN0(NZ,NY,NX)=0.0
      TCO2T(NZ,NY,NX)=0.0
      TCO2A(NZ,NY,NX)=0.0
      TCUPTK(NZ,NY,NX)=0.0
      TCSNC(NZ,NY,NX)=0.0
      TZUPTK(NZ,NY,NX)=0.0
      TZSNC(NZ,NY,NX)=0.0
      TPUPTK(NZ,NY,NX)=0.0
      TPSNC(NZ,NY,NX)=0.0
      TZUPFX(NZ,NY,NX)=0.0
      RNH3C(NZ,NY,NX)=0.0
      TNH3C(NZ,NY,NX)=0.0
      VCOXF(NZ,NY,NX)=0.0
      VNOXF(NZ,NY,NX)=0.0
      VPOXF(NZ,NY,NX)=0.0
      THVSTC(NZ,NY,NX)=0.0
      THVSTN(NZ,NY,NX)=0.0
      THVSTP(NZ,NY,NX)=0.0
      HVSTC(NZ,NY,NX)=0.0
      HVSTN(NZ,NY,NX)=0.0
      HVSTP(NZ,NY,NX)=0.0
      RSETC(NZ,NY,NX)=0.0
      RSETN(NZ,NY,NX)=0.0
      RSETP(NZ,NY,NX)=0.0
      CTRAN(NZ,NY,NX)=0.0
C
C     INITIALIZE STANDING DEAD VARIABLES
C
      IF(IFLGN(NZ,NY,NX).EQ.0)THEN
      WTSTG(NZ,NY,NX)=0.0
      WTSTGN(NZ,NY,NX)=0.0
      WTSTGP(NZ,NY,NX)=0.0
      WTSTDX=WTSTDI(NZ,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      DO 155 M=1,4
      WTSTDG(M,NZ,NY,NX)=WTSTDX*CFOPC(5,M,NZ,NY,NX)
      WTSTDN(M,NZ,NY,NX)=WTSTDX*CNSTK(NZ,NY,NX)
     2*CFOPN(5,M,NZ,NY,NX)
      WTSTDP(M,NZ,NY,NX)=WTSTDX*CPSTK(NZ,NY,NX)
     2*CFOPP(5,M,NZ,NY,NX)
      WTSTG(NZ,NY,NX)=WTSTG(NZ,NY,NX)+WTSTDG(M,NZ,NY,NX) 
      WTSTGN(NZ,NY,NX)=WTSTGN(NZ,NY,NX)+WTSTDN(M,NZ,NY,NX) 
      WTSTGP(NZ,NY,NX)=WTSTGP(NZ,NY,NX)+WTSTDP(M,NZ,NY,NX)
155   CONTINUE
      WTSTDG(5,NZ,NY,NX)=0.0
      WTSTDN(5,NZ,NY,NX)=0.0
      WTSTDP(5,NZ,NY,NX)=0.0
C
C     INITIALIZE STALK VARIABLES
C
      ZG(NZ,NY,NX)=0.0
      ARSTG(NZ,NY,NX)=0.0
      ARLSS(NY,NX)=ARLSS(NY,NX)+ARSTG(NZ,NY,NX)
      DO 6 L=1,JC
      XJC=JC
      ARSTD(L,NZ,NY,NX)=ARSTG(NZ,NY,NX)/XJC
6     CONTINUE
C
C     INITIALIZE PLANT HEAT AND WATER STATUS
C
C     TKQC,VPQC=canopy aerodynamic temperature (K), 
C        vapor pressure (kPa)
C     TKQD,VPQD=standing dead aerodynamic temperature (K), 
C        vapor pressure (kPa)
C     TCC,TKC=canopy surface temperature (oC,K)
C     TCD,TKD=standing dead surface temperature (oC,K)
C     TCG,TKG=canopy surface temperature for phenology (oC,K)
C     PSILT,PSILO,PSILG=canopy total,osmotic,turgor water potl(MPa)
C
      ENGYX(NZ,NY,NX)=0.0
      TKQC(NZ,NY,NX)=ATKA(NY,NX)
      VPQC(NZ,NY,NX)=2.173E-03/TKQC(NZ,NY,NX)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKQC(NZ,NY,NX)))
      TKQD(NZ,NY,NX)=ATKA(NY,NX)
      VPQD(NZ,NY,NX)= 2.173E-03/TKQD(NZ,NY,NX)
     2*0.61*EXP(5360.0*(3.661E-03-1.0/TKQD(NZ,NY,NX)))
      TCC(NZ,NY,NX)=ATCA(NY,NX)
      TKC(NZ,NY,NX)=ATKA(NY,NX)
      TCG(NZ,NY,NX)=ATCA(NY,NX)
      TKG(NZ,NY,NX)=ATKA(NY,NX)
      TCD(NZ,NY,NX)=ATCA(NY,NX)
      TKD(NZ,NY,NX)=ATKA(NY,NX)
      TFN3(NZ,NY,NX)=1.0
      PSILT(NZ,NY,NX)=-1.0E-03
      PSILO(NZ,NY,NX)=OSMO(NZ,NY,NX)+PSILT(NZ,NY,NX)
      PSILG(NZ,NY,NX)=AMAX1(0.0,PSILT(NZ,NY,NX)-PSILO(NZ,NY,NX))
      EP(NZ,NY,NX)=0.0
      ENDIF
      ENDIF
C
C     INITIALIZE ROOT(N=1),MYCORRHIZAL(N=2) MORPHOLOGY AND BIOMASS
C
C     PSIRT,PSIRO,PSIRG=root,myco total,osmotic,turgor water potl(MPa)
C     CO2A,CO2P=root,myco gaseous,aqueous CO2 content (g)
C     OXYA,OXYP=root,myco gaseous,aqueous O2 content (g)
C
      NRT(NZ,NY,NX)=0
      UPNH4(NZ,NY,NX)=0.0
      UPNO3(NZ,NY,NX)=0.0
      UPH2P(NZ,NY,NX)=0.0
      UPH1P(NZ,NY,NX)=0.0
      UPNF(NZ,NY,NX)=0.0
      DO 40 N=1,2
      DO 20 L=1,NL(NY,NX)
      UPWTR(N,L,NZ,NY,NX)=0.0
      PSIRT(N,L,NZ,NY,NX)=-0.01
      PSIRO(N,L,NZ,NY,NX)=OSMO(NZ,NY,NX)+PSIRT(N,L,NZ,NY,NX)
      PSIRG(N,L,NZ,NY,NX)=AMAX1(0.0,PSIRT(N,L,NZ,NY,NX)
     2-PSIRO(N,L,NZ,NY,NX))
      CPOOLR(N,L,NZ,NY,NX)=0.0
      ZPOOLR(N,L,NZ,NY,NX)=0.0
      PPOOLR(N,L,NZ,NY,NX)=0.0
      CCPOLR(N,L,NZ,NY,NX)=0.0
      CZPOLR(N,L,NZ,NY,NX)=0.0
      CPPOLR(N,L,NZ,NY,NX)=0.0
      CWSRTL(N,L,NZ,NY,NX)=CWSRT(NZ,NY,NX)
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
      RTLGA(N,L,NZ,NY,NX)=1.0E-03
      RUPNH4(N,L,NZ,NY,NX)=0.0
      RUPNO3(N,L,NZ,NY,NX)=0.0
      RUPH2P(N,L,NZ,NY,NX)=0.0
      RUPH1P(N,L,NZ,NY,NX)=0.0
      RUPNHB(N,L,NZ,NY,NX)=0.0
      RUPNOB(N,L,NZ,NY,NX)=0.0
      RUPH2B(N,L,NZ,NY,NX)=0.0
      RUPH1B(N,L,NZ,NY,NX)=0.0
      ROXYP(N,L,NZ,NY,NX)=0.0
      RUNNHP(N,L,NZ,NY,NX)=0.0
      RUNNBP(N,L,NZ,NY,NX)=0.0
      RUNNOP(N,L,NZ,NY,NX)=0.0
      RUNNXP(N,L,NZ,NY,NX)=0.0
      RUPP2P(N,L,NZ,NY,NX)=0.0
      RUPP1P(N,L,NZ,NY,NX)=0.0
      RUPP2B(N,L,NZ,NY,NX)=0.0
      RUPP1B(N,L,NZ,NY,NX)=0.0
      CCO2A=CCO2EI(NY,NX)
      CCO2P=0.030*EXP(-2.621-0.0317*ATCA(NY,NX))*CO2EI(NY,NX)
      CO2A(N,L,NZ,NY,NX)=CCO2A*RTVLP(N,L,NZ,NY,NX)
      CO2P(N,L,NZ,NY,NX)=CCO2P*RTVLW(N,L,NZ,NY,NX)
      RCOFLA(N,L,NZ,NY,NX)=0.0
      RCODFA(N,L,NZ,NY,NX)=0.0
      RCO2S(N,L,NZ,NY,NX)=0.0
      RCO2P(N,L,NZ,NY,NX)=0.0
      COXYA=COXYE(NY,NX)
      COXYP=0.032*EXP(-6.175-0.0211*ATCA(NY,NX))*OXYE(NY,NX)
      OXYA(N,L,NZ,NY,NX)=COXYA*RTVLP(N,L,NZ,NY,NX)
      OXYP(N,L,NZ,NY,NX)=COXYP*RTVLW(N,L,NZ,NY,NX)
      CH4A(N,L,NZ,NY,NX)=0.0
      CH4P(N,L,NZ,NY,NX)=0.0
      Z2OA(N,L,NZ,NY,NX)=0.0
      Z2OP(N,L,NZ,NY,NX)=0.0
      ZH3A(N,L,NZ,NY,NX)=0.0
      ZH3P(N,L,NZ,NY,NX)=0.0
      H2GA(N,L,NZ,NY,NX)=0.0
      H2GP(N,L,NZ,NY,NX)=0.0
      WFR(N,L,NZ,NY,NX)=1.0
      DO 30 NR=1,10
      RTN2(N,L,NR,NZ,NY,NX)=0.0
      RTLG1(N,L,NR,NZ,NY,NX)=0.0
      WTRT1(N,L,NR,NZ,NY,NX)=0.0
      WTRT1N(N,L,NR,NZ,NY,NX)=0.0
      WTRT1P(N,L,NR,NZ,NY,NX)=0.0
      RTLG2(N,L,NR,NZ,NY,NX)=0.0
      WTRT2(N,L,NR,NZ,NY,NX)=0.0
      WTRT2N(N,L,NR,NZ,NY,NX)=0.0
      WTRT2P(N,L,NR,NZ,NY,NX)=0.0
      RTDP1(N,NR,NZ,NY,NX)=SDPTH(NZ,NY,NX)
      RTWT1(N,NR,NZ,NY,NX)=0.0
      RTWT1N(N,NR,NZ,NY,NX)=0.0
      RTWT1P(N,NR,NZ,NY,NX)=0.0
30    CONTINUE
      IF(N.EQ.1)THEN
      DO 6400 K=0,1
      DO 6400 M=1,5
      CSNC(M,K,L,NZ,NY,NX)=0.0
      ZSNC(M,K,L,NZ,NY,NX)=0.0
      PSNC(M,K,L,NZ,NY,NX)=0.0
6400  CONTINUE
      CPOOLN(L,NZ,NY,NX)=0.0
      ZPOOLN(L,NZ,NY,NX)=0.0
      PPOOLN(L,NZ,NY,NX)=0.0
      WTNDL(L,NZ,NY,NX)=0.0
      WTNDLN(L,NZ,NY,NX)=0.0
      WTNDLP(L,NZ,NY,NX)=0.0
      RUPNF(L,NZ,NY,NX)=0.0
      ENDIF
20    CONTINUE
40    CONTINUE
C
C     INITIALIZE SEED MORPHOLOGY AND BIOMASS
C
C     WTRVC,WTRVN,WTRVP=C,N,P in storage reserves (g)
C     WTLFB,WTLFBN,WTLFBP=C,N,P in leaves (g)
C     WTLSB=C in leaves+petioles (g)
C     FDM-dry matter fraction (g DM C g FM C-1)
C     VOLWP,VOLWC=water volume in,on canopy (m3)
C     CPOOL,ZPOOL,PPOOL=C,N,P in canopy nonstructural pools (g)
C     WTRT1,WTRT1N,WTRT1P=C,N,P in primary root layer (g)
C     RTWT1,RTWT1N,RTWT1P=total C,N,P in primary root (g)
C     WTRTL,WTRTD=total root C mass (g)
C     WSRTL=total root protein C mass (g)
C     CPOOLR,ZPOOLR,PPOOLR=C,N,P in root,myco nonstructural pools (g)
C
      WTLFBN(1,NZ,NY,NX)=CNGR(NZ,NY,NX)*WTLFB(1,NZ,NY,NX)
      WTLFBP(1,NZ,NY,NX)=CPGR(NZ,NY,NX)*WTLFB(1,NZ,NY,NX)
      WTLSB(1,NZ,NY,NX)=WTLFB(1,NZ,NY,NX)+WTSHEB(1,NZ,NY,NX)
      WTLS(NZ,NY,NX)=WTLS(NZ,NY,NX)+WTLSB(1,NZ,NY,NX)
      WVPLT(NZ,NY,NX)=AMAX1(0.0,WTLS(NZ,NY,NX)
     2+AMIN1(WTSTK(NZ,NY,NX),FARS*ARSTP(NZ,NY,NX)/VSTK))
      APSILT=ABS(PSILT(NZ,NY,NX))
      FDMP=FDMPM+0.10*APSILT/(0.05*APSILT+2.0)
      VOLWP(NZ,NY,NX)=1.0E-06*WVPLT(NZ,NY,NX)/FDMP
      VOLWC(NZ,NY,NX)=0.0
      ZPOOL(1,NZ,NY,NX)=CNGR(NZ,NY,NX)*CPOOL(1,NZ,NY,NX)
      PPOOL(1,NZ,NY,NX)=CPGR(NZ,NY,NX)*CPOOL(1,NZ,NY,NX)
      WTRT1N(1,NG(NZ,NY,NX),1,NZ,NY,NX)=CNGR(NZ,NY,NX)
     2*WTRT1(1,NG(NZ,NY,NX),1,NZ,NY,NX)
      WTRT1P(1,NG(NZ,NY,NX),1,NZ,NY,NX)=CPGR(NZ,NY,NX)
     2*WTRT1(1,NG(NZ,NY,NX),1,NZ,NY,NX)
      RTWT1N(1,1,NZ,NY,NX)=CNGR(NZ,NY,NX)*RTWT1(1,1,NZ,NY,NX)
      RTWT1P(1,1,NZ,NY,NX)=CPGR(NZ,NY,NX)*RTWT1(1,1,NZ,NY,NX)
      WTRTL(1,NG(NZ,NY,NX),NZ,NY,NX)=WTRT1(1,NG(NZ,NY,NX),1,NZ,NY,NX)
      WTRTD(1,NG(NZ,NY,NX),NZ,NY,NX)=WTRT1(1,NG(NZ,NY,NX),1,NZ,NY,NX)
      WSRTL(1,NG(NZ,NY,NX),NZ,NY,NX)=WTRTL(1,NG(NZ,NY,NX),NZ,NY,NX)
     2*CWSRT(NZ,NY,NX)
      ZPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=CNGR(NZ,NY,NX)
     2*CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)
      PPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)=CPGR(NZ,NY,NX)
     2*CPOOLR(1,NG(NZ,NY,NX),NZ,NY,NX)
C     ENDIF
      ENDIF
C     IF(IGTYP(NZ,NY,NX).NE.0)THEN
C     ZEROP(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)
C     ZEROQ(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
C     ZEROP2(NZ,NY,NX)=ZERO2*PP(NZ,NY,NX)
C     ELSE
      ZEROP(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)
      ZEROQ(NZ,NY,NX)=ZERO*PP(NZ,NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      ZEROP2(NZ,NY,NX)=ZERO2*PP(NZ,NY,NX)
C     ENDIF
9985  CONTINUE
C
C     FILL OUT UNUSED ARRAYS
C
      DO 9986 NZ=NP(NY,NX)+1,5
      TCSN0(NZ,NY,NX)=0.0
      TZSN0(NZ,NY,NX)=0.0
      TPSN0(NZ,NY,NX)=0.0
      TCSNC(NZ,NY,NX)=0.0
      TZSNC(NZ,NY,NX)=0.0
      TPSNC(NZ,NY,NX)=0.0
      WTSTG(NZ,NY,NX)=0.0
      WTSTGN(NZ,NY,NX)=0.0
      WTSTGP(NZ,NY,NX)=0.0
      ARSTG(NZ,NY,NX)=0.0
      DO 6401 L=1,NL(NY,NX)
      DO 6401 K=0,1
      DO 6401 M=1,5
      CSNC(M,K,L,NZ,NY,NX)=0.0
      ZSNC(M,K,L,NZ,NY,NX)=0.0
      PSNC(M,K,L,NZ,NY,NX)=0.0
6401  CONTINUE
9986  CONTINUE
      DO 9980 NZ=1,NP(NY,NX)
      HCBFCZ(NZ,NY,NX)=0.0 
      HCBFDZ(NZ,NY,NX)=0.0 
9980  CONTINUE
9990  CONTINUE
9995  CONTINUE
      RETURN
      END


