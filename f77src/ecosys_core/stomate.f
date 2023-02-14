
      SUBROUTINE stomate(I,J,NFZ,NZ,NY,NX)
C
C     THIS SUBROUTINE CALCULATES CANOPY STOMATAL RESISTANCE AT MAXIMUM
C     CANOPY TURGOR FOR USE IN ENERGY BALANCE EQUATIONS IN 'UPTAKE'
C
      include "parameters.h"
      include "blkc.h"
      include "blk1cp.h"
      include "blk1g.h"
      include "blk1n.h"
      include "blk1p.h"
      include "blk2a.h"
      include "blk3.h"
      include "blk5.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk9a.h"
      include "blk9b.h"
      include "blk9c.h"
      include "blk1u.h"
      DIMENSION FLG4Y(0:5)
C
C     QNTM=quantum efficiency (umol e- umol-1 PAR)
C     CURV=shape parameter for e- transport response to PAR
C     ELEC3,ELEC4=e- requirement for CO2 fixn by rubisco,PEP 
C     carboxylase (umol e- umol CO2)
C     CNKI,CPKI=nonstruct N,P inhibition constant on rubisco 
C        (g N,P g-1 C)
C     RSMY=minimum stomatal resistance for CO2 uptake (h m-1)
C     ATRPZ=hours to full dehardening of conifers in spring (h)
C     COMP4=C4 CO2 compensation point (uM)
C     FDML=leaf water content (g H2O g-1 C)
C     FBS,FMP=leaf water content in bundle sheath, mesophyll in C4 CO2
C        fixation
C     C4KI=nonstructural C inhibition constant on PEP carboxylase (uM) 
C     FLG4Y=number of hours with no grain fill to terminate annuals
C
      PARAMETER (QNTM=0.45,CURV=0.70,CURV2=2.0*CURV,CURV4=4.0*CURV
     2,ELEC3=4.5,ELEC4=3.0)
      PARAMETER (CNKI=1.0E+02,CPKI=1.0E+03)
      PARAMETER (RSMY=2.78E-03,ATRPZ=276.9)
      PARAMETER (COMP4=0.5,FDML=6.0,FBS=0.2*FDML,FMP=0.8*FDML
     2,C4KI=5.0E+06)
      DATA FLG4Y/336.0,672.0,672.0,672.0,672.0,672.0/
C
C     CANOPY TEMPERATURE + OFFSET FOR THERMAL ADAPTATION FROM 'READQ'
C
C     FMOL=number of moles of air per m3
C     CO2I=intercellular CO2 concentration
C     FCO2=intercellular:atmospheric CO2 concn ratio from PFT file 
C     SSIN=sine of solar angle
C     ARLFP=PFT leaf area
C
      FMOL(NZ,NY,NX)=1.2194E+04/TKC(NZ,NY,NX)
      CO2I(NZ,NY,NX)=FCO2(NZ,NY,NX)*CO2Q(NY,NX)
      IF(SSIN(NY,NX).GT.0.0.AND.ARLFP(NZ,NY,NX)
     2.GT.ZEROP(NZ,NY,NX))THEN
C
C     CO2 AND O2 AQUEOUS SOLUBILITY
C
C     TCC=canopy temperature
C     SCO2,SO2=solubility of CO2,O2 (uM/(umol mol-1))
C     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
C     DCO2=difference between atmospheric and intercellular CO2
C        concentration (umol m-3)
C     FMOL=number of moles of air per m3
C
      SCO2(NZ,NY,NX)=EXP(-2.621-0.0317*TCC(NZ,NY,NX))
      SO2(NZ,NY,NX)=EXP(-6.175-0.0211*TCC(NZ,NY,NX))
      CO2L(NZ,NY,NX)=CO2I(NZ,NY,NX)*SCO2(NZ,NY,NX)
      O2L(NZ,NY,NX)=O2I(NZ,NY,NX)*SO2(NZ,NY,NX)
      DCO2(NZ,NY,NX)=FMOL(NZ,NY,NX)*(CO2Q(NY,NX)-CO2I(NZ,NY,NX))
C
C     ARRHENIUS FUNCTIONS FOR CARBOXYLATION AND OXYGENATION
C
C     TKC,TKCO=canopy temperature,canopy temperature function
C     OFFST=shift in temperature function for thermal adaptation
C     TFN1,TFN2,TFNE=temperature function for carboxylation,
C        oxygenation,e- transport (25 oC =1)
C     8.313,710.0=gas constant,enthalpy
C     197500,222500=energy of high,low temp inactivn(KJ mol-1)
C     65000,60000,43000=activation energy for carboxylation,
C        oxygenation,e- transport 
C
      CH2O=0.0
      TKCO=TKC(NZ,NY,NX)+OFFST(NZ,NY,NX)
      RTK=8.3143*TKCO
      STK=710.0*TKCO
      ACTV=1.0+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFN1=EXP(26.237-65000/RTK)/ACTV
      TFN2=EXP(24.220-60000/RTK)/ACTV
      TFNE=EXP(17.362-43000/RTK)/ACTV
C
C     M-M CONSTANT FOR CARBOXYLATION FROM 'READQ' ADJUSTED FOR
C     TEMPERATURE
C
C     XKCO2L,XKCO2O=Km for rubisco carboxylation without,with O2
C     XKO2L=Km for rubisco oxygenation
C
      XKCO2L(NZ,NY,NX)=XKCO2(NZ,NY,NX)*EXP(16.136-40000/RTK) 
      XKO2L=XKO2(NZ,NY,NX)*EXP(8.067-20000/RTK) 
      XKCO2O(NZ,NY,NX)=XKCO2L(NZ,NY,NX)*(1.0+O2L(NZ,NY,NX)/XKO2L)
C
C     FOR EACH BRANCH
C
      DO 2900 NB=1,NBR(NZ,NY,NX)
C
C     FEEDBACK ON CO2 FIXATION
C
C     IWTYP=phenology type from PFT file
C     VRNS,VRNL=leafout hours,hours required for leafout 
C     VRNF,VRNX=leafoff hours,hours required for leafoff 
C
      IF(IWTYP(NZ,NY,NX).EQ.0.OR.ISTYP(NZ,NY,NX).EQ.0
     2.OR.VRNS(NB,NZ,NY,NX).GE.VRNL(NB,NZ,NY,NX)
     3.OR.VRNF(NB,NZ,NY,NX).LT.VRNX(NB,NZ,NY,NX))THEN
C
C     FEEDBACK ON C3 CARBOXYLATION FROM NON-STRUCTURAL C:N:P
C
C     CCPOLB,CZPOLB,CPPOLB=nonstructural C,N,P concentration in branch
C     FDBK=N,P feedback inhibition on C3 CO2 fixation
C     CNKI,CPKI=nonstructural N,P inhibition constant on rubisco
C        activity 
C
      IF(CCPOLB(NB,NZ,NY,NX).GT.ZERO)THEN
      FDBK(NB,NZ,NY,NX)=AMIN1(CZPOLB(NB,NZ,NY,NX)
     3/(CZPOLB(NB,NZ,NY,NX)+CCPOLB(NB,NZ,NY,NX)/CNKI)
     4,CPPOLB(NB,NZ,NY,NX)
     5/(CPPOLB(NB,NZ,NY,NX)+CCPOLB(NB,NZ,NY,NX)/CPKI))
      ELSE
      FDBK(NB,NZ,NY,NX)=1.0
      ENDIF
C
C     CHILLING 
C
C     CHILL,HEAT=accumulated chilling,heating hours used to limit 
C        CO2 fixation
C
C     FDBK(NB,NZ,NY,NX)=FDBK(NB,NZ,NY,NX)/(1.0+0.25*CHILL(NZ,NY,NX))
      FDBK(NB,NZ,NY,NX)=FDBK(NB,NZ,NY,NX)/(1.0+HEAT(NZ,NY,NX))
C
C     DEHARDENING OF EVERGREENS IN SPRING
C
C     IWTYP=phenology type:0=evergreen,1=cold deciduous,2=drought
C        deciduous,3=cold and drought deciduous
C     IBTYP=turnover rate of above-ground biomass
C              :0,1=fully deciduous,
C              :2=needleleaf evergreen
C              :3=broadleaf evergreen
C              :4=semi-deciduous
C              :5=semi-evergreen
C     ATRP=hours above threshold temperature for dehardening 
C        since leafout
C     ATRPZ=hours to full dehardening of conifers in spring 
C
      IF(IWTYP(NZ,NY,NX).NE.0.AND.IBTYP(NZ,NY,NX).GE.2)THEN
      FDBK(NB,NZ,NY,NX)=FDBK(NB,NZ,NY,NX)*AMAX1(0.0,AMIN1(1.0
     2,ATRP(NB,NZ,NY,NX)/(0.9*ATRPZ)))
      ENDIF
C
C     TERMINATION OF ANNUALS 
C
C     ISTYP=growth habit:0=annual,1=perennial from PFT file
C     FDBKX=termination feedback inhibition on C3 CO2
C     FLG4=number of hours with no grain fill after start of 
C        grain fill
C     FLG4Y=number of hours with no grain fill to terminate annuals
C
      IF(ISTYP(NZ,NY,NX).EQ.0.AND.FLG4(NB,NZ,NY,NX).GT.0.0)THEN
      FDBKX(NB,NZ,NY,NX)=AMAX1(0.0
     2,1.0-FLG4(NB,NZ,NY,NX)/FLG4Y(IWTYP(NZ,NY,NX)))
      ELSE
      FDBKX(NB,NZ,NY,NX)=1.0
      ENDIF
      FDBK(NB,NZ,NY,NX)=FDBK(NB,NZ,NY,NX)*FDBKX(NB,NZ,NY,NX)
C     IF(NFZ.EQ.1.AND.NZ.EQ.1.AND.J.EQ.15)THEN
C     WRITE(*,4242)'FDBK',I,J,NFZ,NZ,NB,IDTHB(NB,NZ,NY,NX)
C    2,FDBK(NB,NZ,NY,NX),VRNS(NB,NZ,NY,NX),VRNF(NB,NZ,NY,NX)
C    3,VRNL(NB,NZ,NY,NX),VRNX(NB,NZ,NY,NX)
C    3,CCPOLB(NB,NZ,NY,NX),CZPOLB(NB,NZ,NY,NX),CPPOLB(NB,NZ,NY,NX)
C    3,FDBKX(NB,NZ,NY,NX),ATRP(NB,NZ,NY,NX),FLG4(NB,NZ,NY,NX)
C    4,FLG4Y(IWTYP(NZ,NY,NX)),ARLFP(NZ,NY,NX)
C    5,TCC(NZ,NY,NX),HEAT(NZ,NY,NX)
4242  FORMAT(A8,6I4,20E12.4)
C     ENDIF
C
C     FOR EACH NODE
C
C     IDTHB=branch life flag:0=living,1=dead
C     ARLF,WGLF,WSLF=leaf area,C,protein C
C     WSDN=leaf protein surficial density
C
      IF(IDTHB(NB,NZ,NY,NX).EQ.0)THEN
      DO 2800 K=1,25
      IF(ARLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WGLF(K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      WSDN=WSLF(K,NB,NZ,NY,NX)/ARLF(K,NB,NZ,NY,NX)
C    2*AMIN1(1.0,100.0*ARLF(K,NB,NZ,NY,NX)/WGLF(K,NB,NZ,NY,NX))
C     IF((I/60)*60.EQ.I.AND.NFZ.EQ.1.AND.J.EQ.15)THEN
C     WRITE(*,2125)'WSDN',I,J,NFZ,NX,NY,NZ,NB,K,WSDN
C    2,WGLF(K,NB,NZ,NY,NX),WSLF(K,NB,NZ,NY,NX)
C    3,ARLF(K,NB,NZ,NY,NX),WGLFN(K,NB,NZ,NY,NX),WGLFP(K,NB,NZ,NY,NX)
C    4,WGLF(K,NB,NZ,NY,NX)/ARLF(K,NB,NZ,NY,NX)
C    5,WSLF(K,NB,NZ,NY,NX)/ARLF(K,NB,NZ,NY,NX)
C    2,AMIN1(1.0,100.0*ARLF(K,NB,NZ,NY,NX)/WGLF(K,NB,NZ,NY,NX))
2125  FORMAT(A8,8I4,12E12.4)
C     ENDIF
      ELSE
      WSDN=0.0
      ENDIF
      IF(WSDN.GT.ZERO)THEN
C
C     C4 PHOTOSYNTHESIS
C
C     ICTYP=photosynthesis type:3=C3,4=C4 from PFT file
C
      IF(ICTYP(NZ,NY,NX).EQ.4)THEN
C
C     FEEDBACK ON C4 CARBOXYLATION FROM C4 NON-STRUCTURAL C
C
C     CC4M,CCBS=C4 nonstructural C concentration in mesophyll,
C        bundle sheath (uM)
C     CPOOL4,CO2B=C4 nonstructural C mass in mesophyll,
C        bundle sheath
C     WGLF=leaf C mass
C     FBS,FMP=leaf water content fractions in bundle sheath, mesophyll
C     FDBK4=N,P feedback inhibition on C4 CO2 fixation 
C
      CC4M=AMAX1(0.0,0.021E+09*CPOOL4(K,NB,NZ,NY,NX)
     2/(WGLF(K,NB,NZ,NY,NX)*FMP))
      CCBS=AMAX1(0.0,0.083E+09*CO2B(K,NB,NZ,NY,NX)
     2/(WGLF(K,NB,NZ,NY,NX)*FBS))
      FDBK4(K,NB,NZ,NY,NX)=1.0/(1.0+CC4M/C4KI)
      FDBK4(K,NB,NZ,NY,NX)=FDBK4(K,NB,NZ,NY,NX)*FDBKX(NB,NZ,NY,NX)
C
C     SURFICIAL DENSITY OF PEPC AND ITS CHLOROPHYLL
C
C     VCDN4=surficial density of PEP carboxylase in mesophyll
C     ETDN4=surficial density of chlorophyll in mesophyll
C     PEPC=fraction of leaf protein in PEP carboxylase
C        from PFT file
C     CHL4=fraction of leaf protein in mesophyll chlorophyll
C        from PFT file
C     WSDN=leaf protein surficial density
C
      VCDN4=PEPC(NZ,NY,NX)*WSDN
      ETDN4=CHL4(NZ,NY,NX)*WSDN
C
C     CO2-LIMITED C4 CARBOXYLATION RATES
C
C     VCGR4,VGRO4=PEP carboxylation rate unlimited,limited by CO2
C     VCMX4=specific PEP carboxylase activity from PFT file
C     TFN1=temperature function for carboxylation
C     VCDN4=surficial density of PEP carboxylase in mesophyll
C     CO2L=intercellular CO2 concentrations (uM)
C     COMP4=C4 CO2 compensation point (uM)
C     XKCO24=Km for VCMX4 from PFT file (uM)
C
      VCGR4(K,NB,NZ,NY,NX)=VCMX4(NZ,NY,NX)*TFN1*VCDN4
      VGRO4(K,NB,NZ,NY,NX)=AMAX1(0.0,VCGR4(K,NB,NZ,NY,NX)
     2*(CO2L(NZ,NY,NX)-COMP4)/(CO2L(NZ,NY,NX)+XKCO24(NZ,NY,NX)))
C
C     C4 ELECTRON TRANSFER RATES
C
C     ETGR4=light saturated e- transport rate
C     ETMX=specific chlorophyll activity from PFT file
C     TFNE=temperature function for e- transport 
C     ETDN4=surficial density of chlorophyll in mesophyll
C     CBXN4=PEP caboxylation efficiency
C     CO2L=intercellular CO2 concentrations (uM)
C     COMP4=C4 CO2 compensation point (uM)
C     ELEC4=e- requirement for CO2 fixation by PEP carboxylase
C
      ETGR4(K,NB,NZ,NY,NX)=ETMX(NZ,NY,NX)*TFNE*ETDN4
      CBXN4(K,NB,NZ,NY,NX)=AMAX1(0.0,(CO2L(NZ,NY,NX)-COMP4)
     2/(ELEC4*CO2L(NZ,NY,NX)+10.5*COMP4))
C
C     FOR EACH CANOPY LAYER
C
C     ARLFL=leaf area
C     SURFX=unself-shaded leaf surface area
C
      DO 2700 L=JC,1,-1
      IF(ARLFL(L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C
C     FOR EACH INCLINATION AND AZIMUTH CLASS
C
      DO 2600 N=1,4
      DO 2500 M=1,4
C
C     SURFX=unself-shaded leaf surface area
C     PAR=direct PAR flux from ‘hour1.f’
C
      IF(SURFX(N,L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C
C     SUNLIT LEAVES
C
      IF(PAR(N,M,L,NZ,NY,NX).GT.0.0)THEN
C
C     LIGHT-LIMITED CARBOXYLATION RATES
C
C     QNTM=quantum efficiency 
C     PAR=direct PAR flux from ‘hour1.f’
C     ETGR4=light saturated e- transport rate
C     ETLF4=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO4=light-limited PEP carboxylation rate
C
      PARX=QNTM*PAR(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGR4(K,NB,NZ,NY,NX)
      ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ,NY,NX)))/CURV2
      EGRO4=ETLF4*CBXN4(K,NB,NZ,NY,NX)
C
C     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
C
C     VL=PEP carboxylation rate limited by light,CO2,N,P
C     VGRO4=PEP carboxylation rate limited by CO2
C     EGRO4=light-limited PEP carboxylation rate
C     FDBK4=N,P feedback inhibition on C4 CO2 fixation
C     CH2O=total PEP carboxylation rate  
C     SURFX=unself-shaded leaf surface area
C     TAUS=fraction of direct radiation transmitted from layer above
C
      VL=AMIN1(VGRO4(K,NB,NZ,NY,NX),EGRO4)*FDBK4(K,NB,NZ,NY,NX)
      CH2O=CH2O+VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAUS(L+1,NY,NX)
C     IF(L.GT.NC-4.AND.NB.EQ.1.AND.M.EQ.1.AND.N.EQ.3)THEN
C     WRITE(*,6789)'STO',I,J,L,M,N,K,L,VL
C    2,PAR(N,M,L,NZ,NY,NX),RAPS(NY,NX)
C    2,TKC(NZ,NY,NX),CO2Q(NY,NX),ETGR4(K,NB,NZ,NY,NX)
C    3,CBXN4(K,NB,NZ,NY,NX),VGRO4(K,NB,NZ,NY,NX),EGRO4
C    3,FDBK4(K,NB,NZ,NY,NX),CH2O,VGRO4(K,NB,NZ,NY,NX),EGRO4
C    3,SURFX(N,L,K,NB,NZ,NY,NX)
C    3,VCGR4(K,NB,NZ,NY,NX),CO2I(NZ,NY,NX),CO2L(NZ,NY,NX),TFN1,TFN2
C    4,TFNE,WSDN,VCDN4
6789  FORMAT(A8,7I4,40E12.4)
C     ENDIF
      ENDIF
C
C     SHADED LEAVES
C
C     PARDIF=diffuse PAR flux from ‘hour1.f’
C
      IF(PARDIF(N,M,L,NZ,NY,NX).GT.0.0)THEN
C
C     LIGHT-LIMITED CARBOXYLATION RATES
C
C     QNTM=quantum efficiency
C     PARDIF=diffuse PAR flux from ‘hour1.f’
C     ETGR4=light saturated e- transport rate
C     ETLF4=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO4=light-limited PEP carboxylation rate
C
      PARX=QNTM*PARDIF(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGR4(K,NB,NZ,NY,NX)
      ETLF4=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGR4(K,NB,NZ,NY,NX)))/CURV2
      EGRO4=ETLF4*CBXN4(K,NB,NZ,NY,NX)
C
C     C4 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
C
C     VL=PEP carboxylation rate limited by light,CO2,N,P
C     VGRO4=PEP carboxylation rate limited by CO2
C     EGRO4=light-limited PEP carboxylation rate
C     CH2O=total PEP carboxylation rate  
C     FDBK4=N,P feedback inhibition on C4 CO2 fixation 
C     SURFX=unself-shaded leaf surface area from ‘hour1.f’
C     TAU0=fraction of diffuse radiation transmitted from layer above
C        from ‘hour1.f’
C
      VL=AMIN1(VGRO4(K,NB,NZ,NY,NX),EGRO4)*FDBK4(K,NB,NZ,NY,NX)
      CH2O=CH2O+VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAU0(L+1,NY,NX)
C     WRITE(*,6799)'STB',I,J,L,M,N,K
C    2,VL,PAR(N,M,L,NZ,NY,NX),RAPS(NY,NX)
C    2,TKC(NZ,NY,NX),CO2Q(NY,NX),ETGR4(K,NB,NZ,NY,NX)
C    3,CBXN4(K,NB,NZ,NY,NX),VGRO4(K,NB,NZ,NY,NX),EGRO4
C    3,FDBK4(K,NB,NZ,NY,NX),CH2O,VGRO4(K,NB,NZ,NY,NX),EGRO4
C    3,VCGR4(K,NB,NZ,NY,NX),CO2I(NZ,NY,NX),CO2L(NZ,NY,NX)
6799  FORMAT(A8,6I4,40E12.4)
      ENDIF
      ENDIF
2500  CONTINUE
2600  CONTINUE
      ENDIF
2700  CONTINUE
C
C     VARIABLES FOR C3 PHOTOSYNTHESIS DRIVEN BY C4
C
C     VCDN=surficial density of rubisco in bundle sheath
C     ETDN=surficial density of chlorophyll in bundle sheath
C     RUBP=fraction of leaf protein in rubisco
C        from PFT file
C     CHL=fraction of leaf protein in bundle sheath chlorophyll
C        from PFT file
C     WSDN=leaf protein surficial density
C
      VCDN=RUBP(NZ,NY,NX)*WSDN
      ETDN=CHL(NZ,NY,NX)*WSDN
C
C     CO2-LIMITED C3 CARBOXYLATION RATES
C
C     VCGRO=rubisco carboxylation rate unlimited by CO2
C     VCMX=specific rubisco carboxylation activity from PFT file
C     TFN1=temperature function for carboxylation
C     VCDN=surficial density of rubisco in bundle sheath
C     VOGRO=rubisco oxygenation rate
C     TFN2=temperature function for oxygenation
C     COMPL=C3 CO2 compensation point (uM)
C     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
C     XKCO2L,XKCO2O=Km for rubisco carboxylation without,with O2
C     XKO2L=Km for rubisco oxygenation
C     VGRO=rubisco carboxylation rate limited by CO2
C     CCBS=C4 nonstructural C concentration in bundle sheath (uM)
C
      VCGRO(K,NB,NZ,NY,NX)=VCMX(NZ,NY,NX)*TFN1*VCDN
      VOGRO=VOMX(NZ,NY,NX)*TFN2*VCDN
      COMPL(K,NB,NZ,NY,NX)=0.5*O2L(NZ,NY,NX)*VOGRO*XKCO2L(NZ,NY,NX)
     2/(VCGRO(K,NB,NZ,NY,NX)*XKO2L)
      VGRO(K,NB,NZ,NY,NX)=AMAX1(0.0,VCGRO(K,NB,NZ,NY,NX)
     2*(CCBS-COMPL(K,NB,NZ,NY,NX))/(CCBS+XKCO2O(NZ,NY,NX)))
C
C     C3 ELECTRON TRANSFER RATES
C
C     ETGRO=light-limited rubisco carboxylation rate
C     ETMX=specific chlorophyll activity from PFT file
C     TFNE=temperature function for e- transport 
C     ETDN=surficial density of chlorophyll in bundle sheath
C     CBXN=rubisco caboxylation efficiency
C     CO2L=intercellular CO2 concentrations (uM)
C     COMPL=C3 CO2 compensation point (uM)
C     ELEC3=e- requirement for CO2 fixation by rubisco
C
      ETGRO(K,NB,NZ,NY,NX)=ETMX(NZ,NY,NX)*TFNE*ETDN
      CBXN(K,NB,NZ,NY,NX)=AMAX1(0.0,(CCBS-COMPL(K,NB,NZ,NY,NX))
     2/(ELEC3*CCBS+10.5*COMPL(K,NB,NZ,NY,NX)))
C
C     C3 PHOTOSYNTHESIS
C
      ELSE
C
C     SURFICIAL DENSITY OF RUBISCO AND ITS CHLOROPHYLL
C
C     VCDN=surficial density of rubisco in mesophyll
C     ETDN=surficial density of chlorophyll in esophyll
C     RUBP=fraction of leaf protein in rubisco
C        from PFT file
C     CHL=fraction of leaf protein in mesophyll chlorophyll 
C        from PFT file
C     WSDN=leaf protein surficial density
C
      VCDN=RUBP(NZ,NY,NX)*WSDN
      ETDN=CHL(NZ,NY,NX)*WSDN
C
C     CO2-LIMITED C3 CARBOXYLATION RATES
C
C     VCGRO=rubisco carboxylation rate unlimited by CO2
C     VCMX=specific rubisco carboxylation activity from PFT file
C     TFN1=temperature function for carboxylation
C     VCDN=surficial density of rubisco in mesophyll
C     VOGRO=rubisco oxygenation rate
C     TFN2=temperature function for oxygenation
C     COMPL=C3 CO2 compensation point (uM)
C     CO2L,O2L=intercellular CO2,O2 concentrations (uM)
C     XKCO2L,XKCO2O=Km for rubisco carboxylation without,with O2
C     XKO2L=Km for rubisco oxygenation
C     VGRO=rubisco carboxylation rate limited by CO2
C
      VCGRO(K,NB,NZ,NY,NX)=VCMX(NZ,NY,NX)*TFN1*VCDN
      VOGRO=VOMX(NZ,NY,NX)*TFN2*VCDN
      COMPL(K,NB,NZ,NY,NX)=0.5*O2L(NZ,NY,NX)*VOGRO*XKCO2L(NZ,NY,NX)
     2/(VCGRO(K,NB,NZ,NY,NX)*XKO2L)
      VGRO(K,NB,NZ,NY,NX)=AMAX1(0.0,VCGRO(K,NB,NZ,NY,NX)
     2*(CO2L(NZ,NY,NX)-COMPL(K,NB,NZ,NY,NX))
     5/(CO2L(NZ,NY,NX)+XKCO2O(NZ,NY,NX)))
C
C     C3 ELECTRON TRANSFER RATES
C
C     ETGRO=light-limited rubisco carboxylation rate
C     ETMX=specific chlorophyll activity from PFT file
C     TFNE=temperature function for e- transport 
C     ETDN=surficial density of chlorophyll in mesophyll
C     CBXN=rubisco caboxylation efficiency
C     CO2L=intercellular CO2 concentrations (uM)
C     COMPL=C3 CO2 compensation point (uM)
C     ELEC3=e- requirement for CO2 fixation by rubisco
C
      ETGRO(K,NB,NZ,NY,NX)=ETMX(NZ,NY,NX)*TFNE*ETDN
      CBXN(K,NB,NZ,NY,NX)=AMAX1(0.0,(CO2L(NZ,NY,NX)
     2-COMPL(K,NB,NZ,NY,NX))/(ELEC3*CO2L(NZ,NY,NX)
     3+10.5*COMPL(K,NB,NZ,NY,NX)))
C
C     FOR EACH CANOPY LAYER
C
C     ARLFL=node leaf area in canopy layer
C     SURFX=unself-shaded leaf surface area from ‘hour1.f’
C
      DO 3700 L=JC,1,-1
      IF(ARLFL(L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C
C     FOR EACH INCLINATION AND AZIMUTH CLASS
C
      DO 3600 N=1,4
      DO 3500 M=1,4
      IF(SURFX(N,L,K,NB,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
C     PAR=direct PAR flux from ‘hour1.f’
C
C     SUNLIT LEAVES
C
      IF(PAR(N,M,L,NZ,NY,NX).GT.0.0)THEN
C
C     LIGHT-LIMITED CARBOXYLATION RATES
C
C     QNTM=quantum efficiency 
C     PAR=direct PAR flux from ‘hour1.f’
C     ETGRO=light saturated e- transport rate
C     ETLF=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO=light-limited rubisco carboxylation rate
C
      PARX=QNTM*PAR(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
      ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
      EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
C
C     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
C
C     VL=rubisco carboxylation rate limited by light,CO2,N,P
C     VGRO=rubisco carboxylation rate limited by CO2
C     EGRO=light-limited rubisco carboxylation rate
C     FDBK=N,P feedback inhibition on C3 CO2 fixation
C     CH2O=total rubisco carboxylation rate  
C     SURFX=unself-shaded leaf surface area
C     TAUS=fraction of direct radiation transmitted from layer above
C        from ‘hour1.f’
C
      VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*FDBK(NB,NZ,NY,NX)
      CH2O=CH2O+VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAUS(L+1,NY,NX)
C     IF(NB.EQ.1.AND.M.EQ.1.AND.N.EQ.1.AND.K.EQ.KLEAF(NB,NZ,NY,NX)-1
C    2.AND.J.EQ.14)THEN
C     WRITE(*,6798)'STD',I,J,NFZ,L,M,N,K,NB
C    2,VL,CH2O,PAR(N,M,L,NZ,NY,NX),RAPS(NY,NX)
C    2,TKC(NZ,NY,NX),CO2Q(NY,NX),ETGRO(K,NB,NZ,NY,NX)
C    3,CBXN(K,NB,NZ,NY,NX),VGRO(K,NB,NZ,NY,NX),EGRO
C    3,FDBK(NB,NZ,NY,NX),TFN1,TFN2,TFNE,WSDN
C    3,VCGRO(K,NB,NZ,NY,NX),VCDN,CO2I(NZ,NY,NX),CO2L(NZ,NY,NX)
6798  FORMAT(A8,8I4,40E12.4)
C     ENDIF
      ENDIF
C
C     SHADED LEAVES
C
      IF(PARDIF(N,M,L,NZ,NY,NX).GT.0.0)THEN
C
C     LIGHT-LIMITED CARBOXYLATION RATES
C
C     QNTM=quantum efficiency
C     PARDIF=diffuse PAR flux from ‘hour1.f’
C     ETGR=light saturated e- transport rate
C     ETLF=light-limited e- transport rate
C     CURV=shape parameter for e- transport response to PAR
C     EGRO=light-limited rubisco carboxylation rate
C
      PARX=QNTM*PARDIF(N,M,L,NZ,NY,NX)
      PARJ=PARX+ETGRO(K,NB,NZ,NY,NX)
      ETLF=(PARJ-SQRT(PARJ**2-CURV4*PARX*ETGRO(K,NB,NZ,NY,NX)))/CURV2
      EGRO=ETLF*CBXN(K,NB,NZ,NY,NX)
C
C     C3 CARBOXYLATION RATE AND ACCUMULATED PRODUCT
C
C     VL=rubisco carboxylation rate limited by light,CO2,N,P
C     VGRO=rubisco carboxylation rate limited by CO2
C     EGRO=light-limited rubisco carboxylation rate
C     CH2O=total rubisco carboxylation rate  
C     FDBK=N,P feedback inhibition on C3 CO2 fixation 
C     SURFX=unself-shaded leaf surface area
C     TAU0=fraction of diffuse radiation transmitted from layer above
C        from ‘hour1.f’
C
      VL=AMIN1(VGRO(K,NB,NZ,NY,NX),EGRO)*FDBK(NB,NZ,NY,NX)
      CH2O=CH2O+VL*SURFX(N,L,K,NB,NZ,NY,NX)*TAU0(L+1,NY,NX)
      ENDIF
      ENDIF
3500  CONTINUE
3600  CONTINUE
      ENDIF
3700  CONTINUE
      ENDIF
      ELSE
      VCGR4(K,NB,NZ,NY,NX)=0.0
      VCGRO(K,NB,NZ,NY,NX)=0.0
      ENDIF
2800  CONTINUE
      ENDIF
      ELSE
      FDBK(NB,NZ,NY,NX)=0.0
      FDBKX(NB,NZ,NY,NX)=1.0
      DO 2805 K=1,25
      VCGR4(K,NB,NZ,NY,NX)=0.0
      VCGRO(K,NB,NZ,NY,NX)=0.0
2805  CONTINUE
      ENDIF
2900  CONTINUE
C
C     MINIMUM CANOPY STOMATAL RESISTANCE FROM CO2 CONCENTRATION
C     DIFFERENCE DIVIDED BY TOTAL CO2 FIXATION
C
C     RSX,RSMN=minimum canopy stomatal resistance to CO2,H2O (h m-1)
C     CH2O=total PEP(C4) or rubisco(C3) carboxylation rate  
C     FRADP=fraction of radiation received by each PFT canopy
C        received by each PFT canopy from ‘hour1.f’
C     DCO2=difference between atmosph and intercellular CO2
C        concentration (umol m-3)
C     AREA=area of grid cell 
C     RSMY=minimum stomatal resistance for CO2 uptake (h m-1)
C
      IF(CH2O.GT.ZEROP(NZ,NY,NX))THEN
      RSX=FRADP(NZ,NY,NX)*DCO2(NZ,NY,NX)
     2*AREA(3,NU(NY,NX),NY,NX)/(CH2O*3600.0)
      ELSE
      RSX=RSMH(NZ,NY,NX)*1.56
      ENDIF
      RSMN(NZ,NY,NX)=AMIN1(RSMH(NZ,NY,NX),AMAX1(RSMY,RSX*0.641))
      ELSE
      RSMN(NZ,NY,NX)=RSMH(NZ,NY,NX)
      ENDIF
C     IF(ICTYP(NZ,NY,NX).EQ.3)THEN
C     WRITE(*,3010)'CH2O',I,J,NFZ,CH2O
C     ELSEIF(ICTYP(NZ,NY,NX).EQ.4)THEN
C     WRITE(*,3010)'CH2O',I,J,NFZ,CH2O
C     ENDIF
3010  FORMAT(A8,3I4,1E12.4)
      RETURN
      END

