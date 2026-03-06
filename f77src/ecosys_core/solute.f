      SUBROUTINE solute(I,J,NFZ,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE CALCULATES ALL SOLUTE TRANSFORMATIONS
C     FROM THERMODYNAMIC EQUILIBRIA
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
      include "blk15a.h"
      include "blk15b.h"
      include "blk18a.h"
      include "blk18b.h"
      include "blk19a.h"
      include "blk19b.h"
      include "blk19c.h"
      include "blk19d.h"
      include "blk21a.h"
      include "blk21b.h"
C
C     EQUILIBRIUM CONSTANTS
C
      DIMENSION RNHUI(0:2)
C
C     DISSOLUTION AND DISSOCIATION CONSTANTS 
C        (soluble concentrations in mol m-3)
C
C     DPH20=water,SPALO=AlOH3s,SPFEO=FeOH3s,SPCAC=CaCO3s,SPCAS=CaSO4s
C     SPALP=AlPO4s,SPFEP=FePO4s,SPCAM=Ca(H2PO4)2s,SPCAD=CAHPO4s
C     SPCAH=hydroxyapatite,SXOX2=R-OH2,SXOX1=R-OH,SXH2P=R-H2PO4
C     SXH1P=R-HPO4,DPCO2=CO2,DPHCO=HCO3,DPN4=NH4,DPAL1=ALOH
C     DPAL2=ALOH2,DPAL3=AlOH3,DPAL4=ALOH4,DPALS=AlSO4,DPFE1=FeOH
C     DPFE2=FeOH2,DPFE3=FeOH3,DPFE4=FeOH4,DPFES=FeSO4
C     DPCAO=CaOH,DPCAC=CaCO3,DPCAH=CaHCO3,DPCAS=CaSO4
C     DPMGO=MgOH,DPMGC=MgCO3,DPMGH=MgHCO3,DPMGS=MgSO4
C     DPNAC=NaCO3,DPNAS=NaSO4,DPKAS=KSO4,DPH1P=HPO4
C     DPH2P=H2PO4,DPH3P=H3PO4,DPF1P=FeHPO4,DPF2P=FeH2PO4
C     DPC0P=CaPO4,DPC1P=CaHPO4,DPC2P=CaH2PO4
C     DPM1P=MgHPO4,DPCOH=R-COO,COOH=carboxyl sites (mol Mg C-1)
C     SPALSI=Al silicate,SPFESI=Fe silicate,SPCASI=Ca silicate
C     SPMGSI=Mg silicate,SPNASI=Na silicate,SPKASI=K silicate
C
      PARAMETER (DPH2O=1.0E-08,SPALO=1.9E-21,SPFEO=6.3E-26
     2,SPCAC=3.3E-03,SPCAS=1.4E+01,SPALP=9.8E-15,SPFEP=2.5E-19
     3,SPCAM=7.0E+07,SPCAD=1.3E-01,SPCAH=4.0E-31,SXOH2=4.5E-01
     4,SXOH1=8.1E-04,SXH2P=0.5E+06,SXH1P=0.5E+04
     5,DPCO2=4.2E-04,DPHCO=5.6E-08,DPN4=5.5E-07
     6,DPAL1=4.6E-07,DPAL2=7.3E-07,DPAL3=1.8E-05
     7,DPAL4=1.2E-05,DPALS=0.16,DPFE1=2.7E-08,DPFE2=4.5E-07
     8,DPFE3=2.5E-05,DPFE4=1.2E-05,DPFES=7.1E-02,DPCAO=12.5
     9,DPCAC=4.2E-02,DPCAH=13.5,DPCAS=1.2,DPMGO=0.7,DPMGC=0.3
     1,DPMGH=67.0,DPMGS=2.1,DPNAC=0.45,DPNAS=3.3E+02,DPKAS=5.0E+01
     2,DPH1P=4.8E-10,DPH2P=6.2E-05,DPH3P=7.5,DPF1P=4.5E-02
     3,DPF2P=3.7E-03,DPC0P=3.5E-03,DPC1P=1.82,DPC2P=40.0
     4,DPM1P=1.23,DPCOH=1.0E-02,COOH=2.5E+02
     5,SPALSI=1.0E+05,SPFESI=1.0E+05,SPCASI=1.0E+03,SPMGSI=1.0E+03
     6,SPNASI=1.0E+02,SPKASI=1.0E+02)
C
C     DISSOCIATION CONSTANTS CALCULATED FROM OTHER CONSTANTS
C
C     DPCO3=CO3
C
      PARAMETER (DPCO3=DPCO2*DPHCO 
     8,SHCAH1=SPCAH/(DPH2O*DPH1P**3),SHCAH2=SHCAH1/DPH2P**3)
C
C     MRXN=number of cycles for solving reaction equilibria
C     TPDH,TADAH,TADCH,TSLH,TRWH=hourly rate constants for 
C        precipitation,adsorption,solute,silicate weathwering
C        equilibria (h-1)
C     TPDXH,TADAXH,TADCXH,TSLXH,TRWXH=rate constants for 
C        precipitation,adsorption,solute,silicate weathwering
C        equilibria (MRXN-1)
C     FION*=limits on reaction rates imposed by substrates, products
C        (t-1)
C     SPNH4H,SPNH3H,SPNO3H =hourly specific rate constants for 
C        NH4,NH3,NO3 fertilizer dissolution (h-1)
C     SPNHUH=specific hourly urea hydrolysis rate from microbial
C         activity (mol N g C-1 h-1)
C
      PARAMETER (MRXN=60
     2,TPDH=2.5E-02,TPDXH=TPDH/MRXN
     3,TPZH=2.5E-02,TPZXH=TPZH/MRXN
     4,TADAH=5.0E-03,TADAXH=TADAH/MRXN
     5,TADCH=5.0E-03,TADCXH=TADCH/MRXN
     6,TSLH=2.5E-02,TSLXH=TSLH/MRXN
     7,TSZH=2.5E-01,TSZXH=TSZH/MRXN 
     4,TRWH=1.0E-05,TRWXH=TRWH/MRXN
     5,FIONH=1.00,FIONX=0.2,FIONZ=0.2,OHKI=1.0E-06)
      PARAMETER (SPNH4H=1.0E-00,SPNH3H=1.0E-00
     2,SPNO3H=1.0E-00)
      PARAMETER (SPNHUH=0.3E-01)
C
C     DUKM=minimum Km for urea hydrolysis (mol N Mg-1) 
C     DUKI=Ki for microbial activity effects on urea
C        hydrolysis as in ‘nitro.f’ (SOIL SCI 136:56) (g C m-3 h-1)
C     RNHUI=rate constant for decline in urea hydrolysis inhibition 
C        :0=urine,1=normal release,2=slow release (h-1)
C 
      PARAMETER (DUKM=0.05,DUKI=50.0)
      PARAMETER (ZEROC=1.0E-32)
      DATA RNHUI/5.0E-02,1.0E-02,0.5E-02/
C
C     INITIALIZE RATE CONSTANTS FOR SOLUTE TRANSFORMATIONS
C
C     SPNH4,SPNH3,SPNO3,SPNHU=specific rate constants for 
C        NH4,NH3,NO3,urea fertilizer dissolution (t-1)
C     XNFH=time step for solute transformations from ‘wthr.f’ (h t-1)
C
      NPI=INT(NPH/2)
      FIONN=FIONZ 
      FION0=FIONZ 
      SPNH4=SPNH4H*XNFH
      SPNH3=SPNH3H*XNFH
      SPNO3=SPNO3H*XNFH
      SPNHU=SPNHUH*XNFH
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
C
C     RNHUX=rate constant for decline in urea hydrolysis inhibition 
C        (t-1)
C     RNHUI=rate constant for decline in urea hydrolysis inhibition 
C        :0=urine,1=normal release,2=slow release (h-1)
C     IUTYP=urea formulation (0=fast release (urine),
C        1=normal release,2=slow release)
C
      RNHUX=RNHUI(IUTYP(NY,NX))*XNFH
      DO 9985 L=NU(NY,NX),NL(NY,NX)
      IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX)
     2.AND.VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
C
C     WATER VOLUME IN NON-BAND AND BAND SOIL ZONES
C
C     VOLX,VOLWM=soil layer,water volume (m3)
C     VLNH4,VLNHB=fractions of soil volume in NH4 non-band,band
C     VLNO3,VLNOB=fractions of soil volume in NO3 non-band,band
C     VLPO4,VLPOB=fractions of soil volume in H2PO4 non-band,band
C     BKVL=soil mass (Mg)
C
      VOLWNH=VOLW(L,NY,NX)*VLNH4(L,NY,NX)
      VOLWNB=VOLW(L,NY,NX)*VLNHB(L,NY,NX)
      VOLWNO=VOLW(L,NY,NX)*VLNO3(L,NY,NX)
      VOLWNZ=VOLW(L,NY,NX)*VLNOB(L,NY,NX)
      VOLWPO=VOLW(L,NY,NX)*VLPO4(L,NY,NX)
      VOLWPB=VOLW(L,NY,NX)*VLPOB(L,NY,NX)
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      BKVLX=BKVL(L,NY,NX)
      BKVLNH=BKVL(L,NY,NX)*VLNH4(L,NY,NX)
      BKVLNB=BKVL(L,NY,NX)*VLNHB(L,NY,NX)
      BKVLNO=BKVL(L,NY,NX)*VLNO3(L,NY,NX)
      BKVLNZ=BKVL(L,NY,NX)*VLNOB(L,NY,NX)
      BKVLPO=BKVL(L,NY,NX)*VLPO4(L,NY,NX)
      BKVLPB=BKVL(L,NY,NX)*VLPOB(L,NY,NX)
      ELSE
      BKVLX=VOLA(L,NY,NX)
      BKVLNH=VOLWNH
      BKVLNB=VOLWNB
      BKVLNO=VOLWNO
      BKVLNZ=VOLWNZ
      BKVLPO=VOLWPO
      BKVLPB=VOLWPB
      ENDIF
C
C     UREA HYDROLYSIS IN BAND AND NON-BAND SOIL ZONES
C
C     VOLQ=biologically active soil water volume from ‘nitro.f’ (m3)
C     COQCK=concentration of active biomass activity (g C m-3 t-1)
C     TOQCK=total active biomass respiration activity from ‘nitro.f’
C        (g C t-1)
C     DUKD,DUKM=effective,minimum Km for urea hydrolysis (mol N Mg-1) 
C     DUKI=Ki for microbial activity effects on urea
C        hydrolysis as in ‘nitro.f’ (SOIL SCI 136:56) (g C m-3 t-1)
C
      IF(VOLQ(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      COQCK=AMIN1(0.1E+06,TOQCK(L,NY,NX)/(VOLQ(L,NY,NX)*XNFH))
      ELSE
      COQCK=0.1E+06
      ENDIF
      DUKD=DUKM*(1.0+COQCK/DUKI)
C
C     UREA HYDROLYSIS INHIBITION
C
C     ZNHU0,ZNHUI=initial,current urea hydrolysis inhibition activity 
C        from ‘hour1.f’
C     RNHUX=rate constant for decline in urea hydrolysis inhibition 
C        (t-1)
C
      IF(ZNHU0(L,NY,NX).GT.ZEROS(NY,NX)
     2.AND.ZNHUI(L,NY,NX).GT.ZEROS(NY,NX))THEN
      ZNHUI(L,NY,NX)=ZNHUI(L,NY,NX)
     2-RNHUX*ZNHUI(L,NY,NX) 
     3*AMAX1(RNHUX,1.0-ZNHUI(L,NY,NX)/ZNHU0(L,NY,NX))
      ELSE
      ZNHUI(L,NY,NX)=0.0
      ENDIF
C
C     UREA CONCENTRATION AND HYDROLYSIS IN NON-BAND
C
C     ZNHUFA=urea fertilizer in non-band (mol N)
C     BKVL,VOLW=soil mass,water content (Mg,m3)
C     CNHUA=concentration of urea fertilizer in non-band (mol N Mg-1)
C     DFNSA=effect of microbial activity on urea hydrolysis in
C        non-band
C     DUKD=effective Km for urea hydrolysis (mol N Mg-1) 
C     RSNUA=rate of urea hydrolysis in non-band (mol N t-1)
C     SPNHU=specific urea hydrolysis rate from microbial activity 
C        (mol N g C-1)
C     TOQCK=total active biomass respiration activity from ‘nitro.f’
C        (g C t-1)
C     TFNQ=temperature effect on microbial activity from ‘nitro.f’
C     ZNHUI=current urea hydrolysis inhibition activity
C
      IF(ZNHUFA(L,NY,NX).GT.ZEROS(NY,NX)
     2.AND.BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNHUA=ZNHUFA(L,NY,NX)/BKVL(L,NY,NX)
      ELSEIF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      CNHUA=ZNHUFA(L,NY,NX)/VOLW(L,NY,NX)
      ELSE
      CNHUA=0.0
      ENDIF
      DFNSA=CNHUA/(CNHUA+DUKD)
      RSNUA=AMIN1(ZNHUFA(L,NY,NX)
     2,SPNHU*TOQCK(L,NY,NX)*DFNSA*TFNQ(L,NY,NX)*(1.0-ZNHUI(L,NY,NX)))
C
C     UREA CONCENTRATION AND HYDROLYSIS IN BAND
C
C     ZNHUFB=urea fertilizer in band (mol N)
C     BKVL,VOLW=soil mass,water content (Mg,m3)
C     CNHUB=concentration of urea fertilizer in band (mol N Mg-1)
C     DFNSB=effect of microbial activity on urea hydrolysis 
C        in band
C     DUKD=effective Km for urea hydrolysis (mol N Mg-1) 
C     RSNUB=rate of urea hydrolysis in non-band (mol N t-1)
C     SPNHU=specific urea hydrolysis rate from microbial activity 
C        (mol N g C-1)
C     TOQCK=total active biomass respiration activity from ‘nitro.f’
C        (g C t-1)
C     TFNQ=temperature effect on microbial activity from ‘nitro.f’
C
      IF(ZNHUFB(L,NY,NX).GT.ZEROS(NY,NX)
     2.AND.BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CNHUB=ZNHUFB(L,NY,NX)/BKVL(L,NY,NX)
      ELSEIF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      CNHUB=ZNHUFB(L,NY,NX)/VOLW(L,NY,NX)
      ELSE
      CNHUB=0.0
      ENDIF 
      DFNSB=CNHUB/(CNHUB+DUKD)
      RSNUB=AMIN1(ZNHUFB(L,NY,NX)
     2,SPNHU*TOQCK(L,NY,NX)*DFNSB*TFNQ(L,NY,NX)*(1.0-ZNHUI(L,NY,NX)))
C     IF(ZNHUFA(L,NY,NX).GT.ZEROS(NY,NX))THEN
C     WRITE(*,8888)'UREA',I,J,NFZ,L,IUTYP(NY,NX)
C    2,ZNHUI(L,NY,NX),ZNHUFA(L,NY,NX),ZNHUFB(L,NY,NX),RSNUA,RSNUB
C    2,DFNSA,DFNSB,TFNQ(L,NY,NX),CNHUA,DUKD,DUKM,DUKI,TOQCK(L,NY,NX)
C    3,BKVL(L,NY,NX),SPNHU,ZNHU0(L,NY,NX)
C    4,RNHUX,VLNH4(L,NY,NX),VLNHB(L,NY,NX)
C    5,THETWL
8888  FORMAT(A8,5I4,40E12.4)
C     ENDIF
C
C     NH4, NH3, UREA, NO3 DISSOLUTION IN BAND AND NON-BAND
C     SOIL ZONES FROM FIRST-ORDER FUNCTIONS OF REMAINING
C     FERTILIZER (NOTE: SUPERPHOSPHATE AND ROCK PHOSPHATE
C     ARE REPRESENTED AS MONOCALCIUM PHOSPHATE AND HYDROXYAPATITE
C     MODELLED IN PHOSPHORUS REACTIONS BELOW)
C
C     RSN4AA,RSN4BA=rate of broadcast NH4 fertilizer dissolution 
C        in non-band,band (mol N t-1)
C     RSN3AA,RSN3BA=rate of broadcast NH3 fertilizer dissolution 
C        in non-band,band (mol N t-1)
C     RSNUAA,RSNUBA=rate of broadcast urea fertr dissolution 
C        in non-band,band (mol N t-1)
C     RSNOAA,RSNOBA=rate of broadcast NO3 fertilizer dissolution 
C        in non-band,band (mol N t-1)
C     RSN4BB=rate of banded NH4 fertilizer dissolution in band 
C        (mol N t-1)
C     RSN3BB=rate of banded NH3 fertilizer dissolution in band
C        (mol N t-1)
C     RSNUBB=rate of banded urea fertilizer dissolution in band
C        (mol N t-1)
C     RSNOBB=rate of banded NO3 fertilizer dissolution in band
C        (mol N t-1)
C     SPNH4,SPNH3,SPNO3 =specific rate constants for 
C        NH4,NH3,NO3 fertilizer dissolution (t-1)
C     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=broadcast NH4,NH3,urea,NO3
C        fertilizer (mol N)
C     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
C        (mol N)
C     RSNUA,RSNUB=rate of urea hydrolysis in non-band,band (mol N t-1)
C     VLNH4,VLNHB=fractions of soil volume in NH4 non-band,band
C     VLNO3,VLNOB=fractions of soil volume in N03 non-band,band
C     VLPO4,VLPOB=fractions of soil volume in H2PO4 non-band,band
C     THETW=soil water concentration (m3 m-3)
C
      THETWL=THETW(L,NY,NX)
      RSN4AA=SPNH4*ZNH4FA(L,NY,NX)*VLNH4(L,NY,NX)*THETWL
      RSN3AA=SPNH3*ZNH3FA(L,NY,NX)*VLNH4(L,NY,NX)
      RSNUAA=RSNUA*VLNH4(L,NY,NX) 
      RSNOAA=SPNO3*ZNO3FA(L,NY,NX)*VLNO3(L,NY,NX)*THETWL
      RSN4BA=SPNH4*ZNH4FA(L,NY,NX)*VLNHB(L,NY,NX)*THETWL
      RSN3BA=SPNH3*ZNH3FA(L,NY,NX)*VLNHB(L,NY,NX)
      RSNUBA=RSNUA*VLNHB(L,NY,NX) 
      RSNOBA=SPNO3*ZNO3FA(L,NY,NX)*VLNOB(L,NY,NX)*THETWL
      RSN4BB=SPNH4*ZNH4FB(L,NY,NX)*THETWL
      RSN3BB=SPNH3*ZNH3FB(L,NY,NX) 
      RSNUBB=RSNUB*VLNHB(L,NY,NX) 
      RSNOBB=SPNO3*ZNO3FB(L,NY,NX)*THETWL
C
C     SOLUBLE AND EXCHANGEABLE NH4 CONCENTRATIONS
C     IN NON-BAND AND BAND SOIL ZONES
C
C     VOLWNH,VOLWNB=water volume in NH4 non-band,band (m3)
C     RN4X,RN3X=NH4,NH3 input from uptake,mineralization,dissolution
C        in non-band (g N t-1) 
C     RNBX,R3BX=NH4,NH3 input from uptake,mineralization,dissolution
C        in band (g N t-1) 
C     TUPNH4,TUPNH3=soil-root exchange of NH4,NH3 in non-band 
C        from ‘uptake.f’ (g N t-1) 
C     TUPNHB,TUPN3B=soil-root exchange of NH4,NH3 in band 
C        from ‘uptake.f’ (g N t-1) 
C     RSN4AA,RSN4BA=rate of broadcast NH4 fertilizer dissolution 
C        in non-band,band (mol N t-1)
C     RSNUAA,RSNUBA=rate of broadcast urea fertr dissolution 
C        in non-band,band (mol N t-1)
C     RSN4BB=rate of banded NH4 fertilizer dissolution in band 
C        (mol N t-1)
C     RSNUBB=rate of banded urea fertilizer dissolution in band
C        (mol N t-1)
C     XNH4S,XNH4B=net change in NH4 in band,non-band from ‘nitro.f’ 
C        (g N t-1) 
C     CN41,CN31,CN4B,CN3B=total NH4,NH3 concentration in non-band,band
C        (mol N m-3)
C     XN41,XNB1=adsorbed NH4 concentration in non-band,band
C        (mol N Mg-1)
C
      IF(ISALTG.NE.0)THEN
      XMRXN=MRXN
      ELSE
      XMRXN=1.0
      ENDIF
      IF(VOLWNH.GT.ZEROS2(NY,NX))THEN
      VOLWNX=14.0*VOLWNH 
      RN4X=(XNH4S(L,NY,NX)-TUPNH4(L,NY,NX)+14.0*RSN4AA)
     2/(XMRXN*VOLWNX)
      RN3X=(-TUPN3S(L,NY,NX)+14.0*RSNUAA)
     2/(XMRXN*VOLWNX)
      CN41=AMAX1(ZEROC,ZNH4S(L,NY,NX)/VOLWNX+RN4X) 
      CN31=AMAX1(ZEROC,ZNH3S(L,NY,NX)/VOLWNX+RN3X) 
      XN41=AMAX1(ZEROC,XN4(L,NY,NX)/BKVLNH)
      BKVLWH=BKVLNH/VOLWNH 
C     IF((I/30)*30.EQ.I.AND.J.EQ.24)THEN
C     WRITE(*,25)'RN4X',I,J,NFZ,NX,NY,L
C    2,RN4X,XNH4S(L,NY,NX),TUPNH4(L,NY,NX),14.0*RSN4AA
C    2,RN3X,TUPN3S(L,NY,NX),14.0*RSNUAA,ZNHUI(L,NY,NX) 
C    2,XMRXN,VOLWNX,CN41,CN31 
25    FORMAT(A8,6I4,20E12.4)
C     ENDIF
      ELSE
      RN4X=0.0
      RN3X=0.0
      CN41=0.0
      CN31=0.0
      XN41=0.0
      BKVLWH=0.0 
      ENDIF
      IF(VOLWNB.GT.ZEROS2(NY,NX))THEN
      VOLWNX=14.0*VOLWNB 
      RNBX=(XNH4B(L,NY,NX)-TUPNHB(L,NY,NX)+14.0*(RSN4BA+RSN4BB))
     2/(XMRXN*VOLWNX)
      R3BX=(-TUPN3B(L,NY,NX)+14.0*(RSNUBA+RSNUBB))
     2/(XMRXN*VOLWNX)
      CN4B=AMAX1(ZEROC,ZNH4B(L,NY,NX)/VOLWNX+RNBX)
      CN3B=AMAX1(ZEROC,ZNH3B(L,NY,NX)/VOLWNX+R3BX)
      XNB1=AMAX1(ZEROC,XNB(L,NY,NX)/BKVLNB)
      BKVLWB=BKVLNB/VOLWNB 
      ELSE
      RNBX=0.0
      R3BX=0.0
      CN4B=0.0
      CN3B=0.0
      XNB1=0.0
      BKVLWB=0.0 
      ENDIF
C     IF(IYRC.EQ.2012.AND.I.EQ.151.AND.NX.EQ.1)THEN
C     WRITE(*,4141)'RN4X',I,J,NX,NY,L,RN4X,RN3X,RNBX,R3BX
C    2,CN41,CN31,CN4B,CN3B,ZNH4S(L,NY,NX),ZNH3S(L,NY,NX)
C    3,RSN4AA,TUPN3S(L,NY,NX),RSNUAA,TUPNHB(L,NY,NX)
C    4,XNH4B(L,NY,NX),RSN4BA,RSN4BB,TUPN3B(L,NY,NX)
C    5,RSNUBA,RSNUBB,ZNH4S(L,NY,NX),ZNH3S(L,NY,NX)
C    6,VOLWNX,BKVLNH 
4141  FORMAT(A8,5I4,30E12.4)
C     ENDIF
C
C     SOLUBLE, EXCHANGEABLE AND PRECIPITATED PO4 CONCENTRATIONS IN
C     NON-BAND AND BAND SOIL ZONES
C
C     VOLWPO,VOLWPB=water volume in H2PO4 non-band,band (m3)
C     RH1PX,RH2PX=HPO4,H2PO4 inputs from mineralization,uptake 
C        in non-band (mol P t-1)
C     RH1BX,RH2BX=HPO4,H2PO4 inputs from mineralization,uptake 
C        in band (mol P t-1)
C     XH1PS,XH1BS=net change in HPO4 in band,non-band from ‘nitro.f’ 
C        (g P t-1)
C     TUPH1P,TUPH2P=soil-root exchange of HPO4,H2PO4 in non-band 
C        from ‘uptake.f’ (g P t-1)
C     TUPH1B,TUPH2B=soil-root exchange of HPO4,H2PO4 in band 
C        from ‘uptake.f’ (g P t-1)
C     CH1P1,CH2P1=HPO4,H2PO4 concentrations in non-band (mol P m-3)
C     CH1PB,CH2PB=HPO4,H2PO4 concentrations in band (mol P m-3)
C     XOH01,XOH11,XOH21=concentration of adsorption sites R-,R-OH,
C        R-OH2 in non-band (mol Mg-1)
C     XOH1B,XH11B,XH21B= concentration of adsorption sites R-,R-OH,
C        R-OH2 in band (mol Mg-1)
C     XH1P1,XH2P1=concentration of adsorbed HPO4,H2PO4 in non-band
C        (mol Mg-1)
C     XH11B,X2P1B=concentration of adsorbed HPO4,H2PO4 in band
C        (mol Mg-1)
C     PALPO1,PFEPO1=precipitated AlPO4,FEPO4 in non-band (mol)
C     PALPOB,PFEPOB=precipitated AlPO4,FEPO4 in band (mol)
C     PCAPM1,PCAPD1,PCAPH1=precipitated CaH2PO4,CaHPO4,apatite 
C        in non-band (mol)
C     PCAPMB,PCAPDB,PCAPHB=precipitated CaH2PO4,CaHPO4,apatite
C        in band (mol)
C
      IF(VOLWPO.GT.ZEROS2(NY,NX))THEN
      VOLWPX=31.0*VOLWPO
      RH1PX=(XH1PS(L,NY,NX)-TUPH1P(L,NY,NX))/(XMRXN*VOLWPX) 
      RH2PX=(XH2PS(L,NY,NX)-TUPH2P(L,NY,NX))/(XMRXN*VOLWPX)
      CH1P1=AMAX1(ZEROC,H1PO4(L,NY,NX)/VOLWPX+RH1PX)
      CH2P1=AMAX1(ZEROC,H2PO4(L,NY,NX)/VOLWPX+RH2PX)
      XOH01=AMAX1(ZEROC,XOH0(L,NY,NX))/BKVLPO
      XOH11=AMAX1(ZEROC,XOH1(L,NY,NX))/BKVLPO
      XOH21=AMAX1(ZEROC,XOH2(L,NY,NX))/BKVLPO
      XH1P1=AMAX1(ZEROC,XH1P(L,NY,NX))/BKVLPO
      XH2P1=AMAX1(ZEROC,XH2P(L,NY,NX))/BKVLPO
      PALPO1=AMAX1(0.0,PALPO(L,NY,NX))/VOLWPO
      PFEPO1=AMAX1(0.0,PFEPO(L,NY,NX))/VOLWPO 
      PCAPM1=AMAX1(0.0,PCAPM(L,NY,NX))/VOLWPO 
      PCAPD1=AMAX1(0.0,PCAPD(L,NY,NX))/VOLWPO 
      PCAPH1=AMAX1(0.0,PCAPH(L,NY,NX))/VOLWPO 
C     IF(I.GT.140.AND.L.LE.3)THEN
C     WRITE(*,8642)'CH2P1',I,J,L,CH2P1,H2PO4(L,NY,NX)
C    2,VOLWPX,RH2PX,XH2PS(L,NY,NX),TUPH2P(L,NY,NX)
8642  FORMAT(A8,3I4,20E12.4)
C     ENDIF
      ELSE
      RH1PX=0.0
      RH2PX=0.0
      CH1P1=0.0
      CH2P1=0.0
      XOH01=0.0
      XOH11=0.0
      XOH21=0.0
      XH1P1=0.0
      XH2P1=0.0
      PALPO1=0.0
      PFEPO1=0.0
      PCAPM1=0.0
      PCAPD1=0.0
      PCAPH1=0.0
      ENDIF
      IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
      VOLWBX=31.0*VOLWPB
      RH1BX=(XH1BS(L,NY,NX)-TUPH1B(L,NY,NX))/(XMRXN*VOLWBX) 
      RH2BX=(XH2BS(L,NY,NX)-TUPH2B(L,NY,NX))/(XMRXN*VOLWBX) 
      CH1PB=AMAX1(ZEROC,H1POB(L,NY,NX)/VOLWBX+RH1BX)
      CH2PB=AMAX1(ZEROC,H2POB(L,NY,NX)/VOLWBX+RH2BX)
      XH01B=AMAX1(ZEROC,XOH0B(L,NY,NX))/BKVLPB
      XH11B=AMAX1(ZEROC,XOH1B(L,NY,NX))/BKVLPB
      XH21B=AMAX1(ZEROC,XOH2B(L,NY,NX))/BKVLPB
      X1P1B=AMAX1(ZEROC,XH1PB(L,NY,NX))/BKVLPB
      X2P1B=AMAX1(ZEROC,XH2PB(L,NY,NX))/BKVLPB
      PALPOB=AMAX1(0.0,PALPB(L,NY,NX))/VOLWPB 
      PFEPOB=AMAX1(0.0,PFEPB(L,NY,NX))/VOLWPB 
      PCAPMB=AMAX1(0.0,PCPMB(L,NY,NX))/VOLWPB 
      PCAPDB=AMAX1(0.0,PCPDB(L,NY,NX))/VOLWPB 
      PCAPHB=AMAX1(0.0,PCPHB(L,NY,NX))/VOLWPB 
      ELSE
      RH1BX=0.0
      RH2BX=0.0
      CH1PB=0.0
      CH2PB=0.0
      XH01B=0.0
      XH11B=0.0
      XH21B=0.0
      X1P1B=0.0
      X2P1B=0.0
      PALPOB=0.0
      PFEPOB=0.0
      PCAPMB=0.0
      PCAPDB=0.0
      PCAPHB=0.0
      ENDIF
C
C     SOLUBLE NO3 CONCENTRATIONS
C     IN NON-BAND AND BAND SOIL ZONES
C
C     VOLWNO,VOLWNZ=soil water volume in NO3 non-band,band (m3)
C     ZNO3S,ZNO3B=NO3 mass in non-band,band (g N)
C     CNO1,CNOB=NO3 concentrations in non-band,band (mol N m-3)
C
      IF(VOLWNO.GT.ZEROS2(NY,NX))THEN
      CNO1=AMAX1(ZEROC,ZNO3S(L,NY,NX)/(14.0*VOLWNO))
      ELSE
      CNO1=0.0
      ENDIF
      IF(VOLWNZ.GT.ZEROS2(NY,NX))THEN
      CNOB=AMAX1(ZEROC,ZNO3B(L,NY,NX)/(14.0*VOLWNZ))
      ELSE
      CNOB=0.0
      ENDIF
C
C     CATION CONCENTRATIONS
C
C     CCEC,XCEC=total cation exchange concentration,capacity 
C        (mol Mg-1,mol)
C     BKVLX=soil mass (Mg) 
C     VOLWM=soil water volume (m3)
C     X*=cation exchange concentration (mol Mg-1)
C     cation code:HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C          :COOH=total carboxyl,HC=carboxyl + H+ 
C
      IF(BKVLX.GT.ZEROS(NY,NX).AND.VOLW(L,NY,NX).GT.ZEROS(NY,NX))THEN
      CCEC=AMAX1(ZEROC,XCEC(L,NY,NX)/BKVLX)
      XHY1=AMAX1(ZEROC,XHY(L,NY,NX)/BKVLX)
      XAL1=AMAX1(ZEROC,XAL(L,NY,NX)/BKVLX)
      XFE1=AMAX1(ZEROC,XFE(L,NY,NX)/BKVLX)
      XCA1=AMAX1(ZEROC,XCA(L,NY,NX)/BKVLX)
      XMG1=AMAX1(ZEROC,XMG(L,NY,NX)/BKVLX)
      XNA1=AMAX1(ZEROC,XNA(L,NY,NX)/BKVLX)
      XKA1=AMAX1(ZEROC,XKA(L,NY,NX)/BKVLX)
      XHC1=AMAX1(ZEROC,XHC(L,NY,NX)/BKVLX)
      XCOOH=AMAX1(ZEROC,COOH*1.0E-06*ORGC(L,NY,NX)/BKVLX)
      BKVLW=BKVLX/VOLW(L,NY,NX) 
      ELSE
      CCEC=ZERO
      XHY1=0.0
      XAL1=0.0
      XFE1=0.0
      XCA1=0.0
      XMG1=0.0
      XNA1=0.0
      XKA1=0.0
      XHC1=0.0
      XCOOH=0.0
      BKVLW=1.0
      ENDIF
C
C     IF SALT OPTION SELECTED IN SITE FILE
C     THEN SOLVE FULL SET OF EQUILIBRIA REACTIONS
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C
      IF(ISALTG.NE.0)THEN
C
C     C*=solute concentration (mol m-3)
C     Z*=solute mass (mol)
C     VOLW=soil water content (m3)
C     XZHYS=total H+ production from ‘nitro.f’ and ‘grosub.f’
C        (mol t-1)
C     solute code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+
C          :*CA*=Ca2+,*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-
C          :*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
C          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
C          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
C          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
C          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH
C          :*MGC*=MgCO3,*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-
C          :*NAS*=NaSO4-,*KAS*=KSO4-
C     anion code:OH=OH-,*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C          :*F2P*=F1H2PO4-,*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C          :*M1P*=MgHPO4,*COO*=COO-
C          :*1=non-band,*B=band
C
      IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
C
C     TIME STEPS
C
C     TPDXH,TADAXH,TADCXH,TSLXH,TRWXH=maximum 
C        precipitation,adsorption,solute,silicate weathwering
C        rates from PARAMETER (MRXN-1)
C     TPDX,TADAX,TADCX,TSLX,TRWX=maximum 
C        precipitation,adsorption,solute,silicate weathwering
C        rates (t-1)
C     XNFH=time step for solute transformations from ‘wthr.f’ (h t-1)
C     TFNW=temperature sensitivity of silicate weathering
C     SSA=soil surface area from ‘hour1.f’ (m2 m-3)
C     XZHYS,RZHYS=net H+ production vs consumption from nutrient and 
C        salt uptake from ‘nitro.f’ and ‘grosub.f’ 
C        (mol t-1,mol m-3 t-1)
C     VOLW=soil water content (m3)
C
      TPDX=TPDXH*XNFH 
      TPDZ=TPZXH*XNFH
      TADAX=TADAXH*XNFH 
      TADCX=TADCXH*XNFH 
      TSLX=TSLXH*XNFH 
      TSZX=TSZXH*XNFH 
      TRWX=TRWXH*XNFH 
      TPDA=TRWX 
      RTK=8.3143*TKS(L,NY,NX)
      STK=710.0*TKS(L,NY,NX)
      ACTV=1.0+EXP((197500-STK)/RTK)+EXP((STK-222500)/RTK)
      TFNW=EXP(25.229-62500/RTK)/ACTV
      TRWXA=TRWX*SSA(L,NY,NX)*TFNW
      RZHYS=XZHYS(L,NY,NX)/(XMRXN*VOLW(L,NY,NX))
C
C     INITIALIZE ION CONCENTRATIONS
C
      CHY1=AMAX1(ZEROC,ZHY(L,NY,NX)/VOLW(L,NY,NX))
      COH1=AMAX1(ZEROC,ZOH(L,NY,NX)/VOLW(L,NY,NX))
      CAL1=AMAX1(ZEROC,ZAL(L,NY,NX)/VOLW(L,NY,NX))
      CFE1=AMAX1(ZEROC,ZFE(L,NY,NX)/VOLW(L,NY,NX))
      CCA1=AMAX1(ZEROC,ZCA(L,NY,NX)/VOLW(L,NY,NX))
      CMG1=AMAX1(ZEROC,ZMG(L,NY,NX)/VOLW(L,NY,NX))
      CNA1=AMAX1(ZEROC,ZNA(L,NY,NX)/VOLW(L,NY,NX))
      CKA1=AMAX1(ZEROC,ZKA(L,NY,NX)/VOLW(L,NY,NX))
      CSO41=AMAX1(ZEROC,ZSO4(L,NY,NX)/VOLW(L,NY,NX))
      CCL1=AMAX1(ZEROC,ZCL(L,NY,NX)/VOLW(L,NY,NX))
      CCO31=AMAX1(ZEROC,ZCO3(L,NY,NX)/VOLW(L,NY,NX))
      CHCO31=AMAX1(ZEROC,ZHCO3(L,NY,NX)/VOLW(L,NY,NX))
      CCO21=AMAX1(ZEROC,CO2S(L,NY,NX)/(12.0*VOLW(L,NY,NX)))
      CALO1=AMAX1(ZEROC,ZALOH1(L,NY,NX)/VOLW(L,NY,NX))
      CALO2=AMAX1(ZEROC,ZALOH2(L,NY,NX)/VOLW(L,NY,NX))
      CALO3=AMAX1(ZEROC,ZALOH3(L,NY,NX)/VOLW(L,NY,NX))
      CALO4=AMAX1(ZEROC,ZALOH4(L,NY,NX)/VOLW(L,NY,NX))
      CALS1=AMAX1(ZEROC,ZALS(L,NY,NX)/VOLW(L,NY,NX))
      CFEO1=AMAX1(ZEROC,ZFEOH1(L,NY,NX)/VOLW(L,NY,NX))
      CFEO2=AMAX1(ZEROC,ZFEOH2(L,NY,NX)/VOLW(L,NY,NX))
      CFEO3=AMAX1(ZEROC,ZFEOH3(L,NY,NX)/VOLW(L,NY,NX))
      CFEO4=AMAX1(ZEROC,ZFEOH4(L,NY,NX)/VOLW(L,NY,NX))
      CFES1=AMAX1(ZEROC,ZFES(L,NY,NX)/VOLW(L,NY,NX))
      CCAO1=AMAX1(ZEROC,ZCAO(L,NY,NX)/VOLW(L,NY,NX))
      CCAC1=AMAX1(ZEROC,ZCAC(L,NY,NX)/VOLW(L,NY,NX))
      CCAH1=AMAX1(ZEROC,ZCAH(L,NY,NX)/VOLW(L,NY,NX))
      CCAS1=AMAX1(ZEROC,ZCAS(L,NY,NX)/VOLW(L,NY,NX))
      CMGO1=AMAX1(ZEROC,ZMGO(L,NY,NX)/VOLW(L,NY,NX))
      CMGC1=AMAX1(ZEROC,ZMGC(L,NY,NX)/VOLW(L,NY,NX))
      CMGH1=AMAX1(ZEROC,ZMGH(L,NY,NX)/VOLW(L,NY,NX))
      CMGS1=AMAX1(ZEROC,ZMGS(L,NY,NX)/VOLW(L,NY,NX))
      CNAC1=AMAX1(ZEROC,ZNAC(L,NY,NX)/VOLW(L,NY,NX))
      CNAS1=AMAX1(ZEROC,ZNAS(L,NY,NX)/VOLW(L,NY,NX))
      CKAS1=AMAX1(ZEROC,ZKAS(L,NY,NX)/VOLW(L,NY,NX))
      CHYSI1=AMAX1(ZEROC,ZHYSI(L,NY,NX)/VOLW(L,NY,NX))
      ELSE
      CHY1=10.0**(-PH(L,NY,NX))*1.0E+03
      COH1=DPH2O/CHY1
      CAL1=0.0
      CFE1=0.0
      CCA1=0.0
      CMG1=0.0
      CNA1=0.0
      CKA1=0.0
      CSO41=0.0
      CCL1=0.0
      CCO31=0.0
      CHCO31=0.0
      CCO21=0.0
      CALO1=0.0
      CALO2=0.0
      CALO3=0.0
      CALO4=0.0
      CALS1=0.0
      CFEO1=0.0
      CFEO2=0.0
      CFEO3=0.0
      CFEO4=0.0
      CFES1=0.0
      CCAO1=0.0
      CCAC1=0.0
      CCAH1=0.0
      CCAS1=0.0
      CMGO1=0.0
      CMGC1=0.0
      CMGH1=0.0
      CMGS1=0.0
      CNAC1=0.0
      CNAS1=0.0
      CKAS1=0.0
      CHYSI1=0.0
      ENDIF
      CHY0=CHY1
      COH0=COH1      
C
C     PO4 CONCENTRATIONS IN NON-BAND AND BAND SOIL ZONES
C
C     C*1,C*B=PO4 concentrations in non-band,band (mol m-3)
C     P code:
C        H0P=PO4,H3P=H3PO4,F1P=FeHPO4,F2P=FeH2PO4,C0P=CaPO4
C        C1P=CaHPO4,C2P=CaH2PO4,M1P=MgHPO4
C     VOLWPO,VOLWPB=water volume in PO4 non-band,band (m3)
C
      IF(VOLWPO.GT.ZEROS2(NY,NX))THEN
      CH0P1=AMAX1(ZEROC,H0PO4(L,NY,NX)/VOLWPO)
      CH3P1=AMAX1(ZEROC,H3PO4(L,NY,NX)/VOLWPO)
      CF1P1=AMAX1(ZEROC,ZFE1P(L,NY,NX)/VOLWPO)
      CF2P1=AMAX1(ZEROC,ZFE2P(L,NY,NX)/VOLWPO)
      CC0P1=AMAX1(ZEROC,ZCA0P(L,NY,NX)/VOLWPO)
      CC1P1=AMAX1(ZEROC,ZCA1P(L,NY,NX)/VOLWPO)
      CC2P1=AMAX1(ZEROC,ZCA2P(L,NY,NX)/VOLWPO)
      CM1P1=AMAX1(ZEROC,ZMG1P(L,NY,NX)/VOLWPO)
      ELSE
      CH0P1=0.0
      CH3P1=0.0
      CF1P1=0.0
      CF2P1=0.0
      CC0P1=0.0
      CC1P1=0.0
      CC2P1=0.0
      CM1P1=0.0
      ENDIF
      IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
      CH0PB=AMAX1(ZEROC,H0POB(L,NY,NX)/VOLWPB)
      CH3PB=AMAX1(ZEROC,H3POB(L,NY,NX)/VOLWPB)
      CF1PB=AMAX1(ZEROC,ZFE1PB(L,NY,NX)/VOLWPB)
      CF2PB=AMAX1(ZEROC,ZFE2PB(L,NY,NX)/VOLWPB)
      CC0PB=AMAX1(ZEROC,ZCA0PB(L,NY,NX)/VOLWPB)
      CC1PB=AMAX1(ZEROC,ZCA1PB(L,NY,NX)/VOLWPB)
      CC2PB=AMAX1(ZEROC,ZCA2PB(L,NY,NX)/VOLWPB)
      CM1PB=AMAX1(ZEROC,ZMG1PB(L,NY,NX)/VOLWPB)
      ELSE
      CH0PB=0.0
      CH3PB=0.0
      CF1PB=0.0
      CF2PB=0.0
      CC0PB=0.0
      CC1PB=0.0
      CC2PB=0.0
      CM1PB=0.0
      ENDIF
C
C     EXCHANGEABLE ION CONCENTRATIONS
C
      IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
C
C     PRECIPITATE CONCENTRATIONS
C
C     PALOH,PFEOH,PCACO,PCASO=AL(OH)3,FE(OH)3,CACO3,CASO4
C        mass (mol)
C     Q*SI=silicate content (mol)
C        silicate code: AL=Al,FE=Fe,CA=Ca,MG=Mg,NA=Na,KA=K
C
      PALOH1=AMAX1(0.0,PALOH(L,NY,NX))/VOLW(L,NY,NX)
      PFEOH1=AMAX1(0.0,PFEOH(L,NY,NX))/VOLW(L,NY,NX)
      PCACO1=AMAX1(0.0,PCACO(L,NY,NX))/VOLW(L,NY,NX)
      PCASO1=AMAX1(0.0,PCASO(L,NY,NX))/VOLW(L,NY,NX)
      QALSI1=AMAX1(0.0,QALSI(L,NY,NX))/VOLW(L,NY,NX)
      QFESI1=AMAX1(0.0,QFESI(L,NY,NX))/VOLW(L,NY,NX)
      QCASI1=AMAX1(0.0,QCASI(L,NY,NX))/VOLW(L,NY,NX)
      QMGSI1=AMAX1(0.0,QMGSI(L,NY,NX))/VOLW(L,NY,NX)
      QNASI1=AMAX1(0.0,QNASI(L,NY,NX))/VOLW(L,NY,NX)
      QKASI1=AMAX1(0.0,QKASI(L,NY,NX))/VOLW(L,NY,NX)
      ELSE
      PALOH1=0.0
      PFEOH1=0.0
      PCACO1=0.0
      PCASO1=0.0
      QALSI1=0.0 
      QFESI1=0.0 
      QCASI1=0.0 
      QMGSI1=0.0 
      QNASI1=0.0 
      QKASI1=0.0 
      ENDIF
C
C     ITERATIONS TOWARDS SOLUTE EQUILIBRIA
C
      DO 1000 M=1,MRXN
C
C     RESET PH AT START OF ITERATION
C
C     CHY1,COH1=H+,OH- concentration (mol m-3)
C     AHY2,AOH2=H+,OH- activity (mol m-3)
C     DPH2O=H2O solubility product from PARAMETER (mol m-3)2
C     A1=activity coefficient from ‘hour1.f’
C     RHHX=flux to equilibrate H+, OH- at DPH2O 
C     TRHY,TROH=total H+,OH- transformations used in ‘redist.f’
C     PH=soil pH
C
      IF(CHY1.GT.ZEROC)THEN
      AHY2=CHY1*A1(L,NY,NX)
      ELSE
      AHY2=CHY1
      ENDIF
      IF(COH1.GT.ZEROC)THEN
      AOH2=COH1*A1(L,NY,NX)
      ELSE
      AOH2=COH1
      ENDIF
      S0=AHY2+AOH2
      S1=AMAX1(0.0,S0**2-4.0*(AHY2*AOH2-DPH2O))
      RHHX=0.5*(S0-SQRT(S1)) 
      AHY2=AHY2-RHHX
      AOH2=AOH2-RHHX
      CHY1=AHY2/A1(L,NY,NX)
      COH1=AOH2/A1(L,NY,NX)
      TRHY(L,NY,NX)=TRHY(L,NY,NX)-RHHX
      TROH(L,NY,NX)=TROH(L,NY,NX)-RHHX
      TBH2O(L,NY,NX)=TBH2O(L,NY,NX)+RHHX
      PH(L,NY,NX)=-LOG10(AHY2*1.0E-03)
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     IF(L.EQ.1)THEN
C     WRITE(*,23)'RHHX',I,J,NFZ,NX,NY,L,M 
C    2,RHHX,CHY1,COH1,AHY2,AOH2,AHY2*AOH2 
C    3,DPH2O,ZHY(L,NY,NX),ZOH(L,NY,NX),TROH(L,NY,NX),PH(L,NY,NX)
C     ENDIF
C
C     END RESET PH AT START OF ITERATION
C
C     INITIALIZE ION ACTIVITIES
C
      AHY1=CHY1*A1(L,NY,NX)
      AOH1=COH1*A1(L,NY,NX)
      AAL1=CAL1*A3(L,NY,NX)
      AALO1=CALO1*A2(L,NY,NX)
      AALO2=CALO2*A1(L,NY,NX)
      AALO3=CALO3
      AALO4=CALO4*A1(L,NY,NX)
      AFE1=CFE1*A3(L,NY,NX)
      AFEO1=CFEO1*A2(L,NY,NX)
      AFEO2=CFEO2*A1(L,NY,NX)
      AFEO3=CFEO3
      AFEO4=CFEO4*A1(L,NY,NX)
      ACA1=CCA1*A2(L,NY,NX)
      ACO31=CCO31*A2(L,NY,NX)
      AHCO31=CHCO31*A1(L,NY,NX)
      ACO21=CCO21
      ASO41=CSO41*A2(L,NY,NX)
      AH0P1=CH0P1*A3(L,NY,NX)
      AH1P1=CH1P1*A2(L,NY,NX)
      AH2P1=CH2P1*A1(L,NY,NX)
      AH3P1=CH3P1
      AF1P1=CF1P1*A2(L,NY,NX)
      AF2P1=CF2P1*A1(L,NY,NX)
      AC0P1=CC0P1*A1(L,NY,NX)
      AC1P1=CC1P1
      AC2P1=CC2P1*A1(L,NY,NX)
      AM1P1=CM1P1
      AH0PB=CH0PB*A3(L,NY,NX)
      AH1PB=CH1PB*A2(L,NY,NX)
      AH2PB=CH2PB*A1(L,NY,NX)
      AH3PB=CH3PB
      AF1PB=CF1PB*A2(L,NY,NX)
      AF2PB=CF2PB*A1(L,NY,NX)
      AC0PB=CC0PB*A1(L,NY,NX)
      AC1PB=CC1PB
      AC2PB=CC2PB*A1(L,NY,NX)
      AM1PB=CM1PB
      AN41=CN41*A1(L,NY,NX)
      AN4B=CN4B*A1(L,NY,NX)
      AN31=CN31
      AN3B=CN3B
      AMG1=CMG1*A2(L,NY,NX)
      ANA1=CNA1*A1(L,NY,NX)
      AKA1=CKA1*A1(L,NY,NX)
      AALS1=CALS1*A1(L,NY,NX)
      AFES1=CFES1*A1(L,NY,NX)
      ACAO1=CCAO1*A1(L,NY,NX)
      ACAC1=CCAC1
      ACAS1=CCAS1
      ACAH1=CCAH1*A1(L,NY,NX)
      AMGO1=CMGO1*A1(L,NY,NX)
      AMGC1=CMGC1
      AMGH1=CMGH1*A1(L,NY,NX)
      AMGS1=CMGS1
      ANAC1=CNAC1*A1(L,NY,NX)
      ANAS1=CNAS1*A1(L,NY,NX)
      AKAS1=CKAS1*A1(L,NY,NX)
C
C     PRECIPITATION-DISSOLUTION REACTIONS FOR H2PO4 AND CO-PRECIPITATES
C     IN NON-BAND, BAND ZONES CALCULATED FROM ACTIVITIES OF REACTANTS
C     AND PRODUCTS THROUGH SOLUTIONS FOR THEIR EQUILIBRIUM
C     CONCENTRATIONS 
C
C     for all precipitation-dissolution reactions:
C
C     SP*=solubility product from PARAMETER above
C     A*1,A*B=ion activity in non-band,band (mol m-3)
C     A*Q,A*BQ=ion equilibrium activity at SP*in non-band,band 
C        (mol m-3)
C     A1,A2,A3=activity coefficients from ‘hour1.f’ 
C     R*X=precipitation(+ve) or dissolution (-ve) rate 
C        (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FION*=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     TPDX,TPDA=maximum precipitation rate (mol m-3 t-1)
C     salt solute code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+
C          :*CA*=Ca2+,*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-
C          :*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
C          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
C          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
C          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
C          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH
C          :*MGC*=MgCO3,*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-
C          :*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus solute code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C          :*F2P*=F1H2PO4-,*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C          :*M1P*=MgHPO4,*COO*=COO-
C          :*1=non-band,*B=band
C
C     PHOSPHORUS PRECIPITATION-DISSOLUTION IN NON-BAND SOIL ZONE
C
      IF(VOLWPO.GT.ZEROS2(NY,NX))THEN
C
C     ALUMINUM PHOSPHATE ALPO4 (VARISCITE)
C
      XMINN=FIONX*AH2P1
      XMINP=FIONX*AMIN1(AAL1,AH2P1)
      AAL1Q=SPALP/AH0P1
      RPALPX=AMAX1(-AMAX1(0.0,PALPO1),-TPDX,-XMINN
     2,(AMIN1(TPDX,XMINP,AAL1-AAL1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     WRITE(*,1112)'ALPO4',I,J,NFZ,NX,NY,L,M
C    2,RPALPX,AAL1,AAL1Q,PALPO1 
C    2,AH2P1,AH1P1,AH0P1,AAL1*AH0P1,SPALP,TPDX,XMINN,XMINP 
C    3,PALPO(L,NY,NX),PH(L,NY,NX),VLPO4(L,NY,NX)
1112  FORMAT(A8,7I4,50E12.4)
C     ENDIF
C
C     IRON PHOSPHATE FEPO4 (STRENGITE)
C
      XMINN=FIONX*AH2P1
      XMINP=FIONX*AMIN1(AFE1,AH2P1)
      AFE1Q=SPFEP/AH0P1
      RPFEPX=AMAX1(-AMAX1(0.0,PFEPO1),-TPDX,-XMINN
     2,(AMIN1(TPDX,XMINP,AFE1-AFE1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     WRITE(*,1112)'FEPO4',I,J,NFZ,NX,NY,L,M 
C    2,RPFEPX,AFE1,AFE1Q,PFEPO1 
C    2,AH2P1,AH1P1,AH0P1,AFE1*AH0P1,SPFEP,TPDX,XMINN,XMINP 
C    3,PFEPO(L,NY,NX),PH(L,NY,NX),VLPO4(L,NY,NX)
C     ENDIF
C
C     DICALCIUM PHOSPHATE CAHPO4
C
      XMINN=FIONX*AH2P1
      XMINP=FIONX*AMIN1(ACA1,AH2P1)
      ACA1Q=SPCAD/AH1P1 
      RPCADX=AMAX1(-AMAX1(0.0,PCAPD1),-TPDX,-XMINN
     2,(AMIN1(TPDX,XMINP,ACA1-ACA1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     WRITE(*,1112)'RPCADX',I,J,NFZ,NX,NY,L,M
C    2,RPCADX,ACA1,ACA1Q,AH1P1,ACA1*AH1P1,SPCAD
C    3,TPDX,PCAPD(L,NY,NX) 
C     ENDIF
C
C     HYDROXYAPATITE CA5(PO4)3OH
C
      XMINN=FIONX*AH2P1/3.0
      XMINP=FIONX*AMIN1(ACA1,AH2P1)/3.0
      ACA1Q=(SPCAH/(AH0P1**3*AOH1))**0.2
      RPCAHX=AMAX1(-AMAX1(0.0,PCAPH1),-TPDA,-XMINN
     2,(AMIN1(TPDA,XMINP,ACA1-ACA1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     WRITE(*,1115)'APAT',I,J,NFZ,NX,NY,L,M
C    2,RPCAHX,ACA1,ACA1Q,AH2P1,AH0P1,AOH1,COH1,A1(L,NY,NX) 
C    3,ACA1**5*AH0P1**3*AOH1,SPCAH,PCAPH(L,NY,NX) 
C    4,TPDA,XMINN,XMINP,PH(L,NY,NX),VLPO4(L,NY,NX)
C    3,CH0P1,CH1P1,CH2P1,XOH01,XOH11,XOH21,XH1P1,XH2P1
C    4,H1PO4(L,NY,NX),VOLWPX,RH1PX 
C    5,H2PO4(L,NY,NX),VOLWPX,RH2PX 
1115  FORMAT(A8,7I4,50E12.4)
C     ENDIF
C
C     MONOCALCIUM PHOSPHATE CA(H2PO4)2
C
      XMINP=FIONX*AMIN1(ACA1,AH2P1)
      ACA1Q=SPCAM/AH2P1**2
      RPCAMX=AMAX1(-AMAX1(0.0,PCAPM1),-TPDX 
     2,(AMIN1(TPDX,XMINP,ACA1-ACA1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     WRITE(*,1112)'RPCAMX',I,J,NFZ,NX,NY,L,M 
C    2,RPCAMX,ACA1,ACA1Q,AH2P1
C    3,ACA1*AH2P1**2,SPCAM,PCAPM(L,NY,NX) 
C     ENDIF
      ELSE
      RPALPX=0.0
      RPFEPX=0.0
      RPCADX=0.0
      RPCAHX=0.0
      RPCAMX=0.0
      ENDIF
C
C     PHOSPHORUS PRECIPITATION-DISSOLUTION IN BAND SOIL ZONE
C
      IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
C
C     ALUMINUM PHOSPHATE ALPO4 (VARISCITE)
C
      XMINN=FIONX*AH2PB
      XMINP=FIONX*AMIN1(AAL1,AH2PB)
      AAL1Q=SPALP/AH0PB
      RPALBX=AMAX1(-AMAX1(0.0,PALPOB),-TPDX,-XMINN
     2,(AMIN1(TPDX,XMINP,AAL1-AAL1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     IF(L.EQ.2)THEN
C     WRITE(*,1112)'ALPOB',I,J,NFZ,NX,NY,L,M
C    2,RPALBX,AAL1,AAL1Q,PALPOB 
C    2,AH2PB,AAL1*AH0PB,SPALP,TPDX,XMINN,XMINP 
C    3,PALPB(L,NY,NX),PH(L,NY,NX),VLPOB(L,NY,NX)
C     ENDIF
C
C     IRON PHOSPHATE FEPO4 (STRENGITE)
C
      XMINN=FIONX*AH2PB
      XMINP=FIONX*AMIN1(AFE1,AH2PB)
      AFE1Q=SPFEP/AH0PB
      RPFEBX=AMAX1(-AMAX1(0.0,PFEPOB),-TPDX,-XMINN
     2,(AMIN1(TPDX,XMINP,AFE1-AFE1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     WRITE(*,1112)'FEPOB',I,J,NFZ,NX,NY,L,M 
C     2,RPFEBX,AFE1,AFE1Q,PFEPOB 
C    2,AH2PB,AFE1*AH0PB,SPFEP,TPDX,XMINN,XMINP 
C    3,PFEPB(L,NY,NX),PH(L,NY,NX),VLPOB(L,NY,NX)
C     ENDIF
C
C     DICALCIUM PHOSPHATE CAHPO4
C
      XMINN=FIONX*AH2PB
      XMINP=FIONX*AMIN1(ACA1,AH2PB)
      ACA1Q=SPCAD/AH1PB 
      RPCDBX=AMAX1(-AMAX1(0.0,PCAPDB),-TPDX,-XMINN
     2,(AMIN1(TPDX,XMINP,ACA1-ACA1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     WRITE(*,1112)'RPCDBX',I,J,NFZ,NX,NY,L,M
C    2,RPCDBX,ACA1,ACA1Q,AH1PB,ACA1*AH1PB 
C    3,SPCAD,TPDX,XMINN,XMINP,PCPDB(L,NY,NX)
C     ENDIF
C
C     HYDROXYAPATITE CA5(PO4)3OH
C
      XMINN=FIONX*AH2PB/3.0
      XMINP=FIONX*AMIN1(ACA1,AH2PB)/3.0
      ACA1Q=(SPCAH/(AH0PB**3*AOH1))**0.2
      RPCHBX=AMAX1(-AMAX1(0.0,PCAPHB),-TPDA,-XMINN
     2,(AMIN1(TPDA,XMINP,ACA1-ACA1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     WRITE(*,1112)'APATB',I,J,NFZ,NX,NY,L,M
C    2,RPCHBX,ACA1,ACA1Q,AH2PB,AH0PB,AOH1
C    3,ACA1**5*AH0PB**3*AOH1,SPCAH,PCPHB(L,NY,NX)
C    4,TPDZ,XMINN,XMINP,PH(L,NY,NX)
C     ENDIF 
C
C     MONOCALCIUM PHOSPHATE CA(H2PO4)2
C
      XMINP=FIONX*AMIN1(ACA1,AH2PB)
      ACA1Q=SPCAM/AH2PB**2
      RPCMBX=AMAX1(-AMAX1(0.0,PCAPMB),-TPDX
     2,(AMIN1(TPDX,XMINP,ACA1-ACA1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     WRITE(*,1112)'RPCMBX',I,J,NFZ,NX,NY,L,M
C    2,RPCMBX,ACA1,ACA1Q,AH2PB
C    3,ACA1*AH2PB**2,SPCAM
C    4,RH1BX,RH2BX,TUPH1B(L,NY,NX),TUPH2B(L,NY,NX) 
C    5,TPDX,XMINP,PALPB(L,NY,NX),PFEPB(L,NY,NX)
C    6,PCPDB(L,NY,NX),3.0*PCPHB(L,NY,NX),2.0*PCPMB(L,NY,NX)
C    6,PALPO(L,NY,NX),PFEPO(L,NY,NX)
C    6,PCAPD(L,NY,NX),3.0*PCAPH(L,NY,NX),2.0*PCAPM(L,NY,NX)
C    7,VLPOB(L,NY,NX)
C     ENDIF 
      ELSE
      RPALBX=0.0
      RPFEBX=0.0
      RPCDBX=0.0
      RPCHBX=0.0
      RPCMBX=0.0
      ENDIF
C
C     PHOSPHORUS ANION EXCHANGE IN NON-BAND SOIL ZONE
C     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
C     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
C     EXCHANGE SITES
C
C     XAEC=anion exchange capacity (mol)
C     VOLWPO=soil water volume in non-band (m3)
C     TADAX=maximum anion adsorption rate (mol m-3 t-1)
C     RXOH2,RXOH1=OH2,OH adsorption(+ve) or desorption(-ve) 
C        with R-OH2,R-OH in non-band (mol P m-3 t-1)
C     AHY1=H+ activity in non-band (mol H m-3)
C     SXOH2,SXOH1=equilibrium constant for OH2,OH exchange with 
C        R-OH2,R-OH from PARAMETER
C     XOH01,XOH11,XOH21=concentration of adsorption sites R-,R-OH,
C        R-OH2 in non-band (mol Mg-1)
C     XOH0Q,XOH1Q=equilibrium concentration of adsorption sites 
C        R-,R-OH in non-band at SX* (mol Mg-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONX=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C
      IF(VOLWPO.GT.ZEROS2(NY,NX)
     2.AND.XAEC(L,NY,NX).GT.ZEROS(NY,NX))THEN
      XMINN=FIONX*XOH21
      XMINP=FIONX*AMIN1(AHY1,XOH11)
      XOH1Q=SXOH2*XOH21/AHY1
      RXOH2=AMAX1(-TADAX,-XMINN
     2,(AMIN1(TADAX,XMINP,XOH11-XOH1Q)))
      XMINN=FIONX*XOH11
      XMINP=FIONX*AMIN1(AHY1,XOH01)
      XOH0Q=SXOH1*XOH11/AHY1
      RXOH1=AMAX1(-TADAX,-XMINN
     2,(AMIN1(TADAX,XMINP,XOH01-XOH0Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.L.EQ.2)THEN
C     WRITE(*,23)'RXOH2',I,J,NFZ,NX,NY,L,M
C    2,RXOH2,XOH11,XOH1Q,AHY1,SXOH2,XOH21 
C    3,RXOH1,XOH01,XOH0Q,AHY1,SXOH1,XOH11,TADAX,AHY1
C     ENDIF
C
C     H2PO4 EXCHANGE IN NON-BAND SOIL ZONE FROM EQUILIBRIUM
C     AMONG H2PO4-, H+, OH-, X-OH AND X-H2PO4
C
C     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with 
C        R-OH2,R-OH from PARAMETER
C     RXH2P,RYH2P=H2PO4 adsorption(+ve), desorption(-ve) 
C        with R-OH2,R-OH in non-band (mol P m-3 t-1)
C     AH2P1=H2PO4 activity in non-band (mol P m-3)
C     AH2PQ2,AH2PQ1=H2PO4 equilibrium activity in non-band at SP*,SX*
C        (mol P m-3)
C     XH2P1=exchangeable H2PO4 concentration in non-band (mol Mg-1)
C     XOH11,XOH21=concentration of adsorption sites R-OH,
C        R-OH2 in non-band (mol Mg-1)
C     AHY1,AOH1=H+,OH- activity (mol m-3)
C     TADAX=maximum anion adsorption rate (mol m-3 t-1)
C     XMINP*,XMINN*=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONX=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C
      XMINN2=FIONX*XH2P1 
      XMINP2=FIONX*AMIN1(AH2P1,XOH21) 
      SPH2P=SXH2P*AHY1*AOH1
      AH2PQ2=SPH2P*XH2P1/XOH21
      RXH2P=AMAX1(-TADAX,-XMINN2
     2,(AMIN1(TADAX,XMINP2,AH2P1-AH2PQ2)))
      XMINN1=FIONX*XH2P1 
      XMINP1=FIONX*AMIN1(AH2P1,XOH11) 
      SYH2P=SXH2P*AOH1
      AH2PQ1=SYH2P*XH2P1/XOH11
      RYH2P=AMAX1(-TADAX,-XMINN1
     2,(AMIN1(TADAX,XMINP1,AH2P1-AH2PQ1)))
C
C     HPO4 EXCHANGE IN NON-BAND SOIL ZONE FROM EQUILIBRIUM
C     AMONG HPO4--, H+, OH-, X-OH AND X-HPO4
C
C     SPH1P=equilibrium constant for HPO4 exchange with R-OH 
C        from PARAMETER
C     RXH1P=HPO4 adsorption(+ve), desorption(-ve) with R-OH 
C        in non-band (mol P m-3 t-1)
C     AH1P1=HPO4 activity in non-band (mol P m-3)
C     AH1PQ=HPO4 equilibrium activity in non-band at SP*(mol P m-3)
C     XH1P1=exchangeable HPO4 concentration in non-band (mol P Mg-1)
C     XOH11=concentration of adsorption sites R-OH
C        in non-band (mol Mg-1)
C     TADAX=maximum anion adsorption rate (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONX=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C
      XMINN=FIONX*XH1P1
      XMINP=FIONX*AMIN1(AH1P1,XOH11) 
      SPH1P=SXH1P*DPH2O/DPH2P
      AH1PQ=SPH1P*XH1P1/XOH11
      RXH1P=AMAX1(-TADAX,-XMINN
     2,(AMIN1(TADAX,XMINP,AH1P1-AH1PQ)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN 
C     WRITE(*,23)'RXH2P',I,J,NFZ,NX,NY,L,M
C    2,RXH2P,TADAX,XOH21,AH2P1,AH2PQ2,AH2P1*XOH21/XH2P1,SPH2P 
C    2,XH2P1,BKVLW,CH2P1,A1(L,NY,NX)
C    2,H2PO4(L,NY,NX),RH2PX,VOLWPX,XMINN2,XMINP2,TADAX
C    2,XOH01+XOH11+XOH21+XH2P1+XH1P1
C    2,XOH01,XOH11,XOH21,XH2P1,XH1P1
C     WRITE(*,23)'RYH2P',I,J,NFZ,NX,NY,L,M
C    2,RYH2P,TADAX,XOH11,AH2P1,AH2PQ1,AH2P1*XOH11/XH2P1,SYH2P 
C    2,XH2P1,AOH1,BKVLW,CH2P1,A1(L,NY,NX)
C    2,H2PO4(L,NY,NX),RH2PX,VOLWPX,XMINN1,XMINP1,TADAX
C    2,XOH01+XOH11+XOH21+XH2P1+XH1P1
C     WRITE(*,23)'RXH1P',I,J,NFZ,NX,NY,L,M
C    2,RXH1P,TADAX,XOH11,AH1P1,AH1PQ,AH1P1*XOH11/XH1P1,SPH1P 
C    2,XH1P1,BKVLW,CH1P1,A1(L,NY,NX)
C    2,H1PO4(L,NY,NX),RH1PX,VOLWPX,XMINN,XMINP,TADAX
C    2,XOH01+XOH11+XOH21+XH2P1+XH1P1
C     ENDIF
      ELSE
      RXOH2=0.0
      RXOH1=0.0
      RXH2P=0.0
      RYH2P=0.0
      RXH1P=0.0
      ENDIF
C
C     PHOSPHORUS ANION EXCHANGE IN BAND SOIL ZONE
C     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
C     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
C     EXCHANGE SITES
C
C     XAEC=anion exchange capacity (mol)
C     VOLWPB=soil water volume in band (m3)
C     RXO2B,RXO1B=OH2,OH adsorption(+ve), desorption(-ve) with 
C        R-OH2,R-OH in band (mol OH m-3 t-1)
C     AHY1=H+ activity (mol H m-3)
C     XOH1B,XH11B,XH21B=concentration of adsorption sites R-,R-OH,
C        R-OH2 in band (mol OH Mg-1)
C     XOH1BQ,XH11BQ=equilibrium concentration of adsorption sites 
C        R-,R-OH in band at SP*(mol OH Mg-1)
C     SXOH2,SXOH1=equilibrium constant for OH2,OH exchange with 
C        R-OH2,R-OH from PARAMETER 
C     TADAX=maximum anion adsorption rate (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONX=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C
      IF(VOLWPB.GT.ZEROS2(NY,NX)
     2.AND.XAEC(L,NY,NX).GT.ZEROS(NY,NX))THEN
      XMINN=FIONX*XH21B
      XMINP=FIONX*AMIN1(AHY1,XH11B)
      XH11BQ=SXOH2*XH21B/AHY1
      RXO2B=AMAX1(-TADAX,-XMINN
     2,(AMIN1(TADAX,XMINP,XH11B-XH11BQ)))
      XMINN=FIONX*XH11B
      XMINP=FIONX*AMIN1(AHY1,XH01B)
      XH01BQ=SXOH1*XH11B/AHY1
      RXO1B=AMAX1(-TADAX,-XMINN
     2,(AMIN1(TADAX,XMINP,XH01B-XH01BQ)))
C
C     H2PO4 EXCHANGE IN BAND SOIL ZONE FROM SOLUTION FOR 
C     EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH AND X-H2PO4
C
C     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with 
C        R-OH2,R-OH from PARAMETER 
C     RXH2B,RYH2B=H2PO4 adsorption(+ve), desorption(-ve) with 
C        R-OH2,R-OH in band (mol P m-3 t-1)
C     AH2PB=H2PO4 activity in band (mol P m-3)
C     AH2BQ2,AH2BQ1=H2PO4 equilibrium activity in band at SP*
C        (mol P m-3)
C     X2P1B=exchangeable H2PO4 concentration in band (mol P Mg-1)
C     XH11B,XH21B=concentration of adsorption sites R-OH,
C        R-OH2 in band (mol OH Mg-1)
C     AOH1=OH activity (mol m-3)
C     TADAX=maximum anion adsorption rate (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONX=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C
      XMINN2=FIONX*X2P1B 
      XMINP2=FIONX*AMIN1(AH2PB,XH21B) 
      SPH2P=SXH2P*DPH2O 
      AH2BQ2=SPH2P*X2P1B/XH21B
      RXH2B=AMAX1(-TADAX,-XMINN2
     2,(AMIN1(TADAX,XMINP2,AH2PB-AH2BQ2)))
      XMINN1=FIONX*X2P1B 
      XMINP1=FIONX*AMIN1(AH2PB,XH11B) 
      SYH2P=SXH2P*AOH1
      AH2BQ1=SYH2P*X2P1B/XH11B
      RYH2B=AMAX1(-TADAX,-XMINN1
     2,(AMIN1(TADAX,XMINP1,AH2PB-AH2BQ1)))
C
C     HPO4 EXCHANGE IN BAND SOIL ZONE FROM SOLUTION 
C     FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH AND X-HPO4
C
C     SPH1P=equilibrium constant for HPO4 exchange with R-OH
C     RXH1B=HPO4 adsorption(+ve),desorption(-ve) with R-OH 
C        in band (mol P m-3 t-1)
C     AH1PB=HPO4 activity in band (mol H m-3)
C     AH1BQ=HPO4 equilibrium activity in band at SP* (mol H m-3)
C     XH11B=concentration of adsorption sites R-OH,
C        in band (mol OH Mg-1
C     X1P1B=exchangeable HPO4 concentration in band (mol P Mg-1)
C     TADAX=maximum anion adsorption rate (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONX=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C
      XMINN=FIONX*X1P1B 
      XMINP=FIONX*AMIN1(AH1PB,XH11B) 
      SPH1P=SXH1P*DPH2O/DPH2P
      AH1BQ=SPH1P*X1P1B/XH11B
      RXH1B=AMAX1(-TADAX,-XMINN
     2,(AMIN1(TADAX,XMINP,AH1PB-AH1BQ)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     WRITE(*,3334)'RXH2B',I,J,NFZ,NX,NY,L,M
C    2,RXH2B,TADAX,XH21B,AH2PB,AH2BQ2,SPH2P,X2P1B
C    3,FIONX*AH2PB,BKVLW,VOLWPB,XH2PB(L,NY,NX) 
C    2,XH01B,XH11B,XH21B,X2P1B,X1P1B
C    2,XH01B+XH11B+XH21B+X2P1B+X1P1B,XMINN2,XMINP2,TADAX 
C     WRITE(*,3334)'RYH2B',I,J,NFZ,NX,NY,L,M
C    2,RYH2B,TADAX,XH11B,AH2PB,AH2BQ1,SYH2P,X2P1B 
C    3,FIONX*AH2PB,BKVLW,VOLWPB,XH2PB(L,NY,NX),AOH1 
C    2,XH01B,XH11B,XH21B,X2P1B,X1P1B
C    2,XH01B+XH11B+XH21B+X2P1B+X1P1B,XMINN1,XMINP1,TADAX 
C     WRITE(*,3334)'RXH1B',I,J,NFZ,NX,NY,L,M
C    2,RXH1B,TADAX,X1P1B,AH1PB,AH1BQ,SPH1P,X1P1B 
C    3,FIONX*AH1PB,BKVLW,VOLWPB,XH1PB(L,NY,NX),PH(L,NY,NX) 
C    2,XH01B,XH11B,XH21B,X2P1B,X1P1B
C    2,XH01B+XH11B+XH21B+X2P1B+X1P1B,XMINN,XMINP,TADAX
3334  FORMAT(A8,7I4,40E12.4)
C     ENDIF
      ELSE
      RXO2B=0.0
      RXO1B=0.0
      RXH2B=0.0
      RYH2B=0.0
      RXH1B=0.0
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.1.AND.NFZ.EQ.1)THEN
C     WRITE(*,1112)'RXH2PS',I,J,NFZ,NX,NY,L,M
C    2,RXH2P,RYH2P,RXH1P,RXOH2,RXOH1,RXH2B,RYH2B,RXH1B,RXO2B,RXO1B
C    2,RXH2P+RYH2P+RXH1P+RXOH2+RXOH1,RXH2B+RYH2B+RXH1B+RXO2B+RXO1B
C    3,XH1P1,XH2P1,XOH01,XOH11,XOH21,X1P1B,X2P1B,XH01B,XH11B,XH21B
C    3,XH1P1+XH2P1+XOH01+XOH11+XOH21,X1P1B+X2P1B+XH01B+XH11B+XH21B
C    2,XAEC(L,NY,NX),VLPOB(L,NY,NX)
C     ENDIF
C
C     CATION EXCHANGE FROM GAPON SELECTIVITY COEFFICIENTS
C     FOR CA-NH4, CA-H, CA-AL, CA-MG, CA-NA, CA-K
C
      IF(XCEC(L,NY,NX).GT.ZEROS(NY,NX))THEN
C
C     CCEC,XCEC=cation exchange concentration,capacity (mol Mg-1,mol)
C     A*1=cation solute activity (mol m-3)
C     cation code:AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C          :N4=NH4+ non-band,NB=NH4+ band,HY=H+
C
      AALX=AMAX1(ZEROC,AAL1)**0.333
      AFEX=AMAX1(ZEROC,AFE1)**0.333
      ACAX=AMAX1(ZEROC,ACA1)**0.500
      AMGX=AMAX1(ZEROC,AMG1)**0.500
C
C     EQUILIBRIUM X-CA CONCENTRATION FROM CEC AND CATION
C     CONCENTRATIONS
C
C     CCEC,XCEC=total cation exchange concentration,capacity 
C        (mol Mg-1,mol)
C     XCAX=equilibrium X-Ca concentration (mol+ Mg-1)
C     A*1=cation activity (mol m-3)
C     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
C        CA-NH4,CA-H,CA-AL,CA-MG,CA-NA,CA-K
C     X*Q,X*Y=equilibrium,current exchangeable cation concentrations
C        (mol Mg-1)
C     XTLQ,XTL1=total equilibrium exchangeable concentration 
C        (mol, mol(+) Mg-1
C     cation code:AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ non-band,NB=NH4+ band,HY=H+
C 
      IF(ACAX.GT.ZERO.AND.CCEC.GT.ZERO)THEN
      XCAX=CCEC/(1.0
     2+GKC4(L,NY,NX)*AN41/ACAX*VLNH4(L,NY,NX)
     3+GKC4(L,NY,NX)*AN4B/ACAX*VLNHB(L,NY,NX)
     4+GKCH(L,NY,NX)*AHY1/ACAX
     5+GKCA(L,NY,NX)*AALX/ACAX*3.0
     6+GKCA(L,NY,NX)*AFEX/ACAX*3.0
     7+GKCM(L,NY,NX)*AMGX/ACAX*2.0
     8+GKCN(L,NY,NX)*ANA1/ACAX
     9+GKCK(L,NY,NX)*AKA1/ACAX)
      XN4Q=XCAX*AN41/ACAX*GKC4(L,NY,NX)*VLNH4(L,NY,NX)
      XNBQ=XCAX*AN4B/ACAX*GKC4(L,NY,NX)*VLNHB(L,NY,NX)
      XHYQ=XCAX*AHY1/ACAX*GKCH(L,NY,NX)
      XALQ=XCAX*AALX/ACAX*GKCA(L,NY,NX)*3.0
      XFEQ=XCAX*AFEX/ACAX*GKCA(L,NY,NX)*3.0
      XCAQ=XCAX*ACAX/ACAX*2.0
      XMGQ=XCAX*AMGX/ACAX*GKCM(L,NY,NX)*2.0
      XNAQ=XCAX*ANA1/ACAX*GKCN(L,NY,NX)
      XKAQ=XCAX*AKA1/ACAX*GKCK(L,NY,NX)
      XTLQ=XN4Q+XNBQ+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ 
      XTL1=XN41+XNB1+XHY1+XAL1*3.0+XFE1*3.0
     2+XCA1*2.0+XMG1*2.0+XNA1+XKA1 
      IF(XTLQ.GT.ZERO)THEN
      FX=CCEC/XTLQ
      ELSE
      FX=0.0
      ENDIF
      IF(XTL1.GT.ZERO)THEN
      FY=CCEC/XTL1
      ELSE
      FY=0.0
      ENDIF
      XN4Q=FX*XN4Q
      XNBQ=FX*XNBQ
      XHYQ=FX*XHYQ
      XALQ=FX*XALQ
      XFEQ=FX*XFEQ
      XCAQ=FX*XCAQ
      XMGQ=FX*XMGQ
      XNAQ=FX*XNAQ
      XKAQ=FX*XKAQ
      XN4Y=FY*XN41
      XNBY=FY*XNB1
      XHYY=FY*XHY1
      XALY=FY*XAL1*3.0
      XFEY=FY*XFE1*3.0
      XCAY=FY*XCA1*2.0
      XMGY=FY*XMG1*2.0
      XNAY=FY*XNA1
      XKAY=FY*XKA1
C
C     CATION EXCHANGE IN NON-BAND AND BAND SOIL ZONES
C
C     FIONC=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     BKVLW=soil mass:water (Mg m-3)
C     TADCC=maximum cation adsorption rate (mol Mg t-1)
C     XMIN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C
C     for NH4+
C
C     RXN4Q,RXNBQ=NH4 equilibrium desorption(-ve), adsorption(+ve) 
C        rate in non-band,band (mol N m-3 t-1)
C     RXN4,RXNBQ=NH4 current desorption(-ve), adsorption(+ve)
C        rate in non-band,band (mol N m-3 t-1)
C     XN4Q,XN4Y=equilibrium,current exchangeable NH4 concentration
C        non-band (mol Mg-1)
C     CN41,AN41=aqueous NH4 concentration,activity non-band (mol N m-3)
C     XNBQ,XNBY=equilibrium,current exchangeable NH4 concentration
C        band (mol N Mg-1)
C     CN4B,AN4B=aqueous NH4 concentration,activity band (mol N m-3)
C
C     for all other cations:
C
C     RX*Q,RX*=equilibrium,current cation desorption(-ve),
C        adsorption(+ve) rate (mol m-3 t-1)
C     X*Q,X*Y=equilibrium,current exchangeable cation concentrations
C        (mol Mg-1)
C     C*1,A*1=aqueous cation concentration,activity (mol m-3) 
C     cation code
C        :HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ non-band,NB=NH4+ band
C
      FIONC=FIONX/BKVLW
      TADCC=TADCX/BKVLW
      XMIN=FIONC*AMIN1(XN41,CN41)
      RXN4Q=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XN4Q-XN4Y)))*VLNH4(L,NY,NX)
      XMIN=FIONC*AMIN1(XNB1,CN4B)
      RXNBQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XNBQ-XNBY)))*VLNHB(L,NY,NX)
      XMIN=FIONC*AMIN1(XHY1,CHY1)
      RXHYQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XHYQ-XHYY)))
      XMIN=FIONC*3.0*AMIN1(XAL1,CAL1)
      RXALQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XALQ-XALY)))
      XMIN=FIONC*3.0*AMIN1(XFE1,CFE1)
      RXFEQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XFEQ-XFEY)))
      XMIN=FIONC*2.0*AMIN1(XCA1,CCA1)
      RXCAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XCAQ-XCAY)))
      XMIN=FIONC*2.0*AMIN1(XMG1,CMG1)
      RXMGQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XMGQ-XMGY)))
      XMIN=FIONC*AMIN1(XNA1,CNA1)
      RXNAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XNAQ-XNAY)))
      XMIN=FIONC*AMIN1(XKA1,CKA1)
      RXKAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XKAQ-XKAY)))
      TXXX=RXN4Q+RXNBQ+RXHYQ+RXALQ+RXFEQ+RXCAQ+RXMGQ+RXNAQ+RXKAQ
      TXXY=ABS(RXN4Q)+ABS(RXNBQ)+ABS(RXHYQ)+ABS(RXALQ)+ABS(RXFEQ)
     3+ABS(RXCAQ)+ABS(RXMGQ)+ABS(RXNAQ)+ABS(RXKAQ)
      RXN4=RXN4Q-TXXX*ABS(RXN4Q)/TXXY
      RXNB=RXNBQ-TXXX*ABS(RXNBQ)/TXXY
      RXHY=RXHYQ-TXXX*ABS(RXHYQ)/TXXY
      RXAL=(RXALQ-TXXX*ABS(RXALQ)/TXXY)/3.0
      RXFE=(RXFEQ-TXXX*ABS(RXFEQ)/TXXY)/3.0
      RXCA=(RXCAQ-TXXX*ABS(RXCAQ)/TXXY)/2.0
      RXMG=(RXMGQ-TXXX*ABS(RXMGQ)/TXXY)/2.0
      RXNA=RXNAQ-TXXX*ABS(RXNAQ)/TXXY
      RXKA=RXKAQ-TXXX*ABS(RXKAQ)/TXXY
      ELSE
      RXN4=0.0
      RXNB=0.0
      RXHY=0.0
      RXAL=0.0
      RXFE=0.0
      RXCA=0.0
      RXMG=0.0
      RXNA=0.0
      RXKA=0.0
      ENDIF
      ELSE
      RXN4=0.0
      RXNB=0.0
      RXHY=0.0
      RXAL=0.0
      RXFE=0.0
      RXCA=0.0
      RXMG=0.0
      RXNA=0.0
      RXKA=0.0
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.L.EQ.1)THEN
C     WRITE(*,2221)'RXN4S',I,J,NFZ,L,M
C    2,XN4Q,XHYQ,XALQ,XFEQ,XCAQ,XMGQ,XNAQ,XKAQ,XNBQ
C    2,XN41,XHY1,XAL1,XFE1,XCA1,XMG1,XNA1,XKA1,XNB1
C    2,XN4Y,XHYY,XALY,XFEY,XCAY,XMGY,XNAY,XKAY,XNBY
C    2,AN41,AHY1,AAL1,AFE1,ACA1,AMG1,ANA1,AKA1,AN4B
C    2,RXN4,RXHY,RXAL,RXFE,RXCA,RXMG,RXNA,RXKA,RXNB
C    2,RXN4+RXNB
C    2+RXHY+RXAL*3.0+RXFE*3.0+RXCA*2.0+RXMG*2.0+RXNA+RXKA
C    2,XN4Y*VLNH4(L,NY,NX)+XNBY*VLNHB(L,NY,NX)
C    2+XHYY+XALY+XFEY+XCAY+XMGY+XNAY+XKAY
C    2,CCEC,XCAX,ACAX,XTLQ,FX,XTL1,FY
C    2,TADCX,TXXX,TXXY,XCEC(L,NY,NX),BKVLX,CN41
C    2,ZNH4S(L,NY,NX),VOLWNX,RN4X 
2221  FORMAT(A8,5I4,/,9F16.8,/,9F16.8,/,9F16.8,/,9F16.8,/
     2,9F16.8,/,9F16.8,/,9F16.8,/,9F16.8,/,9F16.8)
C     ENDIF
C
C     DISSOCIATION OF CARBOXYL RADICALS
C
C     TADCC=maximum cation adsorption rate (mol Mg t-1)
C     XCOOH,XHC1=total carboxyl exchange sites,occupied by H+ (mol)
C     RXHC=COOH-COO+H desorption(-ve) or adsorption(+ve)(mol C m-3 t-1)
C     XMIN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONX=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     BKVLW=soil mass:water (Mg m-3)
C
      XMIN=FIONX/BKVLW*XHC1
      TADCC=TADCX/BKVLW
      XCOO=AMAX1(ZEROC,XCOOH-XHC1)
      XCOOQ=AMIN1(XCOOH,DPCOH*XHC1/AHY1)
      RXHC=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XCOO-XCOOQ))) 
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,23)'RXHC',I,J,NFZ,NX,NY,L,M
C    2,RXHC,RXHC*BKVLW,XHC1,XCOO,XCOOQ,XCOOH,CHY1,AHY1,XCOO*AHY1/XHC1 
C    3,DPCOH,COOH,TADCC,XMIN,XHC(L,NY,NX),BKVLW 
C     ENDIF
C
C     ION PAIRING REACTIONS
C
C     for all reactions:
C
C     DP*=dissociation constant from PARAMETER
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONX=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     A*1=ion activity (mol m-3)
C     A*Q=ion equilibrium activity at DP* (mol m-3)
C     TSLX=maximum reaction rate (mol m-3 t-1)
C     R*=dissociation(-ve) or association(+ve) rate (mol m-3 t-1)
C     saltsolute code
C        :*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+
C        :*CA*=Ca2+,*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-
C        :*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
C        :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
C        :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
C        :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
C        :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH
C        :*MGC*=MgCO3,*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-
C        :*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphate anion code
C        :OH=OH-,*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COO-
C        :*1=non-band,*B=band
C
C     NH4=NH3+H dissociation
C    
      IF(VOLWNH.GT.ZEROS2(NY,NX))THEN
      XMINN=FIONH*AN41
      XMINP=FIONH*AMIN1(AHY1,AN31)
      AN3Q=DPN4*AN41/AHY1
      RNH4=AMAX1(-TSZX,-XMINN
     2,(AMIN1(TSZX,XMINP,AN31-AN3Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,23)'RNH4',I,J,NFZ,NX,NY,L,M 
C    2,RNH4,AN41,AN31,AN3Q,AHY1,AN31*AHY1/AN41,DPN4
C    4,CHY1,XMINN,XMINP,TSZX,PH(L,NY,NX) 
C     ENDIF
      ELSE
      RNH4=0.0
      ENDIF
      IF(VOLWNB.GT.ZEROS2(NY,NX))THEN
      XMINN=FIONX*AN4B
      XMINP=FIONX*AMIN1(AHY1,AN3B)
      AN3BQ=DPN4*AN4B/AHY1
      RNHB=AMAX1(-TSZX,-XMINN
     2,(AMIN1(TSZX,XMINP,AN3B-AN3BQ)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,23)'RNHB',I,J,NFZ,NX,NY,L,M 
C    3,RNHB,AN4B,AN3B,AN3BQ,AHY1,AN3B*AHY1/AN4B,DPN4
C    4,CHY1,XMINN,XMINP,TSZX,PH(L,NY,NX) 
C     ENDIF
      ELSE
      RNHB=0.0
      ENDIF
      CHYX=CHY1-RNH4*VLNH4(L,NY,NX)-RNHB*VLNHB(L,NY,NX)
      AHY1=CHYX*A1(L,NY,NX)
C
C     RHCO3=HCO3-CO3+H dissociation
C
      XMINN=FIONX*AHCO31
      XMINP=FIONX*AMIN1(AHY1,ACO31)
      ACO3Q=DPHCO*AHCO31/AHY1
      RHCO3=AMAX1(-TSLX,-XMINN 
     2,(AMIN1(TSLX,XMINP,ACO31-ACO3Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,23)'RHCO3',I,J,NFZ,NX,NY,L,M
C    2,RHCO3,CHY1,CHYX,AHY1,AHCO31,ACO31,ACO3Q,AHY1*ACO31/AHCO31
C    3,DPHCO,RPCACX,CCO31,XMINN,XMINP,TSLX,CCO31,A2(L,NY,NX) 
C     ENDIF
C
C     RCO2Q=CO2-HCO3+H dissociation
C
      XMINN=FIONX*ACO21
      XMINP=FIONX*AMIN1(AHY1,AHCO31)
      AHCO3Q=DPCO2*ACO21/AHY1
      RCO2Q=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AHCO31-AHCO3Q)))
C    2,(AMIN1(TSLX,XMINP,ACO31/AHCO31*(AHCO31-AHCO3Q))))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,23)'RCO2Q',I,J,NFZ,NX,NY,L,M
C    2,RCO2Q,CHY1,CHYX,AHY1,ACO21,AHCO31,AHCO3Q 
C    3,DPCO2,CHCO31,RHCO3,XMINN,XMINP,TSLX 
C     ENDIF
C
C     RALO1=ALOH-AL+OH dissociation
C
      XMINN=FIONX*AALO1
      XMINP=FIONX*AMIN1(AAL1,AOH1)
      AAL1Q=DPAL1*AALO1/AOH1
      RALO1=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AAL1-AAL1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24 
C    2.AND.L.EQ.2)THEN
C     IF(L.EQ.2)THEN
C     WRITE(*,23)'RALO1',I,J,NFZ,NX,NY,L,M
C    2,RALO1,TSLX,AAL1,AAL1Q,DPAL1,AALO1,AOH1,XMINN,XMINP 
C     ENDIF
C
C     RALO2=ALOH2-ALOH+OH dissociation 
C
      XMINN=FIONX*AALO2
      XMINP=FIONX*AMIN1(AALO1,AOH1)
      AALO1Q=DPAL2*AALO2/AOH1
      RALO2=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AALO1-AALO1Q)))
C
C     RALO3=ALOH3-ALOH2+OH dissociation 
C
      XMINN=FIONX*AALO3
      XMINP=FIONX*AMIN1(AALO2,AOH1)
      AALO2Q=DPAL3*AALO3/AOH1
      RALO3=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AALO2-AALO2Q)))
C
C     RALO4=ALOH4-ALOH3+OH dissociation 
C
      XMINN=FIONX*AALO4
      XMINP=FIONX*AMIN1(AALO3,AOH1)
      AALO3Q=DPAL4*AALO4/AOH1
      RALO4=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AALO3-AALO3Q)))
C     IF(L.EQ.11)THEN
C     WRITE(*,23)'RALO4',I,J,NFZ,NX,NY,L,M
C    2,RALO4,TSLX,AALO3,AOH1,DPAL4,AALO4,AALO3*AOH1-DPAL4*AALO4
C    3,COH1,A1(L,NY,NX),ZOH(L,NY,NX),VOLW(L,NY,NX) 
C     ENDIF
C
C     RALS=ALSO4-AL+SO4 dissociation 
C
      XMINN=FIONX*AALS1
      XMINP=FIONX*AMIN1(AAL1,ASO41)
      AAL1Q=DPALS*AALS1/ASO41
      RALS=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AAL1-AAL1Q))) 
C     WRITE(*,23)'RALS',I,J,NFZ,NX,NY,L,M
C    2,RALS,TSLX,AAL1,AAL1Q,ASO41,DPALS,AALS1 
C
C     RFEO1=FEOH-FE+OH dissociation 
C
      XMINN=FIONX*AFEO1
      XMINP=FIONX*AMIN1(AFE1,AOH1)
      AFE1Q=DPFE1*AFEO1/AOH1
      RFEO1=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AFE1-AFE1Q)))
C
C     RFEO2=FEOH2-FEOH+OH dissociation 
C
      XMINN=FIONX*AFEO2
      XMINP=FIONX*AMIN1(AFEO1,AOH1)
      AFEO1Q=DPFE2*AFEO2/AOH1
      RFEO2=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AFEO1-AFEO1Q)))
C
C     RFEO3=FEOH3-FEOH2+OH dissociation 
C
      XMINN=FIONX*AFEO3
      XMINP=FIONX*AMIN1(AFEO2,AOH1)
      AFEO2Q=DPFE3*AFEO3/AOH1
      RFEO3=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AFEO2-AFEO2Q)))
C
C     RFEO4=ALOH4-ALOH3+OH dissociation 
C
      XMINN=FIONX*AFEO4
      XMINP=FIONX*AMIN1(AFEO3,AOH1)
      AFEO3Q=DPFE4*AFEO4/AOH1
      RFEO4=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AFEO3-AFEO3Q)))
C
C     RFES-FE+SO4 dissociation 
C
      XMINN=FIONX*AFES1
      XMINP=FIONX*AMIN1(AFE1,ASO41)
      AFE1Q=DPFES*AFES1/ASO41
      RFES=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AFE1-AFE1Q))) 
C
C     RCAO=CAOH-CA+OH dissociation 
C
      XMINN=FIONX*ACAO1
      XMINP=FIONX*AMIN1(ACA1,AOH1)
      ACA1Q=DPCAO*ACAO1/AOH1
      RCAO=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ACA1-ACA1Q)))
C
C     RCAC=CACO3-CA+CO3 dissociation 
C
      XMINN=FIONX*ACAC1
      XMINP=FIONX*AMIN1(ACA1,ACO31)
      ACA1Q=DPCAC*ACAC1/ACO31
      RCAC=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ACA1-ACA1Q))) 
C
C     RCAH=CAHCO3-CA+HCO3 dissociation 
C
      XMINN=FIONX*ACAH1
      XMINP=FIONX*AMIN1(ACA1,AHCO31)
      ACA1Q=DPCAH*ACAH1/AHCO31
      RCAH=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ACA1-ACA1Q))) 
C
C     RCAS=CASO4-CA+SO4 dissociation 
C
      XMINN=FIONX*ACAS1
      XMINP=FIONX*AMIN1(ACA1,ASO41)
      ACA1Q=DPCAS*ACAS1/ASO41
      RCAS=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ACA1-ACA1Q))) 
C
C     RMGO=MGOH-MG+OH dissociation 
C
      XMINN=FIONX*AMGO1
      XMINP=FIONX*AMIN1(AMG1,AOH1)
      AMG1Q=DPMGO*AMGO1/AOH1
      RMGO=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AMG1-AMG1Q)))
C
C     RMGC=MGCO3-MG+CO3 dissociation 
C
      XMINN=FIONX*AMGC1
      XMINP=FIONX*AMIN1(AMG1,ACO31)
      AMG1Q=DPMGC*AMGC1/ACO31
      RMGC=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AMG1-AMG1Q)))
C
C     RMGH=MGHCO3-MG+HCO3 dissociation 
C
      XMINN=FIONX*AMGH1
      XMINP=FIONX*AMIN1(AMG1,AHCO31)
      AMG1Q=DPMGH*AMGH1/AHCO31
      RMGH=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AMG1-AMG1Q)))
C
C     RMGS=MGSO4-MG+SO4 dissociation 
C
      XMINN=FIONX*AMGS1
      XMINP=FIONX*AMIN1(AMG1,ASO41)
      AMG1Q=DPMGS*AMGS1/ASO41
      RMGS=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AMG1-AMG1Q)))
C
C     RNAC=NACO3-NA+CO3 dissociation 
C
      XMINN=FIONX*ANAC1
      XMINP=FIONX*AMIN1(ANA1,ACO31)
      ANA1Q=DPNAC*ANAC1/ACO31
      RNAC=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ANA1-ANA1Q)))
C
C     RNAS=NASO4-NA+SO4 dissociation 
C
      XMINN=FIONX*ANAS1
      XMINP=FIONX*AMIN1(ANA1,ASO41)
      ANA1Q=DPNAS*ANAS1/ASO41
      RNAS=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ANA1-ANA1Q)))
C
C     RKAS=KSO4-K+SO4 dissociation 
C
      XMINN=FIONX*AKAS1
      XMINP=FIONX*AMIN1(AKA1,ASO41)
      AKA1Q=DPKAS*AKAS1/ASO41
      RKAS=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AKA1-AKA1Q)))
C
C     PHOSPHORUS IN NON-BAND SOIL ZONE
C
      IF(VOLWPO.GT.ZEROS2(NY,NX))THEN
C
C     RH1P=HPO4-H+PO4 dissociation in non-band 
C
      XMINN=FIONX*AH1P1
      XMINP=FIONX*AMIN1(AHY1,AH0P1)
      AH0P1Q=DPH1P*AH1P1/AHY1
      RH1P=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AH0P1-AH0P1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     WRITE(*,22)'RH1P',I,J,NFZ,NX,NY,L,M
C    2,RH1P,AH2P1,AH1P1,AH0P1,AH0P1Q,AHY1*AH0P1/AH1P1,DPH1P
C    3,AHY1,TSLH,TSLX 
C     ENDIF
C
C     RH2P=H2PO4-H+HPO4 dissociation in non-band 
C
      XMINN=FIONX*AH2P1
      XMINP=FIONX*AMIN1(AHY1,AH1P1)
      AH1P1Q=DPH2P*AH2P1/AHY1
      RH2P=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AH1P1-AH1P1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     IF(L.EQ.2)THEN
C     WRITE(*,22)'RH2P',I,J,NFZ,NX,NY,L,M
C    2,RH2P,AH2P1,AH1P1,AH1P1Q,AHY1*AH1P1/AH2P1,DPH2P
C    3,AHY1,TSLH,TSLX 
22    FORMAT(A8,7I4,60E12.4)
C     ENDIF
C
C     RH3P=H3PO4-H+H2PO4 dissociation in non-band 
C
      XMINN=FIONX*AH3P1
      XMINP=FIONX*AMIN1(AHY1,AH2P1)
      AH2P1Q=DPH3P*AH3P1/AHY1
      RH3P=AMAX1(-TSLX,-XMINN 
     2,(AMIN1(TSLX,XMINP,AH2P1-AH2P1Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.L.EQ.2)THEN
C     WRITE(*,22)'RH3P',I,J,NFZ,NX,NY,L,M
C    2,RH3P,AH3P1,AH2P1,AH2P1Q,AHY1*AH2P1/AH3P1
C    3,DPH3P,AHY1,TSLH,TSLX 
C     ENDIF
C
C     RF1P=FEHPO4-FE+HPO4 dissociation in non-band 
C
      XMINN=FIONX*AF1P1
      XMINP=FIONX*AMIN1(AH1P1,AFE1)
      AFE1Q=DPF1P*AF1P1/AH1P1
      RF1P=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AFE1-AFE1Q)))
C
C     RF2P=FEH2PO4-FE+H2PO4 dissociation in non-band 
C
      XMINN=FIONX*AF2P1
      XMINP=FIONX*AMIN1(AH2P1,AFE1)
      AFE1Q=DPF2P*AF2P1/AH2P1
      RF2P=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AFE1-AFE1Q)))
C
C     RC0P=CAPO4-CA+PO4 dissociation in non-band 
C
      RC0P=0.0
C
C     RC1P=CAHPO4-CA+HPO4 dissociation in non-band 
C
      XMINN=FIONX*AC1P1
      XMINP=FIONX*AMIN1(AH1P1,ACA1)
      ACA1Q=DPC1P*AC1P1/AH1P1
      RC1P=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ACA1-ACA1Q)))
C
C     RC2P=CAH2PO4-CA+H2PO4 dissociation in non-band 
C
      XMINN=FIONX*AC2P1
      XMINP=FIONX*AMIN1(AH2P1,ACA1)
      ACA1Q=DPC2P*AC2P1/AH2P1
      RC2P=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ACA1-ACA1Q)))
C
C     RM1P=MGHPO4-MG+HPO4 dissociation in non-band 
C
      XMINN=FIONX*AM1P1
      XMINP=FIONX*AMIN1(AH1P1,AMG1)
      AMG1Q=DPM1P*AM1P1/AH1P1
      RM1P=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AMG1-AMG1Q)))
      ELSE
      RH1P=0.0
      RH2P=0.0
      RH3P=0.0
      RF1P=0.0
      RF2P=0.0
      RC0P=0.0
      RC1P=0.0
      RC2P=0.0
      RM1P=0.0
      ENDIF
C
C     PHOSPHORUS IN BAND SOIL ZONE
C
      IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
C
C     RH1B=HPO4-H+PO4 dissociation in band
C
      XMINN=FIONX*AH1PB
      XMINP=FIONX*AMIN1(AHY1,AH0PB)
      AH0PBQ=DPH1P*AH1PB/AHY1
      RH1B=AMAX1(-TSLX,-XMINN 
     2,(AMIN1(TSLX,XMINP,AH0PB-AH0PBQ)))
C
C     RH2B=H2PO4-H+HPO4 dissociation in band 
C
      XMINN=FIONX*AH2PB
      XMINP=FIONX*AMIN1(AHY1,AH1PB)
      AH1PBQ=DPH2P*AH2PB/AHY1
      RH2B=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AH1PB-AH1PBQ)))
C
C     RH3B=H3PO4-H+H2PO4 dissociation in band 
C
      XMINN=FIONX*AH3PB
      XMINP=FIONX*AMIN1(AHY1,AH2PB)
      AH2PBQ=DPH3P*AH3PB/AHY1
      RH3B=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AH2PB-AH2PBQ)))
C
C     RF1B=FEHPO4-FE+HPO4 dissociation in band 
C
      XMINN=FIONX*AF1PB
      XMINP=FIONX*AMIN1(AH1PB,AFE1)
      AFEBQ=DPF1P*AF1PB/AH1PB
      RF1B=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AFE1-AFEBQ)))
C
C     RF2B=FEH2PO4-FE+H2PO4 dissociation in band 
C
      XMINN=FIONX*AF2PB
      XMINP=FIONX*AMIN1(AH2PB,AFE1)
      AFEBQ=DPF2P*AF2PB/AH2PB
      RF2B=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AFE1-AFEBQ)))
C
C     RC0B=CAPO4-CA+PO4 dissociation in band 
C
      RC0B=0.0
C
C     RC1B=CAHPO4-CA+HPO4 dissociation in band 
C
      XMINN=FIONX*AC1PB
      XMINP=FIONX*AMIN1(AH1PB,ACA1)
      ACABQ=DPC1P*AC1PB/AH1PB
      RC1B=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ACA1-ACABQ)))
C
C     RC2B=CAH2PO4-CA+H2PO4 dissociation in band 
C
      XMINN=FIONX*AC2PB
      XMINP=FIONX*AMIN1(AH2PB,ACA1)
      ACABQ=DPC2P*AC2PB/AH2PB
      RC2B=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,ACA1-ACABQ)))
C
C     RM1B=MGHPO4-MG+HPO4 dissociation in band 
C
      XMINN=FIONX*AM1PB
      XMINP=FIONX*AMIN1(AH1PB,AMG1)
      AMGBQ=DPM1P*AM1PB/AH1PB
      RM1B=AMAX1(-TSLX,-XMINN
     2,(AMIN1(TSLX,XMINP,AMG1-AMGBQ)))
      ELSE
      RH1B=0.0
      RH2B=0.0
      RH3B=0.0
      RF1B=0.0
      RF2B=0.0
      RC0B=0.0
      RC1B=0.0
      RC2B=0.0
      RM1B=0.0
      ENDIF
C
C     SILICATE ROCK WEATHERING
C
C     Q*SI=silicate (mol)
C     A*1=ion activity (mol m-3)
C     A*Q=ion activity in equilibrium with Q*SI at SP* (mol m-3)
C     SP*SI=rock solubility product
C     CHYSI=soluble H silicate (mol m-3)
C     RQ*SI=silicate weathering (mol m-3 t-1)
C     TRWXA=maximum rate for silicate weathering (mol m-3 t-1) 	
C     salt code
C        :AL=Al,FE=Fe,CA=Ca,MG=Mg,NA=Na,KA=K 
C
      AALQ=SPALSI*AHY1**3/CHYSI1**0.75
      RQALSI=AMIN1(0.0,AMAX1(-QALSI1,-TRWXA,(AMIN1(TRWXA
     2,(AAL1-AALQ)))))
      AFEQ=SPFESI*AHY1**3/CHYSI1**0.75
      RQFESI=AMIN1(0.0,AMAX1(-QFESI1,-TRWXA,(AMIN1(TRWXA
     2,(AFE1-AFEQ)))))
      ACAQ=SPCASI*AHY1**2/CHYSI1**0.50
      RQCASI=AMIN1(0.0,AMAX1(-QCASI1,-TRWXA,(AMIN1(TRWXA
     2,(ACA1-ACAQ)))))
      AMGQ=SPMGSI*AHY1**2/CHYSI1**0.50
      RQMGSI=AMIN1(0.0,AMAX1(-QMGSI1,-TRWXA,(AMIN1(TRWXA
     2,(AMG1-AMGQ)))))
      ANAQ=SPNASI*AHY1/CHYSI1**0.25
      RQNASI=AMIN1(0.0,AMAX1(-QNASI1,-TRWXA,(AMIN1(TRWXA
     2,(ANA1-ANAQ)))))
      AKAQ=SPKASI*AHY1/CHYSI1**0.25
      RQKASI=AMIN1(0.0,AMAX1(-QKASI1,-TRWXA,(AMIN1(TRWXA
     2,(AKA1-AKAQ)))))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.M.EQ.MRXN)THEN
C     WRITE(*,1111)'RQALS',I,J,NFZ,NX,NY,L,M
C    2,RQALSI,QALSI1,TRWXA,AAL1,AALQ,SPALSI,AHY1,CHY1,CHYSI1
C    3,ZHYSI(L,NY,NX),SSA(L,NY,NX) 
C     WRITE(*,1111)'RQFES',I,J,NFZ,NX,NY,L,M
C    2,RQFESI,QFESI1,TRWXA,AFE1,AFEQ,SPFESI,AHY1,CHY1,CHYSI1
C    3,ZHYSI(L,NY,NX),SSA(L,NY,NX) 
C     WRITE(*,1111)'RQCAS',I,J,NFZ,NX,NY,L,M
C    2,RQCASI,QCASI1,TRWXA,ACA1,ACAQ,SPCASI,AHY1,CHY1,CHYSI1
C    3,ZHYSI(L,NY,NX),SSA(L,NY,NX)
C     WRITE(*,1111)'RQMGS',I,J,NFZ,NX,NY,L,M
C    2,RQMGSI,QMGSI1,TRWXA,AMG1,AMGQ,SPMGSI,AHY1,CHY1,CHYSI1
C    3,ZHYSI(L,NY,NX),SSA(L,NY,NX) 
C     ENDIF 
C
C     TOTAL ION TRANSFORMATIONS FOR CURRENT ITERATION
C     FROM ALL REACTIONS ABOVE
C     (all transformations in mol m-3 t-1)
C
C     RN4S,RN4B=net NH4 transformation in non-band,band
C     RN3S,RN3B=net NH3 transformation in non-band,band
C     RAL,RFE,RHY,RCA,RMG,RNA,RKA,ROH=net Al,Fe,H,Ca,Mg,Na,K,OH
C        transformation
C     RSO4,RCO3,RHCO,RCO2=net SO4,CO3,HCO3,CO2 transformation
C     RAL1,RAL2,RAL3,RAL4,RALS=net AlOH,AlOH2,AlOH3,AlOH4,AlSO4
C     RFE1,RFE2,RFE3,RFE4,RFES=net FeOH,FeOH2,FeOH3,FeOH4,FeSO4
C     RHP0,RHP1,RHP2,RHP3=net PO4,HPO4,H2PO4,H3PO4 transformation 
C        in non-band
C     RXH0,RXH1,RXH2,RX1P,RX2P=net R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 
C        in non-band
C     RHB0,RHB1,RHB2,RHB3=net PO4,HPO4,H2PO4,H3PO4 transformation 
C        in band
C     RBH0,RBH1,RBH2,RB1P,RB2P=net R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 
C        in band
C     RH2O=net change soil in water content 
C
      RN4S=RNH4-RXN4*BKVLWH
      RN4B=RNHB-RXNB*BKVLWB 
      RN3S=-RNH4
      RN3B=-RNHB
      RAL=-RQALSI-RXAL*BKVLW-RALO1-RALS
     2-RPALPX*VLPO4(L,NY,NX)-RPALBX*VLPOB(L,NY,NX)
      RFE=-RQFESI-RXFE*BKVLW-RFEO1-RFES
     2-(RPFEPX+RF1P+RF2P)*VLPO4(L,NY,NX)
     2-(RPFEBX+RF1B+RF2B)*VLPOB(L,NY,NX)
      RHY=RZHYS+3.0*(RQALSI+RQFESI)+2.0*(RQCASI+RQMGSI)
     2+RQNASI+RQKASI-RXHY*BKVLW-RXHC*BKVLW-RCO2Q-RHCO3
     3-RNH4*VLNH4(L,NY,NX)-RNHB*VLNHB(L,NY,NX)
     4+2.0*(RPALPX+RPFEPX)*VLPO4(L,NY,NX)
     5+2.0*(RPALBX+RPFEBX)*VLPOB(L,NY,NX)
     6+6.0*(RPCAHX*VLPO4(L,NY,NX)+RPCHBX*VLPOB(L,NY,NX))
     7+(RPCADX-RXOH2-RXOH1-RH1P-RH2P-RH3P)*VLPO4(L,NY,NX)
     8+(RPCDBX-RXO2B-RXO1B-RH1B-RH2B-RH3B)*VLPOB(L,NY,NX)
      ROH=(-RPCAHX+RYH2P+RXH1P)*VLPO4(L,NY,NX)
     3+(-RPCHBX+RYH2B+RXH1B)*VLPOB(L,NY,NX)
     3-RCAO-RMGO-RALO1-RALO2-RALO3-RALO4-RFEO1-RFEO2-RFEO3-RFEO4
      RCA=-RQCASI-RXCA*BKVLW-RCAO-RCAC-RCAH-RCAS
     2-(RPCADX+RPCAMX+RC0P+RC1P+RC2P)*VLPO4(L,NY,NX)
     3-(RPCDBX+RPCMBX+RC0B+RC1B+RC2B)*VLPOB(L,NY,NX)
     4-5.0*(RPCAHX*VLPO4(L,NY,NX)+RPCHBX*VLPOB(L,NY,NX))
      RMG=-RQMGSI-RXMG*BKVLW-RMGO-RMGC-RMGH-RMGS
     2-RM1P*VLPO4(L,NY,NX)-RM1B*VLPOB(L,NY,NX)
      RNA=-RQNASI-RXNA*BKVLW-RNAC-RNAS
      RKA=-RQKASI-RXKA*BKVLW-RKAS
      RSO4=-RALS-RFES-RCAS-RMGS-RNAS-RKAS
      RCO3=-RHCO3-RCAC-RMGC-RNAC
      RHCO=-RCO2Q-RCAH-RMGH+RHCO3
      RCO2=RCO2Q
      RAL1=RALO1-RALO2
      RAL2=RALO2-RALO3
      RAL3=RALO3-RALO4
      RAL4=RALO4
      RFE1=RFEO1-RFEO2
      RFE2=RFEO2-RFEO3
      RFE3=RFEO3-RFEO4
      RFE4=RFEO4
      RHP0=-RH1P-RC0P
      RHP1=-RXH1P+RH1P-RH2P-RF1P-RC1P-RM1P
      RHP2=-RPALPX-RPFEPX-RPCADX-3.0*RPCAHX
     4-2.0*RPCAMX-RXH2P-RYH2P+RH2P-RH3P-RF2P-RC2P
      RHP3=RH3P
      RXH0=-RXOH1
      RXH1=RXOH1-RXOH2-RYH2P-RXH1P 
      RXH2=RXOH2-RXH2P
      RX1P=RXH1P 
      RX2P=RXH2P+RYH2P 
      RHB0=-RH1B-RC0B
      RHB1=-RXH1B+RH1B-RH2B-RF1B-RC1B-RM1B
      RHB2=-RPALBX-RPFEBX-RPCDBX-3.0*RPCHBX
     4-2.0*RPCMBX-RXH2B-RYH2B+RH2B-RH3B-RF2B-RC2B
      RHB3=RH3B
      RBH0=-RXO1B
      RBH1=RXO1B-RXO2B-RYH2B-RXH1B 
      RBH2=RXO2B-RXH2B
      RB1P=RXH1B 
      RB2P=RXH2B+RYH2B 
      RH2O=RCO2Q 
     2+RXH2P*VLPO4(L,NY,NX)+RXH2B*VLPOB(L,NY,NX)
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     IF(L.EQ.1)THEN
C     WRITE(*,23)'RHYL',I,J,NFZ,NX,NY,L,M
C    2,CHY1,CHYX,AHY1,RHY,RZHYS,3.0*(RQALSI+RQFESI),2.0*(RQCASI+RQMGSI)
C    2,RQNASI+RQKASI,-RXHY*BKVLW,-RXHC*BKVLW,-RCO2Q,-RHCO3
C    2,-RNH4,VLNH4(L,NY,NX),-RNHB,VLNHB(L,NY,NX)
C    4,2.0*RPALPX*VLPO4(L,NY,NX),2.0*RPALBX*VLPOB(L,NY,NX)
C    6,2.0*RPFEPX*VLPO4(L,NY,NX),2.0*RPFEBX*VLPOB(L,NY,NX)
C    9,6.0*RPCAHX*VLPO4(L,NY,NX),6.0*RPCHBX*VLPOB(L,NY,NX)
C    6,RPCADX,-RXOH2,-RXOH1,-RH1P,-RH2P,-RH3P,VLPO4(L,NY,NX)
C    7,RPCDBX,-RXO2B,-RXO1B,-RH1B,-RH2B,-RH3B,VLPOB(L,NY,NX)
C    8,ZHY(L,NY,NX),PH(L,NY,NX)
C     WRITE(*,23)'ROHL',I,J,NFZ,NX,NY,L,M
C    7,COH1,AOH1,ROH,RHHX,-3.0*RPALOX,-3.0*RPFEOX 
C    3,-RPCAHX,RYH2P,RXH1P,VLPO4(L,NY,NX)
C    5,-RPCHBX,RYH2B,RXH1B,VLPOB(L,NY,NX)
C    2,-RCAO,-RMGO 
C    2,-RALO1,-RALO2,-RALO3,-RALO4
C    4,-RFEO1,-RFEO2,-RFEO3,-RFEO4
C    8,ZOH(L,NY,NX)/VOLW(L,NY,NX),ZOH(L,NY,NX),TROH(L,NY,NX)
C     WRITE(*,23)'RN4L',I,J,NFZ,NX,NY,L,M
C    2,RN4S,RNH4,-RXN4*BKVLWH,-RXN4,CN41,RN4X 
C     WRITE(*,23)'RALL',I,J,NFZ,NX,NY,L,M
C    2,CAL1,AAL1,RAL,-RPALOX,-RQALSI,-RXAL*BKVLW,-RALO1,-RALS
C    2,-RPALPX*VLPO4(L,NY,NX),-RPALBX*VLPOB(L,NY,NX)
C     WRITE(*,23)'RFEL',I,J,NFZ,NX,NY,L,M
C    2,CFE1,AFE1,RFE,-RPFEOX,-RQFESI,-RXFE*BKVLW,-RFEO1,-RFES
C    2,-RPFEPX*VLPO4(L,NY,NX),-RF1P*VLPO4(L,NY,NX),-RF2P*VLPO4(L,NY,NX)
C    3,-RPFEBX*VLPOB(L,NY,NX),-RF1B*VLPOB(L,NY,NX),-RF2B*VLPOB(L,NY,NX)
C     WRITE(*,23)'RCAL',I,J,NFZ,NX,NY,L,M,CCA1,ACA1,RCA
C    2,-RPCACX,-RQCASI,-RPCASO,-RXCA*BKVLW,-RCAO,-RCAC,-RCAH,-RCAS
C    2,-RPCADX,-RPCAMX,-RC0P,-RC1P,-RC2P,VLPO4(L,NY,NX)
C    3,-RPCDBX,-RPCMBX,-RC0B,-RC1B,-RC2B,VLPOB(L,NY,NX)
C    4,-5.0*RPCAHX*VLPO4(L,NY,NX),-5.0*RPCHBX*VLPOB(L,NY,NX)
C     WRITE(*,23)'RKA',I,J,NFZ,NX,NY,L,M
C    2,RKA,RXKA,RQKASI,RKAS
C     WRITE(*,23)'RHP2',I,J,NFZ,NX,NY,L,M 
C    2,CH2P1,AH2P1,RHP2,-RPALPX,-RPFEPX,-RPCADX,-3.0*RPCAHX
C    4,-2.0*RPCAMX,-RXH2P,-RYH2P,+RH2P,-RH3P,-RF2P,-RC2P,RH2PX
C    5,PH(L,NY,NX)
C     WRITE(*,23)'RHP1',I,J,NFZ,NX,NY,L,M 
C    3,CH1P1,AH1P1,RHP1,-RXH1P,RH1P,-RH2P,-RF1P,-RC1P,-RM1P,RH1PX
C     WRITE(*,23)'RHP0',I,J,NFZ,NX,NY,L,M
C    2,CH0P1,AH0P1,RHP0,-RH1P,-RC0P 
C     WRITE(*,23)'RHB2',I,J,NFZ,NX,NY,L,M
C    2,CH2PB,AH2PB,RHB2,-RPALBX,-RPFEBX,-RPCDBX,-3.0*RPCHBX
C    4,-2.0*RPCMBX,-RXH2B,-RYH2B,+RH2B,-RH3B,-RF2B,-RC2B
C    3,VLPOB(L,NY,NX)
C     WRITE(*,23)'X21PB',I,J,NFZ,NX,NY,L,M 
C    2,X21PB,RBH2,RXO2B,-RXH2B 
C     WRITE(*,23)'X1P1B',I,J,NFZ,NX,NY,L,M 
C    2,X1P1B,RBH1 RXO1B,-RXO2B,-RYH2B,-RXH1B 
C     WRITE(*,23)'RHCO',I,J,NFZ,NX,NY,L,M,RHCO
C    2,RCO2Q,RCAH,RMGH,RHCO3
C     WRITE(*,23)'RCO3',I,J,NFZ,NX,NY,L,M
C    2,RCO3,-RPCACX,-RHCO3,-RCAC,-RMGC,-RNAC
C    3,RHCO,-RCO2Q,-RCAH,-RMGH,RHCO3
C    4,RCO2,RCO2Q
C    2,CCO31,CHCO31,CCO21,DPHCO,DPCO2
C    5,AHY1,AHCO31,ACO21,TRCO2(L,NY,NX) 
C     WRITE(*,23)'RSO4',I,J,NFZ,NX,NY,L,M
C    2,CSO41,RSO4,RPCASO,RALS,RFES,RCAS,RMGS,RNAS,RKAS
23    FORMAT(A8,7I4,100E12.4)
C     ENDIF
C
C     UPDATE ION CONCENTRATIONS FOR NEXT ITERATION
C     FROM TOTAL ION TRANSFORMATIONS ABOVE
C
C     C*1=ion concentration (mol m-3)
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+
C          :*CA*=Ca2+,*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-
C          :*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
C          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
C          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
C          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
C          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH
C          :*MGC*=MgCO3,*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-
C          :*NAS*=NaSO4-,*KAS*=KSO4-
C     cation code:AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C          :N4=NH4+ non-band,NB=NH4+ band,HY=H+
C     anion code:OH=OH-,*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C          :*F2P*=F1H2PO4-,*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C          :*M1P*=MgHPO4,*COO*=COOH-
C          :*1=non-band,*B=band
C
C     UPDATE EXCHANGEABLE ION CONCENTRATIONS FOR NEXT
C     ITERATION FROM TOTAL ION TRANSFORMATIONS
C
C     C*1=soluble ion concentrations (mol m-3)
C     X*1=exchangeable ion concentrations (mol Mg-1)
C     R*=ion transformation rate (mol m-3 t-1)
C     RX*=ion adsorption  rate (mol m-3 t-1)
C     cation code:AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C          :N4=NH4+ non-band,NB=NH4+ band,HY=H+
C     anion code:OH=OH-,*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C          :*F2P*=F1H2PO4-,*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C          :*M1P*=MgHPO4,*COO*=COOH-
C          :*1=non-band,*B=band
C     silicate code:HYSI=H silicate
C
      CHY1=CHY1+RHY
      COH1=COH1+ROH
      CAL1=AMAX1(ZEROC,CAL1+RAL)
      CFE1=AMAX1(ZEROC,CFE1+RFE)
      IF(M.NE.MRXN)THEN
      CN41=AMAX1(ZEROC,CN41+RN4S+RN4X)
      CN4B=AMAX1(ZEROC,CN4B+RN4B+RNBX)
      CN31=AMAX1(ZEROC,CN31+RN3S+RN3X)
      CN3B=AMAX1(ZEROC,CN3B+RN3B+R3BX)
      CCA1=AMAX1(ZEROC,CCA1+RCA)
      CMG1=AMAX1(ZEROC,CMG1+RMG)
      CNA1=AMAX1(ZEROC,CNA1+RNA)
      CKA1=AMAX1(ZEROC,CKA1+RKA)
      CSO41=AMAX1(ZEROC,CSO41+RSO4)
      CCO31=AMAX1(ZEROC,CCO31+RCO3) 
      CHCO31=AMAX1(ZEROC,CHCO31+RHCO)
      CCO21=AMAX1(ZEROC,CCO21+RCO2)
      CALO1=AMAX1(ZEROC,CALO1+RAL1)
      CALO2=AMAX1(ZEROC,CALO2+RAL2) 
      CALO3=AMAX1(ZEROC,CALO3+RAL3)
      CALO4=AMAX1(ZEROC,CALO4+RAL4)
      CALS1=AMAX1(ZEROC,CALS1+RALS)
      CFEO1=AMAX1(ZEROC,CFEO1+RFE1)
      CFEO2=AMAX1(ZEROC,CFEO2+RFE2) 
      CFEO3=AMAX1(ZEROC,CFEO3+RFE3)
      CFEO4=AMAX1(ZEROC,CFEO4+RFE4)
      CFES1=AMAX1(ZEROC,CFES1+RFES)
      CCAO1=AMAX1(ZEROC,CCAO1+RCAO)
      CCAC1=AMAX1(ZEROC,CCAC1+RCAC)
      CCAH1=AMAX1(ZEROC,CCAH1+RCAH)
      CCAS1=AMAX1(ZEROC,CCAS1+RCAS)
      CMGO1=AMAX1(ZEROC,CMGO1+RMGO)
      CMGC1=AMAX1(ZEROC,CMGC1+RMGC)
      CMGH1=AMAX1(ZEROC,CMGH1+RMGH)
      CMGS1=AMAX1(ZEROC,CMGS1+RMGS)
      CNAC1=AMAX1(ZEROC,CNAC1+RNAC)
      CNAS1=AMAX1(ZEROC,CNAS1+RNAS)
      CKAS1=AMAX1(ZEROC,CKAS1+RKAS)
      CH0P1=AMAX1(ZEROC,CH0P1+RHP0)
      CH1P1=AMAX1(ZEROC,CH1P1+RHP1+RH1PX)
      CH2P1=AMAX1(ZEROC,CH2P1+RHP2+RH2PX)
      CH3P1=AMAX1(ZEROC,CH3P1+RHP3)
      CF1P1=AMAX1(ZEROC,CF1P1+RF1P)
      CF2P1=AMAX1(ZEROC,CF2P1+RF2P)
      CC0P1=AMAX1(ZEROC,CC0P1+RC0P)
      CC1P1=AMAX1(ZEROC,CC1P1+RC1P)
      CC2P1=AMAX1(ZEROC,CC2P1+RC2P)
      CM1P1=AMAX1(ZEROC,CM1P1+RM1P)
      CH0PB=AMAX1(ZEROC,CH0PB+RHB0)
      CH1PB=AMAX1(ZEROC,CH1PB+RHB1+RH1BX)
      CH2PB=AMAX1(ZEROC,CH2PB+RHB2+RH2BX)
      CH3PB=AMAX1(ZEROC,CH3PB+RHB3)
      CF1PB=AMAX1(ZEROC,CF1PB+RF1B)
      CF2PB=AMAX1(ZEROC,CF2PB+RF2B)
      CC0PB=AMAX1(ZEROC,CC0PB+RC0B)
      CC1PB=AMAX1(ZEROC,CC1PB+RC1B)
      CC2PB=AMAX1(ZEROC,CC2PB+RC2B)
      CM1PB=AMAX1(ZEROC,CM1PB+RM1B)
      CHYSI1=AMAX1(ZEROC,CHYSI1-0.75*(RQALSI+RQFESI)
     2-0.50*(RQCASI+RQMGSI)-0.25*(RQNASI+RQKASI))
      XN41=XN41+RXN4
      XNB1=XNB1+RXNB
      XHY1=XHY1+RXHY
      XAL1=XAL1+RXAL
      XFE1=XFE1+RXFE
      XCA1=XCA1+RXCA
      XMG1=XMG1+RXMG
      XNA1=XNA1+RXNA
      XKA1=XKA1+RXKA
      XHC1=XHC1+RXHC
      XOH01=XOH01+RXH0 
      XOH11=XOH11+RXH1 
      XOH21=XOH21+RXH2 
      XH1P1=XH1P1+RX1P 
      XH2P1=XH2P1+RX2P 
      XH01B=XH01B+RBH0 
      XH11B=XH11B+RBH1 
      XH21B=XH21B+RBH2 
      X1P1B=X1P1B+RB1P 
      X2P1B=X2P1B+RB2P
C
C     UPDATE PRECIPITATE CONCENTRATIONS FOR NEXT
C     ITERATION FROM TOTAL ION TRANSFORMATIONS ABOVE
C
C     PALOH1,PFEOH1,PCACO1,PCASO1=precipitated
C         AL(OH)3,FE(OH)3,CACO3,CASO4 vs water (mol m-3)
C     PALPO1,PFEPO1=precipitated AlPO4,FEPO4 in non-band vs water 
C        (mol m-3)
C     PALPOB,PFEPOB=precipitated AlPO4,FEPO4 in band vs water (mol m-3)
C     PCAPM1,PCAPD1,PCAPH1=precipitated CaH2PO4,CaHPO4,apatite 
C        in non-band vs water (mol m-3)
C     PCAPMB,PCAPDB,PCAPHB=precipitated CaH2PO4,CaHPO4,apatite
C        in band vs water (mol m-3)
C     QALSI1,QFESI1,QCASI1,QMGSI1,QNASI1,QKASI1=Al,Fe,Ca,Mg,Na,K
C        silicate vs water (mol m-3)
C     R*X,RQ*I=precipitation(+ve) or dissolution (-ve) rate 
C           (mol m-3 t-1)
C
      QALSI1=QALSI1+RQALSI 
      QFESI1=QFESI1+RQFESI 
      QCASI1=QCASI1+RQCASI 
      QMGSI1=QMGSI1+RQMGSI 
      QNASI1=QNASI1+RQNASI 
      QKASI1=QKASI1+RQKASI 
      PALPO1=PALPO1+RPALPX
      PFEPO1=PFEPO1+RPFEPX 
      PCAPD1=PCAPD1+RPCADX 
      PCAPH1=PCAPH1+RPCAHX 
      PCAPM1=PCAPM1+RPCAMX 
      PALPOB=PALPOB+RPALBX 
      PFEPOB=PFEPOB+RPFEBX 
      PCAPDB=PCAPDB+RPCDBX 
      PCAPHB=PCAPHB+RPCHBX 
      PCAPMB=PCAPMB+RPCMBX 
      ENDIF
C
C     COMPLETE ITERATION WITH PRECIPITATION-DISSOLUTION REACTIONS 
C     FPR AL, FE AND CA THAT DIRECTLY AFFECT PH
C
C     SP*=solubility product from PARAMETER above
C     A*1=ion activity (mol m-3)
C     A*Q=ion equilibrium activity at SP* (mol m-3)
C     A1,A2,A3=activity coefficients from ‘hour1.f’ 
C     RP*X=precipitation(+ve) or dissolution (-ve) rate 
C        (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FION*=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     TPDX=maximum precipitation rate (mol m-3 t-1)
C
C     ALUMINUM HYDROXIDE AL(OH)3 (GIBBSITE)
C
      AAL1=AAL1+RAL 
      XMINN=FIONX*AAL1
      XMINP=FIONH*AAL1
      AAL1Q=SPALO/AOH1**3
      RPALOX=AMAX1(-AMAX1(0.0,PALOH1),-TPDX,-XMINN 
     2,(AMIN1(TPDX,XMINP,AAL1-AAL1Q)))
      RAL=RAL-RPALOX 
      ROH=ROH-3.0*RPALOX
      CAL1=CAL1-RPALOX
      COH1=COH1-3.0*RPALOX
      PALOH1=PALOH1+RPALOX
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,1112)'ALOH3',I,J,NFZ,NX,NY,L,M 
C    2,RPALOX,AAL1,AAL1Q,AOH1,CAL1*A3(L,NY,NX)*AOH1**3
C    3,SPALO,RAL,ROH,PALOH(L,NY,NX),XMINP,TPDX,PH(L,NY,NX)
C     ENDIF 
C
C     IRON HYDROXIDE FE(OH)3
C
      AFE1=AFE1+RFE
      XMINN=FIONX*AFE1
      XMINP=FIONH*AFE1
      AFE1Q=SPFEO/AOH1**3
      RPFEOX=AMAX1(-AMAX1(0.0,PFEOH1),-TPDX,-XMINN 
     2,(AMIN1(TPDX,XMINP,AFE1-AFE1Q)))
      RFE=RFE-RPFEOX 
      ROH=ROH-3.0*RPFEOX
      CFE1=CFE1-RPFEOX
      COH1=COH1-3.0*RPFEOX
      PFEOH1=PFEOH1+RPFEOX
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,1112)'FEOH3',I,J,NFZ,NX,NY,L,M 
C    2,RPFEOX,AFE1,AFE1Q,AOH1,CFE1*A3(L,NY,NX)*AOH1**3
C    3,SPFEO,RFE,ROH,PFEOH(L,NY,NX),XMINP,TPDX,PH(L,NY,NX)
C     ENDIF
C
C     CALCITE CACO3
C
      TPDZO=TPDZ
     2/(1.0+AOH1/OHKI)
      ACA1=ACA1+RCA
      ACO31=ACO31+RCO3
      XMINP=FIONX*AMIN1(ACA1,ACO31)
      ACA1Q=SPCAC/ACO31
      RPCACX=AMAX1(-AMAX1(0.0,PCACO1),-TPDZO
     2,(AMIN1(XMINP,ACA1-ACA1Q)))
      RCA=RCA-RPCACX 
      RCO3=RCO3-RPCACX 
      CCA1=CCA1-RPCACX
      CCO31=CCO31-RPCACX 
      PCACO1=PCACO1+RPCACX
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.(M/5)*5.EQ.M)THEN
C     IF((I/5)*5.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,1112)'CALC',I,J,NFZ,NX,NY,L,M
C    2,RPCACX,ACA1,ACA1Q,CCO31,ACO31,AHCO31,ACO21,AHY1
C    3,RCA,RCO3,ACA1*ACO31,SPCAC
C    4,PCACO(L,NY,NX),ZCA(L,NY,NX)
C    4,XMINP,TPDZ,TPDZO,1.0+AOH1/OHKI,PH(L,NY,NX)
C     ENDIF
C
C     GYPSUM CASO4
C
      XMINP=FIONX*AMIN1(ACA1,ASO41)
      ACA1Q=SPCAS/ASO41
      RPCASO=AMAX1(-AMAX1(0.0,PCASO1),-TPDX
     2,(AMIN1(TPDX,XMINP,ACA1-ACA1Q)))
      RCA=RCA-RPCASO 
      CCA1=CCA1-RPCASO
      RSO4=RSO4-RPCASO
      CSO41=CSO41-RPCASO
      PCASO1=PCASO1+RPCASO
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH.AND.M.EQ.MRXN
C    2.AND.L.EQ.2)THEN
C     WRITE(*,1112)'GYP',I,J,NFZ,NX,NY,L,M,RPCASO
C    2,PCASO1,ACA1,ACA1Q,ASO41,A1(L,NY,NX),A2(L,NY,NX)
C    3,CCA1*CSO41,ACA1*ASO41,SPCAS
C    4,ZCA(L,NY,NX),VOLW(L,NY,NX),CCA1,XMINP,TPDX
C     ENDIF
C
C     ACCUMULATE TOTAL ION TRANSFORMATIONS FROM ALL ITERATIONS
C     FOR USE IN REDIST.F
C        (all units in mol m-3 t-1)
C
C     TRN4S,TRN4B=total NH4 transformation in non-band,band 
C     TRN3S,TRN3B=total NH3 transformation in non-band,band 
C     TRAL,TRFE,TRHY,TRCA,TRMG,TRNA,TRKA,TROH=totalAl,Fe,H,Ca,Mg,
C        Na,K,OH transformation 
C     TRSO4,TRCO3,TRHCO,TRCO2=total SO4,CO3,HCO3,CO2 transformation
C     TRAL1,TRAL2,TRAL3,TRAL4,TRALS=total AlOH,AlOH2,AlOH3,AlOH4,AlSO4
C     TRFE1,TRFE2,TRFE3,TRFE4,TRFES=total FeOH,FeOH2,FeOH3,FeOH4,FeSO4
C     TRCAO,TRCAC,TRCAH,TRCAS=total CaOH,CaCO3,CaHCO3,CaSO4
C        transformation 
C     TRMGO,TRMGC,TRMGH,TRMGS=total MgOH,MgCO3,MgHCO3,MgSO4
C        transformation
C     TRNAC,TRNAS,TRKAS=total NaCO3,NaSO4,KSO4 transformation
C     TRH0P,TRH1P,TRH2P,TRH3P=net PO4,HPO4,H2PO4,H3PO4 transformation
C        in non-band 
C     TRF1P,TRF2P,TRC0P,TRC1P,TRC2P,TRM1P
C        =total FeHPO4,FeH2PO4,CaPO4,CaHPO4,CaH2PO4,MgHPO4 in non-band
C     TRH0B,TRH1B,TRH2B,TRH3B=net PO4,HPO4,H2PO4,H3PO4 transformation
C        in band 
C     TRF1B,TRF2B,TRC0B,TRC1B,TRC2B,TRM1B
C        =total FeHPO4,FeH2PO4,CaPO4,CaHPO4,CaH2PO4,MgHPO4 in band
C     TRXN4,TRXNB=total NH4 adsorption in non-band,band 
C     TRXHY,TRXAL,TRXFE,TRXCA,TRXMG,TRXNA,TRXKA=total
C        H,Al,Fe,Ca,Mg,Na,K adsorption 
C     TRXHC,TRXAL2,TRXFE2=total HCO3,AlOH2,FeOH2 adsorption 
C     TRXH0,TRXH1,TRXH2,TRX1P,TRX2P
C        =total R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in non-band
C     TRBH0,TRBH1,TRBH2,TRB1P,TRB2P
C        =total R-O,R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in band
C     TRALOH,TRFEOH,TRCACO,TRCASO=total AlOH3,FeOH3,CaCO3,CaSO4
C        precipitation 
C     TRALPO,TRFEPO,TRCAPD,TRCAPH,TRCAPM
C        =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation 
C         in non-band 
C     TRALPB,TRFEPB,TRCPDB,TRCPHB,TRCPMB
C        =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation 
C         in band 
C     TRHYSI=total H precipitation from silicate weathering
C     TRALSI,TRFESI,TRCASI,TRMGSI,TRNASI,TRKASI=Al,Fe,Ca,Mg,Na,K
C        silicate weathering
C
      TRN4S(L,NY,NX)=TRN4S(L,NY,NX)+RN4S
      TRN4B(L,NY,NX)=TRN4B(L,NY,NX)+RN4B
      TRN3S(L,NY,NX)=TRN3S(L,NY,NX)+RN3S
      TRN3B(L,NY,NX)=TRN3B(L,NY,NX)+RN3B
      TRHY(L,NY,NX)=TRHY(L,NY,NX)+RHY 
      TROH(L,NY,NX)=TROH(L,NY,NX)+ROH 
      TRAL(L,NY,NX)=TRAL(L,NY,NX)+RAL 
      TRFE(L,NY,NX)=TRFE(L,NY,NX)+RFE 
      TRCA(L,NY,NX)=TRCA(L,NY,NX)+RCA 
      TRMG(L,NY,NX)=TRMG(L,NY,NX)+RMG
      TRNA(L,NY,NX)=TRNA(L,NY,NX)+RNA
      TRKA(L,NY,NX)=TRKA(L,NY,NX)+RKA
      TRSO4(L,NY,NX)=TRSO4(L,NY,NX)+RSO4 
      TRCO3(L,NY,NX)=TRCO3(L,NY,NX)+RCO3 
      TRHCO(L,NY,NX)=TRHCO(L,NY,NX)+RHCO
      TRCO2(L,NY,NX)=TRCO2(L,NY,NX)+RCO2
      TRH2O(L,NY,NX)=TRH2O(L,NY,NX)+RH2O 
      TRAL1(L,NY,NX)=TRAL1(L,NY,NX)+RAL1
      TRAL2(L,NY,NX)=TRAL2(L,NY,NX)+RAL2 
      TRAL3(L,NY,NX)=TRAL3(L,NY,NX)+RAL3
      TRAL4(L,NY,NX)=TRAL4(L,NY,NX)+RAL4
      TRALS(L,NY,NX)=TRALS(L,NY,NX)+RALS
      TRFE1(L,NY,NX)=TRFE1(L,NY,NX)+RFE1
      TRFE2(L,NY,NX)=TRFE2(L,NY,NX)+RFE2
      TRFE3(L,NY,NX)=TRFE3(L,NY,NX)+RFE3
      TRFE4(L,NY,NX)=TRFE4(L,NY,NX)+RFE4
      TRFES(L,NY,NX)=TRFES(L,NY,NX)+RFES
      TRCAO(L,NY,NX)=TRCAO(L,NY,NX)+RCAO
      TRCAC(L,NY,NX)=TRCAC(L,NY,NX)+RCAC
      TRCAH(L,NY,NX)=TRCAH(L,NY,NX)+RCAH
      TRCAS(L,NY,NX)=TRCAS(L,NY,NX)+RCAS
      TRMGO(L,NY,NX)=TRMGO(L,NY,NX)+RMGO
      TRMGC(L,NY,NX)=TRMGC(L,NY,NX)+RMGC
      TRMGH(L,NY,NX)=TRMGH(L,NY,NX)+RMGH
      TRMGS(L,NY,NX)=TRMGS(L,NY,NX)+RMGS
      TRNAC(L,NY,NX)=TRNAC(L,NY,NX)+RNAC
      TRNAS(L,NY,NX)=TRNAS(L,NY,NX)+RNAS
      TRKAS(L,NY,NX)=TRKAS(L,NY,NX)+RKAS
      TRHYSI(L,NY,NX)=TRHYSI(L,NY,NX)-0.75*(RQALSI+RQFESI)
     2-0.50*(RQCASI+RQMGSI)-0.25*(RQNASI+RQKASI)
      TRH0P(L,NY,NX)=TRH0P(L,NY,NX)+RHP0
      TRH1P(L,NY,NX)=TRH1P(L,NY,NX)+RHP1
      TRH2P(L,NY,NX)=TRH2P(L,NY,NX)+RHP2
      TRH3P(L,NY,NX)=TRH3P(L,NY,NX)+RHP3
      TRF1P(L,NY,NX)=TRF1P(L,NY,NX)+RF1P
      TRF2P(L,NY,NX)=TRF2P(L,NY,NX)+RF2P
      TRC0P(L,NY,NX)=TRC0P(L,NY,NX)+RC0P
      TRC1P(L,NY,NX)=TRC1P(L,NY,NX)+RC1P
      TRC2P(L,NY,NX)=TRC2P(L,NY,NX)+RC2P
      TRM1P(L,NY,NX)=TRM1P(L,NY,NX)+RM1P
      TRH0B(L,NY,NX)=TRH0B(L,NY,NX)+RHB0
      TRH1B(L,NY,NX)=TRH1B(L,NY,NX)+RHB1
      TRH2B(L,NY,NX)=TRH2B(L,NY,NX)+RHB2
      TRH3B(L,NY,NX)=TRH3B(L,NY,NX)+RHB3
      TRF1B(L,NY,NX)=TRF1B(L,NY,NX)+RF1B
      TRF2B(L,NY,NX)=TRF2B(L,NY,NX)+RF2B
      TRC0B(L,NY,NX)=TRC0B(L,NY,NX)+RC0B
      TRC1B(L,NY,NX)=TRC1B(L,NY,NX)+RC1B
      TRC2B(L,NY,NX)=TRC2B(L,NY,NX)+RC2B
      TRM1B(L,NY,NX)=TRM1B(L,NY,NX)+RM1B
      TRXN4(L,NY,NX)=TRXN4(L,NY,NX)+RXN4
      TRXNB(L,NY,NX)=TRXNB(L,NY,NX)+RXNB
      TRXHY(L,NY,NX)=TRXHY(L,NY,NX)+RXHY
      TRXAL(L,NY,NX)=TRXAL(L,NY,NX)+RXAL
      TRXFE(L,NY,NX)=TRXFE(L,NY,NX)+RXFE
      TRXCA(L,NY,NX)=TRXCA(L,NY,NX)+RXCA
      TRXMG(L,NY,NX)=TRXMG(L,NY,NX)+RXMG
      TRXNA(L,NY,NX)=TRXNA(L,NY,NX)+RXNA
      TRXKA(L,NY,NX)=TRXKA(L,NY,NX)+RXKA
      TRXHC(L,NY,NX)=TRXHC(L,NY,NX)+RXHC
      TRXH0(L,NY,NX)=TRXH0(L,NY,NX)+RXH0
      TRXH1(L,NY,NX)=TRXH1(L,NY,NX)+RXH1
      TRXH2(L,NY,NX)=TRXH2(L,NY,NX)+RXH2
      TRX1P(L,NY,NX)=TRX1P(L,NY,NX)+RX1P
      TRX2P(L,NY,NX)=TRX2P(L,NY,NX)+RX2P
      TRBH0(L,NY,NX)=TRBH0(L,NY,NX)+RBH0
      TRBH1(L,NY,NX)=TRBH1(L,NY,NX)+RBH1
      TRBH2(L,NY,NX)=TRBH2(L,NY,NX)+RBH2
      TRB1P(L,NY,NX)=TRB1P(L,NY,NX)+RB1P
      TRB2P(L,NY,NX)=TRB2P(L,NY,NX)+RB2P
      TRALOH(L,NY,NX)=TRALOH(L,NY,NX)+RPALOX
      TRFEOH(L,NY,NX)=TRFEOH(L,NY,NX)+RPFEOX
      TRCACO(L,NY,NX)=TRCACO(L,NY,NX)+RPCACX
      TRCASO(L,NY,NX)=TRCASO(L,NY,NX)+RPCASO
      TRALSI(L,NY,NX)=TRALSI(L,NY,NX)+RQALSI
      TRFESI(L,NY,NX)=TRFESI(L,NY,NX)+RQFESI
      TRCASI(L,NY,NX)=TRCASI(L,NY,NX)+RQCASI
      TRMGSI(L,NY,NX)=TRMGSI(L,NY,NX)+RQMGSI
      TRNASI(L,NY,NX)=TRNASI(L,NY,NX)+RQNASI
      TRKASI(L,NY,NX)=TRKASI(L,NY,NX)+RQKASI
      TRALPO(L,NY,NX)=TRALPO(L,NY,NX)+RPALPX
      TRFEPO(L,NY,NX)=TRFEPO(L,NY,NX)+RPFEPX
      TRCAPD(L,NY,NX)=TRCAPD(L,NY,NX)+RPCADX
      TRCAPH(L,NY,NX)=TRCAPH(L,NY,NX)+RPCAHX
      TRCAPM(L,NY,NX)=TRCAPM(L,NY,NX)+RPCAMX
      TRALPB(L,NY,NX)=TRALPB(L,NY,NX)+RPALBX
      TRFEPB(L,NY,NX)=TRFEPB(L,NY,NX)+RPFEBX
      TRCPDB(L,NY,NX)=TRCPDB(L,NY,NX)+RPCDBX
      TRCPHB(L,NY,NX)=TRCPHB(L,NY,NX)+RPCHBX
      TRCPMB(L,NY,NX)=TRCPMB(L,NY,NX)+RPCMBX
C
C     GO TO NEXT ITERATION
C
1000  CONTINUE
C
C     ITERATIONS COMPLETED
C
C     RESET PH AT END OF ITERATION
C
C     CHY2,COH2=H+,OH- concentration (mol m-3)
C     AHY2,AOH2=H+,OH- activity (mol m-3)
C     TRHY,TROH=total H+,OH- transformation (mol m-3 t-1)
C     DPH2O=H2O solubility product from PARAMETER (mol m-3)2
C     A1=activity coefficient from ‘hour1.f’
C     RHHX=flux to equilibrate H+, OH- at DPH2O
C     TRHY,TROH=total H+,OH- transformations used in ‘redist.f’
C     PH=soil pH
C
      CHY2=CHY0+TRHY(L,NY,NX)
      COH2=COH0+TROH(L,NY,NX)
      IF(CHY2.GT.ZEROC)THEN
      AHY2=CHY2*A1(L,NY,NX)
      ELSE
      AHY2=CHY2
      ENDIF
      IF(COH2.GT.ZEROC)THEN
      AOH2=COH2*A1(L,NY,NX)
      ELSE
      AOH2=COH2
      ENDIF
      S0=AHY2+AOH2
      S1=AMAX1(0.0,S0**2-4.0*(AHY2*AOH2-DPH2O))
      RHHX=0.5*(S0-SQRT(S1))
      AHY2=AHY2-RHHX
      AOH2=AOH2-RHHX
      TRHY(L,NY,NX)=TRHY(L,NY,NX)-RHHX
      TROH(L,NY,NX)=TROH(L,NY,NX)-RHHX
      TBH2O(L,NY,NX)=TBH2O(L,NY,NX)+RHHX
C     IF(L.EQ.1)THEN
C     WRITE(*,23)'RHHR',I,J,NFZ,NX,NY,L,M 
C    2,RHHX,CHY0,COH0,CHY2,COH2,AHY2,AOH2,AHY2*AOH2,DPH2O,TROH(L,NY,NX) 
C    3,ZHY(L,NY,NX),ZOH(L,NY,NX),TROH(L,NY,NX)*VOLW(L,NY,NX)
C     ENDIF
C
C     END RESET PH AT END OF ITERATION
C
C     CONVERT TOTAL ION TRANSFORMATIONS FROM CHANGES IN CONCENTRATION
C     TO CHANGES IN MASS FOR USE IN 'REDIST'
C        (change mol m-3 t-1 to mol t-1)
C
C     VOLW=soil water content (m3)
C     VOLWNH,VOLWNB=water volume in NH4 non-band,band (m3)
C     VOLWPO,VOLWPB=water volume in H2PO4 non-band,band (m3)
C
      TRN4S(L,NY,NX)=TRN4S(L,NY,NX)*VOLWNH
      TRN4B(L,NY,NX)=TRN4B(L,NY,NX)*VOLWNB
      TRN3S(L,NY,NX)=TRN3S(L,NY,NX)*VOLWNH
      TRN3B(L,NY,NX)=TRN3B(L,NY,NX)*VOLWNB
      TRHY(L,NY,NX)=TRHY(L,NY,NX)*VOLW(L,NY,NX)
      TROH(L,NY,NX)=TROH(L,NY,NX)*VOLW(L,NY,NX)
      TRAL(L,NY,NX)=TRAL(L,NY,NX)*VOLW(L,NY,NX)
      TRFE(L,NY,NX)=TRFE(L,NY,NX)*VOLW(L,NY,NX)
      TRCA(L,NY,NX)=TRCA(L,NY,NX)*VOLW(L,NY,NX)
      TRMG(L,NY,NX)=TRMG(L,NY,NX)*VOLW(L,NY,NX)
      TRNA(L,NY,NX)=TRNA(L,NY,NX)*VOLW(L,NY,NX)
      TRKA(L,NY,NX)=TRKA(L,NY,NX)*VOLW(L,NY,NX)
      TRSO4(L,NY,NX)=TRSO4(L,NY,NX)*VOLW(L,NY,NX)
      TRCO3(L,NY,NX)=TRCO3(L,NY,NX)*VOLW(L,NY,NX)
      TRHCO(L,NY,NX)=TRHCO(L,NY,NX)*VOLW(L,NY,NX)
      TRCO2(L,NY,NX)=TRCO2(L,NY,NX)*VOLW(L,NY,NX)
      TBH2O(L,NY,NX)=TBH2O(L,NY,NX)*VOLW(L,NY,NX)
      TRH2O(L,NY,NX)=TRH2O(L,NY,NX)*VOLW(L,NY,NX)
      TRAL1(L,NY,NX)=TRAL1(L,NY,NX)*VOLW(L,NY,NX)
      TRAL2(L,NY,NX)=TRAL2(L,NY,NX)*VOLW(L,NY,NX)
      TRAL3(L,NY,NX)=TRAL3(L,NY,NX)*VOLW(L,NY,NX)
      TRAL4(L,NY,NX)=TRAL4(L,NY,NX)*VOLW(L,NY,NX)
      TRALS(L,NY,NX)=TRALS(L,NY,NX)*VOLW(L,NY,NX)
      TRFE1(L,NY,NX)=TRFE1(L,NY,NX)*VOLW(L,NY,NX)
      TRFE2(L,NY,NX)=TRFE2(L,NY,NX)*VOLW(L,NY,NX)
      TRFE3(L,NY,NX)=TRFE3(L,NY,NX)*VOLW(L,NY,NX)
      TRFE4(L,NY,NX)=TRFE4(L,NY,NX)*VOLW(L,NY,NX)
      TRFES(L,NY,NX)=TRFES(L,NY,NX)*VOLW(L,NY,NX)
      TRCAO(L,NY,NX)=TRCAO(L,NY,NX)*VOLW(L,NY,NX)
      TRCAC(L,NY,NX)=TRCAC(L,NY,NX)*VOLW(L,NY,NX)
      TRCAH(L,NY,NX)=TRCAH(L,NY,NX)*VOLW(L,NY,NX)
      TRCAS(L,NY,NX)=TRCAS(L,NY,NX)*VOLW(L,NY,NX)
      TRMGO(L,NY,NX)=TRMGO(L,NY,NX)*VOLW(L,NY,NX)
      TRMGC(L,NY,NX)=TRMGC(L,NY,NX)*VOLW(L,NY,NX)
      TRMGH(L,NY,NX)=TRMGH(L,NY,NX)*VOLW(L,NY,NX)
      TRMGS(L,NY,NX)=TRMGS(L,NY,NX)*VOLW(L,NY,NX)
      TRNAC(L,NY,NX)=TRNAC(L,NY,NX)*VOLW(L,NY,NX)
      TRNAS(L,NY,NX)=TRNAS(L,NY,NX)*VOLW(L,NY,NX)
      TRKAS(L,NY,NX)=TRKAS(L,NY,NX)*VOLW(L,NY,NX)
      TRH0P(L,NY,NX)=TRH0P(L,NY,NX)*VOLWPO
      TRH1P(L,NY,NX)=TRH1P(L,NY,NX)*VOLWPO
      TRH2P(L,NY,NX)=TRH2P(L,NY,NX)*VOLWPO
      TRH3P(L,NY,NX)=TRH3P(L,NY,NX)*VOLWPO
      TRF1P(L,NY,NX)=TRF1P(L,NY,NX)*VOLWPO
      TRF2P(L,NY,NX)=TRF2P(L,NY,NX)*VOLWPO
      TRC0P(L,NY,NX)=TRC0P(L,NY,NX)*VOLWPO
      TRC1P(L,NY,NX)=TRC1P(L,NY,NX)*VOLWPO
      TRC2P(L,NY,NX)=TRC2P(L,NY,NX)*VOLWPO
      TRM1P(L,NY,NX)=TRM1P(L,NY,NX)*VOLWPO
      TRH0B(L,NY,NX)=TRH0B(L,NY,NX)*VOLWPB
      TRH1B(L,NY,NX)=TRH1B(L,NY,NX)*VOLWPB
      TRH2B(L,NY,NX)=TRH2B(L,NY,NX)*VOLWPB
      TRH3B(L,NY,NX)=TRH3B(L,NY,NX)*VOLWPB
      TRF1B(L,NY,NX)=TRF1B(L,NY,NX)*VOLWPB
      TRF2B(L,NY,NX)=TRF2B(L,NY,NX)*VOLWPB
      TRC0B(L,NY,NX)=TRC0B(L,NY,NX)*VOLWPB
      TRC1B(L,NY,NX)=TRC1B(L,NY,NX)*VOLWPB
      TRC2B(L,NY,NX)=TRC2B(L,NY,NX)*VOLWPB
      TRM1B(L,NY,NX)=TRM1B(L,NY,NX)*VOLWPB
      TRHYSI(L,NY,NX)=TRHYSI(L,NY,NX)*VOLW(L,NY,NX)
      TRXN4(L,NY,NX)=TRXN4(L,NY,NX)*BKVLNH
      TRXNB(L,NY,NX)=TRXNB(L,NY,NX)*BKVLNB
      TRXHY(L,NY,NX)=TRXHY(L,NY,NX)*BKVLX
      TRXAL(L,NY,NX)=TRXAL(L,NY,NX)*BKVLX
      TRXFE(L,NY,NX)=TRXFE(L,NY,NX)*BKVLX
      TRXCA(L,NY,NX)=TRXCA(L,NY,NX)*BKVLX
      TRXMG(L,NY,NX)=TRXMG(L,NY,NX)*BKVLX
      TRXNA(L,NY,NX)=TRXNA(L,NY,NX)*BKVLX
      TRXKA(L,NY,NX)=TRXKA(L,NY,NX)*BKVLX
      TRXHC(L,NY,NX)=TRXHC(L,NY,NX)*BKVLX
C     TRXAL2(L,NY,NX)=TRXAL2(L,NY,NX)*VOLW(L,NY,NX)
C     TRXFE2(L,NY,NX)=TRXFE2(L,NY,NX)*VOLW(L,NY,NX)
      TRXH0(L,NY,NX)=TRXH0(L,NY,NX)*VOLWPO
      TRXH1(L,NY,NX)=TRXH1(L,NY,NX)*VOLWPO
      TRXH2(L,NY,NX)=TRXH2(L,NY,NX)*VOLWPO
      TRX1P(L,NY,NX)=TRX1P(L,NY,NX)*VOLWPO
      TRX2P(L,NY,NX)=TRX2P(L,NY,NX)*VOLWPO
      TRBH0(L,NY,NX)=TRBH0(L,NY,NX)*VOLWPB
      TRBH1(L,NY,NX)=TRBH1(L,NY,NX)*VOLWPB
      TRBH2(L,NY,NX)=TRBH2(L,NY,NX)*VOLWPB
      TRB1P(L,NY,NX)=TRB1P(L,NY,NX)*VOLWPB
      TRB2P(L,NY,NX)=TRB2P(L,NY,NX)*VOLWPB
      TRALOH(L,NY,NX)=TRALOH(L,NY,NX)*VOLW(L,NY,NX)
      TRFEOH(L,NY,NX)=TRFEOH(L,NY,NX)*VOLW(L,NY,NX)
      TRCACO(L,NY,NX)=TRCACO(L,NY,NX)*VOLW(L,NY,NX)
      TRCASO(L,NY,NX)=TRCASO(L,NY,NX)*VOLW(L,NY,NX)
      TRALSI(L,NY,NX)=TRALSI(L,NY,NX)*VOLW(L,NY,NX)
      TRFESI(L,NY,NX)=TRFESI(L,NY,NX)*VOLW(L,NY,NX)
      TRCASI(L,NY,NX)=TRCASI(L,NY,NX)*VOLW(L,NY,NX)
      TRMGSI(L,NY,NX)=TRMGSI(L,NY,NX)*VOLW(L,NY,NX)
      TRNASI(L,NY,NX)=TRNASI(L,NY,NX)*VOLW(L,NY,NX)
      TRKASI(L,NY,NX)=TRKASI(L,NY,NX)*VOLW(L,NY,NX)
      TRALPO(L,NY,NX)=TRALPO(L,NY,NX)*VOLWPO
      TRFEPO(L,NY,NX)=TRFEPO(L,NY,NX)*VOLWPO
      TRCAPD(L,NY,NX)=TRCAPD(L,NY,NX)*VOLWPO
      TRCAPH(L,NY,NX)=TRCAPH(L,NY,NX)*VOLWPO
      TRCAPM(L,NY,NX)=TRCAPM(L,NY,NX)*VOLWPO
      TRALPB(L,NY,NX)=TRALPB(L,NY,NX)*VOLWPB
      TRFEPB(L,NY,NX)=TRFEPB(L,NY,NX)*VOLWPB
      TRCPDB(L,NY,NX)=TRCPDB(L,NY,NX)*VOLWPB
      TRCPHB(L,NY,NX)=TRCPHB(L,NY,NX)*VOLWPB
      TRCPMB(L,NY,NX)=TRCPMB(L,NY,NX)*VOLWPB
C     IF(L.EQ.3)THEN
C     WRITE(*,1121)'TRHY',I,J,NFZ,L,TRHY(L,NY,NX),RHY,VOLW(L,NY,NX) 
C     WRITE(*,1121)'TRCA',I,J,NFZ,L,TRCA(L,NY,NX),TRCO3(L,NY,NX)
C    2,TRHCO(L,NY,NX),TRCO2(L,NY,NX),TRCACO(L,NY,NX),ZCA(L,NY,NX)
C     WRITE(*,1121)'TRFE4',I,J,NFZ,L,TRFE4(L,NY,NX),RFE4
C    2,VOLW(L,NY,NX)
1121  FORMAT(A8,4I4,30F16.8)
C     ENDIF
C
C     IF NO SALT IS SELECTED IN SITE FILE THEN A SUBSET
C     OF SOIL EQUILIBRIA REACTIONS ARE SOLVED: MOSTLY THOSE
C     FOR PHOSPHORUS AND CO-REACTANTS
C
      ELSE
C
C     TIME STEPS
C
C     TPDH,TADAH,TADCH,TSLH,TRWH=maximum  
C        precipitation,adsorption,solute,silicate weathering
C        rates from PARAMETER (mol m-3 h-1)
C     TPD,TADA,TADC,TSL,TRW=maximum 
C        precipitation,adsorption,solute,silicate weathering
C        rates from PARAMETER (mol m-3 h-1)
C     XNFH=time step for solute transformations from ‘wthr.f’ (h t-1)
C
      TPD=TPDH*XNFH 
      TPZ=TPZH*XNFH
      TADA=TADAH*XNFH 
      TADC=TADCH*XNFH 
      TSL=TSLH*XNFH 
      TSZ=TSZH*XNFH 
      TRW=TRWH*XNFH 
      TPA=TRW 
C
C     INITIALIZE ION CONCENTRATIONS, ACTIVITIES
C
C     CCEC,XCEC=cation exchange concentration,capacity (mol Mg-1,mol)
C     BKVLX=soil mass (Mg)
C     C*1,A*1=aqueous cation concentration,activity (mol m-3) 
C     cation code
C        :HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ non-band,NB=NH4+ band
C
      IF(BKVLX.GT.ZEROS2(NY,NX))THEN
      CCEC=AMAX1(ZEROC,XCEC(L,NY,NX)/BKVLX)
      ELSE
      CCEC=ZERO
      ENDIF
      IF(BKVLX.GT.ZEROS2(NY,NX))THEN
      CCEC=AMAX1(ZEROC,XCEC(L,NY,NX)/BKVLX)
      ELSE
      CCEC=ZERO
      ENDIF
      CHY1=AMAX1(ZEROC,ZHY(L,NY,NX)/VOLW(L,NY,NX))
      COH1=AMAX1(ZEROC,ZOH(L,NY,NX)/VOLW(L,NY,NX))
      CAL1=AMAX1(ZEROC,ZAL(L,NY,NX)/VOLW(L,NY,NX))
      CFE1=AMAX1(ZEROC,ZFE(L,NY,NX)/VOLW(L,NY,NX))
      CCA1=AMAX1(ZEROC,ZCA(L,NY,NX)/VOLW(L,NY,NX))
      CMG1=AMAX1(ZEROC,ZMG(L,NY,NX)/VOLW(L,NY,NX))
      CNA1=AMAX1(ZEROC,ZNA(L,NY,NX)/VOLW(L,NY,NX))
      CKA1=AMAX1(ZEROC,ZKA(L,NY,NX)/VOLW(L,NY,NX))
      CSO41=AMAX1(ZEROC,ZSO4(L,NY,NX)/VOLW(L,NY,NX))
      CCL1=AMAX1(ZEROC,ZCL(L,NY,NX)/VOLW(L,NY,NX))
      CCO21=AMAX1(ZEROC,CCO2S(L,NY,NX)/12.0)
      CHCO31=AMAX1(ZEROC,CCO21*DPHCO/CHY1)
      CCO31=AMAX1(ZEROC,CCO21*DPCO3/CHY1**2)
      CALO2=AMAX1(ZEROC,ZALOH2(L,NY,NX)/VOLW(L,NY,NX))
      CFEO2=AMAX1(ZEROC,ZFEOH2(L,NY,NX)/VOLW(L,NY,NX))
      AHY1=CHY1*A1(L,NY,NX)
      AOH1=COH1*A1(L,NY,NX)
      AAL1=CAL1*A3(L,NY,NX)
      AFE1=CFE1*A3(L,NY,NX)
      ACA1=CCA1*A2(L,NY,NX)
      AMG1=CMG1*A2(L,NY,NX)
      ANA1=CNA1*A1(L,NY,NX)
      AKA1=CKA1*A1(L,NY,NX)
      ASO41=CSO41*A2(L,NY,NX)
      ACO21=CCO21
      AHCO31=CHCO31*A1(L,NY,NX)
      ACO31=CCO31*A2(L,NY,NX)
      AALO2=CALO2*A1(L,NY,NX)
      AFEO2=CFEO2*A1(L,NY,NX)
      AH1P1=CH1P1*A2(L,NY,NX)
      AH2P1=CH2P1*A1(L,NY,NX)
      AH1PB=CH1PB*A2(L,NY,NX)
      AH2PB=CH2PB*A1(L,NY,NX)
      AN41=CN41*A1(L,NY,NX)
      AN4B=CN4B*A1(L,NY,NX)
      AN31=CN31
      AN3B=CN3B
C
C     PRECIPITATION-DISSOLUTION REACTIONS FOR H2PO4 
C     IN NON-BAND, BAND ZONES CALCULATED FROM ACTIVITIES OF REACTANTS
C     AND PRODUCTS THROUGH SOLUTIONS FOR THEIR EQUILIBRIUM
C     CONCENTRATIONS 
C
C     for all precipitation-dissolution reactions:
C
C     SP*=solubility product from PARAMETER above
C     A*1=ion activity (mol m-3)
C     A*Q=ion equilibrium activity at SP* (mol m-3)
C     A1,A2,A3=activity coefficients from ‘hour1.f’ 
C     R*X=precipitation(+ve) or dissolution (-ve) rate 
C        (mol m-3 t-1)
C     XMIN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONN=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     TPD,TPA=maximum precipitation rate (mol m-3 t-1)
C     phosphorus solute code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C          :*F2P*=F1H2PO4-,*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C          :*M1P*=MgHPO4,*COO*=COO-
C          :*1=non-band,*B=band
C
C     PHOSPHORUS PRECIPITATION-DISSOLUTION IN NON-BAND SOIL ZONE
C
      IF(VOLWPO.GT.ZEROS2(NY,NX))THEN
C
C     ALUMINUM PHOSPHATE (VARISCITE)
C
      XMIN=FIONN*AMIN1(AAL1,AH2P1)
      AAL1Q=SPALO/AOH1**3
      AH2PQ=SHA0P2*AHY1**2/AAL1Q
      RPALPX=AMAX1(-AMAX1(0.0,PALPO1),-TPD,-XMIN
     2,(AMIN1(TPD,XMIN,AH2P1-AH2PQ)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24 
C    2.AND.L.EQ.2)THEN
C     WRITE(*,1117)'RPALPX',I,J,NFZ,L 
C    2,RPALPX,AH2P1,AH2PQ,SHA0P2,AAL1Q,AHY1,AAL1 
C    2,PALPO1,ZAL(L,NY,NX),SPALO/AOH1**3 
C     ENDIF
C
C     IRON PHOSPHATE (STRENGITE)
C
      XMIN=FIONN*AMIN1(AFE1,AH2P1)
      AFE1Q=SPFEO/AOH1**3
      AH2PQ=SHF0P2*AHY1**2/AFE1Q
      RPFEPX=AMAX1(-AMAX1(0.0,PFEPO1),-TPD,-XMIN
     2,(AMIN1(TPD,XMIN,AH2P1-AH2PQ)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24 
C    2.AND.L.EQ.2)THEN
C     WRITE(*,1117)'RPFEPX',I,J,NFZ,L 
C     2,RPFEPX,AH2P1,AH2PQ,SHF0P2,AFE1Q,AHY1,AFE1 
C    2,PFEPO1,ZFE(L,NY,NX),SPFEO/AOH1**3  
C     ENDIF
C
C     DICALCIUM PHOSPHATE
C
      XMIN=FIONN*AMIN1(ACA1,AH1P1)
      AH1PQ=SPCAD/ACA1 
      RPCADX=AMAX1(-AMAX1(0.0,PCAPD1),-TSL
     2,(AMIN1(TPD,XMIN,AH1P1-AH1PQ)))
C
C     HYDROXYAPATITE
C
      XMIN=FIONN*AMIN1(ACA1,AH2P1)
      AH2PH=(SHCAH2*AHY1**7/ACA1**5)**0.333
      RPCAHX=AMAX1(-AMAX1(0.0,PCAPH1),-XMIN,-TPA
     2,AMIN1(XMIN,TPA,(AH2P1-AH2PH)))
C     IF(L.EQ.9)THEN
C     WRITE(*,1117)'RPCADX',I,J,NFZ,L,RPCADX,AH2P1,SYCAD2 
C    2,ACA1,SPCAC,DPCO3,ACO31,ACO21,AHY1,AOH1,PH(L,NY,NX),PCAPD1,BKVLW
C     WRITE(*,1117)'RPCAHX',I,J,NFZ,L,RPCAHX,AH2P1,AH2PH,SHCAH2 
C    2,ACA1,SPCAC,AHY1,AOH1,PH(L,NY,NX),PCAPH1,TPD
C    3,XMIN,CCA1,A2(L,NY,NX)
C     ENDIF
C
C     MONOCALCIUM PHOSPHATE
C
      XMIN=FIONN*AMIN1(ACA1,AH2P1)
      AH2PQ=(SPCAM/ACA1)**0.5
      RPCAMX=AMAX1(-AMAX1(0.0,PCAPM1),-TPD 
     2,(AMIN1(TPD,XMIN,AH2P1-AH2PQ)))
C     IF(I.GT.315)THEN
C     WRITE(*,1117)'RPCADX',I,J,NFZ,L,RPCADX,AH2P1,PCAPD1,RPCAHX 
C    2,AH2PH,SYA0P2,AAL1,AOH1,SHCAH2,ACA1,ACO21,ACO31,PCAPH1 
C    3,VOLWPO,SPCAC/ACO31,CCA(L,NY,NX),H2PO4(L,NY,X)
C    4,VOLW(L,NY,NX),ZCA(L,NY,NX),CCO2S(L,NY,NX)
1117  FORMAT(A8,4I4,30E12.4)
C     ENDIF
C
C     PHOSPHORUS ANION EXCHANGE IN NON-BAND SOIL ZONE
C     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
C     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
C     EXCHANGE SITES
C
C     XAEC=anion exchange capacity (mol)
C     VOLWPO=soil water volume in non-band (m3)
C     TADA=maximum anion adsorption rate (mol m-3 t-1)
C     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with 
C        R-OH2,R-OH from PARAMETER
C     RXH2P,RYH2P=H2PO4 adsorption(+ve), desorption(-ve) 
C        with R-OH2,R-OH in non-band (mol P m-3 t-1)
C     AH2P1=H2PO4 activity in non-band (mol P m-3)
C     AH2PQ2,AH2PQ1=H2PO4 equilibrium activity in non-band at SP*,SX*
C        (mol P m-3)
C     XH2P1=exchangeable H2PO4 concentration in non-band (mol Mg-1)
C     XOH11,XOH21=concentration of adsorption sites R-OH,
C        R-OH2 in non-band (mol Mg-1)
C     AHY1,AOH1=H+,OH- activity (mol m-3)
C     TADA=maximum anion adsorption rate (mol m-3 t-1)
C     XMINP*,XMINN*=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONN=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C
      IF(XAEC(L,NY,NX).GT.ZEROS(NY,NX))THEN
C
C     H2PO4 EXCHANGE IN NON-BAND SOIL ZONE FROM SOLUTION 
C     FOR EQUILIBRIUM AMONG H2PO4-, H+, OH-, X-OH AND X-H2PO4
C
      XMINN2=FIONN*XH2P1 
      XMINP2=FIONN*AMIN1(AH2P1,XOH21) 
      SPH2P=SXH2P*AHY1*AOH1 
      AH2PQ2=SPH2P*XH2P1/XOH21
      RXH2P=AMAX1(-TADA,-XMINN2
     2,(AMIN1(TADA,XMINP2,AH2P1-AH2PQ2)))
      XMINN1=FIONN*XH2P1 
      XMINP1=FIONN*AMIN1(AH2P1,XOH11) 
      SYH2P=SXH2P*AOH1
      AH2PQ1=SYH2P*XH2P1/XOH11
      RYH2P=AMAX1(-TADA,-XMINN1
     2,(AMIN1(TADA,XMINP1,AH2P1-AH2PQ1)))
C
C     HPO4 EXCHANGE IN NON-BAND SOIL ZONE FROM SOLUTION 
C     FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH AND X-HPO4
C
C     SPH1P=equilibrium constant for HPO4 exchange with R-OH 
C        from PARAMETER
C     RXH1P=HPO4 adsorption(+ve), desorption(-ve) with R-OH 
C        in non-band (mol P m-3 t-1)
C     AH1P1=HPO4 activity in non-band (mol P m-3)
C     AH1PQ=HPO4 equilibrium activity in non-band at SP*(mol P m-3)
C     XH1P1=exchangeable HPO4 concentration in non-band (mol P Mg-1)
C     XOH11=concentration of adsorption sites R-OH
C        in non-band (mol Mg-1)
C     TADA=maximum anion adsorption rate (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONN=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C
      XMINN=FIONN*XH1P1
      XMINP=FIONN*AMIN1(AH1P1,XOH11) 
      SPH1P=SXH1P*DPH2O/DPH2P
      AH1PQ=SPH1P*XH1P1/XOH11
      RXH1P=AMAX1(-TADA,-XMINN
     2,(AMIN1(TADA,XMINP,AH1P1-AH1PQ)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24 
C    2.AND.L.EQ.2)THEN
C     WRITE(*,24)'RXH2P',I,J,NFZ,L
C    2,RXH2P,TADA,XOH21,AH2P1,AH2PQ2,SPH2P,XH2P1
C    2,BKVLW,CH2P1,A1(L,NY,NX)
C    2,H2PO4(L,NY,NX),RH2PX,VOLWPX 
C    2,XOH01+XOH11+XOH21+XH2P1+XH1P1
C     WRITE(*,24)'RYH2P',I,J,NFZ,L
C    2,RYH2P,TADA,XOH11,AH2P1,AH2PQ1,SYH2P,XH2P1,AOH1
C    2,BKVLW,CH2P1,A1(L,NY,NX)
C    2,H2PO4(L,NY,NX),RH2PX,VOLWPX
C    2,XOH01+XOH11+XOH21+XH2P1+XH1P1
C     WRITE(*,24)'RXH1P',I,J,NFZ,L
C    2,RXH1P,TADA,XOH11,AH1P1,AH1PQ,SPH1P,XH1P1
C    2,BKVLW,CH1P1,A1(L,NY,NX)
C    2,H1PO4(L,NY,NX),RH1PX,VOLWPX
C    2,XOH01+XOH11+XOH21+XH2P1+XH1P1
C     ENDIF
      ELSE
      RXH2P=0.0
      RYH2P=0.0
      RXH1P=0.0
      ENDIF
C
C     H2PO4-H+HPO4 IN NON-BAND
C
C     DPH2P=dissociation constant from PARAMETER
C     AH1P1,AH2P1=HPO4,H2PO4 activity in non-band (mol m-3)
C     AH1PQ=HPO4 equilibrium activity in non-band at DPH2P (mol P m-3)
C     RH2P=H2PO4-H+HPO4 dissociation(-ve),association(+ve) in non-band 
C        (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONN=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     TSL=maximum reaction rate (mol m-3 t-1)
C     
      XMINN=FIONN*AH2P1
      XMINP=FIONN*AH1P1
      AH1P1Q=DPH2P*AH2P1/AHY1
      RH2P=AMAX1(-TSL,-XMINN
     2,(AMIN1(TSL,XMINP,AH1P1-AH1P1Q)))
      ELSE
      RPALPX=0.0
      RPFEPX=0.0
      RPCADX=0.0
      RPCAHX=0.0
      RPCAMX=0.0
      RXH2P=0.0
      RYH2P=0.0
      RXH1P=0.0
      RH2P=0.0
      ENDIF
C     IF(J.EQ.1)THEN
C     WRITE(*,24)'PO4',I,J,NFZ,L,AH2P1,PALPO1 
C    2,PFEPO1,PCAPD1,PCAPH1,PCAPM1,AH2PF,AH2PH 
C    2,RPALPX,RPFEPX,RPCADX,RPCAHX,RPCAMX
C    3,XH2P1,RXH2P,RYH2P
C    3,AAL1,AFE1,ACA1,AHY1,AOH1
2222  FORMAT(A8,4I4,/,9F16.8,/,9F16.8,/,9F16.8,/,9F16.8,/
     2,9F16.8,/,9F16.8,/,9F16.8,/,9F16.8,/,9F16.8)
C     ENDIF
C
C     PHOSPHORUS PRECIPITATION-DISSOLUTION IN BAND SOIL ZONE
C
C     SP*=solubility product from PARAMETER above
C     A*B=ion activity in band (mol m-3)
C     A*BQ=ion equilibrium activity at SP* in band (mol m-3)
C     A1,A2,A3=activity coefficients from ‘hour1.f’ 
C     R*X=precipitation(+ve) or dissolution (-ve) rate 
C        (mol m-3 t-1)
C     XMIN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONN=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     TPD,TPA=maximum precipitation rate (mol m-3 t-1)
C     phosphorus solute code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C          :*F2P*=F1H2PO4-,*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C          :*M1P*=MgHPO4,*COO*=COO-
C          :*1=non-band,*B=band
C
      IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
C
C     ALUMINUM PHOSPHATE (VARISCITE)
C
      XMIN=FIONN*AMIN1(AAL1,AH2PB)
      AAL1Q=SPALO/AOH1**3
      AH2PBQ=SHA0P2*AHY1**2/AAL1Q
      RPALBX=AMAX1(-AMAX1(0.0,PALPOB),-TPD,-XMIN
     2,(AMIN1(TPD,XMIN,AH2PB-AH2PBQ)))
C
C     IRON PHOSPHATE (STRENGITE)
C
      XMIN=FIONN*AMIN1(AFE1,AH2PB)
      AFE1Q=SPFEO/AOH1**3
      AH2PBQ=SHF0P2*AHY1**2/AFE1Q
      RPFEBX=AMAX1(-AMAX1(0.0,PFEPOB),-TPD,-XMIN
     2,(AMIN1(TPD,XMIN,AH2PB-AH2PBQ)))
C
C     DICALCIUM PHOSPHATE
C
      XMIN=FIONN*AMIN1(ACA1,AH1PB)
      AH1PQ=SPCAD/ACA1 
      RPCDBX=AMAX1(-AMAX1(0.0,PCAPDB),-TSL
     2,(AMIN1(TPD,XMIN,AH1PB-AH1PQ)))
C
C     HYDROXYAPATITE CA5(PO4)3OH
C
      XMIN=FIONN*AMIN1(ACA1,AH2PB)
      AH2PH=(SHCAH2*AHY1**7/ACA1**5)**0.333
      RPCHBX=AMAX1(-AMAX1(0.0,PCAPHB),-TPZ,-XMIN
     2,AMIN1(TPZ,XMIN,(AH2PB-AH2PH)))
C
C     MONOCALCIUM PHOSPHATE CA(H2PO4)2
C
      XMIN=FIONN*AMIN1(ACA1,AH2PB)
      AH2PQ=(SPCAM/ACA1)**0.5
      RPCMBX=AMAX1(-AMAX1(0.0,PCAPMB),-TPD
     2,(AMIN1(TPD,XMIN,AH2PB-AH2PQ)))
C     IF(I.GT.315)THEN
C     WRITE(*,1117)'RPCMBX',I,J,NFZ,L
C    2,RPCMBX,ACA1,AH2PB,SPCAM,CCA1,A2(L,NY,NX),S0,S1 
C    2,PCAPMB,BKVLW,TPD,PCPMB(L,NY,NX),VLPOB(L,NY,NX)
C    3,PALPOB,PFEPOB,PCAPHB,PCAPDB,PCAPMB
C     ENDIF
C
C     PHOSPHORUS ANION EXCHANGE IN BAND SOIL ZONE
C     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
C     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
C     EXCHANGE SITES
C
C     XAEC=anion exchange capacity (mol)
C     VOLWPB=soil water volume in band (m3)
C     TADA=maximum anion adsorption rate (mol m-3 t-1)
C     SPH2P,SXH2P=equilibrium constant for H2PO4 exchange with 
C        R-OH2,R-OH from PARAMETER
C     RXH2B,RYH2B=H2PO4 adsorption(+ve), desorption(-ve) 
C        with R-OH2,R-OH in band (mol P m-3 t-1)
C     AH2PB=H2PO4 activity in band (mol P m-3)
C     AH2BQ2,AH2BQ1=H2PO4 equilibrium activity in band at SP*,SX*
C        (mol P m-3)
C     XH2PB=exchangeable H2PO4 concentration in band (mol Mg-1)
C     XH11B,XH21B=concentration of adsorption sites R-OH,
C        R-OH2 in band (mol Mg-1)
C     AHY1,AOH1=H+,OH- activity (mol m-3)
C     TADA=maximum anion adsorption rate (mol m-3 t-1)
C     XMINP*,XMINN*=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONN=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C
      IF(XAEC(L,NY,NX).GT.ZEROS(NY,NX))THEN
C
      XMINN2=FIONN*X2P1B
      XMINP2=FIONN*AMIN1(AH2PB,XH21B) 
      SPH2P=SXH2P*DPH2O 
      AH2BQ2=SPH2P*X2P1B/XH21B
      RXH2B=AMAX1(-TADA,-XMINN2
     2,(AMIN1(TADA,XMINP2,AH2PB-AH2BQ2)))
      XMINN1=FIONX*X2P1B 
      XMINP1=FIONX*AMIN1(AH2PB,XH11B) 
      SYH2P=SXH2P*AOH1
      AH2BQ1=SYH2P*X2P1B/XH11B
      RYH2B=AMAX1(-TADA,-XMINN1
     2,(AMIN1(TADA,XMINP1,AH2PB-AH2BQ1)))
C 
C     HPO4 EXCHANGE IN BAND SOIL ZONE FROM SOLUTION 
C     FOR EQUILIBRIUM AMONG HPO4--, H+, OH-, X-OH AND X-HPO4
C
C     SPH1P=equilibrium constant for HPO4 exchange with R-OH
C        from PARAMETER
C     RXH1B=HPO4 adsorption(+ve), desorption(-ve) with R-OH 
C        in band (mol m-3 t-1)
C     AH1PB=HPO4 activity in band (mol P m-3)
C     AH1PBQ=HPO4 equilibrium activity in band at SP*(mol P m-3)
C     X1P1B=exchangeable HPO4 concentration in band (mol Mg-1)
C     XH11B=concentration of adsorption sites R-OH,
C        in band (mol Mg-1)
C     TADA=anion adsorption rate constant (t-1)
C
      XMINN=FIONX*X1P1B 
      XMINP=FIONX*AMIN1(AH1PB,XH11B) 
      SPH1P=SXH1P*DPH2O/DPH2P
      AH1BQ=SPH1P*X1P1B/XH11B
      RXH1B=AMAX1(-TADA,-XMINN
     2,(AMIN1(TADA,XMINP,AH1PB-AH1BQ)))
      ELSE
      RXH2B=0.0
      RYH2B=0.0
      RXH1B=0.0
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,1114)'RXH2B',I,J,NFZ,NX,NY,L,RXH2B,TADA,XMINN2,XMINP2
C     2,XH21B,AH2PB,SPH2P,X2P1B,BKVLW,VOLWPB,XH2PB(L,NY,NX)
C     WRITE(*,1114)'RYH2B',I,J,NFZ,NX,NY,L,RYH2B,TADA,XMINN1,XMINP1
C    2,XH11B,AH2PB,SXH2P,X2P1B,AOH1,BKVLW,VOLWPB 
C     WRITE(*,1114)'RXH1B',I,J,NFZ,NX,NY,L,RXH1B,TADA,XMINN,XMINP 
C    2,XH11B,AH1PB,SPH1P,X1P1B,BKVLW,VOLWPB,XH1PB(L,NY,NX) 
1114  FORMAT(A8,6I4,40E12.4)
C     ENDIF
C
C     H2PO4-H+HPO4 IN BAND
C
C     DPH2P=dissociation constant from PARAMETER
C     AH1PB,AH2PB=HPO4,H2PO4 activity in band (mol m-3)
C     AH1PBQ=HPO4 equilibrium activity in band at DPH2P (mol P m-3)
C     RH2B=H2PO4-H+HPO4 dissociation(-ve),association(+ve)in band 
C        (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONN=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     TSL=maximum reaction rate (mol m-3 t-1)
C
      XMINN=FIONN*AH2PB
      XMINP=FIONN*AH1PB
      AH1PBQ=DPH2P*AH2PB/AHY1
      RH2B=AMAX1(-TSL,-XMINN
     2,(AMIN1(TSL,XMINP,AH1PB-AH1PBQ)))
      ELSE
      RPALBX=0.0
      RPFEBX=0.0
      RPCDBX=0.0
      RPCHBX=0.0
      RPCMBX=0.0
      RXH2B=0.0
      RYH2B=0.0
      RXH1B=0.0
      RH2B=0.0
      ENDIF
C
C     CATION EXCHANGE FROM GAPON SELECTIVITY COEFFICIENTS
C     FOR CA-NH4, CA-H, CA-AL
C
      IF(XCEC(L,NY,NX).GT.ZEROS(NY,NX))THEN
C
C     CATION CONCENTRATIONS
C
C     EQUILIBRIUM X-CA CONCENTRATION FROM CEC, GAPON COEFFICIENTS
C     AND CATION CONCENTRATIONS
C
C     CCEC,XCEC=cation exchange concentration,capacity (mol Mg-1,mol)
C     A*1=cation solute activity (mol m-3)
C     cation code:AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C          :N4=NH4+ non-band,NB=NH4+ band,HY=H+
C
      AALX=AMAX1(ZEROC,AAL1)**0.333
      AFEX=AMAX1(ZEROC,AFE1)**0.333
      ACAX=AMAX1(ZEROC,ACA1)**0.500
      AMGX=AMAX1(ZEROC,AMG1)**0.500
C
C     EQUILIBRIUM X-CA CONCENTRATION FROM CEC AND CATION
C     CONCENTRATIONS
C
C     CCEC,XCEC=total cation exchange concentration,capacity 
C        (mol Mg-1,mol)
C     XCAX=equilibrium X-Ca concentration (mol+ Mg-1)
C     A*1=cation activity (mol m-3)
C     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
C        CA-NH4,CA-H,CA-AL,CA-MG,CA-NA,CA-K
C     X*Q,X*Y=equilibrium,current exchangeable cation concentrations
C        (mol Mg-1)
C     XTLQ,XTL1=total equilibrium exchangeable concentration 
C        (mol, mol(+) Mg-1
C     cation code:AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ non-band,NB=NH4+ band,HY=H+
C
      IF(ACAX.GT.ZERO.AND.CCEC.GT.ZERO)THEN 
      XCAX=CCEC/(1.0
     2+GKC4(L,NY,NX)*AN41/ACAX*VLNH4(L,NY,NX)
     3+GKC4(L,NY,NX)*AN4B/ACAX*VLNHB(L,NY,NX)
     4+GKCH(L,NY,NX)*AHY1/ACAX
     5+GKCA(L,NY,NX)*AALX/ACAX*3.0 
     6+GKCA(L,NY,NX)*AFEX/ACAX*3.0
     7+GKCM(L,NY,NX)*AMGX/ACAX*2.0
     8+GKCN(L,NY,NX)*ANA1/ACAX
     9+GKCK(L,NY,NX)*AKA1/ACAX)
      XN4Q=XCAX*AN41/ACAX*GKC4(L,NY,NX)*VLNH4(L,NY,NX)
      XNBQ=XCAX*AN4B/ACAX*GKC4(L,NY,NX)*VLNHB(L,NY,NX)
      XHYQ=XCAX*AHY1/ACAX*GKCH(L,NY,NX)
      XALQ=XCAX*AALX/ACAX*GKCA(L,NY,NX)*3.0
      XFEQ=XCAX*AFEX/ACAX*GKCA(L,NY,NX)*3.0
      XCAQ=XCAX*ACAX/ACAX*2.0
      XMGQ=XCAX*AMGX/ACAX*GKCM(L,NY,NX)*2.0
      XNAQ=XCAX*ANA1/ACAX*GKCN(L,NY,NX)
      XKAQ=XCAX*AKA1/ACAX*GKCK(L,NY,NX)
      XTLQ=XN4Q+XNBQ+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ 
      XTL1=XN41+XNB1+XHY1+XAL1*3.0+XFE1*3.0
     2+XCA1*2.0+XMG1*2.0+XNA1+XKA1 
      IF(XTLQ.GT.ZERO)THEN
      FX=CCEC/XTLQ
      ELSE
      FX=0.0
      ENDIF
      IF(XTL1.GT.ZERO)THEN
      FY=CCEC/XTL1
      ELSE
      FY=0.0
      ENDIF
      XN4Q=FX*XN4Q
      XNBQ=FX*XNBQ
      XHYQ=FX*XHYQ
      XALQ=FX*XALQ
      XFEQ=FX*XFEQ
      XCAQ=FX*XCAQ
      XMGQ=FX*XMGQ
      XNAQ=FX*XNAQ
      XKAQ=FX*XKAQ
      XN4Y=FY*XN41
      XNBY=FY*XNB1
      XHYY=FY*XHY1
      XALY=FY*XAL1*3.0
      XFEY=FY*XFE1*3.0
      XCAY=FY*XCA1*2.0
      XMGY=FY*XMG1*2.0
      XNAY=FY*XNA1
      XKAY=FY*XKA1
C
C     CATION EXCHANGE IN NON-BAND AND BAND SOIL ZONES
C
C     FIONC=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     BKVLW=soil mass:water (Mg m-3)
C     TADCC=maximum cation adsorption rate (mol Mg t-1)
C     XMIN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C
C     for NH4+
C
C     RXN4Q,RXNBQ=NH4 equilibrium desorption(-ve), adsorption(+ve) 
C        rate in non-band,band (mol N m-3 t-1)
C     RXN4,RXNBQ=NH4 current desorption(-ve), adsorption(+ve)
C        rate in non-band,band (mol N m-3 t-1)
C     XN4Q,XN4Y=equilibrium,current exchangeable NH4 concentration
C        non-band (mol Mg-1)
C     CN41,AN41=aqueous NH4 concentration,activity non-band (mol N m-3)
C     XNBQ,XNBY=equilibrium,current exchangeable NH4 concentration
C        band (mol N Mg-1)
C     CN4B,AN4B=aqueous NH4 concentration,activity band (mol N m-3)
C
C     for all other cations:
C
C     RX*Q,RX*=equilibrium,current cation desorption(-ve),
C        adsorption(+ve) rate (mol m-3 t-1)
C     X*Q,X*Y=equilibrium,current exchangeable cation concentrations
C        (mol Mg-1)
C     C*1,A*1=aqueous cation concentration,activity (mol m-3) 
C     cation code
C        :HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ non-band,NB=NH4+ band
C
      FIONC=FIONN/BKVLW
      TADCC=TADC/BKVLW
      XMIN=FIONC*AMIN1(XN41,CN41)
      RXN4Q=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XN4Q-XN4Y)))*VLNH4(L,NY,NX)
      XMIN=FIONC*AMIN1(XNB1,CN4B)
      RXNBQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XNBQ-XNBY)))*VLNHB(L,NY,NX)
      XMIN=FIONC*AMIN1(XHY1,CHY1)
      RXHYQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XHYQ-XHYY)))
      XMIN=FIONC*3.0*AMIN1(XAL1,CAL1)
      RXALQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XALQ-XALY)))
      XMIN=FIONC*3.0*AMIN1(XFE1,CFE1)
      RXFEQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XFEQ-XFEY)))
      XMIN=FIONC*2.0*AMIN1(XCA1,CCA1)
      RXCAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XCAQ-XCAY)))
      XMIN=FIONC*2.0*AMIN1(XMG1,CMG1)
      RXMGQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XMGQ-XMGY)))
      XMIN=FIONC*AMIN1(XNA1,CNA1)
      RXNAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XNAQ-XNAY)))
      XMIN=FIONC*AMIN1(XKA1,CKA1)
      RXKAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XKAQ-XKAY)))
      TXXX=RXN4Q+RXNBQ+RXHYQ+RXALQ+RXFEQ+RXCAQ+RXMGQ+RXNAQ+RXKAQ
      TXXY=ABS(RXN4Q)+ABS(RXNBQ)+ABS(RXHYQ)+ABS(RXALQ)+ABS(RXFEQ)
     3+ABS(RXCAQ)+ABS(RXMGQ)+ABS(RXNAQ)+ABS(RXKAQ)
      RXN4=RXN4Q-TXXX*ABS(RXN4Q)/TXXY
      RXNB=RXNBQ-TXXX*ABS(RXNBQ)/TXXY
      RXHY=RXHYQ-TXXX*ABS(RXHYQ)/TXXY
      RXAL=(RXALQ-TXXX*ABS(RXALQ)/TXXY)/3.0
      RXFE=(RXFEQ-TXXX*ABS(RXFEQ)/TXXY)/3.0
      RXCA=(RXCAQ-TXXX*ABS(RXCAQ)/TXXY)/2.0
      RXMG=(RXMGQ-TXXX*ABS(RXMGQ)/TXXY)/2.0
      RXNA=RXNAQ-TXXX*ABS(RXNAQ)/TXXY
      RXKA=RXKAQ-TXXX*ABS(RXKAQ)/TXXY
      ELSE
      RXN4=0.0
      RXNB=0.0
      RXHY=0.0
      RXAL=0.0
      RXFE=0.0
      RXCA=0.0
      RXMG=0.0
      RXNA=0.0
      RXKA=0.0
      ENDIF
      ELSE
      RXN4=0.0
      RXNB=0.0
      RXHY=0.0
      RXAL=0.0
      RXFE=0.0
      RXCA=0.0
      RXMG=0.0
      RXNA=0.0
      RXKA=0.0
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.L.EQ.1)THEN
C     IF(NX.EQ.1.AND.L.EQ.1)THEN
C     WRITE(*,2222)'RXN4X',I,J,NFZ,L
C    2,XN4Q,XHYQ,XALQ,XFEQ,XCAQ,XMGQ,XNAQ,XKAQ,XNBQ
C    2,XN41,XHY1,XAL1,XFE1,XCA1,XMG1,XNA1,XKA1,XNB1
C    2,XN4Y,XHYY,XALY,XFEY,XCAY,XMGY,XNAY,XKAY,XNBY
C    2,AN41,AHY1,AAL1,AFE1,ACA1,AMG1,ANA1,AKA1,AN4B
C    2,RXN4,RXHY,RXAL,RXFE,RXCA,RXMG,RXNA,RXKA,RXNB
C    2,RXN4*VLNH4(L,NY,NX)+RXNB*VLNHB(L,NY,NX)
C    2+RXHY+RXAL*3.0+RXFE*3.0+RXCA*2.0+RXMG*2.0+RXNA+RXKA
C    2,XN4Y*VLNH4(L,NY,NX)+XNBY*VLNHB(L,NY,NX)
C    2+XHYY+XALY+XFEY+XCAY+XMGY+XNAY+XKAY
C    2,CCEC,XCEC(L,NY,NX),XCAX,XTLQ,FX,XTL1,FY,A1(L,NY,NX),BKVLW 
C     ENDIF
C
C     NH4-NH3+H IN NON-BAND AND BAND SOIL ZONES
C
C     RNH4,RNHB=NH4-NH3+H dissociation in non-band,band (mol m-3 t-1)
C     DPN4=NH4 dissociation constant from PARAMETER
C     AN41,AN31,AN4B,AN3B=NH4,NH3 activity in non-band,band (mol m-3)
C     AN3Q,AN3BQ=NH3 equilibrium activity at DPN4 in non-band,band 
C        (mol m-3)
C     AHY1=H activity (mol m-3)
C     TSL=maximum reaction rate (mol m-3 t-1)
C
      IF(VOLWNH.GT.ZEROS2(NY,NX))THEN
      XMINN=FIONH*AN41
      XMINP=FIONH*AN31
      AN3Q=DPN4*AN41/AHY1
      RNH4=AMAX1(-TSZ,-XMINN
     2,(AMIN1(TSZ,XMINP,AN31-AN3Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24)THEN
C     IF(L.EQ.3)THEN
C     WRITE(*,27)'RNH4',I,J,NFZ,NX,NY,L
C    2,RNH4,AN41,AN31,AN3Q,AHY1,AN31*AHY1/AN41,DPN4
C    4,CHY1,XMINN,XMINP,TSZ,PH(L,NY,NX)
C    5,ZNH4S(L,NY,NX),ZNH3S(L,NY,NX)
27    FORMAT(A8,6I4,30E12.4)
C     ENDIF 
      ELSE
      RNH4=0.0
      ENDIF
      IF(VOLWNB.GT.ZEROS2(NY,NX))THEN
      XMINN=FIONH*AN4B
      XMINP=FIONH*AN3B
      AN3BQ=DPN4*AN4B/AHY1
      RNHB=AMAX1(-TSZ,-XMIN
     2,(AMIN1(TSZ,XMINP,AN3B-AN3BQ)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24)THEN
C     IF(L.EQ.3)THEN
C     WRITE(*,27)'RNHB',I,J,NFZ,NX,NY,L
C    3,RNHB,AN4B,AN3B,AN3BQ,AHY1,AN3B*AHY1/AN4B,DPN4
C    4,CHY1,XMINN,XMINP,TSZ,PH(L,NY,NX) 
C     ENDIF
      ELSE
      RNHB=0.0
      ENDIF
C
C     TOTAL ION TRANSFORMATIONS FOR ALL REACTIONS ABOVE
C        (all units in mol m-3 t-1)
C
C     RN4S,RN4B=net NH4 transformation in non-band,band
C     RN3S,RN3B=net NH3 transformation in non-band,band
C     RAL,RFE,RHY,RCA,RMG,RNA,RKA,ROH=net Al,Fe,H,Ca,Mg,Na,K,OH
C        transformation
C     RHP1,RHP2=net HPO4,H2PO4 transformation in non-band
C     RXH1,RXH2,RX1P,RX2P=net R-OH,R-OH2,R-HPO4,R-H2PO4 in non-band
C     RHB1,RHB2=net HPO4,H2PO4 transformation in band
C     RBH1,RBH2,RB1P,RB2P=net R-OH,R-OH2,R-HPO4,R-H2PO4 in band
C     BKVLW=soil mass:water (Mg m-3)
C
      RN4S=RNH4-RXN4*BKVLWH
      RN4B=RNHB-RXNB*BKVLWB 
      RN3S=-RNH4
      RN3B=-RNHB
      RHP1=-RH2P-RXH1P 
      RHP2=RH2P-RXH2P-RYH2P
     2-RPALPX-RPFEPX-RPCADX-2.0*RPCAMX-3.0*RPCAHX
      RHB1=-RH2B-RXH1B 
      RHB2=RH2B-RXH2B-RYH2B
     2-RPALBX-RPFEBX-RPCDBX-2.0*RPCMBX-3.0*RPCHBX
      RXH1=-RYH2P-RXH1P
      RXH2=-RXH2P
      RX1P=RXH1P
      RX2P=RXH2P+RYH2P
      RBH1=-RYH2B-RXH1B
      RBH2=-RXH2B
      RB1P=RXH1B
      RB2P=RXH2B+RYH2B
      RHY=-RXHY*BKVLW
      RAL=-RXAL*BKVLW
      RFE=-RXFE*BKVLW
      RCA=-RXCA*BKVLW 
      RMG=-RXMG*BKVLW 
      RNA=-RXNA*BKVLW 
      RKA=-RXKA*BKVLW
      ROH=0.0
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,26)'RALL',I,J,NFZ,NX,NY,L 
C    2,-RXAL*BKVLW,-RXFE*BKVLW,RXAL,RXFE,BKVLW 
C     WRITE(*,26)'RCAL',I,J,NFZ,NX,NY,L,CCA1,ACA1,RCA 
C    2,-RXCA*BKVLW,RAL*VOLW(L,NY,NX) 
C    2,-RPCADX,-RPCAMX,VLPO4(L,NY,NX)
C    3,-RPCDBX,-RPCMBX,VLPOB(L,NY,NX)
C    4,-5.0*RPCAHX*VLPO4(L,NY,NX),-5.0*RPCHBX*VLPOB(L,NY,NX)
26    FORMAT(A8,6I4,100E12.4)
C     ENDIF
C
C     RESET H2O=H+ + OH- AT PH SET IN SOIL FILE
C
C     AHY2,AOH2=H,OH activity (mol m-3)
C     A1=ion activity coefficient from ‘hour1.f’
C     RHHY,RHOH=H,OH equilibration (mol m-3 t-1)
C     PH=soil pH from soil file
C
      AHY2=10.0**(-PH(L,NY,NX))*1.0E+03*A1(L,NY,NX)
      AOH2=DPH2O/AHY2
      RHHY=AHY2-AHY1-RHY
      RHOH=AOH2-AOH1-ROH
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,1109)'RHHY',I,J,NFZ,L,AHY2,AOH2,AHY1,AOH1 
C    2,RHY,ROH,RHHY,RHOH
C    2,CHY1,COH1,A1(L,NY,NX) 
C    3,DPH2O,PH(L,NY,NX),ZHY(L,NY,NX),ZOH(L,NY,NX)
1109  FORMAT(A8,4I4,20E12.4)
C     ENDIF
C
C     CONVERT TOTAL ION TRANSFORMATIONS FROM CHANGES IN CONCENTRATION
C     TO CHANGES IN MASS FOR USE IN 'REDIST'
C        (change mol m-3 t-1 to mol t-1)
C
C     TRN4S,TRN4B=total NH4 transformation in non-band,band
C     TRN3S,TRN3B=total NH3 transformation in non-band,band
C     TRH1P,TRH2P=net HPO4,H2PO4 transformation in non-band
C     TRH1B,TRH2B=net HPO4,H2PO4 transformation in band
C     TRXN4,TRXNB=total NH4 adsorption in non-band,band
C     TRXH1,TRXH2,TRX1P,TRX2P
C        =total R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in non-band
C     TRBH1,TRBH2,TRB1P,TRB2P
C        =total R-OH,R-OH2,R-HPO4,R-H2PO4 adsorption in band
C     TRALPO,TRFEPO,TRCAPD,TRCAPH,TRCAPM
C        =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation 
C         in non-band
C     TRALPB,TRFEPB,TRCPDB,TRCPHB,TRCPMB
C        =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation 
C         in band
C     VOLW=soil water content (m3)
C     VOLWNH,VOLWNB=water volume in NH4 non-band,band (m3)
C     VOLWPO,VOLWPB=water volume in H2PO4 non-band,band (m3)
C
      TRN4S(L,NY,NX)=TRN4S(L,NY,NX)+RN4S*VOLWNH
      TRN4B(L,NY,NX)=TRN4B(L,NY,NX)+RN4B*VOLWNB
      TRN3S(L,NY,NX)=TRN3S(L,NY,NX)+RN3S*VOLWNH
      TRN3B(L,NY,NX)=TRN3B(L,NY,NX)+RN3B*VOLWNB
      TRHY(L,NY,NX)=TRHY(L,NY,NX)+(RHY+RHHY)*VOLW(L,NY,NX) 
      TROH(L,NY,NX)=TROH(L,NY,NX)+(ROH+RHOH)*VOLW(L,NY,NX) 
      TRAL(L,NY,NX)=TRAL(L,NY,NX)+RAL*VOLW(L,NY,NX)
      TRFE(L,NY,NX)=TRFE(L,NY,NX)+RFE*VOLW(L,NY,NX)
      TRCA(L,NY,NX)=TRCA(L,NY,NX)+RCA*VOLW(L,NY,NX)
      TRMG(L,NY,NX)=TRMG(L,NY,NX)+RMG*VOLW(L,NY,NX)
      TRNA(L,NY,NX)=TRNA(L,NY,NX)+RNA*VOLW(L,NY,NX)
      TRKA(L,NY,NX)=TRKA(L,NY,NX)+RKA*VOLW(L,NY,NX)
      TBH2O(L,NY,NX)=TBH2O(L,NY,NX)+RHHY*VOLW(L,NY,NX) 
      TRH1P(L,NY,NX)=TRH1P(L,NY,NX)+RHP1*VOLWPO
      TRH2P(L,NY,NX)=TRH2P(L,NY,NX)+RHP2*VOLWPO
      TRH1B(L,NY,NX)=TRH1B(L,NY,NX)+RHB1*VOLWPB
      TRH2B(L,NY,NX)=TRH2B(L,NY,NX)+RHB2*VOLWPB
      TRXN4(L,NY,NX)=TRXN4(L,NY,NX)+RXN4*BKVLNH
      TRXNB(L,NY,NX)=TRXNB(L,NY,NX)+RXNB*BKVLNB
      TRXHY(L,NY,NX)=TRXHY(L,NY,NX)+RXHY*BKVLX 
      TRXAL(L,NY,NX)=TRXAL(L,NY,NX)+RXAL*BKVLX 
      TRXFE(L,NY,NX)=TRXFE(L,NY,NX)+RXFE*BKVLX
      TRXCA(L,NY,NX)=TRXCA(L,NY,NX)+RXCA*BKVLX
      TRXMG(L,NY,NX)=TRXMG(L,NY,NX)+RXMG*BKVLX
      TRXNA(L,NY,NX)=TRXNA(L,NY,NX)+RXNA*BKVLX
      TRXKA(L,NY,NX)=TRXKA(L,NY,NX)+RXKA*BKVLX
      TRXH1(L,NY,NX)=TRXH1(L,NY,NX)+RXH1*VOLWPO
      TRXH2(L,NY,NX)=TRXH2(L,NY,NX)+RXH2*VOLWPO
      TRX1P(L,NY,NX)=TRX1P(L,NY,NX)+RX1P*VOLWPO
      TRX2P(L,NY,NX)=TRX2P(L,NY,NX)+RX2P*VOLWPO
      TRBH1(L,NY,NX)=TRBH1(L,NY,NX)+RBH1*VOLWPB
      TRBH2(L,NY,NX)=TRBH2(L,NY,NX)+RBH2*VOLWPB
      TRB1P(L,NY,NX)=TRB1P(L,NY,NX)+RB1P*VOLWPB
      TRB2P(L,NY,NX)=TRB2P(L,NY,NX)+RB2P*VOLWPB
      TRALPO(L,NY,NX)=TRALPO(L,NY,NX)+RPALPX*VOLWPO
      TRFEPO(L,NY,NX)=TRFEPO(L,NY,NX)+RPFEPX*VOLWPO
      TRCAPD(L,NY,NX)=TRCAPD(L,NY,NX)+RPCADX*VOLWPO
      TRCAPH(L,NY,NX)=TRCAPH(L,NY,NX)+RPCAHX*VOLWPO
      TRCAPM(L,NY,NX)=TRCAPM(L,NY,NX)+RPCAMX*VOLWPO
      TRALPB(L,NY,NX)=TRALPB(L,NY,NX)+RPALBX*VOLWPB
      TRFEPB(L,NY,NX)=TRFEPB(L,NY,NX)+RPFEBX*VOLWPB
      TRCPDB(L,NY,NX)=TRCPDB(L,NY,NX)+RPCDBX*VOLWPB
      TRCPHB(L,NY,NX)=TRCPHB(L,NY,NX)+RPCHBX*VOLWPB
      TRCPMB(L,NY,NX)=TRCPMB(L,NY,NX)+RPCMBX*VOLWPB
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH 
C    2.AND.L.EQ.2)THEN
C     WRITE(*,24)'AHY1',I,J,NFZ,L,AHY1,AOH1,AHY1*AOH1
C    2,DPH2O,CHY1,COH1,CHY1*COH1,DPH2O
C    3,PH(L,NY,NX)
C     WRITE(*,24)'RN4S',I,J,NFZ,L,RN4S,RN3S,RNH4,RXN4,VOLWNH 
C     WRITE(*,24)'RHP1',I,J,NFZ,L,RHP1,RH2P,RXH1P
C    2,TRX1P(L,NY,NX),TRH2P(L,NY,NX)
C     WRITE(*,24)'TRAL',I,J,NFZ,L,TRAL(L,NY,NX),RAL,-RXAL*BKVLW 
C    2,-RPALPX*VLPO4(L,NY,NX),-RPALBX*VLPOB(L,NY,NX),VOLW(L,NY,NX)
C    3,ZAL(L,NY,NX)
C     WRITE(*,24)'TRFE',I,J,NFZ,L,TRFE(L,NY,NX),RFE,-RXFE*BKVLW 
C    2,-RPFEPX*VLPO4(L,NY,NX),-RPFEBX*VLPOB(L,NY,NX),VOLW(L,NY,NX)
C    3,ZFE(L,NY,NX)
C     WRITE(*,24)'RHP2',I,J,NFZ,L,CH2P1,RHP2,RH2P,-RXH2P,-RYH2P
C    2,-RPALPX,-RPFEPX,-RPCADX,-2.0*RPCAMX,-3.0*RPCAHX 
C     WRITE(*,24)'RCA',I,J,NFZ,L,TRCA(L,NY,NX),RCA*VOLW(L,NY,NX)
C    2,RCA,RXCA,RPCADX,RPCAMX,VLPO4(L,NY,NX)
C    3,RPCDBX,RPCMBX,VLPOB(L,NY,NX),5.0*RPCAHX,5.0*RPCHBX
24    FORMAT(A8,4I4,60E12.4)
C     ENDIF
      ENDIF
C
C     END SOLUTE REACTIONS
C
C     CHANGE IN WIDTHS AND DEPTHS OF FERTILIZER BANDS FROM
C     VERTICAL AND HORIZONTAL DIFFUSION DRIVEN BY CONCENTRATION
C     DIFFERENCES BETWEEN BAND AND NON-BAND SOIL ZONES
C
C     NH4 FERTILIZER BAND
C
C     IFNHB=banded NH4 fertilizer flag (0=no band,1=band)
C     ROWN=NH4 fertilizer band row width from fertilizer file (m)
C     CDPTH=depth to bottom of soil layer (m)
C     DPNH4=NH4 fertilizer band depth (m) 
C
      IF(IFNHB(NY,NX).EQ.1.AND.ROWN(NY,NX).GT.0.0)THEN
C
C     NH4 BAND WIDTH
C
C     DWNH4=change in NH4 fertilizer band width (m t-1)
C     WDNHB=layer NH4 fertilizer band width (m)
C     ZNSGL=NH4 diffusivity from ‘hour1.f’ (m2 h-1)
C     TORT=tortuosity from ‘watsub.f’
C     XNFH=time step for solute fluxes from ‘wthr.f’ (h t-1)
C
C     TRANSFER NH4, NH3 FROM NON-BAND TO BAND
C     DURING BAND GROWTH
C
C     DNH4S,DNH3S,DXNH4=transfer of NH4,NH3,exchangeable NH4 
C        (mol t-1)
C     FVLNH4=relative change in VLNH4 (t-1)
C     TRN4S,TRN4B=total NH4 transformation in non-band,band 
C        (mol t-1) 
C     TRN3S,TRN3B=total NH3 transformation in non-band,band 
C        (mol t-1) 
C     TRXN4,TRXNB=total NH4 adsorption in non-band,band
C        (mol t-1) 
C     ZNH4FA,ZNH3FA,ZNHUF =broadcast NH4,NH3,urea fertilizer (mol N)
C     ZNH4FB,ZNH3FB,ZNHUFB=banded NH4,NH3,urea fertilizer (mol N)
C
      DNH4S=FVLNH4(L,NY,NX)*ZNH4S(L,NY,NX)/14.0
      DNH3S=FVLNH4(L,NY,NX)*ZNH3S(L,NY,NX)/14.0
      DXNH4=FVLNH4(L,NY,NX)*XN4(L,NY,NX)
      TRN4S(L,NY,NX)=TRN4S(L,NY,NX)+DNH4S
      TRN4B(L,NY,NX)=TRN4B(L,NY,NX)-DNH4S
      TRN3S(L,NY,NX)=TRN3S(L,NY,NX)+DNH3S
      TRN3B(L,NY,NX)=TRN3B(L,NY,NX)-DNH3S
      TRXN4(L,NY,NX)=TRXN4(L,NY,NX)+DXNH4
      TRXNB(L,NY,NX)=TRXNB(L,NY,NX)-DXNH4
      DNH4FA=FVLNH4(L,NY,NX)*ZNH4FA(L,NY,NX)
      DNH3FA=FVLNH4(L,NY,NX)*ZNH3FA(L,NY,NX)
      DNHUFA=FVLNH4(L,NY,NX)*ZNHUFA(L,NY,NX)
      ZNH4FA(L,NY,NX)=ZNH4FA(L,NY,NX)+DNH4FA 
      ZNH3FA(L,NY,NX)=ZNH3FA(L,NY,NX)+DNH3FA 
      ZNHUFA(L,NY,NX)=ZNHUFA(L,NY,NX)+DNHUFA 
      ZNH4FB(L,NY,NX)=ZNH4FB(L,NY,NX)-DNH4FA 
      ZNH3FB(L,NY,NX)=ZNH3FB(L,NY,NX)-DNH3FA 
      ZNHUFB(L,NY,NX)=ZNHUFB(L,NY,NX)-DNHUFA 
      ENDIF
C
C     NO3 FERTILIZER BAND
C
C     IFNOB=banded NO3 fertilizer flag (m)
C     ROWO=NO3 fertilizer band row width from fertilizer file (m)
C     CDPTH=depth to bottom of soil layer (m)
C     DPNO3=NO3 fertilizer band depth (m) 
C
      IF(IFNOB(NY,NX).EQ.1.AND.ROWO(NY,NX).GT.0.0)THEN
C
C     TRANSFER NO3 FROM NON-BAND TO BAND
C     DURING BAND GROWTH
C
C     DNO3S,DNO2S=transfer of NO3,NO2 (mol t-1)
C     FVLNO3=relative change in VLNO3 (t-1)
C     TRNO3,TRNOB=NO3 dissolution in non-band,band (mol t-1)
C     TRNO2,TRN2B=NO2 dissolution in non-band,band (mol t-1)
C     ZNO3FA,ZNO3FB=broadcast,banded NO3 fertilizer (mol)
C
      DNO3S=FVLNO3(L,NY,NX)*ZNO3S(L,NY,NX)/14.0
      DNO2S=FVLNO3(L,NY,NX)*ZNO2S(L,NY,NX)/14.0
      TRNO3(L,NY,NX)=TRNO3(L,NY,NX)+DNO3S
      TRNO2(L,NY,NX)=TRNO2(L,NY,NX)+DNO2S
      TRNOB(L,NY,NX)=TRNOB(L,NY,NX)-DNO3S
      TRN2B(L,NY,NX)=TRN2B(L,NY,NX)-DNO2S
      DNO3FA=FVLNO3(L,NY,NX)*ZNO3FA(L,NY,NX)
      ZNO3FA(L,NY,NX)=ZNO3FA(L,NY,NX)+DNO3FA 
      ZNO3FB(L,NY,NX)=ZNO3FB(L,NY,NX)-DNO3FA 
      ENDIF
C
C     PO4 FERTILIZER BAND
C
C     IFPOB=banded H2PO4 fertilizer flag (m)
C     ROWP=H2PO4 fertilizer band row width from fertilizer file (m)
C     CDPTH=depth to bottom of soil layer (m)
C     DPPO4=H2PO4 fertilizer band depth (m) 
C
      IF(IFPOB(NY,NX).EQ.1.AND.ROWP(NY,NX).GT.0.0)THEN
      DZH1P=FVLPO4(L,NY,NX)*H1PO4(L,NY,NX)/31.0
      DZH2P=FVLPO4(L,NY,NX)*H2PO4(L,NY,NX)/31.0
      DXOH0=FVLPO4(L,NY,NX)*XOH0(L,NY,NX)
      DXOH1=FVLPO4(L,NY,NX)*XOH1(L,NY,NX)
      DXOH2=FVLPO4(L,NY,NX)*XOH2(L,NY,NX)
      DXH1P=FVLPO4(L,NY,NX)*XH1P(L,NY,NX)
      DXH2P=FVLPO4(L,NY,NX)*XH2P(L,NY,NX)
      DPALP=FVLPO4(L,NY,NX)*PALPO(L,NY,NX)
      DPFEP=FVLPO4(L,NY,NX)*PFEPO(L,NY,NX)
      DPCDP=FVLPO4(L,NY,NX)*PCAPD(L,NY,NX)
      DPCHP=FVLPO4(L,NY,NX)*PCAPH(L,NY,NX)
      DPCMP=FVLPO4(L,NY,NX)*PCAPM(L,NY,NX)
      TRH1P(L,NY,NX)=TRH1P(L,NY,NX)+DZH1P
      TRH2P(L,NY,NX)=TRH2P(L,NY,NX)+DZH2P
      TRH1B(L,NY,NX)=TRH1B(L,NY,NX)-DZH1P
      TRH2B(L,NY,NX)=TRH2B(L,NY,NX)-DZH2P
      TRXH0(L,NY,NX)=TRXH0(L,NY,NX)+DXOH0
      TRXH1(L,NY,NX)=TRXH1(L,NY,NX)+DXOH1
      TRXH2(L,NY,NX)=TRXH2(L,NY,NX)+DXOH2
      TRX1P(L,NY,NX)=TRX1P(L,NY,NX)+DXH1P
      TRX2P(L,NY,NX)=TRX2P(L,NY,NX)+DXH2P
      TRBH0(L,NY,NX)=TRBH0(L,NY,NX)-DXOH0
      TRBH1(L,NY,NX)=TRBH1(L,NY,NX)-DXOH1
      TRBH2(L,NY,NX)=TRBH2(L,NY,NX)-DXOH2
      TRB1P(L,NY,NX)=TRB1P(L,NY,NX)-DXH1P
      TRB2P(L,NY,NX)=TRB2P(L,NY,NX)-DXH2P
      TRALPO(L,NY,NX)=TRALPO(L,NY,NX)+DPALP
      TRFEPO(L,NY,NX)=TRFEPO(L,NY,NX)+DPFEP
      TRCAPD(L,NY,NX)=TRCAPD(L,NY,NX)+DPCDP
      TRCAPH(L,NY,NX)=TRCAPH(L,NY,NX)+DPCHP
      TRCAPM(L,NY,NX)=TRCAPM(L,NY,NX)+DPCMP
      TRALPB(L,NY,NX)=TRALPB(L,NY,NX)-DPALP
      TRFEPB(L,NY,NX)=TRFEPB(L,NY,NX)-DPFEP
      TRCPDB(L,NY,NX)=TRCPDB(L,NY,NX)-DPCDP
      TRCPHB(L,NY,NX)=TRCPHB(L,NY,NX)-DPCHP
      TRCPMB(L,NY,NX)=TRCPMB(L,NY,NX)-DPCMP
C
C     TRANSFER HPO4,H2PO4 FROM NON-BAND TO BAND
C     DURING BAND GROWTH DEPENDING ON SALT
C     VS. NON-SALT OPTION
C
C     D*=transfer of solute,adsorbed,precipitated HPO4,H2PO4 (mol t-1)
C     FVLPO4=relative change in VLPO4 (t-1)
C     TR*=net P transformation (mol t-1) 
C
      IF(ISALTG.NE.0)THEN
      DZH0P=FVLPO4(L,NY,NX)*H0PO4(L,NY,NX)
      DZH3P=FVLPO4(L,NY,NX)*H3PO4(L,NY,NX)
      DZF1P=FVLPO4(L,NY,NX)*ZFE1P(L,NY,NX)
      DZF2P=FVLPO4(L,NY,NX)*ZFE2P(L,NY,NX)
      DZC0P=FVLPO4(L,NY,NX)*ZCA0P(L,NY,NX)
      DZC1P=FVLPO4(L,NY,NX)*ZCA1P(L,NY,NX)
      DZC2P=FVLPO4(L,NY,NX)*ZCA2P(L,NY,NX)
      DZM1P=FVLPO4(L,NY,NX)*ZMG1P(L,NY,NX)
      TRH0P(L,NY,NX)=TRH0P(L,NY,NX)+DZH0P
      TRH3P(L,NY,NX)=TRH3P(L,NY,NX)+DZH3P
      TRF1P(L,NY,NX)=TRF1P(L,NY,NX)+DZF1P
      TRF2P(L,NY,NX)=TRF2P(L,NY,NX)+DZF2P
      TRC0P(L,NY,NX)=TRC0P(L,NY,NX)+DZC0P
      TRC1P(L,NY,NX)=TRC1P(L,NY,NX)+DZC1P
      TRC2P(L,NY,NX)=TRC2P(L,NY,NX)+DZC2P
      TRM1P(L,NY,NX)=TRM1P(L,NY,NX)+DZM1P
      TRH0B(L,NY,NX)=TRH0B(L,NY,NX)-DZH0P
      TRH3B(L,NY,NX)=TRH3B(L,NY,NX)-DZH3P
      TRF1B(L,NY,NX)=TRF1B(L,NY,NX)-DZF1P
      TRF2B(L,NY,NX)=TRF2B(L,NY,NX)-DZF2P
      TRC0B(L,NY,NX)=TRC0B(L,NY,NX)-DZC0P
      TRC1B(L,NY,NX)=TRC1B(L,NY,NX)-DZC1P
      TRC2B(L,NY,NX)=TRC2B(L,NY,NX)-DZC2P
      TRM1B(L,NY,NX)=TRM1B(L,NY,NX)-DZM1P
      ENDIF
      ENDIF
C     ENDIF
C
C     SUBTRACT FERTILIZER DISSOLUTION FROM FERTILIZER POOLS
C
C     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=broadcast NH4,NH3,urea,NO3
C        fertilizer (mol)
C     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 
C        fertilizer (mol)
C     RSN4AA,RSN4BA=rate of broadcast NH4 fertilizer dissolution 
C        in non-band,band (mol t-1)
C     RSN3AA,RSN3BA=rate of broadcast NH3 fertilizer dissolution 
C        in non-band,band (mol t-1)
C     RSNUAA,RSNUBA=rate of broadcast urea fertilizer dissolution 
C        in non-band,band (mol t-1)
C     RSNOAA,RSNOBA=rate of broadcast NO3 fertilizer dissolution 
C        in non-band,band (mol t-1)
C     RSN4BB=rate of banded NH4 fertilizer dissolution in band 
C        (mol t-1)
C     RSN3BB=rate of banded NH3 fertilizer dissolution in band
C        (mol t-1)
C     RSNUBB=rate of banded urea fertilizer dissolution in band
C        (mol t-1)
C     RSNOBB=rate of banded NO3 fertilizer dissolution in band
C        (mol t-1)
C
      ZNH4FA(L,NY,NX)=ZNH4FA(L,NY,NX)-RSN4AA-RSN4BA
      ZNH3FA(L,NY,NX)=ZNH3FA(L,NY,NX)-RSN3AA-RSN3BA
      ZNHUFA(L,NY,NX)=ZNHUFA(L,NY,NX)-RSNUAA-RSNUBA
      ZNO3FA(L,NY,NX)=ZNO3FA(L,NY,NX)-RSNOAA-RSNOBA
      ZNH4FB(L,NY,NX)=ZNH4FB(L,NY,NX)-RSN4BB
      ZNH3FB(L,NY,NX)=ZNH3FB(L,NY,NX)-RSN3BB
      ZNHUFB(L,NY,NX)=ZNHUFB(L,NY,NX)-RSNUBB
      ZNO3FB(L,NY,NX)=ZNO3FB(L,NY,NX)-RSNOBB
C
C     ADD FERTILIZER DISSOLUTION TO ION TRANSFORMATIONS 
C     AND CONVERT TO MASS FOR USE IN ‘REDIST.F’
C
C     TRN3G=NH3 dissolution (g N t-1)
C     TRN4S,TRN4B=NH4 dissolution in non-band,band (g N t-1) 
C     TRN3S,TRN3B=NH3 dissolution from urea in non-band,band (g N t-1)
C     TRNO3,TRNOB=NO3 dissolution in non-band,band (g N t-1)
C     TRH1P,TRH1B=HPO4 dissolution in non-band,band (g P t-1)
C     TRH2P,TRH2B=H2PO4 dissolution in non-band,band (g P t-1)
C
      TRN3G(L,NY,NX)=TRN3G(L,NY,NX)+RSN3AA+RSN3BA+RSN3BB 
      TRN4S(L,NY,NX)=TRN4S(L,NY,NX)+RSN4AA 
      TRN4B(L,NY,NX)=TRN4B(L,NY,NX)+RSN4BA+RSN4BB 
      TRN3S(L,NY,NX)=TRN3S(L,NY,NX)+RSNUAA
      TRN3B(L,NY,NX)=TRN3B(L,NY,NX)+RSNUBA+RSNUBB
      TRNO3(L,NY,NX)=TRNO3(L,NY,NX)+RSNOAA
      TRNOB(L,NY,NX)=TRNOB(L,NY,NX)+RSNOBA+RSNOBB
C
C     TRCO2=net CO2 transformation (g C t-1)
C     TRN3G=net NH3 gas transformation (g N t-1)
C     TRN4S,TRN3S,TRNO2,TRNO3=net NH4,NH3,NO2,NO3 transformation 
C        in non-band(g N t-1) 
C     TRN4B,TRN3B,TRN2B,TRNOB=net NH4,NH3,NO2,NO3 transformation 
C        in band (g N t-1) 
C     TRH1P,TRH2P=net HPO4,H2PO4 transformation in non-band (g P t-1) 
C     TRH1B,TRH2B=net HPO4,H2PO4 transformation in band (g P t-1) 
C
      TRCO2(L,NY,NX)=TRCO2(L,NY,NX)*12.0
      TRN3G(L,NY,NX)=TRN3G(L,NY,NX)*14.0
      TRN4S(L,NY,NX)=TRN4S(L,NY,NX)*14.0
      TRN4B(L,NY,NX)=TRN4B(L,NY,NX)*14.0
      TRN3S(L,NY,NX)=TRN3S(L,NY,NX)*14.0
      TRN3B(L,NY,NX)=TRN3B(L,NY,NX)*14.0
      TRNO3(L,NY,NX)=TRNO3(L,NY,NX)*14.0
      TRNOB(L,NY,NX)=TRNOB(L,NY,NX)*14.0
      TRNO2(L,NY,NX)=TRNO2(L,NY,NX)*14.0
      TRN2B(L,NY,NX)=TRN2B(L,NY,NX)*14.0
      TRH1P(L,NY,NX)=TRH1P(L,NY,NX)*31.0
      TRH2P(L,NY,NX)=TRH2P(L,NY,NX)*31.0
      TRH1B(L,NY,NX)=TRH1B(L,NY,NX)*31.0
      TRH2B(L,NY,NX)=TRH2B(L,NY,NX)*31.0
      ENDIF
C     IF(I.EQ.116)THEN
C     WRITE(*,9984)'TRN3S',I,J,L,TRN4S(L,NY,NX),TRN3S(L,NY,NX)
9984  FORMAT(A8,3I4,20F14.7)
C     ENDIF
9985  CONTINUE 
C
C     SURFACE LITTER
C
C     TIME STEPS
C
C     TPDH,TADAH,TADCH,TSLH,TRWH=maximum
C        precipitation,adsorption,solute,silicate weathwering
C        rates from PARAMETER (mol m-3 h-1)
C     TPD,TADA,TADC,TSL,TRW=maximum
C        precipitation,adsorption,solute,silicate weathwering
C        rates (mol m-3 t-1)
C     XNFH=time step for solute transformations from ‘wthr.f’ (h t-1)
C
      TPD=TPDH*XNFH 
      TPZ=TPZH*XNFH 
      TADA=TADAH*XNFH 
      TADC=TADCH*XNFH 
      TSL=TSLH*XNFH 
      TSZ=TSZH*XNFH 
      TRW=TRWH*XNFH 
      TPA=TRW 
C
      IF(VOLW(0,NY,NX).GT.ZEROS2(NY,NX))THEN
C
C     BKVL=litter mass (Mg DM)
C     VOLW(0,THETW(0,=litter water volume, concentration (m3,m3 m-3)
C     CCEC0=litter CEC (mol Mg-1)
C     BKVLW=soil mass:water (Mg DM m-3)
C     ORGC=litter mass (g C)
C     COOH=carboxyl sites (mol Mg C-1)
C
      THETW0=THETW(0,NY,NX)
      BKVLX=BKVL(0,NY,NX)
      BKVLW=BKVLX/VOLW(0,NY,NX)
      BKVLN=BKVL(NU(NY,NX),NY,NX)
      IF(BKVLX.GT.ZEROS2(NY,NX))THEN
      CCEC0=AMAX1(ZEROC,COOH*1.0E-06*ORGC(0,NY,NX)/BKVLX)
      ELSE
      CCEC0=ZERO
      ENDIF
C
C     UREA HYDROLYSIS IN SURFACE RESIDUE
C
C     VOLQ=biologically active litter water volume from ‘nitro.f’ (m3)
C     COQCK=aqueous concentration of microbial activity (g C m-3 h-1)
C     TOQCK=total active biomass respiration activity from ‘nitro.f’
C        (g C t-1)
C     DUKD,DUKM=effective,minimum Km for urea hydrolysis from PARAMETER
C        (mol Mg-1) 
C     DUKI=Ki for microbial activity effects on urea hydrolysis as in
C        ‘nitro.f’(SOIL SCI 136:56) from PARAMETER (g C m-3 h-1)
C     XNFH=time step for solute transformations from ‘wthr.f’ (h t-1)
C
      IF(VOLQ(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      COQCK=AMIN1(0.1E+06,TOQCK(0,NY,NX)/(VOLQ(0,NY,NX)*XNFH))
      ELSE
      COQCK=0.1E+06
      ENDIF
      DUKD=DUKM*(1.0+COQCK/DUKI)
C
C     UREA HYDROLYSIS INHIBITION
C
C     ZNHU0,ZNHUI=initial,current urea hydrolysis inhibition activity
C        from ‘hour1.f’
C     RNHUX=rate constant for decline in urea hydrolysis inhibition
C        (t-1)
C
      IF(ZNHU0(0,NY,NX).GT.ZEROS(NY,NX)
     2.AND.ZNHUI(0,NY,NX).GT.ZEROS(NY,NX))THEN
      ZNHUI(0,NY,NX)=ZNHUI(0,NY,NX)
     2-RNHUX*ZNHUI(0,NY,NX)
     3*AMAX1(RNHUX,1.0-ZNHUI(0,NY,NX)/ZNHU0(0,NY,NX))
      ELSE
      ZNHUI(0,NY,NX)=0.0
      ENDIF
C
C     UREA CONCENTRATION AND HYDROLYSIS IN SURFACE RESIDUE
C
C     ZNHUFA=urea fertilizer (mol N) 
C     BKVL,VOLW=litter mass,water content (Mg,m3)
C     CNHUA=concentration of urea fertilizer (mol N Mg-1) 
C     DFNSA=effect of microbial concentration on urea hydrolysis 
C     DUKD=effective Km for urea hydrolysis (mol N Mg-1) 
C     RSNUA=rate of urea hydrolysis (mol N t-1) 
C     SPNHU=specific urea hydrolysis rate for microbial activity 
C        from PARAMETER (mol N g C-1)
C     TOQCK=total active biomass respiration activity from ‘nitro.f’
C        (g C t-1)
C     TFNQ=temperature effect on microbial activity from ‘nitro.f’
C     ZNHUI=current inhibition activity
C
      IF(ZNHUFA(0,NY,NX).GT.ZEROS(NY,NX)
     2.AND.BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
      CNHUA=ZNHUFA(0,NY,NX)/BKVL(0,NY,NX)
      DFNSA=CNHUA/(CNHUA+DUKD)
      RSNUA=AMIN1(ZNHUFA(0,NY,NX)
     2,SPNHU*TOQCK(0,NY,NX)*DFNSA*TFNQ(0,NY,NX)*(1.0-ZNHUI(0,NY,NX)))
      ELSE
      RSNUA=0.0
      ENDIF
C     IF(ZNHUFA(0,NY,NX).GT.ZEROS(NY,NX))THEN
C     WRITE(*,8778)'UREA0',I,J,NFZ,IUTYP(NY,NX)
C    2,ZNHUFA(0,NY,NX),RSNUA
C    2,DFNSA,TFNQ(0,NY,NX),CNHUA,DUKD,DUKM,DUKI,TOQCK(0,NY,NX)
C    3,BKVL(0,NY,NX),TFNQ(0,NY,NX),SPNHU,ZNHU0(0,NY,NX),ZNHUI(0,NY,NX)
C    4,RNHUX 
8778  FORMAT(A8,4I4,40E12.4)
C     ENDIF
C
C     NH4, NH3, UREA, NO3 DISSOLUTION IN SURFACE RESIDUE
C     FROM FIRST-ORDER FUNCTIONS OF REMAINING
C     FERTILIZER (NOTE: SUPERPHOSPHATE AND ROCK PHOSPHATE
C     ARE REPRESENTED AS MONOCALCIUM PHOSPHATE AND HYDROXYAPATITE
C     MODELLED IN PHOSPHORUS REACTIONS BELOW)
C
C     RSN4AA=rate of broadcast NH4 fertilizer dissolution (mol N t-1) 
C     RSN3AA=rate of broadcast NH3 fertilizer dissolution (mol N t-1) 
C     RSNUAA=rate of broadcast urea fertilizer dissolution (mol N t-1) 
C     RSNOAA=rate of broadcast NO3 fertilizer dissolution (mol N t-1)  
C     SPNH4,SPNH3,SPNO3 =specific rate constants for 
C        NH4,NH3,NO3 fertilizer dissolution from PARAMETER (t-1)
C     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=broadcast NH4,NH3,urea,NO3
C        fertilizer from ‘hour1.f’ (mol N)
C     THETW0,=litter water concentration (m3 m-3)
C
      RSN4AA=SPNH4*ZNH4FA(0,NY,NX)*THETW0
      RSN3AA=SPNH3*ZNH3FA(0,NY,NX)
      RSNUAA=RSNUA
      RSNOAA=SPNO3*ZNO3FA(0,NY,NX)*THETW0
C
C     SOLUBLE AND EXCHANGEABLE NH4 CONCENTRATIONS
C
C     VOLW=litter water volume (m3)
C     RN4X,RN3X=NH4,NH3 input from uptake,mineralization,dissolution 
C        (g N t-1) 
C     XNH4S=net change in NH4 from ‘nitro.f’ (g N t-1)
C     RSN4AA=rate of broadcast NH4 fertilizer dissolution (mol N t-1) 
C     RSNUAA=rate of broadcast urea fertilizer dissolution (mol N t-1) 
C     CN41,CN31=total NH4,NH3 concentration (mol N m-3) 
C     XN41=adsorbed NH4 concentration (mol N Mg-1) 
C
      VOLWMX=14.0*VOLW(0,NY,NX)
      RN4X=(XNH4S(0,NY,NX)+14.0*RSN4AA) 
      RN3X=14.0*RSNUAA 
      CN41=AMAX1(ZEROC,ZNH4S(0,NY,NX)+RN4X)/VOLWMX 
      CN31=AMAX1(ZEROC,ZNH3S(0,NY,NX)+RN3X)/VOLWMX 
      IF(BKVLX.GT.ZEROS2(NY,NX))THEN
      XN41=AMAX1(ZEROC,XN4(0,NY,NX)/BKVLX)
      ELSE
      XN41=0.0
      ENDIF
C
C     SOLUBLE, EXCHANGEABLE AND PRECIPITATED PO4 CONCENTRATIONS 
C
C     VOLW=litter water volume (m3) 
C     RH1PX,RH2PX=HPO4,H2PO4 inputs from mineralization, uptake 
C        (mol P t-1)
C     XH1PS,XH2PS=net change in HPO4,H2PO4 from ‘nitro.f’ (g P t-1) 
C     CH1P1,CH2P1=HPO4,H2PO4 concentrations (mol P m-3) 
C
      VOLWMP=31.0*VOLW(0,NY,NX)
      RH1PX=XH1PS(0,NY,NX)/VOLWMP
      RH2PX=XH2PS(0,NY,NX)/VOLWMP
      CH1P1=AMAX1(ZEROC,H1PO4(0,NY,NX)/VOLWMP+RH1PX)
      CH2P1=AMAX1(ZEROC,H2PO4(0,NY,NX)/VOLWMP+RH2PX)
C
C     PHOSPHORUS TRANSFORMATIONS IN SURFACE RESIDUE
C     AND EXCHANGE WITH PHOSPHORUS IN SURFACE LAYER (NU)
C
C     PALOH1,PFEOH1,PCACO1,PCASO1=precipitated
C         AL(OH)3,FE(OH)3,CACO3,CASO4 vs water (mol m-3)
C     PALPO1,PFEPO1=precipitated AlPO4,FEPO4 vs water (mol m-3) 
C     PCAPM1,PCAPD1,PCAPH1=precipitated CaH2PO4,CaHPO4,apatite vs water
C        (mol M-3) 
C
      IF(VOLW(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      PALOH1=AMAX1(0.0,PALOH(0,NY,NX))/VOLW(0,NY,NX)
      PFEOH1=AMAX1(0.0,PFEOH(0,NY,NX))/VOLW(0,NY,NX)
      PCACO1=AMAX1(0.0,PCACO(0,NY,NX))/VOLW(0,NY,NX)
      PCASO1=AMAX1(0.0,PCASO(0,NY,NX))/VOLW(0,NY,NX)
      PALPO1=AMAX1(0.0,PALPO(0,NY,NX))/VOLW(0,NY,NX)
      PFEPO1=AMAX1(0.0,PFEPO(0,NY,NX))/VOLW(0,NY,NX)
      PCAPM1=AMAX1(0.0,PCAPM(0,NY,NX))/VOLW(0,NY,NX)
      PCAPD1=AMAX1(0.0,PCAPD(0,NY,NX))/VOLW(0,NY,NX)
      PCAPH1=AMAX1(0.0,PCAPH(0,NY,NX))/VOLW(0,NY,NX)
      ELSE
      PALOH1=0.0
      PFEOH1=0.0
      PCACO1=0.0
      PCASO1=0.0
      PALPO1=0.0
      PFEPO1=0.0
      PCAPM1=0.0
      PCAPD1=0.0
      PCAPH1=0.0
      ENDIF
C
C     IF SALT OPTION SELECTED IN SITE FILE
C     THEN SOLVE EQUILIBRIA REACTIONS IN LITTER
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C
      IF(ISALTG.NE.0)THEN
C
C     PRECIPITATION-DISSOLUTION REACTIONS FOR H2PO4 AND CO-PRECIPITATES
C     CALCULATED FROM ACTIVITIES OF REACTANTS AND PRODUCTS THROUGH 
C     SOLUTIONS FOR THEIR EQUILIBRIUM CONCENTRATIONS 
C
C     VOLW=litter water content (m3)
C     XZHYS=total H+ production from ‘nitro.f’ and ‘grosub.f’
C        (mol t-1)
C     C*1,A*1=litter ion concentration,activity (mol m-3)
C     Z*=litter ion content (mol)
C        ion code:HY=H+,OH=OH-,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+
C            :NA=Na+,KA=K+,SO4=SO42-,CL=Cl-,CO2=CO2,CO3=CO3
C            :H1P=HPO42-,H2P=H2PO4-,N4=NH4+,N3=NH3 
C     DP*,SP*=dissociation,solubility products from PARAMETER above
C     SPX=equilibrium product concentration (mol m-3)
C     R*X=precipitation(+ve) or dissolution (-ve) rate 
C        (mol m-3 t-1)
C     TPD=maximum precipitation rate (mol m-3 t-1)
C
      RZHYS=XZHYS(0,NY,NX)/VOLW(0,NY,NX)
      CHY1=ZHY(0,NY,NX)/VOLW(0,NY,NX)
      COH1=ZOH(0,NY,NX)/VOLW(0,NY,NX)
      CAL1=AMAX1(ZEROC,ZAL(0,NY,NX)/VOLW(0,NY,NX))
      CFE1=AMAX1(ZEROC,ZFE(0,NY,NX)/VOLW(0,NY,NX))
      CCA1=AMAX1(ZEROC,ZCA(0,NY,NX)/VOLW(0,NY,NX))
      CMG1=AMAX1(ZEROC,ZMG(0,NY,NX)/VOLW(0,NY,NX))
      CNA1=AMAX1(ZEROC,ZNA(0,NY,NX)/VOLW(0,NY,NX))
      CKA1=AMAX1(ZEROC,ZKA(0,NY,NX)/VOLW(0,NY,NX))
      CSO41=AMAX1(ZEROC,ZSO4(0,NY,NX)/VOLW(0,NY,NX))
      CCL1=AMAX1(ZEROC,ZCL(0,NY,NX)/VOLW(0,NY,NX))
      CCO20=AMAX1(ZEROC,CCO2S(0,NY,NX)/12.0)
      CHCO31=AMAX1(ZEROC,ZHCO3(0,NY,NX)/VOLW(0,NY,NX))
      CCO31=AMAX1(ZEROC,ZCO3(0,NY,NX)/VOLW(0,NY,NX))
C
C     RESET PH IN LITTER
C
C     CHY2,COH2=H+,OH- concentrations (mol m-3)
C     AHY2,AOH2=H+,OH- activity (mol m-3)
C     A1=ion activity coefficient from ‘hour1.f’
C     DPH2O=H2O solubility product from PARAMETER (mol m-3)2
C     RHHX=flux to equilibrate H+, OH- at DPH2O
C     TRHY,TROH=total H+,OH- transformations used in ‘redist.f’
C     PH=litter pH 
C
      IF(CHY1.GT.ZEROC)THEN
      AHY2=CHY1*A1(0,NY,NX)
      ELSE
      AHY2=CHY1
      ENDIF
      IF(COH1.GT.ZEROC)THEN
      AOH2=COH1*A1(0,NY,NX)
      ELSE
      AOH2=COH1
      ENDIF
      S0=AHY2+AOH2
      S1=AMAX1(0.0,S0**2-4.0*(AHY2*AOH2-DPH2O))
      RHHX=0.5*(S0-SQRT(S1)) 
      AHY1=AMAX1(ZEROC,AHY2-RHHX)
      AOH1=AMAX1(ZEROC,AOH2-RHHX)
      CHY1=AHY1/A1(0,NY,NX)
      COH1=AOH1/A1(0,NY,NX)
      PH(0,NY,NX)=-LOG10(AHY1*1.0E-03)
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,1110)'RHHX0',I,J,NFZ,NX,NY
C    2,RHHX,AHY2,AOH2,AHY1,AOH1,CHY1,COH1,A1(0,NY,NX),AHY1*AOH1 
C    3,DPH2O,ZHY(0,NY,NX),ZOH(0,NY,NX),PH(0,NY,NX)
C     ENDIF 
C
C     END RESET PH IN LITTER 
C
C     INITIALIZE ION ACTIVITIES
C
      AAL1=CAL1*A3(0,NY,NX)
      AFE1=CFE1*A3(0,NY,NX)
      ACA1=CCA1*A2(0,NY,NX)
      AMG1=CMG1*A2(0,NY,NX)
      ANA1=CNA1*A1(0,NY,NX)
      AKA1=CKA1*A1(0,NY,NX)
      ASO41=CSO41*A2(0,NY,NX)
      ACO21=CCO20
      AHCO31=CHCO31*A1(0,NY,NX) 
      ACO31=CCO31*A2(0,NY,NX)
      AH1P1=CH1P1*A2(0,NY,NX)
      AH2P1=CH2P1*A1(0,NY,NX)
      AN41=CN41*A1(0,NY,NX)
      AN31=CN31
C
C     PRECIPITATION-DISSOLUTION REACTIONS FOR H2PO4 AND CO-PRECIPITATES
C     CALCULATED FROM ACTIVITIES OF REACTANTS
C     AND PRODUCTS THROUGH SOLUTIONS FOR THEIR EQUILIBRIUM
C     CONCENTRATIONS 
C
C     for all precipitation-dissolution reactions:
C
C     SP*=solubility product from PARAMETER above
C     A*1=ion activity (mol m-3)
C     A*Q=ion equilibrium activity at SP* (mol m-3)
C     A1,A2,A3=activity coefficients from ‘hour1.f’ 
C     R*X=precipitation(+ve) or dissolution (-ve) rate 
C        (mol m-3 t-1)
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONN=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     TPD,TPA=maximum precipitation rate (mol m-3 t-1)
C     salt solute code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+
C          :*CA*=Ca2+,*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-
C          :*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
C          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
C          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
C          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
C          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH
C          :*MGC*=MgCO3,*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-
C          :*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus solute code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C          :*F2P*=F1H2PO4-,*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C          :*M1P*=MgHPO4,*COO*=COO-
C          :*1=non-band,*B=band
C
C     ALUMINUM HYDROXIDE (GIBBSITE)
C
      XMIN=FION0*AAL1
      AAL1Q=SPALO/AOH1**3
      RPALOX=AMAX1(-AMAX1(0.0,PALOH1),-TPD,-XMIN
     2,(AMIN1(TPD,XMIN,AAL1-AAL1Q)))
C
C     IRON HYDROXIDE FE(OH)3
C
      XMIN=FION0*AFE1
      AFE1Q=SPFEO/AOH1**3
      RPFEOX=AMAX1(-AMAX1(0.0,PFEOH1),-TPD,-XMIN
     2,(AMIN1(TPD,XMIN,AFE1-AFE1Q)))
C
C     CALCITE (CACO3)
C
      XMIN=FION0*AMIN1(ACA1,ACO31)
      ACA1Q=SPCAC/ACO31
      RPCACX=AMAX1(-AMAX1(0.0,PCACO1),-TPD 
     2,(AMIN1(TPD,XMIN,ACA1-ACA1Q)))
C    3*1.0/(1.0+AOH1/OHKI)
C
C     GYPSUM (CASO4)
C
      XMIN=FION0*AMIN1(ACA1,ASO41)
      ACA1Q=SPCAS/ASO41
      RPCASO=AMAX1(-AMAX1(0.0,PCASO1),-TPD
     2,(AMIN1(TPD,XMIN,ACA1-ACA1Q)))
C
C     ALUMINUM PHOSPHATE ALPO4 (VARISCITE)
C
C     IF(AAL1.GT.ZEROC)THEN
      XMINN=FION0*AH2P1
      XMINP=FION0*AMIN1(AAL1,AH2P1)
      AAL1Q=SPALP/AH0P1
      RPALPX=AMAX1(-AMAX1(0.0,PALPO1),-TPD,-XMINN
     2,(AMIN1(TPD,XMINP,AAL1-AAL1Q)))
C     ELSE
C     RPALPX=0.0
C     ENDIF
C     WRITE(*,1110)'RPALPX',I,J,NFZ,NX,NY
C    2,RPALPX,PALPO1,TPD,AAL1,AAL1Q,AH2P1,SPX
C
C     IRON PHOSPHATE FEPO4 (STRENGITE)
C
C     IF(AFE1.GT.ZEROC)THEN
      XMINN=FION0*AH2P1
      XMINP=FION0*AMIN1(AFE1,AH2P1)
      AFE1Q=SPFEP/AH0P1
      RPFEPX=AMAX1(-AMAX1(0.0,PFEPO1),-TPD,-XMINN
     2,(AMIN1(TPD,XMINP,AFE1-AFE1Q)))
C     ELSE
C     RPFEPX=0.0
C     ENDIF
C     WRITE(*,1110)'RFELPX',I,J,NFZ,NX,NY
C    2,RPFEPX,PFEPO1,TPD,AFE1,AH2P1,SPX
C
C     DICALCIUM PHOSPHATE CAHPO4
C
C     IF(ACA1.GT.ZEROC)THEN
      XMINN=FION0*AH2P1
      XMINP=FION0*AMIN1(ACA1,AH2P1)
      ACA1Q=SPCAD/AH1P1 
      RPCADX=AMAX1(-AMAX1(0.0,PCAPD1),-TPD,-XMINN
     2,(AMIN1(TPD,XMINP,ACA1-ACA1Q)))
C     ELSE
C     RPCADX=0.0
C     ENDIF
C
C     HYDROXYAPATITE CA5(PO4)3OH
C
C     IF(ACA1.GT.ZEROC)THEN
      XMINN=FION0*AH2P1
      XMINP=FION0*AMIN1(ACA1,AH2P1)
      ACA1Q=(SPCAH/(AH0P1**3*AOH1))**0.2
      RPCAHX=AMAX1(-AMAX1(0.0,PCAPH1),-TPA,-XMINN
     2,(AMIN1(TPA,XMINP,ACA1-ACA1Q)))
C     ELSE
C     RPCAHX=0.0
C     ENDIF
C
C     MONOCALCIUM PHOSPHATE CA(H2PO4)2
C
C     IF(ACA1.GT.ZEROC)THEN
      XMINP=FION0*AMIN1(ACA1,AH2P1)
      ACA1Q=SPCAM/AH2P1**2
      RPCAMX=AMAX1(-AMAX1(0.0,PCAPM1),-TPD 
     2,(AMIN1(TPD,XMINP,ACA1-ACA1Q)))
C     ELSE
C     RPCAMX=0.0
C     ENDIF
C     IF(I.GT.315)THEN
C     WRITE(*,2227)'RPCADX0',I,J,NFZ,NX,NY 
C    2,RPCADX,AH2P1,SYCAD2,ACA1,AOH1,ACA1,ACO21,ACO31,PCAPD1
C    3,ORGC(0,NY,NX),SPCAC/ACO31,H2PO4(0,NY,NX)
C    4,PCAPD(0,NY,NX),VOLW(0,NY,NX),BKVLX,BKVLW
C     WRITE(*,2227)'RPCAHX0',I,J,NFZ,NX,NY 
C    2,RPCAHX,AH2P1,AH0P1,ACA1,ACA1Q,PCAPH1,FION0 
C    3,ACA1**5*AH0P1**3*AOH1,SPCAH,PCAPH(0,NY,NX)
2227  FORMAT(A8,5I4,20E12.4)
C     ENDIF
C
C     EQUILIBRIUM X-CA CONCENTRATION FROM CEC, GAPON COEFFICIENTS 
C     AND CATION CONCENTRATIONS
C
C     BKVLX=litter mass (Mg)
C     X*1=cation exchangeable concentration (mol Mg-1)
C     cation code
C        :HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ 
C
      XHY1=AMAX1(ZEROC,XHY(0,NY,NX)/BKVLX)
      XAL1=AMAX1(ZEROC,XAL(0,NY,NX)/BKVLX)
      XFE1=AMAX1(ZEROC,XFE(0,NY,NX)/BKVLX)
      XCA1=AMAX1(ZEROC,XCA(0,NY,NX)/BKVLX)
      XMG1=AMAX1(ZEROC,XMG(0,NY,NX)/BKVLX)
      XNA1=AMAX1(ZEROC,XNA(0,NY,NX)/BKVLX)
      XKA1=AMAX1(ZEROC,XKA(0,NY,NX)/BKVLX)
      XHC1=AMAX1(ZEROC,XHC(0,NY,NX)/BKVLX)
      AALX=AMAX1(ZEROC,AAL1)**0.333
      AFEX=AMAX1(ZEROC,AFE1)**0.333
      ACAX=AMAX1(ZEROC,ACA1)**0.500
      AMGX=AMAX1(ZEROC,AMG1)**0.500
C
C     EQUILIBRIUM X-CA CONCENTRATION FROM CEC AND CATION
C     CONCENTRATIONS
C
C     CCEC0=total cation exchange concentration (mol Mg-1)
C     XCAX=equilibrium R-Ca concentration (mol Mg-1)
C     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
C        CA-NH4,CA-H,CA-AL,CA-MG,CA-NA,CA-K
C     X*Q,X*Y=equilibrium,current exchangeable cation concentrations
C        (mol Mg-1)
C     XTLQ,XTL1=total equilibrium exchangeable concentration 
C        (mol, mol(+) Mg-1
C     cation code
C        :HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ 
C
      IF(ACAX.GT.ZERO.AND.CCEC0.GT.ZERO)THEN
      XCAX=CCEC0/(1.0
     2+GKC4(NU(NY,NX),NY,NX)*AN41/ACAX
     3+GKCH(NU(NY,NX),NY,NX)*AHY1/ACAX
     3+GKCA(NU(NY,NX),NY,NX)*AALX/ACAX*3.0
     4+GKCA(NU(NY,NX),NY,NX)*AFEX/ACAX*3.0
     4+GKCM(NU(NY,NX),NY,NX)*AMGX/ACAX*2.0
     5+GKCN(NU(NY,NX),NY,NX)*ANA1/ACAX
     5+GKCK(NU(NY,NX),NY,NX)*AKA1/ACAX)
      XN4Q=XCAX*AN41/ACAX*GKC4(NU(NY,NX),NY,NX)
      XHYQ=XCAX*AHY1/ACAX*GKCH(NU(NY,NX),NY,NX)
      XALQ=XCAX*AALX/ACAX*GKCA(NU(NY,NX),NY,NX)
      XFEQ=XCAX*AFEX/ACAX*GKCA(NU(NY,NX),NY,NX)
      XCAQ=XCAX*ACAX/ACAX 
      XMGQ=XCAX*AMGX/ACAX*GKCM(NU(NY,NX),NY,NX)
      XNAQ=XCAX*ANA1/ACAX*GKCN(NU(NY,NX),NY,NX)
      XKAQ=XCAX*AKA1/ACAX*GKCK(NU(NY,NX),NY,NX)
      XTLQ=XN4Q+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ 
      XTL1=XN41+XHY1+XAL1*3.0+XFE1*3.0+XCA1*2.0+XMG1*2.0+XNA1+XKA1 
      IF(XTLQ.GT.ZERO)THEN
      FX=CCEC0/XTLQ
      ELSE
      FX=0.0
      ENDIF
      IF(XTL1.GT.ZERO)THEN
      FY=CCEC0/XTL1
      ELSE
      FY=0.0
      ENDIF
      XN4Q=FX*XN4Q
      XHYQ=FX*XHYQ
      XALQ=FX*XALQ
      XFEQ=FX*XFEQ
      XCAQ=FX*XCAQ
      XMGQ=FX*XMGQ
      XNAQ=FX*XNAQ
      XKAQ=FX*XKAQ
      XN4Y=FY*XN41
      XHYY=FY*XHY1
      XALY=FY*XAL1*3.0
      XFEY=FY*XFE1*3.0
      XCAY=FY*XCA1*2.0
      XMGY=FY*XMG1*2.0
      XNAY=FY*XNA1
      XKAY=FY*XKA1
C
C     NH4,H,AL,FE,CA,MG,NA,K EXCHANGE IN SURFACE LITTER
C
C     FIONC=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     BKVLW=soil mass:water (Mg m-3)
C     TADCC=maximum cation adsorption rate (mol Mg t-1)
C     XMIN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     RX*Q,RX*=equilibrium,current cation desorption(-ve),
C        adsorption(+ve) rate (mol m-3 t-1)
C     X*Q,X*Y=equilibrium,current exchangeable cation concentrations
C        (mol Mg-1)
C     C*1,A*1=aqueous cation concentration,activity (mol m-3) 
C     cation code
C        :HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ non-band,NB=NH4+ band
C
      FIONC=FION0/BKVLW
      TADCC=TADC/BKVLW
      XMIN=FIONC*AMIN1(XN41,CN41)
      RXN4Q=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XN4Q-XN4Y)))
      XMIN=FIONC*AMIN1(XHY1,CHY1)
      RXHYQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XHYQ-XHYY)))
      XMIN=FIONC*3.0*AMIN1(XAL1,CAL1)
      RXALQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XALQ-XALY)))
      XMIN=FIONC*3.0*AMIN1(XFE1,CFE1)
      RXFEQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XFEQ-XFEY)))
      XMIN=FIONC*2.0*AMIN1(XCA1,CCA1)
      RXCAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XCAQ-XCAY)))
      XMIN=FIONC*2.0*AMIN1(XMG1,CMG1)
      RXMGQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XMGQ-XMGY)))
      XMIN=FIONC*AMIN1(XNA1,CNA1)
      RXNAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XNAQ-XNAY)))
      XMIN=FIONC*AMIN1(XKA1,CKA1)
      RXKAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XKAQ-XKAY)))
      TXXX=RXN4Q+RXNBQ+RXHYQ+RXALQ+RXFEQ+RXCAQ+RXMGQ+RXNAQ+RXKAQ
      TXXY=ABS(RXN4Q)+ABS(RXNBQ)+ABS(RXHYQ)+ABS(RXALQ)+ABS(RXFEQ)
     3+ABS(RXCAQ)+ABS(RXMGQ)+ABS(RXNAQ)+ABS(RXKAQ)
      RXN4=RXN4Q-TXXX*ABS(RXN4Q)/TXXY
      RXNB=RXNBQ-TXXX*ABS(RXNBQ)/TXXY
      RXHY=RXHYQ-TXXX*ABS(RXHYQ)/TXXY
      RXAL=(RXALQ-TXXX*ABS(RXALQ)/TXXY)/3.0
      RXFE=(RXFEQ-TXXX*ABS(RXFEQ)/TXXY)/3.0
      RXCA=(RXCAQ-TXXX*ABS(RXCAQ)/TXXY)/2.0
      RXMG=(RXMGQ-TXXX*ABS(RXMGQ)/TXXY)/2.0
      RXNA=RXNAQ-TXXX*ABS(RXNAQ)/TXXY
      RXKA=RXKAQ-TXXX*ABS(RXKAQ)/TXXY
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,2223)'RXN40S',I,J,NFZ
C    2,XN4Q,XHYQ,XALQ,XFEQ,XCAQ,XMGQ,XNAQ,XKAQ
C    2,XN41,XHY1,XAL1,XFE1,XCA1,XMG1,XNA1,XKA1
C    2,XN4Y,XHYY,XALY,XFEY,XCAY,XMGY,XNAY,XKAY
C    2,AN41,AHY1,AAL1,AFE1,ACA1,AMG1,ANA1,AKA1
C    2,RXN4,RXHY,RXAL,RXFE,RXCA,RXMG,RXNA,RXKA
C    2,RXN4+RXHY+RXAL*3.0+RXFE*3.0+RXCA*2.0+RXMG*2.0+RXNA+RXKA
C    2,XN4Y+XHYY+XALY*3.0+XFEY*3.0+XCAY*2.0+XMGY*2.0+XNAY+XKAY
C    2,CCEC0,XCAX,XTLQ,FX,FY,A1(0,NY,NX),BKVLX,XN4(0,NY,NX),TXXX,TXXY
C     ENDIF
      ELSE
      RXN4=0.0
      RXHY=0.0
      RXAL=0.0
      RXFE=0.0
      RXCA=0.0
      RXMG=0.0
      RXNA=0.0
      RXKA=0.0
      ENDIF
C
C     DISSOCIATION OF CARBOXYL RADICALS
C
C     XMIN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     TADCC=maximum cation adsorption rate (mol Mg t-1)
C     CCEC0=litter CEC (mol Mg-1)
C     XCOOH,XHC1=total carboxyl exchange sites,occupied by H+ (mol)
C     RXHC=COOH-COO+H desorption(-ve) or adsorption(+ve)(mol C m-3 t-1)
C
      XMIN=FION0/BKVLW*XHC1
      TADCC=TADC/BKVLW
      XCOOH=CCEC0
      XCOO=AMAX1(ZEROC,XCOOH-XHC1)
      XCOOQ=AMIN1(XCOOH,DPCOH*XHC1/AHY1)
      RXHC=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XCOO-XCOOQ))) 
C     WRITE(*,1110)'RXHC0',I,J,NFZ,NX,NY
C    2,RXHC,RXHC*BKVLW,XHC1,XCOO,XCOOQ,XCOOH,CHY1,AHY1,DPCOH,COOH 
C    2,ZHY(0,NY,NX),BKVL(0,NY,NX),BKVLW 
C
C     ION PAIRING REACTIONS IN LITTER
C
C     for all reactions:
C
C     DP*=dissociation constant from PARAMETER
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONH=limits on reaction rates imposed by substrate, product
C        activities (h-1)
C     A*1=ion activity (mol m-3)
C     A*Q=ion equilibrium activity at DP* (mol m-3)
C     TSL,TSZ=maximum reaction rate (mol m-3 t-1)
C     R*=dissociation(-ve) or association(+ve) rate (mol m-3 t-1)
C     solute code
C        :*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+
C        :*CA*=Ca2+,*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-
C        :*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
C
C     NH4-NH3+H IN SURFACE LITTER
C
      XMINN=FIONH*AN41
      XMINP=FIONH*AMIN1(AHY1,AN31)
      AN3Q=DPN4*AN41/AHY1
      RNH4=AMAX1(-TSZ,-XMINN
     2,(AMIN1(TSZ,XMINP,AN31-AN3Q)))
      CHYX=CHY1-RNH4
      AHY1=CHYX*A1(0,NY,NX)
C     WRITE(*,1110)'RNH40S',I,J,NFZ,NX,NY
C    2,RNH4,TSL,AHY1,AN31,DPN4,AN41,CHY1,A1(0,NY,NX)
C
C     RHCO3=HCO3-CO3+H dissociation
C
      XMINN=FION0*AHCO31
      XMINP=FION0*AMIN1(AHY1,ACO31)
      ACO3Q=DPHCO*AHCO31/AHY1
      RHCO3=AMAX1(-TSL,-XMIN
     2,(AMIN1(TSL,XMINP,ACO31-ACO3Q)))
C
C     RCO2Q=CO2-HCO3+H dissociation
C
      XMINN=FION0*ACO21
      XMINP=FION0*AMIN1(AHY1,AHCO31)
      AHCO3Q=DPCO2*ACO21/AHY1
      RCO2Q=AMAX1(-TSL,-XMINN
     2,(AMIN1(TSL,XMINP,AHCO31-AHCO3Q)))
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,1110)'RCO2Q0',I,J,NFZ,NX,NY
C    2,RCO2Q,TSZ,AHY1,ACO21,AHCO31,AHCO3Q,AHY1*AHCO31/ACO21
C    3,DPCO2,CHY1,CHYX,CHCO31,CHCO31X,RHCO3,FION0,XMINN,XMINP 
C     ENDIF
C
C     PHOSPHORUS ANION PAIRING IN SURFACE LITTER
C     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
C     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
C     EXCHANGE SITES (NOT CALCULATED)
C
C     H2PO4-H+HPO4
C    
      XMINN=FION0*AH2P1
      XMINP=FION0*AMIN1(AHY1,AH1P1)
      AH1P1Q=DPH2P*AH2P1/AHY1
      RH2P=AMAX1(-TSL,-XMINN
     2,(AMIN1(TSL,XMINP,AH1P1-AH1P1Q)))
C     WRITE(*,1113)'RH2P0',I,J,NFZ,NX,NY
C    2,RH2P,TSL,AH1P1,AHY1,DPH2P,AH2P1,FION0,PH(0,NY,NX)
      ELSE
C
C     IF NO SALT IS SELECTED IN SITE FILE THEN A SUBSET
C     OF LITTER EQUILIBRIA REACTIONS ARE SOLVED: MOSTLY THOSE
C     FOR PHOSPHORUS AND CO-REACTANTS
C
C     INITIALIZE ION CONCENTRATIONS, ACTIVITIES
C
      CHY1=AMAX1(ZEROC,ZHY(0,NY,NX)/VOLW(0,NY,NX))
      COH1=AMAX1(ZEROC,ZOH(0,NY,NX)/VOLW(0,NY,NX))
      CAL1=AMAX1(ZEROC,ZAL(0,NY,NX)/VOLW(0,NY,NX))
      CFE1=AMAX1(ZEROC,ZFE(0,NY,NX)/VOLW(0,NY,NX))
      CCA1=AMAX1(ZEROC,ZCA(0,NY,NX)/VOLW(0,NY,NX))
      CMG1=AMAX1(ZEROC,ZMG(0,NY,NX)/VOLW(0,NY,NX))
      CNA1=AMAX1(ZEROC,ZNA(0,NY,NX)/VOLW(0,NY,NX))
      CKA1=AMAX1(ZEROC,ZKA(0,NY,NX)/VOLW(0,NY,NX))
      CSO41=AMAX1(ZEROC,ZSO4(0,NY,NX)/VOLW(0,NY,NX))
      CCL1=AMAX1(ZEROC,ZCL(0,NY,NX)/VOLW(0,NY,NX))
      CCO20=AMAX1(ZEROC,CCO2S(0,NY,NX)/12.0)
      CHCO31=AMAX1(ZEROC,CCO20*DPHCO/CHY1)
      CCO31=AMAX1(ZEROC,CCO20*DPCO3/CHY1**2)
      AHY1=CHY1*A1(0,NY,NX)
      AOH1=COH1*A1(0,NY,NX)
      AAL1=CAL1*A3(0,NY,NX)
      AFE1=CFE1*A3(0,NY,NX)
      ACA1=CCA1*A2(0,NY,NX)
      AMG1=CMG1*A2(0,NY,NX)
      ANA1=CNA1*A1(0,NY,NX)
      AKA1=CKA1*A1(0,NY,NX)
      ASO41=CSO41*A2(0,NY,NX)
      ACO21=CCO20
      AHCO31=CHCO31*A1(0,NY,NX)
      ACO31=CCO31*A2(0,NY,NX)
      AH1P1=CH1P1*A2(0,NY,NX)
      AH2P1=CH2P1*A1(0,NY,NX)
      AN41=CN41*A1(0,NY,NX)
      AN31=CN31
C
C     CALCULATE H2PO4 COPRECIPITATES FROM LITTER PH
C
C     VOLW=litter water content (m3)
C     C*1,A*1=litter ion concentration,activity (mol m-3)
C     Z*=litter ion content (mol)
C     ion code
C       :HY=H+,OH=OH-,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+
C       :NA=Na+,KA=K+,SO4=SO42-,CL=Cl-,CO2=CO2,CO3=CO3
C       :H1P=HPO42-,H2P=H2PO4- 
C     DP*,SP*=dissociation,solubility products from PARAMETER above
C     R*X=precipitation(+ve) or dissolution (-ve) rate 
C        (mol m-3 t-1)
C     TPD=maximum precipitation rate (mol m-3 t-1)
C
C     ALUMINUM PHOSPHATE ALPO4 (VARISCITE)
C
      XMIN=FION0*AH2P1
      AAL1Q=SPALO/AOH1**3
      AH2PQ=SHA0P2*AHY1**2/AAL1Q
      RPALPX=AMAX1(-AMAX1(0.0,PALPO1),-TPD,-XMIN
     2,(AMIN1(TPD,XMIN,AH2P1-AH2PQ)))
C
C     IRON PHOSPHATE FEPO4 (STRENGITE)
C
      XMIN=FION0*AH2P1
      AFE1Q=SPFEO/AOH1**3
      AH2PQ=SHF0P2*AHY1**2/AFE1Q
      RPFEPX=AMAX1(-AMAX1(0.0,PFEPO1),-TPD,-XMIN
     2,(AMIN1(TPD,XMIN,AH2P1-AH2PQ)))
C
C     DICALCIUM PHOSPHATE CAHPO4
C
      XMIN=FION0*AMIN1(ACA1,AH1P1)
      AH1PQ=SPCAD/ACA1 
C     RPCADX=AMAX1(-AMAX1(0.0,PCAPD1),-TPD,(AMIN1(TPD
C    2,FION0*AH1P1,ACA1-ACA1Q)))
      RPCADX=AMAX1(-AMAX1(0.0,PCAPD1),-TPD,-XMIN
     2,(AMIN1(TPD,XMIN,AH1P1-AH1PQ)))
C
C     HYDROXYAPATITE CA5(PO4)3OH
C
      XMIN=FION0*AMIN1(ACA1,AH2P1)
      AH2PH=(SHCAH2*AHY1**7/ACA1**5)**0.333
      RPCAHX=AMAX1(-AMAX1(0.0,PCAPH1),-TPA,-XMIN
     2,(AMIN1(TPA,XMIN,AH2P1-AH2PH)))
C
C     MONOCALCIUM PHOSPHATE CA(H2PO4)2
C
      XMIN=FION0*AMIN1(ACA1,AH2P1)
      AH2PQ=(SPCAM/ACA1)**0.5
      RPCAMX=AMAX1(-AMAX1(0.0,PCAPM1),-TPD 
     2,(AMIN1(TPD,XMIN,AH2P1-AH2PQ)))
C     IF(I.GT.315)THEN
C     WRITE(*,2228)'RPCADX0',I,J,NFZ,RPCADX,AH2P1,AH2PH 
C    2,SYCAD2,ACA1,AOH1,CCA1,ACO21,ACO31,PCAPD1
C    3,ORGC(0,NY,NX),SPCAC/ACO31,H2PO4(0,NY,NX)
C    4,PCAPD(0,NY,NX),VOLW(0,NY,NX),BKVLX,BKVLW
C     WRITE(*,2228)'RPCAHX0',I,J,NFZ,RPCAHX,PCAPH1,BKVLW 
C    2,AH2P1,AH2PH,SHCAH2,AHY1,ACA1 
C    3,VOLW(0,NY,NX),SPCAC/ACO31,H2PO4(0,NY,NX)
C    4,CCO20,DPCO3,AHY1,CCO2S(0,NY,NX)
2228  FORMAT(A8,3I4,20E12.4)
C     ENDIF
C
C     EQUILIBRIUM X-CA CONCENTRATION FROM CEC, GAPON COEFFICIENTS 
C     AND CATION CONCENTRATIONS
C
C     CCEC,XCEC=cation exchange concentration,capacity (mol Mg-1,mol)
C     BKVLX=litter mass (Mg)
C     CCEC,XCEC=cation exchange concentration,capacity (mol Mg-1,mol)
C     A*1=cation solute activity (mol m-3)
C     cation code
C        :HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ 
C
      AALX=AMAX1(ZEROC,AAL1)**0.333
      AFEX=AMAX1(ZEROC,AFE1)**0.333
      ACAX=AMAX1(ZEROC,ACA1)**0.500
      AMGX=AMAX1(ZEROC,AMG1)**0.500
C
C     EQUILIBRIUM X-CA CONCENTRATION FROM CEC AND CATION
C     CONCENTRATIONS
C
C     CCEC,XCEC=total cation exchange concentration,capacity 
C        (mol Mg-1,mol)
C     XCAX=equilibrium X-Ca concentration (mol+ Mg-1)
C     A*1=cation activity (mol m-3)
C     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
C        CA-NH4,CA-H,CA-AL,CA-MG,CA-NA,CA-K
C     X*Q,X*Y=equilibrium,current exchangeable cation concentrations
C        (mol Mg-1)
C     XTLQ,XTL1=total equilibrium exchangeable concentration 
C        (mol, mol(+) Mg-1
C     cation code
C        :HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ 
C
      IF(ACAX.GT.ZERO.AND.CCEC0.GT.ZERO)THEN
      XCAX=CCEC0/(1.0
     2+GKC4(NU(NY,NX),NY,NX)*AN41/ACAX
     3+GKCH(NU(NY,NX),NY,NX)*AHY1/ACAX
     3+GKCA(NU(NY,NX),NY,NX)*AALX/ACAX*3.0
     4+GKCA(NU(NY,NX),NY,NX)*AFEX/ACAX*3.0
     4+GKCM(NU(NY,NX),NY,NX)*AMGX/ACAX*2.0
     5+GKCN(NU(NY,NX),NY,NX)*ANA1/ACAX
     5+GKCK(NU(NY,NX),NY,NX)*AKA1/ACAX)
      XN4Q=XCAX*AN41/ACAX*GKC4(NU(NY,NX),NY,NX)
      XHYQ=XCAX*AHY1/ACAX*GKCH(NU(NY,NX),NY,NX)
      XALQ=XCAX*AALX/ACAX*GKCA(NU(NY,NX),NY,NX)
      XFEQ=XCAX*AFEX/ACAX*GKCA(NU(NY,NX),NY,NX)
      XCAQ=XCAX*ACAX/ACAX 
      XMGQ=XCAX*AMGX/ACAX*GKCM(NU(NY,NX),NY,NX)
      XNAQ=XCAX*ANA1/ACAX*GKCN(NU(NY,NX),NY,NX)
      XKAQ=XCAX*AKA1/ACAX*GKCK(NU(NY,NX),NY,NX)
      XTLQ=XN4Q+XHYQ+XALQ+XFEQ+XCAQ+XMGQ+XNAQ+XKAQ 
      XTL1=XN41+XHY1+XAL1*3.0+XFE1*3.0+XCA1*2.0+XMG1*2.0+XNA1+XKA1 
      IF(XTLQ.GT.ZERO)THEN
      FX=CCEC0/XTLQ
      ELSE
      FX=0.0
      ENDIF
      IF(XTL1.GT.ZERO)THEN
      FY=CCEC0/XTL1
      ELSE
      FY=0.0
      ENDIF
      XN4Q=FX*XN4Q
      XHYQ=FX*XHYQ
      XALQ=FX*XALQ
      XFEQ=FX*XFEQ
      XCAQ=FX*XCAQ
      XMGQ=FX*XMGQ
      XNAQ=FX*XNAQ
      XKAQ=FX*XKAQ
      XN4Y=FY*XN41
      XHYY=FY*XHY1
      XALY=FY*XAL1*3.0
      XFEY=FY*XFE1*3.0
      XCAY=FY*XCA1*2.0
      XMGY=FY*XMG1*2.0
      XNAY=FY*XNA1
      XKAY=FY*XKA1
C
C     NH4,H,AL,FE,CA,MG,NA,K EXCHANGE IN SURFACE LITTER
C
C     FIONC=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     BKVLW=soil mass:water (Mg m-3)
C     TADCC=maximum cation adsorption rate (mol Mg t-1)
C     XMIN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     RX*Q,RX*=equilibrium,current cation desorption(-ve),
C        adsorption(+ve) rate (mol m-3 t-1)
C     X*Q,X*Y=equilibrium,current exchangeable cation concentrations
C        (mol Mg-1)
C     C*1,A*1=aqueous cation concentration,activity (mol m-3) 
C     cation code
C        :HY=H+,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+,KA=K+
C        :N4=NH4+ non-band 
C
      FIONC=FION0/BKVLW
      TADCC=TADC/BKVLW
      XMIN=FIONC*AMIN1(XN41,CN41)
      RXN4Q=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XN4Q-XN4Y)))
      XMIN=FIONC*AMIN1(XHY1,CHY1)
      RXHYQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XHYQ-XHYY)))
      XMIN=FIONC*3.0*AMIN1(XAL1,CAL1)
      RXALQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XALQ-XALY)))
      XMIN=FIONC*3.0*AMIN1(XFE1,CFE1)
      RXFEQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XFEQ-XFEY)))
      XMIN=FIONC*2.0*AMIN1(XCA1,CCA1)
      RXCAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XCAQ-XCAY)))
      XMIN=FIONC*2.0*AMIN1(XMG1,CMG1)
      RXMGQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XMGQ-XMGY)))
      XMIN=FIONC*AMIN1(XNA1,CNA1)
      RXNAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XNAQ-XNAY)))
      XMIN=FIONC*AMIN1(XKA1,CKA1)
      RXKAQ=AMAX1(-TADCC,-XMIN
     2,(AMIN1(TADCC,XMIN,XKAQ-XKAY)))
      TXXX=RXN4Q+RXNBQ+RXHYQ+RXALQ+RXFEQ+RXCAQ+RXMGQ+RXNAQ+RXKAQ
      TXXY=ABS(RXN4Q)+ABS(RXNBQ)+ABS(RXHYQ)+ABS(RXALQ)+ABS(RXFEQ)
     3+ABS(RXCAQ)+ABS(RXMGQ)+ABS(RXNAQ)+ABS(RXKAQ)
      RXN4=RXN4Q-TXXX*ABS(RXN4Q)/TXXY
      RXNB=RXNBQ-TXXX*ABS(RXNBQ)/TXXY
      RXHY=RXHYQ-TXXX*ABS(RXHYQ)/TXXY
      RXAL=(RXALQ-TXXX*ABS(RXALQ)/TXXY)/3.0
      RXFE=(RXFEQ-TXXX*ABS(RXFEQ)/TXXY)/3.0
      RXCA=(RXCAQ-TXXX*ABS(RXCAQ)/TXXY)/2.0
      RXMG=(RXMGQ-TXXX*ABS(RXMGQ)/TXXY)/2.0
      RXNA=RXNAQ-TXXX*ABS(RXNAQ)/TXXY
      RXKA=RXKAQ-TXXX*ABS(RXKAQ)/TXXY
C     IF(J.EQ.12)THEN
C     WRITE(*,2223)'RXN40X',I,J,NFZ
C    2,XN4Q,XHYQ,XALQ,XFEQ,XCAQ,XMGQ,XNAQ,XKAQ
C    2,XN41,XHY1,XAL1,XFE1,XCA1,XMG1,XNA1,XKA1
C    2,XN4Y,XHYY,XALY,XFEY,XCAY,XMGY,XNAY,XKAY
C    2,AN41,AHY1,AAL1,AFE1,ACA1,AMG1,ANA1,AKA1
C    2,RXN4,RXHY,RXAL,RXFE,RXCA,RXMG,RXNA,RXKA
C    2,RXN4+RXHY+RXAL*3.0+RXFE*3.0+RXCA*2.0+RXMG*2.0+RXNA+RXKA
C    2,XN4Y+XHYY+XALY*3.0+XFEY*3.0+XCAY*2.0+XMGY*2.0+XNAY+XKAY
C    2,CCEC0,XCAX,XTLQ,FX,FY,A1(0,NY,NX),AALX,AFEX,TXXX,TXXY
2223  FORMAT(A8,3I4,/,8F16.8,/,8F16.8,/,8F16.8,/,8F16.8,/
     2,8F16.8,/,8F16.8,/,8F16.8)
C     ENDIF
      ELSE
      RXN4=0.0
      RXHY=0.0
      RXAL=0.0
      RXFE=0.0
      RXCA=0.0
      RXMG=0.0
      RXNA=0.0
      RXKA=0.0
      ENDIF
C
C     NH4-NH3+H IN SURFACE LITTER
C
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FIONH=limits on reaction rates imposed by substrate, product
C        activities (h-1)
C     RNH4=NH4 dissociation(-ve) or association(+ve) rate (mol m-3 t-1) 
C     DPN4=NH4 dissociation constant from PARAMETER
C     AN41,AN31=total NH4,NH3 activity (mol m-3) 
C     AN3Q=NH3 equilibrium activity at DPNH4 (mol m-3)
C     AHY1=H activity (mol m-3)
C     TSL=maximum reaction rate (mol m-3 t-1)
C
      XMINN=FIONH*AN41
      XMINP=FIONH*AMIN1(AHY1,AN31)
      AN3Q=DPN4*AN41/AHY1
      RNH4=AMAX1(-TSZ,-XMINN
     2,(AMIN1(TSZ,XMINP,AN31-AN3Q)))
C
C     H2PO4-H+HPO4
C
C     XMINP,XMINN=limits on reaction rates imposed by substrate,
C        product activities (mol m-3 t-1)
C     FION0=limits on reaction rates imposed by substrate, product
C        activities (t-1)
C     DPH2P=H2PO4 dissociation constant from PARAMETER
C     RH2P=H2PO4 dissociation(+ve), association(+ve) (mol m-3 t-1)
C     AHY1,AH2P1=H,H2PO4 activity (mol m-3)
C     AH1P1Q=HPO4 equilibrium activity at DPH2P (mol P m-3)
C     TSL=maximum reaction rate (mol m-3 t-1)
C     
      AH1P1Q=DPH2P*AH2P1/AHY1
      RH2P=AMAX1(-TSL,-FION0*AH2P1
     2,(AMIN1(TSL,FION0*AH1P1,AH1P1-AH1P1Q)))
C     WRITE(*,9989)'RH2P',I,J,NFZ,RH2P,TSL,AH1P1,AHY1
C    2,DPH2P,AH2P1,DPH2P,AHY1
      ENDIF
      ELSE
      BKVLX=0.0
      RHHX=0.0
      RZHYS=0.0
      AHY1=10.0**(-PH(0,NY,NX))*1.0E+03
      AOH1=DPH2O/AHY1
      RSN4AA=0.0
      RSN3AA=0.0
      RSNUAA=0.0
      RSNOAA=0.0
      RPALOX=0.0
      RPFEOX=0.0
      RPCACX=0.0
      RPCASO=0.0
      RPALPX=0.0
      RPFEPX=0.0
      RPCADX=0.0
      RPCAHX=0.0
      RPCAMX=0.0
      RN4X=0.0
      RN3X=0.0
      CN41=0.0
      CN31=0.0
      XN41=0.0
      RH1PX=0.0
      RH2PX=0.0
      CH1P1=0.0
      CH2P1=0.0
      RXN4=0.0
      RNH4=0.0
      RHCO3=0.0
      RCO2Q=0.0
      RH2P=0.0
      RXHY=0.0
      RXAL=0.0
      RXFE=0.0
      RXCA=0.0
      RXMG=0.0
      RXNA=0.0
      RXKA=0.0
      RXHC=0.0
      PALOH1=0.0
      PFEOH1=0.0
      PCACO1=0.0
      PCASO1=0.0
      BKVLW=1.0
      ENDIF
C
C     TOTAL ION TRANSFORMATIONS FOR ALL REACTIONS ABOVE
C
C     RN4S=net NH4 transformation in litter (mol m-3 t-1)
C     RN3S=net NH3 transformation in litter (mol m-3 t-1)
C     RHP1,RHP2=net HPO4,H2PO4 transformation in litter (mol m-3 t-1)
C     ROH,RAL,RFE,RCA,RMG,RNA,RKA=net OH,Al,Fe,Ca,Mg,Na,K
C        transformation in litter (mol m-3 t-1) 
C     BKVLW=litter mass:water (Mg m-3)
C
      RN4S=RNH4-RXN4*BKVLW
      RN3S=-RNH4
      RHP1=-RH2P 
      RHP2=RH2P-RPALPX-RPFEPX-RPCADX-2.0*RPCAMX-3.0*RPCAHX 
      RHY=-RXHY*BKVLW 
      ROH=0.0 
      RAL=-RXAL*BKVLW 
      RFE=-RXFE*BKVLW 
      RCA=-RXCA*BKVLW 
      RMG=-RXMG*BKVLW 
      RNA=-RXNA*BKVLW 
      RKA=-RXKA*BKVLW 
C
C     TOTAL SALT ION TRANSFORMATIONS FOR ALL REACTIONS ABOVE
C     IF SALT OPTION IS SELECTED
C
C     ISALTG:0=salt concentrations entered in soil file generate
C              equilibrium concentrations that remain static during
C              model run
C           :1=salt equilibrium concentrations are solved
C              dynamically in ‘solute.f’ and transported in ‘trnsfrs.f’ 
C     RHY,ROH,RAL,RFE,RCA,RHCO,RSO4=net H,OH,Al,Fe,Ca,HCO3,SO4
C        transformation including salt reactions (mol m-3 t-1) 
C
      IF(ISALTG.NE.0)THEN
      RHY=RHY+RZHYS-RXHC*BKVLW-RCO2Q-RHCO3
     2+2.0*(RPALPX+RPFEPX)+RPCADX 
     3+6.0*RPCAHX-RH2P-RNH4 
      ROH=ROH-3.0*(RPALOX+RPFEOX)-RPCAHX
      RAL=RAL-RPALOX-RPALPX 
      RFE=RFE-RPFEOX-RPFEPX 
      RCA=RCA-RPCACX-RPCASO-RPCADX-5.0*RPCAHX-RPCAMX
      RCO3=-RHCO3-RPCACX 
      RHCO=-RCO2Q+RHCO3
      RCO2=RCO2Q
      RSO4=-RPCASO
      RH2O=RCO2Q 
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,1110)'RHY0',I,J,NFZ,NX,NY 
C    2,CHY1,AHY1,RHY,RZHYS,-RXHY*BKVLW,-RXHC*BKVLW,-RCO2Q,-RHCO3 
C    3,RPCACX,2.0*(RPALPX+RPFEPX),RPCADX,6.0*RPCAHX,-RH2P,-RNH4
C    4,ZHY(0,NY,NX)/VOLW(0,NY,NX),ZHY(0,NY,NX),VOLW(0,NY,NX)
C    5,RHHX,THETW0,XZHYS(0,NY,NX),PH(0,NY,NX) 
C     WRITE(*,1110)'ROH0',I,J,NFZ,NX,NY 
C    2,COH1,AOH1,ROH,-3.0*RPALOX,-3.0*RPFEOX,-RPCAHX 
C     ENDIF
      ELSE
C
C     IF NO SALT IS SELECTED THEN SET H+,OH- FROM PH IN SOIL FILE
C
C     PH=soil pH from soil file
C     DPH2O=H2O solubility product from PARAMETER (mol m-3)2
C     AHY2,AOH2=H,OH activity (mol m-3)
C     A1=ion activity coefficient from ‘hour1.f’
C     RHHY,RHOH=H,OH equilibration fluxes at DPH2O (mol m-3 t-1)
C
      AHY2=10.0**(-PH(0,NY,NX))*1.0E+03*A1(0,NY,NX)
      AOH2=DPH2O/AHY2
      RHHY=AHY2-AHY1-RHY
      RHOH=AOH2-AOH1-ROH
C     IF((I/30)*30.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,1110)'RHHY0',I,J,NFZ,NX,NY,AHY2,AOH2,AHY2*AOH2,SPX
C    2,RHHX,RHY,ROH
C    2,CHY1,COH1,AHY2,AOH2,A1(0,NY,NX) 
C    3,DPH2O,PH(0,NY,NX)
1110  FORMAT(A8,5I4,80E12.4)
C     ENDIF
      ENDIF
C
C     END SET H+,OH- 
C
C     CONVERT TOTAL SALT ION TRANSFORMATIONS FROM CHANGES IN 
C     CONCENTRATION TO CHANGES IN MASS (mol m-3 t-1 to mol t-1)
C     FOR USE IN 'REDIST'
C
C     VOLW=litter water content (m3)
C     TRSO4,TRCO3,TRHCO=total SO4,CO3,HCO3 transformation
C     TRXHC=total HCO3 adsorption 
C     TRALOH,TRFEOH,TRCACO,TRCASO=total AlOH3,FeOH3,CaCO3,CaSO4
C        precipitation 
C
      IF(ISALTG.NE.0)THEN
      TRCO3(0,NY,NX)=TRCO3(0,NY,NX)+RCO3*VOLW(0,NY,NX)
      TRCO2(0,NY,NX)=TRCO2(0,NY,NX)+RCO2*VOLW(0,NY,NX)
      TRHCO(0,NY,NX)=TRHCO(0,NY,NX)+RHCO*VOLW(0,NY,NX)
      TRSO4(0,NY,NX)=TRSO4(0,NY,NX)+RSO4*VOLW(0,NY,NX)
      TRXHC(0,NY,NX)=TRXHC(0,NY,NX)+RXHC*BKVLX
      TRALOH(0,NY,NX)=TRALOH(0,NY,NX)+RPALOX*VOLW(0,NY,NX)
      TRFEOH(0,NY,NX)=TRFEOH(0,NY,NX)+RPFEOX*VOLW(0,NY,NX)
      TRCACO(0,NY,NX)=TRCACO(0,NY,NX)+RPCACX*VOLW(0,NY,NX)
      TRCASO(0,NY,NX)=TRCASO(0,NY,NX)+RPCASO*VOLW(0,NY,NX)
      TRHY(0,NY,NX)=TRHY(0,NY,NX)+(RHY-RHHX)*VOLW(0,NY,NX)
      TROH(0,NY,NX)=TROH(0,NY,NX)+(ROH-RHHX)*VOLW(0,NY,NX)
      TRH2O(0,NY,NX)=TRH2O(0,NY,NX)+RH2O*VOLW(0,NY,NX)
      TBH2O(0,NY,NX)=TBH2O(0,NY,NX)+RHHX*VOLW(0,NY,NX)
C     WRITE(*,1113)'RAL0',I,J,NFZ,NX,NY
C    2,RFE,RPALOX,RPALPX,RXAL,CAL1,CH2P1 
C    2,PALPO1,BKVLW,TPD,CAL1*CH2P1,SHF0P2*CHY1**2
C    3,TRAL(0,NY,NX),TRMG(0,NY,NX),VOLW(0,NY,NX)
C     WRITE(*,1113)'RFE0',I,J,NFZ,NX,NY
C    2,RFE,RPFEOX,RPFEPX,RXFE,CFE1,CH2P1 
C    2,PFEPO1,BKVLW,TPD,CFE1*CH2P1,SHF0P2*CHY1**2
C    3,XFEQ,XFE1,CFEX,XFEQ
C    3,TRFE(0,NY,NX),TRMG(0,NY,NX),VOLW(0,NY,NX)
C     WRITE(*,1113)'RCA0',I,J,NFZ,NX,NY
C    2,RCA,CCA1,RPCACX,RPCASO
C    2,RPCADX,5.0*RPCAHX,RPCAMX
C    2,PCAPH1,BKVLW,TPD,CH2P1,CH2PH,SHCAH2,CHY1
C    2,TRCA(0,NY,NX),VOLW(0,NY,NX)
C     WRITE(*,1113)'TRHY0',I,J,NFZ,NX,NY
C    2,TRHY(0,NY,NX),RHY,ROH,CHY1,COH1,AHY1,AOH1
C    2,RPCACX,2.0*(RPALPX+RPFEPX),3.0*RPCAHX,RH2P,RNH4,RXHY,RXHC
C    3,RHHY,RHOH,AHY2,AOH2,TRHY(0,NY,NX),ZHY(0,NY,NX)
C    4,ZHY(0,NY,NX)/VOLW(0,NY,NX),ZOH(0,NY,NX)/VOLW(0,NY,NX)
C     WRITE(*,1113)'RH2O',I,J,NFZ,NX,NY,TRH2O(0,NY,NX),RH2O 
1113  FORMAT(A8,5I4,80E12.4)
      ELSE
      TRHY(0,NY,NX)=TRHY(0,NY,NX)+(RHY+RHHY)*VOLW(0,NY,NX)
      TROH(0,NY,NX)=TROH(0,NY,NX)+(ROH+RHOH)*VOLW(0,NY,NX)
      TBH2O(0,NY,NX)=TBH2O(0,NY,NX)+RHHY*VOLW(0,NY,NX)
      ENDIF
C
C     CONVERT TOTAL NON-SALT ION TRANSFORMATIONS FROM CHANGES IN 
C     CONCENTRATION TO CHANGES IN MASS (mol m-3 t-1 to mol t-1)
C     FOR USE IN 'REDIST'
C
C     VOLW=litter water content (m3)
C     TRN4S=total NH4 transformation 
C     TRN3S=total NH3 transformation 
C     TRH1P,TRH2P=net HPO4,H2PO4 transformation 
C     TRXN4=total NH4 adsorption 
C     TRALPO,TRFEPO,TRCAPD,TRCAPH,TRCAPM
C        =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation 
C     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=broadcast NH4,NH3,urea,NO3
C        fertilizer
C     RSN4AA=rate of broadcast NH4 fertilizer dissolution 
C     RSN3AA=rate of broadcast NH3 fertilizer dissolution 
C     RSNUAA=rate of broadcast urea fertilizer dissolution 
C     RSNOAA=rate of broadcast NO3 fertilizer dissolution 
C
      TRN4S(0,NY,NX)=TRN4S(0,NY,NX)+RN4S*VOLW(0,NY,NX)
      TRN3S(0,NY,NX)=TRN3S(0,NY,NX)+RN3S*VOLW(0,NY,NX)
      TRH1P(0,NY,NX)=TRH1P(0,NY,NX)+RHP1*VOLW(0,NY,NX)
      TRH2P(0,NY,NX)=TRH2P(0,NY,NX)+RHP2*VOLW(0,NY,NX)
      TRXN4(0,NY,NX)=TRXN4(0,NY,NX)+RXN4*BKVLX
      TRXHY(0,NY,NX)=TRXHY(0,NY,NX)+RXHY*BKVLX
      TRXAL(0,NY,NX)=TRXAL(0,NY,NX)+RXAL*BKVLX
      TRXFE(0,NY,NX)=TRXFE(0,NY,NX)+RXFE*BKVLX
      TRXCA(0,NY,NX)=TRXCA(0,NY,NX)+RXCA*BKVLX
      TRXMG(0,NY,NX)=TRXMG(0,NY,NX)+RXMG*BKVLX
      TRXNA(0,NY,NX)=TRXNA(0,NY,NX)+RXNA*BKVLX
      TRXKA(0,NY,NX)=TRXKA(0,NY,NX)+RXKA*BKVLX
      TRALPO(0,NY,NX)=TRALPO(0,NY,NX)+RPALPX*VOLW(0,NY,NX)
      TRFEPO(0,NY,NX)=TRFEPO(0,NY,NX)+RPFEPX*VOLW(0,NY,NX)
      TRCAPD(0,NY,NX)=TRCAPD(0,NY,NX)+RPCADX*VOLW(0,NY,NX)
      TRCAPH(0,NY,NX)=TRCAPH(0,NY,NX)+RPCAHX*VOLW(0,NY,NX)
      TRCAPM(0,NY,NX)=TRCAPM(0,NY,NX)+RPCAMX*VOLW(0,NY,NX)
      TRAL(0,NY,NX)=TRAL(0,NY,NX)+RAL*VOLW(0,NY,NX)
      TRFE(0,NY,NX)=TRFE(0,NY,NX)+RFE*VOLW(0,NY,NX)
      TRCA(0,NY,NX)=TRCA(0,NY,NX)+RCA*VOLW(0,NY,NX)
      TRMG(0,NY,NX)=TRMG(0,NY,NX)+RMG*VOLW(0,NY,NX)
      TRNA(0,NY,NX)=TRNA(0,NY,NX)+RNA*VOLW(0,NY,NX)
      TRKA(0,NY,NX)=TRKA(0,NY,NX)+RKA*VOLW(0,NY,NX)
      ZNH4FA(0,NY,NX)=ZNH4FA(0,NY,NX)-RSN4AA
      ZNH3FA(0,NY,NX)=ZNH3FA(0,NY,NX)-RSN3AA
      ZNHUFA(0,NY,NX)=ZNHUFA(0,NY,NX)-RSNUAA
      ZNO3FA(0,NY,NX)=ZNO3FA(0,NY,NX)-RSNOAA
      TRN4S(0,NY,NX)=TRN4S(0,NY,NX)+RSN4AA 
      TRN3S(0,NY,NX)=TRN3S(0,NY,NX)+RSN3AA+RSNUAA 
      TRNO3(0,NY,NX)=TRNO3(0,NY,NX)+RSNOAA
C
C     CONVERT C,N,TO MASS UNITS
C
C     TRN4S,TRN3S,TRNO3=net NH4,NH3,NO3 transformation (g N t-1) 
C     TRH1P,TRH2P=net HPO4,H2PO4 transformation (g P t-1) 
C
      TRCO2(0,NY,NX)=TRCO2(0,NY,NX)*12.0
      TRN4S(0,NY,NX)=TRN4S(0,NY,NX)*14.0
      TRN3S(0,NY,NX)=TRN3S(0,NY,NX)*14.0
      TRNO3(0,NY,NX)=TRNO3(0,NY,NX)*14.0
      TRH1P(0,NY,NX)=TRH1P(0,NY,NX)*31.0
      TRH2P(0,NY,NX)=TRH2P(0,NY,NX)*31.0
C     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,9989)'TROH0',I,J,NFZ,TROH(0,NY,NX)
C    2,ROH,RHOH,VOLW(0,NY,NX),AOH2,AOH1,DPH2O,A1(0,NY,NX),AHY2
C    3,PH(0,NY,NX)
C     WRITE(*,9989)'TRN3S0',I,J,NFZ,TRN3S(0,NY,NX)
C    2,RN3S,VOLW(0,NY,NX),RSN3AA,RSNUAA 
C     WRITE(*,9989)'TRN4S0',I,J,NFZ,TRN4S(0,NY,NX),TRXN4(0,NY,NX)
C    2,RN4S,RNH4,RXN4,RSN4AA,VOLW(0,NY,NX)
C    3,SPNH4,ZNH4FA(0,NY,NX),THETW0
C     WRITE(*,9989)'TRXHY0',I,J,NFZ,TRXHY(0,NY,NX),RXHY 
C    2,TRHY(0,NY,NX),RHY,RHHY,VOLW(0,NY,NX)
C     WRITE(*,9989)'RHP10',I,J,NFZ,TRH1P(0,NY,NX)
C    2,RHP1,RH2P,RPCADX,3.0*RPCAHX
C     WRITE(*,9989)'TRCA0',I,J,NFZ,TRCA(0,NY,NX)
C    2,RCA,RCA*VOLW(0,NY,NX),RXCA,RPCADX,5.0*RPCAHX,RPCAMX
C     WRITE(*,9989)'TRAL0',I,J,NFZ,TRAL(0,NY,NX)
C    2,RAL,RXAL*BKVLW,RPALPX,RXAL,BKVLW,ZAL(0,NY,NX) 
9989  FORMAT(A8,3I4,20E12.4)
C     ENDIF
9990  CONTINUE
9995  CONTINUE
      RETURN
      END


