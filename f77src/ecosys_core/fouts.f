      SUBROUTINE fouts(NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE OPENS AND LABELS OUTPUT FILES FOR SOIL DATA
C
      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"
      include "blk17.h"
      CHARACTER*16 DATA(30),DATAC(30,250,250),DATAP(JP,JY,JX)
     2,DATAM(JP,JY,JX),DATAX(JP),DATAY(JP),DATAZ(JP,JY,JX)
     3,OUTS(10),OUTP(10),OUTFILS(10,JY,JX),OUTFILP(10,JP,JY,JX)
      CHARACTER*3 CHOICE(102,20)
      CHARACTER*8 CDATE
      CHARACTER*8 CDOY,DATE,HOUR
      CHARACTER*1 CHARM,CHARN,CHARZ,CHARO
      CHARACTER*2 CHARX,CHARY
      CHARACTER*4 CHARR
      CHARACTER*16 HEAD(50)
      CHARACTER*80 PREFIX
      CDOY='DOY     '
      DATE='DATE    '
      HOUR='HOUR    '
      CHARM='0'
      CHARO='0'
      CHARZ='0'
      DO 8995 NX=NHW,NHE
      DO 8995 NY=NVN,NVS
      DO 8995 N=1,10
      OUTFILS(N,NY,NX)= '                '
8995  CONTINUE
C
C     OPEN AND NAME OUTPUT FILES
C
      DO 1010 N=21,30
      IF(DATAC(N,NE,NEX).NE.'NO')THEN
      OPEN(15,FILE=TRIM(PREFIX)//DATAC(N,NE,NEX),STATUS='OLD')
      WRITE(CHARR,'(I4)')IYRC
      OUTS(N-20)=CHARO//CHARR//DATAC(N,NE,NEX)
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      IF(NX.LT.10)THEN
      WRITE(CHARN,'(I1)')NX
      CHARX=CHARM//CHARN
      ELSE
      WRITE(CHARX,'(I2)')NX
      ENDIF
      IF(NY.LT.10)THEN
      WRITE(CHARN,'(I1)')NY
      CHARY=CHARM//CHARN
      ELSE
      WRITE(CHARY,'(I2)')NY
      ENDIF
      OUTFILS(N-20,NY,NX)=CHARX//CHARY//CHARZ//CHARR
     2//DATAC(N,NE,NEX)
9990  CONTINUE
9995  CONTINUE
      DO 4010 M=1,2
      READ(15,*)DY
      IDY1=INT(DY/1.0E+02)
      IDY2=INT(DY/1.0E+00-IDY1*1.0E+02)
      IF(IDY2.GT.2)LPY=1
      IF(IDY2.EQ.1)GO TO 4525
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
      GO TO 4530
4525  IDY=IDY1
4530  IF(M.EQ.1)IDATA(N)=IDY
      IF(M.EQ.2)IDATA(N+20)=IDY
4010  CONTINUE
      M=0
4020  CONTINUE
      M=M+1
      READ(15,'(A3)',END=4030)CHOICE(M,N-20)
      GO TO 4020
4030  CONTINUE
      CLOSE(15)
      LUN=N+10
      OPEN(LUN,FILE=trim(outdir)//OUTS(N-20),STATUS='UNKNOWN')
C
C     WRITE HEADINGS TO OUTPUT FILES
C
      M=0
      IF(N.EQ.21)THEN
      DO 1021 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='SOIL_CO2_FLUX'
      IF(L.EQ.2)HEAD(M)='ECO_CO2_FLUX'
      IF(L.EQ.3)HEAD(M)='CH4_FLUX'
      IF(L.EQ.4)HEAD(M)='O2_FLUX'
      IF(L.EQ.5)HEAD(M)='CO2_1'
      IF(L.EQ.6)HEAD(M)='CO2_2'
      IF(L.EQ.7)HEAD(M)='CO2_3'
      IF(L.EQ.8)HEAD(M)='CO2_4'
      IF(L.EQ.9)HEAD(M)='CO2_5'
      IF(L.EQ.10)HEAD(M)='CO2_6'
      IF(L.EQ.11)HEAD(M)='CO2_7'
      IF(L.EQ.12)HEAD(M)='CO2_8'
      IF(L.EQ.13)HEAD(M)='CO2_9'
      IF(L.EQ.14)HEAD(M)='CO2_10'
      IF(L.EQ.15)HEAD(M)='CO2_11'
      IF(L.EQ.16)HEAD(M)='CO2_12'
      IF(L.EQ.17)HEAD(M)='CO2_13'
      IF(L.EQ.18)HEAD(M)='CO2_14'
      IF(L.EQ.19)HEAD(M)='CO2_LIT'
      IF(L.EQ.20)HEAD(M)='CH4_1'
      IF(L.EQ.21)HEAD(M)='CH4_2'
      IF(L.EQ.22)HEAD(M)='CH4_3'
      IF(L.EQ.23)HEAD(M)='CH4_4'
      IF(L.EQ.24)HEAD(M)='CH4_5'
      IF(L.EQ.25)HEAD(M)='CH4_6'
      IF(L.EQ.26)HEAD(M)='CH4_7'
      IF(L.EQ.27)HEAD(M)='CH4_8'
      IF(L.EQ.28)HEAD(M)='CH4_9'
      IF(L.EQ.29)HEAD(M)='CH4_10'
      IF(L.EQ.30)HEAD(M)='CH4_11'
      IF(L.EQ.31)HEAD(M)='CH4_12'
      IF(L.EQ.32)HEAD(M)='CH4_13'
      IF(L.EQ.33)HEAD(M)='CH4_14'
      IF(L.EQ.34)HEAD(M)='CH4_15'
      IF(L.EQ.35)HEAD(M)='O2_1'
      IF(L.EQ.36)HEAD(M)='O2_2'
      IF(L.EQ.37)HEAD(M)='O2_3'
      IF(L.EQ.38)HEAD(M)='O2_4'
      IF(L.EQ.39)HEAD(M)='O2_5'
      IF(L.EQ.40)HEAD(M)='O2_6'
      IF(L.EQ.41)HEAD(M)='O2_7'
      IF(L.EQ.42)HEAD(M)='O2_8'
      IF(L.EQ.43)HEAD(M)='O2_9'
      IF(L.EQ.44)HEAD(M)='O2_10'
      IF(L.EQ.45)HEAD(M)='O2_11'
      IF(L.EQ.46)HEAD(M)='O2_12'
      IF(L.EQ.47)HEAD(M)='O2_13'
      IF(L.EQ.48)HEAD(M)='O2_14'
      IF(L.EQ.49)HEAD(M)='O2_15'
      IF(L.EQ.50)HEAD(M)='O2_LIT'
      ENDIF
1021  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.EQ.22)THEN
      DO 1022 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='EVAPN'
      IF(L.EQ.2)HEAD(M)='RUNOFF'
      IF(L.EQ.3)HEAD(M)='SEDIMENT'
      IF(L.EQ.4)HEAD(M)='TTL_SWC'
      IF(L.EQ.5)HEAD(M)='DISCHG'
      IF(L.EQ.6)HEAD(M)='SNOWPACK'
      IF(L.EQ.7)HEAD(M)='WTR_1'
      IF(L.EQ.8)HEAD(M)='WTR_2'
      IF(L.EQ.9)HEAD(M)='WTR_3'
      IF(L.EQ.10)HEAD(M)='WTR_4'
      IF(L.EQ.11)HEAD(M)='WTR_5'
      IF(L.EQ.12)HEAD(M)='WTR_6'
      IF(L.EQ.13)HEAD(M)='WTR_7'
      IF(L.EQ.14)HEAD(M)='WTR_8'
      IF(L.EQ.15)HEAD(M)='WTR_9'
      IF(L.EQ.16)HEAD(M)='WTR_10'
      IF(L.EQ.17)HEAD(M)='WTR_11'
      IF(L.EQ.18)HEAD(M)='WTR_12'
      IF(L.EQ.19)HEAD(M)='WTR_13'
      IF(L.EQ.20)HEAD(M)='WTR_14'
      IF(L.EQ.21)HEAD(M)='WTR_15'
      IF(L.EQ.22)HEAD(M)='WTR_16'
      IF(L.EQ.23)HEAD(M)='WTR_17'
      IF(L.EQ.24)HEAD(M)='WTR_18'
      IF(L.EQ.25)HEAD(M)='WTR_19'
      IF(L.EQ.26)HEAD(M)='WTR_20'
      IF(L.EQ.27)HEAD(M)='SURF_WTR'
      IF(L.EQ.28)HEAD(M)='ICE_1'
      IF(L.EQ.29)HEAD(M)='ICE_2'
      IF(L.EQ.30)HEAD(M)='ICE_3'
      IF(L.EQ.31)HEAD(M)='ICE_4'
      IF(L.EQ.32)HEAD(M)='ICE_5'
      IF(L.EQ.33)HEAD(M)='ICE_6'
      IF(L.EQ.34)HEAD(M)='ICE_7'
      IF(L.EQ.35)HEAD(M)='ICE_8'
      IF(L.EQ.36)HEAD(M)='ICE_9'
      IF(L.EQ.37)HEAD(M)='ICE_10'
      IF(L.EQ.38)HEAD(M)='ICE_11'
      IF(L.EQ.39)HEAD(M)='ICE_12'
      IF(L.EQ.40)HEAD(M)='ICE_13'
      IF(L.EQ.41)HEAD(M)='ICE_14'
      IF(L.EQ.42)HEAD(M)='ICE_15'
      IF(L.EQ.43)HEAD(M)='ICE_16'
      IF(L.EQ.44)HEAD(M)='ICE_17'
      IF(L.EQ.45)HEAD(M)='ICE_18'
      IF(L.EQ.46)HEAD(M)='ICE_19'
      IF(L.EQ.47)HEAD(M)='ICE_20'
      IF(L.EQ.48)HEAD(M)='SURF_ICE'
      IF(L.EQ.49)HEAD(M)='ACTV_LYR'
      IF(L.EQ.50)HEAD(M)='WTR_TBL'
      ENDIF
1022  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.EQ.23)THEN
      DO 1023 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='N2O_FLUX'
      IF(L.EQ.2)HEAD(M)='N2G_FLUX'
      IF(L.EQ.3)HEAD(M)='NH3_FLUX'
      IF(L.EQ.4)HEAD(M)='SURF_N_FLUX'
      IF(L.EQ.5)HEAD(M)='SUBS_N_FLUX'
      IF(L.EQ.6)HEAD(M)='N2O_1'
      IF(L.EQ.7)HEAD(M)='N2O_2'
      IF(L.EQ.8)HEAD(M)='N2O_3'
      IF(L.EQ.9)HEAD(M)='N2O_4'
      IF(L.EQ.10)HEAD(M)='N2O_5'
      IF(L.EQ.11)HEAD(M)='N2O_6'
      IF(L.EQ.12)HEAD(M)='N2O_7'
      IF(L.EQ.13)HEAD(M)='N2O_8'
      IF(L.EQ.14)HEAD(M)='N2O_9'
      IF(L.EQ.15)HEAD(M)='N2O_10'
      IF(L.EQ.16)HEAD(M)='N2O_11'
      IF(L.EQ.17)HEAD(M)='N2O_12'
      IF(L.EQ.18)HEAD(M)='N2O_13'
      IF(L.EQ.19)HEAD(M)='N2O_14'
      IF(L.EQ.20)HEAD(M)='N2O_15'
      IF(L.EQ.21)HEAD(M)='NH3_1'
      IF(L.EQ.22)HEAD(M)='NH3_2'
      IF(L.EQ.23)HEAD(M)='NH3_3'
      IF(L.EQ.24)HEAD(M)='NH3_4'
      IF(L.EQ.25)HEAD(M)='NH3_5'
      IF(L.EQ.26)HEAD(M)='NH3_6'
      IF(L.EQ.27)HEAD(M)='NH3_7'
      IF(L.EQ.28)HEAD(M)='NH3_8'
      IF(L.EQ.29)HEAD(M)='NH3_9'
      IF(L.EQ.30)HEAD(M)='NH3_10'
      IF(L.EQ.31)HEAD(M)='NH3_11'
      IF(L.EQ.32)HEAD(M)='NH3_12'
      IF(L.EQ.33)HEAD(M)='NH3_13'
      IF(L.EQ.34)HEAD(M)='NH3_14'
      IF(L.EQ.35)HEAD(M)='NH3_15'
      IF(L.EQ.36)HEAD(M)='N2O_LIT'
      IF(L.EQ.37)HEAD(M)='NH3_LIT'
      ENDIF
1023  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.EQ.24)THEN
      DO 1024 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='SURF_P_FLUX'
      IF(L.EQ.2)HEAD(M)='SUBS_P_FLUX'
      ENDIF
1024  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.EQ.25)THEN
      DO 1025 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='SOL_RADN'
      IF(L.EQ.2)HEAD(M)='AIR_TEMP'
      IF(L.EQ.3)HEAD(M)='HUM'
      IF(L.EQ.4)HEAD(M)='WIND'
      IF(L.EQ.5)HEAD(M)='PREC'
      IF(L.EQ.6)HEAD(M)='SOIL_RN'
      IF(L.EQ.7)HEAD(M)='SOIL_LE'
      IF(L.EQ.8)HEAD(M)='SOIL_H'
      IF(L.EQ.9)HEAD(M)='SOIL_G'
      IF(L.EQ.10)HEAD(M)='ECO_RN'
      IF(L.EQ.11)HEAD(M)='ECO_LE'
      IF(L.EQ.12)HEAD(M)='ECO_H'
      IF(L.EQ.13)HEAD(M)='ECO_G'
      IF(L.EQ.14)HEAD(M)='TEMP_1'
      IF(L.EQ.15)HEAD(M)='TEMP_2'
      IF(L.EQ.16)HEAD(M)='TEMP_3'
      IF(L.EQ.17)HEAD(M)='TEMP_4'
      IF(L.EQ.18)HEAD(M)='TEMP_5'
      IF(L.EQ.19)HEAD(M)='TEMP_6'
      IF(L.EQ.20)HEAD(M)='TEMP_7'
      IF(L.EQ.21)HEAD(M)='TEMP_8'
      IF(L.EQ.22)HEAD(M)='TEMP_9'
      IF(L.EQ.23)HEAD(M)='TEMP_10'
      IF(L.EQ.24)HEAD(M)='TEMP_11'
      IF(L.EQ.25)HEAD(M)='TEMP_12'
      IF(L.EQ.26)HEAD(M)='TEMP_13'
      IF(L.EQ.27)HEAD(M)='TEMP_14'
      IF(L.EQ.28)HEAD(M)='TEMP_15'
      IF(L.EQ.29)HEAD(M)='TEMP_16'
      IF(L.EQ.30)HEAD(M)='TEMP_17'
      IF(L.EQ.31)HEAD(M)='TEMP_18'
      IF(L.EQ.32)HEAD(M)='TEMP_19'
      IF(L.EQ.33)HEAD(M)='TEMP_20'
      IF(L.EQ.34)HEAD(M)='TEMP_LITTER'
      IF(L.EQ.35)HEAD(M)='TEMP_SNOW'
      ENDIF
1025  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.EQ.26)THEN
      DO 1026 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='RESIDUE_C'
      IF(L.EQ.2)HEAD(M)='HUMUS_C'
      IF(L.EQ.3)HEAD(M)='AMENDED_C'
      IF(L.EQ.4)HEAD(M)='LITTER_C'
      IF(L.EQ.5)HEAD(M)='CO2_FLUX'
      IF(L.EQ.6)HEAD(M)='O2_FLUX'
      IF(L.EQ.7)HEAD(M)='AUTO_RESP'
      IF(L.EQ.8)HEAD(M)='MICRO_C'
      IF(L.EQ.9)HEAD(M)='SURF_RES'
      IF(L.EQ.10)HEAD(M)='CH4_FLUX'
      IF(L.EQ.11)HEAD(M)='SUR_DOC_SED_FLX'
      IF(L.EQ.12)HEAD(M)='SUB_DOC_FLX'
      IF(L.EQ.13)HEAD(M)='SUR_DIC_FLX'
      IF(L.EQ.14)HEAD(M)='SUB_DIC_FLX'
      IF(L.EQ.15)HEAD(M)='ATM_CO2'
      IF(L.EQ.16)HEAD(M)='NBP'
      IF(L.EQ.17)HEAD(M)='FIRE_CO2'
      IF(L.EQ.18)HEAD(M)='SOC_1'
      IF(L.EQ.19)HEAD(M)='SOC_2'
      IF(L.EQ.20)HEAD(M)='SOC_3'
      IF(L.EQ.21)HEAD(M)='SOC_4'
      IF(L.EQ.22)HEAD(M)='SOC_5'
      IF(L.EQ.23)HEAD(M)='SOC_6'
      IF(L.EQ.24)HEAD(M)='SOC_7'
      IF(L.EQ.25)HEAD(M)='SOC_8'
      IF(L.EQ.26)HEAD(M)='SOC_9'
      IF(L.EQ.27)HEAD(M)='SOC_10'
      IF(L.EQ.28)HEAD(M)='SOC_11'
      IF(L.EQ.29)HEAD(M)='SOC_12'
      IF(L.EQ.30)HEAD(M)='SOC_13'
      IF(L.EQ.31)HEAD(M)='SOC_14'
      IF(L.EQ.32)HEAD(M)='SOC_15'
      IF(L.EQ.41)HEAD(M)='H2_FLUX'
      IF(L.EQ.42)HEAD(M)='ECO_HVST_C'
      IF(L.EQ.43)HEAD(M)='ECO_LAI'
      IF(L.EQ.44)HEAD(M)='ECO_GPP'
      IF(L.EQ.45)HEAD(M)='ECO_RA'
      IF(L.EQ.46)HEAD(M)='ECO_NPP'
      IF(L.EQ.47)HEAD(M)='ECO_RH'
      IF(L.EQ.48)HEAD(M)='FIRE_CH4'
      IF(L.EQ.49)HEAD(M)='TTL_DIC'
      IF(L.EQ.50)HEAD(M)='STG_DEAD'
      ENDIF
1026  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.EQ.27)THEN
      DO 1027 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='PRECN'
      IF(L.EQ.2)HEAD(M)='ET'
      IF(L.EQ.3)HEAD(M)='RUNOFF'
      IF(L.EQ.4)HEAD(M)='WATER'
      IF(L.EQ.5)HEAD(M)='DISCHG'
      IF(L.EQ.6)HEAD(M)='SNOWPACK'
      IF(L.EQ.7)HEAD(M)='WTR_1'
      IF(L.EQ.8)HEAD(M)='WTR_2'
      IF(L.EQ.9)HEAD(M)='WTR_3'
      IF(L.EQ.10)HEAD(M)='WTR_4'
      IF(L.EQ.11)HEAD(M)='WTR_5'
      IF(L.EQ.12)HEAD(M)='WTR_6'
      IF(L.EQ.13)HEAD(M)='WTR_7'
      IF(L.EQ.14)HEAD(M)='WTR_8'
      IF(L.EQ.15)HEAD(M)='WTR_9'
      IF(L.EQ.16)HEAD(M)='WTR_10'
      IF(L.EQ.17)HEAD(M)='WTR_11'
      IF(L.EQ.18)HEAD(M)='WTR_12'
      IF(L.EQ.19)HEAD(M)='WTR_13'
      IF(L.EQ.20)HEAD(M)='SURF_WTR'
      IF(L.EQ.21)HEAD(M)='ICE_1'
      IF(L.EQ.22)HEAD(M)='ICE_2'
      IF(L.EQ.23)HEAD(M)='ICE_3'
      IF(L.EQ.24)HEAD(M)='ICE_4'
      IF(L.EQ.25)HEAD(M)='ICE_5'
      IF(L.EQ.26)HEAD(M)='ICE_6'
      IF(L.EQ.27)HEAD(M)='ICE_7'
      IF(L.EQ.28)HEAD(M)='ICE_8'
      IF(L.EQ.29)HEAD(M)='ICE_9'
      IF(L.EQ.30)HEAD(M)='ICE_10'
      IF(L.EQ.31)HEAD(M)='ICE_11'
      IF(L.EQ.32)HEAD(M)='ICE_12'
      IF(L.EQ.33)HEAD(M)='ICE_13'
      IF(L.EQ.34)HEAD(M)='SURF_ICE'
      IF(L.EQ.35)HEAD(M)='PSI_1'
      IF(L.EQ.36)HEAD(M)='PSI_2'
      IF(L.EQ.37)HEAD(M)='PSI_3'
      IF(L.EQ.38)HEAD(M)='PSI_4'
      IF(L.EQ.39)HEAD(M)='PSI_5'
      IF(L.EQ.40)HEAD(M)='PSI_6'
      IF(L.EQ.41)HEAD(M)='PSI_7'
      IF(L.EQ.42)HEAD(M)='PSI_8'
      IF(L.EQ.43)HEAD(M)='PSI_9'
      IF(L.EQ.44)HEAD(M)='PSI_10'
      IF(L.EQ.45)HEAD(M)='PSI_11'
      IF(L.EQ.46)HEAD(M)='SEDIMENT'
      IF(L.EQ.47)HEAD(M)='PSI_SURF'
      IF(L.EQ.48)HEAD(M)='SURF_ELEV'
      IF(L.EQ.49)HEAD(M)='ACTV_LYR'
      IF(L.EQ.50)HEAD(M)='WTR_TBL'
      ENDIF
1027  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.EQ.28)THEN
      DO 1028 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='RESIDUE_N'
      IF(L.EQ.2)HEAD(M)='HUMUS_N'
      IF(L.EQ.3)HEAD(M)='FERTZR_N'
      IF(L.EQ.4)HEAD(M)='NET_PL_EXCH_N'
      IF(L.EQ.5)HEAD(M)='NH4'
      IF(L.EQ.6)HEAD(M)='NO3'
      IF(L.EQ.7)HEAD(M)='SUR_DON_SED_FLX'
      IF(L.EQ.8)HEAD(M)='SUB_DON_FLX'
      IF(L.EQ.9)HEAD(M)='SUR_DIN_FLX'
      IF(L.EQ.10)HEAD(M)='SUB_DIN_FLX'
      IF(L.EQ.11)HEAD(M)='N2O_FLUX'
      IF(L.EQ.12)HEAD(M)='NH3_FLUX'
      IF(L.EQ.13)HEAD(M)='N2_FIXN'
      IF(L.EQ.14)HEAD(M)='MICRO_N'
      IF(L.EQ.15)HEAD(M)='NH4_1'
      IF(L.EQ.16)HEAD(M)='NH4_2'
      IF(L.EQ.17)HEAD(M)='NH4_3'
      IF(L.EQ.18)HEAD(M)='NH4_4'
      IF(L.EQ.19)HEAD(M)='NH4_5'
      IF(L.EQ.20)HEAD(M)='NH4_6'
      IF(L.EQ.21)HEAD(M)='NH4_7'
      IF(L.EQ.22)HEAD(M)='NH4_8'
      IF(L.EQ.23)HEAD(M)='NH4_9'
      IF(L.EQ.24)HEAD(M)='NH4_10'
      IF(L.EQ.25)HEAD(M)='NH4_11'
      IF(L.EQ.26)HEAD(M)='NH4_12'
      IF(L.EQ.27)HEAD(M)='NH4_13'
      IF(L.EQ.28)HEAD(M)='NH4_14'
      IF(L.EQ.29)HEAD(M)='NH4_15'
      IF(L.EQ.30)HEAD(M)='NO3_1'
      IF(L.EQ.31)HEAD(M)='NO3_2'
      IF(L.EQ.32)HEAD(M)='NO3_3'
      IF(L.EQ.33)HEAD(M)='NO3_4'
      IF(L.EQ.34)HEAD(M)='NO3_5'
      IF(L.EQ.35)HEAD(M)='NO3_6'
      IF(L.EQ.36)HEAD(M)='NO3_7'
      IF(L.EQ.37)HEAD(M)='NO3_8'
      IF(L.EQ.38)HEAD(M)='NO3_9'
      IF(L.EQ.39)HEAD(M)='NO3_10'
      IF(L.EQ.40)HEAD(M)='NO3_11'
      IF(L.EQ.41)HEAD(M)='NO3_12'
      IF(L.EQ.42)HEAD(M)='NO3_13'
      IF(L.EQ.43)HEAD(M)='NO3_14'
      IF(L.EQ.44)HEAD(M)='NO3_15'
      IF(L.EQ.45)HEAD(M)='NH4_RES'
      IF(L.EQ.46)HEAD(M)='NO3_RES'
      IF(L.EQ.47)HEAD(M)='ECO_HVST_N'
      IF(L.EQ.48)HEAD(M)='NET_N_MIN'
      IF(L.EQ.49)HEAD(M)='FIRE_N'
      IF(L.EQ.50)HEAD(M)='N2_FLUX'
      ENDIF
1028  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.EQ.29)THEN
      DO 1029 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='RESIDUE_P'
      IF(L.EQ.2)HEAD(M)='HUMUS_P'
      IF(L.EQ.3)HEAD(M)='FERTZR_P'
      IF(L.EQ.4)HEAD(M)='NET_PL_EXCH_P'
      IF(L.EQ.5)HEAD(M)='EXCH_PO4'
      IF(L.EQ.6)HEAD(M)='SUR_DOP_SED_FLX'
      IF(L.EQ.7)HEAD(M)='SUB_DOP_FLX'
      IF(L.EQ.8)HEAD(M)='SUR_DIP_FLX'
      IF(L.EQ.9)HEAD(M)='SUB_DIP_FLX'
      IF(L.EQ.10)HEAD(M)='PRECIP_P'
      IF(L.EQ.11)HEAD(M)='MICRO_P'
      IF(L.EQ.12)HEAD(M)='FIRE_P'
      IF(L.EQ.13)HEAD(M)='PO4_1'
      IF(L.EQ.14)HEAD(M)='PO4_2'
      IF(L.EQ.15)HEAD(M)='PO4_3'
      IF(L.EQ.16)HEAD(M)='PO4_4'
      IF(L.EQ.17)HEAD(M)='PO4_5'
      IF(L.EQ.18)HEAD(M)='PO4_6'
      IF(L.EQ.19)HEAD(M)='PO4_7'
      IF(L.EQ.20)HEAD(M)='PO4_8'
      IF(L.EQ.21)HEAD(M)='PO4_9'
      IF(L.EQ.22)HEAD(M)='PO4_10'
      IF(L.EQ.23)HEAD(M)='PO4_11'
      IF(L.EQ.24)HEAD(M)='PO4_12'
      IF(L.EQ.25)HEAD(M)='PO4_13'
      IF(L.EQ.26)HEAD(M)='PO4_14'
      IF(L.EQ.27)HEAD(M)='PO4_15'
      IF(L.EQ.28)HEAD(M)='EXCH_P_1'
      IF(L.EQ.29)HEAD(M)='EXCH_P_2'
      IF(L.EQ.30)HEAD(M)='EXCH_P_3'
      IF(L.EQ.31)HEAD(M)='EXCH_P_4'
      IF(L.EQ.32)HEAD(M)='EXCH_P_5'
      IF(L.EQ.33)HEAD(M)='EXCH_P_6'
      IF(L.EQ.34)HEAD(M)='EXCH_P_7'
      IF(L.EQ.35)HEAD(M)='EXCH_P_8'
      IF(L.EQ.36)HEAD(M)='EXCH_P_9'
      IF(L.EQ.37)HEAD(M)='EXCH_P_10'
      IF(L.EQ.38)HEAD(M)='EXCH_P_11'
      IF(L.EQ.39)HEAD(M)='EXCH_P_12'
      IF(L.EQ.40)HEAD(M)='EXCH_P_13'
      IF(L.EQ.41)HEAD(M)='EXCH_P_14'
      IF(L.EQ.42)HEAD(M)='EXCH_P_15'
      IF(L.EQ.43)HEAD(M)='PO4_RES'
      IF(L.EQ.44)HEAD(M)='EXCH_P_RES'
      IF(L.EQ.47)HEAD(M)='ECO_HVST_P'
      IF(L.EQ.48)HEAD(M)='NET_P_MIN'
      ENDIF
1029  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.EQ.30)THEN
      DO 1030 L=1,50
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.1)HEAD(M)='RADN'
      IF(L.EQ.2)HEAD(M)='TMAX_AIR'
      IF(L.EQ.3)HEAD(M)='TMIN_AIR'
      IF(L.EQ.4)HEAD(M)='HMAX_AIR'
      IF(L.EQ.5)HEAD(M)='HMIN_AIR'
      IF(L.EQ.6)HEAD(M)='WIND'
      IF(L.EQ.7)HEAD(M)='PRECN'
      IF(L.EQ.8)HEAD(M)='TMAX_SOIL_1'
      IF(L.EQ.9)HEAD(M)='TMIN_SOIL_1'
      IF(L.EQ.10)HEAD(M)='TMAX_SOIL_2'
      IF(L.EQ.11)HEAD(M)='TMIN_SOIL_2'
      IF(L.EQ.12)HEAD(M)='TMAX_SOIL_3'
      IF(L.EQ.13)HEAD(M)='TMIN_SOIL_3'
      IF(L.EQ.14)HEAD(M)='TMAX_SOIL_4'
      IF(L.EQ.15)HEAD(M)='TMIN_SOIL_4'
      IF(L.EQ.16)HEAD(M)='TMAX_SOIL_5'
      IF(L.EQ.17)HEAD(M)='TMIN_SOIL_5'
      IF(L.EQ.18)HEAD(M)='TMAX_SOIL_6'
      IF(L.EQ.19)HEAD(M)='TMIN_SOIL_6'
      IF(L.EQ.20)HEAD(M)='TMAX_SOIL_7'
      IF(L.EQ.21)HEAD(M)='TMIN_SOIL_7'
      IF(L.EQ.22)HEAD(M)='TMAX_SOIL_8'
      IF(L.EQ.23)HEAD(M)='TMIN_SOIL_8'
      IF(L.EQ.24)HEAD(M)='TMAX_SOIL_9'
      IF(L.EQ.25)HEAD(M)='TMIN_SOIL_9'
      IF(L.EQ.26)HEAD(M)='TMAX_SOIL_10'
      IF(L.EQ.27)HEAD(M)='TMIN_SOIL_10'
      IF(L.EQ.28)HEAD(M)='TMAX_SOIL_11'
      IF(L.EQ.29)HEAD(M)='TMIN_SOIL_11'
      IF(L.EQ.30)HEAD(M)='TMAX_SOIL_12'
      IF(L.EQ.31)HEAD(M)='TMIN_SOIL_12'
      IF(L.EQ.32)HEAD(M)='TMAX_SOIL_13'
      IF(L.EQ.33)HEAD(M)='TMIN_SOIL_13'
      IF(L.EQ.34)HEAD(M)='TMAX_SOIL_14'
      IF(L.EQ.35)HEAD(M)='TMIN_SOIL_14'
      IF(L.EQ.36)HEAD(M)='TMAX_LITTER'
      IF(L.EQ.37)HEAD(M)='TMIN_LITTER'
      IF(L.EQ.38)HEAD(M)='ECND_1'
      IF(L.EQ.39)HEAD(M)='ECND_2'
      IF(L.EQ.40)HEAD(M)='ECND_3'
      IF(L.EQ.41)HEAD(M)='ECND_4'
      IF(L.EQ.42)HEAD(M)='ECND_5'
      IF(L.EQ.43)HEAD(M)='ECND_6'
      IF(L.EQ.44)HEAD(M)='ECND_7'
      IF(L.EQ.45)HEAD(M)='ECND_8'
      IF(L.EQ.46)HEAD(M)='ECND_9'
      IF(L.EQ.47)HEAD(M)='ECND_10'
      IF(L.EQ.48)HEAD(M)='ECND_11'
      IF(L.EQ.49)HEAD(M)='ECND_12'
      IF(L.EQ.50)HEAD(M)='TTL_SALT_DISCHG'
      ENDIF
1030  CONTINUE
      NOUTS(N-20)=M
      ENDIF
      IF(N.LE.25)WRITE(LUN,'(A12,2A8,50A16)')CDOY,DATE,HOUR
     2,(HEAD(K),K=1,M)
      IF(N.GE.26)WRITE(LUN,'(A12,A8,50A16)')CDOY,DATE
     2,(HEAD(K),K=1,M)
      ENDIF
1010  CONTINUE
      RETURN
      END