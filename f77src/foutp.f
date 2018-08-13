
      SUBROUTINE foutp(NT,NE,NAX,NDX,NTX,NEX,NF,NFX,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE OPENS AND LABELS OUTPUT FILES
C     FOR PLANT DATA
C
      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"
      include "blk3.h"
      SAVE LUN
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
      CHARO='1'
      DO 8995 NX=NHW,NHE
      DO 8995 NY=NVN,NVS
      DO 8995 NZ=1,NP0(NY,NX)
      DO 8995 N=1,10
      OUTFILP(N,NZ,NY,NX)= '                '
8995  CONTINUE
C
C     OPEN AND NAME OUTPUT FILES
C
      DO 1010 N=21,30
      IF(DATAC(N,NE,NEX).NE.'NO')THEN
      WRITE(CHARR,'(I4)')IYRC
      OUTP(N-20)=CHARO//CHARR//DATAC(N,NE,NEX)
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      DO 9985 NZ=1,NP0(NY,NX)
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
      WRITE(CHARZ,'(I1)')NZ
      OUTFILP(N-20,NZ,NY,NX)=CHARX//CHARY//CHARZ//CHARR
     2//DATAC(N,NE,NEX)
9985  CONTINUE
9990  CONTINUE
9995  CONTINUE
      LUN=N+20
C      OPEN(LUN,FILE=OUTP(N-20),STATUS='UNKNOWN')
      OPEN(LUN,FILE=trim(outdir)//OUTP(N-20),STATUS='UNKNOWN')   
C
C     WRITE HEADINGS TO OUTPUT FILES
C
      M=0
      IF(N.EQ.21)THEN
      DO 1021 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='CAN_CO2_FLUX'
      IF(L.EQ.52)HEAD(M)='CAN_GPP'
      IF(L.EQ.53)HEAD(M)='CAN_RA'
      IF(L.EQ.54)HEAD(M)='[TNC]'
      IF(L.EQ.55)HEAD(M)='STOML_RSC'
      IF(L.EQ.56)HEAD(M)='BLYR_RSC'
      IF(L.EQ.57)HEAD(M)='[CAN_CO2]'
      IF(L.EQ.58)HEAD(M)='LAI'
      ENDIF
1021  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.EQ.22)THEN
      DO 1022 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='PSI_CAN'
      IF(L.EQ.52)HEAD(M)='TURG_CAN'
      IF(L.EQ.53)HEAD(M)='STOM_RSC'
      IF(L.EQ.54)HEAD(M)='BLYR_RSC'
      IF(L.EQ.55)HEAD(M)='TRANSPN'
      IF(L.EQ.56)HEAD(M)='O2_STRESS'
      IF(L.EQ.57)HEAD(M)='PSI_RT_1'
      IF(L.EQ.58)HEAD(M)='PSI_RT_2'
      IF(L.EQ.59)HEAD(M)='PSI_RT_3'
      IF(L.EQ.60)HEAD(M)='PSI_RT_4'
      IF(L.EQ.61)HEAD(M)='PSI_RT_5'
      IF(L.EQ.62)HEAD(M)='PSI_RT_6'
      IF(L.EQ.63)HEAD(M)='PSI_RT_7'
      IF(L.EQ.64)HEAD(M)='PSI_RT_8'
      IF(L.EQ.65)HEAD(M)='PSI_RT_9'
      IF(L.EQ.66)HEAD(M)='PSI_RT_10'
      IF(L.EQ.67)HEAD(M)='PSI_RT_11'
      IF(L.EQ.68)HEAD(M)='PSI_RT_12'
      IF(L.EQ.69)HEAD(M)='PSI_RT_13'
      IF(L.EQ.70)HEAD(M)='PSI_RT_14'
      IF(L.EQ.71)HEAD(M)='PSI_RT_15'
      ENDIF
1022  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.EQ.23)THEN
      DO 1023 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='NH4_UPTK'
      IF(L.EQ.52)HEAD(M)='NO3_UPTK'
      IF(L.EQ.53)HEAD(M)='N2_FIXN'
      IF(L.EQ.54)HEAD(M)='[TNN]'
      IF(L.EQ.55)HEAD(M)='NH3_FLUX'
      IF(L.EQ.56)HEAD(M)='UP_NH4_1'
      IF(L.EQ.57)HEAD(M)='UP_NH4_2'
      IF(L.EQ.58)HEAD(M)='UP_NH4_3'
      IF(L.EQ.59)HEAD(M)='UP_NH4_4'
      IF(L.EQ.60)HEAD(M)='UP_NH4_5'
      IF(L.EQ.61)HEAD(M)='UP_NH4_6'
      IF(L.EQ.62)HEAD(M)='UP_NH4_7'
      IF(L.EQ.63)HEAD(M)='UP_NH4_8'
      IF(L.EQ.64)HEAD(M)='UP_NH4_9'
      IF(L.EQ.65)HEAD(M)='UP_NH4_10'
      IF(L.EQ.66)HEAD(M)='UP_NH4_11'
      IF(L.EQ.67)HEAD(M)='UP_NH4_12'
      IF(L.EQ.68)HEAD(M)='UP_NH4_13'
      IF(L.EQ.69)HEAD(M)='UP_NH4_14'
      IF(L.EQ.70)HEAD(M)='UP_NH4_15'
      IF(L.EQ.71)HEAD(M)='UP_NO3_1'
      IF(L.EQ.72)HEAD(M)='UP_NO3_2'
      IF(L.EQ.73)HEAD(M)='UP_NO3_3'
      IF(L.EQ.74)HEAD(M)='UP_NO3_4'
      IF(L.EQ.75)HEAD(M)='UP_NO3_5'
      IF(L.EQ.76)HEAD(M)='UP_NO3_6'
      IF(L.EQ.77)HEAD(M)='UP_NO3_7'
      IF(L.EQ.78)HEAD(M)='UP_NO3_8'
      IF(L.EQ.79)HEAD(M)='UP_NO3_9'
      IF(L.EQ.80)HEAD(M)='UP_NO3_10'
      IF(L.EQ.81)HEAD(M)='UP_NO3_11'
      IF(L.EQ.82)HEAD(M)='UP_NO3_12'
      IF(L.EQ.83)HEAD(M)='UP_NO3_13'
      IF(L.EQ.84)HEAD(M)='UP_NO3_14'
      IF(L.EQ.85)HEAD(M)='UP_NO3_15'
      ENDIF
1023  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.EQ.24)THEN
      DO 1024 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='PO4_UPTK'
      IF(L.EQ.52)HEAD(M)='[TNP]'
      IF(L.EQ.53)HEAD(M)='UP_PO4_1'
      IF(L.EQ.54)HEAD(M)='UP_PO4_2'
      IF(L.EQ.55)HEAD(M)='UP_PO4_3'
      IF(L.EQ.56)HEAD(M)='UP_PO4_4'
      IF(L.EQ.57)HEAD(M)='UP_PO4_5'
      IF(L.EQ.58)HEAD(M)='UP_PO4_6'
      IF(L.EQ.59)HEAD(M)='UP_PO4_7'
      IF(L.EQ.60)HEAD(M)='UP_PO4_8'
      IF(L.EQ.61)HEAD(M)='UP_PO4_9'
      IF(L.EQ.62)HEAD(M)='UP_PO4_10'
      IF(L.EQ.63)HEAD(M)='UP_PO4_11'
      IF(L.EQ.64)HEAD(M)='UP_PO4_12'
      IF(L.EQ.65)HEAD(M)='UP_PO4_13'
      IF(L.EQ.66)HEAD(M)='UP_PO4_14'
      IF(L.EQ.67)HEAD(M)='UP_PO4_15'
      ENDIF
1024  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.EQ.25)THEN
      DO 1025 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='CAN_RN'
      IF(L.EQ.52)HEAD(M)='CAN_LE'
      IF(L.EQ.53)HEAD(M)='CAN_H'
      IF(L.EQ.54)HEAD(M)='CAN_G'
      IF(L.EQ.55)HEAD(M)='CAN_TEMP'
      IF(L.EQ.56)HEAD(M)='TEMP_FN'
      ENDIF
1025  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.EQ.26)THEN
      DO 1026 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='SHOOT_C'
      IF(L.EQ.52)HEAD(M)='LEAF_C'
      IF(L.EQ.53)HEAD(M)='SHTH_C'
      IF(L.EQ.54)HEAD(M)='STALK_C'
      IF(L.EQ.55)HEAD(M)='RESERVE_C'
      IF(L.EQ.56)HEAD(M)='REPRO_C'
      IF(L.EQ.57)HEAD(M)='GRAIN_C'
      IF(L.EQ.58)HEAD(M)='ROOT_C'
      IF(L.EQ.59)HEAD(M)='NDL_C'
      IF(L.EQ.60)HEAD(M)='STORED_C'
      IF(L.EQ.61)HEAD(M)='GRAIN_NO.'
      IF(L.EQ.62)HEAD(M)='LAI'
      IF(L.EQ.63)HEAD(M)='GPP'
      IF(L.EQ.64)HEAD(M)='EXUD_C'
      IF(L.EQ.65)HEAD(M)='LITTER_C'
      IF(L.EQ.66)HEAD(M)='SF_LIT_C'
      IF(L.EQ.67)HEAD(M)='AUTO_RESP'
      IF(L.EQ.68)HEAD(M)='ABV_GRD_RESP'
      IF(L.EQ.69)HEAD(M)='[SOL._C]'
      IF(L.EQ.70)HEAD(M)='HVST_C'
      IF(L.EQ.71)HEAD(M)='DNS_RT_1'
      IF(L.EQ.72)HEAD(M)='DNS_RT_2'
      IF(L.EQ.73)HEAD(M)='DNS_RT_3'
      IF(L.EQ.74)HEAD(M)='DNS_RT_4'
      IF(L.EQ.75)HEAD(M)='DNS_RT_5'
      IF(L.EQ.76)HEAD(M)='DNS_RT_6'
      IF(L.EQ.77)HEAD(M)='DNS_RT_7'
      IF(L.EQ.78)HEAD(M)='DNS_RT_8'
      IF(L.EQ.79)HEAD(M)='DNS_RT_9'
      IF(L.EQ.80)HEAD(M)='DNS_RT_10'
      IF(L.EQ.81)HEAD(M)='DNS_RT_11'
      IF(L.EQ.82)HEAD(M)='DNS_RT_12'
      IF(L.EQ.83)HEAD(M)='DNS_RT_13'
      IF(L.EQ.84)HEAD(M)='DNS_RT_14'
      IF(L.EQ.85)HEAD(M)='DNS_RT_15'
      IF(L.EQ.86)HEAD(M)='C_BALANCE'
      IF(L.EQ.87)HEAD(M)='STG_DEAD_C'
      IF(L.EQ.88)HEAD(M)='FIRE_CO2'
      IF(L.EQ.89)HEAD(M)='FIRE_CH4'
      IF(L.EQ.90)HEAD(M)='NPP'
      IF(L.EQ.91)HEAD(M)='CAN_HT'
      IF(L.EQ.92)HEAD(M)='POPN'
      ENDIF
1026  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.EQ.27)THEN
      DO 1027 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='TRANSPN'
      IF(L.EQ.52)HEAD(M)='WTR_STRESS'
      IF(L.EQ.53)HEAD(M)='OXY_STRESS'
      ENDIF
1027  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.EQ.28)THEN
      DO 1028 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='SHOOT_N'
      IF(L.EQ.52)HEAD(M)='LEAF_N'
      IF(L.EQ.53)HEAD(M)='SHTH_N'
      IF(L.EQ.54)HEAD(M)='STALK_N'
      IF(L.EQ.55)HEAD(M)='RESERVE_N'
      IF(L.EQ.56)HEAD(M)='HUSK_N'
      IF(L.EQ.57)HEAD(M)='GRAIN_N'
      IF(L.EQ.58)HEAD(M)='ROOT_N'
      IF(L.EQ.59)HEAD(M)='NDL_N'
      IF(L.EQ.60)HEAD(M)='STORED_N'
      IF(L.EQ.61)HEAD(M)='UPTK_N'
      IF(L.EQ.62)HEAD(M)='LITTER_N'
      IF(L.EQ.63)HEAD(M)='TL_N_FIXED'
      IF(L.EQ.64)HEAD(M)='[SOL_N]'
      IF(L.EQ.65)HEAD(M)='N_STRESS'
      IF(L.EQ.66)HEAD(M)='[SOL_P]'
      IF(L.EQ.67)HEAD(M)='P_STRESS'
      IF(L.EQ.68)HEAD(M)='NH3_FLUX'
      IF(L.EQ.69)HEAD(M)='HVST_N'
      IF(L.EQ.70)HEAD(M)='N_BALANCE'
      IF(L.EQ.71)HEAD(M)='STG_DEAD_N'
      IF(L.EQ.72)HEAD(M)='FIRE_N'
      IF(L.EQ.73)HEAD(M)='SF_LIT_N'
      ENDIF
1028  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.EQ.29)THEN
      DO 1029 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='SHOOT_P'
      IF(L.EQ.52)HEAD(M)='LEAF_P'
      IF(L.EQ.53)HEAD(M)='SHTH_P'
      IF(L.EQ.54)HEAD(M)='STALK_P'
      IF(L.EQ.55)HEAD(M)='RESERVE_P'
      IF(L.EQ.56)HEAD(M)='HUSK_P'
      IF(L.EQ.57)HEAD(M)='GRAIN_P'
      IF(L.EQ.58)HEAD(M)='ROOT_P'
      IF(L.EQ.59)HEAD(M)='NDL_P'
      IF(L.EQ.60)HEAD(M)='STORED_P'
      IF(L.EQ.61)HEAD(M)='UPTK_P'
      IF(L.EQ.62)HEAD(M)='LITTER_P'
      IF(L.EQ.63)HEAD(M)='[SOL_P]'
      IF(L.EQ.64)HEAD(M)='P_STRESS'
      IF(L.EQ.65)HEAD(M)='HVST_P'
      IF(L.EQ.66)HEAD(M)='P_BALANCE'
      IF(L.EQ.67)HEAD(M)='STG_DEAD_P'
      IF(L.EQ.68)HEAD(M)='FIRE_P'
      IF(L.EQ.69)HEAD(M)='SF_LIT_P'
      ENDIF
1029  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.EQ.30)THEN
      DO 1030 L=51,100
      IF(CHOICE(L,N-20).EQ.'YES')THEN
      M=M+1
      IF(L.EQ.51)HEAD(M)='GROWTH_STG'
      IF(L.EQ.52)HEAD(M)='BRANCH_NO.'
      IF(L.EQ.53)HEAD(M)='NODE_NO.'
      IF(L.EQ.54)HEAD(M)='RUB_ACTVN'
      IF(L.EQ.55)HEAD(M)='LEAF_N/C'
      IF(L.EQ.56)HEAD(M)='LEAF_P/C'
      IF(L.EQ.57)HEAD(M)='MIN_LWP'
      IF(L.EQ.58)HEAD(M)='O2_STRESS'
      IF(L.EQ.59)HEAD(M)='TEMP_STRESS'
      ENDIF
1030  CONTINUE
      NOUTP(N-20)=M
      ENDIF
      IF(N.LE.25)WRITE(LUN,'(A12,2A8,50A16)')CDOY,DATE,HOUR
     2,(HEAD(K),K=1,M)
      IF(N.GE.26)WRITE(LUN,'(A12,A8,50A16)')CDOY,DATE
     2,(HEAD(K),K=1,M)
      ENDIF
1010  CONTINUE
      RETURN
      END
