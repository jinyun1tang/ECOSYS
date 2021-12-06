      SUBROUTINE exec(I)
C
C     THIS SUBROUTINE TAKES MASS BALANCE VARIABLES CALCULATED
C     IN 'REDIST' AND PERFORMS MASS BALANCE CHECKS AT THE END
C     OF EACH DAY OF THE MODEL RUN, AND ERRORS OF > 1UG ARE FLAGGED.
C
      include "parameters.h"
      include "blkc.h"
      include "blk2a.h"
      include "blk2b.h"
      include "blk2c.h"
      include "blk16.h"
      SAVE TLW,TLH,TLO,TLC,TLN,TLP,TLI
C
C     CALCULATE MASS BALANCES FOR WATER, HEAT, O2, C, N, P AND SOLUTES
C
      IF(I.EQ.IBEGIN.OR.I.EQ.ISTART.OR.I.EQ.ILAST+1)THEN
      TLW=VOLWSO-CRAIN+CRUN+CEVAP+VOLWOU
      TLH=HEATSO-HEATIN+HEATOU
      TLO=OXYGSO-OXYGIN+OXYGOU
      TLC=TLRSDC+TLORGC+TLCO2G-CO2GIN+TCOU-TORGF-XCSN
      TLN=TLRSDN+TLORGN+TLN2G+TLNH4+TLNO3-ZN2GIN-TZIN+TZOU-TORGN-XZSN
      TLP=TLRSDP+TLORGP+TLPO4-TPIN+TPOU-TORGP-XPSN
      TLI=TION-TIONIN+TIONOU
      ENDIF
C
C     CALCULATE DEVIATION SINCE MASS BALANCE WAS LAST RESET
C
      IF((I/IOUT)*IOUT.EQ.I)THEN
      DIFFQ=(VOLWSO-CRAIN+CRUN+CEVAP+VOLWOU-TLW)/TAREA
      DIFFH=(HEATSO-HEATIN+HEATOU-TLH)/TAREA
      DIFFO=(OXYGSO-OXYGIN+OXYGOU-TLO)/TAREA
      DIFFC=(TLRSDC+TLORGC+TLCO2G-CO2GIN+TCOU-TORGF-XCSN-TLC)/TAREA
      DIFFN=(TLRSDN+TLORGN+TLN2G+TLNH4+TLNO3-ZN2GIN-TZIN+TZOU
     2-TORGN-XZSN-TLN)/TAREA
      DIFFP=(TLRSDP+TLORGP+TLPO4-TPIN+TPOU-TORGP-XPSN-TLP)/TAREA
      DIFFI=(TION-TIONIN+TIONOU-TLI)/TAREA
      WRITE(*,212)I,IYRC
      WRITE(18,213)I,IYRC,DIFFQ,DIFFH,DIFFO,DIFFC,DIFFN
     2,DIFFP,DIFFI
212   FORMAT('NOW EXECUTING DAY',I6,'   OF YEAR',I6)
213   FORMAT(2I6,10F16.6)
C
C     FLAG DEVIATIONS > 1UG OR 1 J IN ENTIRE MODEL LANDSCAPE,
C     RESET MASS BALANCE
C
      IF(ABS(DIFFC).GT.1.0E-06)THEN
      WRITE(18,191)I,IYRC
191   FORMAT('CARBON BALANCE LOST ON DAY, YEAR',2I4)
      TLC=TLRSDC+TLORGC+TLCO2G-CO2GIN+TCOU-TORGF-XCSN
      ENDIF
      IF(ABS(DIFFN).GT.1.0E-06)THEN
      WRITE(18,192)I,IYRC
192   FORMAT('NITROGEN BALANCE LOST ON DAY, YEAR',2I4)
      TLN=TLRSDN+TLORGN+TLN2G+TLNH4+TLNO3-ZN2GIN-TZIN+TZOU-TORGN-XZSN
      ENDIF
      IF(ABS(DIFFP).GT.1.0E-06)THEN
      WRITE(18,193)I,IYRC
193   FORMAT('PHOSPHORUS BALANCE LOST ON DAY, YEAR',2I4)
      TLP=TLRSDP+TLORGP+TLPO4-TPIN+TPOU-TORGP-XPSN
      ENDIF
      IF(ABS(DIFFQ).GT.1.0E-06)THEN
      WRITE(18,194)I,IYRC
194   FORMAT('WATER BALANCE LOST ON DAY, YEAR',2I4)
      TLW=VOLWSO-CRAIN+CRUN+CEVAP+VOLWOU
      ENDIF
      IF(ABS(DIFFH).GT.1.0E-06)THEN
      WRITE(18,195)I,IYRC
195   FORMAT('THERMAL BALANCE LOST ON DAY, YEAR',2I4)
      TLH=HEATSO-HEATIN+HEATOU
      ENDIF
      IF(ABS(DIFFO).GT.1.0E-06)THEN
      WRITE(18,196)I,IYRC
196   FORMAT('OXYGEN BALANCE LOST ON DAY, YEAR',2I4)
      TLO=OXYGSO-OXYGIN+OXYGOU
      ENDIF
      IF(ABS(DIFFI).GT.1.0E-06)THEN
      WRITE(18,197)I,IYRC
197   FORMAT('ION BALANCE LOST ON DAY, YEAR',2I4)
      TLI=TION-TIONIN+TIONOU
      ENDIF
      ENDIF
      IF(IDAYR.LT.0)THEN
      IDAYR=LYRX+IDAYR
      ELSE
      IDAYR=I
      ENDIF
      IOLD=I
      IMNG=0
      NYR=0
      RETURN
      END