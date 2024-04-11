      SUBROUTINE erosion(I,J,NFZ,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE CALCULATES DETACHMENT AND OVERLAND TRANSPORT
C     OF SURFACE SEDIMENT FROM PRECIPITATION IN WEATHER FILE AND
C     FROM RUNOFF IN 'WATSUB'
C
      include "parameters.h"
      include "blkc.h"
      include "blk2a.h"
      include "blk5.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk10.h"
      include "blk11a.h"
      include "blk13a.h"
      include "blk13b.h"
      include "blk13c.h"
      include "blk19a.h"
      include "blk19b.h"
      include "blk19c.h"
      include "blk20f.h"
      DIMENSION RERSED(2,2,JV,JH),TERSED(JY,JX),RDTSED(JY,JX)
     2,FVOLIM(JY,JX),FVOLWM(JY,JX),FERSNM(JY,JX),RERSED0(JY,JX)
C
C     INTERNAL TIME STEP AT WHICH SEDIMENT DETACHMENT AND TRANSPORT
C     IS CALCULATED. DETACHMENT IS THE SUM OF THAT BY RAINFALL AND
C     OVERLAND FLOW
C
      DO 30 M=1,NPH
      DO 9895 NX=NHW,NHE
      DO 9890 NY=NVN,NVS
C
C     IERSNG=options for disturbance effects on soil profile layer
C        depths and contents from site file:
C           :-1=no effects
C           :0=freeze-thaw
C           :1=freeze-thaw+erosion
C           :2=freeze-thaw+SOM gain or loss
C           :3=freeze-thaw+erosion+SOM gain or loss 
C     TERSED=net sediment erosion (Mg t-1)
C     RDTSED=sediment detachment rate (Mg t-1)
C     FVOLWM,FVOLIM=fraction of surface ponding capacity 
C        occupied by water,ice
C     FERSM=fraction of surface ponding capacity 
C        not occupied by water,ice 
C
      IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
      TERSED(NY,NX)=0.0
      RDTSED(NY,NX)=0.0
      FVOLWM(NY,NX)=AMIN1(1.0,AMAX1(0.0,XVOLWM(M,NY,NX)/VOLWG(NY,NX)))
      FVOLIM(NY,NX)=AMIN1(1.0,AMAX1(0.0,XVOLIM(M,NY,NX)/VOLWG(NY,NX)))
      FERSNM(NY,NX)=(1.0-FVOLIM(NY,NX))*FVOLWM(NY,NX) 
C
C     DETACHMENT BY RAINFALL WHEN SURFACE WATER IS PRESENT
C
C     BKDS=surface soil bulk density (0=water,>0=soil) (Mg m-3)
C     ENGYPM=total rainfall energy impact from �watsub.f� (J t-1)
C     XVOLWM=surface water in excess of litter water 
C        retention capacity (m3) 
C
      IF(BKDS(NU(NY,NX),NY,NX).GT.ZERO
     2.AND.ENGYPM(M,NY,NX).GT.0.0
     3.AND.XVOLWM(M,NY,NX).GT.ZEROS(NY,NX))THEN 
C
C     DETACHMENT OF SEDIMENT FROM SURFACE SOIL DEPENDS ON RAINFALL
C     KINETIC ENERGY AND FROM DETACHMENT COEFFICIENT IN 'HOUR1'
C     ATTENUATED BY DEPTH OF SURFACE WATER
C
C     DETS=soil detachability from rainfall impact from �hour1.f�
C        (g J-1)  
C     DETW=DETS accounting for soil surface water content (g J-1)
C     VOLWM,VOLA=soil surface water content,porous volume (m3)
C     DETR=sediment detachment from rainfall impact (Mg t-1)
C     BKVL=soil mass (Mg) 
C     ENGYPM=total rainfall energy impact from �watsub.f� (J t-1)
C     AREA=grid cell surface area (m2)
C     FMPR=1.0-(coarse fragment+macropore) fraction
C     FSNX=fraction of snow-free cover
C     FVOLIM=fraction of surface ponding capacity occupied by ice 
C     RDTSED=sediment detachment rate (Mg t-1)
C
      DETW=DETS(NY,NX)*(1.0+2.0*VOLWM(M,NU(NY,NX),NY,NX)
     2/VOLA(NU(NY,NX),NY,NX))
      DETR=AMIN1(BKVL(NU(NY,NX),NY,NX)*XNPXX 
     2,DETW*ENGYPM(M,NY,NX)*AREA(3,NU(NY,NX),NY,NX)
     3*FMPR(NU(NY,NX),NY,NX)*FSNX(NY,NX)*(1.0-FVOLIM(NY,NX)))
      RDTSED(NY,NX)=RDTSED(NY,NX)+DETR
C     WRITE(*,1117)'DETR',I,J,NFZ,M,NX,NY,RDTSED(NY,NX)
C    2,SED(NY,NX),DETR
C    2,PRECA(NY,NX)*1000.0,PRECD(NY,NX),PRECB(NY,NX),DETS(NY,NX) 
C    3,FSNX(NY,NX),FVOLIM(NY,NX),XVOLWM(M,NY,NX),ENGYPM(M,NY,NX)
C    4,DETW,VOLWM(M,NU(NY,NX),NY,NX),VOLA(NU(NY,NX),NY,NX)
1117  FORMAT(A8,6I4,20E12.4)
      ENDIF
C
C     DEPOSITION OF SEDIMENT TO SOIL SURFACE FROM IMMOBILE SURFACE
C     WATER
C
C     SED=sediment content of surface water (Mg)
C     RDTSED=sediment detachment rate (Mg t-1)
C     XVOLTM=surface water+ice in excess of litter water 
C        retention capacity (m3) 
C     FERSM=fraction of surface ponding capacity 
C        not occupied by water,ice 
C     FVOLWM=fraction of surface ponding capacity occupied by water
C     CSEDD=sediment concentration in surface water (Mg m-3)
C     DEPI=sediment deposition rate (Mg t-1)
C     VLS=sediment sinking rate from �hour1.f� (m h-1) 
C     AREA=grid cell surface area (m2)
C     FMPR=1.0-(coarse fragment+macropore) fraction
C     XNPHX=time step for water fluxes from �wthr.f� (h t-1)
C
      SEDX=SED(NY,NX)+RDTSED(NY,NX)
      IF(BKDS(NU(NY,NX),NY,NX).GT.ZERO
     2.AND.SEDX.GT.ZEROS(NY,NX)
     3.AND.XVOLTM(M,NY,NX).LE.VOLWG(NY,NX)
     4.AND.FERSNM(NY,NX).GT.ZERO)THEN
      CSEDD=AMAX1(0.0,SEDX/XVOLWM(M,NY,NX))
      DEPI=AMAX1(-SEDX 
     2,VLS(NY,NX)*(0.0-CSEDD)*AREA(3,NU(NY,NX),NY,NX) 
     3*FERSNM(NY,NX)*FMPR(NU(NY,NX),NY,NX)*XNPHX)
      RDTSED(NY,NX)=RDTSED(NY,NX)+DEPI
C     WRITE(*,1117)'DEPI',I,J,NFZ,M,NX,NY,RDTSED(NY,NX),SEDX,DEPI
C    2,CSEDD,FVOLIM(NY,NX),FVOLWM(NY,NX),FERSNM(NY,NX)
C    3,VLS(NY,NX),XVOLWM(M,NY,NX)
      ENDIF
C
C     DETACHMENT IN SURFACE WATER FROM OVERLAND WATER
C     VELOCITY FROM 'WATSUB' USED TO CALCULATE STREAM POWER,
C     AND FROM SEDIMENT TRANSPORT CAPACITY VS. CURRENT SEDIMENT
C     CONCENTRATION IN SURFACE WATER, MODIFIED BY SOIL COHESION
C     FROM 'HOUR1'
C
C     BKDS=surface soil bulk density (0=water,>0=soil) (Mg m-3)
C     XVOLTM=surface water+ice in excess of litter water 
C        retention capacity (m3)
C     FERSM=fraction of surface ponding capacity 
C        not occupied by water,ice 
C     STPR=stream power of runoff (cm s-1)
C     QRV=runoff velocity from �watsub.f� (m s-1)
C     SLOPE=sin(slope) from site file
C     PTDSNU=particle density in soil surface layer from �hour1.f� 
C        (Mg m-3)
C     CER,XER=parameters for runoff transport capacity from �hour1.f�
C     CSEDX=sediment concentration holding capacity of runoff (Mg m-3) 
C     CSEDD=sediment concentration in surface water (Mg m-3)
C     DETI=sediment detachment(+ve) or deposition(-ve) rate 
C        from runoff (Mg t-1)
C     BKVL=soil mass (Mg)
C     DETE=soil detachability coefficient from �hour1.f�    
C     AREA=grid cell surface area (m2)
C     FMPR=1.0-(coarse fragment+macropore) fraction
C     VLS=sediment sinking rate (m h-1) 
C     RDTSED=sediment detachment rate (Mg t-1)
C     XNPHX=time step for water fluxes from �wthr.f� (h t-1)
C
      IF(BKDS(NU(NY,NX),NY,NX).GT.ZERO
     3.AND.XVOLTM(M,NY,NX).GT.VOLWG(NY,NX)
     3.AND.FERSNM(NY,NX).GT.ZERO)THEN
      STPR=1.0E+02*QRV(M,NY,NX)*ABS(SLOPE(0,NY,NX))
      CSEDX=PTDSNU(NY,NX)*CER(NY,NX)*AMAX1(0.0,STPR-0.4)**XER(NY,NX)
      CSEDD=AMAX1(0.0,SEDX/XVOLWM(M,NY,NX))
      IF(CSEDX.GT.CSEDD)THEN
      DETI=AMIN1(BKVL(NU(NY,NX),NY,NX)*XNPXX 
     2,DETE(NY,NX)*VLS(NY,NX)*(CSEDX-CSEDD)*AREA(3,NU(NY,NX),NY,NX) 
     3*FERSNM(NY,NX)*FMPR(NU(NY,NX),NY,NX)*XNPHX)
      ELSE
      IF(SEDX.GT.ZEROS(NY,NX))THEN
      DETI=AMAX1(-SEDX 
     2,VLS(NY,NX)*(CSEDX-CSEDD)*AREA(3,NU(NY,NX),NY,NX) 
     3*FERSNM(NY,NX)*FMPR(NU(NY,NX),NY,NX)*XNPHX)
      ELSE
      DETI=0.0
      ENDIF
      ENDIF
      RDTSED(NY,NX)=RDTSED(NY,NX)+DETI
C     WRITE(*,1112)'DETI',I,J,NFZ,M,NX,NY,RDTSED(NY,NX)
C    2,SED(NY,NX),DETI
C    2,QRM(M,NY,NX),QRV(M,NY,NX),SEDX,XVOLWM(M,NY,NX),XVOLTM(M,NY,NX)  
C    4,VOLWG(NY,NX),STPR,CSEDX,CSEDD,VLS(NY,NX),PTDSNU(NY,NX)
C    5,DETE(NY,NX),SLOPE(0,NY,NX),ZM(NY,NX),FERSNM(NY,NX)
1112  FORMAT(A8,6I4,30E12.4)
      ENDIF
C
C     TRANSPORT OF SEDIMENT IN OVERLAND FLOW FROM SEDIMENT
C     CONCENTRATION TIMES OVERLAND WATER FLUX FROM 'WATSUB'
C
C     N2,N1=NY,NX of source grid cell
C     QRM=downslope runoff rate from �watsub.f� (m3 t-1)
C     BKDS=surface soil bulk density (0=water,>0=soil) (Mg m-3)
C     RERSED0=sediment transport down hillslope (Mg t-1)
C     XVOLWM=surface water in excess of litter water 
C        retention capacity (m3) 
C     SED=sediment content of surface water (Mg)
C     RDTSED=sediment detachment rate (Mg t-1)
C     CSEDE=sediment concentration in surface water (Mg m-3)
C     FVOLIM=fraction of surface ponding capacity occupied by ice
C
      N1=NX
      N2=NY
      IF(QRM(M,N2,N1).LE.0.0
     2.OR.BKDS(NU(N2,N1),N2,N1).LE.ZERO)THEN
      RERSED0(N2,N1)=0.0
      ELSE
      IF(XVOLWM(M,N2,N1).GT.ZEROS2(N2,N1))THEN
      SEDX=SED(N2,N1)+RDTSED(N2,N1)
      CSEDE=AMAX1(0.0,SEDX/XVOLWM(M,N2,N1))
      RERSED0(N2,N1)=AMIN1(SEDX,CSEDE*QRM(M,N2,N1)
     2*(1.0-FVOLIM(N2,N1)))
      ELSE
      RERSED0(N2,N1)=0.0
      ENDIF
      ENDIF
C     IF(RERSED0(N2,N1).GT.ZEROS(N2,N1))THEN
C     WRITE(*,1121)'RERSED0',I,J,NFZ,M,N1,N2
C    2,RERSED0(N2,N1),QRM(M,N2,N1),XVOLWM(M,N2,N1)
C    3,SED(N2,N1),RDTSED(N2,N1),CSEDE,FVOLIM(N2,N1)
1121  FORMAT(A8,6I4,12E12.4)
C     ENDIF
C
C     LOCATE INTERNAL BOUNDARIES
C
C     RERSED0=sediment transport down hillslope (Mg t-1)
C     N2,N1=NY,NX of source grid cell
C     N5,N4=NY,NX of destination grid cell E or S
C     N5B,N4B=NY,NX of destination grid cell W or N
C     NN=boundary:N=1:NN=1 east,NN=2 west, N=2:NN=1 south,NN=2 north
C
      IF(RERSED0(N2,N1).GT.0.0)THEN
      DO 4310 N=1,2
      DO 4305 NN=1,2
      IF(N.EQ.1)THEN
      IF(NX.EQ.NHE.AND.NN.EQ.1
     2.OR.NX.EQ.NHW.AND.NN.EQ.2)THEN
      GO TO 4305
      ELSE
      N4=NX+1
      N5=NY
      N4B=NX-1
      N5B=NY
      ENDIF
      ELSEIF(N.EQ.2)THEN
      IF(NY.EQ.NVS.AND.NN.EQ.1
     2.OR.NY.EQ.NVN.AND.NN.EQ.2)THEN
      GO TO 4305
      ELSE
      N4=NX
      N5=NY+1
      N4B=NX
      N5B=NY-1
      ENDIF
      ENDIF
C
C     PARTITION DOWNSLOPE EROSION INTO WE AND NS DIRECTIONS
C
C     QRM=runoff from �watsub.f� (m3 t-1)
C     QRMN=runoff in EW(N=1),NS(N=2) directions from �watsub.f� 
C        (m3 t-1)
C     FERM=fraction of QRM partitioned to WE(N=1),NS(N=2) directions
C     RERSED=sediment transport in WE,NS directions (Mg t-1)
C     XSEDER=cumulative sediment transport in WE,NS directions (Mg t-1)
C
      IF(RERSED0(N2,N1).GT.ZEROS(N2,N1))THEN
      IF(NN.EQ.1)THEN
      FERM=QRMN(M,N,2,N5,N4)/QRM(M,N2,N1)
      RERSED(N,2,N5,N4)=RERSED0(N2,N1)*FERM  
      XSEDER(N,2,N5,N4)=XSEDER(N,2,N5,N4)+RERSED(N,2,N5,N4)
C     IF(N1.EQ.1)THEN
C     WRITE(*,1113)'INTF',I,J,NFZ,M,N1,N2,N4,N5,N
C    2,RERSED0(N2,N1),RERSED(N,2,N5,N4),XSEDER(N,2,N5,N4) 
C    3,SED(N2,N1),XVOLWM(M,N2,N1) 
1113  FORMAT(A8,9I4,30E12.4)
C     ENDIF
      ELSE
      RERSED(N,2,N5,N4)=0.0  
      ENDIF
      IF(NN.EQ.2)THEN
      IF(N4B.GT.0.AND.N5B.GT.0)THEN
      FERM=QRMN(M,N,1,N5B,N4B)/QRM(M,N2,N1)
      RERSED(N,1,N5B,N4B)=RERSED0(N2,N1)*FERM 
      XSEDER(N,1,N5B,N4B)=XSEDER(N,1,N5B,N4B)+RERSED(N,1,N5B,N4B)
C     IF(N1.EQ.1)THEN
C     WRITE(*,1113)'INTB',I,J,NFZ,M,N1,N2,N4B,N5B,N
C    2,RERSED0(N2,N1),RERSED(N,1,N5B,N4B),XSEDER(N,1,N5B,N4B) 
C    3,SED(N2,N1),XVOLWM(M,N2,N1) 
C     ENDIF
      ELSE
      RERSED(N,1,N5B,N4B)=0.0 
      ENDIF
      ENDIF
      ELSE
      RERSED(N,2,N5,N4)=0.0
      IF(N4B.GT.0.AND.N5B.GT.0)THEN
      RERSED(N,1,N5B,N4B)=0.0 
      ENDIF
      ENDIF 
4305  CONTINUE
4310  CONTINUE
      ENDIF
      ENDIF
9890  CONTINUE
9895  CONTINUE
C
C     BOUNDARY SEDIMENT FLUXES
C
C     N2,N1=NY,NX of source grid cell
C     IERSNG=options for disturbance effects on soil profile layer
C        depths and contents:
C           :-1=no effects
C           :0=freeze-thaw
C           :1=freeze-thaw+erosion
C           :2=freeze-thaw+SOM gain or loss
C           :3=freeze-thaw+erosion+SOM gain or loss 
C     QRM=downslope runoff rate from �watsub.f� (m3 t-1)
C     BKDS=surface soil bulk density (0=water,>0=soil) (Mg m-3)
C     RERSED0=sediment transport down hillslope (Mg t-1)
C     XVOLWM=surface water in excess of litter water 
C        retention capacity (m3) 
C     SED=sediment content of surface water (Mg)
C     RDTSED=sediment detachment rate (Mg t-1)
C     CSEDE=sediment concentration in surface water (Mg m-3)
C
      DO 9595 NX=NHW,NHE
      DO 9590 NY=NVN,NVS
      IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
      N1=NX
      N2=NY
      IF(QRM(M,N2,N1).LE.0.0
     2.OR.BKDS(NU(N2,N1),N2,N1).LE.ZERO)THEN
      RERSED0(N2,N1)=0.0
      ELSE
      IF(XVOLWM(M,N2,N1).GT.ZEROS2(N2,N1))THEN
      SEDX=SED(NY,NX)+RDTSED(NY,NX)
      CSEDE=AMAX1(0.0,SEDX/XVOLWM(M,N2,N1))
      RERSED0(N2,N1)=AMIN1(SEDX,CSEDE*QRM(M,N2,N1))
      ELSE
      RERSED0(N2,N1)=0.0
      ENDIF
      ENDIF
C
C     LOCATE EXTERNAL BOUNDARIES
C
C     M5,M4=NY,NX of destination grid cell
C     N5,N4=NY,NX of destination grid cell E or S
C     N5B,N4B=NY,NX of destination grid cell W or N
C     NN=boundary:N=1:NN=1 east,NN=2 west, N=2:NN=1 south,NN=2 north
C     RCHQE,RCHQW,RCHQN,RCHQS=flux conditions at E,W,N,S boundaries
C        from site file
C
      DO 9580 N=1,2
      DO 9575 NN=1,2
      IF(N.EQ.1)THEN
      N4=NX+1
      N5=NY
      N4B=NX-1
      N5B=NY
      IF(NN.EQ.1)THEN
      IF(NX.EQ.NHE)THEN
      M1=NX
      M2=NY
      M4=NX+1
      M5=NY
      XN=-1.0
      RCHQF=RCHQE(M2,M1)
      ELSE
      GO TO 9575
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NX.EQ.NHW)THEN
      M1=NX
      M2=NY
      M4=NX
      M5=NY
      XN=1.0
      RCHQF=RCHQW(M5,M4)
      ELSE
      GO TO 9575
      ENDIF
      ENDIF
      ELSEIF(N.EQ.2)THEN
      N4=NX
      N5=NY+1
      N4B=NX
      N5B=NY-1
      IF(NN.EQ.1)THEN
      IF(NY.EQ.NVS)THEN
      M1=NX
      M2=NY
      M4=NX
      M5=NY+1
      XN=-1.0
      RCHQF=RCHQS(M2,M1)
      ELSE
      GO TO 9575
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NY.EQ.NVN)THEN
      M1=NX
      M2=NY
      M4=NX
      M5=NY
      XN=1.0
      RCHQF=RCHQN(M5,M4)
      ELSE
      GO TO 9575
      ENDIF
      ENDIF
      ENDIF
C
C     SEDIMENT TRANSPORT ACROSS BOUNDARY FROM BOUNDARY RUNOFF
C     IN 'WATSUB' TIMES BOUNDARY SEDIMENT CONCENTRATION IN
C     SURFACE WATER
C
C     IRCHG=topographic constraints on runoff from �starts.f�
C     NN=boundary:N=1:NN=1 east,NN=2 west, N=2:NN=1 south,NN=2 north
C     RERSED0=sediment transport down hillslope (Mg t-1)
C     RERSED=sediment transport in WE,NS directions (Mg t-1)
C     QRM=runoff from �watsub.f� (m3 t-1)
C     QRMN=runoff in EW(N=1),NS(N=2) directions from �watsub.f� 
C        (m3 t-1)
C     FERM=fraction of QRM partitioned to WE,NS directions
C     RERSED=sediment transport in WE,NS directions (Mg t-1)
C     XSEDER=hourly sediment transport in WE,NS directions (Mg t-1)
C
      IF(IRCHG(NN,N,N2,N1).EQ.0.OR.RCHQF.EQ.0.0
     2.OR.RERSED0(N2,N1).LE.ZEROS(N2,N1))THEN
      RERSED(N,NN,M5,M4)=0.0
      ELSE
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
      IF((NN.EQ.1.AND.QRMN(M,N,NN,M5,M4).GT.ZEROS(N2,N1))
     2.OR.(NN.EQ.2.AND.QRMN(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      FERM=QRMN(M,N,NN,M5,M4)/QRM(M,N2,N1)
      RERSED(N,NN,M5,M4)=RERSED0(N2,N1)*FERM
      XSEDER(N,NN,M5,M4)=XSEDER(N,NN,M5,M4)+RERSED(N,NN,M5,M4)
      ELSEIF((NN.EQ.2.AND.QRMN(M,N,NN,M5,M4).GT.ZEROS(N2,N1))
     2.OR.(NN.EQ.1.AND.QRMN(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
      RERSED(N,NN,M5,M4)=0.0
      ELSE
      RERSED(N,NN,M5,M4)=0.0
      ENDIF
      ENDIF
      ENDIF
C     IF(RERSED0(N2,N1).GT.ZEROS(N2,N1))THEN
C     WRITE(*,1114)'BNDY',I,J,NFZ,M,N1,N2,M4,M5,N,NN,IRCHG(NN,N,N2,N1)
C    2,RCHQF,RERSED0(N2,N1),RERSED(N,NN,M5,M4),XSEDER(N,NN,M5,M4) 
C    3,SED(M2,M1),QRM(M,N2,N1),QRMN(M,N,NN,M5,M4),XVOLWM(M,M2,M1)
1114  FORMAT(A8,11I4,30E12.4)
C     ENDIF
9575  CONTINUE
C
C     TOTAL SEDIMENT FLUXES
C
C     TERSED=net sediment transport (Mg t-1)
C     RERSED=sediment transport in WE,NS directions (Mg t-1)
C
      DO 1202 NN=1,2
      TERSED(N2,N1)=TERSED(N2,N1)+RERSED(N,NN,N2,N1)
      IF(IFLBM(M,N,NN,N5,N4).EQ.0)THEN
      TERSED(N2,N1)=TERSED(N2,N1)-RERSED(N,NN,N5,N4)
      ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      TERSED(N2,N1)=TERSED(N2,N1)-RERSED(N,NN,N5B,N4B)
      ENDIF
1202  CONTINUE
9580  CONTINUE
      ENDIF
9590  CONTINUE
9595  CONTINUE
C
C     UPDATE STATE VARIABLES FOR SEDIMENT TRANSPORT
C
C     IERSNG=options for disturbance effects on soil profile layer
C        depths and contents:
C           :-1=no effects
C           :0=freeze-thaw
C           :1=freeze-thaw+erosion
C           :2=freeze-thaw+SOM gain or loss
C           :3=freeze-thaw+erosion+SOM gain or loss 
C     SED=sediment content of surface water (Mg)
C     TERSED=net sediment transport (Mg t-1)
C     RDTSED=sediment detachment rate (Mg t-1)
C
      DO 9695 NX=NHW,NHE
      DO 9690 NY=NVN,NVS
      IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
      SED(NY,NX)=SED(NY,NX)+TERSED(NY,NX)+RDTSED(NY,NX)
C     WRITE(*,1118)'SED2',I,J,NFZ,M,NX,NY,SED(NY,NX)
C    2,TERSED(NY,NX),RDTSED(NY,NX) 
1118  FORMAT(A8,6I4,12E12.4)
      ENDIF
9690  CONTINUE
9695  CONTINUE
30    CONTINUE
C
C     INTERNAL SEDIMENT FLUXES
C
C     IERSNG=options for disturbance effects on soil profile layer
C        depths and contents:
C           :-1=no effects
C           :0=freeze-thaw
C           :1=freeze-thaw+erosion
C           :2=freeze-thaw+SOM gain or loss
C           :3=freeze-thaw+erosion+SOM gain or loss 
C     N2,N1=NY,NX of source grid cell
C     N5,N4=NY,NX of destination grid cell E or S
C     N5B,N4B=NY,NX of destination grid cell W or N
C     NN=boundary:N=1:NN=1 east,NN=2 west, N=2:NN=1 south,NN=2 north
C
      DO 9495 NX=NHW,NHE
      DO 9490 NY=NVN,NVS
      IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
      N1=NX
      N2=NY
      DO 9485 N=1,2
      DO 9450 NN=1,2
      IF(N.EQ.1)THEN
      IF(NX.EQ.NHE.AND.NN.EQ.1
     2.OR.NX.EQ.NHW.AND.NN.EQ.2)THEN
      GO TO 9450
      ELSE
      N4=NX+1
      N5=NY
      N4B=NX-1
      N5B=NY
      ENDIF
      ELSEIF(N.EQ.2)THEN
      IF(NY.EQ.NVS.AND.NN.EQ.1
     2.OR.NY.EQ.NVN.AND.NN.EQ.2)THEN
      GO TO 9450
      ELSE
      N4=NX
      N5=NY+1
      N4B=NX
      N5B=NY-1
      ENDIF
      ENDIF
C
C     FLUXES OF ALL SOLID MATERIALS IN SEDIMENT ARE CALCULATED
C     FROM VALUES OF THEIR CURRENT STATE VARIABLES MULTIPLIED
C     BY THE FRACTION OF THE TOTAL SURFACE LAYER MASS THAT IS
C     TRANSPORTED IN SEDIMENT
C
C     XSEDER=cumulative sediment transport in WE,NS directions (Mg t-1)
C     BKVLNU=mass of surface soil layer from �hour1.f� (Mg)
C     FSEDER=fraction of soil surface layer mass transported 
C        by erosion (t-1)
C     X*ER,X*EB= sediment flux in non-band,band used in �redist.f�
C        (Mg,g or mol t-1) 
C
      IF(NN.EQ.1)THEN
      FSEDER=AMIN1(1.0,XSEDER(N,2,N5,N4)/BKVLNU(N2,N1))
C
C     SOIL MINERAL EROSION
C
C     sediment code:SAN=sand,SIL=silt,CLA=clay (Mg t-1)
C                  :CEC=cation exchange capacity (mol t-1)
C                  :AEC=anion exchange capacity (mol t-1)
C
      XSANER(N,2,N5,N4)=FSEDER*SAND(NU(N2,N1),N2,N1)
      XSILER(N,2,N5,N4)=FSEDER*SILT(NU(N2,N1),N2,N1)
      XCLAER(N,2,N5,N4)=FSEDER*CLAY(NU(N2,N1),N2,N1)
      XCECER(N,2,N5,N4)=FSEDER*XCEC(NU(N2,N1),N2,N1)
      XAECER(N,2,N5,N4)=FSEDER*XAEC(NU(N2,N1),N2,N1)
C
C     EROSION FROM FERTILIZER POOLS (mol t-1)
C
C     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3 
C        in non-band *R and band *B
C
      XNH4ER(N,2,N5,N4)=FSEDER*ZNH4FA(NU(N2,N1),N2,N1)
      XNH3ER(N,2,N5,N4)=FSEDER*ZNH3FA(NU(N2,N1),N2,N1)
      XNHUER(N,2,N5,N4)=FSEDER*ZNHUFA(NU(N2,N1),N2,N1)
      XNO3ER(N,2,N5,N4)=FSEDER*ZNO3FA(NU(N2,N1),N2,N1)
      XNH4EB(N,2,N5,N4)=FSEDER*ZNH4FB(NU(N2,N1),N2,N1)
      XNH3EB(N,2,N5,N4)=FSEDER*ZNH3FB(NU(N2,N1),N2,N1)
      XNHUEB(N,2,N5,N4)=FSEDER*ZNHUFB(NU(N2,N1),N2,N1)
      XNO3EB(N,2,N5,N4)=FSEDER*ZNO3FB(NU(N2,N1),N2,N1)
C
C     EXCHANGEABLE CATION AND ANION EROSION (mol t-1)
C
C     sediment code
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC
C           =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3 
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
C       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C
      XN4ER(N,2,N5,N4)=FSEDER*XN4(NU(N2,N1),N2,N1)
      XNBER(N,2,N5,N4)=FSEDER*XNB(NU(N2,N1),N2,N1)
      XHYER(N,2,N5,N4)=FSEDER*XHY(NU(N2,N1),N2,N1)
      XALER(N,2,N5,N4)=FSEDER*XAL(NU(N2,N1),N2,N1)
      XFEER(N,2,N5,N4)=FSEDER*XFE(NU(N2,N1),N2,N1)
      XCAER(N,2,N5,N4)=FSEDER*XCA(NU(N2,N1),N2,N1)
      XMGER(N,2,N5,N4)=FSEDER*XMG(NU(N2,N1),N2,N1)
      XNAER(N,2,N5,N4)=FSEDER*XNA(NU(N2,N1),N2,N1)
      XKAER(N,2,N5,N4)=FSEDER*XKA(NU(N2,N1),N2,N1)
      XHCER(N,2,N5,N4)=FSEDER*XHC(NU(N2,N1),N2,N1)
      XOH0ER(N,2,N5,N4)=FSEDER*XOH0(NU(N2,N1),N2,N1)
      XOH1ER(N,2,N5,N4)=FSEDER*XOH1(NU(N2,N1),N2,N1)
      XOH2ER(N,2,N5,N4)=FSEDER*XOH2(NU(N2,N1),N2,N1)
      XH1PER(N,2,N5,N4)=FSEDER*XH1P(NU(N2,N1),N2,N1)
      XH2PER(N,2,N5,N4)=FSEDER*XH2P(NU(N2,N1),N2,N1)
      XOH0EB(N,2,N5,N4)=FSEDER*XOH0B(NU(N2,N1),N2,N1)
      XOH1EB(N,2,N5,N4)=FSEDER*XOH1B(NU(N2,N1),N2,N1)
      XOH2EB(N,2,N5,N4)=FSEDER*XOH2B(NU(N2,N1),N2,N1)
      XH1PEB(N,2,N5,N4)=FSEDER*XH1PB(NU(N2,N1),N2,N1)
      XH2PEB(N,2,N5,N4)=FSEDER*XH2PB(NU(N2,N1),N2,N1)
C
C     EROSION OF PRECIPITATES (mol t-1)
C
C     sediment code
C       :PALO,PFEO=precipitated AlOH,FeOH 
C       :PCAC,PCAS=precipitated CaCO3,CaSO4
C       :PALP,PFEP=precipitated AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precipitated AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precipitated CaH2PO4,CaHPO4,apatite
C          in non-band
C       :PCPMB,PCPDB,PCPHB= precipitated CaH2PO4,CaHPO4,apatite 
C          in band
C
      PALOER(N,2,N5,N4)=FSEDER*PALOH(NU(N2,N1),N2,N1)
      PFEOER(N,2,N5,N4)=FSEDER*PFEOH(NU(N2,N1),N2,N1)
      PCACER(N,2,N5,N4)=FSEDER*PCACO(NU(N2,N1),N2,N1)
      PCASER(N,2,N5,N4)=FSEDER*PCASO(NU(N2,N1),N2,N1)
      QALSER(N,2,N5,N4)=FSEDER*QALSI(NU(N2,N1),N2,N1)
      QFESER(N,2,N5,N4)=FSEDER*QFESI(NU(N2,N1),N2,N1)
      QCASER(N,2,N5,N4)=FSEDER*QCASI(NU(N2,N1),N2,N1)
      QMGSER(N,2,N5,N4)=FSEDER*QMGSI(NU(N2,N1),N2,N1)
      QNASER(N,2,N5,N4)=FSEDER*QNASI(NU(N2,N1),N2,N1)
      QKASER(N,2,N5,N4)=FSEDER*QKASI(NU(N2,N1),N2,N1)
      PALPER(N,2,N5,N4)=FSEDER*PALPO(NU(N2,N1),N2,N1)
      PFEPER(N,2,N5,N4)=FSEDER*PFEPO(NU(N2,N1),N2,N1)
      PCPDER(N,2,N5,N4)=FSEDER*PCAPD(NU(N2,N1),N2,N1)
      PCPHER(N,2,N5,N4)=FSEDER*PCAPH(NU(N2,N1),N2,N1)
      PCPMER(N,2,N5,N4)=FSEDER*PCAPM(NU(N2,N1),N2,N1)
      PALPEB(N,2,N5,N4)=FSEDER*PALPB(NU(N2,N1),N2,N1)
      PFEPEB(N,2,N5,N4)=FSEDER*PFEPB(NU(N2,N1),N2,N1)
      PCPDEB(N,2,N5,N4)=FSEDER*PCPDB(NU(N2,N1),N2,N1)
      PCPHEB(N,2,N5,N4)=FSEDER*PCPHB(NU(N2,N1),N2,N1)
      PCPMEB(N,2,N5,N4)=FSEDER*PCPMB(NU(N2,N1),N2,N1)
C
C     EROSION OF ORGANIC MATTER (g t-1)
C
C     sediment code
C        :OMC,OMN,OMP=microbial C,N,P
C        :ORC,ORN,ORP=microbial residue C,N,P
C        :OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate 
C        :OSC,OSA,OSN,OSP=SOC,colonized SOC,SON,SOP 
C           (K=0:woody litter, K=1:non-woody litter,
C            K=2:manure, K=3:POC, K=4:humus)
C
      DO 9480 K=0,5
      DO 9480 NO=1,7
      DO 9480 M=1,3
      OMCER(M,NO,K,N,2,N5,N4)=FSEDER*OMC(M,NO,K,NU(N2,N1),N2,N1)
      OMNER(M,NO,K,N,2,N5,N4)=FSEDER*OMN(M,NO,K,NU(N2,N1),N2,N1)
      OMPER(M,NO,K,N,2,N5,N4)=FSEDER*OMP(M,NO,K,NU(N2,N1),N2,N1)
9480  CONTINUE
      DO 9475 K=0,4
      DO 9470 M=1,2
      ORCER(M,K,N,2,N5,N4)=FSEDER*ORC(M,K,NU(N2,N1),N2,N1)
      ORNER(M,K,N,2,N5,N4)=FSEDER*ORN(M,K,NU(N2,N1),N2,N1)
      ORPER(M,K,N,2,N5,N4)=FSEDER*ORP(M,K,NU(N2,N1),N2,N1)
9470  CONTINUE
      OHCER(K,N,2,N5,N4)=FSEDER*OHC(K,NU(N2,N1),N2,N1)
      OHNER(K,N,2,N5,N4)=FSEDER*OHN(K,NU(N2,N1),N2,N1)
      OHPER(K,N,2,N5,N4)=FSEDER*OHP(K,NU(N2,N1),N2,N1)
      OHAER(K,N,2,N5,N4)=FSEDER*OHA(K,NU(N2,N1),N2,N1)
      DO 9465 M=1,5
      OSCER(M,K,N,2,N5,N4)=FSEDER*OSC(M,K,NU(N2,N1),N2,N1)
      OSAER(M,K,N,2,N5,N4)=FSEDER*OSA(M,K,NU(N2,N1),N2,N1)
      OSNER(M,K,N,2,N5,N4)=FSEDER*OSN(M,K,NU(N2,N1),N2,N1)
      OSPER(M,K,N,2,N5,N4)=FSEDER*OSP(M,K,NU(N2,N1),N2,N1)
9465  CONTINUE
9475  CONTINUE
      ELSE
C
C     SEDIMENT POOLS
C
      XSANER(N,2,N5,N4)=0.0
      XSILER(N,2,N5,N4)=0.0
      XCLAER(N,2,N5,N4)=0.0
      XCECER(N,2,N5,N4)=0.0
      XAECER(N,2,N5,N4)=0.0
C
C     FERTILIZER POOLS
C
      XNH4ER(N,2,N5,N4)=0.0
      XNH3ER(N,2,N5,N4)=0.0
      XNHUER(N,2,N5,N4)=0.0
      XNO3ER(N,2,N5,N4)=0.0
      XNH4EB(N,2,N5,N4)=0.0
      XNH3EB(N,2,N5,N4)=0.0
      XNHUEB(N,2,N5,N4)=0.0
      XNO3EB(N,2,N5,N4)=0.0
C
C     EXCHANGEABLE CATIONS AND ANIONS
C
      XN4ER(N,2,N5,N4)=0.0
      XNBER(N,2,N5,N4)=0.0
      XHYER(N,2,N5,N4)=0.0
      XALER(N,2,N5,N4)=0.0
      XFEER(N,2,N5,N4)=0.0
      XCAER(N,2,N5,N4)=0.0
      XMGER(N,2,N5,N4)=0.0
      XNAER(N,2,N5,N4)=0.0
      XKAER(N,2,N5,N4)=0.0
      XHCER(N,2,N5,N4)=0.0
      XAL2ER(N,2,N5,N4)=0.0
      XFE2ER(N,2,N5,N4)=0.0
      XOH0ER(N,2,N5,N4)=0.0
      XOH1ER(N,2,N5,N4)=0.0
      XOH2ER(N,2,N5,N4)=0.0
      XH1PER(N,2,N5,N4)=0.0
      XH2PER(N,2,N5,N4)=0.0
      XOH0EB(N,2,N5,N4)=0.0
      XOH1EB(N,2,N5,N4)=0.0
      XOH2EB(N,2,N5,N4)=0.0
      XH1PEB(N,2,N5,N4)=0.0
      XH2PEB(N,2,N5,N4)=0.0
C
C     PRECIPITATES
C
      PALOER(N,2,N5,N4)=0.0
      PFEOER(N,2,N5,N4)=0.0
      PCACER(N,2,N5,N4)=0.0
      PCASER(N,2,N5,N4)=0.0
      PALPER(N,2,N5,N4)=0.0
      PFEPER(N,2,N5,N4)=0.0
      PCPDER(N,2,N5,N4)=0.0
      PCPHER(N,2,N5,N4)=0.0
      PCPMER(N,2,N5,N4)=0.0
      PALPEB(N,2,N5,N4)=0.0
      PFEPEB(N,2,N5,N4)=0.0
      PCPDEB(N,2,N5,N4)=0.0
      PCPHEB(N,2,N5,N4)=0.0
      PCPMEB(N,2,N5,N4)=0.0
C
C     ORGANIC MATTER
C
      DO 8480 K=0,5
      DO 8480 NO=1,7
      DO 8480 M=1,3
      OMCER(M,NO,K,N,2,N5,N4)=0.0
      OMNER(M,NO,K,N,2,N5,N4)=0.0
      OMPER(M,NO,K,N,2,N5,N4)=0.0
8480  CONTINUE
      DO 8475 K=0,4
      DO 8470 M=1,2
      ORCER(M,K,N,2,N5,N4)=0.0
      ORNER(M,K,N,2,N5,N4)=0.0
      ORPER(M,K,N,2,N5,N4)=0.0
8470  CONTINUE
      OHCER(K,N,2,N5,N4)=0.0
      OHNER(K,N,2,N5,N4)=0.0
      OHPER(K,N,2,N5,N4)=0.0
      OHAER(K,N,2,N5,N4)=0.0
      DO 8465 M=1,5
      OSCER(M,K,N,2,N5,N4)=0.0
      OSAER(M,K,N,2,N5,N4)=0.0
      OSNER(M,K,N,2,N5,N4)=0.0
      OSPER(M,K,N,2,N5,N4)=0.0
8465  CONTINUE
8475  CONTINUE
      ENDIF
      IF(NN.EQ.2)THEN
      IF(N4B.GT.0.AND.N5B.GT.0)THEN
      FSEDER=AMIN1(1.0,XSEDER(N,1,N5B,N4B)/BKVLNU(N2,N1))
C
C     SOIL MINERAL EROSION
C
C     sediment code:SAN=sand,SIL=silt,CLA=clay (Mg t-1)
C                  :CEC=cation exchange capacity (mol t-1)
C                  :AEC=anion exchange capacity (mol t-1)
C
      XSANER(N,1,N5B,N4B)=FSEDER*SAND(NU(N2,N1),N2,N1)
      XSILER(N,1,N5B,N4B)=FSEDER*SILT(NU(N2,N1),N2,N1)
      XCLAER(N,1,N5B,N4B)=FSEDER*CLAY(NU(N2,N1),N2,N1)
      XCECER(N,1,N5B,N4B)=FSEDER*XCEC(NU(N2,N1),N2,N1)
      XAECER(N,1,N5B,N4B)=FSEDER*XAEC(NU(N2,N1),N2,N1)
C
C     EROSION FROM FERTILIZER POOLS (mol t-1)
C
C     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3 in non-band
C                  :NH4B,NH3B,NHUB,NO3B=NH4,NH3,urea,NO3 in band
C
      XNH4ER(N,1,N5B,N4B)=FSEDER*ZNH4FA(NU(N2,N1),N2,N1)
      XNH3ER(N,1,N5B,N4B)=FSEDER*ZNH3FA(NU(N2,N1),N2,N1)
      XNHUER(N,1,N5B,N4B)=FSEDER*ZNHUFA(NU(N2,N1),N2,N1)
      XNO3ER(N,1,N5B,N4B)=FSEDER*ZNO3FA(NU(N2,N1),N2,N1)
      XNH4EB(N,1,N5B,N4B)=FSEDER*ZNH4FB(NU(N2,N1),N2,N1)
      XNH3EB(N,1,N5B,N4B)=FSEDER*ZNH3FB(NU(N2,N1),N2,N1)
      XNHUEB(N,1,N5B,N4B)=FSEDER*ZNHUFB(NU(N2,N1),N2,N1)
      XNO3EB(N,1,N5B,N4B)=FSEDER*ZNO3FB(NU(N2,N1),N2,N1)
C
C     EXCHANGEABLE CATION AND ANION EROSION (mol t-1)
C
C     sediment code
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC
C           =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3 
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
C       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C
      XN4ER(N,1,N5B,N4B)=FSEDER*XN4(NU(N2,N1),N2,N1)
      XNBER(N,1,N5B,N4B)=FSEDER*XNB(NU(N2,N1),N2,N1)
      XHYER(N,1,N5B,N4B)=FSEDER*XHY(NU(N2,N1),N2,N1)
      XALER(N,1,N5B,N4B)=FSEDER*XAL(NU(N2,N1),N2,N1)
      XFEER(N,1,N5B,N4B)=FSEDER*XFE(NU(N2,N1),N2,N1)
      XCAER(N,1,N5B,N4B)=FSEDER*XCA(NU(N2,N1),N2,N1)
      XMGER(N,1,N5B,N4B)=FSEDER*XMG(NU(N2,N1),N2,N1)
      XNAER(N,1,N5B,N4B)=FSEDER*XNA(NU(N2,N1),N2,N1)
      XKAER(N,1,N5B,N4B)=FSEDER*XKA(NU(N2,N1),N2,N1)
      XHCER(N,1,N5B,N4B)=FSEDER*XHC(NU(N2,N1),N2,N1)
      XOH0ER(N,1,N5B,N4B)=FSEDER*XOH0(NU(N2,N1),N2,N1)
      XOH1ER(N,1,N5B,N4B)=FSEDER*XOH1(NU(N2,N1),N2,N1)
      XOH2ER(N,1,N5B,N4B)=FSEDER*XOH2(NU(N2,N1),N2,N1)
      XH1PER(N,1,N5B,N4B)=FSEDER*XH1P(NU(N2,N1),N2,N1)
      XH2PER(N,1,N5B,N4B)=FSEDER*XH2P(NU(N2,N1),N2,N1)
      XOH0EB(N,1,N5B,N4B)=FSEDER*XOH0B(NU(N2,N1),N2,N1)
      XOH1EB(N,1,N5B,N4B)=FSEDER*XOH1B(NU(N2,N1),N2,N1)
      XOH2EB(N,1,N5B,N4B)=FSEDER*XOH2B(NU(N2,N1),N2,N1)
      XH1PEB(N,1,N5B,N4B)=FSEDER*XH1PB(NU(N2,N1),N2,N1)
      XH2PEB(N,1,N5B,N4B)=FSEDER*XH2PB(NU(N2,N1),N2,N1)
C
C     EROSION OF PRECIPITATES (mol t-1)
C
C     sediment code
C       :PALO,PFEO=precipitated AlOH,FeOH 
C       :PCAC,PCAS=precipitated CaCO3,CaSO4
C       :PALP,PFEP=precipitated AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precipitated AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precipitated CaH2PO4,CaHPO4,apatite 
C           in non-band
C       :PCPMB,PCPDB,PCPHB= precipitated CaH2PO4,CaHPO4,apatite 
C           in band
C
      PALOER(N,1,N5B,N4B)=FSEDER*PALOH(NU(N2,N1),N2,N1)
      PFEOER(N,1,N5B,N4B)=FSEDER*PFEOH(NU(N2,N1),N2,N1)
      PCACER(N,1,N5B,N4B)=FSEDER*PCACO(NU(N2,N1),N2,N1)
      PCASER(N,1,N5B,N4B)=FSEDER*PCASO(NU(N2,N1),N2,N1)
      QALSER(N,1,N5B,N4B)=FSEDER*QALSI(NU(N2,N1),N2,N1)
      QFESER(N,1,N5B,N4B)=FSEDER*QFESI(NU(N2,N1),N2,N1)
      QCASER(N,1,N5B,N4B)=FSEDER*QCASI(NU(N2,N1),N2,N1)
      QMGSER(N,1,N5B,N4B)=FSEDER*QMGSI(NU(N2,N1),N2,N1)
      QNASER(N,1,N5B,N4B)=FSEDER*QNASI(NU(N2,N1),N2,N1)
      QKASER(N,1,N5B,N4B)=FSEDER*QKASI(NU(N2,N1),N2,N1)
      PALPER(N,1,N5B,N4B)=FSEDER*PALPO(NU(N2,N1),N2,N1)
      PFEPER(N,1,N5B,N4B)=FSEDER*PFEPO(NU(N2,N1),N2,N1)
      PCPDER(N,1,N5B,N4B)=FSEDER*PCAPD(NU(N2,N1),N2,N1)
      PCPHER(N,1,N5B,N4B)=FSEDER*PCAPH(NU(N2,N1),N2,N1)
      PCPMER(N,1,N5B,N4B)=FSEDER*PCAPM(NU(N2,N1),N2,N1)
      PALPEB(N,1,N5B,N4B)=FSEDER*PALPB(NU(N2,N1),N2,N1)
      PFEPEB(N,1,N5B,N4B)=FSEDER*PFEPB(NU(N2,N1),N2,N1)
      PCPDEB(N,1,N5B,N4B)=FSEDER*PCPDB(NU(N2,N1),N2,N1)
      PCPHEB(N,1,N5B,N4B)=FSEDER*PCPHB(NU(N2,N1),N2,N1)
      PCPMEB(N,1,N5B,N4B)=FSEDER*PCPMB(NU(N2,N1),N2,N1)
C
C     EROSION OF ORGANIC MATTER (g t-1)
C
C     sediment code:OMC=microbial biomass
C        :ORC=microbial residue
C        :OQC,OQCH=DOC in micropores,macropores
C        :OQA,OQAH=acetate in micropores,macropores
C        :OHC,OHA=adsorbed SOC,acetate 
C        :OSC=SOC(K=0:woody litter, K=1:non-woody litter,
C        :K=2:manure, K=3:POC, K=4:humus)
C
      DO 7480 K=0,5
      DO 7480 NO=1,7
      DO 7480 M=1,3
      OMCER(M,NO,K,N,1,N5B,N4B)=FSEDER*OMC(M,NO,K,NU(N2,N1),N2,N1)
      OMNER(M,NO,K,N,1,N5B,N4B)=FSEDER*OMN(M,NO,K,NU(N2,N1),N2,N1)
      OMPER(M,NO,K,N,1,N5B,N4B)=FSEDER*OMP(M,NO,K,NU(N2,N1),N2,N1)
7480  CONTINUE
      DO 7475 K=0,4
      DO 7470 M=1,2
      ORCER(M,K,N,1,N5B,N4B)=FSEDER*ORC(M,K,NU(N2,N1),N2,N1)
      ORNER(M,K,N,1,N5B,N4B)=FSEDER*ORN(M,K,NU(N2,N1),N2,N1)
      ORPER(M,K,N,1,N5B,N4B)=FSEDER*ORP(M,K,NU(N2,N1),N2,N1)
7470  CONTINUE
      OHCER(K,N,1,N5B,N4B)=FSEDER*OHC(K,NU(N2,N1),N2,N1)
      OHNER(K,N,1,N5B,N4B)=FSEDER*OHN(K,NU(N2,N1),N2,N1)
      OHPER(K,N,1,N5B,N4B)=FSEDER*OHP(K,NU(N2,N1),N2,N1)
      OHAER(K,N,1,N5B,N4B)=FSEDER*OHA(K,NU(N2,N1),N2,N1)
      DO 7465 M=1,5
      OSCER(M,K,N,1,N5B,N4B)=FSEDER*OSC(M,K,NU(N2,N1),N2,N1)
      OSAER(M,K,N,1,N5B,N4B)=FSEDER*OSA(M,K,NU(N2,N1),N2,N1)
      OSNER(M,K,N,1,N5B,N4B)=FSEDER*OSN(M,K,NU(N2,N1),N2,N1)
      OSPER(M,K,N,1,N5B,N4B)=FSEDER*OSP(M,K,NU(N2,N1),N2,N1)
7465  CONTINUE
7475  CONTINUE
      ELSE
      XSANER(N,1,N5B,N4B)=0.0
      XSILER(N,1,N5B,N4B)=0.0
      XCLAER(N,1,N5B,N4B)=0.0
      XCECER(N,1,N5B,N4B)=0.0
      XAECER(N,1,N5B,N4B)=0.0
C
C     FERTILIZER POOLS
C
      XNH4ER(N,1,N5B,N4B)=0.0
      XNH3ER(N,1,N5B,N4B)=0.0
      XNHUER(N,1,N5B,N4B)=0.0
      XNO3ER(N,1,N5B,N4B)=0.0
      XNH4EB(N,1,N5B,N4B)=0.0
      XNH3EB(N,1,N5B,N4B)=0.0
      XNHUEB(N,1,N5B,N4B)=0.0
      XNO3EB(N,1,N5B,N4B)=0.0
C
C     EXCHANGEABLE CATIONS AND ANIONS
C
      XN4ER(N,1,N5B,N4B)=0.0
      XNBER(N,1,N5B,N4B)=0.0
      XHYER(N,1,N5B,N4B)=0.0
      XALER(N,1,N5B,N4B)=0.0
      XFEER(N,1,N5B,N4B)=0.0
      XCAER(N,1,N5B,N4B)=0.0
      XMGER(N,1,N5B,N4B)=0.0
      XNAER(N,1,N5B,N4B)=0.0
      XKAER(N,1,N5B,N4B)=0.0
      XHCER(N,1,N5B,N4B)=0.0
      XAL2ER(N,1,N5B,N4B)=0.0
      XFE2ER(N,1,N5B,N4B)=0.0
      XOH0ER(N,1,N5B,N4B)=0.0
      XOH1ER(N,1,N5B,N4B)=0.0
      XOH2ER(N,1,N5B,N4B)=0.0
      XH1PER(N,1,N5B,N4B)=0.0
      XH2PER(N,1,N5B,N4B)=0.0
      XOH0EB(N,1,N5B,N4B)=0.0
      XOH1EB(N,1,N5B,N4B)=0.0
      XOH2EB(N,1,N5B,N4B)=0.0
      XH1PEB(N,1,N5B,N4B)=0.0
      XH2PEB(N,1,N5B,N4B)=0.0
C
C     PRECIPITATES
C
      PALOER(N,1,N5B,N4B)=0.0
      PFEOER(N,1,N5B,N4B)=0.0
      PCACER(N,1,N5B,N4B)=0.0
      PCASER(N,1,N5B,N4B)=0.0
      PALPER(N,1,N5B,N4B)=0.0
      PFEPER(N,1,N5B,N4B)=0.0
      PCPDER(N,1,N5B,N4B)=0.0
      PCPHER(N,1,N5B,N4B)=0.0
      PCPMER(N,1,N5B,N4B)=0.0
      PALPEB(N,1,N5B,N4B)=0.0
      PFEPEB(N,1,N5B,N4B)=0.0
      PCPDEB(N,1,N5B,N4B)=0.0
      PCPHEB(N,1,N5B,N4B)=0.0
      PCPMEB(N,1,N5B,N4B)=0.0
C
C     ORGANIC MATTER
C
      DO 6480 K=0,5
      DO 6480 NO=1,7
      DO 6480 M=1,3
      OMCER(M,NO,K,N,1,N5B,N4B)=0.0
      OMNER(M,NO,K,N,1,N5B,N4B)=0.0
      OMPER(M,NO,K,N,1,N5B,N4B)=0.0
6480  CONTINUE
      DO 6475 K=0,4
      DO 6470 M=1,2
      ORCER(M,K,N,1,N5B,N4B)=0.0
      ORNER(M,K,N,1,N5B,N4B)=0.0
      ORPER(M,K,N,1,N5B,N4B)=0.0
6470  CONTINUE
      OHCER(K,N,1,N5B,N4B)=0.0
      OHNER(K,N,1,N5B,N4B)=0.0
      OHPER(K,N,1,N5B,N4B)=0.0
      OHAER(K,N,1,N5B,N4B)=0.0
      DO 6465 M=1,5
      OSCER(M,K,N,1,N5B,N4B)=0.0
      OSAER(M,K,N,1,N5B,N4B)=0.0
      OSNER(M,K,N,1,N5B,N4B)=0.0
      OSPER(M,K,N,1,N5B,N4B)=0.0
6465  CONTINUE
6475  CONTINUE
      ENDIF
      ENDIF
9450  CONTINUE
9485  CONTINUE
      ENDIF
9490  CONTINUE
9495  CONTINUE
C
C     EXTERNAL BOUNDARY SEDIMENT FLUXES
C
C     IERSNG=options for disturbance effects on soil profile layer
C        depths and contents:
C           :-1=no effects
C           :0=freeze-thaw
C           :1=freeze-thaw+erosion
C           :2=freeze-thaw+SOM gain or loss
C           :3=freeze-thaw+erosion+SOM gain or loss 
C     N2,N1=NY,NX of source grid cell
C     N5,N4=NY,NX of destination grid cell E or S
C     N5B,N4B=NY,NX of destination grid cell W or N
C     M5,M4=NY,NX of destination grid cell 
C     NN=boundary:N=1:NN=1 east,NN=2 west, N=2:NN=1 south,NN=2 north
C
      DO 8995 NX=NHW,NHE
      DO 8990 NY=NVN,NVS
      IF((IERSNG.EQ.1.OR.IERSNG.EQ.3)
     2.AND.BKDS(NU(NY,NX),NY,NX).GT.ZERO)THEN
      N1=NX
      N2=NY
      DO 8980 N=1,2
      DO 8975 NN=1,2
      IF(N.EQ.1)THEN
      N4=NX+1
      N5=NY
      IF(NN.EQ.1)THEN
      IF(NX.EQ.NHE)THEN
      M1=NX
      M2=NY
      M4=NX+1
      M5=NY
      XN=-1.0
      RCHQF=RCHQE(M2,M1)
      ELSE
      GO TO 8975
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NX.EQ.NHW)THEN
      M1=NX
      M2=NY
      M4=NX
      M5=NY
      XN=1.0
      RCHQF=RCHQW(M5,M4)
      ELSE
      GO TO 8975
      ENDIF
      ENDIF
      ELSEIF(N.EQ.2)THEN
      N4=NX
      N5=NY+1
      IF(NN.EQ.1)THEN
      IF(NY.EQ.NVS)THEN
      M1=NX
      M2=NY
      M4=NX
      M5=NY+1
      XN=-1.0
      RCHQF=RCHQS(M2,M1)
      ELSE
      GO TO 8975
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NY.EQ.NVN)THEN
      M1=NX
      M2=NY
      M4=NX
      M5=NY
      XN=1.0
      RCHQF=RCHQN(M5,M4)
      ELSE
      GO TO 8975
      ENDIF
      ENDIF
      ENDIF
C
C     IF NO EROSION
C
      IF(IRCHG(NN,N,N2,N1).EQ.0.OR.RCHQF.EQ.0.0
     2.OR.ABS(XSEDER(N,NN,M5,M4)).LE.ZEROS(N2,N1))THEN
C
C     SOIL MINERALS
C 
      XSANER(N,NN,M5,M4)=0.0
      XSILER(N,NN,M5,M4)=0.0
      XCLAER(N,NN,M5,M4)=0.0
      XCECER(N,NN,M5,M4)=0.0
      XAECER(N,NN,M5,M4)=0.0
C
C     FERTILIZER POOLS
C
      XNH4ER(N,NN,M5,M4)=0.0
      XNH3ER(N,NN,M5,M4)=0.0
      XNHUER(N,NN,M5,M4)=0.0
      XNO3ER(N,NN,M5,M4)=0.0
      XNH4EB(N,NN,M5,M4)=0.0
      XNH3EB(N,NN,M5,M4)=0.0
      XNHUEB(N,NN,M5,M4)=0.0
      XNO3EB(N,NN,M5,M4)=0.0
C
C     EXCHANGEABLE CATIONS AND ANIONS
C
      XN4ER(N,NN,M5,M4)=0.0
      XNBER(N,NN,M5,M4)=0.0
      XHYER(N,NN,M5,M4)=0.0
      XALER(N,NN,M5,M4)=0.0
      XFEER(N,NN,M5,M4)=0.0
      XCAER(N,NN,M5,M4)=0.0
      XMGER(N,NN,M5,M4)=0.0
      XNAER(N,NN,M5,M4)=0.0
      XKAER(N,NN,M5,M4)=0.0
      XHCER(N,NN,M5,M4)=0.0
      XAL2ER(N,NN,M5,M4)=0.0
      XFE2ER(N,NN,M5,M4)=0.0
      XOH0ER(N,NN,M5,M4)=0.0
      XOH1ER(N,NN,M5,M4)=0.0
      XOH2ER(N,NN,M5,M4)=0.0
      XH1PER(N,NN,M5,M4)=0.0
      XH2PER(N,NN,M5,M4)=0.0
      XOH0EB(N,NN,M5,M4)=0.0
      XOH1EB(N,NN,M5,M4)=0.0
      XOH2EB(N,NN,M5,M4)=0.0
      XH1PEB(N,NN,M5,M4)=0.0
      XH2PEB(N,NN,M5,M4)=0.0
C
C     PRECIPITATES
C
      PALOER(N,NN,M5,M4)=0.0
      PFEOER(N,NN,M5,M4)=0.0
      PCACER(N,NN,M5,M4)=0.0
      PCASER(N,NN,M5,M4)=0.0
      PALPER(N,NN,M5,M4)=0.0
      PFEPER(N,NN,M5,M4)=0.0
      PCPDER(N,NN,M5,M4)=0.0
      PCPHER(N,NN,M5,M4)=0.0
      PCPMER(N,NN,M5,M4)=0.0
      PALPEB(N,NN,M5,M4)=0.0
      PFEPEB(N,NN,M5,M4)=0.0
      PCPDEB(N,NN,M5,M4)=0.0
      PCPHEB(N,NN,M5,M4)=0.0
      PCPMEB(N,NN,M5,M4)=0.0
C
C     ORGANIC MATTER
C
      DO 5480 K=0,5
      DO 5480 NO=1,7
      DO 5480 M=1,3
      OMCER(M,NO,K,N,NN,M5,M4)=0.0
      OMNER(M,NO,K,N,NN,M5,M4)=0.0
      OMPER(M,NO,K,N,NN,M5,M4)=0.0
5480  CONTINUE
      DO 5475 K=0,4
      DO 5470 M=1,2
      ORCER(M,K,N,NN,M5,M4)=0.0
      ORNER(M,K,N,NN,M5,M4)=0.0
      ORPER(M,K,N,NN,M5,M4)=0.0
5470  CONTINUE
      OHCER(K,N,NN,M5,M4)=0.0
      OHNER(K,N,NN,M5,M4)=0.0
      OHPER(K,N,NN,M5,M4)=0.0
      OHAER(K,N,NN,M5,M4)=0.0
      DO 5465 M=1,5
      OSCER(M,K,N,NN,M5,M4)=0.0
      OSAER(M,K,N,NN,M5,M4)=0.0
      OSNER(M,K,N,NN,M5,M4)=0.0
      OSPER(M,K,N,NN,M5,M4)=0.0
5465  CONTINUE
5475  CONTINUE
C
C     CALCULATE FRACTION OF SURACE MATERIAL ERODED
C
      ELSE
C
C     XSEDER=cumulative sediment transport in WE,NS directions (Mg t-1)
C     BKVLNU=mass of surface soil layer from �hour1.f� (Mg)
C     FSEDER=fraction of soil surface layer mass transported 
C        by erosion (t-1)
C     X*ER,X*EB= sediment flux in non-band,band used in �redist.f�
C        (Mg,g or mol t-1) 
C
      FSEDER=AMIN1(1.0,XSEDER(N,NN,N5,N4)/BKVLNU(N2,N1))
C
C     SOIL MINERAL EROSION
C
C     sediment code:SAN=sand,SIL=silt,CLA=clay (Mg t-1)
C                  :CEC=cation exchange capacity (mol t-1)
C                  :AEC=anion exchange capacity (mol t-1)
C
      XSANER(N,NN,M5,M4)=FSEDER*SAND(NU(N2,N1),N2,N1)
      XSILER(N,NN,M5,M4)=FSEDER*SILT(NU(N2,N1),N2,N1)
      XCLAER(N,NN,M5,M4)=FSEDER*CLAY(NU(N2,N1),N2,N1)
      XCECER(N,NN,M5,M4)=FSEDER*XCEC(NU(N2,N1),N2,N1)
      XAECER(N,NN,M5,M4)=FSEDER*XAEC(NU(N2,N1),N2,N1)
C
C     EROSION FROM FERTILIZER POOLS (mol t-1)
C
C     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3 
C        in non-band *R and band *B
C
      XNH4ER(N,NN,M5,M4)=FSEDER*ZNH4FA(NU(N2,N1),N2,N1)
      XNH3ER(N,NN,M5,M4)=FSEDER*ZNH3FA(NU(N2,N1),N2,N1)
      XNHUER(N,NN,M5,M4)=FSEDER*ZNHUFA(NU(N2,N1),N2,N1)
      XNO3ER(N,NN,M5,M4)=FSEDER*ZNO3FA(NU(N2,N1),N2,N1)
      XNH4EB(N,NN,M5,M4)=FSEDER*ZNH4FB(NU(N2,N1),N2,N1)
      XNH3EB(N,NN,M5,M4)=FSEDER*ZNH3FB(NU(N2,N1),N2,N1)
      XNHUEB(N,NN,M5,M4)=FSEDER*ZNHUFB(NU(N2,N1),N2,N1)
      XNO3EB(N,NN,M5,M4)=FSEDER*ZNO3FB(NU(N2,N1),N2,N1)
C
C     EXCHANGEABLE CATION AND ANION EROSION (mol t-1)
C
C     sediment code
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC
C           =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3 
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
C       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C
      XN4ER(N,NN,M5,M4)=FSEDER*XN4(NU(N2,N1),N2,N1)
      XNBER(N,NN,M5,M4)=FSEDER*XNB(NU(N2,N1),N2,N1)
      XHYER(N,NN,M5,M4)=FSEDER*XHY(NU(N2,N1),N2,N1)
      XALER(N,NN,M5,M4)=FSEDER*XAL(NU(N2,N1),N2,N1)
      XFEER(N,NN,M5,M4)=FSEDER*XFE(NU(N2,N1),N2,N1)
      XCAER(N,NN,M5,M4)=FSEDER*XCA(NU(N2,N1),N2,N1)
      XMGER(N,NN,M5,M4)=FSEDER*XMG(NU(N2,N1),N2,N1)
      XNAER(N,NN,M5,M4)=FSEDER*XNA(NU(N2,N1),N2,N1)
      XKAER(N,NN,M5,M4)=FSEDER*XKA(NU(N2,N1),N2,N1)
      XHCER(N,NN,M5,M4)=FSEDER*XHC(NU(N2,N1),N2,N1)
      XOH0ER(N,NN,M5,M4)=FSEDER*XOH0(NU(N2,N1),N2,N1)
      XOH1ER(N,NN,M5,M4)=FSEDER*XOH1(NU(N2,N1),N2,N1)
      XOH2ER(N,NN,M5,M4)=FSEDER*XOH2(NU(N2,N1),N2,N1)
      XH1PER(N,NN,M5,M4)=FSEDER*XH1P(NU(N2,N1),N2,N1)
      XH2PER(N,NN,M5,M4)=FSEDER*XH2P(NU(N2,N1),N2,N1)
      XOH0EB(N,NN,M5,M4)=FSEDER*XOH0B(NU(N2,N1),N2,N1)
      XOH1EB(N,NN,M5,M4)=FSEDER*XOH1B(NU(N2,N1),N2,N1)
      XOH2EB(N,NN,M5,M4)=FSEDER*XOH2B(NU(N2,N1),N2,N1)
      XH1PEB(N,NN,M5,M4)=FSEDER*XH1PB(NU(N2,N1),N2,N1)
      XH2PEB(N,NN,M5,M4)=FSEDER*XH2PB(NU(N2,N1),N2,N1)
C
C     EROSION OF PRECIPITATES (mol t-1)
C
C     sediment code
C       :PALO,PFEO=precipitated AlOH,FeOH 
C       :PCAC,PCAS=precipitated CaCO3,CaSO4
C       :PALP,PFEP=precipitated AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precipitated AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precipitated CaH2PO4,CaHPO4,apatite 
C           in non-band
C       :PCPMB,PCPDB,PCPHB=precipitated CaH2PO4,CaHPO4,apatite 
C           in band
C
      PALOER(N,NN,M5,M4)=FSEDER*PALOH(NU(N2,N1),N2,N1)
      PFEOER(N,NN,M5,M4)=FSEDER*PFEOH(NU(N2,N1),N2,N1)
      PCACER(N,NN,M5,M4)=FSEDER*PCACO(NU(N2,N1),N2,N1)
      PCASER(N,NN,M5,M4)=FSEDER*PCASO(NU(N2,N1),N2,N1)
      QALSER(N,NN,M5,M4)=FSEDER*QALSI(NU(N2,N1),N2,N1)
      QFESER(N,NN,M5,M4)=FSEDER*QFESI(NU(N2,N1),N2,N1)
      QCASER(N,NN,M5,M4)=FSEDER*QCASI(NU(N2,N1),N2,N1)
      QMGSER(N,NN,M5,M4)=FSEDER*QMGSI(NU(N2,N1),N2,N1)
      QNASER(N,NN,M5,M4)=FSEDER*QNASI(NU(N2,N1),N2,N1)
      QKASER(N,NN,M5,M4)=FSEDER*QKASI(NU(N2,N1),N2,N1)
      PALPER(N,NN,M5,M4)=FSEDER*PALPO(NU(N2,N1),N2,N1)
      PFEPER(N,NN,M5,M4)=FSEDER*PFEPO(NU(N2,N1),N2,N1)
      PCPDER(N,NN,M5,M4)=FSEDER*PCAPD(NU(N2,N1),N2,N1)
      PCPHER(N,NN,M5,M4)=FSEDER*PCAPH(NU(N2,N1),N2,N1)
      PCPMER(N,NN,M5,M4)=FSEDER*PCAPM(NU(N2,N1),N2,N1)
      PALPEB(N,NN,M5,M4)=FSEDER*PALPB(NU(N2,N1),N2,N1)
      PFEPEB(N,NN,M5,M4)=FSEDER*PFEPB(NU(N2,N1),N2,N1)
      PCPDEB(N,NN,M5,M4)=FSEDER*PCPDB(NU(N2,N1),N2,N1)
      PCPHEB(N,NN,M5,M4)=FSEDER*PCPHB(NU(N2,N1),N2,N1)
      PCPMEB(N,NN,M5,M4)=FSEDER*PCPMB(NU(N2,N1),N2,N1)
C
C     EROSION OF ORGANIC MATTER (g t-1)
C
C     sediment code:OMC=microbial biomass
C        :ORC=microbial residue
C        :OQC,OQCH=DOC in micropores,macropores
C        :OQA,OQAH=acetate in micropores,macropores
C        :OHC,OHA=adsorbed SOC,acetate 
C        :OSC=SOC(K=0:woody litter, K=1:non-woody litter,
C           :K=2:manure, K=3:POC, K=4:humus)
C
      DO 4880 K=0,5
      DO 4880 NO=1,7
      DO 4880 M=1,3
      OMCER(M,NO,K,N,NN,M5,M4)=FSEDER*OMC(M,NO,K,NU(N2,N1),N2,N1)
      OMNER(M,NO,K,N,NN,M5,M4)=FSEDER*OMN(M,NO,K,NU(N2,N1),N2,N1)
      OMPER(M,NO,K,N,NN,M5,M4)=FSEDER*OMP(M,NO,K,NU(N2,N1),N2,N1)
4880  CONTINUE
      DO 4875 K=0,4
      DO 4870 M=1,2
      ORCER(M,K,N,NN,M5,M4)=FSEDER*ORC(M,K,NU(N2,N1),N2,N1)
      ORNER(M,K,N,NN,M5,M4)=FSEDER*ORN(M,K,NU(N2,N1),N2,N1)
      ORPER(M,K,N,NN,M5,M4)=FSEDER*ORP(M,K,NU(N2,N1),N2,N1)
4870  CONTINUE
      OHCER(K,N,NN,M5,M4)=FSEDER*OHC(K,NU(N2,N1),N2,N1)
      OHNER(K,N,NN,M5,M4)=FSEDER*OHN(K,NU(N2,N1),N2,N1)
      OHPER(K,N,NN,M5,M4)=FSEDER*OHP(K,NU(N2,N1),N2,N1)
      OHAER(K,N,NN,M5,M4)=FSEDER*OHA(K,NU(N2,N1),N2,N1)
      DO 4865 M=1,5
      OSCER(M,K,N,NN,M5,M4)=FSEDER*OSC(M,K,NU(N2,N1),N2,N1)
      OSAER(M,K,N,NN,M5,M4)=FSEDER*OSA(M,K,NU(N2,N1),N2,N1)
      OSNER(M,K,N,NN,M5,M4)=FSEDER*OSN(M,K,NU(N2,N1),N2,N1)
      OSPER(M,K,N,NN,M5,M4)=FSEDER*OSP(M,K,NU(N2,N1),N2,N1)
4865  CONTINUE
4875  CONTINUE
      ENDIF
C     IF(ABS(XSEDER(N,NN,M5,M4)).GT.ZEROS(M5,M4))THEN
C     WRITE(*,1116)'EDGE',I,J,NFZ,N1,N2,N,XSEDER(N,NN,M5,M4),FSEDER
C    2,BKVLNU(N2,N1),ORGC(NU(N2,N1),N2,N1)
C    2,XCLAER(N,NN,M5,M4),CLAY(NU(N2,N1),N2,N1) 
C    3,ORGC(NU(N2,N1),N2,N1),DLYR(3,NU(N2,N1),N2,N1)
C    4,BKVL(NU(N2,N1),N2,N1)
1116  FORMAT(A8,6I4,30E12.4)
C     ENDIF
8975  CONTINUE
8980  CONTINUE
      ENDIF
8990  CONTINUE
8995  CONTINUE
      RETURN
      END

