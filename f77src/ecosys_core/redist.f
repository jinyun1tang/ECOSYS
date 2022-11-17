
      SUBROUTINE redist(I,J,NFZ,NHW,NHE,NVN,NVS)
C
C     THIS SUBROUTINE UPDATES SOIL STATE VARIABLES WITH WATER, HEAT,
C     C, N, P, SOLUTE FLUXES CALCULATED IN EARLIER SUBROUTINES
C
      include "parameters.h"
      include "blkc.h"
      include "blk1g.h"
      include "blk1n.h"
      include "blk1p.h"
      include "blk1cr.h"
      include "blk2a.h"
      include "blk2b.h"
      include "blk2c.h"
      include "blk3.h"
      include "blk5.h"
      include "blk8a.h"
      include "blk8b.h"
      include "blk9b.h"
      include "blk11a.h"
      include "blk11b.h"
      include "blk12a.h"
      include "blk12b.h"
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
      include "blk22a.h"
      include "blk22b.h"
      include "blk22c.h"
      include "blktest.h"
      DIMENSION TFLW(JZ,JY,JX),TFLWX(JZ,JY,JX),THFLW(JZ,JY,JX)
     1,TFLWH(JZ,JY,JX),TOCFLS(0:4,JZ,JY,JX),TONFLS(0:4,JZ,JY,JX)
     2,TOPFLS(0:4,JZ,JY,JX),TOAFLS(0:4,JZ,JY,JX),TCOFLS(JZ,JY,JX)
     3,TCHFLS(JZ,JY,JX),TOXFLS(JZ,JY,JX),TNXFLB(JZ,JY,JX)
     4,TNGFLS(JZ,JY,JX),TN2FLS(JZ,JY,JX),TN4FLS(JZ,JY,JX)
     5,TN4FLB(JZ,JY,JX),TN3FLS(JZ,JY,JX),TN3FLB(JZ,JY,JX)
     6,TNOFLS(JZ,JY,JX),TNOFLB(JZ,JY,JX),TPOFLS(JZ,JY,JX)
     7,TH2BFB(JZ,JY,JX),TNXFLS(JZ,JY,JX),TOCFHS(0:4,JZ,JY,JX)
     8,TONFHS(0:4,JZ,JY,JX),TOPFHS(0:4,JZ,JY,JX),TOAFHS(0:4,JZ,JY,JX)
     9,TCOFHS(JZ,JY,JX),TCHFHS(JZ,JY,JX),TNXFHB(JZ,JY,JX)
     2,TOXFHS(JZ,JY,JX),TNGFHS(JZ,JY,JX),TN2FHS(JZ,JY,JX)
     2,TN4FHS(JZ,JY,JX),TN4FHB(JZ,JY,JX),TN3FHS(JZ,JY,JX)
     3,TN3FHB(JZ,JY,JX),TNOFHS(JZ,JY,JX),TNOFHB(JZ,JY,JX)
     4,TPOFHS(JZ,JY,JX),TH2BHB(JZ,JY,JX),TNXFHS(JZ,JY,JX)
     5,TCOFLG(JZ,JY,JX),TCHFLG(JZ,JY,JX),TOXFLG(JZ,JY,JX)
     6,TNGFLG(JZ,JY,JX),TN2FLG(JZ,JY,JX),TNHFLG(JZ,JY,JX)
     7,TWFLFL(JZ,JY,JX),THFLFL(JZ,JY,JX),TWFLFH(JZ,JY,JX)
     8,TP1FLS(JZ,JY,JX),TP1FHS(JZ,JY,JX),TH1BFB(JZ,JY,JX)
     9,TH1BHB(JZ,JY,JX),VOLW1(JZ,JY,JX),VOLI1(JZ,JY,JX)
     1,VOLWH1(JZ,JY,JX),VOLIH1(JZ,JY,JX),CDPTHX(JZ,JY,JX)
     2,TWFLVL(JZ,JY,JX),THFLVL(JZ,JY,JX),CDPTHY(0:JZ,JY,JX)
     3,TFLV(JZ,JY,JX),TVPFLG(JZ,JY,JX)
      DIMENSION TQR(JY,JX),THQR(JY,JX),TQS(JY,JX),TQW(JY,JX)
     2,TQI(JY,JX),THQS(JY,JX),TFLWS(JS,JY,JX),TFLWW(JS,JY,JX)
     3,TFLWI(JS,JY,JX),THFLWW(JS,JY,JX),TOCQRS(0:4,JY,JX)
     1,TONQRS(0:4,JY,JX),TOPQRS(0:4,JY,JX),TOAQRS(0:4,JY,JX)
     2,TCOQRS(JY,JX),TCHQRS(JY,JX),TOXQRS(JY,JX),TQRH1P(JY,JX)
     3,TNGQRS(JY,JX),TN2QRS(JY,JX),TN4QRS(JY,JX),TN3QRS(JY,JX)
     4,TNOQRS(JY,JX),TPOQRS(JY,JX),TNXQRS(JY,JX),TQRAL(JY,JX)
     6,TQRFE(JY,JX),TQRHY(JY,JX),TQRCA(JY,JX),TQRMG(JY,JX)
     7,TQRNA(JY,JX),TQRKA(JY,JX),TQROH(JY,JX),TQRSO(JY,JX)
     8,TQRCL(JY,JX),TQRC3(JY,JX),TQRHC(JY,JX),TQRAL1(JY,JX)
     9,TQRAL2(JY,JX),TQRAL3(JY,JX),TQRAL4(JY,JX),TQRALS(JY,JX)
     1,TQRFE1(JY,JX),TQRFE2(JY,JX),TQRFE3(JY,JX),TQRFE4(JY,JX)
     2,TQRFES(JY,JX),TQRCAO(JY,JX),TQRCAC(JY,JX),TQRCAH(JY,JX)
     3,TQRCAS(JY,JX),TQRMGO(JY,JX),TQRMGC(JY,JX),TQRMGH(JY,JX)
     4,TQRMGS(JY,JX),TQRNAC(JY,JX),TQRNAS(JY,JX),TQRKAS(JY,JX)
     5,TQRH0P(JY,JX),TQRH3P(JY,JX),TQRF1P(JY,JX),TP1QRS(JY,JX)
     6,TQRF2P(JY,JX),TQRC0P(JY,JX),TQRC1P(JY,JX),TQRC2P(JY,JX)
     7,TQRM1P(JY,JX),TCOQSS(JY,JX),TCHQSS(JY,JX),TOXQSS(JY,JX)
     3,TNGQSS(JY,JX),TN2QSS(JY,JX),TN4QSS(JY,JX),TN3QSS(JY,JX)
     4,TNOQSS(JY,JX),TPOQSS(JY,JX),TP1QSS(JY,JX),TQSAL(JY,JX)
     6,TQSFE(JY,JX),TQSHY(JY,JX),TQSCA(JY,JX),TQSMG(JY,JX)
     7,TQSNA(JY,JX),TQSKA(JY,JX),TQSOH(JY,JX),TQSSO(JY,JX)
     8,TQSCL(JY,JX),TQSC3(JY,JX),TQSHC(JY,JX),TQSAL1(JY,JX)
     9,TQSAL2(JY,JX),TQSAL3(JY,JX),TQSAL4(JY,JX),TQSALS(JY,JX)
     1,TQSFE1(JY,JX),TQSFE2(JY,JX),TQSFE3(JY,JX),TQSFE4(JY,JX)
     2,TQSFES(JY,JX),TQSCAO(JY,JX),TQSCAC(JY,JX),TQSCAH(JY,JX)
     3,TQSCAS(JY,JX),TQSMGO(JY,JX),TQSMGC(JY,JX),TQSMGH(JY,JX)
     4,TQSMGS(JY,JX),TQSNAC(JY,JX),TQSNAS(JY,JX),TQSKAS(JY,JX)
     5,TQSH0P(JY,JX),TQSH1P(JY,JX),TQSH3P(JY,JX),TQSF1P(JY,JX)
     6,TQSF2P(JY,JX),TQSC0P(JY,JX),TQSC1P(JY,JX),TQSC2P(JY,JX)
     7,TQSM1P(JY,JX),TFLWV(JS,JY,JX)
      DIMENSION TALFLS(JZ,JY,JX),TFEFLS(JZ,JY,JX)
     1,TCAFLS(JZ,JY,JX),THYFLS(JZ,JY,JX),TMGFLS(JZ,JY,JX)
     2,TNAFLS(JZ,JY,JX),TKAFLS(JZ,JY,JX),TOHFLS(JZ,JY,JX)
     3,TSOFLS(JZ,JY,JX),TCLFLS(JZ,JY,JX),TC3FLS(JZ,JY,JX)
     4,THCFLS(JZ,JY,JX),TAL1FS(JZ,JY,JX),TAL2FS(JZ,JY,JX)
     5,TAL3FS(JZ,JY,JX),TAL4FS(JZ,JY,JX),TALSFS(JZ,JY,JX)
     6,TFE1FS(JZ,JY,JX),TFE2FS(JZ,JY,JX)
     7,TFE3FS(JZ,JY,JX),TFE4FS(JZ,JY,JX),TFESFS(JZ,JY,JX)
     8,TCAOFS(JZ,JY,JX),TCACFS(JZ,JY,JX),TCAHFS(JZ,JY,JX)
     9,TCASFS(JZ,JY,JX),TMGOFS(JZ,JY,JX),TMGCFS(JZ,JY,JX)
     1,TMGHFS(JZ,JY,JX),TMGSFS(JZ,JY,JX),TNACFS(JZ,JY,JX)
     2,TNASFS(JZ,JY,JX),TKASFS(JZ,JY,JX),TH0PFS(JZ,JY,JX)
     3,TH1PFS(JZ,JY,JX),TH3PFS(JZ,JY,JX),TF1PFS(JZ,JY,JX)
     4,TF2PFS(JZ,JY,JX),TC0PFS(JZ,JY,JX),TC1PFS(JZ,JY,JX)
     5,TC2PFS(JZ,JY,JX),TM1PFS(JZ,JY,JX),TH0BFB(JZ,JY,JX)
     6,TH3BFB(JZ,JY,JX),TF1BFB(JZ,JY,JX)
     7,TF2BFB(JZ,JY,JX),TC0BFB(JZ,JY,JX),TC1BFB(JZ,JY,JX)
     8,TC2BFB(JZ,JY,JX),TM1BFB(JZ,JY,JX)
      DIMENSION TALFHS(JZ,JY,JX),TFEFHS(JZ,JY,JX)
     1,THYFHS(JZ,JY,JX),TCAFHS(JZ,JY,JX),TMGFHS(JZ,JY,JX)
     2,TNAFHS(JZ,JY,JX),TKAFHS(JZ,JY,JX),TOHFHS(JZ,JY,JX)
     3,TSOFHS(JZ,JY,JX),TCLFHS(JZ,JY,JX),TC3FHS(JZ,JY,JX)
     4,THCFHS(JZ,JY,JX),TAL1HS(JZ,JY,JX),TAL2HS(JZ,JY,JX)
     5,TAL3HS(JZ,JY,JX),TAL4HS(JZ,JY,JX),TALSHS(JZ,JY,JX)
     6,TFE1HS(JZ,JY,JX),TFE2HS(JZ,JY,JX)
     7,TFE3HS(JZ,JY,JX),TFE4HS(JZ,JY,JX),TFESHS(JZ,JY,JX)
     8,TCAOHS(JZ,JY,JX),TCACHS(JZ,JY,JX),TCAHHS(JZ,JY,JX)
     9,TCASHS(JZ,JY,JX),TMGOHS(JZ,JY,JX),TMGCHS(JZ,JY,JX)
     1,TMGHHS(JZ,JY,JX),TMGSHS(JZ,JY,JX),TNACHS(JZ,JY,JX)
     2,TNASHS(JZ,JY,JX),TKASHS(JZ,JY,JX),TH0PHS(JZ,JY,JX)
     3,TH3PHS(JZ,JY,JX),TF1PHS(JZ,JY,JX)
     4,TF2PHS(JZ,JY,JX),TC0PHS(JZ,JY,JX),TC1PHS(JZ,JY,JX)
     5,TC2PHS(JZ,JY,JX),TM1PHS(JZ,JY,JX),TH0BHB(JZ,JY,JX)
     6,TH3BHB(JZ,JY,JX),TF1BHB(JZ,JY,JX)
     7,TF2BHB(JZ,JY,JX),TC0BHB(JZ,JY,JX),TC1BHB(JZ,JY,JX)
     8,TC2BHB(JZ,JY,JX),TM1BHB(JZ,JY,JX)
      DIMENSION TSANER(JY,JX),TSILER(JY,JX),TCLAER(JY,JX)
     2,TCECER(JY,JX),TAECER(JY,JX),TNH4ER(JY,JX),TNH3ER(JY,JX)
     3,TNHUER(JY,JX),TNO3ER(JY,JX),TNH4EB(JY,JX),TNH3EB(JY,JX)
     4,TNHUEB(JY,JX),TNO3EB(JY,JX),TN4ER(JY,JX),TNBER(JY,JX)
     5,THYER(JY,JX),TALER(JY,JX),TCAER(JY,JX),TMGER(JY,JX)
     6,TNAER(JY,JX),TKAER(JY,JX),THCER(JY,JX),TAL2ER(JY,JX)
     7,TOH0ER(JY,JX),TOH1ER(JY,JX),TOH2ER(JY,JX),TH1PER(JY,JX)
     8,TH2PER(JY,JX),TOH0EB(JY,JX),TOH1EB(JY,JX),TOH2EB(JY,JX)
     9,TH1PEB(JY,JX),TH2PEB(JY,JX),TALOER(JY,JX),TFEOER(JY,JX)
     1,TCACER(JY,JX),TCASER(JY,JX),TALPER(JY,JX),TFEPER(JY,JX)
     2,TCPDER(JY,JX),TCPHER(JY,JX),TCPMER(JY,JX),TALPEB(JY,JX)
     3,TFEPEB(JY,JX),TCPDEB(JY,JX),TCPHEB(JY,JX),TCPMEB(JY,JX)
     4,TOMCER(3,7,0:5,JY,JX),TOMNER(3,7,0:5,JY,JX)
     4,TOMPER(3,7,0:5,JY,JX),TFEER(JY,JX),TFE2ER(JY,JX)
     5,TORCER(2,0:4,JY,JX),TORNER(2,0:4,JY,JX),TORPER(2,0:4,JY,JX)
     6,TOHCER(0:4,JY,JX),TOHNER(0:4,JY,JX),TOHPER(0:4,JY,JX)
     7,TOHAER(0:4,JY,JX),TOSCER(5,0:4,JY,JX),TOSAER(5,0:4,JY,JX)
     8,TOSNER(5,0:4,JY,JX),TOSPER(5,0:4,JY,JX),TSEDER(JY,JX)
      DIMENSION TOMC(3,7,0:5),TOMN(3,7,0:5),TOMP(3,7,0:5),TORC(2,0:4)
     2,TORN(2,0:4),TORP(2,0:4),TOQC(0:4),TOQN(0:4),TOQP(0:4)
     3,TOQA(0:4),TOHC(0:4),TOHN(0:4),TOHP(0:4),TOHA(0:4),TOSC(5,0:4)
     4,TOSA(5,0:4),TOSN(5,0:4),TOSP(5,0:4),TOSGC(5,0:2),TOSGA(5,0:2)
     5,TOSGN(5,0:2),TOSGP(5,0:2),TOMGC(3,7,0:5),TOMGN(3,7,0:5)
     6,TOMGP(3,7,0:5),TORXC(2,0:2),TORXN(2,0:2),TORXP(2,0:2)
     7,TOQGC(0:2),TOQGN(0:2),TOQGP(0:2),TOQHC(0:2),TOQHN(0:2)
     8,TOQHP(0:2),TOHGC(0:2),TOHGN(0:2),TOHGP(0:2),TOHGA(0:2)
     9,TOQGA(0:2),TOQHA(0:2),THGQRS(JY,JX),THGFHS(JZ,JY,JX)
     1,THGFLG(JZ,JY,JX),THGFLS(JZ,JY,JX),TXCO2(JY,JX),DORGC(JZ,JY,JX)
     2,DORGE(JY,JX),OMCL(0:JZ,JY,JX),OMNL(0:JZ,JY,JX)
     3,DDLYR(0:JZ,6),DDLYX(0:JZ,6)
     4,DVOLW(JZ,JY,JX),DVOLI(JZ,JY,JX),VHCPWZ(JZ,JY,JX)
     5,DDLYRX(3),DDLYRY(JZ),TSEDSK(JY,JX),IFLGL(0:JZ,6),NUL(JY,JX)
      DIMENSION TCOBLS(JS,JY,JX),TCHBLS(JS,JY,JX),TOXBLS(JS,JY,JX)
     2,TNGBLS(JS,JY,JX),TN2BLS(JS,JY,JX),TN4BLW(JS,JY,JX)
     3,TN3BLW(JS,JY,JX),TNOBLW(JS,JY,JX),TH1PBS(JS,JY,JX)
     4,TH2PBS(JS,JY,JX),TALBLS(JS,JY,JX),TFEBLS(JS,JY,JX)
     5,THYBLS(JS,JY,JX),TCABLS(JS,JY,JX),TMGBLS(JS,JY,JX)
     6,TNABLS(JS,JY,JX),TKABLS(JS,JY,JX),TOHBLS(JS,JY,JX)
     7,TSOBLS(JS,JY,JX),TCLBLS(JS,JY,JX),TC3BLS(JS,JY,JX)
     8,THCBLS(JS,JY,JX),TAL1BS(JS,JY,JX),TAL2BS(JS,JY,JX)
     9,TAL3BS(JS,JY,JX),TAL4BS(JS,JY,JX),TALSBS(JS,JY,JX)
     1,TFE1BS(JS,JY,JX),TFE2BS(JS,JY,JX),TFE3BS(JS,JY,JX)
     2,TFE4BS(JS,JY,JX),TFESBS(JS,JY,JX),TCAOBS(JS,JY,JX)
     3,TCACBS(JS,JY,JX),TCAHBS(JS,JY,JX),TCASBS(JS,JY,JX)
     4,TMGOBS(JS,JY,JX),TMGCBS(JS,JY,JX),TMGHBS(JS,JY,JX)
     5,TMGSBS(JS,JY,JX),TNACBS(JS,JY,JX),TNASBS(JS,JY,JX)
     6,TKASBS(JS,JY,JX),TH0PBS(JS,JY,JX),TH3PBS(JS,JY,JX)
     7,TF1PBS(JS,JY,JX),TF2PBS(JS,JY,JX),TC0PBS(JS,JY,JX)
     8,TC1PBS(JS,JY,JX),TC2PBS(JS,JY,JX),TM1PBS(JS,JY,JX)
C
C     FCHMC,FCHGC=fraction of anaerobic C combustion to charcoal,CO2
C     FCHMN,FCHGN=fraction of anaerobic N combustion to NH4,
C        gaseous NOX
C     FCHMP,FCHGP=fraction of anaerobic P combustion to H2PO4,
C        gaseous POX
C     FCOMC,FCOGC=fraction of aerobic C combustion to charcoal,CO2
C     FCOMN,FGOMN=fraction of aerobic N combustion to NH4,
C        gaseous NOX
C     FCOMP,FCOGP=fraction of aerobic P combustion to H2PO4,
C        gaseous POX
C
      PARAMETER (FCHMC=0.3,FCHGC=1.0-FCHMC
     2,FCHMN=0.1,FCHGN=1.0-FCHMN
     3,FCHMP=0.9,FCHGP=1.0-FCHMP
     4,FCOMC=0.0,FCOGC=1.0-FCOMC
     5,FCOMN=0.1,FCOGN=1.0-FCOMN
     6,FCOMP=0.9,FCOGP=1.0-FCOMP)
      DATA SG/0.0/
      VOLISO=0.0
      UDVOLI=0.0
      UDLYXF=0.0
      TFLWT=0.0
      VOLPT=0.0
      VOLTT=0.0
      DO 9995 NX=NHW,NHE
      DO 9990 NY=NVN,NVS
      TXCO2(NY,NX)=0.0
      DORGE(NY,NX)=0.0
C
C     ADD MANURE LITTERFALL FROM EXTRACT.F TO SURFACE RESIDUE
C
C     OSC,OSN,OSP=SOC,SON,SOP
C     CSNMT,ZSNMT,PSNMT=organic C,N,P in manure litterfall
C        from extract.f
C     ZSNIT,PSNIT=inorganic N,P in manure litterfall
C        from extract.f
C
      IF(NFZ.EQ.1.AND.CSNMT(1,NY,NX).GT.ZEROS(NY,NX))THEN
      DO 6955 M=1,4
      OSC(M,2,0,NY,NX)=OSC(M,2,0,NY,NX)+CSNMT(M,NY,NX)
      OSA(M,2,0,NY,NX)=OSA(M,2,0,NY,NX)+CSNMT(M,NY,NX)
      OSN(M,2,0,NY,NX)=OSN(M,2,0,NY,NX)+ZSNMT(M,NY,NX)
      OSP(M,2,0,NY,NX)=OSP(M,2,0,NY,NX)+PSNMT(M,NY,NX)
C     WRITE(*,8487)'OSC0',I,J,NFZ,M,OSC(M,2,0,NY,NX)
C    2,OSN(M,2,0,NY,NX),OSP(M,2,0,NY,NX),CSNMT(M,NY,NX)
C    3,ZSNMT(M,NY,NX),PSNMT(M,NY,NX),ZSNIT(NY,NX),PSNIT(NY,NX)
8487  FORMAT(A8,4I4,20E12.4)
6955  CONTINUE
      ZNHUFA(0,NY,NX)=ZNHUFA(0,NY,NX)+ZSNIT(NY,NX)/14.0
      PCAPM(0,NY,NX)=PCAPM(0,NY,NX)+PSNIT(NY,NX)/62.0
      IUTYP(NY,NX)=0
      ENDIF
C
C     ADD PLANT LITTERFALL FROM EXTRACT.F TO SURFACE RESIDUE
C
C     OSC,OSN,OSP=SOC,SON,SOP
C     CSNT,ZSNT,PSNT=total C,N,P in plant litterfall from extract.f
C
      DO 6965 K=0,1
      DO 6965 M=1,5
      OSC(M,K,0,NY,NX)=OSC(M,K,0,NY,NX)+CSNT(M,K,0,NY,NX)
C     OSA(M,K,0,NY,NX)=OSA(M,K,0,NY,NX)+CSNT(M,K,0,NY,NX)*OMCI(1,K)
      OSN(M,K,0,NY,NX)=OSN(M,K,0,NY,NX)+ZSNT(M,K,0,NY,NX)
      OSP(M,K,0,NY,NX)=OSP(M,K,0,NY,NX)+PSNT(M,K,0,NY,NX)
      IF(M.LE.4)THEN
      ORGC(0,NY,NX)=ORGC(0,NY,NX)+CSNT(M,K,0,NY,NX)
      ELSE
      ORGCC(0,NY,NX)=ORGCC(0,NY,NX)+CSNT(M,K,0,NY,NX)
      ENDIF
C     IF(ICHKF.EQ.1)THEN
C     WRITE(*,8486)'OSC0',I,J,NFZ,K,M,OSC(M,K,0,NY,NX)
C    2,OSN(M,K,0,NY,NX),OSP(M,K,0,NY,NX),CSNT(M,K,0,NY,NX)
C    3,ZSNT(M,K,0,NY,NX),PSNT(M,K,0,NY,NX)
8486  FORMAT(A8,5I4,20E12.4)
C     ENDIF
6965  CONTINUE
C
C     GAS AND SOLUTE EXCHANGE WITHIN SURFACE LITTER ADDED TO ECOSYSTEM
C     TOTALS FOR CALCULATING COMPETITION CONSTRAINTS ON MICROBIAL,
C     ROOT AND MYCORRHIZAL POPULATIONS IN NITRO.F AND UPTAKE.F
C
C     ROXYX=O2 demand by all microbial,root,mycorrhizal populations
C     RNH4X=NH4 demand in non-band by all microbial,root,mycorrhizal
C        populations
C     RNO3X=NO3 demand in non-band by all microbial,root,mycorrhizal
C        populations
C     RPO4X=H2PO4 demand in non-band by all microbial,root,mycorrhizal
C        populations
C     RP14X=HPO4 demand in non-band by all microbial,root,mycorrhizal
C        populations
C     ROXYS=O2 demand from DOC,DOA oxidation by all microbial
C        populations
C     RVMX4=demand for NH4 oxidation by all microbial populations
C     RVMX3=demand for NO3 reduction by all microbial populations
C     RVMX2=demand for NO2 oxidation by all microbial populations
C     RVMX1=demand for N2O reduction by all microbial populations
C     RINHO,RINHOR=substrate-unlimited NH4 immobilization by all
C        microbial,root,mycorrhizal populations
C     RINOO,RINOOR=substrate-unlimited NO3 immobilization by all
C        microbial,root,mycorrhizal populations
C     RIPOO,RIPOOR=substrate-unlimited H2PO4 immobilization by all
C        microbial,root,mycorrhizal populations
C     RIPO1,RIPO1R=substrate-unlimited HPO4 immobilization by all
C        microbial,root,mycorrhizal populations
C     ROQCX,ROQAX=total DOC,DOA demand from DOC,DOA oxidation by all
C        microbial populations
C     ROQCS,ROQAS=DOC,DOA demand from DOC,DOA oxidation by all
C        microbial populations
C     RNO2X=total demand for NO2 reduction by all microbial
C        populations
C     RVMXC=demand for NO2 reduction by all microbial populations
C
      DO 8990 K=0,5
      IF(K.NE.3.AND.K.NE.4)THEN
      DO 8980 N=1,7
      ROXYX(0,NY,NX)=ROXYX(0,NY,NX)+ROXYS(N,K,0,NY,NX)
      RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RVMX4(N,K,0,NY,NX)
      RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RVMX3(N,K,0,NY,NX)
      RNO2X(0,NY,NX)=RNO2X(0,NY,NX)+RVMX2(N,K,0,NY,NX)
      RN2OX(0,NY,NX)=RN2OX(0,NY,NX)+RVMX1(N,K,0,NY,NX)
      RNH4X(0,NY,NX)=RNH4X(0,NY,NX)+RINHO(N,K,0,NY,NX)
      RNO3X(0,NY,NX)=RNO3X(0,NY,NX)+RINOO(N,K,0,NY,NX)
      RPO4X(0,NY,NX)=RPO4X(0,NY,NX)+RIPOO(N,K,0,NY,NX)
      RP14X(0,NY,NX)=RP14X(0,NY,NX)+RIPO1(N,K,0,NY,NX)
      RNH4X(NU(NY,NX),NY,NX)=RNH4X(NU(NY,NX),NY,NX)+RINHOR(N,K,NY,NX)
      RNO3X(NU(NY,NX),NY,NX)=RNO3X(NU(NY,NX),NY,NX)+RINOOR(N,K,NY,NX)
      RPO4X(NU(NY,NX),NY,NX)=RPO4X(NU(NY,NX),NY,NX)+RIPOOR(N,K,NY,NX)
      RP14X(NU(NY,NX),NY,NX)=RP14X(NU(NY,NX),NY,NX)+RIPO1R(N,K,NY,NX)
      IF(K.LE.4)THEN
      ROQCX(K,0,NY,NX)=ROQCX(K,0,NY,NX)+ROQCS(N,K,0,NY,NX)
      ROQAX(K,0,NY,NX)=ROQAX(K,0,NY,NX)+ROQAS(N,K,0,NY,NX)
      ENDIF
8980  CONTINUE
      ENDIF
8990  CONTINUE
      RNO2X(0,NY,NX)=RNO2X(0,NY,NX)+RVMXC(0,NY,NX)
C
C     SINK ALL SOLID C,N,P IN POND
C
C     BKDS=bulk density (0=water,>0=soil)
C     DLYR=layer depth
C     VLS=sinking rate from hour1.f
C     FSINK=sediment sinking fraction
C     XNFH=time step from wthr.f
C     L,LL=pond layer from,to which material sinks
C
      TSEDSK(NY,NX)=0.0
      DO 9885 L=NL(NY,NX)-1,0,-1
      IF((BKDS(L,NY,NX).LE.ZERO
     2.OR.(L.EQ.0.AND.BKDS(NU(NY,NX),NY,NX).LE.ZERO))
     2.AND.DLYR(3,L,NY,NX).GT.ZERO)THEN
      DO 9880 LL=L+1,NL(NY,NX)
      IF(DLYR(3,LL,NY,NX).GT.ZEROS(NY,NX))GO TO 9881
9880  CONTINUE
9881  CONTINUE
      FSINK=AMIN1(1.0,VLS(NY,NX)/DLYR(3,L,NY,NX))*XNFH
C
C     SINK SOIL MINERALS
C
C     SAND,SILT,CLAY=sand,silt,clay content
C     XCEC,XAEC=cation,anion exchange capacity
C
      IF(L.GT.0)THEN
      FSAN=FSINK*SAND(L,NY,NX)
      FSIL=FSINK*SILT(L,NY,NX)
      FCLA=FSINK*CLAY(L,NY,NX)
      FCEC=FSINK*XCEC(L,NY,NX)
      FAEC=FSINK*XAEC(L,NY,NX)
      SAND(L,NY,NX)=SAND(L,NY,NX)-FSAN
      SILT(L,NY,NX)=SILT(L,NY,NX)-FSIL
      CLAY(L,NY,NX)=CLAY(L,NY,NX)-FCLA
      XCEC(L,NY,NX)=XCEC(L,NY,NX)-FCEC
      XAEC(L,NY,NX)=XAEC(L,NY,NX)-FAEC
      SAND(LL,NY,NX)=SAND(LL,NY,NX)+FSAN
      SILT(LL,NY,NX)=SILT(LL,NY,NX)+FSIL
      CLAY(LL,NY,NX)=CLAY(LL,NY,NX)+FCLA
      XCEC(LL,NY,NX)=XCEC(LL,NY,NX)+FCEC
      XAEC(LL,NY,NX)=XAEC(LL,NY,NX)+FAEC
      IF(BKDS(LL,NY,NX).GT.ZERO)THEN
      NUL(NY,NX)=LL
      TSEDSK(NY,NX)=TSEDSK(NY,NX)+FSAN+FSIL+FCLA
      ENDIF
C
C     SINK PRECIPITATES
C
C     sediment code
C       :PALO,PFEO=precipitated AlOH,FeOH
C       :PCAC,PCAS=precipitated CaCO3,CaSO4
C
      FALO=FSINK*PALOH(L,NY,NX)
      FFEO=FSINK*PFEOH(L,NY,NX)
      FCAC=FSINK*PCACO(L,NY,NX)
      FCAS=FSINK*PCASO(L,NY,NX)
      PALOH(L,NY,NX)=PALOH(L,NY,NX)-FALO
      PFEOH(L,NY,NX)=PFEOH(L,NY,NX)-FFEO
      PCACO(L,NY,NX)=PCACO(L,NY,NX)-FCAC
      PCASO(L,NY,NX)=PCASO(L,NY,NX)-FCAS
      PALOH(LL,NY,NX)=PALOH(LL,NY,NX)+FALO
      PFEOH(LL,NY,NX)=PFEOH(LL,NY,NX)+FFEO
      PCACO(LL,NY,NX)=PCACO(LL,NY,NX)+FCAC
      PCASO(LL,NY,NX)=PCASO(LL,NY,NX)+FCAS
      ENDIF
C
C     SINK FERTILIZER POOLS
C
C     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3
C
      FNH4=FSINK*ZNH4FA(L,NY,NX)
      FNH3=FSINK*ZNH3FA(L,NY,NX)
      FNHU=FSINK*ZNHUFA(L,NY,NX)
      FNO3=FSINK*ZNO3FA(L,NY,NX)
      ZNH4FA(L,NY,NX)=ZNH4FA(L,NY,NX)-FNH4
      ZNH3FA(L,NY,NX)=ZNH3FA(L,NY,NX)-FNH3
      ZNHUFA(L,NY,NX)=ZNHUFA(L,NY,NX)-FNHU
      ZNO3FA(L,NY,NX)=ZNO3FA(L,NY,NX)-FNO3
      ZNH4FA(LL,NY,NX)=ZNH4FA(LL,NY,NX)+FNH4
      ZNH3FA(LL,NY,NX)=ZNH3FA(LL,NY,NX)+FNH3
      ZNHUFA(LL,NY,NX)=ZNHUFA(LL,NY,NX)+FNHU
      ZNO3FA(LL,NY,NX)=ZNO3FA(LL,NY,NX)+FNO3
C
C     SINK EXCHANGEABLE CATIONS AND ANIONS
C
C     sediment code
C       :XN4=adsorbed NH4
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
C           =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2
C       :XH1P,XH2P=adsorbed HPO4,H2PO4
C
      FN4=FSINK*XN4(L,NY,NX)
      FOH0=FSINK*XOH0(L,NY,NX)
      FOH1=FSINK*XOH1(L,NY,NX)
      FOH2=FSINK*XOH2(L,NY,NX)
      FH1P=FSINK*XH1P(L,NY,NX)
      FH2P=FSINK*XH2P(L,NY,NX)
      XN4(L,NY,NX)=XN4(L,NY,NX)-FN4
      XOH0(L,NY,NX)=XOH0(L,NY,NX)-FOH0
      XOH1(L,NY,NX)=XOH1(L,NY,NX)-FOH1
      XOH2(L,NY,NX)=XOH2(L,NY,NX)-FOH2
      XH1P(L,NY,NX)=XH1P(L,NY,NX)-FH1P
      XH2P(L,NY,NX)=XH2P(L,NY,NX)-FH2P
      XN4(LL,NY,NX)=XN4(LL,NY,NX)+FN4
      XOH0(LL,NY,NX)=XOH0(LL,NY,NX)+FOH0
      XOH1(LL,NY,NX)=XOH1(LL,NY,NX)+FOH1
      XOH2(LL,NY,NX)=XOH2(LL,NY,NX)+FOH2
      XH1P(LL,NY,NX)=XH1P(LL,NY,NX)+FH1P
      XH2P(LL,NY,NX)=XH2P(LL,NY,NX)+FH2P
      IF(L.GT.0)THEN
      FHY=FSINK*XHY(L,NY,NX)
      FAL=FSINK*XAL(L,NY,NX)
      FFE=FSINK*XFE(L,NY,NX)
      FCA=FSINK*XCA(L,NY,NX)
      FMG=FSINK*XMG(L,NY,NX)
      FNA=FSINK*XNA(L,NY,NX)
      FKA=FSINK*XKA(L,NY,NX)
      FHC=FSINK*XHC(L,NY,NX)
      FAL2=FSINK*XALO2(L,NY,NX)
      FFE2=FSINK*XFEO2(L,NY,NX)
      XHY(L,NY,NX)=XHY(L,NY,NX)-FHY
      XAL(L,NY,NX)=XAL(L,NY,NX)-FAL
      XFE(L,NY,NX)=XFE(L,NY,NX)-FFE
      XCA(L,NY,NX)=XCA(L,NY,NX)-FCA
      XMG(L,NY,NX)=XMG(L,NY,NX)-FMG
      XNA(L,NY,NX)=XNA(L,NY,NX)-FNA
      XKA(L,NY,NX)=XKA(L,NY,NX)-FKA
      XHC(L,NY,NX)=XHC(L,NY,NX)-FHC
      XALO2(L,NY,NX)=XALO2(L,NY,NX)-FAL2
      XFEO2(L,NY,NX)=XFEO2(L,NY,NX)-FFE2
      XHY(LL,NY,NX)=XHY(LL,NY,NX)+FHY
      XAL(LL,NY,NX)=XAL(LL,NY,NX)+FAL
      XFE(LL,NY,NX)=XFE(LL,NY,NX)+FFE
      XCA(LL,NY,NX)=XCA(LL,NY,NX)+FCA
      XMG(LL,NY,NX)=XMG(LL,NY,NX)+FMG
      XNA(LL,NY,NX)=XNA(LL,NY,NX)+FNA
      XKA(LL,NY,NX)=XKA(LL,NY,NX)+FKA
      XHC(LL,NY,NX)=XHC(LL,NY,NX)+FHC
      XALO2(LL,NY,NX)=XALO2(LL,NY,NX)+FAL2
      XFEO2(LL,NY,NX)=XFEO2(LL,NY,NX)+FFE2
      ENDIF
C
C     SINK PHOSPHORUS PRECIPITATES
C
C     sediment code
C       :PALP,PFEP=precipitated AlPO4,FEPO4
C       :PCPM,PCPD,PCPH=precipitated CaH2PO4,CaHPO4,apatite
C
      FALP=FSINK*PALPO(L,NY,NX)
      FFEP=FSINK*PFEPO(L,NY,NX)
      FCPD=FSINK*PCAPD(L,NY,NX)
      FCPH=FSINK*PCAPH(L,NY,NX)
      FCPM=FSINK*PCAPM(L,NY,NX)
      PALPO(L,NY,NX)=PALPO(L,NY,NX)-FALP
      PFEPO(L,NY,NX)=PFEPO(L,NY,NX)-FFEP
      PCAPD(L,NY,NX)=PCAPD(L,NY,NX)-FCPD
      PCAPH(L,NY,NX)=PCAPH(L,NY,NX)-FCPH
      PCAPM(L,NY,NX)=PCAPM(L,NY,NX)-FCPM
      PALPO(LL,NY,NX)=PALPO(LL,NY,NX)+FALP
      PFEPO(LL,NY,NX)=PFEPO(LL,NY,NX)+FFEP
      PCAPD(LL,NY,NX)=PCAPD(LL,NY,NX)+FCPD
      PCAPH(LL,NY,NX)=PCAPH(LL,NY,NX)+FCPH
      PCAPM(LL,NY,NX)=PCAPM(LL,NY,NX)+FCPM
C
C     SINK MICROBIAL C,N,P
C
C     OMC,OMN,OMP=microbial C,N,P
C
      DO 1970 K=0,5
      DO 1960 N=1,7
      DO 1960 M=1,3
      FOMC=FSINK*OMC(M,N,K,L,NY,NX)
      FOMN=FSINK*OMN(M,N,K,L,NY,NX)
      FOMP=FSINK*OMP(M,N,K,L,NY,NX)
      OMC(M,N,K,LL,NY,NX)=OMC(M,N,K,LL,NY,NX)+FOMC
      OMN(M,N,K,LL,NY,NX)=OMN(M,N,K,LL,NY,NX)+FOMN
      OMP(M,N,K,LL,NY,NX)=OMP(M,N,K,LL,NY,NX)+FOMP
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)-FOMC
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)-FOMN
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)-FOMP
C     IF(BKDS(LL,NY,NX).GT.ZERO)THEN
C     TSEDSK(NY,NX)=TSEDSK(NY,NX)+FOMC*1.82E-06
C     ENDIF
1960  CONTINUE
1970  CONTINUE
C
C     SINK MICROBIAL RESIDUE C,N,P
C
C     ORC,ORN,ORP=microbial residue C,N,P
C
      DO 1900 K=0,4
      DO 1940 M=1,2
      FORC=FSINK*ORC(M,K,L,NY,NX)
      FORN=FSINK*ORN(M,K,L,NY,NX)
      FORP=FSINK*ORP(M,K,L,NY,NX)
      ORC(M,K,LL,NY,NX)=ORC(M,K,LL,NY,NX)+FORC
      ORN(M,K,LL,NY,NX)=ORN(M,K,LL,NY,NX)+FORN
      ORP(M,K,LL,NY,NX)=ORP(M,K,LL,NY,NX)+FORP
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-FORC
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-FORN
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-FORP
C     IF(BKDS(LL,NY,NX).GT.ZERO)THEN
C     TSEDSK(NY,NX)=TSEDSK(NY,NX)+FORC*1.82E-06
C     ENDIF
1940  CONTINUE
C
C     SINK ADSORBED C,N,P
C
C     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
C
      FOHC=FSINK*OHC(K,L,NY,NX)
      FOHN=FSINK*OHN(K,L,NY,NX)
      FOHP=FSINK*OHP(K,L,NY,NX)
      FOHA=FSINK*OHA(K,L,NY,NX)
      OHC(K,LL,NY,NX)=OHC(K,LL,NY,NX)+FOHC
      OHN(K,LL,NY,NX)=OHN(K,LL,NY,NX)+FOHN
      OHP(K,LL,NY,NX)=OHP(K,LL,NY,NX)+FOHP
      OHA(K,LL,NY,NX)=OHA(K,LL,NY,NX)+FOHA
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-FOHC
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-FOHN
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-FOHP
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-FOHA
C     IF(BKDS(LL,NY,NX).GT.ZERO)THEN
C     TSEDSK(NY,NX)=TSEDSK(NY,NX)+FOHC*1.82E-06
C     ENDIF
C
C     SINK SOC,SON,SOP
C
C     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
C
      DO 1930 M=1,5
      FOSC=FSINK*OSC(M,K,L,NY,NX)
      FOSA=FSINK*OSA(M,K,L,NY,NX)
      FOSN=FSINK*OSN(M,K,L,NY,NX)
      FOSP=FSINK*OSP(M,K,L,NY,NX)
      OSC(M,K,LL,NY,NX)=OSC(M,K,LL,NY,NX)+FOSC
      OSA(M,K,LL,NY,NX)=OSA(M,K,LL,NY,NX)+FOSA
      OSN(M,K,LL,NY,NX)=OSN(M,K,LL,NY,NX)+FOSN
      OSP(M,K,LL,NY,NX)=OSP(M,K,LL,NY,NX)+FOSP
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-FOSC
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-FOSA
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-FOSN
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-FOSP
C     IF(BKDS(LL,NY,NX).GT.ZERO)THEN
C     TSEDSK(NY,NX)=TSEDSK(NY,NX)+FOSC*1.82E-06
C     ENDIF
1930  CONTINUE
1900  CONTINUE
      ENDIF
9885  CONTINUE
C
C     RUNOFF AND SUBSURFACE BOUNDARY FLUXES
C
      DO 9985 L=NU(NY,NX),NL(NY,NX)
C
C     LOCATE EXTERNAL BOUNDARIES
C
C     N2,N1=NY,NX of source grid cell
C     N5,N4=NY,NX of destination grid cell
C
      N1=NX
      N2=NY
      DO 9980 N=1,3
      DO 9975 NN=1,2
      IF(N.EQ.1)THEN
      IF(NN.EQ.1)THEN
      IF(NX.EQ.NHE)THEN
      N4=NX+1
      N5=NY
      N6=L
      XN=-1.0
      ELSE
      GO TO 9975
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NX.EQ.NHW)THEN
      N4=NX
      N5=NY
      N6=L
      XN=1.0
      ELSE
      GO TO 9975
      ENDIF
      ENDIF
      ELSEIF(N.EQ.2)THEN
      IF(NN.EQ.1)THEN
      IF(NY.EQ.NVS)THEN
      N4=NX
      N5=NY+1
      N6=L
      XN=-1.0
      ELSE
      GO TO 9975
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      IF(NY.EQ.NVN)THEN
      N4=NX
      N5=NY
      N6=L
      XN=1.0
      ELSE
      GO TO 9975
      ENDIF
      ENDIF
      ELSEIF(N.EQ.3)THEN
      IF(NN.EQ.1)THEN
      IF(L.EQ.NL(NY,NX))THEN
      N4=NX
      N5=NY
      N6=L+1
      XN=-1.0
      ELSE
      GO TO 9975
      ENDIF
      ELSEIF(NN.EQ.2)THEN
      GO TO 9975
      ENDIF
      ENDIF
C
C     RUNOFF BOUNDARY FLUXES OF WATER AND HEAT
C
C     QR,QS,QW,QI=runoff from surface water, snowpack
C        snow,water,ice from watsub.f
C     CRUN,URUN=cumulative water+snow runoff
C     HEATOU=cumulative heat loss through lateral and lower boundaries
C
      IF(N.NE.3.AND.L.EQ.NU(NY,NX))THEN
      WQRN=XN*QR(N,NN,N5,N4)
      WQRH(N2,N1)=WQRH(N2,N1)+WQRN
      IF(ABS(WQRN).GT.ZEROS(N5,N4))THEN
      CRUN=CRUN-WQRN
      URUN(NY,NX)=URUN(NY,NX)-WQRN
      HQRN=XN*HQR(N,NN,N5,N4)
      HEATOU=HEATOU-HQRN
C     WRITE(*,2229)'WQRN',I,J,NFZ,N,NN,N5,N4
C    2,HEATOU,HQRN,HQR(N,NN,N5,N4),CRUN,WQRN,QR(N,NN,N5,N4)
2229  FORMAT(A8,7I4,20F20.6)
C
C     RUNOFF BOUNDARY FLUXES OF C, N AND P
C
C     X*QRS(W),X*QSS=solute flux in runoff, snow drift from trnsfr.f
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,P1=HPO4,P4=H2PO4
C     XN=direction indicator
C     TCOU,TZOU,TPOU=cumulative C,N,P loss through lateral and
C        lower boundaries
C     UDOCQ,UDICQ=dissolved organic,inorganic C loss through lateral
C        and lower boundaries
C     UDONQ,UDINQ=dissolved organic,inorganic N loss through lateral
C        and lower boundaries
C     UDOPQ,UDIPQ=dissolved organic,inorganic P loss through lateral
C        and lower boundaries
C
      CXR=XN*(XCOQRS(N,NN,N5,N4)+XCHQRS(N,NN,N5,N4))
      ZXR=XN*(XN4QRW(N,NN,N5,N4)+XN3QRW(N,NN,N5,N4)
     2+XNOQRW(N,NN,N5,N4)+XNXQRS(N,NN,N5,N4))
      ZGR=XN*(XN2QRS(N,NN,N5,N4)+XNGQRS(N,NN,N5,N4))
      PXR=XN*(XP4QRW(N,NN,N5,N4)+XP1QRW(N,NN,N5,N4))
      COR=0.0
      ZOR=0.0
      POR=0.0
      DO 2575 K=0,4
      COR=COR+XN*(XOCQRS(K,N,NN,N5,N4)+XOAQRS(K,N,NN,N5,N4))
      ZOR=ZOR+XN*XONQRS(K,N,NN,N5,N4)
      POR=POR+XN*XOPQRS(K,N,NN,N5,N4)
2575  CONTINUE
      TCOU=TCOU-CXR-COR
      TZOU=TZOU-ZXR-ZOR-ZGR
      TPOU=TPOU-PXR-POR
      UDOCQ(NY,NX)=UDOCQ(NY,NX)-COR
      UDICQ(NY,NX)=UDICQ(NY,NX)-CXR
      UDONQ(NY,NX)=UDONQ(NY,NX)-ZOR
      UDINQ(NY,NX)=UDINQ(NY,NX)-ZXR
      UDOPQ(NY,NX)=UDOPQ(NY,NX)-POR
      UDIPQ(NY,NX)=UDIPQ(NY,NX)-PXR
      OXR=XN*XOXQRS(N,NN,N5,N4)
      OXYGOU=OXYGOU-OXR
      HGR=XN*XHGQRS(N,NN,N5,N4)
      H2GOU=H2GOU+HGR
C     WRITE(*,6636)'XPO',I,J,NFZ,N4,N5,N,XN
C    3,UDICQ(NY,NX),UDICQ(NY,NX),UDICQ(NY,NX)
C    2,COR,CXR,ZOR,ZXR,ZGR,POR,PXR
C    2,XP4QRW(N,NN,N5,N4),XP1QRW(N,NN,N5,N4)
6636  FORMAT(A8,6I4,40E12.4)
C
C     RUNOFF BOUNDARY FLUXES OF SOLUTES
C
C     ISALTG=salt flag
C     XQR*,XQS*=solute loss in runoff,snow drift from trnsfrs.f
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C        :*1=non-band,*B=band
C     XN=direction indicator
C     TPOU,TIONOU=accumulated P,salt loss through lateral and lower
C         boundaries
C
      IF(ISALTG.NE.0)THEN
      PSS=XN*31.0*(XQRH0P(N,NN,N5,N4)
     2+XQRC0P(N,NN,N5,N4)+XQRF1P(N,NN,N5,N4)+XQRC1P(N,NN,N5,N4)
     3+XQRM1P(N,NN,N5,N4)+XQRH3P(N,NN,N5,N4)+XQRF2P(N,NN,N5,N4)
     4+XQRC2P(N,NN,N5,N4))
      SS1=XN*(XQRAL(N,NN,N5,N4)+XQRFE(N,NN,N5,N4)+XQRHY(N,NN,N5,N4)
     2+XQRCA(N,NN,N5,N4)+XQRMG(N,NN,N5,N4)+XQRNA(N,NN,N5,N4)
     3+XQRKA(N,NN,N5,N4)+XQROH(N,NN,N5,N4)+XQRSO(N,NN,N5,N4)
     4+XQRCL(N,NN,N5,N4)+XQRC3(N,NN,N5,N4)+XQRH0P(N,NN,N5,N4))
      SS2=XN*2.0*(XQRHC(N,NN,N5,N4)+XQRAL1(N,NN,N5,N4)
     2+XQRALS(N,NN,N5,N4)+XQRFE1(N,NN,N5,N4)+XQRFES(N,NN,N5,N4)
     3+XQRCAO(N,NN,N5,N4)+XQRCAC(N,NN,N5,N4)+XQRCAS(N,NN,N5,N4)
     4+XQRMGO(N,NN,N5,N4)+XQRMGC(N,NN,N5,N4)+XQRMGS(N,NN,N5,N4)
     5+XQRNAC(N,NN,N5,N4)+XQRNAS(N,NN,N5,N4)+XQRKAS(N,NN,N5,N4)
     6+XQRC0P(N,NN,N5,N4))
      SS3=XN*3.0*(XQRAL2(N,NN,N5,N4)+XQRFE2(N,NN,N5,N4)
     2+XQRCAH(N,NN,N5,N4)+XQRMGH(N,NN,N5,N4)+XQRF1P(N,NN,N5,N4)
     3+XQRC1P(N,NN,N5,N4)+XQRM1P(N,NN,N5,N4))
      SS4=XN*4.0*(XQRAL3(N,NN,N5,N4)+XQRFE3(N,NN,N5,N4)
     2+XQRH3P(N,NN,N5,N4)+XQRF2P(N,NN,N5,N4)+XQRC2P(N,NN,N5,N4))
     5+XN*5.0*(XQRAL4(N,NN,N5,N4)+XQRFE4(N,NN,N5,N4))
      PSS=PSS+XN*31.0*XQSH0P(N,NN,N5,N4)
      TPOU=TPOU-PSS
      SSR=SS1+SS2+SS3+SS4
      TIONOU=TIONOU-SSR
      UIONOU(NY,NX)=UIONOU(NY,NX)-SSR
C     WRITE(20,3336)'SSR',I,J,N,N5,N4,SSR,SS1,SS2,SS3,SS4,TIONOU
3336  FORMAT(A8,5I6,20F16.9)
C
C     SURFACE RUNOFF ELECTRICAL CONDUCTIVITY
C
C     QR=surface runoff from watsub.f
C     XQR*=solute flux in runoff from trnsfrs.f
C     EC*=electrical conductivity of ion
C     ECNDQ=total electrical conductivity of runoff
C     ion code:HY=H,OH=OH,AL=Al,FE=Fe,CA=Ca,MG=Mg,NA=Na,KA=K
C             :CO=CO3,HC=HCO3,SO=SO4,CL=Cl
C
      WX=QR(N,NN,N5,N4)
      IF(ABS(WX).GT.ZEROS(N5,N4))THEN
      ECHY=0.337*AMAX1(0.0,XQRHY(N,NN,N5,N4)/WX)
      ECOH=0.192*AMAX1(0.0,XQROH(N,NN,N5,N4)/WX)
      ECAL=0.056*AMAX1(0.0,XQRAL(N,NN,N5,N4)*3.0/WX)
      ECFE=0.051*AMAX1(0.0,XQRFE(N,NN,N5,N4)*3.0/WX)
      ECCA=0.060*AMAX1(0.0,XQRCA(N,NN,N5,N4)*2.0/WX)
      ECMG=0.053*AMAX1(0.0,XQRMG(N,NN,N5,N4)*2.0/WX)
      ECNA=0.050*AMAX1(0.0,XQRNA(N,NN,N5,N4)/WX)
      ECKA=0.070*AMAX1(0.0,XQRKA(N,NN,N5,N4)/WX)
      ECCO=0.072*AMAX1(0.0,XQRC3(N,NN,N5,N4)*2.0/WX)
      ECHC=0.044*AMAX1(0.0,XQRHC(N,NN,N5,N4)/WX)
      ECSO=0.080*AMAX1(0.0,XQRSO(N,NN,N5,N4)*2.0/WX)
      ECCL=0.076*AMAX1(0.0,XQRCL(N,NN,N5,N4)/WX)
      ECNO=0.071*AMAX1(0.0,XNOQRW(N,NN,N5,N4)/(WX*14.0))
      ECNDQ=ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA
     2+ECCO+ECHC+ECSO+ECCL+ECNO
C     WRITE(*,9991)'ECNDQ',IYRC,I,J,N4,N5,N,NN,WX,ECNDQ
9991  FORMAT(A8,7I4,2E12.4)
      ELSE
      ECNDQ=0.0
      ENDIF
      ENDIF
      ENDIF
C
C     WATER, HEAT, C,N,P, SOLUTES IN SNOW DRIFT
C
C     QS,QW,QI=aggregated snow,water,ice transfer from watsub.f
C     CRUN,URUN=cumulative water+snow runoff
C     HQS=aggregated convective heat from snow,water,ice transfer
C        from watsub.f�
C     HEATOU=cumulative heat loss through lateral and lower boundaries
C     X*QSS=aggregated solute in snow flux from trnsfr.f
C     solute code:CO=CO2,CH=CH4,N4=NH4,N3=NH3,NO=NO3,N2=N2O,NG=N2
C                :P4=H2PO4,P1=HPO4
C     XN=direction indicator
C     TCOU,TZOU,TPOU,OXYGOU=cumulative C,N,P,O2 loss through lateral
C        and lower boundaries
C     UDICQ=dissolved inorganic C loss through lateral
C        and lower boundaries
C     UDINQ=dissolved inorganic N loss through lateral
C        and lower boundaries
C     UDIPQ=dissolved inorganic P loss through lateral
C        and lower boundaries
C
      WQRS=XN*(QS(N,NN,N5,N4)+QW(N,NN,N5,N4)+QI(N,NN,N5,N4)*DENSI)
      IF(ABS(WQRS).GT.ZEROS(N5,N4))THEN
      CRUN=CRUN-WQRS
      URUN(NY,NX)=URUN(NY,NX)-WQRS
      HQRS=XN*HQS(N,NN,N5,N4)
      HEATOU=HEATOU-HQRS
C     WRITE(*,2229)'HQRS',I,J,NFZ,N,NN,N5,N4
C    2,HEATOU,HQRS,HQS(N,NN,N5,N4),CRUN,WQRS,QS(N,NN,N5,N4)
C    3,QW(N,NN,N5,N4),QI(N,NN,N5,N4)
      CXS=XN*(XCOQSS(N,NN,N5,N4)+XCHQSS(N,NN,N5,N4))
      ZXS=XN*(XN4QSS(N,NN,N5,N4)+XN3QSS(N,NN,N5,N4)
     2+XNOQSS(N,NN,N5,N4))
      ZGS=XN*(XN2QSS(N,NN,N5,N4)+XNGQSS(N,NN,N5,N4))
      PXS=XN*(XP4QSS(N,NN,N5,N4)+XP1QSS(N,NN,N5,N4))
      TCOU=TCOU-CXS
      TZOU=TZOU-ZXS-ZGS
      TPOU=TPOU-PXS
      UDICQ(NY,NX)=UDICQ(NY,NX)-CXR
      UDINQ(NY,NX)=UDINQ(NY,NX)-ZXR-ZGR
      UDIPQ(NY,NX)=UDIPQ(NY,NX)-PXR
      OXS=XN*XOXQSS(N,NN,N5,N4)
      OXYGOU=OXYGOU-OXS
C
C     SALT SOLUTES IN SNOW DRIFT
C
C     QSS*=aggregated solute in snow flux from trnsfrs.f
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C     TPOU,TIONOU=accumulated P,salt loss through lateral and lower
C         boundaries
C
      IF(ISALTG.NE.0)THEN
      PSS=XN*31.0*(XQSC0P(N,NN,N5,N4)+XQSF1P(N,NN,N5,N4)
     2+XQSC1P(N,NN,N5,N4)+XQSM1P(N,NN,N5,N4)+XQSH3P(N,NN,N5,N4)
     3+XQSF2P(N,NN,N5,N4)+XQSC2P(N,NN,N5,N4))
      TPOU=TPOU-PSS
      SS1=XN*(XQSAL(N,NN,N5,N4)+XQSFE(N,NN,N5,N4)
     2+XQSHY(N,NN,N5,N4)+XQSCA(N,NN,N5,N4)+XQSMG(N,NN,N5,N4)
     3+XQSNA(N,NN,N5,N4) +XQSKA(N,NN,N5,N4)+XQSOH(N,NN,N5,N4)
     4+XQSSO(N,NN,N5,N4)+XQSCL(N,NN,N5,N4)+XQSC3(N,NN,N5,N4)
     5+XQSH0P(N,NN,N5,N4))
      SS2=XN*2.0*(XQSHC(N,NN,N5,N4)+XQSAL1(N,NN,N5,N4)
     2+XQSALS(N,NN,N5,N4)+XQSFE1(N,NN,N5,N4)+XQSFES(N,NN,N5,N4)
     3+XQSCAO(N,NN,N5,N4)+XQSCAC(N,NN,N5,N4)+XQSCAS(N,NN,N5,N4)
     4+XQSMGO(N,NN,N5,N4)+XQSMGC(N,NN,N5,N4)+XQSMGS(N,NN,N5,N4)
     5+XQSNAC(N,NN,N5,N4)+XQSNAS(N,NN,N5,N4)+XQSKAS(N,NN,N5,N4)
     6+XQSC0P(N,NN,N5,N4))
      SS3=XN*3.0*(XQSAL2(N,NN,N5,N4)+XQSFE2(N,NN,N5,N4)
     2+XQSCAH(N,NN,N5,N4)+XQSMGH(N,NN,N5,N4)+XQSF1P(N,NN,N5,N4)
     3+XQSC1P(N,NN,N5,N4)+XQSM1P(N,NN,N5,N4))
      SS4=XN*4.0*(XQSAL3(N,NN,N5,N4)+XQSFE3(N,NN,N5,N4)
     2+XQSH3P(N,NN,N5,N4)+XQSF2P(N,NN,N5,N4)+XQSC2P(N,NN,N5,N4))
     3+XN*5.0*(XQSAL4(N,NN,N5,N4)+XQSFE4(N,NN,N5,N4))
      SSR=SS1+SS2+SS3+SS4
      TIONOU=TIONOU-SSR
      UIONOU(NY,NX)=UIONOU(NY,NX)-SSR
      ENDIF
      ENDIF
C
C     RUNOFF BOUNDARY FLUXES OF SEDIMENT FROM EROSION
C
C     IERSNG=erosion flag from site file
C     XSEDER=total sediment flux from erosion.f
C     TSEDOU,USEDOU=cumulative sediment loss through lateral and
C        lower boundaries
C
      IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
      IF(ABS(XSEDER(N,NN,N5,N4)).GT.ZEROS(N5,N4))THEN
      ER=XN*XSEDER(N,NN,N5,N4)
      TSEDOU=TSEDOU-ER
      USEDOU(NY,NX)=USEDOU(NY,NX)-ER
C
C     EXCHANGEABLE AND PRECIPITATED N,P IN RUNOFF SEDIMENT
C
C     *ER=sediment flux from erosion.f
C     sediment code
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XNH4,XNH3,XNHU,XNO3=NH4,NH3,urea,NO3 in non-band
C       :XNH4B,XNH3B,XNHUB,XNO3B=fertilizer NH4,NH3,urea,NO3 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C       :PALO,PFEO=precipitated AlOH,FeOH
C       :PCAC,PCAS=precipitated CaCO3,CaSO4
C       :PALP,PFEP=precipitated AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precipitated AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precipitated CaH2PO4,CaHPO4,apatite
C          in non-band
C       :PCPMB,PCPDB,PCPHB= precipitated CaH2PO4,CaHPO4,apatite
C          in band
C
      CXE=0.0
      ZXE=XN*14.0*(XN4ER(N,NN,N5,N4)+XNBER(N,NN,N5,N4))
      ZPE=XN*14.0*(XNH4ER(N,NN,N5,N4)+XNH3ER(N,NN,N5,N4)
     2+XNHUER(N,NN,N5,N4)+XNO3ER(N,NN,N5,N4)+XNH4EB(N,NN,N5,N4)
     3+XNH3EB(N,NN,N5,N4)+XNHUEB(N,NN,N5,N4)+XNO3EB(N,NN,N5,N4))
      PXE=XN*31.0*(XH1PER(N,NN,N5,N4)+XH2PER(N,NN,N5,N4)
     4+XH1PEB(N,NN,N5,N4)+XH2PEB(N,NN,N5,N4))
      PPE=XN*(31.0*(PALPER(N,NN,N5,N4)+PFEPER(N,NN,N5,N4)
     2+PCPDER(N,NN,N5,N4)+PALPEB(N,NN,N5,N4)
     3+PFEPEB(N,NN,N5,N4)+PCPDEB(N,NN,N5,N4))
     7+62.0*(PCPMER(N,NN,N5,N4)+PCPMEB(N,NN,N5,N4))
     8+93.0*(PCPHER(N,NN,N5,N4)+PCPHEB(N,NN,N5,N4)))
C
C     MICROBIAL RESIDUE C IN RUNOFF SEDIMENT
C
C     *ER=sediment flux from erosion.f
C     sediment code:OMC,OMN,OMP=microbial C,N,P
C
      COE=0.0
      ZOE=0.0
      POE=0.0
      DO 3580 K=0,5
      DO 3580 NO=1,7
      DO 3580 M=1,3
      COE=COE+XN*OMCER(M,NO,K,N,NN,N5,N4)
      ZOE=ZOE+XN*OMNER(M,NO,K,N,NN,N5,N4)
      POE=POE+XN*OMPER(M,NO,K,N,NN,N5,N4)
3580  CONTINUE
C
C     MICROBIAL RESIDUE C IN RUNOFF SEDIMENT
C
C     *ER=sediment flux from erosion.f
C     sediment code:ORC,ORN,ORP=microbial residue C,N,P
C
      DO 3575 K=0,4
      DO 3570 M=1,2
      COE=COE+XN*ORCER(M,K,N,NN,N5,N4)
      ZOE=ZOE+XN*ORNER(M,K,N,NN,N5,N4)
      POE=POE+XN*ORPER(M,K,N,NN,N5,N4)
3570  CONTINUE
C
C     DOC, ADSORBED AND HUMUS C IN RUNOFF SEDIMENT
C
C     *ER=sediment flux from erosion.f
C     sediment code:OHC,OHA,OHN,OHP=adsorbed C,acetate,N,P
C                  :OSC,OSA,OSN,OSP=SOC,acetate,SON,SOP
C
      COE=COE+XN*(OHCER(K,N,NN,N5,N4)+OHAER(K,N,NN,N5,N4))
      ZOE=ZOE+XN*OHNER(K,N,NN,N5,N4)
      POE=POE+XN*OHPER(K,N,NN,N5,N4)
      DO 3565 M=1,5
      COE=COE+XN*OSCER(M,K,N,NN,N5,N4)
      ZOE=ZOE+XN*OSNER(M,K,N,NN,N5,N4)
      POE=POE+XN*OSPER(M,K,N,NN,N5,N4)
3565  CONTINUE
3575  CONTINUE
C
C     TOTAL EROSION C,N,P LOSSES
C
C     TCOU,TZOU,TPOU=cumulative C,N,P loss through lateral and
C        lower boundaries
C     UDOCQ,UDICQ=cumulative dissolved organic,inorganic C loss
C        through lateral and lower boundaries
C     UDONQ,UDINQ=cumulative dissolved organic,inorganic N loss
C        through lateral and lower boundaries
C     UDOPQ,UDIPQ=cumulative dissolved organic,inorganic P loss
C        through lateral and lower boundaries
C     COE,CXE=organic,exchangeable C
C     ZOE,ZXE,ZPE=organic,exchangeable,fertilizer N
C     POE,PXE,PPE=organic,exchangeable,precipitated N
C
      TCOU=TCOU-COE-CXE
      TZOU=TZOU-ZOE-ZXE-ZPE
      TPOU=TPOU-POE-PXE-PPE
      UDOCQ(NY,NX)=UDOCQ(NY,NX)-COE
      UDICQ(NY,NX)=UDICQ(NY,NX)-CXE
      UDONQ(NY,NX)=UDONQ(NY,NX)-ZOE
      UDINQ(NY,NX)=UDINQ(NY,NX)-ZXE-ZPE
      UDOPQ(NY,NX)=UDOPQ(NY,NX)-POE
      UDIPQ(NY,NX)=UDIPQ(NY,NX)-PXE-PPE
C     WRITE(*,6635)'POE',I,J,N4,N5,N,NN
C    2,COE,CXE,ZOE,ZXE,ZPE
C    3,POE,PXE,PPE,TPOU,XSEDER(N,NN,N5,N4)
C    3,XN,TCOU,TZOU,TPOU
6635  FORMAT(A8,6I4,20F17.8)
C
C     ADSORBED AND PRECIPITATED SALTS IN RUNOFF SEDIMENTS
C
C     ISALTG=salt flag
C     *ER=sediment flux from erosion.f
C     sediment code
C       :NH4,NH3,NHU,NO3=fertilizer NH4,NH3,urea,NO3 in non-band
C       :NH4B,NH3B,NHUB,NO3B=fertilizer NH4,NH3,urea,NO3 in band
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
C           =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
C       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C       :PALO,PFEO=precip AlOH,FeOH
C       :PCAC,PCAS=precip CaCO3,CaSO4
C       :PALP,PFEP=precip AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precip AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite in non-band
C       :PCPMB,PCPDB,PCPHB=precip CaH2PO4,CaHPO4,apatite in band
C     TIONOU,UIONOU=total salt loss through lateral and lower
C        boundaries
C
      IF(ISALTG.NE.0)THEN
      SEF=XN*(XNH3ER(N,NN,N5,N4)+XNHUER(N,NN,N5,N4)+XNO3ER(N,NN,N5,N4)
     5+XNH3EB(N,NN,N5,N4)+XNHUEB(N,NN,N5,N4)+XNO3EB(N,NN,N5,N4))
     2+2.0*(XNH4ER(N,NN,N5,N4)+XNH4EB(N,NN,N5,N4))
      SEX=XN*(XHYER(N,NN,N5,N4)+XALER(N,NN,N5,N4)
     2+XFEER(N,NN,N5,N4)+XCAER(N,NN,N5,N4)+XMGER(N,NN,N5,N4)
     3+XNAER(N,NN,N5,N4)+XKAER(N,NN,N5,N4)+XHCER(N,NN,N5,N4)
     4+XOH0ER(N,NN,N5,N4)+XOH0EB(N,NN,N5,N4))
     5+XN*2.0*(XN4ER(N,NN,N5,N4)+XNBER(N,NN,N5,N4)
     6+XOH1ER(N,NN,N5,N4)+XOH1EB(N,NN,N5,N4))
     7+XN*3.0*(XAL2ER(N,NN,N5,N4)+XFE2ER(N,NN,N5,N4)
     8+XOH2ER(N,NN,N5,N4)+XOH2EB(N,NN,N5,N4)
     9+XH1PER(N,NN,N5,N4)+XH1PEB(N,NN,N5,N4))
     1+XN*4.0*(XH2PER(N,NN,N5,N4)+XH2PEB(N,NN,N5,N4))
      SEP=XN*2.0*(PCACER(N,NN,N5,N4)+PCASER(N,NN,N5,N4)
     2+PALPER(N,NN,N5,N4)+PFEPER(N,NN,N5,N4)
     3+PALPEB(N,NN,N5,N4)+PFEPEB(N,NN,N5,N4))
     4+XN*3.0*(PCPDER(N,NN,N5,N4)+PCPDEB(N,NN,N5,N4))
     5+XN*4.0*(PALOER(N,NN,N5,N4)+PFEOER(N,NN,N5,N4))
     6+XN*7.0*(PCPMER(N,NN,N5,N4)+PCPMEB(N,NN,N5,N4))
     7+XN*9.0*(PCPHER(N,NN,N5,N4)+PCPHEB(N,NN,N5,N4))
      SET=SEF+SEX+SEP
      TIONOU=TIONOU-SET
      UIONOU(NY,NX)=UIONOU(NY,NX)-SET
C     WRITE(*,3342)'SET',I,J,N4,N5,NN,N,SET,SEF,SEX,SEP,TIONOU
3342  FORMAT(A8,6I4,12F18.6)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
C     SUBSURFACE BOUNDARY FLUXES OF WATER AND HEAT
C
C     FLW,FLWH,HFLW=micropore,macropore water,heat flux through
C        lateral and lower boundaries from watsub.f
C     FLWY,FLWHY=micropore,macropore discharge through lateral and
C        lower boundaries from artificial drainage from watsub.f
C     VOLWOU,HEATOU=cumulative water, heat transfer through lateral
C        and lower boundaries
C     UVOLO,HVOLO=water transfer through lateral and lower boundaries
C     UVOLY=water loss through lateral and lower boundaries from
C        artificial drainage
C
      IF(NCN(NY,NX).NE.3.OR.N.EQ.3)THEN
      HEATOU=HEATOU-XN*HFLW(N,N6,N5,N4)
      WO=XN*(FLW(N,N6,N5,N4)+FLWH(N,N6,N5,N4))
      IF(WO.NE.0.0)THEN
      VOLWOU=VOLWOU-WO
      HVOLO(N2,N1)=HVOLO(N2,N1)-WO
      UVOLO(N2,N1)=UVOLO(N2,N1)-WO
      WY=XN*(FLWY(N,N6,N5,N4)+FLWHY(N,N6,N5,N4))
      UVOLY(N2,N1)=UVOLY(N2,N1)-WY
C     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,3488)'UVOLO',I,J,NFZ,N6,N5,N4,N,XN,WO,WY
C    2,UVOLO(NY,NX),FLW(N,N6,N5,N4),FLWH(N,N6,N5,N4)
C    2,UVOLY(NY,NX),FLWY(N,N6,N5,N4),FLWHY(N,N6,N5,N4)
3488  FORMAT(A8,7I4,12E12.4)
C     ENDIF
C
C     SUBSURFACE BOUNDARY FLUXES OF C,N,P
C
C     X*FLS,X*FHS=solute flux in macropores,micropores from trnsfr.f
C     X*FLG=convective+diffusive gas flux from trnsfr.f
C     TCOU=cumulative C loss through lateral and lower boundaries
C     UDOCD,UDICD=dissolved organic,inorganic C loss through
C        subsurface boundaries
C     TZOU=cumulative N loss through lateral and lower boundaries
C     UDOND,UDIND=dissolved organic,inorganic N loss through
C        subsurface boundaries
C     TPOU=cumulative P loss through lateral and lower boundaries
C     UDOPD,UDIPD=dissolved organic,inorganic P loss through
C        subsurface boundaries
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,H1P=HPO4,H2P=H2PO4 in non-band
C        :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,H1B=HPO4,H2B=H2PO4 in band
C     TCOU,TZOU,TPOU=cumulative C,N,P loss through lateral and
C        lower boundaries
C     UDOCQ,UDICQ=cumulative dissolved organic,inorganic C loss
C        through lateral and lower boundaries
C     UDONQ,UDINQ=cumulative dissolved organic,inorganic N loss
C        through lateral and lower boundaries
C     UDOPQ,UDIPQ=cumulative dissolved organic,inorganic P loss
C        through lateral and lower boundaries
C
      COD=0.0
      ZOD=0.0
      POD=0.0
      DO 450 K=0,4
      COD=COD+XN*(XOCFLS(K,N,N6,N5,N4)+XOAFLS(K,N,N6,N5,N4)
     4+XOCFHS(K,N,N6,N5,N4)+XOAFHS(K,N,N6,N5,N4))
      ZOD=ZOD+XN*(XONFLS(K,N,N6,N5,N4)+XONFHS(K,N,N6,N5,N4))
      POD=POD+XN*(XOPFLS(K,N,N6,N5,N4)+XOPFHS(K,N,N6,N5,N4))
450   CONTINUE
      CXD=XN*(XCOFLS(N,N6,N5,N4)+XCOFHS(N,N6,N5,N4)
     2+XCOFLG(N,N6,N5,N4)+XCHFLS(N,N6,N5,N4)
     3+XCHFHS(N,N6,N5,N4)+XCHFLG(N,N6,N5,N4))
      ZXD=XN*(XN4FLW(N,N6,N5,N4)+XN3FLW(N,N6,N5,N4)+XNOFLW(N,N6,N5,N4)
     2+XN4FLB(N,N6,N5,N4)+XN3FLB(N,N6,N5,N4)+XNOFLB(N,N6,N5,N4)
     3+XNXFLS(N,N6,N5,N4)+XNXFLB(N,N6,N5,N4)
     5+XN4FHW(N,N6,N5,N4)+XN3FHW(N,N6,N5,N4)+XNOFHW(N,N6,N5,N4)
     6+XN4FHB(N,N6,N5,N4)+XN3FHB(N,N6,N5,N4)+XNOFHB(N,N6,N5,N4)
     7+XNXFHS(N,N6,N5,N4)+XNXFHB(N,N6,N5,N4))
      ZGD=XN*(XNGFLS(N,N6,N5,N4)+XNGFLG(N,N6,N5,N4)+XNGFHS(N,N6,N5,N4)
     2+XN2FLS(N,N6,N5,N4)+XN2FLG(N,N6,N5,N4)+XN2FHS(N,N6,N5,N4)
     3+XN3FLG(N,N6,N5,N4))
      PXD=XN*(XH2PFS(N,N6,N5,N4)+XH2BFB(N,N6,N5,N4)
     2+XH2PHS(N,N6,N5,N4)+XH2BHB(N,N6,N5,N4)+XH1PFS(N,N6,N5,N4)
     3+XH1BFB(N,N6,N5,N4)+XH1PHS(N,N6,N5,N4)+XH1BHB(N,N6,N5,N4))
      TCOU=TCOU-COD-CXD
      TZOU=TZOU-ZOD-ZXD-ZGD
      TPOU=TPOU-POD-PXD
      UDOCD(N2,N1)=UDOCD(N2,N1)-COD
      UDICD(N2,N1)=UDICD(N2,N1)-CXD
      UDOND(N2,N1)=UDOND(N2,N1)-ZOD
      UDIND(N2,N1)=UDIND(N2,N1)-ZXD
      UDOPD(N2,N1)=UDOPD(N2,N1)-POD
      UDIPD(N2,N1)=UDIPD(N2,N1)-PXD
C     IF(J.EQ.15.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,3489)'UDICD',I,J,NFZ,N,N2,N1,N6,N5,N4,UDICD(N2,N1),CXD
C    2,XN,FLW(N,N6,N5,N4),XCOFLS(N,N6,N5,N4),XCOFHS(N,N6,N5,N4)
C    2,XCOFLG(N,N6,N5,N4),XCHFLS(N,N6,N5,N4)
C    3,XCHFHS(N,N6,N5,N4),XCHFLG(N,N6,N5,N4)
C     WRITE(*,3489)'UDIND',I,J,NFZ,N,N2,N1,N6,N5,N4,UDIND(N2,N1),ZXD
C    2,XN,XN4FLW(N,N6,N5,N4),XN3FLW(N,N6,N5,N4),XNOFLW(N,N6,N5,N4)
C    2,XN4FLB(N,N6,N5,N4),XN3FLB(N,N6,N5,N4),XNOFLB(N,N6,N5,N4)
C    3,XNXFLS(N,N6,N5,N4),XNXFLB(N,N6,N5,N4)
C    5,XN4FHW(N,N6,N5,N4),XN3FHW(N,N6,N5,N4),XNOFHW(N,N6,N5,N4)
C    6,XN4FHB(N,N6,N5,N4),XN3FHB(N,N6,N5,N4),XNOFHB(N,N6,N5,N4)
C    7,XNXFHS(N,N6,N5,N4),XNXFHB(N,N6,N5,N4)
3489  FORMAT(A8,9I4,30E12.4)
C     ENDIF
C
C     SUBSURFACE BOUNDARY FLUXES OF O2
C
C     X*FLS,X*FHS=solute flux in macropores,micropores from trnsfr.f
C     X*FLG=convective+diffusive gas flux from trnsfr.f
C     solute code:OX=O2,HG=H2
C     OXYGOU,H2GOU=cumulative O2,H2 loss through lateral and lower
C        boundaries
C
      OOD=XN*(XOXFLS(N,N6,N5,N4)+XOXFHS(N,N6,N5,N4)
     2+XOXFLG(N,N6,N5,N4))
      OXYGOU=OXYGOU-OOD
      HOD=XN*(XHGFLS(N,N6,N5,N4)+XHGFHS(N,N6,N5,N4)
     2+XHGFLG(N,N6,N5,N4))
      H2GOU=H2GOU-HOD
C     WRITE(*,6643)'OXYGOU',I,J,NFZ,N,N6,N5,N4
C    2,OXYGOU,OOD,XN,XOXFLS(N,N6,N5,N4),XOXFHS(N,N6,N5,N4)
C    3,XOXFLG(N,N6,N5,N4)
6643  FORMAT(A8,7I4,20E12.4)
C
C     SUBSURFACE BOUNDARY FLUXES OF SOLUTES
C
C     X*FLS=convective+diffusive solute flux through micropores
C        from trnsfrs.f
C     X*FLW,X*FLB=convective+diffusive solute flux through micropores
C        in non-band,band from trnsfrs.f
C     X*FHS=convective+diffusive solute flux through macropores
C        from trnsfrs.f
C     X*FHW,X*FHB=convective+diffusive solute flux through macropores
C        in non-band,band from trnsfrs.f
C     TIONOU,UIONOU=cumulative salt loss through lateral and lower
C        boundaries
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C        :*1=non-band,*B=band
C
      IF(ISALTG.NE.0)THEN
      PQD=XN*31.0*(XH0PFS(N,N6,N5,N4)+XH0BFB(N,N6,N5,N4)
     2+XC0PFS(N,N6,N5,N4)+XC0BFB(N,N6,N5,N4)+XF1PFS(N,N6,N5,N4)
     3+XC1PFS(N,N6,N5,N4)+XM1PFS(N,N6,N5,N4)+XF1BFB(N,N6,N5,N4)
     4+XC1BFB(N,N6,N5,N4)+XM1BFB(N,N6,N5,N4)+XH3PFS(N,N6,N5,N4)
     5+XF2PFS(N,N6,N5,N4)+XC2PFS(N,N6,N5,N4)+XH3BFB(N,N6,N5,N4)
     6+XF2BFB(N,N6,N5,N4)+XC2BFB(N,N6,N5,N4))
      PHD=XN*31.0*(XH0PHS(N,N6,N5,N4)+XH0BHB(N,N6,N5,N4)
     2+XC0PHS(N,N6,N5,N4)+XC0BHB(N,N6,N5,N4)+XF1PHS(N,N6,N5,N4)
     3+XC1PHS(N,N6,N5,N4)+XM1PHS(N,N6,N5,N4)+XF1BHB(N,N6,N5,N4)
     4+XC1BHB(N,N6,N5,N4)+XM1BHB(N,N6,N5,N4)+XH3PHS(N,N6,N5,N4)
     7+XF2PHS(N,N6,N5,N4)+XC2PHS(N,N6,N5,N4)+XH3BHB(N,N6,N5,N4)
     8+XF2BHB(N,N6,N5,N4)+XC2BHB(N,N6,N5,N4))
      TPOU=TPOU-PQD-PHD
      SSD=XN*(XALFLS(N,N6,N5,N4)+XFEFLS(N,N6,N5,N4)+XHYFLS(N,N6,N5,N4)
     2+XCAFLS(N,N6,N5,N4)+XMGFLS(N,N6,N5,N4)+XNAFLS(N,N6,N5,N4)
     3+XKAFLS(N,N6,N5,N4)+XOHFLS(N,N6,N5,N4)+XSOFLS(N,N6,N5,N4)
     4+XCLFLS(N,N6,N5,N4)+XC3FLS(N,N6,N5,N4)+XH0PFS(N,N6,N5,N4)
     5+XH0BFB(N,N6,N5,N4)
     6+2.0*(XHCFLS(N,N6,N5,N4)+XAL1FS(N,N6,N5,N4)
     6+XALSFS(N,N6,N5,N4)+XFE1FS(N,N6,N5,N4)+XFESFS(N,N6,N5,N4)
     7+XCAOFS(N,N6,N5,N4)+XCACFS(N,N6,N5,N4)
     8+XCASFS(N,N6,N5,N4)+XMGOFS(N,N6,N5,N4)+XMGCFS(N,N6,N5,N4)
     9+XMGSFS(N,N6,N5,N4)+XNACFS(N,N6,N5,N4)+XNASFS(N,N6,N5,N4)
     1+XKASFS(N,N6,N5,N4)+XC0PFS(N,N6,N5,N4)+XC0BFB(N,N6,N5,N4))
     3+3.0*(XAL2FS(N,N6,N5,N4)
     3+XFE2FS(N,N6,N5,N4)+XCAHFS(N,N6,N5,N4)+XMGHFS(N,N6,N5,N4)
     4+XF1PFS(N,N6,N5,N4)+XC1PFS(N,N6,N5,N4)+XM1PFS(N,N6,N5,N4)
     5+XF1BFB(N,N6,N5,N4)+XC1BFB(N,N6,N5,N4)+XM1BFB(N,N6,N5,N4))
     6+4.0*(XAL3FS(N,N6,N5,N4)+XFE3FS(N,N6,N5,N4)+XH3PFS(N,N6,N5,N4)
     7+XF2PFS(N,N6,N5,N4)+XC2PFS(N,N6,N5,N4)+XH3BFB(N,N6,N5,N4)
     8+XF2BFB(N,N6,N5,N4)+XC2BFB(N,N6,N5,N4))
     9+5.0*(XAL4FS(N,N6,N5,N4)+XFE4FS(N,N6,N5,N4)))
      SHD=XN*(XALFHS(N,N6,N5,N4)+XFEFHS(N,N6,N5,N4)+XHYFHS(N,N6,N5,N4)
     2+XCAFHS(N,N6,N5,N4)+XMGFHS(N,N6,N5,N4)+XNAFHS(N,N6,N5,N4)
     3+XKAFHS(N,N6,N5,N4)+XOHFHS(N,N6,N5,N4)+XSOFHS(N,N6,N5,N4)
     4+XCLFHS(N,N6,N5,N4)+XC3FHS(N,N6,N5,N4)+XH0PHS(N,N6,N5,N4)
     5+XH0BHB(N,N6,N5,N4)
     5+2.0*(XHCFHS(N,N6,N5,N4)+XAL1HS(N,N6,N5,N4)
     6+XALSHS(N,N6,N5,N4)+XFE1HS(N,N6,N5,N4)+XFESHS(N,N6,N5,N4)
     7+XCAOHS(N,N6,N5,N4)+XCACHS(N,N6,N5,N4)
     8+XCASHS(N,N6,N5,N4)+XMGOHS(N,N6,N5,N4)+XMGCHS(N,N6,N5,N4)
     9+XMGSHS(N,N6,N5,N4)+XNACHS(N,N6,N5,N4)+XNASHS(N,N6,N5,N4)
     1+XKASHS(N,N6,N5,N4)+XC0PHS(N,N6,N5,N4)+XC0BHB(N,N6,N5,N4))
     3+3.0*(XAL2HS(N,N6,N5,N4)
     3+XFE2HS(N,N6,N5,N4)+XCAHHS(N,N6,N5,N4)+XMGHHS(N,N6,N5,N4)
     4+XF1PHS(N,N6,N5,N4)+XC1PHS(N,N6,N5,N4)+XM1PHS(N,N6,N5,N4)
     5+XF1BHB(N,N6,N5,N4)+XC1BHB(N,N6,N5,N4)+XM1BHB(N,N6,N5,N4))
     6+4.0*(XAL3HS(N,N6,N5,N4)+XFE3HS(N,N6,N5,N4)+XH3PHS(N,N6,N5,N4)
     7+XF2PHS(N,N6,N5,N4)+XC2PHS(N,N6,N5,N4)+XH3BHB(N,N6,N5,N4)
     8+XF2BHB(N,N6,N5,N4)+XC2BHB(N,N6,N5,N4))
     9+5.0*(XAL4HS(N,N6,N5,N4)+XAL4HS(N,N6,N5,N4)))
      SO=SSD+SHD
      TIONOU=TIONOU-SO
      UIONOU(N2,N1)=UIONOU(N2,N1)-SO
C     IF(I.EQ.180.AND.J.EQ.12)THEN
C     WRITE(*,3337)'SSD',I,J,N,N6,N5,N4,SSD
C    2,XALFLS(N,N6,N5,N4),XFEFLS(N,N6,N5,N4),XHYFLS(N,N6,N5,N4)
C    2,XCAFLS(N,N6,N5,N4),XMGFLS(N,N6,N5,N4),XNAFLS(N,N6,N5,N4)
C    3,XKAFLS(N,N6,N5,N4),XOHFLS(N,N6,N5,N4),XSOFLS(N,N6,N5,N4)
C    4,XCLFLS(N,N6,N5,N4),XC3FLS(N,N6,N5,N4),XH0PFS(N,N6,N5,N4)
C    5,XH0BFB(N,N6,N5,N4)
C    6,XHCFLS(N,N6,N5,N4),XAL1FS(N,N6,N5,N4)
C    6,XALSFS(N,N6,N5,N4),XFE1FS(N,N6,N5,N4),XFESFS(N,N6,N5,N4)
C    7,XCAOFS(N,N6,N5,N4),XCACFS(N,N6,N5,N4)
C    8,XCASFS(N,N6,N5,N4),XMGOFS(N,N6,N5,N4),XMGCFS(N,N6,N5,N4)
C    9,XMGSFS(N,N6,N5,N4),XNACFS(N,N6,N5,N4),XNASFS(N,N6,N5,N4)
C    1,XKASFS(N,N6,N5,N4),XC0PFS(N,N6,N5,N4),XC0BFB(N,N6,N5,N4)
C    3,XAL2FS(N,N6,N5,N4)
C    3,XFE2FS(N,N6,N5,N4),XCAHFS(N,N6,N5,N4),XMGHFS(N,N6,N5,N4)
C    4,XF1PFS(N,N6,N5,N4),XC1PFS(N,N6,N5,N4),XM1PFS(N,N6,N5,N4)
C    5,XF1BFB(N,N6,N5,N4),XC1BFB(N,N6,N5,N4),XM1BFB(N,N6,N5,N4)
C    6,XAL3FS(N,N6,N5,N4),XFE3FS(N,N6,N5,N4),XH3PFS(N,N6,N5,N4)
C    7,XF2PFS(N,N6,N5,N4),XC2PFS(N,N6,N5,N4),XH3BFB(N,N6,N5,N4)
C    8,XF2BFB(N,N6,N5,N4),XC2BFB(N,N6,N5,N4)
C    9,XAL4FS(N,N6,N5,N4),XFE4FS(N,N6,N5,N4)
3337  FORMAT(A8,6I4,80E12.4)
C     ENDIF
C
C     SUBSURFACE FLUX ELECTRICAL CONDUCTIVITY
C
C     FLW,FLWH=micropore,macropore flux through lateral and lower
C        boundaries from watsub.f
C     X*FLS=convective+diffusive solute flux through micropores
C        from trnsfrs.f
C     X*FLW,X*FLB=convective+diffusive solute flux through micropores
C        in non-band,band from trnsfrs.f
C     X*FHS=convective+diffusive solute flux through macropores
C        from trnsfrs.f
C     X*FHW,X*FHB=convective+diffusive solute flux through macropores
C        in non-band,band from trnsfrs.f
C     EC*=electrical conductivity of ion
C     ECNDQ=total electrical conductivity of water flux
C     ion code:HY=H,OH=OH,AL=Al,FE=Fe,CA=Ca,MG=Mg,NA=Na,KA=K
C             :CO=CO3,HC=HCO3,SO=SO4,CL=Cl
C
      WX=FLW(N,N6,N5,N4)+FLWH(N,N6,N5,N4)
      IF(ABS(WX).GT.ZEROS(N2,N1))THEN
      ECHY=0.337*AMAX1(0.0,(XHYFLS(N,N6,N5,N4)
     2+XHYFHS(N,N6,N5,N4))/WX)
      ECOH=0.192*AMAX1(0.0,(XOHFLS(N,N6,N5,N4)
     2+XOHFHS(N,N6,N5,N4))/WX)
      ECAL=0.056*AMAX1(0.0,(XALFLS(N,N6,N5,N4)
     2+XCAFHS(N,N6,N5,N4))*3.0/WX)
      ECFE=0.051*AMAX1(0.0,(XFEFLS(N,N6,N5,N4)
     2+XFEFHS(N,N6,N5,N4))*3.0/WX)
      ECCA=0.060*AMAX1(0.0,(XCAFLS(N,N6,N5,N4)
     2+XCAFHS(N,N6,N5,N4))*2.0/WX)
      ECMG=0.053*AMAX1(0.0,(XMGFLS(N,N6,N5,N4)
     2+XMGFHS(N,N6,N5,N4))*2.0/WX)
      ECNA=0.050*AMAX1(0.0,(XNAFLS(N,N6,N5,N4)
     2+XNAFHS(N,N6,N5,N4))/WX)
      ECKA=0.070*AMAX1(0.0,(XKAFLS(N,N6,N5,N4)
     2+XKAFHS(N,N6,N5,N4))/WX)
      ECCO=0.072*AMAX1(0.0,(XC3FLS(N,N6,N5,N4)
     2+XC3FHS(N,N6,N5,N4))*2.0/WX)
      ECHC=0.044*AMAX1(0.0,(XHCFLS(N,N6,N5,N4)
     2+XHCFHS(N,N6,N5,N4))/WX)
      ECSO=0.080*AMAX1(0.0,(XSOFLS(N,N6,N5,N4)
     2+XSOFHS(N,N6,N5,N4))*2.0/WX)
      ECCL=0.076*AMAX1(0.0,(XCLFLS(N,N6,N5,N4)
     2+XCLFHS(N,N6,N5,N4))/WX)
      ECNO=0.071*AMAX1(0.0,(XNOFLW(N,N6,N5,N4)
     2+XNOFHW(N,N6,N5,N4))/(WX*14.0))
      ECNDX=ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA
     2+ECCO+ECHC+ECSO+ECCL+ECNO
C     IF((I/10)*10.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,9992)'ECNDX',IYRC,I,J,N4,N5,N6,N,WX,ECNDX
C    2,FLW(N,N6,N5,N4),FLWH(N,N6,N5,N4)
9992  FORMAT(A8,7I4,4E12.4)
C     ENDIF
      ELSE
      ECNDX=0.0
      ENDIF
      ENDIF
      SG=SG+XHGFLS(N,N6,N5,N4)+XHGFLG(N,N6,N5,N4)
      ENDIF
      ENDIF
9975  CONTINUE
9980  CONTINUE
9985  CONTINUE
C
C     INITIALIZE NET WATER AND HEAT FLUXES FOR RUNOFF
C
      TQR(NY,NX)=0.0
      THQR(NY,NX)=0.0
      TQS(NY,NX)=0.0
      TQW(NY,NX)=0.0
      TQI(NY,NX)=0.0
      THQS(NY,NX)=0.0
C
C     INITIALIZE NET SOLUTE AND GAS FLUXES FOR RUNOFF
C
      DO 9960 K=0,2
      TOCQRS(K,NY,NX)=0.0
      TONQRS(K,NY,NX)=0.0
      TOPQRS(K,NY,NX)=0.0
      TOAQRS(K,NY,NX)=0.0
9960  CONTINUE
      TCOQRS(NY,NX)=0.0
      TCHQRS(NY,NX)=0.0
      TOXQRS(NY,NX)=0.0
      TNGQRS(NY,NX)=0.0
      TN2QRS(NY,NX)=0.0
      THGQRS(NY,NX)=0.0
      TN4QRS(NY,NX)=0.0
      TN3QRS(NY,NX)=0.0
      TNOQRS(NY,NX)=0.0
      TNXQRS(NY,NX)=0.0
      TP1QRS(NY,NX)=0.0
      TPOQRS(NY,NX)=0.0
      TCOQSS(NY,NX)=0.0
      TCHQSS(NY,NX)=0.0
      TOXQSS(NY,NX)=0.0
      TNGQSS(NY,NX)=0.0
      TN2QSS(NY,NX)=0.0
      TN4QSS(NY,NX)=0.0
      TN3QSS(NY,NX)=0.0
      TNOQSS(NY,NX)=0.0
      TP1QSS(NY,NX)=0.0
      TPOQSS(NY,NX)=0.0
      DO 9955 L=1,JS
      TCOBLS(L,NY,NX)=0.0
      TCHBLS(L,NY,NX)=0.0
      TOXBLS(L,NY,NX)=0.0
      TNGBLS(L,NY,NX)=0.0
      TN2BLS(L,NY,NX)=0.0
      TN4BLW(L,NY,NX)=0.0
      TN3BLW(L,NY,NX)=0.0
      TNOBLW(L,NY,NX)=0.0
      TH1PBS(L,NY,NX)=0.0
      TH2PBS(L,NY,NX)=0.0
9955  CONTINUE
      IF(ISALTG.NE.0)THEN
      TQRAL(NY,NX)=0.0
      TQRFE(NY,NX)=0.0
      TQRHY(NY,NX)=0.0
      TQRCA(NY,NX)=0.0
      TQRMG(NY,NX)=0.0
      TQRNA(NY,NX)=0.0
      TQRKA(NY,NX)=0.0
      TQROH(NY,NX)=0.0
      TQRSO(NY,NX)=0.0
      TQRCL(NY,NX)=0.0
      TQRC3(NY,NX)=0.0
      TQRHC(NY,NX)=0.0
      TQRAL1(NY,NX)=0.0
      TQRAL2(NY,NX)=0.0
      TQRAL3(NY,NX)=0.0
      TQRAL4(NY,NX)=0.0
      TQRALS(NY,NX)=0.0
      TQRFE1(NY,NX)=0.0
      TQRFE2(NY,NX)=0.0
      TQRFE3(NY,NX)=0.0
      TQRFE4(NY,NX)=0.0
      TQRFES(NY,NX)=0.0
      TQRCAO(NY,NX)=0.0
      TQRCAC(NY,NX)=0.0
      TQRCAH(NY,NX)=0.0
      TQRCAS(NY,NX)=0.0
      TQRMGO(NY,NX)=0.0
      TQRMGC(NY,NX)=0.0
      TQRMGH(NY,NX)=0.0
      TQRMGS(NY,NX)=0.0
      TQRNAC(NY,NX)=0.0
      TQRNAS(NY,NX)=0.0
      TQRKAS(NY,NX)=0.0
      TQRH0P(NY,NX)=0.0
      TQRH3P(NY,NX)=0.0
      TQRF1P(NY,NX)=0.0
      TQRF2P(NY,NX)=0.0
      TQRC0P(NY,NX)=0.0
      TQRC1P(NY,NX)=0.0
      TQRC2P(NY,NX)=0.0
      TQRM1P(NY,NX)=0.0
      TQSAL(NY,NX)=0.0
      TQSFE(NY,NX)=0.0
      TQSHY(NY,NX)=0.0
      TQSCA(NY,NX)=0.0
      TQSMG(NY,NX)=0.0
      TQSNA(NY,NX)=0.0
      TQSKA(NY,NX)=0.0
      TQSOH(NY,NX)=0.0
      TQSSO(NY,NX)=0.0
      TQSCL(NY,NX)=0.0
      TQSC3(NY,NX)=0.0
      TQSHC(NY,NX)=0.0
      TQSAL1(NY,NX)=0.0
      TQSAL2(NY,NX)=0.0
      TQSAL3(NY,NX)=0.0
      TQSAL4(NY,NX)=0.0
      TQSALS(NY,NX)=0.0
      TQSFE1(NY,NX)=0.0
      TQSFE2(NY,NX)=0.0
      TQSFE3(NY,NX)=0.0
      TQSFE4(NY,NX)=0.0
      TQSFES(NY,NX)=0.0
      TQSCAO(NY,NX)=0.0
      TQSCAC(NY,NX)=0.0
      TQSCAH(NY,NX)=0.0
      TQSCAS(NY,NX)=0.0
      TQSMGO(NY,NX)=0.0
      TQSMGC(NY,NX)=0.0
      TQSMGH(NY,NX)=0.0
      TQSMGS(NY,NX)=0.0
      TQSNAC(NY,NX)=0.0
      TQSNAS(NY,NX)=0.0
      TQSKAS(NY,NX)=0.0
      TQSH0P(NY,NX)=0.0
      TQSH3P(NY,NX)=0.0
      TQSF1P(NY,NX)=0.0
      TQSF2P(NY,NX)=0.0
      TQSC0P(NY,NX)=0.0
      TQSC1P(NY,NX)=0.0
      TQSC2P(NY,NX)=0.0
      TQSM1P(NY,NX)=0.0
C
C     INITIALIZE NET SOLUTE AND GAS FLUXES FROM SNOWPACK DRIFT
C
      DO 9950 L=1,JS
      TALBLS(L,NY,NX)=0.0
      TFEBLS(L,NY,NX)=0.0
      THYBLS(L,NY,NX)=0.0
      TCABLS(L,NY,NX)=0.0
      TMGBLS(L,NY,NX)=0.0
      TNABLS(L,NY,NX)=0.0
      TKABLS(L,NY,NX)=0.0
      TOHBLS(L,NY,NX)=0.0
      TSOBLS(L,NY,NX)=0.0
      TCLBLS(L,NY,NX)=0.0
      TC3BLS(L,NY,NX)=0.0
      THCBLS(L,NY,NX)=0.0
      TAL1BS(L,NY,NX)=0.0
      TAL2BS(L,NY,NX)=0.0
      TAL3BS(L,NY,NX)=0.0
      TAL4BS(L,NY,NX)=0.0
      TALSBS(L,NY,NX)=0.0
      TFE1BS(L,NY,NX)=0.0
      TFE2BS(L,NY,NX)=0.0
      TFE3BS(L,NY,NX)=0.0
      TFE4BS(L,NY,NX)=0.0
      TFESBS(L,NY,NX)=0.0
      TCAOBS(L,NY,NX)=0.0
      TCACBS(L,NY,NX)=0.0
      TCAHBS(L,NY,NX)=0.0
      TCASBS(L,NY,NX)=0.0
      TMGOBS(L,NY,NX)=0.0
      TMGCBS(L,NY,NX)=0.0
      TMGHBS(L,NY,NX)=0.0
      TMGSBS(L,NY,NX)=0.0
      TNACBS(L,NY,NX)=0.0
      TNASBS(L,NY,NX)=0.0
      TKASBS(L,NY,NX)=0.0
      TH0PBS(L,NY,NX)=0.0
      TH3PBS(L,NY,NX)=0.0
      TF1PBS(L,NY,NX)=0.0
      TF2PBS(L,NY,NX)=0.0
      TC0PBS(L,NY,NX)=0.0
      TC1PBS(L,NY,NX)=0.0
      TC2PBS(L,NY,NX)=0.0
      TM1PBS(L,NY,NX)=0.0
9950  CONTINUE
      ENDIF
C
C     INITIALIZE NET SEDIMENT FLUXES FROM EROSION
C
      IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
      TSEDER(NY,NX)=0.0
      TSANER(NY,NX)=0.0
      TSILER(NY,NX)=0.0
      TCLAER(NY,NX)=0.0
      TCECER(NY,NX)=0.0
      TAECER(NY,NX)=0.0
      TNH4ER(NY,NX)=0.0
      TNH3ER(NY,NX)=0.0
      TNHUER(NY,NX)=0.0
      TNO3ER(NY,NX)=0.0
      TNH4EB(NY,NX)=0.0
      TNH3EB(NY,NX)=0.0
      TNHUEB(NY,NX)=0.0
      TNO3EB(NY,NX)=0.0
      TN4ER(NY,NX)=0.0
      TNBER(NY,NX)=0.0
      THYER(NY,NX)=0.0
      TALER(NY,NX)=0.0
      TFEER(NY,NX)=0.0
      TCAER(NY,NX)=0.0
      TMGER(NY,NX)=0.0
      TNAER(NY,NX)=0.0
      TKAER(NY,NX)=0.0
      THCER(NY,NX)=0.0
      TAL2ER(NY,NX)=0.0
      TFE2ER(NY,NX)=0.0
      TOH0ER(NY,NX)=0.0
      TOH1ER(NY,NX)=0.0
      TOH2ER(NY,NX)=0.0
      TH1PER(NY,NX)=0.0
      TH2PER(NY,NX)=0.0
      TOH0EB(NY,NX)=0.0
      TOH1EB(NY,NX)=0.0
      TOH2EB(NY,NX)=0.0
      TH1PEB(NY,NX)=0.0
      TH2PEB(NY,NX)=0.0
      TALOER(NY,NX)=0.0
      TFEOER(NY,NX)=0.0
      TCACER(NY,NX)=0.0
      TCASER(NY,NX)=0.0
      TALPER(NY,NX)=0.0
      TFEPER(NY,NX)=0.0
      TCPDER(NY,NX)=0.0
      TCPHER(NY,NX)=0.0
      TCPMER(NY,NX)=0.0
      TALPEB(NY,NX)=0.0
      TFEPEB(NY,NX)=0.0
      TCPDEB(NY,NX)=0.0
      TCPHEB(NY,NX)=0.0
      TCPMEB(NY,NX)=0.0
      DO 9480 K=0,5
      DO 9480 NO=1,7
      DO 9480 M=1,3
      TOMCER(M,NO,K,NY,NX)=0.0
      TOMNER(M,NO,K,NY,NX)=0.0
      TOMPER(M,NO,K,NY,NX)=0.0
9480  CONTINUE
      DO 9475 K=0,4
      DO 9470 M=1,2
      TORCER(M,K,NY,NX)=0.0
      TORNER(M,K,NY,NX)=0.0
      TORPER(M,K,NY,NX)=0.0
9470  CONTINUE
      TOHCER(K,NY,NX)=0.0
      TOHNER(K,NY,NX)=0.0
      TOHPER(K,NY,NX)=0.0
      TOHAER(K,NY,NX)=0.0
      DO 9465 M=1,5
      TOSCER(M,K,NY,NX)=0.0
      TOSAER(M,K,NY,NX)=0.0
      TOSNER(M,K,NY,NX)=0.0
      TOSPER(M,K,NY,NX)=0.0
9465  CONTINUE
9475  CONTINUE
      ENDIF
C
C     INITIALIZE NET SNOWPACK FLUXES WITHIN SNOWPACK
C
      DO 8475 L=1,JS
      TFLWS(L,NY,NX)=0.0
      TFLWW(L,NY,NX)=0.0
      TFLWV(L,NY,NX)=0.0
      TFLWI(L,NY,NX)=0.0
      THFLWW(L,NY,NX)=0.0
8475  CONTINUE
      LG=0
      LX=0
      DO 8575 L=NU(NY,NX),NL(NY,NX)
C
C     IDENTIFY LAYERS FOR BUBBLE FLUX TRANSFER
C
C     LG=lowest soil layer with gas phase
C     V*G2=molar gas content
C     *G=soil gas content
C     VOLP=soil air-filled porosity
C     VTATM=molar gas content at atmospheric pressure
C     VTGAS=total molar gas contest
C     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
C             :*ZN3*=NH3,*H2G*=H2
C
      VCO2G2=CO2G(L,NY,NX)/12.0
      VCH4G2=CH4G(L,NY,NX)/12.0
      VOXYG2=OXYG(L,NY,NX)/32.0
      VZ2GG2=Z2GG(L,NY,NX)/28.0
      VZ2OG2=Z2OG(L,NY,NX)/28.0
      VNH3G2=ZNH3G(L,NY,NX)/14.0
      VH2GG2=H2GG(L,NY,NX)/2.0
      VTATM=AMAX1(0.0,1.2194E+04*VOLP(L,NY,NX)/TKS(L,NY,NX))
      VTGAS=VCO2G2+VCH4G2+VOXYG2+VZ2GG2+VZ2OG2+VNH3G2+VH2GG2
      IF(THETP(L,NY,NX).LT.THETX.OR.VTGAS.GT.VTATM)LX=1
      IF(THETP(L,NY,NX).GE.THETX.AND.LX.EQ.0)LG=L
C     WRITE(*,5431)'LG',I,J,NX,NY,L,LG,LX,THETP(L,NY,NX),THETX
C    2,VOLP(L,NY,NX),TKS(L,NY,NX),VTGAS,VTATM
5431  FORMAT(A8,7I4,12E12.4)
      VOLW1(L,NY,NX)=VOLW(L,NY,NX)
      VOLI1(L,NY,NX)=VOLI(L,NY,NX)
      VOLWH1(L,NY,NX)=VOLWH(L,NY,NX)
      VOLIH1(L,NY,NX)=VOLIH(L,NY,NX)
C
C     INITIALIZE WATER AND HEAT NET FLUX ACCUMULATORS WITHIN SOIL
C
      TFLW(L,NY,NX)=0.0
      TFLV(L,NY,NX)=0.0
      TFLWX(L,NY,NX)=0.0
      TFLWH(L,NY,NX)=0.0
      THFLW(L,NY,NX)=0.0
      TWFLFL(L,NY,NX)=0.0
      TWFLFH(L,NY,NX)=0.0
      THFLFL(L,NY,NX)=0.0
      TWFLVL(L,NY,NX)=0.0
      THFLVL(L,NY,NX)=0.0
C
C     INITIALIZE GAS AND SOLUTE NET FLUX ACCUMULATORS WITHIN SOIL
C
      DO 8595 K=0,4
      TOCFLS(K,L,NY,NX)=0.0
      TONFLS(K,L,NY,NX)=0.0
      TOPFLS(K,L,NY,NX)=0.0
      TOAFLS(K,L,NY,NX)=0.0
      TOCFHS(K,L,NY,NX)=0.0
      TONFHS(K,L,NY,NX)=0.0
      TOPFHS(K,L,NY,NX)=0.0
      TOAFHS(K,L,NY,NX)=0.0
8595  CONTINUE
      TCOFLS(L,NY,NX)=0.0
      TCHFLS(L,NY,NX)=0.0
      TOXFLS(L,NY,NX)=0.0
      TNGFLS(L,NY,NX)=0.0
      TN2FLS(L,NY,NX)=0.0
      THGFLS(L,NY,NX)=0.0
      TN4FLS(L,NY,NX)=0.0
      TN3FLS(L,NY,NX)=0.0
      TNOFLS(L,NY,NX)=0.0
      TNXFLS(L,NY,NX)=0.0
      TP1FLS(L,NY,NX)=0.0
      TPOFLS(L,NY,NX)=0.0
      TN4FLB(L,NY,NX)=0.0
      TN3FLB(L,NY,NX)=0.0
      TNOFLB(L,NY,NX)=0.0
      TNXFLB(L,NY,NX)=0.0
      TH1BFB(L,NY,NX)=0.0
      TH2BFB(L,NY,NX)=0.0
      TCOFHS(L,NY,NX)=0.0
      TCHFHS(L,NY,NX)=0.0
      TOXFHS(L,NY,NX)=0.0
      TNGFHS(L,NY,NX)=0.0
      TN2FHS(L,NY,NX)=0.0
      THGFHS(L,NY,NX)=0.0
      TN4FHS(L,NY,NX)=0.0
      TN3FHS(L,NY,NX)=0.0
      TNOFHS(L,NY,NX)=0.0
      TNXFHS(L,NY,NX)=0.0
      TP1FHS(L,NY,NX)=0.0
      TPOFHS(L,NY,NX)=0.0
      TN4FHB(L,NY,NX)=0.0
      TN3FHB(L,NY,NX)=0.0
      TNOFHB(L,NY,NX)=0.0
      TNXFHB(L,NY,NX)=0.0
      TH1BHB(L,NY,NX)=0.0
      TH2BHB(L,NY,NX)=0.0
      TVPFLG(L,NY,NX)=0.0
      TCOFLG(L,NY,NX)=0.0
      TCHFLG(L,NY,NX)=0.0
      TOXFLG(L,NY,NX)=0.0
      TNGFLG(L,NY,NX)=0.0
      TN2FLG(L,NY,NX)=0.0
      TNHFLG(L,NY,NX)=0.0
      THGFLG(L,NY,NX)=0.0
      IF(ISALTG.NE.0)THEN
      TALFLS(L,NY,NX)=0.0
      TFEFLS(L,NY,NX)=0.0
      THYFLS(L,NY,NX)=0.0
      TCAFLS(L,NY,NX)=0.0
      TMGFLS(L,NY,NX)=0.0
      TNAFLS(L,NY,NX)=0.0
      TKAFLS(L,NY,NX)=0.0
      TOHFLS(L,NY,NX)=0.0
      TSOFLS(L,NY,NX)=0.0
      TCLFLS(L,NY,NX)=0.0
      TC3FLS(L,NY,NX)=0.0
      THCFLS(L,NY,NX)=0.0
      TAL1FS(L,NY,NX)=0.0
      TAL2FS(L,NY,NX)=0.0
      TAL3FS(L,NY,NX)=0.0
      TAL4FS(L,NY,NX)=0.0
      TALSFS(L,NY,NX)=0.0
      TFE1FS(L,NY,NX)=0.0
      TFE2FS(L,NY,NX)=0.0
      TFE3FS(L,NY,NX)=0.0
      TFE4FS(L,NY,NX)=0.0
      TFESFS(L,NY,NX)=0.0
      TCAOFS(L,NY,NX)=0.0
      TCACFS(L,NY,NX)=0.0
      TCAHFS(L,NY,NX)=0.0
      TCASFS(L,NY,NX)=0.0
      TMGOFS(L,NY,NX)=0.0
      TMGCFS(L,NY,NX)=0.0
      TMGHFS(L,NY,NX)=0.0
      TMGSFS(L,NY,NX)=0.0
      TNACFS(L,NY,NX)=0.0
      TNASFS(L,NY,NX)=0.0
      TKASFS(L,NY,NX)=0.0
      TH0PFS(L,NY,NX)=0.0
      TH3PFS(L,NY,NX)=0.0
      TF1PFS(L,NY,NX)=0.0
      TF2PFS(L,NY,NX)=0.0
      TC0PFS(L,NY,NX)=0.0
      TC1PFS(L,NY,NX)=0.0
      TC2PFS(L,NY,NX)=0.0
      TM1PFS(L,NY,NX)=0.0
      TH0BFB(L,NY,NX)=0.0
      TH3BFB(L,NY,NX)=0.0
      TF1BFB(L,NY,NX)=0.0
      TF2BFB(L,NY,NX)=0.0
      TC0BFB(L,NY,NX)=0.0
      TC1BFB(L,NY,NX)=0.0
      TC2BFB(L,NY,NX)=0.0
      TM1BFB(L,NY,NX)=0.0
      TALFHS(L,NY,NX)=0.0
      TFEFHS(L,NY,NX)=0.0
      THYFHS(L,NY,NX)=0.0
      TCAFHS(L,NY,NX)=0.0
      TMGFHS(L,NY,NX)=0.0
      TNAFHS(L,NY,NX)=0.0
      TKAFHS(L,NY,NX)=0.0
      TOHFHS(L,NY,NX)=0.0
      TSOFHS(L,NY,NX)=0.0
      TCLFHS(L,NY,NX)=0.0
      TC3FHS(L,NY,NX)=0.0
      THCFHS(L,NY,NX)=0.0
      TAL1HS(L,NY,NX)=0.0
      TAL2HS(L,NY,NX)=0.0
      TAL3HS(L,NY,NX)=0.0
      TAL4HS(L,NY,NX)=0.0
      TALSHS(L,NY,NX)=0.0
      TFE1HS(L,NY,NX)=0.0
      TFE2HS(L,NY,NX)=0.0
      TFE3HS(L,NY,NX)=0.0
      TFE4HS(L,NY,NX)=0.0
      TFESHS(L,NY,NX)=0.0
      TCAOHS(L,NY,NX)=0.0
      TCACHS(L,NY,NX)=0.0
      TCAHHS(L,NY,NX)=0.0
      TCASHS(L,NY,NX)=0.0
      TMGOHS(L,NY,NX)=0.0
      TMGCHS(L,NY,NX)=0.0
      TMGHHS(L,NY,NX)=0.0
      TMGSHS(L,NY,NX)=0.0
      TNACHS(L,NY,NX)=0.0
      TNASHS(L,NY,NX)=0.0
      TKASHS(L,NY,NX)=0.0
      TH0PHS(L,NY,NX)=0.0
      TH3PHS(L,NY,NX)=0.0
      TF1PHS(L,NY,NX)=0.0
      TF2PHS(L,NY,NX)=0.0
      TC0PHS(L,NY,NX)=0.0
      TC1PHS(L,NY,NX)=0.0
      TC2PHS(L,NY,NX)=0.0
      TM1PHS(L,NY,NX)=0.0
      TH0BHB(L,NY,NX)=0.0
      TH3BHB(L,NY,NX)=0.0
      TF1BHB(L,NY,NX)=0.0
      TF2BHB(L,NY,NX)=0.0
      TC0BHB(L,NY,NX)=0.0
      TC1BHB(L,NY,NX)=0.0
      TC2BHB(L,NY,NX)=0.0
      TM1BHB(L,NY,NX)=0.0
      ENDIF
C
C     NET WATER, HEAT, GAS, SOLUTE, SEDIMENT FLUX
C
C     N3,N2,N1=L,NY,NX of source grid cell
C     N6,N5,N4=L,NY,NX of destination grid cell
C
      N1=NX
      N2=NY
      N3=L
      DO 8580 N=1,3
      IF(N.EQ.1)THEN
      N4=NX+1
      N5=NY
      N4B=NX-1
      N5B=NY
      N6=L
      ELSEIF(N.EQ.2)THEN
      N4=NX
      N5=NY+1
      N4B=NX
      N5B=NY-1
      N6=L
      ELSEIF(N.EQ.3)THEN
      N4=NX
      N5=NY
      N6=L+1
      ENDIF
C
C     NET WATER, SNOW AND HEAT FLUXES FROM RUNOFF
C
      IF(L.EQ.NUM(N2,N1))THEN
      IF(N.NE.3)THEN
C
C     TQR,TQS,TQW,TQI=net water and snowpack snow,water,ice runoff
C     THQR,THQS=net convective heat from surface water and snow runoff
C     QR,HQR=runoff, convective heat from runoff from watsub.f
C     QS,QW,QI=snow,water,ice transfer from watsub.f
C     HQS=convective heat transfer from snow,water,ice transfer
C        from watsub.f
C     T*QRS=net overland solute flux from runoff
C     X*QRS=solute in flux runoff from trnsfr.f
C     T*QSS=net overland solute flux from snowpack
C     X*QSS=solute flux in snowpack flux from trnsfr.f
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,H1P=HPO4,H2P=H2PO4 in non-band
C        :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,H1B=HPO4,H2B=H2PO4 in band
C
      DO 1202 NN=1,2
      TQR(N2,N1)=TQR(N2,N1)+QR(N,NN,N2,N1)
      THQR(N2,N1)=THQR(N2,N1)+HQR(N,NN,N2,N1)
      TQS(N2,N1)=TQS(N2,N1)+QS(N,NN,N2,N1)
      TQW(N2,N1)=TQW(N2,N1)+QW(N,NN,N2,N1)
      TQI(N2,N1)=TQI(N2,N1)+QI(N,NN,N2,N1)
      THQS(N2,N1)=THQS(N2,N1)+HQS(N,NN,N2,N1)
      DO 8590 K=0,2
      TOCQRS(K,N2,N1)=TOCQRS(K,N2,N1)+XOCQRS(K,N,NN,N2,N1)
      TONQRS(K,N2,N1)=TONQRS(K,N2,N1)+XONQRS(K,N,NN,N2,N1)
      TOPQRS(K,N2,N1)=TOPQRS(K,N2,N1)+XOPQRS(K,N,NN,N2,N1)
      TOAQRS(K,N2,N1)=TOAQRS(K,N2,N1)+XOAQRS(K,N,NN,N2,N1)
8590  CONTINUE
      TCOQRS(N2,N1)=TCOQRS(N2,N1)+XCOQRS(N,NN,N2,N1)
      TCHQRS(N2,N1)=TCHQRS(N2,N1)+XCHQRS(N,NN,N2,N1)
      TOXQRS(N2,N1)=TOXQRS(N2,N1)+XOXQRS(N,NN,N2,N1)
      TNGQRS(N2,N1)=TNGQRS(N2,N1)+XNGQRS(N,NN,N2,N1)
      TN2QRS(N2,N1)=TN2QRS(N2,N1)+XN2QRS(N,NN,N2,N1)
      THGQRS(N2,N1)=THGQRS(N2,N1)+XHGQRS(N,NN,N2,N1)
      TN4QRS(N2,N1)=TN4QRS(N2,N1)+XN4QRW(N,NN,N2,N1)
      TN3QRS(N2,N1)=TN3QRS(N2,N1)+XN3QRW(N,NN,N2,N1)
      TNOQRS(N2,N1)=TNOQRS(N2,N1)+XNOQRW(N,NN,N2,N1)
      TNXQRS(N2,N1)=TNXQRS(N2,N1)+XNXQRS(N,NN,N2,N1)
      TP1QRS(N2,N1)=TP1QRS(N2,N1)+XP1QRW(N,NN,N2,N1)
      TPOQRS(N2,N1)=TPOQRS(N2,N1)+XP4QRW(N,NN,N2,N1)
      TCOQSS(N2,N1)=TCOQSS(N2,N1)+XCOQSS(N,NN,N2,N1)
      TCHQSS(N2,N1)=TCHQSS(N2,N1)+XCHQSS(N,NN,N2,N1)
      TOXQSS(N2,N1)=TOXQSS(N2,N1)+XOXQSS(N,NN,N2,N1)
      TNGQSS(N2,N1)=TNGQSS(N2,N1)+XNGQSS(N,NN,N2,N1)
      TN2QSS(N2,N1)=TN2QSS(N2,N1)+XN2QSS(N,NN,N2,N1)
      TN4QSS(N2,N1)=TN4QSS(N2,N1)+XN4QSS(N,NN,N2,N1)
      TN3QSS(N2,N1)=TN3QSS(N2,N1)+XN3QSS(N,NN,N2,N1)
      TNOQSS(N2,N1)=TNOQSS(N2,N1)+XNOQSS(N,NN,N2,N1)
      TP1QSS(N2,N1)=TP1QSS(N2,N1)+XP1QSS(N,NN,N2,N1)
      TPOQSS(N2,N1)=TPOQSS(N2,N1)+XP4QSS(N,NN,N2,N1)
      IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
      TQR(N2,N1)=TQR(N2,N1)-QR(N,NN,N5,N4)
      THQR(N2,N1)=THQR(N2,N1)-HQR(N,NN,N5,N4)
      TQS(N2,N1)=TQS(N2,N1)-QS(N,NN,N5,N4)
      TQW(N2,N1)=TQW(N2,N1)-QW(N,NN,N5,N4)
      TQI(N2,N1)=TQI(N2,N1)-QI(N,NN,N5,N4)
      THQS(N2,N1)=THQS(N2,N1)-HQS(N,NN,N5,N4)
      DO 8591 K=0,2
      TOCQRS(K,N2,N1)=TOCQRS(K,N2,N1)-XOCQRS(K,N,NN,N5,N4)
      TONQRS(K,N2,N1)=TONQRS(K,N2,N1)-XONQRS(K,N,NN,N5,N4)
      TOPQRS(K,N2,N1)=TOPQRS(K,N2,N1)-XOPQRS(K,N,NN,N5,N4)
      TOAQRS(K,N2,N1)=TOAQRS(K,N2,N1)-XOAQRS(K,N,NN,N5,N4)
8591  CONTINUE
      TCOQRS(N2,N1)=TCOQRS(N2,N1)-XCOQRS(N,NN,N5,N4)
      TCHQRS(N2,N1)=TCHQRS(N2,N1)-XCHQRS(N,NN,N5,N4)
      TOXQRS(N2,N1)=TOXQRS(N2,N1)-XOXQRS(N,NN,N5,N4)
      TNGQRS(N2,N1)=TNGQRS(N2,N1)-XNGQRS(N,NN,N5,N4)
      TN2QRS(N2,N1)=TN2QRS(N2,N1)-XN2QRS(N,NN,N5,N4)
      THGQRS(N2,N1)=THGQRS(N2,N1)-XHGQRS(N,NN,N5,N4)
      TN4QRS(N2,N1)=TN4QRS(N2,N1)-XN4QRW(N,NN,N5,N4)
      TN3QRS(N2,N1)=TN3QRS(N2,N1)-XN3QRW(N,NN,N5,N4)
      TNOQRS(N2,N1)=TNOQRS(N2,N1)-XNOQRW(N,NN,N5,N4)
      TNXQRS(N2,N1)=TNXQRS(N2,N1)-XNXQRS(N,NN,N5,N4)
      TP1QRS(N2,N1)=TP1QRS(N2,N1)-XP1QRW(N,NN,N5,N4)
      TPOQRS(N2,N1)=TPOQRS(N2,N1)-XP4QRW(N,NN,N5,N4)
      ENDIF
C     WRITE(*,6631)'TQR',I,J,N1,N2,N4,N5,N,NN
C    2,IFLBH(N,NN,N5,N4)
C    2,TQR(N2,N1),QR(N,NN,N2,N1),QR(N,NN,N5,N4)
      IF(IFLBHS(N,NN,N5,N4).EQ.0)THEN
      TCOQSS(N2,N1)=TCOQSS(N2,N1)-XCOQSS(N,NN,N5,N4)
      TCHQSS(N2,N1)=TCHQSS(N2,N1)-XCHQSS(N,NN,N5,N4)
      TOXQSS(N2,N1)=TOXQSS(N2,N1)-XOXQSS(N,NN,N5,N4)
      TNGQSS(N2,N1)=TNGQSS(N2,N1)-XNGQSS(N,NN,N5,N4)
      TN2QSS(N2,N1)=TN2QSS(N2,N1)-XN2QSS(N,NN,N5,N4)
      TN4QSS(N2,N1)=TN4QSS(N2,N1)-XN4QSS(N,NN,N5,N4)
      TN3QSS(N2,N1)=TN3QSS(N2,N1)-XN3QSS(N,NN,N5,N4)
      TNOQSS(N2,N1)=TNOQSS(N2,N1)-XNOQSS(N,NN,N5,N4)
      TP1QSS(N2,N1)=TP1QSS(N2,N1)-XP1QSS(N,NN,N5,N4)
      TPOQSS(N2,N1)=TPOQSS(N2,N1)-XP4QSS(N,NN,N5,N4)
      ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      TQR(N2,N1)=TQR(N2,N1)-QR(N,NN,N5B,N4B)
      THQR(N2,N1)=THQR(N2,N1)-HQR(N,NN,N5B,N4B)
      TQS(N2,N1)=TQS(N2,N1)-QS(N,NN,N5B,N4B)
      TQW(N2,N1)=TQW(N2,N1)-QW(N,NN,N5B,N4B)
      TQI(N2,N1)=TQI(N2,N1)-QI(N,NN,N5B,N4B)
      THQS(N2,N1)=THQS(N2,N1)-HQS(N,NN,N5B,N4B)
      DO 8592 K=0,2
      TOCQRS(K,N2,N1)=TOCQRS(K,N2,N1)-XOCQRS(K,N,NN,N5B,N4B)
      TONQRS(K,N2,N1)=TONQRS(K,N2,N1)-XONQRS(K,N,NN,N5B,N4B)
      TOPQRS(K,N2,N1)=TOPQRS(K,N2,N1)-XOPQRS(K,N,NN,N5B,N4B)
      TOAQRS(K,N2,N1)=TOAQRS(K,N2,N1)-XOAQRS(K,N,NN,N5B,N4B)
8592  CONTINUE
      TCOQRS(N2,N1)=TCOQRS(N2,N1)-XCOQRS(N,NN,N5B,N4B)
      TCHQRS(N2,N1)=TCHQRS(N2,N1)-XCHQRS(N,NN,N5B,N4B)
      TOXQRS(N2,N1)=TOXQRS(N2,N1)-XOXQRS(N,NN,N5B,N4B)
      TNGQRS(N2,N1)=TNGQRS(N2,N1)-XNGQRS(N,NN,N5B,N4B)
      TN2QRS(N2,N1)=TN2QRS(N2,N1)-XN2QRS(N,NN,N5B,N4B)
      THGQRS(N2,N1)=THGQRS(N2,N1)-XHGQRS(N,NN,N5B,N4B)
      TN4QRS(N2,N1)=TN4QRS(N2,N1)-XN4QRW(N,NN,N5B,N4B)
      TN3QRS(N2,N1)=TN3QRS(N2,N1)-XN3QRW(N,NN,N5B,N4B)
      TNOQRS(N2,N1)=TNOQRS(N2,N1)-XNOQRW(N,NN,N5B,N4B)
      TNXQRS(N2,N1)=TNXQRS(N2,N1)-XNXQRS(N,NN,N5B,N4B)
      TP1QRS(N2,N1)=TP1QRS(N2,N1)-XP1QRW(N,NN,N5B,N4B)
      TPOQRS(N2,N1)=TPOQRS(N2,N1)-XP4QRW(N,NN,N5B,N4B)
      TCOQSS(N2,N1)=TCOQSS(N2,N1)-XCOQSS(N,NN,N5B,N4B)
      TCHQSS(N2,N1)=TCHQSS(N2,N1)-XCHQSS(N,NN,N5B,N4B)
      TOXQSS(N2,N1)=TOXQSS(N2,N1)-XOXQSS(N,NN,N5B,N4B)
      TNGQSS(N2,N1)=TNGQSS(N2,N1)-XNGQSS(N,NN,N5B,N4B)
      TN2QSS(N2,N1)=TN2QSS(N2,N1)-XN2QSS(N,NN,N5B,N4B)
      TN4QSS(N2,N1)=TN4QSS(N2,N1)-XN4QSS(N,NN,N5B,N4B)
      TN3QSS(N2,N1)=TN3QSS(N2,N1)-XN3QSS(N,NN,N5B,N4B)
      TNOQSS(N2,N1)=TNOQSS(N2,N1)-XNOQSS(N,NN,N5B,N4B)
      TP1QSS(N2,N1)=TP1QSS(N2,N1)-XP1QSS(N,NN,N5B,N4B)
      TPOQSS(N2,N1)=TPOQSS(N2,N1)-XP4QSS(N,NN,N5B,N4B)
C     WRITE(*,6631)'TQRB',I,J,N1,N2,N4B,N5B,N,NN
C    2,IFLBH(N,NN,N5B,N4B)
C    2,TQR(N2,N1),QR(N,NN,N5B,N4B)
6631  FORMAT(A8,9I4,12E12.4)
      ENDIF
1202  CONTINUE
C
C     NET SALT FLUXES FROM RUNOFF AND SNOWPACK
C
C     TQR*=net overland solute flux in runoff
C     XQR*=solute in runoff from trnsfrs.f
C     TQS*=net overland solute flux in snow drift
C     XQS*=solute in snow drift from trnsfrs.f
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C        :*1=non-band,*B=band
C
      IF(ISALTG.NE.0)THEN
      DO 1203 NN=1,2
      TQRAL(N2,N1)=TQRAL(N2,N1)+XQRAL(N,NN,N2,N1)
      TQRFE(N2,N1)=TQRFE(N2,N1)+XQRFE(N,NN,N2,N1)
      TQRHY(N2,N1)=TQRHY(N2,N1)+XQRHY(N,NN,N2,N1)
      TQRCA(N2,N1)=TQRCA(N2,N1)+XQRCA(N,NN,N2,N1)
      TQRMG(N2,N1)=TQRMG(N2,N1)+XQRMG(N,NN,N2,N1)
      TQRNA(N2,N1)=TQRNA(N2,N1)+XQRNA(N,NN,N2,N1)
      TQRKA(N2,N1)=TQRKA(N2,N1)+XQRKA(N,NN,N2,N1)
      TQROH(N2,N1)=TQROH(N2,N1)+XQROH(N,NN,N2,N1)
      TQRSO(N2,N1)=TQRSO(N2,N1)+XQRSO(N,NN,N2,N1)
      TQRCL(N2,N1)=TQRCL(N2,N1)+XQRCL(N,NN,N2,N1)
      TQRC3(N2,N1)=TQRC3(N2,N1)+XQRC3(N,NN,N2,N1)
      TQRHC(N2,N1)=TQRHC(N2,N1)+XQRHC(N,NN,N2,N1)
      TQRAL1(N2,N1)=TQRAL1(N2,N1)+XQRAL1(N,NN,N2,N1)
      TQRAL2(N2,N1)=TQRAL2(N2,N1)+XQRAL2(N,NN,N2,N1)
      TQRAL3(N2,N1)=TQRAL3(N2,N1)+XQRAL3(N,NN,N2,N1)
      TQRAL4(N2,N1)=TQRAL4(N2,N1)+XQRAL4(N,NN,N2,N1)
      TQRALS(N2,N1)=TQRALS(N2,N1)+XQRALS(N,NN,N2,N1)
      TQRFE1(N2,N1)=TQRFE1(N2,N1)+XQRFE1(N,NN,N2,N1)
      TQRFE2(N2,N1)=TQRFE2(N2,N1)+XQRFE2(N,NN,N2,N1)
      TQRFE3(N2,N1)=TQRFE3(N2,N1)+XQRFE3(N,NN,N2,N1)
      TQRFE4(N2,N1)=TQRFE4(N2,N1)+XQRFE4(N,NN,N2,N1)
      TQRFES(N2,N1)=TQRFES(N2,N1)+XQRFES(N,NN,N2,N1)
      TQRCAO(N2,N1)=TQRCAO(N2,N1)+XQRCAO(N,NN,N2,N1)
      TQRCAC(N2,N1)=TQRCAC(N2,N1)+XQRCAC(N,NN,N2,N1)
      TQRCAH(N2,N1)=TQRCAH(N2,N1)+XQRCAH(N,NN,N2,N1)
      TQRCAS(N2,N1)=TQRCAS(N2,N1)+XQRCAS(N,NN,N2,N1)
      TQRMGO(N2,N1)=TQRMGO(N2,N1)+XQRMGO(N,NN,N2,N1)
      TQRMGC(N2,N1)=TQRMGC(N2,N1)+XQRMGC(N,NN,N2,N1)
      TQRMGH(N2,N1)=TQRMGH(N2,N1)+XQRMGH(N,NN,N2,N1)
      TQRMGS(N2,N1)=TQRMGS(N2,N1)+XQRMGS(N,NN,N2,N1)
      TQRNAC(N2,N1)=TQRNAC(N2,N1)+XQRNAC(N,NN,N2,N1)
      TQRNAS(N2,N1)=TQRNAS(N2,N1)+XQRNAS(N,NN,N2,N1)
      TQRKAS(N2,N1)=TQRKAS(N2,N1)+XQRKAS(N,NN,N2,N1)
      TQRH0P(N2,N1)=TQRH0P(N2,N1)+XQRH0P(N,NN,N2,N1)
      TQRH3P(N2,N1)=TQRH3P(N2,N1)+XQRH3P(N,NN,N2,N1)
      TQRF1P(N2,N1)=TQRF1P(N2,N1)+XQRF1P(N,NN,N2,N1)
      TQRF2P(N2,N1)=TQRF2P(N2,N1)+XQRF2P(N,NN,N2,N1)
      TQRC0P(N2,N1)=TQRC0P(N2,N1)+XQRC0P(N,NN,N2,N1)
      TQRC1P(N2,N1)=TQRC1P(N2,N1)+XQRC1P(N,NN,N2,N1)
      TQRC2P(N2,N1)=TQRC2P(N2,N1)+XQRC2P(N,NN,N2,N1)
      TQRM1P(N2,N1)=TQRM1P(N2,N1)+XQRM1P(N,NN,N2,N1)
      TQSAL(N2,N1)=TQSAL(N2,N1)+XQSAL(N,NN,N2,N1)
      TQSFE(N2,N1)=TQSFE(N2,N1)+XQSFE(N,NN,N2,N1)
      TQSHY(N2,N1)=TQSHY(N2,N1)+XQSHY(N,NN,N2,N1)
      TQSCA(N2,N1)=TQSCA(N2,N1)+XQSCA(N,NN,N2,N1)
      TQSMG(N2,N1)=TQSMG(N2,N1)+XQSMG(N,NN,N2,N1)
      TQSNA(N2,N1)=TQSNA(N2,N1)+XQSNA(N,NN,N2,N1)
      TQSKA(N2,N1)=TQSKA(N2,N1)+XQSKA(N,NN,N2,N1)
      TQSOH(N2,N1)=TQSOH(N2,N1)+XQSOH(N,NN,N2,N1)
      TQSSO(N2,N1)=TQSSO(N2,N1)+XQSSO(N,NN,N2,N1)
      TQSCL(N2,N1)=TQSCL(N2,N1)+XQSCL(N,NN,N2,N1)
      TQSC3(N2,N1)=TQSC3(N2,N1)+XQSC3(N,NN,N2,N1)
      TQSHC(N2,N1)=TQSHC(N2,N1)+XQSHC(N,NN,N2,N1)
      TQSAL1(N2,N1)=TQSAL1(N2,N1)+XQSAL1(N,NN,N2,N1)
      TQSAL2(N2,N1)=TQSAL2(N2,N1)+XQSAL2(N,NN,N2,N1)
      TQSAL3(N2,N1)=TQSAL3(N2,N1)+XQSAL3(N,NN,N2,N1)
      TQSAL4(N2,N1)=TQSAL4(N2,N1)+XQSAL4(N,NN,N2,N1)
      TQSALS(N2,N1)=TQSALS(N2,N1)+XQSALS(N,NN,N2,N1)
      TQSFE1(N2,N1)=TQSFE1(N2,N1)+XQSFE1(N,NN,N2,N1)
      TQSFE2(N2,N1)=TQSFE2(N2,N1)+XQSFE2(N,NN,N2,N1)
      TQSFE3(N2,N1)=TQSFE3(N2,N1)+XQSFE3(N,NN,N2,N1)
      TQSFE4(N2,N1)=TQSFE4(N2,N1)+XQSFE4(N,NN,N2,N1)
      TQSFES(N2,N1)=TQSFES(N2,N1)+XQSFES(N,NN,N2,N1)
      TQSCAO(N2,N1)=TQSCAO(N2,N1)+XQSCAO(N,NN,N2,N1)
      TQSCAC(N2,N1)=TQSCAC(N2,N1)+XQSCAC(N,NN,N2,N1)
      TQSCAH(N2,N1)=TQSCAH(N2,N1)+XQSCAH(N,NN,N2,N1)
      TQSCAS(N2,N1)=TQSCAS(N2,N1)+XQSCAS(N,NN,N2,N1)
      TQSMGO(N2,N1)=TQSMGO(N2,N1)+XQSMGO(N,NN,N2,N1)
      TQSMGC(N2,N1)=TQSMGC(N2,N1)+XQSMGC(N,NN,N2,N1)
      TQSMGH(N2,N1)=TQSMGH(N2,N1)+XQSMGH(N,NN,N2,N1)
      TQSMGS(N2,N1)=TQSMGS(N2,N1)+XQSMGS(N,NN,N2,N1)
      TQSNAC(N2,N1)=TQSNAC(N2,N1)+XQSNAC(N,NN,N2,N1)
      TQSNAS(N2,N1)=TQSNAS(N2,N1)+XQSNAS(N,NN,N2,N1)
      TQSKAS(N2,N1)=TQSKAS(N2,N1)+XQSKAS(N,NN,N2,N1)
      TQSH0P(N2,N1)=TQSH0P(N2,N1)+XQSH0P(N,NN,N2,N1)
      TQSH3P(N2,N1)=TQSH3P(N2,N1)+XQSH3P(N,NN,N2,N1)
      TQSF1P(N2,N1)=TQSF1P(N2,N1)+XQSF1P(N,NN,N2,N1)
      TQSF2P(N2,N1)=TQSF2P(N2,N1)+XQSF2P(N,NN,N2,N1)
      TQSC0P(N2,N1)=TQSC0P(N2,N1)+XQSC0P(N,NN,N2,N1)
      TQSC1P(N2,N1)=TQSC1P(N2,N1)+XQSC1P(N,NN,N2,N1)
      TQSC2P(N2,N1)=TQSC2P(N2,N1)+XQSC2P(N,NN,N2,N1)
      TQSM1P(N2,N1)=TQSM1P(N2,N1)+XQSM1P(N,NN,N2,N1)
      IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
      TQRAL(N2,N1)=TQRAL(N2,N1)-XQRAL(N,NN,N5,N4)
      TQRFE(N2,N1)=TQRFE(N2,N1)-XQRFE(N,NN,N5,N4)
      TQRHY(N2,N1)=TQRHY(N2,N1)-XQRHY(N,NN,N5,N4)
      TQRCA(N2,N1)=TQRCA(N2,N1)-XQRCA(N,NN,N5,N4)
      TQRMG(N2,N1)=TQRMG(N2,N1)-XQRMG(N,NN,N5,N4)
      TQRNA(N2,N1)=TQRNA(N2,N1)-XQRNA(N,NN,N5,N4)
      TQRKA(N2,N1)=TQRKA(N2,N1)-XQRKA(N,NN,N5,N4)
      TQROH(N2,N1)=TQROH(N2,N1)-XQROH(N,NN,N5,N4)
      TQRSO(N2,N1)=TQRSO(N2,N1)-XQRSO(N,NN,N5,N4)
      TQRCL(N2,N1)=TQRCL(N2,N1)-XQRCL(N,NN,N5,N4)
      TQRC3(N2,N1)=TQRC3(N2,N1)-XQRC3(N,NN,N5,N4)
      TQRHC(N2,N1)=TQRHC(N2,N1)-XQRHC(N,NN,N5,N4)
      TQRAL1(N2,N1)=TQRAL1(N2,N1)+XQRAL1(N,NN,N2,N1)
      TQRAL2(N2,N1)=TQRAL2(N2,N1)+XQRAL2(N,NN,N2,N1)
      TQRAL3(N2,N1)=TQRAL3(N2,N1)+XQRAL3(N,NN,N2,N1)
      TQRAL4(N2,N1)=TQRAL4(N2,N1)+XQRAL4(N,NN,N2,N1)
      TQRALS(N2,N1)=TQRALS(N2,N1)+XQRALS(N,NN,N2,N1)
      TQRFE1(N2,N1)=TQRFE1(N2,N1)+XQRFE1(N,NN,N2,N1)
      TQRFE2(N2,N1)=TQRFE2(N2,N1)+XQRFE2(N,NN,N2,N1)
      TQRFE3(N2,N1)=TQRFE3(N2,N1)+XQRFE3(N,NN,N2,N1)
      TQRFE4(N2,N1)=TQRFE4(N2,N1)+XQRFE4(N,NN,N2,N1)
      TQRFES(N2,N1)=TQRFES(N2,N1)+XQRFES(N,NN,N2,N1)
      TQRCAO(N2,N1)=TQRCAO(N2,N1)+XQRCAO(N,NN,N2,N1)
      TQRCAC(N2,N1)=TQRCAC(N2,N1)+XQRCAC(N,NN,N2,N1)
      TQRCAH(N2,N1)=TQRCAH(N2,N1)+XQRCAH(N,NN,N2,N1)
      TQRCAS(N2,N1)=TQRCAS(N2,N1)+XQRCAS(N,NN,N2,N1)
      TQRMGO(N2,N1)=TQRMGO(N2,N1)+XQRMGO(N,NN,N2,N1)
      TQRMGC(N2,N1)=TQRMGC(N2,N1)+XQRMGC(N,NN,N2,N1)
      TQRMGH(N2,N1)=TQRMGH(N2,N1)+XQRMGH(N,NN,N2,N1)
      TQRMGS(N2,N1)=TQRMGS(N2,N1)+XQRMGS(N,NN,N2,N1)
      TQRNAC(N2,N1)=TQRNAC(N2,N1)+XQRNAC(N,NN,N2,N1)
      TQRNAS(N2,N1)=TQRNAS(N2,N1)+XQRNAS(N,NN,N2,N1)
      TQRKAS(N2,N1)=TQRKAS(N2,N1)+XQRKAS(N,NN,N2,N1)
      TQRH0P(N2,N1)=TQRH0P(N2,N1)+XQRH0P(N,NN,N2,N1)
      TQRH3P(N2,N1)=TQRH3P(N2,N1)+XQRH3P(N,NN,N2,N1)
      TQRF1P(N2,N1)=TQRF1P(N2,N1)+XQRF1P(N,NN,N2,N1)
      TQRF2P(N2,N1)=TQRF2P(N2,N1)+XQRF2P(N,NN,N2,N1)
      TQRC0P(N2,N1)=TQRC0P(N2,N1)+XQRC0P(N,NN,N2,N1)
      TQRC1P(N2,N1)=TQRC1P(N2,N1)+XQRC1P(N,NN,N2,N1)
      TQRC2P(N2,N1)=TQRC2P(N2,N1)+XQRC2P(N,NN,N2,N1)
      TQRM1P(N2,N1)=TQRM1P(N2,N1)+XQRM1P(N,NN,N2,N1)
      ENDIF
      IF(IFLBHS(N,NN,N5,N4).EQ.0)THEN
      TQSAL(N2,N1)=TQSAL(N2,N1)-XQSAL(N,NN,N5,N4)
      TQSFE(N2,N1)=TQSFE(N2,N1)-XQSFE(N,NN,N5,N4)
      TQSHY(N2,N1)=TQSHY(N2,N1)-XQSHY(N,NN,N5,N4)
      TQSCA(N2,N1)=TQSCA(N2,N1)-XQSCA(N,NN,N5,N4)
      TQSMG(N2,N1)=TQSMG(N2,N1)-XQSMG(N,NN,N5,N4)
      TQSNA(N2,N1)=TQSNA(N2,N1)-XQSNA(N,NN,N5,N4)
      TQSKA(N2,N1)=TQSKA(N2,N1)-XQSKA(N,NN,N5,N4)
      TQSOH(N2,N1)=TQSOH(N2,N1)-XQSOH(N,NN,N5,N4)
      TQSSO(N2,N1)=TQSSO(N2,N1)-XQSSO(N,NN,N5,N4)
      TQSCL(N2,N1)=TQSCL(N2,N1)-XQSCL(N,NN,N5,N4)
      TQSC3(N2,N1)=TQSC3(N2,N1)-XQSC3(N,NN,N5,N4)
      TQSHC(N2,N1)=TQSHC(N2,N1)-XQSHC(N,NN,N5,N4)
      TQSAL1(N2,N1)=TQSAL1(N2,N1)-XQSAL1(N,NN,N5,N4)
      TQSAL2(N2,N1)=TQSAL2(N2,N1)-XQSAL2(N,NN,N5,N4)
      TQSAL3(N2,N1)=TQSAL3(N2,N1)-XQSAL3(N,NN,N5,N4)
      TQSAL4(N2,N1)=TQSAL4(N2,N1)-XQSAL4(N,NN,N5,N4)
      TQSALS(N2,N1)=TQSALS(N2,N1)-XQSALS(N,NN,N5,N4)
      TQSFE1(N2,N1)=TQSFE1(N2,N1)-XQSFE1(N,NN,N5,N4)
      TQSFE2(N2,N1)=TQSFE2(N2,N1)-XQSFE2(N,NN,N5,N4)
      TQSFE3(N2,N1)=TQSFE3(N2,N1)-XQSFE3(N,NN,N5,N4)
      TQSFE4(N2,N1)=TQSFE4(N2,N1)-XQSFE4(N,NN,N5,N4)
      TQSFES(N2,N1)=TQSFES(N2,N1)-XQSFES(N,NN,N5,N4)
      TQSCAO(N2,N1)=TQSCAO(N2,N1)-XQSCAO(N,NN,N5,N4)
      TQSCAC(N2,N1)=TQSCAC(N2,N1)-XQSCAC(N,NN,N5,N4)
      TQSCAH(N2,N1)=TQSCAH(N2,N1)-XQSCAH(N,NN,N5,N4)
      TQSCAS(N2,N1)=TQSCAS(N2,N1)-XQSCAS(N,NN,N5,N4)
      TQSMGO(N2,N1)=TQSMGO(N2,N1)-XQSMGO(N,NN,N5,N4)
      TQSMGC(N2,N1)=TQSMGC(N2,N1)-XQSMGC(N,NN,N5,N4)
      TQSMGH(N2,N1)=TQSMGH(N2,N1)-XQSMGH(N,NN,N5,N4)
      TQSMGS(N2,N1)=TQSMGS(N2,N1)-XQSMGS(N,NN,N5,N4)
      TQSNAC(N2,N1)=TQSNAC(N2,N1)-XQSNAC(N,NN,N5,N4)
      TQSNAS(N2,N1)=TQSNAS(N2,N1)-XQSNAS(N,NN,N5,N4)
      TQSKAS(N2,N1)=TQSKAS(N2,N1)-XQSKAS(N,NN,N5,N4)
      TQSH0P(N2,N1)=TQSH0P(N2,N1)-XQSH0P(N,NN,N5,N4)
      TQSH3P(N2,N1)=TQSH3P(N2,N1)-XQSH3P(N,NN,N5,N4)
      TQSF1P(N2,N1)=TQSF1P(N2,N1)-XQSF1P(N,NN,N5,N4)
      TQSF2P(N2,N1)=TQSF2P(N2,N1)-XQSF2P(N,NN,N5,N4)
      TQSC0P(N2,N1)=TQSC0P(N2,N1)-XQSC0P(N,NN,N5,N4)
      TQSC1P(N2,N1)=TQSC1P(N2,N1)-XQSC1P(N,NN,N5,N4)
      TQSC2P(N2,N1)=TQSC2P(N2,N1)-XQSC2P(N,NN,N5,N4)
      TQSM1P(N2,N1)=TQSM1P(N2,N1)-XQSM1P(N,NN,N5,N4)
      ENDIF
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      TQRAL(N2,N1)=TQRAL(N2,N1)-XQRAL(N,NN,N5B,N4B)
      TQRFE(N2,N1)=TQRFE(N2,N1)-XQRFE(N,NN,N5B,N4B)
      TQRHY(N2,N1)=TQRHY(N2,N1)-XQRHY(N,NN,N5B,N4B)
      TQRCA(N2,N1)=TQRCA(N2,N1)-XQRCA(N,NN,N5B,N4B)
      TQRMG(N2,N1)=TQRMG(N2,N1)-XQRMG(N,NN,N5B,N4B)
      TQRNA(N2,N1)=TQRNA(N2,N1)-XQRNA(N,NN,N5B,N4B)
      TQRKA(N2,N1)=TQRKA(N2,N1)-XQRKA(N,NN,N5B,N4B)
      TQROH(N2,N1)=TQROH(N2,N1)-XQROH(N,NN,N5B,N4B)
      TQRSO(N2,N1)=TQRSO(N2,N1)-XQRSO(N,NN,N5B,N4B)
      TQRCL(N2,N1)=TQRCL(N2,N1)-XQRCL(N,NN,N5B,N4B)
      TQRC3(N2,N1)=TQRC3(N2,N1)-XQRC3(N,NN,N5B,N4B)
      TQRHC(N2,N1)=TQRHC(N2,N1)-XQRHC(N,NN,N5B,N4B)
      TQRAL1(N2,N1)=TQRAL1(N2,N1)-XQRAL1(N,NN,N5B,N4B)
      TQRAL2(N2,N1)=TQRAL2(N2,N1)-XQRAL2(N,NN,N5B,N4B)
      TQRAL3(N2,N1)=TQRAL3(N2,N1)-XQRAL3(N,NN,N5B,N4B)
      TQRAL4(N2,N1)=TQRAL4(N2,N1)-XQRAL4(N,NN,N5B,N4B)
      TQRALS(N2,N1)=TQRALS(N2,N1)-XQRALS(N,NN,N5B,N4B)
      TQRFE1(N2,N1)=TQRFE1(N2,N1)-XQRFE1(N,NN,N5B,N4B)
      TQRFE2(N2,N1)=TQRFE2(N2,N1)-XQRFE2(N,NN,N5B,N4B)
      TQRFE3(N2,N1)=TQRFE3(N2,N1)-XQRFE3(N,NN,N5B,N4B)
      TQRFE4(N2,N1)=TQRFE4(N2,N1)-XQRFE4(N,NN,N5B,N4B)
      TQRFES(N2,N1)=TQRFES(N2,N1)-XQRFES(N,NN,N5B,N4B)
      TQRCAO(N2,N1)=TQRCAO(N2,N1)-XQRCAO(N,NN,N5B,N4B)
      TQRCAC(N2,N1)=TQRCAC(N2,N1)-XQRCAC(N,NN,N5B,N4B)
      TQRCAH(N2,N1)=TQRCAH(N2,N1)-XQRCAH(N,NN,N5B,N4B)
      TQRCAS(N2,N1)=TQRCAS(N2,N1)-XQRCAS(N,NN,N5B,N4B)
      TQRMGO(N2,N1)=TQRMGO(N2,N1)-XQRMGO(N,NN,N5B,N4B)
      TQRMGC(N2,N1)=TQRMGC(N2,N1)-XQRMGC(N,NN,N5B,N4B)
      TQRMGH(N2,N1)=TQRMGH(N2,N1)-XQRMGH(N,NN,N5B,N4B)
      TQRMGS(N2,N1)=TQRMGS(N2,N1)-XQRMGS(N,NN,N5B,N4B)
      TQRNAC(N2,N1)=TQRNAC(N2,N1)-XQRNAC(N,NN,N5B,N4B)
      TQRNAS(N2,N1)=TQRNAS(N2,N1)-XQRNAS(N,NN,N5B,N4B)
      TQRKAS(N2,N1)=TQRKAS(N2,N1)-XQRKAS(N,NN,N5B,N4B)
      TQRH0P(N2,N1)=TQRH0P(N2,N1)-XQRH0P(N,NN,N5B,N4B)
      TQRH3P(N2,N1)=TQRH3P(N2,N1)-XQRH3P(N,NN,N5B,N4B)
      TQRF1P(N2,N1)=TQRF1P(N2,N1)-XQRF1P(N,NN,N5B,N4B)
      TQRF2P(N2,N1)=TQRF2P(N2,N1)-XQRF2P(N,NN,N5B,N4B)
      TQRC0P(N2,N1)=TQRC0P(N2,N1)-XQRC0P(N,NN,N5B,N4B)
      TQRC1P(N2,N1)=TQRC1P(N2,N1)-XQRC1P(N,NN,N5B,N4B)
      TQRC2P(N2,N1)=TQRC2P(N2,N1)-XQRC2P(N,NN,N5B,N4B)
      TQRM1P(N2,N1)=TQRM1P(N2,N1)-XQRM1P(N,NN,N5B,N4B)
      TQSAL(N2,N1)=TQSAL(N2,N1)-XQSAL(N,NN,N5B,N4B)
      TQSFE(N2,N1)=TQSFE(N2,N1)-XQSFE(N,NN,N5B,N4B)
      TQSHY(N2,N1)=TQSHY(N2,N1)-XQSHY(N,NN,N5B,N4B)
      TQSCA(N2,N1)=TQSCA(N2,N1)-XQSCA(N,NN,N5B,N4B)
      TQSMG(N2,N1)=TQSMG(N2,N1)-XQSMG(N,NN,N5B,N4B)
      TQSNA(N2,N1)=TQSNA(N2,N1)-XQSNA(N,NN,N5B,N4B)
      TQSKA(N2,N1)=TQSKA(N2,N1)-XQSKA(N,NN,N5B,N4B)
      TQSOH(N2,N1)=TQSOH(N2,N1)-XQSOH(N,NN,N5B,N4B)
      TQSSO(N2,N1)=TQSSO(N2,N1)-XQSSO(N,NN,N5B,N4B)
      TQSCL(N2,N1)=TQSCL(N2,N1)-XQSCL(N,NN,N5B,N4B)
      TQSC3(N2,N1)=TQSC3(N2,N1)-XQSC3(N,NN,N5B,N4B)
      TQSHC(N2,N1)=TQSHC(N2,N1)-XQSHC(N,NN,N5B,N4B)
      TQSAL1(N2,N1)=TQSAL1(N2,N1)-XQSAL1(N,NN,N5B,N4B)
      TQSAL2(N2,N1)=TQSAL2(N2,N1)-XQSAL2(N,NN,N5B,N4B)
      TQSAL3(N2,N1)=TQSAL3(N2,N1)-XQSAL3(N,NN,N5B,N4B)
      TQSAL4(N2,N1)=TQSAL4(N2,N1)-XQSAL4(N,NN,N5B,N4B)
      TQSALS(N2,N1)=TQSALS(N2,N1)-XQSALS(N,NN,N5B,N4B)
      TQSFE1(N2,N1)=TQSFE1(N2,N1)-XQSFE1(N,NN,N5B,N4B)
      TQSFE2(N2,N1)=TQSFE2(N2,N1)-XQSFE2(N,NN,N5B,N4B)
      TQSFE3(N2,N1)=TQSFE3(N2,N1)-XQSFE3(N,NN,N5B,N4B)
      TQSFE4(N2,N1)=TQSFE4(N2,N1)-XQSFE4(N,NN,N5B,N4B)
      TQSFES(N2,N1)=TQSFES(N2,N1)-XQSFES(N,NN,N5B,N4B)
      TQSCAO(N2,N1)=TQSCAO(N2,N1)-XQSCAO(N,NN,N5B,N4B)
      TQSCAC(N2,N1)=TQSCAC(N2,N1)-XQSCAC(N,NN,N5B,N4B)
      TQSCAH(N2,N1)=TQSCAH(N2,N1)-XQSCAH(N,NN,N5B,N4B)
      TQSCAS(N2,N1)=TQSCAS(N2,N1)-XQSCAS(N,NN,N5B,N4B)
      TQSMGO(N2,N1)=TQSMGO(N2,N1)-XQSMGO(N,NN,N5B,N4B)
      TQSMGC(N2,N1)=TQSMGC(N2,N1)-XQSMGC(N,NN,N5B,N4B)
      TQSMGH(N2,N1)=TQSMGH(N2,N1)-XQSMGH(N,NN,N5B,N4B)
      TQSMGS(N2,N1)=TQSMGS(N2,N1)-XQSMGS(N,NN,N5B,N4B)
      TQSNAC(N2,N1)=TQSNAC(N2,N1)-XQSNAC(N,NN,N5B,N4B)
      TQSNAS(N2,N1)=TQSNAS(N2,N1)-XQSNAS(N,NN,N5B,N4B)
      TQSKAS(N2,N1)=TQSKAS(N2,N1)-XQSKAS(N,NN,N5B,N4B)
      TQSH0P(N2,N1)=TQSH0P(N2,N1)-XQSH0P(N,NN,N5B,N4B)
      TQSH3P(N2,N1)=TQSH3P(N2,N1)-XQSH3P(N,NN,N5B,N4B)
      TQSF1P(N2,N1)=TQSF1P(N2,N1)-XQSF1P(N,NN,N5B,N4B)
      TQSF2P(N2,N1)=TQSF2P(N2,N1)-XQSF2P(N,NN,N5B,N4B)
      TQSC0P(N2,N1)=TQSC0P(N2,N1)-XQSC0P(N,NN,N5B,N4B)
      TQSC1P(N2,N1)=TQSC1P(N2,N1)-XQSC1P(N,NN,N5B,N4B)
      TQSC2P(N2,N1)=TQSC2P(N2,N1)-XQSC2P(N,NN,N5B,N4B)
      TQSM1P(N2,N1)=TQSM1P(N2,N1)-XQSM1P(N,NN,N5B,N4B)
      ENDIF
1203  CONTINUE
      ENDIF
C
C     NET WATER AND HEAT FLUXES THROUGH SNOWPACK
C
C     VHCPW,VHCPWX=current, minimum snowpack heat capacities
C     TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
C     THFLWW=convective heat fluxes of snow,water,ice in snowpack
C     XFLWS,XFLWW,XFLWI=snow,water,ice transfer from watsub.f
C     XHFLWW=convective heat flux from snow,water,ice transfer
C        from watsub.f
C     FLSW,FLSWH,FLSWR=water flux from lowest snow layer to soil
C        macropore,micropore,litter
C     HFLSW,HFLSWR=heat flux from lowest snow layer to soil,litter
C
      ELSEIF(N.EQ.3)THEN
      DO 1205 LS=1,JS
      IF(VHCPW(LS,NY,NX).GT.VHCPWX(NY,NX))THEN
      LS2=MIN(JS,LS+1)
C
C     IF LOWER LAYER IS IN THE SNOWPACK
C
      IF(LS.LT.JS.AND.VHCPW(LS2,N2,N1).GT.VHCPWX(N2,N1))THEN
      TFLWS(LS,N2,N1)=TFLWS(LS,N2,N1)+XFLWS(LS,N2,N1)
     2-XFLWS(LS2,N2,N1)
      TFLWW(LS,N2,N1)=TFLWW(LS,N2,N1)+XFLWW(LS,N2,N1)
     2-XFLWW(LS2,N2,N1)
     2-FLSWR(LS,N2,N1)-FLSW(LS,N2,N1)-FLSWH(LS,N2,N1)
      TFLWV(LS,N2,N1)=TFLWV(LS,N2,N1)+XFLWV(LS,N2,N1)
     2-XFLWV(LS2,N2,N1)
     2-FLSVR(LS,N2,N1)-FLSV(LS,N2,N1)
      TFLWI(LS,N2,N1)=TFLWI(LS,N2,N1)+XFLWI(LS,N2,N1)
     2-XFLWI(LS2,N2,N1)
      THFLWW(LS,N2,N1)=THFLWW(LS,N2,N1)+XHFLWW(LS,N2,N1)
     2-XHFLWW(LS2,N2,N1)
     2-HFLSWR(LS,N2,N1)-HFLSW(LS,N2,N1)
C     IF(LS.EQ.5)THEN
C     WRITE(*,7754)'LS',I,J,N1,N2,LS,LS2,TFLWW(LS,N2,N1),XFLWW(LS,N2,N1)
C    2,XFLWW(LS2,N2,N1),TFLWV(LS,N2,N1),XFLWV(LS,N2,N1)
C    2,XFLWV(LS2,N2,N1),FLSWR(LS,N2,N1),FLSW(LS,N2,N1),FLSWH(LS,N2,N1)
C    2,FLSVR(LS,N2,N1),FLSV(LS,N2,N1)
7754  FORMAT(A8,6I4,100E14.6)
C     ENDIF
C
C     NET SOLUTE FLUXES THROUGH SNOWPACK
C
C     T*BLS=net solute flux in snowpack
C     X*BLS=solute flux in snowpack from trnsfr.f
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,H1P=HPO4,H2P=H2PO4 in non-band
C        :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,H1B=HPO4,H2B=H2PO4 in band
C
      TCOBLS(LS,N2,N1)=TCOBLS(LS,N2,N1)+XCOBLS(LS,N2,N1)
     2-XCOBLS(LS2,N2,N1)
      TCHBLS(LS,N2,N1)=TCHBLS(LS,N2,N1)+XCHBLS(LS,N2,N1)
     2-XCHBLS(LS2,N2,N1)
      TOXBLS(LS,N2,N1)=TOXBLS(LS,N2,N1)+XOXBLS(LS,N2,N1)
     2-XOXBLS(LS2,N2,N1)
      TNGBLS(LS,N2,N1)=TNGBLS(LS,N2,N1)+XNGBLS(LS,N2,N1)
     2-XNGBLS(LS2,N2,N1)
      TN2BLS(LS,N2,N1)=TN2BLS(LS,N2,N1)+XN2BLS(LS,N2,N1)
     2-XN2BLS(LS2,N2,N1)
      TN4BLW(LS,N2,N1)=TN4BLW(LS,N2,N1)+XN4BLW(LS,N2,N1)
     2-XN4BLW(LS2,N2,N1)
      TN3BLW(LS,N2,N1)=TN3BLW(LS,N2,N1)+XN3BLW(LS,N2,N1)
     2-XN3BLW(LS2,N2,N1)
      TNOBLW(LS,N2,N1)=TNOBLW(LS,N2,N1)+XNOBLW(LS,N2,N1)
     2-XNOBLW(LS2,N2,N1)
      TH1PBS(LS,N2,N1)=TH1PBS(LS,N2,N1)+XH1PBS(LS,N2,N1)
     2-XH1PBS(LS2,N2,N1)
      TH2PBS(LS,N2,N1)=TH2PBS(LS,N2,N1)+XH2PBS(LS,N2,N1)
     2-XH2PBS(LS2,N2,N1)
C
C     NET SALT FLUXES THROUGH SNOWPACK
C
C     T*BLS=net solute flux in snowpack
C     X*BLS=solute flux in snowpack from trnsfrs.f
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C        :*1=non-band,*B=band
C
      IF(ISALTG.NE.0)THEN
      TALBLS(LS,N2,N1)=TALBLS(LS,N2,N1)+XALBLS(LS,N2,N1)
     2-XALBLS(LS2,N2,N1)
      TFEBLS(LS,N2,N1)=TFEBLS(LS,N2,N1)+XFEBLS(LS,N2,N1)
     2-XFEBLS(LS2,N2,N1)
      THYBLS(LS,N2,N1)=THYBLS(LS,N2,N1)+XHYBLS(LS,N2,N1)
     2-XHYBLS(LS2,N2,N1)
      TCABLS(LS,N2,N1)=TCABLS(LS,N2,N1)+XCABLS(LS,N2,N1)
     2-XCABLS(LS2,N2,N1)
      TMGBLS(LS,N2,N1)=TMGBLS(LS,N2,N1)+XMGBLS(LS,N2,N1)
     2-XMGBLS(LS2,N2,N1)
      TNABLS(LS,N2,N1)=TNABLS(LS,N2,N1)+XNABLS(LS,N2,N1)
     2-XNABLS(LS2,N2,N1)
      TKABLS(LS,N2,N1)=TKABLS(LS,N2,N1)+XKABLS(LS,N2,N1)
     2-XKABLS(LS2,N2,N1)
      TOHBLS(LS,N2,N1)=TOHBLS(LS,N2,N1)+XOHBLS(LS,N2,N1)
     2-XOHBLS(LS2,N2,N1)
      TSOBLS(LS,N2,N1)=TSOBLS(LS,N2,N1)+XSOBLS(LS,N2,N1)
     2-XSOBLS(LS2,N2,N1)
      TCLBLS(LS,N2,N1)=TCLBLS(LS,N2,N1)+XCLBLS(LS,N2,N1)
     2-XCLBLS(LS2,N2,N1)
      TC3BLS(LS,N2,N1)=TC3BLS(LS,N2,N1)+XC3BLS(LS,N2,N1)
     2-XC3BLS(LS2,N2,N1)
      THCBLS(LS,N2,N1)=THCBLS(LS,N2,N1)+XHCBLS(LS,N2,N1)
     2-XHCBLS(LS2,N2,N1)
      TAL1BS(LS,N2,N1)=TAL1BS(LS,N2,N1)+XAL1BS(LS,N2,N1)
     2-XAL1BS(LS2,N2,N1)
      TAL2BS(LS,N2,N1)=TAL2BS(LS,N2,N1)+XAL2BS(LS,N2,N1)
     2-XAL2BS(LS2,N2,N1)
      TAL3BS(LS,N2,N1)=TAL3BS(LS,N2,N1)+XAL3BS(LS,N2,N1)
     2-XAL3BS(LS2,N2,N1)
      TAL4BS(LS,N2,N1)=TAL4BS(LS,N2,N1)+XAL4BS(LS,N2,N1)
     2-XAL4BS(LS2,N2,N1)
      TALSBS(LS,N2,N1)=TALSBS(LS,N2,N1)+XALSBS(LS,N2,N1)
     2-XALSBS(LS2,N2,N1)
      TFE1BS(LS,N2,N1)=TFE1BS(LS,N2,N1)+XFE1BS(LS,N2,N1)
     2-XFE1BS(LS2,N2,N1)
      TFE2BS(LS,N2,N1)=TFE2BS(LS,N2,N1)+XFE2BS(LS,N2,N1)
     2-XFE2BS(LS2,N2,N1)
      TFE3BS(LS,N2,N1)=TFE3BS(LS,N2,N1)+XFE3BS(LS,N2,N1)
     2-XFE3BS(LS2,N2,N1)
      TFE4BS(LS,N2,N1)=TFE4BS(LS,N2,N1)+XFE4BS(LS,N2,N1)
     2-XFE4BS(LS2,N2,N1)
      TFESBS(LS,N2,N1)=TFESBS(LS,N2,N1)+XFESBS(LS,N2,N1)
     2-XFESBS(LS2,N2,N1)
      TCAOBS(LS,N2,N1)=TCAOBS(LS,N2,N1)+XCAOBS(LS,N2,N1)
     2-XCAOBS(LS2,N2,N1)
      TCACBS(LS,N2,N1)=TCACBS(LS,N2,N1)+XCACBS(LS,N2,N1)
     2-XCACBS(LS2,N2,N1)
      TCAHBS(LS,N2,N1)=TCAHBS(LS,N2,N1)+XCAHBS(LS,N2,N1)
     2-XCAHBS(LS2,N2,N1)
      TCASBS(LS,N2,N1)=TCASBS(LS,N2,N1)+XCASBS(LS,N2,N1)
     2-XCASBS(LS2,N2,N1)
      TMGOBS(LS,N2,N1)=TMGOBS(LS,N2,N1)+XMGOBS(LS,N2,N1)
     2-XMGOBS(LS2,N2,N1)
      TMGCBS(LS,N2,N1)=TMGCBS(LS,N2,N1)+XMGCBS(LS,N2,N1)
     2-XMGCBS(LS2,N2,N1)
      TMGHBS(LS,N2,N1)=TMGHBS(LS,N2,N1)+XMGHBS(LS,N2,N1)
     2-XMGHBS(LS2,N2,N1)
      TMGSBS(LS,N2,N1)=TMGSBS(LS,N2,N1)+XMGSBS(LS,N2,N1)
     2-XMGSBS(LS2,N2,N1)
      TNACBS(LS,N2,N1)=TNACBS(LS,N2,N1)+XNACBS(LS,N2,N1)
     2-XNACBS(LS2,N2,N1)
      TNASBS(LS,N2,N1)=TNASBS(LS,N2,N1)+XNASBS(LS,N2,N1)
     2-XNASBS(LS2,N2,N1)
      TKASBS(LS,N2,N1)=TKASBS(LS,N2,N1)+XKASBS(LS,N2,N1)
     2-XKASBS(LS2,N2,N1)
      TH0PBS(LS,N2,N1)=TH0PBS(LS,N2,N1)+XH0PBS(LS,N2,N1)
     2-XH0PBS(LS2,N2,N1)
      TH3PBS(LS,N2,N1)=TH3PBS(LS,N2,N1)+XH3PBS(LS,N2,N1)
     2-XH3PBS(LS2,N2,N1)
      TF1PBS(LS,N2,N1)=TF1PBS(LS,N2,N1)+XF1PBS(LS,N2,N1)
     2-XF1PBS(LS2,N2,N1)
      TF2PBS(LS,N2,N1)=TF2PBS(LS,N2,N1)+XF2PBS(LS,N2,N1)
     2-XF2PBS(LS2,N2,N1)
      TC0PBS(LS,N2,N1)=TC0PBS(LS,N2,N1)+XC0PBS(LS,N2,N1)
     2-XC0PBS(LS2,N2,N1)
      TC1PBS(LS,N2,N1)=TC1PBS(LS,N2,N1)+XC1PBS(LS,N2,N1)
     2-XC1PBS(LS2,N2,N1)
      TC2PBS(LS,N2,N1)=TC2PBS(LS,N2,N1)+XC2PBS(LS,N2,N1)
     2-XC2PBS(LS2,N2,N1)
      TM1PBS(LS,N2,N1)=TM1PBS(LS,N2,N1)+XM1PBS(LS,N2,N1)
     2-XM1PBS(LS2,N2,N1)
      ENDIF
C
C     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
C
      ELSE
      TFLWS(LS,N2,N1)=TFLWS(LS,N2,N1)+XFLWS(LS,N2,N1)
      TFLWW(LS,N2,N1)=TFLWW(LS,N2,N1)+XFLWW(LS,N2,N1)
     2-FLSWR(LS,N2,N1)-FLSW(LS,N2,N1)-FLSWH(LS,N2,N1)
      TFLWV(LS,N2,N1)=TFLWV(LS,N2,N1)+XFLWV(LS,N2,N1)
     2-FLSVR(LS,N2,N1)-FLSV(LS,N2,N1)
      TFLWI(LS,N2,N1)=TFLWI(LS,N2,N1)+XFLWI(LS,N2,N1)
      THFLWW(LS,N2,N1)=THFLWW(LS,N2,N1)+XHFLWW(LS,N2,N1)
     2-HFLSWR(LS,N2,N1)-HFLSW(LS,N2,N1)
C     IF(LS.EQ.5)THEN
C     WRITE(*,7755)'LS',I,J,N1,N2,LS,LS2
C    2,TFLWW(LS,N2,N1),XFLWW(LS,N2,N1)
C    2,XFLWW(LS2,N2,N1),TFLWV(LS,N2,N1),XFLWV(LS,N2,N1)
C    2,XFLWV(LS2,N2,N1),FLSWR(LS,N2,N1),FLSW(LS,N2,N1),FLSWH(LS,N2,N1)
C    2,FLSVR(LS,N2,N1),FLSV(LS,N2,N1)
7755  FORMAT(A8,6I4,100E14.6)
C     ENDIF
      TCOBLS(LS,N2,N1)=TCOBLS(LS,N2,N1)+XCOBLS(LS,N2,N1)
     2-XCOFLS(3,0,N2,N1)-XCOFLS(3,NUM(N2,N1),N2,N1)
     3-XCOFHS(3,NUM(N2,N1),N2,N1)
      TCHBLS(LS,N2,N1)=TCHBLS(LS,N2,N1)+XCHBLS(LS,N2,N1)
     2-XCHFLS(3,0,N2,N1)-XCHFLS(3,NUM(N2,N1),N2,N1)
     3-XCHFHS(3,NUM(N2,N1),N2,N1)
      TOXBLS(LS,N2,N1)=TOXBLS(LS,N2,N1)+XOXBLS(LS,N2,N1)
     2-XOXFLS(3,0,N2,N1)-XOXFLS(3,NUM(N2,N1),N2,N1)
     3-XOXFHS(3,NUM(N2,N1),N2,N1)
      TNGBLS(LS,N2,N1)=TNGBLS(LS,N2,N1)+XNGBLS(LS,N2,N1)
     2-XNGFLS(3,0,N2,N1)-XNGFLS(3,NUM(N2,N1),N2,N1)
     3-XNGFHS(3,NUM(N2,N1),N2,N1)
      TN2BLS(LS,N2,N1)=TN2BLS(LS,N2,N1)+XN2BLS(LS,N2,N1)
     2-XN2FLS(3,0,N2,N1)-XN2FLS(3,NUM(N2,N1),N2,N1)
     3-XN2FHS(3,NUM(N2,N1),N2,N1)
      TN4BLW(LS,N2,N1)=TN4BLW(LS,N2,N1)+XN4BLW(LS,N2,N1)
     2-XN4FLW(3,0,N2,N1)-XN4FLW(3,NUM(N2,N1),N2,N1)
     3-XN4FHW(3,NUM(N2,N1),N2,N1)-XN4FLB(3,NUM(N2,N1),N2,N1)
     3-XN4FHB(3,NUM(N2,N1),N2,N1)
      TN3BLW(LS,N2,N1)=TN3BLW(LS,N2,N1)+XN3BLW(LS,N2,N1)
     2-XN3FLW(3,0,N2,N1)-XN3FLW(3,NUM(N2,N1),N2,N1)
     3-XN3FHW(3,NUM(N2,N1),N2,N1)-XN3FLB(3,NUM(N2,N1),N2,N1)
     3-XN3FHB(3,NUM(N2,N1),N2,N1)
      TNOBLW(LS,N2,N1)=TNOBLW(LS,N2,N1)+XNOBLW(LS,N2,N1)
     2-XNOFLW(3,0,N2,N1)-XNOFLW(3,NUM(N2,N1),N2,N1)
     3-XNOFHW(3,NUM(N2,N1),N2,N1)-XNOFLB(3,NUM(N2,N1),N2,N1)
     3-XNOFHB(3,NUM(N2,N1),N2,N1)
      TH1PBS(LS,N2,N1)=TH1PBS(LS,N2,N1)+XH1PBS(LS,N2,N1)
     2-XH1PFS(3,0,N2,N1)-XH1PFS(3,NUM(N2,N1),N2,N1)
     3-XH1PHS(3,NUM(N2,N1),N2,N1)-XH1BFB(3,NUM(N2,N1),N2,N1)
     3-XH1BHB(3,NUM(N2,N1),N2,N1)
      TH2PBS(LS,N2,N1)=TH2PBS(LS,N2,N1)+XH2PBS(LS,N2,N1)
     2-XH2PFS(3,0,N2,N1)-XH2PFS(3,NUM(N2,N1),N2,N1)
     3-XH2PHS(3,NUM(N2,N1),N2,N1)-XH2BFB(3,NUM(N2,N1),N2,N1)
     3-XH2BHB(3,NUM(N2,N1),N2,N1)
      IF(ISALTG.NE.0)THEN
      TALBLS(LS,NY,NX)=TALBLS(LS,NY,NX)+XALBLS(LS,NY,NX)
     2-XALFLS(3,0,N2,N1)-XALFLS(3,NUM(N2,N1),N2,N1)
     3-XALFHS(3,NUM(N2,N1),N2,N1)
      TFEBLS(LS,NY,NX)=TFEBLS(LS,NY,NX)+XFEBLS(LS,NY,NX)
     2-XFEFLS(3,0,N2,N1)-XFEFLS(3,NUM(N2,N1),N2,N1)
     3-XFEFHS(3,NUM(N2,N1),N2,N1)
      THYBLS(LS,NY,NX)=THYBLS(LS,NY,NX)+XHYBLS(LS,NY,NX)
     2-XHYFLS(3,0,N2,N1)-XHYFLS(3,NUM(N2,N1),N2,N1)
     3-XHYFHS(3,NUM(N2,N1),N2,N1)
      TCABLS(LS,NY,NX)=TCABLS(LS,NY,NX)+XCABLS(LS,NY,NX)
     2-XCAFLS(3,0,N2,N1)-XCAFLS(3,NUM(N2,N1),N2,N1)
     3-XCAFHS(3,NUM(N2,N1),N2,N1)
      TMGBLS(LS,NY,NX)=TMGBLS(LS,NY,NX)+XMGBLS(LS,NY,NX)
     2-XMGFLS(3,0,N2,N1)-XMGFLS(3,NUM(N2,N1),N2,N1)
     3-XMGFHS(3,NUM(N2,N1),N2,N1)
      TNABLS(LS,NY,NX)=TNABLS(LS,NY,NX)+XNABLS(LS,NY,NX)
     2-XNAFLS(3,0,N2,N1)-XNAFLS(3,NUM(N2,N1),N2,N1)
     3-XNAFHS(3,NUM(N2,N1),N2,N1)
      TKABLS(LS,NY,NX)=TKABLS(LS,NY,NX)+XKABLS(LS,NY,NX)
     2-XKAFLS(3,0,N2,N1)-XKAFLS(3,NUM(N2,N1),N2,N1)
     3-XKAFHS(3,NUM(N2,N1),N2,N1)
      TOHBLS(LS,NY,NX)=TOHBLS(LS,NY,NX)+XOHBLS(LS,NY,NX)
     2-XOHFLS(3,0,N2,N1)-XOHFLS(3,NUM(N2,N1),N2,N1)
     3-XOHFHS(3,NUM(N2,N1),N2,N1)
      TSOBLS(LS,NY,NX)=TSOBLS(LS,NY,NX)+XSOBLS(LS,NY,NX)
     2-XSOFLS(3,0,N2,N1)-XSOFLS(3,NUM(N2,N1),N2,N1)
     3-XSOFHS(3,NUM(N2,N1),N2,N1)
      TCLBLS(LS,NY,NX)=TCLBLS(LS,NY,NX)+XCLBLS(LS,NY,NX)
     2-XCLFLS(3,0,N2,N1)-XCLFLS(3,NUM(N2,N1),N2,N1)
     3-XCLFHS(3,NUM(N2,N1),N2,N1)
      TC3BLS(LS,NY,NX)=TC3BLS(LS,NY,NX)+XC3BLS(LS,NY,NX)
     2-XC3FLS(3,0,N2,N1)-XC3FLS(3,NUM(N2,N1),N2,N1)
     3-XC3FHS(3,NUM(N2,N1),N2,N1)
      THCBLS(LS,NY,NX)=THCBLS(LS,NY,NX)+XHCBLS(LS,NY,NX)
     2-XHCFLS(3,0,N2,N1)-XHCFLS(3,NUM(N2,N1),N2,N1)
     3-XHCFHS(3,NUM(N2,N1),N2,N1)
      TAL1BS(LS,NY,NX)=TAL1BS(LS,NY,NX)+XAL1BS(LS,NY,NX)
     2-XAL1FS(3,0,N2,N1)-XAL1FS(3,NUM(N2,N1),N2,N1)
     3-XAL1HS(3,NUM(N2,N1),N2,N1)
      TAL2BS(LS,NY,NX)=TAL2BS(LS,NY,NX)+XAL2BS(LS,NY,NX)
     2-XAL2FS(3,0,N2,N1)-XAL2FS(3,NUM(N2,N1),N2,N1)
     3-XAL2HS(3,NUM(N2,N1),N2,N1)
      TAL3BS(LS,NY,NX)=TAL3BS(LS,NY,NX)+XAL3BS(LS,NY,NX)
     2-XAL3FS(3,0,N2,N1)-XAL3FS(3,NUM(N2,N1),N2,N1)
     3-XAL3HS(3,NUM(N2,N1),N2,N1)
      TAL4BS(LS,NY,NX)=TAL4BS(LS,NY,NX)+XAL4BS(LS,NY,NX)
     2-XAL4FS(3,0,N2,N1)-XAL4FS(3,NUM(N2,N1),N2,N1)
     3-XAL4HS(3,NUM(N2,N1),N2,N1)
      TALSBS(LS,NY,NX)=TALSBS(LS,NY,NX)+XALSBS(LS,NY,NX)
     2-XALSFS(3,0,N2,N1)-XALSFS(3,NUM(N2,N1),N2,N1)
     3-XALSHS(3,NUM(N2,N1),N2,N1)
      TFE1BS(LS,NY,NX)=TFE1BS(LS,NY,NX)+XFE1BS(LS,NY,NX)
     2-XFE1FS(3,0,N2,N1)-XFE1FS(3,NUM(N2,N1),N2,N1)
     3-XFE1HS(3,NUM(N2,N1),N2,N1)
      TFE2BS(LS,NY,NX)=TFE2BS(LS,NY,NX)+XFE2BS(LS,NY,NX)
     2-XFE2FS(3,0,N2,N1)-XFE2FS(3,NUM(N2,N1),N2,N1)
     3-XFE2HS(3,NUM(N2,N1),N2,N1)
      TFE3BS(LS,NY,NX)=TFE3BS(LS,NY,NX)+XFE3BS(LS,NY,NX)
     2-XFE3FS(3,0,N2,N1)-XFE3FS(3,NUM(N2,N1),N2,N1)
     3-XFE3HS(3,NUM(N2,N1),N2,N1)
      TFE4BS(LS,NY,NX)=TFE4BS(LS,NY,NX)+XFE4BS(LS,NY,NX)
     2-XFE4FS(3,0,N2,N1)-XFE4FS(3,NUM(N2,N1),N2,N1)
     3-XFE4HS(3,NUM(N2,N1),N2,N1)
      TFESBS(LS,NY,NX)=TFESBS(LS,NY,NX)+XFESBS(LS,NY,NX)
     2-XFESFS(3,0,N2,N1)-XFESFS(3,NUM(N2,N1),N2,N1)
     3-XFESHS(3,NUM(N2,N1),N2,N1)
      TCAOBS(LS,NY,NX)=TCAOBS(LS,NY,NX)+XCAOBS(LS,NY,NX)
     2-XCAOFS(3,0,N2,N1)-XCAOFS(3,NUM(N2,N1),N2,N1)
     3-XCAOHS(3,NUM(N2,N1),N2,N1)
      TCACBS(LS,NY,NX)=TCACBS(LS,NY,NX)+XCACBS(LS,NY,NX)
     2-XCACFS(3,0,N2,N1)-XCACFS(3,NUM(N2,N1),N2,N1)
     3-XCACHS(3,NUM(N2,N1),N2,N1)
      TCAHBS(LS,NY,NX)=TCAHBS(LS,NY,NX)+XCAHBS(LS,NY,NX)
     2-XCAHFS(3,0,N2,N1)-XCAHFS(3,NUM(N2,N1),N2,N1)
     3-XALFHS(3,NUM(N2,N1),N2,N1)
      TCASBS(LS,NY,NX)=TCASBS(LS,NY,NX)+XCASBS(LS,NY,NX)
     2-XCASFS(3,0,N2,N1)-XCASFS(3,NUM(N2,N1),N2,N1)
     3-XCASHS(3,NUM(N2,N1),N2,N1)
      TMGOBS(LS,NY,NX)=TMGOBS(LS,NY,NX)+XMGOBS(LS,NY,NX)
     2-XMGOFS(3,0,N2,N1)-XMGOFS(3,NUM(N2,N1),N2,N1)
     3-XMGOHS(3,NUM(N2,N1),N2,N1)
      TMGCBS(LS,NY,NX)=TMGCBS(LS,NY,NX)+XMGCBS(LS,NY,NX)
     2-XMGCFS(3,0,N2,N1)-XMGCFS(3,NUM(N2,N1),N2,N1)
     3-XMGCHS(3,NUM(N2,N1),N2,N1)
      TMGHBS(LS,NY,NX)=TMGHBS(LS,NY,NX)+XMGHBS(LS,NY,NX)
     2-XMGHFS(3,0,N2,N1)-XMGHFS(3,NUM(N2,N1),N2,N1)
     3-XMGHHS(3,NUM(N2,N1),N2,N1)
      TMGSBS(LS,NY,NX)=TMGSBS(LS,NY,NX)+XMGSBS(LS,NY,NX)
     2-XMGSFS(3,0,N2,N1)-XMGSFS(3,NUM(N2,N1),N2,N1)
     3-XMGSHS(3,NUM(N2,N1),N2,N1)
      TNACBS(LS,NY,NX)=TNACBS(LS,NY,NX)+XNACBS(LS,NY,NX)
     2-XNACFS(3,0,N2,N1)-XNACFS(3,NUM(N2,N1),N2,N1)
     3-XNACHS(3,NUM(N2,N1),N2,N1)
      TNASBS(LS,NY,NX)=TNASBS(LS,NY,NX)+XNASBS(LS,NY,NX)
     2-XNASFS(3,0,N2,N1)-XNASFS(3,NUM(N2,N1),N2,N1)
     3-XNASHS(3,NUM(N2,N1),N2,N1)
      TKASBS(LS,NY,NX)=TKASBS(LS,NY,NX)+XKASBS(LS,NY,NX)
     2-XKASFS(3,0,N2,N1)-XKASFS(3,NUM(N2,N1),N2,N1)
     3-XKASHS(3,NUM(N2,N1),N2,N1)
      TH0PBS(LS,NY,NX)=TH0PBS(LS,NY,NX)+XH0PBS(LS,NY,NX)
     2-XH0PFS(3,0,N2,N1)-XH0PFS(3,NUM(N2,N1),N2,N1)
     3-XH0PHS(3,NUM(N2,N1),N2,N1)-XH0BFB(3,NUM(N2,N1),N2,N1)
     3-XH0BHB(3,NUM(N2,N1),N2,N1)
      TH3PBS(LS,NY,NX)=TH3PBS(LS,NY,NX)+XH3PBS(LS,NY,NX)
     2-XH3PFS(3,0,N2,N1)-XH3PHS(3,NUM(N2,N1),N2,N1)
     3-XH3PHS(3,NUM(N2,N1),N2,N1)-XH3BFB(3,NUM(N2,N1),N2,N1)
     3-XH3BHB(3,NUM(N2,N1),N2,N1)
      TF1PBS(LS,NY,NX)=TF1PBS(LS,NY,NX)+XF1PBS(LS,NY,NX)
     2-XF1PFS(3,0,N2,N1)-XF1PFS(3,NUM(N2,N1),N2,N1)
     3-XF1PHS(3,NUM(N2,N1),N2,N1)-XF1BFB(3,NUM(N2,N1),N2,N1)
     3-XF1BHB(3,NUM(N2,N1),N2,N1)
      TF2PBS(LS,NY,NX)=TF2PBS(LS,NY,NX)+XF2PBS(LS,NY,NX)
     2-XF2PFS(3,0,N2,N1)-XF2PFS(3,NUM(N2,N1),N2,N1)
     3-XF2PHS(3,NUM(N2,N1),N2,N1)-XF2BFB(3,NUM(N2,N1),N2,N1)
     3-XF2BHB(3,NUM(N2,N1),N2,N1)
      TC0PBS(LS,NY,NX)=TC0PBS(LS,NY,NX)+XC0PBS(LS,NY,NX)
     2-XC0PFS(3,0,N2,N1)-XC0PFS(3,NUM(N2,N1),N2,N1)
     3-XC0PHS(3,NUM(N2,N1),N2,N1)-XC0BFB(3,NUM(N2,N1),N2,N1)
     3-XC0BHB(3,NUM(N2,N1),N2,N1)
      TC1PBS(LS,NY,NX)=TC1PBS(LS,NY,NX)+XC1PBS(LS,NY,NX)
     2-XC1PFS(3,0,N2,N1)-XC1PFS(3,NUM(N2,N1),N2,N1)
     3-XC1PHS(3,NUM(N2,N1),N2,N1)-XC1BFB(3,NUM(N2,N1),N2,N1)
     3-XC1BHB(3,NUM(N2,N1),N2,N1)
      TC2PBS(LS,NY,NX)=TC2PBS(LS,NY,NX)+XC2PBS(LS,NY,NX)
     2-XC2PFS(3,0,N2,N1)-XC2PFS(3,NUM(N2,N1),N2,N1)
     3-XC2PHS(3,NUM(N2,N1),N2,N1)-XC2BFB(3,NUM(N2,N1),N2,N1)
     3-XC2BHB(3,NUM(N2,N1),N2,N1)
      TM1PBS(LS,NY,NX)=TM1PBS(LS,NY,NX)+XM1PBS(LS,NY,NX)
     2-XM1PFS(3,0,N2,N1)-XM1PFS(3,NUM(N2,N1),N2,N1)
     3-XM1PHS(3,NUM(N2,N1),N2,N1)-XM1BFB(3,NUM(N2,N1),N2,N1)
     3-XM1BHB(3,NUM(N2,N1),N2,N1)
      ENDIF
      ENDIF
C
C     WATER,GAS,SOLUTE FLUXES INTO SNOWPACK SURFACE
C
C     TFLWS,TFLWW,TFLWV,TFLWI=net fluxes of snow,water,vapor,ice at
C        snowpack surface
C     THFLWW=convective heat fluxes of snow,water,ice at snowpack
C        surface
C     XFLWS,XFLWW,XFLWV,XFLWI=snow,water,vapor,ice transfer
C        from watsub.f
C     T*BLS=net solute flux in snowpack
C     X*BLS=solute flux in snowpack from trnsfr.f
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,H1P=HPO4,H2P=H2PO4
C
      ELSEIF(LS.EQ.1)THEN
      IF(XFLWS(LS,N2,N1).NE.0.0)THEN
      TFLWS(LS,N2,N1)=TFLWS(LS,N2,N1)+XFLWS(LS,N2,N1)
      TFLWW(LS,N2,N1)=TFLWW(LS,N2,N1)+XFLWW(LS,N2,N1)
      TFLWV(LS,N2,N1)=TFLWV(LS,N2,N1)+XFLWV(LS,N2,N1)
      TFLWI(LS,N2,N1)=TFLWI(LS,N2,N1)+XFLWI(LS,N2,N1)
      THFLWW(LS,N2,N1)=THFLWW(LS,N2,N1)+XHFLWW(LS,N2,N1)
      TCOBLS(LS,N2,N1)=TCOBLS(LS,N2,N1)+XCOBLS(LS,N2,N1)
      TCHBLS(LS,N2,N1)=TCHBLS(LS,N2,N1)+XCHBLS(LS,N2,N1)
      TOXBLS(LS,N2,N1)=TOXBLS(LS,N2,N1)+XOXBLS(LS,N2,N1)
      TNGBLS(LS,N2,N1)=TNGBLS(LS,N2,N1)+XNGBLS(LS,N2,N1)
      TN2BLS(LS,N2,N1)=TN2BLS(LS,N2,N1)+XN2BLS(LS,N2,N1)
      TN4BLW(LS,N2,N1)=TN4BLW(LS,N2,N1)+XN4BLW(LS,N2,N1)
      TN3BLW(LS,N2,N1)=TN3BLW(LS,N2,N1)+XN3BLW(LS,N2,N1)
      TNOBLW(LS,N2,N1)=TNOBLW(LS,N2,N1)+XNOBLW(LS,N2,N1)
      TH1PBS(LS,N2,N1)=TH1PBS(LS,N2,N1)+XH1PBS(LS,N2,N1)
      TH2PBS(LS,N2,N1)=TH2PBS(LS,N2,N1)+XH2PBS(LS,N2,N1)
C
C     SALT FLUXES INTO SNOWPACK SURFACE
C
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C
      IF(ISALTG.NE.0)THEN
      TALBLS(LS,N2,N1)=TALBLS(LS,N2,N1)+XALBLS(LS,N2,N1)
      TFEBLS(LS,N2,N1)=TFEBLS(LS,N2,N1)+XFEBLS(LS,N2,N1)
      THYBLS(LS,N2,N1)=THYBLS(LS,N2,N1)+XHYBLS(LS,N2,N1)
      TCABLS(LS,N2,N1)=TCABLS(LS,N2,N1)+XCABLS(LS,N2,N1)
      TMGBLS(LS,N2,N1)=TMGBLS(LS,N2,N1)+XMGBLS(LS,N2,N1)
      TNABLS(LS,N2,N1)=TNABLS(LS,N2,N1)+XNABLS(LS,N2,N1)
      TKABLS(LS,N2,N1)=TKABLS(LS,N2,N1)+XKABLS(LS,N2,N1)
      TOHBLS(LS,N2,N1)=TOHBLS(LS,N2,N1)+XOHBLS(LS,N2,N1)
      TSOBLS(LS,N2,N1)=TSOBLS(LS,N2,N1)+XSOBLS(LS,N2,N1)
      TCLBLS(LS,N2,N1)=TCLBLS(LS,N2,N1)+XCLBLS(LS,N2,N1)
      TC3BLS(LS,N2,N1)=TC3BLS(LS,N2,N1)+XC3BLS(LS,N2,N1)
      THCBLS(LS,N2,N1)=THCBLS(LS,N2,N1)+XHCBLS(LS,N2,N1)
      TAL1BS(LS,N2,N1)=TAL1BS(LS,N2,N1)+XAL1BS(LS,N2,N1)
      TAL2BS(LS,N2,N1)=TAL2BS(LS,N2,N1)+XAL2BS(LS,N2,N1)
      TAL3BS(LS,N2,N1)=TAL3BS(LS,N2,N1)+XAL3BS(LS,N2,N1)
      TAL4BS(LS,N2,N1)=TAL4BS(LS,N2,N1)+XAL4BS(LS,N2,N1)
      TALSBS(LS,N2,N1)=TALSBS(LS,N2,N1)+XALSBS(LS,N2,N1)
      TFE1BS(LS,N2,N1)=TFE1BS(LS,N2,N1)+XFE1BS(LS,N2,N1)
      TFE2BS(LS,N2,N1)=TFE2BS(LS,N2,N1)+XFE2BS(LS,N2,N1)
      TFE3BS(LS,N2,N1)=TFE3BS(LS,N2,N1)+XFE3BS(LS,N2,N1)
      TFE4BS(LS,N2,N1)=TFE4BS(LS,N2,N1)+XFE4BS(LS,N2,N1)
      TFESBS(LS,N2,N1)=TFESBS(LS,N2,N1)+XFESBS(LS,N2,N1)
      TCAOBS(LS,N2,N1)=TCAOBS(LS,N2,N1)+XCAOBS(LS,N2,N1)
      TCACBS(LS,N2,N1)=TCACBS(LS,N2,N1)+XCACBS(LS,N2,N1)
      TCAHBS(LS,N2,N1)=TCAHBS(LS,N2,N1)+XCAHBS(LS,N2,N1)
      TCASBS(LS,N2,N1)=TCASBS(LS,N2,N1)+XCASBS(LS,N2,N1)
      TMGOBS(LS,N2,N1)=TMGOBS(LS,N2,N1)+XMGOBS(LS,N2,N1)
      TMGCBS(LS,N2,N1)=TMGCBS(LS,N2,N1)+XMGCBS(LS,N2,N1)
      TMGHBS(LS,N2,N1)=TMGHBS(LS,N2,N1)+XMGHBS(LS,N2,N1)
      TMGSBS(LS,N2,N1)=TMGSBS(LS,N2,N1)+XMGSBS(LS,N2,N1)
      TNACBS(LS,N2,N1)=TNACBS(LS,N2,N1)+XNACBS(LS,N2,N1)
      TNASBS(LS,N2,N1)=TNASBS(LS,N2,N1)+XNASBS(LS,N2,N1)
      TKASBS(LS,N2,N1)=TKASBS(LS,N2,N1)+XKASBS(LS,N2,N1)
      TH0PBS(LS,N2,N1)=TH0PBS(LS,N2,N1)+XH0PBS(LS,N2,N1)
      TH3PBS(LS,N2,N1)=TH3PBS(LS,N2,N1)+XH3PBS(LS,N2,N1)
      TF1PBS(LS,N2,N1)=TF1PBS(LS,N2,N1)+XF1PBS(LS,N2,N1)
      TF2PBS(LS,N2,N1)=TF2PBS(LS,N2,N1)+XF2PBS(LS,N2,N1)
      TC0PBS(LS,N2,N1)=TC0PBS(LS,N2,N1)+XC0PBS(LS,N2,N1)
      TC1PBS(LS,N2,N1)=TC1PBS(LS,N2,N1)+XC1PBS(LS,N2,N1)
      TC2PBS(LS,N2,N1)=TC2PBS(LS,N2,N1)+XC2PBS(LS,N2,N1)
      TM1PBS(LS,N2,N1)=TM1PBS(LS,N2,N1)+XM1PBS(LS,N2,N1)
      ENDIF
      ENDIF
      ENDIF
1205  CONTINUE
      ENDIF
C
C     INCOMING FLUXES FROM EROSION
C
C     NN=boundary:N=1:NN=1 east,NN=2 west, N=2:NN=1 south,NN=2 north
C
      IF(N.NE.3.AND.(IERSNG.EQ.1.OR.IERSNG.EQ.3))THEN
      DO 9350 NN=1,2
C
C     T*ER=net sediment flux
C     X*ER=sediment flux from erosion.f
C
      IF(ABS(XSEDER(N,NN,N2,N1)).GT.ZEROS(N2,N1)
     3.OR.ABS(XSEDER(N,NN,N5,N4)).GT.ZEROS(N5,N4))THEN
C
C     SEDIMENT FROM EROSION
C
C     sediment code:XSED=total,XSAN=sand,XSIL=silt,XCLA=clay
C                  :CEC=cation exchange capacity
C                  :AEC=anion exchange capacity
C
      TSEDER(N2,N1)=TSEDER(N2,N1)+XSEDER(N,NN,N2,N1)
      TSANER(N2,N1)=TSANER(N2,N1)+XSANER(N,NN,N2,N1)
      TSILER(N2,N1)=TSILER(N2,N1)+XSILER(N,NN,N2,N1)
      TCLAER(N2,N1)=TCLAER(N2,N1)+XCLAER(N,NN,N2,N1)
      TCECER(N2,N1)=TCECER(N2,N1)+XCECER(N,NN,N2,N1)
      TAECER(N2,N1)=TAECER(N2,N1)+XAECER(N,NN,N2,N1)
C     IF(N1.EQ.1)THEN
C     WRITE(*,1139)'TSEDER1',I,J,NFZ,N,NN,N1,N2
C    2,TSEDER(N2,N1),XSEDER(N,NN,N2,N1)
1139  FORMAT(A8,7I4,12E12.4)
C     ENDIF
C
C     FERTILIZER FROM EROSION
C
C     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3 in non-band
C                  :NH4B,NH3B,NHUB,NO3B=NH4,NH3,urea,NO3 in band
C
      TNH4ER(N2,N1)=TNH4ER(N2,N1)+XNH4ER(N,NN,N2,N1)
      TNH3ER(N2,N1)=TNH3ER(N2,N1)+XNH3ER(N,NN,N2,N1)
      TNHUER(N2,N1)=TNHUER(N2,N1)+XNHUER(N,NN,N2,N1)
      TNO3ER(N2,N1)=TNO3ER(N2,N1)+XNO3ER(N,NN,N2,N1)
      TNH4EB(N2,N1)=TNH4EB(N2,N1)+XNH4EB(N,NN,N2,N1)
      TNH3EB(N2,N1)=TNH3EB(N2,N1)+XNH3EB(N,NN,N2,N1)
      TNHUEB(N2,N1)=TNHUEB(N2,N1)+XNHUEB(N,NN,N2,N1)
      TNO3EB(N2,N1)=TNO3EB(N2,N1)+XNO3EB(N,NN,N2,N1)
C
C     EXCHANGEABLE CATIONS AND ANIONS FROM EROSION
C
C     sediment code
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
C           =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
C       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C
      TN4ER(N2,N1)=TN4ER(N2,N1)+XN4ER(N,NN,N2,N1)
      TNBER(N2,N1)=TNBER(N2,N1)+XNBER(N,NN,N2,N1)
      THYER(N2,N1)=THYER(N2,N1)+XHYER(N,NN,N2,N1)
      TALER(N2,N1)=TALER(N2,N1)+XALER(N,NN,N2,N1)
      TFEER(N2,N1)=TFEER(N2,N1)+XFEER(N,NN,N2,N1)
      TCAER(N2,N1)=TCAER(N2,N1)+XCAER(N,NN,N2,N1)
      TMGER(N2,N1)=TMGER(N2,N1)+XMGER(N,NN,N2,N1)
      TNAER(N2,N1)=TNAER(N2,N1)+XNAER(N,NN,N2,N1)
      TKAER(N2,N1)=TKAER(N2,N1)+XKAER(N,NN,N2,N1)
      THCER(N2,N1)=THCER(N2,N1)+XHCER(N,NN,N2,N1)
      TAL2ER(N2,N1)=TAL2ER(N2,N1)+XAL2ER(N,NN,N2,N1)
      TFE2ER(N2,N1)=TFE2ER(N2,N1)+XFE2ER(N,NN,N2,N1)
      TOH0ER(N2,N1)=TOH0ER(N2,N1)+XOH0ER(N,NN,N2,N1)
      TOH1ER(N2,N1)=TOH1ER(N2,N1)+XOH1ER(N,NN,N2,N1)
      TOH2ER(N2,N1)=TOH2ER(N2,N1)+XOH2ER(N,NN,N2,N1)
      TH1PER(N2,N1)=TH1PER(N2,N1)+XH1PER(N,NN,N2,N1)
      TH2PER(N2,N1)=TH2PER(N2,N1)+XH2PER(N,NN,N2,N1)
      TOH0EB(N2,N1)=TOH0EB(N2,N1)+XOH0EB(N,NN,N2,N1)
      TOH1EB(N2,N1)=TOH1EB(N2,N1)+XOH1EB(N,NN,N2,N1)
      TOH2EB(N2,N1)=TOH2EB(N2,N1)+XOH2EB(N,NN,N2,N1)
      TH1PEB(N2,N1)=TH1PEB(N2,N1)+XH1PEB(N,NN,N2,N1)
      TH2PEB(N2,N1)=TH2PEB(N2,N1)+XH2PEB(N,NN,N2,N1)
C
C     PRECIPITATES FROM EROSION
C
C     sediment code
C       :PALO,PFEO=precip AlOH,FeOH
C       :PCAC,PCAS=precip CaCO3,CaSO4
C       :PALP,PFEP=precip AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precip AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite in non-band
C       :PCPMB,PCPDB,PCPHB=precip CaH2PO4,CaHPO4,apatite in band
C
      TALOER(N2,N1)=TALOER(N2,N1)+PALOER(N,NN,N2,N1)
      TFEOER(N2,N1)=TFEOER(N2,N1)+PFEOER(N,NN,N2,N1)
      TCACER(N2,N1)=TCACER(N2,N1)+PCACER(N,NN,N2,N1)
      TCASER(N2,N1)=TCASER(N2,N1)+PCASER(N,NN,N2,N1)
      TALPER(N2,N1)=TALPER(N2,N1)+PALPER(N,NN,N2,N1)
      TFEPER(N2,N1)=TFEPER(N2,N1)+PFEPER(N,NN,N2,N1)
      TCPDER(N2,N1)=TCPDER(N2,N1)+PCPDER(N,NN,N2,N1)
      TCPHER(N2,N1)=TCPHER(N2,N1)+PCPHER(N,NN,N2,N1)
      TCPMER(N2,N1)=TCPMER(N2,N1)+PCPMER(N,NN,N2,N1)
      TALPEB(N2,N1)=TALPEB(N2,N1)+PALPEB(N,NN,N2,N1)
      TFEPEB(N2,N1)=TFEPEB(N2,N1)+PFEPEB(N,NN,N2,N1)
      TCPDEB(N2,N1)=TCPDEB(N2,N1)+PCPDEB(N,NN,N2,N1)
      TCPHEB(N2,N1)=TCPHEB(N2,N1)+PCPHEB(N,NN,N2,N1)
      TCPMEB(N2,N1)=TCPMEB(N2,N1)+PCPMEB(N,NN,N2,N1)
C
C     ORGANIC MATTER FROM EROSION
C
C     sediment code
C        :OMC,OMN,OMP=microbial C,N,P
C        :ORC,ORN,ORP=microbial residue C,N,P
C        :OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
C        :OSC,OSA,OSN,OSP=SOC,colonized SOC,SON,SOP
C           (K=0:woody litter, K=1:non-woody litter,
C            K=2:manure, K=3:POC, K=4:humus)
C
      DO 9380 K=0,5
      DO 9380 NO=1,7
      DO 9380 M=1,3
      TOMCER(M,NO,K,N2,N1)=TOMCER(M,NO,K,N2,N1)
     2+OMCER(M,NO,K,N,NN,N2,N1)
      TOMNER(M,NO,K,N2,N1)=TOMNER(M,NO,K,N2,N1)
     2+OMNER(M,NO,K,N,NN,N2,N1)
      TOMPER(M,NO,K,N2,N1)=TOMPER(M,NO,K,N2,N1)
     2+OMPER(M,NO,K,N,NN,N2,N1)
9380  CONTINUE
      DO 9375 K=0,4
      DO 9370 M=1,2
      TORCER(M,K,N2,N1)=TORCER(M,K,N2,N1)+ORCER(M,K,N,NN,N2,N1)
      TORNER(M,K,N2,N1)=TORNER(M,K,N2,N1)+ORNER(M,K,N,NN,N2,N1)
      TORPER(M,K,N2,N1)=TORPER(M,K,N2,N1)+ORPER(M,K,N,NN,N2,N1)
9370  CONTINUE
      TOHCER(K,N2,N1)=TOHCER(K,N2,N1)+OHCER(K,N,NN,N2,N1)
      TOHNER(K,N2,N1)=TOHNER(K,N2,N1)+OHNER(K,N,NN,N2,N1)
      TOHPER(K,N2,N1)=TOHPER(K,N2,N1)+OHPER(K,N,NN,N2,N1)
      TOHAER(K,N2,N1)=TOHAER(K,N2,N1)+OHAER(K,N,NN,N2,N1)
      DO 9365 M=1,5
      TOSCER(M,K,N2,N1)=TOSCER(M,K,N2,N1)+OSCER(M,K,N,NN,N2,N1)
      TOSAER(M,K,N2,N1)=TOSAER(M,K,N2,N1)+OSAER(M,K,N,NN,N2,N1)
      TOSNER(M,K,N2,N1)=TOSNER(M,K,N2,N1)+OSNER(M,K,N,NN,N2,N1)
      TOSPER(M,K,N2,N1)=TOSPER(M,K,N2,N1)+OSPER(M,K,N,NN,N2,N1)
9365  CONTINUE
9375  CONTINUE
C     IF(NN.EQ.2)THEN
C
C     OUTGOING FLUXES FROM EROSION E AND S
C
C     NN=boundary:N=1:NN=1 east,NN=2 west, N=2:NN=1 south,NN=2 north
C
C
C     SEDIMENT FROM EROSION
C
C     sediment code:XSED=total,XSAN=sand,XSIL=silt,XCLA=clay
C                  :CEC=cation exchange capacity
C                  :AEC=anion exchange capacity
C
      TSEDER(N2,N1)=TSEDER(N2,N1)-XSEDER(N,NN,N5,N4)
      TSANER(N2,N1)=TSANER(N2,N1)-XSANER(N,NN,N5,N4)
      TSILER(N2,N1)=TSILER(N2,N1)-XSILER(N,NN,N5,N4)
      TCLAER(N2,N1)=TCLAER(N2,N1)-XCLAER(N,NN,N5,N4)
      TCECER(N2,N1)=TCECER(N2,N1)-XCECER(N,NN,N5,N4)
      TAECER(N2,N1)=TAECER(N2,N1)-XAECER(N,NN,N5,N4)
C     IF(N1.EQ.1)THEN
C     WRITE(*,1138)'TSEDER2',I,J,NFZ,N,NN,N1,N2,N4,N5
C    2,TSEDER(N2,N1),XSEDER(N,NN,N4,N5)
1138  FORMAT(A8,9I4,12E12.4)
C     ENDIF
C
C     FERTILIZER FROM EROSION
C
C     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3 in non-band
C                  :NH4B,NH3B,NHUB,NO3B=NH4,NH3,urea,NO3 in band
C
      TNH4ER(N2,N1)=TNH4ER(N2,N1)-XNH4ER(N,NN,N5,N4)
      TNH3ER(N2,N1)=TNH3ER(N2,N1)-XNH3ER(N,NN,N5,N4)
      TNHUER(N2,N1)=TNHUER(N2,N1)-XNHUER(N,NN,N5,N4)
      TNO3ER(N2,N1)=TNO3ER(N2,N1)-XNO3ER(N,NN,N5,N4)
      TNH4EB(N2,N1)=TNH4EB(N2,N1)-XNH4EB(N,NN,N5,N4)
      TNH3EB(N2,N1)=TNH3EB(N2,N1)-XNH3EB(N,NN,N5,N4)
      TNHUEB(N2,N1)=TNHUEB(N2,N1)-XNHUEB(N,NN,N5,N4)
      TNO3EB(N2,N1)=TNO3EB(N2,N1)-XNO3EB(N,NN,N5,N4)
C
C     EXCHANGEABLE CATIONS AND ANIONS FROM EROSION
C
C     sediment code
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
C           =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
C       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C
      TN4ER(N2,N1)=TN4ER(N2,N1)-XN4ER(N,NN,N5,N4)
      TNBER(N2,N1)=TNBER(N2,N1)-XNBER(N,NN,N5,N4)
      THYER(N2,N1)=THYER(N2,N1)-XHYER(N,NN,N5,N4)
      TALER(N2,N1)=TALER(N2,N1)-XALER(N,NN,N5,N4)
      TFEER(N2,N1)=TFEER(N2,N1)-XFEER(N,NN,N5,N4)
      TCAER(N2,N1)=TCAER(N2,N1)-XCAER(N,NN,N5,N4)
      TMGER(N2,N1)=TMGER(N2,N1)-XMGER(N,NN,N5,N4)
      TNAER(N2,N1)=TNAER(N2,N1)-XNAER(N,NN,N5,N4)
      TKAER(N2,N1)=TKAER(N2,N1)-XKAER(N,NN,N5,N4)
      THCER(N2,N1)=THCER(N2,N1)-XHCER(N,NN,N5,N4)
      TAL2ER(N2,N1)=TAL2ER(N2,N1)-XAL2ER(N,NN,N5,N4)
      TFE2ER(N2,N1)=TFE2ER(N2,N1)-XFE2ER(N,NN,N5,N4)
      TOH0ER(N2,N1)=TOH0ER(N2,N1)-XOH0ER(N,NN,N5,N4)
      TOH1ER(N2,N1)=TOH1ER(N2,N1)-XOH1ER(N,NN,N5,N4)
      TOH2ER(N2,N1)=TOH2ER(N2,N1)-XOH2ER(N,NN,N5,N4)
      TH1PER(N2,N1)=TH1PER(N2,N1)-XH1PER(N,NN,N5,N4)
      TH2PER(N2,N1)=TH2PER(N2,N1)-XH2PER(N,NN,N5,N4)
      TOH0EB(N2,N1)=TOH0EB(N2,N1)-XOH0EB(N,NN,N5,N4)
      TOH1EB(N2,N1)=TOH1EB(N2,N1)-XOH1EB(N,NN,N5,N4)
      TOH2EB(N2,N1)=TOH2EB(N2,N1)-XOH2EB(N,NN,N5,N4)
      TH1PEB(N2,N1)=TH1PEB(N2,N1)-XH1PEB(N,NN,N5,N4)
      TH2PEB(N2,N1)=TH2PEB(N2,N1)-XH2PEB(N,NN,N5,N4)
C
C     PRECIPITATES FROM EROSION
C
C     sediment code
C       :PALO,PFEO=precip AlOH,FeOH
C       :PCAC,PCAS=precip CaCO3,CaSO4
C       :PALP,PFEP=precip AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precip AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite in non-band
C       :PCPMB,PCPDB,PCPHB=precip CaH2PO4,CaHPO4,apatite in band
C
      TALOER(N2,N1)=TALOER(N2,N1)-PALOER(N,NN,N5,N4)
      TFEOER(N2,N1)=TFEOER(N2,N1)-PFEOER(N,NN,N5,N4)
      TCACER(N2,N1)=TCACER(N2,N1)-PCACER(N,NN,N5,N4)
      TCASER(N2,N1)=TCASER(N2,N1)-PCASER(N,NN,N5,N4)
      TALPER(N2,N1)=TALPER(N2,N1)-PALPER(N,NN,N5,N4)
      TFEPER(N2,N1)=TFEPER(N2,N1)-PFEPER(N,NN,N5,N4)
      TCPDER(N2,N1)=TCPDER(N2,N1)-PCPDER(N,NN,N5,N4)
      TCPHER(N2,N1)=TCPHER(N2,N1)-PCPHER(N,NN,N5,N4)
      TCPMER(N2,N1)=TCPMER(N2,N1)-PCPMER(N,NN,N5,N4)
      TALPEB(N2,N1)=TALPEB(N2,N1)-PALPEB(N,NN,N5,N4)
      TFEPEB(N2,N1)=TFEPEB(N2,N1)-PFEPEB(N,NN,N5,N4)
      TCPDEB(N2,N1)=TCPDEB(N2,N1)-PCPDEB(N,NN,N5,N4)
      TCPHEB(N2,N1)=TCPHEB(N2,N1)-PCPHEB(N,NN,N5,N4)
      TCPMEB(N2,N1)=TCPMEB(N2,N1)-PCPMEB(N,NN,N5,N4)
C
C     ORGANIC MATTER FROM EROSION
C
C     sediment code
C        :OMC,OMN,OMP=microbial C,N,P
C        :ORC,ORN,ORP=microbial residue C,N,P
C        :OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
C        :OSC,OSA,OSN,OSP=SOC,colonized SOC,SON,SOP
C           (K=0:woody litter, K=1:non-woody litter,
C            K=2:manure, K=3:POC, K=4:humus)
C
      DO 7380 K=0,5
      DO 7380 NO=1,7
      DO 7380 M=1,3
      TOMCER(M,NO,K,N2,N1)=TOMCER(M,NO,K,N2,N1)
     2-OMCER(M,NO,K,N,NN,N5,N4)
      TOMNER(M,NO,K,N2,N1)=TOMNER(M,NO,K,N2,N1)
     2-OMNER(M,NO,K,N,NN,N5,N4)
      TOMPER(M,NO,K,N2,N1)=TOMPER(M,NO,K,N2,N1)
     2-OMPER(M,NO,K,N,NN,N5,N4)
7380  CONTINUE
      DO 7375 K=0,4
      DO 7370 M=1,2
      TORCER(M,K,N2,N1)=TORCER(M,K,N2,N1)-ORCER(M,K,N,NN,N5,N4)
      TORNER(M,K,N2,N1)=TORNER(M,K,N2,N1)-ORNER(M,K,N,NN,N5,N4)
      TORPER(M,K,N2,N1)=TORPER(M,K,N2,N1)-ORPER(M,K,N,NN,N5,N4)
7370  CONTINUE
      TOHCER(K,N2,N1)=TOHCER(K,N2,N1)-OHCER(K,N,NN,N5,N4)
      TOHNER(K,N2,N1)=TOHNER(K,N2,N1)-OHNER(K,N,NN,N5,N4)
      TOHPER(K,N2,N1)=TOHPER(K,N2,N1)-OHPER(K,N,NN,N5,N4)
      TOHAER(K,N2,N1)=TOHAER(K,N2,N1)-OHAER(K,N,NN,N5,N4)
      DO 7365 M=1,5
      TOSCER(M,K,N2,N1)=TOSCER(M,K,N2,N1)-OSCER(M,K,N,NN,N5,N4)
      TOSAER(M,K,N2,N1)=TOSAER(M,K,N2,N1)-OSAER(M,K,N,NN,N5,N4)
      TOSNER(M,K,N2,N1)=TOSNER(M,K,N2,N1)-OSNER(M,K,N,NN,N5,N4)
      TOSPER(M,K,N2,N1)=TOSPER(M,K,N2,N1)-OSPER(M,K,N,NN,N5,N4)
7365  CONTINUE
7375  CONTINUE
C     ENDIF
      ENDIF
C
C     OUTGOING FLUXES FROM EROSION W AND N
C
      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      IF(ABS(XSEDER(N,NN,N5B,N4B)).GT.ZEROS(N5,N4))THEN
C
C     SEDIMENT FROM EROSION
C
C     sediment code:XSED=total,XSAN=sand,XSIL=silt,XCLA=clay
C                  :CEC=cation exchange capacity
C                  :AEC=anion exchange capacity
C
      TSEDER(N2,N1)=TSEDER(N2,N1)-XSEDER(N,NN,N5B,N4B)
      TSANER(N2,N1)=TSANER(N2,N1)-XSANER(N,NN,N5B,N4B)
      TSILER(N2,N1)=TSILER(N2,N1)-XSILER(N,NN,N5B,N4B)
      TCLAER(N2,N1)=TCLAER(N2,N1)-XCLAER(N,NN,N5B,N4B)
      TCECER(N2,N1)=TCECER(N2,N1)-XCECER(N,NN,N5B,N4B)
      TAECER(N2,N1)=TAECER(N2,N1)-XAECER(N,NN,N5B,N4B)
C     IF(N1.EQ.1)THEN
C     WRITE(*,1137)'TSEDER3',I,J,NFZ,N,NN,N1,N2,N4B,N5B
C    2,TSEDER(N2,N1),XSEDER(N,NN,N5B,N4B)
1137  FORMAT(A8,9I4,12E12.4)
C     ENDIF
C
C     FERTILIZER FROM EROSION
C
C     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3 in non-band
C                  :NH4B,NH3B,NHUB,NO3B=NH4,NH3,urea,NO3 in band
C
      TNH4ER(N2,N1)=TNH4ER(N2,N1)-XNH4ER(N,NN,N5B,N4B)
      TNH3ER(N2,N1)=TNH3ER(N2,N1)-XNH3ER(N,NN,N5B,N4B)
      TNHUER(N2,N1)=TNHUER(N2,N1)-XNHUER(N,NN,N5B,N4B)
      TNO3ER(N2,N1)=TNO3ER(N2,N1)-XNO3ER(N,NN,N5B,N4B)
      TNH4EB(N2,N1)=TNH4EB(N2,N1)-XNH4EB(N,NN,N5B,N4B)
      TNH3EB(N2,N1)=TNH3EB(N2,N1)-XNH3EB(N,NN,N5B,N4B)
      TNHUEB(N2,N1)=TNHUEB(N2,N1)-XNHUEB(N,NN,N5B,N4B)
      TNO3EB(N2,N1)=TNO3EB(N2,N1)-XNO3EB(N,NN,N5B,N4B)
C
C     EXCHANGEABLE CATIONS AND ANIONS FROM EROSION
C
C     sediment code
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
C           =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
C       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C
      TN4ER(N2,N1)=TN4ER(N2,N1)-XN4ER(N,NN,N5B,N4B)
      TNBER(N2,N1)=TNBER(N2,N1)-XNBER(N,NN,N5B,N4B)
      THYER(N2,N1)=THYER(N2,N1)-XHYER(N,NN,N5B,N4B)
      TALER(N2,N1)=TALER(N2,N1)-XALER(N,NN,N5B,N4B)
      TFEER(N2,N1)=TFEER(N2,N1)-XFEER(N,NN,N5B,N4B)
      TCAER(N2,N1)=TCAER(N2,N1)-XCAER(N,NN,N5B,N4B)
      TMGER(N2,N1)=TMGER(N2,N1)-XMGER(N,NN,N5B,N4B)
      TNAER(N2,N1)=TNAER(N2,N1)-XNAER(N,NN,N5B,N4B)
      TKAER(N2,N1)=TKAER(N2,N1)-XKAER(N,NN,N5B,N4B)
      THCER(N2,N1)=THCER(N2,N1)-XHCER(N,NN,N5B,N4B)
      TAL2ER(N2,N1)=TAL2ER(N2,N1)-XAL2ER(N,NN,N5B,N4B)
      TFE2ER(N2,N1)=TFE2ER(N2,N1)-XFE2ER(N,NN,N5B,N4B)
      TOH0ER(N2,N1)=TOH0ER(N2,N1)-XOH0ER(N,NN,N5B,N4B)
      TOH1ER(N2,N1)=TOH1ER(N2,N1)-XOH1ER(N,NN,N5B,N4B)
      TOH2ER(N2,N1)=TOH2ER(N2,N1)-XOH2ER(N,NN,N5B,N4B)
      TH1PER(N2,N1)=TH1PER(N2,N1)-XH1PER(N,NN,N5B,N4B)
      TH2PER(N2,N1)=TH2PER(N2,N1)-XH2PER(N,NN,N5B,N4B)
      TOH0EB(N2,N1)=TOH0EB(N2,N1)-XOH0EB(N,NN,N5B,N4B)
      TOH1EB(N2,N1)=TOH1EB(N2,N1)-XOH1EB(N,NN,N5B,N4B)
      TOH2EB(N2,N1)=TOH2EB(N2,N1)-XOH2EB(N,NN,N5B,N4B)
      TH1PEB(N2,N1)=TH1PEB(N2,N1)-XH1PEB(N,NN,N5B,N4B)
      TH2PEB(N2,N1)=TH2PEB(N2,N1)-XH2PEB(N,NN,N5B,N4B)
C
C     PRECIPITATES FROM EROSION
C
C     sediment code
C       :PALO,PFEO=precip AlOH,FeOH
C       :PCAC,PCAS=precip CaCO3,CaSO4
C       :PALP,PFEP=precip AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precip AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite in non-band
C       :PCPMB,PCPDB,PCPHB=precip CaH2PO4,CaHPO4,apatite in band
C
      TALOER(N2,N1)=TALOER(N2,N1)-PALOER(N,NN,N5B,N4B)
      TFEOER(N2,N1)=TFEOER(N2,N1)-PFEOER(N,NN,N5B,N4B)
      TCACER(N2,N1)=TCACER(N2,N1)-PCACER(N,NN,N5B,N4B)
      TCASER(N2,N1)=TCASER(N2,N1)-PCASER(N,NN,N5B,N4B)
      TALPER(N2,N1)=TALPER(N2,N1)-PALPER(N,NN,N5B,N4B)
      TFEPER(N2,N1)=TFEPER(N2,N1)-PFEPER(N,NN,N5B,N4B)
      TCPDER(N2,N1)=TCPDER(N2,N1)-PCPDER(N,NN,N5B,N4B)
      TCPHER(N2,N1)=TCPHER(N2,N1)-PCPHER(N,NN,N5B,N4B)
      TCPMER(N2,N1)=TCPMER(N2,N1)-PCPMER(N,NN,N5B,N4B)
      TALPEB(N2,N1)=TALPEB(N2,N1)-PALPEB(N,NN,N5B,N4B)
      TFEPEB(N2,N1)=TFEPEB(N2,N1)-PFEPEB(N,NN,N5B,N4B)
      TCPDEB(N2,N1)=TCPDEB(N2,N1)-PCPDEB(N,NN,N5B,N4B)
      TCPHEB(N2,N1)=TCPHEB(N2,N1)-PCPHEB(N,NN,N5B,N4B)
      TCPMEB(N2,N1)=TCPMEB(N2,N1)-PCPMEB(N,NN,N5B,N4B)
C
C     ORGANIC MATTER FROM EROSION
C
C     sediment code
C        :OMC,OMN,OMP=microbial C,N,P
C        :ORC,ORN,ORP=microbial residue C,N,P
C        :OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
C        :OSC,OSA,OSN,OSP=SOC,colonized SOC,SON,SOP
C           (K=0:woody litter, K=1:non-woody litter,
C            K=2:manure, K=3:POC, K=4:humus)
C
      DO 8380 K=0,5
      DO 8380 NO=1,7
      DO 8380 M=1,3
      TOMCER(M,NO,K,N2,N1)=TOMCER(M,NO,K,N2,N1)
     2-OMCER(M,NO,K,N,NN,N5B,N4B)
      TOMNER(M,NO,K,N2,N1)=TOMNER(M,NO,K,N2,N1)
     2-OMNER(M,NO,K,N,NN,N5B,N4B)
      TOMPER(M,NO,K,N2,N1)=TOMPER(M,NO,K,N2,N1)
     2-OMPER(M,NO,K,N,NN,N5B,N4B)
8380  CONTINUE
      DO 8375 K=0,4
      DO 8370 M=1,2
      TORCER(M,K,N2,N1)=TORCER(M,K,N2,N1)-ORCER(M,K,N,NN,N5B,N4B)
      TORNER(M,K,N2,N1)=TORNER(M,K,N2,N1)-ORNER(M,K,N,NN,N5B,N4B)
      TORPER(M,K,N2,N1)=TORPER(M,K,N2,N1)-ORPER(M,K,N,NN,N5B,N4B)
8370  CONTINUE
      TOHCER(K,N2,N1)=TOHCER(K,N2,N1)-OHCER(K,N,NN,N5B,N4B)
      TOHNER(K,N2,N1)=TOHNER(K,N2,N1)-OHNER(K,N,NN,N5B,N4B)
      TOHPER(K,N2,N1)=TOHPER(K,N2,N1)-OHPER(K,N,NN,N5B,N4B)
      TOHAER(K,N2,N1)=TOHAER(K,N2,N1)-OHAER(K,N,NN,N5B,N4B)
      DO 8365 M=1,5
      TOSCER(M,K,N2,N1)=TOSCER(M,K,N2,N1)-OSCER(M,K,N,NN,N5B,N4B)
      TOSAER(M,K,N2,N1)=TOSAER(M,K,N2,N1)-OSAER(M,K,N,NN,N5B,N4B)
      TOSNER(M,K,N2,N1)=TOSNER(M,K,N2,N1)-OSNER(M,K,N,NN,N5B,N4B)
      TOSPER(M,K,N2,N1)=TOSPER(M,K,N2,N1)-OSPER(M,K,N,NN,N5B,N4B)
8365  CONTINUE
8375  CONTINUE
      ENDIF
      ENDIF
9350  CONTINUE
      ENDIF
      ENDIF
C
C     NET HEAT, WATER FLUXES BETWEEN ADJACENT
C     GRID CELLS
C
C     TFLW,TFLV,TFLWH,THFLW=net micropore,macropore water,vapor flux,
C        heat flux
C     FLW,FLV,FLWH,HFLW=micropore,macropore water,vapor flux,
C        heat flux from watsub.f
C     FLWNU,FLWHNU,HFLWNU=lake surface water flux, heat flux
C        from watsub.f if lake surface disappears
C
      IF(NCN(N2,N1).NE.3.OR.N.EQ.3)THEN
      DO 1200 LL=N6,NL(N5,N4)
      IF(VOLX(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
      N6=LL
      GO TO 1201
      ENDIF
1200  CONTINUE
1201  CONTINUE
      IF(VOLX(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      IF(N3.EQ.NU(N2,N1).AND.N.EQ.3)THEN
      TFLW(N3,N2,N1)=TFLW(N3,N2,N1)+FLW(N,N3,N2,N1)-FLWNU(N5,N4)
      TFLV(N3,N2,N1)=TFLV(N3,N2,N1)+FLV(N,N3,N2,N1)-FLVNU(N5,N4)
      TFLWX(N3,N2,N1)=TFLWX(N3,N2,N1)+FLWX(N,N3,N2,N1)-FLWXNU(N5,N4)
      TFLWH(N3,N2,N1)=TFLWH(N3,N2,N1)+FLWH(N,N3,N2,N1)-FLWHNU(N5,N4)
      THFLW(N3,N2,N1)=THFLW(N3,N2,N1)+HFLW(N,N3,N2,N1)-HFLWNU(N5,N4)
      ELSE
      TFLW(N3,N2,N1)=TFLW(N3,N2,N1)+FLW(N,N3,N2,N1)-FLW(N,N6,N5,N4)
      TFLV(N3,N2,N1)=TFLV(N3,N2,N1)+FLV(N,N3,N2,N1)-FLV(N,N6,N5,N4)
      TFLWX(N3,N2,N1)=TFLWX(N3,N2,N1)+FLWX(N,N3,N2,N1)
     2-FLWX(N,N6,N5,N4)
      TFLWH(N3,N2,N1)=TFLWH(N3,N2,N1)+FLWH(N,N3,N2,N1)
     2-FLWH(N,N6,N5,N4)
      THFLW(N3,N2,N1)=THFLW(N3,N2,N1)+HFLW(N,N3,N2,N1)
     2-HFLW(N,N6,N5,N4)
      ENDIF
C     IF(N1.EQ.1.AND.N3.EQ.1)THEN
C     WRITE(*,6632)'TFLW',I,J,N,N1,N2,N3,N4,N5,N6,NU(N2,N1)
C    2,TFLW(N3,N2,N1),FLW(N,N3,N2,N1),FLW(N,N6,N5,N4),FLWNU(N5,N4)
C    3,THFLW(N3,N2,N1),HFLW(N,N3,N2,N1),HFLW(N,N6,N5,N4)
C    2,HFLWNU(N5,N4),VOLW(N3,N2,N1)
6632  FORMAT(A8,10I4,12E16.8)
C     ENDIF
C
C     NET SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
C
C     T*FLS,T*FHS=net convective+diffusive solute flux through
C        micropores,macropores
C     X*FLS,X*FHS=convective+diffusive solute flux through
C        micropores, macropores from trnsfr.f
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,H1P=HPO4,H2P=H2PO4 in non-band
C        :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,H1B=HPO4,H2B=H2PO4 in band
C
      DO 8585 K=0,4
      TOCFLS(K,N3,N2,N1)=TOCFLS(K,N3,N2,N1)+XOCFLS(K,N,N3,N2,N1)
     2-XOCFLS(K,N,N6,N5,N4)
      TONFLS(K,N3,N2,N1)=TONFLS(K,N3,N2,N1)+XONFLS(K,N,N3,N2,N1)
     2-XONFLS(K,N,N6,N5,N4)
      TOPFLS(K,N3,N2,N1)=TOPFLS(K,N3,N2,N1)+XOPFLS(K,N,N3,N2,N1)
     2-XOPFLS(K,N,N6,N5,N4)
      TOAFLS(K,N3,N2,N1)=TOAFLS(K,N3,N2,N1)+XOAFLS(K,N,N3,N2,N1)
     2-XOAFLS(K,N,N6,N5,N4)
      TOCFHS(K,N3,N2,N1)=TOCFHS(K,N3,N2,N1)+XOCFHS(K,N,N3,N2,N1)
     2-XOCFHS(K,N,N6,N5,N4)
      TONFHS(K,N3,N2,N1)=TONFHS(K,N3,N2,N1)+XONFHS(K,N,N3,N2,N1)
     2-XONFHS(K,N,N6,N5,N4)
      TOPFHS(K,N3,N2,N1)=TOPFHS(K,N3,N2,N1)+XOPFHS(K,N,N3,N2,N1)
     2-XOPFHS(K,N,N6,N5,N4)
      TOAFHS(K,N3,N2,N1)=TOAFHS(K,N3,N2,N1)+XOAFHS(K,N,N3,N2,N1)
     2-XOAFHS(K,N,N6,N5,N4)
C     WRITE(*,2629)'TOCFLS',I,J,NFZ,N1,N2,N3,N4,N5,N6,N,K
C    2,TOCFLS(K,N3,N2,N1),XOCFLS(K,N,N3,N2,N1),XOCFLS(K,N,N6,N5,N4)
2629  FORMAT(A8,11I4,20E12.4)
8585  CONTINUE
      TCOFLS(N3,N2,N1)=TCOFLS(N3,N2,N1)+XCOFLS(N,N3,N2,N1)
     2-XCOFLS(N,N6,N5,N4)
      TCHFLS(N3,N2,N1)=TCHFLS(N3,N2,N1)+XCHFLS(N,N3,N2,N1)
     2-XCHFLS(N,N6,N5,N4)
      TOXFLS(N3,N2,N1)=TOXFLS(N3,N2,N1)+XOXFLS(N,N3,N2,N1)
     2-XOXFLS(N,N6,N5,N4)
      TNGFLS(N3,N2,N1)=TNGFLS(N3,N2,N1)+XNGFLS(N,N3,N2,N1)
     2-XNGFLS(N,N6,N5,N4)
      TN2FLS(N3,N2,N1)=TN2FLS(N3,N2,N1)+XN2FLS(N,N3,N2,N1)
     2-XN2FLS(N,N6,N5,N4)
      THGFLS(N3,N2,N1)=THGFLS(N3,N2,N1)+XHGFLS(N,N3,N2,N1)
     2-XHGFLS(N,N6,N5,N4)
      TN4FLS(N3,N2,N1)=TN4FLS(N3,N2,N1)+XN4FLW(N,N3,N2,N1)
     2-XN4FLW(N,N6,N5,N4)
      TN3FLS(N3,N2,N1)=TN3FLS(N3,N2,N1)+XN3FLW(N,N3,N2,N1)
     2-XN3FLW(N,N6,N5,N4)
      TNOFLS(N3,N2,N1)=TNOFLS(N3,N2,N1)+XNOFLW(N,N3,N2,N1)
     2-XNOFLW(N,N6,N5,N4)
      TNXFLS(N3,N2,N1)=TNXFLS(N3,N2,N1)+XNXFLS(N,N3,N2,N1)
     2-XNXFLS(N,N6,N5,N4)
      TP1FLS(N3,N2,N1)=TP1FLS(N3,N2,N1)+XH1PFS(N,N3,N2,N1)
     2-XH1PFS(N,N6,N5,N4)
      TPOFLS(N3,N2,N1)=TPOFLS(N3,N2,N1)+XH2PFS(N,N3,N2,N1)
     2-XH2PFS(N,N6,N5,N4)
      TN4FLB(N3,N2,N1)=TN4FLB(N3,N2,N1)+XN4FLB(N,N3,N2,N1)
     2-XN4FLB(N,N6,N5,N4)
      TN3FLB(N3,N2,N1)=TN3FLB(N3,N2,N1)+XN3FLB(N,N3,N2,N1)
     2-XN3FLB(N,N6,N5,N4)
      TNOFLB(N3,N2,N1)=TNOFLB(N3,N2,N1)+XNOFLB(N,N3,N2,N1)
     2-XNOFLB(N,N6,N5,N4)
      TNXFLB(N3,N2,N1)=TNXFLB(N3,N2,N1)+XNXFLB(N,N3,N2,N1)
     2-XNXFLB(N,N6,N5,N4)
      TH1BFB(N3,N2,N1)=TH1BFB(N3,N2,N1)+XH1BFB(N,N3,N2,N1)
     2-XH1BFB(N,N6,N5,N4)
      TH2BFB(N3,N2,N1)=TH2BFB(N3,N2,N1)+XH2BFB(N,N3,N2,N1)
     2-XH2BFB(N,N6,N5,N4)
      TCOFHS(N3,N2,N1)=TCOFHS(N3,N2,N1)+XCOFHS(N,N3,N2,N1)
     2-XCOFHS(N,N6,N5,N4)
      TCHFHS(N3,N2,N1)=TCHFHS(N3,N2,N1)+XCHFHS(N,N3,N2,N1)
     2-XCHFHS(N,N6,N5,N4)
      TOXFHS(N3,N2,N1)=TOXFHS(N3,N2,N1)+XOXFHS(N,N3,N2,N1)
     2-XOXFHS(N,N6,N5,N4)
      TNGFHS(N3,N2,N1)=TNGFHS(N3,N2,N1)+XNGFHS(N,N3,N2,N1)
     2-XNGFHS(N,N6,N5,N4)
      TN2FHS(N3,N2,N1)=TN2FHS(N3,N2,N1)+XN2FHS(N,N3,N2,N1)
     2-XN2FHS(N,N6,N5,N4)
      THGFHS(N3,N2,N1)=THGFHS(N3,N2,N1)+XHGFHS(N,N3,N2,N1)
     2-XHGFHS(N,N6,N5,N4)
      TN4FHS(N3,N2,N1)=TN4FHS(N3,N2,N1)+XN4FHW(N,N3,N2,N1)
     2-XN4FHW(N,N6,N5,N4)
      TN3FHS(N3,N2,N1)=TN3FHS(N3,N2,N1)+XN3FHW(N,N3,N2,N1)
     2-XN3FHW(N,N6,N5,N4)
      TNOFHS(N3,N2,N1)=TNOFHS(N3,N2,N1)+XNOFHW(N,N3,N2,N1)
     2-XNOFHW(N,N6,N5,N4)
      TNXFHS(N3,N2,N1)=TNXFHS(N3,N2,N1)+XNXFHS(N,N3,N2,N1)
     2-XNXFHS(N,N6,N5,N4)
      TP1FHS(N3,N2,N1)=TP1FHS(N3,N2,N1)+XH1PHS(N,N3,N2,N1)
     2-XH1PHS(N,N6,N5,N4)
      TPOFHS(N3,N2,N1)=TPOFHS(N3,N2,N1)+XH2PHS(N,N3,N2,N1)
     2-XH2PHS(N,N6,N5,N4)
      TN4FHB(N3,N2,N1)=TN4FHB(N3,N2,N1)+XN4FHB(N,N3,N2,N1)
     2-XN4FHB(N,N6,N5,N4)
      TN3FHB(N3,N2,N1)=TN3FHB(N3,N2,N1)+XN3FHB(N,N3,N2,N1)
     2-XN3FHB(N,N6,N5,N4)
      TNOFHB(N3,N2,N1)=TNOFHB(N3,N2,N1)+XNOFHB(N,N3,N2,N1)
     2-XNOFHB(N,N6,N5,N4)
      TNXFHB(N3,N2,N1)=TNXFHB(N3,N2,N1)+XNXFHB(N,N3,N2,N1)
     2-XNXFHB(N,N6,N5,N4)
      TH1BHB(N3,N2,N1)=TH1BHB(N3,N2,N1)+XH1BHB(N,N3,N2,N1)
     2-XH1BHB(N,N6,N5,N4)
      TH2BHB(N3,N2,N1)=TH2BHB(N3,N2,N1)+XH2BHB(N,N3,N2,N1)
     2-XH2BHB(N,N6,N5,N4)
C
C     NET GAS FLUXES BETWEEN ADJACENT GRID CELLS
C
C     T*FLG=net convective+diffusive gas flux
C     X*FLG=convective+diffusive gas flux from trnsfr.f
C     gas code:*VP*=vapor,*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2
C             :*Z2O*=N2O,*ZN3*=NH3,*H2G*=H2
C
      TVPFLG(N3,N2,N1)=TVPFLG(N3,N2,N1)+XVPFLG(N,N3,N2,N1)
     2-XVPFLG(N,N6,N5,N4)
      TCOFLG(N3,N2,N1)=TCOFLG(N3,N2,N1)+XCOFLG(N,N3,N2,N1)
     2-XCOFLG(N,N6,N5,N4)
      TCHFLG(N3,N2,N1)=TCHFLG(N3,N2,N1)+XCHFLG(N,N3,N2,N1)
     2-XCHFLG(N,N6,N5,N4)
      TOXFLG(N3,N2,N1)=TOXFLG(N3,N2,N1)+XOXFLG(N,N3,N2,N1)
     2-XOXFLG(N,N6,N5,N4)
      TNGFLG(N3,N2,N1)=TNGFLG(N3,N2,N1)+XNGFLG(N,N3,N2,N1)
     2-XNGFLG(N,N6,N5,N4)
      TN2FLG(N3,N2,N1)=TN2FLG(N3,N2,N1)+XN2FLG(N,N3,N2,N1)
     2-XN2FLG(N,N6,N5,N4)
      TNHFLG(N3,N2,N1)=TNHFLG(N3,N2,N1)+XN3FLG(N,N3,N2,N1)
     2-XN3FLG(N,N6,N5,N4)
      THGFLG(N3,N2,N1)=THGFLG(N3,N2,N1)+XHGFLG(N,N3,N2,N1)
     2-XHGFLG(N,N6,N5,N4)
C
C     NET SALT FLUXES BETWEEN ADJACENT GRID CELLS
C
C     T*FLS,T*FHS=net convective+diffusive solute flux through
C        micropores,macropores
C     X*FLS,X*FHS=convective+diffusive solute flux through
C        micropores, macropores from trnsfrs.f
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C        :*1=non-band,*B=band
C
      IF(ISALTG.NE.0)THEN
      TALFLS(N3,N2,N1)=TALFLS(N3,N2,N1)+XALFLS(N,N3,N2,N1)
     2-XALFLS(N,N6,N5,N4)
      TFEFLS(N3,N2,N1)=TFEFLS(N3,N2,N1)+XFEFLS(N,N3,N2,N1)
     2-XFEFLS(N,N6,N5,N4)
      THYFLS(N3,N2,N1)=THYFLS(N3,N2,N1)+XHYFLS(N,N3,N2,N1)
     2-XHYFLS(N,N6,N5,N4)
      TCAFLS(N3,N2,N1)=TCAFLS(N3,N2,N1)+XCAFLS(N,N3,N2,N1)
     2-XCAFLS(N,N6,N5,N4)
      TMGFLS(N3,N2,N1)=TMGFLS(N3,N2,N1)+XMGFLS(N,N3,N2,N1)
     2-XMGFLS(N,N6,N5,N4)
      TNAFLS(N3,N2,N1)=TNAFLS(N3,N2,N1)+XNAFLS(N,N3,N2,N1)
     2-XNAFLS(N,N6,N5,N4)
      TKAFLS(N3,N2,N1)=TKAFLS(N3,N2,N1)+XKAFLS(N,N3,N2,N1)
     2-XKAFLS(N,N6,N5,N4)
      TOHFLS(N3,N2,N1)=TOHFLS(N3,N2,N1)+XOHFLS(N,N3,N2,N1)
     2-XOHFLS(N,N6,N5,N4)
      TSOFLS(N3,N2,N1)=TSOFLS(N3,N2,N1)+XSOFLS(N,N3,N2,N1)
     2-XSOFLS(N,N6,N5,N4)
      TCLFLS(N3,N2,N1)=TCLFLS(N3,N2,N1)+XCLFLS(N,N3,N2,N1)
     2-XCLFLS(N,N6,N5,N4)
      TC3FLS(N3,N2,N1)=TC3FLS(N3,N2,N1)+XC3FLS(N,N3,N2,N1)
     2-XC3FLS(N,N6,N5,N4)
      THCFLS(N3,N2,N1)=THCFLS(N3,N2,N1)+XHCFLS(N,N3,N2,N1)
     2-XHCFLS(N,N6,N5,N4)
      TAL1FS(N3,N2,N1)=TAL1FS(N3,N2,N1)+XAL1FS(N,N3,N2,N1)
     2-XAL1FS(N,N6,N5,N4)
      TAL2FS(N3,N2,N1)=TAL2FS(N3,N2,N1)+XAL2FS(N,N3,N2,N1)
     2-XAL2FS(N,N6,N5,N4)
      TAL3FS(N3,N2,N1)=TAL3FS(N3,N2,N1)+XAL3FS(N,N3,N2,N1)
     2-XAL3FS(N,N6,N5,N4)
      TAL4FS(N3,N2,N1)=TAL4FS(N3,N2,N1)+XAL4FS(N,N3,N2,N1)
     2-XAL4FS(N,N6,N5,N4)
      TALSFS(N3,N2,N1)=TALSFS(N3,N2,N1)+XALSFS(N,N3,N2,N1)
     2-XALSFS(N,N6,N5,N4)
      TFE1FS(N3,N2,N1)=TFE1FS(N3,N2,N1)+XFE1FS(N,N3,N2,N1)
     2-XFE1FS(N,N6,N5,N4)
      TFE2FS(N3,N2,N1)=TFE2FS(N3,N2,N1)+XFE2FS(N,N3,N2,N1)
     2-XFE2FS(N,N6,N5,N4)
      TFE3FS(N3,N2,N1)=TFE3FS(N3,N2,N1)+XFE3FS(N,N3,N2,N1)
     2-XFE3FS(N,N6,N5,N4)
      TFE4FS(N3,N2,N1)=TFE4FS(N3,N2,N1)+XFE4FS(N,N3,N2,N1)
     2-XFE4FS(N,N6,N5,N4)
      TFESFS(N3,N2,N1)=TFESFS(N3,N2,N1)+XFESFS(N,N3,N2,N1)
     2-XFESFS(N,N6,N5,N4)
      TCAOFS(N3,N2,N1)=TCAOFS(N3,N2,N1)+XCAOFS(N,N3,N2,N1)
     2-XCAOFS(N,N6,N5,N4)
      TCACFS(N3,N2,N1)=TCACFS(N3,N2,N1)+XCACFS(N,N3,N2,N1)
     2-XCACFS(N,N6,N5,N4)
      TCAHFS(N3,N2,N1)=TCAHFS(N3,N2,N1)+XCAHFS(N,N3,N2,N1)
     2-XCAHFS(N,N6,N5,N4)
      TCASFS(N3,N2,N1)=TCASFS(N3,N2,N1)+XCASFS(N,N3,N2,N1)
     2-XCASFS(N,N6,N5,N4)
      TMGOFS(N3,N2,N1)=TMGOFS(N3,N2,N1)+XMGOFS(N,N3,N2,N1)
     2-XMGOFS(N,N6,N5,N4)
      TMGCFS(N3,N2,N1)=TMGCFS(N3,N2,N1)+XMGCFS(N,N3,N2,N1)
     2-XMGCFS(N,N6,N5,N4)
      TMGHFS(N3,N2,N1)=TMGHFS(N3,N2,N1)+XMGHFS(N,N3,N2,N1)
     2-XMGHFS(N,N6,N5,N4)
      TMGSFS(N3,N2,N1)=TMGSFS(N3,N2,N1)+XMGSFS(N,N3,N2,N1)
     2-XMGSFS(N,N6,N5,N4)
      TNACFS(N3,N2,N1)=TNACFS(N3,N2,N1)+XNACFS(N,N3,N2,N1)
     2-XNACFS(N,N6,N5,N4)
      TNASFS(N3,N2,N1)=TNASFS(N3,N2,N1)+XNASFS(N,N3,N2,N1)
     2-XNASFS(N,N6,N5,N4)
      TKASFS(N3,N2,N1)=TKASFS(N3,N2,N1)+XKASFS(N,N3,N2,N1)
     2-XKASFS(N,N6,N5,N4)
      TH0PFS(N3,N2,N1)=TH0PFS(N3,N2,N1)+XH0PFS(N,N3,N2,N1)
     2-XH0PFS(N,N6,N5,N4)
      TH3PFS(N3,N2,N1)=TH3PFS(N3,N2,N1)+XH3PFS(N,N3,N2,N1)
     2-XH3PFS(N,N6,N5,N4)
      TF1PFS(N3,N2,N1)=TF1PFS(N3,N2,N1)+XF1PFS(N,N3,N2,N1)
     2-XF1PFS(N,N6,N5,N4)
      TF2PFS(N3,N2,N1)=TF2PFS(N3,N2,N1)+XF2PFS(N,N3,N2,N1)
     2-XF2PFS(N,N6,N5,N4)
      TC0PFS(N3,N2,N1)=TC0PFS(N3,N2,N1)+XC0PFS(N,N3,N2,N1)
     2-XC0PFS(N,N6,N5,N4)
      TC1PFS(N3,N2,N1)=TC1PFS(N3,N2,N1)+XC1PFS(N,N3,N2,N1)
     2-XC1PFS(N,N6,N5,N4)
      TC2PFS(N3,N2,N1)=TC2PFS(N3,N2,N1)+XC2PFS(N,N3,N2,N1)
     2-XC2PFS(N,N6,N5,N4)
      TM1PFS(N3,N2,N1)=TM1PFS(N3,N2,N1)+XM1PFS(N,N3,N2,N1)
     2-XM1PFS(N,N6,N5,N4)
      TH0BFB(N3,N2,N1)=TH0BFB(N3,N2,N1)+XH0BFB(N,N3,N2,N1)
     2-XH0BFB(N,N6,N5,N4)
      TH3BFB(N3,N2,N1)=TH3BFB(N3,N2,N1)+XH3BFB(N,N3,N2,N1)
     2-XH3BFB(N,N6,N5,N4)
      TF1BFB(N3,N2,N1)=TF1BFB(N3,N2,N1)+XF1BFB(N,N3,N2,N1)
     2-XF1BFB(N,N6,N5,N4)
      TF2BFB(N3,N2,N1)=TF2BFB(N3,N2,N1)+XF2BFB(N,N3,N2,N1)
     2-XF2BFB(N,N6,N5,N4)
      TC0BFB(N3,N2,N1)=TC0BFB(N3,N2,N1)+XC0BFB(N,N3,N2,N1)
     2-XC0BFB(N,N6,N5,N4)
      TC1BFB(N3,N2,N1)=TC1BFB(N3,N2,N1)+XC1BFB(N,N3,N2,N1)
     2-XC1BFB(N,N6,N5,N4)
      TC2BFB(N3,N2,N1)=TC2BFB(N3,N2,N1)+XC2BFB(N,N3,N2,N1)
     2-XC2BFB(N,N6,N5,N4)
      TM1BFB(N3,N2,N1)=TM1BFB(N3,N2,N1)+XM1BFB(N,N3,N2,N1)
     2-XM1BFB(N,N6,N5,N4)
      TALFHS(N3,N2,N1)=TALFHS(N3,N2,N1)+XALFHS(N,N3,N2,N1)
     2-XALFHS(N,N6,N5,N4)
      TFEFHS(N3,N2,N1)=TFEFHS(N3,N2,N1)+XFEFHS(N,N3,N2,N1)
     2-XFEFHS(N,N6,N5,N4)
      THYFHS(N3,N2,N1)=THYFHS(N3,N2,N1)+XHYFHS(N,N3,N2,N1)
     2-XHYFHS(N,N6,N5,N4)
      TCAFHS(N3,N2,N1)=TCAFHS(N3,N2,N1)+XCAFHS(N,N3,N2,N1)
     2-XCAFHS(N,N6,N5,N4)
      TMGFHS(N3,N2,N1)=TMGFHS(N3,N2,N1)+XMGFHS(N,N3,N2,N1)
     2-XMGFHS(N,N6,N5,N4)
      TNAFHS(N3,N2,N1)=TNAFHS(N3,N2,N1)+XNAFHS(N,N3,N2,N1)
     2-XNAFHS(N,N6,N5,N4)
      TKAFHS(N3,N2,N1)=TKAFHS(N3,N2,N1)+XKAFHS(N,N3,N2,N1)
     2-XKAFHS(N,N6,N5,N4)
      TOHFHS(N3,N2,N1)=TOHFHS(N3,N2,N1)+XOHFHS(N,N3,N2,N1)
     2-XOHFHS(N,N6,N5,N4)
      TSOFHS(N3,N2,N1)=TSOFHS(N3,N2,N1)+XSOFHS(N,N3,N2,N1)
     2-XSOFHS(N,N6,N5,N4)
      TCLFHS(N3,N2,N1)=TCLFHS(N3,N2,N1)+XCLFHS(N,N3,N2,N1)
     2-XCLFHS(N,N6,N5,N4)
      TC3FHS(N3,N2,N1)=TC3FHS(N3,N2,N1)+XC3FHS(N,N3,N2,N1)
     2-XC3FHS(N,N6,N5,N4)
      THCFHS(N3,N2,N1)=THCFHS(N3,N2,N1)+XHCFHS(N,N3,N2,N1)
     2-XHCFHS(N,N6,N5,N4)
      TAL1HS(N3,N2,N1)=TAL1HS(N3,N2,N1)+XAL1HS(N,N3,N2,N1)
     2-XAL1HS(N,N6,N5,N4)
      TAL2HS(N3,N2,N1)=TAL2HS(N3,N2,N1)+XAL2HS(N,N3,N2,N1)
     2-XAL2HS(N,N6,N5,N4)
      TAL3HS(N3,N2,N1)=TAL3HS(N3,N2,N1)+XAL3HS(N,N3,N2,N1)
     2-XAL3HS(N,N6,N5,N4)
      TAL4HS(N3,N2,N1)=TAL4HS(N3,N2,N1)+XAL4HS(N,N3,N2,N1)
     2-XAL4HS(N,N6,N5,N4)
      TALSHS(N3,N2,N1)=TALSHS(N3,N2,N1)+XALSHS(N,N3,N2,N1)
     2-XALSHS(N,N6,N5,N4)
      TFE1HS(N3,N2,N1)=TFE1HS(N3,N2,N1)+XFE1HS(N,N3,N2,N1)
     2-XFE1HS(N,N6,N5,N4)
      TFE2HS(N3,N2,N1)=TFE2HS(N3,N2,N1)+XFE2HS(N,N3,N2,N1)
     2-XFE2HS(N,N6,N5,N4)
      TFE3HS(N3,N2,N1)=TFE3HS(N3,N2,N1)+XFE3HS(N,N3,N2,N1)
     2-XFE3HS(N,N6,N5,N4)
      TFE4HS(N3,N2,N1)=TFE4HS(N3,N2,N1)+XFE4HS(N,N3,N2,N1)
     2-XFE4HS(N,N6,N5,N4)
      TFESHS(N3,N2,N1)=TFESHS(N3,N2,N1)+XFESHS(N,N3,N2,N1)
     2-XFESHS(N,N6,N5,N4)
      TCAOHS(N3,N2,N1)=TCAOHS(N3,N2,N1)+XCAOHS(N,N3,N2,N1)
     2-XCAOHS(N,N6,N5,N4)
      TCACHS(N3,N2,N1)=TCACHS(N3,N2,N1)+XCACHS(N,N3,N2,N1)
     2-XCACHS(N,N6,N5,N4)
      TCAHHS(N3,N2,N1)=TCAHHS(N3,N2,N1)+XCAHHS(N,N3,N2,N1)
     2-XCAHHS(N,N6,N5,N4)
      TCASHS(N3,N2,N1)=TCASHS(N3,N2,N1)+XCASHS(N,N3,N2,N1)
     2-XCASHS(N,N6,N5,N4)
      TMGOHS(N3,N2,N1)=TMGOHS(N3,N2,N1)+XMGOHS(N,N3,N2,N1)
     2-XMGOHS(N,N6,N5,N4)
      TMGCHS(N3,N2,N1)=TMGCHS(N3,N2,N1)+XMGCHS(N,N3,N2,N1)
     2-XMGCHS(N,N6,N5,N4)
      TMGHHS(N3,N2,N1)=TMGHHS(N3,N2,N1)+XMGHHS(N,N3,N2,N1)
     2-XMGHHS(N,N6,N5,N4)
      TMGSHS(N3,N2,N1)=TMGSHS(N3,N2,N1)+XMGSHS(N,N3,N2,N1)
     2-XMGSHS(N,N6,N5,N4)
      TNACHS(N3,N2,N1)=TNACHS(N3,N2,N1)+XNACHS(N,N3,N2,N1)
     2-XNACHS(N,N6,N5,N4)
      TNASHS(N3,N2,N1)=TNASHS(N3,N2,N1)+XNASHS(N,N3,N2,N1)
     2-XNASHS(N,N6,N5,N4)
      TKASHS(N3,N2,N1)=TKASHS(N3,N2,N1)+XKASHS(N,N3,N2,N1)
     2-XKASHS(N,N6,N5,N4)
      TH0PHS(N3,N2,N1)=TH0PHS(N3,N2,N1)+XH0PHS(N,N3,N2,N1)
     2-XH0PHS(N,N6,N5,N4)
      TH3PHS(N3,N2,N1)=TH3PHS(N3,N2,N1)+XH3PHS(N,N3,N2,N1)
     2-XH3PHS(N,N6,N5,N4)
      TF1PHS(N3,N2,N1)=TF1PHS(N3,N2,N1)+XF1PHS(N,N3,N2,N1)
     2-XF1PHS(N,N6,N5,N4)
      TF2PHS(N3,N2,N1)=TF2PHS(N3,N2,N1)+XF2PHS(N,N3,N2,N1)
     2-XF2PHS(N,N6,N5,N4)
      TC0PHS(N3,N2,N1)=TC0PHS(N3,N2,N1)+XC0PHS(N,N3,N2,N1)
     2-XC0PHS(N,N6,N5,N4)
      TC1PHS(N3,N2,N1)=TC1PHS(N3,N2,N1)+XC1PHS(N,N3,N2,N1)
     2-XC1PHS(N,N6,N5,N4)
      TC2PHS(N3,N2,N1)=TC2PHS(N3,N2,N1)+XC2PHS(N,N3,N2,N1)
     2-XC2PHS(N,N6,N5,N4)
      TM1PHS(N3,N2,N1)=TM1PHS(N3,N2,N1)+XM1PHS(N,N3,N2,N1)
     2-XM1PHS(N,N6,N5,N4)
      TH0BHB(N3,N2,N1)=TH0BHB(N3,N2,N1)+XH0BHB(N,N3,N2,N1)
     2-XH0BHB(N,N6,N5,N4)
      TH3BHB(N3,N2,N1)=TH3BHB(N3,N2,N1)+XH3BHB(N,N3,N2,N1)
     2-XH3BHB(N,N6,N5,N4)
      TF1BHB(N3,N2,N1)=TF1BHB(N3,N2,N1)+XF1BHB(N,N3,N2,N1)
     2-XF1BHB(N,N6,N5,N4)
      TF2BHB(N3,N2,N1)=TF2BHB(N3,N2,N1)+XF2BHB(N,N3,N2,N1)
     2-XF2BHB(N,N6,N5,N4)
      TC0BHB(N3,N2,N1)=TC0BHB(N3,N2,N1)+XC0BHB(N,N3,N2,N1)
     2-XC0BHB(N,N6,N5,N4)
      TC1BHB(N3,N2,N1)=TC1BHB(N3,N2,N1)+XC1BHB(N,N3,N2,N1)
     2-XC1BHB(N,N6,N5,N4)
      TC2BHB(N3,N2,N1)=TC2BHB(N3,N2,N1)+XC2BHB(N,N3,N2,N1)
     2-XC2BHB(N,N6,N5,N4)
      TM1BHB(N3,N2,N1)=TM1BHB(N3,N2,N1)+XM1BHB(N,N3,N2,N1)
     2-XM1BHB(N,N6,N5,N4)
      ENDIF
      ELSE
      TFLW(N3,N2,N1)=0.0
      TFLV(N3,N2,N1)=0.0
      TFLWX(N3,N2,N1)=0.0
      TFLWH(N3,N2,N1)=0.0
      THFLW(N3,N2,N1)=0.0
      TWFLFL(N3,N2,N1)=0.0
      TWFLFH(N3,N2,N1)=0.0
      THFLFL(N3,N2,N1)=0.0
      TWFLVL(N3,N2,N1)=0.0
      THFLVL(N3,N2,N1)=0.0
      DO 8596 K=0,4
      TOCFLS(K,N3,N2,N1)=0.0
      TONFLS(K,N3,N2,N1)=0.0
      TOPFLS(K,N3,N2,N1)=0.0
      TOAFLS(K,N3,N2,N1)=0.0
      TOCFHS(K,N3,N2,N1)=0.0
      TONFHS(K,N3,N2,N1)=0.0
      TOPFHS(K,N3,N2,N1)=0.0
      TOAFHS(K,N3,N2,N1)=0.0
8596  CONTINUE
      TCOFLS(N3,N2,N1)=0.0
      TCHFLS(N3,N2,N1)=0.0
      TOXFLS(N3,N2,N1)=0.0
      TNGFLS(N3,N2,N1)=0.0
      TN2FLS(N3,N2,N1)=0.0
      THGFLS(N3,N2,N1)=0.0
      TN4FLS(N3,N2,N1)=0.0
      TN3FLS(N3,N2,N1)=0.0
      TNOFLS(N3,N2,N1)=0.0
      TNXFLS(N3,N2,N1)=0.0
      TP1FLS(N3,N2,N1)=0.0
      TPOFLS(N3,N2,N1)=0.0
      TN4FLB(N3,N2,N1)=0.0
      TN3FLB(N3,N2,N1)=0.0
      TNOFLB(N3,N2,N1)=0.0
      TNXFLB(N3,N2,N1)=0.0
      TH1BFB(N3,N2,N1)=0.0
      TH2BFB(N3,N2,N1)=0.0
      TCOFHS(N3,N2,N1)=0.0
      TCHFHS(N3,N2,N1)=0.0
      TOXFHS(N3,N2,N1)=0.0
      TNGFHS(N3,N2,N1)=0.0
      TN2FHS(N3,N2,N1)=0.0
      THGFHS(N3,N2,N1)=0.0
      TN4FHS(N3,N2,N1)=0.0
      TN3FHS(N3,N2,N1)=0.0
      TNOFHS(N3,N2,N1)=0.0
      TNXFHS(N3,N2,N1)=0.0
      TP1FHS(N3,N2,N1)=0.0
      TPOFHS(N3,N2,N1)=0.0
      TN4FHB(N3,N2,N1)=0.0
      TN3FHB(N3,N2,N1)=0.0
      TNOFHB(N3,N2,N1)=0.0
      TNXFHB(N3,N2,N1)=0.0
      TH1BHB(N3,N2,N1)=0.0
      TH2BHB(N3,N2,N1)=0.0
      TVPFLG(N3,N2,N1)=0.0
      TCOFLG(N3,N2,N1)=0.0
      TCHFLG(N3,N2,N1)=0.0
      TOXFLG(N3,N2,N1)=0.0
      TNGFLG(N3,N2,N1)=0.0
      TN2FLG(N3,N2,N1)=0.0
      TNHFLG(N3,N2,N1)=0.0
      THGFLG(N3,N2,N1)=0.0
      IF(ISALTG.NE.0)THEN
      TALFLS(N3,N2,N1)=0.0
      TFEFLS(N3,N2,N1)=0.0
      THYFLS(N3,N2,N1)=0.0
      TCAFLS(N3,N2,N1)=0.0
      TMGFLS(N3,N2,N1)=0.0
      TNAFLS(N3,N2,N1)=0.0
      TKAFLS(N3,N2,N1)=0.0
      TOHFLS(N3,N2,N1)=0.0
      TSOFLS(N3,N2,N1)=0.0
      TCLFLS(N3,N2,N1)=0.0
      TC3FLS(N3,N2,N1)=0.0
      THCFLS(N3,N2,N1)=0.0
      TAL1FS(N3,N2,N1)=0.0
      TAL2FS(N3,N2,N1)=0.0
      TAL3FS(N3,N2,N1)=0.0
      TAL4FS(N3,N2,N1)=0.0
      TALSFS(N3,N2,N1)=0.0
      TFE1FS(N3,N2,N1)=0.0
      TFE2FS(N3,N2,N1)=0.0
      TFE3FS(N3,N2,N1)=0.0
      TFE4FS(N3,N2,N1)=0.0
      TFESFS(N3,N2,N1)=0.0
      TCAOFS(N3,N2,N1)=0.0
      TCACFS(N3,N2,N1)=0.0
      TCAHFS(N3,N2,N1)=0.0
      TCASFS(N3,N2,N1)=0.0
      TMGOFS(N3,N2,N1)=0.0
      TMGCFS(N3,N2,N1)=0.0
      TMGHFS(N3,N2,N1)=0.0
      TMGSFS(N3,N2,N1)=0.0
      TNACFS(N3,N2,N1)=0.0
      TNASFS(N3,N2,N1)=0.0
      TKASFS(N3,N2,N1)=0.0
      TH0PFS(N3,N2,N1)=0.0
      TH3PFS(N3,N2,N1)=0.0
      TF1PFS(N3,N2,N1)=0.0
      TF2PFS(N3,N2,N1)=0.0
      TC0PFS(N3,N2,N1)=0.0
      TC1PFS(N3,N2,N1)=0.0
      TC2PFS(N3,N2,N1)=0.0
      TM1PFS(N3,N2,N1)=0.0
      TH0BFB(N3,N2,N1)=0.0
      TH3BFB(N3,N2,N1)=0.0
      TF1BFB(N3,N2,N1)=0.0
      TF2BFB(N3,N2,N1)=0.0
      TC0BFB(N3,N2,N1)=0.0
      TC1BFB(N3,N2,N1)=0.0
      TC2BFB(N3,N2,N1)=0.0
      TM1BFB(N3,N2,N1)=0.0
      TALFHS(N3,N2,N1)=0.0
      TFEFHS(N3,N2,N1)=0.0
      THYFHS(N3,N2,N1)=0.0
      TCAFHS(N3,N2,N1)=0.0
      TMGFHS(N3,N2,N1)=0.0
      TNAFHS(N3,N2,N1)=0.0
      TKAFHS(N3,N2,N1)=0.0
      TOHFHS(N3,N2,N1)=0.0
      TSOFHS(N3,N2,N1)=0.0
      TCLFHS(N3,N2,N1)=0.0
      TC3FHS(N3,N2,N1)=0.0
      THCFHS(N3,N2,N1)=0.0
      TAL1HS(N3,N2,N1)=0.0
      TAL2HS(N3,N2,N1)=0.0
      TAL3HS(N3,N2,N1)=0.0
      TAL4HS(N3,N2,N1)=0.0
      TALSHS(N3,N2,N1)=0.0
      TFE1HS(N3,N2,N1)=0.0
      TFE2HS(N3,N2,N1)=0.0
      TFE3HS(N3,N2,N1)=0.0
      TFE4HS(N3,N2,N1)=0.0
      TFESHS(N3,N2,N1)=0.0
      TCAOHS(N3,N2,N1)=0.0
      TCACHS(N3,N2,N1)=0.0
      TCAHHS(N3,N2,N1)=0.0
      TCASHS(N3,N2,N1)=0.0
      TMGOHS(N3,N2,N1)=0.0
      TMGCHS(N3,N2,N1)=0.0
      TMGHHS(N3,N2,N1)=0.0
      TMGSHS(N3,N2,N1)=0.0
      TNACHS(N3,N2,N1)=0.0
      TNASHS(N3,N2,N1)=0.0
      TKASHS(N3,N2,N1)=0.0
      TH0PHS(N3,N2,N1)=0.0
      TH3PHS(N3,N2,N1)=0.0
      TF1PHS(N3,N2,N1)=0.0
      TF2PHS(N3,N2,N1)=0.0
      TC0PHS(N3,N2,N1)=0.0
      TC1PHS(N3,N2,N1)=0.0
      TC2PHS(N3,N2,N1)=0.0
      TM1PHS(N3,N2,N1)=0.0
      TH0BHB(N3,N2,N1)=0.0
      TH3BHB(N3,N2,N1)=0.0
      TF1BHB(N3,N2,N1)=0.0
      TF2BHB(N3,N2,N1)=0.0
      TC0BHB(N3,N2,N1)=0.0
      TC1BHB(N3,N2,N1)=0.0
      TC2BHB(N3,N2,N1)=0.0
      TM1BHB(N3,N2,N1)=0.0
      ENDIF
      ENDIF
      ENDIF
8580  CONTINUE
C
C     NET FREEZE-THAW
C
C     XWFLFL,XWFLFH=total freeze-thaw flux in micropores,macropores
C        from watsub.f
C     XHFLFL=total freeze-thaw latent heat flux
C        from watsub.f
C     XWFLVL=total evaporation-condensation
C        from watsub.f
C     XHFLVL=total evaporation-condensation latent heat flux
C        from watsub.f
C     TWFLFL,TWFLFH=net freeze-thaw in micropores,macropores
C     THFLFL=net freeze-thaw latent heat flux
C     TWFLVL=net evaporation-condensation
C     THFLVL=net evaporation-condensation latent heat flux
C
      TWFLFL(N3,N2,N1)=TWFLFL(N3,N2,N1)+XWFLFL(N3,N2,N1)
      TWFLFH(N3,N2,N1)=TWFLFH(N3,N2,N1)+XWFLFH(N3,N2,N1)
      THFLFL(N3,N2,N1)=THFLFL(N3,N2,N1)+XHFLFL(N3,N2,N1)
      TWFLVL(N3,N2,N1)=TWFLVL(N3,N2,N1)+XWFLVL(N3,N2,N1)
      THFLVL(N3,N2,N1)=THFLVL(N3,N2,N1)+XHFLVL(N3,N2,N1)
8575  CONTINUE
C
C     CALCULATE SNOWPACK TEMPERATURE FROM ITS CHANGE
C     IN HEAT STORAGE
C
      VOLSS(NY,NX)=0.0
      VOLWS(NY,NX)=0.0
      VOLIS(NY,NX)=0.0
      VOLS(NY,NX)=0.0
      DPTHS(NY,NX)=0.0
      VOLSWI=0.0
      DO 9780 L=1,JS
C
C     ADD CHANGES IN SNOW, WATER AND ICE
C
C     VOLSSL,VOLWSL,VOLVSL,VOLISL=snow water equivalent,water,vapor,
C        ice volume in snowpack layer
C     TFLWS,TFLWW,TFLWV,TFLWI=net fluxes of snow,water,vapor,ice
C        in snowpack
C     XWFLFS,XWFLFI,XWFLVS=freeze-thaw,vapor flux from watsub.f
C     DENSI=ice density from starts.f
C
      VOLSSL(L,NY,NX)=VOLSSL(L,NY,NX)+TFLWS(L,NY,NX)-XWFLFS(L,NY,NX)
     2+XWFLVS(L,NY,NX)
      VOLWSL(L,NY,NX)=VOLWSL(L,NY,NX)+TFLWW(L,NY,NX)+XWFLFS(L,NY,NX)
     2+XWFLFI(L,NY,NX)+XWFLVW(L,NY,NX)
      VOLVSL(L,NY,NX)=VOLVSL(L,NY,NX)+TFLWV(L,NY,NX)-XWFLVS(L,NY,NX)
     2-XWFLVW(L,NY,NX)
      VOLISL(L,NY,NX)=VOLISL(L,NY,NX)+TFLWI(L,NY,NX)
     2-XWFLFI(L,NY,NX)/DENSI
C
C     ACCUMULATE SNOW MASS FOR CALCULATING COMPRESSION
C
C     VOLWSI=accumulated water equivalent volume in snowpack
C     XFLWS=snow transfer from watsub.f
C     TCASF,DENSF,VOLSF=snowfall temperature,density,volume
C     DENSS=snow density in layer
C
      IF(L.EQ.1)THEN
      VOLSWI=VOLSWI+0.5*(VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX)
     2+VOLISL(L,NY,NX)*DENSI)
C
C     RESET SNOW SURFACE DENSITY FOR SNOWFALL
C
      IF(XFLWS(L,NY,NX).GT.0.0)THEN
      DENSX=DENSS(L,NY,NX)
      TCASF=AMAX1(-15.0,AMIN1(2.0,TCA(NY,NX)))
      DENSF=0.05+1.7E-03*(TCASF+15.0)**1.5
      VOLSF=XFLWS(L,NY,NX)/DENSF
     2+(VOLSSL(L,NY,NX)-XFLWS(L,NY,NX))/DENSS(L,NY,NX)
      DENSS(L,NY,NX)=VOLSSL(L,NY,NX)/VOLSF
      ENDIF
      ELSE
      VOLSWI=VOLSWI+0.5*(VOLSSL(L-1,NY,NX)+VOLWSL(L-1,NY,NX)
     2+VOLISL(L-1,NY,NX)*DENSI+VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX)
     2+VOLISL(L,NY,NX)*DENSI)
      ENDIF
C
C     SNOWPACK COMPRESSION
C
C     DDENS1,DDENS2=temperature, compression effect on snow density
C     DENSS=snow density in snowpack layer
C     VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume
C        in snowpack layer
C     VOLSL=snowpack layer volume
C     DLYRS=snowpack layer depth
C     CDPTHS=cumulative depth to bottom of snowpack layer
C     VHCPW=snowpack layer heat capacity
C     TKW,TCW=snowpack layer temperature K,oC
C     THFLWW=convective heat fluxes of snow,water,ice in snowpack
C     XTHAWW=latent heat flux from freeze-thaw from watsub.f
C     HEATIN=cumulative net surface heat transfer
C     VOLSS,VOLWS,VOLVS,VOLIS=total snow water equivalent,water,vapor
C        ice content of snowpack
C     VOLS,DPTHS=total snowpack volume, depth
C
      IF(DENSS(L,NY,NX).LT.0.25)THEN
      DDENS1=DENSS(L,NY,NX)*1.0E-05*EXP(0.04*TCW(L,NY,NX))
      ELSE
      DDENS1=0.0
      ENDIF
      CVISC=0.25*EXP(-0.08*TCW(L,NY,NX)+23.0*DENSS(L,NY,NX))
      DDENS2=DENSS(L,NY,NX)*VOLSWI/(AREA(3,NU(NY,NX),NY,NX)*CVISC)
      DENSS(L,NY,NX)=DENSS(L,NY,NX)+(DDENS1+DDENS2)*XNFH
      IF(VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX)
     2.GT.ZEROS2(NY,NX))THEN
      VOLSL(L,NY,NX)=VOLSSL(L,NY,NX)/DENSS(L,NY,NX)
     2+VOLWSL(L,NY,NX)+VOLISL(L,NY,NX)
      DLYRS(L,NY,NX)=AMAX1(0.0,VOLSL(L,NY,NX))
     2/AREA(3,NU(NY,NX),NY,NX)
      CDPTHS(L,NY,NX)=CDPTHS(L-1,NY,NX)+DLYRS(L,NY,NX)
      VHCPWZ(L,NY,NX)=VHCPW(L,NY,NX)
      TKWX=TKW(L,NY,NX)
      ENGYW=VHCPW(L,NY,NX)*TKW(L,NY,NX)
      VHCPW(L,NY,NX)=2.095*VOLSSL(L,NY,NX)
     2+4.19*(VOLWSL(L,NY,NX)+VOLVSL(L,NY,NX))
     2+1.9274*VOLISL(L,NY,NX)
      IF(VHCPWZ(L,NY,NX).GT.VHCPWX(NY,NX)
     2.AND.VHCPW(L,NY,NX).GT.ZEROS(NY,NX))THEN
      TKW(L,NY,NX)=(ENGYW+THFLWW(L,NY,NX)+XHFLF0(L,NY,NX)
     2+XHFLV0(L,NY,NX))/VHCPW(L,NY,NX)
      ELSE
      IF(L.EQ.1)THEN
      TKW(L,NY,NX)=TKAM(NY,NX)
      ELSE
      TKW(L,NY,NX)=TKW(L-1,NY,NX)
      ENDIF
      IF(VHCPWZ(L,NY,NX).GT.VHCPWX(NY,NX))THEN
      HEATIN=HEATIN+(TKW(L,NY,NX)-TKWX)*VHCPW(L,NY,NX)
      ENDIF
      ENDIF
      VOLSS(NY,NX)=VOLSS(NY,NX)+VOLSSL(L,NY,NX)
      VOLWS(NY,NX)=VOLWS(NY,NX)+VOLWSL(L,NY,NX)
      VOLIS(NY,NX)=VOLIS(NY,NX)+VOLISL(L,NY,NX)
      VOLS(NY,NX)=VOLS(NY,NX)+VOLSL(L,NY,NX)
      DPTHS(NY,NX)=DPTHS(NY,NX)+DLYRS(L,NY,NX)
C     IF(L.EQ.5)THEN
C     WRITE(*,7753)'VOLSS',I,J,NFZ,NX,NY,L,VOLSSL(L,NY,NX)
C    3,VOLWSL(L,NY,NX),VOLVSL(L,NY,NX),VOLISL(L,NY,NX)
C    3,VOLSL(L,NY,NX),XFLWS(L,NY,NX),XFLWW(L,NY,NX),XFLWI(L,NY,NX)
C    3,XFLWS(L+1,NY,NX),XFLWW(L+1,NY,NX),XFLWI(L+1,NY,NX)
C    4,TFLWS(L,NY,NX),TFLWW(L,NY,NX),TFLWV(L,NY,NX),TFLWI(L,NY,NX)
C    5,XWFLFS(L,NY,NX),XWFLFI(L,NY,NX),XWFLVS(L,NY,NX),XWFLVW(L,NY,NX)
C    3,VOLSS(NY,NX),VOLWS(NY,NX),VOLIS(NY,NX),VOLS(NY,NX),DPTHS(NY,NX)
C    5,DLYRS(L,NY,NX),DENSS(L,NY,NX),DDENS1,DDENS2,TCW(L,NY,NX)
C    6,CVISC,VOLSWI,TCA(NY,NX)
C     WRITE(*,7753)'TKW',I,J,NFZ,NX,NY,L,TKW(L,NY,NX),TCW(L,NY,NX)
C    3,THFLWW(L,NY,NX),XHFLF0(L,NY,NX),XHFLV0(L,NY,NX),XHFLWW(L,NY,NX)
C    4,XHFLWW(L+1,NY,NX),HFLSWR(L,NY,NX),HFLSW(L,NY,NX)
C    4,XHFLV0(L,NY,NX)
C    4,VHCPWX(NY,NX),VHCPW(L,NY,NX),VHCPWZ(L,NY,NX)
7753  FORMAT(A8,6I4,100E14.6)
C     ENDIF
      ELSE
      VOLSSL(L,NY,NX)=0.0
      VOLWSL(L,NY,NX)=0.0
      VOLVSL(L,NY,NX)=0.0
      VOLISL(L,NY,NX)=0.0
      VOLSL(L,NY,NX)=0.0
      DLYRS(L,NY,NX)=0.0
      CDPTHS(L,NY,NX)=CDPTHS(L-1,NY,NX)
      VHCPW(L,NY,NX)=0.0
      IF(L.EQ.1)THEN
      TKW(L,NY,NX)=TKAM(NY,NX)
      DENSS(L,NY,NX)=DENS0(NY,NX)
      ELSE
      TKW(L,NY,NX)=TKW(L-1,NY,NX)
      DENSS(L,NY,NX)=DENSS(L-1,NY,NX)
      ENDIF
      ENDIF
      TCW(L,NY,NX)=TKW(L,NY,NX)-273.15
C
C     SNOWPACK SOLUTE CONTENT
C
C     *W2=solute content of snowpack
C     T*BLS=net solute flux in snowpack
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C             :N4=NH4,N3=NH3,NO=NO3,1P=HPO4,HP=H2PO4
C
      CO2W(L,NY,NX)=CO2W(L,NY,NX)+TCOBLS(L,NY,NX)
      CH4W(L,NY,NX)=CH4W(L,NY,NX)+TCHBLS(L,NY,NX)
      OXYW(L,NY,NX)=OXYW(L,NY,NX)+TOXBLS(L,NY,NX)
      ZNGW(L,NY,NX)=ZNGW(L,NY,NX)+TNGBLS(L,NY,NX)
      ZN2W(L,NY,NX)=ZN2W(L,NY,NX)+TN2BLS(L,NY,NX)
      ZN4W(L,NY,NX)=ZN4W(L,NY,NX)+TN4BLW(L,NY,NX)
      ZN3W(L,NY,NX)=ZN3W(L,NY,NX)+TN3BLW(L,NY,NX)
      ZNOW(L,NY,NX)=ZNOW(L,NY,NX)+TNOBLW(L,NY,NX)
      Z1PW(L,NY,NX)=Z1PW(L,NY,NX)+TH1PBS(L,NY,NX)
      ZHPW(L,NY,NX)=ZHPW(L,NY,NX)+TH2PBS(L,NY,NX)
C
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C
      IF(ISALTG.NE.0)THEN
      ZALW(L,NY,NX)=ZALW(L,NY,NX)+TALBLS(L,NY,NX)
      ZFEW(L,NY,NX)=ZFEW(L,NY,NX)+TFEBLS(L,NY,NX)
      ZHYW(L,NY,NX)=ZHYW(L,NY,NX)+THYBLS(L,NY,NX)
      ZCAW(L,NY,NX)=ZCAW(L,NY,NX)+TCABLS(L,NY,NX)
      ZMGW(L,NY,NX)=ZMGW(L,NY,NX)+TMGBLS(L,NY,NX)
      ZNAW(L,NY,NX)=ZNAW(L,NY,NX)+TNABLS(L,NY,NX)
      ZKAW(L,NY,NX)=ZKAW(L,NY,NX)+TKABLS(L,NY,NX)
      ZOHW(L,NY,NX)=ZOHW(L,NY,NX)+TOHBLS(L,NY,NX)
      ZSO4W(L,NY,NX)=ZSO4W(L,NY,NX)+TSOBLS(L,NY,NX)
      ZCLW(L,NY,NX)=ZCLW(L,NY,NX)+TCLBLS(L,NY,NX)
      ZCO3W(L,NY,NX)=ZCO3W(L,NY,NX)+TC3BLS(L,NY,NX)
      ZHCO3W(L,NY,NX)=ZHCO3W(L,NY,NX)+THCBLS(L,NY,NX)
      ZALH1W(L,NY,NX)=ZALH1W(L,NY,NX)+TAL1BS(L,NY,NX)
      ZALH2W(L,NY,NX)=ZALH2W(L,NY,NX)+TAL2BS(L,NY,NX)
      ZALH3W(L,NY,NX)=ZALH3W(L,NY,NX)+TAL3BS(L,NY,NX)
      ZALH4W(L,NY,NX)=ZALH4W(L,NY,NX)+TAL4BS(L,NY,NX)
      ZALSW(L,NY,NX)=ZALSW(L,NY,NX)+TALSBS(L,NY,NX)
      ZFEH1W(L,NY,NX)=ZFEH1W(L,NY,NX)+TFE1BS(L,NY,NX)
      ZFEH2W(L,NY,NX)=ZFEH2W(L,NY,NX)+TFE2BS(L,NY,NX)
      ZFEH3W(L,NY,NX)=ZFEH3W(L,NY,NX)+TFE3BS(L,NY,NX)
      ZFEH4W(L,NY,NX)=ZFEH4W(L,NY,NX)+TFE4BS(L,NY,NX)
      ZFESW(L,NY,NX)=ZFESW(L,NY,NX)+TFESBS(L,NY,NX)
      ZCAOW(L,NY,NX)=ZCAOW(L,NY,NX)+TCAOBS(L,NY,NX)
      ZCACW(L,NY,NX)=ZCACW(L,NY,NX)+TCACBS(L,NY,NX)
      ZCAHW(L,NY,NX)=ZCAHW(L,NY,NX)+TCAHBS(L,NY,NX)
      ZCASW(L,NY,NX)=ZCASW(L,NY,NX)+TCASBS(L,NY,NX)
      ZMGOW(L,NY,NX)=ZMGOW(L,NY,NX)+TMGOBS(L,NY,NX)
      ZMGCW(L,NY,NX)=ZMGCW(L,NY,NX)+TMGCBS(L,NY,NX)
      ZMGHW(L,NY,NX)=ZMGHW(L,NY,NX)+TMGHBS(L,NY,NX)
      ZMGSW(L,NY,NX)=ZMGSW(L,NY,NX)+TMGSBS(L,NY,NX)
      ZNACW(L,NY,NX)=ZNACW(L,NY,NX)+TNACBS(L,NY,NX)
      ZNASW(L,NY,NX)=ZNASW(L,NY,NX)+TNASBS(L,NY,NX)
      ZKASW(L,NY,NX)=ZKASW(L,NY,NX)+TKASBS(L,NY,NX)
      H0PO4W(L,NY,NX)=H0PO4W(L,NY,NX)+TH0PBS(L,NY,NX)
      H3PO4W(L,NY,NX)=H3PO4W(L,NY,NX)+TH3PBS(L,NY,NX)
      ZFE1PW(L,NY,NX)=ZFE1PW(L,NY,NX)+TF1PBS(L,NY,NX)
      ZFE2PW(L,NY,NX)=ZFE2PW(L,NY,NX)+TF2PBS(L,NY,NX)
      ZCA0PW(L,NY,NX)=ZCA0PW(L,NY,NX)+TC0PBS(L,NY,NX)
      ZCA1PW(L,NY,NX)=ZCA1PW(L,NY,NX)+TC1PBS(L,NY,NX)
      ZCA2PW(L,NY,NX)=ZCA2PW(L,NY,NX)+TC2PBS(L,NY,NX)
      ZMG1PW(L,NY,NX)=ZMG1PW(L,NY,NX)+TM1PBS(L,NY,NX)
      ENDIF
9780  CONTINUE
C
C     SNOW RUNOFF
C
C     VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume
C        in snowpack layer
C     TQS,TQW,TQI,THQS=snowpack snow,water,ice,heat runoff
C     VHCPW=snowpack layer heat capacity
C     TKW,TCW=snowpack layer temperature K,oC
C
      VOLSSL(1,NY,NX)=VOLSSL(1,NY,NX)+TQS(NY,NX)
      VOLWSL(1,NY,NX)=VOLWSL(1,NY,NX)+TQW(NY,NX)
      VOLISL(1,NY,NX)=VOLISL(1,NY,NX)+TQI(NY,NX)
      ENGYW=VHCPW(1,NY,NX)*TKW(1,NY,NX)
      VHCPW(1,NY,NX)=2.095*VOLSSL(1,NY,NX)
     2+4.19*(VOLWSL(1,NY,NX)+VOLVSL(1,NY,NX))
     2+1.9274*VOLISL(1,NY,NX)
      IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
      TKW(1,NY,NX)=(ENGYW+THQS(NY,NX))
     2/VHCPW(1,NY,NX)
      ELSE
      TKW(1,NY,NX)=TKAM(NY,NX)
      ENDIF
C     WRITE(*,7752)'TQS',I,J,NX,NY,TKW(1,NY,NX)
C    3,VOLSSL(1,NY,NX),VOLWSL(1,NY,NX),VOLISL(1,NY,NX),VOLSL(1,NY,NX)
C    3,TQS(NY,NX),FLWS(NY,NX),TQW(NY,NX),FLWW(NY,NX)
C    4,TQI(NY,NX),FLWI(NY,NX),THQS(NY,NX),VHCPW(1,NY,NX)
7752  FORMAT(A8,4I4,20E12.4)
C
C     IF SNOWPACK DISAPPEARS
C
C     VHCPW=snowpack layer heat capacity
C     TKW,TCW,TKA=snowpack layer,air temperature K,oC
C     VOLSSL,VOLWSL,VOLVSL,VOLISL=snow water equivalent,water,vapor,
C        ice volume in snowpack layer
C     VOLW,VOLV,VOLI=surface water,vapor,ice content
C
      IF(XHFLWX(NY,NX).GT.ZEROS(NY,NX))THEN
      VOLSSL(1,NY,NX)=VOLSSL(1,NY,NX)-XFLWSX(NY,NX)
      VOLWSL(1,NY,NX)=VOLWSL(1,NY,NX)-XFLWWX(NY,NX)
      VOLVSL(1,NY,NX)=VOLVSL(1,NY,NX)-XFLWVX(NY,NX)
      VOLISL(1,NY,NX)=VOLISL(1,NY,NX)-XFLWIX(NY,NX)
      VOLW(0,NY,NX)=VOLW(0,NY,NX)+XFLWWX(NY,NX)
      VOLV(0,NY,NX)=VOLV(0,NY,NX)+XFLWVX(NY,NX)
      VOLI(0,NY,NX)=VOLI(0,NY,NX)+XFLWIX(NY,NX)+XFLWSX(NY,NX)/DENSI
      ENGYW=TKW(1,NY,NX)*VHCPW(1,NY,NX)
      VHCPW(1,NY,NX)=2.095*VOLSSL(1,NY,NX)
     2+4.19*(VOLWSL(1,NY,NX)+VOLVSL(1,NY,NX))
     2+1.9274*VOLISL(1,NY,NX)
      IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
      TKW(1,NY,NX)=(ENGYW-XHFLWX(NY,NX))/VHCPW(1,NY,NX)
      ELSE
      TKW(1,NY,NX)=TKQ(NY,NX)
      ENDIF
      ENGY0=TKS(0,NY,NX)*VHCP(0,NY,NX)
      VHCP(0,NY,NX)=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))
     2+4.19*(VOLW(0,NY,NX)+VOLV(0,NY,NX))
     3+1.9274*VOLI(0,NY,NX)
      IF(VHCP(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      TKS(0,NY,NX)=(ENGY0+XHFLWX(NY,NX))/VHCP(0,NY,NX)
      ELSE
      TKS(0,NY,NX)=TKS(NUM(NY,NX),NY,NX)
      ENDIF
      VOLSS(NY,NX)=VOLSSL(1,NY,NX)
      VOLWS(NY,NX)=VOLWSL(1,NY,NX)
      VOLIS(NY,NX)=VOLISL(1,NY,NX)
      DO 9770 L=1,JS
      DENSS(L,NY,NX)=DENS0(NY,NX)
9770  CONTINUE
      VOLS(NY,NX)=VOLSSL(1,NY,NX)/DENSS(1,NY,NX)
     2+VOLWSL(1,NY,NX)+VOLISL(1,NY,NX)
      DPTHS(NY,NX)=0.0
C     WRITE(*,2122)'SNW',I,J,NFZ,NX,NY
C    2,XFLWSX(NY,NX),XFLWWX(NY,NX),XFLWIX(NY,NX),XHFLWX(NY,NX)
C    3,VOLSSL(1,NY,NX),VOLWSL(1,NY,NX),VOLISL(1,NY,NX)
C    4,VOLVSL(1,NY,NX),VOLW(0,NY,NX),VOLV(0,NY,NX),VOLI(0,NY,NX)
C    3,VHCP(0,NY,NX),TKS(0,NY,NX),VHCPW(1,NY,NX),TKW(1,NY,NX)
2122  FORMAT(A8,5I4,20E12.4)
      ENDIF
C
C     CALCULATE SURFACE RESIDUE TEMPERATURE FROM ITS CHANGE
C     IN HEAT STORAGE
C
C     VHCP=litter heat capacity
C     VOLW,VOLV,VOLI=litter water,vapor,ice content
C     ORGC,ORGCC=litter SOC,charcoal content
C     FLWR,XWFLVR,XWFLFR,TQR=litter water,evaporation-condensation,
C        freeze-thaw,runoff flux from watsub.f
C     DENSI=ice dnsity from starts.f
C     TKS=surface litter temperature
C     HFLWR,XHFLFR,XHFLVR,THQR=litter water,evaporation-condensation,
C        freeze-thaw,runoff heat flux from watsub.f
C     HCBFX=heat released by litter combustion
C     HEATIN=cumulative surface heat exchange
C
      VHCPZ=VHCP(0,NY,NX)
      VHCPY=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))
     2+4.19*(VOLW(0,NY,NX)+VOLV(0,NY,NX))
     2+1.9274*VOLI(0,NY,NX)
      VHCPO=VHCPY-VHCPZ
      HFLXO=VHCPO*TKS(0,NY,NX)
      VOLW(0,NY,NX)=VOLW(0,NY,NX)+FLWR(NY,NX)+XWFLVR(NY,NX)
     2+XWFLFR(NY,NX)+TQR(NY,NX)
      VOLV(0,NY,NX)=VOLV(0,NY,NX)+FLVR(NY,NX)-XWFLVR(NY,NX)
      VOLI(0,NY,NX)=VOLI(0,NY,NX)-XWFLFR(NY,NX)/DENSI
      ENGYZ=VHCPZ*TKS(0,NY,NX)
      VHCP(0,NY,NX)=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))
     2+4.19*(VOLW(0,NY,NX)+VOLV(0,NY,NX))
     2+1.9274*VOLI(0,NY,NX)
      IF(VHCP(0,NY,NX).GT.VHCPRX(NY,NX))THEN
      TKS(0,NY,NX)=(ENGYZ+HFLWR(NY,NX)+XHFLFR(NY,NX)+XHFLVR(NY,NX)
     2+THQR(NY,NX)+HFLXO+HCBFX(0,NY,NX))/VHCP(0,NY,NX)
      HEATIN=HEATIN+HFLXO+HCBFX(0,NY,NX)
      ELSE
      HEATIN=HEATIN+HFLXO+HCBFX(0,NY,NX)
     2+(TKS(NUM(NY,NX),NY,NX)-TKS(0,NY,NX))*VHCP(0,NY,NX)
      TKS(0,NY,NX)=TKS(NUM(NY,NX),NY,NX)
      ENDIF
      ENGYR=VHCP(0,NY,NX)*TKS(0,NY,NX)
      HEATSO=HEATSO+ENGYR
      HEATIN=HEATIN+XHFLFR(NY,NX)+XHFLVR(NY,NX)
      TCS(0,NY,NX)=TKS(0,NY,NX)-273.15
      TSMX(0,NY,NX)=AMAX1(TSMX(0,NY,NX),TCS(0,NY,NX))
      TSMN(0,NY,NX)=AMIN1(TSMN(0,NY,NX),TCS(0,NY,NX))
C     IF(ICHKF.EQ.1)THEN
C     WRITE(*,6634)'TKR',I,J,NFZ,NX,NY
C    2,VOLW(0,NY,NX),VOLV(0,NY,NX)
C    2,VOLI(0,NY,NX),VOLP(0,NY,NX),FLWR(NY,NX),FLVR(NY,NX)
C    3,XWFLFR(NY,NX),XWFLVR(NY,NX),TQR(NY,NX)
C    4,TVOLWC(NY,NX),TVOLWP(NY,NX)
C    3,TKS(0,NY,NX),ENGYZ,HFLWR(NY,NX),XHFLFR(NY,NX),XHFLVR(NY,NX)
C    4,THQR(NY,NX),HCBFX(0,NY,NX),HFLXO,VHCPO,VHCP(0,NY,NX)
C    5,ORGC(0,NY,NX),ORGCC(0,NY,NX)
C    5,TKS(NUM(NY,NX),NY,NX),VHCPRX(NY,NX)
6634  FORMAT(A8,5I4,30E12.4)
C     ENDIF
C
C     CALCULATE CANOPY AIR TEMPERATURE, VAPOR CONCENTRATION
C
C     TKQ=bulk canopy,ground air temperature
C     TKQG=ground surface air temperature
C     TKQT=bulk canopy air temperature
C     VPQ=bulk canopy,ground vapor pressure
C     VPQG=ground surface vapor pressure
C     VPQT=bulk canopy vapor pressure
C     FRADG=fraction of radiation received at ground surface
C        from hour1.f
C
      TKQ(NY,NX)=TKQG(NY,NX)*FRADG(NY,NX)
     2+TKQT(NY,NX)*(1.0-FRADG(NY,NX))
      TCQ(NY,NX)=TKQ(NY,NX)-273.15
      VPQ(NY,NX)=VPQG(NY,NX)*FRADG(NY,NX)
     2+VPQT(NY,NX)*(1.0-FRADG(NY,NX))
      IF(ICHKF.EQ.1)THEN
      WRITE(*,3114)'TKQ',I,J,NFZ,NX,NY
     2,TKQ(NY,NX),TKAM(NY,NX),VPQ(NY,NX),VPAM(NY,NX)
     3,TKQG(NY,NX),TKQT(NY,NX),(TKQC(NZ,NY,NX),NZ=1,NP(NY,NX))
     4,(TKQD(NZ,NY,NX),NZ=1,NP(NY,NX))
     3,TKS(0,NY,NX),TKCT(NY,NX),(TKC(NZ,NY,NX),NZ=1,NP(NY,NX))
     4,(TKD(NZ,NY,NX),NZ=1,NP(NY,NX))
     3,VPQG(NY,NX),VPQT(NY,NX),(VPQC(NZ,NY,NX),NZ=1,NP(NY,NX))
     4,(VPQD(NZ,NY,NX),NZ=1,NP(NY,NX))
     4,TKS(NU(NY,NX),NY,NX),CO2Q(NY,NX),CH4Q(NY,NX),OXYQ(NY,NX)
3114  FORMAT(A8,5I4,50E12.4)
      ENDIF
C
C     SURFACE BOUNDARY WATER FLUXES
C
C     PRECQ,PRECI,PRECU=rain+snow,surface,subsurface irrigation
C     TEVAPG,TEVAPP=total evaporation,transpiration
C     CRAIN,CEVAP=cumulative precipitation, evapotranspiration
C     VOLWOU=cumulative water transfer through lateral and
C        lower boundaries
C     XVPFLG=aggregated convective+diffusive vapor flux
C        from trnsfr.f
C
      WI=(PRECQ(NY,NX)+PRECI(NY,NX))*XNFH
      CRAIN=CRAIN+WI
      URAIN(NY,NX)=URAIN(NY,NX)+WI
      WO=TEVAPG(NY,NX)+TEVAPP(NY,NX)+XVPFLG(3,NU(NY,NX),NY,NX)
      CEVAP=CEVAP-WO
      UEVAP(NY,NX)=UEVAP(NY,NX)-WO
      VOLWOU=VOLWOU-PRECU(NY,NX)*XNFH
      HVOLO(NY,NX)=HVOLO(NY,NX)-PRECU(NY,NX)*XNFH
      UVOLO(NY,NX)=UVOLO(NY,NX)-PRECU(NY,NX)*XNFH
      UDRAIN(NY,NX)=UDRAIN(NY,NX)+FLW(3,NK(NY,NX),NY,NX)
      TEVPGH(NY,NX)=TEVPGH(NY,NX)+TEVAPG(NY,NX)
     2+XVPFLG(3,NU(NY,NX),NY,NX)
      TEVPPH(NY,NX)=TEVPPH(NY,NX)+TEVAPP(NY,NX)
C
C     SURFACE BOUNDARY HEAT FLUXES
C
C     HEATIN,HEATOU=cumulative net surface,subsurface heat transfer
C     PRECA,PRECW=rain+irrigation,snowfall
C     HEATH,THFLXC=net surface,canopy heat exchange
C     XTHAWW=latent snowpack heat flux
C     PRECU=subsurface irrigation
C     TKA=air temperature
C     XHFLF0,XHFLV0=snowpack latent heat flux from frezze-thaw,
C        evaporation-consensation
C
      HEATIN=HEATIN+(4.19*TKAM(NY,NX)*PRECA(NY,NX)
     2+2.095*TKAM(NY,NX)*PRECW(NY,NX))*XNFH
      HEATIN=HEATIN+HEATH(NY,NX)+THFLXC(NY,NX)
      DO 5150 L=1,JS
      HEATIN=HEATIN+XHFLF0(L,NY,NX)+XHFLV0(L,NY,NX)
5150  CONTINUE
      HEATOU=HEATOU-4.19*TKAM(NY,NX)*PRECU(NY,NX)*XNFH
C
C     SURFACE BOUNDARY CO2, CH4 AND DOC FLUXES
C
C     X*DFS=soil surface gas exchange from trnsfr.f
C     X*FLG=soil surface convective+diffusive gas flux from trnsfr.f
C     FLQGQ,FLQGI=water flux to snowpack or soil surface from
C        rain,irrigation
C     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
C     C*R,C*Q=precipitation,irrigation solute concentrations
C     X*DFG=soil surface gas volatilization-dissolution from trnsfr.f
C     X*DFR=litter gas volatilization from trnsfr.f
C     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
C             :*ZN3*=NH3,*H2G*=H2
C     H*G,U*G=current,cumulative gas exchange
C     CO2GIN,TCOU=cumulative surface,subsurface gas C exchange
C     XCNET,ZCNET=net,cumulative CO2 flux in canopy air from
C        soil+plants
C     XHNET,ZHNET=net,cumulative CH4 flux in canopy air from
C        soil+plants
C     XONET,ZONET=net O2,cumulative flux in canopy air from
C        soil+plants
C
      CI=XCODFS(NY,NX)+XCOFLG(3,NU(NY,NX),NY,NX)+TCO2Z(NY,NX)
     2+XCOFLG(3,0,NY,NX)+XCODFR(NY,NX)
     3+(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCOR(NY,NX)
     4+(FLQGI(NY,NX)+FLQRI(NY,NX))*CCOQ(NY,NX)
      CH=XCHDFS(NY,NX)+XCHFLG(3,NU(NY,NX),NY,NX)+TCH4Z(NY,NX)
     2+XCHFLG(3,0,NY,NX)+XCHDFR(NY,NX)
     3+(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCHR(NY,NX)
     4+(FLQGI(NY,NX)+FLQRI(NY,NX))*CCHQ(NY,NX)
      CO=-PRECU(NY,NX)*CCOQ(NY,NX)*XNFH
      CX=-PRECU(NY,NX)*CCHQ(NY,NX)*XNFH
      XCNET(NY,NX)=XCNET(NY,NX)+CI
      XHNET(NY,NX)=XHNET(NY,NX)+CH
      UCO2G(NY,NX)=UCO2G(NY,NX)+CI
      HCO2G(NY,NX)=HCO2G(NY,NX)+CI
      UCH4G(NY,NX)=UCH4G(NY,NX)+CH
      HCH4G(NY,NX)=HCH4G(NY,NX)+CH
      CO2GIN=CO2GIN+CI+CH
      TCOU=TCOU+CO+CX
C
C     SURFACE BOUNDARY O2 FLUXES
C
C     OXYGIN,OXYGOU=cumulative surface,subsurface gas O2 exchange
C     H2GIN,H2GOU=cumulative surface,subsurface gas H2 exchange
C     RUPOXO=microbial O2 uptake in litter from nitro.f
C     ROGOX,RC4OX=O2-limited litter O2 uptake,O2 and CH4-limited
C        litter CH4 combustion
C
      OI=XOXDFS(NY,NX)+XOXFLG(3,NU(NY,NX),NY,NX)+TOXYZ(NY,NX)
     2+XOXFLG(3,0,NY,NX)+XOXDFR(NY,NX)
     3+(FLQGQ(NY,NX)+FLQRQ(NY,NX))*COXR(NY,NX)
     4+(FLQGI(NY,NX)+FLQRI(NY,NX))*COXQ(NY,NX)
      XONET(NY,NX)=XONET(NY,NX)+OI
      OXYGIN=OXYGIN+OI
      OO=RUPOXO(0,NY,NX)-PRECU(NY,NX)*COXQ(NY,NX)*XNFH
     2+ROGOX(0,NY,NX)+RC4OX(0,NY,NX)*2.667
      OXYGOU=OXYGOU+OO
      UOXYG(NY,NX)=UOXYG(NY,NX)+OI
      HOXYG(NY,NX)=HOXYG(NY,NX)+OI
      HI=XHGDFS(NY,NX)+XHGFLG(3,NU(NY,NX),NY,NX)+TH2GZ(NY,NX)
     2+XHGFLG(3,0,NY,NX)+XHGDFR(NY,NX)
      H2GIN=H2GIN+HI
      HO=RH2GO(0,NY,NX)
      H2GOU=H2GOU+HO
C     IF(I.EQ.256)THEN
C     WRITE(*,6646)'UCO2G',I,J,NFZ,NX,NY
C    2,XCNET(NY,NX),XONET(NY,NX),XHNET(NY,NX)
C    2,CI,OI,CH,HCO2G(NY,NX)
C    2,XCODFS(NY,NX),XCOFLG(3,NU(NY,NX),NY,NX),TCO2Z(NY,NX)
C    4,XCOFLG(3,0,NY,NX),XCODFR(NY,NX),XNFH
C    2,FLQGQ(NY,NX),FLQRQ(NY,NX),CCOR(NY,NX)
C    3,FLQGI(NY,NX),FLQRI(NY,NX),CCOQ(NY,NX)
C    5,(TLCO2P(L,NY,NX),L=1,10)
C     WRITE(*,6646)'UOXYG',I,J,NFZ,NX,NY
C    5,HOXYG(NY,NX),OI
C    6,XOXDFS(NY,NX),XOXFLG(3,NU(NY,NX),NY,NX)
C    4,TOXYZ(NY,NX),XOXDFG(0,NY,NX),XOXDFR(NY,NX)
C    2,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*COXR(NY,NX)
C    3,(FLQGI(NY,NX)+FLQRI(NY,NX))*COXQ(NY,NX)
C    6,(TLOXYP(L,NY,NX),L=1,10)
C     WRITE(*,6646)'UCH4G',I,J,NFZ,NX,NY,UCH4G(NY,NX),CH
C    2,XCHDFS(NY,NX),XCHFLG(3,NU(NY,NX),NY,NX),TCH4Z(NY,NX)
C    2,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCHR(NY,NX)
C    3,(FLQGI(NY,NX)+FLQRI(NY,NX))*CCHQ(NY,NX)
C    4,XCHDFG(0,NY,NX),XCHDFR(NY,NX)
C    5,(TLCH4P(L,NY,NX),L=1,10)
6646  FORMAT(A8,5I4,60E12.4)
C     ENDIF
C
C     SURFACE BOUNDARY N2, N2O, NH3, NH4, NO3, AND DON FLUXES
C
C     X*DFS=soil surface gas exchange from trnsfr.f
C     X*FLG=soil surface convective+diffusive gas flux from trnsfr.f
C     FLQGQ,FLQGI=water flux to snowpack or soil surface from
C        rain,irrigation
C     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
C     C*R,C*Q=precipitation,irrigation solute concentrations
C     X*DFG=soil surface gas volatilization-dissolution from trnsfr.f
C     X*DFR=litter gas volatilization from trnsfr.f
C     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
C             :*ZN3*=NH3,*H2G*=H2
C     TZIN,TZOU=cumulative surface,subsurface aqueous NH4,NH3,NO3
C        exchange
C     ZN2GIN=cumulative surface gas N2,N2O,NH3 exchange
C     ZNGGIN,ZN2OIN,ZNH3IN=cumulative surface gas N2,N2O,NH3 exchange
C
      ZSI=((FLQGQ(NY,NX)+FLQRQ(NY,NX))
     2*(CN4R(NY,NX)+CN3R(NY,NX)+CNOR(NY,NX))
     3+(FLQGI(NY,NX)+FLQRI(NY,NX))
     4*(CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX)))*14.0
      ZXB=(-PRECU(NY,NX)*(CNNQ(NY,NX)+CN2Q(NY,NX))-PRECU(NY,NX)
     2*(CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX))*14.0)*XNFH
      TZIN=TZIN+ZSI
      TZOU=TZOU+ZXB
      ZGI=(FLQGQ(NY,NX)+FLQRQ(NY,NX))*(CNNR(NY,NX)+CN2R(NY,NX))
     2+(FLQGI(NY,NX)+FLQRI(NY,NX))*(CNNQ(NY,NX)+CN2Q(NY,NX))
     3+XNGDFS(NY,NX)+XN2DFS(NY,NX)+XN3DFS(NY,NX)
     2+XNBDFS(NY,NX)+XNGFLG(3,NU(NY,NX),NY,NX)
     2+XN2FLG(3,NU(NY,NX),NY,NX)+XN3FLG(3,NU(NY,NX),NY,NX)
     3+TN2OZ(NY,NX)+TNH3Z(NY,NX)
     6+XN2FLG(3,0,NY,NX)+XNGFLG(3,0,NY,NX)+XN3FLG(3,0,NY,NX)
     7+XNGDFR(NY,NX)+XN2DFR(NY,NX)+XN3DFR(NY,NX)
      ZN2GIN=ZN2GIN+ZGI
      ZDRAIN(NY,NX)=ZDRAIN(NY,NX)+XN4FLW(3,NK(NY,NX),NY,NX)
     2+XN3FLW(3,NK(NY,NX),NY,NX)+XNOFLW(3,NK(NY,NX),NY,NX)
     3+XNXFLS(3,NK(NY,NX),NY,NX)+XN4FLB(3,NK(NY,NX),NY,NX)
     4+XN3FLB(3,NK(NY,NX),NY,NX)+XNOFLB(3,NK(NY,NX),NY,NX)
     5+XNXFLB(3,NK(NY,NX),NY,NX)
      ZNGGIN=XNGDFS(NY,NX)+XNGFLG(3,NU(NY,NX),NY,NX)+XNGDFG(0,NY,NX)
      ZN2OIN=XN2DFS(NY,NX)+XN2FLG(3,NU(NY,NX),NY,NX)+XN2DFG(0,NY,NX)
      ZNH3IN=XN3DFS(NY,NX)+XNBDFS(NY,NX)+XN3FLG(3,NU(NY,NX),NY,NX)
     2+XN3DFG(0,NY,NX)
C     UN2GG(NY,NX)=UN2GG(NY,NX)+ZNGGIN
C     HN2GG(NY,NX)=HN2GG(NY,NX)+ZNGGIN
      UN2OG(NY,NX)=UN2OG(NY,NX)+ZN2OIN
      HN2OG(NY,NX)=HN2OG(NY,NX)+ZN2OIN
      UNH3G(NY,NX)=UNH3G(NY,NX)+ZNH3IN
      HNH3G(NY,NX)=HNH3G(NY,NX)+ZNH3IN
      UN2GS(NY,NX)=UN2GS(NY,NX)+XN2GS(0,NY,NX)
      UH2GG(NY,NX)=UH2GG(NY,NX)+HI
C     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,6644)'HNH3G',I,J,NFZ,NX,NY,NU(NY,NX),HNH3G(NY,NX),ZNH3IN
C    2,XN3DFS(NY,NX),XNBDFS(NY,NX),XN3FLG(3,NU(NY,NX),NY,NX)
C    2,XN3DFG(0,NY,NX)
C     WRITE(*,6644)'ZN2GIN',I,J,NFZ,NX,NY,NU(NY,NX)
C    2,ZN2GIN,XNGDFS(NY,NX)
C    3,XN2DFS(NY,NX),XN3DFS(NY,NX)
C    2,XNBDFS(NY,NX),XNGFLG(3,NU(NY,NX),NY,NX)
C    2,XN2FLG(3,NU(NY,NX),NY,NX)
C    3,XN3FLG(3,NU(NY,NX),NY,NX),TN2OZ(NY,NX),TNH3Z(NY,NX)
C    4,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*(CNNR(NY,NX)+CN2R(NY,NX))
C    5,(FLQGI(NY,NX)+FLQRI(NY,NX))*(CNNQ(NY,NX)+CN2Q(NY,NX))
C    6,XN2DFG(0,NY,NX)+XNGDFG(0,NY,NX),XN3DFG(0,NY,NX)
C    7,XNGDFR(NY,NX)+XN2DFR(NY,NX),XN3DFR(NY,NX)
C     ENDIF
C
C     SURFACE BOUNDARY PO4 AND DOP FLUXES
C
C     TPIN,TPOU=cumulative surface,subsurface aqueous H2PO4,HPO4
C        exchange
C
      PI=31.0*((FLQGQ(NY,NX)+FLQRQ(NY,NX))
     2*(CPOR(NY,NX)+CH1PR(NY,NX))
     3+(FLQGI(NY,NX)+FLQRI(NY,NX))
     4*(CPOQ(I,NY,NX)+CH1PQ(I,NY,NX)))
      PXB=-31.0*PRECU(NY,NX)*(CPOQ(I,NY,NX)+CH1PQ(I,NY,NX))*XNFH
      TPIN=TPIN+PI
      TPOU=TPOU+PXB
      PDRAIN(NY,NX)=PDRAIN(NY,NX)+XH2PFS(3,NK(NY,NX),NY,NX)
     2+XH2BFB(3,NK(NY,NX),NY,NX)+XH1PFS(3,NK(NY,NX),NY,NX)
     2+XH1BFB(3,NK(NY,NX),NY,NX)
C
C     SURFACE BOUNDARY ION FLUXES
C
C     X*DFS=soil surface gas exchange from trnsfr.f
C     X*FLG=soil surface convective+diffusive gas flux from trnsfr.f
C     FLQGQ,FLQGI=water flux to snowpack or soil surface from
C        rain,irrigation
C     FLQRQ,FLQRI=water flux to surface litter from rain,irrigation
C     C*R,C*Q=precipitation,irrigation solute concentrations
C     X*DFG=soil surface gas volatilization-dissolution from trnsfr.f
C     X*DFR=litter gas volatilization from trnsfr.f
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C
      SIN=((FLQGQ(NY,NX)+FLQRQ(NY,NX))
     2*(2.0*CN4R(NY,NX)+CN3R(NY,NX)+CNOR(NY,NX))
     3+(FLQGI(NY,NX)+FLQRI(NY,NX))
     4*(2.0*CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX)))
      SGN=(2.0*(FLQGQ(NY,NX)+FLQRQ(NY,NX))*(CNNR(NY,NX)+CN2R(NY,NX))
     2+2.0*(FLQGI(NY,NX)+FLQRI(NY,NX))*(CNNQ(NY,NX)+CN2Q(NY,NX))
     3+2.0*(XNGDFS(NY,NX)+XN2DFS(NY,NX))+XN3DFS(NY,NX)
     2+XNBDFS(NY,NX)+2.0*(XNGFLG(3,NU(NY,NX),NY,NX)
     2+XN2FLG(3,NU(NY,NX),NY,NX))+XN3FLG(3,NU(NY,NX),NY,NX)
     3+2.0*TN2OZ(NY,NX)+TNH3Z(NY,NX)
     6+2.0*(XN2DFG(0,NY,NX)+XNGDFG(0,NY,NX))+XN3DFG(0,NY,NX)
     7+2.0*(XNGDFR(NY,NX)+XN2DFR(NY,NX))+XN3DFR(NY,NX))/14.0
      SIP=((FLQGQ(NY,NX)+FLQRQ(NY,NX))
     2*(3.0*CPOR(NY,NX)+2.0*CH1PR(NY,NX))
     3+(FLQGI(NY,NX)+FLQRI(NY,NX))
     4*(3.0*CPOQ(I,NY,NX)+2.0*CH1PQ(I,NY,NX)))
      SNB=(-PRECU(NY,NX)*(CNNQ(NY,NX)+CN2Q(NY,NX))-PRECU(NY,NX)
     2*(2.0*CN4Q(I,NY,NX)+CN3Q(I,NY,NX)+CNOQ(I,NY,NX)))*XNFH
      SPB=(-PRECU(NY,NX)*(3.0*CPOQ(I,NY,NX)+2.0*CH1PQ(I,NY,NX)))*XNFH
      SNM0=(2.0*XNH4S(0,NY,NX)+XNO3S(0,NY,NX)+XNO2S(0,NY,NX)
     2-2.0*XN2GS(0,NY,NX))/14.0
      SPM0=(2.0*XH1PS(0,NY,NX)+3.0*XH2PS(0,NY,NX))/31.0
C
C     ACCUMULATE PLANT LITTERFALL FLUXES
C
C     ZCSNC,ZZSNC,ZPSNC=total plant C,N,P litterfall from extract.f
C     UXCSN,UXZFN,UXPSN=cumulative plant C,N,P litterfall
C
      XCSN=XCSN+ZCSNC(NY,NX)
      XZSN=XZSN+ZZSNC(NY,NX)
      XPSN=XPSN+ZPSNC(NY,NX)
      UXCSN(NY,NX)=UXCSN(NY,NX)+ZCSNC(NY,NX)
      UXZSN(NY,NX)=UXZSN(NY,NX)+ZZSNC(NY,NX)
      UXPSN(NY,NX)=UXPSN(NY,NX)+ZPSNC(NY,NX)
C
C     SURFACE BOUNDARY SALT FLUXES FROM RAINFALL AND SURFACE
C     IRRIGATION
C
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C        :*1=non-band,*B=band
C     TIONIN=total salt surface flux
C
      IF(ISALTG.NE.0)THEN
      SIR=PRECQ(NY,NX)*(CALR(NY,NX)+CFER(NY,NX)+CHYR(NY,NX)+CCAR(NY,NX)
     2+CMGR(NY,NX)+CNAR(NY,NX)+CKAR(NY,NX)+COHR(NY,NX)+CSOR(NY,NX)
     3+CCLR(NY,NX)+CC3R(NY,NX)+CH0PR(NY,NX)
     4+2.0*(CHCR(NY,NX)+CAL1R(NY,NX)+CALSR(NY,NX)+CFE1R(NY,NX)
     5+CFESR(NY,NX)+CCAOR(NY,NX)+CCACR(NY,NX)+CCASR(NY,NX)+CMGOR(NY,NX)
     6+CMGCR(NY,NX)+CMGSR(NY,NX)+CNACR(NY,NX)+CNASR(NY,NX)
     7+CKASR(NY,NX)+CC0PR(NY,NX))
     8+3.0*(CAL2R(NY,NX)+CFE2R(NY,NX)+CCAHR(NY,NX)+CMGHR(NY,NX)
     9+CF1PR(NY,NX)+CC1PR(NY,NX)+CM1PR(NY,NX))
     1+4.0*(CAL3R(NY,NX)+CFE3R(NY,NX)+CH3PR(NY,NX)+CF2PR(NY,NX)
     2+CC2PR(NY,NX))
     3+5.0*(CAL4R(NY,NX)+CFE4R(NY,NX)))*XNFH
      SII=PRECI(NY,NX)*(CALQ(I,NY,NX)+CFEQ(I,NY,NX)+CHYQ(I,NY,NX)
     2+CCAQ(I,NY,NX)+CMGQ(I,NY,NX)+CNAQ(I,NY,NX)+CKAQ(I,NY,NX)
     3+COHQ(I,NY,NX)+CSOQ(I,NY,NX)+CCLQ(I,NY,NX)+CC3Q(I,NY,NX)
     4+CH0PQ(I,NY,NX)
     5+2.0*(CHCQ(I,NY,NX)+CAL1Q(I,NY,NX)+CALSQ(I,NY,NX)
     5+CFE1Q(I,NY,NX)+CFESQ(I,NY,NX)+CCAOQ(I,NY,NX)+CCACQ(I,NY,NX)
     6+CCASQ(I,NY,NX)+CMGOQ(I,NY,NX)+CMGCQ(I,NY,NX)+CMGSQ(I,NY,NX)
     7+CNACQ(I,NY,NX)+CNASQ(I,NY,NX)+CKASQ(I,NY,NX)+CC0PQ(I,NY,NX))
     9+3.0*(CAL2Q(I,NY,NX)+CFE2Q(I,NY,NX)+CCAHQ(I,NY,NX)
     9+CMGHQ(I,NY,NX)+CF1PQ(I,NY,NX)+CC1PQ(I,NY,NX)+CM1PQ(I,NY,NX))
     2+4.0*(CAL3Q(I,NY,NX)+CFE3Q(I,NY,NX)
     2+CH3PQ(I,NY,NX)+CF2PQ(I,NY,NX)+CC2PQ(I,NY,NX))
     3+5.0*(CAL4Q(I,NY,NX)+CFE4Q(I,NY,NX)))*XNFH
      TIONIN=TIONIN+SIR+SII
C     WRITE(*,3338)'SIR',I,J,SIR,PRECQ(NY,NX)*XNFH
C    2,CALR(NY,NX),CFER(NY,NX),CHYR(NY,NX),CCAR(NY,NX)
C    2,CMGR(NY,NX),CNAR(NY,NX),CKAR(NY,NX),COHR(NY,NX),CSOR(NY,NX)
C    3,CCLR(NY,NX),CC3R(NY,NX),CH0PR(NY,NX)
C    4,CHCR(NY,NX),CAL1R(NY,NX),CALSR(NY,NX),CFE1R(NY,NX)
C    5,CFESR(NY,NX),CCAOR(NY,NX),CCACR(NY,NX),CCASR(NY,NX),CMGOR(NY,NX)
C    6,CMGCR(NY,NX),CMGSR(NY,NX),CNACR(NY,NX),CNASR(NY,NX)
C    7,CKASR(NY,NX),CC0PR(NY,NX)
C    8,CAL2R(NY,NX),CFE2R(NY,NX),CCAHR(NY,NX),CMGHR(NY,NX)
C    9,CF1PR(NY,NX),CC1PR(NY,NX),CM1PR(NY,NX)
C    1,CAL3R(NY,NX),CFE3R(NY,NX),CH3PR(NY,NX),CF2PR(NY,NX)
C    2,CC2PR(NY,NX),CAL4R(NY,NX),CFE4R(NY,NX)
C
C     SUBSURFACE BOUNDARY SALT FLUXES FROM SUBSURFACE IRRIGATION
C
C     TIONOU=total salt subsurface flux
C
      SBU=-PRECU(NY,NX)*(CALQ(I,NY,NX)+CFEQ(I,NY,NX)+CHYQ(I,NY,NX)
     2+CCAQ(I,NY,NX)+CMGQ(I,NY,NX)+CNAQ(I,NY,NX)+CKAQ(I,NY,NX)
     3+COHQ(I,NY,NX)+CSOQ(I,NY,NX)+CCLQ(I,NY,NX)+CC3Q(I,NY,NX)
     4+CH0PQ(I,NY,NX)
     5+2.0*(CHCQ(I,NY,NX)+CAL1Q(I,NY,NX)+CALSQ(I,NY,NX)
     5+CFE1Q(I,NY,NX)+CFESQ(I,NY,NX)+CCAOQ(I,NY,NX)+CCACQ(I,NY,NX)
     6+CCASQ(I,NY,NX)+CMGOQ(I,NY,NX)+CMGCQ(I,NY,NX)+CMGSQ(I,NY,NX)
     7+CNACQ(I,NY,NX)+CNASQ(I,NY,NX)+CKASQ(I,NY,NX)+CC0PQ(I,NY,NX))
     9+3.0*(CAL2Q(I,NY,NX)+CFE2Q(I,NY,NX)+CCAHQ(I,NY,NX)+CMGHQ(I,NY,NX)
     9+CF1PQ(I,NY,NX)+CC1PQ(I,NY,NX)+CM1PQ(I,NY,NX))
     4+4.0*(CAL3Q(I,NY,NX)+CFE3Q(I,NY,NX)
     2+CH3PQ(I,NY,NX)+CF2PQ(I,NY,NX)+CC2PQ(I,NY,NX))
     3+5.0*(CAL4Q(I,NY,NX)+CFE4Q(I,NY,NX)))*XNFH
      TIONOU=TIONOU+SBU
      ENDIF
C
C     GAS EXCHANGE FROM LITTER SURFACE VOLATILIZATION-DISSOLUTION
C
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in litter
C     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2
C        in litter
C     X*FLS=convective + diffusive solute flux from trnsfr.f
C     X*DFR=gas exchange between atmosphere and litter from trnsfr.f
C     X*DFG=litter surface gas volatilization from trnsfr.f
C     R*O=net gas transformation in surface litter from nitro.f
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,H1P=HPO4,H2P=H2PO4 in non-band
C        :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,H1B=HPO4,H2B=H2PO4 in band
C     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
C             :*ZN3*=NH3,*H2G*=H2
C
      DO 9680 K=0,2
      OQC(K,0,NY,NX)=OQC(K,0,NY,NX)+XOCFLS(K,3,0,NY,NX)
      OQN(K,0,NY,NX)=OQN(K,0,NY,NX)+XONFLS(K,3,0,NY,NX)
      OQP(K,0,NY,NX)=OQP(K,0,NY,NX)+XOPFLS(K,3,0,NY,NX)
      OQA(K,0,NY,NX)=OQA(K,0,NY,NX)+XOAFLS(K,3,0,NY,NX)
C     IF(NX.EQ.1.AND.NY.EQ.6)THEN
C     WRITE(*,2626)'OQCR',I,J,NFZ,NX,NY,K,OQC(K,0,NY,NX)
C    2,XOCFLS(K,3,0,NY,NX),OQN(K,0,NY,NX),TONQRS(K,NY,NX)
2626  FORMAT(A8,6I4,20E12.4)
C     ENDIF
9680  CONTINUE
      CO2G(0,NY,NX)=CO2G(0,NY,NX)+XCOFLG(3,0,NY,NX)
     2-XCODFG(0,NY,NX)+RCGOX(0,NY,NX)+RC4OX(0,NY,NX)
      CH4G(0,NY,NX)=CH4G(0,NY,NX)+XCHFLG(3,0,NY,NX)
     2-XCHDFG(0,NY,NX)+RCHOX(0,NY,NX)-RC4OX(0,NY,NX)
      OXYG(0,NY,NX)=OXYG(0,NY,NX)+XOXFLG(3,0,NY,NX)
     2-XOXDFG(0,NY,NX)-ROGOX(0,NY,NX)-RC4OX(0,NY,NX)*2.667
C     WRITE(*,1189)'CO2G0',I,J,NFZ,NX,NY
C    2,CO2G(0,NY,NX),XCOFLG(3,0,NY,NX)
C    2,XCODFG(0,NY,NX),RCGOX(0,NY,NX),RC4OX(0,NY,NX)
1189  FORMAT(A8,5I4,30F16.8)
      Z2GG(0,NY,NX)=Z2GG(0,NY,NX)+XNGFLG(3,0,NY,NX)
     2-XNGDFG(0,NY,NX)
      Z2OG(0,NY,NX)=Z2OG(0,NY,NX)+XN2FLG(3,0,NY,NX)
     2-XN2DFG(0,NY,NX)
      ZNH3G(0,NY,NX)=ZNH3G(0,NY,NX)+XN3FLG(3,0,NY,NX)
     2-XN3DFG(0,NY,NX)
      H2GG(0,NY,NX)=H2GG(0,NY,NX)+XHGFLG(3,0,NY,NX)
     2-XHGDFG(0,NY,NX)
      CO2S(0,NY,NX)=CO2S(0,NY,NX)+XCODFR(NY,NX)+XCOFLS(3,0,NY,NX)
     2+XCODFG(0,NY,NX)-RCO2O(0,NY,NX)
      CH4S(0,NY,NX)=CH4S(0,NY,NX)+XCHDFR(NY,NX)+XCHFLS(3,0,NY,NX)
     2+XCHDFG(0,NY,NX)-RCH4O(0,NY,NX)
      OXYS(0,NY,NX)=OXYS(0,NY,NX)+XOXDFR(NY,NX)+XOXFLS(3,0,NY,NX)
     2+XOXDFG(0,NY,NX)-RUPOXO(0,NY,NX)
      Z2GS(0,NY,NX)=Z2GS(0,NY,NX)+XNGDFR(NY,NX)+XNGFLS(3,0,NY,NX)
     2+XNGDFG(0,NY,NX)-RN2G(0,NY,NX)-XN2GS(0,NY,NX)
      Z2OS(0,NY,NX)=Z2OS(0,NY,NX)+XN2DFR(NY,NX)+XN2FLS(3,0,NY,NX)
     2+XN2DFG(0,NY,NX)-RN2O(0,NY,NX)
      H2GS(0,NY,NX)=H2GS(0,NY,NX)+XHGDFR(NY,NX)+XHGFLS(3,0,NY,NX)
     2+XHGDFG(0,NY,NX)-RH2GO(0,NY,NX)
      ZNH4S(0,NY,NX)=ZNH4S(0,NY,NX)+XN4FLW(3,0,NY,NX)
     2+XNH4S(0,NY,NX)+TRN4S(0,NY,NX)
      ZNH3S(0,NY,NX)=ZNH3S(0,NY,NX)+XN3DFR(NY,NX)+XN3FLW(3,0,NY,NX)
     2+XN3DFG(0,NY,NX)+TRN3S(0,NY,NX)
      ZNO3S(0,NY,NX)=ZNO3S(0,NY,NX)+XNOFLW(3,0,NY,NX)
     2+XNO3S(0,NY,NX)+TRNO3(0,NY,NX)
      ZNO2S(0,NY,NX)=ZNO2S(0,NY,NX)+XNXFLS(3,0,NY,NX)
     2+XNO2S(0,NY,NX)
      H1PO4(0,NY,NX)=H1PO4(0,NY,NX)+TRH1P(0,NY,NX)+XH1PFS(3,0,NY,NX)
     2+XH1PS(0,NY,NX)
      H2PO4(0,NY,NX)=H2PO4(0,NY,NX)+TRH2P(0,NY,NX)+XH2PFS(3,0,NY,NX)
     2+XH2PS(0,NY,NX)
C
C     GAS EXCHANGE FROM VOLATILIZTION-DISSOLUTION AT SOIL SURFACE,
C     SURFACE LITTER
C
C     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS
C        =aqueous CO2,CH4,O2,N2,N2O,H2 in soil surface
C     X*DFS= gas exchange between atmosphere and soil surface
C        from trnsfr.f
C     R*O=net gas transformation in soil surface from nitro.f
C     R*L=total gas flux for use in next hour
C     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
C             :*ZN3*=NH3,*H2G*=H2
C
      CO2S(NU(NY,NX),NY,NX)=CO2S(NU(NY,NX),NY,NX)+XCODFS(NY,NX)
      CH4S(NU(NY,NX),NY,NX)=CH4S(NU(NY,NX),NY,NX)+XCHDFS(NY,NX)
      OXYS(NU(NY,NX),NY,NX)=OXYS(NU(NY,NX),NY,NX)+XOXDFS(NY,NX)
      Z2GS(NU(NY,NX),NY,NX)=Z2GS(NU(NY,NX),NY,NX)+XNGDFS(NY,NX)
      Z2OS(NU(NY,NX),NY,NX)=Z2OS(NU(NY,NX),NY,NX)+XN2DFS(NY,NX)
      ZNH3S(NU(NY,NX),NY,NX)=ZNH3S(NU(NY,NX),NY,NX)+XN3DFS(NY,NX)
      ZNH3B(NU(NY,NX),NY,NX)=ZNH3B(NU(NY,NX),NY,NX)+XNBDFS(NY,NX)
      H2GS(NU(NY,NX),NY,NX)=H2GS(NU(NY,NX),NY,NX)+XHGDFS(NY,NX)
      THRE(NY,NX)=THRE(NY,NX)+RCO2O(0,NY,NX)+RCH4O(0,NY,NX)
      UN2GG(NY,NX)=UN2GG(NY,NX)+RN2G(0,NY,NX)
      HN2GG(NY,NX)=HN2GG(NY,NX)+RN2G(0,NY,NX)
C     ROXYF(0,NY,NX)=XOXDFG(0,NY,NX)
C     RCO2F(0,NY,NX)=XCODFG(0,NY,NX)
C     RCH4F(0,NY,NX)=XCHDFG(0,NY,NX)
      ROXYL(0,NY,NX)=XOXDFR(NY,NX)+XOXFLS(3,0,NY,NX)+XOXDFG(0,NY,NX)
     2-(FLQRQ(NY,NX)*COXR(NY,NX)+FLQRI(NY,NX)*COXQ(NY,NX))
      RCH4L(0,NY,NX)=XCHDFR(NY,NX)+XCHFLS(3,0,NY,NX)+XCHDFG(0,NY,NX)
     2-(FLQRQ(NY,NX)*CCHR(NY,NX)+FLQRI(NY,NX)*CCHQ(NY,NX))
      ROXYL(NU(NY,NX),NY,NX)=ROXYL(NU(NY,NX),NY,NX)+XOXDFS(NY,NX)
      RCH4L(NU(NY,NX),NY,NX)=RCH4L(NU(NY,NX),NY,NX)+XCHDFS(NY,NX)
C     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,6644)'OXYS0',I,J,NFZ,NX,NY,NU(NY,NX)
C    2,OXYS(0,NY,NX),XOXDFR(NY,NX),XOXFLS(3,0,NY,NX)
C    3,OXYG(0,NY,NX),XOXFLG(3,0,NY,NX),XOXDFG(0,NY,NX)
C    2,ROGOX(0,NY,NX),RCGSK(0,NY,NX)*2.667,COXYG(0,NY,NX)
C    3,COXYE(NY,NX)
C    3,RUPOXO(0,NY,NX),ROXYL(0,NY,NX)
C    3,FLQRQ(NY,NX)*COXR(NY,NX)+FLQRI(NY,NX)*COXQ(NY,NX)
C     WRITE(*,6644)'NH3S0',I,J,NFZ,NX,NY,NU(NY,NX)
C    2,ZNH3S(0,NY,NX),XN3DFR(NY,NX),XN3FLW(3,0,NY,NX)
C    2,XN3DFG(0,NY,NX),TRN3S(0,NY,NX)
C     WRITE(*,6644)'CO2',I,J,NFZ,NX,NY,NU(NY,NX)
C    2,CO2S(NU(NY,NX),NY,NX),HCO2G(NY,NX)
C    2,CI,XCODFS(NY,NX),XCOFLG(3,NU(NY,NX),NY,NX),TCO2Z(NY,NX)
C    3,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCOR(NY,NX)
C    4,(FLQGI(NY,NX)+FLQRI(NY,NX))*CCOQ(NY,NX)
C    5,CO2S(0,NY,NX),XCODFR(NY,NX),XCOFLS(3,0,NY,NX)
C    2,XCODFG(0,NY,NX),RCO2O(0,NY,NX)
C    5,XCODFG(0,NY,NX),XCODFR(NY,NX),VOLP(0,NY,NX)
C    6,VOLP(NU(NY,NX),NY,NX)
C     WRITE(*,6644)'OXYS1',I,J,NFZ,NX,NY,NU(NY,NX)
C    2,OXYS(NU(NY,NX),NY,NX),HOXYG(NY,NX)
C    2,XOXDFS(NY,NX),XOXFLG(3,NU(NY,NX),NY,NX),TOXYZ(NY,NX)
C    5,XOXDFG(NU(NY,NX),NY,NX),VOLP(NU(NY,NX),NY,NX)
C    3,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*CCOR(NY,NX)
C    4,(FLQGI(NY,NX)+FLQRI(NY,NX))*CCOQ(NY,NX)
C     WRITE(*,6644)'CH4',I,J,NFZ,NX,NY,NU(NY,NX),CH,XCHDFS(NY,NX)
C    2,XCHFLG(3,NU(NY,NX),NY,NX),TCH4Z(NY,NX),FLQGQ(NY,NX)
C    3,FLQRQ(NY,NX),FLQGI(NY,NX),FLQRI(NY,NX),CCHR(NY,NX),CCHQ(NY,NX)
C    4,XCHDFG(0,NY,NX),XCHDFR(NY,NX),CH4S(NU(NY,NX),NY,NX)
C     WRITE(*,6644)'NH3',I,J,NFZ,NX,NY,NU(NY,NX)
C    2,ZNH3S(0,NY,NX),XN3DFR(NY,NX),XN3FLW(3,0,NY,NX)
C    2,XN3DFG(0,NY,NX),TRN3S(0,NY,NX)
6644  FORMAT(A8,6I4,30E12.4)
C     ENDIF
C
C     SURFACE LITTER ION EXCHANGE AND PRECIPITATION
C
C     TR*=total exchange+precipitation from solute.f
C     XN4,XH1P,XH2P=exchangeable NH4,HPO4,H2PO4
C     XOH0,XOH1,XOH2=adsorption sites R-,R-OH,R-OH2
C     PALPO,PFEPO=precipitated AlPO4,FEPO4
C     PCAPM,PCAPD,PCAPH=precipitated CaH2PO4,CaHPO4,apatite
C
      XN4(0,NY,NX)=XN4(0,NY,NX)+TRXN4(0,NY,NX)
      XOH0(0,NY,NX)=XOH0(0,NY,NX)+TRXH0(0,NY,NX)
      XOH1(0,NY,NX)=XOH1(0,NY,NX)+TRXH1(0,NY,NX)
      XOH2(0,NY,NX)=XOH2(0,NY,NX)+TRXH2(0,NY,NX)
      XH1P(0,NY,NX)=XH1P(0,NY,NX)+TRX1P(0,NY,NX)
      XH2P(0,NY,NX)=XH2P(0,NY,NX)+TRX2P(0,NY,NX)
      PALPO(0,NY,NX)=PALPO(0,NY,NX)+TRALPO(0,NY,NX)
      PFEPO(0,NY,NX)=PFEPO(0,NY,NX)+TRFEPO(0,NY,NX)
      PCAPD(0,NY,NX)=PCAPD(0,NY,NX)+TRCAPD(0,NY,NX)
      PCAPH(0,NY,NX)=PCAPH(0,NY,NX)+TRCAPH(0,NY,NX)
      PCAPM(0,NY,NX)=PCAPM(0,NY,NX)+TRCAPM(0,NY,NX)
C     IF(PCAPD(0,NY,NX).LT.ZEROS(NY,NX))THEN
C     WRITE(*,6639)'PCAPD0',I,J,NFZ,PCAPD(0,NY,NX),TRCAPD(0,NY,NX)
6639  FORMAT(A8,3I4,12E12.4)
C     ENDIF
C
C     SURFACE LITTER OUTPUTS
C
C     IF(I.GE.350)THEN
C     WRITE(*,1119)'CO2S0',I,J,NX,NY,CO2S(0,NY,NX),XCODFS(NY,NX)
C    2,XCODFR(NY,NX),XCOFLS(3,0,NY,NX),XCODFG(0,NY,NX),RCO2O(0,NY,NX)
C    3,VOLT(0,NY,NX),CVRD(NY,NX)
C     WRITE(*,1119)'CH4S0',I,J,NX,NY,CH4S(0,NY,NX),XCHDFS(NY,NX)
C    2,XCHDFR(NY,NX),XCHFLS(3,0,NY,NX),XCHDFG(0,NY,NX),RCH4O(0,NY,NX)
C    3,RCH4L(0,NY,NX)
C     WRITE(*,1119)'OXYS0',I,J,NX,NY,OXYS(0,NY,NX),XOXDFR(NY,NX)
C    2,XOXFLS(3,0,NY,NX),XOXDFG(0,NY,NX),RUPOXO(0,NY,NX)
C    3,ROXYL(0,NY,NX),TOXQRS(NY,NX),COXYS(0,NY,NX)
1119  FORMAT(A8,4I4,12E12.4)
C     ENDIF
C     IF(NX.EQ.5)THEN
C     WRITE(*,5533)'NH30',I,J,NX,NY,ZNH4S(0,NY,NX),XN4FLW(3,0,NY,NX)
C    2,XNH4S(0,NY,NX),XN3FLW(3,0,NY,NX),TRN4S(0,NY,NX),TRXN4(0,NY,NX)
C    3,ZNH3S(0,NY,NX),TRN3S(0,NY,NX),XN3DFG(0,NY,NX),XN3DFR(NY,NX)
C    4,ZNHUFA(0,NY,NX),XNO2S(0,NY,NX),XN4(0,NY,NX)*14.0
C     WRITE(*,5533)'ZNO3S0',I,J,NX,NY,ZNO3S(0,NY,NX),XNOFLW(3,0,NY,NX)
C    2,XNO3S(0,NY,NX),TRNO3(0,NY,NX),ZNO2S(0,NY,NX),XNXFLS(3,0,NY,NX)
C    3,XNO2S(0,NY,NX)
C     WRITE(*,5533)'H2PO40',I,J,NX,NY,H2PO4(0,NY,NX)
C    2,XH2PFS(3,0,NY,NX),XH2PS(0,NY,NX),TRH2P(0,NY,NX)
5533  FORMAT(A8,4I4,20F12.4)
C     ENDIF
C     WRITE(*,5544)'HP140',I,J,NX,NY,H1PO4(0,NY,NX)
C    2,XH1P(0,NY,NX),TRH1P(0,NY,NX),XH1PFS(3,0,NY,NX)
C    2,XH1PS(0,NY,NX),TP1QRS(NY,NX)
C     WRITE(*,5544)'HP240',I,J,NX,NY,H2PO4(0,NY,NX)
C    2,XH2P(0,NY,NX),TRH2P(0,NY,NX),XH2PFS(3,0,NY,NX)
C    2,XH2PS(0,NY,NX),TPOQRS(NY,NX)
5544  FORMAT(A8,4I4,40E12.4)
C
C     OVERLAND FLOW
C
C     TQR=net water runoff from watsub.f
C
      IF(ABS(TQR(NY,NX)).GT.ZEROS(NY,NX))THEN
C     DOC, DON, DOP
C
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate
C     T*QRS=net overland solute flux from runoff
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,H1P=HPO4,H2P=H2PO4 in non-band
C
      DO 8570 K=0,2
      OQC(K,0,NY,NX)=OQC(K,0,NY,NX)+TOCQRS(K,NY,NX)
      OQN(K,0,NY,NX)=OQN(K,0,NY,NX)+TONQRS(K,NY,NX)
      OQP(K,0,NY,NX)=OQP(K,0,NY,NX)+TOPQRS(K,NY,NX)
      OQA(K,0,NY,NX)=OQA(K,0,NY,NX)+TOAQRS(K,NY,NX)
C     IF(NX.EQ.1.AND.NY.EQ.6)THEN
C     WRITE(*,2628)'OQCS',I,J,NFZ,NX,NY,K,OQC(K,0,NY,NX)
C    2,TOCQRS(K,NY,NX),OQN(K,0,NY,NX),TONQRS(K,NY,NX)
2628  FORMAT(A8,6I4,20E12.4)
C     ENDIF
8570  CONTINUE
C
C     LITTER SOLUTES
C
      CO2S(0,NY,NX)=CO2S(0,NY,NX)+TCOQRS(NY,NX)
      CH4S(0,NY,NX)=CH4S(0,NY,NX)+TCHQRS(NY,NX)
      OXYS(0,NY,NX)=OXYS(0,NY,NX)+TOXQRS(NY,NX)
      Z2GS(0,NY,NX)=Z2GS(0,NY,NX)+TNGQRS(NY,NX)
      Z2OS(0,NY,NX)=Z2OS(0,NY,NX)+TN2QRS(NY,NX)
      H2GS(0,NY,NX)=H2GS(0,NY,NX)+THGQRS(NY,NX)
      ZNH4S(0,NY,NX)=ZNH4S(0,NY,NX)+TN4QRS(NY,NX)
      ZNH3S(0,NY,NX)=ZNH3S(0,NY,NX)+TN3QRS(NY,NX)
      ZNO3S(0,NY,NX)=ZNO3S(0,NY,NX)+TNOQRS(NY,NX)
      ZNO2S(0,NY,NX)=ZNO2S(0,NY,NX)+TNXQRS(NY,NX)
      H1PO4(0,NY,NX)=H1PO4(0,NY,NX)+TP1QRS(NY,NX)
      H2PO4(0,NY,NX)=H2PO4(0,NY,NX)+TPOQRS(NY,NX)
C
C     TQR*=net overland solute flux in runoff
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C        :*1=non-band,*B=band
C
      IF(ISALTG.NE.0)THEN
      ZAL(0,NY,NX)=ZAL(0,NY,NX)+TQRAL(NY,NX)
      ZFE(0,NY,NX)=ZFE(0,NY,NX)+TQRFE(NY,NX)
      ZHY(0,NY,NX)=ZHY(0,NY,NX)+TQRHY(NY,NX)
      ZCA(0,NY,NX)=ZCA(0,NY,NX)+TQRCA(NY,NX)
      ZMG(0,NY,NX)=ZMG(0,NY,NX)+TQRMG(NY,NX)
      ZNA(0,NY,NX)=ZNA(0,NY,NX)+TQRNA(NY,NX)
      ZKA(0,NY,NX)=ZKA(0,NY,NX)+TQRKA(NY,NX)
      ZOH(0,NY,NX)=ZOH(0,NY,NX)+TQROH(NY,NX)
      ZSO4(0,NY,NX)=ZSO4(0,NY,NX)+TQRSO(NY,NX)
      ZCL(0,NY,NX)=ZCL(0,NY,NX)+TQRCL(NY,NX)
      ZCO3(0,NY,NX)=ZCO3(0,NY,NX)+TQRC3(NY,NX)
      ZHCO3(0,NY,NX)=ZHCO3(0,NY,NX)+TQRHC(NY,NX)
      ZALOH1(0,NY,NX)=ZALOH1(0,NY,NX)+TQRAL1(NY,NX)
      ZALOH2(0,NY,NX)=ZALOH2(0,NY,NX)+TQRAL2(NY,NX)
      ZALOH3(0,NY,NX)=ZALOH3(0,NY,NX)+TQRAL3(NY,NX)
      ZALOH4(0,NY,NX)=ZALOH4(0,NY,NX)+TQRAL4(NY,NX)
      ZALS(0,NY,NX)=ZALS(0,NY,NX)+TQRALS(NY,NX)
      ZFEOH1(0,NY,NX)=ZFEOH1(0,NY,NX)+TQRFE1(NY,NX)
      ZFEOH2(0,NY,NX)=ZFEOH2(0,NY,NX)+TQRFE2(NY,NX)
      ZFEOH3(0,NY,NX)=ZFEOH3(0,NY,NX)+TQRFE3(NY,NX)
      ZFEOH4(0,NY,NX)=ZFEOH4(0,NY,NX)+TQRFE4(NY,NX)
      ZFES(0,NY,NX)=ZFES(0,NY,NX)+TQRFES(NY,NX)
      ZCAO(0,NY,NX)=ZCAO(0,NY,NX)+TQRCAO(NY,NX)
      ZCAC(0,NY,NX)=ZCAC(0,NY,NX)+TQRCAC(NY,NX)
      ZCAH(0,NY,NX)=ZCAH(0,NY,NX)+TQRCAH(NY,NX)
      ZCAS(0,NY,NX)=ZCAS(0,NY,NX)+TQRCAS(NY,NX)
      ZMGO(0,NY,NX)=ZMGO(0,NY,NX)+TQRMGO(NY,NX)
      ZMGC(0,NY,NX)=ZMGC(0,NY,NX)+TQRMGC(NY,NX)
      ZMGH(0,NY,NX)=ZMGH(0,NY,NX)+TQRMGH(NY,NX)
      ZMGS(0,NY,NX)=ZMGS(0,NY,NX)+TQRMGS(NY,NX)
      ZNAC(0,NY,NX)=ZNAC(0,NY,NX)+TQRNAC(NY,NX)
      ZNAS(0,NY,NX)=ZNAS(0,NY,NX)+TQRNAS(NY,NX)
      ZKAS(0,NY,NX)=ZKAS(0,NY,NX)+TQRKAS(NY,NX)
      H0PO4(0,NY,NX)=H0PO4(0,NY,NX)+TQRH0P(NY,NX)
      H3PO4(0,NY,NX)=H3PO4(0,NY,NX)+TQRH3P(NY,NX)
      ZFE1P(0,NY,NX)=ZFE1P(0,NY,NX)+TQRF1P(NY,NX)
      ZFE2P(0,NY,NX)=ZFE2P(0,NY,NX)+TQRF2P(NY,NX)
      ZCA0P(0,NY,NX)=ZCA0P(0,NY,NX)+TQRC0P(NY,NX)
      ZCA1P(0,NY,NX)=ZCA1P(0,NY,NX)+TQRC1P(NY,NX)
      ZCA2P(0,NY,NX)=ZCA2P(0,NY,NX)+TQRC2P(NY,NX)
      ZMG1P(0,NY,NX)=ZMG1P(0,NY,NX)+TQRM1P(NY,NX)
      ENDIF
      ENDIF
C
C     INTERNAL SURFACE SEDIMENT TRANSPORT
C
C     IERSNG=erosion flag from site file
C     T*ER=net sediment flux
C     sediment code:SED=total sediment,SAN=sand,SIL=silt,CLA=clay
C       :OMC,OMN,OMP=microbial C,N,P; ORC=microbial residue C,N,P
C       :OHC,OHN,OHP=adsorbed C,N,P; OSC,OSN,OSP=humus C,N,P
C       :NH4,NH3,NHU,NO3=fertilizer NH4,NH3,urea,NO3 in non-band
C       :NH4B,NH3B,NHUB,NO3B=fertilizer NH4,NH3,urea,NO3 in band
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
C        =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
C       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C       :PALO,PFEO=precip AlOH,FeOH
C       :PCAC,PCAS=precip CaCO3,CaSO4
C       :PALP,PFEP=precip AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precip AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite in non-band
C       :PCPMB,PCPDB,PCPHB=precip CaH2PO4,CaHPO4,apatite in band
C
      IF((IERSNG.EQ.1.OR.IERSNG.EQ.3)
     2.AND.ABS(TSEDER(NY,NX)).GT.ZEROS(NY,NX))THEN
      TSED(NY,NX)=TSED(NY,NX)+TSEDER(NY,NX)
C
C     SOIL MINERAL FRACTIONS
C
C     sediment code:XSED=total,XSAN=sand,XSIL=silt,XCLA=clay
C                  :CEC=cation exchange capacity
C                  :AEC=anion exchange capacity
C
      SAND(NU(NY,NX),NY,NX)=SAND(NU(NY,NX),NY,NX)+TSANER(NY,NX)
      SILT(NU(NY,NX),NY,NX)=SILT(NU(NY,NX),NY,NX)+TSILER(NY,NX)
      CLAY(NU(NY,NX),NY,NX)=CLAY(NU(NY,NX),NY,NX)+TCLAER(NY,NX)
      XCEC(NU(NY,NX),NY,NX)=XCEC(NU(NY,NX),NY,NX)+TCECER(NY,NX)
      XAEC(NU(NY,NX),NY,NX)=XAEC(NU(NY,NX),NY,NX)+TAECER(NY,NX)
C
C     FERTILIZER POOLS
C
C     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3 in non-band
C                  :NH4B,NH3B,NHUB,NO3B=NH4,NH3,urea,NO3 in band
C
      ZNH4FA(NU(NY,NX),NY,NX)=ZNH4FA(NU(NY,NX),NY,NX)+TNH4ER(NY,NX)
      ZNH3FA(NU(NY,NX),NY,NX)=ZNH3FA(NU(NY,NX),NY,NX)+TNH3ER(NY,NX)
      ZNHUFA(NU(NY,NX),NY,NX)=ZNHUFA(NU(NY,NX),NY,NX)+TNHUER(NY,NX)
      ZNO3FA(NU(NY,NX),NY,NX)=ZNO3FA(NU(NY,NX),NY,NX)+TNO3ER(NY,NX)
      ZNH4FB(NU(NY,NX),NY,NX)=ZNH4FB(NU(NY,NX),NY,NX)+TNH4EB(NY,NX)
      ZNH3FB(NU(NY,NX),NY,NX)=ZNH3FB(NU(NY,NX),NY,NX)+TNH3EB(NY,NX)
      ZNHUFB(NU(NY,NX),NY,NX)=ZNHUFB(NU(NY,NX),NY,NX)+TNHUEB(NY,NX)
      ZNO3FB(NU(NY,NX),NY,NX)=ZNO3FB(NU(NY,NX),NY,NX)+TNO3EB(NY,NX)
C
C     EXCHANGEABLE CATIONS AND ANIONS
C
C     sediment code
C       :XN4,XNB=adsorbed NH4 in non-band,band
C       :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
C           =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
C       :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2 in non-band
C       :XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C       :XH1P,XH2P=adsorbed HPO4,H2PO4 in non-band
C       :XH1PB,XP2PB=adsorbed HPO4,H2PO4 in band
C
      XN4(NU(NY,NX),NY,NX)=XN4(NU(NY,NX),NY,NX)+TN4ER(NY,NX)
      XNB(NU(NY,NX),NY,NX)=XNB(NU(NY,NX),NY,NX)+TNBER(NY,NX)
      XHY(NU(NY,NX),NY,NX)=XHY(NU(NY,NX),NY,NX)+THYER(NY,NX)
      XAL(NU(NY,NX),NY,NX)=XAL(NU(NY,NX),NY,NX)+TALER(NY,NX)
      XFE(NU(NY,NX),NY,NX)=XFE(NU(NY,NX),NY,NX)+TFEER(NY,NX)
      XCA(NU(NY,NX),NY,NX)=XCA(NU(NY,NX),NY,NX)+TCAER(NY,NX)
      XMG(NU(NY,NX),NY,NX)=XMG(NU(NY,NX),NY,NX)+TMGER(NY,NX)
      XNA(NU(NY,NX),NY,NX)=XNA(NU(NY,NX),NY,NX)+TNAER(NY,NX)
      XKA(NU(NY,NX),NY,NX)=XKA(NU(NY,NX),NY,NX)+TKAER(NY,NX)
      XHC(NU(NY,NX),NY,NX)=XHC(NU(NY,NX),NY,NX)+THCER(NY,NX)
      XALO2(NU(NY,NX),NY,NX)=XALO2(NU(NY,NX),NY,NX)+TAL2ER(NY,NX)
      XFEO2(NU(NY,NX),NY,NX)=XFEO2(NU(NY,NX),NY,NX)+TFE2ER(NY,NX)
      XOH0(NU(NY,NX),NY,NX)=XOH0(NU(NY,NX),NY,NX)+TOH0ER(NY,NX)
      XOH1(NU(NY,NX),NY,NX)=XOH1(NU(NY,NX),NY,NX)+TOH1ER(NY,NX)
      XOH2(NU(NY,NX),NY,NX)=XOH2(NU(NY,NX),NY,NX)+TOH2ER(NY,NX)
      XH1P(NU(NY,NX),NY,NX)=XH1P(NU(NY,NX),NY,NX)+TH1PER(NY,NX)
      XH2P(NU(NY,NX),NY,NX)=XH2P(NU(NY,NX),NY,NX)+TH2PER(NY,NX)
      XOH0B(NU(NY,NX),NY,NX)=XOH0B(NU(NY,NX),NY,NX)+TOH0EB(NY,NX)
      XOH1B(NU(NY,NX),NY,NX)=XOH1B(NU(NY,NX),NY,NX)+TOH1EB(NY,NX)
      XOH2B(NU(NY,NX),NY,NX)=XOH2B(NU(NY,NX),NY,NX)+TOH2EB(NY,NX)
      XH1PB(NU(NY,NX),NY,NX)=XH1PB(NU(NY,NX),NY,NX)+TH1PEB(NY,NX)
      XH2PB(NU(NY,NX),NY,NX)=XH2PB(NU(NY,NX),NY,NX)+TH2PEB(NY,NX)
C
C     PRECIPITATES
C
C     sediment code
C       :PALO,PFEO=precip AlOH,FeOH
C       :PCAC,PCAS=precip CaCO3,CaSO4
C       :PALP,PFEP=precip AlPO4,FEPO4 in non-band
C       :PALPB,PFEPB=precip AlPO4,FEPO4 in band
C       :PCPM,PCPD,PCPH=precip CaH2PO4,CaHPO4,apatite in non-band
C       :PCPMB,PCPDB,PCPHB=precip CaH2PO4,CaHPO4,apatite in band
C
      PALOH(NU(NY,NX),NY,NX)=PALOH(NU(NY,NX),NY,NX)+TALOER(NY,NX)
      PFEOH(NU(NY,NX),NY,NX)=PFEOH(NU(NY,NX),NY,NX)+TFEOER(NY,NX)
      PCACO(NU(NY,NX),NY,NX)=PCACO(NU(NY,NX),NY,NX)+TCACER(NY,NX)
      PCASO(NU(NY,NX),NY,NX)=PCASO(NU(NY,NX),NY,NX)+TCASER(NY,NX)
      PALPO(NU(NY,NX),NY,NX)=PALPO(NU(NY,NX),NY,NX)+TALPER(NY,NX)
      PFEPO(NU(NY,NX),NY,NX)=PFEPO(NU(NY,NX),NY,NX)+TFEPER(NY,NX)
      PCAPD(NU(NY,NX),NY,NX)=PCAPD(NU(NY,NX),NY,NX)+TCPDER(NY,NX)
      PCAPH(NU(NY,NX),NY,NX)=PCAPH(NU(NY,NX),NY,NX)+TCPHER(NY,NX)
      PCAPM(NU(NY,NX),NY,NX)=PCAPM(NU(NY,NX),NY,NX)+TCPMER(NY,NX)
      PALPB(NU(NY,NX),NY,NX)=PALPB(NU(NY,NX),NY,NX)+TALPEB(NY,NX)
      PFEPB(NU(NY,NX),NY,NX)=PFEPB(NU(NY,NX),NY,NX)+TFEPEB(NY,NX)
      PCPDB(NU(NY,NX),NY,NX)=PCPDB(NU(NY,NX),NY,NX)+TCPDEB(NY,NX)
      PCPHB(NU(NY,NX),NY,NX)=PCPHB(NU(NY,NX),NY,NX)+TCPHEB(NY,NX)
      PCPMB(NU(NY,NX),NY,NX)=PCPMB(NU(NY,NX),NY,NX)+TCPMEB(NY,NX)
C
C     ORGANIC MATTER
C
C     sediment code
C        :OMC,OMN,OMP=microbial C,N,P
C        :ORC,ORN,ORP=microbial residue C,N,P
C        :OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
C        :OSC,OSA,OSN,OSP=SOC,colonized SOC,SON,SOP
C           (K=0:woody litter, K=1:non-woody litter,
C            K=2:manure, K=3:POC, K=4:humus)
C
      DORGP=0.0
      DO 9280 K=0,5
      DO 9280 NO=1,7
      DO 9280 M=1,3
      OMC(M,NO,K,NU(NY,NX),NY,NX)=OMC(M,NO,K,NU(NY,NX),NY,NX)
     2+TOMCER(M,NO,K,NY,NX)
      OMN(M,NO,K,NU(NY,NX),NY,NX)=OMN(M,NO,K,NU(NY,NX),NY,NX)
     2+TOMNER(M,NO,K,NY,NX)
      OMP(M,NO,K,NU(NY,NX),NY,NX)=OMP(M,NO,K,NU(NY,NX),NY,NX)
     2+TOMPER(M,NO,K,NY,NX)
      DORGE(NY,NX)=DORGE(NY,NX)+TOMCER(M,NO,K,NY,NX)
      DORGP=DORGP+TOMPER(M,NO,K,NY,NX)
9280  CONTINUE
      DO 9275 K=0,4
      DO 9270 M=1,2
      ORC(M,K,NU(NY,NX),NY,NX)=ORC(M,K,NU(NY,NX),NY,NX)
     2+TORCER(M,K,NY,NX)
      ORN(M,K,NU(NY,NX),NY,NX)=ORN(M,K,NU(NY,NX),NY,NX)
     2+TORNER(M,K,NY,NX)
      ORP(M,K,NU(NY,NX),NY,NX)=ORP(M,K,NU(NY,NX),NY,NX)
     2+TORPER(M,K,NY,NX)
      DORGE(NY,NX)=DORGE(NY,NX)+TORCER(M,K,NY,NX)
      DORGP=DORGP+TORPER(M,K,NY,NX)
9270  CONTINUE
      OHC(K,NU(NY,NX),NY,NX)=OHC(K,NU(NY,NX),NY,NX)+TOHCER(K,NY,NX)
      OHN(K,NU(NY,NX),NY,NX)=OHN(K,NU(NY,NX),NY,NX)+TOHNER(K,NY,NX)
      OHP(K,NU(NY,NX),NY,NX)=OHP(K,NU(NY,NX),NY,NX)+TOHPER(K,NY,NX)
      OHA(K,NU(NY,NX),NY,NX)=OHA(K,NU(NY,NX),NY,NX)+TOHAER(K,NY,NX)
      DORGE(NY,NX)=DORGE(NY,NX)+TOHCER(K,NY,NX)+TOHAER(K,NY,NX)
      DORGP=DORGP+TOHPER(K,NY,NX)
      DO 9265 M=1,5
      OSC(M,K,NU(NY,NX),NY,NX)=OSC(M,K,NU(NY,NX),NY,NX)
     2+TOSCER(M,K,NY,NX)
      OSA(M,K,NU(NY,NX),NY,NX)=OSA(M,K,NU(NY,NX),NY,NX)
     2+TOSAER(M,K,NY,NX)
      OSN(M,K,NU(NY,NX),NY,NX)=OSN(M,K,NU(NY,NX),NY,NX)
     2+TOSNER(M,K,NY,NX)
      OSP(M,K,NU(NY,NX),NY,NX)=OSP(M,K,NU(NY,NX),NY,NX)
     2+TOSPER(M,K,NY,NX)
      DORGE(NY,NX)=DORGE(NY,NX)+TOSCER(M,K,NY,NX)
      DORGP=DORGP+TOSPER(M,K,NY,NX)
9265  CONTINUE
9275  CONTINUE
C     WRITE(*,6637)'DORGP',I,J,NY,NX,DORGP,TSEDER(NY,NX)
6637  FORMAT(A8,4I4,12E12.4)
      ENDIF
C
C     OVERLAND SNOW REDISTRIBUTION
C
C     *W=solute content of snowpack
C     T*QSS=net overland solute flux from snowpack
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C             :OC=DOC,ON=DON,OP=DOP,OA=acetate
C             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4
C
      IF(TQS(NY,NX).NE.0.0)THEN
      CO2W(1,NY,NX)=CO2W(1,NY,NX)+TCOQSS(NY,NX)
      CH4W(1,NY,NX)=CH4W(1,NY,NX)+TCHQSS(NY,NX)
      OXYW(1,NY,NX)=OXYW(1,NY,NX)+TOXQSS(NY,NX)
      ZNGW(1,NY,NX)=ZNGW(1,NY,NX)+TNGQSS(NY,NX)
      ZN2W(1,NY,NX)=ZN2W(1,NY,NX)+TN2QSS(NY,NX)
      ZN4W(1,NY,NX)=ZN4W(1,NY,NX)+TN4QSS(NY,NX)
      ZN3W(1,NY,NX)=ZN3W(1,NY,NX)+TN3QSS(NY,NX)
      ZNOW(1,NY,NX)=ZNOW(1,NY,NX)+TNOQSS(NY,NX)
      Z1PW(1,NY,NX)=Z1PW(1,NY,NX)+TP1QSS(NY,NX)
      ZHPW(1,NY,NX)=ZHPW(1,NY,NX)+TPOQSS(NY,NX)
C
C     NET SALT FLUXES FROM SNOWPACK
C
C     TQS*=net overland solute flux in snow drift
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C        :*1=non-band,*B=band
C
      IF(ISALTG.NE.0)THEN
      ZALW(1,NY,NX)=ZALW(1,NY,NX)+TQSAL(NY,NX)
      ZFEW(1,NY,NX)=ZFEW(1,NY,NX)+TQSFE(NY,NX)
      ZHYW(1,NY,NX)=ZHYW(1,NY,NX)+TQSHY(NY,NX)
      ZCAW(1,NY,NX)=ZCAW(1,NY,NX)+TQSCA(NY,NX)
      ZMGW(1,NY,NX)=ZMGW(1,NY,NX)+TQSMG(NY,NX)
      ZNAW(1,NY,NX)=ZNAW(1,NY,NX)+TQSNA(NY,NX)
      ZKAW(1,NY,NX)=ZKAW(1,NY,NX)+TQSKA(NY,NX)
      ZOHW(1,NY,NX)=ZOHW(1,NY,NX)+TQSOH(NY,NX)
      ZSO4W(1,NY,NX)=ZSO4W(1,NY,NX)+TQSSO(NY,NX)
      ZCLW(1,NY,NX)=ZCLW(1,NY,NX)+TQSCL(NY,NX)
      ZCO3W(1,NY,NX)=ZCO3W(1,NY,NX)+TQSC3(NY,NX)
      ZHCO3W(1,NY,NX)=ZHCO3W(1,NY,NX)+TQSHC(NY,NX)
      ZALH1W(1,NY,NX)=ZALH1W(1,NY,NX)+TQSAL1(NY,NX)
      ZALH2W(1,NY,NX)=ZALH2W(1,NY,NX)+TQSAL2(NY,NX)
      ZALH3W(1,NY,NX)=ZALH3W(1,NY,NX)+TQSAL3(NY,NX)
      ZALH4W(1,NY,NX)=ZALH4W(1,NY,NX)+TQSAL4(NY,NX)
      ZALSW(1,NY,NX)=ZALSW(1,NY,NX)+TQSALS(NY,NX)
      ZFEH1W(1,NY,NX)=ZFEH1W(1,NY,NX)+TQSFE1(NY,NX)
      ZFEH2W(1,NY,NX)=ZFEH2W(1,NY,NX)+TQSFE2(NY,NX)
      ZFEH3W(1,NY,NX)=ZFEH3W(1,NY,NX)+TQSFE3(NY,NX)
      ZFEH4W(1,NY,NX)=ZFEH4W(1,NY,NX)+TQSFE4(NY,NX)
      ZFESW(1,NY,NX)=ZFESW(1,NY,NX)+TQSFES(NY,NX)
      ZCAOW(1,NY,NX)=ZCAOW(1,NY,NX)+TQSCAO(NY,NX)
      ZCACW(1,NY,NX)=ZCACW(1,NY,NX)+TQSCAC(NY,NX)
      ZCAHW(1,NY,NX)=ZCAHW(1,NY,NX)+TQSCAH(NY,NX)
      ZCASW(1,NY,NX)=ZCASW(1,NY,NX)+TQSCAS(NY,NX)
      ZMGOW(1,NY,NX)=ZMGOW(1,NY,NX)+TQSMGO(NY,NX)
      ZMGCW(1,NY,NX)=ZMGCW(1,NY,NX)+TQSMGC(NY,NX)
      ZMGHW(1,NY,NX)=ZMGHW(1,NY,NX)+TQSMGH(NY,NX)
      ZMGSW(1,NY,NX)=ZMGSW(1,NY,NX)+TQSMGS(NY,NX)
      ZNACW(1,NY,NX)=ZNACW(1,NY,NX)+TQSNAC(NY,NX)
      ZNASW(1,NY,NX)=ZNASW(1,NY,NX)+TQSNAS(NY,NX)
      ZKASW(1,NY,NX)=ZKASW(1,NY,NX)+TQSKAS(NY,NX)
      H0PO4W(1,NY,NX)=H0PO4W(1,NY,NX)+TQSH0P(NY,NX)
      H3PO4W(1,NY,NX)=H3PO4W(1,NY,NX)+TQSH3P(NY,NX)
      ZFE1PW(1,NY,NX)=ZFE1PW(1,NY,NX)+TQSF1P(NY,NX)
      ZFE2PW(1,NY,NX)=ZFE2PW(1,NY,NX)+TQSF2P(NY,NX)
      ZCA0PW(1,NY,NX)=ZCA0PW(1,NY,NX)+TQSC0P(NY,NX)
      ZCA1PW(1,NY,NX)=ZCA1PW(1,NY,NX)+TQSC1P(NY,NX)
      ZCA2PW(1,NY,NX)=ZCA2PW(1,NY,NX)+TQSC2P(NY,NX)
      ZMG1PW(1,NY,NX)=ZMG1PW(1,NY,NX)+TQSM1P(NY,NX)
      ENDIF
      ENDIF
C
C     UPDATE STATE VARIABLES WITH TOTAL FLUXES CALCULATED ABOVE
C
C     IF(J.EQ.24)THEN
C
C     SNOWPACK VARIABLES NEEDED FOR WATER, C, N, P, O, SOLUTE AND
C     ENERGY BALANCES INCLUDING SUM OF ALL CURRENT STATE VARIABLES,
C     CUMULATIVE SUMS OF ALL ADDITIONS AND REMOVALS
C
C     VOLSSL,VOLWSL,VOLVSL,VOLISL=snow water equivalent,water,
C        vapor,ice volume in snowpack layer
C     VOLWSO,UVOLW=total landscape, grid cell water content
C     DENSI=ice density from starts.f
C     HEATSO=total landscape heat content
C     VHCPW=snowpack layer heat capacity
C     TKW=snowpack layer temperature
C     TLCO2G,UCO2S=total landscape, grid cell CO2 content
C     OXYGSO=total landscape O2 content
C     TLN2G,TLNH4,TLNO3,TLPO4=total landscape N2,NH4,NO3,PO4 content
C     TION=total landscape ion content
C     *W=solute content of snowpack
C
      DO 9785 L=1,JS
      WS=VOLSSL(L,NY,NX)+VOLWSL(L,NY,NX)+VOLVSL(L,NY,NX)
     2+VOLISL(L,NY,NX)*DENSI
      VOLWSO=VOLWSO+WS
      UVOLW(NY,NX)=UVOLW(NY,NX)+WS
      ENGYW=VHCPW(L,NY,NX)*TKW(L,NY,NX)
      HEATSO=HEATSO+ENGYW
      TLCO2G=TLCO2G+CO2W(L,NY,NX)+CH4W(L,NY,NX)
      UCO2S(NY,NX)=UCO2S(NY,NX)+CO2W(L,NY,NX)+CH4W(L,NY,NX)
      OXYGSO=OXYGSO+OXYW(L,NY,NX)
      TLN2G=TLN2G+ZNGW(L,NY,NX)+ZN2W(L,NY,NX)
      TLNH4=TLNH4+ZN4W(L,NY,NX)+ZN3W(L,NY,NX)
      TLNO3=TLNO3+ZNOW(L,NY,NX)
      TLPO4=TLPO4+Z1PW(L,NY,NX)+ZHPW(L,NY,NX)
      IF(ISALTG.NE.0)THEN
      SSW=ZALW(L,NY,NX)+ZFEW(L,NY,NX)+ZHYW(L,NY,NX)+ZCAW(L,NY,NX)
     2+ZMGW(L,NY,NX)+ZNAW(L,NY,NX)+ZKAW(L,NY,NX)+ZOHW(L,NY,NX)
     3+ZSO4W(L,NY,NX)+ZCLW(L,NY,NX)+ZCO3W(L,NY,NX)+H0PO4W(L,NY,NX)
     4+2.0*(ZHCO3W(L,NY,NX)+ZALH1W(L,NY,NX)
     5+ZALSW(L,NY,NX)+ZFEH1W(L,NY,NX)+ZFESW(L,NY,NX)+ZCAOW(L,NY,NX)
     6+ZCACW(L,NY,NX)+ZCASW(L,NY,NX)+ZMGOW(L,NY,NX)+ZMGCW(L,NY,NX)
     7+ZMGSW(L,NY,NX)+ZNACW(L,NY,NX)+ZNASW(L,NY,NX)+ZKASW(L,NY,NX)
     8+ZCA0PW(L,NY,NX))
     9+3.0*(ZALH2W(L,NY,NX)+ZFEH2W(L,NY,NX)+ZCAHW(L,NY,NX)
     1+ZMGHW(L,NY,NX)+ZFE1PW(L,NY,NX)+ZCA1PW(L,NY,NX)
     2+ZMG1PW(L,NY,NX))
     2+4.0*(ZALH3W(L,NY,NX)+ZFEH3W(L,NY,NX)+H3PO4W(L,NY,NX)
     3+ZFE2PW(L,NY,NX)+ZCA2PW(L,NY,NX))
     5+5.0*(ZALH4W(L,NY,NX)+ZFEH4W(L,NY,NX))
      TION=TION+SSW
C     WRITE(*,3335)'SSW',I,J,L,SSW
C    2,ZALW(L,NY,NX),ZFEW(L,NY,NX),ZHYW(L,NY,NX),ZCAW(L,NY,NX)
C    2,ZMGW(L,NY,NX),ZNAW(L,NY,NX),ZKAW(L,NY,NX),ZOHW(L,NY,NX)
C    3,ZSO4W(L,NY,NX),ZCLW(L,NY,NX),ZCO3W(L,NY,NX),H0PO4W(L,NY,NX)
C    4,ZHCO3W(L,NY,NX),ZALH1W(L,NY,NX)
C    5,ZALSW(L,NY,NX),ZFEH1W(L,NY,NX),ZFESW(L,NY,NX),ZCAOW(L,NY,NX)
C    6,ZCACW(L,NY,NX),ZCASW(L,NY,NX),ZMGOW(L,NY,NX),ZMGCW(L,NY,NX)
C    7,ZMGSW(L,NY,NX),ZNACW(L,NY,NX),ZNASW(L,NY,NX),ZKASW(L,NY,NX)
C    8,ZCA0PW(L,NY,NX)
C    9,ZALH2W(L,NY,NX),ZFEH2W(L,NY,NX),ZCAHW(L,NY,NX)
C    1,ZMGHW(L,NY,NX),ZFE1PW(L,NY,NX),ZCA1PW(L,NY,NX)
C    2,ZMG1PW(L,NY,NX)
C    2,ZALH3W(L,NY,NX),ZFEH3W(L,NY,NX),H3PO4W(L,NY,NX)
C    3,ZFE2PW(L,NY,NX),ZCA2PW(L,NY,NX)
C    5,ZALH4W(L,NY,NX),ZFEH4W(L,NY,NX)
      ENDIF
9785  CONTINUE
C
C     TOTAL C,N,P, SALTS IN SURFACE LITTER
C
      DC=0.0
      DN=0.0
      DP=0.0
      DCC=0.0
      DNC=0.0
      DPC=0.0
      DO 6975 K=0,5
      RC0(K,NY,NX)=0.0
6975  CONTINUE
      OMCL(0,NY,NX)=0.0
      OMNL(0,NY,NX)=0.0
      DO 6970 K=0,5
      IF(K.NE.3.AND.K.NE.4)THEN
C
C     TOTAL MICROBIAL C,N,P
C
C     RC0=surface litter C content
C     TOMT,TONT,TOPT=total microbial C,N,P
C     OMC,OMN,OMP=microbial C,N,P
C
      DO 6960 N=1,7
      DO 6960 M=1,3
      DC=DC+OMC(M,N,K,0,NY,NX)
      DN=DN+OMN(M,N,K,0,NY,NX)
      DP=DP+OMP(M,N,K,0,NY,NX)
      RC0(K,NY,NX)=RC0(K,NY,NX)+OMC(M,N,K,0,NY,NX)
      TOMT(NY,NX)=TOMT(NY,NX)+OMC(M,N,K,0,NY,NX)
      TONT(NY,NX)=TONT(NY,NX)+OMN(M,N,K,0,NY,NX)
      TOPT(NY,NX)=TOPT(NY,NX)+OMP(M,N,K,0,NY,NX)
      OMCL(0,NY,NX)=OMCL(0,NY,NX)+OMC(M,N,K,0,NY,NX)
      OMNL(0,NY,NX)=OMNL(0,NY,NX)+OMN(M,N,K,0,NY,NX)
6960  CONTINUE
      ENDIF
6970  CONTINUE
C
C     TOTAL MICROBIAL RESIDUE C,N,P
C
C     ORC,ORN,ORP=microbial residue C,N,P
C
      DO 6900 K=0,2
      DO 6940 M=1,2
      DC=DC+ORC(M,K,0,NY,NX)
      DN=DN+ORN(M,K,0,NY,NX)
      DP=DP+ORP(M,K,0,NY,NX)
      RC0(K,NY,NX)=RC0(K,NY,NX)+ORC(M,K,0,NY,NX)
6940  CONTINUE
C
C     TOTAL DOC, DON, DOP
C
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
C     OQCH,OQNH,OQPH,OQAH=DOC,DON,DOP,acetate in macropores
C     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
C
      DC=DC+OQC(K,0,NY,NX)+OQCH(K,0,NY,NX)+OHC(K,0,NY,NX)
     2+OQA(K,0,NY,NX)+OQAH(K,0,NY,NX)+OHA(K,0,NY,NX)
      DN=DN+OQN(K,0,NY,NX)+OQNH(K,0,NY,NX)+OHN(K,0,NY,NX)
      DP=DP+OQP(K,0,NY,NX)+OQPH(K,0,NY,NX)+OHP(K,0,NY,NX)
      RC0(K,NY,NX)=RC0(K,NY,NX)+OQC(K,0,NY,NX)+OQCH(K,0,NY,NX)
     2+OHC(K,0,NY,NX)+OQA(K,0,NY,NX)+OQAH(K,0,NY,NX)+OHA(K,0,NY,NX)
C
C     TOTAL PLANT RESIDUE C,N,P
C
C     OSC,OSN,OSP=SOC,SON,SOP
C
      DO 6930 M=1,5
      IF(M.LE.4)THEN
      DC=DC+OSC(M,K,0,NY,NX)
      DN=DN+OSN(M,K,0,NY,NX)
      DP=DP+OSP(M,K,0,NY,NX)
      ELSE
      DCC=DCC+OSC(M,K,0,NY,NX)
      DNC=DNC+OSN(M,K,0,NY,NX)
      DPC=DPC+OSP(M,K,0,NY,NX)
      ENDIF
      RC0(K,NY,NX)=RC0(K,NY,NX)+OSC(M,K,0,NY,NX)
6930  CONTINUE
6900  CONTINUE
C
C     TOTAL ORGANIC MATTER
C
C     ORGC,ORGN,ORGP=total organic C,N,P
C     TLRSDC,TLRSDN,TLRSDP=total landscape litter C,N,P
C     URSDC,URSDN,URSDP=total grid cell litter C,N,P
C     TVOLWC,TVOLWP=canopy surface,internal water content
C        from uptake.f
C     VOLW,VOLV,VOLI=litter water,vapor,ice content
C     VOLWSO,UVOLW=total landscape, grid cell water content
C     HEATSO=total landscape heat content
C     TLCO2G,UCO2S=total landscape, grid cell CO2 content
C     OXYGSO=total landscape O2 content
C     TLN2G,TLNH4,TLNO3,TLPO4=total landscape N2,NH4,NO3,PO4 content
C     TION=total landscape ion content
C     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2
C        in litter
C     ZNH4S,ZNH3S,ZNO3S,ZNO2S=NH4,NH3,NO3,NO2 litter content
C     ZNH4FA,ZNHUFA,ZNH3FA,ZNO3FA=fertilizer NH4,urea,NH3,NO3
C        in litter
C     H2PO4,H1PO4=H2PO4,HPO4 litter content
C     XH2P,XH1P=exchangeable H2PO4,HPO4 in litter
C     PALPO,PFEPO=precip AlPO4,FEPO4
C     PCAPM,PCAPD,PCAPH=precip CaH2PO4,CaHPO4,apatite
C
      ORGC(0,NY,NX)=DC
      ORGN(0,NY,NX)=DN
      ORGR(0,NY,NX)=DC
      ORGCC(0,NY,NX)=DCC
      ORGNC(0,NY,NX)=DNC
      TLRSDC=TLRSDC+DC+DCC
      URSDC(NY,NX)=URSDC(NY,NX)+DC+DCC
      TLRSDN=TLRSDN+DN+DNC
      URSDN(NY,NX)=URSDN(NY,NX)+DN+DNC
      TLRSDP=TLRSDP+DP+DPC
      URSDP(NY,NX)=URSDP(NY,NX)+DP+DPC
      WSS=VOLW(0,NY,NX)+VOLV(0,NY,NX)+VOLI(0,NY,NX)*DENSI
      WSP=TVOLWC(NY,NX)+TVOLWP(NY,NX)
      VOLWSO=VOLWSO+WSS+WSP
      UVOLW(NY,NX)=UVOLW(NY,NX)+WSS+WSP
      HEATSO=HEATSO+TENGYC(NY,NX)
      CS=CO2S(0,NY,NX)+CH4S(0,NY,NX)+CO2G(0,NY,NX)+CH4G(0,NY,NX)
      TLCO2G=TLCO2G+CS
      UCO2S(NY,NX)=UCO2S(NY,NX)+CS
      HS=H2GS(0,NY,NX)+H2GG(0,NY,NX)
      TLH2G=TLH2G+HS
      OS=OXYS(0,NY,NX)+OXYG(0,NY,NX)
      OXYGSO=OXYGSO+OS
      ZGS=Z2GS(0,NY,NX)+Z2OS(0,NY,NX)+Z2GG(0,NY,NX)+Z2OG(0,NY,NX)
      TLN2G=TLN2G+ZGS
      Z4S=ZNH4S(0,NY,NX)+ZNH3S(0,NY,NX)+ZNH3G(0,NY,NX)
      Z4X=14.0*XN4(0,NY,NX)
      Z4F=14.0*(ZNH4FA(0,NY,NX)+ZNHUFA(0,NY,NX)+ZNH3FA(0,NY,NX))
      TLNH4=TLNH4+Z4S+Z4X+Z4F
      UNH4(NY,NX)=UNH4(NY,NX)+Z4S+Z4X
C     IF((I/30)*30.EQ.I.AND.J.EQ.24)THEN
C     DO 4342 K=0,2
C     WRITE(*,4341)'ORGC0',I,J,NFZ,NX,NY,K,ORGC(0,NY,NX),DN,DNC
C    2,((OMC(M,N,K,0,NY,NX),M=1,3),N=1,7)
C    3,(ORC(M,K,0,NY,NX),M=1,2)
C    4,OQC(K,0,NY,NX),OQCH(K,0,NY,NX),OHC(K,0,NY,NX)
C    2,OQA(K,0,NY,NX),OQAH(K,0,NY,NX),OHA(K,0,NY,NX)
C    5,(OSC(M,K,0,NY,NX),M=1,4)
C     WRITE(*,4341)'ORGN0',I,J,NFZ,NX,NY,K,DN,DNC,ON,ONC
C    2,ORGN(0,NY,NX),ORGNC(0,NY,NX)
C    3,((OMN(M,N,K,0,NY,NX),M=1,3),N=1,7)
C    3,(ORN(M,K,0,NY,NX),M=1,2)
C    4,(OSN(M,K,0,NY,NX),M=1,5)
C    4,OQN(K,0,NY,NX),OQNH(K,0,NY,NX),OHN(K,0,NY,NX)
4341  FORMAT(A8,6I4,120E14.6)
4342  CONTINUE
C     WRITE(*,5456)'TLCO2G0',I,J,NX,NY,TLCO2G
C    2,CS,CO2S(0,NY,NX),CH4S(0,NY,NX)
C     WRITE(*,5456)'TLN2G0',I,J,NX,NY,TLN2G
C    2,ZG,Z2GS(0,NY,NX),Z2OS(0,NY,NX)
C     WRITE(*,5456)'TLNH40',I,J,NX,NY,TLNH4,UNH4(NY,NX)
C    2,Z4S,Z4X,Z4F,XN4(0,NY,NX)
C    2,ZNH4S(0,NY,NX),ZNH3S(0,NY,NX)
5456  FORMAT(A8,4I4,30E16.6)
C     ENDIF
      ZOS=ZNO3S(0,NY,NX)+ZNO2S(0,NY,NX)
      ZOF=14.0*ZNO3FA(0,NY,NX)
      TLNO3=TLNO3+ZOS+ZOF
      UNO3(NY,NX)=UNO3(NY,NX)+ZOS
      P4S=H1PO4(0,NY,NX)+H2PO4(0,NY,NX)
      P4X=31.0*(XH1P(0,NY,NX)+XH2P(0,NY,NX))
      P4P=31.0*(PALPO(0,NY,NX)+PFEPO(0,NY,NX)
     2+PCAPD(0,NY,NX))
     3+62.0*PCAPM(0,NY,NX)
     4+93.0*PCAPH(0,NY,NX)
      TLPO4=TLPO4+P4S+P4X+P4P
      UPO4(NY,NX)=UPO4(NY,NX)+P4S
      UPX4(NY,NX)=UPX4(NY,NX)+P4X
      UPP4(NY,NX)=UPP4(NY,NX)+P4P
C
C     SURFCE LITTER SALT CONTENT
C
C     Z*=litter salt content
C     X*S=litter salt flux from trnsfrs.f
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C     TION,UION=total landscape, grid cell ion content
C
      IF(ISALTG.NE.0)THEN
      ZAL(0,NY,NX)=ZAL(0,NY,NX)+XALFLS(3,0,NY,NX)
      ZFE(0,NY,NX)=ZFE(0,NY,NX)+XFEFLS(3,0,NY,NX)
      ZHY(0,NY,NX)=ZHY(0,NY,NX)+XHYFLS(3,0,NY,NX)
      ZCA(0,NY,NX)=ZCA(0,NY,NX)+XCAFLS(3,0,NY,NX)
      ZMG(0,NY,NX)=ZMG(0,NY,NX)+XMGFLS(3,0,NY,NX)
      ZNA(0,NY,NX)=ZNA(0,NY,NX)+XNAFLS(3,0,NY,NX)
      ZKA(0,NY,NX)=ZKA(0,NY,NX)+XKAFLS(3,0,NY,NX)
      ZOH(0,NY,NX)=ZOH(0,NY,NX)+XOHFLS(3,0,NY,NX)
C     WRITE(20,3338)'ZHY0',I,J,ZHY(0,NY,NX),TQRHY(NY,NX)
C    2,XHYFLS(3,0,NY,NX),ZOH(0,NY,NX),TQROH(NY,NX)
C    4,XOHFLS(3,0,NY,NX)
      ZSO4(0,NY,NX)=ZSO4(0,NY,NX)+XSOFLS(3,0,NY,NX)
      ZCL(0,NY,NX)=ZCL(0,NY,NX)+XCLFLS(3,0,NY,NX)
      ZCO3(0,NY,NX)=ZCO3(0,NY,NX)+XC3FLS(3,0,NY,NX)
      ZHCO3(0,NY,NX)=ZHCO3(0,NY,NX)+XHCFLS(3,0,NY,NX)
      ZALOH1(0,NY,NX)=ZALOH1(0,NY,NX)+XAL1FS(3,0,NY,NX)
      ZALOH2(0,NY,NX)=ZALOH2(0,NY,NX)+XAL2FS(3,0,NY,NX)
      ZALOH3(0,NY,NX)=ZALOH3(0,NY,NX)+XAL3FS(3,0,NY,NX)
      ZALOH4(0,NY,NX)=ZALOH4(0,NY,NX)+XAL4FS(3,0,NY,NX)
      ZALS(0,NY,NX)=ZALS(0,NY,NX)+XALSFS(3,0,NY,NX)
      ZFEOH1(0,NY,NX)=ZFEOH1(0,NY,NX)+XFE1FS(3,0,NY,NX)
      ZFEOH2(0,NY,NX)=ZFEOH2(0,NY,NX)+XFE2FS(3,0,NY,NX)
      ZFEOH3(0,NY,NX)=ZFEOH3(0,NY,NX)+XFE3FS(3,0,NY,NX)
      ZFEOH4(0,NY,NX)=ZFEOH4(0,NY,NX)+XFE4FS(3,0,NY,NX)
      ZFES(0,NY,NX)=ZFES(0,NY,NX)+XFESFS(3,0,NY,NX)
      ZCAO(0,NY,NX)=ZCAO(0,NY,NX)+XCAOFS(3,0,NY,NX)
      ZCAC(0,NY,NX)=ZCAC(0,NY,NX)+XCACFS(3,0,NY,NX)
      ZCAH(0,NY,NX)=ZCAH(0,NY,NX)+XCAHFS(3,0,NY,NX)
      ZCAS(0,NY,NX)=ZCAS(0,NY,NX)+XCASFS(3,0,NY,NX)
      ZMGO(0,NY,NX)=ZMGO(0,NY,NX)+XMGOFS(3,0,NY,NX)
      ZMGC(0,NY,NX)=ZMGC(0,NY,NX)+XMGCFS(3,0,NY,NX)
      ZMGH(0,NY,NX)=ZMGH(0,NY,NX)+XMGHFS(3,0,NY,NX)
      ZMGS(0,NY,NX)=ZMGS(0,NY,NX)+XMGSFS(3,0,NY,NX)
      ZNAC(0,NY,NX)=ZNAC(0,NY,NX)+XNACFS(3,0,NY,NX)
      ZNAS(0,NY,NX)=ZNAS(0,NY,NX)+XNASFS(3,0,NY,NX)
      ZKAS(0,NY,NX)=ZKAS(0,NY,NX)+XKASFS(3,0,NY,NX)
      H0PO4(0,NY,NX)=H0PO4(0,NY,NX)+XH0PFS(3,0,NY,NX)
      H3PO4(0,NY,NX)=H3PO4(0,NY,NX)+XH3PFS(3,0,NY,NX)
      ZFE1P(0,NY,NX)=ZFE1P(0,NY,NX)+XF1PFS(3,0,NY,NX)
      ZFE2P(0,NY,NX)=ZFE2P(0,NY,NX)+XF2PFS(3,0,NY,NX)
      ZCA0P(0,NY,NX)=ZCA0P(0,NY,NX)+XC0PFS(3,0,NY,NX)
      ZCA1P(0,NY,NX)=ZCA1P(0,NY,NX)+XC1PFS(3,0,NY,NX)
      ZCA2P(0,NY,NX)=ZCA2P(0,NY,NX)+XC2PFS(3,0,NY,NX)
      ZMG1P(0,NY,NX)=ZMG1P(0,NY,NX)+XM1PFS(3,0,NY,NX)
      PSS=31.0*(H0PO4(0,NY,NX)+H3PO4(0,NY,NX)+ZFE1P(0,NY,NX)
     2+ZFE2P(0,NY,NX)+ZCA0P(0,NY,NX)+ZCA1P(0,NY,NX)
     3+ZCA2P(0,NY,NX)+ZMG1P(0,NY,NX))
      TLPO4=TLPO4+PSS
      SSS=ZAL(0,NY,NX)+ZFE(0,NY,NX)+ZHY(0,NY,NX)+ZCA(0,NY,NX)
     2+ZMG(0,NY,NX)+ZNA(0,NY,NX)+ZKA(0,NY,NX)+ZOH(0,NY,NX)
     3+ZSO4(0,NY,NX)+ZCL(0,NY,NX)+ZCO3(0,NY,NX)+H0PO4(0,NY,NX)
     4+2.0*(ZHCO3(0,NY,NX)+ZALOH1(0,NY,NX)+ZALS(0,NY,NX)
     5+ZFEOH1(0,NY,NX)+ZFES(0,NY,NX)+ZCAO(0,NY,NX)+ZCAC(0,NY,NX)
     6+ZCAS(0,NY,NX)+ZMGO(0,NY,NX)+ZMGC(0,NY,NX)+ZMGS(0,NY,NX)
     7+ZNAC(0,NY,NX)+ZNAS(0,NY,NX)+ZKAS(0,NY,NX)+ZCA0P(0,NY,NX))
     9+3.0*(ZALOH2(0,NY,NX)+ZFEOH2(0,NY,NX)+ZCAH(0,NY,NX)
     1+ZMGH(0,NY,NX)+ZFE1P(0,NY,NX)+ZCA1P(0,NY,NX)+ZMG1P(0,NY,NX))
     2+4.0*(ZALOH3(0,NY,NX)+ZFEOH3(0,NY,NX)+H3PO4(0,NY,NX)
     3+ZFE2P(0,NY,NX)+ZCA2P(0,NY,NX))
     4+5.0*(ZALOH4(0,NY,NX)+ZFEOH4(0,NY,NX))
      TION=TION+SSS
      UION(NY,NX)=UION(NY,NX)+SSS
      ENDIF
C     WRITE(20,3338)'SBN',I,J,TLNH4,TLNO3,TZIN,TZOU
C    2,Z4S,Z4X,Z4F,ZOS,ZOF,ZG
C    2,ZSI,ZGI,ZGB,Z2B,ZHB
C    3,ZXR,ZGR,ZOR,ZXB,TRN3S(0,NY,NX)
C    3,ZN4W(NY,NX),ZN3W(NY,NX),ZNOW(NY,NX),ZNGW(NY,NX),ZN2W(NY,NX)
C     WRITE(20,3338)'SBP',I,J,P4S,P4X,P4P,PSS,PXR,PXB,POR
C    2,TLPO4,TPIN,TPOU
C    2,XH1PS(0,NY,NX),XH2PS(0,NY,NX),H1PO4(0,NY,NX),H2PO4(0,NY,NX)
C    3,XH1P(0,NY,NX),XH2P(0,NY,NX),PALPO(0,NY,NX),PFEPO(0,NY,NX)
C    6,PCAPD(0,NY,NX),PCAPM(0,NY,NX),PCAPH(0,NY,NX),TRH1P(0,NY,NX)
C    2,TRH2P(0,NY,NX),XH1PFS(3,0,NY,NX),XH2PFS(3,0,NY,NX)
C    3,Z1PW(NY,NX),ZHPW(NY,NX),XH1PBS(NY,NX),XH2PBS(NY,NX)
C    4,FLQGQ(NY,NX),FLQRQ(NY,NX)
C     WRITE(20,3338)'SBS',I,J,TION,TIONIN,TIONOU
C    2,SSW,SSS,SIR,SII,SSR,SQE,SBU
3338  FORMAT(A8,2I4,50F20.8)
3335  FORMAT(A8,3I4,50F20.8)
C     ENDIF
C
C     UPDATE SOIL LAYER VARIABLES WITH TOTAL FLUXES
C
      TVHCP=0.0
      TVHCM=0.0
      TVOLW=0.0
      TVOLV=0.0
      TVOLWH=0.0
      TVOLI=0.0
      TVOLIH=0.0
      TENGY=0.0
      DO 125 L=NU(NY,NX),NL(NY,NX)
C
C     SOIL WATER, ICE, HEAT, TEMPERATURE
C
C     VOLW,VOLV,VOLI,VOLWH,VOLIH=water,vapor,ice content
C        in micropores,macropores
C     TFLW,TFLV,TFLWH,TFLWH=net micropore,macropore water,vapor flux,
C        heat flux
C     TWFLFL,TWFLFH=net freeze-thaw in micropores,macropores
C     TWFLVL=net evaporation-condensation
C     TVPFLG=net convective+diffusive vapor flux
C     FINH=micropore-macropore flux
C     TTHAW,TTHAWH=net freeze-thaw flux in micropores,macropores
C     TUPWTR=total root water uptake from extract.f
C     FLU=subsurface water input
C     DENSI=ice density from starts.f
C     DVOLW,DVOLI=change in soil water,ice content
C     VOLP=air-filled porosity
C     VHCP=soil heat capacity
C     TKS=soil temperature
C     THFLW=net soil conductive, convective heat flux
C     THFLFL=net freeze-thaw latent heat flux
C     THFLVL=net evaporation-condensation latent heat flux
C     TUPHT=total convective heat in root water uptake from extract.f
C     HWFLU=subsurface convective heat flux
C     HCBFX=heat released by soil combustion
C
      TKSX=TKS(L,NY,NX)
      VHCPX=VHCP(L,NY,NX)
      VOLWXX=VOLW(L,NY,NX)
      VOLIXX=VOLI(L,NY,NX)
      VOLW(L,NY,NX)=VOLW(L,NY,NX)+TFLW(L,NY,NX)+TWFLVL(L,NY,NX)
     2+TWFLFL(L,NY,NX)+FINH(L,NY,NX)+TUPWTR(L,NY,NX)+FLU(L,NY,NX)
      VOLV(L,NY,NX)=VOLV(L,NY,NX)+TFLV(L,NY,NX)-TWFLVL(L,NY,NX)
     2+TVPFLG(L,NY,NX)
      VOLWX(L,NY,NX)=VOLWX(L,NY,NX)+TFLWX(L,NY,NX)+FINH(L,NY,NX)
     2+TWFLVL(L,NY,NX)+TUPWTR(L,NY,NX)+FLU(L,NY,NX)
      VOLWX(L,NY,NX)=AMIN1(VOLW(L,NY,NX)
     2,VOLWX(L,NY,NX)+0.01*(VOLW(L,NY,NX)-VOLWX(L,NY,NX)))
      VOLI(L,NY,NX)=VOLI(L,NY,NX)-TWFLFL(L,NY,NX)/DENSI
      VOLWH(L,NY,NX)=VOLWH(L,NY,NX)+TFLWH(L,NY,NX)-FINH(L,NY,NX)
     2+TWFLFH(L,NY,NX)
      VOLIH(L,NY,NX)=VOLIH(L,NY,NX)-TWFLFH(L,NY,NX)/DENSI
      DVOLW(L,NY,NX)=VOLW1(L,NY,NX)+VOLWH1(L,NY,NX)
     2-VOLW(L,NY,NX)-VOLWH(L,NY,NX)
      DVOLI(L,NY,NX)=VOLI1(L,NY,NX)+VOLIH1(L,NY,NX)
     2-VOLI(L,NY,NX)-VOLIH(L,NY,NX)
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      VOLP(L,NY,NX)=AMAX1(0.0,VOLA(L,NY,NX)-VOLW(L,NY,NX)
     2-VOLI(L,NY,NX)+VOLAH(L,NY,NX)-VOLWH(L,NY,NX)-VOLIH(L,NY,NX))
      ELSE
      VOLP(L,NY,NX)=0.0
C     VOLA(L,NY,NX)=VOLW(L,NY,NX)+VOLI(L,NY,NX)
C    2+DVOLW(L,NY,NX)+DVOLI(L,NY,NX)
C     VOLX(L,NY,NX)=VOLA(L,NY,NX)
C     VOLT(L,NY,NX)=VOLA(L,NY,NX)
      ENDIF
      ENGY=VHCPX*TKSX
      VHCP(L,NY,NX)=VHCM(L,NY,NX)
     2+4.19*(VOLW(L,NY,NX)+VOLV(L,NY,NX)+VOLWH(L,NY,NX))
     3+1.9274*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
      TVHCP=TVHCP+VHCP(L,NY,NX)
      TVHCM=TVHCM+VHCM(L,NY,NX)
      TVOLW=TVOLW+VOLW(L,NY,NX)
      TVOLV=TVOLV+VOLV(L,NY,NX)
      TVOLWH=TVOLWH+VOLWH(L,NY,NX)
      TVOLI=TVOLI+VOLI(L,NY,NX)
      TVOLIH=TVOLIH+VOLIH(L,NY,NX)
      TENGY=TENGY+ENGY
C
C     ADD HEAT FLUX FOR ARTIFICIAL SOIL WARMING
C
C     IF(NX.EQ.3.AND.NY.EQ.2.AND.L.GT.NU(NY,NX)
C    3.AND.L.LE.17.AND.I.GE.152.AND.I.LE.304)THEN
C     THFLW(L,NY,NX)=THFLW(L,NY,NX)
C    2+(TKSZ(I,J,L)-TKS(L,NY,NX))*VHCP(L,NY,NX)*XNFH
C     WRITE(*,3379)'REDIST',I,J,NX,NY,L,TKSZ(I,J,L)
C    2,TKS(L,NY,NX),VHCP(L,NY,NX),THFLW(L,NY,NX)
3379  FORMAT(A8,6I4,12E12.4)
C     ENDIF
C
C     END ARTIFICIAL SOIL WARMING
C
      IF(VHCP(L,NY,NX).GT.VHCPNX(NY,NX))THEN
      TKS(L,NY,NX)=(ENGY+THFLW(L,NY,NX)+THFLFL(L,NY,NX)
     2+THFLVL(L,NY,NX)+TUPHT(L,NY,NX)+HWFLU(L,NY,NX)
     2+HCBFX(L,NY,NX))/VHCP(L,NY,NX)
      ELSE
      TKS(L,NY,NX)=TKS(L+1,NY,NX)
      ENDIF
      HEATIN=HEATIN+HCBFX(L,NY,NX)
      TCS(L,NY,NX)=TKS(L,NY,NX)-273.15
      TSMX(L,NY,NX)=AMAX1(TSMX(L,NY,NX),TCS(L,NY,NX))
      TSMN(L,NY,NX)=AMIN1(TSMN(L,NY,NX),TCS(L,NY,NX))
C     IF(L.EQ.9)THEN
C     WRITE(*,6547)'VOLW',I,J,NFZ,NX,NY,L,VOLW(L,NY,NX),VOLI(L,NY,NX)
C    2,VOLV(L,NY,NX),VOLP(L,NY,NX)
C    3,VOLWH(L,NY,NX),VOLIH(L,NY,NX),VOLAH(L,NY,NX),DVOLI(L,NY,NX)
C    4,TFLW(L,NY,NX),TWFLVL(L,NY,NX),TWFLFL(L,NY,NX)
C    3,FINH(L,NY,NX),TUPWTR(L,NY,NX),FLU(L,NY,NX),TQR(NY,NX)
C    5,TFLV(L,NY,NX),TWFLVL(L,NY,NX),TVPFLG(L,NY,NX)
C    6,XVPFLG(3,L,NY,NX),TEVAPG(NY,NX)
C    5,XVPFLG(3,L,NY,NX),XVPFLG(3,L+1,NY,NX)
C    2,FLW(3,L,NY,NX),FLW(3,L+1,NY,NX+1)
C    5,FLV(3,L,NY,NX),FLV(3,L+1,NY,NX)
C     WRITE(*,6547)'VOLWH',I,J,NFZ,NX,NY,L,VOLWH(L,NY,NX)
C    5,FLWH(3,L,NY,NX),FLWH(3,L+1,NY,NX)
C    2,TFLWH(L,NY,NX),FINH(L,NY,NX),TWFLFH(L,NY,NX)
C    3,VOLIH(L,NY,NX),VOLAH(L,NY,NX)
6547  FORMAT(A8,6I4,40E16.8)
C     WRITE(*,6633)'TKS',I,J,NFZ,NX,NY,L
C    2,TKS(L,NY,NX),ENGY,THFLW(L,NY,NX)
C    2,THFLFL(L,NY,NX),THFLVL(L,NY,NX),HWFLU(L,NY,NX),VHCP(L,NY,NX)
C    3,TVHCP,TVHCM,TVOLW,TVOLWH,TVOLI,TVOLIH,TENGY,TKSX,VHCPX
C    3,VOLWXX,VOLIXX,VOLW(L,NY,NX),VOLWH(L,NY,NX),VOLI(L,NY,NX)
C    4,VOLIH(L,NY,NX),TFLW(L,NY,NX),FINH(L,NY,NX),TWFLFL(L,NY,NX)
C    5,TUPWTR(L,NY,NX),FLU(L,NY,NX),TQR(NY,NX)
C    6,HFLW(3,L,NY,NX),HFLW(3,L+1,NY,NX)
C    7,ENGY+THFLW(L,NY,NX)+THFLFL(L,NY,NX)+TUPHT(L,NY,NX)
C    8+HWFLU(L,NY,NX)
6633  FORMAT(A8,6I4,40E14.6)
C     ENDIF
C
C     SOIL LITTER FROM PLANT ROOT LITTERFALL
C
C     OMC,OMN,OMP=microbial C,N,P
C     CSNT,ZSNT,PSNT=total C,N,P litterfall from extract.f
C
      DO 8565 K=0,1
      DO 8565 M=1,5
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)+CSNT(M,K,L,NY,NX)
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)+ZSNT(M,K,L,NY,NX)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)+PSNT(M,K,L,NY,NX)
C     IF(ZSNT(M,K,L,NY,NX).LT.0.0)THEN
C     WRITE(*,8484)'CSNT',I,J,NFZ,NX,NY,L,K,M,OSC(M,K,L,NY,NX)
C    2,OSN(M,K,L,NY,NX),OSP(M,K,L,NY,NX),CSNT(M,K,L,NY,NX)
C    3,ZSNT(M,K,L,NY,NX),PSNT(M,K,L,NY,NX)
8484  FORMAT(A8,8I4,12E14.6)
C     ENDIF
8565  CONTINUE
C
C     DOC,DON,DOP FROM AQUEOUS TRANSPORT
C
C     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores
C     OQCH,OQNH,OQPH,OQAH=DOC,DON,DOP,acetate in macropores
C     T*FLS,T*FHS=net solute flux in micropores,macropores
C     X*FXS=convective+diffusive solute flux between macropores
C        and micropores from trnsfr.f
C
C
      DO 8560 K=0,4
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+TOCFLS(K,L,NY,NX)
     2+XOCFXS(K,L,NY,NX)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+TONFLS(K,L,NY,NX)
     2+XONFXS(K,L,NY,NX)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+TOPFLS(K,L,NY,NX)
     2+XOPFXS(K,L,NY,NX)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+TOAFLS(K,L,NY,NX)
     2+XOAFXS(K,L,NY,NX)
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)+TOCFHS(K,L,NY,NX)
     2-XOCFXS(K,L,NY,NX)
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)+TONFHS(K,L,NY,NX)
     2-XONFXS(K,L,NY,NX)
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)+TOPFHS(K,L,NY,NX)
     2-XOPFXS(K,L,NY,NX)
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)+TOAFHS(K,L,NY,NX)
     2-XOAFXS(K,L,NY,NX)
C     IF(L.EQ.1)THEN
C     WRITE(*,2627)'OQCL',I,J,NFZ,NX,NY,L,K
C    2,OQC(K,L,NY,NX),TOCFLS(K,L,NY,NX),XOCFXS(K,L,NY,NX)
C    4,OQCH(K,L,NY,NX),TOCFHS(K,L,NY,NX),XOCFXS(K,L,NY,NX)
C    3,OQN(K,L,NY,NX),TONFLS(K,L,NY,NX),XONFXS(K,L,NY,NX)
C    4,OQNH(K,L,NY,NX),TONFHS(K,L,NY,NX),XONFXS(K,L,NY,NX)
C    5,OQA(K,L,NY,NX),TOAFLS(K,L,NY,NX),XOAFXS(K,L,NY,NX)
C    5,OQAH(K,L,NY,NX),TOAFHS(K,L,NY,NX),XOAFXS(K,L,NY,NX)
2627  FORMAT(A8,7I4,20E12.4)
C     ENDIF
8560  CONTINUE
C
C     SOIL DOC,DON,DOP FROM PLANT EXUDATION
C
C     OQC,OQN,OQP=DOC,DON,DOP in micropores
C     TDFOMC,TDFOMN,TDFOMP=root-soil nonstructural C,N,P exchange
C        from extract.f
C
      DO 195 K=0,4
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+TDFOMC(K,L,NY,NX)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+TDFOMN(K,L,NY,NX)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+TDFOMP(K,L,NY,NX)
195   CONTINUE
C
C     SOIL SOLUTES FROM AQUEOUS TRANSPORT, MICROBIAL AND ROOT
C     EXCHANGE, EQUILIBRIUM REACTIONS, GAS EXCHANGE,
C     MICROPORE-MACROPORE EXCHANGE
C
C
C     SOIL GAS SOLUTES
C
C     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2
C        in micropores
C     T*FLS=net solute flux in micropores
C     X*DFG=soil gas volatilization-dissolution from trnsfr.f
C     R*=net gas transformation in soil from nitro.f
C     TUP*S,TUP*B=total root-soil gas, solute exchange in
C        non-band,band from extract.f
C     X*BBL=bubble flux from trnsfr.f
C     TR*S=net solute transformation from solute.f
C     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
C             :*ZN3*=NH3,*H2G*=H2
C
      CO2S(L,NY,NX)=CO2S(L,NY,NX)+TCOFLS(L,NY,NX)+XCODFG(L,NY,NX)
     2-RCO2O(L,NY,NX)-TCO2S(L,NY,NX)+RCOFLU(L,NY,NX)+XCOFXS(L,NY,NX)
     3+TRCO2(L,NY,NX)+XCOBBL(L,NY,NX)
      CH4S(L,NY,NX)=CH4S(L,NY,NX)+TCHFLS(L,NY,NX)+XCHDFG(L,NY,NX)
     2-RCH4O(L,NY,NX)-TUPCHS(L,NY,NX)+RCHFLU(L,NY,NX)
     3+XCHFXS(L,NY,NX)+XCHBBL(L,NY,NX)
      OXYS(L,NY,NX)=OXYS(L,NY,NX)+TOXFLS(L,NY,NX)+XOXDFG(L,NY,NX)
     2-RUPOXO(L,NY,NX)-TUPOXS(L,NY,NX)+ROXFLU(L,NY,NX)
     3+XOXFXS(L,NY,NX)+XOXBBL(L,NY,NX)
      Z2GS(L,NY,NX)=Z2GS(L,NY,NX)+TNGFLS(L,NY,NX)+XNGDFG(L,NY,NX)
     2-RN2G(L,NY,NX)-TUPNF(L,NY,NX)+RNGFLU(L,NY,NX)+XNGFXS(L,NY,NX)
     3-XN2GS(L,NY,NX)+XNGBBL(L,NY,NX)
      Z2OS(L,NY,NX)=Z2OS(L,NY,NX)+TN2FLS(L,NY,NX)+XN2DFG(L,NY,NX)
     2-RN2O(L,NY,NX)-TUPN2S(L,NY,NX)+RN2FLU(L,NY,NX)+XN2FXS(L,NY,NX)
     3+XN2BBL(L,NY,NX)
      H2GS(L,NY,NX)=H2GS(L,NY,NX)+THGFLS(L,NY,NX)+XHGDFG(L,NY,NX)
     2-RH2GO(L,NY,NX)-TUPHGS(L,NY,NX)+RHGFLU(L,NY,NX)
     3+XHGFXS(L,NY,NX)+XHGBBL(L,NY,NX)
C     IF(L.EQ.5)THEN
C     WRITE(*,5432)'CO2SL',I,J,NFZ,NX,NY,L
C    2,CO2S(L,NY,NX),TCOFLS(L,NY,NX)
C    2,XCODFG(L,NY,NX),RCO2O(L,NY,NX),TCO2S(L,NY,NX)
C    3,RCOFLU(L,NY,NX),XCOFXS(L,NY,NX),TRCO2(L,NY,NX)
C    4,XCOBBL(L,NY,NX),CO2G(L,NY,NX),ORGC(L,NY,NX)
C     WRITE(*,5432)'CH4SL',I,J,NFZ,NX,NY,L,CH4S(L,NY,NX),CH4G(L,NY,NX)
C    2,TCHFLS(L,NY,NX),XCHDFG(L,NY,NX),RCH4O(L,NY,NX),TUPCHS(L,NY,NX)
C    3,RCHFLU(L,NY,NX),XCHFXS(L,NY,NX),XCHBBL(L,NY,NX)
C    4,XCOBBL(L,NY,NX),XCHFLS(3,L,NY,NX),XCHFLS(3,L+1,NY,NX)
C    5,CCH4S(L,NY,NX)
C     WRITE(*,5432)'OXYSL',I,J,NFZ,NX,NY,L
C    2,OXYS(L,NY,NX),OXYG(L,NY,NX)
C    4,XOXFLS(3,L,NY,NX),XOXFLS(3,L+1,NY,NX),TOXFLS(L,NY,NX)
C    2,XOXDFG(L,NY,NX),RUPOXO(L,NY,NX),TUPOXS(L,NY,NX)
C    3,ROXFLU(L,NY,NX),XOXFXS(L,NY,NX),XOXBBL(L,NY,NX)
C    5,XOXDFS(NY,NX),COXYS(L,NY,NX)
C     WRITE(*,5432)'Z2OS',I,J,NFZ,NX,NY,L
C    2,Z2OS(L,NY,NX),TN2FLS(L,NY,NX)
C    2,XN2DFG(L,NY,NX),RN2O(L,NY,NX),TUPN2S(L,NY,NX),RN2FLU(L,NY,NX)
C    3,XN2FXS(L,NY,NX),XN2BBL(L,NY,NX),XN2FLS(3,L,NY,NX)
C    4,XN2FLS(3,L+1,NY,NX),XN2DFS(NY,NX)
C    3,Z2GS(L,NY,NX),TNGFLS(L,NY,NX),XNGDFG(L,NY,NX)
C    4,RN2G(L,NY,NX),TUPNF(L,NY,NX),RNGFLU(L,NY,NX),XNGFXS(L,NY,NX)
C    5,XN2GS(L,NY,NX),XNGBBL(L,NY,NX),XNGDFS(NY,NX)
C     WRITE(*,5432)'H2GS',I,J,NFZ,NX,NY,L
C    2,H2GS(L,NY,NX),THGFLS(L,NY,NX)
C    2,XHGDFG(L,NY,NX),RH2GO(L,NY,NX),TUPHGS(L,NY,NX),RHGFLU(L,NY,NX)
C    3,XHGFXS(L,NY,NX),XHGBBL(L,NY,NX),XHGDFS(NY,NX)
5432  FORMAT(A8,6I4,20E12.4)
C     ENDIF
C
C     SOIL MINERAL N,P SOLUTES NON-BAND
C
C     ZNH4S,ZNH3S,ZNO3S,ZNO2S=NH4,NH3,NO3,NO2 in non-band micropores
C     ZNH4B,ZNH3B,ZNO3B,ZNO2B=NH4,NH3,NO3,NO2 in band micropores
C     H1PO4,H2PO4=HPO4,H2PO4 in non-band micropores
C     H1POB,H2POB=HPO4,H2PO4 in band micropores
C     T*FLS=net solute flux in micropores
C     TUP*S,TUP*B=total root-soil gas, solute exchange in
C        non-band,band from extract.f
C     R*FLU,R*FBU=subsurface solute flux in non-band,band
C     X*FXS,X*FXB=convective+diffusive solute flux between macropores
C        and micropores in non-band,band from trnsfr.f
C     TR*S=net solute transformation from solute.f
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,H1P=HPO4,H2P=H2PO4 in non-band
C        :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,H1B=HPO4,H2B=H2PO4 in band
C
      ZNH3S(L,NY,NX)=ZNH3S(L,NY,NX)+TN3FLS(L,NY,NX)+XN3DFG(L,NY,NX)
     2+TRN3S(L,NY,NX)-TUPN3S(L,NY,NX)+RN3FLU(L,NY,NX)
     3+XN3FXW(L,NY,NX)+XN3BBL(L,NY,NX)
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+TN4FLS(L,NY,NX)+XNH4S(L,NY,NX)
     2+TRN4S(L,NY,NX)-TUPNH4(L,NY,NX)+RN4FLU(L,NY,NX)
     3+XN4FXW(L,NY,NX)
      ZNO3S(L,NY,NX)=ZNO3S(L,NY,NX)+TNOFLS(L,NY,NX)+XNO3S(L,NY,NX)
     2+TRNO3(L,NY,NX)-TUPNO3(L,NY,NX)+RNOFLU(L,NY,NX)
     3+XNOFXW(L,NY,NX)
      ZNO2S(L,NY,NX)=ZNO2S(L,NY,NX)+TNXFLS(L,NY,NX)+XNO2S(L,NY,NX)
     2+TRNO2(L,NY,NX)+XNXFXS(L,NY,NX)
C     IF(L.EQ.NU(NY,NX))THEN
C     WRITE(*,4444)'NH3',I,J,NX,NY,L,ZNH3S(L,NY,NX),TN3FLS(L,NY,NX)
C    2,XN3DFG(L,NY,NX),TRN3S(L,NY,NX),TUPN3S(L,NY,NX)
C    3,RN3FLU(L,NY,NX),XN3FXW(L,NY,NX),XN3BBL(L,NY,NX),XN3DFS(NY,NX)
C    4,ZNH4S(L,NY,NX)
C    4,TN4FLS(L,NY,NX),XNH4S(L,NY,NX),TRN4S(L,NY,NX),TUPNH4(L,NY,NX)
C    5,RN4FLU(L,NY,NX),XN4FXW(L,NY,NX),TN4QRS(NY,NX),TN3QRS(NY,NX)
C    6,ZNH3SH(L,NY,NX),ZNH4SH(L,NY,NX),14.0*XN4(L,NY,NX)
4443  FORMAT(A8,5I4,30F16.8)
4444  FORMAT(A8,5I4,30E12.4)
C     WRITE(*,5545)'NO3',I,J,NX,NY,L,ZNO3S(L,NY,NX),TNOFLS(L,NY,NX)
C    2,XNO3S(L,NY,NX),TRNO3(L,NY,NX),TUPNO3(L,NY,NX),RNOFLU(L,NY,NX)
C    3,XNOFXW(L,NY,NX),ZNO2S(L,NY,NX),TNXFLS(L,NY,NX)
C    4,XNO2S(L,NY,NX),TRNO2(L,NY,NX),XNXFXS(L,NY,NX),TNXQRS(NY,NX)
5545  FORMAT(A8,5I4,40E12.4)
C     ENDIF
      H1PO4(L,NY,NX)=H1PO4(L,NY,NX)+TP1FLS(L,NY,NX)+XH1PS(L,NY,NX)
     2+TRH1P(L,NY,NX)-TUPH1P(L,NY,NX)+RH1PFU(L,NY,NX)+XH1PXS(L,NY,NX)
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+TPOFLS(L,NY,NX)+XH2PS(L,NY,NX)
     2+TRH2P(L,NY,NX)-TUPH2P(L,NY,NX)+RH2PFU(L,NY,NX)+XH2PXS(L,NY,NX)
C
C     SOIL MINERAL N,P SOLUTES BAND
C
      ZNH3B(L,NY,NX)=ZNH3B(L,NY,NX)+TN3FLB(L,NY,NX)+XNBDFG(L,NY,NX)
     2+TRN3B(L,NY,NX)-TUPN3B(L,NY,NX)+RN3FBU(L,NY,NX)
     3+XN3FXB(L,NY,NX)+XNBBBL(L,NY,NX)
      ZNH4B(L,NY,NX)=ZNH4B(L,NY,NX)+TN4FLB(L,NY,NX)+XNH4B(L,NY,NX)
     2+TRN4B(L,NY,NX)-TUPNHB(L,NY,NX)+RN4FBU(L,NY,NX)
     3+XN4FXB(L,NY,NX)
      ZNO3B(L,NY,NX)=ZNO3B(L,NY,NX)+TNOFLB(L,NY,NX)+XNO3B(L,NY,NX)
     2+TRNOB(L,NY,NX)-TUPNOB(L,NY,NX)+RNOFBU(L,NY,NX)
     3+XNOFXB(L,NY,NX)
      ZNO2B(L,NY,NX)=ZNO2B(L,NY,NX)+TNXFLB(L,NY,NX)+XNO2B(L,NY,NX)
     2+TRN2B(L,NY,NX)+XNXFXB(L,NY,NX)
      H1POB(L,NY,NX)=H1POB(L,NY,NX)+TH1BFB(L,NY,NX)+XH1BS(L,NY,NX)
     2+TRH1B(L,NY,NX)-TUPH1B(L,NY,NX)+RH1BBU(L,NY,NX)
     3+XH1BXB(L,NY,NX)
      H2POB(L,NY,NX)=H2POB(L,NY,NX)+TH2BFB(L,NY,NX)+XH2BS(L,NY,NX)
     2+TRH2B(L,NY,NX)-TUPH2B(L,NY,NX)+RH2BBU(L,NY,NX)
     3+XH2BXB(L,NY,NX)
      THRE(NY,NX)=THRE(NY,NX)+RCO2O(L,NY,NX)+RCH4O(L,NY,NX)
      UN2GS(NY,NX)=UN2GS(NY,NX)+XN2GS(L,NY,NX)
      UN2GG(NY,NX)=UN2GG(NY,NX)+RN2G(L,NY,NX)
      HN2GG(NY,NX)=HN2GG(NY,NX)+RN2G(L,NY,NX)
C
C     SOIL EXCHANGEABLE CATIONS AND ANIONS FROM EXCHANGE REACTIONS
C
C     TR*=total exchange+precipitation from solute.f
C     XN4,XH1P,XH2P=exchangeable NH4,HPO4,H2PO4 in non-band
C     XOH0,XOH1,XOH2=adsorption sites R-,R-OH,R-OH2 in non-band
C     XNB,XH1PB,XH2PB=exchangeable NH4,HPO4,H2PO4 in band
C     XOH0B,XOH1B,XOH2B=adsorption sites R-,R-OH,R-OH2 in band
C
      XN4(L,NY,NX)=XN4(L,NY,NX)+TRXN4(L,NY,NX)
      XNB(L,NY,NX)=XNB(L,NY,NX)+TRXNB(L,NY,NX)
      XOH0(L,NY,NX)=XOH0(L,NY,NX)+TRXH0(L,NY,NX)
      XOH1(L,NY,NX)=XOH1(L,NY,NX)+TRXH1(L,NY,NX)
      XOH2(L,NY,NX)=XOH2(L,NY,NX)+TRXH2(L,NY,NX)
      XH1P(L,NY,NX)=XH1P(L,NY,NX)+TRX1P(L,NY,NX)
      XH2P(L,NY,NX)=XH2P(L,NY,NX)+TRX2P(L,NY,NX)
      XOH0B(L,NY,NX)=XOH0B(L,NY,NX)+TRBH0(L,NY,NX)
      XOH1B(L,NY,NX)=XOH1B(L,NY,NX)+TRBH1(L,NY,NX)
      XOH2B(L,NY,NX)=XOH2B(L,NY,NX)+TRBH2(L,NY,NX)
      XH1PB(L,NY,NX)=XH1PB(L,NY,NX)+TRB1P(L,NY,NX)
      XH2PB(L,NY,NX)=XH2PB(L,NY,NX)+TRB2P(L,NY,NX)
C     IF(J.EQ.12.AND.L.LE.4)THEN
C     WRITE(*,4445)'NHB',I,J,NX,NY,L,ZNH3B(L,NY,NX),TN3FLB(L,NY,NX)
C    2,XNBDFG(L,NY,NX),TRN3B(L,NY,NX),TUPN3B(L,NY,NX)
C    3,RN3FBU(L,NY,NX),XN3FXB(L,NY,NX),XNBBBL(L,NY,NX),TUPNHB(L,NY,NX)
C    4,ZNH4B(L,NY,NX),TN4FLB(L,NY,NX),XNH4B(L,NY,NX)
C    5,TRN4B(L,NY,NX),TUPNHB(L,NY,NX),RN4FBU(L,NY,NX)
C    6,XNB(L,NY,NX)*14.0
C     WRITE(*,4445)'NOB',I,J,NX,NY,L,ZNO2B(L,NY,NX),TNXFLB(L,NY,NX)
C    2,XNO2B(L,NY,NX),TRN2B(L,NY,NX),XNXFXB(L,NY,NX)
4445  FORMAT(A8,5I4,20E12.4)
C     ENDIF
C
C     SOIL PRECIPITATES FROM PRECIPITATION-DISSOLUTION REACTIONS
C
C     TR*=total exchange+precipitation from solute.f
C     PALPO,PFEPO=precipitated AlPO4,FEPO4 in non-band
C     PCAPM,PCAPD,PCAPH=precipitated CaH2PO4,CaHPO4,apatite
C        in non-band
C     PALPB,PFEPB=precipitated AlPO4,FEPO4 in band
C     PCPMB,PCPDB,PCPHB=precipitated CaH2PO4,CaHPO4,apatite in band
C
      PALPO(L,NY,NX)=PALPO(L,NY,NX)+TRALPO(L,NY,NX)
      PFEPO(L,NY,NX)=PFEPO(L,NY,NX)+TRFEPO(L,NY,NX)
      PCAPD(L,NY,NX)=PCAPD(L,NY,NX)+TRCAPD(L,NY,NX)
      PCAPH(L,NY,NX)=PCAPH(L,NY,NX)+TRCAPH(L,NY,NX)
      PCAPM(L,NY,NX)=PCAPM(L,NY,NX)+TRCAPM(L,NY,NX)
      PALPB(L,NY,NX)=PALPB(L,NY,NX)+TRALPB(L,NY,NX)
      PFEPB(L,NY,NX)=PFEPB(L,NY,NX)+TRFEPB(L,NY,NX)
      PCPDB(L,NY,NX)=PCPDB(L,NY,NX)+TRCPDB(L,NY,NX)
      PCPHB(L,NY,NX)=PCPHB(L,NY,NX)+TRCPHB(L,NY,NX)
      PCPMB(L,NY,NX)=PCPMB(L,NY,NX)+TRCPMB(L,NY,NX)
C     IF(PCAPD(L,NY,NX).LT.ZEROS(NY,NX))THEN
C     WRITE(*,6638)'PCAPD',I,J,NFZ,L,PCAPD(L,NY,NX),TRCAPD(L,NY,NX)
6638  FORMAT(A8,4I4,12E12.4)
C     ENDIF
C
C     MACROPORE SOLUTES FROM MACROPORE-MICROPORE EXCHANGE
C
C     CO2SH,CH4SH,OXYSH,Z2GSH,Z2OSH,H2GSH=aqueous CO2,CH4,O2,N2,N2O,H2
C        in macropores
C     ZNH4SH,ZNH3SH,ZNO3SH,ZNO2SH=NH4,NH3,NO3,NO2
C        in non-band macropores
C     ZNH4BH,ZNH3BH,ZNO3BH,ZNO2BH=NH4,NH3,NO3,NO2 in band macropores
C     H1PO4H,H2PO4H=HPO4,H2PO4 in non-band macropores
C     H1POBH,H2POBH=HPO4,H2PO4 in band macropores
C     T*FHS=net solute flux in macropores
C     X*FXS,X*FXB=convective+diffusive solute flux between macropores
C        and micropores in non-band,band from trnsfr.f
C
      CO2SH(L,NY,NX)=CO2SH(L,NY,NX)+TCOFHS(L,NY,NX)-XCOFXS(L,NY,NX)
      CH4SH(L,NY,NX)=CH4SH(L,NY,NX)+TCHFHS(L,NY,NX)-XCHFXS(L,NY,NX)
      OXYSH(L,NY,NX)=OXYSH(L,NY,NX)+TOXFHS(L,NY,NX)-XOXFXS(L,NY,NX)
      Z2GSH(L,NY,NX)=Z2GSH(L,NY,NX)+TNGFHS(L,NY,NX)-XNGFXS(L,NY,NX)
      Z2OSH(L,NY,NX)=Z2OSH(L,NY,NX)+TN2FHS(L,NY,NX)-XN2FXS(L,NY,NX)
      H2GSH(L,NY,NX)=H2GSH(L,NY,NX)+THGFHS(L,NY,NX)-XHGFXS(L,NY,NX)
      ZNH4SH(L,NY,NX)=ZNH4SH(L,NY,NX)+TN4FHS(L,NY,NX)-XN4FXW(L,NY,NX)
      ZNH3SH(L,NY,NX)=ZNH3SH(L,NY,NX)+TN3FHS(L,NY,NX)-XN3FXW(L,NY,NX)
      ZNO3SH(L,NY,NX)=ZNO3SH(L,NY,NX)+TNOFHS(L,NY,NX)-XNOFXW(L,NY,NX)
      ZNO2SH(L,NY,NX)=ZNO2SH(L,NY,NX)+TNXFHS(L,NY,NX)-XNXFXS(L,NY,NX)
      H1PO4H(L,NY,NX)=H1PO4H(L,NY,NX)+TP1FHS(L,NY,NX)-XH1PXS(L,NY,NX)
      H2PO4H(L,NY,NX)=H2PO4H(L,NY,NX)+TPOFHS(L,NY,NX)-XH2PXS(L,NY,NX)
      ZNH4BH(L,NY,NX)=ZNH4BH(L,NY,NX)+TN4FHB(L,NY,NX)-XN4FXB(L,NY,NX)
      ZNH3BH(L,NY,NX)=ZNH3BH(L,NY,NX)+TN3FHB(L,NY,NX)-XN3FXB(L,NY,NX)
      ZNO3BH(L,NY,NX)=ZNO3BH(L,NY,NX)+TNOFHB(L,NY,NX)-XNOFXB(L,NY,NX)
      ZNO2BH(L,NY,NX)=ZNO2BH(L,NY,NX)+TNXFHB(L,NY,NX)-XNXFXB(L,NY,NX)
      H1POBH(L,NY,NX)=H1POBH(L,NY,NX)+TH1BHB(L,NY,NX)-XH1BXB(L,NY,NX)
      H2POBH(L,NY,NX)=H2POBH(L,NY,NX)+TH2BHB(L,NY,NX)-XH2BXB(L,NY,NX)
C     IF(L.EQ.4)THEN
C     WRITE(*,4747)'ZNH4SH',I,J,NX,NY,L
C    2,ZNH4SH(L,NY,NX),TN4FHS(L,NY,NX),XN4FXW(L,NY,NX)
C    2,ZNH3SH(L,NY,NX),TN3FHS(L,NY,NX),XN3FXW(L,NY,NX)
C     WRITE(*,4747)'ZNO3SH',I,J,NX,NY,L,ZNO3SH(L,NY,NX)
C    2,TNOFHS(L,NY,NX),XNOFXW(L,NY,NX)
C    3,ZNO2SH(L,NY,NX),TNXFHS(L,NY,NX),XNXFXS(L,NY,NX)
4747  FORMAT(A8,5I4,12E12.4)
C     IF((I/30)*30.EQ.I.AND.J.EQ.24)THEN
C     WRITE(*,5545)'HP14',I,J,NX,NY,L,H1PO4(L,NY,NX),TP1FLS(L,NY,NX)
C    2,XH1PS(L,NY,NX),TRH1P(L,NY,NX),TUPH1P(L,NY,NX),RH1PFU(L,NY,NX)
C    3,XH1PXS(L,NY,NX),XH1P(L,NY,NX),H1POB(L,NY,NX)
C    4,TH1BFB(L,NY,NX),XH1BS(L,NY,NX),TRH1B(L,NY,NX),TUPH1B(L,NY,NX)
C    2,RH1BBU(L,NY,NX),XH1BXB(L,NY,NX),XH1PB(L,NY,NX)
C    2,H1PO4H(L,NY,NX),TP1FHS(L,NY,NX),XH1PXS(L,NY,NX)
C    2,H1POBH(L,NY,NX),TH1BHB(L,NY,NX),XH1BXB(L,NY,NX)
C     WRITE(*,5545)'HP24',I,J,NX,NY,L,H2PO4(L,NY,NX),TPOFLS(L,NY,NX)
C    2,XH2PS(L,NY,NX),TRH2P(L,NY,NX),TUPH2P(L,NY,NX),RH2PFU(L,NY,NX)
C    3,XH2PXS(L,NY,NX),XH2P(L,NY,NX),H2POB(L,NY,NX)
C    4,TH2BFB(L,NY,NX),XH2BS(L,NY,NX),TRH2B(L,NY,NX),TUPH2B(L,NY,NX)
C    5,RH2BBU(L,NY,NX),XH2BXB(L,NY,NX),XH2PB(L,NY,NX)
C    2,H2PO4H(L,NY,NX),TPOFHS(L,NY,NX),XH2PXS(L,NY,NX)
C    2,H2POBH(L,NY,NX),TH2BHB(L,NY,NX),XH2BXB(L,NY,NX)
C     ENDIF
C     ENDIF
C
C     SOIL GASES FROM VOLATILIZATION-DISSOLUTION AND GAS TRANSFER
C
C     *G=soil gas content
C     T*FLG=net convective+diffusive gas flux from trnsfr.f
C     X*DFG=net water-air gas flux from trnsfr.f
C     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
C             :*ZN3*=NH3,*H2G*=H2
C     TRN3G=NH3 dissolution from solute.f
C     R*F=net gaseous flux for use in nitro.f, uptake.f
C     R*L=net aqueous flux for use in nitro.f, uptake.f
C     X*BBL=bubble flux from trnsfr.f
C
      CO2G(L,NY,NX)=CO2G(L,NY,NX)+TCOFLG(L,NY,NX)-XCODFG(L,NY,NX)
     2+RCGOX(L,NY,NX)+RC4OX(L,NY,NX)
      CH4G(L,NY,NX)=CH4G(L,NY,NX)+TCHFLG(L,NY,NX)-XCHDFG(L,NY,NX)
     2+RCHOX(L,NY,NX)-RC4OX(L,NY,NX)
      OXYG(L,NY,NX)=OXYG(L,NY,NX)+TOXFLG(L,NY,NX)-XOXDFG(L,NY,NX)
     2-ROGOX(L,NY,NX)-RC4OX(L,NY,NX)*2.667
      Z2GG(L,NY,NX)=Z2GG(L,NY,NX)+TNGFLG(L,NY,NX)-XNGDFG(L,NY,NX)
      Z2OG(L,NY,NX)=Z2OG(L,NY,NX)+TN2FLG(L,NY,NX)-XN2DFG(L,NY,NX)
      ZNH3G(L,NY,NX)=ZNH3G(L,NY,NX)+TNHFLG(L,NY,NX)-XN3DFG(L,NY,NX)
     2-XNBDFG(L,NY,NX)+TRN3G(L,NY,NX)
      H2GG(L,NY,NX)=H2GG(L,NY,NX)+THGFLG(L,NY,NX)-XHGDFG(L,NY,NX)
      ROXYF(L,NY,NX)=TOXFLG(L,NY,NX)
      RCO2F(L,NY,NX)=TCOFLG(L,NY,NX)
      RCH4F(L,NY,NX)=TCHFLG(L,NY,NX)
      ROXYL(L,NY,NX)=TOXFLS(L,NY,NX)+ROXFLU(L,NY,NX)+XOXFXS(L,NY,NX)
     2+XOXBBL(L,NY,NX)
      RCH4L(L,NY,NX)=TCHFLS(L,NY,NX)+RCHFLU(L,NY,NX)+XCHFXS(L,NY,NX)
     2+XCHBBL(L,NY,NX)
C     IF((I/30)*30.EQ.I.AND.J.EQ.15.AND.L.EQ.1)THEN
C     WRITE(*,5432)'CO2GL',I,J,NFZ,NX,NY,L
C    2,CO2G(L,NY,NX),TCOFLG(L,NY,NX)
C    2,XCODFG(L,NY,NX),THETP(L,NY,NX)
C     WRITE(*,5432)'OXYGL',I,J,NFZ,NX,NY,L
C    2,OXYG(L,NY,NX),TOXFLG(L,NY,NX)
C    2,XOXDFG(L,NY,NX),COXYG(L,NY,NX),XOXFLG(3,L,NY,NX)
C    3,XOXFLG(3,L+1,NY,NX),XOXFLG(1,L,NY,NX+1)
C     WRITE(*,5432)'CH4GL',I,J,NFZ,NX,NY,L
C    2,CH4G(L,NY,NX),TCHFLG(L,NY,NX)
C    2,XCHDFG(L,NY,NX),CCH4G(L,NY,NX),XCHFLG(3,L,NY,NX)
C    3,XCHFLG(3,L+1,NY,NX),XCHDFS(NY,NX),RCH4F(L,NY,NX)
C    4,RCH4L(L,NY,NX),TCHFLS(L,NY,NX),RCHFLU(L,NY,NX)
C    5,XCHFXS(L,NY,NX),XCHBBL(L,NY,NX)
C     ENDIF
C
C     GRID CELL BOUNDARY FLUXES FROM ROOT GAS TRANSFER
C
C     HEATIN=cumulative net surface heat transfer
C     THTHAW=net freeze-thaw latent heat flux
C     TUPHT=total convective heat in root water uptake from extract.f
C     T*FLA=total root gaseous-atmosphere gas exchange from extract.f
C     X*BBL=bubble flux from trnsfr.f
C     *G=soil gas content
C     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
C     CO2GIN,TCOU=cumulative surface,subsurface gas C exchange
C     OXYGIN,OXYGOU=cumulative surface,subsurface gas O2 exchange
C     H2GIN,H2GOU=cumulative surface,subsurface gas H2 exchange
C     TZIN,TZOU=cumulative surface,subsurface aqueous NH4,NH3,NO3
C        exchange
C     H*G,U*G=current,cumulative gas exchange
C     T*S,T*P=total root gas exchange with soil, root
C     TXCO2=CO2 net change from all solute equilibria
C
      HEATIN=HEATIN+THFLFL(L,NY,NX)+THFLVL(L,NY,NX)+TUPHT(L,NY,NX)
      CIB=TCOFLA(L,NY,NX)
      CHB=TCHFLA(L,NY,NX)
      OIB=TOXFLA(L,NY,NX)
C     HGB=THGFLA(L,NY,NX)
      HGB=0.0
      ZGB=0.0
      Z2B=TN2FLA(L,NY,NX)
      ZHB=TNHFLA(L,NY,NX)
C
C     GRID CELL BOUNDARY FLUXES BUBBLING
C
      IF(LG.EQ.0)THEN
      LL=0
      CIB=CIB+XCOBBL(L,NY,NX)
      CHB=CHB+XCHBBL(L,NY,NX)
      OIB=OIB+XOXBBL(L,NY,NX)
      ZGB=ZGB+XNGBBL(L,NY,NX)
      Z2B=Z2B+XN2BBL(L,NY,NX)
      ZHB=ZHB+XN3BBL(L,NY,NX)+XNBBBL(L,NY,NX)
      HGB=HGB+XHGBBL(L,NY,NX)
      ELSE
      LL=MIN(L,LG)
      CO2G(LL,NY,NX)=CO2G(LL,NY,NX)-XCOBBL(L,NY,NX)
      CH4G(LL,NY,NX)=CH4G(LL,NY,NX)-XCHBBL(L,NY,NX)
      OXYG(LL,NY,NX)=OXYG(LL,NY,NX)-XOXBBL(L,NY,NX)
      Z2GG(LL,NY,NX)=Z2GG(LL,NY,NX)-XNGBBL(L,NY,NX)
      Z2OG(LL,NY,NX)=Z2OG(LL,NY,NX)-XN2BBL(L,NY,NX)
      ZNH3G(LL,NY,NX)=ZNH3G(LL,NY,NX)-XN3BBL(L,NY,NX)-XNBBBL(L,NY,NX)
      H2GG(LL,NY,NX)=H2GG(LL,NY,NX)-XHGBBL(L,NY,NX)
      IF(LG.LT.L)THEN
      TLCO2G=TLCO2G-XCOBBL(L,NY,NX)-XCHBBL(L,NY,NX)
      UCO2S(NY,NX)=UCO2S(NY,NX)-XCOBBL(L,NY,NX)-XCHBBL(L,NY,NX)
      OXYGSO=OXYGSO-XOXBBL(L,NY,NX)
      TLN2G=TLN2G-XNGBBL(L,NY,NX)-XN2BBL(L,NY,NX)
     2-XN3BBL(L,NY,NX)-XNBBBL(L,NY,NX)
      TLH2G=TLH2G-XHGBBL(L,NY,NX)
      ENDIF
      ENDIF
      CO2GIN=CO2GIN+CIB+CHB
      COB=TCO2P(L,NY,NX)+TCO2S(L,NY,NX)-TRCO2(L,NY,NX)
      TCOU=TCOU+COB
      XCNET(NY,NX)=XCNET(NY,NX)+CIB
      XHNET(NY,NX)=XHNET(NY,NX)+CHB
      HCO2G(NY,NX)=HCO2G(NY,NX)+CIB
      UCO2G(NY,NX)=UCO2G(NY,NX)+CIB
      HCH4G(NY,NX)=HCH4G(NY,NX)+CHB
      UCH4G(NY,NX)=UCH4G(NY,NX)+CHB
      UCOP(NY,NX)=UCOP(NY,NX)+TCO2P(L,NY,NX)+TCO2S(L,NY,NX)
      UDICD(NY,NX)=UDICD(NY,NX)-12.0*TBCO2(L,NY,NX)
      TXCO2(NY,NX)=TXCO2(NY,NX)+12.0*TBCO2(L,NY,NX)
      XONET(NY,NX)=XONET(NY,NX)+OIB
      OXYGIN=OXYGIN+OIB
      OOB=RUPOXO(L,NY,NX)+TUPOXP(L,NY,NX)+TUPOXS(L,NY,NX)
     2+ROGOX(L,NY,NX)+RC4OX(L,NY,NX)*2.667
      OXYGOU=OXYGOU+OOB
      UOXYG(NY,NX)=UOXYG(NY,NX)+OIB
      HOXYG(NY,NX)=HOXYG(NY,NX)+OIB
      H2GIN=H2GIN+HGB
      HOB=RH2GO(L,NY,NX)+TUPHGS(L,NY,NX)
      H2GOU=H2GOU+HOB
      ZN2GIN=ZN2GIN+ZGB+Z2B+ZHB
C     UN2GG(NY,NX)=UN2GG(NY,NX)+ZGB
C     HN2GG(NY,NX)=HN2GG(NY,NX)+ZGB
      UN2OG(NY,NX)=UN2OG(NY,NX)+Z2B
      HN2OG(NY,NX)=HN2OG(NY,NX)+Z2B
      UNH3G(NY,NX)=UNH3G(NY,NX)+ZHB
      HNH3G(NY,NX)=HNH3G(NY,NX)+ZHB
      UH2GG(NY,NX)=UH2GG(NY,NX)+HGB
C     IF(NY.EQ.5)THEN
C     WRITE(*,6645)'PLT',I,J,NFZ,NX,NY,L,LG,LL
C    2,XCNET(NY,NX),XONET(NY,NX),XHNET(NY,NX)
C    2,TXCO2(NY,NX),TBCO2(L,NY,NX)
C    2,HCH4G(NY,NX),CHB,TCHFLA(L,NY,NX),XCHBBL(L,NY,NX)
C    2,HOXYG(NY,NX),OIB
C    3,XOXBBL(L,NY,NX),TUPOXP(L,NY,NX),TUPOXS(L,NY,NX)
C    4,TOXFLA(L,NY,NX),OXYG(L,NY,NX),SOXYL(L,NY,NX)
C    4,HCO2G(NY,NX),CIB,TCOFLA(L,NY,NX),XCOBBL(L,NY,NX)
C    4,TRCO2(L,NY,NX)
C    2,UN2OG(NY,NX),ZGI,XN2BBL(L,NY,NX)
C    5,TN2FLA(L,NY,NX),TNHFLA(L,NY,NX),THGFLA(L,NY,NX)
C    2,UN2GG(NY,NX),ZGI,XNGBBL(L,NY,NX)
C    5,TN2FLA(L,NY,NX),TNHFLA(L,NY,NX),THGFLA(L,NY,NX)
C    6,CH4G(LL,NY,NX)
6645  FORMAT(A8,8I4,30E12.4)
C      ENDIF
C
C     SOIL SOLUTE CHANGES FROM SOIL TRANSFORMATION, ROOT UPTAKE,
C     ION REACTIONS
C
C     X*S=soil solute exchange from nitro.f, solute.f
C     TUP*=root solute uptake from uptake.f
C     TIONOU=total solute salt subsurface change
C
      SNM=(2.0*(XNH4S(L,NY,NX)+XNH4B(L,NY,NX)-TUPNH4(L,NY,NX)
     2-TUPNHB(L,NY,NX)-XN2GS(L,NY,NX))-TUPN3S(L,NY,NX)-TUPN3B(L,NY,NX)
     3+XNO3S(L,NY,NX)+XNO3B(L,NY,NX)-TUPNO3(L,NY,NX)-TUPNOB(L,NY,NX)
     5+XNO2S(L,NY,NX)+XNO2B(L,NY,NX))/14.0
      SPM=(2.0*(XH1PS(L,NY,NX)+XH1BS(L,NY,NX)-TUPH1P(L,NY,NX)
     2-TUPH1B(L,NY,NX))+3.0*(XH2PS(L,NY,NX)+XH2BS(L,NY,NX)
     3-TUPH2P(L,NY,NX)-TUPH2B(L,NY,NX)))/31.0
      SSB=TRH2O(L,NY,NX)+TBCO2(L,NY,NX)+XZHYS(L,NY,NX)
     2+TBION(L,NY,NX)
      TIONOU=TIONOU-SSB
C     UIONOU(NY,NX)=UIONOU(NY,NX)-SSB
C     WRITE(20,3339)'SSB',I,J,L,SSB,TRH2O(L,NY,NX)
C    2,TBCO2(L,NY,NX),XZHYS(L,NY,NX),TBION(L,NY,NX)
C
C     SOIL GAS AND SOLUTE EXCHANGE WITHIN GRID CELL ADDED TO ECOSYSTEM
C     TOTALS FOR CALCULATING COMPETITION CONSTRAINTS ON MICROBIAL
C     AND ROOT POPULATIONS
C
C     R*X=total substrate demand from all substrate unlimited by
C        uptake calculated in nitro.f, uptake.f for use in next
C        nitro.f, uptake.f
C     substrate code: OXY=O2, NH4=NH4 non-band, NB4=NH4 band
C        NO3=NO3 non-band, NB3=NO3 band, PO4=H2PO4 non-band
C        POB=H2PO4 band,P14=HPO4 non-band, P1B=HPO4 band, OQC=DOC
C        oxidation, OQA=acetate oxidation
C
      DO 7990 K=0,5
      DO 7980 N=1,7
      ROXYX(L,NY,NX)=ROXYX(L,NY,NX)+ROXYS(N,K,L,NY,NX)
      RNH4X(L,NY,NX)=RNH4X(L,NY,NX)+RVMX4(N,K,L,NY,NX)
     2+RINHO(N,K,L,NY,NX)
      RNO3X(L,NY,NX)=RNO3X(L,NY,NX)+RVMX3(N,K,L,NY,NX)
     2+RINOO(N,K,L,NY,NX)
      RNO2X(L,NY,NX)=RNO2X(L,NY,NX)+RVMX2(N,K,L,NY,NX)
      RN2OX(L,NY,NX)=RN2OX(L,NY,NX)+RVMX1(N,K,L,NY,NX)
      RPO4X(L,NY,NX)=RPO4X(L,NY,NX)+RIPOO(N,K,L,NY,NX)
      RP14X(L,NY,NX)=RP14X(L,NY,NX)+RIPO1(N,K,L,NY,NX)
      RNHBX(L,NY,NX)=RNHBX(L,NY,NX)+RVMB4(N,K,L,NY,NX)
     2+RINHB(N,K,L,NY,NX)
      RN3BX(L,NY,NX)=RN3BX(L,NY,NX)+RVMB3(N,K,L,NY,NX)
     2+RINOB(N,K,L,NY,NX)
      RN2BX(L,NY,NX)=RN2BX(L,NY,NX)+RVMB2(N,K,L,NY,NX)
      RPOBX(L,NY,NX)=RPOBX(L,NY,NX)+RIPBO(N,K,L,NY,NX)
      RP1BX(L,NY,NX)=RP1BX(L,NY,NX)+RIPB1(N,K,L,NY,NX)
      IF(K.LE.4)THEN
      ROQCX(K,L,NY,NX)=ROQCX(K,L,NY,NX)+ROQCS(N,K,L,NY,NX)
      ROQAX(K,L,NY,NX)=ROQAX(K,L,NY,NX)+ROQAS(N,K,L,NY,NX)
      ENDIF
7980  CONTINUE
7990  CONTINUE
      RNO2X(L,NY,NX)=RNO2X(L,NY,NX)+RVMXC(L,NY,NX)
      RN2BX(L,NY,NX)=RN2BX(L,NY,NX)+RVMBC(L,NY,NX)
C
C     GRID CELL VARIABLES NEEDED FOR WATER, C, N, P, O, SOLUTE AND
C     ENERGY BALANCES INCLUDING SUM OF ALL CURRENT STATE VARIABLES,
C     CUMULATIVE SUMS OF ALL ADDITIONS AND REMOVALS SINCE START OF RUN
C
C     VOLW,VOLI,VOLWH,VOLIH=water,ice content in micropores,macropores
C     VOLWSO,UVOLW=total landscape, grid cell water content
C     VHCP=soil heat capacity
C     TKS=soil temperature
C     HEATSO=total landscape heat content
C     SAND,SILT,CLAY=soil sand, silt, clay content
C     TSEDSO=total sediment content
C     TL*,U*=total landscape, grid cell gas, solute content
C     *G,*S,*SH=soil gaseous gas, soil aqueous gas in
C        micropores,macropores
C     T*P=total root gas content
C     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
C             :*ZN3*=NH3,*H2G*=H2
C     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
C        :OC=DOC,ON=DON,OP=DOP,OA=acetate
C        :N4=NH4,N3=NH3,NO=NO3,NX=NO2,H1P=HPO4,H2P=H2PO4 in non-band
C
C     IF(J.EQ.24)THEN
      WS=VOLW(L,NY,NX)+VOLV(L,NY,NX)+VOLWH(L,NY,NX)
     2+(VOLI(L,NY,NX)+VOLIH(L,NY,NX))*DENSI
      VOLWSO=VOLWSO+WS
      VOLISO=VOLISO+VOLI(L,NY,NX)+VOLIH(L,NY,NX)
      UVOLW(NY,NX)=UVOLW(NY,NX)+WS
C    2-WP(L,NY,NX)*VOLX(L,NY,NX)
      HEATSO=HEATSO+VHCP(L,NY,NX)*TKS(L,NY,NX)
      SD=SAND(L,NY,NX)+SILT(L,NY,NX)+CLAY(L,NY,NX)
      TSEDSO=TSEDSO+SD
      CS=CO2G(L,NY,NX)+CO2S(L,NY,NX)+CO2SH(L,NY,NX)+TLCO2P(L,NY,NX)
     2+CH4G(L,NY,NX)+CH4S(L,NY,NX)+CH4SH(L,NY,NX)+TLCH4P(L,NY,NX)
      TLCO2G=TLCO2G+CS
      UCO2S(NY,NX)=UCO2S(NY,NX)+CS
      HS=H2GG(L,NY,NX)+H2GS(L,NY,NX)+H2GSH(L,NY,NX)+TLH2GP(L,NY,NX)
      TLH2G=TLH2G+HS
C     IF(NX.EQ.1.AND.NY.EQ.1)THEN
C     WRITE(*,8642)'TLCO2G',I,J,NFZ,L
C    2,TLCO2G,CS,CO2G(L,NY,NX),CO2S(L,NY,NX)
C    2,CO2SH(L,NY,NX),TLCO2P(L,NY,NX),CH4G(L,NY,NX),CH4S(L,NY,NX)
C    3,CH4SH(L,NY,NX),TLCH4P(L,NY,NX)
8642  FORMAT(A8,4I4,20F16.6)
C     ENDIF
      OS=OXYG(L,NY,NX)+OXYS(L,NY,NX)+OXYSH(L,NY,NX)+TLOXYP(L,NY,NX)
      OXYGSO=OXYGSO+OS
      ZGS=Z2GG(L,NY,NX)+Z2GS(L,NY,NX)+Z2GSH(L,NY,NX)+TLN2OP(L,NY,NX)
     2+Z2OG(L,NY,NX)+Z2OS(L,NY,NX)+Z2OSH(L,NY,NX)+TLNH3P(L,NY,NX)
     3+ZNH3G(L,NY,NX)
      TLN2G=TLN2G+ZGS
      Z4S=ZNH4S(L,NY,NX)+ZNH4SH(L,NY,NX)+ZNH4B(L,NY,NX)
     2+ZNH4BH(L,NY,NX)+ZNH3S(L,NY,NX)+ZNH3SH(L,NY,NX)
     3+ZNH3B(L,NY,NX)+ZNH3BH(L,NY,NX)
      Z4X=14.0*(XN4(L,NY,NX)+XNB(L,NY,NX))
      Z4F=14.0*(ZNH4FA(L,NY,NX)+ZNHUFA(L,NY,NX)+ZNH3FA(L,NY,NX)
     2+ZNH4FB(L,NY,NX)+ZNHUFB(L,NY,NX)+ZNH3FB(L,NY,NX))
      TLNH4=TLNH4+Z4S+Z4X+Z4F
      UNH4(NY,NX)=UNH4(NY,NX)+Z4S+Z4X
C     IF(I.EQ.168)THEN
C     WRITE(*,5455)'TLN2GL',I,J,NX,NY,L,TLN2G
C    2,ZG,Z2GG(L,NY,NX),Z2GS(L,NY,NX),Z2GSH(L,NY,NX),TLN2OP(L,NY,NX)
C    2,Z2OG(L,NY,NX),Z2OS(L,NY,NX),Z2OSH(L,NY,NX),TLNH3P(L,NY,NX)
C    3,ZNH3G(L,NY,NX)
C     WRITE(*,5455)'TLNH4L',I,J,NX,NY,L,TLNH4,UNH4(NY,NX)
C    2,Z4S,Z4X,Z4F,XN4(L,NY,NX)
C    2,XNB(L,NY,NX),ZNH4S(L,NY,NX),ZNH4SH(L,NY,NX)
C    3,ZNH4B(L,NY,NX),ZNH4BH(L,NY,NX),ZNH3S(L,NY,NX),ZNH3SH(L,NY,NX)
C    4,ZNH3B(L,NY,NX),ZNH3BH(L,NY,NX),TN4FHB(L,NY,NX),XN4FXB(L,NY,NX)
5455  FORMAT(A8,5I4,30E12.4)
C     ENDIF
      ZOS=ZNO3S(L,NY,NX)+ZNO3SH(L,NY,NX)+ZNO3B(L,NY,NX)
     2+ZNO3BH(L,NY,NX)+ZNO2S(L,NY,NX)+ZNO2SH(L,NY,NX)
     3+ZNO2B(L,NY,NX)+ZNO2BH(L,NY,NX)
      ZOF=14.0*(ZNO3FA(L,NY,NX)+ZNO3FA(L,NY,NX))
      TLNO3=TLNO3+ZOS+ZOF
      UNO3(NY,NX)=UNO3(NY,NX)+ZOS
      P4S=H2PO4(L,NY,NX)+H2PO4H(L,NY,NX)+H2POB(L,NY,NX)
     2+H2POBH(L,NY,NX)+H1PO4(L,NY,NX)+H1PO4H(L,NY,NX)
     3+H1POB(L,NY,NX)+H1POBH(L,NY,NX)
      P4X=31.0*(XH1P(L,NY,NX)+XH2P(L,NY,NX)
     4+XH1PB(L,NY,NX)+XH2PB(L,NY,NX))
      P4P=31.0*(PALPO(L,NY,NX)+PFEPO(L,NY,NX)+PCAPD(L,NY,NX)
     6+PALPB(L,NY,NX)+PFEPB(L,NY,NX)+PCPDB(L,NY,NX))
     7+62.0*(PCAPM(L,NY,NX)+PCPMB(L,NY,NX))
     8+93.0*(PCAPH(L,NY,NX)+PCPHB(L,NY,NX))
      TLPO4=TLPO4+P4S+P4X+P4P
      UPO4(NY,NX)=UPO4(NY,NX)+P4S
      UPX4(NY,NX)=UPX4(NY,NX)+P4X
      UPP4(NY,NX)=UPP4(NY,NX)+P4P
C     IF((I/10)*10.EQ.I.AND.J.EQ.24.AND.NFZ.EQ.1)THEN
C     WRITE(*,2233)'TLPO4',I,J,NFZ,L,UPP4(NY,NX)
C    2,P4S,P4X,P4P,TLPO4,TSEDER(NY,NX)
C    3,PALPO(L,NY,NX),PFEPO(L,NY,NX),PCAPD(L,NY,NX)
C    6,PALPB(L,NY,NX),PFEPB(L,NY,NX),PCPDB(L,NY,NX)
C    7,PCAPM(L,NY,NX),PCPMB(L,NY,NX)
C    8,PCAPH(L,NY,NX),PCPHB(L,NY,NX)
C    3,TH1PEB(NY,NX),TH2PEB(NY,NX),XH1PB(L,NY,NX),XH2PB(L,NY,NX)
2233  FORMAT(A8,4I4,30E12.4)
C     ENDIF
C     ENDIF
C
C     TOTAL SOC,SON,SOP
C
C     OMC,OMN,OMP=microbial C,N,P
C     OQC,OQN,OPQ,OQCH,OQNH,OQPH=DOC,DON,DOP in micropores,macropores
C     ORC,ORN,ORP=microbial residue C,N,P
C     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
C     OSC,OSA,OSN,OSP=SOC,colonized SOC,SON,SOP
C     TOMT,TONT,TOPT=total microbial C,N,P
C     DORGC=change in total SOC
C     ORGC,ORGCC,ORGN,ORGNC=total SOC,charcoal,SON,charcoal N
C     TLRSDC,TLRSDN,TLRSDP=total landscape litter C,N,P
C     TLORGC,TLORGN,TLORGP=total landscape humus C,N,P
C     URSDC,URSDN,URSDP=total grid cell litter C,N,P
C     UORGC,UORGN,UORGP=total grid cell humus C,N,P
C
      DC=0.0
      DN=0.0
      DP=0.0
      OC=0.0
      ON=0.0
      OP=0.0
      DCC=0.0
      DNC=0.0
      DPC=0.0
      OCC=0.0
      ONC=0.0
      OPC=0.0
      OMCL(L,NY,NX)=0.0
      OMNL(L,NY,NX)=0.0
      DO 7970 K=0,5
      IF(K.LE.2)THEN
      DO 7960 N=1,7
      DO 7960 M=1,3
      DC=DC+OMC(M,N,K,L,NY,NX)
      DN=DN+OMN(M,N,K,L,NY,NX)
      DP=DP+OMP(M,N,K,L,NY,NX)
      TOMT(NY,NX)=TOMT(NY,NX)+OMC(M,N,K,L,NY,NX)
      TONT(NY,NX)=TONT(NY,NX)+OMN(M,N,K,L,NY,NX)
      TOPT(NY,NX)=TOPT(NY,NX)+OMP(M,N,K,L,NY,NX)
      OMCL(L,NY,NX)=OMCL(L,NY,NX)+OMC(M,N,K,L,NY,NX)
      OMNL(L,NY,NX)=OMNL(L,NY,NX)+OMN(M,N,K,L,NY,NX)
7960  CONTINUE
      ELSE
      DO 7950 N=1,7
      DO 7950 M=1,3
      OC=OC+OMC(M,N,K,L,NY,NX)
      ON=ON+OMN(M,N,K,L,NY,NX)
      OP=OP+OMP(M,N,K,L,NY,NX)
      TOMT(NY,NX)=TOMT(NY,NX)+OMC(M,N,K,L,NY,NX)
      TONT(NY,NX)=TONT(NY,NX)+OMN(M,N,K,L,NY,NX)
      TOPT(NY,NX)=TOPT(NY,NX)+OMP(M,N,K,L,NY,NX)
      OMCL(L,NY,NX)=OMCL(L,NY,NX)+OMC(M,N,K,L,NY,NX)
      OMNL(L,NY,NX)=OMNL(L,NY,NX)+OMN(M,N,K,L,NY,NX)
7950  CONTINUE
      ENDIF
7970  CONTINUE
      DO 7900 K=0,4
      IF(K.LE.2)THEN
      DO 7940 M=1,2
      DC=DC+ORC(M,K,L,NY,NX)
      DN=DN+ORN(M,K,L,NY,NX)
      DP=DP+ORP(M,K,L,NY,NX)
7940  CONTINUE
      DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
     2+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      DP=DP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
      DO 7930 M=1,5
      IF(M.LE.4)THEN
      DC=DC+OSC(M,K,L,NY,NX)
      DN=DN+OSN(M,K,L,NY,NX)
      DP=DP+OSP(M,K,L,NY,NX)
      ELSE
      DCC=DCC+OSC(M,K,L,NY,NX)
      DNC=DNC+OSN(M,K,L,NY,NX)
      DPC=DPC+OSP(M,K,L,NY,NX)
      ENDIF
7930  CONTINUE
      ELSE
      DO 7920 M=1,2
      OC=OC+ORC(M,K,L,NY,NX)
      ON=ON+ORN(M,K,L,NY,NX)
      OP=OP+ORP(M,K,L,NY,NX)
7920  CONTINUE
      OC=OC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
     2+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      ON=ON+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      OP=OP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
      DO 7910 M=1,5
      IF(M.LE.4)THEN
      OC=OC+OSC(M,K,L,NY,NX)
      ON=ON+OSN(M,K,L,NY,NX)
      OP=OP+OSP(M,K,L,NY,NX)
      ELSE
      OCC=OCC+OSC(M,K,L,NY,NX)
      ONC=ONC+OSN(M,K,L,NY,NX)
      OPC=OPC+OSP(M,K,L,NY,NX)
      ENDIF
7910  CONTINUE
      ENDIF
7900  CONTINUE
      ORGC(L,NY,NX)=DC+OC
      ORGN(L,NY,NX)=DN+ON
      ORGCC(L,NY,NX)=DCC+OCC
      ORGNC(L,NY,NX)=DNC+ONC
      ORGR(L,NY,NX)=DC
      IF(IERSNG.EQ.2.OR.IERSNG.EQ.3)THEN
      DORGC(L,NY,NX)=ORGCX(L,NY,NX)-ORGC(L,NY,NX)-ORGCC(L,NY,NX)
      IF(L.EQ.NU(NY,NX))THEN
      DORGC(L,NY,NX)=DORGC(L,NY,NX)+DORGE(NY,NX)
      ENDIF
      ELSE
      DORGC(L,NY,NX)=0.0
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.24)THEN
C     DO 4344 K=0,4
C     WRITE(*,4343)'ORGC',I,J,NFZ,NX,NY,L,K
C    2,ORGC(L,NY,NX),ORGCC(L,NY,NX),DC,OC,DCC,OCC
C    2,((OMC(M,N,K,L,NY,NX),M=1,3),N=1,7)
C    3,(ORC(M,K,L,NY,NX),M=1,2)
C    4,OQC(K,L,NY,NX),OQCH(K,L,NY,NX),OHC(K,L,NY,NX)
C    2,OQA(K,L,NY,NX),OQAH(K,L,NY,NX),OHA(K,L,NY,NX)
C    5,(OSC(M,K,L,NY,NX),M=1,5)
C     WRITE(*,4343)'ORGN',I,J,NFZ,NX,NY,L,K,DN,DNC,ON,ONC
C    2,ORGN(L,NY,NX),ORGNC(L,NY,NX)
C    2,((OMN(M,N,K,L,NY,NX),M=1,3),N=1,7)
C    3,(ORN(M,K,L,NY,NX),M=1,2)
C    4,(OSN(M,K,L,NY,NX),M=1,5)
C    4,OQN(K,L,NY,NX),OQNH(K,L,NY,NX),OHN(K,L,NY,NX)
4343  FORMAT(A8,7I4,120E14.6)
4344  CONTINUE
C     ENDIF
      TLRSDC=TLRSDC+DC+DCC
      URSDC(NY,NX)=URSDC(NY,NX)+DC+DCC
      TLRSDN=TLRSDN+DN+DNC
      URSDN(NY,NX)=URSDN(NY,NX)+DN+DNC
      TLRSDP=TLRSDP+DP+DPC
      URSDP(NY,NX)=URSDP(NY,NX)+DP+DPC
      TLORGC=TLORGC+OC+OCC
      UORGC(NY,NX)=UORGC(NY,NX)+OC+OCC
      TLORGN=TLORGN+ON+ONC
      UORGN(NY,NX)=UORGN(NY,NX)+ON+ONC
      TLORGP=TLORGP+OP+OPC
      UORGP(NY,NX)=UORGP(NY,NX)+OP+OPC
      TSEDSO=TSEDSO+(DC+OC)*1.0E-06
C     IF(L.EQ.NU(NY,NX))THEN
C     WRITE(*,2234)'TLORGP',I,J,NX,NY,L,TLRSDC,TLORGC,DC,OC
2234  FORMAT(A8,5I4,20F16.6)
C     ENDIF
C
C     TOTAL SALT IONS
C
C     Z*,Z*H=soil salt contents in micropores,macrpores
C     TR*=total solute transformation from solute.f
C     T*FLS,T*FHS=net convective+diffusive solute flux through
C        micropores
C       ,macropores in non-band from trnsfrs.f
C     T*FLB,T*FHB=net convective+diffusive solute flux through
C        micropores,macropores in band from trnsfrs.f
C     R*FLU,R*FBU=subsurface solute flux in non-band,band
C     X*FXS,X*FXB=convective+diffusive solute flux between macropores
C        and micropores in non-band,band from trnsfr.f
C     salt code:*HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+
C        :*MG*=Mg2+,*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-
C        :*CO3*=CO32-,*HCO3*=HCO3-,*CO2*=CO2,*ALO1*=AlOH2-
C        :*ALOH2=AlOH2-,*ALOH3*=AlOH3,*ALOH4*=AlOH4+,*ALS*=AlSO4+
C        :*FEO1*=FeOH2-,*FEOH2=F3OH2-,*FEOH3*=FeOH3,*FEOH4*=FeOH4+
C        :*FES*=FeSO4+,*CAO*=CaOH,*CAC*=CaCO3,*CAH*=CaHCO3-
C        :*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3,*MHG*=MgHCO3-
C        :*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
C     phosphorus code:*H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-
C        :*F2P*=F1H2PO4-,*C0P*=CaPO4,*C1P*=CaHPO4,*C2P*=CaH2PO4+
C        :*M1P*=MgHPO4,*COO*=COOH-
C        :*1=non-band,*B=band
C     TION=total landscape ion content
C
      IF(ISALTG.NE.0)THEN
C
C     MICROPORES
C
      ZAL(L,NY,NX)=ZAL(L,NY,NX)+TRAL(L,NY,NX)+TALFLS(L,NY,NX)
     2+RALFLU(L,NY,NX)+XALFXS(L,NY,NX)
      ZFE(L,NY,NX)=ZFE(L,NY,NX)+TRFE(L,NY,NX)+TFEFLS(L,NY,NX)
     2+RFEFLU(L,NY,NX)+XFEFXS(L,NY,NX)
      ZHY(L,NY,NX)=ZHY(L,NY,NX)+TRHY(L,NY,NX)+THYFLS(L,NY,NX)
     2+RHYFLU(L,NY,NX)+XHYFXS(L,NY,NX)+XZHYS(L,NY,NX)
      ZCA(L,NY,NX)=ZCA(L,NY,NX)+TRCA(L,NY,NX)+TCAFLS(L,NY,NX)
     2+RCAFLU(L,NY,NX)+XCAFXS(L,NY,NX)
      ZMG(L,NY,NX)=ZMG(L,NY,NX)+TRMG(L,NY,NX)+TMGFLS(L,NY,NX)
     2+RMGFLU(L,NY,NX)+XMGFXS(L,NY,NX)
      ZNA(L,NY,NX)=ZNA(L,NY,NX)+TRNA(L,NY,NX)+TNAFLS(L,NY,NX)
     2+RNAFLU(L,NY,NX)+XNAFXS(L,NY,NX)
      ZKA(L,NY,NX)=ZKA(L,NY,NX)+TRKA(L,NY,NX)+TKAFLS(L,NY,NX)
     2+RKAFLU(L,NY,NX)+XKAFXS(L,NY,NX)
      ZOH(L,NY,NX)=ZOH(L,NY,NX)+TROH(L,NY,NX)+TOHFLS(L,NY,NX)
     2+ROHFLU(L,NY,NX)+XOHFXS(L,NY,NX)
C     IF(L.EQ.1)THEN
C     WRITE(*,5545)'ZOH',I,J,NX,NY,L
C    2,ZOH(L,NY,NX),TROH(L,NY,NX)
C    2,TOHFLS(L,NY,NX),ROHFLU(L,NY,NX),XOHFXS(L,NY,NX)
C    2,ZHY(L,NY,NX),TRHY(L,NY,NX),THYFLS(L,NY,NX)
C    2,RHYFLU(L,NY,NX),XHYFXS(L,NY,NX),XZHYS(L,NY,NX)
C    3,XHY(L,NY,NX),TRXHY(L,NY,NX)
C    4,ZOH(L,NY,NX)/VOLW(L,NY,NX),ZHY(L,NY,NX)/VOLW(L,NY,NX)
C     WRITE(*,5545)'ZAL',I,J,NX,NY,L,ZAL(L,NY,NX),TRAL(L,NY,NX)
C    2,TALFLS(L,NY,NX),RALFLU(L,NY,NX),XALFXS(L,NY,NX)
C     ENDIF
      ZSO4(L,NY,NX)=ZSO4(L,NY,NX)+TRSO4(L,NY,NX)+TSOFLS(L,NY,NX)
     2+RSOFLU(L,NY,NX)+XSOFXS(L,NY,NX)
      ZCL(L,NY,NX)=ZCL(L,NY,NX)+TCLFLS(L,NY,NX)+RCLFLU(L,NY,NX)
     2+XCLFXS(L,NY,NX)
      ZCO3(L,NY,NX)=ZCO3(L,NY,NX)+TRCO3(L,NY,NX)+TC3FLS(L,NY,NX)
     2+XC3FXS(L,NY,NX)
      ZHCO3(L,NY,NX)=ZHCO3(L,NY,NX)+TRHCO(L,NY,NX)+THCFLS(L,NY,NX)
     2+XHCFXS(L,NY,NX)
      ZALOH1(L,NY,NX)=ZALOH1(L,NY,NX)+TRAL1(L,NY,NX)+TAL1FS(L,NY,NX)
     2+XAL1XS(L,NY,NX)
      ZALOH2(L,NY,NX)=ZALOH2(L,NY,NX)+TRAL2(L,NY,NX)+TAL2FS(L,NY,NX)
     2+XAL2XS(L,NY,NX)-TRXAL2(L,NY,NX)
      ZALOH3(L,NY,NX)=ZALOH3(L,NY,NX)+TRAL3(L,NY,NX)+TAL3FS(L,NY,NX)
     2+XAL3XS(L,NY,NX)
      ZALOH4(L,NY,NX)=ZALOH4(L,NY,NX)+TRAL4(L,NY,NX)+TAL4FS(L,NY,NX)
     2+XAL4XS(L,NY,NX)
      ZALS(L,NY,NX)=ZALS(L,NY,NX)+TRALS(L,NY,NX)+TALSFS(L,NY,NX)
     2+XALSXS(L,NY,NX)
      ZFEOH1(L,NY,NX)=ZFEOH1(L,NY,NX)+TRFE1(L,NY,NX)+TFE1FS(L,NY,NX)
     2+XFE1XS(L,NY,NX)
      ZFEOH2(L,NY,NX)=ZFEOH2(L,NY,NX)+TRFE2(L,NY,NX)+TFE2FS(L,NY,NX)
     2+XFE2XS(L,NY,NX)-TRXFE2(L,NY,NX)
      ZFEOH3(L,NY,NX)=ZFEOH3(L,NY,NX)+TRFE3(L,NY,NX)+TFE3FS(L,NY,NX)
     2+XFE3XS(L,NY,NX)
      ZFEOH4(L,NY,NX)=ZFEOH4(L,NY,NX)+TRFE4(L,NY,NX)+TFE4FS(L,NY,NX)
     2+XFE4XS(L,NY,NX)
      ZFES(L,NY,NX)=ZFES(L,NY,NX)+TRFES(L,NY,NX)+TFESFS(L,NY,NX)
     2+XFESXS(L,NY,NX)
      ZCAO(L,NY,NX)=ZCAO(L,NY,NX)+TRCAO(L,NY,NX)+TCAOFS(L,NY,NX)
     2+XCAOXS(L,NY,NX)
      ZCAC(L,NY,NX)=ZCAC(L,NY,NX)+TRCAC(L,NY,NX)+TCACFS(L,NY,NX)
     2+XCACXS(L,NY,NX)
      ZCAH(L,NY,NX)=ZCAH(L,NY,NX)+TRCAH(L,NY,NX)+TCAHFS(L,NY,NX)
     2+XCAHXS(L,NY,NX)
      ZCAS(L,NY,NX)=ZCAS(L,NY,NX)+TRCAS(L,NY,NX)+TCASFS(L,NY,NX)
     2+XCASXS(L,NY,NX)
      ZMGO(L,NY,NX)=ZMGO(L,NY,NX)+TRMGO(L,NY,NX)+TMGOFS(L,NY,NX)
     2+XMGOXS(L,NY,NX)
      ZMGC(L,NY,NX)=ZMGC(L,NY,NX)+TRMGC(L,NY,NX)+TMGCFS(L,NY,NX)
     2+XMGCXS(L,NY,NX)
      ZMGH(L,NY,NX)=ZMGH(L,NY,NX)+TRMGH(L,NY,NX)+TMGHFS(L,NY,NX)
     2+XMGHXS(L,NY,NX)
      ZMGS(L,NY,NX)=ZMGS(L,NY,NX)+TRMGS(L,NY,NX)+TMGSFS(L,NY,NX)
     2+XMGSXS(L,NY,NX)
      ZNAC(L,NY,NX)=ZNAC(L,NY,NX)+TRNAC(L,NY,NX)+TNACFS(L,NY,NX)
     2+XNACXS(L,NY,NX)
      ZNAS(L,NY,NX)=ZNAS(L,NY,NX)+TRNAS(L,NY,NX)+TNASFS(L,NY,NX)
     2+XNASXS(L,NY,NX)
      ZKAS(L,NY,NX)=ZKAS(L,NY,NX)+TRKAS(L,NY,NX)+TKASFS(L,NY,NX)
     2+XKASXS(L,NY,NX)
      H0PO4(L,NY,NX)=H0PO4(L,NY,NX)+TRH0P(L,NY,NX)+TH0PFS(L,NY,NX)
     2+XH0PXS(L,NY,NX)
      H3PO4(L,NY,NX)=H3PO4(L,NY,NX)+TRH3P(L,NY,NX)+TH3PFS(L,NY,NX)
     2+XH3PXS(L,NY,NX)
      ZFE1P(L,NY,NX)=ZFE1P(L,NY,NX)+TRF1P(L,NY,NX)+TF1PFS(L,NY,NX)
     2+XF1PXS(L,NY,NX)
      ZFE2P(L,NY,NX)=ZFE2P(L,NY,NX)+TRF2P(L,NY,NX)+TF2PFS(L,NY,NX)
     2+XF2PXS(L,NY,NX)
      ZCA0P(L,NY,NX)=ZCA0P(L,NY,NX)+TRC0P(L,NY,NX)+TC0PFS(L,NY,NX)
     2+XC0PXS(L,NY,NX)
      ZCA1P(L,NY,NX)=ZCA1P(L,NY,NX)+TRC1P(L,NY,NX)+TC1PFS(L,NY,NX)
     2+XC1PXS(L,NY,NX)
      ZCA2P(L,NY,NX)=ZCA2P(L,NY,NX)+TRC2P(L,NY,NX)+TC2PFS(L,NY,NX)
     2+XC2PXS(L,NY,NX)
      ZMG1P(L,NY,NX)=ZMG1P(L,NY,NX)+TRM1P(L,NY,NX)+TM1PFS(L,NY,NX)
     2+XM1PXS(L,NY,NX)
      H0POB(L,NY,NX)=H0POB(L,NY,NX)+TRH0B(L,NY,NX)+TH0BFB(L,NY,NX)
     2+XH0BXB(L,NY,NX)
      H3POB(L,NY,NX)=H3POB(L,NY,NX)+TRH3B(L,NY,NX)+TH3BFB(L,NY,NX)
     2+XH3BXB(L,NY,NX)
      ZFE1PB(L,NY,NX)=ZFE1PB(L,NY,NX)+TRF1B(L,NY,NX)+TF1BFB(L,NY,NX)
     2+XF1BXB(L,NY,NX)
      ZFE2PB(L,NY,NX)=ZFE2PB(L,NY,NX)+TRF2B(L,NY,NX)+TF2BFB(L,NY,NX)
     2+XF2BXB(L,NY,NX)
      ZCA0PB(L,NY,NX)=ZCA0PB(L,NY,NX)+TRC0B(L,NY,NX)+TC0BFB(L,NY,NX)
     2+XC0BXB(L,NY,NX)
      ZCA1PB(L,NY,NX)=ZCA1PB(L,NY,NX)+TRC1B(L,NY,NX)+TC1BFB(L,NY,NX)
     2+XC1BXB(L,NY,NX)
      ZCA2PB(L,NY,NX)=ZCA2PB(L,NY,NX)+TRC2B(L,NY,NX)+TC2BFB(L,NY,NX)
     2+XC2BXB(L,NY,NX)
      ZMG1PB(L,NY,NX)=ZMG1PB(L,NY,NX)+TRM1B(L,NY,NX)+TM1BFB(L,NY,NX)
     2+XM1BXB(L,NY,NX)
C
C     MACROPORES
C
      ZALH(L,NY,NX)=ZALH(L,NY,NX)+TALFHS(L,NY,NX)-XALFXS(L,NY,NX)
      ZFEH(L,NY,NX)=ZFEH(L,NY,NX)+TFEFHS(L,NY,NX)-XFEFXS(L,NY,NX)
      ZHYH(L,NY,NX)=ZHYH(L,NY,NX)+THYFHS(L,NY,NX)-XHYFXS(L,NY,NX)
      ZCCH(L,NY,NX)=ZCCH(L,NY,NX)+TCAFHS(L,NY,NX)-XCAFXS(L,NY,NX)
      ZMAH(L,NY,NX)=ZMAH(L,NY,NX)+TMGFHS(L,NY,NX)-XMGFXS(L,NY,NX)
      ZNAH(L,NY,NX)=ZNAH(L,NY,NX)+TNAFHS(L,NY,NX)-XNAFXS(L,NY,NX)
      ZKAH(L,NY,NX)=ZKAH(L,NY,NX)+TKAFHS(L,NY,NX)-XKAFXS(L,NY,NX)
      ZOHH(L,NY,NX)=ZOHH(L,NY,NX)+TOHFHS(L,NY,NX)-XOHFXS(L,NY,NX)
      ZSO4H(L,NY,NX)=ZSO4H(L,NY,NX)+TSOFHS(L,NY,NX)-XSOFXS(L,NY,NX)
      ZCLH(L,NY,NX)=ZCLH(L,NY,NX)+TCLFHS(L,NY,NX)-XCLFXS(L,NY,NX)
      ZCO3H(L,NY,NX)=ZCO3H(L,NY,NX)+TC3FHS(L,NY,NX)-XC3FXS(L,NY,NX)
      ZHCO3H(L,NY,NX)=ZHCO3H(L,NY,NX)+THCFHS(L,NY,NX)-XHCFXS(L,NY,NX)
      ZALO1H(L,NY,NX)=ZALO1H(L,NY,NX)+TAL1HS(L,NY,NX)-XAL1XS(L,NY,NX)
      ZALO2H(L,NY,NX)=ZALO2H(L,NY,NX)+TAL2HS(L,NY,NX)-XAL2XS(L,NY,NX)
      ZALO3H(L,NY,NX)=ZALO3H(L,NY,NX)+TAL3HS(L,NY,NX)-XAL3XS(L,NY,NX)
      ZALO4H(L,NY,NX)=ZALO4H(L,NY,NX)+TAL4HS(L,NY,NX)-XAL4XS(L,NY,NX)
      ZALSH(L,NY,NX)=ZALSH(L,NY,NX)+TALSHS(L,NY,NX)-XALSXS(L,NY,NX)
      ZFEO1H(L,NY,NX)=ZFEO1H(L,NY,NX)+TFE1HS(L,NY,NX)-XFE1XS(L,NY,NX)
      ZFEO2H(L,NY,NX)=ZFEO2H(L,NY,NX)+TFE2HS(L,NY,NX)-XFE2XS(L,NY,NX)
      ZFEO3H(L,NY,NX)=ZFEO3H(L,NY,NX)+TFE3HS(L,NY,NX)-XFE3XS(L,NY,NX)
      ZFEO4H(L,NY,NX)=ZFEO4H(L,NY,NX)+TFE4HS(L,NY,NX)-XFE4XS(L,NY,NX)
      ZFESH(L,NY,NX)=ZFESH(L,NY,NX)+TFESHS(L,NY,NX)-XFESXS(L,NY,NX)
      ZCAOH(L,NY,NX)=ZCAOH(L,NY,NX)+TCAOHS(L,NY,NX)-XCAOXS(L,NY,NX)
      ZCACH(L,NY,NX)=ZCACH(L,NY,NX)+TCACHS(L,NY,NX)-XCACXS(L,NY,NX)
      ZCAHH(L,NY,NX)=ZCAHH(L,NY,NX)+TCAHHS(L,NY,NX)-XCAHXS(L,NY,NX)
      ZCASH(L,NY,NX)=ZCASH(L,NY,NX)+TCASHS(L,NY,NX)-XCASXS(L,NY,NX)
      ZMGOH(L,NY,NX)=ZMGOH(L,NY,NX)+TMGOHS(L,NY,NX)-XMGOXS(L,NY,NX)
      ZMGCH(L,NY,NX)=ZMGCH(L,NY,NX)+TMGCHS(L,NY,NX)-XMGCXS(L,NY,NX)
      ZMGHH(L,NY,NX)=ZMGHH(L,NY,NX)+TMGHHS(L,NY,NX)-XMGHXS(L,NY,NX)
      ZMGSH(L,NY,NX)=ZMGSH(L,NY,NX)+TMGSHS(L,NY,NX)-XMGSXS(L,NY,NX)
      ZNACH(L,NY,NX)=ZNACH(L,NY,NX)+TNACHS(L,NY,NX)-XNACXS(L,NY,NX)
      ZNASH(L,NY,NX)=ZNASH(L,NY,NX)+TNASHS(L,NY,NX)-XNASXS(L,NY,NX)
      ZKASH(L,NY,NX)=ZKASH(L,NY,NX)+TKASHS(L,NY,NX)-XKASXS(L,NY,NX)
      H0PO4H(L,NY,NX)=H0PO4H(L,NY,NX)+TH0PHS(L,NY,NX)-XH0PXS(L,NY,NX)
      H3PO4H(L,NY,NX)=H3PO4H(L,NY,NX)+TH3PHS(L,NY,NX)-XH3PXS(L,NY,NX)
      ZFE1PH(L,NY,NX)=ZFE1PH(L,NY,NX)+TF1PHS(L,NY,NX)-XF1PXS(L,NY,NX)
      ZFE2PH(L,NY,NX)=ZFE2PH(L,NY,NX)+TF2PHS(L,NY,NX)-XF2PXS(L,NY,NX)
      ZCA0PH(L,NY,NX)=ZCA0PH(L,NY,NX)+TC0PHS(L,NY,NX)-XC0PXS(L,NY,NX)
      ZCA1PH(L,NY,NX)=ZCA1PH(L,NY,NX)+TC1PHS(L,NY,NX)-XC1PXS(L,NY,NX)
      ZCA2PH(L,NY,NX)=ZCA2PH(L,NY,NX)+TC2PHS(L,NY,NX)-XC2PXS(L,NY,NX)
      ZMG1PH(L,NY,NX)=ZMG1PH(L,NY,NX)+TM1PHS(L,NY,NX)-XM1PXS(L,NY,NX)
      H0POBH(L,NY,NX)=H0POBH(L,NY,NX)+TH0BHB(L,NY,NX)-XH0BXB(L,NY,NX)
      H3POBH(L,NY,NX)=H3POBH(L,NY,NX)+TH3BHB(L,NY,NX)-XH3BXB(L,NY,NX)
      ZFE1BH(L,NY,NX)=ZFE1BH(L,NY,NX)+TF1BHB(L,NY,NX)-XF1BXB(L,NY,NX)
      ZFE2BH(L,NY,NX)=ZFE2BH(L,NY,NX)+TF2BHB(L,NY,NX)-XF2BXB(L,NY,NX)
      ZCA0BH(L,NY,NX)=ZCA0BH(L,NY,NX)+TC0BHB(L,NY,NX)-XC0BXB(L,NY,NX)
      ZCA1BH(L,NY,NX)=ZCA1BH(L,NY,NX)+TC1BHB(L,NY,NX)-XC1BXB(L,NY,NX)
      ZCA2BH(L,NY,NX)=ZCA2BH(L,NY,NX)+TC2BHB(L,NY,NX)-XC2BXB(L,NY,NX)
      ZMG1BH(L,NY,NX)=ZMG1BH(L,NY,NX)+TM1BHB(L,NY,NX)-XM1BXB(L,NY,NX)
      XHY(L,NY,NX)=XHY(L,NY,NX)+TRXHY(L,NY,NX)
      XAL(L,NY,NX)=XAL(L,NY,NX)+TRXAL(L,NY,NX)
      XFE(L,NY,NX)=XFE(L,NY,NX)+TRXFE(L,NY,NX)
      XCA(L,NY,NX)=XCA(L,NY,NX)+TRXCA(L,NY,NX)
      XMG(L,NY,NX)=XMG(L,NY,NX)+TRXMG(L,NY,NX)
      XNA(L,NY,NX)=XNA(L,NY,NX)+TRXNA(L,NY,NX)
      XKA(L,NY,NX)=XKA(L,NY,NX)+TRXKA(L,NY,NX)
      XHC(L,NY,NX)=XHC(L,NY,NX)+TRXHC(L,NY,NX)
      XALO2(L,NY,NX)=XALO2(L,NY,NX)+TRXAL2(L,NY,NX)
      XFEO2(L,NY,NX)=XFEO2(L,NY,NX)+TRXFE2(L,NY,NX)
      PALOH(L,NY,NX)=PALOH(L,NY,NX)+TRALOH(L,NY,NX)
      PFEOH(L,NY,NX)=PFEOH(L,NY,NX)+TRFEOH(L,NY,NX)
      PCACO(L,NY,NX)=PCACO(L,NY,NX)+TRCACO(L,NY,NX)
      PCASO(L,NY,NX)=PCASO(L,NY,NX)+TRCASO(L,NY,NX)
      PSS=31.0*(H0PO4(L,NY,NX)+H3PO4(L,NY,NX)+ZFE1P(L,NY,NX)
     2+ZFE2P(L,NY,NX)+ZCA0P(L,NY,NX)+ZCA1P(L,NY,NX)
     3+ZCA2P(L,NY,NX)+ZMG1P(L,NY,NX)+H0POB(L,NY,NX)
     4+H3POB(L,NY,NX)+ZFE1PB(L,NY,NX)+ZFE2PB(L,NY,NX)
     5+ZCA0PB(L,NY,NX)+ZCA1PB(L,NY,NX)+ZCA2PB(L,NY,NX)
     6+ZMG1PB(L,NY,NX)+H0PO4H(L,NY,NX)+H3PO4H(L,NY,NX)
     7+ZFE1PH(L,NY,NX)+ZFE2PH(L,NY,NX)+ZCA0PH(L,NY,NX)
     8+ZCA1PH(L,NY,NX)+ZCA2PH(L,NY,NX)+ZMG1PH(L,NY,NX)
     9+H0POBH(L,NY,NX)+H3POBH(L,NY,NX)+ZFE1BH(L,NY,NX)
     1+ZFE2BH(L,NY,NX)+ZCA0BH(L,NY,NX)+ZCA1BH(L,NY,NX)
     2+ZCA2BH(L,NY,NX)+ZMG1BH(L,NY,NX))
      TLPO4=TLPO4+PSS
      SSS=ZAL(L,NY,NX)+ZFE(L,NY,NX)+ZHY(L,NY,NX)+ZCA(L,NY,NX)
     2+ZMG(L,NY,NX)+ZNA(L,NY,NX)+ZKA(L,NY,NX)+ZOH(L,NY,NX)
     3+ZSO4(L,NY,NX)+ZCL(L,NY,NX)+ZCO3(L,NY,NX)+H0PO4(L,NY,NX)
     4+H0POB(L,NY,NX)
     5+2.0*(ZHCO3(L,NY,NX)+ZALOH1(L,NY,NX)
     5+ZALS(L,NY,NX)+ZFEOH1(L,NY,NX)+ZFES(L,NY,NX)+ZCAO(L,NY,NX)
     6+ZCAC(L,NY,NX)+ZCAS(L,NY,NX)+ZMGO(L,NY,NX)+ZMGC(L,NY,NX)
     7+ZMGS(L,NY,NX)+ZNAC(L,NY,NX)+ZNAS(L,NY,NX)+ZKAS(L,NY,NX)
     8+ZCA0P(L,NY,NX)+ZCA0PB(L,NY,NX))
     9+3.0*(ZALOH2(L,NY,NX)+ZFEOH2(L,NY,NX)+ZCAH(L,NY,NX)
     1+ZMGH(L,NY,NX)+ZFE1P(L,NY,NX)+ZCA1P(L,NY,NX)+ZMG1P(L,NY,NX)
     2+ZFE1PB(L,NY,NX)+ZCA1PB(L,NY,NX)+ZMG1PB(L,NY,NX))
     3+4.0*(ZALOH3(L,NY,NX)+ZFEOH3(L,NY,NX)+H3PO4(L,NY,NX)
     3+ZFE2P(L,NY,NX)+ZCA2P(L,NY,NX)+H3POB(L,NY,NX)+ZFE2PB(L,NY,NX)
     5+ZCA2PB(L,NY,NX))
     6+5.0*(ZALOH4(L,NY,NX)+ZFEOH4(L,NY,NX))
      SSH=ZALH(L,NY,NX)+ZFEH(L,NY,NX)+ZHYH(L,NY,NX)+ZCCH(L,NY,NX)
     2+ZMAH(L,NY,NX)+ZNAH(L,NY,NX)+ZKAH(L,NY,NX)+ZOHH(L,NY,NX)
     3+ZSO4H(L,NY,NX)+ZCLH(L,NY,NX)+ZCO3H(L,NY,NX) +H0PO4H(L,NY,NX)
     4+H0POBH(L,NY,NX)
     5+2.0*(ZHCO3H(L,NY,NX)+ZALO1H(L,NY,NX)
     5+ZALSH(L,NY,NX)+ZFEO1H(L,NY,NX)+ZFESH(L,NY,NX)+ZCAOH(L,NY,NX)
     6+ZCACH(L,NY,NX)+ZCASH(L,NY,NX)+ZMGOH(L,NY,NX)+ZMGCH(L,NY,NX)
     7+ZMGSH(L,NY,NX)+ZNACH(L,NY,NX)+ZNASH(L,NY,NX)+ZKASH(L,NY,NX)
     8+ZCA0PH(L,NY,NX)+ZCA0BH(L,NY,NX))
     9+3.0*(ZALO2H(L,NY,NX)+ZFEO2H(L,NY,NX)+ZCAHH(L,NY,NX)
     1+ZMGHH(L,NY,NX)+ZFE1PH(L,NY,NX)+ZCA1PH(L,NY,NX)+ZMG1PH(L,NY,NX)
     2+ZFE1BH(L,NY,NX)+ZCA1BH(L,NY,NX)+ZMG1BH(L,NY,NX))
     3+4.0*(ZALO3H(L,NY,NX)+ZFEO3H(L,NY,NX)+H3PO4H(L,NY,NX)
     4+ZFE2PH(L,NY,NX)+ZCA2PH(L,NY,NX)+H3POBH(L,NY,NX)
     5+ZFE2BH(L,NY,NX)+ZCA2BH(L,NY,NX))
     6+5.0*(ZALO4H(L,NY,NX)+ZFEO4H(L,NY,NX))
C
C     TOTAL FERILIZER,EXCHANGEABLE CATIONS AND ANIONS, PRECIPITATES
C
      SSF=ZNH3FA(L,NY,NX)+ZNHUFA(L,NY,NX)+ZNO3FA(L,NY,NX)
     5+ZNH3FB(L,NY,NX)+ZNHUFB(L,NY,NX)+ZNO3FB(L,NY,NX)
     2+2.0*(ZNH4FA(L,NY,NX)+ZNH4FB(L,NY,NX))
      SSX=XHY(L,NY,NX)+XAL(L,NY,NX)
     2+XFE(L,NY,NX)+XCA(L,NY,NX)+XMG(L,NY,NX)
     3+XNA(L,NY,NX)+XKA(L,NY,NX)+XHC(L,NY,NX)
     4+XOH0(L,NY,NX)+XOH0B(L,NY,NX)
     5+2.0*(XN4(L,NY,NX)+XNB(L,NY,NX)
     6+XOH1(L,NY,NX)+XOH1B(L,NY,NX))
     7+3.0*(XALO2(L,NY,NX)+XFEO2(L,NY,NX)
     8+XOH2(L,NY,NX)+XOH2B(L,NY,NX)
     9+XH1P(L,NY,NX)+XH1PB(L,NY,NX))
     1+4.0*(XH2P(L,NY,NX)+XH2PB(L,NY,NX))
      SSP=2.0*(PCACO(L,NY,NX)+PCASO(L,NY,NX)
     2+PALPO(L,NY,NX)+PFEPO(L,NY,NX)
     3+PALPB(L,NY,NX)+PFEPB(L,NY,NX))
     4+3.0*(PCAPD(L,NY,NX)+PCPDB(L,NY,NX))
     5+4.0*(PALOH(L,NY,NX)+PFEOH(L,NY,NX))
     6+7.0*(PCAPM(L,NY,NX)+PCPMB(L,NY,NX))
     7+9.0*(PCAPH(L,NY,NX)+PCPHB(L,NY,NX))
      SST=SSS+SSH+SSF+SSX+SSP
      TION=TION+SST
      UION(NY,NX)=UION(NY,NX)+SST
C     IF(I.EQ.180.AND.J.EQ.12)THEN
C     WRITE(*,3341)'SSS',I,J,NX,NY,L,SST,SSS,SSH,SSF,SSX,SSP,TION
C    2 ZHY(L,NY,NX),ZAL(L,NY,NX),ZFE(L,NY,NX),ZCA(L,NY,NX)
C    2,ZMG(L,NY,NX),ZNA(L,NY,NX),ZKA(L, NY,NX),ZOH(L,NY,NX)
C    3,ZSO4(L,NY,NX),ZCL(L,NY,NX),ZCO3(L,NY,NX),H0PO4(L,NY,NX)
C    4,H0POB(L,NY,NX)
C    5,ZHCO3(L,NY,NX),ZALOH1(L,NY,NX)
C    5,ZALS(L,NY,NX),ZFEOH1(L,NY,NX),ZFES(L,NY,NX),ZCAO(L,NY,NX)
C    6,ZCAC(L,NY,NX),ZCAS(L,NY,NX),ZMGO(L,NY,NX),ZMGC(L,NY,NX)
C    7,ZMGS(L,NY,NX),ZNAC(L,NY,NX),ZNAS(L,NY,NX),ZKAS(L,NY,NX)
C    8,ZCA0P(L,NY,NX),ZCA0PB(L,NY,NX)
C    9,ZALOH2(L,NY,NX),ZFEOH2(L,NY,NX),ZCAH(L,NY,NX)
C    1,ZMGH(L,NY,NX),ZFE1P(L,NY,NX),ZCA1P(L,NY,NX),ZMG1P(L,NY,NX)
C    2,ZFE1PB(L,NY,NX),ZCA1PB(L,NY,NX),ZMG1PB(L,NY,NX)
C    3,ZALOH3(L,NY,NX),ZFEOH3(L,NY,NX),H3PO4(L,NY,NX)
C    3,ZFE2P(L,NY,NX),ZCA2P(L,NY,NX),H3POB(L,NY,NX),ZFE2PB(L,NY,NX)
C    5,ZCA2PB(L,NY,NX)
C    6,ZALOH4(L,NY,NX),ZFEOH4(L,NY,NX)
C     WRITE(*,3341)'SSX',I,J,NX,NY,L,SSX,XHY(L,NY,NX),XAL(L,NY,NX)
C    2,XFE(L,NY,NX),XCA(L,NY,NX),XMG(L,NY,NX),XNA(L,NY,NX)
C    3,XKA(L,NY,NX),XHC(L,NY,NX),XALO2(L,NY,NX),XFEO2(L,NY,NX)
C    4,PCACO(L,NY,NX),PCASO(L,NY,NX),PALOH(L,NY,NX),PFEOH(L,NY,NX)
3341  FORMAT(A8,5I4,20F14.6)
C     ENDIF
C
C     SOIL ELECTRICAL CONDUCTIVITY
C
C     VOLW=soil water content
C     E*=electrical conductivity of each solute
C     ECND=total electrical conductivity of water flux
C     solute code:HY=H+,OH=OH-,AL=Al3+,FE=Fe3+,CA=Ca2+,MG=Mg2+,NA=Na+
C            :KA=K+,CO3=CO32-,HCO3=HCO3-,SO4=SO42-,CL=Cl-,NO3S=NO3-
C
      IF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      ECHY=0.337*AMAX1(0.0,ZHY(L,NY,NX)/VOLW(L,NY,NX))
      ECOH=0.192*AMAX1(0.0,ZOH(L,NY,NX)/VOLW(L,NY,NX))
      ECAL=0.056*AMAX1(0.0,ZAL(L,NY,NX)*3.0/VOLW(L,NY,NX))
      ECFE=0.051*AMAX1(0.0,ZFE(L,NY,NX)*3.0/VOLW(L,NY,NX))
      ECCA=0.060*AMAX1(0.0,ZCA(L,NY,NX)*2.0/VOLW(L,NY,NX))
      ECMG=0.053*AMAX1(0.0,ZMG(L,NY,NX)*2.0/VOLW(L,NY,NX))
      ECNA=0.050*AMAX1(0.0,ZNA(L,NY,NX)/VOLW(L,NY,NX))
      ECKA=0.070*AMAX1(0.0,ZKA(L,NY,NX)/VOLW(L,NY,NX))
      ECCO=0.072*AMAX1(0.0,ZCO3(L,NY,NX)*2.0/VOLW(L,NY,NX))
      ECHC=0.044*AMAX1(0.0,ZHCO3(L,NY,NX)/VOLW(L,NY,NX))
      ECSO=0.080*AMAX1(0.0,ZSO4(L,NY,NX)*2.0/VOLW(L,NY,NX))
      ECCL=0.076*AMAX1(0.0,ZCL(L,NY,NX)/VOLW(L,NY,NX))
      ECNO=0.071*AMAX1(0.0,ZNO3S(L,NY,NX)/(VOLW(L,NY,NX)*14.0))
      ECND(L,NY,NX)=ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA
     2+ECCO+ECHC+ECSO+ECCL+ECNO
C     IF(I.EQ.180.AND.J.EQ.12)THEN
C     WRITE(*,5656)'ECND',I,J,L
C    2,ECND(L,NY,NX),VOLW(L,NY,NX),ECHY,ECOH,ECAL,ECFE,ECCA
C    3,ECMG,ECNA,ECKA,ECCO,ECHC,ECSO,ECCL,ECNO
5656  FORMAT(A8,3I4,30E12.4)
C     ENDIF
      ELSE
      ECND(L,NY,NX)=0.0
      ENDIF
      ENDIF
C     ENDIF
C     WRITE(20,3339)'LBN',I,J,L,TLNH4,TLNO3,TZIN,TZOU
C    2,Z4S,Z4X,Z4F,ZOS,ZOF,ZG
C    2,ZOD,ZXD,ZGD
C    3,ZGB,Z2B,ZHB
C    3,XNH4S(L,NY,NX),ZNH4S(L,NY,NX)
C    3,ZNH4SH(L,NY,NX),ZNH4B(L,NY,NX),ZNH4BH(L,NY,NX)
C    2,ZNH3S(L,NY,NX),ZNH3SH(L,NY,NX),ZNH3B(L,NY,NX),ZNH3BH(L,NY,NX)
C     WRITE(20,3339)'LBP',I,J,L,TLPO4,TPIN,TPOU,POD,PXD,PQD,PHD
C    2,P4S,P4X,P4P,PSS
C    2,XH1PS(L,NY,NX),XH2PS(L,NY,NX),H1PO4(L,NY,NX),H2PO4(L,NY,NX)
C    3,XH1P(L,NY,NX),XH2P(L,NY,NX),PALPO(L,NY,NX),PFEPO(L,NY,NX)
C    6,PCAPD(L,NY,NX),PCAPM(L,NY,NX),PCAPH(L,NY,NX)
C     WRITE(*,3339)'LBS',I,J,L,TION,TIONIN,TIONOU
C    2,SSS,SSH,SSX,SSP,SSD,SHD,SSB
3339  FORMAT(A8,3I4,80E16.8)
125   CONTINUE
C
C     SNOWPACK RELAYERING WITH CHANGES IN SNOWPACK CONTENT
C
C     VHCPW,VHCPWX=snowpack layer,minimum heat capacity
C     VOLSI=snowpack layer volume set in starts.f
C     VOLSL=current snowpack layer volume
C     VOLSSL,VOLWSL,VOLVSL,VOLISL=snow water equivalent,water,vapor,
C        ice volume in snowpack layer
C     DLYRS=snowpack layer depth
C     DENSS=snow density in snowpack layer
C     DDLYRS=change in DLYRS to maintain snowpack layer depth set
C        in starts.f
C     AREA=grid cell surface area
C
      IF(VHCPW(1,NY,NX).GT.VHCPWX(NY,NX))THEN
      DO 325 L=1,JS-1
      VOLSLX=VOLSL(L,NY,NX)
      IF(VOLSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DDLYXS=(VOLSI(L,NY,NX)-VOLSSL(L,NY,NX)/DENSS(L,NY,NX)
     2-VOLWSL(L,NY,NX)-VOLISL(L,NY,NX))/AREA(3,L,NY,NX)
      DDLYXX=DDLYXS
      IF(DDLYXS.LT.-ZERO.OR.DLYRS(L+1,NY,NX).GT.ZERO)THEN
      DDLYRS=AMIN1(DDLYXS,DLYRS(L+1,NY,NX))
      IFLGLS=1
      ELSE
      DDLYXS=(VOLSL(L,NY,NX)-VOLSSL(L,NY,NX)/DENSS(L,NY,NX)
     2-VOLWSL(L,NY,NX)-VOLISL(L,NY,NX))/AREA(3,L,NY,NX)
      DDLYRS=DDLYXS
      IFLGLS=2
      ENDIF
      ELSE
      DDLYRS=0.0
      IFLGLS=0
      ENDIF
C
C     RESET SNOW LAYER DEPTHS
C
C     CDPTHS=cumulative depth to bottom of snowpack layer
C     DDLYRS=change in DLYRS to maintain snowpack layer depth
C        set in starts.f
C     DLYRS=snowpack layer depth
C
      CDPTHS(L,NY,NX)=CDPTHS(L,NY,NX)+DDLYRS
      DLYRS(L,NY,NX)=CDPTHS(L,NY,NX)-CDPTHS(L-1,NY,NX)
C
C     TRANSFER STATE VARIABLES BETWEEN SNOWPACK LAYERS
C
C     L0,L1=snowpack layer from, to which snowpack contents are
C        redistributed
C     FX=fraction of snowpack contents redistributed from L0 to L1
C     VOLSI=snowpack layer volume set in starts.f
C     VOLSL=current snowpack layer volume
C
      IF(ABS(DDLYRS).GT.ZERO)THEN
C     WRITE(*,1113)'DDLYRS',I,J,NX,NY,L,IFLGLS,DDLYRS,DDLYXS
C    2,CDPTHS(L,NY,NX),DLYRS(L,NY,NX),VOLSI(L,NY,NX),VOLSL(L,NY,NX)
C    3,VOLSLX,DDLYXX,CDPTHS(L+1,NY,NX),DLYRS(L+1,NY,NX)
1113  FORMAT(A8,6I4,30E14.6)
      IF(DDLYRS.GT.0.0)THEN
      L1=L
      L0=L+1
      IF(DDLYRS.LT.DDLYXS)THEN
      FX=1.0
      ELSE
      FX=AMIN1(1.0,DDLYRS*AREA(3,L0,NY,NX)/VOLSL(L0,NY,NX))
      ENDIF
      ELSE
      L1=L+1
      L0=L
      IF(VOLSL(L0,NY,NX).LT.VOLSI(L0,NY,NX))THEN
      FX=0.0
      ELSE
      FX=AMIN1(1.0,-DDLYRS*AREA(3,L0,NY,NX)/VOLSL(L0,NY,NX))
      ENDIF
      ENDIF
C
C     REDISTRIBUTE SNOWPACK CONTENTS
C
C     FX=fraction of snowpack contents redistributed from L0 to L1
C
      IF(FX.GT.0.0)THEN
      FY=1.0-FX
C     IF(IYRC.EQ.2006.AND.I.EQ.361.AND.NX.EQ.1)THEN
C     WRITE(*,5596)'SNOW1',I,J,NX,NY,L,NU(NY,NX),L1,L0,FX,FY
C    3,DDLYRS,VOLSI(L0,NY,NX),VOLSL(L0,NY,NX),VOLSSL(L0,NY,NX)
C    3,VOLWSL(L0,NY,NX),VOLISL(L0,NY,NX),VOLSI(L1,NY,NX)
C    4,VOLSL(L1,NY,NX),VOLSSL(L1,NY,NX),VOLWSL(L1,NY,NX)
C    4,VOLISL(L1,NY,NX),CDPTHS(L0,NY,NX),CDPTHS(L1,NY,NX)
C    5,DENSS(L1,NY,NX),DENSS(L0,NY,NX)
C    5,TKW(L0,NY,NX),TKW(L1,NY,NX)
C    5,VHCPW(L0,NY,NX),VHCPW(L1,NY,NX)
C    5,TKW(L0,NY,NX)*VHCPW(L0,NY,NX)
C    5,VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
C    5,TKW(L0,NY,NX)*VHCPW(L0,NY,NX)
C    5+VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
5596  FORMAT(A8,8I4,100E14.6)
C     ENDIF
C
C     REDISTRIBUTE SNOWPACK CONTENTS TO DESTINATION LAYER
C
C     SNOWPACK SNOW, WATER, VAPOR, ICE, HEAT
C
      VOLSSL(L1,NY,NX)=VOLSSL(L1,NY,NX)+FX*VOLSSL(L0,NY,NX)
      VOLWSL(L1,NY,NX)=VOLWSL(L1,NY,NX)+FX*VOLWSL(L0,NY,NX)
      VOLVSL(L1,NY,NX)=VOLVSL(L1,NY,NX)+FX*VOLVSL(L0,NY,NX)
      VOLISL(L1,NY,NX)=VOLISL(L1,NY,NX)+FX*VOLISL(L0,NY,NX)
      VOLSL(L1,NY,NX)=VOLSSL(L1,NY,NX)/DENSS(L1,NY,NX)
     2+VOLWSL(L1,NY,NX)+VOLISL(L1,NY,NX)
      ENGY1X=VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
      ENGY0X=VHCPW(L0,NY,NX)*TKW(L0,NY,NX)
      ENGY1=ENGY1X+FX*ENGY0X
      VHCPW(L1,NY,NX)=2.095*VOLSSL(L1,NY,NX)
     2+4.19*(VOLWSL(L1,NY,NX)+VOLVSL(L1,NY,NX))
     2+1.9274*VOLISL(L1,NY,NX)
      IF(VHCPW(L1,NY,NX).GT.VHCPWX(NY,NX))THEN
      TKW(L1,NY,NX)=ENGY1/VHCPW(L1,NY,NX)
      ELSE
      TKW(L1,NY,NX)=TKW(L0,NY,NX)
      ENDIF
      TCW(L1,NY,NX)=TKW(L1,NY,NX)-273.15
C
C     SNOWPACK SOLUTES
C
      CO2W(L1,NY,NX)=CO2W(L1,NY,NX)+FX*CO2W(L0,NY,NX)
      CH4W(L1,NY,NX)=CH4W(L1,NY,NX)+FX*CH4W(L0,NY,NX)
      OXYW(L1,NY,NX)=OXYW(L1,NY,NX)+FX*OXYW(L0,NY,NX)
      ZNGW(L1,NY,NX)=ZNGW(L1,NY,NX)+FX*ZNGW(L0,NY,NX)
      ZN2W(L1,NY,NX)=ZN2W(L1,NY,NX)+FX*ZN2W(L0,NY,NX)
      ZN4W(L1,NY,NX)=ZN4W(L1,NY,NX)+FX*ZN4W(L0,NY,NX)
      ZN3W(L1,NY,NX)=ZN3W(L1,NY,NX)+FX*ZN3W(L0,NY,NX)
      ZNOW(L1,NY,NX)=ZNOW(L1,NY,NX)+FX*ZNOW(L0,NY,NX)
      Z1PW(L1,NY,NX)=Z1PW(L1,NY,NX)+FX*Z1PW(L0,NY,NX)
      ZHPW(L1,NY,NX)=ZHPW(L1,NY,NX)+FX*ZHPW(L0,NY,NX)
C
C     SNOWPACK SALTS
C
      IF(ISALTG.NE.0)THEN
      ZALW(L1,NY,NX)=ZALW(L1,NY,NX)+FX*ZALW(L0,NY,NX)
      ZFEW(L1,NY,NX)=ZFEW(L1,NY,NX)+FX*ZFEW(L0,NY,NX)
      ZHYW(L1,NY,NX)=ZHYW(L1,NY,NX)+FX*ZHYW(L0,NY,NX)
      ZCAW(L1,NY,NX)=ZCAW(L1,NY,NX)+FX*ZCAW(L0,NY,NX)
      ZMGW(L1,NY,NX)=ZMGW(L1,NY,NX)+FX*ZMGW(L0,NY,NX)
      ZNAW(L1,NY,NX)=ZNAW(L1,NY,NX)+FX*ZNAW(L0,NY,NX)
      ZKAW(L1,NY,NX)=ZKAW(L1,NY,NX)+FX*ZKAW(L0,NY,NX)
      ZOHW(L1,NY,NX)=ZOHW(L1,NY,NX)+FX*ZOHW(L0,NY,NX)
      ZSO4W(L1,NY,NX)=ZSO4W(L1,NY,NX)+FX*ZSO4W(L0,NY,NX)
      ZCLW(L1,NY,NX)=ZCLW(L1,NY,NX)+FX*ZCLW(L0,NY,NX)
      ZCO3W(L1,NY,NX)=ZCO3W(L1,NY,NX)+FX*ZCO3W(L0,NY,NX)
      ZHCO3W(L1,NY,NX)=ZHCO3W(L1,NY,NX)+FX*ZHCO3W(L0,NY,NX)
      ZALH1W(L1,NY,NX)=ZALH1W(L1,NY,NX)+FX*ZALH1W(L0,NY,NX)
      ZALH2W(L1,NY,NX)=ZALH2W(L1,NY,NX)+FX*ZALH2W(L0,NY,NX)
      ZALH3W(L1,NY,NX)=ZALH3W(L1,NY,NX)+FX*ZALH3W(L0,NY,NX)
      ZALH4W(L1,NY,NX)=ZALH4W(L1,NY,NX)+FX*ZALH4W(L0,NY,NX)
      ZALSW(L1,NY,NX)=ZALSW(L1,NY,NX)+FX*ZALSW(L0,NY,NX)
      ZFEH1W(L1,NY,NX)=ZFEH1W(L1,NY,NX)+FX*ZFEH1W(L0,NY,NX)
      ZFEH2W(L1,NY,NX)=ZFEH2W(L1,NY,NX)+FX*ZFEH2W(L0,NY,NX)
      ZFEH3W(L1,NY,NX)=ZFEH3W(L1,NY,NX)+FX*ZFEH3W(L0,NY,NX)
      ZFEH4W(L1,NY,NX)=ZFEH4W(L1,NY,NX)+FX*ZFEH4W(L0,NY,NX)
      ZFESW(L1,NY,NX)=ZFESW(L1,NY,NX)+FX*ZFESW(L0,NY,NX)
      ZCAOW(L1,NY,NX)=ZCAOW(L1,NY,NX)+FX*ZCAOW(L0,NY,NX)
      ZCACW(L1,NY,NX)=ZCACW(L1,NY,NX)+FX*ZCACW(L0,NY,NX)
      ZCAHW(L1,NY,NX)=ZCAHW(L1,NY,NX)+FX*ZCAHW(L0,NY,NX)
      ZCASW(L1,NY,NX)=ZCASW(L1,NY,NX)+FX*ZCASW(L0,NY,NX)
      ZMGOW(L1,NY,NX)=ZMGOW(L1,NY,NX)+FX*ZMGOW(L0,NY,NX)
      ZMGCW(L1,NY,NX)=ZMGCW(L1,NY,NX)+FX*ZMGCW(L0,NY,NX)
      ZMGHW(L1,NY,NX)=ZMGHW(L1,NY,NX)+FX*ZMGHW(L0,NY,NX)
      ZMGSW(L1,NY,NX)=ZMGSW(L1,NY,NX)+FX*ZMGSW(L0,NY,NX)
      ZNACW(L1,NY,NX)=ZNACW(L1,NY,NX)+FX*ZNACW(L0,NY,NX)
      ZNASW(L1,NY,NX)=ZNASW(L1,NY,NX)+FX*ZNASW(L0,NY,NX)
      ZKASW(L1,NY,NX)=ZKASW(L1,NY,NX)+FX*ZKASW(L0,NY,NX)
      H0PO4W(L1,NY,NX)=H0PO4W(L1,NY,NX)+FX*H0PO4W(L0,NY,NX)
      H3PO4W(L1,NY,NX)=H3PO4W(L1,NY,NX)+FX*H3PO4W(L0,NY,NX)
      ZFE1PW(L1,NY,NX)=ZFE1PW(L1,NY,NX)+FX*ZFE1PW(L0,NY,NX)
      ZFE2PW(L1,NY,NX)=ZFE2PW(L1,NY,NX)+FX*ZFE2PW(L0,NY,NX)
      ZCA0PW(L1,NY,NX)=ZCA0PW(L1,NY,NX)+FX*ZCA0PW(L0,NY,NX)
      ZCA1PW(L1,NY,NX)=ZCA1PW(L1,NY,NX)+FX*ZCA1PW(L0,NY,NX)
      ZCA2PW(L1,NY,NX)=ZCA2PW(L1,NY,NX)+FX*ZCA2PW(L0,NY,NX)
      ZMG1PW(L1,NY,NX)=ZMG1PW(L1,NY,NX)+FX*ZMG1PW(L0,NY,NX)
      ENDIF
C
C     REDISTRIBUTE ALL SNOWPACK CONTENTS FROM SOURCE LAYER
C
C     SNOWPACK SNOW,WATER,ICE,HEAT
C
      VOLSSL(L0,NY,NX)=FY*VOLSSL(L0,NY,NX)
      VOLWSL(L0,NY,NX)=FY*VOLWSL(L0,NY,NX)
      VOLVSL(L0,NY,NX)=FY*VOLVSL(L0,NY,NX)
      VOLISL(L0,NY,NX)=FY*VOLISL(L0,NY,NX)
      VOLSL(L0,NY,NX)=VOLSSL(L0,NY,NX)/DENSS(L0,NY,NX)
     2+VOLWSL(L0,NY,NX)+VOLISL(L0,NY,NX)
      ENGY0=FY*ENGY0X
      VHCPW(L0,NY,NX)=2.095*VOLSSL(L0,NY,NX)
     2+4.19*(VOLWSL(L0,NY,NX)+VOLVSL(L0,NY,NX))
     2+1.9274*VOLISL(L0,NY,NX)
      IF(VHCPW(L0,NY,NX).GT.VHCPWX(NY,NX))THEN
      TKW(L0,NY,NX)=ENGY0/VHCPW(L0,NY,NX)
      ELSE
      TKW(L0,NY,NX)=TKW(L1,NY,NX)
      ENDIF
      TCW(L0,NY,NX)=TKW(L0,NY,NX)-273.15
C
C     SNOWPACK SOLUTES
C
      CO2W(L0,NY,NX)=FY*CO2W(L0,NY,NX)
      CH4W(L0,NY,NX)=FY*CH4W(L0,NY,NX)
      OXYW(L0,NY,NX)=FY*OXYW(L0,NY,NX)
      ZNGW(L0,NY,NX)=FY*ZNGW(L0,NY,NX)
      ZN2W(L0,NY,NX)=FY*ZN2W(L0,NY,NX)
      ZN4W(L0,NY,NX)=FY*ZN4W(L0,NY,NX)
      ZN3W(L0,NY,NX)=FY*ZN3W(L0,NY,NX)
      ZNOW(L0,NY,NX)=FY*ZNOW(L0,NY,NX)
      Z1PW(L0,NY,NX)=FY*Z1PW(L0,NY,NX)
      ZHPW(L0,NY,NX)=FY*ZHPW(L0,NY,NX)
C
C     SNOWPACK SALTS
C
      IF(ISALTG.NE.0)THEN
      ZALW(L0,NY,NX)=FY*ZALW(L0,NY,NX)
      ZFEW(L0,NY,NX)=FY*ZFEW(L0,NY,NX)
      ZHYW(L0,NY,NX)=FY*ZHYW(L0,NY,NX)
      ZCAW(L0,NY,NX)=FY*ZCAW(L0,NY,NX)
      ZMGW(L0,NY,NX)=FY*ZMGW(L0,NY,NX)
      ZNAW(L0,NY,NX)=FY*ZNAW(L0,NY,NX)
      ZKAW(L0,NY,NX)=FY*ZKAW(L0,NY,NX)
      ZOHW(L0,NY,NX)=FY*ZOHW(L0,NY,NX)
      ZSO4W(L0,NY,NX)=FY*ZSO4W(L0,NY,NX)
      ZCLW(L0,NY,NX)=FY*ZCLW(L0,NY,NX)
      ZCO3W(L0,NY,NX)=FY*ZCO3W(L0,NY,NX)
      ZHCO3W(L0,NY,NX)=FY*ZHCO3W(L0,NY,NX)
      ZALH1W(L0,NY,NX)=FY*ZALH1W(L0,NY,NX)
      ZALH2W(L0,NY,NX)=FY*ZALH2W(L0,NY,NX)
      ZALH3W(L0,NY,NX)=FY*ZALH3W(L0,NY,NX)
      ZALH4W(L0,NY,NX)=FY*ZALH4W(L0,NY,NX)
      ZALSW(L0,NY,NX)=FY*ZALSW(L0,NY,NX)
      ZFEH1W(L0,NY,NX)=FY*ZFEH1W(L0,NY,NX)
      ZFEH2W(L0,NY,NX)=FY*ZFEH2W(L0,NY,NX)
      ZFEH3W(L0,NY,NX)=FY*ZFEH3W(L0,NY,NX)
      ZFEH4W(L0,NY,NX)=FY*ZFEH4W(L0,NY,NX)
      ZFESW(L0,NY,NX)=FY*ZFESW(L0,NY,NX)
      ZCAOW(L0,NY,NX)=FY*ZCAOW(L0,NY,NX)
      ZCACW(L0,NY,NX)=FY*ZCACW(L0,NY,NX)
      ZCAHW(L0,NY,NX)=FY*ZCAHW(L0,NY,NX)
      ZCASW(L0,NY,NX)=FY*ZCASW(L0,NY,NX)
      ZMGOW(L0,NY,NX)=FY*ZMGOW(L0,NY,NX)
      ZMGCW(L0,NY,NX)=FY*ZMGCW(L0,NY,NX)
      ZMGHW(L0,NY,NX)=FY*ZMGHW(L0,NY,NX)
      ZMGSW(L0,NY,NX)=FY*ZMGSW(L0,NY,NX)
      ZNACW(L0,NY,NX)=FY*ZNACW(L0,NY,NX)
      ZNASW(L0,NY,NX)=FY*ZNASW(L0,NY,NX)
      ZKASW(L0,NY,NX)=FY*ZKASW(L0,NY,NX)
      H0PO4W(L0,NY,NX)=FY*H0PO4W(L0,NY,NX)
      H3PO4W(L0,NY,NX)=FY*H3PO4W(L0,NY,NX)
      ZFE1PW(L0,NY,NX)=FY*ZFE1PW(L0,NY,NX)
      ZFE2PW(L0,NY,NX)=FY*ZFE2PW(L0,NY,NX)
      ZCA0PW(L0,NY,NX)=FY*ZCA0PW(L0,NY,NX)
      ZCA1PW(L0,NY,NX)=FY*ZCA1PW(L0,NY,NX)
      ZCA2PW(L0,NY,NX)=FY*ZCA2PW(L0,NY,NX)
      ZMG1PW(L0,NY,NX)=FY*ZMG1PW(L0,NY,NX)
      ENDIF
C     IF(IYRC.EQ.2006.AND.I.EQ.361.AND.NX.EQ.1)THEN
C     WRITE(*,5596)'SNOW2',I,J,NX,NY,L,NU(NY,NX),L1,L0,FX,FY
C    3,DDLYRS,VOLSI(L0,NY,NX),VOLSL(L0,NY,NX),VOLSSL(L0,NY,NX)
C    3,VOLWSL(L0,NY,NX),VOLISL(L0,NY,NX),VOLSI(L1,NY,NX)
C    4,VOLSL(L1,NY,NX),VOLSSL(L1,NY,NX),VOLWSL(L1,NY,NX)
C    4,VOLISL(L1,NY,NX),CDPTHS(L0,NY,NX),CDPTHS(L1,NY,NX)
C    5,DENSS(L1,NY,NX),DENSS(L0,NY,NX)
C    5,TKW(L0,NY,NX),TKW(L1,NY,NX)
C    5,VHCPW(L0,NY,NX),VHCPW(L1,NY,NX)
C    5,TKW(L0,NY,NX)*VHCPW(L0,NY,NX)
C    5,VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
C    5,TKW(L0,NY,NX)*VHCPW(L0,NY,NX)
C    5+VHCPW(L1,NY,NX)*TKW(L1,NY,NX)
C     ENDIF
      ENDIF
      ENDIF
325   CONTINUE
      ENDIF
C
C     PROFILE RELAYERING WITH EVAPORATION, FREEZE-THAW, EROSION, SOC
C
C     IERSNG=erosion flag from site file
C     BKDS=bulk density (0=water,>0=soil)
C     CDPTH=cumulative depth to bottom of soil layer
C     CDPTHY=CDPTH excluding freeze-thaw effects
C
      IF(IERSNG.GE.0)THEN
      IF(BKDS(NU(NY,NX),NY,NX).LE.ZERO)THEN
      ICHKLX=0
      ELSE
      ICHKLX=1
      ENDIF
      DO 225 LX=NL(NY,NX),NU(NY,NX),-1
      CDPTHX(LX,NY,NX)=CDPTH(LX,NY,NX)
      CDPTHY(LX,NY,NX)=CDPTH(LX,NY,NX)
C
C     IF LAYER IS IN A POND
C
C     DLYR=layer depth
C     DDLYXP=change in DLYR to maintain pond layer depth
C        set in soil file
C     DDLYX,DDLYR=cumulative DDLYXP from bottom to top of profile
C     VOLW,VOLI=pond water,ice content
C     AREA=grid cell surface area
C     IFLGL=pond boundary flag:layer underneath is 1:pond, 2:soil
C
      IF(BKDS(LX,NY,NX).LE.ZERO)THEN
      IF(BKDS(LX+1,NY,NX).GT.ZERO.OR.ICHKLX.EQ.1)THEN
      DDLYXP=DLYR(3,LX,NY,NX)-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))
     2/AREA(3,LX,NY,NX)
      DDLYX(LX,1)=DDLYXP+DDLYX(LX+1,1)
      DDLYR(LX,1)=DDLYX(LX+1,1)
      IFLGL(LX,1)=2
      ELSE
      DDLYXP=DLYRI(3,LX,NY,NX)-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))
     2/AREA(3,LX,NY,NX)
      DPTWI=(VOLW(LX+1,NY,NX)+VOLI(LX+1,NY,NX))/AREA(3,LX,NY,NX)
      IF(DDLYXP.LT.-ZERO.OR.DPTWI.GT.ZERO)THEN
      DDLYX(LX,1)=DDLYXP+DDLYX(LX+1,1)
      DDLYR(LX,1)=AMIN1(DDLYX(LX+1,1),DPTWI)
      IF(DPTWI.GT.ZERO)THEN
      IFLGL(LX,1)=1
      ELSE
      IFLGL(LX,1)=2
      ENDIF
      ELSE
      DDLYXP=DLYR(3,LX,NY,NX)-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))
     2/AREA(3,LX,NY,NX)
      DDLYX(LX,1)=DDLYXP+DDLYX(LX+1,1)
      DDLYR(LX,1)=DDLYX(LX+1,1)
      IFLGL(LX,1)=2
      ENDIF
      ENDIF
      IF(LX.EQ.NU(NY,NX).OR.BKDS(LX-1,NY,NX).GT.ZERO)THEN
      DDLYX(LX-1,1)=DDLYX(LX,1)
      DDLYR(LX-1,1)=DDLYX(LX,1)
      IFLGL(LX-1,1)=1
      ENDIF
      DDLYX(LX,4)=0.0
      DDLYR(LX,4)=0.0
      IFLGL(LX,4)=0
      DDLYX(LX,5)=0.0
      DDLYR(LX,5)=0.0
      IFLGL(LX,5)=0
      DDLYX(LX,6)=0.0
      DDLYR(LX,6)=0.0
      IFLGL(LX,6)=0
C     WRITE(*,1123)'LAKE',I,J,NX,NY,LX,IFLGL(LX,1),ICHKLX
C    2,DDLYXP,DDLYX(LX,1),DDLYR(LX,1),VOLAI(LX,NY,NX),VOLY(LX,NY,NX)
C    3,VOLT(LX,NY,NX),VOLA(LX,NY,NX),VOLW(LX,NY,NX),VOLI(LX,NY,NX)
C    3,DPTWI,POROS(LX,NY,NX),DLYR(3,LX,NY,NX),BKDS(LX,NY,NX)
C    4,CDPTH(LX-1,NY,NX),CDPTH(LX,NY,NX)
C    4,CDPTH(LX,NY,NX)-CDPTH(LX-1,NY,NX)
C    5,VOLW(LX,NY,NX)+VOLI(LX,NY,NX),DLYR(3,LX,NY,NX)
C    6-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))/AREA(3,LX,NY,NX)
C    7,DVOLI(LX,NY,NX)
1123  FORMAT(A8,7I4,20E16.8)
C
C     IF LAYER IS IN SOIL
C
      ELSE
C
C     CHANGE IN LAYER DEPTH WITH FREEZE-THAW
C
C     DVOLI=change in soil ice content
C     DDLYXF=change in DLYR to maintain soil layer depth
C        set in soil file
C     DENSJ=1.0-DENSI:ice density
C     AREA=grid cell surface area
C     DDLYX,DDLYR=cumulative DDLYXF from bottom to top of profile
C
      IF(ABS(DVOLI(LX,NY,NX)).GT.ZEROS(NY,NX))THEN
      DDLYXF=DVOLI(LX,NY,NX)*DENSJ/AREA(3,LX,NY,NX)
      IF(LX.EQ.NL(NY,NX))THEN
      DDLYX(LX,4)=DDLYXF
      DDLYR(LX,4)=0.0
      IFLGL(LX,4)=0
      ELSE
      DDLYX(LX,4)=DDLYXF+DDLYX(LX+1,4)
      DDLYR(LX,4)=DDLYX(LX+1,4)
C    2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
      IFLGL(LX,4)=0
      IF(LX.EQ.NU(NY,NX).OR.BKDS(LX-1,NY,NX).LE.ZERO)THEN
      DDLYX(LX-1,4)=DDLYX(LX,4)
      DDLYR(LX-1,4)=DDLYX(LX,4)
C    2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
      IFLGL(LX-1,4)=0
      ENDIF
      ENDIF
      ELSE
      DDLYXF=0.0
      IF(LX.EQ.NL(NY,NX))THEN
      DDLYX(LX,4)=0.0
      DDLYR(LX,4)=0.0
      IFLGL(LX,4)=0
      ELSE
      DDLYX(LX,4)=DDLYX(LX+1,4)
      DDLYR(LX,4)=DDLYX(LX+1,4)
      IFLGL(LX,4)=0
      IF(LX.EQ.NU(NY,NX))THEN
      DDLYX(LX-1,4)=DDLYX(LX,4)
      DDLYR(LX-1,4)=DDLYX(LX,4)
      IFLGL(LX-1,4)=0
      ENDIF
      ENDIF
      ENDIF
      IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,1126)'FREEZ',I,J,NFZ,NX,NY,LX
C    2,DDLYX(LX,4),DDLYR(LX,4),DDLYXF
C    2,VOLP(LX,NY,NX),VOLT(LX,NY,NX),DVOLI(LX,NY,NX),DVOLW(LX,NY,NX)
C    5,TFLW(LX,NY,NX),FINH(LX,NY,NX)
C    2,TWFLFL(LX,NY,NX),TUPWTR(LX,NY,NX),FLU(LX,NY,NX)
C    2,VOLA(LX,NY,NX)-VOLI(LX,NY,NX)-VOLW(LX,NY,NX)
C    4+VOLAH(LX,NY,NX)-VOLIH(LX,NY,NX)-VOLWH(LX,NY,NX)
C    3,BKDSI(LX,NY,NX),BKDS(LX,NY,NX)
C    4,DLYR(3,LX,NY,NX),CDPTH(LX,NY,NX)
C     IF(LX.EQ.NU(NY,NX))THEN
C     WRITE(*,1126)'FREEZ0',I,J,NFZ,NX,NY,LX
C    2,DDLYX(LX-1,4),DDLYR(LX-1,4)
C    2,DDLYX(LX,4),DDLYR(LX,4),DDLYXF,DDLYX0,DVOLI(LX,NY,NX)*DENSJ
C    2,BKDS(LX,NY,NX),BKDSI(LX,NY,NX)
C    4,DLYR(3,LX,NY,NX),CDPTH(LX-1,NY,NX),CDPTH(LX,NY,NX)
C    5,DVOLI(LX,NY,NX),DVOLW(LX,NY,NX),VOLI(LX,NY,NX),VOLW(LX,NY,NX)
1126  FORMAT(A8,6I4,30E14.6)
C     ENDIF
      ENDIF
      UDVOLI=UDVOLI+DVOLI(L,NY,NX)
      UDLYXF=UDLYXF+DDLYXF
C     IF(I.EQ.365.AND.J.EQ.24.AND.LX.EQ.NU(NY,NX))THEN
C     WRITE(*,1111)'VOLISO',I,J,NX,NY,LX,IFLGL(LX,NN)
C    5,VOLISO,UDVOLI,UDLYXF
C     ENDIF
C
C     CHANGE IN LAYER DEPTH WITH EROSION
C
C     IERSNG=erosion flag from site file
C     TSEDER=net sediment flux
C     DDLYXE=change in DLYR to maintain soil layer depth
C        set in soil file
C     BKVLNU,VOLX=bulk density,volume of soil surface layer
C     DDLYX,DDLYR=DDLYXE from bottom to top of profile
C
      IF((IERSNG.EQ.1.OR.IERSNG.EQ.3)
     2.AND.(ABS(TSEDER(NY,NX)).GT.ZEROS(NY,NX)
     3.OR.TSEDSK(NY,NX).GT.ZEROS(NY,NX)))THEN
      IF(LX.EQ.NL(NY,NX))THEN
      IF(BKDS(NU(NY,NX),NY,NX).GT.ZERO)THEN
      DDLYXE=-TSEDER(NY,NX)/(AREA(3,LX,NY,NX)
     2*BKVLNU(NY,NX)/VOLX(NU(NY,NX),NY,NX))
      ELSEIF(TSEDSK(NY,NX).GT.ZEROS(NY,NX))THEN
      DDLYXE=-TSEDSK(NY,NX)/(AREA(3,LX,NY,NX)
     2*BKDSI(NUL(NY,NX),NY,NX))
      ELSE
      DDLYXE=0.0
      ENDIF
      ENDIF
      DDLYX(LX,5)=DDLYXE
      DDLYR(LX,5)=DDLYXE
      IFLGL(LX,5)=1
      ELSE
      DDLYX(LX,5)=0.0
      DDLYR(LX,5)=0.0
      IFLGL(LX,5)=0
      ENDIF
C     IF(NX.EQ.1)THEN
C     WRITE(*,1121)'SED',I,J,NFZ,NX,NY,LX,NUL(NY,NX)
C    2,DDLYXE,DDLYX(LX,5),TSEDER(NY,NX),TSEDSK(NY,NX)
C    3,BKDS(NU(NY,NX),NY,NX),DLYR(3,LX,NY,NX)
C    3,BKVLNU(NY,NX),VOLX(NU(NY,NX),NY,NX)
C    4,BKVL(NUL(NY,NX),NY,NX),VOLX(NUL(NY,NX),NY,NX)
C    4,SAND(NUL(NY,NX),NY,NX),SILT(NUL(NY,NX),NY,NX)
C    4,CLAY(NUL(NY,NX),NY,NX),ORGC(NUL(NY,NX),NY,NX)
1121  FORMAT(A8,7I4,20E14.6)
C     ENDIF
C
C     CHANGE IN LAYER DEPTH WITH SOC GAIN OR LOSS
C
C     IERSNG=erosion flag from site file
C     DORGC=change in total SOC
C     DDLYXC=change in DLYR to maintain soil layer depth
C        set in soil file
C     FHOL=macropore fraction
C     BKDSI=bulk density as read from soil file
C     AREA=grid cell surface area
C     DDLYX,DDLYR=cumulative DDLYXC from bottom to top of profile
C     DLYRI,DLYR=initial layer depth from soil file,
C        current layer depth
C
      IF((IERSNG.EQ.2.OR.IERSNG.EQ.3)
     2.AND.ABS(DORGC(LX,NY,NX)).GT.ZEROS(NY,NX))THEN
      DDLYXC=1.82E-06*DORGC(LX,NY,NX)/AREA(3,LX,NY,NX)
     2/((1.0-FHOL(LX,NY,NX))*BKDSI(LX,NY,NX))
      IF(LX.EQ.NL(NY,NX).OR.BKDS(LX+1,NY,NX).LE.ZERO)THEN
      DDLYX(LX,6)=DDLYXC
      DDLYR(LX,6)=0.0
      IFLGL(LX,6)=1
      ELSE
      DDLYX(LX,6)=DDLYXC+DDLYX(LX+1,6)
      DDLYR(LX,6)=DDLYX(LX+1,6)
     2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
      IFLGL(LX,6)=1
      IF(LX.EQ.NU(NY,NX).OR.BKDS(LX-1,NY,NX).LE.ZERO)THEN
      DDLYX(LX-1,6)=DDLYX(LX,6)
      DDLYR(LX-1,6)=DDLYX(LX,6)
C    2+DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
      IFLGL(LX-1,6)=1
      ENDIF
      ENDIF
      ELSE
      IF(LX.EQ.NL(NY,NX))THEN
      DDLYX(LX,6)=0.0
      DDLYR(LX,6)=0.0
      IFLGL(LX,6)=0
      ELSE
      DDLYX(LX,6)=DDLYX(LX+1,6)
      DDLYR(LX,6)=DDLYX(LX+1,6)
      IFLGL(LX,6)=0
      ENDIF
      ENDIF
      DDLYX(LX,1)=0.0
      DDLYR(LX,1)=0.0
      IFLGL(LX,1)=0
C     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,1111)'SOC',I,J,NFZ,NX,NY,LX,IFLGL(LX,6),DDLYX(LX,6)
C    2,DDLYR(LX,6),DDLYXC,DORGC(LX,NY,NX),BKDS(LX,NY,NX)
C    3,DLYRI(3,LX,NY,NX),DLYR(3,LX,NY,NX)
C    4,DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
C    5,CDPTH(0,NY,NX)
C    6,ORGCX(LX,NY,NX),ORGC(LX,NY,NX),ORGCC(LX,NY,NX)
1111  FORMAT(A8,7I4,40E14.6)
C     IF(LX.EQ.NU(NY,NX).OR.BKDS(LX-1,NY,NX).LE.ZERO)THEN
C     WRITE(*,1111)'SOC0',I,J,NFZ,NX,NY,LX,IFLGL(LX-1,6),DDLYX(LX-1,6)
C    2,DDLYR(LX-1,6),DDLYXC,DORGC(LX,NY,NX),BKDS(LX,NY,NX)
C    3,DLYRI(3,LX,NY,NX),DLYR(3,LX,NY,NX)
C    4,DLYRI(3,LX,NY,NX)-DLYR(3,LX,NY,NX)
C    5,CDPTH(0,NY,NX)
C     ENDIF
C     ENDIF
      ENDIF
C
C     RESET LAYER DEPTHS
C
      DO 200 NN=1,6
C     IF(ABS(DDLYX(LX,NN)).GT.ZERO)IFLGS(NY,NX)=1
      IF(NN.NE.2.AND.NN.NE.3)THEN
C
C     POND LAYER DEPTHS
C
C     BKDS=bulk density (0=water,>0=soil)
C     CDPTH=cumulative depth to bottom of soil layer
C     DDLYX,DDLYR=cumulative change in layer depth from bottom
C        to top of profile
C     IFLGL=pond boundary flag:layer underneath is 1:pond, 2:soil
C
      IF(BKDS(LX,NY,NX).LE.ZERO)THEN
      IF(IFLGL(LX,NN).NE.0)THEN
      CDPTH(LX,NY,NX)=CDPTH(LX,NY,NX)+DDLYR(LX,NN)
      CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)
      IF(LX.NE.NU(NY,NX).AND.IFLGL(LX,1).EQ.2)THEN
      DO 201 LL=LX-1,0,-1
      CDPTH(LL,NY,NX)=CDPTH(LL,NY,NX)+DDLYX(LX,NN)
      CDPTHY(LL,NY,NX)=CDPTHY(LL,NY,NX)+DDLYX(LX,NN)
C     WRITE(*,1116)'DDLYRZ',I,J,NX,NY,LX,LL,IFLGL(LX,1),IFLGL(LX-1,1)
C    2,DDLYX(LX,NN),CDPTH(LL,NY,NX),DLYR(3,LX,NY,NX)
C    2,VOLW(LX,NY,NX),DVOLW(LX,NY,NX),AREA(3,LX,NY,NX)
1116  FORMAT(A8,8I4,40E16.8)
201   CONTINUE
      DDLYX(LX,NN)=0.0
      ENDIF
C     ENDIF
C     IF(J.EQ.24)THEN
C     WRITE(*,1117)'POND',I,J,NX,NY,LX,NN,IFLGL(LX,NN),ICHKLX
C    2,DDLYXP,DDLYX(LX,NN),DDLYR(LX,NN),DLYR(3,LX,NY,NX)
C    3,VOLW(LX,NY,NX),VOLI(LX,NY,NX)
C    4,CDPTH(LX-1,NY,NX),CDPTH(LX,NY,NX)
C    4,CDPTH(LX,NY,NX)-CDPTH(LX-1,NY,NX)
C    5,VOLW(LX,NY,NX)+VOLI(LX,NY,NX)
1117  FORMAT(A8,8I4,12E16.8)
C     ENDIF
      IF(LX.EQ.NU(NY,NX))THEN
      CDPTH(LX-1,NY,NX)=CDPTH(LX,NY,NX)
     2-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))/AREA(3,LX,NY,NX)
      CDPTHY(LX-1,NY,NX)=CDPTHY(LX,NY,NX)
     2-(VOLW(LX,NY,NX)+VOLI(LX,NY,NX))/AREA(3,LX,NY,NX)
C     WRITE(*,1128)'POND0',I,J,NX,NY,LX,NN,IFLGL(LX-1,NN)
C    2,DDLYR(LX-1,NN),DDLYR(LX,NN),CDPTH(LX-1,NY,NX),CDPTH(LX,NY,NX)
C    4,DLYR(3,LX,NY,NX)
      ENDIF
      ENDIF
C
C     SOIL LAYER DEPTHS
C
      ELSE
C     IF(DDLYR(L,NN).NE.0.OR.DDLYR(LX,NN).NE.0)THEN
C
C     CDPTH=cumulative depth to bottom of soil layer
C     DDLYX,DDLYR=cumulative change in layer depth from bottom
C        to top of profile
C
C     FREEZE-THAW
C
      IF(NN.EQ.4)THEN
      CDPTH(LX,NY,NX)=CDPTH(LX,NY,NX)+DDLYR(LX,NN)
C     CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)
C     WRITE(*,1127)'DFREEZ',I,J,NX,NY,LX,IFLGL(LX,NN)
C    2,DDLYR(LX,NN),CDPTH(LX,NY,NX),DLYR(3,LX,NY,NX),DDLYXF
C    5,CDPTH(LX,NY,NX)-CDPTH(LX-1,NY,NX)
      IF(LX.EQ.NU(NY,NX))THEN
      CDPTH(LX-1,NY,NX)=CDPTH(LX-1,NY,NX)+DDLYR(LX-1,NN)
C     CDPTHY(LX-1,NY,NX)=CDPTHY(LX-1,NY,NX)+DDLYR(LX-1,NN)
C     WRITE(*,1127)'DFREEZ0',I,J,NX,NY,LX,IFLGL(LX-1,NN)
C    2,DDLYR(LX-1,NN),DDLYR(LX,NN),CDPTH(LX-1,NY,NX)
C    3,CDPTH(LX,NY,NX),DLYR(3,LX,NY,NX),DDLYXF
C    5,CDPTH(LX,NY,NX)-CDPTH(LX-1,NY,NX)
1127  FORMAT(A8,6I4,30E16.8)
      ENDIF
      ENDIF
C
C     RESET SURFACE ELEVATION FOR SOIL EROSION
C
      IF(NN.EQ.5.AND.IFLGL(LX,NN).EQ.1)THEN
      CDPTH(LX,NY,NX)=CDPTH(LX,NY,NX)+DDLYR(LX,NN)
      CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)
C     WRITE(*,1122)'CDSED',I,J,NX,NY,LX,IFLGL(LX,5),DDLYR(LX,5)
C    2,CDPTH(LX,NY,NX),DLYR(3,LX,NY,NX)
1122  FORMAT(A8,6I4,12E16.8)
      IF(LX.EQ.NU(NY,NX))THEN
      CDPTH(LX-1,NY,NX)=CDPTH(LX-1,NY,NX)+DDLYR(LX,NN)
      CDPTHY(LX-1,NY,NX)=CDPTHY(LX-1,NY,NX)+DDLYR(LX,NN)
C     WRITE(*,1122)'CDSED0',I,J,NX,NY,LX,IFLGL(LX,5),DDLYR(LX,5)
C    2,CDPTH(LX-1,NY,NX),CDPTH(LX,NY,NX),DLYR(3,LX,NY,NX)
      ENDIF
      ENDIF
C
C     SET SOIL LAYER DEPTHS FOR CHANGES IN SOC
C
      IF(NN.EQ.6.AND.IFLGL(LX,NN).EQ.1)THEN
      CDPTH(LX,NY,NX)=CDPTH(LX,NY,NX)+DDLYR(LX,NN)
      CDPTHY(LX,NY,NX)=CDPTHY(LX,NY,NX)+DDLYR(LX,NN)
C     IF(NX.EQ.1)THEN
C     WRITE(*,1128)'DSOC',I,J,NX,NY,LX,NN,IFLGL(LX,NN)
C    2,DDLYR(LX,NN),DDLYXC,CDPTH(LX,NY,NX),DLYR(3,LX,NY,NX)
C    3,CDPTH(0,NY,NX)
1128  FORMAT(A8,7I4,30E16.8)
C     ENDIF
      IF(LX.EQ.NU(NY,NX).OR.BKDS(LX-1,NY,NX).LE.ZERO)THEN
      CDPTH(LX-1,NY,NX)=CDPTH(LX-1,NY,NX)+DDLYR(LX-1,NN)
      CDPTHY(LX-1,NY,NX)=CDPTHY(LX-1,NY,NX)+DDLYR(LX-1,NN)
C     IF(NX.EQ.1)THEN
C     WRITE(*,1128)'DSOC0',I,J,NX,NY,LX,NN,IFLGL(LX-1,NN)
C    2,DDLYR(LX-1,NN),DDLYXC,CDPTH(LX-1,NY,NX)
C    3,DLYR(3,LX-1,NY,NX),CDPTH(0,NY,NX)
C     ENDIF
      IF(BKDS(LX-1,NY,NX).LE.ZERO)THEN
      DO 255 LY=LX-2,0,-1
      IF(BKDS(LY+1,NY,NX).LE.ZERO)THEN
      CDPTH(LY,NY,NX)=CDPTH(LY,NY,NX)+DDLYR(LX-1,NN)
      CDPTHY(LY,NY,NX)=CDPTHY(LY,NY,NX)+DDLYR(LX-1,NN)
C     WRITE(*,1129)'DCDPTH',I,J,NX,NY,LX,LY,NN,IFLGL(LX,NN)
C    2,DDLYR(LX-1,NN),CDPTH(LY,NY,NX),DLYR(3,LY,NY,NX)
C    3,CDPTH(0,NY,NX)
1129  FORMAT(A8,8I4,20E16.8)
      ENDIF
255   CONTINUE
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C     ENDIF
      ENDIF
200   CONTINUE
      VOLY(LX,NY,NX)=VOLX(LX,NY,NX)
225   CONTINUE
      IF(BKDS(NU(NY,NX),NY,NX).LE.ZERO)THEN
      VOLY(0,NY,NX)=VOLW(0,NY,NX)+VOLI(0,NY,NX)
      ELSE
      VOLY(0,NY,NX)=VOLX(0,NY,NX)
      ENDIF
C
C     RECALCULATE SOIL LAYER THICKNESS
C
C     DLYR=layer depth
C     CDPTH=cumulative depth to bottom of soil layer
C     IFLGL=pond boundary flag:layer underneath is 1:pond, 2:soil
C     BKDS=bulk density (0=water,>0=soil)
C     DDLYRX=total change in DLYR from all disturbances used to
C        calculate redistribution of soil material
C     DLYRI=initial layer depth from soil file
C     DDLYRY=change in DLYR to restore DLYRI
C     DPTH=depth to midpoint of soil layer
C     CDPTHY=CDPTH used to calculate DDLYRX excluding
C        freeze-thaw effects
C     CDPTHZ=cumulative depth to bottom of soil layer from current
C        soil surface depth
C
      ICHKL=0
      TDVOLI=0.0
      DO 245 L=NU(NY,NX),NL(NY,NX)-1
      DO 230 NN=1,3
      IF(NN.EQ.1)THEN
      DLYR(3,L,NY,NX)=CDPTH(L,NY,NX)-CDPTH(L-1,NY,NX)
      DLYRXX=DLYR(3,L,NY,NX)
      IF(IFLGL(L,1).EQ.0.AND.IFLGL(L+1,1).NE.0)THEN
      DDLYRX(NN)=0.0
      IF(BKDS(L,NY,NX).LE.ZERO)THEN
      DDLYRY(L)=DLYRI(3,L,NY,NX)-DLYR(3,L,NY,NX)
      ELSE
      DDLYRY(L)=0.0
      ENDIF
      ICHKL=1
      ELSEIF(IFLGL(L,1).EQ.2.AND.(IFLGL(L+1,1).EQ.0
     3.OR.DLYR(3,L,NY,NX).LE.DLYRI(3,L,NY,NX)))THEN
      DDLYRX(NN)=0.0
      IF(L.EQ.NU(NY,NX).OR.ICHKL.EQ.0)THEN
      DDLYRY(L)=0.0
      ELSE
      DDLYRY(L)=DDLYRY(L-1)
      ENDIF
      IF(IFLGL(L,1).EQ.2.AND.IFLGL(L+1,1).EQ.0)ICHKL=0
      ELSE
      IF(ICHKL.EQ.0)THEN
      DDLYRX(NN)=DLYRI(3,L,NY,NX)-DLYR(3,L,NY,NX)
      DDLYRY(L)=DDLYRX(NN)
      ELSE
      DDLYRX(NN)=0.0
      DDLYRY(L)=DDLYRY(L-1)
      ENDIF
      ENDIF
      CDPTH(L,NY,NX)=CDPTH(L,NY,NX)+DDLYRY(L)
      CDPTHY(L,NY,NX)=CDPTHY(L,NY,NX)+DDLYRY(L)
      DLYR(3,L,NY,NX)=DLYR(3,L,NY,NX)+DDLYRY(L)
      DPTH(L,NY,NX)=0.5*(CDPTH(L,NY,NX)+CDPTH(L-1,NY,NX))
      CDPTHZ(L,NY,NX)=CDPTH(L,NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX)
      IF(L.EQ.NL(NY,NX)-1)THEN
      DLYR(3,L+1,NY,NX)=CDPTH(L+1,NY,NX)-CDPTH(L,NY,NX)
      DPTH(L+1,NY,NX)=0.5*(CDPTH(L+1,NY,NX)+CDPTH(L,NY,NX))
      CDPTHZ(L+1,NY,NX)=CDPTH(L+1,NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX)
      ENDIF
      IF(L.EQ.NU(NY,NX))THEN
      DPTHZ(L,NY,NX)=0.5*CDPTHZ(L,NY,NX)
      ELSE
      DPTHZ(L,NY,NX)=0.5*(CDPTHZ(L,NY,NX)+CDPTHZ(L-1,NY,NX))
      ENDIF
      IF(BKDS(L,NY,NX).GT.ZERO)THEN
      DDLYRX(NN)=CDPTHY(L,NY,NX)-CDPTHX(L,NY,NX)
      TDVOLI=TDVOLI+DVOLI(L,NY,NX)
      DDLYRX(NN)=DDLYRX(NN)-TDVOLI*DENSJ/AREA(3,L,NY,NX)
C     IF(L.EQ.NU(NY,NX))THEN
C     DDLYRX(NN)=DDLYRX(NN)-DDLYR(L-1,6)
C     ENDIF
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,1114)'DDLYR',I,J,NFZ,NX,NY,L,NN,IFLGL(L,1),IFLGL(L+1,1)
C    2,ICHKL,DDLYRX(NN),DDLYRY(L),(DDLYR(L,NM),NM=1,6)
C    3,TSEDSK(NY,NX),DVOLI(L,NY,NX),TDVOLI
C    2,TSEDER(NY,NX),DORGC(L,NY,NX),DORGE(NY,NX)
C    4,VOLTI(L,NY,NX),VOLT(L,NY,NX),VOLA(L,NY,NX),VOLW(L,NY,NX)
C    4,VOLI(L,NY,NX),VOLP(L,NY,NX),DVOLW(L,NY,NX),DVOLI(L,NY,NX)
C    5,DLYR(3,L,NY,NX),CDPTHY(L,NY,NX),CDPTHX(L,NY,NX)
C    4,VOLI(L,NY,NX),VOLP(L,NY,NX),DVOLW(L,NY,NX)
C    5,VOLW(L,NY,NX)+VOLI(L,NY,NX)
C    5,BKDS(L,NY,NX),BKDSI(L,NY,NX)
C    6,DTBLZ(NY,NX),DTBLX(NY,NX)
C    7,VOLAH(L,NY,NX),VOLWH(L,NY,NX),VOLIH(L,NY,NX)
C    8,ORGCX(L,NY,NX),ORGC(L,NY,NX)
1114  FORMAT(A8,10I4,40E16.8)
C     ENDIF
C
C     RESET POND SURFACE LAYER NUMBER IF SURFACE LAYER
C     LOST TO EVAPORATION
C
C     VHCP,VHCPNX=current, minimum soil heat capacity
C     VOLT,VOLX=total, micropore soil volume
C     DDLYRX=total change in DLYR from all disturbances used to
C        calculate redistribution of soil material
C     IFLGL=pond boundary flag:layer underneath is 1:pond, 2:soil
C     DLYR=layer depth
C     BKDS=bulk density (0=water,>0=soil)
C     FMPR=1.0-(coarse fragment+macropore) fraction
C
      ELSEIF(NN.EQ.2)THEN
      IF((L.EQ.NU(NY,NX).AND.BKDS(NU(NY,NX),NY,NX).LE.ZERO)
     2.AND.(VHCP(NU(NY,NX),NY,NX).LE.VHCPNX(NY,NX)
     3.OR.NUM(NY,NX).GT.NU(NY,NX)))THEN
      NUX=NU(NY,NX)
      DO 9970 LL=NUX+1,NL(NY,NX)
      IF(VOLX(LL,NY,NX).GT.ZEROS2(NY,NX))THEN
      NU(NY,NX)=LL
      DDLYRX(NN)=DLYR(3,NUX,NY,NX)
      IFLGL(L,NN)=1
      DLYR(3,NUX,NY,NX)=0.0
      IF(BKDS(NUX,NY,NX).LE.ZERO)THEN
      VOLT(NUX,NY,NX)=AREA(3,NUX,NY,NX)*DLYR(3,NUX,NY,NX)
      VOLX(NUX,NY,NX)=VOLT(NUX,NY,NX)*FMPR(NUX,NY,NX)
      ENDIF
C     WRITE(*,5598)'SURFX',I,J,NX,NY,L,LL,NUX,NU(NY,NX),NUM(NY,NX)
C    2,DDLYRX(NN),VOLX(LL,NY,NX),VOLW(LL,NY,NX),VOLI(LL,NY,NX)
C    2,VHCP(LL,NY,NX),VHCPNX(NY,NX)
5598  FORMAT(A8,9I4,12E16.8)
      GO TO 9971
      ENDIF
9970  CONTINUE
      ELSE
      DDLYRX(NN)=0.0
      IFLGL(L,NN)=0
      ENDIF
9971  CONTINUE
C
C     RESET POND SURFACE LAYER NUMBER IF SURFACE LAYER
C     RESTORED BY PRECIPITATION
C
C     VOLW,VOLI=soil surface layer water,ice content
C     VOLWD=soil surface water ponding capacity
C     BKDS=bulk density (0=water,>0=soil)
C     DDLYRX=total change in DLYR from all disturbances used to
C        calculate redistribution of soil material
C     IFLGL=pond boundary flag:layer underneath is 1:pond, 2:soil
C     DLYR=layer depth: 0=pond, NU=soil surface layer
C     VOLR=dry litter volume
C     VOLWRX=litter water holding capacity
C     CDPTH=cumulative depth to bottom of soil layer
C     CDPTHY=CDPTH excluding freeze-thaw effects
C     CDPTHZ=cumulative depth to bottom of soil layer from current
C        soil surface depth
C     DPTH=depth to midpoint of soil layer
C
      ELSEIF(NN.EQ.3)THEN
      XVOLWP=AMAX1(0.0,VOLW(0,NY,NX)-VOLWD(NY,NX))
      IF(L.EQ.NU(NY,NX).AND.CDPTH(0,NY,NX).GT.CDPTHI(NY,NX)
     2.AND.XVOLWP.GT.VOLWD(NY,NX)+VHCPNX(NY,NX)/4.19)THEN
C     IF((BKDS(L,NY,NX).GT.ZERO.AND.NU(NY,NX).GT.NUI(NY,NX))
C    2.OR.(BKDS(L,NY,NX).LE.ZERO))THEN
      IF(BKDS(L,NY,NX).GT.ZERO.AND.NU(NY,NX).GT.NUI(NY,NX))THEN
      NU(NY,NX)=NUI(NY,NX)
      NUM(NY,NX)=NUI(NY,NX)
      DDLYRX(NN)=(VOLWD(NY,NX)-XVOLWP)/AREA(3,0,NY,NX)
      IFLGL(L,NN)=1
      DLYR0=(AMAX1(0.0,VOLW(0,NY,NX)+VOLI(0,NY,NX)-VOLWRX(NY,NX))
     2+VOLR(NY,NX))/AREA(3,0,NY,NX)
      DLYR(3,0,NY,NX)=DLYR0+DDLYRX(NN)
      DLYR(3,NU(NY,NX),NY,NX)=DLYR(3,NU(NY,NX),NY,NX)-DDLYRX(NN)
      IF(L.GT.2)THEN
      DO 260 LL=L-2,NU(NY,NX),-1
      CDPTH(LL,NY,NX)=CDPTH(L-1,NY,NX)
      CDPTHY(LL,NY,NX)=CDPTHY(L-1,NY,NX)
260   CONTINUE
      ENDIF
      CDPTH(0,NY,NX)=CDPTH(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
      CDPTHY(0,NY,NX)=CDPTHY(NU(NY,NX),NY,NX)-DLYR(3,NU(NY,NX),NY,NX)
      DPTH(NU(NY,NX),NY,NX)=0.5*(CDPTH(NU(NY,NX),NY,NX)
     2+CDPTH(0,NY,NX))
      CDPTHZ(NU(NY,NX),NY,NX)=DLYR(3,NU(NY,NX),NY,NX)
      DPTHZ(NU(NY,NX),NY,NX)=0.5*CDPTHZ(NU(NY,NX),NY,NX)
C     WRITE(*,5597)'SURFY',I,J,NX,NY,L,NUI(NY,NX),NU(NY,NX)
C    2,NUM(NY,NX),IFLGL(L,NN),VOLWD(NY,NX),XVOLWP,DDLYRX(NN)
C    2,DLYR0,DLYR(3,0,NY,NX),DLYR(3,NU(NY,NX),NY,NX),ORGC(0,NY,NX)
C    2,VOLW(0,NY,NX),VOLX(0,NY,NX),VOLW(NU(NY,NX),NY,NX)
C    3,CDPTH(0,NY,NX),CDPTH(NU(NY,NX),NY,NX),VOLR(NY,NX)
5597  FORMAT(A8,9I4,20E16.8)
      ELSE
      DDLYRX(NN)=0.0
      IFLGL(L,NN)=0
      ENDIF
      ELSE
      DDLYRX(NN)=0.0
      IFLGL(L,NN)=0
      ENDIF
      ENDIF
C
C     TRANSFER STATE VARIABLES BETWEEN LAYERS
C
C     IFLGL=pond boundary flag:layer underneath is 1:pond, 2:soil
C     BKDS=bulk density (0=water,>0=soil)
C     DDLYRX=total change in DLYR from all disturbances used to
C        calculate redistribution of soil material
C     L0,L1=soil layer from, to which soil contents are redistributed
C     FX=fraction of soil contents redistributed from L0 to L1
C     VOLW,VOLI=soil water,ice content
C     AREA=grid cell surface area
C     DLYR=layer depth
C
C     IF(IFLGL(L,NN).EQ.1)THEN
      IF(ABS(DDLYRX(NN)).GT.ZERO)THEN
      IF(DDLYRX(NN).GT.ZERO)THEN
C
C     IF LAYER GETS DEEPER, MATERIALS IN LAYER MOVE UP
C
      IF(IFLGL(L,2).EQ.0)THEN
      L1=L
      L0=L+1
      ELSE
      L1=NU(NY,NX)
      L0=NUX
      ENDIF
      IF((BKDS(L,NY,NX).LE.ZERO.AND.IFLGL(L,1).EQ.2)
     2.OR.(DLYR(3,L0,NY,NX).LE.ZERO2.AND.IFLGL(L,6).EQ.1))THEN
      FX=1.0
      ELSE
      IF(BKDS(L0,NY,NX).LE.ZERO)THEN
      DPTWI=(VOLW(L0,NY,NX)+VOLI(L0,NY,NX))/AREA(3,L0,NY,NX)
      IF(DPTWI.GT.ZERO)THEN
      FX=AMIN1(1.0,DDLYRX(NN)/DPTWI)
      ELSE
      FX=1.0
      ENDIF
      ELSE
      IF(DLYR(3,L0,NY,NX).GT.ZERO)THEN
      FX=AMIN1(1.0,DDLYRX(NN)/DLYR(3,L0,NY,NX))
      ELSE
      FX=1.0
      ENDIF
      ENDIF
      ENDIF
      ELSE
C
C     IF LAYER GETS SHALLOWER, MATERIALS IN LAYER MOVE DOWN
C
      IF(IFLGL(L,3).EQ.0)THEN
      L1=L+1
      L0=L
      ELSE
      L1=NU(NY,NX)
      L0=0
      ENDIF
      IF(BKDS(L0,NY,NX).LE.ZERO)THEN
      DPTWI=(VOLW(L0,NY,NX)+VOLI(L0,NY,NX))/AREA(3,L0,NY,NX)
      IF(DPTWI.GT.ZERO)THEN
      FX=AMIN1(1.0,-DDLYRX(NN)/DPTWI)
      ELSE
      FX=1.0
      ENDIF
      ELSE
      IF(DLYR(3,L0,NY,NX).GT.ZERO)THEN
      FX=AMIN1(1.0,-DDLYRX(NN)/DLYR(3,L0,NY,NX))
      ELSE
      FX=1.0
      ENDIF
      ENDIF
      ENDIF
C     IF(NX.EQ.1)THEN
C     WRITE(*,5601)'FX',I,J,NX,NY,L,L0,L1,NN,FO,FX
C    2,DDLYRX(NN),VOLW(L0,NY,NX),BKDS(L0,NY,NX),DPTWI
C    3,DLYR(3,L1,NY,NX),DLYR(3,L0,NY,NX)
5601  FORMAT(A8,8I4,20E12.4)
C     ENDIF
C
C     REDISTRIBUTE POND, SOIL MATERIAL
C
C     FX=fraction of pond or soil contents redistributed from L0 to L1
C     IFLGS=disturbance flag
C     BKDS=bulk density (0=water,>0=soil)
C
      IF(FX.GT.ZERO)THEN
      IFLGS(NY,NX)=1
      FY=1.0-FX
C     IF(FY.LE.ZERO2)FY=0.0
C
C     REDISTRIBUTE POND MATERIAL
C
      IF(BKDS(L0,NY,NX).LE.ZERO)THEN
C     IF(I.EQ.178)THEN
C     WRITE(*,5599)'POND1',I,J,NFZ,NX,NY,L,L0,L1,NU(NY,NX),NN,FX,FY
C    2,DDLYRX(NN),VOLY(L0,NY,NX),VOLX(L0,NY,NX),VOLW(L0,NY,NX)
C    3,VOLI(L0,NY,NX),VOLY(L1,NY,NX),VOLX(L1,NY,NX),VOLW(L1,NY,NX)
C    4,VOLI(L1,NY,NX),CDPTH(L0,NY,NX),CDPTH(L1,NY,NX)
C    5,VOLW(L0,NY,NX)+VOLW(L1,NY,NX)
C    5,(OSC(M,4,L0,NY,NX),M=1,5),(OSC(M,4,L1,NY,NX),M=1,5)
C    5,(OQC(K,L0,NY,NX),K=0,4),(OQC(K,L1,NY,NX),K=0,4)
C    6,DLYR(3,L0,NY,NX),DLYR(3,L1,NY,NX)
C    5,TKS(L0,NY,NX),TKS(L1,NY,NX)
C    5,VHCP(L0,NY,NX),VHCP(L1,NY,NX)
C    6,ZNH4S(L0,NY,NX),ZNH4B(L0,NY,NX),ZNH3S(L0,NY,NX),ZNH3B(L0,NY,NX)
C    6,ZNH4S(L1,NY,NX),ZNH4B(L1,NY,NX),ZNH3S(L1,NY,NX),ZNH3B(L1,NY,NX)
C    6,(WTRT1(1,L1,NR,1,NY,NX),NR=1,NRT(1,NY,NX))
C    6,(WTRT2(1,L1,NR,1,NY,NX),NR=1,NRT(1,NY,NX))
C    6,WTRTL(1,L1,1,NY,NX),WTNDL(L1,1,NY,NX)
C    6,CPOOLR(1,L1,1,NY,NX),ZPOOLR(1,L1,1,NY,NX)
C    6,OXYA(1,L1,1,NY,NX),OXYP(1,L1,1,NY,NX)
C    6,(WTRT1(1,L0,NR,1,NY,NX),NR=1,NRT(1,NY,NX))
C    6,(WTRT2(1,L0,NR,1,NY,NX),NR=1,NRT(1,NY,NX))
C    6,WTRTL(1,L0,1,NY,NX),WTNDL(L0,1,NY,NX)
C    6,CPOOLR(1,L0,1,NY,NX),ZPOOLR(1,L0,1,NY,NX)
C    6,OXYA(1,L0,1,NY,NX),OXYP(1,L0,1,NY,NX)
5599  FORMAT(A8,10I4,60E16.8)
C     ENDIF
C
C     REDISTRIBUTE POND CONTENTS TO DESTINATION LAYER
C
C     POND WATER, ICE, HEAT
C
      VOLW(L1,NY,NX)=VOLW(L1,NY,NX)
     2+FX*VOLW(L0,NY,NX)
      VOLI(L1,NY,NX)=VOLI(L1,NY,NX)
     2+FX*VOLI(L0,NY,NX)
      VOLP(L1,NY,NX)=VOLP(L1,NY,NX)
     2+FX*VOLP(L0,NY,NX)
      VOLA(L1,NY,NX)=VOLA(L1,NY,NX)
     2+FX*VOLA(L0,NY,NX)
      VOLY(L1,NY,NX)=VOLY(L1,NY,NX)
     2+FX*VOLY(L0,NY,NX)
      VOLWX(L1,NY,NX)=VOLW(L1,NY,NX)
      ENGY1=VHCP(L1,NY,NX)*TKS(L1,NY,NX)
      ENGY0=VHCP(L0,NY,NX)*TKS(L0,NY,NX)
      ENGY1=ENGY1+FX*ENGY0
      VHCM(L1,NY,NX)=VHCM(L1,NY,NX)
     2+FX*VHCM(L0,NY,NX)
      VHCP(L1,NY,NX)=VHCM(L1,NY,NX)
     2+4.19*(VOLW(L1,NY,NX)+VOLWH(L1,NY,NX))
     3+1.9274*(VOLI(L1,NY,NX)+VOLIH(L1,NY,NX))
      IF(VHCP(L1,NY,NX).GT.VHCPNX(NY,NX))THEN
      TKS(L1,NY,NX)=ENGY1/VHCP(L1,NY,NX)
      ELSE
      TKS(L1,NY,NX)=TKS(L0,NY,NX)
      ENDIF
      TCS(L1,NY,NX)=TKS(L1,NY,NX)-273.15
C
C     POND N,P SOLUTES IN BAND
C
      ZNH4S(L1,NY,NX)=ZNH4S(L1,NY,NX)
     2+FX*ZNH4S(L0,NY,NX)
      ZNH4B(L1,NY,NX)=ZNH4B(L1,NY,NX)
     2+FX*ZNH4B(L0,NY,NX)
      ZNH3S(L1,NY,NX)=ZNH3S(L1,NY,NX)
     2+FX*ZNH3S(L0,NY,NX)
      ZNH3B(L1,NY,NX)=ZNH3B(L1,NY,NX)
     2+FX*ZNH3B(L0,NY,NX)
      ZNO3S(L1,NY,NX)=ZNO3S(L1,NY,NX)
     2+FX*ZNO3S(L0,NY,NX)
      ZNO3B(L1,NY,NX)=ZNO3B(L1,NY,NX)
     2+FX*ZNO3B(L0,NY,NX)
      ZNO2S(L1,NY,NX)=ZNO2S(L1,NY,NX)
     2+FX*ZNO2S(L0,NY,NX)
      ZNO2B(L1,NY,NX)=ZNO2B(L1,NY,NX)
     2+FX*ZNO2B(L0,NY,NX)
      H1PO4(L1,NY,NX)=H1PO4(L1,NY,NX)
     2+FX*H1PO4(L0,NY,NX)
      H2PO4(L1,NY,NX)=H2PO4(L1,NY,NX)
     2+FX*H2PO4(L0,NY,NX)
C
C     POND SALTS IN NON-BAND
C
      IF(ISALTG.NE.0)THEN
      ZAL(L1,NY,NX)=ZAL(L1,NY,NX)
     2+FX*ZAL(L0,NY,NX)
      ZFE(L1,NY,NX)=ZFE(L1,NY,NX)
     2+FX*ZFE(L0,NY,NX)
      ZHY(L1,NY,NX)=ZHY(L1,NY,NX)
     2+FX*ZHY(L0,NY,NX)
      ZCA(L1,NY,NX)=ZCA(L1,NY,NX)
     2+FX*ZCA(L0,NY,NX)
      ZMG(L1,NY,NX)=ZMG(L1,NY,NX)
     2+FX*ZMG(L0,NY,NX)
      ZNA(L1,NY,NX)=ZNA(L1,NY,NX)
     2+FX*ZNA(L0,NY,NX)
      ZKA(L1,NY,NX)=ZKA(L1,NY,NX)
     2+FX*ZKA(L0,NY,NX)
      ZOH(L1,NY,NX)=ZOH(L1,NY,NX)
     2+FX*ZOH(L0,NY,NX)
      ZSO4(L1,NY,NX)=ZSO4(L1,NY,NX)
     2+FX*ZSO4(L0,NY,NX)
      ZCL(L1,NY,NX)=ZCL(L1,NY,NX)
     2+FX*ZCL(L0,NY,NX)
      ZCO3(L1,NY,NX)=ZCO3(L1,NY,NX)
     2+FX*ZCO3(L0,NY,NX)
      ZHCO3(L1,NY,NX)=ZHCO3(L1,NY,NX)
     2+FX*ZHCO3(L0,NY,NX)
      ZALOH1(L1,NY,NX)=ZALOH1(L1,NY,NX)
     2+FX*ZALOH1(L0,NY,NX)
      ZALOH2(L1,NY,NX)=ZALOH2(L1,NY,NX)
     2+FX*ZALOH2(L0,NY,NX)
      ZALOH3(L1,NY,NX)=ZALOH3(L1,NY,NX)
     2+FX*ZALOH3(L0,NY,NX)
      ZALOH4(L1,NY,NX)=ZALOH4(L1,NY,NX)
     2+FX*ZALOH4(L0,NY,NX)
      ZALS(L1,NY,NX)=ZALS(L1,NY,NX)
     2+FX*ZALS(L0,NY,NX)
      ZFEOH1(L1,NY,NX)=ZFEOH1(L1,NY,NX)
     2+FX*ZFEOH1(L0,NY,NX)
      ZFEOH2(L1,NY,NX)=ZFEOH2(L1,NY,NX)
     2+FX*ZFEOH2(L0,NY,NX)
      ZFEOH3(L1,NY,NX)=ZFEOH3(L1,NY,NX)
     2+FX*ZFEOH3(L0,NY,NX)
      ZFEOH4(L1,NY,NX)=ZFEOH4(L1,NY,NX)
     2+FX*ZFEOH4(L0,NY,NX)
      ZFES(L1,NY,NX)=ZFES(L1,NY,NX)
     2+FX*ZFES(L0,NY,NX)
      ZCAO(L1,NY,NX)=ZCAO(L1,NY,NX)
     2+FX*ZCAO(L0,NY,NX)
      ZCAC(L1,NY,NX)=ZCAC(L1,NY,NX)
     2+FX*ZCAC(L0,NY,NX)
      ZCAH(L1,NY,NX)=ZCAH(L1,NY,NX)
     2+FX*ZCAH(L0,NY,NX)
      ZCAS(L1,NY,NX)=ZCAS(L1,NY,NX)
     2+FX*ZCAS(L0,NY,NX)
      ZMGO(L1,NY,NX)=ZMGO(L1,NY,NX)
     2+FX*ZMGO(L0,NY,NX)
      ZMGC(L1,NY,NX)=ZMGC(L1,NY,NX)
     2+FX*ZMGC(L0,NY,NX)
      ZMGH(L1,NY,NX)=ZMGH(L1,NY,NX)
     2+FX*ZMGH(L0,NY,NX)
      ZMGS(L1,NY,NX)=ZMGS(L1,NY,NX)
     2+FX*ZMGS(L0,NY,NX)
      ZNAC(L1,NY,NX)=ZNAC(L1,NY,NX)
     2+FX*ZNAC(L0,NY,NX)
      ZNAS(L1,NY,NX)=ZNAS(L1,NY,NX)
     2+FX*ZNAS(L0,NY,NX)
      ZKAS(L1,NY,NX)=ZKAS(L1,NY,NX)
     2+FX*ZKAS(L0,NY,NX)
      H0PO4(L1,NY,NX)=H0PO4(L1,NY,NX)
     2+FX*H0PO4(L0,NY,NX)
      H3PO4(L1,NY,NX)=H3PO4(L1,NY,NX)
     2+FX*H3PO4(L0,NY,NX)
      ZFE1P(L1,NY,NX)=ZFE1P(L1,NY,NX)
     2+FX*ZFE1P(L0,NY,NX)
      ZFE2P(L1,NY,NX)=ZFE2P(L1,NY,NX)
     2+FX*ZFE2P(L0,NY,NX)
      ZCA0P(L1,NY,NX)=ZCA0P(L1,NY,NX)
     2+FX*ZCA0P(L0,NY,NX)
      ZCA1P(L1,NY,NX)=ZCA1P(L1,NY,NX)
     2+FX*ZCA1P(L0,NY,NX)
      ZCA2P(L1,NY,NX)=ZCA2P(L1,NY,NX)
     2+FX*ZCA2P(L0,NY,NX)
      ZMG1P(L1,NY,NX)=ZMG1P(L1,NY,NX)
     2+FX*ZMG1P(L0,NY,NX)
      ENDIF
C
C     POND SALTS IN BAND
C
      IF(L0.NE.0)THEN
      H1POB(L1,NY,NX)=H1POB(L1,NY,NX)
     2+FX*H1POB(L0,NY,NX)
      H2POB(L1,NY,NX)=H2POB(L1,NY,NX)
     2+FX*H2POB(L0,NY,NX)
      IF(ISALTG.NE.0)THEN
      H0POB(L1,NY,NX)=H0POB(L1,NY,NX)
     2+FX*H0POB(L0,NY,NX)
      H3POB(L1,NY,NX)=H3POB(L1,NY,NX)
     2+FX*H3POB(L0,NY,NX)
      ZFE1PB(L1,NY,NX)=ZFE1PB(L1,NY,NX)
     2+FX*ZFE1PB(L0,NY,NX)
      ZFE2PB(L1,NY,NX)=ZFE2PB(L1,NY,NX)
     2+FX*ZFE2PB(L0,NY,NX)
      ZCA0PB(L1,NY,NX)=ZCA0PB(L1,NY,NX)
     2+FX*ZCA0PB(L0,NY,NX)
      ZCA1PB(L1,NY,NX)=ZCA1PB(L1,NY,NX)
     2+FX*ZCA1PB(L0,NY,NX)
      ZCA2PB(L1,NY,NX)=ZCA2PB(L1,NY,NX)
     2+FX*ZCA2PB(L0,NY,NX)
      ZMG1PB(L1,NY,NX)=ZMG1PB(L1,NY,NX)
     2+FX*ZMG1PB(L0,NY,NX)
      ENDIF
C
C     POND GASEOUS GASES
C
      CO2G(L1,NY,NX)=CO2G(L1,NY,NX)
     2+FX*CO2G(L0,NY,NX)
      CH4G(L1,NY,NX)=CH4G(L1,NY,NX)
     2+FX*CH4G(L0,NY,NX)
      OXYG(L1,NY,NX)=OXYG(L1,NY,NX)
     2+FX*OXYG(L0,NY,NX)
      Z2GG(L1,NY,NX)=Z2GG(L1,NY,NX)
     2+FX*Z2GG(L0,NY,NX)
      Z2OG(L1,NY,NX)=Z2OG(L1,NY,NX)
     2+FX*Z2OG(L0,NY,NX)
      ZNH3G(L1,NY,NX)=ZNH3G(L1,NY,NX)
     2+FX*ZNH3G(L0,NY,NX)
      H2GG(L1,NY,NX)=H2GG(L1,NY,NX)
     2+FX*H2GG(L0,NY,NX)
      ENDIF
C
C     POND AQUEOUS GASES
C
      CO2S(L1,NY,NX)=CO2S(L1,NY,NX)
     2+FX*CO2S(L0,NY,NX)
      CH4S(L1,NY,NX)=CH4S(L1,NY,NX)
     2+FX*CH4S(L0,NY,NX)
      OXYS(L1,NY,NX)=OXYS(L1,NY,NX)
     2+FX*OXYS(L0,NY,NX)
      Z2GS(L1,NY,NX)=Z2GS(L1,NY,NX)
     2+FX*Z2GS(L0,NY,NX)
      Z2OS(L1,NY,NX)=Z2OS(L1,NY,NX)
     2+FX*Z2OS(L0,NY,NX)
      H2GS(L1,NY,NX)=H2GS(L1,NY,NX)
     2+FX*H2GS(L0,NY,NX)
C
C     POND ORGANIC MATTER
C
      IF(IFLGL(L,3).EQ.0)THEN
      DO 7780 K=0,4
      OQC(K,L1,NY,NX)=OQC(K,L1,NY,NX)
     2+FX*OQC(K,L0,NY,NX)
      OQN(K,L1,NY,NX)=OQN(K,L1,NY,NX)
     2+FX*OQN(K,L0,NY,NX)
      OQP(K,L1,NY,NX)=OQP(K,L1,NY,NX)
     2+FX*OQP(K,L0,NY,NX)
      OQA(K,L1,NY,NX)=OQA(K,L1,NY,NX)
     2+FX*OQA(K,L0,NY,NX)
      OQCH(K,L1,NY,NX)=OQCH(K,L1,NY,NX)
     2+FX*OQCH(K,L0,NY,NX)
      OQNH(K,L1,NY,NX)=OQNH(K,L1,NY,NX)
     2+FX*OQNH(K,L0,NY,NX)
      OQPH(K,L1,NY,NX)=OQPH(K,L1,NY,NX)
     2+FX*OQPH(K,L0,NY,NX)
      OQAH(K,L1,NY,NX)=OQAH(K,L1,NY,NX)
     2+FX*OQAH(K,L0,NY,NX)
7780  CONTINUE
      ENDIF
C
C     POND PLANT MATTER
C
      IF(L0.NE.0)THEN
      DO 8900 NZ=1,NP(NY,NX)
      IF(WTRTL(1,L0,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WTRTL(1,L1,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      DO 8895 N=1,MY(NZ,NY,NX)
C
C     POND ROOT GASES
C
      CO2A(N,L1,NZ,NY,NX)=CO2A(N,L1,NZ,NY,NX)
     2+FX*CO2A(N,L0,NZ,NY,NX)
      OXYA(N,L1,NZ,NY,NX)=OXYA(N,L1,NZ,NY,NX)
     2+FX*OXYA(N,L0,NZ,NY,NX)
      CH4A(N,L1,NZ,NY,NX)=CH4A(N,L1,NZ,NY,NX)
     2+FX*CH4A(N,L0,NZ,NY,NX)
      Z2OA(N,L1,NZ,NY,NX)=Z2OA(N,L1,NZ,NY,NX)
     2+FX*Z2OA(N,L0,NZ,NY,NX)
      ZH3A(N,L1,NZ,NY,NX)=ZH3A(N,L1,NZ,NY,NX)
     2+FX*ZH3A(N,L0,NZ,NY,NX)
      H2GA(N,L1,NZ,NY,NX)=H2GA(N,L1,NZ,NY,NX)
     2+FX*H2GA(N,L0,NZ,NY,NX)
      CO2P(N,L1,NZ,NY,NX)=CO2P(N,L1,NZ,NY,NX)
     2+FX*CO2P(N,L0,NZ,NY,NX)
      OXYP(N,L1,NZ,NY,NX)=OXYP(N,L1,NZ,NY,NX)
     2+FX*OXYP(N,L0,NZ,NY,NX)
      CH4P(N,L1,NZ,NY,NX)=CH4P(N,L1,NZ,NY,NX)
     2+FX*CH4P(N,L0,NZ,NY,NX)
      Z2OP(N,L1,NZ,NY,NX)=Z2OP(N,L1,NZ,NY,NX)
     2+FX*Z2OP(N,L0,NZ,NY,NX)
      ZH3P(N,L1,NZ,NY,NX)=ZH3P(N,L1,NZ,NY,NX)
     2+FX*ZH3P(N,L0,NZ,NY,NX)
      H2GP(N,L1,NZ,NY,NX)=H2GP(N,L1,NZ,NY,NX)
     2+FX*H2GP(N,L0,NZ,NY,NX)
C
C     POND ROOTS
C
      DO 8870 NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L1,NR,NZ,NY,NX)=WTRT1(N,L1,NR,NZ,NY,NX)
     2+FX*WTRT1(N,L0,NR,NZ,NY,NX)
      WTRT1N(N,L1,NR,NZ,NY,NX)=WTRT1N(N,L1,NR,NZ,NY,NX)
     2+FX*WTRT1N(N,L0,NR,NZ,NY,NX)
      WTRT1P(N,L1,NR,NZ,NY,NX)=WTRT1P(N,L1,NR,NZ,NY,NX)
     2+FX*WTRT1P(N,L0,NR,NZ,NY,NX)
      WTRT2(N,L1,NR,NZ,NY,NX)=WTRT2(N,L1,NR,NZ,NY,NX)
     2+FX*WTRT2(N,L0,NR,NZ,NY,NX)
      WTRT2N(N,L1,NR,NZ,NY,NX)=WTRT2N(N,L1,NR,NZ,NY,NX)
     2+FX*WTRT2N(N,L0,NR,NZ,NY,NX)
      WTRT2P(N,L1,NR,NZ,NY,NX)=WTRT2P(N,L1,NR,NZ,NY,NX)
     2+FX*WTRT2P(N,L0,NR,NZ,NY,NX)
      RTLG1(N,L1,NR,NZ,NY,NX)=RTLG1(N,L1,NR,NZ,NY,NX)
     2+FX*RTLG1(N,L0,NR,NZ,NY,NX)
      RTLG2(N,L1,NR,NZ,NY,NX)=RTLG2(N,L1,NR,NZ,NY,NX)
     2+FX*RTLG2(N,L0,NR,NZ,NY,NX)
      RTN2(N,L1,NR,NZ,NY,NX)=RTN2(N,L1,NR,NZ,NY,NX)
     2+FX*RTN2(N,L0,NR,NZ,NY,NX)
8870  CONTINUE
      CPOOLR(N,L1,NZ,NY,NX)=CPOOLR(N,L1,NZ,NY,NX)
     2+FX*CPOOLR(N,L0,NZ,NY,NX)
      ZPOOLR(N,L1,NZ,NY,NX)=ZPOOLR(N,L1,NZ,NY,NX)
     2+FX*ZPOOLR(N,L0,NZ,NY,NX)
      PPOOLR(N,L1,NZ,NY,NX)=PPOOLR(N,L1,NZ,NY,NX)
     2+FX*PPOOLR(N,L0,NZ,NY,NX)
      WTRTL(N,L1,NZ,NY,NX)=WTRTL(N,L1,NZ,NY,NX)
     2+FX*WTRTL(N,L0,NZ,NY,NX)
      WTRTD(N,L1,NZ,NY,NX)=WTRTD(N,L1,NZ,NY,NX)
     2+FX*WTRTD(N,L0,NZ,NY,NX)
      WSRTL(N,L1,NZ,NY,NX)=WSRTL(N,L1,NZ,NY,NX)
     2+FX*WSRTL(N,L0,NZ,NY,NX)
      RTN1(N,L1,NZ,NY,NX)=RTN1(N,L1,NZ,NY,NX)
     2+FX*RTN1(N,L0,NZ,NY,NX)
      RTNL(N,L1,NZ,NY,NX)=RTNL(N,L1,NZ,NY,NX)
     2+FX*RTNL(N,L0,NZ,NY,NX)
      RTLGP(N,L1,NZ,NY,NX)=RTLGP(N,L1,NZ,NY,NX)
     2+FX*RTLGP(N,L0,NZ,NY,NX)
      RTDNP(N,L1,NZ,NY,NX)=RTDNP(N,L1,NZ,NY,NX)
     2+FX*RTDNP(N,L0,NZ,NY,NX)
      RTVLP(N,L1,NZ,NY,NX)=RTVLP(N,L1,NZ,NY,NX)
     2+FX*RTVLP(N,L0,NZ,NY,NX)
      RTVLW(N,L1,NZ,NY,NX)=RTVLW(N,L1,NZ,NY,NX)
     2+FX*RTVLW(N,L0,NZ,NY,NX)
      RRAD1(N,L1,NZ,NY,NX)=RRAD1(N,L1,NZ,NY,NX)
     2+FX*RRAD1(N,L0,NZ,NY,NX)
      RRAD2(N,L1,NZ,NY,NX)=RRAD2(N,L1,NZ,NY,NX)
     2+FX*RRAD2(N,L0,NZ,NY,NX)
      RTARP(N,L1,NZ,NY,NX)=RTARP(N,L1,NZ,NY,NX)
     2+FX*RTARP(N,L0,NZ,NY,NX)
      RTLGA(N,L1,NZ,NY,NX)=RTLGA(N,L1,NZ,NY,NX)
     2+FX*RTLGA(N,L0,NZ,NY,NX)
8895  CONTINUE
C
C     POND ROOT NODULES
C
      WTNDL(L1,NZ,NY,NX)=WTNDL(L1,NZ,NY,NX)
     2+FX*WTNDL(L0,NZ,NY,NX)
      WTNDLN(L1,NZ,NY,NX)=WTNDLN(L1,NZ,NY,NX)
     2+FX*WTNDLN(L0,NZ,NY,NX)
      WTNDLP(L1,NZ,NY,NX)=WTNDLP(L1,NZ,NY,NX)
     2+FX*WTNDLP(L0,NZ,NY,NX)
      CPOOLN(L1,NZ,NY,NX)=CPOOLN(L1,NZ,NY,NX)
     2+FX*CPOOLN(L0,NZ,NY,NX)
      ZPOOLN(L1,NZ,NY,NX)=ZPOOLN(L1,NZ,NY,NX)
     2+FX*ZPOOLN(L0,NZ,NY,NX)
      PPOOLN(L1,NZ,NY,NX)=PPOOLN(L1,NZ,NY,NX)
     2+FX*PPOOLN(L0,NZ,NY,NX)
      ENDIF
8900  CONTINUE
      ENDIF
C
C     POND SEDIMENT
C
      IF(FX.EQ.1.0)THEN
      IF(L0.NE.0)THEN
      SAND(L1,NY,NX)=SAND(L1,NY,NX)
     2+FX*SAND(L0,NY,NX)
      SILT(L1,NY,NX)=SILT(L1,NY,NX)
     2+FX*SILT(L0,NY,NX)
      CLAY(L1,NY,NX)=CLAY(L1,NY,NX)
     2+FX*CLAY(L0,NY,NX)
      XCEC(L1,NY,NX)=XCEC(L1,NY,NX)
     2+FX*XCEC(L0,NY,NX)
      XAEC(L1,NY,NX)=XAEC(L1,NY,NX)
     2+FX*XAEC(L0,NY,NX)
      ENDIF
C
C     POND FERTILIZER IN BAND, NON-BAND
C
      ZNH4FA(L1,NY,NX)=ZNH4FA(L1,NY,NX)
     2+FX*ZNH4FA(L0,NY,NX)
      ZNH3FA(L1,NY,NX)=ZNH3FA(L1,NY,NX)
     2+FX*ZNH3FA(L0,NY,NX)
      ZNHUFA(L1,NY,NX)=ZNHUFA(L1,NY,NX)
     2+FX*ZNHUFA(L0,NY,NX)
      ZNO3FA(L1,NY,NX)=ZNO3FA(L1,NY,NX)
     2+FX*ZNO3FA(L0,NY,NX)
      ZNH4FB(L1,NY,NX)=ZNH4FB(L1,NY,NX)
     2+FX*ZNH4FB(L0,NY,NX)
      ZNH3FB(L1,NY,NX)=ZNH3FB(L1,NY,NX)
     2+FX*ZNH3FB(L0,NY,NX)
      ZNHUFB(L1,NY,NX)=ZNHUFB(L1,NY,NX)
     2+FX*ZNHUFB(L0,NY,NX)
      ZNO3FB(L1,NY,NX)=ZNO3FB(L1,NY,NX)
     2+FX*ZNO3FB(L0,NY,NX)
C
C     POND ADSORBED CATIONS IN BAND, NON-BAND
C
      XN4(L1,NY,NX)=XN4(L1,NY,NX)
     2+FX*XN4(L0,NY,NX)
      XNB(L1,NY,NX)=XNB(L1,NY,NX)
     2+FX*XNB(L0,NY,NX)
      XHY(L1,NY,NX)=XHY(L1,NY,NX)
     2+FX*XHY(L0,NY,NX)
      XAL(L1,NY,NX)=XAL(L1,NY,NX)
     2+FX*XAL(L0,NY,NX)
      XFE(L1,NY,NX)=XFE(L1,NY,NX)
     2+FX*XFE(L0,NY,NX)
      XCA(L1,NY,NX)=XCA(L1,NY,NX)
     2+FX*XCA(L0,NY,NX)
      XMG(L1,NY,NX)=XMG(L1,NY,NX)
     2+FX*XMG(L0,NY,NX)
      XNA(L1,NY,NX)=XNA(L1,NY,NX)
     2+FX*XNA(L0,NY,NX)
      XKA(L1,NY,NX)=XKA(L1,NY,NX)
     2+FX*XKA(L0,NY,NX)
      XHC(L1,NY,NX)=XHC(L1,NY,NX)
     2+FX*XHC(L0,NY,NX)
      XALO2(L1,NY,NX)=XALO2(L1,NY,NX)
     2+FX*XALO2(L0,NY,NX)
      XFEO2(L1,NY,NX)=XFEO2(L1,NY,NX)
     2+FX*XFEO2(L0,NY,NX)
C
C     POND ADSORBED ANIONS IN BAND, NON-BAND
C
      XOH0(L1,NY,NX)=XOH0(L1,NY,NX)
     2+FX*XOH0(L0,NY,NX)
      XOH1(L1,NY,NX)=XOH1(L1,NY,NX)
     2+FX*XOH1(L0,NY,NX)
      XOH2(L1,NY,NX)=XOH2(L1,NY,NX)
     2+FX*XOH2(L0,NY,NX)
      XH1P(L1,NY,NX)=XH1P(L1,NY,NX)
     2+FX*XH1P(L0,NY,NX)
      XH2P(L1,NY,NX)=XH2P(L1,NY,NX)
     2+FX*XH2P(L0,NY,NX)
      XOH0B(L1,NY,NX)=XOH0B(L1,NY,NX)
     2+FX*XOH0B(L0,NY,NX)
      XOH1B(L1,NY,NX)=XOH1B(L1,NY,NX)
     2+FX*XOH1B(L0,NY,NX)
      XOH2B(L1,NY,NX)=XOH2B(L1,NY,NX)
     2+FX*XOH2B(L0,NY,NX)
      XH1PB(L1,NY,NX)=XH1PB(L1,NY,NX)
     2+FX*XH1PB(L0,NY,NX)
      XH2PB(L1,NY,NX)=XH2PB(L1,NY,NX)
     2+FX*XH2PB(L0,NY,NX)
C
C     POND PRECIPITATES IN BAND, NON-BAND
C
      PALOH(L1,NY,NX)=PALOH(L1,NY,NX)
     2+FX*PALOH(L0,NY,NX)
      PFEOH(L1,NY,NX)=PFEOH(L1,NY,NX)
     2+FX*PFEOH(L0,NY,NX)
      PCACO(L1,NY,NX)=PCACO(L1,NY,NX)
     2+FX*PCACO(L0,NY,NX)
      PCASO(L1,NY,NX)=PCASO(L1,NY,NX)
     2+FX*PCASO(L0,NY,NX)
      PALPO(L1,NY,NX)=PALPO(L1,NY,NX)
     2+FX*PALPO(L0,NY,NX)
      PFEPO(L1,NY,NX)=PFEPO(L1,NY,NX)
     2+FX*PFEPO(L0,NY,NX)
      PCAPD(L1,NY,NX)=PCAPD(L1,NY,NX)
     2+FX*PCAPD(L0,NY,NX)
      PCAPH(L1,NY,NX)=PCAPH(L1,NY,NX)
     2+FX*PCAPH(L0,NY,NX)
      PCAPM(L1,NY,NX)=PCAPM(L1,NY,NX)
     2+FX*PCAPM(L0,NY,NX)
      PALPB(L1,NY,NX)=PALPB(L1,NY,NX)
     2+FX*PALPB(L0,NY,NX)
      PFEPB(L1,NY,NX)=PFEPB(L1,NY,NX)
     2+FX*PFEPB(L0,NY,NX)
      PCPDB(L1,NY,NX)=PCPDB(L1,NY,NX)
     2+FX*PCPDB(L0,NY,NX)
      PCPHB(L1,NY,NX)=PCPHB(L1,NY,NX)
     2+FX*PCPHB(L0,NY,NX)
      PCPMB(L1,NY,NX)=PCPMB(L1,NY,NX)
     2+FX*PCPMB(L0,NY,NX)
C
C     POND ORGANIC MATTER
C
      IF(IFLGL(L,3).EQ.0)THEN
      ORGC(L1,NY,NX)=ORGC(L1,NY,NX)+FX*ORGC(L0,NY,NX)
      DO 7965 K=0,5
      DO 7965 N=1,7
      DO 7965 M=1,3
      OMC(M,N,K,L1,NY,NX)=OMC(M,N,K,L1,NY,NX)
     2+FX*OMC(M,N,K,L0,NY,NX)
      OMN(M,N,K,L1,NY,NX)=OMN(M,N,K,L1,NY,NX)
     2+FX*OMN(M,N,K,L0,NY,NX)
      OMP(M,N,K,L1,NY,NX)=OMP(M,N,K,L1,NY,NX)
     2+FX*OMP(M,N,K,L0,NY,NX)
7965  CONTINUE
      DO 7785 K=0,4
      DO 7775 M=1,2
      ORC(M,K,L1,NY,NX)=ORC(M,K,L1,NY,NX)
     2+FX*ORC(M,K,L0,NY,NX)
      ORN(M,K,L1,NY,NX)=ORN(M,K,L1,NY,NX)
     2+FX*ORN(M,K,L0,NY,NX)
      ORP(M,K,L1,NY,NX)=ORP(M,K,L1,NY,NX)
     2+FX*ORP(M,K,L0,NY,NX)
7775  CONTINUE
      OHC(K,L1,NY,NX)=OHC(K,L1,NY,NX)
     2+FX*OHC(K,L0,NY,NX)
      OHN(K,L1,NY,NX)=OHN(K,L1,NY,NX)
     2+FX*OHN(K,L0,NY,NX)
      OHP(K,L1,NY,NX)=OHP(K,L1,NY,NX)
     2+FX*OHP(K,L0,NY,NX)
      OHA(K,L1,NY,NX)=OHA(K,L1,NY,NX)
     2+FX*OHA(K,L0,NY,NX)
      DO 7770 M=1,5
      OSC(M,K,L1,NY,NX)=OSC(M,K,L1,NY,NX)
     2+FX*OSC(M,K,L0,NY,NX)
      OSA(M,K,L1,NY,NX)=OSA(M,K,L1,NY,NX)
     2+FX*OSA(M,K,L0,NY,NX)
      OSN(M,K,L1,NY,NX)=OSN(M,K,L1,NY,NX)
     2+FX*OSN(M,K,L0,NY,NX)
      OSP(M,K,L1,NY,NX)=OSP(M,K,L1,NY,NX)
     2+FX*OSP(M,K,L0,NY,NX)
7770  CONTINUE
7785  CONTINUE
      ENDIF
      ENDIF
C
C     REDISTRIBUTE POND CONTENTS FROM SOURCE LAYER
C
C     POND WATER, ICE, HEAT
C
      VOLW(L0,NY,NX)=FY*VOLW(L0,NY,NX)
      VOLI(L0,NY,NX)=FY*VOLI(L0,NY,NX)
      VOLP(L0,NY,NX)=FY*VOLP(L0,NY,NX)
      VOLA(L0,NY,NX)=FY*VOLA(L0,NY,NX)
      VOLY(L0,NY,NX)=FY*VOLY(L0,NY,NX)
      VOLWX(L0,NY,NX)=VOLW(L0,NY,NX)
      ENGY0=FY*ENGY0
      VHCM(L0,NY,NX)=FY*VHCM(L0,NY,NX)
      IF(L0.NE.0)THEN
      VHCP(L0,NY,NX)=VHCM(L0,NY,NX)
     2+4.19*(VOLW(L0,NY,NX)+VOLWH(L0,NY,NX))
     3+1.9274*(VOLI(L0,NY,NX)+VOLIH(L0,NY,NX))
      ELSE
      VHCP(L0,NY,NX)=VHCM(L0,NY,NX)
     2+4.19*VOLW(L0,NY,NX)+1.9274*VOLI(L0,NY,NX)
      ENDIF
      IF(VHCP(L0,NY,NX).GT.VHCPNX(NY,NX))THEN
      TKS(L0,NY,NX)=ENGY0/VHCP(L0,NY,NX)
      ELSE
      TKS(L0,NY,NX)=TKS(L1,NY,NX)
      ENDIF
      TCS(L0,NY,NX)=TKS(L0,NY,NX)-273.15
C
C     POND N,P SOLUTES IN BAND, NON-BAND
C
      ZNH4S(L0,NY,NX)=FY*ZNH4S(L0,NY,NX)
      ZNH4B(L0,NY,NX)=FY*ZNH4B(L0,NY,NX)
      ZNH3S(L0,NY,NX)=FY*ZNH3S(L0,NY,NX)
      ZNH3B(L0,NY,NX)=FY*ZNH3B(L0,NY,NX)
      ZNO3S(L0,NY,NX)=FY*ZNO3S(L0,NY,NX)
      ZNO3B(L0,NY,NX)=FY*ZNO3B(L0,NY,NX)
      ZNO2S(L0,NY,NX)=FY*ZNO2S(L0,NY,NX)
      ZNO2B(L0,NY,NX)=FY*ZNO2B(L0,NY,NX)
      H1PO4(L0,NY,NX)=FY*H1PO4(L0,NY,NX)
      H2PO4(L0,NY,NX)=FY*H2PO4(L0,NY,NX)
C
C     POND SALTS IN NON-BAND
C
      IF(ISALTG.NE.0)THEN
      ZAL(L0,NY,NX)=FY*ZAL(L0,NY,NX)
      ZFE(L0,NY,NX)=FY*ZFE(L0,NY,NX)
      ZHY(L0,NY,NX)=FY*ZHY(L0,NY,NX)
      ZCA(L0,NY,NX)=FY*ZCA(L0,NY,NX)
      ZMG(L0,NY,NX)=FY*ZMG(L0,NY,NX)
      ZNA(L0,NY,NX)=FY*ZNA(L0,NY,NX)
      ZKA(L0,NY,NX)=FY*ZKA(L0,NY,NX)
      ZOH(L0,NY,NX)=FY*ZOH(L0,NY,NX)
      ZSO4(L0,NY,NX)=FY*ZSO4(L0,NY,NX)
      ZCL(L0,NY,NX)=FY*ZCL(L0,NY,NX)
      ZCO3(L0,NY,NX)=FY*ZCO3(L0,NY,NX)
      ZHCO3(L0,NY,NX)=FY*ZHCO3(L0,NY,NX)
      ZALOH1(L0,NY,NX)=FY*ZALOH1(L0,NY,NX)
      ZALOH2(L0,NY,NX)=FY*ZALOH2(L0,NY,NX)
      ZALOH3(L0,NY,NX)=FY*ZALOH3(L0,NY,NX)
      ZALOH4(L0,NY,NX)=FY*ZALOH4(L0,NY,NX)
      ZALS(L0,NY,NX)=FY*ZALS(L0,NY,NX)
      ZFEOH1(L0,NY,NX)=FY*ZFEOH1(L0,NY,NX)
      ZFEOH2(L0,NY,NX)=FY*ZFEOH2(L0,NY,NX)
      ZFEOH3(L0,NY,NX)=FY*ZFEOH3(L0,NY,NX)
      ZFEOH4(L0,NY,NX)=FY*ZFEOH4(L0,NY,NX)
      ZFES(L0,NY,NX)=FY*ZFES(L0,NY,NX)
      ZCAO(L0,NY,NX)=FY*ZCAO(L0,NY,NX)
      ZCAC(L0,NY,NX)=FY*ZCAC(L0,NY,NX)
      ZCAH(L0,NY,NX)=FY*ZCAH(L0,NY,NX)
      ZCAS(L0,NY,NX)=FY*ZCAS(L0,NY,NX)
      ZMGO(L0,NY,NX)=FY*ZMGO(L0,NY,NX)
      ZMGC(L0,NY,NX)=FY*ZMGC(L0,NY,NX)
      ZMGH(L0,NY,NX)=FY*ZMGH(L0,NY,NX)
      ZMGS(L0,NY,NX)=FY*ZMGS(L0,NY,NX)
      ZNAC(L0,NY,NX)=FY*ZNAC(L0,NY,NX)
      ZNAS(L0,NY,NX)=FY*ZNAS(L0,NY,NX)
      ZKAS(L0,NY,NX)=FY*ZKAS(L0,NY,NX)
      H0PO4(L0,NY,NX)=FY*H0PO4(L0,NY,NX)
      H3PO4(L0,NY,NX)=FY*H3PO4(L0,NY,NX)
      ZFE1P(L0,NY,NX)=FY*ZFE1P(L0,NY,NX)
      ZFE2P(L0,NY,NX)=FY*ZFE2P(L0,NY,NX)
      ZCA0P(L0,NY,NX)=FY*ZCA0P(L0,NY,NX)
      ZCA1P(L0,NY,NX)=FY*ZCA1P(L0,NY,NX)
      ZCA2P(L0,NY,NX)=FY*ZCA2P(L0,NY,NX)
      ZMG1P(L0,NY,NX)=FY*ZMG1P(L0,NY,NX)
      ENDIF
C
C     POND SALTS IN BAND
C
      IF(L0.NE.0)THEN
      H1POB(L0,NY,NX)=FY*H1POB(L0,NY,NX)
      H2POB(L0,NY,NX)=FY*H2POB(L0,NY,NX)
      IF(ISALTG.NE.0)THEN
      H0POB(L0,NY,NX)=FY*H0POB(L0,NY,NX)
      H3POB(L0,NY,NX)=FY*H3POB(L0,NY,NX)
      ZFE1PB(L0,NY,NX)=FY*ZFE1PB(L0,NY,NX)
      ZFE2PB(L0,NY,NX)=FY*ZFE2PB(L0,NY,NX)
      ZCA0PB(L0,NY,NX)=FY*ZCA0PB(L0,NY,NX)
      ZCA1PB(L0,NY,NX)=FY*ZCA1PB(L0,NY,NX)
      ZCA2PB(L0,NY,NX)=FY*ZCA2PB(L0,NY,NX)
      ZMG1PB(L0,NY,NX)=FY*ZMG1PB(L0,NY,NX)
      ENDIF
C
C     POND GASEOUS GASES
C
      CO2G(L0,NY,NX)=FY*CO2G(L0,NY,NX)
      CH4G(L0,NY,NX)=FY*CH4G(L0,NY,NX)
      OXYG(L0,NY,NX)=FY*OXYG(L0,NY,NX)
      Z2GG(L0,NY,NX)=FY*Z2GG(L0,NY,NX)
      Z2OG(L0,NY,NX)=FY*Z2OG(L0,NY,NX)
      ZNH3G(L0,NY,NX)=FY*ZNH3G(L0,NY,NX)
      H2GG(L0,NY,NX)=FY*H2GG(L0,NY,NX)
      ENDIF
C
C     POND AQUEOUS GASES
C
      CO2S(L0,NY,NX)=FY*CO2S(L0,NY,NX)
      CH4S(L0,NY,NX)=FY*CH4S(L0,NY,NX)
      OXYS(L0,NY,NX)=FY*OXYS(L0,NY,NX)
      Z2GS(L0,NY,NX)=FY*Z2GS(L0,NY,NX)
      Z2OS(L0,NY,NX)=FY*Z2OS(L0,NY,NX)
      H2GS(L0,NY,NX)=FY*H2GS(L0,NY,NX)
C
C     POND ORGANIC MATTER
C
      IF(IFLGL(L,3).EQ.0)THEN
      DO 7880 K=0,4
      OQC(K,L0,NY,NX)=FY*OQC(K,L0,NY,NX)
      OQN(K,L0,NY,NX)=FY*OQN(K,L0,NY,NX)
      OQP(K,L0,NY,NX)=FY*OQP(K,L0,NY,NX)
      OQA(K,L0,NY,NX)=FY*OQA(K,L0,NY,NX)
      OQCH(K,L0,NY,NX)=FY*OQCH(K,L0,NY,NX)
      OQNH(K,L0,NY,NX)=FY*OQNH(K,L0,NY,NX)
      OQPH(K,L0,NY,NX)=FY*OQPH(K,L0,NY,NX)
      OQAH(K,L0,NY,NX)=FY*OQAH(K,L0,NY,NX)
7880  CONTINUE
      ENDIF
C
C     POND PLANT MATTER
C
      IF(L0.NE.0)THEN
      DO 8910 NZ=1,NP(NY,NX)
      IF(WTRTL(1,L0,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WTRTL(1,L1,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      DO 8905 N=1,MY(NZ,NY,NX)
C
C     POND ROOT GASES
C
      CO2A(N,L0,NZ,NY,NX)=FY*CO2A(N,L0,NZ,NY,NX)
      OXYA(N,L0,NZ,NY,NX)=FY*OXYA(N,L0,NZ,NY,NX)
      CH4A(N,L0,NZ,NY,NX)=FY*CH4A(N,L0,NZ,NY,NX)
      Z2OA(N,L0,NZ,NY,NX)=FY*Z2OA(N,L0,NZ,NY,NX)
      ZH3A(N,L0,NZ,NY,NX)=FY*ZH3A(N,L0,NZ,NY,NX)
      H2GA(N,L0,NZ,NY,NX)=FY*H2GA(N,L0,NZ,NY,NX)
      CO2P(N,L0,NZ,NY,NX)=FY*CO2P(N,L0,NZ,NY,NX)
      OXYP(N,L0,NZ,NY,NX)=FY*OXYP(N,L0,NZ,NY,NX)
      CH4P(N,L0,NZ,NY,NX)=FY*CH4P(N,L0,NZ,NY,NX)
      Z2OP(N,L0,NZ,NY,NX)=FY*Z2OP(N,L0,NZ,NY,NX)
      ZH3P(N,L0,NZ,NY,NX)=FY*ZH3P(N,L0,NZ,NY,NX)
      H2GP(N,L0,NZ,NY,NX)=FY*H2GP(N,L0,NZ,NY,NX)
C
C     POND ROOTS
C
      DO 8970 NR=1,NRT(NZ,NY,NX)
      WTRT1(N,L0,NR,NZ,NY,NX)=FY*WTRT1(N,L0,NR,NZ,NY,NX)
      WTRT1N(N,L0,NR,NZ,NY,NX)=FY*WTRT1N(N,L0,NR,NZ,NY,NX)
      WTRT1P(N,L0,NR,NZ,NY,NX)=FY*WTRT1P(N,L0,NR,NZ,NY,NX)
      WTRT2(N,L0,NR,NZ,NY,NX)=FY*WTRT2(N,L0,NR,NZ,NY,NX)
      WTRT2N(N,L0,NR,NZ,NY,NX)=FY*WTRT2N(N,L0,NR,NZ,NY,NX)
      WTRT2P(N,L0,NR,NZ,NY,NX)=FY*WTRT2P(N,L0,NR,NZ,NY,NX)
      RTLG1(N,L0,NR,NZ,NY,NX)=FY*RTLG1(N,L0,NR,NZ,NY,NX)
      RTLG2(N,L0,NR,NZ,NY,NX)=FY*RTLG2(N,L0,NR,NZ,NY,NX)
      RTN2(N,L0,NR,NZ,NY,NX)=FY*RTN2(N,L0,NR,NZ,NY,NX)
8970  CONTINUE
      CPOOLR(N,L0,NZ,NY,NX)=FY*CPOOLR(N,L0,NZ,NY,NX)
      ZPOOLR(N,L0,NZ,NY,NX)=FY*ZPOOLR(N,L0,NZ,NY,NX)
      PPOOLR(N,L0,NZ,NY,NX)=FY*PPOOLR(N,L0,NZ,NY,NX)
      WTRTL(N,L0,NZ,NY,NX)=FY*WTRTL(N,L0,NZ,NY,NX)
      WTRTD(N,L0,NZ,NY,NX)=FY*WTRTD(N,L0,NZ,NY,NX)
      WSRTL(N,L0,NZ,NY,NX)=FY*WSRTL(N,L0,NZ,NY,NX)
      RTN1(N,L0,NZ,NY,NX)=FY*RTN1(N,L0,NZ,NY,NX)
      RTNL(N,L0,NZ,NY,NX)=FY*RTNL(N,L0,NZ,NY,NX)
      RTLGP(N,L0,NZ,NY,NX)=FY*RTLGP(N,L0,NZ,NY,NX)
      RTDNP(N,L0,NZ,NY,NX)=FY*RTDNP(N,L0,NZ,NY,NX)
      RTVLP(N,L0,NZ,NY,NX)=FY*RTVLP(N,L0,NZ,NY,NX)
      RTVLW(N,L0,NZ,NY,NX)=FY*RTVLW(N,L0,NZ,NY,NX)
      RRAD1(N,L0,NZ,NY,NX)=FY*RRAD1(N,L0,NZ,NY,NX)
      RRAD2(N,L0,NZ,NY,NX)=FY*RRAD2(N,L0,NZ,NY,NX)
      RTARP(N,L0,NZ,NY,NX)=FY*RTARP(N,L0,NZ,NY,NX)
      RTLGA(N,L0,NZ,NY,NX)=FY*RTLGA(N,L0,NZ,NY,NX)
8905  CONTINUE
C
C     POND ROOT NODULES
C
      WTNDL(L0,NZ,NY,NX)=FY*WTNDL(L0,NZ,NY,NX)
      WTNDLN(L0,NZ,NY,NX)=FY*WTNDLN(L0,NZ,NY,NX)
      WTNDLP(L0,NZ,NY,NX)=FY*WTNDLP(L0,NZ,NY,NX)
      CPOOLN(L0,NZ,NY,NX)=FY*CPOOLN(L0,NZ,NY,NX)
      ZPOOLN(L0,NZ,NY,NX)=FY*ZPOOLN(L0,NZ,NY,NX)
      PPOOLN(L0,NZ,NY,NX)=FY*PPOOLN(L0,NZ,NY,NX)
      ENDIF
8910  CONTINUE
      ENDIF
C
C     POND SEDIMENT
C
      IF(FX.EQ.1.0)THEN
      IF(L0.NE.0)THEN
      SAND(L0,NY,NX)=FY*SAND(L0,NY,NX)
      SILT(L0,NY,NX)=FY*SILT(L0,NY,NX)
      CLAY(L0,NY,NX)=FY*CLAY(L0,NY,NX)
      XCEC(L0,NY,NX)=FY*XCEC(L0,NY,NX)
      XAEC(L0,NY,NX)=FY*XAEC(L0,NY,NX)
      ENDIF
C
C     POND FERTILIZER IN BAND, NON-BAND
C
      ZNH4FA(L0,NY,NX)=FY*ZNH4FA(L0,NY,NX)
      ZNH3FA(L0,NY,NX)=FY*ZNH3FA(L0,NY,NX)
      ZNHUFA(L0,NY,NX)=FY*ZNHUFA(L0,NY,NX)
      ZNO3FA(L0,NY,NX)=FY*ZNO3FA(L0,NY,NX)
      ZNH4FB(L0,NY,NX)=FY*ZNH4FB(L0,NY,NX)
      ZNH3FB(L0,NY,NX)=FY*ZNH3FB(L0,NY,NX)
      ZNHUFB(L0,NY,NX)=FY*ZNHUFB(L0,NY,NX)
      ZNO3FB(L0,NY,NX)=FY*ZNO3FB(L0,NY,NX)
C
C     POND ADSORBED CATIONS IN BAND, NON-BAND
C
      XN4(L0,NY,NX)=FY*XN4(L0,NY,NX)
      XNB(L0,NY,NX)=FY*XNB(L0,NY,NX)
      XHY(L0,NY,NX)=FY*XHY(L0,NY,NX)
      XAL(L0,NY,NX)=FY*XAL(L0,NY,NX)
      XFE(L0,NY,NX)=FY*XFE(L0,NY,NX)
      XCA(L0,NY,NX)=FY*XCA(L0,NY,NX)
      XMG(L0,NY,NX)=FY*XMG(L0,NY,NX)
      XNA(L0,NY,NX)=FY*XNA(L0,NY,NX)
      XKA(L0,NY,NX)=FY*XKA(L0,NY,NX)
      XHC(L0,NY,NX)=FY*XHC(L0,NY,NX)
      XALO2(L0,NY,NX)=FY*XALO2(L0,NY,NX)
      XFEO2(L0,NY,NX)=FY*XFEO2(L0,NY,NX)
C
C     POND ADSORBED ANIONS IN BAND, NON-BAND
C
      XOH0(L0,NY,NX)=FY*XOH0(L0,NY,NX)
      XOH1(L0,NY,NX)=FY*XOH1(L0,NY,NX)
      XOH2(L0,NY,NX)=FY*XOH2(L0,NY,NX)
      XH1P(L0,NY,NX)=FY*XH1P(L0,NY,NX)
      XH2P(L0,NY,NX)=FY*XH2P(L0,NY,NX)
      XOH0B(L0,NY,NX)=FY*XOH0B(L0,NY,NX)
      XOH1B(L0,NY,NX)=FY*XOH1B(L0,NY,NX)
      XOH2B(L0,NY,NX)=FY*XOH2B(L0,NY,NX)
      XH1PB(L0,NY,NX)=FY*XH1PB(L0,NY,NX)
      XH2PB(L0,NY,NX)=FY*XH2PB(L0,NY,NX)
C
C     POND PRECIPITATES IN BAND, NON-BAND
C
      PALOH(L0,NY,NX)=FY*PALOH(L0,NY,NX)
      PFEOH(L0,NY,NX)=FY*PFEOH(L0,NY,NX)
      PCACO(L0,NY,NX)=FY*PCACO(L0,NY,NX)
      PCASO(L0,NY,NX)=FY*PCASO(L0,NY,NX)
      PALPO(L0,NY,NX)=FY*PALPO(L0,NY,NX)
      PFEPO(L0,NY,NX)=FY*PFEPO(L0,NY,NX)
      PCAPD(L0,NY,NX)=FY*PCAPD(L0,NY,NX)
      PCAPH(L0,NY,NX)=FY*PCAPH(L0,NY,NX)
      PCAPM(L0,NY,NX)=FY*PCAPM(L0,NY,NX)
      PALPB(L0,NY,NX)=FY*PALPB(L0,NY,NX)
      PFEPB(L0,NY,NX)=FY*PFEPB(L0,NY,NX)
      PCPDB(L0,NY,NX)=FY*PCPDB(L0,NY,NX)
      PCPHB(L0,NY,NX)=FY*PCPHB(L0,NY,NX)
      PCPMB(L0,NY,NX)=FY*PCPMB(L0,NY,NX)
C
C     POND ORGANIC MATTER
C
      IF(IFLGL(L,3).EQ.0)THEN
      ORGC(L0,NY,NX)=FY*ORGC(L0,NY,NX)
      DO 7865 K=0,5
      DO 7865 N=1,7
      DO 7865 M=1,3
      OMC(M,N,K,L0,NY,NX)=FY*OMC(M,N,K,L0,NY,NX)
      OMN(M,N,K,L0,NY,NX)=FY*OMN(M,N,K,L0,NY,NX)
      OMP(M,N,K,L0,NY,NX)=FY*OMP(M,N,K,L0,NY,NX)
7865  CONTINUE
      DO 7885 K=0,4
      DO 7875 M=1,2
      ORC(M,K,L0,NY,NX)=FY*ORC(M,K,L0,NY,NX)
      ORN(M,K,L0,NY,NX)=FY*ORN(M,K,L0,NY,NX)
      ORP(M,K,L0,NY,NX)=FY*ORP(M,K,L0,NY,NX)
7875  CONTINUE
      OHC(K,L0,NY,NX)=FY*OHC(K,L0,NY,NX)
      OHN(K,L0,NY,NX)=FY*OHN(K,L0,NY,NX)
      OHP(K,L0,NY,NX)=FY*OHP(K,L0,NY,NX)
      OHA(K,L0,NY,NX)=FY*OHA(K,L0,NY,NX)
      DO 7870 M=1,5
      OSC(M,K,L0,NY,NX)=FY*OSC(M,K,L0,NY,NX)
      OSA(M,K,L0,NY,NX)=FY*OSA(M,K,L0,NY,NX)
      OSN(M,K,L0,NY,NX)=FY*OSN(M,K,L0,NY,NX)
      OSP(M,K,L0,NY,NX)=FY*OSP(M,K,L0,NY,NX)
7870  CONTINUE
7885  CONTINUE
      ENDIF
      ENDIF
C
C     CLOSE OUT EMPTY LAYERS
C
      IF(NN.EQ.1)THEN
      IF(BKDS(L0,NY,NX).LE.ZERO.AND.BKDS(L1,NY,NX).LE.ZERO
     3.AND.VOLW(L0,NY,NX)+VOLI(L0,NY,NX).LE.ZEROS(NY,NX))THEN
      CDPTH(L1,NY,NX)=CDPTH(L0,NY,NX)
      CDPTHY(L1,NY,NX)=CDPTHY(L0,NY,NX)
      ENDIF
      ENDIF
C     IF(I.EQ.178)THEN
C     WRITE(*,5599)'POND2',I,J,NFZ,NX,NY,L,L0,L1,NU(NY,NX),NN,FX,FY
C    2,DDLYRX(NN),VOLY(L0,NY,NX),VOLX(L0,NY,NX),VOLW(L0,NY,NX)
C    3,VOLI(L0,NY,NX),VOLY(L1,NY,NX),VOLX(L1,NY,NX),VOLW(L1,NY,NX)
C    4,VOLI(L1,NY,NX),CDPTH(L0,NY,NX),CDPTH(L1,NY,NX)
C    5,VOLW(L0,NY,NX)+VOLW(L1,NY,NX)
C    5,(OSC(M,4,L0,NY,NX),M=1,5),(OSC(M,4,L1,NY,NX),M=1,5)
C    5,(OQC(K,L0,NY,NX),K=0,4),(OQC(K,L1,NY,NX),K=0,4)
C    6,DLYR(3,L0,NY,NX),DLYR(3,L1,NY,NX)
C    5,TKS(L0,NY,NX),TKS(L1,NY,NX)
C    5,VHCP(L0,NY,NX),VHCP(L1,NY,NX)
C    6,ZNH4S(L0,NY,NX),ZNH4B(L0,NY,NX),ZNH3S(L0,NY,NX),ZNH3B(L0,NY,NX)
C    6,ZNH4S(L1,NY,NX),ZNH4B(L1,NY,NX),ZNH3S(L1,NY,NX),ZNH3B(L1,NY,NX)
C    6,(WTRT1(1,L1,NR,1,NY,NX),NR=1,NRT(1,NY,NX))
C    6,(WTRT2(1,L1,NR,1,NY,NX),NR=1,NRT(1,NY,NX))
C    6,WTRTL(1,L1,1,NY,NX),WTNDL(L1,1,NY,NX)
C    6,CPOOLR(1,L1,1,NY,NX),ZPOOLR(1,L1,1,NY,NX)
C    6,OXYA(1,L1,1,NY,NX),OXYP(1,L1,1,NY,NX)
C    6,(WTRT1(1,L0,NR,1,NY,NX),NR=1,NRT(1,NY,NX))
C    6,(WTRT2(1,L0,NR,1,NY,NX),NR=1,NRT(1,NY,NX))
C    6,WTRTL(1,L0,1,NY,NX),WTNDL(L0,1,NY,NX)
C    6,CPOOLR(1,L0,1,NY,NX),ZPOOLR(1,L0,1,NY,NX)
C    6,OXYA(1,L0,1,NY,NX),OXYP(1,L0,1,NY,NX)
C     ENDIF
C
C     REDISTRIBUTE SOIL MATERIAL
C
      ELSE
C     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,5591)'SOIL1',I,J,NFZ,NX,NY,L,L0,L1,NN,DDLYRX(NN),FX,FO
C    1,ZNH4SH(L0,NY,NX),ZNH4SH(L1,NY,NX)
C    1,BKDS(L0,NY,NX),BKDS(L1,NY,NX),FXBKDS,FXVOLW
C    5,ORGC(L0,NY,NX),ORGC(L1,NY,NX),CORGCI(L0,NY,NX),CORGCI(L1,NY,NX)
C    5,VHCM(L0,NY,NX),VHCM(L1,NY,NX),VHCP(L0,NY,NX),VHCP(L1,NY,NX)
C    5,CLAY(L0,NY,NX),CLAY(L1,NY,NX),SILT(L0,NY,NX),SILT(L1,NY,NX)
C    5,SAND(L0,NY,NX),SAND(L1,NY,NX),ORGC(L0,NY,NX),ORGC(L1,NY,NX)
C    2,VOLA(L0,NY,NX),VOLW(L0,NY,NX),VOLI(L0,NY,NX),VOLP(L0,NY,NX)
C    3,VOLA(L1,NY,NX),VOLW(L1,NY,NX),VOLI(L1,NY,NX),VOLP(L1,NY,NX)
C    3,VOLX(L0,NY,NX),VOLX(L1,NY,NX),VOLY(L0,NY,NX)+VOLY(L1,NY,NX)
C    6,DLYR(3,L0,NY,NX),DLYR(3,L1,NY,NX)
C    4,VOLA(L0,NY,NX)+VOLA(L1,NY,NX),VOLW(L0,NY,NX)+VOLW(L1,NY,NX)
C    4,VOLI(L0,NY,NX)+VOLI(L1,NY,NX),VOLP(L0,NY,NX)+VOLP(L1,NY,NX)
C    4,CDPTH(L0,NY,NX),CDPTH(L1,NY,NX)
C    5,PCAPD(L0,NY,NX),PCAPM(L1,NY,NX),FXPCAPD
C    5,((OMC(M,N,1,L0,NY,NX),M=1,3),N=1,7)
C    6,((OMC(M,N,1,L1,NY,NX),M=1,3),N=1,7)
C    5,(OSC(M,4,L0,NY,NX),M=1,2),(OSC(M,4,L1,NY,NX),M=1,2)
C    5,(ORC(M,4,L0,NY,NX),M=1,2),(ORC(M,4,L1,NY,NX),M=1,2)
C    6,((OSC(M,4,L0,NY,NX)+OSC(M,4,L1,NY,NX)),M=1,2)
C    5,(OQC(K,L0,NY,NX),K=0,4),(OQC(K,L1,NY,NX),K=0,4)
C    6,DLYR(3,L0,NY,NX),DLYR(3,L1,NY,NX)
C    5,TKS(L0,NY,NX),TKS(L1,NY,NX),VHCP(L0,NY,NX),VHCP(L1,NY,NX)
C    6,TKS(L0,NY,NX)*VHCP(L0,NY,NX)+TKS(L1,NY,NX)*VHCP(L1,NY,NX)
C    7,VHCM(L0,NY,NX),VHCM(L1,NY,NX)
C    6,ZNH4S(L0,NY,NX),ZNH4B(L0,NY,NX),ZNH3S(L0,NY,NX),ZNH3B(L0,NY,NX)
C    6,ZNH4S(L1,NY,NX),ZNH4B(L1,NY,NX),ZNH3S(L1,NY,NX),ZNH3B(L1,NY,NX)
C    6,OQAH(K,L,NY,NX),OQAH(K,L0,NY,NX),OQAH(K,L1,NY,NX),FXOQAH
C    7,CH4G(L0,NY,NX),CH4G(L1,NY,NX),FXCH4G
C    6,(WTRT1(1,L1,NR,2,NY,NX),NR=1,NRT(1,NY,NX))
C    6,(WTRT2(1,L1,NR,2,NY,NX),NR=1,NRT(1,NY,NX))
C    6,WTRTL(1,L1,2,NY,NX),WTNDL(L1,2,NY,NX)
C    6,CPOOLR(1,L1,2,NY,NX),ZPOOLR(1,L1,2,NY,NX)
C    6,OXYA(1,L1,2,NY,NX),OXYP(1,L1,2,NY,NX)
C    6,(WTRT1(1,L0,NR,2,NY,NX),NR=1,NRT(1,NY,NX))
C    6,(WTRT2(1,L0,NR,2,NY,NX),NR=1,NRT(1,NY,NX))
C    6,WTRTL(1,L0,2,NY,NX),WTNDL(L0,2,NY,NX)
C    6,CPOOLR(1,L0,2,NY,NX),ZPOOLR(1,L0,2,NY,NX)
C    6,OXYG(L0,NY,NX),OXYS(L0,NY,NX)
C    6,OXYG(L1,NY,NX),OXYS(L1,NY,NX)
5591  FORMAT(A8,9I4,60E18.10)
C     ENDIF
C
C     REDISTRIBUTE SOIL CONTENTS FROM SOURCE TO DESTINATION LAYER
C
      IF(L0.NE.0)THEN
      IF(DLYR(3,L1,NY,NX).GT.ZERO
     2.AND.DLYR(3,L0,NY,NX).GT.ZERO)THEN
C
C     RESET SOIL FERTILIZER NH4, NO3, PO4 BAND DIMNSIONS
C
      IF(IFNHB(NY,NX).EQ.1.AND.ROWN(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CDPTH(L-1,NY,NX).LT.DPNH4(NY,NX))THEN
      WDNHBDL=WDNHB(L,NY,NX)*DLYR(3,L,NY,NX)
      WDNHBD0=WDNHB(L0,NY,NX)*DLYR(3,L0,NY,NX)
      WDNHBD1=WDNHB(L1,NY,NX)*DLYR(3,L1,NY,NX)
      FXWDNHB=AMIN1(FX*WDNHBDL,WDNHBD0)
      WDNHBD1=WDNHBD1+FXWDNHB
      WDNHBD0=WDNHBD0-FXWDNHB
      WDNHB(L1,NY,NX)=WDNHBD1/DLYR(3,L1,NY,NX)
      WDNHB(L0,NY,NX)=WDNHBD0/DLYR(3,L0,NY,NX)
      IF(CDPTH(L,NY,NX).GE.DPNH4(NY,NX))THEN
      FXDPNHB=AMIN1(FX*DPNHB(L,NY,NX),DPNHB(L0,NY,NX))
      DPNHB(L1,NY,NX)=DPNHB(L1,NY,NX)+FXDPNHB
      DPNHB(L0,NY,NX)=DPNHB(L0,NY,NX)-FXDPNHB
      ENDIF
      VLNHB(L1,NY,NX)=AMAX1(0.0,AMIN1(0.999,WDNHB(L1,NY,NX)
     2/ROWN(NY,NX)*DPNHB(L1,NY,NX)/DLYR(3,L1,NY,NX)))
      VLNHB(L0,NY,NX)=AMAX1(0.0,AMIN1(0.999,WDNHB(L0,NY,NX)
     2/ROWN(NY,NX)*DPNHB(L0,NY,NX)/DLYR(3,L0,NY,NX)))
      VLNH4(L1,NY,NX)=1.0-VLNHB(L1,NY,NX)
      VLNH4(L0,NY,NX)=1.0-VLNHB(L0,NY,NX)
C     WRITE(*,6601)'VLNHB',I,J,NX,NY,L,L0,L1,FX
C    2,WDNHB(L0,NY,NX),WDNHB(L1,NY,NX)
C    3,DPNHB(L0,NY,NX),DPNHB(L1,NY,NX)
C    2,VLNHB(L0,NY,NX),DLYR(3,L0,NY,NX)
C    3,VLNHB(L1,NY,NX),DLYR(3,L1,NY,NX)
C    2,VLNHB(L0,NY,NX)*DLYR(3,L0,NY,NX)
C    3+VLNHB(L1,NY,NX)*DLYR(3,L1,NY,NX)
C    3,WDNHBDL,WDNHBD0,WDNHBD1,FXWDNHB,FXDPNHB,DPNH4(NY,NX)
6601  FORMAT(A8,7I4,30E16.8)
      ENDIF
      ENDIF
      IF(IFNOB(NY,NX).EQ.1.AND.ROWO(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CDPTH(L-1,NY,NX).LT.DPNO3(NY,NX))THEN
      WDNOBDL=WDNOB(L,NY,NX)*DLYR(3,L,NY,NX)
      WDNOBD0=WDNOB(L0,NY,NX)*DLYR(3,L0,NY,NX)
      WDNOBD1=WDNOB(L1,NY,NX)*DLYR(3,L1,NY,NX)
      FXWDNOB=AMIN1(FX*WDNOBDL,WDNOBD0)
      WDNOBD1=WDNOBD1+FXWDNOB
      WDNOBD0=WDNOBD0-FXWDNOB
      WDNOB(L1,NY,NX)=WDNOBD1/DLYR(3,L1,NY,NX)
      WDNOB(L0,NY,NX)=WDNOBD0/DLYR(3,L0,NY,NX)
      IF(CDPTH(L,NY,NX).GE.DPNO3(NY,NX))THEN
      FXDPNOB=AMIN1(FX*DPNOB(L,NY,NX),DPNOB(L0,NY,NX))
      DPNOB(L1,NY,NX)=DPNOB(L1,NY,NX)+FXDPNOB
      DPNOB(L0,NY,NX)=DPNOB(L0,NY,NX)-FXDPNOB
      ENDIF
      VLNOB(L1,NY,NX)=AMAX1(0.0,AMIN1(0.999,WDNOB(L1,NY,NX)
     2/ROWO(NY,NX)*DPNOB(L1,NY,NX)/DLYR(3,L1,NY,NX)))
      VLNOB(L0,NY,NX)=AMAX1(0.0,AMIN1(0.999,WDNOB(L0,NY,NX)
     2/ROWO(NY,NX)*DPNOB(L0,NY,NX)/DLYR(3,L0,NY,NX)))
      VLNO3(L1,NY,NX)=1.0-VLNOB(L1,NY,NX)
      VLNO3(L0,NY,NX)=1.0-VLNOB(L0,NY,NX)
C     WRITE(*,6601)'VLNOB',I,J,NX,NY,L,L0,L1,FX
C    2,WDNOB(L0,NY,NX),WDNOB(L1,NY,NX)
C    2,DPNOB(L0,NY,NX),DPNOB(L1,NY,NX)
C    2,VLNOB(L0,NY,NX),DLYR(3,L0,NY,NX)
C    3,VLNOB(L1,NY,NX),DLYR(3,L1,NY,NX)
C    2,VLNOB(L0,NY,NX)*DLYR(3,L0,NY,NX)
C    3+VLNOB(L1,NY,NX)*DLYR(3,L1,NY,NX)
C    3,WDNOBDL,WDNOBD0,WDNOBD1,FXWDNOB,FXDPNOB,DPNO3(NY,NX)
      ENDIF
      ENDIF
      IF(IFPOB(NY,NX).EQ.1.AND.ROWP(NY,NX).GT.0.0)THEN
      IF(L.EQ.NU(NY,NX).OR.CDPTH(L-1,NY,NX).LT.DPPO4(NY,NX))THEN
      WDPOBDL=WDPOB(L,NY,NX)*DLYR(3,L,NY,NX)
      WDPOBD0=WDPOB(L0,NY,NX)*DLYR(3,L0,NY,NX)
      WDPOBD1=WDPOB(L1,NY,NX)*DLYR(3,L1,NY,NX)
      FXWDPOB=AMIN1(FX*WDPOBDL,WDPOBD0)
      WDPOBD1=WDPOBD1+FXWDPOB
      WDPOBD0=WDPOBD0-FXWDPOB
      WDPOB(L1,NY,NX)=WDPOBD1/DLYR(3,L1,NY,NX)
      WDPOB(L0,NY,NX)=WDPOBD0/DLYR(3,L0,NY,NX)
      IF(CDPTH(L,NY,NX).GE.DPPO4(NY,NX))THEN
      FXDPPOB=AMIN1(FX*DPPOB(L,NY,NX),DPPOB(L0,NY,NX))
      DPPOB(L1,NY,NX)=DPPOB(L1,NY,NX)+FXDPPOB
      DPPOB(L0,NY,NX)=DPPOB(L0,NY,NX)-FXDPPOB
      ENDIF
      VLPOB(L1,NY,NX)=AMAX1(0.0,AMIN1(0.999,WDPOB(L1,NY,NX)
     2/ROWP(NY,NX)*DPPOB(L1,NY,NX)/DLYR(3,L1,NY,NX)))
      VLPOB(L0,NY,NX)=AMAX1(0.0,AMIN1(0.999,WDPOB(L0,NY,NX)
     2/ROWP(NY,NX)*DPPOB(L0,NY,NX)/DLYR(3,L0,NY,NX)))
      VLPO4(L1,NY,NX)=1.0-VLPOB(L1,NY,NX)
      VLPO4(L0,NY,NX)=1.0-VLPOB(L0,NY,NX)
C     WRITE(*,6601)'VLPOB',I,J,NX,NY,L,L0,L1,FX
C    2,WDPOB(L0,NY,NX),WDPOB(L1,NY,NX)
C    3,DPPOB(L0,NY,NX),DPPOB(L1,NY,NX)
C    2,VLPOB(L0,NY,NX),DLYR(3,L0,NY,NX)
C    3,VLPOB(L1,NY,NX),DLYR(3,L1,NY,NX)
C    2,VLPOB(L0,NY,NX)*DLYR(3,L0,NY,NX)
C    3+VLPOB(L1,NY,NX)*DLYR(3,L1,NY,NX)
C    3,WDPOBDL,WDPOBD0,WDPOBD1,FXWDPOB,FXDPPOB,DPPO4(NY,NX)
      ENDIF
      ENDIF
      ENDIF
C
C     SOIL SEDIMENT
C
C     IF(L0.EQ.L.OR.BKDSI(L0,NY,NX).LE.ZERO)THEN
      FBO=FX
C     ELSE
C     FBO=AMIN1(1.0,FX*BKDSI(L1,NY,NX)/BKDSI(L0,NY,NX))
C     ENDIF
C     BKDS(L1,NY,NX)=(1.0-FX)*BKDS(L1,NY,NX)+FX*BKDSI(L0,NY,NX)
      PH(L1,NY,NX)=(1.0-FX)*PH(L1,NY,NX)+FX*PH(L0,NY,NX)
      FXSAND=FBO*SAND(L0,NY,NX)
      SAND(L1,NY,NX)=SAND(L1,NY,NX)+FXSAND
      SAND(L0,NY,NX)=SAND(L0,NY,NX)-FXSAND
      FXSILT=FBO*SILT(L0,NY,NX)
      SILT(L1,NY,NX)=SILT(L1,NY,NX)+FXSILT
      SILT(L0,NY,NX)=SILT(L0,NY,NX)-FXSILT
      FXCLAY=FBO*CLAY(L0,NY,NX)
      CLAY(L1,NY,NX)=CLAY(L1,NY,NX)+FXCLAY
      CLAY(L0,NY,NX)=CLAY(L0,NY,NX)-FXCLAY
      FXROCK=FBO*ROCK(L0,NY,NX)
      ROCK(L1,NY,NX)=ROCK(L1,NY,NX)+FXROCK
      ROCK(L0,NY,NX)=ROCK(L0,NY,NX)-FXROCK
C
C     SOIL WATER, VAPOR, ICE AND HEAT
C
      IF(FHOL(L1,NY,NX).GT.ZERO.AND.FHOL(L0,NY,NX).GT.ZERO)THEN
C     IF(L0.EQ.L.OR.FHOLI(L0,NY,NX).LE.ZERO)THEN
      FHO=FX
C     ELSE
C     FHO=AMIN1(1.0,FX*FHOLI(L1,NY,NX)/FHOLI(L0,NY,NX))
C     ENDIF
      FHOL(L1,NY,NX)=(1.0-FX)*FHOL(L1,NY,NX)+FX*FHOL(L0,NY,NX)
      FXVOLWH=FHO*VOLWH(L0,NY,NX)
      VOLWH(L1,NY,NX)=VOLWH(L1,NY,NX)+FXVOLWH
      VOLWH(L0,NY,NX)=VOLWH(L0,NY,NX)-FXVOLWH
      FXVOLIH=FHO*VOLIH(L0,NY,NX)
      VOLIH(L1,NY,NX)=VOLIH(L1,NY,NX)+FXVOLIH
      VOLIH(L0,NY,NX)=VOLIH(L0,NY,NX)-FXVOLIH
      FXVOLAH=FHO*VOLAH(L0,NY,NX)
      VOLAH(L1,NY,NX)=VOLAH(L1,NY,NX)+FXVOLAH
      VOLAH(L0,NY,NX)=VOLAH(L0,NY,NX)-FXVOLAH
      ENDIF
      ENDIF
C     IF(L0.EQ.L.OR.POROSI(L0,NY,NX).LE.ZERO)THEN
      FWO=FX
C     ELSE
C     FWO=AMIN1(1.0,FX*POROSI(L1,NY,NX)/POROSI(L0,NY,NX))
C     ENDIF
      IF(L0.EQ.0)THEN
      FXVOLW=FX*AMAX1(0.0,XVOLWP-VOLWD(NY,NX))
      FXVOLV=0.0
      ELSE
      FXVOLW=FWO*VOLW(L0,NY,NX)
      FXVOLV=FWO*VOLV(L0,NY,NX)
      ENDIF
      VOLW(L1,NY,NX)=VOLW(L1,NY,NX)+FXVOLW
      VOLW(L0,NY,NX)=VOLW(L0,NY,NX)-FXVOLW
      VOLV(L1,NY,NX)=VOLV(L1,NY,NX)+FXVOLV
      VOLV(L0,NY,NX)=VOLV(L0,NY,NX)-FXVOLV
      FXVOLI=FWO*VOLI(L0,NY,NX)
      VOLI(L1,NY,NX)=VOLI(L1,NY,NX)+FXVOLI
      VOLI(L0,NY,NX)=VOLI(L0,NY,NX)-FXVOLI
      FXVOLY=FWO*VOLY(L0,NY,NX)
      VOLY(L1,NY,NX)=VOLY(L1,NY,NX)+FXVOLY
      VOLY(L0,NY,NX)=VOLY(L0,NY,NX)-FXVOLY
      FXVOLWX=FWO*VOLWX(L0,NY,NX)
      VOLWX(L1,NY,NX)=VOLWX(L1,NY,NX)+FXVOLWX
      VOLWX(L0,NY,NX)=VOLWX(L0,NY,NX)-FXVOLWX
      FXVHCM=FWO*VHCM(L0,NY,NX)
      VHCM(L1,NY,NX)=VHCM(L1,NY,NX)+FXVHCM
      VHCM(L0,NY,NX)=VHCM(L0,NY,NX)-FXVHCM
      FXENGY=TKS(L0,NY,NX)*(FXVHCM+4.19*FXVOLW+1.9274*FXVOLI)
      ENGY1=VHCP(L1,NY,NX)*TKS(L1,NY,NX)+FXENGY
      ENGY0=VHCP(L0,NY,NX)*TKS(L0,NY,NX)-FXENGY
      VHCP(L1,NY,NX)=VHCP(L1,NY,NX)
     2+FXVHCM+4.19*(FXVOLW+FXVOLV)+1.9274*FXVOLI
      VHCP(L0,NY,NX)=VHCP(L0,NY,NX)
     2-FXVHCM-4.19*(FXVOLW+FXVOLV)-1.9274*FXVOLI
      IF(VHCP(L1,NY,NX).GT.VHCPNX(NY,NX))THEN
      TKS(L1,NY,NX)=ENGY1/VHCP(L1,NY,NX)
      ELSE
      TKS(L1,NY,NX)=TKS(L,NY,NX)
      ENDIF
      TCS(L1,NY,NX)=TKS(L1,NY,NX)-273.15
      IF(VHCP(L0,NY,NX).GT.VHCPNX(NY,NX))THEN
      TKS(L0,NY,NX)=ENGY0/VHCP(L0,NY,NX)
      ELSE
      TKS(L0,NY,NX)=TKS(L,NY,NX)
      ENDIF
      TCS(L0,NY,NX)=TKS(L0,NY,NX)-273.15
C
C     SOIL FERTILIZER IN BAND, NON-BAND
C
      FXZNH4FA=AMIN1(FX*ZNH4FA(L,NY,NX),ZNH4FA(L0,NY,NX))
      ZNH4FA(L1,NY,NX)=ZNH4FA(L1,NY,NX)+FXZNH4FA
      ZNH4FA(L0,NY,NX)=ZNH4FA(L0,NY,NX)-FXZNH4FA
      FXZNH3FA=AMIN1(FX*ZNH3FA(L,NY,NX),ZNH3FA(L0,NY,NX))
      ZNH3FA(L1,NY,NX)=ZNH3FA(L1,NY,NX)+FXZNH3FA
      ZNH3FA(L0,NY,NX)=ZNH3FA(L0,NY,NX)-FXZNH3FA
      FXZNHUFA=AMIN1(FX*ZNHUFA(L,NY,NX),ZNHUFA(L0,NY,NX))
      ZNHUFA(L1,NY,NX)=ZNHUFA(L1,NY,NX)+FXZNHUFA
      ZNHUFA(L0,NY,NX)=ZNHUFA(L0,NY,NX)-FXZNHUFA
      FXZNO3FA=AMIN1(FX*ZNO3FA(L,NY,NX),ZNO3FA(L0,NY,NX))
      ZNO3FA(L1,NY,NX)=ZNO3FA(L1,NY,NX)+FXZNO3FA
      ZNO3FA(L0,NY,NX)=ZNO3FA(L0,NY,NX)-FXZNO3FA
      FXZNH4FB=AMIN1(FX*ZNH4FB(L,NY,NX),ZNH4FB(L0,NY,NX))
      ZNH4FB(L1,NY,NX)=ZNH4FB(L1,NY,NX)+FXZNH4FB
      ZNH4FB(L0,NY,NX)=ZNH4FB(L0,NY,NX)-FXZNH4FB
      FXZNH3FB=AMIN1(FX*ZNH3FB(L,NY,NX),ZNH3FB(L0,NY,NX))
      ZNH3FB(L1,NY,NX)=ZNH3FB(L1,NY,NX)+FXZNH3FB
      ZNH3FB(L0,NY,NX)=ZNH3FB(L0,NY,NX)-FXZNH3FB
      FXZNHUFB=AMIN1(FX*ZNHUFB(L,NY,NX),ZNHUFB(L0,NY,NX))
      ZNHUFB(L1,NY,NX)=ZNHUFB(L1,NY,NX)+FXZNHUFB
      ZNHUFB(L0,NY,NX)=ZNHUFB(L0,NY,NX)-FXZNHUFB
      FXZNO3FB=AMIN1(FX*ZNO3FB(L,NY,NX),ZNO3FB(L0,NY,NX))
      ZNO3FB(L1,NY,NX)=ZNO3FB(L1,NY,NX)+FXZNO3FB
      ZNO3FB(L0,NY,NX)=ZNO3FB(L0,NY,NX)-FXZNO3FB
C
C     SOIL N,P SOLUTES IN BAND, NON-BAND
C
      FXZNH4S=FWO*ZNH4S(L0,NY,NX)
      ZNH4S(L1,NY,NX)=ZNH4S(L1,NY,NX)+FXZNH4S
      ZNH4S(L0,NY,NX)=ZNH4S(L0,NY,NX)-FXZNH4S
      FXZNH4B=FWO*ZNH4B(L0,NY,NX)
      ZNH4B(L1,NY,NX)=ZNH4B(L1,NY,NX)+FXZNH4B
      ZNH4B(L0,NY,NX)=ZNH4B(L0,NY,NX)-FXZNH4B
      FXZNH3S=FWO*ZNH3S(L0,NY,NX)
      ZNH3S(L1,NY,NX)=ZNH3S(L1,NY,NX)+FXZNH3S
      ZNH3S(L0,NY,NX)=ZNH3S(L0,NY,NX)-FXZNH3S
      FXZNH3B=FWO*ZNH3B(L0,NY,NX)
      ZNH3B(L1,NY,NX)=ZNH3B(L1,NY,NX)+FXZNH3B
      ZNH3B(L0,NY,NX)=ZNH3B(L0,NY,NX)-FXZNH3B
      FXZNO3S=FWO*ZNO3S(L0,NY,NX)
      ZNO3S(L1,NY,NX)=ZNO3S(L1,NY,NX)+FXZNO3S
      ZNO3S(L0,NY,NX)=ZNO3S(L0,NY,NX)-FXZNO3S
      FXZNO3B=FWO*ZNO3B(L0,NY,NX)
      ZNO3B(L1,NY,NX)=ZNO3B(L1,NY,NX)+FXZNO3B
      ZNO3B(L0,NY,NX)=ZNO3B(L0,NY,NX)-FXZNO3B
      FXZNO2S=FWO*ZNO2S(L0,NY,NX)
      ZNO2S(L1,NY,NX)=ZNO2S(L1,NY,NX)+FXZNO2S
      ZNO2S(L0,NY,NX)=ZNO2S(L0,NY,NX)-FXZNO2S
      FXZNO2B=FWO*ZNO2B(L0,NY,NX)
      ZNO2B(L1,NY,NX)=ZNO2B(L1,NY,NX)+FXZNO2B
      ZNO2B(L0,NY,NX)=ZNO2B(L0,NY,NX)-FXZNO2B
      FXH1PO4=FWO*H1PO4(L0,NY,NX)
      H1PO4(L1,NY,NX)=H1PO4(L1,NY,NX)+FXH1PO4
      H1PO4(L0,NY,NX)=H1PO4(L0,NY,NX)-FXH1PO4
      FXH2PO4=FWO*H2PO4(L0,NY,NX)
      H2PO4(L1,NY,NX)=H2PO4(L1,NY,NX)+FXH2PO4
      H2PO4(L0,NY,NX)=H2PO4(L0,NY,NX)-FXH2PO4
C
C     SOIL SALTS IN NON-BAND
C
      IF(ISALTG.NE.0)THEN
      FXZAL=FWO*ZAL(L0,NY,NX)
      ZAL(L1,NY,NX)=ZAL(L1,NY,NX)+FXZAL
      ZAL(L0,NY,NX)=ZAL(L0,NY,NX)-FXZAL
      FXZFE=FWO*ZFE(L0,NY,NX)
      ZFE(L1,NY,NX)=ZFE(L1,NY,NX)+FXZFE
      ZFE(L0,NY,NX)=ZFE(L0,NY,NX)-FXZFE
      FXZHY=FWO*ZHY(L0,NY,NX)
      ZHY(L1,NY,NX)=ZHY(L1,NY,NX)+FXZHY
      ZHY(L0,NY,NX)=ZHY(L0,NY,NX)-FXZHY
      FXZCA=FWO*ZCA(L0,NY,NX)
      ZCA(L1,NY,NX)=ZCA(L1,NY,NX)+FXZCA
      ZCA(L0,NY,NX)=ZCA(L0,NY,NX)-FXZCA
      FXZMG=FWO*ZMG(L0,NY,NX)
      ZMG(L1,NY,NX)=ZMG(L1,NY,NX)+FXZMG
      ZMG(L0,NY,NX)=ZMG(L0,NY,NX)-FXZMG
      FXZNA=FWO*ZNA(L0,NY,NX)
      ZNA(L1,NY,NX)=ZNA(L1,NY,NX)+FXZNA
      ZNA(L0,NY,NX)=ZNA(L0,NY,NX)-FXZNA
      FXZKA=FWO*ZKA(L0,NY,NX)
      ZKA(L1,NY,NX)=ZKA(L1,NY,NX)+FXZKA
      ZKA(L0,NY,NX)=ZKA(L0,NY,NX)-FXZKA
      FXZOH=FWO*ZOH(L0,NY,NX)
      ZOH(L1,NY,NX)=ZOH(L1,NY,NX)+FXZOH
      ZOH(L0,NY,NX)=ZOH(L0,NY,NX)-FXZOH
      FXZSO4=FWO*ZSO4(L0,NY,NX)
      ZSO4(L1,NY,NX)=ZSO4(L1,NY,NX)+FXZSO4
      ZSO4(L0,NY,NX)=ZSO4(L0,NY,NX)-FXZSO4
      FXZCL=FWO*ZCL(L0,NY,NX)
      ZCL(L1,NY,NX)=ZCL(L1,NY,NX)+FXZCL
      ZCL(L0,NY,NX)=ZCL(L0,NY,NX)-FXZCL
      FXZCO3=FWO*ZCO3(L0,NY,NX)
      ZCO3(L1,NY,NX)=ZCO3(L1,NY,NX)+FXZCO3
      ZCO3(L0,NY,NX)=ZCO3(L0,NY,NX)-FXZCO3
      FXZHCO3=FWO*ZHCO3(L0,NY,NX)
      ZHCO3(L1,NY,NX)=ZHCO3(L1,NY,NX)+FXZHCO3
      ZHCO3(L0,NY,NX)=ZHCO3(L0,NY,NX)-FXZHCO3
      FXZALOH1=FWO*ZALOH1(L0,NY,NX)
      ZALOH1(L1,NY,NX)=ZALOH1(L1,NY,NX)+FXZALOH1
      ZALOH1(L0,NY,NX)=ZALOH1(L0,NY,NX)-FXZALOH1
      FXZALOH2=FWO*ZALOH2(L0,NY,NX)
      ZALOH2(L1,NY,NX)=ZALOH2(L1,NY,NX)+FXZALOH2
      ZALOH2(L0,NY,NX)=ZALOH2(L0,NY,NX)-FXZALOH2
      FXZALOH3=FWO*ZALOH3(L0,NY,NX)
      ZALOH3(L1,NY,NX)=ZALOH3(L1,NY,NX)+FXZALOH3
      ZALOH3(L0,NY,NX)=ZALOH3(L0,NY,NX)-FXZALOH3
      FXZALOH4=FWO*ZALOH4(L0,NY,NX)
      ZALOH4(L1,NY,NX)=ZALOH4(L1,NY,NX)+FXZALOH4
      ZALOH4(L0,NY,NX)=ZALOH4(L0,NY,NX)-FXZALOH4
      FXZALS=FWO*ZALS(L0,NY,NX)
      ZALS(L1,NY,NX)=ZALS(L1,NY,NX)+FXZALS
      ZALS(L0,NY,NX)=ZALS(L0,NY,NX)-FXZALS
      FXZFEOH1=FWO*ZFEOH1(L0,NY,NX)
      ZFEOH1(L1,NY,NX)=ZFEOH1(L1,NY,NX)+FXZFEOH1
      ZFEOH1(L0,NY,NX)=ZFEOH1(L0,NY,NX)-FXZFEOH1
      FXZFEOH2=FWO*ZFEOH2(L0,NY,NX)
      ZFEOH2(L1,NY,NX)=ZFEOH2(L1,NY,NX)+FXZFEOH2
      ZFEOH2(L0,NY,NX)=ZFEOH2(L0,NY,NX)-FXZFEOH2
      FXZFEOH3=FWO*ZFEOH3(L0,NY,NX)
      ZFEOH3(L1,NY,NX)=ZFEOH3(L1,NY,NX)+FXZFEOH3
      ZFEOH3(L0,NY,NX)=ZFEOH3(L0,NY,NX)-FXZFEOH3
      FXZFEOH4=FWO*ZFEOH4(L0,NY,NX)
      ZFEOH4(L1,NY,NX)=ZFEOH4(L1,NY,NX)+FXZFEOH4
      ZFEOH4(L0,NY,NX)=ZFEOH4(L0,NY,NX)-FXZFEOH4
      FXZFES=FWO*ZFES(L0,NY,NX)
      ZFES(L1,NY,NX)=ZFES(L1,NY,NX)+FXZFES
      ZFES(L0,NY,NX)=ZFES(L0,NY,NX)-FXZFES
      FXZCAO=FWO*ZCAO(L0,NY,NX)
      ZCAO(L1,NY,NX)=ZCAO(L1,NY,NX)+FXZCAO
      ZCAO(L0,NY,NX)=ZCAO(L0,NY,NX)-FXZCAO
      FXZCAC=FWO*ZCAC(L0,NY,NX)
      ZCAC(L1,NY,NX)=ZCAC(L1,NY,NX)+FXZCAC
      ZCAC(L0,NY,NX)=ZCAC(L0,NY,NX)-FXZCAC
      FXZCAH=FWO*ZCAH(L0,NY,NX)
      ZCAH(L1,NY,NX)=ZCAH(L1,NY,NX)+FXZCAH
      ZCAH(L0,NY,NX)=ZCAH(L0,NY,NX)-FXZCAH
      FXZCAS=FWO*ZCAS(L0,NY,NX)
      ZCAS(L1,NY,NX)=ZCAS(L1,NY,NX)+FXZCAS
      ZCAS(L0,NY,NX)=ZCAS(L0,NY,NX)-FXZCAS
      FXZMGO=FWO*ZMGO(L0,NY,NX)
      ZMGO(L1,NY,NX)=ZMGO(L1,NY,NX)+FXZMGO
      ZMGO(L0,NY,NX)=ZMGO(L0,NY,NX)-FXZMGO
      FXZMGC=FWO*ZMGC(L0,NY,NX)
      ZMGC(L1,NY,NX)=ZMGC(L1,NY,NX)+FXZMGC
      ZMGC(L0,NY,NX)=ZMGC(L0,NY,NX)-FXZMGC
      FXZMGH=FWO*ZMGH(L0,NY,NX)
      ZMGH(L1,NY,NX)=ZMGH(L1,NY,NX)+FXZMGH
      ZMGH(L0,NY,NX)=ZMGH(L0,NY,NX)-FXZMGH
      FXZMGS=FWO*ZMGS(L0,NY,NX)
      ZMGS(L1,NY,NX)=ZMGS(L1,NY,NX)+FXZMGS
      ZMGS(L0,NY,NX)=ZMGS(L0,NY,NX)-FXZMGS
      FXZNAC=FWO*ZNAC(L0,NY,NX)
      ZNAC(L1,NY,NX)=ZNAC(L1,NY,NX)+FXZNAC
      ZNAC(L0,NY,NX)=ZNAC(L0,NY,NX)-FXZNAC
      FXZNAS=FWO*ZNAS(L0,NY,NX)
      ZNAS(L1,NY,NX)=ZNAS(L1,NY,NX)+FXZNAS
      ZNAS(L0,NY,NX)=ZNAS(L0,NY,NX)-FXZNAS
      FXZKAS=FWO*ZKAS(L0,NY,NX)
      ZKAS(L1,NY,NX)=ZKAS(L1,NY,NX)+FXZKAS
      ZKAS(L0,NY,NX)=ZKAS(L0,NY,NX)-FXZKAS
      FXH0PO4=FWO*H0PO4(L0,NY,NX)
      H0PO4(L1,NY,NX)=H0PO4(L1,NY,NX)+FXH0PO4
      H0PO4(L0,NY,NX)=H0PO4(L0,NY,NX)-FXH0PO4
      FXH3PO4=FWO*H3PO4(L0,NY,NX)
      H3PO4(L1,NY,NX)=H3PO4(L1,NY,NX)+FXH3PO4
      H3PO4(L0,NY,NX)=H3PO4(L0,NY,NX)-FXH3PO4
      FXZFE1P=FWO*ZFE1P(L0,NY,NX)
      ZFE1P(L1,NY,NX)=ZFE1P(L1,NY,NX)+FXZFE1P
      ZFE1P(L0,NY,NX)=ZFE1P(L0,NY,NX)-FXZFE1P
      FXZFE2P=FWO*ZFE2P(L0,NY,NX)
      ZFE2P(L1,NY,NX)=ZFE2P(L1,NY,NX)+FXZFE2P
      ZFE2P(L0,NY,NX)=ZFE2P(L0,NY,NX)-FXZFE2P
      FXZCA0P=FWO*ZCA0P(L0,NY,NX)
      ZCA0P(L1,NY,NX)=ZCA0P(L1,NY,NX)+FXZCA0P
      ZCA0P(L0,NY,NX)=ZCA0P(L0,NY,NX)-FXZCA0P
      FXZCA1P=FWO*ZCA1P(L0,NY,NX)
      ZCA1P(L1,NY,NX)=ZCA1P(L1,NY,NX)+FXZCA1P
      ZCA1P(L0,NY,NX)=ZCA1P(L0,NY,NX)-FXZCA1P
      FXZCA2P=FWO*ZCA2P(L0,NY,NX)
      ZCA2P(L1,NY,NX)=ZCA2P(L1,NY,NX)+FXZCA2P
      ZCA2P(L0,NY,NX)=ZCA2P(L0,NY,NX)-FXZCA2P
      FXZMG1P=FWO*ZMG1P(L0,NY,NX)
      ZMG1P(L1,NY,NX)=ZMG1P(L1,NY,NX)+FXZMG1P
      ZMG1P(L0,NY,NX)=ZMG1P(L0,NY,NX)-FXZMG1P
      ENDIF
C
C     SOIL SOLUTES IN BAND
C
      IF(L0.NE.0)THEN
      FXH1POB=FWO*H1POB(L0,NY,NX)
      H1POB(L1,NY,NX)=H1POB(L1,NY,NX)+FXH1POB
      H1POB(L0,NY,NX)=H1POB(L0,NY,NX)-FXH1POB
      FXH2POB=FWO*H2POB(L0,NY,NX)
      H2POB(L1,NY,NX)=H2POB(L1,NY,NX)+FXH2POB
      H2POB(L0,NY,NX)=H2POB(L0,NY,NX)-FXH2POB

      IF(ISALTG.NE.0)THEN
      FXH0POB=FWO*H0POB(L0,NY,NX)
      H0POB(L1,NY,NX)=H0POB(L1,NY,NX)+FXH0POB
      H0POB(L0,NY,NX)=H0POB(L0,NY,NX)-FXH0POB
      FXH3POB=FWO*H3POB(L0,NY,NX)
      H3POB(L1,NY,NX)=H3POB(L1,NY,NX)+FXH3POB
      H3POB(L0,NY,NX)=H3POB(L0,NY,NX)-FXH3POB
      FXZFE1PB=FWO*ZFE1PB(L0,NY,NX)
      ZFE1PB(L1,NY,NX)=ZFE1PB(L1,NY,NX)+FXZFE1PB
      ZFE1PB(L0,NY,NX)=ZFE1PB(L0,NY,NX)-FXZFE1PB
      FXZFE2PB=FWO*ZFE2PB(L0,NY,NX)
      ZFE2PB(L1,NY,NX)=ZFE2PB(L1,NY,NX)+FXZFE2PB
      ZFE2PB(L0,NY,NX)=ZFE2PB(L0,NY,NX)-FXZFE2PB
      FXZCA0PB=FWO*ZCA0PB(L0,NY,NX)
      ZCA0PB(L1,NY,NX)=ZCA0PB(L1,NY,NX)+FXZCA0PB
      ZCA0PB(L0,NY,NX)=ZCA0PB(L0,NY,NX)-FXZCA0PB
      FXZCA1PB=FWO*ZCA1PB(L0,NY,NX)
      ZCA1PB(L1,NY,NX)=ZCA1PB(L1,NY,NX)+FXZCA1PB
      ZCA1PB(L0,NY,NX)=ZCA1PB(L0,NY,NX)-FXZCA1PB
      FXZCA2PB=FWO*ZCA2PB(L0,NY,NX)
      ZCA2PB(L1,NY,NX)=ZCA2PB(L1,NY,NX)+FXZCA2PB
      ZCA2PB(L0,NY,NX)=ZCA2PB(L0,NY,NX)-FXZCA2PB
      FXZMG1PB=FWO*ZMG1PB(L0,NY,NX)
      ZMG1PB(L1,NY,NX)=ZMG1PB(L1,NY,NX)+FXZMG1PB
      ZMG1PB(L0,NY,NX)=ZMG1PB(L0,NY,NX)-FXZMG1PB
      ENDIF
C
C     SOIL ADSORBED CATIONS IN BAND, NON-BAND
C
C     IF(L0.EQ.L.OR.CEC(L0,NY,NX).LE.ZERO)THEN
      FCO=FX
C     ELSE
C     FCO=AMIN1(1.0,FX*CEC(L1,NY,NX)/CEC(L0,NY,NX))
C     ENDIF
      FXXCEC=FCO*XCEC(L0,NY,NX)
      XCEC(L1,NY,NX)=XCEC(L1,NY,NX)+FXXCEC
      XCEC(L0,NY,NX)=XCEC(L0,NY,NX)-FXXCEC
      FXXN4=FCO*XN4(L0,NY,NX)
      XN4(L1,NY,NX)=XN4(L1,NY,NX)+FXXN4
      XN4(L0,NY,NX)=XN4(L0,NY,NX)-FXXN4
      FXXNB=FCO*XNB(L0,NY,NX)
      XNB(L1,NY,NX)=XNB(L1,NY,NX)+FXXNB
      XNB(L0,NY,NX)=XNB(L0,NY,NX)-FXXNB
      FXXHY=FCO*XHY(L0,NY,NX)
      XHY(L1,NY,NX)=XHY(L1,NY,NX)+FXXHY
      XHY(L0,NY,NX)=XHY(L0,NY,NX)-FXXHY
      FXXAL=FCO*XAL(L0,NY,NX)
      XAL(L1,NY,NX)=XAL(L1,NY,NX)+FXXAL
      XAL(L0,NY,NX)=XAL(L0,NY,NX)-FXXAL
      FXXFE=FCO*XFE(L0,NY,NX)
      XFE(L1,NY,NX)=XFE(L1,NY,NX)+FXXFE
      XFE(L0,NY,NX)=XFE(L0,NY,NX)-FXXFE
      FXXCA=FCO*XCA(L0,NY,NX)
      XCA(L1,NY,NX)=XCA(L1,NY,NX)+FXXCA
      XCA(L0,NY,NX)=XCA(L0,NY,NX)-FXXCA
      FXXMG=FCO*XMG(L0,NY,NX)
      XMG(L1,NY,NX)=XMG(L1,NY,NX)+FXXMG
      XMG(L0,NY,NX)=XMG(L0,NY,NX)-FXXMG
      FXXNA=FCO*XNA(L0,NY,NX)
      XNA(L1,NY,NX)=XNA(L1,NY,NX)+FXXNA
      XNA(L0,NY,NX)=XNA(L0,NY,NX)-FXXNA
      FXXKA=FCO*XKA(L0,NY,NX)
      XKA(L1,NY,NX)=XKA(L1,NY,NX)+FXXKA
      XKA(L0,NY,NX)=XKA(L0,NY,NX)-FXXKA
      FXXHC=FCO*XHC(L0,NY,NX)
      XHC(L1,NY,NX)=XHC(L1,NY,NX)+FXXHC
      XHC(L0,NY,NX)=XHC(L0,NY,NX)-FXXHC
      FXXALO2=FCO*XALO2(L0,NY,NX)
      XALO2(L1,NY,NX)=XALO2(L1,NY,NX)+FXXALO2
      XALO2(L0,NY,NX)=XALO2(L0,NY,NX)-FXXALO2
      FXXFEO2=FCO*XFEO2(L0,NY,NX)
      XFEO2(L1,NY,NX)=XFEO2(L1,NY,NX)+FXXFEO2
      XFEO2(L0,NY,NX)=XFEO2(L0,NY,NX)-FXXFEO2
C
C     SOIL ADSORBED ANIONS IN BAND, NON-BAND
C
C     IF(L0.EQ.L.OR.AEC(L0,NY,NX).LE.ZERO)THEN
      FAO=FX
C     ELSE
C     FAO=AMIN1(1.0,FX*AEC(L1,NY,NX)/AEC(L0,NY,NX))
C     ENDIF
      FXXAEC=FAO*XAEC(L0,NY,NX)
      XAEC(L1,NY,NX)=XAEC(L1,NY,NX)+FXXAEC
      XAEC(L0,NY,NX)=XAEC(L0,NY,NX)-FXXAEC
      FXXOH0=FAO*XOH0(L0,NY,NX)
      XOH0(L1,NY,NX)=XOH0(L1,NY,NX)+FXXOH0
      XOH0(L0,NY,NX)=XOH0(L0,NY,NX)-FXXOH0
      FXXOH1=FAO*XOH1(L0,NY,NX)
      XOH1(L1,NY,NX)=XOH1(L1,NY,NX)+FXXOH1
      XOH1(L0,NY,NX)=XOH1(L0,NY,NX)-FXXOH1
      FXXOH2=FAO*XOH2(L0,NY,NX)
      XOH2(L1,NY,NX)=XOH2(L1,NY,NX)+FXXOH2
      XOH2(L0,NY,NX)=XOH2(L0,NY,NX)-FXXOH2
      FXXH1P=FAO*XH1P(L0,NY,NX)
      XH1P(L1,NY,NX)=XH1P(L1,NY,NX)+FXXH1P
      XH1P(L0,NY,NX)=XH1P(L0,NY,NX)-FXXH1P
      FXXH2P=FAO*XH2P(L0,NY,NX)
      XH2P(L1,NY,NX)=XH2P(L1,NY,NX)+FXXH2P
      XH2P(L0,NY,NX)=XH2P(L0,NY,NX)-FXXH2P
      FXXOH0B=FAO*XOH0B(L0,NY,NX)
      XOH0B(L1,NY,NX)=XOH0B(L1,NY,NX)+FXXOH0B
      XOH0B(L0,NY,NX)=XOH0B(L0,NY,NX)-FXXOH0B
      FXXOH1B=FAO*XOH1B(L0,NY,NX)
      XOH1B(L1,NY,NX)=XOH1B(L1,NY,NX)+FXXOH1B
      XOH1B(L0,NY,NX)=XOH1B(L0,NY,NX)-FXXOH1B
      FXXOH2B=FAO*XOH2B(L0,NY,NX)
      XOH2B(L1,NY,NX)=XOH2B(L1,NY,NX)+FXXOH2B
      XOH2B(L0,NY,NX)=XOH2B(L0,NY,NX)-FXXOH2B
      FXXH1PB=FAO*XH1PB(L0,NY,NX)
      XH1PB(L1,NY,NX)=XH1PB(L1,NY,NX)+FXXH1PB
      XH1PB(L0,NY,NX)=XH1PB(L0,NY,NX)-FXXH1PB
      FXXH2PB=FAO*XH2PB(L0,NY,NX)
      XH2PB(L1,NY,NX)=XH2PB(L1,NY,NX)+FXXH2PB
      XH2PB(L0,NY,NX)=XH2PB(L0,NY,NX)-FXXH2PB
C
C     SOIL PRECIPITATES IN BAND, NON-BAND
C
      FXPALOH=AMIN1(FX*PALOH(L,NY,NX),PALOH(L0,NY,NX))
      PALOH(L1,NY,NX)=PALOH(L1,NY,NX)+FXPALOH
      PALOH(L0,NY,NX)=PALOH(L0,NY,NX)-FXPALOH
      FXPFEOH=AMIN1(FX*PFEOH(L,NY,NX),PFEOH(L0,NY,NX))
      PFEOH(L1,NY,NX)=PFEOH(L1,NY,NX)+FXPFEOH
      PFEOH(L0,NY,NX)=PFEOH(L0,NY,NX)-FXPFEOH
      FXPCACO=AMIN1(FX*PCACO(L,NY,NX),PCACO(L0,NY,NX))
      PCACO(L1,NY,NX)=PCACO(L1,NY,NX)+FXPCACO
      PCACO(L0,NY,NX)=PCACO(L0,NY,NX)-FXPCACO
      FXPCASO=AMIN1(FX*PCASO(L,NY,NX),PCASO(L0,NY,NX))
      PCASO(L1,NY,NX)=PCASO(L1,NY,NX)+FXPCASO
      PCASO(L0,NY,NX)=PCASO(L0,NY,NX)-FXPCASO
      FXPALPO=AMIN1(FX*PALPO(L,NY,NX),PALPO(L0,NY,NX))
      PALPO(L1,NY,NX)=PALPO(L1,NY,NX)+FXPALPO
      PALPO(L0,NY,NX)=PALPO(L0,NY,NX)-FXPALPO
      FXPFEPO=AMIN1(FX*PFEPO(L,NY,NX),PFEPO(L0,NY,NX))
      PFEPO(L1,NY,NX)=PFEPO(L1,NY,NX)+FXPFEPO
      PFEPO(L0,NY,NX)=PFEPO(L0,NY,NX)-FXPFEPO
      FXPCAPD=AMIN1(FX*PCAPD(L,NY,NX),PCAPD(L0,NY,NX))
      PCAPD(L1,NY,NX)=PCAPD(L1,NY,NX)+FXPCAPD
      PCAPD(L0,NY,NX)=PCAPD(L0,NY,NX)-FXPCAPD
      FXPCAPH=AMIN1(FX*PCAPH(L,NY,NX),PCAPH(L0,NY,NX))
      PCAPH(L1,NY,NX)=PCAPH(L1,NY,NX)+FXPCAPH
      PCAPH(L0,NY,NX)=PCAPH(L0,NY,NX)-FXPCAPH
      FXPCAPM=AMIN1(FX*PCAPM(L,NY,NX),PCAPM(L0,NY,NX))
      PCAPM(L1,NY,NX)=PCAPM(L1,NY,NX)+FXPCAPM
      PCAPM(L0,NY,NX)=PCAPM(L0,NY,NX)-FXPCAPM
      FXPALPB=AMIN1(FX*PALPB(L,NY,NX),PALPB(L0,NY,NX))
      PALPB(L1,NY,NX)=PALPB(L1,NY,NX)+FXPALPB
      PALPB(L0,NY,NX)=PALPB(L0,NY,NX)-FXPALPB
      FXPFEPB=AMIN1(FX*PFEPB(L,NY,NX),PFEPB(L0,NY,NX))
      PFEPB(L1,NY,NX)=PFEPB(L1,NY,NX)+FXPFEPB
      PFEPB(L0,NY,NX)=PFEPB(L0,NY,NX)-FXPFEPB
      FXPCPDB=AMIN1(FX*PCPDB(L,NY,NX),PCPDB(L0,NY,NX))
      PCPDB(L1,NY,NX)=PCPDB(L1,NY,NX)+FXPCPDB
      PCPDB(L0,NY,NX)=PCPDB(L0,NY,NX)-FXPCPDB
      FXPCPHB=AMIN1(FX*PCPHB(L,NY,NX),PCPHB(L0,NY,NX))
      PCPHB(L1,NY,NX)=PCPHB(L1,NY,NX)+FXPCPHB
      PCPHB(L0,NY,NX)=PCPHB(L0,NY,NX)-FXPCPHB
      FXPCPMB=AMIN1(FX*PCPMB(L,NY,NX),PCPMB(L0,NY,NX))
      PCPMB(L1,NY,NX)=PCPMB(L1,NY,NX)+FXPCPMB
      PCPMB(L0,NY,NX)=PCPMB(L0,NY,NX)-FXPCPMB
C
C     SOIL GASEOUS GASES
C
      FXCO2G=FWO*CO2G(L0,NY,NX)
      CO2G(L1,NY,NX)=CO2G(L1,NY,NX)+FXCO2G
      CO2G(L0,NY,NX)=CO2G(L0,NY,NX)-FXCO2G
      FXCH4G=FWO*CH4G(L0,NY,NX)
      CH4G(L1,NY,NX)=CH4G(L1,NY,NX)+FXCH4G
      CH4G(L0,NY,NX)=CH4G(L0,NY,NX)-FXCH4G
      FXOXYG=FWO*OXYG(L0,NY,NX)
      OXYG(L1,NY,NX)=OXYG(L1,NY,NX)+FXOXYG
      OXYG(L0,NY,NX)=OXYG(L0,NY,NX)-FXOXYG
      FXZ2GG=FWO*Z2GG(L0,NY,NX)
      Z2GG(L1,NY,NX)=Z2GG(L1,NY,NX)+FXZ2GG
      Z2GG(L0,NY,NX)=Z2GG(L0,NY,NX)-FXZ2GG
      FXZ2OG=FWO*Z2OG(L0,NY,NX)
      Z2OG(L1,NY,NX)=Z2OG(L1,NY,NX)+FXZ2OG
      Z2OG(L0,NY,NX)=Z2OG(L0,NY,NX)-FXZ2OG
      FXZNH3G=FWO*ZNH3G(L0,NY,NX)
      ZNH3G(L1,NY,NX)=ZNH3G(L1,NY,NX)+FXZNH3G
      ZNH3G(L0,NY,NX)=ZNH3G(L0,NY,NX)-FXZNH3G
      FXH2GG=FWO*H2GG(L0,NY,NX)
      H2GG(L1,NY,NX)=H2GG(L1,NY,NX)+FXH2GG
      H2GG(L0,NY,NX)=H2GG(L0,NY,NX)-FXH2GG
      ENDIF
C
C     SOIL AQUEOUS GASES
C
      IF(L0.NE.0)THEN
      FXCO2S=FWO*CO2S(L0,NY,NX)
      CO2S(L1,NY,NX)=CO2S(L1,NY,NX)+FXCO2S
      CO2S(L0,NY,NX)=CO2S(L0,NY,NX)-FXCO2S
      FXCH4S=FWO*CH4S(L0,NY,NX)
      CH4S(L1,NY,NX)=CH4S(L1,NY,NX)+FXCH4S
      CH4S(L0,NY,NX)=CH4S(L0,NY,NX)-FXCH4S
      FXOXYS=FWO*OXYS(L0,NY,NX)
      OXYS(L1,NY,NX)=OXYS(L1,NY,NX)+FXOXYS
      OXYS(L0,NY,NX)=OXYS(L0,NY,NX)-FXOXYS
      FXZ2GS=FWO*Z2GS(L0,NY,NX)
      Z2GS(L1,NY,NX)=Z2GS(L1,NY,NX)+FXZ2GS
      Z2GS(L0,NY,NX)=Z2GS(L0,NY,NX)-FXZ2GS
      FXZ2OS=FWO*Z2OS(L0,NY,NX)
      Z2OS(L1,NY,NX)=Z2OS(L1,NY,NX)+FXZ2OS
      Z2OS(L0,NY,NX)=Z2OS(L0,NY,NX)-FXZ2OS
      FXH2GS=FWO*H2GS(L0,NY,NX)
      H2GS(L1,NY,NX)=H2GS(L1,NY,NX)+FXH2GS
      H2GS(L0,NY,NX)=H2GS(L0,NY,NX)-FXH2GS
C
C     SOIL MACROPORE N,P SOLUTES
C
      IF(FHOL(L1,NY,NX).GT.ZERO.AND.FHOL(L0,NY,NX).GT.ZERO)THEN
      FXZNH4SH=FHO*ZNH4SH(L0,NY,NX)
      ZNH4SH(L1,NY,NX)=ZNH4SH(L1,NY,NX)+FXZNH4SH
      ZNH4SH(L0,NY,NX)=ZNH4SH(L0,NY,NX)-FXZNH4SH
      FXZNH3SH=FHO*ZNH3SH(L0,NY,NX)
      ZNH3SH(L1,NY,NX)=ZNH3SH(L1,NY,NX)+FXZNH3SH
      ZNH3SH(L0,NY,NX)=ZNH3SH(L0,NY,NX)-FXZNH3SH
      FXZNO3SH=FHO*ZNO3SH(L0,NY,NX)
      ZNO3SH(L1,NY,NX)=ZNO3SH(L1,NY,NX)+FXZNO3SH
      ZNO3SH(L0,NY,NX)=ZNO3SH(L0,NY,NX)-FXZNO3SH
      FXZNO2SH=FHO*ZNO2SH(L0,NY,NX)
      ZNO2SH(L1,NY,NX)=ZNO2SH(L1,NY,NX)+FXZNO2SH
      ZNO2SH(L0,NY,NX)=ZNO2SH(L0,NY,NX)-FXZNO2SH
      FXZNH4BH=FHO*ZNH4BH(L0,NY,NX)
      ZNH4BH(L1,NY,NX)=ZNH4BH(L1,NY,NX)+FXZNH4BH
      ZNH4BH(L0,NY,NX)=ZNH4BH(L0,NY,NX)-FXZNH4BH
      FXZNH3BH=FHO*ZNH3BH(L0,NY,NX)
      ZNH3BH(L1,NY,NX)=ZNH3BH(L1,NY,NX)+FXZNH3BH
      ZNH3BH(L0,NY,NX)=ZNH3BH(L0,NY,NX)-FXZNH3BH
      FXZNO3BH=FHO*ZNO3BH(L0,NY,NX)
      ZNO3BH(L1,NY,NX)=ZNO3BH(L1,NY,NX)+FXZNO3BH
      ZNO3BH(L0,NY,NX)=ZNO3BH(L0,NY,NX)-FXZNO3BH
      FXZNO2BH=FHO*ZNO2BH(L0,NY,NX)
      ZNO2BH(L1,NY,NX)=ZNO2BH(L1,NY,NX)+FXZNO2BH
      ZNO2BH(L0,NY,NX)=ZNO2BH(L0,NY,NX)-FXZNO2BH
      FXH1PO4H=FHO*H1PO4H(L0,NY,NX)
      H1PO4H(L1,NY,NX)=H1PO4H(L1,NY,NX)+FXH1PO4H
      H1PO4H(L0,NY,NX)=H1PO4H(L0,NY,NX)-FXH1PO4H
      FXH2PO4H=FHO*H2PO4H(L0,NY,NX)
      H2PO4H(L1,NY,NX)=H2PO4H(L1,NY,NX)+FXH2PO4H
      H2PO4H(L0,NY,NX)=H2PO4H(L0,NY,NX)-FXH2PO4H
      FXH1POBH=FHO*H1POBH(L0,NY,NX)
      H1POBH(L1,NY,NX)=H1POBH(L1,NY,NX)+FXH1POBH
      H1POBH(L0,NY,NX)=H1POBH(L0,NY,NX)-FXH1POBH
      FXH2POBH=FHO*H2POBH(L0,NY,NX)
      H2POBH(L1,NY,NX)=H2POBH(L1,NY,NX)+FXH2POBH
      H2POBH(L0,NY,NX)=H2POBH(L0,NY,NX)-FXH2POBH
C
C     SOIL MACROPORE SOLUBLE SALTS
C
      IF(ISALTG.NE.0)THEN
      FXZALH=FHO*ZALH(L0,NY,NX)
      ZALH(L1,NY,NX)=ZALH(L1,NY,NX)+FXZALH
      ZALH(L0,NY,NX)=ZALH(L0,NY,NX)-FXZALH
      FXZFEH=FHO*ZFEH(L0,NY,NX)
      ZFEH(L1,NY,NX)=ZFEH(L1,NY,NX)+FXZFEH
      ZFEH(L0,NY,NX)=ZFEH(L0,NY,NX)-FXZFEH
      FXZHYH=FHO*ZHYH(L0,NY,NX)
      ZHYH(L1,NY,NX)=ZHYH(L1,NY,NX)+FXZHYH
      ZHYH(L0,NY,NX)=ZHYH(L0,NY,NX)-FXZHYH
      FXZCCH=FHO*ZCCH(L0,NY,NX)
      ZCCH(L1,NY,NX)=ZCCH(L1,NY,NX)+FXZCCH
      ZCCH(L0,NY,NX)=ZCCH(L0,NY,NX)-FXZCCH
      FXZMAH=FHO*ZMAH(L0,NY,NX)
      ZMAH(L1,NY,NX)=ZMAH(L1,NY,NX)+FXZMAH
      ZMAH(L0,NY,NX)=ZMAH(L0,NY,NX)-FXZMAH
      FXZNAH=FHO*ZNAH(L0,NY,NX)
      ZNAH(L1,NY,NX)=ZNAH(L1,NY,NX)+FXZNAH
      ZNAH(L0,NY,NX)=ZNAH(L0,NY,NX)-FXZNAH
      FXZKAH=FHO*ZKAH(L0,NY,NX)
      ZKAH(L1,NY,NX)=ZKAH(L1,NY,NX)+FXZKAH
      ZKAH(L0,NY,NX)=ZKAH(L0,NY,NX)-FXZKAH
      FXZOHH=FHO*ZOHH(L0,NY,NX)
      ZOHH(L1,NY,NX)=ZOHH(L1,NY,NX)+FXZOHH
      ZOHH(L0,NY,NX)=ZOHH(L0,NY,NX)-FXZOHH
      FXZSO4H=FHO*ZSO4H(L0,NY,NX)
      ZSO4H(L1,NY,NX)=ZSO4H(L1,NY,NX)+FXZSO4H
      ZSO4H(L0,NY,NX)=ZSO4H(L0,NY,NX)-FXZSO4H
      FXZCLH=FHO*ZCLH(L0,NY,NX)
      ZCLH(L1,NY,NX)=ZCLH(L1,NY,NX)+FXZCLH
      ZCLH(L0,NY,NX)=ZCLH(L0,NY,NX)-FXZCLH
      FXZCO3H=FHO*ZCO3H(L0,NY,NX)
      ZCO3H(L1,NY,NX)=ZCO3H(L1,NY,NX)+FXZCO3H
      ZCO3H(L0,NY,NX)=ZCO3H(L0,NY,NX)-FXZCO3H
      FXZHCO3H=FHO*ZHCO3H(L0,NY,NX)
      ZHCO3H(L1,NY,NX)=ZHCO3H(L1,NY,NX)+FXZHCO3H
      ZHCO3H(L0,NY,NX)=ZHCO3H(L0,NY,NX)-FXZHCO3H
      FXZALO1H=FHO*ZALO1H(L0,NY,NX)
      ZALO1H(L1,NY,NX)=ZALO1H(L1,NY,NX)+FXZALO1H
      ZALO1H(L0,NY,NX)=ZALO1H(L0,NY,NX)-FXZALO1H
      FXZALO2H=FHO*ZALO2H(L0,NY,NX)
      ZALO2H(L1,NY,NX)=ZALO2H(L1,NY,NX)+FXZALO2H
      ZALO2H(L0,NY,NX)=ZALO2H(L0,NY,NX)-FXZALO2H
      FXZALO3H=FHO*ZALO3H(L0,NY,NX)
      ZALO3H(L1,NY,NX)=ZALO3H(L1,NY,NX)+FXZALO3H
      ZALO3H(L0,NY,NX)=ZALO3H(L0,NY,NX)-FXZALO3H
      FXZALO4H=FHO*ZALO4H(L0,NY,NX)
      ZALO4H(L1,NY,NX)=ZALO4H(L1,NY,NX)+FXZALO4H
      ZALO4H(L0,NY,NX)=ZALO4H(L0,NY,NX)-FXZALO4H
      FXZALSH=FHO*ZALSH(L0,NY,NX)
      ZALSH(L1,NY,NX)=ZALSH(L1,NY,NX)+FXZALSH
      ZALSH(L0,NY,NX)=ZALSH(L0,NY,NX)-FXZALSH
      FXZFEO1H=FHO*ZFEO1H(L0,NY,NX)
      ZFEO1H(L1,NY,NX)=ZFEO1H(L1,NY,NX)+FXZFEO1H
      ZFEO1H(L0,NY,NX)=ZFEO1H(L0,NY,NX)-FXZFEO1H
      FXZFEO2H=FHO*ZFEO2H(L0,NY,NX)
      ZFEO2H(L1,NY,NX)=ZFEO2H(L1,NY,NX)+FXZFEO2H
      ZFEO2H(L0,NY,NX)=ZFEO2H(L0,NY,NX)-FXZFEO2H
      FXZFEO3H=FHO*ZFEO3H(L0,NY,NX)
      ZFEO3H(L1,NY,NX)=ZFEO3H(L1,NY,NX)+FXZFEO3H
      ZFEO3H(L0,NY,NX)=ZFEO3H(L0,NY,NX)-FXZFEO3H
      FXZFEO4H=FHO*ZFEO4H(L0,NY,NX)
      ZFEO4H(L1,NY,NX)=ZFEO4H(L1,NY,NX)+FXZFEO4H
      ZFEO4H(L0,NY,NX)=ZFEO4H(L0,NY,NX)-FXZFEO4H
      FXZFESH=FHO*ZFESH(L0,NY,NX)
      ZFESH(L1,NY,NX)=ZFESH(L1,NY,NX)+FXZFESH
      ZFESH(L0,NY,NX)=ZFESH(L0,NY,NX)-FXZFESH
      FXZCAOH=FHO*ZCAOH(L0,NY,NX)
      ZCAOH(L1,NY,NX)=ZCAOH(L1,NY,NX)+FXZCAOH
      ZCAOH(L0,NY,NX)=ZCAOH(L0,NY,NX)-FXZCAOH
      FXZCACH=FHO*ZCACH(L0,NY,NX)
      ZCACH(L1,NY,NX)=ZCACH(L1,NY,NX)+FXZCACH
      ZCACH(L0,NY,NX)=ZCACH(L0,NY,NX)-FXZCACH
      FXZCAHH=FHO*ZCAHH(L0,NY,NX)
      ZCAHH(L1,NY,NX)=ZCAHH(L1,NY,NX)+FXZCAHH
      ZCAHH(L0,NY,NX)=ZCAHH(L0,NY,NX)-FXZCAHH
      FXZCASH=FHO*ZCASH(L0,NY,NX)
      ZCASH(L1,NY,NX)=ZCASH(L1,NY,NX)+FXZCASH
      ZCASH(L0,NY,NX)=ZCASH(L0,NY,NX)-FXZCASH
      FXZMGOH=FHO*ZMGOH(L0,NY,NX)
      ZMGOH(L1,NY,NX)=ZMGOH(L1,NY,NX)+FXZMGOH
      ZMGOH(L0,NY,NX)=ZMGOH(L0,NY,NX)-FXZMGOH
      FXZMGCH=FHO*ZMGCH(L0,NY,NX)
      ZMGCH(L1,NY,NX)=ZMGCH(L1,NY,NX)+FXZMGCH
      ZMGCH(L0,NY,NX)=ZMGCH(L0,NY,NX)-FXZMGCH
      FXZMGHH=FHO*ZMGHH(L0,NY,NX)
      ZMGHH(L1,NY,NX)=ZMGHH(L1,NY,NX)+FXZMGHH
      ZMGHH(L0,NY,NX)=ZMGHH(L0,NY,NX)-FXZMGHH
      FXZMGSH=FHO*ZMGSH(L0,NY,NX)
      ZMGSH(L1,NY,NX)=ZMGSH(L1,NY,NX)+FXZMGSH
      ZMGSH(L0,NY,NX)=ZMGSH(L0,NY,NX)-FXZMGSH
      FXZNACH=FHO*ZNACH(L0,NY,NX)
      ZNACH(L1,NY,NX)=ZNACH(L1,NY,NX)+FXZNACH
      ZNACH(L0,NY,NX)=ZNACH(L0,NY,NX)-FXZNACH
      FXZNASH=FHO*ZNASH(L0,NY,NX)
      ZNASH(L1,NY,NX)=ZNASH(L1,NY,NX)+FXZNASH
      ZNASH(L0,NY,NX)=ZNASH(L0,NY,NX)-FXZNASH
      FXZKASH=FHO*ZKASH(L0,NY,NX)
      ZKASH(L1,NY,NX)=ZKASH(L1,NY,NX)+FXZKASH
      ZKASH(L0,NY,NX)=ZKASH(L0,NY,NX)-FXZKASH
      FXH0PO4H=FHO*H0PO4H(L0,NY,NX)
      H0PO4H(L1,NY,NX)=H0PO4H(L1,NY,NX)+FXH0PO4H
      H0PO4H(L0,NY,NX)=H0PO4H(L0,NY,NX)-FXH0PO4H
      FXH3PO4H=FHO*H3PO4H(L0,NY,NX)
      H3PO4H(L1,NY,NX)=H3PO4H(L1,NY,NX)+FXH3PO4H
      H3PO4H(L0,NY,NX)=H3PO4H(L0,NY,NX)-FXH3PO4H
      FXZFE1PH=FHO*ZFE1PH(L0,NY,NX)
      ZFE1PH(L1,NY,NX)=ZFE1PH(L1,NY,NX)+FXZFE1PH
      ZFE1PH(L0,NY,NX)=ZFE1PH(L0,NY,NX)-FXZFE1PH
      FXZFE2PH=FHO*ZFE2PH(L0,NY,NX)
      ZFE2PH(L1,NY,NX)=ZFE2PH(L1,NY,NX)+FXZFE2PH
      ZFE2PH(L0,NY,NX)=ZFE2PH(L0,NY,NX)-FXZFE2PH
      FXZCA0PH=FHO*ZCA0PH(L0,NY,NX)
      ZCA0PH(L1,NY,NX)=ZCA0PH(L1,NY,NX)+FXZCA0PH
      ZCA0PH(L0,NY,NX)=ZCA0PH(L0,NY,NX)-FXZCA0PH
      FXZCA1PH=FHO*ZCA1PH(L0,NY,NX)
      ZCA1PH(L1,NY,NX)=ZCA1PH(L1,NY,NX)+FXZCA1PH
      ZCA1PH(L0,NY,NX)=ZCA1PH(L0,NY,NX)-FXZCA1PH
      FXZCA2PH=FHO*ZCA2PH(L0,NY,NX)
      ZCA2PH(L1,NY,NX)=ZCA2PH(L1,NY,NX)+FXZCA2PH
      ZCA2PH(L0,NY,NX)=ZCA2PH(L0,NY,NX)-FXZCA2PH
      FXZMG1PH=FHO*ZMG1PH(L0,NY,NX)
      ZMG1PH(L1,NY,NX)=ZMG1PH(L1,NY,NX)+FXZMG1PH
      ZMG1PH(L0,NY,NX)=ZMG1PH(L0,NY,NX)-FXZMG1PH
      FXH0POBH=FHO*H0POBH(L0,NY,NX)
      H0POBH(L1,NY,NX)=H0POBH(L1,NY,NX)+FXH0POBH
      H0POBH(L0,NY,NX)=H0POBH(L0,NY,NX)-FXH0POBH
      FXH3POBH=FHO*H3POBH(L0,NY,NX)
      H3POBH(L1,NY,NX)=H3POBH(L1,NY,NX)+FXH3POBH
      H3POBH(L0,NY,NX)=H3POBH(L0,NY,NX)-FXH3POBH
      FXZFE1BH=FHO*ZFE1BH(L0,NY,NX)
      ZFE1BH(L1,NY,NX)=ZFE1BH(L1,NY,NX)+FXZFE1BH
      ZFE1BH(L0,NY,NX)=ZFE1BH(L0,NY,NX)-FXZFE1BH
      FXZFE2BH=FHO*ZFE2BH(L0,NY,NX)
      ZFE2BH(L1,NY,NX)=ZFE2BH(L1,NY,NX)+FXZFE2BH
      ZFE2BH(L0,NY,NX)=ZFE2BH(L0,NY,NX)-FXZFE2BH
      FXZCA0BH=FHO*ZCA0BH(L0,NY,NX)
      ZCA0BH(L1,NY,NX)=ZCA0BH(L1,NY,NX)+FXZCA0BH
      ZCA0BH(L0,NY,NX)=ZCA0BH(L0,NY,NX)-FXZCA0BH
      FXZCA1BH=FHO*ZCA1BH(L0,NY,NX)
      ZCA1BH(L1,NY,NX)=ZCA1BH(L1,NY,NX)+FXZCA1BH
      ZCA1BH(L0,NY,NX)=ZCA1BH(L0,NY,NX)-FXZCA1BH
      FXZCA2BH=FHO*ZCA2BH(L0,NY,NX)
      ZCA2BH(L1,NY,NX)=ZCA2BH(L1,NY,NX)+FXZCA2BH
      ZCA2BH(L0,NY,NX)=ZCA2BH(L0,NY,NX)-FXZCA2BH
      FXZMG1BH=FHO*ZMG1BH(L0,NY,NX)
      ZMG1BH(L1,NY,NX)=ZMG1BH(L1,NY,NX)+FXZMG1BH
      ZMG1BH(L0,NY,NX)=ZMG1BH(L0,NY,NX)-FXZMG1BH
      ENDIF
C
C     SOIL MACROPORE AQUEOUS GASES
C
      FXCO2SH=FHO*CO2SH(L0,NY,NX)
      CO2SH(L1,NY,NX)=CO2SH(L1,NY,NX)+FXCO2SH
      CO2SH(L0,NY,NX)=CO2SH(L0,NY,NX)-FXCO2SH
      FXCH4SH=FHO*CH4SH(L0,NY,NX)
      CH4SH(L1,NY,NX)=CH4SH(L1,NY,NX)+FXCH4SH
      CH4SH(L0,NY,NX)=CH4SH(L0,NY,NX)-FXCH4SH
      FXOXYSH=FHO*OXYSH(L0,NY,NX)
      OXYSH(L1,NY,NX)=OXYSH(L1,NY,NX)+FXOXYSH
      OXYSH(L0,NY,NX)=OXYSH(L0,NY,NX)-FXOXYSH
      FXZ2GSH=FHO*Z2GSH(L0,NY,NX)
      Z2GSH(L1,NY,NX)=Z2GSH(L1,NY,NX)+FXZ2GSH
      Z2GSH(L0,NY,NX)=Z2GSH(L0,NY,NX)-FXZ2GSH
      FXZ2OSH=FHO*Z2OSH(L0,NY,NX)
      Z2OSH(L1,NY,NX)=Z2OSH(L1,NY,NX)+FXZ2OSH
      Z2OSH(L0,NY,NX)=Z2OSH(L0,NY,NX)-FXZ2OSH
      FXH2GSH=FHO*H2GSH(L0,NY,NX)
      H2GSH(L1,NY,NX)=H2GSH(L1,NY,NX)+FXH2GSH
      H2GSH(L0,NY,NX)=H2GSH(L0,NY,NX)-FXH2GSH

      ENDIF
      ENDIF
C
C     SOIL ORGANIC MATTER
C
      IF(IFLGL(L,3).EQ.0.AND.L0.NE.0
     2.AND.VOLX(L0,NY,NX).GT.ZEROS(NY,NX)
     3.AND.VOLX(L1,NY,NX).GT.ZEROS(NY,NX))THEN
C     IF(L0.EQ.L.OR.CORGCI(L0,NY,NX).LE.ZERO)THEN
      FXO=FX
C     ELSE
C     FXO=AMIN1(1.0,FX*AMIN1(10.0,CORGCI(L1,NY,NX)
C    2/CORGCI(L0,NY,NX)))
C     ENDIF
      DO 7966 K=0,5
      DO 7966 N=1,7
      DO 7966 M=1,3
      FXOMC=FXO*OMC(M,N,K,L0,NY,NX)
      OMC(M,N,K,L1,NY,NX)=OMC(M,N,K,L1,NY,NX)+FXOMC
      OMC(M,N,K,L0,NY,NX)=OMC(M,N,K,L0,NY,NX)-FXOMC
      FXOMN=FXO*OMN(M,N,K,L0,NY,NX)
      OMN(M,N,K,L1,NY,NX)=OMN(M,N,K,L1,NY,NX)+FXOMN
      OMN(M,N,K,L0,NY,NX)=OMN(M,N,K,L0,NY,NX)-FXOMN
      FXOMP=FXO*OMP(M,N,K,L0,NY,NX)
      OMP(M,N,K,L1,NY,NX)=OMP(M,N,K,L1,NY,NX)+FXOMP
      OMP(M,N,K,L0,NY,NX)=OMP(M,N,K,L0,NY,NX)-FXOMP
7966  CONTINUE
      DO 7781 K=0,4
      DO 7776 M=1,2
      FXORC=FXO*ORC(M,K,L0,NY,NX)
      ORC(M,K,L1,NY,NX)=ORC(M,K,L1,NY,NX)+FXORC
      ORC(M,K,L0,NY,NX)=ORC(M,K,L0,NY,NX)-FXORC
      FXORN=FXO*ORN(M,K,L0,NY,NX)
      ORN(M,K,L1,NY,NX)=ORN(M,K,L1,NY,NX)+FXORN
      ORN(M,K,L0,NY,NX)=ORN(M,K,L0,NY,NX)-FXORN
      FXORP=FXO*ORP(M,K,L0,NY,NX)
      ORP(M,K,L1,NY,NX)=ORP(M,K,L1,NY,NX)+FXORP
      ORP(M,K,L0,NY,NX)=ORP(M,K,L0,NY,NX)-FXORP
7776  CONTINUE
      FXOQC=FXO*OQC(K,L0,NY,NX)
      OQC(K,L1,NY,NX)=OQC(K,L1,NY,NX)+FXOQC
      OQC(K,L0,NY,NX)=OQC(K,L0,NY,NX)-FXOQC
      FXOQN=FXO*OQN(K,L0,NY,NX)
      OQN(K,L1,NY,NX)=OQN(K,L1,NY,NX)+FXOQN
      OQN(K,L0,NY,NX)=OQN(K,L0,NY,NX)-FXOQN
      FXOQP=FXO*OQP(K,L0,NY,NX)
      OQP(K,L1,NY,NX)=OQP(K,L1,NY,NX)+FXOQP
      OQP(K,L0,NY,NX)=OQP(K,L0,NY,NX)-FXOQP
      FXOQA=FXO*OQA(K,L0,NY,NX)
      OQA(K,L1,NY,NX)=OQA(K,L1,NY,NX)+FXOQA
      OQA(K,L0,NY,NX)=OQA(K,L0,NY,NX)-FXOQA
      IF(FHOL(L1,NY,NX).GT.ZERO.AND.FHOL(L0,NY,NX).GT.ZERO)THEN
      FXOQCH=FXO*OQCH(K,L0,NY,NX)
      OQCH(K,L1,NY,NX)=OQCH(K,L1,NY,NX)+FXOQCH
      OQCH(K,L0,NY,NX)=OQCH(K,L0,NY,NX)-FXOQCH
      FXOQNH=FXO*OQNH(K,L0,NY,NX)
      OQNH(K,L1,NY,NX)=OQNH(K,L1,NY,NX)+FXOQNH
      OQNH(K,L0,NY,NX)=OQNH(K,L0,NY,NX)-FXOQNH
      FXOQPH=FXO*OQPH(K,L0,NY,NX)
      OQPH(K,L1,NY,NX)=OQPH(K,L1,NY,NX)+FXOQPH
      OQPH(K,L0,NY,NX)=OQPH(K,L0,NY,NX)-FXOQPH
      FXOQAH=FXO*OQAH(K,L0,NY,NX)
      OQAH(K,L1,NY,NX)=OQAH(K,L1,NY,NX)+FXOQAH
      OQAH(K,L0,NY,NX)=OQAH(K,L0,NY,NX)-FXOQAH
      ENDIF
      FXOHC=FXO*OHC(K,L0,NY,NX)
      OHC(K,L1,NY,NX)=OHC(K,L1,NY,NX)+FXOHC
      OHC(K,L0,NY,NX)=OHC(K,L0,NY,NX)-FXOHC
      FXOHN=FXO*OHN(K,L0,NY,NX)
      OHN(K,L1,NY,NX)=OHN(K,L1,NY,NX)+FXOHN
      OHN(K,L0,NY,NX)=OHN(K,L0,NY,NX)-FXOHN
      FXOHP=FXO*OHP(K,L0,NY,NX)
      OHP(K,L1,NY,NX)=OHP(K,L1,NY,NX)+FXOHP
      OHP(K,L0,NY,NX)=OHP(K,L0,NY,NX)-FXOHP
      FXOHA=FXO*OHA(K,L0,NY,NX)
      OHA(K,L1,NY,NX)=OHA(K,L1,NY,NX)+FXOHA
      OHA(K,L0,NY,NX)=OHA(K,L0,NY,NX)-FXOHA
      DO 7771 M=1,5
      FXOSC=FXO*OSC(M,K,L0,NY,NX)
      OSC(M,K,L1,NY,NX)=OSC(M,K,L1,NY,NX)+FXOSC
      OSC(M,K,L0,NY,NX)=OSC(M,K,L0,NY,NX)-FXOSC
      FXOSA=FXO*OSA(M,K,L0,NY,NX)
      OSA(M,K,L1,NY,NX)=OSA(M,K,L1,NY,NX)+FXOSA
      OSA(M,K,L0,NY,NX)=OSA(M,K,L0,NY,NX)-FXOSA
      FXOSN=FXO*OSN(M,K,L0,NY,NX)
      OSN(M,K,L1,NY,NX)=OSN(M,K,L1,NY,NX)+FXOSN
      OSN(M,K,L0,NY,NX)=OSN(M,K,L0,NY,NX)-FXOSN
      FXOSP=FXO*OSP(M,K,L0,NY,NX)
      OSP(M,K,L1,NY,NX)=OSP(M,K,L1,NY,NX)+FXOSP
      OSP(M,K,L0,NY,NX)=OSP(M,K,L0,NY,NX)-FXOSP
7771  CONTINUE
7781  CONTINUE
C
C     SOIL PLANT MATTER
C
      GO TO 6005
      DO 8901 NZ=1,NP(NY,NX)
      IF(WTRTL(1,L0,NZ,NY,NX).GT.ZEROP(NZ,NY,NX)
     2.AND.WTRTL(1,L1,NZ,NY,NX).GT.ZEROP(NZ,NY,NX))THEN
      IF(L0.EQ.L.OR.DPTHZ(L1,NY,NX).LE.ZERO)THEN
      FRO=FX
      ELSE
      FRO=AMIN1(1.0,FX*DPTHZ(L0,NY,NX)/DPTHZ(L1,NY,NX))
      ENDIF
      DO 8896 N=1,MY(NZ,NY,NX)
C
C     SOIL ROOT GASES
C
      FXCO2A=FRO*CO2A(N,L0,NZ,NY,NX)
      CO2A(N,L1,NZ,NY,NX)=CO2A(N,L1,NZ,NY,NX)+FXCO2A
      CO2A(N,L0,NZ,NY,NX)=CO2A(N,L0,NZ,NY,NX)-FXCO2A
      FXOXYA=FRO*OXYA(N,L0,NZ,NY,NX)
      OXYA(N,L1,NZ,NY,NX)=OXYA(N,L1,NZ,NY,NX)+FXOXYA
      OXYA(N,L0,NZ,NY,NX)=OXYA(N,L0,NZ,NY,NX)-FXOXYA
      FXCH4A=FRO*CH4A(N,L0,NZ,NY,NX)
      CH4A(N,L1,NZ,NY,NX)=CH4A(N,L1,NZ,NY,NX)+FXCH4A
      CH4A(N,L0,NZ,NY,NX)=CH4A(N,L0,NZ,NY,NX)-FXCH4A
      FXZ2OA=FRO*Z2OA(N,L0,NZ,NY,NX)
      Z2OA(N,L1,NZ,NY,NX)=Z2OA(N,L1,NZ,NY,NX)+FXZ2OA
      Z2OA(N,L0,NZ,NY,NX)=Z2OA(N,L0,NZ,NY,NX)-FXZ2OA
      FXZH3A=FRO*ZH3A(N,L0,NZ,NY,NX)
      ZH3A(N,L1,NZ,NY,NX)=ZH3A(N,L1,NZ,NY,NX)+FXZH3A
      ZH3A(N,L0,NZ,NY,NX)=ZH3A(N,L0,NZ,NY,NX)-FXZH3A
      FXH2GA=FRO*H2GA(N,L0,NZ,NY,NX)
      H2GA(N,L1,NZ,NY,NX)=H2GA(N,L1,NZ,NY,NX)+FXH2GA
      H2GA(N,L0,NZ,NY,NX)=H2GA(N,L0,NZ,NY,NX)-FXH2GA
      FXCO2P=FRO*CO2P(N,L0,NZ,NY,NX)
      CO2P(N,L1,NZ,NY,NX)=CO2P(N,L1,NZ,NY,NX)+FXCO2P
      CO2P(N,L0,NZ,NY,NX)=CO2P(N,L0,NZ,NY,NX)-FXCO2P
      FXOXYP=FRO*OXYP(N,L0,NZ,NY,NX)
      OXYP(N,L1,NZ,NY,NX)=OXYP(N,L1,NZ,NY,NX)+FXOXYP
      OXYP(N,L0,NZ,NY,NX)=OXYP(N,L0,NZ,NY,NX)-FXOXYP
      FXCH4P=FRO*CH4P(N,L0,NZ,NY,NX)
      CH4P(N,L1,NZ,NY,NX)=CH4P(N,L1,NZ,NY,NX)+FXCH4P
      CH4P(N,L0,NZ,NY,NX)=CH4P(N,L0,NZ,NY,NX)-FXCH4P
      FXZ2OP=FRO*Z2OP(N,L0,NZ,NY,NX)
      Z2OP(N,L1,NZ,NY,NX)=Z2OP(N,L1,NZ,NY,NX)+FXZ2OP
      Z2OP(N,L0,NZ,NY,NX)=Z2OP(N,L0,NZ,NY,NX)-FXZ2OP
      FXZH3P=FRO*ZH3P(N,L0,NZ,NY,NX)
      ZH3P(N,L1,NZ,NY,NX)=ZH3P(N,L1,NZ,NY,NX)+FXZH3P
      ZH3P(N,L0,NZ,NY,NX)=ZH3P(N,L0,NZ,NY,NX)-FXZH3P
      FXH2GP=FRO*H2GP(N,L0,NZ,NY,NX)
      H2GP(N,L1,NZ,NY,NX)=H2GP(N,L1,NZ,NY,NX)+FXH2GP
      H2GP(N,L0,NZ,NY,NX)=H2GP(N,L0,NZ,NY,NX)-FXH2GP
C
C     SOIL ROOTS
C
      DO 8871 NR=1,NRT(NZ,NY,NX)
      FXWTRT1=FRO*WTRT1(N,L0,NR,NZ,NY,NX)
      WTRT1(N,L1,NR,NZ,NY,NX)=WTRT1(N,L1,NR,NZ,NY,NX)+FXWTRT1
      WTRT1(N,L0,NR,NZ,NY,NX)=WTRT1(N,L0,NR,NZ,NY,NX)-FXWTRT1
      FXWTR1N=FRO*WTRT1N(N,L0,NR,NZ,NY,NX)
      WTRT1N(N,L1,NR,NZ,NY,NX)=WTRT1N(N,L1,NR,NZ,NY,NX)+FXWTR1N
      WTRT1N(N,L0,NR,NZ,NY,NX)=WTRT1N(N,L0,NR,NZ,NY,NX)-FXWTR1N
      FXWTR1P=FRO*WTRT1P(N,L0,NR,NZ,NY,NX)
      WTRT1P(N,L1,NR,NZ,NY,NX)=WTRT1P(N,L1,NR,NZ,NY,NX)+FXWTR1P
      WTRT1P(N,L0,NR,NZ,NY,NX)=WTRT1P(N,L0,NR,NZ,NY,NX)-FXWTR1P
      FXWTRT2=FRO*WTRT2(N,L0,NR,NZ,NY,NX)
      WTRT2(N,L1,NR,NZ,NY,NX)=WTRT2(N,L1,NR,NZ,NY,NX)+FXWTRT2
      WTRT2(N,L0,NR,NZ,NY,NX)=WTRT2(N,L0,NR,NZ,NY,NX)-FXWTRT2
      FXWTR2N=FRO*WTRT2N(N,L0,NR,NZ,NY,NX)
      WTRT2N(N,L1,NR,NZ,NY,NX)=WTRT2N(N,L1,NR,NZ,NY,NX)+FXWTR2N
      WTRT2N(N,L0,NR,NZ,NY,NX)=WTRT2N(N,L0,NR,NZ,NY,NX)-FXWTR2N
      FXWTR2P=FRO*WTRT2P(N,L0,NR,NZ,NY,NX)
      WTRT2P(N,L1,NR,NZ,NY,NX)=WTRT2P(N,L1,NR,NZ,NY,NX)+FXWTR2P
      WTRT2P(N,L0,NR,NZ,NY,NX)=WTRT2P(N,L0,NR,NZ,NY,NX)-FXWTR2P
      FXRTLG1=FRO*RTLG1(N,L0,NR,NZ,NY,NX)
      RTLG1(N,L1,NR,NZ,NY,NX)=RTLG1(N,L1,NR,NZ,NY,NX)+FXRTLG1
      RTLG1(N,L0,NR,NZ,NY,NX)=RTLG1(N,L0,NR,NZ,NY,NX)-FXRTLG1
      FXRTLG2=FRO*RTLG2(N,L0,NR,NZ,NY,NX)
      RTLG2(N,L1,NR,NZ,NY,NX)=RTLG2(N,L1,NR,NZ,NY,NX)+FXRTLG2
      RTLG2(N,L0,NR,NZ,NY,NX)=RTLG2(N,L0,NR,NZ,NY,NX)-FXRTLG2
      FXRTN2=FRO*RTN2(N,L0,NR,NZ,NY,NX)
      RTN2(N,L1,NR,NZ,NY,NX)=RTN2(N,L1,NR,NZ,NY,NX)+FXRTN2
      RTN2(N,L0,NR,NZ,NY,NX)=RTN2(N,L0,NR,NZ,NY,NX)-FXRTN2
8871  CONTINUE
      FXCPOOLR=FRO*CPOOLR(N,L0,NZ,NY,NX)
      CPOOLR(N,L1,NZ,NY,NX)=CPOOLR(N,L1,NZ,NY,NX)+FXCPOOLR
      CPOOLR(N,L0,NZ,NY,NX)=CPOOLR(N,L0,NZ,NY,NX)-FXCPOOLR
      FXZPOOLR=FRO*ZPOOLR(N,L0,NZ,NY,NX)
      ZPOOLR(N,L1,NZ,NY,NX)=ZPOOLR(N,L1,NZ,NY,NX)+FXZPOOLR
      ZPOOLR(N,L0,NZ,NY,NX)=ZPOOLR(N,L0,NZ,NY,NX)-FXZPOOLR
      FXPPOOLR=FRO*PPOOLR(N,L0,NZ,NY,NX)
      PPOOLR(N,L1,NZ,NY,NX)=PPOOLR(N,L1,NZ,NY,NX)+FXPPOOLR
      PPOOLR(N,L0,NZ,NY,NX)=PPOOLR(N,L0,NZ,NY,NX)-FXPPOOLR
      FXWTRTL=FRO*WTRTL(N,L0,NZ,NY,NX)
      WTRTL(N,L1,NZ,NY,NX)=WTRTL(N,L1,NZ,NY,NX)+FXWTRTL
      WTRTL(N,L0,NZ,NY,NX)=WTRTL(N,L0,NZ,NY,NX)-FXWTRTL
      FXWTRTD=FRO*WTRTD(N,L0,NZ,NY,NX)
      WTRTD(N,L1,NZ,NY,NX)=WTRTD(N,L1,NZ,NY,NX)+FXWTRTD
      WTRTD(N,L0,NZ,NY,NX)=WTRTD(N,L0,NZ,NY,NX)-FXWTRTD
      FXWSRTL=FRO*WSRTL(N,L0,NZ,NY,NX)
      WSRTL(N,L1,NZ,NY,NX)=WSRTL(N,L1,NZ,NY,NX)+FXWSRTL
      WSRTL(N,L0,NZ,NY,NX)=WSRTL(N,L0,NZ,NY,NX)-FXWSRTL
      FXRTN1=FRO*RTN1(N,L0,NZ,NY,NX)
      RTN1(N,L1,NZ,NY,NX)=RTN1(N,L1,NZ,NY,NX)+FXRTN1
      RTN1(N,L0,NZ,NY,NX)=RTN1(N,L0,NZ,NY,NX)-FXRTN1
      FXRTNL=FRO*RTNL(N,L0,NZ,NY,NX)
      RTNL(N,L1,NZ,NY,NX)=RTNL(N,L1,NZ,NY,NX)+FXRTNL
      RTNL(N,L0,NZ,NY,NX)=RTNL(N,L0,NZ,NY,NX)-FXRTNL
      FXRTLGP=FRO*RTLGP(N,L0,NZ,NY,NX)
      RTLGP(N,L1,NZ,NY,NX)=RTLGP(N,L1,NZ,NY,NX)+FXRTLGP
      RTLGP(N,L0,NZ,NY,NX)=RTLGP(N,L0,NZ,NY,NX)-FXRTLGP
      FXRTDNP=FRO*RTDNP(N,L0,NZ,NY,NX)
      RTDNP(N,L1,NZ,NY,NX)=RTDNP(N,L1,NZ,NY,NX)+FXRTDNP
      RTDNP(N,L0,NZ,NY,NX)=RTDNP(N,L0,NZ,NY,NX)-FXRTDNP
      FXRTVLP=FRO*RTVLP(N,L0,NZ,NY,NX)
      RTVLP(N,L1,NZ,NY,NX)=RTVLP(N,L1,NZ,NY,NX)+FXRTVLP
      RTVLP(N,L0,NZ,NY,NX)=RTVLP(N,L0,NZ,NY,NX)-FXRTVLP
      FXRTVLW=FRO*RTVLW(N,L0,NZ,NY,NX)
      RTVLW(N,L1,NZ,NY,NX)=RTVLW(N,L1,NZ,NY,NX)+FXRTVLW
      RTVLW(N,L0,NZ,NY,NX)=RTVLW(N,L0,NZ,NY,NX)-FXRTVLW
      FXRRAD1=FRO*RRAD1(N,L0,NZ,NY,NX)
      RRAD1(N,L1,NZ,NY,NX)=RRAD1(N,L1,NZ,NY,NX)+FXRRAD1
      RRAD1(N,L0,NZ,NY,NX)=RRAD1(N,L0,NZ,NY,NX)-FXRRAD1
      FXRRAD2=FRO*RRAD2(N,L0,NZ,NY,NX)
      RRAD2(N,L1,NZ,NY,NX)=RRAD2(N,L1,NZ,NY,NX)+FXRRAD2
      RRAD2(N,L0,NZ,NY,NX)=RRAD2(N,L0,NZ,NY,NX)-FXRRAD2
      FXRTARP=FRO*RTARP(N,L0,NZ,NY,NX)
      RTARP(N,L1,NZ,NY,NX)=RTARP(N,L1,NZ,NY,NX)+FXRTARP
      RTARP(N,L0,NZ,NY,NX)=RTARP(N,L0,NZ,NY,NX)-FXRTARP
      FXRTLGA=FRO*RTLGA(N,L0,NZ,NY,NX)
      RTLGA(N,L1,NZ,NY,NX)=RTLGA(N,L1,NZ,NY,NX)+FXRTLGA
      RTLGA(N,L0,NZ,NY,NX)=RTLGA(N,L0,NZ,NY,NX)-FXRTLGA
8896  CONTINUE
C
C     SOIL ROOT NODULES
C
      FXWTNDL=FRO*WTNDL(L0,NZ,NY,NX)
      WTNDL(L1,NZ,NY,NX)=WTNDL(L1,NZ,NY,NX)+FXWTNDL
      WTNDL(L0,NZ,NY,NX)=WTNDL(L0,NZ,NY,NX)-FXWTNDL
      FXWTNDLN=FRO*WTNDLN(L0,NZ,NY,NX)
      WTNDLN(L1,NZ,NY,NX)=WTNDLN(L1,NZ,NY,NX)+FXWTNDLN
      WTNDLN(L0,NZ,NY,NX)=WTNDLN(L0,NZ,NY,NX)-FXWTNDLN
      FXWTNDLP=FRO*WTNDLP(L0,NZ,NY,NX)
      WTNDLP(L1,NZ,NY,NX)=WTNDLP(L1,NZ,NY,NX)+FXWTNDLP
      WTNDLP(L0,NZ,NY,NX)=WTNDLP(L0,NZ,NY,NX)-FXWTNDLP
      FXCPOOLN=FRO*CPOOLN(L0,NZ,NY,NX)
      CPOOLN(L1,NZ,NY,NX)=CPOOLN(L1,NZ,NY,NX)+FXCPOOLN
      CPOOLN(L0,NZ,NY,NX)=CPOOLN(L0,NZ,NY,NX)-FXCPOOLN
      FXZPOOLN=FRO*ZPOOLN(L0,NZ,NY,NX)
      ZPOOLN(L1,NZ,NY,NX)=ZPOOLN(L1,NZ,NY,NX)+FXZPOOLN
      ZPOOLN(L0,NZ,NY,NX)=ZPOOLN(L0,NZ,NY,NX)-FXZPOOLN
      FXPPOOLN=FRO*PPOOLN(L0,NZ,NY,NX)
      PPOOLN(L1,NZ,NY,NX)=PPOOLN(L1,NZ,NY,NX)+FXPPOOLN
      PPOOLN(L0,NZ,NY,NX)=PPOOLN(L0,NZ,NY,NX)-FXPPOOLN
      ENDIF
8901  CONTINUE
      ENDIF
6005  CONTINUE
      IF(NN.EQ.1)THEN
      IF(BKDS(L0,NY,NX).LE.ZERO.AND.BKDS(L1,NY,NX).LE.ZERO
     3.AND.VOLW(L0,NY,NX)+VOLI(L0,NY,NX).LE.ZEROS(NY,NX))THEN
      CDPTH(L1,NY,NX)=CDPTH(L0,NY,NX)
      CDPTHY(L1,NY,NX)=CDPTHY(L0,NY,NX)
      ENDIF
      ENDIF
C     IF((I/30)*30.EQ.I.AND.J.EQ.15)THEN
C     WRITE(*,5591)'SOIL2',I,J,NFZ,NX,NY,L,L0,L1,NN,DDLYRX(NN),FX,FO
C    1,ZNH4SH(L0,NY,NX),ZNH4SH(L1,NY,NX),FXZNH4SH
C    1,BKDS(L0,NY,NX),BKDS(L1,NY,NX),FXBKDS,FXVOLW
C    5,ORGC(L0,NY,NX),ORGC(L1,NY,NX),CORGCI(L0,NY,NX),CORGCI(L1,NY,NX)
C    5,VHCM(L0,NY,NX),VHCM(L1,NY,NX),VHCP(L0,NY,NX),VHCP(L1,NY,NX)
C    5,CLAY(L0,NY,NX),CLAY(L1,NY,NX),SILT(L0,NY,NX),SILT(L1,NY,NX)
C    5,SAND(L0,NY,NX),SAND(L1,NY,NX),ORGC(L0,NY,NX),ORGC(L1,NY,NX)
C    2,VOLA(L0,NY,NX),VOLW(L0,NY,NX),VOLI(L0,NY,NX),VOLP(L0,NY,NX)
C    3,VOLA(L1,NY,NX),VOLW(L1,NY,NX),VOLI(L1,NY,NX),VOLP(L1,NY,NX)
C    3,VOLY(L0,NY,NX),VOLY(L1,NY,NX),VOLY(L0,NY,NX)+VOLY(L1,NY,NX)
C    6,DLYR(3,L0,NY,NX),DLYR(3,L1,NY,NX)
C    4,VOLA(L0,NY,NX)+VOLA(L1,NY,NX),VOLW(L0,NY,NX)+VOLW(L1,NY,NX)
C    4,VOLI(L0,NY,NX)+VOLI(L1,NY,NX),VOLP(L0,NY,NX)+VOLP(L1,NY,NX)
C    4,CDPTH(L0,NY,NX),CDPTH(L1,NY,NX)
C    5,PCAPD(L0,NY,NX),PCAPM(L1,NY,NX),FXPCAPD
C    5,((OMC(M,N,1,L0,NY,NX),M=1,3),N=1,7)
C    6,((OMC(M,N,1,L1,NY,NX),M=1,3),N=1,7)
C    5,(OSC(M,4,L0,NY,NX),M=1,2),(OSC(M,4,L1,NY,NX),M=1,2)
C    5,(FXO*OSC(M,4,L0,NY,NX),M=1,2),(FXO*OSC(M,4,L1,NY,NX),M=1,2)
C    5,(ORC(M,4,L0,NY,NX),M=1,2),(ORC(M,4,L1,NY,NX),M=1,2)
C    6,((OSC(M,4,L0,NY,NX)+OSC(M,4,L1,NY,NX)),M=1,2)
C    5,(OQC(K,L0,NY,NX),K=0,4),(OQC(K,L1,NY,NX),K=0,4)
C    6,DLYR(3,L0,NY,NX),DLYR(3,L1,NY,NX)
C    5,TKS(L0,NY,NX),TKS(L1,NY,NX),VHCP(L0,NY,NX),VHCP(L1,NY,NX)
C    6,TKS(L0,NY,NX)*VHCP(L0,NY,NX)+TKS(L1,NY,NX)*VHCP(L1,NY,NX)
C    7,VHCM(L0,NY,NX),VHCM(L1,NY,NX),ENGY0,ENGY1,FXENGY
C    6,ZNH4S(L0,NY,NX),ZNH4B(L0,NY,NX),ZNH3S(L0,NY,NX),ZNH3B(L0,NY,NX)
C    6,ZNH4S(L1,NY,NX),ZNH4B(L1,NY,NX),ZNH3S(L1,NY,NX),ZNH3B(L1,NY,NX)
C    6,OQAH(K,L,NY,NX),OQAH(K,L0,NY,NX),OQAH(K,L1,NY,NX),FXOQAH
C    7,CH4G(L0,NY,NX),CH4G(L1,NY,NX),FXCH4G
C    6,(WTRT1(1,L1,NR,2,NY,NX),NR=1,NRT(1,NY,NX))
C    6,(WTRT2(1,L1,NR,2,NY,NX),NR=1,NRT(1,NY,NX))
C    6,WTRTL(1,L1,2,NY,NX),WTNDL(L1,2,NY,NX)
C    6,CPOOLR(1,L1,2,NY,NX),ZPOOLR(1,L1,2,NY,NX)
C    6,OXYA(1,L1,2,NY,NX),OXYP(1,L1,2,NY,NX)
C    6,(WTRT1(1,L0,NR,2,NY,NX),NR=1,NRT(1,NY,NX))
C    6,(WTRT2(1,L0,NR,2,NY,NX),NR=1,NRT(1,NY,NX))
C    6,WTRTL(1,L0,2,NY,NX),WTNDL(L0,2,NY,NX)
C    6,CPOOLR(1,L0,2,NY,NX),ZPOOLR(1,L0,2,NY,NX)
C    6,OXYG(L0,NY,NX),OXYS(L0,NY,NX)
C    6,OXYG(L1,NY,NX),OXYS(L1,NY,NX)
C     ENDIF
      ENDIF
C
C     RESET LOWER LAYER NUMBER WITH EROSION
C
C     IERSNG=erosion flag from site file
C     DLYR,DLYRI=layer depth,initial value set in soil file
C     CDPTH=cumulative depth to bottom of soil layer
C     CDPTHY=CDPTH used to calculate DDLYRX excluding
C        freeze-thaw effects
C
      IF(IERSNG.EQ.1.OR.IERSNG.EQ.3)THEN
      IF(L.EQ.NL(NY,NX).AND.DLYR(3,L,NY,NX)
     2.GT.DLYRI(3,L,NY,NX))THEN
      NL(NY,NX)=MIN(NLI(NY,NX),NL(NY,NX)+1)
      ENDIF
      IF(L.EQ.NL(NY,NX)-1.AND.CDPTH(NL(NY,NX),NY,NX)
     2-CDPTH(L,NY,NX).LE.ZERO2)THEN
      CDPTH(L,NY,NX)=CDPTH(L,NY,NX)+DLYR(3,NL(NY,NX),NY,NX)
      CDPTHY(L,NY,NX)=CDPTHY(L,NY,NX)+DLYR(3,NL(NY,NX),NY,NX)
      CDPTH(NL(NY,NX),NY,NX)=CDPTH(L,NY,NX)
      CDPTHY(NL(NY,NX),NY,NX)=CDPTHY(L,NY,NX)
      DLYR(3,NL(NY,NX),NY,NX)=0.0
      NL(NY,NX)=L
C     WRITE(*,5595)'ERSNX2',I,J,NX,NY,L,NLI(NY,NX),NL(NY,NX)
C    2,DLYR(3,NL(NY,NX)+1,NY,NX),CDPTH(NL(NY,NX),NY,NX)
C    2,CDPTH(NL(NY,NX)+1,NY,NX)
5595  FORMAT(A8,7I4,12E14.6)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
230   CONTINUE
245   CONTINUE
      ENDIF
C
C     ECOSYSTEM ENERGY, CO2 FLUXES
C
C     TRN=total net SW+LW radiation absorbed by ecosystem
C     HEATI=net SW+LW radiation absorbed by snowpack+litter+soil
C        surface
C     TLE=total ecosystem latent heat flux
C     HEATE=latent heat flux at snowpack+litter+soil surface
C     TSH=total ecosystem sensible heat flux
C     SFLXC=sensible heat flux at snowpack+litter+soil surface
C     TGH=total ecosystem storage heat flux
C     HFLXC=storage heat flux at snowpack+litter+soil surface
C     TCAN=total ecosystem net CO2 flux
C     TCCAN=total PFT net CO2 fixation
C     TNBP=ecosystem net biome productivity
C     UCO2G,UCH4G=cumulative CO2,CH4 exchange by snowpack+litter+soil
C     UDOCQ,UDICQ=dissolved organic,inorganic C loss through lateral
C        and lower boundaries
C     UDOCD,UDICD=dissolved organic,inorganic C loss through
C        subsurface boundaries
C     TXCO2=CO2 net change from all solute equilibria
C
      TRN(NY,NX)=TRN(NY,NX)+HEATI(NY,NX)
      TLE(NY,NX)=TLE(NY,NX)+HEATE(NY,NX)
      TSH(NY,NX)=TSH(NY,NX)+HEATS(NY,NX)
      TGH(NY,NX)=TGH(NY,NX)-(HEATH(NY,NX)-HEATV(NY,NX))
      TRNS(NY,NX)=TRNS(NY,NX)+HEATI(NY,NX)
      TLES(NY,NX)=TLES(NY,NX)+HEATE(NY,NX)
      TSHS(NY,NX)=TSHS(NY,NX)+HEATS(NY,NX)
      TGHS(NY,NX)=TGHS(NY,NX)-(HEATH(NY,NX)-HEATV(NY,NX))
      TNPP(NY,NX)=TGPP(NY,NX)+TRAU(NY,NX)
C
C     CANOPY CO2 CONCENTRATION FROM CO2 INFLUXES AND EFFLUXES
C
C     CO2Q,CO2E=CO2 concentrations in canopy air,atmosphere
C     CH4Q,CH4E=CH4 concentrations in canopy air,atmosphere
C     RAB=aerodynamic boundary layer resistance from watsub.f
C     OXYQ,OXYE=O2 concentrations in canopy air,atmosphere
C     CNETX=net CO2 flux in canopy air from soil+plants
C     HNETX=net CH4 flux in canopy air from soil+plants
C     ONETX=net O2 flux in canopy air from soil+plants
C     OXYC=canopy air O2 content
C
      CNETX=XCNET(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      HNETX=XHNET(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      ONETX=XONET(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
      FMOLQ=1.2194E+04/TKQ(NY,NX)
      CO2Q(NY,NX)=CO2E(NY,NX)-8.333E+04*CNETX/FMOLQ*RAB(NY,NX)/XNFH
      CH4Q(NY,NX)=CH4E(NY,NX)-8.333E+04*HNETX/FMOLQ*RAB(NY,NX)/XNFH
      OXYQ(NY,NX)=OXYE(NY,NX)-3.125E+04*ONETX/FMOLQ*RAB(NY,NX)/XNFH
      OXYC(NY,NX)=OXYQ(NY,NX)*1.429E-03*273.15/TKQ(NY,NX)
     2*ZT(NY,NX)*AREA(3,NU(NY,NX),NY,NX)
      ZCNET(NY,NX)=ZCNET(NY,NX)+XCNET(NY,NX)
      ZHNET(NY,NX)=ZHNET(NY,NX)+XHNET(NY,NX)
      ZONET(NY,NX)=ZONET(NY,NX)+XONET(NY,NX)
C     IF(ICHKF.EQ.1)THEN
C     WRITE(*,6647)'CO2Q',I,J,NFZ,NX,NY
C    2,CO2Q(NY,NX),CO2E(NY,NX),CNETX
C    2,OXYQ(NY,NX),OXYE(NY,NX),ONETX
C    2,CH4Q(NY,NX),CH4E(NY,NX),HNETX
C    2,OXYC(NY,NX),ZT(NY,NX)
C    3,RAB(NY,NX),FMOLQ,TKQ(NY,NX),TKAM(NY,NX)
C    4,ROGOX(0,NY,NX),RCGOX(0,NY,NX),RCHOX(0,NY,NX),RC4OX(0,NY,NX)
C    4,XOXDFS(NY,NX),XOXFLG(3,NU(NY,NX),NY,NX),TOXYZ(NY,NX)
C    2,XOXFLG(3,0,NY,NX),XOXDFR(NY,NX)
C    3,(FLQGQ(NY,NX)+FLQRQ(NY,NX))*COXR(NY,NX)
C    4,(FLQGI(NY,NX)+FLQRI(NY,NX))*COXQ(NY,NX)
C     ENDIF
C
C     ECOSYSTEM C FLUXES
C
C     TCNET=cumulative ecosystem ne CO2 exchange
C     TCCAN,TCAN=total,cumulative net CO2 fixation
C     RECO=ecosystem respiration
C     HCO2G=ground surface CO2 exchange
C     TNBP=cumulative ecosystem net biome productivity
C     UCO2G,UCH4G=cumulative CO2,CH4 exchange by snowpack+litter+soil
C     UDOCQ,UDICQ=dissolved organic,inorganic C loss through lateral
C        and lower boundaries
C     UDOCD,UDICD=dissolved organic,inorganic C loss through
C        subsurface boundaries
C     TXCO2=CO2 net change from all solute equilibria
C
      IF(NFZ.EQ.NFH)THEN
      TCNET(NY,NX)=TCCAN(NY,NX)+HCO2G(NY,NX)
      RECO(NY,NX)=RECO(NY,NX)+HCO2G(NY,NX)
      TCAN(NY,NX)=TCAN(NY,NX)+TCCAN(NY,NX)
      TNBP(NY,NX)=TCAN(NY,NX)+UCO2G(NY,NX)+UCH4G(NY,NX)
     2-UDOCQ(NY,NX)-UDICQ(NY,NX)-UDOCD(NY,NX)-UDICD(NY,NX)
     3+TXCO2(NY,NX)
C     IF(J.EQ.15.AND.NFZ.EQ.NFH)THEN
C     WRITE(*,6647)'TNBP',I,J,NFZ,NX,NY
C    2,HEATE(NY,NX),RAC(NY,NX)
C    2,TNBP(NY,NX),TCAN(NY,NX),TCCAN(NY,NX),TCNET(NY,NX)
C    2,TNPP(NY,NX),THRE(NY,NX),UCO2G(NY,NX),UCH4G(NY,NX)
C    2,UDOCQ(NY,NX),UDICQ(NY,NX),UDOCD(NY,NX),UDICD(NY,NX)
C    3,TXCO2(NY,NX)
6647  FORMAT(A8,5I4,30E12.4)
C     ENDIF
      ENDIF
      IF(NU(NY,NX).GT.NUI(NY,NX))THEN
      DO 235 L=NUI(NY,NX),NU(NY,NX)-1
      IF(VOLX(L,NY,NX).LE.ZEROS2(NY,NX))THEN
      TKS(L,NY,NX)=TKS(NU(NY,NX),NY,NX)
      TCS(L,NY,NX)=TKS(L,NY,NX)-273.15
      ENDIF
235   CONTINUE
      ENDIF
C
C     FIRE
C
C     ICHKF=fire flag
C
      IF(ICHKF.EQ.1)THEN
C
C     RCGSK=total soil combustion from nitro.f
C     FCBOX,FCBCH=fraction of C combusted aerobically,anaerobically
C     RCGOX=O2-limited soil CO2 emission from trnsfr.f
C     RCHOX=O2,CH4-unlimited soil CH4 combustion from trnsfr.f
C     RC4OX=O2,CH4-limited soil CH4 combustion from trnsfr.f
C     FCHMC,FCHGC=fraction of anaerobic C combustion to charcoal,CO2
C     FCHMN,FCHGN=fraction of anaerobic N combustion to NH4,
C        gaseous NOX
C     FCHMP,FCHGP=fraction of anaerobic P combustion to H2PO4,
C        gaseous POX
C     FCOMC,FCOGC=fraction of aerobic C combustion to charcoal,CO2
C     FCOMN,FGOMN=fraction of aerobic N combustion to NH4,
C        gaseous NOX
C     FCOMP,FCOGP=fraction of aerobic P combustion to H2PO4,
C        gaseous POX
C     HCBFL=heat released by soil combustion
C     GCBCO,GCBCH,GCBC4=combustion energy of organic C (aerobic,
C        anaerobic), and CH4 from starts.f
C
      ZOX=0.0
      POX=0.0
      COM=0.0
      Z4M=0.0
      P4M=0.0
      DO 2850 L=0,NL(NY,NX)
      IF(RCGSK(L,NY,NX).GT.ZEROS(NY,NX))THEN
      FCBOX=AMIN1(1.0,RCGOX(L,NY,NX)/RCGSK(L,NY,NX))
      FCBCH=AMIN1(1.0,RCHOX(L,NY,NX)/(RCGSK(L,NY,NX)*FCHGC))
      HCBFL(L,NY,NX)=HCBFL(L,NY,NX)+RCGOX(L,NY,NX)*GCBCO
     2+RCHOX(L,NY,NX)/FCHGC*GCBCH+RC4OX(L,NY,NX)*GCBC4
C
C     COMBUST MICROBIAL BIOMASS
C
C     RCBOMC,RCBOMN,RCBOMP=combustion rate of OMC,OMN,OMP from nitro.f
C     CO2G,CH4G=soil gaseous CO2,CH4
C     ZNH4S,H2PO4=soil NH4,H2PO4
C     OSC(5,OSA(5=soil charcoal C,N
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C     COM,Z4M,P4M=C organic,N,P mineral additions to litter
C
      DO 2870 K=0,5
      DO 2860 N=1,7
      DO 2860 M=1,3
C     CO2G(L,NY,NX)=CO2G(L,NY,NX)+RCBOMC(M,N,K,L,NY,NX)
C    2*FCBOX*FCOGC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOMN(M,N,K,L,NY,NX)
     2*FCBOX*FCOMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOMP(M,N,K,L,NY,NX)
     2*FCBOX*FCOMP
C     CH4G(L,NY,NX)=CH4G(L,NY,NX)+RCBOMC(M,N,K,L,NY,NX)
C    2*FCBCH*FCHGC
      IF(K.LE.4)THEN
      KK=K
      ELSE
      KK=1
      ENDIF
      OSC(5,KK,L,NY,NX)=OSC(5,KK,L,NY,NX)+RCBOMC(M,N,K,L,NY,NX)
     2*FCBCH*FCHMC
      OSA(5,KK,L,NY,NX)=OSA(5,KK,L,NY,NX)+RCBOMC(M,N,K,L,NY,NX)
     2*FCBCH*FCHMC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOMN(M,N,K,L,NY,NX)
     2*FCBCH*FCHMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOMP(M,N,K,L,NY,NX)
     2*FCBCH*FCHMP
      ZOX=ZOX+RCBOMN(M,N,K,L,NY,NX)*FCBOX*FCOGN
      POX=POX+RCBOMP(M,N,K,L,NY,NX)*FCBOX*FCOGP
      ZOX=ZOX+RCBOMN(M,N,K,L,NY,NX)*FCBCH*FCHGN
      POX=POX+RCBOMP(M,N,K,L,NY,NX)*FCBCH*FCHGP
      Z4M=Z4M+RCBOMN(M,N,K,L,NY,NX)*FCBOX*FCOMN
      P4M=P4M+RCBOMP(M,N,K,L,NY,NX)*FCBOX*FCOMP
      COM=COM+RCBOMC(M,N,K,L,NY,NX)*FCBCH*FCHMC
      Z4M=Z4M+RCBOMN(M,N,K,L,NY,NX)*FCBCH*FCHMN
      P4M=P4M+RCBOMP(M,N,K,L,NY,NX)*FCBCH*FCHMP
2860  CONTINUE
2870  CONTINUE
C
C     COMBUST MICROBIAL RESIDUE
C
C     RCBORC,RCBORN,RCBORP=combustion rate of ORC,ORN,ORP from nitro.f
C     CO2G,CH4G=soil gaseous CO2,CH4
C     ZNH4S,H2PO4=soil NH4,H2PO4
C     OSC(5,OSA(5=soil charcoal C,N
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C     COM,Z4M,P4M=C organic,N,P mineral additions to litter
C
      DO 2800 K=0,4
      DO 2840 M=1,2
C     CO2G(L,NY,NX)=CO2G(L,NY,NX)+RCBORC(M,K,L,NY,NX)
C    2*FCBOX*FCOGC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBORN(M,K,L,NY,NX)
     2*FCBOX*FCOMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBORP(M,K,L,NY,NX)
     2*FCBOX*FCOMP
C     CH4G(L,NY,NX)=CH4G(L,NY,NX)+RCBORC(M,K,L,NY,NX)
C    2*FCBCH*FCHGC
      OSC(5,K,L,NY,NX)=OSC(5,K,L,NY,NX)+RCBORC(M,K,L,NY,NX)
     2*FCBCH*FCHMC
      OSA(5,K,L,NY,NX)=OSA(5,K,L,NY,NX)+RCBORC(M,K,L,NY,NX)
     2*FCBCH*FCHMC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBORN(M,K,L,NY,NX)
     2*FCBCH*FCHMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBORP(M,K,L,NY,NX)
     2*FCBCH*FCHMP
      ZOX=ZOX+RCBORN(M,K,L,NY,NX)*FCBOX*FCOGN
      POX=POX+RCBORP(M,K,L,NY,NX)*FCBOX*FCOGP
      ZOX=ZOX+RCBORN(M,K,L,NY,NX)*FCBCH*FCHGN
      POX=POX+RCBORP(M,K,L,NY,NX)*FCBCH*FCHGP
      Z4M=Z4M+RCBORN(M,K,L,NY,NX)*FCBOX*FCOMN
      P4M=P4M+RCBORP(M,K,L,NY,NX)*FCBOX*FCOMP
      COM=COM+RCBORC(M,K,L,NY,NX)*FCBCH*FCHMC
      Z4M=Z4M+RCBORN(M,K,L,NY,NX)*FCBCH*FCHMN
      P4M=P4M+RCBORP(M,K,L,NY,NX)*FCBCH*FCHMP
2840  CONTINUE
C
C     COMBUST DOC, DOA, DON, DOP
C
C     RCBOQC,RCBOQA,RCBOQN,RCBOQP=combustion rate of OQC,OQA,OQN,OQP
C        from nitro.f
C     CO2G,CH4G=soil gaseous CO2,CH4
C     ZNH4S,H2PO4=soil NH4,H2PO4
C     OSC(5,OSA(5=soil charcoal C,N
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C     COM,Z4M,P4M=C organic,N,P mineral additions to litter
C
C     CO2G(L,NY,NX)=CO2G(L,NY,NX)+(RCBOQC(K,L,NY,NX)
C    2+RCBOQA(K,L,NY,NX))*FCBOX*FCOGC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOQN(K,L,NY,NX)
     2*FCBOX*FCOMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOQP(K,L,NY,NX)
     2*FCBOX*FCOMP
C     CH4G(L,NY,NX)=CH4G(L,NY,NX)+(RCBOQC(K,L,NY,NX)
C    2+RCBOQA(K,L,NY,NX))*FCBCH*FCHGC
      OSC(5,K,L,NY,NX)=OSC(5,K,L,NY,NX)+(RCBOQC(K,L,NY,NX)
     2+RCBOQA(K,L,NY,NX))*FCBCH*FCHMC
      OSA(5,K,L,NY,NX)=OSA(5,K,L,NY,NX)+(RCBOQC(K,L,NY,NX)
     2+RCBOQA(K,L,NY,NX))*FCBCH*FCHMC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOQN(K,L,NY,NX)
     2*FCBCH*FCHMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOQP(K,L,NY,NX)
     2*FCBCH*FCHMP
      ZOX=ZOX+RCBOQN(K,L,NY,NX)*FCBOX*FCOGN
      POX=POX+RCBOQP(K,L,NY,NX)*FCBOX*FCOGP
      ZOX=ZOX+RCBOQN(K,L,NY,NX)*FCBCH*FCHGN
      POX=POX+RCBOQP(K,L,NY,NX)*FCBCH*FCHGP
      Z4M=Z4M+RCBOQN(K,L,NY,NX)*FCBOX*FCOMN
      P4M=P4M+RCBOQP(K,L,NY,NX)*FCBOX*FCOMP
      COM=COM+(RCBOQC(K,L,NY,NX)+RCBOQA(K,L,NY,NX))*FCBCH*FCHMC
      Z4M=Z4M+RCBOQN(K,L,NY,NX)*FCBCH*FCHMN
      P4M=P4M+RCBOQP(K,L,NY,NX)*FCBCH*FCHMP
C
C     COMBUST ADSORBED OM
C
C     RCBOHC,RCBOHA,RCBOHN,RCBOHP=combustion rate of OHC,OHA,OHN,OHP
C        from nitro.f
C     CO2G,CH4G=soil gaseous CO2,CH4
C     ZNH4S,H2PO4=soil NH4,H2PO4
C     OSC(5,OSA(5=soil charcoal C,N
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C     COM,Z4M,P4M=C organic,N,P mineral additions to litter
C
C     CO2G(L,NY,NX)=CO2G(L,NY,NX)+(RCBOHC(K,L,NY,NX)
C    2+RCBOHA(K,L,NY,NX))*FCBOX*FCOGC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOHN(K,L,NY,NX)
     2*FCBOX*FCOMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOHP(K,L,NY,NX)
     2*FCBOX*FCOMP
C     CH4G(L,NY,NX)=CH4G(L,NY,NX)+(RCBOHC(K,L,NY,NX)
C    2+RCBOHA(K,L,NY,NX))*FCBCH*FCHGC
      OSC(5,K,L,NY,NX)=OSC(5,K,L,NY,NX)+(RCBOHC(K,L,NY,NX)
     2+RCBOHA(K,L,NY,NX))*FCBCH*FCHMC
      OSA(5,K,L,NY,NX)=OSA(5,K,L,NY,NX)+(RCBOHC(K,L,NY,NX)
     2+RCBOHA(K,L,NY,NX))*FCBCH*FCHMC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOHN(K,L,NY,NX)
     2*FCBCH*FCHMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOHP(K,L,NY,NX)
     2*FCBCH*FCHMP
      ZOX=ZOX+RCBOHN(K,L,NY,NX)*FCBOX*FCOGN
      POX=POX+RCBOHP(K,L,NY,NX)*FCBOX*FCOGP
      ZOX=ZOX+RCBOHN(K,L,NY,NX)*FCBCH*FCHGN
      POX=POX+RCBOHP(K,L,NY,NX)*FCBCH*FCHGP
      Z4M=Z4M+RCBOHN(K,L,NY,NX)*FCBOX*FCOMN
      P4M=P4M+RCBOHP(K,L,NY,NX)*FCBOX*FCOMP
      COM=COM+(RCBOHC(K,L,NY,NX)+RCBOHA(K,L,NY,NX))*FCBCH*FCHMC
      Z4M=Z4M+RCBOHN(K,L,NY,NX)*FCBCH*FCHMN
      P4M=P4M+RCBOHP(K,L,NY,NX)*FCBCH*FCHMP
C
C     COMBUST SOM
C
C     RCBOSC,RCBOSN,RCBOSP=combustion rate of OSC,OSN,OSP from nitro.f
C     CO2G,CH4G=soil gaseous CO2,CH4
C     ZNH4S,H2PO4=soil NH4,H2PO4
C     OSC(5,OSA(5=soil charcoal C,N
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C     COM,Z4M,P4M=C organic,N,P mineral additions to litter
C
      DO 2830 M=1,4
C     CO2G(L,NY,NX)=CO2G(L,NY,NX)+RCBOSC(M,K,L,NY,NX)
C    2*FCBOX*FCOGC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOSN(M,K,L,NY,NX)
     2*FCBOX*FCOMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOSP(M,K,L,NY,NX)
     2*FCBOX*FCOMP
C     CH4G(L,NY,NX)=CH4G(L,NY,NX)+RCBOSC(M,K,L,NY,NX)
C    2*FCBCH*FCHGC
      OSC(5,K,L,NY,NX)=OSC(5,K,L,NY,NX)+RCBOSC(M,K,L,NY,NX)
     2*FCBCH*FCHMC
      OSA(5,K,L,NY,NX)=OSA(5,K,L,NY,NX)+RCBOSA(M,K,L,NY,NX)
     2*FCBCH*FCHMC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOSN(M,K,L,NY,NX)
     2*FCBCH*FCHMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOSP(M,K,L,NY,NX)
     2*FCBCH*FCHMP
      ZOX=ZOX+RCBOSN(M,K,L,NY,NX)*FCBOX*FCOGN
      POX=POX+RCBOSP(M,K,L,NY,NX)*FCBOX*FCOGP
      ZOX=ZOX+RCBOSN(M,K,L,NY,NX)*FCBCH*FCHGN
      POX=POX+RCBOSP(M,K,L,NY,NX)*FCBCH*FCHGP
      Z4M=Z4M+RCBOSN(M,K,L,NY,NX)*FCBOX*FCOMN
      P4M=P4M+RCBOSP(M,K,L,NY,NX)*FCBOX*FCOMP
      COM=COM+RCBOSC(M,K,L,NY,NX)*FCBCH*FCHMC
      Z4M=Z4M+RCBOSN(M,K,L,NY,NX)*FCBCH*FCHMN
      P4M=P4M+RCBOSP(M,K,L,NY,NX)*FCBCH*FCHMP
2830  CONTINUE
C
C     COMBUST CHAR
C
C     RCBOSC(5,RCBOSN(5,RCBOSP(5=combustion rate of OSC(5,OSN(5,OSP(5
C        from nitro.f
C     CO2G,CH4G=soil gaseous CO2,CH4
C     ZNH4S,H2PO4=soil NH4,H2PO4
C     OSC(5,OSA(5=soil charcoal C,N
C     COX,CHX,ZOX,POX=gaseous CO2,CH4,NOX,POX emissions to atmosphere
C     COM,Z4M,P4M=C organic,N,P mineral additions to litter
C
C     CO2G(L,NY,NX)=CO2G(L,NY,NX)+RCBOSC(5,K,L,NY,NX)*FCBOX
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOSN(5,K,L,NY,NX)
     2*FCBOX*FCOMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOSP(5,K,L,NY,NX)
     2*FCBOX*FCOMP
      CH4G(L,NY,NX)=CH4G(L,NY,NX)+RCBOSC(5,K,L,NY,NX)*FCBCH*FCHMC
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+RCBOSN(5,K,L,NY,NX)
     2*FCBCH*FCHMN
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+RCBOSP(5,K,L,NY,NX)
     2*FCBCH*FCHMP
      ZOX=ZOX+RCBOSN(5,K,L,NY,NX)*FCBOX*FCOGN
      POX=POX+RCBOSP(5,K,L,NY,NX)*FCBOX*FCOGP
      ZOX=ZOX+RCBOSN(5,K,L,NY,NX)*FCBCH*FCHGN
      POX=POX+RCBOSP(5,K,L,NY,NX)*FCBCH*FCHGP
      Z4M=Z4M+RCBOSN(5,K,L,NY,NX)*FCBOX*FCOMN
      P4M=P4M+RCBOSP(5,K,L,NY,NX)*FCBOX*FCOMP
      Z4M=Z4M+RCBOSN(5,K,L,NY,NX)*FCBCH*FCHMN
      P4M=P4M+RCBOSP(5,K,L,NY,NX)*FCBCH*FCHMP
      RCHOX(L,NY,NX)=RCHOX(L,NY,NX)+RCBOSC(5,K,L,NY,NX)*FCBCH*FCHMC
      TLCO2G=TLCO2G+RCBOSC(5,K,L,NY,NX)*FCBCH*FCHMC
C     WRITE(*,1188)'RCGSKL',I,J,NFZ,NX,NY,L,K,RCGSK(L,NY,NX)
C    2,RCGSK(L,NY,NX)*2.667,ROGOX(L,NY,NX),FCBOX,FCBCH,FCBOX+FCBCH
C    3,HCBFL(L,NY,NX),(RCBOSC(M,K,L,NY,NX),M=1,4),RCBOSC(5,K,L,NY,NX)
C    4,OSC(5,K,L,NY,NX),RCBOSN(5,K,L,NY,NX),OSN(5,K,L,NY,NX)
C    5,RCGOX(L,NY,NX),RCHOX(L,NY,NX),COM,ZOX,Z4M
C    6,ROGOX(L,NY,NX)/2.667+RCHOX(L,NY,NX)/FCHGC
C    4,VNOXG(NY,NX),VNOXR(NY,NX)
1188  FORMAT(A8,7I4,30F16.8)
2800  CONTINUE
C
C     FIRE BOUNDARY FLUXES
C
C     VCO2G,VCH4G=fire CO2,CH4 emission
C     RCGOX=O2-limited soil CO2 emission from trnsfr.f
C     RCHOX=O2,CH4-unlimited soil CH4 combustion from trnsfr.f
C     ZN2GIN=cumulative surface gas N2,N2O,NH3 exchange
C     TPIN=cumulative surface H2PO4,HPO4 exchange
C     TLRSDC=total landscape litter C
C     TLNH4=total landscape NH4
C     TLPO4=total landscape PO4
C     VCOXFS=fire SOC loss
C     VNOXG,VPOXG=fire NOX,POX emission
C
      VCO2G(NY,NX)=VCO2G(NY,NX)-RCGOX(L,NY,NX)
      VCH4G(NY,NX)=VCH4G(NY,NX)-RCHOX(L,NY,NX)
      ENDIF
2850  CONTINUE
      ZN2GIN=ZN2GIN-ZOX
      TPIN=TPIN-POX
      TLRSDC=TLRSDC+COM
      TLNH4=TLNH4+Z4M
      TLPO4=TLPO4+P4M
      VCOXFS(NY,NX)=VCOXFS(NY,NX)+COM
      VNOXG(NY,NX)=VNOXG(NY,NX)-ZOX
      VPOXG(NY,NX)=VPOXG(NY,NX)-POX
      ENDIF
C
C     END FIRE
C
C     OTHER DISTURBANCES
C
      IF(NFZ.EQ.1.AND.J.EQ.INT(ZNOON(NY,NX)))THEN
C
C     CHANGE EXTERNAL WATER TABLE DEPTH
C
C     ZNOON=hour of solar noon from weather file
C     ITILL=23 or 24:natural or artificial drainage from site
C        or soil management file
C     DCORP=reset external water table depth from site
C        or soil management file
C     CDPTH(NU(NY,NX)-1=soil surface elevation
C     DTBLI,DTBLDI=depth of natural,artificial water table
C        from site file
C     DTBLX,DTBLZ=current,initial natural water table depth
C     DTBLY,DTBLD=current,initial artificial water table depth
C     IDTBL=water table flag from site file
C        :0=none
C        :1,2=natural stationary,mobile
C        :3,4=artificial stationary,mobile
C     HVOLO=water loss through lateral and lower boundaries
C
C     CHANGE NATURAL DRAINAGE
C
      IF(ITILL(I,NY,NX).EQ.23)THEN
      DCORPW=DCORP(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      DTBLI(NY,NX)=DCORPW
      DTBLZ(NY,NX)=DTBLI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))
     2*(1.0-DTBLG(NY,NX))
      DTBLX(NY,NX)=DTBLZ(NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      ENDIF
C
C     CHANGE OR ADD ARTIFICIAL DRAINAGE
C
      IF(ITILL(I,NY,NX).EQ.24)THEN
      DCORPW=DCORP(I,NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      IF(IDTBL(NY,NX).EQ.1)THEN
      IDTBL(NY,NX)=3
      ELSEIF(IDTBL(NY,NX).EQ.2)THEN
      IDTBL(NY,NX)=4
      ENDIF
      DTBLDI(NY,NX)=DCORPW
      DTBLD(NY,NX)=AMAX1(0.0,DTBLDI(NY,NX)-(ALTZ(NY,NX)-ALT(NY,NX))
     2*(1.0-DTBLG(NY,NX)))
      DTBLY(NY,NX)=DTBLD(NY,NX)
      ENDIF
C
C     SET DEPTH OF MOBILE EXTERNAL WATER TABLE
C
      IF(IDTBL(NY,NX).EQ.2.OR.IDTBL(NY,NX).EQ.4)THEN
      DTBLX(NY,NX)=DTBLX(NY,NX)
     2-HVOLO(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
     2-0.00167*(DTBLX(NY,NX)-DTBLZ(NY,NX)-CDPTH(NU(NY,NX)-1,NY,NX))
      DTBLX(NY,NX)=DTBLZ(NY,NX)+CDPTH(NU(NY,NX)-1,NY,NX)
      ENDIF
      IF(IDTBL(NY,NX).EQ.4)THEN
      DTBLY(NY,NX)=DTBLY(NY,NX)
     2-HVOLO(NY,NX)/AREA(3,NU(NY,NX),NY,NX)
     3-0.00167*(DTBLY(NY,NX)-DTBLD(NY,NX))
      ENDIF
C
C     IF LITTER REMOVAL EVENT IS ENTERED IN DISTURBANCE FILE
C
C     ITILL=21:litter removal
C     DCORP=fraction of surface litter removed
C
      IF(ITILL(I,NY,NX).EQ.21)THEN
      DCORPC=AMIN1(0.999,DCORP(I,NY,NX))
      DC=0.0
      DN=0.0
      DP=0.0
      OC=0.0
      ON=0.0
      OP=0.0
      DCC=0.0
      DNC=0.0
      DPC=0.0
      OCC=0.0
      ONC=0.0
      OPC=0.0
      DO 2970 K=0,2
C
C     REMOVE LITTER MICROBIAL BIOMASS
C
      DO 2960 N=1,7
      DO 2960 M=1,3
      OCH=DCORPC*OMC(M,N,K,L,NY,NX)
      ONH=DCORPC*OMN(M,N,K,L,NY,NX)
      OPH=DCORPC*OMP(M,N,K,L,NY,NX)
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)-OCH
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)-ONH
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)-OPH
      DC=DC+OMC(M,N,K,L,NY,NX)
      DN=DN+OMN(M,N,K,L,NY,NX)
      DP=DP+OMP(M,N,K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONH
      OP=OP+OPH
2960  CONTINUE
2970  CONTINUE
C
C     REMOVE LITTER MICROBIAL RESIDUE
C
      DO 2900 K=0,2
      DO 2940 M=1,2
      OCH=DCORPC*ORC(M,K,L,NY,NX)
      ONH=DCORPC*ORN(M,K,L,NY,NX)
      OPH=DCORPC*ORP(M,K,L,NY,NX)
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)-OCH
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)-ONH
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)-OPH
      DC=DC+ORC(M,K,L,NY,NX)
      DN=DN+ORN(M,K,L,NY,NX)
      DP=DP+ORP(M,K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONH
      OP=OP+OPH
2940  CONTINUE
C
C     REMOVE LITTER DISSOLVED DOC, DOA, DON, DOP
C
      OCH=DCORPC*OQC(K,L,NY,NX)
      OAH=DCORPC*OQA(K,L,NY,NX)
      ONH=DCORPC*OQN(K,L,NY,NX)
      OPH=DCORPC*OQP(K,L,NY,NX)
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)-OCH
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)-OAH
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)-ONH
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)-OPH
      OC=OC+OCH+OAH
      ON=ON+ONH
      OP=OP+OPH
      OCH=DCORPC*OQCH(K,L,NY,NX)
      OAH=DCORPC*OQAH(K,L,NY,NX)
      ONH=DCORPC*OQNH(K,L,NY,NX)
      OPH=DCORPC*OQPH(K,L,NY,NX)
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)-OCH
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)-OAH
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)-ONH
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)-OPH
      OC=OC+OCH+OAH
      ON=ON+ONH
      OP=OP+OPH
C
C     REMOVE LITTER ADSORBED OHC, OHA, OHN, OHP
C
      OCH=DCORPC*OHC(K,L,NY,NX)
      OAH=DCORPC*OHA(K,L,NY,NX)
      ONH=DCORPC*OHN(K,L,NY,NX)
      OPH=DCORPC*OHP(K,L,NY,NX)
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)-OCH
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)-OAH
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)-ONH
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)-OPH
      DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
     2+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      DP=DP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
      OC=OC+OCH
      ON=ON+ONH
      OP=OP+OPH
C
C     REMOVE LITTER SOC, colonized SOC, SON, SOP
C
      DO 2930 M=1,5
      OCH=DCORPC*OSC(M,K,L,NY,NX)
      OCA=DCORPC*OSA(M,K,L,NY,NX)
      ONH=DCORPC*OSN(M,K,L,NY,NX)
      OPH=DCORPC*OSP(M,K,L,NY,NX)
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)-OCH
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)-OCA
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)-ONH
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)-OPH
      IF(M.LE.4)THEN
      DC=DC+OSC(M,K,L,NY,NX)
      DN=DN+OSN(M,K,L,NY,NX)
      DP=DP+OSP(M,K,L,NY,NX)
      ELSE
      DCC=DCC+OSC(M,K,L,NY,NX)
      DNC=DNC+OSN(M,K,L,NY,NX)
      DPC=DPC+OSP(M,K,L,NY,NX)
      ENDIF
      OC=OC+OCH
      ON=ON+ONH
      OP=OP+OPH
2930  CONTINUE
2900  CONTINUE
C
C     REMOVE FERTILIZER IN LITTER
C
      ON=ON+DCORPC*(ZNH4S(0,NY,NX)+ZNH3S(0,NY,NX)
     2+ZNO3S(0,NY,NX)+ZNO2S(0,NY,NX))
      OP=OP+DCORPC*(H1PO4(0,NY,NX)+H2PO4(0,NY,NX))
      ZNH4S(0,NY,NX)=(1.0-DCORPC)*ZNH4S(0,NY,NX)
      ZNH3S(0,NY,NX)=(1.0-DCORPC)*ZNH3S(0,NY,NX)
      ZNO3S(0,NY,NX)=(1.0-DCORPC)*ZNO3S(0,NY,NX)
      ZNO2S(0,NY,NX)=(1.0-DCORPC)*ZNO2S(0,NY,NX)
      H1PO4(0,NY,NX)=(1.0-DCORPC)*H1PO4(0,NY,NX)
      H2PO4(0,NY,NX)=(1.0-DCORPC)*H2PO4(0,NY,NX)
      XN4(0,NY,NX)=(1.0-DCORPC)*XN4(0,NY,NX)
      PALPO(0,NY,NX)=(1.0-DCORPC)*PALPO(0,NY,NX)
      PFEPO(0,NY,NX)=(1.0-DCORPC)*PFEPO(0,NY,NX)
      PCAPD(0,NY,NX)=(1.0-DCORPC)*PCAPD(0,NY,NX)
      PCAPH(0,NY,NX)=(1.0-DCORPC)*PCAPH(0,NY,NX)
      PCAPM(0,NY,NX)=(1.0-DCORPC)*PCAPM(0,NY,NX)
      ZNH4FA(0,NY,NX)=(1.0-DCORPC)*ZNH4FA(0,NY,NX)
      ZNH3FA(0,NY,NX)=(1.0-DCORPC)*ZNH3FA(0,NY,NX)
      ZNHUFA(0,NY,NX)=(1.0-DCORPC)*ZNHUFA(0,NY,NX)
      ZNO3FA(0,NY,NX)=(1.0-DCORPC)*ZNO3FA(0,NY,NX)
      ORGC(0,NY,NX)=DC
      ORGN(0,NY,NX)=DN
      ORGCC(0,NY,NX)=DCC
      ORGNC(0,NY,NX)=DNC
      HFLXD=2.496E-06*(ORGCX(0,NY,NX)-ORGC(0,NY,NX)
     2-ORGCC(0,NY,NX))*TKS(0,NY,NX)
      HEATOU=HEATOU+HFLXD
C     VHCP(0,NY,NX)=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))
C    2+4.19*(VOLW(0,NY,NX)+VOLV(0,NY,NX))
C    3+1.9274*VOLI(0,NY,NX)
      TCOU=TCOU+OC
      TZOU=TZOU+ON
      TPOU=TPOU+OP
      UDOCQ(NY,NX)=UDOCQ(NY,NX)+OC
      UDONQ(NY,NX)=UDONQ(NY,NX)+ON
      UDOPQ(NY,NX)=UDOPQ(NY,NX)+OP
      TNBP(NY,NX)=TNBP(NY,NX)-OC
      ENDIF
C
C     TILLAGE: MIX ALL SOIL MATERIAL AND INCORPORATE ALL SURFACE
C     RESIDUE STATE VARIABLES WITHIN THE TILLAGE ZONE TO THE EXTENT
C     ASSOCIATED IN 'DAY' WITH EACH TILLAGE EVENT ENTERED IN THE
C     TILLAGE FILE
C
C     IFLGS=disturbance flag
C     ITILL=1 - 20:tillage options from soil management file
C     XCORP=soil mixing fraction from day.f
C     DCORP=tillage depth from soil management file
C
      IF(ITILL(I,NY,NX).GE.1.AND.ITILL(I,NY,NX).LE.20)THEN
C
C     EXTENT OF MIXING
C
      IFLGS(NY,NX)=1
      CORP=1.0-XCORP(NY,NX)
      ENGYP(NY,NX)=0.0
C
C     TEMPORARY ACCUMULATORS FOR CALCULATING SOIL MATERIAL MIXING
C
      TBKDX=0.0
      TFC=0.0
      TWP=0.0
      TSCNV=0.0
      TSCNH=0.0
      TSAND=0.0
      TSILT=0.0
      TCLAY=0.0
      TXCEC=0.0
      TXAEC=0.0
      TGKC4=0.0
      TGKCA=0.0
      TGKCM=0.0
      TGKCN=0.0
      TGKCK=0.0
      TVOLW=0.0
      TVOLV=0.0
      TVOLI=0.0
C     TVOLP=0.0
C     TVOLA=0.0
      TENGY=0.0
      TVHCM=0.0
      TNFNIH=0.0
      TNH4FA=0.0
      TNH3FA=0.0
      TNHUFA=0.0
      TNO3FA=0.0
      TNH4FB=0.0
      TNH3FB=0.0
      TNHUFB=0.0
      TNO3FB=0.0
      TNH4S=0.0
      TNH4B=0.0
      TNH3S=0.0
      TNH3B=0.0
      TNO3S=0.0
      TNO3B=0.0
      TNO2S=0.0
      TNO2B=0.0
      TZAL=0.0
      TZFE=0.0
      TZHY=0.0
      TZCA=0.0
      TZMG=0.0
      TZNA=0.0
      TZKA=0.0
      TZOH=0.0
      TZSO4=0.0
      TZCL=0.0
      TZCO3=0.0
      TZHCO3=0.0
      TZALO1=0.0
      TZALO2=0.0
      TZALO3=0.0
      TZALO4=0.0
      TZALS=0.0
      TZFEO1=0.0
      TZFEO2=0.0
      TZFEO3=0.0
      TZFEO4=0.0
      TZFES=0.0
      TZCAO=0.0
      TZCAC=0.0
      TZCAH=0.0
      TZCAS=0.0
      TZMGO=0.0
      TZMGC=0.0
      TZMGH=0.0
      TZMGS=0.0
      TZNAC=0.0
      TZNAS=0.0
      TZKAS=0.0
      TH0PO4=0.0
      TH1PO4=0.0
      TH2PO4=0.0
      TH3PO4=0.0
      TZFE1P=0.0
      TZFE2P=0.0
      TZCA0P=0.0
      TZCA1P=0.0
      TZCA2P=0.0
      TZMG1P=0.0
      TH0POB=0.0
      TH1POB=0.0
      TH2POB=0.0
      TH3POB=0.0
      TFE1PB=0.0
      TFE2PB=0.0
      TCA0PB=0.0
      TCA1PB=0.0
      TCA2PB=0.0
      TMG1PB=0.0
      TXNH4=0.0
      TXNHB=0.0
      TXHY=0.0
      TXAL=0.0
      TXFE=0.0
      TXCA=0.0
      TXMG=0.0
      TXNA=0.0
      TXKA=0.0
      TXHC=0.0
      TXAL2=0.0
      TXFE2=0.0
      TXOH0=0.0
      TXOH1=0.0
      TXOH2=0.0
      TXH1P=0.0
      TXH2P=0.0
      TXOH0B=0.0
      TXOH1B=0.0
      TXOH2B=0.0
      TXH1PB=0.0
      TXH2PB=0.0
      TPALOH=0.0
      TPFEOH=0.0
      TPCACO=0.0
      TPCASO=0.0
      TPALPO=0.0
      TPFEPO=0.0
      TPCAPD=0.0
      TPCAPH=0.0
      TPCAPM=0.0
      TPALPB=0.0
      TPFEPB=0.0
      TPCPDB=0.0
      TPCPHB=0.0
      TPCPMB=0.0
      TCO2G=0.0
      TCH4G=0.0
      TCOZS=0.0
      TCHFS=0.0
      TOXYG=0.0
      TOXYS=0.0
      TZ2GG=0.0
      TZ2GS=0.0
      TZ2OG=0.0
      TZ2OS=0.0
      TZNH3G=0.0
      TH2GG=0.0
      TH2GS=0.0
      DO 3990 K=0,5
      DO 3990 N=1,7
      DO 3990 M=1,3
      TOMC(M,N,K)=0.0
      TOMN(M,N,K)=0.0
      TOMP(M,N,K)=0.0
3990  CONTINUE
      DO 3980 K=0,4
      DO 3975 M=1,2
      TORC(M,K)=0.0
      TORN(M,K)=0.0
      TORP(M,K)=0.0
3975  CONTINUE
      TOQC(K)=0.0
      TOQN(K)=0.0
      TOQP(K)=0.0
      TOQA(K)=0.0
      TOHC(K)=0.0
      TOHN(K)=0.0
      TOHP(K)=0.0
      TOHA(K)=0.0
      DO 3970 M=1,5
      TOSC(M,K)=0.0
      TOSA(M,K)=0.0
      TOSN(M,K)=0.0
      TOSP(M,K)=0.0
3970  CONTINUE
3980  CONTINUE
      TZNFN2=0.0
      TZNFNI=0.0
      ZNHUX0=0.0
      ZNHUXI=0.0
      ZNFNX0=0.0
C
C     ACCUMULATE STATE VARIABLES IN SURFACE LITTER FOR ADDITION
C     TO SOIL TO TILLAGE MIXING DEPTH
C
      XCORP0=AMAX1(0.001,XCORP(NY,NX))
      CORP0=1.0-XCORP0
      DC=0.0
      DN=0.0
      DP=0.0
      DCC=0.0
      DNC=0.0
      DPC=0.0
      DO 3950 K=0,5
      IF(K.NE.3.AND.K.NE.4)THEN
      DO 3945 N=1,7
      DO 3945 M=1,3
      TOMGC(M,N,K)=OMC(M,N,K,0,NY,NX)*CORP0
      TOMGN(M,N,K)=OMN(M,N,K,0,NY,NX)*CORP0
      TOMGP(M,N,K)=OMP(M,N,K,0,NY,NX)*CORP0
      OMC(M,N,K,0,NY,NX)=OMC(M,N,K,0,NY,NX)*XCORP0
      OMN(M,N,K,0,NY,NX)=OMN(M,N,K,0,NY,NX)*XCORP0
      OMP(M,N,K,0,NY,NX)=OMP(M,N,K,0,NY,NX)*XCORP0
      DC=DC+OMC(M,N,K,0,NY,NX)
      DN=DN+OMN(M,N,K,0,NY,NX)
      DP=DP+OMP(M,N,K,0,NY,NX)
3945  CONTINUE
      ENDIF
3950  CONTINUE
      DO 3940 K=0,2
      DO 3935 M=1,2
      TORXC(M,K)=ORC(M,K,0,NY,NX)*CORP0
      TORXN(M,K)=ORN(M,K,0,NY,NX)*CORP0
      TORXP(M,K)=ORP(M,K,0,NY,NX)*CORP0
      ORC(M,K,0,NY,NX)=ORC(M,K,0,NY,NX)*XCORP0
      ORN(M,K,0,NY,NX)=ORN(M,K,0,NY,NX)*XCORP0
      ORP(M,K,0,NY,NX)=ORP(M,K,0,NY,NX)*XCORP0
      DC=DC+ORC(M,K,0,NY,NX)
      DN=DN+ORN(M,K,0,NY,NX)
      DP=DP+ORP(M,K,0,NY,NX)
3935  CONTINUE
      TOQGC(K)=OQC(K,0,NY,NX)*CORP0
      TOQGN(K)=OQN(K,0,NY,NX)*CORP0
      TOQGP(K)=OQP(K,0,NY,NX)*CORP0
      TOQGA(K)=OQA(K,0,NY,NX)*CORP0
      TOQHC(K)=OQCH(K,0,NY,NX)*CORP0
      TOQHN(K)=OQNH(K,0,NY,NX)*CORP0
      TOQHP(K)=OQPH(K,0,NY,NX)*CORP0
      TOQHA(K)=OQAH(K,0,NY,NX)*CORP0
      TOHGC(K)=OHC(K,0,NY,NX)*CORP0
      TOHGN(K)=OHN(K,0,NY,NX)*CORP0
      TOHGP(K)=OHP(K,0,NY,NX)*CORP0
      TOHGA(K)=OHA(K,0,NY,NX)*CORP0
C
C     REDUCE SURFACE RESIDUE STATE VARIABLES FOR INCORPORATION
C
      OQC(K,0,NY,NX)=OQC(K,0,NY,NX)*XCORP0
      OQN(K,0,NY,NX)=OQN(K,0,NY,NX)*XCORP0
      OQP(K,0,NY,NX)=OQP(K,0,NY,NX)*XCORP0
      OQA(K,0,NY,NX)=OQA(K,0,NY,NX)*XCORP0
      OQCH(K,0,NY,NX)=OQCH(K,0,NY,NX)*XCORP0
      OQNH(K,0,NY,NX)=OQNH(K,0,NY,NX)*XCORP0
      OQPH(K,0,NY,NX)=OQPH(K,0,NY,NX)*XCORP0
      OQAH(K,0,NY,NX)=OQAH(K,0,NY,NX)*XCORP0
      OHC(K,0,NY,NX)=OHC(K,0,NY,NX)*XCORP0
      OHN(K,0,NY,NX)=OHN(K,0,NY,NX)*XCORP0
      OHP(K,0,NY,NX)=OHP(K,0,NY,NX)*XCORP0
      OHA(K,0,NY,NX)=OHA(K,0,NY,NX)*XCORP0
      DC=DC+OQC(K,0,NY,NX)+OQCH(K,0,NY,NX)+OHC(K,0,NY,NX)
     2+OQA(K,0,NY,NX)+OQAH(K,0,NY,NX)+OHA(K,0,NY,NX)
      DN=DN+OQN(K,0,NY,NX)+OQNH(K,0,NY,NX)+OHN(K,0,NY,NX)
      DP=DP+OQP(K,0,NY,NX)+OQPH(K,0,NY,NX)+OHP(K,0,NY,NX)
      DO 3965 M=1,5
      TOSGC(M,K)=OSC(M,K,0,NY,NX)*CORP0
      TOSGA(M,K)=OSA(M,K,0,NY,NX)*CORP0
      TOSGN(M,K)=OSN(M,K,0,NY,NX)*CORP0
      TOSGP(M,K)=OSP(M,K,0,NY,NX)*CORP0
      OSC(M,K,0,NY,NX)=OSC(M,K,0,NY,NX)*XCORP0
      OSA(M,K,0,NY,NX)=OSA(M,K,0,NY,NX)*XCORP0
      OSN(M,K,0,NY,NX)=OSN(M,K,0,NY,NX)*XCORP0
      OSP(M,K,0,NY,NX)=OSP(M,K,0,NY,NX)*XCORP0
      IF(M.LE.4)THEN
      DC=DC+OSC(M,K,0,NY,NX)
      DN=DN+OSN(M,K,0,NY,NX)
      DP=DP+OSP(M,K,0,NY,NX)
      ELSE
      DCC=DCC+OSC(M,K,0,NY,NX)
      DNC=DNC+OSN(M,K,0,NY,NX)
      DPC=DPC+OSP(M,K,0,NY,NX)
      ENDIF
3965  CONTINUE
3940  CONTINUE
      TCO2GS=CO2S(0,NY,NX)*CORP0
      TCH4GS=CH4S(0,NY,NX)*CORP0
      TOXYGS=OXYS(0,NY,NX)*CORP0
      TZ2GSG=Z2GS(0,NY,NX)*CORP0
      TZ2OGS=Z2OS(0,NY,NX)*CORP0
      TH2GGS=H2GS(0,NY,NX)*CORP0
      TNH4GS=ZNH4S(0,NY,NX)*CORP0
      TNH3GS=ZNH3S(0,NY,NX)*CORP0
      TNO3GS=ZNO3S(0,NY,NX)*CORP0
      TNO2GS=ZNO2S(0,NY,NX)*CORP0
      TP14GS=H1PO4(0,NY,NX)*CORP0
      TPO4GS=H2PO4(0,NY,NX)*CORP0
      TXN4G=XN4(0,NY,NX)*CORP0
      TXOH0G=XOH0(0,NY,NX)*CORP0
      TXOH1G=XOH1(0,NY,NX)*CORP0
      TXOH2G=XOH2(0,NY,NX)*CORP0
      TXH1PG=XH1P(0,NY,NX)*CORP0
      TXH2PG=XH2P(0,NY,NX)*CORP0
      TALPOG=PALPO(0,NY,NX)*CORP0
      TFEPOG=PFEPO(0,NY,NX)*CORP0
      TCAPDG=PCAPD(0,NY,NX)*CORP0
      TCAPHG=PCAPH(0,NY,NX)*CORP0
      TCAPMG=PCAPM(0,NY,NX)*CORP0
      TNH4FG=ZNH4FA(0,NY,NX)*CORP0
      TNH3FG=ZNH3FA(0,NY,NX)*CORP0
      TNHUFG=ZNHUFA(0,NY,NX)*CORP0
      TNO3FG=ZNO3FA(0,NY,NX)*CORP0
      TZNFNG=ZNFNI(0,NY,NX)*CORP0
      TVOLWR=VOLW(0,NY,NX)*CORP0
      HFLXD=2.496E-06*ORGC(0,NY,NX)*CORP0*TKS(0,NY,NX)
      HEATIN=HEATIN-HFLXD
      HEATSO=HEATSO-HFLXD
      TENGYR=4.19*TVOLWR*TKS(0,NY,NX)
      ORGC(0,NY,NX)=DC
      ORGN(0,NY,NX)=DN
      ORGCC(0,NY,NX)=DCC
      ORGNC(0,NY,NX)=DNC
      ORGR(0,NY,NX)=DC
      CO2S(0,NY,NX)=CO2S(0,NY,NX)*XCORP0
      CH4S(0,NY,NX)=CH4S(0,NY,NX)*XCORP0
      OXYS(0,NY,NX)=OXYS(0,NY,NX)*XCORP0
      Z2GS(0,NY,NX)=Z2GS(0,NY,NX)*XCORP0
      Z2OS(0,NY,NX)=Z2OS(0,NY,NX)*XCORP0
      H2GS(0,NY,NX)=H2GS(0,NY,NX)*XCORP0
      ZNH4S(0,NY,NX)=ZNH4S(0,NY,NX)*XCORP0
      ZNH3S(0,NY,NX)=ZNH3S(0,NY,NX)*XCORP0
      ZNO3S(0,NY,NX)=ZNO3S(0,NY,NX)*XCORP0
      ZNO2S(0,NY,NX)=ZNO2S(0,NY,NX)*XCORP0
      H1PO4(0,NY,NX)=H1PO4(0,NY,NX)*XCORP0
      H2PO4(0,NY,NX)=H2PO4(0,NY,NX)*XCORP0
      XN4(0,NY,NX)=XN4(0,NY,NX)*XCORP0
      XOH0(0,NY,NX)=XOH0(0,NY,NX)*XCORP0
      XOH1(0,NY,NX)=XOH1(0,NY,NX)*XCORP0
      XOH2(0,NY,NX)=XOH2(0,NY,NX)*XCORP0
      XH1P(0,NY,NX)=XH1P(0,NY,NX)*XCORP0
      XH2P(0,NY,NX)=XH2P(0,NY,NX)*XCORP0
      PALPO(0,NY,NX)=PALPO(0,NY,NX)*XCORP0
      PFEPO(0,NY,NX)=PFEPO(0,NY,NX)*XCORP0
      PCAPD(0,NY,NX)=PCAPD(0,NY,NX)*XCORP0
      PCAPH(0,NY,NX)=PCAPH(0,NY,NX)*XCORP0
      PCAPM(0,NY,NX)=PCAPM(0,NY,NX)*XCORP0
      ZNH4FA(0,NY,NX)=ZNH4FA(0,NY,NX)*XCORP0
      ZNH3FA(0,NY,NX)=ZNH3FA(0,NY,NX)*XCORP0
      ZNHUFA(0,NY,NX)=ZNHUFA(0,NY,NX)*XCORP0
      ZNO3FA(0,NY,NX)=ZNO3FA(0,NY,NX)*XCORP0
      VOLW(0,NY,NX)=VOLW(0,NY,NX)*XCORP0
      VOLV(0,NY,NX)=VOLV(0,NY,NX)*XCORP0
      VHCP(0,NY,NX)=2.496E-06*(ORGC(0,NY,NX)+ORGCC(0,NY,NX))
     2+4.19*(VOLW(0,NY,NX)+VOLV(0,NY,NX))
     2+1.9274*VOLI(0,NY,NX)
      VOLR(NY,NX)=VOLR(NY,NX)*XCORP0
      VOLT(0,NY,NX)=VOLT(0,NY,NX)*XCORP0
      ZNHUX0=AMAX1(ZNHUX0,ZNHU0(0,NY,NX))
      ZNHUXI=AMAX1(ZNHUXI,ZNHUI(0,NY,NX))
      ZNFNX0=AMAX1(ZNFNX0,ZNFN0(0,NY,NX))
      LL=NU(NY,NX)
C
C     REDISTRIBUTE SOIL MATERIAL DURING TILLAGE
C
      DCORPZ=AMIN1(DCORP(I,NY,NX),CDPTHZ(NL(NY,NX),NY,NX))
C
C     ACCUMULATE STATE VARIABLES FOR EACH SOIL LAYER WITHIN MIXING
C     ZONE
C
C     CDPTHZ=cumulative depth to bottom of soil layer from current
C        soil surface depth
C     DLYR=layer depth
C
      DO 1000 L=NU(NY,NX),NL(NY,NX)
      IF(CDPTHZ(L,NY,NX)-DLYR(3,L,NY,NX).LT.DCORPZ
     2.AND.DLYR(3,L,NY,NX).GT.ZERO)THEN
C
C     TL=soil layer depth within tillage mixing zone
C     FI=fraction of TL within tillage mixing zone
C     TI=fraction of TL within soil layer
C
      TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CDPTHZ(L,NY,NX)
     2-DLYR(3,L,NY,NX)))
      FI=TL/DCORPZ
      TI=TL/DLYR(3,L,NY,NX)
C
C     ACCUMULATE SOIL MATERIAL IN EACH SOIL LATER OF TILLAGE MIXING
C     ZONE
C
C     SEDIMENT, WATER, ICE, HEAT, GAPON COEFFICIENTS
C
      TBKDX=TBKDX+FI*BKDSI(L,NY,NX)
      TFC=TFC+FI*FC(L,NY,NX)
      TWP=TWP+FI*WP(L,NY,NX)
      TSCNV=TSCNV+FI*SCNV(L,NY,NX)
      TSCNH=TSCNH+FI*SCNH(L,NY,NX)
      TSAND=TSAND+TI*SAND(L,NY,NX)
      TSILT=TSILT+TI*SILT(L,NY,NX)
      TCLAY=TCLAY+TI*CLAY(L,NY,NX)
      TXCEC=TXCEC+TI*XCEC(L,NY,NX)
      TXAEC=TXAEC+TI*XAEC(L,NY,NX)
      TGKC4=TGKC4+FI*GKC4(L,NY,NX)
      TGKCA=TGKCA+FI*GKCA(L,NY,NX)
      TGKCM=TGKCM+FI*GKCM(L,NY,NX)
      TGKCN=TGKCN+FI*GKCN(L,NY,NX)
      TGKCK=TGKCK+FI*GKCK(L,NY,NX)
      TVOLW=TVOLW+TI*VOLW(L,NY,NX)
      TVOLV=TVOLV+TI*VOLV(L,NY,NX)
      TVOLI=TVOLI+TI*VOLI(L,NY,NX)
C     TVOLP=TVOLP+TI*VOLP(L,NY,NX)
C     TVOLA=TVOLA+TI*VOLA(L,NY,NX)
      TENGY=TENGY+TI*(4.19*(VOLW(L,NY,NX)+VOLWH(L,NY,NX))
     2+1.9274*(VOLI(L,NY,NX)+VOLIH(L,NY,NX)))*TKS(L,NY,NX)
C
C     FERTILIZER IN BAND, NON-BAND
C
      TNH4FA=TNH4FA+TI*ZNH4FA(L,NY,NX)
      TNH3FA=TNH3FA+TI*ZNH3FA(L,NY,NX)
      TNHUFA=TNHUFA+TI*ZNHUFA(L,NY,NX)
      TNO3FA=TNO3FA+TI*ZNO3FA(L,NY,NX)
      TNH4FB=TNH4FB+TI*ZNH4FB(L,NY,NX)
      TNH3FB=TNH3FB+TI*ZNH3FB(L,NY,NX)
      TNHUFB=TNHUFB+TI*ZNHUFB(L,NY,NX)
      TNO3FB=TNO3FB+TI*ZNO3FB(L,NY,NX)
C
C     C,N SOLUTES IN BAND, NON-BAND
C
      TNH4S=TNH4S+TI*ZNH4S(L,NY,NX)
      TNH4B=TNH4B+TI*ZNH4B(L,NY,NX)
      TNH3S=TNH3S+TI*ZNH3S(L,NY,NX)
      TNH3B=TNH3B+TI*ZNH3B(L,NY,NX)
      TNO3S=TNO3S+TI*ZNO3S(L,NY,NX)
      TNO3B=TNO3B+TI*ZNO3B(L,NY,NX)
      TNO2S=TNO2S+TI*ZNO2S(L,NY,NX)
      TNO2B=TNO2B+TI*ZNO2B(L,NY,NX)
C
C     SALTS IN NON-BAND
C
      TZAL=TZAL+TI*ZAL(L,NY,NX)
      TZFE=TZFE+TI*ZFE(L,NY,NX)
      TZHY=TZHY+TI*ZHY(L,NY,NX)
      TZCA=TZCA+TI*ZCA(L,NY,NX)
      TZMG=TZMG+TI*ZMG(L,NY,NX)
      TZNA=TZNA+TI*ZNA(L,NY,NX)
      TZKA=TZKA+TI*ZKA(L,NY,NX)
      TZOH=TZOH+TI*ZOH(L,NY,NX)
      TZSO4=TZSO4+TI*ZSO4(L,NY,NX)
      TZCL=TZCL+TI*ZCL(L,NY,NX)
      TZCO3=TZCO3+TI*ZCO3(L,NY,NX)
      TZHCO3=TZHCO3+TI*ZHCO3(L,NY,NX)
      TZALO1=TZALO1+TI*ZALOH1(L,NY,NX)
      TZALO2=TZALO2+TI*ZALOH2(L,NY,NX)
      TZALO3=TZALO3+TI*ZALOH3(L,NY,NX)
      TZALO4=TZALO4+TI*ZALOH4(L,NY,NX)
      TZALS=TZALS+TI*ZALS(L,NY,NX)
      TZFEO1=TZFEO1+TI*ZFEOH1(L,NY,NX)
      TZFEO2=TZFEO2+TI*ZFEOH2(L,NY,NX)
      TZFEO3=TZFEO3+TI*ZFEOH3(L,NY,NX)
      TZFEO4=TZFEO4+TI*ZFEOH4(L,NY,NX)
      TZFES=TZFES+TI*ZFES(L,NY,NX)
      TZCAO=TZCAO+TI*ZCAO(L,NY,NX)
      TZCAC=TZCAC+TI*ZCAC(L,NY,NX)
      TZCAH=TZCAH+TI*ZCAH(L,NY,NX)
      TZCAS=TZCAS+TI*ZCAS(L,NY,NX)
      TZMGO=TZMGO+TI*ZMGO(L,NY,NX)
      TZMGC=TZMGC+TI*ZMGC(L,NY,NX)
      TZMGH=TZMGH+TI*ZMGH(L,NY,NX)
      TZMGS=TZMGS+TI*ZMGS(L,NY,NX)
      TZNAC=TZNAC+TI*ZNAC(L,NY,NX)
      TZNAS=TZNAS+TI*ZNAS(L,NY,NX)
      TZKAS=TZKAS+TI*ZKAS(L,NY,NX)
      TH0PO4=TH0PO4+TI*H0PO4(L,NY,NX)
      TH1PO4=TH1PO4+TI*H1PO4(L,NY,NX)
      TH2PO4=TH2PO4+TI*H2PO4(L,NY,NX)
      TH3PO4=TH3PO4+TI*H3PO4(L,NY,NX)
      TZFE1P=TZFE1P+TI*ZFE1P(L,NY,NX)
      TZFE2P=TZFE2P+TI*ZFE2P(L,NY,NX)
      TZCA0P=TZCA0P+TI*ZCA0P(L,NY,NX)
      TZCA1P=TZCA1P+TI*ZCA1P(L,NY,NX)
      TZCA2P=TZCA2P+TI*ZCA2P(L,NY,NX)
      TZMG1P=TZMG1P+TI*ZMG1P(L,NY,NX)
C
C     SALTS IN BAND
C
      TH0POB=TH0POB+TI*H0POB(L,NY,NX)
      TH1POB=TH1POB+TI*H1POB(L,NY,NX)
      TH2POB=TH2POB+TI*H2POB(L,NY,NX)
      TH3POB=TH3POB+TI*H3POB(L,NY,NX)
      TFE1PB=TFE1PB+TI*ZFE1PB(L,NY,NX)
      TFE2PB=TFE2PB+TI*ZFE2PB(L,NY,NX)
      TCA0PB=TCA0PB+TI*ZCA0PB(L,NY,NX)
      TCA1PB=TCA1PB+TI*ZCA1PB(L,NY,NX)
      TCA2PB=TCA2PB+TI*ZCA2PB(L,NY,NX)
      TMG1PB=TMG1PB+TI*ZMG1PB(L,NY,NX)
C
C     ADSORBED CATIONS AND ANIONS IN BAND, NON-BAND
C
      TXNH4=TXNH4+TI*XN4(L,NY,NX)
      TXNHB=TXNHB+TI*XNB(L,NY,NX)
      TXHY=TXHY+TI*XHY(L,NY,NX)
      TXAL=TXAL+TI*XAL(L,NY,NX)
      TXFE=TXFE+TI*XFE(L,NY,NX)
      TXCA=TXCA+TI*XCA(L,NY,NX)
      TXMG=TXMG+TI*XMG(L,NY,NX)
      TXNA=TXNA+TI*XNA(L,NY,NX)
      TXKA=TXKA+TI*XKA(L,NY,NX)
      TXHC=TXHC+TI*XHC(L,NY,NX)
      TXAL2=TXAL2+TI*XALO2(L,NY,NX)
      TXFE2=TXFE2+TI*XFEO2(L,NY,NX)
      TXOH0=TXOH0+TI*XOH0(L,NY,NX)
      TXOH1=TXOH1+TI*XOH1(L,NY,NX)
      TXOH2=TXOH2+TI*XOH2(L,NY,NX)
      TXH1P=TXH1P+TI*XH1P(L,NY,NX)
      TXH2P=TXH2P+TI*XH2P(L,NY,NX)
      TXOH0B=TXOH0B+TI*XOH0B(L,NY,NX)
      TXOH1B=TXOH1B+TI*XOH1B(L,NY,NX)
      TXOH2B=TXOH2B+TI*XOH2B(L,NY,NX)
      TXH1PB=TXH1PB+TI*XH1PB(L,NY,NX)
      TXH2PB=TXH2PB+TI*XH2PB(L,NY,NX)
C
C     PRECIPITATED CATIONS AND ANIONS IN BAND, NON-BAND
C
      TPALOH=TPALOH+TI*PALOH(L,NY,NX)
      TPFEOH=TPFEOH+TI*PFEOH(L,NY,NX)
      TPCACO=TPCACO+TI*PCACO(L,NY,NX)
      TPCASO=TPCASO+TI*PCASO(L,NY,NX)
      TPALPO=TPALPO+TI*PALPO(L,NY,NX)
      TPFEPO=TPFEPO+TI*PFEPO(L,NY,NX)
      TPCAPD=TPCAPD+TI*PCAPD(L,NY,NX)
      TPCAPH=TPCAPH+TI*PCAPH(L,NY,NX)
      TPCAPM=TPCAPM+TI*PCAPM(L,NY,NX)
      TPALPB=TPALPB+TI*PALPB(L,NY,NX)
      TPFEPB=TPFEPB+TI*PFEPB(L,NY,NX)
      TPCPDB=TPCPDB+TI*PCPDB(L,NY,NX)
      TPCPHB=TPCPHB+TI*PCPHB(L,NY,NX)
      TPCPMB=TPCPMB+TI*PCPMB(L,NY,NX)
C
C     GASEOUS AND AQUEOUS GASES
C
      TCO2G=TCO2G+TI*CO2G(L,NY,NX)
      TCH4G=TCH4G+TI*CH4G(L,NY,NX)
      TCOZS=TCOZS+TI*CO2S(L,NY,NX)
      TCHFS=TCHFS+TI*CH4S(L,NY,NX)
      TOXYG=TOXYG+TI*OXYG(L,NY,NX)
      TOXYS=TOXYS+TI*OXYS(L,NY,NX)
      TZ2GG=TZ2GG+TI*Z2GG(L,NY,NX)
      TZ2GS=TZ2GS+TI*Z2GS(L,NY,NX)
      TZ2OG=TZ2OG+TI*Z2OG(L,NY,NX)
      TZ2OS=TZ2OS+TI*Z2OS(L,NY,NX)
      TZNH3G=TZNH3G+TI*ZNH3G(L,NY,NX)
      TH2GG=TH2GG+TI*H2GG(L,NY,NX)
      TH2GS=TH2GS+TI*H2GS(L,NY,NX)
C
C     ORGANIC MATTER
C
      DO 4985 K=0,5
      DO 4985 N=1,7
      DO 4985 M=1,3
      TOMC(M,N,K)=TOMC(M,N,K)+TI*OMC(M,N,K,L,NY,NX)
      TOMN(M,N,K)=TOMN(M,N,K)+TI*OMN(M,N,K,L,NY,NX)
      TOMP(M,N,K)=TOMP(M,N,K)+TI*OMP(M,N,K,L,NY,NX)
4985  CONTINUE
      DO 4980 K=0,4
      DO 4975 M=1,2
      TORC(M,K)=TORC(M,K)+TI*ORC(M,K,L,NY,NX)
      TORN(M,K)=TORN(M,K)+TI*ORN(M,K,L,NY,NX)
      TORP(M,K)=TORP(M,K)+TI*ORP(M,K,L,NY,NX)
4975  CONTINUE
      TOQC(K)=TOQC(K)+TI*OQC(K,L,NY,NX)
      TOQN(K)=TOQN(K)+TI*OQN(K,L,NY,NX)
      TOQP(K)=TOQP(K)+TI*OQP(K,L,NY,NX)
      TOQA(K)=TOQA(K)+TI*OQA(K,L,NY,NX)
      TOHC(K)=TOHC(K)+TI*OHC(K,L,NY,NX)
      TOHN(K)=TOHN(K)+TI*OHN(K,L,NY,NX)
      TOHP(K)=TOHP(K)+TI*OHP(K,L,NY,NX)
      TOHA(K)=TOHA(K)+TI*OHA(K,L,NY,NX)
      DO 4970 M=1,5
      TOSC(M,K)=TOSC(M,K)+TI*OSC(M,K,L,NY,NX)
      TOSA(M,K)=TOSA(M,K)+TI*OSA(M,K,L,NY,NX)
      TOSN(M,K)=TOSN(M,K)+TI*OSN(M,K,L,NY,NX)
      TOSP(M,K)=TOSP(M,K)+TI*OSP(M,K,L,NY,NX)
4970  CONTINUE
4980  CONTINUE
      ZNHUX0=AMAX1(ZNHUX0,ZNHU0(L,NY,NX))
      ZNHUXI=AMAX1(ZNHUXI,ZNHUI(L,NY,NX))
      ZNFNX0=AMAX1(ZNFNX0,ZNFN0(L,NY,NX))
      TZNFNI=TZNFNI+ZNFNI(L,NY,NX)
      LL=L
      ENDIF
1000  CONTINUE
C
C     CHANGE SOIL STATE VARIABLES IN TILLAGE MIXING ZONE
C     TO ACCOUNT FOR REDISTRIBUTION FROM MIXING
C
C     DLYR=layer depth
C     TL=soil layer depth within tillage mixing zone
C     FI=fraction of TL within tillage mixing zone
C     TI=fraction of TL within soil layer
C     CORP=1.0-mixing fraction
C
      DO 2000 L=NU(NY,NX),LL
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
      TL=AMIN1(DLYR(3,L,NY,NX),DCORPZ-(CDPTHZ(L,NY,NX)
     2-DLYR(3,L,NY,NX)))
      FI=TL/DCORPZ
      TI=TL/DLYR(3,L,NY,NX)
      TX=1.0-TI
C
C     SEDIMENT, WATER, ICE, HEAT, GAPON COEFFICIENTS
C
      BKDSI(L,NY,NX)=TI*(BKDSI(L,NY,NX)+CORP*(TBKDX-BKDSI(L,NY,NX)))
     2+TX*BKDSI(L,NY,NX)
      FC(L,NY,NX)=TI*(FC(L,NY,NX)+CORP*(TFC-FC(L,NY,NX)))
     2+TX*FC(L,NY,NX)
      WP(L,NY,NX)=TI*(WP(L,NY,NX)+CORP*(TWP-WP(L,NY,NX)))
     2+TX*WP(L,NY,NX)
      SCNV(L,NY,NX)=TI*(SCNV(L,NY,NX)+CORP*(TSCNV-SCNV(L,NY,NX)))
     2+TX*SCNV(L,NY,NX)
      SCNH(L,NY,NX)=TI*(SCNH(L,NY,NX)+CORP*(TSCNH-SCNH(L,NY,NX)))
     2+TX*SCNH(L,NY,NX)
      SAND(L,NY,NX)=TI*SAND(L,NY,NX)+CORP*(FI*TSAND-TI*SAND(L,NY,NX))
     2+TX*SAND(L,NY,NX)
      SILT(L,NY,NX)=TI*SILT(L,NY,NX)+CORP*(FI*TSILT-TI*SILT(L,NY,NX))
     2+TX*SILT(L,NY,NX)
      CLAY(L,NY,NX)=TI*CLAY(L,NY,NX)+CORP*(FI*TCLAY-TI*CLAY(L,NY,NX))
     2+TX*CLAY(L,NY,NX)
      XCEC(L,NY,NX)=TI*XCEC(L,NY,NX)+CORP*(FI*TXCEC-TI*XCEC(L,NY,NX))
     2+TX*XCEC(L,NY,NX)
      XAEC(L,NY,NX)=TI*XAEC(L,NY,NX)+CORP*(FI*TXAEC-TI*XAEC(L,NY,NX))
     2+TX*XAEC(L,NY,NX)
      GKC4(L,NY,NX)=TI*(GKC4(L,NY,NX)+CORP*(TGKC4-GKC4(L,NY,NX)))
     2+TX*GKC4(L,NY,NX)
      GKCA(L,NY,NX)=TI*(GKCA(L,NY,NX)+CORP*(TGKCA-GKCA(L,NY,NX)))
     2+TX*GKCA(L,NY,NX)
      GKCM(L,NY,NX)=TI*(GKCM(L,NY,NX)+CORP*(TGKCM-GKCM(L,NY,NX)))
     2+TX*GKCM(L,NY,NX)
      GKCN(L,NY,NX)=TI*(GKCN(L,NY,NX)+CORP*(TGKCN-GKCN(L,NY,NX)))
     2+TX*GKCN(L,NY,NX)
      GKCK(L,NY,NX)=TI*(GKCK(L,NY,NX)+CORP*(TGKCK-GKCK(L,NY,NX)))
     2+TX*GKCK(L,NY,NX)
      ENGYM=VHCM(L,NY,NX)*TKS(L,NY,NX)
      ENGYV=(4.19*(VOLW(L,NY,NX)+VOLWH(L,NY,NX))
     2+1.9274*(VOLI(L,NY,NX)+VOLIH(L,NY,NX)))*TKS(L,NY,NX)
      VOLW(L,NY,NX)=TI*VOLW(L,NY,NX)+CORP*(FI*TVOLW-TI*VOLW(L,NY,NX))
     2+TX*VOLW(L,NY,NX)+FI*TVOLWR
      VOLV(L,NY,NX)=TI*VOLV(L,NY,NX)+CORP*(FI*TVOLV-TI*VOLV(L,NY,NX))
     2+TX*VOLV(L,NY,NX)
      VOLI(L,NY,NX)=TI*VOLI(L,NY,NX)+CORP*(FI*TVOLI-TI*VOLI(L,NY,NX))
     2+TX*VOLI(L,NY,NX)
C     VOLP(L,NY,NX)=TI*VOLP(L,NY,NX)+CORP*(FI*TVOLP-TI*VOLP(L,NY,NX))
C    2+TX*VOLP(L,NY,NX)
C     VOLA(L,NY,NX)=TI*VOLA(L,NY,NX)+CORP*(FI*TVOLA-TI*VOLA(L,NY,NX))
C    2+TX*VOLA(L,NY,NX)
      VOLWX(L,NY,NX)=VOLW(L,NY,NX)
C     VOLW(L,NY,NX)=VOLW(L,NY,NX)+CORP*VOLWH(L,NY,NX)
C     VOLI(L,NY,NX)=VOLI(L,NY,NX)+CORP*VOLIH(L,NY,NX)
C     VOLA(L,NY,NX)=VOLA(L,NY,NX)+CORP*VOLAH(L,NY,NX)
C     VOLWH(L,NY,NX)=XCORP(NY,NX)*VOLWH(L,NY,NX)
C     VOLIH(L,NY,NX)=XCORP(NY,NX)*VOLIH(L,NY,NX)
C     VOLAH(L,NY,NX)=XCORP(NY,NX)*VOLAH(L,NY,NX)
C     FHOL(L,NY,NX)=XCORP(NY,NX)*FHOL(L,NY,NX)
      ENGYL=TI*ENGYV+CORP*(FI*TENGY-TI*ENGYV)+TX*ENGYV+FI*TENGYR
      VHCP(L,NY,NX)=VHCM(L,NY,NX)
     2+4.19*(VOLW(L,NY,NX)+VOLV(L,NY,NX)+VOLWH(L,NY,NX))
     2+1.9274*(VOLI(L,NY,NX)+VOLIH(L,NY,NX))
      TKS(L,NY,NX)=(ENGYM+ENGYL)/VHCP(L,NY,NX)
      TCS(L,NY,NX)=TKS(L,NY,NX)-273.15
C
C     FERTILIZER IN BAND, NON-BAND
C
      ZNH4FA(L,NY,NX)=TI*ZNH4FA(L,NY,NX)+CORP*(FI*TNH4FA
     2-TI*ZNH4FA(L,NY,NX))+TX*ZNH4FA(L,NY,NX)
      ZNH3FA(L,NY,NX)=TI*ZNH3FA(L,NY,NX)+CORP*(FI*TNH3FA
     2-TI*ZNH3FA(L,NY,NX))+TX*ZNH3FA(L,NY,NX)
      ZNHUFA(L,NY,NX)=TI*ZNHUFA(L,NY,NX)+CORP*(FI*TNHUFA
     2-TI*ZNHUFA(L,NY,NX))+TX*ZNHUFA(L,NY,NX)
      ZNO3FA(L,NY,NX)=TI*ZNO3FA(L,NY,NX)+CORP*(FI*TNO3FA
     2-TI*ZNO3FA(L,NY,NX))+TX*ZNO3FA(L,NY,NX)
      ZNH4FB(L,NY,NX)=TI*ZNH4FB(L,NY,NX)+CORP*(FI*TNH4FB
     2-TI*ZNH4FB(L,NY,NX))+TX*ZNH4FB(L,NY,NX)
      ZNH3FB(L,NY,NX)=TI*ZNH3FB(L,NY,NX)+CORP*(FI*TNH3FB
     2-TI*ZNH3FB(L,NY,NX))+TX*ZNH3FB(L,NY,NX)
      ZNHUFB(L,NY,NX)=TI*ZNHUFB(L,NY,NX)+CORP*(FI*TNHUFB
     2-TI*ZNHUFB(L,NY,NX))+TX*ZNHUFB(L,NY,NX)
      ZNO3FB(L,NY,NX)=TI*ZNO3FB(L,NY,NX)+CORP*(FI*TNO3FB
     2-TI*ZNO3FB(L,NY,NX))+TX*ZNO3FB(L,NY,NX)
C
C     C,N SOLUTES IN BAND, NON-BAND
C
      ZNH4S(L,NY,NX)=TI*ZNH4S(L,NY,NX)+CORP*(FI*TNH4S-TI*ZNH4S(L,NY,NX))
     2+TX*ZNH4S(L,NY,NX)+CORP*ZNH4SH(L,NY,NX)
      ZNH4B(L,NY,NX)=TI*ZNH4B(L,NY,NX)+CORP*(FI*TNH4B-TI*ZNH4B(L,NY,NX))
     2+TX*ZNH4B(L,NY,NX)+CORP*ZNH4BH(L,NY,NX)
      ZNH3S(L,NY,NX)=TI*ZNH3S(L,NY,NX)+CORP*(FI*TNH3S-TI*ZNH3S(L,NY,NX))
     2+TX*ZNH3S(L,NY,NX)+CORP*ZNH3SH(L,NY,NX)
      ZNH3B(L,NY,NX)=TI*ZNH3B(L,NY,NX)+CORP*(FI*TNH3B-TI*ZNH3B(L,NY,NX))
     2+TX*ZNH3B(L,NY,NX)+CORP*ZNH3BH(L,NY,NX)
      ZNO3S(L,NY,NX)=TI*ZNO3S(L,NY,NX)+CORP*(FI*TNO3S-TI*ZNO3S(L,NY,NX))
     2+TX*ZNO3S(L,NY,NX)+CORP*ZNO3SH(L,NY,NX)
      ZNO3B(L,NY,NX)=TI*ZNO3B(L,NY,NX)+CORP*(FI*TNO3B-TI*ZNO3B(L,NY,NX))
     2+TX*ZNO3B(L,NY,NX)+CORP*ZNO3BH(L,NY,NX)
      ZNO2S(L,NY,NX)=TI*ZNO2S(L,NY,NX)+CORP*(FI*TNO2S-TI*ZNO2S(L,NY,NX))
     2+TX*ZNO2S(L,NY,NX)+CORP*ZNO2SH(L,NY,NX)
      ZNO2B(L,NY,NX)=TI*ZNO2B(L,NY,NX)+CORP*(FI*TNO2B-TI*ZNO2B(L,NY,NX))
     2+TX*ZNO2B(L,NY,NX)+CORP*ZNO2BH(L,NY,NX)
C
C     SALTS IN NON-BAND
C
      ZAL(L,NY,NX)=TI*ZAL(L,NY,NX)+CORP*(FI*TZAL-TI*ZAL(L,NY,NX))
     2+TX*ZAL(L,NY,NX)+CORP*ZALH(L,NY,NX)
      ZFE(L,NY,NX)=TI*ZFE(L,NY,NX)+CORP*(FI*TZFE-TI*ZFE(L,NY,NX))
     2+TX*ZFE(L,NY,NX)+CORP*ZFEH(L,NY,NX)
      ZHY(L,NY,NX)=TI*ZHY(L,NY,NX)+CORP*(FI*TZHY-TI*ZHY(L,NY,NX))
     2+TX*ZHY(L,NY,NX)+CORP*ZHYH(L,NY,NX)
      ZCA(L,NY,NX)=TI*ZCA(L,NY,NX)+CORP*(FI*TZCA-TI*ZCA(L,NY,NX))
     2+TX*ZCA(L,NY,NX)+CORP*ZCCH(L,NY,NX)
      ZMG(L,NY,NX)=TI*ZMG(L,NY,NX)+CORP*(FI*TZMG-TI*ZMG(L,NY,NX))
     2+TX*ZMG(L,NY,NX)+CORP*ZMAH(L,NY,NX)
      ZNA(L,NY,NX)=TI*ZNA(L,NY,NX)+CORP*(FI*TZNA-TI*ZNA(L,NY,NX))
     2+TX*ZNA(L,NY,NX)+CORP*ZNAH(L,NY,NX)
      ZKA(L,NY,NX)=TI*ZKA(L,NY,NX)+CORP*(FI*TZKA-TI*ZKA(L,NY,NX))
     2+TX*ZKA(L,NY,NX)+CORP*ZKAH(L,NY,NX)
      ZOH(L,NY,NX)=TI*ZOH(L,NY,NX)+CORP*(FI*TZOH-TI*ZOH(L,NY,NX))
     2+TX*ZOH(L,NY,NX)+CORP*ZOHH(L,NY,NX)
      ZSO4(L,NY,NX)=TI*ZSO4(L,NY,NX)+CORP*(FI*TZSO4-TI*ZSO4(L,NY,NX))
     2+TX*ZSO4(L,NY,NX)+CORP*ZSO4H(L,NY,NX)
      ZCL(L,NY,NX)=TI*ZCL(L,NY,NX)+CORP*(FI*TZCL-TI*ZCL(L,NY,NX))
     2+TX*ZCL(L,NY,NX)+CORP*ZCLH(L,NY,NX)
      ZCO3(L,NY,NX)=TI*ZCO3(L,NY,NX)+CORP*(FI*TZCO3-TI*ZCO3(L,NY,NX))
     2+TX*ZCO3(L,NY,NX)+CORP*ZCO3H(L,NY,NX)
      ZHCO3(L,NY,NX)=TI*ZHCO3(L,NY,NX)+CORP*(FI*TZHCO3
     2-TI*ZHCO3(L,NY,NX))+TX*ZHCO3(L,NY,NX)+CORP*ZHCO3H(L,NY,NX)
      ZALOH1(L,NY,NX)=TI*ZALOH1(L,NY,NX)+CORP*(FI*TZALO1
     2-TI*ZALOH1(L,NY,NX))+TX*ZALOH1(L,NY,NX)+CORP*ZALO1H(L,NY,NX)
      ZALOH2(L,NY,NX)=TI*ZALOH2(L,NY,NX)+CORP*(FI*TZALO2
     2-TI*ZALOH2(L,NY,NX))+TX*ZALOH2(L,NY,NX)+CORP*ZALO2H(L,NY,NX)
      ZALOH3(L,NY,NX)=TI*ZALOH3(L,NY,NX)+CORP*(FI*TZALO3
     2-TI*ZALOH3(L,NY,NX))+TX*ZALOH3(L,NY,NX)+CORP*ZALO3H(L,NY,NX)
      ZALOH4(L,NY,NX)=TI*ZALOH4(L,NY,NX)+CORP*(FI*TZALO4
     2-TI*ZALOH4(L,NY,NX))+TX*ZALOH4(L,NY,NX)+CORP*ZALO4H(L,NY,NX)
      ZALS(L,NY,NX)=TI*ZALS(L,NY,NX)+CORP*(FI*TZALS-TI*ZALS(L,NY,NX))
     2+TX*ZALS(L,NY,NX)+CORP*ZALSH(L,NY,NX)
      ZFEOH1(L,NY,NX)=TI*ZFEOH1(L,NY,NX)+CORP*(FI*TZFEO1
     2-TI*ZFEOH1(L,NY,NX))+TX*ZFEOH1(L,NY,NX)+CORP*ZFEO1H(L,NY,NX)
      ZFEOH2(L,NY,NX)=TI*ZFEOH2(L,NY,NX)+CORP*(FI*TZFEO2
     2-TI*ZFEOH2(L,NY,NX))+TX*ZFEOH2(L,NY,NX)+CORP*ZFEO2H(L,NY,NX)
      ZFEOH3(L,NY,NX)=TI*ZFEOH3(L,NY,NX)+CORP*(FI*TZFEO3
     2-TI*ZFEOH3(L,NY,NX))+TX*ZFEOH3(L,NY,NX)+CORP*ZFEO3H(L,NY,NX)
      ZFEOH4(L,NY,NX)=TI*ZFEOH4(L,NY,NX)+CORP*(FI*TZFEO4
     2-TI*ZFEOH4(L,NY,NX))+TX*ZFEOH4(L,NY,NX)+CORP*ZFEO4H(L,NY,NX)
      ZFES(L,NY,NX)=TI*ZFES(L,NY,NX)+CORP*(FI*TZFES-TI*ZFES(L,NY,NX))
     2+TX*ZFES(L,NY,NX)+CORP*ZFESH(L,NY,NX)
      ZCAO(L,NY,NX)=TI*ZCAO(L,NY,NX)+CORP*(FI*TZCAO-TI*ZCAO(L,NY,NX))
     2+TX*ZCAO(L,NY,NX)+CORP*ZCAOH(L,NY,NX)
      ZCAC(L,NY,NX)=TI*ZCAC(L,NY,NX)+CORP*(FI*TZCAC-TI*ZCAC(L,NY,NX))
     2+TX*ZCAC(L,NY,NX)+CORP*ZCACH(L,NY,NX)
      ZCAH(L,NY,NX)=TI*ZCAH(L,NY,NX)+CORP*(FI*TZCAH-TI*ZCAH(L,NY,NX))
     2+TX*ZCAH(L,NY,NX)+CORP*ZCAHH(L,NY,NX)
      ZCAS(L,NY,NX)=TI*ZCAS(L,NY,NX)+CORP*(FI*TZCAS-TI*ZCAS(L,NY,NX))
     2+TX*ZCAS(L,NY,NX)+CORP*ZCASH(L,NY,NX)
      ZMGO(L,NY,NX)=TI*ZMGO(L,NY,NX)+CORP*(FI*TZMGO-TI*ZMGO(L,NY,NX))
     2+TX*ZMGO(L,NY,NX)+CORP*ZMGOH(L,NY,NX)
      ZMGC(L,NY,NX)=TI*ZMGC(L,NY,NX)+CORP*(FI*TZMGC-TI*ZMGC(L,NY,NX))
     2+TX*ZMGC(L,NY,NX)+CORP*ZMGCH(L,NY,NX)
      ZMGH(L,NY,NX)=TI*ZMGH(L,NY,NX)+CORP*(FI*TZMGH-TI*ZMGH(L,NY,NX))
     2+TX*ZMGH(L,NY,NX)+CORP*ZMGHH(L,NY,NX)
      ZMGS(L,NY,NX)=TI*ZMGS(L,NY,NX)+CORP*(FI*TZMGS-TI*ZMGS(L,NY,NX))
     2+TX*ZMGS(L,NY,NX)+CORP*ZMGSH(L,NY,NX)
      ZNAC(L,NY,NX)=TI*ZNAC(L,NY,NX)+CORP*(FI*TZNAC-TI*ZNAC(L,NY,NX))
     2+TX*ZNAC(L,NY,NX)+CORP*ZNACH(L,NY,NX)
      ZNAS(L,NY,NX)=TI*ZNAS(L,NY,NX)+CORP*(FI*TZNAS-TI*ZNAS(L,NY,NX))
     2+TX*ZNAS(L,NY,NX)+CORP*ZNASH(L,NY,NX)
      ZKAS(L,NY,NX)=TI*ZKAS(L,NY,NX)+CORP*(FI*TZKAS-TI*ZKAS(L,NY,NX))
     2+TX*ZKAS(L,NY,NX)+CORP*ZKASH(L,NY,NX)
      H0PO4(L,NY,NX)=TI*H0PO4(L,NY,NX)+CORP*(FI*TH0PO4
     2-TI*H0PO4(L,NY,NX))+TX*H0PO4(L,NY,NX)+CORP*H0PO4H(L,NY,NX)
      H1PO4(L,NY,NX)=TI*H1PO4(L,NY,NX)+CORP*(FI*TH1PO4
     2-TI*H1PO4(L,NY,NX))+TX*H1PO4(L,NY,NX)+CORP*H1PO4H(L,NY,NX)
      H2PO4(L,NY,NX)=TI*H2PO4(L,NY,NX)+CORP*(FI*TH2PO4
     2-TI*H2PO4(L,NY,NX))+TX*H2PO4(L,NY,NX)+CORP*H2PO4H(L,NY,NX)
      H3PO4(L,NY,NX)=TI*H3PO4(L,NY,NX)+CORP*(FI*TH3PO4
     2-TI*H3PO4(L,NY,NX))+TX*H3PO4(L,NY,NX)+CORP*H3PO4H(L,NY,NX)
      ZFE1P(L,NY,NX)=TI*ZFE1P(L,NY,NX)+CORP*(FI*TZFE1P
     2-TI*ZFE1P(L,NY,NX))+TX*ZFE1P(L,NY,NX)+CORP*ZFE1PH(L,NY,NX)
      ZFE2P(L,NY,NX)=TI*ZFE2P(L,NY,NX)+CORP*(FI*TZFE2P
     2-TI*ZFE2P(L,NY,NX))+TX*ZFE2P(L,NY,NX)+CORP*ZFE2PH(L,NY,NX)
      ZCA0P(L,NY,NX)=TI*ZCA0P(L,NY,NX)+CORP*(FI*TZCA0P
     2-TI*ZCA0P(L,NY,NX))+TX*ZCA0P(L,NY,NX)+CORP*ZCA0PH(L,NY,NX)
      ZCA1P(L,NY,NX)=TI*ZCA1P(L,NY,NX)+CORP*(FI*TZCA1P
     2-TI*ZCA1P(L,NY,NX))+TX*ZCA1P(L,NY,NX)+CORP*ZCA1PH(L,NY,NX)
      ZCA2P(L,NY,NX)=TI*ZCA2P(L,NY,NX)+CORP*(FI*TZCA2P
     2-TI*ZCA2P(L,NY,NX))+TX*ZCA2P(L,NY,NX)+CORP*ZCA2PH(L,NY,NX)
      ZMG1P(L,NY,NX)=TI*ZMG1P(L,NY,NX)+CORP*(FI*TZMG1P
     2-TI*ZMG1P(L,NY,NX))+TX*ZMG1P(L,NY,NX)+CORP*ZMG1PH(L,NY,NX)
C
C     SALTS IN BAND
C
      H0POB(L,NY,NX)=TI*H0POB(L,NY,NX)+CORP*(FI*TH0POB
     2-TI*H0POB(L,NY,NX))+TX*H0POB(L,NY,NX)+CORP*H0POBH(L,NY,NX)
      H1POB(L,NY,NX)=TI*H1POB(L,NY,NX)+CORP*(FI*TH1POB
     2-TI*H1POB(L,NY,NX))+TX*H1POB(L,NY,NX)+CORP*H1POBH(L,NY,NX)
      H2POB(L,NY,NX)=TI*H2POB(L,NY,NX)+CORP*(FI*TH2POB
     2-TI*H2POB(L,NY,NX))+TX*H2POB(L,NY,NX)+CORP*H2POBH(L,NY,NX)
      H3POB(L,NY,NX)=TI*H3POB(L,NY,NX)+CORP*(FI*TH3POB
     2-TI*H3POB(L,NY,NX))+TX*H3POB(L,NY,NX)+CORP*H3POBH(L,NY,NX)
      ZFE1PB(L,NY,NX)=TI*ZFE1PB(L,NY,NX)+CORP*(FI*TFE1PB
     2-TI*ZFE1PB(L,NY,NX))+TX*ZFE1PB(L,NY,NX)+CORP*ZFE1BH(L,NY,NX)
      ZFE2PB(L,NY,NX)=TI*ZFE2PB(L,NY,NX)+CORP*(FI*TFE2PB
     2-TI*ZFE2PB(L,NY,NX))+TX*ZFE2PB(L,NY,NX)+CORP*ZFE2BH(L,NY,NX)
      ZCA0PB(L,NY,NX)=TI*ZCA0PB(L,NY,NX)+CORP*(FI*TCA0PB
     2-TI*ZCA0PB(L,NY,NX))+TX*ZCA0PB(L,NY,NX)+CORP*ZCA0BH(L,NY,NX)
      ZCA1PB(L,NY,NX)=TI*ZCA1PB(L,NY,NX)+CORP*(FI*TCA1PB
     2-TI*ZCA1PB(L,NY,NX))+TX*ZCA1PB(L,NY,NX)+CORP*ZCA1BH(L,NY,NX)
      ZCA2PB(L,NY,NX)=TI*ZCA2PB(L,NY,NX)+CORP*(FI*TCA2PB
     2-TI*ZCA2PB(L,NY,NX))+TX*ZCA2PB(L,NY,NX)+CORP*ZCA2BH(L,NY,NX)
      ZMG1PB(L,NY,NX)=TI*ZMG1PB(L,NY,NX)+CORP*(FI*TMG1PB
     2-TI*ZMG1PB(L,NY,NX))+TX*ZMG1PB(L,NY,NX)+CORP*ZMG1BH(L,NY,NX)
C
C     ADSORBED CATIONS AND ANIONS IN BAND, NON-BAND
C
      XN4(L,NY,NX)=TI*XN4(L,NY,NX)+CORP*(FI*TXNH4-TI*XN4(L,NY,NX))
     2+TX*XN4(L,NY,NX)
      XNB(L,NY,NX)=TI*XNB(L,NY,NX)+CORP*(FI*TXNHB-TI*XNB(L,NY,NX))
     2+TX*XNB(L,NY,NX)
      XHY(L,NY,NX)=TI*XHY(L,NY,NX)+CORP*(FI*TXHY-TI*XHY(L,NY,NX))
     2+TX*XHY(L,NY,NX)
      XAL(L,NY,NX)=TI*XAL(L,NY,NX)+CORP*(FI*TXAL-TI*XAL(L,NY,NX))
     2+TX*XAL(L,NY,NX)
      XFE(L,NY,NX)=TI*XFE(L,NY,NX)+CORP*(FI*TXFE-TI*XFE(L,NY,NX))
     2+TX*XFE(L,NY,NX)
      XCA(L,NY,NX)=TI*XCA(L,NY,NX)+CORP*(FI*TXCA-TI*XCA(L,NY,NX))
     2+TX*XCA(L,NY,NX)
      XMG(L,NY,NX)=TI*XMG(L,NY,NX)+CORP*(FI*TXMG-TI*XMG(L,NY,NX))
     2+TX*XMG(L,NY,NX)
      XNA(L,NY,NX)=TI*XNA(L,NY,NX)+CORP*(FI*TXNA-TI*XNA(L,NY,NX))
     2+TX*XNA(L,NY,NX)
      XKA(L,NY,NX)=TI*XKA(L,NY,NX)+CORP*(FI*TXKA-TI*XKA(L,NY,NX))
     2+TX*XKA(L,NY,NX)
      XHC(L,NY,NX)=TI*XHC(L,NY,NX)+CORP*(FI*TXHC-TI*XHC(L,NY,NX))
     2+TX*XHC(L,NY,NX)
      XALO2(L,NY,NX)=TI*XALO2(L,NY,NX)+CORP*(FI*TXAL2-TI*XALO2(L,NY,NX))
     2+TX*XALO2(L,NY,NX)
      XFEO2(L,NY,NX)=TI*XFEO2(L,NY,NX)+CORP*(FI*TXFE2-TI*XFEO2(L,NY,NX))
     2+TX*XFEO2(L,NY,NX)
      XOH0(L,NY,NX)=TI*XOH0(L,NY,NX)+CORP*(FI*TXOH0-TI*XOH0(L,NY,NX))
     2+TX*XOH0(L,NY,NX)
      XOH1(L,NY,NX)=TI*XOH1(L,NY,NX)+CORP*(FI*TXOH1-TI*XOH1(L,NY,NX))
     2+TX*XOH1(L,NY,NX)
      XOH2(L,NY,NX)=TI*XOH2(L,NY,NX)+CORP*(FI*TXOH2-TI*XOH2(L,NY,NX))
     2+TX*XOH2(L,NY,NX)
      XH1P(L,NY,NX)=TI*XH1P(L,NY,NX)+CORP*(FI*TXH1P-TI*XH1P(L,NY,NX))
     2+TX*XH1P(L,NY,NX)
      XH2P(L,NY,NX)=TI*XH2P(L,NY,NX)+CORP*(FI*TXH2P-TI*XH2P(L,NY,NX))
     2+TX*XH2P(L,NY,NX)
      XOH0B(L,NY,NX)=TI*XOH0B(L,NY,NX)+CORP*(FI*TXOH0B
     2-TI*XOH0B(L,NY,NX))+TX*XOH0B(L,NY,NX)
      XOH1B(L,NY,NX)=TI*XOH1B(L,NY,NX)+CORP*(FI*TXOH1B
     2-TI*XOH1B(L,NY,NX))+TX*XOH1B(L,NY,NX)
      XOH2B(L,NY,NX)=TI*XOH2B(L,NY,NX)+CORP*(FI*TXOH2B
     2-TI*XOH2B(L,NY,NX))+TX*XOH2B(L,NY,NX)
      XH1PB(L,NY,NX)=TI*XH1PB(L,NY,NX)+CORP*(FI*TXH1PB
     2-TI*XH1PB(L,NY,NX))+TX*XH1PB(L,NY,NX)
      XH2PB(L,NY,NX)=TI*XH2PB(L,NY,NX)+CORP*(FI*TXH2PB
     2-TI*XH2PB(L,NY,NX))+TX*XH2PB(L,NY,NX)
C
C     PRECIPITATED CATIONS AND ANIONS IN BAND, NON-BAND
C
      PALOH(L,NY,NX)=TI*PALOH(L,NY,NX)+CORP*(FI*TPALOH
     2-TI*PALOH(L,NY,NX))+TX*PALOH(L,NY,NX)
      PFEOH(L,NY,NX)=TI*PFEOH(L,NY,NX)+CORP*(FI*TPFEOH
     2-TI*PFEOH(L,NY,NX))+TX*PFEOH(L,NY,NX)
      PCACO(L,NY,NX)=TI*PCACO(L,NY,NX)+CORP*(FI*TPCACO
     2-TI*PCACO(L,NY,NX))+TX*PCACO(L,NY,NX)
      PCASO(L,NY,NX)=TI*PCASO(L,NY,NX)+CORP*(FI*TPCASO
     2-TI*PCASO(L,NY,NX))+TX*PCASO(L,NY,NX)
      PALPO(L,NY,NX)=TI*PALPO(L,NY,NX)+CORP*(FI*TPALPO
     2-TI*PALPO(L,NY,NX))+TX*PALPO(L,NY,NX)
      PFEPO(L,NY,NX)=TI*PFEPO(L,NY,NX)+CORP*(FI*TPFEPO
     2-TI*PFEPO(L,NY,NX))+TX*PFEPO(L,NY,NX)
      PCAPD(L,NY,NX)=TI*PCAPD(L,NY,NX)+CORP*(FI*TPCAPD
     2-TI*PCAPD(L,NY,NX))+TX*PCAPD(L,NY,NX)
      PCAPH(L,NY,NX)=TI*PCAPH(L,NY,NX)+CORP*(FI*TPCAPH
     2-TI*PCAPH(L,NY,NX))+TX*PCAPH(L,NY,NX)
      PCAPM(L,NY,NX)=TI*PCAPM(L,NY,NX)+CORP*(FI*TPCAPM
     2-TI*PCAPM(L,NY,NX))+TX*PCAPM(L,NY,NX)
      PALPB(L,NY,NX)=TI*PALPB(L,NY,NX)+CORP*(FI*TPALPB
     2-TI*PALPB(L,NY,NX))+TX*PALPB(L,NY,NX)
      PFEPB(L,NY,NX)=TI*PFEPB(L,NY,NX)+CORP*(FI*TPFEPB
     2-TI*PFEPB(L,NY,NX))+TX*PFEPB(L,NY,NX)
      PCPDB(L,NY,NX)=TI*PCPDB(L,NY,NX)+CORP*(FI*TPCPDB
     2-TI*PCPDB(L,NY,NX))+TX*PCPDB(L,NY,NX)
      PCPHB(L,NY,NX)=TI*PCPHB(L,NY,NX)+CORP*(FI*TPCPHB
     2-TI*PCPHB(L,NY,NX))+TX*PCPHB(L,NY,NX)
      PCPMB(L,NY,NX)=TI*PCPMB(L,NY,NX)+CORP*(FI*TPCPMB
     2-TI*PCPMB(L,NY,NX))+TX*PCPMB(L,NY,NX)
C
C     GASEOUS AND AQUEOUS GASES
C
      CO2G(L,NY,NX)=TI*CO2G(L,NY,NX)+CORP*(FI*TCO2G-TI*CO2G(L,NY,NX))
     2+TX*CO2G(L,NY,NX)
      CH4G(L,NY,NX)=TI*CH4G(L,NY,NX)+CORP*(FI*TCH4G-TI*CH4G(L,NY,NX))
     2+TX*CH4G(L,NY,NX)
      CO2S(L,NY,NX)=TI*CO2S(L,NY,NX)+CORP*(FI*TCOZS-TI*CO2S(L,NY,NX))
     2+TX*CO2S(L,NY,NX)+CORP*CO2SH(L,NY,NX)
      CH4S(L,NY,NX)=TI*CH4S(L,NY,NX)+CORP*(FI*TCHFS-TI*CH4S(L,NY,NX))
     2+TX*CH4S(L,NY,NX)+CORP*CH4SH(L,NY,NX)
      OXYG(L,NY,NX)=TI*OXYG(L,NY,NX)+CORP*(FI*TOXYG-TI*OXYG(L,NY,NX))
     2+TX*OXYG(L,NY,NX)
      OXYS(L,NY,NX)=TI*OXYS(L,NY,NX)+CORP*(FI*TOXYS-TI*OXYS(L,NY,NX))
     2+TX*OXYS(L,NY,NX)+CORP*OXYSH(L,NY,NX)
      Z2GG(L,NY,NX)=TI*Z2GG(L,NY,NX)+CORP*(FI*TZ2GG-TI*Z2GG(L,NY,NX))
     2+TX*Z2GG(L,NY,NX)
      Z2GS(L,NY,NX)=TI*Z2GS(L,NY,NX)+CORP*(FI*TZ2GS-TI*Z2GS(L,NY,NX))
     2+TX*Z2GS(L,NY,NX)+CORP*Z2GSH(L,NY,NX)
      Z2OG(L,NY,NX)=TI*Z2OG(L,NY,NX)+CORP*(FI*TZ2OG-TI*Z2OG(L,NY,NX))
     2+TX*Z2OG(L,NY,NX)
      Z2OS(L,NY,NX)=TI*Z2OS(L,NY,NX)+CORP*(FI*TZ2OS-TI*Z2OS(L,NY,NX))
     2+TX*Z2OS(L,NY,NX)+CORP*Z2OSH(L,NY,NX)
      ZNH3G(L,NY,NX)=TI*ZNH3G(L,NY,NX)+CORP*(FI*TZNH3G
     2-TI*ZNH3G(L,NY,NX))+TX*ZNH3G(L,NY,NX)
      H2GG(L,NY,NX)=TI*H2GG(L,NY,NX)+CORP*(FI*TH2GG-TI*H2GG(L,NY,NX))
     2+TX*H2GG(L,NY,NX)
      H2GS(L,NY,NX)=TI*H2GS(L,NY,NX)+CORP*(FI*TH2GS-TI*H2GS(L,NY,NX))
     2+TX*H2GS(L,NY,NX)+CORP*H2GSH(L,NY,NX)
C
C     SOLUTES IN MACROPORES
C
      ZNH4SH(L,NY,NX)=XCORP(NY,NX)*ZNH4SH(L,NY,NX)
      ZNH3SH(L,NY,NX)=XCORP(NY,NX)*ZNH3SH(L,NY,NX)
      ZNO3SH(L,NY,NX)=XCORP(NY,NX)*ZNO3SH(L,NY,NX)
      ZNO2SH(L,NY,NX)=XCORP(NY,NX)*ZNO2SH(L,NY,NX)
      H1PO4H(L,NY,NX)=XCORP(NY,NX)*H1PO4H(L,NY,NX)
      H2PO4H(L,NY,NX)=XCORP(NY,NX)*H2PO4H(L,NY,NX)
      ZNH4BH(L,NY,NX)=XCORP(NY,NX)*ZNH4BH(L,NY,NX)
      ZNH3BH(L,NY,NX)=XCORP(NY,NX)*ZNH3BH(L,NY,NX)
      ZNO3BH(L,NY,NX)=XCORP(NY,NX)*ZNO3BH(L,NY,NX)
      ZNO2BH(L,NY,NX)=XCORP(NY,NX)*ZNO2BH(L,NY,NX)
      H1POBH(L,NY,NX)=XCORP(NY,NX)*H1POBH(L,NY,NX)
      H2POBH(L,NY,NX)=XCORP(NY,NX)*H2POBH(L,NY,NX)
      ZALH(L,NY,NX)=XCORP(NY,NX)*ZALH(L,NY,NX)
      ZFEH(L,NY,NX)=XCORP(NY,NX)*ZFEH(L,NY,NX)
      ZHYH(L,NY,NX)=XCORP(NY,NX)*ZHYH(L,NY,NX)
      ZCCH(L,NY,NX)=XCORP(NY,NX)*ZCCH(L,NY,NX)
      ZMAH(L,NY,NX)=XCORP(NY,NX)*ZMAH(L,NY,NX)
      ZNAH(L,NY,NX)=XCORP(NY,NX)*ZNAH(L,NY,NX)
      ZKAH(L,NY,NX)=XCORP(NY,NX)*ZKAH(L,NY,NX)
      ZOHH(L,NY,NX)=XCORP(NY,NX)*ZOHH(L,NY,NX)
      ZSO4H(L,NY,NX)=XCORP(NY,NX)*ZSO4H(L,NY,NX)
      ZCLH(L,NY,NX)=XCORP(NY,NX)*ZCLH(L,NY,NX)
      ZCO3H(L,NY,NX)=XCORP(NY,NX)*ZCO3H(L,NY,NX)
      ZHCO3H(L,NY,NX)=XCORP(NY,NX)*ZHCO3H(L,NY,NX)
      ZALO1H(L,NY,NX)=XCORP(NY,NX)*ZALO1H(L,NY,NX)
      ZALO2H(L,NY,NX)=XCORP(NY,NX)*ZALO2H(L,NY,NX)
      ZALO3H(L,NY,NX)=XCORP(NY,NX)*ZALO3H(L,NY,NX)
      ZALO4H(L,NY,NX)=XCORP(NY,NX)*ZALO4H(L,NY,NX)
      ZALSH(L,NY,NX)=XCORP(NY,NX)*ZALSH(L,NY,NX)
      ZFEO1H(L,NY,NX)=XCORP(NY,NX)*ZFEO1H(L,NY,NX)
      ZFEO2H(L,NY,NX)=XCORP(NY,NX)*ZFEO2H(L,NY,NX)
      ZFEO3H(L,NY,NX)=XCORP(NY,NX)*ZFEO3H(L,NY,NX)
      ZFEO4H(L,NY,NX)=XCORP(NY,NX)*ZFEO4H(L,NY,NX)
      ZFESH(L,NY,NX)=XCORP(NY,NX)*ZFESH(L,NY,NX)
      ZCAOH(L,NY,NX)=XCORP(NY,NX)*ZCAOH(L,NY,NX)
      ZCACH(L,NY,NX)=XCORP(NY,NX)*ZCACH(L,NY,NX)
      ZCAHH(L,NY,NX)=XCORP(NY,NX)*ZCAHH(L,NY,NX)
      ZCASH(L,NY,NX)=XCORP(NY,NX)*ZCASH(L,NY,NX)
      ZMGOH(L,NY,NX)=XCORP(NY,NX)*ZMGOH(L,NY,NX)
      ZMGCH(L,NY,NX)=XCORP(NY,NX)*ZMGCH(L,NY,NX)
      ZMGHH(L,NY,NX)=XCORP(NY,NX)*ZMGHH(L,NY,NX)
      ZMGSH(L,NY,NX)=XCORP(NY,NX)*ZMGSH(L,NY,NX)
      ZNACH(L,NY,NX)=XCORP(NY,NX)*ZNACH(L,NY,NX)
      ZNASH(L,NY,NX)=XCORP(NY,NX)*ZNASH(L,NY,NX)
      ZKASH(L,NY,NX)=XCORP(NY,NX)*ZKASH(L,NY,NX)
      H0PO4H(L,NY,NX)=XCORP(NY,NX)*H0PO4H(L,NY,NX)
      H3PO4H(L,NY,NX)=XCORP(NY,NX)*H3PO4H(L,NY,NX)
      ZFE1PH(L,NY,NX)=XCORP(NY,NX)*ZFE1PH(L,NY,NX)
      ZFE2PH(L,NY,NX)=XCORP(NY,NX)*ZFE2PH(L,NY,NX)
      ZCA0PH(L,NY,NX)=XCORP(NY,NX)*ZCA0PH(L,NY,NX)
      ZCA1PH(L,NY,NX)=XCORP(NY,NX)*ZCA1PH(L,NY,NX)
      ZCA2PH(L,NY,NX)=XCORP(NY,NX)*ZCA2PH(L,NY,NX)
      ZMG1PH(L,NY,NX)=XCORP(NY,NX)*ZMG1PH(L,NY,NX)
      H0POBH(L,NY,NX)=XCORP(NY,NX)*H0POBH(L,NY,NX)
      H1POBH(L,NY,NX)=XCORP(NY,NX)*H1POBH(L,NY,NX)
      H3POBH(L,NY,NX)=XCORP(NY,NX)*H3POBH(L,NY,NX)
      ZFE1BH(L,NY,NX)=XCORP(NY,NX)*ZFE1BH(L,NY,NX)
      ZFE2BH(L,NY,NX)=XCORP(NY,NX)*ZFE2BH(L,NY,NX)
      ZCA0BH(L,NY,NX)=XCORP(NY,NX)*ZCA0BH(L,NY,NX)
      ZCA1BH(L,NY,NX)=XCORP(NY,NX)*ZCA1BH(L,NY,NX)
      ZCA2BH(L,NY,NX)=XCORP(NY,NX)*ZCA2BH(L,NY,NX)
      ZMG1BH(L,NY,NX)=XCORP(NY,NX)*ZMG1BH(L,NY,NX)
      CO2SH(L,NY,NX)=XCORP(NY,NX)*CO2SH(L,NY,NX)
      CH4SH(L,NY,NX)=XCORP(NY,NX)*CH4SH(L,NY,NX)
      OXYSH(L,NY,NX)=XCORP(NY,NX)*OXYSH(L,NY,NX)
      Z2GSH(L,NY,NX)=XCORP(NY,NX)*Z2GSH(L,NY,NX)
      Z2OSH(L,NY,NX)=XCORP(NY,NX)*Z2OSH(L,NY,NX)
      H2GSH(L,NY,NX)=XCORP(NY,NX)*H2GSH(L,NY,NX)
C
C     ORGANIC MATTER
C
      DO 5965 K=0,5
      DO 5965 N=1,7
      DO 5965 M=1,3
      OMC(M,N,K,L,NY,NX)=TI*OMC(M,N,K,L,NY,NX)+CORP*(FI*TOMC(M,N,K)
     2-TI*OMC(M,N,K,L,NY,NX))+TX*OMC(M,N,K,L,NY,NX)
      OMN(M,N,K,L,NY,NX)=TI*OMN(M,N,K,L,NY,NX)+CORP*(FI*TOMN(M,N,K)
     2-TI*OMN(M,N,K,L,NY,NX))+TX*OMN(M,N,K,L,NY,NX)
      OMP(M,N,K,L,NY,NX)=TI*OMP(M,N,K,L,NY,NX)+CORP*(FI*TOMP(M,N,K)
     2-TI*OMP(M,N,K,L,NY,NX))+TX*OMP(M,N,K,L,NY,NX)
5965  CONTINUE
      DO 5980 K=0,4
      DO 5975 M=1,2
      ORC(M,K,L,NY,NX)=TI*ORC(M,K,L,NY,NX)+CORP*(FI*TORC(M,K)
     2-TI*ORC(M,K,L,NY,NX))+TX*ORC(M,K,L,NY,NX)
      ORN(M,K,L,NY,NX)=TI*ORN(M,K,L,NY,NX)+CORP*(FI*TORN(M,K)
     2-TI*ORN(M,K,L,NY,NX))+TX*ORN(M,K,L,NY,NX)
      ORP(M,K,L,NY,NX)=TI*ORP(M,K,L,NY,NX)+CORP*(FI*TORP(M,K)
     2-TI*ORP(M,K,L,NY,NX))+TX*ORP(M,K,L,NY,NX)
5975  CONTINUE
      OQC(K,L,NY,NX)=TI*OQC(K,L,NY,NX)+CORP*(FI*TOQC(K)
     2-TI*OQC(K,L,NY,NX))+TX*OQC(K,L,NY,NX)+CORP*OQCH(K,L,NY,NX)
      OQN(K,L,NY,NX)=TI*OQN(K,L,NY,NX)+CORP*(FI*TOQN(K)
     2-TI*OQN(K,L,NY,NX))+TX*OQN(K,L,NY,NX)+CORP*OQNH(K,L,NY,NX)
      OQP(K,L,NY,NX)=TI*OQP(K,L,NY,NX)+CORP*(FI*TOQP(K)
     2-TI*OQP(K,L,NY,NX))+TX*OQP(K,L,NY,NX)+CORP*OQPH(K,L,NY,NX)
      OQA(K,L,NY,NX)=TI*OQA(K,L,NY,NX)+CORP*(FI*TOQA(K)
     2-TI*OQA(K,L,NY,NX))+TX*OQA(K,L,NY,NX)+CORP*OQAH(K,L,NY,NX)
      OQCH(K,L,NY,NX)=XCORP(NY,NX)*OQCH(K,L,NY,NX)
      OQNH(K,L,NY,NX)=XCORP(NY,NX)*OQNH(K,L,NY,NX)
      OQPH(K,L,NY,NX)=XCORP(NY,NX)*OQPH(K,L,NY,NX)
      OQAH(K,L,NY,NX)=XCORP(NY,NX)*OQAH(K,L,NY,NX)
      OHC(K,L,NY,NX)=TI*OHC(K,L,NY,NX)+CORP*(FI*TOHC(K)
     2-TI*OHC(K,L,NY,NX))+TX*OHC(K,L,NY,NX)
      OHN(K,L,NY,NX)=TI*OHN(K,L,NY,NX)+CORP*(FI*TOHN(K)
     2-TI*OHN(K,L,NY,NX))+TX*OHN(K,L,NY,NX)
      OHP(K,L,NY,NX)=TI*OHP(K,L,NY,NX)+CORP*(FI*TOHP(K)
     2-TI*OHP(K,L,NY,NX))+TX*OHP(K,L,NY,NX)
      OHA(K,L,NY,NX)=TI*OHA(K,L,NY,NX)+CORP*(FI*TOHA(K)
     2-TI*OHA(K,L,NY,NX))+TX*OHA(K,L,NY,NX)
      DO 5970 M=1,5
      OSC(M,K,L,NY,NX)=TI*OSC(M,K,L,NY,NX)+CORP*(FI*TOSC(M,K)
     2-TI*OSC(M,K,L,NY,NX))+TX*OSC(M,K,L,NY,NX)
      OSA(M,K,L,NY,NX)=TI*OSA(M,K,L,NY,NX)+CORP*(FI*TOSA(M,K)
     2-TI*OSA(M,K,L,NY,NX))+TX*OSA(M,K,L,NY,NX)
      OSN(M,K,L,NY,NX)=TI*OSN(M,K,L,NY,NX)+CORP*(FI*TOSN(M,K)
     2-TI*OSN(M,K,L,NY,NX))+TX*OSN(M,K,L,NY,NX)
      OSP(M,K,L,NY,NX)=TI*OSP(M,K,L,NY,NX)+CORP*(FI*TOSP(M,K)
     2-TI*OSP(M,K,L,NY,NX))+TX*OSP(M,K,L,NY,NX)
5970  CONTINUE
5980  CONTINUE
C
C     ADD STATE VARIABLES IN SURFACE RESIDUE INCORPORATED
C     WITHIN TILLAGE MIXING ZONE
C
      DO 5910 K=0,5
      IF(K.NE.3.AND.K.NE.4)THEN
      DO 5915 N=1,7
      DO 5915 M=1,3
      OMC(M,N,K,L,NY,NX)=OMC(M,N,K,L,NY,NX)+FI*TOMGC(M,N,K)
      OMN(M,N,K,L,NY,NX)=OMN(M,N,K,L,NY,NX)+FI*TOMGN(M,N,K)
      OMP(M,N,K,L,NY,NX)=OMP(M,N,K,L,NY,NX)+FI*TOMGP(M,N,K)
5915  CONTINUE
      ENDIF
5910  CONTINUE
      DO 5920 K=0,2
      DO 5925 M=1,2
      ORC(M,K,L,NY,NX)=ORC(M,K,L,NY,NX)+FI*TORXC(M,K)
      ORN(M,K,L,NY,NX)=ORN(M,K,L,NY,NX)+FI*TORXN(M,K)
      ORP(M,K,L,NY,NX)=ORP(M,K,L,NY,NX)+FI*TORXP(M,K)
5925  CONTINUE
      OQC(K,L,NY,NX)=OQC(K,L,NY,NX)+FI*TOQGC(K)
      OQN(K,L,NY,NX)=OQN(K,L,NY,NX)+FI*TOQGN(K)
      OQP(K,L,NY,NX)=OQP(K,L,NY,NX)+FI*TOQGP(K)
      OQA(K,L,NY,NX)=OQA(K,L,NY,NX)+FI*TOQGA(K)
      OQCH(K,L,NY,NX)=OQCH(K,L,NY,NX)+FI*TOQHC(K)
      OQNH(K,L,NY,NX)=OQNH(K,L,NY,NX)+FI*TOQHN(K)
      OQPH(K,L,NY,NX)=OQPH(K,L,NY,NX)+FI*TOQHP(K)
      OQAH(K,L,NY,NX)=OQAH(K,L,NY,NX)+FI*TOQHA(K)
      OHC(K,L,NY,NX)=OHC(K,L,NY,NX)+FI*TOHGC(K)
      OHN(K,L,NY,NX)=OHN(K,L,NY,NX)+FI*TOHGN(K)
      OHP(K,L,NY,NX)=OHP(K,L,NY,NX)+FI*TOHGP(K)
      OHA(K,L,NY,NX)=OHA(K,L,NY,NX)+FI*TOHGA(K)
      DO 5930 M=1,5
      OSC(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)+FI*TOSGC(M,K)
      OSA(M,K,L,NY,NX)=OSA(M,K,L,NY,NX)+FI*TOSGA(M,K)
      OSN(M,K,L,NY,NX)=OSN(M,K,L,NY,NX)+FI*TOSGN(M,K)
      OSP(M,K,L,NY,NX)=OSP(M,K,L,NY,NX)+FI*TOSGP(M,K)
5930  CONTINUE
5920  CONTINUE
C
C     RECALCULATE SOIL ORGANIC MATTER IN SOIL LAYERS WITHIN THE MIXING ZONE
C
      DC=0.0
      DN=0.0
      DP=0.0
      OC=0.0
      ON=0.0
      OP=0.0
      DCC=0.0
      DNC=0.0
      DPC=0.0
      OCC=0.0
      ONC=0.0
      OPC=0.0
      DO 5985 K=0,5
      DO 5985 N=1,7
      DO 5985 M=1,3
      OC=OC+OMC(M,N,K,L,NY,NX)
      ON=ON+OMN(M,N,K,L,NY,NX)
      OP=OP+OMP(M,N,K,L,NY,NX)
      IF(K.LE.2)THEN
      DC=DC+OMC(M,N,K,L,NY,NX)
      DN=DN+OMN(M,N,K,L,NY,NX)
      DP=DP+OMP(M,N,K,L,NY,NX)
      ENDIF
5985  CONTINUE
      DO 6995 K=0,4
      DO 6985 M=1,2
      OC=OC+ORC(M,K,L,NY,NX)
      ON=ON+ORN(M,K,L,NY,NX)
      OP=OP+ORP(M,K,L,NY,NX)
      IF(K.LE.2)THEN
      DC=DC+ORC(M,K,L,NY,NX)
      DN=DN+ORN(M,K,L,NY,NX)
      DP=DP+ORP(M,K,L,NY,NX)
      ENDIF
6985  CONTINUE
      OC=OC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
     2+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      ON=ON+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      OP=OP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
      IF(K.LE.2)THEN
      DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
     2+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      DN=DN+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
      DC=DC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX)
      ENDIF
      DO 6980 M=1,5
      IF(K.LE.2)THEN
      IF(M.LE.4)THEN
      DC=DC+OSC(M,K,L,NY,NX)
      DN=DN+OSN(M,K,L,NY,NX)
      DP=DP+OSP(M,K,L,NY,NX)
      ELSE
      DCC=DCC+OSC(M,K,L,NY,NX)
      DNC=DNC+OSN(M,K,L,NY,NX)
      DPC=DP+OSP(M,K,L,NY,NX)
      ENDIF
      ELSE
      IF(M.LE.4)THEN
      OC=OC+OSC(M,K,L,NY,NX)
      ON=ON+OSN(M,K,L,NY,NX)
      OP=OP+OSP(M,K,L,NY,NX)
      ELSE
      OCC=OCC+OSC(M,K,L,NY,NX)
      ONC=ONC+OSN(M,K,L,NY,NX)
      OPC=OPC+OSP(M,K,L,NY,NX)
      ENDIF
      ENDIF
6980  CONTINUE
6995  CONTINUE
      ORGC(L,NY,NX)=DC+OC
      ORGN(L,NY,NX)=DN+ON
      ORGCC(L,NY,NX)=DCC+OCC
      ORGNC(L,NY,NX)=DNC+ONC
      ORGR(L,NY,NX)=DC
      CO2S(L,NY,NX)=CO2S(L,NY,NX)+FI*TCO2GS
      CH4S(L,NY,NX)=CH4S(L,NY,NX)+FI*TCH4GS
      OXYS(L,NY,NX)=OXYS(L,NY,NX)+FI*TOXYGS
      Z2GS(L,NY,NX)=Z2GS(L,NY,NX)+FI*TZ2GSG
      Z2OS(L,NY,NX)=Z2OS(L,NY,NX)+FI*TZ2OGS
      H2GS(L,NY,NX)=H2GS(L,NY,NX)+FI*TH2GGS
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+FI*TNH4GS
      ZNH3S(L,NY,NX)=ZNH3S(L,NY,NX)+FI*TNH3GS
      ZNO3S(L,NY,NX)=ZNO3S(L,NY,NX)+FI*TNO3GS
      ZNO2S(L,NY,NX)=ZNO2S(L,NY,NX)+FI*TNO2GS
      H1PO4(L,NY,NX)=H1PO4(L,NY,NX)+FI*TP14GS
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+FI*TPO4GS
      XN4(L,NY,NX)=XN4(L,NY,NX)+FI*TXN4G
      XOH0(L,NY,NX)=XOH0(L,NY,NX)+FI*TXOH0G
      XOH1(L,NY,NX)=XOH1(L,NY,NX)+FI*TXOH1G
      XOH2(L,NY,NX)=XOH2(L,NY,NX)+FI*TXOH2G
      XH1P(L,NY,NX)=XH1P(L,NY,NX)+FI*TXH1PG
      XH2P(L,NY,NX)=XH2P(L,NY,NX)+FI*TXH2PG
      PALPO(L,NY,NX)=PALPO(L,NY,NX)+FI*TALPOG
      PFEPO(L,NY,NX)=PFEPO(L,NY,NX)+FI*TFEPOG
      PCAPD(L,NY,NX)=PCAPD(L,NY,NX)+FI*TCAPDG
      PCAPH(L,NY,NX)=PCAPH(L,NY,NX)+FI*TCAPHG
      PCAPM(L,NY,NX)=PCAPM(L,NY,NX)+FI*TCAPMG
      ZNH4FA(L,NY,NX)=ZNH4FA(L,NY,NX)+FI*TNH4FG
      ZNH3FA(L,NY,NX)=ZNH3FA(L,NY,NX)+FI*TNH3FG
      ZNHUFA(L,NY,NX)=ZNHUFA(L,NY,NX)+FI*TNHUFG
      ZNO3FA(L,NY,NX)=ZNO3FA(L,NY,NX)+FI*TNO3FG
      ZNHU0(L,NY,NX)=ZNHUX0
      ZNHUI(L,NY,NX)=ZNHUXI
      ZNFN0(L,NY,NX)=ZNFNX0
      ZNFNI(L,NY,NX)=(TI*ZNFNI(L,NY,NX)+CORP*(FI*TZNFNI
     2-TI*ZNFNI(L,NY,NX))+TX*ZNFNI(L,NY,NX)+FI*TZNFNG)/FI
      TZNFN2=TZNFN2+ZNFNI(L,NY,NX)
      ENDIF
2000  CONTINUE
      ZNFN0(0,NY,NX)=ZNFNX0
      ZNFNI(0,NY,NX)=ZNFNI(0,NY,NX)*XCORP0
      TZNFN2=TZNFN2+TZNFNG
      TZNFNI=TZNFNI+TZNFNG
      DO 2001 L=NU(NY,NX),LL
      IF(TZNFN2.GT.ZERO)THEN
      ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)*TZNFNI/TZNFN2
      ZNFNI(L,NY,NX)=ZNFNI(L,NY,NX)
     2+0.5*(ZNFN0(L,NY,NX)-ZNFNI(L,NY,NX))
      ENDIF
2001  CONTINUE
      ENDIF
      ENDIF
C
C     CHECK MATERIAL BALANCES FROM TOTAL ECOSYSTEM CONTENTS
C     LESS CUMULATIVE INPUTS PLUS CUMULATIVE OUTPUTS
C
      IF(I.EQ.365.AND.J.EQ.24.AND.NFZ.EQ.NFH)THEN
      WRITE(19,2221)'ORGC',I,J,IYRC,NX,NY
     2,(ORGC(L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
      WRITE(20,2221)'ORGN',I,J,IYRC,NX,NY
     2,(ORGN(L,NY,NX)/AREA(3,L,NY,NX),L=0,NL(NY,NX))
2221  FORMAT(A8,5I6,21E14.6)
      ENDIF
C     IF(J.LE.2)THEN
C     WRITE(20,2221)'OMCL',I,J,IYRC,NX,NY
C    2,(OMCL(L,NY,NX),L=0,NL(NY,NX))
C     WRITE(20,2221)'OMNL',I,J,IYRC,NX,NY
C    2,(OMNL(L,NY,NX),L=0,NL(NY,NX))
C
C     CHECK C BALANCE
C
C     WRITE(20,2222)'TLC',I,J,NFZ,NX,NY
C    2,TLRSDC+TLORGC+TLCO2G-CO2GIN+TCOU-TORGF-XCSN
C    2,TLRSDC,TLORGC,TLCO2G,CO2GIN,TCOU,TORGF,XCSN
C    5,XCNET(NY,NX),XHNET(NY,NX)
C    3,(RCGOX(L,NY,NX),L=0,NL(NY,NX))
C    3,(RCHOX(L,NY,NX),L=0,NL(NY,NX))
C    3,(RC4OX(L,NY,NX),L=0,NL(NY,NX))
C    3,(ORGC(L,NY,NX),L=0,NL(NY,NX))
C    3,(ORGCC(L,NY,NX),L=0,NL(NY,NX))
C    5,XCODFS(NY,NX),XCOFLG(3,NU(NY,NX),NY,NX),TCO2Z(NY,NX)
C    3,XCODFG(0,NY,NX)
C    2,FLQGQ(NY,NX)*CCOR(NY,NX),FLQGI(NY,NX)*CCOQ(NY,NX)
C    3,XCODFR(NY,NX),XCHDFS(NY,NX),XCHFLG(3,NU(NY,NX),NY,NX)
C    3,XCHDFG(0,NY,NX)
C    2,FLQGQ(NY,NX)*CCHR(NY,NX),FLQGI(NY,NX)*CCHQ(NY,NX)
C    3,XCHDFR(NY,NX),PRECU(NY,NX)*CCOQ(NY,NX)*XNFH
C    4,PRECU(NY,NX)*CCHQ(NY,NX)*XNFH
C    6,TCOQRS(NY,NX),TCHQRS(NY,NX),XCOFLS(1,0,NY,NX+1)
C    7,XCOFLS(2,0,NY+1,NX)
C    3,UCOP(NY,NX),UDOCQ(NY,NX),UDICQ(NY,NX),UDOCD(NY,NX),UDICD(NY,NX)
C    2,(((CSNT(M,K,L,NY,NX),M=1,4),K=0,1),L=0,NJ(NY,NX))
C    3,(TCO2P(L,NY,NX),L=1,NJ(NY,NX)),(TCO2S(L,NY,NX),L=1,NJ(NY,NX))
C    4,CQ,ZCSNC(NY,NX)
C
C     CHECK WATER BALANCE
C
C     WRITE(20,2222)'TLW',I,J,NFZ,NX,NY
C    2,VOLWSO-CRAIN+CRUN+CEVAP+VOLWOU
C    2,VOLWSO,CRAIN,CRUN,CEVAP,VOLWOU
C    3,TVOLWC(NY,NX),TVOLWP(NY,NX),UVOLO(NY,NX)
C    3,UVOLW(NY,NX),UEVAP(NY,NX),URUN(NY,NX),FLWR(NY,NX),TFLWC(NY,NX)
C    4,(VOLWC(NZ,NY,NX),NZ=1,3)
C    4,(VOLWP(NZ,NY,NX),NZ=1,3)
C    4,(VOLWQ(NZ,NY,NX),NZ=1,3)
C    4,VOLW(0,NY,NX),VOLV(0,NY,NX),VOLI(0,NY,NX)*DENSI
C    4,TEVAPG(NY,NX),TEVAPC(NY,NX),TEVAPP(NY,NX)
C    4,TEVAPP(NY,NX)-TEVAPC(NY,NX),XVPFLG(3,NU(NY,NX),NY,NX)
C    5,VOLSS(NY,NX),VOLWS(NY,NX),VOLIS(NY,NX)*DENSI
C    6,(VOLW(L,NY,NX),L=0,NL(NY,NX))
C    6,(VOLV(L,NY,NX),L=0,NL(NY,NX))
C    6,(VOLWH(L,NY,NX),L=NU(NY,NX),NL(NY,NX))
C    6,(VOLI(L,NY,NX)*DENSI,L=0,NL(NY,NX))
C    6,(VOLIH(L,NY,NX)*DENSI,L=NU(NY,NX),NL(NY,NX))
C    6,TQS(NY,NX),TQW(NY,NX),TQI(NY,NX)
C    8,(TUPWTR(L,NY,NX),L=1,NL(NY,NX))
C
C     CHECK HEAT BALANCE
C
C     WRITE(20,2222)'TLH',I,J,NFZ,NX,NY
C    2,HEATSO-HEATIN+HEATOU
C    2,HEATSO,HEATIN,HEATOU
C    I,(HCBFX(L,2,NX),L=0,NL(NY,NX))
C    I,HFLXD,HFLXO,HFLWR(NY,NX)
C    I,4.19*TKAM(NY,NX)*PRECA(NY,NX)*XNFH
C    I,2.095*TKAM(NY,NX)*PRECW(NY,NX)*XNFH
C    I,XHFLVR(NY,NX),XHFLFR(NY,NX)
C    I,HEATH(NY,NX),THFLXC(NY,NX)
C    I,(XHFLF0(L,NY,NX),L=1,JS)
C    I,(XHFLV0(L,NY,NX),L=1,JS)
C    I,(THFLVL(L,NY,NX),L=NU(NY,NX),NL(NY,NX))
C    I,(THFLFL(L,NY,NX),L=NU(NY,NX),NL(NY,NX))
C    I,(TUPHT(L,NY,NX),L=NU(NY,NX),NL(NY,NX))
C    O,4.19*TKAM(NY,NX)*PRECU(NY,NX)*XNFH,THQS(NY,NX)
C    S,((VHCP(L,NY,NX)*TKS(L,NY,NX)),L=0,NL(NY,NX))
C    S,TENGYC(NY,NX)
C    S,((VHCPW(L,NY,NX)*TKW(L,NY,NX)),L=1,JS)
C    S,VHCP(0,NY,NX),VHCPRX(NY,NX),ORGC(0,NY,NX),ORGCC(0,NY,NX)
C    S,VOLW(0,NY,NX),VOLV(0,NY,NX),VOLI(0,NY,NX)
C
C     CHECK O2 BALANCE
C
C     WRITE(19,2222)'TLO',I,J,NFZ,NX,NY
C    2,OXYGSO-OXYGIN+OXYGOU
C    2,OXYGSO,OXYGIN,OXYGOU
C    3,XOXDFS(NY,NX),XOXDFR(NY,NX),XOXFLG(3,NU(NY,NX),NY,NX)
C    3,XOXFLG(3,0,NY,NX),TOXYZ(NY,NX)
C    3,FLQGQ(NY,NX)*COXR(NY,NX)
C    4,FLQGI(NY,NX)*COXQ,PRECU(NY,NX)*COXQ*XNFH
C    3,(ROGOX(L,NY,NX),L=0,NL(NY,NX))
C    3,(RC4OX(L,NY,NX),L=0,NL(NY,NX))
C    3,(OXYG(L,NY,NX),L=0,NL(NY,NX))
C    4,(OXYS(L,NY,NX),L=0,NL(NY,NX))
C    4,(OXYSH(L,NY,NX),L=1,JZ)
C    2,(RUPOXO(L,NY,NX),L=1,NJ(NY,NX))
C    3,(TUPOXP(L,NY,NX),L=1,NJ(NY,NX)),(TOXFLA(L,NY,NX),L=1,NJ(NY,NX))
C
C     CHECK N BALANCE
C
C     WRITE(20,2222)'TLN',I,J,NFZ,NX,NY
C    2,TLRSDN+TLORGN+TLN2G+TLNH4+TLNO3-ZN2GIN-TZIN+TZOU-TORGN-XZSN
C    2,TLRSDN,TLORGN,TLN2G,TLNH4,TLNO3,ZN2GIN,TZIN,TZOU,TORGN,XZSN
C    5,(XN4(L,NY,NX),L=0,NL(NY,NX))
C    4,(((ZSNT(M,K,L,NY,NX),M=1,5),K=0,1),L=0,JZ)
C    5,(TUPNH4(L,NY,NX),L=1,JZ)
C    6,(TUPNO3(L,NY,NX),L=1,JZ)
C    7,(TNHFLA(L,NY,NX),L=1,JZ)
C    7,XN3DFS(NY,NX),XNBDFS(NY,NX)
C    8,XN3FLG(3,NU(NY,NX),NY,NX),TNH3Z(NY,NX),UN2GS(NY,NX)
C    9,(XN2GS(L,NY,NX),L=0,JZ)
C    4,PRECQ(NY,NX)*XNFH,PRECA(NY,NX)*XNFH
C    4,PRECW(NY,NX)*XNFH,PRECI(NY,NX)*XNFH,FLQGM(NY,NX),FLQRM(NY,NX)
C
C     CHECK P BALANCE
C
C     WRITE(20,2222)'TLP',I,J,NFZ,NX,NY
C    2,TLRSDP+TLORGP+TLPO4-TPIN+TPOU-TORGP-XPSN
C    2,TLRSDP,TLORGP,TLPO4,TPIN,TPOU,TORGP,XPSN
C    2,(Z1PW(L,NY,NX),L=1,JS),(ZHPW(L,NY,NX),L=1,JS)
C    3,(H1PO4(L,NY,NX),L=0,JZ),(H2PO4(L,NY,NX),L=0,JZ)
C    4,(H1PO4H(L,NY,NX),L=1,JZ),(H2PO4H(L,NY,NX),L=1,JZ)
C    6,FLQGQ(NY,NX),FLQRQ(NY,NX),CPOR(NY,NX),CH1PR(NY,NX)
C    7,FLQGI(NY,NX),FLQRI(NY,NX),CPOQ(I,NY,NX),CH1PQ(I,NY,NX)
C
C     CHECK SALT BALANCE
C
C     WRITE(20,2222)'TLI',I,J,NFZ,NX,NY
C    2,TION-TIONIN+TIONOU
C    2,TION,TIONIN,TIONOU
C    3,PRECQ(NY,NX)*XNFH,XHGDFS(NY,NX),XHGFLG(3,NU(NY,NX),NY,NX)
C    4,TH2GZ(NY,NX)
C    4,(XHGQRS(N,NY,NX),N=1,2),(RH2GO(L,NY,NX),L=1,JZ)
C    5,(THGFLA(L,NY,NX),L=1,JZ)
C    6,(H2GG(L,NY,NX),L=1,JZ),(TLH2GP(L,NY,NX),L=1,JZ)
C
C     CHECK H2 BALANCE
C
C     WRITE(20,2222)'TLG',I,J,NFZ,NX,NY
C    2,TLH2G-H2GIN+H2GOU
C    3,TLH2G,H2GIN,H2GOU
C    4,(H2GG(L,NY,NX),L=0,NJ(NY,NX)),(H2GS(L,NY,NX),L=0,NJ(NY,NX))
C    3,(H2GSH(L,NY,NX),L=1,NJ(NY,NX)),(TLH2GP(L,NY,NX),L=1,NJ(NY,NX))
C    4,XHGDFS(NY,NX),TH2GZ(NY,NX),(THGFLA(L,NY,NX),L=1,NJ(NY,NX))
C    2,XHGDFG(0,NY,NX),XHGDFR(NY,NX),(XHGBBL(L,NY,NX),L=1,NJ(NY,NX))
C    3,(RH2GO(L,NY,NX),L=0,NJ(NY,NX)),(TUPHGS(L,NY,NX),L=1,NJ(NY,NX))
C    4,(XHGQRS(N,NY,NX),N=1,2)
C    5,((XHGFLS(N,L,NY,NX),N=1,3),L=0,NJ(NY,NX))
C    5,((XHGFHS(N,L,NY,NX),N=1,3),L=0,NJ(NY,NX))
C    6,((XHGFLG(N,L,NY,NX),N=1,3),L=0,NJ(NY,NX))
C
C     CHECK SEDIMENT BALANCE
C
C     WRITE(*,2222)'TLS',I,J,NFZ,NX,NY
C    2,TSEDSO+TSEDOU
C    2,TSEDSO,TSEDOU,USEDOU(NY,NX),DLYR(3,NU(NY,NX),NY,NX)
C    3,BKVL(NU(NY,NX),NY,NX),SAND(NU(NY,NX),NY,NX)
C    4,SILT(NU(NY,NX),NY,NX),CLAY(NU(NY,NX),NY,NX)
C    5,ORGC(NU(NY,NX),NY,NX)
2222  FORMAT(A8,5I4,360F20.6)
C     ENDIF
9990  CONTINUE
9995  CONTINUE
      RETURN
      END
