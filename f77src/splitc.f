
      SUBROUTINE splitc(NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
      include "parameters.h"
      include "filec.h"
      include "files.h"
      include "blkc.h"
      CHARACTER*16 DATA(30),DATAC(30,250,250),DATAP(JP,JY,JX)
     2,DATAM(JP,JY,JX),DATAX(JP),DATAY(JP),DATAZ(JP,JY,JX)
     3,OUTS(10),OUTP(10),OUTFILS(10,JY,JX),OUTFILP(10,JP,JY,JX)
      CHARACTER*3 CHOICE(102,20)
      CHARACTER*8 CDATE
      character*1024 str
      integer nz,nx,ny,n
      integer :: failure
      character(len=*), parameter :: modfile=__FILE__
      nz=1
      do nx=nhw,nhe
         do ny=nvn,nvs
	    if(nz .lt. np0(ny,nx))nz=np0(ny,nx)
         enddo
      enddo
      do N=1,10
	if(datac(N+20,NE,NEX) .NE. 'NO')then
	  close((N+30))
          close((N+40))
C          call splits(NHW,NHE,NVN,NVS,OUTS(N))
          call splits(NHW,NHE,NVN,NVS,outdir,OUTS(N), failure)
          if(failure==1)call endrun('Fail to process file '//
     2trim(outdir)//trim(OUTS(N))//' in '//trim(modfile),__LINE__)
          str='rm -f ' // OUTS(N)
          call system (str)
C          call splitp(NHW,NHE,NVN,NVS,nz,OUTP(N))
          call splitp(NHW,NHE,NVN,NVS,nz,outdir, OUTP(N), failure)
          if(failure==1)call endrun('Fail to process file '//
     2trim(outdir)//trim(OUTP(N))//' in '//trim(modfile),__LINE__)
          str = 'rm -f ' // OUTP(N)
          call system (str)
        endif
      enddo
      RETURN
      END
