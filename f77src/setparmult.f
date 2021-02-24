      SUBROUTINE setparmult
      include "parmult.h"

C     4 multipliers for SPOSC,
C     one for each of protein, carbonhydrate, cellulose, lignin
C      open(unit=99, file='ecosys_parameters.dat', action='read')
C      read(99,*) sposc_mult_1, sposc_mult_2, sposc_mult_3
C     2,sposc_mult_4, dckm_mult, dcki_mult, tsorp_mult, spoha_mult
C     3,spohc_mult, oqkm_mult, vmxo_mult
C      close(99)

C      vcmx_mul=1.0
C      rmplt_mul=1.0
C      vmxsoil_mul=1.0
C      okmsoil_mul=1.0
C      pscpl_mul=1.0
C      rexu_mul = 1.0
C     11 multipliers in all
      open(unit=111,file='ecosyspar.dat',status='old',action='read')
      read(111,*)vcmx_mul,rmplt_mul,vmxsoil_mul,okmsoil_mul
     2,rexu_mul,pscpl_mul
      close(111)

      print*,'vcmx_mul=',vcmx_mul
      print*,'rmplt_mul=',rmplt_mul
      print*,'vmxsoil_mul=',vmxsoil_mul
      print*,'okmsoil_mul=',okmsoil_mul
      print*,'rexu_mul=',rexu_mul
      print*,'pscpl_mul=',pscpl_mul
      END
