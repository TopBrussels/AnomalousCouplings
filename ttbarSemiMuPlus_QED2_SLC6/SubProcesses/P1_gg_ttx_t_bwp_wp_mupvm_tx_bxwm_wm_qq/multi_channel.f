C+-----------------------------------------------------------------------+
C|                MULTI-CHANNEL WEIGHT FOR MADWEIGHT                     |
C|                                                                       |
C|     Author: Pierre Artoisenet (UCL-CP3)                               |
C|             Olivier Mattelaer (UCL-CP3)                               |
C+-----------------------------------------------------------------------+
C|     This file is generated automaticly by MADWEIGHT-ANALYZER          |
C+-----------------------------------------------------------------------+

       double precision function multi_channel_weight(config)
C+-----------------------------------------------------------------------+
C|     routine returnings the multi channel weight linked to the         |
C|       change of variable 'config'                                     |
C+-----------------------------------------------------------------------+
       	implicit none

       	integer config
       	include 'coupl.inc'
           include 'd_choices.inc'
c    double precision prov1,prov2,prov3,prov4,prov5,prov6,prov7,prov8,prov9
c    double precision prov10,prov11,prov12,prov13,prov14,prov15,prov16
c    double precision prov17,prov18,prov19,prov20,prov21,prov22,prov23 
           double precision num,den
           double precision zero,one
           parameter (zero=0d0)
           parameter (one=1d0)
        double precision local_1
        double precision  Breit_Wigner_for_part
        external  Breit_Wigner_for_part
        double precision local_2

        if (config.eq.1) then
        local_1 = Breit_Wigner_for_part( -4, MT, WT)
        local_2 = Breit_Wigner_for_part( -3, MW, WW)

        num = 1d0 * 1d0 * local_2
        den = 0d0 + 1d0 * local_2 + 1d0 * local_1
        multi_channel_weight = num/den

        elseif (config.eq.2) then
        local_1 = Breit_Wigner_for_part( -4, MT, WT)
        local_2 = Breit_Wigner_for_part( -3, MW, WW)

        num = 1d0 * 1d0 * local_1
        den = 0d0 + 1d0 * local_2 + 1d0 * local_1
        multi_channel_weight = num/den
       	   endif
       	   return
       	   end

