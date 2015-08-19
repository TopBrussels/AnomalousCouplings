ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'

      INCLUDE 'model_functions.inc'
      GC_17 = (MDL_CKM22*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_22 = (MDL_ELECC*MDL_COMPLEXI*MDL_GL)/(MDL_MW*MDL_SW*MDL_SQRT__
     $ 2)
      GC_23 = (MDL_ELECC*MDL_COMPLEXI*MDL_GR)/(MDL_MW*MDL_SW*MDL_SQRT__
     $ 2)
      GC_31 = -((MDL_ELECC*MDL_COMPLEXI*MDL_VL)/(MDL_SW*MDL_SQRT__2))
      GC_32 = -((MDL_ELECC*MDL_COMPLEXI*MDL_VR)/(MDL_SW*MDL_SQRT__2))
      GC_42 = -((MDL_ELECC*MDL_COMPLEXI*MDL_CONJG__GL)/(MDL_MW*MDL_SW
     $ *MDL_SQRT__2))
      GC_43 = -((MDL_ELECC*MDL_COMPLEXI*MDL_CONJG__GR)/(MDL_MW*MDL_SW
     $ *MDL_SQRT__2))
      GC_44 = (MDL_ELECC*MDL_COMPLEXI*MDL_CONJG__VL)/(MDL_SW*MDL_SQRT__
     $ 2)
      GC_45 = (MDL_ELECC*MDL_COMPLEXI*MDL_CONJG__VR)/(MDL_SW*MDL_SQRT__
     $ 2)
      END
