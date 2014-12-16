ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE

      DOUBLE PRECISION PI
      PARAMETER  (PI=3.141592653589793D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'

      INCLUDE 'model_functions.inc'
      GC_13 = (EE*COMPLEXI)/(SW*SQRT__2)
      GC_22 = (ELECC*COMPLEXI*GL)/(MW*SW*SQRT__2)
      GC_23 = (ELECC*COMPLEXI*GR)/(MW*SW*SQRT__2)
      GC_31 = -((ELECC*COMPLEXI*VL)/(SW*SQRT__2))
      GC_32 = -((ELECC*COMPLEXI*VR)/(SW*SQRT__2))
      GC_38 = (EE*COMPLEXI*CONJG__CKM11)/(SW*SQRT__2)
      GC_39 = (EE*COMPLEXI*CONJG__CKM12)/(SW*SQRT__2)
      GC_40 = (EE*COMPLEXI*CONJG__CKM21)/(SW*SQRT__2)
      GC_41 = (EE*COMPLEXI*CONJG__CKM22)/(SW*SQRT__2)
      GC_42 = -((ELECC*COMPLEXI*CONJG__GL)/(MW*SW*SQRT__2))
      GC_43 = -((ELECC*COMPLEXI*CONJG__GR)/(MW*SW*SQRT__2))
      GC_44 = (ELECC*COMPLEXI*CONJG__VL)/(SW*SQRT__2)
      GC_45 = (ELECC*COMPLEXI*CONJG__VR)/(SW*SQRT__2)
      END
