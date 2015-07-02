
C+-----------------------------------------------------------------------+
C|                  TRANSFER FUNCTION FOR MADWEIGHT                      |
C|                                                                       |
C|     Author: Pierre Artoisenet (UCL-CP3)                               |
C|             Olivier Mattelaer (UCL-CP3)                               |
C+-----------------------------------------------------------------------+
C|     This file is generated automaticly by MADWEIGHT-TF_BUILDER        |
C+-----------------------------------------------------------------------+     


C+-----------------------------------------------------------------------+
C|    Transfer function for tf_PT_lepton
C+-----------------------------------------------------------------------+
      subroutine tf_PT_lepton(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'nb_tf.inc'
      include 'TF_param.inc'




        prov1=(tf_lepton_PT_7(curr_tf)+tf_lepton_PT_8(curr_tf)*dsqrt(pt(p))+
     &tf_lepton_PT_9(curr_tf)*pt(p))
        prov2=(tf_lepton_PT_10(curr_tf)+tf_lepton_PT_11(curr_tf)*dsqrt(pt(p))+
     &tf_lepton_PT_12(curr_tf)*pt(p))

        tf=(exp(-(pt(p)-pt(pexp)-prov1)**2/2d0/prov2**2))                !first gaussian
        tf=tf*((1d0/dsqrt(2d0*pi))/(prov2))            !normalisation



      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_PT_lepton
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_PT_lepton(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

       	  include 'nb_tf.inc'
      include 'TF_param.inc'





        width=(tf_lepton_PT_10(curr_tf)+tf_lepton_PT_11(curr_tf)*dsqrt(pt(pexp))+
     &tf_lepton_PT_12(curr_tf)*pt(pexp))




      width_PT_lepton= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_THETA_lepton
C+-----------------------------------------------------------------------+
      subroutine tf_THETA_lepton(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'nb_tf.inc'
      include 'TF_param.inc'



        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_THETA_lepton
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_THETA_lepton(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

       	  include 'nb_tf.inc'
      include 'TF_param.inc'




        width=0d0


      width_THETA_lepton= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_PHI_lepton
C+-----------------------------------------------------------------------+
      subroutine tf_PHI_lepton(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'nb_tf.inc'
      include 'TF_param.inc'



        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_PHI_lepton
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_PHI_lepton(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

       	  include 'nb_tf.inc'
      include 'TF_param.inc'




        width=0d0


      width_PHI_lepton= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_PT_jet
C+-----------------------------------------------------------------------+
      subroutine tf_PT_jet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'nb_tf.inc'
      include 'TF_param.inc'




        prov1=(tf_jet_PT_1(curr_tf)+tf_jet_PT_2(curr_tf)*dsqrt(pt(p))+
     &tf_jet_PT_3(curr_tf)*pt(p))
        prov2=(tf_jet_PT_4(curr_tf)+tf_jet_PT_5(curr_tf)*dsqrt(pt(p))+
     &tf_jet_PT_6(curr_tf)*pt(p))

        tf=(exp(-(pt(p)-pt(pexp)-prov1)**2/2d0/prov2**2))                !first gaussian
        tf=tf*((1d0/dsqrt(2d0*pi))/(prov2))            !normalisation 	



      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_PT_jet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_PT_jet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

       	  include 'nb_tf.inc'
      include 'TF_param.inc'





        width=(tf_jet_PT_4(curr_tf)+tf_jet_PT_5(curr_tf)*dsqrt(pt(pexp))+
     &tf_jet_PT_6(curr_tf)*pt(pexp))




      width_PT_jet= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_THETA_jet
C+-----------------------------------------------------------------------+
      subroutine tf_THETA_jet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'nb_tf.inc'
      include 'TF_param.inc'



        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_THETA_jet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_THETA_jet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

       	  include 'nb_tf.inc'
      include 'TF_param.inc'




        width=0d0


      width_THETA_jet= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_PHI_jet
C+-----------------------------------------------------------------------+
      subroutine tf_PHI_jet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'nb_tf.inc'
      include 'TF_param.inc'



        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_PHI_jet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_PHI_jet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

       	  include 'nb_tf.inc'
      include 'TF_param.inc'




        width=0d0


      width_PHI_jet= width

      return
      end



