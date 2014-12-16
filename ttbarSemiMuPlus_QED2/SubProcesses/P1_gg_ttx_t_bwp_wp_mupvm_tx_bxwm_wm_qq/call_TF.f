
C+-----------------------------------------------------------------------+
C|                  TRANSFER FUNCTION FOR MADWEIGHT                      |
C|                                                                       |
C|     Author: Pierre Artoisenet (UCL-CP3)                               |
C|             Olivier Mattelaer (UCL-CP3)                               |
C+-----------------------------------------------------------------------+
C|     This file is generated automaticly by MADWEIGHT-TF_BUILDER        |
C+-----------------------------------------------------------------------+     


C+-----------------------------------------------------------------------+
C|     Subroutine: GET_CENTRAL_POINT                                     |
C|                                                                       |
C|     purpose: Define the central point of integration c_point(P,I,J,1) |
C|              and the width associated c_point(P,I,J,2)                |
C|                 P: Permutation                                        | 
C|                 I: MadGraph Position                                  |
C|                 J: 1,2,3 for THETA,PHI,RHO variable                   |
C+-----------------------------------------------------------------------+
      subroutine get_central_point
      implicit none
c                                                                        
c     parameter                                                          
c                                                                        
      include 'run.inc'
      include 'nexternal.inc'
      include 'permutation.inc'
      include 'TF_param.inc'

c                                                                        
c     global                                                             
c                                                                        
      double precision c_point(NPERM,1:max_particles,3,2)
      common/ph_sp_init/c_point
c                                                                        
c      double precision p_exp(0:3,nexternal)
c      common /to_pexp/p_exp
      double precision pexp_init(0:3,nexternal)  !impulsion in original configuration
      common/to_pexp_init/pexp_init
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c                                                                        
c     local                                                              
c                                                                       
      integer perm_id(nexternal-2)
      integer i,perm
      integer n_lhco



        external width_PT_bjet, width_THETA_bjet, width_PHI_bjet
        double precision width_PT_bjet, width_THETA_bjet, width_PHI_bjet
        external width_PT_muon, width_THETA_muon, width_PHI_muon
        double precision width_PT_muon, width_THETA_muon, width_PHI_muon
        external width_PT_nonbjet, width_THETA_nonbjet, width_PHI_nonbjet
        double precision width_PT_nonbjet, width_THETA_nonbjet, width_PHI_nonbjet
      do perm =1,NPERM
        call get_perm(perm, perm_id)

C+-----------------------------------------------------------------------+
C|     Start the definition                                              |
C+-----------------------------------------------------------------------+


        c_point(perm,1,1,2)=-1d0
        c_point(perm,1,2,2)=-1d0
        c_point(perm,1,3,2)=-1d0

        c_point(perm,2,1,2)=-1d0
        c_point(perm,2,2,2)=-1d0
        c_point(perm,2,3,2)=-1d0

        c_point(perm, 3,1,1)=theta(pexp_init(0,2+perm_id(1)))
        c_point(perm, 3,2,1)=phi(pexp_init(0,2+perm_id(1)))
        c_point(perm, 3,3,1)=rho(pexp_init(0,2+perm_id(1)))
        c_point(perm,3,1,2)=width_THETA_bjet(pexp_init(0,2+perm_id(1)),tag_lhco(3))
        c_point(perm,3,2,2)=width_PHI_bjet(pexp_init(0,2+perm_id(1)),tag_lhco(3))
        c_point(perm,3,3,2)=width_PT_bjet(pexp_init(0,2+perm_id(1)),tag_lhco(3))

        c_point(perm, 4,1,1)=theta(pexp_init(0,2+perm_id(2)))
        c_point(perm, 4,2,1)=phi(pexp_init(0,2+perm_id(2)))
        c_point(perm, 4,3,1)=rho(pexp_init(0,2+perm_id(2)))
        c_point(perm,4,1,2)=width_THETA_muon(pexp_init(0,2+perm_id(2)),tag_lhco(4))
        c_point(perm,4,2,2)=width_PHI_muon(pexp_init(0,2+perm_id(2)),tag_lhco(4))
        c_point(perm,4,3,2)=width_PT_muon(pexp_init(0,2+perm_id(2)),tag_lhco(4))

        c_point(perm,5,1,2)=-1d0
        c_point(perm,5,2,2)=-1d0
        c_point(perm,5,3,2)=-1d0

        c_point(perm, 6,1,1)=theta(pexp_init(0,2+perm_id(4)))
        c_point(perm, 6,2,1)=phi(pexp_init(0,2+perm_id(4)))
        c_point(perm, 6,3,1)=rho(pexp_init(0,2+perm_id(4)))
        c_point(perm,6,1,2)=width_THETA_bjet(pexp_init(0,2+perm_id(4)),tag_lhco(6))
        c_point(perm,6,2,2)=width_PHI_bjet(pexp_init(0,2+perm_id(4)),tag_lhco(6))
        c_point(perm,6,3,2)=width_PT_bjet(pexp_init(0,2+perm_id(4)),tag_lhco(6))

        c_point(perm, 7,1,1)=theta(pexp_init(0,2+perm_id(5)))
        c_point(perm, 7,2,1)=phi(pexp_init(0,2+perm_id(5)))
        c_point(perm, 7,3,1)=rho(pexp_init(0,2+perm_id(5)))
        c_point(perm,7,1,2)=width_THETA_nonbjet(pexp_init(0,2+perm_id(5)),tag_lhco(7))
        c_point(perm,7,2,2)=width_PHI_nonbjet(pexp_init(0,2+perm_id(5)),tag_lhco(7))
        c_point(perm,7,3,2)=width_PT_nonbjet(pexp_init(0,2+perm_id(5)),tag_lhco(7))

        c_point(perm, 8,1,1)=theta(pexp_init(0,2+perm_id(6)))
        c_point(perm, 8,2,1)=phi(pexp_init(0,2+perm_id(6)))
        c_point(perm, 8,3,1)=rho(pexp_init(0,2+perm_id(6)))
        c_point(perm,8,1,2)=width_THETA_nonbjet(pexp_init(0,2+perm_id(6)),tag_lhco(8))
        c_point(perm,8,2,2)=width_PHI_nonbjet(pexp_init(0,2+perm_id(6)),tag_lhco(8))
        c_point(perm,8,3,2)=width_PT_nonbjet(pexp_init(0,2+perm_id(6)),tag_lhco(8))


        enddo
        return
        end

C+-----------------------------------------------------------------------+
C|     Subroutine: transfer_fct(P,weight)                                |
C|                                                                       |
C|     purpose: returns the weight by the coeficient comming from the    |
C|              transfer_fct                                             |
C+-----------------------------------------------------------------------+
            subroutine transfer_fct(P,weight)
      implicit none
c                                                 
      integer    maxexternal
      parameter (maxexternal=15)
c                               
c     ARGUMENTS                 
c                               
      DOUBLE PRECISION P(0:3,maxexternal)
      DOUBLE PRECISION weight
c                                        
c     INCLUDE and COMMON                 
c                                        
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c                            
c                            
      double precision pexp(0:3,nexternal)
      integer tag_lhco(3:nexternal)
      common/to_pexp/pexp
      common/lhco_order/tag_lhco
C                               
C     SPECIAL FCT              
C                               
      DOUBLE PRECISION R2,DOT,ET,ETA,DJ,SumDot,PT
c                
c     LOCAL                                      
c                
      integer i,k
      integer n_lhco


        weight=1d0
        n_lhco=tag_lhco(3)
        call tf_PT_bjet(pexp(0,3),p(0,3),n_lhco,weight)
        call tf_THETA_bjet(pexp(0,3),p(0,3),n_lhco,weight)
        call tf_PHI_bjet(pexp(0,3),p(0,3),n_lhco,weight)

        n_lhco=tag_lhco(4)
        call tf_PT_muon(pexp(0,4),p(0,4),n_lhco,weight)
        call tf_THETA_muon(pexp(0,4),p(0,4),n_lhco,weight)
        call tf_PHI_muon(pexp(0,4),p(0,4),n_lhco,weight)

        n_lhco=tag_lhco(6)
        call tf_PT_bjet(pexp(0,6),p(0,6),n_lhco,weight)
        call tf_THETA_bjet(pexp(0,6),p(0,6),n_lhco,weight)
        call tf_PHI_bjet(pexp(0,6),p(0,6),n_lhco,weight)

        n_lhco=tag_lhco(7)
        call tf_PT_nonbjet(pexp(0,7),p(0,7),n_lhco,weight)
        call tf_THETA_nonbjet(pexp(0,7),p(0,7),n_lhco,weight)
        call tf_PHI_nonbjet(pexp(0,7),p(0,7),n_lhco,weight)

        n_lhco=tag_lhco(8)
        call tf_PT_nonbjet(pexp(0,8),p(0,8),n_lhco,weight)
        call tf_THETA_nonbjet(pexp(0,8),p(0,8),n_lhco,weight)
        call tf_PHI_nonbjet(pexp(0,8),p(0,8),n_lhco,weight)


        call check_nan(weight)
        return
        end

C+-----------------------------------------------------------------------+
C|     Subroutine: tf_E_for_part(MG_num)                                 |
C|                                                                       |
C|     purpose: returns the value of the transfert function (in energy)  |
C|              for the particles associated to the MG_number            |
C+-----------------------------------------------------------------------+
      double precision function  tf_PT_for_part(MG_num)
      implicit none
c
c     ARGUMENTS
c
      INTEGER MG_num
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer n_lhco
c


        if(MG_num.eq.1) then
        tf_PT_for_part=1d0
        return
        endif
        if(MG_num.eq.2) then
        tf_PT_for_part=1d0
        return
        endif
        if(MG_num.eq.3) then
        tf_PT_for_part=1d0
        n_lhco=tag_lhco(3)
        call tf_PT_bjet(pexp(0,3),momenta(0,3),n_lhco,tf_PT_for_part)

        return
        endif
        if(MG_num.eq.4) then
        tf_PT_for_part=1d0
        n_lhco=tag_lhco(4)
        call tf_PT_muon(pexp(0,4),momenta(0,4),n_lhco,tf_PT_for_part)

        return
        endif
        if(MG_num.eq.5) then
        tf_PT_for_part=1d0
        return
        endif
        if(MG_num.eq.6) then
        tf_PT_for_part=1d0
        n_lhco=tag_lhco(6)
        call tf_PT_bjet(pexp(0,6),momenta(0,6),n_lhco,tf_PT_for_part)

        return
        endif
        if(MG_num.eq.7) then
        tf_PT_for_part=1d0
        n_lhco=tag_lhco(7)
        call tf_PT_nonbjet(pexp(0,7),momenta(0,7),n_lhco,tf_PT_for_part)

        return
        endif
        if(MG_num.eq.8) then
        tf_PT_for_part=1d0
        n_lhco=tag_lhco(8)
        call tf_PT_nonbjet(pexp(0,8),momenta(0,8),n_lhco,tf_PT_for_part)

        return
        endif

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_PT_for_XX()                                   |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+



        double precision function tf_PT_for_1()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_PT_for_1=1d0

        return
        end



C+-----------------------------------------------------------------------+
C|          Subroutine: tf_PT_for_XX()                                   |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+



        double precision function tf_PT_for_2()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_PT_for_2=1d0

        return
        end



C+-----------------------------------------------------------------------+
C|          Subroutine: tf_PT_for_XX()                                   |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+



        double precision function tf_PT_for_3()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_PT_for_3=1d0
        n_lhco=tag_lhco(3)
        call tf_PT_bjet(pexp(0,3),momenta(0,3),n_lhco,tf_PT_for_3)

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_PT_for_XX()                                   |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+



        double precision function tf_PT_for_4()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_PT_for_4=1d0
        n_lhco=tag_lhco(4)
        call tf_PT_muon(pexp(0,4),momenta(0,4),n_lhco,tf_PT_for_4)

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_PT_for_XX()                                   |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+




        double precision function tf_PT_for_5()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_PT_for_5=1d0

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_PT_for_XX()                                   |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+




        double precision function tf_PT_for_6()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_PT_for_6=1d0
        n_lhco=tag_lhco(6)
        call tf_PT_bjet(pexp(0,6),momenta(0,6),n_lhco,tf_PT_for_6)

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_PT_for_XX()                                   |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+




        double precision function tf_PT_for_7()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_PT_for_7=1d0
        n_lhco=tag_lhco(7)
        call tf_PT_nonbjet(pexp(0,7),momenta(0,7),n_lhco,tf_PT_for_7)

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_PT_for_XX()                                   |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+



        double precision function tf_PT_for_8()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_PT_for_8=1d0
        n_lhco=tag_lhco(8)
        call tf_PT_nonbjet(pexp(0,8),momenta(0,8),n_lhco,tf_PT_for_8)

        return
        end



