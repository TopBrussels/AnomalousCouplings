      subroutine GET_PERM_WEIGHT()
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      include 'permutation.inc'
      integer perm
      integer perm_id(nexternal-2)
      integer content(nexternal)
      integer i
      double precision WEIGHT
      double precision weight_perm_global, weight_perm_BW
      external weight_perm_global, weight_perm_BW


      do perm = 1, NPERM
         curr_perm = perm
         call get_perm(perm, perm_id)
         call assign_perm(perm_id)

        
        weight = weight_perm_global(perm,perm_id)
                  content(1) = 8
          content(2) = 7
          content(3) = 0
                weight = weight * weight_perm_BW(perm, perm_id, mdl_MW, mdl_WW, content, -24)
                
          content(1) = 8
          content(2) = 7
          content(3) = 6
          content(4) = 0
                weight = weight * weight_perm_BW(perm, perm_id, mdl_MT, mdl_WT, content, -6)
                
          content(1) = 5
          content(2) = 4
          content(3) = 0
                weight = weight * weight_perm_BW(perm, perm_id, mdl_MW, mdl_WW, content, 24)
                
          content(1) = 5
          content(2) = 4
          content(3) = 3
          content(4) = 0
                weight = weight * weight_perm_BW(perm, perm_id, mdl_MT, mdl_WT, content, 6)
                
        perm_value(perm, 1) = weight
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                    PERMUTATION MODULE                                   cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Need to contain TWO functions:                                                       cc                   
c           - weight_perm_global                                                           cc
c           - weight_perm_BW  (associate to a given propagator                             cc
c    The total wieght is the product of those weight. The permuation below a given         cc
c           relative cutoff are not consider for the integration.                          cc
c    The default function consist in a simple breit-wigner with an effective width dependingc
c            on the width of the transfer functions.                                        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      double precision function weight_perm_global(perm, perm_id)
c      implicit none
c      include 'nexternal.inc'
c      integer perm ! perm unique identifier
c      integer perm_id(nexternal-2) !order of the particles.
c      !      
c      !  DO nothing
c      !
c      weight_perm_global = 1d0
c
c      return 
c      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ANOTHER OPTION:
c     Forbid the permutation between the leptonic and hadronic b-jet.
c     b have position (tag_lhco or tag_init) 3 and 6 in MG (p p > t t~, (t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > j j)
c      --> This because counting only starts at number 3!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function weight_perm_global(perm, perm_id)
      implicit none
      include 'nexternal.inc'       
      integer perm ! perm unique identifier
      integer perm_id(nexternal-2) !order of the particles.
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco   ! get the id of the line of the LHCO information
c
c LHCO input
c
      include 'maxparticles.inc'
      integer tag_init(3:max_particles),type(max_particles),run_number,trigger
      double precision eta_init(max_particles),phi_init(max_particles),
     &pt_init(max_particles),j_mass(max_particles),ntrk(max_particles),
     &btag(max_particles),had_em(max_particles),dummy1(max_particles),
     &dummy2(max_particles)
      common/LHCO_input/eta_init,phi_init,pt_init,
     &j_mass,ntrk,btag,had_em,dummy1,dummy2,tag_init,type,run_number,
     &trigger
      integer lhco_b1, lhco_b2, lhco_bLept, lhco_bHadr 

      lhco_bLept = tag_init(3) !Always points to initial numbering
      lhco_bHadr = tag_init(6)
      lhco_b1 = tag_lhco(3)    !Numbering depends on the considered permutation
      lhco_b2 = tag_lhco(6)
      if (lhco_b1.ne.lhco_bLept.or.lhco_b2.ne.lhco_bHadr) then
	 !If initial numbering is different from numbering in considered permutation, the permutation should be discarded!
         weight_perm_global = 0d0
      else
         weight_perm_global = 1d0
      endif
      write(*,*) 'Global weight of perm',perm,'=',weight_perm_global
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function weight_perm_BW(perm, perm_id, Mass, Width, content, id)
      implicit none
      include 'nexternal.inc'
      include 'maxparticles.inc'
      include 'permutation.inc'
      integer perm_id(nexternal-2)
      integer perm
      integer content(nexternal)
      integer id
      double precision Mass, Width, InvMass2
      integer i,j,k
      double precision E_error(3:nexternal), eta_error(3:nexternal), phi_error(3:nexternal)
      double precision p_temp(0:3)
      double precision dot
      external dot
      double precision pexp(0:3,nexternal)
      integer tag_lhco(3:nexternal)
      common/to_pexp/pexp
      common/lhco_order/tag_lhco
      double precision c_point(NPERM,1:max_particles,3,2)
      common/ph_sp_init/c_point
c
c     User define parameter
c      
      double precision pi
      data pi/3.14159265359d0/ 
      double precision eff_width
c
c     Generic function
c
      do i=3, nexternal
         E_error(i) = c_point(perm,i,1,2)
         eta_error(i) =c_point(perm,i,2,2)
         phi_error(i) = c_point(perm,i,3,2)
      enddo
      do i=0,3
         p_temp(i) = 0d0
      enddo
      do i =1, nexternal
         j = content(i)
         if (j.gt.0) then
            !pexp information changes for the different permutations
            !So even if content keeps pointing to number 8, the
            !kinematic information gets changed, as desired!!
            p_temp(0) = p_temp(0) + pexp(0,j)
            p_temp(1) = p_temp(1) + pexp(1,j)
            p_temp(2) = p_temp(2) + pexp(2,j)
            p_temp(3) = p_temp(3) + pexp(3,j)
         else
            EXIT
         endif
      enddo
      InvMass2 = DOT(p_temp, p_temp) 
c
c     weighting by the Breit-Wigner
c
      eff_width = width**2
      do i = 1, nexternal
         if (content(i).gt.0) then
            if (E_error(content(i)).eq.-1d0)then
               !Don't want any weight calculation when considering a
               !neutrino, this because it doesn't have any x/y/z
               !kinematic component 
               weight_perm_bw = 1d0 ! no weight since neutrino!
               return
            endif
            eff_width = eff_width + E_error(content(i))**2
         else
            EXIT
         endif
      enddo
      eff_width = DSQRT(eff_width)
      weight_perm_bw = (InvMass2-Mass**2)**2+Mass**2*Width**2
      weight_perm_bw = Mass*Width/weight_perm_bw/pi
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
