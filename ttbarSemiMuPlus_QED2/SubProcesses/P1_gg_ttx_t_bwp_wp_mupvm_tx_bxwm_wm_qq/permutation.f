C**************************************************************************************
      subroutine assign_perm(perm_id)
C
C
      include 'nexternal.inc'
      include "phasespace.inc"
      integer perm_id(nexternal-2)     !permutation of 1,2,...,nexternal-2
C
      integer i,j

      double precision pexp_init(0:3,nexternal)  !impulsion in original configuration
      common/to_pexp_init/pexp_init
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp

      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
      integer tag_init(3:max_particles),type(max_particles),run_number,trigger
      double precision eta_init(max_particles),phi_init(max_particles),
     &pt_init(max_particles),j_mass(max_particles),ntrk(max_particles),
     &btag(max_particles),had_em(max_particles),dummy1(max_particles),
     &dummy2(max_particles)
      common/LHCO_input/eta_init,phi_init,pt_init,
     &j_mass,ntrk,btag,had_em,dummy1,dummy2,tag_init,type,run_number,
     &trigger

      integer matching_type_part(3:max_particles) !modif/link 
     &!between our order by type for permutation
      integer inv_matching_type_part(3:max_particles)
      common/madgraph_order_type/matching_type_part,
     & inv_matching_type_part

      do j=3,nexternal
         do i=0,3
            pexp(i,j)=pexp_init(i, 2+perm_id(j-2))
         enddo
         tag_lhco(j)=tag_init(2+perm_id(j-2))
      enddo

      end

      subroutine get_perm(nb, perm)
      implicit none
      integer i,j
      include 'nexternal.inc'
      INTEGER    NB
      INTEGER    PERM(NEXTERNAL-2)
      include 'permutation.inc'
      INTEGER PERMS(NPERM, NEXTERNAL-2)
      DATA (PERMS(1,I),I=1,6) /1,2,3,4,5,6/
      DATA (PERMS(2,I),I=1,6) /1,2,3,4,6,5/
      DATA (PERMS(3,I),I=1,6) /1,2,3,5,4,6/
      DATA (PERMS(4,I),I=1,6) /1,2,3,5,6,4/
      DATA (PERMS(5,I),I=1,6) /1,2,3,6,4,5/
      DATA (PERMS(6,I),I=1,6) /1,2,3,6,5,4/
      DATA (PERMS(7,I),I=1,6) /4,2,3,1,5,6/
      DATA (PERMS(8,I),I=1,6) /4,2,3,1,6,5/
      DATA (PERMS(9,I),I=1,6) /4,2,3,5,1,6/
      DATA (PERMS(10,I),I=1,6) /4,2,3,5,6,1/
      DATA (PERMS(11,I),I=1,6) /4,2,3,6,1,5/
      DATA (PERMS(12,I),I=1,6) /4,2,3,6,5,1/
      DATA (PERMS(13,I),I=1,6) /5,2,3,1,4,6/
      DATA (PERMS(14,I),I=1,6) /5,2,3,1,6,4/
      DATA (PERMS(15,I),I=1,6) /5,2,3,4,1,6/
      DATA (PERMS(16,I),I=1,6) /5,2,3,4,6,1/
      DATA (PERMS(17,I),I=1,6) /5,2,3,6,1,4/
      DATA (PERMS(18,I),I=1,6) /5,2,3,6,4,1/
      DATA (PERMS(19,I),I=1,6) /6,2,3,1,4,5/
      DATA (PERMS(20,I),I=1,6) /6,2,3,1,5,4/
      DATA (PERMS(21,I),I=1,6) /6,2,3,4,1,5/
      DATA (PERMS(22,I),I=1,6) /6,2,3,4,5,1/
      DATA (PERMS(23,I),I=1,6) /6,2,3,5,1,4/
      DATA (PERMS(24,I),I=1,6) /6,2,3,5,4,1/
        do i=1, NEXTERNAL-2
            perm(i) = PERMS(nb, i)
        enddo
        return
        end

