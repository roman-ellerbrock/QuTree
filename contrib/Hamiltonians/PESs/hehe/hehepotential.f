      subroutine hehepotential(r,e)
      use potential_interface
      implicit none
      
      logical ret !dont use the retarded energy
      parameter (ret = .false.)

      real*8 r !distance of the helium atoms
      real*8 e !potential energy
      
      e = V_tot(ret,r) !calculates the totale potential
                       !energy

      end subroutine
