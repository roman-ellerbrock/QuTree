! subroutine to determin PES matrix for dvrms using diabatic model from genetic
! written by W. Eisfeld
! Mar. 6. 2014

      subroutine vpot(v, mu, q, p, nm, ns)
!      subroutine vpot(Ed, mu, q, p, nm, ns)
      implicit none

!--------------params.incl----------------------------------------------------
      include 'params.incl'

!--------------states.incl----------------------------------------------------
      include 'states.incl'

!--------------vwzprec.incl---------------------------------------------------
      include 'vwzprec.incl'

!--------------common.incl---------------------------------------------------
      include 'common.incl'



!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      integer i,j,k,ii

!dimensions: modes, states:
      integer nm, ns

!switch for type of coordinate transformation:
      integer symm

!input coordinates and internally used coordiates
      real*8 q(nm), lq(6), qcart(12), qint(7)

!potential parameters:
      real*8 p(*)

!energies:
      real*8 Ed(nstat,nstat), Ea(nstat)      !diabatic matrix and adiabatic energies
      real*8 v(ns,ns)

!dipole moment:
      real*8 mu(3)

!geometry parameters
      real*8 reno, reoo, rad, refen
      parameter(reno=2.344419d0,reoo=4.0606528d0,rad=0.017453293d0)
      parameter(refen=-0.279586867D+03)    !jump parameters
!      parameter(refen=-0.279589924D+03)    !no-jump parameters
      real*8 ev
      parameter (ev=27.2114d0)
!      parameter (ev=1.d0)


!tuning parameters for potential
      real*8 dum, rho, thr, fac, alpha
      parameter (thr=12.d0, fac=1.0d0, alpha=0.05d0)

!switches for various stuff
      logical dbg, skip, angle, valence, cut


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


!.....INITIALIZE PARAMETERS
!      call pinit(p)
      angle = .true.          !use bending angles for O-N-O rather than O-O distances


!.....generate cartesian coordinates:
      do i=1,12
        qcart(i)=qref(i)             !initialize with reference geometry
      enddo
      do i=1,nm
        lq(i)=q(i)/dl(i)         !remove factor for dimensionless coordinate
	do j=1,12
	  qcart(j)=qcart(j)+lq(i)*mode(j,i)
	enddo
      enddo
!      write(6,'(6f12.6)') q
!      write(6,'(6f12.6)') dl
!      write(6,'(6f12.6)') lq
!      write(6,*)
!      write(6,'(3f12.6)') qcart
!      write(6,*)


!.....TRANSFORM FROM CARTESIAN COORDINATES TO INTERNALLY USED SET
      call internal(qcart,qint)             !this returns N-O distances in bohr and projected O-N-O angles
!      qint(7)=lq(6)/rad                           !this is supposed to be in degrees!!!
      cut=.false.
      if (abs(qint(1)).gt.1.d0) cut=.true.
      if (abs(qint(2)).gt.1.d0) cut=.true.
      if (abs(qint(3)).gt.1.d0) cut=.true.


!      call adia(Ed,mu,qint,p,pst,skip)

      call adia(Ea,mu,qint,p,pst,skip)

      v(1,1)=Ea(1)-refen

!make sure PES is bound
      if (cut) v(1,1)=0.25d0

      if (v(1,1).gt.0.25d0) v(1,1)=0.25d0
      do i=1,nstat
        Ed(i,i)=Ed(i,i) - refen
      enddo

      return

!---------------------------------------------------------------------
!..   ensure potential is bound
      rho=0.d0
      do j=1,nm
        rho=rho+abs(q(j))
      enddo
!      rho=rho-abs(q(2))
      if (rho.gt.thr) then
        dum=fac*(exp(alpha*(rho-thr))-1.d0)
!	write(6,'(10f12.4)') q, rho, dum
	v(1,1)=v(1,1)+dum
      endif

100   format(2f12.6,2f16.7)


      end

