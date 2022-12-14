c changed for usage with generic kinetic engery operator
c DS May 2018

c!-->  program to get potential energy for a given geometry after NN fitting
c!-->  global variables are declared in this module
c       module nnparam
c       implicit none
c       real*8,parameter::alpha=1.0d0,PI=3.1415926d0,radian=PI/180.0d0
c       integer,parameter::nbasis=1053,ndim=15,natom=6
c       real*8,parameter::vpesmin=-500.148240637577d0
c       integer ninput,noutput,nhid,nlayer,ifunc,nwe,nodemax
c       integer nscale
c       real*8,allocatable:: pesnn(:),pes(:),delta1(:)
c       integer, allocatable::nodes(:)
c       real*8, allocatable::weighta(:,:,:),biasa(:,:)
c       real*8, allocatable::pdela(:),pavga(:)
c
c       end module nnparam
! -------------------------------------------------------------------

      subroutine ch4clpot(v,q,m)
      implicit none
      integer    natom,ndimq,ndimx,i
      real*8     m(6)
      real*8     v,vmax
      parameter (natom=6)
      parameter (ndimx=natom*3)
      parameter (ndimq=ndimx-6)
      real*8     x(ndimx),q(ndimq),x0(ndimx)
      parameter (vmax=5.0d0/27.2114d0)

      call IntToCartStereo(q,x0,m)

      call ARRANGE(x0,x)

C     change from bohr to angstrom
      x = 0.5291772 * x

      call pipNN(x,v)
      v=v/27.2114d0

C     thresholds
      if (v.gt.vmax) v=vmax

C     just in case: if r < 0.5 A, then v -> vmax
c      if (q(7).lt.51) v=vmax

      return

      end
!  ***************************************************************************
!  *
!  *   Program    :  arrange
!  *   Function   :  re-arrange cartesian coordinates for ClCH4 potential    *
!  *                 the fourth H is reactive
!  *
!  ***************************************************************************
      subroutine    ARRANGE(cart,cart1)
      implicit      none
      real*8        cart(3,6),cart1(3,6)

      call DCOPY1D(3,cart(1,2),cart1(1,1))  ! H
      call DCOPY1D(3,cart(1,3),cart1(1,2))  ! H
      call DCOPY1D(3,cart(1,4),cart1(1,3))  ! H
      call DCOPY1D(3,cart(1,5),cart1(1,4))  ! reactive H
      call DCOPY1D(3,cart(1,1),cart1(1,5))  ! C
      call DCOPY1D(3,cart(1,6),cart1(1,6))  ! Cl
      
      return
      end

c ------------------------------------------------------------------------

      subroutine IntToCartStereo(q,x,m)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine transforms the internal coordinates       c
c to cartesian ones.                                        c
c Subroutine koordread must be called in main program       c
c Don't change anything below this line unless you          c
c know what you do. This part is written to work with the   c
c generic kinetic energy operator                           c
c DS / 4 May 2018                                           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer    dimch3
      parameter  (dimch3=12)

      real*8  q(12),x(18)
      real*8  r(3),xcm(3), snew, tnew
      real*8  radau(dimch3)
      real*8  mch3,mch4
      real*8     hmass, pi
      parameter (hmass=1836.109d0)
      parameter (pi=3.14159265359d0)

      real*8   rmch3
      real*8   dist

      integer i,j,k

      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      character*20                      coord(5:maxN)
      integer                           blockstart(9:20,5:maxN,5:maxN)

      real*8    m(6)
      real*8    work(4,4)

      real*8    redm(maxN), mabs(maxN), sredm(maxN)
      integer   start

      real*8 alpha(4)


c hardcode N
      N = 6
c hardcode stereographic
      do i=5, N
          coord(i)='stereo'
      enddo

c mass of CH3-group
      mch3=0.0d0
      do i=1,4
         mch3=mch3+m(i)
      end do
      mabs(4) = mch3

      do i = 5, N
        mabs(i) = 0.d0
        mabs(i) = mabs(i-1) + m(i)
        redm(i) = mabs(i-1)*m(i)/(mabs(i))
        sredm(i) = 1.d0/sqrt(redm(i))
      enddo

c hardcode rctrans (now called work)
      do i=1,4
        alpha(i) = sqrt(m(i)/mch3)
      enddo
      do i=1,4
        work(i,1)=alpha(i)
        work(1,i)=alpha(i)
      enddo
      do i=2,4
        do j=2,4
          work(j,i)=alpha(j)*alpha(i)/(alpha(1)+1.0d0)
        enddo
      enddo
      do i=2,4
        work(i,i)=work(i,i) - 1.0d0
      enddo

c order of q's:
c  1 rho_2        4 theta     7 R           10 r
c  2 phi_2        5 phi       8 theta_R     11 theta_r
c  3 theta_2      6 chi       9 phi_R       12 phi_r

c order of radau's
c  1 center of mass    4 x H_1     7 x H_2   10 x H_3
c  2 "                 5 y H_1     8 y H_2   11 y H_3
c  3 "                 6 z H_1     9 z H_2   12 z H_3

c order of cartesians
c 1-3 C (x,y,z)
c 4-6, 7-9, 10-12 methyl-H
c 13- 3*N coordinates for additional atoms

c METHYL

c First step: calculate r1,r2,r3 from rho_2, phi_2, theta_2

      r(1) = q(1)*cos(q(2))
      r(2) = q(1)*sin(q(2))*cos(q(3))
      r(3) = q(1)*sin(q(2))*sin(q(3))
c Second step: calculate mass weighted Radau's
c              from r1,r2,r3,theta,phi,chi

c     print*,"r's in Ruecktransformation:"
c     print'(F10.5)',(r(i),i=1,3)

      radau( 1)=0.0d0
      radau( 2)=0.0d0
      radau( 3)=0.0d0

      radau( 4)=r(1)*sin(q(4))
      radau( 5)=0.0d0
      radau( 6)=r(1)*cos(q(4))

      radau( 7)=r(2)*sin(q(4))*cos(q(5)+q(6))
      radau( 8)=r(2)*sin(q(4))*sin(q(5)+q(6))
      radau( 9)=r(2)*cos(q(4))

      radau(10)=r(3)*sin(q(4))*cos(-q(5)+q(6))
      radau(11)=r(3)*sin(q(4))*sin(-q(5)+q(6))
      radau(12)=r(3)*cos(q(4))

      !print*, r,q

c     print*,"Radau in Ruecktransformation"
      !do i=1,4
      !  print'(3F12.6)',(radau(j+(i-1)*3),j=1,3)
      !end do


c Third step: calculate cartesians for methyl group
c             from Radau's
      do i = 1,18
        x(i) = 0.0d0
      enddo
      do i=1,dimch3,3
         do j=1,dimch3,3
!            print*, 'koord', i, j
!            print*, i/3, j/3
!            print*, rctrans(i/3+1,j/3+1)
! hannes: hier ist work vorher rctrans gewesen:
!           x(i  )=x(i  )+rctrans(i/3+1,j/3+1)*radau(j  )
!           x(i+1)=x(i+1)+rctrans(i/3+1,j/3+1)*radau(j+1)
!           x(i+2)=x(i+2)+rctrans(i/3+1,j/3+1)*radau(j+2)
           x(i  )=x(i  )+work(i/3+1,j/3+1)*radau(j  )
           x(i+1)=x(i+1)+work(i/3+1,j/3+1)*radau(j+1)
           x(i+2)=x(i+2)+work(i/3+1,j/3+1)*radau(j+2)

         end do
      end do

c  unmassweight cartesians, correct sign
      do i=1,dimch3,3
         x(i  )=-x(i  )/sqrt(m(i/3+1))
         x(i+1)=-x(i+1)/sqrt(m(i/3+1))
         x(i+2)=-x(i+2)/sqrt(m(i/3+1))
      end do

c ADDITIONAL ATOMS 5 to N

      do i = 5, N

        start = 6+(i-5)*3
c     unmassweight r
      dist=q(start+1)*sredm(i)

      if (coord(i).eq.'stereo') then
c     inttocart for STEREOGRAPHIC COORDINATES (r, s, t)

        snew = 0.5d0*(q(start+2)*(-cos(2.d0/3.d0*q(6))
     .                      +sqrt(3.d0)*sin(2.d0/3.d0*q(6)))
     .               +q(start+3)*( sin(2.d0/3.d0*q(6))
     .                      +sqrt(3.d0)*cos(2.d0/3.d0*q(6))))


        tnew = 0.5d0*(q(start+3)*(-cos(2.d0/3.d0*q(6))
     .                      +sqrt(3.d0)*sin(2.d0/3.d0*q(6)))
     .               -q(start+2)*( sin(2.d0/3.d0*q(6))
     .                      +sqrt(3.d0)*cos(2.d0/3.d0*q(6))))

c     Center of mass for atoms 1 to i-1
        xcm(1)=0.0d0
        xcm(2)=0.0d0
        xcm(3)=0.0d0
        do k=1, i -1
                xcm(1)=xcm(1)+x(3*(k-1)+1)*m(k)
                xcm(2)=xcm(2)+x(3*(k-1)+2)*m(k)
                xcm(3)=xcm(3)+x(3*(k-1)+3)*m(k)
        end do
        do k=1,3
                xcm(k)=xcm(k)/(mabs(i-1))
        end do

c  calculate cartesians
        x((i-1)*3+1)=xcm(1)+2*dist*snew/(1.d0+snew**2+tnew**2)
        x((i-1)*3+2)=xcm(2)+2*dist*tnew/(1.d0+snew**2+tnew**2)
        x((i-1)*3+3)=xcm(3)
     .       +dist*(1.d0-snew**2-tnew**2)/(1.d0+snew**2+tnew**2)

      else if (coord(i).eq.'spherical') then
c     inttocart for SPHERICAL COORDINATES (r, theta, phi)

c    Center of mass for atoms 1 to i-1
        xcm(1)=0.0d0
        xcm(2)=0.0d0
        xcm(3)=0.0d0
        do k=1, i -1
                xcm(1)=xcm(1)+x(3*(k-1)+1)*m(k)
                xcm(2)=xcm(2)+x(3*(k-1)+2)*m(k)
                xcm(3)=xcm(3)+x(3*(k-1)+3)*m(k)
        end do
        do k=1,3
                xcm(k)=xcm(k)/(mabs(i-1))
        end do

        x((i-1)*3+1) = xcm(1) +
     .          dist*sin(q(start+2))*cos(q(start+3)+2.d0/3.d0*(q(6)-pi))
        x((i-1)*3+2) = xcm(2) +
     .          dist*sin(q(start+2))*sin(q(start+3)+2.d0/3.d0*(q(6)-pi))
        x((i-1)*3+3) = xcm(3) + dist*cos(q(start+2))

      else
                print*, 'ERROR: not supported: ', coord(i)
                print*, 'coordinates for atom:', i
                stop
      endif
      enddo

      end
c ------------------------------------------------------------------------
c=============================================================================
c     copy 1-dimen real array
c=============================================================================
       subroutine dcopy1d(ndim,xdat,ydat)
       implicit none
       integer :: ndim,idim
       real*8  :: xdat(ndim),ydat(ndim)
       do 1200 idim=1,ndim
          ydat(idim)=xdat(idim)
 1200  continue
       return
       end
