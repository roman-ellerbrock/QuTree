      subroutine potentialch5p (v,q)
      
      implicit none
      
      real*8 vbed,vmax
      parameter (vbed=12.d0/219.47463068d0)
      parameter (vmax=5.d0/27.212d0)
      real*8    q(*), v
      real*8 en
      
      integer N
      parameter (N=6)
      real*8 cart(N,3)
      integer i, j
      logical*4 vm
      real*8 mh
      real*8 ar
      
c     Massenparameter, Handbook of Atomic, Molecular,
c     and Optical Physics, Gordon W. F. Drake
      parameter (ar=1822.88848d0)

      real*8 r1, r2, theta
      real*8 x(0:N-1,3), rad(0:N-1,3)
      real*8 Cadj (0:N-1,0:N-1)
      real*8 m(0:N-1), alpha(0:N-1), mg
      real*8 dist(N-1,N-1)
      
      r1=q(1)
      r2=q(2)
c     Adjustment for c++-Version
c      q(3)=acos(q(3))
      theta=q(3)

c     12C-Isotopenmasse      
      
      m(0)=12.0d0*ar

      mg=m(0)
      
c     Isotopenmasse, Handbook of Atomic, Molecular,
c     and Optical Physics, Gordon W. F. Drake
      mh=1.007825047d0*ar      
       
      do i=1,N-1
         m(i)=mh
         mg=mg+m(i)
      enddo

      do i=0,N-1
         alpha(i)=dsqrt(m(i)/mg)
      enddo

c     Umrechnung interne -> kart. Radaus

c     Massenschwerpunkt -> (0,0,0)T
      
      rad(0,1)=0.0d0
      rad(0,2)=0.0d0
      rad(0,3)=0.0d0

      rad(1,1)=0.0d0
      rad(1,2)=0.0d0
      rad(1,3)=r1

      rad(2,1)=r2*dsin(theta)
      rad(2,2)=0.0d0
      rad(2,3)=r2*dcos(theta)
      
      do i=3,N-1
         rad(i,1)=2.0d0*q(3*(i-2)+1)*q(3*(i-2)+2)
     .            /(1+q(3*(i-2)+2)**2+q(3*(i-2)+3)**2)

         rad(i,2)=2.0d0*q(3*(i-2)+1)*q(3*(i-2)+3)
     .            /(1+q(3*(i-2)+2)**2+q(3*(i-2)+3)**2)

c     negatives Vorzeichen z-Komponente bei Projektion vom Nordpol       
       
       rad(i,3)=-1.0d0*q(3*(i-2)+1)*(1-q(3*(i-2)+2)**2-q(3*(i-2)+3)**2)
     .            /(1+q(3*(i-2)+2)**2+q(3*(i-2)+3)**2)
      enddo

c     Transformationsmatrix aufbauen
            
      do i=0,N-1
         Cadj(0,i)=alpha(i)
         Cadj(i,0)=alpha(i)
      enddo

      do i=1,N-1
         do j=1, N-1
           Cadj(i,j)=alpha(i)*alpha(j)/(alpha(0)+1.0d0)
         enddo
      enddo

c     Diagonalelemente korrigieren

      do i=1,N-1
         Cadj(i,i)=Cadj(i,i)-1.0d0              
      enddo
      
c     Transformation der Radau-Vektoren

      do i=0,N-1
         x(i,1)=0.0d0
         x(i,2)=0.0d0
         x(i,3)=0.0d0
      enddo

      do i=0,N-1
         do j=0,N-1
            x(i,1)=x(i,1)+Cadj(i,j)*rad(j,1)
            x(i,2)=x(i,2)+Cadj(i,j)*rad(j,2)
            x(i,3)=x(i,3)+Cadj(i,j)*rad(j,3)
         enddo
      enddo

c     Massengewichtung rueckgaengig, Vorzeichenwechsel

      do i=0,N-1
         x(i,1)=-1.0d0*x(i,1)/dsqrt(m(i))
         x(i,2)=-1.0d0*x(i,2)/dsqrt(m(i))
         x(i,3)=-1.0d0*x(i,3)/dsqrt(m(i))
      enddo
 
c     Koordinaten Potentialroutine      
      
      do i=0,N-1
         cart(i+1,1)=x(i,1)
         cart(i+1,2)=x(i,2)
         cart(i+1,3)=x(i,3)
      enddo

c     test for unphysical interproton distances
      call checkdistch5p (cart,N,vm)
     
      call getpot(cart,en)
      
      v=en
      v=v+40.652764870344082086480940942556d0
      
      if (vm) then
c         if (v.lt.vbed) v=2.d0*vbed-v
         if (v.lt.vbed) v=vbed
      endif
       
      if (v.gt.vmax) v=vmax
      
c      v=v-8000.d0/219474.63068d0


c      if (vm) then
c         do i=1,6
c            write(98,*) cart(i,1)
c            write(98,*) cart(i,2)
c            write(98,*) cart(i,3)
c         enddo
c         write(98,*) v*219474.63068d0
c         write(98,*)
c      endif

      return 
      end
      
      subroutine checkdistch5p (x,N,vm)

      implicit none
      
      real*8 hh2,hh3,hh4,hh5
      parameter (hh2=1.12d0,hh3=1.7d0,hh4=1.9d0,hh5=3.d0)
      integer nh
      parameter (nh=10)

      integer N,i,j,l
      logical*4 vm
      real*8 x(N,3),dist(nh),h

c calculate interproton distances
      l=0
      do i=2,N
         do j=i+1,N
            l=l+1
            dist(l)=(x(j,1)-x(i,1))**2
     .             +(x(j,2)-x(i,2))**2
     .             +(x(j,3)-x(i,3))**2        
            dist(l)=dsqrt(dist(l))
         enddo
      enddo

c     sort interproton distances
      do i=1,nh
         do j=i+1,nh
            if (dist(j).lt.dist(i)) then
               h=dist(i)
               dist(i)=dist(j)
               dist(j)=h
            endif
         enddo
      enddo

c     check interproton distances
      vm=.true.
c      print*, 'x='
c      do i=1,6
c      print*, x(i,:)
c      enddo
c      print*, 'dist='
c      print*, dist(:)
c      read(5,*)
      if (dist(2).le.hh2) return 
      if (dist(3).le.hh3) return 
      if (dist(4).le.hh4) return 
      if (dist(5).le.hh5) return 
      vm=.false.

      return
      end

!***********************************************************************
!CH5+ Potential Energy Surface, version 7a, 2005-10-18
!   by Zhong Jin, Bastiaan J. Braams, Joel M. Bowman
!
!reference:
!   Zhong Jin, Bastiaan J. Braams, Joel M. Bowman
!   Journal of Physical Chemistry A  accepted(2005)
!
!notes:
!
!   x(6,3) is the Cartesian coordinates for six atoms in order of
!   C H H H H H. (in bohr)
!
!   r is defined as the distance between the center of H2 and Carbon 
!   in a.u.
!
!   Potential energy epot is returned in a.u.
!
!   CCSD(T)/aug-cc-pVTZ based
!
!   When compiling, "-r8" option is required for accurate potential
!   energy calculation
!   
!All rights reserved. Contact bowman@euch4e.chem.emory.edu for
!details or any latest updates.
! 
************************************************************************
      SUBROUTINE getpot(x,epot)
!
      implicit none
!
      integer, parameter :: natoms = 6
      integer :: i, j, k, l, m

      double precision :: sw, a, b, theta
      double precision :: vfit, epot
      double precision, dimension(natoms,3) :: x
      double precision, dimension(3,natoms) :: xn
      double precision, dimension(4,3) :: xch3
      double precision, dimension(2,3) :: xh2
      double precision, parameter :: de = -40.57833477d0
      double precision :: epotch3, epoth2, elr
      double precision :: r, rh2, rp
      double precision, parameter :: rmin = 11.d0, rmax = 15.d0

      do i = 1, natoms
         do j = 1, 3
            xn(j,i) = x(i,j)
         end do
      end do

      call calcr(xn,xch3,xh2,r,rh2,theta)

c      write(*,'a2,2x,f9.2') 'r=', r
      a = 0.d0
      b = 0.d0

      if (r<=rmax) then
         call getfit(xn,vfit)
      end if

      rp = (r - rmin)/(rmax - rmin)
!      write(*,*) 'rp=', rp
      if (r<rmin) then
         epot = vfit
      else
         call getlr(r,rh2,theta,elr)
         call getch3pot(xch3,epotch3)
         call geth2pot(xh2,epoth2)

c         write(*,'a8,f9.3')'epotch3=', epotch3
c         write(*,'a7,1x,f9.3')'epoth2=', epoth2
c         write(*,'a10,e20.10')'longrange=',elr
         if (r>rmax) then
            epot = elr + epoth2 + epotch3 + de
         else
            epot = (1.d0-sw(a,b,rp))*vfit + sw(a,b,rp)*
     *             (elr + epoth2 + epotch3 + de)
         end if
      end if

      END SUBROUTINE getpot

      double precision FUNCTION sw(a,b,x)
      implicit none
      double precision :: a, b, x
      sw = 0.5d0*a*x**2 + (10.d0 + b/2.d0 - 3.d0*a/2.d0)*x**3 + 
     &    (3.d0*a/2.d0 - b - 15.d0)*x**4 + (6.d0 + b/2.d0 - a/2.d0)*x**5
      END FUNCTION sw

******************************************************************************
!
!.... Subroutine to obtain the long-range interaction
!
!.... 02/14/05   Zhong      SUBROUTINE created
!.... 02/22/05   Zhong      Modified for the long range effect 

      SUBROUTINE getlr(r,rh2,theta,elr)

      implicit none

      double precision :: r, rh2, theta
      double precision :: alpha, alpha0, alpha90
      double precision :: u3, u4, p2, q
      double precision :: epot, elr
      double precision :: pi

      pi = 4.d0 * datan(1.d0)
      call cspline(rh2,alpha0,alpha90,q) 

c      write(*,'a4,1x,f9.5')'rh2=',rh2
      theta = theta / 180.d0 * pi
      alpha = 1.d0/3.d0*(alpha0 + 2*alpha90)
      p2 = 0.5d0*(3*(cos(theta))**2 - 1.d0)
      u3 = q*p2/r**3
      u4 = -(0.5d0*alpha + 1.d0/3.d0*(alpha0 - alpha90)*p2)/r**4
      elr = u3 + u4 

      END SUBROUTINE getlr     

      SUBROUTINE geth2pot(xyz1,v)
!...  program is to calculate the potential energy of H2
!...  Based on pes6n4
!...
      implicit none
      integer, parameter :: natoms = 6
      integer, parameter :: natoms1 = 2
      integer :: i, j

      double precision, dimension(5:6,3) :: xyz1
      double precision, dimension(natoms,3) :: xyz
      double precision, dimension(3,natoms) :: xyzt, f
      double precision :: v, ch3pot, ch3potab, h20
      double precision, parameter :: ang = 219474.6d0

      character(len=1), dimension(natoms) :: xname

      do i = 5, 6
         do j = 1, 3
            xyz(i,j) = xyz1(i,j)
         end do
      end do

      xyz(1,1) = 0.d0
      xyz(1,2) = 0.d0
      xyz(1,3) = 0.d0
      xyz(2,1) = 0.d0
      xyz(2,2) = 2.0657297d0
      xyz(2,3) = 0.d0
      xyz(3,1) = 1.7889744d0
      xyz(3,2) = -1.0328649d0
      xyz(3,3) = 0.d0
      xyz(4,1) = -1.7889744d0
      xyz(4,2) = -1.0328649d0
      xyz(4,3) = 0.d0

      do i = 1, natoms
         do j = 1, 3
            xyzt(j,i) = xyz(i,j)
         end do
      end do

      call getfit(xyzt,v)

      ch3pot = -39.4046592d0 + 0.0006415d0/2.d0

!      ch3potab = -39.40570072d0
      h20 = -1.1745554d0 + 0.0006728d0/2.d0

      v = v - ch3pot - h20
      if ((v*ang.le.5.d0).and.(v*ang.ge.-5.d0)) then
         v = 0.d0
      end if
      END SUBROUTINE geth2pot

***********************************************************************
      SUBROUTINE getfit(xn,f0)
!
      implicit none
!

!    m = numbers of coefficients for the fit

      integer, parameter :: m = 2239

      integer, parameter :: mr = 32

      double precision :: dc0 = 0.30547300513
      double precision :: dc1 = 0.19476072488
      double precision :: dw0 = 0.09846745640
      double precision :: dw1 = 0.12149234254

!   coef = the coefficients of the fit

      double precision,dimension(0:m+2*mr-1) :: coef 
     
      include "fit.2.inc"
!
      double precision, dimension(0:2,0:5) :: xn
      double precision, dimension(0:2,0:5) :: xn1
      double precision, dimension(0:2,0:5) :: gf0
      double precision, dimension(0:5,0:5) :: d0, r0
      double precision, dimension(0:2) :: rvec
      double precision, dimension(0:m+2*mr-1) :: vec,vec0,vec1
      double precision :: f0 
      double precision :: t0,t1
      double precision, parameter :: dd = 1.0d-6
!
      integer :: i, j, l, k
!
!     -----------------------------------------------------------------
! The following lines are the code that evaluates my fitted function.
! The xnuc are the nuclear coordinates, as found in your table.
! Routine getd0 evaluates my internal coordinates d0(0:5,0:5), which are
! functions of the interparticle distances.  I've defined d0 both above
! and below the diagonal, but it is symmetric (d0(i,j)=d0(j,i)), and d0
! vanishes on the diagonal.  Then routine vec evaluates the (m) basis
! functions at the configuration described by d0, and finally f0 is the
! fitted potential.  The basis functions are polynomials of the entries
! in d0, obeying all the expected symmetries.
!

      call getd0 (xn, dc0, dc1, dw0, dw1, d0, r0)
      call getvec (m, d0, vec(0:m-1))
      call getrvec (3, r0, rvec)
      do l = 0, mr-1
         do k = 0, 1
            vec(m+2*l+k) = rvec(k+1)*vec(l)
         enddo
      enddo

      f0 = dot_product(coef,vec)

!.... Now the code is serving for PE calculation, gradient
!...  is ignored for decreasing computational cost

!       do i = 0, 2
!          do j = 0, 5
!             xn1 = xn ; xn1(i,j) = xn1(i,j)-dd
!             call getd0 (xn1, dc0, dc1, dw0, dw1, d0, r0)
!             call getvec (m, d0, vec(0:m-1))
!             call getrvec (3, r0, rvec)
!             do l = 0, mr-1
!                do k = 0, 1
!                vec(m+2*l+k) = rvec(k+1)*vec(l)
!             enddo
!          enddo
!           t0 = dot_product(coef,vec)

!           xn1 = xn ; xn1(i,j) = xn1(i,j)+dd
!           call getd0 (xn1, dc0, dc1, dw0, dw1, d0, r0)
!           call getvec (m, d0, vec(0:m-1))
!           call getrvec (3, r0, rvec)
!          do l = 0, mr-1
!             do k = 0, 1
!                vec(m+2*l+k) = rvec(k+1)*vec(l)
!             enddo
!          enddo
!           t1 = dot_product(coef,vec)

!
! Note that gf0 will be the negative gradient
!
!           gf0(i,j) = -(t1-t0)/(2*dd)
!         enddo
!       enddo
!
      END SUBROUTINE getfit

!********************************************************************
      SUBROUTINE getd0 (xn, dc0, dc1, dw0, dw1, d0, r0)

      implicit none

      double precision, dimension(0:2,0:5) :: xn
      double precision :: dc0, dc1, dw0, dw1
      double precision, dimension(0:5,0:5) :: d0, r0
      double precision, parameter :: p = 2.0d0
      integer :: i, j
      double precision :: t0

! Note: CH5+.  Nuclei 0..5; i=0 for the C, i in 1..5 for the H.

      d0(0,0) = 0.d0
      r0(0,0) = 0.d0

      do j = 1, 5
         t0 = sqrt((xn(0,j)-xn(0,0))**2+(xn(1,j)-xn(1,0))**2+ 
     $     (xn(2,j)-xn(2,0))**2)
!         d0(0,j) = (log(t0)-dc0)/dw0
!         d0(0,j) = (dexp(-t0/3.d0)-dc0)/dw0
         d0(0,j) = (dexp(-t0/p)-dc0)/dw0
         d0(j,0) = d0(0,j)
         r0(0,j) = t0
         r0(j,0) = t0
      end do

      do i = 1, 5
         d0(i,i) = 0.d0
         r0(i,i) = 0.d0
         do j = i+1, 5
            t0 = sqrt((xn(0,j)-xn(0,i))**2+(xn(1,j)-xn(1,i))**2+ 
     $        (xn(2,j)-xn(2,i))**2)
!           d0(i,j) = (log(t0)-dc1)/dw1
!           d0(i,j) = (dexp(-t0/p)-dc1)/dw1
           d0(i,j) = (dexp(-t0/p)-dc1)/dw1
            d0(j,i) = d0(i,j)
            r0(i,j) = t0
            r0(j,i) = t0
         end do
      end do
      END SUBROUTINE getd0

!**********************************************************************
      SUBROUTINE getvec (m, d, vec)

      implicit none

      integer :: m

      double precision, dimension(0:5,0:5) :: d
      double precision, dimension(0:m-1) :: vec

      integer :: j0, j1, j2, j3, j4, j5, j, k, k0
      double precision,dimension(0:1) :: x
      double precision,dimension(0:3) :: y
      double precision,dimension(0:9) :: z
      double precision,dimension(0:21) :: u
      double precision,dimension(0:41) :: v
      double precision,dimension(0:59) :: w
      double precision,dimension(0:5,0:5) :: d2,d3,d4,d5
      double precision :: t0
      double precision :: her2, her3, her4, her5, her6, her7

      her2(t0) = (4*t0**2-2)/sqrt(dble(8*2))
      her3(t0) = (8*t0**2-12)*t0/sqrt(dble(16*6))
      her4(t0) = ((16*t0**2-48)*t0**2+12)/sqrt(dble(32*24))
      her5(t0) = ((32*t0**2-160)*t0**2+120)*t0/sqrt(dble(64*120))
      her6(t0) = (((64*t0**2-480)*t0**2+720)*t0**2-120)/sqrt(dble(
     $    128*720))
      her7(t0) = (((128*t0**2-1344)*t0**2+3360)*t0**2-1680)*t0/ 
     $    sqrt(dble(256*5040))
!     ------------------------------------------------------------------
!     Test for compatibility
      if (.not.(m.eq.1.or.m.eq.3.or.m.eq.10.or.m.eq.32.or. 
     $   m.eq.101.or.m.eq.299.or.m.eq.849.or.m.eq.2239)) then
      read(5,*)
       stop 'getvec - wrong dimension'
      endif

!     Computation
      x = 0 ; y = 0 ; z = 0 ; u = 0 ; v = 0 ; w = 0
      do j0 = 0, 5
       do j1 = j0+1, 5
        d2(j0,j1) = her2(d(j0,j1))
        d2(j1,j0) = d2(j0,j1)
        d3(j0,j1) = her3(d(j0,j1))
        d3(j1,j0) = d3(j0,j1)
        d4(j0,j1) = her4(d(j0,j1))
        d4(j1,j0) = d4(j0,j1)
        d5(j0,j1) = her5(d(j0,j1))
        d5(j1,j0) = d5(j0,j1)
       enddo
      enddo

      do j1 = 1, 5
       x(0) = x(0)+d(0,j1)/5
       y(0) = y(0)+d2(0,j1)/5
       z(0) = z(0)+d3(0,j1)/5
       u(0) = u(0)+d4(0,j1)/5
       v(0) = v(0)+d5(0,j1)/5
c       w() = w()+her6(d(0,j1))/5
       do j2 = 1, 5
       if (j2.ne.j1) then
        y(1) = y(1)+d(0,j1)*d(j1,j2)/20
        z(1) = z(1)+d2(0,j1)*d(j1,j2)/20
        z(2) = z(2)+d(0,j1)*d2(j1,j2)/20
        z(3) = z(3)+d(0,j1)*d(0,j2)*d(j1,j2)/20
        u(1) = u(1)+d3(0,j1)*d(j1,j2)/20
        u(2) = u(2)+d2(0,j1)*d2(j1,j2)/20
        u(3) = u(3)+d(0,j1)*d3(j1,j2)/20
        u(4) = u(4)+d2(0,j1)*d(0,j2)*d(j1,j2)/20
        u(5) = u(5)+d(0,j1)*d(0,j2)*d2(j1,j2)/20
        v(1) = v(1)+d4(0,j1)*d(j1,j2)/20
        v(2) = v(2)+d3(0,j1)*d2(j1,j2)/20
        v(3) = v(3)+d2(0,j1)*d3(j1,j2)/20
        v(4) = v(4)+d(0,j1)*d4(j1,j2)/20
        v(5) = v(5)+d3(0,j1)*d(0,j2)*d(j1,j2)/20
        v(6) = v(6)+d2(0,j1)*d(0,j2)*d2(j1,j2)/20
        v(7) = v(7)+d(0,j1)*d(0,j2)*d3(j1,j2)/20
        w(0) = w(0)+d4(0,j1)*d2(j1,j2)/20
        w(1) = w(1)+d3(0,j1)*d3(j1,j2)/20
        w(2) = w(2)+d2(0,j1)*d4(j1,j2)/20
        w(3) = w(3)+d(0,j1)*d5(j1,j2)/20
        w(4) = w(4)+d3(0,j1)*d2(0,j2)*d(j1,j2)/20
        w(5) = w(5)+d3(0,j1)*d(0,j2)*d2(j1,j2)/20
        w(6) = w(6)+d2(0,j1)*d(0,j2)*d3(j1,j2)/20
        w(7) = w(7)+d(0,j1)*d(0,j2)*d4(j1,j2)/20
c        w() = w()+d5(0,j1)*d(j1,j2)/20
c        w() = w()+d4(0,j1)*d(0,j2)*d(j1,j2)/20
c        w() = w()+d2(0,j1)*d2(0,j2)*d2(j1,j2)/20
        do j3 = 1, 5
        if (j3.ne.j2.and.j3.ne.j1) then
         z(4) = z(4)+d(0,j1)*d(j1,j2)*d(j2,j3)/60
         z(5) = z(5)+d(0,j1)*d(j1,j2)*d(j1,j3)/60
         u(6) = u(6)+d(0,j1)*d2(j1,j2)*d(j2,j3)/60
         u(7) = u(7)+d(0,j1)*d(j1,j2)*d2(j2,j3)/60
         u(8) = u(8)+d(0,j1)*d2(j1,j2)*d(j1,j3)/60
         u(9) = u(9)+d(0,j1)*d(j1,j2)*d(j1,j3)*d(j2,j3)/60
         u(10) = u(10)+d2(0,j1)*d(j1,j2)*d(j2,j3)/60
         u(11) = u(11)+d2(0,j1)*d(j1,j2)*d(j1,j3)/60
         u(12) = u(12)+d(0,j1)*d(0,j2)*d(j1,j3)*d(j2,j3)/60
         v(8) = v(8)+d3(0,j1)*d(j1,j2)*d(j2,j3)/60
         v(9) = v(9)+d2(0,j1)*d2(j1,j2)*d(j2,j3)/60
         v(10) = v(10)+d2(0,j1)*d(j1,j2)*d2(j2,j3)/60
         v(11) = v(11)+d(0,j1)*d3(j1,j2)*d(j2,j3)/60
         v(12) = v(12)+d(0,j1)*d2(j1,j2)*d2(j2,j3)/60
         v(13) = v(13)+d(0,j1)*d(j1,j2)*d3(j2,j3)/60
         v(14) = v(14)+d3(0,j1)*d(j1,j2)*d(j1,j3)/60
         v(15) = v(15)+d2(0,j1)*d2(j1,j2)*d(j1,j3)/60
         v(16) = v(16)+d(0,j1)*d3(j1,j2)*d(j1,j3)/60
         v(17) = v(17)+d(0,j1)*d2(j1,j2)*d2(j1,j3)/60
         v(18) = v(18)+d2(0,j1)*d(j1,j2)*d(j1,j3)*d(j2,j3)/60
         v(19) = v(19)+d(0,j1)*d2(j1,j2)*d(j1,j3)*d(j2,j3)/60
         v(20) = v(20)+d(0,j1)*d(j1,j2)*d(j1,j3)*d2(j2,j3)/60
         v(21) = v(21)+d2(0,j1)*d(0,j2)*d(j1,j3)*d(j2,j3)/60
         v(22) = v(22)+d(0,j1)*d(0,j2)*d2(j1,j3)*d(j2,j3)/60
         v(23) = v(23)+d(0,j1)*d(0,j2)*d(j1,j3)*d(j2,j3)*d(j1,j2)/60
         v(24) = v(24)+d(0,j1)*d(0,j2)*d(j1,j2)*d2(j1,j3)/60
c         w() = w()+d4(0,j1)*d(j1,j2)*d(j2,j3)/60
         w(8) = w(8)+d(0,j1)*d(j1,j2)*d4(j2,j3)/60
         w(9) = w(9)+d(0,j1)*d4(j1,j2)*d(j2,j3)/60
         w(10) = w(10)+d3(0,j1)*d2(j1,j2)*d(j2,j3)/60
         w(11) = w(11)+d(0,j1)*d2(j1,j2)*d3(j2,j3)/60
         w(12) = w(12)+d2(0,j1)*d3(j1,j2)*d(j2,j3)/60
         w(13) = w(13)+d(0,j1)*d3(j1,j2)*d2(j2,j3)/60
         w(14) = w(14)+d3(0,j1)*d(j1,j2)*d2(j2,j3)/60
         w(15) = w(15)+d2(0,j1)*d(j1,j2)*d3(j2,j3)/60
         w(16) = w(16)+d2(0,j1)*d2(j1,j2)*d2(j2,j3)/60
c         w() = w()+d4(0,j1)*d(j1,j2)*d(j1,j3)/60
         w(17) = w(17)+d(0,j1)*d4(j1,j2)*d(j1,j3)/60
         w(18) = w(18)+d3(0,j1)*d2(j1,j2)*d(j1,j3)/60
         w(19) = w(19)+d2(0,j1)*d3(j1,j2)*d(j1,j3)/60
         w(20) = w(20)+d(0,j1)*d3(j1,j2)*d2(j1,j3)/60
         w(21) = w(21)+d2(0,j1)*d2(j1,j2)*d2(j1,j3)/60
         w(22) = w(22)+d3(0,j1)*d(0,j2)*d(j1,j3)*d(j2,j3)/60
         w(23) = w(23)+d(0,j1)*d(0,j2)*d3(j1,j3)*d(j2,j3)/60
         w(44) = w(44)+d2(0,j1)*d2(0,j2)*d(j1,j3)*d(j2,j3)/60
         w(24) = w(24)+d2(0,j1)*d(0,j2)*d2(j1,j3)*d(j2,j3)/60
         w(25) = w(25)+d(0,j1)*d(0,j2)*d2(j1,j3)*d2(j2,j3)/60
         w(26) = w(26)+d2(0,j1)*d(0,j2)*d(j1,j3)*d2(j2,j3)/60
         w(45) = w(45)+d3(0,j1)*d(0,j2)*d(j1,j2)*d(j2,j3)/60
         w(27) = w(27)+d(0,j1)*d(j1,j2)*d(j1,j3)*d3(j2,j3)/60
         w(28) = w(28)+d(0,j1)*d(0,j2)*d3(j1,j2)*d(j2,j3)/60
         w(46) = w(46)+d3(0,j1)*d(0,j2)*d(j1,j2)*d(j1,j3)/60
c         w() = w()+d3(0,j1)*d(0,j2)*d(0,j3)*d(j1,j2)/60
         w(29) = w(29)+d(0,j1)*d3(j1,j2)*d(j1,j3)*d(j2,j3)/60
c         w() = w()+d(0,j1)*d(0,j2)*d(j1,j2)*d3(j2,j3)/60
         w(30) = w(30)+d3(0,j1)*d(j1,j2)*d(j1,j3)*d(j2,j3)/60
         w(31) = w(31)+d2(0,j1)*d(0,j2)*d2(j1,j2)*d(j2,j3)/60
c         w() = w()+d2(0,j1)*d2(0,j2)*d(j1,j2)*d(j2,j3)/60
c         w() = w()+d2(0,j1)*d(0,j2)*d(0,j3)*d2(j1,j2)/60
         w(32) = w(32)+d(0,j1)*d2(j1,j2)*d(j1,j3)*d2(j2,j3)/60
c         w() = w()+d2(0,j1)*d(0,j2)*d(j1,j2)*d2(j2,j3)/60
         w(33) = w(33)+d2(0,j1)*d(j1,j2)*d(j1,j3)*d2(j2,j3)/60
c         w() = w()+d2(0,j1)*d(0,j2)*d2(j1,j2)*d(j1,j3)/60
c         w() = w()+d2(0,j1)*d2(0,j2)*d(0,j3)*d(j1,j2)/60
         w(34) = w(34)+d(0,j1)*d2(j1,j2)*d2(j1,j3)*d(j2,j3)/60
c         w() = w()+d(0,j1)*d(0,j2)*d2(j1,j2)*d2(j2,j3)/60
c         w() = w()+d(0,j1)*d2(0,j2)*d(j1,j2)*d2(j2,j3)/60
         w(35) = w(35)+d2(0,j1)*d(j1,j2)*d2(j1,j3)*d(j2,j3)/60
         w(36) = w(36)+d2(0,j1)*d(0,j2)*d(0,j3)*d(j1,j2)*d(j2,j3)/60
c         w() = w()+d2(0,j1)*d(0,j2)*d(j1,j2)*d(j1,j3)*d(j2,j3)/60
c         w() = w()+d(0,j1)*d(0,j2)*d(0,j3)*d2(j1,j2)*d(j2,j3)/60
         w(37) = w(37)+d(0,j1)*d(0,j2)*d(j1,j2)*d2(j1,j3)*d(j2,j3)/60
         w(38) = w(38)+d2(0,j1)*d(0,j2)*d(0,j3)*d(j1,j2)*d(j1,j3)/60
         w(39) = w(39)+d(0,j1)*d(0,j2)*d2(j1,j2)*d(j1,j3)*d(j2,j3)/60
c       w() = w()+d(0,j1)*d(0,j2)*d(0,j3)*d(j1,j2)*d(j1,j3)*d(j2,j3)/60
         do j4 = 1, 5
         if (j4.ne.j3.and.j4.ne.j2.and.j4.ne.j1) then
          u(13) = u(13)+d(0,j1)*d(j1,j2)*d(j1,j3)*d(j3,j4)/120
          u(14) = u(14)+d(0,j1)*d(j1,j2)*d(j1,j3)*d(j1,j4)/120
          v(25) = v(25)+d(0,j1)*d2(j1,j2)*d(j2,j3)*d(j3,j4)/120
          v(26) = v(26)+d2(0,j1)*d(j1,j2)*d(j1,j3)*d(j3,j4)/120
          v(27) = v(27)+d(0,j1)*d2(j1,j2)*d(j1,j3)*d(j3,j4)/120
          v(28) = v(28)+d(0,j1)*d(j1,j2)*d2(j1,j3)*d(j3,j4)/120
          v(29) = v(29)+d(0,j1)*d(j1,j2)*d(j1,j3)*d2(j3,j4)/120
          v(30) = v(30)+d2(0,j1)*d(j1,j2)*d(j1,j3)*d(j1,j4)/120
          v(31) = v(31)+d(0,j1)*d2(j1,j2)*d(j1,j3)*d(j1,j4)/120
          w(40) = w(40)+d(0,j1)*d3(j1,j2)*d(j1,j3)*d(j1,j4)/120
          w(41) = w(41)+d2(0,j1)*d(j1,j2)*d2(j2,j3)*d(j3,j4)/120
          w(42) = w(42)+d(0,j1)*d(j1,j2)*d2(j2,j3)*d2(j3,j4)/120
          w(43) = w(43)+d(0,j1)*d2(j1,j2)*d2(j1,j3)*d(j1,j4)/120
          w(44) = w(44)+d(0,j1)*d(j1,j2)*d(j2,j3)*d3(j3,j4)/120
          w(45) = w(45)+d(0,j1)*d2(j1,j2)*d(j2,j3)*d2(j3,j4)/120
          w(46) = w(46)+d(0,j1)*d(j1,j2)*d(j1,j3)*d2(j2,j3)*d(j1,j4)/120
c          w() = w()+d3(0,j1)*d(j1,j2)*d(j2,j3)*d(j3,j4)/120
c          w() = w()+d2(0,j1)*d2(j1,j2)*d(j2,j3)*d(j3,j4)/120
c          w() = w()+d3(0,j1)*d(j1,j2)*d(j1,j3)*d(j1,j4)/120
c          w() = w()+d2(0,j1)*d2(j1,j2)*d(j1,j3)*d(j1,j4)/120
c          w() = w()+d(0,j1)*d(j1,j2)*d(j1,j3)*d(j2,j3)*d2(j1,j4)/120
c          w() = w()+d2(0,j1)*d(j1,j2)*d(j1,j3)*d(j2,j3)*d(j1,j4)/120
c          w() = w()+d(0,j1)*d2(j1,j2)*d(j1,j3)*d(j2,j3)*d(j1,j4)/120
         endif
         enddo
        endif
        enddo
       endif
       enddo
      enddo
      do j0 = 1, 5
       do j1 = 1, 5
       if (j1.ne.j0) then
        x(1) = x(1)+d(j0,j1)/4
        y(2) = y(2)+d2(j0,j1)/4
        z(6) = z(6)+d3(j0,j1)/4
        u(15) = u(15)+d4(j0,j1)/4
        v(32) = v(32)+d5(j0,j1)/4
c        w() = w()+her6(d(j0,j1))/4
        do j2 = 1, 5
        if (j2.ne.j1.and.j2.ne.j0) then
         y(3) = y(3)+d(j0,j1)*d(j1,j2)/12
         z(7) = z(7)+d2(j0,j1)*d(j1,j2)/12
         z(8) = z(8)+d(j0,j1)*d(j0,j2)*d(j1,j2)/12
         u(16) = u(16)+d3(j0,j1)*d(j1,j2)/12
         u(17) = u(17)+d2(j0,j1)*d2(j1,j2)/12
         u(18) = u(18)+d2(j0,j1)*d(j0,j2)*d(j1,j2)/12
         v(33) = v(33)+d4(j0,j1)*d(j1,j2)/12
         v(34) = v(34)+d3(j0,j1)*d2(j1,j2)/12
         v(35) = v(35)+d3(j0,j1)*d(j0,j2)*d(j1,j2)/12
         v(36) = v(36)+d2(j0,j1)*d2(j0,j2)*d(j1,j2)/12
         w(47) = w(47)+d4(j0,j1)*d2(j1,j2)/12
         w(48) = w(48)+d3(j0,j1)*d3(j1,j2)/12
         w(49) = w(49)+d2(j0,j1)*d2(j0,j2)*d2(j1,j2)/12
         w(50) = w(50)+d4(j0,j1)*d(j0,j2)*d(j1,j2)/12
c         w() = w()+d3(j0,j1)*d2(j0,j2)*d(j1,j2)/12
c         w() = w()+d5(j0,j1)*d(j1,j2)/12
         do j3 = 1, 5
         if (j3.ne.j2.and.j3.ne.j1.and.j3.ne.j0) then
          z(9) = z(9)+d(j0,j1)*d(j0,j2)*d(j0,j3)/24
          u(19) = u(19)+d2(j0,j1)*d(j0,j2)*d(j0,j3)/24
          u(20) = u(20)+d(j0,j1)*d(j1,j2)*d(j1,j3)*d(j2,j3)/24
          u(21) = u(21)+d(j0,j1)*d(j0,j2)*d(j1,j3)*d(j2,j3)/24
          v(37) = v(37)+d3(j0,j1)*d(j0,j2)*d(j0,j3)/24
          v(38) = v(38)+d2(j0,j1)*d2(j0,j2)*d(j0,j3)/24
          v(39) = v(39)+d2(j0,j1)*d(j1,j2)*d(j1,j3)*d(j2,j3)/24
          v(40) = v(40)+d(j0,j1)*d2(j1,j2)*d(j1,j3)*d(j2,j3)/24
          v(41) = v(41)+d(j0,j1)*d(j1,j2)*d(j1,j3)*d2(j2,j3)/24
          w(51) = w(51)+d4(j0,j1)*d(j1,j2)*d(j2,j3)/24
          w(52) = w(52)+d(j0,j1)*d4(j1,j2)*d(j2,j3)/24
          w(53) = w(53)+d3(j0,j1)*d2(j1,j2)*d(j2,j3)/24
          w(54) = w(54)+d3(j0,j1)*d(j1,j2)*d2(j2,j3)/24
          w(55) = w(55)+d2(j0,j1)*d3(j1,j2)*d(j2,j3)/24
          w(56) = w(56)+d2(j0,j1)*d2(j1,j2)*d2(j2,j3)/24
          w(57) = w(57)+d4(j0,j1)*d(j0,j2)*d(j0,j3)/24
          w(58) = w(58)+d2(j0,j1)*d2(j0,j2)*d2(j0,j3)/24
          w(59) = w(59)+d3(j0,j1)*d2(j0,j2)*d(j0,j3)/24
         endif
         enddo
        endif
        enddo
       endif
       enddo
      enddo
      vec(0) = 1
      if (3.le.m) then
       vec(1) = x(0)
       vec(2) = x(1)
      endif
      if (10.le.m) then
       vec(3) = her2(x(0))
       vec(4) = x(0)*x(1)
       vec(5) = her2(x(1))
       do k = 0, 3
        vec(6+k) = y(k)
       enddo
      endif
!! third order terms
      if (32.le.m) then
       k0 = 10
       vec(k0) = her3(x(0))
       vec(k0+1) = her2(x(0))*x(1)
       vec(k0+2) = x(0)*her2(x(1))
       vec(k0+3) = x(0)*y(0)
       vec(k0+4) = x(0)*y(1)
       vec(k0+5) = x(0)*y(2)
       vec(k0+6) = x(0)*y(3)
       vec(k0+7) = her3(x(1))
       vec(k0+8) = x(1)*y(0)
       vec(k0+9) = x(1)*y(1)
       vec(k0+10) = x(1)*y(2)
       vec(k0+11) = x(1)*y(3)
       k0 = 22
       do k = 0, 9
        vec(k0+k) = z(k)
       enddo
      endif
!! fourth order terms
      if (101.le.m) then
       k0 = 32
       vec(k0) = her4(x(0))
       vec(k0+1) = her3(x(0))*x(1)
       vec(k0+2) = her2(x(0))*her2(x(1))
       do k = 0, 3
        vec(k0+3+k) = her2(x(0))*y(k)
       enddo
       vec(k0+7) = x(0)*her3(x(1))
       do k = 0, 3
        vec(k0+8+k) = x(0)*x(1)*y(k)
       enddo
       do k = 0, 9
        vec(k0+12+k) = x(0)*z(k)
       enddo
       vec(k0+22) = her4(x(1))
       do k = 0, 3
        vec(k0+23+k) = her2(x(1))*y(k)
       enddo
       do k = 0, 9
        vec(k0+27+k) = x(1)*z(k)
       enddo
       k0 = k0+37
       do j = 0, 3
        do k = j, 3
         vec(k0) = y(j)*y(k)
         k0 = k0+1
        enddo
       enddo
       if (k0.ne.79) then
       read(5,*)
        stop 'getvec: counting error'
       endif
       do k = 0, 21
        vec(k0+k) = u(k)
       enddo
      endif
!! fifth order terms
      if (299.le.m) then
       k0 = 101
       vec(k0) = her5(x(0))
       vec(k0+1) = her4(x(0))*x(1)
       vec(k0+2) = her3(x(0))*her2(x(1))
       do k = 0, 3
        vec(k0+3+k) = her3(x(0))*y(k)
       enddo
       vec(k0+7) = her2(x(0))*her3(x(1))
       do k = 0, 3
        vec(k0+8+k) = her2(x(0))*x(1)*y(k)
       enddo
       do k = 0, 9
        vec(k0+12+k) = her2(x(0))*z(k)
       enddo
       vec(k0+22) = x(0)*her4(x(1))
       do k = 0, 3
        vec(k0+23+k) = x(0)*her2(x(1))*y(k)
       enddo
       do k = 0, 9
        vec(k0+27+k) = x(0)*x(1)*z(k)
       enddo
       k0 = k0+37
       do j = 0, 3
        do k = j, 3
         vec(k0) = x(0)*y(j)*y(k)
         k0 = k0+1
        enddo
       enddo
       if (k0.ne.148) then
       read(5,*)
        stop 'getvec: counting error'
       endif
       do k = 0, 21
        vec(k0+k) = x(0)*u(k)
       enddo
       vec(k0+22) = her5(x(1))
       do k = 0, 3
        vec(k0+23+k) = her3(x(1))*y(k)
       enddo
       do k = 0, 9
        vec(k0+27+k) = her2(x(1))*z(k)
       enddo
       k0 = k0+37
       do j = 0, 3
        do k = j, 3
         vec(k0) = x(1)*y(j)*y(k)
         k0 = k0+1
        enddo
       enddo
       if (k0.ne.195) then
       read(5,*)
        stop 'getvec: counting error'
       endif
       do k = 0, 21
        vec(k0+k) = x(1)*u(k)
       enddo
       do j = 0, 3
        do k = 0, 9
         vec(k0+22+10*j+k) = y(j)*z(k)
        enddo
       enddo
       k0 = 257
       do k = 0, 41
        vec(k0+k) = v(k)
       enddo
      endif
!! sixth order terms
      if (849.le.m) then
       k0 = 299
       do k = 0, 197
        vec(k0+k) = x(0)*vec(101+k)
       enddo
       do k = 0, 128
        vec(k0+198+k) = x(1)*vec(170+k)
       enddo
       do k = 0, 31
        vec(k0+327+k) = y(0)*vec(69+k)
       enddo
       do k = 0, 27
        vec(k0+359+k) = y(1)*vec(73+k)
       enddo
       do k = 0, 24
        vec(k0+387+k) = y(2)*vec(76+k)
       enddo
       do k = 0, 22
        vec(k0+412+k) = y(3)*vec(78+k)
       enddo
       k0 = k0+435
       do j = 0, 9
        do k = j, 9
         vec(k0) = z(j)*z(k)
         k0 = k0+1
        enddo
       enddo
       if (k0.ne.789) then
       read(5,*)
        stop 'getvec: counting error'
       endif
       do k = 0, 59
        vec(k0+k) = w(k)
       enddo
      endif
!! seventh order terms (incomplete)
      if (2239.le.m) then
       k0 = 849
       do k = 0, 549 ! 849-299-1
        vec(k0+k) = x(0)*vec(299+k)
       enddo
       do k = 0, 351 ! 849-497-1
        vec(k0+550+k) = x(1)*vec(497+k)
       enddo
       do k = 0, 81 ! 299-217-1
        vec(k0+902+k) = y(0)*vec(217+k)
       enddo
       do k = 0, 71
        vec(k0+984+k) = y(1)*vec(227+k)
       enddo
       do k = 0, 61
        vec(k0+1056+k) = y(2)*vec(237+k)
       enddo
       do k = 0, 51
        vec(k0+1118+k) = y(3)*vec(247+k)
       enddo
       k0 = k0+1170
       do j = 0, 9
        do k = 0, 21
         vec(k0+22*j+k) = z(j)*u(k)
        enddo
       enddo
      endif
      return
      END SUBROUTINE getvec
!********************************************************
      SUBROUTINE getrvec (m, r, vec)

      implicit none

      integer m

      double precision, dimension(0:5,0:5) :: r, r1
      double precision, dimension(0:m-1) :: vec
      integer :: j0, j1, j2, i, j, k
      double precision, dimension(0:1) :: x
      double precision, dimension(0:3) :: y
!     ------------------------------------------------------------------
!     Test for compatibility
      if (.not.(m.eq.1.or.m.eq.3.or.m.eq.10)) then
       read(5,*)
       stop 'getrvec - wrong dimension'
      endif
!     Computation
      x = 0 ; y = 0
      do i = 0, 5
       do j = 0, 5
        if (i.eq.j) then
         r1(i,j) = 0
        else
         r1(i,j) = exp(-r(i,j))/r(i,j)
        endif
       enddo
      enddo
      do j1 = 1, 5
       x(0) = x(0)+r1(0,j1)/5
       y(0) = y(0)+r1(0,j1)**2/5
       do j2 = 1, 5
       if (j2.ne.j1) then
        y(1) = y(1)+r1(0,j1)*r1(j1,j2)/20
       endif
       enddo
      enddo
      do j0 = 1, 5
       do j1 = 1, 5
       if (j1.ne.j0) then
        x(1) = x(1)+r1(j0,j1)/4
        y(2) = y(2)+r1(j0,j1)**2/4
        do j2 = 1, 5
        if (j2.ne.j1.and.j2.ne.j0) then
         y(3) = y(3)+r1(j0,j1)*r1(j1,j2)/12
        endif
        enddo
       endif
       enddo
      enddo
      vec(0) = 1
      if (3.le.m) then
       vec(1) = x(0)
       vec(2) = x(1)
      endif
      if (10.le.m) then
       vec(3) = x(0)**2
       vec(4) = x(0)*x(1)
       vec(5) = x(1)**2
       do k = 0, 3
        vec(6+k) = y(k)
       enddo
      endif
      END SUBROUTINE getrvec     
      
      SUBROUTINE getch3pot(xyz1,v)
!...  program is to calculate the potential energy of CH3+
!...  Based on pes6n4
!...
      implicit none
      integer, parameter :: natoms = 6
      integer, parameter :: natoms1 = 4
      integer :: i, j

      double precision, dimension(natoms1,3) :: xyz1
      double precision, dimension(natoms,3) :: xyz
      double precision, dimension(3,natoms) :: xyzt, f
      double precision :: v, h2pot, h2potab, ch30
      double precision, parameter :: ang = 219474.6d0

      character(len=1), dimension(natoms) :: xname

      do i = 1, natoms1
         do j = 1, 3
            xyz(i,j) = xyz1(i,j)
         end do
      end do

      xyz(5,1) = 0.7035d0
      xyz(5,2) = 0.d0
      xyz(5,3) = 14.0750143d0
      xyz(6,1) = -0.7035d0
      xyz(6,2) = 0.d0
      xyz(6,3) = 14.0750143d0

      do i = 1, natoms
         do j = 1, 3
            xyzt(j,i) = xyz(i,j)
         end do
      end do

      call getfit(xyzt,v)

      h2pot = -1.1745554d0 + 0.0006728d0/2.d0
!      h2potab = -1.17263405d0
      ch30 = -39.4046592d0 + 0.0006415d0/2.d0
!      ch30 = -39.40570072d0

      v = v - h2pot - ch30

      if ((v*ang.le.5.d0).and.(v*ang.ge.-5.d0)) then
         v = 0.d0
      end if
      END SUBROUTINE getch3pot

*****************************************************************************
!...  subroutine for cubic spline
!...  driver for routine splint, which calls spline
!...  x ....       H-H bond distance
!     y0 ...       alpha parallel
!     y90 ...      alpha perpendicular
!     yq ...       quadrupole      
*****************************************************************************
      SUBROUTINE cspline(x,y0,y90,yq)

      implicit none

      integer, parameter :: NP = 21, NP2 = 257
      integer :: i, nfunc
      double precision :: x, y0, yp01, yp0n, y90, yp901, yp90n
      double precision :: ypq1, ypqn, yq
      double precision, dimension(NP) :: xa, ya0, ya90, y2s0, y2s90
      double precision, dimension(NP2) :: xa1, ya1, yq2

      open(unit=70,file='../data/CH5+/alpha_1967.dat',status='old')
      open(unit=71,file='../data/CH5+/q.dat',status='old')

      read(70,*)
      do i = 1, NP
         read(70,*) xa(i), ya0(i), ya90(i)
      end do
         yp01 = 0.d0
         yp0n = 0.d0
         yp901 = 0.d0
         yp90n = 0.d0

C     call SPLINE to get second derivatives
         call spline(xa,ya0,NP,yp01,yp0n,y2s0)
         call spline(xa,ya90,NP,yp901,yp90n,y2s90)

C     call SPLINT for interpolations
         call splint(xa,ya0,y2s0,NP,x,y0)
         call splint(xa,ya90,y2s90,NP,x,y90)

      do i = 1, NP2
         read(71,*) xa1(i), ya1(i)
      end do
         ypq1 = 0.d0
         ypqn = 0.d0

         call spline(xa1,ya1,NP2,ypq1,ypqn,yq2)
         call splint(xa1,ya1,yq2,NP2,x,yq)

      close(70)
      close(71)         
      END SUBROUTINE cspline
***********************************************************************
!...  Cubic spline code
!...  Original cubic spline code from Numerial Recipe has been modified
!...  to adapt double precision
!...  02/23/05    Zhong
!
***********************************************************************
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit none

      integer, parameter :: NMAX=500
      integer :: i, k, n

      double precision :: yp1,ypn
      double precision :: p, qn, sig, un
      double precision, dimension(n) :: x, y, y2
      double precision, dimension(NMAX) :: u
 
      if (yp1.gt..99e30) then
         y2(1)=0.5d0
         u(1)=0.5d0
      else
         y2(1)=-.5d0
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i = 2, n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/
     $ (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      if (ypn.gt..99e30) then
         qn=0.d0
         un=0.d0
      else
         qn=.5d0
         un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

      do k = n-1, 1, -1
         y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      END SUBROUTINE spline
****************************************************************************
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      implicit none
      integer :: n
      double precision :: x, y
      double precision, dimension(n) :: xa, y2a, ya

      integer :: k, khi, klo
      double precision :: a, b, h

      klo = 1
      khi = n

1     if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if (xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
        goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     $2)/6.d0
      END SUBROUTINE splint

!***********************************************************************
!
!.... Subroutine to calculate R, which is the distance between C and the
!.... center of H2 in CH5+
!
!.....History:
!.... Date     Modified by  Comment
!.... ----     -----------  -------
!.... 02/20/05   Zhong      created
!.... 03/01/05   Zhong      refine it

!.... dist ....  R
!.... h2dist ... bond length of H2
!.... theta ...  the angle between H2 and axis of C-center of H2
! 
************************************************************************
      SUBROUTINE calcr(x,xch3,xh2a,dist,h2dist,theta)

      implicit none
      integer :: i, j, npp1, npp2, jr
      integer, parameter :: natoms = 6

      double precision, dimension(3,natoms) :: x
      double precision, dimension(natoms-1) :: rch
      double precision, dimension(3) :: xh2m, xh1, xh2
      double precision, dimension(2) :: rchmax
      double precision :: dist, h2dist, h2dist5, theta
      double precision, dimension(3,4) :: xch3t
      double precision, dimension(3,2) :: xh2at
      double precision, dimension(4,3) :: xch3
      double precision, dimension(2,3) :: xh2a
       
      do j = 1, natoms - 1
         rch(j) = sqrt((x(1,j+1)-x(1,1))**2+(x(2,j+1)-x(2,1))**2
     $            +(x(3,j+1)-x(3,1))**2)
      end do
      rchmax(1) = 0.d0
      rchmax(2) = 0.d0
      npp1 = 0
      npp2 = 0

      do j = 1, natoms - 1
         if (rchmax(1)<rch(j)) then
            rchmax(1) = rch(j)
            npp1 = j
         end if
      end do

      do j = 1, natoms - 1
         if (npp1==j) then
         else
            if (rchmax(2)<rch(j)) then
               rchmax(2) = rch(j)
               npp2 = j
            end if
         end if
      end do

!.... Find the Cartisian coordinate of the center of H2
      do i = 1, 3
         xh2m(i) = 0.5d0*(x(i,npp1+1)+x(i,npp2+1))
      end do

!.... Obtain the R -- dist

      call calcdist(x,xh2m,dist)

!.... The Cartisian coordinate of H2
      do i = 1, 3
         xh1(i) = x(i,npp1+1)
         xh2(i) = x(i,npp2+1)
      end do

      do i = 1, 3
         xh2at(i,1) = x(i,npp1+1)
         xh2at(i,2) = x(i,npp2+1)
      end do

      jr = 0
      do j = 1, natoms
         if (j==npp1+1.or.j==npp2+1) then
         else
            jr = jr + 1
            do i = 1, 3
               xch3t(i,jr) = x(i,j)
            end do
         end if
      end do

         
      do i = 1, 3
         do j = 1, 4
            xch3(j,i) = xch3t(i,j)
         end do
      end do

      do i = 1, 3
         do j = 1, 2
            xh2a(j,i) = xh2at(i,j)
         end do
      end do

!.... Calculate the bond distance of H2
      call calcdist(xh1,xh2,h2dist)

      h2dist5 = h2dist/2.d0
      call calcangle(rch(npp2),dist,h2dist5,theta)
      END SUBROUTINE calcr
      
      SUBROUTINE calcdist(x1,x2,dist)
!.... a subroutine to calculate the bond length
!.... x1, x2 are Cartisian coordinate
!.... 
      implicit none

      double precision, dimension(3) :: x1, x2
      double precision :: dist

      dist = dsqrt((x1(1) - x2(1))**2 + (x1(2) - x2(2))**2
     $     + (x1(3) - x2(3))**2)

      END SUBROUTINE calcdist
      
      SUBROUTINE calcangle(a,b,c,theta)
!.... program to calculate theta
!.... a, b, c are the length of sides of triangle

      implicit none

      double precision :: a, b, c, theta, pi, d

      pi = 4.d0*datan(1.d0)

      if (abs(b-c-a).le.1e-5) then
         theta = 0.d0
      else
         d = b**2 + c**2 - a**2
         theta = 180.d0*dacos(d/(2.d0*b*c))/pi
      end if

      END SUBROUTINE calcangle
