! Version 15.07.2016
! input:
! 3 NO distances in bohr
! 2 ONO angles in degrees
! 1 umbrella in degrees
!
! tmc (or morse) for 3 distances
! scalling angles and umbrella
!
!  On output, S is corresponding coordinate vector in symmetry
!  coordinates in tmc Morse coordinates.
!  T currently contains only the hyperradius as first entry.
!
!cccccccccccccccccccccccccccccccc

      subroutine ctrans(q,s,t,skip,p)
      implicit none
      include 'morse.incl'
      include 'states.f'
      integer i,symm
      double precision q(*), s(*), t(*)                      !symmetry coordinates
      double precision x(6)
      double precision no1,no2,no3  !, oo1,oo2,oo3              !primitive coordinates
      double precision a1,a2,a3                              !O-N-O angles (3-1-4, 2-1-4, 2-1-3)
      double precision theta
      double precision  reno, reoo, rad, refangle    ! reference values
      parameter(reno=2.344419d0,reoo=4.0606528d0,rad=0.017453293d0  &
     &     ,refangle=120.d0)
      double precision sq2, sq3, sq6
      double precision f , p(20)                        ! exponent prefactor for tunable Morse coord.
      double precision  pi
      logical skip, angle

      pi=dacos(-1.d0)

      symm=2

      sq2=1.d0/dsqrt(2.d0)
      sq3=1.d0/dsqrt(3.d0)
      sq6=1.d0/dsqrt(6.d0)


!..   do no transformation if symm=0:
      if (symm.eq.0) then
      do i=1,qn
        s(i) = q(i)
      enddo
      return
      endif

!..   get primitive no displacement coordinates:
        no1 = q(1)   !- reno
        no2 = q(2)   !- reno
        no3 = q(3)   !- reno

      skip=.false.
      if (symm.eq.1) then
         write(6,*) 'symm .eq. 1 not coded'
         stop
      elseif (symm.eq.2) then
!..     transform distances to Morse or TMC  coordinates:
        if (f(no1,p).lt.-3.d0) skip=.true.
        if (f(no2,p).lt.-3.d0) skip=.true.
        if (f(no3,p).lt.-3.d0) skip=.true.
        x(1) = 1.d0 - dexp(-f(no1,p)*no1)
        x(2) = 1.d0 - dexp(-f(no2,p)*no2)
        x(3) = 1.d0 - dexp(-f(no3,p)*no3)
      endif
      do i=1,3
        if (x(i).lt.-8.d0) skip=.true.
      enddo
      if (skip) then
        write(6,*) 'tmc transformation unreasonable'
        return
      endif

!      write(6,100)(q(i), i=1,6),(x(i), i=1,6)
! 100  format('Q',6f10.3,' X ',6f10.3)
! 101  format('Q',6f10.3,' Fv ',6f10.3)

!--------------------------------------------------------------------------
!...  trisector angle:
      theta=q(7)*rad     !convert to radians

! projected angles directly from input
! correcting to get displacement
      a1=q(4)
      a2=q(5)
      a3=q(6)
      a1=a1*rad
      a2=a2*rad
      a3=a3*rad

!...  transform angles so they properly dissociate:
      no1=q(1)+reno
      no2=q(2)+reno
      no3=q(3)+reno
      a1=(a1)/(no2*no3)
      a2=(a2)/(no3*no1)
      a3=(a3)/(no1*no2)

!..   scale umbrella angle:
      theta=theta/(no1*no2*no3)      ! scale umbrella angle for proper dissociation
      theta=theta*reno**3            ! factor to get scaled angle to resonable value around reference


!..   project on symmetry basis to get new coordinate sets:
      do i=1,6
         s(i)=0.d0
      enddo
      s(1)=sq3*(x(1)+x(2)+x(3))
      s(2)=theta                     ! scaled trisector angle
      s(3)=sq6*(2.d0*x(1)-x(2)-x(3))
      s(4)=sq2*(         +x(2)-x(3))
      s(5)=sq6*(2.d0*a1-a2-a3)
      s(6)=sq2*(       +a2-a3)

!      write(6,123) (q(i), i=1,3), (s(i), i=1,3)
!      write(6,124) (q(i), i=1,6), (s(i), i=1,6)
123   format('old:  ',3f10.3,10x,'new:  ',3f10.5)
124   format('old:  ',6f10.3,2x,'new:  ',6f10.5)

!     create t set of additional coordinates:
      do i=1,qn
         t(i)=0.d0
      enddo
      t(1)=no1 +no2 +no3      !symmetrized NO distances for "hyper radius"
!
!      do i=1,6
!       Write(6,*)'s ctran',i,s(i)
!      enddo

      Return
      end



!------------------------------------------------------------------------------
! returns exponent of tunable Morse coordinate
! exponent is polynomial * gaussian (skewed)
      function f(x,p)

      implicit none
      integer ii,i, pn, nmax
      double precision x, r, f, gaus, poly, skew
      integer npoly(2)
      double precision a(2,20), p(20)

      common /tmc/a,npoly

! a(1):        position_ of equilibrium
! a(2):        constant of exponent
! a(3):        constant for skewing the gaussian
! a(4):        tuning for skewing the gaussian
! a(5):        Gaussian exponent
! a(6):        Shift of Gaussian maximum
! a(7)...:     polynomial coefficients
! a(8+n)...:   coefficients of Morse Power series

! Tunable Morse function
! Power series in Tunable Morse coordinates of order m
! exponent is polynomial of order npoly * gaussian + switching function

!.....set r  r-r_e
      r=x            !/0.529177d0
      r=r-p(1)

      p(2)=abs(p(2))
      p(5)=abs(p(5))

!..   set up_ skewing function:
      skew=0.5d0 * p(3) *(tanh(abs(p(4))*(r-p(6))) + 1.d0)

!..   set up_ gaussian function:
      gaus=exp(-p(5)*(r)**2)

!..   set up_ power series:
      poly=0.d0
      do i=0,npoly(1)-1
        poly=poly+p(7+i)*r**i
      enddo

!..   set up_ full exponent function:
      f=p(2)  + skew + gaus*poly
!      write(6,*) 'p(2):', p(2)
!      write(6,*) 'skew:', skew
!      write(6,*) 'gaus:', gaus
!      write(6,*) 'ploy:', poly
!      write(6,*)

      end

!--------------------------------------------------------------------------
!routine to transform valence coordiantes q to cartesian coordinates qcart
      subroutine cartesian(q,qcart)
      implicit none
      include 'states.incl'
      real*8 r1,r2,r3,a1,a2,th
      real*8 q(qn), qcart(12)
      real*8 reno, reoo, rad
      parameter(reno=2.344419d0,reoo=4.0606528d0,rad=0.017453293d0)

      r1=q(1)
      r2=q(2)
      r3=q(3)
      a1=q(4)*rad
      a2=q(5)*rad
      th=q(6)*rad

! N position_ vector, xyz
      qcart(1)=0.d0
      qcart(2)=0.d0
      qcart(3)=0.d0

      qcart(4)=r1*sin(th)
      qcart(5)=0.d0
      qcart(6)=r1*cos(th)

      qcart(7)=r2*sin(th)
      qcart(8)=r2*sin(a1)*cos(th)
      qcart(9)=r2*cos(a1)*cos(th)

      qcart(10)=r3*sin(th)
      qcart(11)=-r3*sin(a2)*cos(th)
      qcart(12)=r3*cos(a2)*cos(th)

!      geometry={
!      a1,n,xx1,yy1,zz1       qcart(1),qcart(5)...
!      a2,o,xx2,yy2,zz2       qcart(2)
!      a3,o,xx3,yy3,zz3       qcart(3)
!      a4,o,xx4,yy4,zz4}      qcart(4)


      end

!------------------------------------------------------------------
! routine to get N-O and O-O distances:
      subroutine internal(xx,qint)
      implicit none
      integer i, j, nat
      real*8 qint(7), xx(3,4), dum(3), norm, sc, scalar
      real*8 no(3), oo(3), r1, r2, r3
      real*8 v1(3), v2(3), v3(3)
      real*8 n1(3), n2(3), n3(3), tr(3)
      real*8 reno, reoo, rad
      parameter(reno=2.344419d0,reoo=4.0606528d0)

      rad=acos(-1.d0)/180.d0

777   format(3f12.6)
!      do i=1,4
!        write(6,777) (xx(j,i), j=1,3)
!      enddo

!..   get N-O vectors and distances:
      do i=1,3
        v1(i)=xx(i,2)-xx(i,1)
	v2(i)=xx(i,3)-xx(i,1)
	v3(i)=xx(i,4)-xx(i,1)
      enddo
      no(1)=norm(v1,3)
      no(2)=norm(v2,3)
      no(3)=norm(v3,3)

!..   write results onto internal coordiante array:
      qint(1)=no(1)-reno
      qint(2)=no(2)-reno
      qint(3)=no(3)-reno

!..   determine trisector:
      do i=1,3
        v1(i)=v1(i)/no(1)
        v2(i)=v2(i)/no(2)
        v3(i)=v3(i)/no(3)
      enddo

!..   compute three normal vectors for the ONO planes:
      call xprod(n1,v1,v2)
      call xprod(n2,v2,v3)
      call xprod(n3,v3,v1)
!      write(6,100) n1
!      write(6,100) n2
!      write(6,100) n3
!      write(6,*)

!..   compute trisector:
      do i=1,3
        tr(i)=(n1(i)+n2(i)+n3(i))/3.d0
      enddo
      r1=norm(tr,3)
      do i=1,3
        tr(i)=tr(i)/r1
      enddo
!      write(6,100) tr
!      write(6,*)
!      write(6,100) scal(v1,tr,3), scal(v2,tr,3), scal(v3,tr,3)
!      write(6,*)

!..   determine trisector angle:
      sc=scalar(v1,tr,3)
      qint(7)=90.d0 - acos(sc)/rad
      qint(7)=sign(qint(7),xx(1,2))

!..   molecule lies in yz plane, compute projected ONO angles:
      v1(1)=0.d0
      v2(1)=0.d0
      v3(1)=0.d0
      r1=norm(v1,3)
      r2=norm(v2,3)
      r3=norm(v3,3)
      do i=2,3
        v1(i)=v1(i)/r1
        v2(i)=v2(i)/r2
        v3(i)=v3(i)/r3
      enddo
      qint(4)=acos(scalar(v2,v3,3))/rad-120.d0     !projected ONO angles in radians
      qint(5)=acos(scalar(v3,v1,3))/rad-120.d0
      qint(6)=acos(scalar(v1,v2,3))/rad-120.d0


100   format(3f14.9)

      end

!-------------------------------------------------------------------
! compute vector product n1 of vectors v1 x v2
      subroutine xprod(n1,v1,v2)
      implicit none

      real*8 n1(3), v1(3), v2(3)

      n1(1) = v1(2)*v2(3) - v1(3)*v2(2)
      n1(2) = v1(3)*v2(1) - v1(1)*v2(3)
      n1(3) = v1(1)*v2(2) - v1(2)*v2(1)

      end


!-------------------------------------------------------------------
! compute scalar product of vectors v1 and v2:
      real*8 function scalar(v1,v2,n)
      implicit none
      integer i, n
      real*8 v1(*), v2(*)

      scalar=0.d0
      do i=1,n
        scalar=scalar+v1(i)*v2(i)
      enddo

      end

!-------------------------------------------------------------------
! compute norm of vector:
      real*8 function norm(x,n)
      implicit none
      integer i, n
      real*8 x(*)

      norm=0.d0
      do i=1,n
        norm=norm+x(i)**2
      enddo
      norm=sqrt(norm)

      end

