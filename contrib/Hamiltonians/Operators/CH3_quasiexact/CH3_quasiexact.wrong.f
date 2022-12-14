c######################################################################
c MCTDH-Modul Operatoren
c
c Enthaelt die Unterprogramme "h init" und "h" sowie dazu lokale
c Unterprogramme. "h" berechnet die Wirkung einer Komponente des
c Hamiltonoperators dargestellt in der Primaerbasis auf eine Funktion
c der Sekundaerbasis. Der Uebergabeparameter "mode" spezifiziert den
c Freiheitsgrad waehrend "teil" die Nummer der Kompentente des
c Hamilton-Operators angibt. Die Matrizen "trafo" und "matrix" werden
c gegebenenfalls benoetigt. "h init" initialiert das MCTDH-Programm
c fuer die Benutzung des Unterprogramms "h" sowie die Felder
c "koeff" und "diag".
c
c V1.0, 26.2.1991
c######################################################################

      subroutine h init chhh(koeff,diag,dim)

      implicit none

      complex*16  koeff(*)
      integer     dim, k zahl, i
      logical*4   diag(3)

      k zahl=58
      do i=1,k zahl
         koeff(i)=1
      enddo
      koeff(k zahl+1)=0

c Bestimmung der diagonalen Komponenten des Hamilton-Operators

      call diagonchhh(diag,dim,k zahl)
      return
      end

c--------------------------------------------------------------

      subroutine diagonchhh(diag,dim,k zahl)

      implicit none

      integer     dim, k zahl, i, k
      logical*4   diag(dim,k zahl)

      do k=1,k zahl
         do i=1,dim
            diag(i,k)=.false.
         enddo
      enddo

      return
      end

c##############################################################

      subroutine hchhh (mode,teil,h psi,psi,dim,matrix,trafo,ort)

      implicit none

      integer       mode, teil, dim, i, j
      complex*16    h psi(dim), psi(dim), work(1024)
      real*8        matrix(dim,dim,*), trafo(dim,dim), ort(dim)
      real*8        vol, G

      real*8        pi
      parameter     (pi=3.1415926535897931d0)

c structure
c row 1:  11, 22, 33
c row 2:  44[1, 2, 3, 4]
c row 3:  45[1, 2, 3, 4, 5]
c row 4:  54[1, 2, 3, 4, 5]
c row 5:  46[1, 2, 3, 4, 5]
c row 6:  64[1, 2, 3, 4, 5]
c row 7:  55[1, 2, 3, 4, 5, 6]
c row 8:  56[1, 2, 3, 4, 5, 6, 7]
c row 9:  65[1, 2, 3, 4, 5, 6, 7]
c row 10: 66[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
      goto (1000,2000,3000,4000,5000,6000) mode
 1000 goto (
     .      2,8,8,
     .      8,8,8,8,
c    .      8,
     .      8,8,8,8,8,
     .      8,8,8,8,8,
     .      8,8,8,8,8,
     .      8,8,8,8,8,
     .      8,8,8,8,8,8,
c    .      8,
     .      8,8,8,8,8,8,8,
     .      8,8,8,8,8,8,8,
c    .      8
     .      8,8,8,8,8,8,8,8,8,8,8
     .      ) teil
 2000 goto (
     .      1,2,8,
     .      8,8,8,8,
c    .      1,
     .      8,8,8,8,8,
     .      8,8,8,8,8,
     .      8,8,8,8,8,
     .      8,8,8,8,8,
     .      8,8,8,8,8,8,
c    .      1,
     .      8,8,8,8,8,8,8,
     .      8,8,8,8,8,8,8,
c    .      1
     .      8,8,8,8,8,8,8,8,8,8,8
     .      ) teil
 3000 goto (
     .      1,1,2,
     .      8,8,8,1,
c    .      1,
     .      8,8,8,8,1,
     .      8,8,8,8,1,
     .      8,8,8,8,1,
     .      8,8,8,8,1,
     .      8,8,8,8,8,1,
c    .      1,
     .      8,8,8,8,8,8,1,
     .      8,8,8,8,8,8,1,
c    .      1
     .      8,8,8,8,8,8,8,8,8,8,1
     .      ) teil
 4000 goto (
     .      1,1,1,
     .      2,2,2,2,
c    .      2,
     .      5,5,5,5,5,
     .      7,7,7,7,7,
     .      5,5,5,5,5,
     .      7,7,7,7,7,
     .      8,8,1,1,8,8,
c    .      8,
     .      8,8,1,1,1,8,8,
     .      8,8,1,1,1,8,8,
c    .      8
     .      8,8,8,8,1,1,8,8,8,8,1
     .      ) teil
 5000 goto (
     .      1,1,1,
     .      8,8,8,8,
c    .      1,
     .      7,7,7,7,7,
     .      5,5,5,5,5,
     .      8,8,8,8,8,
     .      8,8,8,8,8,
     .      3,3,3,3,3,3,
c    .      2,
     .      5,5,5,5,5,5,5,
     .      7,7,7,7,7,7,7,
c    .      1
     .      8,8,8,8,8,8,8,8,8,8,8
     .      ) teil
 6000 goto (
     .      1,1,1,
     .      8,8,8,8,
c    .      1,
     .      8,8,8,8,8,
     .      8,8,8,8,8,
     .      7,7,7,7,7,
     .      5,5,5,5,5,
     .      8,8,8,8,8,8,
c    .      1,
     .      7,7,7,7,7,7,7,
     .      5,5,5,5,5,5,5,
c    .      3
     .      3,3,3,3,3,3,3,3,3,3,3
     .      ) teil

c      H = 1

 1    do i=1,dim
         h psi(i)=psi(i)
      enddo
      return

c      H = -0.5 1/sqrt(v) d/dx v d/dx 1/sqrt(v)

 2    do i=1,dim
         h psi(i)=psi(i)/sqrt(vol(ort(i),mode))
      enddo
      call ddx (mode, work, h psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=work(i)*vol(ort(i),mode)
      enddo
      call ddx (mode, work, h psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-0.5d0*work(i)/sqrt(vol(ort(i),mode))
      enddo
      return

c      H = -0.5 1/sqrt(v) d/dx v G d/dx 1/sqrt(v)

 3    do i=1,dim
         h psi(i)=psi(i)/sqrt(vol(ort(i),mode))
      enddo
      call ddx (mode, work, h psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=work(i)*vol(ort(i),mode)*G(ort(i),mode,teil)
      enddo
      call ddx (mode, work, h psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-0.5d0*work(i)/sqrt(vol(ort(i),mode))
      enddo
      return

c      H =  sqrt(v) d/dx 1/sqrt(v)

 4    do i=1,dim
         h psi(i)=psi(i)/sqrt(vol(ort(i),mode))
      enddo
      call ddx (mode, work, h psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=work(i)*sqrt(vol(ort(i),mode))
      enddo
      return

c      H = G sqrt(v) d/dx 1/sqrt(v)

 5    do i=1,dim
         h psi(i)=psi(i)/sqrt(vol(ort(i),mode))
      enddo
      call ddx (mode, work, h psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=work(i)*sqrt(vol(ort(i),mode))*G(ort(i),mode,teil)
      enddo
      return

c      H = -0.5 /sqrt(v) d/dx sqrt(v)

 6    do i=1,dim
         h psi(i)=psi(i)*sqrt(vol(ort(i),mode))
      enddo
      call ddx (mode, work, h psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-0.5*work(i)/sqrt(vol(ort(i),mode))
      enddo
      return

c      H = -0.5 /sqrt(v) d/dx sqrt(v) G

 7    do i=1,dim
         h psi(i)=psi(i)*sqrt(vol(ort(i),mode))*G(ort(i),mode,teil)
      enddo
      call ddx (mode, work, h psi, dim, matrix, trafo)
      do i=1,dim
         h psi(i)=-0.5*work(i)/sqrt(vol(ort(i),mode))
      enddo
      return

c      H = G

 8    do i=1,dim
         h psi(i)=psi(i)*G(ort(i),mode,teil)
      enddo
      return

      end

c--------------------------------------------------------------

      real*8 function vol(x,mode)

      real*8 pi,pi3,x
      integer mode
      parameter     (pi=3.1415926535897931d0)
      pi3=pi/3.d0

      select case (mode)
         case(1)
            vol=x**2
         case(2)
            vol=dsin(x)
         case(3)
            vol=1.d0
         case(4)
            vol=dsin(x)**3
         case(5)
            vol=  1.d0/(1.d0
     .                +(x-pi3)**2
     .                -(x-pi3)**3/(3.d0*dsqrt(3.d0))
     .                +(x-pi3)**4*3.d0/4.d0)
c           vol=(dcos(x)+1.d0)/dsin(x)
         case(6)
            vol=  1.d0/(1.d0
     .                +(x-pi )**2/3.d0
     .                +(x-pi )**4/12.d0)
c           vol=(0.5d0-dcos(x))
      end select

      return

      end

c--------------------------------------------------------------

      real*8 function G(x,mode,teil)

      real*8 x,rhoeq,pi,pi3
      integer mode,teil
      parameter     (pi=3.1415926535897931d0)
      pi3=pi/3.d0


      rhoeq=dsqrt(3.d0*1836.109d0)*2.04035

      G=1.d0

      if(mode.eq.1.and.(teil.ge.4.and.teil.le.90))then
         G=16.d0/27.d0
      endif

c gamma phi
      if(mode.eq.5.and.(teil.ge.4.and.teil.le.90))then
         G=( 1.d0
     .      +(x-pi3)**2
     .      -(x-pi3)**3/(3.d0*dsqrt(3.d0))
     .      +(x-pi3)**4*3.d0/4.d0)**2
c        G=9.d0/4.d0/(dcos(x)+1.d0)**2/dsin(x)**2
      endif

c gamma chi
      if(mode.eq.6.and.(teil.ge.4.and.teil.le.90))then
         G=( 1.d0
     .      +(x-pi)**2/ 3.d0
     .      +(x-pi)**4/12.d0)**2
c        G=1.d0/(1.d0+dcos(x))**2
      endif


      goto (1000,2000,3000,4000,5000,6000) mode
 1000 goto (
     .         1, 10, 10,                         ! 11,22,33 terms
     .        14, 14, 17, 14,                     ! the 44 term
c    .        10,
     .        17, 16, 15, 16, 13,                 ! the 45 term
     .        17, 16, 15, 16, 13,                 ! the 54 term
     .        15, 14, 17, 16, 13,                 ! the 46 term
     .        15, 14, 17, 16, 13,                 ! the 64 term
     .        12, 15, 15, 15, 13, 10,             ! the 55 term
c    .        21,                                 ! c3v 55 term
     .        12, 17, 14, 15, 14, 13, 10,         ! the 56 term
     .        12, 17, 14, 15, 14, 13, 10,         ! the 65 term
c    .        10
     .10, 13, 12, 19, 14, 15, 14, 13, 11, 10, 11  ! the 66 term
c    .                                            10 ! c3v 66 term
     .      ) teil
 2000 goto (
     .         1,  1, 32,                         ! 11,22,33 terms
     .        32, 32, 32, 33,                     ! the 44 term
c    .         1,
     .        32, 32, 32, 32, 33,                 ! the 45 term
     .        32, 32, 32, 32, 33,                 ! the 54 term
     .        32, 32, 32, 32, 33,                 ! the 46 term
     .        32, 32, 32, 32, 33,                 ! the 64 term
     .        32, 32, 32, 32, 32, 33,             ! the 55 term
c    .         1,                                 ! c3v 55 term
     .        32, 32, 32, 32, 32, 32, 33,         ! the 56 term
     .        32, 32, 32, 32, 32, 32, 33,         ! the 65 term
c    .         1
     .32, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33  ! the 66 term
c    .                                            1 ! c3v 56 term
     .      ) teil
 3000 goto (
     .         1,  1,  1,                         ! 11,22,33 terms
     .        34, 34, 35,  1,                     ! the 44 term
c    .         1,
     .        34, 35, 34, 35,  1,                 ! the 45 term
     .        34, 35, 34, 35,  1,                 ! the 54 term
     .        34, 35, 34, 35,  1,                 ! the 46 term
     .        34, 35, 34, 35,  1,                 ! the 64 term
     .        34, 34, 34, 34, 35,  1,             ! the 55 term
c    .         1,                                 ! c3v 55 term
     .        34, 34, 35, 35, 35, 35,  1,         ! the 56 term
     .        34, 34, 35, 35, 35, 35,  1,         ! the 65 term
c    .         1
     .34, 34, 34, 34, 34, 34, 35, 35,  1,  1,  1  ! the 66 term
c    .                                            1  ! c3v 66 term
     .      ) teil
 4000 goto (
     .         1,  1,  1,                         ! 11,22,33 terms
     .         1,  1,  1,  1,                     ! the 44 term
c    .         1,
     .        40, 40, 40, 40, 40,
     .        40, 40, 40, 40, 40,                 ! the 45 term
     .        40, 40, 40, 40, 40,                 ! the 54 term
     .        40, 40, 40, 40, 40,                 ! the 46 term
     .        43, 43,  1,  1, 41, 41,             ! the 55 term
c    .        45,                                 ! c3v 55 term
     .        41, 41,  1,  1,  1, 43, 41,         ! the 56 term
     .        41, 41,  1,  1,  1, 43, 41,         ! the 65 term
c    .        46
     .41, 41, 41,  1,  1,  1, 41, 41, 43, 43,  1  ! the 66 term
c    .                                            46 ! c3v 66 term
     ,      ) teil
 5000 goto (
     .         1,  1,  1,                         ! 11,22,33 terms
     .        52, 53, 58, 58,                     ! the 44 term
c    .         1,
     .        56, 52, 54, 60, 61,                 ! the 45 term
     .        56, 52, 54, 60, 61,                 ! the 54 term
     .        66, 65, 71, 56, 60,                 ! the 46 term
     .        66, 65, 71, 56, 60,                 ! the 64 term
     .        52, 60, 52, 55, 54, 55,             ! the 55 term
c    .         1,                                 ! c3v 55 term
     .        50, 56, 67, 60, 52, 52, 54,         ! the 56 term
     .        50, 56, 67, 60, 52, 52, 54,         ! the 65 term
c    .         1
     .53, 76, 64, 58, 52, 60, 56, 65, 60, 67, 52  ! the 66 term
c    .                                            1 ! c3v 66 term
     .      ) teil
 6000 goto (
     .         1,  1,  1,                         ! 11,22,33 terms
     .        53, 52, 58,  1,                     ! the 44 term
c    .         1,
     .        52, 56, 51, 50,  1,                 ! the 45 term
     .        52, 56, 51, 50,  1,                 ! the 54 term
     .        50, 51, 56, 71, 50,                 ! the 46 term
     .        50, 51, 56, 71, 50,                 ! the 64 term
     .         1, 51, 52,  1, 50,  1,             ! the 55 term
c    .         1,                                 ! c3v 55 term
     .        50, 56,  1, 51, 52, 52, 50,         ! the 56 term
     .        50, 56,  1, 51, 52, 52, 50,         ! the 65 term
c    .         1
     .53, 51,  1,  1, 53, 51, 56, 50, 51,  1, 52  ! the 66 term
c    .                                            1 ! c3v 66 term
     .      ) teil

1     return

10    G=+G*1.d0/x**2        ! 1/x**2
      return
11    G=G*(-1.d0)/x**2
      return
12    G=+G/2.d0/x**2
      return
13    G=-G/2.d0/x**2
      return
14    G=+G/4.d0/x**2
      return
15    G=-G/4.d0/x**2
      return
16    G=+G/8.d0/x**2
      return
17    G=-G/8.d0/x**2
      return
18    G=+G*3.d0/16.d0/x**2
      return
19    G=+G     /16.d0/x**2
      return
20    G=+G*3.d0/2.d0/x**2
      return
21    G=+G*3.d0/4.d0/x**2
      return
22    G=+G*9.d0/4.d0/x**2
      return
23    G=+G     /9.d0/x**2
      return
24    G=+G*4.d0/9.d0/x**2
      return

32    G=G/dsin(x)**2                        ! 1/sin^2 x
      return
33    G=G/dcos(x)**2                        ! 1/cos^2 x
      return
34    G=G/dsin(x)**2/dcos(x)**2             ! 1/sin^2 x 1/cos^2 x
      return
35    G=G*(1.d0/dcos(x)**2-1.d0/dsin(x)**2) ! (1/cos^2 x - 1/sin^2 x)
      return
36    G=G/dsin(2.d0*x)**2                   ! ( 1/sin^2 2x)
      return

40    G=G/dtan(x)         ! cot x
      return
41    G=G/dtan(x)**2      ! cot^2 x
      return
42    G=G/dsin(x)         ! csc x
      return
43    G=G/dsin(x)**2      ! csc^2 x
      return
45    G=G*(3.d0+dcos(2.d0*x))/dsin(x)**2    ! (3 + cos 2x) csc^2 x
      return
46    G=G*9.d0/4.d0*(3.d0+dcos(2.d0*x))/dsin(x)**2
      return


50    G=G*dsin(x)         ! sin x
      return
51    G=G*dcos(x)         ! cos x
      return
52    G=G*dsin(x)**2      ! sin^2 x
      return
53    G=G*dcos(x)**2      ! cos^2 x
      return
54    G=G*dsin(x)**3      ! sin^3 x
      return
55    G=G*dsin(x)**4      ! sin^4 x
      return
56    G=G*dsin(2.d0*x)    ! sin 2x
      return
57    G=G*dcos(2.d0*x)    ! cos 2x
      return
58    G=G*dsin(2.d0*x)**2 ! sin^2 2x
      return
59    G=G*dcos(2.d0*x)**2 ! cos^2 2x
      return

60    G=G*dsin(x)   *dsin(2.d0*x)             ! sin x   sin 2x
      return
61    G=G*dsin(x)**2*dsin(2.d0*x)             ! sin^2 x sin 2x
      return
62    G=G           *(1.d0+dsin(x)**2)        !         (1 + sin^2 x)
      return
63    G=G*dsin(x)   *(1.d0+dsin(x)**2)        ! sin x   (1 + sin^2 x)
      return
64    G=G           *(1.d0+dcos(x)**2)        !         (1 + cos^2 x)
      return
65    G=G*dsin(x)   *(1.d0+dcos(x)**2)        ! sin x   (1 + cos^2 x)
      return
66    G=G*dcos(x)   *(1.d0+dcos(x)**2)        ! cos x   (1 + cos^2 x)
      return
67    G=G*dsin(x)**2*(1.d0+dcos(x)**2)        ! sin^2 x (1 + cos^2 x)
      return


71    G=G*        (1.d0-3.d0*dcos(x)**2)      !        1 - 3 cos^2 x
      return
72    G=G*dsin(x)*(2.d0+     dcos(x)**2)      ! sin x (2 -   cos^2 x)
      return
74    G=G*        (3.d0+dcos(x)**2)           !       (3 +   cos^2 x)
      return
75    G=G*dsin(x)*(3.d0+dcos(x)**2)           ! sin x (3 +   cos^2 x)
      return
76    G=G*dcos(x)*(3.d0+dcos(x)**2)           ! cos x (3 +   cos^2 x)
      return

      return

      end

	subroutine IntToCartch3qe(q, cart)
	use iso_c_binding
	implicit none
      real*8	q(*)
      integer	sys(4)
      real*8	para(4,6), mass(4),M
      real*8	rad(3,4),cart(3,4),alpha(4),C(4,4)
	real*8	hmass
	integer	i, j, n, l
      parameter (hmass=1836.1526838917d0)
      common /sys/ sys, para, mass

      mass(1)=12.d0
      mass(2)=1.d0
      mass(3)=1.d0
      mass(4)=1.d0

      M=mass(1)+mass(2)+mass(3)+mass(4)
      do i=1,4
         alpha(i)=dsqrt(mass(i)/M)
      enddo

c the mass weighted cartesian to mass weighted radua matrix
      do i=1,4
         C(1,i)=alpha(i)
         C(i,1)=alpha(i)
      enddo
      do i=2,4
      do j=2,4
         if(i.eq.j)then
            C(i,j)=alpha(i)*alpha(j)/(alpha(1)+1)-1
         else
            C(i,j)=alpha(i)*alpha(j)/(alpha(1)+1)
         endif
      enddo
      enddo

c Cartesian form for the radau vectors
      rad(1,1)=0.d0
      rad(2,1)=0.d0
      rad(3,1)=0.d0

      rad(1,2)=q(1)*dcos(q(2))           *dsin(q(4))
      rad(2,2)=0.d0
      rad(3,2)=q(1)*dcos(q(2))           *dcos(q(4))

      rad(1,3)=q(1)*dsin(q(2))*dcos(q(3))*dsin(q(4))*dcos(+q(5)+q(6))
      rad(2,3)=q(1)*dsin(q(2))*dcos(q(3))*dsin(q(4))*dsin(+q(5)+q(6))
      rad(3,3)=q(1)*dsin(q(2))*dcos(q(3))*dcos(q(4))

      rad(1,4)=q(1)*dsin(q(2))*dsin(q(3))*dsin(q(4))*dcos(-q(5)+q(6))
      rad(2,4)=q(1)*dsin(q(2))*dsin(q(3))*dsin(q(4))*dsin(-q(5)+q(6))
      rad(3,4)=q(1)*dsin(q(2))*dsin(q(3))*dcos(q(4))

c convert from mass weighted raduas to mass weighted cartestians
      do n=1,4
      do i=1,3
         cart(i,n)=0.d0
         do l=1,4
            cart(i,n)=cart(i,n)+C(l,n)*rad(i,l)
         enddo
      enddo
      enddo

c un-mass weight the coordinates
      do n=1,4
      do i=1,3
         cart(i,n)=cart(i,n)/dsqrt(mass(n)*hmass)
      enddo
      enddo

      end subroutine

