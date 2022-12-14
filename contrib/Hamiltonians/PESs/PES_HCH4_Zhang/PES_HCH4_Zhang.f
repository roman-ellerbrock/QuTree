! -------------------------------------------------------------------
! updated 2013/5/17
!
! overall : 47783 points, rmse= 7.1689 meV, rmse(<1.36)= 3.2870 meV
!
! ==> n15.txt <==
! Sep_04_f12a_15_add1, nn15_b
! 16835 points for H+CH4 asy, rmse= 1.9646 meV, rmse(<1.36)= 1.2314 meV
! 
! ==> n24.txt <==
! Aug_30_f12a_24, nn24_k
! 10186 points for H2+CH3 asy, rmse= 3.1841 meV, rmse(<1.36)= 2.1204 meV
! 
! ==> n60.txt <==
! Sep_04_f12a_60_add1, nn60_o
! 20762 points for CH5 INT, rmse= 10.4964 meV, rmse(<1.36)= 4.6293 meV
!
! -------------------------------------------------------------------
! input cartesian or distances in bohr, atom order:  c h h h h h
! output energy in hartree, h+ch4 asy = 0.d0
! -------------------------------------------------------------------
!
      function f(x_in)
      implicit none
      real*8 x_in(3,6),x(3,6),r(15),v,f
      real*8 a15,a24,a60
      x=x_in
      call cart2dist(x,r)
      call nninter1(r,v,a15,a24,a60)
      f=v+40.95573198d0 ! RHF-UCCSD(T)-F12A/AVTZ
      return
      end function f
!
! -------------------------------------------------------------------
!
      subroutine surf(v,x_in)
      implicit none
      real*8 x_in(3,6),x(3,6),r(15),v
      real*8 a15,a24,a60
      x=x_in
      call cart2dist(x,r)
      call nninter1(r,v,a15,a24,a60)
      v=v+40.95573198d0 ! RHF-UCCSD(T)-F12A/AVTZ
      return
      end subroutine surf
!
! -------------------------------------------------------------------
      subroutine potdriva(r_in,v)
      implicit none
      real*8 r_in(15),r(15),v
      real*8 a15,a24,a60
      r=r_in
      call nninter1(r,v,a15,a24,a60)
      v=v+40.95573198d0 ! RHF-UCCSD(T)-F12A/AVTZ
      return
      end subroutine potdriva
!
! -------------------------------------------------------------------
      module nnfunctions
      implicit none

      interface nsim
              module procedure nsim_1 ! one hidden layer_
              module procedure nsim_2 ! two hidden layers
      end interface

      contains

      subroutine nsim_1(rin,vx,ndim,neu,nxin,nxw1,nxb1,nxw2,nxb2,nxout)
      implicit none
      integer :: ndim,neu,i,j,dx
      double precision :: rin(ndim),vrange,rrange(ndim),dv(ndim)
      double precision :: r(ndim),vx,nxin(2,ndim),nxout(2),nxb2,rtmp
      double precision :: nxb1(neu),nxw1(ndim,neu),nxw2(neu),ax(neu)

      dx=0
      vx=0.0
      dv=0.0
      r=rin
      vrange=nxout(2)-nxout(1)
      ! mapminmax [-1,1]
      do i=1,ndim
        rrange(i)=nxin(2,i)-nxin(1,i)
        r(i)=2.d0*(r(i)-nxin(1,i))/rrange(i)-1.d0
c        write (6, *) i-1, nxin(2,i), nxin(1,i)
c        write (6, *) i-1, rrange(i), r(i)
      end do
      ! 1st layer_

      do i=1,neu
        rtmp=nxb1(i)
        do j=1,ndim
          rtmp=rtmp+nxw1(j,i)*r(j)
        end do
        ax(i)=tanh(rtmp)
c        write (6, *) i-1, ax(i)
      end do

      ! output layer_
      vx=nxb2
      do i=1,neu
        vx=vx+nxw2(i)*ax(i) 
      end do
c      write (6, *) vx, vrange, nxout(1);
      !reverse map
      vx=vrange*(vx+1.d0)/2+nxout(1)

      return
      if (dx.eq.0) return

      ! calculate first derivatives, dv/dr(i), dv(i)
      do i=1,ndim
        do j=1,neu
          dv(i)=dv(i)+nxw2(j)*nxw1(i,j)*(1-ax(j)**2)
        enddo
        dv(i)=dv(i)*vrange/rrange(i)
      enddo

      return
      end subroutine nsim_1

      subroutine nsim_2(rin,vx,ndim,neu1,neu2,trfcn,nxin,nxw1,nxb1,
     &nxw2,nxb2,nxw3,nxb3,nxout)
      implicit none
      integer :: ndim,neu1,neu2,i,j,k,dx
      real*8 :: r(ndim),rin(ndim),vx,nxin(2,ndim),nxout(2)
      real*8 :: nxw1(ndim,neu1),nxb1(neu1),nxw2(neu1,neu2),nxb2(neu2)
      real*8 :: nxw3(neu2),nxb3,ax(neu1),bx(neu2),rtmp
      real*8 :: dv(ndim),dvtm,vrange,rrange(ndim)
      character*8 :: trfcn(2)

      dx=0
      vx=0.0
      dv=0.0
      r=rin
      vrange=nxout(2)-nxout(1)
      ! mapminmax [-1,1]
      do i=1,ndim
        rrange(i)=nxin(2,i)-nxin(1,i)
        r(i)=2.0*(r(i)-nxin(1,i))/rrange(i)-1.0
      end do
      ! 1st layer_
      if (trfcn(1).eq.'tansig') then
        do j=1,neu1
          rtmp=nxb1(j)
          do i=1,ndim
            rtmp=rtmp+nxw1(i,j)*r(i)
          end do
          ax(j)=tanh(rtmp)
        end do
      else if (trfcn(1).eq.'purelin') then
        do j=1,neu1
          ax(j)=nxb1(j)
          do i=1,ndim
            ax(j)=ax(j)+nxw1(i,j)*r(i)
          end do
        end do
      else
        write(*,"('transferfunction',a8,'not defined')") trfcn(1)
        stop
      endif

      ! 2nd layer_
      if (trfcn(1).eq.'tansig') then
        do k=1,neu2
          rtmp=nxb2(k)
          do j=1,neu1
            rtmp=rtmp+nxw2(j,k)*ax(j)
          end do
          bx(k)=tanh(rtmp)
        end do
      else
        write(*,"('transferfunction',a8,'not defined')") trfcn(1)
        stop
      endif

      ! output layer_
      vx=nxb3
      do k=1,neu2
        vx=vx+nxw3(k)*bx(k)
      end do

      !reverse map
      vx=vrange*(vx+1.d0)/2+nxout(1)

      return
      if (dx.eq.0) return

      ! calculate first derivatives, dv/dr(i), dv(i)
      do i=1,ndim
        do k=1,neu2
          dvtm=0.d0
          do j=1,neu1
            dvtm=dvtm+nxw2(j,k)*nxw1(i,j)*(1-ax(j)**2)
          enddo
          dv(i)=dv(i)+nxw3(k)*dvtm*(1-bx(k)**2)
        enddo
        dv(i)=dv(i)*vrange/rrange(i)
      enddo

      return
      end subroutine nsim_2

      end module nnfunctions

! --------------------------------------------------------------------
      module nnparam1
      implicit none
      save

      integer,parameter :: ndim=15

      character*20 :: n01f='../data/HCH4/n60.txt'
      integer,parameter :: n01s1=50,n01s2=60
      character*20 :: n13f='../data/HCH4/n15.txt'
      integer,parameter :: n13s=80
      character*20 :: n22f='../data/HCH4/n24.txt'
      integer,parameter :: n22s=80

c double hidden layers
      real*8 :: n01w1(ndim,n01s1),n01b1(n01s1)
      real*8 :: n01w2(n01s1,n01s2),n01b2(n01s2)
      real*8 :: n01w3(n01s2),n01b3
      real*8 :: n01in(2,ndim),n01out(2)
      character*8 :: n01fcn(2)
c single hidden layer_
      real*8 :: n13w1(ndim,n13s),n13b1(n13s),n13w2(n13s),n13b2
      real*8 :: n13in(2,ndim),n13out(2)
c single hidden layer_
      real*8 :: n22w1(ndim,n22s),n22b1(n22s),n22w2(n22s),n22b2
      real*8 :: n22in(2,ndim), n22out(2)
 
      contains

      subroutine nninit
!init nn01
      open(100,file=trim(n01f),status='old', err=200)
      read(100,*)n01w1,n01b1,n01w2,n01b2,n01w3,n01b3,n01in,n01out,n01fcn
      close(100)
!init nn13
      open(100,file=trim(n13f),status='old', err=200)
      read(100,*) n13w1,n13b1,n13w2,n13b2,n13in,n13out
      close(100)
!init nn22
      open(100,file=trim(n22f),status='old', err=200)
      read(100,*) n22w1,n22b1,n22w2,n22b2,n22in,n22out
      close(100)

      return

200   print*, 'Cannot open file.'
      stop
      end subroutine nninit

      end module nnparam1
c---------------------------------------------------------------------
      subroutine nninter1(r0,vx,a13,a22,a40)
      use nnparam1
      use nnfunctions
      integer,save :: init_nn=0
      real*8 :: r0(ndim),vx,r(ndim)
      real*8 :: v13,v22,v40,a13,a22,a40,s13,s22
      real*8 :: logsig,x
      logsig(x)=1.d0/(1.d0+exp(-x))

      r=r0
      call sortdist_ch5(r)
      a13=0.0;a22=0.0;a40=0.0
      v13=0.0;v22=0.0;v40=0.0
      s13=20.0;s22=20.0

c     a13=logsig(s13*(r(5)-r(4)-3.0))
c     a22=(1-a13)*logsig(s22*(r(4)-4.0))
      a22=logsig(s22*(r(4)-4.0))
      a13=(1-a22)*logsig(s13*(r(5)-r(4)-3.0))
      a40=1-a13-a22

! h+ch4
      if( a13.gt.1e-7 ) then
c       if(r(4).gt.4.d0) then
c       v13=-40.7937d0
c       else
        call nsim(r,v13,ndim,n13s,n13in,n13w1,n13b1,n13w2,n13b2,n13out)
c       endif
      endif

! h2+ch3
      if( a22.gt.1e-7 ) then
        if (r(15).gt.4.d0) then
        v22=-40.7937d0
        else
c       if (r(15).gt.4.d0) r(15)=atan(r(15)-4.d0)/3.1415926d0+4.d0
c       if (r(15).gt.4.d0) r(15)=4.d0
        call nsim(r,v22,ndim,n22s,n22in,n22w1,n22b1,n22w2,n22b2,n22out)
        endif
      endif

! ch5
      if( a40.gt.1e-7 ) then
        call nsim(r,v40,ndim,n01s1,n01s2,n01fcn,n01in,
     &  n01w1,n01b1,n01w2,n01b2,n01w3,n01b3,n01out)
      endif

c      write (6, *) a13, v13, r(15)
      vx=a13*v13+a22*v22+a40*v40
      return

      end subroutine nninter1
! --------------------------------------------------------------------
c$$$      subroutine echo_time
c$$$      implicit none
c$$$      integer m1,m2,m3,m4,m5,m6,m7
c$$$      call gettim(m1,m2,m3,m4)
c$$$      call getdat(m5,m6,m7)
c$$$      write(*,777)m5,m6,m7,m1,m2,m3,m4
c$$$777   format(1x,i4.4,'/',i2.2,'/',i2.2,2x,
c$$$     &       i2.2,':',i2.2,':',i2.2,'.',i2.2)
c$$$      return
c$$$      end
! --------------------------------------------------------------------
! for CH5 or SiH5
      subroutine sortdist_ch5(r)
      implicit none
      real*8 r(15),tmp
      integer id(6,6)
      integer i,j,k
      parameter( id=reshape([ -1,  1,  2,  3,  4,  5, 
     &                         1, -1,  6,  7,  8,  9, 
     &                         2,  6, -1, 10, 11, 12, 
     &                         3,  7, 10, -1, 13, 14, 
     &                         4,  8, 11, 13, -1, 15, 
     &                         5,  9, 12, 14, 15, -1] , (/ 6,6 /)))

      ! sort to r1<=r2<=r3<=r4<=r5
      do i=2,5
      do j=i+1,6
      if(r(i-1).gt.r(j-1)) then
        ! exchange atom i with j
        do k=1,6
          if(id(i,k).eq.-1 .or. id(j,k).eq.-1) cycle
          tmp=r(id(i,k))
          r(id(i,k))=r(id(j,k))
          r(id(j,k))=tmp
        enddo
      endif
      enddo
      enddo

      return
      end
! --------------------------------------------------------------------
! cartesian to distances, 6 atom
      subroutine cart2dist(x,r)
      implicit none
      real*8 :: x(3,6),r(15),d(3)
      integer :: i,j,k,m

      r=0.d0
      k=0
      do i=1,5
      do j=i+1,6
        k=k+1
        d(1:3)=x(1:3,i)-x(1:3,j)
        r(k)=dsqrt(dot_product(d,d))
      end do
      end do

      return
      end

! -------------------------------------------------------------------
      subroutine potinihch4zhang
	use iso_c_binding
      use nnparam1
      use nnfunctions
      implicit none
      integer,save :: init_nn=0

      if (init_nn.eq.0) then
        init_nn=1
        call nninit
      endif
!      call koordreadf

      end
! -------------------------------------------------------------------

      subroutine potentialhch4zhang(v, x, npart)
      use iso_c_binding
      implicit none

      integer    natom,ndimx,ndimq,i, npart
      parameter (natom=6)
      parameter (ndimx=3*natom,ndimq=natom*3-6)
      
      real*8     v,vmax
      real*8     x(ndimx)
      parameter (vmax=5.0d0/27.2114d0)

c      call IntToCartf(q,x)
c      print*, 'q='
c      print*, q(:)
c      print*, 'x='
c      print*, x(:)
      call surf(v,x)
c      v=v-0.04344984d0+0.115d0
c      write(6,*) v*27.2116d0*8065

      if (v.gt.vmax) v=vmax

      return

      end 

