c     potential by schwenke for ch_4
c     modified by twestmann for OMP use
      subroutine potentialch4 (v,x)

      implicit none

      integer N
      parameter (N=5)
      real*8 vmax, q(9), qq(9), x(15)
      parameter (vmax=5.d0/27.212d0)
      real*8 v0
      parameter (v0=28.11193d0/219474.63d0)

      integer i,j
      real*8 v
      real*8 coor(3*N)
      real*8 coor2(0:N-1,3)
c      real*8 coor(3*N-6)

c      do i = 1, 9
c         qq(i) = q(i)
c      enddo

c      call IntToCartfch4(qq,x)

      do i=0,N-1
        do j=1,3
           coor2(i,j)=x(3*i+j)
        enddo
      enddo

      call potentialcart (v,coor2,N)
!      call potentialint (v,coor,N)
      v=v+v0

      if (v.gt.vmax) v=vmax

      return
      end

      subroutine potentialcart (v,x,N)

      implicit none

      integer N,i
      real*8 rad(0:N-1,3)
      real*8 rij(10)
      real*8 x(0:N-1,3),v, v_arr(10)

      call cartrad (x,rad,N)
      call radpotint (rad,rij,N)

c      do i=1,10
c         write (6,*) rij(i)
c      enddo
c      stop

      call vibpot (rij,v_arr,1)
	v = v_arr(1)

      return
      end

      subroutine potentialint (v,q,N)

      implicit none

      real*8    q(*), v, v_arr(10)
      real*8 en

      integer N
      real*8 cart(3*N)
      integer i, j
      real*8 mh,mc,ar
      real*8 r1, r2, theta
      real*8 x(0:N-1,3), rad(0:N-1,3)
      real*8 Cadj (0:N-1,0:N-1)
      real*8 m(0:N-1), alpha(0:N-1), mg
      real*8 dist(N-1,N-1)
      real*8 rij(10)

c     Masseninitialisierung

      call massen (mh,mc,ar)

      r1=q(1)
      r2=q(2)
      theta=q(3)

      m(0)=mc*ar
      mg=m(0)

      do i=1,N-1
         m(i)=mh*ar
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
         rad(i,1)=2*q(3*(i-2)+1)*q(3*(i-2)+2)
     .            /(1+q(3*(i-2)+2)**2+q(3*(i-2)+3)**2)

         rad(i,2)=2*q(3*(i-2)+1)*q(3*(i-2)+3)
     .            /(1+q(3*(i-2)+2)**2+q(3*(i-2)+3)**2)

c     negatives Vorzeichen z-Komponente bei Projektion vom Nordpol

       rad(i,3)=-1.0d0*q(3*(i-2)+1)*(1-q(3*(i-2)+2)**2-q(3*(i-2)+3)**2)
     .            /(1+q(3*(i-2)+2)**2+q(3*(i-2)+3)**2)
c H4 Suedpol !
      if (i.eq.(N-1)) rad(i,3)=-rad(i,3)
      enddo

c      write(*,*) 'Radau-Vektoren'
c      do i=0,N-1
c         write(*,20) rad(i,1),rad(i,2),rad(i,3)
c      enddo

c     Abstaende der H-Atome

c      write(*,*)'Abstaende der H-Atome, Radau-Vektoren'
c      do i=1,N-1
c         do j=i+1,N-1
c         dist(i,j)=(rad(i,1)-rad(j,1))**2+(rad(i,2)-rad(j,2))**2
c     .               +(rad(i,3)-rad(j,3))**2
c            dist(i,j)=sqrt(dist(i,j))/sqrt(m(i))
c            write(*,21) i,j,dist(i,j)
c      enddo
c      enddo

c     Berechnung der internen Koordinaten fuer Schwenke-pes aus Radaus

      call radpotint (rad,rij,N)

      call vibpot (rij,v_arr,1)
	v = v_arr(1)

c     Transformationsmatrix aufbauen
c
c     do i=0,N-1
c        Cadj(0,i)=alpha(i)
c        Cadj(i,0)=alpha(i)
c     enddo

c     do i=1,N-1
c        do j=1, N-1
c          Cadj(i,j)=alpha(i)*alpha(j)/(alpha(0)+1.0d0)
c        enddo
c     enddo

c     Diagonalelemente korrigieren

c     do i=1,N-1
c        Cadj(i,i)=Cadj(i,i)-1.0d0
c     enddo
c
c     Transformation der Radau-Vektoren

c     do i=0,N-1
c        x(i,1)=0.0d0
c        x(i,2)=0.0d0
c        x(i,3)=0.0d0
c     enddo

c     do i=0,N-1
c        do j=0,N-1
c           x(i,1)=x(i,1)+Cadj(i,j)*rad(j,1)
c           x(i,2)=x(i,2)+Cadj(i,j)*rad(j,2)
c           x(i,3)=x(i,3)+Cadj(i,j)*rad(j,3)
c        enddo
c     enddo

c     Massengewichtung rueckgaengig, Vorzeichenwechsel

c     do i=0,N-1
c        x(i,1)=-1.0d0*x(i,1)/sqrt(m(i))
c        x(i,2)=-1.0d0*x(i,2)/sqrt(m(i))
c        x(i,3)=-1.0d0*x(i,3)/sqrt(m(i))
c     enddo
c
c 20    format(3F11.6)
c 21    format(I3,I3,F11.6)

c     Koordinaten nach Umrechnung

c     x(0,1-3) xyz C
c     x(1,1-3) xyz H1
c     x(2,1-3) xyz H2
c     x(3,1-3) xyz H3
c     x(4,1-3) xyz H4

c     Koordinaten Potentialroutine
c
c x(1-3) xyz H1; x(7-9) xyz H3; x(13-15) xyz C
c x(4-6) xyz H2; x(10-12) xyz H4; x(16-18) xyz Hb

c     do i=1,N-1
c        do j=1,3
c           cart(3*(i-1)+j)=ang*x(i,j)
c        enddo
c     enddo
c
c     cart(13)=ang*x(0,1)
c     cart(14)=ang*x(0,2)
c     cart(15)=ang*x(0,3)
c
c     if (v.gt.vmax) v=vmax

      return
      end

      subroutine radpotint (rad,rij,N)

      implicit none

      integer N
      integer i,j
      real*8 mh,mc,ar
      real*8 rad(0:N-1,3)
      real*8 rij(10)

      call massen (mh,mc,ar)

c     Laengen der Radau-Vektoren

      do i=1,N-1
       rij(i)=0.d0
       do j=1,3
            rij(i)=rij(i)+rad(i,j)**2
         enddo
         rij(i)=dsqrt(rij(i))
      enddo


c     Winkel zwischen Radau-Vektoren

      rij( 5)=rad(1,1)*rad(2,1)+rad(1,2)*rad(2,2)+rad(1,3)*rad(2,3)
      rij( 6)=rad(1,1)*rad(3,1)+rad(1,2)*rad(3,2)+rad(1,3)*rad(3,3)
      rij( 7)=rad(1,1)*rad(4,1)+rad(1,2)*rad(4,2)+rad(1,3)*rad(4,3)
      rij( 8)=rad(2,1)*rad(3,1)+rad(2,2)*rad(3,2)+rad(2,3)*rad(3,3)
      rij( 9)=rad(2,1)*rad(4,1)+rad(2,2)*rad(4,2)+rad(2,3)*rad(4,3)
      rij(10)=rad(3,1)*rad(4,1)+rad(3,2)*rad(4,2)+rad(3,3)*rad(4,3)

      rij( 5)=rij( 5)/(rij(1)*rij(2))
      rij( 6)=rij( 6)/(rij(1)*rij(3))
      rij( 7)=rij( 7)/(rij(1)*rij(4))
      rij( 8)=rij( 8)/(rij(2)*rij(3))
      rij( 9)=rij( 9)/(rij(2)*rij(4))
      rij(10)=rij(10)/(rij(3)*rij(4))

      rij( 5)=dacos(rij( 5))
      rij( 6)=dacos(rij( 6))
      rij( 7)=dacos(rij( 7))
      rij( 8)=dacos(rij( 8))
      rij( 9)=dacos(rij( 9))
      rij(10)=dacos(rij(10))

c     ohne Massenwichtung

      do i=1,N-1
         rij(i)=rij(i)/dsqrt(mh*ar)
      enddo

      return
      end

      subroutine cartrad (x,rad,N)

      implicit none

      integer N
      integer i,j

      real*8 x(0:N-1,3),rad(0:N-1,3)
      real*8 mh,mc,ar
      real*8 r1, r2, theta
      real*8 C(0:N-1,0:N-1)
      real*8 m(0:N-1), alpha(0:N-1), mg

      call massen (mh,mc,ar)

      m(0)=mc*ar
      mg=m(0)

      do i=1,N-1
         m(i)=mh*ar
         mg=mg+m(i)
      enddo

      do i=0,N-1
         alpha(i)=sqrt(m(i)/mg)
      enddo

c     Transformationsmatrix aufbauen

      do i=0,N-1
         C(0,i)=alpha(i)
         C(i,0)=alpha(i)
      enddo

      do i=1,N-1
         do j=1, N-1
           C(i,j)=alpha(i)*alpha(j)/(alpha(0)+1.0d0)
         enddo
      enddo

c     Diagonalelemente korrigieren

      do i=1,N-1
         C(i,i)=C(i,i)-1.0d0
      enddo

c     Massenwichtung

      do i=0,N-1
         x(i,1)=x(i,1)*dsqrt(m(i))
         x(i,2)=x(i,2)*dsqrt(m(i))
         x(i,3)=x(i,3)*dsqrt(m(i))
      enddo

c     Transformation der kart. Vektoren

      do i=0,N-1
         rad(i,1)=0.0d0
       rad(i,2)=0.0d0
       rad(i,3)=0.0d0
      enddo

      do i=0,N-1
         do j=0,N-1
          rad(i,1)=rad(i,1)+C(i,j)*x(j,1)
          rad(i,2)=rad(i,2)+C(i,j)*x(j,2)
          rad(i,3)=rad(i,3)+C(i,j)*x(j,3)
       enddo
      enddo

c     Vorzeichenwechsel

      do i=0,N-1
         rad(i,1)=-1.0d0*rad(i,1)
         rad(i,2)=-1.0d0*rad(i,2)
       rad(i,3)=-1.0d0*rad(i,3)
      enddo

c      write(*,*) 'Radau-Vektoren'
c      do i=0,N-1
c         write(*,20) rad(i,1),rad(i,2),rad(i,3)
c      enddo
c      stop

20    format(3F11.6)

      return
      end

      subroutine massen (mh,mc,ar)

      implicit none

      real*8 mh,mc,ar

c     atomare Masseneinheit

      ar=1822.88848d0

c     isotopenreine Atommassen

c      mh=1.007825047d0
c      mc=12.d0

c     isotopenreine Kernmassen

      mh=1.00727646688d0
      mc=11.99341704d0

      return
      end

      subroutine vibpot(rij,v,n)
      implicit real*8 (a-h,o-z)
       character*80 title
       real*8 tocm,bohr,bohri
       common /pes1/  tocm,bohr,bohri
       real*8 re0,ex0,beta0,r01,r02,r03,r04,a02,a03,a04,a06,a08
       common /pes2/re0,ex0,beta0,r01,r02,r03,r04,a02,a03,a04,a06,a08
       real*8  ccc
       common /pes3/ ccc
       real*8 reau,ex,thrd,exa,srh,sr12,cs,cb,beta,const
       common /pes4/ reau,ex,thrd,exa,srh,sr12,cs,cb,beta,const
       real*8 cands
	 integer nexpr
       common /pes5/ cands,nexpr
       integer nexpa,ixa,ixr,nx,ixp,ixw,ixam
       common /pesint/ nexpa,ixa,ixr,nx,ixp,ixw
       integer id,n
       real*8 rij,v,fmat,rmat
      parameter (id=2000)
      dimension rij(n,10),v(n,1),cs(6),cb(6),                           7d10s01
     $          fmat(9,9),nx(8),ixp(id,8,8),ixw(id,8,8),ccc(id,8)       5d30s00
      dimension ixa(id,5),ixr(id,4),cands(id,id),rmat(id),              7d27s00
     $          ixam(id,5)                                              7d27s00

      integer intmul,maxbc
      real*8 tovr
      common/iszcm/intmul                                               12d8s88
      common/timrcm/tovr
      parameter (maxbc=1 000 000)
      integer ibcoff,ibcmax,ibca,mnbcv,ibc
      real*8 bc(maxbc)
      common ibcoff,ibcmax,ibca,mnbcv,bc
      dimension ibc(maxbc)
      equivalence (ibc,bc)
      integer i,j

!       if(ifirst.eq.0)then
!       tocm=1d0/4.556335d-6                                              5d23s91
!       bohr=0.529177249d0
!       bohri=1d0/bohr
! c       write(6,111)
!   111  format(/1x,'pes for ch4: fast version')
!        open(unit=1,file='pes',form='formatted',status='old')
!       read(23,1492)title                                                 8d1s00
!  1492 format(a80)                                                       8d1s00
! c      write(6,*)('title for PES:')                                      8d1s00
! c      write(6,1492)title                                                8d1s00
!       read(23,*)re0,ex0,beta0
!       read(23,*)r01,r02,r03,r04
!       read(23,*)a02,a03,a04,a06,a08                                      8d3s00
!       re0=re0*bohri                                                     8d1s00
!       r01=r01*6d0
!       r02=r02*6d0
!       r03=r03*6d0
!       r04=r04*6d0
!       a02=a02*4d0
!       a03=a03*4d0
!       a04=a04*4d0
! c      write(6,*)('zero order fit parameters: ')
! c      write(6,*)('re0,ex0,beta0 '),re0,ex0,beta0
! c      write(6,*)('r0i '),r01,r02,r03,r04
! c      write(6,*)('a0i '),a02,a03,a04,a06,a08                            8d3s00
! c      call flush(6)
!        call setupvp(nx,ixp,ixw,ccc,id,reau,ex,thrd,exa,srh,sr12,cs,cb,
!      $              beta,const,ixa,nexpa,ixr,cands,nexpr,r01,r02,       7d27s00
!      $              r03,r04,a02,a03,a04,a06,a08)                        10d6s00
! c       write(6,*)('in vibpot, ex= '),ex
!        do 510 j=1,5                                                     7d27s00
!         do 1676 i=1,nexpa                                               7d27s00
!          ixam(i,j)=ixa(i,j)-1
!  1676   continue                                                        7d27s00
!   510  continue                                                         7d27s00
!        ifirst=1
!       end if

      do 100 ig=1,n
       t1=cos(rij(ig,5))+thrd                                           7d24s00
       t2=cos(rij(ig,6))+thrd                                           7d24s00
       t3=cos(rij(ig,7))+thrd                                           7d24s00
       t4=cos(rij(ig,8))+thrd                                           7d24s00
       t5=cos(rij(ig,9))+thrd                                           7d24s00
       t6=cos(rij(ig,10))+thrd                                          7d24s00
       rr1=rij(ig,1)-re0                                                7d24s00
       rr2=rij(ig,2)-re0                                                7d24s00
       rr3=rij(ig,3)-re0                                                7d24s00
       rr4=rij(ig,4)-re0                                                7d24s00
       d1=exp(-beta0*(rr1**2+rr2**2))
       d2=exp(-beta0*(rr1**2+rr3**2))
       d3=exp(-beta0*(rr1**2+rr4**2))
       d4=exp(-beta0*(rr2**2+rr3**2))
       d5=exp(-beta0*(rr2**2+rr4**2))
       d6=exp(-beta0*(rr3**2+rr4**2))
       y1=1d0-exp(-ex0*rr1)
       y2=1d0-exp(-ex0*rr2)
       y3=1d0-exp(-ex0*rr3)
       y4=1d0-exp(-ex0*rr4)
       t12=t1*t1
       t22=t2*t2
       t32=t3*t3
       t42=t4*t4
       t52=t5*t5
       t62=t6*t6
       vs=y1*(r01+y1*(r02+y1*(r03+y1*r04)))
     $            +y2*(r01+y2*(r02+y2*(r03+y2*r04)))
     $            +y3*(r01+y3*(r02+y3*(r03+y3*r04)))
     $            +y4*(r01+y4*(r02+y4*(r03+y4*r04)))
     $            +d1*t12*(a02+t1*(a03+t1*(a04+t12*(a06+t12*a08))))     8d3s00
     $            +d2*t22*(a02+t2*(a03+t2*(a04+t22*(a06+t22*a08))))     8d3s00
     $            +d3*t32*(a02+t3*(a03+t3*(a04+t32*(a06+t32*a08))))     8d3s00
     $            +d4*t42*(a02+t4*(a03+t4*(a04+t42*(a06+t42*a08))))     8d3s00
     $            +d5*t52*(a02+t5*(a03+t5*(a04+t52*(a06+t52*a08))))     8d3s00
     $            +d6*t62*(a02+t6*(a03+t6*(a04+t62*(a06+t62*a08))))     8d3s00
              d1=rij(ig,1)-reau                                                   5d8s00
              d2=rij(ig,2)-reau                                                   5d8s00
              d3=rij(ig,3)-reau                                                   5d8s00
              d4=rij(ig,4)-reau                                                   5d8s00
       if(ex.gt.0d0)then                                                10d6s00
        d1=1d0-exp(-ex*d1)                                              10d11s00
        d2=1d0-exp(-ex*d2)                                              10d11s00
        d3=1d0-exp(-ex*d3)                                              10d11s00
        d4=1d0-exp(-ex*d4)                                              10d11s00
       end if                                                           10d6s00
              damp=exp(-beta*(d1*d1+d2*d2+d3*d3+d4*d4))                 7d24s00
              s1=0.5d0*(d1+d2+d3+d4)
              s2a=sr12*(2d0*t1-t2-t3-t4-t5+2d0*t6)
              s2b=0.5d0*(t2-t3-t4+t5)
              s3x=0.5d0*(d1-d2+d3-d4)
              s3y=0.5d0*(d1-d2-d3+d4)
              s3z=0.5d0*(d1+d2-d3-d4)
              s4x=srh*(t5-t2)
              s4y=srh*(t4-t3)
              s4z=srh*(t6-t1)
              fmat(1,1)=1d0
              fmat(1,2)=1d0
              fmat(1,3)=1d0
              fmat(1,4)=1d0
              fmat(1,5)=1d0
              fmat(1,6)=1d0
              fmat(1,7)=1d0
              fmat(1,8)=1d0
              fmat(1,9)=1d0
              fmat(2,1)=s1
              fmat(2,2)=s2a
              fmat(2,3)=s2b
              fmat(2,4)=s3x
              fmat(2,5)=s3y
              fmat(2,6)=s3z
              fmat(2,7)=s4x
              fmat(2,8)=s4y
              fmat(2,9)=s4z
              do 4353 i=3,9
               fmat(i,1)=fmat(i-1,1)*fmat(2,1)
               fmat(i,2)=fmat(i-1,2)*fmat(2,2)
               fmat(i,3)=fmat(i-1,3)*fmat(2,3)
               fmat(i,4)=fmat(i-1,4)*fmat(2,4)
               fmat(i,5)=fmat(i-1,5)*fmat(2,5)
               fmat(i,6)=fmat(i-1,6)*fmat(2,6)
               fmat(i,7)=fmat(i-1,7)*fmat(2,7)
               fmat(i,8)=fmat(i-1,8)*fmat(2,8)
               fmat(i,9)=fmat(i-1,9)*fmat(2,9)
 4353         continue
              do 67 i=1,nx(1)                                            5d30s00
               vs=vs+ccc(i,1)*fmat(ixp(i,1,1),ixw(i,1,1))*damp               5d30s00
   67         continue                                                  5d30s00
              do 68 i=1,nx(2)                                            5d30s00
               vs=vs+ccc(i,2)*fmat(ixp(i,1,2),ixw(i,1,2))*damp               5d30s00
     $                       *fmat(ixp(i,2,2),ixw(i,2,2))               5d30s00
   68         continue                                                  5d30s00
              do 69 i=1,nx(3)                                            5d30s00
               vs=vs+ccc(i,3)*fmat(ixp(i,1,3),ixw(i,1,3))*damp               5d30s00
     $                       *fmat(ixp(i,2,3),ixw(i,2,3))               5d30s00
     $                       *fmat(ixp(i,3,3),ixw(i,3,3))               5d30s00
   69         continue                                                  5d30s00
              do 70 i=1,nx(4)                                            5d30s00
               vs=vs+ccc(i,4)*fmat(ixp(i,1,4),ixw(i,1,4))*damp               5d30s00
     $                       *fmat(ixp(i,2,4),ixw(i,2,4))               5d30s00
     $                       *fmat(ixp(i,3,4),ixw(i,3,4))               5d30s00
     $                       *fmat(ixp(i,4,4),ixw(i,4,4))               5d30s00
   70         continue                                                  5d30s00
              do 71 i=1,nx(5)                                            5d30s00
               vs=vs+ccc(i,5)*fmat(ixp(i,1,5),ixw(i,1,5))*damp               5d30s00
     $                       *fmat(ixp(i,2,5),ixw(i,2,5))               5d30s00
     $                       *fmat(ixp(i,3,5),ixw(i,3,5))               5d30s00
     $                       *fmat(ixp(i,4,5),ixw(i,4,5))               5d30s00
     $                       *fmat(ixp(i,5,5),ixw(i,5,5))               5d30s00
   71         continue                                                  5d30s00
              do 72 i=1,nx(6)                                            5d30s00
               vs=vs+ccc(i,6)*fmat(ixp(i,1,6),ixw(i,1,6))*damp               5d30s00
     $                       *fmat(ixp(i,2,6),ixw(i,2,6))               5d30s00
     $                       *fmat(ixp(i,3,6),ixw(i,3,6))               5d30s00
     $                       *fmat(ixp(i,4,6),ixw(i,4,6))               5d30s00
     $                       *fmat(ixp(i,5,6),ixw(i,5,6))               5d30s00
     $                       *fmat(ixp(i,6,6),ixw(i,6,6))               5d30s00
   72         continue                                                  5d30s00
              do 73 i=1,nx(7)                                            5d30s00
               vs=vs+ccc(i,7)*fmat(ixp(i,1,7),ixw(i,1,7))*damp               5d30s00
     $                       *fmat(ixp(i,2,7),ixw(i,2,7))               5d30s00
     $                       *fmat(ixp(i,3,7),ixw(i,3,7))               5d30s00
     $                       *fmat(ixp(i,4,7),ixw(i,4,7))               5d30s00
     $                       *fmat(ixp(i,5,7),ixw(i,5,7))               5d30s00
     $                       *fmat(ixp(i,6,7),ixw(i,6,7))               5d30s00
     $                       *fmat(ixp(i,7,7),ixw(i,7,7))               5d30s00
   73         continue                                                  5d30s00
              do 74 i=1,nx(8)                                            5d30s00
               vs=vs+ccc(i,8)*fmat(ixp(i,1,8),ixw(i,1,8))*damp               5d30s00
     $                       *fmat(ixp(i,2,8),ixw(i,2,8))               5d30s00
     $                       *fmat(ixp(i,3,8),ixw(i,3,8))               5d30s00
     $                       *fmat(ixp(i,4,8),ixw(i,4,8))               5d30s00
     $                       *fmat(ixp(i,5,8),ixw(i,5,8))               5d30s00
     $                       *fmat(ixp(i,6,8),ixw(i,6,8))               5d30s00
     $                       *fmat(ixp(i,7,8),ixw(i,7,8))               5d30s00
     $                       *fmat(ixp(i,8,8),ixw(i,8,8))               5d30s00
   74         continue                                                  5d30s00
              vs=vs+const*damp                                          8d3s00
c$$$              write(6,1492)(rij(ig,j),j=1,10),v(ig),vs
c$$$ 1492         format(1p12e15.7)
              v(ig,1)=vs                                                7d25s00
  100 continue
      return
      end
      subroutine setupvp(nx,ixp,ixw,ccc,idx,reau,ex,thrd,exa,srh,sr12,
     $           cs,cb,beta,const,ixa,nexpa,ixr,cands,nexpr,r01,        7d27s00
     $           r02,r03,r04,a02,a03,a04)                               10d6s00
      implicit real*8 (a-h,o-z)
      parameter (id=6400,idu=20,idf=100)
      dimension nx(8),ixp(idx,8,8),ixw(idx,8,8),ccc(idx,8),             7d25s00
     $          ixa(idx,5),ixr(idx,4),ixar(id,2),cands(idx,1)           7d27s00
      dimension cs(6),cb(6),
     $          ix(id,9),cc(id),iperm(4,24),scart(3,4),
     $          iaperm(6,24),ca(6),stry(8,24),rx(4),try(24),ihit(24),
     $          ixs(id,9),ccs(id),jnu(idu),inux(idf,9),cnu(idu,idf),    6d10s00
     $          cnua(idu,idu),cnub(idu,idf),cnuc(idu),ipvt(idu),        6d10s00
     $          cnuu(idu,idf),                                          6d10s00
     $          icnuc(idu),cnuaa(idu,idu),cnud(idu,idf),icnud(idu)      6d10s00
      data iperm/1,2,3,4, 1,2,4,3, 1,3,2,4, 1,3,4,2, 1,4,2,3, 1,4,3,2,  5d8s00
     $           2,1,3,4, 2,1,4,3, 2,3,1,4, 2,3,4,1, 2,4,1,3, 2,4,3,1,  5d8s00
     $           3,1,2,4, 3,1,4,2, 3,2,1,4, 3,2,4,1, 3,4,1,2, 3,4,2,1,  5d8s00
     $           4,1,2,3, 4,1,3,2, 4,2,1,3, 4,2,3,1, 4,3,1,2, 4,3,2,1/  5d8s00
      data iaperm/1,2,3,4,5,6, 1,3,2,5,4,6, 2,1,3,4,6,5, 2,3,1,6,4,5,
     $            3,1,2,5,6,4, 3,2,1,6,5,4, 1,4,5,2,3,6, 1,5,4,3,2,6,
     $            4,1,5,2,6,3, 4,5,1,6,2,3, 5,1,4,3,6,2, 5,4,1,6,3,2,
     $            2,4,6,1,3,5, 2,6,4,3,1,5, 4,2,6,1,5,3, 4,6,2,5,1,3,
     $            6,2,4,3,5,1, 6,4,2,5,3,1, 3,5,6,1,2,4, 3,6,5,2,1,4,
     $            5,3,6,1,4,2, 5,6,3,4,1,2, 6,3,5,2,4,1, 6,5,3,4,2,1/
       read(23,*)re,ex,exa,beta,iradau,lm
c       write(6,*)('re = '),re
c       write(6,*)('stretch morse parameter '),ex
c       write(6,*)('bend extrapolation parameter '),exa
c       write(6,*)('damping parameter = '),beta                          6d15s00
c       write(6,*)('iradau = '),iradau                                   6d19s00
c       write(6,*)('lm = '),lm                                           6d19s00
       tocm=1d0/4.556335d-6                                              5d23s91
       sr12=sqrt(1d0/12d0)
       srh=sqrt(0.5d0)
       nexp=0
       ix1=0
       ix2=0
c       call flush(6)
    2  continue
        read(23,*,end=4)i1,i2,i3,i4,i5,i6,i7,i8,i9,coef
        in=i3+i4+i5+i6+i7+i8+i9
        j1=in+i2
        j2=in+i1
        if(j1.eq.0.and.lm.ne.0)then                                     6d23s00
c         write(6,*)('pure stretch power: '),i1,coef
         cs(i1)=coef*6d0
         ix1=max(ix1,i1)
         coef=0d0                                                       6d13s00
        else if(j2.eq.0.and.lm.ne.0)then                                6d23s00
c         write(6,*)('pure bend power: '),i2,coef
         cb(i2)=coef*4d0
         ix2=max(ix2,i2)
         coef=0d0                                                       6d13s00
        end if                                                          6d13s00
         nexp=nexp+1
         if(nexp.gt.id)stop 'id'
         ix(nexp,1)=i1
         ix(nexp,2)=i2
         ix(nexp,3)=i3
         ix(nexp,4)=i4
         ix(nexp,5)=i5
         ix(nexp,6)=i6
         ix(nexp,7)=i7
         ix(nexp,8)=i8
         ix(nexp,9)=i9
         cc(nexp)=coef
c         write(6,3)nexp,(ix(nexp,j),j=1,9),cc(nexp)
c$$$         cc(nexp)=cc(nexp)*4.556335d-6                                  7d24s00
    3    format(i5,9i2,1pe22.14)
c$$$        end if
        go to 2
    4  continue
       energy=1d0/(1d-18*6.0221367d23*1d-3*3.808798d-4)
       reau=re*1.889725989d0
       do 5 i=1,nexp
c        write(6,3)i,(ix(i,j),j=1,9),cc(i)
c$$$        cc(i)=cc(i)*energy
        do 6 j=1,9
         ix(i,j)=ix(i,j)+1
    6   continue
    5  continue
       thrd=1d0/3d0
       close(unit=23)
c       write(6,*)('analyzing fcns for numbers of multiplications ')
       nx(1)=0                                                             5d30s00
       nx(2)=0                                                             5d30s00
       nx(3)=0                                                             5d30s00
       nx(4)=0                                                             5d30s00
       nx(5)=0
       nx(6)=0
       nx(7)=0
       nx(8)=0
       nxzero=0                                                         7d24s00
       do 600 i=1,nexp                                                  5d30s00
        nnz=0                                                           5d30s00
        do 601 j=1,9                                                    5d30s00
         if(ix(i,j).gt.1)nnz=nnz+1                                      5d30s00
  601   continue                                                        5d30s00
        if(nnz.eq.0)then                                                7d24s00
         nxzero=nxzero+1                                                7d24s00
         const=cc(i)                                                    7d24s00
        else                                                            7d24s00
         if(nnz.eq.7)then
c          write(6,1492)(ix(i,j),j=1,9)
 1492     format(9i2)
         end if
        nx(nnz)=nx(nnz)+1                                                 5d30s00
        if(nx(nnz).gt.idx)stop 'idx'                                    7d24s00
        ccc(nx(nnz),nnz)=cc(i)                                           5d30s00
        ii=1                                                            5d30s00
        do 602 j=1,9                                                    5d30s00
         if(ix(i,j).gt.1)then                                           5d30s00
          ixp(nx(nnz),ii,nnz)=ix(i,j)                                    5d30s00
          ixw(nx(nnz),ii,nnz)=j                                          5d30s00
          ii=ii+1                                                       5d30s00
         end if                                                         5d30s00
  602   continue                                                        5d30s00
        end if                                                          7d24s00
  600  continue                                                         5d30s00
c       write(6,*)('number of fcns with zero nonzero: '),nxzero          7d24s00
c       write(6,*)('number of fcns with one nonzero: '),nx(1)
c       write(6,*)('number of fcns with two nonzero: '),nx(2)
c       write(6,*)('number of fcns with three nonzero: '),nx(3)
c       write(6,*)('number of fcns with four nonzero: '),nx(4)
c       write(6,*)('number of fcns with five nonzero: '),nx(5)
c       write(6,*)('number of fcns with six nonzero: '),nx(6)
c       write(6,*)('number of fcns with seven nonzero: '),nx(7)
c       write(6,*)('number of fcns with eight nonzero: '),nx(8)
       return
       end

       subroutine potinich4
       implicit none
       character*80 title
       real*8 tocm,bohr,bohri
       common /pes1/  tocm,bohr,bohri
       real*8 re0,ex0,beta0,r01,r02,r03,r04,a02,a03,a04,a06,a08
       common /pes2/re0,ex0,beta0,r01,r02,r03,r04,a02,a03,a04,a06,a08
       real*8  ccc
       common /pes3/ ccc
       real*8 reau,ex,thrd,exa,srh,sr12,cs,cb,beta,const
       common /pes4/ reau,ex,thrd,exa,srh,sr12,cs,cb,beta,const
       real*8 cands
	 integer nexpr
       common /pes5/ cands,nexpr
       integer nexpa,ixa,ixr,nx,ixp,ixw,ixam
       common /pesint/ nexpa,ixa,ixr,nx,ixp,ixw
       integer id,n
       parameter (n=1)
       real*8 rij,v,fmat,rmat
      parameter (id=2000)
      dimension rij(n,10),v(n,1),cs(6),cb(6),                           7d10s01
     $          fmat(9,9),nx(8),ixp(id,8,8),ixw(id,8,8),ccc(id,8)       5d30s00
      dimension ixa(id,5),ixr(id,4),cands(id,id),rmat(id),              7d27s00
     $          ixam(id,5)                                              7d27s00

      integer intmul,maxbc
      real*8 tovr
      common/iszcm/intmul                                               12d8s88
      common/timrcm/tovr
      parameter (maxbc=1 000 000)
      integer ibcoff,ibcmax,ibca,mnbcv,ibc
      real*8 bc(maxbc)
      common ibcoff,ibcmax,ibca,mnbcv,bc
      dimension ibc(maxbc)
      equivalence (ibc,bc)
      integer i,j

      call koordreadfch4

      tocm=1d0/4.556335d-6                                              5d23s91
      bohr=0.529177249d0
      bohri=1d0/bohr
c       write(6,111)
  111  format(/1x,'pes for ch4: fast version')
       open(unit=23,file='../data/CHD3/pes',
     . form='formatted',status='old')
      read(23,1492)title                                                 8d1s00
 1492 format(a80)                                                       8d1s00
c      write(6,*)('title for PES:')                                      8d1s00
c      write(6,1492)title                                                8d1s00
      read(23,*)re0,ex0,beta0
      read(23,*)r01,r02,r03,r04
      read(23,*)a02,a03,a04,a06,a08                                      8d3s00
      re0=re0*bohri                                                     8d1s00
      r01=r01*6d0
      r02=r02*6d0
      r03=r03*6d0
      r04=r04*6d0
      a02=a02*4d0
      a03=a03*4d0
      a04=a04*4d0
c      write(6,*)('zero order fit parameters: ')
c      write(6,*)('re0,ex0,beta0 '),re0,ex0,beta0
c      write(6,*)('r0i '),r01,r02,r03,r04
c      write(6,*)('a0i '),a02,a03,a04,a06,a08                            8d3s00
c      call flush(6)
       call setupvp(nx,ixp,ixw,ccc,id,reau,ex,thrd,exa,srh,sr12,cs,cb,
     $              beta,const,ixa,nexpa,ixr,cands,nexpr,r01,r02,       7d27s00
     $              r03,r04,a02,a03,a04)                                10d6s00
!       call setupvp(nx,ixp,ixw,ccc,id,reau,ex,thrd,exa,srh,sr12,cs,cb,
!     $              beta,const,ixa,nexpa,ixr,cands,nexpr,r01,r02,       7d27s00
!     $              r03,r04,a02,a03,a04,a06,a08)                        10d6s00
c       write(6,*)('in vibpot, ex= '),ex
       do 510 j=1,5                                                     7d27s00
        do 1676 i=1,nexpa                                               7d27s00
         ixam(i,j)=ixa(i,j)-1
 1676   continue                                                        7d27s00
  510  continue                                                         7d27s00
       close (23)
       end

c -------------------------------------------------------------------
      subroutine koordreadfch4
      implicit none

      integer    na,na3,dimch3
      parameter (na=5,na3=3*na,dimch3=12)

      real*8     mch3,hmass
      parameter (hmass=1836.109d0)

      integer i,j

      real*8 work(4,4),m(na),sm1(na),alpha(4)
      common /rctrans/ work,m,sm1

c masses
      m(1)=11.907d0
      m(2)=1.9984d0
      m(3)=1.9984d0
      m(4)=1.9984d0
      m(5)=1.0d0
!      m(1)=11.907d0
!      m(2)=1.0d0
!      m(3)=1.0d0
!      m(4)=1.0d0
!      m(5)=1.0d0
      do i=1,na
         m(i)=m(i)*hmass
         sm1(i)=1.0d0/sqrt(m(i))
      end do

      mch3=0.0d0
      do i=1,4
         mch3=mch3+m(i)
      end do

c calculate matrix for transformation between mass weighted
c cartesians and Radaus

      do i=1,4
         alpha(i)=sqrt(m(i)/mch3)
      end do

      do i=1,4
         work(i,1)=alpha(i)
         work(1,i)=alpha(i)
      end do

      do i=2,4
         do j=2,4
            work(j,i)=alpha(j)*alpha(i)/(alpha(1)+1.0d0)
         end do
      end do

      do i=2,4
         work(i,i)=work(i,i)-1.0d0
      end do

      end
c ------------------------------------------------------------------------

      subroutine IntToCartfch4(q,x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine transforms the internal coordinates       c
c to cartesian ones.                                        c
c Subroutine koordread must be called in main program       c
c GS / 9 Feb 2006                                           c
c rst, RST                                                  c
c GS / 25 July 2008                                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer    na,na3,dimq,dimch3
      parameter  (na=5,na3=3*na,dimq=na3-6,dimch3=12)

      real*8  q(dimq),x(na3)
      real*8  r(3),xcm(3)
      real*8  radau(dimch3)
      real*8  mch3,mch4
      real*8     hmass
      parameter (hmass=1836.109d0)

      real*8   rmch3,rmch4,srmch31,srmch41
      real*8   dist1,dist2,nenn1,nenn2

      integer i,j,k

      real*8 rctrans(4,4),m(na),sm1(na)
      common /rctrans/ rctrans,m,sm1

c mass of CH3-group
      mch3=0.0d0
      do i=1,4
         mch3=mch3+m(i)
      end do

c mass of CH4-group
      mch4=0.0d0
      do i=1,5
         mch4=mch4+m(i)
      end do

c reduced mass between CH3 and 4th H

      rmch3=mch3*m(5)/(mch3+m(5))
      srmch31=1.0d0/sqrt(rmch3)


c reduced mass between CH4 and 5th H

c      rmch4=mch4*m(6)/(mch4+m(6))
c      srmch41=1.0d0/sqrt(rmch4)


c order of q's:
c  1 rho          4 theta     7 r           10 R
c  2 theta_rho    5 phi       8 s           11 S
c  3 phi_rho      6 chi       9 t           12 T

c order of radau's
c  1 center of mass    4 x H_1     7 x H_2   10 x H_3
c  2 "                 5 y H_1     8 y H_2   11 y H_3
c  3 "                 6 z H_1     9 z H_2   12 z H_3

c order of cartesians
c 1-3 C (x,y,z)
c 4-6, 7-9, 10-12 methyl-H
c 13-15, 16-18 H2-H

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

c     print*,"Radau in Ruecktransformation"
c     do i=1,4
c        print'(3F12.6)',(radau(j+(i-1)*3),j=1,3)
c     end do

c Third step: calculate cartesians for methyl group
c             from Radau's
      do i=1,dimch3,3
         x(i  )=0.0d0
         x(i+1)=0.0d0
         x(i+2)=0.0d0
         do j=1,dimch3,3
           x(i  )=x(i  )+rctrans(i/3+1,j/3+1)*radau(j  )
           x(i+1)=x(i+1)+rctrans(i/3+1,j/3+1)*radau(j+1)
           x(i+2)=x(i+2)+rctrans(i/3+1,j/3+1)*radau(j+2)
         end do
      end do

c  unmassweight cartesians, correct sign
      do i=1,dimch3,3
         x(i  )=-x(i  )*sm1(i/3+1)
         x(i+1)=-x(i+1)*sm1(i/3+1)
         x(i+2)=-x(i+2)*sm1(i/3+1)
      end do

c Fourth step: add 4th hydrogen atom

c  unmassweight r
      dist1=q(7)*srmch31
      nenn1=1.0d0/(1.0d0+q(8)**2+q(9)**2)

c  calculate cartesians
      x(13)=2*dist1*q(8)*nenn1
      x(14)=2*dist1*q(9)*nenn1
      x(15)=dist1*(1.0d0-q(8)**2-q(9)**2)*nenn1

      end
c ------------------------------------------------------------------------
