C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% FUNCTION EXTPOT(R,COSTH)
C% Author: J. M. M. Howson and J. D. Hutson 2001
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION EXTPOT(R,COSTH)
c
c     Program to project out in legendre polys
c     in the angular co-ordinate.
c
      implicit double precision (a-h,o-z)
      integer mxparm, mxdata, mxsite
      save
      
      PARAMETER (MXSITE=10)
      parameter (mxparm=500,mxdata=500)
      
      double precision parm0(mxparm), rdata(mxdata),RLAM(50),BLAM(50)
      double precision p0(mxparm),vx(10),ax(10),sumleg
      double precision alpha(19,13),x(19),rad_kernel,rnw,cthnw
      double precision x_k,r,costh,pi,plm,pot,vpp,ra,ro,rs,rc,rb
      double precision eps_scal, r_scal, extpot, extint,ang,aur,aue
      
      integer n_r,n_lam, np,k
      integer i,j,n,npts
c      common npts
C     
      DATA RLAM/50*0.D0/, BLAM/50*0.D0/
      DATA AUE/219474.6354D0/
      DATA AUR/0.52917706D0/
C     
      ro = 1.68029d0
      rs = 1.03731d0

      ra = (r*r+rs*rs+2.0d0*r*rs*costh)**0.5d0
      rb = (r*r+ro*ro-2.0d0*r*ro*costh)**0.5d0

      rnw = (ra + rb) /(ro+rs)
      cthnw = (ra - rb) /(ro+rs)

      eps_scal = sumleg(vx,np,cthnw)
      
      r_scal = rnw*sumleg(ax,np,cthnw)
c     write(*,*)'r=',r_scal,'costh=',costh,'eps_scal=',eps_scal
      
      pi=dacos(-1.d0)
      
      pot=0.0d0

      do j = 1,n_lam
         
         vpp = 0.0d0
         do i=1,npts 
            x_k=x(i)
            vpp=vpp+alpha(i,j)*rad_kernel(x_k,r_scal)
         enddo

         pot = pot+ vpp*plm(j-1,0,cthnw)
c         pot = pot+ vpp*legp(j-1,cthnw)
      enddo
*
      extpot=eps_scal*pot
*
      return

      ENTRY EXTINT
c     read number of scaling paras
      read(5,*) np
c     read scaling paras
      do i=1,np
	 read(5,*) vx(i), ax(i)
      enddo
c     read number of alphas(=no. of rvalues and theta values)
      read(5,*) npts, n_lam
c     NB:x(i)=radial pt at which have abinitio pt
c     NB:y(i)=(1-cos(theta))/2.0
c     where theta is angular pt at which have abinitio pt
c     read x,y and alphas for  RK
      
      do j = 1,n_lam
         do i=1, npts
            read(5,*) x(i),alpha(i,j)
         enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% FUNCTION RAD_KERNEL(X,R)
C% Author: J. M. M. Howson and J. D. Hutson 2001
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function rad_kernel(x,r)
      implicit none
      real*8 rad_kernel,x,r,x_l,x_s,R0
      
      x_l=r
      x_s=x
      if(r.lt.x) then
         x_s=r
         x_l=x
      endif
c     kernel for protonated systems (long range prop to R-4)
c     rad_kernel=2.d0/(15.d0*x_l**5)*(1.d0-5.d0*x_s/(7.d0*x_l))
c     kernel for protonated systems (long range prop to R-6)m=5,n=2
c     rad_kernel=2.d0/(21.d0*x_l**6.0d0)*(1.d0-3.d0*x_s/(4.d0*x_l))
c     m=6, n=2
      rad_kernel=1.d0/(14.d0*x_l**7.0d0)*(1.d0-7.d0*x_s/(9.d0*x_l))
c     Increase order of smoothness.n=3
c     rad_kernel=1.0d0/(28.0d0*x_l**7.0d0)*(1.d0-7.d0*x_s/(5.d0*x_l)
c     . + 28.d0*x_s*x_s/(55.d0*x_l*x_l))
c     n=3 m=5
c     rad_kernel=3.d0/(56.0d0*x_l**6)*(1.d0-4.d0*x_s/(3.d0*x_l)
c     . + 7.d0*x_s*x_s/(15.d0*x_l*x_l))
      
c     rad_kernel=1.d0/(14.0d0*(R0+x_l)**7)*
c     .    (1.d0-7.d0*(x_s+R0)/(9.d0*(x_l+R0)))
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


