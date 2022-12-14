C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE CTRANS(...)
C% 
C% A. Viel 21 09 2010
C%
C% Modified by T. Weike 03.2016
C% 
C% Routine to transform symmetryinput coordinates  to symmetrized 
C% coordinates. Distances Are discribet by Morse coordinates or
C% TMC depending on Set Parameters in the Genetic Input.
C% 
C% Input variables:
C% q:    primitive coordinates (double[qn])
C%       q(1): CCl-distance
C%       q(2): CCl-angel-in x direction
C%       q(3): CCl-angle-in y direction
C%       q(4): CH1-distance
C%       q(5): CH2-distance
C%       q(6): CH3-distance
C%       q(7): umbrella angle
C%       q(8): CH2-bending angle (direction away from yaxis)
C%       q(9): CH3-bending angle (direction to yaxis)
C% symm: 0 no transformation, 1 symmetrize
C% t:    dummy (double[qn])
C% skip: dummy (logical)
C% p:    parameter vector
C% npar: length of parameter vector
C% 
C% Output variables
C% s: symmetrized coordinates (double[qn])
C%    s(1): CCl-Streck, 2:
C%    s(2): CH-symetric streatch
C%    s(3): CH-umbrella
C%    s(4): CCl-angle-ex
C%    s(5): CCl-angle-ey
C%    s(6): CH-asymetric streatch-ex
C%    s(7): CH-asymetric streatch-ey
C%    s(8): CH-bend-ex
C%    s(9): CH-bend-ey
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine ctrans(q,s,t,symm,skip,p,npar)
      implicit none

      include 'common.incl'
      include 'states.incl'

c     Funktions
      double precision func_f !morse exponent
      double precision damp

c     Variables
      integer i !running indices

      integer symm !trasform or not

      double precision q(qn) !given coordinates
      double precision s(qn) !output coordinates
      double precision t(qn) !dummy
      
      integer npar             !length of parameter vector
      double precision p(npar) !parameter vector

      double precision ang2rad !conversion factor
      parameter (ang2rad=0.017453292519943295769236907684886127134428718
     &     885417254560d0)

      double precision pi  !pi
      double precision pi2 !2*pi
      parameter (pi=3.14159265358979323846264338327950288419716939937510
     &     5820974d0)
      parameter (pi2=6.2831853071795864769252867665590057683943387987502
     &     11641949d0)

      double precision dchequi !ch-equilibrium distance
      parameter (dchequi=1.08485600d0)
      double precision dcclequi !ccl-equilibrium distance
      parameter (dcclequi=1.80969700d0)
      double precision umbequi !umbralla equilibrium angle
      parameter (umbequi=15.012025d0)

      double precision sq2,sq3,sq6 !precalculatet roots

      !primitive coordinates
      double precision dccl           !ccl-distance
      double precision axcl,aycl      !cl-angeles (x,y)
      double precision dch1,dch2,dch3 !ch-distances
      double precision umb            !Umbrella Angle from xy-plane 
      double precision ah2,ah3        !Ang. from neg. x-Achsis to h2/h3

      !tuned distances
      double precision tmcccl                 !C-Cl-Distances
      double precision tmcch1, tmcch2, tmcch3 !C-H-Distances

      !scaledangles
      double precision alpha1,alpha2,alpha3 !HCH-Angeles
      double precision theta                !Umbrella
      double precision betax,betay          !Cl-Angles

      !Symmetry coordinates
      double precision  aCl, aR, aAng !a1-modes Cl, H-Dist., H-Ang.
      double precision exCl,exR,exAng !ex components Cl, H-Dist., H-Ang.
      double precision eyCl,eyR,eyAng !ey components Cl, H-Dist., H-Ang.

      logical skip,dbg

      dbg=.true.
      dbg=.false.

c     subpressing warnings
      t(1)=0.d0

      sq2=1.d0/sqrt(2.d0)
      sq3=1.d0/sqrt(3.d0)
      sq6=1.d0/sqrt(6.d0)


c     do no transformation if symm=0:
      if (symm.eq.0) then
         do i=1,qn
            s(i) = q(i)
         enddo
         return
      endif 

c     get the coordinates for better readebility
      dccl = q(1)
      axcl = q(2)*ang2rad
      aycl =-q(3)*ang2rad
      dch1 = q(4)
      dch2 = q(5)               !Warning: Changed -> mathematical positiv
      dch3 = q(6)
      umb  =(q(7)-umbequi)*ang2rad
      ah2  =(-q(8))*ang2rad  !Warning: Changed -> mathematical positiv
      ah3  =( q(9))*ang2rad

ct      dccl = q(1)
ct      axcl = q(2)
ct      aycl =-q(3)
ct      dch1 = q(4)
ct      dch2 = q(5)         !Warning: Changed -> mathematical positiv
ct      dch3 = q(6)
ct      umb  =(q(7)-umbequi)
ct      ah2  =(120+q(8))  !Warning:Changed -> mathematical positiv
ct      ah3  =(120+q(9))

c     transforming distances to tunable morse coordiantes
      tmcccl = 1.d0 - dexp(-func_f(dccl,p(pst(1,2)),1)*(dccl-dcclequi)) 
      tmcch1 = 1.d0 - dexp(-func_f(dch1,p(pst(1,3)),2)*(dch1-dchequi))
      tmcch2 = 1.d0 - dexp(-func_f(dch2,p(pst(1,3)),2)*(dch2-dchequi))
      tmcch3 = 1.d0 - dexp(-func_f(dch3,p(pst(1,3)),2)*(dch3-dchequi))
c      tmcch1 = dch1
c      tmcch2 = dch2
c      tmcch3 = dch3

      !call switch(damp,-tmcccl,0.3d0,2.d0)
      !tmcccl=damp*tmcccl-1.9d0*(1.d0-damp)
      
c     transformation of angles for better representation of dissociation
      betax  = axcl/dccl
      betay  = aycl/dccl
c      alpha1 =(pi2-ah2-ah3)-120!/(dch2*dch3)
      alpha1 =(pi2-ah2-ah3)-120!/(dch2*dch3)
      alpha2 = ah3-120!/(dch1*dch3)  !Change: the angel is now on the 
                                     !opsite side of the atom
      alpha3 = ah2-120!/(dch1*dch2)
      theta  = umb!/(dch1*dch2*dch3) !Does not help for dissociation
      
      if(dbg) write(6,*) 'given coordinates'
      if(dbg) write(6,*) 'dccl: ', dccl
      if(dbg) write(6,*) 'axcl: ', axcl
      if(dbg) write(6,*) 'aycl: ', aycl
      if(dbg) write(6,*) 'dch1: ', dch1
      if(dbg) write(6,*) 'dch2: ', dch2
      if(dbg) write(6,*) 'dch3: ', dch3
      if(dbg) write(6,*) 'umb : ', umb 
      if(dbg) write(6,*) 'ah2 : ', ah2 
      if(dbg) write(6,*) 'ah3 : ', ah3 

      if(dbg) write(6,*) 'tmcs'
      if(dbg) write(6,*) 'tmcccl: ', tmcccl
      if(dbg) write(6,*) 'tmcch1: ', tmcch1
      if(dbg) write(6,*) 'tmcch2: ', tmcch2
      if(dbg) write(6,*) 'tmcch3: ', tmcch3

      if(dbg) write(6,*) 'scaled angles'
      if(dbg) write(6,*) 'betax : ', betax 
      if(dbg) write(6,*) 'betay : ', betay 
      if(dbg) write(6,*) 'alpha1: ', alpha1
      if(dbg) write(6,*) 'alpha2: ', alpha2
      if(dbg) write(6,*) 'alpha3: ', alpha3

c     symmetriesation of coordinates
CT       aCl  = tmcccl
CT      exCl  = betax
CT      eyCl  =-betay
CT       aR   = sq3*(tmcch1+tmcch2+tmcch3)
CT      exR   = sq6*(2*tmcch1-tmcch2-tmcch3)
CT      eyR   = sq2*(tmcch2-tmcch3)
CT       aAng = theta
CT      eyAng = sq6*(2*alpha1-alpha2-alpha3)
CT      exAng =-sq2*(alpha2-alpha3) 

       aCl  = tmcccl
      exCl  = betax
      eyCl  =-betay
       aR   = sq3*(tmcch1+tmcch2+tmcch3)
      exR   = sq6*(2*tmcch1-tmcch2-tmcch3)
      eyR   =-sq2*(tmcch2-tmcch3)
       aAng = theta
      exAng =-sq6*(2*alpha1-alpha2-alpha3)
      eyAng = sq2*(alpha2-alpha3) 
c      exAng = sq2*(alpha2+alpha3)
c      eyAng = sq2*(alpha2-alpha3) 


      if(dbg) write(6,*) 'symmetrized coordinates'
      if(dbg) write(6,*) ' aCl : ',  aCl 
      if(dbg) write(6,*) 'exCl : ', exCl 
      if(dbg) write(6,*) 'eyCl : ', eyCl 
      if(dbg) write(6,*) ' aR  : ',  aR  
      if(dbg) write(6,*) 'exR  : ', exR  
      if(dbg) write(6,*) 'eyR  : ', eyR  
      if(dbg) write(6,*) ' aAng: ',  aAng
      if(dbg) write(6,*) 'exAng: ', exAng
      if(dbg) write(6,*) 'eyAng: ', eyAng

      s(1) =  aCl
      s(2) =  aR
      s(3) =  aAng
      s(4) = exCl
      s(5) = eyCl
      s(6) = exR
      s(7) = eyR
      s(8) = exAng
      s(9) = eyAng

      end


C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% FUNCTION F(...)
C% 
C% Returns exponent of tunable Morse coordinate
C% exponent is polynomial * gaussian (skewed)
C% ilabel = 1 or 2 selects the parameters a and sfac to be used
C% 
C% Background: better representation of the prefector in the
C%             exponend of the morse function.
C% Formular: f(r) = latest no3 paper
C%
C% Variables:
C% x:  distance of atoms (double)
C% p:  parameter vector (double[20])
C% ii: 1 for CCl and 2 for CCH (int)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function func_f(x,p,ii)
      implicit none

      integer ii !1 for CCL and 2 for CCH 

      integer i !running index

      double precision x !coordinate

      double precision r    !equilibrium distance
      double precision gaus !gaus part of f
      double precision poly !polynom part of f
      double precision skew !tanh part of f

      double precision func_f !prefactor of exponent and returned value

      integer npoly(2) !order of polynom

      double precision p(11)   !parameter-vector

c     Maximum polynom order
      npoly(1)=5
      npoly(2)=5

! p(1):        position_ of equilibrium
! p(2):        constant of exponent
! p(3):        constant for skewing the gaussian
! p(4):        tuning for skewing the gaussian
! p(5):        Gaussian exponent
! p(6):        Shift of Gaussian maximum
! p(7)...:     polynomial coefficients
! p(8+n)...:   coefficients of Morse Power series

!1-exp{[p(2)+exp{-p(5)[x-p(6)]^2}[Taylor{p(7+n)}(x-p(6))]][x-p(1)]}

! Tunable Morse function
! Power series in Tunable Morse coordinates of order m
! exponent is polynomial of order npoly * gaussian + switching function

c     set r  r-r_e
      r=x
      r=r-p(1)

      p(2)=abs(p(2))
      p(5)=abs(p(5))

c     set up_ skewing function:
      skew=0.5d0 * p(3) *(tanh(abs(p(4))*(r-p(6))) + 1.d0)

c     set up_ gaussian function:
      gaus=exp(-p(5)*(r-p(6))**2)

c     set up_ power series:
      poly=0.d0
      do i=0,npoly(ii)-1
        poly=poly+p(7+i)*(r-p(6))**i
      enddo

c     set up_ full exponent function:
      func_f=p(2)  + skew + gaus*poly

      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE EXTTOINTCOORD(EXT,INT)
C%
C% Calculates the primitive internal coordinates from cartisian
C% external coordinates.
C%
C% Input variable:
C% ext: cartisian external coordinates (double[15])
C%      ext( 1): Cl-x
C%      ext( 2): Cl-y
C%      ext( 3): Cl-z
C%      ext( 4): H1-x
C%      ext( 5): H1-y
C%      ext( 6): H1-z
C%      ext( 7): H2-x
C%      ext( 8): H2-y
C%      ext( 9): H2-z
C%      ext(10): H3-x
C%      ext(11): H3-y
C%      ext(12): H3-z
C%      ext(13): C-x
C%      ext(14): C-y
C%      ext(15): C-z 
C%
C% Output variable:
C% in:  primitive internal coordinates (double[qn])
C%      int(1): CCl-distance
C%      int(2): CCl-angel-in x direction
C%      int(3): CCl-angle-in y direction
C%      int(4): CH1-distance
C%      int(5): CH2-distance
C%      int(6): CH3-distance
C%      int(7): umbrella angle
C%      int(8): CH2-bending angle (direction away from yaxis)
C%      int(9): CH3-bending angle (direction to yaxis)       
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine extToIntCoord(ext,in)
      implicit none

      include 'states.incl'

c     functions
      double precision asin
      double precision atan2

c     variables
      integer          exn      !number of external coordinates
      parameter (exn=15)
      double precision ext(exn) !external cartisian coordinates
      double precision in(qn)  !internal coordinates

      double precision ang2rad !conversion factor
      parameter (ang2rad=0.017453292519943295769236907684886127134428718
     &     885417254560)

      logical dbg !more Output for debugging

      dbg=.true.
      dbg=.false.

      if(dbg) then
         write(6,'('' 1: Cl-x'')')
         write(6,'('' 2: Cl-y'')')
         write(6,'('' 3: Cl-z'')')
         write(6,'('' 4: H1-x'')')
         write(6,'('' 5: H1-y'')')
         write(6,'('' 6: H1-z'')')
         write(6,'('' 7: H2-x'')')
         write(6,'('' 8: H2-y'')')
         write(6,'('' 9: H2-z'')')
         write(6,'(''10: H3-x'')')
         write(6,'(''11: H3-y'')')
         write(6,'(''12: H3-z'')')
         write(6,'(''13: C-x '')')
         write(6,'(''14: C-y '')')
         write(6,'(''15: C-z '')')
      endif
      if(dbg) write(6,'(''External coordinates:'')')
      if(dbg) write(6,'(3f9.5)') ext

c     Setting C to the origin
      ext( 1) = ext( 1)-ext(13)
      ext( 2) = ext( 2)-ext(14)
      ext( 3) = ext( 3)-ext(15)

      ext( 4) = ext( 4)-ext(13)
      ext( 5) = ext( 5)-ext(14)
      ext( 6) = ext( 6)-ext(15)

      ext( 7) = ext( 7)-ext(13)
      ext( 8) = ext( 8)-ext(14)
      ext( 9) = ext( 9)-ext(15)

      ext(10) = ext(10)-ext(13)
      ext(11) = ext(11)-ext(14)
      ext(12) = ext(12)-ext(15)

c     calculatin distances
      call dista(in(1),ext( 1)) !C-Cl
      call dista(in(4),ext( 4)) !C-H1
      call dista(in(5),ext( 7)) !C-H2
      call dista(in(6),ext(10)) !C-H3

c     calculate cl-angles
      in(2)=asin(ext(1)/in(1))/ang2rad
      in(3)=asin(ext(2)/in(1))/ang2rad

c     calculate umbrella-angle
      in(7)=-asin(ext(12)/in(6))/ang2rad

c     calculate CH-Bending-angles
      in(8)=60.d0-atan2(ext(8),ext(7))/ang2rad-180.d0
      in(9)=120.d0-atan2(ext(11),ext(10))/ang2rad
      if(dbg) then
         write(6,'(''1: CCl-distance'')')
         write(6,'(''2: CCl-angel-in x direction'')')
         write(6,'(''3: CCl-angle-in y direction'')')
         write(6,'(''4: CH1-distance'')')
         write(6,'(''5: CH2-distance'')')
         write(6,'(''6: CH3-distance'')')
         write(6,'(''7: umbrella angle'')')
         write(6,'(''8: CH2-bending angle'')')
         write(6,'(''9: CH3-bending angle'')')
      endif
      if(dbg) write(6,'(''Internal coordinates:'')')
      if(dbg) write(6,'(3f9.5)') in

      end subroutine

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE DIST(DIST,X1,X2)
C%
C% Calculates the distance to the origin.
C% 
C% Input variables:
C% x1: Point one (double[3])
C% x2: Point two (double[3])
C%
C% Output variable:
C% dist: distance between the points
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine dista(d,x1)
      implicit none
      
      double precision d !distance

      double precision x1(3) !point one

      integer i !running index

      d=0.d0
      do i=1,3
         d=d+x1(i)**2
      enddo
      d=dsqrt(d)

      end subroutine

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE EXTTOINTCOORD(EXT,INT)
C%
C% Calculates the primitive internal coordinates from cartisian
C% external coordinates.
C%
C% Input variable:
C% ext: cartisian external coordinates (double[15])
C%      ext( 1): Cl-x
C%      ext( 2): Cl-y
C%      ext( 3): Cl-z
C%      ext( 4): H1-x
C%      ext( 5): H1-y
C%      ext( 6): H1-z
C%      ext( 7): H2-x
C%      ext( 8): H2-y
C%      ext( 9): H2-z
C%      ext(10): H3-x
C%      ext(11): H3-y
C%      ext(12): H3-z
C%      ext(13): C-x
C%      ext(14): C-y
C%      ext(15): C-z 
C%
C% Output variable:
C% in:  primitive internal coordinates (double[qn])
C%      int(1): CCl-distance
C%      int(2): CCl-angel-in x direction
C%      int(3): CCl-angle-in y direction
C%      int(4): CH1-distance
C%      int(5): CH2-distance
C%      int(6): CH3-distance
C%      int(7): umbrella angle
C%      int(8): CH2-bending angle (direction away from yaxis)
C%      int(9): CH3-bending angle (direction to yaxis)       
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine katToIntCoord(ext,in)

c     coordinate vectors
      double precision ext(15)
      double precision in(9)

c     kartisian koordinates
      double precision IX
      double precision IY
      double precision IZ
      double precision H1X
      double precision H1Y
      double precision H1Z
      double precision H2X
      double precision H2Y
      double precision H2Z
      double precision H3X
      double precision H3Y
      double precision H3Z

c     internal koordinates
      double precision DCCL
      double precision AXCL
      double precision AYCL
      double precision DCH1
      double precision DCH2
      double precision DCH3
      double precision UMB
      double precision AH2
      double precision AH3

c     normation of dist. for projection in plane
      double precision alpha

      integer j,i

      double precision ang2rad !conversion factor
      parameter (ang2rad=0.017453292519943295769236907684886127134428718
     &     885417254560d0)

      IX  = ext( 1) - ext(13)
      IY  = ext( 2) - ext(14)
      IZ  = ext( 3) - ext(15)
      H1X = ext( 4) - ext(13) 
      H1Y = ext( 5) - ext(14) 
      H1Z = ext( 6) - ext(15) 
      H2X = ext( 7) - ext(13) 
      H2Y = ext( 8) - ext(14) 
      H2Z = ext( 9) - ext(15) 
      H3X = ext(10) - ext(13) 
      H3Y = ext(11) - ext(14) 
      H3Z = ext(12) - ext(15)

c     calculate internal coordinates
c     distances
      DCCL=sqrt(IX**2+IY**2+IZ**2)
      DCH1=sqrt(H1X**2+H1Y**2+H1Z**2)
      DCH2=sqrt(H2X**2+H2Y**2+H2Z**2)
      DCH3=sqrt(H3X**2+H3Y**2+H3Z**2)

c     normation
      alpha=atan2(H1Y,H1X)
      H1X= H1X*cos(alpha)+H1Y*sin(alpha)
      H1Y=-H1X*sin(alpha)+H1Y*cos(alpha)
      H2X= H2X*cos(alpha)+H2Y*sin(alpha)
      H2Y=-H2X*sin(alpha)+H2Y*cos(alpha)
      H3X= H3X*cos(alpha)+H3Y*sin(alpha)
      H3Y=-H3X*sin(alpha)+H3Y*cos(alpha)
      IX = IX *cos(alpha)+IY *sin(alpha)
      IY =-IX *sin(alpha)+IY *cos(alpha)
      
c     umbrella angle
      UMB=atan2(-H1Z,H1X)/ang2rad

c     ch bending angles
      AH2=atan2(H2Y,H2X)/ang2rad
      AH3=atan2(H3Y,H3X)/ang2rad

c     cl angles
      AXCL=atan2(IX,IZ)/ang2rad
      AYCL=atan2(IY,IZ)/ang2rad

      in(1) = DCCL
      in(2) = AXCL
      in(3) = AYCL
      in(4) = DCH1
      in(5) = DCH2
      in(6) = DCH3
      in(7) = UMB
      in(8) = AH2
      in(9) = AH3

      end subroutine
