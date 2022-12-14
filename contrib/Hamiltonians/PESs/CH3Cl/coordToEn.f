			subroutine potentialch3cl(ex,w,part)
			implicit none

			double precision ex(15)		!cartesian coordiantes
			double precision in(9)		!internal coordinates
			double precision en(3,3)	  !energy
			double precision w

			double precision p(1000)	!parameter
c			integer pst(2,400)				!pointers for parameters

			integer npar							!number of parameters
			integer part ! which element of the diabatic matrix

			logical skip							!dummy

			double precision s(9),t(9) !transformed coordiantes

			double precision cmtoau

	            double precision ev(3,3), ew(3)
			parameter (cmtoau=4.55633d-6)

			include 'common.incl'

			skip=.false.
			include 'fort.66'

			p(2) = p(2) - p(1)
			p(1) = 0.

c     calculate number of parameters
			npar=pst(1,304)+pst(2,304)-1

c     calculate internal coordinates
c			call extToIntCoord(ex,in)
			call katToIntCoord(ex,in)
c     transform coordinates
			call ctrans(in,s,t,1,skip,p,npar)

c     calculate monomials
			call vwzprec(1,s)

c     get energy of groundstate
			call calculatePotential(en,in,p,npar,pst,skip)
c test diagonalizing
!      	call rs(3,3,en,ew,1,ev)
!		w = ew(part+1)

c			if(en(1,1)/cmtoau.lt.-600000.d0) then
c				 write(6,'(''--------------------------------------'')')
c				 write(6,'(''found hole'')')
c				 write(6,'(''Energy is:''f12.4)') en/cmtoau
c				 write(6,*)
c				 write(6,'(''catisian Coordinates are:'')')
c				 write(6,'(3f9.4)') ex
c				 write(6,*)
c				 write(6,'(''primitive Coordinates are:'')')
c				 write(6,'(3f9.4)') in
c				 write(6,*)
c				 write(6,'(''transformed Coordinates are:'')')
c				 write(6,'(3f9.4)') s
c				 write(6,*)
c				 write(6,'(''kill calculation'')')
c				 stop
c			endif

       if (part.eq.0) then
       	w = en(1,1)
      else if (part.eq.1) then
       	w = en(2,1)
      else if (part.eq.2) then
       	w = en(3,1)
      else if (part.eq.4) then
       	w = en(2,2)
      else if (part.eq.5) then
       	w = en(3,2)
      else if (part.eq.8) then
       	w = en(3,3)
      else
            stop "gar nicht gut : ("
      endif
!			if (part.eq.0) w = en(1,1)
!			if (part.eq.1) w = en(2,2)
!			if (part.eq.2) w = en(3,3)
!			if (part.eq.3) w = en(1,2)
!			if (part.eq.4) w = en(1,3)
!			if (part.eq.5) w = en(2,3)
			end subroutine

			
