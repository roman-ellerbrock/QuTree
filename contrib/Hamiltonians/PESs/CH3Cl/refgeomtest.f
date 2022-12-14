      program testScan
      implicit none

      double precision qstart(9), qend(9), qstep(9), q(9)

      integer nsteps
      parameter(nsteps = 100)

      integer i,j

			double precision en				!energy of groundstate

			double precision p(3000)	!parameter
c			integer pst(2,400)				!pointers for parameters

			integer npar							!number of parameters

			logical skip							!dummy

			double precision s(9),t(9) !transformed coordiantes

			double precision y(12)		!energies and cis

      include 'common.incl'

      skip=.false.
      include 'fort.66'

      p(2) = p(2) - p(1)
      p(1) = 0.

c     calculate number of parameters
      npar=pst(1,304)+pst(2,304)-1

      qstart(1) =  1.51250586081310     
      qstart(2) = -1.54890360362008     
      qstart(3) =  20.8905057397179       
      qstart(4) = 0.492453504006022          
      qstart(5) =  1.90381655585127     
      qstart(6) =  3.20534143377021      
      qstart(7) =  9.64612521292266       
      qstart(8) = -129.810995409334       
      qstart(9) =  101.360208931371     
 
      qend(1) =  1.26689688014270     
      qend(2) = -15.0022901054448            
      qend(3) = -2.34040972015036       
      qend(4) = 0.804978018430861          
      qend(5) =  2.13422356806596            
      qend(6) =  3.17893378465997        
      qend(7) =  50.1326524159092         
      qend(8) = -31.2210477514816         
      qend(9) = -166.357992153139     

      qstep(1) = (qend(1)-qstart(1))/(1.*nsteps) 
      qstep(2) = (qend(2)-qstart(2))/(1.*nsteps) 
      qstep(3) = (qend(3)-qstart(3))/(1.*nsteps) 
      qstep(4) = (qend(4)-qstart(4))/(1.*nsteps) 
      qstep(5) = (qend(5)-qstart(5))/(1.*nsteps) 
      qstep(6) = (qend(6)-qstart(6))/(1.*nsteps) 
      qstep(7) = (qend(7)-qstart(7))/(1.*nsteps) 
      qstep(8) = (qend(8)-qstart(8))/(1.*nsteps) 
      qstep(9) = (qend(9)-qstart(9))/(1.*nsteps) 
      

      do i=1,nsteps*10
         do j=1,9
            q(j) = qstep(j)*i + qstart(j) - qstep(j)*nsteps*5
         enddo

         write(6,*) 
         write(6,*) 'STEP:', i
         write(6,*) 'COORDINATES:'
         write(6,*) q
         
         call vwzprec(1,s)

c     get energy of groundstate
         call adia(y,1,q,p,npar,pst,skip)
         en=y(1)
         write(6,*) 'ENERGY:', en
      enddo
      
      end
