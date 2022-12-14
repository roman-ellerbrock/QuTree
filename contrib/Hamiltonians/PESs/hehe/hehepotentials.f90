function damp(n,eta,r)
implicit none
integer,parameter :: prec = kind(1.d0)
real(prec) :: damp
integer    :: n
real(prec) :: eta,r
real(prec) :: br,suma,term
integer    :: i
br = eta*r
if(br>1.d0) then
! 1 - exp(-x) sum(x^i/i!,i=0..n)
   term = 1.d0
   suma = 1.d0
   do i=1,n
      term = term*br/i
      suma = suma + term
   enddo
   damp = 1.d0 - exp(-br)*suma
else
! exp(-x) (exp(x) - sum(x^i/i!,i=0..n)) = exp(-x) sum(x^i/i!,i=n+1..inf)
   term = 1.d0
   do i=1,n
      term = term*br/i
   enddo
   suma = 0.d0
   do i=n+1,n+40
      term = term*br/i
      suma = suma + term
   enddo
   damp = suma*exp(-br)
endif
end function damp

function damp_mod(n,eta,r)
implicit none
integer,parameter :: prec = kind(1.d0)
real(prec) :: damp_mod
integer    :: n
real(prec) :: eta,r
real(prec) :: br,suma,term
integer    :: i
br = eta*r
! - exp(-x) sum(x^i/i!,i=0..n)
term = 1.d0
suma = 1.d0
do i=1,n
   term = term*br/i
   suma = suma + term
enddo
damp_mod = -exp(-br)*suma
end function damp_mod

!===============================================================================

module retardation
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: fsalpha = 1.d0/137.035999679d0
real(prec),parameter :: pi = 3.14159265358979324d0

real(prec),parameter :: polarizability = 1.38319217440d0
real(prec),parameter :: C6  = 1.460977837725d0
real(prec),parameter :: W4  = 0.35322d-04/fsalpha**2
real(prec),parameter :: AS3 = 0.577235d-06/fsalpha**3
real(prec),parameter :: K7  = 23.d0/(4.d0*pi)*polarizability**2/fsalpha
real(prec),parameter :: ratio = fsalpha*K7/C6

real(prec),parameter :: B1 =  8.454943177941253d0, &
                        B2 = 16.006586066260556d0, &
			B3 = 10.378373954734820d0, &
			B4 =  3.515803817223855d0, &
			B5 =  0.591502377533792d0, &
			B6 =  0.059455768329599d0

real(prec),parameter :: A1 = B1, &
                        A2 = B2 - W4/C6, &
			A3 = B3 - B1*W4/C6 + AS3/C6, &
			A4 = ratio*B5, &
			A5 = ratio*B6

private
public damp_ret

contains

real(prec) function damp_ret(r)
implicit none
real(prec) :: r,x,tmp1,tmp2
x = r*fsalpha
tmp1 = ((((A5*x + A4)*x + A3)*x + A2)*x + A1)*x + 1.d0
tmp2 =(((((B6*x + B5)*x + B4)*x + B3)*x + B2)*x + B1)*x + 1.d0
damp_ret = tmp1/tmp2
end function damp_ret

end module retardation

!===============================================================================

module Total_Fit
use retardation
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: a   = 3.64890303652830d0,&
                        b   = 2.36824871743591d0,&
                        eta = 4.09423805117871d0
real(prec),parameter :: P0  = -25.4701669416621d0,&
                        P1  = 269.244425630616d0,&
			P2  = -56.3879970402079d0,&
			Q0  =  38.7957487310071d0,&
			Q1  =  -2.76577136772754d0

real(prec),parameter :: C6BO = 1.460977837725d0
real(prec),parameter,dimension(3:16) :: Cn = &
(/     0.577235d-06,&
      -0.35322d-04,&
       0.1377841d-05,&
       1.46183d0,&
       0.d0,&
      14.1235d0,&
       0.d0,&
     183.7497d0,&
     -76.74d0,&
    3372.d0,&
   -3806.d0,&
   85340.d0,&
 -170700.d0,&
 2860000.d0&
/)

private
public total

contains

real(prec) function total(ret6,r)
implicit none
logical    :: ret6
real(prec) :: r,term,asy
integer    :: n
real(prec),external :: damp,damp_mod
term = (P0 + P1*r + P2*r**2)*exp(-a*r)
term = term + (Q0 + Q1*r)*exp(-b*r)
if(ret6) then
    asy = - damp_mod(3,eta,r)*Cn(3)/r**3 &
          - damp_mod(4,eta,r)*Cn(4)/r**4 &
	  - (damp_ret(r) + damp_mod(6,eta,r))*C6BO/r**6 &
	  - damp(6,eta,r)*(Cn(6) - C6BO)/r**6
else
    asy = - damp(3,eta,r)*Cn(3)/r**3 &
          - damp(4,eta,r)*Cn(4)/r**4 &
	  - damp(6,eta,r)*Cn(6)/r**6
endif
asy = asy - damp(5,eta,r)*Cn(5)/r**5 &
          - damp(8,eta,r)*Cn(8)/r**8
do n=10,16
    asy = asy - damp(n,eta,r)*Cn(n)/r**n
enddo
total = term + asy        ! in hartree
end function total

end module Total_Fit

!===============================================================================

module Total_Fit_sigma
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: C(3) = (/ 0.16702d-3, 0.4524d-5, 0.1843d-7 /)
real(prec),parameter :: A(3) = (/ 2.456d0,    1.100d0,   0.4381d0  /)

private
public sigma_total

contains

real(prec) function sigma_total(r)
implicit none
real(prec) :: r,term
integer    :: i
term = 0.d0
do i=1,3
   term = term + C(i)*exp(-A(i)*r)
enddo
sigma_total = term              ! in hartree
end function sigma_total

end module Total_Fit_sigma

!===============================================================================

module Born_Oppenheimer
use retardation
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: alpha = 3.65271949356113d0,&
                        beta  = 2.36720871471273d0,&
                        eta   = 4.09707982651218d0
real(prec),parameter :: P0  = -25.2616315711638d0,&
                        P1  = 272.551324570603d0,&
                        P2  = -56.5506587837287d0,&
			Q0  =  38.6339627675476d0,&
			Q1  =  -2.76537930016271d0

real(prec),parameter,dimension(6:16) :: Cn = &
(/     1.460977837725d0,&
       0.d0,&
      14.11785737d0,&
       0.d0,&
     183.691075d0,&
     -76.74d0,&
    3372.d0,&
   -3806.d0,&
   85340.d0,&
 -170700.d0,&
 2860000.d0/)

private
public BO

contains

real(prec) function BO(ret6,r)
implicit none
real(prec) :: r,term,asy
logical    :: ret6
integer    :: n
real(prec),external :: damp,damp_mod
term = (P0 + P1*r + P2*r**2)*exp(-alpha*r)
term = term + (Q0 + Q1*r)*exp(-beta*r)
if(ret6) then
   asy = - (damp_ret(r) + damp_mod(6,eta,r))*Cn(6)/r**6
else
   asy = - damp(6,eta,r)*Cn(6)/r**6
endif
asy = asy - damp(8,eta,r)*Cn(8)/r**8
do n=10,16
   asy = asy - damp(n,eta,r)*Cn(n)/r**n
enddo
BO = term + asy        ! in hartree
end function BO

end module Born_Oppenheimer

!===============================================================================

module adiabatic_correction
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: a   = 1.93303496128277d0,&
                        b   = 3.66494461616141d0,&
                        eta = 1.43603556683987d0
real(prec),parameter :: P0  = 0.377267595818441d-02,&
                        P1  =-0.668389341060405d-03,&
			Q0  = 0.241186947243731d-01,&
			Q1  =-0.219794150747149d-01

real(prec),parameter :: A6  = 0.0011445d0,&
                        A8  = 0.006519d0,&
			A10 = 0.0668d0

private
public ad

contains

real(prec) function ad(r)
implicit none
real(prec) :: r,term,asy
real(prec),external :: damp
term = (P0 + P1*r)*exp(-a*r)
term = term + (Q0 + Q1*r)*exp(-b*r)
asy = - damp(6,eta,r)*A6/r**6
asy = asy - damp(8,eta,r)*A8/r**8
asy = asy - damp(10,eta,r)*A10/r**10
ad = term + asy        ! in hartree
end function ad

end module adiabatic_correction

!===============================================================================

module Cowan_Griffin
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: a   = 2.90006031878770d0,&
                        b   = 2.24256968191719d0,&
                        c   = 0.425026783288716d-01,&
			eta = 2.95110060835781d0
real(prec),parameter :: P0  =-0.225884444742344d-02,&
                        P1  =-0.340134134669477d-02,&
                        Q0  = 0.369060490502093d-02,&
			Q1  =-0.163513409459839d-02

real(prec),parameter :: A6  =-0.000257d0,&
                        A8  =-0.00286d0,&
			A10 =-0.0401d0

private
public CG

contains

real(prec) function CG(r)
implicit none
real(prec) :: r,term,asy
real(prec),external :: damp
term = (P0 + P1*r)*exp(-a*r)
term = term + (Q0 + Q1*r)*exp(-b*r-c*r**2)
asy = - damp(6,eta,r)*A6/r**6
asy = asy - damp(8,eta,r)*A8/r**8
asy = asy - damp(10,eta,r)*A10/r**10
CG = term + asy        ! in hartree
end function CG

end module Cowan_Griffin

!===============================================================================

module Darwin_1el
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: fsalpha = 1.d0/137.035999679d0
real(prec),parameter :: pi = 3.14159265358979324d0

real(prec),parameter :: a   = 2.01042647863543d0,&
                        b   = 2.11313790518324d0,&
                        c   = 0.186163623408221d0,&
			eta = 4.07421274788735d0
real(prec),parameter :: P0  = 0.120025637135018d-01,&
                        P1  =-0.136621851503271d-02,&
                        Q0  = 0.601999646519627d-02,&
			Q1  =-0.445269531689116d-02

real(prec),parameter :: A6  = 0.001426d0,&
                        A8  = 0.01753d0,&
			A10 = 0.2735d0

real(prec),parameter :: ln_k0 = 4.3701602220d0
real(prec),parameter :: a2_a3 = fsalpha*(4.d0/3.d0)*&
                                (19.d0/30.d0 - 2.d0*log(fsalpha) - ln_k0)/&
				(pi/2.d0)

private
public a3D1

contains

real(prec) function D1(r)
implicit none
real(prec) :: r,term,asy
real(prec),external :: damp
term = (P0 + P1*r)*exp(-a*r)
term = term + (Q0 + Q1*r)*exp(-b*r-c*r**2)
asy = - damp(6,eta,r)*A6/r**6
asy = asy - damp(8,eta,r)*A8/r**8
asy = asy - damp(10,eta,r)*A10/r**10
D1 = term + asy        ! in hartree
end function D1

real(prec) function a3D1(r)
implicit none
real(prec) :: r
a3D1 = a2_a3*D1(r)
end function a3D1

end module Darwin_1el

!===============================================================================

module Darwin_2el
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: fsalpha = 1.d0/137.035999679d0
real(prec),parameter :: pi = 3.14159265358979324d0

real(prec),parameter :: a   = 1.85295073126638d0,&
                        b   = 2.52729337296464d0,&
                        c   = 0.351932113973596d0,&
			eta = 4.44314823119677d0
real(prec),parameter :: P0  = 0.544670094065963d-03,&
                        P1  =-0.699315570914574d-04,&
                        Q0  = 0.248652720153703d-02,&
			Q1  =-0.283371150273082d-04

real(prec),parameter :: A6  = 0.000103d0,&
                        A8  = 0.00136d0,&
			A10 = 0.0222d0

real(prec),parameter :: a2_a3 = fsalpha*&
                                (164.d0/15.d0 + 14.d0/3.d0*log(fsalpha))/&
                                pi

private
public D2,a3D2

contains

real(prec) function D2(r)
implicit none
real(prec) :: r,term,asy
real(prec),external :: damp
term = (P0 + P1*r)*exp(-a*r)
term = term + (Q0 + Q1*r)*exp(-b*r-c*r**2)
asy = - damp(6,eta,r)*A6/r**6
asy = asy - damp(8,eta,r)*A8/r**8
asy = asy - damp(10,eta,r)*A10/r**10
D2 = term + asy        ! in hartree
end function D2

real(prec) function a3D2(r)
implicit none
real(prec) :: r
a3D2 = a2_a3*D2(r)
end function a3D2

end module Darwin_2el

!===============================================================================

module Breit
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: a   = 1.74515046839322d0,&
                        b   = 5.00000000000000d0,&
                        c   = 0.263295834449133d0,&
                        eta = 3.45686578155164d0
real(prec),parameter :: P0  = 0.203310117594328d-03,&
                        P1  =-0.705614471514519d-04,&
                        P2  = 0.721004448825699d-05,&
                        Q0  =-0.270357188415393d-02,&
			Q1  = 0.516061855783347d-03

real(prec),parameter :: A4  =-0.000035322d0,&
                        A6  =-0.0001894d0

private
public Br

contains

real(prec) function Br(ret6,r)
implicit none
real(prec) :: r,term,asy
logical    :: ret6
real(prec),external :: damp,damp_mod
term = (P0 + P1*r + P2*r**2)*exp(-a*r)
term = term + (Q0 + Q1*r)*exp(-b*r - c*r**2)
if(ret6) then
    asy = - damp_mod(4,eta,r)*A4/r**4
else
    asy = - damp(4,eta,r)*A4/r**4
endif
asy = asy - damp(6,eta,r)*A6/r**6
Br = term + asy        ! in hartree
end function Br

end module Breit

!===============================================================================

module Araki_Sucher
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: a   = 0.262672224977970d0,&
                        b   = 2.02714404018130d0,&
                        eta = 2.88735558042463d0
real(prec),parameter :: P0  = 0.897489244934804d-09,&
                        P1  =-0.186471329166686d-09,&
                        Q0  = 0.186184585055017d-05,&
			Q1  =-0.152618822249290d-05

real(prec),parameter :: A3  = 0.577235d-06,&
                        A5  = 0.1377841d-05

private
public AS

contains

real(prec) function AS(ret6,r)
implicit none
real(prec) :: r,term,asy
logical    :: ret6
real(prec),external :: damp,damp_mod
term = (P0 + P1*r)*exp(-a*r)
term = term + (Q0 + Q1*r)*exp(-b*r)
if(ret6) then
    asy = - damp_mod(3,eta,r)*A3/r**3
else
    asy = - damp(3,eta,r)*A3/r**3
endif
asy = asy - damp(5,eta,r)*A5/r**5
AS = term + asy        ! in hartree
end function AS

end module Araki_Sucher

!===============================================================================

module potential_interface
use Total_Fit
use Total_Fit_sigma
use Born_Oppenheimer
use adiabatic_correction
use Cowan_Griffin
use Darwin_1el
use Darwin_2el
use Breit
use Araki_Sucher
implicit none
integer,parameter :: prec = kind(1.d0)

real(prec),parameter :: Mnuc4 = 7294.2995365d0
real(prec),parameter :: Mnuc3 = 5495.8852765d0
real(prec),parameter :: MU44 = Mnuc4/2.d0,&
                        MU34 = Mnuc3*Mnuc4/(Mnuc3+Mnuc4),&
                        MU33 = Mnuc3/2.d0

real(prec),parameter :: mult44 = 1.d0,&
                        mult34 = MU44/MU34,&
                        mult33 = MU44/MU33

private
public mult44,mult34,mult33
public V_tot,V_sigma,V_BO,V_ad,V_rel,V_QED

contains

real(prec) function V_tot(ret6,r)
implicit none
real(prec) :: r
logical    :: ret6
V_tot = total(ret6,r)
end function V_tot

real(prec) function V_sigma(r)
implicit none
real(prec) :: r
logical    :: ret6
V_sigma = sigma_total(r)
end function V_sigma

real(prec) function V_BO(ret6,r)
implicit none
real(prec) :: r
logical    :: ret6
V_BO = BO(ret6,r)
end function V_BO

real(prec) function V_ad(ret6,r)
implicit none
real(prec) :: r
logical    :: ret6
V_ad = ad(r)
end function V_ad

real(prec) function V_rel(ret6,r)
implicit none
real(prec) :: r
logical    :: ret6
V_rel = CG(r) + D2(r) + Br(ret6,r)
end function V_rel

real(prec) function V_QED(ret6,r)
implicit none
real(prec) :: r
logical    :: ret6
V_QED = a3D1(r) + a3D2(r) + AS(ret6,r)
end function V_QED

end module potential_interface
