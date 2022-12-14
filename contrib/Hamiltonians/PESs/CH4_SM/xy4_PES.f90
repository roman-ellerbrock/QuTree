program poten
  implicit none

  integer, parameter :: ark         = selected_real_kind(25,32)
  real(ark) :: deg, pi
  integer ipar, ieq, parmax, info, term(289)
  double precision  f, local(10), par_(2)
  double precision  force(289)
  double precision  alpha12, alpha13, alpha14, alpha23, alpha24, alpha34, r1, r2, r3, r4
  character(14) buf30

  interface
  function poten_xy4(local,parmax,par_,force,term)
    integer parmax, term(289)
    double precision local(10), par_(2),force(289), poten_xy4
  end function poten_xy4
  end interface

  pi = 4.0d0 * datan2(1.0d0,1.0d0)

  deg = pi/180.0d0

  !read parameters

  read*, parmax

  !read eq parameters
  do ieq =1, 2
    read*, buf30, par_(ieq)
  enddo
  !read potential parameters
  do ipar=1, parmax
    read*, buf30, term(ipar),force(ipar)
  enddo

  !read grid, print potential energy

  do
    read(*,*,end=4,err=4) local(1:10)

    r1  = local(1)
    r2  = local(2)
    r3  = local(3)
    r4  = local(4)

    alpha12 = local(5)*deg
    alpha13 = local(6)*deg
    alpha14 = local(7)*deg
    alpha23 = local(8)*deg
    alpha24 = local(9)*deg
    alpha34 = local(10)*deg

    f = poten_xy4(local, parmax, par_, force, term)
    write(*,'(10(1x,f12.7),2x,f14.6)') local(1:10), f
    cycle
4  exit
  enddo


end program poten


function poten_xy4(local,parmax,par_,force,term) result (f)
  implicit none

  integer, parameter :: ark         = selected_real_kind(25,32)
  integer, parameter :: ik          = selected_int_kind(8)

  integer parmax, ipar, term(289), k
  double precision local(10), force(289), par_(2), dF(289)
  double precision f

  integer(ik)          :: N

  real(ark) :: y1, y2, y3, y4, y5, y6, y7, y8, y9, r1e, alphae, a0, deg, pi
  real(ark) :: s1,s2,s3
  real(ark) :: r1,r2,r3,r4,alpha12,alpha13,alpha14,alpha23,alpha24,alpha34

  pi = 4.0d0 * datan2(1.0d0,1.0d0)

  deg = pi/180.0d0

  r1e    = par_(1)
  a0     = par_(2)

  r1  = local(1)
  r2  = local(2)
  r3  = local(3)
  r4  = local(4)

  alpha12 = local(5)*deg
  alpha13 = local(6)*deg
  alpha14 = local(7)*deg
  alpha23 = local(8)*deg
  alpha24 = local(9)*deg
  alpha34 = local(10)*deg

  y1=1.0_ark-exp(-a0*(r1-r1e))
  y2=1.0_ark-exp(-a0*(r2-r1e))
  y3=1.0_ark-exp(-a0*(r3-r1e))
  y4=1.0_ark-exp(-a0*(r4-r1e))
      
  y5=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
  y6=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
  y7=(alpha24-alpha13)/sqrt(2.0_ark)
  y8=(alpha23-alpha14)/sqrt(2.0_ark)
  y9=(alpha34-alpha12)/sqrt(2.0_ark)


      dF(1) = 0._ark
      dF(2) = 0._ark
      dF(3) = 1.0_ark
      dF(4) = y2+y3+y4+y1
      dF(5) = y8**2+y7**2+y9**2
      dF(6) = y6**2+y5**2
      dF(7) = (-y7-y8-y9)*y1+(y7-y9+y8)*y2+(y8+y9-y7)*y3+(y9+y7-y8)*y4
      dF(8) = (y4+y3+y2)*y1+(y4+y3)*y2+y3*y4
      dF(9) = y2**2+y3**2+y4**2+y1**2
      dF(10) = y7*y8*y9
      dF(11) = (-sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6
      dF(12) = y5**3-3._ark*y5*y6**2
      dF(13) = ((y8-2._ark*y9+y7)*y5+(sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y6)*y1+((-y8-2._ark*y9- &
          y7)*y5+(sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y6)*y2+((2._ark*y9+y7-y8)*y5+(-sqrt(3._ark)*y8- &
          sqrt(3._ark)*y7)*y6)*y3+((2._ark*y9+y8-y7)*y5+(sqrt(3._ark)*y7+sqrt(3._ark)*y8)*y6)*y4
      dF(14) = ((y9+y8)*y7+y8*y9)*y1+((y8-y9)*y7-y8*y9)*y2+((-y9-y8)*y7+y8*y9)*y3+((- &
          y8+y9)*y7-y8*y9)*y4
      dF(15) = (y8**2+y7**2+y9**2)*y1+(y8**2+y7**2+y9**2)*y2+(y8**2+y7**2+y9**2)*y3+ &
          (y8**2+y7**2+y9**2)*y4
      dF(16) = (y6**2+y5**2)*y1+(y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4
      dF(17) = (y3*y7+y4*y8+y2*y9)*y1+(-y3*y8-y4*y7)*y2-y3*y4*y9
      dF(18) = (y2*y5+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark-sqrt(3._ark)*y6/ &
          2._ark)*y4)*y1+((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark+sqrt(3._ark)*y6/ &
          2._ark)*y4)*y2+y3*y4*y5
      dF(19) = ((y4+y3)*y2+y3*y4)*y1+y2*y3*y4
      dF(20) = (y9+y7+y8)*y1**2+(-y7+y9-y8)*y2**2+(-y8-y9+y7)*y3**2+(-y7-y9+y8)*y4**2
      dF(21) = (y4+y3+y2)*y1**2+(y3**2+y2**2+y4**2)*y1+(y4+y3)*y2**2+(y3**2+y4**2)*y2+ &
          y3**2*y4+y3*y4**2
      dF(22) = y3**3+y1**3+y4**3+y2**3
      dF(23) = (y9**2+y8**2)*y7**2+y8**2*y9**2
      dF(24) = y9**4+y8**4+y7**4
      dF(25) = -sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/ &
          6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2
      dF(26) = (y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2
      dF(27) = y6**4+y5**4+2._ark*y5**2*y6**2
      dF(28) = (y3**2+y2**2+y4**2)*y1**2+(y3**2+y4**2)*y2**2+y3**2*y4**2
      dF(29) = (-y7-y8-y9)*y1**3+(y7-y9+y8)*y2**3+(y8+y9-y7)*y3**3+(y9+y7-y8)*y4**3
      dF(30) = y4**4+y3**4+y2**4+y1**4
      dF(31) = y3*y7*y8*y9+y2*y7*y8*y9+y4*y7*y8*y9+y1*y7*y8*y9
      dF(32) = ((y9+y8)*y7**2+(y9**2+y8**2)*y7+y8*y9**2+y8**2*y9)*y1+((-y8+y9)*y7**2+ &
          (-y9**2-y8**2)*y7+y8**2*y9-y8*y9**2)*y2+((-y9-y8)*y7**2+(y9**2+y8**2)*y7- &
          y8**2*y9-y8*y9**2)*y3+((y8-y9)*y7**2+(-y9**2-y8**2)*y7+y8*y9**2-y8**2*y9)*y4
      dF(33) = (y9**3+y8**3+y7**3)*y1+(y9**3-y8**3-y7**3)*y2+(-y8**3-y9**3+y7**3)*y3+ &
          (-y7**3-y9**3+y8**3)*y4
      dF(34) = ((-sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6)*y1+((-sqrt(3._ark)*y7**2/3._ark- &
          sqrt(3._ark)*y8**2/3._ark+2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6)*y2+((- &
          sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+ &
          (y7**2-y8**2)*y6)*y3+((-sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6)*y4
      dF(35) = (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y6)*y1+(((y8+y9/2._ark)*y7+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y2+(((y9/2._ark-y8)*y7-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3+(((-y9/2._ark- &
          y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4
      dF(36) = (-sqrt(3._ark)*y5**2*y9/2._ark+(y7-y8)*y6*y5+(sqrt(3._ark)*y9/6._ark- &
          sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y6**2)*y1+(-sqrt(3._ark)*y5**2*y9/2._ark+(y8- &
          y7)*y6*y5+(sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y6**2)*y2+ &
          (sqrt(3._ark)*y5**2*y9/2._ark+(y7+y8)*y6*y5+(sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark- &
          sqrt(3._ark)*y9/6._ark)*y6**2)*y3+(sqrt(3._ark)*y5**2*y9/2._ark+(-y8-y7)*y6*y5+(- &
          sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y6**2)*y4
      dF(37) = ((y9+y7+y8)*y5**2+(y9+y7+y8)*y6**2)*y1+((-y7+y9-y8)*y5**2+(-y7+y9- &
          y8)*y6**2)*y2+((-y8-y9+y7)*y5**2+(-y8-y9+y7)*y6**2)*y3+((-y7-y9+y8)*y5**2+(-y7- &
          y9+y8)*y6**2)*y4
      dF(38) = (y5**3-3._ark*y5*y6**2)*y1+(y5**3-3._ark*y5*y6**2)*y2+(y5**3- &
          3._ark*y5*y6**2)*y3+(y5**3-3._ark*y5*y6**2)*y4
      dF(39) = (y3*y7**2+y2*y9**2+y4*y8**2)*y1+(y4*y7**2+y3*y8**2)*y2+y3*y4*y9**2
      dF(40) = (y4*y7*y9+y3*y8*y9+y2*y7*y8)*y1+(-y4*y8*y9-y3*y7*y9)*y2-y3*y4*y7*y8
      dF(41) = ((y8**2+y7**2)*y2+(y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y1+((y7**2+ &
          y9**2)*y3+(y9**2+y8**2)*y4)*y2+(y8**2+y7**2)*y4*y3
      dF(42) = ((y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y1+((y6**2+ &
          y5**2)*y3+(y6**2+y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3
      dF(43) = ((sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark)*y2+(y5*y6+ &
          sqrt(3._ark)*y5**2/3._ark)*y3+(sqrt(3._ark)*y5**2/3._ark-y5*y6)*y4)*y1+ &
          ((sqrt(3._ark)*y5**2/3._ark-y5*y6)*y3+(y5*y6+sqrt(3._ark)*y5**2/3._ark)*y4)*y2+ &
          (sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark)*y4*y3
      dF(44) = (y2*y5*y9+(-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y3+(-y5*y8/2._ark- &
          sqrt(3._ark)*y6*y8/2._ark)*y4)*y1+((y5*y8/2._ark+sqrt(3._ark)*y6*y8/2._ark)*y3+(y5*y7/ &
          2._ark-sqrt(3._ark)*y6*y7/2._ark)*y4)*y2-y3*y4*y5*y9
      dF(45) = (((y9+y7-y8)*y3+(y8+y9-y7)*y4)*y2+(y7-y9+y8)*y4*y3)*y1+(-y7-y8- &
          y9)*y4*y3*y2
      dF(46) = y1*y2*y3*y4
      dF(47) = (y3*y7+y4*y8+y2*y9)*y1**2+(y4**2*y8+y3**2*y7+y2**2*y9)*y1+(-y3*y8- &
          y4*y7)*y2**2+(-y3**2*y8-y4**2*y7)*y2-y3*y4**2*y9-y3**2*y4*y9
      dF(48) = ((-y8-y7)*y2+(-y9-y8)*y3+(-y9-y7)*y4)*y1**2+((y7+y8)*y2**2+(y9+ &
          y8)*y3**2+(y7+y9)*y4**2)*y1+((y7-y9)*y3+(y8-y9)*y4)*y2**2+((y9-y7)*y3**2+(-y8+ &
          y9)*y4**2)*y2+(y8-y7)*y4*y3**2+(y7-y8)*y4**2*y3
      dF(49) = (y2*y5+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark-sqrt(3._ark)*y6/ &
          2._ark)*y4)*y1**2+(y2**2*y5+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4**2)*y1+((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark+ &
          sqrt(3._ark)*y6/2._ark)*y4)*y2**2+((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/2._ark+ &
          sqrt(3._ark)*y6/2._ark)*y4**2)*y2+y3*y4**2*y5+y3**2*y4*y5
      dF(50) = ((y4+y3)*y2+y3*y4)*y1**2+((y4+y3)*y2**2+(y3**2+y4**2)*y2+y3*y4**2+ &
          y3**2*y4)*y1+y2**2*y3*y4+(y3*y4**2+y3**2*y4)*y2
      dF(51) = (y4+y3+y2)*y1**3+(y4**3+y3**3+y2**3)*y1+(y4+y3)*y2**3+(y4**3+y3**3)*y2+ &
          y3*y4**3+y3**3*y4
      dF(52) = ((y9+y8)*y7+y8*y9)*y1**2+((y8-y9)*y7-y8*y9)*y2**2+((-y9-y8)*y7+ &
          y8*y9)*y3**2+((-y8+y9)*y7-y8*y9)*y4**2
      dF(53) = (y8**2+y7**2+y9**2)*y1**2+(y8**2+y7**2+y9**2)*y2**2+(y8**2+y7**2+ &
          y9**2)*y3**2+(y8**2+y7**2+y9**2)*y4**2
      dF(54) = (y6**2+y5**2)*y1**2+(y6**2+y5**2)*y2**2+(y6**2+y5**2)*y3**2+(y6**2+ &
          y5**2)*y4**2
      dF(55) = ((y9-y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y1**2+((y9+y7/2._ark+y8/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
          2._ark)*y6)*y2**2+((-y7/2._ark-y9+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
          2._ark)*y6)*y3**2+((y7/2._ark-y8/2._ark-y9)*y5+(-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y4**2
      dF(56) = y7**3*y8*y9+(y8**3*y9+y8*y9**3)*y7
      dF(57) = (2._ark/3._ark*sqrt(3._ark)*y9**4-sqrt(3._ark)*y8**4/3._ark-sqrt(3._ark)*y7**4/ &
          3._ark)*y5+(-y8**4+y7**4)*y6
      dF(58) = sqrt(3._ark)*y5**3*y9**2+(y7**2-y8**2)*y6*y5**2+(-sqrt(3._ark)*y9**2/3._ark- &
         4._ark/3._ark*sqrt(3._ark)*y7**2-4._ark/3._ark*sqrt(3._ark)*y8**2)*y6**2*y5+(y7**2- &
         y8**2)*y6**3
      dF(59) = ((-y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+ &
          (sqrt(3._ark)*y8**2*y9**2/2._ark-sqrt(3._ark)*y7**2*y9**2/2._ark)*y6
      dF(60) = y5**2*y7*y8*y9+y6**2*y7*y8*y9
      dF(61) = (y8**2+y7**2+y9**2)*y5**3+(-3._ark*y8**2-3._ark*y7**2-3._ark*y9**2)*y6**2*y5
      dF(62) = -3._ark*y5*y6**4+y5**5-2._ark*y5**3*y6**2
      s1 = (((y9/2._ark-y8)*y7**2+(y9**2/2._ark-y8**2)*y7+y8*y9**2/2._ark+y8**2*y9/2._ark)*y5+ &
          (-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y8**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y7**2*y9/2._ark)*y6)*y1+(((y8+y9/2._ark)*y7**2+(-y9**2/2._ark+y8**2)*y7- &
          y8*y9**2/2._ark+y8**2*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8*y9**2/ &
          2._ark+sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark)*y6)*y2
      dF(63) = s1+(((y8-y9/2._ark)*y7**2+(y9**2/2._ark-y8**2)*y7-y8**2*y9/2._ark-y8*y9**2/ &
          2._ark)*y5+(-sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark+ &
          sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8**2*y9/2._ark)*y6)*y3+(((-y9/2._ark-y8)*y7**2+ &
          (-y9**2/2._ark+y8**2)*y7-y8**2*y9/2._ark+y8*y9**2/2._ark)*y5+(sqrt(3._ark)*y8**2*y9/ &
          2._ark-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y6)*y4
      dF(64) = ((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(y8*y9- &
          y7*y9)*y6*y5+((sqrt(3._ark)*y9/6._ark+2._ark/3._ark*sqrt(3._ark)*y8)*y7+sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y1+((-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5**2+(-y8*y9+ &
          y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark+2._ark/3._ark*sqrt(3._ark)*y8)*y7-sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y2+((-sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/2._ark)*y5**2+(y8*y9+ &
          y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark-2._ark/3._ark*sqrt(3._ark)*y8)*y7+sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y3+((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(-y8*y9- &
          y7*y9)*y6*y5+((-2._ark/3._ark*sqrt(3._ark)*y8+sqrt(3._ark)*y9/6._ark)*y7-sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y4
      dF(65) = (-sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/ &
          6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y1+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y2+(-sqrt(3._ark)*y5**2*y9**2/ &
          2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/3._ark- &
          sqrt(3._ark)*y7**2/3._ark)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2- &
          y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/ &
          3._ark)*y6**2)*y4
      dF(66) = (((y9+y8)*y7+y8*y9)*y5**2+((y9+y8)*y7+y8*y9)*y6**2)*y1+(((y8-y9)*y7- &
          y8*y9)*y5**2+((y8-y9)*y7-y8*y9)*y6**2)*y2+(((-y9-y8)*y7+y8*y9)*y5**2+((-y9- &
          y8)*y7+y8*y9)*y6**2)*y3+(((-y8+y9)*y7-y8*y9)*y5**2+((-y8+y9)*y7-y8*y9)*y6**2)*y4
      dF(67) = ((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y1+((y8**2+y7**2+ &
          y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y2+((y8**2+y7**2+y9**2)*y5**2+(y8**2+ &
          y7**2+y9**2)*y6**2)*y3+((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y4
      dF(68) = ((sqrt(3._ark)*y7+sqrt(3._ark)*y8)*y5**3+(y8-y7)*y6*y5**2+(-5._ark/ &
          3._ark*sqrt(3._ark)*y7-5._ark/3._ark*sqrt(3._ark)*y8-8._ark/3._ark*sqrt(3._ark)*y9)*y6**2*y5+ &
          (y8-y7)*y6**3)*y1+((-sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y5**3+(y7-y8)*y6*y5**2+(-8._ark/ &
          3._ark*sqrt(3._ark)*y9+5._ark/3._ark*sqrt(3._ark)*y7+5._ark/3._ark*sqrt(3._ark)*y8)*y6**2*y5+ &
          (y7-y8)*y6**3)*y2+((sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y5**3+(-y8-y7)*y6*y5**2+(8._ark/ &
          3._ark*sqrt(3._ark)*y9+5._ark/3._ark*sqrt(3._ark)*y8-5._ark/3._ark*sqrt(3._ark)*y7)*y6**2*y5+(- &
          y8-y7)*y6**3)*y3+((sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y5**3+(y7+y8)*y6*y5**2+(5._ark/ &
          3._ark*sqrt(3._ark)*y7-5._ark/3._ark*sqrt(3._ark)*y8+8._ark/3._ark*sqrt(3._ark)*y9)*y6**2*y5+ &
          (y7+y8)*y6**3)*y4
      dF(69) = ((y9+y7+y8)*y5**3+(-3._ark*y8-3._ark*y9-3._ark*y7)*y6**2*y5)*y1+((-y7+y9- &
          y8)*y5**3+(3._ark*y7+3._ark*y8-3._ark*y9)*y6**2*y5)*y2+((-y8-y9+y7)*y5**3+(-3._ark*y7+ &
          3._ark*y8+3._ark*y9)*y6**2*y5)*y3+((-y7-y9+y8)*y5**3+(3._ark*y7+3._ark*y9- &
          3._ark*y8)*y6**2*y5)*y4
      dF(70) = (y6**4+y5**4+2._ark*y5**2*y6**2)*y1+(y6**4+y5**4+2._ark*y5**2*y6**2)*y2+ &
          (y6**4+y5**4+2._ark*y5**2*y6**2)*y3+(y6**4+y5**4+2._ark*y5**2*y6**2)*y4
      dF(71) = (y4*y7*y8*y9+y3*y7*y8*y9+y2*y7*y8*y9)*y1+(y4*y7*y8*y9+y3*y7*y8*y9)*y2+ &
          y3*y4*y7*y8*y9
      dF(72) = ((-y8**2*y9-y7**2*y9)*y2+(-y9**2-y8**2)*y7*y3+(-y8*y9**2- &
          y7**2*y8)*y4)*y1+((y7**2*y8+y8*y9**2)*y3+(y9**2+y8**2)*y7*y4)*y2+(y7**2*y9+ &
          y8**2*y9)*y4*y3
      dF(73) = (-y2*y9**3-y4*y8**3-y3*y7**3)*y1+(y4*y7**3+y3*y8**3)*y2+y3*y4*y9**3
      dF(74) = (((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5+(y8**2/2._ark-y7**2/ &
          2._ark)*y6)*y2+(-sqrt(3._ark)*y5*y9**2/2._ark+(y9**2/2._ark+y8**2)*y6)*y3+(- &
          sqrt(3._ark)*y5*y9**2/2._ark+(-y7**2-y9**2/2._ark)*y6)*y4)*y1+((-sqrt(3._ark)*y5*y9**2/ &
          2._ark+(-y7**2-y9**2/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y9**2/2._ark+(y9**2/2._ark+ &
          y8**2)*y6)*y4)*y2+((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5+(y8**2/2._ark- &
          y7**2/2._ark)*y6)*y4*y3
      dF(75) = (2._ark*y2*y5*y7*y8+(-y5*y8*y9+sqrt(3._ark)*y6*y8*y9)*y3+(- &
          sqrt(3._ark)*y6*y7*y9-y5*y7*y9)*y4)*y1+((y5*y7*y9+sqrt(3._ark)*y6*y7*y9)*y3+(- &
          sqrt(3._ark)*y6*y8*y9+y5*y8*y9)*y4)*y2-2._ark*y3*y4*y5*y7*y8
      dF(76) = (((-y8**2/2._ark-y7**2/2._ark)*y5+(sqrt(3._ark)*y8**2/2._ark-sqrt(3._ark)*y7**2/ &
          2._ark)*y6)*y2+((-y9**2/2._ark+y8**2)*y5-sqrt(3._ark)*y6*y9**2/2._ark)*y3+((y7**2-y9**2/ &
          2._ark)*y5+sqrt(3._ark)*y6*y9**2/2._ark)*y4)*y1+(((y7**2-y9**2/2._ark)*y5+ &
          sqrt(3._ark)*y6*y9**2/2._ark)*y3+((-y9**2/2._ark+y8**2)*y5-sqrt(3._ark)*y6*y9**2/ &
          2._ark)*y4)*y2+((-y8**2/2._ark-y7**2/2._ark)*y5+(sqrt(3._ark)*y8**2/2._ark- &
          sqrt(3._ark)*y7**2/2._ark)*y6)*y4*y3
      dF(77) = (-2._ark*y2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y3+(y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y4)*y1+((y5*y8**2+sqrt(3._ark)*y6*y8**2)*y3+(- &
          sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4)*y2-2._ark*y3*y4*y5*y9**2
      dF(78) = ((-sqrt(3._ark)*y6**2*y9/6._ark+sqrt(3._ark)*y5**2*y9/2._ark)*y2+(-y5*y6*y7+ &
          sqrt(3._ark)*y6**2*y7/3._ark)*y3+(sqrt(3._ark)*y6**2*y8/3._ark+y5*y6*y8)*y4)*y1+((- &
          sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y3+(-sqrt(3._ark)*y6**2*y7/3._ark+ &
          y5*y6*y7)*y4)*y2+(sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y4*y3
      dF(79) = ((-y5**3/3._ark+y5*y6**2)*y2+(-y5**3/3._ark+y5*y6**2)*y3+(-y5**3/3._ark+ &
          y5*y6**2)*y4)*y1+((-y5**3/3._ark+y5*y6**2)*y3+(-y5**3/3._ark+y5*y6**2)*y4)*y2+(- &
          y5**3/3._ark+y5*y6**2)*y4*y3
      dF(80) = ((-y9*y6**2-y9*y5**2)*y2+(-y5**2*y7-y7*y6**2)*y3+(-y5**2*y8- &
          y8*y6**2)*y4)*y1+((y8*y6**2+y5**2*y8)*y3+(y7*y6**2+y5**2*y7)*y4)*y2+(y9*y6**2+ &
          y9*y5**2)*y4*y3
      dF(81) = ((5._ark/9._ark*sqrt(3._ark)*y5**3+sqrt(3._ark)*y5*y6**2)*y2+(y6**3+y5**2*y6- &
          4._ark/9._ark*sqrt(3._ark)*y5**3)*y3+(-y6**3-y5**2*y6-4._ark/ &
          9._ark*sqrt(3._ark)*y5**3)*y4)*y1+((-y6**3-y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5**3)*y3+ &
          (y6**3+y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5**3)*y4)*y2+(5._ark/9._ark*sqrt(3._ark)*y5**3+ &
          sqrt(3._ark)*y5*y6**2)*y4*y3
      dF(82) = (-y3*y8*y9-y2*y7*y8-y4*y7*y9)*y1**2+(-y4**2*y7*y9-y2**2*y7*y8- &
          y3**2*y8*y9)*y1+(y3*y7*y9+y4*y8*y9)*y2**2+(y4**2*y8*y9+y3**2*y7*y9)*y2+ &
          y3**2*y4*y7*y8+y3*y4**2*y7*y8
      dF(83) = ((y8**2+y7**2)*y2+(y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y1**2+((y8**2+ &
          y7**2)*y2**2+(y9**2+y8**2)*y3**2+(y7**2+y9**2)*y4**2)*y1+((y7**2+y9**2)*y3+ &
          (y9**2+y8**2)*y4)*y2**2+((y7**2+y9**2)*y3**2+(y9**2+y8**2)*y4**2)*y2+(y8**2+ &
          y7**2)*y4*y3**2+(y8**2+y7**2)*y4**2*y3
      dF(84) = ((-y8*y9-y7*y9)*y2+(-y9-y8)*y7*y3+(-y7*y8-y8*y9)*y4)*y1**2+((y8*y9+ &
          y7*y9)*y2**2+(y9+y8)*y7*y3**2+(y8*y9+y7*y8)*y4**2)*y1+((-y7*y8+y8*y9)*y3+(-y8+ &
          y9)*y7*y4)*y2**2+((-y8*y9+y7*y8)*y3**2+(y8-y9)*y7*y4**2)*y2+(-y8*y9+ &
          y7*y9)*y4*y3**2+(y8*y9-y7*y9)*y4**2*y3
      dF(85) = (y3*y7**2+y2*y9**2+y4*y8**2)*y1**2+(y4**2*y8**2+y3**2*y7**2+ &
          y2**2*y9**2)*y1+(y4*y7**2+y3*y8**2)*y2**2+(y4**2*y7**2+y3**2*y8**2)*y2+ &
          y3*y4**2*y9**2+y3**2*y4*y9**2
      s1 = (((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/2._ark)*y6)*y2+(- &
          sqrt(3._ark)*y5*y9/2._ark+(y8+y9/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark- &
          y7)*y6)*y4)*y1**2+(((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark-y8/ &
          2._ark)*y6)*y2**2+(sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark-y8)*y6)*y3**2+ &
          (sqrt(3._ark)*y5*y9/2._ark+(y7+y9/2._ark)*y6)*y4**2)*y1+((-sqrt(3._ark)*y5*y9/2._ark+(-y9/ &
         2._ark+y7)*y6)*y3+(-sqrt(3._ark)*y5*y9/2._ark+(y9/2._ark-y8)*y6)*y4)*y2**2
      dF(86) = s1+((sqrt(3._ark)*y5*y9/2._ark+(-y7+y9/2._ark)*y6)*y3**2+(sqrt(3._ark)*y5*y9/ &
          2._ark+(y8-y9/2._ark)*y6)*y4**2)*y2+((sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(- &
          y8/2._ark-y7/2._ark)*y6)*y4*y3**2+((-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y7/ &
          2._ark+y8/2._ark)*y6)*y4**2*y3
      dF(87) = ((y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y1**2+((y6**2+ &
          y5**2)*y2**2+(y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y1+((y6**2+y5**2)*y3+ &
          (y6**2+y5**2)*y4)*y2**2+((y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y2+(y6**2+ &
          y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3
      s1 = (((-y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6)*y2+ &
          ((y8-y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y3+((-y9/2._ark+y7)*y5+sqrt(3._ark)*y6*y9/ &
          2._ark)*y4)*y1**2+(((y7/2._ark+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y2**2+((y9/2._ark-y8)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3**2+((-y7+y9/2._ark)*y5- &
          sqrt(3._ark)*y6*y9/2._ark)*y4**2)*y1+(((-y9/2._ark-y7)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3+ &
          ((-y9/2._ark-y8)*y5-sqrt(3._ark)*y6*y9/2._ark)*y4)*y2**2
      dF(88) = s1+(((y7+y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y3**2+((y8+y9/2._ark)*y5+ &
          sqrt(3._ark)*y6*y9/2._ark)*y4**2)*y2+((y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y6)*y4*y3**2+((y7/2._ark-y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6)*y4**2*y3
      dF(89) = ((((-y8+y9)*y7-y8*y9)*y3+((-y9-y8)*y7+y8*y9)*y4)*y2+((y8-y9)*y7- &
          y8*y9)*y4*y3)*y1+((y9+y8)*y7+y8*y9)*y4*y3*y2
      dF(90) = (((y8**2+y7**2+y9**2)*y3+(y8**2+y7**2+y9**2)*y4)*y2+(y8**2+y7**2+ &
          y9**2)*y4*y3)*y1+(y8**2+y7**2+y9**2)*y4*y3*y2
      dF(91) = (((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3)*y1+(y6**2+ &
          y5**2)*y4*y3*y2
      dF(92) = ((-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y2+(-y5*y6- &
          sqrt(3._ark)*y5**2/3._ark)*y3+(y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4)*y1**2+((- &
          sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y2**2+(-y5*y6-sqrt(3._ark)*y5**2/ &
          3._ark)*y3**2+(y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4**2)*y1+((y5*y6-sqrt(3._ark)*y5**2/ &
          3._ark)*y3+(-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4)*y2**2+((y5*y6-sqrt(3._ark)*y5**2/ &
          3._ark)*y3**2+(-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4**2)*y2+(-sqrt(3._ark)*y6**2/2._ark+ &
          sqrt(3._ark)*y5**2/6._ark)*y4*y3**2+(-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/ &
          6._ark)*y4**2*y3
      dF(93) = ((y3*y8+y4*y7)*y2+y3*y4*y9)*y1**2+((-y3*y7-y4*y8)*y2**2+(-y4**2*y9- &
          y3**2*y9)*y2-y3*y4**2*y7-y3**2*y4*y8)*y1+y2**2*y3*y4*y9+(y3*y4**2*y8+ &
          y3**2*y4*y7)*y2
      dF(94) = (((-sqrt(3._ark)*y5/3._ark-y6)*y3+(-sqrt(3._ark)*y5/3._ark+y6)*y4)*y2+2._ark/ &
          3._ark*sqrt(3._ark)*y3*y4*y5)*y1**2+(((-sqrt(3._ark)*y5/3._ark+y6)*y3+(-sqrt(3._ark)*y5/ &
          3._ark-y6)*y4)*y2**2+(2._ark/3._ark*sqrt(3._ark)*y4**2*y5+2._ark/ &
          3._ark*sqrt(3._ark)*y3**2*y5)*y2+(-sqrt(3._ark)*y5/3._ark-y6)*y4*y3**2+(-sqrt(3._ark)*y5/ &
          3._ark+y6)*y4**2*y3)*y1+2._ark/3._ark*sqrt(3._ark)*y2**2*y3*y4*y5+((-sqrt(3._ark)*y5/3._ark+ &
          y6)*y4*y3**2+(-sqrt(3._ark)*y5/3._ark-y6)*y4**2*y3)*y2
      dF(95) = ((y4+y3)*y2**2+(y3**2+y4**2)*y2+y3*y4**2+y3**2*y4)*y1**2+((y3**2+ &
          y4**2)*y2**2+y3**2*y4**2)*y1+(y3*y4**2+y3**2*y4)*y2**2+y2*y3**2*y4**2
      dF(96) = (-y4*y8-y3*y7-y2*y9)*y1**3+(-y4**3*y8-y3**3*y7-y2**3*y9)*y1+(y3*y8+ &
          y4*y7)*y2**3+(y3**3*y8+y4**3*y7)*y2+y3**3*y4*y9+y3*y4**3*y9
      dF(97) = ((y7+y8)*y2+(y9+y8)*y3+(y7+y9)*y4)*y1**3+((-y8-y7)*y2**3+(-y9- &
          y8)*y3**3+(-y9-y7)*y4**3)*y1+((y9-y7)*y3+(-y8+y9)*y4)*y2**3+((y7-y9)*y3**3+(y8- &
          y9)*y4**3)*y2+(y7-y8)*y4*y3**3+(y8-y7)*y4**3*y3
      dF(98) = (-2._ark/3._ark*sqrt(3._ark)*y2*y5+(-y6+sqrt(3._ark)*y5/3._ark)*y3+(y6+ &
          sqrt(3._ark)*y5/3._ark)*y4)*y1**3+(-2._ark/3._ark*sqrt(3._ark)*y2**3*y5+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y3**3+(y6+sqrt(3._ark)*y5/3._ark)*y4**3)*y1+((y6+sqrt(3._ark)*y5/ &
          3._ark)*y3+(-y6+sqrt(3._ark)*y5/3._ark)*y4)*y2**3+((y6+sqrt(3._ark)*y5/3._ark)*y3**3+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**3)*y2-2._ark/3._ark*sqrt(3._ark)*y3**3*y4*y5-2._ark/ &
          3._ark*sqrt(3._ark)*y3*y4**3*y5
      dF(99) = ((y4+y3)*y2+y3*y4)*y1**3+((y4+y3)*y2**3+(y4**3+y3**3)*y2+y3**3*y4+ &
          y3*y4**3)*y1+y2**3*y3*y4+(y3**3*y4+y3*y4**3)*y2
      dF(100) = (y4+y3+y2)*y1**4+(y4**4+y2**4+y3**4)*y1+(y4+y3)*y2**4+(y4**4+ &
          y3**4)*y2+y3*y4**4+y3**4*y4
      dF(101) = y2**2*y7*y8*y9+y3**2*y7*y8*y9+y1**2*y7*y8*y9+y4**2*y7*y8*y9
      dF(102) = ((-y9-y8)*y7**2+(-y9**2-y8**2)*y7-y8*y9**2-y8**2*y9)*y1**2+((y8- &
          y9)*y7**2+(y9**2+y8**2)*y7-y8**2*y9+y8*y9**2)*y2**2+((y9+y8)*y7**2+(-y9**2- &
          y8**2)*y7+y8**2*y9+y8*y9**2)*y3**2+((-y8+y9)*y7**2+(y9**2+y8**2)*y7-y8*y9**2+ &
          y8**2*y9)*y4**2
      dF(103) = (-y8**3-y7**3-y9**3)*y1**2+(-y9**3+y8**3+y7**3)*y2**2+(y8**3+y9**3- &
          y7**3)*y3**2+(y7**3+y9**3-y8**3)*y4**2
      dF(104) = (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y6)*y1**2+(((y8+y9/2._ark)*y7+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y2**2+(((y9/2._ark-y8)*y7-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3**2+(((-y9/2._ark- &
          y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4**2
      dF(105) = ((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y1**2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y2**2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y3**2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4**2
      dF(106) = ((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5**2+(y8-y7)*y6*y5+(- &
          sqrt(3._ark)*y8/6._ark-2._ark/3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark)*y6**2)*y1**2+ &
          ((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5**2+(y7-y8)*y6*y5+(-2._ark/ &
          3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y7/6._ark)*y6**2)*y2**2+((- &
          sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5**2+(-y8-y7)*y6*y5+(-sqrt(3._ark)*y7/ &
          6._ark+2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/6._ark)*y6**2)*y3**2+((sqrt(3._ark)*y7/ &
          2._ark-sqrt(3._ark)*y8/2._ark)*y5**2+(y7+y8)*y6*y5+(2._ark/3._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/6._ark)*y6**2)*y4**2
      dF(107) = ((y9+y7+y8)*y5**2+(y9+y7+y8)*y6**2)*y1**2+((-y7+y9-y8)*y5**2+(-y7+y9- &
          y8)*y6**2)*y2**2+((-y8-y9+y7)*y5**2+(-y8-y9+y7)*y6**2)*y3**2+((-y7-y9+y8)*y5**2+ &
          (-y7-y9+y8)*y6**2)*y4**2
      dF(108) = (y5**3-3._ark*y5*y6**2)*y1**2+(y5**3-3._ark*y5*y6**2)*y2**2+(y5**3- &
          3._ark*y5*y6**2)*y3**2+(y5**3-3._ark*y5*y6**2)*y4**2
      dF(109) = (2._ark/3._ark*sqrt(3._ark)*y2*y5*y9+(-sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3+(- &
          sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4)*y1**2+(2._ark/3._ark*sqrt(3._ark)*y2**2*y5*y9+(- &
          sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3**2+(-sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4**2)*y1+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4)*y2**2+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3**2+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4**2)*y2- &
          2._ark/3._ark*sqrt(3._ark)*y3*y4**2*y5*y9-2._ark/3._ark*sqrt(3._ark)*y3**2*y4*y5*y9
      dF(110) = (-y3**2*y7-y4**2*y8-y2**2*y9)*y1**2+(y4**2*y7+y3**2*y8)*y2**2+ &
          y3**2*y4**2*y9
      dF(111) = (-2._ark/3._ark*sqrt(3._ark)*y2**2*y5+(-y6+sqrt(3._ark)*y5/3._ark)*y3**2+(y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**2)*y1**2+((y6+sqrt(3._ark)*y5/3._ark)*y3**2+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**2)*y2**2-2._ark/3._ark*sqrt(3._ark)*y3**2*y4**2*y5
      dF(112) = ((y9+y8)*y7+y8*y9)*y1**3+((y8-y9)*y7-y8*y9)*y2**3+((-y9-y8)*y7+ &
          y8*y9)*y3**3+((-y8+y9)*y7-y8*y9)*y4**3
      dF(113) = (y8**2+y7**2+y9**2)*y1**3+(y8**2+y7**2+y9**2)*y2**3+(y8**2+y7**2+ &
          y9**2)*y3**3+(y8**2+y7**2+y9**2)*y4**3
      dF(114) = ((-2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5+ &
          (y8-y7)*y6)*y1**3+((-sqrt(3._ark)*y7/3._ark-2._ark/3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y8/ &
          3._ark)*y5+(y7-y8)*y6)*y2**3+((sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5+(-y8-y7)*y6)*y3**3+((2._ark/3._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark)*y5+(y7+y8)*y6)*y4**3
      dF(115) = (y6**2+y5**2)*y1**3+(y6**2+y5**2)*y2**3+(y6**2+y5**2)*y3**3+(y6**2+ &
          y5**2)*y4**3
      dF(116) = (y3**2+y2**2+y4**2)*y1**3+(y4**3+y3**3+y2**3)*y1**2+(y3**2+ &
          y4**2)*y2**3+(y4**3+y3**3)*y2**2+y3**2*y4**3+y3**3*y4**2
      dF(117) = (-y7-y8-y9)*y1**4+(y7-y9+y8)*y2**4+(y8+y9-y7)*y3**4+(y9+y7-y8)*y4**4
      dF(118) = y4**5+y3**5+y2**5+y1**5
      dF(119) = (y7**2*y8*y9+(y8*y9**2+y8**2*y9)*y7)*y1+(-y7**2*y8*y9+(y8*y9**2- &
          y8**2*y9)*y7)*y2+(y7**2*y8*y9+(-y8**2*y9-y8*y9**2)*y7)*y3+(-y7**2*y8*y9+(- &
          y8*y9**2+y8**2*y9)*y7)*y4
      dF(120) = ((y9**2+y8**2)*y7**2+y8**2*y9**2)*y1+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y2+((y9**2+y8**2)*y7**2+y8**2*y9**2)*y3+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y4
      dF(121) = ((y9+y8)*y7**3+(y9**3+y8**3)*y7+y8**3*y9+y8*y9**3)*y1+((y8-y9)*y7**3+ &
          (y8**3-y9**3)*y7-y8*y9**3-y8**3*y9)*y2+((-y9-y8)*y7**3+(-y9**3-y8**3)*y7+ &
          y8**3*y9+y8*y9**3)*y3+((-y8+y9)*y7**3+(y9**3-y8**3)*y7-y8*y9**3-y8**3*y9)*y4
      dF(122) = (y9**4+y8**4+y7**4)*y1+(y9**4+y8**4+y7**4)*y2+(y9**4+y8**4+y7**4)*y3+ &
          (y9**4+y8**4+y7**4)*y4
      s1 = ((sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**2/2._ark)*y5+((y8+y9/2._ark)*y7**2+(-y9**2/2._ark-y8**2)*y7- &
          y8**2*y9/2._ark+y8*y9**2/2._ark)*y6)*y1+((-sqrt(3._ark)*y7*y9**2/2._ark- &
          sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark)*y5+ &
          ((y9/2._ark-y8)*y7**2+(y9**2/2._ark+y8**2)*y7-y8*y9**2/2._ark-y8**2*y9/2._ark)*y6)*y2
      dF(123) = s1+((sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y8*y9**2/2._ark)*y5+((-y9/2._ark-y8)*y7**2+(- &
          y9**2/2._ark-y8**2)*y7+y8**2*y9/2._ark-y8*y9**2/2._ark)*y6)*y3+((sqrt(3._ark)*y8**2*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y5+((y8-y9/2._ark)*y7**2+(y9**2/2._ark+y8**2)*y7+y8**2*y9/2._ark+y8*y9**2/ &
          2._ark)*y6)*y4
      dF(124) = ((-sqrt(3._ark)*y7**3/3._ark-sqrt(3._ark)*y8**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3)*y5+(-y8**3+y7**3)*y6)*y1+((sqrt(3._ark)*y7**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3+y8**3)*y6)*y2+((-2._ark/ &
          3._ark*sqrt(3._ark)*y9**3-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/3._ark)*y5+(y8**3+ &
          y7**3)*y6)*y3+((-2._ark/3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y7**3/3._ark- &
          sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3-y8**3)*y6)*y4
      dF(125) = ((((-sqrt(3._ark)*y8/3._ark-2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/ &
          3._ark)*y5+(-y8-y7)*y6)*y3+((-sqrt(3._ark)*y7/3._ark+sqrt(3._ark)*y8/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5+(y7+y8)*y6)*y4)*y2+((sqrt(3._ark)*y8/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/3._ark)*y5+(y8-y7)*y6)*y4*y3)*y1+((2._ark/ &
          3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y5+(y7- &
          y8)*y6)*y4*y3*y2
      dF(126) = (((y7+y9)*y3+(y9+y8)*y4)*y2+(y7+y8)*y4*y3)*y1**2+(((-y8+y9)*y3+(y9- &
          y7)*y4)*y2**2+((y7-y8)*y3**2+(y8-y7)*y4**2)*y2+(y7-y9)*y4*y3**2+(y8- &
          y9)*y4**2*y3)*y1+(-y8-y7)*y4*y3*y2**2+((-y9-y8)*y4*y3**2+(-y9-y7)*y4**2*y3)*y2
      dF(127) = y1**2*y2*y3*y4+(y2**2*y3*y4+(y3*y4**2+y3**2*y4)*y2)*y1
      dF(128) = (y9**2+y8**2)*y7**4+(y8**4+y9**4)*y7**2+y8**4*y9**2+y8**2*y9**4
      dF(129) = y7**6+y9**6+y8**6
      dF(130) = y7**2*y8**2*y9**2
      dF(131) = ((y9**2+y8**2)*y7**2+y8**2*y9**2)*y5**2+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y6**2
      dF(132) = -6._ark*y5**2*y6**4+9._ark*y5**4*y6**2+y6**6
      dF(133) = (y7**3*y8*y9+(-2._ark*y8*y9**3+y8**3*y9)*y7)*y5+(sqrt(3._ark)*y7*y8**3*y9- &
          sqrt(3._ark)*y7**3*y8*y9)*y6
      dF(134) = ((2._ark/3._ark*sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2/6._ark)*y7**2+ &
          sqrt(3._ark)*y8**2*y9**2/6._ark)*y5**2+(-y8**2*y9**2+y7**2*y9**2)*y6*y5+ &
          (sqrt(3._ark)*y7**2*y9**2/2._ark+sqrt(3._ark)*y8**2*y9**2/2._ark)*y6**2
      dF(135) = -sqrt(3._ark)*y5**2*y9**4/2._ark+(-y8**4+y7**4)*y6*y5+(-sqrt(3._ark)*y7**4/ &
          3._ark-sqrt(3._ark)*y8**4/3._ark+sqrt(3._ark)*y9**4/6._ark)*y6**2
      dF(136) = -y5**3*y7*y8*y9/3._ark+y5*y6**2*y7*y8*y9
      dF(137) = (y9**4+y8**4+y7**4)*y5**2+(y9**4+y8**4+y7**4)*y6**2
      dF(138) = 3._ark/4._ark*y5**4*y9**2+(-y9**2/2._ark+y7**2+y8**2)*y6**2*y5**2+(-2._ark/ &
          3._ark*sqrt(3._ark)*y7**2+2._ark/3._ark*sqrt(3._ark)*y8**2)*y6**3*y5+(y8**2/3._ark+y9**2/ &
          12._ark+y7**2/3._ark)*y6**4
      dF(139) = -sqrt(3._ark)*y5**4*y9**2/4._ark+(y7**2-y8**2)*y6*y5**3- &
          sqrt(3._ark)*y5**2*y6**2*y9**2/2._ark+(-y8**2/3._ark+y7**2/3._ark)*y6**3*y5+(-2._ark/ &
          9._ark*sqrt(3._ark)*y8**2+7._ark/36._ark*sqrt(3._ark)*y9**2-2._ark/ &
          9._ark*sqrt(3._ark)*y7**2)*y6**4
      dF(140) = (-y9**2/2._ark+y7**2+y8**2)*y5**4+3._ark*y5**2*y6**2*y9**2+(4._ark/ &
          3._ark*sqrt(3._ark)*y7**2-4._ark/3._ark*sqrt(3._ark)*y8**2)*y6**3*y5+(y8**2/3._ark+5._ark/ &
          6._ark*y9**2+y7**2/3._ark)*y6**4
      dF(141) = 9._ark*y5**2*y6**4+y5**6-6._ark*y5**4*y6**2
      dF(142) = ((y9**2+y8**2)*y7**3+(y9**3+y8**3)*y7**2+y8**3*y9**2+y8**2*y9**3)*y1+ &
          ((-y9**2-y8**2)*y7**3+(y9**3-y8**3)*y7**2+y8**2*y9**3-y8**3*y9**2)*y2+((y9**2+ &
          y8**2)*y7**3+(-y9**3-y8**3)*y7**2-y8**2*y9**3-y8**3*y9**2)*y3+((-y9**2- &
          y8**2)*y7**3+(y8**3-y9**3)*y7**2+y8**3*y9**2-y8**2*y9**3)*y4
      dF(143) = ((y8*y9**2+y8**2*y9)*y7**2+y7*y8**2*y9**2)*y1+((-y8*y9**2+ &
          y8**2*y9)*y7**2-y7*y8**2*y9**2)*y2+((-y8**2*y9-y8*y9**2)*y7**2+ &
          y7*y8**2*y9**2)*y3+((y8*y9**2-y8**2*y9)*y7**2-y7*y8**2*y9**2)*y4
      dF(144) = (y7**3*y8*y9+(y8**3*y9+y8*y9**3)*y7)*y1+(y7**3*y8*y9+(y8**3*y9+ &
          y8*y9**3)*y7)*y2+(y7**3*y8*y9+(y8**3*y9+y8*y9**3)*y7)*y3+(y7**3*y8*y9+(y8**3*y9+ &
          y8*y9**3)*y7)*y4
      dF(145) = ((y9+y8)*y7**4+(y8**4+y9**4)*y7+y8**4*y9+y8*y9**4)*y1+((-y8+y9)*y7**4+ &
          (-y8**4-y9**4)*y7+y8**4*y9-y8*y9**4)*y2+((-y9-y8)*y7**4+(y8**4+y9**4)*y7- &
          y8**4*y9-y8*y9**4)*y3+((y8-y9)*y7**4+(-y8**4-y9**4)*y7+y8*y9**4-y8**4*y9)*y4
      dF(146) = (-y7**5-y8**5-y9**5)*y1+(y7**5-y9**5+y8**5)*y2+(-y7**5+y8**5+ &
          y9**5)*y3+(y7**5+y9**5-y8**5)*y4
      s1 = ((-sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y7*y9**3/2._ark+sqrt(3._ark)*y8*y9**3/ &
          2._ark-sqrt(3._ark)*y8**3*y9/2._ark)*y5+((y8+y9/2._ark)*y7**3+(-y9**3/2._ark-y8**3)*y7- &
          y8**3*y9/2._ark+y8*y9**3/2._ark)*y6)*y1+((sqrt(3._ark)*y8**3*y9/2._ark+ &
          sqrt(3._ark)*y7**3*y9/2._ark-sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y8*y9**3/2._ark)*y5+ &
          ((y8-y9/2._ark)*y7**3+(y9**3/2._ark-y8**3)*y7+y8**3*y9/2._ark-y8*y9**3/2._ark)*y6)*y2
      dF(147) = s1+((-sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y8**3*y9/2._ark+ &
          sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y8*y9**3/2._ark)*y5+((-y9/2._ark-y8)*y7**3+ &
          (y8**3+y9**3/2._ark)*y7-y8**3*y9/2._ark+y8*y9**3/2._ark)*y6)*y3+((sqrt(3._ark)*y8**3*y9/ &
          2._ark-sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y8*y9**3/ &
          2._ark)*y5+((y9/2._ark-y8)*y7**3+(-y9**3/2._ark+y8**3)*y7+y8**3*y9/2._ark-y8*y9**3/ &
          2._ark)*y6)*y4
      dF(148) = ((y7**2*y8*y9/2._ark+(y8**2*y9/2._ark-y8*y9**2)*y7)*y5+ &
          (sqrt(3._ark)*y7*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y8*y9/2._ark)*y6)*y1+((-y7**2*y8*y9/ &
          2._ark+(-y8**2*y9/2._ark-y8*y9**2)*y7)*y5+(-sqrt(3._ark)*y7*y8**2*y9/2._ark+ &
          sqrt(3._ark)*y7**2*y8*y9/2._ark)*y6)*y2+((y7**2*y8*y9/2._ark+(-y8**2*y9/2._ark+ &
          y8*y9**2)*y7)*y5+(-sqrt(3._ark)*y7*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y8*y9/ &
          2._ark)*y6)*y3+((-y7**2*y8*y9/2._ark+(y8*y9**2+y8**2*y9/2._ark)*y7)*y5+ &
          (sqrt(3._ark)*y7*y8**2*y9/2._ark+sqrt(3._ark)*y7**2*y8*y9/2._ark)*y6)*y4
      dF(149) = ((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y1+((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y2+((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y3+((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y4
      dF(150) = ((y9+y8)*y7+y8*y9)*y1**4+((y8-y9)*y7-y8*y9)*y2**4+((-y9-y8)*y7+ &
          y8*y9)*y3**4+((-y8+y9)*y7-y8*y9)*y4**4
      dF(151) = (((-y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+ &
          (sqrt(3._ark)*y8**2*y9**2/2._ark-sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y1+(((-y9**2/2._ark+ &
          y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+(sqrt(3._ark)*y8**2*y9**2/2._ark- &
          sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y2+(((-y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/ &
          2._ark)*y5+(sqrt(3._ark)*y8**2*y9**2/2._ark-sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y3+(((- &
          y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+(sqrt(3._ark)*y8**2*y9**2/2._ark- &
          sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y4
      s1 = (((y9/2._ark-y8)*y7**3+(y9**3/2._ark-y8**3)*y7+y8**3*y9/2._ark+y8*y9**3/2._ark)*y5+ &
          (-sqrt(3._ark)*y8*y9**3/2._ark+sqrt(3._ark)*y7*y9**3/2._ark+sqrt(3._ark)*y7**3*y9/2._ark- &
          sqrt(3._ark)*y8**3*y9/2._ark)*y6)*y1+(((-y9/2._ark-y8)*y7**3+(-y9**3/2._ark-y8**3)*y7- &
          y8**3*y9/2._ark-y8*y9**3/2._ark)*y5+(-sqrt(3._ark)*y7*y9**3/2._ark+sqrt(3._ark)*y8**3*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**3/2._ark-sqrt(3._ark)*y7**3*y9/2._ark)*y6)*y2
      dF(152) = s1+(((y8-y9/2._ark)*y7**3+(-y9**3/2._ark+y8**3)*y7+y8**3*y9/2._ark+y8*y9**3/ &
          2._ark)*y5+(-sqrt(3._ark)*y8*y9**3/2._ark-sqrt(3._ark)*y7*y9**3/2._ark- &
          sqrt(3._ark)*y8**3*y9/2._ark-sqrt(3._ark)*y7**3*y9/2._ark)*y6)*y3+(((y8+y9/2._ark)*y7**3+ &
          (y8**3+y9**3/2._ark)*y7-y8**3*y9/2._ark-y8*y9**3/2._ark)*y5+(sqrt(3._ark)*y7*y9**3/2._ark+ &
          sqrt(3._ark)*y8**3*y9/2._ark+sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y8*y9**3/ &
          2._ark)*y6)*y4
      s1 = (((-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark)*y7**2+(sqrt(3._ark)*y8**2/3._ark+ &
          sqrt(3._ark)*y9**2/3._ark)*y7+sqrt(3._ark)*y8*y9**2/3._ark-sqrt(3._ark)*y8**2*y9/ &
          6._ark)*y5**2+(-y7**2*y8+(y9**2+y8**2)*y7-y8*y9**2)*y6*y5+(sqrt(3._ark)*y7**2*y9/ &
          2._ark+sqrt(3._ark)*y8**2*y9/2._ark)*y6**2)*y1+(((-sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/ &
          3._ark)*y7**2+(-sqrt(3._ark)*y9**2/3._ark-sqrt(3._ark)*y8**2/3._ark)*y7- &
          sqrt(3._ark)*y8*y9**2/3._ark-sqrt(3._ark)*y8**2*y9/6._ark)*y5**2+(y7**2*y8+(-y9**2- &
          y8**2)*y7+y8*y9**2)*y6*y5+(sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y8**2*y9/ &
          2._ark)*y6**2)*y2
      dF(153) = s1+(((sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/3._ark)*y7**2+(sqrt(3._ark)*y8**2/ &
          3._ark+sqrt(3._ark)*y9**2/3._ark)*y7+sqrt(3._ark)*y8**2*y9/6._ark-sqrt(3._ark)*y8*y9**2/ &
          3._ark)*y5**2+(y7**2*y8+(y9**2+y8**2)*y7+y8*y9**2)*y6*y5+(-sqrt(3._ark)*y8**2*y9/ &
          2._ark-sqrt(3._ark)*y7**2*y9/2._ark)*y6**2)*y3+(((sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y9/ &
          6._ark)*y7**2+(-sqrt(3._ark)*y9**2/3._ark-sqrt(3._ark)*y8**2/3._ark)*y7+ &
          sqrt(3._ark)*y8*y9**2/3._ark+sqrt(3._ark)*y8**2*y9/6._ark)*y5**2+(-y7**2*y8+(-y9**2- &
          y8**2)*y7-y8*y9**2)*y6*y5+(-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y6**2)*y4
      s1 = ((sqrt(3._ark)*y8/6._ark-5._ark/24._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/6._ark)*y5**4+(- &
          sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/6._ark+7._ark/12._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(y7- &
          y8)*y6**3*y5+sqrt(3._ark)*y6**4*y9/8._ark)*y1+((-sqrt(3._ark)*y8/6._ark-5._ark/ &
          24._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark)*y5**4+(7._ark/12._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y8/6._ark)*y6**2*y5**2+(y8-y7)*y6**3*y5+ &
          sqrt(3._ark)*y6**4*y9/8._ark)*y2
      dF(154) = s1+((-sqrt(3._ark)*y8/6._ark+5._ark/24._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/ &
          6._ark)*y5**4+(-7._ark/12._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y8/ &
          6._ark)*y6**2*y5**2+(y7+y8)*y6**3*y5-sqrt(3._ark)*y6**4*y9/8._ark)*y3+((sqrt(3._ark)*y8/ &
          6._ark-sqrt(3._ark)*y7/6._ark+5._ark/24._ark*sqrt(3._ark)*y9)*y5**4+(sqrt(3._ark)*y7/6._ark- &
          sqrt(3._ark)*y8/6._ark-7._ark/12._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(-y8-y7)*y6**3*y5- &
          sqrt(3._ark)*y6**4*y9/8._ark)*y4
      dF(155) = (((-y9-y8)*y7-y8*y9)*y5**2+((-y9-y8)*y7-y8*y9)*y6**2)*y1**2+(((-y8+ &
          y9)*y7+y8*y9)*y5**2+((-y8+y9)*y7+y8*y9)*y6**2)*y2**2+(((y9+y8)*y7-y8*y9)*y5**2+ &
          ((y9+y8)*y7-y8*y9)*y6**2)*y3**2+(((y8-y9)*y7+y8*y9)*y5**2+((y8-y9)*y7+ &
          y8*y9)*y6**2)*y4**2
      dF(156) = (y6**4+y5**4+2._ark*y5**2*y6**2)*y1**2+(y6**4+y5**4+ &
          2._ark*y5**2*y6**2)*y2**2+(y6**4+y5**4+2._ark*y5**2*y6**2)*y3**2+(y6**4+y5**4+ &
          2._ark*y5**2*y6**2)*y4**2
      dF(157) = ((-y7**3-y8**3)*y2+(-y9**3-y8**3)*y3+(-y7**3-y9**3)*y4)*y1**2+((y8**3+ &
          y7**3)*y2**2+(y9**3+y8**3)*y3**2+(y9**3+y7**3)*y4**2)*y1+((-y9**3+y7**3)*y3+ &
          (y8**3-y9**3)*y4)*y2**2+((-y7**3+y9**3)*y3**2+(y9**3-y8**3)*y4**2)*y2+(-y7**3+ &
          y8**3)*y4*y3**2+(-y8**3+y7**3)*y4**2*y3
      dF(158) = ((y8*y9**2+y7*y9**2)*y2+(y9+y8)*y7**2*y3+(y8**2*y9+ &
          y7*y8**2)*y4)*y1**2+((-y8*y9**2-y7*y9**2)*y2**2+(-y9-y8)*y7**2*y3**2+(-y8**2*y9- &
          y7*y8**2)*y4**2)*y1+((-y7*y8**2+y8**2*y9)*y3+(-y8+y9)*y7**2*y4)*y2**2+((- &
          y8**2*y9+y7*y8**2)*y3**2+(y8-y9)*y7**2*y4**2)*y2+(-y8*y9**2+y7*y9**2)*y4*y3**2+ &
          (-y7*y9**2+y8*y9**2)*y4**2*y3
      dF(159) = (y4*y7*y8*y9+y3*y7*y8*y9+y2*y7*y8*y9)*y1**2+(y4**2*y7*y8*y9+ &
          y3**2*y7*y8*y9+y2**2*y7*y8*y9)*y1+(y4*y7*y8*y9+y3*y7*y8*y9)*y2**2+ &
          (y4**2*y7*y8*y9+y3**2*y7*y8*y9)*y2+y3**2*y4*y7*y8*y9+y3*y4**2*y7*y8*y9
      dF(160) = ((y7*y8**2+y7**2*y8)*y2+(y8*y9**2+y8**2*y9)*y3+(y7*y9**2+ &
          y7**2*y9)*y4)*y1**2+((-y7**2*y8-y7*y8**2)*y2**2+(-y8**2*y9-y8*y9**2)*y3**2+(- &
          y7*y9**2-y7**2*y9)*y4**2)*y1+((-y7*y9**2+y7**2*y9)*y3+(-y8*y9**2+ &
          y8**2*y9)*y4)*y2**2+((y7*y9**2-y7**2*y9)*y3**2+(y8*y9**2-y8**2*y9)*y4**2)*y2+ &
          (y7*y8**2-y7**2*y8)*y4*y3**2+(-y7*y8**2+y7**2*y8)*y4**2*y3
      dF(161) = ((-y8**2*y9-y7**2*y9)*y2+(-y9**2-y8**2)*y7*y3+(-y8*y9**2- &
          y7**2*y8)*y4)*y1**2+((-y8**2*y9-y7**2*y9)*y2**2+(-y9**2-y8**2)*y7*y3**2+(- &
          y8*y9**2-y7**2*y8)*y4**2)*y1+((y7**2*y8+y8*y9**2)*y3+(y9**2+y8**2)*y7*y4)*y2**2+ &
          ((y7**2*y8+y8*y9**2)*y3**2+(y9**2+y8**2)*y7*y4**2)*y2+(y7**2*y9+ &
          y8**2*y9)*y4*y3**2+(y7**2*y9+y8**2*y9)*y4**2*y3
      dF(162) = (y4**2*y8**2+y3**2*y7**2+y2**2*y9**2)*y1**2+(y4**2*y7**2+ &
          y3**2*y8**2)*y2**2+y3**2*y4**2*y9**2
      dF(163) = ((y8**2+y7**2)*y2**2+(y9**2+y8**2)*y3**2+(y7**2+y9**2)*y4**2)*y1**2+ &
          ((y7**2+y9**2)*y3**2+(y9**2+y8**2)*y4**2)*y2**2+(y8**2+y7**2)*y4**2*y3**2
      dF(164) = ((-y9-y8)*y7**2+(-y9**2-y8**2)*y7-y8*y9**2-y8**2*y9)*y1**3+((y8- &
          y9)*y7**2+(y9**2+y8**2)*y7-y8**2*y9+y8*y9**2)*y2**3+((y9+y8)*y7**2+(-y9**2- &
          y8**2)*y7+y8**2*y9+y8*y9**2)*y3**3+((-y8+y9)*y7**2+(y9**2+y8**2)*y7-y8*y9**2+ &
          y8**2*y9)*y4**3
      dF(165) = (-y8**3-y7**3-y9**3)*y1**3+(-y9**3+y8**3+y7**3)*y2**3+(y8**3+y9**3- &
          y7**3)*y3**3+(y7**3+y9**3-y8**3)*y4**3
      dF(166) = y4**3*y7*y8*y9+y2**3*y7*y8*y9+y1**3*y7*y8*y9+y3**3*y7*y8*y9
      dF(167) = ((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(-y7**2+y8**2)*y6)*y1**3+((sqrt(3._ark)*y8**2/3._ark+ &
          sqrt(3._ark)*y7**2/3._ark-2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+(-y7**2+y8**2)*y6)*y2**3+ &
          ((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark-2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+(- &
          y7**2+y8**2)*y6)*y3**3+((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(-y7**2+y8**2)*y6)*y4**3
      dF(168) = ((2._ark/3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y5+ &
          (y7-y8)*y6)*y1**4+((sqrt(3._ark)*y8/3._ark+2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/ &
          3._ark)*y5+(y8-y7)*y6)*y2**4+((-sqrt(3._ark)*y7/3._ark+sqrt(3._ark)*y8/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5+(y7+y8)*y6)*y3**4+((-sqrt(3._ark)*y8/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/3._ark)*y5+(-y8-y7)*y6)*y4**4
      s1 = (((-sqrt(3._ark)*y9/3._ark-sqrt(3._ark)*y8/3._ark)*y7**2+(-sqrt(3._ark)*y8**2/3._ark+ &
          sqrt(3._ark)*y9**2/6._ark)*y7+sqrt(3._ark)*y8*y9**2/6._ark-sqrt(3._ark)*y8**2*y9/ &
          3._ark)*y5**2+((-y9-y8)*y7**2+y7*y8**2+y8**2*y9)*y6*y5+(-sqrt(3._ark)*y8*y9**2/2._ark- &
          sqrt(3._ark)*y7*y9**2/2._ark)*y6**2)*y1+(((-sqrt(3._ark)*y9/3._ark+sqrt(3._ark)*y8/ &
          3._ark)*y7**2+(sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y9**2/6._ark)*y7- &
          sqrt(3._ark)*y8**2*y9/3._ark-sqrt(3._ark)*y8*y9**2/6._ark)*y5**2+((y8-y9)*y7**2- &
          y7*y8**2+y8**2*y9)*y6*y5+(sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8*y9**2/ &
          2._ark)*y6**2)*y2
      dF(169) = s1+(((sqrt(3._ark)*y9/3._ark+sqrt(3._ark)*y8/3._ark)*y7**2+(-sqrt(3._ark)*y8**2/ &
          3._ark+sqrt(3._ark)*y9**2/6._ark)*y7-sqrt(3._ark)*y8*y9**2/6._ark+sqrt(3._ark)*y8**2*y9/ &
          3._ark)*y5**2+((y9+y8)*y7**2+y7*y8**2-y8**2*y9)*y6*y5+(sqrt(3._ark)*y8*y9**2/2._ark- &
          sqrt(3._ark)*y7*y9**2/2._ark)*y6**2)*y3+(((-sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y9/ &
          3._ark)*y7**2+(sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y9**2/6._ark)*y7+ &
          sqrt(3._ark)*y8*y9**2/6._ark+sqrt(3._ark)*y8**2*y9/3._ark)*y5**2+((-y8+y9)*y7**2- &
          y7*y8**2-y8**2*y9)*y6*y5+(-sqrt(3._ark)*y8*y9**2/2._ark+sqrt(3._ark)*y7*y9**2/ &
          2._ark)*y6**2)*y4
      dF(170) = ((-sqrt(3._ark)*y9**3/6._ark+sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/ &
          3._ark)*y5**2+(-y8**3+y7**3)*y6*y5+sqrt(3._ark)*y6**2*y9**3/2._ark)*y1+((- &
          sqrt(3._ark)*y9**3/6._ark-sqrt(3._ark)*y8**3/3._ark-sqrt(3._ark)*y7**3/3._ark)*y5**2+(- &
          y7**3+y8**3)*y6*y5+sqrt(3._ark)*y6**2*y9**3/2._ark)*y2+((sqrt(3._ark)*y7**3/3._ark- &
          sqrt(3._ark)*y8**3/3._ark+sqrt(3._ark)*y9**3/6._ark)*y5**2+(y8**3+y7**3)*y6*y5- &
          sqrt(3._ark)*y6**2*y9**3/2._ark)*y3+((-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/3._ark+ &
          sqrt(3._ark)*y9**3/6._ark)*y5**2+(-y7**3-y8**3)*y6*y5-sqrt(3._ark)*y6**2*y9**3/ &
          2._ark)*y4
      dF(171) = ((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y1+((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y2+((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y3+((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y4
      dF(172) = (((y8/3._ark+y9/3._ark)*y7+y8*y9/3._ark)*y5**3+((-y9-y8)*y7- &
          y8*y9)*y6**2*y5)*y1+(((y8/3._ark-y9/3._ark)*y7-y8*y9/3._ark)*y5**3+((-y8+y9)*y7+ &
          y8*y9)*y6**2*y5)*y2+(((-y8/3._ark-y9/3._ark)*y7+y8*y9/3._ark)*y5**3+((y9+y8)*y7- &
          y8*y9)*y6**2*y5)*y3+(((y9/3._ark-y8/3._ark)*y7-y8*y9/3._ark)*y5**3+((y8-y9)*y7+ &
          y8*y9)*y6**2*y5)*y4
      dF(173) = ((y7**2*y9**2+y8**2*y9**2)*y2+(y9**2+y8**2)*y7**2*y3+(y8**2*y9**2+ &
          y7**2*y8**2)*y4)*y1+((y8**2*y9**2+y7**2*y8**2)*y3+(y9**2+y8**2)*y7**2*y4)*y2+ &
          (y7**2*y9**2+y8**2*y9**2)*y4*y3
      dF(174) = (((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5**2+(- &
          sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y2+((-y9**2-y8**2)*y6*y5+ &
          (sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y3+((y7**2+y9**2)*y6*y5+ &
          (sqrt(3._ark)*y7**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y4)*y1+(((y7**2+ &
          y9**2)*y6*y5+(sqrt(3._ark)*y7**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y3+((-y9**2- &
          y8**2)*y6*y5+(sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y4)*y2+ &
          ((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5**2+(-sqrt(3._ark)*y8**2/6._ark- &
          sqrt(3._ark)*y7**2/6._ark)*y6**2)*y4*y3
      dF(175) = (((y8**2+y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y2+((y9**2+y8**2)*y5**2+ &
          (y9**2+y8**2)*y6**2)*y3+((y7**2+y9**2)*y5**2+(y7**2+y9**2)*y6**2)*y4)*y1+ &
          (((y7**2+y9**2)*y5**2+(y7**2+y9**2)*y6**2)*y3+((y9**2+y8**2)*y5**2+(y9**2+ &
          y8**2)*y6**2)*y4)*y2+((y8**2+y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y4*y3
      s1 = (((-sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y5+(-y7**2+y8**2)*y6)*y2+ &
          (sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y3+(sqrt(3._ark)*y5*y7**2+ &
          (2._ark*y9**2+y7**2)*y6)*y4)*y1**2+(((-sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y5+(- &
          y7**2+y8**2)*y6)*y2**2+(sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y3**2+ &
          (sqrt(3._ark)*y5*y7**2+(2._ark*y9**2+y7**2)*y6)*y4**2)*y1+((sqrt(3._ark)*y5*y7**2+ &
          (2._ark*y9**2+y7**2)*y6)*y3+(sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y4)*y2**2
      dF(176) = s1+((sqrt(3._ark)*y5*y7**2+(2._ark*y9**2+y7**2)*y6)*y3**2+ &
          (sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y4**2)*y2+((-sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y5+(-y7**2+y8**2)*y6)*y4*y3**2+((-sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y5+(-y7**2+y8**2)*y6)*y4**2*y3
      dF(177) = (y4+y3+y2)*y1**5+(y2**5+y4**5+y3**5)*y1+(y4+y3)*y2**5+(y4**5+ &
          y3**5)*y2+y3**5*y4+y3*y4**5
      s2 = (((-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y5**2+(y7-y8)*y6*y5)*y2+ &
          ((sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/3._ark)*y5**2+y5*y6*y9-sqrt(3._ark)*y6**2*y8/ &
          2._ark)*y3+((sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y9/3._ark)*y5**2-y5*y6*y9- &
          sqrt(3._ark)*y6**2*y7/2._ark)*y4)*y1**2
      s3 = (((sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5**2+(y8-y7)*y6*y5)*y2**2+((- &
          sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2-y5*y6*y9+sqrt(3._ark)*y6**2*y8/ &
          2._ark)*y3**2+((-sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2+y5*y6*y9+ &
          sqrt(3._ark)*y6**2*y7/2._ark)*y4**2)*y1+(((-sqrt(3._ark)*y9/3._ark-sqrt(3._ark)*y7/ &
          6._ark)*y5**2-y5*y6*y9+sqrt(3._ark)*y6**2*y7/2._ark)*y3+((-sqrt(3._ark)*y9/3._ark- &
          sqrt(3._ark)*y8/6._ark)*y5**2+y5*y6*y9+sqrt(3._ark)*y6**2*y8/2._ark)*y4)*y2**2
      s1 = s2+s3
      dF(178) = s1+(((sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2+y5*y6*y9- &
          sqrt(3._ark)*y6**2*y7/2._ark)*y3**2+((sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2- &
          y5*y6*y9-sqrt(3._ark)*y6**2*y8/2._ark)*y4**2)*y2+((sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/ &
          3._ark)*y5**2+(y7+y8)*y6*y5)*y4*y3**2+((sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/ &
          3._ark)*y5**2+(-y8-y7)*y6*y5)*y4**2*y3
      s2 = (((-sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/6._ark)*y5**2+(y7-y8)*y6*y5+(- &
          sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y6**2)*y2+((-2._ark/3._ark*sqrt(3._ark)*y9- &
          sqrt(3._ark)*y8/6._ark)*y5**2-y5*y6*y8-sqrt(3._ark)*y6**2*y8/2._ark)*y3+((-2._ark/ &
          3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark)*y5**2+y5*y6*y7-sqrt(3._ark)*y6**2*y7/ &
          2._ark)*y4)*y1**2
      s3 = (((sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y7/6._ark)*y5**2+(y8-y7)*y6*y5+ &
          (sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6**2)*y2**2+((2._ark/3._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y8/6._ark)*y5**2+y5*y6*y8+sqrt(3._ark)*y6**2*y8/2._ark)*y3**2+ &
          ((sqrt(3._ark)*y7/6._ark+2._ark/3._ark*sqrt(3._ark)*y9)*y5**2-y5*y6*y7+ &
          sqrt(3._ark)*y6**2*y7/2._ark)*y4**2)*y1+(((sqrt(3._ark)*y7/6._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5**2-y5*y6*y7+sqrt(3._ark)*y6**2*y7/2._ark)*y3+((sqrt(3._ark)*y8/ &
          6._ark-2._ark/3._ark*sqrt(3._ark)*y9)*y5**2+y5*y6*y8+sqrt(3._ark)*y6**2*y8/2._ark)*y4)*y2**2
      s1 = s2+s3
      dF(179) = s1+(((-sqrt(3._ark)*y7/6._ark+2._ark/3._ark*sqrt(3._ark)*y9)*y5**2+y5*y6*y7- &
          sqrt(3._ark)*y6**2*y7/2._ark)*y3**2+((-sqrt(3._ark)*y8/6._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5**2-y5*y6*y8-sqrt(3._ark)*y6**2*y8/2._ark)*y4**2)*y2+((- &
          sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y8/6._ark)*y5**2+(y7+y8)*y6*y5+(-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6**2)*y4*y3**2+((sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/ &
          6._ark)*y5**2+(-y8-y7)*y6*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6**2)*y4**2*y3
      dF(180) = ((y4*y7**2+y3*y8**2)*y2+y3*y4*y9**2)*y1**2+((y4*y8**2+y3*y7**2)*y2**2+ &
          (y4**2*y9**2+y3**2*y9**2)*y2+y3**2*y4*y8**2+y3*y4**2*y7**2)*y1+ &
          y2**2*y3*y4*y9**2+(y3*y4**2*y8**2+y3**2*y4*y7**2)*y2
      dF(181) = (((-y7*y8-y8*y9)*y3+(-y9-y8)*y7*y4)*y2+(-y8*y9-y7*y9)*y4*y3)*y1**2+ &
          (((-y8+y9)*y7*y3+(-y7*y8+y8*y9)*y4)*y2**2+((-y8*y9+y7*y9)*y3**2+(y8*y9- &
          y7*y9)*y4**2)*y2+(-y8*y9+y7*y8)*y4*y3**2+(y8-y9)*y7*y4**2*y3)*y1+(y8*y9+ &
          y7*y9)*y4*y3*y2**2+((y9+y8)*y7*y4*y3**2+(y8*y9+y7*y8)*y4**2*y3)*y2
      dF(182) = (y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y1+(y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y2+ &
          (y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y3+(y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y4
      dF(183) = (((y9+y8)*y7**2+(y9**2+y8**2)*y7+y8*y9**2+y8**2*y9)*y5**2+((y9+ &
          y8)*y7**2+(y9**2+y8**2)*y7+y8*y9**2+y8**2*y9)*y6**2)*y1+(((-y8+y9)*y7**2+(- &
          y9**2-y8**2)*y7+y8**2*y9-y8*y9**2)*y5**2+((-y8+y9)*y7**2+(-y9**2-y8**2)*y7+ &
          y8**2*y9-y8*y9**2)*y6**2)*y2+(((-y9-y8)*y7**2+(y9**2+y8**2)*y7-y8**2*y9- &
          y8*y9**2)*y5**2+((-y9-y8)*y7**2+(y9**2+y8**2)*y7-y8**2*y9-y8*y9**2)*y6**2)*y3+ &
          (((y8-y9)*y7**2+(-y9**2-y8**2)*y7+y8*y9**2-y8**2*y9)*y5**2+((y8-y9)*y7**2+(- &
          y9**2-y8**2)*y7+y8*y9**2-y8**2*y9)*y6**2)*y4
      dF(184) = ((y9**3+y8**3+y7**3)*y5**2+(y9**3+y8**3+y7**3)*y6**2)*y1+((y9**3- &
          y8**3-y7**3)*y5**2+(y9**3-y8**3-y7**3)*y6**2)*y2+((-y8**3-y9**3+y7**3)*y5**2+(- &
          y8**3-y9**3+y7**3)*y6**2)*y3+((-y7**3-y9**3+y8**3)*y5**2+(-y7**3-y9**3+ &
          y8**3)*y6**2)*y4
      dF(185) = (((4._ark/9._ark*sqrt(3._ark)*y9-5._ark/9._ark*sqrt(3._ark)*y8)*y7+4._ark/ &
          9._ark*sqrt(3._ark)*y8*y9)*y5**3+(-y8*y9+y7*y9)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y7*y8+ &
          (-y8*y9+y7*y9)*y6**3)*y1+(((-5._ark/9._ark*sqrt(3._ark)*y8-4._ark/ &
          9._ark*sqrt(3._ark)*y9)*y7-4._ark/9._ark*sqrt(3._ark)*y8*y9)*y5**3+(y8*y9-y7*y9)*y6*y5**2- &
          sqrt(3._ark)*y5*y6**2*y7*y8+(y8*y9-y7*y9)*y6**3)*y2+(((-4._ark/9._ark*sqrt(3._ark)*y9+ &
          5._ark/9._ark*sqrt(3._ark)*y8)*y7+4._ark/9._ark*sqrt(3._ark)*y8*y9)*y5**3+(-y8*y9- &
          y7*y9)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y7*y8+(-y8*y9-y7*y9)*y6**3)*y3+(((4._ark/ &
          9._ark*sqrt(3._ark)*y9+5._ark/9._ark*sqrt(3._ark)*y8)*y7-4._ark/ &
          9._ark*sqrt(3._ark)*y8*y9)*y5**3+(y8*y9+y7*y9)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y7*y8+ &
          (y8*y9+y7*y9)*y6**3)*y4
      dF(186) = ((-4._ark/9._ark*sqrt(3._ark)*y8**2+5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/ &
          9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2-y8**2)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9**2+ &
          (y7**2-y8**2)*y6**3)*y1+((-4._ark/9._ark*sqrt(3._ark)*y8**2+5._ark/ &
          9._ark*sqrt(3._ark)*y9**2-4._ark/9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2-y8**2)*y6*y5**2+ &
          sqrt(3._ark)*y5*y6**2*y9**2+(y7**2-y8**2)*y6**3)*y2+((-4._ark/9._ark*sqrt(3._ark)*y8**2+ &
          5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2- &
          y8**2)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9**2+(y7**2-y8**2)*y6**3)*y3+((-4._ark/ &
          9._ark*sqrt(3._ark)*y8**2+5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/ &
          9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2-y8**2)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9**2+ &
          (y7**2-y8**2)*y6**3)*y4
      dF(187) = ((y9+y7+y8)*y5**4+(2._ark*y8+2._ark*y9+2._ark*y7)*y6**2*y5**2+(y9+y7+ &
          y8)*y6**4)*y1+((-y7+y9-y8)*y5**4+(-2._ark*y7+2._ark*y9-2._ark*y8)*y6**2*y5**2+(-y7+y9- &
          y8)*y6**4)*y2+((-y8-y9+y7)*y5**4+(-2._ark*y8+2._ark*y7-2._ark*y9)*y6**2*y5**2+(-y8-y9+ &
          y7)*y6**4)*y3+((-y7-y9+y8)*y5**4+(-2._ark*y9+2._ark*y8-2._ark*y7)*y6**2*y5**2+(-y7-y9+ &
          y8)*y6**4)*y4
      s1 = ((sqrt(3._ark)*y9/24._ark+sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y7/6._ark)*y5**4+(y7- &
          y8)*y6*y5**3+(sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/4._ark+sqrt(3._ark)*y7/ &
          2._ark)*y6**2*y5**2+3._ark/8._ark*sqrt(3._ark)*y6**4*y9)*y1+((-sqrt(3._ark)*y8/6._ark- &
          sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/24._ark)*y5**4+(y8-y7)*y6*y5**3+(-sqrt(3._ark)*y8/ &
          2._ark-sqrt(3._ark)*y9/4._ark-sqrt(3._ark)*y7/2._ark)*y6**2*y5**2+3._ark/ &
          8._ark*sqrt(3._ark)*y6**4*y9)*y2
      dF(188) = s1+((-sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/24._ark+sqrt(3._ark)*y7/ &
          6._ark)*y5**4+(y7+y8)*y6*y5**3+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark+ &
          sqrt(3._ark)*y9/4._ark)*y6**2*y5**2-3._ark/8._ark*sqrt(3._ark)*y6**4*y9)*y3+ &
          ((sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/24._ark-sqrt(3._ark)*y7/6._ark)*y5**4+(-y8- &
          y7)*y6*y5**3+(sqrt(3._ark)*y9/4._ark-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
          2._ark)*y6**2*y5**2-3._ark/8._ark*sqrt(3._ark)*y6**4*y9)*y4
      dF(189) = (-3._ark*y5*y6**4+y5**5-2._ark*y5**3*y6**2)*y1+(-3._ark*y5*y6**4+y5**5- &
          2._ark*y5**3*y6**2)*y2+(-3._ark*y5*y6**4+y5**5-2._ark*y5**3*y6**2)*y3+(-3._ark*y5*y6**4+ &
          y5**5-2._ark*y5**3*y6**2)*y4
      dF(190) = ((sqrt(3._ark)*y5**2*y9**2/2._ark-sqrt(3._ark)*y6**2*y9**2/6._ark)*y2+(- &
          y5*y6*y7**2+sqrt(3._ark)*y6**2*y7**2/3._ark)*y3+(y5*y6*y8**2+sqrt(3._ark)*y6**2*y8**2/ &
          3._ark)*y4)*y1+((y5*y6*y8**2+sqrt(3._ark)*y6**2*y8**2/3._ark)*y3+(-y5*y6*y7**2+ &
          sqrt(3._ark)*y6**2*y7**2/3._ark)*y4)*y2+(sqrt(3._ark)*y5**2*y9**2/2._ark- &
          sqrt(3._ark)*y6**2*y9**2/6._ark)*y4*y3
      dF(191) = ((sqrt(3._ark)*y5**2*y7*y8/2._ark-sqrt(3._ark)*y6**2*y7*y8/6._ark)*y2+ &
          (sqrt(3._ark)*y6**2*y8*y9/3._ark-y5*y6*y8*y9)*y3+(sqrt(3._ark)*y6**2*y7*y9/3._ark+ &
          y5*y6*y7*y9)*y4)*y1+((-sqrt(3._ark)*y6**2*y7*y9/3._ark-y5*y6*y7*y9)*y3+(- &
          sqrt(3._ark)*y6**2*y8*y9/3._ark+y5*y6*y8*y9)*y4)*y2+(sqrt(3._ark)*y6**2*y7*y8/6._ark- &
          sqrt(3._ark)*y5**2*y7*y8/2._ark)*y4*y3
      dF(192) = ((y8*y7*y5**2+y8*y7*y6**2)*y2+(y5**2*y8*y9+y9*y8*y6**2)*y3+ &
          (y5**2*y7*y9+y9*y7*y6**2)*y4)*y1+((-y9*y7*y6**2-y5**2*y7*y9)*y3+(-y5**2*y8*y9- &
          y9*y8*y6**2)*y4)*y2+(-y8*y7*y5**2-y8*y7*y6**2)*y4*y3
      dF(193) = ((2._ark/9._ark*sqrt(3._ark)*y6**4+2._ark*sqrt(3._ark)*y5**2*y6**2)*y2+(5._ark/ &
          3._ark*y5*y6**3+7._ark/18._ark*sqrt(3._ark)*y6**4-y5**3*y6+sqrt(3._ark)*y5**4/2._ark)*y3+ &
          (7._ark/18._ark*sqrt(3._ark)*y6**4+sqrt(3._ark)*y5**4/2._ark+y5**3*y6-5._ark/ &
          3._ark*y5*y6**3)*y4)*y1+((7._ark/18._ark*sqrt(3._ark)*y6**4+sqrt(3._ark)*y5**4/2._ark+ &
          y5**3*y6-5._ark/3._ark*y5*y6**3)*y3+(5._ark/3._ark*y5*y6**3+7._ark/18._ark*sqrt(3._ark)*y6**4- &
          y5**3*y6+sqrt(3._ark)*y5**4/2._ark)*y4)*y2+(2._ark/9._ark*sqrt(3._ark)*y6**4+ &
          2._ark*sqrt(3._ark)*y5**2*y6**2)*y4*y3
      dF(194) = (y4*y8**3+y2*y9**3+y3*y7**3)*y1**2+(y2**2*y9**3+y4**2*y8**3+ &
          y3**2*y7**3)*y1+(-y4*y7**3-y3*y8**3)*y2**2+(-y3**2*y8**3-y4**2*y7**3)*y2- &
          y3**2*y4*y9**3-y3*y4**2*y9**3
      s1 = (((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5+(-y7*y9/2._ark+y8*y9/ &
          2._ark)*y6)*y2+(-sqrt(3._ark)*y5*y7*y9/2._ark+(y8+y9/2._ark)*y7*y6)*y3+(- &
          sqrt(3._ark)*y5*y8*y9/2._ark+(-y8*y9/2._ark-y7*y8)*y6)*y4)*y1**2+(((-sqrt(3._ark)*y8*y9/ &
          2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5+(y7*y9/2._ark-y8*y9/2._ark)*y6)*y2**2+ &
          (sqrt(3._ark)*y5*y7*y9/2._ark+(-y9/2._ark-y8)*y7*y6)*y3**2+(sqrt(3._ark)*y5*y8*y9/2._ark+ &
          (y7*y8+y8*y9/2._ark)*y6)*y4**2)*y1+((sqrt(3._ark)*y5*y8*y9/2._ark+(y8*y9/2._ark- &
          y7*y8)*y6)*y3+(sqrt(3._ark)*y5*y7*y9/2._ark+(y8-y9/2._ark)*y7*y6)*y4)*y2**2
      dF(195) = s1+((-sqrt(3._ark)*y5*y8*y9/2._ark+(y7*y8-y8*y9/2._ark)*y6)*y3**2+(- &
          sqrt(3._ark)*y5*y7*y9/2._ark+(y9/2._ark-y8)*y7*y6)*y4**2)*y2+((-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y5+(y8*y9/2._ark+y7*y9/2._ark)*y6)*y4*y3**2+((- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5+(-y8*y9/2._ark-y7*y9/ &
          2._ark)*y6)*y4**2*y3
      dF(196) = (((-y8-y7)*y5**2+(-y8-y7)*y6**2)*y2+((-y9-y8)*y5**2+(-y9- &
          y8)*y6**2)*y3+((-y9-y7)*y5**2+(-y9-y7)*y6**2)*y4)*y1**2+(((y7+y8)*y5**2+(y7+ &
          y8)*y6**2)*y2**2+((y9+y8)*y5**2+(y9+y8)*y6**2)*y3**2+((y7+y9)*y5**2+(y7+ &
          y9)*y6**2)*y4**2)*y1+(((y7-y9)*y5**2+(y7-y9)*y6**2)*y3+((y8-y9)*y5**2+(y8- &
          y9)*y6**2)*y4)*y2**2+(((y9-y7)*y5**2+(y9-y7)*y6**2)*y3**2+((-y8+y9)*y5**2+(-y8+ &
          y9)*y6**2)*y4**2)*y2+((y8-y7)*y5**2+(y8-y7)*y6**2)*y4*y3**2+((y7-y8)*y5**2+(y7- &
          y8)*y6**2)*y4**2*y3
      s1 = (((y8**2+y7**2)*y5+(sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y6)*y2+((y9**2- &
          2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y3+((-2._ark*y7**2+y9**2)*y5- &
          sqrt(3._ark)*y6*y9**2)*y4)*y1**2+(((y8**2+y7**2)*y5+(sqrt(3._ark)*y7**2- &
          sqrt(3._ark)*y8**2)*y6)*y2**2+((y9**2-2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y3**2+ &
          ((-2._ark*y7**2+y9**2)*y5-sqrt(3._ark)*y6*y9**2)*y4**2)*y1+(((-2._ark*y7**2+y9**2)*y5- &
          sqrt(3._ark)*y6*y9**2)*y3+((y9**2-2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y4)*y2**2
      dF(197) = s1+(((-2._ark*y7**2+y9**2)*y5-sqrt(3._ark)*y6*y9**2)*y3**2+((y9**2- &
          2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y4**2)*y2+((y8**2+y7**2)*y5+ &
          (sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y6)*y4*y3**2+((y8**2+y7**2)*y5+ &
          (sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y6)*y4**2*y3
      dF(198) = (-2._ark*y2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y3+(y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y4)*y1**2+(-2._ark*y2**2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+ &
          y5*y7**2)*y3**2+(y5*y8**2+sqrt(3._ark)*y6*y8**2)*y4**2)*y1+((y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y3+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4)*y2**2+((y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y3**2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4**2)*y2- &
          2._ark*y3*y4**2*y5*y9**2-2._ark*y3**2*y4*y5*y9**2
      s1 = (((y8*y9/2._ark+y7*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/ &
          2._ark)*y6)*y2+((y9/2._ark-y8)*y7*y5+sqrt(3._ark)*y6*y7*y9/2._ark)*y3+((y8*y9/2._ark- &
          y7*y8)*y5-sqrt(3._ark)*y6*y8*y9/2._ark)*y4)*y1**2+(((-y8*y9/2._ark-y7*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/2._ark)*y6)*y2**2+((y8-y9/2._ark)*y7*y5- &
          sqrt(3._ark)*y6*y7*y9/2._ark)*y3**2+((y7*y8-y8*y9/2._ark)*y5+sqrt(3._ark)*y6*y8*y9/ &
          2._ark)*y4**2)*y1+(((-y8*y9/2._ark-y7*y8)*y5+sqrt(3._ark)*y6*y8*y9/2._ark)*y3+((-y9/ &
          2._ark-y8)*y7*y5-sqrt(3._ark)*y6*y7*y9/2._ark)*y4)*y2**2
      dF(199) = s1+(((y7*y8+y8*y9/2._ark)*y5-sqrt(3._ark)*y6*y8*y9/2._ark)*y3**2+((y8+y9/ &
          2._ark)*y7*y5+sqrt(3._ark)*y6*y7*y9/2._ark)*y4**2)*y2+((-y7*y9/2._ark+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4*y3**2+((y7*y9/2._ark-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4**2*y3
      dF(200) = ((y9*y6**2+y9*y5**2)*y2+(y7*y6**2+y5**2*y7)*y3+(y8*y6**2+ &
          y5**2*y8)*y4)*y1**2+((y9*y6**2+y9*y5**2)*y2**2+(y7*y6**2+y5**2*y7)*y3**2+ &
          (y8*y6**2+y5**2*y8)*y4**2)*y1+((-y5**2*y8-y8*y6**2)*y3+(-y5**2*y7- &
          y7*y6**2)*y4)*y2**2+((-y5**2*y8-y8*y6**2)*y3**2+(-y5**2*y7-y7*y6**2)*y4**2)*y2+ &
          (-y9*y6**2-y9*y5**2)*y4*y3**2+(-y9*y6**2-y9*y5**2)*y4**2*y3
      dF(201) = ((y5**3-3._ark*y5*y6**2)*y2+(y5**3-3._ark*y5*y6**2)*y3+(y5**3- &
          3._ark*y5*y6**2)*y4)*y1**2+((y5**3-3._ark*y5*y6**2)*y2**2+(y5**3- &
          3._ark*y5*y6**2)*y3**2+(y5**3-3._ark*y5*y6**2)*y4**2)*y1+((y5**3-3._ark*y5*y6**2)*y3+ &
          (y5**3-3._ark*y5*y6**2)*y4)*y2**2+((y5**3-3._ark*y5*y6**2)*y3**2+(y5**3- &
          3._ark*y5*y6**2)*y4**2)*y2+(y5**3-3._ark*y5*y6**2)*y4*y3**2+(y5**3- &
          3._ark*y5*y6**2)*y4**2*y3
      dF(202) = ((y8**2+y7**2)*y2+(y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y1**3+((y8**2+ &
          y7**2)*y2**3+(y9**2+y8**2)*y3**3+(y7**2+y9**2)*y4**3)*y1+((y7**2+y9**2)*y3+ &
          (y9**2+y8**2)*y4)*y2**3+((y7**2+y9**2)*y3**3+(y9**2+y8**2)*y4**3)*y2+(y8**2+ &
          y7**2)*y4*y3**3+(y8**2+y7**2)*y4**3*y3
      dF(203) = ((-y8*y9-y7*y9)*y2+(-y9-y8)*y7*y3+(-y7*y8-y8*y9)*y4)*y1**3+((y8*y9+ &
          y7*y9)*y2**3+(y9+y8)*y7*y3**3+(y8*y9+y7*y8)*y4**3)*y1+((-y7*y8+y8*y9)*y3+(-y8+ &
          y9)*y7*y4)*y2**3+((-y8*y9+y7*y8)*y3**3+(y8-y9)*y7*y4**3)*y2+(-y8*y9+ &
          y7*y9)*y4*y3**3+(y8*y9-y7*y9)*y4**3*y3
      dF(204) = (y3*y7**2+y2*y9**2+y4*y8**2)*y1**3+(y3**3*y7**2+y4**3*y8**2+ &
          y2**3*y9**2)*y1+(y4*y7**2+y3*y8**2)*y2**3+(y3**3*y8**2+y4**3*y7**2)*y2+ &
          y3**3*y4*y9**2+y3*y4**3*y9**2
      dF(205) = (y4*y7*y9+y3*y8*y9+y2*y7*y8)*y1**3+(y3**3*y8*y9+y4**3*y7*y9+ &
          y2**3*y7*y8)*y1+(-y4*y8*y9-y3*y7*y9)*y2**3+(-y3**3*y7*y9-y4**3*y8*y9)*y2- &
          y3**3*y4*y7*y8-y3*y4**3*y7*y8
      s1 = (((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark-y8/2._ark)*y6)*y2+(- &
          sqrt(3._ark)*y5*y8/2._ark+(y9+y8/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y7/2._ark+(-y9-y7/ &
          2._ark)*y6)*y4)*y1**3+(((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/ &
          2._ark)*y6)*y2**3+(sqrt(3._ark)*y5*y8/2._ark+(-y8/2._ark-y9)*y6)*y3**3+ &
          (sqrt(3._ark)*y5*y7/2._ark+(y9+y7/2._ark)*y6)*y4**3)*y1+((sqrt(3._ark)*y5*y7/2._ark+(y7/ &
          2._ark-y9)*y6)*y3+(sqrt(3._ark)*y5*y8/2._ark+(-y8/2._ark+y9)*y6)*y4)*y2**3
      dF(206) = s1+((-sqrt(3._ark)*y5*y7/2._ark+(-y7/2._ark+y9)*y6)*y3**3+(- &
          sqrt(3._ark)*y5*y8/2._ark+(-y9+y8/2._ark)*y6)*y4**3)*y2+((sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark+y8/2._ark)*y6)*y4*y3**3+((-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark-y7/2._ark)*y6)*y4**3*y3
      s1 = (((y7/2._ark+y8/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6)*y2+((- &
          y9+y8/2._ark)*y5+sqrt(3._ark)*y6*y8/2._ark)*y3+((y7/2._ark-y9)*y5-sqrt(3._ark)*y6*y7/ &
          2._ark)*y4)*y1**3+(((-y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y2**3+((-y8/2._ark+y9)*y5-sqrt(3._ark)*y6*y8/2._ark)*y3**3+((-y7/2._ark+ &
          y9)*y5+sqrt(3._ark)*y6*y7/2._ark)*y4**3)*y1+(((-y9-y7/2._ark)*y5+sqrt(3._ark)*y6*y7/ &
          2._ark)*y3+((-y8/2._ark-y9)*y5-sqrt(3._ark)*y6*y8/2._ark)*y4)*y2**3
      dF(207) = s1+(((y9+y7/2._ark)*y5-sqrt(3._ark)*y6*y7/2._ark)*y3**3+((y9+y8/2._ark)*y5+ &
          sqrt(3._ark)*y6*y8/2._ark)*y4**3)*y2+((y7/2._ark-y8/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y6)*y4*y3**3+((y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6)*y4**3*y3
      dF(208) = ((-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y2+(-y5*y6+ &
          sqrt(3._ark)*y6**2/3._ark)*y3+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4)*y1**3+((- &
          sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y2**3+(-y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y3**3+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4**3)*y1+((y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y3+(-y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4)*y2**3+((y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y3**3+(-y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4**3)*y2+(-sqrt(3._ark)*y6**2/6._ark+ &
          sqrt(3._ark)*y5**2/2._ark)*y4*y3**3+(-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/ &
          2._ark)*y4**3*y3
      dF(209) = ((y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y1**3+((y6**2+ &
          y5**2)*y2**3+(y6**2+y5**2)*y3**3+(y6**2+y5**2)*y4**3)*y1+((y6**2+y5**2)*y3+ &
          (y6**2+y5**2)*y4)*y2**3+((y6**2+y5**2)*y3**3+(y6**2+y5**2)*y4**3)*y2+(y6**2+ &
          y5**2)*y4*y3**3+(y6**2+y5**2)*y4**3*y3
      dF(210) = (y3*y7+y4*y8+y2*y9)*y1**4+(y2**4*y9+y4**4*y8+y3**4*y7)*y1+(-y3*y8- &
          y4*y7)*y2**4+(-y4**4*y7-y3**4*y8)*y2-y3**4*y4*y9-y3*y4**4*y9
      dF(211) = ((-y8-y7)*y2+(-y9-y8)*y3+(-y9-y7)*y4)*y1**4+((y7+y8)*y2**4+(y9+ &
          y8)*y3**4+(y7+y9)*y4**4)*y1+((y7-y9)*y3+(y8-y9)*y4)*y2**4+((y9-y7)*y3**4+(-y8+ &
          y9)*y4**4)*y2+(y8-y7)*y4*y3**4+(y7-y8)*y4**4*y3
      dF(212) = (-2._ark/3._ark*sqrt(3._ark)*y2*y5+(-y6+sqrt(3._ark)*y5/3._ark)*y3+(y6+ &
          sqrt(3._ark)*y5/3._ark)*y4)*y1**4+(-2._ark/3._ark*sqrt(3._ark)*y2**4*y5+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y3**4+(y6+sqrt(3._ark)*y5/3._ark)*y4**4)*y1+((y6+sqrt(3._ark)*y5/ &
          3._ark)*y3+(-y6+sqrt(3._ark)*y5/3._ark)*y4)*y2**4+((y6+sqrt(3._ark)*y5/3._ark)*y3**4+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**4)*y2-2._ark/3._ark*sqrt(3._ark)*y3*y4**4*y5-2._ark/ &
          3._ark*sqrt(3._ark)*y3**4*y4*y5
      dF(213) = ((y8**4+y7**4)*y2+(y8**4+y9**4)*y3+(y9**4+y7**4)*y4)*y1+((y9**4+ &
          y7**4)*y3+(y8**4+y9**4)*y4)*y2+(y8**4+y7**4)*y4*y3
      dF(214) = ((y7**3*y8+y7*y8**3)*y2+(y8**3*y9+y8*y9**3)*y3+(y7*y9**3+ &
          y7**3*y9)*y4)*y1+((-y7*y9**3-y7**3*y9)*y3+(-y8**3*y9-y8*y9**3)*y4)*y2+(- &
          y7**3*y8-y7*y8**3)*y4*y3
      dF(215) = (y4*y7**2*y9**2+y2*y7**2*y8**2+y3*y8**2*y9**2)*y1+(y3*y7**2*y9**2+ &
          y4*y8**2*y9**2)*y2+y3*y4*y7**2*y8**2
      dF(216) = (y3*y7**2*y8*y9+y4*y7*y8**2*y9+y2*y7*y8*y9**2)*y1+(-y4*y7**2*y8*y9- &
          y3*y7*y8**2*y9)*y2-y3*y4*y7*y8*y9**2
      dF(217) = (y4*y8**4+y3*y7**4+y2*y9**4)*y1+(y3*y8**4+y4*y7**4)*y2+y3*y4*y9**4
      dF(218) = (8._ark/3._ark*sqrt(3._ark)*y2*y5*y6**2*y9+(-sqrt(3._ark)*y5**3*y7+ &
          y5**2*y6*y7+5._ark/3._ark*sqrt(3._ark)*y5*y6**2*y7+y6**3*y7)*y3+(-y5**2*y6*y8+5._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2*y8-sqrt(3._ark)*y5**3*y8-y6**3*y8)*y4)*y1+((y6**3*y8+ &
          y5**2*y6*y8-5._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8+sqrt(3._ark)*y5**3*y8)*y3+(- &
          y5**2*y6*y7-y6**3*y7+sqrt(3._ark)*y5**3*y7-5._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2*y7)*y4)*y2-8._ark/3._ark*sqrt(3._ark)*y3*y4*y5*y6**2*y9
      dF(219) = ((y5**2*y9**2+y6**2*y9**2)*y2+(y6**2*y7**2+y5**2*y7**2)*y3+ &
          (y6**2*y8**2+y5**2*y8**2)*y4)*y1+((y6**2*y8**2+y5**2*y8**2)*y3+(y6**2*y7**2+ &
          y5**2*y7**2)*y4)*y2+(y5**2*y9**2+y6**2*y9**2)*y4*y3
      dF(220) = ((4._ark/3._ark*y6**4+4._ark*y5**2*y6**2)*y2+(4._ark/3._ark*sqrt(3._ark)*y5*y6**3+ &
          y5**2*y6**2+5._ark/6._ark*y6**4+3._ark/2._ark*y5**4)*y3+(3._ark/2._ark*y5**4-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2+5._ark/6._ark*y6**4)*y4)*y1+((3._ark/2._ark*y5**4- &
          4._ark/3._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2+5._ark/6._ark*y6**4)*y3+(4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2+5._ark/6._ark*y6**4+3._ark/2._ark*y5**4)*y4)*y2+ &
          (4._ark/3._ark*y6**4+4._ark*y5**2*y6**2)*y4*y3
      dF(221) = ((y4*y7*y8*y9+y3*y7*y8*y9)*y2+y3*y4*y7*y8*y9)*y1+y2*y3*y4*y7*y8*y9
      dF(222) = (((((-y9/2._ark-y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark- &
          sqrt(3._ark)*y7*y9/2._ark)*y6)*y3+(((y9/2._ark-y8)*y7-y8*y9/2._ark)*y5+ &
          (sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4)*y2+(((y8+y9/2._ark)*y7+ &
          y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4*y3)*y1+ &
          (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/ &
          2._ark)*y6)*y4*y3*y2
      dF(223) = ((((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y3+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4)*y2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4*y3)*y1+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4*y3*y2
      dF(224) = ((((y9+y7-y8)*y5**2+(y9+y7-y8)*y6**2)*y3+((y8+y9-y7)*y5**2+(y8+y9- &
          y7)*y6**2)*y4)*y2+((y7-y9+y8)*y5**2+(y7-y9+y8)*y6**2)*y4*y3)*y1+((-y7-y8- &
          y9)*y5**2+(-y7-y8-y9)*y6**2)*y4*y3*y2
      dF(225) = (((y5**3-3._ark*y5*y6**2)*y3+(y5**3-3._ark*y5*y6**2)*y4)*y2+(y5**3- &
          3._ark*y5*y6**2)*y4*y3)*y1+(y5**3-3._ark*y5*y6**2)*y4*y3*y2
      dF(226) = (((-sqrt(3._ark)*y6*y8-y5*y8)*y3+(sqrt(3._ark)*y6*y7-y5*y7)*y4)*y2+ &
          2._ark*y3*y4*y5*y9)*y1**2+(((-sqrt(3._ark)*y6*y7+y5*y7)*y3+(sqrt(3._ark)*y6*y8+ &
          y5*y8)*y4)*y2**2+(-2._ark*y3**2*y5*y9-2._ark*y4**2*y5*y9)*y2+(sqrt(3._ark)*y6*y8+ &
          y5*y8)*y4*y3**2+(-sqrt(3._ark)*y6*y7+y5*y7)*y4**2*y3)*y1+2._ark*y2**2*y3*y4*y5*y9+ &
          ((sqrt(3._ark)*y6*y7-y5*y7)*y4*y3**2+(-sqrt(3._ark)*y6*y8-y5*y8)*y4**2*y3)*y2
      dF(227) = (((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3)*y1**2+ &
          (((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2**2+((y6**2+y5**2)*y3**2+(y6**2+ &
          y5**2)*y4**2)*y2+(y6**2+y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3)*y1+(y6**2+ &
          y5**2)*y4*y3*y2**2+((y6**2+y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3)*y2
      dF(228) = ((sqrt(3._ark)*y5**3-sqrt(3._ark)*y5*y6**2/3._ark)*y2+(y6**3-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6)*y3+(-y6**3-y5**2*y6-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2)*y4)*y1**2+((sqrt(3._ark)*y5**3-sqrt(3._ark)*y5*y6**2/ &
          3._ark)*y2**2+(y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6)*y3**2+(-y6**3- &
          y5**2*y6-4._ark/3._ark*sqrt(3._ark)*y5*y6**2)*y4**2)*y1+((-y6**3-y5**2*y6-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2)*y3+(y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2+ &
          y5**2*y6)*y4)*y2**2+((-y6**3-y5**2*y6-4._ark/3._ark*sqrt(3._ark)*y5*y6**2)*y3**2+ &
          (y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6)*y4**2)*y2+(sqrt(3._ark)*y5**3- &
          sqrt(3._ark)*y5*y6**2/3._ark)*y4*y3**2+(sqrt(3._ark)*y5**3-sqrt(3._ark)*y5*y6**2/ &
          3._ark)*y4**2*y3
      dF(229) = ((-y8**2*y9+y7**2*y9)*y6*y2+((-sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y9**2/ &
          2._ark)*y7*y5+(y9**2/2._ark-y8**2/2._ark)*y7*y6)*y3+((-sqrt(3._ark)*y7**2*y8/2._ark+ &
          sqrt(3._ark)*y8*y9**2/2._ark)*y5+(-y8*y9**2/2._ark+y7**2*y8/2._ark)*y6)*y4)*y1+(((- &
          sqrt(3._ark)*y8*y9**2/2._ark+sqrt(3._ark)*y7**2*y8/2._ark)*y5+(-y7**2*y8/2._ark+y8*y9**2/ &
          2._ark)*y6)*y3+((sqrt(3._ark)*y8**2/2._ark-sqrt(3._ark)*y9**2/2._ark)*y7*y5+(y8**2/2._ark- &
          y9**2/2._ark)*y7*y6)*y4)*y2+(-y7**2*y9+y8**2*y9)*y6*y4*y3
      dF(230) = (y2*y5*y9**3+(sqrt(3._ark)*y6*y7**3/2._ark-y5*y7**3/2._ark)*y3+(- &
          sqrt(3._ark)*y6*y8**3/2._ark-y5*y8**3/2._ark)*y4)*y1+((y5*y8**3/2._ark+ &
          sqrt(3._ark)*y6*y8**3/2._ark)*y3+(-sqrt(3._ark)*y6*y7**3/2._ark+y5*y7**3/2._ark)*y4)*y2- &
          y3*y4*y5*y9**3
      dF(231) = (y2*y5*y7*y8*y9+(sqrt(3._ark)*y6*y7*y8*y9/2._ark-y5*y7*y8*y9/2._ark)*y3+(- &
          y5*y7*y8*y9/2._ark-sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y4)*y1+((-y5*y7*y8*y9/2._ark- &
          sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y3+(sqrt(3._ark)*y6*y7*y8*y9/2._ark-y5*y7*y8*y9/ &
          2._ark)*y4)*y2+y3*y4*y5*y7*y8*y9
      dF(232) = ((y7**2*y9+y8**2*y9)*y5*y2+((-y8**2/2._ark-y9**2/2._ark)*y7*y5+ &
          (sqrt(3._ark)*y9**2/2._ark+sqrt(3._ark)*y8**2/2._ark)*y7*y6)*y3+((-y8*y9**2/2._ark- &
          y7**2*y8/2._ark)*y5+(-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y8/ &
          2._ark)*y6)*y4)*y1+(((y8*y9**2/2._ark+y7**2*y8/2._ark)*y5+(sqrt(3._ark)*y8*y9**2/2._ark+ &
          sqrt(3._ark)*y7**2*y8/2._ark)*y6)*y3+((y8**2/2._ark+y9**2/2._ark)*y7*y5+(- &
          sqrt(3._ark)*y9**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y7*y6)*y4)*y2+(-y8**2*y9- &
          y7**2*y9)*y5*y4*y3
      dF(233) = (((-sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y5**2+(y7**2- &
          y8**2)*y6*y5+(-sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y9**2+(-sqrt(3._ark)*y9**2/6._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y8**2)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9**2/2._ark-y5*y6*y9**2+(- &
          sqrt(3._ark)*y9**2/6._ark-2._ark/3._ark*sqrt(3._ark)*y7**2)*y6**2)*y4)*y1+((- &
          sqrt(3._ark)*y5**2*y9**2/2._ark-y5*y6*y9**2+(-sqrt(3._ark)*y9**2/6._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y7**2)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y9**2+(- &
          sqrt(3._ark)*y9**2/6._ark-2._ark/3._ark*sqrt(3._ark)*y8**2)*y6**2)*y4)*y2+((- &
          sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y5**2+(y7**2-y8**2)*y6*y5+(- &
          sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y4*y3
      dF(234) = ((-3._ark*y5*y6**2*y9+y5**3*y9)*y2+(y5**3*y7-3._ark*y5*y6**2*y7)*y3+ &
          (y5**3*y8-3._ark*y5*y6**2*y8)*y4)*y1+((-y5**3*y8+3._ark*y5*y6**2*y8)*y3+(-y5**3*y7+ &
          3._ark*y5*y6**2*y7)*y4)*y2+(-y5**3*y9+3._ark*y5*y6**2*y9)*y4*y3
      dF(235) = ((-5._ark/3._ark*y6**4-6._ark*y5**2*y6**2+y5**4)*y2+(-8._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**3-2._ark/3._ark*y6**4-2._ark*y5**4)*y3+(-2._ark*y5**4-2._ark/ &
          3._ark*y6**4+8._ark/3._ark*sqrt(3._ark)*y5*y6**3)*y4)*y1+((-2._ark*y5**4-2._ark/3._ark*y6**4+ &
          8._ark/3._ark*sqrt(3._ark)*y5*y6**3)*y3+(-8._ark/3._ark*sqrt(3._ark)*y5*y6**3-2._ark/ &
          3._ark*y6**4-2._ark*y5**4)*y4)*y2+(-5._ark/3._ark*y6**4-6._ark*y5**2*y6**2+y5**4)*y4*y3
      dF(236) = (((y7**3+y9**3-y8**3)*y3+(y8**3+y9**3-y7**3)*y4)*y2+(-y9**3+y8**3+ &
          y7**3)*y4*y3)*y1+(-y8**3-y7**3-y9**3)*y4*y3*y2
      dF(237) = ((((-y8+y9)*y7**2+(y9**2+y8**2)*y7-y8*y9**2+y8**2*y9)*y3+((y9+ &
          y8)*y7**2+(-y9**2-y8**2)*y7+y8**2*y9+y8*y9**2)*y4)*y2+((y8-y9)*y7**2+(y9**2+ &
          y8**2)*y7-y8**2*y9+y8*y9**2)*y4*y3)*y1+((-y9-y8)*y7**2+(-y9**2-y8**2)*y7- &
          y8*y9**2-y8**2*y9)*y4*y3*y2
      dF(238) = (((-sqrt(3._ark)*y5**2*y9/2._ark+(y7+y8)*y6*y5+(sqrt(3._ark)*y9/6._ark+ &
          sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9/2._ark+(- &
          y8-y7)*y6*y5+(sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/ &
          3._ark)*y6**2)*y4)*y2+(sqrt(3._ark)*y5**2*y9/2._ark+(y7-y8)*y6*y5+(-sqrt(3._ark)*y7/ &
          3._ark-sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y9/6._ark)*y6**2)*y4*y3)*y1+ &
          (sqrt(3._ark)*y5**2*y9/2._ark+(y8-y7)*y6*y5+(sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark- &
          sqrt(3._ark)*y9/6._ark)*y6**2)*y4*y3*y2
      dF(239) = (y8**2+y7**2+y9**2)*y4*y3*y2*y1
      dF(240) = (y6**2+y5**2)*y4*y3*y2*y1
      dF(241) = (y9**4+y8**4+y7**4)*y1**2+(y9**4+y8**4+y7**4)*y2**2+(y9**4+y8**4+ &
          y7**4)*y3**2+(y9**4+y8**4+y7**4)*y4**2
      dF(242) = (y7**2*y8*y9+(y8*y9**2+y8**2*y9)*y7)*y1**2+(-y7**2*y8*y9+(y8*y9**2- &
          y8**2*y9)*y7)*y2**2+(y7**2*y8*y9+(-y8**2*y9-y8*y9**2)*y7)*y3**2+(-y7**2*y8*y9+(- &
          y8*y9**2+y8**2*y9)*y7)*y4**2
      dF(243) = ((y9**2+y8**2)*y7**2+y8**2*y9**2)*y1**2+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y2**2+((y9**2+y8**2)*y7**2+y8**2*y9**2)*y3**2+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y4**2
      dF(244) = ((y9+y8)*y7**3+(y9**3+y8**3)*y7+y8**3*y9+y8*y9**3)*y1**2+((y8- &
          y9)*y7**3+(y8**3-y9**3)*y7-y8*y9**3-y8**3*y9)*y2**2+((-y9-y8)*y7**3+(-y9**3- &
          y8**3)*y7+y8**3*y9+y8*y9**3)*y3**2+((-y8+y9)*y7**3+(y9**3-y8**3)*y7-y8*y9**3- &
          y8**3*y9)*y4**2
      s1 = ((sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**2/2._ark)*y5+((y8+y9/2._ark)*y7**2+(-y9**2/2._ark-y8**2)*y7- &
          y8**2*y9/2._ark+y8*y9**2/2._ark)*y6)*y1**2+((-sqrt(3._ark)*y7*y9**2/2._ark- &
          sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark)*y5+ &
          ((y9/2._ark-y8)*y7**2+(y9**2/2._ark+y8**2)*y7-y8*y9**2/2._ark-y8**2*y9/2._ark)*y6)*y2**2
      dF(245) = s1+((sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y8*y9**2/2._ark)*y5+((-y9/2._ark-y8)*y7**2+(- &
          y9**2/2._ark-y8**2)*y7+y8**2*y9/2._ark-y8*y9**2/2._ark)*y6)*y3**2+ &
          ((sqrt(3._ark)*y8**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y7**2*y9/2._ark)*y5+((y8-y9/2._ark)*y7**2+(y9**2/2._ark+y8**2)*y7+y8**2*y9/ &
          2._ark+y8*y9**2/2._ark)*y6)*y4**2
      dF(246) = ((-sqrt(3._ark)*y7**3/3._ark-sqrt(3._ark)*y8**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3)*y5+(-y8**3+y7**3)*y6)*y1**2+((sqrt(3._ark)*y7**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3+y8**3)*y6)*y2**2+((- &
          2._ark/3._ark*sqrt(3._ark)*y9**3-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/3._ark)*y5+ &
          (y8**3+y7**3)*y6)*y3**2+((-2._ark/3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y7**3/3._ark- &
          sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3-y8**3)*y6)*y4**2
      s1 = (((y8-y9/2._ark)*y7**2+(-y9**2/2._ark+y8**2)*y7-y8*y9**2/2._ark-y8**2*y9/ &
          2._ark)*y5+(-sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark+ &
          sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7*y9**2/2._ark)*y6)*y1**2+(((-y9/2._ark- &
          y8)*y7**2+(y9**2/2._ark-y8**2)*y7+y8*y9**2/2._ark-y8**2*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8**2*y9/2._ark- &
          sqrt(3._ark)*y8*y9**2/2._ark)*y6)*y2**2
      dF(247) = s1+(((y9/2._ark-y8)*y7**2+(-y9**2/2._ark+y8**2)*y7+y8**2*y9/2._ark+y8*y9**2/ &
          2._ark)*y5+(sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y8*y9**2/ &
          2._ark-sqrt(3._ark)*y7*y9**2/2._ark)*y6)*y3**2+(((y8+y9/2._ark)*y7**2+(y9**2/2._ark- &
          y8**2)*y7+y8**2*y9/2._ark-y8*y9**2/2._ark)*y5+(sqrt(3._ark)*y7*y9**2/2._ark- &
          sqrt(3._ark)*y8**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark+sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y6)*y4**2
      dF(248) = ((-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5**2+(-y8*y9+ &
          y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark-2._ark/3._ark*sqrt(3._ark)*y8)*y7-sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y1**2+((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(y8*y9- &
          y7*y9)*y6*y5+((-2._ark/3._ark*sqrt(3._ark)*y8+sqrt(3._ark)*y9/6._ark)*y7+sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y2**2+((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(- &
          y8*y9-y7*y9)*y6*y5+((sqrt(3._ark)*y9/6._ark+2._ark/3._ark*sqrt(3._ark)*y8)*y7- &
          sqrt(3._ark)*y8*y9/6._ark)*y6**2)*y3**2+((-sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/ &
          2._ark)*y5**2+(y8*y9+y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y8)*y7+sqrt(3._ark)*y8*y9/6._ark)*y6**2)*y4**2
      dF(249) = (-sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/ &
          6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y1**2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y2**2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y3**2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y4**2
      dF(250) = ((-y9/3._ark-y8/3._ark-y7/3._ark)*y5**3+(y9+y7+y8)*y6**2*y5)*y1**2+((y8/ &
          3._ark-y9/3._ark+y7/3._ark)*y5**3+(-y7+y9-y8)*y6**2*y5)*y2**2+((y9/3._ark+y8/3._ark-y7/ &
          3._ark)*y5**3+(-y8-y9+y7)*y6**2*y5)*y3**2+((y7/3._ark+y9/3._ark-y8/3._ark)*y5**3+(-y7- &
          y9+y8)*y6**2*y5)*y4**2
      dF(251) = ((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y1**2+((y8**2+ &
          y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y2**2+((y8**2+y7**2+y9**2)*y5**2+ &
          (y8**2+y7**2+y9**2)*y6**2)*y3**2+((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+ &
          y9**2)*y6**2)*y4**2
      dF(252) = ((-4._ark/9._ark*sqrt(3._ark)*y8-4._ark/9._ark*sqrt(3._ark)*y7+5._ark/ &
          9._ark*sqrt(3._ark)*y9)*y5**3+(y7-y8)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9+(y7- &
          y8)*y6**3)*y1**2+((4._ark/9._ark*sqrt(3._ark)*y7+5._ark/9._ark*sqrt(3._ark)*y9+4._ark/ &
          9._ark*sqrt(3._ark)*y8)*y5**3+(y8-y7)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9+(y8- &
          y7)*y6**3)*y2**2+((4._ark/9._ark*sqrt(3._ark)*y8-5._ark/9._ark*sqrt(3._ark)*y9-4._ark/ &
          9._ark*sqrt(3._ark)*y7)*y5**3+(y7+y8)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y9+(y7+ &
          y8)*y6**3)*y3**2+((-5._ark/9._ark*sqrt(3._ark)*y9+4._ark/9._ark*sqrt(3._ark)*y7-4._ark/ &
          9._ark*sqrt(3._ark)*y8)*y5**3+(-y8-y7)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y9+(-y8- &
          y7)*y6**3)*y4**2
      dF(253) = ((sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y2+(- &
          sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3+(-sqrt(3._ark)*y6**2*y8/3._ark- &
          y5*y6*y8)*y4)*y1**2+((sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y2**2+ &
          (-sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3**2+(-sqrt(3._ark)*y6**2*y8/3._ark- &
          y5*y6*y8)*y4**2)*y1+((sqrt(3._ark)*y6**2*y8/3._ark+y5*y6*y8)*y3+(-y5*y6*y7+ &
          sqrt(3._ark)*y6**2*y7/3._ark)*y4)*y2**2+((sqrt(3._ark)*y6**2*y8/3._ark+y5*y6*y8)*y3**2+ &
          (-y5*y6*y7+sqrt(3._ark)*y6**2*y7/3._ark)*y4**2)*y2+(-sqrt(3._ark)*y6**2*y9/6._ark+ &
          sqrt(3._ark)*y5**2*y9/2._ark)*y4*y3**2+(-sqrt(3._ark)*y6**2*y9/6._ark+ &
          sqrt(3._ark)*y5**2*y9/2._ark)*y4**2*y3
      dF(254) = (2._ark/3._ark*sqrt(3._ark)*y2**2*y6**2+(y5*y6+sqrt(3._ark)*y5**2/2._ark+ &
          sqrt(3._ark)*y6**2/6._ark)*y3**2+(-y5*y6+sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/ &
          2._ark)*y4**2)*y1**2+((-y5*y6+sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y3**2+ &
          (y5*y6+sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/6._ark)*y4**2)*y2**2+2._ark/ &
          3._ark*sqrt(3._ark)*y3**2*y4**2*y6**2
      dF(255) = (y2*y5*y7*y8+(sqrt(3._ark)*y6*y8*y9/2._ark-y5*y8*y9/2._ark)*y3+(- &
          sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y4)*y1**2+(y2**2*y5*y7*y8+ &
          (sqrt(3._ark)*y6*y8*y9/2._ark-y5*y8*y9/2._ark)*y3**2+(-sqrt(3._ark)*y6*y7*y9/2._ark- &
          y5*y7*y9/2._ark)*y4**2)*y1+((sqrt(3._ark)*y6*y7*y9/2._ark+y5*y7*y9/2._ark)*y3+(- &
          sqrt(3._ark)*y6*y8*y9/2._ark+y5*y8*y9/2._ark)*y4)*y2**2+((sqrt(3._ark)*y6*y7*y9/2._ark+ &
          y5*y7*y9/2._ark)*y3**2+(-sqrt(3._ark)*y6*y8*y9/2._ark+y5*y8*y9/2._ark)*y4**2)*y2- &
          y3**2*y4*y5*y7*y8-y3*y4**2*y5*y7*y8
      dF(256) = ((y3*y7*y9+y4*y8*y9)*y2+y3*y4*y7*y8)*y1**2+((-y4*y7*y9- &
          y3*y8*y9)*y2**2+(-y3**2*y7*y8-y4**2*y7*y8)*y2-y3*y4**2*y8*y9-y3**2*y4*y7*y9)*y1+ &
          y2**2*y3*y4*y7*y8+(y3**2*y4*y8*y9+y3*y4**2*y7*y9)*y2
      dF(257) = (((y7**2+y9**2)*y3+(y9**2+y8**2)*y4)*y2+(y8**2+y7**2)*y4*y3)*y1**2+ &
          (((y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y2**2+((y8**2+y7**2)*y3**2+(y8**2+ &
          y7**2)*y4**2)*y2+(y7**2+y9**2)*y4*y3**2+(y9**2+y8**2)*y4**2*y3)*y1+(y8**2+ &
          y7**2)*y4*y3*y2**2+((y9**2+y8**2)*y4*y3**2+(y7**2+y9**2)*y4**2*y3)*y2
      s1 = (((sqrt(3._ark)*y5*y9/2._ark+(y7+y9/2._ark)*y6)*y3+(sqrt(3._ark)*y5*y9/2._ark+(-y9/ &
          2._ark-y8)*y6)*y4)*y2+((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark-y8/ &
          2._ark)*y6)*y4*y3)*y1**2+(((sqrt(3._ark)*y5*y9/2._ark+(y8-y9/2._ark)*y6)*y3+ &
          (sqrt(3._ark)*y5*y9/2._ark+(-y7+y9/2._ark)*y6)*y4)*y2**2+(((-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark+y8/2._ark)*y6)*y3**2+((sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark-y7/2._ark)*y6)*y4**2)*y2+(-sqrt(3._ark)*y5*y9/2._ark+ &
          (-y9/2._ark+y7)*y6)*y4*y3**2+(-sqrt(3._ark)*y5*y9/2._ark+(y9/2._ark-y8)*y6)*y4**2*y3)*y1
      dF(258) = s1+((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/ &
          2._ark)*y6)*y4*y3*y2**2+((-sqrt(3._ark)*y5*y9/2._ark+(y8+y9/2._ark)*y6)*y4*y3**2+(- &
          sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark-y7)*y6)*y4**2*y3)*y2
      s1 = ((((-y9/2._ark+y7)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3+((y8-y9/2._ark)*y5- &
          sqrt(3._ark)*y6*y9/2._ark)*y4)*y2+((-y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6)*y4*y3)*y1**2+((((-y9/2._ark-y8)*y5-sqrt(3._ark)*y6*y9/ &
          2._ark)*y3+((-y9/2._ark-y7)*y5+sqrt(3._ark)*y6*y9/2._ark)*y4)*y2**2+(((y8/2._ark-y7/ &
          2._ark)*y5+(-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y6)*y3**2+((y7/2._ark-y8/ &
          2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6)*y4**2)*y2+((y7+y9/2._ark)*y5- &
          sqrt(3._ark)*y6*y9/2._ark)*y4*y3**2+((y8+y9/2._ark)*y5+sqrt(3._ark)*y6*y9/ &
          2._ark)*y4**2*y3)*y1
      dF(259) = s1+((y7/2._ark+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y4*y3*y2**2+(((y9/2._ark-y8)*y5+sqrt(3._ark)*y6*y9/2._ark)*y4*y3**2+((-y7+ &
          y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y4**2*y3)*y2
      dF(260) = (((y5*y6+sqrt(3._ark)*y6**2/3._ark)*y3+(-y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y4)*y2+(-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y4*y3)*y1**2+(((- &
          y5*y6+sqrt(3._ark)*y6**2/3._ark)*y3+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4)*y2**2+((- &
          sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y3**2+(-sqrt(3._ark)*y6**2/6._ark+ &
          sqrt(3._ark)*y5**2/2._ark)*y4**2)*y2+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4*y3**2+(-y5*y6+ &
          sqrt(3._ark)*y6**2/3._ark)*y4**2*y3)*y1+(-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/ &
          2._ark)*y4*y3*y2**2+((-y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4*y3**2+(y5*y6+ &
          sqrt(3._ark)*y6**2/3._ark)*y4**2*y3)*y2
      dF(261) = (y9+y7+y8)*y4*y3*y2*y1**2+((-y7+y9-y8)*y4*y3*y2**2+((-y8-y9+ &
          y7)*y4*y3**2+(-y7-y9+y8)*y4**2*y3)*y2)*y1
      dF(262) = (y4**2*y7*y9+y2**2*y7*y8+y3**2*y8*y9)*y1**2+(-y4**2*y8*y9- &
          y3**2*y7*y9)*y2**2-y3**2*y4**2*y7*y8
      dF(263) = (y2**2*y5*y9+(-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y3**2+(-y5*y8/2._ark- &
          sqrt(3._ark)*y6*y8/2._ark)*y4**2)*y1**2+((y5*y8/2._ark+sqrt(3._ark)*y6*y8/2._ark)*y3**2+ &
          (y5*y7/2._ark-sqrt(3._ark)*y6*y7/2._ark)*y4**2)*y2**2-y3**2*y4**2*y5*y9
      dF(264) = ((y6**2+y5**2)*y2**2+(y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y1**2+ &
          ((y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y2**2+(y6**2+y5**2)*y4**2*y3**2
      dF(265) = ((y3*y9+y4*y9)*y2**2+(y4**2*y8+y3**2*y7)*y2+y3*y4**2*y8+ &
          y3**2*y4*y7)*y1**2+((-y3**2*y8-y4**2*y7)*y2**2-y3**2*y4**2*y9)*y1+(-y3*y4**2*y7- &
          y3**2*y4*y8)*y2**2-y2*y3**2*y4**2*y9
      dF(266) = (((y7-y8)*y3+(y8-y7)*y4)*y2**2+((-y8+y9)*y3**2+(y9-y7)*y4**2)*y2+(y8- &
          y9)*y4*y3**2+(y7-y9)*y4**2*y3)*y1**2+(((y7+y9)*y3**2+(y9+y8)*y4**2)*y2**2+(y7+ &
          y8)*y4**2*y3**2)*y1+((-y9-y7)*y4*y3**2+(-y9-y8)*y4**2*y3)*y2**2+(-y8- &
          y7)*y4**2*y3**2*y2
      dF(267) = ((y4*y5+y3*y5)*y2**2+((-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4**2)*y2+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y4*y3**2+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4**2*y3)*y1**2+(((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/ &
          2._ark+sqrt(3._ark)*y6/2._ark)*y4**2)*y2**2+y3**2*y4**2*y5)*y1+((-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4*y3**2+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y4**2*y3)*y2**2+ &
          y2*y3**2*y4**2*y5
      dF(268) = (y2**2*y3*y4+(y3*y4**2+y3**2*y4)*y2)*y1**2+((y3*y4**2+y3**2*y4)*y2**2+ &
          y2*y3**2*y4**2)*y1
      dF(269) = ((y3**2+y4**2)*y2**2+y3**2*y4**2)*y1**2+y2**2*y3**2*y4**2
      dF(270) = (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y6)*y1**3+(((y8+y9/2._ark)*y7+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y2**3+(((y9/2._ark-y8)*y7-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3**3+(((-y9/2._ark- &
          y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4**3
      dF(271) = (-sqrt(3._ark)*y5**2*y9/2._ark+(y7-y8)*y6*y5+(sqrt(3._ark)*y9/6._ark- &
          sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y6**2)*y1**3+(-sqrt(3._ark)*y5**2*y9/2._ark+ &
          (y8-y7)*y6*y5+(sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/ &
          3._ark)*y6**2)*y2**3+(sqrt(3._ark)*y5**2*y9/2._ark+(y7+y8)*y6*y5+(sqrt(3._ark)*y8/3._ark- &
          sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y9/6._ark)*y6**2)*y3**3+(sqrt(3._ark)*y5**2*y9/2._ark+(- &
          y8-y7)*y6*y5+(-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/ &
          3._ark)*y6**2)*y4**3
      dF(272) = ((y9+y7+y8)*y5**2+(y9+y7+y8)*y6**2)*y1**3+((-y7+y9-y8)*y5**2+(-y7+y9- &
          y8)*y6**2)*y2**3+((-y8-y9+y7)*y5**2+(-y8-y9+y7)*y6**2)*y3**3+((-y7-y9+y8)*y5**2+ &
          (-y7-y9+y8)*y6**2)*y4**3
      dF(273) = (y5**3-3._ark*y5*y6**2)*y1**3+(y5**3-3._ark*y5*y6**2)*y2**3+(y5**3- &
          3._ark*y5*y6**2)*y3**3+(y5**3-3._ark*y5*y6**2)*y4**3
      dF(274) = (y4**3+y3**3+y2**3)*y1**3+(y4**3+y3**3)*y2**3+y3**3*y4**3
      dF(275) = (2._ark/3._ark*sqrt(3._ark)*y2*y5*y9+(-sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3+(- &
          sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4)*y1**3+(2._ark/3._ark*sqrt(3._ark)*y2**3*y5*y9+(- &
          sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3**3+(-sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4**3)*y1+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4)*y2**3+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3**3+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4**3)*y2- &
          2._ark/3._ark*sqrt(3._ark)*y3*y4**3*y5*y9-2._ark/3._ark*sqrt(3._ark)*y3**3*y4*y5*y9
      dF(276) = ((y3*y8+y4*y7)*y2+y3*y4*y9)*y1**3+((-y3*y7-y4*y8)*y2**3+(-y3**3*y9- &
          y4**3*y9)*y2-y3**3*y4*y8-y3*y4**3*y7)*y1+y2**3*y3*y4*y9+(y3*y4**3*y8+ &
          y3**3*y4*y7)*y2
      dF(277) = (((y7+y9)*y3+(y9+y8)*y4)*y2+(y7+y8)*y4*y3)*y1**3+(((-y8+y9)*y3+(y9- &
          y7)*y4)*y2**3+((y7-y8)*y3**3+(y8-y7)*y4**3)*y2+(y7-y9)*y4*y3**3+(y8- &
          y9)*y4**3*y3)*y1+(-y8-y7)*y4*y3*y2**3+((-y9-y8)*y4*y3**3+(-y9-y7)*y4**3*y3)*y2
      dF(278) = (((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark+sqrt(3._ark)*y6/ &
          2._ark)*y4)*y2+y3*y4*y5)*y1**3+(((-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4)*y2**3+(y3**3*y5+y4**3*y5)*y2+(-y5/2._ark-sqrt(3._ark)*y6/ &
          2._ark)*y4*y3**3+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y4**3*y3)*y1+y2**3*y3*y4*y5+((-y5/ &
          2._ark+sqrt(3._ark)*y6/2._ark)*y4*y3**3+(-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y4**3*y3)*y2
      dF(279) = ((y7+y8)*y2**2+(y9+y8)*y3**2+(y7+y9)*y4**2)*y1**3+((-y8-y7)*y2**3+(- &
          y9-y8)*y3**3+(-y9-y7)*y4**3)*y1**2+((y9-y7)*y3**2+(-y8+y9)*y4**2)*y2**3+((y7- &
          y9)*y3**3+(y8-y9)*y4**3)*y2**2+(y7-y8)*y4**2*y3**3+(y8-y7)*y4**3*y3**2
      dF(280) = (y4**2*y8+y3**2*y7+y2**2*y9)*y1**3+(y3**3*y7+y2**3*y9+y4**3*y8)*y1**2+ &
          (-y3**2*y8-y4**2*y7)*y2**3+(-y4**3*y7-y3**3*y8)*y2**2-y3**3*y4**2*y9- &
          y3**2*y4**3*y9
      dF(281) = (2._ark/3._ark*sqrt(3._ark)*y2**2*y5+(-sqrt(3._ark)*y5/3._ark+y6)*y3**2+(- &
          sqrt(3._ark)*y5/3._ark-y6)*y4**2)*y1**3+(2._ark/3._ark*sqrt(3._ark)*y2**3*y5+(- &
          sqrt(3._ark)*y5/3._ark+y6)*y3**3+(-sqrt(3._ark)*y5/3._ark-y6)*y4**3)*y1**2+((- &
          sqrt(3._ark)*y5/3._ark-y6)*y3**2+(-sqrt(3._ark)*y5/3._ark+y6)*y4**2)*y2**3+((- &
          sqrt(3._ark)*y5/3._ark-y6)*y3**3+(-sqrt(3._ark)*y5/3._ark+y6)*y4**3)*y2**2+2._ark/ &
          3._ark*sqrt(3._ark)*y3**3*y4**2*y5+2._ark/3._ark*sqrt(3._ark)*y3**2*y4**3*y5
      dF(282) = ((y4+y3)*y2**2+(y3**2+y4**2)*y2+y3*y4**2+y3**2*y4)*y1**3+((y4+ &
          y3)*y2**3+(y4**3+y3**3)*y2+y3**3*y4+y3*y4**3)*y1**2+((y3**2+y4**2)*y2**3+(y4**3+ &
          y3**3)*y2**2+y3**3*y4**2+y3**2*y4**3)*y1+(y3*y4**2+y3**2*y4)*y2**3+(y3**3*y4+ &
          y3*y4**3)*y2**2+(y3**3*y4**2+y3**2*y4**3)*y2
      dF(283) = y1**3*y2*y3*y4+(y2**3*y3*y4+(y3**3*y4+y3*y4**3)*y2)*y1
      dF(284) = (y8**2+y7**2+y9**2)*y1**4+(y8**2+y7**2+y9**2)*y2**4+(y8**2+y7**2+ &
          y9**2)*y3**4+(y8**2+y7**2+y9**2)*y4**4
      dF(285) = (y6**2+y5**2)*y1**4+(y6**2+y5**2)*y2**4+(y6**2+y5**2)*y3**4+(y6**2+ &
          y5**2)*y4**4
      dF(286) = (y3**2+y2**2+y4**2)*y1**4+(y4**4+y2**4+y3**4)*y1**2+(y3**2+ &
          y4**2)*y2**4+(y4**4+y3**4)*y2**2+y3**4*y4**2+y3**2*y4**4
      dF(287) = ((y4+y3)*y2+y3*y4)*y1**4+((y4+y3)*y2**4+(y4**4+y3**4)*y2+y3*y4**4+ &
          y3**4*y4)*y1+y2**4*y3*y4+(y3*y4**4+y3**4*y4)*y2
      dF(288) = (y9+y7+y8)*y1**5+(-y7+y9-y8)*y2**5+(-y8-y9+y7)*y3**5+(-y7-y9+y8)*y4**5
      dF(289) = y3**6+y4**6+y1**6+y2**6

     f = 0.0_ark

     do ipar = 1,parmax
       !
       k = term(ipar)
       !
       f = f + force(ipar)*dF(k)
       !
     enddo      

end function poten_xy4
