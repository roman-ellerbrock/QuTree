
program poten
  use pesparam
  implicit none

  integer, parameter :: ark         = selected_real_kind(25,32)
  real(ark) :: deg, pi
  integer ipar, ieq, info
  double precision  f, local(10)
  double precision  alpha12, alpha13, alpha14, alpha23, alpha24, alpha34, r1, r2, r3, r4
  character(14) buf30
  integer i
  real*8 x(3, 5)
  real*8 poten_xy4

!  interface
!  function poten_xy4(local,parmax,par_,force,term)
!    integer parmax, term(289)
!    double precision local(10), par_(2),force(289), poten_xy4
!  end function poten_xy4
!  end interface

  pi = 4.0d0 * datan2(1.0d0,1.0d0)

  deg = pi/180.0d0

  !read parameters

  call potini_ch4_xy4()

  !read grid, print potential energy

!  do
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

    f = poten_xy4(local)
    write(*,'(10(1x,f12.7),2x,f14.6)') local(1:10), f
!    cycle
!4  exit
!  enddo

4  CONTINUE
  do i=1, 5
    read*, x(1, i), x(2, i), x(3, i)
    print*, x(1, i), x(2, i), x(3, i)
  enddo
  call carttoint(local, x)
  f = poten_xy4(local)
  write(*,'(10(1x,f12.7),2x,f14.6)') local(1:10), f

end program poten


