
      subroutine ddx (mode,hpsi,psi,dim,matrix,trafo)
      use iso_c_binding
      implicit none
      interface
         subroutine cddx(mode, hpsi, psi, dim, mat, trafo)
     .        bind(c,name="DDXCPP")
         use, intrinsic :: iso_c_binding
          integer (kind = c_int):: mode
          complex (kind = c_double), dimension(*) :: psi
          complex (kind = c_double), dimension(*) :: hpsi
          integer (kind = c_int) :: dim
          real (kind = c_double), dimension(*) :: mat
          real (kind = c_double), dimension(*) :: trafo
         end subroutine cddx
      end interface
      integer       mode, teil, dim, i, j
      complex*16    h psi(dim), psi(dim)
      real*8        matrix(dim,dim,*), trafo(dim,dim), ort(dim)

      call cddx(mode, hpsi, psi, dim, matrix, trafo)

      end

      subroutine kin(mode,hpsi,psi,dim,matrix,trafo)
      implicit none
      interface
         subroutine ckin(mode, hpsi, psi, dim, mat, trafo)
     .        bind(c,name="KINCPP")
         use, intrinsic :: iso_c_binding
          integer (kind = c_int):: mode
          complex (kind = c_double), dimension(*) :: psi
          complex (kind = c_double), dimension(*) :: hpsi
          integer (kind = c_int) :: dim
          real (kind = c_double), dimension(*) :: mat
          real (kind = c_double), dimension(*) :: trafo
         end subroutine ckin
      end interface
      integer       mode, teil, dim, i, j
      complex*16    h psi(dim), psi(dim)
      real*8        matrix(dim,dim,*), trafo(dim,dim), ort(dim)

      call ckin(mode, hpsi, psi, dim, matrix, trafo)

      end

c feed sysblock for mctdh
      subroutine fillpara(f,para_cpp)
      implicit none
      integer    f
      real*8     para(4, 4048), para_cpp(4, f)
      integer    sys(4), i, j

      common /sys/ sys, para

      do i=1, f
          do j=1,4
              para(j, i)=para_cpp(j, i)
          enddo
      enddo
      end

