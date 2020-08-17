c######################################################################
c Libary module LA_lib
c
c Contains subroutines for linear algebra:
c "matvec", "rmatvec", "tmatvec", "ctmatvec", "trmatvec", 
c "vecvec", "rvecvec", "copyvec", "rcopyvec", 
c "nullvec", "rnullvec", "unitymatrix", "runitymatrix",
c "rhomat", "rrhomat", "realrhomat", "rhodiag", "orthonormal"
c "rorthonormal", "multiorthonormal", "rmultiorthonormal"
c
c V1.0, 13.4.1992
c######################################################################
c----------------------------------------------------------------------
c Libary subroutine matvec
c
c Multiplication of a complex tensor of third order with a complex
c matrix: 
c           mulpsi(i,j,k) = matrix(j,l)*psi(i,l,k) .
c
c If the label add is set to be ".true.", the initial tensor
c "mulpsi" is added to the result of the above multiplication.
c
c Input-variables:  psi    - complex tensor of third order
c                   matrix - complex matrix
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c                   add    - label as defined above
c                   mulpsi - complex tensor of third order 
c                            added to the result, if add=.true.
c Output-variables: mulpsi - resulting complex tensor
c----------------------------------------------------------------------

      subroutine matvec (mulpsi, psi, matrix, a, b, c, add)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, k, l
      complex*16    psi(b,a,c), mulpsi(b,a,c), matrix(a,a)
      logical*4     add

      if (.not.add) then
         do i=1,b*a*c
            mulpsi(i,1,1)=0
         enddo
      endif
      if (b.ne.1) then
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k,i)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  do i=1,b
                     mulpsi(i,k,j)=mulpsi(i,k,j)
     .                             +matrix(k,l)*psi(i,l,j)
                  enddo
               enddo
            enddo
         enddo
      else
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  mulpsi(1,k,j)=mulpsi(1,k,j)
     .                          +matrix(k,l)*psi(1,l,j)
               enddo
            enddo
         enddo
      endif
      return
      end

c----------------------------------------------------------------------
c Libary subroutine tmatvec
c
c Multiplication of a complex tensor of third order with a transposed
c complex matrix: 
c           mulpsi(i,j,k) = matrix(l,j)*psi(i,l,k) .
c
c If the label add is set to be ".true.", the initial tensor
c "mulpsi" is added to the result of the above multiplication.
c
c Input-variables:  psi    - complex tensor of third order
c                   matrix - complex matrix
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c                   add    - label as defined above
c                   mulpsi - complex tensor of third order 
c                            added to the result, if add=.true.
c Output-variables: mulpsi - resulting complex tensor
c----------------------------------------------------------------------

      subroutine tmatvec (mulpsi, psi, matrix, a, b, c, add)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, k, l
      complex*16    psi(b,a,c), mulpsi(b,a,c), matrix(a,a)
      logical*4     add

      if (.not.add) then
         do i=1,b*a*c
            mulpsi(i,1,1)=0
         enddo
      endif
      if (b.ne.1) then
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k,i)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  do i=1,b
                     mulpsi(i,k,j)=mulpsi(i,k,j)
     .                             +matrix(l,k)*psi(i,l,j)
                  enddo
               enddo
            enddo
         enddo
      else
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  mulpsi(1,k,j)=mulpsi(1,k,j)
     .                          +matrix(l,k)*psi(1,l,j)
               enddo
            enddo
         enddo
      endif
      return
      end

c----------------------------------------------------------------------
c Libary subroutine cmatvec
c
c Multiplication of a complex tensor of third order with a
c complex conjugate matrix: 
c           mulpsi(i,j,k) = conjg(matrix(j,l))*psi(i,l,k) .
c
c If the label add is set to be ".true.", the initial tensor
c "mulpsi" is added to the result of the above multiplication.
c
c Input-variables:  psi    - complex tensor of third order
c                   matrix - complex matrix
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c                   add    - label as defined above
c                   mulpsi - complex tensor of third order 
c                            added to the result, if add=.true.
c Output-variables: mulpsi - resulting complex tensor
c----------------------------------------------------------------------

      subroutine cmatvec (mulpsi, psi, matrix, a, b, c, add)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, k, l
      complex*16    psi(b,a,c), mulpsi(b,a,c), matrix(a,a)
      logical*4     add

      if (.not.add) then
         do i=1,b*a*c
            mulpsi(i,1,1)=0
         enddo
      endif
      if (b.ne.1) then
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k,i)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  do i=1,b
                     mulpsi(i,k,j)=mulpsi(i,k,j)
     .                             +conjg(matrix(k,l))*psi(i,l,j)
                  enddo
               enddo
            enddo
         enddo
      else
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  mulpsi(1,k,j)=mulpsi(1,k,j)
     .                          +conjg(matrix(k,l))*psi(1,l,j)
               enddo
            enddo
         enddo
      endif
      return
      end

c----------------------------------------------------------------------
c Libary subroutine ctmatvec
c
c Multiplication of a complex tensor of third order with a
c complex conjugate transposed matrix: 
c           mulpsi(i,j,k) = conjg(matrix(l,j))*psi(i,l,k) .
c
c If the label add is set to be ".true.", the initial tensor
c "mulpsi" is added to the result of the above multiplication.
c
c Input-variables:  psi    - complex tensor of third order
c                   matrix - complex matrix
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c                   add    - label as defined above
c                   mulpsi - complex tensor of third order 
c                            added to the result, if add=.true.
c Output-variables: mulpsi - resulting complex tensor
c----------------------------------------------------------------------

      subroutine ctmatvec (mulpsi, psi, matrix, a, b, c, add)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, k, l
      complex*16    psi(b,a,c), mulpsi(b,a,c), matrix(a,a)
      logical*4     add

      if (.not.add) then
         do i=1,b*a*c
            mulpsi(i,1,1)=0
         enddo
      endif
      if (b.ne.1) then
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k,i)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  do i=1,b
                     mulpsi(i,k,j)=mulpsi(i,k,j)
     .                             +conjg(matrix(l,k))*psi(i,l,j)
                  enddo
               enddo
            enddo
         enddo
      else
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  mulpsi(1,k,j)=mulpsi(1,k,j)
     .                          +conjg(matrix(l,k))*psi(1,l,j)
               enddo
            enddo
         enddo
      endif
      return
      end

c----------------------------------------------------------------------
c Libary subroutine rmatvec
c
c Multiplication of a real tensor of third order with a real
c matrix: 
c           mulpsi(i,j,k) = matrix(j,l)*psi(i,l,k) .
c
c If the label add is set to be ".true.", the initial tensor
c "mulpsi" is added to the result of the above multiplication.
c
c Input-variables:  psi    - real tensor of third order
c                   matrix - real matrix
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c                   add    - label as defined above
c                   mulpsi - real tensor of third order 
c                            added to the result, if add=.true.
c Output-variables: mulpsi - resulting real tensor
c----------------------------------------------------------------------

      subroutine rmatvec (mulpsi, psi, matrix, a, b, c, add)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, k, l
      real*8        psi(b,a,c), mulpsi(b,a,c), matrix(a,a)
      logical*4     add

      if (.not.add) then
         do i=1,b*a*c
            mulpsi(i,1,1)=0
         enddo
      endif
      if (b.gt.2) then
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k,i)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  do i=1,b
                     mulpsi(i,k,j)=mulpsi(i,k,j)
     .                             +matrix(k,l)*psi(i,l,j)
                  enddo
               enddo
            enddo
         enddo
      else
         if (b.eq.2) then
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
            do j=1,c
               do k=1,a
                  do l=1,a
                     mulpsi(1,k,j)=mulpsi(1,k,j)
     .                             +matrix(k,l)*psi(1,l,j)
                     mulpsi(2,k,j)=mulpsi(2,k,j)
     .                             +matrix(k,l)*psi(2,l,j)
                  enddo
               enddo
            enddo
         else
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
            do j=1,c
               do k=1,a
                  do l=1,a
                     mulpsi(1,k,j)=mulpsi(1,k,j)
     .                             +matrix(k,l)*psi(1,l,j)
                  enddo
               enddo
            enddo
         endif
      endif
      return
      end

c----------------------------------------------------------------------
c Libary subroutine trmatvec
c
c Multiplication of a real tensor of third order with a transposed real
c matrix: 
c           mulpsi(i,j,k) = matrix(l,j)*psi(i,l,k) .
c
c If the label add is set to be ".true.", the initial tensor
c "mulpsi" is added to the result of the above multiplication.
c
c Input-variables:  psi    - real tensor of third order
c                   matrix - real matrix
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c                   add    - label as defined above
c                   mulpsi - real tensor of third order 
c                            added to the result, if add=.true.
c Output-variables: mulpsi - resulting real tensor
c----------------------------------------------------------------------

      subroutine trmatvec (mulpsi, psi, matrix, a, b, c, add)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, k, l
      real*8        psi(b,a,c), mulpsi(b,a,c), matrix(a,a)
      logical*4     add

      if (.not.add) then
         do i=1,b*a*c
            mulpsi(i,1,1)=0
         enddo
      endif
      if (b.gt.2) then
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k,i)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
         do j=1,c
            do k=1,a
               do l=1,a
                  do i=1,b
                     mulpsi(i,k,j)=mulpsi(i,k,j)
     .                             +matrix(l,k)*psi(i,l,j)
                  enddo
               enddo
            enddo
         enddo
      else
         if (b.eq.2) then
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
            do j=1,c
               do k=1,a
                  do l=1,a
                     mulpsi(1,k,j)=mulpsi(1,k,j)
     .                             +matrix(l,k)*psi(1,l,j)
                     mulpsi(2,k,j)=mulpsi(2,k,j)
     .                             +matrix(l,k)*psi(2,l,j)
                  enddo
               enddo
            enddo
         else
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,l,k)
!$omp& schedule(static)
!$omp& if((c*a*a*b.ge.effort).and.(c.gt.1))
            do j=1,c
               do k=1,a
                  do l=1,a
                     mulpsi(1,k,j)=mulpsi(1,k,j)
     .                             +matrix(l,k)*psi(1,l,j)
                  enddo
               enddo
            enddo
         endif
      endif
      return
      end

c----------------------------------------------------------------------
c Libary subroutine vecvec
c
c Multiplication of a complex tensor of third order with a complex
c vector: 
c           mulpsi(i,j,k) = vector(j)*psi(i,j,k) .
c
c Input-variables:  psi    - complex tensor of third order
c                   vector - complex vector
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c Output-variables: mulpsi - resulting complex tensor
c----------------------------------------------------------------------

      subroutine vecvec (mulpsi, psi, vector, a, b, c)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, k
      complex*16    psi(b,a,c), mulpsi(b,a,c), vector(a)

      if (b.ne.1) then
         do j=1,c
            do k=1,a
               do i=1,b
                  mulpsi(i,k,j)=vector(k)*psi(i,k,j)
               enddo
            enddo
         enddo
      else
         do j=1,c
            do k=1,a
               mulpsi(1,k,j)=vector(k)*psi(1,k,j)
            enddo
         enddo
      endif
      return
      end

c----------------------------------------------------------------------
c Libary subroutine cvecvec
c
c Multiplication of a complex tensor of third order with a complex
c conjugate vector: 
c           mulpsi(i,j,k) = conjg(vector(j))*psi(i,j,k) .
c
c Input-variables:  psi    - complex tensor of third order
c                   vector - complex vector
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c Output-variables: mulpsi - resulting complex tensor
c----------------------------------------------------------------------

      subroutine cvecvec (mulpsi, psi, vector, a, b, c)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, k
      complex*16    psi(b,a,c), mulpsi(b,a,c), vector(a)

      if (b.ne.1) then
         do j=1,c
            do k=1,a
               do i=1,b
                  mulpsi(i,k,j)=conjg(vector(k))*psi(i,k,j)
               enddo
            enddo
         enddo
      else
         do j=1,c
            do k=1,a
               mulpsi(1,k,j)=conjg(vector(k))*psi(1,k,j)
            enddo
         enddo
      endif
      return
      end

c----------------------------------------------------------------------
c Libary subroutine rvecvec
c
c Multiplication of a real tensor of third order with a real
c vector: 
c           mulpsi(i,j,k) = vector(j)*psi(i,j,k) .
c
c Input-variables:  psi    - real tensor of third order
c                   vector - real vector
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c Output-variables: mulpsi - resulting real tensor
c----------------------------------------------------------------------

      subroutine rvecvec (mulpsi, psi, vector, a, b, c)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, k
      real*8        psi(b,a,c), mulpsi(b,a,c), vector(a)

      if (b.gt.2) then
         do j=1,c
            do k=1,a
               do i=1,b
                  mulpsi(i,k,j)=vector(k)*psi(i,k,j)
               enddo
            enddo
         enddo
      else
         if (b.eq.2) then
            do j=1,c
               do k=1,a
                  mulpsi(1,k,j)=vector(k)*psi(1,k,j)
                  mulpsi(2,k,j)=vector(k)*psi(2,k,j)
               enddo
            enddo
         else
            do j=1,c
               do k=1,a
                  mulpsi(1,k,j)=vector(k)*psi(1,k,j)
               enddo
            enddo
         endif
      endif
      return
      end


c----------------------------------------------------------------------
c Libary subroutine copyvec
c----------------------------------------------------------------------

      subroutine copyvec (copypsi, psi, dim)

      implicit none

      integer       dim, i
      complex*16    psi(dim), copypsi(dim)

      do i=1,dim
         copypsi(i)=psi(i)
      enddo
      end

c----------------------------------------------------------------------
c Libary subroutine rcopyvec
c----------------------------------------------------------------------

      subroutine rcopyvec (copypsi, psi, dim)

      implicit none

      integer       dim, i
      real*8        psi(dim), copypsi(dim)

      do i=1,dim
         copypsi(i)=psi(i)
      enddo
      end

c----------------------------------------------------------------------
c Libary subroutine nullvec
c----------------------------------------------------------------------

      subroutine nullvec (nullpsi, dim)

      implicit none

      integer       dim, i
      complex*16    nullpsi(dim)

      do i=1,dim
         nullpsi(i)=0
      enddo
      end

c----------------------------------------------------------------------
c Libary subroutine rnullvec
c----------------------------------------------------------------------

      subroutine rnullvec (nullpsi, dim)

      implicit none

      integer       dim, i
      real*8        nullpsi(dim)

      do i=1,dim
         nullpsi(i)=0
      enddo
      end

c----------------------------------------------------------------------
c Libary subroutine unitymatrix
c----------------------------------------------------------------------

      subroutine unitymatrix (matrix,dim)

      implicit none

      integer       dim, i
      complex*16    matrix(dim,dim)

      call nullvec(matrix,dim*dim)
      do i=1,dim
         matrix(i,i)=(1.d0,0.d0)
      enddo
      end

c----------------------------------------------------------------------
c Libary subroutine runitymatrix
c----------------------------------------------------------------------

      subroutine runitymatrix (matrix,dim)

      implicit none

      integer       dim, i
      real*8        matrix(dim,dim)

      call rnullvec(matrix,dim*dim)
      do i=1,dim
         matrix(i,i)=1.d0
      enddo
      end

c----------------------------------------------------------------------
c Libary subroutine rhomat
c
c Multiplication of two complex tensors of third order:
c
c           matrix(j,l) = conjg(bra(i,j,k))*ket(i,l,k) .
c
c Input-variables:  bra, ket  - complex tensors of third order
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c Output-variables: matrix - complex matrix
c----------------------------------------------------------------------


      subroutine rhomat (bra,ket,matrix,a,b,c)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, n, m
      complex*16    bra(b,a,c), ket(b,a,c), matrix(a,a)

      do i=1,a
         do j=1,a
            matrix(j,i)=0
         enddo
      enddo
      if (b.ne.1) then
         do m=1,c
         do i=1,a
         do j=1,a
            do n=1,b
               matrix(j,i)=matrix(j,i)+conjg(bra(n,j,m))*ket(n,i,m)
            enddo
         enddo
         enddo
         enddo
      else
         if (c.gt.a) then
            do i=1,a
            do j=1,a
               do m=1,c
                  matrix(j,i)=matrix(j,i)+conjg(bra(1,j,m))*ket(1,i,m)
               enddo
            enddo
            enddo
         else
            do m=1,c
               do i=1,a
               do j=1,a
                  matrix(j,i)=matrix(j,i)+conjg(bra(1,j,m))*ket(1,i,m)
               enddo
            enddo
            enddo
         endif
      endif
      return
      end


c----------------------------------------------------------------------
c Libary subroutine realrhomat
c
c Multiplication of two complex tensors of third order:
c
c           matrix(j,l) = bra(i,j,k)*ket(i,l,k) .
c
c Input-variables:  bra, ket  - complex tensors of third order
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c Output-variables: matrix - complex matrix
c----------------------------------------------------------------------


      subroutine realrhomat (bra,ket,matrix,a,b,c)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, n, m
      complex*16    bra(b,a,c), ket(b,a,c), matrix(a,a)

      do i=1,a
         do j=1,a
            matrix(j,i)=0
         enddo
      enddo
      if (b.ne.1) then
         do m=1,c
         do i=1,a
         do j=1,a
            do n=1,b
               matrix(j,i)=matrix(j,i)+bra(n,j,m)*ket(n,i,m)
            enddo
         enddo
         enddo
         enddo
      else
         if (c.gt.a) then
            do i=1,a
            do j=1,a
               do m=1,c
                  matrix(j,i)=matrix(j,i)+bra(1,j,m)*ket(1,i,m)
               enddo
            enddo
            enddo
         else
            do m=1,c
               do i=1,a
               do j=1,a
                  matrix(j,i)=matrix(j,i)+bra(1,j,m)*ket(1,i,m)
               enddo
            enddo
            enddo
         endif
      endif
      return
      end


c----------------------------------------------------------------------
c Libary subroutine rrhomat
c
c Multiplication of two real tensors of third order:
c
c           matrix(j,l) = bra(i,j,k))*ket(i,l,k) .
c
c Input-variables:  bra, ket  - real tensors of third order
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c Output-variables: matrix - real matrix
c----------------------------------------------------------------------


      subroutine rrhomat (bra,ket,matrix,a,b,c)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, j, n, m
      real*8        bra(b,a,c), ket(b,a,c), matrix(a,a)

      do i=1,a
         do j=1,a
            matrix(j,i)=0
         enddo
      enddo
      if (b.ne.1) then
         do m=1,c
         do i=1,a
         do j=1,a
            do n=1,b
               matrix(j,i)=matrix(j,i)+bra(n,j,m)*ket(n,i,m)
            enddo
         enddo
         enddo
         enddo
      else
         if (c.gt.a) then
            do i=1,a
            do j=1,a
               do m=1,c
                  matrix(j,i)=matrix(j,i)+bra(1,j,m)*ket(1,i,m)
               enddo
            enddo
            enddo
         else
            do m=1,c
               do i=1,a
               do j=1,a
                  matrix(j,i)=matrix(j,i)+bra(1,j,m)*ket(1,i,m)
               enddo
            enddo
            enddo
         endif
      endif
      return
      end


c----------------------------------------------------------------------
c Libary subroutine rhodiag
c
c Multiplication of two complex tensors of third order:
c
c           vector(j) = conjg(bra(i,j,k))*ket(i,j,k) .
c
c Input-variables:  bra, ket  - complex tensors of third order
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c Output-variables: vector - complex vector
c----------------------------------------------------------------------


      subroutine rhodiag (bra,ket,vector,a,b,c)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, n, m
      complex*16    bra(b,a,c), ket(b,a,c), vector(a)

      do i=1,a
         vector(i)=0
      enddo
      if (b.ne.1) then
         do m=1,c
            do i=1,a
               do n=1,b
                  vector(i)=vector(i)+conjg(bra(n,i,m))*ket(n,i,m)
               enddo
            enddo
         enddo
      else
         if (c.gt.a) then
            do i=1,a
               do m=1,c
                  vector(i)=vector(i)+conjg(bra(1,i,m))*ket(1,i,m)
               enddo
            enddo
         else
            do m=1,c
               do i=1,a
                  vector(i)=vector(i)+conjg(bra(1,i,m))*ket(1,i,m)
               enddo
            enddo
         endif
      endif
      return
      end


c----------------------------------------------------------------------
c Libary subroutine rrhodiag
c
c Multiplication of two real tensors of third order:
c
c           vector(j) = bra(i,j,k)*ket(i,j,k) .
c
c Input-variables:  bra, ket  - real tensors of third order
c                   a      - dimension of second index
c                   b      - dimension of first index
c                   c      - dimension of third index
c Output-variables: vector - real vector
c----------------------------------------------------------------------


      subroutine rrhodiag (bra,ket,vector,a,b,c)

      implicit none

      integer       effort
      parameter     (effort=1000)
      integer       a, b, c, i, n, m
      real*8        bra(b,a,c), ket(b,a,c), vector(a)

      do i=1,a
         vector(i)=0
      enddo
      if (b.ne.1) then
         do m=1,c
            do i=1,a
               do n=1,b
                  vector(i)=vector(i)+bra(n,i,m)*ket(n,i,m)
               enddo
            enddo
         enddo
      else
         if (c.gt.a) then
            do i=1,a
               do m=1,c
                  vector(i)=vector(i)+bra(1,i,m)*ket(1,i,m)
               enddo
            enddo
         else
            do m=1,c
               do i=1,a
                  vector(i)=vector(i)+bra(1,i,m)*ket(1,i,m)
               enddo
            enddo
         endif
      endif
      return
      end

c----------------------------------------------------------------------
c Libary subroutine orthonormal
c
c Gram-Schmidt orthogonalisation and normalisation of a set of 
c complex vectors vector(dim,zahl).
c
c Input-variables:  vector - set of real vectors to be normalized
c                   dim    - dimension of individual vectors
c                   zahl   - number of vectors
c Output-variables: vector - normalized set of real vectors
c----------------------------------------------------------------------

      subroutine orthonormal (vector, dim, zahl)

      integer       effort, maxstep
      parameter     (effort=1000, maxstep=10)
      integer       dim, zahl, n, i, j, step
      complex*16    vector(dim,zahl), x
      real*8        y, z, ratio
      
      ratio=0
      step=1
      do while ((ratio.lt.1d-4).and.(step.lt.maxstep))
         ratio=1
         do n=1,zahl
            z=0
            do i=1,dim
               z=z+dconjg(vector(i,n))*vector(i,n)
            enddo
            do j=1,n-1
               x=0
               do i=1,dim
                  x=x+dconjg(vector(i,j))*vector(i,n)
               enddo
               do i=1,dim
                  vector(i,n)=vector(i,n)-x*vector(i,j)
               enddo
            enddo
            y=0
            do i=1,dim
               y=y+dconjg(vector(i,n))*vector(i,n)
            enddo
            if (y/z.lt.ratio) ratio=y/z
            y=1/dsqrt(y)
            do i=1,dim
               vector(i,n)=vector(i,n)*y
            enddo
         enddo
         step=step+1
      enddo
      if (step.ge.maxstep) then
         write (6,*) 'orthonormalization failed'
         stop
      endif

      return
      end

c----------------------------------------------------------------------
c Libary subroutine rorthonormal
c
c Gram-Schmidt orthogonalisation and normalisation of a set of vectors
c vector(dim,zahl).
c
c Input-variables:  vector - set of real vectors to be normalized
c                   dim    - dimension of individual vectors
c                   zahl   - number of vectors
c Output-variables: vector - normalized set of real vectors
c----------------------------------------------------------------------

      subroutine rorthonormal (vector, dim, zahl)
      
      integer       effort, maxstep
      parameter     (effort=1000, maxstep=10)
      integer       dim, zahl, n, i, j, step
      real*8        vector(dim,zahl), x
      real*8        y, z, ratio
      
      ratio=0
      step=1
      do while ((ratio.lt.1d-4).and.(step.lt.maxstep))
         ratio=1
         do n=1,zahl
            z=0
            do i=1,dim
               z=z+vector(i,n)*vector(i,n)
            enddo
            do j=1,n-1
               x=0
               do i=1,dim
                  x=x+vector(i,j)*vector(i,n)
               enddo
               do i=1,dim
                  vector(i,n)=vector(i,n)-x*vector(i,j)
               enddo
            enddo
            y=0
            do i=1,dim
               y=y+vector(i,n)*vector(i,n)
            enddo
            if (y/z.lt.ratio) ratio=y/z
            y=1/dsqrt(y)
            do i=1,dim
               vector(i,n)=vector(i,n)*y
            enddo
         enddo
         step=step+1
      enddo
      if (step.ge.maxstep) then
         write (6,*) 'orthonormalization failed'
         stop
      endif

      return
      end

c----------------------------------------------------------------------
c Libary subroutine multiorthonormal
c
c Gram-Schmidt orthogonalisation and normalisation of a set of 
c complex vectors vector(dim,zahl,1) and corresponding transformation
c of the following sets of complex vectors vectors(dim,zahl,nr).
c
c Input-variables:  vector - set of real vectors to be normalized
c                   dim    - dimension of individual vectors
c                   zahl   - number of vectors
c Output-variables: vector - normalized set of real vectors
c----------------------------------------------------------------------

      subroutine multiorthonormal (vector, dim, zahl, nr)

      integer       effort
      parameter     (effort=1000)
      integer       dim, zahl, n, i, j, k
      complex*16    vector(dim,zahl,nr), x
      real*8        y
      
      do n=1,zahl
         do j=1,n-1
            x=0
            do i=1,dim
               x=x+dconjg(vector(i,j,1))*vector(i,n,1)
            enddo
            do k=1,nr
               do i=1,dim
                  vector(i,n,k)=vector(i,n,k)-x*vector(i,j,k)
               enddo
            enddo
         enddo
         y=0
         do i=1,dim
            y=y+dconjg(vector(i,n,1))*vector(i,n,1)
         enddo
         y=1/dsqrt(y)
         do k=1,nr
            do i=1,dim
               vector(i,n,k)=vector(i,n,k)*y
            enddo
         enddo
      enddo

      return
      end

c----------------------------------------------------------------------
c Libary subroutine rmultiorthonormal
c
c Gram-Schmidt orthogonalisation and normalisation of a set of vectors
c vector(dim,zahl) and corresponding transformation of the following 
c sets of vectors vectors(dim,zahl,nr).
c
c Input-variables:  vector - set of real vectors to be normalized
c                   dim    - dimension of individual vectors
c                   zahl   - number of vectors
c Output-variables: vector - normalized set of real vectors
c----------------------------------------------------------------------

      subroutine rmultiorthonormal (vector, dim, zahl, nr)

      integer       effort
      parameter     (effort=1000)
      integer       dim, zahl, n, i, j, k
      real*8        vector(dim,zahl,nr), x
      
      do n=1,zahl
         do j=1,n-1
            x=0
            do i=1,dim
               x=x+vector(i,j,1)*vector(i,n,1)
            enddo
            do k=1,nr
               do i=1,dim
                  vector(i,n,k)=vector(i,n,k)-x*vector(i,j,k)
               enddo
            enddo
         enddo
         x=0
         do i=1,dim
            x=x+vector(i,n,1)*vector(i,n,1)
         enddo
         x=1/dsqrt(x)
         do k=1,nr
            do i=1,dim
               vector(i,n,k)=vector(i,n,k)*x
            enddo
         enddo
      enddo

      return
      end

!c----------------------------------------------------------------------
!c Libary subroutine symorthonormal
!c
!c Symmetric (Loewdin) orthogonalisation and normalisation of a set of
!c complex vectors vector(dim,zahl).
!c
!c Input-variables:  vector    - set of complex vectors to be normalized
!c                   dim       - dimension of individual vectors
!c                   copytrans - if .true., copy transformation matrix
!c                                  S^(1/2) into s array
!c Output-variables: vector    - normalized set of complex vectors
!c                   s         - transformation matrix S^(1/2)
!c----------------------------------------------------------------------
!
!      subroutine symorthonormal(vector, dim, zahl, copytrans, s)
!
!      implicit none
!
!      integer       i,j,k,m,n
!      integer       dim, zahl
!      complex*16    vector(dim,zahl), ev(zahl,zahl)
!      complex*16    v2(dim,zahl)
!      complex*16    s(zahl,zahl), s1(zahl,zahl), s2(zahl,zahl)
!      real*8        ew(zahl), ew2(zahl)
!      logical*4     copytrans
!
!c Diagonalize overlap matrix s
!      call rhomat(vector,vector,s1,zahl,dim,1)
!      call cDiag(s1,zahl,zahl,zahl,ew,ev,zahl,zahl,.false.)
!
!c Square root (inverse) matrix
!      do i=1,zahl
!         ew2(i) = dsqrt(ew(i))
!         ew(i) = 1/dsqrt(ew(i))
!      enddo
!
!c Backtransform s matrix
!      do i=1,zahl
!         do j=1,zahl
!            s1(j,i) = 0d0
!            s2(j,i) = 0d0
!            do n=1,zahl
!               s1(j,i) = s1(j,i) + ev(i,n) * ew(n) * dconjg(ev(j,n))
!               s2(j,i) = s2(j,i) + ev(i,n) * ew2(n) * dconjg(ev(j,n))
!            enddo
!         enddo
!      enddo
!
!c Transform vectors
!      call matvec(v2,vector,s1,zahl,dim,1,.false.)
!      do i=1,zahl
!         do j=1,dim
!            vector(j,i) = v2(j,i)
!         enddo
!      enddo
!
!c Copy transformation matrix
!      if( copytrans ) then
!         do i=1,zahl
!            do j=1,zahl
!               s(j,i) = s2(i,j)
!            enddo
!         enddo
!      endif
!
!      end
!