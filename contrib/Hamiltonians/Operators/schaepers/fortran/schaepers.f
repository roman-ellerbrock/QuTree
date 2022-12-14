c######################################################################
c MCTDH-Modul Operatoren
c
c Contains subroutines "h init" and "h".
c "h init" initializes the arrays "coeff" and "diag" which
c specify the coefficients employed in the Hamiltonian series.
c "h" calculates the action of a single-particle operator (on the
c primitive grid). "mode" and "teil" define the number of the
c coordinate and of the term in the series, respectively.
c
c generic KEO
c DS May 2018
c######################################################################


      subroutine hinitschaepers(coeff,diag,dim,na,mass,coupling)

      implicit none

      integer                           dim, number_of_coefficients
      integer                           N, k, maxN,na
      parameter                         (maxN=20)
      integer                           i, coeff_per_part(4:maxN)
      character*20                      coord(5:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      logical*4                         diag(3)
      complex*16                        coeff(*)
      integer                           coupling(na-4)
      
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart
      SAVE      /coordinate_choice/
      integer   njx(5:maxN), njy(5:maxN), njz(5:maxN)
      common    /nj/                    njx, njy, njz
      SAVE      /nj/
      real*8    mass(maxN)
      
        N = na     !number of atoms in your system
        
        if (N.gt.maxN) then
        print*, 'ERROR: maximum number of atoms is 20'
        print*, 'yours is: ', N
        stop
        endif

      do i = 1,na-4
        if(coupling(1).eq.0) then
            coord(4+i) = 'stereo'
        elseif (coupling(1).eq.1) then
            coord(4+i) = 'spherical'
        else
            write(*,*) 'coupling type not supported'
        endif
      enddo

c       generic part starts here
c       don't change anything of the rest of the file unless you want to
c       change the definition of the kinetic energy operator
        call koordreadschaepers(mass,na)

        do i=5, N
                if (coord(i).eq.'stereo') then
                        print*, 'stereographic coordinates for atom ', i
                else if (coord(i).eq.'spherical') then
                        print*, 'spherical coordinates for atom', i
                else
                        print*, 'ERROR: not supported: ', coord(i)
                        print*, 'coordinates for atom:', i
                        stop
                endif
        enddo
        

        number_of_coefficients = 0

c      Methyl-Potential
       call set_coefficients_methyl(coeff, number_of_coefficients)
       
       coeff_per_part(4) = number_of_coefficients

c      potential for atom 5 (no mixed terms between atoms not belonging
c      to ch3)

       call set_coefficients_atom(coeff, number_of_coefficients, 5)

       coeff_per_part(5) = number_of_coefficients - coeff_per_part(4)
c      potential for atom i (i=6 to N) (mixed terms to atoms 5 to i-1
c      included)
        
       do i = 6, N
                call set_coefficients_atom_mix(coeff,
     .                  number_of_coefficients, i)
                coeff_per_part(i) = number_of_coefficients
               do k = 4, i-1       
               coeff_per_part(i) = coeff_per_part(i) - coeff_per_part(k)
               enddo
       enddo

       coeff(number_of_coefficients+1)=0.d0

c Bestimmung der diagonalen Komponenten des Hamilton-Operators
      call diagon(diag,dim,number_of_coefficients)
      
      end

c----------------------------------------------------------------
      subroutine diagon(diag,dim,number_of_coefficients)

      implicit none

      integer     dim, number_of_coefficients, i, k, N, part
      logical*4   diag(dim,number_of_coefficients)


        do i = 1, number_of_coefficients
           do k= 1, dim
                diag(k,i)=.true.
           enddo
        enddo
        
      N = (dim - 6)/3 + 4


      call set_diag_methyl(diag, dim, i)  
      call set_diag_atom(diag, dim, i, 5)

      do part = 6, N
        call set_diag_atom_mix(diag, dim, i, part)
      enddo

      end


c#####################################################################

      subroutine h(mode,teilin,hpsi,psi,dim,matrix,trafo,ort)

      implicit none

      integer       mode, teilin, dim, i
      complex*16    hpsi(dim), psi(dim)
      complex*16    work(1024),work2(1024),work3(1024)
      real*8        matrix(dim,dim), trafo(dim,dim), ort(dim)

      real*8        xphi,xchi,r0

      integer       term, teil, part, index

      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      character*20                      coord(5:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

      integer      dummy(4)
      real*8       baspar(4,12)
      common /sys/ dummy, baspar, r0

      
        i = 4
        teil = teilin
        
! h=theta(r-r0)  dividing surface
      if ((teilin.eq.-1)) then
        if (mode.eq.7) then
          do i=1,dim
            if (ort(i).lt.r0) then
              hpsi(i)=0.0d0
            else
              hpsi(i)=psi(i)
            endif
          enddo
        else
          do i=1,dim
            hpsi(i)=psi(i)
          enddo
        endif
        return
      endif

! [H,h(r-r0)] flux operator
      if ((teilin.eq.-2)) then
        if (mode.eq.7) then
          call kin(mode,hpsi,psi,dim,matrix,trafo)
          do i=1,dim
            if (ort(i).lt.r0) hpsi(i)=0.0d0
          enddo
          do i=1,dim
            if (ort(i).lt.r0) then
               work(i)=0.0d0
            else
               work(i)=psi(i)
            endif
          enddo
          call kin(mode,work2,work,dim,matrix,trafo)
          do i=1,dim
            hpsi(i)=work2(i)-hpsi(i)
          enddo
        else
          do i=1,dim
            hpsi(i)=psi(i)
          enddo
        endif
        return
      endif
        
        do while (teil.gt.0)
                teil = teil-coeff_per_part(i)
                i = i +1
        enddo

        part = i-1
        term = teil + coeff_per_part(part)

        if (part.eq.4) then
                call get_index_methyl(index, term, mode)
        else if (part.eq.5) then
                call get_index_atom(index, term, mode, part)
        else if (part.gt.5) then
                call get_index_atom_mix(index, term, mode, part)
        else
                stop
        endif
        
c        print*, 'Teil: ', teilin, ' Mode: ', mode, ' Index: ', index
        

        SELECT CASE(index)
                CASE( 2) ! H = 1/x**2
                        do i=1,dim
                                hpsi(i)=psi(i)/ort(i)**2
                        end do
                        return
                CASE( 3) ! H = sin(x)
                        do i=1,dim
                                hpsi(i)=psi(i)*sin(ort(i))
                        end do
                        return
                CASE( 4) ! H = cos(x)
                        do i=1,dim
                                hpsi(i)=psi(i)*cos(ort(i))
                        end do
                        return
                CASE( 5) ! H = cot(x)
                        do i=1,dim
                                hpsi(i)=psi(i)*cos(ort(i))/sin(ort(i))
                        end do
                        return
                CASE( 6) ! H = sin(x)**2
                        do i=1,dim
                                hpsi(i)=psi(i)*sin(ort(i))**2
                        end do
                        return
                CASE( 7) ! H = cos(x)**2
                        do i=1,dim
                                hpsi(i)=psi(i)*cos(ort(i))**2
                        end do
                        return
                CASE( 8) ! H = cot(x)**2
                        do i=1,dim
                                hpsi(i)=psi(i)/tan(ort(i))**2
                        end do
                        return
                CASE( 9) ! H = cos(x)**2 - sin(x)**2
                        do i=1,dim
                        hpsi(i)=psi(i)*(cos(ort(i))**2-sin(ort(i))**2)
                        end do
                        return
                CASE(10) ! H = 1/sin(x)**2
                        do i=1,dim
                                hpsi(i)=psi(i)/sin(ort(i))**2
                        end do
                        return
                CASE(11) ! H = 5 - 3/sin(x)**2
                        do i=1,dim
                        hpsi(i)=psi(i)*(1.0d0+1.0d0/sin(ort(i))**2)
                        end do
                        return
                CASE(12) ! H = 3 - 1/sin(x)**2
                        do i=1,dim
                        hpsi(i)=psi(i)*(3.0d0-1.0d0/sin(ort(i))**2)
                        end do
                        return
                CASE(13) ! H = 1/sin(x)**2 - 3/2
                        do i=1,dim
                        hpsi(i)=psi(i)*(1.0d0/sin(ort(i))**2-1.5d0)
                        end do
                        return
                CASE(14) ! H = d/dx 
                        call ddx (mode,hpsi,psi,dim,matrix,trafo)
                        return
                CASE(15) ! H = -1/2 * (d/dx)**2
                        call kin(mode,hpsi,psi,dim,matrix,trafo)
                        return
                CASE(16) ! H = cot(x) d/dx + 1/sin(x) d/dx cos(x)
                         ! Legendre
                        do i=1,dim
                                work(i)=psi(i)*cos(ort(i))/sin(ort(i))
                        end do
                        call ddx (mode,work2,work,dim,matrix,trafo)
                        do i=1,dim
                                work2(i)=work2(i)/sin(ort(i))
                        end do
                        do i=1,dim
                                work(i)=psi(i)/sin(ort(i))
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                        hpsi(i)=work2(i)+hpsi(i)*cos(ort(i))/sin(ort(i))
                        end do
                        return
                CASE(17) ! H = 1/sin(x) d/dx sin(x)  +  d/dx
                         !  Legendre 
                        do i=1,dim
                                work(i)=psi(i)/sin(ort(i))
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)/sin(ort(i))
                        end do
                        return
                CASE(18) ! H = d/dx sin(x)  +  sin(x) d/dx
                        do i=1,dim
                                work(i)=sin(ort(i))*psi(i)
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*sin(ort(i))
                        end do
                        return
                CASE(19) ! H = d/dx cos(x)  +  cos(x) d/dx
                        do i=1,dim
                                work(i)=cos(ort(i))*psi(i)
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*cos(ort(i))
                        end do
                        return
                CASE(20) ! H =  d/dx cot(x)  +  cot(x) d/dx
                        do i=1,dim
                                work(i)=psi(i)*cos(ort(i))/sin(ort(i))
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                        hpsi(i)=hpsi(i)+work(i)*cos(ort(i))/sin(ort(i))
                        end do
                        return
                CASE(21) ! H = d/dx sin(x) cos(x) + sin(x) cos(x) d/dx
                        do i=1,dim
                                work(i)=cos(ort(i))*sin(ort(i))*psi(i)
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                        hpsi(i)=hpsi(i)+work(i)*cos(ort(i))*sin(ort(i))
                        end do
                        return
                CASE(22) ! H = 1/sqrt( y(x) ) * d/dx * y(x) * d/dx * 1/sqrt( y(x) ) , y=x(phi)
                        do i=1,dim
                                work(i)=psi(i)/sqrt(xphi(ort(i)))
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                                work(i)=hpsi(i)*xphi(ort(i))
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)/sqrt(xphi(ort(i)))
                        end do
                        return 
                CASE(23) ! H = 1/sqrt( y(x) ) * d/dx * y(x) * d/dx * 1/sqrt( y(x) ) , y=x(chi)
                        do i=1,dim
                                work(i)=psi(i)/sqrt(xchi(ort(i)))
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                                work(i)=hpsi(i)*xchi(ort(i))
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)/sqrt(xchi(ort(i)))
                        end do
                        return 
                CASE(24) ! H =  1/sin(x) d/dx sin(x) d/dx  + 2 cot(x)**2  - 1/sin(x) (d/dx)**2 sin(x)  - (d/dx)**2
                         ! Legendre
                    call kin(mode,hpsi,psi,dim,matrix,trafo)
                    do i=1,dim
                    hpsi(i)=-2.0d0*hpsi(i)+2.0d0*psi(i)/(tan(ort(i))**2)
                    end do
                    call ddx(mode,work,psi,dim,matrix,trafo)
                    do i=1,dim
                        work(i)=work(i)/sin(ort(i))
                    end do
                    call ddx(mode,work2,work,dim,matrix,trafo)
                    do i=1,dim
                        hpsi(i)=hpsi(i)-work2(i)/sin(ort(i))
                    end do
                    do i=1,dim
                        work(i)=psi(i)/sin(ort(i))
                    end do
                    call ddx(mode,work2,work,dim,matrix,trafo)
                    do i=1,dim
                        work(i)=work2(i)/sin(ort(i))
                    end do
                    call ddx(mode,work2,work,dim,matrix,trafo)
                    do i=1,dim
                        hpsi(i)=hpsi(i)-work2(i)
                    end do
                    return    
                CASE(25) ! H =  -1/2 (d/dx)**2  +  1/2 (d/dx)**2 sin(x)**2  +  1/2 sin(x)**2 (d/dx)**2
                        do i=1,dim
                                work(i)=-sin(ort(i))**2*psi(i)
                        end do
                        call kin(mode,hpsi,work,dim,matrix,trafo)
                        call kin(mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)-work(i)*sin(ort(i))**2
                        end do
                        call kin(mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)
                        end do
                        return
                CASE(26) ! H =  -1/2 (d/dx)**2 cos(x)  -1/4 cos(x)  -1/2 cos(x) (d/dx)**2
                        do i=1,dim
                                work(i)=cos(ort(i))*psi(i)
                        end do
                        call kin(mode,hpsi,work,dim,matrix,trafo)
                        call kin(mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*cos(ort(i))
     .                                  -0.25d0*cos(ort(i))*psi(i)
                        end do
                        return
               CASE(27) ! H = sqrt( x(phi) * d/dx * 1/sqrt( x(phi) )  -  1/sqrt( x(phi) * d/dx * sqrt( x(phi) ) 
                        do i=1,dim
                                work(i)=psi(i)/sqrt(xphi(ort(i)))
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)*sqrt(xphi(ort(i)))
                        end do
                        do i=1,dim
                                work(i)=sqrt(xphi(ort(i)))*psi(i)
                        end do
                        call ddx (mode,work2,work,dim,matrix,trafo)
                        do i=1,dim
                           hpsi(i)=hpsi(i)-work2(i)/sqrt(xphi(ort(i)))
                        enddo
                        return
                CASE(28) ! H = x
                        do i=1,dim
                                hpsi(i)=psi(i)*ort(i)
                        enddo
                        return
                CASE(29) ! H = x^2
                        do i=1,dim
                                hpsi(i)=psi(i)*ort(i)*ort(i)
                        enddo
                        return
                CASE(30) ! H = (1-x^2)
                        do i=1,dim
                                hpsi(i)=psi(i)*(-ort(i)*ort(i)+1.0d0)
                        enddo
                        return
                CASE(31) ! H = (1-x^2)^2
                        do i=1,dim
                                hpsi(i)=psi(i)*(ort(i)**2-1.0d0)**2
                        enddo
                        return
                CASE(32) ! H = x*(1-x^2)
                        do i=1,dim
                           hpsi(i)=-psi(i)*(ort(i)*ort(i)-1.0d0)*ort(i)
                        enddo
                        return
                CASE(33) ! H = 1+x^2
                        do i=1,dim
                                hpsi(i)=psi(i)*(ort(i)*ort(i)+1.0d0)
                        enddo
                        return

                CASE(35) ! H = (1+x^2)^2
                        do i=1,dim
                                hpsi(i)=psi(i)*(ort(i)**2+1.0d0)**2
                        enddo
                        return
                CASE(36) ! H = x d/dx x
                        do i=1,dim
                                work(i)=ort(i)*psi(i)
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)*ort(i)
                        end do
                        return
                CASE(37) ! H =  d/dx x  +  x d/dx
                        do i=1,dim
                                work(i)=psi(i)*ort(i)
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*ort(i)
                        end do
                        return
                CASE(38) ! H =  d/dx x^3  +  x^3 d/dx
                        do i=1,dim
                                work(i)=psi(i)*ort(i)**3
                        end do
                        call ddx (mode,hpsi,work,dim,matrix,trafo)
                        call ddx (mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*ort(i)**3
                        end do
                        return
                CASE(39) ! H =  -1/2 x (d/dx)^2 x
                        do i=1,dim
                                work(i)=ort(i)*psi(i)
                        end do
                        call kin(mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)*ort(i)
                        end do
                        return
                CASE(40) ! H =  -1/2 x^2 (d/dx)^2 x^2
                        do i=1,dim
                                work(i)=ort(i)**2*psi(i)
                        end do
                        call kin(mode,hpsi,work,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)*ort(i)**2
                        end do
                        return
                CASE(41) ! H =  -1/2 (d/dx)^2 x  -1/2 x (d/dx)^2
                        do i=1,dim
                                work(i)=ort(i)*psi(i)
                        end do
                        call kin(mode,hpsi,work,dim,matrix,trafo)
                        call kin(mode,work,psi,dim,matrix,trafo)
                        do i=1,dim
                                hpsi(i)=hpsi(i)+work(i)*ort(i)
                        end do
                        return
                CASE DEFAULT
                        print*, 'ERROR: unknown index in h'
                        print*, 'teilin: ', teilin
                        print*, 'part: ', part
                        print*, 'term: ', term
                        print*, 'mode: ', mode
                        print*, 'index: ', index
                        stop
        END SELECT
                        


      end subroutine
c -------------------------------------------------------------------
      real*8 function xphi(phi)
      implicit none
 
      real*8     phi
      real*8     pi,pi3
      parameter (pi=3.14159265358979323844d0,pi3=pi/3.0d0)

      xphi=(1.0d0/( 1.0d0 +                   (phi-pi3)**2
     .                    - sqrt(3.0d0)/9.0d0*(phi-pi3)**3
     .                    +       3.0d0/4.0d0*(phi-pi3)**4  ))

      end 
c -------------------------------------------------------------------
      real*8 function xchi(chi)
      implicit none

      real*8     chi
      real*8     pi
      parameter (pi=3.14159265358979323844d0)
 
      xchi=( 1.0d0/( 1.0d0 + ((chi-pi)**2)/3.0d0
     .                     + ((chi-pi)**4)/12.0d0 ))

      end
c-------------------------------------------------------
        subroutine print_error(name, error, term, mode)

        implicit none
        integer         term, mode
        character*30    name,error

        print*, 'ERROR in KEO: ', name
        print*, error
        print*, 'term: ', term
        print*, 'mode: ', mode

        stop

        end subroutine
c -------------------------------------------------------------------

        subroutine get_index_methyl(index, term, mode)

        implicit none

        integer index, term, mode
        character*30    name, error2
        parameter       (name = 'get_index_methyl')
        parameter       (error2='mode unknown in this part')
        
        index = 0
 
        SELECT CASE (term)
                CASE(1)
                        if (mode.eq.1) index = 15
                CASE(2)
                        if (mode.eq.1) index = 2
                        if (mode.eq.2) index = 10
                        if (mode.eq.3) index = 15
                CASE(3)
                        if (mode.eq.1) index = 2
                        if (mode.eq.2) index = 15
                CASE(4)
                        if (mode.eq.1) index = 2
                        if (mode.eq.2) index = 11
                CASE(5)
                        if (mode.eq.1) index = 2
                        if (mode.eq.4) index = 12
                CASE(6)
                        if (mode.eq.1) index = 2
                        if (mode.eq.4) index = 15
                CASE(7)
                        if (mode.eq.1) index = 2
                        if (mode.eq.4) index = 10
                        if (mode.eq.5) index = 22
                CASE(8)
                        if (mode.eq.1) index = 2
                        if (mode.eq.5) index = 22
                CASE(9)
                        if (mode.eq.1) index = 2
                        if (mode.eq.4) index = 10
                        if (mode.eq.6) index = 23
                CASE(10)
                        if (mode.eq.1) index = 2
                        if (mode.eq.6) index = 23
                CASE(11)
                        if (mode.eq.1) index = 2
                        if (mode.eq.4) index = 20
                        if (mode.eq.5) index = 14
                CASE(12)
                        if (mode.eq.1) index = 2
                        if (mode.eq.4) index = 13
                        if (mode.eq.5) index = 27
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT
        
        if (index.eq.0) call print_error(name, error2, term, mode)
        end

c -------------------------------------------------------------------

        subroutine get_index_atom(index, term, mode, part)

        implicit none

        integer         index, term, mode, part, vecindex, tmp, i
        character*30    name, error2
        parameter       (name = 'get_index_atom')
        parameter       (error2='mode unknown in this part')

        integer                           N, maxN, block
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        
        index = 0
        block = 0

        i = 10       

        do while ((term.ge.blockstart(i,part,part)).and.(i.lt.19))
                i = i + 1
        enddo
        
        block = i - 1 
        

        if (block.le.0) then
                print*, 'Block: ', block
                call print_error(name, 'cannot find block',
     .                                  term, mode)
        endif    
      
        term = term - blockstart(block,part,part) + 1

        tmp = mode-6-(part-5)*3-1
        vecindex = mod(tmp,3)

        if ((block.ge.9).and.(block.lt.18)) then
                if (mode.eq.1) then
                        index = 2
                        return
                endif
        endif


        SELECT CASE (block)
                CASE(9)
                        if (mode.eq.4) then
                                 index = 14
                                 return
                        endif
                        call get_index_jy(index,term,vecindex,mode,part)
                CASE(10)
                        if (mode.eq.4) then
                                 index = 5
                                 return
                        endif
                        if (mode.eq.5) then
                                 index = 14
                                 return
                        endif
                        call get_index_jy(index,term,vecindex,mode,part)
                CASE(11)
                        if (mode.eq.4) then
                                 index = 5
                                 return
                        endif
                        if (mode.eq.6) then
                                 index = 14
                                 return
                        endif
                        call get_index_jx(index,term,vecindex,mode,part)
                CASE(12)
                        if (mode.eq.4) then
                                 index = 10
                                 return
                        endif
                        if (mode.eq.6) then
                                 index = 14
                                 return
                        endif
                        call get_index_jz(index,term,vecindex,mode,part)
                CASE(13)
                       call get_index_jx2(index,term,vecindex,mode,part)
                CASE(14)
                       call get_index_jy2(index,term,vecindex,mode,part)
                CASE(15)
                        if (mode.eq.4) then
                                index = 10
                                return
                        endif
                       call get_index_jz2(index,term,vecindex,mode,part)
                CASE(16)
                       call get_index_jz2(index,term,vecindex,mode,part)
                CASE(17)
                        if (mode.eq.4) then
                                index = 5
                                return
                        endif
                      call get_index_jxjz(index,term,vecindex,mode,part)
                CASE(18)
                        call get_index_t(index,term,vecindex,mode,part)
                CASE DEFAULT
                    call print_error(name,'block unknown', term, mode)
        END SELECT
        
        if (index.eq.0) call print_error(name, error2, term, mode)
        end

c -------------------------------------------------------------------

        subroutine get_index_atom_mix(index, term, mode, part)

        implicit none

        integer         index, mode, part, tmp, term, i
        integer         vib, vip, b
        logical*4       bmode, pmode
        character*30    name, error2
        parameter       (name = 'get_index_atom_mix')
        parameter       (error2='mode unknown in this part')

        integer                           N, maxN, block
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        integer         njx(5:maxN), njy(5:maxN), njz(5:maxN)
        common  /nj/    njx, njy, njz      
        
        index = 0
        block = 0
        b = 6
        
        i = 13
       
c       find out if mixed or part for this atom

        if (term.lt.blockstart(13,5,part)) then
                call get_index_atom(index, term, mode, part)
                return            
        endif

c       mixed
c       mixed with which atom?

        do while ((term.ge.blockstart(i,b,part)).and.(b.lt.part))
                b = b + 1
        enddo

        b = b - 1
c       which block?
        do while ((term.ge.blockstart(i,b,part)).and.(i.lt.20))
                i = i + 1
        enddo

        block = i - 1

        if (block.eq.18) block = 19

        if (block.le.12) then
                print*, 'Block: ', block
                call print_error(name, 'cannot find block',
     .                                  term, mode)    
        endif
       
        term = term - blockstart(block,b,part) + 1

        vip = -1
        vib = -1
        pmode = .false.
        bmode = .false.        

        tmp = mode-6-(part-5)*3-1
        if ((tmp.ge.0).and.(tmp.le.2)) then
                vip = mod(tmp,3)
                pmode = .true.
        endif
                
        tmp = mode-6-(b-5)*3-1
        if ((tmp.ge.0).and.(tmp.le.2)) then
                vib = mod(tmp,3)
                bmode = .true.
        endif

        if ((block.ge.13).and.(block.lt.18)) then
                if (mode.eq.1) then
                        index = 2
                        return
                endif
        endif

        if (block.eq.19) then
                if (mode.eq.1) then
                        index = 2
                        return
                endif
        endif                

        SELECT CASE (block)
                CASE(13)
                        tmp = mod(term, njx(b))
                        if (tmp.eq.0) tmp = njx(b)
                        if (bmode) then
                                term = tmp
                                call get_index_jx(index,term,vib,mode,b)
                                return
                        endif
                        if (pmode) then
                                term = (term - tmp) /njx(b) + 1
                             call get_index_jx(index,term,vip,mode,part)
                                return
                        endif
                CASE(14)
                        tmp = mod(term, njy(b))
                        if (tmp.eq.0) tmp = njy(b)
                        if (bmode) then
                                term = tmp
                                call get_index_jy(index,term,vib,mode,b)
                                return
                        endif
                        if (pmode) then
                                term = (term - tmp) /njy(b) + 1
                             call get_index_jy(index,term,vip,mode,part)
                                return
                        endif
                        return
                CASE(15)
                        if (mode.eq.4) then
                                index = 10
                                return
                        endif
                        tmp = mod(term, njz(b))
                        if (tmp.eq.0) tmp = njz(b)
                        if (bmode) then
                                term = tmp
                                call get_index_jz(index,term,vib,mode,b)
                                return
                        endif
                        if (pmode) then
                                term = (term - tmp) /njz(b) + 1
                             call get_index_jz(index,term,vip,mode,part)
                                return
                        endif
                        return
                CASE(16)
                        tmp = mod(term, njz(b))
                        if (tmp.eq.0) tmp = njz(b)
                        if (bmode) then
                                term = tmp
                                call get_index_jz(index,term,vib,mode,b)
                                return
                        endif
                        if (pmode) then
                                term = (term - tmp) /njz(b) + 1
                             call get_index_jz(index,term,vip,mode,part)
                                return
                        endif
                        return
                CASE(17)
                        if (mode.eq.4) then
                                index = 5
                                return
                        endif
                        tmp = mod(term, njx(b))
                        if (tmp.eq.0) tmp = njx(b)
                        if (bmode) then
                                term = tmp
                                call get_index_jx(index,term,vib,mode,b)
                                return
                        endif
                        if (pmode) then
                                term = (term - tmp) /njx(b) + 1
                             call get_index_jz(index,term,vip,mode,part)
                                return
                        endif
                        return
                CASE(19)
                        if (mode.eq.4) then
                                index = 5
                                return
                        endif
                        tmp = mod(term, njz(b))
                        if (tmp.eq.0) tmp = njz(b)
                        if (bmode) then
                                term = tmp
                                call get_index_jz(index,term,vib,mode,b)
                                return
                        endif
                        if (pmode) then
                                term = (term - tmp) /njz(b) + 1
                             call get_index_jx(index,term,vip,mode,part)
                                return
                        endif
                        return
                CASE DEFAULT
                    call print_error(name,'block unknown', term, mode)
        END SELECT
        
        if (index.eq.0) call print_error(name, error2, term, mode)
        end

c-----------------------------------------------------------
        subroutine get_index_jx(index, term, vecindex, mode, part)

        implicit none

        integer         index, term, vecindex, part, mode
        character*30    name, error2
        parameter       (name = 'get_index_jx')
        parameter       (error2='term unknown in this part')

        integer                           N, maxN
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        if (coord(part).eq.'stereo') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 37
                        if (vecindex.eq.2) index = 28
                CASE(2)
                        if (vecindex.eq.1) index = 30
                        if (vecindex.eq.2) index = 14
                CASE(3)
                        if (vecindex.eq.2) index = 36
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else if (coord(part).eq.'spherical') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 17
                        if (vecindex.eq.2) index = 3
                CASE(2)
                        if (vecindex.eq.1) index = 5
                        if (vecindex.eq.2) index = 19
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                
        end subroutine
c-----------------------------------------------------------
        subroutine get_index_jy(index, term, vecindex, mode, part)

        implicit none

        integer         index, term, vecindex, mode, part
        character*30    name, error2
        parameter       (name = 'get_index_jy')
        parameter       (error2='term unknown in this part')

        integer                           N, maxN
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        if (coord(part).eq.'stereo') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 28
                        if (vecindex.eq.2) index = 37
                CASE(2)
                        if (vecindex.eq.1) index = 14
                        if (vecindex.eq.2) index = 30
                CASE(3)
                        if (vecindex.eq.1) index = 36
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else if (coord(part).eq.'spherical') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 17
                        if (vecindex.eq.2) index = 4
                CASE(2)
                        if (vecindex.eq.1) index = 5
                        if (vecindex.eq.2) index = 18
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                
        end subroutine
c-----------------------------------------------------------
        subroutine get_index_jz(index, term, vecindex, mode, part)

        implicit none

        integer         index, term, mode, vecindex, part
        character*30    name, error2
        parameter       (name = 'get_index_jz')
        parameter       (error2='term unknown in this part')

        integer                           N, maxN
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        if (coord(part).eq.'stereo') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 28
                        if (vecindex.eq.2) index = 14
                CASE(2)
                        if (vecindex.eq.1) index = 14
                        if (vecindex.eq.2) index = 28
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else if (coord(part).eq.'spherical') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.2) index = 14
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                
        end subroutine
c-----------------------------------------------------------
        subroutine get_index_jx2(index, term, vecindex, mode, part)

        implicit none

        integer         index, term, mode, vecindex, part
        character*30    name, error2
        parameter       (name = 'get_index_jx2')
        parameter       (error2='term unknown in this part')

        integer                           N, maxN
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        if (coord(part).eq.'stereo') then
        SELECT CASE (term)
                CASE(2)
                        if (vecindex.eq.1) index = 29
                CASE(3)
                        if (vecindex.eq.2) index = 29
                CASE(4)
                        if (vecindex.eq.1) index = 39
                        if (vecindex.eq.2) index = 29
                CASE(5)
                        if (vecindex.eq.1) index = 38
                        if (vecindex.eq.2) index = 37
                CASE(6)
                        if (vecindex.eq.1) index = 37
                        if (vecindex.eq.2) index = 37
                CASE(7)
                        if (vecindex.eq.1) index = 37
                        if (vecindex.eq.2) index = 38
                CASE(8)
                        if (vecindex.eq.1) index = 31
                        if (vecindex.eq.2) index = 15
                CASE(9)
                        if (vecindex.eq.1) index = 30
                        if (vecindex.eq.2) index = 39
                CASE(10)
                        if (vecindex.eq.2) index = 40
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else if (coord(part).eq.'spherical') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 15
                CASE(2)
                        if (vecindex.eq.1) index = 8
                        if (vecindex.eq.2) index = 15
                CASE(3)
                        if (vecindex.eq.1) index = 24
                        if (vecindex.eq.2) index = 9
                CASE(4)
                        if (vecindex.eq.1) index = 16
                        if (vecindex.eq.2) index = 21
                CASE(5)
                        if (vecindex.eq.1) index = 8
                        if (vecindex.eq.2) index = 25
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                
        end subroutine

c-----------------------------------------------------------
        subroutine get_index_jy2(index, term, vecindex, mode, part)

        implicit none

        integer         index, term, vecindex, part, mode
        character*30    name, error2
        parameter       (name = 'get_index_jy2')
        parameter       (error2='term unknown in this part')

        integer                           N, maxN
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        if (coord(part).eq.'stereo') then
        SELECT CASE (term)
                CASE(2)
                        if (vecindex.eq.1) index = 29
                CASE(3)
                        if (vecindex.eq.2) index = 29
                CASE(4)
                        if (vecindex.eq.1) index = 40
                CASE(5)
                        if (vecindex.eq.1) index = 39
                        if (vecindex.eq.2) index = 30
                CASE(6)
                        if (vecindex.eq.1) index = 15
                        if (vecindex.eq.2) index = 31
                CASE(7)
                        if (vecindex.eq.1) index = 38
                        if (vecindex.eq.2) index = 37
                CASE(8)
                        if (vecindex.eq.1) index = 37
                        if (vecindex.eq.2) index = 37
                CASE(9)
                        if (vecindex.eq.1) index = 37
                        if (vecindex.eq.2) index = 38
                CASE(10)
                        if (vecindex.eq.1) index = 29
                        if (vecindex.eq.2) index = 39
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else if (coord(part).eq.'spherical') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 15
                CASE(2)
                        if (vecindex.eq.1) index = 8
                        if (vecindex.eq.2) index = 15
                CASE(3)
                        if (vecindex.eq.1) index = 24
                        if (vecindex.eq.2) index = 9
                CASE(4)
                        if (vecindex.eq.1) index = 16
                        if (vecindex.eq.2) index = 21
                CASE(5)
                        if (vecindex.eq.1) index = 8
                        if (vecindex.eq.2) index = 25
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                
        end subroutine

c-----------------------------------------------------------
        subroutine get_index_jz2(index, term, vecindex, mode, part)

        implicit none

        integer         index, term, vecindex, part, mode
        character*30    name, error2
        parameter       (name = 'get_index_jz2')
        parameter       (error2='term unknown in this part')

        integer                           N, maxN
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        if (coord(part).eq.'stereo') then
        SELECT CASE (term)
                CASE(2)
                        if (vecindex.eq.1) index = 15
                        if (vecindex.eq.2) index = 29
                CASE(3)
                        if (vecindex.eq.1) index = 29
                        if (vecindex.eq.2) index = 15
                CASE(4)
                        if (vecindex.eq.1) index = 37
                        if (vecindex.eq.2) index = 37
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else if (coord(part).eq.'spherical') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.2) index = 15
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                
        end subroutine

c-----------------------------------------------------------
        subroutine get_index_jxjz(index, term, vecindex, mode, part)

        implicit none

        integer         index, term, vecindex, part, mode
        character*30    name, error2
        parameter       (name = 'get_index_jxjz')
        parameter       (error2='term unknown in this part')

        integer                           N, maxN
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        if (coord(part).eq.'stereo') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 28
                CASE(2)
                        if (vecindex.eq.1) index = 32
                        if (vecindex.eq.2) index = 15
                CASE(3)
                        if (vecindex.eq.1) index = 41
                        if (vecindex.eq.2) index = 29
                CASE(4)
                        if (vecindex.eq.1) index = 28
                        if (vecindex.eq.2) index = 39
                CASE(5)
                        if (vecindex.eq.1) index = 36
                        if (vecindex.eq.2) index = 37
                CASE(6)
                        if (vecindex.eq.1) index = 14
                        if (vecindex.eq.2) index = 37
                CASE(7)
                        if (vecindex.eq.1) index = 14
                        if (vecindex.eq.2) index = 38
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else if (coord(part).eq.'spherical') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 17
                        if (vecindex.eq.2) index = 18
                CASE(2)
                        if (vecindex.eq.1) index = 5
                        if (vecindex.eq.2) index = 26
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                
        end subroutine

c-----------------------------------------------------------
        subroutine get_index_j2(index, term, vecindex, mode, part)

        implicit none

        integer         index, term, vecindex, part, mode
        character*30    name, error2
        parameter       (name = 'get_index_j2')
        parameter       (error2='term unknown in this part')

        integer                           N, maxN
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
        if (coord(part).eq.'stereo') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 15
                        if (vecindex.eq.2) index = 35
                CASE(2)
                        if (vecindex.eq.1) index = 35
                        if (vecindex.eq.2) index = 15
                CASE(3)
                        if (vecindex.eq.1) index = 39
                        if (vecindex.eq.2) index = 33
                CASE(4)
                        if (vecindex.eq.1) index = 33
                        if (vecindex.eq.2) index = 39
                CASE(5)
                        if (vecindex.eq.1) index = 40
                CASE(6)
                        if (vecindex.eq.2) index = 40
                CASE(7)
                        if (vecindex.eq.1) index = 33
                CASE(8)
                        if (vecindex.eq.2) index = 33
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else if (coord(part).eq.'spherical') then
        SELECT CASE (term)
                CASE(1)
                        if (vecindex.eq.1) index = 15
                CASE(2)
                        if (vecindex.eq.1) index = 10
                        if (vecindex.eq.2) index = 15
                CASE DEFAULT
                       call print_error(name,'term unknown', term, mode)
        END SELECT

        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                
        end subroutine

c-----------------------------------------------------------
        subroutine get_index_t(index, term, vecindex, mode, part)

        implicit none

        integer         index, term, vecindex, part, mode
        character*30    name, error2
        parameter       (name = 'get_index_t')
        parameter       (error2='term unknown in this part')

        integer                           N, maxN
        parameter                         (maxN=20)
        integer                           blockstart(9:19,5:maxN,5:maxN)
        integer                           coeff_per_part(4:maxN)
        character*20                      coord(5:maxN)
        common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                    blockstart
                
        if (vecindex.eq.0) then
                if (term.eq.1) then
                        index = 15
                else
                        index = 2
                endif
        else
                if (term.gt.1) then
                    term = term - 1
                    call get_index_j2(index, term, vecindex, mode, part)
                endif
        endif

        end subroutine

c-------------------------------------------------------------
      subroutine set_coefficients_methyl(coeff, number_of_coefficients)

      implicit none

      complex*16                        coeff(*), value
      integer                           number_of_coefficients
      
      real*8    a_th, a_phi, c_phi, a_chi, c_chi, a_th_phi
      integer   l, i


      call get_methyl_coefficients(a_th, a_phi, c_phi, a_chi, c_chi,
     .          a_th_phi)

      l = 0

      call set_coeff((1.d0,0.d0), l, coeff)
      call set_coeff((1.d0,0.d0), l, coeff)
      call set_coeff((1.d0,0.d0), l, coeff)
      call set_coeff((-0.125d0,0.d0), l, coeff)
      value = -3.d0/8.d0*a_th
      call set_coeff(value, l, coeff)
      value = a_th*(1.d0,0.d0)
      call set_coeff(value, l, coeff)
      value = -0.5d0*c_phi
      call set_coeff(value, l, coeff)
      value = -0.5d0*(a_phi-c_phi)
      call set_coeff(value, l, coeff)
      value = -0.5d0*c_chi
      call set_coeff(value, l, coeff)
      value = -0.5d0*(a_chi-c_chi)
      call set_coeff(value, l, coeff)
      value = -0.5d0*a_th_phi
      call set_coeff(value, l, coeff)
      value = -0.5d0*a_th_phi
      call set_coeff(value, l, coeff)
 
      number_of_coefficients = l

      end subroutine
c-------------------------------------------------------------
      subroutine set_coefficients_atom(coeff, number_of_coefficients,
     .                                  part)

      implicit none

      complex*16                        coeff(*), value
      integer                           number_of_coefficients
      integer                           N, maxN, part
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart
      
      real*8    a_th, a_phi, c_phi, a_chi, c_chi, a_th_phi
      real*8   a_th_y, b_phi_y, b_chi_x, a_chi_z, c_chi_z, a_x
      real*8   a_y, a_z, c_z, b_x_z 
      integer   l, nb, i
      complex*16        ca, cb(1:10)

        call get_methyl_coefficients(a_th, a_phi, c_phi, a_chi, c_chi,
     .          a_th_phi)

        call get_atom_coefficients(a_th_y, b_phi_y, b_chi_x,
     .          a_chi_z, c_chi_z, a_x, a_y, a_z, c_z, b_x_z)

      l = number_of_coefficients

      call get_jy_coefficients(coord(part),cb,nb,part)
      ca = 1.d0/2.d0*2.d0*(0.d0,1.d0)*a_th_y
      blockstart(9,part,part) = l -number_of_coefficients + 1
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo


      blockstart(10,part,part) = l -number_of_coefficients + 1
      ca = 1.d0/2.d0*2.d0*(0.d0,1.d0)*b_phi_y
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo

      call get_jx_coefficients(coord(part),cb,nb,part)
      blockstart(11,part,part) = l -number_of_coefficients + 1
      ca = 1.d0/2.d0*2.d0*(0.d0,1.d0)*b_chi_x
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo

      call get_jz_coefficients(coord(part),cb,nb,part)
      blockstart(12,part,part) = l -number_of_coefficients + 1
      ca = 1.d0/2.d0*2.d0/3.d0*(0.d0,1.d0)*(2.d0*c_chi + 3.d0*c_chi_z)
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo

      call get_jx2_coefficients(coord(part),cb,nb,part)
      blockstart(13,part,part) = l -number_of_coefficients + 1
      ca = 1.d0/2.d0*a_x
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo

      call get_jy2_coefficients(coord(part),cb,nb,part)
      blockstart(14,part,part) = l -number_of_coefficients + 1
      ca = 1.d0/2.d0*a_y
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo

      call get_jz2_coefficients(coord(part),cb,nb,part)
      blockstart(15,part,part) = l -number_of_coefficients + 1
      ca = 1.d0/(2.d0*9.d0)*(4.d0*c_chi+12.d0*c_chi_z+9.d0*c_z)
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo

      blockstart(16,part,part) = l -number_of_coefficients + 1
      ca = 1.d0/(2.d0*9.d0)*(4.d0*(a_chi-c_chi)+12.d0*(a_chi_z-c_chi_z)
     .   + 9.d0*(a_z-c_z))
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo

      call get_jxjz_coefficients(coord(part),cb,nb,part)
      blockstart(17,part,part) = l -number_of_coefficients + 1
      ca = 1.d0/(6.d0)*(2.d0*b_chi_x+3.d0*b_x_z)
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo

      call get_t_coefficients(coord(part),cb,nb,part)
      blockstart(18,part,part) = l -number_of_coefficients + 1
      ca = 1.d0
      do i = 1, nb
        value = ca*cb(i)
        call set_coeff(value, l, coeff)
      enddo


      blockstart(19,part,part) = l -number_of_coefficients + 1

      coeff_per_part(part) = l - number_of_coefficients 
      number_of_coefficients = l

      end subroutine

c-------------------------------------------------------------
      subroutine set_coefficients_atom_mix(coeff,number_of_coefficients,
     .                                  part)

      implicit none

      complex*16                        coeff(*), value
      integer                           number_of_coefficients, save_nof
      integer                           N, maxN, part, b
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart
      real*8    a_th, a_phi, c_phi, a_chi, c_chi, a_th_phi
      real*8   a_th_y, b_phi_y, b_chi_x, a_chi_z, c_chi_z, a_x
      real*8   a_y, a_z, c_z, b_x_z 
      integer   l, np, nb, ip, ib, i
      complex*16        ca, cb(1:10), cp(1:10)

     

        save_nof = number_of_coefficients
        call set_coefficients_atom(coeff,number_of_coefficients,part)
      
        call get_methyl_coefficients(a_th, a_phi, c_phi, a_chi, c_chi,
     .          a_th_phi)

        call get_atom_coefficients(a_th_y, b_phi_y, b_chi_x,
     .          a_chi_z, c_chi_z, a_x, a_y, a_z, c_z, b_x_z)

      l = number_of_coefficients
      blockstart(13,5,part) = l +1      

      do b= 5, part-1

      call get_jx_coefficients(coord(b),cb,nb,b)
      call get_jx_coefficients(coord(part),cp,np,part)
      blockstart(13,b,part) = l - save_nof + 1
      ca =  1.d0/2.d0*2*a_x
      do ip = 1, np
        do ib = 1, nb
                value = ca*cb(ib)*cp(ip)
                call set_coeff(value, l, coeff)
        enddo
      enddo

      call get_jy_coefficients(coord(b),cb,nb,b)
      call get_jy_coefficients(coord(part),cp,np,part)
      blockstart(14,b,part) = l - save_nof + 1
      ca = 1.d0/2.d0*2.d0*a_y
      do ip = 1, np
        do ib = 1, nb
                value = ca*cb(ib)*cp(ip)
                call set_coeff(value, l, coeff)
        enddo
      enddo

      call get_jz_coefficients(coord(b),cb,nb,b)
      call get_jz_coefficients(coord(part),cp,np,part)
      blockstart(15,b,part) = l - save_nof + 1
      ca = 2.d0*1.d0/(2.d0*9.d0)*(4.d0*c_chi+12.d0*c_chi_z+9.d0*c_z)
      do ip = 1, np
        do ib = 1, nb
                value = ca*cb(ib)*cp(ip)
                call set_coeff(value, l, coeff)
        enddo
      enddo

      blockstart(16,b,part) = l - save_nof + 1
      ca = 2.d0*1.d0/(2.d0*9.d0)*(4.d0*(a_chi-c_chi)
     .              +12.d0*(a_chi_z-c_chi_z)+9.d0*(a_z-c_z))
      do ip = 1, np
        do ib = 1, nb
                value = ca*cb(ib)*cp(ip)
                call set_coeff(value, l, coeff)
        enddo
      enddo

      call get_jx_coefficients(coord(b),cb,nb,b)
      call get_jz_coefficients(coord(part),cp,np,part)
      blockstart(17,b,part) = l - save_nof + 1
      ca = 2.d0/(6.d0)*(2.d0*b_chi_x+3.d0*b_x_z)
      do ip = 1, np
        do ib = 1, nb
                value = ca*cb(ib)*cp(ip)
                call set_coeff(value, l, coeff)
        enddo
      enddo

      blockstart(18,b,part) = l - save_nof + 1
      call get_jz_coefficients(coord(b),cb,nb,b)
      call get_jx_coefficients(coord(part),cp,np,part)
      blockstart(19,b,part) = l - save_nof + 1
      ca = 2.d0/(6.d0)*(2.d0*b_chi_x+3.d0*b_x_z)
      do ip = 1, np
        do ib = 1, nb
                value = ca*cb(ib)*cp(ip)
                call set_coeff(value, l, coeff)
        enddo
      enddo
        
      enddo

      coeff_per_part(part) = l - save_nof
      number_of_coefficients = l

      end subroutine

c-------------------------------------------------------------
        subroutine get_jx_coefficients(type, c, n, part)
        implicit none
        integer         n
        complex*16      c(1:10)
        character*20    type, name
        parameter       (name = 'get_jx_coefficients')
        integer         maxN, part
        parameter       (maxN=20)
        integer         njx(5:maxN), njy(5:maxN), njz(5:maxN)
        common  /nj/    njx, njy, njz      


        if (type.eq.'stereo') then
                n = 3
                c(1) = -1.d0/2.d0*(0.d0,-1.d0)
                c(2) = -1.d0/2.d0*(0.d0,-1.d0)
                c(3) = -1.d0/2.d0*(0.d0,-1.d0)
        else if (type.eq.'spherical') then
                n = 2
                c(1) =  1.d0/2.d0*(0.d0, 1.d0)
                c(2) =  1.d0/2.d0*(0.d0, 1.d0)
        else
                print*, 'ERROR in ', name
                print*, 'coordinates not supported: ', type
                stop
        endif
        
        njx(part) = n

        end subroutine
c------------------------------------------------------------
        subroutine get_jy_coefficients(type, c, n, part)
        implicit none
        integer         n
        complex*16      c(1:10)
        character*20    type, name
        parameter       (name = 'get_jy_coefficients')
        integer         maxN, part
        parameter       (maxN=20)
        integer         njx(5:maxN), njy(5:maxN), njz(5:maxN)
        common  /nj/    njx, njy, njz      


        if (type.eq.'stereo') then
                n = 3
                c(1) = 1.d0/2.d0*(0.d0,-1.d0)
                c(2) = 1.d0/2.d0*(0.d0,-1.d0)
                c(3) = 1.d0/2.d0*(0.d0,-1.d0)
        else if (type.eq.'spherical') then
                n = 2
                c(1) =   1.d0/2.d0*(0.d0,-1.d0)
                c(2) = - 1.d0/2.d0*(0.d0,-1.d0)
        else
                print*, 'ERROR in ', name
                print*, 'coordinates not supported: ', type
                stop
        endif
        
        njy(part) = n

        end subroutine
c------------------------------------------------------------
        subroutine get_jz_coefficients(type, c, n, part)
        implicit none
        integer         n
        complex*16      c(1:10)
        character*20    type, name
        parameter       (name = 'get_jz_coefficients')
        integer         maxN, part
        parameter       (maxN=20)
        integer         njx(5:maxN), njy(5:maxN), njz(5:maxN)
        common  /nj/    njx, njy, njz      


        if (type.eq.'stereo') then
                n = 2
                c(1) = (0.d0,-1.d0)
                c(2) = (0.d0, 1.d0)
        else if (type.eq.'spherical') then
                n = 1
                c(1) = (0.d0,-1.d0)
        else
                print*, 'ERROR in ', name
                print*, 'coordinates not supported: ', type
                stop
        endif
        
        njz(part) = n

        end subroutine
c------------------------------------------------------------
        subroutine get_jx2_coefficients(type, c, n, part)
        implicit none
        integer         n
        complex*16      c(1:10)
        character*20    type, name
        parameter       (name = 'get_jx2_coefficients')
        integer         part


        if (type.eq.'stereo') then
                n = 10
                c( 1) = -1.d0/4.d0 
                c( 2) = -1.d0/4.d0
                c( 3) = -1.d0/4.d0
                c( 4) =  2.d0
                c( 5) =  1.d0/4.d0
                c( 6) = -1.d0/4.d0
                c( 7) = -1.d0/4.d0
                c( 8) =  1.d0/2.d0
                c( 9) =  1.d0
                c(10) =  1.d0/2.d0
        else if (type.eq.'spherical') then
                n = 5
                c(1) =  1.d0
                c(2) =  1.d0
                c(3) = -1.d0/2.d0
                c(4) = -1.d0/2.d0
                c(5) =  1.d0
        else
                print*, 'ERROR in ', name
                print*, 'coordinates not supported: ', type
                stop
        endif
        
        end subroutine
c------------------------------------------------------------
        subroutine get_jy2_coefficients(type, c, n, part)
        implicit none
        integer         n
        complex*16      c(1:10)
        character*20    type, name
        parameter       (name = 'get_jy2_coefficients')
        integer         part


        if (type.eq.'stereo') then
                n = 10
                c( 1) = -1.d0/4.d0 
                c( 2) = -1.d0/4.d0
                c( 3) = -1.d0/4.d0
                c( 4) =  1.d0/2.d0
                c( 5) =  1.d0
                c( 6) =  1.d0/2.d0
                c( 7) = -1.d0/4.d0
                c( 8) = -1.d0/4.d0
                c( 9) =  1.d0/4.d0
                c(10) =  2.d0
        else if (type.eq.'spherical') then
                n = 5
                c(1) =  1.d0
                c(2) =  1.d0
                c(3) =  1.d0/2.d0
                c(4) =  1.d0/2.d0
                c(5) = -1.d0
        else
                print*, 'ERROR in ', name
                print*, 'coordinates not supported: ', type
                stop
        endif
        
        end subroutine
c------------------------------------------------------------
        subroutine get_jz2_coefficients(type, c, n, part)
        implicit none
        integer         n
        complex*16      c(1:10)
        character*20    type, name
        parameter       (name = 'get_jz2_coefficients')
        integer         part


        if (type.eq.'stereo') then
                n = 4
                c(1) = -1.d0/2.d0 
                c(2) =  2.d0
                c(3) =  2.d0
                c(4) =  1.d0/2.d0
        else if (type.eq.'spherical') then
                n = 1
                c(1) =  2.d0
        else
                print*, 'ERROR in ', name
                print*, 'coordinates not supported: ', type
                stop
        endif
        
        end subroutine
c------------------------------------------------------------
        subroutine get_jxjz_coefficients(type, c, n, part)
c       JxJz+JzJx
        implicit none
        integer         n
        complex*16      c(1:10)
        character*20    type, name
        parameter       (name = 'get_jxjz_coefficients')
        integer         part


        if (type.eq.'stereo') then
                n = 7
                c(1) = -1.d0/2.d0 
                c(2) = -2.d0
                c(3) =  2.d0
                c(4) = -2.d0
                c(5) =  3.d0/2.d0
                c(6) = -1.d0/2.d0
                c(7) = -1.d0/2.d0
        else if (type.eq.'spherical') then
                n = 2
                c(1) =  1.d0/2.d0
                c(2) = -2.d0
        else
                print*, 'ERROR in ', name
                print*, 'coordinates not supported: ', type
                stop
        endif
        
        end subroutine
c------------------------------------------------------------
        subroutine get_j2_coefficients(type, c, n, part)
        implicit none
        integer         n
        complex*16      c(1:10)
        character*20    type, name
        parameter       (name = 'get_j2_coefficients')
        integer         part


        if (type.eq.'stereo') then
                n = 8
                c(1) =  1.d0/2.d0 
                c(2) =  1.d0/2.d0
                c(3) =  1.d0
                c(4) =  1.d0
                c(5) =  1.d0/2.d0
                c(6) =  1.d0/2.d0
                c(7) = -1.d0/2.d0
                c(8) = -1.d0/2.d0
        else if (type.eq.'spherical') then
                n = 2
                c(1) =  2.d0
                c(2) =  2.d0
        else
                print*, 'ERROR in ', name
                print*, 'coordinates not supported: ', type
                stop
        endif
        
        end subroutine
c------------------------------------------------------------
        subroutine get_t_coefficients(type, c, n, part)
c       T = -1/2 d^2/dr^2 + 1/2 1/r^2 J^2
        implicit none
        integer         n, i
        complex*16      c(1:10)
        character*20    type, name
        parameter       (name = 'get_j2_coefficients')
        integer         part

        call get_j2_coefficients(type, c, n, part)
        do i = n, 1, -1
         c(i+1) = c(i)*0.5d0
        enddo
        c(1) = 1.d0
        n = n + 1
        
        end subroutine
c-------------------------------------------------------------
        subroutine get_methyl_coefficients(a_th, a_phi, c_phi, a_chi,
     .          c_chi, a_th_phi)
        
       implicit none

       real*8    a_th, a_phi, c_phi, a_chi, c_chi, a_th_phi
       real*8      dth, dph1, dph2 ,dchi

       dth  = 0.0030d0
       dph1 = 0.0045d0
       dph2 = 0.0100d0
       dchi = 0.0300d0

       a_th = 1.d0 + 2.d0*dth + 4.d0/3.d0*dph1 + 1.d0/3.d0*dph2 
     .      + 1.d0/9.d0*dchi
       
       a_phi = 3.d0/2.d0 + 15.d0/8.d0*dth + 3.d0*dph1

       c_phi = 3.d0 + 6.d0*dth + 4.d0*dph1 + 3.d0/2.d0*dph2
     .       + 5.d0/6.d0*dchi

       a_chi = 9.d0/2.d0 + 99.d0/8.d0*dth + 3.d0*dph1

       c_chi = 9.d0 + 18.d0*dth + 12.d0*dph1 + 15.d0/2.d0*dph2
     .       + 6.d0/4.d0*dchi

       a_th_phi = 1.d0/(4.d0*sqrt(3.d0))*(-9.d0*dth + 8.d0*dph1
     .          + 7.d0*dph2 - 7.d0/3.d0*dchi)

        end subroutine
c-------------------------------------------------------------
        subroutine get_atom_coefficients(a_th_y, b_phi_y, b_chi_x,
     .          a_chi_z, c_chi_z, a_x, a_y, a_z, c_z, b_x_z)
        
       implicit none

       real*8   a_th_y, b_phi_y, b_chi_x, a_chi_z, c_chi_z, a_x
       real*8   a_y, a_z, c_z, b_x_z 
       real*8   dth, dph1, dph2 ,dchi

       dth  = 0.0030d0
       dph1 = 0.0045d0
       dph2 = 0.0100d0
       dchi = 0.0300d0

        a_th_y  = 3.d0/2.d0*dth - 4.d0/3.d0*dph1 - 5.d0/6.d0*dph2
     .          + 1.d0/18.d0*dchi

        b_phi_y = sqrt(3.d0)*(-1.d0 - 11.d0/4.d0*dth - 2.d0/3.d0*dph1
     .          - 3.d0/4.d0*dph2 - 1.d0/36.d0*dchi)

        b_chi_x = -3.d0 - 15.d0/4.d0*dth -6.d0*dph1 - 15.d0/4.d0*dph2
     .          - 17.d0/12.d0*dchi

        a_chi_z = -3.d0 - 21.d0/2.d0*dth

        c_chi_z = -6.d0 - 57.d0/4.d0*dth -6.d0*dph1 -15.d0/4.d0*dph2
     .          - 17.d0/12.d0*dchi

        a_x     = 2.d0 + 5.d0/2.d0*dth + 4.d0*dph1 + 2.d0*dph2
     .          + 4.d0/3.d0*dchi

        a_y     = 2.d0 + 11.d0/2.d0*dth + 4.d0/3.d0*dph1
     .          + 4.d0/3.d0*dph2 - 2.d0/9.d0*dchi

        a_z     = 3.d0 + 21.d0/2.d0*dth

        c_z     = 5.d0 + 13.d0*dth + 4.d0*dph1 + 2.d0*dph2 + 4.d0*dchi

        b_x_z   = 2.d0 + 5.d0/2.d0*dth + 4.d0*dph1 + 2.d0*dph2
     .          + 4.d0/3.d0*dchi


        end subroutine
c--------------------------------------------------------------
        subroutine set_coeff(value, index, coeff)
        
        implicit none

        integer index
        complex*16  value, coeff(*)
        

        index = index +1
        coeff(index) = value

        end subroutine


c-------------------------------------------------------------------
        subroutine set_diag_methyl(diag, dim, i)
        
        implicit none

        integer         dim, i
        logical*4       diag(dim,*)
        
        diag(1, 1)=.false.
        diag(1, 2)=.false.
        diag(2, 2)=.false.
        diag(3, 2)=.false.
        diag(1, 3)=.false.
        diag(2, 3)=.false.
        diag(1, 4)=.false.
        diag(2, 4)=.false.
        diag(1, 5)=.false.
        diag(4, 5)=.false.
        
        diag(1, 6)=.false.
        diag(4, 6)=.false.
        diag(1, 7)=.false.
        diag(4, 7)=.false.
        diag(5, 7)=.false.
        diag(1, 8)=.false.
        diag(5, 8)=.false.
        diag(1, 9)=.false.
        diag(4, 9)=.false.
        diag(6, 9)=.false.
        diag(1,10)=.false.
        diag(6,10)=.false.

        diag(1,11)=.false.
        diag(4,11)=.false.
        diag(5,11)=.false.
        diag(1,12)=.false.
        diag(4,12)=.false.
        diag(5,12)=.false.

        i = 12

        end subroutine

c-------------------------------------------------------------------
        subroutine set_diag_atom(diag, dim, i, part)
        
        implicit none

        integer         dim, i, k, term, part, c
        logical*4       diag(dim,*)
        
      integer                           N, maxN, block
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart


        block = 9
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                term = k - blockstart(block, part, part) +1
                i = i + 1
                diag(1,i) = .false.
                diag(4,i) = .false.
                call set_diag_jy(diag, dim, part, term, i)
        enddo

        block = 10
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                i = i + 1
                term = k - blockstart(block, part, part) +1
                diag(1, i) = .false.
                diag(4, i) = .false.
                diag(5, i) = .false.
                call set_diag_jy(diag, dim, part, term, i)
        enddo
        
        block = 11
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                term = k - blockstart(block, part, part) +1
                i = i + 1
                diag(1, i) = .false.
                diag(4, i) = .false.
                diag(6, i) = .false.
                call set_diag_jx(diag, dim, part, term, i)
        enddo

        block = 12
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                term = k - blockstart(block, part, part) +1
                i = i + 1
                diag(1, i) = .false.
                diag(4, i) = .false.
                diag(6, i) = .false.
                call set_diag_jz(diag, dim, part, term, i)
        enddo

        block = 13
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                i = i + 1
                term = k - blockstart(block, part, part) +1
                diag(1, i) = .false.
                call set_diag_jx2(diag, dim, part, term, i)
        enddo

        block = 14
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                i = i + 1
                term = k - blockstart(block, part, part) +1
                diag(1, i) = .false.
                call set_diag_jy2(diag, dim, part, term, i)
        enddo

        block = 15
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                i = i + 1
                term = k - blockstart(block, part, part) +1
                diag(1, i) = .false.
                diag(4, i) = .false.
                call set_diag_jz2(diag, dim, part, term, i)
        enddo

        block = 16
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                i = i + 1
                term = k - blockstart(block, part, part) +1
                diag(1, i) = .false.
                call set_diag_jz2(diag, dim, part, term, i)
        enddo

        block = 17
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                i = i + 1
                term = k - blockstart(block, part, part) +1
                diag(1, i) = .false.
                diag(4, i) = .false.
                call set_diag_jxjz(diag, dim, part, term, i)
        enddo

        block = 18
        do k = blockstart(block, part, part),
     .  blockstart(block+1, part, part)-1
                i = i + 1
                term = k - blockstart(block, part, part) +1
                call set_diag_t(diag, dim, part, term, i)
        enddo


        end subroutine
c-------------------------------------------------------------------
        subroutine set_diag_atom_mix(diag, dim, i, part)
        
        implicit none

        integer         dim, i, kb, kp, b, part
        logical*4       diag(dim,*)
        
      integer                           N, maxN, block
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart
        integer         njx(5:maxN), njy(5:maxN), njz(5:maxN)
        common  /nj/    njx, njy, njz      

        call set_diag_atom(diag,dim,i,part)


        do b = 5, part-1

        block = 13
        do kp = 1, njx(part)
                do kb = 1, njx(b) 
                        i = i + 1
                        diag(1, i) = .false.
                        call set_diag_jx(diag, dim, b, kb, i)
                        call set_diag_jx(diag, dim, part, kp, i)
                enddo
        enddo
        

        block = 14
        do kp = 1, njy(part)
                do kb = 1, njy(b) 
                        i = i + 1
                        diag(1, i) = .false.
                        call set_diag_jy(diag, dim, b, kb, i)
                        call set_diag_jy(diag, dim, part, kp, i)
                enddo
        enddo
        

        block = 15
        do kp = 1, njz(part)
                do kb = 1, njz(b) 
                        i = i + 1
                        diag(1, i) = .false.
                        diag(4, i) = .false.
                        call set_diag_jz(diag, dim, b, kb, i)
                        call set_diag_jz(diag, dim, part, kp, i)
                enddo
        enddo

        block = 16
        do kp = 1, njz(part)
                do kb = 1, njz(b) 
                        i = i + 1
                        diag(1, i) = .false.
                        call set_diag_jz(diag, dim, b, kb, i)
                        call set_diag_jz(diag, dim, part, kp, i)
                enddo
        enddo
        

        block = 17
        do kp = 1, njz(part)
                do kb = 1, njx(b) 
                        i = i + 1
                        diag(1, i) = .false.
                        diag(4, i) = .false.
                        call set_diag_jx(diag, dim, b, kb, i)
                        call set_diag_jz(diag, dim, part, kp, i)
                enddo
        enddo

        block = 19
        do kp = 1, njx(part)
                do kb = 1, njz(b) 
                        i = i + 1
                        diag(1, i) = .false.
                        diag(4, i) = .false.
                        call set_diag_jz(diag, dim, b, kb, i)
                        call set_diag_jx(diag, dim, part, kp, i)
                enddo
        enddo

        enddo

        end subroutine
c------------------------------------------------------------------ 
        subroutine set_diag_jx(diag, dim, part, term, k)
        
        implicit none

        integer         dim, part, term, k
        logical*4       diag(dim,*)
        character*20    name
        parameter       (name='set_diag_jx')
        
      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

c     hannes fix
      coord(part) = 'stereo'

        if (coord(part).eq.'stereo') then
                SELECT CASE (term)
                       CASE(1:2)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE(3)
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else if (coord(part).eq.'spherical') then
                SELECT CASE (term)
                       CASE(1:2)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                                 


        end subroutine
c------------------------------------------------------------------ 
        subroutine set_diag_jy(diag, dim, part, term, k)
        
        implicit none

        integer         dim, part, term, k
        logical*4       diag(dim,*)
        character*20    name
        parameter       (name='set_diag_jy')
        
      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

        if (coord(part).eq.'stereo') then
                SELECT CASE (term)
                       CASE(1:2)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE(3)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else if (coord(part).eq.'spherical') then
                SELECT CASE (term)
                       CASE(1:2)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                                 


        end subroutine
c------------------------------------------------------------------ 
        subroutine set_diag_jz(diag, dim, part, term, k)
        
        implicit none

        integer         dim, part, term, k
        logical*4       diag(dim,*)
        character*30    name
        parameter       (name='set_diag_jz')
        
      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart
c     hannes fix
      coord(part) = 'stereo'

        if (coord(part).eq.'stereo') then
                SELECT CASE (term)
                       CASE(1:2)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else if (coord(part).eq.'spherical') then
                SELECT CASE (term)
                       CASE(1)
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                                 


        end subroutine
c------------------------------------------------------------------ 
        subroutine set_diag_jx2(diag, dim, part, term, k)
        
        implicit none

        integer         dim, part, term, k
        logical*4       diag(dim,*)
        character*30    name
        parameter       (name='set_diag_jx2')
        
      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

        if (coord(part).eq.'stereo') then
                SELECT CASE (term)
                       CASE(1)
                                return
                       CASE(2)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE(3)
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE(4:9)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE(10)
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else if (coord(part).eq.'spherical') then
                SELECT CASE (term)
                       CASE(1)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE(2:5)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif
                                 
        end subroutine
c------------------------------------------------------------------ 
        subroutine set_diag_jy2(diag, dim, part, term, k)
        
        implicit none

        integer         dim, part, term, k
        logical*4       diag(dim,*)
        character*20    name
        parameter       (name='set_diag_jy2')
        
      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

        if (coord(part).eq.'stereo') then
                SELECT CASE (term)
                       CASE(1)
                                return
                       CASE(2)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE(3)
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE(4)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE(5:10)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else if (coord(part).eq.'spherical') then
                SELECT CASE (term)
                       CASE(1)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE(2:5)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif

        end subroutine
c------------------------------------------------------------------ 
        subroutine set_diag_jz2(diag, dim, part, term, k)
        
        implicit none

        integer         dim, part, term, k
        logical*4       diag(dim,*)
        character*20    name
        parameter       (name='set_diag_jz2')
        
      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

        if (coord(part).eq.'stereo') then
                SELECT CASE (term)
                       CASE(1)
                                return
                       CASE(2:4)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else if (coord(part).eq.'spherical') then
                SELECT CASE (term)
                       CASE(1)
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif

        end subroutine
c------------------------------------------------------------------ 
        subroutine set_diag_jxjz(diag, dim, part, term, k)
        
        implicit none

        integer         dim, part, term, k
        logical*4       diag(dim,*)
        character*20    name
        parameter       (name='set_diag_jxjz')
        
      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

        if (coord(part).eq.'stereo') then
                SELECT CASE (term)
                       CASE(1)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE(2:7)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else if (coord(part).eq.'spherical') then
                SELECT CASE (term)
                       CASE(1:2)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif

        end subroutine
c------------------------------------------------------------------ 
        subroutine set_diag_j2(diag, dim, part, term, k)
        
        implicit none

        integer         dim, part, term, k
        logical*4       diag(dim,*)
        character*20    name
        parameter       (name='set_diag_j2')
        
      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

        if (coord(part).eq.'stereo') then
                SELECT CASE (term)
                       CASE(1:4)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE(5)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE(6)
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE(7)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE(8)
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else if (coord(part).eq.'spherical') then
                SELECT CASE (term)
                       CASE(1)
                                diag(6+(part-5)*3+2,k) = .false.
                       CASE(2)
                                diag(6+(part-5)*3+2,k) = .false.
                                diag(6+(part-5)*3+3,k) = .false.
                       CASE DEFAULT
                       call print_error(name,'term unknown', term, part)
                END SELECT
        else
            print*, 'ERROR in ', name
            print*, 'type not supported: ', coord(part)
            stop
        endif

        end subroutine
c------------------------------------------------------------------ 
        subroutine set_diag_t(diag, dim, part, term, k)
        
        implicit none

        integer         dim, part, term, k
        logical*4       diag(dim,*)
        character*20    name
        parameter       (name='set_diag_t')
        
      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      integer                           blockstart(9:19,5:maxN,5:maxN)
      character*20                      coord(5:maxN)
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

        diag(6+(part-5)*3+1,k) = .false.
        if (term.gt.1) then
                term = term -1
                call set_diag_j2(diag, dim, part, term, k)
        endif

        end subroutine

      subroutine koordreadschaepers(mass,na)
      implicit none

      integer    na,na3,dimch3

      real*8     mch3,hmass,mass(20)
      parameter (hmass=1836.109d0)

      integer i,j

      real*8 work(4,4),m(20),alpha(4)
      common /rctrans/ work,m
      SAVE /rctrans/

      na3 = 3 * na
      dimch3 = 12

c masses
      write(*,*) 'mass of non-reactive H is ', mass(1)
      m(1)=mass(2)/hmass
      m(2)=mass(1)/hmass
      m(3)=mass(1)/hmass
      m(4)=mass(1)/hmass

      do i = 3,na
        m(2+i) = mass(i)/hmass
      enddo

c      don't touch stuff below this line
c      generic part starts here
      do i=1,na
         m(i)=m(i)*hmass
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

      subroutine IntToCart(q,x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine transforms the internal coordinates       c
c to cartesian ones.                                        c
c Subroutine koordread must be called in main program       c
c Don't change anything below this line unless you          c
c know what you do. This part is written to work with the   c
c generic kinetic energy operator                           c
c DS / 4 May 2018                                           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer    dimch3
      parameter  (dimch3=12)

      real*8  q(60),x(60)
      real*8  r(3),xcm(3), snew, tnew
      real*8  radau(dimch3)
      real*8  mch3,mch4
      real*8     hmass, pi
      parameter (hmass=1836.109d0)
      parameter (pi=3.14159265359d0)

      real*8   rmch3
      real*8   dist

      integer i,j,k

      integer                           N, maxN
      parameter                         (maxN=20)
      integer                           coeff_per_part(4:maxN)
      character*20                      coord(5:maxN)
      integer                           blockstart(9:20,5:maxN,5:maxN)

      real*8 rctrans(4,4),m(maxN)
      common /rctrans/ rctrans,m
      common    /coordinate_choice/     N, coord, coeff_per_part,
     .                                  blockstart

      real*8    redm(maxN), mabs(maxN), sredm(maxN)
      integer   start

c mass of CH3-group
      mch3=0.0d0
      do i=1,4
         mch3=mch3+m(i)
      end do
      mabs(4) = mch3

      do i = 5, N
        mabs(i) = 0.d0
        mabs(i) = mabs(i-1) + m(i)
        redm(i) = mabs(i-1)*m(i)/(mabs(i))
        sredm(i) = 1.d0/sqrt(redm(i))
      enddo


c order of q's:
c  1 rho_2        4 theta     7 R           10 r
c  2 phi_2        5 phi       8 theta_R     11 theta_r
c  3 theta_2      6 chi       9 phi_R       12 phi_r

c order of radau's
c  1 center of mass    4 x H_1     7 x H_2   10 x H_3
c  2 "                 5 y H_1     8 y H_2   11 y H_3
c  3 "                 6 z H_1     9 z H_2   12 z H_3

c order of cartesians
c 1-3 C (x,y,z)
c 4-6, 7-9, 10-12 methyl-H
c 13- 3*N coordinates for additional atoms


c METHYL

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

      !print*, r,q

c     print*,"Radau in Ruecktransformation"
      !do i=1,4
      !  print'(3F12.6)',(radau(j+(i-1)*3),j=1,3)
      !end do

c Third step: calculate cartesians for methyl group
c             from Radau's
      do i=1,dimch3,3
         x(i  )=0.0d0
         x(i+1)=0.0d0
         x(i+2)=0.0d0
         do j=1,dimch3,3
!            print*, 'koord', i, j
!            print*, i/3, j/3
!            print*, rctrans(i/3+1,j/3+1)
           x(i  )=x(i  )+rctrans(i/3+1,j/3+1)*radau(j  )
           x(i+1)=x(i+1)+rctrans(i/3+1,j/3+1)*radau(j+1)
           x(i+2)=x(i+2)+rctrans(i/3+1,j/3+1)*radau(j+2)
         end do
      end do
        !print*, 'x',x
c  unmassweight cartesians, correct sign
      do i=1,dimch3,3
         x(i  )=-x(i  )/sqrt(m(i/3+1))
         x(i+1)=-x(i+1)/sqrt(m(i/3+1))
         x(i+2)=-x(i+2)/sqrt(m(i/3+1))
      end do

c ADDITIONAL ATOMS 5 to N

      do i = 5, N

        start = 6+(i-5)*3
c     unmassweight r
      dist=q(start+1)*sredm(i)

      if (coord(i).eq.'stereo') then
c     inttocart for STEREOGRAPHIC COORDINATES (r, s, t)

        snew = 0.5d0*(q(start+2)*(-cos(2.d0/3.d0*q(6))
     .                      +sqrt(3.d0)*sin(2.d0/3.d0*q(6)))
     .               +q(start+3)*( sin(2.d0/3.d0*q(6))
     .                      +sqrt(3.d0)*cos(2.d0/3.d0*q(6))))


        tnew = 0.5d0*(q(start+3)*(-cos(2.d0/3.d0*q(6))
     .                      +sqrt(3.d0)*sin(2.d0/3.d0*q(6)))
     .               -q(start+2)*( sin(2.d0/3.d0*q(6))
     .                      +sqrt(3.d0)*cos(2.d0/3.d0*q(6))))

c     Center of mass for atoms 1 to i-1
        xcm(1)=0.0d0
        xcm(2)=0.0d0
        xcm(3)=0.0d0
        do k=1, i -1
                xcm(1)=xcm(1)+x(3*(k-1)+1)*m(k)
                xcm(2)=xcm(2)+x(3*(k-1)+2)*m(k)
                xcm(3)=xcm(3)+x(3*(k-1)+3)*m(k)
        end do
        do k=1,3
                xcm(k)=xcm(k)/(mabs(i-1))
        end do

c  calculate cartesians
        x((i-1)*3+1)=xcm(1)+2*dist*snew/(1.d0+snew**2+tnew**2)
        x((i-1)*3+2)=xcm(2)+2*dist*tnew/(1.d0+snew**2+tnew**2)
        x((i-1)*3+3)=xcm(3)
     .       +dist*(1.d0-snew**2-tnew**2)/(1.d0+snew**2+tnew**2)

      else if (coord(i).eq.'spherical') then
c     inttocart for SPHERICAL COORDINATES (r, theta, phi)

c    Center of mass for atoms 1 to i-1
        xcm(1)=0.0d0
        xcm(2)=0.0d0
        xcm(3)=0.0d0
        do k=1, i -1
                xcm(1)=xcm(1)+x(3*(k-1)+1)*m(k)
                xcm(2)=xcm(2)+x(3*(k-1)+2)*m(k)
                xcm(3)=xcm(3)+x(3*(k-1)+3)*m(k)
        end do
        do k=1,3
                xcm(k)=xcm(k)/(mabs(i-1))
        end do

        x((i-1)*3+1) = xcm(1) +
     .          dist*sin(q(start+2))*cos(q(start+3)+2.d0/3.d0*(q(6)-pi))
        x((i-1)*3+2) = xcm(2) +
     .          dist*sin(q(start+2))*sin(q(start+3)+2.d0/3.d0*(q(6)-pi))
        x((i-1)*3+3) = xcm(3) + dist*cos(q(start+2))

      else
                print*, 'ERROR: not supported: ', coord(i)
                print*, 'coordinates for atom:', i
                stop
      endif
      enddo

      end