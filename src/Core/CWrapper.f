
      subroutine openfortran(num, filename)
      implicit none
      integer	num
      character*200 filename
      
      open(num, file=filename, form='unformatted')

      end
      
      subroutine freadcomplexfortran(id, array, dim)
      implicit none
      integer	id, dim, i
	complex*16  array(dim)
      
      do i=1, dim
      read(id) array(i)
      enddo

      end
      
      
