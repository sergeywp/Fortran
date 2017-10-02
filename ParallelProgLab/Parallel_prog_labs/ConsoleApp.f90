      program prog
	  use omp_lib

      implicit none

	  integer,parameter:: count=100
      real*8 :: e=1,fact=1,t1,t2
	  integer :: i
	  t1=omp_get_wtime()
	  !$omp parallel reduction(+:e)
	  !$omp do
	  do i=1,count/4
		e=1/fact
	  enddo
	  !$omp end parallel
	  t2=omp_get_wtime()
	  write (*,*) "e: ",e, "Time: ",t2-t1
	  read(*,*)
      end program prog
