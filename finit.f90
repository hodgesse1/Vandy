subroutine finit(ierror)
  implicit none
 include 'mpif.h'
  
  integer :: ierror

  call MPI_INIT(ierror)

end subroutine finit

