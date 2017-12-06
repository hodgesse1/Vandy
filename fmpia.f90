subroutine fmpia(comm,size)
  implicit none

  include 'mpif.h'



  integer rank,size,ierror,tag,status(MPI_STATUS_SIZE),i,np
  integer comm

  ierror = 0
  
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierror)

  
    
end subroutine fmpia





  
  
