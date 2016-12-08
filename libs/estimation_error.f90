module estimation_error
  implicit none

  contains

!************************************************************

subroutine error(usol,alpha,yk,p,err,error_max)
  use functions
  use linear_system
!********************************************************
! Put error in error.dat file and display the maximal error
!   --inputs--
! usol: solution of the linear system
! alpha: set by the user
! yk:
! p: the number of point of the solution to display
!   --outputs--
! err:
! error_max:
! error.dat is created in the same directory
!*********************************************************
  integer, intent(in) :: p
  integer :: j
  double precision, intent(in) :: alpha
  double precision :: error_max
  double precision, dimension(p+1), intent(in) :: yk
  double precision, dimension(p+1) :: err
  double precision, dimension(p) :: usol
!
  error_max = 0.d0
  do j = 1_8,p
     err(j) = abs(usol(j) - sol(yk(j),alpha))
     if (error_max < err(j)) then
       error_max = err(j)
     end if
  end do
  err(p+1) = 0.d0
!
  open(unit=11, file="error.dat", form="formatted")
  do j = 1_8,p+1
     write (11,*) yk(j), err(j)
  end do
  close(11)
!
end subroutine error

!************************************************************
end module estimation_error
