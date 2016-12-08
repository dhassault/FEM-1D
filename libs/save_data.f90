module save_data
  implicit none

  contains

!************************************************************

subroutine read_sol(u,alpha,n,p,yk,usol)
  use functions
  use linear_system
!******************************************
! Put the solution in a file named sol.dat
! Analytical and fem solutions are written
!   --inputs--
! u: the solution of Au=b
! alpha: set by the user
! n: set by the user
! p: the number of point of the solution to display
! yk: vector of points where usol is calculated
! usol: solution of the linear system
!   --outputs--
! sol.dat is created in the same directory
!******************************************
  integer, intent(in) :: n,p
  integer :: i,j
  double precision, intent(in) :: alpha
  double precision, dimension(p+1), intent(in) :: yk
  double precision, dimension(n), intent(in) :: u
  double precision, dimension(p) :: usol

  usol(:) = 0.d0
  do j = 1_8,p
     do i = 1_8,n
        usol(j) = usol(j) + u(i)*fbase(i,n,yk(j))
     end do
  end do

  if (p<20) then
     write (6,*) 'U (solution):'
     do j = 1_8,p
        write (6,'(f12.8)') usol(j)
     end do
     write (6,'(f12.8)') 0.d0
     write (6,*) ''
  end if

  open(unit = 11,file="sol.dat", form="formatted")
  do j = 1_8,p
     write (11,*) yk(j), usol(j), sol(yk(j),alpha)
  end do
  write (11,*) 1.d0, 0.d0, 0.d0
  close(11)

end subroutine read_sol

!************************************************************
end module save_data
