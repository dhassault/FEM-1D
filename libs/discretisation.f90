module discretisation
  implicit none

  contains
    
!************************************************************

subroutine mesh(n,xk)
!******************************************
! Discretisation of ]0,1[ in n elements.
!   --inputs--
! n: set by the user
!   --outputs--
! xk: vector of n+1 elements
!******************************************
  integer :: i,n
  double precision :: h
  double precision, dimension(*) :: xk

  h = 1.d0/n
  do i = 0_8,n
     xk(i+1) = i*h
  end do
end subroutine mesh

!************************************************************
end module discretisation
