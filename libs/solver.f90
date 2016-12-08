module solver
  implicit none

  contains
!************************************************************

subroutine solve(a,b,n,u)
!******************************************
! Resolving Au=b with Thomas algorithm
! A have to be tridiagonal and positive
!   --inputs--
! a: matrix A
! b: matrix b
! n: set by the user
!   --outputs--
! u: the solution of Au=b
!******************************************
  integer, intent(in) :: n
  integer :: i
  double precision, dimension(n) :: b,c,u
  double precision, dimension(n,n) :: a

  c(1) = 0.d0
  do i = 2_8,n
     c(i) = a(i,i-1)/a(i-1,i-1)
     a(i,i) = a(i,i) - c(i)*a(i-1,i)
     b(i) = b(i) - c(i)*b(i-1)
  end do

  u(n) = b(n)/a(n,n)
  do i = n-1,1_8,-1
     u(i) = (b(i) - a(i,i+1)*u(i+1))/a(i,i)
  end do

  if (n < 20) then
     write(6,*) 'U'
     do i = 1_8,n
        write(6,'(f12.8)') u(i)
     end do
     write (6,*) ''
  end if

end subroutine solve

!************************************************************
end module solver
