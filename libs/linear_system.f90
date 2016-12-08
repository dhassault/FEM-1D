module linear_system
  implicit none

  contains
!************************************************************

subroutine integ(K,l1,l2,n,xk,s)
  use functions
!******************************************
! Calculating the integral s on the K element of phi_l1*phi_l2
! K: K element (global)
! l1: local numbering of the k element
! l2: local numbering of the k element
! n: set by the user
! xk: vector of n+1 elements
! s: the value of the integral
!******************************************
  integer, intent(in) :: K,l1,l2,n
  integer :: i,j
  double precision, parameter :: cm = (0.5d0*(1.d0 - 1.d0/sqrt(3.d0))) , cp =  (0.5d0*(1.d0 + 1.d0/sqrt(3.d0)))
  double precision :: h,s,y11,y12,y21,y22
  double precision, dimension(*), intent(in) :: xk
!
  h = 1.d0/n
  i = local2global(K,l1)
  j = local2global(K,l2)

  y11 = fbase(i,n,xk(K)+cm*h) ! y11 = phi_l1(nm)
  y12 = fbase(i,n,xk(K)+cp*h) ! y12 = phi_l1(np)
  y21 = fbase(j,n,xk(K)+cm*h) ! y21 = phi_l2(nm)
  y22 = fbase(j,n,xk(K)+cp*h) ! y22 = phi_l2(np)

  s = (1/(2.d0*n))*(y11*y21 + y12*y22)

end subroutine integ

subroutine integd(K,l1,l2,n,ds)
  use functions
!******************************************
! Calculating the integral ds on the K element of dphi_l1*dphi_l2
! K: K element (global)
! l1: local numbering of the k element
! l2: local numbering of the k element
! n: set by the user
! ds: the value of the integral
!******************************************
  integer, intent(in) :: K,l1,l2,n
  integer :: i,j
  double precision :: ds, h, y1,y2

  h = 1.d0/n
  i = local2global(K,l1)
  j = local2global(K,l2)
  y1 = dfbase(i,n,(K-1)*h)
  y2 = dfbase(j,n,(K-1)*h)
  ds = h*y1*y2

end subroutine integd

subroutine integf(K,l,n,xk,sf)
  use functions
!******************************************
! Calculating the integral sf on the K element of f*dphi_l
! K: K element (global)
! l: local numbering of the k element
! n: set by the user
! xk: vector of n+1 elements
! sf: the value of the integral
!******************************************
  integer, intent(in) :: K,l,n
  integer :: i
  double precision, parameter :: cm = (0.5d0*(1.d0 - 1.d0/sqrt(3.d0))) , cp =  (0.5d0*(1.d0 + 1.d0/sqrt(3.d0)))
  double precision :: sf,h, y11,y12,y21,y22
  double precision, dimension(*), intent(in) :: xk

  h = 1.d0/n
  i = local2global(K,l)

  y11 = f(xk(K)+cm*h)       ! y11 = f(n1)
  y12 = f(xk(K)+cp*h)       ! y12 = f(n2)
  y21 = fbase(i,n,xk(K)+cm*h) ! y21 = phi_l(n1)
  y22 = fbase(i,n,xk(K)+cp*h) ! y22 = phi_l(n2)

  sf = 0.5d0*h*(y11*y21+y12*y22)

end subroutine integf

subroutine calel(K,n,xk,ak,bk)
!******************************************
! Calculating the elementary matrix A_K and the second one b_K of the linear system
! K: K element (global)
! xk: vector of n+1 elements
! ak: matrix A_K
! bk: matrix b_K
!******************************************
  integer, intent(in) :: K,n
  integer :: l,l1,l2
  double precision :: s,ds,sf
  double precision, dimension(2) :: bk
  double precision, dimension(2,2) :: ak
  double precision, dimension(*) :: xk

  do l1 = 1_8,2_8
     do l2 = 1_8,2_8
        call integ(K,l1,l2,n,xk,s)
        call integd(K,l1,l2,n,ds)
        ak(l1,l2) = s + ds
     end do
  end do

  do l = 1_8,2_8
     call integf(K,l,n,xk,sf)
     bk(l) = sf
  end do

end subroutine calel

subroutine calAB(alpha,n,xk,a,b)
  use functions
!******************************************
! Gathering the matrix A and the vector b
! alpha: set by the user
! n: set by the user
! xk: vector of n+1 elements
! a: matrix A
! b: matrix b
!******************************************
  integer, intent(in) :: n
  integer :: i,j,K,l,l1,l2
  double precision, intent(in) :: alpha
  double precision, dimension(*) :: xk
  double precision, dimension(n) :: b
  double precision, dimension(n,n) :: a
  double precision, dimension(2) :: bk
  double precision, dimension(2,2) :: ak

  a(:,:) = 0.d0
  b(:) = 0.d0
  b(1) = -alpha

  do K = 1_8,n-1
     call calel(K,n,xk,ak,bk)

     do l1 = 1_8,2_8
        do l2 = 1_8,2_8
           i = local2global(K,l1)
           j = local2global(K,l2)
           a(i,j) = a(i,j) + ak(l1,l2)
        end do
     end do

     do l = 1_8,2_8
        i = local2global(K,l)
        b(i) = b(i) + bk(l)
     end do
  end do

  call calel(n,n,xk,ak,bk)
  a(n,n) = a(n,n) + ak(1,1)
  b(n) = b(n) + bk(1)

  if (n < 20) then
     write (6,*) 'A'
     do i = 1_8,n
        do j = 1_8,n
           write(6,'(f12.5,$)') a(i,j)
        end do
        write (6,*) ''
     end do
     write (6,*) ''

     write(6,*) 'B'
     do i = 1_8,n
        write(6,'(f12.5)') b(i)
     end do
     write (6,*) ''
  end if

end subroutine calAB

!***********************************************************
end module linear_system
