module functions
  implicit none

  contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function f(x)
  double precision, intent(in) :: x
  double precision :: f
  f = x
end function f

function sol(x,alpha)
! Fonction solution pour f = x.
!******************************************
! solutions of f=x
! x: the point where f is calculated
! alpha: set by the user
!******************************************
  double precision, intent(in) :: x,alpha
  double precision :: beta, sol
  beta = alpha -1.d0
  sol = -((beta*sinh(1.d0)+1.d0)/cosh(1.d0))*cosh(x) + beta*sinh(x) + x
end function sol

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function fbase(i,n,x)
!******************************************
! For phi_i functions
! n: set by the user
! i: local2global(K,l1) with l1 the local numbering
! x: the point where phi_i is calculated
! phi: value of phi_i
!******************************************
  integer :: i,n
  double precision :: h, fbase, x

  h = 1.d0/n
  if ((x > (i-2.d0)*h) .AND. (x < (i-1.d0)*h)) then
     fbase = n*x - i + 2.d0
  else if ((x >= (i-1.d0)*h) .AND. (x < i*h)) then
     fbase = -n*x + i
  else
     fbase = 0.d0
  end if

end function fbase

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function dfbase(i,n,x)
!******************************************
! Differentiation of phi_i previously calculated
! n: set by the user
! i: local2global(K,l1) with l1 the local numbering
! x: the point where phi_i is calculated
! phi: value of phi_i
!******************************************
  integer, intent(in) :: i,n
  double precision :: h, dfbase, x
  h = 1.d0/n
  if ((x >= (i-2.d0)*h) .AND. (x < (i-1.d0)*h)) then
     dfbase = n
  else if ((x >= (i-1.d0)*h) .AND. (x <= i*h)) then
     dfbase = -n
  else
     dfbase = 0.d0
  end if

end function dfbase

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function local2global(K,l)
!******************************************
! local to global numbering
! K: K element (global)
! l: local numbering
!******************************************
  integer, intent(in) :: K,l
  integer :: local2global
  local2global = K + l - 1

end function local2global

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module functions
