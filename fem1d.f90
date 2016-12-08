program fem1d
  use linear_system
  use discretisation
  use solver
  use save_data
  use estimation_error
  implicit none

  integer :: n, p
  double precision :: alpha, error_max

  double precision, allocatable :: b(:), err(:), u(:), usol(:), xk(:), yk(:)
  double precision, allocatable :: a(:,:)

  write (6,*) ''
  write (6,*) '1D Finite element method'
  write (6,*) ''
  write (6,*) '-u'''' + u = f     in ]0,1['
  write (6,*) 'u''(0) = alpha'
  write (6,*) 'u(1) = 0'
  write (6,*) ''
  write (6,'(a,$)') 'Enter alpha: '
  read (5,*) alpha
  write (6,'(a,$)') 'Enter the number of elements n: '
  read (5,*) n

  allocate(xk(n+1))
  allocate(a(n,n))
  allocate(b(n))
  allocate(u(n))

! We are first shaping the linear system and we are resolving it

  call mesh(n,xk)             !creating the mesh from the space ]0,1[
  call calAB(alpha,n,xk,a,b)  !shaping of the matrix A and vector b
  call solve(a,b,n,u)         !Resolving the linear system Au=b

! Then we put the solution on a file sol.dat

  write (6,*) ''
  write(6,'(a,$)') 'Enter the number of point of the solution to show: '
  read (5,*) p
  p = p-1
  p = p/2*2
  write (6,*) ''

  allocate(yk(p+1))
  allocate(usol(p))

  call mesh(p,yk)
  call read_sol(u,alpha,n,p,yk,usol)

! Now we can estimate the error

  allocate(err(p))

  call error(usol,alpha,yk,p,err,error_max)
  write (6,*) 'error_max = ', error_max

! Job done, we deallocate everything

  deallocate(a)
  deallocate(b)
  deallocate(err)
  deallocate(u)
  deallocate(usol)
  deallocate(xk)
  deallocate(yk)


end program fem1d
