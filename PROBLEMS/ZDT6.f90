
! ###############################################################################
! #
! #   As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
! #   Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary 
! #   Computation 8(2): 173-195, 2000.
! #
! #   Example T6.
! #
! #   This file is part of a collection of problems developed for
! #   derivative-free multiobjective optimization in
! #   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
! #   Direct Multisearch for Multiobjective Optimization, 2010.
! #
! #   Written by the authors in June 1, 2010.
! #
! ###############################################################################
!
! # parameters
! param m := 10, >=2;
! param pi := 4*atan(1);
!

! # x variable
! var x{1..m};

! # functions
! var ff1 = 1-exp(-4*x[1])*sin(6*pi*x[1])^6;

! # g(x)
! var gx = 1 + 9 * ((sum {i in 2..m} (x[i]))/(m-1))^0.25;

! var h = 1-(ff1/gx)^2;

! minimize f1:
!	ff1;
! minimize f2:
!    gx*h;
!
! subject to bounds {i in 1..m}:
!    0.0 <= x[i] <= 1.0;




subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 10
	m = 0
	q = 2

	return
end subroutine setdim

subroutine startp(n,x)
	implicit none
	integer	:: n
	real*8		:: x(n), l(n), u(n)

	!call setbounds(n,l,u)

	x = 0.d0
	
	return
end subroutine startp

subroutine functs(n,x,M,f)
	implicit none
	integer	:: n, M, i
	real*8		:: x(n), f(M), h
    real*8     :: pi



 pi = 4.d0*atan(1.d0)


! # functions
! var ff1 = 1-exp(-4*x[1])*sin(6*pi*x[1])^6;

! # g(x)
! var gx = 1 + 9 * ((sum {i in 2..m} (x[i]))/(m-1))^0.25;

! var h = 1-(ff1/gx)^2;

! minimize f1:
!	ff1;
! minimize f2:
!    gx*h;

f=0.0d0

f(1)= 1- dexp(-4.0d0*x(1))*dsin(6.0d0*pi*x(1))**6

f(2)= 0.0d0 

do i=2, n
   f(2)= f(2)+ x(i)
end do
   
f(2)= ( f(2)/dble(n-1) )**0.25d0

f(2) = 9.0d0 * f(2)

f(2) = f(2) + 1.0d0 

h= 1.0d0 - (f(1)/f(2))**2

f(2) = f(2)*h

return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.0d0;	    ub = 1.0d0

	return
end subroutine setbounds

