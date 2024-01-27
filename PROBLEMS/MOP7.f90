!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example MOP7, Van Valedhuizen's test suit.
!#
!#   This file is part of a collection of problems developed for
!#   derivative-free multiobjective optimization in
!#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
!#   Direct Multisearch for Multiobjective Optimization, 2010.
!#
!#   Written by the authors in June 1, 2010.
!#
!###############################################################################
!
!# Define objective function variables
!var x{1..2} >= -400, <=400;
!
!# The objective functions
!minimize fobj1:
!    (x[1]-2)^2/2+(x[2]+1)^2/13+3;
!
!minimize fobj2:
!    (x[1]+x[2]-3)^2/36+(-x[1]+x[2]+2)^2/8-17;
!
!minimize fobj3:
!    (x[1]+2*x[2]-1)^2/175+(-x[1]+2*x[2])^2/17-13;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 2
	m = 0
	q = 3

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
	real*8		:: x(n), f(M)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

!minimize fobj1:
!    (x[1]-2)^2/2+(x[2]+1)^2/13+3;
	f(1) = (x(1)-2.d0)**2.d0/2.d0+(x(2)+1.d0)**2.d0/13.d0+3.d0
!minimize fobj2:
!    (x[1]+x[2]-3)^2/36+(-x[1]+x[2]+2)^2/8-17;
	f(2) = (x(1)+x(2)-3.d0)**2.d0/36.d0+(-x(1)+x(2)+2.d0)**2.d0/8.d0-17.d0
!minimize fobj3:
!    (x[1]+2*x[2]-1)^2/175+(-x[1]+2*x[2])^2/17-13;
	f(3) = (x(1)+2.d0*x(2)-1.d0)**2.d0/175.d0+(-x(1)+2.d0*x(2))**2.d0/17.d0-13.d0

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -400.d0; 	    ub = 400.d0

	return
end subroutine setbounds
