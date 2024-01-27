!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example IM1, see the previous cited paper for the original reference.
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
!# x variable
!var x{1..2};
!
!minimize f1:
!    2*sqrt(x[1]);
!minimize f2:
!    x[1]*(1-x[2])+5;
!
!subject to bound1:
!    1 <= x[1] <= 4;
!subject to bound2:
!    1 <= x[2] <= 2;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 2
	m = 0
	q = 2

	return
end subroutine setdim

subroutine startp(n,x)
	implicit none
	integer	:: n
	real*8		:: x(n), l(n), u(n)

	call setbounds(n,l,u)

	x = (l+u)/2.d0
	
	return
end subroutine startp

subroutine functs(n,x,M,f)
	implicit none
	integer	:: n, M
	real*8		:: x(n), f(M)

	f(1) = 2.d0*sqrt(x(1))
	f(2) = x(1)*(1.d0-x(2))+5.d0

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb(1) = 1.d0; 	    ub(1) = 4.d0
	lb(2) = 1.d0; 	    ub(2) = 2.d0

	return
end subroutine setbounds
