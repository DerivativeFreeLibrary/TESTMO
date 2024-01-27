!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example TKLY1, see the previous cited paper for the original reference.
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
!var x{1..4};
!
!minimize f1:
!    x[1];
!minimize f2:
!    (prod {i in 2..4} (2.0-exp(-((x[i]-0.1)/0.004)^2)-0.8*exp(-((x[i]-0.9)/0.4)^2)))/x[1];
!    
!subject to bound:
!    0.1<=x[1]<=1;
!
!subject to bounds{i in 2..4}:
!    0<=x[i]<=1;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 4
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
	integer	:: n, M, i
	real*8		:: x(n), f(M)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

!minimize f1:
!    x[1];
	f(1) = x(1)
!minimize f2:
!    (prod {i in 2..4} (2.0-exp(-((x[i]-0.1)/0.004)^2)-0.8*exp(-((x[i]-0.9)/0.4)^2)))/x[1];
	f(2) = 1.d0
	do i = 2,n
		f(2) = f(2)*(2.d0-exp(-((x(i)-0.1d0)/0.004d0)**2.d0)-0.8d0*exp(-((x(i)-0.9d0)/0.4d0)**2.d0))
	enddo
	f(2) = f(2) / x(1)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

	lb = 0.d0;	    ub = 1.d0
	lb(1) = 0.1d0

	return
end subroutine setbounds
