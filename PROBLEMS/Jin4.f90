!###############################################################################
!#
!#   As described by Y. Jin, M. Olhofer and B. Sendhoff. "Dynamic weighted
!#   aggregation for evolutionary multi-objective optimization: Why does it
!#   work and how?", in Proceedings of Genetic and Evolutionary Computation 
!#   Conference, pp.1042-1049, San Francisco, USA, 2001.
!#
!#   Test function 4, F4.
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
!# parameters
!param n := 2;
!
!# x variable
!var x{1..n} >=0, <=1;
!
!# g(x) function
!var gx = 1+(9*sum{i in 2..n} x[i])/(n-1);
!
!minimize f1:
!    x[1];
!minimize f2:
!    gx*(1-(x[1]/gx)^0.25-(x[1]/gx)^4);

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

	!call setbounds(n,l,u)

	x = 0.d0
	
	return
end subroutine startp

subroutine functs(n,x,M,f)
	implicit none
	integer	:: n, M, i
	real*8		:: x(n), f(M), gx

!# g(x) function
!var gx = 1+(9*sum{i in 2..n} x[i])/(n-1);
	gx = 0.d0
	do i = 2,n
		gx = gx + x(i)
	enddo
	gx = 1.d0 + (9.d0*gx)/dble(n-1)
!
!minimize f1:
!    x[1];
	f(1) = x(1)
!minimize f2:
!    gx*(1-(x[1]/gx)^0.25-(x[1]/gx)^4);
	f(2) = gx*(1.d0-(x(1)/gx)**0.25d0-(x(1)/gx)**4.d0)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.d0; 	    ub = 1.d0

	return
end subroutine setbounds
