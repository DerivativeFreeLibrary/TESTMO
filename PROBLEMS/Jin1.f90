!###############################################################################
!#
!#   As described by Y. Jin, M. Olhofer and B. Sendhoff. "Dynamic weighted
!#   aggregation for evolutionary multi-objective optimization: Why does it
!#   work and how?", in Proceedings of Genetic and Evolutionary Computation 
!#   Conference, pp.1042-1049, San Francisco, USA, 2001.
!#
!#   Test function 1, F1.
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
!minimize f1:
!    (sum {i in 1..n} (x[i]^2))/n;
!minimize f2:
!    (sum {i in 1..n} ((x[i]-2)^2))/n;

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
	real*8		:: x(n), f(M)

	f(1) = 0.d0
	f(2) = 0.d0
	do i = 1,n
		f(1) = f(1) + x(1)**2.d0
		f(2) = f(2) + (x(1)-2.d0)**2.d0
	enddo
	f(1) = f(1)/dble(n)
	f(2) = f(2)/dble(n)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.d0; 	    ub = 1.d0

	return
end subroutine setbounds
