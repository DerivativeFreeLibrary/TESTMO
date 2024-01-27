!###############################################################################
!#
!#   As described by C.M. Fonseca and P.J. Fleming in "Multiobjective
!#   Optimization and Multiple Constraint Handling with Evolutionary
!#   Algorithms Part I: A Unified Formulation", in IEEE Transactions 
!#   on Systems, Man, and Cybernetics Part A: Systems and Humans, 
!#   vol. 28, no. 1, January 1998.
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
!    1-exp(-(x[1]-1)^2-(x[2]+1)^2);
!minimize f2:
!    1-exp(-(x[1]+1)^2-(x[2]-1)^2);
!
!
!subject to bounds {i in 1..2}:
!    -4.0 <= x[i] <= 4.0;

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

subroutine functs(n,x,q,f)
	implicit none
	integer	:: n, q, i
	real*8		:: x(n), f(q)
	real*8		:: pi

	pi = 4.d0*atan(1.d0)

	f(1) = 1.d0-exp(-(x(1)-1.d0)**2.d0-(x(2)+1.d0)**2.d0)
	f(2) = 1.d0-exp(-(x(1)+1.d0)**2.d0-(x(2)-1.d0)**2.d0)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -4.d0; 	    ub = 4.d0

	return
end subroutine setbounds
