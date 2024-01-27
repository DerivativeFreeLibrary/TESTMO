!###############################################################################
!#
!#
!#   As described by A. Lovison in "A synthetic approach to multiobjective
!#   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
!#
!#   Example 1.
!#
!#   In the above paper/papers the variables bounds were not set.
!#   We considered 0<=x[i]<=3, i=1,2.
!#
!#   This file is part of a collection of problems developed for
!#   derivative-free multiobjective optimization in
!#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
!#   Direct Multisearch for Multiobjective Optimization, 2010.
!#
!#
!#   Written by the authors in June 1, 2010.
!#
!###############################################################################

!# x variable
!var x{1..2} >=0, <=3;

!maximize f1:
!	-1.05*x[1]^2-0.98*x[2]^2;
!maximize f2:
!    -0.99*(x[1]-3)^2-1.03*(x[2]-2.5)^2;


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
	integer	:: n, M
	real*8		:: x(n), f(M)

	f(1) = 1.05d0*x(1)**2.d0+0.98d0*x(2)**2.d0
	f(2) = 0.99d0*(x(1)-3.d0)**2.d0+1.03d0*(x(2)-2.5d0)**2.d0


	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.d0; 	    ub = 3.d0

	return
end subroutine setbounds
