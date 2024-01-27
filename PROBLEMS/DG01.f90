!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example DG01, see the previous cited paper for the original reference.
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
!var x >=-10, <=13;
!
!minimize f1:
!	sin(x);
!minimize f2:
!    sin(x+0.7);

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 1
	m = 0
	q = 2

	return
end subroutine setdim

subroutine startp(n,x)
	implicit none
	integer	:: n
	real*8		:: x(n), l(n), u(n)

	!call setbounds(n,l,u)

	x = 0.0d0
	
	return
end subroutine startp

subroutine functs(n,x,qq,f)
	implicit none
	integer	:: n,qq
	real*8		:: x(n), f(qq)

	f(1) = sin(x(1))
	f(2) = sin(x(1)+0.7d0)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -10.d0; 	    ub = 13.d0

	return
end subroutine setbounds
