!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example MHHM2, see the previous cited paper for the original reference.
!#
!#   This file is part of a collection of problems developed for
!#   derivative-free multiobjective optimization in
!#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
!#   Direct Multisearch for Multiobjective Optimization, 2010.
!#
!#   Written by the authors in June 1, 2010.
!#
!###############################################################################

!# x variable
!var x{1..2} >=0, <=1;

!minimize f1:
!    (x[1]-0.8)^2+(x[2]-0.6)^2;
!minimize f2:
!    (x[1]-0.85)^2+(x[2]-0.7)^2;
!minimize f3:
!    (x[1]-0.9)^2+(x[2]-0.6)^2;




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
	integer	:: n, M
	real*8		:: x(n), f(M)


	f(1) =(x(1)-0.8d0)**2.d0+(x(2)-0.6d0)**2.d0

	f(2) = (x(1)-0.85d0)**2.d0+(x(2)-0.7d0)**2.d0

	f(3) = (x(1)-0.9d0)**2.d0+(x(2)-0.6d0)**2.d0

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.d0; 	    ub = 1.d0
	

	return
end subroutine setbounds
