!###############################################################################
!#
!#   As described by A. Lovison in "A synthetic approach to multiobjective
!#   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
!#
!#   Example 4.
!#
!#   In the above paper/papers the variables bounds were not set.
!#   We considered 0<=x[1]<=6 and -1<=x[2]<=1.
!#
!#   This file is part of a collection of problems developed for
!#   derivative-free multiobjective optimization in
!#   A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
!#   Direct Multisearch for Multiobjective Optimization, 2010.
!#
!#   Written by the authors in June 1, 2010.
!#
!###############################################################################

!# x variable
!var x{1..2};

!maximize f1:
!	-x[1]^2-x[2]^2-4*(exp(-(x[1]+2)^2-x[2]^2)+exp(-(x[1]-2)^2-x[2]^2));
!maximize f2:
!    -(x[1]-6)^2-(x[2]+0.5)^2;

!subject to xbound1:
!    0 <= x[1] <= 6;

!subject to xbound2:
!    -1 <= x[2] <= 1;





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

	f(1) =x(1)**2.d0+x(2)**2.d0+4.d0*(exp(-(x(1)+2.d0)**2.d0-x(2)**2.d0)-exp(-(x(1)-2.d0)**2.d0-x(2)**2.d0));

	f(2) = (x(1)-6.d0)**2.d0+(x(2)+0.5d0)**2.d0


	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb(1) = 0.d0; 	    ub(1) = 6.d0
	lb(2) = -1.d0; 	    ub(2) = 1.d0

	return
end subroutine setbounds
