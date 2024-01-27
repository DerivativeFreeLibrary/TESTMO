!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example MLF2, See the previous cited paper for the original reference.
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
!var x{1..2} >=-2, <=2;
!var y{i in 1..2} = 2*x[i];
!
!maximize f1:
!    5-((x[1]^2+x[2]-11)^2+(x[1]+x[2]^2-7)^2)/200;
!maximize f2:
!    5-((y[1]^2+y[2]-11)^2+(y[1]+y[2]^2-7)^2)/200;

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
	real*8		:: x(n), f(M), y(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

!var y{i in 1..2} = 2*x[i];
	y = 2.d0*x

!maximize f1:
!    5-((x[1]^2+x[2]-11)^2+(x[1]+x[2]^2-7)^2)/200;
	f(1) = 5.d0-((x(1)**2.d0+x(2)-11.d0)**2.d0+(x(1)+x(2)**2.d0-7.d0)**2.d0)/200.d0
	f(1) = -f(1)

!maximize f2:
!    5-((y[1]^2+y[2]-11)^2+(y[1]+y[2]^2-7)^2)/200;
	f(2) = 5.d0-((y(1)**2.d0+y(2)-11.d0)**2.d0+(y(1)+y(2)**2.d0-7.d0)**2.d0)/200.d0;
	f(2) = -f(2)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -2.d0; 	    ub = 2.d0

	return
end subroutine setbounds
