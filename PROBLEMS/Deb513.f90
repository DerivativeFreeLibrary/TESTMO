!###############################################################################
!#
!#   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
!#   Difficulties and Construction of Test Problems", Evolutionary Computation 
!#   7(3): 205-230, 1999.
!#
!#   Example 5.1.3 (Discontinuous Pareto-optimal Front).
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
!param beta := 1;
!param alpha := 2;
!param q := 4;
!param pi := 4*atan(1);
!
!# x variable
!var x{1..2};
!
!# functions
!var ff1 = x[1];
!
!# g(x[2])
!var gx = 1+10*x[2];
!
!var h = 1-(ff1/gx)^alpha-(ff1/gx)*sin(2*pi*q*ff1);
!
!minimize f1:
!	ff1;
!minimize f2:
!    gx*h;
!
!subject to bounds1:
!	0.0 <= x[1] <= 1.0;
!subject to bounds2:
!	0.0 <= x[2] <= 1.0;

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

	x = 0.0d0
	
	return
end subroutine startp

subroutine functs(n,x,qq,f)
	implicit none
	integer	:: n,qq
	real*8		:: x(n), f(qq)
	real*8		:: gx, h
	real*8, parameter :: beta  = 1.0d0
	real*8, parameter :: alpha = 2.0d0
	real*8, parameter :: q     = 4.0d0
	real*8		:: pi

	pi = 4.d0*atan(1.d0);

	gx = 1.d0+10.d0*x(2);

	f(1) = x(1);

	h = 1.d0-(f(1)/gx)**alpha-(f(1)/gx)*sin(2.d0*pi*q*f(1));

	f(2) = gx*h

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.d0; 	    ub = 1.d0

	return
end subroutine setbounds
