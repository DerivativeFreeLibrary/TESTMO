!###############################################################################
!#
!#   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
!#   Difficulties and Construction of Test Problems", Evolutionary Computation 
!#   7(3): 205-230, 1999.
!#              
!#   Example 5.1.2 (Non-convex local and convex global Pareto-optimal Front).
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
!
!
!# x variable
!var x{1..2};
!
!# functions
!var ff1 = 4*x[1];
!
!# g(x[2])
!var gx = if (x[2]<=0.4) then 4-3*exp(-((x[2]-0.2)/0.02)^2)
!            else 4-2*exp(-((x[2]-0.7)/0.2)^2);
!
!var alpha = 0.25+3.75*(gx-1);
!
!var h = if (ff1<=beta*gx) then (1-(ff1/(beta*gx))^alpha) else 0;
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

subroutine functs(n,x,q,f)
	implicit none
	integer	:: n,q
	real*8		:: x(n), f(q)
	real*8		:: gx, h, alpha
	real*8, parameter :: beta  = 1.0d0

	if (x(2) <= 0.4d0) then 
		gx = 4.d0-3.d0*exp(-((x(2)-0.2d0)/0.02d0)**2.d0)
	else 
		gx = 4.d0-2.d0*exp(-((x(2)-0.7d0)/0.2d0)**2.d0);
	endif 

	f(1) = 4.d0*x(1);
	alpha = 0.25d0+3.75d0*(gx-1.d0);

	if (f(1) <= beta*gx) then 
		h = (1.d0-(f(1)/(beta*gx))**alpha) 
	else 
		h = 0.d0
	endif

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
