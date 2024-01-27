!###############################################################################
!#
!#   As described by T. Okabe, Y. Jin, M. Olhofer, and B. Sendhoff. "On test
!#   functions for evolutionary multi-objective optimization.", Parallel
!#   Problem Solving from Nature, VIII, LNCS 3242, Springer, pp.792-802,
!#   September 2004.
!#
!#   Test function OKA2.
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
!param pi := 4*tan(1);
!
!# x variable
!var x{1..3};
!
!minimize f1:
!    x[1];
!minimize f2:
!    1-(x[1]+pi)^2/(4*pi^2)+abs(x[2]-5*cos(x[1]))^(1/3)+abs(x[3]-5*sin(x[1]))^(1/3);
!
!# Simple bounds
!subject to b1:
!    -pi<= x[1] <= pi;
!
!subject to b2 {i in 2..3}:
!    -5<= x[i] <=5;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 3
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
	real*8		:: x(n), f(M), y(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

!minimize f1:
!    x[1];
	f(1) = x(1)
!minimize f2:
!    1-(x[1]+pi)^2/(4*pi^2)+abs(x[2]-5*cos(x[1]))^(1/3)+abs(x[3]-5*sin(x[1]))^(1/3);
	f(2) = 1.d0-(x(1)+pi)**2.d0/(4.d0*pi**2.d0)+abs(x(2)-5.d0*cos(x(1)))**(1.d0/3.d0)+ &
		abs(x(3)-5.d0*sin(x(1)))**(1.d0/3.d0)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

	lb(1) = -pi; 	    ub(1) = pi
	lb(2) = -5.d0;	    ub(2) = 5.d0
	lb(3) = -5.d0;	    ub(3) = 5.d0

	return
end subroutine setbounds
