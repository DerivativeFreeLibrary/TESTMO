!###############################################################################
!#
!#   As described by T. Okabe, Y. Jin, M. Olhofer, and B. Sendhoff. "On test
!#   functions for evolutionary multi-objective optimization.", Parallel
!#   Problem Solving from Nature, VIII, LNCS 3242, Springer, pp.792-802,
!#   September 2004.
!#
!#   Test function OKA1.
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
!var x{1..2};
!
!# y variable
!var y{i in 1..2} = if i<=1 then cos(pi/12)*x[1]-sin(pi/12)*x[2]
!    else sin(pi/12)*x[1]+cos(pi/12)*x[2];
!
!minimize f1:
!    y[1];
!minimize f2:
!    sqrt(2*pi)-sqrt(abs(y[1]))+2*abs(y[2]-3*cos(y[1])-3)^(1/3);
!
!# Simple bounds
!subject to b1:
!    6*sin(pi/12)<= x[1] <= 6*sin(pi/12)+2*pi*cos(pi/12);
!    
!subject to b2:
!    -2*pi*sin(pi/12)<= x[2] <=6*cos(pi/12);

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

	call setbounds(n,l,u)

	x = (l+u)/2.d0
	
	return
end subroutine startp

subroutine functs(n,x,M,f)
	implicit none
	integer	:: n, M, i
	real*8		:: x(n), f(M), y(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

!# y variable
!var y{i in 1..2} = if i<=1 then cos(pi/12)*x[1]-sin(pi/12)*x[2]
!    else sin(pi/12)*x[1]+cos(pi/12)*x[2];
	y(1) = cos(pi/12.d0)*x(1)-sin(pi/12.d0)*x(2)
	y(2) = sin(pi/12.d0)*x(1)+cos(pi/12.d0)*x(2)

!minimize f1:
!    y[1];
	f(1) = y(1)

!minimize f2:
!    sqrt(2*pi)-sqrt(abs(y[1]))+2*abs(y[2]-3*cos(y[1])-3)^(1/3);
	f(2) = sqrt(2.d0*pi)-sqrt(abs(y(1)))+2.d0*abs(y(2)-3.d0*cos(y(1))-3.d0)**(1.d0/3.d0)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

	lb(1) = 6.d0*sin(pi/12.d0); 	    ub(1) = 6.d0*sin(pi/12.d0)+2.d0*pi*cos(pi/12.d0)
	lb(2) = -2.d0*pi*sin(pi/12.d0);	    ub(2) = 6.d0*cos(pi/12.d0)

	return
end subroutine setbounds
