!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example MOP6, Van Valedhuizen's test suit.
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
!# Pi
!param pi := 4*atan(1);
!
!# Define objective function variables
!var x{1..2} >= 0, <=1;
!
!# The objective functions
!minimize fobj1:
!    x[1];
!
!minimize fobj2:
!    (1+10*x[2])*(1-(x[1]/(1+10*x[2]))^2-x[1]/(1+10*x[2])*sin(8*pi*x[1]));

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
	integer	:: n, M, i
	real*8		:: x(n), f(M)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

!minimize fobj1:
!    x[1];
	f(1) = x(1)
!minimize fobj2:
!    (1+10*x[2])*(1-(x[1]/(1+10*x[2]))^2-x[1]/(1+10*x[2])*sin(8*pi*x[1]));
	f(2) = (1.d0+10.d0*x(2))*(1.d0-(x(1)/(1.d0+10.d0*x(2)))**2.d0	&
		-x(1)/(1.d0+10.d0*x(2))*sin(8*pi*x(1)))

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.d0; 	    ub = 1.d0

	return
end subroutine setbounds
