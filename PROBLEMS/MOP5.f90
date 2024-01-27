!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example MOP5, Van Valedhuizen's test suit.
!#
!#   This file is part of a collection of problems developed for
!#   derivative-free multiobjective optimization in
!#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
!#   Direct Multisearch for Multiobjective Optimization, 2010.
!#
!#   Written by the authors in June 1, 2010.
!#
!###############################################################################
!# Define objective function variables
!var x{1..2} >= -30, <=30;
!
!# The objective functions
!minimize fobj1:
!    0.5*(x[1]^2+x[2]^2)+sin(x[1]^2+x[2]^2);
!
!minimize fobj2:
!    (3*x[1]-2*x[2]+4)^2/8+(x[1]-x[2]+1)^2/27+15;
!
!minimize fobj3:
!    1/(x[1]^2+x[2]^2+1)-1.1*exp(-x[1]^2-x[2]^2);

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
	integer	:: n, M, i
	real*8		:: x(n), f(M)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

!minimize fobj1:
!    0.5*(x[1]^2+x[2]^2)+sin(x[1]^2+x[2]^2);
	f(1) = 0.5d0*(x(1)**2.d0+x(2)**2.d0)+sin(x(1)**2.d0+x(2)**2.d0)
!minimize fobj2:
!    (3*x[1]-2*x[2]+4)^2/8+(x[1]-x[2]+1)^2/27+15;
	f(2) = (3.d0*x(1)-2.d0*x(2)+4.d0)**2.d0/8.d0+(x(1)-x(2)+1.d0)**2.d0/27.d0+15.d0
!minimize fobj3:
!    1/(x[1]^2+x[2]^2+1)-1.1*exp(-x[1]^2-x[2]^2);
	f(3) = 1.d0/(x(1)**2.d0+x(2)**2.d0+1.d0)-1.1d0*exp(-x(1)**2.d0-x(2)**2.d0)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -30.d0; 	    ub = 30.d0

	return
end subroutine setbounds
