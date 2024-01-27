!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example QV1, see the previous cited paper for the original reference.
!#
!#   In the original reference the number of variables was n=16. 
!#   We selected n=10 as default.
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
!param pi := 4*atan(1);
!param n  := 10;
!
!# x variable
!var x{1..n} >=-5.12, <= 5.12;
!
!minimize f1:
!    (sum {i in 1..n} (x[i]^2-10*cos(2*pi*x[i])+10)/n)^0.25;
!minimize f2:
!    (sum {i in 1..n} ((x[i]-1.5)^2-10*cos(2*pi*(x[i]-1.5))+10)/n)^0.25;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 10
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

!minimize f1:
!    (sum {i in 1..n} (x[i]^2-10*cos(2*pi*x[i])+10)/n)^0.25;
	f(1) = 0.d0
	do i = 1,n
		f(1) = f(1) + (x(i)**2.d0-10.d0*cos(2.d0*pi*x(i))+10.d0)
	enddo
	f(1) = (f(1)/dble(n))**0.25d0

!minimize f2:
!    (sum {i in 1..n} ((x[i]-1.5)^2-10*cos(2*pi*(x[i]-1.5))+10)/n)^0.25;
	f(2) = 0.d0
	do i = 1,n
		f(2) = f(2) + ((x(i)-1.5d0)**2.d0-10.d0*cos(2.d0*pi*(x(i)-1.5d0))+10.d0)
	enddo
	f(2) = (f(2)/dble(n))**0.25d0

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

	lb = -5.12d0;	    ub = 5.12d0

	return
end subroutine setbounds
