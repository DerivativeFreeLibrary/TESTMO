!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example FES2, see the previous cited paper for the original reference.
!#
!#   In the above paper the number of variables was left undefined. 
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
!# parameter
!param pi := 4*atan(1);
!param n := 10;
!
!# x variable
!var x{1..n} >=0, <=1;
!
!minimize f1:
!    sum {i in 1..n} ((x[i]-0.5*cos(10*pi*i/n)-0.5)^2);
!minimize f2:
!    sum {i in 1..n} (abs(x[i]-sin(i-1)^2*cos(i-1)^2)^0.5);
!minimize f3:
!    sum {i in 1..n} (abs(x[i]-0.25*cos(i-1)*cos(2*i-2)-0.5)^0.5);

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 10
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

subroutine functs(n,x,q,f)
	implicit none
	integer	:: n, q, i
	real*8		:: x(n), f(q)
	real*8		:: pi

	pi = 4.d0*atan(1.d0)

	f = 0.d0
	do i = 1,n
		f(1) = f(1) + ((x(i)-0.5d0*cos(10.d0*pi*(dble(i)/dble(n)))-0.5d0)**2.d0)
		f(2) = f(2) + (abs(x(i)-sin(dble(i-1))**2.d0*cos(dble(i-1))**2.d0)**0.5d0)
		f(3) = f(3) + (abs(x(i)-0.25d0*cos(dble(i-1))*cos(dble(2*i-2))-0.5d0)**0.5d0)
	enddo

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.d0; 	    ub = 1.d0

	return
end subroutine setbounds
