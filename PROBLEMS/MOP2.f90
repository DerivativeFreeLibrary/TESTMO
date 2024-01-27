!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example MOP2, Van Valedhuizen's test suit.
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
!# Number of variables
!param n := 4; #??
!# Define objective function variables
!var x{1..n} >= -4, <=4;
!
!# The objective functions
!minimize fobj1:
!    1-exp(-sum {i in 1..n} (x[i]-1/sqrt(n))^2);
!
!minimize fobj2:
!    1-exp(-sum {i in 1..n} (x[i]+1/sqrt(n))^2);

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 4
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
!    1-exp(-sum {i in 1..n} (x[i]-1/sqrt(n))^2);
	f(1) = 0.d0
	do i = 1,n
		f(1) = f(1) + (x(i)-1.d0/sqrt(dble(n)))**2.d0
	enddo
	f(1) = 1.d0 - exp(-f(1))

!minimize fobj2:
!    1-exp(-sum {i in 1..n} (x[i]+1/sqrt(n))^2);
	f(2) = 0.d0
	do i = 1,n
		f(2) = f(2) + (x(i)+1.d0/sqrt(dble(n)))**2.d0
	enddo
	f(2) = 1.d0 - exp(-f(2))

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -4.d0; 	    ub = 4.d0

	return
end subroutine setbounds
