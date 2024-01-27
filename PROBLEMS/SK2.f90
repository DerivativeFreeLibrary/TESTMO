!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example SK2, see the previous cited paper for the original reference.
!#
!#   In the above paper/papers the variables bounds were not set.
!#   We considered -10<=x[i]<=10, i=1,2,3,4.
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
!var x{1..4} >=-10, <= 10;
!
!maximize f1:
!    -(x[1]-2)^2-(x[2]+3)^2-(x[3]-5)^2-(x[4]-4)^2+5;
!maximize f2:
!    (sum {i in 1..4} sin(x[i]))/(1+(sum {i in 1..4} x[i]^2)/100);

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
	real*8		:: a, b

!maximize f1:
!    -(x[1]-2)^2-(x[2]+3)^2-(x[3]-5)^2-(x[4]-4)^2+5;
	f(1) = -(x(1)-2.d0)**2.d0-(x(2)+3.d0)**2.d0-(x(3)-5.d0)**2.d0-(x(4)-4.d0)**2.d0+5.d0 

!maximize f2:
!    (sum {i in 1..4} sin(x[i]))/(1+(sum {i in 1..4} x[i]^2)/100);
	a = 0.d0
	b = 0.d0
	do i = 1,n
		a = a + sin(x(i))
		b = b + x(i)**2.d0
	enddo
	f(2) = a / (1.d0 + b/100.d0)

	f = -f

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

	lb = -10.d0;	    ub = 10.d0

	return
end subroutine setbounds
