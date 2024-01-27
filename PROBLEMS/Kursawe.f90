!###############################################################################
!#
!#   As described by F. Kursawe in "A variant of evolution strategies for
!#   vector optimization", in H. P. Schwefel and R. Männer, editors, Parallel 
!#   Problem Solving from Nature, 1st Workshop, PPSN I, volume 496 of Lecture 
!#   Notes in Computer Science, pages 193-197, Berlin, Germany, Oct 1991, 
!#   Springer-Verlag.
!#
!#   In the above paper the variables bounds were not set.
!#   We considered -5.0 <= x[i] <= 5.0, i=1,2,3.
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
!param n := 3, >= 0; # Number of variables
!param pi := 4*atan(1);
!
!# x variable
!var x{1..n};
!
!minimize f1:
!    sum {i in 1..n-1} (-10*exp(-0.2*sqrt(x[i]^2+x[i+1]^2)));
!minimize f2:
!    sum {i in 1..n} (abs(x[i])^0.8+5*sin(x[i])^3);
!
!subject to bounds {i in 1..n}:
!    -5.0 <= x[i] <= 5.0;

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
	real*8		:: x(n), f(M), gx
	real*8, parameter :: pi = 4.d0*atan(1.d0);

!minimize f1:
!    sum {i in 1..n-1} (-10*exp(-0.2*sqrt(x[i]^2+x[i+1]^2)));
	f(1) = 0.d0
	do i = 1,n-1
		f(1) = f(1) + (-10.d0*exp(-0.2d0*sqrt(x(i)**2.d0+x(i+1)**2.d0)))
	enddo
!minimize f2:
!    sum {i in 1..n} (abs(x[i])^0.8+5*sin(x[i])^3);
	f(2) = 0.d0
	do i = 1,n
		f(2) = f(2) + (abs(x(i))**0.8d0+5.d0*sin(x(i))**3.d0)
	enddo

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -5.d0; 	    ub = 5.d0

	return
end subroutine setbounds
