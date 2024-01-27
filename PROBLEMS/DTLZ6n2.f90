!###############################################################################
!#
!#   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
!#   Multi-Objective Optimization Test Problems", Congress on Evolutionary 
!#   Computation (CEC2002): 825-830, 2002.
!#
!#   Example DTLZ6 with M=2 and n=2.
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
!param M := 2, >=2;  # Number of objectives
!param n := 2, >= M; # Number of variables
!param k := n-M+1;
!param pi := 4*atan(1);
!
!# x variable
!var x{1..n};
!
!# g(x)
!var gx = 1 + 9/k * sum {i in M..n} (x[i]);
!
!# functions
!var ffM = (1+gx) * (M - sum{i in 1..M-1} (x[i]/(1+gx)*(1+sin(3*pi*x[i]))));
!
!minimize f {i in 1..M-1}:
!    x[i];
!minimize fM:
!    ffM;
!
!subject to bounds {i in 1..n}:
!    0.0 <= x[i] <= 1.0;

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

subroutine functs(n,x,qq,f)
	implicit none
	integer	:: n, qq, i, j
	real*8		:: x(n), f(qq)
	real*8		:: gx, pi, k
	real*8, parameter :: alpha = 100.d0

	pi = 4.d0*atan(1.d0)
	k  = dble(n-qq+1)

	gx = 0.d0
	do i = qq,n
		gx = gx + x(i)
	enddo
	gx = 1.d0 + (9.d0/k)*gx

	do i = 1,qq-1
		f(i) = x(i)
	enddo

	f(qq) = 0.d0
	do i = 1,qq-1 
		f(qq) = f(qq) + (x(i)/(1.d0+gx)*(1.d0+sin(3.d0*pi*x(i))))
	enddo
	f(qq) = (1.d0+gx) * (dble(qq) - f(qq))

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.d0; 	    ub = 1.d0

	return
end subroutine setbounds
