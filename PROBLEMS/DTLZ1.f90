!###############################################################################
!#
!#   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
!#   Multi-Objective Optimization Test Problems", Congress on Evolutionary 
!#   Computation (CEC2002): 825-830, 2002.
!#
!#   Example DTLZ1.
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
!param M := 3, >=2;  # Number of objectives
!param n := 7, >= M; # Number of variables
!param k := n-M+1;
!param pi := 4*atan(1);
!
!# x variable
!var x{1..n};
!
!# g(x)
!var gx = 100*(k + (sum {i in M..n} ((x[i]-0.5)^2-cos(20*pi*(x[i]-0.5)))));
!
!# functions
!var ff1 = 0.5*(1+gx)*prod {j in 1..M-1} (x[j]);
!var ff {i in 2..M} = 0.5*(1+gx) * (prod {j in 1..M-i} (x[j]))*(1-x[M-i+1]);
!
!minimize f1:
!    ff1;
!minimize f {i in 2..M}:
!	ff[i];
!
!subject to bounds {i in 1..n}:
!    0.0 <= x[i] <= 1.0;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 7
	m = 0
	q = 3

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

	pi = 4.d0*atan(1.d0)
	k  = dble(n-qq+1)

	gx = 0.d0
	do i = qq,n
		gx = gx + ((x(i)-0.5d0)**2.d0-cos(20.d0*pi*(x(i)-0.5d0)))
	enddo
	gx = 100.d0*(k + gx)

	f(1) = 1.d0
	do i = 1,qq-1
		f(1) = f(1)*x(i)
	enddo
	f(1) = 0.5d0*(1.d0+gx)*f(1)

	do i = 2,qq
		f(i) = 1.d0
		do j = 1,qq-i
			f(i) = f(i)*x(j)
		enddo
		f(i) = 0.5d0*(1.d0+gx)*f(i)*(1.d0-x(qq-i+1))
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
