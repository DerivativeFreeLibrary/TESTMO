!###############################################################################
!#
!#   As described by A. Lovison in "A synthetic approach to multiobjective
!#   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
!#
!#   Example 5.
!#
!#   In the above paper/papers the variables bounds were not set.
!#   We considered -1<=x[i]<=4, i=1,2,3.
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
!param n := 3;           # number of variables
!param m := 4;           # number of function used in the objectives
!param C{1..n,1..m};     #:=Uniform(-1,1); # distinct non collinear points
!param alpha{1..m,1..n}; #:=Uniform(0,1) # negative definite matrix diagonal
!param beta{1..m};       #:=Uniform(-1,1);
!param gamma{1..m};      #:=Uniform(-1,1);
!param pi := 4*atan(1);
!
!# x variable
!var x{1..n} >= -1, <= 4;
!var f{j in 1..m} = sum{i in 1..n} (-alpha[j,i]*(x[i]-C[i,j])^2);
!
!maximize u1:
!	f[1]+beta[1]*exp(f[4]/gamma[1]);
!maximize u2:
!    f[2]+beta[2]*sin(pi*(x[1]+x[2])/gamma[2]);
!maximize u3:
!    f[3]+beta[3]*cos(pi*(x[1]-x[2])/gamma[3]);

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 3
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

subroutine functs(n,x,M,u)
	implicit none
	integer	:: n, M, i, j
	real*8		:: x(n), f(4), u(M)
	real*8, parameter :: pi = 4.d0*atan(1.d0);
	real*8 	:: C(4,3), alpha(3,4)
	real*8		:: beta(4), gamma(4)
	save C, alpha, beta, gamma
	data C / &
	0.218418, -0.620254, 0.843784, 0.914311, &
	-0.788548, 0.428212, 0.103064, -0.47373, &
	-0.300792, -0.185507, 0.330423, 0.151614 /
	data alpha / &
	0.942022, 0.363525, 0.00308876, &
	0.755598, 0.450103, 0.170122, &
	0.787748, 0.837808, 0.590166, &
	0.203093, 0.253639, 0.532339 /
	data beta / -0.666503, -0.945716, -0.334582, 0.611894 /
	data gamma / 0.281032, 0.508749, -0.0265389, -0.920133 /

!var f{j in 1..m} = sum{i in 1..n} (-alpha[j,i]*(x[i]-C[i,j])^2);
	do j = 1,4
		f(j) = 0.d0
		do i = 1,n
			f(j) = f(j) - alpha(i,j)*(x(i)-C(j,i))**2.d0
		enddo
	enddo

!
!maximize u1:
!	f[1]+beta[1]*exp(f[4]/gamma[1]);
	u(1) = -f(1)-beta(1)*exp(f(4)/gamma(1))

!maximize u2:
!    f[2]+beta[2]*sin(pi*(x[1]+x[2])/gamma[2]);
	u(2) = -f(2)-beta(2)*sin(pi*(x(1)+x(2))/gamma(2))

!maximize u3:
!    f[3]+beta[3]*cos(pi*(x[1]-x[2])/gamma[3]);
	u(3) = -f(3)-beta(3)*cos(pi*(x(1)-x(2))/gamma(3))

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -1.d0; 	    ub = 4.d0

	return
end subroutine setbounds
