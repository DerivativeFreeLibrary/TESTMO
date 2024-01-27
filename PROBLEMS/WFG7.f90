!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#           
!#   Example WFG7.
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
!param M := 3;
!param k := 4;
!param l := 4;
!param n := k+l;
!
!param pi := 4*atan(1);
!param pi2:= 2*atan(1);
!
!param S {m in 1..M} := 2*m;
!
!# neq WFG3
!param A {i in 1..M-1} := 1;
!
!# problem variables
!param zmax {i in 1..n} := 2*i;
!var z{i in 1..n} >=0, <= zmax[i];
!
!# transform z into [0,1] set
!var y{i in 1..n} = z[i]/zmax[i];
!
!# first level mapping
!param w{1..n}, default 1.0;
!param AA := 0.98/49.98;
!param BB := 0.02;
!param CC := 50;
!var r_sum{i in 1..k} = (sum {j in i+1..n} (w[j]*y[j]))/(sum {j in i+1..n} w[j]);
!var t1{i in 1..n} = if i<=k then y[i]^(BB+(CC-BB)*(AA-(1-2*r_sum[i])*abs(floor(0.5-r_sum[i])+AA)))
!    else y[i];
!
!# second level mapping
!var t2{i in 1..n} = if i <= k then t1[i] 
!      else abs(t1[i]-0.35)/abs(floor(0.35-t1[i])+0.35);
!      
!# third level mapping
!# w already defined
!var t3{i in 1..M} = if i<=M-1 then (sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (w[j]*t2[j]))/(sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} w[j])
!    else (sum {j in k+1..n} (w[j]*t2[j]))/(sum {j in k+1..n} w[j]);
!
!# Define objective function variables
!var x{i in 1..M} = if i<=M-1 then max(t3[M],A[i])*(t3[i]-0.5)+0.5
!    else t3[M];
!
!# Define objective function function h
!var h{m in 1..M} = if m==1 then prod {i in 1..M-1} sin(x[i]*pi2)
!    else if m<=M-1 then (prod {i in 1..M-m} sin(x[i]*pi2))*cos(x[M-m+1]*pi2)
!        else cos(x[1]*pi2);
!
!# The objective functions
!minimize fobj {m in 1..M}:
!    x[M]+S[m]*h[m];

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 8
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

subroutine functs(n,z,M,f)
	implicit none
	integer	:: n, M, i, j
	real*8		:: z(n), f(M)
	integer, parameter :: k = 4
	integer, parameter :: l = 4
	real*8, parameter :: pi = 4.d0*atan(1.d0)
	real*8, parameter :: pi2 = 2.d0*atan(1.d0)
	real*8		:: S(3), zmax(8), A(2), y(8), h(3), rsum(4)
	real*8		:: t1(8), t2(8), t3(3), w(8), t4(3), x(3)
!	real*8		:: S(M), zmax(n), A(M-1), y(n), h(M), rsum(k)
!	real*8		:: t1(n), t2(n), t3(M), w(n), t4(M), x(M)
	real*8		:: gsum1, gsum2
	real*8, parameter :: AA = 0.98d0/49.98d0;
	real*8, parameter :: BB = 0.02d0;
	real*8, parameter :: CC = 50.d0;
	real*8, parameter :: AAA = 0.02d0;
	real*8, parameter :: alpha = 1.d0
	real*8, parameter :: AAAA = 5.d0

!param S {m in 1..M} := 2*m;
	do i = 1,M
		S(i) = 2.d0*dble(i)
	enddo
!param zmax {i in 1..n} := 2*i;
	do i = 1,n
		zmax(i) = 2.d0*dble(i)
	enddo
!param A {i in 1..M-1} := 1;
	do i = 1,M-1
		A(i) = 1.d0
	enddo

!# problem variables
!# transform z into [0,1] set
!var y{i in 1..n} = z[i]/zmax[i];
	do i = 1,n
		y(i) = z(i)/zmax(i)
	enddo

!# first level mapping
!var r_sum{i in 1..k} = (sum {j in i+1..n} (w[j]*y[j]))/(sum {j in i+1..n} w[j]);
!var t1{i in 1..n} = if i<=k then y[i]^(BB+(CC-BB)*(AA-(1-2*r_sum[i])*abs(floor(0.5-r_sum[i])+AA)))
!    else y[i];
	do i = 1,n
		w(i) = 1.d0
	enddo
	do i = 1,k
		gsum1 = 0.d0
		gsum2 = 0.d0
		do j = i+1,n
			gsum1 = gsum1 + w(j)*y(j)
			gsum2 = gsum2 + w(j)
		enddo
		rsum(i) = gsum1 / gsum2
	enddo
	do i = 1,k
		t1(i) = y(i)**(BB+(CC-BB)*(AA-(1.d0-2.d0*rsum(i))*abs(floor(0.5d0-rsum(i))+AA)))
	enddo
	do i = k+1,n
		t1(i) = y(i)
	enddo

!# second level mapping
!var t2{i in 1..n} = if i <= k then t1[i] 
!      else abs(t1[i]-0.35)/abs(floor(0.35-t1[i])+0.35);
	do i = 1,k
		t2(i) = t1(i)
	enddo
	do i = k+1,n
		t2(i) = abs(t1(i)-0.35d0)/abs(floor(0.35d0-t1(i))+0.35d0)
	enddo

!# third level mapping
!# w already defined
!var t3{i in 1..M} = if i<=M-1 then (sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (w[j]*t2[j]))/(sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} w[j])
!    else (sum {j in k+1..n} (w[j]*t2[j]))/(sum {j in k+1..n} w[j]);
	do i = 1,M-1
		gsum1 = 0.d0
		gsum2 = 0.d0
		do j = ((i-1)*k/(M-1)+1),(i*k/(M-1))
			gsum1 = gsum1 + w(j)*t2(j)
			gsum2 = gsum2 + w(j)
		enddo
		t3(i) = gsum1 / gsum2
	enddo
	gsum1 = 0.d0
	gsum2 = 0.d0
	do j = k+1,n
		gsum1 = gsum1 + w(j)*t2(j)
		gsum2 = gsum2 + w(j)
	enddo
	t3(M) = gsum1 / gsum2

!# Define objective function variables
!var x{i in 1..M} = if i<=M-1 then max(t3[M],A[i])*(t3[i]-0.5)+0.5
!    else t3[M];
	do i = 1,M-1
		x(i) = max(t3(M),A(i))*(t3(i)-0.5d0)+0.5d0
	enddo
	x(M) = t3(M)

!# Define objective function function h
!var h{m in 1..M} = if m==1 then prod {i in 1..M-1} sin(x[i]*pi2)
!    else if m<=M-1 then (prod {i in 1..M-m} sin(x[i]*pi2))*cos(x[M-m+1]*pi2)
!        else cos(x[1]*pi2);
	h(1) = 1.d0
	do i = 1,M-1
		h(1) = h(1)*sin(x(i)*pi2)
	enddo
	do j = 2,M-1
		h(j) = 1.d0
		do i = 1,M-j
			h(j) = h(j)*sin(x(i)*pi2)
		enddo
		h(j) = h(j)*cos(x(M-j+1)*pi2)
	enddo
	h(M) = cos(x(1)*pi2)

!# The objective functions
!minimize fobj {m in 1..M}:
!    x[M]+S[m]*h[m];
	do i = 1,M
		f(i) = x(M)+S(i)*h(i)
	enddo

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n, i
	real*8		:: lb(n), ub(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

	lb = 0.d0
	do i = 1,n
		ub(i) = 2.d0*dble(i)
	enddo

	return
end subroutine setbounds
