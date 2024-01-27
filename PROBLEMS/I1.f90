!###############################################################################
!#
!#   As described by Huband et al. in "A Scalable Multi-objective Test Problem
!#   Toolkit", in C. A. Coello Coello et al. (Eds.): EMO 2005, LNCS 3410, 
!#   pp. 280295, 2005, Springer-Verlag Berlin Heidelberg 2005.
!#
!#   Example I1.
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
!param S {m in 1..M} := 1;
!
!# neq WFG3
!param A {i in 1..M-1} := 1;
!
!# problem variables
!param zmax {i in 1..n} := 1;
!var z{i in 1..n} >=0, <= zmax[i];
!
!# transform z into [0,1] set
!var y{i in 1..n} = z[i]/zmax[i];
!
!# first level mapping
!var t1{i in 1..n} = y[i];
!
!# second level mapping
!var t2{i in 1..n} = if i<=k then t1[i]
!    else abs(t1[i]-0.35)/abs(floor(0.35-t1[i])+0.35);
!
!# third level mapping
!param w{i in 1..n} := 1;
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

	n = 4+4
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
	real*8		:: pi, pi2
	real*8		:: t1(8), t2(8), zmax(8), y(8), w(8)
!	real*8		:: t1(n), t2(n), zmax(n), y(n), w(n)
	real*8		:: t3(3), x(3), h(3), S(3)
!	real*8		:: t3(M), x(M), h(M), S(M)
	real*8		:: Av(2)
!	real*8		:: Av(M-1)
	real*8		:: a, b
	integer, parameter :: k = 4
	integer, parameter :: l = 4
	real*8, parameter  :: kr = 4.d0
	real*8, parameter  :: lr = 4.d0

	pi  = 4.d0*atan(1.d0)
	pi2 = 2.d0*atan(1.d0)

	do i = 1,M
		S(i) = 1.d0
	enddo
	do i = 1,M-1
		Av(i) = 1.d0
	enddo
	do i = 1,n
		zmax(i) = 1.d0
		y(i) = z(i)/zmax(i)
		t1(i) = y(i)
	enddo

!# second level mapping
!var t2{i in 1..n} = if i<=k then t1[i]
!    else abs(t1[i]-0.35)/abs(floor(0.35-t1[i])+0.35);
	do i = 1,k
		t2(i) = t1(i)
	enddo
	do i = k+1,n
		t2(i) = abs(t1(i)-0.35d0)/abs(floor(0.35d0-t1(i))+0.35d0)
	enddo

!# third level mapping
!param w{i in 1..n} := 1;
!var t3{i in 1..M} = if i<=M-1 then (sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (w[j]*t2[j]))/(sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} w[j])
!    else (sum {j in k+1..n} (w[j]*t2[j]))/(sum {j in k+1..n} w[j]);
	do i = 1,n
		w(i) = 1.d0
	enddo
	do i = 1,M-1
		a = 0.d0
		b = 0.d0
		do j = ((i-1)*k/(M-1)+1), (i*k/(M-1)) 
			a = a + (w(j)*t2(j))
			b = b + w(j)
		enddo
		t3(i) = a/b
	enddo
	a = 0.d0
	b = 0.d0
	do j = k+1,n 
		a = a + (w(j)*t2(j))
		b = b + w(j)
	enddo
	t3(M) = a/b

!# Define objective function variables
!var x{i in 1..M} =	 if i<=M-1 then max(t3[M],A[i])*(t3[i]-0.5)+0.5
!    else t3[M];
	do i = 1,M-1
		x(i) = max(t3(M),Av(i))*(t3(i)-0.5d0)+0.5d0
	enddo
	x(M) = t3(M)

!# Define objective function function h
!var h{m in 1..M} = if m==1 then prod {i in 1..M-1} sin(x[i]*pi2)
!    else if m<=M-1 then (prod {i in 1..M-m} sin(x[i]*pi2))*cos(x[M-m+1]*pi2)
!        else cos(x[1]*pi2);
	h(1) = 1.d0
	do j = 1,M-1
		h(1) = h(1)*sin(x(j)*pi2)
	enddo
	do i = 2,M-1
		h(i) = 1.d0
		do j = 1,M-i
			h(i) = h(i)*sin(x(j)*pi2)
		enddo
		h(i) = h(i)*cos(x(M-i+1)*pi2)
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
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.d0; 	    ub = 1.d0

	return
end subroutine setbounds
