!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#           
!#   Example WFG6.
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
!var t1{i in 1..n} = if i <= k then y[i] 
!      else abs(y[i]-0.35)/abs(floor(0.35-y[i])+0.35);
!
!# second level mapping
!var t2{i in 1..M} = if i<=M-1 then sum {ii in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (t1[ii]+sum {jj in 0..(k/(M-1)-2)} abs(t1[ii]-t1[((i-1)*k/(M-1)+1)+((ii+jj-((i-1)*k/(M-1)+1)+1) mod ((i*k/(M-1))-((i-1)*k/(M-1)+1)+1))]))/(((i*k/(M-1))-((i-1)*k/(M-1)+1)+1)/(k/(M-1))*ceil(k/(M-1)/2)*(1+2*k/(M-1)-2*ceil(k/(M-1)/2)))
!    else sum {ii in k+1..n} (t1[ii]+sum {jj in 0..(l-2)} abs(t1[ii]-t1[k+1+((ii+jj-(k+1)+1) mod (n-k))]))/(((n-k)/l)*ceil(l/2)*(1+2*l-2*ceil(l/2)));
!
!# Define objective function variables
!var x{i in 1..M} = if i<=M-1 then max(t2[M],A[i])*(t2[i]-0.5)+0.5
!    else t2[M];
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
	integer	:: n, M, i, j, ii, jj
	real*8		:: z(n), f(M)
	integer, parameter :: k = 4
	integer, parameter :: l = 4
	real*8, parameter :: pi = 4.d0*atan(1.d0)
	real*8, parameter :: pi2 = 2.d0*atan(1.d0)
	real*8		:: S(3), zmax(8), A(2), y(8), h(3)
	real*8		:: t1(8), t2(3), t3(8), w(8), t4(3), x(3)
!	real*8		:: S(M), zmax(n), A(M-1), y(n), h(M)
!	real*8		:: t1(n), t2(M), t3(n), w(n), t4(M), x(M)
	real*8		:: gsum1, gsum2
	real*8, parameter :: AA = 0.35d0;
	real*8, parameter :: BB = 0.001d0;
	real*8, parameter :: CC = 0.05d0;
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
!var t1{i in 1..n} = if i <= k then y[i] 
!      else abs(y[i]-0.35)/abs(floor(0.35-y[i])+0.35);
	do i = 1,k
		t1(i) = y(i)
	enddo
	do i = k+1,n
		t1(i) = abs(y(i)-0.35d0)/abs(floor(0.35d0-y(i))+0.35d0)
	enddo

!# second level mapping
!var t2{i in 1..M} = if i<=M-1 then sum {ii in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (t1[ii]+sum {jj in 0..(k/(M-1)-2)} abs(t1[ii]-t1[((i-1)*k/(M-1)+1)+((ii+jj-((i-1)*k/(M-1)+1)+1) mod ((i*k/(M-1))-((i-1)*k/(M-1)+1)+1))]))/(((i*k/(M-1))-((i-1)*k/(M-1)+1)+1)/(k/(M-1))*ceil(k/(M-1)/2)*(1+2*k/(M-1)-2*ceil(k/(M-1)/2)))
!    else sum {ii in k+1..n} (t1[ii]+sum {jj in 0..(l-2)} abs(t1[ii]-t1[k+1+((ii+jj-(k+1)+1) mod (n-k))]))/(((n-k)/l)*ceil(l/2)*(1+2*l-2*ceil(l/2)));
	do i = 1,M-1
		gsum2 = 0.d0
		do ii = ((i-1)*k/(M-1)+1), (i*k/(M-1))
			gsum1 = 0.d0
			do jj = 0, (k/(M-1)-2)
				gsum1 = gsum1 + abs(t1(ii)-t1(((i-1)*k/(M-1)+1)+mod((ii+jj-((i-1)*k/(M-1)+1)+1),((i*k/(M-1))-((i-1)*k/(M-1)+1)+1))))
			enddo
			gsum2 = gsum2 + t1(ii) + gsum1
		enddo
		t2(i) = gsum2 / (((dble(i)*dble(k)/dble(M-1))-(dble(i-1)*dble(k)/dble(M-1)+1.d0)+1.d0)/  &
				 (dble(k)/dble(M-1))*ceiling(dble(k)/dble(M-1)/2.d0)*(1.d0+2.d0*dble(k)/ &
					dble(M-1)-2*ceiling(dble(k)/dble(M-1)/2.d0)))
	enddo
	t2(M) = 0.d0
	do ii = k+1,n
		gsum1 = 0.d0
		do jj = 0, l-2
			gsum1 = gsum1 + abs(t1(ii)-t1(k+1+mod((ii+jj-(k+1)+1), (n-k))))
		enddo
		t2(M) = t2(M) + t1(ii) + gsum1
	enddo
	t2(M) = t2(M) / ((dble(n-k)/dble(l))*2.d0*(1.d0+2.d0*dble(l)-4.d0))

!# Define objective function variables
!var x{i in 1..M} = if i<=M-1 then max(t2[M],A[i])*(t2[i]-0.5)+0.5
!    else t2[M];
	do i = 1,M-1
		x(i) = max(t2(M),A(i))*(t2(i)-0.5d0)+0.5d0
	enddo
	x(M) = t2(M)

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
