!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#           
!#   Example WFG3.
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
!# WFG3
!param A {i in 1..M-1} := if i<=1 then 1 else 0;
!
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
!param AA := 2;
!var t2{i in 1..k+l/2} = if i<=k then t1[i]
!    else sum {ii in (k+2*(i-k)-1)..(k+2*(i-k))} (t1[ii]+sum {jj in 0..AA-2} abs(t1[ii]-t1[(k+2*(i-k)-1)+((ii+jj-(k+2*(i-k)-1)+1) mod (k+2*(i-k)-(k+2*(i-k)-1)+1))]))/((k+2*(i-k)-(k+2*(i-k)-1)+1)/AA*ceil(AA/2)*(1+2*AA-2*ceil(AA/2)));
!
!# third level mapping
!param w{i in 1..n} := 1;
!var t3{i in 1..M} = if i<=M-1 then (sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (w[j]*t2[j]))/(sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} w[j])
!    else (sum {j in k+1..k+l/2} (w[j]*t2[j]))/(sum {j in k+1..k+l/2} w[j]);
!
!# Define objective function variables
!var x{i in 1..M} = if i<=M-1 then max(t3[M],A[i])*(t3[i]-0.5)+0.5
!    else t3[M];
!
!# Define objective function function h
!var h{m in 1..M} = if m==1 then prod {i in 1..M-1} x[i]
!    else if m<=M-1 then (prod {i in 1..M-m} x[i])*(1-x[M-m+1])
!        else 1-x[1];
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
	real*8		:: t1(8), t2(6), t3(3), w(8), t4(3), x(3)
!	real*8		:: S(M), zmax(n), A(M-1), y(n), h(M)
!	real*8		:: t1(n), t2(k+l/2), t3(M), w(n), t4(M), x(M)
	real*8		:: gsum1, gsum2
	integer, parameter :: AA = 2;
	real*8, parameter :: alpha = 1.d0
	real*8, parameter :: beta = 1.d0
	real*8, parameter :: AAAA = 5.d0

!param S {m in 1..M} := 2*m;
	do i = 1,M
		S(i) = 2.d0*dble(i)
	enddo

!param zmax {i in 1..n} := 2*i;
	do i = 1,n
		zmax(i) = 2.d0*dble(i)
	enddo

!param A {i in 1..M-1} := if i<=1 then 1 else 0;
	A(1) = 1.d0
	do i = 2,M-1
		A(i) = 0.d0
	enddo

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
!var t2{i in 1..k+l/2} = if i<=k then t1[i]
!    else sum {ii in (k+2*(i-k)-1)..(k+2*(i-k))} (t1[ii]+sum {jj in 0..AA-2} abs(t1[ii]-t1[(k+2*(i-k)-1)+((ii+jj-(k+2*(i-k)-1)+1) mod (k+2*(i-k)-(k+2*(i-k)-1)+1))]))/
!						((k+2*(i-k)-(k+2*(i-k)-1)+1)/AA*ceil(AA/2)*(1+2*AA-2*ceil(AA/2)));    
	do i = 1,k
		t2(i) = t1(i)
	enddo
	do i = k+1,k+l/2
		gsum2 = 0.d0
		do ii = (k+2*(i-k)-1), (k+2*(i-k))
			gsum1 = 0.d0
			do jj = 0, AA-2
				gsum1 = gsum1 + abs(t1(ii)-t1((k+2*(i-k)-1)+mod(((ii+jj-(k+2*(i-k)-1)+1)),(k+2*(i-k)-(k+2*(i-k)-1)+1))))
			enddo
			gsum2 = gsum2 + t1(ii) + gsum1
		enddo
		t2(i) = gsum2 / (dble(k+2*(i-k)-(k+2*(i-k)-1)+1)/dble(AA)*(AA/2)*(1.d0+2.d0*dble(AA)-2.d0*(AA/2))) 
	enddo

!# third level mapping
!param w{i in 1..n} := 1;
!var t3{i in 1..M} = if i<=M-1 then (sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (w[j]*t2[j]))/(sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} w[j])
!    else (sum {j in k+1..k+l/2} (w[j]*t2[j]))/(sum {j in k+1..k+l/2} w[j]);
	do i = 1,n
		w(i) = 1.d0
	enddo
	do i = 1,M-1
		gsum1 = 0.d0
		gsum2 = 0.d0
		do j = ((i-1)*k/(M-1)+1), (i*k/(M-1))
			gsum1 = gsum1 + w(j)*t2(j)
			gsum2 = gsum2 + w(j)
		enddo
		t3(i) = gsum1 / gsum2
	enddo
	gsum1 = 0.d0
	gsum2 = 0.d0
	do j = k+1,k+l/2
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
!var h{m in 1..M} = if m==1 then prod {i in 1..M-1} x[i]
!    else if m<=M-1 then (prod {i in 1..M-m} x[i])*(1-x[M-m+1])
!        else 1-x[1];
	h(1) = 1.d0
	do i = 1,M-1
		h(1) = h(1)*x(i)
	enddo
	do j = 2,M-1
		h(j) = 1.d0
		do i = 1,M-j
			h(j) = h(j)*x(i)
		enddo
		h(j) = h(j)*(1.d0-x(M-j+1))
	enddo
	h(M) = 1.d0-x(1)

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
