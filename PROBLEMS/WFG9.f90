!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#           
!#   Example WFG9.
!#
!#   This file is part of a collection of problems developed for
!#   derivative-free multiobjective optimization in
!#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
!#   Direct Multisearch for Multiobjective Optimization, 2010.
!#
!#   Written by the authors in June 1, 2010.
!#
!###############################################################################

!param M := 3;
!param k := 4;
!param l := 4;
!param n := k+l;

!param pi := 4*atan(1);
!param pi2:= 2*atan(1);

!param S {m in 1..M} := 2*m;

!# neq WFG3
!param A {i in 1..M-1} := 1;


!# problem variables
!param zmax {i in 1..n} := 2*i;
!var z{i in 1..n} >=0, <= zmax[i];

!# transform z into [0,1] set
!var y{i in 1..n} = z[i]/zmax[i];

!# first level mapping
!param w{1..n}, default 1.0;
!param AA := 0.98/49.98;
!param BB := 0.02;
!param CC := 50;
!var r_sum{i in 1..n-1} = (sum {j in i+1..n} (w[j]*y[j]))/(sum {j in i+1..n} w[j]);
!var t1{i in 1..n} = if i<=n-1 then y[i]^(BB+(CC-BB)*(AA-(1-2*r_sum[i])*abs(floor(0.5-r_sum[i])+AA)))
!    else y[i];

!# second level mapping
!param AAA := 0.35;
!param BBB := 0.001;
!param CCC := 0.05;
!param AAAA := 30;
!param BBBB := 95;
!param CCCC := 0.35;
!var t2{i in 1..n} = if i<=k then 1+(abs(t1[i]-AAA)-BBB)*((floor(t1[i]-AAA+BBB)*(1-CCC+(AAA-BBB)/BBB))/(AAA-BBB)+(floor(AAA+BBB-t1[i])*(1-CCC+(1-AAA-BBB)/BBB))/(1-AAA-BBB)+1/BBB)
!    else (1+cos((4*AAAA+2)*pi*(0.5-abs(t1[i]-CCCC)/(2*(floor(CCCC-t1[i])+CCCC))))+4*BBBB*(abs(t1[i]-CCCC)/(2*(floor(CCCC-t1[i])+CCCC)))^2)/(BBBB+2);

!# third level mapping
!var t3{i in 1..M} = if i<=M-1 then sum {ii in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (t2[ii]+sum {jj in 0..(k/(M-1)-2)} abs(t2[ii]-t2[((i-1)*k/(M-1)+1)+((ii+jj-((i-1)*k/(M-1)+1)+1) mod ((i*k/(M-1))-((i-1)*k/(M-1)+1)+1))]))/(((i*k/(M-1))-((i-1)*k/(M-1)+1)+1)/(k/(M-1))*ceil(k/(M-1)/2)*(1+2*k/(M-1)-2*ceil(k/(M-1)/2)))
!    else sum {ii in k+1..n} (t2[ii]+sum {jj in 0..(l-2)} abs(t2[ii]-t2[k+1+((ii+jj-(k+1)+1) mod (n-k))]))/(((n-k)/l)*ceil(l/2)*(1+2*l-2*ceil(l/2)));

!# Define objective function variables
!var x{i in 1..M} = if i<=M-1 then max(t3[M],A[i])*(t3[i]-0.5)+0.5
!    else t3[M];

!# Define objective function function h
!var h{m in 1..M} = if m==1 then prod {i in 1..M-1} sin(x[i]*pi2)
!    else if m<=M-1 then (prod {i in 1..M-m} sin(x[i]*pi2))*cos(x[M-m+1]*pi2)
!        else cos(x[1]*pi2);


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
	real*8		:: S(3), zmax(8), A(2), y(8), h(3), rsum(4)
	real*8		:: t1(8), t2(8), t3(3), w(8), t4(3), x(3)
!	real*8		:: S(M), zmax(n), A(M-1), y(n), h(M), rsum(k)
!	real*8		:: t1(n), t2(n), t3(M), w(n), t4(M), x(M)
	real*8		:: gsum1, gsum2
	real*8, parameter :: AA = 0.98d0/49.98d0;
	real*8, parameter :: BB = 0.02d0;
	real*8, parameter :: CC = 50.d0;

	!param AAA := 0.35;
    !param BBB := 0.001;
    !param CCC := 0.05;
    !param AAAA := 30;
    !param BBBB := 95;
    !param CCCC := 0.35;

	real*8, parameter :: AAA = 0.35d0;
	real*8, parameter :: BBB = 0.001d0;
	real*8, parameter :: CCC = 0.05d0;
	real*8, parameter :: AAAA = 30.d0
	real*8, parameter :: BBBB = 95.d0
	real*8, parameter :: CCCC = 0.35d0

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

!param w{1..n}, default 1.0;
    do i = 1,n
		w(i) = 1.d0
	enddo


!# first level mapping
!var r_sum{i in 1..n-1} = (sum {j in i+1..n} (w[j]*y[j]))/(sum {j in i+1..n} w[j]);
!var t1{i in 1..n} = if i<=n-1 then y[i]^(BB+(CC-BB)*(AA-(1-2*r_sum[i])*abs(floor(0.5-r_sum[i])+AA)))
!    else y[i];
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
!var t2{i in 1..n} = if i<=k then 
!1+(abs(t1[i]-AAA)-BBB)*((floor(t1[i]-AAA+BBB)*(1-CCC+(AAA-BBB)/BBB))/(AAA-BBB)+(floor(AAA+BBB-t1[i])*(1-CCC+(1-AAA-BBB)/BBB))/(1-AAA-BBB)+1/BBB)
!else 
! (1+cos((4*AAAA+2)*pi*(0.5-abs(t1[i]-CCCC)/(2*(floor(CCCC-t1[i])+CCCC))))+4*BBBB*(abs(t1[i]-CCCC)/(2*(floor(CCCC-t1[i])+CCCC)))^2)/(BBBB+2);
	do i = 1,k
           t2(i) = 1+(dabs(t1(i)-AAA)-BBB)*((floor(t1(i)-AAA+BBB)*(1.0d0-CCC+(AAA-BBB)/BBB))/(AAA-BBB))
           t2(i) = t2(i)+(dabs(t1(i)-AAA)-BBB)*((floor(AAA+BBB-t1(i))*(1.0d0-CCC+(1.0d0-AAA-BBB)/BBB))/(1.0d0-AAA-BBB)+1.0d0/BBB)
	enddo
	do i = k+1,n
		t2(i) = 1.0d0+dcos((4.0d0*AAAA+2.0d0)*pi*(0.50d0-dabs(t1(i)-CCCC)/(2.0d0*(floor(CCCC-t1(i))+CCCC))))
	        t2(i) = t2(i)+4.0d0*BBBB*(dabs(t1(i)-CCCC)/(2.0d0*(floor(CCCC-t1(i))+CCCC)))**2
                t2(i) = t2(i)/(BBBB+2.0d0)
        enddo

!# third level mapping
!# w already defined

!# third level mapping
!var t3{i in 1..M} = if i<=M-1 
! then 
! sum {ii in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (t2[ii]+sum {jj in 0..(k/(M-1)-2)} abs(t2[ii]-t2[((i-1)*k/(M-1)+1)+((ii+jj-((i-1)*k/(M-1)+1)+1) mod ((i*k/(M-1))-((i-1)*k/(M-1)+1)+1))]))/(((i*k/(M-1))-((i-1)*k/(M-1)+1)+1)/(k/(M-1))*ceil(k/(M-1)/2)*(1+2*k/(M-1)-2*ceil(k/(M-1)/2)))
! else 
! sum {ii in k+1..n} (t2[ii]+sum {jj in 0..(l-2)} abs(t2[ii]-t2[k+1+((ii+jj-(k+1)+1) mod (n-k))] ) )/(((n-k)/l)*ceil(l/2)*(1+2*l-2*ceil(l/2)));

	do i = 1,M-1
		gsum1 = 0.d0
		do ii = ((i-1)*k/(M-1)+1),(i*k/(M-1))
			
		   gsum1=gsum1+t2(ii)

 		   do jj = 0, (k/(M-1)-2)
		     gsum1= gsum1 +  dabs(t2(ii)-t2( ( (i-1)*k/(M-1)+1 ) +mod((ii+jj-((i-1)*k/(M-1)+1)+1),((i*k/(M-1))-((i-1)*k/(M-1)+1)+1))))
	           end do
			
		
		enddo
                 
      gsum1=gsum1/(((i*k/(M-1))-((i-1)*k/(M-1)+1)+1)/(k/(M-1))*ceiling(dble(k/(M-1)/2))*(1+2*k/(M-1)-2*ceiling(dble(k/(M-1)/2))))

	
	
		t3(i) = gsum1
	enddo


	gsum1 = 0.d0
	do ii = k+1,n
		gsum1=gsum1+t2(ii)

 			do jj = 0, l-2
		               gsum1= gsum1 +  dabs(t2(ii)-t2( k+1+mod((ii+jj-(k+1)+1), (n-k)) ) )
			end do	
		
	enddo
	gsum1=gsum1/(((n-k)/l)*ceiling(dble(l/2))*(1+2*l-2*ceiling(dble(l/2))))
        
	t3(M) = gsum1 


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
