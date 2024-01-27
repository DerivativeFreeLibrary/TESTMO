!###############################################################################
!#
!#   As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
!#   test problems, linkages, and evolutionary methodologies", GECCO'06}: 
!#   Proceedings of the 8th Annual Conference on Genetic and Evolutionary 
!#   Computation, 1141-1148, 2006.
!#
!#   Example T4, with linkage L1.
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
!param m := 10, >=2;
!param pi := 4*atan(1);
!param A{1..m,1..m};
!
!# x variable
!var x{1..m};
!
!# y variable
!var y{i in 1..m} = sum {j in 1..m} A[i,j] * x[j];
!
!# functions
!var ff1 = y[1]^2;
!
!# g(x)
!var gx = 1 + 10*(m-1) + sum {i in 2..m} (y[i]^2-10*cos(4*pi*y[i]));
!
!var h = 1-sqrt(ff1/gx);
!
!minimize f1:
!	ff1;
!minimize f2:
!    gx*h;
!
!subject to bounds1:
!	0.0 <= x[1] <= 1.0;
!
!subject to bounds {i in 2..m}:
!	-5.0 <= x[i] <= 5.0;
!
!data;
!
!param A :        1            2           3            4           5            6
! :=
!1     1     0           0            0           0            0
!2     0     0.884043   -0.272951    -0.993822    0.511197    -0.0997948
!3     0    -0.492722    0.0646786   -0.666503   -0.945716    -0.334582
!4     0     0.308861   -0.0437502   -0.374203    0.207359    -0.219433
!5     0    -0.708948   -0.37902      0.576578    0.0194674   -0.470262
!6     0    -0.827302    0.669248     0.494475    0.691715    -0.198585
!7     0    -0.715997    0.220772     0.692356    0.646453    -0.401724
!8     0     0.613732   -0.525712    -0.995728    0.389633    -0.064173
!9     0    -0.160446   -0.394585    -0.167581    0.0679849    0.449799
!10    0     0.162711    0.294454    -0.563345   -0.114993     0.549589
!
!:        7             8             9            10        :=
!1     0            0             0             0
!2    -0.659756     0.575496      0.675617      0.180332
!3     0.611894     0.281032      0.508749     -0.0265389
!4     0.914104     0.184408      0.520599     -0.88565
!5     0.572576     0.351245     -0.480477      0.238261
!6     0.0492812    0.959669      0.884086     -0.218632
!7     0.615443    -0.0601957    -0.748176     -0.207987
!8     0.662131    -0.707048     -0.340423      0.60624
!9     0.733505    -0.00918638    0.00446808    0.404396
!10   -0.775141     0.677726      0.610715      0.0850755;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 10
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
	integer	:: n, M, i, j
	real*8		:: x(n), f(M), gx, y(n), h
	real*8, parameter :: pi = 4.d0*atan(1.d0);
	real*8		:: A(10,10)
	save A
	data A /    1.d0  ,   0.d0        ,   0.d0         ,	   0.d0  ,   0.d0          ,  0.d0,	&
		     0.d0  ,   0.d0        ,     0.d0       ,      0.d0	,&
		     0.d0  ,   0.884043d0  , -0.272951d0    ,	-0.993822d0  ,  0.511197d0 ,   -0.0997948d0,	&
		-0.659756d0,     0.575496d0,      0.675617d0,      0.180332d0,	&
		     0.d0  ,  -0.492722d0  ,  0.0646786d0   ,	-0.666503d0  , -0.945716d0 ,   -0.334582d0,	&
		0.611894d0 ,    0.281032d0 ,     0.508749d0 ,    -0.0265389d0,	&
		     0.d0  ,   0.308861d0  , -0.0437502d0   ,	-0.374203d0  ,  0.207359d0 ,   -0.219433d0,	&
		0.914104d0 ,    0.184408d0 ,     0.520599d0 ,    -0.88565d0,	&
		     0.d0  ,  -0.708948d0  , -0.37902d0     , 	 0.576578d0 ,   0.0194674d0,   -0.470262d0,	&
		0.572576d0 ,    0.351245d0 ,    -0.480477d0 ,     0.238261d0,	&
		     0.d0  ,  -0.827302d0  ,  0.669248d0    ,	 0.494475d0  ,  0.691715d0 ,   -0.198585d0,	&
		0.0492812d0,    0.959669d0 ,     0.884086d0 ,    -0.218632d0,	&
		     0.d0  ,  -0.715997d0  ,  0.220772d0    ,	 0.692356d0  ,  0.646453d0 ,   -0.401724d0,	&
		0.615443d0 ,   -0.0601957d0,    -0.748176d0 ,    -0.207987d0,	&
		     0.d0  ,   0.613732d0  , -0.525712d0    ,	-0.995728d0  ,  0.389633d0 ,   -0.064173d0,	&
		0.662131d0 ,   -0.707048d0 ,    -0.340423d0 ,     0.60624d0,	&
		     0.d0  ,  -0.160446d0  , -0.394585d0    ,	-0.167581d0 ,   0.0679849d0,    0.449799d0,	&
		0.733505d0 ,   -0.00918638d0,    0.00446808d0,    0.404396d0,	&
		     0.d0  ,   0.162711d0   , 0.294454d0    ,	-0.563345d0  , -0.114993d0 ,    0.549589d0,	&
		-0.775141d0,     0.677726d0 ,     0.610715d0,    0.0850755d0 /

!# y variable
!var y{i in 1..m} = sum {j in 1..m} A[i,j] * x[j];
	do i = 1,n
		y(i) = 0.d0
		do j = 1,n
			y(i) = y(i) + A(j,i)*x(j)
		enddo
	enddo
!# functions
!var ff1 = y[1]^2;
	f(1) = y(1)**2.d0
!# g(x)
!var gx = 1 + 10*(m-1) + sum {i in 2..m} (y[i]^2-10*cos(4*pi*y[i]));
	gx = 0.d0
	do i = 2,n
		gx = gx + (y(i)**2.d0-10.d0*cos(4.d0*pi*y(i)))
	enddo
	gx = 1.d0 + 10.d0*dble(n-1) + gx
!
!var h = 1-sqrt(ff1/gx);
	h = 1.d0 - sqrt(f(1)/gx)
!minimize f2:
!    gx*h;
	f(2) = gx*h

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -5.d0; 	    ub = 5.d0
	lb(1) = 0.d0;       ub(1) = 1.d0

	return
end subroutine setbounds
