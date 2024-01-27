!###############################################################################
!#
!#   As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
!#   Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary 
!#   Computation 8(2): 173-195, 2000.
!#
!#   Example T1.
!#
!#   This file is part of a collection of problems developed for
!#   derivative-free multiobjective optimization in
!#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
!#   Direct Multisearch for Multiobjective Optimization, 2010.
!#
!#   Written by the authors in June 1, 2010.
!#
!###############################################################################

!# parameters
!param m := 30, >=2;

!# x variable
!var x{1..m};

!# functions
!var ff1 = x[1];

!# g(x)
!var gx = 1 + 9/(m-1) * sum {i in 2..m} (x[i]);

!var h = 1-sqrt(ff1/gx);

!minimize f1:
!	ff1;
!minimize f2:
!    gx*h;

!subject to bounds {i in 1..m}:
!	0.0 <= x[i] <= 1.0;





subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 30
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
	real*8		:: x(n), f(M), h
    real*8     :: pi



 pi = 4.d0*atan(1.d0)


!# functions
!var ff1 = x[1];

!# g(x)
!var gx = 1 + 9/(m-1) * sum {i in 2..m} (x[i]);

!var h = 1-(ff1/gx)^2;

!minimize f1:
!	ff1;
!minimize f2:
!    gx*h;


f=0.0d0

f(1)= x(1)


f(2)= 0.0d0 

do i=2, n
   f(2)= f(2)+ x(i)
end do
   
f(2)= 1.0d0 + 9.0d0/dble(n-1) * f(2)

h= 1.0d0 - dsqrt(f(1)/f(2))

f(2) = f(2)*h

return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = 0.0d0;	    ub = 1.0d0

	return
end subroutine setbounds

