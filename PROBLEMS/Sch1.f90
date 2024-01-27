!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example Sch1, see the previous cited paper for the original reference.
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
!# x variable
!var x >=0, <= 5;
!
!minimize f1:
!    if x<=1 then -x else
!        if x<= 3 then -2+x else
!            if x<=4 then 4-x else
!                -4+x;
!minimize f2:
!    (x-5)^2;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 1
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
	real*8		:: x(n), f(M)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

!minimize f1:
!    if x<=1 then -x else
!        if x<= 3 then -2+x else
!            if x<=4 then 4-x else
!                -4+x;
	if(x(1) <= 1.d0) then
		f(1) = -x(1)
	elseif(x(1) <= 3.d0) then
		f(1) = -2.d0 + x(1)
	elseif(x(1) <= 4.d0) then
		f(1) = 4.d0 - x(1)
	else
		f(1) = -4.d0 + x(1)
	endif

!minimize f2:
!    (x-5)^2;
	f(2) = (x(1)-5.d0)**2.d0

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

	lb = 0.d0;	    ub = 5.d0

	return
end subroutine setbounds
