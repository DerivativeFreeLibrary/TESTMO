!###############################################################################
!#
!#   As described by F.Y. Cheng and X.S. Li, "Generalized center method for
!#   multiobjective engineering optimization", Engineering Optimization,31:5,
!#   641-661, 1999.
!#
!#   Example 2, four bar truss.
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
!param F :=  10;
!param E := 2*10^5;
!param L := 200;
!param sigma := 10;
!
!# x variable
!var x{1..4};
!
!
!minimize f1:
!    L*(2*x[1]+sqrt(2)*x[2]+sqrt(x[3])+x[4]);
!minimize f2:
!	F*L/E*(2/x[1]+(2*sqrt(2))/x[2]-(2*sqrt(2))/x[3]+2/x[4]);
!
!
!subject to bound1:
!    F/sigma <= x[1] <= 3*F/sigma;
!subject to bound2:
!    sqrt(2)*F/sigma <= x[2] <= 3*F/sigma;
!subject to bound3:
!    sqrt(2)*F/sigma <= x[3] <= 3*F/sigma;
!subject to bound4:
!    F/sigma <= x[4] <= 3*F/sigma;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 4
	m = 0
	q = 2

	return
end subroutine setdim

subroutine startp(n,x)
	implicit none
	integer	:: n
	real*8		:: x(n), l(n), u(n)

	call setbounds(n,l,u)

	x = (u+l)/2.d0
	
	return
end subroutine startp

subroutine functs(n,x,q,ff)
	implicit none
	integer	:: n,q
	real*8		:: x(n), ff(q)
	real*8, parameter :: F = 10.d0
	real*8, parameter :: E = 2.d+5
	real*8, parameter :: L = 200.d0
	real*8, parameter :: sigma = 10.d0

	ff(1) = L*(2.d0*x(1) + sqrt(2.d0)*x(2) + sqrt(x(3)) + x(4))
	ff(2) = F*L/E*(2.d0/x(1) + (2.d0*sqrt(2.d0))/x(2) - (2.d0*sqrt(2.d0))/x(3) + 2.d0/x(4))

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)
	real*8, parameter :: F = 10.d0
	real*8, parameter :: E = 2.d+5
	real*8, parameter :: L = 200.d0
	real*8, parameter :: sigma = 10.d0

	lb(1) = F/sigma; 	    ub(1) = 3*F/sigma
	lb(2) = sqrt(2.d0)*F/sigma; ub(2) = 3*F/sigma
	lb(3) = sqrt(2.d0)*F/sigma; ub(3) = 3*F/sigma
	lb(4) = F/sigma; 	    ub(4) = 3*F/sigma

	return
end subroutine setbounds
