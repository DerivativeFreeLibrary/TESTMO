!# ex005.mod
!# Original AMPL coding by Sven Leyffer, Argonne Natl. Lab.
!#
!# A simple multi-objective optimization problem (p. 281):
!# C.-L. Hwang & A. S. Md. Masud, Multiple Objective
!# Decision Making - Methods and Applications, No. 164 in 
!# Lecture Notes in Economics and Mathematical Systems,
!# Springer, 1979.
!
!# ... variables
!var x1 >= -1, <= 2;
!var x2 >=  1, <= 2;
!
!# ... objective functions
!minimize  f1: x1^2 - x2^2;
!minimize  f2: x1/x2;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 2
	m = 0
	q = 2

	return
end subroutine setdim

subroutine startp(n,x)
	implicit none
	integer	:: n
	real*8		:: x(n), l(n), u(n)

	call setbounds(n,l,u)

	x = (l+u)/2.d0
	
	return
end subroutine startp

subroutine functs(n,x,q,f)
	implicit none
	integer	:: n, q
	real*8		:: x(n), f(q)

	f(1) = x(1)**2.d0 - x(2)**2.d0
	f(2) = x(1)/x(2)

	return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb(1) = -1.d0; 	    ub(1) = 2.d0
	lb(2) =  1.d0; 	    ub(2) = 2.d0

	return
end subroutine setbounds
