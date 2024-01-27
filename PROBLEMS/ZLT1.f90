!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example ZLT1, see the previous cited paper for the original reference.
!#
!#   In the above paper the number of variables was set equal to 100. 
!#   We selected n=10 as default.
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
!# parameters for problem dimension
!param n := 10;
!param M := 3;
!
!# x variable
!var x{1..n} >= -1000, <=1000;
!
!
!minimize f{m in 1..M}:
!    (x[m]-1)^2 + sum {i in 1..n} (if i<>m then x[i]^2 else 0);


subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 10
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

subroutine functs(n,x,M,f)
	implicit none
	integer	:: n, M, i, k
	real*8		:: x(n), f(M)


!minimize f{m in 1..M}:
!    (x[m]-1)^2 + sum {i in 1..n} (if i<>m then x[i]^2 else 0);


do k=1,M
    f(k)=(x(k)-1)**2

    do i=1, n
	   if (i/=k) then 
	       f(k)=f(k)+x(i)**2
	   end if
    end do
	
end do


return
end subroutine functs

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n
	real*8		:: lb(n), ub(n)

	lb = -1000.0d0;	    ub = 1000.0d0

	return
end subroutine setbounds
