c	Vector		Version 1 8/2/1990		Patrice Koehl
c
c	This file contains several subroutine that can be used for any
c	vector operations (vector in 3D cartesian space)
c
c	This includes :	crossvect	: cross vector of two vectors
c			dotvect		: dot product of two vectors
c			normvect	: norm of a vector
c			detvect		: determinant of three vectors
c			diffvect	: substract two vectors
c			addvect		: add two vectors
c
c	For each subroutine : u is for vectors (arrays of size 3)
c			      all other value are scalar
c			      calculations are done in double precision
c
c	1 . crossvect :
c
	subroutine crossvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	u3(1) = u1(2)*u2(3) - u1(3)*u2(2)
	u3(2) = -u1(1)*u2(3) + u1(3)*u2(1)
	u3(3) = u1(1)*u2(2) - u1(2)*u2(1)
c
	return
	end
c
c	2. dotvect :
c
	subroutine dotvect(u1,u2,dot)
c
	integer	i
c
	real*8	u1(3),u2(3),dot
c
	dot = 0.d0
	do 100 i = 1,3
		dot = dot + u1(i)*u2(i)
100	continue
c
	return
	end
c
c	3. normvect :
c
	subroutine normvect(u1,norm)
c
	real*8	u1(3),norm
c
	call dotvect(u1,u1,norm)
	norm = dsqrt(norm)
c
	return
	end
c
c	4. detvect :
c
	subroutine detvect(u1,u2,u3,det)
c
	real*8	u1(3),u2(3),u3(3),det,u4(3)
c
	call crossvect(u2,u3,u4)
	call dotvect(u1,u4,det)
c
	return
	end
c
c	5. diffvect :
c
	subroutine diffvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	integer i
c
	do 100 i = 1,3
		u3(i) = u2(i) - u1(i)
100	continue
c
	return
	end
c
c	6. addvect :
c
	subroutine addvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	integer i
c
	do 100 i = 1,3
		u3(i) = u1(i) + u2(i)
100	continue
c
	return
	end
c
c	7. Normalise a vector : given a vector u1, output u1/norm(u1) :
c
	subroutine unitvector(u1,u2)
c
	real*8  u1(3),u2(3),norm
c
	integer i
c
	call normvect(u1,norm)
c
	do 100 i = 1,3
		u2(i) = u1(i)/norm
100	continue
c
	return
	end
