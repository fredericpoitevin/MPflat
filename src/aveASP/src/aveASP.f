c       Frederic Poitevin - 2018
c	
	include 'param.h'
c	
	character*80 fparam
c	
	call getarg(1,fparam)
	call readinput(fparam)
	call cubature
	call atomparam
	call readpdb
	call surftens
c
	stop
	end
