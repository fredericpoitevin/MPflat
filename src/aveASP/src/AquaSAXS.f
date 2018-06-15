c	
c	!!!!!!!!!!!!!!!!!!!!!! AQUASAXS v2.2 !!!!!!!!!!!!!!!!!!!!!!!!
c	
c	Frederic Poitevin	November, 10th 2010
c	July 2011: upgrade (original)
c	August 2011: new handling for user-defined atomic types.
c	
	include 'param.h'
c	
	character*80 fparam
c	
c	-------------------------------------------------------------
c		Parameters					
c	-------------------------------------------------------------
	call getarg(1,fparam)
	call readinput(fparam)
	call cubature
	call atomparam
c	-------------------------------------------------------------
c		Read and store input data			
c	-------------------------------------------------------------
	call readpdb
	call profilexp
c	-------------------------------------------------------------
c		Compute configuration independent terms (fi,gi)	
c	-------------------------------------------------------------
	call atomicsf
c	-------------------------------------------------------------
c		Add solvation layer and compute profile		
c	-------------------------------------------------------------
	call solvterm
	call saxs
	call writeprof
c	-------------------------------------------------------------
c		End of the program				
c	-------------------------------------------------------------
	stop
	end
