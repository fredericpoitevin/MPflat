c
c	
	subroutine writeprof
c	
	include'param.h'
c	
	integer ulog,isample,nsample
	integer flagsolv,flagfit,flagmeth,fsolvfmt
c	
	real*8 qrange,qstep,profil(nqmax),expq(nqmax)
        real*8 sigexp(nqmax),profexp(nqmax)
	real*8 bulkdens,cbest,c1best,c2best,c1,c2,chisqmin
	real*8 c1tmp
	real*8 redres
c	
	character*64 fpdb,fasexp,fsolvent,fprof,flog
c	
1	format(e12.7,1(1x,e12.7))
2	format(e12.7,3(1x,e12.7),1x,f12.7)
3	format("# ... Here is your AquaSaxs profile ...")
4	format("# PDB model: ",a30)
5	format("# I exp:     ",a30)
6	format("# Fitting parameters: ")
7	format("#                 C1: ",f6.3)
8	format("#                 C2: ",f6.3)
9	format("#                 C : ",e15.7)
11	format("#                Chi: ",e15.7)
12	format("#   q(A-1)   .   Icalc    .   Iexp     .   sigma    . resid")
13	format("#   q(A-1)   .   Icalc   ")
800	format(">> Write computed SAXS profile in: ",a30)
801	format("# Outputs #")
c	
	common /files/ fpdb,fasexp,fsolvent,fprof,flog,ulog
	common /profinfo/ qrange,qstep,profil,expq,
     &                    profexp,sigexp,nsample
	common /fitparam/ c1best,c2best,cbest,chisqmin
	common /flags/ flagsolv,flagfit,flagmeth,fsolvfmt
	common /add/ c1tmp
	common /params/ bulkdens
c	
	open(unit=ulog,file=flog,access='append',status='unknown')
	  write(ulog,*) " "
c	  write(ulog,801)
	  write(ulog,*) " "
c	
	c1 = c1tmp
	c2 = c2best/bulkdens
	open(unit=1,file=fprof,status='unknown')
	  write(1,3)
	  write(1,4) fpdb
	  if(fasexp.ne.'') write(1,5) fasexp
	  if(flagfit.eq.1) then
	    write(1,6)
	    write(1,7) c1
	    write(1,8) c2
	    write(1,9) cbest
	    write(1,11) chisqmin
	    write(1,12)
	  else
	    write(1,13)
	  endif
	  do 10 isample = 1,nsample
	    if(fasexp.eq.'') then
	      write(1,1) expq(isample),profil(isample)
	    else
	      redres = profexp(isample)-profil(isample)
	      redres = redres/sigexp(isample)
	      write(1,2) expq(isample),profil(isample),
     1	            profexp(isample),sigexp(isample),redres
	    endif
10	  continue
	close(unit=1)
c	
	close(unit=ulog)
c	
	return
	end
