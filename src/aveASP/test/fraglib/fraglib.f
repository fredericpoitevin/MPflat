c	
c	This program fragments "table.out" in as many files
c	as needed
c	
	integer flag,natom,natomtmp
c	
	character record*20,fnew*7
c	
1	format(a)
c	
	flag = 0
	open(unit=1,file="table.out",status='unknown')
100	  read(1,1,end=200) record
	  if(flag.eq.0.and.record(1:3).eq."...") then
	    flag = 1
	    natomtmp = 0
	    goto 100
	  endif
	  if(flag.eq.1) then
	    fnew = record(1:3)//".lib"
	    read(record(4:7),*) natom
	    flag = 2
	  endif
	  if(flag.eq.2.and.natomtmp.le.natom) then
	    natomtmp = natomtmp + 1
	    open(unit=2,file=fnew,access='append',status='unknown')
	      write(2,1) record
	    close(unit=2)
	    if(natomtmp.eq.natom+1) flag = 0
	    if(natom.eq.0) flag = 0
	  endif
	  goto 100
200	  continue
	close(unit=1)
c	
	stop
	end
