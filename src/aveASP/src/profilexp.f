c	
c	30/06/10	FP
c	This subroutines reads and store a .dat SAXS profile file
c	
	subroutine profilexp
c	
	include 'param.h'
c	
	integer i,j,k,isample,nsample
	integer ulog,ilog,ulogtmp
	integer imin,imax,nsample2
	integer ntest,nsigma
	integer iline1,iline2,iline3
	integer flagunit
	integer flagauth
	integer lenline,flag,icol,indx(200),istart(3),iend(3)
	integer i1,i2,i3,i4,i5,i6
c	
	real*8 qrange,qstep,profil(nqmax),expq(nqmax)
	real*8 sigexp(nqmax),profexp(nqmax)
	real*8 sigexptmp(nqmax),profexptmp(nqmax),expqtmp(nqmax)
	real*8 diff1,diff2,pi
	real*8 moy
	real*8 thresh
c	
	character record*200,line1*16,line2*16,line3*16
	character inum1*1,auth(16)*1,testauth*1,test*1
	character*64 fpdb,fasexp,fsolvent,fprof,flog
c	
	common /files/ fpdb,fasexp,fsolvent,fprof,flog,ulog
	common /profinfo/ qrange,qstep,profil,expq,
     &	                  profexp,sigexp,nsample
	common /qunit/ flagunit
c	
1	format(a)
2	format("# Reading experimental profile #")
3	format(">> Number of computed points: ",i4)
4	format(">>   ranging from: ",f6.3," to: ",f6.3," A-1")
5	format("ERR there is a format problem here: ",a50)
6	format("ERR Program stopped. File needs to be edited.")
7	format("ERR at line ",i4,", position ",i4)
c	
	data (auth(j),j=1,16) /'1','2','3','4','5','6','7','8','9',
     1  '.','E','e','-','+',' ','0'/
c	
	pi = acos(-1.d0)
	thresh = 1.d-4
	flagauth = 1
c	
	do 301 ilog = 1,2
	  if(ilog.eq.1) then
	    ulogtmp = 6
	  else
	    ulogtmp = ulog
	open(unit=ulogtmp,file=flog,access='append',status='unknown')
	  endif
c	
	if(fasexp.ne.'') then
	  write(ulogtmp,*) ' '
	  write(ulogtmp,*) ' '
	  write(ulogtmp,2)
	  write(ulogtmp,*) ' '
	endif
301	continue
c	
	  nsample = int(qrange/qstep)+1
	  do 100 isample = 1,nsample
	    expq(isample) = (isample-1)*qstep
100	  continue
	if(fasexp.ne.'') then
	  open(unit=3,file="filtered.dat",status='unknown')
	  nsample2 = 0
	  open(unit=4,file=fasexp,status='unknown')
200	    read(4,1,end=300) record
c	
	    lenline = len_trim(record)
	    do 201 i = 1,lenline
	      test = record(i:i)
	      flag = 0
	      do 202 j = 1,16
	        if(test.eq.auth(j)) then
	          flag = 1
	        endif
202	      continue
	      if(flag.eq.0) goto 200
201	    continue
c	
	    nsample2 = nsample2 + 1
	    write(3,1) record
	    icol = 0
	    do 203 i = 1,lenline
	      test = record(i:i)
	      if(test.ne.' ') then
	        if(i.gt.1) then
	          if(indx(i-1).eq.0) then
	            icol = icol + 1
	            indx(i) = icol
	            istart(icol) = i
	          else
	            indx(i) = icol
	          endif
	        else
	          icol = icol + 1
	          indx(i) = icol
	          istart(icol) = i
	        endif
	      else
	        indx(i) = 0
	        if(i.gt.1) then
	          if(indx(i-1).gt.0) then
	            iend(icol) = i-1
	          endif
	        endif
	      endif
203	    continue
c	
	    i1 = istart(1)
	    i2 = iend(1)
	    i3 = istart(2)
	    if(icol.eq.2) then
	      nsigma = 0
	      i4 = lenline
	    else
	      nsigma = 1
	      i4 = iend(2)
	      i5 = istart(3)
	      i6 = lenline
	    endif
c	
	    read(record(i1:i2),*) expqtmp(nsample2)
	    read(record(i3:i4),*) profexptmp(nsample2)
	    if(nsigma.eq.1) then
	      read(record(i5:i6),*) sigexptmp(nsample2)
	    endif
c	
	    if(flagunit.eq.2) then
	      expqtmp(nsample2) = 1.d1*expqtmp(nsample2)
	    elseif(flagunit.eq.3) then
	      expqtmp(nsample2) = 2*pi*expqtmp(nsample2)
	    elseif(flagunit.eq.4) then
	      expqtmp(nsample2) = 2.d1*pi*expqtmp(nsample2)
	    endif
c	
	    goto 200
300	    continue
	  close(unit=4)
	  close(unit=3)
c	
	  if(nsample2.eq.0) then
	    open(unit=5,file="data.err",status='unknown')
	      write(5,*) "No data was recognized in the data file."
	    close(unit=5)
	    stop
	  endif
c	
	if(nsigma.eq.0) then
	  do 500 isample = 1,nsample2
	    if(isample.eq.1) then
	      moy = (profexptmp(1)+profexptmp(2))/2
	      sigexptmp(1) = sqrt((moy-profexptmp(1))**2)
	    elseif(isample.eq.nsample2) then
	      moy = (profexptmp(nsample2-1)+profexptmp(nsample2))/2
	      sigexptmp(isample) = sqrt((moy-profexptmp(isample))**2)
	    else
	      moy = (profexptmp(isample-1)+profexptmp(isample+1))/2
	      sigexptmp(isample) = sqrt((moy-profexptmp(isample))**2)
	    endif
	    if(sigexptmp(isample).lt.thresh) then
	       sigexptmp(isample) = thresh
	    endif
cc	    write(6,*) sigexptmp(isample)
cc	    sigexptmp(isample) = sigexptmp(isample)*10.d2
500	  continue
	endif
c	
	  imin = 0
	  imax = nsample
	  do 10 isample = 1,nsample
	    if(expq(isample).lt.expqtmp(1)) then
	       imin = imin + 1
	    elseif(expq(isample).gt.expqtmp(nsample2)) then
	       imax = isample - 1
	       goto 11
	    endif
10	  continue
11	  continue
	  nsample = imax-imin
	  do 302 ilog = 1,2
	    if(ilog.eq.1) then
	      ulogtmp = 6
	    else
	      ulogtmp = ulog
	    endif
	    write(ulogtmp,3) nsample
	    write(ulogtmp,4) expq(1+imin),expq(nsample)
302	  continue
	  do 20 isample = 1,nsample
	    expq(isample) = expq(isample+imin)
	    do 21 i = 1,nsample2
	      if(expqtmp(i).le.expq(isample).and.
     1	         expqtmp(i+1).ge.expq(isample)) then
	        diff1 = expq(isample) - expqtmp(i)
	        diff2 = expqtmp(i+1) - expq(isample)
	        if(diff1.lt.diff2) then
	          profexp(isample) = profexptmp(i)
	          sigexp(isample) = sigexptmp(i)
	          expq(isample) = expqtmp(i)
	        else
	          profexp(isample) = profexptmp(i+1)
                  sigexp(isample) = sigexptmp(i+1)
                  expq(isample) = expqtmp(i+1)
	        endif
c	 write(6,*) expq(isample),profexp(isample),sigexp(isample)
	        goto 20
	      endif
21	    continue
20	  continue
	endif
c	
999	continue
	close(unit=ulog)
c	
	return
	end
