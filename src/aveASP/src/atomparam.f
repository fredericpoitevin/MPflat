c	
c	This subroutine contains data block with useful info
c	on atom types...
c	
c	We bound H to their heavy atoms...
c	
c	For now, we do not differentiate protonated and deprotonated
c	states.
c	
	subroutine atomparam
c	
	include 'param.h'
c	
	integer i,j,inum,iline1,iline2
	integer iatom(nresdef),atomtype(nresdef,24)
	integer iaddatom(naddrmax),atomaddtype(naddrmax,100)
	integer ntype
	integer naddatyp,naddrtyp
	integer ulog,ulogtmp,ilog
	integer addHb(naddamax,1)
c	
	real*8 asf(sftyp,4),bsf(sftyp,4),csf(sftyp,1),Xvol(sftyp,1)
        real*8 Xrad(sftyp,1)
	real*8 addasf(naddamax,4),addbsf(naddamax,4)
	real*8 addcsf(naddamax,1),addXvol(naddamax,1)
	real*8 addXrad(naddamax,1)
c	
	character nameres(nresdef)*3,nameatom(nresdef,24)*3
	character*64 faddrestype,faddatomtype
	character record*100,nameaddat(naddamax)*4
	character nameaddres(naddrmax)*3,inum1*1,gtest*2,line*7
	character nameaddatom(naddrmax,100)*4,name2num(14)*4
	character*64 fpdb,fasexp,fsolvent,fprof,flog
c	
1	format(a)
2       format("> Number of newly defined atomic types : ",i5)
3       format("> Number of newly defined residue types: ",i5)
4       format(">> ",a4)
5	format(">> ",a3)
c	
	common /sfparam/ asf,bsf,csf,Xvol,Xrad
	common /addfile/ faddrestype,faddatomtype
	common /addtype/ naddrtyp,naddatyp,nameaddres,iaddatom,
     1	nameaddatom,atomaddtype
	common /addsfparam/ addasf,addbsf,addcsf,addXvol,addXrad,addHb
	common /files/ fpdb,fasexp,fsolvent,fprof,flog,ulog
c	
c	atomtypes:
c	1	C		9	O
c	2	CH		10	OH
c	3	CH2		11	S
c	4	CH3		12	SH
c	5	N		13	P
c	6	NH		14	H
c	7	NH2
c	8	NH3
c	
	data (name2num(i),i=1,14)/'C  ','CH ','CH2','CH3','N  ',
     1	'NH ','NH2','NH3','O  ','OH ','S  ','SH ','P  ','H  '/
c	
c       Hydrogen        
        data (asf(14,i),i=1,4) /0.493002,0.322912,0.140191,0.040810/
        data (bsf(14,i),i=1,4) /10.510900,26.125700,3.142360,57.799698/
        data (csf(14,i),i=1,1) /0.003038/
        data (Xvol(14,i),i=1,1) /5.150000/
        data (Xrad(14,i),i=1,1) /1.07/
c       
c       Carbon          
        data (asf(1,i),i=1,4) /2.310000,1.020000,1.588600,0.865000/
        data (bsf(1,i),i=1,4) /20.843899,10.207500,0.568700,51.651199/
        data (csf(1,i),i=1,1) /0.215600/
        data (Xvol(1,i),i=1,1) /16.44000/
        data (Xrad(1,i),i=1,1) /1.58/
c       
c	CH1
	data (Xvol(2,i),i=1,1) /21.59000/
	data (Xrad(2,i),i=1,1) /1.73/
c	
c	CH2
	data (Xvol(3,i),i=1,1) /26.74/
	data (Xrad(3,i),i=1,1) /1.85/
c	
c	CH3
	data (Xvol(4,i),i=1,1) /31.89/
	data (Xrad(4,i),i=1,1) /1.97/
c	
c       Nitrogen        
        data (asf(5,i),i=1,4) /12.212600,3.132200,2.012500,1.166300/
        data (bsf(5,i),i=1,4) /0.005700,9.893300,28.997499,0.582600/
        data (csf(5,i),i=1,1) /-11.528999/
        data (Xvol(5,i),i=1,1) /2.49000000/
        data (Xrad(5,i),i=1,1) /0.84/
c	
c	NH
	data (Xvol(6,i),i=1,1) /7.64/
	data (Xrad(6,i),i=1,1) /1.22/
c	
c	NH2
	data (Xvol(7,i),i=1,1) /12.79/
	data (Xrad(7,i),i=1,1) /1.45/
c	
c	NH3
	data (Xvol(8,i),i=1,1) /17.94/
	data (Xrad(8,i),i=1,1) /1.62/
c       
c       Oxygen          
        data (asf(9,i),i=1,4) /3.048500,2.286800,1.546300,0.867000/
        data (bsf(9,i),i=1,4) /13.277100,5.701100,0.323900,32.908897/
        data (csf(9,i),i=1,1) /0.250800/
        data (Xvol(9,i),i=1,1) /9.130000/
        data (Xrad(9,i),i=1,1) /1.30/
c	
c	OH
	data (Xvol(10,i),i=1,1) /14.28/
	data (Xrad(10,i),i=1,1) /1.50/
c       
c       Sulfur          
        data (asf(11,i),i=1,4) /6.905300,5.203400,1.437900,1.586300/
        data (bsf(11,i),i=1,4) /1.467900,22.215099,0.253600,56.172001/
        data (csf(11,i),i=1,1) /0.866900/
        data (Xvol(11,i),i=1,1) /19.86000/
        data (Xrad(11,i),i=1,1) /1.68/
c	
c	SH
	data (Xvol(12,i),i=1,1) /25.10/
	data (Xrad(12,i),i=1,1) /1.81/
c       
c       Phosphorus      
        data (asf(13,i),i=1,4) /6.434500,4.179100,1.780000,1.490800/
        data (bsf(13,i),i=1,4) /1.906700,27.157000,0.526000,68.164497/
        data (csf(13,i),i=1,1) /1.114900/
        data (Xvol(13,i),i=1,1) /5.730000/
        data (Xrad(13,i),i=1,1) /1.11/
c	
c	-------------------------------------------------------------
c		Any user-defined residue or atom types?		
c	-------------------------------------------------------------
c	
	naddatyp = 0
	if(faddatomtype.ne.'') then
	  open(unit=1,file=faddatomtype,status='unknown')
100	    read(1,1,end=200) record
	    if(record(1:1).eq."#") goto 100
	    if(record(1:3).eq."...") then
	      naddatyp = naddatyp +1
	      goto 100
	    endif
	    if(record(1:4).eq."NAME") then
	      i = 4
101	      i = i+1
	      if(record(i:i).eq.' ') goto 101
	      nameaddat(naddatyp) = record(i:i+2)
	      goto 100
	    endif
	    do 99 inum = 1,5
	      write(inum1,'(i1)') inum
	      gtest = "G"//inum1
	      if(record(1:2).eq.gtest) then
	        i = 2
102	        i = i+1
	        if(record(i:i).eq.' ') goto 102
	        iline1 = i
103	        i = i+1
	        if(record(i:i).ne.' ') goto 103
	        if(inum.le.4) then
	          read(record(iline1:i-1),*) addasf(naddatyp,inum)
	        else
	          read(record(iline1:i-1),*) addcsf(naddatyp,1)
c	          write(6,*)  naddatyp,addcsf(naddatyp,1)
	          goto 106
	        endif
104	        i = i+1
	        if(record(i:i).eq.' ') goto 104
	        iline2 = i
105	        i = i+1
	        if(record(i:i).ne.' ') goto 105
	        read(record(iline2:i-1),*) addbsf(naddatyp,inum)
	      endif
99	    continue
106	    continue
	    if(record(1:2).eq."Xv") then
	      i = 2
107	      i = i+1
	      if(record(i:i).eq.' ') goto 107
	      iline1 = i
108	      i = i+1
	      if(record(i:i).ne.' ') goto 108
	      read(record(iline1:i-1),*) addXvol(naddatyp,1)
	    elseif(record(1:2).eq."Xr") then
	      i = 2
109	      i = i+1
	      if(record(i:i).eq.' ') goto 109
	      iline1 = i
110	      i = i+1
	      if(record(i:i).ne.' ') goto 110
	      read(record(iline1:i-1),*) addXrad(naddatyp,1)
	    elseif(record(1:2).eq."Hb") then
	      i = 2
111	      i = i+1
	      if(record(i:i).eq.' ') goto 111
	      iline1 = i
112	      i = i+1
	      if(record(i:i).ne.' ') goto 112
	      read(record(iline1:i-1),*) addHb(naddatyp,1)
	    endif
	    goto 100
200	    continue
	  close(unit=1)
	do 202 ilog = 1,2
	  if(ilog.eq.1) then
	    ulogtmp = 6
	  else
	    ulogtmp = ulog
	    open(unit=ulogtmp,file=flog,access='append',status='unknown')
	  endif
	    write(ulogtmp,2) naddatyp
	    do 201 i = 1,naddatyp
	      write(ulogtmp,4) nameaddat(i)
201	    continue
	  if(ilog.eq.2) then
	    close(unit=ulogtmp)
	  endif
202	continue
	endif
c	
	naddrtyp = 0
	if(faddrestype.ne.'') then
	  open(unit=1,file=faddrestype,status='unknown')
300	    read(1,1,end=400) record
	    if(record(1:1).eq."#") goto 300
	    if(record(1:3).eq."...") then
	      naddrtyp = naddrtyp + 1
	      goto 300
	    endif
	    read(record(1:3),*) nameaddres(naddrtyp)
	    read(record(5:10),*) iaddatom(naddrtyp)
	    do 301 i = 1,iaddatom(naddrtyp)
	      read(1,1) line
	      nameaddatom(naddrtyp,i) = line(1:4)
	      do 302 j = 1,14
	        if(line(5:7).eq.name2num(j)) then
	          atomaddtype(naddrtyp,i) = j
	        endif
302	      continue
	      do 303 j = 1,naddatyp
	        if(line(5:7).eq.nameaddat(j)) then
	          atomaddtype(naddrtyp,i) = 14+j
	        endif
303	      continue
301	    continue
	    goto 300
400	    continue
	  close(unit=1)
	  do 402 ilog = 1,2
	    if(ilog.eq.1) then
	      ulogtmp = 6
	    else
	      ulogtmp = ulog
	      open(unit=ulogtmp,file=flog,access='append',status='unknown')
	    endif
	    write(ulogtmp,3) naddrtyp
	    do 401 i = 1,naddrtyp
	      write(ulogtmp,5) nameaddres(i)
401	    continue
	    if(ilog.eq.2) then
	      close(unit=ulogtmp)
	    endif
402	  continue
	endif
c	
c ======= End of the subroutine =====================================
c	
	return
	end
