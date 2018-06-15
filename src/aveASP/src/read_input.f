c	
	subroutine readinput(fparam)
c	
	include 'param.h'
c	
	integer nkeys
	integer i,j,imin,imax
	integer flagsolv,flagfit,flagmeth,fsolvfmt
	integer flagunit
	integer ulog,nsample
	integer naddrtyp,naddatyp
	integer restrictchain,nch
	integer flagcuser
c	
	real*8 bulkdens,qrange,qstep,profil(nqmax)
	real*8 expq(nqmax),sigexp(nqmax),profexp(nqmax)
	real*8 c1user,c2user
c	
	character keys(13)*35,fparam*80,jnum*1
	character*64 fpdb,fasexp,fsolvent,fprof,flog
	character*64 faddrestype,faddatomtype
	character line*80,param*45,keyword*35,value*45
	character okchain(20)*1
c	
1	format(a)
2	format("#### Welcome to AquaSAXS ####")
3	format("# Input/Parameter Summary #")
4	format(">> Name of input PDB/PQR file: ",a30)
5	format(">> Fitting will be performed...")
6	format(">> Name of exp. SAXS profile : ",a30)
7	format(">> Solvation accounted for through SAS")
8	format(">> Solvation accounted for through solvation map")
9	format(">>     Name of input Solvent map: ",a30)
c	
	logical file_exists
c	
	common /files/		fpdb,fasexp,fsolvent,fprof,flog,ulog
	common /flags/		flagsolv,flagfit,flagmeth,fsolvfmt
	common /params/		bulkdens
	common /profinfo/	qrange,qstep,profil,expq,
     &	                 	profexp,sigexp,nsample
	common /qunit/		flagunit
	common /addfile/	faddrestype,faddatomtype
	common /chains/		restrictchain,nch,okchain
	common /usercs/		flagcuser,c1user,c2user
c	
	data nkeys /13/
	data (keys(i),i=1,13) /
     &	'SOLUTE  pdb                       :',
     &	'EXPER.  saxs profile filename     :',
     &	'SOLVATION SAS(0) or Dipolar(1)  ? :',
     &	'  if Dipolar, solvent map name    :',
     &	'BULK    Density (default:0.334)   :',
     &	'PROFILE q-max                     :',
     &	'PROFILE q-step                    :',
     &	'PROFILE unit                      :',
     &	'ADDITIONAL res/monomer type file  :',
     &	'ADDITIONAL atom types             :',
     &	'Chain(s) to read ?                :',
     &	'If no fit and no default, c1 value:',
     &	'If no fit and no default, c2 value:'
     &	/
c	
c ======= Defaults ==================================================
c	
	bulkdens = 0.334
	flagsolv = 0
	fasexp = ''
	flagfit = 0
	faddrestype = ''
	faddatomtype = ''
	flagunit = 1
	restrictchain = 0
	flagcuser = 0
c	
c ======= Read input file ===========================================
c	
	open(unit=1,file=fparam,status='old')
10	  read(1,1,end=20) line
	    param = '                              '
	    if(line(1:1).eq.'#') goto 10
	    keyword = line(1:35)
	    value   = line(36:80)
c	
	    do 30 i = 1,45
	      if(value(i:i).eq.'!') then
	        imax = i-1
	        goto 40
	      endif
30	    continue
	    imax = 45
40	    continue
	    do 50 i = imax,1,-1
	      if(value(i:i).ne.' ') goto 60
50	    continue
60	    continue
	    imax = i
	    do 70 i = 1,imax
	      if(value(i:i).ne.' ') goto 80
70	    continue
80	    continue
	    imin = i
	    param(1:imax-imin+1) = value(imin:imax)
	    if(param.eq.'') goto 10
c	
	    do 100 i = 1,nkeys
	      if(keyword.eq.keys(i)) goto 110
100	    continue
110	    continue
c	
	    if(i.eq.1) then
	      read(param,*) fpdb
	    elseif(i.eq.2) then
	      read(param,*) fasexp
	    elseif(i.eq.3) then
	      read(param,*) flagsolv
	    elseif(i.eq.4) then
	      read(param,*) fsolvent
	    elseif(i.eq.5) then
	      read(param,*) bulkdens
	    elseif(i.eq.6) then
	      read(param,*) qrange
	    elseif(i.eq.7) then
	      read(param,*) qstep
	    elseif(i.eq.8) then
	      read(param,*) flagunit
	    elseif(i.eq.9) then
	      read(param,*) faddrestype
	    elseif(i.eq.10) then
	      read(param,*) faddatomtype
	    elseif(i.eq.11) then
	      restrictchain = 1
	      do 200 j = 1,imax-imin+1
	        okchain(j) = param(j:j)
200	      continue
	      nch = imax-imin+1
	    elseif(i.eq.12) then
	      read(param,*) c1user
	      flagcuser = 1
	    elseif(i.eq.13) then
	      read(param,*) c2user
	      flagcuser = 1
	    endif
c	
	  goto 10
20	  continue
	close(unit=1)
c	
c ======= Check map format, if needed ===============================
c	
        if(flagsolv.eq.1) then
	  fsolvfmt = 0
          i = 0
400       i = i+1
          if(fsolvent(i:i+2).eq.'.cn') then
            fsolvfmt = 1
            goto 410
          else
            goto 400
          endif
          if(fsolvfmt.eq.0) stop
	endif
410     continue
c	
c ======= Define output names =======================================
c	
	i = 0
500	i = i+1
	if(fpdb(i:i+3).eq.'.pdb'.or.fpdb(i:i+3).eq.'.pqr') then
	  j = -1
501	  j = j+1
	  write(jnum,'(i1)') j
	  fprof = fpdb(1:i-1)//jnum//".saxs"
	  flog  = fpdb(1:i-1)//jnum//".log"
	  inquire(file=flog,exist=file_exists)
	  if(file_exists) then
	    goto 501
	  else
	    goto 510
	  endif
	else
	  goto 500
	endif
510	continue
c	
	do 511 i = 1,2
	  if(i.eq.1) then
	    ulog = 6
	  else
	    ulog = 15
	    open(unit=ulog,file=flog,status='unknown')
	  endif
	  write(ulog,2)
	  write(ulog,*)" "
	  write(ulog,3)
c	  write(ulog," ")
	  write(ulog,4) fpdb
	  if(fasexp.ne.'') then
	    flagfit = 1
	    write(ulog,5)
	    write(ulog,6) fasexp
	  endif
	  if(flagsolv.eq.0) then
	    write(ulog,7)
	  else
	    write(ulog,8)
	    write(ulog,9) fsolvent
	  endif
	  if(i.eq.2) then
	    close(unit=ulog)
	  endif
511	continue
c	
c ======= End of the subroutine =====================================
c	
	return
	end
c____________________________________________________________________
c	
c	This subroutine only retains heavy atoms,
c	and attributes atom type...
c	
	subroutine readpdb
c	
	include 'param.h'
c	
	integer i,j,k,jmin,jmax
	integer ulog,ilog,ulogtmp,natom,natom2,resnum(natot)
	integer iatom(nresdef),atomtype(nresdef,24)
	integer ntype,foundres,found(natot),atomtyp(natot)
	integer naddatyp,naddrtyp,iaddatom(naddrmax)
	integer atomaddtype(naddrmax,100)
	integer errexistr,errexista
	integer restrictchain,nch
	integer testch,ich,readok,nwat
	integer flageisenb,atmsolvtyp(natot)
	integer addHb(naddamax,1)
c	
	real*8 atom_crd(3*natot)
	real*8 asf(sftyp,4),bsf(sftyp,4),csf(sftyp,1),Xvol(sftyp,1)
        real*8 Xrad(sftyp,1),averad
	real*8 addasf(naddamax,4),addbsf(naddamax,4)
        real*8 addcsf(naddamax,1),addXvol(naddamax,1)
        real*8 addXrad(naddamax,1)
	real*8 occ,beta
c	
	character*64 fpdb,fasexp,fsolvent,fprof,flog,outpdb
	character atom*6,record*80,chain(natot)*1
	character atname(natot)*4,resname(natot)*3
	character nameres(nresdef)*3,nameatom(nresdef,24)*3
	character nameaddres(naddrmax)*3,nameaddatom(naddrmax,100)*4
	character okchain(20)*1,hetatm*6,wat*3
	character col3*3,col4*4,col9*4,col13*3,test*4
	character flib*60,restest*3
	character test2*1
c	
	logical file_exists
c	
1	format(a)
2	format(12x,a4,1x,a3,1x,a1,i4,4x,3(f8.3))
3	format("ERR > Did not recognize residue: ",
     1	i5,1x,a3,1x,a3,1x,a1,1x,i4)
4	format("ERR >    Did not recognize atom: ",
     1  i5,1x,a4,1x,a3,1x,a1,1x,i4)
5	format(">> Number of atoms read: ",i6)
6	format("# Reading PDB/PQR file #")
7	format(">> Number of atoms recognized: ",i6)
8	format("ERR > Consider defining a new residue/atom type?")
9	format(">> Average atomic radius for this molecule: ",f5.3," A")
10	format(">> Number of water molecule discarded: ",i5)
11	format(a6,i6,a4,1x,a3,1x,a1,i4,4x,3(f8.3),2(f6.2))
12	format(">> Chain read: ",a1)
c	
	common /files/ fpdb,fasexp,fsolvent,fprof,flog,ulog
	common /atominfo/ natom2,atname,resname,chain,resnum
	common /atomtyp/ natom,atomtyp,atom_crd
	common /sfparam/ asf,bsf,csf,Xvol,Xrad
	common /radave/ averad
	common /addtype/ naddrtyp,naddatyp,nameaddres,iaddatom,
     1	nameaddatom,atomaddtype
	common /addsfparam/ addasf,addbsf,addcsf,addXvol,addXrad,addHb
	common /chains/ restrictchain,nch,okchain
        common /atomsolventtype/ atmsolvtyp
c	
c====================================================================
c	Parameters						
c====================================================================
	atom = "ATOM  "
	hetatm = "HETATM"
	wat = "HOH"
	natom = 0
	averad = 0.d0
	errexistr = 0
	errexista = 0
	outpdb = "filtered.pdb"
	occ = 1.d0
	beta = 2.d1
	nwat = 0
c	
c====================================================================
c	Read/Edit PDB/PQR file					
c====================================================================
	open(unit=1,file=fpdb,status='old')
100	  read(1,1,end=200) record
c	
c	    ### Only read lines starting with ATOM ot HETATM ###
	    readok = 0
	    if(record(1:6).eq.atom) then
	      readok = 1
	    elseif(record(1:6).eq.hetatm) then
	      readok = 1
	    endif
	    if(readok.eq.0) goto 100
c	
c	    ### Discard hydrogens ###
	    if(record(13:13).eq.'H'.or.record(14:14).eq.'H') goto 100
c	
c	    ### Discard water ###
	    if(record(18:20).eq.wat) then
	      nwat = nwat+1
	      goto 100
	    endif
c	
	    natom = natom + 1
	    read(record,2) atname(natom),resname(natom),
     1	        chain(natom),resnum(natom),atom_crd(3*(natom-1)+1),
     2	            atom_crd(3*(natom-1)+2),atom_crd(3*(natom-1)+3)
	    if(restrictchain.eq.1) then
	      testch = 0
	      do 120 ich = 1,nch
	        if(chain(natom).eq.okchain(ich)) then
	          testch = 1
	        endif
120	      continue
	      if(testch.eq.0) then
	        natom = natom - 1
	      endif
	    endif
	  goto 100
200	  continue
	close(unit=1)
c	
	do 202 ilog = 1,2
	  if(ilog.eq.1) then
	    ulogtmp = 6
	  else
	    ulogtmp = ulog
	open(unit=ulogtmp,file=flog,access='append',status='unknown')
	  endif
c	  write(6,5) natom
	  write(ulogtmp,*) ' '
	  write(ulogtmp,*) ' '
	  write(ulogtmp,6)
	  write(ulogtmp,*) ' '
	  write(ulogtmp,5) natom
	  write(ulogtmp,10) nwat
	  do 201 ich = 1,nch
	    write(ulogtmp,12) okchain(ich)
201	  continue
	  if(ilog.eq.2) then
	close(unit=ulogtmp)
	  endif
202	continue
c	
c====================================================================
c	Get atomic types					
c====================================================================
c	
c       First, the solvation type
c       N (neutral)     1
c       O (neutral)     2
c       N (charged)     3
c       O (charged)     4
c       C               5
c       S               6
c
        do 250 i = 1,natom
          test=atname(i)
          restest=resname(i)
          if(test(2:2).eq.'N') then
            if(restest.eq.'ARG'.and.test(2:4).eq.'NH1') then
              atmsolvtyp(i) = 3
            elseif(restest.eq.'LYS'.and.test(2:4).eq.'NZ') then
              atmsolvtyp(i) = 3
            else
              atmsolvtyp(i) = 1
            endif
          elseif(test(2:2).eq.'O') then
            if(restest.eq.'ASP'.and.test(2:4).eq.'OD1') then
              atmsolvtyp(i) = 4
            elseif(restest.eq.'GLU'.and.test(2:4).eq.'OE1') then
              atmsolvtyp(i) = 4
            else
              atmsolvtyp(i) = 2
            endif
          elseif(test(2:2).eq.'C') then
            atmsolvtyp(i) = 5
          elseif(test(2:2).eq.'S') then
            atmsolvtyp(i) = 6
          else
            atmsolvtyp(i) = 0
          endif
250     continue       
c
c       Second, the scattering type
c       
	natom2 = 0
	do 300 i = 1,natom
	  foundres = 0
	  found(i) = 0
	  restest = resname(i)
	  j=0
801	  j=j+1
	  if(restest(j:j).eq." ") goto 801
	  jmin=j
	  j=4
803	  j=j-1
	  if(restest(j:j).eq." ") goto 803
	  jmax=j
c	  flib = "../../aquasaxs/fraglib/"//restest(jmin:jmax)//".lib"
	  flib = "./fraglib/"//restest(jmin:jmax)//".lib"
	  inquire(file=flib,exist=file_exists)
          if(file_exists) then
	    foundres = 1
c	
c	  ### Scan the standard equivalence table ###
	   open(unit=2,file=flib,status='unknown')
110	    read(2,1,end=210) record
	    col3 = record(1:3)
	    col4 = record(1:4)
	    col9 = record(6:9)
	    col13= record(11:13)
	    if(record(13:13).eq." ") goto 110
	    if(foundres.eq.1) then
c	      ### careful with 3 or 4-letter nomenclature ###
	      test = atname(i)
	      if(test(1:1).eq.' ') then
	        col4(1:4) = " "//col4(1:3)
	      endif
	      if(col4.eq.atname(i)) then
	        read(col13,*) atomtyp(i)
	        found(i) = 1
	        goto 210
	      endif
	      if(found(i).eq.0.and.col9.eq.atname(i)) then
	        read(col13,*) atomtyp(i)
	        found(i) = 1
	        goto 210
	      endif
	    endif
	    goto 110
210	    continue
	   close(unit=2)
	  endif
c	
c	  ### Scan the user-defined equivalence table ###
c	  ### Note that it over-rides the precedent   ###
	  do 303 j = 1,naddrtyp
	    if(resname(i).eq.nameaddres(j)) then
	      foundres = 1
c	      ### careful with 3 or 4-letter nomenclature ###
	      test = atname(i)
	      if(test(1:1).eq.' ') then
	        test(1:4)=test(2:4)//' '
	      endif
	      do 304 k = 1,iaddatom(j)
	        if(test.eq.nameaddatom(j,k)) then
c	        if(atname(i).eq.nameaddatom(j,k)) then
	          found(i) = 1
	          atomtyp(i) = atomaddtype(j,k)
	        endif
304	      continue
	    endif
303	  continue
	  if(found(i).eq.0) then
	    errexista = 1
	    do 305 ilog = 1,2
	      if(ilog.eq.1) then
	        ulogtmp = 6
	      else
	        ulogtmp = ulog
	    open(unit=ulogtmp,file=flog,access='append',status='unknown')
	      endif
             write(ulogtmp,4) i,atname(i),resname(i),chain(i),resnum(i)
	      if(ilog.eq.2) then
            close(unit=ulogtmp)
	      endif
305	    continue
	  else
c	
c	    ### Reorder the arrays without gaps ###
	    natom2 = natom2 + 1
	    atname(natom2) = atname(i)
	    resname(natom2) = resname(i)
	    chain(natom2) = chain(i)
	    resnum(natom2) = resnum(i)
	    atomtyp(natom2) = atomtyp(i)
	    if(atomtyp(natom2).le.14) then
	      averad = averad + Xrad(atomtyp(i),1)
	    else
	      averad = averad + addXrad(atomtyp(i)-14,1)
	    endif
            atmsolvtyp(natom2) = atmsolvtyp(i)
c	
	    do 130 k = 1,3
	      atom_crd(3*(natom2-1)+k) = atom_crd(3*(i-1)+k)
130	    continue
	    open(unit=2,file=outpdb,access='append',status='unknown')
	      write(2,11) atom,i,atname(i),resname(i),chain(i),
     1	resnum(i),atom_crd(3*(i-1)+1),atom_crd(3*(i-1)+2),
     1	atom_crd(3*(i-1)+3),occ,beta
	    close(unit=2)
	  endif
300	continue
c	
	averad = averad/real(natom2)
	do 306 ilog = 1,2
	  if(ilog.eq.1) then
	    ulogtmp = 6
	  else
	    ulogtmp = ulog
	open(unit=ulogtmp,file=flog,access='append',status='unknown')
	  endif
	  if(errexista.eq.1) then
	    write(ulogtmp,8)
	  endif
	  write(ulogtmp,9) averad
	  write(ulogtmp,7) natom2
	  natom = natom2
	  if(ilog.eq.2) then
	close(unit=ulog)
	  endif
306	continue
c	
c ======= End of the subroutine =====================================
c	
	return
	end
