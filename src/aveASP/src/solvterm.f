c	
c	subroutine solvterm
c	
c	Computes the hydration shell contribution
c	
	subroutine solvterm
c	
	include 'param.h'
c	
	integer i,j,k,j2
	integer natom2
	integer flagsolv,flagfit,flagmeth,fsolvfmt
	integer natom,atomtyp(natot)
	integer nneigh(natot),neighbour(natot,1000)
	integer ncub,icub,ncount
	integer resnum(natot),npoint(natot)
	integer ibeg(3),iend(3),ipas(3)
	integer ulog,ilog,ulogtmp
	integer addHb(naddamax,1)
c	
	real*8 pi
	real*8 asf(sftyp,4),bsf(sftyp,4),csf(sftyp,1),Xvol(sftyp,1)
        real*8 Xrad(sftyp,1)
        real*8 addasf(naddamax,4),addbsf(naddamax,4)
        real*8 addcsf(naddamax,1),addXvol(naddamax,1)
        real*8 addXrad(naddamax,1)
	real*8 atom_crd(3*natot),rad(natot),maxrad
	real*8 Rprobe,cutoff,dij
	real*8 qpos(3,ncubmax),qweight(1,ncubmax)
	real*8 pos(3),sas(natot),sastot
	real*8 map(ix,iy,iz),cell(6),sigma,ave
	real*8 resol
c	
	character*3 atname(natot),resname(natot),chain(natot)*1
	character*64 fpdb,fasexp,fsolvent,fprof,flog
c	
	logical qform
c	
1	format(i5,1x,a3,1x,a3,1x,a1,1x,i4,1x,i5,1x,i5,1x,f12.7)
2	format("ATOM    188  CG2 VAL    13    ",3f8.3)
3	format("# Computing/Storing info on the hydration term #")
4	format(">> Approx. total Solvent Accessible Surface: ",f10.3," A2")
5	format(">> Solvent map resolution: ",f10.3," A")
c	
	common /flags/ flagsolv,flagfit,flagmeth,fsolvfmt
	common /files/ fpdb,fasexp,fsolvent,fprof,flog,ulog
	common /atominfo/ natom2,atname,resname,chain,resnum
	common /atomtyp/ natom,atomtyp,atom_crd
	common /cubparam/ qpos,qweight,ncub
	common /sfparam/ asf,bsf,csf,Xvol,Xrad
	common /asas/ sas
	common /addsfparam/ addasf,addbsf,addcsf,addXvol,addXrad,addHb
	common /smap/ ibeg,iend,ipas,cell,map
c	
c	-------------------------------------------------------------
c	Parameters						
c	-------------------------------------------------------------
	do 100 ilog = 1,2
	  if(ilog.eq.1) then
	    ulogtmp = 6
	  else
	    ulogtmp = ulog
	open(unit=ulogtmp,file=flog,access='append',status='unknown')
	  endif
	  write(ulogtmp,*)' '
	  write(ulogtmp,*)' '
	  write(ulogtmp,3)
	  write(ulogtmp,*)' '
100	continue
c	
	maxrad = 0.d0
	Rprobe = 1.4
	sastot = 0.d0
	pi = acos(-1.d0)
	qform = .true.
	if(flagsolv.eq.0) then
c	-------------------------------------------------------------
c	SAS approach						
c	-------------------------------------------------------------
c	=== Find maxrad ===
	  do 10 i = 1,natom
	    if(atomtyp(i).le.14) then
	      rad(i) = Xrad(atomtyp(i),1)
	    else
	      rad(i) = addXrad(atomtyp(i)-14,1)
	    endif
	    if(rad(i).gt.maxrad) maxrad = rad(i)
10	  continue
c	=== List neighbours ===
	  do 20 i = 1,natom
	    nneigh(i) = 0
	    do 21 j = 1,natom
	      if(i.eq.j) goto 21
	      cutoff = 4*(Rprobe+maxrad)*(Rprobe+maxrad) 
	      dij = 0.d0
	      do 22 k = 1,3
	        dij = dij + 
     1	              (atom_crd(3*(j-1)+k)-atom_crd(3*(i-1)+k))**2
22	      continue
c	      dij = sqrt(dij)
	      if(dij.le.cutoff) then
	        nneigh(i) = nneigh(i) + 1
	        neighbour(i,nneigh(i)) = j
	      endif
21	    continue
20	  continue
c	=== Compute atomic SA surface ===
	  do 30 i = 1,natom
	    ncount = ncub
	    do 31 icub = 1,ncub
	      do 32 k = 1,3
	        pos(k) = atom_crd(3*(i-1)+k) +
     1	            (Rprobe+rad(i))*qpos(k,icub)
32	      continue
	      do 33 j = 1,nneigh(i)
	        j2 = neighbour(i,j)
	        cutoff = (rad(j2) + Rprobe)**2
	        dij = 0.d0
	        do 34 k = 1,3
	          dij = dij +
     1	           (atom_crd(3*(j2-1)+k)-pos(k))**2
34	        continue
	        if(dij.lt.cutoff) then
	          ncount = ncount - 1
	          goto 31
	        endif
33	      continue
31	    continue
	    sas(i) = 4*pi*(rad(i)+Rprobe)*(rad(i)+Rprobe)
	    sas(i) = sas(i)*real(ncount)/real(ncub)
	    sastot = sastot + sas(i)
	    npoint(i) = ncount
30	  continue
c	... Check: write output ...
	do 101 ilog = 1,2
	  if(ilog.eq.1) then
	    ulogtmp = 6
	  else
	    ulogtmp = ulog
	  endif
	  write(ulogtmp,4) sastot
101	continue
	else
c	-------------------------------------------------------------
c	Map approach						
c	-------------------------------------------------------------
c	==== Read the map (only CNS allowed for now...) ===
	  call readmap_cns(ibeg,iend,ipas,map,ave,sigma,cell,
     1	                   fsolvent,qform)
	  resol = cell(1)/ipas(1)
	  do 102 ilog = 1,2
	    if(ilog.eq.1) then
	      ulogtmp = 6
	    else
	      ulogtmp = ulog
	    endif
	    write(ulogtmp,5) resol
102	  continue
	  call DefGrid
	endif
c	
c	-------------------------------------------------------------
c	End of the subroutine					
c	-------------------------------------------------------------
	close(unit=ulog)
999	continue
	return
	end
