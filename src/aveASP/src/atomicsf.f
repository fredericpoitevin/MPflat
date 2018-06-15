c	
	subroutine atomicsf
c	
c	NOTE: all this might be tabulated... do it later...
c	
	include 'param.h'
c	
	integer i,j,k,j2
	integer natom,atomtyp(natot)
	integer isample,nsample
	integer loopstart
	integer ncount(14)
	integer addHb(naddamax,1),j3,j4,nHb
c	
	real*8 atom_crd(3*natot)
	real*8 asf(sftyp,4),bsf(sftyp,4),csf(sftyp,1),Xvol(sftyp,1)
        real*8 Xrad(sftyp,1)
	real*8 qrange,qstep,profil(nqmax),expq(nqmax)
        real*8 sigexp(nqmax),profexp(nqmax)
	real*8 qnorm,pi,pifour,fac,sum1
	real*8 fi(natot,nqmax),gi(natot,nqmax),fw(nqmax)
	real*8 wfi(natot,nqmax),wgi(natot,nqmax),wfw(nqmax)
	real*8 fityp(14,nqmax),gityp(14,nqmax)
	real*8 addasf(naddamax,4),addbsf(naddamax,4)
        real*8 addcsf(naddamax,1),addXvol(naddamax,1)
        real*8 addXrad(naddamax,1)
c	
c	character
c	
1	format(f12.7,14(1x,f12.7))
2	format(f12.7,1x,f12.7)
c	
	common /atomtyp/ natom,atomtyp,atom_crd
	common /sfparam/ asf,bsf,csf,Xvol,Xrad
	common /profinfo/ qrange,qstep,profil,expq,
     &                    profexp,sigexp,nsample
	common /atomsf/ fi,gi,fw
	common /addsfparam/ addasf,addbsf,addcsf,addXvol,addXrad,addHb
c	
	pi = acos(-1.d0)
	pifour = 4*pi
c	
	do 10 isample = 1,nsample
	  qnorm = expq(isample)
	  fac = qnorm*qnorm/pifour
c	-------------------------------------------------------------
c		Compute the atomic in vacuo form factors	
c		as well as the atomic excl.vol. one		
c	-------------------------------------------------------------
	  do 20 i = 1,natom
	      fi(i,isample) = 0.d0
	      gi(i,isample) = 0.d0
	      if(atomtyp(i).le.4) then
	        loopstart = 1
	      elseif(atomtyp(i).gt.4.and.atomtyp(i).le.8) then
	        loopstart = 5
	      elseif(atomtyp(i).gt.8.and.atomtyp(i).le.10) then
	        loopstart = 9
	      elseif(atomtyp(i).gt.10.and.atomtyp(i).le.12) then
	        loopstart = 11
	      else
	        loopstart = atomtyp(i)
	      endif
	      do 22 j = loopstart,atomtyp(i)
	        if(j.eq.loopstart) then
	          j2 = j
	        else
	          j2 = 14
	        endif
	        if(atomtyp(i).le.14) then
	         fi(i,isample) = fi(i,isample) + csf(j2,1)
	         do 21 k = 1,4
	          sum1 = -fac*bsf(j2,k)/pifour
	          sum1 = asf(j2,k)*exp(sum1)
	          fi(i,isample) = fi(i,isample) + sum1
21	         continue
	        else
	         fi(i,isample) = fi(i,isample) + addcsf(j2-14,1)
	         do 23 k = 1,4
	          sum1 = -fac*addbsf(j2-14,k)/pifour
                  sum1 = addasf(j2-14,k)*exp(sum1)
                  fi(i,isample) = fi(i,isample) + sum1
23	         continue
	         j3 = atomtyp(i)
	         nHb = addHb(j3-14,1)
	         if(nHb.ne.0) then
	           do 24 j4 = 1,nHb
	             fi(i,isample) = fi(i,isample) + csf(14,1)
	             do 25 k = 1,4
	               sum1 = -fac*bsf(14,k)/pifour
	               sum1 = asf(14,k)*exp(sum1)
	               fi(i,isample) = fi(i,isample) + sum1
25	             continue
24	           continue
	         endif
	        endif
22	      continue
	      j = atomtyp(i)
	      if(j.le.14) then
	       sum1 = -fac*(Xvol(j,1))**(2/3)
	       gi(i,isample) = gi(i,isample) + Xvol(j,1)*exp(sum1)
	      else
               sum1 = -fac*(addXvol(j-14,1))**(2/3)
               gi(i,isample) = gi(i,isample) + addXvol(j-14,1)*exp(sum1)
	      endif
20	  continue
c	-------------------------------------------------------------
c		Compute the water form factor			
c	-------------------------------------------------------------
	fw(isample) = 0.d0
	do 30 i = 1,3
	  if(i.eq.1) then
	    j = 9
	  else
	    j = 14
	  endif
	  fw(isample) = fw(isample) + csf(j,1)
	  do 31 k = 1,4
	    sum1 = -fac*bsf(j,k)/pifour
	    sum1 = asf(j,k)*exp(sum1)
	    fw(isample) = fw(isample) + sum1
31	  continue
30	continue
10	continue
c	
c	-------------------------------------------------------------
c		Check: Write for each atom type			
c	-------------------------------------------------------------
	if(1.eq.0) then
	do 400 j = 1,14
	  ncount(j) = 0
	  do 401 i = 1,natom
	    if(atomtyp(i).eq.j) then
	      ncount(j) = 1
	      do 402 isample = 1,nsample
	        fityp(j,isample) = fi(i,isample)
	        gityp(j,isample) = gi(i,isample)
402	      continue
	      goto 400
	    endif
401	  continue
400	continue
	do 403 j = 1,14
	  if(ncount(j).eq.0) then
	    do 404 isample = 1,nsample
	      fityp(j,isample) = 0.d0
	      gityp(j,isample) = 0.d0
404	    continue
	  endif
403	continue
	open(unit=2,file="fi.dat",status='unknown')
	open(unit=3,file="gi.dat",status='unknown')
	open(unit=4,file="fw.dat",status='unknown')
	  do 500 isample = 1,nsample
	    write(2,1) expq(isample),fityp(1,isample),
     &	               fityp(2,isample),
     1	               fityp(3,isample),fityp(4,isample),
     2	               fityp(5,isample),fityp(6,isample),
     3	               fityp(7,isample),fityp(8,isample),
     4	               fityp(9,isample),fityp(10,isample),
     5	               fityp(11,isample),fityp(12,isample),
     6	               fityp(13,isample),fityp(14,isample)
	    write(3,1) expq(isample),gityp(1,isample),
     &	               gityp(2,isample),
     1	               gityp(3,isample),gityp(4,isample),
     2	               gityp(5,isample),gityp(6,isample),
     3	               gityp(7,isample),gityp(8,isample),
     4	               gityp(9,isample),gityp(10,isample),
     5	               gityp(11,isample),gityp(12,isample),
     6	               gityp(13,isample),gityp(14,isample)
	    write(4,2) expq(isample),fw(isample)
500	  continue
	close(unit=4)
	close(unit=3)
	close(unit=2)
	endif
c	-------------------------------------------------------------
c		End of the subroutine				
c	-------------------------------------------------------------
	return
	end
