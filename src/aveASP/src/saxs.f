c	
c	subroutine saxs.f
c	
	subroutine saxs
c	
	include 'param.h'
c	
	integer i,j,k,iatom,ulog,ilog,ulogtmp
	integer isample,nsample,natom,atomtyp(natot)
	integer flagsolv,flagfit,flagmeth,fsolvfmt
	integer iexcl,ishell
	integer ma,iafit(1)
	integer ncub
	integer ibeg(3),iend(3),ipas(3)
	integer nx,ny,nz,npoint,nsize_byte
	integer ia,ib,ic
	integer flagcuser
c	
	real*8 qrange,qstep,profil(nqmax),expq(nqmax)
        real*8 sigexp(nqmax),profexp(nqmax)
	real*8 fi(natot,nqmax),gi(natot,nqmax),fw(nqmax)
	real*8 sas(natot)
	real*8 bulkdens,c1,c2,c1best,c2best
	real*8 chisq,chisqmin,covar(1,1)
	real*8 qnorm,intens,faci,facj,dij
	real*8 atom_crd(3*natot)
	real*8 proftmp(nqmax)
	real*8 afit(1),cbest
	real*8 qpos(3,ncubmax),qweight(1,ncubmax)
	real*8 wj,a,b,dot,ri(3),qj(3),cqjri,sqjri
	real*8 dotu(nqmax,natot)
	real*8 vc(nqmax,ncubmax),vs(nqmax,ncubmax)
	real*8 xc(nqmax,ncubmax),xs(nqmax,ncubmax)
	real*8 wc(nqmax,ncubmax),ws(nqmax,ncubmax)
	real*8 u(nqmax,1),v(1,1),w(1)
	real*8 map(ix,iy,iz),cell(6),gridvol
	real*8 hx(ix),hy(iy),hz(iz)
	real*8 xsmap(*),posx(*),posy(*),posz(*)
	real*8 cutmin,cutmax,tol
	real*8 averad,rzero,pi,rzerobest,facexp
	real*8 c1tmp
	real*8 c1user,c2user
	real*8 wsum
c	
	character charstart*8,charend*8
	character*64 fpdb,fasexp,fsolvent,fprof,flog
c	
	pointer (ptr_xsmap,xsmap)
        pointer (ptr_posx,posx)
        pointer (ptr_posy,posy)
        pointer (ptr_posz,posz)
c	
1	format(a1,$)
2	format(" >>> Percent done: ",$)
3	format(f5.0,$)
4	format(" New chisqmin: ",e12.7)
5	format("# Actual computation of the SAXS profile #")
6	format(">> Profile computation started at: ",a8)
7	format(">> Profile computation ended at: ",a8)
8	format(">> Compressing solvent excess density map...")
9	format(">>   initial size:    ",i12)
800	format(">>   compressed size: ",i12)
801	format(">> C1 scanned range: ",f5.3," to ",f5.3)
802	format(">> C2 scanned range: ",f5.3," to ",f5.3)
804	format(">> Integration volume (solvent region): ",f6.3," A3")
805	format(">> C1 used value: ",f5.3)
806	format(">> C2 used value: ",f5.3)
807	format(">>    ie an average atomic radius of: ",f5.3," A")
808	format(3(f12.5))
809	format("# Goodness-of-fit as a function of C1 and C2 #")
810	format("#     C1          C2          Chi    #")
811	format(">> Chi value    : ",e10.5)
c	
	common /files/ fpdb,fasexp,fsolvent,fprof,flog,ulog
	common /flags/ flagsolv,flagfit,flagmeth,fsolvfmt
	common /atomtyp/ natom,atomtyp,atom_crd
	common /profinfo/ qrange,qstep,profil,expq,
     &                    profexp,sigexp,nsample
        common /atomsf/ fi,gi,fw
	common /asas/ sas
	common /cubparam/ qpos,qweight,ncub
	common /fitparam/ c1best,c2best,cbest,chisqmin
	common /smap/ ibeg,iend,ipas,cell,map
	common /posmap/ hx,hy,hz,nx,ny,nz
	common /radave/ averad
	common /add/ c1tmp
	common /usercs/ flagcuser,c1user,c2user
	common /params/ bulkdens
c	
	chisqmin = 1.d10
	tol = 1.d-4
	cutmin = 1.d0-tol
	cutmax = 1.d0+tol
	pi = acos(-1.d0)
	wsum = 12.5663706143165
c	
	do 706 ilog = 1,2
	  if(ilog.eq.1) then
	    ulogtmp = 6
	  else
	    ulogtmp = ulog
	open(unit=ulogtmp,file=flog,access='append',status='unknown')
	  endif
	  write(ulogtmp,*) ' '
	  write(ulogtmp,*) ' '
	  write(ulogtmp,5)
	  write(ulogtmp,*) ' '
706	continue
c	
	call time(charstart)
c	=============================================================
c		Spherical averaging using the cubature formula	
c	=============================================================
c	
c	 ### If dipolar option allocate compressed memory ###
	  if(flagsolv.eq.1) then
	    gridvol = (cell(1)/ipas(1))**3
	    nsize_byte = 8*nx*ny*nz
	    ptr_xsmap  = malloc(nsize_byte)
	    ptr_posx   = malloc(nsize_byte)
	    ptr_posy   = malloc(nsize_byte)
	    ptr_posz   = malloc(nsize_byte)
	    do 707 ilog = 1,2
	      if(ilog.eq.1) then
	        ulogtmp = 6
	      else
	        ulogtmp = ulog
	      endif
	      write(ulogtmp,804) gridvol
              write(ulogtmp,8)
707	    continue
            npoint = 0
            do 500 ic = 1,nz
              do 499 ib = 1,ny
                do 498 ia = 1,nx
	          if(map(ia,ib,ic).eq.0) goto 498
                  if(map(ia,ib,ic).gt.cutmin.and.
     1	            map(ia,ib,ic).lt.cutmax) goto 498
                    npoint = npoint + 1
                    xsmap(npoint) = map(ia,ib,ic) - 1.d0
                    posx(npoint) = hx(ia)
                    posy(npoint) = hy(ib)
                    posz(npoint) = hz(ic)
498             continue
499           continue
500         continue
	    do 708 ilog = 1,2
	      if(ilog.eq.1) then
	        ulogtmp = 6
	      else
	        ulogtmp = ulog
	      endif
              write(ulogtmp,9) nx*ny*nz
              write(ulogtmp,800) npoint
708	    continue
	  else
	    nsize_byte = 0
	    ptr_xsmap = malloc(nsize_byte)
	    ptr_posx   = malloc(nsize_byte)
            ptr_posy   = malloc(nsize_byte)
            ptr_posz   = malloc(nsize_byte)
	  endif
c	
c	  ### Compute individual terms first ###
	  do 199 j = 1,ncub
	    do 198 iatom = 1,natom
	      dotu(j,iatom) = 0.d0
	      do 197 k = 1,3
	        dotu(j,iatom) = dotu(j,iatom)
     1         + atom_crd(3*(iatom-1)+k)*qpos(k,j)
197	      continue
198	    continue
199	  continue
	  do 200 isample = 1,nsample
	    do 201 j = 1,ncub
	      vc(isample,j) = 0.d0
	      vs(isample,j) = 0.d0
	      xc(isample,j) = 0.d0
	      xs(isample,j) = 0.d0
	      wc(isample,j) = 0.d0
	      ws(isample,j) = 0.d0
	      do 202 iatom = 1,natom
	        dot = expq(isample)*dotu(j,iatom)
	        cqjri = cos(dot)
	        sqjri = sin(dot)
	        vc(isample,j) = vc(isample,j)
     1	          + fi(iatom,isample)*cqjri
	        vs(isample,j) = vs(isample,j)
     1	          + fi(iatom,isample)*sqjri
	        xc(isample,j) = xc(isample,j)
     1	          + gi(iatom,isample)*cqjri
	        xs(isample,j) = xs(isample,j)
     1	          + gi(iatom,isample)*sqjri
202	      continue
c	
	      if(flagsolv.eq.0) then
	        do 203 iatom = 1,natom
	          dot = expq(isample)*dotu(j,iatom)
                  cqjri = cos(dot)
                  sqjri = sin(dot)
	          wc(isample,j) = wc(isample,j)
     1	            + sas(iatom)*fw(isample)*cqjri
	          ws(isample,j) = ws(isample,j)
     1	            + sas(iatom)*fw(isample)*sqjri
203	        continue
	      else
	        do 204 ic = 1,npoint
	          ri(1) = posx(ic)
	          ri(2) = posy(ic)
	          ri(3) = posz(ic)
	          dot = 0.d0
	          do 205 k = 1,3
	            dot = dot + ri(k)*qpos(k,j)*expq(isample)
205	          continue
	          cqjri = cos(dot)
	          sqjri = sin(dot)
c	
	          wc(isample,j) = wc(isample,j)
     1	            + xsmap(ic)*cqjri*gridvol
	          ws(isample,j) = ws(isample,j)
     1	            - xsmap(ic)*sqjri*gridvol
204	        continue
	      endif
201	    continue
200	  continue
c	
	  if(flagfit.eq.1) then
c	    ### Starting fitting mode ###
c	    ### Warn the user about the fitting range ###
	    do 709 ilog = 1,2
	      if(ilog.eq.1) then
	        ulogtmp = 6
	      else
	        ulogtmp = ulog
	      endif
	      write(ulogtmp,801) 9.d-1,1.12d0
	      if(flagsolv.eq.0) then
	        write(ulogtmp,802) 0.d0,4.d0
	      else
	        write(ulogtmp,802) 0.d0,1.4d0
	      endif
709	    continue
	    open(unit=1,file="isochi.txt",access='append',
     1               status='unknown')
	      write(1,809)
	      write(1,810)
	    close(unit=1)
c	
c	    ### Now mix terms. You warned them ###
	    do 300 iexcl = 1,101
c	
c	      ### Scan the SAV weight ###
	      c1 = 9.d-1 + (real(iexcl)-1)*2.d-3
	      c1tmp = c1
c	
c	      ### SEV term is scaled with a gaussian ###
	      rzero = c1*averad
	      c1 = bulkdens*c1*c1*c1
c	
	      do 301 ishell = 1,101
c	
c	        ### Scan the hydration shell weight ###
c	        ### Careful: range is model-dpdt ###
	        if(flagsolv.eq.0) then
	          c2 = 0.d0 + (real(ishell)-1)*4.d0/1.d2
	        else
	          c2 = 0.d-1 + (real(ishell)-1)*1.4/1.d2
	        endif
	        c2 = bulkdens*c2
c	
	        do 302 isample=1,nsample
	          intens = 0.d0
	          facexp = -pi*((4*pi)/3)**(3/2)
	          facexp = facexp*expq(isample)*expq(isample)
	          facexp = facexp*(rzero**2-averad**2)
	          do 303 j = 1,ncub
	            wj = qweight(1,j)
	            a = vc(isample,j) - c1*exp(facexp)*xc(isample,j)
     1	              + c2*wc(isample,j)
	            b = vs(isample,j) - c1*exp(facexp)*xs(isample,j)
     1	              + c2*ws(isample,j)
	            intens = intens + wj*(a*a+b*b)
303	          continue
	          proftmp(isample) = intens/wsum
302	        continue
c	
c	        ### Compute the deviation to exp ###
	        ma = 1
	        call svdfit(expq,profexp,sigexp,nsample,afit,ma,u,v,
     1	                    w,nsample,ma,chisq,proftmp)
c	
	        open(unit=1,file="isochi.txt",access='append',
     1	             status='unknown')
	          write(1,808) c1tmp,c2/bulkdens,chisq
	        close(unit=1)
c	
	        if(chisq.lt.chisqmin) then
c	
c	          ### Remember this one if better ###
	          chisqmin = chisq
	          c1best = c1tmp
	          c2best = c2
	          cbest = afit(1)
	          rzerobest = rzero
                endif
301	      continue
300	    continue
c	    write(6,*) " " 
c	
c	    ### Done with fitting mode ###
	  else
c	    ### Starting the non-fitting mode ###
c	
	    if(flagcuser.eq.1) then
	      rzerobest = c1user*averad
	      c1best = c1user
	      c2best = c2user*bulkdens
	      cbest = 1.d0
	    else
	      rzerobest = averad
	      c1best = 1.d0
	      c2best = bulkdens
	      cbest = 1.d0
	    endif
c	
c	    ### Done with non-fitting mode ###
	  endif
c	
c	  ### Write some info to the outside world ###
	  do 710 ilog = 1,2
	    if(ilog.eq.1) then
	      ulogtmp = 6
	    else
	      ulogtmp = ulog
	    endif
	    write(ulogtmp,805) c1best
	    write(ulogtmp,807) rzerobest
	    write(ulogtmp,806) c2best/bulkdens
	    if(flagfit.eq.1) then
	      write(ulogtmp,811) chisqmin
	    endif
710	  continue
c	
	  c1tmp = c1best
	  c1best = c1best*c1best*c1best*bulkdens
	  do 400 isample = 1,nsample
	    intens = 0.d0
	    facexp = -pi*((4*pi)/3)**(3/2)
            facexp = facexp*expq(isample)*expq(isample)
            facexp = facexp*(rzerobest**2-averad**2)
	    do 401 j = 1,ncub
	      wj = qweight(1,j)
              a = vc(isample,j) - c1best*exp(facexp)*xc(isample,j)
     1            + c2best*wc(isample,j)
              b = vs(isample,j) - c1best*exp(facexp)*xs(isample,j)
     1            + c2best*ws(isample,j)
              intens = intens + wj*(a*a+b*b)
401	    continue
	    profil(isample) = intens*cbest/wsum
400	  continue
c	
c	=============================================================
c		End of the subroutine				
c	=============================================================
	call time(charend)
	do 711 ilog = 1,2
	  if(ilog.eq.1) then
	    ulogtmp = 6
	  else
	    ulogtmp = ulog
	  endif
	  write(ulogtmp,6) charstart
	  write(ulogtmp,7) charend
711	continue
	call free(ptr_xsmap)
        call free(ptr_posx)
        call free(ptr_posy)
        call free(ptr_posz)
	close(unit=ulog)
c	
	return
	end
