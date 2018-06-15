c	All routines about Surface definition are here
c	04/06/10	FP (adapted from P.Koehl)
c	
c____________________________________________________________________
	subroutine DefGrid
c	
	include 'param.h'
c	
	integer i,j,k
	integer nx,ny,nz,ix2,iy2,iz2
	integer ibeg(3),iend(3),ipas(3)
c	
	real*8 bulkdens,solvxs
	real*8 hx(ix),hy(iy),hz(iz)
	real*8 sizex,sizey,sizez
	real*8 map(ix,iy,iz),cell(6)
c	
1	format("ATOM    188  CG2 VAL    13    ",3f8.3)
c	
	common /smap/ ibeg,iend,ipas,cell,map
	common /posmap/ hx,hy,hz,nx,ny,nz
c	
c	
c ===================================================================
c		Set coordinates of grid points			
c ===================================================================
c	
	nx = ipas(1)-1
	ny = ipas(2)-1
	nz = ipas(3)-1
	sizex = cell(1)/ipas(1)
	sizey = cell(2)/ipas(2)
	sizez = cell(3)/ipas(3)
	do 100 i = 1,nx
	  hx(i) = (cell(1)/ipas(1))*(i-1) - cell(1)/2. + sizex/2
100	continue
	do 200 i = 1,ny
          hy(i) = (cell(2)/ipas(2))*(i-1) - cell(2)/2. + sizey/2
200     continue
	do 300 i = 1,nz
          hz(i) = (cell(3)/ipas(3))*(i-1) - cell(3)/2. + sizez/2
300     continue
c	do 400 i = 1,10
c	  do 401 j = 1,10
c	    do 402 k=  1,10
c	write(6,1) hx(k),hy(j),hz(i)
c402	    continue
c401	  continue
c400	continue
c ===================================================================
c		End of DefGrid					
c ===================================================================
	return
	end
c	
c____________________________________________________________________
	subroutine GenSurf(surf_type,R_probe,R_cush,grid_partition,
     1	                   grid_cushion,hx,hy,hz,nx,ny,nz)
c	
	include 'param.h'
c	
	integer natom
	integer i,j,k,kx,ky,kz
	integer ix2,iy2,iz2,ipoint
	integer istepx,istepy,istepz
	integer nx,ny,nz
	integer surf_type
	integer*1 grid_partition(nx,ny,nz)
	integer*1 grid_cushion(nx,ny,nz)
c	
	real*8 x,y,z
	real*8 atom_coord(ncortot),atom_radius(natot)
	real*8 R_probe,R_cush,rad,rad2,radc,radc2
	real*8 val,valx,valy,valz
	real*8 hx(nx),hy(ny),hz(nz)
c
	common /charges/ natom,atom_coord,atom_radius
c	
c ===================================================================
c		Initialise					
c ===================================================================
c	
	do 300 k = 1,nz
	  do 200 j = 1,ny
	    do 100 i = 1,nx
	      grid_partition(i,j,k) = 0
	      grid_cushion(i,j,k) = 1
100	    continue
200	  continue
300	continue
c	
c ===================================================================
c		Surface type 1: solvent accessible surface	
c ===================================================================
c	
	do 700 ipoint = 1,natom
	  x = atom_coord(3*(ipoint-1)+1)
	  y = atom_coord(3*(ipoint-1) +2)
	  z = atom_coord(3*(ipoint-1) +3)
	  call locate_cell(x,hx,ix2,nx)
          call locate_cell(y,hy,iy2,ny)
          call locate_cell(z,hz,iz2,nz)
	  rad = atom_radius(ipoint)
	  rad = rad + R_probe
	  radc = rad + R_cush
	  rad2 = rad*rad
	  radc2 = radc*radc
c	  call defstep(rad,hx,ix2,istepx,nx)
c          call defstep(rad,hy,iy2,istepy,ny)
c          call defstep(rad,hz,iz2,istepz,nz)
	  call defstep(radc,hx,ix2,istepx,nx)
          call defstep(radc,hy,iy2,istepy,ny)
          call defstep(radc,hz,iz2,istepz,nz)
	  do 600 k = 1,2*istepz + 1
c
                kz = iz2 - istepz + k -1
                kz=max(1,kz)
                kz=min(kz,nz)
                valz = (z - hz(kz))**2
c
                do 500 j = 1,2*istepy + 1
c
                   ky = iy2 - istepy + j -1
                   ky=max(1,ky)
                   ky=min(ky,ny)
                   valy=(y-hy(ky))**2
c
                   do 400 i = 1,2*istepx + 1
c
                           kx = ix2 - istepx + i -1
                           kx=max(1,kx)
                           kx=min(kx,nx)
                           if(grid_partition(kx,ky,kz).eq.0) then
                                valx = (x - hx(kx))**2
                                val=valx+valy+valz
c
c                               if vertex inside solvated ball ipoint,
c                               set bit zero to 1
c                               (all vertices with bit 0 set to 1
c                               belongs to the solvated volume of the molecule)
c
			        if(val.le.radc2) then
			            grid_cushion(kx,ky,kz) = 0
			        endif
                                if(val.le.rad2) then
                                    grid_partition(kx,ky,kz) = 1
			            grid_cushion(kx,ky,kz) = 1
                                endif
c
                           endif
c
400                continue
500             continue
600        continue
c
700     continue
c
c       We are done if we wanted to define the accessible surface area
c
        if(surf_type.eq.0) return
c	
c	
	return
	end
c____________________________________________________________________
c	
c       Version 1:      7/05/07
c       Author:         Patrice Koehl
c                       in collaboration with Marc Delarue
c
c       This subroutine uses bisection to search an ordered array
c
        subroutine locate_cell(x,array,idx,n)
c
c ====================================================================
c = Define variables
c ====================================================================
c
        integer i_low,i_high,i_mid
        integer n,idx
c
        real*8  x
        real*8  array(n)
c
c ====================================================================
c = Perform bisection search
c ====================================================================
c
        i_low = 0
        i_high = n+1
c
100     if(i_high-i_low.gt.1) then
c
                i_mid = (i_high+i_low)/2
                if(x.ge.array(i_mid)) then
                        i_low = i_mid
                else
                        i_high = i_mid
                endif
c
                goto 100
c
        endif
c
        if(x.eq.array(1)) then
                idx = 1
        elseif(x.eq.array(n)) then
                idx = n-1
        else
                idx = i_low
        endif
c	
	return
	end
c____________________________________________________________________
c       Version 1:      7/05/07
c       Author:         Patrice Koehl
c                       in collaboration with Marc Delarue
c
c       This subroutine finds the number of positions that need to be tested around a
c       given position to cover a given range
c
        subroutine defstep(val,array,ix,irange,n)
c
c ===================================================================
c = Define variables
c ===================================================================
c
        integer ix,irange,i,n,iup,ilow
c
        real*8  val,array(n)
c
        iup = 0
        ilow = 0
c
c ===================================================================
c = Find upper range
c ===================================================================
c
        do 100 i = ix+1,n
                if(array(i)-array(ix).gt.val) then
                        iup = i-ix
                        goto 200
                endif
100     continue
200     continue
c
c ===================================================================
c = Find lower range
c ===================================================================
c
        do 300 i = ix-1,1,-1
                if(array(ix)-array(i).gt.val) then
                        ilow = ix-i
                        goto 400
                endif
300     continue
400     continue
c
        irange=max(ilow,iup)
c
        return
        end
