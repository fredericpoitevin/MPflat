c	
c	28/06/10	FP
c	This file regroups some routines for map handling...
c	
c____________________________________________________________________
	subroutine writemap_cns(ibeg,iend,ipas,rmap,rave,rsigma,
     1				cell,nom_map,qform)
c	
	include 'param.h'
c	
	real*8	a,b,c,alphad,betad,gammad
	real*8	rmap(ix,iy,iz)
	real*8	cell(6),rsigma,rave
c
	integer i,na,nb,nc
	integer	ibeg(3),iend(3),ipas(3),unit,ksect,ntitle
	integer ia,amin,amax,ib,bmin,bmax,ic,cmin,cmax
c
	logical qform
c
	common	/unit_cell/a,b,c,alphad,betad,gammad
c
	character*80	ligne,title(100)
	character*64	nom_map
	character*3  	mode
c
c	write(6,*)' '
c	write(6,*)' Entering writemap subr.'
c	write(6,*)' name_file_map=',nom_map
c	write(6,*)' qform=',qform
	open(unit=1,file=nom_map,form='formatted',status='unknown')
C 
C write title 
C
	 NTITLE=1
         WRITE(1,'(/I8)') NTITLE 
         WRITE(1, '(A)',ERR=7) ' Calculated CNS Formatted Map' 
C 
C write sectioning information 
C
	amin=ibeg(1)
	amax=iend(1)
	bmin=ibeg(2)
	bmax=iend(2)
	cmin=ibeg(3)
	cmax=iend(3)
	na=ipas(1)
	nb=ipas(2)
	nc=ipas(3)
C	
        WRITE(1,'(9I8)',ERR=7) 
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX 
        WRITE(6,'(9I8)') 
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX 
C 
C write unit cell constants in angstroms and degrees 
C
         write(1,'(6E12.5)',ERR=7) (CELL(I),I=1,6) 
C
C write matrix mode 
C
	MODE='ZYX'
        WRITE(1,'(3A)',ERR=7) MODE 
C           
C write density matrix, c is slowest ("z-sections"). 
C
	ksect = 0
	do ic = 1,cmax-cmin+1
	  write(1,'(i8)',err=7) ksect
	  write(1,'(6e12.5)',err=7)
     &	  ((rmap(ia,ib,ic),ia=1,amax-amin+1),ib=1,bmax-bmin+1)
	  ksect = ksect + 1
	end do
c	
C write average and scale factor for map
C
        WRITE(1,'(I8)',ERR=7) KSECT
        WRITE(1,'(2(E12.4,1X))',ERR=7) RAVE, RSIGMA
C
	goto 8
c
7  	write(6,*)' Error in writing CNS map...'
C
8  	write(6,*)' Writing of CNS map completed...'
C
        RETURN
        END
c	
c____________________________________________________________________
	subroutine readmap_cns(ibeg,iend,ipas,rmap,rave,rsigma,
     1				cell,nom_map,qform)
c
	include 'param.h'
c
	real*8	a,b,c,alphad,betad,gammad
	real*8	rmap(ix,iy,iz)
	real*8	cell(6),rsigma,rave
c
	integer i,j,k,na,nb,nc
	integer	ibeg(3),iend(3),ipas(3),unit,ksect,ntitle
	integer ia,amin,amax,ib,bmin,bmax,ic,cmin,cmax
c
	logical qform
c
	common	/unit_cell/a,b,c,alphad,betad,gammad
c
	character*80	ligne,title(100)
	character*64	nom_map
	character*3  	mode
c
c	write(6,*)' '
c	write(6,*)' Entering readmap subr.'
c	write(6,*)' name_file_map=',nom_map
c	write(6,*)' qform=',qform
	open(unit=1,file=nom_map,form='formatted',status='old')
C 
C read title 
C
      IF (QFORM) THEN 
         READ(1,'(/I8)') NTITLE 
C         WRITE(6,'(a,I8)')' Ntitle=',NTITLE 
         IF (NTITLE .LE. 0) THEN 
            READ(1, '(A)',END=6,ERR=7) 
         ELSE 
            DO J = 1, NTITLE 
               TITLE(J) = ' ' 
               READ(1, '(A)',END=6,ERR=7) TITLE(J) 
            ENDDO 
         ENDIF 
C            write(6,'(A)')(title(k),k=1,ntitle)
      ELSE 
         DO J=1,100
            TITLE(J)=' ' 
         END DO 
         READ(1,END=6,ERR=7) NTITLE,(TITLE(J)(1:80),J=1,NTITLE) 
      END IF 
C 
C read sectioning information 
C
      IF (QFORM) THEN 
        READ(1,'(9I8)',END=6,ERR=7) 
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX 
      ELSE 
        READ(1,END=6,ERR=7) 
     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX 
      END IF 
C        WRITE(6,'(9I8)') 
C     &  NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX 
C 
C read unit cell constants in angstroms and degrees 
C
      IF (QFORM) THEN 
         READ(1,'(6E12.5)',END=6,ERR=7) (CELL(I),I=1,6) 
      ELSE 
         READ(1,END=6,ERR=7) (CELL(I),I=1,6) 
      END IF 
C
	a=cell(1)
	b=cell(2)
	c=cell(3)
	alphad=cell(4)
	betad =cell(5)
	gammad=cell(6)
C 
C read matrix mode 
C
      IF (QFORM) THEN 
        READ(1,'(3A)',END=6,ERR=7) MODE 
      ELSE 
        READ(1,END=6,ERR=7) MODE 
      END IF 
      IF (MODE.NE.'ZYX') THEN 
         WRITE(6,*)' ERROR IN MODE OF SECTIONING:',MODE
         STOP
      END IF 
C           
C read density matrix, c is slowest ("z-sections"). 
C
	do ic = 1,cmax-cmin+1
	  if(qform) then
	    read(1,'(i8)',end=6,err=7) ksect
	    read(1,'(6e12.5)',end=6,err=7)
     &	    ((rmap(ia,ib,ic),ia=1,amax-amin+1),ib=1,bmax-bmin+1)
	  else
	    read(1,end=6,err=7) ksect
	    read(1,end=6,err=7)
     &	    ((rmap(ia,ib,ic),ia=1,amax-amin+1),ib=1,bmax-bmin+1)
	  endif
	end do
c	
C read average and scale factor for map
C
      IF (QFORM) THEN 
        READ(1,'(I8)',END=6,ERR=7) KSECT
        READ(1,'(2(E12.4,1X))',END=6,ERR=7) RAVE, RSIGMA
      ELSE 
        READ(1,END=6,ERR=7) KSECT
        READ(1,END=6,ERR=7) RAVE, RSIGMA
      END IF    
C
6	ibeg(1)=amin
	iend(1)=amax
	ibeg(2)=bmin
	iend(2)=bmax
	ibeg(3)=cmin
	iend(3)=cmax
C
	if(ibeg(1).lt.-ix/2.or.iend(1).gt.ix/2) then
		write(6,*)' Error in dimensioning map (nx)',ix
		stop
	elseif(ibeg(2).lt.-iy/2.or.iend(2).gt.iy/2) then
		write(6,*)' Error in dimensioning map (ny)',iy
		stop
	elseif(ibeg(3).lt.-iz/2.or.iend(3).gt.iz/2) then
		write(6,*)' Error in dimensioning map (nz)',iz
		stop
	endif
C
	ipas(1)=na
	ipas(2)=nb
	ipas(3)=nc
C	
	goto 8
c
7  	write(6,*)' Error in reading CNS map...'
C
8  	continue
C	write(6,*)' Reading of CNS map completed...'
C
        RETURN
        END
c
c____________________________________________________________________
c       ReadMap_OpenDX.f
c
c	22/06/20	FP
c	adapted from P.Koehl, in collaboration with M.Delarue
c	
c       Purpose:        This subroutine read a map in OpenDX format
c
c
        subroutine readmap_dx(ibeg,iend,ipas,rmap,cell,nom_map)
c	
	include 'param.h'
c	
	integer ibeg(3),iend(3),ipas(3)
        integer i,j,k
        integer nx,ny,nz,nvol
	integer tmp1,tmp2
c	
        character title*72,dieze*1,nom_map*64
c	
        real*8  origin(3),dx(3),dy(3),dz(3),cell(6)
        real*8  rmap(ix,iy,iz)
        real*8  sizex,sizey,sizez
c	
1       format('# Data from GenSolv V 0.1',/,'#')
2       format('#')
3       format('object 1 class gridpositions counts',5x,i3,5x,i3,5x,i3)
c4       format('origin ',f15.6,1x,f15.6,1x,f15.6)
4	format('origin ',1x,f14.8,1x,f14.8,1x,f14.8)
c5       format('delta ',f15.6,1x,f15.6,1x,f15.6)
5	format('delta ',f16.8,1x,i1,1x,i1)
6       format('object 2 class gridconnections counts',5x,i3,5x,i3,
     1	5x,i3)
7       format('object 3 class array type double rank 0 items',i9,
     1           ' follows')
c8       format(1x,1pe15.6,1x,1pe15.6,1x,1pe15.6)
8	format(e16.5,e16.5,e16.53)
9       format('attribute "dep" string "positions"')
10      format(
     1  'object "regular positions regular connections" class field')
11      format('component "positions" value 1')
12      format('component "connections" value 2')
13      format('component "data" value 3')
14      format('# ',a72)
15	format('delta ',1x,i1,f16.8,1x,i1)
16	format('delta ',1x,i1,1x,i1,f16.8)
c	
c	Skip eventual header...
c	
	open(unit=11,file=nom_map,status='unknown')
c100	  read(1,'(a1)') dieze
c	  if(dieze.eq.' '.or.dieze.eq.'#') goto 100
c200	  continue
c	
	  read(11,3) nx,ny,nz
	  read(11,4) origin(1),origin(2),origin(3)
c	  read(11,5) dx(1),dx(2),dx(3)
	  read(11,5) dx(1),tmp1,tmp2
c	  read(11,15) dy(1),dy(2),dy(3)
	  read(11,15) tmp1,dy(2),tmp2
c	  read(11,16) dz(1),dz(2),dz(3)
	  read(11,16) tmp1,tmp2,dz(3)
c	
	  read(11,6) nx,ny,nz
	  read(11,7) nvol
c	
	  if(nvol.ne.nx*ny*nz) then
	    write(6,*)'... Problem reading OpenDX map'
	    stop
	  endif
c	
	  ibeg(1) = origin(1)/dx(1)
	  ibeg(2) = origin(2)/dy(2)
	  ibeg(3) = origin(3)/dz(3)
	  iend(1) = ibeg(1)+nx
	  iend(2) = ibeg(2)+ny
	  iend(3) = ibeg(3)+nz
	  ipas(1) = nx
	  ipas(2) = ny
	  ipas(3) = nz
	  cell(1) = nx*dx(1)
	  cell(2) = ny*dy(2)
	  cell(3) = nz*dz(3)
	  cell(4) = 90.
	  cell(5) = 90.
	  cell(6) = 90.
c	
	  read(11,8) (((rmap(i,j,k),k=1,nz),j=1,ny),i=1,nx)
c	  read(11,8) (((rmap(i,j,k),k=ibeg(3),iend(3)),
c     1	j=ibeg(2),iend(2)),i=ibeg(1),iend(1))
	close(unit=11)
	write(6,*) rmap(1,1,1)
c	
	write(6,*)'     Done reading OpenDX map'
        return
        end
c	
