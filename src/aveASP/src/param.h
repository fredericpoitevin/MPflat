C	param.h	Version 1 28/06/2010	Frederic Poitevin
c
c	This include file contains all the dimensions of the
c	arrays needed for all routines in Profmodel
c
c	Grid parameters
c	===============
c	
	integer ix,iy,iz
c	
	parameter	(ix = 200, iy = 200, iz	= 200)
c	
c	Protein parameters
c	==================
c	
	integer nresdef,naddrmax
c	
	parameter	(nresdef	= 24)
	parameter	(naddrmax	= 100)
c	
c	Atom parameters
c	===============
c	
	integer natot,ncortot,sftyp,naddamax
c	
	parameter	(natot		= 30000)
	parameter	(ncortot	= 3*natot)
	parameter	(sftyp		= 14)
	parameter	(naddamax	= 100)
c	
c	Scattering parameters
c	=====================
c	
	integer nqmax,ncubmax
c	
	parameter	(nqmax		= 1500)
	parameter	(ncubmax	= 900)
c	

