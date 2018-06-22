c
c       Based on the code 'WaterProfile' by Patrice Koehl
c       
        integer natom
        real*8 coord(*)
        character label(*)*30
        character fname_pdb*100
        character fname_pdb*100
c
        pointer (ptr_coord,coord)
        pointer (ptr_label,label)
        common /arrays/ ptr_coord,ptr_label,natom
c       
        data atm/'N ','O ','N+','O-','C ','S '/
c
        call getarg(1,fname_pdb)
        call read_pdb(fname_pdb)
        do 100 i = 1,natom
          atom_type(i) = atm_solvtype(label(i))
100     continue
c
        stop
        end
c       ___________________
c
        subroutine read_pdb(fname_pdb)
c       
        integer l_in
        integer natom
        real*8 coord(*)
        character record*80
        character label(*)*30
        character fname_pdb*100
c
        pointer (ptr_coord,coord)
        pointer (ptr_label,label)
        common /arrays/ ptr_coord,ptr_label,natom
c
        l_in=11
        open(unit=l_in,file=fname_pdb,status='old',err=1000)
c       First count number of atoms
        natom = 0
50      read(l_in,1,end=75,err=1000) record
          if(record(1:6).ne.'ATOM  ') goto 50
          natom = natom + 1
          goto 50
75      continue
        write(6,*) "Number of atoms read",natom
c       Allocate space
        ptr_coord  = malloc(24*natom)
        ptr_label  = malloc(30*natom)
        rewind(l_in)
c       Second pass: store arrays
        natom = 0
100     read(l_in,1,end=200,err=1000) record
          if(record(1:6).ne.'ATOM  ') goto 100
          natom = natom + 1
          read(record,'(a30)') atom_label(natom)
          read(record(31:54),'(a24)') (atom_coord(3*(natom-1)+i),i=1,3)
        goto 100
200     continue
        close(unit=11)
c
        return
        end
c==========================================================================================
c==========================================================================================
c       atm_solvtype.f
c==========================================================================================
c==========================================================================================
c
c       This utility defines the type of an atom based on its name and the name
c       of the amino acid it belongs to
c
c       N (neutral)     1
c       O (neutral)     2
c       N (charged)     3
c       O (charged)     4
c       C               5
c       S               6
c
        function atm_solvtype(label)
c
        integer atm_solvtype
c
        character       label*30,name*1
        character       residue*3,longname*3
c
        name = label(14:14)
        residue=label(18:20)
        longname=label(14:16)
c
        if(name.eq.'N') then
                if(residue.eq.'ARG'.and.longname.eq.'NH1') then
                        atm_solvtype = 3
                elseif(residue.eq.'LYS'.and.longname.eq.'NZ') then
                        atm_solvtype = 3
                else
                        atm_solvtype = 1
                endif
        elseif(name.eq.'O') then
                if(residue.eq.'ASP'.and.longname.eq.'OD1') then
                        atm_solvtype = 4
                elseif(residue.eq.'GLU'.and.longname.eq.'OE1') then
                        atm_solvtype = 4
                else
                        atm_solvtype = 2
                endif
        elseif(name.eq.'C') then
                atm_solvtype = 5
        elseif(name.eq.'S') then
                atm_solvtype = 6
        else
                atm_solvtype = 0
        endif
c
        return
        end
