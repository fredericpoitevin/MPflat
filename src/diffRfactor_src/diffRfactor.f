c
        integer iq,nq
        real*8 q(200),R(200),dRfact
        character*100 finput,line
c
        call getarg(1,finput)
c
        iq=0
        open(unit=1,file=finput,status='unknown')
100       read(1,'(a)',end=200) line
          iq = iq + 1
          read(line,*) q(iq),R(iq)
          goto 100
200       continue
        close(unit=1)
        nq=iq
c
        dRfact=0.d0
        do 300 iq = 2,nq
          dRfact = dRfact 
     1           + (q(iq)-q(iq-1))*R(iq)
300     continue
c       
        write(6,'(e12.5)') dRfact
c       
        stop
        end
