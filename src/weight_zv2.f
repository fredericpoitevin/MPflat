C CALCULATES DIFFERENCE STRUCTURE FACTORS AND WEIGHTS.
C SORT ORDER IN OUTPUT IS THE SAME AS IN FILE1
C CALC. MEAN SQUARE (AND NOT SQUARED MEAN) 
C INPUT: FILE1 REFERENCE: NORMALLY LIGHT (FOR LIGHT MINUS DARK AMPLITUDES)
C                         COLUNS H,K,L,F,SIGF REQUIRED, 6th COLUMN IGNORED 
C        FILE2          : IF HKL MATCH DIFFERENCE WILL BE WRITTEN OUT
C                         (NORMALLY DARK FOR LIGHT MINUS DARK)
C                         COLUNS H,K,L,F,SIGF REQUIRED, 6th COLUMN IGNORED 
C        FILE3          : PHASES MUST BE IN 6th COLUMN, WILL BE ONLY TRANSFERED TO
C                         OUTPUT FOR IMMEDIATE FFT
C        OUTFIL         : OUTPUT FILE NAME
C
C WEIGHT USED BY THIS PROGRAM:
C 
C                       1
C  W = ----------------------------------- 
C               F**2         SIG(F)**2
C       1  +  --------  +  ------------
C             < F**2 >     < SIG(F)**2 >
C 
C NOTE ! ZHONG REN USES THE SQUARE OF THE MEAN OF THE ABSOLUTE VALUES (E.G. <|F|>**2)
C        WHEREAS THIS PROGRAM USES THE MEAN OF SQUARES (OR THE MEAN SQUARE
C        WHICH IS MORE CORRECT IN A STATISTICAL MANNER) 
C STRUCTURE FACTORS WILL BE NORMALIZED BY THE MEAN WEIGHT IN ORDER TO PRESERVE
C ABSOLUTE SCALE OF A WEIGHTED MAP. ELECTRON DENSITY OF UNWEIGHTED MAPS CALCULATED 
C BY THE OUTPUT OF THIS PROGRAM HAS TO BE DIVIDED BY THE MEAN WEIGHT TO BE ON 
C ABSOLUTE SCALE
C BUG WITH OUTPUT DATA REMOVED
C MARIUS SCHMIDT SEPT/OCT 2000

           CHARACTER*60 FIL1,FIL2,FIL3,OFIL 
           DIMENSION FDARK(-100:100,-100:100,-100:100,3) 
C
C INITIALIZE FDARKSs PHASE FIELD (PHASE MUST BE DEFINED BY SFCALC PROGRAM
C PHASE ANGLE 0 ALLOWED THEREFORE INITALIZE WITH 999.9 
           DO 50 I=-100,100
            DO 50 J=-100,100
             DO 50 K=-100,100
 50           FDARK(I,J,K,3)=999.9
           WRITE(*,100) 'FIRST  DATA SET (REFERENCE)    : '
           READ(*,200) FIL1
           WRITE(*,*) FIL1
           WRITE(*,100) 'SECOND DATA SET                : '
           READ(*,200) FIL2
           WRITE(*,*) FIL2
           WRITE(*,100) 'PHASES DATA SET                : '
           READ(*,200) FIL3
           WRITE(*,*) FIL3
           WRITE(*,100) 'DIFFERENCE F AND WEIGHT OUTPUT : '
           READ(*,200) OFIL
           WRITE(*,*) OFIL
 100       FORMAT(A,$)
 200       FORMAT(A)

           OPEN(UNIT=20,FILE=FIL3,STATUS='OLD',ERR=1000)
           GOTO 1200
 1000      WRITE(*,*) ' PHASE FILE DOES NOT EXIST '
           GOTO 1003
 1001      WRITE(*,*) ' DATA SET 1 DOES NOT EXIST ' 
           GOTO 1003
 1002      WRITE(*,*) ' DATA SET 2 DOES NOT EXIST '
 1003      CONTINUE
           STOP 

 1200      CONTINUE
C PHASES ARE READ IN
           IP=0
 1300      CONTINUE
           READ(20,*,END=1400) IH,IK,IL,DM,DM,PHASE
           FDARK(IH,IK,IL,3)=PHASE
           IP=IP+1
           GOTO 1300
 1400      CONTINUE
C DEBUG    write(*,*) IP,IHMIN,IHMAX,KMIN,KMAX,LMIN,LMAX
C 
           I2=0 
           OPEN(UNIT=10,FILE=FIL1,STATUS='OLD',ERR=1001)
           OPEN(UNIT=11,FILE=FIL2,STATUS='OLD',ERR=1002)
 1500      CONTINUE
           READ(11,*,END=2000) IH,K,L,F1,SIGMA 
           FDARK(IH,K,L,1)=F1
           FDARK(IH,K,L,2)=SIGMA
           I2=I2+1
           GOTO 1500

 2000      CONTINUE
 
C LOOP OVER AMPLITUDES IN FILE 1 AND SEARCH MATCHING HKLs IN FILE 2.
C THEN CALCULATE MEAN DIFFERENCE AMPLITUDE AND THE MEAN OF THE
C SUM OF SIGMA1 AND SIGMA2.
           I1=0
           IC=0
           DFSQ=0.0
           S12SQ=0.0
           DF=0.0
           S12=0.0
           ADF=0.0
 2200      CONTINUE          
           READ(10,*,END=2500) IH,K,L,F1,SIGF1
           I1=I1+1
           IF (FDARK(IH,K,L,1).NE.0.0.AND.FDARK(IH,K,L,3).NE.999.9) THEN
                DFSQ=(F1-FDARK(IH,K,L,1))**2+DFSQ 
                S12SQ=SIGF1**2+FDARK(IH,K,L,2)**2+S12SQ
C BE CONSISTENT WITH ZHONG EASY TO CHANGE (SEE BELOW)
                DF=F1-FDARK(IH,K,L,1)+DF
                ADF=ABS(F1-FDARK(IH,K,L,1))+ADF
                S12=SQRT(SIGF1**2+FDARK(IH,K,L,2)**2)+S12
                IC=IC+1
             ENDIF
           GOTO 2200
 2500      CONTINUE           
           DFSQM=DFSQ/IC 
           S12SQM=S12SQ/IC
           DFM=DF/IC
           S12M=S12/IC
           ADFM=ADF/IC
C DEBUG    write(*,*) 'matching: ', IC

C DETERMINE MEAN WEIGHT TO PRESERVE ABSOLUTE SCALE
           IO=0
           WMEAN=0
           WZMEAN=0
           REWIND (10)
 2700      CONTINUE
           READ(10,*,END=2900) IH,K,L,F1,SIGF1
           IF (FDARK(IH,K,L,1).NE.0.0.AND.FDARK(IH,K,L,3).NE.999.9) THEN
C MARIUS WEIGHT MEAN SQUARE
                DF=F1-FDARK(IH,K,L,1)
C DON'T TAKE SQRT, IN NEXT LINE WILL BE SQUARED AGAIN 
                S12=SIGF1**2+FDARK(IH,K,L,2)**2
                W=1/(1+(DF**2/DFSQM)+(S12/S12SQM))
C ZHONGS WEIGTH SQUARED MEAN
                WZ=1/(1+(DF**2/ADFM**2)+(S12/S12M**2))
                WMEAN=WMEAN+W
                WZMEAN=WZMEAN+WZ
                IO=IO+1
            ENDIF
           GOTO 2700
 2900      CONTINUE
C           write(*,*)IO,WMEAN,WZMEAN
           WMEAN=WMEAN/IO
           WZMEAN=WZMEAN/IO
C 
           OPEN(UNIT=12,FILE=OFIL)
C LOOP AGAIN OVER ALL HKLs in FILE1 (PRESERVES SORT ORDER !!!!!), 
C CALCULATE WEIGHTS AND WRITE OUT DATA
           REWIND(10)
           DFDW=0
 3200      CONTINUE
           READ(10,*,END=3500) IH,K,L,F1,SIGF1 
           IF (FDARK(IH,K,L,1).NE.0.0.AND.FDARK(IH,K,L,3).NE.999.9) THEN
                DF=(F1-FDARK(IH,K,L,1))
                DFDW=DFDW+DF/WMEAN
C DON'T TAKE SQRT, WILL BE SQUARED IN NEXT LINE AGAIN
                S12=SIGF1**2+FDARK(IH,K,L,2)**2
C MARIUS WEIGHT MEAN SQUARE 
                W=1/(1+(DF**2/DFSQM)+(S12/S12SQM))
C ZHONGS WEIGTH SQUARED MEAN
                WZ=1/(1+(DF**2/ADFM**2)+(S12/S12M**2))
C DIVIDE BY AVERAGE WEIGHT TO PRESERVE SCALE 
                PHASE=FDARK(IH,K,L,3)
                IF (DF.LT.0.0) THEN
                  DF=ABS(DF)
                  PHASE=PHASE+180.0
                  IF (PHASE.GT.180.0) PHASE=PHASE-360.0
                ENDIF
                WRITE(12,5600) IH,K,L,DF/WMEAN,W,PHASE
            ENDIF
           GOTO 3200
 3500      CONTINUE
 5600      FORMAT(3I5,3F10.4)
           DFDW=DFDW/IO
C
           WRITE(*,6100) I1,I2,IP,IC,IO
 6100     FORMAT(' NUMBER OF HKL IN FILE1             : ',I10,/,
     *           ' NUMBER OF HKL IN FILE2             : ',I10,/,
     *           ' NUMBER OF PHASES IN FILE3          : ',I10,/,
     *           ' NUMBER OF MATCHING HKL             : ',I10,/,
     *           ' NUMBER OF DIFFERENCE-F WRITTEN OUT : ',I10 )
        WRITE(*,6200) DFM,DFDW,DFM**2,ADFM,ADFM**2,DFSQM,S12M,S12M**2,
     *                S12SQM,WMEAN,WZMEAN
 6200 FORMAT(' MEAN AMPLITUDE DIFFERENCE                      :  ',
     *                                                          G14.4,/,
     *    ' MEAN AMPLITUDE DIFFERENCE/<WEIGHT>          (M):  ',F10.4,/,
     *    ' MEAN AMPLITUDE DIFFERENCE SQUARED              :  ',G14.4,/,
     *    ' MEAN ABSOLUTE AMPLITUDE DIFFERENCE             :  ',F10.4,/,
     *    ' MEAN ABSOLUTE AMPLITUDE DIFFERENCE SQUARED  (Z):  ',F10.4,/,
     *    ' MEAN SQUARE AMPLITUDE DIFFERENCE            (M):  ',F10.4,/,
     *    ' MEAN SIGMA OF DIFFERENCE AMPLITUDES            :  ',F10.4,/,
     *    ' MEAN SIGMA OF DIFFERENCE AMPLITUDES SQUARED (Z):  ',F10.4,/,
     *    ' MEAN SQUARE SIGMA OF DIFFERENCE AMPLITUDES  (M):  ',F10.4,/,
     *    ' AVERAGE WEIGHT                              (M):  ',F10.4,/,
     *    ' AVERAGE ZHONG WEIGHT                        (Z):  ',F10.4 )    

           CLOSE(10)
           CLOSE(11)
           CLOSE(12)
           CLOSE(20)
           END


           
            
          

          
           
          
 
         
