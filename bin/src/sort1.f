        SUBROUTINE SORT1(QQ,I1,IN)
        SAVE
!
!       Finds ordering for data in QQ. No. of entries in QQ is I1.
!       IN(1),IN(2),.....are indices of QQ from lowest to highest.
!       NB Identical entries in QQ not ordered randomly.
!
!       Single precision
!
        
        DIMENSION QQ(I1),IN(I1),NB(I1+1),NX(I1+1),LH(I1+1),QB(I1+1),IB(I1+1)
!
        IF(I1 .EQ. 1) THEN
        IN(1)=1
        GO TO 99
        END IF
!
        QMN=QQ(1)
        QMX=QQ(1)
        DO 20 I=2,I1
        QMN=AMIN1(QMN,QQ(I))
20      QMX=AMAX1(QMX,QQ(I))
!
        DO 30 I=1,I1
        NB(I)=0
30      LH(I)=0
!
        Q1=QMN-1.E-5*(QMX-QMN)
        Q2=QMX+1.E-5*(QMX-QMN)
        DQ=(Q2-Q1)/I1
        DO 40 I=1,I1
        LB=(QQ(I)-Q1)/DQ+1
        NB(LB)=NB(LB)+1                      ! NO. IN BIN LB
        NX(I)=LH(LB)                         ! INDEX OF PREVIOUS ENTRY
40      LH(LB)=I                             ! INDEX OF LATEST ENTRY  
!
        J=0
        DO 60 LB=1,I1
        N1=NB(LB)
        IF(N1 .EQ. 0) GO TO 60

        N=0
        I=LH(LB)
50      N=N+1
        QB(N)=QQ(I)
        IB(N)=I
        I=NX(I)
        IF(I .GT. 0) GO TO 50
!
        M=0
52      QBMN=Q2
        DO 55 N=1,N1
        IF(QB(N) .LT. QBMN) THEN
        QBMN=QB(N)
        NM=N
        END IF
55      CONTINUE
        M=M+1
        J=J+1
        IN(J)=IB(NM)
        QB(NM)=Q2
        IF(M .LT. N1) GO TO 52
60      CONTINUE
!
99      RETURN 
        END
