c     Classical 0(N3) hierarchical clustering based on R's hclust

      SUBROUTINE LLAHCLUSTN3(N,LEN,IOPT,IA,IB,CRIT,MEMBR,
     X     FLAG,SIMI,SIMI2)
c     Args
      INTEGER N, LEN, IOPT
      INTEGER IA(N), IB(N)
      LOGICAL FLAG(N)
      DOUBLE PRECISION CRIT(N), MEMBR(N), SIMI(LEN), SIMI2(LEN)
c     Var
      INTEGER IM, JM, I, NCL, J, IND, I2, J2, K, IND1, IND2, IND3
      DOUBLE PRECISION INF, SMAX, D12, SOUT
c     External function
      INTEGER IOFFST
c     
c     was 1D+20
      DATA INF/-1.D+300/;
c
c     unnecessary initialization of im jj jm to keep g77 -Wall happy
c
      IM = 0
      JM = 0
C
C     Initializations
C
      DO 10 I=1,N
C     We do not initialize MEMBR in order to be able to restart the
C     algorithm from a cut.
C     MEMBR(I)=1.
 10      FLAG(I)=.TRUE.
         NCL=N
C     

 800  CONTINUE
C
C      call dblepr("SIMI", 4, SIMI, N*(N-1)/2) 
C
C     Determine highest sim.
      SMAX=INF
      DO 80 I=1,N-1 
         IF (.NOT.FLAG(I)) GOTO 80
         DO 90 J=I+1,N
            IF (.NOT.FLAG(J)) GOTO 90
            IND=IOFFST(N,I,J)
            IF (SIMI(IND).GE.SMAX) THEN 
               SMAX=SIMI(IND)
               JM=J
               IM=I
            ENDIF
 90      CONTINUE
 80   CONTINUE     
C       
C      call dblepr("SMAX", 4, SMAX, 1) 
C      call intpr("IM", 2, IM, 1)
C      call intpr("JM", 2, JM, 1)
C
       NCL=NCL-1
C     
C     This allows an agglomeration to be carried out.
C     
       I2=MIN0(IM,JM)
       J2=MAX0(IM,JM)
       IA(N-NCL)=I2
       IB(N-NCL)=J2
       CRIT(N-NCL)=SMAX
       FLAG(J2)=.FALSE.
C     
C     Update similarities from new cluster.
C     
       SMAX=INF
       DO 50 K=1,N
C     WRITE(*,*)(MEMBR(III) , III=1,N)
          IF (.NOT.FLAG(K)) GOTO 50
          IF (K.EQ.I2) GOTO 50
C     WRITE(*,*)I2
          IF (I2.LT.K) THEN
             IND1=IOFFST(N,I2,K)
          ELSE
             IND1=IOFFST(N,K,I2)
          ENDIF
C     WRITE(*,*)J2
          IF (J2.LT.K) THEN
             IND2=IOFFST(N,J2,K)
          ELSE
             IND2=IOFFST(N,K,J2)
          ENDIF
          IND3=IOFFST(N,I2,J2)
          D12=SIMI(IND3)
C     
C     FISHER METHOD - IOPT=5.
C     
         IF (IOPT.EQ.5) THEN
            call fisher(MEMBR(K), MEMBR(I2), MEMBR(J2), 
     X           SIMI(IND1), SIMI(IND2), SOUT)
            SIMI(IND1) = SOUT
         ENDIF
C     
C     UNIFORM METHOD - IOPT=6.
C     
         IF (IOPT.EQ.6) THEN
            call uniform(K, IA, IB,  
     X           N - NCL, N, SIMI2, SOUT)
            SIMI(IND1) = SOUT
         ENDIF   
C     
C     NORMAL METHOD - IOPT=7.
C     
         IF (IOPT.EQ.7) THEN
            call normal(MEMBR(I2), MEMBR(J2), 
     X           SIMI(IND1), SIMI(IND2),SOUT)
            SIMI(IND1) = SOUT
         ENDIF    
C     
C     MAXIMUM METHOD - IOPT=8.
C     
         IF (IOPT.EQ.8) THEN
            call complete(MEMBR(I2), MEMBR(J2), 
     X           SIMI(IND1), SIMI(IND2), SOUT)      
            SIMI(IND1) = SOUT
         ENDIF
C
 50   CONTINUE
      MEMBR(I2)=MEMBR(I2)+MEMBR(J2)
C     
C     Repeat previous steps until N-1 agglomerations carried out.
C     
      IF (NCL.GT.1) GOTO 800
C     
C     Reverse criteria      
C     Criteria is a p-value (1 - similarity) 
C     
         DO 600 I=1,N-1
            CRIT(I)=1-CRIT(I)
 600     CONTINUE
      
      RETURN
      END
C
