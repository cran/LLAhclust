C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                            C
C  HIERARCHICAL CLUSTERING using (user-specified) criterion. C
C                                                            C
C  Modifications of the original R version (package stats) C
C  Here we work with similarities instead of dissimilarities C 
C  and therefore we maximize instead of minimizing C
C  Parameters:                                               C
C                                                            C
C  N                 the number of points being clustered    C
C  SIMI(LEN)         similarities in lower half diagonal  C
C                    storage; LEN = N.N-1/2,                 C
C  IOPT              clustering criterion to be used,        C
C  IA, IB, CRIT      history of agglomerations; dimensions   C
C                    N, first N-1 locations only used,       C
C  MEMBR, NN, SIMNN  vectors of length N, used to store      C
C                    cluster cardinalities, current nearest  C
C                    neighbour, and the similarity assoc. C
C                    with the latter.                        C
C                    MEMBR must be initialized by R to the   C
C                    default of  rep(1, N)                   C
C  FLAG              boolean indicator of agglomerable obj./ C
C                    clusters.                               C
C                                                            C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       C
C  Modifications for R: Ross Ihaka, Dec 1996                 C
C                       Fritz Leisch, Jun 2000               C
C  all vars declared:   Martin Maechler, Apr 2001            C
C  PR#4195 fixed by BDR Nov 2003                             C
C                                                                               C
C  LLA modification: Ivan Kojadinovic, Mar 2007. C
C------------------------------------------------------------C

      SUBROUTINE LLAHCLUSTNN(N,LEN,IOPT,IA,IB,CRIT,MEMBR,NN,SIMNN,
     X     FLAG,SIMI,EPSILON)
c     Args
      INTEGER N, LEN, IOPT
      INTEGER IA(N), IB(N), NN(N)
      LOGICAL FLAG(N)
      DOUBLE PRECISION CRIT(N), MEMBR(N), SIMI(LEN), SIMNN(N)
      DOUBLE PRECISION EPSILON
c     Var
      INTEGER IM, JJ, JM, I, NCL, J, IND, I2, J2, K, IND1, IND2, IND3
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
      JJ = 0
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
C     Carry out an agglomeration - first create list of NNs
C     Note NN and SIMNN are nearest neighbours and its "distance"
C     TO THE RIGHT of I.
C
      DO 80 I=1,N-1
         SMAX=INF
         DO 90 J=I+1,N
            IND=IOFFST(N,I,J)
            IF (SIMI(IND).LE.SMAX) GOTO 90
            SMAX=SIMI(IND)
            JM=J
 90      CONTINUE
         NN(I)=JM
         SIMNN(I)=SMAX
 80   CONTINUE
C     
 800  CONTINUE
C     Next, determine highest sim. using list of NNs
      SMAX=INF
      DO 700 I=1,N-1
          IF (.NOT.FLAG(I)) GOTO 700
          IF (SIMNN(I).LE.SMAX) GOTO 700
          SMAX=SIMNN(I)
          IM=I
          JM=NN(I)
 700   CONTINUE
       
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
C     LLA METHOD - IOPT=1.
C     
         IF (IOPT.EQ.1) THEN
            SIMI(IND1)=-EPSILON*LOG(MEMBR(I2)+MEMBR(J2))
     X           +MAX(SIMI(IND1)+EPSILON*LOG(MEMBR(I2)),SIMI(IND2)
     X           +EPSILON*LOG(MEMBR(J2)))
         ENDIF
C     
C     TIPPETT METHOD - IOPT=2.
C
         IF (IOPT.EQ.2) THEN
            call tippett(MEMBR(I2), MEMBR(J2), 
     X           SIMI(IND1), SIMI(IND2), EPSILON, SOUT)
            SIMI(IND1) = SOUT
         ENDIF  
C     
C     AVERAGE LINK (OR GROUP AVERAGE) METHOD - IOPT=3.
C     
         IF (IOPT.EQ.3) THEN
            SIMI(IND1)=(MEMBR(I2)*SIMI(IND1)+MEMBR(J2)*SIMI(IND2))/
     X           (MEMBR(I2)+MEMBR(J2))
         ENDIF
C     
C     COMPLETE LINK METHOD - IOPT=4.
C     
         IF (IOPT.EQ.4) THEN
            SIMI(IND1)=MIN(SIMI(IND1),SIMI(IND2))
         ENDIF
C
C     
 50   CONTINUE
      MEMBR(I2)=MEMBR(I2)+MEMBR(J2)
C     
C     Update list of NNs
C     
      DO 900 I=1,N-1
         IF (.NOT.FLAG(I)) GOTO 900
C     (Redetermine NN of I:)
         SMAX=INF
         DO 870 J=I+1,N
            IF (.NOT.FLAG(J)) GOTO 870
            IND=IOFFST(N,I,J)
            IF (SIMI(IND).LE.SMAX) GOTO 870
            SMAX=SIMI(IND)
            JJ=J
 870     CONTINUE
         NN(I)=JJ
         SIMNN(I)=SMAX
 900  CONTINUE
C     
C     Repeat previous steps until N-1 agglomerations carried out.
C     
      IF (NCL.GT.1) GOTO 800
C     
C     Reverse criteria   
C     
      
      IF (IOPT.EQ.1) THEN
         SMAX = CRIT(1) + 1
         DO 500 I=1,N-1
            CRIT(I)=SMAX-CRIT(I)
 500     CONTINUE
      ELSE
C     
C     Criteria is a p-value (1 - similarity) 
C     
         DO 600 I=1,N-1
            CRIT(I)=1-CRIT(I)
 600     CONTINUE
      ENDIF
      
      RETURN
      END
C
C
      INTEGER FUNCTION IOFFST(N,I,J)
C     Map row I and column J of upper half diagonal symmetric matrix
C     onto vector.
      INTEGER N,I,J
      IOFFST=J+(I-1)*N-(I*(I+1))/2
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C  Given a HIERARCHIC CLUSTERING, described as a sequence of    C
C  agglomerations, prepare the seq. of aggloms. and "horiz."    C
C  order of objects for plotting the dendrogram using S routine C
C  'plclust'.                                                   C
C                                                               C
C  Parameters:                                                  C
C                                                               C
C  IA, IB:       vectors of dimension N defining the agglomer-  C
C                 ations.                                       C
C  IIA, IIB:     used to store IA and IB values differently     C
C                (in form needed for S command `plclust`        C
C  IORDER:       "horiz." order of objects for dendrogram       C
C                                                               C
C  F. Murtagh, ESA/ESO/STECF, Garching, June 1991               C
C                                                               C
C  HISTORY                                                      C
C                                                               C
C  Adapted from routine HCASS, which additionally determines    C
C   cluster assignments at all levels, at extra comput. expense C
C                                                               C
C---------------------------------------------------------------C
      SUBROUTINE HCASS2(N,IA,IB,IORDER,IIA,IIB)
c Args
      INTEGER N,IA(N),IB(N),IORDER(N),IIA(N),IIB(N)
c Var
      INTEGER I, J, K, K1, K2, LOC
C
C     Following bit is to get seq. of merges into format acceptable to plclust
C     I coded clusters as lowest seq. no. of constituents; S's `hclust' codes
C     singletons as -ve numbers, and non-singletons with their seq. nos.
C
      DO 912 I=1,N
         IIA(I)=IA(I)
         IIB(I)=IB(I)
  912 CONTINUE
      DO 915 I=1,N-2
C        In the following, smallest (+ve or -ve) seq. no. wanted
         K=MIN(IA(I),IB(I))
         DO 913 J=I+1, N-1
            IF(IA(J).EQ.K) IIA(J)=-I
            IF(IB(J).EQ.K) IIB(J)=-I
  913    CONTINUE
  915 CONTINUE
      DO 916 I=1,N-1
         IIA(I)=-IIA(I)
         IIB(I)=-IIB(I)
  916 CONTINUE
      DO 917 I=1,N-1
         IF (IIA(I).GT.0.AND.IIB(I).LT.0) THEN
            K = IIA(I)
            IIA(I) = IIB(I)
            IIB(I) = K
         ENDIF
         IF (IIA(I).GT.0.AND.IIB(I).GT.0) THEN
            K1 = MIN(IIA(I),IIB(I))
            K2 = MAX(IIA(I),IIB(I))
            IIA(I) = K1
            IIB(I) = K2
         ENDIF
  917 CONTINUE
C
C
C     NEW PART FOR `ORDER'
C
      IORDER(1) =IIA(N-1)
      IORDER(2) =IIB(N-1)
      LOC=2
      DO 175 I=N-2,1,-1
        DO 169 J=1,LOC
          IF(IORDER(J).EQ.I) THEN
C           REPLACE IORDER(J) WITH IIA(I) AND IIB(I)
            IORDER(J)=IIA(I)
            IF (J.EQ.LOC) THEN
                LOC=LOC+1
                IORDER(LOC)=IIB(I)
                GOTO 171
            ENDIF
            LOC=LOC+1
            DO 95 K=LOC,J+2,-1
               IORDER(K)=IORDER(K-1)
  95        CONTINUE
            IORDER(J+1)=IIB(I)
            GOTO 171
          ENDIF
 169    CONTINUE
C       SHOULD NEVER REACH HERE
 171    CONTINUE
 175  CONTINUE
C
C
      DO 181 I=1,N
         IORDER(I) = -IORDER(I)
 181  CONTINUE
C
C
      RETURN
      END
