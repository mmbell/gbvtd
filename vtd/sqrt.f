      SUBROUTINE CALCULATE_XY(IRAYPROC,NUMCOEFF,AZ,VE,XLLS,YLLS,
     $NDATA,WEIGHT)

C  This subroutine calculates the independent variables used in the
C  regression for all good data on a ring and stores them in arrays.

      IMPLICIT NONE
      INCLUDE 'vad_analysis_memory.for'
      INCLUDE 'vad_analysis_constants.for'
      INTEGER IRAYPROC,NUMCOEFF,NDATA,I,J,N
      REAL AZ(IRAYPROC),VE(IRAYPROC),XLLS(MAXCOEFFS,IRAYPROC),
     $YLLS(IRAYPROC),A,WEIGHT(IRAYPROC)
cc -- Changed 1-24-98 by wcl, NDATA is not initialized the first time
cc      DO I=1,NDATA
      do i=1,irayproc
      WEIGHT(I)=1.
      ENDDO
      NDATA=0
      N=NUMCOEFF/2
      DO 100 I=1,IRAYPROC
         IF(VE(I).LE.60..AND.VE(I).GE.-70.)THEN
            NDATA=NDATA+1
c            A=AZ(I)*RAD_CONVERSION
            A=AZ(I)
            XLLS(1,NDATA)=1.
            DO 50 J=1,N
               XLLS(2*J,NDATA)=SIN(FLOAT(J)*A)
               XLLS(2*J+1,NDATA)=COS(FLOAT(J)*A)
50          CONTINUE
            YLLS(NDATA)=VE(I)
         ENDIF
100   CONTINUE
      RETURN
      END
C==============================================================================

      SUBROUTINE LLS(NUMCOEFF,NDATA,EFFECTIVE_NDATA,X,Y,WEIGHT,S,
     $STAND_ERR,C,FLAG)

C  This subroutine performs a weighted linear least squares regression.

C  Input:
C  NDATA is the number of data points to fit.
C  EFFECTIVE_NDATA is the effective number of data points to use when
C  estimating the standard deviation about the regression.
C  It should equal NDATA except when groups of data have been weighted
C  to appear to be single data points.
C  Then it should be the number of groups of data.
C  NUMCOEFF is the number of independent variables (including the constant
C  term if used) in the fit.
C  Each independent variable has a coefficient to be determined.
C  X(I,J) is the value of the Ith independent variable for the Jth data point.
C  Y(J) is the value of the dependent variable for the Jth data point.
C  WEIGHT(J) is the weight assigned to the Jth data point.

C  Output:
C  C(I) is the coefficient for the Ith independent variable.
C  STAND_ERR(I) is the standard error of Ith coefficient.
C  S is the estimated standard deviation about the regression.
C  FLAG is returned .TRUE. if and only if there was sufficient data for
C  the regression to be performed.

      IMPLICIT NONE
      INCLUDE 'vad_analysis_memory.for'
      LOGICAL FLAG
      INTEGER NDATA,NUMCOEFF,IPOINT,ICOEFF,IROW,ICOL,EFFECTIVE_NDATA
      REAL S,SSRES,YREG
CC Change the array size from numcoeff to MAXCOEFFS 12/14/2000 by WCL
      REAL C(MAXCOEFFS),STAND_ERR(MAXCOEFFS)
      REAL Y(NDATA),WEIGHT(NDATA)
      REAL B(MAXCOEFFS)
      REAL X(MAXCOEFFS,NDATA)
      REAL A(MAXCOEFFS,MAXCOEFFS),
     $AINV(MAXCOEFFS,MAXCOEFFS)

C  Check whether there are enough data.
C  Since the standard error of the coefficients is being estimated, need
C  at least one more data than coefficients.
      IF(NDATA.LE.NUMCOEFF.OR.EFFECTIVE_NDATA.LE.NUMCOEFF)THEN
         FLAG=.FALSE.
         RETURN
      ELSE
         FLAG=.TRUE.
      ENDIF

C  Initialize the regression matrices.
      DO 10 ICOL=1,MAXCOEFFS
         DO 11 IROW=1,MAXCOEFFS
            A(IROW,ICOL)=0.
11       CONTINUE
10    CONTINUE
      DO 12 IROW=1,MAXCOEFFS
         B(IROW)=0.
         C(IROW)=0.  !reset C to 0.
12    CONTINUE

C  Loop through the data and accumulate the covariances in the regression
C  matrices.
      DO 20 IPOINT=1,NDATA
         DO 21 ICOL=1,NUMCOEFF
            DO 22 IROW=1,NUMCOEFF
               A(IROW,ICOL)=A(IROW,ICOL)+
     $         X(IROW,IPOINT)*X(ICOL,IPOINT)*WEIGHT(IPOINT)
22          CONTINUE
21       CONTINUE
         DO 23 IROW=1,NUMCOEFF
            B(IROW)=B(IROW)+X(IROW,IPOINT)*Y(IPOINT)*WEIGHT(IPOINT)
23       CONTINUE
20    CONTINUE

C  Find the inverse of A.
      CALL MINV(A,MAXCOEFFS,MAXCOEFFS,AINV,MAXCOEFFS,MAXCOEFFS,
     $NUMCOEFF)

C  Multiply AINV and B to get the coefficients C.
      CALL MAT_MULT(AINV,MAXCOEFFS,MAXCOEFFS,NUMCOEFF,NUMCOEFF,
     $B,MAXCOEFFS,1,NUMCOEFF,1,C,NUMCOEFF,1)

C  Calculate the sum of squares about the regression.
      SSRES=0.
      DO 200 IPOINT=1,NDATA
         YREG=0.
         DO 180 ICOEFF=1,NUMCOEFF
            YREG=YREG+C(ICOEFF)*X(ICOEFF,IPOINT)
180      CONTINUE
         SSRES=SSRES+WEIGHT(IPOINT)*(Y(IPOINT)-YREG)**2
200   CONTINUE

C  Calculate the estimated standard deviation about the regression.
      S=SQRT(SSRES/FLOAT(EFFECTIVE_NDATA-NUMCOEFF))

C  Calculate the standard error of the coefficients.
      DO 300 ICOEFF=1,NUMCOEFF
         STAND_ERR(ICOEFF)=S*SQRT(AINV(ICOEFF,ICOEFF))
300   CONTINUE

      RETURN
      END
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE MINV(A,ISIZE1,ISIZE2,B,ISIZE3,ISIZE4,M)

C  This subroutine calculates the inverse of a matrix.

C  Input:
C  A is the square matrix whose inverse is to be found.
C  M is the number of rows, as well as the number of columns, of A.
C  It must be less than or equal to MAXCOEFFS.
C  ISIZE1 is the first dimension of A in the calling program.
C  It must be greater than or equal to M.
C  ISIZE2 is the second dimension of A in the calling program.
C  It must be greater than or equal to M.
C  ISIZE3 is the first dimension of B in the calling program.
C  It must be greater than or equal to M.
C  ISIZE4 is the second dimension of B in the calling program.
C  It must be greater than or equal to M.

C  Output:
C  B is the inverse matrix of A.
C  It also has M rows and M columns.

      IMPLICIT NONE
      INCLUDE 'vad_analysis_memory.for'
      INCLUDE 'vad_analysis_constants.for'
      INTEGER ISIZE1,ISIZE2,ISIZE3,ISIZE4,M,I,J
      REAL A(ISIZE1,ISIZE2)
      REAL B(ISIZE3,ISIZE4)
      DOUBLEPRECISION DUMDP(MAXCOEFFS,2*MAXCOEFFS)

      DO 200 I=1,MAXCOEFFS
         DO 201 J=1,2*MAXCOEFFS
            DUMDP(I,J)=0.D0
201      CONTINUE
200   CONTINUE

      DO 101 I=1,M
         DO 102 J=1,M
            DUMDP(I,J)=DBLE(A(I,J))
102      CONTINUE
101   CONTINUE

      DO 103 I=1,M
         DO 104 J=M+1,2*M
            IF(I.EQ.J-M)THEN
               DUMDP(I,J)=1.D0
            ELSE
               DUMDP(I,J)=0.D0
            ENDIF
104      CONTINUE
103   CONTINUE

      CALL ROWRED_PP(DUMDP,MAXCOEFFS,2*MAXCOEFFS,M,2*M,0,
     $LOG_OUTPUT_UNIT)

      DO 105 I=1,M
         DO 106 J=1,M
            B(I,J)=SNGL(DUMDP(I,J+M))
106      CONTINUE
105   CONTINUE

      RETURN
      END
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE ROWRED_PP(A,ISIZE1,ISIZE2,M,N,IWRITE,LU)

C  This subroutine row reduces a matrix using partial pivoting.

C  Input and output:
C  A is the double precision matrix to be row reduced.
C  The row reduced form replaces the original matrix.

C  Input:
C  M is the number of rows in A.
C  N is the number of columns in A.
C  ISIZE1 is the first dimension of A in the calling program.
C  It must be greater than or equal to M.
C  ISIZE2 is the second dimension of A in the calling program.
C  It must be greater than or equal to N.
C  LU is the logical unit number to which written output is sent.
C  IWRITE controls what is written to LU as follows.
C     If IWRITE=0, nothing is written.
C     If IWRITE=1, the initial and row reduced matrices are written.
C     If IWRITE=2, the initial, final, and all intermediate matrices during
C     the row reduction are written.

      IMPLICIT NONE
      INTEGER IWRITE,LU,M,N,ISIZE1,ISIZE2,I,J,II,JJ,IP
      DOUBLEPRECISION A(ISIZE1,ISIZE2),W,X,Y,Z,STORE

      IF(IWRITE.GE.1)THEN
         CALL ROWRED_PP_WRITE(LU,A,ISIZE1,ISIZE2,M,N,1)
      ENDIF

C  Initialize the pivot row and the pivot column.
      I=1
      J=1

C  Check whether the pivot column is all zero.
104   DO 101 II=I,M
         IF(A(II,J).NE.0.D0)THEN
            GOTO 102
         ENDIF
101   CONTINUE

C  The pivot column is all zero.
C  Advance the pivot column.
      IF(J.LT.N)THEN
         J=J+1
         GOTO 104
      ELSE
         GOTO 103
      ENDIF

C  The pivot column is not all zero.
C  Interchange rows so that the term in the pivot column and pivot row
C  divided by the largest term in the row is maximum.
102   X=0.D0
      DO 105 II=I,M
         W=0.D0
         DO 1002 JJ=J,N
            IF(DABS(A(II,JJ)).GT.W)THEN
               W=DABS(A(II,JJ))
            ENDIF
1002     CONTINUE
         IF(W.NE.0.D0)THEN
            IF(DABS(A(II,J))/W.GT.X)THEN
               X=DABS(A(II,J)/W)
               IP=II
            ENDIF
         ENDIF
105   CONTINUE
110   IF(IP.NE.I)THEN
         DO 107 JJ=1,N
            STORE=A(I,JJ)
            A(I,JJ)=A(IP,JJ)
            A(IP,JJ)=STORE
107      CONTINUE
         IF(IWRITE.GE.2)THEN
            CALL ROWRED_PP_WRITE(LU,A,ISIZE1,ISIZE2,M,N,2)
         ENDIF
      ENDIF

C  Divide the terms in the pivot row by the term in the pivot row and
C  column.
106   Y=A(I,J)
      A(I,J)=1.D0
      IF(J+1.LE.N)THEN
         DO 108 JJ=J+1,N
            A(I,JJ)=A(I,JJ)/Y
108      CONTINUE
      ENDIF
      IF(IWRITE.GE.2)THEN
         CALL ROWRED_PP_WRITE(LU,A,ISIZE1,ISIZE2,M,N,2)
      ENDIF

C  Eliminate the terms in the pivot column in rows other than the pivot
C  row.
      DO 109 II=1,M
         IF(II.NE.I)THEN
            Z=A(II,J)
            A(II,J)=0.D0
            IF(J+1.LE.N)THEN
               DO 111 JJ=J+1,N
                  A(II,JJ)=A(II,JJ)-A(I,JJ)*Z
111            CONTINUE
            ENDIF
         ENDIF
109   CONTINUE
      IF(IWRITE.GE.2)THEN
         CALL ROWRED_PP_WRITE(LU,A,ISIZE1,ISIZE2,M,N,2)
      ENDIF
 
C  Row reduction for the current pivot row and pivot column is complete.
C  Advance the pivot row and pivot column.
      IF(I.LT.M.AND.J.LT.N)THEN
         I=I+1
         J=J+1
         GOTO 104
      ENDIF

C  The matrix is row reduced.
103   IF(IWRITE.EQ.1)THEN
         CALL ROWRED_PP_WRITE(LU,A,ISIZE1,ISIZE2,M,N,3)
      ENDIF

      RETURN
      END
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE ROWRED_PP_WRITE(LU,A,ISIZE1,ISIZE2,M,N,KODE)

C  This subroutine performs the writing of matrices for subroutine
C  ROWRED_PP.

      IMPLICIT NONE
      INTEGER LU,ISIZE1,ISIZE2,M,N,KODE,I,J
      DOUBLEPRECISION A(ISIZE1,ISIZE2)

      IF(KODE.EQ.1)THEN
         WRITE(LU,64)
64       FORMAT(' ORIGINAL MATRIX')
      ELSEIF(KODE.EQ.2)THEN
         WRITE(LU,65)
65       FORMAT(' MATRIX AFTER ROW REDUCTION OPERATION')
      ELSEIF(KODE.EQ.3)THEN
         WRITE(LU,66)
66       FORMAT(' ROW REDUCED MATRIX')
      ENDIF

      WRITE(LU,63)
63    FORMAT(' ')

      DO 101 I=1,M
         WRITE(LU,61)(A(I,J),J=1,N)
61       FORMAT(' ',8E16.8)
         WRITE(LU,63)
101   CONTINUE

      WRITE(LU,62)
62    FORMAT(//)

      RETURN
      END
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE MAT_MULT(A,ISIZE1,ISIZE2,MA,NA,B,ISIZE3,ISIZE4,MB,NB,
     $C,ISIZE5,ISIZE6)

C  This subroutine multiples two matrices.

C  Input:
C  A is the first matrix to be multiplied.
C  B is the second matrix to be multiplied.
C  MA is the number of rows in A.
C  NA is the number of columns in A.
C  MB is the number of rows in B.
C  MB must be equal to NA for matrix multiplication to be defined.
C  NB is the number of columns in B.
C  ISIZE1 is the first dimension of A in the calling program.
C  It must be greater than or equal to MA.
C  ISIZE2 is the second dimension of A in the calling program.
C  It must be greater than or equal to NA.
C  ISIZE3 is the first dimension of B in the calling program.
C  It must be greater than or equal to MB.
C  ISIZE4 is the second dimension of B in the calling program.
C  It must be greater than or equal to NB.
C  ISIZE5 is the first dimension of C in the calling program.
C  It must be greater than or equal to MA.
C  ISIZE6 is the second dimension of C in the calling program.
C  It must be greater than or equal to NB.

C  Output:
C  C is the resultant matrix.
C  It has MA rows and NB columns.

      IMPLICIT NONE
      INTEGER ISIZE1,ISIZE2,MA,NA,ISIZE3,ISIZE4,MB,NB,ISIZE5,ISIZE6,
     $IROW,ICOL,I
      REAL SUM
      REAL A(ISIZE1,ISIZE2),B(ISIZE3,ISIZE4),C(ISIZE5,ISIZE6)

      DO 2 ICOL=1,NB
         DO 3 IROW=1,MA
            SUM=0.
            DO 1 I=1,NA
               SUM=SUM+A(IROW,I)*B(I,ICOL)
1           CONTINUE
            C(IROW,ICOL)=SUM
3        CONTINUE
2     CONTINUE
      RETURN
      END
