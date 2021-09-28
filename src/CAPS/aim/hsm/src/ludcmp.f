
      SUBROUTINE LUDCMP(NSIZ,N,A)
C     *******************************************************
C     *                                                     *
C     *   Factors a full NxN matrix into an LU form.        *
C     *   Subr. BAKSUB can back-substitute it with some RHS.*
C     *   Assumes matrix is non-singular...                 *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   A is the matrix...                                *
C     *     ...replaced with its LU factors.                *
C     *                                                     *
C     *                              Mark Drela  1988       *
C     *******************************************************
      IMPLICIT NONE

      INTEGER NSIZ,N
      REAL A(NSIZ,NSIZ)

      INTEGER I,J,K
      REAL SUM, DUM
C
      DO J=1, N
        DO I=1, J-1
          SUM = A(I,J)
          DO K=1, I-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
        ENDDO
C
        DO I=J, N
          SUM = A(I,J)
          DO K=1, J-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
        ENDDO
C
        DUM = 1.0/A(J,J)
        DO I=J+1, N
          A(I,J) = A(I,J)*DUM
        ENDDO
      ENDDO
C
      RETURN
      END ! LUDCMP



      SUBROUTINE LUDCMPT(NSIZ,N,A)
      IMPLICIT NONE

      INTEGER NSIZ,N
      REAL A(NSIZ,NSIZ)

      INTEGER I,J,K
      REAL SUM, DUM
C
      DO J=1, N
        DO I=1, J-1
          SUM = A(J,I)
          DO K=1, I-1
            SUM = SUM - A(K,I)*A(J,K)
          ENDDO
          A(J,I) = SUM
        ENDDO
C
        DO I=J, N
          SUM = A(J,I)
          DO K=1, J-1
            SUM = SUM - A(K,I)*A(J,K)
          ENDDO
          A(J,I) = SUM
        ENDDO
C
        DUM = 1.0/A(J,J)
        DO I=J+1, N
          A(J,I) = A(J,I)*DUM
        ENDDO
      ENDDO
C
      RETURN
      END ! LUDCMPT



      SUBROUTINE BAKSUB(NSIZ,N,A,B)
      IMPLICIT NONE

      INTEGER NSIZ,N
      REAL A(NSIZ,NSIZ), B(NSIZ)
C
      INTEGER I,J,K,II
      REAL SUM, DUM

      II = 0
      DO I=1, N
        SUM = B(I)
        IF(II.NE.0) THEN
         DO J=II, I-1
           SUM = SUM - A(I,J)*B(J)
         ENDDO
        ELSE IF(SUM.NE.0.0) THEN
         II = I
        ENDIF
        B(I) = SUM
      ENDDO
C
      DO I=N, 1, -1
        SUM = B(I)
        IF(I.LT.N) THEN
         DO J=I+1, N
           SUM = SUM - A(I,J)*B(J)
         ENDDO
        ENDIF
        B(I) = SUM/A(I,I)
      ENDDO
C
      RETURN
      END ! BAKSUB



      SUBROUTINE BAKSUBB(NSIZ,N,A,B)
      IMPLICIT NONE

      INTEGER NSIZ,N
      REAL A(NSIZ,NSIZ), B(NSIZ,NSIZ)

      INTEGER I,J,K,L
      REAL SUM, DUM

      DO L = 1, N

      DO I = 1, N
        SUM = B(I,L)
        DO J = 1, I-1
           SUM = SUM - A(I,J)*B(J,L)
        ENDDO
        B(I,L) = SUM
      ENDDO

      DO I = N, 1, -1
        SUM = B(I,L) 
        DO J = I+1, N
           SUM = SUM - A(I,J)*B(J,L)
        ENDDO
        B(I,L) = SUM/A(I,I)
      ENDDO

      ENDDO
C
      RETURN
      END ! BAKSUBB



      SUBROUTINE BAKSUBBT(NSIZ,N,A,B)
      IMPLICIT NONE

      INTEGER NSIZ,N
      REAL A(NSIZ,NSIZ), B(NSIZ,NSIZ)

      INTEGER I,J,K,L
      REAL SUM, DUM

      DO L = 1, N

      DO I = 1, N
        SUM = B(L,I)
        DO J = 1, I-1
           SUM = SUM - A(J,I)*B(L,J)
        ENDDO
        B(L,I) = SUM
      ENDDO

      DO I = N, 1, -1
        SUM = B(L,I) 
        DO J = I+1, N
           SUM = SUM - A(J,I)*B(L,J)
        ENDDO
        B(L,I) = SUM/A(I,I)
      ENDDO

      ENDDO

      RETURN
      END ! BAKSUBBT



      SUBROUTINE LUDCMPI(NSIZ,N,A,INDX)
C     *******************************************************
C     *                                                     *
C     *   Factors a full NxN matrix into an LU form.        *
C     *   Subr. BAKSUB can back-substitute it with some RHS.*
C     *   Assumes matrix is non-singular...                 *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   A is the matrix...                                *
C     *     ...replaced with its LU factors.                *
C     *                                                     *
C     *                              Mark Drela  1988       *
C     *******************************************************
      IMPLICIT NONE

      INTEGER NSIZ,N
      REAL A(NSIZ,NSIZ)
      INTEGER INDX(NSIZ)
C
      INTEGER NVX
      PARAMETER (NVX=1000)
      REAL VV(NVX)

      INTEGER I,J,K,IMAX
      REAL AAMAX, SUM, DUM
C
      IF(N.GT.NVX) STOP 'LUDCMPI: Array overflow. Increase NVX.'
C
      DO I=1, N
        AAMAX = 0.
        DO J=1, N
          AAMAX = MAX( ABS(A(I,J)) , AAMAX )
        ENDDO
        VV(I) = 1.0/AAMAX
      ENDDO
C
      DO J=1, N
        DO I=1, J-1
          SUM = A(I,J)
          DO K=1, I-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
        ENDDO
C
        AAMAX = 0.
        DO I=J, N
          SUM = A(I,J)
          DO K=1, J-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
C
          DUM = VV(I)*ABS(SUM)
          IF(DUM.GE.AAMAX) THEN
           IMAX = I
           AAMAX = DUM
          ENDIF
        ENDDO
C
        IF(J.NE.IMAX) THEN
         DO K=1, N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
         ENDDO
         VV(IMAX) = VV(J)
        ENDIF
C
        INDX(J) = IMAX

        IF(J.NE.N) THEN
         DUM = 1.0/A(J,J)
         DO I=J+1, N
           A(I,J) = A(I,J)*DUM
         ENDDO
        ENDIF
C
      ENDDO
C
      RETURN
      END ! LUDCMPI



      SUBROUTINE BAKSUBI(NSIZ,N,A,INDX,B)
      IMPLICIT NONE

      INTEGER NSIZ,N
      REAL A(NSIZ,NSIZ), B(NSIZ)
      INTEGER INDX(NSIZ)

      INTEGER I,J,II,LL
      REAL SUM
C
      II = 0
      DO I=1, N
        LL = INDX(I)
        SUM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
         DO J=II, I-1
           SUM = SUM - A(I,J)*B(J)
         ENDDO
        ELSE IF(SUM.NE.0.0) THEN
         II = I
        ENDIF
        B(I) = SUM
      ENDDO
C
      DO I=N, 1, -1
        SUM = B(I)
        IF(I.LT.N) THEN
         DO J=I+1, N
           SUM = SUM - A(I,J)*B(J)
         ENDDO
        ENDIF
        B(I) = SUM/A(I,I)
      ENDDO
C
      RETURN
      END ! BAKSUBI




      SUBROUTINE BAKSUBBI(NSIZ,N,A,INDX,B)
      IMPLICIT NONE

      INTEGER NSIZ,N
      REAL A(NSIZ,NSIZ), B(NSIZ,NSIZ)
      INTEGER INDX(NSIZ)

      INTEGER I,J,II,LL, L
      REAL SUM
C
      DO L = 1, N

      II = 0
      DO I=1, N
        LL = INDX(I)
        SUM = B(LL,L)
        B(LL,L) = B(I,L)
        IF(II.NE.0) THEN
         DO J=II, I-1
           SUM = SUM - A(I,J)*B(J,L)
         ENDDO
        ELSE IF(SUM.NE.0.0) THEN
         II = I
        ENDIF
        B(I,L) = SUM
      ENDDO
C
      DO I=N, 1, -1
        SUM = B(I,L)
        IF(I.LT.N) THEN
         DO J=I+1, N
           SUM = SUM - A(I,J)*B(J,L)
         ENDDO
        ENDIF
        B(I,L) = SUM/A(I,I)
      ENDDO
C
      ENDDO

      RETURN
      END ! BAKSUBBI



      SUBROUTINE ZLUDCMP(NSIZ,N,A)
C     *******************************************************
C     *                                                     *
C     *   Factors a full NxN matrix into an LU form.        *
C     *   Subr. BAKSUB can back-substitute it with some RHS.*
C     *   Assumes matrix is non-singular...                 *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   A is the matrix...                                *
C     *     ...replaced with its LU factors.                *
C     *                                                     *
C     *                              Mark Drela  1988       *
C     *******************************************************
      IMPLICIT NONE

      INTEGER NSIZ,N
      COMPLEX*16 A(NSIZ,NSIZ)

      INTEGER I,J,K
      COMPLEX*16 SUM, AINV
C
      DO J=1, N
        DO I=1, J-1
          SUM = A(I,J)
          DO K=1, I-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
        ENDDO
C
        DO I=J, N
          SUM = A(I,J)
          DO K=1, J-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
        ENDDO
C
        AINV = 1.0/A(J,J)

        DO I=J+1, N
          A(I,J) = A(I,J)*AINV
        ENDDO
C
      ENDDO
C
      RETURN
      END ! ZLUDCMP



      SUBROUTINE ZBAKSUB(NSIZ,N,A,B)
      IMPLICIT NONE

      INTEGER NSIZ,N
      COMPLEX*16 A(NSIZ,NSIZ), B(NSIZ)

      INTEGER I,J,II
      COMPLEX*16 SUM
C
      II = 0
      DO I=1, N
        SUM = B(I)
        IF(II.NE.0) THEN
         DO J=II, I-1
           SUM = SUM - A(I,J)*B(J)
         ENDDO
        ELSE IF(SUM.NE.0.0) THEN
         II = I
        ENDIF
        B(I) = SUM
      ENDDO
C
      DO I=N, 1, -1
        SUM = B(I)
        IF(I.LT.N) THEN
         DO J=I+1, N
           SUM = SUM - A(I,J)*B(J)
         ENDDO
        ENDIF
        B(I) = SUM/A(I,I)
      ENDDO
C
      RETURN
      END ! ZBAKSUB



      SUBROUTINE ZLUDCMPI(NSIZ,N,A,INDX)
C     *******************************************************
C     *                                                     *
C     *   Factors a full NxN matrix into an LU form.        *
C     *   Subr. BAKSUB can back-substitute it with some RHS.*
C     *   Assumes matrix is non-singular...                 *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   A is the matrix...                                *
C     *     ...replaced with its LU factors.                *
C     *                                                     *
C     *                              Mark Drela  1988       *
C     *******************************************************
      IMPLICIT NONE

      INTEGER NSIZ,N
      COMPLEX*16 A(NSIZ,NSIZ)
      INTEGER INDX(NSIZ)
C
      INTEGER NVX
      PARAMETER (NVX=1000)
      COMPLEX*16 VV(NVX)
C
      INTEGER I,J,K,IMAX
      COMPLEX*16 CAMAX,SUM,ATMP

      REAL RAMAX, RDUM
C
      IF(N.GT.NVX) STOP 'ZLUDCMPI: Array overflow. Increase NVX.'
C
      DO I=1, N
        CAMAX = (0.0,0.0)
        DO J=1, N
          IF(ABS(A(I,J)) .GT. ABS(CAMAX)) CAMAX = A(I,J)
        ENDDO
        VV(I) = 1.0/CAMAX
      ENDDO
C
      DO J=1, N
        DO I=1, J-1
          SUM = A(I,J)
          DO K=1, I-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
        ENDDO
C
        RAMAX = 0.
        DO I=J, N
          SUM = A(I,J)
          DO K=1, J-1
            SUM = SUM - A(I,K)*A(K,J)
          ENDDO
          A(I,J) = SUM
C
          RDUM = ABS(VV(I)*SUM)
          IF(RDUM.GE.RAMAX) THEN
           IMAX = I
           RAMAX = RDUM
          ENDIF
        ENDDO
C
        IF(J.NE.IMAX) THEN
         DO K=1, N
           ATMP = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = ATMP
         ENDDO
         VV(IMAX) = VV(J)
        ENDIF
C
        INDX(J) = IMAX
        IF(J.NE.N) THEN

cc        write(*,1200) j, indx(j), a(j,j)
cc 1200   format(1x,2i4,f18.8)

         ATMP = 1.0/A(J,J)
         DO I=J+1, N
           A(I,J) = A(I,J)*ATMP
         ENDDO
        ENDIF
C
      ENDDO
C
      RETURN
      END ! ZLUDCMPI



      SUBROUTINE ZBAKSUBI(NSIZ,N,A,INDX,B)
      IMPLICIT NONE

      INTEGER NSIZ,N
      COMPLEX*16 A(NSIZ,NSIZ), B(NSIZ)
      INTEGER INDX(NSIZ)
C
      INTEGER I,J,II,LL
      COMPLEX*16 SUM

      II = 0
      DO I=1, N
        LL = INDX(I)
        SUM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
         DO J=II, I-1
           SUM = SUM - A(I,J)*B(J)
         ENDDO
        ELSE IF( SUM .NE. (0.0,0.0) ) THEN
         II = I
        ENDIF
        B(I) = SUM
      ENDDO
C
      DO I=N, 1, -1
        SUM = B(I)
        IF(I.LT.N) THEN
         DO J=I+1, N
           SUM = SUM - A(I,J)*B(J)
         ENDDO
        ENDIF
        B(I) = SUM/A(I,I)
      ENDDO
C
      RETURN
      END ! ZBAKSUBI


