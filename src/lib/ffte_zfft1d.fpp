!C
!C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!C
!C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
!C                BY
!C         DAISUKE TAKAHASHI
!C         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!C         UNIVERSITY OF TSUKUBA
!C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!C         E-MAIL: daisuke@cs.tsukuba.ac.jp
!C
!C
!C     1-D COMPLEX FFT ROUTINE
!C
!C     FORTRAN77 SOURCE PROGRAM
!C
!C     CALL ZFFT1D(A,N,IOPT,B)
!C
!C     A(N) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX(WP))
!C     B(N*2) IS WORK/COEFFICIENT VECTOR (COMPLEX(WP))
!C     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*4)
!C       -----------------------------------
!C         N = (2**IP) * (3**IQ) * (5**IR)
!C       -----------------------------------
!C     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
!C          = -1 FOR FORWARD TRANSFORM
!C          = +1 FOR INVERSE TRANSFORM
!C
!C     WRITTEN BY DAISUKE TAKAHASHI
!C

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

      SUBROUTINE ZFFT1D(A,N,IOPT,B)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
      COMPLEX(WP) A(*),B(*)
      DIMENSION IP(3)
!C
      CALL FACTOR(N,IP)
!C
      IF (IOPT .EQ. 1) THEN
        DO 10 I=1,N
          A(I)=CONJG(A(I))
   10   CONTINUE
      END IF
!CJJS
      IF (IOPT .EQ. 2) THEN
        DO I=1,N/2+1
          A(I)=CONJG(A(I))
        ENDDO
      ENDIF
!CJJE
!C
        IF (IOPT .EQ. 0) THEN
          CALL SETTBL(B(N+1),N)
          RETURN
        END IF
!C
!CJJS
      IF (IOPT .EQ. 2) THEN
        A(1)=DCMPLX(REAL(A(1),WP),0.0_WP)
        DO I=2,N/2+1
          A(N-I+2)=CONJG(A(I))
        ENDDO
      ENDIF
        CALL FFT235(A,B,B(N+1),N,IP)
!C
      IF (IOPT .EQ. 1 .OR. IOPT .EQ. 2) THEN
        DN=1.0_WP
        DO 20 I=1,N
          A(I)=CONJG(A(I))*DN
   20   CONTINUE
      END IF
      RETURN
      END
