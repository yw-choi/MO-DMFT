        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 19 02:02:59 2016
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LINMIN__genmod
          INTERFACE 
            SUBROUTINE LINMIN(P,XI,N,FRET,NS,NWLOC,NXSIZE,OMEGA,D_EV,   &
     &NBATH,NORB,COMM)
              INTEGER(KIND=4) :: NORB
              INTEGER(KIND=4) :: NWLOC
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: P(N)
              REAL(KIND=8) :: XI(N)
              REAL(KIND=8) :: FRET
              INTEGER(KIND=4) :: NS
              INTEGER(KIND=4) :: NXSIZE
              REAL(KIND=8) :: OMEGA(NWLOC)
              COMPLEX(KIND=8) :: D_EV(NORB,NWLOC)
              INTEGER(KIND=4) :: NBATH
              INTEGER(KIND=4) :: COMM
            END SUBROUTINE LINMIN
          END INTERFACE 
        END MODULE LINMIN__genmod
