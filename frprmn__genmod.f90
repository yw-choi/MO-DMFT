        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 19 02:02:59 2016
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FRPRMN__genmod
          INTERFACE 
            SUBROUTINE FRPRMN(P,N,FTOL,ITER,NXSIZE,NWLOC,NS,OMEGA,D_EV, &
     &NBATH,NORB,FRET,COMM)
              INTEGER(KIND=4) :: NORB
              INTEGER(KIND=4) :: NWLOC
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: P(N)
              REAL(KIND=8) :: FTOL
              INTEGER(KIND=4) :: ITER
              INTEGER(KIND=4) :: NXSIZE
              INTEGER(KIND=4) :: NS
              REAL(KIND=8) :: OMEGA(NWLOC)
              COMPLEX(KIND=8) :: D_EV(NORB,NWLOC)
              INTEGER(KIND=4) :: NBATH
              REAL(KIND=8) :: FRET
              INTEGER(KIND=4) :: COMM
            END SUBROUTINE FRPRMN
          END INTERFACE 
        END MODULE FRPRMN__genmod
