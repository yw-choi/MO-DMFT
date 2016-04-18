        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 19 02:02:59 2016
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE F1DIM__genmod
          INTERFACE 
            FUNCTION F1DIM(X,NS,NWLOC,NXSIZE,OMEGA,D_EV,NBATH,NORB,PCOM,&
     &XICOM,COMM)
              INTEGER(KIND=4) :: NORB
              INTEGER(KIND=4) :: NXSIZE
              INTEGER(KIND=4) :: NWLOC
              REAL(KIND=8) :: X
              INTEGER(KIND=4) :: NS
              REAL(KIND=8) :: OMEGA(NWLOC)
              COMPLEX(KIND=8) :: D_EV(NORB,NWLOC)
              INTEGER(KIND=4) :: NBATH
              REAL(KIND=8) :: PCOM(NXSIZE)
              REAL(KIND=8) :: XICOM(NXSIZE)
              INTEGER(KIND=4) :: COMM
              REAL(KIND=8) :: F1DIM
            END FUNCTION F1DIM
          END INTERFACE 
        END MODULE F1DIM__genmod
