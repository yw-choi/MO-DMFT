        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 19 02:02:59 2016
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BRENT__genmod
          INTERFACE 
            FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN,NS,NWLOC,NXSIZE,OMEGA,   &
     &D_EV,NBATH,NORB,PCOM,XICOM,COMM)
              INTEGER(KIND=4) :: NORB
              INTEGER(KIND=4) :: NWLOC
              REAL(KIND=8) :: AX
              REAL(KIND=8) :: BX
              REAL(KIND=8) :: CX
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: XMIN
              INTEGER(KIND=4) :: NS
              INTEGER(KIND=4) :: NXSIZE
              REAL(KIND=8) :: OMEGA(NWLOC)
              COMPLEX(KIND=8) :: D_EV(NORB,NWLOC)
              INTEGER(KIND=4) :: NBATH
              REAL(KIND=8) :: PCOM
              REAL(KIND=8) :: XICOM
              INTEGER(KIND=4) :: COMM
              REAL(KIND=8) :: BRENT
            END FUNCTION BRENT
          END INTERFACE 
        END MODULE BRENT__genmod
