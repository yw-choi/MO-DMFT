        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 19 02:02:59 2016
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MNBRAK__genmod
          INTERFACE 
            SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC,NS,NWLOC,NXSIZE,   &
     &OMEGA,D_EV,NBATH,NORB,PCOM,XICOM,COMM)
              INTEGER(KIND=4) :: NORB
              INTEGER(KIND=4) :: NWLOC
              REAL(KIND=8) :: AX
              REAL(KIND=8) :: BX
              REAL(KIND=8) :: CX
              REAL(KIND=8) :: FA
              REAL(KIND=8) :: FB
              REAL(KIND=8) :: FC
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              INTEGER(KIND=4) :: NS
              INTEGER(KIND=4) :: NXSIZE
              REAL(KIND=8) :: OMEGA(NWLOC)
              COMPLEX(KIND=8) :: D_EV(NORB,NWLOC)
              INTEGER(KIND=4) :: NBATH
              REAL(KIND=8) :: PCOM
              REAL(KIND=8) :: XICOM
              INTEGER(KIND=4) :: COMM
            END SUBROUTINE MNBRAK
          END INTERFACE 
        END MODULE MNBRAK__genmod
