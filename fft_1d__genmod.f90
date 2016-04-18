        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 19 02:01:55 2016
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FFT_1D__genmod
          INTERFACE 
            SUBROUTINE FFT_1D(N,IN_FFT,OUT_FFT)
              INTEGER(KIND=8) :: N
              REAL(KIND=8) :: IN_FFT(N)
              COMPLEX(KIND=8) :: OUT_FFT(N/2+1)
            END SUBROUTINE FFT_1D
          END INTERFACE 
        END MODULE FFT_1D__genmod
