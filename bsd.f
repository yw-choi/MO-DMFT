      SUBROUTINE CPUTIM (TIME)

      DOUBLE PRECISION TIME
      REAL TIMES(2)
C
      TIME = ETIME(TIMES)
      END
      subroutine flush_(i)
      call flush(i)
      end
