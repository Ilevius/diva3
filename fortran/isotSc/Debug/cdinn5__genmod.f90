        !COMPILER-GENERATED INTERFACE MODULE: Sun Nov 20 12:34:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CDINN5__genmod
          INTERFACE 
            SUBROUTINE CDINN5(CF,A,B,EPS,PR,N,SB,RD)
              INTEGER(KIND=4) :: N
              EXTERNAL CF
              COMPLEX(KIND=8) :: A
              COMPLEX(KIND=8) :: B
              REAL(KIND=8) :: EPS
              REAL(KIND=8) :: PR
              COMPLEX(KIND=8) :: SB(N)
              COMPLEX(KIND=8) :: RD(N)
            END SUBROUTINE CDINN5
          END INTERFACE 
        END MODULE CDINN5__genmod
