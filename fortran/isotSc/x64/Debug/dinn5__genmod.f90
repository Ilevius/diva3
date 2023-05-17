        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 10 10:59:16 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DINN5__genmod
          INTERFACE 
            SUBROUTINE DINN5(CF,T1,T2,T3,T4,TM,TP,EPS,PR,GR,N,RD)
              INTEGER(KIND=4) :: N
              EXTERNAL CF
              REAL(KIND=8) :: T1
              REAL(KIND=8) :: T2
              REAL(KIND=8) :: T3
              REAL(KIND=8) :: T4
              REAL(KIND=8) :: TM
              REAL(KIND=8) :: TP
              REAL(KIND=8) :: EPS
              REAL(KIND=8) :: PR
              REAL(KIND=8) :: GR
              COMPLEX(KIND=8) :: RD(N)
            END SUBROUTINE DINN5
          END INTERFACE 
        END MODULE DINN5__genmod
