        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 10 10:59:15 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STAR5__genmod
          INTERFACE 
            SUBROUTINE STAR5(A,B,C,R,ND,N,M,KORT)
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ND
              COMPLEX(KIND=8) :: A(ND,N)
              COMPLEX(KIND=8) :: B(ND,M)
              COMPLEX(KIND=8) :: C(ND,N)
              COMPLEX(KIND=8) :: R(ND,M+3)
              INTEGER(KIND=4) :: KORT
            END SUBROUTINE STAR5
          END INTERFACE 
        END MODULE STAR5__genmod
