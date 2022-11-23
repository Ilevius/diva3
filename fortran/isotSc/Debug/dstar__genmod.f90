        !COMPILER-GENERATED INTERFACE MODULE: Sun Nov 20 12:34:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSTAR__genmod
          INTERFACE 
            SUBROUTINE DSTAR(A,C,S,SM,DET,ND,N,KORT)
              INTEGER(KIND=4) :: ND
              COMPLEX(KIND=8) :: A(ND,ND)
              COMPLEX(KIND=8) :: C(ND,ND)
              COMPLEX(KIND=8) :: S(ND)
              COMPLEX(KIND=8) :: SM(ND)
              COMPLEX(KIND=8) :: DET
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: KORT
            END SUBROUTINE DSTAR
          END INTERFACE 
        END MODULE DSTAR__genmod
