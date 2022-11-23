
!STAR5 is a solver for a complex linear algebraic system Ax=b 
!       by the double-orthogonalization method 
! N - system size; M - number of the right-hand parts;
! Nd - dimension (number of the lines) declared in fact; 
! A(Nd,N) - N*N complex matrix of the system; 
! B(Nd,M) - N*M complex matrix of the right-hand parts;
! C(Nd,N),R(Nd,M+3) - complex work arrays
! KORT - number of the ortogonalization cycles 
!(2 - 3 are usually quite enough)
 
      SUBROUTINE STAR5(A,B,C,R,ND,N,M,kort)
      implicit none
      
	  integer Nd,N,M,kort,i,k,k1,kw,jr,im,Ms,Msm,Mt
	  
	  real*8 ra,rt

	  COMPLEX*16 A(Nd,N),B(Nd,M),C(Nd,N),R(Nd,M+3),TK,ST

	  Ms=M+1; Msm=M+2; Mt=M+3

! orthogonolization

      DO 1 K=1,N
        K1=K-1;  RA=0D0;  KW=0

        DO 2 I=1,N
          ST=A(I,K); R(I,Ms)=ST;  R(I,Msm)=ST
  2     RA=RA+ST*CONJG(ST)
        RA=SQRT(RA)

        IF(K.NE.1)GOTO 12
        R(1,Mt)=RA

        DO 13 I=1,N
 13     C(I,1)=R(I,Ms)/RA
        GOTO 1

 12     DO 4 JR=1,K1
          ST=0D0
          DO 5 I=1,N
  5       ST=ST+R(I,Ms)*CONJG(C(I,JR))

          DO 7 I=1,N
  7       R(I,Msm)=R(I,Msm)-ST*C(I,JR)

  4     CONTINUE

        RT=0D0
        DO 14 I=1,N
 14     RT=RT+ABS(R(I,Msm))**2
        
        RT=SQRT(RT)
        IF(RT.GT.1D-13)GOTO 9
     
          DO 8 I=1,N
  8       R(I,Ms)=0D0

          ST=0D0
        GOTO 28
 
  9     KW=KW+1
        DO 10 I=1,N
          ST=R(I,Msm)/RT;  R(I,Ms)=ST
 10     R(I,Msm)=ST

        IF(KW.LT.KORT)GOTO 12
 
 34     ST=0D0
        DO 21 I=1,N
 21     ST=ST+R(I,Ms)*CONJG(A(I,K))
  
 28     R(K,Mt)=ST
 
        DO 23 I=1,N
 23     C(I,K)=R(I,Ms)

  1   CONTINUE
      
! backward steps
	
	K=N
      DO 44 IM=1,M
  3     DO 44 I=1,N
 44     R(I,IM)=B(I,IM)

 45     DO 49 I=1,N
 49     R(I,Ms)=C(I,K)

        TK=R(K,Mt)
        IF(ABS(TK).GT.1D-12)GOTO 50
          DO 55 IM=1,M
 55       B(K,IM)=0D0
        GOTO 56
 
 50     DO 51 IM=1,M
          ST=0D0
          DO 52 I=1,N
 52       ST=ST+R(I,Ms)*CONJG(R(I,IM))

 51     B(K,IM)=CONJG(ST/TK)
        
		IF(K.EQ.1) RETURN
      
	    DO 54 I=1,N
          TK=A(I,K)
          DO 54 IM=1,M
 54     R(I,IM)=R(I,IM)-B(K,IM)*TK

 56     K=K-1
        IF(K.GE.1)GOTO 45
      
	    RETURN
        END

!-----------------------

      SUBROUTINE DGAUS(A,ND,N,CD)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 A(ND,ND),CD,CR
     
      CD=(1D0,0d0)
      DO 1 I=1,N
        J=I
        JS=J
        IF(I.EQ.N)GOTO 2
        I1=I+1
        AM=ABS(A(I,J))

        DO 3 L=I1,N
          AT=ABS(A(I,L))
          IF(AT.LE.AM)GOTO 3
          JS=L
          AM=AT
  3     CONTINUE
      
        IF(AM.GT.1D-20)GOTO 6
        CD=(0D0,0d0)
      RETURN

  6     IF(JS.EQ.J)GOTO 4
        DO 5 K=1,N
          CR=A(K,J)
          A(K,J)=A(K,JS)
  5     A(K,JS)=CR
        CD=-CD
  
  4     DO 7 K=I1,N
          CR=A(K,J)/A(I,J)
          DO 7 L=I1,N
  7     A(K,L)=A(K,L)-CR*A(I,L)
  
  2     CD=CD*A(I,J)
      
  1   CONTINUE                    
  
      RETURN
      END

!-------------------------

! DET - determinant of a complex matrix A(N,N)

      SUBROUTINE DSTAR(A,C,S,SM,DET,ND,N,KORT)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 A(ND,ND),DET,S(ND),SM(ND),ST,C(ND,ND)
      
      DO 1 K=1,N
        K1=K-1
        KW=0
        IF(K.NE.1)GO TO 15
        RA=0D0
        DO 2 I=1,N
          ST=A(I,K)
          S(I)=ST
          SM(I)=ST
  2     RA=RA+ST*CONJG(ST)
        RA=SQRT(RA)
        IF(RA.GT.1D-15)GO TO 50
        DET=0D0
      RETURN

50      DO 13 I=1,N
 13     C(I,1)=S(I)/RA
        GOTO 1
        
 15     DO 16 I=1,N
          ST=A(I,K)
          S(I)=ST
 16     SM(I)=ST
 
 12     DO 4 JR=1,K1
          ST=0D0
          DO 5 I=1,N
  5       ST=ST+S(I)*CONJG(C(I,JR))
  
          DO 7 I=1,N
  7       SM(I)=SM(I)-ST*C(I,JR)
  4     CONTINUE
      
        RT=0D0
        DO 14 I=1,N
 14     RT=RT+ABS(SM(I))**2
        RT=SQRT(RT)
        RA=RA*RT                        
              
        IF(RT.GT.1D-15)GOTO 9
        DET=0D0
      RETURN
  
  9     KW=KW+1
        DO 10 I=1,N
          ST=SM(I)/RT
          S(I)=ST
 10     SM(I)=ST
        IF(KW.LT.KORT)GOTO 12

        DO 23 I=1,N
 23     C(I,K)=S(I)
  1   CONTINUE
     
      CALL DGAUS(C,ND,N,DET)
     
      DET=DET*RA
      RETURN
      END
