MODULE functions
    real*8 pi
    complex*16 ci
    parameter (pi=3.141592653589793d0)
    parameter (ci = (0d0,1d0))
    
    CONTAINS
    
        FUNCTION matrixA(alfa, kappa, kappaCap, lambda, mu)
        real*8 kappa(2), kappaCap(2), lambda(2), mu(2)
        complex*16 matrixA(4,4), sigma(2), sigmaCap(2), alfa
        !alfa=1d0;
            sigma = makeSigma(kappa, alfa)
            sigmaCap = makeSigma(kappaCap, alfa)
            ! elements defining
            matrixA(1,1) = -ci*alfa; 
            matrixA(1,2) = sigma(2); 
            matrixA(1,3) = ci*alfa; 
            matrixA(1,4) = sigmaCap(2);
            
            matrixA(2,1) = -sigma(1); 
            matrixA(2,2) = -ci*alfa; 
            matrixA(2,3) = -sigmaCap(1); 
            matrixA(2,4) = ci*alfa;
            
            matrixA(3,1) = 2d0*mu(1)*ci*alfa*sigma(1); 
            matrixA(3,2) = -mu(1)*(sigma(2)**2+alfa**2); 
            matrixA(3,3) = 2d0*mu(2)*ci*alfa*sigmaCap(1); 
            matrixA(3,4) = mu(2)*(sigmaCap(2)**2+alfa**2);
            
            matrixA(4,1) = mu(1)*(alfa**2+sigma(2)**2); 
            matrixA(4,2) = 2d0*mu(1)*ci*alfa*sigma(2); 
            matrixA(4,3) = -mu(2)*(alfa**2+sigmaCap(2)**2); 
            matrixA(4,4) = 2d0*mu(2)*ci*alfa*sigmaCap(2);                     
        END FUNCTION matrixA
        
        
        

        
        
        FUNCTION delta(alfa, mu, kappa)
        ! single mu!!!
        real*8 mu, kappa(2)
        complex*16 alfa, delta, sigma(2)
            sigma = makeSigma(kappa, alfa)
            delta = 2d0*mu*(-(alfa**2 - 0.5d0*kappa(2)**2)**2 + alfa**2*sigma(1)*sigma(2)) 
        END
        
        
        
        FUNCTION makeSigma(kappa, alfa)
        real*8 kappa(2)
        complex*16 ci, alfa, makeSigma(2)
        parameter (ci = (0d0,1d0))
            if (Imag(alfa) == 0d0)  then 
                    if (abs(alfa) < Kappa(1)) then
                        MakeSigma(1) = -ci*sqrt(Kappa(1)**2 - real(alfa)**2)
                    else 
                        MakeSigma(1) = sqrt(real(alfa)**2 - Kappa(1)**2)
                    endif    
                    
                    if (abs(alfa) < Kappa(2)) then
                        MakeSigma(2) = -ci*sqrt(Kappa(2)**2 - real(alfa)**2)
                    else 
                        MakeSigma(2) = sqrt(real(alfa)**2 - Kappa(2)**2)
                    endif  
                    
                
            else 
                    MakeSigma(1) = sqrt(alfa**2 - cmplx(Kappa(1)**2))
                    MakeSigma(2) = sqrt(alfa**2 - cmplx(Kappa(2)**2))
            endif  
        END
    
        FUNCTION U1(alfa, kappa, mu)
        ! single mu!!!
            IMPLICIT NONE;
            real*8 kappa(2), mu
            complex*16 alfa, U1, Sigma(2)
                Sigma = makeSigma(kappa, alfa)
                U1 = ci*alfa*(alfa**2 - 0.5d0*(kappa(2))**2)/Delta(alfa, mu, kappa)    
        END FUNCTION U1
    
        FUNCTION U2(alfa, kappa, mu)
        ! single mu!!!
            IMPLICIT NONE;
            real*8 kappa(2), mu
            complex*16 alfa, U2, Sigma(2)
                Sigma = makeSigma(kappa, alfa)
                U2 = -ci*alfa*Sigma(1)*Sigma(2)/Delta(alfa, mu, kappa)    
        END FUNCTION U2
    
        FUNCTION W1(alfa, kappa, mu)
        ! single mu!!!
            IMPLICIT NONE;
            complex*16 alfa, W1, Sigma(2)
            real*8 kappa(2), mu
                Sigma = makeSigma(kappa, alfa)
                W1 = -Sigma(1)*(alfa**2 - 0.5d0*Kappa(2)**2)/Delta(alfa, mu, kappa)
        END FUNCTION W1
        
        FUNCTION W2(alfa, kappa, mu)
        ! single mu!!!
            IMPLICIT NONE;
            complex*16 alfa, W2, Sigma(2)
            real*8 kappa(2), mu
                Sigma = makeSigma(kappa, alfa)
                W2 = Sigma(1)*alfa**2/Delta(alfa, mu, kappa)
        END FUNCTION W2
        
        FUNCTION partB(i, mu, lambda, kappa, alfa)
        ! single mu!!!
            IMPLICIT NONE;
            integer i
            real*8 kappa(2), mu, lambda
            complex*16 alfa, partB(4,1), Sigma(2), theU1, theU2, theW1, theW2         
                sigma = makeSigma(kappa, alfa)   
                if (i == 1) then
                    theU1 = U1(alfa, kappa, mu); theW1 = W1(alfa, kappa, mu);
                    partB(1,1) = -theU1
                    partB(2,1) = -theW1
                    partB(3,1) = -mu*(sigma(1)*theU1-ci*alfa*theW1)
                    partB(4,1) = lambda*ci*alfa*theU1-(lambda+2d0*mu)*sigma(1)*theW1
                else 
                    theU2 = U2(alfa, kappa, mu); theW2 = W2(alfa, kappa, mu);
                    partB(1,1) = -theU2
                    partB(2,1) = -theW2
                    partB(3,1) = -mu*(sigma(2)*theU2-ci*alfa*theW2)
                    partB(4,1) = lambda*ci*alfa*theU2-(lambda+2d0*mu)*sigma(2)*theW2
                endif
        END FUNCTION partB  
        
        
        FUNCTION CramDelta(i, j, kappa, kappaCap, lambda, mu, alfa)
        implicit none
        integer i, j
        real*8 kappa(2), kappaCap(2), lambda(2), mu(2)
        complex*16 CramDelta, alfa, matrix(4,4), column(4,1)
        complex*16 C(4,4), S(4), SM(4)
     
            matrix = matrixA(alfa, kappa, kappaCap, lambda, mu)         
            if (i>0 .AND. j>0) then 
                column = partB(j, mu(1), lambda(1), kappa, alfa)
                matrix(:,i) = column(:,1)
            endif          
            call DSTAR(matrix,C,S,SM,CramDelta,4,4,2)      
        END FUNCTION CramDelta
        
    
    
END MODULE functions