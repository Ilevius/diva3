MODULE test
    CONTAINS
    
        SUBROUTINE source_field(testName, kappa, mu, lambda)
        use functions
        IMPLICIT NONE
        character (len=*) :: testName
        real*8 alfa, kappa(2), mu, lambda, theH
        complex*16 alfa_c, sigma(2), object1, object2
        
!                                           Here we test if mu(U'-i*alfa*W) equals Q1, as Q1=0 it's the same to U'=i*alfa*W
            if (testName == 'first boundary') then
                open(1, file='test.txt', FORM='FORMATTED')
                write(1, *) "%Here we test if U'=i*alfa*W at z = 0 "
                write(1, *) "% Objects 1, 2 should coinside "
                write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
                do alfa=-100, 100d0, 1d0
                    alfa_c = cmplx(alfa)
                    sigma = makeSigma(kappa, alfa_c)
                    object1 = sigma(1)*U1(alfa_c, kappa, mu)+sigma(2)*U2(alfa_c, kappa, mu)
                    object2 = ci*alfa*(W1(alfa_c, kappa, mu)+W2(alfa_c, kappa, mu))
                    write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
                enddo 
                close(1)
                
!                                           Here we test if -lambda*i*alfa*U + (lambda+2*mu)*W' equals Q2 or the same  (lambda+2*mu)*W' = lambda*i*alfa*U+1  

            else if (testName == 'second boundary') then
                open(1, file='test.txt', FORM='FORMATTED')
                write(1, *) "%Here we test if (lambda+2*mu)*W' = lambda*i*alfa*U+1 at z = 0 "
                write(1, *) "% Objects 1, 2 should coinside "
                write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
                do alfa=-100, 100d0, 1d0
                    alfa_c = cmplx(alfa)
                    sigma = makeSigma(kappa, alfa_c)
                    object1 = (lambda+2d0*mu)*(sigma(1)*W1(alfa_c, kappa, mu)+sigma(2)*W2(alfa_c, kappa, mu))
                    object2 = lambda*ci*alfa*(U1(alfa_c, kappa, mu)+U2(alfa_c, kappa, mu)) + 1d0
                    write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
                enddo 
                close(1)
                
                
!                       Here we test if (kappa2^2*mu-alfa^2(lam+2mu))U = i*alfa*(lam+mu)*W'-mu*U''
            else if (testName == 'first Lame') then
                open(1, file='test.txt', FORM='FORMATTED')
                write(1, *) "%Here we test if (lambda+2*mu)*W' = lambda*i*alfa*U+1 at z = 0 "
                write(1, *) "% Objects 1, 2 should coinside "
                write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
                do alfa=-100, 100d0, 1d0
                    alfa_c = cmplx(alfa)
                    theH = -0.9d0
                    sigma = makeSigma(kappa, alfa_c)
                    object1 = (kappa(2)**2*mu-alfa**2*(lambda+2d0*mu))*(U1(alfa_c, kappa, mu)*exp(sigma(1)*theH)+U2(alfa_c, kappa, mu)*exp(sigma(2)*theH))
                    object2 = ci*alfa*(lambda+mu)*(W1(alfa_c, kappa, mu)*sigma(1)*exp(sigma(1)*theH)+W2(alfa_c, kappa, mu)*sigma(2)*exp(sigma(2)*theH))
                    object2 = object2 - mu*(U1(alfa_c, kappa, mu)*sigma(1)**2*exp(sigma(1)*theH)+U2(alfa_c, kappa, mu)*sigma(2)**2*exp(sigma(2)*theH)) 
                    write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
                enddo 
                close(1)
  
                
!                       Here we test if (kappa2^2*mu-alfa^2*mu)W = i*alfa*(lam+mu)*U'-(lambda+2mu)*W''                
            else if (testName == 'second Lame') then
                open(1, file='test.txt', FORM='FORMATTED')
                write(1, *) "%Here we test if (lambda+2*mu)*W' = lambda*i*alfa*U+1 at z = 0 "
                write(1, *) "% Objects 1, 2 should coinside "
                write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
                do alfa=-100, 100d0, 1d0
                    alfa_c = cmplx(alfa)
                    theH = -0.3d0
                    sigma = makeSigma(kappa, alfa_c)
                    object1 = (kappa(2)**2*mu-alfa**2*mu)*(W1(alfa_c, kappa, mu)*exp(sigma(1)*theH)+W2(alfa_c, kappa, mu)*exp(sigma(2)*theH))
                    object2 = ci*alfa*(lambda+mu)*(U1(alfa_c, kappa, mu)*sigma(1)*exp(sigma(1)*theH)+U2(alfa_c, kappa, mu)*sigma(2)*exp(sigma(2)*theH))
                    object2 = object2 - (lambda+2d0*mu)*(W1(alfa_c, kappa, mu)*sigma(1)**2*exp(sigma(1)*theH)+W2(alfa_c, kappa, mu)*sigma(2)**2*exp(sigma(2)*theH)) 
                    write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
                enddo 
                close(1)
                
            endif    
        
            !   осталось протестировать: 
            !           уравнения для образов и праобразов
            !           граничные условия для праобразов
            
        END SUBROUTINE source_field
        
        SUBROUTINE xzPointsTest(x, z, psi, R, h, n)
            integer n
            real*8 x(n), z(n), psi(n), R(n), h
            open(1, file='points.txt', FORM='FORMATTED')
            write(1,*) '% x, z, xRe, zRe, psi, R, h'
            do i = 1, n          
                write(1,'(7E15.6E3)') x(i), z(i), R(i)*cosd(psi(i)), R(i)*sind(psi(i)) - 2d0*h, psi(i), R(i), h
            enddo
            close(1)
        END SUBROUTINE xzPointsTest
    
        
        
        !SUBROUTINE reflected_field(testName, kappa, kappaCap, mu, lambda, h)
        !use functions
        !IMPLICIT NONE
        !character (len=*) :: testName
        !real*8 alfa, kappa(2), kappaCap(2), mu(2), lambda(2), h
        !complex*16 alfa_c, sigma(2), sigmaCap(2), object1, object2, t(4)
        !complex*16 output
        !
        !    !                                           Here we test if U0+Uminus equals Uplus at z = -h
        !    if (testName == 'U') then
        !        open(1, file='test.txt', FORM='FORMATTED')
        !        write(1, *) "%Here we test if U0+Uminus equals Uplus at z = -h "
        !        write(1, *) "% Objects 1, 2 should coinside "
        !        write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
        !        do alfa=-100, 100d0, 0.5d0
        !            alfa_c = cmplx(alfa)
        !            sigma = makeSigma(kappa, alfa_c)
        !            sigmaCap = makeSigma(kappaCap, alfa_c)
        !            t(1) = CramDelta(1, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
        !            t(1) = t(1) + CramDelta(1, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
        !            
        !            t(2) = CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
        !            t(2) = t(2) + CramDelta(2, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
        !            
        !            t(3) = CramDelta(3, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
        !            t(3) = t(3) + CramDelta(3, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
        !            
        !            t(4) = CramDelta(4, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
        !            t(4) = t(4) + CramDelta(4, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
        !
        !            ! U0 + Uminus
        !            object1 = -ci*alfa*t(1)+sigma(2)*t(2) + U1(alfa_c, kappa, mu(1))*exp(-sigma(1)*h) + U2(alfa_c, kappa, mu(1))*exp(-sigma(2)*h)
        !            ! Uplus
        !            object2 = -ci*alfa*t(3)-sigmaCap(2)*t(4)
        !
        !
        !            write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
        !        enddo 
        !        close(1)
        !        
        !    !                                           Here we test if W0+Wminus equals Wplus at z = -h   
        !
        !    else if (testName == 'first Lame') then
        !        open(1, file='test.txt', FORM='FORMATTED')
        !        write(1, *) "%Here we test if W0+Wminus equals Wplus at z = -h "
        !        write(1, *) "% Objects 1, 2 should coinside "
        !        write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
        !        do alfa=-100, 100d0, 0.5d0
        !            alfa_c = cmplx(alfa)
        !            sigma = makeSigma(kappa, alfa_c)
        !            sigmaCap = makeSigma(kappaCap, alfa_c)
        !
        !            ! U0 + Uminus
        !            call ups_integrand_horis(alfa, output, 1)
        !            object1 = output
        !            ! Uplus
        !            object2 = sigmaCap(1)*t(3)-ci*alfa*t(4)
        !
        !
        !            write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
        !        enddo 
        !        close(1)
        !        
        !        
        !    else if (testName == 'W') then
        !        open(1, file='test.txt', FORM='FORMATTED')
        !        write(1, *) "%Here we test if W0+Wminus equals Wplus at z = -h "
        !        write(1, *) "% Objects 1, 2 should coinside "
        !        write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
        !        do alfa=-100, 100d0, 0.5d0
        !            alfa_c = cmplx(alfa)
        !            sigma = makeSigma(kappa, alfa_c)
        !            sigmaCap = makeSigma(kappaCap, alfa_c)
        !            t(1) = CramDelta(1, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
        !            t(1) = t(1) + CramDelta(1, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
        !            
        !            t(2) = CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
        !            t(2) = t(2) + CramDelta(2, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
        !            
        !            t(3) = CramDelta(3, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
        !            t(3) = t(3) + CramDelta(3, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
        !            
        !            t(4) = CramDelta(4, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
        !            t(4) = t(4) + CramDelta(4, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
        !            ! U0 + Uminus
        !            object1 = -sigma(1)*t(1)-ci*alfa*t(2) + W1(alfa_c, kappa, mu(1))*exp(-sigma(1)*h) + W2(alfa_c, kappa, mu(1))*exp(-sigma(2)*h)
        !            ! Uplus
        !            object2 = sigmaCap(1)*t(3)-ci*alfa*t(4)
        !            write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
        !        enddo 
        !        close(1)
        !
        !        
        !        
        !
        !        
        !    endif  
        !    
        !    !   осталось протестировать: 
        !    !           уравнения для образов и праобразов
        !    !           граничное условие для праобразов
        !
        !END SUBROUTINE reflected_field
        !
        

END MODULE test