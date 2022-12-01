MODULE test
    CONTAINS
    
        SUBROUTINE source_field(testName, kappa, mu)
        use functions
        IMPLICIT NONE
        character (len=*) :: testName
        real*8 alfa, kappa(2), mu
        complex*16 alfa_c, sigma(2), object1, object2
            ! Here we test if mu(U'-i*alfa*W) equals Q1
            if (testName == 'first boundary') then
                open(1, file='test.txt', FORM='FORMATTED')
                write(1, *) "%Here we test if mu(U'-i*alfa*W) equals Q1 at z = 0 "
                write(1, *) "% Objects 1, 2 should coinside "
                write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
                do alfa=-1000, 1000d0, 2d0
                    alfa_c = cmplx(alfa)
                    sigma = makeSigma(kappa, alfa_c)
                    object1 = sigma(1)*U1(alfa_c, kappa, mu)+sigma(2)*U2(alfa_c, kappa, mu)
                    object2 = ci*alfa*(W1(alfa_c, kappa, mu)+W2(alfa_c, kappa, mu))
                    write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
                enddo 
                close(1)
            ! Here we test if -lambda*i*alfa*U + (lambda+2*mu)*W' equals Q2   
                                                                                    ! T O   B E   D O N E  !!!!!!!!!!!!
            else if (testName == 'second boundary') then

                
            endif    
        
        END SUBROUTINE source_field
        
        
        SUBROUTINE scattered_field(testName, kappa, kappaCap, mu, lambda, h)
        use functions
        IMPLICIT NONE
        character (len=*) :: testName
        real*8 alfa, kappa(2), kappaCap(2), mu(2), lambda(2), h
        complex*16 alfa_c, sigma(2), sigmaCap(2), object1, object2, t(4)
            ! Here we test if U0+Uminus equals Uplus at z = -h
            if (testName == 'U') then
                open(1, file='test.txt', FORM='FORMATTED')
                write(1, *) "%Here we test if U0+Uminus equals Uplus at z = -h "
                write(1, *) "% Objects 1, 2 should coinside "
                write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
                do alfa=-100, 100d0, 0.5d0
                    alfa_c = cmplx(alfa)
                    sigma = makeSigma(kappa, alfa_c)
                    sigmaCap = makeSigma(kappaCap, alfa_c)
                    t(1) = CramDelta(1, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
                    t(1) = t(1) + CramDelta(1, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
                    
                    t(2) = CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
                    t(2) = t(2) + CramDelta(2, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
                    
                    t(3) = CramDelta(3, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
                    t(3) = t(3) + CramDelta(3, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)
                    
                    t(4) = CramDelta(4, 1, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(1)*h)
                    t(4) = t(4) + CramDelta(4, 2, kappa, kappaCap, lambda, mu, alfa_c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa_c)*exp(-sigma(2)*h)

                    ! U0 + Uminus
                    object1 = -ci*alfa*t(1)+sigma(2)*t(2) + U1(alfa_c, kappa, mu(1))*exp(-sigma(1)*h) + U2(alfa_c, kappa, mu(1))*exp(-sigma(2)*h)
                    ! Uplus
                    object2 = -ci*alfa*t(3)-sigmaCap(2)*t(4)


                    write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
                enddo 
                close(1)
            ! Here we test if -lambda*i*alfa*U + (lambda+2*mu)*W' equals Q2   
                                                                                    ! T O   B E   D O N E  !!!!!!!!!!!!
            else if (testName == 'W') then

                
            endif    
        
        END SUBROUTINE scattered_field
        
        
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
    
END MODULE test