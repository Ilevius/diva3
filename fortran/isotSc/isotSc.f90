!  isotSc.f90 
!
!  FUNCTIONS:
!  isotSc - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: isotSc
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

PROGRAM isotSc
    use functions
    use test
    implicit none
    integer pointsNumber, i
    real*8 cp(2), cs(2), rho(2), h
    namelist /media/ cp, cs, rho, h
    real*8 singleF, singlePsi, psiMin, psiNumber, psiStep, singleR, Rmin, Rstep, Rmax
    character(len=5) :: mode
    namelist /study/mode, singleF, singlePsi, psiMin, psiNumber, psiStep, singleR, Rmin, Rstep, Rmax
    
    real*8 t1, t2, t3, t4, tm, tp, eps, step, IntLimit
    namelist /dinn5Settings/ t1, t2, t3, t4, tm, tp, eps, step, IntLimit
    real*8 mu(2), lambda(2), kappa(2), kappaCap(2), w
    real*8, allocatable:: x(:), z(:), psis(:), Rs(:)
    complex*16, allocatable:: integral_vert(:), integral_horis(:), stPhase(:,:)
    
        open(unit=1,file='input.txt',status='old')
        read(1, media) 
        read(1, study)
        read(1, dinn5Settings)
        close(1)
        
        kappa(1) = singleF/cp(1); kappa(2) = singleF/cs(1);
        kappaCap(1) = singleF/cp(2); kappaCap(2) = singleF/cs(2);
        
        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0+1d0);
        
        mu(1) = rho(1)*cs(1)**2; mu(2) = rho(2)*cs(2)**2;
        lambda(1) = rho(1)*(cp(1)**2 - 2d0*cs(1)**2); lambda(2) = rho(2)*(cp(2)**2 - 2d0*cs(2)**2); 
        
        !                        AS                                  T E S T S                  BEGIN
        !                         YOU'LL SEE JUST WHO YOU ARE           * * * 
                    
        !call source_field('first boundary', kappa, mu(1))  
        
        !call scattered_field('U', kappa, kappaCap, mu, lambda, h)
        
        !                                                           * * *
        !                                                   END OF  T E S T S
            
   
        if (mode == 'polar') then
            pointsNumber = psiNumber
            allocate(x(pointsNumber), z(pointsNumber), psis(pointsNumber), Rs(pointsNumber), integral_vert(pointsNumber), integral_horis(pointsNumber), stPhase(2,pointsNumber));
            open(1, file='points.txt', FORM='FORMATTED')
            do i = 1, pointsNumber
                psis(i) = (i - 1)*psiStep + psiMin
                Rs(i) = singleR
                x(i) = singleR*cosd(psis(i)) 
                z(i) = singleR*sind(psis(i)) - 2d0*h
                write(1,*) x(i), z(i)
            enddo
            close(1)
            
        else if (mode == 'flat') then
            pointsNumber = psiNumber
            allocate(x(pointsNumber), z(pointsNumber));
            do i = 1, pointsNumber
                x(i) = singleR*cosd((i - 1)*psiStep + psiMin) 
                z(i) = singleR*sind((i - 1)*psiStep + psiMin) - 2d0*h
            enddo    
        endif
            

                 
        call dinn5(uss_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
        call dinn5(uss_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
        
        call uss_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)
        
        open(1, file='integral.txt', FORM='FORMATTED')
        write(1,*) "%psi, re(u), im(u), abs(u), re(w), im(w), abs(w)"
        open(2, file='stPhase.txt', FORM='FORMATTED')
        write(2,*) "%psi, re(u), im(u), abs(u), re(w), im(w), abs(w)"
        do i = 1, pointsNumber
            write(1, '(F7.2, 6E15.6E3)') psis(i), real(integral_horis(i)), imag(integral_horis(i)), abs(integral_horis(i)), real(integral_vert(i)), imag(integral_vert(i)), abs(integral_vert(i))
            write(2, '(F7.2, 6E15.6E3)') psis(i), real(stPhase(1,i)), imag(stPhase(1,i)), abs(stPhase(1,i)), real(stPhase(2,i)), imag(stPhase(2,i)), abs(stPhase(2,i))
        
        enddo  
        close(1); close(2);
        
    CONTAINS
    
        SUBROUTINE uss_integrand_horis(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus 
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(1, 1, kappa, kappaCap, lambda, mu, alfa)*(-ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(1, 1, kappa, kappaCap, lambda, mu, -alfa)*(ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
            s(i) = no_x_part_plus*exp(-sigma(1)*(z(i)+2d0*h)-ci*alfa*x(i))   
            s(i) = s(i) + no_x_part_minus*exp(-sigma(1)*(z(i)+2d0*h)+ci*alfa*x(i))
            s(i) = s(i)/(2d0*pi)
            enddo   
            
        END SUBROUTINE uss_integrand_horis
    
        SUBROUTINE uss_integrand_vert(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus  
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(1, 1, kappa, kappaCap, lambda, mu, alfa)*(-sigma(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa)
            no_x_part_minus = CramDelta(1, 1, kappa, kappaCap, lambda, mu, -alfa)*(-sigma(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
            s(i) = no_x_part_plus*exp(-sigma(1)*(z(i)+2d0*h)-ci*alfa*x(i))   
            s(i) = s(i) + no_x_part_minus*exp(-sigma(1)*(z(i)+2d0*h)+ci*alfa*x(i))
            s(i) = s(i)/(2d0*pi)
            enddo    
        END SUBROUTINE uss_integrand_vert
        
    

        SUBROUTINE uss_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda, lambdaCap, mu, muCap, theResult)
        implicit none
        integer pointsNumber
        real*8 psis(pointsNumber), Rs(pointsNumber), kappa(2), kappaCap(2), lambda, lambdaCap, mu, muCap
        complex*16 theResult(2, pointsNumber)  
        integer i
        complex*16 alfa0, sigma(2), commonRes
            do i = 1, pointsNumber
                alfa0 = -cosd(psis(i))*kappa(1)
                sigma = makeSigma(kappa, alfa0)
                commonRes = sqrt(kappa(1)*sind(psis(i))**2/(2d0*pi*Rs(i)))*CramDelta(1, 1, kappa, kappaCap, lambda, mu, alfa0)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0)*exp(ci*Rs(i)*kappa(1) - ci*pi/4)
                theResult(1,i) = commonRes*(-ci*alfa0)
                theResult(2,i) = commonRes*(-sigma(1))    
            enddo      
        END SUBROUTINE uss_stPhase       
            
    END PROGRAM isotSc

