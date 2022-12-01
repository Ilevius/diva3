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
    real*8 singleW, singlePsi, psiMin, psiNumber, psiStep, singleR, Rmin, Rstep, Rmax, zSingle, xMin, xMax, xStep
    character(len=3) :: field
    character (len=5) :: mode
    namelist /study/field, mode, singleW, singlePsi, psiMin, psiNumber, psiStep, singleR, Rmin, Rstep, Rmax, zSingle, xMin, xMax, xStep
    
    real*8 t1, t2, t3, t4, tm, tp, eps, step, IntLimit
    namelist /dinn5Settings/ t1, t2, t3, t4, tm, tp, eps, step, IntLimit
    real*8 mu(2), lambda(2), kappa(2), kappaCap(2), w, currentPsi, currentR, currentX, currentZ
    real*8, allocatable:: x(:), z(:), psis(:), Rs(:)
    complex*16, allocatable:: integral_vert(:), integral_horis(:), stPhase(:,:)
    
        open(unit=1,file='input.txt',status='old')
        read(1, media) 
        read(1, study)
        read(1, dinn5Settings)
        close(1)
        
        kappa(1) = singleW/cp(1); kappa(2) = singleW/cs(1);
        kappaCap(1) = singleW/cp(2); kappaCap(2) = singleW/cs(2);
        mu(1) = rho(1)*cs(1)**2; mu(2) = rho(2)*cs(2)**2;
        lambda(1) = rho(1)*(cp(1)**2 - 2d0*cs(1)**2); lambda(2) = rho(2)*(cp(2)**2 - 2d0*cs(2)**2); 
        
        
        
        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0+1d0);
        

        
        
        !                                                            T E S T S                  
        !                                                              * * * 
                    
        !call source_field('first boundary', kappa, mu(1))  
        
        !call scattered_field('U', kappa, kappaCap, mu, lambda, h)
        
        !                                                           * * *
        !                                                   END OF  T E S T S
              
        call makeStudy
        
        
        open(1, file='integral.txt', FORM='FORMATTED')
        
        write(1,*) "% the field is ", field
        write(1,*) "%psi, R, x, z, re(u), im(u), abs(u), re(w), im(w), abs(w)"
        
        open(2, file='stPhase.txt', FORM='FORMATTED')
        write(2,*) "% the field is ", field
        write(2,*) "%psi, R, x, z, re(u), im(u), abs(u), re(w), im(w), abs(w)"
        do i = 1, pointsNumber
            write(1, '(10E15.6E3)') psis(i), Rs(i), x(i), z(i), real(integral_horis(i)), imag(integral_horis(i)), abs(integral_horis(i)), real(integral_vert(i)), imag(integral_vert(i)), abs(integral_vert(i))
            write(2, '(10E15.6E3)') psis(i), Rs(i), x(i), z(i), real(stPhase(1,i)),      imag(stPhase(1,i)),      abs(stPhase(1,i)),      real(stPhase(2,i)),     imag(stPhase(2,i)),     abs(stPhase(2,i))
        enddo  
        close(1); 
        close(2);
     
        
        
    CONTAINS
        SUBROUTINE makeStudy   
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
            
            else if (mode == 'flat_') then
                pointsNumber = ceiling((xMax-xMin)/xStep)
                allocate(x(pointsNumber), z(pointsNumber), psis(pointsNumber), Rs(pointsNumber), integral_vert(pointsNumber), integral_horis(pointsNumber), stPhase(2,pointsNumber));
                do i = 1, pointsNumber
                    x(i) = xMin+(i-1)*xStep
                    z(i) = zSingle
                    Rs(i) = sqrt(x(i)**2+4d0*h**2)
                    if (x(i)>0d0) then
                        psis(i) = atan(2d0*h/x(i))
                    else if (x(i)<0) then
                        psis(i) = atan(2d0*h/x(i)) + pi
                    else 
                        psis(i) = pi/2d0
                    endif
                    psis(i) = psis(i)*180/pi
                enddo 
                call xzPointsTest(x, z, psis, Rs, h, pointsNumber)
            endif
        
        
            if (field == "upp") then          
                call dinn5(upp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                call dinn5(upp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                call upp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)            
            else if (field == "ups") then
                call dinn5(ups_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                call dinn5(ups_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                call ups_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)
            else if (field == "usp") then 
                call dinn5(usp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                call dinn5(usp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                call usp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase) 
            else if (field == "uss") then
                call dinn5(uss_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                call dinn5(uss_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                call uss_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)
            else 
                print*, "Check field name in input file!"; pause;
            end if    
            
        END SUBROUTINE makeStudy
    
    
    !                                                   *   *   *   
    !                                                      u_pp
        SUBROUTINE upp_integrand_horis(alfa, s, n)
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
        END SUBROUTINE upp_integrand_horis
        
    
        SUBROUTINE upp_integrand_vert(alfa, s, n)
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
        END SUBROUTINE upp_integrand_vert
    

        SUBROUTINE upp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda, lambdaCap, mu, muCap, theResult)
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
        END SUBROUTINE upp_stPhase  
        
        
        !                                                   *   *   *   
        !                                                      u_ps
        SUBROUTINE ups_integrand_horis(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus 
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa)*(sigma(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(2, 1, kappa, kappaCap, lambda, mu, -alfa)*(sigma(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*exp(-sigma(2)*(z(i)+h)-sigma(1)*h-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(-sigma(2)*(z(i)+h)-sigma(1)*h+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE ups_integrand_horis
        
    
        SUBROUTINE ups_integrand_vert(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus
            sigma = makeSigma(kappa, alfa)
            
            
            no_x_part_plus = CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa)*(-ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa)
            no_x_part_minus = CramDelta(2, 1, kappa, kappaCap, lambda, mu, -alfa)*(ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            
            
            do i = 1, n
                s(i) = no_x_part_plus*exp(-sigma(2)*(z(i)+h)-sigma(1)*h-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(-sigma(2)*(z(i)+h)-sigma(1)*h+ci*alfa*x(i))
               s(i) = s(i)/(2d0*pi)
            enddo    
        END SUBROUTINE ups_integrand_vert
    

        SUBROUTINE ups_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda, lambdaCap, mu, muCap, theResult)
        implicit none
        integer pointsNumber
        real*8 psis(pointsNumber), Rs(pointsNumber), kappa(2), kappaCap(2), lambda, lambdaCap, mu, muCap
        real*8 stPoints(10), theta, D2theta, alfa0
        complex*16 theResult(2, pointsNumber), test  
        integer i, stPointsNumber
        complex*16 alfa0c, sigma(2), commonRes
            do i = 1, pointsNumber
                currentPsi = psis(i)
                currentR = Rs(i)
                ! Находим стац точки, сохраняем первую в дей и компл форме
                call Halfc(DThetaPSHalfC, -kappa(1), kappa(1), 2d-3, 1d-8, 10, stPoints, stPointsNumber)
                alfa0 = stPoints(1)
                alfa0c = cmplx(stPoints(1))
                sigma = makeSigma(kappa, alfa0c)
                
                if (stPointsNumber > 1 ) then
                    print*, 'oops'
                    pause
                endif  
                
                if (abs(alfa0)>kappa(1)) then
                    print*, 'oops'
                    pause
                endif
                
                ! находим значение фазовой функции и ее производных в стац точке
                theta = ThetaPS(alfa0, h, Rs(i), kappa, psis(i))
                D2theta = D2ThetaPS(alfa0, h, Rs(i), kappa, psis(i))
                
                ! obtaining asymptotics 
                test = CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)
                commonRes = sqrt(1d0/(2d0*pi*Rs(i)))*sqrt(1d0/abs(D2theta))*CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)*exp( ci*Rs(i)*theta + ci*pi/4d0*sign(1d0, D2Theta) )
                theResult(1,i) = commonRes*(sigma(2))
                theResult(2,i) = commonRes*(-ci*alfa0)    
            enddo      
        END SUBROUTINE ups_stPhase  
        
        FUNCTION DThetaPSHalfC(alfa)
            implicit none
            real*8 alfa
            complex*16 DThetaPSHalfC
                DThetaPSHalfC = cmplx(DThetaPS(alfa, h, singleR, kappa, currentPsi))
        END FUNCTION DThetaPSHalfC
        
        
        
        !                                                     *  *  *    
        !                                                       u_sp 
        
        SUBROUTINE usp_integrand_horis(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus 
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(1, 2, kappa, kappaCap, lambda, mu, alfa)*(-ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(1, 2, kappa, kappaCap, lambda, mu, -alfa)*(ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
            s(i) = no_x_part_plus*exp(-sigma(1)*(z(i)+h)-sigma(2)*h-ci*alfa*x(i))   
            s(i) = s(i) + no_x_part_minus*exp(-sigma(1)*(z(i)+h)-sigma(2)*h+ci*alfa*x(i))
            s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE usp_integrand_horis
        
    
        SUBROUTINE usp_integrand_vert(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus  
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(1, 2, kappa, kappaCap, lambda, mu, alfa)*(-sigma(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa)
            no_x_part_minus = CramDelta(1, 2, kappa, kappaCap, lambda, mu, -alfa)*(-sigma(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
            s(i) = no_x_part_plus*exp(-sigma(1)*(z(i)+h)-sigma(2)*h-ci*alfa*x(i))   
            s(i) = s(i) + no_x_part_minus*exp(-sigma(1)*(z(i)+h)-sigma(2)*h+ci*alfa*x(i))
            s(i) = s(i)/(2d0*pi)
            enddo    
        END SUBROUTINE usp_integrand_vert
        
    
        SUBROUTINE usp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda, lambdaCap, mu, muCap, theResult)
        implicit none
        integer pointsNumber
        real*8 psis(pointsNumber), Rs(pointsNumber), kappa(2), kappaCap(2), lambda, lambdaCap, mu, muCap
        real*8 stPoints(10), theta, D2theta, alfa0, test
        complex*16 theResult(2, pointsNumber)  
        integer i, stPointsNumber
        complex*16 alfa0c, sigma(2), commonRes
            do i = 1, pointsNumber
                currentPsi = psis(i)
                
                ! Находим стац точки, сохраняем первую в дей и компл форме
                call Halfc(DThetaSPHalfC, -100d0, 100d0, 2d-3, 1d-6, 10, stPoints, stPointsNumber)
                alfa0 = stPoints(1)
                alfa0c = cmplx(stPoints(1))
                sigma = makeSigma(kappa, alfa0c)
                
                if (stPointsNumber > 1 ) then
                    print*, 'oops'
                    pause
                endif  
                
                if (abs(alfa0)>kappa(1)) then
                    print*, 'oops'
                    pause
                endif
                
                ! находим значение фазовой функции и ее производных в стац точке
                theta = ThetaSP(alfa0, h, Rs(i), kappa, psis(i))
                D2theta = D2ThetaSP(alfa0, h, Rs(i), kappa, psis(i))
                
                ! obtaining asymptotics 
                commonRes = sqrt(1d0/(2d0*pi*Rs(i)))*sqrt(1d0/abs(D2theta))*CramDelta(1, 2, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)*exp( ci*Rs(i)*theta + ci*pi/4d0*sign(1d0, D2Theta) )
                theResult(1,i) = commonRes*(-ci*alfa0)
                theResult(2,i) = commonRes*(-sigma(1))    
            enddo      
        END SUBROUTINE usp_stPhase  
        
        FUNCTION DThetaSPHalfC(alfa)
            implicit none
            real*8 alfa
            complex*16 DThetaSPHalfC
                DThetaSPHalfC = cmplx(DThetaSP(alfa, h, singleR, kappa, currentPsi))
        END FUNCTION DThetaSPHalfC
        
        
        
        !                                                     *  *  *    
        !                                                       u_ss 
        
        SUBROUTINE uss_integrand_horis(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus 
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(2, 2, kappa, kappaCap, lambda, mu, alfa)*(sigma(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(2, 2, kappa, kappaCap, lambda, mu, -alfa)*(sigma(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
            s(i) = no_x_part_plus*exp(-sigma(2)*(z(i)+2d0*h)-ci*alfa*x(i))   
            s(i) = s(i) + no_x_part_minus*exp(-sigma(2)*(z(i)+2d0*h)+ci*alfa*x(i))
            s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE uss_integrand_horis
        
    
        SUBROUTINE uss_integrand_vert(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus  
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(2, 2, kappa, kappaCap, lambda, mu, alfa)*(-ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa)
            no_x_part_minus = CramDelta(2, 2, kappa, kappaCap, lambda, mu, -alfa)*(ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
            s(i) = no_x_part_plus*exp(-sigma(2)*(z(i)+2d0*h)-ci*alfa*x(i))   
            s(i) = s(i) + no_x_part_minus*exp(-sigma(2)*(z(i)+2d0*h)+ci*alfa*x(i))
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
                alfa0 = -cosd(psis(i))*kappa(2)
                sigma = makeSigma(kappa, alfa0)
                commonRes = sqrt(kappa(2)*sind(psis(i))**2/(2d0*pi*Rs(i)))*CramDelta(2, 2, kappa, kappaCap, lambda, mu, alfa0)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0)*exp(ci*Rs(i)*kappa(2) - ci*pi/4)
                theResult(1,i) = commonRes*(sigma(2))
                theResult(2,i) = commonRes*(-ci*alfa0)    
            enddo      
        END SUBROUTINE uss_stPhase 
            
    END PROGRAM isotSc

