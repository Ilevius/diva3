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
    real*8 cp(2), cs(2), rho(2), h, fieldCode
    namelist /media/ cp, cs, rho, h
    real*8 singleW, wMin, wStep, wMax, singlePsi, psiMin, psiNumber, psiStep, singleR, Rmin, Rstep, Rmax, zSingle, xMin, xMax, xStep, xSingle
    character(len=3) :: field
    character (len=5) :: mode
    namelist /study/field, mode, singleW, wMin, wStep, wMax, singlePsi, psiMin, psiNumber, psiStep, singleR, Rmin, Rstep, Rmax, zSingle, xMin, xMax, xStep, xSingle
    
    
    real*8 t1, t2, t3, t4, tm, tp, eps, step, IntLimit
    namelist /dinn5Settings/ t1, t2, t3, t4, tm, tp, eps, step, IntLimit
    real*8 mu(2), lambda(2), kappa(2), kappaCap(2), w, currentPsi, currentR, currentX, currentZ
    real*8 sWaveLen1, sWaveLen2, pWaveLen1, pWaveLen2
    real*8, allocatable:: x(:), z(:), psis(:), Rs(:), ws(:)
    complex*16, allocatable:: integral_vert(:), integral_horis(:), stPhase(:,:), integral_vert_all(:,:), integral_horis_all(:,:), stPhase_all(:,:,:)
    
        open(unit=1,file='input.txt',status='old')
        read(1, media) 
        read(1, study)
        read(1, dinn5Settings)
        close(1)
        
        kappa(1) = singleW/cp(1); kappa(2) = singleW/cs(1);
        kappaCap(1) = singleW/cp(2); kappaCap(2) = singleW/cs(2);
        mu(1) = rho(1)*cs(1)**2; mu(2) = rho(2)*cs(2)**2;
        lambda(1) = rho(1)*(cp(1)**2 - 2d0*cs(1)**2); lambda(2) = rho(2)*(cp(2)**2 - 2d0*cs(2)**2); 
        sWaveLen1 = 2d0*pi*cs(1)/singleW
        sWaveLen2 = 2d0*pi*cs(2)/singleW
        pWaveLen1= 2d0*pi*cp(1)/singleW
        pWaveLen2 = 2d0*pi*cp(2)/singleW
        !IntLimit = Kappa(2)*1.4d0+1d0
        
        
        
        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0);
        
        call makePoints
        
        
        !                                                            T E S T S                  
        !                                                              * * * 
        
        !call source_field('second boundary', kappa, mu(1), lambda(1))        
        !call reflected_field('first Lame', kappa, kappaCap, mu, lambda, h)
        
        !                                                           * * *
        !                                                   END OF  T E S T S
              
        call makeStudy
        
        call saveResults
        
        

     
        
        
    CONTAINS
    !                               D A T A S E T   P R E P A R I N G
    
        SUBROUTINE makePoints   
            if (mode == 'polar') then
                pointsNumber = psiNumber
                allocate(x(pointsNumber), z(pointsNumber), psis(pointsNumber), Rs(pointsNumber), integral_vert(pointsNumber), integral_horis(pointsNumber), stPhase(2,pointsNumber));
                allocate(integral_vert_all(4,pointsNumber), integral_horis_all(4,pointsNumber), stPhase_all(4,2,pointsNumber));

                ! строим точки по заданным в input полярным координатам
                do i = 1, pointsNumber
                    x(i) = singleR*cosd((i - 1)*psiStep + psiMin) 
                    z(i) = singleR*sind((i - 1)*psiStep + psiMin) + zSingle
                enddo
                ! получаем координаты точек во внутренней полярной системе координат
                do i = 1, pointsNumber
                    Rs(i) = sqrt(x(i)**2 + (z(i)+2d0*h)**2)
                    psis(i) = atan2((z(i)+2d0*h), x(i)) + 2d0*pi
                    psis(i) = psis(i)*180/pi
                enddo 
            
            else if (mode == 'decar') then
                pointsNumber = ceiling((xMax-xMin)/xStep)
                allocate(x(pointsNumber), z(pointsNumber), psis(pointsNumber), Rs(pointsNumber), integral_vert(pointsNumber), integral_horis(pointsNumber), stPhase(2,pointsNumber));
                allocate(integral_vert_all(4,pointsNumber), integral_horis_all(4,pointsNumber), stPhase_all(4,2,pointsNumber));
                
                !!Выражаем геометрические параметры в длинах поперечной волны соответствующего слоя
                !if (zSingle<-h) then
                !    zSingle = zSingle*sWaveLen2
                !    xMin = xMin*sWaveLen2
                !    xStep = xStep*sWaveLen2
                !else
                !    zSingle = zSingle*sWaveLen1
                !    xMin = xMin*sWaveLen1
                !    xStep = xStep*sWaveLen1
                !endif
                !

                

                
                do i = 1, pointsNumber
                    x(i) = xMin+(i-1)*xStep
                    z(i) = zSingle
                    Rs(i) = sqrt(x(i)**2 + (z(i)+2d0*h)**2)
                    psis(i) = atan2((z(i)+2d0*h), x(i)) + 2d0*pi
                    psis(i) = psis(i)*180/pi
                enddo 
                
                
                do i = 1, pointsNumber
                    Rs(i) = sqrt(x(i)**2 + (z(i)+2d0*h)**2)
                    psis(i) = atan2((z(i)+2d0*h), x(i)) + 2d0*pi
                    psis(i) = psis(i)*180/pi
                enddo 
                
                
            else if (mode == 'freqs') then
                pointsNumber = ceiling((wMax-wMin)/wStep)
                allocate(x(1), z(1), psis(1), Rs(1), integral_vert(pointsNumber), integral_horis(pointsNumber), stPhase(2,pointsNumber), ws(pointsNumber));
                allocate(integral_vert_all(4,pointsNumber), integral_horis_all(4,pointsNumber), stPhase_all(4,2,pointsNumber));
                do i = 1, pointsNumber
                    ws(i) = wMin+(i-1)*wStep
                enddo
                
                !Выражаем геометрические параметры в длинах поперечной волны соответствующего слоя
                if (zSingle<-h) then
                    x(1) = xSingle!*sWaveLen2
                    z(1) = zSingle!*sWaveLen2
                else
                    x(1) = xSingle!*sWaveLen1
                    z(1) = zSingle!*sWaveLen1
                endif
                
                
                Rs(1) = sqrt(x(1)**2 + (z(1)+2d0*h)**2)
                psis(1) = atan2((z(1)+2d0*h), x(1)) + 2d0*pi
                psis(1) = psis(1)*180/pi
                
            else
                print*, 'Wrong mode key!'
                pause
            endif
            
            
            call xzPointsTest(x, z, psis, Rs, h, pointsNumber)
        END SUBROUTINE makePoints

    
    
    !                               T H E   S T U D Y   I T S E L F
        SUBROUTINE makeStudy           
        
            if (field == "upp") then
                fieldCode = 11
                if (mode == "freqs") then
                    do i = 1, pointsNumber
                        singleW = ws(i)
                        kappa(1) = ws(i)/cp(1); kappa(2) = ws(i)/cs(1);
                        kappaCap(1) = ws(i)/cp(2); kappaCap(2) = ws(i)/cs(2);
                        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0);
                        call dinn5(upp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_horis(i))
                        call dinn5(upp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_vert(i))
                        call upp_stPhase(1, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase(:,i))  
                     enddo   
                else
                    call dinn5(upp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                    call dinn5(upp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                    call upp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)  
                endif
                
            else if (field == "usp") then
                fieldCode = 12
                if (mode == "freqs") then
                    do i = 1, pointsNumber
                        singleW = ws(i)
                        kappa(1) = ws(i)/cp(1); kappa(2) = ws(i)/cs(1);
                        kappaCap(1) = ws(i)/cp(2); kappaCap(2) = ws(i)/cs(2);
                        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0);
                        call dinn5(ups_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_horis(i))
                        call dinn5(ups_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_vert(i))
                        call upp_stPhase(1, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase(:,i))  
                     enddo   
                else
                    call dinn5(ups_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                    call dinn5(ups_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                    call ups_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase) 
                endif
                
                

            else if (field == "ups") then 
                fieldCode = 21
                
                if (mode == "freqs") then
                    do i = 1, pointsNumber
                        singleW = ws(i)
                        kappa(1) = ws(i)/cp(1); kappa(2) = ws(i)/cs(1);
                        kappaCap(1) = ws(i)/cp(2); kappaCap(2) = ws(i)/cs(2);
                        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0);
                        call dinn5(usp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_horis(i))
                        call dinn5(usp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_vert(i))
                        call usp_stPhase(1, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase(:,i))  
                     enddo   
                else
                    call dinn5(usp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                    call dinn5(usp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                    call usp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase) 
                endif
                
                
            else if (field == "uss") then
                fieldCode = 22
                
                if (mode == "freqs") then
                    do i = 1, pointsNumber
                        singleW = ws(i)
                        kappa(1) = ws(i)/cp(1); kappa(2) = ws(i)/cs(1);
                        kappaCap(1) = ws(i)/cp(2); kappaCap(2) = ws(i)/cs(2);
                        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0);
                        call dinn5(uss_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_horis(i))
                        call dinn5(uss_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_vert(i))
                        call uss_stPhase(1, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase(:,i))  
                     enddo   
                else
                    call dinn5(uss_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                    call dinn5(uss_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                    call uss_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)
                endif

                
            else if (field == "all") then
                call dinn5(upp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert_all(1,:))
                call dinn5(upp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis_all(1,:))
                call upp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase_all(1,:,:))
                
                call dinn5(ups_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert_all(2,:))
                call dinn5(ups_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis_all(2,:))
                call ups_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase_all(2,:,:))
                
                call dinn5(usp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert_all(3,:))
                call dinn5(usp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis_all(3,:))
                call usp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase_all(3,:,:))
                
                call dinn5(uss_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert_all(4,:))
                call dinn5(uss_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis_all(4,:))
                call uss_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase_all(4,:,:))
                
                integral_vert = integral_vert_all(1,:) + integral_vert_all(2,:) + integral_vert_all(3,:) + integral_vert_all(4,:)
                integral_horis = integral_horis_all(1,:)+integral_horis_all(2,:)+integral_horis_all(3,:)+integral_horis_all(4,:)
                stPhase = stPhase_all(1,:,:)+stPhase_all(2,:,:)+stPhase_all(3,:,:)+stPhase_all(4,:,:)
            
            
            else if (field == "ppp") then
                fieldCode = 11
                
                if (mode == "freqs") then
                    do i = 1, pointsNumber
                        singleW = ws(i)
                        kappa(1) = ws(i)/cp(1); kappa(2) = ws(i)/cs(1);
                        kappaCap(1) = ws(i)/cp(2); kappaCap(2) = ws(i)/cs(2);
                        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0);
                        call dinn5(u_plus_pp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_horis(i))
                        call dinn5(u_plus_pp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_vert(i))
                        call u_plus_pp_stPhase(1, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase(:,i))  
                     enddo   
                else
                    call dinn5(u_plus_pp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                    call dinn5(u_plus_pp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                    call u_plus_pp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)
                endif
 
                
                
            else if (field == "pps") then
                fieldCode = 12
                
                if (mode == "freqs") then
                    do i = 1, pointsNumber
                        singleW = ws(i)
                        kappa(1) = ws(i)/cp(1); kappa(2) = ws(i)/cs(1);
                        kappaCap(1) = ws(i)/cp(2); kappaCap(2) = ws(i)/cs(2);
                        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0);
                        call dinn5(u_plus_ps_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_horis(i))
                        call dinn5(u_plus_ps_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_vert(i))
                        call u_plus_ps_stPhase(1, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase(:,i))  
                     enddo   
                else
                    call dinn5(u_plus_ps_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                    call dinn5(u_plus_ps_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                    call u_plus_ps_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)
                endif

                
               
            else if (field == "psp") then
                fieldCode = 21
                
                if (mode == "freqs") then
                    do i = 1, pointsNumber
                        singleW = ws(i)
                        kappa(1) = ws(i)/cp(1); kappa(2) = ws(i)/cs(1);
                        kappaCap(1) = ws(i)/cp(2); kappaCap(2) = ws(i)/cs(2);
                        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0);
                        call dinn5(u_plus_sp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_horis(i))
                        call dinn5(u_plus_sp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_vert(i))
                        call u_plus_sp_stPhase(1, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase(:,i))  
                     enddo   
                else
                    call dinn5(u_plus_sp_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                    call dinn5(u_plus_sp_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                    call u_plus_sp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)    
                endif
                                
            else if (field == "pss") then
                fieldCode = 22
                
                if (mode == "freqs") then
                    do i = 1, pointsNumber
                        singleW = ws(i)
                        kappa(1) = ws(i)/cp(1); kappa(2) = ws(i)/cs(1);
                        kappaCap(1) = ws(i)/cp(2); kappaCap(2) = ws(i)/cs(2);
                        t1 = Kappa(1)*0.5; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0);
                        call dinn5(u_plus_ss_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_horis(i))
                        call dinn5(u_plus_ss_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,1,integral_vert(i))
                        call u_plus_ss_stPhase(1, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase(:,i))  
                     enddo   
                else
                    call dinn5(u_plus_ss_integrand_horis,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_horis)
                    call dinn5(u_plus_ss_integrand_vert,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,pointsNumber,integral_vert)
                    call u_plus_ss_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda(1), lambda(2), mu(1), mu(2), stPhase)    
                endif

                

                
            else 
                print*, "Check field name in input file!"; pause;
            end if           
        END SUBROUTINE makeStudy
        
        
        
        SUBROUTINE saveResults
            if (mode == 'polar') then
                open(1, file='integral_polar.txt', FORM='FORMATTED')
                open(2, file='stPhase_polar.txt', FORM='FORMATTED')
        
                write(1,*) "% the field is ", field
                write(1,'(A)') "% Dinn5 settings , t1, t2, t3, t4, tm, tp, eps, step, IntLimit"
                write(1,'((A),9E15.6E3)') "% ", t1, t2, t3, t4, tm, tp, eps, step, IntLimit
                write(1,'(A)') "% 1)psi, 2) R, 3) x, 4) z, 5) re(u), 6) im(u), 7) abs(u), 8) re(w), 9) im(w), 10) abs(w), 11) fieldCode, 12) cp(1), 13) cp(2), 14) cs(1), 15) cs(2), 16) rho(1), 17) rho(2), 18) h, 19) w"
        
                write(2,*) "% the field is ", field
                write(2,'(A)') "% Dinn5 settings , t1, t2, t3, t4, tm, tp, eps, step, IntLimit"
                write(2,'((A),9E15.6E3)') "% ", t1, t2, t3, t4, tm, tp, eps, step, IntLimit
                write(2,'(A)') "%psi, R, x, z, re(u), im(u), abs(u), re(w), im(w), abs(w), fieldCode, cp(1), cp(2), cs(1), cs(2), rho(1), rho(2), h, w"
        
                do i = 1, pointsNumber
                    write(1, '(19E15.6E3)') (i - 1)*psiStep + psiMin, singleR, x(i), z(i), real(integral_horis(i)), imag(integral_horis(i)), abs(integral_horis(i)), real(integral_vert(i)), imag(integral_vert(i)), abs(integral_vert(i)), fieldCode, cp(1), cp(2), cs(1), cs(2), rho(1), rho(2), h, singleW
                    write(2, '(19E15.6E3)') (i - 1)*psiStep + psiMin, singleR, x(i), z(i), real(stPhase(1,i)),      imag(stPhase(1,i)),      abs(stPhase(1,i)),      real(stPhase(2,i)),     imag(stPhase(2,i)),     abs(stPhase(2,i)), fieldCode, cp(1), cp(2), cs(1), cs(2), rho(1), rho(2), h, singleW
                enddo  
                close(1); 
                close(2);
                
            else if (mode == 'decar') then
                !           перевод в длины волн x, z
                !       ОСТОРОЖНО! Может не сработать корректно, если глубина будет задана меньше 5-10 длинн волны
                !
                !
                !if (zSingle<-h) then
                !    x = x/sWaveLen2
                !    z = z/sWaveLen2
                !else
                !    x = x/sWaveLen1
                !    z = z/sWaveLen1
                !endif
                
                
                
                open(1, file='integral.txt', FORM='FORMATTED')
                open(2, file='stPhase.txt', FORM='FORMATTED')
        
                write(1,*) "% the field is ", field
                write(1,'(A)') "% Dinn5 settings , t1, t2, t3, t4, tm, tp, eps, step, IntLimit"
                write(1,'((A),9E15.6E3)') "% ", t1, t2, t3, t4, tm, tp, eps, step, IntLimit
                write(1,'(A)') "% 1)psi, 2) R, 3) x, 4) z, 5) re(u), 6) im(u), 7) abs(u), 8) re(w), 9) im(w), 10) abs(w), 11) fieldCode, 12) cp(1), 13) cp(2), 14) cs(1), 15) cs(2), 16) rho(1), 17) rho(2), 18) h, 19) w"
        
                write(2,*) "% the field is ", field
                write(2,'(A)') "% Dinn5 settings , t1, t2, t3, t4, tm, tp, eps, step, IntLimit"
                write(2,'((A),9E15.6E3)') "% ", t1, t2, t3, t4, tm, tp, eps, step, IntLimit
                write(2,'(A)') "%psi, R, x, z, re(u), im(u), abs(u), re(w), im(w), abs(w), fieldCode, cp(1), cp(2), cs(1), cs(2), rho(1), rho(2), h, w"
        
                do i = 1, pointsNumber
                    write(1, '(19E15.6E3)') psis(i), Rs(i), x(i), z(i), real(integral_horis(i)), imag(integral_horis(i)), abs(integral_horis(i)), real(integral_vert(i)), imag(integral_vert(i)), abs(integral_vert(i)), fieldCode, cp(1), cp(2), cs(1), cs(2), rho(1), rho(2), h, singleW
                    write(2, '(19E15.6E3)') psis(i), Rs(i), x(i), z(i), real(stPhase(1,i)),      imag(stPhase(1,i)),      abs(stPhase(1,i)),      real(stPhase(2,i)),     imag(stPhase(2,i)),     abs(stPhase(2,i)), fieldCode, cp(1), cp(2), cs(1), cs(2), rho(1), rho(2), h, singleW
                enddo  
                close(1); 
                close(2);
                
            else if (mode == 'freqs') then
                open(1, file='integral.txt', FORM='FORMATTED')
                open(2, file='stPhase.txt', FORM='FORMATTED')
        
                write(1,*) "% the field is ", field
                write(1,'(A)') "% Dinn5 settings , t1, t2, t3, t4, tm, tp, eps, step, IntLimit"
                write(1,'((A),9E15.6E3)') "% ", t1, t2, t3, t4, tm, tp, eps, step, IntLimit
                write(1,'(A)') "% 1)psi, 2) R, 3) x, 4) z, 5) re(u), 6) im(u), 7) abs(u), 8) re(w), 9) im(w), 10) abs(w), 11) fieldCode, 12) cp(1), 13) cp(2), 14) cs(1), 15) cs(2), 16) rho(1), 17) rho(2), 18) h, 19) w"
        
                write(2,*) "% the field is ", field
                write(2,'(A)') "% Dinn5 settings , t1, t2, t3, t4, tm, tp, eps, step, IntLimit"
                write(2,'((A),9E15.6E3)') "% ", t1, t2, t3, t4, tm, tp, eps, step, IntLimit
                write(2,'(A)') "%psi, R, x, z, re(u), im(u), abs(u), re(w), im(w), abs(w), fieldCode, cp(1), cp(2), cs(1), cs(2), rho(1), rho(2), h, w"
        
                do i = 1, pointsNumber
                    write(1, '(19E15.6E3)') psis(1), Rs(1), x(1), z(1), real(integral_horis(i)), imag(integral_horis(i)), abs(integral_horis(i)), real(integral_vert(i)), imag(integral_vert(i)), abs(integral_vert(i)), fieldCode, cp(1), cp(2), cs(1), cs(2), rho(1), rho(2), h, ws(i)
                    write(2, '(19E15.6E3)') psis(1), Rs(1), x(1), z(1), real(stPhase(1,i)),      imag(stPhase(1,i)),      abs(stPhase(1,i)),      real(stPhase(2,i)),     imag(stPhase(2,i)),     abs(stPhase(2,i)), fieldCode, cp(1), cp(2), cs(1), cs(2), rho(1), rho(2), h, ws(i)
                enddo  
                close(1); 
                close(2);

            endif
        END SUBROUTINE saveResults
    
    
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
                call Halfc(DThetaPSHalfC, -kappa(1), kappa(1), 1d-3, 1d-8, 10, stPoints, stPointsNumber)
                alfa0 = stPoints(1)
                alfa0c = cmplx(stPoints(1))
                sigma = makeSigma(kappa, alfa0c)              
                if (stPointsNumber > 1 ) then
                    print*, 'oops, we have a few stationary points'
                    pause
                endif           
                if (abs(alfa0)>kappa(1)) then
                    print*, 'oops, we have alfa>kappa1'
                    pause
                endif
        ! находим значение фазовой функции и ее производных в стац точке
                theta = ThetaPS(alfa0, h, Rs(i), kappa, psis(i))
                D2theta = D2ThetaPS(alfa0, h, Rs(i), kappa, psis(i))  
        ! obtaining asymptotics 
                commonRes = sqrt(1d0/(2d0*pi*Rs(i)))*sqrt(1d0/abs(D2theta))*CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)*exp( ci*Rs(i)*theta + ci*pi/4d0*sign(1d0, D2Theta) )
                theResult(1,i) = commonRes*(sigma(2))
                theResult(2,i) = commonRes*(-ci*alfa0)    
            enddo      
        END SUBROUTINE ups_stPhase  
        
        FUNCTION DThetaPSHalfC(alfa)
            implicit none
            real*8 alfa
            complex*16 DThetaPSHalfC
                DThetaPSHalfC = cmplx(DThetaPS(alfa, h, currentR, kappa, currentPsi))
        END FUNCTION DThetaPSHalfC
        
        
        !                           N E W     ST   P H A S E       PS
        
        SUBROUTINE ups_stPhaseNew(pointsNumber, psis, Rs, kappa, kappaCap, lambda, lambdaCap, mu, muCap, theResult)
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
                call Halfc(DThetaPSHalfCNew, -kappa(1), kappa(1), 2d-3, 1d-8, 10, stPoints, stPointsNumber)
                alfa0 = stPoints(1)
                alfa0c = cmplx(stPoints(1))
                sigma = makeSigma(kappa, alfa0c)              
                if (stPointsNumber > 1 ) then
                    print*, 'oops, we have a few stationary points'
                    pause
                endif           
                if (abs(alfa0)>kappa(1)) then
                    print*, 'oops, we have alfa>kappa1'
                    pause
                endif
        ! находим значение фазовой функции и ее производных в стац точке
                theta = ThetaPSNew(alfa0, h, Rs(i), 40d0, kappa, psis(i))
                D2theta = D2ThetaPSNew(alfa0, h, Rs(i), 40d0, kappa, psis(i))  
        ! obtaining asymptotics 
                commonRes = sqrt(1d0/(2d0*pi*Rs(i)))*sqrt(1d0/abs(D2theta))*CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)*exp( ci*Rs(i)*theta + ci*pi/4d0*sign(1d0, D2Theta) )
                theResult(1,i) = commonRes*(sigma(2))
                theResult(2,i) = commonRes*(-ci*alfa0)    
            enddo      
        END SUBROUTINE ups_stPhaseNew  
        
        FUNCTION DThetaPSHalfCNew(alfa)
            implicit none
            real*8 alfa
            complex*16 DThetaPSHalfCNew
                DThetaPSHalfCNew = cmplx(DThetaPSNew(alfa, h, currentR, 40d0, kappa, currentPsi))
        END FUNCTION DThetaPSHalfCNew
        
        
        
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
        integer i, stPointsNumber, j
        complex*16 alfa0c, sigma(2), commonRes
            do i = 1, pointsNumber
                currentPsi = psis(i)
                currentR = Rs(i)
                ! Находим стац точки, сохраняем первую в дей и компл форме
                call Halfc(DThetaSPHalfC, -100d0, 100d0, 1d-3, 1d-6, 10, stPoints, stPointsNumber)
                if (stPointsNumber > 1 ) then
                    print*, 'oops, we have a few stationary points: ', stPointsNumber, stPoints(1), stPoints(2)
                    pause
                endif             
                theResult(1,i) = 0d0
                theResult(2,i) = 0d0       
                do j = 1, stPointsNumber
                    alfa0 = stPoints(j)
                    alfa0c = cmplx(stPoints(j))
                    sigma = makeSigma(kappa, alfa0c)
                    if (abs(alfa0)>kappa(1)) then
                        print*, 'oops, we have alfa>kappa1'
                        pause
                    endif
                    ! находим значение фазовой функции и ее производных в стац точке
                    theta = ThetaSP(alfa0, h, Rs(i), kappa, psis(i))
                    D2theta = D2ThetaSP(alfa0, h, Rs(i), kappa, psis(i))
                    ! obtaining asymptotics 
                    commonRes = sqrt(1d0/(2d0*pi*Rs(i)))*sqrt(1d0/abs(D2theta))*CramDelta(1, 2, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)*exp( ci*Rs(i)*theta + ci*pi/4d0*sign(1d0, D2Theta) )
                    theResult(1,i) = theResult(1,i) + commonRes*(-ci*alfa0)
                    theResult(2,i) = theResult(2,i) + commonRes*(-sigma(1))    
                enddo   
            enddo      
        END SUBROUTINE usp_stPhase  
        
        
        
        
        FUNCTION DThetaSPHalfC(alfa)
            implicit none
            real*8 alfa
            complex*16 DThetaSPHalfC
                DThetaSPHalfC = cmplx(DThetaSP(alfa, h, currentR, kappa, currentPsi))
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
        
        
        
 !                                        ______ ____   ___  __  __  __  ___  ___ __ ______ ______  ____ ____  
 !                                        | || | || \\ // \\ ||\ || (( \ ||\\//|| || | || | | || | ||    || \\ 
 !                                          ||   ||_// ||=|| ||\\||  \\  || \/ || ||   ||     ||   ||==  ||  ))
 !                                          ||   || \\ || || || \|| \_)) ||    || ||   ||     ||   ||___ ||_//              f i e l d
 !                                                      https://patorjk.com/software/taag/                                                                   


        
!                                                                             __ __ ____  ____ 
 !                                                                            || || || \\ || \\
!                                                                             || || ||_// ||_//
 !                                                                            \\_// ||    ||   
!                  
        
        SUBROUTINE u_plus_pp_integrand_horis(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigmaCap(2), sigma(2), no_x_part_minus, no_x_part_plus 
            sigmaCap = makeSigma(kappaCap, alfa)
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(3, 1, kappa, kappaCap, lambda, mu, alfa)*(-ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(3, 1, kappa, kappaCap, lambda, mu, -alfa)*(ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*exp(sigmaCap(1)*(z(i)+h)-h*sigma(1)-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(sigmaCap(1)*(z(i)+h)-h*sigma(1)+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE u_plus_pp_integrand_horis
        
        SUBROUTINE u_plus_pp_integrand_vert(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigmaCap(2), sigma(2), no_x_part_minus, no_x_part_plus 

            sigmaCap = makeSigma(kappaCap, alfa)
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(3, 1, kappa, kappaCap, lambda, mu, alfa)*(sigmaCap(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(3, 1, kappa, kappaCap, lambda, mu, -alfa)*(sigmaCap(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*exp(sigmaCap(1)*(z(i)+h)-h*sigma(1)-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(sigmaCap(1)*(z(i)+h)-h*sigma(1)+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE u_plus_pp_integrand_vert
        
        
        SUBROUTINE u_plus_pp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda, lambdaCap, mu, muCap, theResult)
        implicit none
        integer pointsNumber
        real*8 psis(pointsNumber), Rs(pointsNumber), kappa(2), kappaCap(2), lambda, lambdaCap, mu, muCap
        real*8 stPoints(10), theta, D2theta, alfa0
        complex*16 theResult(2, pointsNumber), test  
        integer i, stPointsNumber
        complex*16 alfa0c, sigma(2), sigmaCap(2), commonRes
            do i = 1, pointsNumber
                currentPsi = psis(i)
                currentR = Rs(i)
        ! Находим стац точки, сохраняем первую в дей и компл форме
                call Halfc(DEtaPPHalfC, -100d0, 100d0, 2d-3, 1d-8, 10, stPoints, stPointsNumber)
                alfa0 = stPoints(1)
                alfa0c = cmplx(stPoints(1))
                sigma = makeSigma(kappa, alfa0c)
                sigmaCap = makeSigma(kappaCap, alfa0c)
                if (stPointsNumber > 1 ) then
                    print*, 'oops, we have a few stationary points'
                    pause
                endif           
                if (abs(alfa0)>kappa(1)) then
                    print*, 'oops, we have alfa>kappa1'
                    pause
                endif
        ! находим значение фазовой функции и ее производных в стац точке
                theta = EtaPP(alfa0, h, Rs(i), kappa, kappaCap, psis(i))
                D2theta = D2EtaPP(alfa0, h, Rs(i), kappa, kappaCap, psis(i))  
        ! obtaining asymptotics 
                commonRes = sqrt(1d0/(2d0*pi*Rs(i)))*sqrt(1d0/abs(D2theta))*CramDelta(3, 1, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)*exp( ci*Rs(i)*theta + ci*pi/4d0*sign(1d0, D2Theta) )
                theResult(1,i) = commonRes*(-ci*alfa0)
                theResult(2,i) = commonRes*(sigmaCap(1))    
            enddo      
        END SUBROUTINE u_plus_pp_stPhase  
        
        FUNCTION DEtaPPHalfC(alfa)
            implicit none
            real*8 alfa
            complex*16 DEtaPPHalfC
                DEtaPPHalfC = cmplx(DEtaPP(alfa, h, currentR, kappa, kappaCap, currentPsi))
        END FUNCTION DEtaPPHalfC
        
        

!                                                                         __ __ ____   __ 
!                                                                         || || || \\ (( \
!                                                                         || || ||_//  \\ 
!                                                                         \\_// ||    \_))

        
        SUBROUTINE u_plus_ps_integrand_horis(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigmaCap(2), sigma(2), no_x_part_minus, no_x_part_plus 
            sigmaCap = makeSigma(kappaCap, alfa)
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(3, 2, kappa, kappaCap, lambda, mu, alfa)*(-ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(3, 2, kappa, kappaCap, lambda, mu, -alfa)*(ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*exp(sigmaCap(1)*(z(i)+h)-h*sigma(2)-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(sigmaCap(1)*(z(i)+h)-h*sigma(2)+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE u_plus_ps_integrand_horis
        
        SUBROUTINE u_plus_ps_integrand_vert(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigmaCap(2), sigma(2), no_x_part_minus, no_x_part_plus 
            sigmaCap = makeSigma(kappaCap, alfa)
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(3, 2, kappa, kappaCap, lambda, mu, alfa)*(sigmaCap(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(3, 2, kappa, kappaCap, lambda, mu, -alfa)*(sigmaCap(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*exp(sigmaCap(1)*(z(i)+h)-h*sigma(2)-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(sigmaCap(1)*(z(i)+h)-h*sigma(2)+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE u_plus_ps_integrand_vert
        
        
        SUBROUTINE u_plus_ps_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda, lambdaCap, mu, muCap, theResult)
        implicit none
        integer pointsNumber
        real*8 psis(pointsNumber), Rs(pointsNumber), kappa(2), kappaCap(2), lambda, lambdaCap, mu, muCap
        real*8 stPoints(10), theta, D2theta, alfa0
        complex*16 theResult(2, pointsNumber), test  
        integer i, stPointsNumber
        complex*16 alfa0c, sigma(2), sigmaCap(2), commonRes
            do i = 1, pointsNumber
                currentPsi = psis(i)
                currentR = Rs(i)
        ! Находим стац точки, сохраняем первую в дей и компл форме
                call Halfc(DEtaPSHalfC, -100d0, 100d0, 2d-3, 1d-8, 10, stPoints, stPointsNumber)
                alfa0 = stPoints(1)
                alfa0c = cmplx(stPoints(1))
                sigma = makeSigma(kappa, alfa0c)
                sigmaCap = makeSigma(kappaCap, alfa0c)
                if (stPointsNumber > 1 ) then
                    print*, 'oops, we have a few stationary points'
                    pause
                endif           
                if (abs(alfa0)>kappaCap(1)) then
                    print*, 'oops, we have alfa>kappa1'
                    !pause
                endif
        ! находим значение фазовой функции и ее производных в стац точке
                theta = EtaPS(alfa0, h, Rs(i), kappa, kappaCap, psis(i))
                D2theta = D2EtaPS(alfa0, h, Rs(i), kappa, kappaCap, psis(i))  
        ! obtaining asymptotics 
                commonRes = sqrt(1d0/(2d0*pi*Rs(i)))*sqrt(1d0/abs(D2theta))*CramDelta(3, 2, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)*exp( ci*Rs(i)*theta + ci*pi/4d0*sign(1d0, D2Theta) )
                theResult(1,i) = commonRes*(-ci*alfa0)
                theResult(2,i) = commonRes*(sigmaCap(1))    
            enddo      
        END SUBROUTINE u_plus_ps_stPhase  
        
        FUNCTION DEtaPSHalfC(alfa)
            implicit none
            real*8 alfa
            complex*16 DEtaPSHalfC
                DEtaPSHalfC = cmplx(DEtaPS(alfa, h, currentR, kappa, kappaCap, currentPsi))
        END FUNCTION DEtaPSHalfC
              
       
!                                                                     __ __  __  ____ 
 !                                                                    || || (( \ || \\
!                                                                     || ||  \\  ||_//
 !                                                                    \\_// \_)) ||   
                 
        SUBROUTINE u_plus_sp_integrand_horis(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigmaCap(2), sigma(2), no_x_part_minus, no_x_part_plus 
            sigmaCap = makeSigma(kappaCap, alfa)
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(4, 1, kappa, kappaCap, lambda, mu, alfa)*(-sigmaCap(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(4, 1, kappa, kappaCap, lambda, mu, -alfa)*(-sigmaCap(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*exp(sigmaCap(2)*(z(i)+h)-h*sigma(1)-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(sigmaCap(2)*(z(i)+h)-h*sigma(1)+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE u_plus_sp_integrand_horis
        
        SUBROUTINE u_plus_sp_integrand_vert(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigmaCap(2), sigma(2), no_x_part_minus, no_x_part_plus 
            sigmaCap = makeSigma(kappaCap, alfa)
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(4, 1, kappa, kappaCap, lambda, mu, alfa)*(-ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(4, 1, kappa, kappaCap, lambda, mu, -alfa)*(ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*exp(sigmaCap(2)*(z(i)+h)-h*sigma(1)-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(sigmaCap(2)*(z(i)+h)-h*sigma(1)+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE u_plus_sp_integrand_vert
        
        
        SUBROUTINE u_plus_sp_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda, lambdaCap, mu, muCap, theResult)
        implicit none
        integer pointsNumber
        real*8 psis(pointsNumber), Rs(pointsNumber), kappa(2), kappaCap(2), lambda, lambdaCap, mu, muCap
        real*8 stPoints(10), theta, D2theta, alfa0
        complex*16 theResult(2, pointsNumber), test  
        integer i, stPointsNumber
        complex*16 alfa0c, sigma(2), sigmaCap(2), commonRes
            do i = 1, pointsNumber
                currentPsi = psis(i)
                currentR = Rs(i)
        ! Находим стац точки, сохраняем первую в дей и компл форме
                call Halfc(DEtaSPHalfC, -100d0, 100d0, 2d-3, 1d-8, 10, stPoints, stPointsNumber)
                alfa0 = stPoints(1)
                alfa0c = cmplx(stPoints(1))
                sigma = makeSigma(kappa, alfa0c)
                sigmaCap = makeSigma(kappaCap, alfa0c)
                if (stPointsNumber > 1 ) then
                    print*, 'oops, we have a few stationary points'
                    pause
                endif           
                if (abs(alfa0)>kappa(1)) then
                    print*, 'oops, we have alfa>kappa1'
                    pause
                endif
        ! находим значение фазовой функции и ее производных в стац точке
                theta = EtaSP(alfa0, h, Rs(i), kappa, kappaCap, psis(i))
                D2theta = D2EtaSP(alfa0, h, Rs(i), kappa, kappaCap, psis(i))  
        ! obtaining asymptotics 
                commonRes = sqrt(1d0/(2d0*pi*Rs(i)))*sqrt(1d0/abs(D2theta))*CramDelta(4, 1, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)*exp( ci*Rs(i)*theta + ci*pi/4d0*sign(1d0, D2Theta) )
                theResult(1,i) = commonRes*(sigmaCap(2))
                theResult(2,i) = commonRes*(-ci*alfa0)    
            enddo      
        END SUBROUTINE u_plus_sp_stPhase  
        
        FUNCTION DEtaSPHalfC(alfa)
            implicit none
            real*8 alfa
            complex*16 DEtaSPHalfC
                DEtaSPHalfC = cmplx(DEtaSP(alfa, h, currentR, kappa, kappaCap, currentPsi))
        END FUNCTION DEtaSPHalfC
        
        
!                                                                    __ __  __   __ 
 !                                                                   || || (( \ (( \
!                                                                    || ||  \\   \\ 
 !                                                                   \\_// \_)) \_))
!                
  
         
        SUBROUTINE u_plus_ss_integrand_horis(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigmaCap(2), sigma(2), no_x_part_minus, no_x_part_plus 
            sigmaCap = makeSigma(kappaCap, alfa)
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(4, 2, kappa, kappaCap, lambda, mu, alfa)*(-sigmaCap(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(4, 2, kappa, kappaCap, lambda, mu, -alfa)*(-sigmaCap(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*exp(sigmaCap(2)*(z(i)+h)-h*sigma(2)-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(sigmaCap(2)*(z(i)+h)-h*sigma(2)+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE u_plus_ss_integrand_horis
        
        SUBROUTINE u_plus_ss_integrand_vert(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigmaCap(2), sigma(2), no_x_part_minus, no_x_part_plus 
            sigmaCap = makeSigma(kappaCap, alfa)
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(4, 2, kappa, kappaCap, lambda, mu, alfa)*(-ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(4, 2, kappa, kappaCap, lambda, mu, -alfa)*(ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*exp(sigmaCap(2)*(z(i)+h)-h*sigma(2)-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*exp(sigmaCap(2)*(z(i)+h)-h*sigma(2)+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE u_plus_ss_integrand_vert
        
        
        SUBROUTINE u_plus_ss_stPhase(pointsNumber, psis, Rs, kappa, kappaCap, lambda, lambdaCap, mu, muCap, theResult)
        implicit none
        integer pointsNumber
        real*8 psis(pointsNumber), Rs(pointsNumber), kappa(2), kappaCap(2), lambda, lambdaCap, mu, muCap
        real*8 stPoints(10), theta, D2theta, alfa0
        complex*16 theResult(2, pointsNumber), test  
        integer i, stPointsNumber
        complex*16 alfa0c, sigma(2), sigmaCap(2), commonRes
            do i = 1, pointsNumber
                currentPsi = psis(i)
                currentR = Rs(i)
        ! Находим стац точки, сохраняем первую в дей и компл форме
                call Halfc(DEtaSSHalfC, -200d0, 200d0, 2d-3, 1d-8, 10, stPoints, stPointsNumber)
                alfa0 = stPoints(1)
                alfa0c = cmplx(stPoints(1))
                sigma = makeSigma(kappa, alfa0c)
                sigmaCap = makeSigma(kappaCap, alfa0c)
                if (stPointsNumber > 1 ) then
                    print*, 'oops, we have a few stationary points'
                    pause
                endif           
                if (abs(alfa0)>kappa(1)) then
                    print*, 'oops, we have alfa>kappa1'
                    !pause
                endif
        ! находим значение фазовой функции и ее производных в стац точке
                theta = EtaSS(alfa0, h, Rs(i), kappa, kappaCap, psis(i))
                D2theta = D2EtaSS(alfa0, h, Rs(i), kappa, kappaCap, psis(i))  
        ! obtaining asymptotics 
                commonRes = sqrt(1d0/(2d0*pi*Rs(i)))*sqrt(1d0/abs(D2theta))*CramDelta(4, 2, kappa, kappaCap, lambda, mu, alfa0c)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa0c)*exp( ci*Rs(i)*theta + ci*pi/4d0*sign(1d0, D2Theta) )
                theResult(1,i) = commonRes*(sigmaCap(2))
                theResult(2,i) = commonRes*(-ci*alfa0)    
            enddo      
        END SUBROUTINE u_plus_ss_stPhase  
        
            FUNCTION DEtaSSHalfC(alfa)
            implicit none
            real*8 alfa
            complex*16 DEtaSSHalfC
                DEtaSSHalfC = cmplx(DEtaSS(alfa, h, currentR, kappa, kappaCap, currentPsi))
        END FUNCTION DEtaSSHalfC
        
        

!                                                _________ _______  _______ _________ _______ 
!                                                \__   __/(  ____ \(  ____ \\__   __/(  ____ \
!                                                   ) (   | (    \/| (    \/   ) (   | (    \/
!                                                   | |   | (__    | (_____    | |   | (_____ 
!                                                   | |   |  __)   (_____  )   | |   (_____  )
!                                                   | |   | (            ) |   | |         ) |
!                                                   | |   | (____/\/\____) |   | |   /\____) |
!                                                   )_(   (_______/\_______)   )_(   \_______)
                                             


                                                         
        
        !                                                           T E S T S    
        
        SUBROUTINE reflected_field(testName, kappa, kappaCap, mu, lambda, h)
        use functions
        IMPLICIT NONE
        character (len=*) :: testName
        real*8 alfa, kappa(2), kappaCap(2), mu(2), lambda(2), h
        complex*16 alfa_c, sigma(2), sigmaCap(2), object1, object2, t(4)
        complex*16 output1(1), output2(1), output3(1)
                
            !                                           Here we test if first Lame equation holds for Uminus   

            if (testName == 'first Lame') then
                open(1, file='test.txt', FORM='FORMATTED')
                write(1, *) "%Here we test if W0+Wminus equals Wplus at z = -h "
                write(1, *) "% Objects 1, 2 should coinside "
                write(1, *) "%alfa, real(object1), imag(object1), real(object2), imag(object2)"
                x = 0d0
                z = 0d0
                do alfa=-6d0, 6d0, 0.05d0
                    alfa_c = cmplx(alfa)
                    sigma = makeSigma(kappa, alfa_c)
                    sigmaCap = makeSigma(kappaCap, alfa_c)
                    
                    call upp_integrand_horis(alfa_c, output1, 1)                   
                    call upp_integrand_vert_d(alfa_c, output2, 1)
                    call upp_integrand_horis_d2(alfa_c, output3, 1)
                    
                    object1 = output1(1)*(kappa(2)**2*mu(1)-alfa**2*(lambda(1)+2d0*mu(1)))
                    object2 =  ci*alfa*(lambda(1)+mu(1))*output2(1) - mu(1)*output3(1)
                    write(1, '(5E15.6E3)') alfa, real(object1), imag(object1), real(object2), imag(object2)
                enddo 
                close(1)     
            endif   
            
        END SUBROUTINE reflected_field
        
        
        
        SUBROUTINE upp_integrand_horis_d2(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus 
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(1, 1, kappa, kappaCap, lambda, mu, alfa)*(-ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(1, 1, kappa, kappaCap, lambda, mu, -alfa)*(ci*alfa)/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*sigma(1)**2*exp(-sigma(1)*(z(i)+2d0*h)-ci*alfa*x(i))   
                !s(i) = s(i) + no_x_part_minus*sigma(1)**2*exp(-sigma(1)*(z(i)+2d0*h)+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE upp_integrand_horis_d2
        
        
        SUBROUTINE upp_integrand_vert_d(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus  
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(1, 1, kappa, kappaCap, lambda, mu, alfa)*(-sigma(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa)
            no_x_part_minus = CramDelta(1, 1, kappa, kappaCap, lambda, mu, -alfa)*(-sigma(1))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
            s(i) = no_x_part_plus*(-sigma(1))*exp(-sigma(1)*(z(i)+2d0*h)-ci*alfa*x(i))   
            !s(i) = s(i) + no_x_part_minus*(-sigma(1))*exp(-sigma(1)*(z(i)+2d0*h)+ci*alfa*x(i))
            s(i) = s(i)/(2d0*pi)
            enddo    
        END SUBROUTINE upp_integrand_vert_d
        
        
        SUBROUTINE ups_integrand_horis_d2(alfa, s, n)
        implicit none;
        integer n, i
        complex*16 alfa, s(n), sigma(2), no_x_part_minus, no_x_part_plus 
            sigma = makeSigma(kappa, alfa)
            no_x_part_plus = CramDelta(2, 1, kappa, kappaCap, lambda, mu, alfa)*(sigma(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, alfa) 
            no_x_part_minus = CramDelta(2, 1, kappa, kappaCap, lambda, mu, -alfa)*(sigma(2))/CramDelta(0, 0, kappa, kappaCap, lambda, mu, -alfa)
            do i = 1, n
                s(i) = no_x_part_plus*sigma(2)**2*exp(-sigma(2)*(z(i)+h)-sigma(1)*h-ci*alfa*x(i))   
                s(i) = s(i) + no_x_part_minus*sigma(2)**2*exp(-sigma(2)*(z(i)+h)-sigma(1)*h+ci*alfa*x(i))
                s(i) = s(i)/(2d0*pi)
            enddo     
        END SUBROUTINE ups_integrand_horis_d2

        

        
        
            
    END PROGRAM isotSc

