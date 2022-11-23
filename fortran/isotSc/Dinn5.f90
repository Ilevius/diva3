! ѕрограмма вычислени€ N интегралов по полубесконечному контуру 
!           в случае обратной волны 
!
!  subroutine —F(u,s,n) - подынтегральные функции; u - аргумент (complex*16),
!                   s(n) - массив значений функций в точке u (complex*16),
!                   n - число интегралов (integer)
! [t1,t2],[t3,t4] - участки отклонени€ контура вниз (real*8)
!         [t2,t3] - участок отклонени€ контура вверх (real*8)
! tm,tp > 0 - величины отклонени€ контура вниз и вверх (real*8)
! (если нет обратной волны, то следует положить t2=t3=t1, tp=0 
!  обход полюсов при этом будет только снизу на участке [t1,t4]
!  с отклонением на tm)
! eps -  отн. погрешность,  pr - начальный шаг, 
! gr - верхний предел (все real*8)
! N- число интегралов (integer) 
! Rd(N) - выход: массив значений интегралов

    subroutine DINN5(CF,t1,t2,t3,t4,tm,tp,eps,pr,gr,N,Rd)
    implicit none
	
	integer N,ib,inf,i

	real*8 t1,t2,t3,t4,tm,tp,eps,pr,gr
	
	complex*16 Rd(N),a,b,sb(N),h

    external CF
    

    common/simp5/h,ib,inf

! initial constants

    h=pr/4; ib=1; inf=-1
	do i=1,N;  Rd(i)=0;	end do

! [0, t1]
    a=0; b=t1
    call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

	if(t3-t2 < eps)then
! no inverse poles case

! [t1, t1-i*tm]    

	  a=b; b=cmplx(t1,-tm)
      call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

!  [t1-i*tm, t4-i*tm]    

      a=b; b=cmplx(t4,-tm)
      call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

!  [t4-i*tm, t4]    

      a=b; b=t4
      call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

    else
! t2 < t3 - an inverse pole case

      if(t2-t1 > eps)then
! first deviation from below

!   [t1, t1-i*tm]    

	    a=b; b=cmplx(t1,-tm)
        call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

!   [t1-i*tm,t2-i*tm]    

	    a=b; b=cmplx(t2,-tm)
        call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

	  end if ! t2 > t1
! diviation from above

!   [b,t2+i*tp]    

      a=b; b=cmplx(t2,tp)
      call CDINN5(CF,a,b,eps,pr,N,sb,Rd)
	   
!  [t2+i*tp, t3+i*tp]    

      a=b; b=cmplx(t3,tp)
      call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

!  [t3+i*tp, t3-i*tm]    

	  a=b; b=cmplx(t3,-tm)
      call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

! second diviation from below

!  [t3-i*tm, t4-i*tm] 

	  a=b; b=cmplx(t4,-tm)
      call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

!  [t4-i*tm, t4] 

	  a=b; b=t4
      call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

	end if ! t3 > t2

! [t4, inf.] 

	a=b; b=gr; inf=1
    call CDINN5(CF,a,b,eps,pr,N,sb,Rd)

    end