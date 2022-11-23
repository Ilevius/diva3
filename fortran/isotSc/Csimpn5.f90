! Программа счета N интегралов методом Симпсона 
! по комплексному отрезку [a,b] 

! Вход:
!  subroutine СF(u,s,N) - подынтегральные функции: 
!         u - аргумент (complex*16),
!         s(N) - массив значений функций в точке u (complex*16),
!         N - число интегралов (integer);
! eps -  отн. погрешность, real*8
! N- число интегралов (integer) 

! Bыход:
! R(N) -  массив значений интегралов
 
    subroutine CSIMPN5(CF,a,b,eps,N,sb,R)
    implicit none
    
	integer N,Nh,i,k,ib,inf
	
	real*8 eps,pm,t

	complex*16 a,b,R(N),s(N),s1(N),s2(N),s3(N), &
	           h,h1,h3,h43,sb(N),ri,u

    common/simp5/h,ib,inf
	 
! first step

    h=(b-a)/2; h3=h/3; h43=4*h3; Nh=1
    
	if(ib > 0)then
	  call CF(a,s1,N)
    else
	  do i=1,N
	    s1(i)=sb(i)
      end do
	end if ! ib

	call CF(b  ,sb,N)
	call CF(a+h,s3,N)

	do i=1,N
	  s1(i)=(s1(i)+sb(i))*h3
	  s2(i)=0
	  s3(i)=s3(i)*h43
	end do ! i

! next steps

 1  h1=h; h=h/2; h43=4*h/3; Nh=2*Nh; u=a+h

	do i=1,N
	  s1(i)=s1(i)/2
	  s2(i)=s2(i)/2+s3(i)/4
	  s3(i)=0
	end do ! i

    do k=1,Nh
	  call CF(u,s,N)
	  
	  do i=1,N
	    s3(i)=s3(i)+s(i)*h43
	  end do ! i

	  u=u+h1
	end do ! k

	pm=0
    do i=1,N
	  ri=s1(i)+s2(i)+s3(i)
      if(abs(ri) > 1d-12)then
	    t=abs((ri-R(i))/ri)
		if(t > pm) pm=t
	  end if ! ri > 1d-12
	  R(i)=ri   
	end do ! i

! accuracy control

    if((pm < eps))return
	if(abs(h) > eps) go to 1

    end
