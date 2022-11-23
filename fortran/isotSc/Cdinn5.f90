! ѕрограмма счета N интегралов по комплексному отрезку [a,b] 
! в том числе и дл€ полуоси (b = inf.);
! переменный шаг интегрировани€ в зависимости от быстроты
! получени€ требуемой точности на каждом предыдущем шаге 
! (считает CSIMPN5) 

	subroutine CDINN5(CF,a,b,eps,pr,N,sb,rd)
    implicit none

	integer N,ib,inf,ipr,i,it
	
	real*8 eps,pr,f(N),t,t1,t2,eps10,th,pm,pt

	complex*16 a,b,rd(N),r(N),h,sb(N),h1,a1,b1,s1

    common/simp5/h,ib,inf

	external CF

! initial constants

    ipr=0; t=abs(b-a); t2=(1-1d-12)*t
    eps10=10*eps;      if(t < 1d-13)return

    s1=(b-a)/t; h1=s1*pr; 

! beginning of integration
	
	b1=a; it=1

 	do while(it > 0)
	  a1=b1;  b1=b1+h1;  t1=abs(b1-a)
	  if(t1 > t)then
	    b1=b; it=-1
	  end if

      call CSIMPN5(CF,a1,b1,eps,N,sb,R); ib=-1
	  if (minval(abs(R)).lt.1d-8.and.abs(a1).gt.110.0d0.and.maxval(abs(R)).lt.1d-6) then
!		print*,minval(abs(R)),maxval(abs(R)),abs(a1)
		goto 10
	  endif

! for cheking-up the convergence at infinity

      if(inf > 0)then
	    th=abs(b1-a1)
		do i=1,N
		  f(i)=abs(r(i))/th
		end do
	  end if ! inf > 0

! summing

      do i=1,N
        rd(i)=rd(i)+r(i)
	  end do

      if(abs(h1) < 10*abs(h)) then
	    h1=4*h1;  else;   h1=4*h
	  end if

! at infinity

      if(inf > 0)then
	    pm=0
        do i=1,N
          if(abs(rd(i)) > 1d-15)then
		    pt=abs(f(i)*t1/rd(i))
		    if(pt > pm) pm=pt
		  end if
	    end do

        if(pm < eps10)then
          ipr=ipr+1
		  if(ipr > 4)return
		else
		  ipr=0
	    end if ! pm < eps10
	  end if ! inf >0

    end do ! while(it > 0)

10    end
