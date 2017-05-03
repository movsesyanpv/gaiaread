
MODULE Projection

USE DFLIB

IMPLICIT NONE


REAL(8), PARAMETER :: PI = 3.141592653589793238462643383279502884197_8


CONTAINS
  

SUBROUTINE Aitoff(l,b,x,y)
 REAL(8), INTENT(IN) :: l,b  ! Сферические координаты в радианах 
 REAL(8), INTENT(OUT) :: x,y ! Декартовы координаты 
 REAL(8) :: s, l1  ! Вспомогательные переменные
 real(8) :: n = 900.0 !300.0 на 1366х768
 LOGICAL(4) status
 TYPE (windowconfig) wc
 IF (l>PI) THEN ! Приведение l в диапазон -Pi до +Pi 
   l1=l-2*Pi 
  ELSE 
   l1=l
 END IF
 status = GETWINDOWCONFIG(wc)
 S=sqrt(1.0+cos(b)*cos(l1/2)) ! Знаменатель формул (4.1)
 x=n*-2*cos(b)*sin(l1/2)/s + wc%numxpixels/2.0             !стало n*orig+960
 y=-n*sin(b)/s + wc%numypixels/2.0                      !n*orig+540
 END SUBROUTINE Aitoff


REAL(8) FUNCTION radi(x) ! Перевод градусов в радианы 
INTEGER, INTENT(IN) :: x
 radi=x/180.0*PI
END FUNCTION radi

REAL(8) FUNCTION rad(x) ! Перевод градусов в радианы 
REAL(8), INTENT(IN) :: x
 rad=x/180.0*PI
END FUNCTION rad

REAL(8) FUNCTION grad(x) ! Перевод радиан в градусы
REAL(8), INTENT(IN) :: x
 grad=x*180.0/PI
END FUNCTION grad


SUBROUTINE AitoffGrid (Step,Gr)

 INTEGER, INTENT(IN) :: Step  ! Шаг сетки в градусах 
 LOGICAL, INTENT(IN) :: Gr ! Флаг - в градусах или в часах разметка долготы

 INTEGER :: i,j ! Переменные циклов for 
 REAL(8) :: l,b ! Галактические координаты 
 REAL(8) :: x,y ! Декартовы координаты 
 CHARACTER(8)  s ! Строка для подписей 
 INTEGER ::  h  ! Для разметки осей 
 
 TYPE (wxycoord) wxy
 INTEGER(2) :: status2
 INTEGER(4) :: status4
 
 ! Нанесение сетки меридианов 
 
 status4 = SetColorRGB(#0000ff)

 DO i=-180,+180,Step
  l=radi(i) !  Перевод в радианы }
  DO j=-90,+90,5 ! Цикл построения вдоль меридиана 
   ! Вычисление точки меридиана 
   b=radi(j) ! Перевод в радианы широты 
   CALL Aitoff(l,b,x,y) ! Перевод в декартовы координаты 
   ! Если точка первая (j=-90), то помещаем графический курсор
   !  в точку (x,y) функцией MoveTo_W, если точка не первая, то
   ! "прочерчиваем" курсором линию из предыдущей точки 
   !  в точку (u,v) функцией LineTo_W.
   IF (j==-90) THEN
     CALL MoveTo_W(x,y,wxy) 
	ELSE 
	 status2=LineTo_W(x,y);
   END IF
  END DO ! j
 END DO ! i

 ! Нанесение сетки параллелей - аналогично предыдущему 
 DO j=-90,+90,Step
  b=radi(j)
  DO i=-180,+180,5  ! цикл построения вдоль параллели 
   l=radi(i)
   CALL Aitoff(l,b,x,y);
   if (i==-180) then 
     CALL MoveTo_W(x,y,wxy)
    else 
	 status2=LineTo_W(x,y)
   end if
  END DO ! i
 END DO ! j

 status2=SetFont('t''Arial''h48')
 status4 = SetColorRGB(#FF0000)

 ! Подписи меридианов вдоль экватора 
 DO i=-180,+180,Step
  ! Вычисление координаты точки вывода надписи 
  l=Radi(i);
  CALL Aitoff(l,0.0_8,x,y)
  ! Если Gr истина, то разметка в градусах, иначе в - часах 
  if (Gr) then  
    h=i
  else 
    h=i/15 
	if (h<0) h=h+24;
  end if
  write(s,'(I4)') h ! Преобразование значения h  в текстовую строку 
  Call MoveTo_W(x, y, wxy)
  Call OUTGTEXT(s)
 END DO


 ! Подписи параллелей вдоль нулевого меридиана - аналогично 
 DO j=-90,+90,Step
  if (j /= 0) then ! Экватор не подписываем 
   b=Radi(j);
   CALL Aitoff(0.0_8,b,x,y)
   write(s,'(I4)') j 
   Call MoveTo_W(x, y, wxy)
   Call OUTGTEXT(s)
  end if
 END DO

END SUBROUTINE AitoffGrid




END MODULE