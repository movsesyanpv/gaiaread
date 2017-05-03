MODULE GALACTIC
USE DFLIB

IMPLICIT NONE



 real(8), parameter, private :: Leo = 282.85948083 ! °
 real(8), parameter, private :: AL  = Leo-90.0_8 ! °
 real(8), parameter, private :: L0  =  32.931918056 ! ° 
 real(8), parameter, private :: si  = 0.88998807641_8 ! sin 62.871748611° 
 real(8), parameter, private :: ci  = 0.45598379779_8 ! cos 62.871748611° 
 
 real(8), parameter :: Pi  = 3.1415926535897932384626433832795_8


 real(8), parameter, private :: Leo1 = 282.25 ! °
 real(8), parameter, private :: L01  =  33 ! ° 
 real(8), parameter, private :: si1  = 0.88781538514_8 ! sin 62.6° 
 real(8), parameter, private :: ci1  = 0.46019978478_8 ! cos 62.6° 

 ! Прецессионные константы
REAL(8), PARAMETER :: M_1975 = 3.07374_8*15.0_8/3600.0_8, N_1975 = 20.0405_8/3600.0_8
REAL(8), PARAMETER :: M_1950 = 3.07327_8*15.0_8/3600.0_8, N_1950 = 20.0426_8/3600.0_8



CONTAINS





INTEGER(4) FUNCTION Round(x)
REAL(8), INTENT(IN) :: x

 Round=INT(ANINT(x))
 
END FUNCTION



 Subroutine Galaxy(a,d,l,b)
 real(8),intent(in) :: a,d
 real(8),intent(out) :: l,b
 real(8) :: al,sa,ca,sd,cd 

  al=a-Leo
  sa=dsind(al);  ca=dcosd(al)
  sd=dsind(d);   cd=dcosd(d)
  b=dasind(sd*ci-cd*si*sa)
  l=datan2d(sd*si+cd*ci*sa,cd*ca)+L0
  if (l<0) then 
   l=l+360.0_8
  endif
 end subroutine


Subroutine Galaxy1950(a,d,l,b)
 real(8),intent(in) :: a,d
 real(8),intent(out) :: l,b
 real(8) :: al,sa,ca,sd,cd 

  al=a-Leo1
  sa=dsind(al);  ca=dcosd(al)
  sd=dsind(d);   cd=dcosd(d)
  b=dasind(sd*ci1-cd*si1*sa)
  l=datan2d(sd*si1+cd*ci1*sa,cd*ca)+L01
  if (l<0) then 
   l=l+360.0_8
  endif
 end subroutine

 Subroutine Equatorial(l,b,a,d)
 real(8),intent(in) :: l,b
 real(8),intent(out) :: a,d
 real(8) cb,sb,l1,cl,sl
 l1=l-L0
 cb=dcosd(b);  sb=dsind(b)
 cl=dcosd(l1); sl=dsind(l1)

 d=dasind(cb*si*sl+sb*ci)
 a=datan2d(cb*cl,sb*si-cb*ci*sl)+AL  
 if (a>360.0_8) then
   a=a-360.0_8
 end if
 end subroutine



 Subroutine GalaxMu(mua,mud,l,b,d, mul,mub)
 real(8),intent(in) :: mua,mud,l,b,d 
 real(8),intent(out) :: mul,mub
 real(8) :: cd,sfi,cfi
   cd=dcosd(d);
   sfi=si*dcosd(l-L0)/cd
   cfi=(dcosd(b)*ci-dsind(b)*si*dsind(l-L0))/cd
   mul= cfi*mua+sfi*mud
   mub=-sfi*mua+cfi*mud
 end subroutine


 REAL(4) function Ranorm(s)
 ! s - дисперсия
 real(4), intent(in) :: s
 real(4) :: a
 integer, save :: z=37843827
 integer(4) ::   i;
  a=0.0
  DO i=1,12 
   a=a+ran(z)
  END DO
  Ranorm=(a-6.0)*s
 end FUNCTION RANORM

 REAL(8) FUNCTION PRECRA(RA,DE, Dt)
 ! Прецессия по прямому восхождению в градусах
 ! RA, DE - в градусах
 ! Dt - в годах
 REAL(8), INTENT(IN) :: RA,DE,Dt
  PRECRA=(M_1975+N_1975*DSIND(RA)*DTAND(DE))*Dt
 END FUNCTION 

 REAL(8) FUNCTION PRECDE(RA,DE, Dt)
 ! Прецессия по склонению восхождению в градусах
 REAL(8), INTENT(IN) :: RA,DE,Dt
  PRECDE=N_1975*DCOSD(RA)*Dt
 END FUNCTION 


 REAL(8) FUNCTION RHO(L1,B1,L2,B2)
 ! Расстояние между объектами в градусах
 REAL(8), INTENT(IN) :: L1,B1,L2,B2
  RHO=DACOSD( DSIND(B1)*DSIND(B2)+ DCOSD(B1)*DCOSD(B2)*DCOSD(L1-L2) )
 END FUNCTION


END MODULE