Module HealPix

  IMPLICIT none

  INTEGER(4), private, PARAMETER :: ns_max=8192 ! 2^13 : largest nside available


  ! Numerical Constant (Double precision)
  REAL(8), PARAMETER, public :: QUARTPI=0.785398163397448309615660845819875721049_8
  REAL(8), PARAMETER, public :: HALFPI= 1.570796326794896619231321691639751442099_8
  REAL(8), PARAMETER, public :: PI    = 3.141592653589793238462643383279502884197_8
  REAL(8), PARAMETER, public :: TWOPI = 6.283185307179586476925286766559005768394_8
  REAL(8), PARAMETER, public :: FOURPI=12.56637061435917295385057353311801153679_8
  REAL(8), PARAMETER, public :: SQRT2 = 1.41421356237309504880168872420969807856967_8
  REAL(8), PARAMETER, public :: EULER = 0.5772156649015328606065120900824024310422_8
  REAL(8), PARAMETER, public :: SQ4PI_INV = 0.2820947917738781434740397257803862929220_8
  REAL(8), PARAMETER, public :: TWOTHIRD = 0.6666666666666666666666666666666666666666_8


   integer(4), private, save, dimension(128) :: x2pix=0,y2pix=0

CONTAINS
  !=======================================================================
  subroutine ang2pix_ring(nside, phi, theta, ipix)
    !=======================================================================
    !     renders the pixel number ipix (RING scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map
    !     resolution parameter nside
    !=======================================================================
    INTEGER(4), INTENT(IN) :: nside
    INTEGER(4), INTENT(OUT) :: ipix
    REAL(8), INTENT(IN) ::  theta, phi

    INTEGER(4) ::  nl4, jp, jm
    REAL(8)    ::  z, za, tt, tp, tmp, temp1, temp2
    INTEGER(4) ::  ir, ip, kshift

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) STOP "ang2pix_ring: nside out of range"
    if (theta<-HALFPI .or. theta>HALFPI)  then
       print *,"ANG2PIX_RING: theta : ",theta," is out of range [-Pi/2, +Pi/2]"
       STOP
    endif

    z = DCOS(HALFPI-theta)
    za = DABS(z)
    tt = MODULO( phi, twopi) / halfpi  ! in [0,4)


    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(.5_8+tt)
       temp2 = nside*.75_8*z
       jp = int(temp1-temp2) ! index of  ascending edge line
       jm = int(temp1+temp2) ! index of descending edge line

       ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
       kshift = 1 - modulo(ir,2) ! kshift=1 if ir even, 0 otherwise

       nl4 = 4*nside
       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*nside*(nside-1) + nl4*(ir-1) + ip

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_8)
       tmp = nside * DSQRT( 3.0_8*(1.0_8 - za) )

       jp = INT(tp          * tmp ) ! increasing edge line index
       jm = INT((1.0_8 - tp) * tmp ) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir )     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0._8) then
          ipix = 2*ir*(ir-1) + ip
       else
          ipix = 12*nside**2 - 2*ir*(ir+1) + ip
       endif

    endif

    return
  end subroutine ang2pix_ring

 

  !=======================================================================
  subroutine pix2ang_ring(nside, ipix, phi, theta)
    !=======================================================================
    !     renders theta and phi coordinates of the nominal pixel center
    !     for the pixel number ipix (RING scheme)
    !     given the map resolution parameter nside
    !=======================================================================
    INTEGER(KIND=4), INTENT(IN) :: ipix, nside
    REAL(KIND=8), INTENT(OUT) ::  theta, phi

    INTEGER(KIND=4) ::  nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
    REAL(KIND=8) ::  fodd, hip, fihip
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) STOP "nside out of range"
    npix = 12*nside**2       ! total number of points
    if (ipix <0 .or. ipix>npix-1) STOP "ipix out of range"

    ipix1 = ipix + 1 ! in {1, npix}
    nl2 = 2*nside
    ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1

    if (ipix1 <= ncap) then ! North Polar cap -------------

       hip   = ipix1*0.5_8
       fihip = AINT ( hip , kind=8)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
       iphi  = ipix1 - 2*iring*(iring - 1)

       theta = ACOS( 1.0_8 - iring**2 / (3.0_8*nside**2) )
       phi   = (real(iphi,8) - 0.5_8) * PI/(2.0_8*iring)

    elseif (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------

       ip    = ipix1 - ncap - 1
       nl4   = 4*nside
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = MODULO(ip,nl4) + 1

       fodd  = 0.5_8 * (1 + MODULO(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
       theta = ACOS( (nl2 - iring) / (1.5_8*nside) )
       phi   = (real(iphi,8) - fodd) * PI /(2.0_8*nside)

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix1 + 1
       hip   = ip*0.5_8
       fihip = AINT ( hip , kind=8)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

       theta = ACOS( -1.0_8 + iring**2 / (3.0_8*nside**2) )
       phi   = (real(iphi,8) - 0.5_8) * PI/(2.0_8*iring)

    endif
    
	theta=HALFPI-theta 

    return
  end subroutine pix2ang_ring ! pix2ang_ring


  END MODULE HEALPIX