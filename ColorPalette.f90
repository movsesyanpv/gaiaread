module ColorPalette
    
    USE DFLIB
    use Projection
    use HealPix
    
    integer, parameter :: aR = 0       ! RGB for our 1st color (blue in this case).
    integer, parameter :: aG = 0
    integer, parameter :: aB = 0  
    integer, parameter :: bR = 255     ! RGB for our 2nd color (yellow in this case).
    integer, parameter :: bG = 255
    integer, parameter :: bB = 255
    
    contains

    function GetColorGradVect(value,valueexp) result(res)
        
        real(8), dimension(:) :: value
        integer, dimension(size(value)) :: red, green, blue
        integer, dimension(size(value)) :: res
        real(8) :: valueexp
        
        value = value**valueexp
        
        red   = (bR - aR) * value + aR      ! Evaluated as -255*value + 255.
        green = (bG - aG) * value + aG      ! Evaluates as 0.
        blue  = (bB - aB) * value + aB      ! Evaluates as 255*value + 0.
        
        res = red + 256*green + 65536*blue
        
    end function GetColorGradVect
  
    
    
    function GetColorGrad(value,valueexp) result(res)
        
        real(8) :: value, valueexp
        integer :: red, green, blue
        integer :: res 
        
        value = value**valueexp
         
        red   = (bR - aR) * value + aR      ! Evaluated as -255*value + 255.
        green = (bG - aG) * value + aG      ! Evaluates as 0.
        blue  = (bB - aB) * value + aB      ! Evaluates as 255*value + 0.
        
        res = red + 256*green + 65536*blue
        
    end function GetColorGrad
    
    function GetWColorGrad(value) result(res)
        
        real(8) :: value, valueexp
        integer :: red, green, blue
        integer :: res 
         
        red   = 255 * value      ! Evaluated as -255*value + 255.
        green = 255 * value      ! Evaluates as 0.
        blue  = 255 * value      ! Evaluates as 255*value + 0.
        
        res = red + 256*green + 65536*blue
        
    end function GetWColorGrad
    
    subroutine DrawHealpix(ipixcolor,nside)
    
    integer, dimension(0:) :: ipixcolor
    real(8) :: l,b,x,y
    integer :: i,j,status,k
    
    !!$omp parallel private(b,l,x,y)
    !!$omp do
    do i = 0, 3500!*3600
        do j = 0, 4000!*3600
            l = -halfpi+halfpi*i/(1750)!*3600)
            b = -4.0*atan(1.0)+4.0*atan(1.0)*j/(2000)!*3600)
            !
            !call ang2pix_ring(nside, b, l, k)
            !status = setcolorrgb(ipixcolor(k))
            !call Aitoff(b,l,x,y)
            !print*, l, b
            !if((x.lt.0).or.(y.lt.0))then
            !    print*, l, b
            !    print*, x, y
            !    stop
            !endif
            !l = -l
            call ang2pix_ring(nside, b, l, k)
            call Aitoff(b,l,x,y)
            !!$omp critical
            status = setcolorrgb(ipixcolor(k))
            status = rectangle($GFILLINTERIOR,x,y,x,y)
            !!$omp end critical
            !l = -l
            !call ang2pix_ring(nside, b, l, k)
            !status = setcolorrgb(ipixcolor(k))
            !call Aitoff(b,l,x,y)
            !status = rectangle($GFILLINTERIOR,x,y,x,y)
            !b = -b
            !call ang2pix_ring(nside, b, l, k)
            !status = setcolorrgb(ipixcolor(k))
            !call Aitoff(b,l,x,y)
            !status = rectangle($GFILLINTERIOR,x,y,x,y)
            !l = -l
            !call ang2pix_ring(nside, b, l, k)
            !status = setcolorrgb(ipixcolor(k))
            !call Aitoff(b,l,x,y)
            !status = rectangle($GFILLINTERIOR,x,y,x,y)
        enddo
    enddo
    !!$omp end do
    !!$omp end parallel
    
    end subroutine DrawHealpix
    
    subroutine DrawLegend(name,meddist,pos)
    
    character(32) :: name
    character(8) :: legnum
    real(8), dimension(:) :: meddist
    integer :: i
    logical(4) :: status
    TYPE (windowconfig) wc
    real(8) :: colorgrad, valueexp
    TYPE (xycoord) pos
    
    valueexp = 1.0
    
    status = getwindowconfig(wc)
    
    do i = 0, 255
        colorgrad = i/255.0
        status = setcolorrgb(#000000)
        if((i/64.0-i/64).eq.0)then
            status = rectangle($GFILLINTERIOR,2*i+20,wc%numypixels-125,2*i+20,wc%numypixels-100)
        endif
        status = setcolorrgb(GetColorGrad(colorgrad,valueexp))
        status = rectangle($GFILLINTERIOR,2*i+20,wc%numypixels-100,2*(i+1)+20,wc%numypixels)
    enddo
    status = setcolorrgb(#000000)
    call moveto(10,wc%numypixels-215,pos)
    call outgtext(trim(name))
    call moveto(10,wc%numypixels-170,pos)
    call outgtext("0")
    call moveto(2*64,wc%numypixels-170,pos)
    write(legnum,"(F3.1)") maxval(meddist)/4000.0
    call outgtext(legnum)
    call moveto(2*128,wc%numypixels-170,pos)
    write(legnum,"(F3.1)") maxval(meddist)/2000.0
    call outgtext(legnum)
    call moveto(2*192,wc%numypixels-170,pos)
    write(legnum,"(F3.1)") maxval(meddist)/4000.0*3
    call outgtext(legnum)
    call moveto(2*i,wc%numypixels-170,pos)
    write(legnum,"(F3.1)") maxval(meddist)/1000.0
    call outgtext(legnum)
    status = rectangle($GFILLINTERIOR,2*i+20,wc%numypixels-125,2*i+20,wc%numypixels-100)
    
    end subroutine DrawLegend
    
end module ColorPalette