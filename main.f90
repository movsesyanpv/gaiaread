    use CustomDataTypes
    use GALACTIC
    use HealPix
    use LSQ
    use OGOROD
    
    implicit none
    
    type(TgasEntry), allocatable, dimension(:) :: GaiaData
    type(MedTable), allocatable, dimension(:,:) :: medall
    real(8), allocatable, dimension(:,:) :: a,r
    real(8), allocatable, dimension(:) :: ym,w,dv,v
    integer, allocatable, dimension(:,:) :: nall
    character(766) :: readStr
    integer :: line_num, i, j, k
    integer, allocatable, dimension(:) :: comp_count
    character(9) :: file_num = ""
    character(256) :: csv_file_name = ""
    real(8) :: dist, maxdist,cond,s
    real(8), dimension(400) :: med_err 
    integer, dimension(400) :: n
    
    integer, allocatable, dimension(:) :: ipixcolor, nhealp
    real(8), allocatable, dimension(:) :: meddist,nhealpr
    integer :: ipix
    integer :: nside = 12 
    integer :: Nh
    integer :: m = 11
    real(8) :: k_mod = 4.738
    character(2) :: mchar
    character(8) :: legnum = ""
    real(8) :: RAD_TO_DEG = 180/4.0/atan(1.0)
    real(8) :: DEG_TO_RAD = 4.0*atan(1.0)/180
    
    real(8) :: l,b

    real(8) :: mua, mud, d, ra
    
    integer :: mindst = 900
    integer :: maxdst = 1000
    character(6) :: mindch
    character(6) :: maxdch
    
    type(tuvw) :: Vs
    type(trot) :: Ws
    type(tdef) :: Ms
    
    real(8) :: phi, psi
    
    open(20,file='D:\gaiaread\testout.csv')
    write(20,*)"#solution_id,source_id,random_index,ref_epoch,ra,ra_error,dec,dec_error&
    ,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,ra_dec_corr,ra_parallax_corr&
    ,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr&
    ,parallax_pmdec_corr,pmra_pmdec_corr,astrometric_n_obs_al,astrometric_n_obs_ac&
    ,astrometric_n_good_obs_al,astrometric_n_good_obs_ac,astrometric_n_bad_obs_al,astrometric_n_bad_obs_ac&
    ,astrometric_delta_q,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_primary_flag&
    ,astrometric_relegation_factor,astrometric_weight_al,astrometric_weight_ac,astrometric_priors_used&
    ,matched_observations,duplicated_source,scan_direction_strength_k1,scan_direction_strength_k2&
    ,scan_direction_strength_k3,scan_direction_strength_k4,scan_direction_mean_k1,scan_direction_mean_k2&
    ,scan_direction_mean_k3,scan_direction_mean_k4,phot_g_n_obs,phot_g_mean_flux,phot_g_mean_flux_error&
    ,phot_g_mean_mag,phot_variable_flag,l,b,ecl_lon,ecl_lat"
    
    Nh = 12*nside**2
    
    allocate(comp_count(0:15), medall(0:Nh-1,9),meddist(0:Nh-1), nall(0:Nh-1,9) ,nhealp(0:Nh-1))
    allocate(nhealpr(0:Nh-1),a(2*Nh,m),ym(2*Nh),w(2*Nh),v(m),dv(m),r(m,m))
    
    comp_count = 0
    
    n = 0
    nhealp = 0
    meddist = 0
    
    do j = 0, 15
    
      line_num = 0
    
      write(file_num, "(I9.9)") j
    
      csv_file_name = 'D:\gaiaread\TgasSource_'//file_num(1:3)//'-'//file_num(4:6)//'-'//file_num(7:9)//'.csv'
      call csv_file_line_count(csv_file_name,line_num)
    
      allocate(GaiaData(line_num-1))
        
      open(10,file=csv_file_name)
      read(10,*) readStr
      do i = 1, size(GaiaData)
        read(10,*) GaiaData(i)
      enddo
      close(10)
    
      maxdist = 0
      
      nall = 0
      
      !!$omp parallel private(l,b,x,y)
      !!$omp do
      do i = 1, size(GaiaData)
        if((1.0/GaiaData(i)%parallax .ge. mindst/1000.0).and.((1.0/GaiaData(i)%parallax) .le. maxdst/1000.0))then!.and.(GaiaData(i)%pmra.lt.1290000000)) then
            dist = 1000.0/GaiaData(i)%parallax
            l = gaiadata(i)%l*DEG_TO_RAD
            b = gaiadata(i)%b*DEG_TO_RAD
            call ang2pix_ring(nside, l, b, ipix)
            meddist(ipix) = (meddist(ipix) * nhealp(ipix) + dist)/(nhealp(ipix)+1)
            maxdist = max(dist,maxdist)
            medall(ipix,1)%ihealp = ipix
            call pix2ang_ring(nside, ipix, medall(ipix,1)%lcen, medall(ipix,1)%bcen)
            if(GaiaData(i)%l .lt. 0)then
                GaiaData(i)%l = GaiaData(i)%l + 360.0
            endif
            medall(ipix,1)%l = (medall(ipix,1)%l+GaiaData(i)%l)
            medall(ipix,1)%b = (medall(ipix,1)%b+GaiaData(i)%b)
            medall(ipix,1)%ra = (medall(ipix,1)%ra+GaiaData(i)%ra)
            medall(ipix,1)%dec = (medall(ipix,1)%dec+GaiaData(i)%dec)
            medall(ipix,1)%pml = (medall(ipix,1)%pml+GaiaData(i)%pmra)
            medall(ipix,1)%pmb = (medall(ipix,1)%pmb+GaiaData(i)%pmdec)
            medall(ipix,1)%parallax = (medall(ipix,1)%parallax+GaiaData(i)%parallax)
            medall(ipix,1)%nhealp = medall(ipix,1)%nhealp + 1
            comp_count(j) = comp_count(j) + 1
        endif
      enddo
      !!$omp end do
      !!$omp end parallel
    
      deallocate(GaiaData)

      print*, trim(csv_file_name), ' contains', comp_count(j), ' compatible entries, max distance is', maxdist, 'pc'
    enddo
    
    do i=0,Nh
        medall(i,1)%lcen     = medall(i,1)%lcen     * RAD_TO_DEG
        medall(i,1)%bcen     = medall(i,1)%bcen     * RAD_TO_DEG
        medall(i,1)%l        = medall(i,1)%l        / medall(i,1)%nhealp
        medall(i,1)%b        = medall(i,1)%b        / medall(i,1)%nhealp
        medall(i,1)%ra       = medall(i,1)%ra       / medall(i,1)%nhealp
        medall(i,1)%dec      = medall(i,1)%dec      / medall(i,1)%nhealp
        medall(i,1)%pml      = medall(i,1)%pml      / medall(i,1)%nhealp
        medall(i,1)%pmb      = medall(i,1)%pmb      / medall(i,1)%nhealp
        medall(i,1)%parallax = medall(i,1)%parallax / medall(i,1)%nhealp
    enddo
        
    
    if(any((medall(:,1)%nhealp.eq.0)))then
        print*, "There is an empty HEALPIX sector at ", minloc(abs(medall(:,1)%nhealp))  
        stop
    endif

    open(30,file="D:\gaiaread\out.csv",CARRIAGECONTROL='NONE',ENCODING='UTF-8')
    write(30,*) "ihealp,lcen,bcen,l,b,pml,pmb,parallax,nhealp",new_line(" ")
    do i = 0,Nh-1
        mua = medall(i,1)%pml
        mud = medall(i,1)%pmb
        call Equatorial(medall(i,1)%lcen,medall(i,1)%bcen,ra,d)
        call GalaxMu(mua,mud,medall(i,1)%lcen,medall(i,1)%bcen,d,medall(i,1)%pml,medall(i,1)%pmb)
        write(30,*) medall(i,1)%ihealp,",",medall(i,1)%lcen,",",medall(i,1)%bcen,",",medall(i,1)%l,",",medall(i,1)%b,",",&
                    medall(i,1)%pml,",",medall(i,1)%pmb,",",medall(i,1)%parallax,",",medall(i,1)%nhealp,new_line(" ")
    enddo
    close(30)
    
    do i = 0,Nh-1
        do j = 1,m
            a(i+1,j) = kmul_base_prim(j,medall(i,1)%lcen,medall(i,1)%bcen,medall(i,1)%parallax)
            a(i+1+12*nside**2,j) = kmub_base_prim(j,medall(i,1)%lcen,medall(i,1)%bcen,medall(i,1)%parallax)
        enddo
        ym(i+1) = medall(i,1)%pml*k_mod
        ym(i+1+Nh) = medall(i,1)%pmb*k_mod
        !ym(i+1) = kmul(Vs,Ws,Ms,medall(i,1)%lcen,medall(i,1)%bcen,medall(i,1)%parallax)        !Генерация по модели Огородникова-Милна
        !ym(i+1+Nh) = kmub(Vs,Ws,Ms,medall(i,1)%lcen,medall(i,1)%bcen,medall(i,1)%parallax)
    enddo
    
    open(20, file='D:\gaiaread\matrix.dat',CARRIAGECONTROL='NONE')
    do i = 1, 2*Nh
        do j=1,m
            write(20,*) a(i,j),","
        enddo
        write(20,*)ym(i),new_line(" ")
    enddo
    close(20)

    v(:)= 0
    dv(:) = 0
    w = 1
    
    call LSQM(a,ym,w,v,dv,s,r,cond)
    
    write(mindch,"(I6)") mindst
    write(maxdch,"(I6)") maxdst
    
    write(mchar,"(i2)") m
    open(20,file='D:\gaiaread\corrMatrix.dat',carriagecontrol='none')
    do i=1,m
        write(20,"("//trim(adjustl(mchar))//"f7.3,a1)") r(i,:),new_line("")
    enddo
    close(20)
    
    open(10,file='D:\gaiaread\og'//trim(adjustl(mindch))//'-'//trim(adjustl(maxdch))//'.dat')
    print*, 'total compatible entries', sum(comp_count)
    do i = 1, m
        write(10,*) v(i), dv(i)
    enddo
    close(10)

contains

    subroutine csv_file_line_count ( csv_file_name, line_num )

          character ( len = * ) csv_file_name
          integer ( kind = 4 ) ierror
          integer ( kind = 4 ) input_status
          integer ( kind = 4 ) input_unit
          character ( len = 1023 ) line
          integer ( kind = 4 ) line_num

          line_num = -1

          call get_unit ( input_unit )

          open ( unit = input_unit, file = csv_file_name, status = 'old', &
            iostat = input_status )

          if ( input_status /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'CSV_FILE_LINE_COUNT - Fatal error!'
            write ( *, '(a,i8)' ) '  Could not open "' // trim ( csv_file_name ) // '".'
            stop
          end if

          line_num = 0

          do

            read ( input_unit, '(a)', iostat = input_status ) line

            if ( input_status /= 0 ) then
              ierror = line_num
              exit
            end if

            line_num = line_num + 1

          end do

          close ( unit = input_unit )

          return

        end subroutine csv_file_line_count

    subroutine get_unit ( iunit )

          integer ( kind = 4 ) i
          integer ( kind = 4 ) ios
          integer ( kind = 4 ) iunit
          logical lopen

          iunit = 0

          do i = 1, 99

            if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

              inquire ( unit = i, opened = lopen, iostat = ios )

              if ( ios == 0 ) then
                if ( .not. lopen ) then
                  iunit = i
                  return
                end if
              end if

            end if

          end do

          return

    end

end