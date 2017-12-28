    use CustomDataTypes
    !use Projection
    use GALACTIC
    use HealPix
    !use ColorPalette
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
    
    real(8) :: colorgrad, valueexp
    integer, allocatable, dimension(:) :: ipixcolor, nhealp
    real(8), allocatable, dimension(:) :: meddist,nhealpr
    integer :: ipix
    integer :: nside = 10
    character(8) :: legnum = ""
    
    real(8) :: x,y,l,b
    
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
    
    allocate(comp_count(0:15), medall(0:12*nside**2-1,9),meddist(0:12*nside**2-1), nall(0:12*nside**2-1,9) ,nhealp(0:12*nside**2-1))
    allocate(nhealpr(0:12*nside**2),a(24*nside**2,3),ym(24*nside**2),w(24*nside**2),v(3),dv(3),r(3,3))
    
    comp_count = 0
    
    !med_err = 0
    n = 0
    nhealp = 0
    meddist = 0
    
     do j = 0, 15
    
      line_num = 0
    
      write(file_num, "(I9.9)") j
    
      csv_file_name = 'D:\gaiaread\TgasSource_'//file_num(1:3)//'-'//file_num(4:6)//'-'//file_num(7:9)//'.csv'
      call csv_file_line_count(csv_file_name,line_num)
    
      allocate(GaiaData(line_num-1))
    
      !print*,csv_file_name, ' contains ', line_num-1, ' entries'
    
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
        if((GaiaData(i)%parallax .ge. 0).and.((1.0/GaiaData(i)%parallax) .le. 0.1))then!.and.(GaiaData(i)%pmra.lt.1290000000)) then
            dist = 1000.0/GaiaData(i)%parallax
            !call Galaxy(GaiaData(i)%ra,GaiaData(i)%dec,l,b)
            l = gaiadata(i)%ra/180*4*atan(1.0)
            b = gaiadata(i)%dec/180*4*atan(1.0)
            call ang2pix_ring(nside, l, b, ipix)
            meddist(ipix) = (meddist(ipix) * nhealp(ipix) + dist)/(nhealp(ipix)+1)
          maxdist = max(dist,maxdist)
          k = int(aint(dist/100))+1
          !print*, k
          medall(ipix,1)%ihealp = ipix
          call pix2ang_ring(nside, ipix, medall(ipix,1)%lcen, medall(ipix,1)%bcen)
          medall(ipix,1)%l = (medall(ipix,1)%l*nhealp(ipix)+GaiaData(i)%ra)/(nhealp(ipix)+1)
          medall(ipix,1)%b = (medall(ipix,1)%b*nhealp(ipix)+GaiaData(i)%dec)/(nhealp(ipix)+1)
          medall(ipix,1)%pml = (medall(ipix,1)%pml*nhealp(ipix)+GaiaData(i)%pmra)/(nhealp(ipix)+1)
          medall(ipix,1)%pmb = (medall(ipix,1)%pmb*nhealp(ipix)+GaiaData(i)%pmdec)/(nhealp(ipix)+1)
          medall(ipix,1)%parallax = (medall(ipix,1)%parallax*nhealp(ipix)+GaiaData(i)%parallax)/(nhealp(ipix)+1)
          medall(ipix,1)%nhealp = medall(ipix,1)%nhealp + 1
          !medall(ipix,k)%ra = (medall(ipix,k)%ra*nall(ipix,k)+GaiaData(i)%ra)/(nall(ipix,k)+1)
          !medall(ipix,k)%dec = (medall(ipix,k)%dec*nall(ipix,k)+GaiaData(i)%dec)/(nall(ipix,k)+1)
          !medall(ipix,k)%pmra = (medall(ipix,k)%pmra*nall(ipix,k)+GaiaData(i)%pmra)/(nall(ipix,k)+1)
          !medall(ipix,k)%pmdec = (medall(ipix,k)%pmdec*nall(ipix,k)+GaiaData(i)%pmdec)/(nall(ipix,k)+1)
          !medall(ipix,k)%parallax = (medall(ipix,k)%parallax*nall(ipix,k)+GaiaData(i)%parallax)/(nall(ipix,k)+1)
          !nall(ipix,k) = nall(ipix,k) + 1
          med_err(k) = (med_err(k) * n(k) + GaiaData(i)%parallax_error/GaiaData(i)%parallax)/(n(k)+1)
          n(k) = n(k) + 1
          !write(20,*) GaiaData(i)%parallax, GaiaData(i)%parallax_error/GaiaData(i)%parallax
          !write(20,*) log10(dist), GaiaData(i)%parallax_error*dist/1000.0
          comp_count(j) = comp_count(j) + 1
        !else
        !  !print*,GaiaData(i)
        endif
      enddo
      !!$omp end do
      !!$omp end parallel
    
      deallocate(GaiaData)
    !  
      print*, trim(csv_file_name), ' contains', comp_count(j), ' compatible entries, max distance is', maxdist, 'pc'
    !  print*, 'max distance in ', trim(csv_file_name), ' is', maxdist, 'pc'
    !
    enddo
    
    open(30,file="D:\gaiaread\out.csv",CARRIAGECONTROL='NONE',ENCODING='UTF-8')
    write(30,*) "ihealp,lcen,bcen,ra,dec,pmra,pmdec,parallax,nhealp",new_line(" ")
    do i = 0,12*nside**2-1
        write(30,*) medall(i,1)%ihealp,",",medall(i,1)%lcen*180/4.0/atan(1.0),",",medall(i,1)%bcen*180/4.0/atan(1.0),",",medall(i,1)%l,",",medall(i,1)%b,",",medall(i,1)%pml,",",medall(i,1)%pmb,",",medall(i,1)%parallax,",",medall(i,1)%nhealp,new_line(" ")
    enddo
    close(30)
    
    do i = 0, 12*nside**2-1
        a(i+1,1) = kmul_base(1,medall(i,1)%l,medall(i,1)%b,medall(i,1)%parallax)           !TODO: перевести mura mudec в mul mub; составить систему и решить ее при помощи МНК
        a(i+1,2) = kmul_base(2,medall(i,1)%l,medall(i,1)%b,medall(i,1)%parallax)
        a(i+1,3) = 0!kmul_base(3,medall(i,1)%ra,medall(i,1)%dec,medall(i,1)%parallax)
        ym(i+1) = medall(i,1)%pml/1000.0/3600.0 * cosd(medall(i,1)%b)
        a(i+1+12*nside**2,1) = kmub_base(1,medall(i,1)%l,medall(i,1)%b,medall(i,1)%parallax)
        a(i+1+12*nside**2,2) = kmub_base(2,medall(i,1)%l,medall(i,1)%b,medall(i,1)%parallax)
        a(i+1+12*nside**2,3) = kmub_base(3,medall(i,1)%l,medall(i,1)%b,medall(i,1)%parallax)
        ym(i+1+12*nside**2) = medall(i,1)%pmb/1000.0/3600.0
    enddo
    
    open(20, file='D:\gaiaread\matrix.dat',CARRIAGECONTROL='NONE')
    do i = 1, 2*12*nside**2
        write(20,*) a(i,1),"v1+", a(i,2),"v2+", a(i,3),"v3=", ym(i),new_line(" ")
    enddo
    close(20)

    v = 0
    dv = 0
    
    call LSQM(a,ym,w,v,dv,s,r,cond)
    
    print*, 'total compatible entries', sum(comp_count)
    do i = 1, 3
        print*, v(i), dv(i)
    enddo

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