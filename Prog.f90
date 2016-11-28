    program GaiaRead

    use CustomDataTypes
    
    implicit none
    
    type(TgasEntry), allocatable, dimension(:) :: GaiaData
    character(80) :: readStr
    integer :: line_num, i, j, k
    integer(16), allocatable, dimension(:) :: comp_count
    character(9) :: file_num = ""
    character(26) :: csv_file_name = ""
    real(8) :: dist, maxdist
    real(8), dimension(400) :: med_err 
    integer, dimension(400) :: n


    open(20,file='testout.csv')
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

    allocate(comp_count(0:15))

    comp_count = 0

    med_err = 0
    n = 0

    do j = 0, 15

      line_num = 0

      write(file_num, "(I9.9)") j

      csv_file_name = 'TgasSource_'//file_num(1:3)//'-'//file_num(4:6)//'-'//file_num(7:9)//'.csv'
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

      do i = 1, size(GaiaData)
        if((GaiaData(i)%parallax .ge. 0).and.((aint(1000.0/GaiaData(i)%parallax/100)+1).le.400)) then
          dist = 1000.0/GaiaData(i)%parallax
          maxdist = max(dist,maxdist)
          k = int(aint(dist/100))+1
          !print*, k
          med_err(k) = (med_err(k) * n(k) + GaiaData(i)%parallax_error/GaiaData(i)%parallax)/(n(k)+1)
          n(k) = n(k) + 1
          !write(20,*) GaiaData(i)%parallax, GaiaData(i)%parallax_error/GaiaData(i)%parallax
          !write(20,*) log10(dist), GaiaData(i)%parallax_error*dist/1000.0
          comp_count(j) = comp_count(j) + 1
        else
          !print*,GaiaData(i)
        endif
      enddo

      deallocate(GaiaData)
      
      !print*, csv_file_name, ' contains', comp_count(j), ' compatible entries, max distance is', maxdist, 'pc'
      print*, 'max distance in ', csv_file_name, ' is', maxdist, 'pc'

    enddo

    print*, 'total compatible entries', sum(comp_count)

    do i = 1, 400
      write(20,*) i*100, med_err(i)
    enddo

    close(20)

    line_num = 0
    call csv_file_line_count('testout.csv',line_num)
    print*,'output lines',line_num

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

