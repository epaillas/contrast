program tpcf
    use OMP_LIB
    use procedures
    implicit none
    
    real*8 :: rgrid, vol, mean_density
    real*8 :: boxsize
    real*8 :: disx, disy, disz, dis, dis2
    real*8 :: rwidth, dim1_max, dim1_min, dim1_max2, dim1_min2
    real*8 :: pi = 4.*atan(1.)
    
    integer*8 :: ndata1, ndata2, dim1_nbin, rind
    integer*8 :: i, ii, ix, iy, iz, ix2, iy2, iz2
    integer*8 :: ipx, ipy, ipz, ndif
    integer*8 :: ngrid
    integer*8 :: end, beginning, rate
    integer*4 :: nthreads, use_weights
    
    integer*8, dimension(:, :, :), allocatable :: lirst
    integer*8, dimension(:), allocatable :: ll
    
    real*8, allocatable, dimension(:,:)  :: data1, data2
    real*8, dimension(:), allocatable :: D1D2, weight1, weight2, xi_r
    real*8, dimension(:), allocatable :: rbin, rbin_edges
  
    character(20), external :: str
    character(len=500) :: data_filename1, data_filename2, output_filename, nthreads_char
    character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char, ngrid_char, box_char
    character(len=10) :: use_weights_char, data_fileformat

    logical :: debug = .true.
    
    if (debug) then
      if (iargc() .lt. 11) then
          write(*,*) 'Some arguments are missing.'
          write(*,*) '1) data_filename'
          write(*,*) '2) data_filename_2'
          write(*,*) '3) output_filename'
          write(*,*) '4) boxsize'
          write(*,*) '5) dim1_min'
          write(*,*) '6) dim1_max'
          write(*,*) '7) dim1_nbin'
          write(*,*) '8) ngrid'
          write(*,*) '9) nthreads'
          write(*,*) '10) use_weights'
          write(*,*) '11) data_fileformat'
          write(*,*) ''
          stop
        end if
      end if

    call system_clock(beginning, rate)
      
    call get_command_argument(number=1, value=data_filename1)
    call get_command_argument(number=2, value=data_filename2)
    call get_command_argument(number=3, value=output_filename)
    call get_command_argument(number=4, value=box_char)
    call get_command_argument(number=5, value=dim1_min_char)
    call get_command_argument(number=6, value=dim1_max_char)
    call get_command_argument(number=7, value=dim1_nbin_char)
    call get_command_argument(number=8, value=ngrid_char)
    call get_command_argument(number=9, value=nthreads_char)
    call get_command_argument(number=10, value=use_weights_format)
    call get_command_argument(number=11, value=data_fileformat)
    
    read(box_char, *) boxsize
    read(dim1_min_char, *) dim1_min
    read(dim1_max_char, *) dim1_max
    read(dim1_nbin_char, *) dim1_nbin
    read(ngrid_char, *) ngrid
    read(nthreads_char, *) nthreads
    read(use_weights_char, *) use_weights

    if (debug) then
      write(*,*) '-----------------------'
      write(*,*) 'Running tpcf.exe'
      write(*,*) 'input parameters:'
      write(*,*) ''
      write(*, *) 'data_filename: ', trim(data_filename1)
      write(*, *) 'data_filename_2: ', trim(data_filename2)
      write(*, *) 'boxsize: ', trim(box_char)
      write(*, *) 'output_filename: ', trim(output_filename)
      write(*, *) 'dim1_min: ', trim(dim1_min_char), ' Mpc'
      write(*, *) 'dim1_max: ', trim(dim1_max_char), ' Mpc'
      write(*, *) 'dim1_nbin: ', trim(dim1_nbin_char)
      write(*, *) 'ngrid: ', trim(ngrid_char)
      write(*, *) 'nthreads: ', trim(nthreads_char)
      write(*, *) 'use_weights: ', trim(use_weights_char)
      write(*, *) 'data_fileformat: ', trim(dat_fileformat)
      write(*,*) ''
    end if

    ! read tracers file
        if (trim(tracers_fileformat) == 'ascii') then
            if (use_weights == 1) then
                call read_catalogue_type2(data_filename1, data1, weight1, ndata1)
                call read_catalogue_type2(data_filename2, data2, weight2, ndata2)
            else
                call read_catalogue_type1(data_filename1, data1, weight1, ndata1)
                call read_catalogue_type1(data_filename2, data2, weight2, ndata2)
            end if
        else
            if (use_weights == 1) then
                call read_catalogue_type6(data_filename1, data1, weight1, ndata1)
                call read_catalogue_type6(data_filename2, data2, weight2, ndata2)
            else
                call read_catalogue_type5(data_filename1, data1, weight1, ndata1)
                call read_catalogue_type5(data_filename2, data2, weight2, ndata2)
            end if
        end if

    call linked_list(data2, boxsize, ngrid, ll, lirst, rgrid)
    call binning(dim1_min, dim1_max, dim1_nbin, rbin, rbin_edges, rwidth)

    allocate(D1D2(dim1_nbin))
    allocate(xi_r(dim1_nbin))
    mean_density = ndata2 / (boxsize **3)
    D1D2 = 0
    xi_r = 0
    dim1_min2 = dim1_min ** 2
    dim1_max2 = dim1_max ** 2
    ndif = int(dim1_max / rgrid + 1.)

    call OMP_SET_NUM_THREADS(nthreads)
    if (debug) then
      write(*, *) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()
    end if
      
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, ipz, &
    !$OMP ix, iy, iz, ix2, iy2, iz2, disx, disy, disz, dis, dis2, rind) &
    !$OMP REDUCTION(+:D1D2)
    do i = 1, ndata1
      ipx = int(data1(1, i) / rgrid + 1.)
      ipy = int(data1(2, i) / rgrid + 1.)
      ipz = int(data1(3, i) / rgrid + 1.)
  
      do ix = ipx - ndif, ipx + ndif
        do iy = ipy - ndif, ipy + ndif
          do iz = ipz - ndif, ipz + ndif
            if ((ix - ipx)**2 + (iy - ipy)**2 + (iz - ipz)**2 .gt. (ndif + 1)**2) cycle
    
            ix2 = ix
            iy2 = iy
            iz2 = iz
    
            if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
            if (ix2 .lt. 1) ix2 = ix2 + ngrid
            if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
            if (iy2 .lt. 1) iy2 = iy2 + ngrid
            if (iz2 .gt. ngrid) iz2 = iz2 - ngrid
            if (iz2 .lt. 1) iz2 = iz2 + ngrid
    
            ii = lirst(ix2, iy2, iz2)
            if(ii .ne. 0) then
              do
                ii = ll(ii)
                disx = data2(1, ii) - data1(1, i)
                disy = data2(2, ii) - data1(2, i)
                disz = data2(3, ii) - data1(3, i)
  
                if (disx .lt. -boxsize/2) disx = disx + boxsize
                if (disx .gt. boxsize/2) disx = disx - boxsize
                if (disy .lt. -boxsize/2) disy = disy + boxsize
                if (disy .gt. boxsize/2) disy = disy - boxsize
                if (disz .lt. -boxsize/2) disz = disz + boxsize
                if (disz .gt. boxsize/2) disz = disz - boxsize

                dis2 = disx ** 2 + disy ** 2 + disz ** 2

                if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                  dis = sqrt(dis2)
                  rind = int((dis - dim1_min) / rwidth + 1)
                  D1D2(rind) = D1D2(rind) + weight1(i) * weight2(ii)
                end if
    
                if (ii .eq. lirst(ix2, iy2, iz2)) exit
    
              end do
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  
    do i = 1, dim1_nbin
      vol = 4./3 * pi * (rbin_edges(i + 1) ** 3 - rbin_edges(i) ** 3)
      xi_r(i) = D1D2(i) / (vol * mean_density * ndata1) - 1
    end do
    
    if (debug) then
      write(*,*) ''
      write(*,*) 'Calculation finished. Writing output...'
    end if
    
    open(12, file=output_filename, status='replace')
    do i = 1, dim1_nbin
      write(12, fmt='(3E15.5)') rbin(i), xi_r(i), D1D2(i)
    end do
  
    call system_clock(end)
    if (debug) then
      print *, "elapsed time: ", real(end - beginning) / real(rate)
    end if

    end program tpcf
    
