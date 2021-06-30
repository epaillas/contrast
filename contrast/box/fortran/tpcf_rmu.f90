program tpcf_rmu
  use OMP_LIB
  use procedures
  implicit none
  
  real*8 :: rgrid, boxsize, vol, mean_density
  real*8 :: disx, disy, disz, dis, dis2, mu
  real*8 :: comx, comy, comz
  real*8 :: rwidth, irwidth, dim1_max, dim1_min, dim1_min2, dim1_max2
  real*8 :: dim2_max, dim2_min, muwidth, imuwidth
  real*8 :: pi = 4.*atan(1.)

  integer*8 :: ndata2, ndata1, dim1_nbin, dim2_nbin, rind, muind
  integer*8 :: i, ii, ix, iy, iz, ix2, iy2, iz2
  integer*8 :: ipx, ipy, ipz, ndif
  integer*8 :: ngrid, end, beginning, rate
  integer*4 :: nthreads
  
  integer*8, dimension(:, :, :), allocatable :: lirst
  integer*8, dimension(:), allocatable :: ll
  
  real*8, allocatable, dimension(:,:)  :: data1, data2
  real*8, dimension(:), allocatable :: weight1, weight2
  real*8, dimension(:, :), allocatable :: D1D2, xi_rmu
  real*8, dimension(:), allocatable :: rbin, rbin_edges, mubin, mubin_edges

  logical :: has_velocity
  logical :: debug = .true.
  
  character(20), external :: str
  character(len=500) :: data_filename1, data_filename2, output_filename, nthreads_char
  character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char, ngrid_char, box_char
  character(len=10) :: dim2_max_char, dim2_min_char, dim2_nbin_char
  
  if (debug) then
    if (iargc() .lt. 9) then
        write(*,*) 'Some arguments are missing.'
        write(*,*) '1) data_filename'
        write(*,*) '2) data_filename2'
        write(*,*) '3) output_filename'
        write(*,*) '4) boxsize'
        write(*,*) '5) dim1_min'
        write(*,*) '6) dim1_max'
        write(*,*) '7) dim1_nbin'
        write(*,*) '8) dim2_min'
        write(*,*) '9) dim2_max'
        write(*,*) '10) dim2_nbin'
        write(*,*) '11) ngrid'
        write(*,*) '12) nthreads'
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
  call get_command_argument(number=8, value=dim2_min_char)
  call get_command_argument(number=9, value=dim2_max_char)
  call get_command_argument(number=10, value=dim2_nbin_char)
  call get_command_argument(number=11, value=ngrid_char)
  call get_command_argument(number=12, value=nthreads_char)
  
  read(box_char, *) boxsize
  read(dim1_min_char, *) dim1_min
  read(dim1_max_char, *) dim1_max
  read(dim1_nbin_char, *) dim1_nbin
  read(dim2_min_char, *) dim2_min
  read(dim2_max_char, *) dim2_max
  read(dim2_nbin_char, *) dim2_nbin
  read(ngrid_char, *) ngrid
  read(nthreads_char, *) nthreads

  if (debug) then
    write(*,*) '-----------------------'
    write(*,*) 'Running tpcf_rmu.exe'
    write(*,*) 'input parameters:'
    write(*,*) ''
    write(*, *) 'data_filename: ', trim(data_filename1)
    write(*, *) 'data_filename_2: ', trim(data_filename2)
    write(*, *) 'boxsize: ', trim(box_char)
    write(*, *) 'output_filename: ', trim(output_filename)
    write(*, *) 'dim1_min: ', trim(dim1_min_char), ' Mpc'
    write(*, *) 'dim1_max: ', trim(dim1_max_char), ' Mpc'
    write(*, *) 'dim1_nbin: ', trim(dim1_nbin_char)
    write(*, *) 'dim2_min: ', trim(dim2_min_char)
    write(*, *) 'dim2_max: ', trim(dim2_max_char)
    write(*, *) 'dim2_nbin: ', trim(dim2_nbin_char)
    write(*, *) 'ngrid: ', trim(ngrid_char)
    write(*, *) 'nthreads: ', trim(nthreads_char)
    write(*,*) ''
  end if

  call read_unformatted(data_filename1, data1, weight1, ndata1, has_velocity)
  call read_unformatted(data_filename2, data2, weight2, ndata2, has_velocity)
  call linked_list(data2, boxsize, ngrid, ll, lirst, rgrid)
  call binning(dim1_min, dim1_max, dim1_nbin, rbin, rbin_edges, rwidth)
  call binning(dim2_min, dim2_max, dim2_nbin, mubin, mubin_edges, muwidth)

  allocate(D1D2(dim1_nbin, dim2_nbin))
  allocate(xi_rmu(dim1_nbin, dim2_nbin))
  D1D2 = 0
  dim1_min2 = dim1_min ** 2
  dim1_max2 = dim1_max ** 2
  ndif = int(dim1_max / rgrid + 1.)
  irwidth = 1 / rwidth
  imuwidth = 1 / muwidth
  comx = 0
  comy = 0
  comz = 1 ! assume LOS along the z axis
  mean_density = ndata2 / boxsize ** 3

  call OMP_SET_NUM_THREADS(nthreads)
  if (debug) then
    write(*, *) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()
  end if
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, ipz, &
  !$OMP ix, iy, iz, ix2, iy2, iz2, disx, disy, disz, dis, dis2, rind, &
  !$OMP mu, muind) REDUCTION(+:D1D2)
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
  
          ii = lirst(ix2,iy2,iz2)
          if(ii.ne.0) then
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

              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                dis = sqrt(dis2)

                mu = (disx * comx + disy * comy + disz * comz) &
                & / (dis * sqrt(comx * comx + comy * comy + comz * comz))

                if (mu .eq. 1) cycle ! might have to deal with this later

                rind = int((dis - dim1_min) * irwidth + 1)
                muind = int((mu - dim2_min) * imuwidth + 1)

                D1D2(rind, muind) = D1D2(rind, muind) + weight1(i) * weight2(ii)
              end if
  
              if(ii.eq.lirst(ix2,iy2,iz2)) exit
  
            end do
          end if
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  do i = 1, dim1_nbin
    do ii = 1, dim2_nbin
      vol = 4./3 * pi * (rbin_edges(i + 1) ** 3 &
        & - rbin_edges(i) ** 3) / (dim2_nbin)
      xi_rmu(i, ii) = D1D2(i, ii)  / (vol * mean_density * ndata1) - 1 
    end do
  end do
  
  write(*,*) ''
  write(*,*) 'Calculation finished. Writing output...'
  
  open(12, file=output_filename, status='replace')
  do i = 1, dim1_nbin
    do ii = 1, dim2_nbin
      write(12, fmt='(4E15.5)') rbin(i), mubin(ii), xi_rmu(i, ii), D1D2(i, ii)
    end do
  end do

  call system_clock(end)
  if (debug) then
    print *, "elapsed time: ", real(end - beginning) / real(rate)
  end if

  end program tpcf_rmu
  
