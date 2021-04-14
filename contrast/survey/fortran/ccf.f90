module procedures
  implicit none
contains

  subroutine linked_list(pos, ngrid, gridmin, gridmax, ll, lirst, rgrid)
    implicit none
    integer*8 :: i, ng, ipx, ipy, ipz
    integer*8, intent(in) :: ngrid
    real*8, intent(out) :: rgrid
    real*8, intent(in) :: gridmin, gridmax
    real*8, dimension(:,:), intent(in) :: pos
    integer*8, dimension(:,:,:), intent(out) :: lirst
    integer*8, dimension(:), intent(out) :: ll

    rgrid = (gridmax - gridmin) / ngrid

    ng = size(pos, dim=2)
    lirst = 0
    ll = 0
    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)
      ipz = int((pos(3, i) - gridmin) / rgrid + 1.)
      if(ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.&
      ipz.gt.0.and.ipz.le.ngrid) lirst(ipx, ipy, ipz) = i

    end do

    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)
      ipz = int((pos(3, i) - gridmin) / rgrid + 1.)

      if (ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.ipz&
      &.gt.0.and.ipz.le.ngrid) then
        ll(lirst(ipx, ipy, ipz)) = i
        lirst(ipx, ipy, ipz) = i
      endif
    end do

  end subroutine linked_list

end module procedures



program tpcf
    use procedures
    use OMP_LIB
    implicit none
    
    real*8 :: rgrid, disx, disy, disz, dis, dis2, gridmin, gridmax
    real*8 :: rwidth, dim1_max, dim1_min, dim1_max2, dim1_min2, norm
    
    integer*8 :: ng, nc, nr, dim1_nbin, rind
    integer*8 :: i, ii, ix, iy, iz
    integer*8 :: nrows, ncols
    integer*8 :: ngrid, ipx, ipy, ipz, ndif
    integer*8 :: end, beginning, rate
    integer*4 :: nthreads
    
    integer*8, dimension(:, :, :), allocatable :: lirst_tracers, lirst_randoms
    integer*8, dimension(:), allocatable :: ll_tracers, ll_randoms
    
    real*8, allocatable, dimension(:,:)  :: tracers, centres, randoms
    real*8, dimension(:), allocatable :: DD, DR, delta
    real*8, dimension(:), allocatable :: weights_tracers, weights_centres, weights_randoms
    real*8, dimension(:), allocatable :: rbin, rbin_edges
    real*8, dimension(:, :), allocatable :: DD_i, DR_i
  
    character(20), external :: str
    character(len=500) :: data_filename, data_filename_2, output_filename, randoms_filename
    character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char, ngrid_char
    character(len=10) :: nthreads_char, gridmin_char, gridmax_char
    character(len=10) :: estimator = 'DP'

    logical :: debug = .true.
    
    if (debug) then
      if (iargc() .lt. 11) then
          write(*,*) 'Some arguments are missing.'
          write(*,*) '1) data_filename'
          write(*,*) '2) data_filename_2'
          write(*,*) '3) random_filename'
          write(*,*) '4) output_filename'
          write(*,*) '5) dim1_min'
          write(*,*) '6) dim1_max'
          write(*,*) '7) dim1_nbin'
          write(*,*) '8) ngrid'
          write(*,*) '9) gridmin'
          write(*,*) '10) gridmax'
          write(*,*) '11) nthreads'
          write(*,*) ''
          stop
        end if
      end if

    call system_clock(beginning, rate)
      
    call get_command_argument(number=1, value=data_filename)
    call get_command_argument(number=2, value=data_filename_2)
    call get_command_argument(number=3, value=randoms_filename)
    call get_command_argument(number=4, value=output_filename)
    call get_command_argument(number=5, value=dim1_min_char)
    call get_command_argument(number=6, value=dim1_max_char)
    call get_command_argument(number=7, value=dim1_nbin_char)
    call get_command_argument(number=8, value=ngrid_char)
    call get_command_argument(number=9, value=gridmin_char)
    call get_command_argument(number=10, value=gridmax_char)
    call get_command_argument(number=11, value=nthreads_char)
    
    read(dim1_min_char, *) dim1_min
    read(dim1_max_char, *) dim1_max
    read(dim1_nbin_char, *) dim1_nbin
    read(ngrid_char, *) ngrid
    read(gridmin_char, *) gridmin
    read(gridmax_char, *) gridmax
    read(nthreads_char, *) nthreads

    if (debug) then
      write(*,*) '-----------------------'
      write(*,*) 'Running tpcf.exe'
      write(*,*) 'input parameters:'
      write(*,*) ''
      write(*, *) 'data_filename: ', trim(data_filename)
      write(*, *) 'data_filename_2: ', trim(data_filename_2)
      write(*, *) 'randoms_filename: ', trim(randoms_filename)
      write(*, *) 'output_filename: ', trim(output_filename)
      write(*, *) 'dim1_min: ', trim(dim1_min_char), ' Mpc'
      write(*, *) 'dim1_max: ', trim(dim1_max_char), ' Mpc'
      write(*, *) 'dim1_nbin: ', trim(dim1_nbin_char)
      write(*, *) 'gridmin: ', trim(gridmin_char), ' Mpc'
      write(*, *) 'gridmax: ', trim(gridmax_char), ' Mpc'
      write(*, *) 'ngrid: ', trim(ngrid_char)
      write(*, *) 'nthreads: ', trim(nthreads_char)
      write(*,*) ''
    end if

    ! read tracers catalogue
    open(10, file=data_filename, status='old', form='unformatted')
    read(10) nrows
    read(10) ncols
    allocate(tracers(ncols, nrows))
    allocate(weights_tracers(nrows))
    read(10) tracers
    close(10)
    ng = nrows
    if (ncols .eq. 4) then
      weights_tracers = tracers(4, :)
      if (debug) write(*,*) 'Tracer file has weight information.'
    else
      weights_tracers = 1.0
    end if
    if (debug) then
      write(*,*) 'ntracers dim: ', size(tracers, dim=1), size(tracers, dim=2)
      write(*,*) 'tracers(min, max) = ', minval(tracers(:,:)), maxval(tracers(:,:))
      write(*,*) 'weights_tracers(min, max) = ', minval(weights_tracers), maxval(weights_tracers)
    end if

    open(11, file=data_filename_2, status='old', form='unformatted')
    read(11) nrows
    read(11) ncols
    allocate(centres(ncols, nrows))
    allocate(weights_centres(nrows))
    read(11) centres
    close(11)
    nc = nrows
    if (ncols .eq. 4) then
      weights_centres = centres(4, :)
      if (debug) write(*,*) 'Centres file has weight information.'
    else
      weights_centres = 1.0
    end if
    if (debug) then
      write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)
      write(*,*) 'centres(min), tracers(max) = ', minval(centres(:,:)), maxval(centres(:,:))
      write(*,*) 'weights_centres(min, max) = ', minval(weights_centres), maxval(weights_centres)
    end if

    ! read random catalogue
  open(11, file=randoms_filename, status='old', form='unformatted')
  read(11) nrows
  read(11) ncols
  allocate(randoms(ncols, nrows))
  allocate(weights_randoms(nrows))
  read(11) randoms
  close(11)
  nr = nrows
  if (ncols .eq. 4) then
    weights_randoms = randoms(4, :)
    if (debug) write(*,*) 'Tracer file has weight information.'
  else
    weights_randoms = 1.0
  end if
  if (debug) then 
    write(*,*) 'nrandoms dim: ', size(randoms, dim=1), size(randoms, dim=2)
    write(*,*) 'randoms(min), randoms(max) = ', minval(randoms(:,:)), maxval(randoms(:,:))
    write(*,*) 'weights_randoms(min, max) = ', minval(weights_randoms), maxval(weights_randoms)
  end if

  ! construct linked lists for tracers and randoms
  allocate(ll_tracers(ng))
  allocate(ll_randoms(nr))
  allocate(lirst_tracers(ngrid, ngrid, ngrid))
  allocate(lirst_randoms(ngrid, ngrid, ngrid))
  call linked_list(tracers, ngrid, gridmin, gridmax, ll_tracers, lirst_tracers, rgrid)
  call linked_list(randoms, ngrid, gridmin, gridmax, ll_randoms, lirst_randoms, rgrid)

  allocate(rbin(dim1_nbin))
  allocate(rbin_edges(dim1_nbin + 1))
  allocate(DD(dim1_nbin))
  allocate(DD_i(nc, dim1_nbin))
  allocate(DR(dim1_nbin))
  allocate(DR_i(nc, dim1_nbin))
  allocate(delta(dim1_nbin))
  
  rwidth = (dim1_max - dim1_min) / dim1_nbin
  do i = 1, dim1_nbin + 1
    rbin_edges(i) = dim1_min+(i-1)*rwidth
  end do
  do i = 1, dim1_nbin
    rbin(i) = rbin_edges(i+1)-rwidth/2.
  end do
  
  DD = 0
  DD_i = 0
  DR = 0
  DR_i = 0
  delta = 0
  dim1_min2 = dim1_min ** 2
  dim1_max2 = dim1_max ** 2
  ndif = int(dim1_max / rgrid + 1.)

  call OMP_SET_NUM_THREADS(nthreads)
  if (debug) then
    write(*, *) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()
  end if
    
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, ipz, &
  !$OMP ix, iy, iz, disx, disy, disz, dis, dis2, rind)
  do i = 1, nc
    ipx = int((centres(1, i) - gridmin) / rgrid + 1.)
    ipy = int((centres(2, i) - gridmin) / rgrid + 1.)
    ipz = int((centres(3, i) - gridmin) / rgrid + 1.)

    do ix = ipx - ndif, ipx + ndif, 1
      do iy = ipy - ndif, ipy + ndif, 1
        do iz = ipz - ndif, ipz + ndif, 1 
          if ((ix - ipx)**2 + (iy - ipy)**2 + (iz - ipz)**2 .gt. (ndif+ 1)**2) cycle
  
          ii = lirst_tracers(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_tracers(ii)
              disx = tracers(1, ii) - centres(1, i)
              disy = tracers(2, ii) - centres(2, i)
              disz = tracers(3, ii) - centres(3, i)

              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                dis = sqrt(dis2)
                rind = int((dis - dim1_min) / rwidth + 1)
                DD_i(i, rind) = DD_i(i, rind) + weights_centres(i) * weights_tracers(ii)
              end if
  
              if(ii .eq. lirst_tracers(ix, iy, iz)) exit
  
            end do
          end if

          ii = lirst_randoms(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_randoms(ii)
              disx = randoms(1, ii) - centres(1, i)
              disy = randoms(2, ii) - centres(2, i)
              disz = randoms(3, ii) - centres(3, i)

              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                dis = sqrt(dis2)
                rind = int((dis - dim1_min) / rwidth + 1)
                DR_i(i, rind) = DR_i(i, rind) + weights_centres(i) * weights_randoms(ii)
              end if
  
              if (ii .eq. lirst_randoms(ix, iy, iz)) exit
  
            end do
          end if
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  norm = SUM(weights_randoms) / SUM(weights_tracers)

  ! Calculate density contrast
  if (estimator .eq. 'DP') then
    delta = norm * (DD / DR) - 1
  else
    write(*,*) 'Estimator for the correlation function was not recognized.'
    stop
  end if
  
  if (debug) then
    write(*,*) ''
    write(*,*) 'Calculation finished. Writing output...'
  end if
  
  open(12, file=output_filename, status='replace')
  do i = 1, dim1_nbin
    write(12, fmt='(2f15.5)') rbin(i), delta(i)
  end do

  call system_clock(end)
  if (debug) then
    print *, "elapsed time: ", real(end - beginning) / real(rate)
  end if

  end program tpcf
    