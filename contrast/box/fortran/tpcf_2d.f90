program tpcf
  implicit none
  
  real*8 :: rgrid, area, rhomean
  real*8 :: boxsize
  real*8 :: disx, disy, dis, dis2
  real*8 :: rwidth, dim1_max, dim1_min, dim1_max2, dim1_min2
  real*8 :: pi = 4.*atan(1.)
  
  integer*8 :: ntracers, ncentres, dim1_nbin, rind
  integer*8 :: i, ii, ix, iy, ix2, iy2
  integer*8 :: indx, indy, nrows, ncols
  integer*8 :: ipx, ipy, ndif
  integer*8 :: ngrid
  integer*8 :: end, beginning, rate
  
  integer*8, dimension(:, :), allocatable :: lirst, nlirst
  integer*8, dimension(:), allocatable :: ll
  
  real*8, allocatable, dimension(:,:)  :: tracers, centres
  real*8, dimension(:), allocatable :: DD, delta
  real*8, dimension(:), allocatable :: rbin, rbin_edges

  character(20), external :: str
  character(len=500) :: data_filename, data_filename_2, output_filename
  character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char, ngrid_char, box_char

  logical :: debug = .true.
  
  if (debug) then
    if (iargc() .lt. 8) then
        write(*,*) 'Some arguments are missing.'
        write(*,*) '1) data_filename'
        write(*,*) '2) data_filename_2'
        write(*,*) '3) output_filename'
        write(*,*) '4) boxsize'
        write(*,*) '5) dim1_min'
        write(*,*) '6) dim1_max'
        write(*,*) '7) dim1_nbin'
        write(*,*) '8) ngrid'
        write(*,*) ''
        stop
    end if
  end if

  call system_clock(beginning, rate)
    
  call get_command_argument(number=1, value=data_filename)
  call get_command_argument(number=2, value=data_filename_2)
  call get_command_argument(number=3, value=output_filename)
  call get_command_argument(number=4, value=box_char)
  call get_command_argument(number=5, value=dim1_min_char)
  call get_command_argument(number=6, value=dim1_max_char)
  call get_command_argument(number=7, value=dim1_nbin_char)
  call get_command_argument(number=8, value=ngrid_char)
  
  read(box_char, *) boxsize
  read(dim1_min_char, *) dim1_min
  read(dim1_max_char, *) dim1_max
  read(dim1_nbin_char, *) dim1_nbin
  read(ngrid_char, *) ngrid

  if (debug) then
    write(*,*) '-----------------------'
    write(*,*) 'Running tpcf_2d.exe'
    write(*,*) 'input parameters:'
    write(*,*) ''
    write(*, *) 'data_filename: ', trim(data_filename)
    write(*, *) 'data_filename_2: ', trim(data_filename_2)
    write(*, *) 'boxsize: ', trim(box_char)
    write(*, *) 'output_filename: ', trim(output_filename)
    write(*, *) 'dim1_min: ', trim(dim1_min_char), ' Mpc'
    write(*, *) 'dim1_max: ', trim(dim1_max_char), ' Mpc'
    write(*, *) 'dim1_nbin: ', trim(dim1_nbin_char)
    write(*, *) 'ngrid: ', trim(ngrid_char)
    write(*,*) ''
  end if

  open(10, file=data_filename, status='old', form='unformatted')
  read(10) nrows
  read(10) ncols
  allocate(tracers(ncols, nrows))
  read(10) tracers
  close(10)
  ntracers = nrows
  if (debug) then
    write(*,*) 'ntracers dim: ', size(tracers, dim=1), size(tracers, dim=2)
    write(*,*) 'pos(min), pos(max) = ', minval(tracers(1,:)), maxval(tracers(1,:))
  end if

  open(11, file=data_filename_2, status='old', form='unformatted')
  read(11) nrows
  read(11) ncols
  allocate(centres(ncols, nrows))
  read(11) centres
  close(11)
  ncentres = nrows
  if (debug) then
    write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)
  end if

  allocate(rbin(dim1_nbin))
  allocate(rbin_edges(dim1_nbin + 1))
  allocate(DD(dim1_nbin))
  allocate(delta(dim1_nbin))
  
  rwidth = (dim1_max - dim1_min) / dim1_nbin
  do i = 1, dim1_nbin + 1
    rbin_edges(i) = dim1_min+(i-1)*rwidth
  end do
  do i = 1, dim1_nbin
    rbin(i) = rbin_edges(i+1)-rwidth/2.
  end do
  
  ! Mean density inside the box
  rhomean = ntracers / (boxsize ** 2)
  
  ! Construct linked list for tracers
  allocate(lirst(ngrid, ngrid))
  allocate(nlirst(ngrid, ngrid))
  allocate(ll(ntracers))
  rgrid = (boxsize) / real(ngrid)
  
  lirst = 0
  ll = 0
  
  do i = 1, ntracers
    indx = int((tracers(1, i)) / rgrid + 1.)
    indy = int((tracers(2, i)) / rgrid + 1.)
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid) lirst(indx,indy)=i
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid) nlirst(indx,indy) = &
    nlirst(indx, indy) + 1
  end do
  
  do i = 1, ntracers
    indx = int((tracers(1, i))/ rgrid + 1.)
    indy = int((tracers(2, i))/ rgrid + 1.)
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid) then
      ll(lirst(indx,indy)) = i
      lirst(indx,indy) = i
    endif
  end do
  
  if (debug) then
    write(*,*) 'Linked list successfully constructed'
    write(*,*) ''
    write(*,*) 'Starting loop over tracers...'
  end if
  
  DD = 0
  delta = 0
  dim1_min2 = dim1_min ** 2
  dim1_max2 = dim1_max ** 2
  ndif = int(dim1_max / rgrid + 1.)
  
  do i = 1, ncentres
    ipx = int(centres(1, i) / rgrid + 1.)
    ipy = int(centres(2, i) / rgrid + 1.)

    do ix = ipx - ndif, ipx + ndif
      do iy = ipy - ndif, ipy + ndif
        if ((ix - ipx)**2 + (iy - ipy)**2 .gt. (ndif + 1)**2) cycle

        ix2 = ix
        iy2 = iy

        if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
        if (ix2 .lt. 1) ix2 = ix2 + ngrid
        if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
        if (iy2 .lt. 1) iy2 = iy2 + ngrid

        ii = lirst(ix2,iy2)
        if(ii.ne.0) then
          do
            ii = ll(ii)
            disx = tracers(1, ii) - centres(1, i)
            disy = tracers(2, ii) - centres(2, i)

            if (disx .lt. -boxsize/2) disx = disx + boxsize
            if (disx .gt. boxsize/2) disx = disx - boxsize
            if (disy .lt. -boxsize/2) disy = disy + boxsize
            if (disy .gt. boxsize/2) disy = disy - boxsize

            dis2 = disx ** 2 + disy ** 2

            if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
              dis = sqrt(dis2)
              rind = int((dis - dim1_min) / rwidth + 1)
              DD(rind) = DD(rind) + 1
            end if

            if(ii.eq.lirst(ix2,iy2)) exit

          end do
        end if
      end do
    end do
  end do

  do i = 1, dim1_nbin
    area = pi * (rbin_edges(i + 1) ** 2 - rbin_edges(i) ** 2)
    delta(i) = DD(i) / (area * rhomean * ncentres) - 1
  end do
  
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
    write(*,*) 'Elapsed time: ', real(end - beginning) / real(rate)
  end if

end program tpcf
