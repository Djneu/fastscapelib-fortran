subroutine Fastscape_Write_Visualization (vex, istep, foldername, k, model_height, model_dim, adjustment, &
                                time, output_basement, output_sealevel)

  ! Writes FastScape output for one timestep as VTS files that ParaView
  ! can plot alongside ASPECT's solution.pvd:
  !   Topography<step_str>.vts  - the surface with all FastScape field_data
  !   Basement<step_str>.vts    - the bedrock surface        (if requested)
  !   SeaLevel<step_str>.vts    - a flat plane at sea level  (if requested)
  ! plus a matching .pvd per surface listing every step_str's physical time.
  !
  ! Each .vts stores its numbers as raw binary appended after the XML
  ! header (big-endian Float32). write_vts_file() handles that format;
  ! this routine just builds the point coordinates and field arrays and
  ! hands them over. update_pvd() maintains the .pvd time series.
  !
  ! model_dim picks the vertical axis: in 2D ASPECT height runs in the Y-direction (FastScape
  ! spans X-Z), and in 3D in the Z-direction (FastScape spans X-Y). vex exaggerates
  ! elevation (only |vex| is used); model_height and adjustment shift the
  ! surface to sit correctly relative to ASPECT's mesh.

  use FastScapeContext
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer,          intent(in) :: k, istep, model_dim
  double precision, intent(in) :: vex, model_height, adjustment, time
  character(len=k), intent(in) :: foldername
  character(len=20), allocatable :: field_names(:)
  logical(c_bool),  intent(in) :: output_basement, output_sealevel
  double precision, target, dimension(:), allocatable :: kfint, kdint

  character(len=7)  :: step_str
  double precision  :: dx, dy
  integer           :: node, n_fields

  ! Point coordinates (3 per node) and field values, single precision.
  real(c_float), allocatable :: xyz(:,:)          ! (3, nn)
  real(c_float), allocatable :: field_data(:,:)       ! (nn, n_fields)
  allocate(kfint(nn))
  allocate(kdint(nn))

  dx = xl/(nx - 1)
  dy = yl/(ny - 1)
  write (step_str,'(i7.7)') istep

  allocate(xyz(3, nn))

  n_fields = 12
  if (runMarine) n_fields = 13

  ! ----- Topography: surface at h, carrying every FastScape field -----
  call build_points(h, xyz)
  allocate(field_data(nn, n_fields))
  field_data(:, 1)  = real(h(1:nn),      c_float) ! topography
  field_data(:, 2)  = real(b(1:nn),      c_float) ! basement
  field_data(:, 3)  = real(erate(1:nn),  c_float) ! erosion_rate
  field_data(:, 4)  = real(etot(1:nn),   c_float) ! total erosion
  field_data(:, 5)  = real(a(1:nn),      c_float) ! drainage area
  field_data(:, 6)  = real(catch(1:nn),  c_float) ! catchment
  field_data(:, 7)  = real(precip(1:nn), c_float) ! precipitation
  field_data(:, 8)  = real(u(1:nn),      c_float) ! uplift
  field_data(:, 9)  = real(vx(1:nn),     c_float) ! velocity_0
  field_data(:, 10) = real(vy(1:nn),     c_float) ! velocity_1

  ! Find kf and kd including sediment values. 
  kfint(1:nn) = kf(1:nn)
  where (kfsed > 0.0d0 .and. (h(1:nn) - b(1:nn)) > 1.0d0)    kfint(1:nn) = kfsed

  kdint(1:nn) = kd(1:nn)
  where (kdsed > 0.0d0 .and. (h(1:nn) - b(1:nn)) > 1.0d-6) kdint(1:nn) = kdsed

  field_data(:, 11) = real(kdint(1:nn), c_float)   ! diffusivity
  field_data(:, 12) = real(kfint(1:nn), c_float)   ! river incision

  ! Silt fraction: only exists as a field when marine is active.
  if (runMarine) then
    where (h(1:nn) < sealevel)
      ! Marine diffusivity, matching how sea diffusivity is calculated in ASPECT.
      ! where Fmix is the silt_fraction, kdsea1 is the silt diffusivity, and ksea2 is the sand diffusivity.
      field_data(:, 11) = real(kdsea1(1:nn)*Fmix(1:nn) + (1.0_c_float - Fmix(1:nn))*kdsea2(1:nn), c_float)  ! diffusivity
      field_data(:, 12) = ieee_value(0.0_c_float, ieee_quiet_nan)
      field_data(:, 13) = real(Fmix(1:nn), c_float)              ! silt_fraction
    elsewhere
      field_data(:, 13) = ieee_value(0.0_c_float, ieee_quiet_nan)
    end where
  end if

  allocate(field_names(n_fields))
  field_names(1:12) = [character(len=20) :: 'topography','basement','erosion_rate','total_erosion', &
                     'drainage_area','catchment','precipitation','uplift','velocity_0','velocity_1', &
                     'diffusivity','river_incision_rate']
  if (runMarine) field_names(13) = 'silt_fraction'

  call write_vts_file(trim(foldername)//'/Topography'//step_str//'.vts', xyz, field_data, field_names)

  deallocate(field_names)
  deallocate(field_data)

  ! ----- Basement: surface at b, with its depth below topography -----
  if (output_basement) then
     call build_points(b, xyz)
     allocate(field_data(nn, 2))
     field_data(:, 1) = real(b(1:nn),        c_float)  ! basement
     field_data(:, 2) = real(h(1:nn)-b(1:nn), c_float) ! topography - basement
     call write_vts_file(trim(foldername)//'/Basement'//step_str//'.vts', xyz, field_data, &
          [character(len=20) :: 'basement','basement_depth'])
     deallocate(field_data)
  end if

  ! ----- SeaLevel: flat plane at the current sea level -----
  if (output_sealevel) then
     call build_flat_points(real(sealevel*abs(vex)+model_height, c_float), xyz)
     allocate(field_data(nn, 1))
     field_data(:, 1) = real(sealevel, c_float) ! sea level
     call write_vts_file(trim(foldername)//'/SeaLevel'//step_str//'.vts', xyz, field_data, &
          [character(len=20) :: 'sea_level'])
     deallocate(field_data)
  end if

  ! Function to update the pvd files.
  call update_pvd(step_str, time, foldername, output_basement, output_sealevel)

  deallocate(xyz)
  return

contains

  ! Fill xyz(1:3, node) for every grid node in VTS order (i fastest,
  ! then j). The elevation array sets the height on whichever axis
  ! ASPECT treats as vertical; the other two axes are the FastScape grid.
  subroutine build_points(elevation, xyz)
    double precision, intent(in), dimension(*) :: elevation
    real(c_float),    intent(out) :: xyz(:,:)
    integer :: ii, jj, node
    node = 0
    do jj = 1, ny
       do ii = 1, nx
          node = node + 1
          if (model_dim == 2) then
             xyz(1, node) = real(dx*(ii-1)-adjustment, c_float)
             xyz(2, node) = real(elevation(ii+(jj-1)*nx)*abs(vex)+model_height, c_float)
             xyz(3, node) = real(dy*(jj-1)-yl-adjustment, c_float)
          else
             xyz(1, node) = real(dx*(ii-1)-adjustment, c_float)
             xyz(2, node) = real(dy*(jj-1)-adjustment, c_float)
             xyz(3, node) = real(elevation(ii+(jj-1)*nx)*abs(vex)+model_height, c_float)
          end if
       end do
    end do
  end subroutine build_points

  ! Same layout as build_points but every node sits at one height.
  subroutine build_flat_points(height, xyz)
    real(c_float), intent(in)  :: height
    real(c_float), intent(out) :: xyz(:,:)
    integer :: ii, jj, node
    node = 0
    do jj = 1, ny
       do ii = 1, nx
          node = node + 1
          if (model_dim == 2) then
             xyz(1, node) = real(dx*(ii-1)-adjustment, c_float)
             xyz(2, node) = height
             xyz(3, node) = real(dy*(jj-1)-yl-adjustment, c_float)
          else
             xyz(1, node) = real(dx*(ii-1)-adjustment, c_float)
             xyz(2, node) = real(dy*(jj-1)-adjustment, c_float)
             xyz(3, node) = height
          end if
       end do
    end do
  end subroutine build_flat_points

  ! Write one complete .vts file: XML header naming each field, then a
  ! raw binary block holding the points and every field. Each array in
  ! the block is prefixed by its length as an 8-byte integer, and the
  ! header gives each array's byte offset into the block. Data is
  ! big-endian Float32; convert='big_endian' forces that regardless of
  ! the GFORTRAN_CONVERT_UNIT environment variable.
  subroutine write_vts_file(path, xyz, field_data, names)
    character(len=*), intent(in) :: path
    real(c_float),    intent(in) :: xyz(:,:)          ! (3, nn)
    real(c_float),    intent(in) :: field_data(:,:)       ! (nn, n_fields)
    character(len=*), intent(in) :: names(:)

    integer, parameter :: u = 77
    integer :: f, n_fields, offset
    character(len=64)  :: extent
    character(len=32)  :: offc
    integer(c_int64_t) :: n_bytes

    n_fields = size(names)
    write(extent,'(I0,1x,I0,1x,I0,1x,I0,1x,I0,1x,I0)') 0, nx-1, 0, ny-1, 0, 0

    open(unit=u, file=path, status='unknown', form='unformatted', &
         access='stream', convert='big_endian')

    ! --- XML header ---
    call line(u, '<?xml version="1.0"?>')
    call line(u, '<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian" header_type="UInt64">')
    call line(u, '  <StructuredGrid WholeExtent="'//trim(extent)//'">')
    call line(u, '    <Piece Extent="'//trim(extent)//'">')
    call line(u, '      <Points>')
    call line(u, '        <DataArray type="Float32" NumberOfComponents="3" format="appended" offset="0"/>')
    call line(u, '      </Points>')
    call line(u, '      <PointData Scalars="'//trim(names(1))//'">')
    ! points block occupies 8 + 4*3*nn bytes; each field 8 + 4*nn.
    offset = 8 + 4*3*nn
    do f = 1, n_fields
       write(offc,'(I0)') offset
       call line(u, '        <DataArray type="Float32" Name="'//trim(names(f))// &
                 '" format="appended" offset="'//trim(adjustl(offc))//'"/>')
       offset = offset + 8 + 4*nn
    end do
    call line(u, '      </PointData>')
    call line(u, '    </Piece>')
    call line(u, '  </StructuredGrid>')

    ! --- raw binary block: points, then each field ---
    call line(u, '  <AppendedData encoding="raw">')
    write(u) '_'                              ! marks offset origin; no newline
    n_bytes = int(3*nn*4, c_int64_t)
    write(u) n_bytes
    write(u) xyz
    do f = 1, n_fields
       n_bytes = int(nn*4, c_int64_t)
       write(u) n_bytes
       write(u) field_data(:, f)
    end do
    write(u) char(10)
    call line(u, '  </AppendedData>')
    call line(u, '</VTKFile>')
    close(u)
  end subroutine write_vts_file

  ! Write one text line plus a newline byte to a stream unit.
  subroutine line(u, text)
    integer,          intent(in) :: u
    character(len=*), intent(in) :: text
    write(u) text//char(10)
  end subroutine line

  ! Append this step_str to the hidden .pvd_index file (dropping any entry at
  ! or after the current time, so a restart replaces rather than
  ! duplicates), then rebuild each .pvd from the full index. The index
  ! lives on disk so the .pvd survives checkpoint restarts. Times are
  ! written to match ASPECT's solution.pvd precision so ParaView lines
  ! the two series up on the same steps.
  subroutine update_pvd(step_str, time, foldername, output_basement, output_sealevel)
    character(len=7), intent(in) :: step_str
    double precision, intent(in) :: time
    character(len=*), intent(in) :: foldername
    logical(c_bool),  intent(in) :: output_basement, output_sealevel

    integer, parameter :: max_records = 100000
    character(len=7)  :: steps(max_records), s
    double precision  :: times(max_records), t
    integer           :: count, ios, r
    logical           :: exists
    character(len=1024) :: index_path

    index_path = trim(foldername)//'/.pvd_index'

    count = 0
    inquire(file=trim(index_path), exist=exists)
    if (exists) then
       open(unit=79, file=trim(index_path), status='old', form='formatted', action='read')
       do
          ! * is used to here to match the G0.15 format we use to write to match ASPECT.
          read(79,*,iostat=ios) s, t
          if (ios /= 0) exit
          if (t .lt. time - 1.d-10*max(1.d0,abs(time))) then
             count = count + 1
             steps(count) = s
             times(count) = t
          end if
       end do
       close(79)
    end if
    count = count + 1
    steps(count) = step_str
    times(count) = time

    open(unit=79, file=trim(index_path), status='replace', form='formatted')
    do r = 1, count
       write(79,'(A7,1x,G0.15)') steps(r), times(r)   ! G0.15 matches ASPECT
    end do
    close(79)

    call write_one_pvd(foldername, index_path, 'Topography')
    if (output_sealevel) call write_one_pvd(foldername, index_path, 'SeaLevel')
    if (output_basement) call write_one_pvd(foldername, index_path, 'Basement')
  end subroutine update_pvd

  ! Rebuild one <base>.pvd from every entry in the index file.
  subroutine write_one_pvd(foldername, index_path, base)
    character(len=*), intent(in) :: foldername, index_path, base
    integer :: ios2
    character(len=7) :: s2
    double precision :: t2
    open(unit=79, file=trim(index_path), status='old', form='formatted', action='read')
    open(unit=78, file=trim(foldername)//'/'//trim(base)//'.pvd', status='unknown', form='formatted')
    write(78,'(A)') '<?xml version="1.0"?>'
    write(78,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="BigEndian">'
    write(78,'(A)') '  <Collection>'
    do
       read(79,*,iostat=ios2) s2, t2
       if (ios2 /= 0) exit
       write(78,'(A,G0.12,A)') '    <DataSet timestep="', t2, &   ! G0.12 matches ASPECT
            '" group="" part="0" file="'//trim(base)//s2//'.vts"/>'
    end do
    write(78,'(A)') '  </Collection>'
    write(78,'(A)') '</VTKFile>'
    close(78)
    close(79)
  end subroutine write_one_pvd

end subroutine Fastscape_Write_Visualization
