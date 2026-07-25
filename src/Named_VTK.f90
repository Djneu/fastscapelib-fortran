subroutine Fastscape_Named_VTK (vex, istep, foldername, k, model_height, model_dim, adjustment, &
                                time, output_basement, output_sealevel)

  ! Writes an XML StructuredGrid (.vts) file per timestep into
  ! the fastscape folder of the ASPECT output folder. This will include
  ! a .pvd file to write the times and plot them along-side ASPECT.
  !
  ! We only get one step per call and don't remember the others, so a
  ! hidden .pvd file lists them all. We save each step's time to this file
  ! and rebuild the .pvd from it. This allows the .pvd to work through restarts.
  !
  ! Adjustment and model_height are small factors used to correctly plot
  ! FastScape surface in relation to the ASPECT surface, with adjustment doing the shift
  ! to account for the ghost nodes if they are turned on.
  !
  ! Model dim is used to switch the surface plot. In 2D ASPECT simulations,
  ! depth is called Y in ASPECT and FastScape dimensions are called X and Z,
  ! in 3D ASPECT Z is depth in ASPECT and FastScape is X-Y.

  use FastScapeContext
  use, intrinsic :: iso_c_binding
  implicit none

  integer, intent(in) :: k, istep, model_dim
  double precision, intent(in) :: vex, model_height, adjustment, time
  character(len=k), intent(in) :: foldername
  logical(c_bool), intent(in) :: output_basement, output_sealevel

  character(len=7) :: cstep
  integer :: i, j
  double precision :: dx, dy
  character(len=1024) :: fname
  character(len=64)   :: extent

  ! Points and per-node scratch buffer (single precision, VTK ordering)
  real(c_float), allocatable :: pts(:), buf(:)

  dx = xl/(nx - 1)
  dy = yl/(ny - 1)

  ! zero-padded step string, built once and reused everywhere
  write (cstep,'(i7.7)') istep

  ! StructuredGrid extent: 0-based index range in i, j, k
  write(extent,'(I0,1x,I0,1x,I0,1x,I0,1x,I0,1x,I0)') 0, nx-1, 0, ny-1, 0, 0

  allocate(pts(3*nn), buf(nn))

  ! =========================== Topography ===========================
  ! Build the point coordinates (i fastest, then j) into pts(:).
  call build_points(h, pts)

  fname = trim(foldername)//'/Topography'//cstep//'.vts'
  call open_vts(77, fname, extent)

  ! Offsets are assigned in write order; open_piece_points writes the
  ! Points DataArray header with offset 0, then each scalar advances the
  ! running offset by 4 (length header) + 4*nn (data).
  call write_header_points(77, extent)

  ! Here are the different fields we will write into the 
  ! Topography.vts files.
  call put_line(77, '      <PointData Scalars="topography">')
  call decl_scalar(77, 'topography',          off_after(0))
  call decl_scalar(77, 'basement',            off_after(1))
  call decl_scalar(77, 'erosion_rate',        off_after(2))
  call decl_scalar(77, 'total_erosion',       off_after(3))
  call decl_scalar(77, 'drainage_area',       off_after(4))
  call decl_scalar(77, 'catchment',           off_after(5))
  call decl_scalar(77, 'precipitation',       off_after(6))
  call decl_scalar(77, 'uplift',              off_after(7))
  call decl_scalar(77, 'velocity_0',          off_after(8))
  call decl_scalar(77, 'velocity_1',          off_after(9))
  call decl_scalar(77, 'diffusivity',         off_after(10))
  call decl_scalar(77, 'river_incision_rate', off_after(11))
  call put_line(77, '      </PointData>')
  call close_piece(77)

  ! ---- Appended raw binary block ----
  call open_appended(77)
  call append_real(77, pts)          ! points (3*nn)
  call append_real_from(77, h,      buf)
  call append_real_from(77, b,      buf)
  call append_real_from(77, erate,  buf)
  call append_real_from(77, etot,   buf)
  call append_real_from(77, a,      buf)
  call append_real_from(77, catch,  buf)
  call append_real_from(77, precip, buf)
  call append_real_from(77, u,      buf)
  call append_real_from(77, vx,     buf)
  call append_real_from(77, vy,     buf)
  call append_real_from(77, kd,     buf)
  call append_real_from(77, kf,     buf)
  call close_appended(77)
  close(77)

  ! ============================ Basement ============================
  if (output_basement) then
     call build_points(b, pts)

     fname = trim(foldername)//'/Basement'//cstep//'.vts'
     call open_vts(77, fname, extent)
     call write_header_points(77, extent)
     call put_line(77, '      <PointData Scalars="basement">')
     call decl_scalar(77, 'basement',            off_after(0))
     call decl_scalar(77, 'basement_difference', off_after(1))
     call put_line(77, '      </PointData>')
     call close_piece(77)

     call open_appended(77)
     call append_real(77, pts)               ! points
     call append_real_from(77, b, buf)       ! basement (buf as scratch)
     do i = 1, nn
        buf(i) = real(h(i) - b(i), c_float)   ! basement_difference into buf
     end do
     call append_real(77, buf)               ! basement_difference
     call close_appended(77)
     close(77)
  end if

  ! ============================ SeaLevel ============================
  if (output_sealevel) then
     ! Flat plane: elevation is the constant sealevel at every node.
     call build_points_const(real(sealevel*abs(vex)+model_height, c_float), pts)
     do i = 1, nn
        buf(i) = real(sealevel, c_float)
     end do

     fname = trim(foldername)//'/SeaLevel'//cstep//'.vts'
     call open_vts(77, fname, extent)
     call write_header_points(77, extent)
     call put_line(77, '      <PointData Scalars="SL">')
     call decl_scalar(77, 'SL', off_after(0))
     call put_line(77, '      </PointData>')
     call close_piece(77)

     call open_appended(77)
     call append_real(77, pts)
     call append_real(77, buf)
     call close_appended(77)
     close(77)
  end if

  ! This will call the function to write the PVD files.
  call Fastscape_PVD_Collection (cstep, time, foldername, k, output_basement, output_sealevel)

  deallocate(pts, buf)
  return

contains

  ! ----- byte offset of the array that comes AFTER `m` scalars have
  !       already been declared (points array is index -1 at offset 0) -----
  integer function off_after(m) result(o)
    integer, intent(in) :: m
    ! points block: 8 (len) + 4*3*nn (data); then m scalar blocks each
    ! 8 + 4*nn.  (8-byte UInt64 length prefix.)
    o = (8 + 4*3*nn) + m*(8 + 4*nn)
  end function off_after

  subroutine build_points(elev, p)
    double precision, intent(in), dimension(*) :: elev
    real(c_float), intent(out) :: p(:)
    integer :: ii, jj, idx
    idx = 0
    do jj = 1, ny
       do ii = 1, nx
          if (model_dim == 2) then
             p(idx+1) = real(dx*(ii-1)-adjustment, c_float)
             p(idx+2) = real(elev(ii+(jj-1)*nx)*abs(vex)+model_height, c_float)
             p(idx+3) = real(dy*(jj-1)-yl-adjustment, c_float)
          else
             p(idx+1) = real(dx*(ii-1)-adjustment, c_float)
             p(idx+2) = real(dy*(jj-1)-adjustment, c_float)
             p(idx+3) = real(elev(ii+(jj-1)*nx)*abs(vex)+model_height, c_float)
          end if
          idx = idx + 3
       end do
    end do
  end subroutine build_points

  subroutine build_points_const(zconst, p)
    real(c_float), intent(in) :: zconst
    real(c_float), intent(out) :: p(:)
    integer :: ii, jj, idx
    idx = 0
    do jj = 1, ny
       do ii = 1, nx
          if (model_dim == 2) then
             p(idx+1) = real(dx*(ii-1)-adjustment, c_float)
             p(idx+2) = zconst
             p(idx+3) = real(dy*(jj-1)-yl-adjustment, c_float)
          else
             p(idx+1) = real(dx*(ii-1)-adjustment, c_float)
             p(idx+2) = real(dy*(jj-1)-adjustment, c_float)
             p(idx+3) = zconst
          end if
          idx = idx + 3
       end do
    end do
  end subroutine build_points_const

  ! ---- XML header helpers (all written as formatted text) ----

  subroutine open_vts(unit, path, ext)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: path, ext
    ! stream access so we can mix the text header with the raw binary
    ! appended block in one file.
    open(unit=unit, file=trim(path), status='unknown', form='unformatted', &
         access='stream', convert='big_endian')
    call put_line(unit, '<?xml version="1.0"?>')
    call put_line(unit, '<VTKFile type="StructuredGrid" version="0.1" '// &
                        'byte_order="BigEndian" header_type="UInt64">')
    call put_line(unit, '  <StructuredGrid WholeExtent="'//trim(ext)//'">')
    call put_line(unit, '    <Piece Extent="'//trim(ext)//'">')
  end subroutine open_vts

  subroutine write_header_points(unit, ext)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: ext
    call put_line(unit, '      <Points>')
    call put_line(unit, '        <DataArray type="Float32" '// &
                        'NumberOfComponents="3" format="appended" offset="0"/>')
    call put_line(unit, '      </Points>')
  end subroutine write_header_points

  subroutine decl_scalar(unit, name, offset)
    integer, intent(in) :: unit, offset
    character(len=*), intent(in) :: name
    character(len=32) :: offc
    write(offc,'(I0)') offset
    call put_line(unit, '        <DataArray type="Float32" Name="'//trim(name)// &
                        '" format="appended" offset="'//trim(adjustl(offc))//'"/>')
  end subroutine decl_scalar

  subroutine close_piece(unit)
    integer, intent(in) :: unit
    call put_line(unit, '    </Piece>')
    call put_line(unit, '  </StructuredGrid>')
  end subroutine close_piece

  subroutine open_appended(unit)
    integer, intent(in) :: unit
    ! The underscore marks the start of the raw data; VTK offsets are
    ! measured from the byte immediately after it.
    call put_line(unit, '  <AppendedData encoding="raw">')
    call put_str (unit, '_')
  end subroutine open_appended

  subroutine close_appended(unit)
    integer, intent(in) :: unit
    call put_line(unit, '  </AppendedData>')
    call put_line(unit, '</VTKFile>')
  end subroutine close_appended

  ! ---- raw byte writers (stream unformatted) ----

  subroutine put_line(unit, s)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: s
    write(unit) s//char(10)
  end subroutine put_line

  subroutine put_str(unit, s)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: s
    write(unit) s
  end subroutine put_str

  ! append a real(c_float) array: 4-byte length header then the raw data
  subroutine append_real(unit, arr)
    integer, intent(in) :: unit
    real(c_float), intent(in) :: arr(:)
    integer(c_int64_t) :: nbytes
    nbytes = int(size(arr)*4, c_int64_t)
    write(unit) nbytes
    write(unit) arr
  end subroutine append_real

  ! append a double-precision node field, converting to single first
  subroutine append_real_from(unit, arr, scratch)
    integer, intent(in) :: unit
    double precision, intent(in), dimension(*) :: arr
    real(c_float), intent(inout) :: scratch(:)
    integer :: ii
    integer(c_int64_t) :: nbytes
    do ii = 1, nn
       scratch(ii) = real(arr(ii), c_float)
    end do
    nbytes = int(nn*4, c_int64_t)
    write(unit) nbytes
    write(unit) scratch(1:nn)
  end subroutine append_real_from

!--------------------------------------------------------------------

  subroutine Fastscape_PVD_Collection (cstep, time, foldername, k, output_basement, output_sealevel)

    ! Rewrites the hidden file with this step appended, dropping any
    ! existing entries at or after the current time so that a replayed
    ! step (e.g. after a checkpoint restart) replaces its old record
    ! instead of duplicating it. Then recreates each .pvd from the
    ! index.

    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in) :: k
    character(len=7), intent(in) :: cstep
    double precision, intent(in) :: time
    character(len=k), intent(in) :: foldername
    logical(c_bool), intent(in) :: output_basement, output_sealevel

    character(len=7), allocatable :: cs(:)
    double precision, allocatable :: t(:)
    character(len=7) :: cs1
    double precision :: t1
    integer :: n, ios, r
    logical :: exists

    ! read existing index, keeping only entries strictly before this time
    allocate(cs(100000), t(100000))
    n = 0
    inquire(file=trim(foldername)//'/.pvd_index', exist=exists)
    if (exists) then
       open(unit=79, file=trim(foldername)//'/.pvd_index', status='old', &
            form='formatted', action='read')
       do
          read(79,*,iostat=ios) cs1, t1
          if (ios /= 0) exit
          if (t1 .lt. time - 1.d-10*max(1.d0,abs(time))) then
             n = n + 1
             cs(n) = cs1
             t(n)  = t1
          end if
       end do
       close(79)
    end if

    ! append this step and rewrite the index
    n = n + 1
    cs(n) = cstep
    t(n)  = time

    open(unit=79, file=trim(foldername)//'/.pvd_index', status='replace', &
         form='formatted')
    do r = 1, n
       ! We write at G0.15 to match ASPECT significant figures in solution.pvd
       write(79,'(A7,1x,G0.15)') cs(r), t(r)
    end do
    close(79)

    deallocate(cs, t)

    ! Always write the topography pvd.
    call write_pvd (foldername, k, 'Topography')

    ! Only write sea level and basement pvd if they are requested.
    if (output_sealevel) call write_pvd (foldername, k, 'SeaLevel')
    if (output_basement) call write_pvd (foldername, k, 'Basement')

    return
  end subroutine Fastscape_PVD_Collection

!--------------------------------------------------------------------

  subroutine write_pvd (foldername, k, base)

    ! Rewrites <base>.pvd in full from the hidden file index, so the
    ! collection is always valid even if the run stops early.

    implicit none
    integer, intent(in) :: k
    character(len=k), intent(in) :: foldername
    character(len=*), intent(in) :: base

    character(len=1024) :: fname
    character(len=7) :: cs
    double precision :: t
    integer :: ios

    open(unit=79, file=trim(foldername)//'/.pvd_index', status='old', &
         form='formatted', action='read')

    fname = trim(foldername)//'/'//trim(base)//'.pvd'
    open(unit=78, file=trim(fname), status='unknown', form='formatted')
    write(78,'(A)') '<?xml version="1.0"?>'
    write(78,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="BigEndian">'
    write(78,'(A)') '  <Collection>'
    do
       read(79,*,iostat=ios) cs, t
       if (ios /= 0) exit
       write(78,'(A,G0.12,A)') '    <DataSet timestep="', t, &
            '" group="" part="0" file="'//trim(base)//cs//'.vts"/>'
    end do
    write(78,'(A)') '  </Collection>'
    write(78,'(A)') '</VTKFile>'
    close(78)
    close(79)

    return
  end subroutine write_pvd

end subroutine Fastscape_Named_VTK