subroutine Fastscape_Named_VTK (vex, istep, foldername, k, model_height, model_dim, adjustment, &
                                time, output_basement, output_sealevel)

  ! Writes an XML StructuredGrid (.vts) file per timestep into
  ! <foldername>/.
  !
  ! Companion routine Fastscape_PVD_Collection maintains a .pvd that
  ! indexes each .vts with its physical time, giving ParaView a proper
  ! time series that plots directly alongside ASPECT's solution.pvd.
  !
  ! Adjustment and model_height are small factors used to correctly plot
  ! FastScape surface in relation to the ASPECT surface, with adjustment doing the shift
  ! to account for the ghost nodes if they are turned on.
  !
  ! Model dim is used to switch the surface plot. In 2D ASPECT Y is depth in ASPECT and FastScape is X-Z,
  ! in 3D ASPECT Z is depth in ASPET and FastScape is X-Y.

  use FastScapeContext
  use, intrinsic :: iso_c_binding
  implicit none

  integer, intent(in) :: k, istep, model_dim
  double precision, intent(in) :: vex, model_height, adjustment, time
  character(len=k), intent(in) :: foldername
  character(len=7) :: cstep
  integer :: i, j
  double precision :: dx, dy 
  character(len=1024) :: fname
  character(len=64)   :: extent
  logical(c_bool), intent(in) :: output_basement, output_sealevel

  dx = xl/(nx - 1)
  dy = yl/(ny - 1)

  ! zero-padded step string, built once and reused everywhere
  write (cstep,'(i7.7)') istep

  ! StructuredGrid extent: 0-based index range in i, j, k
  write(extent,'(I0,1x,I0,1x,I0,1x,I0,1x,I0,1x,I0)') 0, nx-1, 0, ny-1, 0, 0

  fname = trim(foldername)//'/Topography'//cstep//'.vts'
  open(unit=77, file=trim(fname), status='unknown', form='formatted')

  write(77,'(A)') '<?xml version="1.0"?>'
  write(77,'(A)') '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(77,'(A)') '  <StructuredGrid WholeExtent="'//trim(extent)//'">'
  write(77,'(A)') '    <Piece Extent="'//trim(extent)//'">'

  ! ---- Points (x, y, z = elevation * vertical exaggeration) ----
  ! Order must be i fastest, then j, then k (VTK structured ordering).
  write(77,'(A)') '      <Points>'
  write(77,'(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
  do j = 1, ny
     do i = 1, nx
      if (model_dim == 2) then
         ! 2D ASPECT: elevation in Y, adjust for ghost nodes and extent in 2D.
         write(77,'(3(1x,ES14.6))') &
              sngl(dx*(i-1)-adjustment), &
              sngl((h(i+(j-1)*nx))*abs(vex)+model_height), &
              sngl(dy*(j-1)-yl-adjustment)
      else
         ! 3D ASPECT: X, Y as-is, elevation in Z
         write(77,'(3(1x,ES14.6))') &
              sngl(dx*(i-1)-adjustment), &
              sngl(dy*(j-1)-adjustment), &
              sngl((h(i+(j-1)*nx))*abs(vex)+model_height)
      end if
     end do
  end do
  write(77,'(A)') '        </DataArray>'
  write(77,'(A)') '      </Points>'

  ! ---- Point data (all your original fields, same node ordering) ----
  write(77,'(A)') '      <PointData Scalars="topography">'
  call write_scalar(77, 'topography',    h,      nn)
  call write_scalar(77, 'basement',      b,      nn)
  call write_scalar(77, 'erosion_rate',  erate,  nn)
  call write_scalar(77, 'total_erosion', etot,   nn)
  call write_scalar(77, 'drainage_area', a,      nn)
  call write_scalar(77, 'catchment',     catch,  nn)
  call write_scalar(77, 'precipitation', precip, nn)

  ! Should these be renamed based on orientation, should z always read uplift?
  call write_scalar(77, 'uplift',          u,   nn)
  call write_scalar(77, 'velocity_0',          vx,  nn)
  call write_scalar(77, 'velocity_1',          vy,  nn)
  call write_scalar(77, 'diffusivity',         kd,  nn)
  call write_scalar(77, 'river_incision_rate', kf,  nn)
  write(77,'(A)') '      </PointData>'

  write(77,'(A)') '    </Piece>'
  write(77,'(A)') '  </StructuredGrid>'
  write(77,'(A)') '</VTKFile>'
  close(77)

  ! ---- Optionally output the basement ----
  if (output_basement) then

     fname = trim(foldername)//'/Basement'//cstep//'.vts'
     open(unit=77, file=trim(fname), status='unknown', form='formatted')
     write(77,'(A)') '<?xml version="1.0"?>'
     write(77,'(A)') '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
     write(77,'(A)') '  <StructuredGrid WholeExtent="'//trim(extent)//'">'
     write(77,'(A)') '    <Piece Extent="'//trim(extent)//'">'
     write(77,'(A)') '      <Points>'
     write(77,'(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
     do j = 1, ny
        do i = 1, nx
          if (model_dim == 2) then
            write(77,'(3(1x,ES14.6))') sngl(dx*(i-1)-adjustment), &
                  sngl(b(i+(j-1)*nx)*abs(vex)+model_height), &
                  sngl(dy*(j-1)-yl-adjustment)
          else
            write(77,'(3(1x,ES14.6))') sngl(dx*(i-1)-adjustment), sngl(dy*(j-1)-adjustment), &
                  sngl(b(i+(j-1)*nx)*abs(vex)+model_height)
          endif
        end do
     end do
     write(77,'(A)') '        </DataArray>'
     write(77,'(A)') '      </Points>'
     write(77,'(A)') '      <PointData Scalars="B">'
     call write_scalar(77, 'basement',     b, nn)
     call write_scalar(77, 'basement_difference', h-b, nn)
     write(77,'(A)') '      </PointData>'
     write(77,'(A)') '    </Piece>'
     write(77,'(A)') '  </StructuredGrid>'
     write(77,'(A)') '</VTKFile>'
     close(77)

   end if 

   ! ---- Optionally output the sealevel ----
   if (output_sealevel) then

     fname = trim(foldername)//'/SeaLevel'//cstep//'.vts'
     open(unit=77, file=trim(fname), status='unknown', form='formatted')
     write(77,'(A)') '<?xml version="1.0"?>'
     write(77,'(A)') '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
     write(77,'(A)') '  <StructuredGrid WholeExtent="'//trim(extent)//'">'
     write(77,'(A)') '    <Piece Extent="'//trim(extent)//'">'
     write(77,'(A)') '      <Points>'
     write(77,'(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
     do j = 1, ny
        do i = 1, nx
          if (model_dim == 2) then
            write(77,'(3(1x,ES14.6))') sngl(dx*(i-1)-adjustment), sngl(sealevel*abs(vex)+model_height), &
                  sngl(dy*(j-1)-yl-adjustment)
          else
            write(77,'(3(1x,ES14.6))') sngl(dx*(i-1)-adjustment), sngl(dy*(j-1)-adjustment), &
                  sngl(sealevel*abs(vex)+model_height)
          endif
        end do
     end do
     write(77,'(A)') '        </DataArray>'
     write(77,'(A)') '      </Points>'
     write(77,'(A)') '      <PointData Scalars="SL">'
     write(77,'(A)') '        <DataArray type="Float32" Name="SL" format="ascii">'
     do i = 1, nn
        write(77,'(1x,ES14.6)') sngl(sealevel)
     end do
     write(77,'(A)') '        </DataArray>'
     write(77,'(A)') '      </PointData>'
     write(77,'(A)') '    </Piece>'
     write(77,'(A)') '  </StructuredGrid>'
     write(77,'(A)') '</VTKFile>'
     close(77)

  end if

  call Fastscape_PVD_Collection (cstep, time, foldername, k, vex)

  return

contains

  subroutine write_scalar(unit, name, arr, n)
    integer, intent(in) :: unit, n
    character(len=*), intent(in) :: name
    double precision, intent(in), dimension(*) :: arr
    integer :: ii
    write(unit,'(A)') '        <DataArray type="Float32" Name="'//trim(name)// &
         '" format="ascii">'
    do ii = 1, n
       write(unit,'(1x,ES14.6)') sngl(arr(ii))
    end do
    write(unit,'(A)') '        </DataArray>'
  end subroutine write_scalar

!--------------------------------------------------------------------

  subroutine Fastscape_PVD_Collection (cstep, time, foldername, k, vex)

    ! Rewrites the sidecar index with this step appended, dropping any
    ! existing entries at or after the current time so that a replayed
    ! step (e.g. after a checkpoint restart) replaces its old record
    ! instead of duplicating it. Then regenerates each .pvd from the
    ! index. Holds no state between calls: the index file on disk is
    ! the only record, so restarts pick up where they left off.

    implicit none
    integer, intent(in) :: k
    character(len=7), intent(in) :: cstep
    double precision, intent(in) :: time, vex
    character(len=k), intent(in) :: foldername

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
          read(79,'(A7,1x,ES16.8)',iostat=ios) cs1, t1
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
       write(79,'(A7,1x,ES16.8)') cs(r), t(r)
    end do
    close(79)

    deallocate(cs, t)

    call write_pvd (foldername, k, 'Topography')

    if (vex.lt.0.d0) then
       call write_pvd (foldername, k, 'SeaLevel')
       call write_pvd (foldername, k, 'Basement')
    end if

    return
  end subroutine Fastscape_PVD_Collection

!--------------------------------------------------------------------

  subroutine write_pvd (foldername, k, base)

    ! Rewrites <base>.pvd in full from the sidecar index, so the
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
    write(78,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
    write(78,'(A)') '  <Collection>'
    do
       read(79,'(A7,1x,ES16.8)',iostat=ios) cs, t
       if (ios /= 0) exit
       write(78,'(A,ES16.8,A)') '    <DataSet timestep="', t, &
            '" group="" part="0" file="'//trim(base)//cs//'.vts"/>'
    end do
    write(78,'(A)') '  </Collection>'
    write(78,'(A)') '</VTKFile>'
    close(78)
    close(79)

    return
  end subroutine write_pvd

end subroutine Fastscape_Named_VTK