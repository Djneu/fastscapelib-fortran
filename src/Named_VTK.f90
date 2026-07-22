subroutine Fastscape_Named_VTK (f, vex, istep, foldername, k, model_height, model_dim, adjustment)

  ! Writes an XML StructuredGrid (.vts) file per timestep into
  ! <foldername>/, preserving the regular nx*ny lattice (connectivity
  ! is implicit, so no Cells block is needed).
  !
  ! Companion routine Fastscape_PVD_Collection maintains a .pvd that
  ! indexes each .vts with its physical time, giving ParaView a proper
  ! time series that plots directly alongside ASPECT's solution.pvd.
  !
  ! This version writes ASCII (inline) data: simplest to get correct and
  ! trivially inspectable. Move to appended-binary later if you want to
  ! match the file size / speed of the original binary .vtk output.

  use FastScapeContext
  implicit none

  integer, intent(in) :: k, istep, model_dim
  double precision, intent(in) :: vex, model_height, adjustment
  double precision, intent(in), dimension(*) :: f
  character(len=k), intent(in) :: foldername
  character cstep*7

  integer :: i, j
  double precision :: dx, dy
  character(len=1024) :: fname
  character(len=64)   :: extent

  dx = xl/(nx - 1)
  dy = yl/(ny - 1)

  ! zero-padded step string (matches your original convention)
  write (cstep,'(i7)') istep
  if (istep.lt.10)      cstep(1:6)='000000'
  if (istep.lt.100)     cstep(1:5)='00000'
  if (istep.lt.1000)    cstep(1:4)='0000'
  if (istep.lt.10000)   cstep(1:3)='000'
  if (istep.lt.100000)  cstep(1:2)='00'
  if (istep.lt.1000000) cstep(1:1)='0'

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
         ! 2D ASPECT: elevation in Y, Z = 0, adjust for ghost nodes and extent in 2D.
         write(77,'(3(1x,ES14.6))') &
              sngl(dx*(i-1)-adjustment), &
              sngl((h(i+(j-1)*nx)+model_height)*abs(vex)), &
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
  call write_scalar(77, 'HHHHH',         f,      nn)
  call write_scalar(77, 'basement',      b,      nn)
  call write_scalar(77, 'erosion_rate',  erate,  nn)
  call write_scalar(77, 'total_erosion', etot,   nn)
  call write_scalar(77, 'drainage_area', a,      nn)
  call write_scalar(77, 'catchment',     catch,  nn)
  call write_scalar(77, 'precipitation',     precip,  nn)

  ! Should these be renamed based on orientation, should z always read uplift?
  call write_scalar(77, 'velocity_z',  u  ,  nn)
  call write_scalar(77, 'velocity_x',  vx  ,  nn)
  call write_scalar(77, 'velocity_y',  vy  ,  nn)
  call write_scalar(77, 'diffusivity',  kd  ,  nn)
  call write_scalar(77, 'river_incision_rate',  kf  ,  nn)
  call write_scalar(77, 'velocity_y',  vy  ,  nn)
  write(77,'(A)') '      </PointData>'

  write(77,'(A)') '    </Piece>'
  write(77,'(A)') '  </StructuredGrid>'
  write(77,'(A)') '</VTKFile>'
  close(77)

  ! ---- Optional basement / sea-level files, mirroring the original ----
  if (vex.lt.0.d0) then

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
     call write_scalar(77, 'B',     b, nn)
     call write_scalar(77, 'HHHHH', f, nn)
     write(77,'(A)') '      </PointData>'
     write(77,'(A)') '    </Piece>'
     write(77,'(A)') '  </StructuredGrid>'
     write(77,'(A)') '</VTKFile>'
     close(77)

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

end subroutine Fastscape_Named_VTK


subroutine Fastscape_PVD_Collection (istep, time, foldername, k)

  ! Maintains <foldername>/Topography.pvd indexing each per-step .vts
  ! with its physical time. Rewrites the whole (small) file each call so
  ! the collection is always valid even if the run stops early.

  implicit none
  integer, intent(in) :: k, istep
  double precision, intent(in) :: time
  character(len=k), intent(in) :: foldername

  integer, parameter :: MAXSTEPS = 100000
  integer, save :: nrec = 0
  integer, save :: steps(MAXSTEPS)
  double precision, save :: times(MAXSTEPS)
  character cstep*7
  integer :: r
  character(len=1024) :: fname

  ! record this step
  nrec = nrec + 1
  if (nrec > MAXSTEPS) then
     write(*,*) 'Fastscape_PVD_Collection: MAXSTEPS exceeded, increase it.'
     nrec = MAXSTEPS
     return
  end if
  steps(nrec) = istep
  times(nrec) = time

  fname = trim(foldername)//'/Topography.pvd'
  open(unit=78, file=trim(fname), status='unknown', form='formatted')
  write(78,'(A)') '<?xml version="1.0"?>'
  write(78,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
  write(78,'(A)') '  <Collection>'
  do r = 1, nrec
     write (cstep,'(i7)') steps(r)
     if (steps(r).lt.10)      cstep(1:6)='000000'
     if (steps(r).lt.100)     cstep(1:5)='00000'
     if (steps(r).lt.1000)    cstep(1:4)='0000'
     if (steps(r).lt.10000)   cstep(1:3)='000'
     if (steps(r).lt.100000)  cstep(1:2)='00'
     if (steps(r).lt.1000000) cstep(1:1)='0'
     write(78,'(A,ES16.8,A)') '    <DataSet timestep="', times(r), &
          '" group="" part="0" file="Topography'//cstep//'.vts"/>'
  end do
  write(78,'(A)') '  </Collection>'
  write(78,'(A)') '</VTKFile>'
  close(78)

  return
end subroutine Fastscape_PVD_Collection
