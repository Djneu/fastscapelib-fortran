program MarineDiffusionRotation

  implicit none

  integer :: nx, ny, i, j, ij
  integer :: istep, nstep

  logical, parameter :: rotate = .true.

  double precision :: xl, yl, dt
  double precision :: sealevel
  double precision :: poro, zporo, ratio, L
  double precision :: kfsed, kdsed, g, m, n

  double precision, allocatable :: x(:),y(:)
  double precision, allocatable :: h(:),u(:)
  double precision, allocatable :: kf(:),kd(:)
  double precision, allocatable :: kdsea1(:),kdsea2(:)

  call FastScape_Init()

!----------------------------------------------------------
! Grid
!----------------------------------------------------------

  if (.not.rotate) then
     nx=201
     ny=101
     xl=400.d3
     yl=200.d3
  else
     nx=101
     ny=201
     xl=200.d3
     yl=400.d3
  endif

  call FastScape_Set_NX_NY(nx,ny)
  call FastScape_Setup()
  call FastScape_Set_XL_YL(xl,yl)

  allocate(x(nx*ny),y(nx*ny))
  allocate(h(nx*ny))
  allocate(u(nx*ny))
  allocate(kf(nx*ny),kd(nx*ny))
  allocate(kdsea1(nx*ny),kdsea2(nx*ny))

!----------------------------------------------------------
! Coordinates
!----------------------------------------------------------

  x=(/((xl*dble(i-1)/dble(nx-1),i=1,nx),j=1,ny)/)
  y=(/((yl*dble(j-1)/dble(ny-1),i=1,nx),j=1,ny)/)

!----------------------------------------------------------
! Initial marine domain
!----------------------------------------------------------

  h=-100.d0

! Gaussian bump

  if (.not.rotate) then

     h = h + 50.d0*exp( &
          -((x-xl/2.d0)/(0.15d0*xl))**2 &
          -((y-yl/2.d0)/(0.15d0*yl))**2 )

  else

     h = h + 50.d0*exp( &
          -((x-xl/2.d0)/(0.15d0*xl))**2 &
          -((y-yl/2.d0)/(0.15d0*yl))**2 )

  endif

  call FastScape_Init_H(h)

!----------------------------------------------------------
! No uplift
!----------------------------------------------------------

  u=0.d0
  call FastScape_Set_U(u)

!----------------------------------------------------------
! Turn hillslope diffusion off
!----------------------------------------------------------

  kf=0.d0
  kd=0.d0

  kfsed=0.d0
  kdsed=0.d0

  m=0.5d0
  n=1.d0
  g=1.d0

  call FastScape_Set_Erosional_Parameters( &
       kf,kfsed,m,n,kd,kdsed,g,g,1.d0)

!----------------------------------------------------------
! Variable marine diffusivities
!----------------------------------------------------------

  if (.not.rotate) then

     kdsea1 = 200.d0 + 300.d0*x/xl
     kdsea2 = 100.d0 + 150.d0*x/xl

  else

     kdsea1 = 200.d0 + 300.d0*y/yl
     kdsea2 = 100.d0 + 150.d0*y/yl

  endif

  sealevel=0.d0
  poro=0.d0
  zporo=1000.d0
  ratio=0.5d0
  L=100.d0

  call FastScape_Set_Marine_Parameters( &
       sealevel,poro,poro,zporo,zporo, &
       ratio,L,kdsea1,kdsea2)

!----------------------------------------------------------
! Time step
!----------------------------------------------------------

  dt=100.d0
  call FastScape_Set_DT(dt)

  call FastScape_Set_BC(1111)

!----------------------------------------------------------
! Run
!----------------------------------------------------------

  nstep=200

  do

     call FastScape_Get_Step(istep)
     if (istep>=nstep) exit

     call FastScape_Execute_Step()

  enddo

!----------------------------------------------------------
! Output
!----------------------------------------------------------

  call FastScape_Copy_H(h)

  if (.not.rotate) then
     open(10,file='topography_A.dat')
  else
     open(10,file='topography_B.dat')
  endif

  do ij=1,nx*ny
     write(10,'(3ES20.10)') x(ij),y(ij),h(ij)
  enddo

  close(10)

  call FastScape_Destroy()

end program MarineDiffusionRotation