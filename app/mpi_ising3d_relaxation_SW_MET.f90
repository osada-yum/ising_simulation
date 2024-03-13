program mpi_ising3d_simulation
  use, intrinsic :: iso_fortran_env
  use mpi
  use ising3d_m
  use variance_covariance_kahan_m
  use mpi_variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 1000
  integer(int32), parameter :: mcs_sw = 1
  integer(int32), parameter :: nsample = 25
  integer(int64), parameter :: nx = 31
  integer(int64), parameter :: ny = 31
  integer(int64), parameter :: nz = 30
  real(real64), parameter :: kbt = 4.511_real64
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny * nz, real64)
  type(ising3d) :: system
  type(variance_covariance_kahan), allocatable :: order_parameter(:)
  integer(int32) :: i, j
  integer(int32) :: myrank, num_proc, ierr
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierr)
  call system%init(nx, ny, nz, kbt)
  if (myrank == 0) then
     write(output_unit, '(a,i0)'    ) "# Nsize: ", system%nall()
     write(output_unit, '(3(a, i0))') "# nx: ", system%nx(), " ny: ", system%ny(), " nz: ", system%nz()
     write(output_unit, '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
     write(output_unit, '(a, g0)' ) "# 温度: ", system%kbt()
     write(output_unit, '(a, i0)' ) "# method: SW: ", mcs_sw
     write(output_unit, '(a, i0)' ) "# method: Metropolis: ", mcs - mcs_sw
     write(output_unit, '(a, i0)' ) "# the number of processors: ", num_proc

     write(error_unit, '(a,i0)'    ) "# Nsize: ", system%nall()
     write(error_unit, '(3(a, i0))') "# nx: ", system%nx(), " ny: ", system%ny(), " nz: ", system%nz()
     write(error_unit, '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
     write(error_unit, '(a, g0)' ) "# 温度: ", system%kbt()
     write(error_unit, '(a, i0)' ) "# method: SW: ", mcs_sw
     write(error_unit, '(a, i0)' ) "# method: Metropolis: ", mcs - mcs_sw
     write(error_unit, '(a, i0)' ) "# the number of processors: ", num_proc
  end if
  allocate(order_parameter(mcs), source = variance_covariance_kahan())
  do j = 1, nsample
     if (myrank == 0) &
          write(error_unit, '(2(a, i0))') "sample: ", j, " / ", nsample
     call system%set_ising_allup()
     do i = 1, mcs_sw
        call system%update_swendsen_wang()
        associate(m => system%calc_magne_summ(), &
             & e => system%calc_energy_summ())
          if (m < 0) &
               & call system%flip_all()
          call order_parameter(i)%add_data(abs(m) * n_inv_r64, e * n_inv_r64)
        end associate
     end do
     do i = mcs_sw + 1, mcs
        call system%update()
        associate(m => system%calc_magne_summ(), &
             & e => system%calc_energy_summ())
          call order_parameter(i)%add_data(m * n_inv_r64, e * n_inv_r64)
        end associate
     end do
  end do
  block
    type(variance_covariance_kahan), allocatable :: all_order_params(:)
    allocate(all_order_params(mcs), source = variance_covariance_kahan())
    call vck_mpi_multi_gather(mcs, order_parameter, all_order_params, 0, myrank, num_proc, ierr)
    if (myrank == 0) then
       write(output_unit, '(a)') "# Nsize, Nsample, mcs, <m>, <e>, <m^2>, <e^2>, χ, C, m'"
       do i = 1, mcs
          write(output_unit, '(*(g0, 1x))') system%nall(), all_order_params(i)%num_sample(), i, &
               & all_order_params(i)%mean1(), all_order_params(i)%mean2(), &
               & all_order_params(i)%square_mean1(), all_order_params(i)%square_mean2(), &
               & system%nall() * all_order_params(i)%var1(), system%nall() * all_order_params(i)%var2(), &
               & system%nall() * all_order_params(i)%cov()
       end do
    end if
  end block
  call MPI_Finalize(ierr)
end program mpi_ising3d_simulation
