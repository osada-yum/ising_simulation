program mpi_ising2d_simulation
  use, intrinsic :: iso_fortran_env
  use mpi_f08
  use ising2d_m
  use variance_covariance_kahan_m
  use mpi_variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 1000, nsample = 25
  integer(int64), parameter :: nx = 1001, ny = 1000
  real(real64), parameter :: kbt = 2.26918531421302_real64, n_inv_r64 = 1 / real(nx * ny, real64)
  type(ising2d) :: system
  type(variance_covariance_kahan) :: order_parameter(mcs)
  integer(int32) :: i, j
  integer(int32) :: myrank, num_proc, ierr
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierr)
  call system%init(nx, ny, kbt)
  if (myrank == 0) then
     write(output_unit, '(a,i0)'    ) "# Nsize: ", system%nall()
     write(output_unit, '(2(a, i0))') "# nx: ", system%nx(), " ny: ", system%ny()
     write(output_unit, '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
     write(output_unit, '(a, g0)' ) "# 温度: ", system%kbt()
     write(output_unit, '(a)' ) "# method: Metropolis"
     write(output_unit, '(a, i0)' ) "# the number of processors: ", num_proc
  end if
  do j = 1, nsample
     if (myrank == 0) &
          write(error_unit, '(a, i0)') "sample: ", j
     call system%set_ising_allup()
     do i = 1, mcs
        call system%update()
        associate(m => system%calc_magne_summ(), &
             & e => system%calc_energy_summ())
          call order_parameter(i)%add_data(m * n_inv_r64, e * n_inv_r64)
        end associate
     end do
  end do
  block
    type(variance_covariance_kahan) :: all_order_params(mcs)
    do i = 1, mcs
       call vck_mpi_gather(order_parameter(i), all_order_params(i), 0, myrank, num_proc, ierr)
    end do
    if (myrank == 0) then
       write(output_unit, '(a)') "# Nsize, Nsample, mcs, <m>, <e>, <m^2>, <e^2>, χ, C, m'"
       do i = 1, mcs
          write(output_unit, '(*(g0, 1x))') system%nall(), all_order_params(i)%num_sample(), i, &
               & all_order_params(i)%mean1(), all_order_params(i)%mean2(), &
               & all_order_params(i)%square_mean1(), all_order_params(i)%square_mean2(), &
               & all_order_params(i)%var1(), all_order_params(i)%var2(), &
               & all_order_params(i)%cov()
       end do
    end if
  end block
  call MPI_Finalize(ierr)
end program mpi_ising2d_simulation
