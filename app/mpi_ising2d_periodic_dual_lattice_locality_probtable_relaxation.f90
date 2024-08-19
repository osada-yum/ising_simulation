program mpi_ising2d_periodic_dual_lattice_locality_probtable_relaxation
  use, intrinsic :: iso_fortran_env
  use mpi
  use gf2xe
  use msmt19937
  use ising2d_periodic_dual_lattice_locality_probtable_m
  use variance_covariance_kahan_m
  use mpi_variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: outs(*) = [output_unit, error_unit]

  character(len=*), parameter :: version = "dual_lattice_locality_probtable"

  integer(int32), parameter :: mcs = 1000
  integer(int32), parameter :: nsample = 25

  integer(int64), parameter :: iseed = 42_int64
  integer(int32), parameter :: n_skip = 0_int32
  integer(int32) :: expo

  type(variance_covariance_kahan), allocatable :: order_parameter(:)

  integer(int32) :: i, j

  integer(int32) :: myrank, num_proc, ierr
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierr)

  call init_genrand(iseed)
  !> Skip random numbers. (num_proc * n_skip + myrank) * 2^e
  !> 2^e must be larger than (nx * ny * mcs * nsample).
  expo = ceiling(log(real(2 * nx * ny * (mcs + 1) * nsample, real64)) / log(2.0d0)) + 1
  if (num_proc * n_skip + myrank /= 0) &
       & call mt_jumpahead(num_proc * n_skip + myrank, expo)

  call init_ising2d()
  if (myrank == 0) then
     call print_version()
     do i = 1, size(outs)
        write(outs(i), '(a)') "#"//version
        write(outs(i), '(a,i0)'    ) "# Nsize: ", nall
        write(outs(i), '(2(a, i0))') "# nx: ", nx, " ny: ", ny
        write(outs(i), '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
        write(outs(i), '(a, g0)' ) "# 温度: ", kbt
        write(outs(i), '(a)' ) "# method: Metropolis"
        write(outs(i), '(a, i0)' ) "# the number of processors: ", num_proc
     end do
  end if

  allocate(order_parameter(mcs), source = variance_covariance_kahan())
  do j = 1, nsample
     if (myrank == 0) &
          write(error_unit, '(a, i0)') "sample: ", j
     call init_ising2d_order()
     do i = 1, mcs
        call update_metropolis()
        associate(m => calc_magne(), &
             & e => calc_energy())
          call order_parameter(i)%add_data(m, e)
        end associate
     end do
  end do
  block
    type(variance_covariance_kahan), allocatable :: all_order_params(:)
    allocate(all_order_params(mcs), source = variance_covariance_kahan())
    call vck_mpi_multi_gather(mcs, order_parameter, all_order_params, 0, myrank, num_proc, ierr)
    if (myrank == 0) then
       write(output_unit, '(a)') "# Nsize, Nsample, mcs, <m>, <e>, <me>, <m^2>, <e^2>, χ, C, m'"
       do i = 1, mcs
          write(output_unit, '(*(g0, 1x))') nall, all_order_params(i)%num_sample(), i, &
               & all_order_params(i)%mean1(), all_order_params(i)%mean2(), all_order_params(i)%mean_v1v2(), &
               & all_order_params(i)%square_mean1(), all_order_params(i)%square_mean2(), &
               & nall * all_order_params(i)%var1(), nall * all_order_params(i)%var2(), &
               & nall * all_order_params(i)%cov()
       end do
    end if
  end block
  call MPI_Finalize(ierr)
end program mpi_ising2d_periodic_dual_lattice_locality_probtable_relaxation
