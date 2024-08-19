module ising2d_periodic_dual_lattice_probtable_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  character(len=*), parameter :: version = "dual_lattice_probtable"
  public :: print_version

  integer(int64), parameter :: nx = 1000_int64, ny = nx, nall = nx * ny
  real(real64), parameter :: nall_inv = 1d0 / nall
  real(real64), parameter :: kbt = 2.269185314213022d0, beta = 1 / kbt

  integer(int32), allocatable :: ising2d_even(:, :)
  integer(int32), allocatable :: ising2d_odd(:, :)

  real(real64), allocatable :: rnds(:, :)

  integer(int64) :: i_c
  real(real64), parameter :: probability_table(1:8) = [(exp(- beta * i_c), i_c = 1, 8)]

  integer(int64), parameter :: nd = 4

  public :: nx, ny, nall, kbt, beta
  public :: init_ising2d, init_ising2d_order, update_metropolis, calc_energy, calc_magne

contains
  impure subroutine print_version()
    write(output_unit, '(a)') "#"//version
    write(error_unit, '(a)') "#"//version
  end subroutine print_version
  !> init_ising2d: Initialize 6-state clock model.
  impure subroutine init_ising2d()
    allocate(ising2d_even(nx, ny / 2))
    allocate(ising2d_odd(nx, ny / 2))
    allocate(rnds(nx, ny))
  end subroutine init_ising2d
  !> init_ising2d_order: Set spins with the all-alinged state.
  impure subroutine init_ising2d_order()
    ising2d_even(:, :) = 1_int32
    ising2d_odd(:, :) = 1_int32
  end subroutine init_ising2d_order
  !> update_metropolis: Update the lattice with 1MCS.
  impure subroutine update_metropolis()
    integer(int64) :: x, y
    do y = 1, ny
       do x = 1, nx
          rnds(x, y) = grnd()
          rnds(x, y) = grnd()
       end do
    end do
    !> even.
    do y = 1, ny / 2
       do x = 1, nx
          call local_flip_even(ising2d_even, ising2d_odd, x, y)
       end do
    end do
    !> odd.
    do y = 1, ny / 2
       do x = 1, nx
          call local_flip_odd(ising2d_odd, ising2d_even, x, y)
       end do
    end do
  end subroutine update_metropolis
  !> local_flip_even: Flip a spin of (x, y) position of lattice in even lattice.
  !> @param ising2d_even A array of spins in even lattice.
  !> @param ising2d_odd A array of spins in odd lattice.
  !> @param x A x-position of the spin.
  !> @param y A y-position of the spin.
  pure subroutine local_flip_even(ising2d_even, ising2d_odd, x, y)
    integer(int32), intent(inout) :: ising2d_even(nx, ny / 2)
    integer(int32), intent(in) :: ising2d_odd(nx, ny / 2)
    integer(int64), intent(in) :: x, y
    integer(int32) :: nearest_states(nd)
    integer(int64) :: rx, lx, uy, dy
    integer(int32) :: delta_e
    real(real64) :: prob
    !> left
    lx = x - 1
    if (lx < 1) lx = nx
    !> right
    rx = x + 1
    if (rx > nx) rx = 1_int64
    !> down
    dy = y - iand(x, b'1')
    if (dy < 1) dy = ny / 2
    !> up
    uy = y + iand(x - 1, b'1')
    if (uy > ny / 2) uy = 1_int64
    nearest_states(1:4) = [ising2d_odd(lx, y), ising2d_odd(rx, y), ising2d_odd(x, dy), ising2d_odd(x, uy)]

    delta_e = 2 * ising2d_even(x, y) * sum(nearest_states(1:4))
    if (delta_e > 0) then
       prob = probability_table(delta_e)
       if (.not. (rnds(x, 2 * y - iand(x, b'1')) < prob)) &
            & return
    end if
    !> delta_e <= 0 .or. rnds(x, y) < exp(- beta * delta_e)
    ising2d_even(x, y) = - ising2d_even(x, y)
  end subroutine local_flip_even
  !> local_flip_odd: Flip a spin of (x, y) position of lattice in odd lattice.
  !> @param ising2d_odd A array of spins in odd lattice.
  !> @param ising2d_even A array of spins in even lattice.
  !> @param x A x-position of the spin.
  !> @param y A y-position of the spin.
  pure subroutine local_flip_odd(ising2d_odd, ising2d_even, x, y)
    integer(int32), intent(inout) :: ising2d_odd(nx, ny / 2)
    integer(int32), intent(in) :: ising2d_even(nx, ny / 2)
    integer(int64), intent(in) :: x, y
    integer(int32) :: nearest_states(nd)
    integer(int64) :: rx, lx, uy, dy
    integer(int32) :: delta_e
    real(real64) :: prob
    !> left
    lx = x - 1
    if (lx < 1) lx = nx
    !> right
    rx = x + 1
    if (rx > nx) rx = 1_int64
    !> down
    dy = y - iand(x - 1, b'1')
    if (dy < 1) dy = ny / 2
    !> up
    uy = y + iand(x, b'1')
    if (uy > ny / 2) uy = 1_int64
    nearest_states(1:4) = [ising2d_even(lx, y), ising2d_even(rx, y), ising2d_even(x, dy), ising2d_even(x, uy)]

    delta_e = 2 * ising2d_odd(x, y) * sum(nearest_states(1:4))
    if (delta_e > 0) then
       prob = probability_table(delta_e)
       if (.not. (rnds(x, 2 * y - iand(x - 1, b'1')) < prob)) &
            & return
    end if
    !> delta_e <= 0 .or. rnds(x, y) < exp(- beta * delta_e)
    ising2d_odd(x, y) = - ising2d_odd(x, y)
  end subroutine local_flip_odd

  !> calc_energy: Calculate the energy density.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: summ
    integer(int64) :: x, y
    integer(int64) :: rx, uy
    summ = 0_int64
    do y = 1, ny / 2
       do x = 1, nx
          rx = x + 1
          if (rx > nx) rx = 1_int64

          uy = y + iand(x - 1, b'1')
          if (uy > ny / 2) uy = 1_int64
          summ = summ - ising2d_even(x, y) * (ising2d_odd(x, uy) + ising2d_odd(rx, y))

          uy = y + iand(x, b'1')
          if (uy > ny / 2) uy = 1_int64
          summ = summ - ising2d_odd(x, y) * (ising2d_even(x, uy) + ising2d_even(rx, y))
       end do
    end do
    res = summ * nall_inv
  end function calc_energy
  !> calc_magne: Calculate the magnetism density.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: summ
    integer(int64) :: x, y
    summ = 0_int64
    do y = 1, ny / 2
       do x = 1, nx
          summ = summ + ising2d_even(x, y) + ising2d_odd(x, y)
       end do
    end do
    res = summ * nall_inv
  end function calc_magne
end module ising2d_periodic_dual_lattice_probtable_m
