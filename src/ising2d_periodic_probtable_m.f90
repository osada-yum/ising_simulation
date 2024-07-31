module ising2d_periodic_probtable_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  character(len=*), parameter :: version = "probtable"
  public :: print_version

  integer(int64), parameter :: nx = 1000_int64, ny = nx, nall = nx * ny
  real(real64), parameter :: nall_inv = 1d0 / nall
  real(real64), parameter :: kbt = 2.269185314213022d0, beta = 1 / kbt

  integer(int32), allocatable :: ising2d(:, :)

  real(real64), allocatable :: rnds(:, :)

  integer(int64) :: i_c
  real(real64), parameter :: probability_table(1:8) = [(exp(- beta * i_c), i_c = 1, 8)]

  integer(int64), parameter :: nd = 4
  integer(int64), parameter :: dy(nd) = [0_int64, 1_int64, 0_int64, -1_int64]
  integer(int64), parameter :: dx(nd) = [1_int64, 0_int64, -1_int64, 0_int64]

  public :: nx, ny, nall, kbt, beta
  public :: init_ising2d, init_ising2d_order, update_metropolis, calc_energy, calc_magne

contains
  impure subroutine print_version()
    write(output_unit, '(a)') "#"//version
    write(error_unit, '(a)') "#"//version
  end subroutine print_version
  !> init_ising2d: Initialize 6-state clock model.
  impure subroutine init_ising2d()
    allocate(ising2d(nx, ny))
    allocate(rnds(nx, ny))
  end subroutine init_ising2d
  !> init_ising2d_order: Set spins with the all-alinged state.
  impure subroutine init_ising2d_order()
    ising2d(:, :) = 1_int32
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
    do y = 1, ny
       do x = 1 + iand(y - 1, b'1'), nx, 2
          call local_flip(x, y)
       end do
    end do
    !> odd.
    do y = 1, ny
       do x = 1 + iand(y, b'1'), nx, 2
          call local_flip(x, y)
       end do
    end do
  end subroutine update_metropolis
  !> local_flip: Flip a spin of (x, y) position of lattice.
  !> @param x A x-position of the spin.
  !> @param y A y-position of the spin.
  impure subroutine local_flip(x, y)
    integer(int64), intent(in) :: x, y
    integer(int32) :: nearest_states(nd)
    integer(int64) :: near_x, near_y
    integer(int32) :: delta_e
    real(real64) :: prob
    integer(int32) :: d
    do d = 1, nd
       near_x = x + dx(d)
       if (near_x > nx) then
          near_x = 1_int64
       else if (near_x < 1_int64) then
          near_x = nx
       end if
       near_y = y + dy(d)
       if (near_y > ny) then
          near_y = 1_int64
       else if (near_y < 1_int64) then
          near_y = ny
       end if
       nearest_states(d) = ising2d(near_x, near_y)
    end do
    delta_e = 2 * ising2d(x, y) * sum(nearest_states(1:4))
    if (delta_e > 0) then
       prob = probability_table(delta_e)
       if (.not. (rnds(x, y) < prob)) &
            & return
    end if
    !> delta_e <= 0 .or. rnds(x, y) < exp(- beta * delta_e)
    ising2d(x, y) = - ising2d(x, y)
  end subroutine local_flip
  !> calc_energy: Calculate the energy density.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: summ
    integer(int64) :: x, y
    integer(int64) :: rx, uy
    summ = 0_int64
    do y = 1, ny
       uy = y + 1
       if (uy > ny) uy = 1_int64
       do x = 1, nx
          rx = x + 1
          if (rx > nx) rx = 1_int64
          summ = summ - ising2d(x, y) * (ising2d(x, uy) + ising2d(rx, y))
       end do
    end do
    res = summ * nall_inv
  end function calc_energy
  !> calc_magne: Calculate the magnetism density.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: summ
    integer(int64) :: x, y
    summ = 0_int64
    do y = 1, ny
       do x = 1, nx
          summ = summ + ising2d(x, y)
       end do
    end do
    res = summ * nall_inv
  end function calc_magne
end module ising2d_periodic_probtable_m
