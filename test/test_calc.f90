program test_calc
  use, intrinsic :: iso_fortran_env
  use ising2d_m
  implicit none
  integer(int64), parameter :: nx = 11, ny = 10, nall = nx * ny, mcs = 100
  real(real64), parameter :: kbt = 2.0_real64
  type(ising2d) :: system
  integer(int32) :: i
  call system%init(nx, ny, kbt)
  if (system%calc_energy_summ() /= -2 * nall) then
     error stop 1
  end if
  if (system%calc_magne_summ() /= nall) then
     error stop 2
  end if
  do i = 1, mcs
     call system%update()
  end do
  write(error_unit, *) system%calc_energy_summ()
  write(error_unit, *) system%calc_magne_summ()
end program test_calc
