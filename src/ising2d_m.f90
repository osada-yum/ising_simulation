module ising2d_m
  use, intrinsic :: iso_fortran_env
  implicit none
  private
  integer(int32), parameter :: lb_exparr = -8, ub_exparr = 8
  public :: ising2d
  type :: ising2d
     private
     integer(int64) :: nx_, ny_, nall_
     real(real64) :: beta_
     integer(int32), allocatable :: spins_(:)
     real(real64), allocatable :: exparr_(:)
   contains
     !> initializer.
     procedure, pass :: init => init_ising2d
     !> setter.
     procedure, pass :: set_ising_allup => set_ising_allup_ising2d
     procedure, pass :: set_ising_random => set_ising_random_ising2d
     procedure, pass :: set_kbt => set_kbt_ising2d
     procedure, pass :: set_beta => set_beta_ising2d
     !> updater.
     procedure, pass :: update => update_ising2d
     procedure, pass, private :: update_onesite => update_onesite_ising2d
     procedure, pass, private :: update_norishiro => update_norishiro_ising2d
     !> calculator.
     procedure, pass :: calc_energy_summ => calc_total_energy_ising2d
     procedure, pass :: calc_magne_summ => calc_total_magne_ising2d
     !> getter.
     procedure, pass :: nx => nx_ising2d
     procedure, pass :: ny => ny_ising2d
     procedure, pass :: nall => nall_ising2d
     procedure, pass :: kbt => kbt_ising2d
     procedure, pass :: beta => beta_ising2d
     procedure, pass, private :: norishiro_begin => norishiro_begin_ising2d
     procedure, pass, private :: norishiro_end => norishiro_end_ising2d
  end type ising2d
contains
  !> init_ising2d: Initialize ising2d once.
  pure subroutine init_ising2d(this, nx, ny, kbt)
    class(ising2d), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny
    real(real64), intent(in) :: kbt
    this%nx_ = nx
    this%ny_ = ny
    this%nall_ = nx * ny
    allocate(this%spins_(this%norishiro_begin() : this%norishiro_end()))
    call this%set_ising_allup()
    allocate(this%exparr_(lb_exparr:ub_exparr))
    call this%set_kbt(kbt)
  end subroutine init_ising2d
  !> set_ising_allup_ising2d: Set spins `1`.
  pure subroutine set_ising_allup_ising2d(this)
    class(ising2d), intent(inout) :: this
    this%spins_(:) = 1_int32
  end subroutine set_ising_allup_ising2d
  !> set_ising_random_ising2d: Set spins `1` or `-1` randomly.
  impure subroutine set_ising_random_ising2d(this)
    class(ising2d), intent(inout) :: this
    real(real64), allocatable :: r(:)
    allocate(r(1:this%nall_))
    call random_number(r)
    this%spins_(1:this%nall_) = merge(1_int32, -1_int32, r < 0.5_real64)
    call this%update_norishiro()
  end subroutine set_ising_random_ising2d
  !> set_kbt_ising2d: Set parameter `beta` as `1 / kbt`.
  pure subroutine set_kbt_ising2d(this, kbt)
    class(ising2d), intent(inout) :: this
    real(real64), intent(in) :: kbt
    call this%set_beta(1 / kbt)
  end subroutine set_kbt_ising2d
  !> set_beta_ising2d: Set parameter `beta` and update `this%exparr_`.
  pure subroutine set_beta_ising2d(this, beta)
    class(ising2d), intent(inout) :: this
    real(real64), intent(in) :: beta
    integer(int32) :: delta_e
    this%beta_ = beta
    this%exparr_(:) = 1.0_real64
    do delta_e = 1, ub_exparr
       this%exparr_(delta_e) = exp(- beta * delta_e)
    end do
  end subroutine set_beta_ising2d

  !> update_ising2d: Update the system by Metropolis method.
  impure subroutine update_ising2d(this)
    class(ising2d), intent(inout) :: this
    real(real64), allocatable :: r(:)
    integer(int64) :: i, j
    allocate(r(this%nall_))
    call random_number(r)
    do j = 1, 2
       do i = j, this%nall_, 2
          call this%update_onesite(i, r(i))
       end do
       call this%update_norishiro()
    end do
  end subroutine update_ising2d
  !> update_onesite_ising2d: Update a spin of the system.
  pure subroutine update_onesite_ising2d(this, idx, r)
    class(ising2d), intent(inout) :: this
    integer(int64), intent(in) :: idx
    real(real64), intent(in) :: r
    integer(int32) :: delta_e
    delta_e = 2 * this%spins_(idx) * (this%spins_(idx + 1) + this%spins_(idx + this%nx_) + &
         & this%spins_(idx - 1) + this%spins_(idx - this%nx_))
    if (r < this%exparr_(delta_e)) &
         this%spins_(idx) = - this%spins_(idx)
  end subroutine update_onesite_ising2d
  !> update_norishiro_ising2d: Update norishiro.
  pure subroutine update_norishiro_ising2d(this)
    class(ising2d), intent(inout) :: this
    integer(int64) :: i
    do i = 1_int64, this%nx_
       this%spins_(this%norishiro_begin() + i - 1) = this%spins_(this%nall_ - this%nx_ + i)
       this%spins_(this%norishiro_end() - this%nx_ + i) = this%spins_(i)
    end do
  end subroutine update_norishiro_ising2d

  !> calc_total_energy_ising2d: Calculate the total energy.
  pure integer(int64) function calc_total_energy_ising2d(this) result(res)
    class(ising2d), intent(in) :: this
    integer(int64) :: i
    res = 0_int64
    do i = 1_int64, this%nall_
       res = res - this%spins_(i) * (this%spins_(i + 1) + this%spins_(i + this%nx_))
    end do
  end function calc_total_energy_ising2d
  !> calc_total_magne_ising2d: Calculate the total magne.
  pure integer(int64) function calc_total_magne_ising2d(this) result(res)
    class(ising2d), intent(in) :: this
    integer(int64) :: i
    res = 0_int64
    do i = 1_int64, this%nall_
       res = res + this%spins_(i)
    end do
  end function calc_total_magne_ising2d

  !> nx_ising2d: Return size of `x` of the system.
  pure integer(int64) function nx_ising2d(this) result(res)
    class(ising2d), intent(in) :: this
    res = this%nx_
  end function nx_ising2d
  !> ny_ising2d: Return size of `y` of the system.
  pure integer(int64) function ny_ising2d(this) result(res)
    class(ising2d), intent(in) :: this
    res = this%ny_
  end function ny_ising2d
  !> nall_ising2d: Return size of the system.
  pure integer(int64) function nall_ising2d(this) result(res)
    class(ising2d), intent(in) :: this
    res = this%nall_
  end function nall_ising2d
  !> kbt_ising2d: Return temperature of the system.
  pure real(real64) function kbt_ising2d(this) result(res)
    class(ising2d), intent(in) :: this
    res = 1 / this%beta_
  end function kbt_ising2d
  !> beta_ising2d: Return inverse temperature of the system.
  pure real(real64) function beta_ising2d(this) result(res)
    class(ising2d), intent(in) :: this
    res = this%beta_
  end function beta_ising2d
  !> norishiro_begin_ising2d: Return start index of `this%spins_(:)`.
  pure integer(int64) function norishiro_begin_ising2d(this) result(res)
    class(ising2d), intent(in) :: this
    res = 1 - this%nx_
  end function norishiro_begin_ising2d
  !> norishiro_end_ising2d: Return end index of `this%spins_(:)`.
  pure integer(int64) function norishiro_end_ising2d(this) result(res)
    class(ising2d), intent(in) :: this
    res = this%nall_ + this%nx_
  end function norishiro_end_ising2d
end module ising2d_m
