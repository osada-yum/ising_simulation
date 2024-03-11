module ising3d_m
  use, intrinsic :: iso_fortran_env
  implicit none
  private
  public :: ising3d
  integer(int32), parameter :: spin_map(0:1) = [-1, 1]
  type :: ising3d
     private
     integer(int64) :: nx_, ny_, nz_, nxy_, nall_
     real(real64) :: beta_
     integer(int8), allocatable :: spins_(:)
     integer(int64), allocatable :: energy_table_(:, :, :, :, :, :, :)
     real(real64), allocatable :: ws_(:, :, :, :, :, :, :)
   contains
     !> initializer.
     procedure, pass :: init => init_ising3d
     !> setter.
     procedure, pass :: set_ising_allup => set_ising_allup_ising3d
     procedure, pass :: set_ising_random => set_ising_random_ising3d
     procedure, pass :: set_kbt => set_kbt_ising3d
     procedure, pass :: set_beta => set_beta_ising3d
     !> updater.
     procedure, pass :: update => update_ising3d
     procedure, pass, private :: update_onesite => update_onesite_ising3d
     procedure, pass, private :: update_norishiro => update_norishiro_ising3d
     procedure, pass, private :: update_energy_table => update_energy_table_ising3d
     procedure, pass, private :: update_ws => update_ws_ising3d
     !> calculator.
     procedure, pass :: calc_energy_summ => calc_total_energy_ising3d
     procedure, pass :: calc_magne_summ => calc_total_magne_ising3d
     !> getter.
     procedure, pass :: nx => nx_ising3d
     procedure, pass :: ny => ny_ising3d
     procedure, pass :: nz => nz_ising3d
     procedure, pass :: nall => nall_ising3d
     procedure, pass :: kbt => kbt_ising3d
     procedure, pass :: beta => beta_ising3d
     procedure, pass, private :: norishiro_begin => norishiro_begin_ising3d
     procedure, pass, private :: norishiro_end => norishiro_end_ising3d
  end type ising3d
contains
  !> init_ising3d: Initialize ising3d once.
  pure subroutine init_ising3d(this, nx, ny, nz, kbt)
    class(ising3d), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny, nz
    real(real64), intent(in) :: kbt
    if (is_even(nx) .or. is_even(ny) .or. (.not. is_even(nz))) then
       error stop "The parity of size must be (x, y, z) == (odd, odd, even)."
    end if
    this%nx_ = nx
    this%ny_ = ny
    this%nz_ = nz
    this%nxy_ = nx * ny
    this%nall_ = nx * ny * nz
    allocate(this%spins_(this%norishiro_begin() : this%norishiro_end()))
    call this%set_ising_allup()
    call this%set_kbt(kbt)
  contains
    pure logical function is_even(v) result(res)
      integer(int64), intent(in) :: v
      res = iand(v, b'1') == 0_int64
    end function is_even
  end subroutine init_ising3d
  !> set_ising_allup_ising3d: Set spins `1`.
  pure subroutine set_ising_allup_ising3d(this)
    class(ising3d), intent(inout) :: this
    this%spins_(:) = 1_int32
  end subroutine set_ising_allup_ising3d
  !> set_ising_random_ising3d: Set spins `1` or `0` randomly.
  impure subroutine set_ising_random_ising3d(this)
    class(ising3d), intent(inout) :: this
    real(real64), allocatable :: r(:)
    allocate(r(1:this%nall_))
    call random_number(r)
    this%spins_(1:this%nall_) = merge(1_int32, 0_int32, r < 0.5_real64)
    call this%update_norishiro()
  end subroutine set_ising_random_ising3d
  !> set_kbt_ising3d: Set parameter `beta` as `1 / kbt`.
  pure subroutine set_kbt_ising3d(this, kbt)
    class(ising3d), intent(inout) :: this
    real(real64), intent(in) :: kbt
    call this%set_beta(1 / kbt)
  end subroutine set_kbt_ising3d
  !> set_beta_ising3d: Set parameter `beta` and update `this%exparr_`.
  pure subroutine set_beta_ising3d(this, beta)
    class(ising3d), intent(inout) :: this
    real(real64), intent(in) :: beta
    integer(int32) :: delta_e
    this%beta_ = beta
    call this%update_energy_table()
    call this%update_ws()
  end subroutine set_beta_ising3d
  pure subroutine update_energy_table_ising3d(this)
    class(ising3d), intent(inout) :: this
    integer(int32) :: u, d, l, r, f, b, center
    if (.not. allocated(this%energy_table_)) &
         & allocate(this%energy_table_(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1))
    do center = 0, 1
       do b = 0, 1
          do f = 0, 1
             do r = 0, 1
                do l = 0, 1
                   do d = 0, 1
                      do u = 0, 1
                         this%energy_table_(u, d, l, r, f, b, center) = &
                              & - spin_map(center) * &
                              & sum(spin_map([u, d, l, r, f, b]))
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine update_energy_table_ising3d
  pure subroutine update_ws_ising3d(this)
    class(ising3d), intent(inout) :: this
    integer(int32) :: u, d, l, r, f, b, center
    if (.not. allocated(this%ws_)) &
         & allocate(this%ws_(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1))
    do center = 0, 1
       do b = 0, 1
          do f = 0, 1
             do r = 0, 1
                do l = 0, 1
                   do d = 0, 1
                      do u = 0, 1
                         associate(delta_e => - 2 * this%energy_table_(u, d, l, r, f, b, center))
                           if (delta_e <= 0.0_real64) then
                              this%ws_(u, d, l, r, f, b, center) = 1.0_real64
                           else
                              this%ws_(u, d, l, r, f, b, center) = exp(- this%beta_ * delta_e)
                           end if
                         end associate
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine update_ws_ising3d
  !> update_ising3d: Update the system by Metropolis method.
  impure subroutine update_ising3d(this)
    class(ising3d), intent(inout) :: this
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
  end subroutine update_ising3d
  !> update_onesite_ising3d: Update a spin of the system.
  pure subroutine update_onesite_ising3d(this, idx, r)
    class(ising3d), intent(inout) :: this
    integer(int64), intent(in) :: idx
    real(real64), intent(in) :: r
    if (r < this%ws_(&
         & this%spins_(idx - 1), this%spins_(idx + 1), &
         & this%spins_(idx - this%nx_), this%spins_(idx + this%nx_), &
         & this%spins_(idx - this%nxy_), this%spins_(idx + this%nxy_), &
         & this%spins_(idx))) then
       this%spins_(idx) = 1_int8 - this%spins_(idx)
    end if
  end subroutine update_onesite_ising3d
  !> update_norishiro_ising3d: Update norishiro.
  pure subroutine update_norishiro_ising3d(this)
    class(ising3d), intent(inout) :: this
    integer(int64) :: i
    do i = 1_int64, this%nxy_
       this%spins_(this%norishiro_begin() + i - 1) = this%spins_(this%nall_ - this%nxy_ + i)
       this%spins_(this%norishiro_end() - this%nxy_ + i) = this%spins_(i)
    end do
  end subroutine update_norishiro_ising3d

  !> calc_total_energy_ising3d: Calculate the total energy.
  pure integer(int64) function calc_total_energy_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    integer(int64) :: i
    res = 0_int64
    do i = 1_int64, this%nall_
       res = res + this%energy_table_(&
            & this%spins_(i - 1), this%spins_(i + 1), &
            & this%spins_(i - this%nx_), this%spins_(i + this%nx_), &
            & this%spins_(i - this%nxy_), this%spins_(i + this%nxy_), &
            & this%spins_(i))
    end do
    res = res / 2
  end function calc_total_energy_ising3d
  !> calc_total_magne_ising3d: Calculate the total magne.
  pure integer(int64) function calc_total_magne_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    integer(int64) :: i
    res = 0_int64
    do i = 1_int64, this%nall_
       res = res + spin_map(this%spins_(i))
    end do
  end function calc_total_magne_ising3d

  !> nx_ising3d: Return size of `x` of the system.
  pure integer(int64) function nx_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    res = this%nx_
  end function nx_ising3d
  !> ny_ising3d: Return size of `y` of the system.
  pure integer(int64) function ny_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    res = this%ny_
  end function ny_ising3d
  !> nz_ising3d: Return size of `z` of the system.
  pure integer(int64) function nz_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    res = this%nz_
  end function nz_ising3d
  !> nall_ising3d: Return size of the system.
  pure integer(int64) function nall_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    res = this%nall_
  end function nall_ising3d
  !> kbt_ising3d: Return temperature of the system.
  pure real(real64) function kbt_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    res = 1 / this%beta_
  end function kbt_ising3d
  !> beta_ising3d: Return inverse temperature of the system.
  pure real(real64) function beta_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    res = this%beta_
  end function beta_ising3d
  !> norishiro_begin_ising3d: Return start index of `this%spins_(:)`.
  pure integer(int64) function norishiro_begin_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    res = 1 - this%nxy_
  end function norishiro_begin_ising3d
  !> norishiro_end_ising3d: Return end index of `this%spins_(:)`.
  pure integer(int64) function norishiro_end_ising3d(this) result(res)
    class(ising3d), intent(in) :: this
    res = this%nall_ + this%nxy_
  end function norishiro_end_ising3d
end module ising3d_m
