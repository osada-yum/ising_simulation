module ising3d_m
  use, intrinsic :: iso_fortran_env
  use union_find_m
  implicit none
  private
  integer(int32), parameter :: spin_map(0:1) = [-1, 1]
  public :: ising3d
  type :: ising3d
     private
     integer(int64) :: nx_, ny_, nz_, nxy_, nall_
     real(real64) :: beta_
     integer(int32), allocatable :: spins_(:)
     integer(int64), allocatable :: energy_table_(:, :)
     real(real64), allocatable :: ws_(:, :)
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
     procedure, pass :: update_swendsen_wang => update_swendsen_wang_ising3d
     procedure, pass :: flip_all => flip_all_ising3d
     procedure, pass, private :: update_onesite => update_onesite_ising3d
     procedure, pass, private :: update_norishiro => update_norishiro_ising3d
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
  !> set_ising_random_ising3d: Set spins `1` or `-1` randomly.
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
    integer(int32) :: i1, i2, i3, i4, i5, i6
    this%beta_ = beta
    if (.not. allocated(this%energy_table_)) &
         & allocate(this%energy_table_(0:3, 0:1))
    if (.not. allocated(this%ws_)) &
         & allocate(this%ws_(0:6, 0:1))
    do i1 = 0, 1
       do i2 = 0, 1
          do i3 = 0, 1
             associate(s => i1 + i2 + i3)
               this%energy_table_(s, 0) = - spin_map(0) * sum(spin_map([i1, i2, i3]))
               this%energy_table_(s, 1) = - spin_map(1) * sum(spin_map([i1, i2, i3]))
           end associate
          end do
       end do
    end do
    do i1 = 0, 1
       do i2 = 0, 1
          do i3 = 0, 1
             associate(s1 => i1 + i2 + i3)
               do i4 = 0, 1
                  do i5 = 0, 1
                     do i6 = 0, 1
                        associate(s2 => i4 + i5 + i6)
                          associate(e1 => this%energy_table_(s1, 0) + this%energy_table_(s2, 0), &
                               & e2 => this%energy_table_(s1, 1) + this%energy_table_(s2, 1))
                            this%ws_(s1 + s2, 0) = min(1.0_real64, exp(- this%beta_ * (e2 - e1)))
                            this%ws_(s1 + s2, 1) = min(1.0_real64, exp(- this%beta_ * (e1 - e2)))
                          end associate
                        end associate
                     end do
                  end do
               end do
             end associate
          end do
       end do
    end do
  end subroutine set_beta_ising3d

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
  !> update_swendsen_wang_ising3d: Update the system by Metropolis method.
  impure subroutine update_swendsen_wang_ising3d(this)
    class(ising3d), intent(inout) :: this
    real(real64), allocatable :: r(:, :)
    type(union_find) :: uf
    logical, allocatable :: is_flip(:)
    integer(int64) :: i
    allocate(r(4, this%nall_))
    allocate(is_flip(this%nall_), source = .false.)
    call random_number(r)
    call uf%init(this%nall_)
    do i = 1, this%nall_
       block
         integer(int64) :: right, up, top
         right = i + 1; up = i + this%nx_; top = i + this%nxy_
         if (right > this%nall_) right = 1
         if (up > this%nall_) up = up - this%nall_
         if (top > this%nall_) top = top - this%nall_
         call maybe_connect(uf, i, right, r(1, i))
         call maybe_connect(uf, i, up, r(2, i))
         call maybe_connect(uf, i, top, r(3, i))
         if (r(4, i) < 0.5_real64) &
              & is_flip(i) = .true.
       end block
    end do
    do i = 1, this%nall_
       if (is_flip(uf%root(i))) &
            & this%spins_(i) = 1 - this%spins_(i)
    end do
    call this%update_norishiro()
  contains
    impure subroutine maybe_connect(uf, i1, i2, r)
      type(union_find), intent(inout) :: uf
      integer(int64), intent(in) :: i1, i2
      real(real64), intent(in) :: r
      if (this%spins_(i1) /= this%spins_(i2)) return
      if (uf%same(i1, i2)) return
      if (r >= 1 - exp(-2 * this%beta_)) return
      call uf%union(i1, i2)
    end subroutine maybe_connect
  end subroutine update_swendsen_wang_ising3d
  pure subroutine flip_all_ising3d(this)
    class(ising3d), intent(inout) :: this
    this%spins_(:) = 1 - this%spins_(:)
  end subroutine flip_all_ising3d
  !> update_onesite_ising3d: Update a spin of the system.
  pure subroutine update_onesite_ising3d(this, idx, r)
    class(ising3d), intent(inout) :: this
    integer(int64), intent(in) :: idx
    real(real64), intent(in) :: r
    associate(s => this%spins_(idx + 1) + this%spins_(idx - 1) + &
         & this%spins_(idx + this%nx_) + this%spins_(idx - this%nx_) + &
         & this%spins_(idx + this%nxy_) + this%spins_(idx - this%nxy_))
      if (r < this%ws_(s, this%spins_(idx))) &
           & this%spins_(idx) = 1 - this%spins_(idx)
    end associate
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
       associate(s => this%spins_(i + 1) + this%spins_(i + this%nx_) + this%spins_(i + this%nxy_))
         res = res + this%energy_table_(s, this%spins_(i))
       end associate
    end do
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
