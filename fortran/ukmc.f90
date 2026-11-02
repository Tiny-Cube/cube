module ukmc
  use iso_fortran_env, only: int64, real64, error_unit
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: i8, dp, max_count, model_t, graph_t, reaction_network_t, edge_t, reaction_t
  public :: engine_t, record_t, rng_t, require, valid_count, nonnegative
  public :: init_graph, init_reactions, marcus, arrhenius, miller_abrahams, infinity
  integer, parameter :: i8 = int64, dp = real64
  integer(i8), parameter :: max_count = 9007199254740991_i8
  real(dp), parameter :: kb = 8.617333262145e-5_dp, hbar = 6.582119569e-16_dp
  real(dp), parameter :: pi = 3.14159265358979323846_dp

  type, abstract :: model_t
  contains
    procedure(propensity_interface), deferred :: propensities
    procedure(fire_interface), deferred :: fire
    procedure :: state => default_state
    procedure :: position => default_position
  end type
  type :: rng_t
    integer(i8) :: state = 2_i8
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
  end type
  type :: record_t
    integer(i8) :: step = 0_i8, event = -1_i8
    real(dp) :: time = 0.0_dp, xyz(3) = 0.0_dp
    integer(i8), allocatable :: state(:)
  end type
  type :: engine_t
    type(rng_t) :: rng
    real(dp) :: time = 0.0_dp
    integer(i8) :: steps = 0_i8
    character(16) :: status = 'ready'
    procedure(uniform_interface), pointer, nopass :: random => null()
  contains
    procedure :: seed => engine_seed
    procedure :: uniform => engine_uniform
    procedure :: snapshot => engine_snapshot
    procedure :: step => engine_step
    procedure :: run => engine_run
  end type
  type :: edge_t
    integer(i8) :: source = 0_i8, target = 0_i8
    real(dp) :: rate = 0.0_dp, displacement(3) = 0.0_dp
  end type
  type, extends(model_t) :: graph_t
    private
    type(edge_t), allocatable :: edges(:)
    integer(i8) :: site = 0_i8
    real(dp) :: xyz(3) = 0.0_dp
  contains
    procedure :: propensities => graph_propensities
    procedure :: fire => graph_fire
    procedure :: state => graph_state
    procedure :: position => graph_position
  end type
  type :: reaction_t
    real(dp) :: rate = 0.0_dp
    integer(i8), allocatable :: reactants(:), products(:)
  end type
  type, extends(model_t) :: reaction_network_t
    private
    integer(i8), allocatable :: populations(:)
    type(reaction_t), allocatable :: reactions(:)
  contains
    procedure :: propensities => reaction_propensities
    procedure :: fire => reaction_fire
    procedure :: state => reaction_state
    procedure :: position => reaction_position
  end type
  abstract interface
    function propensity_interface(self) result(rates)
      import :: model_t, dp
      class(model_t), intent(in) :: self
      real(dp), allocatable :: rates(:)
    end function
    subroutine fire_interface(self, event)
      import :: model_t, i8
      class(model_t), intent(inout) :: self
      integer(i8), intent(in) :: event
    end subroutine
    function state_interface(self) result(state)
      import :: model_t, i8
      class(model_t), intent(in) :: self
      integer(i8), allocatable :: state(:)
    end function
    function position_interface(self) result(xyz)
      import :: model_t, dp
      class(model_t), intent(in) :: self
      real(dp) :: xyz(3)
    end function
    function uniform_interface() result(u)
      import :: dp
      real(dp) :: u
    end function
    subroutine observer_interface(row)
      import :: record_t
      type(record_t), intent(in) :: row
    end subroutine
  end interface
contains
  function default_state(self) result(state)
    class(model_t), intent(in) :: self
    integer(i8), allocatable :: state(:)
    allocate(state(0))
    ! self is intentionally unused: these are optional diagnostic observations.
  end function
  function default_position(self) result(xyz)
    class(model_t), intent(in) :: self
    real(dp) :: xyz(3)
    xyz = 0.0_dp
  end function
  subroutine require(ok, message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not. ok) then
      write(error_unit, '(a)') 'ukmc: '//message
      error stop 2
    end if
  end subroutine
  subroutine valid_count(n)
    integer(i8), intent(in) :: n
    call require(n >= 0_i8 .and. n <= max_count, 'count outside [0,2^53-1]')
  end subroutine
  subroutine nonnegative(a)
    real(dp), intent(in) :: a
    call require(ieee_is_finite(a) .and. a >= 0.0_dp, 'rate must be finite and nonnegative')
  end subroutine
  function infinity() result(x)
    real(dp) :: x
    x = ieee_value(0.0_dp, ieee_positive_inf)
  end function
  subroutine rng_seed(self, seed)
    class(rng_t), intent(inout) :: self
    integer(i8), intent(in) :: seed
    call valid_count(seed)
    self%state = modulo(seed, 2147483646_i8) + 1_i8
  end subroutine
  function rng_uniform(self) result(u)
    class(rng_t), intent(inout) :: self
    real(dp) :: u
    self%state = modulo(48271_i8 * self%state, 2147483647_i8)
    u = real(self%state, dp) / 2147483647.0_dp
  end function
  subroutine engine_seed(self, seed)
    class(engine_t), intent(inout) :: self
    integer(i8), intent(in) :: seed
    call self%rng%seed(seed)
  end subroutine
  function engine_uniform(self) result(u)
    class(engine_t), intent(inout) :: self
    real(dp) :: u
    if (associated(self%random)) then
      u = self%random()
    else
      u = self%rng%uniform()
    end if
    call require(ieee_is_finite(u) .and. u > 0.0_dp .and. u < 1.0_dp, 'RNG must return a value in (0,1)')
  end function
  function engine_snapshot(self, model, event) result(row)
    class(engine_t), intent(in) :: self
    class(model_t), intent(in) :: model
    integer(i8), intent(in), optional :: event
    type(record_t) :: row
    row%step = self%steps
    row%time = self%time
    if (present(event)) row%event = event
    row%xyz = model%position()
    row%state = model%state()
  end function
  logical function engine_step(self, model, until, event) result(advanced)
    class(engine_t), intent(inout) :: self
    class(model_t), intent(inout) :: model
    real(dp), intent(in) :: until
    integer(i8), intent(out), optional :: event
    real(dp), allocatable :: rates(:)
    real(dp) :: total, dt, threshold, cumulative
    integer :: j, chosen
    call require(.not. ieee_is_nan(until) .and. until >= self%time, 'until must be >= current time')
    advanced = .false.
    if (present(event)) event = -1_i8
    if (self%time == until) then
      self%status = 'time_limit'
      return
    end if
    rates = model%propensities()
    total = 0.0_dp
    do j = 1, size(rates)
      call nonnegative(rates(j))
      total = total + rates(j)
    end do
    call require(ieee_is_finite(total), 'total propensity overflow')
    if (total == 0.0_dp) then
      self%status = 'absorbed'
      if (ieee_is_finite(until)) self%time = until
      return
    end if
    dt = -log(self%uniform()) / total
    if (dt > until - self%time) then
      self%time = until
      self%status = 'time_limit'
      return
    end if
    call require(ieee_is_finite(dt) .and. &
      ieee_is_finite(self%time + dt) .and. self%time + dt > self%time, 'waiting time outside floating-point range')
    threshold = self%uniform() * total
    cumulative = 0.0_dp
    chosen = 1
    do j = 1, size(rates)
      if (rates(j) > 0.0_dp) chosen = j
      cumulative = cumulative + rates(j)
      if (threshold < cumulative) exit
    end do
    ! Public event IDs are zero-based, including the Fortran API.
    call model%fire(int(chosen - 1, i8))
    self%time = self%time + dt
    self%steps = self%steps + 1_i8
    self%status = 'running'
    if (present(event)) event = int(chosen - 1, i8)
    advanced = .true.
  end function
  subroutine engine_run(self, model, max_steps, until, observe)
    class(engine_t), intent(inout) :: self
    class(model_t), intent(inout) :: model
    integer(i8), intent(in) :: max_steps
    real(dp), intent(in) :: until
    procedure(observer_interface) :: observe
    integer(i8) :: j, event
    real(dp) :: before
    call valid_count(max_steps)
    call require(.not. ieee_is_nan(until) .and. until >= self%time, 'until must be >= current time')
    call observe(self%snapshot(model))
    do j = 1_i8, max_steps
      before = self%time
      if (.not. self%step(model, until, event)) then
        if (self%time /= before) call observe(self%snapshot(model))
        return
      end if
      call observe(self%snapshot(model, event))
    end do
    if (self%time == until) then
      self%status = 'time_limit'
    else
      self%status = 'step_limit'
    end if
  end subroutine
  subroutine init_graph(model, sites, edges, initial)
    type(graph_t), intent(out) :: model
    integer(i8), intent(in) :: sites, initial
    type(edge_t), intent(in) :: edges(:)
    integer :: j
    call valid_count(sites)
    call valid_count(initial)
    call require(sites > 0 .and. initial < sites, 'invalid graph size or initial site')
    do j = 1, size(edges)
      call valid_count(edges(j)%source)
      call valid_count(edges(j)%target)
      call require(edges(j)%source < sites .and. edges(j)%target < sites, 'edge endpoint out of bounds')
      call nonnegative(edges(j)%rate)
      call require(all(ieee_is_finite(edges(j)%displacement)), 'non-finite displacement')
    end do
    model%edges = edges
    model%site = initial
  end subroutine
  function graph_propensities(self) result(rates)
    class(graph_t), intent(in) :: self
    real(dp), allocatable :: rates(:)
    integer :: j
    call require(allocated(self%edges), 'graph not initialized')
    allocate(rates(size(self%edges)), source=0.0_dp)
    do j = 1, size(rates)
      if (self%edges(j)%source == self%site) rates(j) = self%edges(j)%rate
    end do
  end function
  subroutine graph_fire(self, event)
    class(graph_t), intent(inout) :: self
    integer(i8), intent(in) :: event
    integer :: j
    real(dp) :: next(3)
    call require(allocated(self%edges), 'graph not initialized')
    call require(event >= 0 .and. event < size(self%edges, kind=i8), 'event out of bounds')
    j = int(event) + 1
    call require(self%edges(j)%source == self%site .and. self%edges(j)%rate > 0, 'inactive event')
    next = self%xyz + self%edges(j)%displacement
    call require(all(ieee_is_finite(next)), 'position overflow')
    self%xyz = next
    self%site = self%edges(j)%target
  end subroutine
  function graph_state(self) result(state)
    class(graph_t), intent(in) :: self
    integer(i8), allocatable :: state(:)
    state = [self%site]
  end function
  function graph_position(self) result(xyz)
    class(graph_t), intent(in) :: self
    real(dp) :: xyz(3)
    xyz = self%xyz
  end function
  subroutine init_reactions(model, populations, reactions)
    type(reaction_network_t), intent(out) :: model
    integer(i8), intent(in) :: populations(:)
    type(reaction_t), intent(in) :: reactions(:)
    integer :: i, j
    call require(size(populations) > 0, 'at least one species required')
    do i = 1, size(populations)
      call valid_count(populations(i))
    end do
    do j = 1, size(reactions)
      call nonnegative(reactions(j)%rate)
      call require(allocated(reactions(j)%reactants) .and. allocated(reactions(j)%products), 'unallocated stoichiometry')
      call require(size(reactions(j)%reactants) == size(populations), 'reactants dimension mismatch')
      call require(size(reactions(j)%products) == size(populations), 'products dimension mismatch')
      do i = 1, size(populations)
        call valid_count(reactions(j)%reactants(i))
        call valid_count(reactions(j)%products(i))
      end do
    end do
    model%populations = populations
    model%reactions = reactions
  end subroutine
  function reaction_propensities(self) result(rates)
    class(reaction_network_t), intent(in) :: self
    real(dp), allocatable :: rates(:)
    real(dp) :: a
    integer :: j, i
    integer(i8) :: n, order, k
    call require(allocated(self%reactions), 'reaction network not initialized')
    allocate(rates(size(self%reactions)), source=0.0_dp)
    do j = 1, size(rates)
      if (self%reactions(j)%rate == 0.0_dp) cycle
      if (any(self%populations < self%reactions(j)%reactants)) cycle
      a = self%reactions(j)%rate
      do i = 1, size(self%populations)
        n = self%populations(i)
        order = min(self%reactions(j)%reactants(i), n - self%reactions(j)%reactants(i))
        do k = 1_i8, order
          a = a * (real(n - k + 1_i8, dp) / real(k, dp))
          call nonnegative(a)
        end do
      end do
      rates(j) = a
    end do
  end function
  subroutine reaction_fire(self, event)
    class(reaction_network_t), intent(inout) :: self
    integer(i8), intent(in) :: event
    integer :: j, i
    integer(i8), allocatable :: next(:)
    call require(allocated(self%reactions), 'reaction network not initialized')
    call require(event >= 0 .and. event < size(self%reactions, kind=i8), 'event out of bounds')
    j = int(event) + 1
    call require(self%reactions(j)%rate > 0, 'inactive reaction')
    call require(all(self%populations >= self%reactions(j)%reactants), 'inactive reaction')
    next = self%populations - self%reactions(j)%reactants + self%reactions(j)%products
    do i = 1, size(next)
      call valid_count(next(i))
    end do
    self%populations = next
  end subroutine
  function reaction_state(self) result(state)
    class(reaction_network_t), intent(in) :: self
    integer(i8), allocatable :: state(:)
    call require(allocated(self%populations), 'reaction network not initialized')
    state = self%populations
  end function
  function reaction_position(self) result(xyz)
    class(reaction_network_t), intent(in) :: self
    real(dp) :: xyz(3)
    call require(allocated(self%populations), 'reaction network not initialized')
    xyz = 0.0_dp
  end function
  function thermal(t) result(kt)
    real(dp), intent(in) :: t
    real(dp) :: kt
    call require(ieee_is_finite(t) .and. t > 0, 'positive finite temperature required')
    kt = kb*t
  end function
  function marcus(coupling, lambda, dg, t) result(rate)
    real(dp), intent(in) :: coupling, lambda, dg, t
    real(dp) :: rate, kt, x
    kt = thermal(t)
    call require(ieee_is_finite(coupling) .and. ieee_is_finite(lambda), 'non-finite Marcus energies')
    call require(lambda > 0 .and. ieee_is_finite(dg), 'invalid Marcus energies')
    if (coupling == 0.0_dp) then
      rate = 0.0_dp
      return
    end if
    x = (dg + lambda) / sqrt(lambda*kt)
    rate = exp(log(2*pi/hbar) + 2*log(abs(coupling)) - 0.5_dp*log(4*pi*lambda*kt) - x*x/4)
    call nonnegative(rate)
  end function
  function arrhenius(attempt, barrier, t) result(rate)
    real(dp), intent(in) :: attempt, barrier, t
    real(dp) :: rate
    call nonnegative(attempt)
    call nonnegative(barrier)
    rate = attempt * exp(-barrier / thermal(t))
  end function
  function miller_abrahams(attempt, distance, localization, dg, t) result(rate)
    real(dp), intent(in) :: attempt, distance, localization, dg, t
    real(dp) :: rate, kt
    kt = thermal(t)
    call nonnegative(attempt)
    call nonnegative(distance)
    call require(ieee_is_finite(localization) .and. localization > 0, 'invalid localization length')
    call require(ieee_is_finite(dg), 'non-finite delta_g')
    rate = attempt * exp(-2*distance/localization - max(dg, 0.0_dp)/kt)
  end function
end module
