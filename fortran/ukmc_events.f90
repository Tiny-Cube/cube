module ukmc_events
  use ukmc
  implicit none
  private
  public :: callback_model_t
  ! No prescribed state layout. Data can be any user-defined derived type.
  ! The propensity callback may generate any number of events; their indexing
  ! must agree with the action callback until that event has been applied.
  type, extends(model_t) :: callback_model_t
    class(*), allocatable :: data
    procedure(rates_callback), pointer, nopass :: rates => null()
    procedure(action_callback), pointer, nopass :: action => null()
    procedure(encode_callback), pointer, nopass :: encode => null()
  contains
    procedure :: propensities => callback_propensities
    procedure :: fire => callback_fire
    procedure :: state => callback_state
  end type
  abstract interface
    function rates_callback(data) result(rates)
      import :: dp
      class(*), intent(in) :: data
      real(dp), allocatable :: rates(:)
    end function
    subroutine action_callback(data, event)
      import :: i8
      class(*), intent(inout) :: data
      integer(i8), intent(in) :: event
    end subroutine
    function encode_callback(data) result(state)
      import :: i8
      class(*), intent(in) :: data
      integer(i8), allocatable :: state(:)
    end function
  end interface
contains
  function callback_propensities(self) result(rates)
    class(callback_model_t), intent(in) :: self
    real(dp), allocatable :: rates(:)
    call require(allocated(self%data), 'callback data not initialized')
    call require(associated(self%rates), 'missing rates callback')
    rates = self%rates(self%data)
  end function
  subroutine callback_fire(self, event)
    class(callback_model_t), intent(inout) :: self
    integer(i8), intent(in) :: event
    call require(allocated(self%data), 'callback data not initialized')
    call require(associated(self%action), 'missing action callback')
    call self%action(self%data, event)
  end subroutine
  function callback_state(self) result(state)
    class(callback_model_t), intent(in) :: self
    integer(i8), allocatable :: state(:)
    if (associated(self%encode)) then
      call require(allocated(self%data), 'callback data not initialized')
      state = self%encode(self%data)
    else
      allocate(state(0))
    end if
  end function
end module
