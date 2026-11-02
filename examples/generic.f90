module queue_callbacks
  use ukmc
  implicit none
  type :: queue_t
    integer(i8), allocatable :: job_ids(:)
    integer(i8) :: completed = 0_i8, next_id = 0_i8
  end type
contains
  function rates(data) result(a)
    class(*), intent(in) :: data
    real(dp), allocatable :: a(:)
    select type (data)
    type is (queue_t)
      a = [2.0_dp, 0.0_dp]
      if (size(data%job_ids) > 0) a(2) = 3.0_dp
    class default
      error stop 'wrong queue type'
    end select
  end function
  subroutine action(data, event)
    class(*), intent(inout) :: data
    integer(i8), intent(in) :: event
    select type (data)
    type is (queue_t)
      select case (event)
      case (0)
        data%job_ids = [data%job_ids, data%next_id]
        data%next_id = data%next_id + 1_i8
      case (1)
        call require(size(data%job_ids) > 0, 'empty queue')
        data%job_ids = data%job_ids(2:)
        data%completed = data%completed + 1_i8
      case default
        error stop 'invalid queue event'
      end select
    class default
      error stop 'wrong queue type'
    end select
  end subroutine
  subroutine observe(row)
    type(record_t), intent(in) :: row
    call require(row%time >= 0, 'invalid time')
  end subroutine
end module
program generic_example
  use ukmc
  use ukmc_events
  use queue_callbacks
  implicit none
  type(callback_model_t) :: model
  type(engine_t) :: engine
  type(queue_t) :: queue
  allocate(queue%job_ids(0))
  allocate(model%data, source=queue)
  model%rates => rates
  model%action => action
  call engine%seed(42_i8)
  call engine%run(model, 10000_i8, 100.0_dp, observe)
  select type (data => model%data)
  type is (queue_t)
    print *, size(data%job_ids), data%completed, engine%steps
  end select
end program
