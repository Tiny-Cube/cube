program benchmark
  use ukmc
  use ukmc_io
  implicit none
  class(model_t), allocatable :: model
  type(engine_t) :: engine
  integer(i8) :: steps, seed, checksum, started, finished, ticks
  character(4096) :: arg
  call get_command_argument(1, arg)
  call load_model(trim(arg), model)
  call get_command_argument(2, arg)
  steps = parse_integer(trim(arg))
  call get_command_argument(3, arg)
  seed = parse_integer(trim(arg))
  call engine%seed(seed)
  checksum = 0_i8
  call system_clock(started, ticks)
  call engine%run(model, steps, infinity(), observe)
  call system_clock(finished)
  write(*,'(*(g0,1x))') real(finished-started,dp)/real(ticks,dp), &
    engine%steps, engine%time, checksum, model%position(), model%state()
contains
  subroutine observe(row)
    type(record_t), intent(in) :: row
    checksum = checksum + row%event + 1_i8
  end subroutine
end program
