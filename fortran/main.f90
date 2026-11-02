program main
  use iso_fortran_env, only: output_unit, error_unit
  use, intrinsic :: ieee_arithmetic
  use ukmc
  use ukmc_io
  implicit none
  class(model_t), allocatable :: model
  type(engine_t) :: engine
  integer(i8) :: steps, seed
  integer(i8), allocatable :: state(:)
  integer :: argc, i, length
  real(dp) :: until
  character(:), allocatable :: arg
  argc = command_argument_count()
  call require(argc >= 1 .and. argc <= 4, 'usage: ukmc-fortran MODEL [STEPS [SEED [UNTIL]]]')
  call argument(1, arg)
  call load_model(arg, model)
  steps = 1000_i8
  seed = 1_i8
  until = infinity()
  if (argc >= 2) then
    call argument(2, arg)
    steps = parse_integer(arg)
  end if
  if (argc >= 3) then
    call argument(3, arg)
    seed = parse_integer(arg)
  end if
  if (argc >= 4) then
    call argument(4, arg)
    until = parse_number(arg)
    call require(ieee_is_finite(until) .and. until >= 0.0_dp, 'invalid time limit')
  end if
  call engine%seed(seed)
  state = model%state()
  write(output_unit, '(a)', advance='no') 'step,time,event,x,y,z'
  do i = 1, size(state)
    write(output_unit, '(a,i0)', advance='no') ',state', i-1
  end do
  write(output_unit, *)
  call engine%run(model, steps, until, output)
  write(error_unit, '(a,a,a,i0,a,es25.17e3)') 'status=',trim(engine%status),' events=',engine%steps,' time=',engine%time
contains
  subroutine argument(index, value)
    integer, intent(in) :: index
    character(:), allocatable, intent(out) :: value
    call get_command_argument(index, length=length)
    allocate(character(length) :: value)
    call get_command_argument(index, value)
  end subroutine
  subroutine output(row)
    type(record_t), intent(in) :: row
    integer :: j
    write(output_unit, '(i0,a,es25.17e3,a,i0)', advance='no') row%step,',',row%time,',',row%event
    do j = 1, 3
      write(output_unit, '(a,es25.17e3)', advance='no') ',',row%xyz(j)
    end do
    do j = 1, size(row%state)
      write(output_unit, '(a,i0)', advance='no') ',',row%state(j)
    end do
    write(output_unit, *)
  end subroutine
end program
