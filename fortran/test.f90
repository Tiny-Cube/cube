module test_callbacks
  use ukmc
  implicit none
  type :: task_state_t
    integer, allocatable :: tasks(:)
    real(dp) :: weight = 1.25_dp
  end type
  integer :: index = 0, observations = 0
contains
  function dynamic_rates(data) result(rates)
    class(*), intent(in) :: data
    real(dp), allocatable :: rates(:)
    select type(data)
    type is(task_state_t)
      allocate(rates(size(data%tasks)),source=1.0_dp)
    class default
      error stop 'wrong task type'
    end select
  end function
  subroutine dynamic_action(data,event)
    class(*), intent(inout) :: data
    integer(i8), intent(in) :: event
    integer :: j
    select type(data)
    type is(task_state_t)
      j = int(event) + 1
      call require(j >= 1 .and. j <= size(data%tasks),'invalid task index')
      data%tasks = [data%tasks(:j-1), data%tasks(j+1:)]
      data%weight = data%weight * 0.5_dp
    class default
      error stop 'wrong task type'
    end select
  end subroutine
  subroutine ignore(row)
    type(record_t),intent(in) :: row
    call require(size(row%state) == 0,'optional observation should be empty')
  end subroutine
  function fixed() result(u)
    real(dp) :: u
    index = index + 1
    u = exp(-1.0_dp)
    if (index == 2) u = 0.8_dp
  end function
  subroutine observe(row)
    type(record_t), intent(in) :: row
    observations = observations + 1
    if (observations == 2) call require(row%event == 1, 'selection failed')
    if (observations == 3) call require(row%event == -1, 'censoring failed')
  end subroutine
  subroutine conserved(row)
    type(record_t), intent(in) :: row
    call require(row%state(1) + 2*row%state(2) == 10, 'conservation failed')
  end subroutine
end module
program test
  use ukmc
  use test_callbacks
  use ukmc_events
  implicit none
  type(graph_t) :: graph
  type(reaction_network_t) :: network
  type(reaction_t) :: reactions(2)
  type(edge_t) :: edges(2)
  type(engine_t) :: engine, re, custom_engine
  type(callback_model_t) :: custom
  type(task_state_t) :: tasks
  type(rng_t) :: rng
  real(dp), allocatable :: rates(:)
  tasks%tasks = [1,2,3]
  allocate(custom%data,source=tasks)
  custom%rates => dynamic_rates
  custom%action => dynamic_action
  call custom_engine%run(custom,10_i8,infinity(),ignore)
  call require(custom_engine%steps == 3 .and. trim(custom_engine%status) == 'absorbed','dynamic catalog failed')
  select type(data => custom%data)
  type is(task_state_t)
    call require(size(data%tasks) == 0 .and. abs(data%weight-.15625_dp)<1e-15_dp,'custom state failed')
  end select
  call rng%seed(0_i8)
  call require(abs(rng%uniform()-48271.0_dp/2147483647.0_dp) < 1e-20_dp, 'RNG failed')
  edges(1)%rate = 2
  edges(2)%rate = 3
  call init_graph(graph,1_i8,edges,0_i8)
  engine%random => fixed
  call engine%run(graph,10_i8,0.3_dp,observe)
  call require(observations == 3 .and. engine%steps == 1, 'horizon failed')
  reactions(1)%rate = 0.1_dp
  reactions(1)%reactants = [2_i8,0_i8]
  reactions(1)%products = [0_i8,1_i8]
  reactions(2)%rate = 0.3_dp
  reactions(2)%reactants = [0_i8,1_i8]
  reactions(2)%products = [2_i8,0_i8]
  call init_reactions(network,[10_i8,0_i8],reactions)
  rates = network%propensities()
  call require(abs(rates(1)-4.5_dp) < 1e-12_dp, 'mass action failed')
  call re%seed(99_i8)
  call re%run(network,2000_i8,infinity(),conserved)
  call require(abs(marcus(.01_dp,.2_dp,.05_dp,300._dp)/marcus(.01_dp,.2_dp,-.05_dp,300._dp) &
    - exp(-.05_dp/(8.617333262145e-5_dp*300._dp))) < 1e-13_dp, 'detailed balance failed')
  print '(3(es25.17e3,1x))',marcus(.01_dp,.2_dp,.05_dp,300._dp),arrhenius(1e12_dp,.1_dp,300._dp), &
    miller_abrahams(1e12_dp,1._dp,2._dp,.1_dp,300._dp)
end program
