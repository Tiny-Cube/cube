module ukmc_io
  use iso_fortran_env, only: iostat_end, iostat_eor
  use ukmc
  implicit none
  private
  public :: load_model, parse_integer, parse_number
  type :: tokens_t
    character(:), allocatable :: text
    integer :: offset = 1
  contains
    procedure :: next => next_token
    procedure :: integer => next_integer
    procedure :: number => next_number
  end type
contains
  function next_token(self) result(word)
    class(tokens_t), intent(inout) :: self
    character(:), allocatable :: word
    integer :: first, n
    n = len(self%text)
    do while (self%offset <= n)
      if (iachar(self%text(self%offset:self%offset)) > 32) exit
      self%offset = self%offset + 1
    end do
    first = self%offset
    do while (self%offset <= n)
      if (iachar(self%text(self%offset:self%offset)) <= 32) exit
      self%offset = self%offset + 1
    end do
    word = self%text(first:self%offset - 1)
  end function
  function parse_integer(word) result(n)
    character(*), intent(in) :: word
    integer(i8) :: n
    integer :: j, status
    call require(len(word) > 0, 'expected integer, got end of input')
    do j = 1, len(word)
      call require(index('0123456789', word(j:j)) > 0, 'expected unsigned decimal integer')
    end do
    read(word, *, iostat=status) n
    call require(status == 0, 'integer overflow')
    call valid_count(n)
  end function
  function parse_number(word) result(x)
    character(*), intent(in) :: word
    real(dp) :: x
    integer :: j, n, digits, status
    n = len(word)
    call require(n > 0, 'expected number, got end of input')
    j = 1
    if (word(j:j) == '+' .or. word(j:j) == '-') j = j + 1
    digits = 0
    do while (j <= n)
      if (index('0123456789', word(j:j)) == 0) exit
      digits = digits + 1
      j = j + 1
    end do
    if (j <= n) then
      if (word(j:j) == '.') then
        j = j + 1
        do while (j <= n)
          if (index('0123456789', word(j:j)) == 0) exit
          digits = digits + 1
          j = j + 1
        end do
      end if
    end if
    call require(digits > 0, 'invalid decimal number')
    if (j <= n) then
      if (word(j:j) == 'e' .or. word(j:j) == 'E') then
        j = j + 1
        if (j <= n) then
          if (word(j:j) == '+' .or. word(j:j) == '-') j = j + 1
        end if
        digits = 0
        do while (j <= n)
          if (index('0123456789', word(j:j)) == 0) exit
          digits = digits + 1
          j = j + 1
        end do
        call require(digits > 0, 'missing exponent')
      end if
    end if
    call require(j == n + 1, 'invalid decimal number')
    read(word, *, iostat=status) x
    call require(status == 0, 'number outside supported range')
  end function
  function next_integer(self) result(n)
    class(tokens_t), intent(inout) :: self
    integer(i8) :: n
    n = parse_integer(self%next())
  end function
  function next_number(self) result(x)
    class(tokens_t), intent(inout) :: self
    real(dp) :: x
    x = parse_number(self%next())
  end function
  subroutine load_model(path, model)
    character(*), intent(in) :: path
    class(model_t), allocatable, intent(out) :: model
    type(tokens_t) :: t
    type(graph_t) :: graph
    type(reaction_network_t) :: network
    type(edge_t), allocatable :: edges(:)
    type(reaction_t), allocatable :: reactions(:)
    integer(i8), allocatable :: populations(:)
    character(4096) :: buffer
    character(:), allocatable :: line, kind, word
    integer :: unit, status, used, comment, i, j, k
    integer(i8) :: n, m, initial
    open(newunit=unit, file=path, status='old', action='read', iostat=status)
    call require(status == 0, 'cannot open model')
    t%text = ''
    line = ''
    do
      read(unit, '(a)', advance='no', size=used, iostat=status) buffer
      if (used > 0) line = line//buffer(:used)
      if (status == iostat_eor .or. status == iostat_end) then
        comment = index(line, '#')
        if (comment > 0) line = line(:comment - 1)
        t%text = t%text//line//' '
        line = ''
      end if
      if (status == iostat_end) exit
      call require(status == 0 .or. status == iostat_eor, 'model read failed')
    end do
    close(unit)
    word = t%next()
    call require(word == 'KMC1', 'expected KMC1 header')
    kind = t%next()
    n = t%integer()
    m = t%integer()
    call require(m <= huge(0), 'too many events for this Fortran build')
    select case (kind)
    case ('graph')
      initial = t%integer()
      allocate(edges(int(m)))
      do j = 1, int(m)
        edges(j)%source = t%integer()
        edges(j)%target = t%integer()
        edges(j)%rate = t%number()
        do k = 1, 3
          edges(j)%displacement(k) = t%number()
        end do
      end do
      call init_graph(graph, n, edges, initial)
      allocate(model, source=graph)
    case ('reaction')
      call require(n <= huge(0), 'too many species for this Fortran build')
      allocate(populations(int(n)), reactions(int(m)))
      do i = 1, int(n)
        populations(i) = t%integer()
      end do
      do j = 1, int(m)
        reactions(j)%rate = t%number()
        allocate(reactions(j)%reactants(int(n)), reactions(j)%products(int(n)))
        do i = 1, int(n)
          reactions(j)%reactants(i) = t%integer()
        end do
        do i = 1, int(n)
          reactions(j)%products(i) = t%integer()
        end do
      end do
      call init_reactions(network, populations, reactions)
      allocate(model, source=network)
    case default
      call require(.false., 'unknown model kind')
    end select
    word = t%next()
    call require(len(word) == 0, 'unexpected trailing tokens')
  end subroutine
end module
