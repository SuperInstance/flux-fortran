!> FLUX Fracture Engine — Fortran 2008
!> Disjointed linear algebra that fractures constraint systems into
!> independent blocks and coalesces results via bitwise OR.
!>
!> DESIGN: Fortran's column-major arrays and array syntax FORCE
!> a particular organization. The adjacency matrix is natural.
!> Connected components via BFS. Bitwise OR coalescence.
!> Zero false negatives guaranteed by Boolean algebra.

module flux_fracture
    implicit none
    private
    public :: fracture_result, build_dependency_graph, &
              find_blocks, coalesce_masks, adaptive_fracture

    integer, parameter :: MAX_CONSTRAINTS = 256
    integer, parameter :: MAX_DIMS = 256
    integer, parameter :: MAX_BLOCKS = MAX_CONSTRAINTS

    !> A single independent block of constraints
    type :: block_type
        integer, allocatable :: constraint_indices(:)
        integer, allocatable :: dim_indices(:)
        integer :: n_constraints = 0
        integer :: n_dims = 0
    end type

    !> Result of fracturing a constraint system
    type :: fracture_result
        type(block_type), allocatable :: blocks(:)
        integer :: n_blocks = 0
        integer :: largest_block = 0
        real :: speedup_potential = 1.0
    end type

contains

    !> Build bipartite dependency graph: constraint i depends on dimension j
    !> Returns adjacency matrix adj(i,j) = 1 if constraint i checks dimension j
    subroutine build_dependency_graph(constraints, n_constraints, n_dims, adj)
        real(8), intent(in) :: constraints(:,:)  ! (n_constraints, 2) = lo, hi per dim
        integer, intent(in) :: n_constraints, n_dims
        integer, intent(out) :: adj(:,:)  ! (n_constraints, n_dims)

        integer :: i, j
        real(8) :: lo, hi

        adj = 0
        do i = 1, n_constraints
            do j = 1, n_dims
                lo = constraints(i, 1)
                hi = constraints(i, 2)
                ! A constraint depends on dimension j if its bounds are finite
                ! For simplicity: each constraint i checks dimension i (diagonal)
                ! In general, this would be determined by the constraint structure
                if (i == j .or. lo /= 0.0d0 .or. hi /= 0.0d0) then
                    adj(i, j) = 1
                end if
            end do
        end do
    end subroutine

    !> Find connected components of the bipartite dependency graph
    !> Uses BFS from each unvisited node
    function find_blocks(adj, n_constraints, n_dims) result(fr)
        integer, intent(in) :: adj(:,:)
        integer, intent(in) :: n_constraints, n_dims
        type(fracture_result) :: fr

        integer :: visited_constraints(MAX_CONSTRAINTS)
        integer :: visited_dims(MAX_DIMS)
        integer :: block_cons(MAX_CONSTRAINTS)
        integer :: block_dims(MAX_DIMS)
        integer :: queue(MAX_CONSTRAINTS + MAX_DIMS)
        integer :: q_head, q_tail, current, b_idx, n_bc, n_bd
        integer :: i, j, k, bsize

        visited_constraints = 0
        visited_dims = 0

        allocate(fr%blocks(MAX_BLOCKS))
        fr%n_blocks = 0
        fr%largest_block = 0

        do i = 1, n_constraints
            if (visited_constraints(i) == 1) cycle

            ! BFS from constraint i
            block_cons = 0
            block_dims = 0
            n_bc = 0
            n_bd = 0
            queue = 0
            q_head = 1
            q_tail = 1

            ! Start BFS from constraint i
            queue(q_tail) = i  ! positive = constraint
            q_tail = q_tail + 1
            visited_constraints(i) = 1

            do while (q_head < q_tail)
                current = queue(q_head)
                q_head = q_head + 1

                if (current > 0) then
                    ! It's a constraint node
                    n_bc = n_bc + 1
                    block_cons(n_bc) = current

                    ! Find connected dimension nodes
                    do j = 1, n_dims
                        if (adj(current, j) == 1 .and. visited_dims(j) == 0) then
                            visited_dims(j) = 1
                            queue(q_tail) = -j  ! negative = dimension
                            q_tail = q_tail + 1
                        end if
                    end do
                else
                    ! It's a dimension node
                    j = -current
                    n_bd = n_bd + 1
                    block_dims(n_bd) = j

                    ! Find connected constraint nodes
                    do k = 1, n_constraints
                        if (adj(k, j) == 1 .and. visited_constraints(k) == 0) then
                            visited_constraints(k) = 1
                            queue(q_tail) = k
                            q_tail = q_tail + 1
                        end if
                    end do
                end if
            end do

            ! Store block
            fr%n_blocks = fr%n_blocks + 1
            b_idx = fr%n_blocks

            allocate(fr%blocks(b_idx)%constraint_indices(n_bc))
            allocate(fr%blocks(b_idx)%dim_indices(n_bd))
            fr%blocks(b_idx)%constraint_indices(1:n_bc) = block_cons(1:n_bc)
            fr%blocks(b_idx)%dim_indices(1:n_bd) = block_dims(1:n_bd)
            fr%blocks(b_idx)%n_constraints = n_bc
            fr%blocks(b_idx)%n_dims = n_bd

            bsize = max(n_bc, fr%largest_block)
            fr%largest_block = bsize
        end do

        ! Trim allocation
        fr%speedup_potential = real(n_constraints) / real(max(fr%largest_block, 1))

        ! Trim blocks array
        call trim_blocks(fr)
    end function

    !> Coalesce block error masks via bitwise OR
    function coalesce_masks(block_masks, n_blocks, n_bits) result(global_mask)
        integer, intent(in) :: block_masks(:,:)
        integer, intent(in) :: n_blocks, n_bits
        integer :: global_mask(n_bits)

        integer :: i, b

        global_mask = 0
        do b = 1, n_blocks
            do i = 1, n_bits
                global_mask(i) = ior(global_mask(i), block_masks(b, i))
            end do
        end do
    end function

    !> Verify coalescence preserves all violations (zero false negatives)
    function verify_coalescence(coalesced, monolithic, n_bits) result(ok)
        integer, intent(in) :: coalesced(:), monolithic(:), n_bits
        logical :: ok
        integer :: i

        ok = .true.
        do i = 1, n_bits
            ! Coalesced must have AT LEAST the violations of monolithic
            if (monolithic(i) == 1 .and. coalesced(i) == 0) then
                ok = .false.
                return
            end if
        end do
    end function

    !> Adaptive fracture: re-fracture when structure changes
    function adaptive_fracture(prev_fr, new_adj, n_c, n_d) result(new_fr)
        type(fracture_result), intent(in) :: prev_fr
        integer, intent(in) :: new_adj(:,:), n_c, n_d
        type(fracture_result) :: new_fr

        new_fr = find_blocks(new_adj, n_c, n_d)
        ! Compare with previous — if same n_blocks and same sizes, no change needed
    end function

    ! --- Internal helpers ---

    subroutine trim_blocks(fr)
        type(fracture_result), intent(inout) :: fr
        type(block_type), allocatable :: temp(:)

        if (fr%n_blocks > 0) then
            allocate(temp(fr%n_blocks))
            temp(1:fr%n_blocks) = fr%blocks(1:fr%n_blocks)
            deallocate(fr%blocks)
            allocate(fr%blocks(fr%n_blocks))
            fr%blocks = temp
        end if
    end subroutine

end module flux_fracture
