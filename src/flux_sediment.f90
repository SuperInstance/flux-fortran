!> FLUX Sediment Engine — Fortran 2008
!> Accumulated correctness through sediment layers.
!> Each layer corrects edge cases. Correctness monotonically increases.
!> COBOL theorem: the frozen hot path + accumulated corrections = production system.

module flux_sediment
    implicit none
    private
    public :: sediment_layer, sediment_stack, add_layer, apply_sediment, correctness_density

    integer, parameter :: MAX_LAYERS = 50
    integer, parameter :: MAX_CONSTRAINTS = 256

    !> A sediment layer: edge-case corrections
    type :: sediment_layer
        integer :: constraint_idx = 0       ! which constraint this corrects
        real(8) :: corrected_lo = 0.0d0     ! corrected lower bound
        real(8) :: corrected_hi = 0.0d0     ! corrected upper bound
        real(8) :: surprise = 0.0d0          ! information-theoretic surprise
        integer :: timestamp = 0             ! when added
        logical :: active = .true.
    end type

    !> Stack of sediment layers
    type :: sediment_stack
        type(sediment_layer) :: layers(MAX_LAYERS)
        integer :: depth = 0
        integer :: active_count = 0
        real(8) :: total_surprise = 0.0d0
    end type

contains

    !> Add a correction layer to the stack
    subroutine add_layer(stack, constraint_idx, corrected_lo, corrected_hi, surprise, timestamp)
        type(sediment_stack), intent(inout) :: stack
        integer, intent(in) :: constraint_idx
        real(8), intent(in) :: corrected_lo, corrected_hi, surprise
        integer, intent(in) :: timestamp

        if (stack%depth >= MAX_LAYERS) then
            ! Stack full — supersede oldest layer for same constraint
            call supersede_oldest(stack, constraint_idx, corrected_lo, corrected_hi, surprise, timestamp)
            return
        end if

        stack%depth = stack%depth + 1
        stack%layers(stack%depth)%constraint_idx = constraint_idx
        stack%layers(stack%depth)%corrected_lo = corrected_lo
        stack%layers(stack%depth)%corrected_hi = corrected_hi
        stack%layers(stack%depth)%surprise = surprise
        stack%layers(stack%depth)%timestamp = timestamp
        stack%layers(stack%depth)%active = .true.
        stack%active_count = stack%active_count + 1
        stack%total_surprise = stack%total_surprise + surprise
    end subroutine

    !> Apply sediment corrections to error mask
    !> For each constraint that has a correction: if the original mask says PASS
    !> but the correction says it should FAIL, flip the bit.
    subroutine apply_sediment(stack, values, original_mask, n_constraints, corrected_mask)
        type(sediment_stack), intent(in) :: stack
        real(8), intent(in) :: values(:)
        integer, intent(in) :: original_mask(:), n_constraints
        integer, intent(out) :: corrected_mask(:)

        integer :: i, c_idx
        real(8) :: lo, hi, val

        corrected_mask = original_mask

        do i = 1, stack%depth
            if (.not. stack%layers(i)%active) cycle

            c_idx = stack%layers(i)%constraint_idx
            if (c_idx < 1 .or. c_idx > n_constraints) cycle

            lo = stack%layers(i)%corrected_lo
            hi = stack%layers(i)%corrected_hi

            if (c_idx <= size(values)) then
                val = values(c_idx)
                ! Check with corrected bounds
                if (val < lo .or. val > hi .or. val /= val) then
                    ! NaN check: val /= val is True only for NaN
                    corrected_mask(c_idx) = 1
                end if
            end if
        end do
    end subroutine

    !> Compute correctness density: fraction of constraints covered by sediment
    function correctness_density(stack, n_constraints) result(density)
        type(sediment_stack), intent(in) :: stack
        integer, intent(in) :: n_constraints
        real(8) :: density

        logical :: covered(MAX_CONSTRAINTS)
        integer :: i, n_covered

        covered = .false.
        do i = 1, stack%depth
            if (stack%layers(i)%active) then
                if (stack%layers(i)%constraint_idx >= 1 .and. &
                    stack%layers(i)%constraint_idx <= n_constraints) then
                    covered(stack%layers(i)%constraint_idx) = .true.
                end if
            end if
        end do

        n_covered = count(covered(1:n_constraints))
        density = real(n_covered) / real(n_constraints)
    end function

    ! --- Internal ---

    subroutine supersede_oldest(stack, constraint_idx, lo, hi, surprise, timestamp)
        type(sediment_stack), intent(inout) :: stack
        integer, intent(in) :: constraint_idx, timestamp
        real(8), intent(in) :: lo, hi, surprise

        integer :: i, oldest_idx, oldest_time

        oldest_idx = 1
        oldest_time = stack%layers(1)%timestamp

        do i = 2, stack%depth
            if (stack%layers(i)%timestamp < oldest_time) then
                oldest_time = stack%layers(i)%timestamp
                oldest_idx = i
            end if
        end do

        ! Supersede the oldest
        stack%total_surprise = stack%total_surprise - stack%layers(oldest_idx)%surprise
        stack%layers(oldest_idx)%constraint_idx = constraint_idx
        stack%layers(oldest_idx)%corrected_lo = lo
        stack%layers(oldest_idx)%corrected_hi = hi
        stack%layers(oldest_idx)%surprise = surprise
        stack%layers(oldest_idx)%timestamp = timestamp
        stack%layers(oldest_idx)%active = .true.
        stack%total_surprise = stack%total_surprise + surprise
    end subroutine

end module flux_sediment
