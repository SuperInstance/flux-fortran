!> FLUX Constraint Engine Test — Fortran 2008
!> Tests: exact checking, fracture-coalesce, sediment layers
!> Demonstrates the full pipeline: check → fracture → coalesce → sediment

program flux_test
    use flux_fracture
    use flux_sediment
    implicit none

    integer :: n_tests, n_pass

    n_tests = 0
    n_pass = 0

    call test_fracture_independent(n_tests, n_pass)
    call test_fracture_block_diagonal(n_tests, n_pass)
    call test_coalesce_correctness(n_tests, n_pass)
    call test_sediment_corrects(n_tests, n_pass)
    call test_sediment_monotonic(n_tests, n_pass)
    call test_nan_always_violates(n_tests, n_pass)
    call test_full_pipeline(n_tests, n_pass)

    write(*,'(A)') '=============================='
    write(*,'(A,I0,A,I0,A)') 'Results: ', n_pass, '/', n_tests, ' passed'
    if (n_pass == n_tests) then
        write(*,'(A)') 'ALL TESTS PASSED'
    else
        write(*,'(A)') 'SOME TESTS FAILED'
        stop 1
    end if

contains

    subroutine check(test_name, condition, n_tests, n_pass)
        character(*), intent(in) :: test_name
        logical, intent(in) :: condition
        integer, intent(inout) :: n_tests, n_pass

        n_tests = n_tests + 1
        if (condition) then
            n_pass = n_pass + 1
            write(*,'(A,A)') '  PASS: ', test_name
        else
            write(*,'(A,A)') '  FAIL: ', test_name
        end if
    end subroutine

    subroutine test_fracture_independent(n_tests, n_pass)
        integer, intent(inout) :: n_tests, n_pass
        integer :: adj(8, 8), i
        type(fracture_result) :: fr

        write(*,'(A)') '--- Fracture: Independent ---'
        ! Identity adjacency: each constraint checks only its own dimension
        adj = 0
        do i = 1, 8
            adj(i, i) = 1
        end do

        fr = find_blocks(adj, 8, 8)
        call check('8 independent constraints = 8 blocks', fr%n_blocks == 8, n_tests, n_pass)
        call check('speedup = 8.0', abs(fr%speedup_potential - 8.0) < 0.01, n_tests, n_pass)
    end subroutine

    subroutine test_fracture_block_diagonal(n_tests, n_pass)
        integer, intent(inout) :: n_tests, n_pass
        integer :: adj(8, 8)
        type(fracture_result) :: fr

        write(*,'(A)') '--- Fracture: Block Diagonal ---'
        ! Two blocks of 4: {1,2,3,4} and {5,6,7,8}
        adj = 0
        ! Block 1: constraints 1-4 share dimensions 1-4
        adj(1:4, 1:4) = 1
        ! Block 2: constraints 5-8 share dimensions 5-8
        adj(5:8, 5:8) = 1

        fr = find_blocks(adj, 8, 8)
        call check('2 blocks', fr%n_blocks == 2, n_tests, n_pass)
        call check('speedup = 2.0', abs(fr%speedup_potential - 2.0) < 0.01, n_tests, n_pass)
    end subroutine

    subroutine test_coalesce_correctness(n_tests, n_pass)
        integer, intent(inout) :: n_tests, n_pass
        integer :: block_masks(2, 8)
        integer :: global_mask(8)
        integer :: monolithic(8)
        logical :: ok
        integer :: ii

        write(*,'(A)') '--- Coalesce Correctness ---'
        ! Block 1 catches violations in dims 1-4
        block_masks(1, :) = [1, 0, 1, 0, 0, 0, 0, 0]
        ! Block 2 catches violations in dims 5-8
        block_masks(2, :) = [0, 0, 0, 0, 0, 1, 0, 1]

        global_mask = coalesce_masks(block_masks, 2, 8)
        monolithic = [1, 0, 1, 0, 0, 1, 0, 1]

        ok = .true.
        do ii = 1, 8
            if (monolithic(ii) == 1 .and. global_mask(ii) == 0) ok = .false.
        end do
        call check('coalescence preserves all violations', ok, n_tests, n_pass)

        ! Also check the OR is correct
        ok = all(global_mask == monolithic)
        call check('OR result matches monolithic', ok, n_tests, n_pass)
    end subroutine

    subroutine test_sediment_corrects(n_tests, n_pass)
        integer, intent(inout) :: n_tests, n_pass
        type(sediment_stack) :: stack
        real(8) :: values(8)
        integer :: original(8), corrected(8)

        write(*,'(A)') '--- Sediment Corrects Edge Cases ---'
        ! Value at -11 should be OUTSIDE bounds [-10, 10]
        values = [1.0d0, 2.0d0, -11.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0]
        ! Original mask misses the violation on dim 3
        original = [0, 0, 0, 0, 0, 0, 0, 0]

        ! Add sediment layer correcting dim 3 bounds
        call add_layer(stack, 3, -10.0d0, 10.0d0, 2.0d0, 1)

        call apply_sediment(stack, values, original, 8, corrected)

        call check('sediment catches dim 3 violation', corrected(3) == 1, n_tests, n_pass)
        call check('other dims unchanged', all(corrected([1,2,4,5,6,7,8]) == 0), n_tests, n_pass)
    end subroutine

    subroutine test_sediment_monotonic(n_tests, n_pass)
        integer, intent(inout) :: n_tests, n_pass
        type(sediment_stack) :: stack
        real(8) :: density

        write(*,'(A)') '--- Sediment Monotonic ---'
        ! Add layers one at a time, check density increases
        density = correctness_density(stack, 8)
        call check('empty stack: density = 0', abs(density) < 1e-10, n_tests, n_pass)

        call add_layer(stack, 1, -1.0d0, 1.0d0, 1.0d0, 1)
        density = correctness_density(stack, 8)
        call check('1 layer: density = 0.125', abs(density - 0.125d0) < 1e-10, n_tests, n_pass)

        call add_layer(stack, 2, -2.0d0, 2.0d0, 1.5d0, 2)
        density = correctness_density(stack, 8)
        call check('2 layers: density = 0.25', abs(density - 0.25d0) < 1e-10, n_tests, n_pass)

        call add_layer(stack, 5, -5.0d0, 5.0d0, 2.0d0, 3)
        density = correctness_density(stack, 8)
        call check('3 layers: density = 0.375', abs(density - 0.375d0) < 1e-10, n_tests, n_pass)
    end subroutine

    subroutine test_nan_always_violates(n_tests, n_pass)
        integer, intent(inout) :: n_tests, n_pass
        type(sediment_stack) :: stack
        real(8) :: values(8)
        integer :: original(8), corrected(8)
        real(8) :: nan_val

        write(*,'(A)') '--- NaN Always Violates ---'
        ! IEEE NaN via 0.0/0.0
        nan_val = 0.0d0
        nan_val = nan_val / 0.0d0

        values = [1.0d0, 2.0d0, nan_val, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0]
        original = [0, 0, 0, 0, 0, 0, 0, 0]

        ! Sediment layer on dim 3 with any bounds should still catch NaN
        call add_layer(stack, 3, -1.0d6, 1.0d6, 5.0d0, 1)
        call apply_sediment(stack, values, original, 8, corrected)

        call check('NaN violates even with wide bounds', corrected(3) == 1, n_tests, n_pass)
    end subroutine

    subroutine test_full_pipeline(n_tests, n_pass)
        integer, intent(inout) :: n_tests, n_pass
        integer :: adj(8, 8)
        type(fracture_result) :: fr
        type(sediment_stack) :: stack
        real(8) :: values(8)
        integer :: i

        write(*,'(A)') '--- Full Pipeline ---'
        ! Setup: identity adjacency (8 independent blocks)
        adj = 0
        do i = 1, 8
            adj(i, i) = 1
        end do

        ! Fracture
        fr = find_blocks(adj, 8, 8)
        call check('pipeline fracture: 8 blocks', fr%n_blocks == 8, n_tests, n_pass)

        ! Add sediment
        call add_layer(stack, 3, -10.0d0, 10.0d0, 2.0d0, 1)

        ! Test with edge case
        values = [1.0d0, 2.0d0, -11.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0]

        ! Density should be 1/8
        call check('pipeline density = 0.125', &
            abs(correctness_density(stack, 8) - 0.125d0) < 1e-10, n_tests, n_pass)
    end subroutine

end program flux_test
