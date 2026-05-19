# flux-fortran

Disjointed linear algebra for constraint systems — Fortran 2008.

## How It Works

Constraint systems often have hidden independence: if constraint A checks dimension 1 and constraint B checks dimension 2, they can be evaluated in parallel. This library detects that independence using BFS on a bipartite constraint×dimension dependency graph, then splits the system into independent blocks.

The error mask is 1 bit per constraint. Eight constraints = one byte. When blocks are independent, you can coalesce their results with bitwise OR — zero false negatives, mathematically guaranteed.

Sediment layers accumulate edge-case corrections over time. Each layer tightens a bound. Correctness is monotonically increasing — adding layers never makes the system less correct.

## Build & Test

```bash
make test
```

Requires `gfortran`. Outputs:

```
Results: 15/15 passed
ALL TESTS PASSED
```

## Modules

### `flux_fracture` — Split independent constraints

```fortran
use flux_fracture
type(DependencyGraph) :: graph
type(FractureResult) :: result

! Build adjacency: 8 constraints, 8 dimensions
call graph%from_identity(8)
result = graph%find_blocks()
! result%n_blocks == 8  →  8 independent blocks, 8× parallelism
```

Key types:
- `DependencyGraph` — bipartite adjacency matrix with flat array storage
- `FractureResult` — blocks, count, speedup potential
- `coalesce_masks(block_masks, n_blocks, n_dims)` — bitwise OR merge

### `flux_sediment` — Accumulated correctness

```fortran
use flux_sediment
type(SedimentStack) :: stack

call add_layer(stack, 3, -10.0d0, 10.0d0, 2.0d0, 1)
call apply_sediment(stack, values, original_mask, 8, corrected_mask)
```

Key types:
- `SedimentLayer` — one correction: constraint index, corrected bounds, surprise, timestamp
- `SedimentStack` — fixed array of 50 layers, supersede oldest when full

## What Fortran Teaches Us

Fortran stores arrays column-major. The adjacency matrix `adj(constraint, dimension)` means iterating over constraints in the inner loop touches contiguous memory. The language chose the layout; you organize around it.

Fortran also forces fixed-size arrays and explicit bounds. There's no `Vec::push()` — you declare `real(8) :: values(8)` and that's what you get. The constraint engine's hot path never allocates. Fortran makes that the default, not something you have to enforce.

Column-major storage, preallocated arrays, subroutines that mutate arguments — this is the shape of every high-performance numeric system. Fortran just refuses to let you pretend otherwise.

## Files

| File | What |
|------|------|
| `src/flux_fracture.f90` | Dependency graph, BFS connected components, coalesce |
| `src/flux_sediment.f90` | Sediment stack, add/apply, correctness density |
| `src/flux_test.f90` | 15 tests: fracture, coalesce, sediment, NaN, full pipeline |

## Where to Go Next

| If you want to... | Go to... |
|---|---|
| See the same math in Rust | [flux-fracture](https://github.com/SuperInstance/flux-fracture) |
| See the same math in C | [flux-fracture-c](https://github.com/SuperInstance/flux-fracture-c) |
| Read what old languages teach | [OLD-LANGUAGE-ARCHITECTURE.md](https://github.com/SuperInstance/constraint-theory-ecosystem/blob/main/docs/OLD-LANGUAGE-ARCHITECTURE.md) |
| See the full ecosystem | [constraint-theory-ecosystem](https://github.com/SuperInstance/constraint-theory-ecosystem) |

## License

MIT
