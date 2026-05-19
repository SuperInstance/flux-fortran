# flux-fortran

Disjointed linear algebra for constraint systems — Fortran 2008.

Fractures constraint systems into independent blocks via BFS connected-component detection on the dependency graph. Coalesces results via bitwise OR. Zero false negatives guaranteed by Boolean algebra.

## Modules

- `flux_fracture` — DependencyGraph, FractureResult, find_blocks(), coalesce_masks()
- `flux_sediment` — SedimentLayer, SedimentStack, add_layer(), apply_sediment()

## Build & Test

```bash
make test
```

## What Fortran's Shape Teaches Us

Fortran forces array-first thinking. The adjacency matrix is natural here — column-major storage makes the bipartite graph traversal cache-friendly. The sediment stack as a fixed-size OCCURS array is pure COBOL-era thinking: preallocate, never grow, supersede the oldest.

The shape is: **preallocated arrays, explicit bounds, no dynamic allocation in the hot path.** This is how production systems that can't fail are built.
