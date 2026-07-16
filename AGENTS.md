# AGENTS.md

This file provides guidance to AIs when working with code in this repository.

## What this is

OpenVort simulates the motion of quantized vortices in superfluid He-4 (He II) using the vortex filament model of K. W. Schwarz. Pure C (C99-ish, GNU extensions), parallelized with OpenMP. Dependencies: libconfig, OpenMP, glibc (`feenableexcept`). GPLv3.

## Building and running

CMake with a `build/` directory (Ninja is currently configured):

```sh
cmake -B build          # configure (once)
cmake --build build     # produces build/vortices
```

Run a simulation:

```sh
./build/vortices -c <config file> -o <output directory>
```

`example_config.cfg` documents the config format (libconfig syntax; floats must be written with a decimal point — libconfig does not coerce ints to floats). The full list of recognized options is the chain of `config_lookup_*` calls in `src/configuration.c` (`load_conf`/`setup_init`); config keys map to `global_*` variables declared in `include/vortex_constants.h`.

There are no automated tests. Verification is done by running short simulations and inspecting output. Note `feenableexcept` in `main.c`: the code traps FP overflow/underflow/invalid/div-by-zero, so numerical bugs crash with SIGFPE rather than propagating NaNs.

Output is `init.dat` plus `frame%08d.dat` files in the output directory, one every `frame_shots` steps. Post-processing/plotting scripts live in `python/` (e.g. `python3 python/animate.py --config <config> <output dir>`; use `--slow` with periodic boundaries).

## Architecture

The entire simulation state lives in `struct tangle_state` (`include/tangle.h`):

- Vortex filaments are discretized into nodes stored in flat parallel arrays (`vnodes`, `vels`, `tangents`, `normals`, precomputed segment arrays `seg_*`, etc.), all of length `N`.
- Connectivity is via `connections[i].forward/.reverse` index pairs — vortices are doubly-linked rings threaded through the arrays, not contiguous. Empty slots have status `EMPTY` and are recycled through `next_free`/`get_tangle_next_free()`; arrays grow via `expand_tangle()`. Node status (`FREE`, `PINNED`, `PINNED_SLIP`) tracks wall pinning.
- Walking a vortex means following `connections`; `step_node()` steps along a filament handling boundaries/pins.

Main loop (`src/main.c`): per time step — `reconnect()` → `eliminate_small_loops()` → optional frame save → `inject_vortices()` → `update_tangle()` (tangents/normals then velocities) → `rk4_step()` → `remesh()` → `eliminate_small_loops()` → `enforce_boundaries()`. Runs until no points remain or `max_steps`.

Key modules:

- `src/tangle.c` — the core: tangent/normal computation, Biot–Savart velocity integration (`calculate_vs*`, `update_velocities`), local induction + nonlocal contribution, remeshing, loop elimination. The largest and most performance-critical file; OpenMP loops are here.
- `src/vortex_dynamics.c` — RK4 time stepping and reconnection logic (distance + minimum-angle criteria, `reconnection_distance` / `reconnection_angle_cutoff`).
- `src/octree.c` — Barnes–Hut approximation of the Biot–Savart integral (optional, `use_BH` config; supports a 2D quadtree mode). Tree nodes carry centre of mass and a circulation tensor.
- `src/boundary_images.c` — boundary conditions implemented as *image tangles*: precomputed tables of shifted/reflected copies of the domain (`struct boundary_images`, e.g. `periodic_6`, `wall_1_6`, `channel_z`), selected by the `boundaries` config string. Velocity sums include contributions from these images via `shifted()`.
- `src/external_velocity.c` / `include/external_velocity.h` — pluggable space/time-dependent external normal (`vn_conf`) and superfluid (`vs_conf`) velocity fields, plus optional moving boundary (`vb_conf`); each type is a named entry with a parameter table, so new fields are added here.
- `src/vortex_constants.c` — definitions of all `global_*` simulation parameters (physical constants, discretization, mutual friction, BH settings) that `configuration.c` populates.
- `src/vec3_maths.c` — vec3/mat3 helpers used everywhere.

Units are CGS (cm, s) throughout; the circulation quantum `KAPPA` and mutual friction coefficients `alpha`/`alpha_p` are the main physical parameters.

## Caveats

- Simulation restarting is only rudimentarily implemented; frame files lack the frame time, which matters for time-varying external velocities.
- Linux is the only reliably supported platform.
