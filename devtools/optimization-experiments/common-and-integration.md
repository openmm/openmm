# Common / integration port report

The implementation is ported into `openmm-pr` against current master commit `0c4bcaba734ea574f52f5e09ffaaf4d6fe6109c2`. It does not modify the running Folding installation. No GPU execution, performance measurement, or numerical simulation validation was performed by this porting task.

All newly contributed performance experiments require explicit opt-in. Boolean flags require exactly `1`; the SETTLE block-size option selects the alternative only for `128`. Trace flags remain independently opt-in. This differs intentionally from the local deployed 8.2 build's default-on paths.

## Feature accounting

| Feature / setting | State | Port details / applicability |
| --- | --- | --- |
| A1 potential-only Monte Carlo barostat energy | Ported, new flag `OPENMM_EXPERIMENT_BAROSTAT_POTENTIAL_ONLY` | Exact `LangevinMiddleIntegrator` type only; fallback performs the original `State::Energy` query. Preserves latest master rigid-molecule/per-atom scaling count. No bypass of force evaluation or acceptance calculation. |
| `OPENMM_EXPERIMENT_BONDED_FORCE_ONLY` | Ported | Mixed precision, one context, recognized exact ordinary Force classes and no parameter derivatives. Keeps the original energy-bearing path for other requests. |
| `OPENMM_EXPERIMENT_BONDED_STATIC_FORCE_ONLY` | Ported | Only specializes the optional force-only clone. ABI argument slots remain but request values become device compile-time constants. |
| `OPENMM_EXPERIMENT_BONDED_FORCE_FENCE_ELISION` | Ported, unvalidated device experiment | Only wrapper-owned fence offsets in the optional force-only clone; empty prefix and no derivatives required. Device guard retains fences unless CUDA mixed precision and non-HIP. No claim that compiler output changes only fences. |
| `OPENMM_EXPERIMENT_WATER_BOND_INTEGER_SUM` | Ported | Topology-validated disjoint harmonic/constraint triangles; preserves original per-bond conversion before unsigned fixed-point summation. Residual bonds retain original indices; full original loop for energy-bearing requests. |
| `OPENMM_EXPERIMENT_WATER_BOND_STATIC_DIRECTION` | Ported | Only when every immutable topology-map direction matches; parameters remain dynamically indexed. |
| `OPENMM_EXPERIMENT_GETPOSITIONS_PINNED_TAIL` | Ported / adapted | Preserves current master's two `allowPeriodic` output branches. Mixed CUDA, one context, N >= 4096, CPU PME disabled, array sizes/layout checked. Uses byte copies from raw tail storage. |
| `OPENMM_EXPERIMENT_GETPOSITIONS_SINGLE_ALLOCATION` | Ported | Lower-priority alternative when pinned-tail path is not selected. Avoids redundant intermediate N-size allocation before device-size download. |
| `OPENMM_EXPERIMENT_PARALLEL_SET_POSITIONS` | Ported / adapted | Works with master's packed xyz buffer; original upload and GPU xyz-copy path retained. MSVC x64 worker MXCSR controls must match; full serial rewrite before upload on mismatch. |
| `OPENMM_EXPERIMENT_FUSED_SET_POSITIONS` | Ported / adapted | Packs xyz into first 12N bytes of pinned memory, correction into aligned 16P-byte tail; fits existing 32P-byte mixed-velocity allocation. Uses master's floatBuffer upload and xyz-only GPU copy. No old posq download restored. Unmasked SSE traps reject fusion. |
| `OPENMM_EXPERIMENT_BOX_POSITION_SCRATCH_REUSE` | Ported | Context-owned temporary Vec3 vector; getter success is required in this call before use. Optional allocation failure falls back to local vector. |
| `OPENMM_EXPERIMENT_REUSE_BOX_POSQ_DOWNLOAD` | Obsolete / already addressed structurally upstream | Latest master setter has no posq download at all. Reintroducing the old layout just to reuse that download would regress master. The setting is intentionally absent. |
| `OPENMM_LANGEVIN_MIDDLE_SETTLE_FUSION` | Ported | Exact Middle integrator, mixed CUDA, one context, no virtual sites, validated disjoint positive-finite-mass SETTLE partition. Random generation/count/index and Part2 launch mapping remain unchanged. |
| `OPENMM_EXPERIMENT_MIDDLE_FUSED_TAIL` | Ported | Requires SETTLE fusion; complement list covers real residual atoms, with SHAKE/CCMA before fused finalization. |
| `OPENMM_EXPERIMENT_MIDDLE_KICK_SETTLE` | Ported | Requires SETTLE fusion plus fused tail; combines Part1 with velocity SETTLE, retaining residual velocity constraints before stochastic Part2. |
| `OPENMM_EXPERIMENT_SETTLE_BLOCK_SIZE` | Ported | Alternative value `128`, bounded by compiled kernel maximum; no global launch-geometry changes. |
| `OPENMM_EXPERIMENT_SETTLE_RELOAD_DELTA` | Ported | Optional volatile reload of original deltas in fused SETTLE. Does not assert register or instruction-count improvements. |
| `OPENMM_EXPERIMENT_VELOCITY_SETTLE_BLOCK128` | Ported | Mixed CUDA, one context, kernel maximum at least 128; ordinary public velocity constraint calls still solve all constraints. |
| `OPENMM_EXPERIMENT_CMM_COMPACT_BUFFER` | Ported | CUDA mixed / one context and at least two blocks; buffer sized to the actual 64-thread producer grid, with reduction result slot when necessary. |
| `OPENMM_EXPERIMENT_CMM_REDUCE_ONCE` | Ported | One 64-thread second-stage reduction followed by independent velocity application; unchanged per-invocation timing and arithmetic tree. |
| `OPENMM_EXPERIMENT_CMM_WARP_REDUCTION` | Ported | Requires reduce-once and CUDA Volta+; one CTA barrier plus warp shuffles retaining original binary32 tree. |
| `OPENMM_EXPERIMENT_PARALLEL_BAROSTAT` / separate legacy kernel | Already upstream | Current master contains large-molecule workgroup scaling from PR #5432 plus newer shear/component support. Its functions and kernel source are retained, without reintroducing the historical duplicate. |

Five trace settings are retained for SETTLE fusion, tail, kick, occupancy and velocity block size. They do not enable the corresponding optimization by themselves.
