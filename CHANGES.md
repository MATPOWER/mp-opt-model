Change history for MP-Opt-Model
===============================


Changes since 0.8
-----------------

#### 4/30/20
  - Add `README.md`, `CHANGES.md`, `AUTHORS`, `CONTRIBUTING.md`, `LICENSE`.
  - Refactor `@opt_model class` to use new abstract `@mp_idx_manager`
    class, which can be used to manage the indexing of other sets of
    parameters, etc. in other contexts.


Version 0.8 - *Apr 29, 2020*
----------------------------

#### 4/29/20
  - Version 0.8.
  - Add `mpomver()` to define MP-Opt-Model version information.
  - **INCOMPATIBLE CHANGE:** Renamed the following functions and
    modified the order of their input args so that the MP-Opt-Model
    object appears first. Ideally, these would be defined as methods
    of the `@opt_model` class, but Octave 4.2 and earlier is not
    able to find them via a function handle (as used in the `solve()`
    method) if they are inherited by a sub-class.
    - `opf_consfcn()` --> `nlp_consfcn()`
    - `opf_costfcn()` --> `nlp_costfcn()`
    - `opf_hessfcn()` --> `nlp_hessfcn()`
  - Add GitHub Actions CI workflow and [Travis-CI][3] configuration.
  - Add `test_mp_opt_model()` to run all tests.
  - Remove MATPOWER dependencies.
  - Move code related to solver interfaces, `@opt_model` and a
    few other functions like `have_fcn()` and `nested_struct_copy()`
    from main [MATPOWER][1] repository to new [MP-Opt-Model][2]
    repository.

----
[1]: https://github.com/MATPOWER/matpower
[2]: https://github.com/MATPOWER/mp-opt-model
[3]: https://travis-ci.org
