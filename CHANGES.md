Change history for MP-Opt-Model
===============================


Version 3.0 - *Oct 8, 2020*
---------------------------

#### 10/8/20
  - Release 3.0.

#### 9/23/20
  - Move `have_feature()` to [MP-Test][8] and respective feature detection
    functions to [MP-Test][8], [MIPS][9], and [MATPOWER][1].
    - MP-Test
      - `have_feature()`
      - `have_feature_matlab()`
      - `have_feature_octave()`
    - MIPS
      - `have_feature_lu_vec()`
      - `have_feature_pardiso_legacy()`
      - `have_feature_pardiso_object()`
      - `have_feature_pardiso()`
    - MATPOWER
      - `have_feature_e4st()`
      - `have_feature_minopf()`
      - `have_feature_most()`
      - `have_feature_pdipmopf()`
      - `have_feature_regexp_split()`
      - `have_feature_scpdipmopf()`
      - `have_feature_sdp_pf()`
      - `have_feature_smartmarket()`
      - `have_feature_syngrid()`
      - `have_feature_tralmopf()`

#### 9/18/20
  - Add `have_feature()` as a modular, extensible alternative
    to `have_fcn()`, where the detection of a feature named
    `<tag>` is implemented by the function `have_feature_<tag>()`.
  - Make `have_fcn()` a simple wrapper around the new `have_feature()`.

#### 9/16/20
  - Add `set_params()` method to modify parameters of existing
    variables, costs and constraints of an MP-Opt-Model object.
  - Calling `params_var()` method with empty `idx` no longer results
    in fatal error.
  - For `opt_model`, fixed incorrect evaluation of constant term in
    vector valued quadratic costs with constant term supplied as a
    vector.
  - Simplified logic to determine whether a quadratic cost for an
    MP-Opt-Model object is vector vs. scalar valued. If the quadratic
    coefficient is supplied as a matrix, the cost is scalar varied,
    otherwise it is vector valued.

#### 9/14/20
  - Allow `v0`, `vl`, and `vu` inputs to `opt_model/add_var()` method,
    and `l` and `u` inputs to `opt_model/add_lin_constraint()` to
    be scalars that get expanded automatically to the appropriate
    vector dimension.

#### 9/11/20
  - Add `get_soln()` method to `opt_model` for extracting solved
    results for a given named set of variables, constraints or costs.
  - Add `parse_soln()` method which returns a struct with a complete
    set of solution vector and shadow price values for a solved model.

#### 9/10/20
  - Add caching of problem_type() return value.

#### 9/1/20
  - Add support for OSQP solver from [https://osqp.org][7] for LP
    and QP problems, including functions `qps_osqp()`, `osqpver()`,
    and `osqp_options()`.

#### 8/31/20
  - Save the results of `solve()` method to the `soln` field of the
    MP-Opt-Model object.

#### 8/28/20
  - Add `eval_lin_constraint()` method to evaluate the constraint
    values for the full set or an individual named subset of linear
    constraints.

#### 8/27/20
  - Starting point supplied to `solve()` via `opt.x0` is no longer
    ignored for nonlinear equations.


Version 2.1 - *Aug 25, 2020*
----------------------------

#### 8/25/20
  - Release 2.1.
  - Add core NLEQ solver function `nleqs_core()` with arbitrary,
    user-defined update function, used to implement Gauss-Seidel and
    Newton solvers, `nleqs_gauss_seidel()` and `nleqs_newton()`.

#### 8/20/20
  - Add linear equation (`'LEQ'`) problem type for models with equal
    number of variables and linear equality constraints, no costs,
    and no inequality or nonlinear equality constraints. Solved via
    `mplinsolve()`.
  - The `solve()` method of `opt_model` can now automatically handle
    mixed systems of equations, with both linear and nonlinear equality
    constraints.

#### 7/30/20
  - Add fast-decoupled Newton's and Gauss-Seidel solvers for nonlinear
    equations. Use `alg = 'FD'` and `alg = 'GS'` with `nleqs_master()`.
    See also `nleqs_fd_newton()` and `nleqs_gauss_seidel()`.

#### 7/29/20
  - Allow `solve()` method to pass along number of requested output
    arguments `*_master()` solver functions.
  - **INCOMPATIBLE CHANGE:** In `output` return value from
    `nleqs_newton()`, changed the `normF` field of `output.hist` to
    `normf`, for consistency in using lowercase `f` everywhere.


Version 2.0 - *Jul 8, 2020*
---------------------------

#### 7/8/20
  - Release 2.0.

#### 7/3/20
  - Add to `eval_nln_constraint()` method the ability to compute
    constraints for a single named set.

#### 7/2/20
  - Skip evaluation of gradient if `eval_nln_constraint()` is called
    with a single output argument.
  - Add `params_nln_constraint()` method, and add documentation to the
    manual for it and `params_nln_cost()`.

#### 7/1/20
  - Add `mpopt2nleqopt()` to create or modify an `nleqs_master()`
    options struct from a MATPOWER options struct.
  - Add table of valid `have_fcn()` input tags to User's Manual.

#### 6/24/20
  - Add support for nonlinear equations (NLEQ) to `opt_model`. For
    problems with only nonlinear equality constraints and no costs, the
    `problem_type()` method returns `'NLEQ'` and the `solve()` method
    calls `nleqs_master()` to solve the problem.
  - Add tests for solving LP/QP, MILP/MIQP, NLP and NLEQ problems via
    `opt_model/solve()`.

#### 6/17/20
  - Add `nleqs_master()` function as unified interface for solving
    nonlinear equations, including implementations for `fsolve` and
    Newton's method in functions `nleqs_fsolve()` and `nleqs_newton()`,
    respectively.

#### 6/16/20
  - Add new `'fsolve'` tag to `have_fcn()` to check for availability of
    `fsolve()` function.

#### 5/11/20
  - Remove redundant MIPS tests from `test_mp_opt_model`.


Version 1.0 - *May 8, 2020*
---------------------------

#### 5/8/20
  - Release 1.0.

#### 5/7/20
  - Add MP-Opt-Model User's Manual in `docs`, with LaTeX source in
    `docs/src`.
  - Add usage examples to `README.md`.

#### 4/30/20
  - Add `README.md`, `CHANGES.md`, `AUTHORS`, `CONTRIBUTING.md`, `LICENSE`.
  - Refactor `opt_model` class to inherit from new abstract base class
    `mp_idx_manager`, which can be used to manage the indexing of other sets of
    parameters, etc. in other contexts.


Version 0.8 - *Apr 29, 2020*
----------------------------

#### 4/29/20
  - Version 0.8.
  - Add `mpomver()` to define MP-Opt-Model version information.
  - **INCOMPATIBLE CHANGE:** Renamed the following functions and
    modified the order of their input args so that the MP-Opt-Model
    object appears first. Ideally, these would be defined as methods
    of the `opt_model` class, but Octave 4.2 and earlier is not
    able to find them via a function handle (as used in the `solve()`
    method) if they are inherited by a sub-class.
    - `opf_consfcn()` --> `nlp_consfcn()`
    - `opf_costfcn()` --> `nlp_costfcn()`
    - `opf_hessfcn()` --> `nlp_hessfcn()`
  - Add GitHub Actions CI workflow and [Travis-CI][3] configuration.
  - Add `test_mp_opt_model()` to run all tests.
  - Remove MATPOWER dependencies.
  - Move code related to solver interfaces, `opt_model` and a
    few other functions like `have_fcn()` and `nested_struct_copy()`
    from main [MATPOWER][1] repository to new [MP-Opt-Model][2]
    repository.


⬆ _In [MP-Opt-Model][2] repository_ ⬆

-----------------------------------

⬇ _In [MATPOWER][1] repository_ ⬇


#### 4/28/20
  - Move deprecated `opt_model` methods and code related to legacy
    user-defined OPF costs from `@opt_model` to `@opf_model`.
  - **INCOMPATIBLE CHANGE:** Modify order of default output arguments of
    `opt_model/get_idx()` (again), removing the one related to legacy
    costs.

#### 3/18/20
  - Add `nlpopf_solver()` based on the new `solver()` method of
    `opt_model`. This single function replaces `mipsopf_solver()`,
    `fmincopf_solver()`, `ipoptopf_solver()`, and `ktropf_solver()`.
  - Convert `dcopf_solver()` to use the new `solver()` method of
    `opt_model` instead of calling `qps_matpower()` manually.
  - Add new top-level wrapper function `nlps_matpower()` to provide
    a standard interface for MATPOWER's nonlinear program (NLP)
    solvers (`fmincon`, IPOPT, Artelys Knitro, and MIPS), with
    calling syntax similar to `mips()`. It includes the ability to
    pass in solver-specific input options.
  - Add `nlps_fmincon()`, `nlps_ipopt()` and `nlps_knitro()`, with
    interface that matches `nlps_matpower()` to handle implementation
    for `fmincon`, IPOPT, and Artelys Knitro solvers, respectively.
  - Add `mpopt2nlpopt()` to set up an options struct for
    `nlps_matpower()` based on a MATPOWER options struct.
  - Add three new methods to `opt_model` class:
    - `is_mixed_integer()` - returns true if the model includes any binary
      or integer variables
    - `problem_type()` - returns one of the following strings, based on
      the characteristics of the variables, costs and constraints in the
      model:
      - `'NLP'` - nonlinear program
      - `'LP'` - linear program
      - `'QP'` - quadratic program
      - `'MILP'` - mixed-integer linear program
      - `'MIQP'` - mixed-integer quadratic program
    - `solve()` - solves the model using `qps_matpower()`,
      `miqps_matpower()`, or `nlps_matpower()`, depending on the problem
      type (`'MINLP'` problems are not yet implemented)

#### 3/12/20
  - Fix bug in `ktropf_solver()` where Artelys Knitro was still using
    `fmincon` options.

#### 3/6/20
  - Fix issue with missing objective function value from `miqps_mosek()`
    and `qps_mosek()` when return status is "Stalled at or near optimal
    solution."

#### 3/4/20
  - Remove unused input arguments from `opf_consfcn()` and `opf_hessfcn()`.

#### 2/27/20
  - Add `copy()` method to `opt_model` class to get around issues
    with inheritance in constructors that was preventing copy constructor
    from working in Octave 5.2 and earlier (see also [Octave bug
    52614](https://savannah.gnu.org/bugs/?52614).

#### 2/26/20
  - Significant performance improvement for CPLEX on small problems by
    eliminating call to `cplexoptimset()`, which was a huge bottleneck.
  - Fix CPLEX 12.10 compatibility [issue #90][6].

#### 2/18/20
  - Artelys Knitro 12.1 compatibility fix.

#### 8/15/19
  - Improve performance of `opt_model/add_named_set()`.
    (See [issue #79][5].)
    *Thanks to Baraa Mohandes.*
  - Refactor code in `opt_model/params_lin_constraint()` and
    `opt_model/params_quad_cost()` to speed up sparse matrix construction
    when there are lots of constraint or cost sets. Results in significant
    speedups for some problems during problem setup in MOST.
    (See [pull request #70][4].)
    *Thanks to Daniel Muldrew.*


Version 0.7.0 - *Jun 20, 2019*
------------------------------

#### 6/20/19
  - This change history begins with the code that was part of the
    MATPOWER 7.0 release, which is tagged as version 0.7.0 in the
    MP-Opt-Model repository.

----
[1]: https://github.com/MATPOWER/matpower
[2]: https://github.com/MATPOWER/mp-opt-model
[3]: https://travis-ci.org
[4]: https://github.com/MATPOWER/matpower/pull/70
[5]: https://github.com/MATPOWER/matpower/issues/79
[6]: https://github.com/MATPOWER/matpower/issues/90
[7]: https://osqp.org
[8]: https://github.com/MATPOWER/mptest
[9]: https://github.com/MATPOWER/mips
