MP-Opt-Model
============

[MP-Opt-Model][1] is a package of MATLAB/Octave M-files for constructing
and solving mathematical optimization problems. It provides a unified
interface for calling numerous LP, QP, mixed-integer and nonlinear
solvers, as well as an object-oriented interface for building and
solving your optimization model.

It is based on code that was originally developed by Ray D. Zimmerman as
part of [MATPOWER][2].


System Requirements
-------------------

*   [MATLAB][3] version 8.6 (R2015b) or later, or
*   [GNU Octave][4] version 4.2 or later
*   [MP-Test][5], for running the MP-Opt-Model test suite
*   [MATPOWER Interior Point Solver (MIPS)][6]


Installation
------------

**Note to [MATPOWER][2] users:** _MP-Opt-Model and its prerequisites, MIPS
and MP-Test, are included when you install [MATPOWER][2]. There is generally
no need to install it separately. You can skip directly to step 3 to verify._

Installation and use of MP-Opt-Model requires familiarity with the basic operation
of MATLAB or Octave, including setting up your MATLAB path.

1.  Clone the repository or download and extract the zip file of the MP-Opt-Model
    distribution from the [MP-Opt-Model project page][1] to the location of your
    choice. The files in the resulting `mp-opt-model` or `mp-opt-modelXXX` directory,
    where `XXX` depends on the version of MP-Opt-Model, should not need to be
    modified, so it is recommended that they be kept separate from your
    own code. We will use `<MPOM>` to denote the path to this directory.

2.  Add the following directories to your MATLAB or Octave path:
    *   `<MPOM>/lib`
    *   `<MPOM>/lib/t`

3.  At the MATLAB/Octave prompt, type `test_mp_opt_model` to run the test suite and
    verify that MP-Opt-Model is properly installed and functioning. (Note: The
    tests require functioning installations of both [MP-Test][5] and
    [MIPS][6]) The result should resemble the following:
```matlab
  >> test_mp_opt_model
  t_nested_struct_copy....ok
  t_have_fcn..............ok
  t_mips..................ok
  t_mips_pardiso..........ok
  t_qps_mips..............ok
  t_qps_matpower..........ok (100 of 396 skipped)
  t_miqps_matpower........ok (102 of 288 skipped)
  t_nlps_matpower.........ok
  t_opt_model.............ok
  All tests successful (1759 passed, 202 skipped of 1961)
  Elapsed time 1.98 seconds.
```

Usage
-----

Until we get some time to write some documentation, there are examples in the
test files in `<MPOM>/lib/t`, as well as in the [`opf_setup()`][12] and
[`opf_execute()`][13] functions in [MATPOWER][2].


Documentation
-------------

There are two primary sources of documentation for MP-Opt-Model.

[ **Note: The MP-Opt-Model User's Manual is not yet available.** The
first is the [MP-Opt-Model User's Manual][7]. It can be found in your
MP-Opt-Model distribution at
`<MPOM>/docs/MP-Opt-Model-manual.pdf` and the latest version is
always available at:
<https://github.com/MATPOWER/mp-opt-model/blob/master/docs/MP-Opt-Model-manual.pdf>. ]

And second is the built-in `help` command. As with the built-in
functions and toolbox routines in MATLAB and Octave, you can type `help`
followed by the name of a command or M-file to get help on that particular
function. Many of the M-files in MP-Opt-Model have such documentation and this
should be considered the main reference for the calling options for each
function, e.g.: `qps_matpower`, `miqps_matpower`, and `nlps_matpower`.


[Citing MP-Opt-Model][10]
-------------------------

**Please ignore the following and simply cite the MP-Opt-Model [GitHub page][1] until the MP-Opt-Model User's Manual
becomes available.**

[ *The following will be updated when the MP-Opt-Model User's Manual becomes available.*

We request that publications derived from the use of the MP-Opt-Model
explicitly acknowledge that fact by citing the [MP-Opt-Model User's Manual][7].
The citation and DOI can be version-specific or general, as appropriate.
For version 0.9, use:

>   R. D. Zimmerman. *MP-Opt-Model User's Manual, Version 0.9*. 2020.
    [Online]. Available: https://matpower.org/docs/MP-Opt-Model-manual-0.9.pdf  
    doi: [10.5281/zenodo.???????](https://doi.org/10.5281/zenodo.???????)

For a version non-specific citation, use the following citation and DOI,
with *\<YEAR\>* replaced by the year of the most recent release:

>   R. D. Zimmerman. *MP-Opt-Model User's Manual*. *\<YEAR\>*.
    [Online]. Available: https://matpower.org/docs/MP-Opt-Model-manual.pdf  
    doi: [10.5281/zenodo.???????][11]

A list of versions of the User's Manual with release dates and
version-specific DOI's can be found via the general DOI at
https://doi.org/10.5281/zenodo.???????.
]

Contributing
------------

Please see our [contributing guidelines][8] for details on how to
contribute to the project or report issues.

License
-------

MP-Opt-Model is distributed under the [3-clause BSD license][9].

----
[1]: https://github.com/MATPOWER/mp-opt-model
[2]: https://matpower.org/
[3]: https://www.mathworks.com/
[4]: https://www.gnu.org/software/octave/
[5]: https://github.com/MATPOWER/mptest
[6]: https://github.com/MATPOWER/mips
[7]: docs/MP-Opt-Model-manual.pdf
[8]: CONTRIBUTING.md
[9]: LICENSE
[10]: CITATION
[11]: https://doi.org/10.5281/zenodo.???????
[12]: https://github.com/MATPOWER/matpower/blob/master/lib/opf_setup.m
[13]: https://github.com/MATPOWER/matpower/blob/master/lib/opf_execute.m