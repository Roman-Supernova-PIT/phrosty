Contributing to phrosty
================================================

Reporting Issues
----------------

When opening an issue to report a problem, please try to provide a minimal code
example that reproduces the issue along with details of the operating
system and the Python, NumPy, and `phrosty` versions you are using.

Contributing Code and Documentation
-----------------------------------

So you are interested in contributing to the phrosty Project?  Excellent!
We love contributions! We are open source, built on open source, and
we'd love to have you hang out in our community.


Anti Imposter Syndrome Reassurance
----------------------------------

We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code or
documentation, you can contribute code to open source.
Contributing to open source projects is a fantastic way to advance one's coding
and open source workflow skills. Writing perfect code isn't the measure of a good
developer (that would disqualify all of us!); it's trying to create something,
making mistakes, and learning from those mistakes. That's how we all improve,
and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This text was originally written by
[Adrienne Lowe](https://github.com/adriennefriend) for a
[PyCon talk](https://www.youtube.com/watch?v=6Uj746j9Heo), and was adapted by
Astropy based on its use in the README file for the
[MetPy project](https://github.com/Unidata/MetPy).


How to Contribute, Best Practices
---------------------------------

Most contributions are done via [pull
requests](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests)
from GitHub users' forks of phrosty If you are new to this
style of development, you will want to read over our development workflow in the docs.

When you open a pull request (which should be opened against the ``main``
branch, not against any of the other branches), please make sure to
include the following:

- **Code**: the code you are adding, which should follow
  our coding guidelines

- **Tests**: these are usually tests to ensure code that previously
  failed now works (regression tests), or tests that cover as much as possible
  of the new functionality to make sure it does not break in the future and
  also returns consistent results on all platforms (since we run these tests on
  many platforms/configurations). For more information about how to write
  tests, see our testing guidelines.

- **Documentation**: if you are adding new functionality, be sure to include a
  description in the main documentation (in ``docs/``).


- **Performance improvements**: if you are making changes that impact package
  performance, consider adding a reference for the performance benchmark 


  If you are opening a pull request you may not know
  the PR number yet, but you can add it once the pull request is open.



Checklist for Contributed Code
------------------------------

Before being merged, a pull request for a new feature will be reviewed to see if
it meets the following requirements. If you are unsure about how to meet all of these
requirements, please submit the PR and ask for help and/or guidance. 

**Scientific Quality** (when applicable)
  * Is the submission relevant to astronomy?
  * Are references included to the origin source for the algorithm?
  * Does the code perform as expected?
  * Has the code been tested against previously existing implementations?

**Code Quality**
  * Are the coding guidelines followed?
  * Is the code compatible with the supported versions of Python?
  * Are there dependencies other than the run-time dependencies listed in pyproject.toml?
    * Is the package importable even if the C-extensions are not built?
    * Are additional dependencies handled appropriately?
    * Do functions that require additional dependencies raise an `ImportError`
      if they are not present?

**Testing**
  * Are the testing guidelines followed?
  * Are the inputs to the functions sufficiently tested?
  * Are there tests for any exceptions raised?
  * Are there tests for the expected performance?
  * Are the sources for the tests documented?
  * Have tests that require an optional dependency been marked as such?
  * Does ``tox -e test`` run without failures?

**Documentation**
  * Are the documentation guidelines followed?
  * Is there a docstring in [numpydoc format](https://numpydoc.readthedocs.io/en/latest/format.html) in the function describing:
    * What the code does?
    * The format of the inputs of the function?
    * The format of the outputs of the function?
    * References to the original algorithms?
    * Any exceptions which are raised?
    * An example of running the code?
  * Is there any information needed to be added to the docs to describe the
    function?
  * Does the documentation build without errors or warnings?

**License**
  * Are there any conflicts with this code and existing codes?

**Package requirements**
  * Do all the GitHub Actions tests pass? If not, are they allowed to fail?
  * If applicable, has an entry been added into the changelog?
  * Can you check out the pull request and repeat the examples and tests?

Other Tips
----------

- Behind the scenes, we conduct a number of tests or checks with new pull requests.
  To prevent the automated tests from running, you can add ``[ci skip]``
  to your commit message. This is useful if your PR is a work in progress (WIP) and
  you are not yet ready for the tests to run. For example:

      $ git commit -m "WIP widget [ci skip]"

  - If you already made the commit without including this string, you can edit
    your existing commit message by running:

        $ git commit --amend

- If your commit makes substantial changes to the documentation but none of
  those changes include code snippets, then you can use ``[ci skip]``,
  which will skip all CI except where the documentation is built.

- When contributing trivial documentation fixes (i.e., fixes to typos, spelling,
  grammar) that don't contain any special markup and are not associated with
  code changes, please include the string ``[ci skip]`` in your commit
  message.

      $ git commit -m "Fixed typo [ci skip]"

- ``[ci skip]`` and ``[skip ci]`` are the same and can be used interchangeably.