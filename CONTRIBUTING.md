## How to Contribute to OpenMM Development

We welcome anyone who wants to contribute to the project, whether by adding a feature,
fixing a bug, or improving documentation.  The process is quite simple.

First, it is always best to begin by opening an issue on GitHub that describes the change you
want to make.  This gives everyone a chance to discuss it before you put in a lot of work.
For bug fixes, we will confirm that the behavior is actually a bug and that the proposed fix
is correct.  For new features, we will decide whether the proposed feature is something we
want and discuss possible designs for it.

Once everyone is in agreement, the next step is to
[create a pull request](https://help.github.com/en/articles/about-pull-requests) with the code changes.
For larger features, feel free to create the pull request even before the implementation is
finished so as to get early feedback on the code.  When doing this, put the letters "WIP" at
the start of the title of the pull request to indicate it is still a work in progress.

> [!IMPORTANT]
> Be sure to read our [AI Policy](AI_POLICY.md), which describes requirements and guidelines for
> use of artificial intelligence in OpenMM development.

Keep the following guidelines and tips in mind when proposing changes to OpenMM:

* Each distinct bug fix, feature implementation, or other improvement should be
  submitted in a separate pull request.  A pull request fixing a particular bug
  should not contain unrelated code changes, which should instead be submitted
  separately if they are relevant to another issue or improvement.

* For new features, consult the [New Feature Checklist](https://github.com/openmm/openmm/wiki/Checklist-for-Adding-a-New-Feature),
which lists various items that need to be included before the feature can be merged (documentation,
tests, serialization, support for all APIs, etc.).  Not every item is necessarily applicable to
every new feature, but usually at least some of them are.

* If you are proposing an optimization to OpenMM intended to improve simulation
  performance, be sure to use OpenMM's [standard benchmark suite](https://github.com/openmm/openmm/blob/master/examples/benchmarks/benchmark.py),
  as applicable, to assess the impact of your changes.

* If you are planning on making substantial changes to OpenMM's C++ code,
  familiarize yourself with the [OpenMM Developer Guide](https://docs.openmm.org/latest/developerguide/index.html)
  which contains useful information about its internal architecture.

Following these guidelines will help the core developers review the pull request
and provide feedback more effectively.  The developers may suggest changes to
your pull request.  After you make the requested changes, simply push them to
the branch that is being pulled from, and they will automatically be added to the
pull request.  In addition, the full test suite is automatically run on every pull request,
and rerun every time a change is added.  Once the tests are passing and everyone is satisfied
with the code, the pull request will be merged.  Congratulations on a successful contribution!
