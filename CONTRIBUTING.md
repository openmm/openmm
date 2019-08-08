## How to Contribute to OpenMM Development

We welcome anyone who wants to contribute to the project, whether by adding a feature,
fixing a bug, or improving documentation.  The process is quite simple.

First, it is always best to begin by opening an issue on Github that describes the change you
want to make.  This gives everyone a chance to discuss it before you put in a lot of work.
For bug fixes, we will confirm that the behavior is actually a bug and that the proposed fix
is correct.  For new features, we will decide whether the proposed feature is something we
want and discuss possible designs for it.

Once everyone is in agreement, the next step is to
[create a pull request](https://help.github.com/en/articles/about-pull-requests) with the code changes.
For larger features, feel free to create the pull request even before the implementaton is
finished so as to get early feedback on the code.  When doing this, put the letters "WIP" at
the start of the title of the pull request to indicate it is still a work in progress.

For new features, consult the [New Feature Checklist](https://github.com/openmm/openmm/wiki/Checklist-when-adding-a-new-feature),
which lists various items that need to be included before the feature can be merged (documentation,
tests, serialization, support for all APIs, etc.).  Not every item is necessarily applicable to
every new feature, but usually at least some of them are.

The core developers will review the pull request and may suggest changes.  Simply push the
changes to the branch that is being pulled from, and they will automatically be added to the
pull request.  In addition, the full test suite is automatically run on every pull request,
and rerun every time a change is added.  Once the tests are passing and everyone is satisfied
with the code, the pull request will be merged.  Congratulations on a succesful contribution!