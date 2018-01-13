Update Conda-forge package
--------------------------------

*After* a PyPi release update the conda-forge package. For this do the following
on a local checkout of the package
[feedstock](https://github.com/conda-forge/griddataformats-feedstock)

1. create a new branch
1. conda smithy rerender
1. update the sha256 in the `meta.yaml`
1. update version number

Afterwards upload the new branch to your **own fork** of the feedstock and
generate a PR. Once all tests pass merge the PR and the package will be
published.


Update package on MDAnalysis channel
---------------------------------------------

Don't. We don't have the man power to update all the dependencies we need in the
channel ourselves. Relying on conda-forge is more reliant.
