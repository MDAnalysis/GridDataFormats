# Maintainer documentation

For additional documentation see the [Developer Guide: Release
Management](https://github.com/MDAnalysis/GridDataFormats/wiki/Developer-Guide#release-management)
page. This file gives a brief reminder of what maintainers need to do
for new releases.

1. create a release on GitHub using tag `<major>.<minor>.<patch>`.
1. release on PyPi under https://pypi.org/project/GridDataFormats
1. release on conda-forge https://anaconda.org/conda-forge/griddataformats

## GitHub release

* We use [semantic versioning](https://semver.org) MAJOR.MINOR.PATCH
  (i.e., briefly, major revision changes whenever the API changes in
  backwards-incompatible manner, MINOR changes for new features, PATCH
  changes for fixes that correct functionality; as long as MAJOR == 0,
  we can break the API with increasing MINOR.)
  
* Releases are cut from the master branch and tagged with
  *MAJOR.MINOR.PATCH* (note: the release tag *determines* the tag
  because we use
  [versioneer](https://github.com/warner/python-versioneer/blob/master/INSTALL.md#post-installation-usage),
  which obtains the release number from the git tag). We do from the
  master branch:
  
    1. `git tag <major>.<minor>.<patch>`
    1. `git push --tags`

* Create a
  [GitHub release](https://github.com/MDAnalysis/GridDataFormats/releases)
  from the tag and name it `v<major>.<minor>.<patch>` and add a short description.

## PyPi release

Upload to PyPi can be done by PyPi maintainers and requires `twine`:

 
1. `python setup.py sdist bdist_wheel`
2. `twine upload dist/*`

## Update Conda-forge package

*After* a PyPi release update the conda-forge package. For this do the following
on a local checkout of the package
[feedstock](https://github.com/conda-forge/griddataformats-feedstock)

1. create a new branch
1. conda smithy rerender
1. update the sha256 in the `meta.yaml` (see the [PyPi downloads
   page](https://pypi.org/project/GridDataFormats/#files) for the
   sha256 of the tar.gz file)
1. update version number

Afterwards upload the new branch to your **own fork** of the feedstock and
generate a PR. Once all tests pass merge the PR and the package will be
published.


### Update package on MDAnalysis channel


Don't. We don't have the man power to update all the dependencies we need in the
channel ourselves. Relying on conda-forge is more reliant.



## Documentation

Documentation is automatically generated on Travis CI and pushed to
the gh-pages branch and appears at https://www.mdanalysis.org/GridDataFormats/.


There is also alternative documentation on ReadTheDocs
https://griddataformats.readthedocs.io/, which automatically rebuilds.
