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

* This will automatically trigger a github action (named 'deploy') to using pypa's `build` tool to create a tarball & pure Python wheel and upload it to https://test.pypi.org/project/GridDataFormats

* Once uploaded to testpypi, please check locally that the testpypi build is working as intended. In a clean environment do:

    1. `pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple GridDataFormats=="version number"`
    2. `pip install pytest`
    3. `pytest --pyargs gridData`
  
* Create a
  [GitHub release](https://github.com/MDAnalysis/GridDataFormats/releases)
  from the tag and name it `<major>.<minor>.<patch>` and add a short description.
  
* The GitHub release triggers the `deploy` action to deploy the tarball and wheel to the standard PyPI repository.


## PyPi release

A GitHub release is automatically deployed to PyPI.

## Update Conda-forge package

*After* a PyPI release update the conda-forge package [feedstock](https://github.com/conda-forge/griddataformats-feedstock).

### Automatic

1. Wait for the *regro-cf-autotick-bot* to create a PR, based on the PyPI release (can take a few hours).
2. review the PR
3. merge the PR
4. conda packages will be built


### Manual

Manual updates are rarely necessary. 

If necessary do the following on a local checkout of the package
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



## Documentation

Documentation is automatically generated in CI and pushed to
the gh-pages branch and appears at https://www.mdanalysis.org/GridDataFormats/.


There is also alternative documentation on ReadTheDocs
https://griddataformats.readthedocs.io/, which automatically rebuilds.
