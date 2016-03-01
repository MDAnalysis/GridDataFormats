How to build package
--------------------
To build and run tests on the package for python 2/3 use the following command

    `conda build . --python 3.5 --python 2.7`

To build the package the `meta.yaml` is building from the current source
directory

If the build fails because `tempdir` is missing add the MDAnalysis channel
to your `.condarc`.

Prepare for a release
---------------------

Update the version string in `meta.yaml`

Pushing to official MDAnalysis channel
--------------------------------------

To upload to your own channel `anaconda upload <where conda build the package>`.
To upload to the MDAnalysis channel add the flag `-u MDAnalysis`
