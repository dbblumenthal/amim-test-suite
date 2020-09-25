
=== DESCRIPTION ===

This directory contains an SVN checkout of the PinnacleZ repository, as retrieved on 2011-04-19.

See the "GetRepositories.sh" script for the retrieval URL, etc.

PinnacleZ depends on the libraries 'oiler' and 'modlab'. It comes with a pre-compiled jar-file
of those. However, we downloaded their source code as well for completeness, although we do not
re-compile them.

The last change to the PinnacleZ repository was on the 2nd of November, 2007, so it is not
under active development.

=== PATCHES ===

In the course of our work we implemented a number of patches:

* patches/AddOptionForNormalizingInput.patch

  PinnacleZ will normally start by Z-normalizing the expression data, i.e., translating and scaling
  the expression values per gene, such that mu = 0, sigma = 1.

  This patch introduces a "-z" option that skips this renormalization step. It is most useful in
  case we want to train using a subset of a previously normalized dataset (e.g., during cross-validation).

* patches/InitializeRandomSeeds.patch

  PinnacleZ doesn't initialize its random seeds. This patch fixes that, in an attempt to make PinnacleZ
  deterministic. To make PinnacleZ fully deterministic, the next patch is also required.

* patches/ForceSingleThreaded.patch

  Even when the random seeds are fixed, PinnacleZ is not deterministic. This is because a single random
  generator is used from multiple threads of execution.

  This patch forces PinnacleZ to run single-threaded. In combination with the previous patch, this makes
  PZ's behavior fully deterministic and predictable.

  Unfortunately, this makes PinnacleZ's network search very slow.

=== ISSUES ===

We discovered some issues in the PinnacleZ code. These are documented in the README.txt that goes with the
"ChuangFeatureExtractor.py". Please look there for more info.

=== BUILDING ===

To build PinnacleZ, please execute the "MakePinnacleZ.sh" script.

This unpacks the source, applies patches, and runs "ant".

The end result is three versions of PinnacleZ:

* pinnaclez-ORIGINAL.jar

  A pristine version of PinnacleZ; no patches applied.

* pinnaclez-DETERMINISTIC.jar

  A deterministic (but slow) version of PinnacleZ, that has these patches applied:

  - patches/InitializeRandomSeeds.patch
  - patches/ForceSingleThreaded.patch

  This version was not used in our work.

* pinnaclez-OPTIONAL_NORMALIZE.jar

  A version of PinnacleZ that adds a "-z" option that forgoes the normalization of expression data.
  It has the following patch applied:

  - patches/AddOptionForNormalizingInput.patch
