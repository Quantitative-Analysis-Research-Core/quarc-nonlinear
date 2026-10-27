#!/bin/zsh
# Copy the C++ core (cpp/) into the R and Python packages.
#
# R and Python build from their own package directory, so each carries a copy
# of the core. Never edit the copies: edit cpp/, then run this. The tests in
# both packages fail if a copy has drifted from cpp/.

set -eu
setopt null_glob
repo=${0:A:h:h}
for dest in $repo/rpkg/src $repo/pypkg/src/core; do
  mkdir -p $dest/quarc
  rm -f $dest/core_*.cpp $dest/quarc/*.hpp
  cp $repo/cpp/include/quarc/*.hpp $dest/quarc/
  for f in $repo/cpp/src/*.cpp; do cp $f $dest/core_${f:t}; done
done
echo "synced cpp/ into rpkg/src and pypkg/src/core"
