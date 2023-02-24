#!/usr/bin/sh

if [ -z "$2" ]; then
  min=1
  max=100
else
  min=$2
  max=$2
fi

echo "will search all dirs with $1 files and mindepth/maxdepth = $min/$max "
find . -maxdepth $max -mindepth $min -type d -exec bash -c "echo -ne '{} '; ls '{}' | wc -l" \; | awk '$NF'==$1
