#!/usr/bin/sh

cat tmp  | sort | uniq | grep -v -e 'c1(ccccc1)' -e 'c1ccccc1' -e 'c1ccc(cc1)' -e 'C1=CC=CC=C1' -e 'C1(=CC=CC=C1)' -e '[C]1=CC=CC=C1' -e '\.'
