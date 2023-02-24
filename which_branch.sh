#!/bin/sh

git branch     2>/dev/null | grep "^*" | colrm 1 2
