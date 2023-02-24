#!/usr/bin/sh

if [ -z "$1" ]; then
    echo "No input was given"
    msg='quick push'
else
    echo "Input was given"
    msg=$1
fi

git pull
git status
git add --all
git commit -m "$msg"
git push -u origin HEAD:main
