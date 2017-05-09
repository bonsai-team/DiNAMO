#!/bin/bash
set -ev
if [ "$TRAVIS_BRANCH" = "unstable" ]; then
	make
fi
