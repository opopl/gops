#!/bin/bash
# SVN VERSION SCRIPT
# To be run on compilation - returns the current version of the code.
# This script should be called from within the Makefile, and the results passed to the preprocessor
# For example:
# DEFS+=-DSVNVERSION="`./version.sh`"

# First, we try to run svnversion to get the version number
svnversion . | sed 's/.*://' | sed 's/M//' > version.tmp
# svnversion returns 'exported' if you have it installed but run it in a non svn directory.
# If svnversion is not installed at all, the file will be empty. In either case, we want to
# use the version number recorded previously. 
if [ "`cat version.tmp`" != "exported" ] && [ "`cat version.tmp`" != "" ]; then 
# If we are working with svn however, we want to overwrite the version number with the current one.
    cp version.tmp SVNREV 
fi
# Remove the temporary file
rm version.tmp
# Return the version from the VERSION file
cat SVNREV
