#!/bin/#!/usr/bin/env bash

s=$(printf "%-80s" "*")
echo "${s// /*}"

### Find COPSS Source folder
if [ "$COPSS_DIR" != "" ]; then
  echo "COPSS directory found at $COPSS_DIR"
else
  echo "COPSS source library not found. Exiting..."
  echo "Suggestion: please add 'export COPSS_DIR=[/path/to/COPSS/DIRECTORY]' to\
your environment, e.g., .bashrc file"
  exit -1
fi
cd $COPSS_DIR/src


action=""
package=""
### Find all flags
while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo "options:"
                        echo "-h, --help:   show brief help"
                        echo "-p, --package [PACKAGE]: Please specify an package to compile, options are: POINTPARTICLE, RIGIDPARTICLE"
                        echo "-a, --action [ACTION]:   specify additional actions, options are: clean_first, clean_only"
                        exit -1
                        ;;
                -p|--packa*)
                        shift
                        if test $# -gt 0; then
                                if [ $1 = "POINTPARTICLE" ] || [ $1 = "RIGIDPARTICLE" ]; then
                                  package=$1
                                else
                                  echo "Package $1 not supported. Exiting ..."
                                  exit -1
                                fi
                        fi
                        shift
                        ;;
                -a|--act*)
                        shift
                        if test $# -gt 0; then
                                if [ $1 = "clean_first" ] || [ $1 = "clean_only" ]; then
                                  action=$1
                                else
                                  echo "Action $1 not supported. Exiting ..."
                                  exit -1
                                fi
                        fi
                        shift
                        ;;
                *)
                        echo "Flag '$1' not supported. Use -h for help"
                        exit -1
                        ;;
        esac
done

function clean() {
  make package=POINTPARTICLE clean
  make package=RIGIDPARTICLE clean
}

function compile(){
  if [[ -z "$package" ]]; then
    echo "Compiling package is not specified, set it to POINTPARTICLE by default"
    package="POINTPARTICLE"
  fi
  echo "Start compiling package $package ..."
  echo "${s// /*}"
  if make package=$package; then
      # export exectuble
      copss_exec="${COPSS_DIR}/src/copss-$package-opt"
      echo "${s// /*}"
      echo "Compilation is done, the compiled executable is at $copss_exec"
      echo "${s// /*}"
  else
      echo "Error: Compilation failed"
      exit -1
  fi
}

if [[ $action = "clean_only" ]]; then
  clean
elif [[ $action = "clean_first" ]]; then
  clean
  compile
else
  compile
fi
