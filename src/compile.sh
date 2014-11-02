#!/bin/bash
set -e

matlab_dir="" # set this if matlab is not on your path

# compile mex
mcmd="compileMex; quit;"
"$matlab_dir"matlab -nodisplay -nojvm -r "${mcmd}" > compile.log

exit
