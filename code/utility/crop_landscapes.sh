#!/bin/bash

CURR_DIR=`pwd`

cd $CURR_DIR
cd ./figures/landscapes/gmt_landscapes

PNG_LIST=`ls  */*.png`

echo "Cropping landscapes"

for f in $PNG_LIST; do
 #   convert ${f} -crop 2200x1200+500+600 ${f}
    convert ${f} -gravity center -crop 1600x1100-140+150 ${f}
done
