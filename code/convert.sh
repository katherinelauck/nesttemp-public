#!/bin/bash

for FILE in *; do
    if [[ $FILE == *.h264 ]]
    then
        MP4Box -fps 30 -add $FILE "$FILE.mp4"
    fi
done;
