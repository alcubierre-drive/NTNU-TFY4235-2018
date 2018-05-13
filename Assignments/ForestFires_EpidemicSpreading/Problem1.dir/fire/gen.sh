#!/bin/bash

PREFIX="$1"
OUT="$2"

ffmpeg -y -framerate 10 -i "$PREFIX"'_%03d.png' \
    -c:v libx264 -f mp4 -crf 0 "$OUT"
