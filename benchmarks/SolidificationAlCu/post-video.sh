#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: post-video.sh
# Purpose: Compile PNG images into a WebM video
# Creator: Raphael Schiedung
# -----------------------------------------------------------------------------
# -----------------------------
# User-configurable parameters
# -----------------------------
FPS=20                            # Frames per second
INPUT_DIR="Images/Cropped"        # Input image folder
OUTPUT_DIR="Videos"               # Output video folder
INPUT_NAME="SolidPhaseFraction"   # Base name of input files
INPUT_PATTERN="${INPUT_DIR}/${INPUT_NAME}.%04d.png"
OUTPUT_PATTERN="${OUTPUT_DIR}/${INPUT_NAME}_fps${FPS}.webm"
# -----------------------------

set -euo pipefail

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Pick the first file to probe dimensions
FIRST_FRAME=$(printf "${INPUT_DIR}/${INPUT_NAME}.%04d.png" 0)

if [[ ! -f "${FIRST_FRAME}" ]]; then
  echo "❌ Error: No input image found at ${FIRST_FRAME}"
  exit 1
fi

# Get width and height using ffprobe
WIDTH=$(ffprobe -v error -select_streams v:0 -show_entries stream=width  -of csv=p=0 "${FIRST_FRAME}")
HEIGHT=$(ffprobe -v error -select_streams v:0 -show_entries stream=height -of csv=p=0 "${FIRST_FRAME}")

# Number of frames
NUM_FRAMES=$(ls "$INPUT_DIR"/*.png | wc -l)
DURATION=$(echo "$NUM_FRAMES / $FPS" | bc -l)

echo "Number of frames: $NUM_FRAMES"
echo "Video duration: $DURATION s"
echo "🖼️  Detected resolution: ${WIDTH}x${HEIGHT}"
echo "🎬 Creating video at ${FPS} FPS → ${OUTPUT_PATTERN}"

# Run ffmpeg with white background of correct size
ffmpeg -y -framerate "${FPS}" -i "${INPUT_PATTERN}" \
  -f lavfi -i "color=white:s=${WIDTH}x${HEIGHT}:d=${DURATION}" \
  -filter_complex "[1][0]overlay=format=auto" \
  -c:v libvpx-vp9 -crf 25 -b:v 0 "${OUTPUT_PATTERN}"

echo "✅ Video created: ${OUTPUT_PATTERN}"
xdg-open ${OUTPUT_PATTERN}
