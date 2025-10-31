#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script: post-crop.sh
# Purpose: Crop PNG images to remove unnecessary borders
# Creator: Raphael Schiedung
# -----------------------------------------------------------------------------
#
# -----------------------------
# User-configurable parameters
# -----------------------------
INPUT_DIR="Images"
OUTPUT_DIR="Images/Cropped"
# -----------------------------

set -euo pipefail

mkdir -p "$OUTPUT_DIR"
# Create output directory if it doesn't exist

# First, find max width and height after trimming
MAX_WIDTH=0
MAX_HEIGHT=0

echo "Determine Maximum cropped size..."

for img in "$INPUT_DIR"/*.png; do
    # Get dimensions after trimming
    read WIDTH HEIGHT <<< $(convert "$img" -trim -format "%w %h" info:)
    if (( WIDTH > MAX_WIDTH )); then MAX_WIDTH=$WIDTH; fi
    if (( HEIGHT > MAX_HEIGHT )); then MAX_HEIGHT=$HEIGHT; fi
done

echo "Maximum cropped size: ${MAX_WIDTH}x${MAX_HEIGHT}"

# Now crop and pad all images to uniform size
for img in "$INPUT_DIR"/*.png; do
    base=$(basename "$img")
    convert "$img" -trim \
        -background transparent -gravity center \
        -extent ${MAX_WIDTH}x${MAX_HEIGHT} \
        "$OUTPUT_DIR/$base"
done

echo "✅ All images cropped and padded to uniform size in $OUTPUT_DIR"
