#!/bin/bash
set -e  # Exit on error

# Pairs to process
PAIRS=(
  "avatar-i.png avatar.png composite.png"
  "avatar-i-dm.png avatar-dm.png composite-dm.png"
)

for pair in "${PAIRS[@]}"; do
    read left right output <<< "$pair"

    echo "Composing: $left + $right → $output"

    # Get height of left image
    height=$(identify -format "%h" "$left")

    # Resize right image to match height
    convert "$right" -resize x${height} "resized_$right"

    # Combine side by side
    convert +append "$left" "resized_$right" "$output"

    # Cleanup
    rm "resized_$right"
done

echo "Done!"
