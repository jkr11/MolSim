#!/bin/bash

SRC_DIR="../src"

CLANG_FORMAT="clang-format"

if ! command -v $CLANG_FORMAT &>/dev/null; then
  echo "Error: clang-format s not installed"
  exit 1
fi

echo "Applying clang-format to all files in $SRC_DIR..."
find "$SRC_DIR" \( -name "*.cpp" -o -name "*.h" -o -name "*.cxx" -o -name "*.hxx" \) -exec $CLANG_FORMAT -i {} +

echo "Formatted"