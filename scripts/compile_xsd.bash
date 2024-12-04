#!/bin/bash

INPUT_XSD="../src/io/file/in/xml/input.xsd"
CHECKPOINT_XSD="../src/io/file/out/checkpoint-schema.xsd"
INPUT_DIR=$(dirname "$INPUT_XSD")
CHECKPOINT_DIR=$(dirname "$CHECKPOINT_XSD")
echo "Running xsdcxx tree compiler on $INPUT_XSD"
xsdcxx cxx-tree --output-dir "$INPUT_DIR" "$INPUT_XSD"

if [ $? -ne 0 ]; then
  echo "Error: xsdcxx failed to compile the input XSD file."
  exit 1
fi

xsdcxx cxx-tree --generate-serialization --output-dir "$CHECKPOINT_DIR" "$CHECKPOINT_XSD"

if [ $? -ne 0 ]; then
  echo "Error: xsdcxx failed to compile the XSD file."
  exit 1
else
  echo "Compiled Checkpoints"
fi

echo "xsdcxx compilation completed successfully."

echo "Applying clang-format to generated files..."
for file in "${INPUT_XSD%.xsd}.hxx" "${INPUT_XSD%.xsd}.cxx"; do
  if [ -f "$file" ]; then
    clang-format -i "$file"
    echo "Formatted: $file"
  else
    echo "Warning: File not found: $file"
  fi
done

for file in "${CHECKPOINT_XSD%.xsd}.hxx" "${CHECKPOINT_XSD%.xsd}.cxx"; do
  if [ -f "$file" ]; then
    clang-format -i "$file"
    echo "Formatted: $file"
  else
    echo "Warning: File not found: $file"
  fi
done

echo "clang-format applied successfully to generated files."

exit 0
