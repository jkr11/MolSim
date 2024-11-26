#!/bin/bash

INPUT_XSD="../src/io/file/in/xml/input.xsd"
INPUT_DIR=$(dirname "$INPUT_XSD")
echo "Running xsdcxx tree compiler on $INPUT_XSD"
xsdcxx cxx-tree --output-dir "$INPUT_DIR" --std c++11 --generate-doxygen "$INPUT_XSD"

if [ $? -ne 0 ]; then
  echo "Error: xsdcxx failed to compile the XSD file."
  exit 1
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

echo "clang-format applied successfully to generated files."

exit 0
