#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <version> [output_folder]"
    echo "Example: $0 5.76-107.0 interproscan_db/"
    exit 1
fi

VERSION="$1"
OUTPUT="${2:-interproscan_db}"

REPO="https://github.com/ebi-pf-team/interproscan.git"
FOLDER="core/jms-implementation/support-mini-x86-32/data"

ORIGINAL_DIR=$(pwd)

git clone --depth 1 --filter=blob:none --sparse "$REPO"
cd interproscan
git sparse-checkout set "$FOLDER"
git checkout "$VERSION"

mkdir -p "$ORIGINAL_DIR/$OUTPUT"
mv "$FOLDER" "$ORIGINAL_DIR/$OUTPUT"
cd "$ORIGINAL_DIR"
rm -rf interproscan
