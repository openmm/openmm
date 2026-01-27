#!/bin/bash
# Download OMol25 dataset from HuggingFace
# 
# Prerequisites:
# 1. Request access at https://huggingface.co/facebook/OMol25
# 2. Wait for approval (check your email)
# 3. Run this script

set -e

echo "==================================="
echo "OMol25 Dataset Download Script"
echo "==================================="
echo ""

# Check if huggingface_hub is installed
if ! python -c "import huggingface_hub" 2>/dev/null; then
    echo "Installing huggingface_hub..."
    pip install huggingface_hub[cli]
fi

# Check if user is logged in
if ! huggingface-cli whoami &>/dev/null; then
    echo "Please login to HuggingFace:"
    huggingface-cli login
fi

# Set download directory
DATA_DIR="/media/extradrive/Trajectories/openmm/data/omol25"
mkdir -p "$DATA_DIR"

echo ""
echo "Downloading OMol25 train_4M subset to $DATA_DIR"
echo "This will download ~100s of GB, may take several hours..."
echo ""

# Download train_4M subset
huggingface-cli download facebook/OMol25 \
  --repo-type dataset \
  --include "train_4M/*" \
  --local-dir "$DATA_DIR" \
  --resume-download

echo ""
echo "Downloading validation set..."
huggingface-cli download facebook/OMol25 \
  --repo-type dataset \
  --include "val/*" \
  --local-dir "$DATA_DIR" \
  --resume-download

echo ""
echo "==================================="
echo "Download complete!"
echo "Data location: $DATA_DIR"
echo ""
echo "Directory structure:"
ls -lh "$DATA_DIR"
echo "==================================="
