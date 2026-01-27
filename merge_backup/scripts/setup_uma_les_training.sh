#!/bin/bash
# Complete setup script for UMA+LES training on OMol25
#
# This script:
# 1. Checks dependencies
# 2. Downloads OMol25 data (if needed)
# 3. Tests configuration
# 4. Provides training commands

set -e

echo "=============================================="
echo "  UMA+LES Training Setup"
echo "=============================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Directories
REPO_DIR="/media/extradrive/Trajectories/openmm"
FAIRCHEM_DIR="$REPO_DIR/fairchem"
DATA_DIR="$REPO_DIR/data/omol25"
CONFIG_FILE="$FAIRCHEM_DIR/configs/uma/lr/uma_sm_omol25_les_train.yaml"

# Step 1: Check dependencies
echo "Step 1: Checking dependencies..."
echo "-------------------------------------------"

# Check Python
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}✗ Python 3 not found${NC}"
    exit 1
fi
echo -e "${GREEN}✓ Python 3: $(python3 --version)${NC}"

# Check PyTorch
if python3 -c "import torch" 2>/dev/null; then
    TORCH_VERSION=$(python3 -c "import torch; print(torch.__version__)")
    echo -e "${GREEN}✓ PyTorch: ${TORCH_VERSION}${NC}"
    
    # Check CUDA
    if python3 -c "import torch; exit(0 if torch.cuda.is_available() else 1)" 2>/dev/null; then
        CUDA_VERSION=$(python3 -c "import torch; print(torch.version.cuda)")
        GPU_COUNT=$(python3 -c "import torch; print(torch.cuda.device_count())")
        echo -e "${GREEN}✓ CUDA: ${CUDA_VERSION} (${GPU_COUNT} GPUs)${NC}"
    else
        echo -e "${YELLOW}! CUDA not available (will use CPU)${NC}"
    fi
else
    echo -e "${RED}✗ PyTorch not installed${NC}"
    echo "  Install with: pip install torch"
    exit 1
fi

# Check fairchem-core
if python3 -c "import fairchem.core" 2>/dev/null; then
    echo -e "${GREEN}✓ fairchem-core installed${NC}"
else
    echo -e "${YELLOW}! fairchem-core not installed${NC}"
    echo "  Installing from ${FAIRCHEM_DIR}..."
    cd "$FAIRCHEM_DIR"
    pip install -e packages/fairchem-core[dev] || {
        echo -e "${RED}✗ Failed to install fairchem-core${NC}"
        exit 1
    }
    echo -e "${GREEN}✓ fairchem-core installed${NC}"
fi

# Check if on les_branch
echo ""
echo "Step 2: Checking git branch..."
echo "-------------------------------------------"
cd "$FAIRCHEM_DIR"
CURRENT_BRANCH=$(git branch --show-current)
if [ "$CURRENT_BRANCH" = "les_branch" ]; then
    echo -e "${GREEN}✓ On les_branch${NC}"
else
    echo -e "${YELLOW}! Currently on branch: ${CURRENT_BRANCH}${NC}"
    echo "  Switching to les_branch..."
    git checkout les_branch || {
        echo -e "${RED}✗ Failed to switch to les_branch${NC}"
        exit 1
    }
    echo -e "${GREEN}✓ Switched to les_branch${NC}"
fi

# Step 3: Check data
echo ""
echo "Step 3: Checking OMol25 dataset..."
echo "-------------------------------------------"
if [ -d "$DATA_DIR/train_4M" ]; then
    FILE_COUNT=$(find "$DATA_DIR/train_4M" -name "*.aselmdb" | wc -l)
    echo -e "${GREEN}✓ OMol25 data found: ${FILE_COUNT} files in train_4M${NC}"
else
    echo -e "${YELLOW}! OMol25 data not found${NC}"
    echo "  Please download the dataset:"
    echo "  1. Request access at: https://huggingface.co/facebook/OMol25"
    echo "  2. Run: $REPO_DIR/scripts/download_omol25.sh"
    echo ""
    echo "  Continuing anyway (can test config without data)..."
fi

# Step 4: Verify configurations
echo ""
echo "Step 4: Verifying configurations..."
echo "-------------------------------------------"

CONFIGS=(
    "$FAIRCHEM_DIR/configs/uma/lr/dataset/omol25.yaml"
    "$FAIRCHEM_DIR/configs/uma/lr/tasks/omol25.yaml"
    "$FAIRCHEM_DIR/configs/uma/lr/element_refs/omol25_element_reference.yaml"
    "$FAIRCHEM_DIR/configs/uma/lr/backbone/uma_sm_les.yaml"
    "$CONFIG_FILE"
)

ALL_CONFIGS_EXIST=true
for config in "${CONFIGS[@]}"; do
    if [ -f "$config" ]; then
        echo -e "${GREEN}✓ $(basename $config)${NC}"
    else
        echo -e "${RED}✗ Missing: $config${NC}"
        ALL_CONFIGS_EXIST=false
    fi
done

if [ "$ALL_CONFIGS_EXIST" = false ]; then
    echo -e "${RED}✗ Some configuration files are missing${NC}"
    exit 1
fi

# Step 5: Test configuration
echo ""
echo "Step 5: Testing configuration (dry run)..."
echo "-------------------------------------------"
cd "$FAIRCHEM_DIR"

echo "Running: fairchem -c $CONFIG_FILE job.debug=True epochs=1 batch_size=2 --help"
echo ""

# Just check if the command would work
if command -v fairchem &> /dev/null; then
    echo -e "${GREEN}✓ fairchem CLI available${NC}"
else
    echo -e "${YELLOW}! fairchem CLI not in PATH${NC}"
    echo "  You may need to reinstall fairchem-core"
fi

# Summary and next steps
echo ""
echo "=============================================="
echo "  Setup Complete!"
echo "=============================================="
echo ""
echo "Configuration file: $CONFIG_FILE"
echo ""
echo "Next steps:"
echo ""
echo "1. Test configuration (debug mode):"
echo "   cd $FAIRCHEM_DIR"
echo "   fairchem -c configs/uma/lr/uma_sm_omol25_les_train.yaml \\"
echo "     job.debug=True epochs=1 batch_size=2"
echo ""
echo "2. Run full training (single GPU):"
echo "   fairchem -c configs/uma/lr/uma_sm_omol25_les_train.yaml"
echo ""
echo "3. Run with multiple GPUs (e.g., 4 GPUs):"
echo "   fairchem -c configs/uma/lr/uma_sm_omol25_les_train.yaml \\"
echo "     job.scheduler.ranks_per_node=4"
echo ""
echo "4. Monitor training:"
echo "   - Logs: $REPO_DIR/training_runs/uma_les_omol25/"
echo "   - W&B: wandb login && check dashboard"
echo ""
echo "5. Validate trained model:"
echo "   python $REPO_DIR/scripts/validate_uma_les.py \\"
echo "     --checkpoint training_runs/uma_les_omol25/checkpoints/final/inference_ckpt.pt \\"
echo "     --test-data data/omol25/test/test.aselmdb"
echo ""
echo "=============================================="
