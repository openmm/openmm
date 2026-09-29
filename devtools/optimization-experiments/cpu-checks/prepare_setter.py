"""Extract actual baseline/current setter bodies; does not load OpenMM or CUDA."""
from pathlib import Path
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
baseline = '0c4bcaba734ea574f52f5e09ffaaf4d6fe6109c2'
name = 'platforms/common/src/CommonKernels.cpp'
old = subprocess.check_output(['git', 'show', baseline+':'+name], cwd=args.source).decode('utf-8')
new = (args.source/name).read_text(encoding='utf-8')
def body(text):
    start = text.index('void CommonUpdateStateDataKernel::setPositions(')
    end = text.index('\nvoid CommonUpdateStateDataKernel::getVelocities(', start)
    return text[start:end]
pieces = {
    'exact_baseline.inc': body(old),
    'exact_candidate.inc': body(new),
    'exact_helper.inc': new[new.index('template <class CopyRange>'):new.index('void CommonUpdateStateDataKernel::setPositions(')],
    'exact_copy_coordinates.inc': (args.source/'platforms/common/src/kernels/copyCoordinateBuffers.cc').read_text(encoding='utf-8'),
}
for name, text in pieces.items():
    (args.output/name).write_text(text, encoding='utf-8')
