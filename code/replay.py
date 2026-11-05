#!/usr/bin/env python3
"""Replay the exact finite calculations in a new, separate output directory."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', required=True, type=Path,
                    help='new output directory; existing nonempty directories are rejected')
args = parser.parse_args()
source = Path(__file__).resolve().parent
out = args.output.resolve()
if out.exists() and any(out.iterdir()):
    parser.error('output directory is nonempty; choose a new directory')
out.mkdir(parents=True, exist_ok=True)
working = out / 'working'
working.mkdir()
names = ['mixed_potential_space.py', 'mixed_candidate.py',
         'mixed_explicit_relations.py', 'candidate_certify.py',
         'finite_multiplication.py', 'center_algebra.py',
         'actual_projection_check.py', 'transport_check.py',
         'verify_quadratic_table.py', 'check_constants.py']
hashes = {}
for name in names:
    src = source / name
    hashes[name] = hashlib.sha256(src.read_bytes()).hexdigest()
    shutil.copy2(src, working / name)
runs = []
for name in ['candidate_certify.py', 'mixed_explicit_relations.py',
             'finite_multiplication.py', 'verify_quadratic_table.py',
             'check_constants.py', 'center_algebra.py',
             'actual_projection_check.py', 'transport_check.py']:
    start = time.monotonic()
    print('Running', name, flush=True)
    cmd = [sys.executable, str(working / name)]
    if name == 'check_constants.py':
        cmd += ['--data', str(working/'finite_multiplication.json'),
                '--tex', str(source.parent/'source'/'constants.tex')]
    proc = subprocess.run(cmd, cwd=working, text=True, capture_output=True)
    stem = Path(name).stem
    (out / (stem + '.stdout.txt')).write_text(proc.stdout)
    (out / (stem + '.stderr.txt')).write_text(proc.stderr)
    runs.append({'script': name, 'command': cmd, 'returncode': proc.returncode,
                 'elapsed_seconds': round(time.monotonic()-start, 3)})
    if proc.returncode:
        (out / 'replay_summary.json').write_text(json.dumps(
            {'status': 'FAIL', 'runs': runs, 'source_sha256': hashes}, indent=2)+'\n')
        raise SystemExit(proc.returncode)
    if name in ('center_algebra.py', 'transport_check.py',
                'verify_quadratic_table.py', 'check_constants.py'):
        report = json.loads(proc.stdout)
        (out / (stem + '.json')).write_text(json.dumps(report, indent=2)+'\n')

table = json.loads((working/'finite_multiplication.json').read_text())
assert table['checks'] == {
    'basis_dimension': 83, 'multiplets': 55, 'nonzero_products': 480,
    'nonzero_channels': 270, 'associativity_composable_triples_degree_le_4': 4584,
    'equivariance_equations': 2834, 'channel_coefficient_equations': 1146}
transport = json.loads((out/'transport_check.json').read_text())
assert transport['tested_channels'] == 623
assert transport['tested_input_monomials'] == 3805 and not transport['failures']
centre = json.loads((out/'center_algebra.json').read_text())
assert centre['centre_dimension'] == 9 and centre['commutator_constraint_rank'] == 74
assert centre['basis_product_checks'] == 81 and centre['socle_dimension'] == 5
direct = json.loads((out/'verify_quadratic_table.json').read_text())
assert direct['ordered_products_checked'] == 6889 and direct['nonzero_products'] == 480
assert direct['quotient_dimensions'] == [5, 20, 33, 20, 5, 0]
constants = json.loads((out/'check_constants.json').read_text())
assert constants['nonunit_rows'] == 165
summary = {'status': 'PASS', 'python': sys.version, 'executable': sys.executable,
           'parallelism': 1, 'arithmetic': 'exact integers and fractions.Fraction',
           'source_sha256': hashes, 'runs': runs,
           'scope': 'Exact finite algebra and bounded transport checks; the VOA comparison and uniform-index bridge are the manuscript proofs.'}
(out/'replay_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print('PASS: results saved in', out)
