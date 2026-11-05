#!/usr/bin/env python3
"""Compare the companion's nonunit table with the generated exact JSON data."""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import re

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data', type=Path, default=root / 'code' / 'finite_multiplication.json')
parser.add_argument('--tex', type=Path, default=root / 'source' / 'constants.tex')
args = parser.parse_args()
tex = args.tex.resolve()
data_path = args.data.resolve()
text = tex.read_text()
assert text.count(r'\begin{longtable}') == text.count(r'\end{longtable}') == 1
table = text.split(r'\begin{longtable}', 1)[1].split(r'\end{longtable}', 1)[0]
data = json.loads(data_path.read_text())
multiplets = data['multiplets']
ids = {m['name']: i for i, m in enumerate(multiplets)}
rows = {}
for line in table.splitlines():
    if not line.startswith(r'\texttt{'):
        continue
    match = re.fullmatch(r'\\texttt\{([^}]+)\} & \\texttt\{([^}]+)\} & \\texttt\{([^}]+)\} & \$([^$]+)\$\\\\', line)
    assert match, line
    key = tuple(ids[n.replace(r'\_', '_')] for n in match.groups()[:3])
    value = match.group(4)
    fraction = re.fullmatch(r'\\frac\{(-?\d+)\}\{(\d+)\}', value)
    value = Fraction(int(fraction.group(1)), int(fraction.group(2))) if fraction else Fraction(value)
    assert key not in rows
    rows[key] = value
expected = {(c['first'], c['second'], c['output']): Fraction(c['coefficient'])
            for c in data['equivariant_constants']
            if multiplets[c['first']]['degree'] and multiplets[c['second']]['degree']}
assert rows == expected and len(rows) == 165
print(json.dumps({'status': 'PASS', 'nonunit_rows': len(rows),
                  'unit_channels': len(data['equivariant_constants']) - len(expected),
                  'arithmetic': 'fractions.Fraction',
                  'source_sha256': hashlib.sha256(tex.read_bytes()).hexdigest(),
                  'table_sha256': hashlib.sha256(data_path.read_bytes()).hexdigest()}, indent=2))
