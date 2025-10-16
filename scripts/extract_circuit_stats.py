#!/usr/bin/env python3
"""Extract circuit statistics from simulator logs."""
from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

ANSI_ESCAPE_RE = re.compile(r"\x1B\[[0-?]*[ -/]*[@-~]")
LOG_DIM_RE = re.compile(r"log_dim\s*=\s*(\d+)")
NON_FREE_DEPTH_RE = re.compile(r"poly circuit non_free_depth\s+(\d+)")
LIMB_RE = re.compile(r"limb_?(\d+)")


@dataclass
class CircuitBlock:
    lines: Tuple[str, str, str]
    log_dim: int
    gate_counts: Dict[str, int]
    non_free_depth: int


def strip_ansi(text: str) -> str:
    return ANSI_ESCAPE_RE.sub('', text)


def parse_gate_counts(line: str) -> Dict[str, int]:
    start = line.find('{')
    end = line.rfind('}')
    if start == -1 or end == -1 or end <= start:
        raise ValueError(f'circuit size line missing gate list: {line!r}')
    gate_section = line[start + 1:end]
    counts: Dict[str, int] = {}
    for part in gate_section.split(','):
        part = part.strip()
        if not part:
            continue
        if ':' not in part:
            raise ValueError(f'malformed gate entry: {part!r}')
        name, value = part.split(':', 1)
        counts[name.strip()] = int(value.strip())
    return counts


def parse_blocks(path: Path) -> List[CircuitBlock]:
    lines = [strip_ansi(line.rstrip('\n')) for line in path.read_text().splitlines()]
    blocks: List[CircuitBlock] = []
    idx = 0
    while idx + 2 < len(lines):
        first, second, third = lines[idx], lines[idx + 1], lines[idx + 2]
        if 'circuit constructed with' in first and 'circuit size' in second and 'poly circuit non_free_depth' in third:
            log_dim_match = LOG_DIM_RE.search(first)
            depth_match = NON_FREE_DEPTH_RE.search(third)
            if log_dim_match and depth_match:
                log_dim = int(log_dim_match.group(1))
                gate_counts = parse_gate_counts(second)
                non_free_depth = int(depth_match.group(1))
                blocks.append(CircuitBlock((first, second, third), log_dim, gate_counts, non_free_depth))
                idx += 3
                continue
        idx += 1
    return blocks


def list_blocks(blocks: Iterable[CircuitBlock]) -> None:
    for number, block in enumerate(blocks, start=1):
        print(f'Block {number}:')
        for line in block.lines:
            print(line)
        print()


GATE_ORDER = ['Input', 'Add', 'Sub', 'SmallScalarMul', 'LargeScalarMul', 'Mul']
PUB_LUT_PREFIX = 'PubLut'
PUB_LUT_RE = re.compile(r"PubLut\((\d+)\)")


def ordered_gate_names(names: Iterable[str]) -> List[str]:
    name_set = set(names)
    ordered: List[str] = [name for name in GATE_ORDER if name in name_set]
    pub_luts = []
    for name in name_set:
        if name.startswith(PUB_LUT_PREFIX):
            match = PUB_LUT_RE.search(name)
            index = int(match.group(1)) if match else -1
            pub_luts.append((index, name))
    ordered.extend(name for _, name in sorted(pub_luts))
    remaining = sorted(name_set - set(ordered))
    ordered.extend(remaining)
    return ordered


def write_csv(rows: List[Tuple[str, int | None, CircuitBlock]], output: Path | None) -> None:
    gate_names = ordered_gate_names(gate for _, _, block in rows for gate in block.gate_counts)
    header = ['prefix', 'limb_bit_size', 'non_free_depth', *gate_names]
    stream = output.open('w', newline='') if output else sys.stdout
    should_close = output is not None
    try:
        writer = csv.writer(stream)
        writer.writerow(header)
        for prefix, limb_bit_size, block in rows:
            limb_value = limb_bit_size if limb_bit_size is not None else ''
            row = [prefix, limb_value, block.non_free_depth]
            row.extend(block.gate_counts.get(gate, 0) for gate in gate_names)
            writer.writerow(row)
    finally:
        if should_close:
            stream.close()


def prefix_from_name(path: Path) -> str:
    name = path.name
    marker = '.params'
    return name.split(marker, 1)[0] if marker in name else name.rsplit('.', 1)[0]


def limb_bit_size_from_prefix(prefix: str) -> int | None:
    match = LIMB_RE.search(prefix)
    if match:
        return int(match.group(1))
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description='Collect circuit statistics from simulator logs.')
    parser.add_argument('--log', type=Path, help='Path to a single log file to inspect.')
    parser.add_argument('--logs-dir', type=Path, default=Path('abe/logs'), help='Directory containing log files (default: abe/logs).')
    parser.add_argument('--output', type=Path, help='Optional CSV output path for directory mode.')
    args = parser.parse_args()

    if args.log:
        if not args.log.is_file():
            parser.error(f'log file not found: {args.log}')
        blocks = parse_blocks(args.log)
        if not blocks:
            print('No matching circuit blocks found.', file=sys.stderr)
            sys.exit(1)
        list_blocks(blocks)
        best = min(blocks, key=lambda block: block.log_dim)
        prefix = prefix_from_name(args.log)
        limb_bit_size = limb_bit_size_from_prefix(prefix)
        print(f'Selected block with minimal log_dim = {best.log_dim}')
        if limb_bit_size is not None:
            print(f'limb_bit_size: {limb_bit_size}')
        print(f'non_free_depth: {best.non_free_depth}')
        for gate in ordered_gate_names(best.gate_counts):
            print(f'{gate}: {best.gate_counts[gate]}')
        return

    logs_dir = args.logs_dir
    if not logs_dir.is_dir():
        parser.error(f'logs directory not found: {logs_dir}')

    rows: List[Tuple[str, int | None, CircuitBlock]] = []
    for path in sorted(logs_dir.glob('*.log')):
        blocks = parse_blocks(path)
        if not blocks:
            continue
        best = min(blocks, key=lambda block: block.log_dim)
        prefix = prefix_from_name(path)
        limb_bit_size = limb_bit_size_from_prefix(prefix)
        rows.append((prefix, limb_bit_size, best))

    if not rows:
        print('No circuit statistics found in logs.', file=sys.stderr)
        sys.exit(1)

    write_csv(rows, args.output)


if __name__ == '__main__':
    main()
