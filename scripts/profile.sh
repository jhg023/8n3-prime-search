#!/bin/bash
# Profile script for 8n3-prime-search
# Automates: debug build, Instruments profiling, and data extraction

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
TRACE_OUTPUT="${PROJECT_DIR}/profile.trace"
TEMPLATE="My PMU Counters"

# Default search range (short for profiling)
START="${1:-1e12}"
END="${2:-1.00001e12}"

echo "=========================================="
echo "Prime Search Profiling Script"
echo "=========================================="
echo "Project: $PROJECT_DIR"
echo "Range: $START to $END"
echo ""

# Step 1: Build with debug symbols + optimizations
echo "[1/4] Building optimized binary with debug symbols..."
cd "$PROJECT_DIR"
make clean >/dev/null 2>&1 || true
gcc -Wall -Wextra -std=c11 -Iinclude \
    -O3 -march=native -mtune=native -flto -funroll-loops \
    -g \
    -o search src/search.c -lm
echo "      Build complete: $(file search | cut -d: -f2)"

# Step 2: Run profiling with Instruments
echo ""
echo "[2/4] Profiling with '$TEMPLATE' template..."
rm -rf "$TRACE_OUTPUT" 2>/dev/null || true

if ! xctrace record --template "$TEMPLATE" --output "$TRACE_OUTPUT" \
    --launch -- ./search "$START" "$END" 2>&1; then
    echo "ERROR: Profiling failed. Trying alternative template..."
    TEMPLATE="My CPU Counters"
    xctrace record --template "$TEMPLATE" --output "$TRACE_OUTPUT" \
        --launch -- ./search "$START" "$END" 2>&1 || {
        echo "ERROR: Both templates failed. Check Instruments configuration."
        exit 1
    }
fi
echo "      Trace saved: $TRACE_OUTPUT"

# Step 3: Symbolicate
echo ""
echo "[3/4] Symbolicating trace..."
xctrace symbolicate --input "$TRACE_OUTPUT" --dsym ./search 2>&1 | grep -v '^\[' || true
echo "      Symbolication complete"

# Step 4: Extract and analyze data
echo ""
echo "[4/4] Extracting profiling data..."
echo ""

# Export counters-profile table
EXPORT_FILE="/tmp/profile_export_$$.xml"
xctrace export --input "$TRACE_OUTPUT" \
    --xpath '/trace-toc/run/data/table[@schema="counters-profile"]' \
    > "$EXPORT_FILE" 2>/dev/null

# Extract function hotspots
echo "=========================================="
echo "FUNCTION HOTSPOTS (by sample count)"
echo "=========================================="
grep -o 'name="[^"]*"' "$EXPORT_FILE" 2>/dev/null | \
    sed 's/name="//g; s/"//g' | \
    grep -v '^0x' | \
    grep -v '^$' | \
    sort | uniq -c | sort -rn | head -20 || echo "(no data)"

# Export PMC data for detailed analysis
PMC_FILE="/tmp/pmc_export_$$.xml"
xctrace export --input "$TRACE_OUTPUT" \
    --xpath '/trace-toc/run/data/table[@schema="kdebug-counters-with-time-sample"]' \
    > "$PMC_FILE" 2>/dev/null

# Analyze PMC counters with Python
echo ""
echo "=========================================="
echo "PMU COUNTER ANALYSIS"
echo "=========================================="
echo "Counters: Cycles, Instructions, INST_BRANCH, BRANCH_MISPRED_NONSPEC, L1D_CACHE_MISS_LD"
echo ""

python3 - "$PMC_FILE" << 'PYTHON_SCRIPT'
import re
import sys

pmc_file = sys.argv[1]

try:
    with open(pmc_file, 'r') as f:
        content = f.read()
except Exception as e:
    print(f"Could not read PMC export file: {e}")
    sys.exit(1)

# Extract pmc-events with unique IDs (not refs)
# Format: pmc-events id="N" fmt="(cycles), (instrs), (branches), (mispred), (l1d_miss)"
id_pattern = r'pmc-events id="(\d+)" fmt="([^"]+)"'
id_matches = re.findall(id_pattern, content)

if not id_matches:
    print("No PMC data found in trace")
    sys.exit(0)

def parse_pmc_fmt(fmt_str):
    """Parse PMC format string into list of 5 counter values."""
    groups = re.findall(r'\(([^)]+)\)', fmt_str)
    if len(groups) < 5:
        return None
    try:
        return [int(g.replace(',', '')) for g in groups[:5]]
    except:
        return None

# Parse samples and sort by ID
samples = []
for id_val, fmt in id_matches:
    vals = parse_pmc_fmt(fmt)
    if vals:
        samples.append((int(id_val), vals))
samples.sort(key=lambda x: x[0])

if len(samples) < 2:
    print("Not enough samples for analysis")
    sys.exit(0)

# Compute deltas, filtering out core migrations
# Core migrations cause: counter decreases OR very large cycle jumps
MAX_CYCLE_DELTA = 20_000_000  # 20M cycles max per ~1ms sample
deltas = []
core_switches = 0

for i in range(1, len(samples)):
    prev = samples[i-1][1]
    curr = samples[i][1]
    delta = [curr[j] - prev[j] for j in range(5)]

    # Skip if any counter decreased (core switch)
    if any(d < 0 for d in delta):
        core_switches += 1
        continue

    # Skip if cycle delta is too large (core switch or scheduling gap)
    if delta[0] > MAX_CYCLE_DELTA or delta[0] < 1000:
        core_switches += 1
        continue

    deltas.append(delta)

if not deltas:
    print("Could not compute valid deltas")
    sys.exit(0)

# Calculate totals
total_cycles = sum(d[0] for d in deltas)
total_instrs = sum(d[1] for d in deltas)
total_branches = sum(d[2] for d in deltas)
total_mispred = sum(d[3] for d in deltas)
total_l1d = sum(d[4] for d in deltas)

print(f"Total samples:      {len(samples):>12,}")
print(f"Valid samples:      {len(deltas):>12,}")
print(f"Core switches:      {core_switches:>12,}")
print(f"Total Cycles:       {total_cycles:>12,}")
print(f"Total Instructions: {total_instrs:>12,}")
if total_cycles > 0:
    print(f"IPC:                {total_instrs/total_cycles:>12.3f}")
print(f"Branches:           {total_branches:>12,}")
print(f"Branch Mispreds:    {total_mispred:>12,}")
if total_branches > 0:
    rate = total_mispred / total_branches * 100
    print(f"Mispredict Rate:    {rate:>11.2f}%")
print(f"L1D Cache Misses:   {total_l1d:>12,}")
if total_instrs > 0:
    print(f"L1D Miss/1K Instr:  {total_l1d/total_instrs*1000:>12.3f}")

# Per-sample statistics
if len(deltas) > 10:
    ipcs = [d[1]/d[0] for d in deltas if d[0] > 0]
    mispreds = [d[3]/d[2]*100 for d in deltas if d[2] > 0]

    print("")
    print("Per-Sample Statistics:")
    print(f"  Avg IPC:          {sum(ipcs)/len(ipcs):>12.3f}")
    print(f"  Min IPC:          {min(ipcs):>12.3f}")
    print(f"  Max IPC:          {max(ipcs):>12.3f}")
    if mispreds:
        print(f"  Avg Mispredict:   {sum(mispreds)/len(mispreds):>11.2f}%")
        print(f"  Max Mispredict:   {max(mispreds):>11.2f}%")

# Show a few sample deltas for verification
print("")
print("Sample Deltas (first 5):")
print("  Cycles       Instructions   Branches    Mispred    L1D Miss")
for d in deltas[:5]:
    print(f"  {d[0]:>10,}   {d[1]:>12,}   {d[2]:>8,}   {d[3]:>8,}   {d[4]:>8,}")
PYTHON_SCRIPT

# Cleanup
rm -f "$EXPORT_FILE" "$PMC_FILE" 2>/dev/null || true

echo ""
echo "=========================================="
echo "Profiling complete!"
echo "=========================================="
echo ""
echo "To view in Instruments GUI:"
echo "  open $TRACE_OUTPUT"
echo ""
