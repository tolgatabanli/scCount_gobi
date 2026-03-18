#!/usr/bin/env bash
set -euo pipefail

# Runs read_tracker_v2.py for quality filters 60,50,...,10 with up to 6 parallel jobs.
# Outputs:
#   1) Combined final statistics from all runs
#   2) Time and memory usage report per run

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

PYTHON_SCRIPT="${SCRIPT_DIR}/read_tracker_v2.py"

BAM_A="../data/data/minimap/full_genome_short_rna_mode_Run/sorted.bam"
BAM_B="../data/data/star_run/Aligned.sortedByCoord.out.bam"
GTF="../data/data/Sus_scrofa.Sscrofa11.1.gtf"
REPEATS="../data/data/repeats/susScr3.fa.out"
PARALOG_TABLE="../data/paralogs.tsv"
SEED_FILE="../data/genes_minimap_spaceranger_70p_discrepancy.csv"

MAX_PARALLEL=6
QUALITIES=(60 50 40 30 20 10)

OUT_DIR="${SCRIPT_DIR}/quality_sweep_$(date +%Y%m%d_%H%M%S)"
RUN_LOG_DIR="${OUT_DIR}/logs"
STAT_SNIPPET_DIR="${OUT_DIR}/stat_snippets"
TIME_SNIPPET_DIR="${OUT_DIR}/time_snippets"
mkdir -p "${RUN_LOG_DIR}" "${STAT_SNIPPET_DIR}" "${TIME_SNIPPET_DIR}"

STATS_FILE="${OUT_DIR}/quality_filter_stats.txt"
TIME_MEM_FILE="${OUT_DIR}/quality_filter_time_memory.txt"

if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
  echo "ERROR: Python script not found: ${PYTHON_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${SCRIPT_DIR}/${SEED_FILE}" ]]; then
  echo "ERROR: Seed file not found: ${SCRIPT_DIR}/${SEED_FILE}" >&2
  exit 1
fi

SEED_IDS="$(tr -d '\n' < "${SCRIPT_DIR}/${SEED_FILE}")"

extract_final_stats() {
  local input_log="$1"
  local output_stats="$2"

  awk '
    /^total visited loci:/ ||
    /^tracked reads:/ ||
    /^mapped else where in a:/ ||
    /^mapped else where in b:/ ||
    /^mapped else where / ||
    /^not present in other bam / ||
    /^filtered by quality \(total\) / ||
    /^  filtered in bam A / ||
    /^  filtered in bam B / ||
    /^paralogous / ||
    /^repeat / ||
    /^warning: [0-9]+ tracked reads were not assigned a terminal status/
  ' "${input_log}" > "${output_stats}"
}

run_one_quality() {
  local q="$1"
  local run_log="${RUN_LOG_DIR}/q${q}.log"
  local run_time="${TIME_SNIPPET_DIR}/q${q}.txt"
  local run_stats="${STAT_SNIPPET_DIR}/q${q}.txt"

  (
    cd "${SCRIPT_DIR}"
    /usr/bin/time -f "elapsed_seconds=%e\tmax_rss_kb=%M\texit_code=%x" -o "${run_time}" \
      python "${PYTHON_SCRIPT}" \
      -a "${BAM_A}" \
      -b "${BAM_B}" \
      -g "${GTF}" \
      -r "${REPEATS}" \
      --paralog_table "${PARALOG_TABLE}" \
      -s "${SEED_IDS}" \
      --quality_filter "${q}" \
      > "${run_log}" 2>&1
  )

  extract_final_stats "${run_log}" "${run_stats}"
}

running_jobs=0
for q in "${QUALITIES[@]}"; do
  run_one_quality "${q}" &
  ((running_jobs+=1))

  if (( running_jobs >= MAX_PARALLEL )); then
    wait -n
    ((running_jobs-=1))
  fi
done

wait

{
  echo "quality_filter\ttotal visited loci\ttracked reads\tmapped else where in a\tmapped else where in b\tmapped else where\tnot present in other bam\tfiltered by quality (total)\tfiltered in bam A\tfiltered in bam B\tparalogous\trepeat"
  for q in "${QUALITIES[@]}"; do
    f="${STAT_SNIPPET_DIR}/q${q}.txt"

    total_visited="$(grep -E '^total visited loci:' "${f}" | sed 's/^total visited loci: //')"
    tracked="$(grep -E '^tracked reads:' "${f}" | sed 's/^tracked reads: //')"
    in_a="$(grep -E '^mapped else where in a:' "${f}" | sed 's/^mapped else where in a: //')"
    in_b="$(grep -E '^mapped else where in b:' "${f}" | sed 's/^mapped else where in b: //')"
    mapped="$(grep -E '^mapped else where ' "${f}" | head -n1 | sed 's/^mapped else where //')"
    missing="$(grep -E '^not present in other bam ' "${f}" | sed 's/^not present in other bam //')"
    filtered_total="$(grep -E '^filtered by quality \(total\) ' "${f}" | sed 's/^filtered by quality (total) //')"
    filtered_a="$(grep -E '^  filtered in bam A ' "${f}" | sed 's/^  filtered in bam A //')"
    filtered_b="$(grep -E '^  filtered in bam B ' "${f}" | sed 's/^  filtered in bam B //')"
    para="$(grep -E '^paralogous ' "${f}" | sed 's/^paralogous //')"
    rep="$(grep -E '^repeat ' "${f}" | sed 's/^repeat //')"

    echo -e "${q}\t${total_visited}\t${tracked}\t${in_a}\t${in_b}\t${mapped}\t${missing}\t${filtered_total}\t${filtered_a}\t${filtered_b}\t${para}\t${rep}"
  done
} > "${STATS_FILE}"

{
  echo -e "quality_filter\telapsed_seconds\tmax_rss_kb\texit_code"
  for q in "${QUALITIES[@]}"; do
    sed -E "s/^elapsed_seconds=([0-9.]+)\tmax_rss_kb=([0-9]+)\texit_code=([0-9-]+)/${q}\t\1\t\2\t\3/" "${TIME_SNIPPET_DIR}/q${q}.txt"
  done
} > "${TIME_MEM_FILE}"

echo "Completed quality sweep."
echo "Stats file: ${STATS_FILE}"
echo "Time/memory file: ${TIME_MEM_FILE}"
