#!/usr/bin/env bash

# Re-download and uniformly process the selected public N. crassa H3K4
# ChIP-seq runs for G3 revision issue #69.
#
# Selected datasets (approved 1/2/4/5/6 plan):
#   Ferraro et al. 2021: H3K4me1/2/3 (SRR12229306-8)
#   Sasaki et al. 2014: H3K4me2 (SRR1295547)
#   Storck et al. 2020: H3K4me3 (SRR12202381)
#
# The script deliberately does not:
#   - use the unmatched Ferraro hH3-3xFLAG input;
#   - include Zhu et al. 2019 light/dark samples;
#   - include the Storck delta-lsd1 sample;
#   - pool reads across studies;
#   - call peaks without a matched input; or
#   - reuse the historical RPKM tracks.

set -Eeuo pipefail
IFS=$'\n\t'

readonly DEFAULT_WORK_ROOT="/Volumes/Garage/Re_analysis/260906_issue69_H3K4"
readonly REFERENCE_ACCESSION="GCA_000182925.2"
readonly REFERENCE_BASENAME="GCA_000182925.2_NC12"
readonly REFERENCE_URL_ROOT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/925/GCA_000182925.2_NC12"
readonly SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORK_ROOT="${WORK_ROOT:-${DEFAULT_WORK_ROOT}}"
THREADS="${THREADS:-8}"
MIN_MAPQ="${MIN_MAPQ:-20}"
BIN_SIZE="${BIN_SIZE:-10}"
PYTHON_VERSION="${PYTHON_VERSION:-3.12.13}"
STAGE="${1:-all}"

SRA_ROOT="${WORK_ROOT}/sra"
RAW_FASTQ_ROOT="${WORK_ROOT}/fastq/raw"
TRIMMED_FASTQ_ROOT="${WORK_ROOT}/fastq/trimmed"
REFERENCE_ROOT="${WORK_ROOT}/reference"
INDEX_ROOT="${WORK_ROOT}/reference/bowtie2_index"
BAM_ROOT="${WORK_ROOT}/bam"
BIGWIG_ROOT="${WORK_ROOT}/bigwig"
QC_ROOT="${WORK_ROOT}/qc"
LOG_ROOT="${WORK_ROOT}/logs"
METADATA_ROOT="${WORK_ROOT}/metadata"
TMP_ROOT="${WORK_ROOT}/tmp"

REFERENCE_FASTA_GZ="${REFERENCE_ROOT}/${REFERENCE_BASENAME}_genomic.fna.gz"
REFERENCE_GFF_GZ="${REFERENCE_ROOT}/${REFERENCE_BASENAME}_genomic.gff.gz"
REFERENCE_FASTA="${REFERENCE_ROOT}/${REFERENCE_BASENAME}_genomic.fna"
REFERENCE_GFF="${REFERENCE_ROOT}/${REFERENCE_BASENAME}_genomic.gff"
INDEX_PREFIX="${INDEX_ROOT}/${REFERENCE_BASENAME}"

# run|study|mark|sample label
readonly RUN_MANIFEST=$(cat <<'EOF'
SRR12229306|Ferraro2021|H3K4me1|Ferraro2021_WT_H3K4me1
SRR12229307|Ferraro2021|H3K4me2|Ferraro2021_WT_H3K4me2
SRR12229308|Ferraro2021|H3K4me3|Ferraro2021_WT_H3K4me3
SRR1295547|Sasaki2014|H3K4me2|Sasaki2014_WT_H3K4me2
SRR12202381|Storck2020|H3K4me3|Storck2020_WT_H3K4me3
EOF
)

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

prepare_directories() {
    mkdir -p \
        "${SRA_ROOT}" \
        "${RAW_FASTQ_ROOT}" \
        "${TRIMMED_FASTQ_ROOT}" \
        "${REFERENCE_ROOT}" \
        "${INDEX_ROOT}" \
        "${BAM_ROOT}" \
        "${BIGWIG_ROOT}" \
        "${QC_ROOT}/fastqc_raw" \
        "${QC_ROOT}/fastqc_trimmed" \
        "${QC_ROOT}/alignment" \
        "${QC_ROOT}/duplicates" \
        "${LOG_ROOT}" \
        "${METADATA_ROOT}" \
        "${TMP_ROOT}/fasterq" \
        "${TMP_ROOT}/matplotlib"
}

resolve_tools() {
    local command_name

    for command_name in \
        prefetch fasterq-dump vdb-validate fastqc trimmomatic bowtie2 Rscript \
        bowtie2-build samtools picard curl gzip shasum; do
        require_command "${command_name}"
    done

    require_command pyenv
    export PYENV_VERSION="${PYTHON_VERSION}"
    export MPLCONFIGDIR="${TMP_ROOT}/matplotlib"

    BAM_COVERAGE=$(pyenv which bamCoverage) || die "bamCoverage is unavailable in pyenv Python ${PYTHON_VERSION}"
    [[ -x "${BAM_COVERAGE}" ]] || die "bamCoverage is not executable: ${BAM_COVERAGE}"
    "${BAM_COVERAGE}" --version >/dev/null
    Rscript -e '
        required <- c("IRanges", "rtracklayer", "digest", "jsonlite")
        missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
        if (length(missing) > 0L) {
            stop("Missing required R packages: ", paste(missing, collapse = ", "))
        }
    ' >/dev/null || die "Required R packages are unavailable"

    if command -v brew >/dev/null 2>&1; then
        TRIMMOMATIC_ADAPTERS="${TRIMMOMATIC_ADAPTERS:-$(brew --prefix trimmomatic)/share/trimmomatic/adapters/TruSeq3-SE.fa}"
    else
        TRIMMOMATIC_ADAPTERS="${TRIMMOMATIC_ADAPTERS:-}"
    fi
    [[ -f "${TRIMMOMATIC_ADAPTERS}" ]] || die "Trimmomatic adapter FASTA not found; set TRIMMOMATIC_ADAPTERS"
}

record_manifest() {
    local manifest_path="${METADATA_ROOT}/selected_runs.tsv"

    {
        printf 'run\tstudy\tmark\tsample_label\n'
        printf '%s\n' "${RUN_MANIFEST}" | tr '|' '\t'
    } > "${manifest_path}"
}

record_versions() {
    local versions_path="${METADATA_ROOT}/software_versions.txt"

    {
        printf 'generated_at\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        printf 'reference_accession\t%s\n' "${REFERENCE_ACCESSION}"
        printf 'python\t%s\n' "$(pyenv which python)"
        "$(pyenv which python)" --version 2>&1
        prefetch --version 2>&1 | sed -n '1p' || true
        fasterq-dump --version 2>&1 | sed -n '1p' || true
        fastqc --version 2>&1 | sed -n '1p' || true
        trimmomatic -version 2>&1 | sed -n '1p' || true
        bowtie2 --version 2>&1 | sed -n '1p' || true
        samtools --version 2>&1 | sed -n '1,2p' || true
        picard MarkDuplicates --version 2>&1 | sed -n '1p' || true
        "${BAM_COVERAGE}" --version 2>&1 | sed -n '1p' || true
        Rscript --version 2>&1 | sed -n '1p' || true
    } > "${versions_path}"
}

check_environment() {
    [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
    [[ "${MIN_MAPQ}" =~ ^[0-9]+$ ]] || die "MIN_MAPQ must be a non-negative integer"
    [[ "${BIN_SIZE}" =~ ^[1-9][0-9]*$ ]] || die "BIN_SIZE must be a positive integer"

    prepare_directories
    resolve_tools
    record_manifest
    record_versions

    log "Work root: ${WORK_ROOT}"
    log "Free space at work root:"
    df -h "${WORK_ROOT}"
    log "pyenv bamCoverage: ${BAM_COVERAGE}"
    log "Environment check completed"
}

download_file() {
    local url="$1"
    local destination="$2"
    local partial="${destination}.part"

    if [[ -s "${destination}" ]]; then
        log "Keeping existing file: ${destination}"
        return
    fi

    log "Downloading ${url}"
    curl --fail --location --retry 5 --retry-delay 5 \
        --output "${partial}" "${url}"
    mv "${partial}" "${destination}"
}

prepare_reference() {
    download_file \
        "${REFERENCE_URL_ROOT}/${REFERENCE_BASENAME}_genomic.fna.gz" \
        "${REFERENCE_FASTA_GZ}"
    download_file \
        "${REFERENCE_URL_ROOT}/${REFERENCE_BASENAME}_genomic.gff.gz" \
        "${REFERENCE_GFF_GZ}"

    if [[ ! -s "${REFERENCE_FASTA}" ]]; then
        gzip -dc "${REFERENCE_FASTA_GZ}" > "${REFERENCE_FASTA}.part"
        mv "${REFERENCE_FASTA}.part" "${REFERENCE_FASTA}"
    fi
    if [[ ! -s "${REFERENCE_GFF}" ]]; then
        gzip -dc "${REFERENCE_GFF_GZ}" > "${REFERENCE_GFF}.part"
        mv "${REFERENCE_GFF}.part" "${REFERENCE_GFF}"
    fi

    shasum -a 256 \
        "${REFERENCE_FASTA_GZ}" "${REFERENCE_GFF_GZ}" \
        "${REFERENCE_FASTA}" "${REFERENCE_GFF}" \
        > "${METADATA_ROOT}/reference_sha256.txt"
    samtools faidx "${REFERENCE_FASTA}"

    if [[ ! -s "${INDEX_PREFIX}.1.bt2" && ! -s "${INDEX_PREFIX}.1.bt2l" ]]; then
        log "Building Bowtie2 index for ${REFERENCE_ACCESSION}"
        bowtie2-build --threads "${THREADS}" "${REFERENCE_FASTA}" "${INDEX_PREFIX}" \
            > "${LOG_ROOT}/bowtie2-build.stdout.log" \
            2> "${LOG_ROOT}/bowtie2-build.stderr.log"
    fi

    log "Reference preparation completed"
}

compress_fastq() {
    local fastq_path="$1"

    if command -v bgzip >/dev/null 2>&1; then
        bgzip --threads "${THREADS}" "${fastq_path}"
    else
        gzip -n "${fastq_path}"
    fi
}

normalize_single_end_fastq() {
    local run="$1"
    local fastq_path="${RAW_FASTQ_ROOT}/${run}.fastq"
    local split_fastq_1="${RAW_FASTQ_ROOT}/${run}_1.fastq"
    local split_fastq_2="${RAW_FASTQ_ROOT}/${run}_2.fastq"

    if [[ -s "${fastq_path}" ]]; then
        [[ ! -e "${split_fastq_1}" && ! -e "${split_fastq_2}" ]] || \
            die "Ambiguous FASTQ outputs for ${run}: canonical and split files coexist"
        return
    fi

    # Some nominally single-end SRA runs contain one biological read followed
    # by a zero-length technical read. With --split-files, fasterq-dump writes
    # the biological reads to *_1.fastq rather than to the canonical filename.
    if [[ -s "${split_fastq_1}" && ! -e "${split_fastq_2}" ]]; then
        log "Normalizing single biological read file for ${run}"
        mv "${split_fastq_1}" "${fastq_path}"
        return
    fi

    if [[ -e "${split_fastq_1}" || -e "${split_fastq_2}" ]]; then
        die "Unexpected paired or incomplete FASTQ output for ${run}; inspect ${RAW_FASTQ_ROOT}"
    fi
}

download_one_run() {
    local run="$1"
    local sra_path="${SRA_ROOT}/${run}/${run}.sra"
    local fastq_path="${RAW_FASTQ_ROOT}/${run}.fastq"
    local fastq_gz_path="${fastq_path}.gz"
    local run_tmp="${TMP_ROOT}/fasterq/${run}"

    if [[ ! -s "${sra_path}" ]]; then
        log "Prefetching ${run}"
        prefetch "${run}" --max-size 100G --output-directory "${SRA_ROOT}" \
            > "${LOG_ROOT}/${run}.prefetch.stdout.log" \
            2> "${LOG_ROOT}/${run}.prefetch.stderr.log"
    fi
    vdb-validate "${sra_path}" \
        > "${LOG_ROOT}/${run}.vdb-validate.stdout.log" \
        2> "${LOG_ROOT}/${run}.vdb-validate.stderr.log"

    if [[ ! -s "${fastq_gz_path}" ]]; then
        mkdir -p "${run_tmp}"
        normalize_single_end_fastq "${run}"
        if [[ ! -s "${fastq_path}" ]]; then
            log "Converting ${run} to FASTQ"
            fasterq-dump "${sra_path}" \
                --threads "${THREADS}" \
                --temp "${run_tmp}" \
                --outdir "${RAW_FASTQ_ROOT}" \
                --split-files \
                > "${LOG_ROOT}/${run}.fasterq-dump.stdout.log" \
                2> "${LOG_ROOT}/${run}.fasterq-dump.stderr.log"
            normalize_single_end_fastq "${run}"
        fi
        [[ -s "${fastq_path}" ]] || die "Expected single-end FASTQ was not created: ${fastq_path}"
        compress_fastq "${fastq_path}"
    fi

    shasum -a 256 "${sra_path}" "${fastq_gz_path}" \
        > "${METADATA_ROOT}/${run}.download_sha256.txt"
    fastqc --threads "${THREADS}" --outdir "${QC_ROOT}/fastqc_raw" "${fastq_gz_path}" \
        > "${LOG_ROOT}/${run}.fastqc_raw.stdout.log" \
        2> "${LOG_ROOT}/${run}.fastqc_raw.stderr.log"
}

download_runs() {
    local run study mark label

    while IFS='|' read -r run study mark label; do
        [[ -n "${run}" ]] || continue
        download_one_run "${run}"
    done <<< "${RUN_MANIFEST}"

    log "SRA and FASTQ download completed"
}

trim_one_run() {
    local run="$1"
    local raw_fastq="${RAW_FASTQ_ROOT}/${run}.fastq.gz"
    local trimmed_fastq="${TRIMMED_FASTQ_ROOT}/${run}.trimmed.fastq.gz"

    [[ -s "${raw_fastq}" ]] || die "Raw FASTQ not found: ${raw_fastq}"

    if [[ ! -s "${trimmed_fastq}" ]]; then
        log "Trimming ${run}"
        trimmomatic SE \
            -threads "${THREADS}" \
            -phred33 \
            "${raw_fastq}" \
            "${trimmed_fastq}" \
            "ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:2:30:10" \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
            > "${LOG_ROOT}/${run}.trimmomatic.stdout.log" \
            2> "${LOG_ROOT}/${run}.trimmomatic.stderr.log"
    fi

    shasum -a 256 "${trimmed_fastq}" \
        > "${METADATA_ROOT}/${run}.trimmed_fastq_sha256.txt"
    fastqc --threads "${THREADS}" --outdir "${QC_ROOT}/fastqc_trimmed" "${trimmed_fastq}" \
        > "${LOG_ROOT}/${run}.fastqc_trimmed.stdout.log" \
        2> "${LOG_ROOT}/${run}.fastqc_trimmed.stderr.log"
}

trim_runs() {
    local run study mark label

    while IFS='|' read -r run study mark label; do
        [[ -n "${run}" ]] || continue
        trim_one_run "${run}"
    done <<< "${RUN_MANIFEST}"

    log "Read trimming completed"
}

align_one_run() {
    local run="$1"
    local label="$2"
    local trimmed_fastq="${TRIMMED_FASTQ_ROOT}/${run}.trimmed.fastq.gz"
    local sorted_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.sorted.bam"
    local marked_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.marked.bam"
    local deduplicated_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.nonduplicate.bam"
    local duplicate_metrics="${QC_ROOT}/duplicates/${label}.picard_markduplicates.txt"

    [[ -s "${trimmed_fastq}" ]] || die "Trimmed FASTQ not found: ${trimmed_fastq}"

    if [[ ! -s "${sorted_bam}" ]]; then
        log "Aligning ${run} (${label})"
        bowtie2 \
            --very-sensitive \
            -p "${THREADS}" \
            -x "${INDEX_PREFIX}" \
            --un-gz "${QC_ROOT}/alignment/${label}.unaligned.fastq.gz" \
            -U "${trimmed_fastq}" \
            2> "${LOG_ROOT}/${label}.bowtie2.stderr.log" \
        | samtools view \
            -@ "${THREADS}" \
            -b \
            -q "${MIN_MAPQ}" \
            -F 4 \
            - \
        | samtools sort \
            -@ "${THREADS}" \
            -O BAM \
            -o "${sorted_bam}" \
            -
        samtools index --threads "${THREADS}" "${sorted_bam}"
    fi

    if [[ ! -s "${marked_bam}" ]]; then
        log "Marking duplicate reads for ${label}"
        picard MarkDuplicates \
            INPUT="${sorted_bam}" \
            OUTPUT="${marked_bam}" \
            METRICS_FILE="${duplicate_metrics}" \
            REMOVE_DUPLICATES=false \
            ASSUME_SORT_ORDER=coordinate \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=SILENT \
            > "${LOG_ROOT}/${label}.picard.stdout.log" \
            2> "${LOG_ROOT}/${label}.picard.stderr.log"
    fi

    if [[ ! -s "${deduplicated_bam}" ]]; then
        samtools view \
            -@ "${THREADS}" \
            -b \
            -F 1024 \
            -o "${deduplicated_bam}" \
            "${marked_bam}"
        samtools index --threads "${THREADS}" "${deduplicated_bam}"
    fi

    samtools flagstat --threads "${THREADS}" "${marked_bam}" \
        > "${QC_ROOT}/alignment/${label}.marked.flagstat.txt"
    samtools flagstat --threads "${THREADS}" "${deduplicated_bam}" \
        > "${QC_ROOT}/alignment/${label}.nonduplicate.flagstat.txt"
    samtools idxstats "${marked_bam}" \
        > "${QC_ROOT}/alignment/${label}.marked.idxstats.tsv"
    samtools stats --threads "${THREADS}" "${marked_bam}" \
        > "${QC_ROOT}/alignment/${label}.marked.stats.txt"
}

align_runs() {
    local run study mark label

    [[ -s "${INDEX_PREFIX}.1.bt2" || -s "${INDEX_PREFIX}.1.bt2l" ]] || die "Bowtie2 index not found: ${INDEX_PREFIX}"

    while IFS='|' read -r run study mark label; do
        [[ -n "${run}" ]] || continue
        align_one_run "${run}" "${label}"
    done <<< "${RUN_MANIFEST}"

    log "Alignment and duplicate marking completed"
}

make_one_track() {
    local label="$1"
    local marked_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.marked.bam"
    local deduplicated_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.nonduplicate.bam"
    local all_reads_bigwig="${BIGWIG_ROOT}/${label}.q${MIN_MAPQ}.all_mapped.cpm.bw"
    local nonduplicate_bigwig="${BIGWIG_ROOT}/${label}.q${MIN_MAPQ}.nonduplicate.cpm.bw"

    [[ -s "${marked_bam}" ]] || die "Marked BAM not found: ${marked_bam}"
    [[ -s "${deduplicated_bam}" ]] || die "Nonduplicate BAM not found: ${deduplicated_bam}"

    if [[ ! -s "${all_reads_bigwig}" ]]; then
        "${BAM_COVERAGE}" \
            --bam "${marked_bam}" \
            --outFileName "${all_reads_bigwig}" \
            --outFileFormat bigwig \
            --normalizeUsing CPM \
            --binSize "${BIN_SIZE}" \
            --numberOfProcessors "${THREADS}"
    fi

    if [[ ! -s "${nonduplicate_bigwig}" ]]; then
        "${BAM_COVERAGE}" \
            --bam "${deduplicated_bam}" \
            --outFileName "${nonduplicate_bigwig}" \
            --outFileFormat bigwig \
            --normalizeUsing CPM \
            --binSize "${BIN_SIZE}" \
            --numberOfProcessors "${THREADS}"
    fi

    shasum -a 256 "${all_reads_bigwig}" "${nonduplicate_bigwig}" \
        > "${METADATA_ROOT}/${label}.bigwig_sha256.txt"
}

make_tracks() {
    local run study mark label

    while IFS='|' read -r run study mark label; do
        [[ -n "${run}" ]] || continue
        make_one_track "${label}"
    done <<< "${RUN_MANIFEST}"

    log "CPM bigWig generation completed"
}

analyze_reporter_loci() {
    log "Analyzing candidate reporter loci"
    Rscript "${SCRIPT_ROOT}/analyze_reporter_loci.R" \
        --work-root "${WORK_ROOT}"
    log "Candidate reporter-locus analysis completed"
}

usage() {
    cat <<EOF
Usage: $(basename "$0") [check|reference|download|trim|align|tracks|reporters|all]

Environment overrides:
  WORK_ROOT       Analysis directory (default: ${DEFAULT_WORK_ROOT})
  THREADS         Worker threads (default: 8)
  MIN_MAPQ        Minimum mapping quality retained (default: 20)
  BIN_SIZE        bigWig bin size in bp (default: 10)
  PYTHON_VERSION  pyenv Python containing deepTools (default: 3.12.13)
  TRIMMOMATIC_ADAPTERS
                  TruSeq single-end adapter FASTA, if Homebrew discovery fails

The all stage executes: check, reference, download, trim, align, tracks, reporters.
Existing non-empty downloads and outputs are retained to make reruns resumable.
EOF
}

main() {
    case "${STAGE}" in
        check)
            check_environment
            ;;
        reference)
            check_environment
            prepare_reference
            ;;
        download)
            check_environment
            download_runs
            ;;
        trim)
            check_environment
            trim_runs
            ;;
        align)
            check_environment
            align_runs
            ;;
        tracks)
            check_environment
            make_tracks
            ;;
        reporters)
            check_environment
            analyze_reporter_loci
            ;;
        all)
            check_environment
            prepare_reference
            download_runs
            trim_runs
            align_runs
            make_tracks
            analyze_reporter_loci
            ;;
        -h|--help|help)
            usage
            ;;
        *)
            usage >&2
            die "Unknown stage: ${STAGE}"
            ;;
    esac
}

main "$@"
