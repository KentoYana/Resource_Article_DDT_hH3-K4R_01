#!/usr/bin/env bash

# Re-download and uniformly process selected public N. crassa datasets for
# reporter-locus nucleosome occupancy and chromatin accessibility analysis.
#
# Approved datasets for G3 revision issue #69:
#   Kamei et al. 2021: WT MNase-seq, three biological replicates
#   Kamei et al. 2021: WT total-H3 ChIP-seq, two biological replicates
#   Ferraro et al. 2021: WT ATAC-seq, two biological replicates
#
# Studies and assays remain separate. This workflow does not pool runs across
# studies, call peaks, or treat ATAC-seq as a direct nucleosome assay.

set -Eeuo pipefail
IFS=$'\n\t'

readonly DEFAULT_WORK_ROOT="/Volumes/Garage/Re_analysis/260907_issue69_nucleosome"
readonly REFERENCE_ACCESSION="GCA_000182925.2"
readonly REFERENCE_BASENAME="GCA_000182925.2_NC12"
readonly REFERENCE_URL_ROOT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/925/GCA_000182925.2_NC12"
readonly SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORK_ROOT="${WORK_ROOT:-${DEFAULT_WORK_ROOT}}"
THREADS="${THREADS:-8}"
MIN_MAPQ="${MIN_MAPQ:-20}"
MNASE_MIN_FRAGMENT="${MNASE_MIN_FRAGMENT:-130}"
MNASE_MAX_FRAGMENT="${MNASE_MAX_FRAGMENT:-200}"
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

readonly NUCLEAR_CONTIGS=(
    "CM002236.1" "CM002237.1" "CM002238.1" "CM002239.1"
    "CM002240.1" "CM002241.1" "CM002242.1"
)

# run|study|assay|replicate|sample label
readonly RUN_MANIFEST=$(cat <<'EOF'
SRR13067301|Kamei2021|MNase|1|Kamei2021_WT_MNase_rep1
SRR13067302|Kamei2021|MNase|2|Kamei2021_WT_MNase_rep2
SRR13067303|Kamei2021|MNase|3|Kamei2021_WT_MNase_rep3
SRR13067307|Kamei2021|total_H3|1|Kamei2021_WT_totalH3_rep1
SRR13067308|Kamei2021|total_H3|2|Kamei2021_WT_totalH3_rep2
SRR12229299|Ferraro2021|ATAC|1|Ferraro2021_WT_ATAC_rep1
SRR12229300|Ferraro2021|ATAC|2|Ferraro2021_WT_ATAC_rep2
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
        "${QC_ROOT}/fragment_size" \
        "${LOG_ROOT}" \
        "${METADATA_ROOT}" \
        "${TMP_ROOT}/fasterq" \
        "${TMP_ROOT}/matplotlib" \
        "${TMP_ROOT}/atac_shift"
}

resolve_tools() {
    local command_name

    for command_name in \
        prefetch fasterq-dump vdb-validate fastqc trimmomatic bowtie2 Rscript \
        bowtie2-build samtools picard curl gzip shasum awk; do
        require_command "${command_name}"
    done

    require_command pyenv
    export PYENV_VERSION="${PYTHON_VERSION}"
    export MPLCONFIGDIR="${TMP_ROOT}/matplotlib"

    BAM_COVERAGE=$(pyenv which bamCoverage) || \
        die "bamCoverage is unavailable in pyenv Python ${PYTHON_VERSION}"
    ALIGNMENT_SIEVE=$(pyenv which alignmentSieve) || \
        die "alignmentSieve is unavailable in pyenv Python ${PYTHON_VERSION}"
    [[ -x "${BAM_COVERAGE}" ]] || die "bamCoverage is not executable: ${BAM_COVERAGE}"
    [[ -x "${ALIGNMENT_SIEVE}" ]] || die "alignmentSieve is not executable: ${ALIGNMENT_SIEVE}"
    "${BAM_COVERAGE}" --version >/dev/null
    "${ALIGNMENT_SIEVE}" --version >/dev/null

    Rscript -e '
        required <- c(
            "tidyverse", "tikzDevice", "patchwork", "here",
            "IRanges", "rtracklayer", "digest", "jsonlite"
        )
        missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
        if (length(missing) > 0L) {
            stop("Missing required R packages: ", paste(missing, collapse = ", "))
        }
    ' >/dev/null || die "Required R packages are unavailable"

    if command -v brew >/dev/null 2>&1; then
        TRIMMOMATIC_ADAPTERS="${TRIMMOMATIC_ADAPTERS:-$(brew --prefix trimmomatic)/share/trimmomatic/adapters/TruSeq3-PE.fa}"
    else
        TRIMMOMATIC_ADAPTERS="${TRIMMOMATIC_ADAPTERS:-}"
    fi
    [[ -f "${TRIMMOMATIC_ADAPTERS}" ]] || \
        die "Paired-end Trimmomatic adapter FASTA not found; set TRIMMOMATIC_ADAPTERS"
}

record_manifest() {
    {
        printf 'run\tstudy\tassay\treplicate\tsample_label\n'
        printf '%s\n' "${RUN_MANIFEST}" | tr '|' '\t'
    } > "${METADATA_ROOT}/selected_runs.tsv"
}

record_versions() {
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
        "${ALIGNMENT_SIEVE}" --version 2>&1 | sed -n '1p' || true
        Rscript --version 2>&1 | sed -n '1p' || true
    } > "${METADATA_ROOT}/software_versions.txt"
}

check_environment() {
    [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
    [[ "${MIN_MAPQ}" =~ ^[0-9]+$ ]] || die "MIN_MAPQ must be a non-negative integer"
    [[ "${MNASE_MIN_FRAGMENT}" =~ ^[1-9][0-9]*$ ]] || die "MNASE_MIN_FRAGMENT must be positive"
    [[ "${MNASE_MAX_FRAGMENT}" =~ ^[1-9][0-9]*$ ]] || die "MNASE_MAX_FRAGMENT must be positive"
    (( MNASE_MIN_FRAGMENT <= MNASE_MAX_FRAGMENT )) || \
        die "MNASE_MIN_FRAGMENT must not exceed MNASE_MAX_FRAGMENT"

    prepare_directories
    resolve_tools
    record_manifest
    record_versions

    log "Work root: ${WORK_ROOT}"
    log "Free space at work root:"
    df -h "${WORK_ROOT}"
    log "pyenv bamCoverage: ${BAM_COVERAGE}"
    log "pyenv alignmentSieve: ${ALIGNMENT_SIEVE}"
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

    local contig
    for contig in "${NUCLEAR_CONTIGS[@]}"; do
        grep -q "^${contig}"$'\t' "${REFERENCE_FASTA}.fai" || \
            die "Expected NC12 nuclear contig absent from FASTA: ${contig}"
    done

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

download_one_run() {
    local run="$1"
    local sra_path="${SRA_ROOT}/${run}/${run}.sra"
    local read1="${RAW_FASTQ_ROOT}/${run}_1.fastq"
    local read2="${RAW_FASTQ_ROOT}/${run}_2.fastq"
    local read1_gz="${read1}.gz"
    local read2_gz="${read2}.gz"
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

    if [[ ! -s "${read1_gz}" || ! -s "${read2_gz}" ]]; then
        [[ ! -e "${read1_gz}" && ! -e "${read2_gz}" ]] || \
            die "Only one compressed FASTQ mate exists for ${run}"
        mkdir -p "${run_tmp}"
        if [[ ! -s "${read1}" || ! -s "${read2}" ]]; then
            log "Converting ${run} to paired-end FASTQ"
            fasterq-dump "${sra_path}" \
                --threads "${THREADS}" \
                --temp "${run_tmp}" \
                --outdir "${RAW_FASTQ_ROOT}" \
                --split-files \
                > "${LOG_ROOT}/${run}.fasterq-dump.stdout.log" \
                2> "${LOG_ROOT}/${run}.fasterq-dump.stderr.log"
        fi
        [[ -s "${read1}" && -s "${read2}" ]] || \
            die "Expected paired FASTQ files were not created for ${run}"
        compress_fastq "${read1}"
        compress_fastq "${read2}"
    fi

    shasum -a 256 "${sra_path}" "${read1_gz}" "${read2_gz}" \
        > "${METADATA_ROOT}/${run}.download_sha256.txt"
    fastqc --threads "${THREADS}" --outdir "${QC_ROOT}/fastqc_raw" \
        "${read1_gz}" "${read2_gz}" \
        > "${LOG_ROOT}/${run}.fastqc_raw.stdout.log" \
        2> "${LOG_ROOT}/${run}.fastqc_raw.stderr.log"
}

download_runs() {
    local run study assay replicate label

    while IFS='|' read -r run study assay replicate label; do
        [[ -n "${run}" ]] || continue
        download_one_run "${run}"
    done <<< "${RUN_MANIFEST}"

    log "SRA and paired FASTQ download completed"
}

trim_one_run() {
    local run="$1"
    local raw1="${RAW_FASTQ_ROOT}/${run}_1.fastq.gz"
    local raw2="${RAW_FASTQ_ROOT}/${run}_2.fastq.gz"
    local paired1="${TRIMMED_FASTQ_ROOT}/${run}_1.paired.fastq.gz"
    local unpaired1="${TRIMMED_FASTQ_ROOT}/${run}_1.unpaired.fastq.gz"
    local paired2="${TRIMMED_FASTQ_ROOT}/${run}_2.paired.fastq.gz"
    local unpaired2="${TRIMMED_FASTQ_ROOT}/${run}_2.unpaired.fastq.gz"

    [[ -s "${raw1}" && -s "${raw2}" ]] || die "Raw paired FASTQ not found for ${run}"

    if [[ ! -s "${paired1}" || ! -s "${paired2}" ]]; then
        log "Trimming ${run}"
        trimmomatic PE \
            -threads "${THREADS}" \
            -phred33 \
            "${raw1}" "${raw2}" \
            "${paired1}" "${unpaired1}" \
            "${paired2}" "${unpaired2}" \
            "ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:2:30:10" \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 \
            > "${LOG_ROOT}/${run}.trimmomatic.stdout.log" \
            2> "${LOG_ROOT}/${run}.trimmomatic.stderr.log"
    fi

    shasum -a 256 "${paired1}" "${paired2}" \
        > "${METADATA_ROOT}/${run}.trimmed_fastq_sha256.txt"
    fastqc --threads "${THREADS}" --outdir "${QC_ROOT}/fastqc_trimmed" \
        "${paired1}" "${paired2}" \
        > "${LOG_ROOT}/${run}.fastqc_trimmed.stdout.log" \
        2> "${LOG_ROOT}/${run}.fastqc_trimmed.stderr.log"
}

trim_runs() {
    local run study assay replicate label

    while IFS='|' read -r run study assay replicate label; do
        [[ -n "${run}" ]] || continue
        trim_one_run "${run}"
    done <<< "${RUN_MANIFEST}"

    log "Paired-read trimming completed"
}

write_fragment_histogram() {
    local bam="$1"
    local output="$2"

    {
        printf 'fragment_length_bp\tpair_count\n'
        samtools view -@ "${THREADS}" -f 64 "${bam}" \
        | awk '
            BEGIN { OFS = "\t" }
            {
                fragment_length = $9
                if (fragment_length < 0) fragment_length = -fragment_length
                if (fragment_length > 0) count[fragment_length]++
            }
            END {
                for (fragment_length in count) print fragment_length, count[fragment_length]
            }
        ' \
        | sort -n -k1,1
    } > "${output}"
}

align_one_run() {
    local run="$1"
    local label="$2"
    local read1="${TRIMMED_FASTQ_ROOT}/${run}_1.paired.fastq.gz"
    local read2="${TRIMMED_FASTQ_ROOT}/${run}_2.paired.fastq.gz"
    local aligned_bam="${BAM_ROOT}/${label}.aligned.sorted.bam"
    local q20_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.sorted.bam"
    local marked_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.marked.bam"
    local deduplicated_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.nonduplicate.bam"
    local duplicate_metrics="${QC_ROOT}/duplicates/${label}.picard_markduplicates.txt"

    [[ -s "${read1}" && -s "${read2}" ]] || die "Trimmed paired FASTQ not found for ${run}"

    if [[ ! -s "${aligned_bam}" ]]; then
        log "Aligning ${run} (${label})"
        bowtie2 \
            --very-sensitive \
            --no-mixed \
            --no-discordant \
            -X 1000 \
            -p "${THREADS}" \
            -x "${INDEX_PREFIX}" \
            -1 "${read1}" \
            -2 "${read2}" \
            2> "${LOG_ROOT}/${label}.bowtie2.stderr.log" \
        | samtools sort \
            -@ "${THREADS}" \
            -O BAM \
            -o "${aligned_bam}" \
            -
        samtools index --threads "${THREADS}" "${aligned_bam}"
    fi

    if [[ ! -s "${q20_bam}" ]]; then
        local q20_without_rg="${q20_bam}.without_read_group.bam"
        samtools view \
            -@ "${THREADS}" \
            -b \
            -q "${MIN_MAPQ}" \
            -f 2 \
            -F 2828 \
            -o "${q20_without_rg}" \
            "${aligned_bam}" \
            "${NUCLEAR_CONTIGS[@]}"
        samtools addreplacerg \
            -@ "${THREADS}" \
            -r "ID:${run}" \
            -r "SM:${label}" \
            -r "LB:${label}" \
            -r "PL:ILLUMINA" \
            -o "${q20_bam}" \
            "${q20_without_rg}"
        rm "${q20_without_rg}"
        samtools index --threads "${THREADS}" "${q20_bam}"
    fi

    if ! samtools view -H "${q20_bam}" | grep -q '^@RG'; then
        local q20_with_rg="${q20_bam}.with_read_group.bam"
        log "Adding read group to existing BAM for ${label}"
        samtools addreplacerg \
            -@ "${THREADS}" \
            -r "ID:${run}" \
            -r "SM:${label}" \
            -r "LB:${label}" \
            -r "PL:ILLUMINA" \
            -o "${q20_with_rg}" \
            "${q20_bam}"
        mv "${q20_with_rg}" "${q20_bam}"
        samtools index --threads "${THREADS}" "${q20_bam}"
    fi

    if [[ ! -s "${marked_bam}" ]]; then
        log "Marking duplicate pairs for ${label}"
        picard MarkDuplicates \
            INPUT="${q20_bam}" \
            OUTPUT="${marked_bam}" \
            METRICS_FILE="${duplicate_metrics}" \
            REMOVE_DUPLICATES=false \
            ASSUME_SORT_ORDER=coordinate \
            CREATE_INDEX=true \
            READ_NAME_REGEX=null \
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
    write_fragment_histogram \
        "${marked_bam}" \
        "${QC_ROOT}/fragment_size/${label}.marked.fragment_lengths.tsv"
}

align_runs() {
    local run study assay replicate label

    [[ -s "${INDEX_PREFIX}.1.bt2" || -s "${INDEX_PREFIX}.1.bt2l" ]] || \
        die "Bowtie2 index not found: ${INDEX_PREFIX}"

    while IFS='|' read -r run study assay replicate label; do
        [[ -n "${run}" ]] || continue
        align_one_run "${run}" "${label}"
    done <<< "${RUN_MANIFEST}"

    log "Paired-end alignment and duplicate marking completed"
}

make_atac_track() {
    local input_bam="$1"
    local output_bigwig="$2"
    local shift_stem="$3"
    local high_resolution_bigwig="${output_bigwig%.cpm.bw}.cutsite1bp.cpm.bw"
    local shifted_unsorted="${TMP_ROOT}/atac_shift/${shift_stem}.unsorted.bam"
    local shifted_bam="${TMP_ROOT}/atac_shift/${shift_stem}.bam"

    if [[ ! -s "${shifted_bam}" ]]; then
        "${ALIGNMENT_SIEVE}" \
            --bam "${input_bam}" \
            --outFile "${shifted_unsorted}" \
            --ATACshift \
            --numberOfProcessors "${THREADS}" \
            --filterMetrics "${QC_ROOT}/alignment/${shift_stem}.atac_shift.metrics.txt"
        samtools sort -@ "${THREADS}" -o "${shifted_bam}" "${shifted_unsorted}"
        samtools index --threads "${THREADS}" "${shifted_bam}"
        rm "${shifted_unsorted}"
    fi

    if [[ ! -s "${high_resolution_bigwig}" ]]; then
        "${BAM_COVERAGE}" \
            --bam "${shifted_bam}" \
            --outFileName "${high_resolution_bigwig}" \
            --outFileFormat bigwig \
            --normalizeUsing CPM \
            --exactScaling \
            --binSize 1 \
            --Offset 1 \
            --numberOfProcessors "${THREADS}"
    fi

    "${BAM_COVERAGE}" \
        --bam "${shifted_bam}" \
        --outFileName "${output_bigwig}" \
        --outFileFormat bigwig \
        --normalizeUsing CPM \
        --exactScaling \
        --binSize 10 \
        --Offset 1 \
        --numberOfProcessors "${THREADS}"
}

make_one_track() {
    local assay="$1"
    local label="$2"
    local variant="$3"
    local input_bam
    local output_bigwig="${BIGWIG_ROOT}/${label}.q${MIN_MAPQ}.${variant}.cpm.bw"

    if [[ "${variant}" == "all_mapped" ]]; then
        input_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.marked.bam"
    else
        input_bam="${BAM_ROOT}/${label}.q${MIN_MAPQ}.nonduplicate.bam"
    fi
    [[ -s "${input_bam}" ]] || die "Input BAM not found: ${input_bam}"

    if [[ -s "${output_bigwig}" ]]; then
        log "Keeping existing track: ${output_bigwig}"
        return
    fi

    log "Generating ${assay} ${variant} CPM track for ${label}"
    case "${assay}" in
        MNase)
            "${BAM_COVERAGE}" \
                --bam "${input_bam}" \
                --outFileName "${output_bigwig}" \
                --outFileFormat bigwig \
                --normalizeUsing CPM \
                --exactScaling \
                --binSize 1 \
                --MNase \
                --minFragmentLength "${MNASE_MIN_FRAGMENT}" \
                --maxFragmentLength "${MNASE_MAX_FRAGMENT}" \
                --samFlagInclude 64 \
                --numberOfProcessors "${THREADS}"
            ;;
        total_H3)
            "${BAM_COVERAGE}" \
                --bam "${input_bam}" \
                --outFileName "${output_bigwig}" \
                --outFileFormat bigwig \
                --normalizeUsing CPM \
                --exactScaling \
                --binSize 10 \
                --extendReads \
                --samFlagInclude 64 \
                --numberOfProcessors "${THREADS}"
            ;;
        ATAC)
            make_atac_track \
                "${input_bam}" \
                "${output_bigwig}" \
                "${label}.${variant}"
            ;;
        *)
            die "Unknown assay for track generation: ${assay}"
            ;;
    esac
}

make_tracks() {
    local run study assay replicate label variant

    while IFS='|' read -r run study assay replicate label; do
        [[ -n "${run}" ]] || continue
        for variant in all_mapped nonduplicate; do
            make_one_track "${assay}" "${label}" "${variant}"
        done
        if [[ "${assay}" == "ATAC" ]]; then
            shasum -a 256 \
                "${BIGWIG_ROOT}/${label}.q${MIN_MAPQ}.all_mapped.cpm.bw" \
                "${BIGWIG_ROOT}/${label}.q${MIN_MAPQ}.nonduplicate.cpm.bw" \
                "${BIGWIG_ROOT}/${label}.q${MIN_MAPQ}.all_mapped.cutsite1bp.cpm.bw" \
                "${BIGWIG_ROOT}/${label}.q${MIN_MAPQ}.nonduplicate.cutsite1bp.cpm.bw" \
                > "${METADATA_ROOT}/${label}.bigwig_sha256.txt"
        else
            shasum -a 256 \
                "${BIGWIG_ROOT}/${label}.q${MIN_MAPQ}.all_mapped.cpm.bw" \
                "${BIGWIG_ROOT}/${label}.q${MIN_MAPQ}.nonduplicate.cpm.bw" \
                > "${METADATA_ROOT}/${label}.bigwig_sha256.txt"
        fi
    done <<< "${RUN_MANIFEST}"

    log "MNase, total-H3, and ATAC CPM bigWig generation completed"
}

analyze_reporter_loci() {
    log "Analyzing candidate reporter loci for nucleosome occupancy"
    Rscript "${SCRIPT_ROOT}/analyze_reporter_loci_nucleosome.R" \
        --work-root "${WORK_ROOT}"
    log "Nucleosome reporter-locus analysis completed"
}

usage() {
    cat <<EOF
Usage: $(basename "$0") [check|reference|download|trim|align|tracks|reporters|all]

Environment overrides:
  WORK_ROOT          Analysis directory (default: ${DEFAULT_WORK_ROOT})
  THREADS            Worker threads (default: 8)
  MIN_MAPQ           Minimum mapping quality retained (default: 20)
  MNASE_MIN_FRAGMENT Minimum MNase fragment length (default: 130)
  MNASE_MAX_FRAGMENT Maximum MNase fragment length (default: 200)
  PYTHON_VERSION     pyenv Python containing deepTools (default: 3.12.13)
  TRIMMOMATIC_ADAPTERS
                     TruSeq paired-end adapter FASTA, if discovery fails

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
