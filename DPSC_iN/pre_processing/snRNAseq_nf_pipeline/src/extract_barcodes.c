/*
 * extract_barcodes.c
 *
 * Extracts cell barcodes, UMIs, and lineage barcodes from
 * paired-end snRNA-seq FASTQ files. Uses kseq.h for FASTQ parsing.
 *
 * Usage: extract_barcodes <sample_id> <fastq_1.fastq.gz> [fastq_2.fastq.gz ...]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

/* ── Constants ──────────────────────────────────────────── */
#define CELL_BC_LEN   16
#define UMI_LEN       12
#define REVP_LEN      18
#define VIRUSBC_LEN   20
#define SLIBID_LEN    6
#define GLIBID_LEN    6
#define MIN_R2_LEN    50
#define EXPECTED_LIB  "GGCTAC"
#define MAX_HAMMING   1

/* ── Output file handles ────────────────────────────────── */
typedef struct {
    FILE *barcodes;     /* {sample_id}_barcodes.tsv  */
    FILE *non_lib;      /* {sample_id}_non_libreads.fa */
} OutputFiles;

/* ── Extracted fields from one read pair ───────────────── */
typedef struct {
    char cellbarcode[CELL_BC_LEN + 1];
    char umi        [UMI_LEN     + 1];
    char revp       [REVP_LEN    + 1];
    char virusbc    [VIRUSBC_LEN + 1];
    char slibids    [SLIBID_LEN  + 1];
    char glibids    [GLIBID_LEN  + 1];
} ReadFields;

/* ── Lane pairing structure ─────────────────────────────── */
typedef struct {
    int lane;
    char *r1;
    char *r2;
} LanePair;

/* ══════════════════════════════════════════════════════════
 * hamming_distance
 * Returns the number of differing characters between two
 * fixed-length strings. Both strings must be >= len chars.
 * ══════════════════════════════════════════════════════════ */
static int hamming_distance(const char *s1, const char *s2, int len) {
    int dist = 0;
    for (int i = 0; i < len; i++)
        if (s1[i] != s2[i]) dist++;
    return dist;
}

/* ══════════════════════════════════════════════════════════
 * get_lane_number
 * Extracts the lane number from a filename containing _L00N_.
 * Returns the lane integer, or -1 if not found.
 * ══════════════════════════════════════════════════════════ */
static int get_lane_number(const char *filename) {
    const char *p = strstr(filename, "_L00");
    if (!p) return -1;
    return atoi(p + 4);
}

/* ══════════════════════════════════════════════════════════
 * get_read_number
 * Extracts whether a filename corresponds to R1 or R2.
 * Returns 1 for R1, 2 for R2, or -1 if not found.
 * ══════════════════════════════════════════════════════════ */
static int get_read_number(const char *filename) {
    if (strstr(filename, "_R1_") || strstr(filename, "_R1."))
        return 1;
    if (strstr(filename, "_R2_") || strstr(filename, "_R2."))
        return 2;
    return -1;
}

/* ══════════════════════════════════════════════════════════
 * build_lane_pairs
 * Groups input files by lane and assigns R1/R2 members based
 * on filename, independent of input position.
 * Returns 1 on success, 0 on error.
 * ══════════════════════════════════════════════════════════ */
static int build_lane_pairs(char **files, int n_files,
                            LanePair **pairs_out, int *n_pairs_out) {
    LanePair *pairs = calloc(n_files, sizeof(LanePair));
    if (!pairs) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        return 0;
    }

    int n_pairs = 0;

    for (int i = 0; i < n_files; i++) {
        int lane = get_lane_number(files[i]);
        int read = get_read_number(files[i]);

        if (lane == -1) {
            fprintf(stderr, "Error: Could not extract lane number from file: %s\n", files[i]);
            free(pairs);
            return 0;
        }

        if (read == -1) {
            fprintf(stderr, "Error: Could not determine whether file is R1 or R2: %s\n", files[i]);
            free(pairs);
            return 0;
        }

        int idx = -1;
        for (int j = 0; j < n_pairs; j++) {
            if (pairs[j].lane == lane) {
                idx = j;
                break;
            }
        }

        if (idx == -1) {
            idx = n_pairs++;
            pairs[idx].lane = lane;
            pairs[idx].r1 = NULL;
            pairs[idx].r2 = NULL;
        }

        if (read == 1) {
            if (pairs[idx].r1 != NULL) {
                fprintf(stderr, "Error: Multiple R1 files found for lane %d:\n  %s\n  %s\n",
                        lane, pairs[idx].r1, files[i]);
                free(pairs);
                return 0;
            }
            pairs[idx].r1 = files[i];
        } else {
            if (pairs[idx].r2 != NULL) {
                fprintf(stderr, "Error: Multiple R2 files found for lane %d:\n  %s\n  %s\n",
                        lane, pairs[idx].r2, files[i]);
                free(pairs);
                return 0;
            }
            pairs[idx].r2 = files[i];
        }
    }

    for (int i = 0; i < n_pairs; i++) {
        if (pairs[i].r1 == NULL || pairs[i].r2 == NULL) {
            fprintf(stderr, "Error: Missing mate for lane %d:\n", pairs[i].lane);
            if (pairs[i].r1)
                fprintf(stderr, "  R1: %s\n", pairs[i].r1);
            if (pairs[i].r2)
                fprintf(stderr, "  R2: %s\n", pairs[i].r2);
            free(pairs);
            return 0;
        }
    }

    *pairs_out = pairs;
    *n_pairs_out = n_pairs;
    return 1;
}

/* ══════════════════════════════════════════════════════════
 * open_output_files
 * Opens barcodes TSV and non-lib FASTA output files, both
 * prefixed with sample_id for uniqueness.
 * Returns 1 on success, 0 on failure.
 * ══════════════════════════════════════════════════════════ */
static int open_output_files(const char *sample_id, OutputFiles *out) {
    char path[512];

    snprintf(path, sizeof(path), "%s_barcodes.tsv", sample_id);
    out->barcodes = fopen(path, "w");
    if (!out->barcodes) {
        fprintf(stderr, "Error: Cannot open '%s' for writing.\n", path);
        return 0;
    }

    snprintf(path, sizeof(path), "%s_non_libreads.fa", sample_id);
    out->non_lib = fopen(path, "w");
    if (!out->non_lib) {
        fprintf(stderr, "Error: Cannot open '%s' for writing.\n", path);
        fclose(out->barcodes);
        return 0;
    }

    /* Write TSV header once */
    fprintf(out->barcodes,
        "CELLBARCODE\tUMI\tREVERSEPRIMER\tVIRUSBC\t"
        "SPECIFICLIBID\tGENERALLIBID\tSAMPLE\tFWDREAD\tREVREAD\n");

    return 1;
}

/* ══════════════════════════════════════════════════════════
 * close_output_files
 * Closes both output file handles.
 * ══════════════════════════════════════════════════════════ */
static void close_output_files(OutputFiles *out) {
    if (out->barcodes) fclose(out->barcodes);
    if (out->non_lib)  fclose(out->non_lib);
}

/* ══════════════════════════════════════════════════════════
 * extract_fields
 * Populates a ReadFields struct from a pair of raw sequences.
 * Returns 1 if reads are long enough, 0 if too short to use.
 * ══════════════════════════════════════════════════════════ */
static int extract_fields(const char *r1_seq, size_t r1_len,
                           const char *r2_seq, size_t r2_len,
                           ReadFields *f) {
    if (r1_len < (size_t)(CELL_BC_LEN + UMI_LEN)) return 0;
    if (r2_len < MIN_R2_LEN)                       return 0;

    strncpy(f->cellbarcode, r1_seq,               CELL_BC_LEN);
    strncpy(f->umi,         r1_seq + CELL_BC_LEN, UMI_LEN);
    strncpy(f->revp,        r2_seq,               REVP_LEN);
    strncpy(f->virusbc,     r2_seq + 19,          VIRUSBC_LEN);
    strncpy(f->slibids,     r2_seq + 38,          SLIBID_LEN);
    strncpy(f->glibids,     r2_seq + 44,          GLIBID_LEN);

    f->cellbarcode[CELL_BC_LEN] = '\0';
    f->umi        [UMI_LEN]     = '\0';
    f->revp       [REVP_LEN]    = '\0';
    f->virusbc    [VIRUSBC_LEN] = '\0';
    f->slibids    [SLIBID_LEN]  = '\0';
    f->glibids    [GLIBID_LEN]  = '\0';

    return 1;
}

/* ══════════════════════════════════════════════════════════
 * write_read
 * Writes one read pair to either the barcodes TSV (match) or
 * the non-lib FASTA (no match), based on Hamming distance.
 * ══════════════════════════════════════════════════════════ */
static void write_read(const ReadFields *f, const char *sample_id,
                        const char *r1_seq, const char *r2_seq,
                        OutputFiles *out) {
    int dist = hamming_distance(f->glibids, EXPECTED_LIB, GLIBID_LEN);

    if (dist <= MAX_HAMMING) {
        fprintf(out->barcodes,
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            f->cellbarcode, f->umi, f->revp, f->virusbc,
            f->slibids, f->glibids, sample_id,
            r1_seq, r2_seq);
    } else {
        /* Trim first 18 bases of R1 before writing to non-lib output */
        fprintf(out->non_lib, ">%s\n%s\n",
            f->cellbarcode,
            r1_seq + 18);
    }
}

/* ══════════════════════════════════════════════════════════
 * process_lane_pair
 * Opens one R1/R2 lane pair with kseq, iterates over all
 * read pairs, and writes results to the output files.
 * Returns 1 on success, 0 on error.
 * ══════════════════════════════════════════════════════════ */
static int process_lane_pair(const char *r1_path, const char *r2_path,
                              const char *sample_id, OutputFiles *out) {
    gzFile fp_r1 = gzopen(r1_path, "r");
    if (!fp_r1) {
        fprintf(stderr, "Error: Cannot open R1 file '%s'\n", r1_path);
        return 0;
    }

    gzFile fp_r2 = gzopen(r2_path, "r");
    if (!fp_r2) {
        fprintf(stderr, "Error: Cannot open R2 file '%s'\n", r2_path);
        gzclose(fp_r1);
        return 0;
    }

    kseq_t *r1 = kseq_init(fp_r1);
    kseq_t *r2 = kseq_init(fp_r2);

    long n_reads = 0, n_matched = 0, n_short = 0;
    int ret_r1, ret_r2;

    while ((ret_r1 = kseq_read(r1)) >= 0 &&
           (ret_r2 = kseq_read(r2)) >= 0) {
        n_reads++;

        ReadFields f;
        if (!extract_fields(r1->seq.s, r1->seq.l,
                             r2->seq.s, r2->seq.l, &f)) {
            n_short++;
            continue;
        }

        int dist = hamming_distance(f.glibids, EXPECTED_LIB, GLIBID_LEN);
        if (dist <= MAX_HAMMING) n_matched++;

        write_read(&f, sample_id, r1->seq.s, r2->seq.s, out);
    }

    /* Warn if R1 and R2 had different read counts */
    if ((ret_r1 >= 0) != (ret_r2 >= 0)) {
        fprintf(stderr,
            "Warning: R1 and R2 have different read counts in:\n"
            "  R1: %s\n  R2: %s\n", r1_path, r2_path);
    }

    fprintf(stderr,
        "  Lane summary: %ld reads | %ld matched | %ld too short\n",
        n_reads, n_matched, n_short);

    kseq_destroy(r1);
    kseq_destroy(r2);
    gzclose(fp_r1);
    gzclose(fp_r2);
    return 1;
}

/* ══════════════════════════════════════════════════════════
 * print_usage
 * ══════════════════════════════════════════════════════════ */
static void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s <sample_id> <fastq_1.fastq.gz> [fastq_2.fastq.gz ...]\n"
        "\n"
        "  Provide all FASTQ files for a sample across lanes.\n"
        "  R1/R2 pairing is determined from filenames, not argument order.\n"
        "\n"
        "Options:\n"
        "  -h  Show this help message\n"
        "\n"
        "Output files (prefixed with sample_id):\n"
        "  {sample_id}_barcodes.tsv     matched reads with all fields\n"
        "  {sample_id}_non_libreads.fa  non-matching reads in FASTA\n",
        prog);
}

/* ══════════════════════════════════════════════════════════
 * main
 * ══════════════════════════════════════════════════════════ */
int main(int argc, char *argv[]) {

    /* ── Argument parsing ─────────────────────────────── */
    if (argc < 4) { print_usage(argv[0]); return 1; }
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0) { print_usage(argv[0]); return 0; }
    }

    const char *sample_id = argv[1];

    int n_files = argc - 2;
    char **files = &argv[2];

    LanePair *pairs = NULL;
    int n_pairs = 0;

    /* ── Build lane pairs from filenames ──────────────── */
    if (!build_lane_pairs(files, n_files, &pairs, &n_pairs))
        return 1;

    /* ── Open output files ────────────────────────────── */
    OutputFiles out;
    if (!open_output_files(sample_id, &out)) {
        free(pairs);
        return 1;
    }

    /* ── Process each lane pair ───────────────────────── */
    fprintf(stderr, "Sample: %s | %d lane pair(s)\n", sample_id, n_pairs);

    int success = 1;
    for (int i = 0; i < n_pairs; i++) {
        fprintf(stderr, "Lane %d: %s + %s\n",
                pairs[i].lane, pairs[i].r1, pairs[i].r2);
        if (!process_lane_pair(pairs[i].r1, pairs[i].r2, sample_id, &out)) {
            success = 0;
            break;
        }
    }

    /* ── Cleanup ──────────────────────────────────────── */
    close_output_files(&out);
    free(pairs);

    if (success)
        fprintf(stderr, "Done. Output: %s_barcodes.tsv, %s_non_libreads.fa\n",
                sample_id, sample_id);

    return success ? 0 : 1;
}