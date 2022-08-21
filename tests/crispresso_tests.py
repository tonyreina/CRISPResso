#!/usr/bin/env python
# -*- coding: utf8 -*-

# Run unit tests on CRISPResso using test data

# coverage run -m pytest -vvv crispresso_tests.py; coverage report -m

# CRISPResso --output_folder TEST_OUTPUT_FOLDER  -r1 test_data/test_L001_R1_001.fastq.gz -r2 test_data/test_L001_R2_001.fastq.gz --amplicon_seq gtcgcccctcaaatcttacagctgctcactcccctgcagggcaacgcccagggaccaagttagccccttaagcctaggcaaaagaatcccgcccataatcgagaagcgactcgacatggaggcgatgacgagatcacgcgaggaggaaaggagggagggcttcttccaggcccagggcggtccttacaagacgggaggcagcagagaactcccataaaggtattgcggcactcccctccccctgcccagaagggtgcggccttctctccacctcctccac --guide_seq aatcgagaagcgactcgaca,taaggggctaacttggtccc

import shutil
import numpy as np

import pytest

import CRISPResso as cr

#  -p 8
# --cleavage_offset -28
# --trim_sequences
# --trimmomatic_options_string
# "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36"
# --window_around_sgrna 70
class args:
    output_folder = "TEST_OUTPUT_FOLDER"
    fastq_r1 = "test_data/test_L001_R1_001.fastq.gz"
    fastq_r2 = "test_data/test_L001_R2_001.fastq.gz"
    amplicon_seq = (
        "gtcgcccctcaaatcttacagctgctcactc"
        "ccctgcagggcaacgcccagggaccaagttag"
        "ccccttaagcctaggcaaaagaatcccgccca"
        "taatcgagaagcgactcgacatggaggcgatg"
        "acgagatcacgcgaggaggaaaggagggaggg"
        "cttcttccaggcccagggcggtccttacaaga"
        "cgggaggcagcagagaactcccataaaggtatt"
        "gcggcactcccctccccctgcccagaagggtgc"
        "ggccttctctccacctcctccac"
    )
    guide_seq = "aatcgagaagcgactcgaca,taaggggctaacttggtccc"
    name = "pytest"
    cleavage_offset = -3
    window_around_sgrna = 1
    expected_hdr_amplicon_seq = ""
    donor_seq = ""
    coding_seq = ""
    min_average_read_quality = 0
    min_single_bp_quality = 0
    min_identity_score = 60.0
    split_paired_end = False
    trim_sequences = True
    trimmomatic_options_string = (
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True " " LEADING:3 TRAILING:3 MINLEN:36"
    )
    min_paired_end_reads_overlap = 4
    max_paired_end_reads_overlap = 100
    hide_mutations_outside_window_NHEJ = False
    exclude_bp_from_left = 15
    exclude_bp_from_right = 15
    hdr_perfect_alignment_threshold = 50.0
    ignore_substitutions = False
    ignore_deletions = False
    ignore_insertions = False
    needle_options_string = "-gapopen=10 -gapextend=0.5  -awidth3=5000"
    dump = False
    save_also_png = True
    n_processes = 1
    offset_around_cut_to_plot = 20
    min_frequency_alleles_around_cut_to_plot = 0.2
    max_rows_alleles_around_cut_to_plot = 50
    debug = False
    keep_intermediate = False


def test_count_reads():

    assert cr.get_n_reads_fastq("test_data/test_L001_R1_001.fastq.gz") == 8906
    assert cr.get_n_reads_fastq("test_data/test_L001_R2_001.fastq.gz") == 8906


def test_check_library():

    assert cr.check_library("pandas")
    assert cr.check_library("numpy")

    library_name = "#$231ddRRF^^&&*(( "
    with pytest.raises(Exception) as exc:
        cr.check_library(library_name)

    value = f"You need to install {library_name} to use CRISPResso!"

    assert value in str(exc.value)


def test_average_read_length():

    assert (
        cr.get_average_read_length_fastq("test_data/test_L001_R1_001.fastq.gz") == 151
    )
    assert (
        cr.get_average_read_length_fastq("test_data/test_L001_R2_001.fastq.gz") == 151
    )


def test_filter_se_fastq_by_qual():

    assert (
        cr.filter_se_fastq_by_qual("test_data/test_L001_R1_001.fastq.gz")
        == "test_data/test_L001_R1_001_filtered.fastq.gz"
    )


def test_filter_pe_fastq_by_qual():

    assert cr.filter_pe_fastq_by_qual(
        "test_data/test_L001_R1_001.fastq.gz", "test_data/test_L001_R2_001.fastq.gz"
    ) == (
        "test_data/test_L001_R1_001_filtered.fastq.gz",
        "test_data/test_L001_R2_001_filtered.fastq.gz",
    )


def test_get_ids_reads_to_remove():

    assert cr.get_ids_reads_to_remove("test_data/test_L001_R1_001.fastq.gz", 23) == set(
        [
            "M06879:15:000000000-DFF22:1:1101:25894:23776",
            "M06879:15:000000000-DFF22:1:1101:24046:20708",
        ]
    )
    assert cr.get_ids_reads_to_remove("test_data/test_L001_R2_001.fastq.gz", 15) == set(
        ["M06879:15:000000000-DFF22:1:1102:22078:15849"]
    )


def test_find_wrong_nt():

    assert (
        np.array(cr.find_wrong_nt("ACBTGCNGRCCACTGFNNC")).sort()
        == np.array(["R", "F", "B"]).sort()
    )


def test_reverse_complement():

    assert cr.reverse_complement("ACTGGT") == "ACCAGT"


def test_check_file_found():

    filename = "test_data/test_L001_R1_001.fastq.gz"
    cr.check_file(filename)

    filename = "test_data/test_L001_R2_001.fastq.gz"
    cr.check_file(filename)


def test_check_file_not_found():

    filename = "123test_dhjata/test_L0016hy_R1_001.hy.fastq.gz.in"
    with pytest.raises(Exception) as exc:
        cr.check_file(filename)

    assert f"I cannot open the file: {filename}" in str(exc.value)


def test_check_program():

    assert cr.check_program("date")


@pytest.mark.parametrize(
    "p,keep_intermediate,trim_sequences",
    [(1, True, True), (2, False, True), (8, True, False)],
)
def test_run_crispresso(p, keep_intermediate, trim_sequences):
    """Run end-to-end command"""

    args.fastq_r1 = "test_data/test_L001_R1_001.fastq.gz"
    args.fastq_r2 = "test_data/test_L001_R2_001.fastq.gz"
    args.n_processes = p
    args.keep_intermediate = keep_intermediate
    args.output_folder = f"pytest_p{p}"
    args.trim_sequences = trim_sequences

    (
        n_total,
        n_reads_input,
        n_unmodified,
        n_mixed_hdr_nhej,
        n_modified,
        n_repaired,
        nhej_inserted,
        nhej_deleted,
        nhej_mutated,
        df_indels,
        df_insertion,
        df_deletion,
        df_substitution,
        df_alleles,
    ) = cr.run_crispresso(args)

    assert n_total == 7058
    assert n_reads_input == 8906
    assert n_unmodified == 6853
    assert n_mixed_hdr_nhej == 0
    assert n_modified == 205
    assert n_repaired == 0
    assert nhej_inserted == 0
    assert nhej_deleted == 12
    assert nhej_mutated == 193

    assert tuple(df_indels["fq"].values[:4]) == (1, 0, 0, 0)
    assert tuple(df_insertion["fq"].values[:4]) == (7058, 0, 0, 0)
    assert tuple(df_deletion["fq"].values[:4]) == (7046, 0, 0, 0)
    assert tuple(df_substitution["fq"].values[:4]) == (6865, 188, 5, 0)
    assert tuple(df_alleles["#Reads"].values[:4]) == (1098, 346, 19, 17)
