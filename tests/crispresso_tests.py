#!/usr/bin/env python
# -*- coding: utf8 -*-

# Run unit tests on CRISPResso using test data

# coverage run -m pytest -vvv crispresso_tests.py; coverage report -m

import sys
import os
import numpy as np

import pytest

import CRISPResso as cr

env_directory = os.getenv("CONDA_PREFIX")
trim_dir = f"{env_directory}/share/trimmomatic-0.39-2/adapters/"

# Next line hijacks the command line parameters so that
# we can mock our own during the test.
sys.argv = [
    "CRISPResso",
    "-r1",
    "test_data/test1_L001_R1_001.fastq.gz",
    "--amplicon_seq",
    "gtcgccccgacttctctccacctcctccac",
]


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

    filename = ".123test_dhjata/test_L0.016hy_R1_001.hy.fastq.gz.in"
    with pytest.raises(Exception) as exc:
        cr.check_file(filename)

    assert f"I cannot open the file: {filename}" in str(exc.value)


def test_check_program():

    assert cr.check_program("date")


@pytest.mark.parametrize(
    "p,keep_intermediate",
    [(1, False), (2, True)],
)
def test_run_crispresso(p, keep_intermediate):
    """Run end-to-end command

    Ground truth values are from the original CRIPResso Docker
    CRISPResso --output_folder TEST_OUTPUT_FOLDER  \
        -r1 test_data/test_L001_R1_001.fastq.gz \
         -r2 test_data/test_L001_R2_001.fastq.gz \
         --amplicon_seq gtcgcccctcaaatcttacagctgctcactcccctgcagggcaacgcccagggaccaagttagccccttaagcctaggcaaaagaatcccgcccataatcgagaagcgactcgacatggaggcgatgacgagatcacgcgaggaggaaaggagggagggcttcttccaggcccagggcggtccttacaagacgggaggcagcagagaactcccataaaggtattgcggcactcccctccccctgcccagaagggtgcggccttctctccacctcctccac \
        --guide_seq aatcgagaagcgactcgaca,taaggggctaacttggtccc
    """

    args = cr.parse_args(sys.argv[1:])
    args.fastq_r1 = "test_data/test_L001_R1_001.fastq.gz"
    args.fastq_r2 = "test_data/test_L001_R2_001.fastq.gz"
    args.amplicon_seq = (
        "gtcgcccctcaaatcttacagctgctcactc"
        "ccctgcagggcaacgcccagggaccaagttag"
        "ccccttaagcctaggcaaaagaatcccgccca"
        "taatcgagaagcgactcgacatggaggcgatg"
        "acgagatcacgcgaggaggaaaggagggaggg"
        "cttcttccaggcccagggcggtccttacaaga"
        "cgggaggcagcagagaactcccataaaggtat"
        "tgcggcactcccctccccctgcccagaagggt"
        "gcggccttctctccacctcctccac"
    )
    args.guide_seq = "aatcgagaagcgactcgaca,taaggggctaacttggtccc"
    args.n_processes = p
    args.keep_intermediate = keep_intermediate
    args.output_folder = f"PYTEST_RESULTS_FOLDER/pytest_p{p}"
    args.trim_sequences = False

    print(args)

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


@pytest.mark.parametrize(
    "p,keep_intermediate",
    [(5, False)],
)
def test1_run_crispresso(p, keep_intermediate):
    """Run end-to-end command
    
    Ground truth values are from the original CRIPResso Docker

    CRISPResso --output_folder TEST1_OUTPUT_FOLDER \
        -r1 test_data/test1_L001_R1_001.fastq.gz \
        -r2 test_data/test1_L001_R2_001.fastq.gz \
        --amplicon_seq gtcgcccctcaaatcttacagctgctcactcccctgcagggcaacgcccagggaccaagttagccccttaagcctaggcaaaagaatcccgcccataatcgagaagcgactcgacatggaggcgatgacgagatcacgcgaggaggaaaggagggagggcttcttccaggcccagggcggtccttacaagacgggaggcagcagagaactcccataaaggtattgcggcactcccctccccctgcccagaagggtgcggccttctctccacctcctccac \
        --guide_seq cgagaagcgactcgacatgg,aaggggctaacttggtccct \
        --min_identity_score 30.0  \
        --window_around_sgrna 23 \
        --trim_sequences


        --trimmomatic_options_string "ILLUMINACLIP:NexteraPE-PE.fa:0:90:10:0:true MINLEN:40"

    """

    args = cr.parse_args(sys.argv[1:])
    args.fastq_r1 = "test_data/test1_L001_R1_001.fastq.gz"
    args.fastq_r2 = "test_data/test1_L001_R2_001.fastq.gz"
    args.amplicon_seq = (
        "gtcgcccctcaaatcttacagctgctcactcccctgcagg"
        "gcaacgcccagggaccaagttagccccttaagcctaggcaa"
        "aagaatcccgcccataatcgagaagcgactcgacatggagg"
        "cgatgacgagatcacgcgaggaggaaaggagggagggcttc"
        "ttccaggcccagggcggtccttacaagacgggaggcagcag"
        "agaactcccataaaggtattgcggcactcccctccccctgc"
        "ccagaagggtgcggccttctctccacctcctccac"
    )
    args.guide_seq = "cgagaagcgactcgacatgg,aaggggctaacttggtccct"
    args.n_processes = p
    args.keep_intermediate = keep_intermediate
    args.output_folder = f"PYTEST_RESULTS_FOLDER/pytest1_p{p}"
    args.window_around_sgrna = 23
    args.min_identity_score = 30.0
    args.trim_sequences = True

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

    assert n_total == 4039
    assert n_reads_input == 4941
    assert n_unmodified == 2647
    assert n_mixed_hdr_nhej == 0
    assert n_modified == 1392
    assert n_repaired == 0
    assert nhej_inserted == 49
    assert nhej_deleted == 680
    assert nhej_mutated == 890

    assert tuple(df_indels["fq"].values[:4]) == (2, 4, 5, 5)
    assert tuple(df_insertion["fq"].values[:4]) == (3990, 6, 1, 0)
    assert tuple(df_deletion["fq"].values[:4]) == (3359, 43, 3, 0)
    assert tuple(df_substitution["fq"].values[:4]) == (3149, 693, 105, 23)
    assert tuple(df_alleles["#Reads"].values[:4]) == (184, 68, 44, 26)
