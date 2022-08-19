#!/usr/bin/env python
# -*- coding: utf8 -*-

# Run unit tests on CRISPResso using test data

# coverage run -m pytest -vvv crispresso_tests.py; coverage report -m

# CRISPResso --output_folder TEST_OUTPUT_FOLDER  -r1 test_data/test_L001_R1_001.fastq.gz -r2 test_data/test_L001_R2_001.fastq.gz --amplicon_seq gtcgcccctcaaatcttacagctgctcactcccctgcagggcaacgcccagggaccaagttagccccttaagcctaggcaaaagaatcccgcccataatcgagaagcgactcgacatggaggcgatgacgagatcacgcgaggaggaaaggagggagggcttcttccaggcccagggcggtccttacaagacgggaggcagcagagaactcccataaaggtattgcggcactcccctccccctgcccagaagggtgcggccttctctccacctcctccac --guide_seq aatcgagaagcgactcgaca,taaggggctaacttggtccc

import CRISPResso as cr


def test_count_reads():

    assert cr.get_n_reads_fastq("test_data/test_L001_R1_001.fastq.gz") == 8906
    assert cr.get_n_reads_fastq("test_data/test_L001_R2_001.fastq.gz") == 8906


def test_check_library():

    assert cr.check_library("pandas")
    assert cr.check_library("numpy")


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

    assert (
        cr.filter_pe_fastq_by_qual(
            "test_data/test_L001_R1_001.fastq.gz", "test_data/test_L001_R2_001.fastq.gz"
        )
        == "test_data/test_L001_R1_001_filtered.fastq.gz",
        "test_data/test_L001_R2_001_filtered.fastq.gz",
    )
