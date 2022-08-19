#!/usr/bin/env python
# -*- coding: utf8 -*-

# Run unit tests on CRISPResso using test data

# coverage run -m pytest -vvv crispresso_tests.py; coverage report -m

# CRISPResso --output_folder TEST_OUTPUT_FOLDER  -r1 test_data/test_L001_R1_001.fastq.gz -r2 test_data/test_L001_R2_001.fastq.gz --amplicon_seq gtcgcccctcaaatcttacagctgctcactcccctgcagggcaacgcccagggaccaagttagccccttaagcctaggcaaaagaatcccgcccataatcgagaagcgactcgacatggaggcgatgacgagatcacgcgaggaggaaaggagggagggcttcttccaggcccagggcggtccttacaagacgggaggcagcagagaactcccataaaggtattgcggcactcccctccccctgcccagaagggtgcggccttctctccacctcctccac --guide_seq aatcgagaagcgactcgaca,taaggggctaacttggtccc

import os
import numpy as np

import pytest

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

    assert  np.array(cr.find_wrong_nt("ACBTGCNGRCCACTGFNNC")).sort() == np.array(["R", "F", "B"]).sort()

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
