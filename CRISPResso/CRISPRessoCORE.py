#!/usr/bin/env python
# -*- coding: utf8 -*-

"""
CRISPResso - Luca Pinello 2015
Software pipeline for the analysis of CRISPR-Cas9 genome
editing outcomes from deep sequencing data
https://github.com/lucapinello/CRISPResso
"""
__version__ = "1.1.0"

from typing import List
import sys
import errno
import os
from shutil import which
import subprocess as sb
import argparse
import re
import gzip
from collections import defaultdict
import multiprocessing as mp
import pickle
import unicodedata
import traceback

import datetime

import logging

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager as fm
from matplotlib import colors as colors_mpl
from matplotlib import gridspec

import plotly.express as px

from Bio import SeqIO, pairwise2

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n",
    datefmt="%a, %d %b %Y %H:%M:%S",
    stream=sys.stderr,
    filemode="w",
)
error = logging.critical
warning = logging.warning
debug = logging.debug
info = logging.info

# ###Support functions###


def get_data(path: str) -> str:
    """Combine the data path with the absolute path

    Args:
        path(str): Relative path
    Returns:
        Global absolute path
    """
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", path)


def check_library(library_name: str):
    """Checks if library can be loaded

    Args:
        Library name to check
    """
    try:
        return __import__(library_name)
    except Exception as exc:
        raise Exception(
            f"You need to install {library_name} to use CRISPResso!"
        ) from exc


def check_program(binary_name: str, download_url: str = None):
    """Checks if program exists

    Args:
        binary_name(str): Name of the program to check
        download_url(str): URL for download if it isn't installed
    """
    if not which(binary_name):
        if download_url:
            error(f"You can download it from here: {download_url}")
        raise Exception(
            "You need to install and have the command"
            f" #####{binary_name}##### in your PATH "
            "variable to use CRISPResso!\n Please read the documentation!"
        )

    return True


def check_file(filename: str):
    """Checks if file can be opened

    Args:
        filename(str): Filename to check
    """
    try:
        with open(filename, "r", encoding="utf-8"):
            pass
    except IOError as exc:
        raise Exception(f"I cannot open the file: {filename}") from exc


def force_symlink(src: str, dst: str) -> None:
    """Symbolic link two directories

    Args:
        src(str): Source directory
        dst(str): Destination directory
    """

    if os.path.exists(dst) and os.path.samefile(src, dst):
        return

    try:
        os.symlink(src, dst)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src, dst)


def reverse_complement(sequence: str) -> str:
    """Reverse complement of sequence

    Calculate the reverse complement of the nucleotide sequence.

    Args:
        sequence(str): Input sequence of nucleotides

    Returns:
        Reverse complement string
    """

    nt_complement = dict(
        {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "_": "_", "-": "-"}
    )
    return "".join([nt_complement[c] for c in sequence.upper()[-1::-1]])


def find_wrong_nt(sequence: str) -> List[str]:
    """Finds any incorrect base pair labels

    Makes sure the nucleotide sequence is valid.

    Args:
        sequence(str): The input nucleotide sequence

    Returns:
        Any invalid nucleotide labels

    """
    return list(set(sequence.upper()).difference(set(["A", "T", "C", "G", "N"])))


def get_ids_reads_to_remove(
    fastq_filename: str, min_bp_quality: int = 20, min_single_bp_quality: int = 0
) -> List[str]:
    """Get list of reads to remove

    Gets a list of the reads that are below the basepair quality specified

    Args:
        fastq_filename(str): Filename of fastq file
        min_bp_quality(int): Minimum average basepair quality score
        min_single_bp_quality(int): Minimum single basepair quality score

    Returns:
        List of read labels to remove due to low quality scores

    """

    ids_to_remove = set()
    if fastq_filename.endswith(".gz"):
        fastq_handle = gzip.open(fastq_filename, "rt")
    else:
        fastq_handle = open(fastq_filename, "r", encoding="utf-8")

    for record in SeqIO.parse(fastq_handle, "fastq"):
        if (
            np.array(record.letter_annotations["phred_quality"]).mean() < min_bp_quality
            or np.array(record.letter_annotations["phred_quality"]).min()
            < min_single_bp_quality
        ):
            ids_to_remove.add(record.id)

    return ids_to_remove


def filter_pe_fastq_by_qual(
    fastq_r1: str,
    fastq_r2: str,
    output_filename_r1: str = None,
    output_filename_r2: str = None,
    min_bp_quality: int = 20,
    min_single_bp_quality: int = 0,
):
    """Filters paired end reads by quality score.

    Args:
        fastq_r1(str): Fastq filename for read 1
        fastq_r2(str): Fastq filename for read 2
        output_filename_r1(str): Output filename for read 1
        output_filename_r2(str): Output filename for read 2
        min_bp_quality(int): Minimum average base quality
        min_single_bp_quality(int): Minimum single base quality

    """

    ids_to_remove_s1 = get_ids_reads_to_remove(
        fastq_r1,
        min_bp_quality=min_bp_quality,
        min_single_bp_quality=min_single_bp_quality,
    )
    ids_to_remove_s2 = get_ids_reads_to_remove(
        fastq_r2,
        min_bp_quality=min_bp_quality,
        min_single_bp_quality=min_single_bp_quality,
    )

    ids_to_remove = ids_to_remove_s1.union(ids_to_remove_s2)

    if fastq_r1.endswith(".gz"):
        fastq_handle_r1 = gzip.open(fastq_r1, "rt")
    else:
        fastq_handle_r1 = open(fastq_r1, "r", encoding="utf-8")

    if fastq_r2.endswith(".gz"):
        fastq_handle_r2 = gzip.open(fastq_r2, "rt")
    else:
        fastq_handle_r2 = open(fastq_r2, "r", encoding="utf-8")

    if not output_filename_r1:
        output_filename_r1 = (
            fastq_r1.replace(".fastq", "").replace(".gz", "") + "_filtered.fastq.gz"
        )

    if not output_filename_r2:
        output_filename_r2 = (
            fastq_r2.replace(".fastq", "").replace(".gz", "") + "_filtered.fastq.gz"
        )

    try:
        fastq_filtered_outfile_r1 = gzip.open(output_filename_r1, "wt")

        for record in SeqIO.parse(fastq_handle_r1, "fastq"):
            if record.id not in ids_to_remove:
                fastq_filtered_outfile_r1.write(record.format("fastq"))
    except Exception as exc:
        raise Exception("Error handling the fastq_filtered_outfile_r1") from exc

    try:
        fastq_filtered_outfile_r2 = gzip.open(output_filename_r2, "wt")

        for record in SeqIO.parse(fastq_handle_r2, "fastq"):
            if record.id not in ids_to_remove:
                fastq_filtered_outfile_r2.write(record.format("fastq"))
    except Exception as exc:
        raise Exception("Error handling the fastq_filtered_outfile_r2") from exc

    return output_filename_r1, output_filename_r2


def filter_se_fastq_by_qual(
    fastq_filename: str,
    output_filename=None,
    min_bp_quality: int = 20,
    min_single_bp_quality: int = 0,
):
    """Filters single end reads by quality score.

    Args:
        fastq_filename(str): Fastq filename for read
        output_filename(str): Output filename for read
        min_bp_quality(int): Minimum average base quality
        min_single_bp_quality(int): Minimum single base quality

    """
    if fastq_filename.endswith(".gz"):
        fastq_handle = gzip.open(fastq_filename, "rt")
    else:
        fastq_handle = open(fastq_filename, "r", encoding="utf-8")

    if not output_filename:
        output_filename = (
            fastq_filename.replace(".fastq", "").replace(".gz", "")
            + "_filtered.fastq.gz"
        )

    try:
        fastq_filtered_outfile = gzip.open(output_filename, "wt")

        for record in SeqIO.parse(fastq_handle, "fastq"):
            if (
                np.array(record.letter_annotations["phred_quality"]).mean()
                >= min_bp_quality
                and np.array(record.letter_annotations["phred_quality"]).min()
                >= min_single_bp_quality
            ):
                fastq_filtered_outfile.write(record.format("fastq"))
    except Exception as exc:
        raise Exception("Error handling the fastq_filtered_outfile") from exc

    return output_filename


def get_average_read_length_fastq(fastq_filename: str) -> int:
    """Get average read length for FastQ

    Args:
        fastq_filename(str): FastQ filename
    Returns:
        Average read length of FastQ

    """
    cmd = (
        ("z" if fastq_filename.endswith(".gz") else "")
        + (f"cat < {fastq_filename}")
        + r""" | awk 'BN {n=0;s=0;} NR%4 == 2 {s+=length($0);n++;} END { printf("%d\n",s/n)}' """
    )
    p = sb.Popen(cmd, shell=True, stdout=sb.PIPE)
    return int(p.communicate()[0].strip())


def get_n_reads_fastq(fastq_filename: str) -> int:
    """Get number of reads in FastQ

    Args:
        fastq_filename(str): FastQ filename
    Returns:
        Average read length of FastQ

    """
    p = sb.Popen(
        ("z" if fastq_filename.endswith(".gz") else "")
        + f"cat < {fastq_filename} | wc -l",
        shell=True,
        stdout=sb.PIPE,
    )
    return int(float(p.communicate()[0]) // 4)


matplotlib = check_library("matplotlib")

font = {"size": 22}
matplotlib.rc("font", **font)
matplotlib.use("Agg")

plt = check_library("pylab")

pd = check_library("pandas")
np = check_library("numpy")
Bio = check_library("Bio")

check_program("java")
check_program("flash")
check_program("needle")

sns = check_library("seaborn")
sns.set_context("poster")
sns.set(font_scale=2.2)
sns.set_style("white")

# ##EXCEPTIONS############################


class FlashException(Exception):  # pragma: no cover
    pass


class TrimmomaticException(Exception):  # pragma: no cover
    pass


class NeedleException(Exception):  # pragma: no cover
    pass


class NoReadsAlignedException(Exception):  # pragma: no cover
    pass


class DonorSequenceException(Exception):  # pragma: no cover
    pass


class AmpliconEqualDonorException(Exception):  # pragma: no cover
    pass


class CoreDonorSequenceNotContainedException(Exception):  # pragma: no cover
    pass


class CoreDonorSequenceNotUniqueException(Exception):  # pragma: no cover
    pass


class SgRNASequenceException(Exception):  # pragma: no cover
    pass


class NTException(Exception):  # pragma: no cover
    pass


class ExonSequenceException(Exception):  # pragma: no cover
    pass


class DuplicateSequenceIdException(Exception):  # pragma: no cover
    pass


class NoReadsAfterQualityFiltering(Exception):  # pragma: no cover
    pass


# ########################################


def process_df_chunk(chunk_input):
    """Process one chunk of the dataframe

    Using multiple processes and imap to
    parallelize the dataframe processing.
    (Not sure if this is necessary)

    """

    df_needle_alignment_chunk = chunk_input[0]
    args = chunk_input[1]

    modified_frameshift = 0
    modified_non_frameshift = 0
    non_modified_non_frameshift = 0
    splicing_sites_modified = 0

    # INITIALIZATIONS
    if args.coding_seq:
        perform_frameshift_analysis = True
    else:
        perform_frameshift_analysis = False

    effect_vector_insertion = np.zeros(LEN_AMPLICON)
    effect_vector_deletion = np.zeros(LEN_AMPLICON)
    effect_vector_mutation = np.zeros(LEN_AMPLICON)
    effect_vector_any = np.zeros(LEN_AMPLICON)

    effect_vector_insertion_mixed = np.zeros(LEN_AMPLICON)
    effect_vector_deletion_mixed = np.zeros(LEN_AMPLICON)
    effect_vector_mutation_mixed = np.zeros(LEN_AMPLICON)

    effect_vector_insertion_hdr = np.zeros(LEN_AMPLICON)
    effect_vector_deletion_hdr = np.zeros(LEN_AMPLICON)
    effect_vector_mutation_hdr = np.zeros(LEN_AMPLICON)

    effect_vector_insertion_noncoding = np.zeros(LEN_AMPLICON)
    effect_vector_deletion_noncoding = np.zeros(LEN_AMPLICON)
    effect_vector_mutation_noncoding = np.zeros(LEN_AMPLICON)

    hist_inframe = defaultdict(lambda: 0)
    hist_frameshift = defaultdict(lambda: 0)

    avg_vector_del_all = np.zeros(LEN_AMPLICON)
    avg_vector_ins_all = np.zeros(LEN_AMPLICON)

    re_find_indels = re.compile("(-*-)")
    re_find_substitutions = re.compile(r"(\.*\.)")

    for idx_row, row in df_needle_alignment_chunk.iterrows():

        # GET THE MUTATIONS POSITIONS
        if row.UNMODIFIED:
            continue

        if perform_frameshift_analysis:
            lenght_modified_positions_exons = []
            current_read_exons_modified = False
            current_read_spliced_modified = False

        # quantify substitution
        substitution_positions = []
        if not args.ignore_substitutions:
            for p in re_find_substitutions.finditer(row.align_str):
                st, en = p.span()
                substitution_positions.append(row.ref_positions[st:en])

            if substitution_positions:
                substitution_positions = list(np.hstack(substitution_positions))

        # quantify deletion
        deletion_positions = []
        deletion_positions_flat = []
        deletion_sizes = []

        if not args.ignore_deletions:
            for p in re_find_indels.finditer(row.align_seq):
                st, en = p.span()
                deletion_positions.append(row.ref_positions[st:en])
                deletion_sizes.append(en - st)

            if deletion_positions:
                deletion_positions_flat = np.hstack(deletion_positions)

        # quantify insertion
        insertion_positions = []
        insertion_sizes = []
        insertion_positions_flat = []

        if not args.ignore_insertions:
            for p in re_find_indels.finditer(row.ref_seq):
                st, en = p.span()
                # ref_st=row.ref_positions[st-1]
                # # we report the base preceding the insertion

                # insertion_positions.append(ref_st)
                insertion_positions.append(
                    [
                        row["ref_positions"][max(0, st - 1)],
                        row["ref_positions"][min(len(row["ref_positions"]) - 1, en)],
                    ]
                )
                insertion_sizes.append(en - st)

            if insertion_positions:
                insertion_positions_flat = np.hstack(insertion_positions)

        # #######CLASSIFY READ
        # WE HAVE THE DONOR SEQUENCE
        if args.expected_hdr_amplicon_seq:

            # HDR
            if (row.score_diff < 0) & (
                row.score_repaired >= args.hdr_perfect_alignment_threshold
            ):
                df_needle_alignment_chunk.loc[idx_row, "HDR"] = True

            # MIXED
            elif (row.score_diff < 0) & (
                row.score_repaired < args.hdr_perfect_alignment_threshold
            ):
                df_needle_alignment_chunk.loc[idx_row, "MIXED"] = True

            else:
                # NHEJ
                if (
                    INCLUDE_IDXS.intersection(substitution_positions)
                    or INCLUDE_IDXS.intersection(insertion_positions_flat)
                    or INCLUDE_IDXS.intersection(deletion_positions_flat)
                ):
                    df_needle_alignment_chunk.loc[idx_row, "NHEJ"] = True

                # UNMODIFIED
                else:
                    df_needle_alignment_chunk.loc[idx_row, "UNMODIFIED"] = True

        # NO DONOR SEQUENCE PROVIDED
        else:
            # NHEJ
            if (
                INCLUDE_IDXS.intersection(substitution_positions)
                or INCLUDE_IDXS.intersection(insertion_positions_flat)
                or INCLUDE_IDXS.intersection(deletion_positions_flat)
            ):
                df_needle_alignment_chunk.loc[idx_row, "NHEJ"] = True

            # UNMODIFIED
            else:
                df_needle_alignment_chunk.loc[idx_row, "UNMODIFIED"] = True

        # ##CREATE AVERAGE SIGNALS, HERE WE SHOW EVERYTHING...
        if df_needle_alignment_chunk.loc[idx_row, "MIXED"]:
            effect_vector_mutation_mixed[substitution_positions] += 1
            effect_vector_deletion_mixed[deletion_positions_flat] += 1
            effect_vector_insertion_mixed[insertion_positions_flat] += 1

        elif df_needle_alignment_chunk.loc[idx_row, "HDR"]:
            effect_vector_mutation_hdr[substitution_positions] += 1
            effect_vector_deletion_hdr[deletion_positions_flat] += 1
            effect_vector_insertion_hdr[insertion_positions_flat] += 1

        elif (
            df_needle_alignment_chunk.loc[idx_row, "NHEJ"]
            and not args.hide_mutations_outside_window_NHEJ
        ):
            effect_vector_mutation[substitution_positions] += 1
            effect_vector_deletion[deletion_positions_flat] += 1
            effect_vector_insertion[insertion_positions_flat] += 1

        any_positions = np.unique(
            np.hstack(
                [
                    deletion_positions_flat,
                    insertion_positions_flat,
                    substitution_positions,
                ]
            )
        ).astype(int)
        effect_vector_any[any_positions] += 1

        # For NHEJ we count only the events that overlap
        # the window specified around
        # the cut site (1bp by default)...
        if df_needle_alignment_chunk.loc[idx_row, "NHEJ"] and args.window_around_sgrna:

            substitution_positions = list(
                INCLUDE_IDXS.intersection(substitution_positions)
            )

            insertion_positions_window = []
            insertion_sizes_window = []

            # count insertions overlapping
            for idx_ins, ins_pos_set in enumerate(insertion_positions):
                # print ref_st, insertion_positions
                if INCLUDE_IDXS.intersection(ins_pos_set):
                    insertion_positions_window.append(ins_pos_set)
                    insertion_sizes_window.append(insertion_sizes[idx_ins])

            insertion_positions = insertion_positions_window
            insertion_sizes = insertion_sizes_window

            deletion_positions_window = []
            deletion_sizes_window = []
            for idx_del, del_pos_set in enumerate(deletion_positions):
                if INCLUDE_IDXS.intersection(del_pos_set):
                    deletion_positions_window.append(del_pos_set)
                    deletion_sizes_window.append(deletion_sizes[idx_del])

            deletion_positions = deletion_positions_window
            deletion_sizes = deletion_sizes_window

            if deletion_positions:
                deletion_positions_flat = np.hstack(deletion_positions)

        if (
            df_needle_alignment_chunk.loc[idx_row, "NHEJ"]
            and args.hide_mutations_outside_window_NHEJ
        ):
            effect_vector_mutation[substitution_positions] += 1
            effect_vector_deletion[deletion_positions_flat] += 1
            effect_vector_insertion[insertion_positions_flat] += 1

        # ###QUANTIFICATION AND FRAMESHIFT ANALYSIS
        if not df_needle_alignment_chunk.loc[idx_row, "UNMODIFIED"]:

            df_needle_alignment_chunk.loc[idx_row, "n_mutated"] = len(
                substitution_positions
            )
            df_needle_alignment_chunk.loc[idx_row, "n_inserted"] = np.sum(
                insertion_sizes
            )
            df_needle_alignment_chunk.loc[idx_row, "n_deleted"] = np.sum(deletion_sizes)

            for idx_ins, ins_pos_set in enumerate(insertion_positions):
                avg_vector_ins_all[ins_pos_set] += insertion_sizes[idx_ins]

                if perform_frameshift_analysis:
                    if set(EXON_POSITIONS).intersection(
                        ins_pos_set
                    ):  # check that we are inserting in one exon
                        lenght_modified_positions_exons.append(insertion_sizes[idx_ins])
                        current_read_exons_modified = True

            for idx_del, del_pos_set in enumerate(deletion_positions):
                avg_vector_del_all[del_pos_set] += deletion_sizes[idx_del]

            if perform_frameshift_analysis:
                del_positions_to_append = sorted(
                    set(EXON_POSITIONS).intersection(set(deletion_positions_flat))
                )
                if del_positions_to_append:
                    # Always use the low include upper not
                    current_read_exons_modified = True
                    lenght_modified_positions_exons.append(
                        -len(del_positions_to_append)
                    )

                if set(EXON_POSITIONS).intersection(substitution_positions):
                    current_read_exons_modified = True

                if set(SPLICING_POSITIONS).intersection(substitution_positions):
                    current_read_spliced_modified = True

                if set(SPLICING_POSITIONS).intersection(deletion_positions_flat):
                    current_read_spliced_modified = True

                if set(SPLICING_POSITIONS).intersection(insertion_positions_flat):
                    current_read_spliced_modified = True

                if current_read_spliced_modified:
                    splicing_sites_modified += 1

                # if modified check if frameshift
                if current_read_exons_modified:

                    if not lenght_modified_positions_exons:
                        # there are no indels
                        modified_non_frameshift += 1
                        hist_inframe[0] += 1
                    else:

                        effetive_length = sum(lenght_modified_positions_exons)

                        if (effetive_length % 3) == 0:
                            modified_non_frameshift += 1
                            hist_inframe[effetive_length] += 1
                        else:
                            modified_frameshift += 1
                            hist_frameshift[effetive_length] += 1

                # the indels and subtitutions are outside the exon/s
                # so we don't care!
                else:
                    non_modified_non_frameshift += 1
                    effect_vector_insertion_noncoding[insertion_positions_flat] += 1
                    effect_vector_deletion_noncoding[deletion_positions_flat] += 1
                    effect_vector_mutation_noncoding[substitution_positions] += 1

    hist_inframe = dict(hist_inframe)
    hist_frameshift = dict(hist_frameshift)

    return (
        df_needle_alignment_chunk,
        effect_vector_insertion,
        effect_vector_deletion,
        effect_vector_mutation,
        effect_vector_any,
        effect_vector_insertion_mixed,
        effect_vector_deletion_mixed,
        effect_vector_mutation_mixed,
        effect_vector_insertion_hdr,
        effect_vector_deletion_hdr,
        effect_vector_mutation_hdr,
        effect_vector_insertion_noncoding,
        effect_vector_deletion_noncoding,
        effect_vector_mutation_noncoding,
        hist_inframe,
        hist_frameshift,
        avg_vector_del_all,
        avg_vector_ins_all,
        modified_frameshift,
        modified_non_frameshift,
        non_modified_non_frameshift,
        splicing_sites_modified,
    )


def add_hist(hist_to_add, hist_global):
    for key, value in hist_to_add.items():
        hist_global[key] += value
    return hist_global


def slugify(value, allow_unicode=False):  # adapted from the Django project
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize("NFKC", value)
    else:
        value = (
            unicodedata.normalize("NFKD", value)
            .encode("ascii", "ignore")
            .decode("ascii")
        )
    value = re.sub(r"[^\w\s-]", "", value.lower())
    return re.sub(r"[-\s]+", "-", value).strip("-_")


def split_paired_end_reads_single_file(
    fastq_filename, output_filename_r1, output_filename_r2
):

    if fastq_filename.endswith(".gz"):
        fastq_handle = gzip.open(fastq_filename)
    else:
        fastq_handle = open(fastq_filename, "r", encoding="utf-8")

    try:
        fastq_splitted_outfile_r1 = gzip.open(output_filename_r1, "wt")
        fastq_splitted_outfile_r2 = gzip.open(output_filename_r2, "wt")
    except:
        raise Exception("Error handling the splitting operation")

    return output_filename_r1, output_filename_r2


def get_row_around_cut(row, cut_point, offset):
    cut_idx = row["ref_positions"].index(cut_point)
    return (
        row["Aligned_Sequence"][cut_idx - offset + 1 : cut_idx + offset + 1],
        row["Reference_Sequence"][cut_idx - offset + 1 : cut_idx + offset + 1],
        row["UNMODIFIED"],
        row["%Reads"],
        row["#Reads"],
    )


def get_dataframe_around_cut(df_alleles, cut_point, offset):
    df_alleles_around_cut = pd.DataFrame(
        list(
            df_alleles.apply(
                lambda row: get_row_around_cut(row, cut_point, offset), axis=1
            ).values
        ),
        columns=[
            "Aligned_Sequence",
            "Reference_Sequence",
            "Unedited",
            "%Reads",
            "#Reads",
        ],
    )
    df_alleles_around_cut = (
        df_alleles_around_cut.groupby(["Aligned_Sequence", "Reference_Sequence"])
        .sum()
        .reset_index()
        .set_index("Aligned_Sequence")
    )

    df_alleles_around_cut.sort_values(by="%Reads", inplace=True, ascending=False)
    df_alleles_around_cut["Unedited"] = df_alleles_around_cut["Unedited"] > 0
    return df_alleles_around_cut


# We need to customize the seaborn heatmap class and function
class Custom_HeatMapper(sns.matrix._HeatMapper):
    def __init__(
        self,
        data,
        vmin,
        vmax,
        cmap,
        center,
        robust,
        annot,
        fmt,
        annot_kws,
        per_element_annot_kws,
        cbar,
        cbar_kws,
        xticklabels=True,
        yticklabels=True,
        mask=None,
    ):

        super(Custom_HeatMapper, self).__init__(
            data,
            vmin,
            vmax,
            cmap,
            center,
            robust,
            annot,
            fmt,
            annot_kws,
            cbar,
            cbar_kws,
            xticklabels,
            yticklabels,
            mask,
        )

        if annot is not None:
            if per_element_annot_kws is None:
                self.per_element_annot_kws = np.empty_like(annot, dtype=object)
                self.per_element_annot_kws[:] = dict()
            else:
                self.per_element_annot_kws = per_element_annot_kws

    # add per element dict to syle the annotatiin
    def _annotate_heatmap(self, ax, mesh):
        """Add textual labels with the value in each cell."""
        mesh.update_scalarmappable()
        xpos, ypos = np.meshgrid(ax.get_xticks(), ax.get_yticks())

        for x, y, m, color, val, per_element_dict in zip(
            xpos.flat,
            ypos.flat,
            mesh.get_array(),
            mesh.get_facecolors(),
            self.annot_data.flat,
            self.per_element_annot_kws.flat,
        ):
            # print per_element_dict
            if m is not np.ma.masked:
                l = sns.utils.relative_luminance(color)
                text_color = ".15" if l > 0.408 else "w"
                annotation = ("{:" + self.fmt + "}").format(val)
                text_kwargs = dict(color=text_color, ha="center", va="center")
                text_kwargs.update(self.annot_kws)
                text_kwargs.update(per_element_dict)

                ax.text(x, y, annotation, **text_kwargs)

    # removed the colobar
    def plot(self, ax, kws):
        """Draw the heatmap on the provided Axes."""
        # Remove all the Axes spines
        sns.utils.despine(ax=ax, left=True, bottom=True)

        # Draw the heatmap
        mesh = ax.pcolormesh(
            self.plot_data, vmin=self.vmin, vmax=self.vmax, cmap=self.cmap, **kws
        )

        # Set the axis limits
        ax.set(xlim=(0, self.data.shape[1]), ylim=(0, self.data.shape[0]))

        # Add row and column labels
        ax.set(xticks=self.xticks, yticks=self.yticks)
        xtl = ax.set_xticklabels(self.xticklabels)
        ytl = ax.set_yticklabels(self.yticklabels, rotation="vertical")

        # Possibly rotate them if they overlap
        plt.draw()
        if sns.utils.axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation="vertical")
        if sns.utils.axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")

        # Add the axis labels
        ax.set(xlabel=self.xlabel, ylabel=self.ylabel)

        # Annotate the cells with the formatted values
        if self.annot:
            self._annotate_heatmap(ax, mesh)


def custom_heatmap(
    data,
    vmin=None,
    vmax=None,
    cmap=None,
    center=None,
    robust=False,
    annot=None,
    fmt=".2g",
    annot_kws=None,
    per_element_annot_kws=None,
    linewidths=0,
    linecolor="white",
    cbar=True,
    cbar_kws=None,
    square=False,
    ax=None,
    xticklabels=True,
    yticklabels=True,
    mask=None,
    **kwargs,
):

    # Initialize the plotter object
    plotter = Custom_HeatMapper(
        data,
        vmin,
        vmax,
        cmap,
        center,
        robust,
        annot,
        fmt,
        annot_kws,
        per_element_annot_kws,
        cbar,
        cbar_kws,
        xticklabels,
        yticklabels,
        mask,
    )

    # Add the pcolormesh kwargs here
    kwargs["linewidths"] = linewidths
    kwargs["edgecolor"] = linecolor

    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()
    if square:
        ax.set_aspect("equal")
    plotter.plot(ax, kwargs)
    return ax


def plot_alleles_table(
    args,
    cut_point,
    df_alleles,
    sg_rna_name,
    output_directory,
    min_frequency: float = 0.5,
    max_n_rows: int = 100,
) -> pd.DataFrame:
    reference_seq = args.amplicon_seq
    min_frequency = args.min_frequency_alleles_around_cut_to_plot
    max_n_rows = args.max_rows_alleles_around_cut_to_plot

    # bp we are plotting on each side
    offset_around_cut_to_plot = len(df_alleles.index[0]) // 2

    # make a color map of fixed colors
    alpha = 0.5

    get_color = lambda x, y, z: (x / 255.0, y / 255.0, z / 255.0, alpha)
    a_color = get_color(127, 201, 127)
    t_color = get_color(190, 174, 212)
    c_color = get_color(253, 192, 134)
    g_color = get_color(255, 255, 153)
    n_color = get_color(255, 255, 255)
    indel_color = get_color(230, 230, 230)

    cmap = colors_mpl.ListedColormap(
        [indel_color, a_color, t_color, c_color, g_color, n_color]
    )

    dna_to_numbers = {"-": 0, "A": 1, "T": 2, "C": 3, "G": 4, "N": 5}
    seq_to_numbers = lambda seq: [dna_to_numbers[x] for x in seq]

    x_labels = []
    annot = []
    y_labels = []
    lines = defaultdict(list)

    re_find_indels = re.compile("(-*-)")

    per_element_annot_kws = []
    idx_row = 0

    for idx, row in df_alleles.loc[df_alleles["%Reads"] >= min_frequency][
        :max_n_rows
    ].iterrows():
        x_labels.append(seq_to_numbers(str.upper(idx)))
        annot.append(list(idx))
        y_labels.append(f"{row['%Reads']}% ({row['#Reads']} reads)")

        for p in re_find_indels.finditer(row["Reference_Sequence"]):
            lines[idx_row].append((p.start(), p.end()))

        lines[idx_row].append(row["Reference_Sequence"])
        idx_row += 1

        idxs_sub = [
            i_sub
            for i_sub in range(len(idx))
            if (row["Reference_Sequence"][i_sub] != idx[i_sub])
            and (row["Reference_Sequence"][i_sub] != "-")
            and (idx[i_sub] != "-")
        ]
        to_append = np.array([{}] * len(idx), dtype=object)
        to_append[idxs_sub] = {"weight": "bold", "color": "black", "size": 16}
        per_element_annot_kws.append(to_append)

    ref_seq_around_cut = reference_seq[
        (cut_point - offset_around_cut_to_plot + 1) : (
            cut_point + offset_around_cut_to_plot + 1
        )
    ]

    per_element_annot_kws = np.vstack(per_element_annot_kws[::-1])
    ref_seq_hm = np.expand_dims(seq_to_numbers(ref_seq_around_cut), 1).T
    ref_seq_annot_hm = np.expand_dims(list(ref_seq_around_cut), 1).T

    annot = annot[::-1]
    x_labels = x_labels[::-1]

    sns.set_context("poster")

    n_rows = len(x_labels)
    n_cols = offset_around_cut_to_plot * 2

    plt.figure(figsize=(offset_around_cut_to_plot * 0.6, (n_rows + 1) * 0.6))
    gs1 = gridspec.GridSpec(n_rows + 1, n_cols)
    gs2 = gridspec.GridSpec(n_rows + 1, n_cols)

    ax_hm_ref = plt.subplot(gs1[0, :])
    ax_hm = plt.subplot(gs2[1:, :])

    custom_heatmap(
        ref_seq_hm,
        annot=ref_seq_annot_hm,
        annot_kws={"size": 16},
        cmap=cmap,
        fmt="s",
        ax=ax_hm_ref,
        vmin=0,
        vmax=5,
        square=True,
    )
    custom_heatmap(
        x_labels,
        annot=np.array(annot),
        annot_kws={"size": 16},
        cmap=cmap,
        fmt="s",
        ax=ax_hm,
        square=True,
        vmin=0,
        vmax=5,
        per_element_annot_kws=per_element_annot_kws,
    )

    ax_hm.yaxis.tick_right()
    ax_hm.yaxis.set_ticklabels(y_labels[::-1], rotation=True)
    ax_hm.xaxis.set_ticks([])

    # print lines

    # cut point vertical line
    ax_hm.vlines([offset_around_cut_to_plot], *ax_hm.get_ylim(), linestyles="dashed")

    # create boxes for ins
    for idx, lss in lines.items():
        for ls in lss:
            for line in ls:
                ax_hm.vlines([line], n_rows - idx - 1, n_rows - idx, color="red", lw=3)

            ax_hm.hlines(n_rows - idx - 1, ls[0], ls[1], color="red", lw=3)
            ax_hm.hlines(n_rows - idx, ls[0], ls[1], color="red", lw=3)

    ax_hm_ref.yaxis.tick_right()
    ax_hm_ref.xaxis.set_ticks([])
    ax_hm_ref.yaxis.set_ticklabels(["Reference"], rotation=True)

    gs2.update(
        left=0, right=1, hspace=0.05, wspace=0, top=1 * (((n_rows) * 1.13)) / (n_rows)
    )
    gs1.update(
        left=0,
        right=1,
        hspace=0.05,
        wspace=0,
    )

    sns.set_context(
        rc={
            "lines.markeredgewidth": 1,
            "mathtext.fontset": "stix",
            "text.usetex": True,
            "text.latex.unicode": True,
        }
    )

    proxies = [
        matplotlib.lines.Line2D(
            [0],
            [0],
            linestyle="none",
            mfc="black",
            mec="none",
            marker=r"$\mathbf{{{}}}$".format("bold"),
            ms=18,
        ),
        matplotlib.lines.Line2D(
            [0],
            [0],
            linestyle="none",
            mfc="none",
            mec="red",
            marker="s",
            ms=8,
            markeredgewidth=2.5,
        ),
        matplotlib.lines.Line2D(
            [0],
            [0],
            linestyle="none",
            mfc="none",
            mec="black",
            marker="_",
            ms=2,
        ),
        matplotlib.lines.Line2D([0], [1], linestyle="--", c="black", ms=6),
    ]  #
    descriptions = [
        "Substitutions",
        "Insertions",
        "Deletions",
        "Predicted cleavage position",
    ]
    ax_hm_ref.legend(
        proxies,
        descriptions,
        numpoints=1,
        markerscale=2,
        loc="center",
        bbox_to_anchor=(0.5, 4),
        ncol=1,
    )

    _jp = lambda filename: os.path.join(output_directory, filename)

    plt.savefig(
        _jp(f"9.Alleles_around_cut_site_for_{sg_rna_name}.pdf"), bbox_inches="tight"
    )
    if args.save_also_png:
        plt.savefig(
            _jp(f"9.Alleles_around_cut_site_for_{sg_rna_name}.png"),
            bbox_inches="tight",
            pad_inches=1,
        )


def run_crispresso(args):
    """Run the CRISPResso script

    Args:
        args: The CRISPResso arguments

        usage: CRISPResso [-h]
        -r1 FASTQ_R1
        [-r2 FASTQ_R2]
        -a AMPLICON_SEQ
        [-g GUIDE_SEQ]
        [-e EXPECTED_HDR_AMPLICON_SEQ]
        [-d DONOR_SEQ]
        [-c CODING_SEQ]
        [-q MIN_AVERAGE_READ_QUALITY]
        [-s MIN_SINGLE_BP_QUALITY]
        [--min_identity_score MIN_IDENTITY_SCORE]
        [-n NAME]
        [-o OUTPUT_FOLDER]
        [--split_paired_end]
        [--trim_sequences]
        [--trimmomatic_options_string TRIMMOMATIC_OPTIONS_STRING]
        [--min_paired_end_reads_overlap MIN_PAIRED_END_READS_OVERLAP]
        [--max_paired_end_reads_overlap MAX_PAIRED_END_READS_OVERLAP]
        [--hide_mutations_outside_window_NHEJ]
        [-w WINDOW_AROUND_SGRNA]
        [--cleavage_offset CLEAVAGE_OFFSET]
        [--exclude_bp_from_left EXCLUDE_BP_FROM_LEFT]
        [--exclude_bp_from_right EXCLUDE_BP_FROM_RIGHT]
        [--hdr_perfect_alignment_threshold HDR_PERFECT_ALIGNMENT_THRESHOLD]
        [--ignore_substitutions]
        [--ignore_insertions]
        [--ignore_deletions]
        [--needle_options_string NEEDLE_OPTIONS_STRING]
        [--keep_intermediate]
        [--dump]
        [--save_also_png]
        [-p N_PROCESSES]
        [--offset_around_cut_to_plot OFFSET_AROUND_CUT_TO_PLOT]
        [--min_frequency_alleles_around_cut_to_plot MIN_FREQUENCY_ALLELES_AROUND_CUT_TO_PLOT]
        [--max_rows_alleles_around_cut_to_plot MAX_ROWS_ALLELES_AROUND_CUT_TO_PLOT]
        [--debug]
    """

    # global variables for the multiprocessing
    global INCLUDE_IDXS
    global LEN_AMPLICON
    global EXON_POSITIONS
    global SPLICING_POSITIONS

    # check files
    check_file(args.fastq_r1)
    if args.fastq_r2:
        check_file(args.fastq_r2)

    # normalize name and remove not allowed characters
    if args.name:
        clean_name = slugify(args.name)
        if args.name != clean_name:
            warning(
                f"The specified name {args.name} contained "
                f"characters not allowed and was changed to: {clean_name}"
            )
            args.name = clean_name

    # amplicon sequence check
    # make everything uppercase!
    args.amplicon_seq = args.amplicon_seq.upper().strip().rstrip("\n")
    wrong_nt = find_wrong_nt(args.amplicon_seq)
    if wrong_nt:
        raise NTException(f"The amplicon sequence contains wrong characters:{wrong_nt}")

    LEN_AMPLICON = len(args.amplicon_seq)

    if args.guide_seq:
        cut_points = []
        sg_rna_intervals = []
        offset_plots = []
        sg_rna_sequences = []

        args.guide_seq = args.guide_seq.strip().upper()

        for current_guide_seq in args.guide_seq.split(","):

            if current_guide_seq in args.amplicon_seq:
                offset_plots.append(1)
            else:
                offset_plots.append(0)

            wrong_nt = find_wrong_nt(current_guide_seq)
            if wrong_nt:
                raise NTException(
                    f"The sgRNA sequence contains wrong characters:{wrong_nt}"
                )

            offset_fw = args.cleavage_offset + len(current_guide_seq) - 1
            offset_rc = (-args.cleavage_offset) - 1
            cut_points += [
                m.start() + offset_fw
                for m in re.finditer(current_guide_seq, args.amplicon_seq)
            ] + [
                m.start() + offset_rc
                for m in re.finditer(
                    reverse_complement(current_guide_seq), args.amplicon_seq
                )
            ]
            sg_rna_intervals += [
                (m.start(), m.start() + len(current_guide_seq) - 1)
                for m in re.finditer(current_guide_seq, args.amplicon_seq)
            ] + [
                (m.start(), m.start() + len(current_guide_seq) - 1)
                for m in re.finditer(
                    reverse_complement(current_guide_seq), args.amplicon_seq
                )
            ]
            sg_rna_sequences.append(current_guide_seq)

        offset_plots = np.array(offset_plots)

        if not cut_points:
            raise SgRNASequenceException(
                "The guide sequence/s provided is(are) "
                "not present in the amplicon sequence! "
                "\n\nPlease check your input!"
            )
        info(f"Cut Points from guide seq:{cut_points}")

    else:
        cut_points = []
        sg_rna_intervals = []
        offset_plots = np.array([])
        sg_rna_sequences = []

    if args.expected_hdr_amplicon_seq:
        args.expected_hdr_amplicon_seq = args.expected_hdr_amplicon_seq.strip().upper()

        if args.expected_hdr_amplicon_seq == args.amplicon_seq:
            raise AmpliconEqualDonorException(
                "The amplicon sequence expected after an HDR "
                "and the reference amplicon cannot be the same! "
                "\n\nPlease check your input!"
            )

        wrong_nt = find_wrong_nt(args.expected_hdr_amplicon_seq)
        if wrong_nt:
            raise NTException(
                "The amplicon sequence expected after an HDR "
                f"contains wrong characters:{wrong_nt}"
            )

        # if len(args.expected_hdr_amplicon_seq)!=len(args.amplicon_seq):
        aligned_ref, aligned_exp = pairwise2.align.globalxx(
            args.amplicon_seq, args.expected_hdr_amplicon_seq
        )[0][:2]
        identity_ref_rep = (
            sum([1.0 for a, b in zip(aligned_ref, aligned_exp) if a == b])
            / len(aligned_ref)
            * 100
        )
        if identity_ref_rep < args.min_identity_score:
            raise DonorSequenceException(
                "The amplicon sequence expected after an HDR "
                "should be provided as the reference amplicon "
                "sequence with the relevant part of the donor "
                "sequence replaced, and not just as the donor "
                "sequence. \n\nPlease check your input!"
            )

    if args.donor_seq:
        args.donor_seq = args.donor_seq.strip().upper()
        wrong_nt = find_wrong_nt(args.donor_seq)
        if wrong_nt:
            raise NTException(
                f"The donor sequence contains wrong characters:{wrong_nt}"
            )

        if args.donor_seq not in args.expected_hdr_amplicon_seq:
            raise CoreDonorSequenceNotContainedException(
                "The donor sequence provided is not present in the "
                "expected HDR amplicon sequence, or the expected HDR "
                "amplicon sequence parameter (-e) is not defined.  "
                "\n\nPlease check your input!"
            )

        positions_core_donor_seq = [
            (m.start(), m.start() + len(args.donor_seq))
            for m in re.finditer(
                f"(?={args.donor_seq})", args.expected_hdr_amplicon_seq
            )
        ]
        if len(positions_core_donor_seq) > 1:
            raise CoreDonorSequenceNotUniqueException(
                "The donor sequence provided is not unique in the "
                "expected HDR amplicon sequence.  \n\nPlease check your input!"
            )
        core_donor_seq_st_en = positions_core_donor_seq[0]

    # ##FRAMESHIFT SUPPORT###
    if args.coding_seq:

        perform_frameshift_analysis = True

        EXON_POSITIONS = set()
        exon_intervals = []
        SPLICING_POSITIONS = []

        for exon_seq in args.coding_seq.strip().upper().split(","):

            # check for wrong NT
            wrong_nt = find_wrong_nt(exon_seq)
            if wrong_nt:
                raise NTException(
                    f"The coding sequence contains wrong characters:{wrong_nt}"
                )

            st_exon = args.amplicon_seq.find(exon_seq)
            if st_exon < 0:
                raise ExonSequenceException(
                    "The coding subsequence/s provided:{exon_seq} is(are) not "
                    "contained in the amplicon sequence."
                )
            en_exon = st_exon + len(
                exon_seq
            )  # this do not include the upper bound as usual in python
            exon_intervals.append((st_exon, en_exon))
            EXON_POSITIONS = EXON_POSITIONS.union(set(range(st_exon, en_exon)))

            # consider 2 base pairs before and after each exon
            SPLICING_POSITIONS += [
                max(0, st_exon - 2),
                max(0, st_exon - 1),
                min(LEN_AMPLICON - 1, en_exon),
                min(LEN_AMPLICON - 1, en_exon + 1),
            ]

        EXON_POSITIONS = sorted(EXON_POSITIONS)

        # protect from the wrong splitting of exons by the
        # users to avoid false splicing sites
        SPLICING_POSITIONS = set(SPLICING_POSITIONS).difference(EXON_POSITIONS)

    else:
        perform_frameshift_analysis = False

    # we have insertions/deletions that change the concatenated
    # exon sequence lenght and the difference between the final sequence
    # and the original sequence lenght is not a multiple of 3
    modified_frameshift = 0

    # we have insertions/deletions that change the concatenated
    # exon sequence lenght and the difference between the final sequence
    # and the original sequence lenght is a multiple of 3. We
    # are in this case also when no indels are present but we have
    # substitutions
    modified_non_frameshift = 0

    # we don't touch the exons at all, the read can be still modified tough..
    non_modified_non_frameshift = 0

    splicing_sites_modified = 0

    ################

    get_name_from_fasta = (
        lambda x: os.path.basename(x).replace(".fastq", "").replace(".gz", "")
    )

    if not args.name:
        if args.fastq_r2 != "":
            database_id = (
                f"{get_name_from_fasta(args.fastq_r1)}_"
                f"{get_name_from_fasta(args.fastq_r2)}"
            )
        else:
            database_id = get_name_from_fasta(args.fastq_r1)

    else:
        database_id = args.name

    output_directory = f"CRISPResso_on_{database_id}"

    if args.output_folder:
        output_directory = os.path.join(
            os.path.abspath(args.output_folder), output_directory
        )

    _jp = lambda filename: os.path.join(
        output_directory, filename
    )  # handy function to put a file in the output directory
    log_filename = _jp("CRISPResso_RUNNING_LOG.txt")

    try:
        os.makedirs(output_directory)
        info(f"Creating Folder {output_directory}")
        info("Done!")
    except Exception:
        warning(f"Folder {output_directory} already exists.")

    finally:
        logging.getLogger().addHandler(logging.FileHandler(log_filename))

        with open(log_filename, "wt", encoding="utf-8") as outfile:
            outfile.write(
                f"[Command used]:\nCRISPResso {sys.argv}\n\n"
                f"Args: {repr(args)}\n\n[Execution log]:\n"
            )

    if args.split_paired_end:

        if args.fastq_r2 != "":
            raise Exception(
                "The option --split_paired_end is available "
                "only when a single fastq file is specified!"
            )

        info("Splitting paired end single fastq file in two files...")
        args.fastq_r1, args.fastq_r2 = split_paired_end_reads_single_file(
            args.fastq_r1,
            output_filename_r1=_jp(
                os.path.basename(args.fastq_r1.replace(".fastq", "")).replace(".gz", "")
                + "_splitted_r1.fastq.gz"
            ),
            output_filename_r2=_jp(
                os.path.basename(args.fastq_r1.replace(".fastq", "")).replace(".gz", "")
                + "_splitted_r2.fastq.gz"
            ),
        )
        splitted_files_to_remove = [args.fastq_r1, args.fastq_r2]

        info("Done!")

    if args.min_average_read_quality > 0 or args.min_single_bp_quality > 0:
        info(
            "Filtering reads with average "
            f"bp quality < {args.min_average_read_quality} "
            f"and single bp quality < {args.min_single_bp_quality} ..."
        )
        if args.fastq_r2 != "":
            args.fastq_r1, args.fastq_r2 = filter_pe_fastq_by_qual(
                args.fastq_r1,
                args.fastq_r2,
                output_filename_r1=_jp(
                    os.path.basename(args.fastq_r1.replace(".fastq", "")).replace(
                        ".gz", ""
                    )
                    + "_filtered.fastq.gz"
                ),
                output_filename_r2=_jp(
                    os.path.basename(args.fastq_r2.replace(".fastq", "")).replace(
                        ".gz", ""
                    )
                    + "_filtered.fastq.gz"
                ),
                min_bp_quality=args.min_average_read_quality,
                min_single_bp_quality=args.min_single_bp_quality,
            )
        else:
            args.fastq_r1 = filter_se_fastq_by_qual(
                args.fastq_r1,
                output_filename=_jp(
                    os.path.basename(args.fastq_r1)
                    .replace(".fastq", "")
                    .replace(".gz", "")
                    + "_filtered.fastq.gz"
                ),
                min_bp_quality=args.min_average_read_quality,
                min_single_bp_quality=args.min_single_bp_quality,
            )

    if args.fastq_r2 == "":  # single end reads

        # check if we need to trim
        if not args.trim_sequences:
            # create a symbolic link
            symlink_filename = _jp(os.path.basename(args.fastq_r1))
            force_symlink(os.path.abspath(args.fastq_r1), symlink_filename)
            output_forward_filename = symlink_filename
        else:
            output_forward_filename = _jp("reads.trimmed.fq.gz")
            options_str = args.trimmomatic_options_string.replace(
                "NexteraPE-PE.fa", "TruSeq3-SE.fa"
            )
            # Trimming with trimmomatic
            cmd = (
                f"java -jar {get_data('trimmomatic-0.33.jar')} "
                f"SE -phred33 {args.fastq_r1}  "
                f"{output_forward_filename} "
                f"{options_str} >>{log_filename} 2>&1"
            )
            # print cmd
            trimmomatic_status = sb.call(cmd, shell=True)

            if trimmomatic_status:
                raise TrimmomaticException(
                    "TRIMMOMATIC failed to run, please check the log file."
                )

        processed_output_filename = output_forward_filename

    else:  # paired end reads case

        if not args.trim_sequences:
            output_forward_paired_filename = args.fastq_r1
            output_reverse_paired_filename = args.fastq_r2
        else:
            info("Trimming sequences with Trimmomatic...")
            output_forward_paired_filename = _jp("output_forward_paired.fq.gz")
            output_forward_unpaired_filename = _jp("output_forward_unpaired.fq.gz")
            output_reverse_paired_filename = _jp("output_reverse_paired.fq.gz")
            output_reverse_unpaired_filename = _jp("output_reverse_unpaired.fq.gz")

            # Trimming with trimmomatic
            cmd = (
                f"java -jar {get_data('trimmomatic-0.33.jar')} "
                f"PE -phred33 {args.fastq_r1}  {args.fastq_r2} "
                f"{output_forward_paired_filename}  {output_forward_unpaired_filename} "
                f" {output_reverse_paired_filename}  {output_reverse_unpaired_filename} "
                f"{args.trimmomatic_options_string} >>{log_filename} 2>&1"
            )
            # print cmd
            trimmomatic_status = sb.call(cmd, shell=True)
            if trimmomatic_status:
                raise TrimmomaticException(
                    "TRIMMOMATIC failed to run, please check the log file."
                )

            info("Done!")

        info("Estimating average read length...")
        if get_n_reads_fastq(output_forward_paired_filename):
            avg_read_length = get_average_read_length_fastq(
                output_forward_paired_filename
            )
            std_fragment_length = int(LEN_AMPLICON * 0.1)
        else:
            raise NoReadsAfterQualityFiltering(
                "No reads survived the average or single bp quality filtering."
            )

        # Merging with Flash
        info("Merging paired sequences with Flash...")
        cmd = (
            f"flash {output_forward_paired_filename} "
            f"{output_reverse_paired_filename} "
            f"--allow-outies --max-overlap {args.max_paired_end_reads_overlap} "
            f"--min-overlap {args.min_paired_end_reads_overlap} "
            f"-f {LEN_AMPLICON} -r {avg_read_length} "
            f"-s {std_fragment_length}  -z -d {output_directory} >>{log_filename} 2>&1"
        )

        flash_status = sb.call(cmd, shell=True)
        if flash_status:
            raise FlashException("Flash failed to run, please check the log file.")

        info("Done!")

        flash_hist_filename = _jp("out.hist")
        flash_histogram_filename = _jp("out.histogram")
        flash_not_combined_1_filename = _jp("out.notCombined_1.fastq.gz")
        flash_not_combined_2_filename = _jp("out.notCombined_2.fastq.gz")

        processed_output_filename = _jp("out.extendedFrags.fastq.gz")

    # count reads
    n_reads_input = get_n_reads_fastq(args.fastq_r1)
    n_reads_after_preprocessing = get_n_reads_fastq(processed_output_filename)
    if n_reads_after_preprocessing == 0:
        raise NoReadsAfterQualityFiltering(
            "No reads in input or no reads survived the average or single bp quality filtering."
        )

    info("Preparing files for the alignment...")
    # parsing flash output and prepare the files for alignment

    database_fasta_filename = _jp(f"{database_id}_database.fa")
    needle_output_filename = _jp(f"needle_output_{database_id}.txt.gz")

    # write .fa file only for amplicon the rest we pipe trough awk on the fly!

    with open(database_fasta_filename, "wt", encoding="utf-8") as outfile:
        outfile.write(f">{database_id}\n{args.amplicon_seq}\n")

    if args.expected_hdr_amplicon_seq:
        database_repair_fasta_filename = _jp(f"{database_id}_database_repair.fa")
        needle_output_repair_filename = _jp(
            f"needle_output_repair_{database_id}.txt.gz"
        )

        with open(database_repair_fasta_filename, "wt", encoding="utf-8") as outfile:
            outfile.write(f">{database_id}\n{args.expected_hdr_amplicon_seq}\n")

    def parse_needle_output(needle_filename, name="seq", just_score=False):
        needle_data = []

        try:
            needle_infile = gzip.open(needle_filename, mode="r")

            line = needle_infile.readline().decode("UTF-8")

            while line:

                while line and ("# Aligned_sequences" not in line):
                    line = needle_infile.readline().decode("UTF-8")

                if line:
                    # print line
                    needle_infile.readline().decode("UTF-8")  # skip another line

                    line = needle_infile.readline().decode("UTF-8")
                    id_seq = line.split()[-1].replace("_", ":")

                    for _ in range(5):
                        needle_infile.readline().decode("UTF-8")

                    line = needle_infile.readline().decode("UTF-8")

                    identity_seq = eval(
                        line.strip()
                        .split(" ")[-1]
                        .replace("%", "")
                        .replace(")", "")
                        .replace("(", "")
                    )

                    if just_score:
                        needle_data.append([id_seq, identity_seq])
                    else:
                        for _ in range(7):
                            needle_infile.readline().decode("UTF-8")

                        line = needle_infile.readline().decode("UTF-8")
                        aln_ref_seq = line.split()[2]

                        aln_str = (
                            needle_infile.readline().decode("UTF-8")[21:].rstrip("\n")
                        )
                        line = needle_infile.readline().decode("UTF-8")
                        aln_query_seq = line.split()[2]
                        aln_query_len = line.split()[3]

                        needle_data.append(
                            [
                                id_seq,
                                identity_seq,
                                aln_query_len,
                                aln_ref_seq,
                                aln_str,
                                aln_query_seq,
                            ]
                        )

            if just_score:
                needle_infile.close()
                return pd.DataFrame(
                    needle_data, columns=["ID", "score_" + name]
                ).set_index("ID")
            else:
                needle_infile.close()
                return pd.DataFrame(
                    needle_data,
                    columns=[
                        "ID",
                        "score_" + name,
                        "length",
                        "ref_seq",
                        "align_str",
                        "align_seq",
                    ],
                ).set_index("ID")
        except Exception as exc:
            raise NeedleException("Failed to parse the output of needle!") from exc

    info("Aligning sequences...")
    # Alignment here

    cmd = (
        (
            f"cat {processed_output_filename} |"
            + (" gunzip |" if processed_output_filename.endswith(".gz") else " ")
        )
        + r""" awk 'NR % 4 == 1 {print ">" $0} NR % 4 ==2 {print $0}' """
        + " | sed 's/:/_/g' | needle "
        f"-asequence={database_fasta_filename} "
        "-bsequence=/dev/stdin -outfile=/dev/stdout "
        f"{args.needle_options_string} 2>> {log_filename}  "
        f"| gzip >{needle_output_filename}"
    )

    needle_output = sb.call(cmd, shell=True)
    if needle_output:
        raise NeedleException("Needle failed to run, please check the log file.")

    # If we have a donor sequence we just compare
    # the fq in the two cases and see which one aligns better
    if args.expected_hdr_amplicon_seq:

        cmd_repair = (
            (
                f"cat {processed_output_filename} |"
                + (" gunzip |" if processed_output_filename.endswith(".gz") else " ")
            )
            + r""" awk 'NR % 4 == 1 {print ">" $0} NR % 4 ==2 {print $0}' """
            + " | sed 's/:/_/g' | needle "
            f"-asequence={database_repair_fasta_filename} "
            "-bsequence=/dev/stdin -outfile=/dev/stdout "
            f"{args.needle_options_string} 2>> {log_filename}  "
            f"| gzip >{needle_output_repair_filename}"
        )
        needle_output = sb.call(cmd_repair, shell=True)

        if needle_output:
            raise NeedleException("Needle failed to run, please " "check the log file.")
        info("Done!")

    # merge the flow
    if args.expected_hdr_amplicon_seq:
        df_database = parse_needle_output(needle_output_filename, "ref")
        df_database_repair = parse_needle_output(
            needle_output_repair_filename, "repaired", just_score=True
        )

        df_database_and_repair = df_database.join(df_database_repair)

        del df_database
        del df_database_repair

        # find reads that failed to align and try on the reverse complement
        sr_not_aligned = df_database_and_repair.loc[
            (df_database_and_repair.score_ref < args.min_identity_score)
            & (df_database_and_repair.score_ref < args.min_identity_score)
        ].align_seq.apply(lambda x: x.replace("_", ""))

        # filter out not aligned reads
        df_database_and_repair = df_database_and_repair.loc[
            (df_database_and_repair.score_ref > args.min_identity_score)
            | (df_database_and_repair.score_repaired > args.min_identity_score)
        ]

        df_database_and_repair["score_diff"] = (
            df_database_and_repair.score_ref - df_database_and_repair.score_repaired
        )

        df_needle_alignment = df_database_and_repair

        del df_database_and_repair

    else:
        df_needle_alignment = parse_needle_output(needle_output_filename, "ref")

        sr_not_aligned = df_needle_alignment.loc[
            (df_needle_alignment.score_ref < args.min_identity_score)
        ].align_seq.apply(lambda x: x.replace("_", ""))
        # filter out not aligned reads
        df_needle_alignment = df_needle_alignment.loc[
            df_needle_alignment.score_ref > args.min_identity_score
        ]

    # check if the not aligned reads are in the reverse complement
    if sr_not_aligned.count():

        # write fastq_not_aligned
        fasta_not_aligned_filename = _jp("not_aligned_amplicon_forward.fa.gz")

        outfile = gzip.open(fasta_not_aligned_filename, "wt")

        for x0, x1 in sr_not_aligned.items():
            outfile.write(f">{x0}\n{x1}\n")

        # write reverse complement of ampl and expected amplicon
        database_rc_fasta_filename = _jp(f"{database_id}_database_rc.fa")
        needle_output_rc_filename = _jp(f"needle_output_rc_{database_id}.txt.gz")

        info("Align sequences to reverse complement of the amplicon...")

        with open(database_rc_fasta_filename, "wt", encoding="utf-8") as outfile:
            outfile.write(f">{database_id}\n{reverse_complement(args.amplicon_seq)}\n")

        if args.expected_hdr_amplicon_seq:
            database_repair_rc_fasta_filename = _jp(
                f"{database_id}_database_repair_rc.fa"
            )
            needle_output_repair_rc_filename = _jp(
                f"needle_output_repair_rc_{database_id}.txt.gz"
            )

            with open(
                database_repair_rc_fasta_filename, "wt", encoding="utf-8"
            ) as outfile:
                outfile.write(
                    f">{database_id}\n"
                    f"{reverse_complement(args.expected_hdr_amplicon_seq)}\n"
                )
        info("Done!")

        # Now we do the alignment
        cmd = (
            f"zcat < {fasta_not_aligned_filename} | sed 's/:/_/g' | "
            f"needle -asequence={database_rc_fasta_filename} "
            f"-bsequence=/dev/stdin -outfile=/dev/stdout "
            f"{args.needle_options_string} 2>> {log_filename}  "
            f"| gzip >{needle_output_rc_filename}"
        )

        needle_output = sb.call(cmd, shell=True)
        if needle_output:
            raise NeedleException("Needle failed to run, please check the log file.")

        if args.expected_hdr_amplicon_seq:
            cmd = (
                f"zcat < {fasta_not_aligned_filename} | sed 's/:/_/g' | "
                f"needle -asequence={database_repair_rc_fasta_filename} "
                "-bsequence=/dev/stdin -outfile=/dev/stdout "
                f"args.needle_options_string 2>> {log_filename}  "
                f"| gzip >{needle_output_repair_rc_filename}"
            )

            needle_output = sb.call(cmd, shell=True)
            if needle_output:
                raise NeedleException(
                    "Needle failed to run, please check the log file."
                )

        # merge the flow rev
        if args.expected_hdr_amplicon_seq:
            df_database_rc = parse_needle_output(needle_output_rc_filename, "ref")
            print("df_database_rc")
            print(df_database_rc)
            df_database_repair_rc = parse_needle_output(
                needle_output_repair_rc_filename, "repaired", just_score=True
            )
            print("df_database_repair_rc")
            print(df_database_repair_rc)

            df_database_and_repair_rc = df_database_rc.join(df_database_repair_rc)

            del df_database_rc
            del df_database_repair_rc

            # filter bad alignments also to rc

            df_database_and_repair_rc = df_database_and_repair_rc.loc[
                (df_database_and_repair_rc.score_ref > args.min_identity_score)
                | (df_database_and_repair_rc.score_repaired > args.min_identity_score)
            ]

            df_database_and_repair_rc["score_diff"] = (
                df_database_and_repair_rc.score_ref
                - df_database_and_repair_rc.score_repaired
            )

            df_needle_alignment_rc = df_database_and_repair_rc

            del df_database_and_repair_rc

        else:
            df_needle_alignment_rc = parse_needle_output(
                needle_output_rc_filename, "ref"
            )

            # filter out not aligned reads
            df_needle_alignment_rc = df_needle_alignment_rc.loc[
                df_needle_alignment_rc.score_ref > args.min_identity_score
            ]

        # reverse complement and invert the align string
        # so we have everything in the positive strand
        df_needle_alignment_rc["ref_seq"] = df_needle_alignment_rc["ref_seq"].apply(
            reverse_complement
        )
        df_needle_alignment_rc["align_seq"] = df_needle_alignment_rc["align_seq"].apply(
            reverse_complement
        )
        df_needle_alignment_rc["align_str"] = df_needle_alignment_rc["align_str"].apply(
            lambda x: x[::-1]
        )

        # fix for duplicates when rc alignment
        df_needle_alignment_rc.index = map(
            lambda x: "_".join([x, "RC"]), df_needle_alignment_rc.index
        )

        # append the RC reads to the aligned reads in the original orientation
        df_needle_alignment = pd.concat([df_needle_alignment, df_needle_alignment_rc])

        del df_needle_alignment_rc

    # check for duplicates
    try:
        assert (
            df_needle_alignment.shape[0] == df_needle_alignment.index.unique().shape[0]
        )
    except Exception as exc:
        raise DuplicateSequenceIdException(
            "The .fastq file/s contain/s duplicate sequence IDs"
        ) from exc

    # Initializations
    info("Quantifying indels/substitutions...")
    df_needle_alignment["UNMODIFIED"] = df_needle_alignment.score_ref == 100

    # the rest we have to look one by one to potentially exclude regions
    df_needle_alignment["MIXED"] = False
    df_needle_alignment["HDR"] = False
    df_needle_alignment["NHEJ"] = False

    df_needle_alignment["n_mutated"] = 0
    df_needle_alignment["n_inserted"] = 0
    df_needle_alignment["n_deleted"] = 0

    n_total = df_needle_alignment.shape[0]

    if n_total == 0:
        raise NoReadsAlignedException(
            "Zero sequences aligned, please check your amplicon sequence"
        )

    # remove the mutations in bp equal to 'N'
    if "N" in args.amplicon_seq:

        info(
            "Your amplicon sequence contains one or more N, "
            "excluding these bp for the indel quantification..."
        )

        def ignore_n_in_alignment(row):
            row["align_str"] = "".join(
                [
                    ("|" if (row["ref_seq"][idx] == "N") else c)
                    for idx, c in enumerate(row["align_str"])
                ]
            )
            if len(set(row["align_str"])) == 1:
                row["UNMODIFIED"] = True

            return row

        df_needle_alignment = df_needle_alignment.apply(ignore_n_in_alignment, axis=1)

    # ####QUANTIFICATION START
    def compute_ref_positions(ref_seq):
        pos_idxs = []
        idx = 0
        for c in ref_seq:
            if c in set(["A", "T", "C", "G", "N"]):
                pos_idxs.append(idx)
                idx += 1
            else:
                if idx == 0:
                    pos_idxs.append(-1)
                else:
                    pos_idxs.append(-idx)
        return np.array(pos_idxs)

    # compute positions relative to alignmnet
    df_needle_alignment["ref_positions"] = df_needle_alignment["ref_seq"].apply(
        compute_ref_positions
    )

    def plot1_indel_size(hlengths, hdensity, center_index, pdf, args):
        """Plot figure 1

        Args:
            hlengths
            hdensity
            center_index

        """

        # Create tag for edit type
        # No indel is the center of the bar plot
        # Rest are indels

        # Bar plot by counts
        indel = ["Indel"] * len(hlengths)
        indel[center_index] = "No indel"

        fig1_df = pd.DataFrame.from_dict(
            {
                "Indel size (bp)": hlengths,
                "Sequences (no.)": hdensity,
                "Edit Type": indel,
            }
        )

        fig1a = px.bar(
            data_frame=fig1_df,
            x="Indel size (bp)",
            y="Sequences (no.)",
            color="Edit Type",
            title="Indel size distribution",
        )

        fig1a.update_layout(legend_title="")

        filename_fig1a = "1a.Indel_size_distribution_n_sequences"
        fig1a.write_image(_jp(f"{filename_fig1a}.pdf"))
        fig1a.write_html(_jp(f"{filename_fig1a}.html"))
        if args.save_also_png:
            fig1a.write_image(_jp(f"{filename_fig1a}.png"))

        # Bar plot by percentages

        fig1_df = pd.DataFrame.from_dict(
            {
                "Indel size (bp)": hlengths,
                "Sequences (%)": 100.0 * hdensity / hdensity.sum(),
                "Edit Type": indel,
            }
        )

        fig1b = px.bar(
            data_frame=fig1_df,
            x="Indel size (bp)",
            y="Sequences (%)",
            color="Edit Type",
            title="Indel size distribution",
        )

        fig1b.update_layout(legend_title="")

        filename_fig1b = "1b.Indel_size_distribution_percentage"
        fig1b.write_image(_jp(f"{filename_fig1b}.pdf"))
        fig1b.write_html(_jp(f"{filename_fig1b}.html"))
        if args.save_also_png:
            fig1b.write_image(_jp(f"{filename_fig1b}.png"))

    def plot2(
        n_unmodified,
        n_mixed_hdr_nhej,
        n_modified,
        n_repaired,
        cut_points,
        sg_rna_intervals,
        length_amplicon,
        pdf,
        args,
    ):
        """Create plot 2"""

        # ###PIE CHARTS FOR HDR/NHEJ/MIXED/EVENTS###

        if args.expected_hdr_amplicon_seq:

            data = [
                [f"Unmodified\n({n_unmodified} reads)", n_unmodified],
                [f"Mixed HDR-NHEJ\n({n_mixed_hdr_nhej} reads)", n_mixed_hdr_nhej],
                [f"NHEJ\n({n_modified} reads)", n_modified],
                [f"HDR\n({n_repaired} reads)", n_repaired],
            ]

            df2 = pd.DataFrame(data, columns=["Type", "Reads"])

            fig2 = px.pie(df2, values="Reads", names="Type")
            fig2.update_traces(textinfo="percent+label")

            # plt.figure(figsize=(12 * 1.5, 14.5 * 1.5))
            # ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=5)
            # _, texts, autotexts = ax1.pie(
            #     [n_unmodified, n_mixed_hdr_nhej, n_modified, n_repaired],
            #     labels=[
            #         f"Unmodified\n({n_unmodified} reads)",
            #         f"Mixed HDR-NHEJ\n({n_mixed_hdr_nhej} reads)",
            #         f"NHEJ\n({n_modified} reads)",
            #         f"HDR\n({n_repaired} reads)",
            #     ],
            #     explode=(0, 0, 0, 0),
            #     colors=[(1, 0, 0, 0.2), (0, 1, 1, 0.2), (0, 0, 1, 0.2), (0, 1, 0, 0.2)],
            #     autopct="%1.1f%%",
            # )

            # if cut_points or args.donor_seq:
            #     ax2 = plt.subplot2grid((6, 3), (5, 0), colspan=3, rowspan=1)
            #     ax2.plot(
            #         [0, length_amplicon], [0, 0], "-k", lw=2, label="Amplicon sequence"
            #     )

            #     if args.donor_seq:
            #         ax2.plot(
            #             core_donor_seq_st_en,
            #             [0, 0],
            #             "-",
            #             lw=10,
            #             c=(0, 1, 0, 0.5),
            #             label="Donor Sequence",
            #         )

            #     if cut_points:
            #         ax2.plot(
            #             cut_points + offset_plots,
            #             np.zeros(len(cut_points)),
            #             "vr",
            #             ms=24,
            #             label="Predicted Cas9 cleavage site/s",
            #         )

            #     for idx, sg_rna_int in enumerate(sg_rna_intervals):
            #         if idx == 0:
            #             ax2.plot(
            #                 [sg_rna_int[0], sg_rna_int[1]],
            #                 [0, 0],
            #                 lw=10,
            #                 c=(0, 0, 0, 0.15),
            #                 label="sgRNA",
            #             )
            #         else:
            #             ax2.plot(
            #                 [sg_rna_int[0], sg_rna_int[1]],
            #                 [0, 0],
            #                 lw=10,
            #                 c=(0, 0, 0, 0.15),
            #                 label="_nolegend_",
            #             )

            #     lgd = plt.legend(
            #         bbox_to_anchor=(0, 0, 1.0, 0),
            #         ncol=1,
            #         mode="expand",
            #         borderaxespad=0.0,
            #         numpoints=1,
            #     )
            #     plt.xlim(0, length_amplicon)
            #     plt.axis("off")

            # proptease = fm.FontProperties()
            # proptease.set_size("xx-large")
            # plt.setp(autotexts, fontproperties=proptease)
            # plt.setp(texts, fontproperties=proptease)
            # plt.savefig(
            #     _jp("2.Unmodified_NHEJ_HDR_pie_chart.pdf"),
            #     pad_inches=1,
            #     bbox_inches="tight",
            # )
            # if args.save_also_png:
            #     plt.savefig(
            #         _jp("2.Unmodified_NHEJ_HDR_pie_chart.png"),
            #         pad_inches=1,
            #         bbox_inches="tight",
            #     )

            filename_fig2 = "2.Unmodified_NHEJ_HDR_pie_chart.png"
            fig2.write_image(_jp(f"{filename_fig2}.pdf"))
            fig2.write_html(_jp(f"{filename_fig2}.html"))
            if args.save_also_png:
                fig2.write_image(_jp(f"{filename_fig2}.png"))

        else:

            data = [
                [f"Unmodified\n({n_unmodified} reads)", n_unmodified],
                [f"NHEJ\n({n_modified} reads)", n_modified],
            ]

            df2 = pd.DataFrame(data, columns=["Type", "Reads"])

            fig2 = px.pie(df2, values="Reads", names="Type")
            fig2.update_traces(textinfo="percent+label")

            # plt.figure(figsize=(12 * 1.5, 14.5 * 1.5))
            # ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=5)
            # _, texts, autotexts = ax1.pie(
            #     [n_unmodified / n_total * 100, n_modified / n_total * 100],
            #     labels=[
            #         f"Unmodified\n({n_unmodified} reads)",
            #         f"NHEJ\n({n_modified} reads)",
            #     ],
            #     explode=(0, 0),
            #     colors=[(1, 0, 0, 0.2), (0, 0, 1, 0.2)],
            #     autopct="%1.1f%%",
            # )

            # if cut_points:
            #     ax2 = plt.subplot2grid((6, 3), (5, 0), colspan=3, rowspan=1)
            #     ax2.plot(
            #         [0, length_amplicon], [0, 0], "-k", lw=2, label="Amplicon sequence"
            #     )

            #     for idx, sg_rna_int in enumerate(sg_rna_intervals):
            #         if idx == 0:
            #             ax2.plot(
            #                 [sg_rna_int[0], sg_rna_int[1]],
            #                 [0, 0],
            #                 lw=10,
            #                 c=(0, 0, 0, 0.15),
            #                 label="sgRNA",
            #                 solid_capstyle="butt",
            #             )
            #         else:
            #             ax2.plot(
            #                 [sg_rna_int[0], sg_rna_int[1]],
            #                 [0, 0],
            #                 lw=10,
            #                 c=(0, 0, 0, 0.15),
            #                 label="_nolegend_",
            #                 solid_capstyle="butt",
            #             )

            #     ax2.plot(
            #         cut_points + offset_plots,
            #         np.zeros(len(cut_points)),
            #         "vr",
            #         ms=12,
            #         label="Predicted Cas9 cleavage site/s",
            #     )
            #     lgd = plt.legend(
            #         bbox_to_anchor=(0, 0, 1.0, 0),
            #         ncol=1,
            #         mode="expand",
            #         borderaxespad=0.0,
            #         numpoints=1,
            #         prop={"size": "large"},
            #     )
            #     plt.xlim(0, length_amplicon)
            #     plt.axis("off")

            # proptease = fm.FontProperties()
            # proptease.set_size("xx-large")
            # plt.setp(autotexts, fontproperties=proptease)
            # plt.setp(texts, fontproperties=proptease)
            # plt.savefig(
            #     _jp("2.Unmodified_NHEJ_pie_chart.pdf"),
            #     pad_inches=1,
            #     bbox_inches="tight",
            # )
            # if args.save_also_png:
            #     plt.savefig(
            #         _jp("2.Unmodified_NHEJ_pie_chart.png"),
            #         pad_inches=1,
            #         bbox_inches="tight",
            #     )

        # pdf.attach_note("Unmodified NEHJ pie chart")
        # pdf.savefig()  # saves the current figure into a pdf page
        # plt.close()

        filename_fig2 = "2.Unmodified_NHEJ_HDR_pie_chart.png"
        fig2.write_image(_jp(f"{filename_fig2}.pdf"))
        fig2.write_html(_jp(f"{filename_fig2}.html"))
        if args.save_also_png:
            fig2.write_image(_jp(f"{filename_fig2}.png"))

    def plot3(df_needle_alignment, n_total, pdf, args):
        """Plot 3"""

        # #############################################
        # (3) a graph of frequency of deletions and
        # insertions of various sizes (deletions
        # could be consider as negative numbers and insertions as positive);

        def calculate_range(df, column_name):
            df_not_zero = df.loc[df[column_name] > 0, column_name]
            try:
                r = max(15, int(np.round(np.percentile(df_not_zero, 99))))
            except Exception:
                r = 15
            return r

        range_mut = calculate_range(df_needle_alignment, "n_mutated")
        range_ins = calculate_range(df_needle_alignment, "n_inserted")
        range_del = calculate_range(df_needle_alignment, "n_deleted")

        y_values_mut, x_bins_mut = plt.histogram(
            df_needle_alignment["n_mutated"], bins=range(0, range_mut)
        )
        y_values_ins, x_bins_ins = plt.histogram(
            df_needle_alignment["n_inserted"], bins=range(0, range_ins)
        )
        y_values_del, x_bins_del = plt.histogram(
            df_needle_alignment["n_deleted"], bins=range(0, range_del)
        )

        fig = plt.figure(figsize=(26, 6.5))

        ax = fig.add_subplot(1, 3, 1)
        ax.bar(
            x_bins_ins[:-1], y_values_ins, align="center", linewidth=0, color=(0, 0, 1)
        )
        barlist = ax.bar(
            x_bins_ins[:-1], y_values_ins, align="center", linewidth=0, color=(0, 0, 1)
        )
        barlist[0].set_color("r")

        plt.title("Insertions")
        plt.xlabel("Size (bp)")
        plt.ylabel("Sequences % (no.)")
        lgd = plt.legend(
            ["Non-insertion", "Insertion"][::-1],
            bbox_to_anchor=(0.82, -0.22),
            ncol=1,
            fancybox=True,
            shadow=True,
        )
        lgd.legendHandles[0].set_height(6)
        lgd.legendHandles[1].set_height(6)
        plt.xlim(xmin=-1)
        y_label_values = np.round(np.linspace(0, min(n_total, max(ax.get_yticks())), 6))
        plt.yticks(
            y_label_values,
            [f"{n_reads / n_total * 100}%% ({n_reads})" for n_reads in y_label_values],
        )

        ax = fig.add_subplot(1, 3, 2)
        ax.bar(
            -x_bins_del[:-1], y_values_del, align="center", linewidth=0, color=(0, 0, 1)
        )
        barlist = ax.bar(
            -x_bins_del[:-1], y_values_del, align="center", linewidth=0, color=(0, 0, 1)
        )
        barlist[0].set_color("r")
        plt.title("Deletions")
        plt.xlabel("Size (bp)")
        plt.ylabel("Sequences % (no.)")
        lgd = plt.legend(
            ["Non-deletion", "Deletion"][::-1],
            bbox_to_anchor=(0.82, -0.22),
            ncol=1,
            fancybox=True,
            shadow=True,
        )
        lgd.legendHandles[0].set_height(6)
        lgd.legendHandles[1].set_height(6)
        plt.xlim(xmax=1)
        y_label_values = np.round(
            np.linspace(0, min(n_total, max(ax.get_yticks())), 6)
        )  # np.arange(0,y_max,y_max/6.0)
        plt.yticks(
            y_label_values,
            [
                "%.1f%% (%d)" % (n_reads / n_total * 100, n_reads)
                for n_reads in y_label_values
            ],
        )

        ax = fig.add_subplot(1, 3, 3)
        ax.bar(
            x_bins_mut[:-1], y_values_mut, align="center", linewidth=0, color=(0, 0, 1)
        )
        barlist = ax.bar(
            x_bins_mut[:-1], y_values_mut, align="center", linewidth=0, color=(0, 0, 1)
        )
        barlist[0].set_color("r")
        plt.title("Substitutions")
        plt.xlabel("Positions substituted (number)")
        plt.ylabel("Sequences % (no.)")
        lgd = plt.legend(
            ["Non-substitution", "Substitution"][::-1],
            bbox_to_anchor=(0.82, -0.22),
            ncol=1,
            fancybox=True,
            shadow=True,
        )
        lgd.legendHandles[0].set_height(6)
        lgd.legendHandles[1].set_height(6)
        plt.xlim(xmin=-1)
        y_label_values = np.round(np.linspace(0, min(n_total, max(ax.get_yticks())), 6))
        plt.yticks(
            y_label_values,
            [
                "%.1f%% (%d)" % (n_reads / n_total * 100, n_reads)
                for n_reads in y_label_values
            ],
        )

        plt.tight_layout()

        plt.savefig(
            _jp("3.Insertion_Deletion_Substitutions_size_hist.pdf"), bbox_inches="tight"
        )
        if args.save_also_png:
            plt.savefig(
                _jp("3.Insertion_Deletion_Substitutions_size_hist.png"),
                bbox_inches="tight",
            )

        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

        return (
            y_values_mut,
            x_bins_mut,
            y_values_ins,
            x_bins_ins,
            y_values_del,
            x_bins_del,
        )

    def plot4a(
        effect_vector_any, sg_rna_intervals, cut_points, length_amplicon, pdf, args
    ):
        """Create plot 4a"""

        # (4) another graph with the frequency that
        # each nucleotide within the amplicon was
        # modified in any way (perhaps would consider
        # insertion as modification of the flanking nucleotides);

        # Indels location Plots

        plt.figure(figsize=(10, 10))

        y_max = max(effect_vector_any) * 1.2

        plt.plot(
            effect_vector_any,
            "r",
            lw=3,
            label="Combined Insertions/Deletions/Substitutions",
        )
        # plt.hold(True)

        if cut_points:

            for idx, cut_point in enumerate(cut_points):
                if idx == 0:
                    plt.plot(
                        [cut_point + offset_plots[idx], cut_point + offset_plots[idx]],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="Predicted cleavage position",
                    )
                else:
                    plt.plot(
                        [cut_point + offset_plots[idx], cut_point + offset_plots[idx]],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="_nolegend_",
                    )

            for idx, sg_rna_int in enumerate(sg_rna_intervals):
                if idx == 0:
                    plt.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="sgRNA",
                        solid_capstyle="butt",
                    )
                else:
                    plt.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="_nolegend_",
                        solid_capstyle="butt",
                    )

        lgd = plt.legend(
            loc="center",
            bbox_to_anchor=(0.5, -0.23),
            ncol=1,
            fancybox=True,
            shadow=True,
        )
        if y_max > 0:
            y_label_values = np.arange(0, y_max, y_max / 6.0)
        plt.yticks(
            y_label_values,
            [
                f"{n_reads / float(n_total) * 100}%% ({n_reads})"
                for n_reads in y_label_values
            ],
        )
        plt.xticks(
            np.arange(
                0,
                length_amplicon,
                max(3, (length_amplicon // 6) - (length_amplicon // 6) % 5),
            )
        )

        plt.title("Mutation position distribution")
        plt.xlabel("Reference amplicon position (bp)")
        plt.ylabel("Sequences % (no.)")
        plt.ylim(0, max(1, y_max))
        plt.xlim(xmax=len(args.amplicon_seq) - 1)
        plt.savefig(
            _jp("4a.Combined_Insertion_Deletion_Substitution_Locations.pdf"),
            bbox_extra_artists=(lgd,),
            bbox_inches="tight",
        )
        if args.save_also_png:
            plt.savefig(
                _jp("4a.Combined_Insertion_Deletion_Substitution_Locations.png"),
                bbox_extra_artists=(lgd,),
                bbox_inches="tight",
                pad_inches=1,
            )

        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

    def plot4b(
        effect_vector_insertion,
        effect_vector_deletion,
        effect_vector_mutation,
        length_amplicon,
        n_total,
        n_modified,
        cut_points,
        sg_rna_intervals,
        pdf,
        args,
    ):
        """Plot 4b"""

        # NHEJ
        plt.figure(figsize=(10, 10))
        plt.plot(effect_vector_insertion, "r", lw=3, label="Insertions")

        plt.plot(effect_vector_deletion, "m", lw=3, label="Deletions")
        plt.plot(effect_vector_mutation, "g", lw=3, label="Substitutions")

        y_max = (
            max(
                max(effect_vector_insertion),
                max(effect_vector_deletion),
                max(effect_vector_mutation),
            )
            * 1.2
        )

        if cut_points:

            for idx, cut_point in enumerate(cut_points):
                if idx == 0:
                    plt.plot(
                        [cut_point + offset_plots[idx], cut_point + offset_plots[idx]],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="Predicted cleavage position",
                    )
                else:
                    plt.plot(
                        [cut_point + offset_plots[idx], cut_point + offset_plots[idx]],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="_nolegend_",
                    )

            for idx, sg_rna_int in enumerate(sg_rna_intervals):
                if idx == 0:
                    plt.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="sgRNA",
                        solid_capstyle="butt",
                    )
                else:
                    plt.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="_nolegend_",
                        solid_capstyle="butt",
                    )

        lgd = plt.legend(
            loc="center",
            bbox_to_anchor=(0.5, -0.28),
            ncol=1,
            fancybox=True,
            shadow=True,
        )
        if y_max > 0:
            y_label_values = np.arange(0, y_max, y_max / 6.0)
        plt.yticks(
            y_label_values,
            [
                "%.1f%% (%.1f%% , %d)"
                % (
                    n_reads / float(n_total) * 100,
                    n_reads / float(n_modified) * 100,
                    n_reads,
                )
                for n_reads in y_label_values
            ],
        )
        plt.xticks(
            np.arange(
                0,
                length_amplicon,
                max(3, (length_amplicon // 6) - (length_amplicon // 6) % 5),
            )
        )

        plt.xlabel("Reference amplicon position (bp)")
        plt.ylabel("Sequences: % Total ( % NHEJ, no. )")
        plt.ylim(0, max(1, y_max))
        plt.xlim(xmax=len(args.amplicon_seq) - 1)

        plt.title("Mutation position distribution of NHEJ")
        plt.savefig(
            _jp("4b.Insertion_Deletion_Substitution_Locations_NHEJ.pdf"),
            bbox_extra_artists=(lgd,),
            bbox_inches="tight",
        )
        if args.save_also_png:
            plt.savefig(
                _jp("4b.Insertion_Deletion_Substitution_Locations_NHEJ.png"),
                bbox_extra_artists=(lgd,),
                bbox_inches="tight",
                pad_inches=1,
            )

        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

    # END PLOTS

    # INITIALIZATIONS

    effect_vector_insertion = np.zeros(LEN_AMPLICON)
    effect_vector_deletion = np.zeros(LEN_AMPLICON)
    effect_vector_mutation = np.zeros(LEN_AMPLICON)
    effect_vector_any = np.zeros(LEN_AMPLICON)

    effect_vector_insertion_mixed = np.zeros(LEN_AMPLICON)
    effect_vector_deletion_mixed = np.zeros(LEN_AMPLICON)
    effect_vector_mutation_mixed = np.zeros(LEN_AMPLICON)

    effect_vector_insertion_hdr = np.zeros(LEN_AMPLICON)
    effect_vector_deletion_hdr = np.zeros(LEN_AMPLICON)
    effect_vector_mutation_hdr = np.zeros(LEN_AMPLICON)

    effect_vector_insertion_noncoding = np.zeros(LEN_AMPLICON)
    effect_vector_deletion_noncoding = np.zeros(LEN_AMPLICON)
    effect_vector_mutation_noncoding = np.zeros(LEN_AMPLICON)

    hist_inframe = defaultdict(lambda: 0)
    hist_frameshift = defaultdict(lambda: 0)

    avg_vector_del_all = np.zeros(LEN_AMPLICON)
    avg_vector_ins_all = np.zeros(LEN_AMPLICON)

    # look around the sgRNA(s) only?
    if cut_points and args.window_around_sgrna > 0:
        INCLUDE_IDXS = []
        half_window = max(1, args.window_around_sgrna // 2)
        for cut_p in cut_points:
            st = max(0, cut_p - half_window + 1)
            en = min(len(args.amplicon_seq) - 1, cut_p + half_window + 1)
            INCLUDE_IDXS.append(range(st, en))
    else:
        INCLUDE_IDXS = range(len(args.amplicon_seq))

    exclude_idxs = []

    if args.exclude_bp_from_left:
        exclude_idxs += range(args.exclude_bp_from_left)

    if args.exclude_bp_from_right:
        exclude_idxs += range(LEN_AMPLICON)[-args.exclude_bp_from_right :]

    # flatten the arrays to avoid errors with old numpy library
    INCLUDE_IDXS = np.ravel(INCLUDE_IDXS)
    exclude_idxs = np.ravel(exclude_idxs)

    INCLUDE_IDXS = set(np.setdiff1d(INCLUDE_IDXS, exclude_idxs))

    # handy generator to split in chunks the dataframe, np.split_array is slow!
    def get_chunk(df_needle_alignment, n_processes, args):
        for _, df in df_needle_alignment.groupby(
            np.arange(len(df_needle_alignment))
            // (len(df_needle_alignment) // (n_processes - 1))
        ):
            yield df, args

    # Use a Pool of processes, or just a single process
    if args.n_processes > 1:

        processes = min(df_needle_alignment.shape[0], args.n_processes)
        info(
            "[CRISPResso quantification is running in "
            f"parallel mode with {processes} processes]"
        )

        pool = mp.Pool(processes=processes)
        chunks_computed = []

        for result in pool.imap(
            process_df_chunk, get_chunk(df_needle_alignment, processes, args)
        ):
            (
                df_needle_alignment_chunk,
                effect_vector_insertion_chunk,
                effect_vector_deletion_chunk,
                effect_vector_mutation_chunk,
                effect_vector_any_chunk,
                effect_vector_insertion_mixed_chunk,
                effect_vector_deletion_mixed_chunk,
                effect_vector_mutation_mixed_chunk,
                effect_vector_insertion_hdr_chunk,
                effect_vector_deletion_hdr_chunk,
                effect_vector_mutation_hdr_chunk,
                effect_vector_insertion_noncoding_chunk,
                effect_vector_deletion_noncoding_chunk,
                effect_vector_mutation_noncoding_chunk,
                hist_inframe_chunk,
                hist_frameshift_chunk,
                avg_vector_del_all_chunk,
                avg_vector_ins_all_chunk,
                modified_frameshift_chunk,
                modified_non_frameshift_chunk,
                non_modified_non_frameshift_chunk,
                splicing_sites_modified_chunk,
            ) = result

            chunks_computed.append(df_needle_alignment_chunk)
            effect_vector_insertion += effect_vector_insertion_chunk
            effect_vector_deletion += effect_vector_deletion_chunk
            effect_vector_mutation += effect_vector_mutation_chunk
            effect_vector_any += effect_vector_any_chunk
            effect_vector_insertion_mixed += effect_vector_insertion_mixed_chunk
            effect_vector_deletion_mixed += effect_vector_deletion_mixed_chunk
            effect_vector_mutation_mixed += effect_vector_mutation_mixed_chunk
            effect_vector_insertion_hdr += effect_vector_insertion_hdr_chunk
            effect_vector_deletion_hdr += effect_vector_deletion_hdr_chunk
            effect_vector_mutation_hdr += effect_vector_mutation_hdr_chunk
            effect_vector_insertion_noncoding += effect_vector_insertion_noncoding_chunk
            effect_vector_deletion_noncoding += effect_vector_deletion_noncoding_chunk
            effect_vector_mutation_noncoding += effect_vector_mutation_noncoding_chunk
            add_hist(hist_inframe_chunk, hist_inframe)
            add_hist(hist_frameshift_chunk, hist_frameshift)
            avg_vector_del_all += avg_vector_del_all_chunk
            avg_vector_ins_all += avg_vector_ins_all_chunk
            modified_frameshift += modified_frameshift_chunk
            modified_non_frameshift += modified_non_frameshift_chunk
            non_modified_non_frameshift += non_modified_non_frameshift_chunk
            splicing_sites_modified += splicing_sites_modified_chunk

        pool.close()
        pool.join()
        df_needle_alignment = pd.concat(chunks_computed)
        del chunks_computed

    else:
        (
            df_needle_alignment,
            effect_vector_insertion,
            effect_vector_deletion,
            effect_vector_mutation,
            effect_vector_any,
            effect_vector_insertion_mixed,
            effect_vector_deletion_mixed,
            effect_vector_mutation_mixed,
            effect_vector_insertion_hdr,
            effect_vector_deletion_hdr,
            effect_vector_mutation_hdr,
            effect_vector_insertion_noncoding,
            effect_vector_deletion_noncoding,
            effect_vector_mutation_noncoding,
            hist_inframe,
            hist_frameshift,
            avg_vector_del_all,
            avg_vector_ins_all,
            modified_frameshift,
            modified_non_frameshift,
            non_modified_non_frameshift,
            splicing_sites_modified,
        ) = process_df_chunk([df_needle_alignment, args])

    n_modified = df_needle_alignment["NHEJ"].sum()
    n_unmodified = df_needle_alignment["UNMODIFIED"].sum()
    n_mixed_hdr_nhej = df_needle_alignment["MIXED"].sum()
    n_repaired = df_needle_alignment["HDR"].sum()

    # disable known division warning
    with np.errstate(divide="ignore", invalid="ignore"):

        effect_vector_combined = 100.0 * effect_vector_any / float(n_total)

        avg_vector_ins_all /= (
            effect_vector_insertion
            + effect_vector_insertion_hdr
            + effect_vector_insertion_mixed
        )
        avg_vector_del_all /= (
            effect_vector_deletion
            + effect_vector_deletion_hdr
            + effect_vector_deletion_mixed
        )

    avg_vector_ins_all[np.isnan(avg_vector_ins_all)] = 0
    avg_vector_del_all[np.isnan(avg_vector_del_all)] = 0
    avg_vector_ins_all[np.isinf(avg_vector_ins_all)] = 0
    avg_vector_del_all[np.isinf(avg_vector_del_all)] = 0

    if perform_frameshift_analysis:
        if not dict(hist_inframe):
            hist_inframe = {0: 0}

        if not dict(hist_frameshift):
            hist_frameshift = {0: 0}

    info("Done!")

    info("Calculating indel distribution based on the length of the reads...")

    df_needle_alignment["effective_len"] = df_needle_alignment.apply(
        lambda row: LEN_AMPLICON + row.n_inserted - row.n_deleted, axis=1
    )

    info("Done!")

    # write alleles table
    info("Calculating alleles frequencies...")

    def get_ref_positions(row, df_alignment):
        ref_positions = list(
            df_alignment.loc[[(row.Aligned_Sequence, row.Reference_Sequence)]]
            .iloc[
                0,
            ]
            .loc["ref_positions"]
        )

        return ref_positions

    df_alleles = df_needle_alignment.groupby(
        [
            "align_seq",
            "ref_seq",
            "NHEJ",
            "UNMODIFIED",
            "HDR",
            "n_deleted",
            "n_inserted",
            "n_mutated",
        ]
    ).size()
    df_alleles = df_alleles.reset_index()
    df_alleles.rename(
        columns={
            0: "#Reads",
            "align_seq": "Aligned_Sequence",
            "ref_seq": "Reference_Sequence",
        },
        inplace=True,
    )
    df_alleles["%Reads"] = df_alleles["#Reads"] / df_alleles["#Reads"].sum() * 100.0

    df_alleles.sort_values(by="#Reads", ascending=False, inplace=True)

    # add ref positions for the plot around the cut sites
    df_needle_alignment.set_index(["align_seq", "ref_seq"], inplace=True)
    df_needle_alignment.sort_index(inplace=True)
    df_alleles["ref_positions"] = df_alleles.apply(
        lambda x: get_ref_positions(x, df_needle_alignment), axis=1
    ).values

    info("Done!")

    info("Making Plots...")

    # plot effective length
    if args.guide_seq:
        min_cut = min(cut_points)
        max_cut = max(cut_points)
        xmin, xmax = -min_cut, LEN_AMPLICON - max_cut
    else:
        min_cut = LEN_AMPLICON // 2
        max_cut = LEN_AMPLICON // 2
        xmin, xmax = -min_cut, +max_cut

    hdensity, hlengths = np.histogram(
        df_needle_alignment.effective_len - LEN_AMPLICON, np.arange(xmin, xmax)
    )
    hlengths = hlengths[:-1]
    center_index = np.nonzero(hlengths == 0)[0][0]

    # Save concatenated report to PDF
    pdf = PdfPages(_jp(f"crispresso_report_for_{database_id}.pdf"))

    # Create plot 1
    plot1_indel_size(hlengths, hdensity, center_index, pdf, args)

    # Create plot 2
    plot2(
        n_unmodified,
        n_mixed_hdr_nhej,
        n_modified,
        n_repaired,
        cut_points,
        sg_rna_intervals,
        LEN_AMPLICON,
        pdf,
        args,
    )

    # Create plot 3
    (
        y_values_mut,
        x_bins_mut,
        y_values_ins,
        x_bins_ins,
        y_values_del,
        x_bins_del,
    ) = plot3(df_needle_alignment, n_total, pdf, args)

    plot4a(effect_vector_any, sg_rna_intervals, cut_points, LEN_AMPLICON, pdf, args)

    plot4b(
        effect_vector_insertion,
        effect_vector_deletion,
        effect_vector_mutation,
        LEN_AMPLICON,
        n_total,
        n_modified,
        cut_points,
        sg_rna_intervals,
        pdf,
        args,
    )

    if args.expected_hdr_amplicon_seq:

        # HDR
        plt.figure(figsize=(10, 10))
        plt.plot(effect_vector_insertion_hdr, "r", lw=3, label="Insertions")
        # plt.hold(True)
        plt.plot(effect_vector_deletion_hdr, "m", lw=3, label="Deletions")
        plt.plot(effect_vector_mutation_hdr, "g", lw=3, label="Substitutions")

        y_max = (
            max(
                max(effect_vector_insertion_hdr),
                max(effect_vector_deletion_hdr),
                max(effect_vector_mutation_hdr),
            )
            * 1.2
        )

        if cut_points:

            for idx, cut_point in enumerate(cut_points):
                if idx == 0:
                    plt.plot(
                        [
                            cut_point + offset_plots[idx],
                            cut_point + offset_plots[idx],
                        ],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="Predicted cleavage position",
                    )
                else:
                    plt.plot(
                        [
                            cut_point + offset_plots[idx],
                            cut_point + offset_plots[idx],
                        ],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="_nolegend_",
                    )

            for idx, sg_rna_int in enumerate(sg_rna_intervals):
                if idx == 0:
                    plt.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="sgRNA",
                        solid_capstyle="butt",
                    )
                else:
                    plt.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="_nolegend_",
                        solid_capstyle="butt",
                    )

        lgd = plt.legend(
            loc="center",
            bbox_to_anchor=(0.5, -0.28),
            ncol=1,
            fancybox=True,
            shadow=True,
        )
        if y_max > 0:
            y_label_values = np.arange(0, y_max, y_max // 6).astype(int)
        plt.yticks(
            y_label_values,
            [
                "%.1f%% (%.1f%% , %d)"
                % (
                    n_reads / float(n_total) * 100.0,
                    n_reads / float(n_repaired) * 100.0,
                    n_reads,
                )
                for n_reads in y_label_values
            ],
        )
        plt.xticks(
            np.arange(
                0,
                LEN_AMPLICON,
                max(3, (LEN_AMPLICON // 6) - (LEN_AMPLICON // 6) % 5),
            ).astype(int)
        )

        plt.xlabel("Reference amplicon position (bp)")
        plt.ylabel("Sequences: % Total ( % HDR, no. )")
        plt.ylim(0, max(1, y_max))
        plt.xlim(xmax=len(args.amplicon_seq) - 1)
        plt.title("Mutation position distribution of HDR")
        plt.savefig(
            _jp("4c.Insertion_Deletion_Substitution_Locations_HDR.pdf"),
            bbox_extra_artists=(lgd,),
            bbox_inches="tight",
        )
        if args.save_also_png:
            plt.savefig(
                _jp("4c.Insertion_Deletion_Substitution_Locations_HDR.png"),
                bbox_extra_artists=(lgd,),
                bbox_inches="tight",
                pad_inches=1,
            )

        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

        # MIXED
        plt.figure(figsize=(10, 10))
        plt.plot(effect_vector_insertion_mixed, "r", lw=3, label="Insertions")
        # plt.hold(True)
        plt.plot(effect_vector_deletion_mixed, "m", lw=3, label="Deletions")
        plt.plot(effect_vector_mutation_mixed, "g", lw=3, label="Substitutions")

        y_max = (
            max(
                max(effect_vector_insertion_mixed),
                max(effect_vector_deletion_mixed),
                max(effect_vector_mutation_mixed),
            )
            * 1.2
        )

        if cut_points:

            for idx, cut_point in enumerate(cut_points):
                if idx == 0:
                    plt.plot(
                        [
                            cut_point + offset_plots[idx],
                            cut_point + offset_plots[idx],
                        ],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="Predicted cleavage position",
                    )
                else:
                    plt.plot(
                        [
                            cut_point + offset_plots[idx],
                            cut_point + offset_plots[idx],
                        ],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="_nolegend_",
                    )

            for idx, sg_rna_int in enumerate(sg_rna_intervals):
                if idx == 0:
                    plt.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="sgRNA",
                        solid_capstyle="butt",
                    )
                else:
                    plt.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="_nolegend_",
                        solid_capstyle="butt",
                    )

        lgd = plt.legend(
            loc="center",
            bbox_to_anchor=(0.5, -0.28),
            ncol=1,
            fancybox=True,
            shadow=True,
        )
        if y_max > 0:
            y_label_values = np.arange(0, y_max, y_max // 6).astype(int)
        plt.yticks(
            y_label_values,
            [
                "%.1f%% (%.1f%% , %d)"
                % (
                    n_reads / float(n_total) * 100.0,
                    n_reads / float(n_mixed_hdr_nhej) * 100.0,
                    n_reads,
                )
                for n_reads in y_label_values
            ],
        )
        plt.xticks(
            np.arange(
                0,
                LEN_AMPLICON,
                max(3, (LEN_AMPLICON // 6) - (LEN_AMPLICON // 6) % 5),
            ).astype(int)
        )

        plt.xlabel("Reference amplicon position (bp)")
        plt.ylabel("Sequences: % Total ( % mixed HDR-NHEJ, no. )")
        plt.ylim(0, max(1, y_max))
        plt.xlim(xmax=len(args.amplicon_seq) - 1)
        plt.title("Mutation position distribution of mixed HDR-NHEJ")
        plt.savefig(
            _jp("4d.Insertion_Deletion_Substitution_Locations_Mixed_HDR_NHEJ.pdf"),
            bbox_extra_artists=(lgd,),
            bbox_inches="tight",
        )
        if args.save_also_png:
            plt.savefig(
                _jp("4d.Insertion_Deletion_Substitution_Locations_Mixed_HDR_NHEJ.png"),
                bbox_extra_artists=(lgd,),
                bbox_inches="tight",
                pad_inches=1,
            )
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

    # Position dependent indels plot
    fig = plt.figure(figsize=(24, 10))
    ax1 = fig.add_subplot(1, 2, 1)

    markerline, stemlines, baseline = ax1.stem(avg_vector_ins_all, markerfmt="s")

    plt.setp(markerline, "markerfacecolor", "r", "markersize", 8)
    plt.setp(baseline, "linewidth", 0)
    plt.setp(stemlines, "color", "r", "linewidth", 3)

    # plt.hold(True)
    y_max = max(avg_vector_ins_all) * 1.2

    if cut_points:

        for idx, cut_point in enumerate(cut_points):

            if idx == 0:
                ax1.plot(
                    [cut_point + offset_plots[idx], cut_point + offset_plots[idx]],
                    [0, y_max],
                    "--k",
                    lw=2,
                    label="Predicted cleavage position",
                )
            else:
                ax1.plot(
                    [cut_point + offset_plots[idx], cut_point + offset_plots[idx]],
                    [0, y_max],
                    "--k",
                    lw=2,
                    label="_nolegend_",
                )

    plt.xticks(
        np.arange(
            0, LEN_AMPLICON, max(3, (LEN_AMPLICON // 6) - (LEN_AMPLICON // 6) % 5)
        ).astype(int)
    )
    plt.xlabel("Reference amplicon position (bp)")
    plt.ylabel("Average insertion length")
    plt.ylim(0, max(1, y_max))
    plt.xlim(xmax=LEN_AMPLICON - 1)
    ax1.set_title("Position dependent insertion size")
    plt.tight_layout()

    ax2 = fig.add_subplot(1, 2, 2)
    markerline, stemlines, baseline = ax2.stem(avg_vector_del_all, markerfmt="s")
    plt.setp(markerline, "markerfacecolor", "m", "markersize", 8)
    plt.setp(baseline, "linewidth", 0)
    plt.setp(stemlines, "color", "m", "linewidth", 3)

    y_max = max(avg_vector_del_all) * 1.2
    if cut_points:

        for idx, cut_point in enumerate(cut_points):
            if idx == 0:
                ax2.plot(
                    [cut_point + offset_plots[idx], cut_point + offset_plots[idx]],
                    [0, y_max],
                    "--k",
                    lw=2,
                    label="Predicted cleavage position",
                )
            else:
                ax2.plot(
                    [cut_point + offset_plots[idx], cut_point + offset_plots[idx]],
                    [0, y_max],
                    "--k",
                    lw=2,
                    label="_nolegend_",
                )

    plt.xticks(
        np.arange(
            0, LEN_AMPLICON, max(3, (LEN_AMPLICON // 6) - (LEN_AMPLICON // 6) % 5)
        ).astype(int)
    )
    plt.xlabel("Reference amplicon position (bp)")
    plt.ylabel("Average deletion length")

    plt.ylim(ymin=0, ymax=max(1, y_max))
    plt.xlim(xmax=LEN_AMPLICON - 1)
    ax2.set_title("Position dependent deletion size")

    plt.tight_layout()

    lgd = plt.legend(
        loc="center",
        bbox_to_anchor=(0.5, -0.28),
        ncol=1,
        fancybox=True,
        shadow=True,
    )

    plt.savefig(
        _jp("4e.Position_dependent_average_indel_size.pdf"),
        bbox_extra_artists=(lgd,),
        bbox_inches="tight",
    )
    if args.save_also_png:
        plt.savefig(
            _jp("4e.Position_dependent_average_indel_size.png"),
            bbox_extra_artists=(lgd,),
            bbox_inches="tight",
        )

    pdf.savefig()  # saves the current figure into a pdf page

    if perform_frameshift_analysis:

        # make frameshift plots
        fig = plt.figure(figsize=(12 * 1.5, 14.5 * 1.5))
        ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=5)
        _, texts, autotexts = ax1.pie(
            [
                modified_frameshift,
                modified_non_frameshift,
                non_modified_non_frameshift,
            ],
            labels=[
                f"Frameshift mutation\n({modified_frameshift} reads)",
                f"In-frame mutation\n({modified_non_frameshift} reads)",
                f"Noncoding mutation\n({non_modified_non_frameshift} reads)",
            ],
            explode=(0.0, 0.0, 0.0),
            colors=[
                (0.89019608, 0.29019608, 0.2, 0.8),
                (0.99215686, 0.73333333, 0.51764706, 0.8),
                (0.99607843, 0.90980392, 0.78431373, 0.8),
            ],
            autopct="%1.1f%%",
        )

        ax2 = plt.subplot2grid((6, 3), (5, 0), colspan=3, rowspan=1)
        ax2.plot([0, LEN_AMPLICON], [0, 0], "-k", lw=2, label="Amplicon sequence")

        for idx, exon_interval in enumerate(exon_intervals):
            if idx == 0:
                ax2.plot(
                    exon_interval,
                    [0, 0],
                    "-",
                    lw=10,
                    c=(0, 0, 1, 0.5),
                    label="Coding sequence/s",
                    solid_capstyle="butt",
                )
            else:
                ax2.plot(
                    exon_interval,
                    [0, 0],
                    "-",
                    lw=10,
                    c=(0, 0, 1, 0.5),
                    label="_nolegend_",
                    solid_capstyle="butt",
                )

        if cut_points:
            ax2.plot(
                cut_points + offset_plots,
                np.zeros(len(cut_points)),
                "vr",
                ms=25,
                label="Predicted Cas9 cleavage site/s",
            )

        lgd = plt.legend(
            bbox_to_anchor=(0, 0, 1.0, 0),
            ncol=1,
            mode="expand",
            borderaxespad=0.0,
            numpoints=1,
        )
        plt.xlim(0, LEN_AMPLICON)
        plt.axis("off")

        proptease = fm.FontProperties()
        proptease.set_size("xx-large")
        plt.setp(autotexts, fontproperties=proptease)
        plt.setp(texts, fontproperties=proptease)
        plt.savefig(
            _jp("5.Frameshift_In-frame_mutations_pie_chart.pdf"),
            pad_inches=1,
            bbox_inches="tight",
        )
        if args.save_also_png:
            plt.savefig(
                _jp("5.Frameshift_In-frame_mutations_pie_chart.png"),
                pad_inches=1,
                bbox_inches="tight",
            )

        pdf.savefig()  # saves the current figure into a pdf page

        # profiles--------------------------------------------
        fig = plt.figure(figsize=(22, 10))
        ax1 = fig.add_subplot(2, 1, 1)
        x, y = map(np.array, zip(*[a for a in hist_frameshift.items()]))
        y = y / float(sum(hist_frameshift.values())) * 100
        ax1.bar(x - 0.5, y)
        ax1.set_xlim(-30.5, 30.5)
        ax1.set_frame_on(False)
        ax1.set_xticks([idx for idx in range(-30, 31) if idx % 3])
        ax1.tick_params(
            which="both",  # both major and minor ticks are affected
            bottom="off",  # ticks along the bottom edge are off
            top="off",  # ticks along the top edge are off
            labelbottom="on",
        )  # labels along the bottom edge are off)
        ax1.yaxis.tick_left()
        xmin, xmax = ax1.get_xaxis().get_view_interval()
        ax1.set_xticklabels(
            [str(idx) for idx in [idx for idx in range(-30, 31) if idx % 3]],
            rotation="vertical",
        )
        plt.title("Frameshift profile")
        ax1.tick_params(axis="both", which="major", labelsize=32)
        ax1.tick_params(axis="both", which="minor", labelsize=32)
        plt.tight_layout()
        plt.ylabel("%")

        ax2 = fig.add_subplot(2, 1, 2)
        x, y = map(np.array, zip(*list(hist_inframe.items())))
        y = y / float(sum(hist_inframe.values())) * 100
        ax2.bar(x - 0.5, y, color=(0, 1, 1, 0.2))
        ax2.set_xlim(-30.5, 30.5)
        ax2.set_frame_on(False)
        ax2.set_xticks([idx for idx in range(-30, 31) if (idx % 3 == 0)])
        ax2.tick_params(
            which="both",  # both major and minor ticks are affected
            bottom="off",  # ticks along the bottom edge are off
            top="off",  # ticks along the top edge are off
            labelbottom="on",
        )  # labels along the bottom edge are off)
        ax2.yaxis.tick_left()
        xmin, xmax = ax2.xaxis.get_view_interval()

        ax2.set_xticklabels(
            [str(idx) for idx in [idx for idx in range(-30, 31) if (idx % 3 == 0)]],
            rotation="vertical",
        )
        plt.title("In-frame profile")
        plt.tight_layout()
        plt.ylabel("%")
        ax2.tick_params(axis="both", which="major", labelsize=32)
        ax2.tick_params(axis="both", which="minor", labelsize=32)
        plt.tight_layout()

        plt.savefig(
            _jp("6.Frameshift_In-frame_mutation_profiles.pdf"),
            pad_inches=1,
            bbox_inches="tight",
        )
        if args.save_also_png:
            plt.savefig(
                _jp("6.Frameshift_In-frame_mutation_profiles.png"),
                pad_inches=1,
                bbox_inches="tight",
            )

        pdf.savefig()  # saves the current figure into a pdf page

        # ---------------------------------------------------
        fig = plt.figure(figsize=(12 * 1.5, 12 * 1.5))
        ax = fig.add_subplot(1, 1, 1)
        _, texts, autotexts = ax.pie(
            [
                splicing_sites_modified,
                (df_needle_alignment.shape[0] - splicing_sites_modified),
            ],
            labels=[
                f"Potential splice sites modified\n({splicing_sites_modified}) reads)",
                f"Unmodified\n({df_needle_alignment.shape[0] - splicing_sites_modified} reads)",
            ],
            explode=(0.0, 0),
            colors=[
                (0.89019608, 0.29019608, 0.2, 0.8),
                (0.99607843, 0.90980392, 0.78431373, 0.8),
            ],
            autopct="%1.1f%%",
        )
        proptease = fm.FontProperties()
        proptease.set_size("xx-large")
        plt.setp(autotexts, fontproperties=proptease)
        plt.setp(texts, fontproperties=proptease)
        plt.savefig(
            _jp("8.Potential_Splice_Sites_pie_chart.pdf"),
            pad_inches=1,
            bbox_inches="tight",
        )
        if args.save_also_png:
            plt.savefig(
                _jp("8.Potential_Splice_Sites_pie_chart.png"),
                pad_inches=1,
                bbox_inches="tight",
            )

        pdf.savefig()  # saves the current figure into a pdf page
        plt.close("all")

        # non coding
        plt.figure(figsize=(10, 10))
        plt.plot(effect_vector_insertion_noncoding, "r", lw=3, label="Insertions")
        plt.plot(effect_vector_deletion_noncoding, "m", lw=3, label="Deletions")
        plt.plot(effect_vector_mutation_noncoding, "g", lw=3, label="Substitutions")

        y_max = (
            max(
                max(effect_vector_insertion_noncoding),
                max(effect_vector_deletion_noncoding),
                max(effect_vector_mutation_noncoding),
            )
            * 1.2
        )

        if cut_points:

            for idx, cut_point in enumerate(cut_points):
                if idx == 0:
                    plt.plot(
                        [
                            cut_point + offset_plots[idx],
                            cut_point + offset_plots[idx],
                        ],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="Predicted cleavage position",
                    )
                else:
                    plt.plot(
                        [
                            cut_point + offset_plots[idx],
                            cut_point + offset_plots[idx],
                        ],
                        [0, y_max],
                        "--k",
                        lw=2,
                        label="_nolegend_",
                    )

                for idx, sg_rna_int in enumerate(sg_rna_intervals):
                    if idx == 0:
                        plt.plot(
                            [sg_rna_int[0], sg_rna_int[1]],
                            [0, 0],
                            lw=10,
                            c=(0, 0, 0, 0.15),
                            label="sgRNA",
                            solid_capstyle="butt",
                        )
                    else:
                        plt.plot(
                            [sg_rna_int[0], sg_rna_int[1]],
                            [0, 0],
                            lw=10,
                            c=(0, 0, 0, 0.15),
                            label="_nolegend_",
                            solid_capstyle="butt",
                        )

        lgd = plt.legend(
            loc="center",
            bbox_to_anchor=(0.5, -0.28),
            ncol=1,
            fancybox=True,
            shadow=True,
        )
        plt.xticks(
            np.arange(
                0,
                LEN_AMPLICON,
                max(3, (LEN_AMPLICON // 6) - (LEN_AMPLICON // 6) % 5),
            ).astype(int)
        )

        plt.xlabel("Reference amplicon position (bp)")
        plt.ylabel("Sequences (no.)")
        plt.ylim(0, max(1, y_max))
        plt.xlim(xmax=len(args.amplicon_seq) - 1)
        plt.title("Noncoding mutation position distribution")
        plt.savefig(
            _jp("7.Insertion_Deletion_Substitution_Locations_Noncoding.pdf"),
            bbox_extra_artists=(lgd,),
            bbox_inches="tight",
        )
        if args.save_also_png:
            plt.savefig(
                _jp("7.Insertion_Deletion_Substitution_Locations_Noncoding.png"),
                bbox_extra_artists=(lgd,),
                bbox_inches="tight",
            )

    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    ##new plots alleles around cut_sites
    for sg_rna, cut_point in zip(sg_rna_sequences, cut_points):
        # print sg_rna,cut_point
        df_allele_around_cut = get_dataframe_around_cut(
            df_alleles, cut_point, args.offset_around_cut_to_plot
        )

        # write alleles table to file
        df_allele_around_cut.to_csv(
            _jp(f"Alleles_frequency_table_around_cut_site_for_{sg_rna}.txt"),
            sep="\t",
            header=True,
        )
        plot_alleles_table(
            args, cut_point, df_allele_around_cut, sg_rna, output_directory
        )

        pdf.savefig()
        plt.close("all")

    # We can also set the file's metadata via the PdfPages object:
    d = pdf.infodict()
    d["Title"] = f"CRISPResso Report for {database_id}"
    d["Subject"] = f"Collated results from CRISPResso {database_id}"
    d["Keywords"] = "CRISPResso Amplicon NGS CRISPR 'Gene Therapy'"
    d["CreationDate"] = datetime.datetime.today()
    pdf.close()
    info("Done!")

    if not args.keep_intermediate:
        info("Removing Intermediate files...")

        if args.fastq_r2 != "":
            files_to_remove = [
                processed_output_filename,
                flash_hist_filename,
                flash_histogram_filename,
                flash_not_combined_1_filename,
                flash_not_combined_2_filename,
                database_fasta_filename,
            ]
        else:
            files_to_remove = [processed_output_filename, database_fasta_filename]

        if args.trim_sequences and args.fastq_r2 != "":
            files_to_remove += [
                output_forward_paired_filename,
                output_reverse_paired_filename,
                output_forward_unpaired_filename,
                output_reverse_unpaired_filename,
            ]

        if not args.dump:
            files_to_remove += [needle_output_filename]
            if args.expected_hdr_amplicon_seq:
                files_to_remove += [needle_output_repair_filename]

        if args.expected_hdr_amplicon_seq:
            files_to_remove += [
                database_repair_fasta_filename,
            ]

        if args.split_paired_end:
            files_to_remove += splitted_files_to_remove

        if args.min_average_read_quality > 0 or args.min_single_bp_quality > 0:

            if args.fastq_r2 != "":
                files_to_remove += [args.fastq_r1, args.fastq_r2]
            else:
                files_to_remove += [args.fastq_r1]

        if sr_not_aligned.count():
            files_to_remove += [
                fasta_not_aligned_filename,
                database_rc_fasta_filename,
                needle_output_rc_filename,
            ]

            if args.expected_hdr_amplicon_seq:
                files_to_remove += [
                    database_repair_rc_fasta_filename,
                    needle_output_repair_rc_filename,
                ]

        for file_to_remove in files_to_remove:
            try:
                if os.path.islink(file_to_remove):
                    os.unlink(file_to_remove)
                else:
                    os.remove(file_to_remove)
            except Exception as exc:
                warning(f"{exc}: Skipping {file_to_remove}")

    # write effect vectors as plain text files
    info("Saving processed data...")

    def save_vector_to_file(vector, name):
        np.savetxt(
            _jp(f"{name}.txt"),
            np.vstack([(np.arange(len(vector)) + 1), vector]).T,
            fmt=["%d", "%.18e"],
            delimiter="\t",
            newline="\n",
            header="amplicon position\teffect",
            footer="",
            comments="# ",
        )

    nhej_inserted = np.sum(
        df_needle_alignment[df_needle_alignment["NHEJ"] == True]["n_inserted"] > 0
    )
    if np.isnan(nhej_inserted):
        nhej_inserted = 0

    nhej_deleted = np.sum(
        df_needle_alignment[df_needle_alignment["NHEJ"] == True]["n_deleted"] > 0
    )
    if np.isnan(nhej_deleted):
        nhej_deleted = 0

    nhej_mutated = np.sum(
        df_needle_alignment[df_needle_alignment["NHEJ"] == True]["n_mutated"] > 0
    )
    if np.isnan(nhej_mutated):
        nhej_mutated = 0

    hdr_inserted = np.sum(
        df_needle_alignment[df_needle_alignment["HDR"] == True]["n_inserted"] > 0
    )
    if np.isnan(hdr_inserted):
        hdr_inserted = 0

    hdr_deleted = np.sum(
        df_needle_alignment[df_needle_alignment["HDR"] == True]["n_deleted"] > 0
    )
    if np.isnan(hdr_deleted):
        hdr_deleted = 0

    hdr_mutated = np.sum(
        df_needle_alignment[df_needle_alignment["HDR"] == True]["n_mutated"] > 0
    )
    if np.isnan(hdr_mutated):
        hdr_mutated = 0

    mixed_inserted = np.sum(
        df_needle_alignment[df_needle_alignment["MIXED"] == True]["n_inserted"] > 0
    )
    if np.isnan(mixed_inserted):
        mixed_inserted = 0

    mixed_deleted = np.sum(
        df_needle_alignment[df_needle_alignment["MIXED"] == True]["n_deleted"] > 0
    )
    if np.isnan(mixed_deleted):
        mixed_deleted = 0

    mixed_mutated = np.sum(
        df_needle_alignment[df_needle_alignment["MIXED"] == True]["n_mutated"] > 0
    )
    if np.isnan(mixed_mutated):
        mixed_mutated = 0

    with open(
        _jp("Quantification_of_editing_frequency.txt"), "wt", encoding="utf-8"
    ) as outfile:
        outfile.write(
            (
                "Quantification of editing frequency:\n\t- "
                f"Unmodified:{n_unmodified} reads\n"
            )
            + (
                f"\t- NHEJ:{n_modified} reads "
                f"({nhej_inserted} reads with insertions, "
                f"{nhej_deleted} reads with deletions, "
                f"{nhej_mutated} reads with substitutions)\n"
            )
            + (
                f"\t- HDR:{n_repaired} reads "
                f"({hdr_inserted} reads with insertions, "
                f"{hdr_deleted} reads with deletions, "
                f"{hdr_mutated} reads with substitutions)\n"
            )
            + (
                f"\t- Mixed HDR-NHEJ:{n_mixed_hdr_nhej} reads "
                f"({mixed_inserted} reads with insertions, "
                f"{mixed_deleted} reads with deletions, "
                f"{mixed_mutated} reads with substitutions)\n\n"
            )
            + (f"Total Aligned:{n_total} reads ")
        )

    # write alleles table
    df_alleles.loc[:, :"%Reads"].to_csv(
        _jp("Alleles_frequency_table.txt"), sep="\t", header=True, index=None
    )

    # write statistics
    with open(_jp("Mapping_statistics.txt"), "wt", encoding="utf-8") as outfile:
        outfile.write(
            f"READS IN INPUTS:{n_reads_input}\n"
            f"READS AFTER PREPROCESSING:{n_reads_after_preprocessing}"
            f"\nREADS ALIGNED:{n_total}"
        )

    if perform_frameshift_analysis:
        with open(_jp("Frameshift_analysis.txt"), "wt", encoding="utf-8") as outfile:
            outfile.write(
                "Frameshift analysis:\n\t"
                f"Noncoding mutation:{non_modified_non_frameshift} reads\n\t"
                f"In-frame mutation:{modified_non_frameshift} reads\n\t"
                f"Frameshift mutation:{modified_frameshift} reads\n"
            )

        with open(_jp("Splice_sites_analysis.txt"), "wt", encoding="utf-8") as outfile:
            unmodified = df_needle_alignment.shape[0] - splicing_sites_modified
            outfile.write(
                "Splice sites analysis:\n\t"
                f"Unmodified:{unmodified} reads\n\t"
                f"Potential splice sites modified:{splicing_sites_modified} reads\n"
            )

        save_vector_to_file(
            effect_vector_insertion_noncoding, "effect_vector_insertion_noncoding"
        )
        save_vector_to_file(
            effect_vector_deletion_noncoding, "effect_vector_deletion_noncoding"
        )
        save_vector_to_file(
            effect_vector_mutation_noncoding, "effect_vector_substitution_noncoding"
        )

    save_vector_to_file(effect_vector_insertion, "effect_vector_insertion_NHEJ")
    save_vector_to_file(effect_vector_deletion, "effect_vector_deletion_NHEJ")
    save_vector_to_file(effect_vector_mutation, "effect_vector_substitution_NHEJ")
    save_vector_to_file(effect_vector_combined, "effect_vector_combined")

    save_vector_to_file(
        avg_vector_ins_all, "position_dependent_vector_avg_insertion_size"
    )
    save_vector_to_file(
        avg_vector_del_all, "position_dependent_vector_avg_deletion_size"
    )

    df_indels = pd.DataFrame(
        np.vstack([hlengths, hdensity]).T, columns=["indel_size", "fq"]
    )
    df_indels.to_csv(_jp("indel_histogram.txt"), index=None, sep="\t")

    df_insertion = pd.DataFrame(
        np.vstack([x_bins_ins[:-1], y_values_ins]).T, columns=["ins_size", "fq"]
    )
    df_insertion.to_csv(_jp("insertion_histogram.txt"), index=None, sep="\t")

    df_deletion = pd.DataFrame(
        np.vstack([-x_bins_del[:-1], y_values_del]).T, columns=["del_size", "fq"]
    )
    df_deletion.to_csv(_jp("deletion_histogram.txt"), index=None, sep="\t")

    df_substitution = pd.DataFrame(
        np.vstack([x_bins_mut[:-1], y_values_mut]).T, columns=["sub_size", "fq"]
    )
    df_substitution.to_csv(_jp("substitution_histogram.txt"), index=None, sep="\t")

    if args.expected_hdr_amplicon_seq:
        save_vector_to_file(
            effect_vector_insertion_mixed, "effect_vector_insertion_mixed_hdr_nhej"
        )
        save_vector_to_file(
            effect_vector_deletion_mixed, "effect_vector_deletion_mixed_hdr_nhej"
        )
        save_vector_to_file(
            effect_vector_mutation_mixed,
            "effect_vector_substitution_mixed_hdr_nhej",
        )
        save_vector_to_file(effect_vector_insertion_hdr, "effect_vector_insertion_HDR")
        save_vector_to_file(effect_vector_deletion_hdr, "effect_vector_deletion_HDR")
        save_vector_to_file(
            effect_vector_mutation_hdr, "effect_vector_substitution_HDR"
        )

    if cut_points:
        pickle.dump(sg_rna_intervals, open(_jp("sg_rna_intervals.pickle"), "wb"))

    if sg_rna_intervals:
        pickle.dump(cut_points, open(_jp("cut_points.pickle"), "wb"))

    if offset_plots.any():
        pickle.dump(offset_plots, open(_jp("offset_plots.pickle"), "wb"))

    if args.dump:
        info("Dumping all the processed data...")
        np.savez(_jp("effect_vector_insertion_NHEJ"), effect_vector_insertion)
        np.savez(_jp("effect_vector_deletion_NHEJ"), effect_vector_deletion)
        np.savez(_jp("effect_vector_substitution_NHEJ"), effect_vector_mutation)

        np.savez(_jp("effect_vector_combined"), effect_vector_combined)

        np.savez(
            _jp("position_dependent_vector_avg_insertion_size"), avg_vector_ins_all
        )
        np.savez(_jp("position_dependent_vector_avg_deletion_size"), avg_vector_del_all)

        df_needle_alignment.to_pickle(_jp("processed_reads_dataframe.pickle"))

        if args.expected_hdr_amplicon_seq:
            np.savez(
                _jp("effect_vector_insertion_mixed_hdr_nhej"),
                effect_vector_insertion_mixed,
            )
            np.savez(
                _jp("effect_vector_deletion_mixed_hdr_nhej"),
                effect_vector_deletion_mixed,
            )
            np.savez(
                _jp("effect_vector_substitution_mixed_hdr_nhej"),
                effect_vector_mutation_mixed,
            )

            np.savez(_jp("effect_vector_insertion_HDR"), effect_vector_insertion_hdr)
            np.savez(_jp("effect_vector_deletion_HDR"), effect_vector_deletion_hdr)
            np.savez(_jp("effect_vector_substitution_HDR"), effect_vector_mutation_hdr)

    info("All Done!")
    print(
        r"""
                )
                (
            __)__
            C\|     |
            \     /
            \___/
            """
    )

    return (
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
    )


def parse_args(args):
    """Parse the arguments passed"""
    parser = argparse.ArgumentParser(
        description="CRISPResso Parameters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-r1",
        "--fastq_r1",
        type=str,
        help="First fastq file",
        required=True,
        default="Fastq filename",
    )

    parser.add_argument(
        "-r2",
        "--fastq_r2",
        type=str,
        help="Second fastq file for paired end reads",
        default="",
    )

    parser.add_argument(
        "-a", "--amplicon_seq", type=str, help="Amplicon Sequence", required=True
    )

    # optional
    parser.add_argument(
        "-g",
        "--guide_seq",
        help=(
            "sgRNA sequence, if more than one, please "
            "separate by comma/s. Note that the sgRNA "
            "needs to be input as the guide RNA sequence "
            "(usually 20 nt) immediately adjacent to but "
            "not including the PAM sequence (5' of NGG for"
            " SpCas9). If the PAM is found on the opposite "
            "strand with respect to the Amplicon Sequence, "
            "ensure the sgRNA sequence is also found on the "
            "opposite strand. The CRISPResso convention is "
            "to depict the expected cleavage position using "
            "the value of the parameter cleavage_offset nt  3' "
            "from the end of the guide. In addition, the "
            "use of alternate nucleases to SpCas9 is "
            "supported. For example, if using the Cpf1 "
            "system, enter the sequence (usually 20 nt) "
            "immediately 3' of the PAM sequence and explicitly "
            "set the cleavage_offset parameter to 1, since "
            "the default setting of -3 is suitable only for SpCas9."
        ),
        default="",
    )
    parser.add_argument(
        "-e",
        "--expected_hdr_amplicon_seq",
        help="Amplicon sequence expected after HDR",
        default="",
    )
    parser.add_argument(
        "-d",
        "--donor_seq",
        help=(
            "Donor Sequence. This optional input comprises "
            "a subsequence of the expected HDR amplicon to "
            "be highlighted in plots."
        ),
        default="",
    )
    parser.add_argument(
        "-c",
        "--coding_seq",
        help=(
            "Subsequence/s of the amplicon sequence "
            "covering one or more coding sequences "
            "for the frameshift analysis.If more than "
            "one (for example, split by intron/s), "
            "please separate by comma."
        ),
        default="",
    )
    parser.add_argument(
        "-q",
        "--min_average_read_quality",
        type=int,
        help="Minimum average quality score (phred33) to keep a read",
        default=0,
    )
    parser.add_argument(
        "-s",
        "--min_single_bp_quality",
        type=int,
        help="Minimum single bp score (phred33) to keep a read",
        default=0,
    )
    parser.add_argument(
        "--min_identity_score",
        type=float,
        help="Minimum identity score for the alignment",
        default=60.0,
    )
    parser.add_argument("-n", "--name", help="Output name", default="")
    parser.add_argument("-o", "--output_folder", help="", default="")
    parser.add_argument(
        "--split_paired_end",
        help=(
            "Splits a single fastq file containing paired "
            "end reads in two files before running CRISPResso"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--trim_sequences",
        help="Enable the trimming of Illumina adapters with Trimmomatic",
        action="store_true",
    )
    parser.add_argument(
        "--trimmomatic_options_string",
        type=str,
        help="Override options for Trimmomatic",
        default=f" ILLUMINACLIP:{get_data('NexteraPE-PE.fa')}"
        ":0:90:10:0:true MINLEN:40",
    )
    parser.add_argument(
        "--min_paired_end_reads_overlap",
        type=int,
        help=(
            "Parameter for the FLASH read merging step. "
            "Minimum required overlap length between two "
            "reads to provide a confident overlap. "
        ),
        default=4,
    )
    parser.add_argument(
        "--max_paired_end_reads_overlap",
        type=int,
        help=(
            "Parameter for the FLASH merging step. "
            "Maximum overlap length expected in "
            "approximately 90%% of read pairs. "
            "Please see the FLASH manual for more information."
        ),
        default=100,
    )
    parser.add_argument(
        "--hide_mutations_outside_window_NHEJ",
        help=(
            "This parameter allows to visualize "
            "only the mutations overlapping the "
            "cleavage site and used to classify "
            "a read as NHEJ. This parameter has "
            "no effect on the quanitification "
            "of the NHEJ. It  may be helpful to "
            "mask a pre-existing and known mutations "
            "or sequencing errors outside the window "
            "used for quantification of NHEJ events."
        ),
        action="store_true",
    )
    parser.add_argument(
        "-w",
        "--window_around_sgrna",
        type=int,
        help=(
            "Window(s) in bp around the cleavage "
            "position (half on on each side) as "
            "determined by the provide guide RNA "
            "sequence to quantify the indels. "
            "Any indels outside this window are "
            "excluded. A value of 0 disables this filter."
        ),
        default=1,
    )
    parser.add_argument(
        "--cleavage_offset",
        type=int,
        help=(
            "Cleavage offset to use within respect "
            "to the 3' end of the provided sgRNA "
            "sequence. Remember that the sgRNA "
            "sequence must be entered without the PAM. "
            "The default is -3 and is suitable for the "
            "SpCas9 system. For alternate nucleases, "
            "other cleavage offsets may be appropriate,"
            "for example, if using Cpf1 this parameter "
            "would be set to 1."
        ),
        default=-3,
    )
    parser.add_argument(
        "--exclude_bp_from_left",
        type=int,
        help=(
            "Exclude bp from the left side of "
            "the amplicon sequence for the "
            "quantification of the indels"
        ),
        default=15,
    )
    parser.add_argument(
        "--exclude_bp_from_right",
        type=int,
        help=(
            "Exclude bp from the right side of "
            "the amplicon sequence for the "
            "quantification of the indels"
        ),
        default=15,
    )
    parser.add_argument(
        "--hdr_perfect_alignment_threshold",
        type=float,
        help="Sequence homology %% for an HDR occurrence",
        default=98.0,
    )
    parser.add_argument(
        "--ignore_substitutions",
        help="Ignore substitutions events for the quantification and visualization",
        action="store_true",
    )
    parser.add_argument(
        "--ignore_insertions",
        help="Ignore insertions events for the quantification and visualization",
        action="store_true",
    )
    parser.add_argument(
        "--ignore_deletions",
        help="Ignore deletions events for the quantification and visualization",
        action="store_true",
    )
    parser.add_argument(
        "--needle_options_string",
        type=str,
        help="Override options for the Needle aligner",
        default="-gapopen=10 -gapextend=0.5  -awidth3=5000",
    )
    parser.add_argument(
        "--keep_intermediate",
        help="Keep all the  intermediate files",
        action="store_true",
    )
    parser.add_argument(
        "--dump",
        help="Dump numpy arrays and pandas dataframes "
        "to file for debugging purposes",
        action="store_true",
    )
    parser.add_argument(
        "--save_also_png",
        help="Save also .png images additionally to .pdf files",
        action="store_true",
    )
    parser.add_argument(
        "-p",
        "--n_processes",
        type=int,
        help=(
            "Specify the number of processes to use for the quantification."
            "Please use with caution since increasing this parameter "
            "will increase significantly the memory required to run CRISPResso."
        ),
        default=1,
    )
    parser.add_argument(
        "--offset_around_cut_to_plot",
        type=int,
        help=(
            "Offset to use to summarize alleles around the cut"
            " site in the alleles table plot."
        ),
        default=20,
    )
    parser.add_argument(
        "--min_frequency_alleles_around_cut_to_plot",
        type=float,
        help="Minimum %% reads required to report an allele in the alleles table plot.",
        default=0.2,
    )
    parser.add_argument(
        "--max_rows_alleles_around_cut_to_plot",
        type=int,
        help="Maximum number of rows to report in the alleles table plot. ",
        default=50,
    )
    parser.add_argument(
        "--debug", action="store_true", help="Print stack trace on error."
    )

    return parser.parse_args()


def main():  # pragma: no cover
    """Main program"""

    def print_stacktrace_if_debug():
        """Print the stack track if in debug mode"""
        debug_flag = False
        if "args" in globals() and "debug" in args:
            debug_flag = args.debug

        if debug_flag:
            traceback.print_exc(file=sys.stdout)

    try:

        print("  \n~~~CRISPResso~~~")
        print("-Analysis of CRISPR/Cas9 outcomes from deep sequencing data-")
        print(
            r"""
                      )
                     (
                    __)__
                 C\|     |
                   \     /
                    \___/
             """
        )
        print(
            "\n[Luca Pinello 2015, send bugs, suggestions or "
            "*green coffee* to lucapinello AT gmail DOT com]\n\n"
        )

        print(f"Version {__version__}\n")

        args = parse_args(sys.argv[1:])

        print(args)

        run_crispresso(args)

    except NTException as exc:
        print_stacktrace_if_debug()
        error(f"Alphabet error, please check your input.\n\nERROR: {exc}")
        sys.exit(1)
    except SgRNASequenceException as exc:
        print_stacktrace_if_debug()
        error(f"sgRNA error, please check your input.\n\nERROR: {exc}")
        sys.exit(2)
    except DonorSequenceException as exc:
        print_stacktrace_if_debug()
        error(
            "Problem with the expected hdr amplicon sequence parameter, "
            f"please check your input.\n\nERROR: {exc}"
        )
        sys.exit(3)
    except TrimmomaticException as exc:
        print_stacktrace_if_debug()
        error(f"Trimming error, please check your input.\n\nERROR: {exc}")
        sys.exit(4)
    except FlashException as exc:
        print_stacktrace_if_debug()
        error(f"Merging error, please check your input.\n\nERROR: {exc}")
        sys.exit(5)
    except NeedleException as exc:
        print_stacktrace_if_debug()
        error(f"Alignment error, please check your input.\n\nERROR: {exc}")
        sys.exit(6)
    except NoReadsAlignedException as exc:
        print_stacktrace_if_debug()
        error(f"Alignment error, please check your input.\n\nERROR: {exc}")
        sys.exit(7)
    except AmpliconEqualDonorException as exc:
        print_stacktrace_if_debug()
        error(
            "Problem with the expected hdr amplicon sequence parameter, "
            f"please check your input.\n\nERROR: {exc}"
        )
        sys.exit(8)
    except CoreDonorSequenceNotContainedException as exc:
        print_stacktrace_if_debug()
        error(f"Donor sequence error, please check your input.\n\nERROR: {exc}")
        sys.exit(9)
    except CoreDonorSequenceNotUniqueException as exc:
        print_stacktrace_if_debug()
        error(f"Donor sequence error, please check your input.\n\nERROR: {exc}")
        sys.exit(10)
    except ExonSequenceException as exc:
        print_stacktrace_if_debug()
        error(f"Coding sequence error, please check your input.\n\nERROR: {exc}")
        sys.exit(11)
    except DuplicateSequenceIdException as exc:
        print_stacktrace_if_debug()
        error(f"Fastq file error, please check your input.\n\nERROR: {exc}")
        sys.exit(12)
    except NoReadsAfterQualityFiltering as exc:
        print_stacktrace_if_debug()
        error(f"Filtering error, please check your input.\n\nERROR:{exc}")
        sys.exit(13)
    except Exception as exc:
        print_stacktrace_if_debug()
        error(f"Unexpected error, please check your input.\n\nERROR: {exc}")
        sys.exit(-1)
