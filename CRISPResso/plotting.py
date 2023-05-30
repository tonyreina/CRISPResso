import os

import numpy as np
import pandas as pd
import plotly.express as px

import matplotlib.pyplot as plt


def plot1_indel_size(hlengths, hdensity, center_index, args):
    """Plot figure 1

    Args:
        hlengths
        hdensity
        center_index

    """

    _jp = lambda filename: os.path.join(
        args.output_directory, filename
    )  # handy function to put a file in the output directory

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
    fig1a.write_image(_jp( f"{filename_fig1a}.pdf"))
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

        fig2_df = pd.DataFrame(data, columns=["Type", "Reads"])

        fig2 = px.pie(fig2_df, values="Reads", names="Type")

        filename_fig2 = "2.Unmodified_NHEJ_HDR_pie_chart"
        fig2.write_image(_jp(f"{filename_fig2}.pdf"))
        fig2.write_html(_jp(f"{filename_fig2}.html"))
        if args.save_also_png:
            fig2.write_image(_jp(f"{filename_fig2}.png"))

    else:
        plt.figure(figsize=(12 * 1.5, 14.5 * 1.5))
        ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=5)
        _, texts, autotexts = ax1.pie(
            [n_unmodified / n_total * 100, n_modified / n_total * 100],
            labels=[
                f"Unmodified\n({n_unmodified} reads)",
                f"NHEJ\n({n_modified} reads)",
            ],
            explode=(0, 0),
            colors=[(1, 0, 0, 0.2), (0, 0, 1, 0.2)],
            autopct="%1.1f%%",
        )

        if cut_points:
            ax2 = plt.subplot2grid((6, 3), (5, 0), colspan=3, rowspan=1)
            ax2.plot(
                [0, length_amplicon], [0, 0], "-k", lw=2, label="Amplicon sequence"
            )

            for idx, sg_rna_int in enumerate(sg_rna_intervals):
                if idx == 0:
                    ax2.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="sgRNA",
                        solid_capstyle="butt",
                    )
                else:
                    ax2.plot(
                        [sg_rna_int[0], sg_rna_int[1]],
                        [0, 0],
                        lw=10,
                        c=(0, 0, 0, 0.15),
                        label="_nolegend_",
                        solid_capstyle="butt",
                    )

            ax2.plot(
                cut_points + offset_plots,
                np.zeros(len(cut_points)),
                "vr",
                ms=12,
                label="Predicted Cas9 cleavage site/s",
            )
            lgd = plt.legend(
                bbox_to_anchor=(0, 0, 1.0, 0),
                ncol=1,
                mode="expand",
                borderaxespad=0.0,
                numpoints=1,
                prop={"size": "large"},
            )
            plt.xlim(0, length_amplicon)
            plt.axis("off")

        proptease = fm.FontProperties()
        proptease.set_size("xx-large")
        plt.setp(autotexts, fontproperties=proptease)
        plt.setp(texts, fontproperties=proptease)
        plt.savefig(
            _jp("2.Unmodified_NHEJ_pie_chart.pdf"),
            pad_inches=1,
            bbox_inches="tight",
        )
        if args.save_also_png:
            plt.savefig(
                _jp("2.Unmodified_NHEJ_pie_chart.png"),
                pad_inches=1,
                bbox_inches="tight",
            )

    pdf.attach_note("Unmodified NEHJ pie chart")
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


def plot3(df_needle_alignment, n_total, pdf, args):
    """Plot 3"""

    # #############################################
    # (3) a graph of frequency of deletions and
    # insertions of various sizes (deletions
    # could be consider as negative numbers and insertions as positive);

    fig = plt.figure(figsize=(26, 6.5))

    ax = fig.add_subplot(1, 3, 1)
    ax.bar(x_bins_ins[:-1], y_values_ins, align="center", linewidth=0, color=(0, 0, 1))
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
    ax.bar(-x_bins_del[:-1], y_values_del, align="center", linewidth=0, color=(0, 0, 1))
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
    ax.bar(x_bins_mut[:-1], y_values_mut, align="center", linewidth=0, color=(0, 0, 1))
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
        offset_plots,
        y_values_mut,
        x_bins_mut,
        y_values_ins,
        x_bins_ins,
        y_values_del,
        x_bins_del,
    )


def plot4a(
    offset_plots,
    effect_vector_any,
    sg_rna_intervals,
    cut_points,
    length_amplicon,
    pdf,
    args,
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
    offset_plots,
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


def plot4c(
    offset_plots,
    effect_vector_insertion_hdr,
    effect_vector_deletion_hdr,
    effect_vector_mutation_hdr,
    LEN_AMPLICON,
    n_total,
    n_modified,
    n_repaired,
    cut_points,
    sg_rna_intervals,
    pdf,
    args,
):
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


def plot4d(
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
):
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


def plot4e(
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
):
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
    plt.close()


def plot5(
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
):
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
    plt.close()


def plot6():
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
    plt.close()


def plot7():
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


def plot8():
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
    plt.close()
