#!/usr/bin/env Rscript
#
# SPDX-FileCopyrightText: 2024, 2025 Marcel Schilling <mschilling@ub.edu>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Read sampling R script for `scr4eam` scRNA-seq read simulation pipeline.
#
# Copyright (C) 2024, 2025  Marcel Schilling
#
# This file is part of `scr4eam`.
#
# `scr4eam` is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#######################
# General information #
#######################

# File:     sample_reads.R
# Created:  2024-12-12
# Modified: 2025-07-11
# Author:   Marcel Schilling <mschilling@ub.edu>
# License:  GNU Affero General Public License Version >= 3.0 (GNU AGPL v3+)
# Purpose:  Sample (relative) read coordinates for `scr4eam` based on given
#           iDGE and read/fragment distributions extracted from reference
#           scRNA-seq data.


#####################################
# Changelog (reverse chronological) #
#####################################

# 2025-07-11: Email: IDIBELL -> UB.
# 2024-12-12: Initial working version.


##############
# Parameters #
##############

# Set default parameters.
unambiguous_fragment_dists.single_isoform_expressed_genes.txt.default <-
  "unambiguous_fragment_dists.single_isoform_expressed_genes.txt.gz"
n_reads.per_fragment.txt.default <- "n_reads.per_fragment.txt"
read_dists.anchors_skipped.txt.default <- "read_dists.anchors_skipped.txt"
simulated_reads.tsv.default <- "simulated_reads.tsv.gz"


#############
# Libraries #
#############

# Load required libraries without cluttering the user-facing output.
suppressWarnings(
  suppressPackageStartupMessages(invisible(lapply(c("argparse",
                                                    "dplyr",
                                                    "magrittr",
                                                    "ggplot2",
                                                    "MASS",
                                                    "mefa4",
                                                    "tidyr",
                                                    "readr"),
                                                  library,
                                                  character.only = TRUE,
                                                  warn.conflicts = FALSE,
                                                  quietly = TRUE))))


###############
# Description #
###############

# Define description for help message and add to argument parser.
paste("Sample (relative) read coordinates for `scr4eam` based on given",
      "iDGE and read/fragment distributions extracted from reference",
      "scRNA-seq data.") %>%
  ArgumentParser(description = .,
                 prog = paste0("[Rscript ]",
                               argparse:::get_Rscript_filename())) %>%


##########################
# Command line arguments #
##########################

  # Add all parameters to argument parser, document for help message, parse
  # from command line, and attach to environment.
  {
    .$add_argument("simulated_idge.rds", metavar = "<simulated-idge-rds>",
                   help = paste("R object (rds) file to read (sparse) Matrix",
                                "object holding the simulated iDGE (cells in",
                                "columns, transcripts in rows) from"))
    .$add_argument(
        "--fragment-dists", "-f",
        metavar =
          "<unambiguous-fragment-dists-single-isoform-expressed-genes-txt>",
        default =
          unambiguous_fragment_dists.single_isoform_expressed_genes.txt.default,
        help = paste("(Compressed) text file to read the fragments' (0-based)",
                     "distances to their corresponding transcript's (3') ends",
                     "from (one integer per line); default:",
                     paste0(
                       "`",
                       unambiguous_fragment_dists.single_isoform_expressed_genes.txt.default,
                       "`")),
        dest = "unambiguous_fragment_dists.single_isoform_expressed_genes.txt")
    .$add_argument("--reads-per-fragment", "-n",
                   metavar = "<n-reads-per-fragment-txt>",
                   default = n_reads.per_fragment.txt.default,
                   help = paste("(Compressed) text file to read the read",
                                "counts per fragment from (one integer per",
                                "line); default:",
                                paste0("`", n_reads.per_fragment.txt.default,
                                       "`")),
                   dest = "n_reads.per_fragment.txt")
    .$add_argument("--read-distances", "-r",
                   metavar = "<read-dists-anchors-skipped-txt>",
                   default = read_dists.anchors_skipped.txt.default,
                   help = paste("(Compressed) text file to read the reads'",
                                "distances (0-based) to their corresponding",
                                "fragment's 3' end for non-anchor reads (i.e.",
                                "excluding the first zero-distance read of",
                                "the fragment) from (one integer per line);",
                                "default:",
                                paste0( "`",
                                       read_dists.anchors_skipped.txt.default,
                                       "`")),
                   dest = "read_dists.anchors_skipped.txt")
    .$add_argument("--simulated-reads", "-o",
                   metavar = "<simulated-reads-tsv>",
                   default = simulated_reads.tsv.default,
                   help = paste("(Compressed) TSV file to write synthetic",
                                "read data (columns: transcript name/ID (from",
                                "input iDGE row names), read ID (cell ID",
                                "(from input iDGE column names), fragment ID",
                                "(`fragment_N` with `N` increasing",
                                "(1-based)), and read ID (`read_M` with `M`",
                                "increasing (1-based), concatenated with",
                                "underscores (`_`) as separators)), and",
                                "(relative) distance (0-based) of the read to",
                                "the corresponding transcript(!)'s (3') end)",
                                "to; default:",
                                paste0( "`", simulated_reads.tsv.default,
                                       "`")),
                   dest = "simulated_reads.tsv")
    .
  } %$%
  parse_args() %>%
  attach


#################
# Input reading #
#################

# Read (unambiguous, single-isoform-expressed genes) fragment distances from
# (compressed) text file, one line at a time.
unambiguous_fragment_dists.single_isoform_expressed_genes.txt %>%
  readLines %>%
  as.integer ->
    unambiguous_fragment_dists.single_isoform_expressed_genes

# Read read counts per fragment from (compressed) text file, one line at a
# time.
n_reads.per_fragment.txt %>%
  readLines %>%
  as.integer ->
    n_reads.per_fragment

# Read (non-anchor) read distances from (compressed) text file, one line at a
# time.
read_dists.anchors_skipped.txt  %>%
  readLines %>%
  as.integer ->
    read_dists.anchors_skipped

# Read (isoform-level) DGE from RDS file.
simulated_idge.rds %>%
  readRDS ->
    simulated_idge


#######################
# Fragment simulation #
#######################

# Sample fragments with corresponding read counts.
simulated_idge %>%
  Melt %>%
  as_tibble %>%
  setNames(c("transcript", "cell", "n_fragments")) %>%
  group_by(across(-n_fragments)) %>%
  mutate(fragment = list(1L:n_fragments)) %>%
  dplyr::select(-n_fragments) %>%
  ungroup %>%
  unnest(fragment) %>%
  group_by(cell) %>%
  mutate(fragment = 1L:n()) %>%
  ungroup %>%
  unite(fragment, cell:fragment, sep = "_fragment") %>%
  mutate(fragment_dist =
           sample(unambiguous_fragment_dists.single_isoform_expressed_genes,
                  size = n(), replace = TRUE),
         n_reads =
           sample(n_reads.per_fragment, size = n(), replace = TRUE)) ->
    simulated_fragments

# Cleanup (isoform-level) DGE to free memory.
rm(simulated_idge)


###################
# Read simulation #
###################

# Generate anchor read for each fragment and write to (compressed) output TSV
# file.
simulated_fragments %>%
  dplyr::select(-n_reads) %>%
  mutate(across(fragment, \(fragment) paste0(fragment, "_read1"))) %>%
  rename(read = fragment, read_dist = fragment_dist) %>%
  write_tsv(simulated_reads.tsv)

# Discard single-read fragments to free memory.
simulated_fragments %<>% filter(n_reads > 1L)

# Sample non-anchor read for multi-read fragments and append to (compressed)
# output TSV file.
simulated_fragments %>%
  mutate(across(fragment, as.factor)) %>%
  group_by(across(-n_reads)) %>%
  mutate(read = list(2L:n_reads)) %>%
  dplyr::select(-n_reads) %>%
  ungroup %>%
  unnest(read) %>%
  unite(read, fragment, read, sep = "_read") %>%
  mutate(read_dist = sample(read_dists.anchors_skipped, size = n(),
                            replace = TRUE)) %>%
  mutate(across(read_dist, \(read_dist) read_dist + fragment_dist)) %>%
  dplyr::select(-fragment_dist) %>%
  write_tsv(simulated_reads.tsv, append = TRUE)
