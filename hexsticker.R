#!/usr/bin/env Rscript
#
# SPDX-FileCopyrightText: 2024, 2025 Marcel Schilling <mschilling@ub.edu>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Hex-sticker generation R script for `scr4eam` scRNA-seq simulation pipeline.
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

# File:     hexsticker.R
# Created:  2024-12-12
# Modified: 2025-07-11
# Author:   Marcel Schilling <mschilling@ub.edu>
# License:  GNU Affero General Public License Version >= 3.0 (GNU AGPL v3+)
# Purpose:  Generate a hex-sticker for the `scr4eam` scRNA-seq read simulation
#           pipeline.


#####################################
# Changelog (reverse chronological) #
#####################################

# 2025-07-11: Email: IDIBELL -> UB.
# 2024-12-12: Initial working version with hard-coded parameters.


########
# Note #
########

# Running this script is not required to use `scr4eam`.
# It is solely provided to allow re-creating the hex-sticker shipped with
# `scr4eam`.
# Thus, this script is neither parameterized nor documented beyond the comments
# in the script itself.


#############
# Libraries #
#############

# Note: Running this script is not required to use `scr4eam`.
# Thus, its dependencies won't be distributed with / required by `scr4eam`.
# Developers will have to manually provide these extra dependencies.
library(magrittr)
library(tibble)
library(dplyr)
library(Gviz)
library(svglite)
library(svgparser)
library(grid)
library(emoGG)
library(hexSticker)


######################
# Read coverage plot #
######################

# Create SVG file with read coverage plot cartoon showing multi-isoform gene
# with reads in different shades of grey (representing different cell barcodes)
# distributed along the 3' ends of the isoforms including a spliced read to
# symbolize the key features of `scr4eam`.
set.seed(37L)
svglite("reads_plot.svg", bg = "transparent", height = 2, width = 4)
data.frame(chromosome = "chr",
           start = rep(c(1000L, 1737L, 1800L), 2L),
           strand = "+",
           feature = rep(c(rep("exon", 2L), "utr3"), 2L),
           gene = "gene",
           exon = paste0("ex", c(1:3, 1:2, 4L)),
           transcript = rep(paste0("tx", 1:2), each = 3)) %>%
  mutate(end = c(1200L, 1799L, 2300L, 1200L, 1799L, 2000L),
         width = end - start + 1L,
         symbol = transcript,
         across(where(is.character), as.factor)) %>%
  relocate(chromosome, start, end, width, -symbol) %>%
  GeneRegionTrack(fill = "grey20", col.line = "grey20", alpha = .65) %>%
  list(AnnotationTrack(chromosome = "chr",
                       start = c(1173L, 1737L,
                                 as.integer(rnorm(6L, 1837, 37)),
                                 as.integer(rnorm(7L, 2073, 73))),
                       width = c(28L, 73L - 28L, rep(73L, 13L)),
                       group = paste0("r", c(1L, 1:14)),
                       col.line = "grey40",
                       fill = c(rep("grey40", 2L),
                                sample(c("blacK", paste0("grey", c(20, 40))),
                                       13L, replace = TRUE)),
                       alpha = .95)) %>%
  rev %>%
  plotTracks(showTitle = FALSE, sizes = c(1, pi / 2), col.border.title = NA,
             col = "transparent")
dev.off()


########################
# SVG artifact removal #
########################

# The SVG generated above contains a tiny dot in the bottom left corner that
# visually impacts the hex-sticker.
# The following function scans through all leaves of a graphical object tree
# and removes those that have exactly one set of coordinates with both `x` and
# `y` being equal to zero.
remove_corner <- function(gtree) {
  for(grob.name in grid.ls(gtree, print = FALSE)$name) {
    grob <- getGrob(gtree, grob.name, strict = FALSE)
    if(is.null(grob)) { next }
    if("gTree" %in% class(grob)) { next }
    coords <- grobCoords(grob)
    if(length(coords) != 1L) { next }
    coords %<>% `[[`(1L)
    if(length(coords) != 2L) { next }
    if(!"x" %in% names(coords)) { next }
    if(!"y" %in% names(coords)) { next }
    if(!all(coords$x == 0)) { next }
    if(!all(coords$y == 0)) { next }
    gtree %<>% removeGrob(grob.name)
  }
 return(gtree)
}


#####################################
# Fragment length distribution plot #
#####################################

# Generate a histogram representing the fragment length distribution of
# scRNA-seq data representing the fact that `scr4eam` uses reference data
# statistics to generate realistically distributed reads.
set.seed(73L)
png("scr4eam.png")
rgamma(10000L, pi, 1/73) %>%
  round %>%
  `-`(.73 * mean(.)) %>%
  as.integer %>%
  `[`(. >= 0) %>%
  enframe(name = NULL, value = "fraglen") %>%
  ggplot(aes(fraglen)) %>%
  `+`(add_emoji(emoji_search("munch")[1L, "code"])) %>%
  `+`(geom_histogram(binwidth = 25L, alpha = .73)) %>%
  `+`(xlab("fragment length")) %>%
  `+`(theme_void()) %>%


########################
# Hex-sticker assembly #
########################

  # Inset read coverage plot into fragment length histogram and assemble
  # the final hex-sticker and export it to SVG.
  `+`(annotation_custom(read_svg("reads_plot.svg") %>%
                          remove_corner,
                        xmin = 150, ymin = 400)) %>%
  sticker(package = "scr⁴eam",
          s_x = 1, s_y = .9, s_width = 1.3, s_height = 1,
          p_y = 1.55, p_size = 8, p_color = "gold1",
          h_fill = "mediumpurple3", h_color = "gold1",
          filename = "scr4eam.svg")


###########
# Cleanup #
###########

# The `hexsticker::sticker` function seems to always open a graphical device,
# even though it is writing the sticker to a file anyhow.
# Thus, instead of cluttering the device with `Rplots.pdf` or the likes, the
# `png` device was explicitly opened above (see Fragment length distribution
# plot section) and is closed here before removing the (empty) PNG file created
# in the process.
dev.off()
file.remove("scr4eam.png")

# Cleanup the read coverage plot SVG file that was created and re-imported
# above, leaving only the hex-sticker SVG file behind.
file.remove("reads_plot.svg")
