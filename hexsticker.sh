#!/bin/sh
#
# SPDX-FileCopyrightText: 2024 Marcel Schilling <mschilling@idibell.cat>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Hex-sticker generation & optimization  shell script for `scr4eam` scRNA-seq
# simulation pipeline.
#
# Copyright (C) 2024  Marcel Schilling
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

# File:    hexsticker.sh
# Created: 2024-12-12
# Author:  Marcel Schilling <mschilling@idibell.cat>
# License: GNU Affero General Public License Version >= 3.0 (GNU AGPL v3+)
# Purpose: Generate a hex-sticker for the `scr4eam` scRNA-seq read simulation
#          pipeline and optimize the resulting SVG file.


#####################################
# Changelog (reverse chronological) #
#####################################

# 2024-12-12: Initial working version with hard-coded parameters.


########
# Note #
########

# Running this script is not required to use `scr4eam`.
# It is solely provided to allow re-creating the hex-sticker shipped with
# `scr4eam`.
# Thus, this script is neither parameterized nor documented beyond the comments
# in the script itself.


##########################
# Hex-sticker generation #
##########################

# Generate the initial SVG file with the hex-sticker using the corresponding R
# script.
"$(dirname "$(realpath "$0")")"/hexsticker.R


####################
# SVG optimization #
####################

# Optimize the SVG file generated by R.
# Note: Running this script is not required to use `scr4eam`.
# Thus, its dependencies won't be distributed with / required by `scr4eam`.
# Developers will have to manually provide these extra dependencies.
svgo scr4eam.svg -o scr4eam.min.svg


###########
# Cleanup #
###########

# Overwrite the non-optimized SVG file with its optimized version leaving only
# that one file behind.
mv scr4eam.min.svg scr4eam.svg