<!--
SPDX-FileCopyrightText: 2024, 2025 Marcel Schilling <mschilling@ub.edu>

SPDX-License-Identifier: FSFAP
-->

[#]: # (README Markdown file for `scr4eam` scRNA-seq read simulation pipeline.)
[#]: #
[#]: # (Copyright [C] 2024, 2025  Marcel Schilling <mschilling@ub.edu>)
[#]: #
[#]: # (This file is part of `scr4eam`.)
[#]: #
[#]: # (`scr4eam` is free software: you can redistribute it and/or)
[#]: # (modify it under the terms of the GNU Affero General Public)
[#]: # (License as published by the Free Software Foundation, either)
[#]: # (version 3 of the License, or [at your option] any later)
[#]: # (version.)
[#]: #
[#]: # (This program is distributed in the hope that it will be useful,)
[#]: # (but WITHOUT ANY WARRANTY; without even the implied warranty of)
[#]: # (MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the)
[#]: # (GNU Affero General Public License for more details.)
[#]: #
[#]: # (You should have received a copy of the GNU Affero General Public)
[#]: # (License along with this program.  If not, see)
[#]: # (<http://www.gnu.org/licenses/>.)
[#]: #
[#]: # (Copying and distribution of this file, with or without)
[#]: # (modification, are permitted in any medium without royalty)
[#]: # (provided the copyright notice and this notice are preserved.)
[#]: # (This file is offered as-is, without any warranty.)


[#]: # (#####################)
[#]: # ( General Information )
[#]: # (#####################)

[#]: # (File:     README.md)
[#]: # (Created:  2024-12-12)
[#]: # (Modified: 2025-07-11)
[#]: # (Author:   Marcel Schilling <mschilling@ub.edu>)
[#]: # (License:  GNU All-Permissive License)
[#]: # (Purpose:  Document `scr4eam` script collection.)


[#]: # (###################################)
[#]: # ( Changelog [reverse chronological] )
[#]: # (###################################)

[#]: # (2025-07-11:)
[#]: # ( * Email: IDIBELL -> UB.)
[#]: # (2024-12-12:)
[#]: # ( * Initial version:)
[#]: # (   * Header comments.)
[#]: # (   * Summary and Table of Contents.)
[#]: # (   * Description, License, Installation, and Citation sections.)


[#]: # (######)
[#]: # ( Logo )
[#]: # (######)

![scr⁴eam logo](scr4eam.svg)


[#]: # (########)
[#]: # ( README )
[#]: # (########)

# scr⁴eam - <u>*s*</u>ingle-<u>*c*</u>ell <u>*R*</u>NA-seq <u>*r*</u>ealistic <u>*r*</u>andom <u>*r*</u>ead <u>*e*</u>mitting <u>*A*</u>wk <u>*m*</u>ess

scr⁴eam is is a small collection of scripts (that prominently
feature the use of the Awk programming language) to generate realistic
scRNA-seq reads based on an (isoform-level) DGE, read and fragment mapping
statistic of a reference dataset, and a FASTA file with spliced transcript
sequences.

---
## Table of contents

* [License](#license)
* [Usage](#usage)
* [Installation](#installation)
  * [Manual](#manual)
  * [Using conda](#using-conda)
* [Citation](#citation)
* [Documentation of individual
  scripts](#documentation-of-individual-scripts)
  * [`extract_distributions.awk`](#extract_distributions-awk)
    * [Usage](#usage-1)
    * [Command line arguments](#command-line-arguments)
  * [`sample_reads.R`](#sample_reads-r)
    * [Usage](#usage-2)
    * [Command line arguments](#command-line-arguments-1)

---


## License

Copyright © 2024 [Marcel Schilling](mailto:mschilling@idibell)

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received [a copy of the GNU Affero General Public
License](LICENSES/AGPL-3.0-or-later.txt) along with this program.
If not, see <http://www.gnu.org/licenses/>.


## Usage

The CLI is still work in progress and will be added soon.

However, the scripts themselves are already documented
[below](#documentation-of-individual-scripts).


## Installation


### Manual

Detailed installation instructions are still work in progress and will be
added soon.

For now, simply download the scripts and try to run them.
Error messages should help you identify missing dependencies.


### Using conda

The conda package is still work in progress and will be added soon.

For now, follow the manual installation instructions [above](#manual).


## Citation

The scr⁴eam paper is still a work in progress and will be published as a
pre-print soon.

Meanwhile, if you use scr⁴eam, please reference this repository and consider
citing the following pre-print:

> **Quantification of transcript isoforms at the single-cell level using SCALPEL**  
> Franz Ake, Sandra M. Fernández-Moya, Marcel Schilling, Akshay Jaya Ganesh, Ana
> Gutiérrez-Franco, Lei Li, Mireya Plass  
> bioRxiv 2024.06.21.600022; doi: https://doi.org/10.1101/2024.06.21.600022


## Documentation of individual scripts

### `extract_distributions.awk`

Extract distributions required for read generation using scr4eam from read and
fragment mapping summary TSV generated by [SCALPEL](scalpel).


#### Usage

```sh
[gawk --file=]extract_distributions.awk \
  [--assign=multi_isoform_genes_txt=<multi-isoform-genes-txt>] \
  [--assign=read_dists_anchors_skipped_txt=<read-dists-anchors-skipped-txt>] \
  [--assign=n_reads_per_fragment_txt=<n-reads-per-fragment-txt>] \
  [--assign=read_distance_col=<read-distance-column>] \
  [--assign=gene_col=<gene-column>] \
  [--assign=isoform_col=<isoform-column>] \
  [--assign=fragment_col=<fragment-column>] \
  [--assign=read_col=<read-column>] \
  [--assign=fragment_distance_col=<fragment-distance-column>] \
  [< <fragments-and-reads-stats-tsv>] \
  [> <unambiguous-fragment-dists-per-gene-tsv>]
```


#### Command line arguments

* `<fragments-and-reads-stats-tsv>`: TSV file/stream to read read and fragment
                                     mapping summary data (*e.g.* as generated
                                     by [SCALPEL](scalpel)) from; default:
                                     `STDIN`
* `<unambiguous-fragment-dists-per-gene-tsv>`: TSV file/stream to write the
                                               assigned gene and the fragment's
                                               distance (0-based) to the
                                               transcript's (3') end for each
                                               unambiguously assigned fragment
                                               to; default: `STDOUT`
* `multi_isoform_genes_txt`: Text file to write genes expressing more than just
                             a single isoform to (one gene per line); default:
                             `genes.multiple_expressed_isoforms.txt`
* `read_dists_anchors_skipped_txt`: Text file to write distances (0-based) to
                                    the fragment's 3' end for all non-anchor
                                    reads (*i.e.* excluding the first
                                    zero-distance read of each fragment) to
                                    (one integer per line); default:
                                    `read_dists.anchors_skipped.txt`
* `n_reads_per_fragment_txt`: File name/path to write total read counts for
                              each fragment to; default:
                              `n_reads.per_fragment.txt`
* `read_distance_col`: Column index (1-based) in the input TSV containing the
                       read's distance (0-based) from the transcript's (3')
                       end; default: `9`
* `gene_col`: Column index (1-based) in the input TSV specifying the gene the
              read was assigned to; default: `11`
* `isoform_col`: Column index (1-based) in the input TSV specifying the isoform
                 the read was assigned to; default: `12`
* `fragment_col`: Column index (1-based) in the input TSV specifying the
                  fragment the read was assigned to; default: `13`
* `read_col`: Column index (1-based) in the input TSV identifying the read the
              input line belongs to; default: `14`
* `fragment_distance_col`: Column index (1-based) in the input TSV containing
                           the fragment's distance (0-based) from the
                           transcript's (3') end; default: `15`


### `sample_reads.R`

Sample (relative) read coordinates for `scr4eam` based on given iDGE and
read/fragment distributions extracted from reference scRNA-seq data.


#### Usage

```sh
[Rscript ]sample_reads.R [-h/--help] \
  [-f/--fragment-dists \
    <unambiguous-fragment-dists-single-isoform-expressed-genes-txt>] \
  [-n/--reads-per-fragment <n-reads-per-fragment-txt>] \
  [-r/--read-distances <read-dists-anchors-skipped-txt>] \
  [-o/--simulated-reads <simulated-reads-tsv>] \
  <simulated-idge-rds>
```


#### Command line arguments

##### Positional

* `<simulated-idge-rds>`: R object (rds) file to read (sparse) Matrix object
                          holding the simulated iDGE (cells in columns,
                          transcripts in rows) from

##### Optional

* `--help`/`-h`: Show help message and exit
* `--fragment-dists`/`-f`: (Compressed) text file to read the fragments'
                           (0-based) distances to their corresponding
                           transcript's (3') ends from (one integer per line);
                           default:
                           `unambiguous_fragment_dists.single_isoform_expressed_genes.txt.gz`
* `--reads-per-fragment`/`-n`: (Compressed) text file to read the read counts
                               per fragment from (one integer per line);
                               default: `n_reads.per_fragment.txt`
* `--read-distances`/`-r`: (Compressed) text file to read the reads'
                           distances (0-based) to their corresponding
                           fragment's 3' end for non-anchor reads (*i.e.*
                           excluding the first zero- distance read of the
                           fragment) from (one integer per line); default:
                           `read_dists.anchors_skipped.txt`
* `--simulated-reads`/`-o`: (Compressed) TSV file to write synthetic read data
                            (columns: transcript name/ID (from input iDGE row
                            names), read ID (cell ID (from input iDGE column
                            names), fragment ID (`fragment_N` with `N`
                            increasing (1-based)), and read ID (`read_M` with
                            `M` increasing (1-based), concatenated with
                            underscores (`_`) as separators)), and (relative)
                            distance (0-based) of the read to the corresponding
                            transcript(!)'s (3') end) to; default:
                            `simulated_reads.tsv.gz`


[#]: # (#######)
[#]: # ( Links )
[#]: # (#######)

[scalpel]: "https://github.com/p-CMRC-LAB/SCALPEL" "SCALPEL GitHub page"
