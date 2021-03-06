# bamwindow

`bamwindow` segments aligned `BAM` files for downstream copy number analysis. Rather than using adaptive binning, `bamwindow` splits the genome into bins of fixed size, thus allowing multiple samples to be directly compared. Windows are generated using 1-based coordinated. Both the start and end nucleotide of a region are included. If a region is supplied, binning is only performed on that specified region; otherwise all chromosomes in the BAM file are processed. Unless the window size is an exact fraction of the sequence (or region) length, the final window will be shorted than the specified width.

## Window Matching

There are three possible modes for matching a read to a window, overlap, midpoint and start:

* In overlap mode (`-m0` also the default), a read is assigned to a window if there is any overlap between the read and the window. This allows for a read to be counted in multiple windows.
* In start mode (`-m2`), a read is assigned to a single window based upon its start position.
* In midpoint mode (`-m2`), a read is assigned to a single window by its midpoint (rounded towards the chromsome start).

In start and midpoint modes, a read is guaranteed to be counted at most once. Read clipped regions are discarded before window assignment.

## Licence

bamwindow is released under the [GNU General Public License version 3](http://www.gnu.org/licenses/gpl.html).

The code in the `cmake` directory for finding HTSlib is taken from bamdb by Aaron Dornbrand-Lo. (see https://github.com/D-Lo/bamdb). These files are released under the BSD 2-Clause "Simplified" License.

## Options & Arguments

bamwindow is called as:

    bamwindow [-hvem] n file [region]

* `-h`: print help
* `-v`: print program and `HTSLib` version information
* `-e`: include all regions in the output, even if they are empty
* `-m`: the read mode (0, 1, or 2). See above
* `n `: the window size 
* `file`: the input BAM file to process
* `region`: optionally, a chromosome region to process

## Output Format

Output is returned in a tab-delimited format, containing the following columns:

* window chromosome name
* window start
* window end
* read count

If the `-e` argument was supplied, all processed windows will be returned in order, allowing for direct comparison with other samples. Otherwise, empty windows are not returned.

## Building & Installation

* bamwindow requires the [htslib](https://github.com/samtools/htslib) library. This should be present before you attempt to build bamwindow.
* bamwindow uses [cmake](https://cmake.org) to simplify the make process. This must be installed for installation.

The bamwindow code can be built from the main folder with:

~~~bash
mkdir build
cd build
cmake ../
make
~~~

You can subsequently install bamwindow with:

~~~bash
sudo make install
~~~
