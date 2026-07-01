# Polikit 0.4
[![GitHub](https://img.shields.io/badge/GitHub-V0.4-C71D23?logo=github&logoColor=white&labelColor=000)](https://github.com/jiahuuui/polikit/)
[![Bitbucket](https://img.shields.io/badge/Bitbucket-V0.4-0052CC?logo=bitbucket&logoColor=white&labelColor=000)](https://bitbucket.org/jiahuijiahui/polikit/src/master/)

## A polyhedral analysis toolkit

This package is originally developed for polyhedral analysis of amorphous structures. Now it has modules for other analysis methods, including bond angle analysis, RDF analysis, TCT analysis and ring statistics analysis. For now, file formats including xyz, lammps data file, lammps dump file can be read, but please carefully check the format when using the package. Analysis can be performed in either static or dynamic way, depends on whether the analysis only involves one file, or also comparison with other files.

Please contact **zjh239@foxmail.com** if you have bugs or issues to report.

#### 0.5

- Read from `.toml` files.

#### 0.4  *(3 Mar. 2025)*

 - Non-affine displacement analysis.
 - Pair-wise cutoff values.
 - Dimension-wise periodic boundary condition.

#### 0.3  *(22 Jan. 2025)*

 - Ring statistics analysis.

#### 0.2
 - Bond angle distribution (BAD) analysis.
 - Radial distribution function.

#### 0.1
 - Polyhedral analysis.

## How to compile
At the root directory of the code:

1. `mkdir build && cd build`

2. `cmake ../.` or `cmake -DDEBUG=on ../.` for debug mode.

3. `make` or `make -j`

## What Polikit can do

### Analysis of a static configuration:

1. Polyhedral analysis

2. Bond angle distribution

3. Coordination distribution/change

4. Rings statistics distribution

5. Radial distribution function

7. Cluster analysis

### Analysis of a series of configurations:

1. Non-affine displacement (D2min)

2. Localized plastic events analysis (LPSE)

3. Cluster inheritence analysis

4. Topological constraint analysis

5. Neighbor change analysis

6. Polyhedral neighbor change analysis

## Usage
For static analysis(`-f`):

**`./polikit -f abc.xyz -p 1 -poly 2.3`**

`-p [int]` can be 1 or 0, decides whether periodic boundary condition will be applied.

`-[key] [parameters] ...` gives computing options and corresponding parameters. Now the availables are: `poly [cutoff]` - polyhedral analysis, `bad [cutoff]` - bond angle analysis, `rdf [cutoff]` - radial distribution function, `ring [cutoff] [max size]` - ring statistics analysis, `d2min [cutoff]` - non-affine displacement analysis, `cluster [cutoff]` - cluster analysis, `lpes [d2min_cutoff] [cutoff]` - LPSE analysis (combined of d2min and cluster analysis).

For dynamic analysis(`-d`):

**`./polikit -d dumpfiles -os 3 -p 1 -d2min 4.6`**

`-d` -- giving a directory name that contains .xyz files. 20 is the frame invertal for dynamic comparison. Other parameters work in the same way as in static analysis.

`-os [int]` -- frame interval, which is applicable when frame-wise comparison is performed.

`-skip [int]` -- skipping the first N frames, and starting from the N+1 frame. Note that sometimes the frame number may differ from the file name.

## Examples

- Polyhedral analysis

`./polikit -f ../test/ga2o3_test.xyz -p 1 -poly 2.3`

- Bond angle analysis

`./polikit -f ../test/ga2o3_test.xyz -p 1 -bad 2.3`

- Radial distribution

`./polikit -f ../test/ga2o3_test.xyz -p 1 -rdf 10`

- Wendt-Abraham parameter calculation

`./polikit -f ../test/ga2o3_test.xyz -p 1 -wa 5`

- Honeycutt-Anderson parameters analysis

`./polikit -f ../test/ga2o3_test.xyz -p 1 -ha 2.3`

- Ring statistics analysis

`./polikit -f ../test/ga2o3_test.xyz -p 1 -ring 2.3 8`

- Dynamic neighbor change analysis

`./polikit -d ../test/test_dir/ -os 1 -p 1 -nc 2.3`

- D2min analysis

`./polikit -d ../test/test_dir/ -os 2 -p 1 -d2min 4.6`

- D2min analysis and cluster analysis on high D2min atoms / LPSE analysis

`./polikit -d ../test/test_dir/ -os 2 -p 1 -d2min 4.6 -cluster 2.3`

`./polikit -d ../test/test_dir/ -os 2 -p 1 -lpse 4.6 2.3`

- LPSE inheritence analysis

`./polikit -d ../test/test_dir/ -os 1 -p 1 -ci 4.6 2.3`

## Extract results

Lines correspond to specific results are given a specific initial letter or letters. Therefore, those lines can be extracted using `grep`, `awk`, and `csplit` command. Suppose the code is run with a dumping option `> out.log`, below are the ways to extract specific analysis results.

1. Data that only takes one line for each structure.

- Coordination number distribution

`grep 'c|' out.log > cn.txt`

- Polyhedral type distribution

`grep 'p|' out.log > poly.txt`

2. Data that takes a fixed number of lines.

- Bond angle distribution

`grep -A 180 'b2|' out.log > bad.txt`

- D2min analysis

`grep -A 400 'd1|' out.log > d2min.txt`

- Split file

`csplit -f mbad_ -b %d.txt bad.txt '/--/' '{*}' -s -k`

This splits the file on every '--' line.

3. Data that takes uncertain number of lines, with a starting and ending initial respectively.

- Cluster analysis results

`awk '/ c1\|/  {inblock = 1} / c2\|/  {inblock = 0}  inblock == 1'  out.log > cluster.txt`

- Split file

`csplit -f mcluster_ -n 3 cluster.txt '/ c1|/' '{*}' -s -k`
