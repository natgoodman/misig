Mistakes of Significance
================
Nathan (Nat) Goodman
January 16, 2019

<!-- README.md is generated from README.Rmd. Please edit that file -->
*A collection of R scripts and documents exploring mistakes made by significance testing. The only content at present is a blog post entitled "When You Select Significant Findings, You’re Selecting Inflated Estimates" whose main point is that significance testing is a biased procedure that overestimates effect size.*

**THE SOFTWARE IS STILL ROUGH and SOFTWARE DOCUMENTATION NONEXISTENT. PLEASE GET IN TOUCH IF YOU NEED HELP**

Overview
--------

The program explores the bias of significance testing and why it inflates effect size. The software simulates *studies* for conditions of interest. The studies are simple two group comparisons parameterized by sample size *n* and population effect size *d*<sub>*p**o**p*</sub> (*d*<sub>*p**o**p*</sub> ≥ 0). For each study, I generate two groups of random numbers, each of size *n*. One group, *group0*, comes from a standard normal distribution with *m**e**a**n* = 0; the other, *group1*, is from a standard normal distribution with *m**e**a**n* = *d*<sub>*p**o**p*</sub>. The effect size statistic is standardized difference, aka *Cohen's d*, defined as the mean of *group1* minus the mean of *group0* divided by the pooled standard deviation of the two groups. P-values are from the t-test.

> **Technical note**: The pairing of Cohen's d with the t-test is mathematically natural because both are *standardized*, meaning both are relative to the sample standard deviation. In fact, Cohen's d and the t-statistic are essentially the same statistic, related by the identities $d=t\\sqrt{2/n}$ and $t=d\\sqrt{n/2}$ (for my simulation scenario).

The program performs two types of simulations.

1.  *s**i**m*<sub>*r**a**n**d*</sub>. Randomly selects a large number of *d*<sub>*p**o**p*</sub>s and simulates one study per *d*<sub>*p**o**p*</sub>.
2.  *s**i**m*<sub>*f**i**x**d*</sub>. Fixes *d*<sub>*p**o**p*</sub> to a few values of interest, sets *n* to a range of values, and simulates many studies for each *d*<sub>*p**o**p*</sub> and *n*.

-   *s**i**m*<sub>*r**a**n**d*</sub>. Randomly selects a large number of *d*<sub>*p**o**p*</sub>s and simulates one study per *d*<sub>*p**o**p*</sub>.
-   *s**i**m*<sub>*f**i**x**d*</sub>. Fixes *d*<sub>*p**o**p*</sub> to a few values of interest, sets *n* to a range of values, and simulates many studies for each *d*<sub>*p**o**p*</sub> and *n*.

The program uses *s**i**m*<sub>*f**i**x**d*</sub> to estimate the mean significant observed effect size for each *d*<sub>*p**o**p*</sub> and *n*. It also caculates these values analytically.

After running the simulations and calculating the mean significant observed effect sizes, the program generates the figures that appear in the blog post, as well as the specific results mentioned in the post. It stores the latter in a table for future reference.

**Notation**

<table>
<colgroup>
<col width="10%" />
<col width="7%" />
<col width="82%" />
</colgroup>
<thead>
<tr class="header">
<th>in code</th>
<th>in text</th>
<th>meaning</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>n</code></td>
<td><span class="math inline"><em>n</em></span></td>
<td>sample size</td>
</tr>
<tr class="even">
<td><code>d.pop</code></td>
<td><span class="math inline"><em>d</em><sub><em>p</em><em>o</em><em>p</em></sub></span></td>
<td>population effect size; this arises in the text but not the software</td>
</tr>
<tr class="odd">
<td><code>d.sdz</code></td>
<td><span class="math inline"><em>d</em><sub><em>s</em><em>d</em><em>z</em></sub></span></td>
<td>standardized observed effect size, aka Cohen's d</td>
</tr>
<tr class="even">
<td><code>d</code></td>
<td><span class="math inline"><em>d</em></span></td>
<td>variously means either <span class="math inline"><em>d</em><sub><em>p</em><em>o</em><em>p</em></sub></span> or <span class="math inline"><em>d</em><sub><em>s</em><em>d</em><em>z</em></sub></span>; hopefully clear from context</td>
</tr>
<tr class="odd">
<td><code>sim.rand</code></td>
<td><span class="math inline"><em>s</em><em>i</em><em>m</em><sub><em>r</em><em>a</em><em>n</em><em>d</em></sub></span></td>
<td>simulation that randomly selects a large number of <span class="math inline"><em>d</em><sub><em>p</em><em>o</em><em>p</em></sub></span>s and simulates one study per <span class="math inline"><em>d</em><sub><em>p</em><em>o</em><em>p</em></sub></span>.</td>
</tr>
<tr class="even">
<td><code>sim.fixd</code></td>
<td><span class="math inline"><em>s</em><em>i</em><em>m</em><sub><em>f</em><em>i</em><em>x</em><em>d</em></sub></span></td>
<td>simulation that fixes <span class="math inline"><em>d</em><sub><em>p</em><em>o</em><em>p</em></sub></span> to a few values of interest, sets <span class="math inline"><em>n</em></span> to a range of values, and simulates many studies for each <span class="math inline"><em>d</em><sub><em>p</em><em>o</em><em>p</em></sub></span> and <span class="math inline"><em>n</em></span></td>
</tr>
<tr class="odd">
<td><code>meand.empi</code></td>
<td><span class="math inline"><em>m</em><em>e</em><em>a</em><em>n</em><em>d</em><sub><em>e</em><em>m</em><em>p</em><em>i</em></sub></span></td>
<td>mean significant observed effect size estimated from simulation (<span class="math inline"><em>s</em><em>i</em><em>m</em><sub><em>f</em><em>i</em><em>x</em><em>d</em></sub></span>)</td>
</tr>
<tr class="even">
<td><code>meand.theo</code></td>
<td><span class="math inline"><em>m</em><em>e</em><em>a</em><em>n</em><em>d</em><sub><em>t</em><em>h</em><em>e</em><em>o</em></sub></span></td>
<td>mean significant observed effect size calculated analytically</td>
</tr>
<tr class="odd">
<td><code>pval</code></td>
<td></td>
<td>p-value</td>
</tr>
<tr class="even">
<td><code>df</code></td>
<td></td>
<td>degrees of freedom</td>
</tr>
<tr class="odd">
<td><code>t</code></td>
<td><span class="math inline"><em>t</em></span></td>
<td>t-statistic</td>
</tr>
<tr class="even">
<td><code>d0</code></td>
<td></td>
<td>noncentrality parameter expressed in standardized <span class="math inline"><em>d</em></span> units; this is the presumed population effect size</td>
</tr>
<tr class="odd">
<td><code>ncp</code></td>
<td></td>
<td>noncentrality parameter expressed in <span class="math inline"><em>t</em></span> units</td>
</tr>
<tr class="even">
<td><code>sig.level</code></td>
<td><span class="math inline"><em>α</em></span></td>
<td>significance level</td>
</tr>
</tbody>
</table>

Installation and Usage
----------------------

The software is **not a package** and cannot be installed by `devtools::install_github` or related. Sorry. The simplest way to get the software is to download or clone the entire repo.

The code mostly uses base R capabilities but has a few dependencies: `RColorBrewer`, `akima`, `knitr`, and `pryr`. Since it's not a package, you have to manually install these packages if you don't already have them.

The recommended way to run the program is to `source` the file `R/run.R` into your R session; `R/run.R will source the rest. Once loaded, you can run the program by executing the statement`run()\` as shown below.

``` r
## This code block assumes your working directory is the root of the distribution.

source('R/run.R');
run();
```

This runs the program in a demo-like mode that quickly generates the data and produces the figures that appear in this README document. The default computation simulates 2,500 replications and produces 22 figures and one table. The data generation takes about 30 seconds on my small Linux server; the figures take about a minute, much of which is spent rendering the plots over a remote X11 connection.

The code that creates the outout (figures and table) for this README document is in `R/doc_readme.R`.

To rerun the program from scratch you must specify `clean=T`, since by default the program reuses existing data.

``` r
## This code block assumes your working directory is the root of the distribution
## and you've already sourced R/readme.R into your session

run(clean=T);                # delete data from previous run and rerun program
```

You can run each part separately by running one of the statements below.

``` r
## This code block assumes your working directory is the root of the distribution
## and you've already sourced R/readme.R into your session

init(doc='readme',clean=T);  # you MUST specify doc to run the program in pieces
dosim();                     # run simulations
domeand();                   # calculate mean significant observed effect sizes
dodoc();                     # generate data and figures
dodoc(save.out=F);           # generate data and figures without saving them
dodoc(figscreen=F);          # generate and save figures without plotting to screen. faster!
```

The program can also generate the data and figures for the blog post associated with the project. To generate these, execute `run()` with a suitable `doc` argument as shown below.

``` r
## This code block assumes your working directory is the root of the distribution.

source('R/run.R');
run(doc='ovrfx');            # generate data and figures for blog post 
```

Figures
-------

The default mode produces figures that illustrate the kinds of graphs the program can produce.

1.  *d* vs. *d* scatter plot colored by p-value; by default *d*<sub>*p**o**p*</sub> vs. *d*<sub>*s**d**z*</sub>
2.  histogram typically of *d*<sub>*s**d**z*</sub> colored by p-value
3.  probability distribution vs. *d*<sub>*s**d**z*</sub> colored by p-value; can compute probability density or cumulative probability internally, or you can pass the distribution into the function
4.  multiple line plot; my adaptation of R's `matplot`

The first block of figures show *d* vs. *d* scatter plots for *n* = 20 and *n* = 100. You may notice that the p-value color legends are of different sizes. This is a limitation of the code that draws these legends.

<img src="figure/readme/figure_01-001a_plotdvsd_scat_pop_sdz.png" width="50%" /><img src="figure/readme/figure_01-001b_plotdvsd_scat_sdz_pop.png" width="50%" /><img src="figure/readme/figure_01-002a_plotdvsd_scat_pop_sdz.png" width="50%" /><img src="figure/readme/figure_01-002b_plotdvsd_scat_sdz_pop.png" width="50%" />

The next block are similar but zoom in to the critical region where p-values switch from nonsignificant to significant.

<img src="figure/readme/figure_01-003a_plotdvsd_zoom.png" width="50%" /><img src="figure/readme/figure_01-003b_plotdvsd_zoom.png" width="50%" />

Next are four histograms for *d* = 0.2 and *d* = 0.8 by *n* = 20 and *n* = 100. You may notice that the p-value color legends are in different places and have different sizes. Again, these are limitations of the code.

<img src="figure/readme/figure_02-001a_plothist_hist.png" width="25%" /><img src="figure/readme/figure_02-001b_plothist_hist.png" width="25%" /><img src="figure/readme/figure_02-002a_plothist_hist.png" width="25%" /><img src="figure/readme/figure_02-002b_plothist_hist.png" width="25%" />

The next block are NULL and sampling distributions for *d* = 0.2, *n* = 20, and *d* = 0.8, *n* = 100.

<img src="figure/readme/figure_03-001a_plotpvsd_distn.png" width="25%" /><img src="figure/readme/figure_03-001b_plotpvsd_distn.png" width="25%" /><img src="figure/readme/figure_03-001c_plotpvsd_distn.png" width="25%" /><img src="figure/readme/figure_03-001d_plotpvsd_distn.png" width="25%" /><img src="figure/readme/figure_03-002a_plotpvsd_distn.png" width="25%" /><img src="figure/readme/figure_03-002b_plotpvsd_distn.png" width="25%" /><img src="figure/readme/figure_03-002c_plotpvsd_distn.png" width="25%" /><img src="figure/readme/figure_03-002d_plotpvsd_distn.png" width="25%" />

The final block are line plots of mean significant observed effect size (*m**e**a**n**d*<sub>*e**m**p**i*</sub>, *m**e**a**n**d*<sub>*t**h**e**o*</sub>) with various smoothing functions applied/

<img src="figure/readme/figure_04-001_plotm_meand.png" width="25%" /><img src="figure/readme/figure_04-002_plotm_meand.png" width="25%" /><img src="figure/readme/figure_04-003_plotm_meand.png" width="25%" /><img src="figure/readme/figure_04-004_plotm_meand.png" width="25%" />

Finally, here is a table of supporting data.

|  d.crit|  md20\_0.2|  md20\_0.5|  md20\_0.8|  ov20\_0.2|  ov20\_0.5|  ov20\_0.8|  nov\_0.2|  nov\_0.5|  nov\_0.8|
|-------:|----------:|----------:|----------:|----------:|----------:|----------:|---------:|---------:|---------:|
|    0.64|       0.79|       0.86|       0.98|       3.96|       1.72|       1.22|    166.52|      46.6|     16.65|

Comments Please!
----------------

Please post comments on [Twitter](https://twitter.com/gnatgoodman) or [Facebook](https://www.facebook.com/nathan.goodman.3367), or contact me by email <natg@shore.net>.

Please report bugs, other software problems, and feature requests using the [GitHub Issue Tracker](https://github.com/natgoodman/effit/issues). I will be notified, and you'll be apprised of progress. As already noted, the software is still rough and software documentation nonexistent.

Copyright & License
-------------------

Copyright (c) 2019 Nathan Goodman

The software is **open source and free**, released under the [MIT License](https://opensource.org/licenses/MIT). The documentation is **open access**, released under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0).
