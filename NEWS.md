Revision history for misig repository
================
Nathan (Nat) Goodman
February 12, 2019

<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->
Release 0.90 2019-02-13
-----------------------

Support second version of ovrfx blog post. Current title "When You Select Significant Findings, You’re Selecting Inflated Estimates". Hopefully near final version.

Major code changes

-   Abandon file caching. Way more trouble than it's worth
-   Handle multiple types of simulation
    -   `rand` - the original siumulation scheme: randomly choose 'd' and run one simulation per choice
    -   `fixd` - more typical simulation scheme: set 'd' to certain fixed values and run many simulations for each value
-   Handle two ways of calculating mean significant effect sizes, 'meand'
    -   `empi` - estimate empirically from simulation data
    -   `theo` - theoretical calculation from sampling distribution

Changed files

-   `R/datman.R`
    -   Add functions to `save` and `load` new file types: `sim_rand`, `sim_fixd`, `meand_empi`, `meand_theo`
    -   Replace generic functions, eg, `save_nd`, by ones specific to each file type. Though seemingly less factored, it's way simpler now that I've abandoned caching, esp. since we have so few file types at present
    -   Remove `keep` functions. These supported now-abandoned file caching
    -   Remove caching related code from `save` and `load` functions
    -   Redefine `get` functions to be synonyms for `load`. Previous distinction was needed for now-abandoned file caching
-   `R/doc_ovrfx.R`
    -   Near complete rewrite for second version of blog post
-   `R/init.R`
    -   Add params for multiple types of simulation: `sim_rand`, `sim_fixd`
    -   Add params for calculation of mean significant effect sizes: `meand`
    -   Remove params for caching: `load`, `keep`, `clean.cache`
    -   Remove code for initializing cache
    -   `datadir` no longer has `m` component. Makes no sense with multiple simulation types
-   `R/ovrfx.R`
    -   Add function for `fixd` simulation and rename function for original `rand` simulation: `dosim_fixd`, `dosim_rand`
    -   Add function for empirical calculation of mean significant effect size and rename function for theoretical calculation: `domeand_empi`, `domeand_theo`
    -   Incorporate new `R/datman.R` functions for saving files
    -   Remove code for getting data from cache if possible
-   `R/plot.R`
    -   Add function to plot histogram: `plothist`
    -   Fiddle with 'extra' line params and code. Removed specialized params like `d.crit` and `d.pop` in favor of the more general `vline` and `hline`
    -   Add `x0` param to `pval_legend` function to adjust location. Crude but adequate for immediate needs
    -   Add convenience function `d2col` that maps `d` values to colors by combining `d2pval` and `pval2col`

Release 0.50 2019-02-03
-----------------------

Support first real version of ovrfx blog post ready for external review. Current title "Significance Testing Makes Bad Choices"

-   Rename repo from 'effit' to 'misig'. Better reflects the content I expect for this repo
-   Rename first blog post from 'siglo' to 'ovrfx' - stands for 'overestimate with fixed effect size'

Release 0.10 2019-01-01
-----------------------

First version. Contains boilerplate files such as this one and the document `sigbi.Rmd` (presently just a stub). Will contain code to support the document when released.

Copyright & License
-------------------

Copyright (c) 2019 Nathan Goodman

The software is **open source and free**, released under the [MIT License](https://opensource.org/licenses/MIT). The documentation is **open access**, released under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0).
