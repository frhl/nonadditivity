# Power Analysis Plots for Reviewer

## Quick Start

Open the script in RStudio and step through it, or run from command line:

```bash
Rscript create_reviewer_plots.R
```

## What it does

Creates 4 QQ plots for reviewer response:
- **Figure R7**: Additive architecture (shows non-additive model is null)
- **Figure R8**: Recessive architecture (shows non-additive model captures signal)
- **Figure R9**: Partially recessive (shows both models detect signal)
- **Figure R10**: Overdominant (shows non-additive is primary driver)

## Modify parameters

At the top of the script, change:
```r
H2 <- 0.01      # heritability (0.01 or 0.20)
N_CAUSAL <- 10  # number of causal variants (10 or 20)
```

## Output

PDFs saved to: `/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/power_simulation_plots/`

## Structure

The script has 7 clear steps:
1. Configuration
2. Load all data
3. Helper function for QQ plots
4. Create Figure R7 (additive)
5. Create Figure R8 (recessive)
6. Create Figure R9 (partially recessive)
7. Create Figure R10 (overdominant)

Easy to debug - just run each section one at a time.
