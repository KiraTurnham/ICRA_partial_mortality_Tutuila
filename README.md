# From stability to decline: partial mortality of *Isopora crateriformis* in Tutuila

## Overview

This repository contains the data, R scripts, and selected outputs supporting:

> Turnham K, Huntington B, Couch CS, Sena E, Tanaka K, Oliver T. From stability to decline: threefold increase in partial mortality of threatened *Isopora crateriformis* in American Samoa following the 2024 Global Coral Bleaching Event. *Coral Reefs*. 2026 Jul 9. [https://doi.org/10.1007/s00338-026-02912-0](https://doi.org/10.1007/s00338-026-02912-0)

The study analyzes *I. crateriformis* colony density and partial mortality around Tutuila, American Samoa, using surveys from 2015–2025. It compares surveys before the 2024 Global Coral Bleaching Event (GCBE) with 2025 observations and evaluates whether heat stress and colony size are associated with partial mortality.

The study found that mean partial mortality was 30% in 2025, compared with 8% in previous survey years. Changes were concentrated in medium and large colonies, while small colonies remained relatively unchanged. Partial mortality in medium and large colonies also increased with mean sea surface temperature, consistent with an impact of the 2024 GCBE. These findings support partial mortality as a way to assess bleaching impacts when surveys cannot be conducted during a heat-stress event, and indicate that spatial differences in heat stress may contribute to local refuges.

## Prerequisites

- R and, optionally, RStudio.
- The R packages used by the analysis scripts: `betareg`, `boot`, `broom`, `broom.mixed`, `corrplot`, `cowplot`, `DHARMa`, `dplyr`, `eds`, `emmeans`, `GGally`, `ggeffects`, `ggplot2`, `ggpubr`, `ggridges`, `ggsignif`, `ggtext`, `glmmTMB`, `hrbrthemes`, `lme4`, `lmtest`, `multcomp`, `patchwork`, `performance`, `purrr`, `readxl`, `rstatix`, `sf`, `statmod`, `svglite`, `tidyr`, `tidyverse`, and `viridis`.

Install packages from CRAN where available, using `install.packages("package-name")`. Some packages may require additional system libraries; see the package installation instructions for your operating system.

## Installing and running

Clone or download this repository, then open `ICRA_partial_mortality_Tutuila.Rproj` in RStudio. The scripts use paths relative to the repository root, so run them with the project root as the working directory.

Run the numbered scripts in order:

1. `scripts/00_Prep_colony_level_data.R` prepares colony-level comparisons for NOAA NCRMP and 2025 survey data, applies size filters, and saves intermediate data.
2. `scripts/01_Clean_density_data_site_level.R` calculates site-level colony density and adjusted southern Tutuila density.
3. `scripts/02_colony_level_PM_size_summary.R` summarizes and visualizes colony partial mortality by year and size class.
4. `scripts/03_density_summary.R` summarizes and visualizes site-level density.
5. `scripts/04_site_map.R` creates maps of survey sites and density.
6. `scripts/05_glmm.R` analyzes the relationship between partial mortality and 1year prior satellite-derived SST.
7. `scripts/06_glmm_6m_supplement.R` runs the supplementary six-month prior satellite-derived SST analysis.

For example, from the repository root in a terminal:

```sh
Rscript scripts/00_Prep_colony_level_data.R
```

Repeat with each script in the order above. The numbered scripts use saved intermediate datasets in `data/`; in particular, the heat-stress analyses use prepared environmental data supplied in that directory. Figures are written to `plots/`, and tabular analysis outputs are written to `outputs/` or `data/`.

## Repository contents

- `scripts/` — numbered R scripts for data preparation, summaries, statistical analyses, and mapping.
- `data/` — colony- and site-level survey data, NOAA NCRMP comparison data, prepared environmental/heat-stress data, intermediate R data files, and Tutuila shapefile components.
- `outputs/` — selected CSV tables produced by the analyses.
- `plots/` — selected figures produced by the analyses, in PNG and PDF formats.

## Survey and comparison data

The 2025 *I. crateriformis* survey covered 63 randomly selected sites at 6–18 m depth, with three 2 × 30 m belt transects per site. Divers counted colonies at least 5 cm in diameter and recorded colony diameter and partial mortality for up to the first 30 colonies encountered per transect. Partial mortality is the percentage of colony surface consisting of old dead skeleton. Site summaries pool observations across the three transects.

Historical comparison data are from NOAA's National Coral Reef Monitoring Program (NCRMP), accessed through the [National Centers for Environmental Information accession 0312067](https://accession.nodc.noaa.gov/0312067). NCRMP benthic surveys were conducted in 2015, 2018, and 2023 at stratified-random shallow, medium, and deep forereef sites. NCRMP and 2025 survey methods differ and were not designed for direct methodological equivalency; interpret temporal comparisons with this limitation in mind.

## Resources

- [Published paper](https://doi.org/10.1007/s00338-026-02912-0)
- [NOAA NCRMP data, NCEI accession 0312067](https://accession.nodc.noaa.gov/0312067)
- [Environmental data summary script and associated satellite-derived data access links are available via] (https://github.com/krtanaka/eds).

## Version control platform

- Git
- GitHub: [KiraTurnham/ICRA_partial_mortality_Tutuila](https://github.com/KiraTurnham/ICRA_partial_mortality_Tutuila)

## License

No license file is currently included in this repository. Reuse and redistribution terms have not been specified.

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration or the United States Department of Commerce. Repository content is provided “as is,” and users assume responsibility for its use. Any claims against the Department of Commerce or its bureaus arising from use of this project are governed by applicable federal law. References to commercial products, processes, or services do not imply endorsement, recommendation, or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal of a Department bureau, may not be used to imply such endorsement.
