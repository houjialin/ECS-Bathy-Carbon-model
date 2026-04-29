# Subseafloor Microbial and organic carbon dynamic model
## Overview
This repository contains the R-based bioenergetic model developed to mechanistically resolve the survival strategies behind subsurface microbial succession and quantitatively evaluate their role in sedimentary organic carbon (OC) degradation over geological timescales.
## Model Features
- 2-G Carbon Mineralization: Partitions Total Organic Carbon (TOC) into a reactive labile pool (OC1) and a complex recalcitrant pool (OC2).
- Ecophysiological Guilds: Simulates community dynamics through four distinct microbial groups, parameterized by their specific life-history strategies and substrate affinities:
  - B1 (Copiotrophic Bacteria): Prioritizes rapid biomass turnover utilizing OC1, characterized by high maximum growth rates, yields, and maintenance costs.
  - B2 (Generalist Bacteria): Targets recalcitrant OC2 with comparatively reduced growth rates and maintenance demands.
  - A1 (Moderate Archaea): Utilizes labile OC1 via a balanced, moderately growth-maintenance strategy.
  - A2 (Oligotrophic Archaea，specifically refers to Bathyarchaeia in this study): Deep-biosphere adapted extremophiles with superior affinity for OC2, characterized by ultra-low mortality and minimal basal energy demands.
- Bioenergetic ODE System: Solves coupled Ordinary Differential Equations (ODEs) to track carbon transfers driven by substrate uptake, growth, state-dependent maintenance (transitioning between exogenous and endogenous power), mortality, and necromass recycling.
## Implementation
- Language/Environment: Built using the R programming language.
- Numerical Solver: Integrates equations using the lsoda solver from the deSolve package.
- Timescale & Constraints: Simulates a 1,000-year burial timescale, with initial conditions and physiological parameters tightly constrained by empirical observations and recent bioenergetic estimations.

## Citation 
Physiological trade-offs drive the archaeal dominance and carbon turnover in deep subsurface, Jialin Hou & Lewen Liang et al., 2026 (submitting to BioRxiv)
