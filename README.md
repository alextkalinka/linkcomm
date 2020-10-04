<!-- badges: start -->
  [![R build status](https://github.com/alextkalinka/linkcomm/workflows/R-CMD-check/badge.svg)](https://github.com/alextkalinka/linkcomm/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/linkcomm)](https://CRAN.R-project.org/package=linkcomm)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/linkcomm)](https://cran.r-project.org/package=linkcomm)
<!-- badges: end -->


# `linkcomm`

## Summary
`linkcomm` is an `R` package that provides tools for generating, visualising, and analysing [Link Communities](https://www.nature.com/articles/nature09182) in networks. See the companion [paper](https://academic.oup.com/bioinformatics/article/27/14/2011/194743) for more information.

## Installation

```r
install.packages("linkcomm")
```

## Usage

```r
# Explore the in-built Les Miserables network:
library(linkcomm)
lm <- getLinkCommunities(lesmiserables)
```

![](./imgs/summ.png)

```r
# Visualize the communities:
plot(lm, type = "graph", layout="spencer.circle")
```

![](./imgs/spenc-circ.png)

```r
# Extract the nodes from the first community:
n1 <- getNodesIn(lm, clusterids = 1)

# Nodes shared by communities 10 and 11:
ns <- get.shared.nodes(lm, comms = 10:11)

# Node community centrality scores:
nc <- getCommunityCentrality(lm)

# Community connectedness scores:
cc <- getCommunityConnectedness(lm)

```