<!-- badges: start -->  [![R build status](https://github.com/alextkalinka/linkcomm/workflows/R-CMD-check/badge.svg)](https://github.com/alextkalinka/linkcomm/actions)  <!-- badges: end -->


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
nodes_c1 <- getNodesIn(lm, clusterids = 1)

# Nodes shared by communities 10 and 11:
nodes_sh <- get.shared.nodes(lm, comms = 10:11)

# Community connectedness score:
comm.conn <- getCommunityConnectedness(lm)

```
