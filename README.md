---
title: "The effect of influenza receptor binding on antigenic drift"
author: "James Hay"
output: 
  html_document
header-includes:
  - \usepackage{amsmath}
  - \usepackage{amssymb}
  - \usepackage{amsmath}
  - \usepackage{amsfonts}
  - \usepackage{amsthm}
  - \usepackage{amsmath}
  - \usepackage{mathtools}
  - \newcommand{\HRule}{\rule{\linewidth}{0.4mm}}
---


[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

This git repository contains all of the code needed to run the binding avidity adaptation project either through command line or the accompanying shiny app.

## Installation
Installation is fairly straightforward. You should 

```r
devtools::install_github("jameshay218/driftSim")
```

## Usage

```r
library(driftSim)
```

## License

GPL-2 © [James Hay &lt;james.hay13@imperial.ac.uk&gt;](https://github.com/).



## Equations
### 1. Probability of evading immune system
After entering a host, the virus must first evade the immune response (Equation 1). Here, the probability of escaping the immune response increases as binding avidity increases. The virus must evade the immunity conferred from all previous infections (k), adjusted by the antigenic distance between the virus and the virus that elicited the host's immunity. This relationship can be seen in plot C. The rate of change of $f$ with respect to binding avidity is shown in plot E.

$$f(k,V_i) = [1-e^{-p(V_i+q)}]^{rk - \delta_{ji}}$$

### 2. Probability of Successful Replication Within Host
Binding avidity also affects how well a virus is able to replicate within the host, as described by Equation 2. This relationship is shown in plot A. The naive case (ie. $k=0$) is shown in plot D.

$$g(V_i) = e^{-aV_i^b}$$

### 3. Within Host Reproductive Number
The within-host reproductive number is therefore given by the product of the probability of successful within host infection, and the number of offspring virions produced per event:

$$R_{in} = n \cdot \phi(H_k,V_i)$$

### 4. Infectiousness
The infectiousness of a particular virus is therefore related to the within-host reproductive number and the number of initially infecting virions as follows:

$$ \rho = 1 - (\frac{1}{R_{in}})^{-v} = 1-(\frac{1}{\phi(H_k, V_i)})^{-nv} $$

### 5. Transmission Rate
The transmission rate between hosts, $\beta$, is therefore given by the product of the infectiousness of that virus and the contact rate between hosts. This relationship is shown in plot B. The rate of change of $\beta$ with respect to binding avidity is shown in plot F.

$$\beta = c \cdot \rho$$

