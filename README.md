# Cartel-detection in Rational Complex Networks
This repository contains the code used to produce the results and plots for my [C5.4 Networks](https://courses-archive.maths.ox.ac.uk/node/49460) final assessment and the paper that came of it. The code in the state it was submitted, and thankfully not marked, was [commit d1f4b0d](https://github.com/jacobusmmsmit/network_cartels/commit/d1f4b0d509d8fa8b8dee779650c30521ff057c6a).

In its current state, a lot of the code went unused as I refactored it just before submitting, but also left the old code up because of the time limit.

To run the code, please install [Julia 1.7.X or above](https://julialang.org/downloads/), open the repository in your text editor of choice, and run the following commands to install the required packages:
```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Now the code from the `src/` folder can be run. Relevant files to reproduce results in paper:
* `scripts/plot_avocado.jl` - Creates Fig. 1
* `scripts/paper_results.jl` - Creates Fig. 2, Table. 1

In the end, I only got an okay mark so maybe it was a bad idea after all.
