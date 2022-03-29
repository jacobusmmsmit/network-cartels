# Cartel-detection in Rational Complex Networks
This repository contains the code used in my [C5.4 Networks](https://courses-archive.maths.ox.ac.uk/node/49460) final assessment, the paper of which can be found [here](https://www.jacobussmit.com) (when released). The code in the state it was submitted, and thankfully not marked, was [commit d1f4b0d](https://github.com/jacobusmmsmit/network_cartels/commit/d1f4b0d509d8fa8b8dee779650c30521ff057c6a).

To run the code, please install [Julia 1.7.X or above](https://julialang.org/downloads/), open the repository in your text editor of choice, and run the following commands to install the required packages:
```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Now the code from the `src/` folder can be run. Relevant files to reproduce results in paper:
* `scripts/plot_avocado.jl` - Creates Fig. 1
* `scripts/paper_results.jl` - Creates Fig. 2, Table. 1
