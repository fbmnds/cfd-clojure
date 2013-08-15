## Computational Fluid Dynamics

Prof. Lorena Barba, Dr. Rio Yokota, [CFD Python: 12 steps to Navier-Stokes](http://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/).


### Scope

- Step 1: Linear Convection
- Step 2: Non-Linear Convection
- Step 3: 1D Diffusion
- Step 4: Burgers' Equation
- Step 5: 2D Linear Convection
- Step 6: Convection 2D
- Step 7: Diffusion 2D

### How to run the models

The project uses [Midje](https://github.com/marick/Midje/) to run the models.

`lein midje` will run all models and generate graphics/animation.

This assumes a valid Gnu environment in place. For Windows, [GoW](https://github.com/bmatzelle/gow/wiki) is recommended.

`lein midje :filter -mpeg` will run all models without generating graphics/animation. Models may be selected by filtering [(Midje filters)](https://github.com/marick/Midje/wiki/Using-metadata-to-filter-facts#lein-midje-command-line-arguments) using `:step1 ...`.

### License

CFD Clojure © 2013 Friedrich Boeckh, distributed under the Eclipse Public License, the same as Clojure.
CFD Python is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
