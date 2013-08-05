# Computational Fluid Dynamics

Prof. Lorena Barba, Dr. Rio Yokota, [CFD Python: 12 steps to Navier-Stokes](http://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/).


## How to run the tests

The project uses [Midje](https://github.com/marick/Midje/).

`lein midje` will run all tests and generate graphics/animation.

This assumes a valid Gnu environment in place. For Windows, [GoW](https://github.com/bmatzelle/gow/wiki) is recommended.

`lein midje :filter -mpeg` will run all tests without generation graphics/animation. The tests support further filtering ([Midje filters])(https://github.com/marick/Midje/wiki/Using-metadata-to-filter-facts#lein-midje-command-line-arguments) by using `:step1 ...`.

## License

CFD Python is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
CFD Clojure © 2013 Friedrich Boeckh, distributed under the Eclipse Public License, the same as Clojure.
