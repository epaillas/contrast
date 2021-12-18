### CONTRAST (COrrelatioN funcTions foR ASTrophysics)

A  package for the calculation of galaxy clustering statistics, including correlation functions and velocity statistics.

The interface of the code is written in Python 3, while the engine that takes care of the intensive calculations is written in the Julia programming language, so you must have it installed along with the following Julia packages:

### Julia requirements

    - `CellListMap`
    - `StaticArrays`
    - `LinearAlgebra`

You will also require the following Python packages:

### Python requirements

  - `pyjulia` >= 0.5.7
  - `numpy` >= 1.20.1

### Installation

Under the main directory, install the package with `python setup.py install`. 

### Examples

Example notebooks can be found under the `examples/` directory.
