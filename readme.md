# XFDTD

The project uses c++ to implement a time domain finite difference (FDTD) algorithm.

## Getting Stared

You will require the the following libraries on your system to compile XFDTD: [BLAS](https://www.netlib.org/blas/), [LAPACK](https://www.netlib.org/lapack/), [FFTW](http://www.fftw.org/),[xtensor](https://github.com/xtensor-stack/xtensor),[xtensor-blas](https://github.com/xtensor-stack/xtensor-blas), [xtensor-fftw](https://github.com/xtensor-stack/xtensor-fftw).

## Examples

You can see some examples in `examples` directory.

Test the example:

```bash
$ chmod +x ./run_all_examples.sh
$ cmake -B build && cmake -DCMAKE_BUILD_TYPE=Release build &&cmake --build build  -j $(nproc) && ./run_all_examples.sh
```

There will be data in `visualizing_data/data` directory if everything goes well.

## TODO
- [ ] Add the plane wave source condition
- [ ] Add concept `Flux`
- [ ] A CLI to complete FDTD simulate with config text
- [ ] Apply to periodic structures
- [ ] Parallel computing
- [ ] Add documentation

## Problems

1. There is a problem with the calculation of near field to far field transformation. The result is not enough accurate.
For example. vCompared with Mie theory, perfect conductor sphere scattering with bistatic radar cross section (RCS) has an average error of 1dB.
2. Disable to calculate dispersion material for failing to compline this lib on Ubuntu 20.04.1 LTS. (If you want to use it, go `src/object/object.cpp` , uncomment the line 111 and 115 and comment the line 110. It will work on Mac OS.)
