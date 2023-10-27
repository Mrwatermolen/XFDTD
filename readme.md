# XFDTD

The project uses c++ to implement a time domain finite difference (FDTD) algorithm.

## Getting Stared

You will require the the following libraries on your system to compile XFDTD: [BLAS](https://www.netlib.org/blas/), [LAPACK](https://www.netlib.org/lapack/), [FFTW](http://www.fftw.org/),[xtensor](https://github.com/xtensor-stack/xtensor),[xtensor-blas](https://github.com/xtensor-stack/xtensor-blas), [xtensor-fftw](https://github.com/xtensor-stack/xtensor-fftw).

You can see some examples in the `test` folder.

## TODO
- [ ] Add the plane wave source condition
- [ ] Add concept `Flux`
- [ ] A CLI to complete FDTD simulate with config text
- [ ] Apply to periodic structures
- [ ] Parallel computing
