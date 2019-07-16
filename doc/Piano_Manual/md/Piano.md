# Piano

Piano is a molecular dynamics program developed by Haoyu Lin, the first two characters of which stands for "path integral".

* Author: Haoyu Lin
* E-mail: vileoy@pku.edu.cn
* Repo: [uovie/Piano](https://github.com/uovie/Piano) `private`

## 1 Releases

The version 1.0.0 beta of Piano (`Piano-1.0.0`) has completed, but not been released. Currently, it only integrates tiny part of functions, including certain classical dynamic simulations in the constant temperature and quantum statistical calculations based on path integral molecular dynamics.

## 2 Instructions

* usage tips

  Piano only accepts `.vie` type input files. Type the following command to run a calculation:

  ```shell
  ./Piano-1.0.0 test.vie
  ```

* compile tips

  If you want to compile it by yourself, just invoke `make` under `src` directory .

## 3 Supports

* Classical Dynamics
  * Langevin dynamics (side; middle)
  * Andersen thermostat (side; middle)
  * Nose-Hoover chain (global, side; global, local; local, side; local, middle)
* Quantum Statistical Mechanics
  * PIMD via LD (side; middle)
  * PIMD via AT (side; middle)
  * PIMD via NHC (local, side; local, middle)