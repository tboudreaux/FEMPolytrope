# Finite Element Modeling Of a Polytrope
This is a very basic project where I am attempting to learn how to use finite element modeling and MFEM.

Note that this code *does not currently solve the lane-emden equation
properly*. Instead it requires that θ(0) = 1 and θ(ξₛ) = 0 along with ξₛ being
manually specificed (this differes from the general formulation where the
Neuman condition θ'(0) = 0 is used as opposed to a essential BC at the
surface). The issue here is that I am still learning how to work with Neuman
BCs and essential BCs at the same DOF in FEM. However, if for example you set n
= 1.5 and ξₛ = 3.6538 you find the correct profile for θ.

## Building
In order to build this you will need meson and MFEM installed. Because this
code is not intended for distribution there are some hardcoded paths you will
need to update. These should all be in the root meson.build file.

once you have those updated it should be as simple as

```bash
./mk
```

You can then run the test by running

```bash
./build/src/solver/laneEmden
```

In order to see command line options you can use `--help`. the option `--root` is used to specify ξₛ.

or if you want to build and run back to back

```bash
./mk --run
```

If you want to build with debugging symbols

```bash
./mk --debug
```
