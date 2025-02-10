# Finite Element Modeling Of a Polytrope
This is a very basic project where I am attempting to learn how to use finite element modeling and MFEM.

Note that this code *does not currently solve the lane-emden equation properly*.

## Building
In order to build this you will need meson and MFEM installed. Because this code is 
not intended for distribution there are some hardcoded paths you will need to update. These should all be in the root meson.build file.

once you have those updated it should be as simple as

```bash
./mk
```

You can then run the test by running

```bash
./build/src/solver/laneEmden
```

or if you want to build and run back to back

```bash
./mk --run
```

If you want to build with debugging symbols

```bash
./mk --debug
```
