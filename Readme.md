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

## Weak form
We start with the strong form of the lane-emden equation

$$
\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^{2}\frac{d\theta}{d\xi}\right) + \theta^{n} = 0
$$

where $n$ is the polytropic index. Note how when $\xi=0$ there is a singularity. We multiply both sides by $\xi^{2}$ to avoid this singularity

$$
\frac{d}{d\xi}\left(\xi^{2}\frac{d\theta}{d\xi}\right) + \xi^{2}\theta^{n} = 0
$$

now to put this in its weak form we multiple by a test function ($v(\xi)$) and integrate over some domain $\Omega$

$$
\int_{\Omega} v(\xi)\frac{d}{d\xi}\left(\xi^{2}\frac{d\theta}{d\xi}\right)d\Omega + \int_{\Omega}v(\xi)\xi^{2}\theta^{n}d\Omega = 0
$$

By integration by parts the first term can be rewritten as

$$
\int_{\Omega} v(\xi)\frac{d}{d\xi}\left(\xi^{2}\frac{d\theta}{d\xi}\right)d\Omega = \left[v(\xi)\xi^{2}\frac{d\theta}{d\xi}\right]_ {\partial \Omega} - \int_{\Omega} r^{2}\frac{dv}{d\xi}\frac{d\theta}{d\xi}
$$

so then the weak form of the lane-emden equation is

$$
\left[v(\xi)\xi^{2}\frac{d\theta}{d\xi}\right]_ {\partial \Omega} - \int_{\Omega} r^{2}\frac{dv}{d\xi}\frac{d\theta}{d\xi} + \int_{\Omega}v(\xi)\xi^{2}\theta^{n}d\Omega = 0
$$

$$
\left[v(\xi)\xi^{2}\frac{d\theta}{d\xi}\right]_ {\partial \Omega} - \int_{\Omega} r^{2}\frac{dv}{d\xi}\frac{d\theta}{d\xi} + \int_{\Omega}v(\xi)\xi^{2}\theta^{n}d\Omega = 0
$$

we know that we want to let $\frac{d\theta}{d\xi}(\xi=0) = 0$ so 

$$
v(\xi)\xi^{2}\frac{d\theta}{d\xi}|_ {\xi=\xi_{s}} - \int_{\Omega} r^{2}\frac{dv}{d\xi}\frac{d\theta}{d\xi} + \int_{\Omega}v(\xi)\xi^{2}\theta^{n}d\Omega = 0
$$

this is the final weak form of our expression

## Representing this in MFEM

We can represent this weak form as a nonlinear form composed of a:

1. Domain `DiffusionIntegrator` with a `VectorFunctionCoefficient` evaluated as $-\xi^{2}$
2. Domain `NonlinearPowerIntegrator` with a `FunctionCoefficient evaulated as $\xi^{2}$. Note that this needs to be implimented as a custom integrator. This is done in `src/utils/mfemUtils.cpp`.
3. Some kind of domain integrator I think (the first term is the one I am still not sure how to deal with.
