# Lifting Directional Fields to Minimal Sections

This code implements the minimal section lifting method detailed in the paper

> David Palmer, Albert Chern, and Justin Solomon, __"Lifting Directional Fields to Minimal Sections."__ _ACM Transactions on Graphics_ 43, 4, Article 60 (July 2024). https://doi.org/10.1145/3658198.

The main field optimization code is in `minsec`. First, run
```
meshData = process_surface(verts, faces, Flat=flat, NormalizeArea=area);
```
to generate a surface data structure encoding various basic information about the surface. `flat = true` indicates that the surface is a flat (planar) domain so that the local frames will be chosen to be identical and parallel transport operators will be identity operators—useful for visualization of generated sections. The parameters of the minimal section optimization are defined relative to the area of the underlying surface, set by the `area` parameter. To compute a minimal section, run
```
[S, C, f, z, admmData] = minsec(meshData, N, lambda, degree, fiberRadius, options, vizOpts);
```
where:
- `N` is the number of discrete subdivisions of the fiber. We find `N = 64` to be a suitable choice for most surfaces.
- `lambda` is the regularization weight on the mass norm of the singular measure, denoted $\lambda$ in the paper. Increasing it will penalize proliferation of singularities.
- `degree` is the degree of the bundle. E.g., set `degree = 1` for unit vector fields, `degree = 2` for line fields, or `degree = 4` for cross fields.
- `fiberRadius` is the radius of the circular fibers, setting the length scale of the minimal section. This is the parameter $r$ in the paper. At higher values of $r$, the mass norm behaves more like total variation energy, while reducing $r$ makes it behave more like Dirichlet energy.
- `options` includes options for fixing or masking singularities, the maximum number of iterations, and a flag for continuous vs. finite-difference vertical derivatives.
- `vizOpts` includes options for online visualization of the current as it converges, as well as for saving the results to a file.

## Visualization

The `plot` directory includes various functions for visualizing currents, sections, and fields.

- `plot_section_volume` displays a (pseudo-)volume-render of a current (only applicable over flat domains).
- `plot_section` displays polar plots of the magnitude of the current over sampled points on the base surface.
- `plot_fourier_components` plots the Fourier components of the scalar function `f` parametrizing the exact part of a current. The phases are only meaningful on a flat domain.
- `plot_section_surface` displays the graph of an extracted section `z` (only applicable over flat domains).
- `plot_intrinsic_field` plots an extracted section in the manner of a quiver plot.
- `plot_integral_curves` generates and plots the integral curves of an
extracted section `z`.
- `plot_singularities` plots the singularities of an extracted section `z` as spheres.

The `figures` directory includes scripts to generate figures in the paper and render them to files. `faces_experiment` must be invoked before `faces_experiment_figures` or `faces_mesh_figures`. For `wavy_disk_experiment`, the values of `eccentricity` used in the paper are `1` and `0.7`. The `dataset_experiment` script was used to process the dataset from Myles et al. (2014), "Robust Field-aligned Global Parametrization." For the figures using faces from the Basel faces dataset, you will need to download the dataset `model2019_face12.h5` to the `figures/basel_faces` directory.

## Dependencies

Some of the code may require functions from `gptoolbox`. `minsec_cvx` requires the `CVX` convex-optimization toolkit and an interior-point optimizer such as `MOSEK`. The `ext` directory contains adapters to external C++ code from the `Directional` and `Geometry Central` projects. Compilation instructions are provided in the `README.md` in that directory.