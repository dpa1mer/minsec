# How to Use

This folder contains various adapters to external code written in C++. To use
them, the MEX `.cpp` files must be compiled with the `mexbuild` script,
after replacing the placeholder paths with their concrete values on your
system.

After compiling the MEX code, do not run the MEX scripts directly. Instead,
invoke the provided helper MATLAB scripts. For example, to compute a field by
the method of Knöppel et al. (2013) with the `Geometry Central` library
and then use it to make a quad mesh with the `Directional` library:
```
z = compute_field_geometry_central(meshData, 4);
[V, D, F] = quad_mesh_directional(meshData, z, options);
```
where `meshData` is the output of `process_surface` and `z` is an optimized
cross field encoded as a complex value in the fourth tensor power of the
tangent bundle (represented here in per-vertex local frames).