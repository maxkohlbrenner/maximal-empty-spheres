# CGAL implementation
`empty_spheres_reconstruction_3.cpp` is a very basic example that constructs the contact points for a set of SDF spheres using the experimental Maximal_empty_spheres CGAL module.

input/output is currently only supported via files.

# Building:
Make sure the cgal/ submodule in this folder is initialized. The Maximal_empty_spheres CGAL module is not available in the main branch as of now.

```
mkdir build
cmake -B build -DCGAL_DIR=./cgal
cmake --build build/ --parallel
```

the build/ folder should now contain the executable.
