# The Macrocirculation submodule

## Directory structure

- apps: contains the demo applications copied to bin/macrocirculation
- data: contains our test geometries  
- scripts: contains utility scripts for creating e.g. the geometry 
- src: the source code of the library behind this submodule 
- unittests: A few catch2 unittests for critical components

## Run unittests

The unittests can be run right after building the project with ctest
```
cmake <cmake-options-here> ..
make 
ctest
```