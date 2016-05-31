# ClusterGeneratorBox
Generates a fiber cluster (box) given its size and fiber properties.

Input
----------------

The program requires a project organized as:

- Project
  - INPUT
    - input.in
    - initial_histogram.in (optional)
  - OUTPUT
    - coords.txt


The file input.in should have the following format:
```
&input
fiber_diameter      =           !Average fiber diameter
fiber_length        =           !Average fiber length
nbr_segments        =           !number of segments per fiber
volumetric_fraction =           !Volumetric Fraction
max_theta           =           !Maximum fiber inclination wtih respect of Z
box_size            =           !Cluster Size
box_location        =           !Cluster center location
/
```

The file initial_histogram.in contains information fo the 2D orientation of the fibers in the XY plane with the 
orientation distribution function. 

The file initial_histogram.in should have the following format:
```
#bins
Distribution function value for bin 1
Distribution function value for bin 2
.
.
.
Distribution function value for bin n
```
