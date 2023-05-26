# Calculate_Moire_Vectors
A GUI for determing the moire vectors produced between two sets of real-space lattice vectors.
Written by Evan Telord (ejt2133@columbia.edu)


1) Download and run the file "Calculate_moire_GUI_rt_vectors.py" or "Calculate_moire_GUI_xy_vectors.py".

Note: There are two files containing in this repository that perform the same function. "Calculate_moire_GUI_rt_vectors.py" accepts vector inputs in the form of vector length 1, vector length 2, and angle between the vectors. For example, a square unit cell input would be something like "1","1","90". "Calculate_moire_GUI_xy_vectors.py" accepts vector inputs in the form of vector 1 x component, vector 1 y component, vector 2 x component, vector 2 y component. For example, a square unit cell input would be something like "1","0","0","1".

2) A GUI will pop up requesting the input vector information. You will need to enter two sets of vectors. Each set corresponds to the real-space lattice vectors of the underlying crystals. You will also need to enter the desired twist angle between the sets of real-space lattice vectors.
3) Click "Calculate Moire Vectors" to calculate the moire vectors.
4) The moire vectors informtion will be printed in the GUI below the "Calculate Moire Vectors" button. Additionally, three plots will be generated consisting of quivers plots of left: the real-space lattice vectors twisted with respect to one another, middle: the reciprocal-space lattice vectors twisted with respect to one another along with the calculate moire reciprocal vectors, right: calculated real-space moire lattice vectors.

Note: the code works perfectly for angles below any possible symmetry operations. For angles above the symmetry operations, it should work for most use cases, but it hasn't been tested thoroughly. To deal with angles larger than the lowest symmetry operation, the code calculates the rotational symmetry of the reciprocal vectors, calculates the symmetry operation angle, and then performs the symmetry operation to shift the angle back to an appropriate value for the moire-vector calculation.


