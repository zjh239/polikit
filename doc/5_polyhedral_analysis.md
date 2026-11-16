## 5. Polyhedral analysis

The connectivity of tetrahedral is computed as:

1. Loop over the full list of polyhedra, if it is a cation, it means a center of polyhedron. Check if the coordination number is 4, which means a tetrahedron. Note the difference between polyhedral neighbor number and center atom coordination number. The former is a coarse-grained coordination, the latter in this case defines the shape of the polyhedron.

2. For a tetrahedron, loop over its corner oxygen atoms to check if it is connected with other tetrahedra by sharing itself. If so, mark this corner oxygen, and +1 to the connectivity counter. In other words, we add 1 to the center cation's T_n value when its neighboring polyhedron meets two conditions:

- It is a tetrahedron, the center atom has coordination number of 4;
- It is corner connected with the tetrahedron being investigated;


