This program is HW3D.
There are 4 targets in this program:
HW3D - main target, for user's calls;
Tests - target, for tests;

HW3D:
The first number is number of triangles.
Next numbers is coordinates of points of triangles (3 numbers per point, 3 points per triangle).

Tests:
In tests.txt:
The first number is the number of tests.
The next numbers are the expected value (true or false) and the coordinates of each triangle.

It is expected that the files will not be moved to the project repository.

As you can see, only the krvlib::triangle::intersect() method is tested, because the performance of the program depends on it.
