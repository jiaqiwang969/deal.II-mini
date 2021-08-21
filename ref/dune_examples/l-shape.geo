// Scale global grid element size
Mesh.CharacteristicLengthFactor = 1.5;
// Corners of the domain
Point(0) = { 0, 0, 0};
Point(1) = { 1, 0, 0};
Point(2) = { 1, 0.5, 0};
Point(3) = {0.5, 0.5, 0};
Point(4) = {0.5, 1, 0};
Point(5) = { 0, 1, 0};
// Boundary edges of the domain
Line(1) = {0,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,0};
// Boundary of the domain
Line Loop(1) = {1,2,3,4,5,6};
// The domain itself
Plane Surface(1) = {1};