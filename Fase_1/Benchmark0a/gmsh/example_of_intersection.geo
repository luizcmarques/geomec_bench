

SetFactory("OpenCASCADE");

Include "macros/IntersectLines.geo";

// Begin for the script

Point(1) = {0, 0, 0};
Point(2) = {20, 0, 0};
Point(3) = {0, 3, 0};
Point(4) = {10, -4, 0};

l1=newl; Line(l1) = {3, 4};
l2=newl; Line(l2) = {1, 2};

Call IntersectLines;

//Line(7) = {5,6};