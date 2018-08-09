

SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0};
Point(2) = {20, 0, 0};
Point(3) = {0, 3, 0};
Point(4) = {10, -4, 0};

Line(1) = {3, 4};
Line(2) = {1, 2};

bunch_of_entities[] = BooleanFragments{ Line{1}; }{ Line{2}; };
Printf("intersection created lines ", bunch_of_entities[]);

//Delete{ Line{1,2}; }
//Line(7) = {5,6};