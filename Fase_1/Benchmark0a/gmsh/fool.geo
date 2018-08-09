

SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0};
Point(2) = {20, 0, 0};
Point(3) = {0, 3, 0};
Point(4) = {10, -4, 0};

Line(1) = {3, 4};
Line(2) = {1, 2};

bunch_of_entities[] = BooleanFragments{ Line{1}; }{ Line{2}; };
Printf("intersection created lines ", bunch_of_entities[]);

  n_lines = #bunch_of_entities[];
  If(n_lines==2)
    Printf("No intersection with provided lines ",{f1,f2});
    f1[] = {bunch_of_entities[0]};
    f2[] = {bunch_of_entities[1]};
  Else
    f1[] = {bunch_of_entities[0],bunch_of_entities[1]};
    f2[] = {bunch_of_entities[2],bunch_of_entities[3]};
  EndIf

//Delete{ Line{1,2}; }
//Line(7) = {5,6};