
Macro IntersectLines;

bunch_of_entities[] = BooleanFragments{ Line{larg1}; }{ Line{larg2}; };
Printf("intersection created lines ", bunch_of_entities[]);

  point_cluster[] = {};

  n_lines = #bunch_of_entities[];
  If(n_lines == 2)
    Printf("No intersection with provided lines ");
    f1[] = {bunch_of_entities[0]};
    f2[] = {bunch_of_entities[1]};
    p_int = -1;
  Else
    f1[] = {bunch_of_entities[0],bunch_of_entities[1]};
    f2[] = {bunch_of_entities[2],bunch_of_entities[3]};
    points[] = Boundary{Line{f1[0]};};
    p_int = points[1];
    //Delete{ Line{larg1,larg2}; }
  EndIf


Return
