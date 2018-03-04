List of things that need to be fixed:

 - Lot of variables defined in common.h. Remove unwanted.
 - Perhaps its not efficient to define pointers to large arrays
   - x,y,mass,etc
   in common.h. Remove them and pass them around.

 - Also need to cleanup control.c 