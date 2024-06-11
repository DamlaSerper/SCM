# Switchable contact model (SCM)
SCM modification of LIGGGHTS-PFM code to represent complex porous boundaries and improved primitive geometry functions for fix wall/gran function. SCM enables the representation of repetitive porous structures without the use of meshing due to the use of new primitive wall ycylinder_finite_porous or yplane_finite_porous. Current version of this model is written for repeating square pores on cylinders and flat rectangles, but the source code could be modified to fit other repetitive porous structures. The SCM code also include modifications for including new primitives types for defining finite cylinders, circles and concentric circles.
Please find the journal article describing the SCM here:
The research data of the study above can be found at: https://github.com/DamlaSerper/SCM_Research_Data/tree/main

**Contributing Authors:**
- MSc. Damla Serper (Aalto Univeristy, Finland) '*corresponding author*'
- Dr. Kevin Hanley (University of Edinburgh, UK)

## How to use the code:
- Clone this directory into the LIGGGHTS-PFM src directory on your local computer.
- Here the written .gitignore file allows you to only stage the required files within your folder.  
- Copy and move fix_wall_gran.cpp, fix_wall_gran.h and primitive_wall_definitions.h files from SCM folder to upper src directory, replacing the original LIGGGHTS-PUBLIC versions.
- Recompile the LIGGGHTS-PFM.
- Here are some examples of how to add the SCM and new primitive types within LIGGGHTS .in script:
    
    ### You are intending to use SCM to represent a porous cylinder (primitive type: ycylinder_finite_porous, 11 arguments)
    - param[0] = radius of the cylinder
    - param[1] = first coordinate of center
    - param[2] = second coordinate of center
    - param[3] = np_ver (integer number of pores in the vertical (y) direction; must be >= 1)
    - param[4] = np_hor (integer number of pores in a horizontal x-z plane; must be >= 1)
    - param[5] = ylow (y coordinate of the bottom of the finite cylinder)
    - param[6] = yhigh (y coordinate of the top of the finite cylinder)
    - param[7] = fvoid_ver (fractional distance in vertical repeating cell which is void(pore))
    - param[8] = fvoid_hor (fractional angular distance in planar circumferential cell which is void)
    - param[9] = ftransition_ver (fractional distance in vertical repeating cell which is void + transition)
    - param[10] = ftransition_hor (fractional angular distance in planar circumferential cell which is void + transition)
    - For further explanation on these variables please refer to the article.
    - Example syntax: ```fix wall_name all wall/gran model hertz tangential history primitive type 1 ycylinder_finite_porous 0.069 0 0 3 18 -0.0628 0 0.020044 0.01742 0.020544 0.01792 shear x 157.0796```

    ### You are intending to use SCM to represent a porous flat rectangle (primitive type: yplane_finite_porous, 12 arguments)
    - param[0] = position of the y-plane (The y coordinate, meaning the rectangle will be on the x, z dimensions)
    - param[1] = first coordinate of the left-bottom end of the rectangle (x)
    - param[2] = other coordinate of the left-bottom end of the rectangle (z)
    - param[3] = np_ver (integer number of pores in the vertical (z) direction; must be >= 1)
    - param[4] = np_hor (integer number of pores in a horizontal (x); must be >= 1)
    - param[5] = A_ver (height of the rectangle)
    - param[6] = A_hor (length of the rectangle)
    - param[7] = fvoid_ver (fractional distance in vertical repeating cell which is void)
    - param[8] = fvoid_hor (fractional distance in horizontal repeating cell which is void)
    - param[9] = ftransition_ver (fractional distance in vertical repeating cell which is void + transition)
    - param[10] = ftransition_hor (fractional distance in horizontal repeating cell which is void + transition)
    - param[11] = shape (dummy double variable, give '0')
    - Example syntax: ```fix filtermesh all wall/gran model hertz tangential history primitive type 2 yplane_finite_porous 0.0 0.0 0.0 30.0 30.0 0.162 0.162 0.16407407407407 0.16407407407407 0.17185185185185 0.17185185185185 0.0```

    ### You are inteding to use finite cylinders as a shape (primitive type: ycylinder_finite, 5 arguments)
    - param[0] = radius of the cylinder
    - param[1] = first coordinate of center
    - param[2] = second coordinate of center
    - param[3] = ylow (y coordinate of the bottom of the finite cylinder)
    - param[4] = yhigh (y coordinate of the top of the finite cylinder)
    - Example syntax: ```fix wall_name all wall/gran model hertz tangential history primitive type 1 ycylinder_finite 0.017 0 0 -0.0033 0.1367```   
   
    ### You are intending to use finite circles as a shape (primitive type: yplane_circle_finite, 2 arguments)
    - param[0] = y coordinate of the center
    - param[1] = radius of the circle
    - Example syntax: ```fix bot_disk_prim_wall all wall/gran model hertz tangential history primitive type 3 yplane_circle_finite -0.0628 0.069```
   
    ### You are intending to use concentric circles as a shape (primitive type: yplane_concencircle_finite, 4 arguments)
    - param[0] = y coordinate of the center
    - param[1] = inner radius of the circle
    - param[2] = outer radius of the circle
    - param[3] = +1 for the region between two circles
    - Example syntax: ```fix top_lid_prim_wall all wall/gran model hertz tangential history primitive type 3 yplane_concencircle_finite 0 0.017 0.069 +1```

**- NOTE: Currently all the new primitives are defined along y axis, for other types please modify the code.**
**- NOTE: Shear option is only available with the cylinder primitives.**

