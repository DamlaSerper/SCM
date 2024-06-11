/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Philippe Seil (JKU Linz)
   Paul Kieckhefen (TUHH)
------------------------------------------------------------------------- */

#ifndef LMP_PRIMITIVE_WALL_DEFINITIONS
#define LMP_PRIMITIVE_WALL_DEFINITIONS

#include "vector_liggghts.h"
#include "math_extra_liggghts.h"

/*
 * Necessary steps to add new primitive walls:
 * (1) add an enum for your primitive to WallType, but insert it before NUM_WTYPE
 * (2) add a string that you want to use in your input script to wallString and
 *     the number of arguments the wall requires to numArgs
 * (3) implement distance and neighbor list build functions
 * (4) add them to the switch statements in chooseContactTemplate() and
 *     chooseNeighlistTemplate() located at the bottom of this file
 */

namespace LAMMPS_NS
{
  namespace PRIMITIVE_WALL_DEFINITIONS
  {

    enum WallType
    {
        XPLANE,
        YPLANE,
        ZPLANE,
        XCYLINDER,
        YCYLINDER,
        ZCYLINDER,
        GENERAL_PLANE,
        GENERAL_CYLINDER,
        GENERAL_CONE,
	YCYLINDER_FINITE_POROUS, // FILTER_MESH_WALLS (this is a height constrained finite cylinder with novel boundary representation method for pore representation)
	YCYLINDER_FINITE, // INLET_CYL, TOP_CYL, BOT_CYL (this is a height constrained finite cylinder)
	YPLANE_CIRCLE_FINITE, // BOT_DISK (inward) and TOP_CYL_TOP_DISK (inward) (this is a plane that is contrained to be the area within a circle if inward, excluding this area if outward) 
	YPLANE_CONCENCIRCLE_FINITE, // TOP_LID (inward) and CONE_TOP_DISK (inward) (this is a plane taht is constrained to be two concentric circles' overlapping area if inward, excluding this area if outward)
	YPLANE_FINITE_POROUS,
        NUM_WTYPE
    };

    static const char *wallString[] =
    {
        "xplane",
        "yplane",
        "zplane",
        "xcylinder",
        "ycylinder",
        "zcylinder",
        "general_plane",
        "general_cylinder",
        "general_cone",
	"ycylinder_finite_porous",
	"ycylinder_finite",
	"yplane_circle_finite",
	"yplane_concencircle_finite",
	"yplane_finite_porous"
    };

    static int numArgs[] =
    {
        1,
        1,
        1,
        3,
        3,
        3,
        6,
        7,
        8,
	11,
	5,
	2,
	4,
	12
    }; //~ The numbers immediately above (11, 5, 2, 4) are the number of ycylinder_finite_porous, ycylinder_finite, yplane_cylinder_finite and yplane_concencircle_finite

    /*
     * templates for the different wall primitives
     * each wall type needs to be coded in a template specialization
     * default implementation: no contact
     */
    template<WallType W>
    double resolveContactTemplate(double *x, double r, double *delta, double *param) {return 1.;}

    /*
     * default neighbor list template returns true --> if neighbor list distance function is not
     * coded explicitly, every particle will be added
     */
    template<WallType W>
    bool resolveNeighlistTemplate(double *x, double r, double treshold, double *param) {return true;}

    /*
     * declaration of choosing functions
     * definitions need to be after ALL definitions of resolveContactTemplate and resolveNeighlistTemplate
     */
    inline double chooseContactTemplate(double *x, double r, double *delta, double *param, WallType wType);
    inline bool chooseNeighlistTemplate(double *x, double r, double treshold, double *param, WallType wType);

/* ---------------------------------------------------------------------- */

    /*
     * x,y,z planes can be handled by a template with dimension as template parameter
     */

    template<int dim>
    struct Dim
    {
      static const int x = dim, y = (dim+1)%3, z = (dim+2)%3;
    };

    template<int dim>
    struct Plane : public Dim<dim>
    {
      typedef Dim<dim> d;
      static double resolveContact(double *pos, double r, double *delta, double *param)
      {
        if(pos[d::x] > *param){
          delta[d::x] = *param - pos[d::x]; delta[d::y] = 0.; delta[d::z] = 0.;
          return pos[d::x] - *param - r;
        } else {
          delta[d::x] = *param - pos[d::x]; delta[d::y] = 0.; delta[d::z] = 0.;
          return *param - pos[d::x] - r;
        }
      }
      static bool resolveSameSide(double *pos0, double *pos1, double *param)
      {
        if((pos0[d::x] > *param && pos1[d::x] > *param) ||
           (pos0[d::x] < *param && pos1[d::x] < *param)){
          return true;
        }
        return false;
      }
      static bool resolveNeighlist(double *pos, double r, double treshold, double *param)
      {
        double dMax = r + treshold;
        double dist = pos[d::x] - *param;
        double absdist = (dist > 0.0) ? dist : -dist;
        return (absdist <= dMax);
      }
    };

/* ---------------------------------------------------------------------- */
    /*
     * same holds for x,y,z cylinders
     * param[0] = radius
     * param[1] = first coordinate of center
     * param[2] = second coordinate of center
     */
    template<int dim>
    struct Cylinder : public Dim<dim>
    {
    public:

      typedef Dim<dim> d;
      static double calcRadialDistance(double *pos, double *param, double &dy, double &dz)
      {
        dy = pos[d::y]-param[1];
        dz = pos[d::z]-param[2];
        return sqrt(dy*dy+dz*dz);
      }
      static double calcRadialDistanceSquared(double *pos, double *param, double &dy, double &dz)
      {
        dy = pos[d::y]-param[1];
        dz = pos[d::z]-param[2];
        return (dy*dy+dz*dz);
      }

      static double resolveContact(double *pos, double r, double *delta, double *param)
      {
        double dx, dy,dz, fact;
        double dist = calcRadialDistance(pos,param,dy,dz);
        if(dist > *param){
          dx = dist - *param - r;
          fact = (dist - *param) / dist;
          delta[d::x] = 0.; delta[d::y] = -dy*fact; delta[d::z] = -dz*fact;

        } else{
          dx = *param - dist - r;
          fact = (*param - dist) / dist;
          delta[d::x] = 0.; delta[d::y] = +dy*fact; delta[d::z] = +dz*fact;
        }
        return dx;
      }
      static bool resolveSameSide(double *pos0, double *pos1, double *param)
      {
        double dy,dz;
        double distsq0 = calcRadialDistanceSquared(pos0,param,dy,dz);
        double distsq1 = calcRadialDistanceSquared(pos1,param,dy,dz);
        double rsq = param[0]*param[0];
        if((distsq0 > rsq && distsq1 > rsq) ||
           (distsq0 < rsq && distsq1 < rsq)){
          return true;
        }
        return false;
      }
      static bool resolveNeighlist(double *pos, double r, double treshold, double *param)
      {
        double dy,dz;
        double dMax = r + treshold;
        double dist = calcRadialDistance(pos,param,dy,dz) - *param;
        return (dMax < dist || -dMax < dist);
      }

    };
    
	/* ---------------------------------------------------------------------- */ 
	// Special template for rectangular filter mesh

       /*
       * param[0] = position of the plane (if y plane, the y coordinate, meaning the rectangle will ve on x,z dimensions)
       * param[1] = first coordinate of the left end of the rectangle
       * param[2] = other coordinate of the left end of the rectangle
       * param[3] = np_ver (integer number of pores in the vertical (z) direction; must be >= 1)
       * param[4] = np_hor (integer number of pores in a horizontal (x); must be >= 1)
       * param[5] = A_ver (height of the rectangle)
       * param[6] = A_hor (length of the rectangle)
       * param[7] = fvoid_ver (fractional distance in vertical repeating cell which is void)
       * param[8] = fvoid_hor (fractional distance in horizontal repeating cell which is void)
       * param[9] = ftransition_ver (fractional distance in vertical repeating cell which is void + transition)
       * param[10] = ftransition_hor (fractional distance in horizontal repeating cell which is void + transition)
       * param[11] = shape (rectangle or square, you can use square for circle as well) 
       */
    
    template<int dim>
    struct RectangularFilterMesh : public Dim<dim>
    {
      typedef Dim<dim> d;
      static double resolveContact(double *pos, double r, double *delta, double *param)
      {
        if(pos[d::x] > *param){
          delta[d::x] = *param - pos[d::x]; delta[d::y] = 0.; delta[d::z] = 0.;
          return pos[d::x] - *param - r;
        } else {
          delta[d::x] = pos[d::x] - *param; delta[d::y] = 0.; delta[d::z] = 0.;
          return pos[d::x] - *param - r;
        }
      }
      static bool resolveSameSide(double *pos0, double *pos1, double *param)
      {
        if((pos0[d::x] > *param && pos1[d::x] > *param) ||
           (pos0[d::x] < *param && pos1[d::x] < *param)){
          return true;
        }
        return false;
      }
      static bool resolveNeighlist(double *pos, double r, double treshold, double *param)
      {
        double dMax = r + treshold;
        double dist = pos[d::x] - *param;
        double absdist = (dist > 0.0) ? dist : -dist;
        return (absdist <= dMax);
      }
    };
// End of rectangular filter mesh
/* -------------------------------------------------------------------------*/
    

    /* ---------------------------------------------------------------------- */
      /*
       * Define a special template for the ycylinder_finite_porous.
       * Effectively the same as the cylinder template as the extra parameters
       aren't used in this file.
       * param[0] = radius
       * param[1] = first coordinate of center
       * param[2] = second coordinate of center
       * param[3] = np_ver (integer number of pores in the vertical (y) direction; must be >= 1)
       * param[4] = np_hor (integer number of pores in a horizontal x-z plane; must be >= 1)
       * param[5] = ylow (y coordinate of the bottom of the finite cylinder)
       * param[6] = yhigh (y coordinate of the top of the finite cylinder)
       * param[7] = fvoid_ver (fractional distance in vertical repeating cell which is void)
       * param[8] = fvoid_hor (fractional angular distance in planar circumferential cell which is void)
       * param[9] = ftransition_ver (fractional distance in vertical repeating cell which is void + transition)
       * param[10] = ftransition_hor (fractional angular distance in planar circumferential cell which is void + transition)
       */
      template<int dim>
      struct Filtermesh : public Dim<dim>
      {
      public:

        typedef Dim<dim> d;
        static double calcRadialDistance(double *pos, double *param, double &dy, double &dz)
        {
          dy = pos[d::y]-param[1];
          dz = pos[d::z]-param[2];
          return sqrt(dy*dy+dz*dz);
        }

	static double calcRadialDistanceSquared(double *pos, double *param, double &dy, double &dz)
	{
	  dy = pos[d::y]-param[1];
	  dz = pos[d::z]-param[2];
	  return (dy*dy+dz*dz);
	}

	static bool resolveSameSide(double *pos0, double *pos1, double *param)
	{
	  double dy,dz;
	  double distsq0 = calcRadialDistanceSquared(pos0,param,dy,dz);
	  double distsq1 = calcRadialDistanceSquared(pos1,param,dy,dz);
	  double rsq = param[0]*param[0];
	  if((distsq0 > rsq && distsq1 > rsq) ||
	     (distsq0 < rsq && distsq1 < rsq)){
	    return true;
	  }
	  return false;
	}
	
        static double resolveContact(double *pos, double r, double *delta, double *param)
        {
          double dx, dy,dz, fact;

          const double dist = calcRadialDistance(pos,param,dy,dz);
	  if (MathExtraLiggghts::compDouble(dist, 0.0)) {
	    delta[d::x] = 0.; delta[d::y] = 0.; delta[d::z] = 0.;
	    dx = 0.0;
	    return dx; // break for zero dist (avoid divide-by-zero)
          }
 
	  /*~ This part of the code needed modification. We assume that 
	    usually a sphere won't overlap a solid region of the cylinder
	    by more than its radius. In the transition regions, it is 
	    possible for overlaps larger than this to appear, and the code 
	    needs to calculate a contact force for the particle if a solid
	    cylinder were present.*/
	  dx = *param - dist - r;
	  fact = (*param - dist) / dist;
	  delta[d::x] = 0.; delta[d::y] = +dy*fact; delta[d::z] = +dz*fact;
          return dx;
        }
	
        static bool resolveNeighlist(double *pos, double r, double treshold, double *param)
        {
          double dy,dz;
          double dMax = r + treshold;
          double dist = calcRadialDistance(pos,param,dy,dz) - *param;
          return (dMax < dist || -dMax < dist);
        }
      };
    
/* ---------------------------------------------------------------------- */
    /*
     * a cone stump/frustrum / cone mantle without the bounding
     * circular areas
     * param[0] = radius 1 (= at base point)
     * param[1,2,3] base point on cylinder axis
     * param[4,5,6] end point on cylinder axis
     * param[7] radius 2 (= at end point)
     */
    struct GeneralCone
    {
      static double calcDistance(double *pos, double *param, double *dx, bool& inside)
      {
        const double small = 1e-20;

        double p[3], tAxis[3], nAxis[3];
        vectorSubtract3D(pos, &param[1], p);
        vectorSubtract3D(&param[4], &param[1], tAxis);

        // geometry
        const double a = vectorMag3D(tAxis);  // frustrum length
        const double b = param[7] - param[0]; // difference of radii
        const double intercept[2] = {param[0], 0};
        double slope[2] = {b, a};
        vectorNormalize2D(slope);

        // 2D plane spanned by cone axis and pos.
        // project problem into this plane
        double P[2];             // pos in 2D plane
        vectorScalarDiv3D(tAxis, a + small);
        P[1] = vectorDot3D(tAxis, p);
        vectorAddMultiply3D(p, tAxis, -P[1], nAxis);
        P[0] = vectorMag3D(nAxis);
        vectorScalarDiv3D(nAxis, P[0] + small);

        // project point onto the cone mantle
        double closest[2] = {0., 0.};
        vectorSubtract2D(P, intercept, closest);
        vectorAddMultiply2D
        (
          intercept,
          slope,
          vectorDot2D
          (
            slope,
            closest
          ),
          closest
        );

        // this is not an infinite cone, bound to rims
        if (closest[1] > a) // "above"
        {
          closest[0] = param[7];
          closest[1] = a;
        }
        else if (closest[1] <= 0) // below
        {
          closest[0] = param[0];
          closest[1] = 0.;
        }

        // difference between point in plane and closest point on cone
        vectorSubtract2D(closest, P, closest);
        inside = closest[0] >= 0.;

        // transform back to 3D for dx
        vectorCopy3D(tAxis, dx);
        vectorScalarMult3D(dx, closest[1]);
        vectorAddMultiply3D(dx, nAxis, closest[0], dx);

        return vectorMag2D(closest);
      }

      static double resolveContact(double *pos, const double r, double *delta, double *param)
      {
        bool inside;
        return calcDistance(pos, param, delta, inside)- r;
      }

      static bool resolveSameSide(double *pos0, double *pos1, double *param)
      {
        double dx[3];
        bool inside0, inside1;
        calcDistance(pos0, param, dx, inside0);
        calcDistance(pos1, param, dx, inside1);

        return inside0 == inside1;
      }

      static bool resolveNeighlist(double *pos, const double r, const double treshold, double *param)
      {
        double dx[3];
        bool inside;
        const double dMax = r + treshold;
        const double dist = calcDistance(pos, param, dx, inside);
        return dist <= dMax;
      }

    };

/* ---------------------------------------------------------------------- */

    /*
     * a plane is defined by a point and a vector perpendicular to the
     * plane. Parameters:
     * param[0,1,2] point on plane
     * param[3,4,5] vector perpendicular to the plane
     */

    struct GeneralPlane {
    public:
      // setParams needs to be called before any
      // computation takes place!!!
      static void setParams(double *param, double *p, double *n)
      {
        p[0] = param[0]; p[1] = param[1]; p[2] = param[2];
        n[0] = param[3]; n[1] = param[4]; n[2] = param[5];

        vectorNormalize3D(n);
      }
      // relies on p, n to be properly set
      // returns signed distance
      static double distToPlane(double *pos, double *p,  double *n)
      {
        double vec[3];
        vectorSubtract3D(pos,p,vec);
        return vectorDot3D(vec,n);
      }

      static double resolveContact(double *pos, double r, double *delta, double *param)
      {
        double p[3],n[3];
        setParams(param,p,n);

        double dist = distToPlane(pos,p,n);
        vectorCopy3D(n,delta);
        vectorScalarMult3D(delta,-dist);
        double absdist = (dist > 0.0) ? dist : -dist;
        return absdist - r;
      }

      static bool resolveSameSide(double *pos0, double *pos1, double *param)
      {
        double p[3],n[3];
        setParams(param,p,n);
        const double dist0 = distToPlane(pos0,p,n);
        const double dist1 = distToPlane(pos1,p,n);

        if((dist0 > 0.0 && dist1 > 0.0) ||
           (dist0 < 0.0 && dist1 < 0.0)){
          return true;
        }
        return false;
      }

      static bool resolveNeighlist(double *pos, double r, double treshold, double *param)
      {
        double p[3],n[3];
        setParams(param,p,n);

        double dMax = r + treshold;
        double dist = distToPlane(pos,p,n);
        double absdist = (dist > 0.0) ? dist : -dist;
        return (absdist <= dMax);
      }

    };

/* ---------------------------------------------------------------------- */
    /*
     * a cylinder is defined by a radius, a point on and a vector parallel to
     * cylinder axis
     * param[0] = radius
     * param[1,2,3] point on cylinder axis
     * param[4,5,6] vector parallel to cylinder axis
     */
    struct GeneralCylinder
    {
      static double calcRadialDistance(double *pos, double *param, double &dx, double &dy, double &dz)
      {
        return sqrt(calcRadialDistanceSquared(pos, param, dx, dy, dz));
      }
      static double calcRadialDistanceSquared(double *pos, double *param, double &dx, double &dy, double &dz)
      {
        double p[3],v[3],w[3],pb[3];
        vectorCopy3D(&param[1],p);
        vectorCopy3D(&param[4],v);

        vectorSubtract3D(pos,p,w);

        const double c1 = vectorDot3D(w,v);
        const double c2 = vectorDot3D(v,v);
        const double b = c1 / c2;

        vectorAddMultiply3D(p,v,b,pb);

        dx = pos[0] - pb[0];
        dy = pos[1] - pb[1];
        dz = pos[2] - pb[2];
        return (dx*dx + dy*dy + dz*dz);
      }

      static double resolveContact(double *pos, double r, double *delta, double *param)
      {
        double dx,dy,dz,fact,dr;
        double dist = calcRadialDistance(pos,param,dx,dy,dz);
        if(dist > *param){
          dr = dist - *param - r;
          fact = (dist - *param) / dist;
          delta[0] = -dx*fact; delta[1] = -dy*fact; delta[2] = -dz*fact;
         } else {
          dr = *param - dist - r;
          fact = (*param - dist) / dist;
          delta[0] = +dx*fact; delta[1] = +dy*fact; delta[2] = +dz*fact;
        }
        return dr;
      }

      static bool resolveSameSide(double *pos0, double *pos1, double *param)
      {
        double dx,dy,dz;
        double distsq0 = calcRadialDistanceSquared(pos0,param,dx,dy,dz);
        double distsq1 = calcRadialDistanceSquared(pos1,param,dx,dy,dz);
        double rsq = param[0]*param[0];
        if((distsq0 > rsq && distsq1 > rsq) ||
           (distsq0 < rsq && distsq1 < rsq)){
          return true;
        }
        return false;
      }

      static bool resolveNeighlist(double *pos, double r, double treshold, double *param)
      {
        double dx,dy,dz;
        double dMax = r + treshold;
        double dist = calcRadialDistance(pos,param,dx,dy,dz) - *param;
        return (dMax < dist || -dMax < dist);
      }

    };

/* ---------------------------------------------------------------------- */

    /*
     * functions to choose the correct template
     */

/* ---------------------------------------------------------------------- */

    inline double chooseContactTemplate(double *x, double r, double *delta, double *param, WallType wType)
    {
      //TODO: find a way to create switch statement automatically
      switch(wType){
      case XPLANE:
        return Plane<0>::resolveContact(x,r,delta,param);
      case YPLANE:
        return Plane<1>::resolveContact(x,r,delta,param);
      case ZPLANE:
        return Plane<2>::resolveContact(x,r,delta,param);
      case XCYLINDER:
        return Cylinder<0>::resolveContact(x,r,delta,param);
      case YCYLINDER:
        return Cylinder<1>::resolveContact(x,r,delta,param);
      case ZCYLINDER:
        return Cylinder<2>::resolveContact(x,r,delta,param);
      case GENERAL_PLANE:
        return GeneralPlane::resolveContact(x,r,delta,param);
      case GENERAL_CYLINDER:
        return GeneralCylinder::resolveContact(x,r,delta,param);
      case GENERAL_CONE:
        return GeneralCone::resolveContact(x,r,delta,param);
      case YCYLINDER_FINITE_POROUS: 
	      return Filtermesh<1>::resolveContact(x,r,delta,param);
      case YCYLINDER_FINITE: 
	      return Cylinder<1>::resolveContact(x,r,delta,param);
      case YPLANE_CIRCLE_FINITE: 
	      return Plane<1>::resolveContact(x,r,delta,param);
      case YPLANE_CONCENCIRCLE_FINITE: 
      	      return Plane<1>::resolveContact(x,r,delta,param);
      case YPLANE_FINITE_POROUS: 
	      return RectangularFilterMesh<1>::resolveContact(x,r,delta,param);
      default: // default: no contact
        return 1.;
      }
    }

    inline bool chooseSameSideTemplate(double *x0, double *x1, double *param, WallType wType)
    {
      //TODO: create switch statement automatically
      switch(wType){
      case XPLANE:
        return Plane<0>::resolveSameSide(x0,x1,param);
      case YPLANE:
        return Plane<1>::resolveSameSide(x0,x1,param);
      case ZPLANE:
        return Plane<2>::resolveSameSide(x0,x1,param);
      case XCYLINDER:
        return Cylinder<0>::resolveSameSide(x0,x1,param);
      case YCYLINDER:
        return Cylinder<1>::resolveSameSide(x0,x1,param);
      case ZCYLINDER:
        return Cylinder<2>::resolveSameSide(x0,x1,param);
      case GENERAL_PLANE:
        return GeneralPlane::resolveSameSide(x0,x1,param);
      case GENERAL_CYLINDER:
        return GeneralCylinder::resolveSameSide(x0,x1,param);
      case GENERAL_CONE:
        return GeneralCone::resolveSameSide(x0,x1,param);
      case YCYLINDER_FINITE_POROUS:
        return Filtermesh<1>::resolveSameSide(x0,x1,param);
      case YCYLINDER_FINITE: 
	return Cylinder<1>::resolveSameSide(x0,x1,param);
      case YPLANE_CIRCLE_FINITE: 
	return Plane<1>::resolveSameSide(x0,x1,param);
      case YPLANE_CONCENCIRCLE_FINITE:
        return Plane<1>::resolveSameSide(x0,x1,param);
      case YPLANE_FINITE_POROUS:
        return RectangularFilterMesh<1>::resolveSameSide(x0,x1,param);
      default: // default: same side
        return true;
      }
    }
    inline bool chooseNeighlistTemplate(double *x, double r, double treshold, double *param, WallType wType)
    {
      //TODO: create switch statement automatically
      switch(wType){
      case XPLANE:
        return Plane<0>::resolveNeighlist(x,r,treshold,param);
      case YPLANE:
        return Plane<1>::resolveNeighlist(x,r,treshold,param);
      case ZPLANE:
        return Plane<2>::resolveNeighlist(x,r,treshold,param);
      case XCYLINDER:
        return Cylinder<0>::resolveNeighlist(x,r,treshold,param);
      case YCYLINDER:
        return Cylinder<1>::resolveNeighlist(x,r,treshold,param);
      case ZCYLINDER:
        return Cylinder<2>::resolveNeighlist(x,r,treshold,param);
      case GENERAL_PLANE:
        return GeneralPlane::resolveNeighlist(x,r,treshold,param);
      case GENERAL_CYLINDER:
        return GeneralCylinder::resolveNeighlist(x,r,treshold,param);
      case GENERAL_CONE:
        return GeneralCone::resolveNeighlist(x,r,treshold,param);
      case YCYLINDER_FINITE_POROUS:
        return Filtermesh<1>::resolveNeighlist(x,r,treshold,param);
      case YCYLINDER_FINITE: 
	return Cylinder<1>::resolveNeighlist(x,r,treshold,param);
      case YPLANE_CIRCLE_FINITE: 
	return Plane<1>::resolveNeighlist(x,r,treshold,param);
      case YPLANE_CONCENCIRCLE_FINITE:
        return Plane<1>::resolveNeighlist(x,r,treshold,param);
      case YPLANE_FINITE_POROUS:
        return RectangularFilterMesh<1>::resolveNeighlist(x,r,treshold,param);
      default: // default value: every particle will be added to neighbor list
        return true;
      }
    }

    inline int chooseAxis(WallType wType)
    {
      //TODO: create switch statement automatically
      switch(wType){
      case XCYLINDER:
        return 0;
      case YCYLINDER:
        return 1;
      case ZCYLINDER:
        return 2;
      case YCYLINDER_FINITE_POROUS:
        return 1;
      case YCYLINDER_FINITE: 
	return 1;
      case YPLANE_CONCENCIRCLE_FINITE:
        return 1;
      case YPLANE_FINITE_POROUS:
        return 1;
      default:
        return -1;
      }
    }

    inline double chooseCalcRadialDistance(double *pos, double *param, double &dx, double &dy, double &dz, WallType wType)
    {
      //TODO: create switch statement automatically
      switch(wType){
      case XCYLINDER:
        dx=0.;
        return Cylinder<0>::calcRadialDistance(pos,param,dy,dz);
      case YCYLINDER:
        dy = 0.;
        return Cylinder<1>::calcRadialDistance(pos,param,dz,dx);
      case ZCYLINDER:
        dz = 0.;
        return Cylinder<2>::calcRadialDistance(pos,param,dx,dy);
      case GENERAL_CYLINDER:
        return GeneralCylinder::calcRadialDistance(pos,param,dx,dy,dz);
      case YCYLINDER_FINITE_POROUS:
        dy = 0.;
        return Filtermesh<1>::calcRadialDistance(pos,param,dx,dz);
      case YCYLINDER_FINITE:
	dy = 0.;
        return Cylinder<1>::calcRadialDistance(pos,param,dx,dz);
	
      default:
        dx = dy = dz = 0.;
        return -1.;
      }
    }

    inline int chooseNumArgs(char *style)
    {
        for(int i = 0; i < NUM_WTYPE; i++)
            if(strcmp(style,wallString[i]) == 0)
                return numArgs[i];
        return 0;
    }

    inline int numArgsPrimitiveWall(char *style)
    {
        return chooseNumArgs(style);
    }
  }
}

#endif /* PRIMITIVEWALLDEFINITIONS_H_ */
