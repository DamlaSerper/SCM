/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations 

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_wall_gran.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair_gran.h"
#include "fix_rigid.h"
#include "fix_mesh.h"
#include "fix_contact_history_mesh.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_contact_property_atom_wall.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "compute_pair_gran_local.h"
#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface_stress.h"
#include "tri_mesh.h"
#include "primitive_wall.h"
#include "primitive_wall_definitions.h"
#include "mpi_liggghts.h"
#include "neighbor.h"
#include "contact_interface.h"
#include "fix_property_global.h"
#include <vector>
#include "granular_wall.h"
#include <assert.h>
#include <string>
#include <sstream>

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  #include "math_extra_liggghts_nonspherical.h"
#endif

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LAMMPS_NS::PRIMITIVE_WALL_DEFINITIONS;
using namespace LIGGGHTS::Walls;
using namespace LIGGGHTS::ContactModels;

/*NL*/ #define DEBUGMODE_LMP_FIX_WALL_GRAN false //(20950 < update->ntimestep) //(comm->me ==1)
/*NL*/ #define DEBUG_LMP_FIX_FIX_WALL_GRAN_M_ID 774
/*NL*/ #define DEBUG_LMP_FIX_FIX_WALL_GRAN_P_ID 206

const double SMALL = 1e-12;
double direc = 0.0; // for YPLANE_CIRCLE_FINITE - if this value is minus 1, the area outside

  // modes for conduction contact area calaculation
  // same as in fix_heat_gran_conduction.cpp

  enum{ CONDUCTION_CONTACT_AREA_OVERLAP,
        CONDUCTION_CONTACT_AREA_CONSTANT,
        CONDUCTION_CONTACT_AREA_PROJECTION};

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    // wall/gran requires gran properties
    // sph not
    if (strncmp(style,"wall/gran",9) == 0 && (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag))
        error->fix_error(FLERR,this,"requires atom attributes radius, omega, torque");

    // defaults
    store_force_ = false;
    store_force_contact_ = false;
    stress_flag_ = false;
    n_FixMesh_ = 0;
    dnum_ = 0;
    skinDistance_ = 0.0;

    r0_ = 0.;

    shear_ = 0;
    shearDim_ = shearAxis_ = -1;
    vectorZeroize3D(shearAxisVec_);

    atom_type_wall_ = 1; // will be overwritten during execution, but other fixes require a value here

    // initializations
    fix_wallforce_ = 0;
    fix_wallforce_contact_ = 0;
    fix_rigid_ = NULL;
    heattransfer_flag_ = false;

    FixMesh_list_ = NULL;

    rebuildPrimitiveNeighlist_ = false;

    addflag_ = 0;
    cwl_ = NULL;

    computeflag_ = 1;

    meshwall_ = -1;

    Temp_wall = -1.;
    fixed_contact_area_ = 0.;
    Q = Q_add = 0.;

    area_calculation_mode_ = CONDUCTION_CONTACT_AREA_OVERLAP;

    // parse args
    //style = new char[strlen(arg[2])+2];
    //strcpy(style,arg[2]);

    iarg_ = 3;
    narg_ = narg;

    int nremaining = narg - 3;
    char ** remaining_args = &arg[3];

    int64_t variant = Factory::instance().selectVariant(style, nremaining, remaining_args);
    impl = Factory::instance().create(style, variant, lmp, this);

    iarg_ = narg - nremaining;

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        /*NL*/// if (screen) fprintf(screen,"arg %s\n",arg[iarg_]);

        if (strcmp(arg[iarg_],"primitive") == 0) {
           iarg_++;
           meshwall_ = 0;

           if (meshwall_ == 1)
             error->fix_error(FLERR,this,"'mesh' and 'primitive' are incompatible, choose either of them");

           if (strcmp(arg[iarg_++],"type"))
             error->fix_error(FLERR,this,"expecting keyword 'type'");
           atom_type_wall_ = force->inumeric(FLERR,arg[iarg_++]);
           if (atom_type_wall_ < 1 || atom_type_wall_ > atom->ntypes)
             error->fix_error(FLERR,this,"1 <= type <= max type as defined in create_box'");

           int n_primitives_ = 1;

           if (strcmp(arg[iarg_],"num") == 0) {
             ++iarg_;
             n_primitives_ = atoi(arg[iarg_++]);
           }

           for(int iWall = 0; iWall < n_primitives_; ++iWall) {

             char *wallstyle = arg[iarg_++];
             nPrimitiveArgs = PRIMITIVE_WALL_DEFINITIONS::numArgsPrimitiveWall(wallstyle); //~ Made this available outside the constructor
             /*NL*/// if (screen) fprintf(screen,"nPrimitiveArgs %d\n",nPrimitiveArgs);

             if(narg-iarg_ < nPrimitiveArgs)
              error->fix_error(FLERR,this,"not enough arguments for primitive wall");

             double * argVec = new double[nPrimitiveArgs];
             for(int i=0;i<nPrimitiveArgs;i++)
             {
               /*NL*/ //if (screen) fprintf(screen,"primitive arg %s\n",arg[iarg_]);
               argVec[i] = force->numeric(FLERR,arg[iarg_++]);
             }

	     //~ Add code for new primitives
	     if (nPrimitiveArgs != 1 && nPrimitiveArgs != 3) {
	       double tol = 1.0e-10; //~ A tiny tolerance for checking floating-point numbers
	       if (atom->superquadric_flag) error->fix_error(FLERR,this,"The switchable contact model implementation isn't complete for superquadrics");
	   
	       if (nPrimitiveArgs > 4) {//~ cylinder primitives 
		 
   
		 if (nPrimitiveArgs != 12) {
		 rcylinder = argVec[0];
		 (nPrimitiveArgs == 11) ? ylow = argVec[5] : ylow = argVec[3];
		 (nPrimitiveArgs == 11) ? yhigh = argVec[6] : yhigh = argVec[4];
		 if (fabs(argVec[1]) > tol || fabs(argVec[2]) > tol) error->fix_error(FLERR,this,"Centreline of cylinder in FixWallGran must be along the y-axis");	     
		 if (rcylinder <= 0.0) error->fix_error(FLERR,this,"Cylinder radius in FixWallGran must be positive");
		 if (ylow >= yhigh) error->fix_error(FLERR,this,"ylow must be less than yhigh in FixWallGran");
		 }
		 
if (nPrimitiveArgs == 12) {//~ yplane_finite_porous
		C_y = argVec[0];
		C_x = argVec[1]; // Left end corner (x, hor)
		C_z = argVec[2]; // Left end corner (z, ver)
		np_ver = static_cast<int> (argVec[3]);
		np_hor = static_cast<int> (argVec[4]);
		A_ver = argVec[5];
		A_hor = argVec[6];
		fvoid_ver = argVec[7];
		fvoid_hor = argVec[8];
		ftransition_ver = argVec[9];
		ftransition_hor = argVec[10];
		
		
		if (np_ver <= 0 || np_hor <= 0) error->fix_error(FLERR,this,"np_ver and np_hor in FixWallGran must be positive integers");
		if (fvoid_ver < 0.0 || fvoid_hor < 0.0) error->fix_error(FLERR,this,"fvoid_ver and fvoid_hor in FixWallGran must be between 0.0 and 1.0");
		if (fvoid_ver > ftransition_ver || fvoid_hor > ftransition_hor) error->fix_error(FLERR,this,"fvoid_X cannot be larger than ftransition_X in FixWallGran");
		if (ftransition_ver > 1.0 || ftransition_hor > 1.0) error->fix_error(FLERR,this,"ftransition_ver and ftransition_hor in FixWallGran must be between 0.0 and 1.0");
		
		fsolid_ver = 1.0 - ftransition_ver + fvoid_ver;
		fsolid_hor = 1.0 - ftransition_hor + fvoid_hor;
		if (fsolid_ver < ftransition_ver || fsolid_hor < ftransition_hor) error->fix_error(FLERR,this,"Error in the specification of f* in FixWallGran");
		
		//~ Calculate some constants for use later in the post_force_primitive function
		Lv = A_ver / static_cast<double> (np_ver); //~ Height of a vertical repeating cell containing {void, transition 1, solid, transition 2}
		Lh = A_hor / static_cast<double> (np_hor); //~ Length of a horizontal repeating cell containing {void, transition 1, solid, transition 2}
		fdv = Lv*(1.0 - fsolid_ver); //~ Vertical height of transition region in distance units
		fdh = Lh*(1.0 - fsolid_hor); //~ Horizontal height of transition region in distance units
		
		/*~ Define a critical particle radius above which a particle can never even contact a
		transition/void region, so the contact is simply taken as contact with a solid cylinder.
		Take this as the total distance spanned by the void + 2 transition regions.*/
		double criticalradiusver = Lv*fvoid_ver + 2.0*fdv;
		double criticalradiushor = Lh*fvoid_hor + 2.0*fdh;
		(criticalradiusver > criticalradiushor) ? criticalradius = criticalradiushor : criticalradius = criticalradiusver;
		} // End of finite porous rectangle

		 if (nPrimitiveArgs == 11) {//~ ycylinder_filtermesh
		   np_ver = static_cast<int> (argVec[3]);
		   np_hor = static_cast<int> (argVec[4]);
		   fvoid_ver = argVec[7];
		   fvoid_hor = argVec[8];
		   ftransition_ver = argVec[9];
		   ftransition_hor = argVec[10];
	     
		   if (np_ver <= 0 || np_hor <= 0) error->fix_error(FLERR,this,"np_ver and np_hor in FixWallGran must be positive integers");
	       
		   if (fvoid_ver < 0.0 || fvoid_hor < 0.0) error->fix_error(FLERR,this,"fvoid_ver and fvoid_hor in FixWallGran must be between 0.0 and 1.0");
		   if (fvoid_ver > ftransition_ver || fvoid_hor > ftransition_hor) error->fix_error(FLERR,this,"fvoid_X cannot be larger than ftransition_X in FixWallGran");
		   if (ftransition_ver > 1.0 || ftransition_hor > 1.0) error->fix_error(FLERR,this,"ftransition_ver and ftransition_hor in FixWallGran must be between 0.0 and 1.0");
	     
		   fsolid_ver = 1.0 - ftransition_ver + fvoid_ver;
		   fsolid_hor = 1.0 - ftransition_hor + fvoid_hor;
		   if (fsolid_ver < ftransition_ver || fsolid_hor < ftransition_hor) error->fix_error(FLERR,this,"Error in the specification of f* in FixWallGran");
	     
		   //~ Calculate some constants for use later in the post_force_primitive function
		   Lv = (yhigh - ylow) / static_cast<double> (np_ver); //~ Height of a vertical repeating cell containing {void, transition, solid, transition}
		   Lh = 2.0*M_PI / static_cast<double> (np_hor); //~ Angle in radians demarcating a planar circumferential cell containing {void, transition, solid, transition}
		   fdv = Lv*(1.0 - fsolid_ver); //~ Vertical height of transition region in distance units

		   /*~ Define a critical particle radius above which a particle can never even contact a
		     transition/void region, so the contact is simply taken as contact with a solid cylinder.
		     Take this as the total distance spanned by the void + 2 transition regions. For the
		     horizontal direction, this is taken as a chord length.*/
		   double criticalradiusver = Lv*fvoid_ver + 2.0*fdv;
		   double criticalradiushor = 2.0*rcylinder*sin(0.5*Lh*(1.0 - fsolid_hor + ftransition_hor));
		   (criticalradiusver > criticalradiushor) ? criticalradius = criticalradiushor : criticalradius = criticalradiusver;
		 }
	       } else if (nPrimitiveArgs == 2) {// finite circle planes
		 rcircle = argVec[1]; // Serving the role of the circle radius
		 if (rcircle < 0.0) {
	           direc = -1; // -1 for outside of circle, 1 for inside the circle
	           rcircle = rcircle*(-1);
		 }
		 else direc = 1;
		 if (rcircle <= 0.0) error->fix_error(FLERR,this,"Inlet pipe radius in FixWallGran must be positive");
	       } else  if (nPrimitiveArgs == 4) { // Area within two concentric circles or the area except their intersection
		 rcircle = argVec[1]; // Inner circle radius
		 rocircle = argVec[2]; // Outer circle radius
		 direc = argVec[3]; // Inwards or outwards
	       }
	       else error->fix_error(FLERR,this,"Error in definition of primitive wall");
	     }
	   
             bool setflag = false;
             for(int w=0;w<(int)PRIMITIVE_WALL_DEFINITIONS::NUM_WTYPE;w++)
             {
               /*NL*///if (screen) fprintf(screen,"WallString %s %s\n",PRIMITIVE_WALL_DEFINITIONS::wallString[w],wallstyle);
               if(strcmp(wallstyle,PRIMITIVE_WALL_DEFINITIONS::wallString[w]) == 0)
               {
                 primitiveWalls_.push_back(new PrimitiveWall(lmp,(PRIMITIVE_WALL_DEFINITIONS::WallType)w,nPrimitiveArgs,argVec));
                 setflag = true;
                 break;
               }
             }
             if(!setflag) error->fix_error(FLERR,this,"unknown primitive wall style");
             delete[] argVec;
           }
           hasargs = true;
        } else if (strcmp(arg[iarg_],"mesh") == 0) {
           hasargs = true;
           meshwall_ = 1;
           iarg_ += 1;
        } else if (strcmp(arg[iarg_],"store_force") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"store_force_contact") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_contact_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_contact_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force_contact_'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"n_meshes") == 0) {
          if (meshwall_ != 1)
             error->fix_error(FLERR,this,"have to use keyword 'mesh' before using 'n_meshes'");
          if (iarg_+2 > narg)
             error->fix_error(FLERR,this,"not enough arguments");
          n_FixMesh_ = atoi(arg[iarg_+1]);
          if(n_FixMesh_ < 1)
              error->fix_error(FLERR,this,"'n_meshes' > 0 required");
          hasargs = true;
          iarg_ += 2;
        } else if (strcmp(arg[iarg_],"meshes") == 0) {
          if (meshwall_ != 1)
             error->fix_error(FLERR,this,"have to use keyword 'mesh' before using 'meshes'");
          if(n_FixMesh_ == 0)
              error->fix_error(FLERR,this,"have to define 'n_meshes' before 'meshes'");
          if (narg < iarg_+1+n_FixMesh_)
              error->fix_error(FLERR,this,"not enough arguments");

          FixMesh_list_ = new FixMeshSurface*[n_FixMesh_];
          for(int i = 1; i <= n_FixMesh_; i++)
          {
              int f_i = modify->find_fix(arg[iarg_+i]);
              if (f_i == -1)
                  error->fix_error(FLERR,this,"could not find fix mesh id you provided");
              if (strncmp(modify->fix[f_i]->style,"mesh/surface",12))
                  error->fix_error(FLERR,this,"the fix belonging to the id you provided is not of type mesh");
              FixMesh_list_[i-1] = static_cast<FixMeshSurface*>(modify->fix[f_i]);

              if(FixMesh_list_[i-1]->trackStress())
                stress_flag_ = true;
              //NP create neighbor list for each mesh
          }
          hasargs = true;
          iarg_ += 1+n_FixMesh_;
        } else if (strcmp(arg[iarg_],"shear") == 0) {
          if (iarg_+3 > narg)
            error->fix_error(FLERR,this,"not enough arguments for 'shear'");
          if(primitiveWalls_.size() == 0)
            error->fix_error(FLERR,this,"have to define primitive wall before 'shear'. For mesh walls, please use fix move/mesh");

          if (strcmp(arg[iarg_+1],"x") == 0) shearDim_ = 0;
          else if (strcmp(arg[iarg_+1],"y") == 0) shearDim_ = 1;
          else if (strcmp(arg[iarg_+1],"z") == 0) shearDim_ = 2;
          else error->fix_error(FLERR,this,"illegal 'shear' dim");
          vshear_ = force->numeric(FLERR,arg[iarg_+2]);
          shear_ = 1;

	  //~ Update this for the moving ycylinder
	  if (nPrimitiveArgs == 11 && shearDim_ != 1) {
	    shearAxis_ = 1; //~ As returned by primitive_wall_definitions for ycylinder*
            shearAxisVec_[shearAxis_] = vshear_;
	  }
	  
          // update axis for cylinder etc if needed
          /* FIXME if(shearDim_ != primitiveWall_->axis())
          {
            shearAxis_ = primitiveWall_->axis();
            shearAxisVec_[shearAxis_] = vshear_;
          }*/

          hasargs = true;
          iarg_ += 3;
        } else if (strcmp(arg[iarg_],"temperature") == 0) {
            if (iarg_+1 >= narg)
              error->fix_error(FLERR,this,"not enough arguments for 'temperature'");
            Temp_wall = force->numeric(FLERR,arg[iarg_+1]);
            hasargs = true;
            iarg_ += 2;
        } else if(strcmp(arg[iarg_],"contact_area") == 0) {

          if(strcmp(arg[iarg_+1],"overlap") == 0)
            area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_OVERLAP;
          else if(strcmp(arg[iarg_+1],"projection") == 0)
            area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_PROJECTION;
          else if(strcmp(arg[iarg_+1],"constant") == 0)
          {
            if (iarg_+3 > narg)
                error->fix_error(FLERR,this,"not enough arguments for keyword 'contact_area constant'");
            area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_CONSTANT;
            fixed_contact_area_ = force->numeric(FLERR,arg[iarg_+2]);
            if (fixed_contact_area_ <= 0.)
                error->fix_error(FLERR,this,"'contact_area constant' value must be > 0");
            iarg_++;
          }
          else error->fix_error(FLERR,this,"expecting 'overlap', 'projection' or 'constant' after 'contact_area'");
          iarg_ += 2;
          hasargs = true;
        }
    }

    if(impl)
      impl->settings(narg - iarg_, &arg[iarg_]);

    // error checks

    if(meshwall_ == -1 && primitiveWalls_.size() == 0)
        error->fix_error(FLERR,this,"Need to use define style 'mesh' or 'primitive'");

    if(meshwall_ == 1 && !FixMesh_list_)
        error->fix_error(FLERR,this,"Need to provide the number and a list of meshes by using 'n_meshes' and 'meshes'");
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_create()
{
    if(strncmp(style,"wall/gran",9) != 0)
    {
      // case non-granular (sph)
      dnum_ = 0;
    }

    // register storage for wall force if required
    if(store_force_)
    {
          char *wallforce_name = new char[strlen(id)+1+6];
          strcpy(wallforce_name,"force_");
          strcat(wallforce_name,id);
          const char *fixarg[11];
          fixarg[0] = wallforce_name;
          fixarg[1] = "all";
          fixarg[2] = "property/atom";
          fixarg[3] = wallforce_name;
          fixarg[4] = "vector";
          fixarg[5] = "no";    // restart
          fixarg[6] = "no";    // communicate ghost
          fixarg[7] = "no";    // communicate rev
          fixarg[8] = "0.";
          fixarg[9] = "0.";
          fixarg[10] = "0.";
          modify->add_fix(11,const_cast<char**>(fixarg));
          fix_wallforce_ =
              static_cast<FixPropertyAtom*>(modify->find_fix_property(wallforce_name,"property/atom","vector",3,0,style));
          delete []wallforce_name;
   }

   //NP case primitve wall - register here
   if(store_force_contact_ && 0 == meshwall_)
   {
        const char *fixarg[19];
        char fixid[200],ownid[200];
        sprintf(fixid,"contactforces_%s",id);
        sprintf(ownid,"%s",id);
        fixarg[0]=fixid;
        fixarg[1]="all";
        fixarg[2]="contactproperty/atom/wall";
        fixarg[3]=fixid;
        fixarg[4]="6";
        fixarg[5]="fx";
        fixarg[6]="0";
        fixarg[7]="fy";
        fixarg[8]="0";
        fixarg[9]="fz";
        fixarg[10]="0";
        fixarg[11]="tx";
        fixarg[12]="0";
        fixarg[13]="ty";
        fixarg[14]="0";
        fixarg[15]="tz";
        fixarg[16]="0";
        fixarg[17]="primitive";
        fixarg[18]=ownid;
        modify->add_fix(19,const_cast<char**>(fixarg));
        fix_wallforce_contact_ = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
   }

   // create neighbor list for each mesh
   //NP also create contact history in case of dnum == 0
   //NP so to be able to exclude duplicates contacts for coplanar faces
   for(int i=0;i<n_FixMesh_;i++)
   {
       //NP order is important!!
       //NP neigh list must be built before contact history is refreshed
       FixMesh_list_[i]->createWallNeighList(igroup);
       FixMesh_list_[i]->createContactHistory(dnum());

       if(store_force_contact_)
         FixMesh_list_[i]->createMeshforceContact();
   }

   // contact history for primitive wall
   if(meshwall_ == 0 && dnum_ > 0)
   {
     assert(primitiveWallsHistory_.size() == 0);

     for(size_t i = 0; i < primitiveWalls_.size(); ++i) {
        ostringstream os;
        os << "history_" << id << "_" << i;
        string hist_name = os.str();
        const char **fixarg = new const char*[8+dnum_];
        fixarg[0] = hist_name.c_str();
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = hist_name.c_str();
        fixarg[4] = "vector";
        fixarg[5] = "yes"; // restart
        fixarg[6] = "no";  // communicate ghost
        fixarg[7] = "no";  // communicate rev
        for(int j = 8; j < 8+dnum_; ++j)
          fixarg[j] = "0.";
        modify->add_fix(8+dnum_,const_cast<char**>(fixarg));
        FixPropertyAtom* history = static_cast<FixPropertyAtom*>(modify->find_fix_property(hist_name.c_str(),"property/atom","vector",dnum_,0,style));
        primitiveWallsHistory_.push_back(history);
        delete []fixarg;
     }
   }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::pre_delete(bool unfixflag)
{
    if(unfixflag && store_force_)
        modify->delete_fix(fix_wallforce_->id);
    if(unfixflag && primitiveWallsHistory_.size() > 0) {
      for(size_t i = 0; i < primitiveWallsHistory_.size(); ++i){
        modify->delete_fix(primitiveWallsHistory_[i]->id);
        primitiveWallsHistory_[i] = NULL;
      }
      primitiveWallsHistory_.clear();
    }

   if(unfixflag && store_force_contact_) {
     modify->delete_fix(fix_wallforce_contact_->id);
   }

    if(unfixflag)
    {
       for(int i=0;i<n_FixMesh_;i++)
       {
           //NP remove contact history as new fix wall gran might have a
           //NP different contact history needs
           FixMesh_list_[i]->deleteWallNeighList();
           FixMesh_list_[i]->deleteContactHistory();
       }
    }
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
    if(primitiveWalls_.size() > 0) {
      for(size_t i = 0; i < primitiveWalls_.size(); ++i) {
        delete primitiveWalls_[i];
        primitiveWalls_[i] = NULL;
      }
      primitiveWalls_.clear();
    }
    if(FixMesh_list_) delete []FixMesh_list_;
    delete impl;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::setmask()
{
    int mask = 0;
    mask |= PRE_NEIGHBOR;
    mask |= PRE_FORCE;
    mask |= POST_FORCE;
    mask |= POST_FORCE_RESPA;
    return mask;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::min_type() const
{
    return atom_type_wall_;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::max_type() const
{
    return atom_type_wall_;
}

/* ---------------------------------------------------------------------- */

PrimitiveWall* FixWallGran::primitiveWall()
{ return primitiveWalls_[0]; } // FIXME

/* ---------------------------------------------------------------------- */

void FixWallGran::init()
{
    dt_ = update->dt;

    // case granular
    if(strncmp(style,"wall/gran",9) == 0)
    {
        // check if a fix rigid is registered - important for damp
        fix_rigid_ = static_cast<FixRigid*>(modify->find_fix_style_strict("rigid",0));

        if (strcmp(update->integrate_style,"respa") == 0)
          nlevels_respa_ = ((Respa *) update->integrate)->nlevels;

        if(impl)
          impl->init_granular();
        else
        {
          // init for derived classes
          init_granular();
        }

        // disallow more than one wall of non-rimitive style
        //NP because contact history is stored in the mesh
        //NP and mesh can hold only one history
        if(is_mesh_wall())
        {
            int nfix = modify->n_fixes_style("wall/gran");
            for (int ifix = 0; ifix < nfix; ifix++)
            {
                FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
                if (fwg == this) continue;
                if (fwg->is_mesh_wall())
                    error->fix_error(FLERR,this,"More than one wall of type 'mesh' is not supported");
            }
        }
    }
    //NP out = fopen("shearhistory","w");
}

/* ---------------------------------------------------------------------- */

void FixWallGran::setup(int vflag)
{
    if (strstr(update->integrate_style,"verlet"))
    {
      pre_neighbor();
      pre_force(vflag);
      post_force(vflag);
    }
    else
    {
      ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa_-1);
      post_force_respa(vflag,nlevels_respa_-1,0);
      ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa_-1);
    }

    //NP doing this here because deltan_ratio is set in init() of fix heat/gran
    init_heattransfer();
}

/* ----------------------------------------------------------------------
   neighbor list via fix wall/gran is only relevant for primitive walls
------------------------------------------------------------------------- */

void FixWallGran::pre_neighbor()
{
    rebuildPrimitiveNeighlist_ = (primitiveWalls_.size() > 0);
}

void FixWallGran::pre_force(int vflag)
{
    const double halfskin = neighbor->skin*0.5;
    const int nlocal = atom->nlocal;

    x_ = atom->x;
    radius_ = atom->radius;
    cutneighmax_ = neighbor->cutneighmax;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    if(atom->superquadric_flag) {
      quat_ = atom->quaternion;
      shape_ = atom->shape;
      blockiness_ = atom->blockiness;
    }
#endif

    /*NL*/// if (screen) fprintf(screen,"cutneighmax_ %f, radius is %s, rebuild primitive %s\n",cutneighmax_,radius_?"notnull":"null",rebuildPrimitiveNeighlist_?"yes":"no");
    /*NL*/// error->all(FLERR,"end");

    // build neighlist for primitive walls
    //NP as in neighlist/mesh - hack for sph
    if(rebuildPrimitiveNeighlist_) {
      for(size_t i = 0; i < primitiveWalls_.size(); ++i) {
        primitiveWalls_[i]->buildNeighList(radius_ ? halfskin:(r0_+halfskin),x_,radius_,nlocal);
      }
    }

    rebuildPrimitiveNeighlist_ = false;
}

/* ----------------------------------------------------------------------
   force on each atom calculated via post_force
   called via verlet
------------------------------------------------------------------------- */

void FixWallGran::post_force(int vflag)
{
    computeflag_ = 1;
    shearupdate_ = 1;
    if (update->setupflag) shearupdate_ = 0;
    addflag_ = 0;

    post_force_wall(vflag);
}

/* ----------------------------------------------------------------------
   force on each atom calculated via post_force
   called via compute wall/gran
------------------------------------------------------------------------- */

void FixWallGran::post_force_pgl()
{
    computeflag_ = 0;
    shearupdate_ = 0;
    addflag_ = 1;

    post_force_wall(0);
}

/* ----------------------------------------------------------------------
   post_force
------------------------------------------------------------------------- */

void FixWallGran::post_force_wall(int vflag)
{
  // set pointers and values appropriately
  nlocal_ = atom->nlocal;
  x_ = atom->x;
  f_ = atom->f;
  radius_ = atom->radius;
  rmass_ = atom->rmass;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  if(atom->superquadric_flag) {
    quat_ = atom->quaternion;
    shape_ = atom->shape;
    blockiness_ = atom->blockiness;
  }
#endif

  if(fix_rigid_)
  {
      body_ = fix_rigid_->body;
      masstotal_ = fix_rigid_->masstotal;
  }

  if(fix_wallforce_)
    wallforce_ = fix_wallforce_->array_atom;

  cutneighmax_ = neighbor->cutneighmax;

  if(nlocal_ && !radius_ && r0_ == 0.)
    error->fix_error(FLERR,this,"need either per-atom radius or r0_ being set");

  if(store_force_)
  {
      for(int i = 0; i < nlocal_; i++)
      {
          vectorZeroize3D(wallforce_[i]);
      }
  }

  if(meshwall_ == 1)
    post_force_mesh(vflag);
  else
    post_force_primitive(vflag);

  if(meshwall_ == 0 && store_force_contact_)
    fix_wallforce_contact_->do_forward_comm();

  if(meshwall_ == 1 && store_force_contact_)
  {
    for(int imesh = 0; imesh < n_FixMesh_; imesh++)
        FixMesh_list_[imesh]->meshforceContact()->do_forward_comm();
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force_respa(int vflag, int ilevel, int iloop)
{
    if (ilevel == nlevels_respa_-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   post_force for mesh wall
------------------------------------------------------------------------- */

void FixWallGran::post_force_mesh(int vflag)
{
    //NP groupbit accounted for via neighlist/mesh

    // contact properties
    double v_wall[3],bary[3];
    double delta[3],deltan;
    const int nlocal = atom->nlocal;
    const double contactDistanceMultiplier = neighbor->contactDistanceFactor - 1.0;

    CollisionData cdata;
    cdata.is_wall = true;
    cdata.computeflag = computeflag_;
    cdata.shearupdate = shearupdate_;

    /*NL*/// if(comm->me == 3 && update->ntimestep == 3735 && screen)
    /*NL*///   fprintf(screen,"proc 3 start\n");

    for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
    {
      FixContactHistoryMesh *fix_contact = FixMesh_list_[iMesh]->contactHistory();
      // mark all contacts for deletion at this point
      //NP all detected contacts will be un-marked by fix_contact->handleContact()

      if(fix_contact) {
        FixNeighlistMesh * meshNeighlist = static_cast<FixNeighlistMesh*>(FixMesh_list_[iMesh]->meshNeighlist());
        fix_contact->resetDeletionPage(0);

        std::vector<int> & particle_indices = meshNeighlist->particle_indices;

        // mark all contacts for deletion at this point
        //NP all detected contacts will be un-marked by fix_contact->handleContact()
        for(std::vector<int>::iterator it = particle_indices.begin(); it != particle_indices.end(); ++it) {
          fix_contact->markForDeletion(0, *it);
        }
      }

      //NP extremely dirty; this is that FixWallGranBase can use fix_wallforce_contact_
      //NP pointer for both primitive and mesh case
      if(store_force_contact_)
        fix_wallforce_contact_ = FixMesh_list_[iMesh]->meshforceContact();
    }
    //if (screen) fprintf(screen, "[%d] Mark All took %g seconds.\n", tid, MPI_Wtime() - markAllTime);

    for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
    {
      TriMesh *mesh = FixMesh_list_[iMesh]->triMesh();
      const int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
      FixContactHistoryMesh * const fix_contact = FixMesh_list_[iMesh]->contactHistory();

      // get neighborList and numNeigh
      FixNeighlistMesh * meshNeighlist = static_cast<FixNeighlistMesh*>(FixMesh_list_[iMesh]->meshNeighlist());
      std::vector<int>::iterator b = meshNeighlist->particle_indices.begin();
      std::vector<int>::iterator e = meshNeighlist->particle_indices.end();
      if(b == e) continue; // nothing to do (no neighbors)

      vectorZeroize3D(v_wall);
      MultiVectorContainer<double,3,3> *vMeshC = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");

      cdata.jtype = FixMesh_list_[iMesh]->atomTypeWall();

      // moving mesh
      if(vMeshC)
      {
        double ***vMesh = vMeshC->begin();

        // loop owned and ghost triangles
        for(int iTri = 0; iTri < nTriAll; iTri++)
        {
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();

          for(int iCont = 0; iCont < numneigh; iCont++)
          {
            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            int idTri = mesh->id(iTri);

#ifdef SUPERQUADRIC_ACTIVE_FLAG
            if(atom->superquadric_flag)
            {
              Superquadric particle(x_[iPart], quat_[iPart], shape_[iPart], blockiness_[iPart]);

              if(mesh->sphereTriangleIntersection(iTri, radius_[iPart], x_[iPart])) //check for Bounding Sphere-triangle intersection
              {
                deltan = mesh->resolveTriSuperquadricContact(iTri, delta, cdata.contact_point, particle, bary);
              }
              else
              {
                deltan = LARGE_TRIMESH;
              }

              cdata.is_non_spherical = true; //by default it is false
            }
            else
            {
              deltan = mesh->resolveTriSphereContactBary(iPart,iTri,radius_ ? radius_[iPart]:r0_,x_[iPart],delta,bary);
            }
#else
            deltan = mesh->resolveTriSphereContactBary(iPart,iTri,radius_ ? radius_[iPart]:r0_ ,x_[iPart],delta,bary);
#endif

            if(deltan > cutneighmax_) continue;

            bool intersectflag = (deltan <= 0);

            if(intersectflag || (radius_ && deltan < contactDistanceMultiplier*radius_[iPart]))
            {
              if(fix_contact && ! fix_contact->handleContact(iPart,idTri,cdata.contact_history)) continue;

              for(int i = 0; i < 3; i++)
                v_wall[i] = (bary[0]*vMesh[iTri][0][i] + bary[1]*vMesh[iTri][1][i] + bary[2]*vMesh[iTri][2][i]);

              cdata.i = iPart;
              cdata.deltan = -deltan;
              cdata.delta[0] = -delta[0];
              cdata.delta[1] = -delta[1];
              cdata.delta[2] = -delta[2];
              post_force_eval_contact(cdata, intersectflag, v_wall,iMesh,FixMesh_list_[iMesh],mesh,iTri);
            }

          }
        }
      }
      // non-moving mesh - do not calculate v_wall, use standard distance function
      else
      {
        // loop owned and ghost particles
        for(int iTri = 0; iTri < nTriAll; iTri++)
        {
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();

          for(int iCont = 0; iCont < numneigh; iCont++)
          {
            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            int idTri = mesh->id(iTri);

#ifdef SUPERQUADRIC_ACTIVE_FLAG
            if(atom->superquadric_flag)
            {
              Superquadric particle(x_[iPart], quat_[iPart], shape_[iPart], blockiness_[iPart]);

              if(mesh->sphereTriangleIntersection(iTri, radius_[iPart], x_[iPart])) //check for Bounding Sphere-triangle intersection
              {
                deltan = mesh->resolveTriSuperquadricContact(iTri, delta, cdata.contact_point, particle);
              }
              else
              {
                deltan = LARGE_TRIMESH;
              }

              cdata.is_non_spherical = true; //by default it is false
            }
            else
            {
              deltan = mesh->resolveTriSphereContact(iPart,iTri,radius_ ? radius_[iPart]:r0_,x_[iPart],delta);
            }
#else
            deltan = mesh->resolveTriSphereContact(iPart,iTri,radius_ ? radius_[iPart]:r0_,x_[iPart],delta);
#endif

            if(deltan > cutneighmax_) continue;

            bool intersectflag = (deltan <= 0);

            //NP hack for SPH
            if(intersectflag || (radius_ && deltan < contactDistanceMultiplier*radius_[iPart]))
            {
              //NP continue in case already have a contact with a coplanar face
              if(fix_contact && ! fix_contact->handleContact(iPart,idTri,cdata.contact_history)) continue;

              cdata.i = iPart;
              cdata.deltan = -deltan;
              cdata.delta[0] = -delta[0];
              cdata.delta[1] = -delta[1];
              cdata.delta[2] = -delta[2];
              post_force_eval_contact(cdata, intersectflag, v_wall,iMesh,FixMesh_list_[iMesh],mesh,iTri);
            }
          }
        }
      }
    }

  // clean-up contacts
  //NP i.e. delete all which have not been detected this time-step
  //NP have to do this here: values cannot be deleted before
  //NP since must be available to copy contact history for coplanar neighs
  //double cleanupTime = MPI_Wtime();
  for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
  {
    FixContactHistoryMesh *fix_contact = FixMesh_list_[iMesh]->contactHistory();
      // clean-up contacts
      //NP i.e. delete all which have not been detected this time-step
      if(fix_contact) {
        FixNeighlistMesh * meshNeighlist = static_cast<FixNeighlistMesh*>(FixMesh_list_[iMesh]->meshNeighlist());
        std::vector<int> & particle_indices = meshNeighlist->particle_indices;

        for(std::vector<int>::iterator it = particle_indices.begin(); it != particle_indices.end(); ++it) {
          fix_contact->cleanUpContact(*it);
        }
      }
  }
}

/* ----------------------------------------------------------------------
   post_force for primitive wall
------------------------------------------------------------------------- */

void FixWallGran::post_force_primitive(int vflag)
{
  int *mask = atom->mask;

  CollisionData cdata;
  cdata.is_wall = true;
  cdata.computeflag = computeflag_;
  cdata.shearupdate = shearupdate_;
  cdata.jtype = atom_type_wall_;
  const double contactDistanceMultiplier = neighbor->contactDistanceFactor - 1.0;

  // contact properties
  double v_wall[] = {0.,0.,0.};

  // if shear, set velocity accordingly
  if (shear_) v_wall[shearDim_] = vshear_; // FIXME

  for(size_t iWall = 0; iWall < primitiveWalls_.size(); ++iWall) {
    PrimitiveWall * primitiveWall = primitiveWalls_[iWall];
    double **c_history = 0;

    if(dnum() > 0) {
      c_history = primitiveWallsHistory_[iWall]->array_atom;
    }

    // loop neighbor list
    int *neighborList;
    const int nNeigh = primitiveWall->getNeighbors(neighborList);

    for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++)
    {
      int iPart = *neighborList;

      if(!(mask[iPart] & groupbit)) continue;

      double delta[3]={};
      double deltan = primitiveWall->resolveContact(x_[iPart],radius_?radius_[iPart]:r0_,delta);

      if(deltan > cutneighmax_) continue;
      
      if(deltan <= 0 || deltan < contactDistanceMultiplier*radius_[iPart])
	{
	  // spheres
	  bool intersectflag = (deltan <= 0);

	  //~ Moved these two lines up from later in the function
	  cdata.i = iPart;
	  cdata.contact_history = c_history ? c_history[iPart] : NULL;
	
	  if (nPrimitiveArgs == 11) {
	    /*~ The spherical particle would be in contact with a solid cylinder. 
	      Now determine, based on its location in the vertical (y) direction
	      and in the horizontal plane, whether the sphere is contacting a 
	      solid region of the mesh or a void. In here, we assume that the 
	      centreline of the centrifuge is the y-axis (though this can
	      be moved easily enough if required).*/
	    int deleteparticle, trueoverlapflag = 0;
	    double radialdist, fvtemp, fv, twopi = 2.0*M_PI;
	  
	    radialdist = sqrt(x_[iPart][0]*x_[iPart][0] + x_[iPart][2]*x_[iPart][2]); //~ Radial distance of particle from basket centreline

	    //~ A new addition to compute v_wall for the contact force calculation(s)
	    if (shear_) {
	      v_wall[0] = vshear_*x_[iPart][2]/radialdist;
	      v_wall[1] = 0.0;
	      v_wall[2] = -vshear_*x_[iPart][0]/radialdist;
	    }
	    
	    fvtemp = (x_[iPart][1] - ylow)/Lv; //~ x_[iPart][1] is the y-coordinate of the particle
	    fv = fmod(fvtemp, 1.0); //~ Fractional height of the particle (0.0-1.0) in the vertical repeating cell

	    /*~ Determine if this fv corresponds to void (0), solid (1) or 
	      transition (2 or 3) regions. 2 is used for the first transition
	      from void to solid with increasing fv; 3 is used for the
	      second transition from solid to void with increasing fv*/
	    int fregion[2]; //~ [fv, fh]
	    if (radius_[iPart] >= criticalradius) { 
	      fregion[0] = 1;
	    } else {
	      fregion[0] = region_finder(fvoid_ver, ftransition_ver, fsolid_ver, fv);
	    }
	    /*~ The location in the horizonal plane is trickier as the particle's
	      angle needs to be rotated to account for the angular rotation of the
	      centrifuge basket from its fixed reference configuration*/
	    double unrotatedangle, rotatedangle, rottemp, fhtemp, fh;
	    if (fregion[0] != 1) {
	      //~ Only need to do these calculations if this location is not definitely solid
	      unrotatedangle = atan2(x_[iPart][2], x_[iPart][0]); //~ Angle in range [-pi, pi]
	      if (unrotatedangle <= 0.0) unrotatedangle += twopi; //~ Angle in range [0, 2*pi]
	      rottemp = unrotatedangle - vshear_ * update->ntimestep * update->dt/rcylinder; //~ vshear_ /rcylinder is the angular velocity of the rotating basket
	      rotatedangle = fmod(rottemp, twopi); //~ In the range [-2*pi, 2*pi]
	      if (rotatedangle <= 0.0) rotatedangle += twopi; //~ Angle in range [0, 2*pi]
	      fhtemp = rotatedangle/Lh;
	      fh = fmod(fhtemp, 1.0); //~ Fractional arc distance of the particle (0.0-1.0) in the planar circumferential cell, measured from the positive x-axis
	      fregion[1] = region_finder(fvoid_hor, ftransition_hor, fsolid_hor, fh);
	    }

	    //~ Example of printing out diagnostic information
	    // fprintf(screen,"fh %1.4f fv %1.4f fregion[0] %i fregion[1] %i on timestep " BIGINT_FORMAT "\n",fh,fv,fregion[0],fregion[1],update->ntimestep);
          
	    /*~ Deal with the case of a particle contacting a void region in BOTH 
	      dimensions: zero shear forces, check for the particle being a distance
	      approximately > its radius outside of the basket, and continue*/
	    if (fregion[0] == 0 && fregion[1] == 0) {
	      /*~ Here is the modification to the existing particle egress concept based solely
		on the location of the particle's centrepoint. Now we also take into consideration
		the location of the lowermost, uppermost, leftmost and rightmost (in a 2D
		projection sense) points on the spherical particle.

		Horizontal lines from the basket centreline through all four of these additional
		points must also lie within the void region in both dimensions for egress to be
		permitted. For any point that does not satisfy this, we apply a force to the
		particle which is for a contact with a 45-degree transition plane; thus free
		passage of the particle through the void is hindered.*/
	      int fregion_fourpoints[4]; //~ [lowermost, uppermost, negative angular direction, positive angular direction]
	      double fvtemp_lower, fv_lower, fvtemp_upper, fv_upper;
	      double incrementangle;
	      double rottemp_negative, rotatedangle_negative, fhtemp_negative, fh_negative;
	      double rottemp_positive, rotatedangle_positive, fhtemp_positive, fh_positive;
	    
	      fvtemp_lower = (x_[iPart][1] - radius_[iPart] - ylow)/Lv;
	      fv_lower = fmod(fvtemp_lower, 1.0); //~ Not a problem if fv_lower is slightly negative because it will still be mapped to void by region_finder
	      fregion_fourpoints[0] = region_finder(fvoid_ver, ftransition_ver, fsolid_ver, fv_lower);
 
	      fvtemp_upper = (x_[iPart][1] + radius_[iPart] - ylow)/Lv;
	      fv_upper = fmod(fvtemp_upper, 1.0);
	      fregion_fourpoints[1] = region_finder(fvoid_ver, ftransition_ver, fsolid_ver, fv_upper);
 
	      /*~ incrementangle is the angular amount that a particle's radius changes rotatedangle (found at
		the particle's centrepoint) by*/
	      incrementangle = atan2(radius_[iPart],radialdist);

	      rottemp_negative = rotatedangle - incrementangle;
	      rotatedangle_negative = fmod(rottemp_negative, twopi);
	      if (rotatedangle_negative <= 0.0) rotatedangle_negative += twopi;
	      fhtemp_negative = rotatedangle_negative/Lh;
	      fh_negative = fmod(fhtemp_negative, 1.0);
	      fregion_fourpoints[2] = region_finder(fvoid_hor, ftransition_hor, fsolid_hor, fh_negative);
	    
	      rottemp_positive = rotatedangle + incrementangle;
	      rotatedangle_positive = fmod(rottemp_positive, twopi);
	      if (rotatedangle_positive <= 0.0) rotatedangle_positive += twopi;
	      fhtemp_positive = rotatedangle_positive/Lh;
	      fh_positive = fmod(fhtemp_positive, 1.0);
	      fregion_fourpoints[3] = region_finder(fvoid_hor, ftransition_hor, fsolid_hor, fh_positive);

	      if (fregion_fourpoints[0] == 0 && fregion_fourpoints[1] == 0 && fregion_fourpoints[2] == 0 && fregion_fourpoints[3] == 0) {//~ If all are zero
		if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
		deleteparticle = outside_basket_check(radialdist, rcylinder, radius_[iPart], iPart);
		continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
	      } else {
		/*~ Else assume contact with a transition region to prevent free passage of the
		  particle out of the basket. The code below sets fregion[*] to denote contact
		  with one or more transition regions. Note the way it is done means that fregion[*]
		  can be 2, 3 or 5 with the '5' value indicating contact with both transition regions
		  simultaneously, i.e., a particle that is larger than a pore that can never egress*/
		if (fregion_fourpoints[0] > 0) fregion[0] += 3;
		if (fregion_fourpoints[1] > 0) fregion[0] += 2;
		if (fregion_fourpoints[2] > 0) fregion[1] += 3;
		if (fregion_fourpoints[3] > 0) fregion[1] += 2;
	      }
	    }

	    /*~ Pass on to a function if the particle contacts a transition 
	      region. Bypass if the particle contacts a solid region ('do 
	      nothing' in my table of actions). The vertical and horizontal 
	      (planar) directions are treated as independent for the purpose 
	      of calculating the modified forces.
     
	      For convenience, I have decided to take the physical form of 
	      these transitions as walls with 45 degree orientations rather
	      than cylinders as calculating the contact force between a sphere
	      and a cylinder is trickier than I anticipated!*/
	    ForceData i_forces;
	    ForceData j_forces;
	  
	    if (fregion[0] != 1 && fregion[1] != 1) {
	      if (fregion[0] >= 2 || fregion[1] >= 2) {
        
		//~ Check for the particle being outside of the basket in this case too
		deleteparticle = outside_basket_check(radialdist, rcylinder, radius_[iPart], iPart);
		if (deleteparticle == 1) continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
        
		double diffforce[3]; //~ The contact force vector due to an overlap of -deltan
		double magdiffforce = -1.0; //~ The magnitude of this contact force vector. Initialise at an unreasonable value
		double trueoverlap; //~ The overlap between the particle and the 45 degree plane
		double desiredmagforce; //~ The force magnitude to apply
              
		if (fregion[0] >= 2) {
		  double numeratorinalpha, transitionfrac;

		  while (fregion[0] >= 2) {
		    fregion[0] == 2 ? numeratorinalpha = fmod(ftransition_ver - fv, 1.0) : numeratorinalpha = fmod(fv - fsolid_ver, 1.0); //~ Start with transition region 2 if fregion[0] == 5
		    if (numeratorinalpha < 0.0) numeratorinalpha += 1.0; //~ Ensure always in the range [0.0, 1.0]
		    transitionfrac = 1.0 - fsolid_ver;
          
		    double alpha = fdv*numeratorinalpha/transitionfrac; //~ Distance units due to multiplication by fdv from constructor
		    double beta = rcylinder + alpha - radialdist;
		    double Lver = beta*M_SQRT1_2; //~ M_SQRT1_2 = 1/sqrt(2) = sin(45 degrees)
		    trueoverlap = radius_[iPart] - Lver;

		    /*~ Correction so that the forces don't become unreasonable due to contact
		      with 45 degree planes extended into the basket region, i.e., the trueoverlap
		      above could identify a spurious contact point which is within the cylindrical
		      basket and should not be present in reality.*/
		    double trueoverlaplimitver = radialdist + sqrt(fabs(radius_[iPart]*radius_[iPart] - alpha*alpha)) - rcylinder;

		    if (trueoverlap > trueoverlaplimitver) trueoverlap = trueoverlaplimitver;
		  
		    if (trueoverlap > 0.0) { //~ There is an overlap so a force must be calculated
		      trueoverlapflag = 1;
		      /*~ The current forces on the particle, prior to adding the 
			sphere-cylinder force, are stored in 'oldvforce'. Then the 
			particle force is updated for the sphere-solid cylinder contact. 
			This force and oldvforce are subtracted to find the contact force
			due to an overlap of (-deltan). Knowing the real overlap of 
			'trueoverlap', the contact force magnitude can be scaled to
			be commensurate with 'trueoverlap'. Finally, this contact force 
			magnitude is applied in the appropriate direction (45 degrees
			to the vertical and towards the y-axis).*/
		      double oldvforce[3] = {f_[iPart][0], f_[iPart][1], f_[iPart][2]};
              
		      /*~ The section enclosed between ~~~ has been copy-pasted from 
			below and simplified to include a sphere-cylinder contact force
			in the particle's net force. The only modification is to
			redefine deltan so that things work sensibly when radialdist >
			cylinder, which does happen in the transition regions.*/
		      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		      if (magdiffforce < 0.0) {
			if (radialdist > rcylinder) deltan = rcylinder - radialdist - radius_[iPart];
			if (!cdata.is_non_spherical) cdata.deltan = -deltan;
			cdata.delta[0] = -delta[0];
			cdata.delta[1] = -delta[1];
			cdata.delta[2] = -delta[2];
			cdata.radi = radius_ ? radius_[iPart] : r0_;
			cdata.r = cdata.radi - cdata.deltan; // sign of corrected, because negative value is passed
			cdata.rsq = cdata.r*cdata.r;
			cdata.meff = rmass_ ? rmass_[iPart] : atom->mass[atom->type[iPart]];
			cdata.area_ratio = 1.;
		      }
		    
		      if (impl) impl->compute_force(this, cdata, intersectflag, v_wall, i_forces, j_forces);
		      else {
			cdata.r =  r0_ - cdata.deltan;
			compute_force(cdata, v_wall); // LEGACY CODE (SPH)
		      }
		      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              
		      MathExtra::sub3(f_[iPart],oldvforce,diffforce); //~ New force minus old force vectors
		      magdiffforce = MathExtra::len3(diffforce);
              
		      //~ Assuming a Hertzian model in which the normal force is proportional to overlap to the power of 1.5
		      desiredmagforce = magdiffforce*pow(-trueoverlap/deltan,1.5); //~ The force magnitude to apply
		    
		      //~ Find the direction in which this is applied
		      double forcevdirection[3] = {-x_[iPart][0], radialdist, -x_[iPart][2]};
		      if (fregion[0] == 2) forcevdirection[1] *= -1.0; //~ Flip the sign for a downwards component of velocity
		
		      MathExtra::norm3(forcevdirection); //~ Normalise this vector
		      MathExtra::scale3(desiredmagforce,forcevdirection); //~ And scale it by the required magnitude
		
		      //~ Now apply this force to the particle, overwriting the force update done by the compute_force function
		      MathExtra::add3(oldvforce, forcevdirection, f_[iPart]);
		    } //~ Closes off 'if (trueoverlap > 0.0)'

		    fregion[0] -= 3; //~ Reduce fregion[0] from 5 to 2 to run another iteration of the while loop, or ensure exiting if fregion[0] == 2 or 3
		  } //~ Closes off 'while (fregion[0] >= 2)'
		} //~ Closes off 'if (fregion[0] >= 2)'
	  
		if (fregion[1] >= 2) {
		  /*~ This is easier in one way: the modified force acts solely in the x-z plane,
		    i.e., has no vertical component.
                  
		    There is a convenient equation for the shortest (perpendicular) distance between
		    a point and a line defined by a point and an angle 
		    (https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line).
		    This is the basis for finding Lhor.
                   
		    The fractional angular distances are all for the particle within the planar 
		    circumferential cell in its reference configuration, i.e., starting from the 
		    positive x-axis. However, the particle's location is taken to be at 'unrotatedangle' 
		    so that the angle of the 45 degree wall is correct.*/
                   
		  /*~ Firstly find the signed angular difference between fh (where the particle's centre 
		    is in the reference configuration) and the point at which the 45 degree wall touches 
		    the solid cylinder. Two transition regions -> two cases*/
		  double angulardiff, fdh;

		  while (fregion[1] >= 2) {
		    fregion[1] == 2 ? fdh = ftransition_hor - fh : fdh = fsolid_hor - fh; //~ Start with transition region 2 if fregion[1] == 5
		    angulardiff = fdh*Lh;

		    //~ Angle to known point on 45 degree wall (line) and the angle of the line itself, theta
		    double angletoknownpointonline, theta, Lhor;
		    angletoknownpointonline = unrotatedangle + angulardiff;
		    double pointonline[2] = {rcylinder*cos(angletoknownpointonline), rcylinder*sin(angletoknownpointonline)};
		    fregion[1] == 2 ? theta = angletoknownpointonline - M_PI/4.0 : theta = angletoknownpointonline + M_PI/4.0; //~ +/- 45 degrees
		    Lhor = fabs(cos(theta)*(pointonline[1]-x_[iPart][2]) - sin(theta)*(pointonline[0]-x_[iPart][0]));
		    trueoverlap = radius_[iPart] - Lhor;

		    /*~ Again, correction so that the forces don't become unreasonable due to contact
		      with 45 degree planes extended into the basket region.*/
		    double equivalpha = 2.0*rcylinder*sin(0.5*fabs(angulardiff)); //~ Chord length which is equivalent to alpha in the vertical direction
		    double trueoverlaplimithor = radialdist + sqrt(fabs(radius_[iPart]*radius_[iPart] - equivalpha*equivalpha)) - rcylinder;
		    if (trueoverlap > trueoverlaplimithor) trueoverlap = trueoverlaplimithor;
		  
		    if (trueoverlap > 0.0) {
		      trueoverlapflag = 1;
		      double oldhforce[3] = {f_[iPart][0], f_[iPart][1], f_[iPart][2]}; //~ Always store this anew
		  
		      if (magdiffforce < 0.0) {
			/*~ The section enclosed between ~~~ has been copy-pasted from 
			  below and simplified to include a sphere-cylinder contact force
			  in the particle's net force. The only modification is to
			  redefine deltan so that things work sensibly when radialdist >
			  cylinder, which does happen in the transition regions.*/
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			if (radialdist > rcylinder) deltan = rcylinder - radialdist - radius_[iPart];
			if (!cdata.is_non_spherical) cdata.deltan = -deltan;
			cdata.delta[0] = -delta[0];
			cdata.delta[1] = -delta[1];
			cdata.delta[2] = -delta[2];
			cdata.radi = radius_ ? radius_[iPart] : r0_;
			cdata.r = cdata.radi - cdata.deltan; // sign of corrected, because negative value is passed
			cdata.rsq = cdata.r*cdata.r;
			cdata.meff = rmass_ ? rmass_[iPart] : atom->mass[atom->type[iPart]];
			cdata.area_ratio = 1.;
		      }
		  
		      if (impl) impl->compute_force(this, cdata, intersectflag, v_wall, i_forces, j_forces);
		      else {
			cdata.r =  r0_ - cdata.deltan;
			compute_force(cdata, v_wall); // LEGACY CODE (SPH)
		      }
		      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  
		      MathExtra::sub3(f_[iPart],oldhforce,diffforce); //~ New force minus old force vectors
		      magdiffforce = MathExtra::len3(diffforce);
		      desiredmagforce = magdiffforce*pow(-trueoverlap/deltan,1.5);
		
		      /*~ Follow a similar procedure to above to find the direction
			of the contact force applied to the particle*/
		      double forcehdirection[3] = {sin(theta), 0.0, cos(theta)}; //~ Already a unit vector
		      if (fregion[1] == 2) forcehdirection[2] *= -1.0; //~ Flip the signs as needed
		      else forcehdirection[0] *= -1.0;
		
		      MathExtra::scale3(desiredmagforce,forcehdirection); //~ Scale it by the required magnitude and apply as above
		      MathExtra::add3(oldhforce, forcehdirection, f_[iPart]);
		    } //~ Closes off 'if (trueoverlap > 0.0)'

		    fregion[1] -= 3; //~ Reduce fregion[1] from 5 to 2 to run another iteration of the while loop, or ensure exiting if fregion[1] == 2 or 3
		  } //~ Closes off 'while (fregion[0] >= 2)'
		} //~ Closes off 'if (fregion[1] >= 2)'

		if (!trueoverlapflag && c_history) vectorZeroizeN(c_history[iPart],dnum_); //~ Delete the shear history if no contact
		continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
	      } //~ Closes off 'if (fregion[0] >= 2 || fregion[1] >= 2)'
	    } //~ Closes off 'if (fregion[0] != 1 && fregion[1] != 1)'
          
	    /*~ Do nothing if the particle contacts a solid region; the calculations
	      continue as for a typical solid cylinder*/
	  } //~ Closes off 'if (nPrimitiveArgs == 11)'
	
	  if (nPrimitiveArgs == 12) {

	    /*~ Here are the function arguments (12th argument is a dummy one used for differentiating
	      this primitive from the finite porous cylinder with 11 arguments):
		
	      param[0] = position of the y-plane (The y coordinate, meaning the rectangle will be on
	      the x, z dimensions)
	      param[1] = first coordinate of the left-bottom end of the rectangle (x)
	      param[2] = other coordinate of the left-bottom end of the rectangle (z)
	      param[3] = np_ver (integer number of pores in the vertical (z) direction; must be >= 1)
	      param[4] = np_hor (integer number of pores in a horizontal (x); must be >= 1)
	      param[5] = A_ver (height of the rectangle)
	      param[6] = A_hor (length of the rectangle)
	      param[7] = fvoid_ver (fractional distance in vertical repeating cell which is void)
	      param[8] = fvoid_hor (fractional distance in horizontal repeating cell which is void)
	      param[9] = ftransition_ver (fractional distance in vertical repeating cell which is void 
	      + transition)
	      param[10] = ftransition_hor (fractional distance in horizontal repeating cell which is 
	      void + transition)
	      param[11] = shape (dummy double variable, give '0')*/

	    /*~ Some definitions needed for the code: */
	    int deleteparticle;
	    const double sin_45 = sin(M_PI/4);
	    const double cos_45 = cos(M_PI/4);
	 
	    /*~ If the particle's x-z coordinates are outside the rectangle dimensions A_ver, A_hor,
	      it does not contact the plane of the finite length. In this case, zero shear forces are
	      applied and code can continue to next particle in the neighbour list.*/

	    if (x_[iPart][0] > C_x + A_hor || x_[iPart][0] < C_x || x_[iPart][2] > C_z + A_ver || x_[iPart][2] < C_z) {
	      if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
	      continue; //~ Shortcut the calculations by skipping to the next particle in the neighbour list
	    }

	    /*~ The spherical particle would be in contact with a ractangular region with repeating
	      pores. Now determine, based on its location in the vertical (z) and in the horizontal
	      (x) directions, whether the sphere is contacting a solid region of the mesh or a void, 
	      given that it is in a contact distance in the y direction. In here, we assume that the 
	      plane is a y plane and the particles are falling on to the porous region from the above.*/

	    int trueoverlapflag = 0;
	    double fvtemp, fv, y_distance;

	    y_distance = x_[iPart][1] - C_y; //~ Particle distance from the plane origin (bottom left corner) in the  y direction 

	    fvtemp = (x_[iPart][2] - C_z)/Lv; //~ x_[iPart][2] is the z-coordinate of the particle
	    fv = fmod(fvtemp, 1.0); //~ Fractional height of the particle (0.0-1.0) in the vertical repeating cell

	    /*~ Determine if this fv corresponds to void (0), solid (1) or 
	      transition (2 or 3) regions. 2 is used for the first transition
	      from void to solid with increasing fv; 3 is used for the
	      second transition from solid to void with increasing fv*/

	    int fregion[2]; //~ [fv, fh]
	
	    if (radius_[iPart] >= criticalradius) {
	      fregion[0] = 1; //~ Set as a solid contact to bypass if the particle is too large to egress
	    } else {
	      fregion[0] = region_finder(fvoid_ver, ftransition_ver, fsolid_ver, fv);
	    }
	
	    double fhtemp, fh;
	
	    if (fregion[0] != 1) {
	      //~ Only need to do these calculations if this location is not definitely solid
	    
	      fhtemp = (x_[iPart][0] - C_x)/Lh; //~ x_[iPart][0] is the x-coordinate of the particle
	      fh = fmod(fhtemp, 1.0); //~ Fractional height of the particle (0.0-1.0) in the horizontal repeating cell
	      fregion[1] = region_finder(fvoid_hor, ftransition_hor, fsolid_hor, fh);
	    }

	    /*~ Deal with the case of a particle contacting a void region in BOTH 
	      dimensions: zero shear forces and continue*/

	    if (fregion[0] == 0 && fregion[1] == 0) {
	      /*~ Here is the modification to the existing particle egress concept based solely
		on the location of the particle's centrepoint. Now we also take into consideration
		the location of the lowermost, uppermost, leftmost and rightmost (in a 2D
		projection sense) points on the spherical particle.

		Lines from the rectangle origin through all four of these additional
		points must also lie within the void region in both dimensions for egress to be
		permitted. For any point that does not satisfy this, we apply a force to the
		particle which is for a contact with a 45-degree transition plane; thus free
		passage of the particle through the void is hindered.*/

	      int fregion_fourpoints[4]; //~ [lowermost (vertical, z), uppermost (vertical, z), negative direction (horizontal, x), positive direction (horizontal, x)]
	      double fvtemp_lower, fv_lower;
	      double fvtemp_upper, fv_upper;
	      double fhtemp_negative, fh_negative;
	      double fhtemp_positive, fh_positive;
	    
	      fvtemp_lower = (x_[iPart][2] - radius_[iPart] - C_z)/Lv;
	      fv_lower = fmod(fvtemp_lower, 1.0); //~ Not a problem if fv_lower is slightly negative because it will still be mapped to void by region_finder
	      fregion_fourpoints[0] = region_finder(fvoid_ver, ftransition_ver, fsolid_ver, fv_lower);
 
	      fvtemp_upper = (x_[iPart][2] + radius_[iPart] - C_z)/Lv;
	      fv_upper = fmod(fvtemp_upper, 1.0);
	      fregion_fourpoints[1] = region_finder(fvoid_ver, ftransition_ver, fsolid_ver, fv_upper);
		
	      fhtemp_negative = (x_[iPart][0] - radius_[iPart] - C_x)/Lh;
	      fh_negative = fmod(fhtemp_negative, 1.0);
	      fregion_fourpoints[2] = region_finder(fvoid_hor, ftransition_hor, fsolid_hor, fh_negative);
	    
	      fhtemp_positive = (x_[iPart][0] + radius_[iPart] - C_x)/Lh;
	      fh_positive = fmod(fhtemp_positive, 1.0);
	      fregion_fourpoints[3] = region_finder(fvoid_hor, ftransition_hor, fsolid_hor, fh_positive);
	    
	      if (fregion_fourpoints[0] == 0 && fregion_fourpoints[1] == 0 && fregion_fourpoints[2] == 0 && fregion_fourpoints[3] == 0) {//~ If all are zero
	      
		if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
		deleteparticle = outside_rectangle_check(C_y, radius_[iPart], iPart);
		continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
	    
	      } else {
	    
		/*~ Else assume contact with a transition region to prevent free passage of the
		  particle out of the basket. The code below sets fregion[*] to denote contact
		  with one or more transition regions. Note the way it is done means that fregion[*]
		  can be 2, 3 or 5 with the '5' value indicating contact with both transition regions
		  simultaneously, i.e., a particle that is larger than a pore that can never egress*/
	      
		if (fregion_fourpoints[0] > 0) fregion[0] += 3;
		if (fregion_fourpoints[1] > 0) fregion[0] += 2;
		if (fregion_fourpoints[2] > 0) fregion[1] += 3;
		if (fregion_fourpoints[3] > 0) fregion[1] += 2;
	      }
	    }

  
	    /*~ Pass on to a function if the particle contacts a transition 
	      region. Bypass if the particle contacts a solid region ('do 
	      nothing' in table of actions). The vertical and horizontal 
	      directions are treated as independent for the purpose 
	      of calculating the modified forces.
     
	      For convenience, it has been decided to take the physical form of 
	      these transitions as walls with 45 degree orientations.*/

	    ForceData i_forces;
	    ForceData j_forces;

	    if (fregion[0] != 1 && fregion[1] != 1) {
	      if (fregion[0] >= 2 || fregion[1] >= 2) {
		
		//~ Check for the particle being outside of the basket in this case too
		deleteparticle = outside_rectangle_check(C_y, radius_[iPart], iPart);
		if (deleteparticle == 1) continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
			
		double diffforcever[3]; //~ The contact force vector due to an overlap of -deltan
		double magdiffforcever = -1.0; //~ The magnitude of this contact force vector. Initialise at an unreasonable value
		double trueoverlapver; //~ The overlap between the particle and the 45 degree plane
		double desiredmagforcever; //~ The force magnitude to apply
			
		double diffforcehor[3]; //~ The contact force vector due to an overlap of -deltan
		double magdiffforcehor = -1.0; //~ The magnitude of this contact force vector. Initialise at an unreasonable value
		double trueoverlaphor; //~ The overlap between the particle and the 45 degree plane
		double desiredmagforcehor; //~ The force magnitude to apply
        
		if (fregion[0] >= 2) {
		  double numeratorinalphaver, transitionfracver;

		  while (fregion[0] >= 2) {
		    fregion[0] == 2 ? numeratorinalphaver = fmod(ftransition_ver - fv, 1.0) : numeratorinalphaver = fmod(fv - fsolid_ver, 1.0); //~ Start with transition region 2 if fregion[0] == 5
		    if (numeratorinalphaver < 0.0) numeratorinalphaver += 1.0; //~ Ensure always in the range [0.0, 1.0]
		    transitionfracver = 1.0 - fsolid_ver;

		    double alphaver = fdv*numeratorinalphaver/transitionfracver; //~ Distance units due to multiplication by fdv from constructor
		    double betaver = alphaver + y_distance;
		    double Lver = betaver*sin_45;
		    trueoverlapver = radius_[iPart] - Lver;
				
		    /*~ Correction so that the forces don't become unreasonable due to contact
		      with 45 degree planes extended into the basket region, i.e., the trueoverlap
		      above could identify a spurious contact point that should not be present in reality.*/
		  			
		    double trueoverlaplimitver = - x_[iPart][1] + sqrt(fabs(radius_[iPart]*radius_[iPart] - alphaver*alphaver)) + C_y;
				
		    if (trueoverlapver > trueoverlaplimitver) trueoverlapver = trueoverlaplimitver;
		  
		    if (trueoverlapver > 0.0) { //~ There is an overlap so a force must be calculated
		    			
		      trueoverlapflag = 1;
		    			
		      /*~ The current forces on the particle, prior to adding the 
		      	new force, are stored in 'oldforce'. Then the particle force is 
			updated for the sphere-solid plane contact. This force and oldforce 
			are subtracted to find the contact force due to an overlap of (-deltan). 
			Knowing the real overlap of 'trueoverlap', the contact force magnitude 
			can be scaled to be commensurate with 'trueoverlap'. Finally, this contact 
			force magnitude is applied in the appropriate direction. */
		    			
		      double oldforcever[3] = {f_[iPart][0], f_[iPart][1], f_[iPart][2]};

		      if (magdiffforcever < 0.0) {

			if (x_[iPart][1] < C_y) deltan = x_[iPart][1] - C_y - radius_[iPart];
		     			
			if (!cdata.is_non_spherical) cdata.deltan = -deltan;
		      			
			cdata.delta[0] = -delta[0];
			cdata.delta[1] = -delta[1];
			cdata.delta[2] = -delta[2];
			cdata.radi = radius_ ? radius_[iPart] : r0_;
			cdata.r = cdata.radi - cdata.deltan; // sign of corrected, because negative value is passed
			cdata.rsq = cdata.r*cdata.r;
			cdata.meff = rmass_ ? rmass_[iPart] : atom->mass[atom->type[iPart]];
			cdata.area_ratio = 1.;

		      }

		      if (impl) impl->compute_force(this, cdata, intersectflag, v_wall, i_forces, j_forces);
		      else {
		      			
			cdata.r =  r0_ - cdata.deltan;
			compute_force(cdata, v_wall); // LEGACY CODE (SPH)
		    			
		      }

		      MathExtra::sub3(f_[iPart],oldforcever,diffforcever); //~ New force minus old force vectors
		      magdiffforcever = MathExtra::len3(diffforcever);
              
		      //~ Assuming a Hertzian model in which the normal force is proportional to overlap to the power of 1.5
		      desiredmagforcever = magdiffforcever*pow(-trueoverlapver/deltan,1.5); //~ The force magnitude to apply

		      //~ Find the direction in which this is applied
		      double forcedirectionver[3] = {0, sin_45, cos_45};
		      if (fregion[0] == 2) forcedirectionver[2] *= -1.0; //~ Flip the sign for a downwards component of velocity
		
		      MathExtra::norm3(forcedirectionver); //~ Normalise this vector
		      MathExtra::scale3(desiredmagforcever,forcedirectionver); //~ And scale it by the required magnitude
		
		      //~ Now apply this force to the particle, overwriting the force update done by the compute_force function
		      MathExtra::add3(oldforcever, forcedirectionver, f_[iPart]);

		      //fprintf(screen,"oldforcever[0] %1.9f oldforcever[1] %1.9f oldforcever[2] %1.9f newforce[0] %1.9f newforce[1] %1.9f newforce[2] %1.9f on timestep " BIGINT_FORMAT "\n",oldforcever[0],oldforcever[1],oldforcever[2],f_[iPart][0],f_[iPart][1],f_[iPart][2],update->ntimestep);
	
		    		 		  				
		    } //~ Closes off 'if (trueoverlapver > 0.0)'
         
		    fregion[0] -= 3; //~ Reduce fregion[0] from 5 to 2 to run another iteration of the while loop, or ensure exiting if fregion[0] == 2 or 3
				
		  } //~ Closes off 'while (fregion[0] >= 2)'
	      	
		} //~ Closes off 'if (fregion[0] >= 2)'
		
		if (fregion[1] >= 2) {

		  double numeratorinalphahor, transitionfrachor;

		  while (fregion[1] >= 2) {

		    fregion[1] == 2 ? numeratorinalphahor = fmod(ftransition_hor - fh, 1.0) : numeratorinalphahor = fmod(fh - fsolid_hor, 1.0); //~ Start with transition region 2 if fregion[0] == 5
		  		
		    if (numeratorinalphahor < 0.0) numeratorinalphahor += 1.0; //~ Ensure always in the range [0.0, 1.0]
		    transitionfrachor = 1.0 - fsolid_hor;

		    double alphahor = fdh*numeratorinalphahor/transitionfrachor; //~ Distance units due to multiplication by fdh from constructor
		    double betahor = alphahor+y_distance;
		    double Lhor = betahor*sin_45;
		    trueoverlaphor = radius_[iPart] - Lhor;

		    /*~ Correction so that the forces don't become unreasonable due to contact
		      with 45 degree planes extended into the basket region, i.e., the trueoverlap
		      above could identify a spurious contact point that should not be present in reality.*/
		  			
		    double trueoverlaplimithor = - x_[iPart][1] + sqrt(fabs(radius_[iPart]*radius_[iPart] - alphahor*alphahor)) + C_y;

		    if (trueoverlaphor > trueoverlaplimithor) trueoverlaphor = trueoverlaplimithor;
		  
		    if (trueoverlaphor > 0.0) { //~ There is an overlap so a force must be calculated
		    			
		      trueoverlapflag = 1;
		    			
		      /*~ The current forces on the particle are stored in 'oldforce'. Then the 
			particle force is updated for the sphere-solid plane contact. 
			This force and oldforce are subtracted to find the contact force
			due to an overlap of (-deltan). Knowing the real overlap of 
			'trueoverlap', the contact force magnitude can be scaled to
			be commensurate with 'trueoverlap'. Finally, this contact force 
			magnitude is applied in the appropriate direction.*/
		    			
		      double oldforcehor[3] = {f_[iPart][0], f_[iPart][1], f_[iPart][2]};

		      if (magdiffforcehor < 0.0) {

			if (x_[iPart][1] < C_y) deltan = x_[iPart][1] - C_y - radius_[iPart];
		     			
			if (!cdata.is_non_spherical) cdata.deltan = -deltan;
		      			
			cdata.delta[0] = -delta[0];
		      	cdata.delta[1] = -delta[1];
		      	cdata.delta[2] = -delta[2];
		      	cdata.radi = radius_ ? radius_[iPart] : r0_;
		      	cdata.r = cdata.radi - cdata.deltan; // sign of corrected, because negative value is passed
		      	cdata.rsq = cdata.r*cdata.r;
		      	cdata.meff = rmass_ ? rmass_[iPart] : atom->mass[atom->type[iPart]];
		      	cdata.area_ratio = 1.;

		      }

		      if (impl) impl->compute_force(this, cdata, intersectflag, v_wall, i_forces, j_forces);
		      else {
			cdata.r =  r0_ - cdata.deltan;
			compute_force(cdata, v_wall); // LEGACY CODE (SPH)
		      }

		      MathExtra::sub3(f_[iPart],oldforcehor,diffforcehor); //~ New force minus old force vectors
		      magdiffforcehor = MathExtra::len3(diffforcehor);
              
		      //~ Assuming a Hertzian model in which the normal force is proportional to overlap to the power of 1.5
		      desiredmagforcehor = magdiffforcehor*pow(-trueoverlaphor/deltan,1.5); //~ The force magnitude to apply
		    
		      //~ Find the direction in which this is applied
		      double forcedirectionhor[3] = {cos_45, sin_45, 0};
		      if (fregion[1] == 2) forcedirectionhor[0] *= -1.0; //~ Flip the sign for a downwards component of velocity
		
		      MathExtra::norm3(forcedirectionhor); //~ Normalise this vector
		      MathExtra::scale3(desiredmagforcehor,forcedirectionhor); //~ And scale it by the required magnitude
		
		      //~ Now apply this force to the particle, overwriting the force update done by the compute_force function
		      MathExtra::add3(oldforcehor, forcedirectionhor, f_[iPart]);
		  				
		    } //~ Closes off 'if (trueoverlaphor > 0.0)'

		    fregion[1] -= 3; //~ Reduce fregion[1] from 5 to 2 to run another iteration of the while loop, or ensure exiting if fregion[1] == 2 or 3
			
		  } //~ Closes off 'while (fregion[0] >= 2)'
	    
		} //~ Closes off 'if (fregion[1] >= 2)'

		if (!trueoverlapflag && c_history) vectorZeroizeN(c_history[iPart],dnum_); //~ Delete the shear history if no contact
	    
		continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
	    
	      } //~ Closes off 'if (fregion[0] >= 2 || fregion[1] >= 2)'

	    } //~ Closes off 'if (fregion[0] != 1 && fregion[1] != 1)'
          
	    /*~ Do nothing if the particle contacts a solid region; the calculations
	      continue as for a typical solid plane.*/

	  } //~ Closes off 'if (nPrimitiveArgs == 12)'

	  if (nPrimitiveArgs == 5) {//~ ycylinder_finite
	    /*~ If the particle's y coordinate is above yhigh or below ylow, it does
	      not contact the cylindrical inlet pipe of finite length. In this case,
	      zero shear forces and continue*/
	    if (x_[iPart][1] > yhigh || x_[iPart][1] < ylow) {//~ x_[iPart][1] is the y-coordinate of the particle
	      if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
	      continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
	    }
	  } //~ Closes off 'if (nPrimitiveArgs == 5)'

	  if (nPrimitiveArgs == 2) {//~ yplane_circle_finite
	    double radialdist = sqrt(x_[iPart][0]*x_[iPart][0] + x_[iPart][2]*x_[iPart][2]); //~ Radial distance of particle from basket centreline

	    /*~ If the particle's radialdist is less than the inlet pipe's radius, it 
	      can pass through the wall. In this case, zero shear forces and continue*/
	    if (direc == -1) {
	      if (radialdist < rcircle) {
		if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
		continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
	      }
	    } else {
	      if (radialdist > rcircle) {
		if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
		continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
	      }
	    }
	  } //~ Closes off 'if (nPrimitiveArgs == 2)'
	  if (nPrimitiveArgs == 4) {
	    double radialdist = sqrt(x_[iPart][0]*x_[iPart][0] + x_[iPart][2]*x_[iPart][2]); //~ Radial distance of particle from basket centreline
	    if (direc == -1) {
	      if (radialdist > rcircle && radialdist < rocircle) {
		if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
		continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list
	      } else if (direc == 1) {
		if (radialdist < rcircle || radialdist > rocircle) {
		  if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
		  continue; //~ Shortcut the calculations by skipping to the next particle in neighbour list	  
		}
	      }
	    } 
	  }
	
#ifdef SUPERQUADRIC_ACTIVE_FLAG
	  if(atom->superquadric_flag) {
	    double sphere_contact_point[3];
	    vectorAdd3D(x_[iPart], delta, sphere_contact_point);
	    double closestPoint[3], closestPointProjection[3], point_of_lowest_potential[3];

	    Superquadric particle(x_[iPart], quat_[iPart], shape_[iPart], blockiness_[iPart]);
	    intersectflag = particle.plane_intersection(delta, sphere_contact_point, closestPoint, point_of_lowest_potential);
	    deltan = -MathExtraLiggghtsNonspherical::point_wall_projection(delta, sphere_contact_point, closestPoint, closestPointProjection);

	    vectorCopy3D(closestPoint, cdata.contact_point);
	    cdata.is_non_spherical = true; // by default it is false
	  }
#endif

	  //~ Disabled this for a particle-porous cylinder contact as it doesn't seem to work
	  //~ Instead use the v_wall values calculated earlier
	  if(shear_ && shearAxis_ >= 0 && nPrimitiveArgs != 11)
	    {
	      double rdist[3];
	      primitiveWall->calcRadialDistance(x_[iPart],rdist);
	      vectorCross3D(shearAxisVec_,rdist,v_wall);
	    }

#ifdef SUPERQUADRIC_ACTIVE_FLAG
	  if(atom->superquadric_flag && deltan > 0.0)
	    {
	      if(c_history)
		vectorZeroizeN(c_history[iPart],dnum_);
	      break;
	    }

	  if(!cdata.is_non_spherical || atom->superquadric_flag)
#endif
	    cdata.deltan = -deltan;

	  cdata.delta[0] = -delta[0];
	  cdata.delta[1] = -delta[1];
	  cdata.delta[2] = -delta[2];
	  post_force_eval_contact(cdata,intersectflag,v_wall);
	}
      else
	{
	  if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
	}
    }
  }
}

void FixWallGran::compute_force(CollisionData &, double *)
{
}

/* ---------------------------------------------------------------------- */

int FixWallGran::outside_basket_check(double radialdist, double rcylinder, double radius, int iPart)
{
  /*~ Function to determine if the particle is far enough outside of the
    basket that it will never reenter, in which case it can be removed from
    the simulation. An interesting quirk is the triggering distance cannot be
    the particle's radius since, in that case, there would be no contact between
    the particle and cylinder so this function would never be invoked. Set it to
    95% of the particle's radius instead.*/
  double rdist = rcylinder + 0.95*radius;
  if (radialdist > rdist) {
    /*~ Particle can be removed from the simulation. Rather than deleting
      it explicitly which is computationally expensive, move the particle
      far from the basket so that it will be lost due to the boundary
      conditions adopted*/
    double multfactor = 2.0; //~ Factor to define how far to move the particle. Typically will ensure immediate deletion.
    double **x_ = atom->x;
    x_[iPart][0] *= multfactor;
    x_[iPart][2] *= multfactor;
    return 1;
  } else return 0;
}

int FixWallGran::outside_rectangle_check(double C_y, double radius, int iPart)
{
  /*~ Function to determine if the particle is far enough outside of the
    plane that it will never reenter, in which case it can be removed from
    the simulation.*/ 
  if (C_y > x_[iPart][1]) {
    /*~ Particle can be removed from the simulation. Rather than deleting
      it explicitly which is computationally expensive, move the particle
      far from the plane so that it will be lost due to the boundary
      conditions adopted*/
    double **x_ = atom->x;
    x_[iPart][1] = x_[iPart][1] - 40*radius;
    return 1;
  } else return 0;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::region_finder(double fvoid, double ftransition, double fsolid, double f)
{
  if (f <= fvoid) return 0;
  else if (f >= ftransition && f <= fsolid) return 1;
  else if (f < ftransition) return 2;
  else return 3;
}

/* ----------------------------------------------------------------------
   actually calculate force, called for both mesh and primitive
   ------------------------------------------------------------------------- */

inline void FixWallGran::post_force_eval_contact(CollisionData & cdata, bool intersectflag, double * v_wall, int iMesh, FixMeshSurface *fix_mesh, TriMesh *mesh, int iTri)
{
  const int iPart = cdata.i;

  // deltan > 0 in compute_force
  // but negative in distance algorithm
  cdata.radi = radius_ ? radius_[iPart] : r0_;
  cdata.r = cdata.radi - cdata.deltan; // sign of corrected, because negative value is passed
  cdata.rsq = cdata.r*cdata.r;
  cdata.meff = rmass_ ? rmass_[iPart] : atom->mass[atom->type[iPart]];
  cdata.area_ratio = 1.;

  ForceData i_forces;
  ForceData j_forces;

  // add to cwl
  if(cwl_ && addflag_)
    {
      double contactPoint[3];
      vectorSubtract3D(x_[cdata.i],cdata.delta,contactPoint);
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(atom->superquadric_flag) {
	vectorCopy3D(cdata.contact_point, contactPoint);
      }
#endif
      cwl_->add_wall_1(iMesh,mesh->id(iTri),iPart,contactPoint,v_wall);
    }
  
  if(impl) {
    impl->compute_force(this, cdata, intersectflag, v_wall, i_forces, j_forces);
  } else {
    double force_old[3]={};

    // if force should be stored - remember old force
    if(store_force_ || stress_flag_)
      vectorCopy3D(f_[iPart],force_old);

    compute_force(cdata, v_wall); // LEGACY CODE (SPH)

    if(store_force_ || stress_flag_)
      {
	vectorSubtract3D(f_[iPart], force_old, i_forces.delta_F);
      }
  }
  
  // if force should be stored or evaluated
  if(cdata.has_force_update && (store_force_ || stress_flag_))
    {
      if(store_force_)
        vectorAdd3D (wallforce_[iPart], i_forces.delta_F, wallforce_[iPart]);

      if(stress_flag_ && fix_mesh->trackStress())
	{
	  double delta[3];
	  delta[0] = -cdata.delta[0];
	  delta[1] = -cdata.delta[1];
	  delta[2] = -cdata.delta[2];
	  static_cast<FixMeshSurfaceStress*>(fix_mesh)->add_particle_contribution
	    (
	     iPart,i_forces.delta_F,delta,iTri,v_wall
	     );
	}
    }

  // add heat flux
  if(heattransfer_flag_)
    addHeatFlux(mesh,iPart,cdata.deltan,1.);
}

/* ---------------------------------------------------------------------- */

int FixWallGran::is_moving()
{
  if(is_mesh_wall())
    {
      for(int i = 0; i < n_FixMesh_; i++) {
	if(FixMesh_list_[i]->mesh()->isMoving())
	  return 1;
      }
      return 0;
    }
  return shear_;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local()
{
  if (!is_mesh_wall() || dnum() == 0) return 0;

  int ncontacts = 0;
  for(int i = 0; i < n_FixMesh_; i++)
    ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts();

  return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all()
{
  int ncontacts = n_contacts_local();
  MPI_Sum_Scalar(ncontacts,world);
  return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local(int contact_groupbit)
{
  if (!is_mesh_wall() || dnum() == 0) return 0;

  int ncontacts = 0;
  for(int i = 0; i < n_FixMesh_; i++)
    ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts(contact_groupbit);

  return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all(int contact_groupbit)
{
  int ncontacts = n_contacts_local(contact_groupbit);
  MPI_Sum_Scalar(ncontacts,world);
  return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_ghosts_local()
{
  if (!is_mesh_wall()) return 0;

  int nghosts = 0;
  for(int iMesh = 0; iMesh < n_FixMesh_; ++iMesh) {
    TriMesh *mesh = FixMesh_list_[iMesh]->triMesh();
    nghosts += mesh->sizeGhost();
  }

  return nghosts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_ghosts_all()
{
    int nghosts = n_ghosts_local();
    MPI_Sum_Scalar(nghosts,world);
    return nghosts;
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void FixWallGran::register_compute_wall_local(ComputePairGranLocal *ptr,int &dnum_compute)
{
   if(cwl_ != NULL)
     error->fix_error(FLERR,this,"Fix wall/gran allows only one compute of type wall/gran/local");
   cwl_ = ptr;
   dnum_compute = dnum_; //history values
}

void FixWallGran::unregister_compute_wall_local(ComputePairGranLocal *ptr)
{
   if(cwl_ != ptr)
     error->fix_error(FLERR,this,"Illegal situation in FixWallGran::unregister_compute_wall_local");
   cwl_ = NULL;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::init_heattransfer()
{
    fppa_T = NULL;
    fppa_hf = NULL;
    deltan_ratio = NULL;

    // decide if heat transfer is to be calculated

    if (!is_mesh_wall() && Temp_wall < 0.) return;
    else if (is_mesh_wall())
    {
        int heatflag = 0;
        for(int imesh = 0; imesh < n_meshes(); imesh++)
        {
            heatflag = heatflag || mesh_list()[imesh]->mesh()->prop().getGlobalProperty<ScalarContainer<double> >("Temp") != NULL;
        }

        if(!heatflag) return;
    }

    // heat transfer is to be calculated - continue with initializations

    // set flag so addHeatFlux function is called
    heattransfer_flag_ = true;

    // if(screen && comm->me == 0) fprintf(screen,"Initializing wall/gran heat transfer model\n");
    fppa_T = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",1,0,style));
    fppa_hf = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",1,0,style));

    th_cond = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalConductivity","property/global","peratomtype",0,0,style))->get_values();

    // if youngsModulusOriginal defined, get deltan_ratio
    Fix* ymo_fix = modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",0,0,style,false);
    // deltan_ratio is defined by heat transfer fix, see if there is one
    int n_htf = modify->n_fixes_style("heat/gran/conduction");
    if (n_htf < 1) n_htf = modify->n_fixes_style_strict("heat/gran");

    // get deltan_ratio set by the heat transfer fix
    if(ymo_fix && n_htf) deltan_ratio = static_cast<FixPropertyGlobal*>(ymo_fix)->get_array_modified();
}

/* ---------------------------------------------------------------------- */

//NP modified C.K.
void FixWallGran::addHeatFlux(TriMesh *mesh,int ip, double delta_n, double area_ratio)
{
    //r is the distance between the sphere center and wall
    double tcop, tcowall, hc, Acont=0.0, r;
    double reff_wall = atom->radius[ip];
    int itype = atom->type[ip];
    double ri = atom->radius[ip];

    if(mesh)
        Temp_wall = (*mesh->prop().getGlobalProperty< ScalarContainer<double> >("Temp"))(0);

    double *Temp_p = fppa_T->vector_atom;
    double *heatflux = fppa_hf->vector_atom;

    /*NL*/ //if (screen) fprintf(screen,"size %d Temp %f\n",(*mesh->prop().getGlobalProperty< ScalarContainer<double> >("Temp")).size(),Temp_wall);

    if(CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    {
        //NP adjust overlap that may be superficially large due to softening
        if(deltan_ratio)
           delta_n *= deltan_ratio[itype-1][atom_type_wall_-1];

        r = ri - delta_n;

        Acont = (reff_wall*reff_wall-r*r)*M_PI*area_ratio; //contact area sphere-wall
    }
    else if (CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
        Acont = fixed_contact_area_;
    else if (CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    {
        Acont = M_PI*ri*ri;
    }

    tcop = th_cond[itype-1]; //types start at 1, array at 0
    tcowall = th_cond[atom_type_wall_-1];

    if ((fabs(tcop) < SMALL) || (fabs(tcowall) < SMALL)) hc = 0.;
    else hc = 4.*tcop*tcowall/(tcop+tcowall)*sqrt(Acont);

    if(computeflag_)
    {
        heatflux[ip] += (Temp_wall-Temp_p[ip]) * hc;
        Q_add += (Temp_wall-Temp_p[ip]) * hc * update->dt;
    }
    if(cwl_ && addflag_)
        cwl_->add_heat_wall(ip,(Temp_wall-Temp_p[ip]) * hc);
    /*NL*/ //if (screen) fprintf(screen,"adding heat flux of %g to particle %d, wall %f Part %f hc %f\n",(Temp_wall-Temp_p[ip]) * hc,ip,Temp_wall,Temp_p[ip],hc);
    /*NL*/ //if (screen) fprintf(screen,"tcop %f tcowall %f Acont %f reff_wall %f r %f delta_n %f\n",tcop,tcowall,Acont,reff_wall,r,delta_n);
}

int64_t FixWallGran::hashcode() {
  return impl->hashcode();
}
