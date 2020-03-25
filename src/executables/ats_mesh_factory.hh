/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! Simple wrapper that takes a ParameterList and generates all needed meshes.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!
All processes are simulated on a domain, which is discretized through a mesh.

Multiple domains and therefore meshes can be used in a single simulation, and
multiple meshes can be constructed on the fly.

The base mesh represents the primary domain of simulation.  Simple, structured
meshes may be generated on the fly, or complex unstructured meshes are
provided as ``Exodus II`` files.  The base *mesh* list includes either a
GeneratedMesh_,  MeshFromFile_, or LogicalMesh_ spec, as described below.

Additionally, a SurfaceMesh_ may be formed by lifting the surface of a
provided mesh and then flattening that mesh to a 2D surface.  ColumnMeshes_
which split a base mesh into vertical columns of cells for use in 1D models
may also be generated automatically.

Finally, mesh generation is hard and error-prone.  A mesh audit is provided,
which checks for many common geometric and topologic errors in mesh
generation.  This is reasonably fast, even for big meshes, and can be done
through providing a "verify mesh" option.

``[mesh-typed-spec]``

* `"mesh type`" ``[string]`` One of `"generate mesh`", `"read mesh file`",
   `"logical`", `"surface`", `"subgrid`", or `"column`".
* `"_mesh_type_ parameters`" ``[_mesh_type_-spec]`` List of parameters
  associated with the type.
* `"verify mesh`" ``[bool]`` **false** Perform a mesh audit.
* `"deformable mesh`" ``[bool]`` **false** Will this mesh be deformed?


GeneratedMesh
==============

Generated mesh are by definition structured, with uniform dx, dy, and dz.
Such a mesh is specified by a bounding box high and low coordinate, and a list
of number of cells in each direction.

Specified by `"mesh type`" of `"generate mesh`".

``[mesh-type-generate-mesh-spec]``

* `"domain low coordinate`" ``[Array(double)]`` Location of low corner of domain
* `"domain high coordinate`" ``[Array(double)]`` Location of high corner of domain
* `"number of cells`" ``[Array(int)]`` the number of uniform cells in each coordinate direction

Example:

.. code-block:: xml

   <ParameterList name="mesh">
     <ParameterList name="domain">
       <Parameter name="mesh type" type="string" value="generate mesh"/>
       <ParameterList name="generate mesh parameters"/>
         <Parameter name="number of cells" type="Array(int)" value="{{100, 1, 100}}"/>
         <Parameter name="domain low coordinate" type="Array(double)" value="{{0.0, 0.0, 0.0}}" />
         <Parameter name="domain high coordinate" type="Array(double)" value="{{100.0, 1.0, 10.0}}" />
       </ParameterList>
     </ParameterList>   
   </ParameterList>   


MeshFromFile
==============

Meshes can be pre-generated in a multitude of ways, then written to "Exodus
II" file format, and loaded in ATS.

Specified by `"mesh type`" of `"read mesh file`".

``[mesh-type-read-mesh-file-spec]``

* `"file`" ``[string]`` name of pre-generated mesh file. Note that in the case of an
   Exodus II mesh file, the suffix of the serial mesh file must be .exo and 
   the suffix of the parallel mesh file must be .par.
   When running in serial the code will read this the indicated file directly.
   When running in parallel and the suffix is .par, the code will instead read
   the partitioned files, that have been generated with a Nemesis tool and
   named as filename.par.N.r where N is the number of processors and r is the rank.
   When running in parallel and the suffix is .exo, the code will partition automatically
   the serial file.
     
* `"format`" ``[string]`` format of pre-generated mesh file (`"MSTK`" or `"Exodus II`")

Example:

.. code-block:: xml

   <ParameterList name="mesh">
     <ParameterList name="domain">
       <Parameter name="mesh type" type="string" value="read mesh file"/>
       <ParameterList name="read mesh file parameters">
         <Parameter name="file" type="string" value="mesh_filename.exo"/>
         <Parameter name="format" type="string" value="Exodus II"/>
       </ParameterList>   
       <Parameter name="verify mesh" type="bool" value="true" />
     </ParameterList>
   </ParameterList>


LogicalMesh
==============

** Document me! **

Specified by `"mesh type`" of `"logical`".



SurfaceMesh
==============

To lift a surface off of the mesh, a side-set specifying all surface faces
must be given.  These faces are lifted locally, so the partitioning of the
surface cells will be identical to the partitioning of the subsurface faces
that correspond to these cells.  All communication and ghost cells are set up.
The mesh is flattened, so all surface faces must have non-zero area when
projected in the z-direction.  No checks for holes are performed.  Surface
meshes may similarly be audited to make sure they are reasonable for
computation.

Specified by `"mesh type`" of `"surface`".

``[mesh-type-surface-spec]``

* `"surface sideset name`" ``[string]`` The Region_ name containing all surface faces.
* `"surface sideset names`" ``[Array(string)]`` A list of Region_ names containing the surface faces.  Either this or the singular version must be specified.
* `"verify mesh`" ``[bool]`` **false** Verify validity of surface mesh.
* `"export mesh to file`" ``[string]`` Export the lifted surface mesh to this filename.

Example:

.. code-block:: xml

    <ParameterList name="mesh" type="ParameterList">
      <ParameterList name="surface" type="ParameterList">
        <Parameter name="mesh type" type="string" value="surface" />
        <ParameterList name="surface parameters" type="ParameterList">
          <Parameter name="surface sideset name" type="string" value="{surface_region}" />
          <Parameter name="verify mesh" type="bool" value="true" />
          <Parameter name="export mesh to file" type="string" value="surface_mesh.exo" />
        </ParameterList>
      </ParameterList>
      <ParameterList name="domain" type="ParameterList">
        <Parameter name="mesh type" type="string" value="read mesh file" />
        <ParameterList name="read mesh file parameters" type="ParameterList">
          <Parameter name="file" type="string" value="../data/open-book-2D.exo" />
          <Parameter name="format" type="string" value="Exodus II" />
        </ParameterList>
      </ParameterList>
    </ParameterList>


SubgridMeshes
==============

A collection of meshes formed by associating a new mesh with each entity of a
region.  Used for a few cases, including generating a 1D column for each
surface face of a semi-structured subsurface mesh, or for hanging logical
meshes off of each surface cell as a subgrid model, etc.

The subgrid meshes are then named `"MESH_NAME_X"` for each X, which is an
entity local ID, in a provided region of the provided entity type.

**DOCUMENT ME How is the subgrid mesh type specified?  Add examples for Columns and Transport Subgrid model!**

Specified by `"mesh type`" of `"subgrid`".

``[mesh-type-subgrid-spec]``

* `"subgrid region name`" ``[string]`` Region on which each subgrid mesh will be associated.
* `"entity kind`" ``[string]`` One of `"cell`", `"face`", etc.  Entity of the region (usually
   `"cell`") on which each subgrid mesh will be associated.
* `"parent domain`" ``[string]`` **domain** Mesh which includes the above region.
* `"flyweight mesh`" ``[bool]`` **False** NOT SUPPORTED?  Allows a single mesh instead of one per entity.

    
ColumnMeshes
==============

Note these are never? created manually by a user.  Instead use SubgridMeshes,
which generate a ColumnMesh_ spec for every face of the surface mesh.

Specified by `"mesh type`" of `"column`".

``[mesh-type-column-spec]``

* `"parent domain`" ``[string]`` The Mesh_ name of the 3D mesh from which columns are generated.
   Note that the `"build columns from set`" parameter must be set in that mesh.
* `"verify mesh`" ``[bool]`` **false** Verify validity of surface mesh.
* `"deformable mesh`" ``[bool]`` **false**  Used for deformation PKs to allow non-const access.
* `"entity LID`" ``[int]`` Local ID of the surface cell that is the top of the column.

Example:

.. code-block:: xml

    <ParameterList name="mesh" type="ParameterList">
      <ParameterList name="column" type="ParameterList">
        <ParameterList name="column parameters" type="ParameterList">
          <Parameter name="parent domain" type="string" value="domain" />
          <Parameter name="entity LID" type="int" value="0" />
        </ParameterList>
      </ParameterList>
      <ParameterList name="domain" type="ParameterList">
        <Parameter name="mesh type" type="string" value="read mesh file" />
        <ParameterList name="read mesh file parameters" type="ParameterList">
          <Parameter name="file" type="string" value="../data/open-book-2D.exo" />
          <Parameter name="format" type="string" value="Exodus II" />
        </ParameterList>
      </ParameterList>
    </ParameterList>

SperryMesh
==============

A mesh based on the Sperry et al 98 and Christoffersen et al 16 papers.  Hard
coded currently for simplicity, this is ongoing work and will change a lot.

Example:

.. code-block:: xml

    <ParameterList name="mesh" type="ParameterList">
      <ParameterList name="domain" type="ParameterList">
        <Parameter name="mesh type" type="string" value="Sperry 1D column" />
      </ParameterList>
    </ParameterList>
    
*/

#ifndef ATS_MESH_FACTORY_HH_
#define ATS_MESH_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"
#include "State.hh"


namespace ATS {

bool
checkVerifyMesh(Teuchos::ParameterList& mesh_plist,
                Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh);

void
createMesh(Teuchos::ParameterList& plist,
           const Amanzi::Comm_ptr_type& comm,
           const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
           Amanzi::State& s);


void
createMeshes(Teuchos::ParameterList& plist,
             const Amanzi::Comm_ptr_type& comm,             
             const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
             Amanzi::State& s);


} // namespace ATS

#endif
