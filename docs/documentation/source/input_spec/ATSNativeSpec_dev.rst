ATS Native XML Input Specification V1
#######################################

.. contents:: **Table of Contents**

  
Syntax of the Specification
#######################################

* Input specification for each ParameterList entry consists of two parts.  
  First, a bulleted list defines the usage syntax and available options.  
  This is followed by example snipets of XML code to demonstrate usage.

* In many cases, the input specifies data for a particular parameterized model, and ATS 
  supports a number of parameterizations.  
  For example, initial data might be uniform (the value is required), or linear in y (the value 
  and its gradient are required).  
  Where ATS supports a number of parameterized models for quantity Z, the available 
  models will be listed by name, and then will be described in the subsequent section.  
  For example, the specification for an `"X`" list might begin with the following:

  * `"Y`" ``[string]`` **"default_value"**, `"other`", `"valid`", `"options`"

  * Z ``[Z-spec]`` Model for Z, choose exactly one of the following: (1) `"z1`", or (2) `"z2`" (see below) 

Here, an `"X`" is defined by a `"Y`" and a `"Z`".  
The `"Y`" is a string parameter but the `"Z`" is given by a model (which will require its own set of parameters).
The options for `"Z`" will then be described as a spec:

 * `"z1`" applies model z1.  Requires `"z1a`" ``[string]``

 * `"z2`" applies model z2.  Requires `"z2a`" ``[double]`` and `"z2b`" ``[int]``

An example of using such a specification:

.. code-block:: xml

    <ParameterList name="X">
      <Parameter name="Y" type="string" value="hello"/>
      <ParameterList name="z2">
        <Parameter name="z2a" type="double" value="0.7"/>
        <Parameter name="z2b" type="int" value="3"/>
      </ParameterList>   
    </ParameterList>   
 
Here, the user is defining X with Y="hello", and Z will be a z2 constructed with z2a=0.7 and z2b=3.

Conventions:

* Reserved keywords and labels are `"quoted and italicized`" -- these
  labels or values of parameters in user-generated input files must
  match (using XML matching rules) the specified or allowable values.

* User-defined labels are indicated with ALL-CAPS, and are meant to
  represent a typical name given by a user - these can be names or
  numbers or whatever serves best the organization of the user input
  data.

* Bold values are default values, and are used if the Parameter
  is not provided.


Symbol Index
#############

:math:`|E|` | volume of a cell :math:`[m^X]` (where :math:`X` is the dimension of the mesh)
:math:`g` | gravitational acceleration vector :math:`[m s^-2]`
:math:`h` | ponded depth, or the water head over the surface :math:`[m]`
:math:`` | alternative, in context of the subsurface, water head :math:`[m]`
:math:`h_{snow}` | snow depth :math:`[m]`
:math:`K` | absolute permeability :math:`[m^2]`
:math:`k_r` | relative permeability :math:`[-]`
:math:`n_X` | molar density of phase X :math:`[mol m^-3]`
:math:`p` | pressure of the liquid phase :math:`[Pa]`
:math:`P_{s,r}` | precipitation of rain or snow, noting that snow is always a precipitation rate in snow-water-equivalent (SWE) basis.  :math:`[m s^-1]`
:math:`Q_w` | mass source of water :math:`[mol s^-1]`
:math:`s_X` | saturation of phase X :math:`[-]`
:math:`t` | time variable :math:`[s]`
:math:`z` | elevation :math:`[m]`
:math:`\nu` | dynamic viscosity of water :math:`[Pa s]`
:math:`\phi` | porosity of the soil :math:`[-]`
:math:`\rho` | mass density of a phase :math:`[kg m^-3]`
:math:`\Theta` | extensive water content of a cell :math:`[mol]`

   

  
Main
#######################################

The `"main`" ParameterList frames the entire input spec, and must contain
one sublist for each of the following sections.

* `"mesh`" ``[mesh-spec]``  See the Mesh_ spec.

* `"regions`" ``[list]``

  List of multiple Region_ specs, each in its own sublist named uniquely by the user.

* `"coordinator`" ``[coordinator-spec]``  See the Coordinator_ spec.

* `"visualization`" ``[visualization-spec]`` A Visualization_ spec for the main mesh/domain.

* `"visualization XX`" ``[visualization-spec]``

  Potentially more than one other Visualization_ specs, one for each domain `"XX`".  e.g. `"surface`"

* `"checkpoint`" ``[checkpoint-spec]`` A Checkpoint_ spec.

* `"observations`" ``[observation-spec]`` An Observation_ spec.

* `"PKs`" ``[list]``

  A list containing exactly one sublist, a PK_ spec with the top level PK.

* `"state`" ``[list]`` A State_ spec.

  
Mesh
#####

 Simple wrapper that takes a ParameterList and generates all needed meshes.
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
    




Region
##########


  Region: a geometric or discrete subdomain (abstract)

Regions are geometrical constructs used in Amanzi to define subsets of
the computational domain in order to specify the problem to be solved, and the
output desired. Regions may represents zero-, one-, two- or three-dimensional
subsets of physical space.  for a three-dimensional problem, the simulation
domain will be a three-dimensional region bounded by a set of two-dimensional
regions.  If the simulation domain is N-dimensional, the boundary conditions
must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled *All*, which is the 
entire simulation domain. Currently, the unstructured framework does
not support the *All* region, but it is expected to do so in the
near future.

Amanzi supports parameterized forms for a number of analytic shapes, as well
as more complex definitions based on triangulated surface files.


ONE OF:
* `"region: box`" ``[region-box-spec]``
OR:
* `"region: plane`" ``[region-plane-spec]``
OR:
* `"region: labeled set`" ``[region-labeled-set-spec]``
OR:
* `"region: color function`" ``[region-color-function-spec]``
OR:
* `"region: point`" ``[region-point-spec]``
OR:
* `"region: logical`" ``[region-logical-spec]``
OR:
* `"region: polygon`" ``[region-polygon-spec]``
OR:
* `"region: enumerated`" ``[region-enumerated-spec]``
OR:
* `"region: all`" ``[list]``
OR:
* `"region: boundary`" ``[region-boundary-spec]``
OR:
* `"region: box volume fractions`" ``[region-box-volume-fractions-spec]``
OR:
* `"region: line segment`" ``[region-line-segment-spec]``
END


Notes:

- Surface files contain labeled triangulated face sets.  The user is
  responsible for ensuring that the intersections with other surfaces
  in the problem, including the boundaries, are *exact* (*i.e.* that
  surface intersections are *watertight* where applicable), and that
  the surfaces are contained within the computational domain.  If
  nodes in the surface fall outside the domain, the elements they
  define are ignored.

  Examples of surface files are given in the *Exodus II* file 
  format here.

- Region names must NOT be repeated.

Example:

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="regions">
       <ParameterList name="TOP SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 5}"/>
           <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 8}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="MIDDLE SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 3}"/>
           <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 5}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BOTTOM SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 0}"/>
           <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 3}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="INFLOW SURFACE">
         <ParameterList name="region: labeled set">
           <Parameter name="label"  type="string" value="sideset_2"/>
           <Parameter name="file"   type="string" value="F_area_mesh.exo"/>
           <Parameter name="format" type="string" value="Exodus II"/>
           <Parameter name="entity" type="string" value="face"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="OUTFLOW PLANE">
         <ParameterList name="region: plane">
           <Parameter name="point" type="Array(double)" value="{0.5, 0.5, 0.5}"/>
           <Parameter name="normal" type="Array(double)" value="{0, 0, 1}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BLOODY SAND">
         <ParameterList name="region: color function">
           <Parameter name="file" type="string" value="F_area_col.txt"/>
           <Parameter name="value" type="int" value="25"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="FLUX PLANE">
         <ParameterList name="region: polygon">
           <Parameter name="number of points" type="int" value="5"/>
           <Parameter name="points" type="Array(double)" value="{-0.5, -0.5, -0.5, 
                                                                  0.5, -0.5, -0.5,
                                                                  0.8, 0.0, 0.0,
                                                                  0.5,  0.5, 0.5,
                                                                 -0.5, 0.5, 0.5}"/>
          </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, *TOP SESCTION*, *MIDDLE SECTION* and *BOTTOM SECTION*
are three box-shaped volumetric regions. *INFLOW SURFACE* is a
surface region defined in an Exodus II-formatted labeled set
file and *OUTFLOW PLANE* is a planar region. *BLOODY SAND* is a volumetric
region defined by the value 25 in color function file.





Point
======
 RegionPoint: a point in space.
List *region: point* defines a point in space. 
This region consists of cells containing this point.

* `"coordinate`" ``[Array(double)]`` Location of point in space.

Example:

.. code-block:: xml

   <ParameterList name="DOWN_WIND150"> <!-- parent list defining the name -->
     <ParameterList name="region: point">
       <Parameter name="coordinate" type="Array(double)" value="{-150.0, 0.0, 0.0}"/>
     </ParameterList>
   </ParameterList>




Box
======
 RegionBox: a rectangular region in space, defined by two corners

List *region: box* defines a region bounded by coordinate-aligned
planes. Boxes are allowed to be of zero thickness in only one
direction in which case they are equivalent to planes.

* `"low coordinate`" ``[Array(double)]`` Location of the boundary point with the lowest coordinates.

* `"high coordinate`" ``[Array(double)]`` Location of the boundary points with the highest coordinates.

Example:

.. code-block:: xml

   <ParameterList name="WELL">  <!-- parent list -->
     <ParameterList name="region: box">
       <Parameter name="low coordinate" type="Array(double)" value="{-5.0,-5.0, -5.0}"/>
       <Parameter name="high coordinate" type="Array(double)" value="{5.0, 5.0,  5.0}"/>
     </ParameterList>
   </ParameterList>
  



Plane
======
 RegionPlane: A planar (infinite) region in space, defined by a point and a normal.
List *region: plane* defines a plane using a point lying on the plane and normal to the plane.

* `"normal`" ``[Array(double)]`` Normal to the plane.

* `"point`" ``[Array(double)]`` Point in space.

Example:

.. code-block:: xml

   <ParameterList name="TOP_SECTION"> <!-- parent list -->
     <ParameterList name="region: plane">
       <Parameter name="point" type="Array(double)" value="{2, 3, 5}"/>
       <Parameter name="normal" type="Array(double)" value="{1, 1, 0}"/>
       <ParameterList name="expert parameters">
         <Parameter name="tolerance" type="double" value="1.0e-05"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>




Labeled Set
============
 RegionLabeledSet: A region defined by a set of mesh entities in a mesh file
The list *region: labeled set* defines a named set of mesh entities
existing in an input mesh file. This is the same file that contains
the computational mesh. The name of the entity set is given
by *label*.  For example, a mesh file in the Exodus II
format can be processed to tag cells, faces and/or nodes with
specific labels, using a variety of external tools. Regions based
on such sets are assigned a user-defined label for Amanzi, which may
or may not correspond to the original label in the exodus file.
Note that the file used to express this labeled set may be in any
Amanzi-supported mesh format (the mesh format is specified in the
parameters for this option).  The *entity* parameter may be
necessary to specify a unique set.  For example, an Exodus file
requires *cell*, *face* or *node* as well as a label (which is
an integer).  The resulting region will have the dimensionality 
associated with the entities in the indicated set. 

* `"label`" ``[string]`` Set per label defined in the mesh file.

* `"file`" ``[string]`` File name.

* `"format`" ``[string]`` Currently, we only support mesh files in the "Exodus II" format.

* `"entity`" ``[string]`` Type of the mesh object (cell, face, etc).

Example:

.. code-block:: xml

   <ParameterList name="AQUIFER">
     <ParameterList name="region: labeled set">
       <Parameter name="entity" type="string" value="cell"/>
       <Parameter name="file" type="string" value="porflow4_4.exo"/>
       <Parameter name="format" type="string" value="Exodus II"/>
       <Parameter name="label" type="string" value="1"/>
     </ParameterList>
   </ParameterList>




Color Function
===============
** DOC GENERATION ERROR: file not found ' RegionColorFunction ' **


Coordinator
############

 Coordinator: Simulation controller and top-level driver

In the `"cycle driver`" sublist, the user specifies global control of the
simulation, including starting and ending times and restart options.
 
* `"start time`" ``[double]`` **0.** Specifies the start of time in model time.
 
* `"start time units`" ``[string]`` **"s"** One of `"s`", `"d`", or `"yr`"

* `"end time`" ``[double]`` Specifies the end of the simulation in model time.
 
* `"end time units`" ``[string]`` **"s"** One of `"s`", `"d`", or `"yr`"

* `"end cycle`" ``[int]`` **optional** If provided, specifies the end of the
   simulation in timestep cycles.

* `"restart from checkpoint file`" ``[string]`` **optional** If provided,
   specifies a path to the checkpoint file to continue a stopped simulation.

* `"wallclock duration [hrs]`" ``[double]`` **optional** After this time, the
   simulation will checkpoint and end.

* `"required times`" ``[io-event-spec]`` **optional** An IOEvent_ spec that
   sets a collection of times/cycles at which the simulation is guaranteed to
   hit exactly.  This is useful for situations such as where data is provided
   at a regular interval, and interpolation error related to that data is to
   be minimized.

* `"PK tree`" ``[pk-type-spec-list]`` List of length one, the top level PK spec.
   
Note: Either `"end cycle`" or `"end time`" are required, and if
both are present, the simulation will stop with whichever arrives
first.  An `"end cycle`" is commonly used to ensure that, in the case
of a time step crash, we do not continue on forever spewing output.

``[pk-type-spec]`` is a pk type and a list of subpks.
* `"PK type`" ``[string]`` One of the registered PK types
* `"sub PKs`" ``[pk-type-spec-list]`` **optional** If there are sub pks, list them.

Example:

.. code-block::xml

   <!-- simulation control -->
   <ParameterList name="coordinator">
     <Parameter  name="end cycle" type="int" value="6000"/>
     <Parameter  name="start time" type="double" value="0."/>
     <Parameter  name="start time units" type="string" value="s"/>
     <Parameter  name="end time" type="double" value="1"/>
     <Parameter  name="end time units" type="string" value="yr"/>
     <ParameterList name="required times">
       <Parameter name="start period stop" type="Array(double)" value="{0,-1,86400}" />
     </ParameterList>
   </ParameterList>



   

Visualization
##############

A user may request periodic writes of field data for the purposes of
visualization in the `"visualization`" sublists.  ATS accepts a visualization
list for each domain/mesh, including surface and column meshes.  These are in
separate ParameterLists, entitled `"visualization`" for the main mesh, and
`"visualization surface`" on the surface mesh.  It is expected that, for any
addition meshes, each will have a domain name and therefore admit a spec of
the form: `"visualization DOMAIN-NAME`".

 Visualization: a class for controlling simulation output.

Each list contains all parameters as in a IOEvent_ spec, and also:

* `"file name base`" ``[string]`` **visdump_DOMAIN_data**

* `"dynamic mesh`" ``[bool]`` **false** Write mesh data for every
  visualization dump; this facilitates visualizing deforming meshes.

* `"time units`" ``[string]`` **s** A valid time unit to convert time
  into for output files.  One of `"s`", `"d`", `"y`", or `"yr 365`"
  
INCLUDES:
* ``[io-event-spec]`` An IOEvent_ spec


Example:

.. code-block:: xml

  <ParameterList name="visualization">
    <Parameter name="file name base" type="string" value="visdump_data"/>
  
    <Parameter name="cycles start period stop" type="Array(int)" value="{{0, 100, -1}}" />
    <Parameter name="cycles" type="Array(int)" value="{{999, 1001}}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{{0.0, 10.0, 100.0}}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{{100.0, 25.0, -1.0}}"/>
    <Parameter name="times" type="Array(double)" value="{{101.0, 303.0, 422.0}}"/>

    <Parameter name="dynamic mesh" type="bool" value="false"/>
  </ParameterList>




  
Checkpoint
##############

A user may request periodic dumps of ATS Checkpoint Data in the
`"checkpoint`" sublist.  The user has no explicit control over the
content of these files, but has the guarantee that the ATS run will be
reproducible (with accuracies determined by machine round errors and
randomness due to execution in a parallel computing environment).
Therefore, output controls for Checkpoint Data are limited to file
name generation and writing frequency, by numerical cycle number.
Unlike `"visualization`", there is only one `"checkpoint`" list for
all domains/meshes.

** DOC GENERATION ERROR: file not found ' checkpoint ' **  


 
Observation
##############

 Observable: Collects, reduces, and writes observations during a simulation.
Observations are a localized-in-space but frequent in time view of
data, designed to get at useful diagnostic quantities such as
hydrographs, total water content, quantities at a point, etc.  These
are designed to allow frequent collection in time without saving huge
numbers of visualization files to do postprocessing.  In fact, these
should be though of as orthogonal data queries to visualization -- vis
is pointwise in time but complete in space, while observations are
pointwise/finite in space but complete in time.

A user may request any number of specific observations from ATS.  Each
observation spec involves a field quantity, a functional reduction
operator, a region from which it will extract its source data, and a
list of discrete times for its evaluation.  The observations are
evaluated during the simulation and written to disk.

``[observation-spec]`` consists of the following quantities:

* `"observation output filename`" ``[string]`` user-defined name for the file that the observations are written to.

* `"variable`" ``[string]`` any ATS variable used by any PK, e.g. `"pressure`" or `"surface-water_content`"

* `"region`" ``[string]`` the label of a user-defined region

* `"location name`" ``[string]`` the mesh location of the thing to be measured, i.e. `"cell`", `"face`", or `"node`"

* `"functional`" ``[string]`` the label of a function to apply to the variable across the region.  Valid functionals include:
 - `"observation data: point`" returns the value of the field quantity at a point.  The region and location name must result in a single entity being selected.
 - `"observation data: extensive integral`" returns the sum of an (extensive) variable over the region.  This should be used for extensive quantities such as `"water_content`" or `"energy`".
 - `"observation data: intensive integral`" returns the volume-weighted average of an (intensive) variable over the region.  This should be used for intensive quantities such as `"temperature`" or `"saturation_liquid`".

For flux observations, additional options are available:

* `"direction normalized flux`" ``[bool]`` **optional** Dots the face-normal flux with a vector to ensure fluxes are integrated pointing the same direction.

* `"direction normalized flux direction`" ``[Array(double)]`` **optional**
  Provides the vector to dot the face normal with.  If this is not provided,
  then it is assumed that the faces integrated over are all boundary faces and
  that the default vector is the outward normal direction for each face.

Additionally, each ``[observation-spec]`` contains all parameters as in a IOEvent_ spec, which are used to specify at which times/cycles the observation is collected.

INCLUDES:
* ``[io-event-spec]`` An IOEvent_ spec

  


Example:

.. code-block:: xml
  
  <ParameterList name="observations" type="ParameterList">
    <!-- This measures the hydrograph out the "east" face of the surface domain -->
    <ParameterList name="surface outlet flux" type="ParameterList">
      <Parameter name="variable" type="string" value="surface-mass_flux" />
      <Parameter name="direction normalized flux" type="bool" value="true" />
      <Parameter name="region" type="string" value="east" />
      <Parameter name="functional" type="string" value="observation data: extensive integral" />
      <Parameter name="delimiter" type="string" value=" " />
      <Parameter name="location name" type="string" value="face" />
      <Parameter name="observation output filename" type="string" value="surface_outlet_flux.dat" />
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1.0}" />
    </ParameterList>
    <!-- This measures the total water, in mols, in the entire subsurface domain -->
    <ParameterList name="subsurface water content" type="ParameterList">
      <Parameter name="variable" type="string" value="water_content" />
      <Parameter name="region" type="string" value="computational domain" />
      <Parameter name="functional" type="string" value="observation data: extensive integral" />
      <Parameter name="delimiter" type="string" value=" " />
      <Parameter name="location name" type="string" value="cell" />
      <Parameter name="observation output filename" type="string" value="water_content.dat" />
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1.0}" />
    </ParameterList>
    <!-- This tracks the temperature at a point -->
    <ParameterList name="temperature_probeA" type="ParameterList">
      <Parameter name="variable" type="string" value="temperature" />
      <Parameter name="region" type="string" value="probeA" />
      <Parameter name="functional" type="string" value="observation data: point" />
      <Parameter name="delimiter" type="string" value=" " />
      <Parameter name="location name" type="string" value="cell" />
      <Parameter name="observation output filename" type="string" value="temperature_probeA.dat" />
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1.0}" />
    </ParameterList>
  </ParameterList>





PK
#####

 The interface for a Process Kernel, an equation or system of equations.
A process kernel represents a single or system of partial/ordinary
differential equation(s) or conservation law(s), and is used as the
fundamental unit for coupling strategies.

Implementations of this interface typically are either an MPC
(multi-process coupler) whose job is to heirarchically couple several
other PKs and represent the system of equations, or a Physical PK,
which represents a single equation.

All PKs have the following parameters in their spec:

* `"PK type`" ``[string]``

  The PK type is a special key-word which corresponds to a given class in the PK factory.  See available PK types listed below.

Example:

.. code-block:: xml

  <ParameterList name="PKs">
    <ParameterList name="my cool PK">
      <Parameter name="PK type" type="string" value="my cool PK"/>
       ...
    </ParameterList>
  </ParameterList>

.. code-block:: xml

  <ParameterList name="PKs">
    <ParameterList name="Top level MPC">
      <Parameter name="PK type" type="string" value="strong MPC"/>
       ...
    </ParameterList>
  </ParameterList>



Base PKs
===============

There are several types of PKs, and each PK has its own valid input spec.  However, there are three main types of PKs, from which nearly all PKs derive.  Note that none of these are true PKs and cannot stand alone.

PKPhysicalBase
----------------

 A base class with default implementations of methods for a leaf of the PK tree (a conservation equation, or similar).

``PKPhysicalBase`` is a base class providing some functionality for PKs which
are defined on a single mesh, and represent a single process model.  Typically
all leaves of the PK tree will inherit from ``PKPhysicalBase``.

* `"domain name`" ``[string]`` e.g. `"surface`".

  Domains and meshes are 1-to-1, and the empty string refers to the main domain or mesh.  PKs defined on other domains must specify which domain/mesh they refer to.

* `"primary variable key`" ``[string]`` Sets the primary variable.

  The primary variable associated with this PK, i.e. `"pressure`", `"temperature`", `"surface_pressure`", etc.

* `"initial condition`" ``[initial-conditions-spec]``  See InitialConditions_.


   Indicates that the primary variable field has both CELL and FACE objects, and the FACE values are calculated as the average of the neighboring cells.


NOTE: ``PKPhysicalBase (v)-->`` PKDefaultBase_





PKBDFBase
----------------

 A base class with default implementations of methods for a PK that can be implicitly integrated in time.

``PKBDFBase`` is a base class from which PKs that want to use the ``BDF``
series of implicit time integrators must derive.  It specifies both the
``BDFFnBase`` interface and implements some basic functionality for ``BDF``
PKs.

* `"initial time step`" ``[double]`` **1.**

  The initial timestep size for the PK, this ensures that the initial timestep
  will not be **larger** than this value.

* `"assemble preconditioner`" ``[bool]`` **true** 

  A flag for the PK to not assemble its preconditioner if it is not needed by
  a controlling PK.  This is usually set by the MPC, not by the user.

In the top-most (in the PK tree) PK that is meant to be integrated implicitly,
several additional specs are included.  For instance, in a strongly coupled
flow and energy problem, these specs are included in the ``StrongMPC`` that
couples the flow and energy PKs, not to the flow or energy PK itself.
  
* `"time integrator`" ``[implicit-time-integrator-typed-spec]`` **optional**
  A TimeIntegrator_.  Note that this is only provided if this PK is not
  strongly coupled to other PKs.

* `"preconditioner`" ``[preconditioner-typed-spec]`` **optional** is a Preconditioner_ spec.
  Note that this is only used if this PK is not strongly coupled to other PKs.

  This spec describes how to form the (approximate) inverse of the preconditioner.
  
NOTE: ``PKBDFBase  (v)-->`` PKDefaultBase_





PKPhysicalBDFBase
-------------------

 Standard base for most implemented PKs, this combines both domains/meshes of PKPhysicalBase and BDF methods of PKBDFBase.

A base class for all PKs that are both physical, in the sense that they
implement an equation and are not couplers, and support the implicit
integration interface.  This largely just supplies a default error norm based
on error in conservation relative to the extent of the conserved quantity.

* `"absolute error tolerance`" ``[double]`` **1.0** Absolute tolerance,
  :math:`a_tol` in the equation below.  Note that this default is often
  overridden by PKs with more physical values, and very rarely are these set
  by the user.

* `"relative error tolerance`" ``[double]`` **1.0** Relative tolerance,
  :math:`r_tol` in the equation below.  Note that this default is often
  overridden by PKs with more physical values, and very rarely are these set
  by the user.

* `"flux error tolerance`" ``[double]`` **1.0** Relative tolerance on the
  flux.  Note that this default is often overridden by PKs with more physical
  values, and very rarely are these set by the user.

By default, the error norm used by solvers is given by:

:math:`ENORM(u, du) = |du| / ( a_tol + r_tol * |u| )`

The defaults here are typically good, or else good defaults are set in the
code, so these need not be supplied.


NOTE: ``PKPhysicalBDFBase -->`` PKBDFBase_
      ``PKPhysicalBDFBase -->`` PKPhysicalBase_
      ``PKPhysicalBDFBase (v)-->`` PKDefaultBase_




Physical PKs
===============

Physical PKs are the physical capability implemented within ATS.

Flow PKs
-----------

Richards PK
^^^^^^^^^^^^^^^

 Two-phase, variable density Richards equation.

Solves Richards equation:

.. math::
  \frac{\partial \Theta}{\partial t} - \nabla \frac{k_r n_l}{\mu} K ( \nabla p + \rho g \cdot \hat{z} ) = Q_w


Includes options from:

* ``[pk-physical-default-spec]`` PKPhysicalDefault_
# ``[pk-bdf-default-spec]`` PKBDFDefault_
* ``[pk-physical-bdf-default-spec]`` PKPhysicalBDFDefault_

Other variable names, typically not set as the default is basically always good:

* `"conserved quantity key`" ``[string]`` **DOMAIN-water_content** Typically not set, default is good. ``[mol]``

* `"mass density key`" ``[string]`` **DOMAIN-mass_density_liquid** liquid water density ``[kg m^-3]``

* `"molar density key`" ``[string]`` **DOMAIN-molar_density_liquid** liquid water density ``[mol m^-3]``

* `"permeability key`" ``[string]`` **DOMAIN-permeability** permeability of the soil medium ``[m^2]``

* `"conductivity key`" ``[string]`` **DOMAIN-relative_permeability** scalar coefficient of the permeability ``[-]``

* `"upwind conductivity key`" ``[string]`` **DOMAIN-upwind_relative_permeability** upwinded (face-based) scalar coefficient of the permeability.  Note the units of this are strange, but this represents :math:`\frac{n_l k_r}{\mu}`  ``[mol kg^-1 s^1 m^-2]``

* `"darcy flux key`" ``[string]`` **DOMAIN-mass_flux** mass flux across a face ``[mol s^-1]``

* `"darcy flux direction key`" ``[string]`` **DOMAIN-mass_flux_direction** direction of the darcy flux (used in upwinding :math:`k_r`) ``[??]``

* `"darcy velocity key`" ``[string]`` **DOMAIN-darcy_velocity** darcy velocity vector, interpolated from faces to cells ``[m s^-1]``

* `"darcy flux key`" ``[string]`` **DOMAIN-mass_flux** mass flux across a face ``[mol s^-1]``

* `"saturation key`" ``[string]`` **DOMAIN-saturation_liquid** volume fraction of the liquid phase ``[-]``

Discretization control:

* `"diffusion`" ``[list]`` An PDE_Diffusion_ spec describing the (forward) diffusion operator

* `"diffusion preconditioner`" ``[list]`` An PDE_Diffusion_ spec describing the diffusive parts of the preconditioner.

* `"linear solver`" ``[linear-solver-typed-spec]`` **optional** is a LinearSolver_ spec.  Note
  that this is only used if this PK is not strongly coupled to other PKs.

Boundary conditions:

//* `"boundary conditions`" ``[subsurface-flow-bc-spec]`` Defaults to Neuman, 0 normal flux.  See `Flow-specific Boundary Conditions`_

Physics control:

* `"permeability rescaling`" ``[double]`` **1** Typically 1e7 or order :math:`sqrt(K)` is about right.  This rescales things to stop from multiplying by small numbers (permeability) and then by large number (:math:`\rho / \mu`).

* `"permeability type`" ``[string]`` **'scalar'** The permeability type can be 'scalar', 'horizontal and vertical', 'diagonal tensor', or 'full tensor'. This key is placed in state->field evaluators->permeability. The 'scalar' option requires 1 permeability value, 'horizontal and vertical' requires 2 values, 'diagonal tensor' requires 2 (2D) or 3 (3D) values, and 'full tensor' requires 3 (2D) or 6 (3D) values. The ordering of the permeability values in the input script is important: 'horizontal and vertical'={xx/yy,zz}, 'diagonal tensor'={xx,yy} or {xx,yy,zz}, 'full tensor'={xx,yy,xy/yx} or {xx,yy,zz,xy/yx,xz/zx,yz/zy}.

* `"water retention evaluator`" ``[wrm-evaluator-spec]`` The WRM.  This needs to go away!

This PK additionally requires the following:

EVALUATORS:
- `"conserved quantity`"
- `"mass density`"
- `"molar density`"
- `"permeability`"
- `"conductivity`"
- `"saturation`"




Permafrost Flow PK
^^^^^^^^^^^^^^^^^^^^

Overland Flow, head primary variable PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Overland Flow, pressure primary variable, PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 Overland flow using the diffusion wave equation.

Solves the diffusion wave equation for overland flow with pressure as a primary variable:

.. math::
  \frac{\partial \Theta}{\partial t} - \nabla n_l k \nabla h(p) = Q_w


Options:

Variable naming:

* `"domain`" ``[string]`` **"surface"**  Defaults to the extracted surface mesh.

* `"primary variable`" ``[string]`` The primary variable associated with this PK, typically `"DOMAIN-pressure`"


Other variable names, typically not set as the default is basically always good:

* `"conserved quantity key`" ``[string]`` **"DOMAIN-water_content"** The conserved quantity name.

Discretization control:

* `"diffusion`" ``[list]`` An PDE_Diffusion_ spec describing the (forward) diffusion operator

* `"diffusion preconditioner`" ``[list]`` An PDE_Diffusion_ spec describing the diffusive parts of the preconditioner.

Time integration and timestep control:

* `"initial time step`" ``[double]`` **1.** Max initial time step size ``[s]``.

* `"time integrator`" ``[implicit-time-integrator-typed-spec]`` **optional** is a TimeIntegrator_ spec.
  Note that this is only used if this PK is not strongly coupled to other PKs.

* `"linear solver`" ``[linear-solver-typed-spec]`` **optional** is a LinearSolver_ spec.  Note
  that this is only used if this PK is not strongly coupled to other PKs.

* `"preconditioner`" ``[preconditioner-typed-spec]`` **optional** is a Preconditioner_ spec.
  Note that this is only used if this PK is not strongly coupled to other PKs.

* `"initial condition`" ``[initial-conditions-spec]`` See InitialConditions_.

Error control:

* `"absolute error tolerance`" ``[double]`` **550.** Defaults to 1 cm of water.  A small, but significant, amount of water.

* `"relative error tolerance`" ``[double]`` **1** Take the error relative to the amount of water present in that cell.

* `"flux tolerance`" ``[double]`` **1** Multiplies the error in flux (on a face)
  relative to the min of water in the neighboring cells.  Typically only
  changed if infiltration is very small and the boundary condition is not
  converging, at which point it can be decreased by an order of magnitude at a
  time until the boundary condition is satisfied.

Boundary conditions:

xx* `"boundary conditions`" ``[surface-flow-bc-spec]`` Defaults to Neuman, 0 normal flux.


May inherit options from PKPhysicalBDFBase_.





Snow Distribution PK
^^^^^^^^^^^^^^^^^^^^


Energy PKs
-----------

Advection Diffusion PK
^^^^^^^^^^^^^^^^^^^^^^^

Energy Base PK
^^^^^^^^^^^^^^^^^^^^^^^

Two-Phase subsurface Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three-Phase subsurface Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three-Phase subsurface Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Surface Ice Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Surface Energy Balance PKs
------------------------------


Surface Energy Balance / Snow -- Monolithic Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Surface Energy Balance -- Generic Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Biogeochemistry
-----------------


Biogeochemistry -- Monolithic Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Deformation
-------------


Volumetric Deformation
^^^^^^^^^^^^^^^^^^^^^^



MPCs
===============

MPCs couple other PKs, and are the non-leaf nodes in the PK tree.

WeakMPC
----------

StrongMPC
----------

Physical MPCs
===============
 coupling is an art, and requires special off-diagonal work.  Physical MPCs can derive from default MPCs to provide special work.

Coupled Water MPC
--------------------

 MPCCoupledWater: coupler which integrates surface and subsurface flow.

Couples Richards equation to surface water through continuity of both pressure and fluxes.

Currently requires that the subsurface discretization is a face-based
discretization, i.e. one of the MFD methods.  Then the surface equations are
directly added into the subsurface discrete equations.

* `"PKs order`" ``[Array(string)]`` Supplies the names of the coupled PKs.
  The order must be {subsurface_flow_pk, surface_flow_pk} (subsurface first).

* `"linear solver`" ``[linear-solver-typed-spec]`` **optional** is a LinearSolver_ spec.  Note
  that this is only used if this PK is not strongly coupled to other PKs.

* `"preconditioner`" ``[preconditioner-typed-spec]`` **optional** is a Preconditioner_ spec.
  Note that this is only used if this PK is not strongly coupled to other PKs.

* `"water delegate`" ``[list]`` 



 Globalization hacks to deal with nonlinearity around the appearance/disappearance of surface water.

 The water delegate works to eliminate discontinuities/strong nonlinearities
 when surface cells shift from dry to wet (i.e. the surface pressure goes
 from < atmospheric pressure to > atmospheric pressure.

 These methods work to alter the predictor around this nonlinearity.

 - `"modify predictor with heuristic`" ``[bool]`` **false** This simply
   limits the prediction to backtrack to just above atmospheric on both the
   first and second timesteps that take us over atmospheric.

 - `"modify predictor damp and cap the water spurt`" ``[bool]`` **false** The
   second both limits (caps) and damps all surface cells to ensure that all
   nearby cells are also not overshooting.  This is the preferred method.
    
 These methods work to alter the preconditioned correction for the same
 reasons described above.

 - `"global water face limiter`" ``[default]`` **INF** This is simply a limit
   to the maximum allowed size of the correction (in [Pa]) on all faces.  Any
   correction larger than this is set to this.

 - `"cap the water spurt`" ``[bool]`` **false** If a correction takes the
   pressure on a surface cell from below atmospheric (dry) to above (wet),
   the correction is set to a value which results in the new iterate to being
   CAP_SIZE over atmospheric.

 - `"damp the water spurt`" ``[bool]`` **false** A damping factor (less than
   one) is calculated to multiply the correction such that the largest
   correction takes a cell to just above atmospheric.  All faces (globally)
   are affected.
  
 - `"damp and cap the water spurt`" ``[bool]`` **false** None of the above
   should really be used.  Capping, when the cap is particularly severe,
   results in faces whose values are very out of equilibrium with their
   neighboring cells which are not capped.  Damping results in a tiny
   timestep in which, globally, at MOST one face can go from wet to dry.
   This looks to do a combination, in which all things are damped, but faces
   that are initially expected to go from dry to wet are pre-scaled to ensure
   that, when damped, they are also (like the biggest change) allowed to go
   from dry to wet (so that multiple cells can wet in the same step).  This
   is the preferred method.

 In these methods, the following parameters are useful:

 - `"cap over atmospheric`" ``[double]`` **100 Pa** This sets the max size over
   atmospheric to which things are capped or damped.
  
 



Subsurface MPC
--------------------

Permafrost MPC
--------------------

State
##############

State consists of two sublists, one for evaluators and the other for
atomic constants.  The latter is currently called `"initial
conditions`", which is a terrible name which must be fixed.

example:

.. code-block:: xml
                
  <ParameterList name="state">
    <ParameterList name="field evaluators">
      ...
    </ParameterList>
    <ParameterList name="initial conditions">
      ...
    </ParameterList>
  </ParameterList>

 

Field Evaluators
=================

Many field evaluators exist, but most derive from one of four base types.

Field Evaluator Base Classes
-------------------------------

PrimaryVariableEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^

SecondaryVariableEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^

SecondaryVariablesEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^

IndependentVariableEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While these provide base functionality, all of the physics are in the
following derived classes.

Water Content
-----------------

Water content is the conserved quantity in most flow equations, including
Richard's equation with and without ice.  A variety of evaluators are provided
for inclusion of multiple phases.

RichardsWaterContentEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 The Richards water content evaluator is an algebraic evaluator for liquid only water content
  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

.. math::
  Theta = n * s * phi * cell volume

``[field-evaluator-type-richards-water-content-spec]``

* `"porosity key`" ``[string]`` **DOMAIN-porosity** 
* `"molar density liquid key`" ``[string]`` **DOMAIN-molar_density_liquid** 
* `"saturation liquid key`" ``[string]`` **DOMAIN-saturation_liquid** 
* `"cell volume key`" ``[string]`` **DOMAIN-cell_volume**

EVALUATORS:
- `"porosity`"
- `"molar density liquid`"
- `"saturation liquid`"
- `"cell volume`"



RichardsWaterContentWithVaporEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
** DOC GENERATION ERROR: file not found ' richards_water_content_with_vapor_evaluator ' **
PermafrostWaterContentEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
** DOC GENERATION ERROR: file not found ' permafrost_water_content ' **


Surface Water potential surfaces
---------------------------------

Evaluators for 

SurfaceElevation
^^^^^^^^^^^^^^^^^^
 MeshedElevationEvaluator: evaluates the elevation (z-coordinate) and slope magnitude of a mesh.
Evaluator type: `"meshed elevation`"

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h z|``

* `"elevation key`" ``[string]`` **elevation** Name the elevation variable. [m]
* `"slope magnitude key`" ``[string]`` **slope_magnitude** Name the elevation variable. [-]
* `"dynamic mesh`" ``[bool]`` **false** Lets the evaluator know that the elevation changes in time, and adds the `"deformation`" dependency.
* `"parent domain name`" ``[string]`` **DOMAIN** Domain name of the parent mesh, which is the 3D version of this domain.  Attempts to generate an intelligent default by stripping "surface" from this domain.

Example:

.. code-block:: xml

  <ParameterList name="elevation">
    <Parameter name="evaluator type" type="string" value="meshed elevation"/>
  </ParameterList>




SurfacePotential
^^^^^^^^^^^^^^^^^^^
 PresElevEvaluator: evaluates h + z
Evaluator type: ""

.. math::
  h + z

* `"my key`" ``[string]`` **pres_elev** Names the surface water potential variable, h + z [m]
* `"height key`" ``[string]`` **ponded_depth** Names the height variable. [m]
* `"elevation key`" ``[string]`` **elevation** Names the elevation variable. [m]


NOTE: This is a legacy evaluator, and is not in the factory, so need not be in
the input spec.  However, we include it here because this could easily be
abstracted for new potential surfaces, kinematic wave, etc, at which point it
would need to be added to the factory and the input spec.

NOTE: This could easily be replaced by a generic AdditiveEvaluator_




SnowSurfacePotential
^^^^^^^^^^^^^^^^^^^^^^
 PresElevEvaluator: evaluates h + z
Evaluator type: "snow skin potential"

.. math::
  h + z + h_{{snow}} + dt * P_{{snow}}

* `"my key`" ``[string]`` **snow_skin_potential** Names the potential variable evaluated [m]
* `"ponded depth key`" ``[string]`` **ponded_depth** Names the surface water depth variable. [m]
* `"snow depth key`" ``[string]`` **snow_depth** Names the snow depth variable. [m]
* `"precipitation snow key`" ``[string]`` **precipitation_snow** Names the snow precipitation key. [m]
* `"elevation key`" ``[string]`` **elevation** Names the elevation variable. [m]
* `"dt factor`" ``[double]`` A free-parameter factor for providing a time scale for diffusion of snow precipitation into low-lying areas.  Typically on the order of 1e4-1e7. This timestep times the wave speed of snow provides an approximate length of how far snow precip can travel.  Extremely tunable! [s]

NOTE: This is equivalent to a generic AdditiveEvaluator_

Example:

.. code-block:: xml

  <ParameterList name="snow_skin_potential" type="ParameterList">
    <Parameter name="field evaluator type" type="string" value="snow skin potential" />
    <Parameter name="dt factor" type="double" value="864000.0" />
  </ParameterList>






Generic Evaluators
---------------------------------

Several generic evaluators are provided.

AdditiveEvaluator
^^^^^^^^^^^^^^^^^^^^^^


MultiplicativeEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^




InitialConditions
=================

Initial condition specs are used in two places:

* within the PK_ spec which describes the initial condition of primary variables (true
  initial conditions), and

* in the `"initial conditions`" sublist of state, in which the value
  of atomic constants are provided (not really initial conditions and
  should be renamed).  These atomic values are not controlled by
  evaluators, and are not included in the DaG.  Likely these should be
  removed entirely.
  
Initialization of constant scalars
------------------------------------

A constant scalar field is the global (with respect to the mesh)
constant.  At the moment, the set of such fields includes atmospheric
pressure.  The initialization requires to provide a named sublist with
a single parameter `"value`".

.. code-block:: xml

  <ParameterList name="fluid_density">
    <Parameter name="value" type="double" value="998.0"/>
  </ParameterList>


Initialization of constant vectors
------------------------------------

A constant vector field is the global (with respect to the mesh)
vector constant.  At the moment, the set of such vector constants
includes gravity.  The initialization requires to provide a named
sublist with a single parameter `"Array(double)`". In two dimensions,
is looks like

.. code-block:: xml

  <ParameterList name="gravity">
    <Parameter name="value" type="Array(double)" value="{0.0, -9.81}"/>
  </ParameterList>


Initialization of scalar fields
------------------------------------

A variable scalar field is defined by a few functions (labeled for instance,
`"Mesh Block i`" with non-overlapping ranges. 
The required parameters for each function are `"region`", `"component`",
and the function itself.

.. code-block:: xml

  <ParameterList name="porosity"> 
    <ParameterList name="function">
      <ParameterList name="Mesh Block 1">
        <Parameter name="region" type="string" value="Computational domain"/>
        <Parameter name="component" type="string" value="cell"/>
        <ParameterList name="function">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="0.2"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
      <ParameterList name="Mesh Block 2">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>


Initialization of tensor fields
------------------------------------
 
A variable tensor (or vector) field is defined similarly to a variable
scalar field.  The difference lies in the definition of the function
which is now a multi-values function.  The required parameters are
`"Number of DoFs`" and `"Function type`".

.. code-block:: xml

  <ParameterList name="function">
    <Parameter name="Number of DoFs" type="int" value="2"/>
    <Parameter name="Function type" type="string" value="composite function"/>
    <ParameterList name="DoF 1 Function">
      <ParameterList name="function-constant">
        <Parameter name="value" type="double" value="1.9976e-12"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="DoF 2 Function">
      <ParameterList name="function-constant">
        <Parameter name="value" type="double" value="1.9976e-13"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Initialization from a file
------------------------------------

Some data can be initialized from files. Additional sublist has to be
added to named sublist of the `"state`" list with the file name and
the name of attribute.  For a serial run, the file extension must be
`".exo`".  For a parallel run, it must be `".par`".  Here is an
example:

.. code-block:: xml

  <ParameterList name="permeability">
    <ParameterList name="exodus file initialization">
      <Parameter name="file" type="string" value="mesh_with_data.exo"/>
      <Parameter name="attribute" type="string" value="perm"/>
    </ParameterList>
  </ParameterList>



example:

.. code-block:: xml

  <ParameterList name="state">
    <ParameterList name="initial conditions">
      <ParameterList name="fluid_density">
        <Parameter name="value" type="double" value="998.0"/>
      </ParameterList>

      <ParameterList name="fluid_viscosity">
        <Parameter name="value" type="double" value="0.001"/>
      </ParameterList>

      <ParameterList name="gravity">
        <Parameter name="value" type="Array(double)" value="{0.0, -9.81}"/>
      </ParameterList>

    </ParameterList>
  </ParameterList>



BoundaryConditions
===================



In general, boundary conditions are provided in a heirarchical list by
boundary condition type, then functional form.  Boundary condition specs are
split between two types -- those which require a user-provided function
(i.e. Dirichlet data, etc) and those which do not (i.e. zero gradient
conditions).

A list of conditions might pull in both Dirichlet and Neumann data on
different regions, or use different functions on different regions.  The
following example illustrates how boundary conditions are prescribed across
the domain for a typical PK:

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="DIRICHLET_TYPE">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="DIRICHLET_FUNCTION_NAME">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
     <ParameterList name="BC east">
       <Parameter name="regions" type="Array(string)" value="{east}"/>
       <ParameterList name="DIRICHLET_FUNCTION_NAME">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="102325."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
   <ParameterList name="mass flux">
     <ParameterList name="BC north">
       <Parameter name="regions" type="Array(string)" value="{north}"/>
       <ParameterList name="outward mass flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
   <ParameterList name="zero gradient">
     <ParameterList name="BC south">
       <Parameter name="regions" type="Array(string)" value="{south}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Different PKs populate this general format with different names, replacing
DIRICHLET_TYPE and DIRICHLET_FUNCTION_NAME.
  
 


Flow-specific Boundary Conditions
----------------------------------



Flow boundary conditions must follow the general format shown in
BoundaryConditions_.  Specific conditions implemented include:

Dirichlet (pressure) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for both surface and subsurface flows, this provides pressure data on
boundaries (in [Pa]).

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="pressure">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary pressure">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Neumann (mass flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for both surface and subsurface flows, this provides mass flux data (in [mol m^-2 s^-1], in the outward normal direction) on boundaries.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="mass flux">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward mass flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="-1.e-3"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

 
Seepage face boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A variety of seepage face boundary conditions are permitted for both surface
and subsurface flow PKs.  Typically seepage conditions are of the form:

  - if :math:`q \cdot \hat{n} < 0`, then :math:`q = 0`
  - if :math:`p > p0`, then :math:`p = p0`

This ensures that flow is only out of the domain, but that the max pressure on
the boundary is specified by :math:`p0`.

Example: pressure (for surface or subsurface)

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="seepage face pressure">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary pressure">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Example: head (for surface)
 
.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="seepage face head">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary head">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Additionally, an infiltration flux may be prescribed, which describes the max
flux.  This is for surface faces on which a typical precipitation rate might
be prescribed, to be enforced until the water table rises to the surface, at
which point the precip is turned off and water seeps into runoff.  This
capability is experimental and has not been well tested.

  - if :math:`q \cdot \hat{n} < q0`, then :math:`q = q0`
  - if :math:`p > p_atm`, then :math:`p = p_atm`

Example: seepage with infiltration

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="seepage face with infiltration">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward mass flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="-1.e-5"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

Note it would be straightforward to add both p0 and q0 in the same condition;
this has simply not had a use case yet.


Dirichlet (head) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Used for surface flows, this provides head data (in [m]) on boundaries.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="head">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary head">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0.01"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Fixed level boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For surface flows only.  This fixes the water table at a constant elevation.
It is a head condition that adapts to the surface elevation such that

.. math::
  h = max( h0 - z, 0 )

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="fixed level">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="fixed level">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Zero head gradient boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for surface flows, this is an "outlet" boundary condition which looks to
enforce the condition that

.. math::
  \div h \cdot \hat{n} = 0

for head :math:`h` and outward normal :math:`\hat{n}`.  Note that this is an
"outlet" boundary, in the sense that it should really not be used on a
boundary in which

.. math::
  \div z \cdot \hat{n} > 0.

This makes it a useful boundary condition for benchmark and 2D problems, where
the elevation gradient is clear, but not so useful for DEM-based meshes.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="zero gradient">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Critical depth boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Also for surface flows, this is an "outlet" boundary condition which looks to
set an outward flux to take away runoff.  This condition is given by:

.. math::
  q = \sqrt{g \hat{z}} n_{liq} h^1.5

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="critical depth">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Dynamic boundary condutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The type of boundary conditions maybe changed in time depending on the switch function of TIME.
<ParameterList name="dynamic">
          
     <Parameter name="regions" type="Array(string)" value="{surface west}"/>
     <ParameterList name="switch function">
       <ParameterList name="function-tabular">
         <Parameter name="file" type="string" value="../data/floodplain2.h5" />
         <Parameter name="x header" type="string" value="Time" />
         <Parameter name="y header" type="string" value="Switch" />
         <Parameter name="form" type="Array(string)" value="{constant}"/>
       </ParameterList>
     </ParameterList>
          
     <ParameterList name="bcs">
       <Parameter name="bc types" type="Array(string)" value="{head, mass flux}"/>
       <Parameter name="bc functions" type="Array(string)" value="{boundary head, outward mass flux}"/>

       <ParameterList name="mass flux">
         <ParameterList name="BC west">
           <Parameter name="regions" type="Array(string)" value="{surface west}"/>
           <ParameterList name="outward mass flux">
             <ParameterList name="function-tabular">
               <Parameter name="file" type="string" value="../data/floodplain2.h5" />
               <Parameter name="x header" type="string" value="Time" />
               <Parameter name="y header" type="string" value="Flux" />
               <Parameter name="form" type="Array(string)" value="{linear}"/>
             </ParameterList>
            </ParameterList>
          </ParameterList>
       </ParameterList>

       <ParameterList name="head">  
          <ParameterList name="BC west">
            <Parameter name="regions" type="Array(string)" value="{surface west}"/>
            <ParameterList name="boundary head">
              <ParameterList name="function-tabular">
                 <Parameter name="file" type="string" value="../data/floodplain2.h5" />
                 <Parameter name="x header" type="string" value="Time" />
                 <Parameter name="y header" type="string" value="Head" />
                 <Parameter name="form" type="Array(string)" value="{linear}"/>
               </ParameterList>
            </ParameterList>
          </ParameterList>
        </ParameterList>
     </ParameterList>
                 
 </ParameterList> 
<!-- dynamic -->

 







Time integrators, solvers, and other mathematical specs
####################################################################################

Common specs for all solvers and time integrators, used in PKs.


TimeIntegrator
=================

 Factory for creating TimestepController objects

A TimestepController object sets what size timestep to take.  This can be a
variety of things, from fixed timestep size, to adaptive based upon error
control, to adapter based upon simple nonlinear iteration counts.


* `"timestep controller type`" ``[string]`` Set the type.  One of the below types.
* `"timestep controller X parameters`" ``[list]`` List of parameters for a timestep controller of type X.

Available types include:

- TimestepControllerFixed_  (type `"fixed`"), a constant timestep
- TimestepControllerStandard_ (type `'standard`"), an adaptive timestep based upon nonlinear iterations
- TimestepControllerSmarter_ (type `'smarter`"), an adaptive timestep based upon nonlinear iterations with more control
- TimestepControllerAdaptive_ (type `"adaptive`"), an adaptive timestep based upon error control.
- TimestepControllerFromFile_ (type `"from file`"), uses a timestep history loaded from an HDF5 file.  (Usually only used for regression testing.)




TimestepControllerFixed
--------------------------
  Timestep controller providing constant timestep size.

``TimestepControllerFixed`` is a simple timestep control mechanism which sets
a constant timestep size.  Note that the actual timestep size is given by the
minimum of PK's initial timestep sizes.





TimestepControllerStandard
----------------------------
 Simple timestep control based upon previous iteration count.

``TimestepControllerStandard`` is a simple timestep control mechanism
which sets the next timestep based upon the previous timestep and how many
nonlinear iterations the previous timestep took to converge.

The timestep for step :math:`k+1`, :math:`\Delta t_{k+1}`, is given by:

- if :math:`N_k > N^{max}` then :math:`\Delta t_{k+1} = f_{reduction} * \Delta t_{k}`
- if :math:`N_k < N^{min}` then :math:`\Delta t_{k+1} = f_{increase} * \Delta t_{k}`
- otherwise :math:`\Delta t_{k+1} = \Delta t_{k}`

where :math:`\Delta t_{k}` is the previous timestep and :math:`N_k` is the number of 
nonlinear iterations required to solve step :math:`k`:.

* `"max iterations`" ``[int]`` :math:`N^{max}`, decrease the timestep if the previous step took more than this.
* `"min iterations`" ``[int]`` :math:`N^{min}`, increase the timestep if the previous step took less than this.
* `"time step reduction factor`" ``[double]`` :math:`f_reduction`, reduce the previous timestep by this multiple.
* `"time step increase factor`" ``[double]`` :math:`f_increase`, increase the previous timestep by this multiple.
* `"max time step`" ``[double]`` The max timestep size allowed.
* `"min time step`" ``[double]`` The min timestep size allowed.  If the step has failed and the new step is below this cutoff, the simulation fails.




TimestepControllerSmarter
----------------------------
  Slightly smarter timestep controller based upon a history of previous timesteps.

``TimestepControllerSmarter`` is based on ``TimestepControllerStandard``, but
also tries to be a bit smarter to avoid repeated increase/decrease loops where
the step size decreases, converges in few iterations, increases, but then
fails again.  It also tries to grow the step geometrically to more quickly
recover from tricky nonlinearities.

* `"max iterations`" ``[int]`` :math:`N^{max}`, decrease the timestep if the previous step took more than this.
* `"min iterations`" ``[int]`` :math:`N^{min}`, increase the timestep if the previous step took less than this.
* `"time step reduction factor`" ``[double]`` :math:`f_reduction`, reduce the previous timestep by this multiple.
* `"time step increase factor`" ``[double]`` :math:`f_increase`, increase the previous timestep by this multiple.  Note that this can be modified geometrically in the case of repeated successful steps.
* `"max time step increase factor`" ``[double]`` **10.** The max :math:`f_increase` will ever get.
* `"growth wait after fail`" ``[int]`` Wait at least this many timesteps before attempting to grow the timestep after a failed timestep.
* `"count before increasing increase factor`" ``[int]`` Require this many successive increasions before multiplying :math:`f_increase` by itself.





TimestepControllerFromFile
----------------------------
  Timestep controller which loads a timestep history from file.

``TimestepControllerFromFile`` loads a timestep history from a file, then
advances the step size with those values.  This is mostly used for testing
purposes, where we need to force the same timestep history as previous runs to
do regression testing.  Otherwise even machine roundoff can eventually alter
number of iterations enough to alter the timestep history, resulting in
solutions which are enough different to cause doubt over their correctness.

* `"file name`" ``[string]`` Path to hdf5 file containing timestep information.
* `"timestep header`" ``[string]`` Name of the dataset containing the history of timestep sizes.








Linear Solver Spec
===================

For each solver, a few parameters are used:

* `"iterative method`" ``[string]`` `"pcg`", `"gmres`", or `"nka`"

  defines which method to use.

* `"error tolerance`" ``[double]`` **1.e-6** is used in the convergence test.

* `"maximum number of iterations`" ``[int]`` **100** is used in the convergence test.

* `"convergence criteria`" ``[Array(string)]``  **{"relative rhs"}** specifies multiple convergence criteria. The list
  may include `"relative residual`", `"relative rhs`", and `"absolute residual`", and `"???? force once????`"

* `"size of Krylov space`" ``[int]`` is used in GMRES iterative method. The default value is 10.

.. code-block:: xml

     <ParameterList name="my solver">
       <Parameter name="iterative method" type="string" value="gmres"/>
       <Parameter name="error tolerance" type="double" value="1e-12"/>
       <Parameter name="maximum number of iterations" type="int" value="400"/>
       <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
       <Parameter name="size of Krylov space" type="int" value="10"/>

       <ParameterList name="VerboseObject">
         <Parameter name="Verbosity Level" type="string" value="high"/>
       </ParameterList>
     </ParameterList>


Preconditioner
===================

These can be used by a process kernel lists to define a preconditioner.  The only common parameter required by all lists is the type:

 * `"preconditioner type`" ``[string]`` **"identity"**, `"boomer amg`", `"trilinos ml`", `"block ilu`" ???
 * `"PC TYPE parameters`" ``[list]`` includes a list of parameters specific to the type of PC.

Example:

.. code-block:: xml

     <ParameterList name="my preconditioner">
       <Parameter name="type" type="string" value="trilinos ml"/>
        <ParameterList name="trilinos ml parameters"> ?????? check me!
            ... 
        </ParameterList>
     </ParameterList>


Hypre's Boomer AMG
-------------------
 PreconditionerBoomerAMG: HYPRE's multigrid preconditioner.
Internal parameters for Boomer AMG include

* `"tolerance`" ``[double]`` if is not zero, the preconditioner is dynamic 
  and approximate the inverse matrix with the prescribed tolerance (in
  the energy norm ???).

* `"smoother sweeps`" ``[int]`` **3** defines the number of smoothing loops. Default is 3.

* `"cycle applications`" ``[int]`` **5** defines the number of V-cycles.

* `"strong threshold`" ``[double]`` **0.5** defines the number of V-cycles. Default is 5.

* `"relaxation type`" ``[int]`` **6** defines the smoother to be used. Default is 6 
  which specifies a symmetric hybrid Gauss-Seidel / Jacobi hybrid method. TODO: add others!

* `"coarsen type`" ``[int]`` **0** defines the coarsening strategy to be used. Default is 0 
  which specifies a Falgout method. TODO: add others!

* `"max multigrid levels`" ``[int]`` optionally defined the maximum number of multigrid levels.

* `"use block indices`" ``[bool]`` **false** If true, uses the `"systems of
    PDEs`" code with blocks given by the SuperMap, or one per DoF per entity
    type.

* `"number of functions`" ``[int]`` **1** Any value > 1 tells Boomer AMG to
  use the `"systems of PDEs`" code with strided block type.  Note that, to use
  this approach, unknowns must be ordered with DoF fastest varying (i.e. not
  the native Epetra_MultiVector order).  By default, it uses the `"unknown`"
  approach in which each equation is coarsened and interpolated independently.
  
* `"nodal strength of connection norm`" ``[int]`` tells AMG to coarsen such
    that each variable has the same coarse grid - sometimes this is more
    "physical" for a particular problem. The value chosen here for nodal
    determines how strength of connection is determined between the
    coupled system.  I suggest setting nodal = 1, which uses a Frobenius
    norm.  This does NOT tell AMG to use nodal relaxation.
    Default is 0.

* `"verbosity`" ``[int]`` **0** prints a summary of run time settings and
  timing information to stdout.  `"1`" prints coarsening info, `"2`" prints
  smoothing info, and `"3`'" prints both.

Example:
  
.. code-block:: xml

  <ParameterList name="boomer amg parameters">
    <Parameter name="tolerance" type="double" value="0.0"/>
    <Parameter name="smoother sweeps" type="int" value="3"/>
    <Parameter name="cycle applications" type="int" value="5"/>
    <Parameter name="strong threshold" type="double" value="0.5"/>
    <Parameter name="coarsen type" type="int" value="0"/>
    <Parameter name="relaxation type" type="int" value="3"/>
    <Parameter name="verbosity" type="int" value="0"/>
    <Parameter name="number of functions" type="int" value="1"/>
  </ParameterList>




Trilinos ML
-------------------
 PreconditionerML: Trilinos ML multigrid.
Internal parameters of Trilinos ML includes

Example:

.. code-block:: xml

   <ParameterList name="ml parameters">
     <Parameter name="ML output" type="int" value="0"/>
     <Parameter name="aggregation: damping factor" type="double" value="1.33"/>
     <Parameter name="aggregation: nodes per aggregate" type="int" value="3"/>
     <Parameter name="aggregation: threshold" type="double" value="0.0"/>
     <Parameter name="aggregation: type" type="string" value="Uncoupled"/>
     <Parameter name="coarse: type" type="string" value="Amesos-KLU"/>
     <Parameter name="coarse: max size" type="int" value="128"/>
     <Parameter name="coarse: damping factor" type="double" value="1.0"/>
     <Parameter name="cycle applications" type="int" value="2"/>
     <Parameter name="eigen-analysis: iterations" type="int" value="10"/>
     <Parameter name="eigen-analysis: type" type="string" value="cg"/>
     <Parameter name="max levels" type="int" value="40"/>
     <Parameter name="prec type" type="string" value="MGW"/>
     <Parameter name="smoother: damping factor" type="double" value="1.0"/>
     <Parameter name="smoother: pre or post" type="string" value="both"/>
     <Parameter name="smoother: sweeps" type="int" value="2"/>
     <Parameter name="smoother: type" type="string" value="Gauss-Seidel"/>
   </ParameterList>

 


Block ILU
-------------------
 PreconditionerBlockILU:   Incomplete LU preconditioner.

The internal parameters for block ILU are as follows:

Example:

.. code-block:: xml

  <ParameterList name="block ilu parameters">
    <Parameter name="fact: relax value" type="double" value="1.0"/>
    <Parameter name="fact: absolute threshold" type="double" value="0.0"/>
    <Parameter name="fact: relative threshold" type="double" value="1.0"/>
    <Parameter name="fact: level-of-fill" type="int" value="0"/>
    <Parameter name="overlap" type="int" value="0"/>
    <Parameter name="schwarz: combine mode" type="string" value="Add"/>
    </ParameterList>
  </ParameterList>




Indentity
-------------------
The default, no PC applied.



NonlinearSolver
===================




Other Common Specs
##########################################

IOEvent
===================

 IOEvent: base time/timestep control determing when in time to do something.

The IOEvent is used for multiple objects that need to indicate simulation times or cycles on which to do something.

* `"cycles start period stop`" ``[Array(int)]`` **optional**

    The first entry is the start cycle, the second is the cycle
    period, and the third is the stop cycle or -1, in which case there
    is no stop cycle. A visualization dump is written at such
    cycles that satisfy cycle = start + n*period, for n=0,1,2,... and
    cycle < stop if stop != -1.0.

* `"cycles start period stop 0`" ``[Array(int)]`` **optional** 

    If multiple cycles start period stop parameters are needed, then use these
    parameters.  If one with 0 is found, then one with 1 is looked for, etc,
    until the Nth one is not found.

* `"cycles`" ``[Array(int)]``  **optional**
  
    An array of discrete cycles that at which a visualization dump is
    written.

* `"times start period stop`" ``[Array(double)]`` **optional** 

    The first entry is the start time, the second is the time period,
    and the third is the stop time or -1, in which case there is no
    stop time. A visualization dump is written at such times that
    satisfy time = start + n*period, for n=0,1,2,... and time < stop
    if stop != -1.0.

* `"times start period stop units`" ``string`` **s** 

    Units corresponding to this spec.  One of `"s`", `"d`", `"yr`", or `"yr 365`"
    
* `"times start period stop 0`" ``[Array(double)]`` **optional**

    If multiple start period stop parameters are needed, then use this these
    parameters with N=0,1,2.  If one with 0 is found, then one with 1 is
    looked for, etc, until the Nth one is not found.

* `"times start period stop 0 units`" ``string`` **s** 

    Units corresponding to this spec.  One of `"s`", `"d`", `"yr`", or `"yr 365`"
    See above for continued integer listings.

* `"times`" ``[Array(double)]`` **optional** 

    An array of discrete times that at which a visualization dump
    shall be written.

* `"times units`" ``string`` **s** 

    Units corresponding to this spec.  One of `"s`", `"d`", `"yr`", or `"yr 365`"
    
 


VerboseObject
===================

 VerboseObject: a controller for writing log files on multiple cores with varying verbosity.

This allows control of log-file verbosity for a wide variety of objects
and physics.

* `"verbosity level`" ``[string]`` **GLOBAL_VERBOSITY**, `"low`", `"medium`", `"high`", `"extreme`"

   The default is set by the global verbosity spec, (fix me!)  Typically,
   `"low`" prints out minimal information, `"medium`" prints out errors and
   overall high level information, `"high`" prints out basic debugging, and
   `"extreme`" prints out local debugging information.

Note: while there are other options, users should typically not need them.
Instead, developers can use them to control output.
   
Example:

.. code-block:: xml

  <ParameterList name="verbose object">
    <Parameter name="verbosity level" type="string" value="medium"/>
    <Parameter name="name" type="string" value="my header"/>
    <Parameter name="hide line prefix" type="bool" value="false"/>
    <Parameter name="write on rank" type="int" value="0"/>
  </ParameterList>



   

Function
===================

 Function: base class for all functions of space and time.
Analytic, algabraic functions of space and time are used for a variety of
purposes, including boundary conditions, initial conditions, and independent
variables.

For initial conditions, functions are prescribed of space only, i.e.

:math:`u = f(x,y,z)`

For boundary conditions and independent variables, functions are also a
function of time:

:math:`u = f(t,x,y,z)`

``[function-spec]``

ONE OF:
* `"function: constant`" ``[constant-function-spec]``
OR:
* `"function: tabular`" ``[tabular-function-spec]``
OR:
* `"function: smooth step`" ``[smooth-step-function-spec]``
OR:
* `"function: polynomial`" ``[polynomial-function-spec]``
OR:
* `"function: monomial`" ``[monomial-function-spec]``
OR:
* `"function: linear`" ``[linear-function-spec]``
OR:
* `"function: separable`" ``[separable-function-spec]``
OR:
* `"function: additive`" ``[additive-function-spec]``
OR:
* `"function: multiplicative`" ``[multiplicative-function-spec]``
OR:
* `"function: composition`" ``[composition-function-spec]``
OR:
* `"function: static head`" ``[static-head-function-spec]``
OR:
* `"function: standard math`" ``[standard-math-function-spec]``
OR:
* `"function: bilinear`" ``[bilinear-function-spec]``
OR:
* `"function: distance`" ``[distance-function-spec]``
#OR:
#* `"function: squared distance`" ``[squared-distance-function-spec]``
END



It is straightforward to add new functions as needed.

Constant Function
-------------------------
** DOC GENERATION ERROR: file not found ' ConstantFunction ' **  

Tabular Function
-------------------------
** DOC GENERATION ERROR: file not found ' TabularFunction ' **

Smooth step Function
-------------------------
** DOC GENERATION ERROR: file not found ' SmoothStepFunction ' **

Polynomial Function
-------------------------
** DOC GENERATION ERROR: file not found ' PolynomialFunction ' **  

Multi-variable linear Function
------------------------------
** DOC GENERATION ERROR: file not found ' LinearFunction ' **  

Separable Function
------------------
** DOC GENERATION ERROR: file not found ' SeparableFunction ' **

Additive Function
------------------
** DOC GENERATION ERROR: file not found ' AdditiveFunction ' **

Multiplicative Function
--------------------------
** DOC GENERATION ERROR: file not found ' MultiplicativeFunction ' **

Composition Function
--------------------------
** DOC GENERATION ERROR: file not found ' CompositionFunction ' **

Piecewise Bilinear Function
---------------------------
** DOC GENERATION ERROR: file not found ' BilinearFunction ' **

Distance Function
-------------------
** DOC GENERATION ERROR: file not found ' DistanceFunction ' **

Monomial Function
-------------------
** DOC GENERATION ERROR: file not found ' MonomialFunction ' **

Standard Math Function
-------------------------
** DOC GENERATION ERROR: file not found ' StandardMathFunction ' **



Operator
===================

 Operator represents a linear map, and typically encapsulates a discretization.
``Operator`` represents a map from linear space X to linear space Y.  Typically,
this map is a linear map, and encapsulates much of the discretization involved
in moving from continuous to discrete equations. The spaces X and Y are described 
by CompositeVectors (CV). A few maps X->Y are supported. 

An ``Operator`` provides an interface for applying both the forward and inverse
linear map (assuming the map is invertible).

Typically the ``Operator`` is never seen by the user; instead the user provides
input information for helper classes based on the continuous mathematical
operator and the desired discretization.  These helpers build the needed
``Operator``, which may include information from multiple helpers (i.e. in the
case of Jacobian Operators for a PDE).

However, one option may be provided by the user, which is related to dealing
with nearly singular operators:

* `"diagonal shift`" ``[double]`` **0.0** Adds a scalar shift to the diagonal
  of the ``Operator``, which can be useful if the ``Operator`` is singular or
  near-singular.



OperatorAccumulation
-------------------------

** DOC GENERATION ERROR: file not found ' OperatorAccumulation ' **

OperatorDiffusion
------------------
** DOC GENERATION ERROR: file not found ' OperatorDiffusionFactory ' **

** DOC GENERATION ERROR: file not found ' OperatorDiffusionMFD ' **

** DOC GENERATION ERROR: file not found ' OperatorDiffusionMFDwithGravity ' **

** DOC GENERATION ERROR: file not found ' OperatorDiffusion ' **


OperatorAdvection
-------------------------

** DOC GENERATION ERROR: file not found ' OperatorAdvection ' **



