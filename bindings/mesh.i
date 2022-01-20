%module mesh

%{
#include <iostream>

#include "simple_math/Defines.h"
using namespace simple_math;

#include "mesh/Mesh.h"
using namespace mesh;
%}

%include <typemaps.i>
%include <attribute.i>
%include <std_string.i>
%include <std_array.i>
%include <stdint.i>

%template(ArrayD3) std::array<double, 3>;

// Eigen matrices into Numpy arrays.
%include "eigen.i"

%include "simple_math/Defines.h"
//%eigen_typemaps(simple_math::Quaternion)
%eigen_typemaps(simple_math::Vector3D)
%eigen_typemaps(simple_math::Matrix3D)
%eigen_typemaps(simple_math::VoigtVector)
%eigen_typemaps(simple_math::VoigtMatrix)
%eigen_typemaps(simple_math::GeneralVector)
%eigen_typemaps(simple_math::DynamicDoubleMatrix)
%eigen_typemaps(simple_math::DynamicUInt32Matrix)

%include "utility/Slice.h"
%template(SliceDouble) utility::Slice<double>;


%include "mesh/Defines.h"

%include "mesh/common/EdgeTypes.h"
%include "mesh/common/FaceTypes.h"
%include "mesh/common/ElementTypes.h"
%include "mesh/common/QuaternionContainer.h"

//////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper for NodeAccessor

%include "mesh/node/NodeAccessor.h"

%ignore mesh::NodeAccessorBase::coord;
%extend mesh::NodeAccessorBase
{
std::string __str__()
{
  return self->str();
}
std::array<double, 3>& get_point()
{
  return self->point();
}
void set_point(const std::array<double, 3>& p)
{
  self->point() = p;
}
mesh::GlobalIndex get_global_index()
{
  return self->global_index();
}
void set_global_index(mesh::GlobalIndex i)
{
  self->global_index() = i;
}
mesh::PartitionIndex get_parent_partition()
{
  return self->parent_partition();
}
void set_parent_partition(mesh::PartitionIndex i)
{
  self->parent_partition() = i;
}
};

//%attribute(mesh::NodeAccessorBase, mesh::GlobalIndex, g, global_index);

%template(NodeAccessor) mesh::NodeAccessorBase<mesh::NodeContainer>;

// End NodeAccessor
//////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper for NodeIterator

%include "mesh/node/NodeIterator.h"

%extend mesh::NodeIteratorBase {
mesh::NodeIteratorBase* __iter__()
{
  return self;
}
NodeAccessor __next__()
{
  if (self->is_end())
  {
    throw swig::stop_iteration();
  }
  else
  {
    NodeAccessor n = self->get();
    ++(*self);
    return n;
  }
}
};
%exception mesh::NodeIteratorBase::__next__
{
  try
  {
    $action // calls %extend function next() below
  }
  catch (swig::stop_iteration)
  {
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return nullptr;
  }
};
%template(NodeIterator) mesh::NodeIteratorBase<mesh::NodeAccessorBase<mesh::NodeContainer>>;

// End NodeAccessor
//////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper for NodeContainer

%include "mesh/node/NodeContainer.h"

%extend mesh::NodeContainer
{
NodeIterator __iter__()
{
  NodeIterator iter(self->begin());
  return iter;
}
NodeAccessor __getitem__(int i)
{
  return (*self)[i];
}
};

%attribute(mesh::NodeContainer, mesh::LocalIndex, s, size);

// End NodeAccessor
//////////////////////////////////////////////////////////////////////////////////////////////

%include "mesh/node/NodeSubrangeIterator.h"
%include "mesh/node/NodeGroup.h"

%include "mesh/element/ElementAccessor.h"
%include "mesh/element/ElementIterator.h"
%include "mesh/element/ElementSubrangeIterator.h"
%include "mesh/element/ElementGroup.h"
%include "mesh/element/ElementContainer.h"

%include "mesh/Mesh.h"



