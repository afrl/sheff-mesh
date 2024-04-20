%module mesh

%{
#include <iostream>

#include "simple_math/Defines.h"
using namespace simple_math;

#include "mesh/Mesh.h"
#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/Groups.h"
#include "mesh/algorithms/Extrusion.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/LCSAssignments.h"
#include "mesh/algorithms/Manipulation.h"
#include "mesh/algorithms/Refinement.h"
#include "mesh/algorithms/Transformation.h"
using namespace mesh;
%}

%include <typemaps.i>
%include <attribute.i>
%include <stdint.i>
%include <std_array.i>
%include <std_pair.i>
%include <std_string.i>
%include <std_vector.i>

%template(ArrayDouble3) std::array<double, 3>;
%template(ArrayInt3) std::array<int, 3>;
%template(ArrayUInt3) std::array<uint32_t, 3>;

%template(PairUIntUInt) std::pair<uint32_t, uint32_t>;

%template(VectorUInt) std::vector<uint32_t>;
%template(VectorDouble) std::vector<double>;
%template(VectorPairUIntUInt) std::vector<std::pair<uint32_t, uint32_t>>;

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

// Ignore the build in reference accessor methods, since they will cause problems.
// Instead, make it more Python by exposing them via attributes.
%ignore mesh::NodeAccessorBase::local_index() const;
%ignore mesh::NodeAccessorBase::global_index() const;
%ignore mesh::NodeAccessorBase::x() const;
%ignore mesh::NodeAccessorBase::y() const;
%ignore mesh::NodeAccessorBase::z() const;
%ignore mesh::NodeAccessorBase::point() const;
%ignore mesh::NodeAccessorBase::point(int8_t) const;
%ignore mesh::NodeAccessorBase::coord() const;
%ignore mesh::NodeAccessorBase::parent_partition() const;

%extend mesh::NodeAccessorBase
{
std::string __str__()
{
  return self->str();
}
};

%template(NodeAccessor) mesh::NodeAccessorBase<mesh::NodeContainer>;

%attribute(mesh::NodeAccessorBase<mesh::NodeContainer>,
           mesh::GlobalIndex,
           local_index,
           local_index);
%attributeref(mesh::NodeAccessorBase<mesh::NodeContainer>,
              mesh::GlobalIndex,
              global_index);
%attributeref(mesh::NodeAccessorBase<mesh::NodeContainer>,
              double,
              x);
%attributeref(mesh::NodeAccessorBase<mesh::NodeContainer>,
              double,
              y);
%attributeref(mesh::NodeAccessorBase<mesh::NodeContainer>,
              double,
              z);
%attribute2ref(mesh::NodeAccessorBase<mesh::NodeContainer>,
              %arg(std::array<double, 3>),
              point);
%attributeref(mesh::NodeAccessorBase<mesh::NodeContainer>,
              uint16_t,
              parent_partition);

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

// End NodeIterator
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
mesh::LocalIndex __len__()
{
    return self->size();
}
};

%attribute(mesh::NodeContainer, mesh::LocalIndex, s, size);

// End NodeContainer
//////////////////////////////////////////////////////////////////////////////////////////////

%include "mesh/node/NodeSubrangeIterator.h"
%include "mesh/node/NodeGroup.h"

//////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper for ElementAccessor

%include "mesh/element/ElementAccessor.h"

// Ignore the build in reference accessor methods, since they will cause problems.
// Instead, make it more Python by exposing them via attributes.
%ignore mesh::ElementAccessorBase::local_index() const;
%ignore mesh::ElementAccessorBase::global_index() const;
%ignore mesh::ElementAccessorBase::type() const;
%ignore mesh::ElementAccessorBase::connectivity() const;
// node/nodes?
%ignore mesh::ElementAccessorBase::parent_partition() const;
%ignore mesh::ElementAccessorBase::component_id() const;
%ignore mesh::ElementAccessorBase::local_coord_sys() const;
%ignore mesh::ElementAccessorBase::node_container() const;

%template(ElementAccessor) mesh::ElementAccessorBase<mesh::NodeContainer, mesh::ElementContainer>;

%extend mesh::ElementAccessor
{
std::string __str__()
{
  return self->str();
}
};

%attribute(%arg(mesh::ElementAccessorBase<mesh::NodeContainer, mesh::ElementContainer>),
           mesh::GlobalIndex,
           local_index,
           local_index);
%attributeref(%arg(mesh::ElementAccessorBase<mesh::NodeContainer, mesh::ElementContainer>),
              mesh::GlobalIndex,
              global_index);
%attributeref(%arg(mesh::ElementAccessorBase<mesh::NodeContainer, mesh::ElementContainer>),
              uint16_t,
              parent_partition);
%attributeref(%arg(mesh::ElementAccessorBase<mesh::NodeContainer, mesh::ElementContainer>),
              uint16_t,
              component_id);
%attributeref2(%arg(mesh::ElementAccessorBase<mesh::NodeContainer, mesh::ElementContainer>),
           %arg(utility::Slice<uint32_t, uint16_t>),
           connectivity,
           connectivity);

// End ElementAccessor
//////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper for ElementIterator

%include "mesh/element/ElementIterator.h"

%extend mesh::ElementIteratorBase {
mesh::ElementIteratorBase* __iter__()
{
  return self;
}
ElementAccessor __next__()
{
  if (self->is_end())
  {
    throw swig::stop_iteration();
  }
  else
  {
    ElementAccessor n = self->get();
    ++(*self);
    return n;
  }
}
};
%exception mesh::ElementIteratorBase::__next__
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
%template(ElementIterator) mesh::ElementIteratorBase<mesh::ElementAccessorBase<mesh::NodeContainer, mesh::ElementContainer>>;

// End ElementIterator
//////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper for ElementContainer

%include "mesh/element/ElementContainer.h"

%extend mesh::ElementContainer
{
ElementIterator __iter__()
{
  ElementIterator iter(self->begin());
  return iter;
}
ElementAccessor __getitem__(int i)
{
  return (*self)[i];
}
};

%attribute(mesh::ElementContainer, mesh::LocalIndex, s, size);

// End ElementContainer
//////////////////////////////////////////////////////////////////////////////////////////////


%include "mesh/element/ElementSubrangeIterator.h"
%include "mesh/element/ElementGroup.h"

%include "mesh/Mesh.h"

%include "mesh/algorithms/Creation.h"


%ignore add_to_group(Mesh&, const std::string, NodeSelection);
%ignore add_to_group(Mesh&, const std::string, ElementSelection);
%include "mesh/algorithms/Groups.h"

%include "mesh/algorithms/Extrusion.h"
%include "mesh/algorithms/LCSAssignments.h"
%include "mesh/algorithms/Manipulation.h"
%include "mesh/algorithms/Refinement.h"
%include "mesh/algorithms/Transformation.h"

%include "mesh/external_interfaces/fea_io/IO.h"
%include "mesh/external_interfaces/vtk/IO.h"
