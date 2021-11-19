%module constitutive

%{
#include "utility/SingletonNamedObjectContainer.h"

#include "simple_math/Defines.h"
using namespace simple_math;

#include "constitutive/isv/ISVSet.h"
#include "constitutive/isv/ISVSetContainer.h"
#include "constitutive/enrichment/EnrichmentProperties.h"
#include "constitutive/enrichment/PointwiseEnrichmentModel.h"
#include "constitutive/enrichment/PointwiseEnrichmentModelInputs.h"
#include "constitutive/enrichment/PointwiseEnrichmentModelWorkspace.h"
#include "constitutive/enrichment/deformation/DeformationPointwiseEnrichmentModel.h"
#include "constitutive/enrichment/deformation/turon_06_cz/InfStrainTuron06EnrichmentModel.h"
#include "constitutive/enrichment/deformation/turon_06_cz/Turon06EnrichmentProperties.h"
#include "constitutive/enrichment/deformation/turon_06_cz/Turon06ISVSet.h"
using namespace constitutive;
%}

%include <typemaps.i>
%include <std_string.i>
%include <std_array.i>
%include <stdint.i>

%template(ArrayD3) std::array<double, 3>;

// Eigen matrices into Numpy arrays.
%include "eigen.i"

%eigen_typemaps(simple_math::Vector3D)
%eigen_typemaps(simple_math::Matrix3D)
%eigen_typemaps(simple_math::VoigtVector)
%eigen_typemaps(simple_math::VoigtMatrix)
%eigen_typemaps(simple_math::GeneralVector)
%eigen_typemaps(simple_math::GeneralMatrix)

%ignore *::clone;

%include "utility/Bitset.h"
%include "utility/SingletonNamedObjectContainer.h"
%include "utility/Slice.h"

%template(SliceDouble) utility::Slice<double>;

%include "simple_math/Defines.h"

%include "constitutive/isv/ISVSet.h"
%include "constitutive/isv/ISVSetContainer.h"

%extend constitutive::ISVSetContainer {
  %template(get_set)    get_set<constitutive::ISVSet>;
};

%include "constitutive/enrichment/PointwiseEnrichmentModelInputs.h"

%template(PointwiseEnrichmentModelInputsBitset) utility::bitset<constitutive::PointwiseEnrichmentModelInputs::Flags>;

%include "constitutive/enrichment/EnrichmentProperties.h"

%template(EnrichmentPropertiesCollection) utility::SingletonNamedObjectContainer<constitutive::EnrichmentProperties>;
%template(EnrichmentPropertiesCollectionRef) utility::SingletonNamedObjectContainerRef<constitutive::EnrichmentProperties>;

%include "constitutive/enrichment/PointwiseEnrichmentModelWorkspace.h"
%include "constitutive/enrichment/PointwiseEnrichmentModel.h"
%include "constitutive/enrichment/deformation/DeformationPointwiseEnrichmentModel.h"

%include "constitutive/enrichment/deformation/turon_06_cz/InfStrainTuron06EnrichmentModel.h"
%include "constitutive/enrichment/deformation/turon_06_cz/Turon06EnrichmentProperties.h"
%include "constitutive/enrichment/deformation/turon_06_cz/Turon06ISVSet.h"
