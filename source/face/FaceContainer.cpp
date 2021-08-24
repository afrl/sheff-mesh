#include "mesh/face/FaceContainer.h"

#include "mesh/face/FaceAccessor.h"
#include "mesh/face/FaceIterator.h"

namespace mesh
{
FaceContainer::FaceContainer(ElementConstAccessor element)
  : element_(element)
{
    auto& face_infos = get_face_info(element.type());
    types_.reserve(face_infos.size());
    connectivities_.reserve(face_infos.size());
    connectivity_offsets_.reserve(face_infos.size() + 1);
    connectivity_offsets_ = { 0 };
    for (auto& face_info : face_infos)
    {
        types_.push_back(face_info.type);
        for (auto& i : face_info.connectivity)
        {
            connectivities_.push_back(i);
        }
        connectivity_offsets_.push_back(static_cast<ConnIndex>(connectivities_.size()));
    }
}

FaceConstIterator FaceContainer::begin() const
{
    return FaceConstIterator(*this);
}
FaceConstIterator FaceContainer::end() const
{
    return FaceConstIterator(*this, size());
}

FaceConstAccessor FaceContainer::operator[](ConnIndex i) const
{
    return FaceConstAccessor(*this, i);
}

FaceConstAccessor FaceContainer::back() const
{
    return FaceConstAccessor(*this, size() - 1);
}

ConnIndex FaceContainer::size() const
{
    return static_cast<ConnIndex>(types_.size());
}
}
