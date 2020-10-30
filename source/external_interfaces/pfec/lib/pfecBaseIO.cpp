#include "pfecBaseIO.h"

PFECNativeInterface::BaseIO::BaseIO()
{
    num_local_parts_ = 0;
    header_end_ = 0;
}

PFECNativeInterface::BaseIO::~BaseIO()
{
    close();
}

void PFECNativeInterface::BaseIO::close()
{
    num_local_parts_ = 0;
    local_part_ids_.clear();
    part_id_map_.clear();
    header_end_ = 0;
}
