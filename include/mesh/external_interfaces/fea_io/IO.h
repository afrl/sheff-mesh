#pragma once

#include "mesh/Mesh.h"
#include "mesh/common/ElementTypes.h"

#include <string>
#include <vector>

namespace mesh
{
#ifdef USE_FEAIO

/// Reads a fea_io JSON/binary files and populates a mesh.
/// @param file_stem - stem for the relative (to cwd) or absolute path (this is basically the file
///        path excluding the .m.json or .m.bin extension)
/// @param partition - the partition to read
Mesh read_fea_io(const std::string& file_stem, const PartitionIndex partition = 0);

/// Reads a fea_io JSON/binary files given the relevant meta information and populates a mesh.  This
/// avoids having to know the path of the particular mesh file and travers.
/// @param meta_filepath - relative (to cwd) or absolute path to the meta file (include the
///        .x.json extension)
/// @param partition - the partition to read
Mesh read_fea_io_from_meta(const std::string& meta_filepath,
                           const PartitionIndex partition = 0);

/// @brief Writes a mesh in the JSON/binary format using the fea_io library.
/// 
/// This works for centralized/distributed meshes and partitioned/unpartitioned meshes.  However,
/// there is the assumption that if it is distributed that each rank possesses a single partition.
/// 
/// @important This method should be used when a single rank holds the entire mesh.  See
///            write_distributed_fea_io() if the mesh is distributed across several ranks.
/// @param mesh - mesh to write to the file. If the mesh is partitioned, ghost elements will
///        be included for each partition.
/// @param base_directory - directory to write the resulting files (there will be multiple).
///        This can be a relative (to cwd) or absolute path to the directory.
/// @param mesh_name - name of the mesh files without extensions. (.m.json and .m.bin will be
///        appended)
/// @param meta_file_name (optional) - name of the meta file without an extention.  (.x.json will be
///        appended)  If this is an empty string, no meta file will be written.
/// @param element_category_to_write (optional) - indicates which type of elements to write to a
///        fea_io file.  During pre-processing, it is generally helpful to write any type (default).
///        However, during post-processing, it is helpful to separate solid and interfacial meshes
///        since the data for the respective types will be fundamentally different.
/// @param component_names (optional) - only used if element_to_component_map is included.  Assigns
///        a name to each component, which may make post-processing easier by keeping track of which
///        material corresponds to a given component ID.
void write_fea_io(Mesh& mesh,
                  const std::string& base_directory,
                  const std::string& mesh_name,
                  const std::string meta_file_name                     = "");

#endif
}