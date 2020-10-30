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
/// @param time_step - the time step for which the mesh should be read
/// @param partition - the partition to read
Mesh read_fea_io_from_meta(const std::string& meta_filepath,
                           const int time_step       = 0,
                           const PartitionIndex partition = 0);

/// Serially writes a mesh in the JSON/binary format, which supports distributed meshes.
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
/// @param element_to_component_map (optional) - a map from the local element indices to the
///        component to which the element belongs.  When visualizing the results with the fea_io
///        files, the component ID's basically control whether quantities should be continuous
///        across element boundaries or not.  It is typically useful to assign each material a
///        different component if you want to avoid averaging of quantities across element
///        boundaries.  If omitted, every element gets assigned 0.
/// @param component_names (optional) - only used if element_to_component_map is included.  Assigns
///        a name to each component, which may make post-processing easier by keeping track of which
///        material corresponds to a given component ID.
void write_fea_io(Mesh& mesh,
                  const std::string& base_directory,
                  const std::string& mesh_name,
                  const std::string meta_file_name                             = "",
                  const std::vector<ComponentIndex>& element_to_component_map  = {},
                  const std::map<ComponentIndex, std::string>& component_names = {});

/// Writes a distributed mesh (not just partitioned but broken up across ranks) in parallel to the
/// JSON/binary format, which supports distributed meshes.
/// @important This method should be used multiple ranks each hold a part of the entire mesh.  See
///            write_fea_io() if one rank holds the entire mesh (even if the mesh is partitioned).
/// @param mesh - part of the mesh that this rank holds.  Ghost elements will not be generated, as
///        it is assumed that the partition already has any relevant ghost information inside of it.
/// @param base_directory - directory to write the resulting files (there will be multiple).
///        This can be a relative (to cwd) or absolute path to the directory.
/// @param mesh_name - name of the mesh files without extensions. (.m.json and .m.bin will be
///        appended)
/// @param meta_file_name (optional) - name of the meta file without an extention.  (.x.json will be
///        appended)  If this is an empty string, no meta file will be written.
/// @param element_to_component_map (optional) - a map from the local element indices (this should
/// only include the elements in the mesh passed into this function) to the
///        component to which the element belongs.  When visualizing the results with the fea_io
///        files, the component ID's basically control whether quantities should be continuous
///        across element boundaries or not.  It is typically useful to assign each material a
///        different component if you want to avoid averaging of quantities across element
///        boundaries.  If omitted, every element gets assigned 0.
/// @param component_names (optional) - only used if element_to_component_map is included.  Assigns
///        a name to each component, which may make post-processing easier by keeping track of which
///        material corresponds to a given component ID.
void write_distributed_fea_io(Mesh& mesh,
                              const std::string& base_directory,
                              const std::string& mesh_name,
                              const std::string meta_file_name                             = "",
                              const std::vector<ComponentIndex>& element_to_component_map  = {},
                              const std::map<ComponentIndex, std::string>& component_names = {});

#endif
}