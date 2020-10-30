#include "pfecMetaWriter.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"

using namespace rapidjson;

PFECNativeInterface::MetaWriter::MetaWriter()
	: MetaReader()
{
	time_step_started_ = false;
	curr_time_meta_index_ = 0;
}

PFECNativeInterface::MetaWriter::MetaWriter(std::string filePath, std::string analysisName)
	: MetaWriter()
{
	open(filePath, analysisName);
}

PFECNativeInterface::MetaWriter::~MetaWriter()
{
    if (time_step_started_)
		write_json_meta();
}

void PFECNativeInterface::MetaWriter::open(std::string filePath, std::string analysisName)
{
	meta_filepath_ = filePath;
	meta_file_dir_ = get_directory(filePath);
}

void PFECNativeInterface::MetaWriter::start_new_time_step(int timeStepID, double timeVal)
{
	if (time_step_started_)
		write_json_meta();

	// If the time ID does not already exist, then create a meta for it
	if (time_step_id_map_.find(timeStepID) == time_step_id_map_.end()) {
		TimeStepMeta newStep;
		newStep.time_step_id_ = timeStepID;
		newStep.time_value_ = timeVal;
		time_step_metas_.push_back(newStep);
		time_step_id_map_[timeStepID] = time_step_metas_.size() - 1;
	}
	else if (timeVal != get_time_step_meta(timeStepID).time_value_)
		throw _MSG_EXCEPT_C(ERR_META_FILE_ALREADY_HAS_ENTRY);

	time_step_started_ = true;
    curr_time_meta_index_ = static_cast<PFEC_INT>(time_step_id_map_[timeStepID]);
}

std::unique_ptr<PFECNativeInterface::MeshWriter> PFECNativeInterface::MetaWriter::get_mesh_writer(std::string fileRelPath)
{
    // Check that there it hasn't already been written for this time step
    if (get_time_step_meta(curr_time_meta_index_).mesh_file_path_ != "")
        throw _MSG_EXCEPT_C(ERR_META_FILE_ALREADY_HAS_ENTRY);

	get_time_step_meta(curr_time_meta_index_).mesh_file_path_ = fileRelPath;
    std::unique_ptr<PFECNativeInterface::MeshWriter> mwPtr = std::unique_ptr<MeshWriter>(new MeshWriter());
    mwPtr->open(meta_file_dir_ + "/" + get_time_step_meta(curr_time_meta_index_).mesh_file_path_);
    return std::move(mwPtr);
}

std::unique_ptr<PFECNativeInterface::DataWriter> PFECNativeInterface::MetaWriter::get_orient_writer(std::string fileRelPath)
{
    // Check that there it hasn't already been written for this time step
    if (get_time_step_meta(curr_time_meta_index_).orient_file_path_ != "")
        throw _MSG_EXCEPT_C(ERR_META_FILE_ALREADY_HAS_ENTRY);

	get_time_step_meta(curr_time_meta_index_).orient_file_path_ = fileRelPath;
    std::unique_ptr<PFECNativeInterface::DataWriter> dwPtr = std::unique_ptr<DataWriter>(new DataWriter(meta_file_dir_ + "/" + fileRelPath));
    return std::move(dwPtr);
}

std::unique_ptr<PFECNativeInterface::DataWriter> PFECNativeInterface::MetaWriter::get_data_writer(std::string fileRelPath)
{
    // Check that there it hasn't already been written for this time step
    if (get_time_step_meta(curr_time_meta_index_).data_file_path_ != "")
        throw _MSG_EXCEPT_C(ERR_META_FILE_ALREADY_HAS_ENTRY);

	get_time_step_meta(curr_time_meta_index_).data_file_path_ = fileRelPath;
    std::unique_ptr<PFECNativeInterface::DataWriter> dwPtr = std::unique_ptr<DataWriter>(new DataWriter(meta_file_dir_ + "/" + fileRelPath));
    return std::move(dwPtr);
}

std::unique_ptr<PFECNativeInterface::GroupWriter> PFECNativeInterface::MetaWriter::get_group_writer(std::string fileRelPath)
{
	get_time_step_meta(curr_time_meta_index_).group_file_paths_.push_back(fileRelPath);
    std::unique_ptr<PFECNativeInterface::GroupWriter> gwPtr = std::unique_ptr<GroupWriter>(new GroupWriter());
    gwPtr->open(meta_file_dir_ + "/" + fileRelPath);
    return std::move(gwPtr);
}

void PFECNativeInterface::MetaWriter::prepend_dir_to_paths(const std::string dir)
{
	for (auto& tsMeta : time_step_metas_) {
		if (tsMeta.mesh_file_path_ != "") tsMeta.mesh_file_path_ = dir + tsMeta.mesh_file_path_;
		if (tsMeta.orient_file_path_ != "") tsMeta.orient_file_path_ = dir + tsMeta.orient_file_path_;
		if (tsMeta.data_file_path_ != "") tsMeta.data_file_path_ = dir + tsMeta.data_file_path_;
		for (auto& groupPath : tsMeta.group_file_paths_) groupPath = dir + groupPath;
	}
}

void PFECNativeInterface::MetaWriter::clear_existing_data_references()
{
	for (auto& tsMeta : time_step_metas_) {
		tsMeta.data_file_path_ = "";
	}
}

void PFECNativeInterface::MetaWriter::write_json_meta()
{
	// Create the JSON structure
	Document analysisJSON(kObjectType);  // Default template parameter uses UTF8 and MemoryPoolAllocator.
	Document::AllocatorType& allocator = analysisJSON.GetAllocator();

	Value analysisItem(kObjectType);
	{
		Value analysisName;
        analysisName.SetString(analysis_name_.c_str(), static_cast<rapidjson::SizeType>(analysis_name_.size()), allocator);
		analysisItem.AddMember("name", analysisName, allocator);

		// Expected time step format:
		// {
		//		"id" : 0,
		//		"time" : 0.0,
		//		"mesh" : "Test.pfecm",
		//		"orientation" : "Test.pfeco",
		//		"data" : {
		//			"file": "Test.pfecd",
		//			"arrays" : [
		//				{"type" : "Nodal", "name" : "Displacement"}
		//			]
		//		},
		//		"groups" : [ "Test.pfecg" ]
		// }
		Value timeStepArray(kArrayType);
		auto dataTypeNameMap = get_data_type_name_map();
		for (auto& timeStepMeta : time_step_metas_) {
			Value timeStepObject(kObjectType);

			timeStepObject.AddMember("id", Value(timeStepMeta.time_step_id_), allocator);
			timeStepObject.AddMember("time", Value(timeStepMeta.time_value_), allocator);

			// Now add paths if they exist
			if (timeStepMeta.mesh_file_path_ != "") {
				Value filepath;
                filepath.SetString(timeStepMeta.mesh_file_path_.c_str(), static_cast<rapidjson::SizeType>(timeStepMeta.mesh_file_path_.size()), allocator);
				timeStepObject.AddMember("mesh", filepath, allocator);
			}
			if (timeStepMeta.orient_file_path_ != "") {
				Value filepath;
                filepath.SetString(timeStepMeta.orient_file_path_.c_str(), static_cast<rapidjson::SizeType>(timeStepMeta.orient_file_path_.size()), allocator);
				timeStepObject.AddMember("orientation", filepath, allocator);
			}
			if (timeStepMeta.data_file_path_ != "") {
				Value dataObject(kObjectType);

				Value filepath;
                filepath.SetString(timeStepMeta.data_file_path_.c_str(), static_cast<rapidjson::SizeType>(timeStepMeta.data_file_path_.size()), allocator);
				dataObject.AddMember("file", filepath, allocator);

				// Extract what data types are in the data file
				timeStepMeta.extract_data_array_info(meta_file_dir_);
				Value dataArrays(kArrayType);
				for (auto i = 0; i < timeStepMeta.data_types_.size(); i++) {
					Value dataArrayEntry(kObjectType);
					{
						Value dataType;
						const std::string dataTypeStr = dataTypeNameMap[timeStepMeta.data_types_[i]];
                        dataType.SetString(dataTypeStr.c_str(), static_cast<rapidjson::SizeType>(dataTypeStr.size()), allocator);
						dataArrayEntry.AddMember("type", dataType, allocator);
						Value dataName;
                        dataName.SetString(timeStepMeta.data_names_[i].c_str(), static_cast<rapidjson::SizeType>(timeStepMeta.data_names_[i].size()), allocator);
						dataArrayEntry.AddMember("name", dataName, allocator);
					}
					dataArrays.PushBack(dataArrayEntry, allocator);
				}
				dataObject.AddMember("arrays", dataArrays, allocator);

				timeStepObject.AddMember("data", dataObject, allocator);
			}
			if (timeStepMeta.group_file_paths_.size() > 0) {
				Value groupArray(kArrayType);
				for (auto& groupFilePath : timeStepMeta.group_file_paths_) {
					Value filepath;
                    filepath.SetString(groupFilePath.c_str(), static_cast<rapidjson::SizeType>(groupFilePath.size()), allocator);
					groupArray.PushBack(filepath, allocator);
				}
				timeStepObject.AddMember("groups", groupArray, allocator);
			}

			timeStepArray.PushBack(timeStepObject, allocator);
		}
		analysisItem.AddMember("timesteps", timeStepArray, allocator);
	}
	analysisJSON.AddMember("analysis", analysisItem, allocator);

	// Stringify the JSON
	StringBuffer buffer;
	PrettyWriter<StringBuffer> writer(buffer);
	analysisJSON.Accept(writer);

	// Write the string
	std::ofstream metaFile(meta_filepath_);
	metaFile << buffer.GetString();
	metaFile.close();
}

PFECNativeInterface::TimeStepMeta & PFECNativeInterface::MetaWriter::get_time_step_meta(const int metaTimeIndex)
{
	return time_step_metas_[metaTimeIndex];
}
