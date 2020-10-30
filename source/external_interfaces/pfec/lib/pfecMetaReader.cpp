#include "pfecMetaReader.h"
#include "pfecUtility.h"
#include "rapidjson/document.h"

#include <streambuf>
#include <utility>

using namespace rapidjson;

void PFECNativeInterface::MetaReader::read_meta_file(std::string filePath)
{
	clear();

	meta_filepath_ = filePath;
	meta_file_dir_ = get_directory(filePath);

    std::ifstream metaFile(filePath);
    if (!metaFile.is_open())
        throw _MSG_EXCEPT_C(ERR_COULDNT_OPEN_STREAM);
	// I do not expect these to be long, so just read it all into a string
	std::string metaFileString((std::istreambuf_iterator<char>(metaFile)), std::istreambuf_iterator<char>());
	metaFile.close();

	// Parse the file string into a JSON structure
	Document analysisJSON;  // Default template parameter uses UTF8 and MemoryPoolAllocator.
	if (analysisJSON.Parse(metaFileString.c_str()).HasParseError())
		throw _MSG_EXCEPT_C(ERR_META_INVALID);

	// Parse the JSON structure into the class

	if (!analysisJSON.HasMember("analysis")) throw _MSG_EXCEPT_C(ERR_META_INVALID);
	const Value& analysisItem = analysisJSON["analysis"];

	// Read name
	analysis_name_ = analysisItem["name"].GetString();

	// Read the time steps
	// Expected time step format:
	// {
	//		"id" : 0,
	//		"time" : 0.0,
	//		"mesh" : "Test.pfecm",
	//		"orientation" : "Test.pfeco",
	//		"data" : {
	//			file: "Test.pfecd",
	//			"arrays" : [
	//				{"type" : "Nodal", "name" : "Displacement"}
	//			]
	//		},
	//		"groups" : [ "Test.pfecg" ]
	// }
	const Value& timeStepsItem = analysisItem["timesteps"];
	auto nameDataTypeMap = get_name_data_type_map();
	for (auto& stepItem : timeStepsItem.GetArray()) {
		TimeStepMeta stepMeta;

		if (stepItem.HasMember("id"))
			stepMeta.time_step_id_ = stepItem["id"].GetInt();
		else throw _MSG_EXCEPT_C(ERR_META_INVALID_TIMESTEP);

		if (stepItem.HasMember("time"))
			stepMeta.time_value_ = stepItem["time"].GetDouble();
		else throw _MSG_EXCEPT_C(ERR_META_INVALID_TIMESTEP);

		if (stepItem.HasMember("mesh"))
			stepMeta.mesh_file_path_ = stepItem["mesh"].GetString();
		else if (time_step_metas_.size() > 0 && time_step_metas_.back().mesh_file_path_ != "")
			stepMeta.mesh_file_path_ = time_step_metas_.back().mesh_file_path_;
		else throw _MSG_EXCEPT_C(ERR_META_INVALID_TIMESTEP);

		if (stepItem.HasMember("orientation"))
			stepMeta.orient_file_path_ = stepItem["orientation"].GetString();
		else if (time_step_metas_.size() > 0 && time_step_metas_.back().orient_file_path_ != "")
			stepMeta.orient_file_path_ = time_step_metas_.back().orient_file_path_;
		else stepMeta.orient_file_path_ = "";

		if (stepItem.HasMember("data")) {
			const Value& dataItem = stepItem["data"];
			if (dataItem.HasMember("file"))
				stepMeta.data_file_path_ = dataItem["file"].GetString();
			else throw _MSG_EXCEPT_C(ERR_META_INVALID_TIMESTEP);

			if (dataItem.HasMember("arrays")) {
				const Value& arrayItem = dataItem["arrays"];
				for (Value::ConstValueIterator dataArrayIter = arrayItem.Begin(); dataArrayIter != arrayItem.End(); ++dataArrayIter) {
					std::string arrayType, arrayName;
					if ((*dataArrayIter).HasMember("type"))
						arrayType = (*dataArrayIter)["type"].GetString();
					else throw _MSG_EXCEPT_C(ERR_META_INVALID_TIMESTEP);
					if ((*dataArrayIter).HasMember("name"))
						arrayName = (*dataArrayIter)["name"].GetString();
					else throw _MSG_EXCEPT_C(ERR_META_INVALID_TIMESTEP);
					stepMeta.data_types_.push_back(nameDataTypeMap[arrayType]);
					stepMeta.data_names_.push_back(arrayName);
				}
			}
			else throw _MSG_EXCEPT_C(ERR_META_INVALID_TIMESTEP);

		}

		if (stepItem.HasMember("groups")) {
			const Value& groupsItem = stepItem["groups"];
			for (Value::ConstValueIterator groupIter = groupsItem.Begin(); groupIter != groupsItem.End(); ++groupIter) {
				if (!(*groupIter).IsString()) throw _MSG_EXCEPT_C(ERR_META_INVALID_TIMESTEP);
				stepMeta.group_file_paths_.push_back((*groupIter).GetString());
			}
		}
		else if (time_step_metas_.size() > 0 && time_step_metas_.back().group_file_paths_.size() > 0)
			stepMeta.group_file_paths_ = time_step_metas_.back().group_file_paths_;

		// Save the meta data and update the id map
		time_step_metas_.push_back(stepMeta);
		time_step_id_map_[stepMeta.time_step_id_] = time_step_metas_.size() - 1;
	}
	// End JSON parse
}

std::vector<int> PFECNativeInterface::MetaReader::get_output_time_step_ids()
{
	std::vector<int> stepIDs;
	for (auto& timeStep : time_step_metas_)
		stepIDs.push_back(timeStep.time_step_id_);
	return stepIDs;
}

std::vector<double> PFECNativeInterface::MetaReader::get_output_time_step_vals()
{
	std::vector<double> stepTimes;
	for (auto& timeStep : time_step_metas_)
		stepTimes.push_back(timeStep.time_value_);
	return stepTimes;
}

std::unique_ptr<PFECNativeInterface::MeshReader> PFECNativeInterface::MetaReader::get_mesh_reader(int timeStepID)
{
	if (time_step_id_map_.find(timeStepID) != time_step_id_map_.end() && time_step_metas_[time_step_id_map_[timeStepID]].mesh_file_path_ != "") {
		std::unique_ptr<MeshReader> reader_ptr = std::unique_ptr<MeshReader>(new MeshReader());
		reader_ptr->open(meta_file_dir_ + "/" + time_step_metas_[time_step_id_map_[timeStepID]].mesh_file_path_);
		return std::move(reader_ptr);
	}
	else return nullptr;
}

std::unique_ptr<PFECNativeInterface::DataReader> PFECNativeInterface::MetaReader::get_orient_reader(int timeStepID)
{
	if (time_step_id_map_.find(timeStepID) != time_step_id_map_.end() && time_step_metas_[time_step_id_map_[timeStepID]].orient_file_path_ != "") {
		std::unique_ptr<DataReader> reader_ptr = std::unique_ptr<DataReader>(new DataReader());
		reader_ptr->open(meta_file_dir_ + "/" + time_step_metas_[time_step_id_map_[timeStepID]].orient_file_path_);
		return std::move(reader_ptr);
	}
	else return nullptr;
}

std::unique_ptr<PFECNativeInterface::DataReader> PFECNativeInterface::MetaReader::get_data_reader(int timeStepID)
{
	if (time_step_id_map_.find(timeStepID) != time_step_id_map_.end() && time_step_metas_[time_step_id_map_[timeStepID]].data_file_path_ != "") {
		std::unique_ptr<DataReader> reader_ptr = std::unique_ptr<DataReader>(new DataReader());
		reader_ptr->open(meta_file_dir_ + "/" + time_step_metas_[time_step_id_map_[timeStepID]].data_file_path_);
		return std::move(reader_ptr);
	}
	else return nullptr;
}

std::vector<std::string> PFECNativeInterface::MetaReader::get_group_file_paths(int timeStepID)
{
	if (time_step_id_map_.find(timeStepID) != time_step_id_map_.end())
		return time_step_metas_[time_step_id_map_[timeStepID]].group_file_paths_;
	else
		return std::vector<std::string>();
}

std::set<std::pair<PFECNativeInterface::DataType, std::string>> PFECNativeInterface::MetaReader::get_common_data_outputs() const
{
	std::set<std::pair<PFECNativeInterface::DataType, std::string>> commonDataOutputs;
	for (auto& timeStep : time_step_metas_) {
		for (auto i = 0; i < timeStep.data_types_.size(); i++) {
			const auto dataPair = std::pair<PFECNativeInterface::DataType, std::string>(timeStep.data_types_[i], timeStep.data_names_[i]);
			commonDataOutputs.insert(dataPair);
		}
	}
	return commonDataOutputs;
}

void PFECNativeInterface::MetaReader::clear()
{
	analysis_name_ = "";
	time_step_id_map_.clear();
	time_step_metas_.clear();
	mesh_is_time_dependent_ = false;
	meta_filepath_ = "";
	meta_file_dir_ = "";
}

void PFECNativeInterface::TimeStepMeta::extract_data_array_info(std::string metaFileDir)
{
	if (data_file_path_ != "") {
		DataReader dr;
		dr.open(metaFileDir + "/" + data_file_path_);
		std::vector<PFECNativeInterface::PFEC_INT> dataIDs;
		std::vector<PFECNativeInterface::PFEC_INT> dataNumComp;
		std::vector< std::vector<std::string> > dataCompNames;
		dr.get_data_info(dataIDs, data_types_, dataNumComp, data_names_, dataCompNames);
	}
}
