#pragma once

#include "NaturalSort.h"
#include "utility/Filesystem.h"
#include "utility/StringUtils.h"

#include <string>

namespace mesh
{
/// Given the file name for a test data file, this returns the full filepath
static const std::string get_test_data_filepath(const std::string& data_filename)
{
    // Assumes test data is one directory up in another directory called "data"
    return fs::path(__FILE__)
      .parent_path()
      .parent_path()
      .append("data")
      .append(data_filename)
      .string();
}

static std::vector<std::string> get_files_with_extension(const std::string& path,
                                                         const std::string& ext,
                                                         const bool absolute_path)
{
    std::vector<std::string> files;

    for (auto& p : fs::directory_iterator(path))
    {
        if (utility::istr_is_equal(p.path().extension().string(), ext))
        {
            if (absolute_path)
            {
                files.push_back(p.path().string());
            }
            else
            {
                files.push_back(p.path().relative_path().string());
            }
        }
    }
    std::sort(files.begin(), files.end(), doj::alphanum_less<std::string>());

    return files;
}
}