#pragma once

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace mesh
{
/// @brief A struct for parsing type, keywords, and options of a block in an Abaqus file
struct OptionBlock
{
    enum class BlockType
    {
        NODE,
        ELEMENT,
        NODESET,
        ELEMENTSET,
        ORIENTATION,
        STOP,
        OTHER
    };

    /// @param file - the istream that the .inp file is being read from.  This parameter is only
    /// used if the keywards are long enough in the .inp file to be on multiple lines.  Otherwise,
    /// the stream is not modified.
    ///
    /// @param line - the current line from the .inp file marking the beginning of a block, which
    /// should start with a '*'
    OptionBlock(std::istream& file, std::string line);

    BlockType block_type;
    /// Keywords are stored in key-value pairs.  Keys are forced to be uppercase in the map (even if
    /// the .inp file has them as mixed or lower case), but values perserve case.
    std::map<std::string, std::string> keywords;
    std::set<std::string> options;

    bool has_option(const std::string& option) const;

  private:
    bool ends_with(std::string str, std::string ending);
    void add_keywords(std::vector<std::string>& tokens);
};
}