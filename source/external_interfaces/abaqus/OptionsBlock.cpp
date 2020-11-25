#include "mesh/external_interfaces/abaqus/OptionsBlock.h"

#include "utility/ASCIIStreamHelpers.h"

#include <cctype>

namespace mesh
{
OptionBlock::OptionBlock(std::istream& file, std::string line)
{
    bool ends_with_comma = ends_with(line, ",");

    line.erase(std::remove_if(line.begin(), line.end(), std::isspace),
               line.end()); // remove whitespace
    std::replace(line.begin(), line.end(), ',', ' ');
    auto tokens = utility::tokenize(line);

    if (tokens.size() == 0)
    {
        block_type = BlockType::OTHER;
        return;
    }
    if (utility::istr_is_equal(tokens[0], "*NODE"))
    {
        block_type = BlockType::NODE;
    }
    else if (utility::istr_is_equal(tokens[0], "*ELEMENT"))
    {
        block_type = BlockType::ELEMENT;
    }
    else if (utility::istr_is_equal(tokens[0], "*NSET"))
    {
        block_type = BlockType::NODESET;
    }
    else if (utility::istr_is_equal(tokens[0], "*ELSET"))
    {
        block_type = BlockType::ELEMENTSET;
    }
    else if (utility::istr_is_equal(tokens[0], "*ORIENTATION"))
    {
        block_type = BlockType::ORIENTATION;
    }
    else
    {
        block_type = BlockType::OTHER;
    }

    tokens.erase(tokens.begin()); // Remove block type token
    add_keywords(tokens);

    while (ends_with_comma)
    {
        getline(file, line);
        line.erase(std::remove_if(line.begin(), line.end(), std::isspace),
                   line.end()); // remove whitespace
        std::replace(line.begin(), line.end(), ',', ' ');
        auto tokens = utility::tokenize(line);
        add_keywords(tokens);

        bool ends_with_comma = ends_with(line, ",");
    }
}

bool OptionBlock::ends_with(std::string str, std::string ending)
{
    // Compare ending
    if (str.length() >= ending.length())
    {
        return (0 == str.compare(str.length() - ending.length(), ending.length(), ending));
    }
    else
    {
        return false;
    }
}

void OptionBlock::add_keywords(std::vector<std::string>& tokens)
{
    for (auto& token : tokens)
    {
        std::replace(token.begin(), token.end(), '=', ' ');
        auto subtokens = utility::tokenize(token);
        if (subtokens.size() == 2)
        {
            utility::str_to_upper(subtokens[0]);
            keywords.insert({ subtokens[0], subtokens[1] });
        }
        else if (subtokens.size() == 1)
        {
            utility::str_to_upper(subtokens[0]);
            options.insert(subtokens[0]);
        }
    }
}
}