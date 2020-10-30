#pragma once

#include <iostream>
#include <string>
#include <vector>

namespace PFECNativeInterface {

    template<typename T>
    bool in_vec(const std::vector<T> &vec, const T &val)
    {
		for (auto i : vec)
			if (i == val) return true;
        return false;
    };

    std::string get_directory(const std::string & filePath);
	/// Replaces getline but handles /r, /n and /r/n correctly across platforms
	std::istream& safe_getline(std::istream& is, std::string& t);
};
