#include "misc.h"

using json = nlohmann::json;
json Miscellaneous::read_metadata(const std::string& filename)
{
    std::stringstream ss;
    std::ifstream in(filename);
    ss << in.rdbuf();
    in.close();
    return json::parse(ss.str());
}

bool Miscellaneous::validate_filename_extension(const std::string& filename, const std::string& expected_extension)
{
    const size_t dot = filename.find_last_of(".");
    const std::string extension = filename.substr(dot + 1);
    //const std::string name = filename.substr(0,dot);
    return (dot != std::string::npos && extension.compare(expected_extension) == 0);
}

void Miscellaneous::write_to_file(const std::vector<std::vector<double>>& xss, const std::vector<std::string>& header, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "couldn't open file: " << filename << std::endl; 
        assert(false);
    }

    assert(!xss.empty());      
    assert(header.size() == xss.size());

    for (auto xs : xss)
    {
        if (xs.size() != xss[0].size())
        {
            std::cerr << "xs.size(): " << xs.size() << std::endl;
            std::cerr << "xss[0].size(): " << xss[0].size() << std::endl;
//            assert(false);
        }
    }

    for (uint64_t ix=0; ix<xss.size(); ++ix)
    {
        if (ix < xss.size() - 1)
            file << header[ix] << ",";
        else
            file << header[ix] << std::endl;
    }


    for (uint64_t i=0; i<xss[0].size(); ++i)
    {
        for (uint64_t ix=0; ix<xss.size(); ++ix)
        {
            if (ix < xss.size() - 1)
                file << (xss[ix])[i] << ",";
            else
                file << (xss[ix])[i] << std::endl;
        }
    }

    file.close();
}

bool Miscellaneous::file_exists(const std::string& filename)
{
    const std::filesystem::path p {filename};
    if (!std::filesystem::exists(p))
    {
        std::cerr << " did not find :  " << filename << std::endl;
        return false;
    }
    else if (!std::filesystem::is_regular_file(p))
    {
        std::cerr << " found :  " << filename << " but it's not a regular file" << std::endl;
        return false;
    }

    return true;
}

bool Miscellaneous::directory_exists(const std::string& directory_name)
{
    const std::filesystem::path p {directory_name};
    if (!std::filesystem::exists(p))
    {
        std::cerr << " did not find :  " << directory_name << std::endl;
        return false;
    }
    else if (!std::filesystem::is_directory(p))
    {
        std::cerr << " found :  " << directory_name << " but it's not a directory" << std::endl;
        return false;
    }
    return true;
}

void Miscellaneous::verify_directory(const std::string& directory)
{
    if (!Miscellaneous::directory_exists(directory))
    {
        std::cerr << "error finding directory: " << directory << std::endl;
        exit(1);
    }
}

void Miscellaneous::initialize_directory(const std::string& directory)
{
    std::filesystem::remove_all(directory);
    std::filesystem::create_directory(directory);
    verify_directory(directory);
}

void Miscellaneous::verify_file(const std::string& filename)
{
    if (!Miscellaneous::file_exists(filename))
    {
        std::cerr << "error finding file: " << filename << std::endl;
        exit(1);
    }
}
