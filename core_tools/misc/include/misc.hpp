

template<typename T>
void Miscellaneous::read_file(std::vector<T>& store, const std::string& filename, uint64_t expected_number_of_elements)
{
	std::ifstream in(filename, std::ios::binary);
	if (!in)
	{
		std::cerr << "file error" << std::endl;
		return;
	}

	in.unsetf(std::ios::skipws); // stop eating new lines in binary mode!

	// get size:
	std::streampos file_size;
	in.seekg(0, std::ios::end);
	auto file_size_in_bytes = in.tellg();
	in.seekg(0, std::ios::beg);

	const uint64_t expected_size_in_bytes = expected_number_of_elements * sizeof(T);
	if (static_cast<uint64_t>(file_size_in_bytes) != expected_size_in_bytes)
	{
		std::cerr << "file error: " << std::endl;
		std::cerr << "read: " << file_size_in_bytes << " bytes, but expected to read " << expected_size_in_bytes << " bytes." << std::endl;
		return;
	}

	// reserve capacity 
	store.resize(expected_number_of_elements);

	// read the data
	in.read((char*) &store[0], expected_size_in_bytes);
	if (static_cast<uint64_t>(in.gcount()) != expected_size_in_bytes)
	{
		std::cerr << "file error: " << std::endl;
		std::cerr << "bytes read: " << in.gcount() << std::endl;
		std::cerr << "expected size in bytes: " << expected_size_in_bytes << std::endl;
		return;
	}

	in.close();
}

template<uint64_t dim>
void Miscellaneous::write_points_to_vtk_file(const std::vector<Vector<dim>>& vertices, const std::string& filename_vtk, bool connect_the_dots)
{
    assert(Miscellaneous::validate_filename_extension(filename_vtk, "vtk"));

    std::ofstream file(filename_vtk);
    if (!file.is_open())
    {
        std::cerr << "file error with vtk file: " << filename_vtk << std::endl;
        exit(1);
    }

    file << "# vtk DataFile Version 4.2" << std::endl;
    file << "definition precedes existence" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET POLYDATA" << std::endl;
    file << "POINTS " << vertices.size() << " " << "float" << std::endl;

    switch (dim)
    {
    case 1:
        for (auto vertex : vertices)
            file << vertex(0) << " " <<  0 << " " <<  0 << " " <<  std::endl;
        break;
    
    case 2:
        for (auto vertex : vertices)
            file << vertex(0) << " " <<  vertex(1) << " " <<  0 << " " <<  std::endl;
        break;
    
    case 3:
        for (auto vertex : vertices)
            file << vertex(0) << " " <<  vertex(1) << " " <<  vertex(2) << " " <<  std::endl;
        break;
    
    default:
        std::cerr << "point dimensions other than 1,2,3 are currently not supported" <<std::endl;
        exit(1);
        break;
    }

    if (connect_the_dots)
    {
        file << "POLYGONS 1 " << vertices.size() << std::endl;
        for (uint64_t i=0; i<vertices.size(); ++i)
            file << vertices.size()-i << ' ';

        // "close the loop"
        file << vertices.size() << std::endl;
        file << std::endl;
    }

    file.close();
}

template<uint64_t dim>
std::vector<Vector<dim>> Miscellaneous::read_points_from_vtk_file(const std::string& filename_vtk)
{
    if (dim != 1 && dim != 2 && dim != 3)
    {
        std::cerr << "error: in file: " << filename_vtk << ", dim = " << dim << ", but point dimensions other than 1,2,3 are currently not supported" << std::endl;
        exit(1);
    }

    assert(Miscellaneous::validate_filename_extension(filename_vtk, "vtk"));

    std::ifstream file(filename_vtk);
    assert(file.is_open());

    // skip to the POINTS section
    std::string line;
    std::string first_word;
    std::string second_word;
    std::string third_word;
    do {
        std::getline(file, line);
        std::stringstream ss(line);
        ss >> first_word >> second_word >> third_word;
    } while (first_word.compare(std::string("POINTS")) != 0);
    uint64_t n = std::stoul(second_word);
    assert(third_word.compare("float") == 0);

    std::vector<Vector<dim>> vertices;

    if (dim == 3)
    {
        for (uint64_t i=0; i<n; ++i)
        {
            Vector<dim> v;
            for (uint64_t j=0; j<3; ++j)
            {
                std::string word;
                file >> word;
                float f = std::stof(word);
                v(j) = static_cast<double>(f);
            }
            vertices.push_back(v);
        }
    }
    else
    {
        for (uint64_t i=0; i<n; ++i)
        {
            Vector<dim> v;
            for (uint64_t j=0; j<3; ++j)
            {
                std::string word;
                file >> word;
                float f = std::stof(word);
                if (j<dim)
                {
                    if (f != 0.0)
                    {
                        std::cerr << "error with file: " << filename_vtk << std::endl;
                        std::cerr << "expected a file with data of dimension : " << dim << std::endl;
                        std::cerr << "component " << j << " was expected to be zero but read " << f << std::endl;
                        exit(1);
                    }
                    v(j) = static_cast<double>(f);
                }
            }
            vertices.push_back(v);
        }
    }

    file.close();

    return vertices;
}

