

namespace
{
    bool validate_filename_extension(const std::string& filename, const std::string& expected_extension)
    {
        const size_t dot = filename.find_last_of(".");
        const std::string extension = filename.substr(dot + 1);
        //const std::string name = filename.substr(0,dot);
        return (dot != std::string::npos && extension.compare(expected_extension) == 0);
    }

    // copied from "map_masked_velocity_to_vtk.cpp", 3/18/2022
    template<typename T>
    void write_scalar_slice_to_vtk(uint64_t nx, uint64_t ny, uint64_t nz, const std::vector<std::array<double, 3>>& points, 
    const std::vector<T>& scalar_data, const std::string& label, const std::string& filename)
    {
        std::string vtk_data_type;
        if (typeid(T) == typeid(bool)) vtk_data_type = "bit";
        else if (typeid(T) == typeid(double)) vtk_data_type = "double";
        else if (typeid(T) == typeid(float)) vtk_data_type = "float";
        else if (typeid(T) == typeid(char)) vtk_data_type = "char";
        else if (typeid(T) == typeid(short)) vtk_data_type = "short";
        else if (typeid(T) == typeid(int)) vtk_data_type = "int";
        else if (typeid(T) == typeid(long)) vtk_data_type = "long";
        else if (typeid(T) == typeid(unsigned char)) vtk_data_type = "unsigned_char";
        else if (typeid(T) == typeid(unsigned short)) vtk_data_type = "unsigned_short";
        else if (typeid(T) == typeid(unsigned int)) vtk_data_type = "unsigned_int";
        else if (typeid(T) == typeid(unsigned long)) vtk_data_type = "unsigned_long";
        else
        {
            std::cerr << "error: unsupported type" << std::endl;	
            exit(1);	
        }

        std::ofstream file(filename, std::ios::binary);
        file << "# vtk DataFile Version 2.3" << std::endl; // prior to this version, numComp was not supported
        file << "DANWUZHERE" << std::endl;
        file << "ASCII" << std::endl;
        file << "DATASET STRUCTURED_GRID" << std::endl;
        file << "DIMENSIONS" << " " << nx << " " << ny << " " << nz << std::endl;
        file << "POINTS" << " " << nx*ny*nz << " " << "double" << std::endl;

        for (auto p : points)
            file << p[0] << " " << p[1] << " " << p[2] << std::endl;

        file << "POINT_DATA" << " " << nx*ny*nz << std::endl;
        file << "SCALARS " << label << " " << vtk_data_type << " " << 1 << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;

        for (auto v : scalar_data) 
            file << v << std::endl;

        file.close();
    }

    void write_to_vtk(const VoxelModel<3>& vmodel, const std::string& filename)
    {
        const uint64_t nx = vmodel.voxel_box.n(0);
        const uint64_t ny = vmodel.voxel_box.n(1);
        const uint64_t nz = vmodel.voxel_box.n(2);

        assert(vmodel.voxel_box.num_cells == nx*ny*nz);

        std::vector<std::array<double, 3>> points;
        for (uint64_t rank = 0; rank < nx*ny*nz; ++rank)
        {
            const Vector<3> cell_center = vmodel.voxel_box.position(rank);
            std::array<double, 3> point({cell_center(0), cell_center(1), cell_center(2)});
            points.push_back(point);
        }

        // apparent glitch in visit prevents proper handling of bool / "bit" data in vtk files
        // here is a hack to sidestep this issue:
        std::vector<int> hack;
        for (auto b : vmodel.voxel_mask)
            hack.push_back(b ? 1 : 0);

        write_scalar_slice_to_vtk(nx,ny,nz, points, hack, "mask", filename);
    }

    VoxelModel<3> generate_phantom(Index<3> n, Vector<3> d)
    {
        const double lx = n(0) * d(0);
        const double ly = n(1) * d(1);
        const double lz = n(2) * d(2);
        const double r = std::min(lx, std::min(ly,lz))/2;
        const double r_sq = r*r;

        VoxelBox<3> vbox(n,d);

        std::vector<bool> mask(vbox.num_cells, false);

        for (uint64_t rank=0; rank < vbox.num_cells; ++rank)
        {
            const Vector<3>& p = vbox.position(rank);
            const Vector<3>& o = vbox.position(0);

            if ((p-o).squaredNorm() < r_sq)
                mask[rank] = true;
        }

        VoxelModel<3> full_phantom(vbox, mask);

        return full_phantom;
    }

    std::vector<bool> uniyon(const std::vector<bool>& mask1, const std::vector<bool>& mask2)
    {
        assert(mask1.size() == mask2.size());
        std::vector<bool> rtn = mask1;
        for (uint64_t i=0; i<mask1.size(); ++i)
            rtn[i] =  rtn[i] || mask2[i];        
        return rtn;
    }

    std::vector<bool> intersection(const std::vector<bool>& mask1, const std::vector<bool>& mask2)
    {
        assert(mask1.size() == mask2.size());
        std::vector<bool> rtn = mask1;
        for (uint64_t i=0; i<mask1.size(); ++i)
            rtn[i] = rtn[i] && mask2[i];        
        return rtn;
    }


}

int main()
{


    Index<3> n ({50,50,50}); // lattice dimensions
    Vector<3> d ({1.0,1.0,1.0}); // voxel dimensions

    VoxelModel<3> floating_phantom = generate_phantom(n,d);

    write_to_vtk(floating_phantom, "../data/floating_phantom.vtk");

    Vector<3> normal(1,1,1);
    Vector<3> point(10,10,10);
    VoxelModel<3> plane_cut = VoxelModel<3>(floating_phantom.voxel_box, floating_phantom.plane_cut_mask(point, normal));

    write_to_vtk(plane_cut, "../data/plane_cut.vtk");


//     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

//     Index<3> radii = 2*Index<3>::Ones();

//     std::vector<bool> mask = floating_phantom.boundary_mask(radii);
  
//     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//     std::cerr << "Time elapsed (seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << std::endl;

//     VoxelModel<3> floating_phantom_boundary(floating_phantom.voxel_box, mask);

//     write_to_vtk(floating_phantom_boundary, "../data/floating_phantom_boundary.vtk");


//     Matrix<3> rotation = Matrix<3>::Identity();
//     // rot(0,0) = 0;
//     // rot(0,1) = 1;
//     // rot(0,2) = 0;
//     // rot(1,0) = 0;
//     // rot(1,1) = 0;
//     // rot(1,2) = 1;
//     // rot(2,0) = 1;
//     // rot(2,1) = 0;
//     // rot(2,2) = 0;
//     Vector<3> translation;
//     translation(0) = 5;
//     translation(1) = 10;
//     translation(2) = -10;

//     const Affinity<3> affinity(rotation, translation);
// //    auto transformed_mask = floating_phantom_boundary.transformed_mask(affinity);
//     VoxelModel<3> fixed_phantom(floating_phantom_boundary.voxel_box, floating_phantom_boundary.transformed_mask(affinity));
//     write_to_vtk(fixed_phantom, "../data/fixed_phantom.vtk");

    // register the two phantoms and hope to match the transformation

    
    return 0;
}