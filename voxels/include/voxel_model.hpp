
template<uint64_t dim>
VoxelModel<dim>::VoxelModel(const VoxelBox<dim>& voxel_box_, const std::vector<bool>& voxel_mask_) : voxel_box(voxel_box_), voxel_mask(voxel_mask_)
{
    assert(voxel_mask.size() == voxel_box.num_cells);
}

template<uint64_t dim>
void VoxelModel<dim>::display_mask_stats(const std::vector<bool>& mask)
{
    std::cerr << "mask.size(): " << mask.size() << std::endl;
    uint64_t count=0;
    for (auto i : mask)
        if (i)
            count++;
    std::cerr << "count: " << count << std::endl;
}

template<uint64_t dim>
std::vector<bool> VoxelModel<dim>::transformed_mask(const Affinity<dim> affinity) const
{
    std::vector<bool> transformed(voxel_box.num_cells, false);
    for (uint64_t rank = 0; rank < voxel_box.num_cells; ++rank)
        if (voxel_mask[rank])
        {
            Vector<3> position = voxel_box.position(rank);
            Vector<3> transformed_position = affinity(position);

            // assert(transformed_position == affinity(position));

            Index<3> index_of_transformed_position = voxel_box.index(transformed_position);

            // assert(voxel_box.to_rank(index_of_transformed_position) == rank);

            if (!voxel_box.out_of_range(index_of_transformed_position))
                transformed[voxel_box.to_rank(index_of_transformed_position)] = true;
        }
    return transformed;
}

template<uint64_t dim>
std::vector<bool> VoxelModel<dim>::boundary_mask(const Index<dim>& radii) const
{
    std::vector<bool> on_boundary(voxel_box.num_cells, false);
    for (uint64_t rank = 0; rank < voxel_box.num_cells; ++rank)
    {
        auto neighbours = voxel_box.neighbour_ranks(rank, radii);
        bool near_interior = false;
        bool near_exterior = false;
        for (auto neighbour : neighbours)
        {
            if (voxel_mask[neighbour])
                near_interior = true;
            else
                near_exterior = true;
            
            if (near_exterior && near_interior)
            {
                on_boundary[rank] = true;
                break;
            }
        }
    }
    return on_boundary;
}

template<uint64_t dim>
std::vector<bool> VoxelModel<dim>::plane_cut_mask(const Vector<dim>& point, const Vector<dim>& normal)
{
    std::vector<bool> rtn = voxel_mask;
    for (uint64_t rank = 0; rank<voxel_box.num_cells; ++rank)
    {
        Vector<dim> position = voxel_box.position(rank);
        if ((position-point).dot(normal) > 0)
            rtn[rank] = false;
    }
    return rtn;
}

template<uint64_t dim>
void VoxelModel<dim>::write_to_vtk(const std::string& filename) const
{
    const uint64_t nx = voxel_box.n(0);
    const uint64_t ny = voxel_box.n(1);
    const uint64_t nz = voxel_box.n(2);

    assert(voxel_box.num_cells == nx*ny*nz);

    std::vector<std::array<double, 3>> points;
    for (uint64_t rank = 0; rank < nx*ny*nz; ++rank)
    {
        const Vector<3> cell_center = voxel_box.position(rank);
        std::array<double, 3> point({cell_center(0), cell_center(1), cell_center(2)});
        points.push_back(point);
    }

    // apparent glitch in visit prevents proper handling of bool / "bit" data in vtk files
    // here is a hack to sidestep this issue:
    std::vector<int> hack;
    for (auto b : voxel_mask)
        hack.push_back(b ? 1 : 0);

    Miscellaneous::write_scalar_slice_to_vtk(nx,ny,nz, points, hack, "mask", filename);
}

