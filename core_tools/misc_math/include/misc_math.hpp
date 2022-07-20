namespace
{
    template<typename T>
    void cartesian_product_recurse(std::vector<std::vector<T>>& accum, std::vector<T>& stack, const std::vector<std::vector<T>>& sequences, uint64_t index)
    {
        std::vector<T> sequence = sequences[index];
        for(uint64_t i : sequence)
        {
            stack.push_back(i);
            if(index ==0)
                accum.push_back(stack);
            else
                cartesian_product_recurse(accum, stack, sequences, index-1);
            stack.pop_back();
        }

    }
}

template<uint64_t dim>
uint64_t MiscellaneousMath::closest_point(const std::vector<Vector<dim>>& points, const Vector<dim>& query_point)
{
    assert(!points.empty());
    uint64_t index_of_closest = 0;
    double min_dist = (query_point-points[0]).squaredNorm();

    for (uint64_t i=0; i<points.size(); ++i)
    {
        const double dist = (query_point-points[i]).squaredNorm();
        if (dist < min_dist)
        {
            min_dist = dist;
            index_of_closest = i;
        }
    }
    return index_of_closest;
}

template<uint64_t dim>
Vector<dim> MiscellaneousMath::projection_to_line(
    const Vector<dim>& point_on_line, 
    const Vector<dim>& direction_of_line, 
    const Vector<dim>& query_point)
{
    const double c = (point_on_line - query_point).dot(direction_of_line) / direction_of_line.squaredNorm();
    return c*direction_of_line + point_on_line;
}

template<uint64_t dim>
std::array<Vector<dim>, 2> MiscellaneousMath::bounding_box(const std::vector<Vector<dim>>& points)
{
    assert(!points.empty());

    // auto mins = std::numeric_limits<double>::max() * Vector<dim>::Ones(dim);
    // auto maxs = -std::numeric_limits<double>::max() * Vector<dim>::Ones(dim);
    // for (auto point : points)
    //     for (uint64_t i=0; i<dim; ++i)
    //     {
    //         if (point(i) < mins(i))
    //             mins(i) = point[i];
    //         if (point(i) > maxs(i))
    //             maxs(i) = point[i];
    //     }

    auto mins = points.front();
    auto maxs = points.front();
    for (auto point : points)
        for (uint64_t i=0; i<dim; ++i)
        {
            if (point(i) < mins(i))
                mins(i) = point[i];
            else if (point(i) > maxs(i))
                maxs(i) = point[i];
        }

    return std::array<Vector<dim>, 2>({mins, maxs});
}

template<typename T>
std::vector<std::vector<T>> MiscellaneousMath::cartesian_product(const std::vector<std::vector<T>>& sequences)
{
    std::vector<std::vector<T>> accum;
    std::vector<T> stack;
    if (sequences.size() > 1)
        cartesian_product_recurse(accum, stack, sequences, sequences.size()-1);
    return accum; 
}

template<uint64_t dim>
std::vector<Index<dim>> MiscellaneousMath::cartesian_product(const Index<dim>& min_corner, const Index<dim>& max_corner)
{
    // TODO: implement this directly (will result in simpler, faster code)

    std::vector<std::vector<uint64_t>> sequences;

    for (uint64_t i=0; i<dim; ++i)
    {
        std::vector<uint64_t> sequence;

        for (uint64_t j=min_corner(i); j<= max_corner(i); ++j)
            sequence.push_back(j);

        sequences.push_back(sequence);
    }

    std::vector<std::vector<uint64_t>> temp = cartesian_product(sequences);

    std::vector<Index<dim>> rtn;    
    for (auto t : temp)
    {
        assert(t.size() == dim);
        Index<dim> index;
        for (uint64_t i=0; i<dim; ++i)
            index(i) = t[i];
        rtn.push_back(index);        
    }
    return rtn;
}


template<uint64_t dim>
Eigen::Matrix<double, dim-1, dim> MiscellaneousMath::projection_matrix(const Vector<dim>& normal)
{
    static_assert(dim > 0, "logic error");
    assert(normal.norm() != 0);

    Eigen::JacobiSVD<Vector<dim>> svd(normal, Eigen::ComputeFullU | Eigen::ComputeFullV);

    auto kernel = svd.matrixU().transpose().bottomRows(dim-1);

    const double tolerance = 1e-8;
    if ((kernel*normal).norm() >= tolerance)
    {
        std::cerr << "error in projection matrix computation";
        assert(false);
    }

    return kernel;
}

template<uint64_t dim>
double MiscellaneousMath::diameter(const std::vector<Vector<dim>>& points)
{
    double rtn = 0; 
    for (uint64_t i=0; i<points.size(); ++i)
        for (uint64_t j=i+1; j< points.size(); ++j)
        {
            double d = (points[i]-points[j]).norm();
            if (d > rtn)
                rtn = d;
        }
    return rtn;
}

template<uint64_t dim>
Vector<dim> MiscellaneousMath::centroid(const std::vector<Vector<dim>>& points)
{
    Vector<dim> sum = Vector<dim>::Zero();
    for(auto x : points)
        sum = sum + x;

    return sum * (1.0 / points.size());
}

template<uint64_t dim>
std::pair<uint64_t, std::array<uint64_t,3>> MiscellaneousMath::rank_and_signature(Matrix<dim> A, double tolerance)
{
    Eigen::FullPivLU<Matrix<dim>> lu_decomp(A);

    uint64_t pos_evalues = 0;
    uint64_t neg_evalues = 0;
    uint64_t zero_evalues = 0;

    auto evalues = A.eigenvalues();        

    for (auto evalue : evalues)
    {
        assert(fabs(evalue.imag()) < tolerance);
        if (evalue.real() > tolerance)
            pos_evalues++;
        else if (evalue.real() < -tolerance)
            neg_evalues++;
        else
            zero_evalues++;
    }
    
    std::array<uint64_t,3> signature({pos_evalues, neg_evalues, zero_evalues});
    uint64_t rank = lu_decomp.rank();



    // std::cerr << "evalues: " << std::endl;
    // for (auto evalue : evalues)
    // {
    //     std::cerr << evalue.real() << " ";
    // }
    // std::cerr << std::endl;


    return make_pair(rank, signature);
}

template <uint64_t dim>
Vector<dim> MiscellaneousMath::gaussian_vector(double mu, double sigma)
{
    assert(sigma >=0);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mu, sigma);

    Vector<dim> v;
    for (uint64_t i=0; i<dim; ++i)
        v[i] = distribution(generator);

    return v;
}


