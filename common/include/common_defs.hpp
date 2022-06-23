#pragma once

template<uint64_t dim>
Vector<dim> index_to_vector(const Index<dim>& index)
{
    Vector<dim> rtn;
    for (uint64_t i=0; i<dim; ++i)
        rtn(i) = static_cast<double>(index(i));

    return rtn;
}

template<uint64_t dim>
Index<dim> vector_to_index(const Vector<dim>& vec)
{
    Index<dim> rtn;
    for (uint64_t i=0; i<dim; ++i)
        rtn(i) = std::floor(vec(i));

    return rtn;
}
