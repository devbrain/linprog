#include <sstream>
#include <iostream>
#include <fstream>
#include "problem_description.hpp"

problem_description::problem_description(const char* path)
{
    std::ifstream ifs(path);
    if (!ifs.good())
    {
        throw std::runtime_error("file not found");
    }
    ifs >> m >> n;
    basic.resize(m);
    for (unsigned int i = 0; i < m; i++)
    {
        ifs >> basic[i];
    }
    non_basic.resize(n);
    for (unsigned int i = 0; i < n; i++)
    {
        ifs >> non_basic[i];
    }
    b.resize(m);

    for (unsigned int i = 0; i < m; i++)
    {
        ifs >> b[i];
    }
    a.resize(m*n);
    for (std::size_t i = 0; i < a.size(); i++)
    {
        ifs >> a[i];
    }
    ifs >> z0;
    c.resize(n);
    for (unsigned int i = 0; i < n; i++)
    {
        ifs >> c[i];
    }
}