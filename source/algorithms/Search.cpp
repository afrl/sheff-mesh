#include "..\..\include\mesh\algorithms\Search.h"

#include "mesh/algorithms/Search.h"

namespace mesh
{
std::pair<Point, Point> find_extent(const Mesh& mesh)
{
    Point neg_point = { INFINITY, INFINITY, INFINITY };
    Point pos_point = { -INFINITY, -INFINITY, -INFINITY };
    for (auto& node : mesh.nodes())
    {
        auto& p = node.point();
        for (auto i = 0; i < neg_point.size(); i++)
        {
            if (p[i] < neg_point[i])
            {
                neg_point[i] = p[i];
            }
            if (p[i] > pos_point[i])
            {
                pos_point[i] = p[i];
            }
        }
    }
    return { neg_point, pos_point };
}
std::pair<Point, Point> find_extent(const std::initializer_list<std::pair<Point, Point>>& extents)
{
    Point neg_point = { INFINITY, INFINITY, INFINITY };
    Point pos_point = { -INFINITY, -INFINITY, -INFINITY };
    for (auto& extent : extents)
    {
        for (auto& p : { extent.first, extent.second })
        {
            for (auto i = 0; i < neg_point.size(); i++)
            {
                if (p[i] < neg_point[i])
                {
                    neg_point[i] = p[i];
                }
                if (p[i] > pos_point[i])
                {
                    pos_point[i] = p[i];
                }
            }
        }
    }
    return { neg_point, pos_point };
}
}