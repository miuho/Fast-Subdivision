/**
 * @file mesh.hpp
 * @brief Mesh class and OBJ loader.
 *
 * @author Eric Butler (edbutler)
 * @author Zeyang Li (zeyangl)
 */

#ifndef _462_SCENE_MESH_HPP_
#define _462_SCENE_MESH_HPP_


#include "math/vector.hpp"
#include <stdlib.h>
#include <vector>
#include <utility>
#include <set>
#include <map>
#include <cassert>

#define BOUNDARY 15462
#define INTERIOR 15662
#define V0_V1 0 
#define V1_V2 1
#define V2_V0 2
#define YES 42
#define MIN(a, b) ((a < b) ? a: b)
#define MAX(a, b) ((a > b) ? a: b)

namespace _462 {

struct MeshVertex
{
    MeshVertex() {
        // initialize the parameters
        type = INTERIOR;
        neighbors_position_sum = Vector3::Zero;
        bn_index = 0;
    };
    
    //Vector4 color;
    Vector3 position;
    Vector3 normal;
    Vector2 tex_coord;
    // index into the vertex list of the neighbor of this vertex
    std::set< int > neighbors;
    // either boundary or interior vertex
    int type;
    // index into the vertex list of the boundary neighbors
    // if the vertex is interior, it should be empty
    // if the vertex is boundary, it should have 2 boundary neighbors
    int boundary_neighbors[2];
    // keep count of the boundary_neighbors stored
    unsigned int bn_index;
    // records the sum of neighbor vertices' position
    Vector3 neighbors_position_sum;
};

struct MeshTriangle
{
    MeshTriangle() {
        //initialize the parameters
        odd_vertices_count = 0;
    };
    
    // index into the vertex list of the 3 vertices
    unsigned int vertices[3];
    // index into the odd vertices list of the 3 new odd vertices
    unsigned int odd_vertices[3];
    // count the number of odd vertices generated
    int odd_vertices_count;
    // assume v0, v1, v2 are counter-clockwise stored in vertices[3]
    // we have the 3 new odd vertices as o0, o1, o2 in odd_vertices[3]
    // we define o0, o1, o2 that represent these specific vertices
    //      v0
    //     / \
    //   o0   o2
    //   /     \
    // v1---o1--v2  
    // V0_V1 = 0, V1_V2 = 1, V2_V0 = 2
    // thus, 3 new odd vertices are also counter-clockwise stored
};

struct MeshEdge
{
    // index into the vertex list of the 2 end points
    unsigned int end_points[2];
    // number of triangles share this edge
    unsigned int triangle_count;
    // index into the vertex list of the vertices that belong to
    // the triangles that share this edge but are not end points
    // a.k.a. the third point in the triangle beside the end points
    // Since we can assume there is no more than 2 triangles sharing
    // 1 edge, there cannot be more than 2 third points for each edge
    unsigned int third_points[2];
    // index into the triangle list of the triangles that share this edge
    unsigned int triangles[2];
    // records the orientation of the edge in the triangle(s)
    // V0_V1 or V1_V2 or V2_V0
    // which corresponds to the odd_vertices[3] in MeshTriangle
    unsigned int orientations[2];
};

/**
 * A mesh of triangles.
 */
class Mesh
{
public:

    Mesh();
    virtual ~Mesh();

    typedef std::vector< MeshTriangle > MeshTriangleList;
    typedef std::vector< MeshVertex > MeshVertexList;
    typedef std::vector< MeshEdge > MeshEdgeList;

    // The list of all triangles in this model.
    MeshTriangleList triangles;

    // The list of all vertices in this model.
    MeshVertexList vertices;

    // scene loader stores the filename of the mesh here
    std::string filename;

    bool has_tcoords;
    bool has_normals;
    int has_colors;

    // Loads the model into a list of triangles and vertices.
    bool load();

    // Creates opengl data for rendering and computes normals if needed
    bool create_gl_data();

    bool subdivide();

    // Renders the mesh using opengl.
    void render() const;
private:
    typedef std::vector< float > FloatList;
    typedef std::vector< unsigned int > IndexList;

    // the vertex data used for GL rendering
    FloatList vertex_data;
    // the index data used for GL rendering
    IndexList index_data;

    // prevent copy/assignment
    Mesh( const Mesh& );
    Mesh& operator=( const Mesh& );

};


} /* _462 */

#endif /* _462_SCENE_MESH_HPP_ */
