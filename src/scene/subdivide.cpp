/**
 * @file subdivide.cpp
 *
 * @author HingOn Miu (hmiu)
 * 
 */

#include "scene/mesh.hpp"


namespace _462 {

// a vector that holds all the edges of the triangles.
typedef std::vector< MeshEdge > MeshEdgeList;
// a hashtable that maps a unique key to the index of the edges array.
// the key is generated from the two index of the vertices array, which
// are the endpoints of the edge. in order words, two endpoints map to
// an edge.
typedef std::map< unsigned long long, int > EdgeHashTable;


/**
 * Zudzik's pairing function
 * http://szudzik.com/ElegantPairing.pdf
 * Provided two 32 bit non-negative integer, the function returns a
 * 64 bit non-negative integer, which is unique to the inputs.
 */
unsigned long long two_to_one_mapping(unsigned long long a, 
                                      unsigned long long b) {
    return a >= b ? a * a + a + b : a + b * b;
}

/**
 * Returns 0 if the creation of new edge is successful, and -1 otherwise.
 */
int create_new_edge(EdgeHashTable& edge_ht,
                    MeshEdgeList& edges, 
                    unsigned int& edges_index, unsigned int a,
                    unsigned int b, int orientation,
                    unsigned int third_point,
                    int triangle_index) {
    std::pair< std::map< unsigned long long, int>::iterator, bool> ret;
    ret = edge_ht.insert(std::make_pair<unsigned long long, int>
                        (two_to_one_mapping((unsigned long long)a, 
                        (unsigned long long)b), edges_index));

    // this edge is a new edge
    if (ret.second == true) {
        MeshEdge edge;
        edge.end_points[0] = a;
        edge.end_points[1] = b;
        edge.triangle_count = 1;
        edge.third_points[0] = third_point;
        edge.triangles[0] = triangle_index;
        edge.orientations[0] = orientation;
        // insert the edge to edges array
        edges[edges_index] = edge;
        edges_index++;
        return 0;
    }
    // this edge was already shared by another triangle
    else {
        // get the inserted edge's index of the edges array
        int edge_index = ret.first->second;
        
        // each edge cannot be shared with more than 2 triangles
        if ((edges[edge_index]).triangle_count != 1) {
            return -1;
        } 
        (edges[edge_index]).triangle_count++;
        (edges[edge_index]).third_points[1] = third_point;
        (edges[edge_index]).triangles[1] = triangle_index;
        (edges[edge_index]).orientations[1] = orientation;
        return 0;
    }
}

/**
 * Returns 0 if add the new edges for a triangles is successful,
 * and -1 otherwise.
 */
int add_new_edges(EdgeHashTable& edge_ht, 
                  MeshEdgeList& edges, 
                  unsigned int& edges_index, unsigned int i0,
                  unsigned int i1, unsigned int i2, 
                  int triangle_index) {
    unsigned long long min, max;
    
    // add the edge v0 to v1
    min = (unsigned long long)MIN(i0, i1);
    max = (unsigned long long)MAX(i0, i1);
    if (create_new_edge(edge_ht, edges, edges_index, 
                        min, max, V0_V1, i2, 
                        triangle_index) == -1) {
        return -1;
    }
    
    // add the edge v1 to t2
    min = (unsigned long long)MIN(i1, i2);
    max = (unsigned long long)MAX(i1, i2);
    if (create_new_edge(edge_ht, edges, edges_index,
                        min, max, V1_V2, i0, 
                        triangle_index) == -1) {
        return -1;
    }

    // add the edge v2 to v0
    min = (unsigned long long)MIN(i2, i0);
    max = (unsigned long long)MAX(i2, i0);
    if (create_new_edge(edge_ht, edges, edges_index, 
                        min, max, V2_V0, i1, 
                        triangle_index) == -1) {
        return -1;
    }

    return 0;
}

/**
 * Calculate beta from the number of neighbors
 */
real_t calculate_beta(real_t N) {
    real_t tmp = (3.0/8.0) + (cos((2.0 * PI)/N)/4.0);
    return ((5.0/8.0) - (tmp * tmp))/N;
}

/**
 * Accumulate the sum of neighbors' coordinates
 */
void sum_neighbor_positions(Mesh::MeshVertexList& vertices, unsigned int i0,
                       unsigned int i1, unsigned int i2) {
    std::pair< std::set< int >::iterator, bool> ret;

    // check whether the neighbors were added for each vertices
    ret = (vertices[i0]).neighbors.insert(i1);
    if (ret.second == true) {
        (vertices[i0]).neighbors_position_sum += (vertices[i1]).position;
    }
    ret = (vertices[i0]).neighbors.insert(i2);
    if (ret.second == true) {
        (vertices[i0]).neighbors_position_sum += (vertices[i2]).position;
    }
    ret = (vertices[i1]).neighbors.insert(i0);
    if (ret.second == true) {
        (vertices[i1]).neighbors_position_sum += (vertices[i0]).position;
    }
    ret = (vertices[i1]).neighbors.insert(i2);
    if (ret.second == true) {
        (vertices[i1]).neighbors_position_sum += (vertices[i2]).position;
    }
    ret = (vertices[i2]).neighbors.insert(i1);
    if (ret.second == true) {
        (vertices[i2]).neighbors_position_sum += (vertices[i1]).position;
    }
    ret = (vertices[i2]).neighbors.insert(i0);
    if (ret.second == true) {
        (vertices[i2]).neighbors_position_sum += (vertices[i0]).position;
    }
}

/**
 * Split each triangles into 4 smaller triangles and place them into
 * the new triangles vector.
 */
void split_triangle(Mesh::MeshTriangleList& triangles, unsigned int i,
                    Mesh::MeshTriangleList& new_triangles,
                    unsigned nvertices) {
    unsigned int o0, o1, o2, i0, i1, i2;
    MeshTriangle t1, t2, t3, t4;
    
    // even vertices index
    i0 = (triangles[i]).vertices[0];
    i1 = (triangles[i]).vertices[1];
    i2 = (triangles[i]).vertices[2];
    //odd vertices index
    o0 = (triangles[i]).odd_vertices[0] + nvertices;
    o1 = (triangles[i]).odd_vertices[1] + nvertices;
    o2 = (triangles[i]).odd_vertices[2] + nvertices;
    
    // construct triangles 
    t1.vertices[0] = i0;
    t1.vertices[1] = o0;
    t1.vertices[2] = o2;
    t2.vertices[0] = o0;
    t2.vertices[1] = i1;
    t2.vertices[2] = o1;
    t3.vertices[0] = o0;
    t3.vertices[1] = o1;
    t3.vertices[2] = o2;
    t4.vertices[0] = o2;
    t4.vertices[1] = o1;
    t4.vertices[2] = i2;
    
    new_triangles.push_back(t1);
    new_triangles.push_back(t2);
    new_triangles.push_back(t3);
    new_triangles.push_back(t4);
}

bool Mesh::subdivide()
{
    unsigned nvertices = vertices.size();
    unsigned ntriangles = triangles.size();
    unsigned nedges, new_ntriangles, new_nvertices;
    MeshVertexList odd_vertices;
    MeshVertexList even_vertices;
    MeshEdgeList edges;
    unsigned int edges_index = 0; 
    MeshTriangleList new_triangles;
    EdgeHashTable edge_ht;
    int i, j;
    real_t N, beta;
    unsigned int a, b, c, d, i0, i1, i2;

    // pre-allocate space to avoid slow resize vector operation
    even_vertices.reserve(nvertices + 3 * ntriangles);
    odd_vertices.reserve(3 * ntriangles);
    edges.reserve(3 * ntriangles);
    new_triangles.reserve(4 * ntriangles);
     
    // create the edges list
    for (i = 0; i < ntriangles; i++) {
        i0 = (triangles[i]).vertices[0];
        i1 = (triangles[i]).vertices[1];
        i2 = (triangles[i]).vertices[2];
        
        if (add_new_edges(edge_ht, edges, edges_index, i0, i1, i2, i) == -1) {
            std::cout << "Failed to create edges list" << std::endl;
            return false;
        }
        
        // sum neighbors' coordinates for later even vertices calculation 
        sum_neighbor_positions(vertices, i0, i1, i2);
    }
    
    // get the number of edges generated
    nedges = edges_index;
    // generate the odd vertices
    for (i = 0; i < nedges; i ++) {
        MeshVertex v;
        v.normal = Vector3::Zero;
        a = (edges[i]).end_points[0];
        b = (edges[i]).end_points[1];
        // this edge is shared by 2 triangles
        if ((edges[i]).triangle_count == 2) {
            c = (edges[i]).third_points[0];
            d = (edges[i]).third_points[1];
            
            v.position = ((vertices[a].position * 3.0)/ 8.0) + 
                         ((vertices[b].position * 3.0)/ 8.0) + 
                         (vertices[c].position / 8.0) + 
                         (vertices[d].position / 8.0);           
            
            // add the odd vertex's index to the corresponding triangle
            (triangles[(edges[i]).triangles[0]]).odd_vertices[
                (edges[i]).orientations[0]] = i;
            (triangles[(edges[i]).triangles[0]]).odd_vertices_count++;
            (triangles[(edges[i]).triangles[1]]).odd_vertices[
                (edges[i]).orientations[1]] = i;
            (triangles[(edges[i]).triangles[1]]).odd_vertices_count++;
            
            // check if 3 odd vertices are generated, ready for triangle split
            if ((triangles[(edges[i]).triangles[0]]).odd_vertices_count == 3) {
                split_triangle(triangles, (edges[i]).triangles[0], 
                               new_triangles, nvertices);
            }
            if ((triangles[(edges[i]).triangles[1]]).odd_vertices_count == 3) {
                split_triangle(triangles, (edges[i]).triangles[1], 
                               new_triangles, nvertices);
            }
        }
        // this edge is a boundary edge
        else {
            v.position = (vertices[a].position / 2.0) + 
                         (vertices[b].position / 2.0); 
            
            // identify the boundary vertices
            (vertices[a]).type = BOUNDARY;
            (vertices[a]).boundary_neighbors[(vertices[a]).bn_index] = b;
            (vertices[a]).bn_index++;
            (vertices[b]).type = BOUNDARY;
            (vertices[b]).boundary_neighbors[(vertices[b]).bn_index] = a;
            (vertices[b]).bn_index++;       
            
            // add the odd vertex's index to the corresponding triangle
            (triangles[(edges[i]).triangles[0]]).odd_vertices[
                (edges[i]).orientations[0]] = i;
            (triangles[(edges[i]).triangles[0]]).odd_vertices_count++;
            
            // check if 3 odd vertices are generated, ready for triangle split
            if ((triangles[(edges[i]).triangles[0]]).odd_vertices_count == 3) {
                split_triangle(triangles, (edges[i]).triangles[0], 
                               new_triangles, nvertices);
            }
        }
        odd_vertices.push_back(v);
    }
    
    // generate the even vertices
    for(i = 0; i < nvertices; i++) {
        MeshVertex v;
        v.normal = Vector3::Zero;
        // case for interior vertices
        if ((vertices[i]).type == INTERIOR) {
            N = (real_t)((vertices[i]).neighbors.size());
            beta = calculate_beta(N);
            
            v.position = (1.0 - beta * N) * vertices[i].position + 
                         beta * vertices[i].neighbors_position_sum;
        }
        // case for boundary vertices
        else {
            //assert((vertices[i]).boundary_neighbors.size() == 2);
            a = (vertices[i]).boundary_neighbors[0];
            b = (vertices[i]).boundary_neighbors[1];
            
            v.position = ((vertices[i].position * 3.0)/4.0) +
                         (vertices[a].position / 8.0) +
                         (vertices[b].position /8.0);           
        }
        even_vertices.push_back(v);
    }

    // merge the even vertices list and odd vertices list
    even_vertices.insert(even_vertices.end(), 
                         odd_vertices.begin(), odd_vertices.end());
    
    new_ntriangles = new_triangles.size();
    // calculate the new normals for each vertices
    for (i = 0; i < new_ntriangles; i++) {
        i0 = (new_triangles[i]).vertices[0];
        i1 = (new_triangles[i]).vertices[1];
        i2 = (new_triangles[i]).vertices[2];
        
        Vector3 a = (even_vertices[i1]).position -
                    (even_vertices[i0]).position;
        Vector3 b = (even_vertices[i2]).position -
                    (even_vertices[i0]).position;

        // accumulate the triangle normals for each vertices
        (even_vertices[i0]).normal += cross(a, b);
        (even_vertices[i1]).normal += cross(a, b);
        (even_vertices[i2]).normal += cross(a, b);
    }
    
    new_nvertices = even_vertices.size();
    // normalize each vertices' normal sum
    for (i = 0; i < new_nvertices; i++) {
        (even_vertices[i]).normal = normalize((even_vertices[i]).normal);
    }
     
    // rebuild vertex_data and index_data for GL rendering
    vertices.clear();
    triangles.clear();
    vertices = even_vertices;
    triangles = new_triangles;
    create_gl_data();
    return true;
}

} /* _462 */
