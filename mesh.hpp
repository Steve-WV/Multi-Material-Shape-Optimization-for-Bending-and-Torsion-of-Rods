#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <vector>
#include <set>

#include "edge.hpp"
#include "utils.hpp"

class Tri;
class Vertex;

// The mesh.
class Mesh  {
	public:
  
  //@{
  /// Typedefs for vertex and triangle containers in the mesh.
  typedef             std::vector<Vertex*> VertexCt;
  typedef             std::vector<Tri*> TriCt;
	
  typedef             TriCt::iterator TriIt;
  typedef             TriCt::const_iterator const_TriIt;
  typedef             VertexCt::iterator   VertexIt;
  typedef             VertexCt::const_iterator   const_VertexIt;

  typedef             std::set<Edge, Edge::edge_comp>  EdgeSet;
  typedef             std::set<Vertex*>  VertexSet;
  	
  /// Constructor initializes empty mesh.
  inline Mesh( void ) : m_vc(), m_tc() {}
  /// Destructor cleans up vertex and triangle containers.
  inline ~Mesh( void );

  void Load(std::string);
  	
  //@{
  /// Iterators for the vertex and triangle containers.
  inline TriIt          TriBegin() { return m_tc.begin(); }
  inline TriIt          TriEnd() { return m_tc.end(); }
  inline const_TriIt    TriBegin() const { return m_tc.begin(); }
  inline const_TriIt    TriEnd() const { return m_tc.end(); }
  inline int            TriSize() const { return static_cast<int>(m_tc.size()); }
  inline VertexIt       VertexBegin() { return m_vc.begin(); }
  inline VertexIt       VertexEnd() { return m_vc.end(); }
  inline VertexIt       BdryVertexBegin() { return m_bv.begin(); }
  inline VertexIt       BdryVertexEnd() { return m_bv.end(); }
  inline const_VertexIt VertexBegin() const { return m_vc.begin(); }
  inline const_VertexIt VertexEnd() const { return m_vc.end(); }
  inline const_VertexIt BdryVertexBegin() const { return m_bv.begin(); }
  inline const_VertexIt BdryVertexEnd() const { return m_bv.end(); }
  
  inline int              VertexSize() const { return static_cast<int>(m_vc.size()); }
  //@}

  //@{
  /// Add vertices and triangles to the respective containers.
  inline Vertex* addVertex( Vertex* v ) { m_vc.push_back( v ); return v; }
  inline Tri* addTri( Tri* t ) { m_tc.push_back( t ); return t; }
  //@}

  //@{
  /// Direct accessors to the vertex and triangle containers.
  inline Vertex* &V( int i ) { return m_vc[i]; }
  inline Vertex* V( int i ) const { return m_vc[i]; }
  inline Tri* &T( int i ) { return m_tc[i]; }
  inline Tri* T( int i ) const { return m_tc[i]; }
  //@}
  
  // Calls a function to set up topological data & also mark boundaries & set valences
  inline void Init() {
      // connect up triangles
      BuildTriangleTopology();
      // intialize structures needed by vertex iterators
      //BuildVertexTopology();
  }

 
private:
  /// All vertices.
  VertexCt    m_vc;
  VertexCt    m_bv;
  /// All triangles.
  TriCt   m_tc;
	
  void AddAndOrMatch( EdgeSet&, Vertex*, Vertex*, Tri* );
  void SetBoundary( EdgeSet& );
  //void BuildRing( Tet* t, Vertex *v );
  //void BuildVertexTopology();
  void BuildTriangleTopology();
  
};


#endif /* __MESH_HPP__ */
