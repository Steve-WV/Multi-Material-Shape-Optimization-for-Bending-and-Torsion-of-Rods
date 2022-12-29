#ifndef __FACE_HPP__
#define __FACE_HPP__

#include <functional>
#include "utils.hpp"

class Tri;
class Vertex;

/// Faces in the triangulation.
class Edge{
public:
	/// Construction from two vertices, and maybe an adjacent triangle.
	inline Edge( Vertex*, Vertex*, Tri* = NULL );
	/// Copy constructor. Would prefer to do without..
	inline Edge( const Edge& );

	/// Topological information
	inline bool         Boundary() const;

	// Information about the vertices
	inline Vertex*      OtherV( const Vertex* ) const;
	inline bool         Contains( const Vertex* ) const;
	inline bool         IsEdge( const Vertex*, const Vertex* ) const;

	// Accessors
	Vertex*&            a() { return m_vi[0]; }
	Vertex*             a() const { return m_vi[0]; }
	Vertex*&            b() { return m_vi[1]; }
	Vertex*             b() const { return m_vi[1]; }
	Tri*                t() const { return m_t; }
	Tri*&               t() { return m_t; }
	
	// Comparison function for use in Face sets, lexicographical order
	struct edge_comp : public std::binary_function<Edge, Edge, bool> {
		edge_comp() {}
		bool operator()( const Edge f1, const Edge f2 ) const {
			const Vertex *min1 = std::min( f1.a(), f1.b() );
			const Vertex *max1 = std::max( f1.a(), f1.b() );
			const Vertex *min2 = std::min( f2.a(), f2.b() );
			const Vertex *max2 = std::max( f2.a(), f2.b() );
			
			assert( min1 != max1 );
            assert( min2 != max2 );
            
            if ( min1 < min2 ) return true;
			if ( min1 == min2 && max1 < max2 ) return true;
			return false;
		}
	};
	
private:
	//@{
	/// Vertices defining end point	
	Vertex*       m_vi[2];
	/// one of the faces incident with this Face	
	Tri*          m_t;
};

#include "tri.hpp"

inline
	Edge::Edge( Vertex* a, Vertex* b, Tri* t )
		: m_t( t )
{
	assert( a && b ); assert( a != b );
	assert( !t || ( t->Contains(a) && t->Contains(b) ) );
  
	m_vi[0] = a;
	m_vi[1] = b;
}

inline
	Edge::Edge( const Edge& e )
	 	: m_t( e.m_t )
{
	assert( e.m_vi[0] && e.m_vi[1] ); assert( e.m_vi[0] != e.m_vi[1] );
	assert( !e.m_t || ( e.m_t->Contains(e.m_vi[0]) && e.m_t->Contains(e.m_vi[1]) ) );
  
	m_vi[0] = e.m_vi[0];
	m_vi[1] = e.m_vi[1];
}


/// Is the Face on a boundary?
inline bool
Edge::Boundary() const{
	// if the adjacent face does not have a neighbor pointer than this is a 
	// boundary.
	return !m_t->Across( m_vi[0], m_vi[1]);
}

/// Given one vertex, return the other one.
inline Vertex*
Edge::OtherV( const Vertex* v1 ) const{
	assert( v1 ); assert( Contains( v1 ));
    if ( m_vi[0] == v1 ) return m_vi[1];
	return m_vi[0]; // else
}

/// Is the given vertex one of this faces's points.
inline bool
Edge::Contains( const Vertex* v ) const{
	assert( v );
	return v == m_vi[0] || v == m_vi[1];
}

/// Does this Face contain the three vertices
inline bool
	Edge::IsEdge( const Vertex* v1, const Vertex* v2) const
{
	assert( v1 && v2 );
	if ( m_vi[0] == v1 && m_vi[1] == v2 ) return true;
	if ( m_vi[1] == v1 && m_vi[0] == v2 ) return true;
	return false;
}

#endif	/* __Face_HPP__ */
