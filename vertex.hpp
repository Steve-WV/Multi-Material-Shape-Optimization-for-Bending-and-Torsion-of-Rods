#ifndef	__VERTEX_HPP__
#define	__VERTEX_HPP__

// include files
#include "names.h"
#include "utils.hpp"
//#include "tet.hpp"
#include "tmv.hpp"

class Mesh;
/// The vertices. They are pretty dumb.
class Vertex {
public:
  
  inline Vertex() : m_boundary(false), m_dirbd(false)
    /*: m_ti( NULL )*/ { }
  
  inline Vertex( Vec2 coord ) : m_coord( coord ), m_boundary(false), m_dirbd(false)
	 /*, m_ti( NULL )*/ { }
  
  inline Vertex( double x, double y) : m_boundary(false), m_dirbd(false)
  	{ m_coord.x() = x; m_coord.y() = y; }

  inline ~Vertex() {}

  // Accessors
  inline Vec2& Coord() { return m_coord; }
  inline Vec2 Coord() const { return m_coord; }
  inline double& x() { return m_coord.x(); }
  inline double x() const { return m_coord.x(); }
  inline double& y() { return m_coord.y(); }
  inline double y() const { return m_coord.y(); }
  
  inline int Index() const { return m_index; }
  inline void SetIndex(int i) { m_index = i; }
  
  // Is this a boundary vertex?
  inline bool Boundary() const { return m_boundary; }
  inline bool DirBdry() const { return m_dirbd; }
  inline bool NeuBdry() const { return m_boundary && !m_dirbd; }
  
  inline void MarkDirBdry() { m_dirbd = true; }
  
  inline double& u() {return *m_u;}
  inline double u() const {return *m_u;}
  inline double& Phi() { return *m_up;}
  inline double Phi() const {return *m_up;}
  
  inline void Attach( double* u_addr ) { m_u = u_addr; }
  inline void Attach1( double* u_addr ) { m_up = u_addr; }   
   
  friend class Mesh;
  
private:
	
  int m_index;
	
  Vec2 m_coord;	
  // One of the incident tets.
  // Clockwise most incident triangle if the one ring is open. Any face if
  // it's closed.
  // not sure yet
  //Tet* m_ti;       

  /// valence of this vertex;
  //int m_valence;

  bool m_boundary;
  bool m_dirbd;
  inline void MarkBoundary() { m_boundary = true; }
  /// is this vertex on the boundary?
  
  double* m_u;
 double* m_up;
  
};

#endif	/* __VERTEX_HPP__ */
