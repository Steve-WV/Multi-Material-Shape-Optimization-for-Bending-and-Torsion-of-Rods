#ifndef __TET_HPP__
#define __TET_HPP__

#include <vector>

#include "names.h"
#include "utils.hpp"
#include "tmv.hpp"


class Vertex;

/// Triangles.
class Tri{
public:

	/// Default constructor
	inline Tri() { 
		m_ts[A] = m_ts[B]= m_ts[C] = NULL; }
		/// Constructor from three vertices.
		inline Tri( Vertex*, Vertex*, Vertex* );
		/// Copy constructor.
		inline Tri( const Tri& t );
		/// Destructor removes the triangle connections.
		inline ~Tri();

		/// Query routines
		inline bool         Contains( const Vertex* ) const;
		inline bool         IsAdj( const Tri* ) const;
		inline VertexName   VName( const Vertex* ) const;
		inline VertexName   OtherVName( const Vertex*, const Vertex* ) const;
		inline VertexName   VNameAcross( const Tri* ) const;
  
		/// Adjacency information
		inline Tri*&   Across( VertexName );
		inline Tri*    Across( VertexName ) const;
		inline Tri*&   Across( const Vertex*, const Vertex* );
		inline Tri*    Across( const Vertex*, const Vertex* ) const;
		inline Tri*&   Across( const Vertex* );
		inline Tri*    Across( const Vertex* ) const;
  
		inline Vertex*&     OtherV( const Vertex*, const Vertex* );
		inline Vertex*      OtherV( const Vertex*, const Vertex* ) const;
		inline Vertex*&     Across( const Tri* );
		inline Vertex*      Across( const Tri* ) const;
  
		/// Connect Bits related stuff
		inline void         Unlink( const Tri* );
		inline void         SetConnected( const Vertex * );
		inline bool         AllVerticesConnected() const;

		/// Accessors
		inline Vertex*             a( void ) const { return m_vs[A]; }
		inline Vertex*             b( void ) const { return m_vs[B]; }
		inline Vertex*             c( void ) const { return m_vs[C]; }
		inline Vertex*             v( const VertexName vn ) const { assert( vn != VNONE ); return m_vs[vn]; }
     
		inline double Area() const { return m_area; }
		inline double Init_area();
	
		inline void Calc_Coord(Vec2*) const; 
		inline void Calc_s(int, double*) const;
		inline void Calc_u(double*) const;
		inline void Calc_Phi(double*) const; 
		inline void Init_grad_s();
		inline void Calc_grad_s( int, Vec2& ) const;
		inline void Calc_grad_u( Vec2& )const;
		inline void Calc_gradPhi( Vec2& ) const; 
	        

	private:
		/// Vertices of the tet	
		Vertex*		m_vs[3];
		/// Triangle across from vs[i]	
		Tri*			m_ts[3];
    
		// some properties of triangles
		double m_area;
		Mat22 Ftinv;
		Vec2 grad_s[3];
	
	};

	//#include "vertex.hpp"

	inline Tri::Tri( Vertex* a, Vertex* b, Vertex* c )
	{
		assert( a && b && c );
		assert( ( a != b ) && ( b != c ) && ( c != a ) );
		// set the vertices
		m_vs[A] = a;  m_vs[B] = b;  m_vs[C] = c;
		// set all the connect bits to false
		/*m_cb[A] = m_cb[B] = m_cb[C] = m_cb[D] = false;*/
		// set the adjacencies to null for now
		m_ts[A] = m_ts[B] = m_ts[C] = NULL;
	}

	inline Tri::Tri( const Tri& t ) {
		m_vs[A] = t.m_vs[A]; m_vs[B] = t.m_vs[B]; m_vs[C] = t.m_vs[C];
		m_ts[A] = t.m_ts[A]; m_ts[B] = t.m_ts[B]; m_ts[C] = t.m_ts[C];
		/*m_cb[A] = t.m_cb[A]; m_cb[B] = t.m_cb[B]; m_cb[C] = t.m_cb[C]; m_cb[D] = t.m_cb[D];*/
  
	}

	inline Tri::~Tri( void ) {
		// unlink from neighbors
		if( m_ts[A] ){ m_ts[A]->Unlink( this ); m_ts[A] = NULL; }
		if( m_ts[B] ){ m_ts[B]->Unlink( this ); m_ts[B] = NULL; }
		if( m_ts[C] ){ m_ts[C]->Unlink( this ); m_ts[C] = NULL; }
	}

	/// is this triangle surrounded by triangles?
	/*inline bool
	Tri::Boundary() const {
	return  m_vs[A]->Boundary() || m_vs[B]->Boundary() || m_vs[C]->Boundary() || m_vs[D]->Boundary() ;
	}*/

	/// Is the given vertex one of the three that make up the face
	inline bool
		Tri::Contains( const Vertex *v ) const
	{
		return ( m_vs[A] == v ) || ( m_vs[B] == v ) || ( m_vs[C] == v );
	}

	inline bool
		Tri::IsAdj( const Tri *t ) const
	{
		return ( m_ts[A] == t ) || ( m_ts[B] == t ) || ( m_ts[C] == t );
	}

	inline VertexName
		Tri::VName( const Vertex *v) const
	{	
		assert( Contains(v) );
		if ( m_vs[A] == v ) return A;
		if ( m_vs[B] == v ) return B;
		return C; //else
	}

	/// Given two of the vertices in the face, return the name of the third vertex
	inline VertexName     
	Tri::OtherVName( const Vertex* v1, const Vertex* v2  ) const {
		assert( v1 && v2 ); assert( v1 != v2 );
		assert( Contains( v1 ) && Contains( v2 ) );
	
		VertexName a = VName( v1 );
		VertexName b = VName( v2 );
  
		return OtherVN(a,b); 
	}

	inline VertexName
	Tri::VNameAcross( const Tri* t) const {
		assert( IsAdj(t) );
  
		if ( m_ts[A] == t ) return A;
		if ( m_ts[B] == t ) return B;
		return C; //else
	}

	inline Tri*&
	Tri::Across( VertexName i ) {
		assert( i<3 );
		return m_ts[i];
	}

	inline Tri*
	Tri::Across( VertexName i ) const {
		assert( i<3 );
		return m_ts[i];
	}

	/// Return the tet that is across the face defined by the three
	/// Vertices
	inline Tri*&
	Tri::Across( const Vertex *v1, const Vertex *v2 ) {
		return m_ts[OtherVName( v1, v2 )];
	}

	/// Return the face that is across the edge defined by the three
	/// Vertices. This one is const so that others that don't need to 
	/// change Tri can use it.
	inline Tri*
	Tri::Across( const Vertex* v1, const Vertex* v2 ) const {
		return m_ts[OtherVName( v1, v2 )];
	}

	inline Tri*&
	Tri::Across( const Vertex* v ) {
		return m_ts[VName( v )];
	}

	inline Tri*
	Tri::Across( const Vertex* v ) const {
		return m_ts[VName( v )];
	}

	inline Vertex*&
	Tri::OtherV( const Vertex* v1, const Vertex* v2 ) {
		return m_vs[OtherVName( v1, v2 )];
	}

	inline Vertex*
	Tri::OtherV( const Vertex* v1, const Vertex* v2 ) const {
		return m_vs[OtherVName( v1, v2)];
	}

	inline Vertex*&
	Tri::Across( const Tri* t ) {
		return m_vs[VNameAcross( t )];
	}
	inline Vertex*
	Tri::Across( const Tri* t ) const {
		return m_vs[VNameAcross( t )];
	}

	/// Remove the connections for this triangle.
	inline void
	Tri::Unlink( const Tri* t ){
		if( m_ts[A] == t ) m_ts[A] = NULL;
		else if( m_ts[B] == t ) m_ts[B] = NULL;
		else if( m_ts[C] == t ) m_ts[C] = NULL;
		else die();
	}


#include "vertex.hpp"
	inline double
	Tri::Init_area() {
   
		Vertex *a = m_vs[0], *b = m_vs[1], *c = m_vs[2];
    
		Vec2 CA = c->Coord() - a->Coord();
		Vec2 CB = c->Coord() - b->Coord();
    
		m_area = 1.0/2.0 * cross( CA, CB );

		return m_area;
	}

	inline void
	Tri::Init_grad_s() {
		Mat22 F (a()->Coord() - c()->Coord(), b()->Coord() - c()->Coord() );
		F = F.inv();
		grad_s[0] = F.row0();
		grad_s[1] = F.row1();
		grad_s[2] = -1.0* ( grad_s[0] + grad_s[1] );
	
	}

	inline void
	Tri::Calc_grad_s(int j, Vec2& v) const {
		v = grad_s[j];
	}

	inline void 
	Tri::Calc_s(int j, double* s) const {
 		for (int k= 0; k<NumIntPts; ++k){
			s[k] = GaussPoints[k][j];
		}
	}

	inline void
	Tri::Calc_Coord(Vec2* coord) const {
		for (int k = 0; k<NumIntPts; ++k) {
			coord[k] = GaussPoints[k][0] * a()->Coord() + 
				GaussPoints[k][1] * b()->Coord() +
					GaussPoints[k][2] * c()->Coord();
		} 	
	
	}
	
	inline void
	Tri::Calc_u(double* u) const {
		for (int k = 0; k<NumIntPts; ++k) {
			u[k] = GaussPoints[k][0] * a()->u() + 
				GaussPoints[k][1] * b()->u() +
					GaussPoints[k][2] * c()->u();
		} 	
	}
	
	
	inline void 
	Tri::Calc_Phi(double* u_p) const {
		for (int k=0; k< NumIntPts; ++k) {
			u_p[k]= GaussPoints[k][0] * a() -> Phi()+
				GaussPoints[k][1] * b() -> Phi()+
					GaussPoints[k][2] * c() -> Phi();
		}
	}
	
	
	inline void
	Tri::Calc_grad_u(Vec2& grad_u) const {
		grad_u = a()->u() * grad_s[0];
		grad_u += b()->u() * grad_s[1];
		grad_u += c()->u() * grad_s[2];
	}
	
	
inline void
	Tri::Calc_gradPhi(Vec2& grad_u) const {
		grad_u = a()->Phi() * grad_s[0];
		grad_u += b()->Phi() * grad_s[1];
		grad_u += c()->Phi() * grad_s[2];
	}
	

	/*
	inline void 
	Tri::Calc_Ftinv() {
	Mat33 F (a()->Coord() - d()->Coord(), b()->Coord() - d()->Coord(), c()->Coord() - d()->Coord() );
	Ftinv = F.trans().inv();
	// do something.
	}
	*/


	/*
	inline void
	Tri::Calc_grad_u_loc(Vec3& grad_u_loc) const {
	grad_u_loc.x() = a()->u() - d()->u();
	grad_u_loc.y() = b()->u() - d()->u();
	grad_u_loc.z() = c()->u() - d()->u();
	}


	inline void
	Tri::Calc_grad_v_loc(Vec3& grad_v_loc) const {
		
	// do something.
		
	}
 

	inline void
	Tri::Calc_grad_u(Vec3& grad_u) const {

	grad_u = a()->u() * grad_s[0];
	grad_u += b()->u() * grad_s[1];
	grad_u += c()->u() * grad_s[2];
	grad_u += d()->u() * grad_s[3];
	}


	inline void
	Tri::Calc_grad_v(Vec3& grad_v) const {

	//do something
	}

	inline void
	Tri::Calc_u(double* u) const {
	
	
	for (int k = 0; k<NumIntPts; ++k) {
		
	u[k] = GaussPoints[k][0] * a()->u() + 
	GaussPoints[k][1] * b()->u() +
	GaussPoints[k][2] * c()->u() +
	GaussPoints[k][3] * d()->u();
		
	} 	
	
	}

	inline void
	Tri::Calc_v(double* v) const {
	
	for (int k = 0; k<NumIntPts; ++k) {
	v[k] = GaussPoints[k][0] * a()->v() + 
	GaussPoints[k][1] * b()->v() +
	GaussPoints[k][2] * c()->v() +
	GaussPoints[k][3] * d()->v();
	} 	
	
	}


	inline void 
	Tri::Calc_s(int j, double* s) const {
    
	for (int k= 0; k<NumIntPts; ++k){
	s[k] = GaussPoints[k][j];
	}
	
	// switch j.. then fill up w/ loop. 
	// die if not correct j
	
	}

	inline void
	Tri::Init_grad_s() {
	
	Mat33 F (a()->Coord() - d()->Coord(), b()->Coord() - d()->Coord(), c()->Coord() - d()->Coord() );
	F = F.inv();
	grad_s[0] = F.row0();
	grad_s[1] = F.row1();
	grad_s[2] = F.row2();
	grad_s[3] = -1.0* ( grad_s[0] + grad_s[1] + grad_s[2] );
	
	}

	inline void
	Tri::Calc_grad_s(int j, Vec3& v) const {
	v = grad_s[j];
	}

	*/
#endif	/* __TET_HPP__ */
