#include <vector>
#include <istream>
#include "boost/tokenizer.hpp"

#include "utils.hpp"
#include "mesh.hpp"
#include "tri.hpp"
#include "vertex.hpp"
#include "edge.hpp"


/// This function sets the adjacency pointers for all triangles.
/// It does so by iterating over triangles and calling a helper function
/// that deals with the edges given by the respective vertex pairs.
void Mesh::BuildTriangleTopology()
{
    // build edge list
    EdgeSet es;
    TriIt ti = TriBegin(), te = TriEnd();
    for( ; ti != te ; ++ti ) {
        Tri* t = *ti;
        AddAndOrMatch( es, t->a(), t->b(), t );
        AddAndOrMatch( es, t->b(), t->c(), t );
        AddAndOrMatch( es, t->c(), t->a(), t );
    }
    
    // now we can check for boundary vertices
    SetBoundary( es );
    
}

void Mesh::SetBoundary( Mesh::EdgeSet& es ) {
	
	VertexSet vs;
    
    auto ei = es.begin(), ee = es.end();
    for (; ei!=ee; ++ei ) {
        Edge e = *ei;
        if ( e.Boundary() ) {
            (*ei).a()->MarkBoundary();
            vs.insert((*ei).a());
			(*ei).b()->MarkBoundary();
			vs.insert((*ei).b());
		}
    }
	
	m_bv = VertexCt(vs.begin(), vs.end() );
    
}

/// Helper function to connect triangles.
void Mesh::AddAndOrMatch( EdgeSet& es, Vertex *a, Vertex *b, Tri* tn )
{
    // For a given edge which has a triangle adjacent try to find the triangle
    // on the other side. The edge set es knows whether the edge we're trying to
    // insert is already contained in the set.
    Edge en = Edge( a, b, tn );
    std::pair<EdgeSet::iterator, bool> i = es.insert( en );
    if( !i.second ){ // the edge was already in the set.
        
        // This means that there is a triangle on the other side,
        // and we can set up the adjacency information.
        Edge eo = *i.first; Tri *to = eo.t();
        
        // first we test if the other tet already has some connection
        // set up. This would be fatal!
        assert( !(to->Across( eo.a(), eo.b()) ) );
        
        // Some more testing to make sure that the two triangles
        // have the same endpoints
        
        // XXX maybe check this somehow... idk.
        //assert( ( en.a() == eo.b() ) && ( en.b() == eo.a() ) );
        // edges must originate from opposite triangles
        assert( to != tn );
        
        // set neighbor information for triangles
        to->Across( eo.a(), eo.b() ) = tn;
        tn->Across( eo.a(), eo.b() ) = to;
        
    }
}

#include <fstream>
void Mesh::Load(std::string fname) {
	
	std::ifstream is;
	is.open(fname.c_str());

	std::string line;
	while( std::getline(is, line) ){
		
		boost::char_separator<char> sep(" ");
		boost::tokenizer< boost::char_separator<char> > tok(line, sep);
		boost::tokenizer< boost::char_separator<char> >::iterator ti = tok.begin(), te = tok.end();
		
		if (ti==te) continue; //empty line.
		
       	if( *ti == "v" ){
            
            double x = std::stod(*(++ti));
            double y = std::stod(*(++ti));
            Vertex* v = new Vertex( x, y );
            assert( v );
            addVertex( v );
            
        } 
		if( *ti == "f" ) {
            
			Tri *t = new Tri( V(std::stoi(*(++ti))-1),
                             V(std::stoi(*(++ti))-1),
                             V(std::stoi(*(++ti))-1) );
            assert( t );
            // get anything the Triangle may want (currently empty)
            addTri( t );
        }
    }
    
    // all the initialization code associated with a fully loaded Mesh, i.e.
    // build triangle and vertex topology, set boundary and set vertex valences
    Init();
}


// Destructor
Mesh::~Mesh( void )
{
    // be a good citizen and clean up after yourself...
    for_each( VertexBegin(), VertexEnd(), Delete<Vertex>() );
    for_each( TriBegin(), TriEnd(), Delete<Tri>() );
}






