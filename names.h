#ifndef	__NAMES_H__
#define	__NAMES_H__

enum VertexName { A = 0, B = 1, C = 2, VNONE = 3 };

// turn one type of number (rightmost letter) into another kind of number
// (left letter)
inline VertexName   OtherVN( const VertexName v1, const VertexName v2){
	assert( v1 != v2 ); return VertexName(3 - v1 - v2);
}

inline VertexName   NextVN( const VertexName v ) { return VertexName((v+1)%3); } 
inline VertexName   PrevVN( const VertexName v ) { return VertexName((v+2)%3); } 

#endif	/* __NAMES_H__ */
