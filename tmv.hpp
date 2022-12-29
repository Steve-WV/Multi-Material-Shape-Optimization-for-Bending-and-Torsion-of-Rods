#ifndef __TMV_HPP__
#define __TMV_HPP__

class Vec2 {

public:
	
	inline Vec2() {}
	inline Vec2( const Vec2& a ) : m_x(a.m_x), m_y(a.m_y) {}
	inline Vec2( double x, double y ) : m_x(x), m_y(y) {}
	
	inline double& x() { return m_x; }
	inline double x() const { return m_x; }
	inline double& y() { return m_y; }
	inline double y() const { return m_y; }

    inline Vec2& operator=(const Vec2 &rhs) {

       // Only do assignment if RHS is a different object from this.
       if (this != &rhs) {
		   m_x = rhs.m_x; m_y = rhs.m_y;
       }

       return *this;
     }

    inline Vec2& operator+=(const Vec2& rhs)
    {
        m_x += rhs.m_x; m_y += rhs.m_y;
        return *this;
    }
    inline Vec2& operator-=(const Vec2& rhs)
    {
        m_x -= rhs.m_x; m_y -= rhs.m_y;
        return *this;
    }
    inline Vec2& operator/=(const double& rhs)
    {
        m_x /= rhs; m_y /= rhs;
        return *this;
    }
    inline Vec2& operator*=(const double& rhs)
    {
        m_x *= rhs; m_y *= rhs;
        return *this;
    }
	inline const double normsqr(){
	    double norm = sqr(m_x) + sqr(m_y);
	    return norm;
	}
	
private:
	double m_x, m_y;

};

inline const Vec2 operator+(Vec2 lhs, const Vec2& rhs) {
    return lhs += rhs;
}
inline const Vec2 operator-(Vec2 lhs, const Vec2& rhs) {
    return lhs -= rhs;
}
inline const double dot(const Vec2& lhs, const Vec2& rhs) {
	return lhs.x()*rhs.x() + lhs.y()*rhs.y();
}

inline const double cross(const Vec2& lhs, const Vec2& rhs) {
	return lhs.x()*rhs.y() - lhs.y()*rhs.x();
}

inline const Vec2 operator*(const double lambda, const Vec2 &v) {
    Vec2 result;
    result.x() = lambda * v.x();
    result.y() = lambda * v.y();
    
    return result;
}

inline const Vec2 operator*(const Vec2 &v, const double lambda) {
    Vec2 result;
    result.x() = lambda * v.x();
    result.y() = lambda * v.y();
    
    return result;
}


class Mat22 {
public:
    inline Mat22() {}
    inline Mat22(const Mat22& M) :  m_11 (M.m_11), m_12(M.m_12), m_21(M.m_21), m_22(M.m_22) {}
    
    inline Mat22(Vec2 v1, Vec2 v2) : m_11(v1.x()), m_21(v1.y()), m_12(v2.x()), m_22(v2.y()) {}
    
    inline Mat22(double a_11, double a_12, double a_21, double a_22 ) :
    m_11(a_11), m_12(a_12), m_21(a_21), m_22(a_22) {}
    
    inline const Mat22 t() const {
        Mat22 B(m_11, m_21, m_12, m_22);
        return B;
    }
    
    inline Mat22& operator=(const Mat22 &rhs) {
        
        // Only do assignment if RHS is a different object from this.
        if (this != &rhs) {
            m_11 = rhs.m_11; m_12 = rhs.m_12;
            m_21 = rhs.m_21; m_22 = rhs.m_22;
      }
        
        return *this;
    }

    inline const Mat22 inv() const {
        double det = m_11*m_22 - m_21*m_12;
        
        double A_11 =  m_22 / det; 
		double A_12 = -m_12 / det;
        double A_21 = -m_21 / det;
		double A_22 =  m_11 / det;
    
    Mat22 A (A_11, A_12, A_21, A_22 );
    return A;
    
    }
    
    inline const Vec2 row0(){
        return Vec2(m_11, m_12);
    }
    
    inline const Vec2 row1(){
        return Vec2(m_21, m_22);
    }
    

	inline const Vec2 operator*( const Vec2 &v) const {
		Vec2 result;
		result.x() = m_11*v.x() + m_12*v.y();
		result.y() = m_21*v.x() + m_22*v.y();
		return result;
	}
	
private:
	double m_11, m_12, m_21, m_22;

};

#endif