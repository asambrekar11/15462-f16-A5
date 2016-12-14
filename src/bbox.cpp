#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CMU462 {

  bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

    // TODO:
    // Implement ray - bounding box intersection test
    // If the ray intersected the bouding box within the range given by
    // t0, t1, update t0 and t1 with the new intersection times.

	
	double tx_min = r.sign[0] == 0 ? (min.x - r.o.x) * r.inv_d.x : (max.x - r.o.x) * r.inv_d.x;
	double tx_max = r.sign[0] == 0 ? (max.x - r.o.x) * r.inv_d.x : (min.x - r.o.x) * r.inv_d.x;
	double ty_min = r.sign[1] == 0 ? (min.y - r.o.y) * r.inv_d.y : (max.y - r.o.y) * r.inv_d.y;
	double ty_max = r.sign[1] == 0 ? (max.y - r.o.y) * r.inv_d.y : (min.y - r.o.y) * r.inv_d.y;
	double tz_min = r.sign[2] == 0 ? (min.z - r.o.z) * r.inv_d.z : (max.z - r.o.z) * r.inv_d.z;
	double tz_max = r.sign[2] == 0 ? (max.z - r.o.z) * r.inv_d.z : (min.z - r.o.z) * r.inv_d.z;
	
	double t_min = 0.f, t_max = 0.f;
	
	t_min = std::max(tx_min, std::max(ty_min, tz_min));
	t_max = std::min(tx_max, std::min(ty_max, tz_max));
	
	if(t_min <= t_max )
	{
		if( t_min >= t0 && t_min <= t1 )
		{
			t0 = t_min;
		}
		if( t_max <= t1 && t_max >= t0  )
		{
			t1 = t_max;
		} 
		
		return true;
	}else 
	{
		return false;
	}
	
    

  }

  void BBox::draw(Color c) const {

    glColor4f(c.r, c.g, c.b, c.a);

    // top
    glBegin(GL_LINE_STRIP);
    glVertex3d(max.x, max.y, max.z);
    glVertex3d(max.x, max.y, min.z);
    glVertex3d(min.x, max.y, min.z);
    glVertex3d(min.x, max.y, max.z);
    glVertex3d(max.x, max.y, max.z);
    glEnd();

    // bottom
    glBegin(GL_LINE_STRIP);
    glVertex3d(min.x, min.y, min.z);
    glVertex3d(min.x, min.y, max.z);
    glVertex3d(max.x, min.y, max.z);
    glVertex3d(max.x, min.y, min.z);
    glVertex3d(min.x, min.y, min.z);
    glEnd();

    // side
    glBegin(GL_LINES);
    glVertex3d(max.x, max.y, max.z);
    glVertex3d(max.x, min.y, max.z);
    glVertex3d(max.x, max.y, min.z);
    glVertex3d(max.x, min.y, min.z);
    glVertex3d(min.x, max.y, min.z);
    glVertex3d(min.x, min.y, min.z);
    glVertex3d(min.x, max.y, max.z);
    glVertex3d(min.x, min.y, max.z);
    glEnd();

  }

  std::ostream& operator<<(std::ostream& os, const BBox& b) {
    return os << "BBOX(" << b.min << ", " << b.max << ")";
  }

} // namespace CMU462

// double tx_min = -INF_D, tx_max = INF_D;
// 	if(fabs(r.d.x) < 1e-6)
// 	{
// 		if(r.d.x > 0.f )
// 		{
// 			tx_min = (min.x - r.o.x) / r.d.x;
// 			tx_max = (max.x - r.o.x) / r.d.x;
// 		}else
// 		{
// 			tx_min = (max.x - r.o.x) / r.d.x;
// 			tx_max = (min.x - r.o.x) / r.d.x;
// 		}
// 		
// 	}
// 	
// 	double ty_min = -INF_D, ty_max = INF_D; 
// 	if(fabs(r.d.y) < 1e-6)
// 	{
// 		if(r.d.y > 0.f )
// 		{
// 			ty_min = (min.y - r.o.y) / r.d.y;
// 			ty_max = (max.y - r.o.y) / r.d.y;
// 		}else
// 		{
// 			ty_min = (max.y - r.o.y) / r.d.y;
// 			ty_max = (min.y - r.o.y) / r.d.y;
// 		}
// 		
// 	}
// 	
// 	double tz_min = -INF_D, tz_max = INF_D; 
// 	if(fabs(r.d.z) < 1e-6)
// 	{
// 		if(r.d.z > 0.f )
// 		{
// 			tz_min = (min.z - r.o.z) / r.d.z;
// 			tz_max = (max.z - r.o.z) / r.d.z;
// 		}else
// 		{
// 			tz_min = (max.z - r.o.z) / r.d.z;
// 			tz_max = (min.z - r.o.z) / r.d.z;			
// 		}
// 		
// 	}
// 	
// 	if( (tx_min > ty_max) || (ty_min > tx_max) && (tx_min > tz_max) || (tz_min > tx_max) && (ty_min > tz_max) || (tz_min > ty_max) )
// 	{
// 		return false
// 	}else
// 	{
// 		if( (tx_min < ty_max ) || (ty_min < tx_max) )
// 		{
// 			double t_min,t_max;
// 			
// 			t_min = std::max(tx_min, ty_min);
// 			t_max = std::min(tx_max, ty_max );
// 			if(t_min > t0)
// 			{
// 				t0 = t_min;
// 			}
// 		
// 			if(t_max < t1)
// 			{
// 				t1 = t_max;
// 			}		
// 		}else if( (tx_min < tz_max ) || (tz_min < tx_max) )
// 		{
// 			t_min = std::max(tx_min, tz_min);
// 			t_max = std::min(tx_max, tz_max );
// 			if(t_min > t0)
// 			{
// 				t0 = t_min;
// 			}
// 		
// 			if(t_max < t1)
// 			{
// 				t1 = t_max;
// 			}	
// 		}else
// 		{
// 			t_min = std::max(ty_min, tz_min);
// 			t_max = std::min(ty_max, tz_max );
// 			if(t_min > t0)
// 			{
// 				t0 = t_min;
// 			}
// 		
// 			if(t_max < t1)
// 			{
// 				t1 = t_max;
// 			}	
// 		}
// 		return true;
// 	}
// 	
