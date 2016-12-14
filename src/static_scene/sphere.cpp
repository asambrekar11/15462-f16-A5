#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CMU462 { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  	double a = dot(r.d, r.d);
    double b = 2.f * dot(r.o - o, r.d);
    double c = dot(r.o - o, r.o - o) - r2;
    

  if( (b * b - 4.f * a * c) < 0 )//&& a < 1e-6)
  {
  	return false;
  }else
  {
  	double t_1 =  ( -b + sqrt( b * b - 4.f * a * c) ) / (2.f * a );
	double t_2 =  ( -b - sqrt( b * b - 4.f * a * c) ) / (2.f * a );
	t1 = std::min(t_1, t_2);
	t2 = std::max(t_1, t_2);
	return true;
  }
  
  

}

bool Sphere::intersect(const Ray& r) const {

  // TODO:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1,t2;
  bool hit = test(r, t1, t2);
  double min_t = std::min(t1,t2);
  if(hit == true && min_t >= r.min_t && min_t <= r.max_t)
  {
	r.max_t = min_t;
  	
  }else
  {
  	hit = false;
  }
	return hit;
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double t1,t2;
  if(test(r,t1,t2))
  {	
	double min_t = std::min(t1,t2);
	double max_t = std::max(t1,t2);
	
	if(i == NULL)
	{
		i = new Intersection();
	}
	
	if(min_t < i->t && min_t >= r.min_t && min_t <= r.max_t)
	{
		i->t = min_t;
		
		i->n = normal(r.o + min_t * r.d);
		
		i->bsdf = object->get_bsdf(); 
		
		i->primitive = this;
		
		r.max_t = min_t;
		
		return true;
		
	}else if(max_t < i->t && max_t >= r.min_t && max_t <= r.max_t)
	{
		i->t = max_t;
		
		i->n = normal(r.o + max_t * r.d);
		
		i->bsdf = object->get_bsdf(); 
		
		i->primitive = this;
		
		r.max_t = max_t;
		
		return true;
	}
	  	
  }

  	return false;

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CMU462
