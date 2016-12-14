#include "triangle.h"

#include "CMU462/CMU462.h"
#include "GL/glew.h"

namespace CMU462 { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, vector<size_t>& v) :
    mesh(mesh), v(v) { }
Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {
  
  // TODO: 
  // compute the bounding box of the triangle
  Vector3D p0(mesh->positions[v1].x, mesh->positions[v1].y, mesh->positions[v1].z);
  Vector3D p1(mesh->positions[v2].x, mesh->positions[v2].y, mesh->positions[v2].z);
  Vector3D p2(mesh->positions[v3].x, mesh->positions[v3].y, mesh->positions[v3].z);
  
  BBox tri;
  
  Vector3D min_p(min(p0.x, min(p1.x, p2.x)), min(p0.y, min(p1.y, p2.y)), min(p0.z, min(p1.z, p2.z)) );
  Vector3D max_p(max(p0.x, max(p1.x, p2.x)), max(p0.y, max(p1.y, p2.y)), max(p0.z, max(p1.z, p2.z)) );
  
  tri.min = min_p;
  tri.max = max_p;
  tri.extent = max_p - min_p;
  
  return tri;
}

bool Triangle::intersect(const Ray& r) const {
  
  // TODO: implement ray-triangle intersection
  Vector3D p0(mesh->positions[v1].x, mesh->positions[v1].y, mesh->positions[v1].z);
  Vector3D p1(mesh->positions[v2].x, mesh->positions[v2].y, mesh->positions[v2].z);
  Vector3D p2(mesh->positions[v3].x, mesh->positions[v3].y, mesh->positions[v3].z);
  
  double u = 0.f, v = 0.f, t = 0.f, w = 0.f;
  
  Vector3D s = r.o - p0;
  Vector3D e1 = p1 - p0;
  Vector3D e2 = p2 - p0;
  
  if(fabs(dot(cross(e1,r.d),e2)) < 1e-6)
  {
  	return false;
  
  }else
  {
  	u =  (double)dot(-cross(s,e2),r.d)  / (double)dot(cross(e1,r.d),e2);
    v =  (double)dot( cross(e1,r.d),s)  / (double)dot(cross(e1,r.d),e2);
    t =  (double)dot(-cross(s,e2),  e1) / (double)dot(cross(e1,r.d),e2);
 	w = 1.f - u - v; 

  	
  	if(u >= 0.f && u <= 1.f && v >= 0.f && v <= 1.f && w >= 0.f && w <= 1.f && t >= r.min_t && t <= r.max_t)
  	{
  		r.max_t = t;
  		return true;	
  			
  	}else
  	{
  		return false;
  	}
  }
  
}

bool Triangle::intersect(const Ray& r, Intersection *isect) const {
  
  // TODO: 
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
  		
  		Vector3D p0(mesh->positions[v1].x, mesh->positions[v1].y, mesh->positions[v1].z);
		Vector3D p1(mesh->positions[v2].x, mesh->positions[v2].y, mesh->positions[v2].z);
		Vector3D p2(mesh->positions[v3].x, mesh->positions[v3].y, mesh->positions[v3].z);

		double u = 0.f, v = 0.f, t = 0.f, w = 0.f;

		Vector3D s = r.o - p0;
		Vector3D e1 = p1 - p0;
		Vector3D e2 = p2 - p0;
		
		if(fabs(dot(cross(e1,r.d),e2)) < 1e-6)
  		{
//   			printf("div by zero\n");
  			return false;
  			
  		}else
  		{
  			u =  (double)dot(-cross(s,e2),r.d)  / (double)dot(cross(e1,r.d),e2);
			v =  (double)dot( cross(e1,r.d),s)  / (double)dot(cross(e1,r.d),e2);
			t =  (double)dot(-cross(s,e2),  e1) / (double)dot(cross(e1,r.d),e2);
 			w = 1.f - u - v; 
			
			
			if(u >= 0.f && u <= 1.f && v >= 0.f && v <= 1.f && w >= 0.f && w <= 1.f && t >= r.min_t && t <= r.max_t)
  			{	
  			
  				if(isect == NULL)
  				{
  					isect = new Intersection();
  				}
				
				if(t < isect->t )
				{
					isect->t = t;
					
					auto normal_t = mesh->normals[v1] * (1.f - u - v) + mesh->normals[v2] * u + mesh->normals[v3] * v;
					
					if(dot(normal_t,r.d) < 0.f)
					{
						isect->n = normal_t;
						isect->is_back = false;
						
					}else
					{
						isect->n = -normal_t;
						isect->is_back = true;
					}
					 	
					isect->bsdf = mesh->get_bsdf(); 
		
					isect->primitive = this; 
					
					r.max_t = t;
					
					return true;

				}
					
  			
  			}else
  			{
  				return false;
  			}
  		}
		
  
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CMU462
