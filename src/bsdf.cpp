#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CMU462 {

  void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
  }

  // Diffuse BSDF //

  Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    return albedo * (1.0 / PI);
  }

  Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  	
  	*wi = sampler.get_sample(pdf);
    return f(wo, *wi);
  }

  // Mirror BSDF //

  Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    
    Vector3D dir(0.f, 0.f, 0.f);
    reflect(wo, &dir);
    if( 1.f - (dot(dir.unit(), wi.unit())) < 1e-6)
    return reflectance * (1.f / wi[2]) ;
    
    return Spectrum();
  }

  Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

    // TODO:
    // Implement MirrorBSDF
    reflect(wo, wi);
    *pdf = 1.f ;
    return f(wo, *wi);
  }

  // Glossy BSDF //

  /*
     Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
     return Spectrum();
     }

     Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
   *pdf = 1.0f;
   return reflect(wo, wi, reflectance);
   }
   */

  // Refraction BSDF //

  Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    printf("Refraction!\n");
    return Spectrum();
  }

  Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

    // TODO:
    // Implement RefractionBSDF

    return Spectrum();
  }

  // Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
   
   double ni = 1.0, nt = ior;
    if(wo[2] < 0.0)
    {
    	std::swap(ni, nt);
    }
    Vector3D dir_refract(0.0,0.0,0.0);
   	bool isRefract = refract(wo, &dir_refract, (nt/ni));
   	{
   		double f_r = 1.0, r_perp = 0.0, r_parl = 0.0;
   		double cos_theta_i = 0.0;
		double cos_theta_t = 0.0;
  
  		cos_theta_i = fabs(wo[2]);
		cos_theta_t = !isRefract ? 0.0 : fabs(dir_refract.z);
		

		r_parl = (ni*cos_theta_i - nt*cos_theta_t) / (ni*cos_theta_i + nt*cos_theta_t);
		r_perp = (ni*cos_theta_t - nt*cos_theta_i) / (ni*cos_theta_t + nt*cos_theta_i);
				
		f_r = std::min((double)1.f, std::max((double)0.0, 0.5 * (pow(r_parl,2) + pow(r_perp,2))));
   		
   		Vector3D dir_reflect(0.0,0.0,0.0);
   	
   		reflect(wo,&dir_reflect);
   		
   		if(isRefract && (1.0 - (dot(dir_refract, wi))) < 1e-4)
   		{
   			return transmittance * pow((nt/ ni),2)*(1.0 - f_r)*(1.0 / fabs(wi.z));
   		}
   		
   		if( (1.0 - (dot(dir_reflect, wi))) < 1e-4)
   		{
   			return reflectance * (1.0 / fabs(wi.z)) * f_r;
   			
   		}
   		

   			return Spectrum();			
   		
   	}
   

  }

  Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

    // TODO:
    // Compute Fresnel coefficient and either reflect or refract based on it.
    double ni = 1.0, nt = ior;
    if(wo[2] < 0.0)
    {
    	std::swap(ni, nt);
    }
   	bool isRefract = refract(wo, wi, (nt/ni));

   	{
   		double f_r = 1.0, r_perp = 0.0, r_parl = 0.0;
   		double cos_theta_i = 0.0;
		double cos_theta_t = 0.0;
  
  		cos_theta_i = fabs(wo[2]);
		cos_theta_t = !isRefract ? 0.f : fabs(wi->z);

		r_parl = (nt*cos_theta_t - ni*cos_theta_i) / (nt*cos_theta_t + ni*cos_theta_i);
		r_perp = (ni*cos_theta_t - nt*cos_theta_i) / (ni*cos_theta_t + nt*cos_theta_i);
		
		f_r = std::min((double)1.0 ,std::max(0.0, 0.5 * (pow(r_parl,2) + pow(r_perp,2))));
		
		float rand_prob = (float)rand() / (float) RAND_MAX;
   		if(rand_prob <= f_r)
   		{
   			reflect(wo, wi);
			*pdf = f_r;
			
   		}else
   		{
			*pdf = 1.0 - f_r;	
   		}
   		
   		return f(wo,*wi);
   		
	}

  }
  



  void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

    // TODO:
    // Implement reflection of wo about normal (0,0,1) and store result in wi.
	wi->x = -wo.x;
	wi->y = -wo.y;
	wi->z =  wo.z;
	
	(*wi).normalize();

  }

  bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

    // TODO:
    // Use Snell's Law to refract wo surface and store result ray in wi.
    // Return false if refraction does not occur due to total internal reflection
    // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
    // ray entering the surface through vacuum.
  
  
  	float sign = wo[2] / fabs(wo[2]);
	double cos_theta_i = fabs(wo[2]);
	double cos_theta_t = 0.0;
	double cos_theta_t2 = 0.0;
	cos_theta_t2 = 1.0 - pow((1.0/ior),2) * (1.0 - pow(cos_theta_i,2));
	if(cos_theta_t2 > 0.0)
	{
		cos_theta_t = sqrt(cos_theta_t2);
		wi->x = -wo.x / ior;
		wi->y = -wo.y / ior;
		wi->z = -sign * cos_theta_t;
    
    	(*wi).normalize();
		return true;
	}
	
	return false;
	


  }

  // Emission BSDF //

  Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    return Spectrum();
  }

  Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    *wi  = sampler.get_sample(pdf);
    return Spectrum();
  }

} // namespace CMU462
