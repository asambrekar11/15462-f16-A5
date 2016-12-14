#include "environment_light.h"
#include "math.h"

namespace CMU462 { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
  // TODO: initialize things here as needed
  double dtheta = PI / (double)(envMap->h);
  double dphi = (2.f * PI) / (double)(envMap->w);
  double sum = 0.f;
  
  for( int i = 0; i < this->envMap->w * this->envMap->h; i++)
  {
  	size_t x = i % envMap->w;
  	size_t y = i / envMap->w;
  	
  	double sintheta = sinf(PI - ( ( (double)y + 0.5f) * (PI)) / (double)this->envMap->h );//Pi - this
  	double p = (double)this->envMap->data[x + envMap->w * y].illum() * sintheta * dphi * dtheta;
	sum += (double)this->envMap->data[x + envMap->w * y].illum() * sintheta * dphi * dtheta;
  	cdf.push_back(sum);
  	pdf.push_back(p);
  }
  
  
  for( int i = 0; i < this->envMap->w * this->envMap->h; i++)
  {
  	size_t x = i % envMap->w;
  	size_t y = i / envMap->w;
  	pdf[x + envMap->w*y] /= sum	;
  	cdf[x + envMap->w*y] /= sum;
  	
  }
	  
}

Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  // TODO: Implement
  *distToLight = INF_D;
  *pdf = 1.f / ( 4 * PI);
  double Xi1 = (double)(std::rand()) / (double)RAND_MAX;
  double Xi2 = (double)(std::rand()) / (double)RAND_MAX;
  
  bool isuniform = true;
  double phi = 0.f;
  double theta = 0.f;
  
  if(!isuniform)
  {
  	//Importance sampling
	int low = 0; int high = envMap->h - 1;
	int mid = low + (high - low) / 2;
	double 	probValue = 0.f;
	int rowNum = 0;
	
	while(low <= high)
	{
		mid = low + (high - low) / 2;
		if(mid == 0)
		{
			probValue = cdf[envMap->w - 1 + envMap->w * (mid)];
		}else if(mid == envMap->h  - 1) 
		{
			probValue = cdf[envMap->w - 1 + envMap->w * (mid)];
		}else
		{
			probValue = cdf[envMap->w - 1 + envMap->w * (mid)];
		}
		if(probValue < Xi1 )
		{
			double nextValue = 0.f;
			if(mid == envMap->h - 1)
			{	
				rowNum = mid;
				break;
			}else
			{
				nextValue = cdf[envMap->w - 1 + envMap->w * (mid + 1)];
			}
			
			if(nextValue >= Xi1 )
			{
				rowNum = mid;
				break;
			}else
			{
				low = mid + 1;
			}
		}else if( probValue > Xi1 )
		{
			double prevValue = 0.f;
			if(mid == 0)
			{
				rowNum = mid;
				break;
			}else 
			{
				prevValue = cdf[envMap->w - 1 + envMap->w * (mid - 1)];
			}
			
			if(prevValue <= Xi1)
			{
				rowNum = mid;
				break;
			}else
			{
				high = mid - 1;
			}
		}else
		{
			rowNum = mid;
			break;
		}
	}
	
	double marg_prob = cdf[envMap->w - 1 + envMap->w * (rowNum)];
	// low = envMap->w * rowNum; high = envMap->w - 1 + envMap->w * rowNum;
	low = 0; high = envMap->w - 1;
	mid = low + (low + high) / 2;
	int columnNum = 0;
	
	while(low <= high)
	{
		
		mid = low + (high - low) / 2;
		if(mid == 0)
		{
			
			probValue = cdf[envMap->w * rowNum];
			probValue /= marg_prob;
			
		}else if(mid == envMap->w - 1)
		{
			probValue = cdf[mid + envMap->w * rowNum];
			probValue /= marg_prob;
			
		}else 
		{
			probValue = cdf[mid + envMap->w * rowNum];
			probValue /= marg_prob;
		}
		
		
		if(probValue < Xi2 )
		{
			double nextValue = 0.f;
			if(mid == envMap->w - 1)
			{
				columnNum = mid;
				break;
			}else if ( mid == 0)
			{
				columnNum = mid;
				break;
			}
			{
				nextValue = cdf[mid + 1 + envMap->w * rowNum];
				nextValue /= marg_prob;
			}
			
			if(nextValue >= Xi2 )
			{
				columnNum = mid;
				break;
			}else
			{
				low = mid + 1;
			}
		}else if( probValue > Xi2 )
		{
			double prevValue = 0.f;
			if(mid ==  0)
			{
				columnNum = mid;
				break;
			}else if(mid == envMap->w - 1)
			{
				columnNum = mid;
				break;
			}else
			{
				prevValue = cdf[mid - 1 + envMap->w * rowNum];
				prevValue /= marg_prob;
			}
			
			if(prevValue <= Xi2)
			{
				columnNum = mid;
				break;
			}else
			{
				high = mid - 1;
			}
		}else
		{
			columnNum = mid;
			break;
		}
		
	}
// 	printf("row: %d, col: %d\n",rowNum, columnNum);
	theta = PI - ((((double)rowNum * (PI)) )/ (double)envMap->h);
	phi =   ((double)columnNum * (2.f * PI)) / (double)envMap->w;
	
  }else
  {
  	//Uniform sampling
  	phi = 2.f * PI * Xi1;
  	theta = PI - acos(1.f - 2.f * Xi2);
  }
  
  Matrix3x3 sampleToWorld;
  sampleToWorld[0] = Vector3D(0,  0,  1);
  sampleToWorld[1] = Vector3D(-1,  0, 0);
  sampleToWorld[2] = Vector3D(0,  -1,  0);
   
  Vector3D sampled_dir = Vector3D(cosf(phi)*sinf(theta), sinf(phi)*sinf(theta), cos(theta) );
  
//   sampled_dir = sampleToWorld*sampled_dir;
  

  
  *wi = sampled_dir.unit();   
   
  return sample_dir(Ray(p, *wi));
}

Spectrum EnvironmentLight::sample_dir(const Ray& r) const {
  // TODO: Implement
  
  double phi = atan2(-r.d[0], r.d[2]) < 0 ? 2 * PI + atan2(-r.d[0], r.d[2]) : atan2(-r.d[0], r.d[2]);
  double theta = acos(r.d[1]);
  
  double u = (phi * (envMap->w - 1.f)) / ( 2.f * PI );
  double v = ((theta) * (envMap->h - 1.f)) / ( PI );
  
  u = std::min((double)envMap->w - 1,(double)u);
  v = std::min((double)envMap->h - 1,(double)v);
  
  int x = floor(u);
  int y = floor(v);
  double xi = floor(u);
  double yi = floor(v);
  double u_ratio = u - xi;
  double v_ratio = v - yi;
  double u_opposite = 1.f - u_ratio;
  double v_opposite = 1.f - v_ratio; 
  int x1 = x == envMap->w - 1 ? x : x+1; //>= envMap->w - 1 ? x : floor(u) + 1;
  int y1 = y == envMap->h - 1 ? y : y+1; //>= envMap->h - 1 ? y : floor(v) + 1;

  return (envMap->data[(x + envMap->w*y)] * u_opposite + envMap->data[(x1 + envMap->w*y)] * u_ratio) * v_opposite + 
		 (envMap->data[(x + envMap->w*(y1))] * u_opposite + envMap->data[(x1 + envMap->w*(y1))] * u_ratio) * v_ratio;  
  
}



} // namespace StaticScene
} // namespace CMU462
