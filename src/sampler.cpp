#include "sampler.h"
#include <random>

namespace CMU462 {

  // Uniform Sampler2D Implementation //

  Vector2D UniformGridSampler2D::get_sample() const {

    // TODO:
    // Implement uniform 2D grid sampler
	
	double a = (double)(std::rand()) / (double)RAND_MAX;
    double b = (double)(std::rand()) / (double)RAND_MAX;
// 	std::default_random_engine generator;
//     std::uniform_real_distribution<double> distr(0.0,1.0);

// std::cout << distr(generator) << '\n';

    return Vector2D(a, b);

  }

  // Uniform Hemisphere Sampler3D Implementation //

  Vector3D UniformHemisphereSampler3D::get_sample() const {

    double Xi1 = (double)(std::rand()) / RAND_MAX;
    double Xi2 = (double)(std::rand()) / RAND_MAX;

    double theta = acos(Xi1);
    double phi = 2.0 * PI * Xi2;

    double xs = sinf(theta) * cosf(phi);
    double ys = sinf(theta) * sinf(phi);
    double zs = cosf(theta);

    return Vector3D(xs, ys, zs);

  }

  Vector3D CosineWeightedHemisphereSampler3D::get_sample() const {
    float f;
    return get_sample(&f);
  }

  Vector3D CosineWeightedHemisphereSampler3D::get_sample(float *pdf) const {
    // You may implement this, but don't have to.
    
	double Xi1 = (double)(std::rand()) / RAND_MAX;
    double Xi2 = (double)(std::rand()) / RAND_MAX;

    double costheta = sqrt(Xi1);
    double sintheta = sqrt(1.0f - Xi1);
    double phi = 2.0 * PI * Xi2;
    double theta = acos(sqrt(Xi1));
    
    *pdf = cosf(theta) / PI;
	
    double xs = sinf(theta) * cosf(phi);
    double ys = sinf(theta) * sinf(phi);
    double zs = cosf(theta);
    return Vector3D(xs, ys, zs);
  }


} // namespace CMU462
