// Given a time between 0 and 1, evaluates a cubic polynomial with
// the given endpoint and tangent values at the beginning (0) and
// end (1) of the interval.  Optionally, one can request a derivative
// of the spline (0=no derivative, 1=first derivative, 2=2nd derivative).
template <class T>
inline T Spline<T>::cubicSplineUnitInterval(
      const T& position0,
      const T& position1,
      const T& tangent0,
      const T& tangent1,
      double normalizedTime,
      int derivative )
{
	//printf("%f\n",normalizedTime);
	
	double t = normalizedTime;
	double h00, h10, h01, h11;
	
	switch(derivative)
	{
		case 0 :
			h00 = 2.0 * pow(t,3) - 3.0 * pow(t,2) + 1;
			h10 = pow(t,3) - 2.0 * pow(t,2) + t;
			h01 = - 2.0 * pow(t,3) + 3.0 * pow(t,2);
			h11 = pow(t,3) - pow(t,2);
			break;
	
		case 1 : 
			h00 = 6.0 * pow(t,2) - 6.0 * t;
			h10 = 3.0 * pow(t,2) - 4.0 * t + 1;
			h01 = - 6.0 * pow(t,2) + 6.0 * t;
			h11 = 3.0 * pow(t,2) - 2.0 * t;
			break;
	
		case 2 :
			h00 = 12.0 * t - 6.0;
			h10 = 6.0 * t - 4.0;
			h01 = - 12.0 * t + 6.0;
			h11 = 6.0 * t - 2.0;
			break;
			
		//default : 
		//	printf("Invalid derivative\n");
		//	return T();

	}
			
   // TODO IMPLEMENT ME (TASK 1A)
   return position0 * h00 + tangent0 * h10 + position1 * h01 + tangent1 * h11;
}
            
// Returns a state interpolated between the values directly before and after the given time.
template <class T>
inline T Spline<T>::evaluate( double time, int derivative )
{
   // TODO IMPLEMENT ME (TASK 1B)

   
   KnotIter lower_knot = knots.begin();
   KnotIter upper_knot = knots.upper_bound(time);
   
   if( knots.size() < 1 )
   {
   		return T();
   		
   }else if(knots.size() == 1)
   {
   		KnotIter knot_1 = knots.begin();
   		return knot_1->second;
   		
   }else if((*lower_knot).first >= time)
   {
   		KnotIter knot_1 = knots.begin();
   		return knot_1->second;
   		
   }else if(upper_knot == knots.end())
   {
   		KnotIter knot_1 = knots.end();
   		knot_1--;
   		return knot_1->second;
   }else
   {
   		KnotIter k2 = knots.upper_bound(time);
   		KnotIter k1 = prev(k2);
   		KnotIter k3,k0;
   		
   		double t1 = (*k1).first;
   		double t2 = (*k2).first;
   		
   		double t3, t0;
   		
   		T m1,m2;
   		//printf("general\n");
   		
   		if(t1 == (*knots.begin()).first)
   		{
   			t0 = (*k1).first - ((*k2).first - (*k1).first);
   			const T p0 = knots[t1] - (knots[t2] - knots[t1]);
   			m1 = (knots[t2] - p0) / (t2 - t0);

   		}else
   		{
   			k0 = prev(k1);
   			t0 = (*k0).first;
   			m1 = (knots[t2] - knots[t0]) / (t2 - t0);
   		}
   		
   		if(k2 == prev(knots.end()))
   		{
   			t3 = (*k2).first + ( (*k2).first - (*k1).first);
   			const T p3 = knots[t2] + ( knots[t2] - knots[t1]);
   			m2 = (p3 - knots[t1]) / (t3 - t1);
   			
   			
   		}else
   		{
   			k3 = next(k2);
   			t3 = (*k3).first;
   			m2 = (knots[t3] - knots[t1]) / (t3 - t1);
   			
   		}
   		
   		const T p1 = knots[t1];
   		const T p2 = knots[t2]; 
   		m1 = m1 * (t2 - t1);
   		m2 = m2 * (t2 - t1);
   		double normalizedTime = (time - t1) / (t2 - t1);
   		
   		T cubicPolynomial = cubicSplineUnitInterval(p1, p2, m1, m2, normalizedTime, derivative);
   		
   		return cubicPolynomial * pow(t2 - t1 , derivative);
   }
   
}

// Removes the knot closest to the given time,
//    within the given tolerance..
// returns true iff a knot was removed.
template <class T>
inline bool Spline<T>::removeKnot(double time, double tolerance )
{
   // Empty maps have no knots.
   if( knots.size() < 1 )
   {
      return false;
   }

   // Look up the first element > or = to time.
   typename std::map<double, T>::iterator t2_iter = knots.lower_bound(time);
   typename std::map<double, T>::iterator t1_iter;
   t1_iter = t2_iter;
   t1_iter--;

   if( t2_iter == knots.end() )
   {
      t2_iter = t1_iter;
   }

   // Handle tolerance bounds,
   // because we are working with floating point numbers.
   double t1 = (*t1_iter).first;
   double t2 = (*t2_iter).first;

   double d1 = fabs(t1 - time);
   double d2 = fabs(t2 - time);


   if(d1 < tolerance && d1 < d2)
   {
      knots.erase(t1_iter);
      return true;
   }

   if(d2 < tolerance && d2 < d1)
   {
      knots.erase(t2_iter);
      return true;//t2;
   }

   return false;
}

// Sets the value of the spline at a given time (i.e., knot),
// creating a new knot at this time if necessary.
template <class T>
inline void Spline<T>::setValue( double time, T value )
{
   knots[ time ] = value;
}

template <class T>
inline T Spline<T>::operator()( double time )
{
   return evaluate( time );
}
