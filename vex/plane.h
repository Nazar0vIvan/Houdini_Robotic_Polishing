struct Plane {
  vector4 coeffs;
	
  function float distance(vector pt) {
    return A*pt[0] + B*pt[1] + C*pt[2] + D;
  }
	
  function int is_above(vector pt) {
    float aa = -this.coeffs[0]/this.coeffs[2],
          bb = -this.coeffs[1]/this.coeffs[2],
          cc = -this.coeffs[3]/this.coeffs[2];
          return (pt[2] > aa*pt[0] + bb*pt[1] + cc) ? 1 : 0;
  }
  function vector intersec(vector pt1; vector pt2) {
    float A = this.coeffs[0],
          B = this.coeffs[1],
          C = this.coeffs[2],
          D = this.coeffs[3];
    float x1 = pt1[0],
          y1 = pt1[1],
          z1 = pt1[2];
		vector dir_vec = pt1 - pt2;
    float vx = dir_vec[0], 
          vy = dir_vec[1],
          vz = dir_vec[2];
    float t = -(A*x1+B*y1+C*z1+D)/(A*vx+B*vy+C*vz);
    return set(x1 + vx*t, y1 + vy*t, z1 + vz*t);
  }
}
