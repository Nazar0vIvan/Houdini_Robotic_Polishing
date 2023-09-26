struct Plane {
	float coeffs[];
	
	function int is_above(vector pt) {
		float aa = -this.coeffs[0]/this.coeffs[2],
					bb = -this.coeffs[1]/this.coeffs[2],
					cc = -this.coeffs[3]/this.coeffs[2];
		return (pt[2] > aa*pt[0] + bb*pt[1] + cc) ? 1 : 0;
	}
}