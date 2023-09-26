struct RuledSurf {
	vector a1, a2, b1, b2;
	
	function vector val(float u; float v) {
		return (1-u)*(1-v)*this.a1 + u*(1-v)*this.a2 + (1-u)*v*this.b1 + u*v*this.b2;
	}
	
	function vector du(float v) {
		return (this.a1 - this.a2)*(v-1) - (this.b1 - this.b2)*v;
	}
	
	function vector dv(float u) {
		return (this.a1 - this.b1)*(u-1) - (this.a2 - this.b2)*u;
	}
}