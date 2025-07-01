struct Circle {
    vector center;
		float radius;
}

function int[] addCirc(Circle circ; int n) {
    int pt_nums[];
    vector pc = circ.center;
		float r = circ.radius;
		vector pt;
		float angle;
		for (int i = 0; i < n; i++) {
        angle = 2*PI*i/n;
        pt = pc + set(cos(angle), sin(angle), 0) * r;
        append(pt_nums, addpoint(0, pt));
    }
    if (inOnlyPoints) addprim(0, "polyline", pts);
    return pt_nums;
}

