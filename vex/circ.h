struct Circle {
    vector center;
		float radius;
}

int[] addCirc(Circle c; int n) {
    int pts[];
    for (int i = 0; i < n; i++) {
        float angle = 2*PI*i/n;
        vector pos = c.center + set(cos(angle), sin(angle), 0) * c.radius;
        append(pts, addpoint(0, pos));
    }
    addprim(0, "polyline", pts);
    return pts;
}
