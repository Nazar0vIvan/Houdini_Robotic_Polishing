function vector pts2vector(vector pt1; vector pt2) {
    // v = {v1, v1, v3} - направляющий вектор прямой
    float vx = pt2[0] - pt1[0];
    float vy = pt2[1] - pt1[1];
    float vz = pt2[2] - pt1[2];
    return set(vx,vy,vz);
}

function vector2 linsolve(float x1; float x2; float y1; float y2) {
    matrix2 A = set(x1, x2, 1, 1);
    vector2 B = set(y1, y2);
    return invert(A)*B;
}

function vector linsolve(float x1; float x2; float x3; float y1; float y2; float y3) {
    matrix3 A = set(x1*x1, x2*x2, x3*x3, x1, x2, x3, 1, 1, 1);
    vector B = set(y1, y2, y3);
    return invert(A)*B;
}