function vector4 plane(vector norm; vector pt) {
    float A = norm[0];
    float B = norm[1];
    float C = norm[2];
    float D = -A*pt[0]-B*pt[1]-C*pt[2];
    
    return set(A,B,C,D);
}

function vector4 plane(vector pt1; vector pt2; vector pt3) {
    float x1 = pt1[0]; float y1 = pt1[1]; float z1 = pt1[2];
		float x2 = pt2[0]; float y2 = pt2[1]; float z2 = pt2[2];
    float x3 = pt3[0]; float y3 = pt3[1]; float z3 = pt3[2];
    
    float A = determinant(set(1,y2,z2,1,y1,z1,1,y2,z2));
    float B = determinant(set(x2,1,z2,x1,1,z1,x3,1,z3));
    float C = determinant(set(x2,y2,1,x1,y1,1,x3,y3,1));
    float D = -determinant(set(x2,y2,z2,x1,y1,z1,x3,y3,z3));
    
    return set(A,B,C,D);
}

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