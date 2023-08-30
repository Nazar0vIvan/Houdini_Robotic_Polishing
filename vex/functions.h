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
    
    float A = determinant(set(1,y2,z2,1,y1,z1,1,y3,z3));
    float B = determinant(set(x2,1,z2,x1,1,z1,x3,1,z3));
    float C = determinant(set(x2,y2,1,x1,y1,1,x3,y3,1));
    float D = -determinant(set(x2,y2,z2,x1,y1,z1,x3,y3,z3));
    
    return set(A,B,C,D);
}

function vector circ(vector pt1; vector pt2; vector pt3){
    float x1 = pt1[0]; float x2 = pt2[0]; float x3 = pt3[0];
    float y1 = pt1[1]; float y2 = pt2[1]; float y3 = pt3[1];
    float xc = -0.5*(y1*(x2*x2+y2*y2-x3*x3-y3*y3)+y2*(x3*x3+y3*y3-x1*x1-y1*y1)+y3*(x1*x1+y1*y1-x2*x2-y2*y2))/(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
    float yc = 0.5*(x1*(x2*x2+y2*y2-x3*x3-y3*y3)+x2*(x3*x3+y3*y3-x1*x1-y1*y1)+x3*(x1*x1+y1*y1-x2*x2-y2*y2))/(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
    float r = sqrt((x1-xc)*(x1-xc) + (y1-yc)*(y1-yc));
    return set(xc, yc, r);
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

function vector poly(vector pt1; vector pt2; vector pt3) {
    return linsolve(pt1[0], pt2[0], pt3[0], pt1[1], pt2[1], pt3[1]);  
}

function vector2 polynorm(vector poly; x0) {
    float a = poly[0]; float b = poly[1]; float c = poly[2];
    return set(-1/2*a*x0+b, a*pow(x0,2)+b*x0*c+x0/2*a*x0+b);
}

function float polyval(a,b,x0) {
    return a*x0+b;
}

function float polyval(a,b,c,x0) {
    return a*pow(x0, 2)+b*x0+c;
}