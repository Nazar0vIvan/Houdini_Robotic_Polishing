function float[] plane(vector norm; vector pt) {
		float A = norm[0],
				  B = norm[1],
			    C = norm[2],
					D = -A*pt[0]-B*pt[1]-C*pt[2];
    return array(A,B,C,D);
}

function float[] plane(vector pt1; vector pt2; vector pt3) {
    float x1 = pt1[0],
					y1 = pt1[1],
					z1 = pt1[2],
					x2 = pt2[0],
					y2 = pt2[1],
					z2 = pt2[2],
					x3 = pt3[0],
					y3 = pt3[1],
					z3 = pt3[2];
    float A = determinant(set(1,y2,z2,1,y1,z1,1,y3,z3)),
					B = determinant(set(x2,1,z2,x1,1,z1,x3,1,z3)),
					C = determinant(set(x2,y2,1,x1,y1,1,x3,y3,1)),
					D = -determinant(set(x2,y2,z2,x1,y1,z1,x3,y3,z3));
    return array(A,B,C,D);
}

function vector4 circ(vector pt1; vector pt2; vector pt3){
    float x1 = pt1[0],
				  x2 = pt2[0], 
					x3 = pt3[0];
    float y1 = pt1[1],
					y2 = pt2[1], 
					y3 = pt3[1];
    float xc = -0.5*(y1*(x2*x2+y2*y2-x3*x3-y3*y3)+y2*(x3*x3+y3*y3-x1*x1-y1*y1)+y3*(x1*x1+y1*y1-x2*x2-y2*y2))/(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
    float yc = 0.5*(x1*(x2*x2+y2*y2-x3*x3-y3*y3)+x2*(x3*x3+y3*y3-x1*x1-y1*y1)+x3*(x1*x1+y1*y1-x2*x2-y2*y2))/(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
    float r = sqrt((x1-xc)*(x1-xc) + (y1-yc)*(y1-yc));
    float zc = pt1[2];
    return set(xc, yc, zc, r);
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

function vector2 polynorm(vector poly; float x0) {
    float a = poly[0],
					b = poly[1],
					c = poly[2];
    return set(-1/(2*a*x0+b), a*x0*x0+b*x0+c + x0/(2*a*x0+b));
}

function float polyval(float a; float b; float x0) {
    return a*x0+b;
}

function float polyval(float a; float b; float c; float x0) {
    return a*pow(x0,2)+b*x0+c;
}

function int[] linspace(int x1; int x2; int is_end) {
    int res[];
    int k = x2 > x1 ? 1 : -1;
    for (int i = 0; i < abs(x2-x1); ++i) push(res, x1+k*i);
    if(is_end) push(res, x2);
    return res;
}

function void removepoints(int geohandle; int point_numbers[]) {
    foreach(int ptnum; point_numbers) removepoint(geohandle, ptnum);
}

function void printDict(dict dic) {
    foreach(string key; keys(dic)) {
        vector2 value = dic[key];
        printf("%s: %i\n", key, value);
    }
}

function void printMatrix(matrix3 m) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            printf("%8.3f ", getcomp(m, i, j));
        }
    printf("\n");
    }
}

// intersection of the line (points u1, u2) and cone, defined by generator (line) kg*x+bg
function vector intersec_cone_line(vector u1; vector u2; float kg; float bg) {
    vector vec = u2 - u1;
    float vx = vec[0],
			    vy = vec[1],
					vz = vec[2];
    float x1 = u1[0],
				  y1 = u1[1],
					z1 = u1[2];
    float a = pow(vy/kg,2) - pow(vx,2) + pow(vz/kg,2);
    float b = 2*vy*y1/(pow(kg,2)) - 2*vx*(x1+bg/kg) + 2*vz*z1/(pow(kg,2));
    float c = pow(y1/kg,2) - pow(x1+bg/kg,2) + pow(z1/kg,2);
    
    float t1, t2;
    solvequadratic(a,b,c,t1,t2);
    
    return set(x1 + vx*t1, y1 + vy*t1, z1 + vz*t1);
}
