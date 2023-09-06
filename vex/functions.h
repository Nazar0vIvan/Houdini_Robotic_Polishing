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

function vector4 circ(vector pt1; vector pt2; vector pt3){
    float x1 = pt1[0]; float x2 = pt2[0]; float x3 = pt3[0];
    float y1 = pt1[1]; float y2 = pt2[1]; float y3 = pt3[1];
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
    float a = poly[0]; float b = poly[1]; float c = poly[2];
    return set(-1/(2*a*x0+b), a*x0*x0+b*x0+c + x0/(2*a*x0+b));
}

function float polyval(float a; float b; float x0) {
    return a*x0+b;
}

function float polyval(float a; float b; float c; float x0) {
    return a*pow(x0,2)+b*x0+c;
}

function vector4 lsqr(int ptnums[]) {
	vector pt;
	float X2 = 0.; float XY = 0.; float X = 0.;
	float Y2 = 0.; float Y = 0.;
	float XZ = 0.; float YZ = 0.; float Z = 0.;
	foreach(int ptnum; ptnums) {
			pt = point(0, "P", ptnum);
			
			X2 += pt[0]*pt[0]; Y2 += pt[1]*pt[1];
			XY += pt[0]*pt[1];
			X  += pt[0]; Y  += pt[1];
			
			XZ += pt[0]*pt[2];
			YZ += pt[1]*pt[2];
			Z  += pt[2];   
	}
	matrix3 U = set(X2, XY, X, XY, Y2, Y, X, Y, len(ptnums));
	vector V = set(XZ, YZ, Z);
	vector P = invert(U)*V;

	float C = 1/sqrt(P[0]*P[0]+P[1]*P[1]+1);
	float A = -P[0]*C;
	float B = -P[1]*C;
	float D = -P[2]*C;
	
	return set(A,B,C,D);
}

function int[] linspace(int x1; int x2) {
    int res[];
    if(x2 < x1) return res;
    for (int i = x1; i <= x2; ++i) push(res, i);    
    return res;
}

function void removepoints(int geohandle; int point_numbers[]) {
    foreach(int ptnum; point_numbers) removepoint(geohandle, ptnum);
}

function matrix3 invMatrix(matrix3 m) {
    float m00 = getcomp(m, 0, 0); float m01 = getcomp(m, 0, 1); float m02 = getcomp(m, 0, 2);
    float m10 = getcomp(m, 1, 0); float m11 = getcomp(m, 1, 1); float m12 = getcomp(m, 1, 2);
    float m20 = getcomp(m, 2, 0); float m21 = getcomp(m, 2, 1); float m22 = getcomp(m, 2, 2);
    
    float invdet = 1/determinant(m);
    
    float im00 =  (m11*m22-m21*m12)*invdet;
    float im10 = -(m01*m22-m02*m21)*invdet;
    float im20 =  (m01*m12-m02*m11)*invdet;
    float im01 = -(m10*m22-m12*m20)*invdet;
    float im11 =  (m00*m22-m02*m20)*invdet;
    float im21 = -(m00*m12-m10*m02)*invdet;
    float im02 =  (m10*m21-m20*m11)*invdet;
    float im12 = -(m00*m21-m20*m01)*invdet;
    float im22 =  (m00*m11-m10*m01)*invdet;
    
    return set(im00, im10, im20, im01, im11, im21, im02, im12, im22);
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

// пересечение прямой (точки u1, u2) и конуса, заданного образующей kg*x+bg
function vector intersec_cone_line(vector u1; vector u2; float kg; float bg) {
    vector vec = u2 - u1;
    float vx = vec[0]; float vy = vec[1]; float vz = vec[2];
    float x1 = u1[0]; float y1 = u1[1]; float z1 = u1[2];
     
    float a = pow(vy/kg,2) - pow(vx,2) + pow(vz/kg,2);
    float b = 2*vy*y1/(pow(kg,2)) - 2*vx*(x1+bg/kg) + 2*vz*z1/(pow(kg,2));
    float c = pow(y1/kg,2) - pow(x1+bg/kg,2) + pow(z1/kg,2);
    
    float t1, t2;
    solvequadratic(a,b,c,t1,t2);
    
    return set(x1 + vx*t1, y1 + vy*t1, z1 + vz*t1);
}