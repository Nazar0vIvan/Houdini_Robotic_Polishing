#include "$HIP/vex/functions.h"

// for cx and cv
function int[] addpoints_radius_cxcv(int ptnum; float r; int npts; int ptnums[]; int first_pf_ptnums[]; matrix3 R) {
    int n = len(ptnums); int i;
    if (ptnum == 0) i = 1;
    if (ptnum == n-1) i = n-2;
    vector pt0, pt1, pt2, pt3, pt4;
    pt1 = point(0, "P", ptnums[i-1]);
    pt2 = point(0, "P", ptnums[i]);
    pt3 = point(0, "P", ptnums[i+1]);
    pt4 = point(0, "P", first_pf_ptnums[ptnum]); 
    
    matrix3 RT = transpose(R);

    vector pt0 = ptnum == 0 ? pt1 : ptnum == 1 ? pt3 : pt2 // target point
     // y axis of attached coor sys
    vector uvy = get_uvy(RT*pt1, RT*pt2, RT*pt3, RT*pt0);
    if ((dir == 1 && uvy[1] < 0) || (dir == 0 && uvy[1] > 0)) uvy*=-1;
    // x axis of attached coor sys 
    vector uvx = get_uvx(pt0 + uvy, pt0, pt4);
    // z axis of attached coor sys 
    vector uvz = cross(uvx, uvy);
    // transformation matrix
    matrix T = get_transform(uvx, uvy, uvz, pt_target);
    
    matrix TI = invert(T);
    vector4 pt0T = TI*set(pt0[0], pt0[1], pt0[2], 1.0);
    vector4 pt4T = TI*set(pt4[0], pt4[1], pt4[2], 1.0);
    linsolve(pt0T[1], pt4T[1], pt0T[2], pt4T[2])
    vector4 circ = get_radius_circ(pt0, pt4, r, T);

    push(ptnums, addpoints_radius_sector(circ, T, npts));
    
    return ptnums;
}

// for edges
function int[] addpoints_radius_edge(int ptnum; float r; int npts; int ptnums[]; int first_pf_ptnums[]; matrix3 R) {
    int n = len(ptnums);
    matrix3 RT = transpose(R);
    vector pt0 = point(0, "P", ptnums[0]);
    vector ptm = point(0, "P", ptnums[floor(n/2)]);
    vector ptn = point(0, "P", ptnums[n-1]);
    
    vector edge_circ = circ(RT*pt0, RT*ptm, RT*ptn);
    vector ptc = R*set(edge_circ[0], edge_circ[1], 0.);
    vector pt1, pt2;
    
    pt1 = point(0, "P", ptnums[i]);
    pt2 = point(0, "P", first_pf_ptnums[i]); 
    // y axis of attached coor sys
    vector uvy = get_uvy(pt1, ptc);
    // x axis of attached coor sys 
    vector uvx = get_uvx(pt1 + uvy, pt1, pt2);
    // z axis of attached cs
    vector uvz = cross(uvx, uvy);
    
    matrix T = get_transform(uvx, uvy, uvz, pt1);
    vector4 circ = get_radius_circ(pt1, pt2, r, T);
    
    push(ptnums, addpoints_radius_sector(circ, T, npts));
    
    return ptnums;
}

// uv - unit vector
// for poly
function vector get_uvy(vector pt1; vector pt2; vector pt3; vector pt0) {
    vector2 norm = polynorm(poly(pt1, pt2, pt3), pt0) // normal сoeffs at pt0                 
    float x0 = pt0[0]; float y0 = pt0[1];
    float x1 = 1;      float y1 = norm[0]*x1+norm[1];
    return normalize(set(x1-x0, y1-y0, 0));
}

// for circ
function vector get_uvy(vector pt1; vector pt2) {
    return normalize(pt1 - pt2);
}

function vector get_uvx(vector pt1; vector pt2; vector pt3) {
    vector4 vplane = plane(pt1, pt2, pt3);
    return normalize(set(vplane[0], vplane[1], vplane[2]));
}

function matrix get_transform(vector uvx; vector uvy; vector uvz; vector pt) {
    return set(uvx[0], uvx[1], uvx[2], 0., 
               uvy[0], uvy[1], uvy[2], 0., 
               uvz[0], uvz[1], uvz[2], 0., 
               pt2[0], pt2[1], pt2[2], 1.);
}

// circ data: (yc, r)  - circle center; (yt, zt) - touch point of the line and the circle
function vector4 get_radius_circ(vector pt1; vector pt2; float r; matrix T) {        
    matrix TI = invert(T);
    vector4 pt1T = TI*set(pt1[0], pt1[1], pt1[2], 1.0);
    vector4 pt2T = TI*set(pt2[0], pt2[1], pt2[2], 1.0);
    vector2 tangent = linsolve(pt1T[1], pt2T[1], pt1T[2], pt2T[2]);
    float k = tangent[0]; float b = tangent[1];

    // abscissa of the circ center
    float yc_max = (-b+r*(1+sqrt(k*k+1)))/k;
    float yc_min = (-b+r*(1-sqrt(k*k+1)))/k;
    
    // abscissa of the touch point
    float yt_max = 1/k*(-b+r*(1+1/(sqrt(k*k+1))));
    float yt_min = 1/k*(-b+r*(1-1/(sqrt(k*k+1))));
    
    float yc = k > 0 ? yc_max : yc_min;
    float yt = k > 0 ? yt_max : yt_min;
    float zt = k*yt+b;
    
    return set(yc,r,yt,zt);
}

function int[] addpoints_radius_sector(vector circ; int npt, matrix T) {
    float yc = circ[0]; float r = circ[1]; float yt = circ[1]; float zt = circ[2];
    vector2 clipline = linsolve(yc,yt,r,zt);
    float k = clipline[0]; float b = clipline[1];
    float begin_angle = k > 0 ? -PI+atan(k) : PI+atan(k);
    float end_angle = PI/2-atan(k);
    float step = end_angle/npt; 
    float t = 0.;
    float y1, z1;
    vector4 ptc;
    int ptnums[];
    for(int i = 0; i < npt; i++) {
        y1 = yc + r*cos(t+begin_angle);
        z1 = r*(1 + sin(t+begin_angle));
        ptc = T*set(0., y1, z1, 1.);
        push(ptnums, addpoint(0, vector(ptc)));
        t += step;
    }  
    return ptnums;
}

// BEGIN //

// all points numbers of the first profile 
int first_pf_ptnums[] = findattribval(0, "point", "pf_num", 0);
int first_pf_ptnums_cx[]; int first_pf_ptnums_cv[];
int first_pf_ptnums_le[]; int first_pf_ptnums_re[];
foreach(int ptnum; first_pf_ptnums) {
    if(point(0, "pf_part", ptnum) == "cx") {
        push(first_pf_ptnums_cx, ptnum); continue;
    }
    if(point(0, "pf_part", ptnum) == "cv") {
        push(first_pf_ptnums_cv, ptnum); continue;
    }
    if(point(0, "pf_part", ptnum) == "le") {
        push(first_pf_ptnums_le, ptnum); continue;
    }
    if(point(0, "pf_part", ptnum) == "re") {
        push(first_pf_ptnums_re, ptnum); continue;
    }
}

// all intersection points numbers
int i_ptnums[] = findattribval(0, "point", "pf_num", -2);
// intersection convex points numbers
int i_ptnums_cx[] = findattribval(0, "point", "pf_part", "i_cx");
// intersection concave points numbers
int i_ptnums_cv[] = findattribval(0, "point", "pf_part", "i_cv");
// intersection left edge points numbers
int i_ptnums_le[] = findattribval(0, "point", "pf_part", "i_le");
// intersection right edge points numbers
int i_ptnums_re[] = findattribval(0, "point", "pf_part", "i_re");

float angleY = detail(0, "plane_angle");
matrix3 R = set(cos(angleY), 0., -sin(angleY), 0., 1., 0., sin(angleY), 0., cos(angleY));
float r = 4.0;
int npts = 50;
int ptnums[];

for(int i = 0; i < len(i_ptnums_cx); ++i) {
    ptnums = addpoints_radius_cvcx(i, r, npts, i_ptnums_cx, first_pf_ptnums_cx, R);
    setprimattrib(0, "type", addprim(0, "polyline", ptnums), "radius");
}
for(int i = 0; i < len(i_ptnums_le); ++i) {
    ptnums = addpoints_radius_edge(i, r, npts, i_ptnums_le, first_pf_ptnums_le, R)
    setprimattrib(0, "type", addprim(0, "polyline", ptnums), "radius");
}
for(int i = 0; i < len(i_ptnums_cv); ++i) {
    ptnums = addpoints_radius_cvcx(i, r, npts, i_ptnums_cv, first_pf_ptnums_cv, R)
    setprimattrib(0, "type", addprim(0, "polyline", ptnums), "radius");
}
for(int i = 0; i < len(i_ptnums_re); ++i) {
    ptnums = addpoints_radius_edge(i, r, npts, i_ptnums_re, first_pf_ptnums_re, R)
    setprimattrib(0, "type", addprim(0, "polyline", ptnums), "radius");
}

