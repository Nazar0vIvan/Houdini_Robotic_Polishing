#include "$HIP/vex/functions.h"

// for edges
function int[] addpoints_radius_edge(float r; int ptnums[]; string curvetype; int first_pf_pt_nums[]; matrix3 R) {
    int n = len(ptnums);
    matrix3 RT = transpose(R);
    vector pt0 = point(0, "P", ptnums[0]);
    vector ptm = point(0, "P", ptnums[floor(n/2)]);
    vector ptn = point(0, "P", ptnums[n-1]);
    
    vector edge_circ = circ(RT*pt0, RT*ptm, RT*ptn);
    vector ptc = set(R*edge_circ[0], R*edge_circ[1], pt0[2]);
    vector pt1, pt2;
    
    for(int i = 0; i < n; ++i) {
        pt1 = point(0, "P", ptnums[i]);
        pt2 = point(0, "P", first_pf_pt_nums[i]); 
        // y axis of attached coor sys
        vector uvy = get_uvy(pt1, ptc);
        // x axis of attached coor sys 
        vector uvx = get_uvx(pt1 + uvy, pt1, pt2);
        // z axis of attached cs
        vector uvz = cross(uvx, uvy);
        
        matrix T = get_transform(uvx, uvy, uvz, pt1);
        
        matrix TI = invert(T);
        vector4 pt1T = TI*set(pt1[0], pt1[1], pt1[2], 1.0);
        vector4 pt2T = TI*set(pt2[0], pt2[1], pt2[2], 1.0);
        vector4 circ = get_radius_circ(linsolve(pt1T[1], pt2T[1], pt1T[2], pt2T[2]), r);
        
        float yc = circ[0]; float yt = circ[2]; float zt = circ[3];
        push(ptnums, addpoints_radius_sector(yc, r, linsolve(yc,yt,r,zt), T, 50));
    }
    
    return ptnums;
}

// for cx and cv
function int[] addpoints_radius_cxcv(float r; int ptnums[]; int dir; int first_pf_pt_nums[]; matrix3 R) {
    vector pt1, pt2, pt3, pt4;
    int n = len(ptnums);
    matrix3 RT = transpose(R);
    for(int i = 1; i < n-1; ++i) {
        pt1 = point(0, "P", ptnums[i-1]);
        pt2 = point(0, "P", ptnums[i]);
        pt3 = point(0, "P", ptnums[i+1]);
        pt4 = point(0, "P", first_pf_pt_nums[i]); 
        // y axis of attached coor sys
        vector uvy = get_uvy(RT*pt1, RT*pt2, RT*pt3);
        if ((dir == 1 && uvy[1] < 0) || (dir == 0 && uvy[1] > 0)) uvy*=-1;
        // x axis of attached coor sys 
        vector uvx = get_uvx(pt2 + uvy, pt2, pt4);
        // z axis of attached cs
        vector uvz = cross(uvx, uvy);
        // transformation matrix
        matrix T = get_transform(uvx, uvy, uvz, pt2);
        
        matrix TI = invert(T);
        vector4 pt2T = TI*set(pt2[0], pt2[1], pt2[2], 1.0);
        vector4 pt4T = TI*set(pt4[0], pt4[1], pt4[2], 1.0);
        vector4 circ = get_radius_circ(linsolve(pt2T[1], pt4T[1], pt2T[2], pt4T[2]), r);

        float yc = circ[0]; float yt = circ[2]; float zt = circ[3];
        push(ptnums, addpoints_radius_sector(yc, r, linsolve(yc,yt,r,zt), T, 50));
    }
    return ptnums;
}

// uv - unit vector
function vector get_uvy(vector pt1; vector pt2; vector pt3; int dir, int pt0) {
    vector poly2 = linsolve(pt1[0], pt2[0], pt3[0], pt1[1], pt2[1], pt3[1]);           
    float a = poly2[0]; float b = poly2[1]; float c = poly2[0];                  
    float x0 = pt2[0];
    float y0 = pt2[1];
    float z0 = pt2[2];
    float x1 = 1;
    float y1 = a*x0*x0+b*x0*c - 1/(2*a*x0+b)*(x1-x0);
    return normalize(set(x1-x0, y1-y0, 0));
}

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

// circ data: (yc, r)  - circle center; (yt, zt) - intersection point of the line and the circle
function vector4 get_radius_circ(vector2 tangent; float r) {        
    float k = tangent[0]; float b = tangent[1];
    
    float yc_max = (-b+r*(1+sqrt(k*k+1)))/k;
    float yc_min = (-b+r*(1-sqrt(k*k+1)))/k;
    
    float yt_max = 1/k*(-b+r*(1+1/(sqrt(k*k+1))));
    float yt_min = 1/k*(-b+r*(1-1/(sqrt(k*k+1))));
    
    float yc = k > 0 ? yc_max : yc_min;
    float yt = k > 0 ? yt_max : yt_min;
    float zt = k*yt+b;
    
    return set(yc,r,yt,zt);
}

function int[] addpoints_radius_sector(float yc; float r; vector2 clipline; matrix T; int npt) {
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
    setprimattrib(0, "type", addprim(0, "polyline", ptnums), "radius");
    return ptnums;
}



// BEGIN

// get all points numbers of the first profile 
int first_pf_pt_nums[] = findattribval(0, "point", "pf_num", 0);
// get all penultimate prim points numbers of convex and concave
int first_pf_pt_nums_cx[]; int first_pf_pt_nums_cv[];
int first_pf_pt_nums_le[]; int first_pf_pt_nums_re[];
foreach(int ptnum; first_pf_pt_nums) {
    if(point(0, "pf_part", ptnum) == "cx") {
        push(first_pf_pt_nums_cx, ptnum); continue;
    }
    if(point(0, "pf_part", ptnum) == "cv") {
        push(first_pf_pt_nums_cv, ptnum); continue;
    }
    if(point(0, "pf_part", ptnum) == "le") {
        push(first_pf_pt_nums_le, ptnum); continue;
    }
    if(point(0, "pf_part", ptnum) == "re") {
        push(first_pf_pt_nums_re, ptnum); continue;
    }
}
// intersection convex points numbers
int i_cx_pt_nums[] = findattribval(0, "point", "pf_part", "i_cx");
// intersection concave points numbers
int i_cv_pt_nums[] = findattribval(0, "point", "pf_part", "i_cv");
// intersection left edge points numbers
int i_le_pt_nums[] = findattribval(0, "point", "pf_part", "i_le");
// intersection right edge points numbers
int i_re_pt_nums[] = findattribval(0, "point", "pf_part", "i_re");

float angleY = detail(0, "plane_angle");
matrix3 R = set(cos(angleY), 0., -sin(angleY), 0., 1., 0., sin(angleY), 0., cos(angleY));
float r = 4.0;

int cx[] = addpoints_radius_cxcv(r, i_cx_pt_nums, "cx", first_pf_pt_nums_cx, R);
int le[] = addpoints_radius_edge(r, i_le_pt_nums, "le", first_pf_pt_nums_le, R);
int cv[] = addpoints_radius_cxcv(r, i_cv_pt_nums, "cv", first_pf_pt_nums_cv, R);
int re[] = addpoints_radius_edge(r, i_re_pt_nums, "re", first_pf_pt_nums_re, R);

