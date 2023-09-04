#include "$HIP/vex/functions.h"

// uv - unit vector
// for poly
function vector get_uvy(vector pt1; vector pt2; vector pt3; vector pt0; int is_cx) {
    float x0 = pt0[0]; float y0 = pt0[1];
    vector2 norm = polynorm(poly(pt1, pt2, pt3), x0); // normal сoeffs at pt0                 
    float x1 = 1;      float y1 = norm[0]*x1+norm[1];
    vector uvy = normalize(set(x1-x0, y1-y0, 0));
    return ((uvy[1] < 0 && is_cx) || (uvy[1] > 0 && !is_cx )) ? uvy *= -1 : uvy;
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
               pt[0],  pt[1],  pt[2],  1.);
}

// circ data: (yc, r)  - circle center; (yt, zt) - touch point of the v line and the circle
// return radius circ center and touch point in attached coor sys 
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

function int[] addpoints_radius_sector(vector4 circ; int npt; matrix T; string blade_part) {
    float yc = circ[0]; float r = circ[1]; float yt = circ[2]; float zt = circ[3];
    vector2 clipline = linsolve(yc,yt,r,zt);
    float k = clipline[0]; float b = clipline[1];
    float begin_angle = k > 0 ? -PI+atan(k) : PI+atan(k);
    float end_angle = PI/2-atan(k);
    float step = end_angle/npt; 
    float t = 0.;
    float y1, z1;
    vector4 ptc;
    int ptnums[], ptnum;
    for(int i = 0; i < npt; i++) {
        y1 = yc + r*cos(t+begin_angle);
        z1 = r*(1 + sin(t+begin_angle));
        ptc = T*set(0., y1, z1, 1.);
        ptnum = addpoint(0, vector(ptc));
        setpointattrib(0, "blade_part", ptnum, "r_"+blade_part);
        //setpointattrib(0, "pf_num", ptnum, sector_num);
        push(ptnums, ptnum);
        t += step;
    }
    return ptnums;
}

// for cx and cv
function int[] addpoints_radius(float r; int npts; int ptnums[]; int first_pf_ptnums[]; matrix3 R; matrix3 RT) {
    string blade_part = point(0, "blade_part", first_pf_ptnums[0]);
    int is_cx = blade_part == "cx" ? 1 : 0; // cx or cv: cx -> 1; cv -> 0

    int n = len(ptnums);
    int j = i == 0 ? 1 : i == n-1 ? n-2 : i;
    vector pt1 = point(0, "P", ptnums[j-1]);
    vector pt2 = point(0, "P", ptnums[j]);
    vector pt3 = point(0, "P", ptnums[j+1]);
    vector pt4 = point(0, "P", first_pf_ptnums[i]);
    
    pt0 = i == 0 ? pt1 : i == n-1 ? pt3 : pt2; // target point    
    vector uvy = R*get_uvy(RT*pt1, RT*pt2, RT*pt3, RT*pt0, is_cx);   
    vector uvx = get_uvx(pt0 + uvy, pt0, pt4);
    vector uvz = cross(uvx, uvy);
    
    matrix T = get_transform(uvx, uvy, uvz, pt0);
    vector4 circ = get_radius_circ(pt0, pt4, r, T);
    return addpoints_radius_sector(circ, npts, T, blade_part);
}

// for edges
function int[] addpoints_radius(float r; int npts; int ptnums[]; int first_pf_ptnums[]; matrix3 R; matrix3 RT; vector ptc) {
    string blade_part = point(0, "blade_part", first_pf_ptnums[0]);
    
    vector pt1 = point(0, "P", ptnums[i]);
    vector pt2 = point(0, "P", first_pf_ptnums[i]); 
    
    vector uvy = get_uvy(pt1, ptc);
    vector uvx = get_uvx(pt1 + uvy, pt1, pt2);
    vector uvz = cross(uvx, uvy);
    
    matrix T = get_transform(uvx, uvy, uvz, pt1);
    vector4 circ = get_radius_circ(pt1, pt2, r, T);

    return addpoints_radius_sector(circ, npts, T, blade_part);
}

// ---------------------------------- BEGIN --------------------------------- //

// all points numbers of the first profile 
int first_pf_ptnums[] = findattribval(0, "point", "pf_num", 0);
int first_pf_cx_ptnums[]; int first_pf_cv_ptnums[];
int first_pf_le_ptnums[]; int first_pf_re_ptnums[];
foreach(int ptnum; first_pf_ptnums) {
    if(point(0, "blade_part", ptnum) == "cx") {
        push(first_pf_cx_ptnums, ptnum); continue;
    }
    if(point(0, "blade_part", ptnum) == "le") {
        push(first_pf_le_ptnums, ptnum); continue;
    }
    if(point(0, "blade_part", ptnum) == "cv") {
        push(first_pf_cv_ptnums, ptnum); continue;
    }
    if(point(0, "blade_part", ptnum) == "re") {
        push(first_pf_re_ptnums, ptnum); continue;
    }
}

// intersection convex points numbers
int i_cx_ptnums[] = findattribval(0, "point", "blade_part", "i_cx");
// intersection concave points numbers
int i_cv_ptnums[] = findattribval(0, "point", "blade_part", "i_cv");
// intersection left edge points numbers
int i_le_ptnums[] = findattribval(0, "point", "blade_part", "i_le");
// intersection right edge points numbers
int i_re_ptnums[] = findattribval(0, "point", "blade_part", "i_re");

float angleY = detail(0, "plane_angle");
matrix3 R = set(cos(angleY), 0., -sin(angleY), 0., 1., 0., sin(angleY), 0., cos(angleY));
matrix3 RT = transpose(R);
float r = 4.0;
int npts = 30;

// cx
for (int i = 0; i < len(i_cx_ptnums)-1; ++i) {
    addpoints_radius(i, r, npts, i_ptnums_cx, first_pf_ptnums_cx, R, RT)
}
// le
vector pt0 = point(0, "P", i_le_ptnums[0]);
vector ptm = point(0, "P", i_le_ptnums[floor(n/2)]);
vector ptn = point(0, "P", i_le_ptnums[n-1]);
vector ptc = R*set(vector(circ(RT*pt0, RT*ptm, RT*ptn))); // circ center in 0 coor sys


for (int i = 0; i < len(i_le_ptnums)-1; ++i) {
    addpoints_radius(i, r, npts, i_ptnums_cx, first_pf_ptnums_cx, R, RT, ptc);
}
// cv
for (int i = 0; i < len(i_cv_ptnums); ++i) {
    addpoints_radius(i, r, npts, i_ptnums_cx, first_pf_ptnums_cx, R, RT)
}
// re
vector pt0 = point(0, "P", i_re_ptnums[0]);
vector ptm = point(0, "P", i_re_ptnums[floor(n/2)]);
vector ptn = point(0, "P", i_re_ptnums[n-1]);
vector ptc = R*set(vector(circ(RT*pt0, RT*ptm, RT*ptn))); // circ center in 0 coor sys
for (int i = 0; i < len(i_re_ptnums)-1; ++i) {
    addpoints_radius(i, r, npts, i_ptnums_cx, first_pf_ptnums_cx, R, RT, ptc);
}

/*
int r_ptnums_cx[] = addpoints_radius_cxcv(r, npts, i_ptnums_cx, first_pf_ptnums_cx, R, RT);
int r_ptnums_le[] = addpoints_radius_edge(r, npts, i_ptnums_le, first_pf_ptnums_le, R, RT);
int r_ptnums_cv[] = addpoints_radius_cxcv(r, npts, i_ptnums_cv, first_pf_ptnums_cv, R, RT);
int r_ptnums_re[] = addpoints_radius_edge(r, npts, i_ptnums_re, first_pf_ptnums_re, R, RT);

setdetailattrib(0, "last_r_num", len(i_ptnums_cx)+len(i_ptnums_cv)+len(i_ptnums_le)+len(i_ptnums_re));
*/

// ---------------------------------- END --------------------------------- //