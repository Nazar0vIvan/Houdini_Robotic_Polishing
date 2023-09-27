#include "$HIP/vex/functions.h"

struct EdgeCirc {
    float xc, yc, zc, R;
		vector2 clipline1;
		vector2 clipline2;
}

function EdgeCirc edgecirc_bytans(vector2 line1; vector2 line2; vector pt1; vector pt2) {  
    float a1 = line1[0]; float b1 = line1[1];
    float a2 = line2[0]; float b2 = line2[1];
    float x1 = pt1[0]; float y1 = pt1[1];
    float x2 = pt2[0]; float y2 = pt2[1];
    float zc = pt1[2];
    
    float xc, yc, r, x2t, y2t, x1t, y1t, alpha1, alpha2;
    vector2 clipline1, clipline2;
    
    // circ1 - circ created from the x1; circ2 - circ created from the x2
    // check position of circ1 touch point (x2t) relative to x2
    xc = (a2*x1*(1+a1*a1)*(a2-a1) + a1*((b2-b1)*(1+a1*a2) + sqrt((a1*a1+1)*(a2*a2+1))*(b1-b2+x1*(a1-a2))))/pow((a1-a2),2);
    x2t = (b1-b2+xc/a2+(x1-xc)/a1+a1*x1)/(a2+1/a2);
    if(abs(x2t) > abs(x2)) {
        // continue create circ1
        r = abs((x1 - xc + a1*a1*x1 + a1*b1 - a1*b2 - a1*a2*xc)/(a1*sqrt(a2*a2+1)));
        yc = a1*x1+b1 - 1/a1*(xc-x1);
        y2t = a2*x2t+b2;
        clipline1 = linsolve(x1, xc, y1, yc);
        clipline2 = linsolve(x2t, xc, y2t, yc); 
    } else {
        // start create circ2
        xc = (a1*x2*(1+a2*a2)*(a1-a2) - a2*((b2-b1)*(1+a1*a2) + sqrt((a1*a1+1)*(a2*a2+1))*(b1-b2+x2*(a1-a2))))/pow((a1-a2),2);
        r = abs((x2 - xc + a2*a2*x2 - a2*b1 + a2*b2 - a1*a2*xc)/(a2*sqrt(a1*a1+1)));
        yc = a2*x2+b2 - 1/a2*(xc-x2);
        x1t = (b2-b1+xc/a1+(x2-xc)/a2+a2*x2)/(a1+1/a1);
        y1t = a1*x1t+b1;
        clipline1 = abs(x1) > abs(x1t) ? linsolve(x1, xc, y1, yc) : linsolve(x1t, xc, y1t, yc);
        clipline2 = linsolve(x2, xc, y2, yc); 
    }
    return EdgeCirc(xc, yc, zc, r, clipline1, clipline2);
}

function int[] addpoints_edgepf(EdgeCirc circ; int n; string edgetype) {
    float k1 = circ.clipline1[0]; float k2 = circ.clipline2[0];
    float z = circ.zc;
    float ba; // begin angle
    float ea; // end angle
    if (abs(k1) < 0.0001) {
        if (edgetype == "le") {
            ba = PI/2;
            ea = k2 > 0 ? PI/2+atan(k2) : 3*PI/2+atan(k2);  
        }
        else {
            ba = k2 > 0 ? -PI+atan(k2) : atan(k2);
            ea = k2 > 0 ? 3*PI/2-atan(k2) : PI/2-atan(k2);
        }
    }
    else {
        if(abs(k2) < 0.0001) {
            if (edgetype == "le") {
                ba = k1 > 0 ? atan(k1) : PI+atan(k1);
                ea = k1 > 0 ? 3*PI/2-atan(k1) : PI/2-atan(k1);
            }
            else {
                ba = -PI/2;
                ea = k1 > 0 ? PI/2+atan(k1) : 3*PI/2+atan(k1);
            }
        }
        else{ // k1 != 0 && k2 != 0
            if(edgetype == "le") {
                ba = k1 > 0 ? atan(k1) : PI + atan(k1);
                ea = k1*k2 > 0 ? PI-atan(k1)+atan(k2) : k1 > 0 ? 2*PI-atan(k1)+atan(k2) : -atan(k1)+atan(k2);
            }
            else {
                ba = k2 < 0 ? atan(k2) : -PI+atan(k2);
                ea = k1*k2 > 0 ? PI+atan(k1)-atan(k2) : k1 < 0 ? 2*PI+atan(k1)-atan(k2) : atan(k1)-atan(k2);  
            }
        }
    }
    float step = ea/n;
    float t = step;
    int ptnums []; int ptnum; 
    for(int i = 0; i < n; i++) {
        float x = circ.xc + circ.R*cos(t + ba);
        float y = circ.yc + circ.R*sin(t + ba);
        ptnum = addpoint(0, set(x,y,z) );
        setpointattrib(0, "blade_part", ptnum, edgetype);
        push(ptnums, ptnum);
        t += step;
     }
     return ptnums;
}