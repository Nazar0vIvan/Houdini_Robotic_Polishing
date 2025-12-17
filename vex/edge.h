#include "linalg.h"

struct Circle {
    vector center;
		float radius;
}

function int[] addCirc(Circle circ; int n; int inOnlyPoints) {
    int pt_nums[];
    vector pc = circ.center;
		float r = circ.radius;
		vector pt;
		float angle;
		for (int i = 0; i < n; i++) {
        angle = 2*PI*i/n;
        pt = pc + set(cos(angle), sin(angle), 0) * r;
        append(pt_nums, addpoint(0, pt));
    }
    if (inOnlyPoints) addprim(0, "polyline", pt_nums);
    return pt_nums;
}

vector pickTangent(
    vector pte1; vector pte2;
    vector t1; vector t2;
    int    xDir;  // +1 = pick larger X, -1 = pick smaller X
    int    yDir   // +1 = pick larger Y, -1 = pick smaller Y
)
{
    if (isOppositePoints(pte1, pte2, t1, t2)) {
        return (t1.x * xDir > t2.x * xDir) ? t1 : t2;
    } else {
        return (t1.y * yDir > t2.y * yDir) ? t1 : t2;
    }
}

function vector[] pointCircleTangents(vector pt; Circle circ)
{               
    float  r = circ.radius;
    vector pc = circ.center;
    vector v    = pt - pc;
    float  d2   = dot(v, v),
           a    = r*r/d2,
           h    = r*sqrt(d2 - r*r)/d2;
    vector vp = set(-v.y, v.x, 0);
    return array(pc + a*v + h*vp, pc + a*v - h*vp);
}

function vector[] tangentsPoints(vector edge_pts[]; Circle circ_le; Circle circ_re)
{
    vector out[]; resize(out, 4); 

    for (int i = 0; i < 4; i++) {
        // decide left/right circle & segment endpoints
        int    isLeft = (i < 2);
        Circle circ   = isLeft ? circ_le : circ_re;
        vector pte1   = vector(edge_pts[isLeft ? 0 : 2]);
        vector pte2   = vector(edge_pts[isLeft ? 1 : 3]);
        vector pte     = vector(edge_pts[i]);

        // compute the two raw tangent candidates
        vector tangs[] = pointCircleTangents(pte, circ);
        vector  cand[] = set(tangs); // cand[0], cand[1]
        vector  t1 = cand[0], t2 = cand[1];

        // flags: X by circle side, Y by even/odd index
        int xDir = isLeft ? -1 : +1;         // left=min-X, right=max-X
        int yDir = ((i % 2) == 0) ? +1 : -1; // even=max-Y, odd=min-Y

        // pick and store the correct tangent
        out[i] = pickTangent(pte1, pte2, t1, t2, xDir, yDir);
    }
    return out;
}

function int[] addArc(vector pt1; vector pt2; Circle circ; int n; string side; int isClockwise) {
    int pt_num, pt_nums[];
    vector pc = circ.center;
    float  r = circ.radius;
    float  a1 = atan2(pt1.y - pc.y, pt1.x - pc.x);
    float  a2 = atan2(pt2.y - pc.y, pt2.x - pc.x);
    int isRight = side == "re" ? 1 : 0;

    float d = a2 - a1;
    if (d < 0) d += 2 * M_PI;

    vector mid = pc + r * set(cos(a1 + d * 0.5), sin(a1 + d * 0.5), 0);
    if ((cross(pt2 - pt1, mid - pt1).z > 0) != (isRight != 0)) {
        d -= 2 * M_PI;
    }

    float step = d / (n + 1);
    float ang;
    append(pt_nums, addpoint(geoself(), isClockwise ? pt1 : pt2));
    for (int i = 1; i <= n; i++) {
        ang = a1 + step * (isClockwise ? i : n + 1 - i);
        pt_num = addpoint(geoself(), pc + r * set(cos(ang), sin(ang), 0));
        push(pt_nums, pt_num);
    }
    append(pt_nums, addpoint(geoself(), isClockwise ? pt2 : pt1)); 
    return pt_nums;
}

