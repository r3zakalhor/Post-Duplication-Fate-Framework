// ClassifyDuplicationFates.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <ostream>
#include <map>
#include <cmath>

#include <limits>
#include <iomanip>

using namespace std;
int counter = 0;

struct point2D {

    int64_t x, y;
public:
    point2D()
    {
        x = 0;
        y = 0;
    };
    point2D(int64_t x1, int64_t y1) {
        x = x1;
        y = y1;
        
    }
    point2D& operator=(point2D const& obj)
    {
        x = obj.x;
        y = obj.y;
        return *this;
    }
    bool operator == (const point2D& Ref) const
    {
        return(this->x == Ref.x && this->y == Ref.y);
    }
};
// note, might as well take p as const-reference
std::ostream& operator<<(std::ostream& o, const point2D& p)
{
    o << "x: " << p.x << "\ty: " << p.y << std::endl;
    return o;
}
const int   N = 10; // clipped (new) polygon size

// check if a point is on the LEFT side of an edge
bool inside(point2D p, point2D p1, point2D p2)
{
    return (p2.y - p1.y) * p.x + (p1.x - p2.x) * p.y + (p2.x * p1.y - p1.x * p2.y) < 0;
}

// calculate intersection point
point2D intersection(point2D cp1, point2D cp2, point2D s, point2D e)
{
    point2D dc = { cp1.x - cp2.x, cp1.y - cp2.y };
    point2D dp = { s.x - e.x, s.y - e.y };

    int64_t n1 = cp1.x * cp2.y - cp1.y * cp2.x;
    int64_t n2 = s.x * e.y - s.y * e.x;
    int64_t n3 = (dc.x * dp.y - dc.y * dp.x);

    return { (n1 * dp.x - n2 * dc.x) / n3, (n1 * dp.y - n2 * dc.y) / n3 };
}

// Sutherland-Hodgman clipping
void SutherlandHodgman(point2D* subjectPolygon, int& subjectPolygonSize, point2D* clipPolygon, int& clipPolygonSize, point2D(&newPolygon)[N], int& newPolygonSize)
{

    int o, k;
    for (o = 0; o < subjectPolygonSize; o++)
        if (subjectPolygon[o].y != 0)
            break;
    for (k = 0; k < clipPolygonSize; k++)
        if (clipPolygon[k].y != 0)
            break;
    bool negativeflag = false;

    //see end of function for a note
    if (o < subjectPolygonSize && k < clipPolygonSize)  //ML, fix the stack around g2d being corrupted
    {
        if (subjectPolygon[o].y < 0 && clipPolygon[k].y < 0) {
            subjectPolygon[o].y *= -1;
            clipPolygon[k].y *= -1;
            negativeflag = true;
        }
    }
    
    
    point2D cp1, cp2, s, e, inputPolygon[N];
    // copy subject polygon to new polygon and set its size
    for (int i = 0; i < subjectPolygonSize; i++)
        newPolygon[i] = subjectPolygon[i];

    newPolygonSize = subjectPolygonSize;

    bool check[3];// = { false, false, false };
    for (int i = 0; i < 3; i++) {
        check[i] = false;
        for (int j = 0; j < 3; j++) {
            if (subjectPolygon[i] == clipPolygon[j])
                check[i] = true;
        }
    }
    if (check[0] && check[1] && check[2]) {
        if (negativeflag) {
            subjectPolygon[2].y *= -1;
            clipPolygon[2].y *= -1;
            for (int i = 0; i < newPolygonSize; i++)
                if (newPolygon[i].y != 0)
                    newPolygon[i].y *= -1;
        }
        return;
    }



    for (int j = 0; j < clipPolygonSize; j++)
    {
        // copy new polygon to input polygon & set counter to 0
        for (int k = 0; k < newPolygonSize; k++) { inputPolygon[k] = newPolygon[k]; }
        int counter = 0;

        // get clipping polygon edge
        cp1 = clipPolygon[j];
        cp2 = clipPolygon[(j + 1) % clipPolygonSize];

        for (int i = 0; i < newPolygonSize; i++)
        {
            // get subject polygon edge
            s = inputPolygon[i];
            e = inputPolygon[(i + 1) % newPolygonSize];

            // Case 1: Both vertices are inside:
            // Only the second vertex is added to the output list
            if (inside(s, cp1, cp2) && inside(e, cp1, cp2))
                newPolygon[counter++] = e;

            // Case 2: First vertex is outside while second one is inside:
            // Both the point of intersection of the edge with the clip boundary
            // and the second vertex are added to the output list
            else if (!inside(s, cp1, cp2) && inside(e, cp1, cp2))
            {
                newPolygon[counter++] = intersection(cp1, cp2, s, e);
                newPolygon[counter++] = e;
            }

            // Case 3: First vertex is inside while second one is outside:
            // Only the point of intersection of the edge with the clip boundary
            // is added to the output list
            else if (inside(s, cp1, cp2) && !inside(e, cp1, cp2))
                newPolygon[counter++] = intersection(cp1, cp2, s, e);

            // Case 4: Both vertices are outside
            else if (!inside(s, cp1, cp2) && !inside(e, cp1, cp2))
            {
                // No vertices are added to the output list
            }
        }
        // set new polygon size
        newPolygonSize = counter;
    }

    //ML: not sure this works as intended.  if negativeFlag is true, it is assumed that o=2 and k=2.
    //What happens if this is not the case?  Also what happens if the passed polygons are not triangles but have a negative y?
    if (negativeflag) {
        subjectPolygon[2].y *= -1;
        clipPolygon[2].y *= -1;
        for (int i = 0; i < newPolygonSize; i++)
            if (newPolygon[i].y != 0)
                newPolygon[i].y *= -1;
    }
}


// find equation of a line from two points (y = mx + z)
void findequation(point2D a, point2D b, long double mz[])
{
    mz[0] = (long double)(b.y - a.y) / (long double)(b.x - a.x);
    mz[1] = (long double)(a.y) - (long double)(mz[0] * a.x);
}

// sort vertices of triangle based on x
point2D* sorttrianglevertices(point2D a[3])
{
    point2D* sorted = a;//new point2D[3];
    /*sorted[0] = a[0];
    sorted[1] = a[1];
    sorted[2] = a[2];*/
    point2D temp, temp1;
    if (a[0].x < a[1].x && a[0].x < a[2].x)
    {
        sorted[0] = a[0];
        if (a[1].x < a[2].x) {
            sorted[1] = a[1];
            sorted[2] = a[2];
        }
        else {
            temp = a[1];
            sorted[1] = a[2];
            sorted[2] = temp;
        }
    }
    else if (a[1].x < a[0].x && a[1].x < a[2].x)
    {
        temp = a[0];
        sorted[0] = a[1];
        if (a[0].x < a[2].x) {
            sorted[1] = temp;
            sorted[2] = a[2];
        }
        else {
            sorted[1] = a[2];
            sorted[2] = temp;
        }
    }
    else if (a[2].x < a[0].x && a[2].x < a[1].x)
    {
        temp = a[0];
        temp1 = a[1];
        sorted[0] = a[2];
        if (a[0].x < a[1].x) {
            sorted[1] = temp;
            sorted[2] = temp1;
        }
        else {
            sorted[1] = temp1;
            sorted[2] = temp;
        }
    }
    //for (int i = 0; i < 3; i++)
        //cout << "(" << sorted[i].x << ", " << sorted[i].y << ")" << endl;
    return sorted;
}

// calculate a + b when a and b have intersection ( in this case we have new polygon with 7 vertices)
void calculate_a_plus_b_polygon(point2D* a, point2D* b, point2D(&a_plus_b_polygon)[6])
{
    int size = 0;
    long double mz[2];
    //point2D a_plus_b_polygon[7]{};
    point2D* sorted_a = sorttrianglevertices(a);
    point2D* sorted_b = sorttrianglevertices(b);
    point2D temp[3]{};
    if (sorted_a[0].x > sorted_b[0].x)
    {
        for (int i = 0; i < 3; i++) {
            temp[i] = sorted_a[i];
            sorted_a[i] = sorted_b[i];
            sorted_b[i] = temp[i];
        }
    }
    //cout << sorted_a[0] << sorted_a[1] << sorted_a[2] << "**" << endl;
    //cout << sorted_b[0] << sorted_b[1] << sorted_b[2];
    //point2D* sorted_intersection_a_b = sorttrianglevertices(intersection_a_b);
    a_plus_b_polygon[0] = sorted_a[0];

    if (sorted_b[0].x >= sorted_a[1].x) {
        a_plus_b_polygon[1] = sorted_a[1];
        if (sorted_b[0].x >= sorted_a[2].x) {
            a_plus_b_polygon[2] = sorted_a[2];
            a_plus_b_polygon[3] = sorted_b[0];
            a_plus_b_polygon[4] = sorted_b[1];
            a_plus_b_polygon[5] = sorted_b[2];
        }
        else {
            findequation(sorted_a[1], sorted_a[2], mz);
            //cout << "m ,z" << mz[0] << mz[1] << endl;
            a_plus_b_polygon[2] = { sorted_b[0].x , (int64_t)(mz[0] * (long double)sorted_b[0].x + mz[1]) };
            if (sorted_b[1].x >= sorted_a[2].x) {
                findequation(sorted_b[0], sorted_b[1], mz);
                //cout << "m ,z" << mz[0] << mz[1] << endl;
                a_plus_b_polygon[3] = { sorted_a[2].x , (int64_t)(mz[0] * (long double)sorted_a[2].x + mz[1]) };
                a_plus_b_polygon[4] = sorted_b[1];
                a_plus_b_polygon[5] = sorted_b[2];
            }
            else {
                findequation(sorted_a[1], sorted_a[2], mz);
                //cout << "m ,z" << mz[0] << mz[1] << endl;
                a_plus_b_polygon[3] = { sorted_b[1].x , (int64_t)(mz[0] * (long double) sorted_b[1].x + mz[1] + sorted_b[1].y) };
                if (sorted_b[2].x >= sorted_a[2].x) {
                    findequation(sorted_b[1], sorted_b[2], mz);
                    //cout << "m ,z" << mz[0] << mz[1] << endl;
                    a_plus_b_polygon[4] = { sorted_a[2].x , (int64_t)(mz[0] * (long double) sorted_a[2].x + mz[1]) };
                    a_plus_b_polygon[5] = sorted_b[2];
                }
                else {
                    findequation(sorted_a[1], sorted_a[2], mz);
                    //cout << "m ,z" << mz[0] << mz[1] << endl;
                    a_plus_b_polygon[4] = { sorted_b[2].x , int64_t(mz[0] * (long double) sorted_b[2].x + mz[1]) };
                    a_plus_b_polygon[5] = sorted_a[2];
                }
            }
        }
    }
    else {
        findequation(sorted_a[0], sorted_a[1], mz);
        //cout << "m ,z" << mz[0] << mz[1] << endl;
        a_plus_b_polygon[1] = { sorted_b[0].x , int64_t(mz[0] * (long double)sorted_b[0].x + mz[1]) };
        if (sorted_b[1].x >= sorted_a[1].x) {
            findequation(sorted_b[0], sorted_b[1], mz);
            //cout << "m ,z" << mz[0] << mz[1] << endl;
            a_plus_b_polygon[2] = { sorted_a[1].x , int64_t(mz[0] * (long double)sorted_a[1].x + mz[1] + sorted_a[1].y) };
            if (sorted_b[1].x >= sorted_a[2].x) {
                findequation(sorted_b[0], sorted_b[1], mz);
                //cout << "m ,z" << mz[0] << mz[1] << endl;
                a_plus_b_polygon[3] = { sorted_a[2].x , int64_t(mz[0] * (long double)sorted_a[2].x + mz[1]) };
                a_plus_b_polygon[4] = sorted_b[1];
                a_plus_b_polygon[5] = sorted_b[2];
            }
            else {
                findequation(sorted_a[1], sorted_a[2], mz);
                a_plus_b_polygon[3] = { sorted_b[1].x , int64_t(mz[0] * (long double)sorted_b[1].x + mz[1] + sorted_b[1].y) };
                if (sorted_b[2].x >= sorted_a[2].x) {
                    findequation(sorted_b[1], sorted_b[2], mz);
                    //cout << "m ,z" << mz[0] << mz[1] << endl;
                    a_plus_b_polygon[4] = { sorted_a[2].x , int64_t(mz[0] * (long double)sorted_a[2].x + mz[1]) };
                    a_plus_b_polygon[5] = sorted_b[2];
                }
                else {
                    findequation(sorted_a[1], sorted_a[2], mz);
                    //cout << "m ,z" << mz[0] << mz[1] << endl;
                    a_plus_b_polygon[4] = { sorted_b[2].x , int64_t(mz[0] * (long double)sorted_b[2].x + mz[1]) };
                    a_plus_b_polygon[5] = sorted_a[2];
                }
            }
        }
        else {
            findequation(sorted_a[0], sorted_a[1], mz);
            //cout << "m ,z" << mz[0] << mz[1] << endl;
            a_plus_b_polygon[2] = { sorted_b[1].x , int64_t(mz[0] * (long double)sorted_b[1].x + mz[1] + sorted_b[1].y) };
            if (sorted_a[1].x >= sorted_b[2].x) {
                findequation(sorted_a[0], sorted_a[1], mz);
                //cout << "m ,z" << mz[0] << mz[1] << endl;
                a_plus_b_polygon[3] = { sorted_b[2].x , int64_t(mz[0] * (long double)sorted_b[2].x + mz[1]) };
                a_plus_b_polygon[4] = sorted_a[1];
                a_plus_b_polygon[5] = sorted_a[2];
            }
            else {
                findequation(sorted_b[1], sorted_b[2], mz);
                a_plus_b_polygon[3] = { sorted_a[1].x , int64_t(mz[0] * (long double)sorted_a[1].x + mz[1] + sorted_a[1].y) };
                if (sorted_b[2].x >= sorted_a[2].x) {
                    findequation(sorted_b[1], sorted_b[2], mz);
                    //cout << "m ,z" << mz[0] << mz[1] << endl;
                    a_plus_b_polygon[4] = { sorted_a[2].x , int64_t(mz[0] * (long double)sorted_a[2].x + mz[1]) };
                    a_plus_b_polygon[5] = sorted_b[2];
                }
                else {
                    findequation(sorted_a[1], sorted_a[2], mz);
                    //cout << "m ,z" << mz[0] << mz[1] << endl;
                    a_plus_b_polygon[4] = { sorted_b[2].x , int64_t(mz[0] * (long double)sorted_b[2].x + mz[1]) };
                    a_plus_b_polygon[5] = sorted_a[2];
                }
            }
        }
    }

}

point2D* ret(point2D* a) {
    point2D* b = a;
    return b;
}

// calculate area of triangle using the coordinates
long double trianglearea(point2D a[3])
{
    long double dArea = ((long double)(a[1].x - a[0].x) * (long double)(a[2].y - a[0].y) - (long double)(a[2].x - a[0].x) * (long double)(a[1].y - a[0].y)) / 2.0;
    return (dArea > 0.0) ? dArea : -dArea;
}

long double polygonArea(point2D* polygon, int n)
{
    vector<long double> X;
    vector<long double> Y;

    for (int i = 0; i < n; i++) {
        X.push_back(polygon[i].x);
        Y.push_back(polygon[i].y);
    }

    // Initialize area
    long double area = 0.0;

    // Calculate value of shoelace formula
    int j = n - 1;
    for (int i = 0; i < n; i++)
    {
        area += (X[j] + X[i]) * (Y[j] - Y[i]);
        j = i;  // j is previous vertex to i
    }
    // Return absolute value
    return abs(area / 2.0);
}

long double min(long double a, long double b) {
    return (a > b) ? b : a;
}

long double max(long double a, long double b) {
    return (a > b) ? a : b;
}

vector<point2D> remove_duplicate_vertices(point2D(&a_plus_b_polygon)[6]) {
    vector<point2D> polygon_a_plus_b_uni;
    polygon_a_plus_b_uni.reserve(6);
    bool flag;
    for (int i = 0; i < 6; i++) {
        flag = false;
        for (int j = 0; j < polygon_a_plus_b_uni.size(); j++) {
            if (polygon_a_plus_b_uni[j] == a_plus_b_polygon[i]) {
                flag = true;
            }
        }
        if (!flag)
            polygon_a_plus_b_uni.push_back(a_plus_b_polygon[i]);
    }
    return polygon_a_plus_b_uni;
}



long double getProbIdeal_1(long double param, long double tolerance = 0.1)
{
    if (param >= 1 - tolerance)
        return 1;
    else if (param <= tolerance)
        return 0;
    else
        return 1 / (1 - tolerance) * param;	//scale the 0-0.9 range to 0-1 for smoother numbers
}


long double getProbIdeal_0(long double param, long double tolerance = 0.1)
{
    return getProbIdeal_1(1 - param, tolerance);
}


long double getProbIdeal_05(long double param, long double tolerance = 0.1)
{
    long double factor = 0;
    if (abs(param - 0.5) <= tolerance)
        return 1;
    else if (param < 0.5)
        factor = 1 / (0.5 - tolerance) * param;
    else
    {
        factor = 1 / (tolerance - 0.5) * param + 1 - ((tolerance + 0.5) / (tolerance - 0.5));
    }

    return factor;
}



map<string, long double> getFateProbabilities(long double ig_a, long double ig_b, long double ia_g, long double ib_g, long double ig_a_plus_b, long double ia_plus_b_g, long double pa, long double pb, long double tolerance = 0.1)
{
    /*
    probabilities are based on the following table of ideal values for each fate
    cons = { "ig_a" : 1, "ig_b" : 1, "ia_g" : 1, "ib_g" : 1, "ia_plus_b_g" : 1/2, "pa" : 0, "pb" : 0 };
    neofunc of a = { "ig_a" : 0, "ig_b" : 1, "ia_g" : 0, "ib_g" : 1, "ia_plus_b_g" : 1/2, "pa" : 0, "pb" : 0 };
    neofunc of b = { "ig_a" : 1, "ig_b" : 0, "ia_g" : 1, "ib_g" : 0, "ia_plus_b_g" : 1/2, "pa" : 0, "pb" : 0 };
    specialization = { "ig_a" : 0, "ig_b" : 0, "ia_g" : 0, "ib_g" : 0, "ia_plus_b_g" : 0, "pa" : 0, "pb" : 0 };
    subfunc = { "ia_g" : 1, "ib_g" : 1, "ig_a" : 1/2, "ig_b" : 1/2, "ig_a_plus_b" : 1, "ia_plus_b_g" : 1, "pa" : 0, "pb" : 0 };
    pseudogene a = { "pa" : 1, "pb" : 0 };
    pseudogene b = { "pa" : 0, "pb" : 1 };
    */


    long double prob_conservation = getProbIdeal_1(ig_a, tolerance) *
        getProbIdeal_1(ig_b, tolerance) *
        getProbIdeal_1(ia_g, tolerance) *
        getProbIdeal_1(ib_g, tolerance) *
        getProbIdeal_05(ia_plus_b_g, tolerance) *
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);


    long double prob_neo_a = getProbIdeal_0(ig_a, tolerance) *
        getProbIdeal_1(ig_b, tolerance) *
        getProbIdeal_0(ia_g, tolerance) *
        getProbIdeal_1(ib_g, tolerance) *
        getProbIdeal_05(ia_plus_b_g, tolerance) *
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);

    long double prob_neo_b = getProbIdeal_1(ig_a, tolerance) *
        getProbIdeal_0(ig_b, tolerance) *
        getProbIdeal_1(ia_g, tolerance) *
        getProbIdeal_0(ib_g, tolerance) *
        getProbIdeal_05(ia_plus_b_g, tolerance) *
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);


    long double prob_spec = getProbIdeal_0(ig_a, tolerance) *
        getProbIdeal_0(ig_b, tolerance) *
        getProbIdeal_0(ia_g, tolerance) *
        getProbIdeal_0(ib_g, tolerance) *
        getProbIdeal_0(ia_plus_b_g, tolerance) *
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);

    long double prob_sub = getProbIdeal_05(ig_a, tolerance) *
        getProbIdeal_05(ig_b, tolerance) *
        getProbIdeal_1(ia_g, tolerance) *
        getProbIdeal_1(ib_g, tolerance) *
        getProbIdeal_1(ig_a_plus_b, tolerance) *
        getProbIdeal_1(ia_plus_b_g, tolerance) *
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);

    long double prob_pseudo_a = getProbIdeal_1(pa, tolerance) * getProbIdeal_0(pb, tolerance);

    long double prob_pseudo_b = getProbIdeal_0(pa, tolerance) * getProbIdeal_1(pb, tolerance);

    map<string, long double> probs = {
        {"cons", prob_conservation},
        {"neo_a", prob_neo_a},
        {"neo_b", prob_neo_b},
        {"spec", prob_spec},
        {"subfunc", prob_sub},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b},
    };

    return probs;

}




long double getLinPos(long double val, long double min, long double max)
{
    long double a = 1 / (max - min);
    long double b = -min / (max - min);

    long double fval = a * val + b;

    //min and max caused compile errors, so screw it
    if (fval < 0)
        fval = 0;
    if (fval > 1)
        fval = 1;

    return fval;
}


long double getLinNeg(long double val, long double min, long double max)
{
    long double a = -1 / (max - min);
    long double b = max / (max - min);

    long double fval = a * val + b;


    //min and max caused compile errors, so screw it
    if (fval < 0)
        fval = 0;
    if (fval > 1)
        fval = 1;


    return fval;
}





//long double g[3]: g[0](m), g[1](h), g[2](w)
map<string, long double> getFateProbabilities_v3(long double g[3], long double a[3], long double b[3], long double ig_a, long double ig_b, long double ia_g, long double ib_g, long double ig_a_plus_b, long double ia_plus_b_g, long double tolerance = 0.1)
{
    if (g[2] <= 0.0007 || (g[1] <= 0.0007 && g[1] >= -0.0007))
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"subfunc", 0},
        {"pseudo_a", 0},
        {"pseudo_b", 0}
        };
        return ptmp;
    }

    bool a_died = false;
    bool b_died = false;

    if (a[2] <= 0.0007 || (a[1] <= 0.0007 && a[1] >= -0.0007))
        a_died = true;
    if (b[2] <= 0.0007 || (b[1] <= 0.0007 && b[1] >= -00007))
        b_died = true;

    if (a_died || b_died)
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"subfunc", 0},
        {"pseudo_a",  (a_died) ? 1 : 0  },
        {"pseudo_b",  (b_died) ? 1 : 0  }
        };
        return ptmp;
    }



    long double pa = getLinPos(ia_g, 0.7, 1) * getLinNeg(ig_a, 0, 0.3);
    long double pb = getLinPos(ib_g, 0.7, 1) * getLinNeg(ig_b, 0, 0.3);


    long double tau = 0.05;


    long double prob_conservation =
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinNeg(ia_plus_b_g, 0.5, 1) *
        (1 - pa) * (1 - pb);



    long double prob_neo_a =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_neo_b =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);


    long double prob_spec =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_sub =
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ig_a_plus_b, 0 + tau, 1 - tau) *
        getLinPos(ia_plus_b_g, 0.5, 1) *
        (1 - pa) * (1 - pb);

    //this is a hack to ensure that both a and b are needed to perform g.  I don't know how to integrate that elegantly in the probabilities product
    if (ig_a > 0.8 || ig_b > 0.8)
        prob_sub = 0;


    long double prob_pseudo_a = pa * (1 - pb);

    long double prob_pseudo_b = pb * (1 - pa);



    map<string, long double> probs = {
        {"cons", prob_conservation},
        {"neo_a", prob_neo_a},
        {"neo_b", prob_neo_b},
        {"spec", prob_spec},
        {"subfunc", prob_sub},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b},
    };



    return probs;
}
















//long double g[3]: g[0](m), g[1](h), g[2](w)
map<string, long double> getFateProbabilities_v4(long double g[3], long double a[3], long double b[3], long double ig_a, long double ig_b, long double ia_g, long double ib_g, long double ig_a_plus_b, long double ia_plus_b_g, long double tolerance = 0.1)
{
    if (g[2] <= 0.0005 || (g[1] <= 0.0005 && g[1] >= -0.0005))
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"dblneo", 0},
        {"subfunc", 0},
        {"pseudo_a", 0},
        {"pseudo_b", 0}
        };
        return ptmp;
    }

    long double pa, pb;
    if (a[2] <= 0.0005 || (a[1] <= 0.0005 && a[1] >= -0.0005))
        pa = 1;
    else
        pa = getLinPos(ia_g, 0.7, 1) * getLinNeg(ig_a, 0, 0.3);


    if (b[2] <= 0.0005 || (b[1] <= 0.0005 && b[1] >= -0.0005))
        pb = 1;
    else
        pb = getLinPos(ib_g, 0.7, 1) * getLinNeg(ig_b, 0, 0.3);


    long double prob_pseudo_a = pa * (1 - pb);
    long double prob_pseudo_b = pb * (1 - pa);




    if (pa == 1 || pb == 1)
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"dblneo", 0},
        {"subfunc", 0},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b}
        };

        return ptmp;

    }


    long double tau = 0.02;


    long double prob_conservation =
        getLinPos(ig_a, 0.5, 1 - tau) *
        getLinPos(ig_b, 0.5, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);



    long double prob_neo_a =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_neo_b =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);


    long double prob_spec =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_dblneo =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_sub =
        getLinPos(ig_a, 0 + tau, 0.25) *
        getLinNeg(ig_a, 0.75, 1 - tau) *
        getLinPos(ig_b, 0 + tau, 0.25) *
        getLinNeg(ig_b, 0.75, 1 - tau) *
        getLinPos(ig_a_plus_b, 0 + tau, 1 - tau) *
        getLinPos(ia_plus_b_g, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);






    long double prob_consa_lossb =
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(pb, 0 + tau, 1 - tau);

    long double prob_consb_lossa =
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(pa, 0 + tau, 1 - tau);

    long double prob_neoa_lossb =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinPos(pb, 0 + tau, 1 - tau);

    long double prob_neob_lossa =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinPos(pa, 0 + tau, 1 - tau);






    map<string, long double> probs = {
        {"cons", prob_conservation},
        {"neo_a", prob_neo_a},
        {"neo_b", prob_neo_b},
        {"spec", prob_spec},
        {"dblneo", prob_dblneo},
        {"subfunc", prob_sub},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b},
        {"consa_lossb", prob_consa_lossb},
        {"consb_lossa", prob_consb_lossa},
        {"neoa_lossb", prob_neoa_lossb},
        {"neob_lossa", prob_neob_lossa},
    };



    return probs;
}



























map<string, long double> getFateProbabilities_v5(long double g[3], long double a[3], long double b[3], long double ig_a, long double ig_b, long double ia_g, long double ib_g, long double ig_a_plus_b, long double ia_plus_b_g, long double tolerance = 0.1)
{
    if (g[2] <= 0.0005 || (g[1] <= 0.0005 && g[1] >= -0.0005))
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"dblneo", 0},
        {"subfunc", 0},
        {"pseudo_a", 0},
        {"pseudo_b", 0}
        };
        return ptmp;
    }

    long double pa, pb;
	long double pseudo_threshold = 0.2;

    if (a[2] <= 0.0005 || (a[1] <= 0.0005 && a[1] >= -0.0005))
	{
        pa = 1;
	}
    else
	{
        pa = getLinPos(ia_g, 1 - pseudo_threshold, 1) * getLinNeg(ig_a, 0, pseudo_threshold);
	}


    if (b[2] <= 0.0005 || (b[1] <= 0.0005 && b[1] >= -0.0005))
        pb = 1;
    else
        pb = getLinPos(ib_g, 1 - pseudo_threshold, 1) * getLinNeg(ig_b, 0, pseudo_threshold);


    long double prob_pseudo_a = pa;	// * (1 - pb);
    long double prob_pseudo_b = pb; //* (1 - pa);




    if (pa == 1 || pb == 1)
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"dblneo", 0},
        {"subfunc", 0},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b}
        };

        return ptmp;

    }


    long double tau = 0.02;


    



    long double prob_neo_a =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_neo_b =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);


    long double prob_spec =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);


    long double prob_dblneo =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

	
	long double a_plus_b_subfunc_cost = -3.0 * ia_plus_b_g * ia_plus_b_g + 13.0/2.0 * ia_plus_b_g - 5.0/2.0;
	if (a_plus_b_subfunc_cost < 0)
		a_plus_b_subfunc_cost = 0.0;
	if (a_plus_b_subfunc_cost > 1)
		a_plus_b_subfunc_cost = 1.0;
	
	long double prob_conservation =
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
		getLinPos(ig_a_plus_b, 0 + tau, 1 - tau) *
        //getLinNeg(ia_plus_b_g, 0.5, 0.833333) *
		(1 - a_plus_b_subfunc_cost) * 
		
        (1 - pa) * (1 - pb);


    long double prob_sub =
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ig_a_plus_b, 0 + tau, 1 - tau) *
		a_plus_b_subfunc_cost * 
        //getLinPos(ia_plus_b_g, 0.5, 0.83333) *				//previously, was 1 instead of 0.83333
        (1 - pa) * (1 - pb);






    long double prob_consa_lossb =
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(pb, 0 + tau, 1 - tau);

    long double prob_consb_lossa =
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(pa, 0 + tau, 1 - tau);

    long double prob_neoa_lossb =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinPos(pb, 0 + tau, 1 - tau);

    long double prob_neob_lossa =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinPos(pa, 0 + tau, 1 - tau);






    map<string, long double> probs = {
        {"cons", prob_conservation},
        {"neo_a", prob_neo_a},
        {"neo_b", prob_neo_b},
        {"spec", prob_spec},
        {"dblneo", prob_dblneo},
        {"subfunc", prob_sub},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b},
        {"consa_lossb", prob_consa_lossb},
        {"consb_lossa", prob_consb_lossa},
        {"neoa_lossb", prob_neoa_lossb},
        {"neob_lossa", prob_neob_lossa},
    };



    return probs;
}









map<string, long double> getFateProbabilities_v45(long double g[3], long double a[3], long double b[3], long double ig_a, long double ig_b, long double ia_g, long double ib_g, long double ig_a_plus_b, long double ia_plus_b_g, long double tolerance = 0.1)
{
    if (g[2] <= 0.0005 || (g[1] <= 0.0005 && g[1] >= -0.0005))
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"dblneo", 0},
        {"subfunc", 0},
        {"pseudo_a", 0},
        {"pseudo_b", 0}
        };
        return ptmp;
    }

    long double pa, pb;
	long double pseudo_threshold = 0.2;

    if (a[2] <= 0.0005 || (a[1] <= 0.0005 && a[1] >= -0.0005))
	{
        pa = 1;
	}
    else
	{
        pa = getLinPos(ia_g, 1 - pseudo_threshold, 1) * getLinNeg(ig_a, 0, pseudo_threshold);
	}


    if (b[2] <= 0.0005 || (b[1] <= 0.0005 && b[1] >= -0.0005))
        pb = 1;
    else
        pb = getLinPos(ib_g, 1 - pseudo_threshold, 1) * getLinNeg(ig_b, 0, pseudo_threshold);


    long double prob_pseudo_a = pa;	// * (1 - pb);
    long double prob_pseudo_b = pb; //* (1 - pa);




    if (pa == 1 || pb == 1)
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"dblneo", 0},
        {"subfunc", 0},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b}
        };

        return ptmp;

    }


    long double tau = 0.000001;


    
    long double prob_neo_a =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0.25, 1 - tau) *
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_neo_b =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0.25, 1 - tau) *
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);


    /*long double prob_neo_a =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinNeg(ig_a, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_neo_b =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinNeg(ig_b, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);*/


    long double prob_spec =
	getLinPos(ig_a_plus_b, 0 + tau, 1 - tau) *
        getLinNeg(ia_g, 0.25, 1) *
        getLinNeg(ib_g, 0.25, 1) *
        (1 - pa) * (1 - pb);


	
    long double prob_dblneo =
	getLinNeg(ig_a, 0 + tau, 0.5) *
        getLinNeg(ig_b, 0 + tau, 0.5) *
	getLinNeg(ia_g, 0, 0.5) *
        getLinNeg(ib_g, 0, 0.5) *
        (1 - pa) * (1 - pb);
    /*long double prob_dblneo = getLinNeg(ig_a_plus_b, 0, 0.5) *
        (1 - pa) * (1 - pb);*/

	
	long double prob_conservation =
        getLinPos(ia_g, 0.25, 1 - tau) *
        getLinPos(ib_g, 0.25, 1 - tau) *
	getLinPos(ig_a_plus_b, 0 + tau, 1 - tau) *
        getLinNeg(ia_plus_b_g, 0.5, 1) *		
        (1 - pa) * (1 - pb);


    long double prob_sub =
        getLinPos(ia_g, 0.25, 1 - tau) *
        getLinPos(ib_g, 0.25, 1 - tau) *
        getLinPos(ig_a_plus_b, 0 + tau, 1 - tau) *
	getLinPos(ia_plus_b_g, 0.5, 1) *				
        (1 - pa) * (1 - pb);






    long double prob_consa_lossb =
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(pb, 0 + tau, 1 - tau);

    long double prob_consb_lossa =
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(pa, 0 + tau, 1 - tau);

    long double prob_neoa_lossb =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinPos(pb, 0 + tau, 1 - tau);

    long double prob_neob_lossa =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinPos(pa, 0 + tau, 1 - tau);






    map<string, long double> probs = {
        {"cons", prob_conservation},
        {"neo_a", prob_neo_a},
        {"neo_b", prob_neo_b},
        {"spec", prob_spec},
        {"dblneo", prob_dblneo},
        {"subfunc", prob_sub},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b},
        {"consa_lossb", prob_consa_lossb},
        {"consb_lossa", prob_consb_lossa},
        {"neoa_lossb", prob_neoa_lossb},
        {"neob_lossa", prob_neob_lossa},
    };



    return probs;
}














/*long double getLineCoverage(long double l1, long double r1, long double l2, long double, r2)
{
	if (l1 < l2 && l2 < r1)
	{
		return (r1 - l2)/(r1 - l1);
	}
	if (l2 < l1 && l1 < r2)
	{
		return (r2 - l1)/(r1 - l1);
	}
	return 0;
}


map<string, long double> getFateProbabilities_v6(long double g[3], long double a[3], long double b[3], long double ig_a, long double ig_b, long double ia_g, long double ib_g, long double ig_a_plus_b, long double ia_plus_b_g, long double tolerance = 0.1)
{
    if (g[2] <= 0.0005 || (g[1] <= 0.0005 && g[1] >= -0.0005))
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"subfunc", 0},
        {"pseudo_a", 0},
        {"pseudo_b", 0}
        };
        return ptmp;
    }
	
	
	long double lg = g[0] - g[2];
	long double rg = g[0] + g[2];
	long double la = a[0] - a[2];
	long double ra = a[0] + a[2];
	long double lb = b[0] - b[2];
	long double rb = b[0] + b[2];
	
	
	ig_a = getLineCoverage(lg, rg, la, ra);
	ig_b = getLineCoverage(lg, rg, lb, rb);
	
	ia_g = getLineCoverage(la, ra, lg, rg);
	ib_g = getLineCoverage(lb, rb, lg, rg);
	
	
	
	
	
	
	
	
	

    long double pa, pb;
    if (a[2] <= 0.0005 || (a[1] <= 0.0005 && a[1] >= -0.0005))
        pa = 1;
    else
        pa = getLinPos(ia_g, 0.7, 1) * getLinNeg(ig_a, 0, 0.3);


    if (b[2] <= 0.0005 || (b[1] <= 0.0005 && b[1] >= -0.0005))
        pb = 1;
    else
        pb = getLinPos(ib_g, 0.7, 1) * getLinNeg(ig_b, 0, 0.3);


    long double prob_pseudo_a = pa * (1 - pb);
    long double prob_pseudo_b = pb * (1 - pa);




    if (pa == 1 || pb == 1)
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"subfunc", 0},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b}
        };

        return ptmp;

    }


    long double tau = 0.02;


    long double prob_conservation =
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinNeg(ia_plus_b_g, 0.5, 1) *
	getLinPos(ig_a_plus_b, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);



    long double prob_neo_a =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_neo_b =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);


    long double prob_spec =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        (1 - pa) * (1 - pb);

    long double prob_sub =
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(ig_a_plus_b, 0 + tau, 1 - tau) *
        getLinPos(ia_plus_b_g, 0.5, 1 - tau) *
        (1 - pa) * (1 - pb);






    long double prob_consa_lossb =
        getLinPos(ig_a, 0 + tau, 1 - tau) *
        getLinPos(ia_g, 0 + tau, 1 - tau) *
        getLinPos(pb, 0 + tau, 1 - tau);

    long double prob_consb_lossa =
        getLinPos(ig_b, 0 + tau, 1 - tau) *
        getLinPos(ib_g, 0 + tau, 1 - tau) *
        getLinPos(pa, 0 + tau, 1 - tau);

    long double prob_neoa_lossb =
        getLinNeg(ia_g, 0 + tau, 1 - tau) *
        getLinPos(pb, 0 + tau, 1 - tau);

    long double prob_neob_lossa =
        getLinNeg(ib_g, 0 + tau, 1 - tau) *
        getLinPos(pa, 0 + tau, 1 - tau);






    map<string, long double> probs = {
        {"cons", prob_conservation},
        {"neo_a", prob_neo_a},
        {"neo_b", prob_neo_b},
        {"spec", prob_spec},
        {"subfunc", prob_sub},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b},
        {"consa_lossb", prob_consa_lossb},
        {"consb_lossa", prob_consb_lossa},
        {"neoa_lossb", prob_neoa_lossb},
        {"neob_lossa", prob_neob_lossa},
    };



    return probs;
}
*/


































map<string, long double> getFateProbabilities_v2(long double g[3], long double a[3], long double b[3], long double ig_a, long double ig_b, long double ia_g, long double ib_g, long double ig_a_plus_b, long double ia_plus_b_g, long double tolerance = 0.1)
{
    if (g[1] <= 0.0007 || g[2] <= 0.0007)
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"subfunc", 0},
        {"pseudo_a", 0},
        {"pseudo_b", 0}
        };
        return ptmp;
    }
    if (a[1] <= 0.0007 || a[2] <= 0.0007 || b[1] <= 0.0007 || b[2] <= 0.0007)
    {
        map<string, long double> ptmp = {
        {"cons", 0},
        {"neo_a", 0},
        {"neo_b", 0},
        {"spec", 0},
        {"subfunc", 0},
        {"pseudo_a",  (a[1] == 0 || a[2] == 0) ? 1 : 0  },
        {"pseudo_b",  (b[1] == 0 || b[2] == 0) ? 1 : 0  }
        };
        return ptmp;
    }


    long double pa = 0;
    long double pb = 0;



    long double prob_conservation =
        getProbIdeal_1(ig_a, tolerance) *
        getProbIdeal_1(ig_b, tolerance) *
        getProbIdeal_1(ia_g, tolerance) *
        getProbIdeal_1(ib_g, tolerance) *
        getProbIdeal_05(ia_plus_b_g, tolerance) *
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);


    long double prob_neo_a =
        getProbIdeal_0(ig_a, tolerance) *
        getProbIdeal_1(ig_b, tolerance) *
        getProbIdeal_0(ia_g, tolerance) *
        getProbIdeal_1(ib_g, tolerance) *
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);

    long double prob_neo_b =
        getProbIdeal_1(ig_a, tolerance) *
        getProbIdeal_0(ig_b, tolerance) *
        getProbIdeal_1(ia_g, tolerance) *
        getProbIdeal_0(ib_g, tolerance) *
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);


    long double prob_spec =
        getProbIdeal_0(ig_a, tolerance) *
        getProbIdeal_0(ig_b, tolerance) *
        getProbIdeal_0(ia_g, tolerance) *
        getProbIdeal_0(ib_g, tolerance) *
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);

    long double prob_sub =
        getProbIdeal_1(ia_g, tolerance) *
        getProbIdeal_1(ib_g, tolerance) *
        getProbIdeal_1(ig_a_plus_b, tolerance) *
        //getProbIdeal_1(ia_plus_b_g, tolerance) *
        max(0, getProbIdeal_0(2 * (1 - ia_plus_b_g), tolerance)) * 	//this is a trick to make ia_plus_b_g be between 0.5 and 1
        getProbIdeal_0(pa, tolerance) *
        getProbIdeal_0(pb, tolerance);

    //this is a hack to ensure that both a and b are needed to perform g.  I don't know how to integrate that elegantly in the probabilities product
    if (ig_a > 0.8 || ig_b > 0.8)
        prob_sub = 0;


    long double prob_pseudo_a = getProbIdeal_1(pa, tolerance) * getProbIdeal_0(pb, tolerance);

    long double prob_pseudo_b = getProbIdeal_0(pa, tolerance) * getProbIdeal_1(pb, tolerance);



    map<string, long double> probs = {
        {"cons", prob_conservation},
        {"neo_a", prob_neo_a},
        {"neo_b", prob_neo_b},
        {"spec", prob_spec},
        {"subfunc", prob_sub},
        {"pseudo_a", prob_pseudo_a},
        {"pseudo_b", prob_pseudo_b},
    };

    return probs;
}



























//long double g[3]: g[0](m), g[1](h), g[2](w)
vector<long double> ProbabilitiesCalculation(long double g[3], long double a[3], long double b[3], int prob_version = 5) {
    // convert m, h, w to three vertices

    point2D g2d[3] = {};
    point2D a2d[3] = {};
    point2D b2d[3] = {};
    
    
    
    long double multiplier = 1000.0;
    long double h_multiplier = multiplier;
    
    //TODO: the following is a hack - I don't know why it's needed, but if all triangles are negative, we get different values!
    if (g[1] < 0 && a[1] < 0 && b[1] < 0)
    {
    	h_multiplier *= -1;
    }
    
    g2d[0] = point2D{ (int64_t(g[0] * multiplier) - int64_t(g[2] * multiplier)) * 100, 0 };
    g2d[1] = point2D{ (int64_t(g[0] * multiplier) + int64_t(g[2] * multiplier)) * 100, 0 };
    g2d[2] = point2D{ int64_t(g[0]  * multiplier) * 100 , int64_t(g[1]  * h_multiplier) * 100 };

    a2d[0] = point2D{ (int64_t(a[0] * multiplier) - int64_t(a[2] * multiplier)) * 100, 0 };
    a2d[1] = point2D{ (int64_t(a[0] * multiplier) + int64_t(a[2] * multiplier)) * 100, 0 };
    a2d[2] = point2D{ int64_t(a[0] * multiplier) * 100, int64_t(a[1] * h_multiplier) * 100 };

    b2d[0] = point2D{ (int64_t(b[0] * multiplier) - int64_t(b[2] * multiplier)) * 100, 0 };
    b2d[1] = point2D{ (int64_t(b[0] * multiplier) + int64_t(b[2] * multiplier)) * 100, 0 };
    b2d[2] = point2D{ int64_t(b[0]  * multiplier) * 100, int64_t(b[1]  * h_multiplier) * 100};




    //sorttrianglevertices(g2d);
    //sorttrianglevertices(a2d);
    //sorttrianglevertices(b2d);

    //area of g, a ,b
    long double area_g = polygonArea(g2d, 3);
    long double area_a = polygonArea(a2d, 3);
    long double area_b = polygonArea(b2d, 3);

    //cout << "area g: " << area_g << endl;
    //cout << "area a: " << area_a << endl;
    //cout << "area b: " << area_b << endl;

    //point2D subjectPolygon[] = { g1, g2, g3 };
    int gPolygonSize = sizeof(g2d) / sizeof(g2d[0]);
    int aPolygonSize = sizeof(a2d) / sizeof(a2d[0]);
    int bPolygonSize = sizeof(b2d) / sizeof(b2d[0]);

    // define the new clipped polygon (empty)
    //int newPolygonSize = 0;
    //point2D newPolygon[N] = { };
    //point2D* newPolygon = new point2D[];


    int newPolygonSize_ga = 0;
    point2D newPolygon_ga[N] = { };

    // apply clipping
    SutherlandHodgman(g2d, gPolygonSize, a2d, aPolygonSize, newPolygon_ga, newPolygonSize_ga);
    long double inersection_area_g_a;
    long double inersection_area_g_b;
    long double area_a_plus_b;
    long double intersection_area_g_a_plus_b;
    // calculate intersection area
    /*point2D a2 = newPolygon[2];
    newPolygon[2] = newPolygon[3];
    newPolygon[3] = a2;*/
    //for (int i = 0; i < 3; i++)
        //cout << "(" << newPolygon[i].x << ", " << newPolygon[i].y << ")" << endl;
    if (newPolygonSize_ga != 0)
        inersection_area_g_a = polygonArea(newPolygon_ga, newPolygonSize_ga);
    else
        inersection_area_g_a = 0;

    int newPolygonSize_gb = 0;
    point2D newPolygon_gb[N] = { };
    

    SutherlandHodgman(g2d, gPolygonSize, b2d, bPolygonSize, newPolygon_gb, newPolygonSize_gb);
    if (newPolygonSize_gb != 0)
        inersection_area_g_b = polygonArea(newPolygon_gb, newPolygonSize_gb);
    else
        inersection_area_g_b = 0;
    
    
    

    int newPolygonSize_ab = 0;
    point2D newPolygon_ab[N] = { };

    SutherlandHodgman(a2d, aPolygonSize, b2d, bPolygonSize, newPolygon_ab, newPolygonSize_ab);
    if (newPolygonSize_ab == 0) {
        //cout << "a and b are separte. " << endl;
        area_a_plus_b = area_a + area_b;
        intersection_area_g_a_plus_b = inersection_area_g_b + inersection_area_g_a;
    }
    else {
        //cout << "a and b have intersection. " << endl;
        point2D polygon_a_plus_b[6]{};
        vector<point2D> polygon_a_plus_b_uni;
        calculate_a_plus_b_polygon(b2d, a2d, polygon_a_plus_b);
        polygon_a_plus_b_uni = remove_duplicate_vertices(polygon_a_plus_b);
        /*for (int i = 0; i < 6; i++) {
            polygon_a_plus_b_uni.push_back(polygon_a_plus_b[i]);
        }*/
        int polygon_a_plus_bSize = polygon_a_plus_b_uni.size();
        for (int k = 0; k < polygon_a_plus_b_uni.size(); k++)
            polygon_a_plus_b[k] = polygon_a_plus_b_uni[k];

        area_a_plus_b = polygonArea(polygon_a_plus_b, polygon_a_plus_bSize);


	//ML's debugging
	/*cout << std::fixed;
	cout << std::setprecision(15);
	cout<<"a+b points = "<<polygon_a_plus_bSize<<endl;
	for (int i = 0; i < polygon_a_plus_bSize; ++i)
	{
		cout<<polygon_a_plus_b[i].x<<" "<<polygon_a_plus_b[i].y<<endl;
	}	
	cout<<"g points"<<endl;
	for (int i = 0; i < 3; ++i)
	{
		cout<<g2d[i].x<<" "<<g2d[i].y<<"   g[i] = "<<g[i]<<"   a[i]="<<a[i]<<"   b[i]="<<b[i]<<endl;
	}*/
		
        /*point2D polygon_a_plus_b_sorted[6]{};
        polygon_a_plus_b_sorted[0] = polygon_a_plus_b[0];
        int k = 1;
        for (int i = polygon_a_plus_bSize-1; i > 0; i--){
          polygon_a_plus_b_sorted[k] = polygon_a_plus_b[i];
          k++;
        }*/
        //sorttrianglevertices(g2d);
        int newPolygonSize_aplusb_g = 0;
        point2D newPolygon_aplusb_g[N] = { };

        SutherlandHodgman(polygon_a_plus_b, polygon_a_plus_bSize, g2d, gPolygonSize, newPolygon_aplusb_g, newPolygonSize_aplusb_g);
        intersection_area_g_a_plus_b = polygonArea(newPolygon_aplusb_g, newPolygonSize_aplusb_g);

            
    }
    
    if (isnan(area_a_plus_b))
    {
    	cout<<"NAN"<<endl;
    	int xx;
    	cin>>xx;
    
    }

    //cout << "area of intersection of g, a + b: " << intersection_area_g_a_plus_b << endl;

    long double threshold_pseudogene = 0.7;
    long double ig_a = 0;
    if (area_g > 0)
    	ig_a = (inersection_area_g_a / area_g);

    long double ig_b;
    if (area_g > 0)
        ig_b = (inersection_area_g_b / area_g);    
    
    
    if (b[2] <= g[2] / 1.5 && ig_b > 0.9)
    //if (ig_b < 0.95)
    {
    cout<<"*** PROLEM"<<endl;
    cout<<"polygon size="<<newPolygonSize_gb<<endl;
    for (int i = 0; i < newPolygonSize_gb; ++i)
    {
    	cout<<newPolygon_gb[i].x<<" "<<newPolygon_gb[i].y<<endl;
    
    }
    cout<<g[0]<<" "<<g[1]<<" "<<g[2]<<endl;
        cout<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
            cout<<b[0]<<" "<<b[1]<<" "<<b[2]<<endl;
    
    cout<<"G coord"<<endl;
    cout<<g2d[0].x<<" "<<g2d[0].y<<endl
        <<g2d[1].x<<" "<<g2d[1].y<<endl
        <<g2d[2].x<<" "<<g2d[2].y<<endl;
            
    cout << "area of intersection of a , g: " << inersection_area_g_a << endl;
    cout << "area of intersection of b , g: " << inersection_area_g_b << endl;
    cout<<"area g: "<<area_g<<endl;
    cout<<"ig_b:"<<ig_b<<endl;

    }
    

    long double ia_g = 0;
    if (area_a > 0)
        ia_g = (inersection_area_g_a / area_a);
        
    long double ib_g = 0;
    if (area_b > 0)
	ib_g = (inersection_area_g_b / area_b);

    //cout << "ib_g: " << ib_g << endl;
    long double ig_a_plus_b = 0;
    if (area_g > 0)
    	ig_a_plus_b = (intersection_area_g_a_plus_b / area_g);

    long double ia_plus_b_g = 0;
    if (area_a_plus_b > 0)
    	ia_plus_b_g = (intersection_area_g_a_plus_b / area_a_plus_b);

    long double pa = 1 - min(abs(a[1] / (g[1] * threshold_pseudogene)), 1);
    long double pb = 1 - min(abs(b[1] / (g[1] * threshold_pseudogene)), 1);





    vector<long double> fates_probablities;
    fates_probablities.resize(20);

    //previous way
    /*long double P_conservation = ig_a * ig_b;
    long double P_newfunctionlizatoin = max(ig_b * (1 - ia_g) * (1 - pa), ig_a * (1 - ib_g) * (1 - pb));
    long double P_specialization = (1 - pa) * (1 - pb) * (1 - ig_a_plus_b);
    long double P_pseudogenization = max(pa, pb);
    long double P_subfunctionlization = ig_a * ib_g * ig_a_plus_b * ia_plus_b_g * (1 - pa) * (1 - pb);

    fates_probablities[0] = P_subfunctionlization;
    fates_probablities[1] = P_conservation;
    fates_probablities[2] = P_newfunctionlizatoin;
    fates_probablities[3] = P_pseudogenization;
    fates_probablities[4] = P_specialization;*/


    //signature is
    //getFateProbabilities(long double ig_a, long double ig_b, long double ia_g, long double ib_g, long double ig_a_plus_b, long double ia_plus_b_g, long double pa, long double pb, long double tolerance = 0.1)
    //map<string, long double> probs = getFateProbabilities(ig_a, ig_b, ia_g, ib_g, ig_a_plus_b, ia_plus_b_g, pa, pb, 0.1);


    //getFateProbabilities_v2(long double g[3], long double a[3], long double b[3], long double ig_a, long double ig_b, long double ia_g, long double ib_g, long double ig_a_plus_b, long double ia_plus_b_g, long double tolerance = 0.1)
    map<string, long double> probs;
    
    if (prob_version == 4)
    	probs = getFateProbabilities_v4(g, a, b, ig_a, ig_b, ia_g, ib_g, ig_a_plus_b, ia_plus_b_g, 0.1);
    else if (prob_version == 5)
        probs = getFateProbabilities_v5(g, a, b, ig_a, ig_b, ia_g, ib_g, ig_a_plus_b, ia_plus_b_g, 0.1);
	else if (prob_version == 45)
		probs = getFateProbabilities_v45(g, a, b, ig_a, ig_b, ia_g, ib_g, ig_a_plus_b, ia_plus_b_g, 0.1);
    else
        cout<<"Error: version does not exist."<<endl;

    fates_probablities[0] = probs["subfunc"];
    fates_probablities[1] = probs["cons"];
    fates_probablities[2] = max(probs["neo_a"], probs["neo_b"]);
    fates_probablities[3] = max(probs["pseudo_a"], probs["pseudo_b"]);
    fates_probablities[4] = probs["spec"];
    fates_probablities[13] = probs["dblneo"];

    fates_probablities[5] = ig_a;
    fates_probablities[6] = ig_b;
    fates_probablities[7] = ia_g;
    fates_probablities[8] = ib_g;

    fates_probablities[9] = ig_a_plus_b;
    fates_probablities[10] = ia_plus_b_g;
    fates_probablities[11] = probs["pseudo_a"];
    fates_probablities[12] = probs["pseudo_b"];


    fates_probablities[14] = probs["consa_lossb"];
    fates_probablities[15] = probs["consb_lossa"];
    fates_probablities[16] = probs["neoa_lossb"];
    fates_probablities[17] = probs["neob_lossa"];

    fates_probablities[18] = max(probs["consa_lossb"], probs["consb_lossa"]);
    fates_probablities[19] = max(probs["neoa_lossb"], probs["neob_lossa"]);


    /*if (counter == 0) {
        static std::ofstream trianglevariables("triangle_fates_variables.csv", std::ofstream::trunc);
        trianglevariables << "g.m, g.h, g.w, a.m, a.h, a.w, b.m, b.h, b.w, P_subfunctionlization, P_conservation, P_newfunctionlizatoin, P_pseudogenization, P_specialization, ig_a, ig_b, ia_g, ib_g, ig_a_plus_b, ia_plus_b_g, pa, pb " << std::endl;
        trianglevariables << g[0] << "," << g[1] << "," << g[2] << ","
            << a[0] << "," << a[1] << "," << a[2] << ","
            << b[0] << "," << b[1] << "," << b[2] << ","
            << fates_probablities[0] << "," << fates_probablities[1] << "," << fates_probablities[2] << "," << fates_probablities[3] << "," << fates_probablities[4] << ","
            << fates_probablities[5] << "," << fates_probablities[6] << "," << fates_probablities[7] << "," << fates_probablities[8] << "," << fates_probablities[9] << ","
            << fates_probablities[10] << "," << fates_probablities[11] << "," << fates_probablities[12] << std::endl;
        counter++;
    }*/
    /* cout << "P_subfunctionlization: " << P_subfunctionlization << " P_conservation: " << P_conservation << " P_newfunctionlizatoin " << P_newfunctionlizatoin << " P_pseudogenization " << P_pseudogenization << " P_specialization " << P_specialization << endl;

     if (P_conservation > P_newfunctionlizatoin && P_conservation > P_specialization && P_conservation > P_pseudogenization && P_conservation > P_subfunctionlization)
     {
         printf("Conservatoin");
     }
     else if (P_newfunctionlizatoin > P_specialization && P_newfunctionlizatoin > P_pseudogenization && P_newfunctionlizatoin > P_subfunctionlization)
     {
         printf("Newfunctionlizatoin");
     }
     else if (P_specialization > P_pseudogenization && P_specialization > P_subfunctionlization)

     {
         printf("Specialization");
     }
     else if (P_pseudogenization > P_subfunctionlization)
     {
         printf("Pseudogenization");
     }
     else
     {

         printf("Subfunctionlization");
     }*/
    return fates_probablities;
}
