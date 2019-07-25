__device__ void bicubic(real4 y, real4 y1i, real4 y2i, real4 y12i, real x1, real x1l, real x1u,
                        real x2, real x2l, real x2u, real* energyOut, real* dang1Out, real* dang2Out) {
    real c[4][4];
    real d1 = x1u - x1l;
    real d2 = x2u - x2l;
    real d12 = d1*d2;
    real4 y1 = d1*y1i;
    real4 y2 = d2*y2i;
    real4 y12 = d12*y12i;

    c[0][0] = y.x;
    c[0][1] = y2.x;
    c[0][2] = 3.0f*(y.w - y.x) - (2.0f*y2.x + y2.w);
    c[0][3] = 2.0f*(y.x - y.w) + y2.x + y2.w;
    c[1][0] = y1.x;
    c[1][1] = y12.x;
    c[1][2] = 3.0f*(y1.w - y1.x) - (2.0f*y12.x + y12.w);
    c[1][3] = 2.0f*(y1.x - y1.w) + y12.x + y12.w;
    c[2][0] = 3.0f*(y.y - y.x) - (2.0f*y1.x + y1.y);
    c[2][1] = 3.0f*(y2.y - y2.x) - (2.0f*y12.x + y12.y);
    c[2][2] = 9.0f*(y.x - y.y + y.z - y.w) + 6.0f* y1.x  +  3.0f* y1.y  -  3.0f* y1.z  -  6.0f* y1.w +
                                             6.0f* y2.x  -  6.0f* y2.y  -  3.0f* y2.z  +  3.0f* y2.w +
                                             4.0f*y12.x  +  2.0f*y12.y  +       y12.z  +  2.0f*y12.w;
    c[2][3] = 6.0f*(y.y - y.x + y.w - y.z) + -4.0f* y1.x  -  2.0f* y1.y  +  2.0f* y1.z  +  4.0f* y1.w +
                                             -3.0f* y2.x  +  3.0f* y2.y  +  3.0f* y2.z  -  3.0f* y2.w
                                             -2.0f*y12.x  -       y12.y  -       y12.z  -  2.0f*y12.w;
    c[3][0] = 2.0f*(y.x - y.y) + y1.x + y1.y;
    c[3][1] = 2.0f*(y2.x - y2.y) + y12.x + y12.y;
    c[3][2] = 6.0f*(y.y -  y.x +  y.w -  y.z) +
              3.0f*(y1.z + y1.w - y1.x - y1.y) +
              2.0f*(2.0f*(y2.y - y2.x) + y2.z - y2.w) +
             -2.0f*(y12.x + y12.y) - y12.z - y12.w;
    c[3][3] = 4.0f*( y.x -    y.y  +   y.z -  y.w)  +
              2.0f*(y1.x +   y1.y  -  y1.z -  y1.w) +
              2.0f*(y2.x -   y2.y  -  y2.z +  y2.w) +
                    y12.x +  y12.y  + y12.z + y12.w;

    real t = (x1-x1l) / (x1u-x1l);
    real u = (x2-x2l) / (x2u-x2l);

    real energy =       ((c[3][3]*u + c[3][2])*u + c[3][1])*u + c[3][0];
    energy = t*energy + ((c[2][3]*u + c[2][2])*u + c[2][1])*u + c[2][0];
    energy = t*energy + ((c[1][3]*u + c[1][2])*u + c[1][1])*u + c[1][0];
    energy = t*energy + ((c[0][3]*u + c[0][2])*u + c[0][1])*u + c[0][0];

    real dang1 =           (3.0f*c[3][3]*t + 2.0f*c[2][3])*t + c[1][3];
         dang1 = u*dang1 + (3.0f*c[3][2]*t + 2.0f*c[2][2])*t + c[1][2];
         dang1 = u*dang1 + (3.0f*c[3][1]*t + 2.0f*c[2][1])*t + c[1][1];
         dang1 = u*dang1 + (3.0f*c[3][0]*t + 2.0f*c[2][0])*t + c[1][0];

    real dang2 =           (3.0f*c[3][3]*u + 2.0f*c[3][2])*u + c[3][1];
         dang2 = t*dang2 + (3.0f*c[2][3]*u + 2.0f*c[2][2])*u + c[2][1];
         dang2 = t*dang2 + (3.0f*c[1][3]*u + 2.0f*c[1][2])*u + c[1][1];
         dang2 = t*dang2 + (3.0f*c[0][3]*u + 2.0f*c[0][2])*u + c[0][1];

    dang1 = dang1 / (x1u-x1l);
    dang2 = dang2 / (x2u-x2l);

    *energyOut = energy;
    *dang1Out = dang1;
    *dang2Out = dang2;
}