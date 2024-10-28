#ifdef GL_ES
precision mediump float;
#endif

varying vec2 v_texCoord;
varying vec4 v_fragmentColor;
uniform int u_npts;           // 贝塞尔曲线段数
uniform float u_pts[400];     // 每段贝塞尔曲线的点和控制点数组，400为最大长度，可根据需求调整

// uniform float u_posx;
// uniform float u_posy;

float cubicRoot(float x) 
{
    return x >= 0.0 ? pow(x, 1.0 / 3.0) : -pow(-x, 1.0 / 3.0);
}

int soluteCubicBezier(float posx, float posy, float x0, float y0, float cx1, float cy1, float cx2, float cy2, float x3, float y3)
{
    float a = y3 - 3.0 * cy2 + 3.0 * cy1 - y0;
    float b = 3.0 * cy2 - 6.0 * cy1 + 3.0 * y0;
    float c = 3.0 * cy1 - 3.0 * y0;
    float d = y0 - posy;

    if (abs(a) < 1e-6) 
    { // 二次或一次方程
        if (abs(b) < 1e-6) 
        {
            if (abs(c) < 1e-6)
                return 0;
            float t = -d / c;
            float u = 1.0 - t;
            float x = u * u * u * x0 + 3.0 * u * u * t * cx1 + 3.0 * u * t * t * cx2 + t * t * t * x3;
            return (x >= posx && t >= 0.0 && t <= 1.0) ? 1 : 0;
        }

        float delta = c * c - 4.0 * b * d;
        if (delta > 0.0) 
        {
            float deltaSqrt = sqrt(delta);
            float t1 = (-c - deltaSqrt) / (2.0 * b);
            float t2 = (-c + deltaSqrt) / (2.0 * b);
            int roots = 0;
            if (t1 >= 0.0 && t1 <= 1.0) 
            {
                float u1 = 1.0 - t1;
                if (u1 * u1 * u1 * x0 + 3.0 * u1 * u1 * t1 * cx1 + 3.0 * u1 * t1 * t1 * cx2 + t1 * t1 * t1 * x3 >= posx)
                    roots++;
            }
            if (t2 >= 0.0 && t2 <= 1.0) 
            {
                float u2 = 1.0 - t2;
                if (u2 * u2 * u2 * x0 + 3.0 * u2 * u2 * t2 * cx1 + 3.0 * u2 * t2 * t2 * cx2 + t2 * t2 * t2 * x3 >= posx)
                    roots++;
            }
            return roots;
        }
        return 0;
    }

    // 三次求解
    //将
    //at^3+bt^2+ct+d=0
    //化简为
    //y^3+py+q=0
    //其中 t=y-b/(3a)
    b /= a;
    c /= a;
    d /= a;
    float p = (3.0 * c - b * b) / 3.0;
    float q = (2.0 * b * b * b - 9.0 * b * c + 27.0 * d) / 27.0;
    
    float half_q = q / 2.0;
    float delta = half_q * half_q + (p * p * p) / 27.0;

    if (delta >= 0.0) 
    {
        float A = cubicRoot(-half_q + sqrt(delta));
        float B = cubicRoot(-half_q - sqrt(delta));
        
        //float y = A + B 
        //float t = y - b/3.0;//这里b==原方程中的b/a
        float t = A + B - b/3.0;
        
        float u = 1.0 - t;
        float x = u * u * u * x0 + 3.0 * u * u * t * cx1 + 3.0 * u * t * t * cx2 + t * t * t * x3;
        return (x >= posx && t >= 0.0 && t <= 1.0) ? 1 : 0;
    } 
    else 
    {

        float r = 2.0 * sqrt(-p / 3.0);
        float theta = acos(3.0 * q / (2.0 * p) * sqrt(-3.0 / p));
        int roots = 0;
        float ts[3];
        ts[0]=r * cos(theta / 3.0) - b / 3.0;
        ts[1]=r * cos((theta + 2.0 * 3.14159265359) / 3.0) - b / 3.0; 
        ts[2]=r * cos((theta + 4.0 * 3.14159265359) / 3.0) - b / 3.0;

        for (int i = 0; i < 3; ++i) 
        {
            float t = ts[i];
            if (t >= 0.0 && t <= 1.0) 
            {
                float u = 1.0 - t;
                float x = u * u * u * x0 + 3.0 * u * u * t * cx1 + 3.0 * u * t * t * cx2 + t * t * t * x3;
                if (x >= posx) roots++;
            }  
        }
        return roots;
    }
}


int mod(int x,int y)
{
    return x-x/y*y;
}

void main() 
{
    float scaleFactor=2.0;
    vec2 origin=vec2(200.0,600.0);
    float scale=4.0;

    vec2 factpos = vec2(gl_FragCoord.x/scaleFactor,gl_FragCoord.y/scaleFactor);
    vec2 pos = factpos-origin;
    pos=vec2 (pos.x/scale,pos.y/scale);

    // gl_FragColor=vec4(pos.x-u_posx,pos.y-u_posy,0.0,0.5);

    int   num=0;
    for(int i=0;i<u_npts*2-7;i+=6)
    {
        num += soluteCubicBezier(pos.x,pos.y,u_pts[i],-u_pts[i+1],u_pts[i+2],-u_pts[i+3],u_pts[i+4],-u_pts[i+5],u_pts[i+6],-u_pts[i+7]);
    }

    if (mod(num,2)==1) {
        gl_FragColor=vec4(1.0,1.0,1.0,1.0);
        // gl_FragColor = v_fragmentColor * texture2D(CC_Texture0,v_texCoord);
    } else {
        discard;
    }
}