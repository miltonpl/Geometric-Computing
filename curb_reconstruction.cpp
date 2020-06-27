
/*
    Milton Palaguachi
    Geometric Computing
    N = max number of points = 132;
    3 main function 
    Triangulation computes  initial triangulation
    Delauney Disc computed Delauney triangulation
    Curb Reconstruction recontructs the curb. 
    Note: My Delauney Disc Function is a bit slow since I did not remove the outer edges 
    and the while loop check for all point of edges all over when finds a violation.
*/
// run this on the terminal :
// 1.  c++ test1.cpp -o test1 -lX11 -lm -L/usr/X11R6/lib
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<iostream>
#include <X11/Xlib.h>
#include<algorithm>
#include<math.h>
#include <stack>
#define N 132// number of pair points
typedef struct Edge Edge;
typedef struct Point{
        float x, y;        
        int n;
        // if point left and point are equal return true else false
         friend bool operator ==(const Point& lhs, const Point& rhs)
        {
            return (lhs.x == rhs.x && lhs.y ==rhs.y) ? true : false;
        }
        } Point;
struct Triangle{
    // from the 3 points of the triangle
    Point pi,pj,pk;
    int n;
    // circumcenter coordinates of the circunference defined
    // float xc , yc;
    Point center;
    // circunference radius
    float radius;
};//edge
struct Edge{
    Point start,dest;
    int n;
    Point left,right;
    Triangle *adj_r,*adj_l;
    bool visited;
    friend bool operator ==(const Edge& lhs, const Edge& rhs){ 
        return ((lhs.start == rhs.start && lhs.dest == rhs.dest) || (lhs.start == rhs.dest && rhs.dest == rhs.start))? true: false;
    }
    
};//edge

// Orientation ABO: two vectors OA and OB 
// returns positive for counter clockwise 
// turn and negative for clockwise turn 
//zero if coliniear
float orientation(Point O, Point A, Point B) { 
	return ((A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x)); 
} 
void swap_p(Point *p1,Point *p2){
    Point temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}
void bubbleSortx(std::vector<Point> &points){
    int i,j,n=points.size();

    for (i=0; i< n-1;i++){
        //last i element are already in place
        for(j= 0; j< n-i-1; j++){
            if(points[j].x > points[j+1].x)
                swap_p(&points[j],&points[j+1]);
        }
    } 
}
void bubbleSorty(std::vector<Point> &points){
    int i,j,n=points.size();

    for (i=0; i< n-1;i++){
        //last i element are already in place
        for(j= 0; j< n-i-1; j++){
            if(points[j].y > points[j+1].y)
                swap_p(&points[j],&points[j+1]);
        }
    } 
}
//Calculates the radius and the circumcenter location
/* 
    Finding Circumcentr Using Linear Equations
    let assume (X,Y) be the coordinates of the Circumcenter.
    According to the cicumcenter properties the distance of(X,Y) form each vertex 
    will be the same
    D1= sqrt((X-A.x)^2 + ((Y-A.y)^2)
    D2= sqrt((X-B.x)^2 + ((Y-B.y)^2)
    D3= sqrt((X-C.x)^2 + ((Y-C.y)^2)
    then we have D1=D2=D3
    solve for X,Y using equation 
    D1= D2 q1
    D1=D3 q2
    where X=cm.x  snd Y = cm.y

*/
Point circumCenter(Point A, Point B, Point C,float &r){
    float a,b,c,d,l,m,dx,dy;
    Point cm;
    //A(x1,y1) B(x2,y2) C(x3,y3)
    //Xa +Yb = l (Eq. 1) // D1==D2
    a = (B.x - A.x); b = (B.y - A.y);
    l = ((B.x*B.x + B.y*B.y) - (A.x*A.x + A.y*A.y))/2.0;
     //Xc +Yd = m (Eq. 2)
    c = (C.x - A.x); d = (C.y - A.y);
    m = ((C.x*C.x + C.y*C.y) - (A.x*A.x + A.y*A.y))/2.0;
    
    //solve for Y Eq2 and Eq1  //D1==D3
    cm.y=(m*a -l*c)/(a*d - b*c);
    // solve for X=(m -Yd)/c
    cm.x = (m -cm.y*d)/c;
     // calculates radius
    dx = A.x - cm.x;
    dy = A.y - cm.y;
    r = sqrt(dx*dx + dy*dy);
    // returns point(x,y) of the circum center
    return cm;
}
bool edgeInTriangle(Edge e, Triangle t){
    Point a= t.pi, b=t.pj, c=t.pk;
    //edg e = e.start -> e.destinatin
    //check t=<a,b,c>
    // (e.s->,e.d) == (a->b) or  (e.s->,e.d) == (b->a)
    if((e.start == a && e.dest == b) || (e.start == b && e.dest == a))
        return true;
    // (e.s->,e.d) == (a->c) or (e.s->,e.d) == (c->a)
    if((e.start == a && e.dest == c) || (e.start == c && e.dest == a))
        return true;
    // (e.s->,e.d) == (b->c) or (e.s->,e.d) == (c->b)
    if((e.start == b && e.dest == c)||(e.start == c && e.dest == b))
        return true;
    
    return false;
}
//get all not reapeated edges of the triangle
std::vector<Edge> edgeReconstruction(std::vector<Triangle> T){
    std::vector<Edge>edgs;
    Edge e;
    e.adj_l = NULL; e.adj_r = NULL; e.visited =false;
    
    //initialization of edges
     for (int j=0; j < T.size();j++)
    {
        e.start=T[j].pi; e.dest = T[j].pj;
        edgs.push_back(e);//edge 1

        e.start=T[j].pj; e.dest = T[j].pk;
        edgs.push_back(e);//edge 2

        e.start = T[j].pi; e.dest = T[j].pk;
        edgs.push_back(e); //edge3
    }
    //Remove repeated edges
    for(int s = 0; s < edgs.size(); s++){
        for (size_t i = s+1; i < edgs.size(); i++)
        {
            if(edgs[s] == edgs[i])
            edgs.erase(edgs.begin()+i);
        }
    }
    // id for each edge
    for (size_t i = 0; i < edgs.size(); i++) edgs[i].n = i;
    return edgs;
}
//get adjacen triangle of the vec_edges
void adjacent_tringle(std::vector<Edge> &vec_edgs,std::vector<Triangle> &T){
    //return all not repeated edges of triangle
    vec_edgs = edgeReconstruction(T);
   
    for (size_t i = 0; i < T.size(); i++) T[i].n=i;    //triangle id    
    //triangles that share the same edge
    for (size_t j = 0; j < vec_edgs.size(); j++)
    {
        int adj_count=0;
        for (size_t i = 0; i < T.size(); i++)
        {
            if(edgeInTriangle(vec_edgs[j], T[i]))
            {
                //An Edge can share max 2 triangles 
                if(adj_count == 0){
                    //right adjacent Triangle
                    adj_count++;
                    vec_edgs[j].adj_r = &T[i];
                }
                else if(adj_count == 1){
                    adj_count++;
                    //left adjacent Triangle
                    vec_edgs[j].adj_l = &T[i];
                }
                else// if counter is greater than 1 that means that there are peated triangles therefor therefor there is a violation
                    printf("Repeated triangles\n not allowed!!%d\n",adj_count);
            } 
    }  }
}
//update triangle funciton
void updateTriangle(std::vector<Triangle> &T,Edge edge,Point pc,Point pd){
    Triangle right,left;
     //left triangle new vertices
    left.n = edge.adj_l->n;// copy the index
    left.pi = pc;
    left.pj = pd;
    left.pk = edge.start;
    //get cirmcumCenter of left triangle and its radius
    left.center = circumCenter(pc, pd, edge.start, left.radius);

    //right triangle new vertices;
    right.n = edge.adj_r->n;//copy the index 
    right.pi = pc;
    right.pj = pd;
    right.pk = edge.dest;
    // get circumCenter of right;
    right.center = circumCenter(pc, pd, edge.dest, right.radius);
    //update new triangle
    T[left.n] =left;
    T[right.n] = right;
}
float get_distance(Point p, Point q){
    float x,y;
    x = p.x - q.x;
    y = p.y - q.y;
    return sqrt(x*x + y*y);
}
std::vector<Triangle> Delauney_Disk(std::vector<Triangle> &T){

    std::vector<Edge>vec_edgs;
    Point center1;
    float radius_left,radius_right;
    Triangle *l,*r;
    Edge e;
    //get the circum cencet (x,y) of all triangles
    for (size_t i = 0; i < T.size(); i++)
        T[i].center = circumCenter(T[i].pi, T[i].pj, T[i].pk,T[i].radius);
        //circumcenter
    adjacent_tringle(vec_edgs, T);
    bool check=false;
    while(!check)
    {
        check=true;
        for(size_t i=0; i < vec_edgs.size();i++)
            {   Point p_c,p_d;
                e = vec_edgs[i];//sharing edge
                r = e.adj_r; l = e.adj_l;
                // right(Triangle) <- edge ->left(Trinagle)
                // if the triangle has a right and left adjacent triangle
                if(r != NULL && l != NULL){
                    // Find point in the right triangle that is opposite to sharing edge
                    //if the point pi or pj or pk is not in edge e then the point is on the opposite  to current edge
                    if(r->pi.n != e.start.n && r->pi.n != e.dest.n ) 
                        p_c = r->pi;
                    else if(r->pj.n != e.start.n && r->pj.n != e.dest.n ) 
                        p_c = r->pj;
                    else 
                        p_c = r->pk;
                        //Find point in left triangle that is opposite to edge
                    if(l->pi.n != e.start.n && l->pi.n != e.dest.n ) 
                        p_d = l->pi;
                    else if(l->pj.n != e.start.n && l->pj.n != e.dest.n ) 
                        p_d = l->pj;
                    else 
                        p_d = l->pk;

                    // check point of the right triagle if it is inside left triangle;
                    radius_left = get_distance(l->center, p_c);
                    //distance between the point of the right triangle and center of the triagle
                    radius_right = get_distance(r->center, p_d);
                    //if point p_c or p_d is inside triangle then there is violation of Delauny Disk
                    if(radius_left < l->radius || radius_right < r->radius){
                        updateTriangle(T, vec_edgs[i],p_c,p_d);
                        vec_edgs.clear();//remove all element from a vector
                        adjacent_tringle(vec_edgs,T);//get new edges from updated Triangle
                        check =false;
                    }
                }
            }
    }
    return T;
}
// the mid point of two points ex. a-----m----b
Point bisector(Point a,Point b,float &m){
    Point M;
    float slope;
    M.x = (a.x + b.x)/2;
    M.y = (a.y + b.y)/2;
    //if Point a and Point be form a orizontal or vertical line than its slope is zero or undinied
    if( (a.y - b.y) == 0 ||(a.x - b.x) == 0 ){
        m = 0;
        return M;
    }
    else{
          slope = (a.y - b.y)/(a.x - b.x);
            m = -1/slope;
        return M;
    }
}
std::vector<Edge> curb_recruction(std::vector<Triangle> T){
    //get edges from set of Trinagles 
    std::vector<Edge> curbe,edges = edgeReconstruction(T);
    int counter =0;
    for (size_t i = 0; i < edges.size(); i++)
    {
        Edge e = edges[i];
        Point p,q,s;
        p = edges[i].start; 
        q = edges[i].dest;

        float old_d,new_d;
        // edge p->q
        //distand p->q
        old_d = get_distance(p,q);

        //index have all edges id/index that contain point p
        std::vector<int>index;
        index.push_back(e.n);
        //find all the edges that share the point p 
        for (size_t j = 0; j < edges.size(); j++)
        {
            if (i != j){
                if(p.n == edges[j].start.n || p.n == edges[j].dest.n){
                    index.push_back(edges[j].n);
                    //calculate distance from p to q
                     new_d = get_distance(edges[j].start,edges[j].dest);
                    //update distance/q if new distance is less than old distance
                    if(new_d < old_d ){
                        old_d = new_d;
                        q =(p.n == edges[j].start.n)? edges[j].dest : edges[j].start;
                    }
                }   
            }
        }
        /*At this point we have the shorstes distance of vertex p that is connected to point q.*/

        //new_d for distance, m for bisector slop; same_orientaion for if al least one point on the same orientation
        float m;// m is the slope
        bool same_orientation=false;
        old_d = INT_MAX;
        Point W,Z,S,midPoint;
        midPoint = bisector(p, q, m);
        //if slope is zero //let W.x=M.x and W.y= and check for orientation.
        //same p.y == q.y then vertical line slope = 0
        if((p.y - q.y) == 0){
            W.y = 0;
            W.x = midPoint.x;
        }
        else if((p.x - q.x) == 0){// slope undifined p.x == p.x vertical line
            W.y = midPoint.y;;
            W.x = 0;
        }
        else{//slope does not equal to zero or undefine
              // y = y1 + m*(0 - x1);
              W.y= midPoint.y -(m*midPoint.x);
              W.x=0;
        }
        //assume s initial distanc d = infinity
        for (size_t j = 0; j < index.size(); j++)
        {
            e = edges[index[j]];
            float orientation1,orientation2;// for orientation
            if(p.n == e.dest.n && q.n != e.start.n) 
                Z = e.start;
            else if(p.n == e.start.n && q.n != e.dest.n) 
                Z = e.dest;
            else 
                continue;// ignor whe p.n == e.start  and q.n == e.dest
            orientation1 =orientation(midPoint,W,p);
            orientation2 = orientation(midPoint,W,Z);
            //same orientation
            if((orientation1 > 0 && orientation2 > 0)||(orientation1 < 0 && orientation2 < 0)){
                same_orientation = true;
                //check distance:p->Z
                new_d = get_distance(p,Z);
                if(new_d < old_d){
                    old_d = new_d;
                    S = Z;
                }
            }
        }
        index.clear();        
        e.start = p; e.dest = q;
        curbe.push_back(e);//push p->q
        if(same_orientation){//find at least one at same orientation.
            e.dest = S;
            curbe.push_back(e);    //push p->s    
            }
    }
    return curbe;

}
std::vector<Triangle> triangulation(std::vector<Point> points){
    std::vector<Point> u_cov,l_cov;
    std::vector<Triangle> triangles;
    Triangle t;
    if(points.size() < 3)
        return triangles;//empty triangles return sine a triangle must have at leat 3 points
    bubbleSorty(points);
    bubbleSortx(points);
   
    if(orientation(points[0],points[1],points[2]) != 0){
         t.pi=points[0]; t.pj = points[1]; t.pk = points[2];
          triangles.push_back(t);
    }
    //CCW positive
    if (orientation(points[0], points[1], points[2]) < 0) {
        u_cov.push_back(points[0]);
        u_cov.push_back(points[1]);
        u_cov.push_back(points[2]);

        l_cov.push_back(points[0]);
        l_cov.push_back(points[2]);
    
    }else //negative CW
    if(orientation(points[0], points[1], points[2]) > 0){
        l_cov.push_back(points[0]);
        l_cov.push_back(points[1]);
        l_cov.push_back(points[2]);

        u_cov.push_back(points[0]);
        u_cov.push_back(points[2]);
    }
    else{
        //Collinear
        u_cov.push_back(points[0]);
        u_cov.push_back(points[1]);
        u_cov.push_back(points[2]);

        l_cov.push_back(points[0]);
        l_cov.push_back(points[1]);
        l_cov.push_back(points[2]);
    }
    Point last_point = points[2];
    for(int k = 3; k < points.size(); k++) {

        //when repeated points then go to the end of loop
        // if(last_point.x == points[k].x && last_point.y == points[k].y) continue;
        if(last_point == points[k])continue;
        int i = u_cov.size() - 1;
       // upper
        while(i > 0 && orientation(points[k], u_cov[i], u_cov[i - 1]) < 0){
            t.pi=points[k]; t.pj = u_cov[i]; t.pk=u_cov[i-1];
            triangles.push_back(t);
            u_cov.pop_back();
            i--;
        }
        i = l_cov.size() - 1;
        //lower
        while(i > 0 && orientation(points[k], l_cov[i], l_cov[i - 1]) > 0){
            t.pi=points[k]; t.pj = l_cov[i]; t.pk = l_cov[i-1];
            triangles.push_back(t);
            l_cov.pop_back();
            i--;
        }
        u_cov.push_back(points[k]);
        l_cov.push_back(points[k]);
        last_point = points[k];
    } 
    return triangles;
}
////////////////////////////////

Display *display;
Window  window;
int border_width;
int win_width, win_height;
int win_x, win_y;

XSetWindowAttributes attributes;
long valuemask = 0;
XGCValues gr_values,gc_yellow_values;
Colormap color_map;
XFontStruct *fontinfo;
XColor tmp_color1, tmp_color2;
GC gc_red,gc_yellow,gc_green;
Visual *visual;
int depth;
int screen;
XEvent event;
XColor    color, dummy;
////////////////////////////////
int main(int argc,char *argv[]){
     int n,finish;
     Point p;
    std::vector<Point> points,inner_points;  
    std::vector<Edge> Eds;  
    std::vector<Triangle> T;
    FILE * inputfile;
    if(argc != 2){
          printf("need file name as command line argument.\n");fflush(stdout);
                 exit(0);
          }
    inputfile= fopen(argv[1],"r");
    n = 0 ;finish = 0;
    while(n < N && !finish){
        if(fscanf(inputfile, "P (%f,%f)\n", &(p.x),&(p.y)) != 2)
            finish = 1;
        else
        {
        // place in 2D vec
            p.n = n;
            points.push_back(p);
            n++;
       }
    }
    T = triangulation(points);//simple triangulation  
    Delauney_Disk(T);// Delauney Triangulation
    std::vector<XPoint> pointss,others;
    std::vector<Edge> curb = curb_recruction(T);//Curb reconstruction.
    std::vector<Edge> edges = edgeReconstruction(T);
    int curb_n = curb.size();
    int edges_n = edges.size();
    XSegment seg_curb[curb_n];
    XSegment seg_edges[edges_n];

    for (size_t i = 0; i < curb_n; i++)
    {
        seg_curb[i].x1 = curb[i].start.x; seg_curb[i].y1 = curb[i].start.y;
        seg_curb[i].x2 = curb[i].dest.x; seg_curb[i].y2 = curb[i].dest.y;
    }
    
    for (size_t i=0; i < edges.size();i++){
        seg_edges[i].x1 = edges[i].start.x; seg_edges[i].y1 = edges[i].start.y;
        seg_edges[i].x2 = edges[i].dest.x;  seg_edges[i].y2 = edges[i].dest.y;
    }


    if((display = XOpenDisplay(NULL))== NULL){
          printf("Could not open display. \n"); exit(-1);}
        screen = DefaultScreen(display);
        visual = DefaultVisual(display,screen);
        depth  = DefaultDepth(display,screen);
        color_map= XDefaultColormap(display,screen);
        /* creating the window */
        border_width = 10;
        win_x = 0; win_y = 0;
        win_width= 500;
        win_height= 500;
                
        attributes.background_pixel = XWhitePixel(display,screen);

        window = XCreateWindow( display,XRootWindow(display,screen),
                                win_x, win_y, win_width, win_width, border_width,
                                depth,  InputOutput,
                                visual ,CWBackPixel, &attributes);
        XSelectInput(display,window,ExposureMask | KeyPressMask) ;
        fontinfo = XLoadQueryFont(display,"6x10");

        XAllocNamedColor(display, DefaultColormap(display, screen),"red",&color,&dummy);

        gr_values.font = fontinfo->fid;
        gr_values.foreground = color.pixel;
        gc_red=XCreateGC(display,window,GCFont+GCForeground, &gr_values);
        gc_yellow = XCreateGC( display, window, valuemask, &gc_yellow_values);
        if( XAllocNamedColor( display, color_map, "yellow", &tmp_color1, &tmp_color2 ) == 0 )
        {printf("failed to get color yellow\n"); exit(-1);} 
        else
        XSetForeground( display, gc_yellow, tmp_color1.pixel );

        gc_green = XCreateGC( display, window, valuemask, &gc_yellow_values);
        if( XAllocNamedColor( display, color_map, "green", &tmp_color1, &tmp_color2 ) == 0 )
        {printf("failed to get color yellow\n"); exit(-1);} 
        else
        XSetForeground( display, gc_green, tmp_color1.pixel );

        XFlush(display);
        XMapWindow(display,window);
        XFlush(display);

        while(1){
            XNextEvent(display,&event);

            switch(event.type){
            case Expose:
                XFlush(display);
                XDrawSegments(display, window, gc_red, seg_edges, edges_n);//Delauney triangulation
                //  XDrawSegments(display, window, gc_green, seg_curb, curb_n);//curb
                XFlush(display);
                break;
            case KeyPress:
                XCloseDisplay(display);
                return 1;

            }
    }
    
}
