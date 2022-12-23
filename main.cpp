#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <graphics.h>

#define max_polygon_size 15
#define sites_cnt 100
#define eps 0.00000001

#define width 640
#define height 480

#define format_point "(%f, %f)"
#define r() ((double)rand() / RAND_MAX)

typedef struct{
    double a, b, c;
}Line;

typedef struct{
    double x, y;
}point2d;

typedef struct{
    int size;
    point2d points[max_polygon_size];
}polygon;

double get_signed_dist(Line l, point2d p){
    return l.a * p.x + l.b * p.y + l.c;
}

double dot2d(point2d p1, point2d p2){
    return p1.x*p2.x + p1.y*p2.y;
}

Line get_bisection(point2d p1, point2d p2){
    double a = -(p2.x - p1.x), b = -(p2.y - p1.y);
    point2d mid_p = {(p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2}; 
    double c = -(a*mid_p.x + b*mid_p.y);

    Line l = {a, b, c};
    double d = get_signed_dist(l, p1);

    if(d < -eps){
        l.a *= -1.0; l.b *= -1.0; l.c *= -1.0;
    }

    return l;
}

point2d seg_to_line_intersection(point2d p1, point2d p2, Line l){
    point2d rd = {p2.x - p1.x, p2.y - p1.y};
    point2d n = {l.a, l.b};
    double t = -(l.c + dot2d(p1, n)) / dot2d(rd, n);

    point2d res = {p1.x + t*rd.x, p1.y + t*rd.y};

    return res;
}

void cut_negative(polygon *polygons, int idx, Line l){
    int sz = polygons[idx].size, new_sz = 0;
    polygon buffer;
    for (int i = 0; i < sz; i++){
        double d1 = get_signed_dist(l, polygons[idx].points[i]), d2 = get_signed_dist(l, polygons[idx].points[(i + 1) % sz]);
        if (d1 > eps)
            buffer.points[new_sz++] = polygons[idx].points[i];
        else if(fabs(d1) < eps)
            continue;

        if (fabs(d2) < eps){
            buffer.points[new_sz++] = polygons[idx].points[(i + 1) % sz];
        }else if(d1 > eps && d2 < -eps || d1 < -eps && d2 > eps){
            point2d intersection_point = seg_to_line_intersection(polygons[idx].points[i], polygons[idx].points[(i + 1) % sz], l);
            buffer.points[new_sz++] = intersection_point;
        }
    }

    for(int i = 0; i < new_sz; i++){
        polygons[idx].points[i] = buffer.points[i];
    } 
    polygons[idx].size = new_sz;
}

point2d on_circle(int cnt, int idx, double offset_x, double offset_y){
    double step = (2.0 * M_PI) / cnt;
    point2d pt = {cos(step * idx) + offset_x, sin(step * idx) + offset_y};
    return pt;
}

void construct_voronoi(polygon *cells, point2d sites[], bool first_call){
    double step = (2.0 * M_PI) / sites_cnt;
    for(int i = 0; i < sites_cnt; i++){
        cells[i].size = 4;

        cells[i].points[0].x = 0.0;  cells[i].points[0].y = 0.0;
        cells[i].points[1].x = 1.0;  cells[i].points[1].y = 0.0;
        cells[i].points[2].x = 1.0;  cells[i].points[2].y = 1.0;
        cells[i].points[3].x = 0.0;  cells[i].points[3].y = 1.0;

        // sites[i] = on_circle(sites_cnt, i, 0.5, 0.5);
        if(first_call){
            sites[i].x = r();
            sites[i].y = r();
        }

    }

    for(int i = 0; i < sites_cnt; i++){
        for(int j = 0; j < sites_cnt; j++){
            if(i == j) continue;
            Line bisec = get_bisection(sites[i], sites[j]);
            cut_negative(cells, i, bisec);
        }
    }

}

void Lloyd_relaxation(polygon *cells, point2d sites[]){
    for(int i = 0; i < sites_cnt; i++){
        point2d centroid = {0.0, 0.0};
        for(int j = 0; j < cells[i].size; j++){
            centroid.x += cells[i].points[j].x; centroid.y += cells[i].points[j].y; 
        }
        centroid.x /= cells[i].size; centroid.y /= cells[i].size; 
        sites[i] = centroid;
    }
}


int main(int argc, char *argv[]){
    srand(time(NULL));
    point2d sites[sites_cnt];
    polygon *Voronoi_cells = (polygon *)malloc(sizeof(polygon) * sites_cnt);
    int colors[sites_cnt]; for (int i = 0; i < sites_cnt; i++){colors[i] = RGB((int)(r() * 255), (int)(r() * 255), (int)(r() * 255));}

    

    construct_voronoi(Voronoi_cells, sites, true); 
 
    int win = initwindow(width, height, "Voronoi");
    bool x = false;
    while(1){
        for(int i = 0; i < sites_cnt; i++){
            int *shape = (int*)malloc((Voronoi_cells[i].size + 1) * 2 * sizeof(int));
            for(int j = 0; j < Voronoi_cells[i].size; j++){
                point2d p1 = Voronoi_cells[i].points[j];
                shape[2 * j + 0] = (int)(p1.x * width);
                shape[2 * j + 1] = (int)(p1.y * height);

            }
            shape[2 * Voronoi_cells[i].size + 0] = shape[0];
            shape[2 * Voronoi_cells[i].size + 1] = shape[1];

            int col = colors[i];
            setcolor(col);
            setfillstyle(SOLID_FILL, col);
            fillpoly(Voronoi_cells[i].size + 1, shape);
            free(shape);
        }
        setcolor(RGB(255, 255, 255));
        for(int i = 0; i < sites_cnt; i++){
            circle((int)(sites[i].x * width), (int)(sites[i].y * height), 1);
        }
        swapbuffers();
        clearviewport();
        Sleep(1);
        Lloyd_relaxation(Voronoi_cells, sites);
        construct_voronoi(Voronoi_cells, sites, false);
    }

    free(Voronoi_cells);
    return 0;
}