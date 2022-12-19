#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define max_polygon_size 100
#define sites_cnt 30
#define eps 0.00000001


#define format_point "(%f, %f)"
#define r() ((double)rand() / RAND_MAX)

typedef struct{
    double a, b, c;
}line;

typedef struct{
    double x, y;
}point2d;

typedef struct{
    int size;
    point2d points[max_polygon_size];
}polygon;

double get_signed_dist(line l, point2d p){
    return l.a * p.x + l.b * p.y + l.c;
}

double dot2d(point2d p1, point2d p2){
    return p1.x*p2.x + p1.y*p2.y;
}

line get_bisection(point2d p1, point2d p2){
    double a = -(p2.x - p1.x), b = -(p2.y - p1.y);
    point2d mid_p = {(p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2}; 
    double c = -(a*mid_p.x + b*mid_p.y);

    line l = {a, b, c};
    double d = get_signed_dist(l, p1);

    if(d < -eps){
        l.a *= -1.0; l.b *= -1.0; l.c *= -1.0;
    }

    return l;
}

point2d seg_to_line_intersection(point2d p1, point2d p2, line l){
    point2d rd = {p2.x - p1.x, p2.y - p1.y};
    point2d n = {l.a, l.b};
    double t = -(l.c + dot2d(p1, n)) / dot2d(rd, n);

    point2d res = {p1.x + t*rd.x, p1.y + t*rd.y};

    return res;
}

void cut_negative(polygon *polygons, int idx, line l){
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

void construct_voronoi(polygon *cells, point2d sites[]){
    for(int i = 0; i < sites_cnt; i++){
        cells[i].size = 4;

        cells[i].points[0].x = 0.0;  cells[i].points[0].y = 0.0;
        cells[i].points[1].x = 1.0;  cells[i].points[1].y = 0.0;
        cells[i].points[2].x = 1.0;  cells[i].points[2].y = 1.0;
        cells[i].points[3].x = 0.0;  cells[i].points[3].y = 1.0;

        sites[i].x = r(); sites[i].y = r();  

    }
    FILE *f = fopen("out.txt", "w");
    for(int i = 0; i < sites_cnt; i++){
        fprintf(f, format_point, sites[i].x, sites[i].y);
        fprintf(f, "\n");
    }

    for(int i = 0; i < sites_cnt; i++){
        for(int j = 0; j < sites_cnt; j++){
            if(i == j) continue;
            line bisec = get_bisection(sites[i], sites[j]);
            cut_negative(cells, i, bisec);
        }
    }

    for(int i = 0; i < sites_cnt; i++){
        fprintf(f, "\\operatorname{polygon} (");
        for(int j = 0; j < cells[i].size; j++){
            fprintf(f, format_point, cells[i].points[j].x, cells[i].points[j].y);
            fprintf(f, ",");
        }
        fprintf(f, format_point, cells[i].points[cells[i].size - 1].x, cells[i].points[cells[i].size - 1].y);
        fprintf(f, ")\n");
    }

    fclose(f);

}

int main(){
    srand(time(NULL));
    point2d sites[sites_cnt];
    polygon *Voronoi_cells = malloc(sizeof(polygon) * sites_cnt);

    construct_voronoi(Voronoi_cells, sites);    
    free(Voronoi_cells);
    return 0;
}