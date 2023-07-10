#include <bits/stdc++.h>
using namespace std;

class Edge;
class Face;

class Vertex{
    public:
    int id;
    double x, y;
    Edge *inc_edge;
};

class Edge{
    public:
    Edge *twin;
    Vertex *org;
    Face *left; //incident face
    Edge *next;
    Edge *prev;


};

class Face{
    public:
    Edge *inc_edge;

    vector<Vertex *> enumerate_vertices(Face *f);
};

class DCEL{
    public:
    vector<Vertex *> vertices;
    vector<Edge *> edges;
    vector<Face *> faces;
};

vector<Vertex *> Face:: enumerate_vertices(Face *f){
    Edge *start = f->inc_edge;
    Edge *e = start;
    vector<Vertex *> v;

    return v;
}

Vertex *next(vector<Vertex *> p, Vertex *v){
    for(int i=0; i<p.size(); i++){
        if(p[i] == v){
            if(i==p.size()-1) return p[0];
            else return p[i+1];
        }
    }
    return NULL;
}

bool check_reflex(Vertex* a,Vertex* b,Vertex* c){
    pair<double,double> ca = {a->x, a->y};
    pair<double,double> cb = {b->x, b->y};
    pair<double,double> cc = {c->x, c->y};
    double x1 = ca.first-cb.first;
    double y1 = ca.second-cb.second;
    double x2 = cc.first-cb.first;
    double y2 = cc.second-cb.second;
    double crossproduct = x1*y2-x2*y1;
    if(crossproduct < 0){
        return true;
    }
    else{
        return false;
    }
}

vector<Vertex *> find_notches(vector<Vertex *> p, vector<Vertex *> l){
    vector<Vertex *> notches;
    for(int i=0; i<p.size(); i++){
        if(find(l.begin(), l.end(), p[i]) == l.end()){
            cout<<"not in p "<<i<<endl;
            if(i==0){
                if(check_reflex(p[p.size()-1], p[i], p[i+1])) {
                    cout<<"notch "<<p[i]->x<<" "<<p[i]->y<<endl;
                    notches.push_back(p[i]);
                }
            }
            else if(i==p.size()-1){
                if(check_reflex(p[i-1], p[i], p[0])) {
                    cout<<"notch "<<p[i]->x<<" "<<p[i]->y<<endl;
                    notches.push_back(p[i]);
                }
            }
            else{
                if(check_reflex(p[i-1], p[i], p[i+1])) {
                    cout<<"notch "<<p[i]->x<<" "<<p[i]->y<<endl;
                    notches.push_back(p[i]);
                }
            }
        }
    }
    return notches;
}

bool v_in_polygon(vector<Vertex*> l, Vertex* v){
    int n = l.size();
    bool flag = false;
    for(int i=0; i<n; i++){
        if(i==n-1){
            if(!check_reflex(l[0], l[i], v)){
                cout<<"angle not reflex: "<<l[i]->x<<" "<<l[i]->y<<endl;
                flag = true;
                return false;
            }
        }
        else{
            if(!check_reflex(l[i]->inc_edge->next->org, l[i], v)){
                cout<<"angle not reflex: "<<l[i]->x<<" "<<l[i]->y<<endl;
                flag = true;
                return false;
            }
        }
    }
    return true;
}

vector<Vertex *> vs_in_semiplane(vector<Vertex *> l, Vertex *v){
    //draw a line through v1 and v
    //check if vn and vx lie on same side
    vector<Vertex *> vtr;
    pair<double,double> v1 = {l[0]->x,l[0]->y};
    pair<double,double> vn = {l[l.size()-1]->x, l[l.size()-1]->y};
    pair<double,double> vp = {v->x, v->y};
    //coefficients of line drawn by vp and v1
    double a = vp.second - v1.second; //a=y2-y1
    double b = vp.first - v1.first; //b=x2-x1
    double c = a*v1.first + b*v1.second; //c=ax1+by1
    for(int i=l.size()-1; i>1; i--){
        pair<double,double> vx = {l[i]->x, l[i]->y};
        double one = a*vn.first + b*vn.second - c;
        double two = a*vx.first + b*vx.second - c;
        if(one*two>=0){
            //same side. remove
            vtr.push_back(l[i]);
        }
        else{
            //not same side. stop. ok doubt
            break;
        }
    }
    for(int i=0; i<vtr.size(); i++){
        cout<<"vertex not in semiplane: "<<vtr[i]->x<<" "<<vtr[i]->y<<endl;
    }
    return vtr;
}

void generate_file(vector<Edge *> edges){
    fstream fp;
    fp.open("out1.txt", ios::out); 
    if(fp.is_open()){
        int n = edges.size();
        for(int i=0; i<n; i++){
            // cout<<dcel->edges[i]->org->x<<" "<<dcel->edges[i]->org->y<<" "<<dcel->edges[i]->twin->org->x<<" "<<dcel->edges[i]->twin->org->y<<endl;
            fp << edges[i]->org->x<<" "<<edges[i]->org->y<<" "<<edges[i]->twin->org->x<<" "<<edges[i]->twin->org->y;
            if(i!=n-1) fp<<endl;
        }
        fp.close();
    }
}

DCEL* make_dcel(vector<Vertex *> v){
    int n = v.size();
    DCEL *dcel = new DCEL();
    for(int i=0; i<n; i++){
        dcel->vertices.push_back(v[i]);
        dcel->edges.push_back(v[i]->inc_edge);
    }
    dcel->edges.pop_back();                     // removing the last edge as it will be different in this case
    Edge *e = new Edge(), *et = new Edge();
    e->org = v[0];
    e->prev = v[0]->inc_edge->prev;
    e->next = v[v.size()-1]->inc_edge;
    e->twin = et;
    et->twin = e;
    et->org = v[v.size()-1];
    et->prev = v[v.size()-1]->inc_edge->twin;
    et->next = v[0]->inc_edge->twin->next;
    dcel->edges.push_back(e);
    return dcel;
}

void decompose(DCEL *dcel){
    vector<DCEL *> dcels;
    dcels.push_back(dcel);
    vector<Vertex *> p;                             //p contains all the vertices in the polygon in a clockwise order
    vector<Edge *> ret;                             //stores all the edges for plotting the polygon
    for(int i=0; i<dcel->vertices.size(); i++){
        p.push_back(dcel->vertices[i]);
        ret.push_back(dcel->edges[i]);
        // cout<<p[i]->x<<" "<<p[i]->y<<endl;
    }
    vector<vector<Vertex *>> l;                     //l contains the partitions of our convex polygons - from the 1st index
    l.push_back({p[0]});
    int m = 1;
    while(p.size() > 3){
        Vertex *v1 = l[m-1][l[m-1].size()-1];
        Vertex *v2 = next(p, v1);
        l.push_back({});
        l[m].push_back(v1);
        l[m].push_back(v2);
        cout<<"hi "<<v1->x<<" "<<v1->y<<endl<<"hi "<<v2->x<<" "<<v2->y<<endl;
        int t = 2;
        Vertex *v3 = next(p, v2);
        Vertex *vim = v1, *vi = v2, *vip = v3;
        //checking if the vertices satisfy all 3 conditions to be in the polygon
        while(!check_reflex(vim, vi, vip) && !check_reflex(vi, vip, v1) && !check_reflex(vip, v1, v2) && l[m].size()<p.size()){ 
            cout<<"hi "<<vip->x<<" "<<vip->y<<endl;
            l[m].push_back(vip);
            t++;
            vim = vi;
            vi = vip;
            vip = next(p, vip);                                 //next(p, vip) to get the next vertex
            // cout<<l[m].size()<<endl;
        }
        cout<<"hi1"<<endl;
        if(l[m].size() != p.size()){                            //check if l[m].size == p.size
            vector<Vertex *> lpvs = find_notches(p, l[m]);      //get the list of notches in the polygon
            cout<<"no. of notches: "<<lpvs.size()<<endl;
            while(lpvs.size() > 0){
                double minx = FLT_MAX, maxx = -FLT_MAX, miny = FLT_MAX, maxy = -FLT_MAX;
                for(int j=0; j<l[m].size(); j++){
                    // cout<<"x:"<<l[m][j]->x<<endl;
                    if(l[m][j]->x > maxx){
                        maxx = l[m][j]->x;
                    }
                    if(l[m][j]->x < minx){
                        minx = l[m][j]->x;
                    }
                    if(l[m][j]->y > maxy){
                        maxy = l[m][j]->y;
                    }
                    if(l[m][j]->y < miny){
                        miny = l[m][j]->y;
                    }
                }
                cout<<"x: "<<minx<<" "<<maxx<<endl<<"y: "<<miny<<" "<<maxy<<endl;
                bool backward = false;
                while(!backward && lpvs.size()>0){
                    Vertex *v = lpvs[0];
                    while(v->x < minx || v->x > maxx || v->y < miny || v->y > maxy){
                        cout<<"removing notch"<<endl;
                        lpvs.erase(lpvs.begin()); // check
                        if(lpvs.size() == 0) break;
                        v = lpvs[0];
                    }
                    // exit(0);
                    if(lpvs.size() > 0){
                        cout<<"notches in polygon.."<<endl;
                        if(v_in_polygon(l[m], v)){
                            cout<<" in notch"<<endl;
                            //exit(0);
                            vector<Vertex *> vtr = vs_in_semiplane(l[m], v);
                            for(int j=0; j<vtr.size(); j++){
                                int d = find(l[m].begin(), l[m].end(), vtr[j]) - l[m].begin();
                                cout<<"removing vertex: "<<l[m][d]->x<<" "<<l[m][d]->y<<endl;
                                l[m].erase(l[m].begin()+d);
                            }
                            backward = true;
                        }
                        lpvs.erase(lpvs.begin());
                    }
                } 
            }
        }    
        if(l[m][l[m].size()-1] != v2){                     //it is not a line
            //add edge joining v1 and last of l
            cout<<"joining edges: "<<l[m][0]->x<<" "<<l[m][0]->y<<" & "<<l[m][l[m].size()-1]->x<<" "<<l[m][l[m].size()-1]->y<<endl;
            // dcel->edges.push_back(e);
            DCEL *nw = make_dcel(l[m]);
            ret.push_back(nw->edges[nw->edges.size()-1]);
            dcels.push_back(nw);
            //remove all vertices of l[m] from p except for the first and last one
            for(int j=1; j<l[m].size()-1; j++){ 
                int d = find(p.begin(), p.end(), l[m][j]) - p.begin();
                p.erase(p.begin()+d);
            }
            // n = n - l[m].size() + 2;
        }       
        m++;         
    }
    cout<<dcels.size()<<endl;
    generate_file(ret);
}

void vertex_dependency(DCEL *dcel){
    vector<DCEL *> os;
    int card = INT_MAX, s = 1, n = dcel->vertices.size();
    vector<DCEL *> dcels;
    while(s <= n){
        dcels.push_back(new DCEL());
        for(int i=s-1; i<n; i++){

        }
    }
}

int main(){
    //get the DCEL struct of our polygon
    DCEL *dcel = new DCEL();
    int n;
    vector<Vertex *> vertices;
    vector<Edge *> edges;
    fstream fp;
    fp.open("input.txt", ios::in); 
    if (fp.is_open()){
        fp >> n;
        cout<<n<<" "<<endl;
        for(int i=0; i<n; i++){
            Vertex *v = new Vertex();
            Edge *e = new Edge(), *et = new Edge();
            double x, y;
            fp >> x;
            fp >> y;
            cout<<x<<" "<<y<<endl;
            v->id = i+1;
            v->x = x;
            v->y = y;
            v->inc_edge = e;
            e->org = v;
            e->twin = et;
            et->twin = e;
            vertices.push_back(v);
            edges.push_back(e);
        }
        fp.close();
    }
    for(int i=0; i<n; i+=1){
        if(i!=n-1){
            edges[i]->next = edges[i+1];
            edges[i]->twin->org = vertices[i+1];
            edges[i]->twin->prev = edges[i+1]->twin;
        }
        else{
            edges[i]->next = edges[0];
            edges[i]->twin->org = vertices[0];
            edges[i]->twin->prev = edges[0]->twin;
        }
        if(i!=0){
            edges[i]->prev = edges[i-1];
            edges[i]->twin->next = edges[i-1]->twin;
        }
        else{
            edges[i]->prev = edges[n-1];
            edges[i]->twin->next = edges[n-1]->twin;
        }
    }
    dcel->vertices = vertices;
    dcel->edges = edges;
    //generate_file(dcel->edges);
    decompose(dcel);
    // decompose1(dcel);
    return 0;
}