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
    while(e!=start){
        // cout<<e->org->x<<" "<<e->org->y<<endl;
        v.push_back(e->org);
        e = e->next;
    }
    return v;
}

Vertex *next(DCEL *d, Vertex *v){
    // cout<<"in next";
    // cout<<d->vertices.size()<<endl;
    for(int i=0; i<d->vertices.size(); i++){
        // cout<<" hi2"<<d->vertices[i]->x<<" "<<d->vertices[i]->y<<endl;
        if(d->vertices[i]->x == v->x && d->vertices[i]->y == v->y){
            // cout<<" hi1"<<endl;
            if(i==d->vertices.size()-1) return d->vertices[0];
            else return d->vertices[i+1];
        }
    }
    return NULL;
}

Vertex *prev(DCEL *d, Vertex *v){
    for(int i=0; i<d->vertices.size(); i++){
        if(d->vertices[i]->x == v->x && d->vertices[i]->y == v->y){
            if(i==0) return d->vertices[d->vertices.size()-1];
            else return d->vertices[i-1];
        }
    }
    return NULL;
}

DCEL* make_dcel(vector<Vertex *> v){
    int n = v.size();
    DCEL *dcel = new DCEL();
    // cout<<n<<" :vertices size"<<endl;
    for(int i=0; i<n; i++){
        Vertex *nv = new Vertex();
        Edge *ne = new Edge(), *net = new Edge();
        nv->id = i+1;
        nv->x = v[i]->x;
        nv->y = v[i]->y;
        nv->inc_edge = ne;
        ne->org = nv;
        ne->twin = net;
        net->twin = ne;
        dcel->vertices.push_back(nv);
        dcel->edges.push_back(ne);
    }
    dcel->edges.pop_back();                     // removing the last edge as it will be different in this case
    Edge *e = new Edge(), *et = new Edge();
    dcel->vertices[n-1]->inc_edge = e;
    e->org = dcel->vertices[n-1];
    e->prev = dcel->vertices[n-1]->inc_edge->prev;
    e->next = dcel->vertices[0]->inc_edge;
    e->twin = et;
    et->twin = e;
    et->org = dcel->vertices[0];
    et->prev = dcel->vertices[0]->inc_edge->twin;
    et->next = dcel->vertices[n-1]->inc_edge->twin->next;
    dcel->edges.push_back(e);                   //last edge connects the diagonal in polygon
    for(int i=0; i<n; i+=1){
        if(i!=n-1){
            dcel->edges[i]->next = dcel->edges[i+1];
            dcel->edges[i]->twin->org = dcel->vertices[i+1];
            dcel->edges[i]->twin->prev = dcel->edges[i+1]->twin;
        }
        else{
            dcel->edges[i]->next = dcel->edges[0];
            dcel->edges[i]->twin->org = dcel->vertices[0];
            dcel->edges[i]->twin->prev = dcel->edges[0]->twin;
        }
        if(i!=0){
            dcel->edges[i]->prev = dcel->edges[i-1];
            dcel->edges[i]->twin->next = dcel->edges[i-1]->twin;
        }
        else{
            dcel->edges[i]->prev = dcel->edges[n-1];
            dcel->edges[i]->twin->next = dcel->edges[n-1]->twin;
        }
    }
    // for(int i=0; i<n; i++){
    //     cout<<dcel->vertices[i]->x<<"&"<<dcel->vertices[i]->y<<"&&"<<dcel->vertices[i]->inc_edge->next->org->x<<"&"<<dcel->vertices[i]->inc_edge->next->org->y<<"  ";
    // }
    // cout<<endl;
    return dcel;
}

DCEL *merge_dcel(DCEL *d1, DCEL *d2, Vertex *vs, Vertex *vt){
    vector<Vertex*> v;
    int n1 = d1->vertices.size(), n2 = d2->vertices.size();
    // for(int i=0; i<n1; i++){
    //     cout<<d1->vertices[i]->x<<"&"<<d1->vertices[i]->y<<"&&"<<d1->vertices[i]->inc_edge->next->org->x<<"&"<<d1->vertices[i]->inc_edge->next->org->y<<"  ";
    // }
    // cout<<endl;
    // for(int i=0; i<n2; i++){
    //     cout<<d2->vertices[i]->x<<"&"<<d2->vertices[i]->y<<"&&"<<d2->vertices[i]->inc_edge->next->org->x<<"&"<<d2->vertices[i]->inc_edge->next->org->y<<"  ";
    // }
    // cout<<"vs and vt: "<<vs->x<<"&"<<vs->y<<"  "<<vt->x<<"&"<<vt->y<<endl;
    for(int i=0; i<n1; i++){
        if(d1->vertices[i]->x == vs->x && d1->vertices[i]->y == vs->y){ //add all vertices till first common point
            break;
        }
        // cout<<d1->vertices[i]->x<<"&"<<d1->vertices[i]->y<<"  ";
        v.push_back(d1->vertices[i]);
    }
    // cout<<v.size()<<" ";
    int pos1 = -1;
    for(int i=0; i<n2; i++){
        if(d2->vertices[i]->x == vs->x && d2->vertices[i]->y == vs->y){ //iterate till first common point in second polygon
            pos1 = i;
            break;
        }
    }
    // cout<<pos1<<endl;
    Vertex *temp1 = d2->vertices[pos1];
    while(temp1->x != vt->x || temp1->y != vt->y){        //add all vertices in second polygon till the second common point
        v.push_back(temp1);
        // cout<<temp1->x<<"&"<<temp1->y<<"  ";
        temp1 = next(d2, temp1);
        // cout<<temp1->x<<"||"<<temp1->y<<"|"<<(temp1->x != vt->x && temp1->y != vt->y)<<"  ";
    }
    // cout<<v.size()<<" ";
    int pos2 = -1;
    for(int i=0; i<n1; i++){
        if(d1->vertices[i]->x == vt->x && d1->vertices[i]->y == vt->y){ //iterate till second point in first polygon
            pos2 = i;
            break;
        }
    }
    // cout<<pos2<<endl;
    Vertex *temp2 = d1->vertices[pos2];
    while(temp2->x != d1->vertices[0]->x || temp2->y != d1->vertices[0]->y){        //add all vertices in first polygon till the start
        v.push_back(temp2);
        // cout<<temp2->x<<"&"<<temp2->y<<"  ";
        temp2 = next(d1, temp2);
    }
    // cout<<v.size()<<" ";
    DCEL *dcel = make_dcel(v);
    // cout<<"merge sizes: "<<d1->vertices.size()<<" "<<d2->vertices.size()<<" "<<dcel->vertices.size()<<endl;
    return dcel;
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
    // cout<<"checking reflex\n";
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
            // cout<<"not in p "<<i<<endl;
            if(i==0){
                if(check_reflex(p[p.size()-1], p[i], p[i+1])) {
                    // cout<<"notch "<<p[i]->x<<" "<<p[i]->y<<endl;
                    notches.push_back(p[i]);
                }
            }
            else if(i==p.size()-1){
                if(check_reflex(p[i-1], p[i], p[0])) {
                    // cout<<"notch "<<p[i]->x<<" "<<p[i]->y<<endl;
                    notches.push_back(p[i]);
                }
            }
            else{
                if(check_reflex(p[i-1], p[i], p[i+1])) {
                    // cout<<"notch "<<p[i]->x<<" "<<p[i]->y<<endl;
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
                // cout<<"angle not reflex: "<<l[i]->x<<" "<<l[i]->y<<endl;
                flag = true;
                return false;
            }
        }
        else{
            if(!check_reflex(l[i]->inc_edge->next->org, l[i], v)){
                // cout<<"angle not reflex: "<<l[i]->x<<" "<<l[i]->y<<endl;
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
    // for(int i=0; i<vtr.size(); i++){
    //     // cout<<"vertex not in semiplane: "<<vtr[i]->x<<" "<<vtr[i]->y<<endl;
    // }
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

vector<pair<int,Edge*>> genLpv(Vertex* v,vector<DCEL*> dcel){
    vector<pair<int,Edge*>> lp;
    for(int i=1;i<dcel.size();i++){             //looping through DCELs
        for(auto e:dcel[i]->edges){         //looping through every edge of DCEL[i]
            if(e->org->x==v->x && e->org->y==v->y){
                bool flag = false;
                //check if edge is in main polygon - dcel[0]
                for(auto e0: dcel[0]->edges){
                    
                }
                //else add to lp
                if(!flag) lp.push_back(make_pair(i,e));
            }
        }
    }          

    return lp;
}

bool isConvex(DCEL *dcel, Vertex* v){
    //call check_reflex function
    // cout<<"checking convex for: "<<v->x<<" "<<v->y<<endl;
    Vertex *n = next(dcel, v), *p = prev(dcel, v);
    return !check_reflex(p, v, n);
}

void merging(vector<Edge *> diags,vector<DCEL*> polygons){
    vector<Edge *> ret;
    int n = polygons[0]->edges.size();
    for(int i=0; i<n; i++){
        ret.push_back(polygons[0]->edges[i]);
    }
    // cout<<diags.size();
    int m = diags.size();
    for(int i=0; i<m; i++){
        ret.push_back(diags[i]);
        // cout<<"hi";
    }
    int np = m+1;
    vector<bool> LDP(np+1,true);  //leaving index 0
    vector<int> LUP;
    
    Vertex* vs, *vt, *j1, *j2, *j3, *i1, *i2, *i3;

    for(int i=0;i<=np;i++){
        LUP.push_back(i);       //leaving index 0 
    }
    
    vector<pair<int,Edge*>> Lpvt;
    vector<pair<int,Edge*>> Lpvs;
    for(int j=1;j<=m;j++){      
        vs = diags[j-1]->org;
        vt = diags[j-1]->twin->org;
        
        Lpvt = genLpv(vt,polygons);
        Lpvs = genLpv(vs,polygons);
        // cout<<Lpvs.size()<<" "<<Lpvt.size()<<endl;
        // cout<<"edge: "<<vs->x<<" "<<vs->y<<" & "<<vt->x<<" "<<vt->y<<endl;
        if((Lpvs.size()>2 && Lpvt.size()>2)||(Lpvs.size()>2 && isConvex(polygons[0], vt))||(isConvex(polygons[0], vs) && Lpvt.size()>2)||(isConvex(polygons[0], vs) && isConvex(polygons[0], vt))){    
            j2 = vt;
            i2 = vs;
            // cout<<"bye "<<LUP[j]<<" "<<polygons.size();
            j3 = next(polygons[LUP[j]], vt); //diags[j]->next->org;
            // cout<<"next";
            i1 = prev(polygons[LUP[j]], vs); //diags[j]->prev->org;
            // cout<<" hi "<<endl;
            int u = -1;
            int i = 0;
            for(i=0;i<Lpvt.size();i++){
                if(Lpvt[i].second->twin->org->x==vs->x && Lpvt[i].second->twin->org->y==vs->y){
                    u = Lpvt[i].first;
                    break;
                }
            }
            if(u!=-1){
                // cout<<"yo "<<LUP[u]<<" "<<polygons[u]->vertices.size()<<" "<<(polygons[LUP[u]]==NULL)<<endl;
                i3 = next(polygons[LUP[u]], vs); 
                j1 = prev(polygons[LUP[u]], vt); 
            }
            else continue;
            // cout<<" lol"<<(i1==NULL)<<" "<<(i2==NULL)<<" "<<(i3==NULL)<<" "<<(j1==NULL)<<" "<<(j2==NULL)<<" "<<(j3==NULL)<<" "<<endl;
            if(!check_reflex(i1,i2,i3) && !check_reflex(j1,j2,j3)){
                //diag can be removed, new polygon possible
                np++;
                //create new polygon
                // cout<<" bye"<<endl;
                DCEL* newPoly = merge_dcel(polygons[LUP[j]], polygons[LUP[u]], vs, vt);
                //adding newPoly to list of polygons
                polygons.push_back(newPoly);
                LDP[j]=false;
                LDP[u]=false;
                LDP.push_back(true);
                
                LUP[j]=np;
                LUP[u]=np;
                
                for(int x=0; x<np; x++){
                    if(LUP[x]==j || LUP[x]==u) LUP[x]=np;
                }
                // cout<<"removing edge: "<<vs->x<<" "<<vs->y<<" & "<<vt->x<<" "<<vt->y<<endl;
                int k = find(ret.begin(), ret.end(), diags[j-1]) - ret.begin();
                ret.erase(ret.begin()+k);
            }  
        }   
    }
    // cout<<ret.size();
    generate_file(ret);
} 

Vertex *nextp(vector<Vertex*> d, Vertex *v){
    for(int i=0; i<d.size(); i++){
        if(d[i] == v){
            if(i==d.size()-1) return d[0];
            else return d[i+1];
        }
    }
    return NULL;
}

void decompose(DCEL *dcel){
    vector<DCEL *> dcels;
    dcels.push_back(dcel);
    vector<Vertex *> p;                             //p contains all the vertices in the polygon in a clockwise order
    vector<Edge *> ret, diags;  
    int v_size = dcel->vertices.size();                           //stores all the edges for plotting the polygon
    for(int i=0; i<dcel->vertices.size(); i++){
        p.push_back(dcel->vertices[i]);
        ret.push_back(dcel->edges[i]);
        // cout<<p[i]->x<<" "<<p[i]->y<<endl;
    }
    vector<vector<Vertex *>> l;                     //l contains the partitions of our convex polygons - from the 1st index
    l.push_back({p[0]});
    int m = 1;
    while(p.size() > 3){
        if(m>2*v_size+1){
            //cout<<"Not a simple polygon or vertices are not in clockwise order.";
            break;
        }
        Vertex *v1 = l[m-1][l[m-1].size()-1];
        Vertex *v2 = nextp(p, v1);
        l.push_back({});
        l[m].push_back(v1);
        l[m].push_back(v2);
        // cout<<"hi "<<v1->x<<" "<<v1->y<<endl<<"hi "<<v2->x<<" "<<v2->y<<endl;
        int t = 2;
        Vertex *v3 = nextp(p, v2);
        Vertex *vim = v1, *vi = v2, *vip = v3;
        //checking if the vertices satisfy all 3 conditions to be in the polygon
        while(!check_reflex(vim, vi, vip) && !check_reflex(vi, vip, v1) && !check_reflex(vip, v1, v2) && l[m].size()<p.size()){ 
            // cout<<"hi "<<vip->x<<" "<<vip->y<<endl;
            l[m].push_back(vip);
            t++;
            vim = vi;
            vi = vip;
            vip = nextp(p, vip);                                 //next(p, vip) to get the next vertex
            // cout<<l[m].size()<<endl;
        }
        // cout<<"hi1"<<endl;
        if(l[m].size() != p.size()){                            //check if l[m].size == p.size
            vector<Vertex *> lpvs = find_notches(p, l[m]);      //get the list of notches in the polygon
            // cout<<"no. of notches: "<<lpvs.size()<<endl;
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
                // cout<<"x: "<<minx<<" "<<maxx<<endl<<"y: "<<miny<<" "<<maxy<<endl;
                bool backward = false;
                while(!backward && lpvs.size()>0){
                    Vertex *v = lpvs[0];
                    while(v->x < minx || v->x > maxx || v->y < miny || v->y > maxy){
                        // cout<<"removing notch"<<endl;
                        lpvs.erase(lpvs.begin()); // check
                        if(lpvs.size() == 0) break;
                        v = lpvs[0];
                    }
                    // exit(0);
                    if(lpvs.size() > 0){
                        // cout<<"notches in polygon.."<<endl;
                        if(v_in_polygon(l[m], v)){
                            // cout<<" in notch"<<endl;
                            //exit(0);
                            vector<Vertex *> vtr = vs_in_semiplane(l[m], v);
                            for(int j=0; j<vtr.size(); j++){
                                int d = find(l[m].begin(), l[m].end(), vtr[j]) - l[m].begin();
                                // cout<<"removing vertex: "<<l[m][d]->x<<" "<<l[m][d]->y<<endl;
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
            // cout<<"joining edges: "<<l[m][0]->x<<" "<<l[m][0]->y<<" & "<<l[m][l[m].size()-1]->x<<" "<<l[m][l[m].size()-1]->y<<endl;
            // dcel->edges.push_back(e);
            DCEL *nw = make_dcel(l[m]);
            ret.push_back(nw->edges[nw->edges.size()-1]);
            diags.push_back(nw->edges[nw->edges.size()-1]);
            dcels.push_back(nw);
            //remove all vertices of l[m] from p except for the first and last one
            for(int j=1; j<l[m].size()-1; j++){ 
                int d = find(p.begin(), p.end(), l[m][j]) - p.begin();
                // cout<<" removing vertex: "<<l[m][j]->x<<" "<<l[m][j]->y<<endl;
                p.erase(p.begin()+d);
            }
            // n = n - l[m].size() + 2;
        }       
        m++;         
    }
    // cout<<"************";
    int w = diags.size();
    if(diags.size()>1){
        if(diags[w-1]->org->x==diags[w-2]->twin->org->x && diags[w-1]->org->y==diags[w-2]->twin->org->y)
            diags.pop_back();
    } 
    if(p.size()>2){
        dcels.push_back(make_dcel(p));
    }
    // cout<<dcels.size()<<endl;
    generate_file(ret);
    merging(diags, dcels);
}

// void vertex_dependency(DCEL *dcel){
//     vector<DCEL *> os;
//     int card = INT_MAX, s = 1, n = dcel->vertices.size();
//     vector<DCEL *> dcels;
//     while(s <= n){
//         dcels.push_back(new DCEL());
//         for(int i=s-1; i<n; i++){

//         }
//     }
// }

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
        // cout<<n<<" "<<endl;
        for(int i=0; i<n; i++){
            Vertex *v = new Vertex();
            Edge *e = new Edge(), *et = new Edge();
            double x, y;
            fp >> x;
            fp >> y;
            // cout<<x<<" "<<y<<endl;
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
    // for(int i=0; i<n; i++){
    //     cout<<dcel->vertices[i]->x<<"&"<<dcel->vertices[i]->y<<"&&"<<dcel->vertices[i]->inc_edge->next->org->x<<"&"<<dcel->vertices[i]->inc_edge->next->org->y<<"  ";
    // }
    generate_file(dcel->edges);
    decompose(dcel);
    // decompose1(dcel);
    return 0;
}
