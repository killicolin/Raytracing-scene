
struct Material{
    float Ka; // ambiant        coefficient
    float Kd; // diffuse        coefficient
    float Ks; // specular       coefficient
    float Kn; // specular power coefficient
};
  
struct ShaderColor{
    vec3 color; 
    float K; // reflection coefficient
};

    
struct LightInfo{
    vec3 pos; // light position
    vec3 col; // light color
    float power; // light power
};
    
//Structure use for csg construction with sphere 
struct sphereContact{
    float minDist; //distance between intersec1 and camera
    float maxDist; //distance between intersec2 and camera
    vec3 intersec1; // first intersection
    vec3 normal1; // normal associate to intersec1
    vec3 intersec2; // last intersection
    vec3 normal2; // normal associate to intersec2
};
//bound of light reflection
const int MaxBound=3;    
//sky parameter
const int skyID=0;
const vec3 skyColor= vec3(0,0,0); //black
//camera parameter
vec3 Camerapos=vec3(6,4,-5);
const vec3 cameraTarget=vec3(3,1,-8);
const float CameraFovY=80.0;
//sphere parameter
const int sphereID=1;
vec3 spherePos=cameraTarget+vec3(0,1,2);
const float sphereRadius=1.0;
const vec3 sphereCol=vec3(1.0,0.0,0.0);
const Material sphereMat=Material(0.2,/*Ka*/0.7,/*Kd*/1.0,/*Ks*/50.0/*Kn*/);
// Light Parameter
const vec3 ambiantCol= vec3(0.0,0.0,1.0); //blue
//light 1
const vec3 light1Col=vec3(1,1,1); //white
      vec3 light1Pos=vec3(8,10,-12);
const float light1Pow=0.8;
//light 2
const vec3 light2Col=vec3(1,1,0.5); //white/yellow
      vec3 light2Pos=vec3(3,10,1);
const float light2Pow=0.5;
//number of light
const int NB_LIGHTS = 2;
//light array
LightInfo lights[NB_LIGHTS];

// Plan
const int planID=2;
const vec3 planOrigin=vec3(0,0,0);
const vec3 dirX=vec3(1,0,0);
const vec3 dirZ=vec3(0,0,1);
const vec3 planColor1=vec3(1.0,1.0,1.0);
const vec3 planColor2=vec3(0.3,0.3,0.3);
const Material planMat=Material(0.2,/*Ka*/1.0,/*Kd*/0.2,/*Ks*/5.0/*Kn*/);
// pixel sampling use in the 
const int PIXEL_SAMPLING_GRID_SIZE=4;
//
const int PIXEL_SAMPLING_SIZE=1;
// CSG object parameters
//CSG sphere 1
const vec3 csgSphere1Pos = cameraTarget+vec3(-1.125,2,0);
const float csgSphere1Radius = 1.4;
const vec3 csgCol1 = vec3(1.0,0.5,0.0); // orange
//CSG sphere 2
const vec3 csgSphere2Pos = cameraTarget+vec3(1.125,2,0);
const float csgSphere2Radius = 1.4;
const vec3 csgCol2 = vec3(0.4,1.0,1.0); // cyan
//CSG sphere 3
const vec3 csgSphere3Pos = cameraTarget+vec3(0,2.5,0);
const float csgSphere3Radius = 0.5;
const vec3 csgCol3 = vec3(1.0,0.0,1.0); // magenta
//CSG sphere 4
const vec3 csgSphere4Pos = cameraTarget+vec3(0,2.4,-0.5);
const float csgSphere4Radius = 0.5;
const vec3 csgCol4 = vec3(5.0,5.0,0.0); // hyper yellow

const Material csgMat= Material(0.2,1.0,0.1,90.0);

const int csgId =3;

// id of the CSG sphere find
int subObjectIdC=0;

//init of the light array
void initLigths(){
    lights[1].pos=light1Pos;
    lights[0].pos=light2Pos;
    lights[1].col=light1Col;
    lights[0].col=light2Col;
    lights[1].power=light1Pow;
    lights[0].power=light2Pow;
}

//associate the color for the CSG construction
vec3 getCSGColorAtPoint(vec3 pt){
    if(subObjectIdC==1)
        return csgCol1; 
    if(subObjectIdC==2)
        return csgCol2; 
    if(subObjectIdC==3)
        return csgCol3; 
    if(subObjectIdC==4)
        return csgCol4; 
    return vec3(3);
}
// rayPos the ray origin
// rayDir the ray direction
void computeCameraRayFromPixel(in vec2 fragCoord,out vec3 rayPos,out vec3 rayDir){
    rayPos=Camerapos;
    // the projection
    vec3 P=vec3(fragCoord.x-(iResolution.x/2.),
                iResolution.y/2.-fragCoord.y,
                iResolution.y/(2.*tan(radians(CameraFovY/2.))));
    vec3 cameraDirection=normalize(cameraTarget-Camerapos);
    vec3 cy=vec3(0,-1,0);
    //cx the camera X-Axis convert in the 3d world space
    vec3 cx=normalize(cross(cy,cameraDirection));
    //cy the camera Y-Axis convert in the 3d world space
    cy=cross(cameraDirection,cx);
    rayDir=normalize(P.x*cx+P.y*cy+P.z*cameraDirection);
}

//this is use as a random function
vec2 noise2(vec2 location,vec2 delta){
    const vec2 c = vec2(12.9898,78.233);
    const float m = 43758.5453;
    return vec2(fract(sin(dot(location+delta   ,c))*m),
                fract(sin(dot(location+delta.yx,c))*m));
}

//mathematic function to know if the ray encounter the sphere
// rayDir the ray direction
// diff distance between rayPos and spherePos
// sphereRadius the sphere radius
void subRaySphere(vec3 rayDir,vec3 diff,float sphereRadius,
                  out float a,out float b ,out float c ,out float det){
    a=dot(rayDir,rayDir);
    b=2.*dot(diff,rayDir);
    c=dot(diff,diff)-pow(sphereRadius,2.);
    det=pow(b,2.)-4.*a*c;
}

// if the ray encouter a sphere return the distance, the intersection point and the normal associate
// rayPos the ray origin
// rayDir the ray direction
// spherePos the sphere position
// sphereRadius the sphere radius
// intersecS position of the first intersection
// normalS the normal associate to intersecS
float raySphere(vec3 rayPos,vec3 rayDir,vec3 spherePos,float sphereRadius,
                out vec3 normalS,out vec3 intersecS){
    float dist=-1.;
    vec3 diff=rayPos-spherePos;
    float a,b,c,det;
	subRaySphere(rayDir,diff,sphereRadius,a,b,c,det);
    //if the ray encounter only one point of the sphere
    if(det==0.){
        dist=(-b)/(2.*a);
    }
    //else if the ray encounter two point
    else if(det>0.){
        //take the first intersection distance
        dist=min((-b+sqrt(det))/(2.*a),(-b-sqrt(det))/(2.*a));
    }
    //the position of the first intersection point encountered along the ray.
    intersecS=rayPos+dist*rayDir;
    //the normal associate
    normalS=normalize(intersecS-spherePos);
    return dist;
}

// if the ray encouter a sphere return sphereContact associate
// rayPos the ray origin
// rayDir the ray direction
// spherePos the sphere position
// sphereRadius the sphere radius
sphereContact raySphereCSG(vec3 rayPos,vec3 rayDir,vec3 spherePos,float sphereRadius){
    vec3 diff=rayPos-spherePos;
    float a,b,c,det;
	subRaySphere(rayDir,diff,sphereRadius,a,b,c,det);
    sphereContact res;
    res.minDist=-1.0;
    res.maxDist=-1.0;
    //if the ray encounter only one point of the sphere
    if(det==0.){
        res.minDist=(-b)/(2.*a);
        res.maxDist=res.minDist;
    }
    //else if the ray encounter two point
    else if(det>0.){
        res.minDist=min((-b+sqrt(det))/(2.*a),(-b-sqrt(det))/(2.*a));
        res.maxDist=max((-b+sqrt(det))/(2.*a),(-b-sqrt(det))/(2.*a));
    }
    //the position of the first intersection point encountered along the ray.
    res.intersec1=rayPos+res.minDist*rayDir;
    res.normal1=normalize(res.intersec1-spherePos);
    //the position of the last intersection point.
    res.intersec2=rayPos+res.maxDist*rayDir;
    res.normal2=normalize(spherePos-res.intersec2);
    return res;
}

// if the ray encouter a plan return the distance, the intersection point and the normal associate
// rayPos the ray origin
// rayDir the ray direction
// planOrigin a point of the plan
// dirX and dirZ two vector in the plan, they have to be not colinear
// intersecP position of the intersection between the ray and the plan
// normalP the normal associate to intersecP
float rayPlan(vec3 rayPos,vec3 rayDir,vec3 planOrigin,vec3 dirX,vec3 dirZ,out vec3 normalP,out vec3 intersecP){
    float dist=-1.;
    normalP=cross(dirZ,dirX);
    float den = dot(normalP, rayDir);
    if (abs(den) <= 0.000001)   // To avoid numerical instabilities we consider the ray to be 
        return -1.0;
    dist=(dot(planOrigin,normalP)-dot(rayPos,normalP))/den;
    normalP=-sign(den) *normalP;
    intersecP=rayPos+dist*rayDir;
    return dist;
}

vec3 computePhongShading(vec3 color, LightInfo light,float shadowfactor,Material sphereMat,vec3 normal,vec3 L,vec3 R,vec3 V){
    vec3 ambiant= sphereMat.Ka * ambiantCol;
    vec3 diffuse= shadowfactor*color*sphereMat.Kd * light.col * max(.0,dot(normal,L));
    //if the pixel is in the shadow of an other object.
    float isInShadow=shadowfactor==sphereMat.Ka?0.:1.;
    vec3 specular=isInShadow*sphereMat.Ks *  light.col * pow(max(0.0,dot(R,V)),sphereMat.Kn);
    return (ambiant+diffuse+specular);
}

//define the grid color 
vec3 definePlanColor(vec3 intersecPlan){
    float u=dot(intersecPlan,dirX);
    float v=dot(intersecPlan,dirZ);
    return mod(floor(u*0.5)+floor(v*0.5),2.0)<1.0 ? planColor1 : planColor2;
}

// soustract c2 to c1
//c1 a sphereContact
//c2 a sphereContact to soustract to c1
//i1 the color ID of c1
//i2 the color ID of c2
//subObjectId is the result of c1 - c2
sphereContact notCSG(sphereContact c1,sphereContact c2,int i1,int i2,out int subObjectId){
    sphereContact res;
    res.minDist=-1.0;
    //if the ray never encounter c1 and c2
    if(c1.minDist==-1.0 && c2.minDist==-1.0){
        res.minDist=-1.0;
    }
    //if the ray encounter c1 only
    else if(c2.minDist==-1.0){
        res=c1;
        subObjectId=i1;
    }
    //if the ray encounter c1 and c2
    else if(c1.minDist>=0.0 && c2.minDist>=0.0){
        //if the sphere are separated or c1 is before c2
        if(c1.minDist<c2.minDist || c1.minDist>=c2.maxDist){
            res=c1;
            subObjectId=i1;
        }
        //if the sphere are join, and c2 ending in the middel of c1
        else if(c2.maxDist<=c1.maxDist){
            res.minDist=c2.maxDist;
            res.maxDist=c1.maxDist;
            res.intersec1=c2.intersec2;
            res.normal1=c2.normal2;
            res.intersec2=c1.intersec2;
            res.normal2=c1.normal2;
            subObjectId=i2; 
        }
    }
    return res;
}

// intersection between c1 and c2
//c1 a sphereContact
//c2 a sphereContact
//i1 the color ID of c1
//i2 the color ID of c2
//subObjectId is the result of c1 & c2
sphereContact andCSG(sphereContact c1,sphereContact c2,int i1,int i2,out int subObjectId){
    sphereContact res;
    res.minDist=-1.0;
    //if the sphere are join
    if(c1.maxDist>=c2.minDist&&c2.maxDist>=c1.minDist){
        res.minDist=max(c1.minDist,c2.minDist);
        res.maxDist=min(c1.maxDist,c2.maxDist);
        //if c1 is before
        if(res.minDist==c1.minDist){
            res.intersec1=c1.intersec1;
            res.normal1=c1.normal1;
            res.intersec2=c1.intersec2;
            res.normal2=c1.normal2;
            subObjectId=i1;
        }
        //if c2 is before
        else{
            res.intersec1=c2.intersec1;
            res.normal1=c2.normal1;
            res.intersec2=c2.intersec2;
            res.normal2=c2.normal2;
            subObjectId=i2;
        }
    }
    return res;
}

// Union between c1 and c2
//c1 a sphereContact
//c2 a sphereContact
//i1 the color ID of c1
//i2 the color ID of c2
//subObjectId is the result of c1 | c2
sphereContact orCSG(sphereContact c1,sphereContact c2,int i1,int i2,out int subObjectId){
    sphereContact res;
    subObjectId=0;
    if(c1.minDist==-1.0 || c2.minDist==-1.0){
        //if the ray never encounter c1 and c2
        if(c1.minDist==-1.0 && c2.minDist==-1.0){
            res.minDist=-1.0;
        }
        //if the ray encounter c2 only
        else if(c1.minDist==-1.0){
            res=c2;
            subObjectId=i2;
        }
        //if the ray encounter c1 only
        else{
            res=c1;
            subObjectId=i1;
        }
    }
    //if the sphere are join
    else{
        res.minDist=min(c1.minDist,c2.minDist);
        res.maxDist=max(c1.maxDist,c2.maxDist);
        //if c1 is before
        if(res.minDist==c1.minDist){
            res.intersec1=c1.intersec1;
            res.normal1=c1.normal1;
            res.intersec2=c1.intersec2;
            res.normal2=c1.normal2;
            subObjectId=i1;
        }
        //if c2 is before
        else{
            res.intersec1=c2.intersec1;
            res.normal1=c2.normal1;
            res.intersec2=c2.intersec2;
            res.normal2=c2.normal2;
            subObjectId=i2;
        }
    }
    return res;
}

float rayCSG(vec3 rayPos,vec3 rayDir, out vec3 intersecPt, out vec3 normal, out int subObjectId){
    sphereContact contact1 = raySphereCSG(rayPos,rayDir,csgSphere1Pos,csgSphere1Radius);
    sphereContact contact2 = raySphereCSG(rayPos,rayDir,csgSphere2Pos,csgSphere2Radius);
    sphereContact contact3 = raySphereCSG(rayPos,rayDir,csgSphere3Pos,csgSphere3Radius);
    sphereContact contact4 = raySphereCSG(rayPos,rayDir,csgSphere4Pos,csgSphere4Radius);
    sphereContact res= andCSG(contact1,contact2,1,2,subObjectId);
    res= orCSG(contact3,res,3,subObjectId,subObjectId);
    res= notCSG(res,contact4,subObjectId,4,subObjectId);
    intersecPt=res.intersec1;
    normal=res.normal1;
    return res.minDist;
}

int definefragObject(float dist1,float dist2,float dist3){
    if(((dist1<= dist2 || dist2<=0.) && (dist1<= dist3 || dist3<=0.)) && dist1>0.){
        return sphereID;
    }
    else if (((dist2<= dist1 || dist1<=0.) && (dist2<= dist3 || dist3<=0.)) && dist2>0.){
        return planID;
    }
    else if (((dist3<= dist1 || dist1<=0.) && (dist3<= dist2 || dist2<=0.)) && dist3>0.){
        return csgId;
    }
    return skyID;
}

float getShadowFactorPoint(vec3 I,vec3 normalPoint,Material objectMat,vec3 L, float Ldist){
    vec3 I2=I+normalPoint*0.001;
    vec3 normal,P,intersecCSG,normalCSG;
    vec3 Ldir=L-I2;
    Ldist=distance(L,I2);
    int tmpTrash;
    float dist1=raySphere(I2,Ldir,spherePos,sphereRadius,normal,P);
    float dist2=rayPlan(I2,Ldir,planOrigin,dirX,dirZ,normal,P);
    float dist3=rayCSG(I2,Ldir,intersecCSG,normalCSG,tmpTrash);
    int obj=definefragObject(dist1,dist2,dist3);
    if(obj==sphereID && dist1<Ldist || obj==planID && dist2<Ldist || obj==csgId && dist3<Ldist)
        return objectMat.Ka;
    return 1.;
    
}

// animate camera, light and shpere
void animateScene(float time){
    const float pi=3.1418;
    const float rs=2.0;
    const float npr=5.0;
    float mn=2.0*pi*time/npr;
    spherePos=cameraTarget+rs+vec3(-sin(mn),-1.9,cos(mn));
    //lights[0].pos=vec3(4.*cos(2.*time),11.5+9.5*cos(time),4.*sin(time));
    float targetDist=length(cameraTarget-Camerapos);
    //Camerapos+=vec3(1.,0.,targetDist);
    Camerapos-=targetDist*vec3(sin(time),min(sin(time*0.3),-0.3),cos(time));
}

float computeNearestIntersection(vec3 rayPos, vec3 rayDir,out int objectId, out vec3 intersecI, out vec3 normalI){
    vec3 normalS,normalP,normalCSG,P,intersecPlan,intersecCSG,I;
    float distRes;
    float dist1=raySphere(rayPos,rayDir,spherePos,sphereRadius,normalS,P);
    float dist2=rayPlan(rayPos,rayDir,planOrigin,dirX,dirZ,normalP,intersecPlan);
    float dist3=rayCSG(rayPos,rayDir,intersecCSG,normalCSG,subObjectIdC);
    objectId=definefragObject(dist1,dist2,dist3);
    switch(objectId){
        case(sphereID):
            intersecI=P;
            normalI=normalS;
            distRes=dist1;
        	break;
        case(planID):
            intersecI=intersecPlan;
            normalI=normalP;
            distRes=dist2;
        	break;
        case(csgId):
            intersecI=intersecCSG;
            normalI=normalCSG;
            distRes=dist3;
        	break;
        default:
            distRes=-1.;
            break;
    }
    return distRes;
}

//define the color of the point encounter
//objectId the id of the object find
//pt the position of the point
//objectMat the material of the object find
//return the color of the object encounter
vec3 getObjectColorAtPoint(int objectId, vec3 pt, out Material objectMat)
{
    switch(objectId){
        case(sphereID):
            objectMat = sphereMat;
            return sphereCol;
        case(planID):
            objectMat = planMat;
            return definePlanColor(pt);
        	break;
        case(csgId):
            objectMat = csgMat;
            return getCSGColorAtPoint(pt);
        default:
            return skyColor;
    }
}

vec3 RaycastAtPixelCoord(vec2 pixCoord){
    ShaderColor tab[MaxBound];
    vec3 rayPos,rayDir,normal,intersec;
    float shadowfactor;
    int obj;
    Material objectMat;
    LightInfo tmpLight;
    computeCameraRayFromPixel(pixCoord,rayPos,rayDir);
    for(int i=0; i<MaxBound;++i){
        float dist=computeNearestIntersection(rayPos,rayDir,obj,intersec,normal);
        if(dist>0.0){
            vec3 tmpcolor=vec3(0,0,0);
            for(int j = 0; j <NB_LIGHTS ;j++){
                tmpLight=lights[j];
                vec3 L= normalize(tmpLight.pos-intersec);
                vec3 R= 2.*normal*dot(L,normal)-L;
                vec3 V= normalize(Camerapos-intersec);
                vec3 objColor=getObjectColorAtPoint(obj,intersec,objectMat);
                shadowfactor=getShadowFactorPoint(intersec,normal,objectMat,tmpLight.pos,dist);
                tmpcolor+=tmpLight.power*computePhongShading(objColor,tmpLight,shadowfactor,objectMat,normal,L,R,V);
            }
            tab[i].color=tmpcolor;
            tab[i].K=objectMat.Ks;
            rayPos=intersec+normal*0.001;
            rayDir=reflect(rayDir,normal);
        }
        else{
            tab[i].color=skyColor.xyz;
            tab[i].K=0.0;
            break;
        }
    }
    vec3 tmp=vec3(.0,.0,.0);
    float Ktmp=1.0;
    for(int i=0; i<MaxBound;++i){
        tmp+=Ktmp *tab[i].color;
        Ktmp= tab[i].K;
        if(Ktmp==0.)
            break;
    }
    return tmp;
}
// Antialiasing function
vec3 regularGridAntialiasing(in vec2 fragCoord){
    vec2 pixCoord;
    vec3 res=vec3(0,0,0);
    const float pas=1.0/float(PIXEL_SAMPLING_GRID_SIZE);
    //for one pixel, send a grid of ray
    for(float i=-0.5; i<0.5;i+=pas){ 
        for(float j=-0.5; j<0.5;j+=pas){ 
            pixCoord=vec2(fragCoord.x+i,fragCoord.y+j);
            res+=RaycastAtPixelCoord(pixCoord);
        }
    }
    //normalise the result
    return res/float(PIXEL_SAMPLING_GRID_SIZE*PIXEL_SAMPLING_GRID_SIZE);
}

/*
vec3 stochasticAntialiasing(in vec2 fragCoord){
    vec2 pixCoord;
    vec3 res=vec3(0,0,0);
    for(int i=0; i<PIXEL_SAMPLING_SIZE;i++){ 
        pixCoord=fragCoord+(noise2(fragCoord,vec2(iTime,iTime)));
        res+=RaycastAtPixelCoord(pixCoord);
    }
    return res/float(PIXEL_SAMPLING_SIZE);
}*/

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    initLigths();
    animateScene(iTime);
    fragColor=vec4(regularGridAntialiasing(fragCoord),1.0);
    //fragColor=vec4(stochasticAntialiasing(fragCoord),1.0);
}
