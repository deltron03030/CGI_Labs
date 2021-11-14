/*
 * SDL2 skeleton for lab assignments 1–3 of the KTH course DH2323,
 * Computer Graphics and Interaction (and also SU DA3001, Datalogi II).
 *
 * See README.md for details.
 */

#include <iostream>
#include "glm/glm.hpp"
#include <vector>
#include "SDL2auxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::vec2;
using glm::mat3;

// --------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL2Aux *sdlAux;
int focalLength = SCREEN_HEIGHT;

int t;
vector<Triangle> triangles;
vec3 cameraPos( 0, 0, -3.001 );
vec3 currentColor;

mat3 Ry;
mat3 Rx;
mat3 R;
float yaw = 0; // Yaw angle controlling camera rotation around y-axis
float pitch = 0;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

vec3 lightPos(0,-0.5,-0.7);
vec3 lightPower = 14.1f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

vec3 currentNormal;
vec3 currentReflectance;


struct Pixel
{
  int x;
  int y;
  float zinv;
  vec3 illumination;
  vec3 pos3d;
  vec3 pos3dZ;
};

struct Vertex
{
  vec3 position;
  vec3 normal;
  vec3 reflectance;
};


// --------------------------------------------------------
// FUNCTION DECLARATIONS

void Draw();
void Update();
  void VertexShader1( const vec3& v, ivec2& p );
  void Interpolate1( ivec2 a, ivec2 b, vector<ivec2>& result );
  void DrawLineSDL1( SDL2Aux* screen, ivec2 a, ivec2 b, vec3 color );
  void DrawPolygonEdges( const vector<vec3>& vertices );
  void ComputePolygonRows1(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels);
  void DrawPolygonRows1( const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels );
  void DrawPolygon1( const vector<vec3>& vertices );

void DrawPolygon( const vector<Vertex>& vertices );
void DrawLineSDL(SDL2Aux* screen, Pixel a, Pixel b, vec3 color);
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void VertexShader( const Vertex& v, Pixel& p );
void PixelShader( const Pixel& p );

vec3 ComputeNormal(vec3& v0, vec3& v1, vec3& v2);
void RotateX();
void RotateY();

// --------------------------------------------------------
// FUNCTION DEFINITIONS

int main(int argc, char* argv[]) {
  LoadTestModel( triangles );
  sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);
  t = SDL_GetTicks();

  while (!sdlAux->quitEvent()) {
    Update();
    Draw();
  }

  sdlAux->saveBMP("screenshot.bmp");
  return 0;
}


void Draw() {
  for( int y=0; y<SCREEN_HEIGHT; ++y )
  {
    for( int x=0; x<SCREEN_WIDTH; ++x )
      depthBuffer[y][x] = 0;
  }

  sdlAux->clearPixels();

  for( int i=0; i<triangles.size(); ++i )
  {
    currentColor = triangles[i].color;
    currentNormal = triangles[i].normal;
    currentReflectance = vec3(1,1,1);
    vector<Vertex> vertices(3);
    vertices[0].position = triangles[i].v0;
    vertices[1].position = triangles[i].v1;
    vertices[2].position = triangles[i].v2;
    //for (int j=0; j<3; j++){
      //vertices[j].normal = triangles[i].normal;
      //vertices[j].reflectance = vec3(1,1,1);
    //}
    DrawPolygon(vertices);
 /*
    currentColor = triangles[i].color;
    vector<Vertex> vertices(3);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;
    DrawPolygon( vertices ); */
  }
  sdlAux->render();
}

vec3 ComputeNormal(vec3&v0, vec3&v1, vec3&v2)
{
    glm::vec3 e1 = v1-v0;
    glm::vec3 e2 = v2-v0;
    return glm::normalize( glm::cross( e2, e1 ) );
}

void Update()
{
  // Compute frame time:
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  cout << "Render time: " << dt << " ms." << endl;
  float amount = 0.005f;

  const Uint8* keystate = SDL_GetKeyboardState(0);

  if( keystate[SDL_SCANCODE_UP] )
    cameraPos.z += amount;

  if( keystate[SDL_SCANCODE_DOWN] )
    cameraPos.z -= amount;

  if( keystate[SDL_SCANCODE_LEFT] )
    yaw += amount;

  if( keystate[SDL_SCANCODE_RIGHT] )
    yaw -= amount;

  if( keystate[SDL_SCANCODE_L] )
    pitch -= amount;

  if( keystate[SDL_SCANCODE_PERIOD] )
    pitch += amount;

  if( keystate[SDL_SCANCODE_W] )
    //pitch -=amount/2;
    lightPos.z += 0.02;

  if( keystate[SDL_SCANCODE_S] )
    //pitch +=amount/2;
    lightPos.z -= 0.02;

  if( keystate[SDL_SCANCODE_D] )
    lightPos.x += 0.02;

  if( keystate[SDL_SCANCODE_A] )
    lightPos.x -= 0.02;

  if( keystate[SDL_SCANCODE_E] )
    lightPos.y += 0.02;

  if( keystate[SDL_SCANCODE_Q] )
    lightPos.y -= 0.02;
  
  Rx = mat3(1, 0, 0, 0, cos(pitch), -sin(pitch), 0, sin(pitch), cos(pitch));
  Ry = mat3(cos(yaw), 0, sin(yaw), 0, 1, 0, -sin(yaw), 0, cos(yaw));
}

void VertexShader( const Vertex& v, Pixel& p )
{
  vec3  newP = (v.position - cameraPos) * Rx * Ry;
  p.zinv = float (1.0 / newP.z);
  p.x = focalLength * newP.x / newP.z + SCREEN_WIDTH/2;
  p.y = focalLength * newP.y / newP.z + SCREEN_HEIGHT/2;
  /*
  vec3 r = glm::normalize(lightPos - v.position);
  float dist = glm::distance(lightPos, v.position);
  float dot = glm::dot(r, v.normal);
  vec3 D = lightPower * fmaxf(dot, 0.0f) / float (4 * M_PI * dist * dist);
  p.illumination = v.reflectance * (D + indirectLightPowerPerArea);
  cout<<p.illumination.x<<" "<<p.illumination.y<<" "<<p.illumination.z<<endl; 
  */
  p.pos3d = v.position;
  p.pos3dZ = v.position / v.position.z;
}

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result )
{
  int N = result.size();
  vec3 step;
  step.x = float(b.x - a.x) / float(glm::max(N - 1, 1));
  step.y = float(b.y - a.y) / float(glm::max(N - 1, 1));
  step.z = float(b.zinv - a.zinv) / float(glm::max(N - 1, 1));
  vec3 lightStep = (b.illumination - a.illumination) / float(glm::max(N - 1, 1));
  vec3 pos3dStep = (b.pos3d - a.pos3d ) / float(glm::max(N - 1, 1));
  vec3 pos3dZstep = (b.pos3dZ - a.pos3dZ ) / float(glm::max(N - 1, 1));
  //vec3 pos3dStep = (b.pos3d / b.pos3d.z - a.pos3d / a.pos3d.z) / float(glm::max(N - 1, 1));
  Pixel current(a);
  for (int i=0; i<N; i++)
  {
    current.x = a.x + i*step.x;
    current.y = a.y + i*step.y;
    current.zinv = a.zinv + i*step.z;
    current.illumination = a.illumination + float(i) * lightStep;
    current.pos3d = a.pos3d + float(i) * pos3dStep ;
    current.pos3dZ = a.pos3dZ + float(i) * pos3dZstep ;
    //current.pos3d = a.pos3d / float(a.pos3d.z) + float(i) * pos3dStep ;
    //current.pos3d = current.pos3d * current.pos3d.z;
    result[i] = current;
    // current.pos3d += pos3dStep;
  }
}


void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels)
{
  // 1. Find max and min y-value of the polygon and compute the number of rows it occupies.
  int ymax = -numeric_limits<int>::max();
  int ymin = +numeric_limits<int>::max();
  for (int i=0; i<vertexPixels.size(); i++){
    if (vertexPixels[i].y > ymax)
      ymax = vertexPixels[i].y;
    if (vertexPixels[i].y < ymin)
      ymin = vertexPixels[i].y;
  }
  int ROWS = ymax - ymin + 1;
  // 2. Resize leftPixels and rightPixels so that they have an element for each row
  leftPixels.resize(ROWS);
  rightPixels.resize(ROWS);
  // 3. Initialize the x-coordinates in leftPixels to some really large value 
  // and the x-coordinates in rightPixels to some really small value.
  for( int j=0; j<ROWS; ++j )
  {
    leftPixels[j].x = +numeric_limits<int>::max();
    leftPixels[j].y = j + ymin;
    rightPixels[j].x = -numeric_limits<int>::max();
    rightPixels[j].y = j + ymin;
  }
  // 4. Loop through all edges of the polygon and use linear interpolation to find the x-coordinate for
  // each row it occupies. Update the corresponding values in rightPixels and leftPixels.
  for (int i=0; i<vertexPixels.size(); i++)
  { 
    int next = (i+1) % vertexPixels.size(); 
    // int xD = abs(vertexPixels[i].x - vertexPixels[next].x);   // difference of x-coordinates between two vertices.        
    int yD = abs(vertexPixels[i].y - vertexPixels[next].y);
    // int pixels = glm::max( xD, yD ) + 1;       // find the number of pixels of each polygon edge.
    int pixels = yD + 1;
    vector<Pixel> line(pixels);
    Interpolate(vertexPixels[i], vertexPixels[next], line);   // interpolate each edge.
    for (int j=0; j<line.size(); j++)          // loop through pixels on each edge
    {         
    /*
      for (int c = 0; c < ROWS; c++)              // loop though each row of the polygon
      {
        if (leftPixels[c].y == line[j].y)
        {
          if (leftPixels[c].x>line[j].x)
          {
            leftPixels[c].x = line[j].x;
            leftPixels[c].zinv = line[j].zinv;
            leftPixels[c].pos3d = line[j].pos3d;
            leftPixels[c].pos3dZ = line[j].pos3dZ;
          }
          if (rightPixels[c].x<line[j].x)
          {
            rightPixels[c].x = line[j].x;
            rightPixels[c].zinv = line[j].zinv;
            rightPixels[c].pos3d = line[j].pos3d;
            rightPixels[c].pos3dZ = line[j].pos3dZ;
          }
          break;
        }
      }
      */
      for (Pixel p : line)
      {
        int index = p.y - ymin;                  // the corresponding index on leftPixels/rightPixels arrays.
        if (p.x < leftPixels[index].x)
        {
          leftPixels[index].x = p.x;
          leftPixels[index].zinv = p.zinv;
          leftPixels[index].pos3d = p.pos3d;
          leftPixels[index].pos3dZ = p.pos3dZ;
          //leftPixels[i].step = p.step;
          leftPixels[index].illumination = p.illumination;
        }

        if (p.x > rightPixels[index].x)
        {
          rightPixels[index].x = p.x;
          rightPixels[index].zinv = p.zinv;
          rightPixels[index].pos3d = p.pos3d;
          rightPixels[index].pos3dZ = p.pos3dZ;
          //rightPixels[i].step = p.step;
          rightPixels[index].illumination = p.illumination;
        } 
      }  
    } 
  } 
}

void DrawLineSDL(SDL2Aux* screen, Pixel a, Pixel b, vec3 color)
{
  Pixel delta;
  delta.x = glm::abs(a.x - b.x);
  delta.y = glm::abs(a.y - b.y);

  int pixels = glm::max(delta.x, delta.y) + 1;

  vector<Pixel> line(pixels);
  Interpolate(a, b, line);
  for (Pixel p : line)
    PixelShader(p);
}

void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels)
{
  for (int i=0; i<leftPixels.size(); i++)
  { 
    DrawLineSDL(sdlAux, leftPixels[i], rightPixels[i], currentColor);
    /* int num = rightPixels[i].x - leftPixels[i].x + 1;    // number of pixels on each line.
    vector<Pixel> line(num);                             // vector array of all the pixels on the line.
    Interpolate(leftPixels[i], rightPixels[i], line);
    for (Pixel p : line){
      PixelShader(p);
    } */
  }  
}

void DrawPolygon( const vector<Vertex>& vertices )
{
  int V = vertices.size();
  // vector<ivec2> vertexPixels( V );
  vector<Pixel> vertexPixels( V );
  for( int i=0; i<V; ++i )
     VertexShader( vertices[i], vertexPixels[i] );    // vertexPixels stores 2D coordinates of the vertices
  // vector<ivec2> leftPixels;
  // vector<ivec2> rightPixels;
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
  // cout << "LeftPixels: (" << leftPixels[0].x << " " << leftPixels[0].y << ")" << endl;
  DrawPolygonRows( leftPixels, rightPixels );
}

void PixelShader( const Pixel& p )
{ 
  int x = p.x;
  int y = p.y;
  if( p.zinv > depthBuffer[y][x] )
  {
    depthBuffer[y][x] = p.zinv;
    
    vec3 r = glm::normalize(lightPos - p.pos3d * p.pos3dZ.z);
    //vec3 r = glm::normalize(lightPos - p.pos3d);
    float dist = glm::distance(lightPos, p.pos3d * p.pos3dZ.z);
    //float dist = glm::distance(lightPos, p.pos3d);
    float dot = glm::dot(r, currentNormal);
    vec3 D = lightPower * fmaxf(dot, 0.0f) / float (4 * M_PI * dist * dist);
    vec3 illumination = currentReflectance * (D + indirectLightPowerPerArea); 
    //cout<<illumination.x<<" "<<illumination.y<<" "<<illumination.z<<endl; 
    sdlAux->putPixel(x, y, illumination * currentColor );
    // sdlAux->putPixel(x, y, p.illumination * currentColor);
  }
}
// ---------------------------------functions for first part of lab---------------------------------

void Interpolate1( ivec2 a, ivec2 b, vector<ivec2>& result )
{ 
  int N = result.size();
  vec2 step = vec2(b-a) / float(max(N-1,1));
  vec2 current(a);
  for (int i=0; i<N; i++)
  {
    result[i] = current;
    current += step;
  }
}


void VertexShader1( const vec3& v, ivec2& p )
{
  vec3  newP = (v - cameraPos) * Rx * Ry;
  p.x = focalLength * newP.x / newP.z + SCREEN_WIDTH/2;
  p.y = focalLength * newP.y / newP.z + SCREEN_HEIGHT/2;
}


void DrawPolygonRows1( const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels )
{
  for (int i=0; i<leftPixels.size(); i++)
    DrawLineSDL1(sdlAux, leftPixels[i], rightPixels[i], currentColor);
}

void DrawPolygon1( const vector<vec3>& vertices )
{
  int V = vertices.size();
  vector<ivec2> vertexPixels( V );
  //vector<Pixel> vertexPixels( V );
  for( int i=0; i<V; ++i )
     VertexShader1( vertices[i], vertexPixels[i] );    // vertexPixels stores 2D coordinates of the vertices
  vector<ivec2> leftPixels;
  vector<ivec2> rightPixels;
  //vector<Pixel> leftPixels;
  //vector<Pixel> rightPixels;
  ComputePolygonRows1( vertexPixels, leftPixels, rightPixels );
  // cout << "LeftPixels: (" << leftPixels[0].x << " " << leftPixels[0].y << ")" << endl;
  DrawPolygonRows1( leftPixels, rightPixels );
}

void DrawLineSDL1( SDL2Aux* screen, ivec2 a, ivec2 b, vec3 color )
{
  ivec2 delta = glm::abs( a - b );
  int pixels = glm::max( delta.x, delta.y ) + 1;
  vector<ivec2> line( pixels );
  Interpolate1( a, b, line );
  for (int i=0; i<line.size(); i++)
  {
    screen -> putPixel(line[i].x, line[i].y, color);
  }
}


void DrawPolygonEdges( const vector<vec3>& vertices )
{
  int V = vertices.size();
  // Transform each vertex from 3D world position to 2D image position: 
  vector<ivec2> projectedVertices( V );
  for( int i=0; i<V; ++i )
  {
    VertexShader1( vertices[i], projectedVertices[i] );
  }
  // Loop over all vertices and draw the edge from it to the next vertex:
  for( int i=0; i<V; ++i )
  {
    int j = (i+1)%V;         // The next vertex
    vec3 color( 1, 1, 1 );
    DrawLineSDL1( sdlAux, projectedVertices[i], projectedVertices[j], color);
  }
} 


void ComputePolygonRows1(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels )
{
  // 1. Find max and min y-value of the polygon and compute the number of rows it occupies.
  int ymax = -numeric_limits<int>::max();
  int ymin = +numeric_limits<int>::max();
  for (int i=0; i<vertexPixels.size(); i++){
    if (vertexPixels[i].y > ymax){
      ymax = vertexPixels[i].y;
    }
    if (vertexPixels[i].y < ymin){
      ymin = vertexPixels[i].y;
    }
  }
  int ROWS = ymax - ymin + 1;
  // 2. Resize leftPixels and rightPixels so that they have an element for each row
  leftPixels.resize(ROWS);
  rightPixels.resize(ROWS);
  // 3. Initialize the x-coordinates in leftPixels to some really large value 
  // and the x-coordinates in rightPixels to some really small value.
  for( int j=0; j<ROWS; ++j )
  {
    leftPixels[j].x = +numeric_limits<int>::max();
    leftPixels[j].y = j + ymin;
    rightPixels[j].x = -numeric_limits<int>::max();
    rightPixels[j].y = j + ymin;
  }
  // 4. Loop through all edges of the polygon and use linear interpolation to find the x-coordinate for
  // each row it occupies. Update the corresponding values in rightPixels and leftPixels.
  for (int i=0; i<vertexPixels.size(); i++)
  {
    ivec2 delta = glm::abs( vertexPixels[i] - vertexPixels[(i+1)%vertexPixels.size()] );
    int pixels = glm::max( delta.x, delta.y ) + 1;    // find the number of pixels of each polygon edge.
    vector<ivec2> line(pixels);
    Interpolate1(vertexPixels[i], vertexPixels[(i+1)%vertexPixels.size()], line);   // interpolate each edge.
    for (ivec2 p : line)          // loop through pixels on each edge
    {         
      /*int index = p.y - ymin;               // the corresponding index on leftPixels/rightPixels arrays.
      if ( p.x < leftPixels[index].x )
      {
        leftPixels[index].x = p.x; 
      }
      if ( p.x > rightPixels[index].x )
      {
        rightPixels[index].x = p.x;
      }  */

      for (int c = 0; c < ROWS; c++)
      {
        if (leftPixels[c].y == p.y)
        {
          if (p.x < leftPixels[c].x)
            leftPixels[c].x = p.x;
          if (p.x > rightPixels[c].x)
            rightPixels[c].x = p.x;
          break;
        }
      } 
    }
  }
}

// --------------------------------------------------------------------------
