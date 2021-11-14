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
#include <numeric>
#include <math.h>

using namespace std;
using glm::vec3;
using glm::mat3;

// --------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL2Aux *sdlAux;
int t;
mat3 R;
float yaw = 0;

vector<Triangle> triangles;

float focalLength = SCREEN_HEIGHT;
vec3 cameraPos( 0,0, -(2*focalLength / SCREEN_HEIGHT+1) );
//vec3 cameraPos( 0,0, -3);

vec3 lightPos( 0, -0.5, -0.7 );
vec3 lightColor = 14.f * vec3( 1, 1, 1 );
vec3 indirectLight = 0.5f*vec3( 1, 1, 1 );



struct Intersection
{
      vec3 position;
      float distance;
      int triangleIndex;
};
// --------------------------------------------------------
// FUNCTION DECLARATIONS

void Draw();
void Update();
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection
);
vec3 DirectLight( const Intersection& i );

// --------------------------------------------------------
// FUNCTION DEFINITIONS

int main(int argc, char* argv[]) {
  sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);
  t = SDL_GetTicks();	// Set start value for timer.
  LoadTestModel( triangles );

  while (!sdlAux->quitEvent()) {
    Update();
    Draw();
  }
  
  sdlAux->saveBMP("screenshot.bmp");
  return 0;
}

void Update(){
    // Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;
  float amount = 0.1;
  float amount2 = M_PI/180;

  const Uint8* keystate = SDL_GetKeyboardState( 0 );
      if( keystate[SDL_SCANCODE_UP] )
      {
            // Move camera forward
            cameraPos.z += amount;
      }
      if( keystate[SDL_SCANCODE_DOWN] )
      {
            // Move camera backward
            cameraPos.z -= amount;
      }
      if( keystate[SDL_SCANCODE_LEFT] )
      {
            // Move camera to the left
            yaw = -amount2;
            R = mat3(cos(yaw), 0, sin(yaw), 0, 1, 0, -sin(yaw), 0, cos(yaw));
            cameraPos = R * cameraPos;
      }
      if( keystate[SDL_SCANCODE_RIGHT] )
      {
            // Move camera to the right
            yaw = amount2;
            R = mat3(cos(yaw), 0, sin(yaw), 0, 1, 0, -sin(yaw), 0, cos(yaw));
            cameraPos = R * cameraPos;
      }

      if( keystate[SDL_SCANCODE_W] )
      lightPos.z += amount;         // move light position forward.
      if( keystate[SDL_SCANCODE_S] )
      lightPos.z -= amount;
      if( keystate[SDL_SCANCODE_A] )
      lightPos.x -= amount;
      if( keystate[SDL_SCANCODE_D] )
      lightPos.x += amount;
}

void Draw() {
  sdlAux->clearPixels();
   
  vec3 dir;
  Intersection closestIntersection;

  for (int y = 0; y < SCREEN_HEIGHT; ++y) {
    for (int x = 0; x < SCREEN_WIDTH; ++x) { 

      dir = glm::normalize (vec3 (x-SCREEN_WIDTH/2, y-SCREEN_HEIGHT/2, focalLength));
      bool intersect = ClosestIntersection(cameraPos, dir, triangles, closestIntersection);

      if (intersect){
        vec3 triangleColor = triangles[closestIntersection.triangleIndex].color;
        vec3 directColor = DirectLight(closestIntersection);
        vec3 color = triangleColor * (directColor + indirectLight);
        sdlAux->putPixel(x, y, color);
        }
      else
        sdlAux->putPixel(x, y, vec3 (0,0,0));
    }
  }
  
  sdlAux->render();
}

bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& closestIntersection)
{
  vec3 v0, v1, v2, e1, e2, b, x;
  mat3 A;
  bool if_intersect = false;
  closestIntersection.distance = std::numeric_limits<float>::max();

  for (size_t i=0; i<triangles.size(); i++){
    v0 = triangles[i].v0;
    v1 = triangles[i].v1;
    v2 = triangles[i].v2;
    e1 = v1 - v0;
    e2 = v2 - v0;
    b = start - v0;
    A = mat3 (-dir, e1, e2);
    x = glm::inverse( A ) * b;

    if (x.x>=0.0001 && x.y>=0 && x.z>=0 && (x.y + x.z)<=1){
      if_intersect = true;
      if (x.x < closestIntersection.distance){
        closestIntersection.position = start + x.x * dir;
        closestIntersection.distance = x.x;
        closestIntersection.triangleIndex = i;
      }
    }
  }
  return if_intersect;
}

vec3 DirectLight( const Intersection& i ) {
  vec3 n = triangles[i.triangleIndex].normal;        // surface normal
  vec3 r = glm::normalize (lightPos - i.position);   // unit directional vector from surface point to light source.
  float dot = glm::dot(n,r);
  float dist = glm::distance(lightPos, i.position);  // distance from surface point to light source.
  vec3 d = lightColor * fmaxf(dot, 0) / float (4 * M_PI * dist * dist);
  
  Intersection newIntersection;
  if (ClosestIntersection(i.position, r, triangles, newIntersection)){  //if there's intersection from point to light source.
     if (i.triangleIndex != newIntersection.triangleIndex){   // if the intersection is not on the same triangle
       if (newIntersection.distance < dist){
       return vec3(0,0,0); 
       }
     }     
  }
  return d; 
}