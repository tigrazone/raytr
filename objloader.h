#ifndef OBJLOADER_H

#define OBJLOADER_H

#include "common.h"
#include "objekt.h"

// Load .obj scenes

int LoadMesh (const string &objPath,
			  const json &modifiers,
			  int interpolate_normal,
			  vector<Primitive*> &prims,
			  vector<Destroyable*> &resPool);

#endif