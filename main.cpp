#include <fstream>
#include <iostream>
#include <mutex>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include "image.h"
#include "scene.h"
#include "pathtracer.h"
#include "random.h"
#include "pm.h"
#include "meshsimp1.h"

using namespace std;

void PT (const json &scene);
void PPM (const json &scene);
void MeshSimplify (int argc, char **argv);

int main (int argc, char **argv)
{
	if (argc < 2)
	{
		cerr << "Usage: ./raytr <scene.json> | ./raytr m <mesh.json> <simp_ratio>" << endl ;
		exit(1);
	}

	if (string(argv[1]) == "m")
	{
		MeshSimplify(argc - 1, argv + 1);
		return 0;
	}

	auto scene_file = json::parse(ifstream(argv[1]));
	if (scene_file.count("pm_params"))
		PPM(scene_file);
	else if (scene_file.count("pt_params"))
		PT(scene_file);

	return 0;
}

void PPM (const json &j)
{
	PMParams params;
	params.LoadFrom(j.at("pm_params"));

	omp_set_num_threads(params.n_threads);
	Scene scene = Scene::Load(j);
	Image image(scene.Width() * 2, scene.Height(), (std::string("ph_last_") + j.at("pm_params").at("algo").get<string>()).c_str(), params.display);
	PM pm(scene, image, params);

	if (j.at("pm_params").at("algo").get<string>() == "ppm")
	{
		printf("algo: PPM\n");
		pm.PPM();
	}
	else
	{
		printf("algo: probabilistic PPM\n");
		pm.PPPM();
	}

	image.Show();
	image.Dump(-1);
	image.Show();
}

void PT (const json &j)
{
	PTParams params;
	params.LoadFrom(j.at("pt_params"));

	omp_set_num_threads(params.n_threads);
	Scene scene = Scene::Load(j);
	Image image(scene.Width() * 2, scene.Height(), "im_last", params.display);
	PathTracer pathTracer(scene, image, params);

	Ray ray;
	rand_data rd;
	rand_init(&rd);
	scene.GetCamera()->Project(&rd, 818, 1800 - 1043, ray);
	pathTracer.Trace(&rd, ray);

	pathTracer.PT();
	image.Show();
	image.Dump(-1);
	image.Show();
}

void MeshSimplify (int argc, char **argv)
{
    MeshSimplifier1 reducer(argv[1]);
    double target_ratio = argc > 2 ? std::stof(argv[2]) : 0.9;
    reducer.Run(target_ratio);
}
