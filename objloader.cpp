/*
 * obj_converter.cc: Load obj & ext
 */

#include "objloader.h"
#include "ext/tiny_obj_loader.h"
#include "primitives.h"
#include "shade.h"

using nlohmann::json;
using namespace tinyobj;
using namespace std;


class VecPool: public Destroyable
{
	vector<Vec3> vec;
	size_t cp;
public:
	VecPool (size_t size)
	{
		vec.resize(size);
		cp = 0;
	}
	~ VecPool () {}
	CPVec3 Add (const Vec3 &v)
	{
		assert(cp < vec.size());
		vec[cp++] = v;
		return &vec[cp - 1];
	}
};


class VecMap
{
	vector<pair<Vec3, size_t>> axis_;
	vector<size_t> orig_to_sorted_;
	vector<Vec3> storage_;

public:

	VecMap (const vector<Vec3> &orig, Vec3 v0)
	{
		size_t n = orig.size();
		axis_.resize(n);
		orig_to_sorted_.resize(n);
		storage_ = vector<Vec3>(n, v0);

		for (size_t i = 0; i < n; ++i)
			axis_[i] = make_pair(orig[i], i);

		auto vecLt = [](const pair<Vec3, size_t> &a_, const pair<Vec3, size_t> &b_) {
			const auto &a = a_.first;
			const auto &b = b_.first;
			if (fabs(a.x - b.x) > eps) return a.x < b.x;
			if (fabs(a.y - b.y) > eps) return a.y < b.y;
			if (fabs(a.z - b.z) > eps) return a.z < b.z;
			return false;
		};

		sort(axis_.begin(), axis_.end(), vecLt);

		for (size_t i = 0; i < n; ++i)
		{
			if (i && !vecLt(axis_[i - 1], axis_[i])) // eq
				orig_to_sorted_[axis_[i].second] = orig_to_sorted_[axis_[i - 1].second];
			else
				orig_to_sorted_[axis_[i].second] = i;
		}
	}

	Vec3 &operator[] (size_t p) { return storage_.at(orig_to_sorted_[p]); }
};


int LoadMesh (const string &objPath,
			  const json &modifiers,
			  int interpolate_normal,
			  vector<Primitive*> &prims,
			  vector<Destroyable*> &resPool)
{
	vector<shape_t> shapes;
	vector<material_t> materials;

	std::string err;
	bool ret = LoadObj(shapes, materials, err, objPath.c_str());

	if (!err.empty())
	{
		// `err` may contain warning message.
		std::cerr << "When loading mesh: " << err << std::endl;
	}

	if (!ret)
	{
		std::cerr << "Mesh failed to load\n";
		return 1;
	}

	std::cout << "# of shapes    : " << shapes.size() << std::endl;
	std::cout << "# of materials : " << materials.size() << std::endl;

	vector<PSurfaceProperty> surfaceProps;
	vector<std::map<string, string>> surfaceTags;
	for (material_t &mtl: materials)
	{
#define BC(arr) RGBSpectrum(arr[0], arr[1], arr[2])
#ifndef SPECTRAL
		auto diffuse = BC(mtl.diffuse);
		auto specular = BC(mtl.specular);
		auto transmittance = BC(mtl.transmittance);
		auto emission = BC(mtl.emission);
#else
#define PropToSpectrum(id, idf, typ) (mtl.unknown_parameter.count("S" #id) ? \
		SampledSpectrum::Parse(json::parse(mtl.unknown_parameter.at("S" #id))["S"], typ) : \
		SampledSpectrum::FromRGB(BC(mtl.idf), typ))
		Spectrum diffuse = PropToSpectrum(Kd, diffuse, 1);
		Spectrum specular = PropToSpectrum(Ks, specular, 1);
		Spectrum transmittance = PropToSpectrum(Ki, transmittance, 1);
		Spectrum emission = PropToSpectrum(Ke, emission, 1);
#endif
		double Ns = mtl.shininess;
		if (Ns == 1)
			Ns = 0;

		if (transmittance.nonzero())
			cerr << "Warning: transmittance not supported in material " << mtl.name << "\n";

		surfaceProps.push_back(make_shared<ShadingProperty>(
				diffuse, emission, specular, Ns, mtl.ior,
				mtl.diffuse_texname, mtl.specular_texname, mtl.bump_texname,
				&mtl.unknown_parameter
		));
		surfaceTags.push_back(mtl.unknown_parameter);

	/*
		printf("material[%ld].name = %s\n", i, materials[i].name.c_str());
		printf("  material.Ka = (%f, %f ,%f)\n", materials[i].ambient[0], materials[i].ambient[1],
			   materials[i].ambient[2]);
		printf("  material.Kd = (%f, %f ,%f)\n", materials[i].diffuse[0], materials[i].diffuse[1],
			   materials[i].diffuse[2]);
		printf("  material.Ks = (%f, %f ,%f)\n", materials[i].specular[0], materials[i].specular[1],
			   materials[i].specular[2]);
		printf("  material.Tr = (%f, %f ,%f)\n", materials[i].transmittance[0], materials[i].transmittance[1],
			   materials[i].transmittance[2]);
		printf("  material.Ke = (%f, %f ,%f)\n", materials[i].emission[0], materials[i].emission[1],
			   materials[i].emission[2]);
		printf("  material.Ns = %f\n", materials[i].shininess);
		printf("  material.Ni = %f\n", materials[i].ior);
		printf("  material.dissolve = %f\n", materials[i].dissolve);
		printf("  material.illum = %d\n", materials[i].illum);
		printf("  material.map_Ka = %s\n", materials[i].ambient_texname.c_str());
		printf("  material.map_Kd = %s\n", materials[i].diffuse_texname.c_str());
		printf("  material.map_Ks = %s\n", materials[i].specular_texname.c_str());
		printf("  material.map_Ns = %s\n", materials[i].specular_highlight_texname.c_str());
		std::map<std::string, std::string>::const_iterator it(materials[i].unknown_parameter.begin());
		std::map<std::string, std::string>::const_iterator itEnd(materials[i].unknown_parameter.end());
		for (; it != itEnd; it++)
		{
			printf("  material.%s = %s\n", it->first.c_str(), it->second.c_str());
		}
		printf("\n");
		 */
	}

	for (size_t i = 0; i < shapes.size(); ++i)
	{
		json jmod;
		if (modifiers.count(shapes[i].name))
			jmod = modifiers[shapes[i].name];

		if (jmod.count("excluded"))
			continue;

		auto transform = jmod.count("Tf") ? VecTf(jmod["Tf"]) : VecTf::Eye();

		printf("shape[%ld].name = %s\n", i, shapes[i].name.c_str());
		printf("Size of shape[%ld].mesh_indices: %ld\n", i, shapes[i].mesh.indices.size());
		printf("Size of shape[%ld].material_ids: %ld\n", i, shapes[i].mesh.material_ids.size());

		// Init temp vertex buffer
		vector<Vec3> vertices;
		printf("shape[%ld].vertices: %ld\n", i, shapes[i].mesh.positions.size());
		for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++)
		{
			Vec3 vv(shapes[i].mesh.positions.at(3 * v + 0),
					shapes[i].mesh.positions.at(3 * v + 1),
					shapes[i].mesh.positions.at(3 * v + 2));
			vertices.push_back(transform(vv));
		}
		size_t n_vertices = vertices.size();

		// Init normal vector & texcoord pool
		vector<CPVec3> normals, texcoords;
		const auto &mesh_indices = shapes[i].mesh.indices;

		auto n_normals = shapes[i].mesh.normals.size();
		auto n_texcrds = shapes[i].mesh.texcoords.size();
		assert(n_normals <= 3 || n_normals == n_vertices * 3);
		assert(n_texcrds <= 2 || n_texcrds == n_vertices * 2);
		VecPool *vecPool = new VecPool(n_normals / 3 + n_texcrds / 2);

		if (jmod.count("re_interpolate"))
		{
			clog << "\tRe-interpolating normals\n";
			// Re-interpolate normals

			// tiny_obj_loader creates copies for same vertex in different faces. So we need a map.
			VecMap Nu(vertices, Vec3(0, 0, 0));
			for (size_t f = 0; f < shapes[i].mesh.indices.size(); f += 3)
			{
				size_t a[3] = {mesh_indices[f + 2], mesh_indices[f + 1], mesh_indices[f + 0]};
				Vec3 v[3] = {vertices[a[0]], vertices[a[1]], vertices[a[2]]};
				Vec3 nu = (v[1] - v[0]).cross(v[2] - v[0]);
				for (int d = 0; d < 3; ++d)
				{
					auto ea = v[(d + 1) % 3] - v[d];
					auto eb = v[(d + 2) % 3] - v[d];
					long double ang = acos(ea.normalized().dot(eb.normalized()));
					assert(!std::isnan(ang));
					Nu[a[d]] += ang * nu; // == ang * area * n
				}
			}
			for (size_t i = 0; i < n_vertices; ++i)
				normals.push_back(vecPool->Add(Nu[i] / Nu[i].norm()));
		}
		else if (n_normals == n_vertices * 3 && interpolate_normal)
		{
			const auto &m_normals = shapes[i].mesh.normals;
			auto inv_tf = transform.Inverse();
			for (size_t v = 0; v < m_normals.size() / 3; v++)
			{
				Vec3 n(-m_normals.at(v * 3 + 0), -m_normals.at(v * 3 + 1), -m_normals.at(v * 3 + 2));
				normals.push_back(vecPool->Add(inv_tf.TVec(n).normalized()));
			}
		}
		else
		{
			// Storage waste is O(1).
			normals = vector<CPVec3>(n_vertices, nullptr);
		}

		if (n_texcrds == n_vertices * 2)
		{
			const auto &m_texcoords = shapes[i].mesh.texcoords;
			for (size_t v = 0; v < m_texcoords.size() / 2; v++)
			{
				Vec3 t(m_texcoords.at(v * 2 + 0), m_texcoords.at(v * 2 + 1), 0);
				texcoords.push_back(vecPool->Add(t));
			}
		}
		else
		{
			texcoords = vector<CPVec3>(n_vertices, nullptr);
		}

		// Set the owner of the pool
		resPool.push_back(vecPool);

		assert((mesh_indices.size() % 3) == 0);
		assert(normals.size() == n_vertices || normals.size() <= 1);
		for (size_t f = 0; f < mesh_indices.size() / 3; f++)
		{
			Triangle *tri;

			// In obj file normals point to the outside of object, but I have defined
			// normal in triangle in the opposite way. I decided to modify obj file instead.

			int a = mesh_indices[3 * f + 2], b = mesh_indices[3 * f + 1], c = mesh_indices[3 * f + 0];
			auto Va = vertices[a];
			auto Vb = vertices[b];
			auto Vc = vertices[c];
			if ((Va - Vb).cross(Vc - Va).norm() < eps)
			{
				cerr << "*DT ";
				//cerr << "Degenerated triangle found\n";
				continue;
			}
			tri = new Triangle(
					vertices[a], vertices[b], vertices[c],
					normals[a], normals[b], normals[c],
					texcoords[a], texcoords[b], texcoords[c],
					surfaceProps.at((size_t) shapes[i].mesh.material_ids[f]));

			if (jmod.count("ldrange"))
			{
				tri->SetUnidirectedLight(
						Vec3::Parse(jmod["ldrange"][0]),
						Vec3::Parse(jmod["ldrange"][1]));
			}
			else if (jmod.count("conc"))
			{
				tri->SetLightConcentration(jmod["conc"].get<double>());
			}

			prims.push_back(tri);
		}
		
		printf("\n");
	}

}
