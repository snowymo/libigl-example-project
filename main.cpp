#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <iostream>
#include <igl/triangle_triangle_adjacency.h>
#include <queue>
#include <set>
#include <igl/find.h>
#include <igl/remove_unreferenced.h>
#include <igl/embree/line_mesh_intersection.h>
#include <igl/decimate.h>
#include <igl/jet.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd N_vertices;
Eigen::MatrixXd N_faces;
Eigen::MatrixXd N_corners;


// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
	switch (key)
	{
	case '1':
		viewer.data().set_normals(N_faces);
		return true;
	case '2':
		viewer.data().set_normals(N_vertices);
		return true;
	case '3':
		viewer.data().set_normals(N_corners);
		return true;
	default: break;
	}
	return false;
}

int findFirstZero(Eigen::VectorXi v) {
	for (int i = 0; i < v.rows(); i++)
		if (v[i] == 0)
			return i;
	return -1;
}


void splitMesh(Eigen::MatrixXd V, Eigen::MatrixXi F, std::vector<Eigen::MatrixXd>& subVs, std::vector<Eigen::MatrixXi>& subFs, int bound = 65000) {

	// Plot the mesh
	//std::vector<igl::opengl::glfw::Viewer> viewers;
	igl::opengl::glfw::Viewer viewer;
	int pieces = 0;

	// triangle triangle adjacency
	Eigen::MatrixXi TT;
	igl::triangle_triangle_adjacency(F, TT);

	// face discover to record if the face is already into the queue
	//int *FD = new int[F.cols()]{ 0 };
	Eigen::VectorXi FD = Eigen::VectorXi::Zero(F.rows());

	// face queue
	std::queue<int> FQ;

	bool isOverSize = false;

	// set of face and vertex for current submesh
	std::set<int> FS;
	std::set<int> VS;

	// if there is still face not being discovered.
	while (!FD.isOnes()) {
		// find the first non-one
		int curIdx = findFirstZero(FD);
		if (curIdx < 0) break;

		// visit it
		FQ.push(curIdx);

		while (!FQ.empty()) {
			int curFace = FQ.front();
			FD(curFace) = 1;
			FQ.pop();
			//std::cout << "\ndealing " << curFace << " pushing neighbours ";

			FS.insert(curFace);
			for (int i = 0; i < F.cols(); i++)
				VS.insert(F(curFace, i));

			// get three adjacencies and push into the queue aka visit
			for (int i = 0; i < F.cols(); i++) {
				int neighbour = TT(curFace, i);

				// check if already discoved
				if ((neighbour != -1) && (FD(neighbour) == 0)) {
					FQ.push(neighbour);
					FD(neighbour) = 2;
					//std::cout << neighbour << "\t";
				}
				//std::cout << "\n";
			}

			// check the size of vertices and faces
			if (FS.size() > bound || VS.size() > bound) {
				isOverSize = true;
				break;
			}
		}

		if (isOverSize) {
			isOverSize = false;

			// deal with current submesh, turn set of faces to MatrixXd
			Eigen::MatrixXd subV;
			Eigen::MatrixXi subF(FS.size(), F.cols());
			int i = 0;
			for (std::set<int>::iterator it = FS.begin(); it != FS.end(); ++it, i++) {
				subF.row(i) = F.row(*it);
			}
			Eigen::VectorXi UJ;
			igl::remove_unreferenced(V, subF, subV, subF, UJ);

			subVs.push_back(subV);
			subFs.push_back(subF);

			FS.clear();
			VS.clear();

			// write down
			 	igl::writeOFF("split" + std::to_string(pieces++) + ".off", subV, subF);
				std::cout << "split " << pieces << "\n";

			 //Plot the mesh
				viewer.append_mesh();
				viewer.data().clear();
			 	viewer.data().set_mesh(subV, subF);
				double r, g, b;
				igl::colormap(igl::COLOR_MAP_TYPE_JET, (double)pieces / 10.0f, r, g, b);
				viewer.data().set_colors(Eigen::RowVector3d(r, g, b));
		}
	}

	if (FS.size() > 0) {
		// deal with current submesh, turn set of faces to MatrixXd
		Eigen::MatrixXd subV;
		Eigen::MatrixXi subF(FS.size(), F.cols());
		int i = 0;
		for (std::set<int>::iterator it = FS.begin(); it != FS.end(); ++it, i++) {
			subF.row(i) = F.row(*it);
		}
		Eigen::VectorXi UJ;
		igl::remove_unreferenced(V, subF, subV, subF, UJ);

		subVs.push_back(subV);
		subFs.push_back(subF);

		FS.clear();
		VS.clear();

		// write down
		igl::writeOFF("split" + std::to_string(pieces++) + ".off", subV, subF);
		std::cout << "split once\n";

		//Plot the mesh
		viewer.append_mesh();
		viewer.data().clear();
		viewer.data().set_mesh(subV, subF);
		double r, g, b;
		igl::colormap(igl::COLOR_MAP_TYPE_JET, (double)pieces / 10.0f, r, g, b);
		viewer.data().set_colors(Eigen::RowVector3d(r, g, b));
	}

	viewer.launch();
}

void splitMesh(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd C, std::vector<Eigen::MatrixXd>& subVs, std::vector<Eigen::MatrixXi>& subFs, std::vector<Eigen::MatrixXd>& subCs, int bound /*= 64995*/)
{
	// triangle triangle adjacency
	Eigen::MatrixXi TT;
	igl::triangle_triangle_adjacency(F, TT);

	// face discover to record if the face is already into the queue
	//int *FD = new int[F.cols()]{ 0 };
	Eigen::VectorXi FD = Eigen::VectorXi::Zero(F.rows());

	// face queue
	std::queue<int> FQ;

	// if there is still face not being discovered.
	while (!FD.isOnes()) {
		// find the first non-one
		int curIdx = findFirstZero(FD);
		if (curIdx < 0) break;

		// visit it
		FQ.push(curIdx);

		// set of face and vertex for current submesh
		std::set<int> FS;
		std::set<int> VS;

		while (!FQ.empty()) {
			int curFace = FQ.front();
			FD(curFace) = 1;
			FQ.pop();
			//std::cout << "\ndealing " << curFace << " pushing neighbours ";

			FS.insert(curFace);
			for (int i = 0; i < F.cols(); i++)
				VS.insert(F(curFace, i));

			// get three adjacencies and push into the queue aka visit
			for (int i = 0; i < F.cols(); i++) {
				int neighbour = TT(curFace, i);

				// check if already discoved
				if ((neighbour != -1) && (FD(neighbour) == 0)) {
					FQ.push(neighbour);
					FD(neighbour) = 2;
					//std::cout << neighbour << "\t";
				}
				//std::cout << "\n";
			}

			// check the size of vertices and faces
			if (FS.size() > bound || VS.size() > bound)
				break;
		}

		// deal with current submesh, turn set of faces to MatrixXd
		Eigen::MatrixXd subV, subC;
		Eigen::MatrixXi subF(FS.size(), F.cols());
		Eigen::MatrixXi resF(FS.size(), F.cols());
		int i = 0;
		for (std::set<int>::iterator it = FS.begin(); it != FS.end(); ++it, i++) {
			subF.row(i) = F.row(*it);
		}
		Eigen::VectorXi UJ;
		igl::remove_unreferenced(V, subF, subV, resF, UJ);
		igl::remove_unreferenced(C, subF, subC, resF, UJ);

		subVs.push_back(subV);
		subCs.push_back(subC);
		subFs.push_back(resF);

	}
}

void splitMeshMergeApart(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &C, std::vector<Eigen::MatrixXd>& subVs, std::vector<Eigen::MatrixXi>& subFs, std::vector<Eigen::MatrixXd>& subCs, int bound /*= 64995*/)
{
	// triangle triangle adjacency
	Eigen::MatrixXi TT;
	igl::triangle_triangle_adjacency(F, TT);

	// face discover to record if the face is already into the queue
	//int *FD = new int[F.cols()]{ 0 };
	Eigen::VectorXi FD = Eigen::VectorXi::Zero(F.rows());

	// face queue
	std::queue<int> FQ;

	bool isOverSize = false;

	// set of face and vertex for current submesh
	std::set<int> FS;
	std::set<int> VS;

	// if there is still face not being discovered.
	while (!FD.isOnes()) {
		// find the first non-one
		int curIdx = findFirstZero(FD);
		if (curIdx < 0) break;

		// visit it
		FQ.push(curIdx);

		while (!FQ.empty()) {
			int curFace = FQ.front();
			FD(curFace) = 1;
			FQ.pop();
			//std::cout << "\ndealing " << curFace << " pushing neighbours ";

			FS.insert(curFace);
			for (int i = 0; i < F.cols(); i++)
				VS.insert(F(curFace, i));

			// get three adjacencies and push into the queue aka visit
			for (int i = 0; i < F.cols(); i++) {
				int neighbour = TT(curFace, i);

				// check if already discovered
				if ((neighbour != -1) && (FD(neighbour) == 0)) {
					FQ.push(neighbour);
					FD(neighbour) = 2;
					//std::cout << neighbour << "\t";
				}
				//std::cout << "\n";
			}

			// check the size of vertices and faces
			if (FS.size() > bound || VS.size() > bound) {
				isOverSize = true;
				break;
			}
				
		}

		if (isOverSize) {
			// deal with current submesh, turn set of faces to MatrixXd
			Eigen::MatrixXd subV, subC;
			Eigen::MatrixXi subF(FS.size(), F.cols());
			Eigen::MatrixXi resF(FS.size(), F.cols());
			int i = 0;
			for (std::set<int>::iterator it = FS.begin(); it != FS.end(); ++it, i++) {
				subF.row(i) = F.row(*it);
			}
			Eigen::VectorXi UJ;
			igl::remove_unreferenced(V, subF, subV, resF, UJ);
			igl::remove_unreferenced(C, subF, subC, resF, UJ);

			subVs.push_back(subV);
			subCs.push_back(subC);
			subFs.push_back(resF);

			FS.clear();
			VS.clear();
		}
	}
	if (FS.size() > 0) {
		// deal with current submesh, turn set of faces to MatrixXd
		Eigen::MatrixXd subV, subC;
		Eigen::MatrixXi subF(FS.size(), F.cols());
		Eigen::MatrixXi resF(FS.size(), F.cols());
		int i = 0;
		for (std::set<int>::iterator it = FS.begin(); it != FS.end(); ++it, i++) {
			subF.row(i) = F.row(*it);
		}
		Eigen::VectorXi UJ;
		igl::remove_unreferenced(V, subF, subV, resF, UJ);
		igl::remove_unreferenced(C, subF, subC, resF, UJ);

		subVs.push_back(subV);
		subCs.push_back(subC);
		subFs.push_back(resF);

		FS.clear();
		VS.clear();
	}
}

void Resample(const Eigen::MatrixXd & V_target, const Eigen::MatrixXi & F_target, const Eigen::MatrixXd &C_target, Eigen::MatrixXd &V_source, Eigen::MatrixXi &F_source, Eigen::MatrixXd & C_result) {
	const size_t threshold = 65000;
	
	Eigen::VectorXi  J;//??
	bool ret = igl::decimate(V_target, F_target, 200, V_source, F_source, J);
	std::cout << "decimate:" << ret << "\n";
	C_result.resize(V_source.rows(), 3);
	if (ret) {
		// Compute per-vertex normals
		Eigen::MatrixXd N_source;
		igl::per_vertex_normals(V_source, F_source, N_source);
		Eigen::MatrixXd FandBary = igl::embree::line_mesh_intersection(V_source, N_source, V_target, F_target);
/*
		for vdist0 in vdist0s :
		new_dist = np.zeros(V.shape[0])
			fdist = vdist0[F0[ProjF0]] # v dist for verts from projected faces
			for i, f in enumerate(fdist) :
				new_dist[i] = np.dot(ProjB[i, :], f) # inner product
				results.append(new_dist)
				*/

		// project the *_source to get the color
		for (int i = 0; i < V_source.rows(); i++) {
			int faceId = FandBary(i, 0);
			double b1 = FandBary(i, 1);
			double b2 = FandBary(i, 2);
			double b0 = 1 - b1 - b2;
			std::cout << "face id " << faceId << "\t b0 " << b0 << "\t b1 " << b1 << "\t b2 " << b2 << "\n";
			Eigen::Vector3d c0 = C_target.row(F_target(faceId, 0));
			Eigen::Vector3d c1 = C_target.row(F_target(faceId, 1));
			Eigen::Vector3d c2 = C_target.row(F_target(faceId, 2));
			Eigen::Vector3d finalc = c0 * b0 + c1 * b1 + c2 * b2;
			C_result.row(i) = finalc;
		}
	}
	
}


int main(int argc, char *argv[])
{
	// Load a mesh in OFF format
	//igl::readOFF( "D:\\Projects\\libigl\\tutorial\\shared\\fandisk.off", V, F);
	//igl::readOBJ("D:\\Projects\\libigl\\tutorial\\shared\\brokenFace.obj", V, F);
	igl::readOBJ("D:\\Projects\\libigl\\tutorial\\shared\\armadillo.obj", V, F);
	std::cout << V.size() << " vertices\t" << F.size() << " faces\n";
	std::vector<Eigen::MatrixXd> subVs; std::vector<Eigen::MatrixXi> subFs;


	//
	Eigen::VectorXd Z = V.col(2);
	Eigen::MatrixXd C;
	igl::jet(Z, true, C);

	Eigen::MatrixXd V_source; Eigen::MatrixXi F_source; Eigen::MatrixXd  C_result;
	//Resample(V, F, C, V_source, F_source, C_result);

	splitMesh(V, F, subVs, subFs, 5000);

	// Compute per-face normals
	
	igl::per_face_normals(V, F, N_faces);

	// Compute per-vertex normals
	igl::per_vertex_normals(V, F, N_vertices);

	// Compute per-corner normals, |dihedral angle| > 20 degrees --> crease
	igl::per_corner_normals(V, F, 20, N_corners);

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.callback_key_down = &key_down;
	viewer.data().show_lines = false;
	viewer.data().set_mesh(V_source, F_source);
	//viewer.data().set_normals(N_faces);
	viewer.data().set_colors(C_result);

// 	viewer.append_mesh();
// 	viewer.data().clear();
// 	viewer.data().set_mesh(V, F);
// 	//viewer.data().set_normals(N_faces);
// 	viewer.data().set_colors(C);
// 	std::cout <<
// 		"Press '1' for per-face normals." << std::endl <<
// 		"Press '2' for per-vertex normals." << std::endl <<
// 		"Press '3' for per-corner normals." << std::endl;
	viewer.launch();
}